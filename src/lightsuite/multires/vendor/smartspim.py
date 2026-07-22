"""Build multiresolution pair manifests from SmartSPIM / ASI exports."""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

from lightsuite.multires.manifest import save_pair_manifest
from lightsuite.multires.models import MANIFEST_FORMAT, ManifestVolumeSpec, MultiresPairManifest
from lightsuite.multires.volume import discover_volume_shape

# ASI stage X/Y/Z table values are hundredths of a micron (0.01 µm).
STAGE_COORD_SCALE_UM = 0.01


@dataclass(frozen=True)
class SmartspimGeometryConfig:
    """How ASI stage coordinates map into physical overview / ROI frames."""

    stage_coord_scale_um: float = STAGE_COORD_SCALE_UM
    stage_xy_is_center: bool = True
    lateral_flip: tuple[int, int] = (1, 1)
    # Map ROI origin Y from overview tile-row stage anchors (recommended for overview↔ROI pairs).
    anchor_roi_y_to_overview_stage: bool = False
    # Optional calibration shift applied to the solved origin (µm), after stage placement.
    stage_origin_offset_um: tuple[float, float, float] = (0.0, 0.0, 0.0)


@dataclass(frozen=True)
class SmartspimScanMeta:
    """Parsed fields from a SmartSPIM ``metadata.txt`` export."""

    objective: str
    hres: int
    vres: int
    um_per_pix: float
    z_step_um: float
    tile_centers_stage: list[tuple[float, float, float]]
    tile_num_images: list[int]


def parse_smartspim_metadata(path: Path) -> SmartspimScanMeta:
    """Parse ASI SmartSPIM ``metadata.txt`` (tab-separated export)."""
    path = path.expanduser().resolve()
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    if len(lines) < 2:
        msg = f"SmartSPIM metadata too short: {path}"
        raise ValueError(msg)

    header = lines[0].split("\t")
    values = lines[1].split("\t")
    field_map = {key.strip(): values[index].strip() for index, key in enumerate(header) if key.strip()}

    try:
        objective = field_map["Obj"]
        hres = int(float(field_map["H_Res"]))
        vres = int(float(field_map["V_Res"]))
        um_per_pix = float(field_map["um/pix"])
        z_step_um = float(field_map["Z step (um)"])
    except KeyError as exc:
        msg = f"Missing SmartSPIM header field {exc!s} in {path}"
        raise ValueError(msg) from exc

    tile_centers_stage: list[tuple[float, float, float]] = []
    tile_num_images: list[int] = []
    in_z_table = False
    for raw in lines[2:]:
        line = raw.strip()
        if not line:
            continue
        if line.startswith("wavelength "):
            in_z_table = False
            continue
        if line == "z\tfocus":
            in_z_table = True
            continue
        parts = line.split("\t")
        if in_z_table:
            continue
        if len(parts) >= 9 and parts[0] not in {"X", "Wavelength"}:
            try:
                tile_centers_stage.append((float(parts[0]), float(parts[1]), float(parts[2])))
                tile_num_images.append(int(float(parts[8])))
            except ValueError:
                continue

    if not tile_centers_stage:
        msg = f"No tile stage positions found in {path}"
        raise ValueError(msg)
    if not tile_num_images:
        msg = f"No NumImages values found in tile table: {path}"
        raise ValueError(msg)

    return SmartspimScanMeta(
        objective=objective,
        hres=hres,
        vres=vres,
        um_per_pix=um_per_pix,
        z_step_um=z_step_um,
        tile_centers_stage=tile_centers_stage,
        tile_num_images=tile_num_images,
    )


def _stage_scale(geometry: SmartspimGeometryConfig) -> float:
    return float(geometry.stage_coord_scale_um)


def _direction_from_geometry(geometry: SmartspimGeometryConfig) -> list[float]:
    f0, f1 = geometry.lateral_flip
    return [float(f0), 0.0, 0.0, 0.0, float(f1), 0.0, 0.0, 0.0, 1.0]


def _stage_z_origin_um(tile_centers_stage: list[tuple[float, float, float]]) -> float:
    """Return mean tile-table Z as physical origin (µm).

    Unlike X/Y stage coordinates, SmartSPIM tile Z is already in microns.
    """
    zs = [coord[2] for coord in tile_centers_stage]
    return float(sum(zs) / len(zs))


def _resolve_num_images(tile_num_images: list[int]) -> int:
    """Return the metadata NumImages value used for Z extent."""
    unique = sorted(set(tile_num_images))
    if len(unique) != 1:
        msg = f"Expected a single NumImages value across tiles, got {unique}"
        raise ValueError(msg)
    return unique[0]


def _stage_center_origin_um(
    tile_centers_stage: list[tuple[float, float, float]],
    *,
    shape_yx: tuple[int, int],
    um_per_pix: float,
    geometry: SmartspimGeometryConfig,
) -> tuple[float, float]:
    """Place volume assuming stage X/Y are tile centers in stage table units."""
    scale = _stage_scale(geometry)
    xs = [coord[0] for coord in tile_centers_stage]
    ys = [coord[1] for coord in tile_centers_stage]
    center_x_um = sum(xs) / len(xs) * scale
    center_y_um = sum(ys) / len(ys) * scale
    ny, nx = shape_yx
    f0, f1 = geometry.lateral_flip
    if geometry.stage_xy_is_center:
        origin_x = center_x_um - f0 * um_per_pix * (nx - 1) / 2.0
        origin_y = center_y_um - f1 * um_per_pix * (ny - 1) / 2.0
    else:
        origin_x = min(xs) * scale
        origin_y = min(ys) * scale
    return origin_x, origin_y


def _overview_stage_y_anchors(
    overview_tiles: list[tuple[float, float, float]],
) -> tuple[float, float]:
    ys = sorted({coord[1] for coord in overview_tiles})
    if len(ys) < 2:
        msg = "Need at least two distinct overview tile Y stage values for ROI Y anchoring"
        raise ValueError(msg)
    return ys[0], ys[-1]


def _physical_y_at_overview_stage_row(
    stage_y: float,
    *,
    overview_tiles: list[tuple[float, float, float]],
    shape_yx: tuple[int, int],
    um_per_pix: float,
    geometry: SmartspimGeometryConfig,
) -> float:
    """Map a stage-table Y coordinate into overview physical Y (µm)."""
    y_min_stage, y_max_stage = _overview_stage_y_anchors(overview_tiles)
    if y_max_stage <= y_min_stage:
        msg = "Overview tile Y stage values must span a positive range"
        raise ValueError(msg)

    _, origin_y = _stage_center_origin_um(
        overview_tiles,
        shape_yx=shape_yx,
        um_per_pix=um_per_pix,
        geometry=geometry,
    )
    ny, _nx = shape_yx
    f1 = geometry.lateral_flip[1]
    y_index_extent = f1 * um_per_pix * (ny - 1)
    y_at_min_stage = origin_y
    y_at_max_stage = origin_y + y_index_extent
    frac = (stage_y - y_min_stage) / (y_max_stage - y_min_stage)
    return y_at_min_stage + frac * (y_at_max_stage - y_at_min_stage)


def _roi_center_y_um(
    roi_tiles: list[tuple[float, float, float]],
    *,
    overview_tiles: list[tuple[float, float, float]] | None,
    overview_shape_yx: tuple[int, int] | None,
    overview_um_per_pix: float | None,
    overview_geometry: SmartspimGeometryConfig | None,
    roi_um_per_pix: float,
    roi_shape_yx: tuple[int, int],
    geometry: SmartspimGeometryConfig,
) -> float:
    ys = [coord[1] for coord in roi_tiles]
    mean_stage_y = sum(ys) / len(ys)
    if (
        geometry.anchor_roi_y_to_overview_stage
        and overview_tiles is not None
        and overview_shape_yx is not None
        and overview_um_per_pix is not None
        and overview_geometry is not None
    ):
        return _physical_y_at_overview_stage_row(
            mean_stage_y,
            overview_tiles=overview_tiles,
            shape_yx=overview_shape_yx,
            um_per_pix=overview_um_per_pix,
            geometry=overview_geometry,
        )

    scale = _stage_scale(geometry)
    center_y_um = mean_stage_y * scale
    if geometry.stage_xy_is_center:
        return center_y_um
    _, ny = roi_shape_yx
    return center_y_um + roi_um_per_pix * (ny - 1) / 2.0


def volume_spec_from_smartspim_export(
    *,
    label: str,
    volume_path: Path,
    metadata_path: Path,
    tile_centers_stage: list[tuple[float, float, float]] | None = None,
    geometry: SmartspimGeometryConfig | None = None,
    overview_tiles_for_roi_y: list[tuple[float, float, float]] | None = None,
    overview_shape_yx: tuple[int, int] | None = None,
    overview_um_per_pix: float | None = None,
    overview_geometry: SmartspimGeometryConfig | None = None,
) -> ManifestVolumeSpec:
    """Build a manifest volume spec from a stitched SmartSPIM stack folder."""
    volume_path = volume_path.expanduser().resolve()
    geometry = geometry or SmartspimGeometryConfig()
    meta = parse_smartspim_metadata(metadata_path)
    nz, ny, nx = discover_volume_shape(volume_path)
    centers = tile_centers_stage or meta.tile_centers_stage
    num_images = _resolve_num_images(meta.tile_num_images)
    if nz != num_images:
        msg = (
            f"Volume plane count ({nz}) does not match metadata NumImages ({num_images}) "
            f"for {volume_path}"
        )
        raise ValueError(msg)

    origin_x, origin_y = _stage_center_origin_um(
        centers,
        shape_yx=(ny, nx),
        um_per_pix=meta.um_per_pix,
        geometry=geometry,
    )
    if label == "roi" and geometry.anchor_roi_y_to_overview_stage:
        center_y_um = _roi_center_y_um(
            centers,
            overview_tiles=overview_tiles_for_roi_y,
            overview_shape_yx=overview_shape_yx,
            overview_um_per_pix=overview_um_per_pix,
            overview_geometry=overview_geometry,
            roi_um_per_pix=meta.um_per_pix,
            roi_shape_yx=(ny, nx),
            geometry=geometry,
        )
        origin_y = center_y_um - (ny - 1) / 2.0 * meta.um_per_pix

    origin_z = _stage_z_origin_um(centers)
    offset_x, offset_y, offset_z = geometry.stage_origin_offset_um
    origin_x += float(offset_x)
    origin_y += float(offset_y)
    origin_z += float(offset_z)
    spacing = (meta.um_per_pix, meta.um_per_pix, meta.z_step_um)
    return ManifestVolumeSpec(
        volume_path=str(volume_path),
        shape_zyx=[nz, ny, nx],
        spacing_um=[float(v) for v in spacing],
        origin_um=[origin_x, origin_y, origin_z],
        direction=_direction_from_geometry(geometry),
    )


def build_smartspim_pair_manifest(
    *,
    sample_name: str,
    pair_label: str,
    overview_path: Path,
    roi_path: Path,
    overview_meta_path: Path,
    roi_meta_path: Path,
    overview_tile_centers_stage: list[tuple[float, float, float]] | None = None,
    roi_tile_centers_stage: list[tuple[float, float, float]] | None = None,
    overview_geometry: SmartspimGeometryConfig | None = None,
    roi_geometry: SmartspimGeometryConfig | None = None,
    output_manifest_path: Path | None = None,
) -> MultiresPairManifest:
    """Convert stitched SmartSPIM exports into a LightSuite multires pair manifest."""
    overview_path = overview_path.expanduser().resolve()
    roi_path = roi_path.expanduser().resolve()
    overview_meta_path = overview_meta_path.expanduser().resolve()
    roi_meta_path = roi_meta_path.expanduser().resolve()
    overview_geometry = overview_geometry or SmartspimGeometryConfig()
    roi_geometry = roi_geometry or SmartspimGeometryConfig()

    overview_meta = parse_smartspim_metadata(overview_meta_path)
    overview_tiles = overview_tile_centers_stage or overview_meta.tile_centers_stage
    _, overview_ny, overview_nx = discover_volume_shape(overview_path)

    overview_spec = volume_spec_from_smartspim_export(
        label="overview",
        volume_path=overview_path,
        metadata_path=overview_meta_path,
        tile_centers_stage=overview_tiles,
        geometry=overview_geometry,
    )
    roi_spec = volume_spec_from_smartspim_export(
        label="roi",
        volume_path=roi_path,
        metadata_path=roi_meta_path,
        tile_centers_stage=roi_tile_centers_stage,
        geometry=roi_geometry,
        overview_tiles_for_roi_y=overview_tiles,
        overview_shape_yx=(overview_ny, overview_nx),
        overview_um_per_pix=overview_meta.um_per_pix,
        overview_geometry=overview_geometry,
    )

    manifest = MultiresPairManifest(
        format=MANIFEST_FORMAT,
        sample_name=sample_name,
        pair_label=pair_label,
        overview=overview_spec,
        roi=roi_spec,
        provenance={
            "microscope": "smartspim",
            "source_overview": str(overview_path),
            "source_roi": str(roi_path),
            "overview_meta": str(overview_meta_path),
            "roi_meta": str(roi_meta_path),
            "conversion": "lightsuite.multires.vendor.smartspim",
            "stage_coord_scale_um": str(overview_geometry.stage_coord_scale_um),
            "overview_lateral_flip": ",".join(str(v) for v in overview_geometry.lateral_flip),
            "roi_lateral_flip": ",".join(str(v) for v in roi_geometry.lateral_flip),
            "roi_y_anchored_to_overview_stage": str(roi_geometry.anchor_roi_y_to_overview_stage),
            "roi_stage_origin_offset_um": ",".join(str(v) for v in roi_geometry.stage_origin_offset_um),
            "z_origin": "stage_table",
            "conversion_timestamp": datetime.now(timezone.utc).isoformat(),
        },
    )

    if output_manifest_path is not None:
        save_pair_manifest(manifest, output_manifest_path)

    return manifest
