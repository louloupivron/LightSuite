"""Build multiresolution pair manifests from SmartSPIM / ASI exports."""

from __future__ import annotations

import math
from dataclasses import dataclass, replace
from datetime import datetime, timezone
from pathlib import Path
from typing import Literal

from lightsuite.multires.manifest import save_pair_manifest
from lightsuite.multires.models import MANIFEST_FORMAT, ManifestVolumeSpec, MultiresPairManifest
from lightsuite.multires.volume import discover_volume_shape

# ASI stage X/Y/Z table values are hundredths of a micron (0.01 µm).
STAGE_COORD_SCALE_UM = 0.01

RoiXyPlacement = Literal["mean_center", "upper_left_tile_corner"]
XyPlacement = RoiXyPlacement  # alias; applies to overview and ROI alike


@dataclass(frozen=True)
class SmartspimGeometryConfig:
    """How ASI stage coordinates map into physical overview / ROI frames."""

    stage_coord_scale_um: float = STAGE_COORD_SCALE_UM
    stage_xy_is_center: bool = True
    stage_z_is_center: bool = False
    lateral_flip: tuple[int, int] = (1, 1)
    # XY frame: mean tile center (default) or upper-left tile corner + stitched extent.
    xy_placement: XyPlacement = "mean_center"
    roi_xy_placement: RoiXyPlacement | None = None  # deprecated alias for xy_placement on ROI builds
    # Map ROI origin Y from overview tile-row stage anchors (recommended for overview↔ROI pairs).
    anchor_roi_y_to_overview_stage: bool = False
    # Map ROI origin X from overview tile-column stage anchors (auto-enabled with Y for mosaics).
    anchor_roi_x_to_overview_stage: bool = False
    # Correct mean_center placement when tile stage pitch is smaller than tile width (overlapping mosaics).
    apply_mosaic_pitch_correction: bool = False
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


def _stage_z_origin_um(
    tile_centers_stage: list[tuple[float, float, float]],
    *,
    nz: int,
    z_step_um: float,
    geometry: SmartspimGeometryConfig,
) -> float:
    """Return physical Z origin (µm) from tile-table Z values.

    Unlike X/Y stage coordinates, SmartSPIM tile Z is already in microns.
    When ``stage_z_is_center`` is true, the table value is the stack center.
    """
    zs = [coord[2] for coord in tile_centers_stage]
    z_mean = float(sum(zs) / len(zs))
    if geometry.stage_z_is_center:
        return z_mean - (nz - 1) / 2.0 * z_step_um
    return z_mean


def _resolve_num_images(tile_num_images: list[int]) -> int:
    """Return the metadata NumImages value used for Z extent."""
    unique = sorted(set(tile_num_images))
    if len(unique) != 1:
        msg = f"Expected a single NumImages value across tiles, got {unique}"
        raise ValueError(msg)
    return unique[0]


def _select_upper_left_tile(
    tile_centers_stage: list[tuple[float, float, float]],
    lateral_flip: tuple[int, int],
) -> tuple[float, float, float]:
    """Return the tile whose center anchors image index (0, 0) for the given lateral flip."""
    xs = sorted({coord[0] for coord in tile_centers_stage})
    ys = sorted({coord[1] for coord in tile_centers_stage})
    if not xs or not ys:
        msg = "Need at least one tile stage position"
        raise ValueError(msg)

    col_x = xs[0]
    _f0, f1 = lateral_flip
    row_y = ys[-1] if f1 < 0 else ys[0]
    for tile in tile_centers_stage:
        if math.isclose(tile[0], col_x, abs_tol=1.0) and math.isclose(tile[1], row_y, abs_tol=1.0):
            return tile

    msg = (
        "Could not find upper-left tile for stage grid "
        f"(col_x={col_x}, row_y={row_y}, tiles={len(tile_centers_stage)})"
    )
    raise ValueError(msg)


def _resolved_xy_placement(geometry: SmartspimGeometryConfig) -> XyPlacement:
    if geometry.roi_xy_placement is not None:
        return geometry.roi_xy_placement
    return geometry.xy_placement


def _roi_origin_from_upper_left_tile_corner(
    tile_centers_stage: list[tuple[float, float, float]],
    *,
    hres: int,
    vres: int,
    um_per_pix: float,
    geometry: SmartspimGeometryConfig,
) -> tuple[float, float]:
    """Place ROI origin at the physical corner of the upper-left mosaic tile."""
    tile = _select_upper_left_tile(tile_centers_stage, geometry.lateral_flip)
    scale = _stage_scale(geometry)
    center_x_um = tile[0] * scale
    center_y_um = tile[1] * scale
    tile_width_um = float(hres) * um_per_pix
    tile_height_um = float(vres) * um_per_pix
    f0, f1 = geometry.lateral_flip
    origin_x = center_x_um - f0 * tile_width_um / 2.0
    origin_y = center_y_um - f1 * tile_height_um / 2.0
    return origin_x, origin_y


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


def _mosaic_pitch_correction_um(
    tile_centers_stage: list[tuple[float, float, float]],
    *,
    hres: int,
    vres: int,
    um_per_pix: float,
    geometry: SmartspimGeometryConfig,
) -> tuple[float, float]:
    """Origin shift for overlapping mosaics under mean_center placement.

  When adjacent tile centers in the stage table are much closer than one tile
  width, the arithmetic mean of tile centers is not the stitched volume center.
  Empirically (Multi_RES_SCANs slice-qc sweep) the residual is ~one tile width
  minus stage pitch per axis, scaled by ``(n_tiles - 1)`` along that axis.
    """
    xs = sorted({coord[0] for coord in tile_centers_stage})
    ys = sorted({coord[1] for coord in tile_centers_stage})
    scale = _stage_scale(geometry)
    f0, f1 = geometry.lateral_flip
    tile_w = float(hres) * um_per_pix
    tile_h = float(vres) * um_per_pix

    dx = 0.0
    if len(xs) > 1:
        pitch_x = (xs[1] - xs[0]) * scale
        if pitch_x < tile_w - 1e-3:
            dx = f0 * (tile_w - pitch_x) * (len(xs) - 1)

    dy = 0.0
    if len(ys) > 1:
        pitch_y = (ys[1] - ys[0]) * scale
        if pitch_y < tile_h - 1e-3:
            dy = f1 * (tile_h - pitch_y) * (len(ys) - 1)

    return dx, dy


def _overview_stage_y_anchors(
    overview_tiles: list[tuple[float, float, float]],
) -> tuple[float, float]:
    ys = sorted({coord[1] for coord in overview_tiles})
    if len(ys) < 2:
        msg = "Need at least two distinct overview tile Y stage values for ROI Y anchoring"
        raise ValueError(msg)
    return ys[0], ys[-1]


def _overview_stage_y_anchors_bracketing_roi(
    overview_tiles: list[tuple[float, float, float]],
    roi_tiles: list[tuple[float, float, float]],
) -> tuple[float, float]:
    """Pick overview tile rows that bracket the ROI stage-Y extent (not global min/max)."""
    ov_ys = sorted({coord[1] for coord in overview_tiles})
    if len(ov_ys) < 2:
        msg = "Need at least two distinct overview tile Y stage values for ROI Y anchoring"
        raise ValueError(msg)

    roi_y_min = min(coord[1] for coord in roi_tiles)
    roi_y_max = max(coord[1] for coord in roi_tiles)
    low = [y for y in ov_ys if y <= roi_y_min]
    high = [y for y in ov_ys if y >= roi_y_max]
    y_min_stage = max(low) if low else ov_ys[0]
    y_max_stage = min(high) if high else ov_ys[-1]
    if y_min_stage >= y_max_stage:
        return ov_ys[0], ov_ys[-1]
    return y_min_stage, y_max_stage


def _overview_stage_x_to_pixel_col(
    stage_x: float,
    overview_tiles: list[tuple[float, float, float]],
    *,
    hres: int,
    geometry: SmartspimGeometryConfig,
) -> float:
    """Map a stage-table X coordinate to the overview pixel column of that tile-column center."""
    xs = sorted({coord[0] for coord in overview_tiles})
    if len(xs) < 1:
        msg = "Need at least one overview tile X stage value"
        raise ValueError(msg)

    f0 = geometry.lateral_flip[0]
    if f0 < 0:
        stage_to_col = {x: float(len(xs) - 1 - index) for index, x in enumerate(xs)}
    else:
        stage_to_col = {x: float(index) for index, x in enumerate(xs)}

    if stage_x <= xs[0]:
        col_frac = stage_to_col[xs[0]]
    elif stage_x >= xs[-1]:
        col_frac = stage_to_col[xs[-1]]
    else:
        col_frac = stage_to_col[xs[-1]]
        for index in range(len(xs) - 1):
            x_lo, x_hi = xs[index], xs[index + 1]
            if x_lo <= stage_x <= x_hi:
                col_lo = stage_to_col[x_lo]
                col_hi = stage_to_col[x_hi]
                frac = (stage_x - x_lo) / (x_hi - x_lo)
                col_frac = col_lo + frac * (col_hi - col_lo)
                break

    return col_frac * hres + (hres - 1) / 2.0


def _physical_x_at_overview_stage_col(
    stage_x: float,
    *,
    overview_tiles: list[tuple[float, float, float]],
    shape_yx: tuple[int, int],
    um_per_pix: float,
    hres: int,
    geometry: SmartspimGeometryConfig,
) -> float:
    """Map a stage-table X coordinate into overview physical X (µm)."""
    origin_x, _origin_y = _stage_center_origin_um(
        overview_tiles,
        shape_yx=shape_yx,
        um_per_pix=um_per_pix,
        geometry=geometry,
    )
    f0 = geometry.lateral_flip[0]
    pixel_col = _overview_stage_x_to_pixel_col(
        stage_x,
        overview_tiles,
        hres=hres,
        geometry=geometry,
    )
    return origin_x + f0 * um_per_pix * pixel_col


def _overview_stage_y_to_pixel_row(
    stage_y: float,
    overview_tiles: list[tuple[float, float, float]],
    *,
    vres: int,
    geometry: SmartspimGeometryConfig,
) -> float:
    """Map a stage-table Y coordinate to the overview pixel row of that tile-row center."""
    ys = sorted({coord[1] for coord in overview_tiles})
    if len(ys) < 2:
        msg = "Need at least two distinct overview tile Y stage values"
        raise ValueError(msg)

    f1 = geometry.lateral_flip[1]
    if f1 < 0:
        stage_to_row = {y: float(len(ys) - 1 - index) for index, y in enumerate(ys)}
    else:
        stage_to_row = {y: float(index) for index, y in enumerate(ys)}

    if stage_y <= ys[0]:
        row_frac = stage_to_row[ys[0]]
    elif stage_y >= ys[-1]:
        row_frac = stage_to_row[ys[-1]]
    else:
        row_frac = stage_to_row[ys[-1]]
        for index in range(len(ys) - 1):
            y_lo, y_hi = ys[index], ys[index + 1]
            if y_lo <= stage_y <= y_hi:
                row_lo = stage_to_row[y_lo]
                row_hi = stage_to_row[y_hi]
                frac = (stage_y - y_lo) / (y_hi - y_lo)
                row_frac = row_lo + frac * (row_hi - row_lo)
                break

    return row_frac * vres + (vres - 1) / 2.0


def _physical_y_at_overview_stage_row(
    stage_y: float,
    *,
    overview_tiles: list[tuple[float, float, float]],
    roi_tiles: list[tuple[float, float, float]] | None = None,
    shape_yx: tuple[int, int],
    um_per_pix: float,
    vres: int,
    geometry: SmartspimGeometryConfig,
) -> float:
    """Map a stage-table Y coordinate into overview physical Y (µm)."""
    del roi_tiles  # bracketing anchors are only used by callers for documentation / tests
    _, origin_y = _stage_center_origin_um(
        overview_tiles,
        shape_yx=shape_yx,
        um_per_pix=um_per_pix,
        geometry=geometry,
    )
    f1 = geometry.lateral_flip[1]
    pixel_row = _overview_stage_y_to_pixel_row(
        stage_y,
        overview_tiles,
        vres=vres,
        geometry=geometry,
    )
    return origin_y + f1 * um_per_pix * pixel_row


def _roi_center_x_um(
    roi_tiles: list[tuple[float, float, float]],
    *,
    overview_tiles: list[tuple[float, float, float]] | None,
    overview_shape_yx: tuple[int, int] | None,
    overview_um_per_pix: float | None,
    overview_hres: int | None,
    overview_geometry: SmartspimGeometryConfig | None,
    roi_um_per_pix: float,
    roi_shape_yx: tuple[int, int],
    geometry: SmartspimGeometryConfig,
) -> float:
    xs = [coord[0] for coord in roi_tiles]
    mean_stage_x = sum(xs) / len(xs)
    if (
        geometry.anchor_roi_x_to_overview_stage
        and overview_tiles is not None
        and overview_shape_yx is not None
        and overview_um_per_pix is not None
        and overview_hres is not None
        and overview_geometry is not None
    ):
        return _physical_x_at_overview_stage_col(
            mean_stage_x,
            overview_tiles=overview_tiles,
            shape_yx=overview_shape_yx,
            um_per_pix=overview_um_per_pix,
            hres=overview_hres,
            geometry=overview_geometry,
        )

    scale = _stage_scale(geometry)
    center_x_um = mean_stage_x * scale
    if geometry.stage_xy_is_center:
        return center_x_um
    _ny, nx = roi_shape_yx
    return center_x_um + roi_um_per_pix * (nx - 1) / 2.0


def _roi_center_y_um(
    roi_tiles: list[tuple[float, float, float]],
    *,
    overview_tiles: list[tuple[float, float, float]] | None,
    overview_shape_yx: tuple[int, int] | None,
    overview_um_per_pix: float | None,
    overview_vres: int | None,
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
        and overview_vres is not None
        and overview_geometry is not None
    ):
        return _physical_y_at_overview_stage_row(
            mean_stage_y,
            overview_tiles=overview_tiles,
            roi_tiles=roi_tiles,
            shape_yx=overview_shape_yx,
            um_per_pix=overview_um_per_pix,
            vres=overview_vres,
            geometry=overview_geometry,
        )

    scale = _stage_scale(geometry)
    center_y_um = mean_stage_y * scale
    if geometry.stage_xy_is_center:
        return center_y_um
    ny, _nx = roi_shape_yx
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
    overview_hres: int | None = None,
    overview_vres: int | None = None,
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

    if _resolved_xy_placement(geometry) == "upper_left_tile_corner":
        origin_x, origin_y = _roi_origin_from_upper_left_tile_corner(
            centers,
            hres=meta.hres,
            vres=meta.vres,
            um_per_pix=meta.um_per_pix,
            geometry=geometry,
        )
    else:
        origin_x, origin_y = _stage_center_origin_um(
            centers,
            shape_yx=(ny, nx),
            um_per_pix=meta.um_per_pix,
            geometry=geometry,
        )
        if geometry.apply_mosaic_pitch_correction:
            pitch_dx, _pitch_dy = _mosaic_pitch_correction_um(
                centers,
                hres=meta.hres,
                vres=meta.vres,
                um_per_pix=meta.um_per_pix,
                geometry=geometry,
            )
            origin_x -= pitch_dx
        if label == "roi" and (
            geometry.anchor_roi_x_to_overview_stage or geometry.anchor_roi_y_to_overview_stage
        ):
            if geometry.anchor_roi_x_to_overview_stage:
                center_x_um = _roi_center_x_um(
                    centers,
                    overview_tiles=overview_tiles_for_roi_y,
                    overview_shape_yx=overview_shape_yx,
                    overview_um_per_pix=overview_um_per_pix,
                    overview_hres=overview_hres,
                    overview_geometry=overview_geometry,
                    roi_um_per_pix=meta.um_per_pix,
                    roi_shape_yx=(ny, nx),
                    geometry=geometry,
                )
                f0 = geometry.lateral_flip[0]
                origin_x = center_x_um - f0 * meta.um_per_pix * (nx - 1) / 2.0
            if geometry.anchor_roi_y_to_overview_stage:
                f1 = geometry.lateral_flip[1]
                center_y_um = _roi_center_y_um(
                    centers,
                    overview_tiles=overview_tiles_for_roi_y,
                    overview_shape_yx=overview_shape_yx,
                    overview_um_per_pix=overview_um_per_pix,
                    overview_vres=overview_vres,
                    overview_geometry=overview_geometry,
                    roi_um_per_pix=meta.um_per_pix,
                    roi_shape_yx=(ny, nx),
                    geometry=geometry,
                )
                origin_y = center_y_um - f1 * meta.um_per_pix * (ny - 1) / 2.0

    origin_z = _stage_z_origin_um(
        centers,
        nz=nz,
        z_step_um=meta.z_step_um,
        geometry=geometry,
    )
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
    # Stitched overview placement always uses the full metadata tile table.
    overview_tiles = overview_meta.tile_centers_stage
    roi_tiles = roi_tile_centers_stage
    if roi_tiles is None:
        roi_tiles = parse_smartspim_metadata(roi_meta_path).tile_centers_stage
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
        tile_centers_stage=roi_tiles,
        geometry=roi_geometry,
        overview_tiles_for_roi_y=overview_tiles,
        overview_shape_yx=(overview_ny, overview_nx),
        overview_um_per_pix=overview_meta.um_per_pix,
        overview_hres=overview_meta.hres,
        overview_vres=overview_meta.vres,
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
            "overview_xy_placement": overview_geometry.xy_placement,
            "roi_xy_placement": _resolved_xy_placement(roi_geometry),
            "roi_y_anchored_to_overview_stage": str(roi_geometry.anchor_roi_y_to_overview_stage),
            "roi_x_anchored_to_overview_stage": str(roi_geometry.anchor_roi_x_to_overview_stage),
            "overview_mosaic_pitch_correction": str(overview_geometry.apply_mosaic_pitch_correction),
            "roi_mosaic_pitch_correction": str(roi_geometry.apply_mosaic_pitch_correction),
            "roi_stage_origin_offset_um": ",".join(str(v) for v in roi_geometry.stage_origin_offset_um),
            "z_origin": "stage_table_center" if overview_geometry.stage_z_is_center else "stage_table",
            "stage_z_is_center": str(overview_geometry.stage_z_is_center),
            "overview_tile_count": str(len(overview_tiles)),
            "roi_tile_count": str(len(roi_tiles)),
            "roi_y_anchor_mode": "bracket_roi_stage_y",
            "conversion_timestamp": datetime.now(timezone.utc).isoformat(),
        },
    )

    if output_manifest_path is not None:
        save_pair_manifest(manifest, output_manifest_path)

    return manifest
