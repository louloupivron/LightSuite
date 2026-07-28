"""Build multiresolution pair manifests from SmartSPIM / ASI exports.

Stage-table conventions (see ``examples/config/smartspim/context.txt`` for the
supporting measurements):

* ``X`` and ``Y`` are tile FOV centers in units of 0.1 µm.
* ``Z`` is the *first* plane of the stack, already in µm; the stack spans
  ``NumImages * Z step`` from there.
* Adjacent tiles overlap by a hard-coded 10%, and the stitcher anchors the
  mosaic on the upper-left tile, so the stitched pixel ``[0, 0]`` sits at that
  tile's leading corner.

Metadata may be either legacy tab-separated ``metadata.txt`` or JSON
``metadata.json`` (``sample_metadata`` + ``tiles``).
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from lightsuite.multires.manifest import save_pair_manifest
from lightsuite.multires.models import MANIFEST_FORMAT, ManifestVolumeSpec, MultiresPairManifest
from lightsuite.multires.volume import discover_volume_shape

# ASI stage X/Y table values are tenths of a micron (0.1 µm).
STAGE_COORD_SCALE_UM = 0.1

# Acquisition software hard-codes the tile overlap; used to validate the stage scale.
NOMINAL_TILE_OVERLAP = 0.10


@dataclass(frozen=True)
class SmartspimGeometryConfig:
    """How ASI stage coordinates map into a physical (µm) frame."""

    stage_coord_scale_um: float = STAGE_COORD_SCALE_UM
    stage_xy_is_center: bool = True
    stage_z_is_center: bool = False
    lateral_flip: tuple[int, int] = (1, 1)


@dataclass(frozen=True)
class SmartspimScanMeta:
    """Parsed fields from a SmartSPIM ``metadata.txt`` / ``metadata.json`` export."""

    objective: str
    hres: int
    vres: int
    um_per_pix: float
    z_step_um: float
    tile_centers_stage: list[tuple[float, float, float]]
    tile_num_images: list[int]

    @property
    def tile_width_um(self) -> float:
        return float(self.hres) * self.um_per_pix

    @property
    def tile_height_um(self) -> float:
        return float(self.vres) * self.um_per_pix


def _finalize_scan_meta(
    *,
    path: Path,
    objective: str,
    hres: int,
    vres: int,
    um_per_pix: float,
    z_step_um: float,
    tile_centers_stage: list[tuple[float, float, float]],
    tile_num_images: list[int],
) -> SmartspimScanMeta:
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


def _parse_smartspim_metadata_txt(path: Path) -> SmartspimScanMeta:
    """Parse ASI SmartSPIM ``metadata.txt`` (tab-separated export)."""
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    if len(lines) < 2:
        msg = f"SmartSPIM metadata too short: {path}"
        raise ValueError(msg)

    header = lines[0].split("\t")
    values = lines[1].split("\t")
    field_map = {
        key.strip(): values[index].strip() for index, key in enumerate(header) if key.strip()
    }

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

    return _finalize_scan_meta(
        path=path,
        objective=objective,
        hres=hres,
        vres=vres,
        um_per_pix=um_per_pix,
        z_step_um=z_step_um,
        tile_centers_stage=tile_centers_stage,
        tile_num_images=tile_num_images,
    )


def _parse_smartspim_metadata_json(path: Path) -> SmartspimScanMeta:
    """Parse ASI SmartSPIM ``metadata.json`` (``sample_metadata`` + ``tiles``)."""
    try:
        payload: dict[str, Any] = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        msg = f"Invalid SmartSPIM metadata JSON: {path}"
        raise ValueError(msg) from exc

    sample = payload.get("sample_metadata")
    tiles = payload.get("tiles")
    if not isinstance(sample, dict):
        msg = f"Missing sample_metadata object in {path}"
        raise ValueError(msg)
    if not isinstance(tiles, list) or not tiles:
        msg = f"Missing tiles array in {path}"
        raise ValueError(msg)

    try:
        objective = str(sample["objective"])
        hres = int(float(sample["horizontal_resolution"]))
        vres = int(float(sample["vertical_resolution"]))
        um_per_pix = float(sample["um_per_pix"])
        z_step_um = float(sample["z_step_um"])
    except KeyError as exc:
        msg = f"Missing SmartSPIM sample_metadata field {exc!s} in {path}"
        raise ValueError(msg) from exc
    except (TypeError, ValueError) as exc:
        msg = f"Invalid SmartSPIM sample_metadata numeric field in {path}"
        raise ValueError(msg) from exc

    tile_centers_stage: list[tuple[float, float, float]] = []
    tile_num_images: list[int] = []
    for index, tile in enumerate(tiles):
        if not isinstance(tile, dict):
            msg = f"Tile {index} is not an object in {path}"
            raise ValueError(msg)
        try:
            tile_centers_stage.append((float(tile["X"]), float(tile["Y"]), float(tile["Z"])))
            tile_num_images.append(int(float(tile["NumImages"])))
        except KeyError as exc:
            msg = f"Missing tile field {exc!s} at index {index} in {path}"
            raise ValueError(msg) from exc
        except (TypeError, ValueError) as exc:
            msg = f"Invalid tile numeric field at index {index} in {path}"
            raise ValueError(msg) from exc

    return _finalize_scan_meta(
        path=path,
        objective=objective,
        hres=hres,
        vres=vres,
        um_per_pix=um_per_pix,
        z_step_um=z_step_um,
        tile_centers_stage=tile_centers_stage,
        tile_num_images=tile_num_images,
    )


def parse_smartspim_metadata(path: Path) -> SmartspimScanMeta:
    """Parse ASI SmartSPIM ``metadata.txt`` or ``metadata.json``."""
    path = path.expanduser().resolve()
    suffix = path.suffix.lower()
    if suffix == ".json":
        return _parse_smartspim_metadata_json(path)
    if suffix in {".txt", ".tsv", ""}:
        return _parse_smartspim_metadata_txt(path)
    # Fall back by content for oddly named exports.
    head = path.read_text(encoding="utf-8", errors="replace").lstrip()[:1]
    if head == "{":
        return _parse_smartspim_metadata_json(path)
    return _parse_smartspim_metadata_txt(path)


def stage_pitch_overlap_fractions(
    tile_centers_stage: list[tuple[float, float, float]],
    *,
    meta: SmartspimScanMeta,
    geometry: SmartspimGeometryConfig,
) -> dict[str, float]:
    """Tile overlap fraction implied by the stage pitch, per mosaic axis.

    Axes with a single tile are omitted.
    """
    scale = float(geometry.stage_coord_scale_um)
    fractions: dict[str, float] = {}
    for axis, index, tile_size in (("x", 0, meta.tile_width_um), ("y", 1, meta.tile_height_um)):
        values = sorted({coord[index] for coord in tile_centers_stage})
        if len(values) < 2:
            continue
        pitch_um = (values[1] - values[0]) * scale
        fractions[axis] = 1.0 - pitch_um / tile_size
    return fractions


def check_stage_scale(
    tile_centers_stage: list[tuple[float, float, float]],
    *,
    meta: SmartspimScanMeta,
    geometry: SmartspimGeometryConfig,
    tolerance: float = 0.02,
) -> None:
    """Fail fast when ``stage_coord_scale_um`` disagrees with the hard-coded overlap.

    A wrong scale silently places volumes millimetres apart, which is not
    recoverable downstream, so mosaics are checked against the known 10%
    overlap before a manifest is written.
    """
    fractions = stage_pitch_overlap_fractions(tile_centers_stage, meta=meta, geometry=geometry)
    bad = {
        axis: value
        for axis, value in fractions.items()
        if abs(value - NOMINAL_TILE_OVERLAP) > tolerance
    }
    if not bad:
        return
    detail = ", ".join(f"{axis}={100 * value:.1f}%" for axis, value in sorted(bad.items()))
    msg = (
        f"Stage pitch implies a tile overlap of {detail}, but the acquisition software "
        f"hard-codes {100 * NOMINAL_TILE_OVERLAP:.0f}%. "
        f"stage_coord_scale_um={geometry.stage_coord_scale_um} is probably wrong "
        f"(tile FOV {meta.tile_width_um:.1f} x {meta.tile_height_um:.1f} µm)."
    )
    raise ValueError(msg)


def _resolve_num_images(tile_num_images: list[int]) -> int:
    """Return the metadata NumImages value used for Z extent."""
    unique = sorted(set(tile_num_images))
    if len(unique) != 1:
        msg = f"Expected a single NumImages value across tiles, got {unique}"
        raise ValueError(msg)
    return unique[0]


def _direction_from_geometry(geometry: SmartspimGeometryConfig) -> list[float]:
    f0, f1 = geometry.lateral_flip
    return [float(f0), 0.0, 0.0, 0.0, float(f1), 0.0, 0.0, 0.0, 1.0]


def _stage_xy_origin_um(
    tile_centers_stage: list[tuple[float, float, float]],
    *,
    meta: SmartspimScanMeta,
    geometry: SmartspimGeometryConfig,
) -> tuple[float, float]:
    """Physical XY of stitched pixel ``[0, 0]``, anchored on the upper-left tile."""
    scale = float(geometry.stage_coord_scale_um)
    xs = [coord[0] for coord in tile_centers_stage]
    ys = [coord[1] for coord in tile_centers_stage]
    f0, f1 = geometry.lateral_flip

    origin_x = (min(xs) if f0 > 0 else max(xs)) * scale
    origin_y = (min(ys) if f1 > 0 else max(ys)) * scale
    if geometry.stage_xy_is_center:
        origin_x -= f0 * meta.tile_width_um / 2.0
        origin_y -= f1 * meta.tile_height_um / 2.0
    return origin_x, origin_y


def _stage_z_origin_um(
    tile_centers_stage: list[tuple[float, float, float]],
    *,
    nz: int,
    z_step_um: float,
    geometry: SmartspimGeometryConfig,
) -> float:
    """Physical Z of the first plane, from the tile-table Z column (already µm)."""
    zs = [coord[2] for coord in tile_centers_stage]
    z_mean = float(sum(zs) / len(zs))
    if geometry.stage_z_is_center:
        return z_mean - (nz - 1) / 2.0 * z_step_um
    return z_mean


def volume_spec_from_smartspim_export(
    *,
    volume_path: Path,
    metadata_path: Path,
    geometry: SmartspimGeometryConfig | None = None,
) -> ManifestVolumeSpec:
    """Build a manifest volume spec from a stitched SmartSPIM stack folder."""
    volume_path = volume_path.expanduser().resolve()
    geometry = geometry or SmartspimGeometryConfig()
    meta = parse_smartspim_metadata(metadata_path)
    centers = meta.tile_centers_stage
    check_stage_scale(centers, meta=meta, geometry=geometry)

    nz, ny, nx = discover_volume_shape(volume_path)
    num_images = _resolve_num_images(meta.tile_num_images)
    if nz != num_images:
        msg = (
            f"Volume plane count ({nz}) does not match metadata NumImages ({num_images}) "
            f"for {volume_path}"
        )
        raise ValueError(msg)

    origin_x, origin_y = _stage_xy_origin_um(centers, meta=meta, geometry=geometry)
    origin_z = _stage_z_origin_um(centers, nz=nz, z_step_um=meta.z_step_um, geometry=geometry)
    return ManifestVolumeSpec(
        volume_path=str(volume_path),
        shape_zyx=[nz, ny, nx],
        spacing_um=[meta.um_per_pix, meta.um_per_pix, meta.z_step_um],
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
    geometry: SmartspimGeometryConfig | None = None,
    output_manifest_path: Path | None = None,
) -> MultiresPairManifest:
    """Convert stitched SmartSPIM exports into a LightSuite multires pair manifest.

    Overview and ROI share one geometry config: both are placed in the same
    stage frame, so no overview↔ROI anchoring is required.
    """
    geometry = geometry or SmartspimGeometryConfig()
    overview_spec = volume_spec_from_smartspim_export(
        volume_path=overview_path,
        metadata_path=overview_meta_path,
        geometry=geometry,
    )
    roi_spec = volume_spec_from_smartspim_export(
        volume_path=roi_path,
        metadata_path=roi_meta_path,
        geometry=geometry,
    )

    manifest = MultiresPairManifest(
        format=MANIFEST_FORMAT,
        sample_name=sample_name,
        pair_label=pair_label,
        overview=overview_spec,
        roi=roi_spec,
        provenance={
            "microscope": "smartspim",
            "source_overview": str(Path(overview_path).expanduser().resolve()),
            "source_roi": str(Path(roi_path).expanduser().resolve()),
            "overview_meta": str(Path(overview_meta_path).expanduser().resolve()),
            "roi_meta": str(Path(roi_meta_path).expanduser().resolve()),
            "conversion": "lightsuite.multires.vendor.smartspim",
            "stage_coord_scale_um": str(geometry.stage_coord_scale_um),
            "stage_xy_is_center": str(geometry.stage_xy_is_center),
            "stage_z_is_center": str(geometry.stage_z_is_center),
            "lateral_flip": ",".join(str(v) for v in geometry.lateral_flip),
            "xy_origin": "upper_left_tile_corner",
            "conversion_timestamp": datetime.now(UTC).isoformat(),
        },
    )

    if output_manifest_path is not None:
        save_pair_manifest(manifest, output_manifest_path)

    return manifest
