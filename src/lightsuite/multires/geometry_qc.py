"""Physical-space overlap QC for multiresolution mesoSPIM geometry."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import numpy as np
import SimpleITK as sitk

from lightsuite.mesospim.config_models import MesospimGeometryConfig
from lightsuite.mesospim.meta import parse_mesospim_meta
from lightsuite.multires.config_models import MultiresPipelineConfig
from lightsuite.multires.geometry import resample_to_reference_grid
from lightsuite.multires.models import ManifestVolumeSpec
from lightsuite.multires.plots import (
    _normalize_panel,
    _resample_to_shape,
    normalized_cross_correlation,
)
from lightsuite.multires.resolve import _merge_mesospim_geometry, resolve_pair_manifest
from lightsuite.multires.spec_geometry import (
    crop_index_range_from_physical_box,
    overlap_physical_bounds_from_specs,
    physical_to_continuous_index_xyz,
)
from lightsuite.multires.vendor.mesospim import (
    _overview_volume_spec,
    _resolve_overview_meta_path,
    _resolve_roi_meta_path,
    _roi_volume_spec,
)


@dataclass(frozen=True)
class GeometryQcSlice:
    """One overlap slice with ROI physically resampled onto the overview grid."""

    overview_display: np.ndarray
    roi_resampled_display: np.ndarray
    physical_ncc: float
    overview_z: int
    roi_z: int
    lateral_flip: tuple[int, int]
    has_overlap: bool


def lateral_flip_from_config(cfg: MultiresPipelineConfig) -> tuple[int, int]:
    """Return the configured lateral_flip, preferring overview then ROI overrides."""
    meso_geom = cfg.multires.mesospim_geometry
    if meso_geom is not None:
        if meso_geom.overview is not None and meso_geom.overview.lateral_flip is not None:
            fx, fy = meso_geom.overview.lateral_flip
            return int(fx), int(fy)
        if meso_geom.roi is not None and meso_geom.roi.lateral_flip is not None:
            fx, fy = meso_geom.roi.lateral_flip
            return int(fx), int(fy)
    return 1, -1


def lateral_flip_tuple(*, flip_x: bool, flip_y: bool) -> tuple[int, int]:
    return (-1 if flip_x else 1, -1 if flip_y else 1)


def mesospim_geometry_yaml_snippet(lateral_flip: tuple[int, int]) -> str:
    fx, fy = lateral_flip
    return (
        "  mesospim_geometry:\n"
        "    overview:\n"
        f"      lateral_flip: [{fx}, {fy}]\n"
        "    roi:\n"
        f"      lateral_flip: [{fx}, {fy}]"
    )


def _geometry_for_lateral_flip(
    cfg: MultiresPipelineConfig,
    lateral_flip: tuple[int, int],
) -> MesospimGeometryConfig:
    base = _merge_mesospim_geometry(
        cfg.multires.mesospim_geometry.overview if cfg.multires.mesospim_geometry else None
    ) or MesospimGeometryConfig()
    return base.model_copy(update={"lateral_flip": [int(lateral_flip[0]), int(lateral_flip[1])]})


def _reference_volume_paths(
    cfg: MultiresPipelineConfig,
    *,
    channel: str | None = None,
) -> tuple[Path, Path, Path | None, Path | None]:
    meso = cfg.multires
    if meso.channels:
        ref = channel or meso.registration.reference_channel
        if ref is None:
            msg = "registration.reference_channel is required"
            raise ValueError(msg)
        if ref not in meso.channels:
            msg = (
                f"channel {ref!r} missing from multires.channels "
                f"(available: {sorted(meso.channels)})"
            )
            raise ValueError(msg)
        channel_paths = meso.channels[ref]
        return (
            channel_paths.overview,
            channel_paths.roi,
            channel_paths.overview_meta_path or meso.overview_meta_path,
            channel_paths.roi_meta_path,
        )

    manifest, _manifest_path = resolve_pair_manifest(cfg, rebuild=False)
    prov = manifest.provenance
    overview = prov.get("source_overview")
    roi = prov.get("source_roi")
    if not overview or not roi:
        msg = (
            "inspect-geometry needs multires.channels in the YAML config, or a pair manifest "
            "built from mesoSPIM exports (with source_overview / source_roi provenance)."
        )
        raise ValueError(msg)
    overview_meta = prov.get("overview_meta")
    roi_meta = prov.get("roi_meta")
    return (
        Path(overview).expanduser().resolve(),
        Path(roi).expanduser().resolve(),
        Path(overview_meta).expanduser().resolve() if overview_meta else None,
        Path(roi_meta).expanduser().resolve() if roi_meta else None,
    )


def build_reference_specs_with_geometry(
    cfg: MultiresPipelineConfig,
    lateral_flip: tuple[int, int],
    *,
    channel: str | None = None,
) -> tuple[ManifestVolumeSpec, ManifestVolumeSpec, Path]:
    """Rebuild overview / ROI manifest specs for one shared lateral_flip."""
    overview_path, roi_path, overview_meta_path, roi_meta_path = _reference_volume_paths(
        cfg,
        channel=channel,
    )
    geometry = _geometry_for_lateral_flip(cfg, lateral_flip)

    overview_meta_path = _resolve_overview_meta_path(overview_path, overview_meta_path)
    roi_meta_path = _resolve_roi_meta_path(roi_path, roi_meta_path)
    overview_spec = _overview_volume_spec(
        overview_path,
        parse_mesospim_meta(overview_meta_path),
        geometry,
    )
    roi_spec = _roi_volume_spec(roi_path, parse_mesospim_meta(roi_meta_path), geometry)
    manifest_dir = cfg.multires.resolved_pair_manifest_path(
        cfg.sample.save_path,
        cfg.sample.name,
    ).parent
    return overview_spec, roi_spec, manifest_dir


def _default_z_at_overlap_center(
    overview_spec: ManifestVolumeSpec,
    roi_spec: ManifestVolumeSpec,
) -> tuple[int, int, np.ndarray]:
    overlap_min, overlap_max = overlap_physical_bounds_from_specs(overview_spec, roi_spec)
    center = 0.5 * (overlap_min + overlap_max)
    overview_z = int(round(physical_to_continuous_index_xyz(overview_spec, center)[2]))
    roi_z = int(round(physical_to_continuous_index_xyz(roi_spec, center)[2]))
    nz_ov = int(overview_spec.shape_zyx[0])
    nz_roi = int(roi_spec.shape_zyx[0])
    overview_z = int(np.clip(overview_z, 0, max(nz_ov - 1, 0)))
    roi_z = int(np.clip(roi_z, 0, max(nz_roi - 1, 0)))
    return overview_z, roi_z, center


def compute_geometry_qc_slice(
    overview_spec: ManifestVolumeSpec,
    roi_spec: ManifestVolumeSpec,
    *,
    manifest_dir: Path,
    overview_z: int | None = None,
    roi_z: int | None = None,
    lateral_flip: tuple[int, int] = (1, -1),
    link_z: bool = True,
) -> GeometryQcSlice:
    """Load one overlap slice; resample ROI onto overview in physical space.

    When ``link_z`` is False, each panel shows its native Z slice cropped to the
    shared XY overlap (2D resize only). Physical NCC is not meaningful in that
    mode and is returned as NaN.
    """
    from lightsuite.multires.volume import load_manifest_xyz_crop

    try:
        overlap_min, overlap_max = overlap_physical_bounds_from_specs(overview_spec, roi_spec)
    except ValueError:
        blank = np.zeros((8, 8), dtype=np.float32)
        return GeometryQcSlice(
            overview_display=blank,
            roi_resampled_display=blank,
            physical_ncc=0.0,
            overview_z=0,
            roi_z=0,
            lateral_flip=lateral_flip,
            has_overlap=False,
        )

    center = 0.5 * (overlap_min + overlap_max)
    if overview_z is None or roi_z is None:
        default_ov_z, default_roi_z, _center = _default_z_at_overlap_center(overview_spec, roi_spec)
        overview_z = default_ov_z if overview_z is None else overview_z
        roi_z = default_roi_z if roi_z is None else roi_z

    overview_start, overview_size = crop_index_range_from_physical_box(
        overview_spec,
        overlap_min,
        overlap_max,
    )
    roi_start, roi_size = crop_index_range_from_physical_box(
        roi_spec,
        overlap_min,
        overlap_max,
    )

    overview_ref = load_manifest_xyz_crop(
        overview_spec,
        start_xyz=[overview_start[0], overview_start[1], int(overview_z)],
        crop_size_xyz=[overview_size[0], overview_size[1], 1],
        manifest_dir=manifest_dir,
    )
    roi_mov = load_manifest_xyz_crop(
        roi_spec,
        start_xyz=[roi_start[0], roi_start[1], int(roi_z)],
        crop_size_xyz=[roi_size[0], roi_size[1], 1],
        manifest_dir=manifest_dir,
    )
    sl_overview = np.asarray(sitk.GetArrayFromImage(overview_ref)[0], dtype=np.float32)
    if link_z:
        resampled = resample_to_reference_grid(roi_mov, overview_ref)
        sl_roi = np.asarray(sitk.GetArrayFromImage(resampled)[0], dtype=np.float32)
    else:
        sl_roi = np.asarray(sitk.GetArrayFromImage(roi_mov)[0], dtype=np.float32)
        if sl_roi.shape != sl_overview.shape:
            sl_roi = _resample_to_shape(sl_roi, sl_overview.shape)
    overview_norm = _normalize_panel(sl_overview)
    roi_norm = _normalize_panel(sl_roi)
    ncc = (
        normalized_cross_correlation(overview_norm, roi_norm)
        if link_z
        else float("nan")
    )

    return GeometryQcSlice(
        overview_display=overview_norm,
        roi_resampled_display=roi_norm,
        physical_ncc=ncc,
        overview_z=int(overview_z),
        roi_z=int(roi_z),
        lateral_flip=lateral_flip,
        has_overlap=True,
    )


def evaluate_geometry_qc(
    cfg: MultiresPipelineConfig,
    *,
    lateral_flip: tuple[int, int] | None = None,
    overview_z: int | None = None,
    roi_z: int | None = None,
    link_z: bool = True,
) -> GeometryQcSlice:
    """Convenience wrapper: rebuild specs from config and score one overlap slice."""
    flip = lateral_flip or lateral_flip_from_config(cfg)
    overview_spec, roi_spec, manifest_dir = build_reference_specs_with_geometry(cfg, flip)
    return compute_geometry_qc_slice(
        overview_spec,
        roi_spec,
        manifest_dir=manifest_dir,
        overview_z=overview_z,
        roi_z=roi_z,
        lateral_flip=flip,
        link_z=link_z,
    )


__all__ = [
    "GeometryQcSlice",
    "build_reference_specs_with_geometry",
    "compute_geometry_qc_slice",
    "evaluate_geometry_qc",
    "lateral_flip_from_config",
    "lateral_flip_tuple",
    "mesospim_geometry_yaml_snippet",
]
