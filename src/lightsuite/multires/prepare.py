"""Registration-pair preparation from multiresolution pair manifests."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import SimpleITK as sitk

from lightsuite.multires.config_models import MultiresGeometryMode, MultiresPipelineConfig
from lightsuite.multires.landmark_session import MultiresLandmarkSession
from lightsuite.multires.landmarks import (
    LandmarkFitResult,
    fit_landmark_transform,
    update_landmark_session_fit,
)
from lightsuite.multires.memory import warn_if_overlap_memory_exceeds_system
from lightsuite.multires.models import ManifestVolumeSpec, MultiresPairManifest
from lightsuite.multires.resolve import resolve_pair_manifest
from lightsuite.multires.spec_geometry import (
    crop_index_range_from_physical_box,
    overlap_box_from_landmark_specs,
    overlap_physical_bounds_from_specs,
    sitk_geometry_from_spec,
)
from lightsuite.multires.volume import load_manifest_xyz_crop, stream_resample_to_reference
from lightsuite.reporter import check_stage_cancelled


@dataclass
class MultiresPreparedPair:
    """Overlap crops ready for elastix — full stacks are never held in memory."""

    fixed_cropped: sitk.Image
    moving: sitk.Image
    overlap_box: tuple[np.ndarray, np.ndarray]
    crop_start_index: list[int]
    overview_spec: ManifestVolumeSpec
    landmark_fit: LandmarkFitResult | None = None
    landmark_session: MultiresLandmarkSession | None = None


def load_landmark_session(path: Path) -> MultiresLandmarkSession:
    if not path.is_file():
        msg = (
            f"Landmark session not found: {path}\n"
            "Create multires_landmarks.json or run multires match-points."
        )
        raise FileNotFoundError(msg)
    return MultiresLandmarkSession.load(path)


def prepare_multires_registration_pair(
    cfg: MultiresPipelineConfig,
    *,
    manifest: MultiresPairManifest | None = None,
    landmark_session: MultiresLandmarkSession | None = None,
    overview_spec: ManifestVolumeSpec | None = None,
    roi_spec: ManifestVolumeSpec | None = None,
    progress_prefix: str | None = None,
) -> MultiresPreparedPair:
    """Build fixed/moving overlap crops by streaming planes from disk.

    Neither the overview nor the ROI full stack is materialised. Only the
    shared physical overlap (plus a small Z chunk of the ROI at a time while
    resampling) is held in RAM.
    """
    if manifest is None:
        manifest, manifest_path = resolve_pair_manifest(cfg)
    else:
        manifest_path = cfg.multires.resolved_pair_manifest_path(
            cfg.sample.save_path,
            cfg.sample.name,
        )
    manifest_dir = manifest_path.parent
    margin_um = cfg.multires.registration.overlap_margin_um
    mode = cfg.multires.geometry_mode
    overview_spec = overview_spec or manifest.overview
    roi_spec = roi_spec or manifest.roi

    landmark_fit: LandmarkFitResult | None = None
    session: MultiresLandmarkSession | None = None
    reference_to_moving: np.ndarray | None = None

    if mode == MultiresGeometryMode.METADATA:
        overlap_min, overlap_max = overlap_physical_bounds_from_specs(
            overview_spec,
            roi_spec,
            margin_um=margin_um,
        )
    else:
        session_path = cfg.multires.resolved_landmark_session_path(cfg.sample.save_path, manifest)
        session = landmark_session or load_landmark_session(session_path)
        # Geometry-only 1-voxel images carry spacing/origin/direction for point maps.
        overview_geo = sitk_geometry_from_spec(overview_spec)
        roi_geo = sitk_geometry_from_spec(roi_spec)
        landmark_fit = fit_landmark_transform(
            overview=overview_geo,
            roi=roi_geo,
            session=session,
            fit_mode=cfg.multires.landmarks.fit_mode,
            min_pairs=cfg.multires.landmarks.min_pairs,
        )
        overlap_min, overlap_max = overlap_box_from_landmark_specs(
            overview_spec,
            roi_spec,
            landmark_fit.roi_to_overview_tform,
            margin_um=margin_um,
        )
        reference_to_moving = np.linalg.inv(landmark_fit.roi_to_overview_tform)
        update_landmark_session_fit(session, landmark_fit)
        session.save(session_path)

    crop_start_index, crop_size = crop_index_range_from_physical_box(
        overview_spec,
        overlap_min,
        overlap_max,
    )
    warn_if_overlap_memory_exceeds_system(
        crop_size,
        max_slab_bytes=cfg.multires.registration.max_slab_bytes,
    )
    check_stage_cancelled()
    overview_label = f"{progress_prefix} overview crop" if progress_prefix else None
    roi_label = f"{progress_prefix} ROI resample" if progress_prefix else None
    fixed_cropped = load_manifest_xyz_crop(
        overview_spec,
        start_xyz=crop_start_index,
        crop_size_xyz=crop_size,
        manifest_dir=manifest_dir,
        progress_label=overview_label,
    )
    check_stage_cancelled()
    moving = stream_resample_to_reference(
        roi_spec,
        fixed_cropped,
        manifest_dir=manifest_dir,
        reference_to_moving=reference_to_moving,
        max_slab_bytes=cfg.multires.registration.max_slab_bytes,
        progress_label=roi_label,
    )

    return MultiresPreparedPair(
        fixed_cropped=fixed_cropped,
        moving=moving,
        overlap_box=(overlap_min, overlap_max),
        crop_start_index=crop_start_index,
        overview_spec=overview_spec,
        landmark_fit=landmark_fit,
        landmark_session=session,
    )
