"""Unified volume loading and registration-pair preparation for mesoSPIM."""

from __future__ import annotations

from dataclasses import dataclass

import SimpleITK as sitk

from lightsuite.mesospim.config_models import (
    MesospimConfig,
    MesospimGeometryMode,
    MesospimPipelineConfig,
)
from lightsuite.mesospim.geometry import (
    apply_image_geometry,
    apply_voxel_geometry,
    prepare_registration_pair,
)
from lightsuite.mesospim.io import read_tiff_as_float
from lightsuite.mesospim.landmark_geometry import (
    LandmarkFitResult,
    prepare_registration_pair_from_landmarks,
    update_landmark_session_fit,
)
from lightsuite.mesospim.landmark_session import MesospimLandmarkSession
from lightsuite.mesospim.meta import parse_mesospim_meta


@dataclass
class MesospimPreparedPair:
    overview: sitk.Image
    roi: sitk.Image
    fixed_cropped: sitk.Image
    moving: sitk.Image
    overlap_box: tuple
    crop_start_index: list[int]
    landmark_fit: LandmarkFitResult | None = None
    landmark_session: MesospimLandmarkSession | None = None


def load_overview_roi_volumes(cfg: MesospimPipelineConfig) -> tuple[sitk.Image, sitk.Image]:
    meso = cfg.mesospim
    overview_path = meso.overview.path
    roi_path = meso.roi.path
    overview = read_tiff_as_float(
        overview_path,
        overview_path=overview_path,
        roi_path=roi_path,
        remap=meso.tiff_remap,
    )
    roi = read_tiff_as_float(
        roi_path,
        overview_path=overview_path,
        roi_path=roi_path,
        remap=meso.tiff_remap,
    )
    return overview, roi


def apply_volume_geometry(
    overview: sitk.Image,
    roi: sitk.Image,
    meso: MesospimConfig,
) -> None:
    """Apply metadata or voxel geometry in-place on both volumes."""
    mode = meso.geometry_mode
    if mode in (MesospimGeometryMode.METADATA, MesospimGeometryMode.HYBRID):
        overview_meta = parse_mesospim_meta(meso.overview.resolved_meta_path())
        roi_meta = parse_mesospim_meta(meso.roi.resolved_meta_path())
        apply_image_geometry(overview, overview_meta, meso.geometry)
        apply_image_geometry(roi, roi_meta, meso.geometry)
        return

    overview_voxel_um = meso.overview.voxel_um
    roi_voxel_um = meso.roi.voxel_um
    if overview_voxel_um is None or roi_voxel_um is None:
        msg = "geometry_mode=landmarks requires voxel_um on overview and roi"
        raise ValueError(msg)
    apply_voxel_geometry(overview, overview_voxel_um)
    apply_voxel_geometry(roi, roi_voxel_um)


def load_landmark_session(cfg: MesospimPipelineConfig) -> MesospimLandmarkSession:
    path = cfg.mesospim.resolved_landmark_session_path(cfg.sample.save_path)
    if not path.is_file():
        msg = (
            f"Landmark session not found: {path}\n"
            "Create mesospim_landmarks.json or run mesospim match-points (future GUI)."
        )
        raise FileNotFoundError(msg)
    return MesospimLandmarkSession.load(path)


def prepare_mesospim_registration_pair(
    cfg: MesospimPipelineConfig,
    *,
    overview: sitk.Image | None = None,
    roi: sitk.Image | None = None,
    landmark_session: MesospimLandmarkSession | None = None,
) -> MesospimPreparedPair:
    """Load volumes, apply geometry, and build the fixed/moving pair for elastix."""
    meso = cfg.mesospim
    if overview is None or roi is None:
        overview, roi = load_overview_roi_volumes(cfg)
    apply_volume_geometry(overview, roi, meso)

    margin_um = meso.registration.overlap_margin_um
    mode = meso.geometry_mode

    if mode == MesospimGeometryMode.METADATA:
        fixed_cropped, moving, overlap_box, crop_start_index = prepare_registration_pair(
            overview,
            roi,
            margin_um=margin_um,
        )
        return MesospimPreparedPair(
            overview=overview,
            roi=roi,
            fixed_cropped=fixed_cropped,
            moving=moving,
            overlap_box=overlap_box,
            crop_start_index=crop_start_index,
        )

    session = landmark_session or load_landmark_session(cfg)
    fit_mode = meso.landmarks.fit_mode
    min_pairs = meso.landmarks.min_pairs

    fixed_cropped, moving, overlap_box, crop_start_index, fit = (
        prepare_registration_pair_from_landmarks(
            overview,
            roi,
            session=session,
            fit_mode=fit_mode,
            min_pairs=min_pairs,
            margin_um=margin_um,
        )
    )
    update_landmark_session_fit(session, fit)

    session_path = cfg.mesospim.resolved_landmark_session_path(cfg.sample.save_path)
    session.save(session_path)

    return MesospimPreparedPair(
        overview=overview,
        roi=roi,
        fixed_cropped=fixed_cropped,
        moving=moving,
        overlap_box=overlap_box,
        crop_start_index=crop_start_index,
        landmark_fit=fit,
        landmark_session=session,
    )


def resolve_meta_paths(cfg: MesospimPipelineConfig) -> tuple[str | None, str | None]:
    """Return overview/roi meta paths when available."""
    overview_meta = (
        str(m.resolved_meta_path())
        if (m := cfg.mesospim.overview).has_meta_sidecar()
        else None
    )
    roi_meta = str(m.resolved_meta_path()) if (m := cfg.mesospim.roi).has_meta_sidecar() else None
    return overview_meta, roi_meta
