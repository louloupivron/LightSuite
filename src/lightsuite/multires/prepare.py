"""Registration-pair preparation from multiresolution pair manifests."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import SimpleITK as sitk

from lightsuite.multires.config_models import MultiresGeometryMode, MultiresPipelineConfig
from lightsuite.multires.geometry import prepare_registration_pair
from lightsuite.multires.landmark_session import MultiresLandmarkSession
from lightsuite.multires.landmarks import (
    LandmarkFitResult,
    prepare_registration_pair_from_landmarks,
    update_landmark_session_fit,
)
from lightsuite.multires.manifest import load_pair_manifest
from lightsuite.multires.models import MultiresPairManifest
from lightsuite.multires.volume import load_manifest_volume


@dataclass
class MultiresPreparedPair:
    overview: sitk.Image
    roi: sitk.Image
    fixed_cropped: sitk.Image
    moving: sitk.Image
    overlap_box: tuple
    crop_start_index: list[int]
    landmark_fit: LandmarkFitResult | None = None
    landmark_session: MultiresLandmarkSession | None = None


def load_pair_volumes(
    manifest: MultiresPairManifest,
    *,
    manifest_path: Path,
) -> tuple[sitk.Image, sitk.Image]:
    manifest_dir = manifest_path.parent
    overview = load_manifest_volume(manifest.overview, manifest_dir=manifest_dir)
    roi = load_manifest_volume(manifest.roi, manifest_dir=manifest_dir)
    return overview, roi


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
    overview: sitk.Image | None = None,
    roi: sitk.Image | None = None,
    landmark_session: MultiresLandmarkSession | None = None,
) -> MultiresPreparedPair:
    """Load manifest volumes and build the fixed/moving pair for elastix."""
    manifest_path = cfg.multires.pair_manifest
    if manifest is None:
        manifest = load_pair_manifest(manifest_path)
    if overview is None or roi is None:
        overview, roi = load_pair_volumes(manifest, manifest_path=manifest_path)

    margin_um = cfg.multires.registration.overlap_margin_um
    mode = cfg.multires.geometry_mode

    if mode == MultiresGeometryMode.METADATA:
        fixed_cropped, moving, overlap_box, crop_start_index = prepare_registration_pair(
            overview,
            roi,
            margin_um=margin_um,
        )
        return MultiresPreparedPair(
            overview=overview,
            roi=roi,
            fixed_cropped=fixed_cropped,
            moving=moving,
            overlap_box=overlap_box,
            crop_start_index=crop_start_index,
        )

    session_path = cfg.multires.resolved_landmark_session_path(cfg.sample.save_path, manifest)
    session = landmark_session or load_landmark_session(session_path)
    fit_mode = cfg.multires.landmarks.fit_mode
    min_pairs = cfg.multires.landmarks.min_pairs

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
    session.save(session_path)

    return MultiresPreparedPair(
        overview=overview,
        roi=roi,
        fixed_cropped=fixed_cropped,
        moving=moving,
        overlap_box=overlap_box,
        crop_start_index=crop_start_index,
        landmark_fit=fit,
        landmark_session=session,
    )
