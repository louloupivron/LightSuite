"""Landmark-driven coarse placement of ROI within overview."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import SimpleITK as sitk

from lightsuite.gui.affine import (
    affine_point_errors,
    fit_affine_transform,
    fit_similarity_transform,
    summarize_point_errors,
    transform_points,
)
from lightsuite.multires.geometry import (
    crop_to_physical_box,
    physical_bounds,
    transformed_bounds_in_target_space,
    zyx_points_to_physical,
)
from lightsuite.multires.landmark_session import LandmarkFitMode, MultiresLandmarkSession


def fit_rigid_transform(source: np.ndarray, target: np.ndarray) -> tuple[np.ndarray, float]:
    """Fit rotation + translation without scale (source -> target)."""
    if source.shape != target.shape:
        msg = f"Point arrays must match shape, got {source.shape} vs {target.shape}"
        raise ValueError(msg)
    if source.shape[0] < 3:
        msg = "Need at least 3 point pairs for rigid fit"
        raise ValueError(msg)

    centroid_source = source.mean(axis=0)
    centroid_target = target.mean(axis=0)
    source_centered = source - centroid_source
    target_centered = target - centroid_target

    h = source_centered.T @ target_centered
    u, _, vt = np.linalg.svd(h)
    rotation = vt.T @ u.T
    if np.linalg.det(rotation) < 0:
        vt[-1, :] *= -1
        rotation = vt.T @ u.T

    translation = centroid_target - rotation @ centroid_source
    matrix = np.eye(4)
    matrix[:3, :3] = rotation
    matrix[:3, 3] = translation
    predicted = transform_points(source, matrix)
    mse = float(np.mean(np.sum((predicted - target) ** 2, axis=1)))
    return matrix, mse


def fit_roi_to_overview_transform(
    roi_points_um: np.ndarray,
    overview_points_um: np.ndarray,
    *,
    fit_mode: LandmarkFitMode,
) -> tuple[np.ndarray, float]:
    """Fit a 4x4 transform mapping ROI physical points to overview physical points."""
    if fit_mode == "similarity":
        return fit_similarity_transform(roi_points_um, overview_points_um)
    if fit_mode == "affine":
        return fit_affine_transform(roi_points_um, overview_points_um)
    if fit_mode == "rigid":
        return fit_rigid_transform(roi_points_um, overview_points_um)
    msg = f"Unsupported fit_mode: {fit_mode!r}"
    raise ValueError(msg)


@dataclass
class LandmarkFitResult:
    roi_to_overview_tform: np.ndarray
    rms_error_um: float
    fit_point_errors_um: np.ndarray
    fit_stats: dict[str, float]


def fit_landmark_transform(
    *,
    overview: sitk.Image,
    roi: sitk.Image,
    session: MultiresLandmarkSession,
    fit_mode: LandmarkFitMode,
    min_pairs: int,
) -> LandmarkFitResult:
    """Fit ROI->overview transform from paired landmarks stored as ZYX indices."""
    overview_pts_zyx, roi_pts_zyx = session.paired_points_zyx()
    n_pairs = overview_pts_zyx.shape[0]
    if n_pairs < min_pairs:
        msg = f"Need at least {min_pairs} landmark pairs, got {n_pairs}"
        raise ValueError(msg)

    overview_pts_um = zyx_points_to_physical(overview, overview_pts_zyx)
    roi_pts_um = zyx_points_to_physical(roi, roi_pts_zyx)
    matrix, _mse = fit_roi_to_overview_transform(
        roi_pts_um,
        overview_pts_um,
        fit_mode=fit_mode,
    )
    _fit_mse, errors = affine_point_errors(roi_pts_um, overview_pts_um, matrix)
    stats = summarize_point_errors(errors)
    return LandmarkFitResult(
        roi_to_overview_tform=matrix,
        rms_error_um=float(np.sqrt(np.mean(errors**2))),
        fit_point_errors_um=errors,
        fit_stats=stats,
    )


def overlap_box_from_landmark_transform(
    overview: sitk.Image,
    roi: sitk.Image,
    roi_to_overview: np.ndarray,
    *,
    margin_um: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Predict overlap crop box from landmark transform and volume bounds."""
    roi_min, roi_max = transformed_bounds_in_target_space(roi, roi_to_overview)
    overview_min, overview_max = physical_bounds(overview)
    overlap_min = np.maximum(roi_min, overview_min) - margin_um
    overlap_max = np.minimum(roi_max, overview_max) + margin_um
    if np.any(overlap_min >= overlap_max):
        msg = (
            "Landmark transform places ROI outside the overview bounds.\n"
            f"  ROI in overview space: {roi_min} .. {roi_max}\n"
            f"  Overview bounds: {overview_min} .. {overview_max}"
        )
        raise ValueError(msg)
    return overlap_min, overlap_max


def sitk_affine_from_matrix(reference_to_moving: np.ndarray) -> sitk.AffineTransform:
    """Build a SimpleITK affine for Resample (reference physical -> moving physical)."""
    transform = sitk.AffineTransform(3)
    transform.SetMatrix(reference_to_moving[:3, :3].reshape(-1).tolist())
    transform.SetTranslation(reference_to_moving[:3, 3].tolist())
    return transform


def resample_roi_onto_reference(
    roi: sitk.Image,
    reference: sitk.Image,
    roi_to_overview: np.ndarray,
) -> sitk.Image:
    """Resample ROI onto a reference overview crop using a ROI->overview landmark transform."""
    overview_to_roi = np.linalg.inv(roi_to_overview)
    transform = sitk_affine_from_matrix(overview_to_roi)
    return sitk.Resample(
        roi,
        reference,
        transform,
        sitk.sitkLinear,
        0.0,
        roi.GetPixelID(),
    )


def prepare_registration_pair_from_landmarks(
    overview: sitk.Image,
    roi: sitk.Image,
    *,
    session: MultiresLandmarkSession,
    fit_mode: LandmarkFitMode,
    min_pairs: int,
    margin_um: float,
) -> tuple[sitk.Image, sitk.Image, tuple[np.ndarray, np.ndarray], list[int], LandmarkFitResult]:
    """Crop overview and resample ROI using landmark-derived placement."""
    fit = fit_landmark_transform(
        overview=overview,
        roi=roi,
        session=session,
        fit_mode=fit_mode,
        min_pairs=min_pairs,
    )
    overlap_min, overlap_max = overlap_box_from_landmark_transform(
        overview,
        roi,
        fit.roi_to_overview_tform,
        margin_um=margin_um,
    )
    fixed_cropped, crop_start_index = crop_to_physical_box(overview, overlap_min, overlap_max)
    moving = resample_roi_onto_reference(roi, fixed_cropped, fit.roi_to_overview_tform)
    overlap_box = (overlap_min, overlap_max)
    return fixed_cropped, moving, overlap_box, crop_start_index, fit


def update_landmark_session_fit(
    session: MultiresLandmarkSession,
    fit: LandmarkFitResult,
) -> MultiresLandmarkSession:
    session.roi_to_overview_tform = fit.roi_to_overview_tform.tolist()
    session.rms_error_um = fit.rms_error_um
    session.fit_point_errors_um = fit.fit_point_errors_um.tolist()
    return session
