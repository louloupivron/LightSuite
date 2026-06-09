"""Correspondence-informed affine correction before B-spline registration."""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import numpy as np

from scipy.spatial.distance import cdist

from lightsuite.gui.affine import fit_affine_transform, transform_points, transform_points_inverse
from lightsuite.gui.slice_correspondence import VOLUME_AXES, SliceAnchor, SliceCorrespondence
from lightsuite.gui.slices import volume_index_to_image
from lightsuite.registration.points_utils import subsample_point_pairs, thin_point_list
from lightsuite.registration.warp import swap_xy_transform


@dataclass
class CorrespondenceAffineStats:
    """Diagnostics for slice-correspondence affine pre-warp in register."""

    enabled: bool
    n_anchor_pairs: int
    n_axes: int
    median_residual_before_vox: float | None
    median_residual_after_vox: float | None
    p95_residual_after_vox: float | None
    applied: bool
    skip_reason: str | None = None

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    def save(self, path: Path) -> None:
        path = path.expanduser()
        path.write_text(json.dumps(self.to_dict(), indent=2), encoding="utf-8")


def _slice_tissue_centroid_yxz_1based(
    volume: np.ndarray,
    chooserow: np.ndarray,
) -> np.ndarray:
    """Return tissue centroid as 1-based (Y, X, Z) volume indices for one chooselist row."""
    chooserow = np.asarray(chooserow, dtype=int)
    sample_slice = volume_index_to_image(volume, chooserow)
    cut_axis = int(chooserow[1]) - 1
    plot_axes = [dim for dim in range(3) if dim != cut_axis]
    tissue = (
        sample_slice > float(np.quantile(sample_slice[sample_slice > 0], 0.05))
        if np.any(sample_slice > 0)
        else np.ones_like(sample_slice, dtype=bool)
    )
    yy, xx = np.where(tissue)
    if yy.size == 0:
        h, w = sample_slice.shape
        yy, xx = np.mgrid[0:h, 0:w]
        yy, xx = yy.ravel(), xx.ravel()
    pts = np.zeros((yy.size, 3), dtype=float)
    pts[:, plot_axes[1]] = yy.astype(float) + 1.0
    pts[:, plot_axes[0]] = xx.astype(float) + 1.0
    pts[:, cut_axis] = float(chooserow[0])
    return np.median(pts, axis=0)


def _atlas_target_from_anchor(
    sample_volume: np.ndarray,
    chooserow: np.ndarray,
    anchor: SliceAnchor,
) -> np.ndarray:
    """Atlas-side target point on the registration atlas grid (1-based Y, X, Z)."""
    centroid = _slice_tissue_centroid_yxz_1based(sample_volume, chooserow)
    cut_axis = int(chooserow[1]) - 1
    centroid[cut_axis] = float(anchor.atlas_plane)
    return centroid


def build_correspondence_anchor_pairs(
    correspondence: SliceCorrespondence,
    *,
    sample_warped: np.ndarray,
    original_trans_vol: np.ndarray,
    downfac: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Build atlas/sample landmark pairs from confirmed slice-correspondence anchors.

    Points are in the same spaces used by ``brain_register._prepare_control_points``:
    atlas coordinates are divided by ``downfac``; sample coordinates are inverse-mapped
    from the warped sample grid into permuted registration-volume space.
    """
    atlas_pts: list[np.ndarray] = []
    sample_pts: list[np.ndarray] = []
    original_trans_vol = np.asarray(original_trans_vol, dtype=float)

    for axis in VOLUME_AXES:
        for anchor in correspondence.confirmed_anchors(axis):
            chooserow = np.array([anchor.sample_index, axis, 1, 1], dtype=int)
            atlas_yxz = _atlas_target_from_anchor(sample_warped, chooserow, anchor)
            sample_centroid = _slice_tissue_centroid_yxz_1based(sample_warped, chooserow)
            sample_yxz = transform_points_inverse(
                sample_centroid.reshape(1, 3),
                original_trans_vol,
            )[0]
            atlas_pts.append(atlas_yxz / downfac)
            sample_pts.append(sample_yxz)

    if not atlas_pts:
        return np.zeros((0, 3)), np.zeros((0, 3))
    return np.vstack(atlas_pts), np.vstack(sample_pts)


def _residuals_vox(
    atlas_pts: np.ndarray,
    sample_pts: np.ndarray,
    tform: np.ndarray,
) -> np.ndarray:
    if atlas_pts.shape[0] == 0:
        return np.zeros(0, dtype=float)
    predicted = transform_points(atlas_pts, tform)
    return np.linalg.norm(predicted - sample_pts, axis=1)


def compose_correspondence_affine_correction(
    tform_aff: np.ndarray,
    corr_atlas: np.ndarray,
    corr_sample: np.ndarray,
) -> tuple[np.ndarray, CorrespondenceAffineStats]:
    """Compose a correspondence delta with the control-point affine fit.

    Fits a small atlas-space correction that maps the current affine prediction at each
    correspondence anchor onto the target atlas coordinate, then returns
    ``tform_aff @ inv(T_delta)`` so anchors align before the B-spline stage.
    """
    tform_aff = np.asarray(tform_aff, dtype=float)
    corr_atlas = np.asarray(corr_atlas, dtype=float)
    corr_sample = np.asarray(corr_sample, dtype=float)
    n_pairs = int(corr_atlas.shape[0])

    if n_pairs < 4:
        return tform_aff, CorrespondenceAffineStats(
            enabled=True,
            n_anchor_pairs=n_pairs,
            n_axes=0,
            median_residual_before_vox=None,
            median_residual_after_vox=None,
            p95_residual_after_vox=None,
            applied=False,
            skip_reason="fewer_than_4_anchor_pairs",
        )

    before = _residuals_vox(corr_atlas, corr_sample, tform_aff)
    before_median = float(np.median(before))
    if before_median < 1e-3:
        return tform_aff, CorrespondenceAffineStats(
            enabled=True,
            n_anchor_pairs=n_pairs,
            n_axes=0,
            median_residual_before_vox=before_median,
            median_residual_after_vox=before_median,
            p95_residual_after_vox=float(np.percentile(before, 95)),
            applied=False,
            skip_reason="already_aligned",
        )

    atlas_pred = transform_points_inverse(corr_sample, tform_aff)
    t_delta, _ = fit_affine_transform(atlas_pred, corr_atlas)
    tform_corrected = tform_aff @ np.linalg.inv(t_delta)
    after = _residuals_vox(corr_atlas, corr_sample, tform_corrected)
    after_median = float(np.median(after))
    if after_median > before_median:
        return tform_aff, CorrespondenceAffineStats(
            enabled=True,
            n_anchor_pairs=n_pairs,
            n_axes=0,
            median_residual_before_vox=before_median,
            median_residual_after_vox=before_median,
            p95_residual_after_vox=float(np.percentile(before, 95)),
            applied=False,
            skip_reason="correction_did_not_improve_fit",
        )

    return tform_corrected, CorrespondenceAffineStats(
        enabled=True,
        n_anchor_pairs=n_pairs,
        n_axes=0,
        median_residual_before_vox=before_median,
        median_residual_after_vox=after_median,
        p95_residual_after_vox=float(np.percentile(after, 95)),
        applied=True,
    )


def apply_slice_correspondence_affine(
    tform_aff: np.ndarray,
    correspondence: SliceCorrespondence | None,
    *,
    sample_warped: np.ndarray,
    original_trans: np.ndarray,
    downfac: float,
    enabled: bool = True,
) -> tuple[np.ndarray, CorrespondenceAffineStats]:
    """Apply correspondence correction when anchors are available."""
    if not enabled:
        return tform_aff, CorrespondenceAffineStats(
            enabled=False,
            n_anchor_pairs=0,
            n_axes=0,
            median_residual_before_vox=None,
            median_residual_after_vox=None,
            p95_residual_after_vox=None,
            applied=False,
            skip_reason="disabled_in_config",
        )
    if correspondence is None or not correspondence.has_confirmed_anchors():
        return tform_aff, CorrespondenceAffineStats(
            enabled=True,
            n_anchor_pairs=0,
            n_axes=0,
            median_residual_before_vox=None,
            median_residual_after_vox=None,
            p95_residual_after_vox=None,
            applied=False,
            skip_reason="no_confirmed_anchors",
        )

    original_trans_vol = swap_xy_transform(np.asarray(original_trans, dtype=float))
    corr_atlas, corr_sample = build_correspondence_anchor_pairs(
        correspondence,
        sample_warped=sample_warped,
        original_trans_vol=original_trans_vol,
        downfac=downfac,
    )
    tform_corrected, stats = compose_correspondence_affine_correction(
        tform_aff,
        corr_atlas,
        corr_sample,
    )
    stats.n_axes = correspondence.confirmed_axis_count()
    return tform_corrected, stats


@dataclass
class CorrespondenceLandmarkStats:
    """Diagnostics for slice-correspondence B-spline landmarks in register."""

    enabled: bool
    n_correspondence_pairs: int
    n_existing_pairs: int
    n_added_pairs: int
    n_merged_pairs: int
    n_skipped_near_existing: int
    applied: bool
    skip_reason: str | None = None

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    def save(self, path: Path) -> None:
        path = path.expanduser()
        path.write_text(json.dumps(self.to_dict(), indent=2), encoding="utf-8")


def _filter_pairs_far_from_existing(
    moving_pts: np.ndarray,
    fixed_pts: np.ndarray,
    existing_fixed: np.ndarray,
    min_distance_vox: float,
) -> tuple[np.ndarray, np.ndarray, int]:
    """Drop correspondence pairs whose fixed (sample) point is near an existing landmark."""
    if moving_pts.shape[0] == 0:
        return moving_pts, fixed_pts, 0
    if existing_fixed.shape[0] == 0:
        return moving_pts, fixed_pts, 0
    distances = cdist(fixed_pts, existing_fixed)
    keep = np.all(distances >= float(min_distance_vox), axis=1)
    skipped = int(np.count_nonzero(~keep))
    return moving_pts[keep], fixed_pts[keep], skipped


def append_correspondence_bspline_landmarks(
    cpaffine: np.ndarray,
    cptshistology: np.ndarray,
    correspondence: SliceCorrespondence | None,
    *,
    sample_warped: np.ndarray,
    original_trans: np.ndarray,
    downfac: float,
    tform_aff: np.ndarray,
    min_distance_vox: float,
    max_landmarks: int = 96,
    enabled: bool = True,
) -> tuple[np.ndarray, np.ndarray, CorrespondenceLandmarkStats]:
    """Add align-slices anchor pairs as extra Elastix B-spline landmarks.

    Moving landmarks are atlas points warped into sample space with ``tform_aff``;
    fixed landmarks are sample-space tissue centroids at each confirmed anchor.
    """
    cpaffine = np.asarray(cpaffine, dtype=float)
    cptshistology = np.asarray(cptshistology, dtype=float)
    n_existing = int(cptshistology.shape[0])

    if not enabled:
        return cpaffine, cptshistology, CorrespondenceLandmarkStats(
            enabled=False,
            n_correspondence_pairs=0,
            n_existing_pairs=n_existing,
            n_added_pairs=0,
            n_merged_pairs=n_existing,
            n_skipped_near_existing=0,
            applied=False,
            skip_reason="disabled_in_config",
        )
    if correspondence is None or not correspondence.has_confirmed_anchors():
        return cpaffine, cptshistology, CorrespondenceLandmarkStats(
            enabled=True,
            n_correspondence_pairs=0,
            n_existing_pairs=n_existing,
            n_added_pairs=0,
            n_merged_pairs=n_existing,
            n_skipped_near_existing=0,
            applied=False,
            skip_reason="no_confirmed_anchors",
        )

    original_trans_vol = swap_xy_transform(np.asarray(original_trans, dtype=float))
    corr_atlas, corr_sample = build_correspondence_anchor_pairs(
        correspondence,
        sample_warped=sample_warped,
        original_trans_vol=original_trans_vol,
        downfac=downfac,
    )
    n_corr = int(corr_atlas.shape[0])
    if n_corr == 0:
        return cpaffine, cptshistology, CorrespondenceLandmarkStats(
            enabled=True,
            n_correspondence_pairs=0,
            n_existing_pairs=n_existing,
            n_added_pairs=0,
            n_merged_pairs=n_existing,
            n_skipped_near_existing=0,
            applied=False,
            skip_reason="no_confirmed_anchors",
        )

    corr_moving = transform_points(corr_atlas, tform_aff)
    corr_fixed = corr_sample
    corr_moving, corr_fixed, skipped = _filter_pairs_far_from_existing(
        corr_moving,
        corr_fixed,
        cptshistology,
        min_distance_vox,
    )
    if corr_moving.shape[0]:
        keep = thin_point_list(corr_fixed, float(min_distance_vox))
        corr_moving = corr_moving[keep]
        corr_fixed = corr_fixed[keep]

    n_added = int(corr_moving.shape[0])
    if n_added == 0:
        return cpaffine, cptshistology, CorrespondenceLandmarkStats(
            enabled=True,
            n_correspondence_pairs=n_corr,
            n_existing_pairs=n_existing,
            n_added_pairs=0,
            n_merged_pairs=n_existing,
            n_skipped_near_existing=skipped,
            applied=False,
            skip_reason="all_pairs_near_existing_landmarks",
        )

    merged_moving = np.vstack([cpaffine, corr_moving]) if cpaffine.size else corr_moving
    merged_fixed = np.vstack([cptshistology, corr_fixed]) if cptshistology.size else corr_fixed
    if merged_moving.shape[0] > max_landmarks:
        # Keep all existing landmarks; subsample only the correspondence additions.
        if n_existing >= max_landmarks:
            merged_moving = cpaffine
            merged_fixed = cptshistology
            n_added = 0
        else:
            room = max_landmarks - n_existing
            corr_moving, corr_fixed = subsample_point_pairs(
                corr_moving,
                corr_fixed,
                max_points=room,
            )
            merged_moving = (
                np.vstack([cpaffine, corr_moving]) if cpaffine.size else corr_moving
            )
            merged_fixed = (
                np.vstack([cptshistology, corr_fixed]) if cptshistology.size else corr_fixed
            )
            n_added = int(corr_moving.shape[0])

    stats = CorrespondenceLandmarkStats(
        enabled=True,
        n_correspondence_pairs=n_corr,
        n_existing_pairs=n_existing,
        n_added_pairs=n_added,
        n_merged_pairs=int(merged_fixed.shape[0]),
        n_skipped_near_existing=skipped,
        applied=n_added > 0,
    )
    if not stats.applied:
        stats.skip_reason = "subsample_left_no_correspondence_pairs"
    return merged_moving, merged_fixed, stats
