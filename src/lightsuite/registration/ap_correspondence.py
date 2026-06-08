"""Slice-correspondence helpers for refining automatic control-point pairs."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import json
import numpy as np

from lightsuite.gui.slice_correspondence import VOLUME_AXES, SliceCorrespondence
from lightsuite.registration.warp import transform_points_affinetform

# cut_axis (1-based Y/X/Z) -> cloud column index [x, y, z]
_CLOUD_AXIS_FOR_CUT: dict[int, int] = {1: 1, 2: 0, 3: 2}


def cloud_coord_1based(point: np.ndarray, cut_axis: int) -> float:
    """1-based coordinate along ``cut_axis`` for a cloud point ``[x, y, z]`` (0-based)."""
    col = _CLOUD_AXIS_FOR_CUT[int(cut_axis)]
    return float(point[col] + 1.0)


def atlas_axis_size(atlas_shape: tuple[int, int, int], cut_axis: int) -> int:
    """Volume extent along ``cut_axis`` (1-based)."""
    return int(atlas_shape[cut_axis - 1])


def axis_residual_vox(
    sample_warped: np.ndarray,
    atlas_point: np.ndarray,
    correspondence: SliceCorrespondence,
    atlas_shape: tuple[int, int, int],
    cut_axis: int,
) -> float:
    """|atlas_coord - expected_coord| along one volume axis."""
    sample_coord = cloud_coord_1based(sample_warped, cut_axis)
    expected = correspondence.interpolate_atlas_plane(
        int(round(sample_coord)),
        cut_axis,
        atlas_axis_size(atlas_shape, cut_axis),
    )
    if expected is None:
        return float("nan")
    atlas_coord = cloud_coord_1based(np.asarray(atlas_point, dtype=float), cut_axis)
    return abs(atlas_coord - expected)


def multi_axis_residuals_vox(
    cpsample: np.ndarray,
    cpatlas: np.ndarray,
    correspondence: SliceCorrespondence,
    original_trans: np.ndarray,
    atlas_shape: tuple[int, int, int],
) -> tuple[np.ndarray, dict[int, np.ndarray]]:
    """Per-pair max residual across axes with confirmed correspondence curves."""
    if cpsample.shape[0] == 0:
        return np.zeros(0, dtype=float), {}
    sample_warped = transform_points_affinetform(
        np.asarray(cpsample, dtype=float),
        np.asarray(original_trans, dtype=float),
    )
    active_axes = [axis for axis in VOLUME_AXES if correspondence.has_confirmed_anchors(axis)]
    if not active_axes:
        return np.full(cpsample.shape[0], np.nan, dtype=float), {}

    per_axis: dict[int, np.ndarray] = {}
    max_residual = np.zeros(cpsample.shape[0], dtype=float)
    for axis in active_axes:
        axis_residuals = np.array(
            [
                axis_residual_vox(s_pt, a_pt, correspondence, atlas_shape, axis)
                for s_pt, a_pt in zip(sample_warped, cpatlas, strict=True)
            ],
            dtype=float,
        )
        per_axis[axis] = axis_residuals
        valid = ~np.isnan(axis_residuals)
        max_residual = np.where(
            valid,
            np.maximum(max_residual, axis_residuals),
            max_residual,
        )
    if not np.any(~np.isnan(max_residual)):
        max_residual[:] = np.nan
    return max_residual, per_axis


def ap_residuals_vox(
    cpsample: np.ndarray,
    cpatlas: np.ndarray,
    correspondence: SliceCorrespondence,
    original_trans: np.ndarray,
    atlas_shape: tuple[int, int, int],
) -> np.ndarray:
    """Per-pair max residual across confirmed axes (registration-grid voxels)."""
    residuals, _ = multi_axis_residuals_vox(
        cpsample,
        cpatlas,
        correspondence,
        original_trans,
        atlas_shape,
    )
    return residuals


@dataclass
class ApFilterStats:
    mode: str
    pairs_before: int
    pairs_after: int
    pairs_removed: int
    tolerance_vox: float
    active_axes: list[int]
    n_anchors_by_axis: dict[str, int]
    median_residual_before_vox: float | None
    median_residual_after_vox: float | None
    median_residual_by_axis_before_vox: dict[str, float | None] | None = None
    median_residual_by_axis_after_vox: dict[str, float | None] | None = None

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    def save(self, path: Path) -> None:
        path = path.expanduser()
        path.write_text(json.dumps(self.to_dict(), indent=2), encoding="utf-8")


def _median_or_none(values: np.ndarray) -> float | None:
    valid = values[~np.isnan(values)]
    if valid.size == 0:
        return None
    return float(np.median(valid))


def filter_pairs_by_ap_correspondence(
    cpsample: np.ndarray,
    cpatlas: np.ndarray,
    correspondence: SliceCorrespondence,
    original_trans: np.ndarray,
    atlas_shape: tuple[int, int, int],
    *,
    tolerance_vox: float = 12.0,
    min_pairs_kept: int = 24,
) -> tuple[np.ndarray, np.ndarray, ApFilterStats]:
    """Drop auto pairs that disagree with slice correspondence on any active axis."""
    sample = np.asarray(cpsample, dtype=float)
    atlas = np.asarray(cpatlas, dtype=float)
    n_before = int(sample.shape[0])
    active_axes = [axis for axis in VOLUME_AXES if correspondence.has_confirmed_anchors(axis)]
    n_anchors = {str(axis): len(correspondence.confirmed_anchors(axis)) for axis in active_axes}

    if n_before == 0:
        stats = ApFilterStats(
            mode="multi_axis_filter",
            pairs_before=0,
            pairs_after=0,
            pairs_removed=0,
            tolerance_vox=tolerance_vox,
            active_axes=active_axes,
            n_anchors_by_axis=n_anchors,
            median_residual_before_vox=None,
            median_residual_after_vox=None,
        )
        return sample, atlas, stats

    residuals, per_axis = multi_axis_residuals_vox(
        sample,
        atlas,
        correspondence,
        original_trans,
        atlas_shape,
    )
    valid = ~np.isnan(residuals)
    before_median = _median_or_none(residuals)
    before_by_axis = {
        str(axis): _median_or_none(axis_vals) for axis, axis_vals in per_axis.items()
    }

    def _apply_mask(mask: np.ndarray) -> tuple[np.ndarray, np.ndarray, ApFilterStats]:
        kept_residuals = residuals[mask & valid]
        after_median = _median_or_none(kept_residuals)
        after_by_axis: dict[str, float | None] = {}
        for axis, axis_vals in per_axis.items():
            kept_axis = axis_vals[mask & ~np.isnan(axis_vals)]
            after_by_axis[str(axis)] = _median_or_none(kept_axis)
        stats = ApFilterStats(
            mode="multi_axis_filter",
            pairs_before=n_before,
            pairs_after=int(np.count_nonzero(mask)),
            pairs_removed=int(n_before - np.count_nonzero(mask)),
            tolerance_vox=tolerance_vox,
            active_axes=active_axes,
            n_anchors_by_axis=n_anchors,
            median_residual_before_vox=before_median,
            median_residual_after_vox=after_median,
            median_residual_by_axis_before_vox=before_by_axis,
            median_residual_by_axis_after_vox=after_by_axis,
        )
        return sample[mask], atlas[mask], stats

    tol = float(tolerance_vox)
    mask = valid & (residuals <= tol)
    if int(np.count_nonzero(mask)) < min_pairs_kept:
        relaxed = valid & (residuals <= tol * 1.5)
        if int(np.count_nonzero(relaxed)) >= min_pairs_kept:
            mask = relaxed
        else:
            mask = valid if int(np.count_nonzero(valid)) >= min_pairs_kept else np.ones(n_before, dtype=bool)
    return _apply_mask(mask)
