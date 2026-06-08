"""AP slice-correspondence helpers for refining automatic control-point pairs."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import json
import numpy as np

from lightsuite.gui.slice_correspondence import SliceCorrespondence
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


def expected_atlas_plane(
    correspondence: SliceCorrespondence,
    sample_coord_1based: float,
    atlas_shape: tuple[int, int, int],
) -> int | None:
    """Interpolate atlas plane from sample position along the correspondence cut axis."""
    return correspondence.interpolate_atlas_plane(
        int(round(sample_coord_1based)),
        correspondence.cut_axis,
        atlas_axis_size(atlas_shape, correspondence.cut_axis),
    )


def ap_residuals_vox(
    cpsample: np.ndarray,
    cpatlas: np.ndarray,
    correspondence: SliceCorrespondence,
    original_trans: np.ndarray,
    atlas_shape: tuple[int, int, int],
) -> np.ndarray:
    """Per-pair |atlas_AP - expected_AP| in voxels (registration grid)."""
    if cpsample.shape[0] == 0:
        return np.zeros(0, dtype=float)
    sample_warped = transform_points_affinetform(
        np.asarray(cpsample, dtype=float),
        np.asarray(original_trans, dtype=float),
    )
    cut_axis = correspondence.cut_axis
    residuals = np.full(cpsample.shape[0], np.nan, dtype=float)
    for i, (s_pt, a_pt) in enumerate(zip(sample_warped, cpatlas, strict=True)):
        expected = expected_atlas_plane(correspondence, cloud_coord_1based(s_pt, cut_axis), atlas_shape)
        if expected is None:
            continue
        a_ap = cloud_coord_1based(np.asarray(a_pt, dtype=float), cut_axis)
        residuals[i] = abs(a_ap - expected)
    return residuals


@dataclass
class ApFilterStats:
    mode: str
    pairs_before: int
    pairs_after: int
    pairs_removed_ap: int
    tolerance_vox: float
    cut_axis: int
    n_anchors: int
    median_ap_residual_before_vox: float | None
    median_ap_residual_after_vox: float | None

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    def save(self, path: Path) -> None:
        path = path.expanduser()
        path.write_text(json.dumps(self.to_dict(), indent=2), encoding="utf-8")


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
    """Drop auto pairs whose atlas AP coordinate disagrees with slice correspondence."""
    sample = np.asarray(cpsample, dtype=float)
    atlas = np.asarray(cpatlas, dtype=float)
    n_before = int(sample.shape[0])
    if n_before == 0:
        stats = ApFilterStats(
            mode="ap_filter",
            pairs_before=0,
            pairs_after=0,
            pairs_removed_ap=0,
            tolerance_vox=tolerance_vox,
            cut_axis=correspondence.cut_axis,
            n_anchors=len(correspondence.confirmed_anchors()),
            median_ap_residual_before_vox=None,
            median_ap_residual_after_vox=None,
        )
        return sample, atlas, stats

    residuals = ap_residuals_vox(sample, atlas, correspondence, original_trans, atlas_shape)
    valid = ~np.isnan(residuals)
    before_median = float(np.median(residuals[valid])) if np.any(valid) else None

    def _apply_mask(mask: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        kept_residuals = residuals[mask & valid]
        after_median = float(np.median(kept_residuals)) if kept_residuals.size else None
        stats = ApFilterStats(
            mode="ap_filter",
            pairs_before=n_before,
            pairs_after=int(np.count_nonzero(mask)),
            pairs_removed_ap=int(n_before - np.count_nonzero(mask)),
            tolerance_vox=tolerance_vox,
            cut_axis=correspondence.cut_axis,
            n_anchors=len(correspondence.confirmed_anchors()),
            median_ap_residual_before_vox=before_median,
            median_ap_residual_after_vox=after_median,
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
