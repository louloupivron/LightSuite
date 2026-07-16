"""Optional warp-canvas helpers (not used by MATLAB-parity brain registration).

Brain registration follows ``multiobjRegistration.m``: atlas volumes are warped onto
the sample grid via ``imwarp`` / :func:`~lightsuite.registration.warp.imwarp_volume`
with ``OutputView = size(volume)``. These utilities remain for tests and a future
explicit crop/pad stage.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class WarpCanvasPadding:
    """Symmetric padding applied around the sample grid for affine/B-spline warps."""

    pad_before: tuple[int, int, int]
    pad_after: tuple[int, int, int]

    @property
    def is_zero(self) -> bool:
        return self.pad_before == (0, 0, 0) and self.pad_after == (0, 0, 0)

    def padded_shape(self, original_shape: tuple[int, int, int]) -> tuple[int, int, int]:
        return tuple(
            original + before + after
            for original, before, after in zip(
                original_shape,
                self.pad_before,
                self.pad_after,
                strict=True,
            )
        )


def compute_vd_warp_canvas_padding(
    sample_shape: tuple[int, int, int],
    atlas_shape: tuple[int, int, int],
    *,
    vd_axis: int = 2,
    ap_axis: int = 1,
) -> WarpCanvasPadding:
    """Return symmetric VD padding so warped atlas labels are not clipped at brain poles.

    Padding is chosen to match the template AP:VD aspect ratio (brain grid convention
    after ``permute_brain_volume``), mirroring the relative VD headroom in the
    Perens/Gubra template volume.
    """
    if vd_axis not in {0, 1, 2} or ap_axis not in {0, 1, 2}:
        msg = f"vd_axis and ap_axis must be 0, 1, or 2, got {vd_axis} and {ap_axis}"
        raise ValueError(msg)

    sample_vd = sample_shape[vd_axis]
    sample_ap = sample_shape[ap_axis]
    atlas_vd = atlas_shape[vd_axis]
    atlas_ap = atlas_shape[ap_axis]
    if sample_vd <= 0 or sample_ap <= 0 or atlas_vd <= 0 or atlas_ap <= 0:
        return WarpCanvasPadding((0, 0, 0), (0, 0, 0))

    target_vd = int(round(sample_ap * atlas_vd / atlas_ap))
    total_pad = max(0, target_vd - sample_vd)
    if total_pad <= 0:
        return WarpCanvasPadding((0, 0, 0), (0, 0, 0))

    pad_before = [0, 0, 0]
    pad_after = [0, 0, 0]
    pad_before[vd_axis] = total_pad // 2
    pad_after[vd_axis] = total_pad - pad_before[vd_axis]
    return WarpCanvasPadding(tuple(pad_before), tuple(pad_after))


def pad_volume_for_warp_canvas(
    volume: np.ndarray,
    padding: WarpCanvasPadding,
    *,
    constant: float = 0.0,
) -> np.ndarray:
    """Embed ``volume`` in a larger array with background on each padded side."""
    if padding.is_zero:
        return volume
    pad_width = tuple(zip(padding.pad_before, padding.pad_after, strict=True))
    return np.pad(volume, pad_width, mode="constant", constant_values=constant)


def crop_from_warp_canvas(
    volume: np.ndarray,
    padding: WarpCanvasPadding,
    original_shape: tuple[int, int, int],
) -> np.ndarray:
    """Extract the original sample region from a padded warp canvas."""
    if padding.is_zero:
        return volume
    slices = tuple(
        slice(before, before + original)
        for before, original in zip(padding.pad_before, original_shape, strict=True)
    )
    return volume[slices]


def offset_volume_indices(
    points_yxz: np.ndarray,
    pad_before: tuple[int, int, int],
) -> np.ndarray:
    """Shift sample-space (Y, X, Z) indices after low-side canvas padding."""
    if points_yxz.size == 0 or pad_before == (0, 0, 0):
        return points_yxz
    out = np.asarray(points_yxz, dtype=float).copy()
    out[:, 0] += pad_before[0]
    out[:, 1] += pad_before[1]
    out[:, 2] += pad_before[2]
    return out
