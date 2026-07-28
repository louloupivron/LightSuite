"""Registration working-canvas helpers (Option A — working grid vs canonical atlas)."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from lightsuite.config.models import RegistrationCanvasMode
from lightsuite.registration.content_bbox import ContentBox, crop_volume_yxz


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


@dataclass(frozen=True)
class RegistrationCanvas:
    """Working grid for elastix relative to the (possibly cropped) sample volume."""

    mode: str
    pad_before: tuple[int, int, int]
    pad_after: tuple[int, int, int]
    sample_crop_start: tuple[int, int, int]
    working_shape: tuple[int, int, int]

    @property
    def is_identity(self) -> bool:
        return (
            self.pad_before == (0, 0, 0)
            and self.pad_after == (0, 0, 0)
            and self.sample_crop_start == (0, 0, 0)
        )

    def to_checkpoint_dict(self) -> dict:
        return {
            "mode": self.mode,
            "pad_before": list(self.pad_before),
            "pad_after": list(self.pad_after),
            "sample_crop_start": list(self.sample_crop_start),
            "working_shape": list(self.working_shape),
        }

    @classmethod
    def from_checkpoint_dict(cls, raw: dict | None) -> RegistrationCanvas | None:
        if raw is None:
            return None
        return cls(
            mode=str(raw.get("mode", "off")),
            pad_before=tuple(int(v) for v in raw["pad_before"]),
            pad_after=tuple(int(v) for v in raw["pad_after"]),
            sample_crop_start=tuple(int(v) for v in raw.get("sample_crop_start", [0, 0, 0])),
            working_shape=tuple(int(v) for v in raw["working_shape"]),
        )


def _intersection_box(a_shape: tuple[int, int, int], b_shape: tuple[int, int, int]) -> ContentBox:
    sy = min(a_shape[0], b_shape[0])
    sx = min(a_shape[1], b_shape[1])
    sz = min(a_shape[2], b_shape[2])
    return ContentBox(0, sy - 1, 0, sx - 1, 0, sz - 1)


def compute_registration_canvas(
    sample_shape: tuple[int, int, int],
    atlas_shape: tuple[int, int, int],
    mode: RegistrationCanvasMode,
) -> RegistrationCanvas:
    """Choose elastix working grid given cropped sample and trimmed atlas shapes."""
    if mode == RegistrationCanvasMode.OFF:
        return RegistrationCanvas(
            mode=mode.value,
            pad_before=(0, 0, 0),
            pad_after=(0, 0, 0),
            sample_crop_start=(0, 0, 0),
            working_shape=sample_shape,
        )

    if mode == RegistrationCanvasMode.PAD:
        pad_before = [0, 0, 0]
        pad_after = [0, 0, 0]
        for axis, (s, a) in enumerate(zip(sample_shape, atlas_shape, strict=True)):
            if a > s:
                total = a - s
                pad_before[axis] = total // 2
                pad_after[axis] = total - pad_before[axis]
        working = tuple(
            s + b + a for s, b, a in zip(sample_shape, pad_before, pad_after, strict=True)
        )
        return RegistrationCanvas(
            mode=mode.value,
            pad_before=tuple(pad_before),
            pad_after=tuple(pad_after),
            sample_crop_start=(0, 0, 0),
            working_shape=working,
        )

    if mode == RegistrationCanvasMode.CROP:
        box = _intersection_box(sample_shape, atlas_shape)
        return RegistrationCanvas(
            mode=mode.value,
            pad_before=(0, 0, 0),
            pad_after=(0, 0, 0),
            sample_crop_start=box.start_yxz,
            working_shape=box.size_yxz,
        )

    if mode == RegistrationCanvasMode.UNION:
        working = tuple(max(s, a) for s, a in zip(sample_shape, atlas_shape, strict=True))
        pad_before = [0, 0, 0]
        pad_after = [0, 0, 0]
        for axis, (s, w) in enumerate(zip(sample_shape, working, strict=True)):
            if w > s:
                total = w - s
                pad_before[axis] = total // 2
                pad_after[axis] = total - pad_before[axis]
        return RegistrationCanvas(
            mode=mode.value,
            pad_before=tuple(pad_before),
            pad_after=tuple(pad_after),
            sample_crop_start=(0, 0, 0),
            working_shape=working,
        )

    msg = f"Unknown canvas mode {mode!r}"
    raise ValueError(msg)


def apply_canvas_sample_crop(volume: np.ndarray, canvas: RegistrationCanvas) -> np.ndarray:
    if canvas.sample_crop_start == (0, 0, 0) and volume.shape == canvas.working_shape:
        return volume
    if canvas.sample_crop_start != (0, 0, 0):
        box = ContentBox(
            canvas.sample_crop_start[0],
            canvas.sample_crop_start[0] + canvas.working_shape[0] - 1,
            canvas.sample_crop_start[1],
            canvas.sample_crop_start[1] + canvas.working_shape[1] - 1,
            canvas.sample_crop_start[2],
            canvas.sample_crop_start[2] + canvas.working_shape[2] - 1,
        )
        volume = crop_volume_yxz(volume, box)
    padding = WarpCanvasPadding(canvas.pad_before, canvas.pad_after)
    return pad_volume_for_warp_canvas(volume, padding)


def undo_canvas_sample_crop(
    volume: np.ndarray,
    canvas: RegistrationCanvas,
    full_shape: tuple[int, int, int],
) -> np.ndarray:
    """Inverse of :func:`apply_canvas_sample_crop` for the unpadded working region."""
    if canvas.is_identity and volume.shape == full_shape:
        return volume

    inner_shape = tuple(
        canvas.working_shape[i] - canvas.pad_before[i] - canvas.pad_after[i]
        for i in range(3)
    )
    unpadded = volume
    padding = WarpCanvasPadding(canvas.pad_before, canvas.pad_after)
    if not padding.is_zero:
        unpadded = crop_from_warp_canvas(volume, padding, inner_shape)

    if canvas.sample_crop_start == (0, 0, 0):
        if unpadded.shape == full_shape:
            return unpadded
        msg = (
            f"Cannot embed unpadded volume shape {unpadded.shape} "
            f"into full_shape {full_shape} without crop offset"
        )
        raise ValueError(msg)

    out = np.zeros(full_shape, dtype=volume.dtype)
    sy, sx, sz = canvas.sample_crop_start
    ey = sy + unpadded.shape[0]
    ex = sx + unpadded.shape[1]
    ez = sz + unpadded.shape[2]
    out[sy:ey, sx:ex, sz:ez] = unpadded
    return out
