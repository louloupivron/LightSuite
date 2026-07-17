"""Foreground bounding boxes for atlas trim and sample content crop."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.ndimage import binary_dilation, label


@dataclass(frozen=True)
class ContentBox:
    """Inclusive axis-aligned crop on a (Y, X, Z) volume."""

    y0: int
    y1: int
    x0: int
    x1: int
    z0: int
    z1: int

    @property
    def start_yxz(self) -> tuple[int, int, int]:
        return (self.y0, self.x0, self.z0)

    @property
    def size_yxz(self) -> tuple[int, int, int]:
        return (self.y1 - self.y0 + 1, self.x1 - self.x0 + 1, self.z1 - self.z0 + 1)

    def slices(self) -> tuple[slice, slice, slice]:
        return (slice(self.y0, self.y1 + 1), slice(self.x0, self.x1 + 1), slice(self.z0, self.z1 + 1))

    def to_dict(self) -> dict[str, list[int]]:
        return {
            "crop_start": list(self.start_yxz),
            "crop_size": list(self.size_yxz),
            "box_yxz": [self.y0, self.y1, self.x0, self.x1, self.z0, self.z1],
        }

    def to_manual_list(self) -> list[int]:
        """Return ``[y0, y1, x0, x1, z0, z1]`` for pipeline YAML."""
        return [self.y0, self.y1, self.x0, self.x1, self.z0, self.z1]

    @classmethod
    def full_volume(cls, shape: tuple[int, int, int]) -> ContentBox:
        sy, sx, sz = shape
        return cls(y0=0, y1=sy - 1, x0=0, x1=sx - 1, z0=0, z1=sz - 1)

    @classmethod
    def from_manual_box(cls, box: list[int]) -> ContentBox:
        """Parse ``[y0, y1, x0, x1, z0, z1]`` inclusive indices."""
        if len(box) != 6:
            msg = f"content_box must have 6 values [y0,y1,x0,x1,z0,z1], got {len(box)}"
            raise ValueError(msg)
        y0, y1, x0, x1, z0, z1 = (int(v) for v in box)
        if y1 < y0 or x1 < x0 or z1 < z0:
            msg = f"Invalid content_box ordering: {box}"
            raise ValueError(msg)
        return cls(y0=y0, y1=y1, x0=x0, x1=x1, z0=z0, z1=z1)

    @classmethod
    def from_dict(cls, raw: dict) -> ContentBox:
        if "box_yxz" in raw:
            return cls.from_manual_box(list(raw["box_yxz"]))
        start = tuple(int(v) for v in raw["crop_start"])
        size = tuple(int(v) for v in raw["crop_size"])
        return cls(
            y0=start[0],
            y1=start[0] + size[0] - 1,
            x0=start[1],
            x1=start[1] + size[1] - 1,
            z0=start[2],
            z1=start[2] + size[2] - 1,
        )


def _expand_box(
    y0: int,
    y1: int,
    x0: int,
    x1: int,
    z0: int,
    z1: int,
    shape: tuple[int, int, int],
    margin_vox: int,
) -> ContentBox:
    sy, sx, sz = shape
    return ContentBox(
        y0=max(0, y0 - margin_vox),
        y1=min(sy - 1, y1 + margin_vox),
        x0=max(0, x0 - margin_vox),
        x1=min(sx - 1, x1 + margin_vox),
        z0=max(0, z0 - margin_vox),
        z1=min(sz - 1, z1 + margin_vox),
    )


def bbox_from_mask(mask: np.ndarray, *, margin_vox: int = 0) -> ContentBox | None:
    """Tight box around ``True`` voxels in a (Y, X, Z) mask."""
    coords = np.argwhere(np.asarray(mask, dtype=bool))
    if coords.size == 0:
        return None
    y0, x0, z0 = coords.min(axis=0)
    y1, x1, z1 = coords.max(axis=0)
    return _expand_box(int(y0), int(y1), int(x0), int(x1), int(z0), int(z1), mask.shape, margin_vox)


def bbox_from_annotation(
    annotation: np.ndarray,
    *,
    margin_vox: int = 0,
) -> ContentBox | None:
    return bbox_from_mask(np.asarray(annotation) > 0, margin_vox=margin_vox)


def bbox_from_template_threshold(
    template: np.ndarray,
    *,
    margin_vox: int = 0,
    quantile: float = 0.05,
) -> ContentBox | None:
    vol = np.asarray(template, dtype=np.float32)
    flat = vol.ravel()
    if flat.size == 0:
        return None
    probe = flat if flat.size <= 50_000 else flat[np.random.default_rng(0).choice(flat.size, 50_000, replace=False)]
    threshold = float(np.quantile(probe, quantile)) * 2.0
    return bbox_from_mask(vol > max(threshold, 0.0), margin_vox=margin_vox)


def crop_volume_yxz(volume: np.ndarray, box: ContentBox) -> np.ndarray:
    return np.asanyarray(volume)[box.slices()]


def sample_foreground_bbox(
    volume: np.ndarray,
    *,
    margin_vox: int = 8,
    trim_z: bool = True,
) -> ContentBox | None:
    """Largest connected foreground component on a registration-resolution (Y, X, Z) volume."""
    vol = np.asarray(volume)
    if vol.ndim != 3:
        msg = f"Expected 3D volume, got shape {vol.shape}"
        raise ValueError(msg)

    flat = vol.ravel()
    probe_n = min(20_000, flat.size)
    if probe_n == 0:
        return None
    rng = np.random.default_rng(1)
    probe = flat[rng.choice(flat.size, probe_n, replace=False)] if flat.size > probe_n else flat
    positive = probe[probe > 0]
    if positive.size == 0:
        return None
    mask = vol > 0
    if not np.any(mask):
        backval = int(np.bincount(positive.astype(np.int64)).argmax())
        filled = vol.copy()
        filled[filled == 0] = backval
        mask = filled > backval
    if not np.any(mask):
        return None
    mask = binary_dilation(mask, structure=np.ones((3, 3, 3), dtype=bool))
    labeled, n = label(mask)
    if n == 0:
        return None
    counts = np.bincount(labeled.ravel())
    counts[0] = 0
    keep = int(np.argmax(counts))
    fg = labeled == keep

    box = bbox_from_mask(fg, margin_vox=0)
    if box is None:
        return None

    if trim_z:
        area = fg.sum(axis=(0, 1)).astype(float)
        med = float(np.median(area))
        rstd = float(np.std(area)) if area.size > 1 else 0.0
        high = area > (med + 3.0 * rstd)
        if np.any(high) and float(high.mean()) > 0.01:
            z_idx = np.where(high)[0]
            z0 = int(z_idx[0])
            z1 = int(z_idx[-1])
            box = ContentBox(box.y0, box.y1, box.x0, box.x1, z0, z1)

    return _expand_box(box.y0, box.y1, box.x0, box.x1, box.z0, box.z1, vol.shape, margin_vox)


def scale_box_to_lower_resolution(box: ContentBox, factor: float) -> ContentBox:
    """Scale a native-resolution box to a downsampled grid (e.g. atlas → tvreg)."""
    if factor <= 0:
        msg = f"factor must be positive, got {factor}"
        raise ValueError(msg)
    inv = factor
    return ContentBox(
        y0=int(np.floor(box.y0 / inv)),
        y1=int(np.ceil(box.y1 / inv)),
        x0=int(np.floor(box.x0 / inv)),
        x1=int(np.ceil(box.x1 / inv)),
        z0=int(np.floor(box.z0 / inv)),
        z1=int(np.ceil(box.z1 / inv)),
    )
