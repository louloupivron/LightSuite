"""Volume utilities for registration."""

from __future__ import annotations

import numpy as np
import tifffile
from pathlib import Path
from skimage.transform import resize


def load_registration_volume(path: Path) -> np.ndarray:
    """Load a multi-page registration TIFF as (Y, X, Z) float32."""
    path = path.expanduser()
    with tifffile.TiffFile(path) as tif:
        stack = tif.asarray()
    if stack.ndim == 2:
        return stack.astype(np.float32)
    if stack.ndim == 3:
        # Pages are Z; convert to Y, X, Z
        return np.moveaxis(stack.astype(np.float32), 0, -1)
    msg = f"Unexpected TIFF shape {stack.shape} in {path}"
    raise ValueError(msg)


def resize_atlas_volume(volume: np.ndarray, scale: float, *, nearest: bool = False) -> np.ndarray:
    """Resize 3D atlas volume by isotropic scale factor (atlasres/registres)."""
    if np.isclose(scale, 1.0):
        return volume
    h, w, z = volume.shape
    new_shape = (
        max(1, int(round(h * scale))),
        max(1, int(round(w * scale))),
        max(1, int(round(z * scale))),
    )
    order = 0 if nearest else 1
    out = resize(volume, new_shape, order=order, preserve_range=True, anti_aliasing=not nearest)
    return out.astype(volume.dtype, copy=False)


def permute_brain_volume(volume: np.ndarray, permvec: list[int]) -> np.ndarray:
    """Permute and flip volume axes (permuteBrainVolume.m)."""
    if len(permvec) != 3:
        msg = f"permvec must have length 3, got {permvec}"
        raise ValueError(msg)
    perm_order = [abs(v) - 1 for v in permvec]
    out = np.transpose(volume, perm_order)
    for dim, val in enumerate(permvec):
        if val < 0:
            out = np.flip(out, axis=dim)
    return np.ascontiguousarray(out)


def normalize_registration_volume(volume: np.ndarray) -> np.ndarray:
    """Scale sample volume to ~[0, 1] using central ROI (initializeRegistration.m)."""
    cent = np.array(volume.shape) // 2
    naround = max(1, int(round(min(cent) / 3)))
    sl = tuple(slice(max(0, c - naround), min(s, c + naround + 1)) for c, s in zip(cent, volume.shape, strict=True))
    center = volume[sl]
    top_val = float(np.quantile(center, 0.999)) * 2.0
    if top_val <= 0:
        top_val = float(volume.max()) or 1.0
    return ((volume.astype(np.float32) - 0.0) / top_val).astype(np.float32)
