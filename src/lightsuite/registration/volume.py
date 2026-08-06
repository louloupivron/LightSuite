"""Volume utilities for registration."""

from __future__ import annotations

import numpy as np
import tifffile
from pathlib import Path
from skimage.transform import resize


def _stack_tiff_pages(tif: tifffile.TiffFile, path: Path) -> np.ndarray:
    """Stack TIFF planes to (Z, Y, X), handling modern tifffile series layout."""
    if len(tif.pages) == 0:
        msg = f"No TIFF pages in {path}"
        raise ValueError(msg)

    if len(tif.pages) == 1:
        return np.asarray(tif.pages[0].asarray(), dtype=np.float32)

    # tifffile >= 2024 often exposes each IFD as its own 2D series; asarray() is then
    # only the first plane while Fiji/ImageJ still show the full stack.
    if len(tif.series) == 1 and tif.series[0].ndim == 3:
        data = np.asarray(tif.series[0].asarray(), dtype=np.float32)
        axes = tif.series[0].axes
        if axes in {"ZYX", "IYX"}:
            return data
        if axes == "XYZ":
            return np.moveaxis(data, -1, 0)
        msg = f"Unsupported TIFF series axes {axes!r} in {path}"
        raise ValueError(msg)

    planes = [np.asarray(page.asarray(), dtype=np.float32) for page in tif.pages]
    if any(p.ndim != 2 for p in planes):
        msg = f"Expected 2D TIFF pages, got shapes {[p.shape for p in planes[:3]]} in {path}"
        raise ValueError(msg)
    return np.stack(planes, axis=0)


def load_registration_volume(path: Path) -> np.ndarray:
    """Load a multi-page registration TIFF as (Y, X, Z) float32."""
    path = path.expanduser()
    with tifffile.TiffFile(path) as tif:
        stack = _stack_tiff_pages(tif, path)
    if stack.ndim == 2:
        return stack
    if stack.ndim == 3:
        return np.moveaxis(stack, 0, -1)
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


def unpermute_brain_volume(volume: np.ndarray, permvec: list[int]) -> np.ndarray:
    """Inverse of :func:`permute_brain_volume`."""
    if len(permvec) != 3:
        msg = f"permvec must have length 3, got {permvec}"
        raise ValueError(msg)
    out = np.asarray(volume)
    for dim, val in enumerate(permvec):
        if val < 0:
            out = np.flip(out, axis=dim)
    perm_order = [abs(v) - 1 for v in permvec]
    inverse_order = [0, 0, 0]
    for new_axis, old_axis in enumerate(perm_order):
        inverse_order[old_axis] = new_axis
    return np.ascontiguousarray(np.transpose(out, inverse_order))


def normalize_registration_volume(volume: np.ndarray) -> np.ndarray:
    """Scale sample volume to ~[0, 1] using central ROI (initializeRegistration.m).

    MATLAB indexes a diagonal through the center (``centind(:,1)``, ``centind(:,2)``,
    ``centind(:,3)``), not a full cubic ROI.
    """
    vol = volume.astype(np.float32, copy=False)
    cent_px = np.round(np.array(vol.shape, dtype=float) / 2.0).astype(int)
    naround = max(1, int(round(float(np.min(cent_px)) / 3.0)))
    offsets = np.arange(-naround, naround + 1, dtype=int)
    idx0 = cent_px[0] + offsets
    idx1 = cent_px[1] + offsets
    idx2 = cent_px[2] + offsets
    center_samples = vol[idx0, idx1, idx2]
    top_val = float(np.quantile(center_samples, 0.999)) * 2.0
    if top_val <= 0:
        top_val = float(vol.max()) or 1.0
    return (vol / top_val).astype(np.float32)
