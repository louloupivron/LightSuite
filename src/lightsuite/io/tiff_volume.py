"""Read Z planes from single-file volumetric TIFF stacks (ImageJ / OME)."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import tifffile

_MEMMAP_CACHE: dict[str, np.memmap] = {}


def read_volumetric_slice(path: str | Path, z_index: int, stack_read_mode: str) -> np.ndarray:
    """Return one YX plane from a single-file Z stack (0-based z_index)."""
    key = str(Path(path).resolve())
    if stack_read_mode == "memmap_zyx":
        vol = _MEMMAP_CACHE.get(key)
        if vol is None:
            vol = tifffile.memmap(key)
            _MEMMAP_CACHE[key] = vol
        return np.asarray(vol[z_index], dtype=np.uint16)
    if stack_read_mode == "memmap_yxz":
        vol = _MEMMAP_CACHE.get(key)
        if vol is None:
            vol = tifffile.memmap(key)
            _MEMMAP_CACHE[key] = vol
        return np.asarray(vol[:, :, z_index], dtype=np.uint16)
    msg = f"Unsupported stack_read_mode: {stack_read_mode}"
    raise ValueError(msg)


def clear_memmap_cache() -> None:
    _MEMMAP_CACHE.clear()
