"""Save microscopy volumes as multi-page TIFF (saveastiff.m equivalent)."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import tifffile


def save_registration_volume(volume: np.ndarray, path: Path, *, compression: str = "lzw") -> None:
    """Write a YXxZ uint volume as a compressed multi-page TIFF."""
    path = path.expanduser()
    path.parent.mkdir(parents=True, exist_ok=True)

    data = np.asarray(volume)
    if data.ndim != 3:
        msg = f"Expected 3D volume (Y, X, Z), got shape {data.shape}"
        raise ValueError(msg)

    if data.dtype != np.uint16:
        data = np.clip(data, 0, np.iinfo(np.uint16).max).astype(np.uint16)

    # tifffile expects (pages, Y, X) for imagej-style stacks
    stack = np.moveaxis(data, -1, 0)
    tifffile.imwrite(
        path,
        stack,
        photometric="minisblack",
        compression=compression,
        metadata={"axes": "ZYX"},
    )
