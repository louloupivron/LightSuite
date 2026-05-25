"""MetaImage (MHD/RAW) I/O for elastix."""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np


def write_mhd(volume: np.ndarray, base_path: Path, spacing_mm: float | list[float]) -> Path:
    """Write a 3D volume as MHD+RAW for elastix (volume indexed Y, X, Z)."""
    base_path = base_path.expanduser()
    base_path.parent.mkdir(parents=True, exist_ok=True)

    data = np.ascontiguousarray(volume.astype(np.float32, copy=False))
    if data.ndim != 3:
        msg = f"Expected 3D volume, got shape {data.shape}"
        raise ValueError(msg)

    sp = np.atleast_1d(np.asarray(spacing_mm, dtype=float))
    if sp.size == 1:
        sp = np.repeat(sp, 3)
    if sp.size != 3:
        msg = f"spacing_mm must be scalar or length-3, got {sp.size}"
        raise ValueError(msg)

    ny, nx, nz = data.shape
    raw_name = f"{base_path.name}.raw"
    raw_path = base_path.parent / raw_name
    flat = np.transpose(data, (1, 0, 2)).ravel(order="C")
    flat.tofile(raw_path)

    spacing_line = " ".join(f"{v:.12g}" for v in sp)
    header = (
        "ObjectType = Image\n"
        "NDims = 3\n"
        "BinaryData = True\n"
        "BinaryDataByteOrderMSB = False\n"
        "CompressedData = False\n"
        "TransformMatrix = 1 0 0 0 1 0 0 0 1\n"
        "Offset = 0 0 0\n"
        "CenterOfRotation = 0 0 0\n"
        f"ElementSpacing = {spacing_line}\n"
        f"DimSize = {nx} {ny} {nz}\n"
        "ElementType = MET_FLOAT\n"
        f"ElementDataFile = {raw_name}\n"
    )
    mhd_path = base_path.with_suffix(".mhd")
    mhd_path.write_text(header, encoding="utf-8")
    return mhd_path


def read_mhd_spacing(mhd_path: Path) -> np.ndarray:
    """Read ElementSpacing from an MHD header (mm)."""
    text = mhd_path.expanduser().read_text(encoding="utf-8")
    match = re.search(r"(?i)ElementSpacing\s*=\s*([^\r\n#]+)", text)
    if not match:
        msg = f"ElementSpacing not found in {mhd_path}"
        raise ValueError(msg)
    values = [float(v) for v in match.group(1).split()]
    if len(values) != 3:
        msg = f"Expected 3 spacing values in {mhd_path}, got {len(values)}"
        raise ValueError(msg)
    return np.asarray(values, dtype=float)
