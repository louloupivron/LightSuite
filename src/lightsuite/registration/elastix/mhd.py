"""MetaImage (MHD/RAW) I/O for elastix."""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np


def scale_volume_for_elastix_mi(volume: np.ndarray, *, quantile: float = 0.999) -> np.ndarray:
    """Scale a volume to uint16 for elastix MI (matches MATLAB ``mhd_write`` ushort inputs).

    Fixed (sample) and moving (affine-warped atlas template) often occupy very different
    raw intensity ranges. Mattes MI uses a fixed bin count; when one image spans ~0–15k
    and the other ~0–500, the moving histogram collapses to a few bins at full resolution
    while pyramid downsampling artificially spreads it (R2 looks fine, R3 MI ≈ 0).
    Per-volume quantile scaling puts both images on a comparable 16-bit range.
    """
    vol = np.asarray(volume, dtype=np.float32)
    positive = vol[vol > 0]
    if positive.size == 0:
        return np.zeros(vol.shape, dtype=np.uint16)
    hi = float(np.quantile(positive, quantile))
    if hi <= 0:
        hi = float(vol.max()) or 1.0
    return np.clip(vol / hi * 65535.0, 0, 65535).astype(np.uint16)


def write_mhd(volume: np.ndarray, base_path: Path, spacing_mm: float | list[float]) -> Path:
    """Write a 3D volume as MHD+RAW for elastix (volume indexed Y, X, Z)."""
    base_path = base_path.expanduser()
    base_path.parent.mkdir(parents=True, exist_ok=True)

    if volume.dtype == np.uint16:
        data = np.ascontiguousarray(volume)
        element_type = "MET_USHORT"
    else:
        data = np.ascontiguousarray(volume.astype(np.float32, copy=False))
        element_type = "MET_FLOAT"
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
        f"ElementType = {element_type}\n"
        f"ElementDataFile = {raw_name}\n"
    )
    mhd_path = base_path.with_suffix(".mhd")
    mhd_path.write_text(header, encoding="utf-8")
    return mhd_path


_MHD_ELEMENT_DTYPES: dict[str, np.dtype] = {
    "MET_FLOAT": np.dtype(np.float32),
    "MET_DOUBLE": np.dtype(np.float64),
    "MET_SHORT": np.dtype(np.int16),
    "MET_USHORT": np.dtype(np.uint16),
    "MET_INT": np.dtype(np.int32),
    "MET_UINT": np.dtype(np.uint32),
    "MET_CHAR": np.dtype(np.int8),
    "MET_UCHAR": np.dtype(np.uint8),
}


def _mhd_header_field(text: str, key: str) -> str | None:
    match = re.search(rf"(?im)^{re.escape(key)}\s*=\s*([^\r\n#]+)", text)
    return match.group(1).strip() if match else None


def _mhd_element_dtype(element_type: str) -> np.dtype:
    key = element_type.strip().upper()
    if key not in _MHD_ELEMENT_DTYPES:
        msg = f"Unsupported MHD ElementType {element_type!r}"
        raise ValueError(msg)
    return _MHD_ELEMENT_DTYPES[key]


def read_mhd_volume(mhd_path: Path) -> np.ndarray:
    """Load a 3D volume from MHD+RAW (returns Y, X, Z numpy indexing)."""
    mhd_path = mhd_path.expanduser()
    text = mhd_path.read_text(encoding="utf-8")
    dim_field = _mhd_header_field(text, "DimSize")
    if dim_field is None:
        msg = f"DimSize missing in {mhd_path}"
        raise ValueError(msg)
    nx, ny, nz = (int(v) for v in dim_field.split())
    raw_name = _mhd_header_field(text, "ElementDataFile")
    if raw_name is None:
        msg = f"ElementDataFile missing in {mhd_path}"
        raise ValueError(msg)
    element_type = _mhd_header_field(text, "ElementType") or "MET_FLOAT"
    dtype = _mhd_element_dtype(element_type)
    raw_path = mhd_path.parent / raw_name
    if not raw_path.is_file():
        msg = f"RAW file missing: {raw_path}"
        raise FileNotFoundError(msg)
    expected_count = nx * ny * nz
    expected_bytes = expected_count * dtype.itemsize
    raw_bytes = raw_path.stat().st_size
    if raw_bytes != expected_bytes:
        msg = (
            f"RAW size mismatch for {raw_path}: {raw_bytes} bytes, "
            f"expected {expected_bytes} for DimSize {nx} {ny} {nz} and {element_type}"
        )
        raise ValueError(msg)
    flat = np.fromfile(raw_path, dtype=dtype)
    vol = flat.reshape((nx, ny, nz), order="C")
    return np.transpose(vol, (1, 0, 2)).astype(np.float32, copy=False)


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
