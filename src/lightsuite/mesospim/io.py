"""TIFF I/O and axis remapping for mesoSPIM stacks."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import SimpleITK as sitk
import tifffile

from lightsuite.mesospim.config_models import MesospimGeometryConfig, MesospimTiffRemapConfig
from lightsuite.mesospim.geometry import apply_image_geometry

_MEMMAP_CACHE: dict[str, np.memmap] = {}


def tiff_shape(path: Path) -> tuple[int, int, int]:
    """Read array shape from TIFF header only (no full load)."""
    with tifffile.TiffFile(path) as tf:
        shape = tf.series[0].shape
    if len(shape) == 2:
        return (1, shape[0], shape[1])
    if len(shape) == 3:
        return shape  # (Z, Y, X) typical for mesoSPIM
    if len(shape) == 4:
        return shape[1:] if shape[0] == 1 else shape[:3]
    msg = f"Unexpected TIFF shape {shape} for {path}"
    raise ValueError(msg)


def empty_image_from_shape(shape_zyx: tuple[int, int, int]) -> sitk.Image:
    """Placeholder volume for geometry-only checks (no pixel I/O)."""
    z, y, x = shape_zyx
    arr = np.zeros((z, y, x), dtype=np.float32)
    return sitk.GetImageFromArray(arr)


def _rot90_k_for_path(
    path: Path,
    *,
    overview_path: Path,
    roi_path: Path,
    remap: MesospimTiffRemapConfig,
) -> int:
    resolved = path.resolve()
    if resolved == roi_path.resolve():
        return remap.rot90_k_roi
    if resolved == overview_path.resolve():
        return remap.rot90_k_overview
    return remap.rot90_k_overview


def remap_tiff_array_zyx(
    arr: np.ndarray,
    path: Path,
    *,
    overview_path: Path,
    roi_path: Path,
    remap: MesospimTiffRemapConfig,
) -> np.ndarray:
    """Remap raw tifffile (Z,Y,X) using configured axis transforms."""
    if arr.ndim != 3:
        msg = f"remap expects 3D ZYX, got {arr.shape}"
        raise ValueError(msg)

    out = arr
    if remap.reverse_z:
        out = out[::-1]

    k = _rot90_k_for_path(path, overview_path=overview_path, roi_path=roi_path, remap=remap)
    if k:
        if out.shape[1] != out.shape[2]:
            msg = f"rot90 k={k} requires square XY; got ZYX {out.shape} ({path.name})"
            raise ValueError(msg)
        out = np.rot90(out, k=k, axes=(1, 2))

    if remap.flip_row:
        out = out[:, ::-1, :]
    if remap.flip_col:
        out = out[:, :, ::-1]
    if remap.swap_xy:
        out = np.swapaxes(out, 1, 2)
    return np.ascontiguousarray(out)


def _apply_path_extra_flips(
    arr: np.ndarray,
    path: Path,
    *,
    overview_path: Path,
    roi_path: Path,
    remap: MesospimTiffRemapConfig,
) -> np.ndarray:
    out = arr
    resolved = path.resolve()
    if resolved == overview_path.resolve():
        if remap.extra_flip_row_overview:
            out = out[:, ::-1, :]
        if remap.extra_flip_col_overview:
            out = out[:, :, ::-1]
    elif resolved == roi_path.resolve():
        if remap.extra_flip_row_roi:
            out = out[:, ::-1, :]
        if remap.extra_flip_col_roi:
            out = out[:, :, ::-1]
    return np.ascontiguousarray(out)


def _normalize_tiff_array(arr: np.ndarray, path: Path) -> np.ndarray:
    if arr.size == 0:
        msg = f"Empty array from {path}"
        raise ValueError(msg)

    while arr.ndim > 3 and arr.shape[0] == 1:
        arr = arr[0]

    if arr.ndim == 2:
        arr = arr[np.newaxis, ...]
    elif arr.ndim == 4:
        shape = arr.shape
        if shape[-1] <= 4 and min(shape[1], shape[2]) > 8:
            arr = arr[..., 0]
        elif shape[1] <= 4 and min(shape[2], shape[3]) > 8:
            arr = arr[:, 0, :, :]
        else:
            msg = f"4D TIFF with ambiguous layout {shape} for {path}"
            raise ValueError(msg)

    if arr.ndim != 3:
        msg = f"Expected 3D array, got shape {arr.shape} for {path}"
        raise ValueError(msg)
    return arr


def _remapped_z_to_raw_page(remapped_z: int, nz: int, *, reverse_z: bool) -> int:
    if reverse_z:
        return nz - 1 - remapped_z
    return remapped_z


def _read_raw_plane_yx(path: Path, raw_z: int) -> np.ndarray:
    """Read one native YX plane (0-based Z) from a mesoSPIM TIFF stack."""
    path = path.expanduser().resolve()
    with tifffile.TiffFile(str(path)) as tf:
        if len(tf.pages) > 1 and tf.pages[0].ndim == 2:
            if raw_z < 0 or raw_z >= len(tf.pages):
                msg = f"Z index {raw_z} out of range for {len(tf.pages)} pages in {path.name}"
                raise IndexError(msg)
            return np.asarray(tf.pages[raw_z].asarray(), dtype=np.float32)

    key = str(path)
    vol = _MEMMAP_CACHE.get(key)
    if vol is None:
        vol = tifffile.memmap(key)
        _MEMMAP_CACHE[key] = vol
    if raw_z < 0 or raw_z >= vol.shape[0]:
        msg = f"Z index {raw_z} out of range for shape {vol.shape} in {path.name}"
        raise IndexError(msg)
    return np.asarray(vol[raw_z], dtype=np.float32)


def read_tiff_xy_slice_at_z_index(
    path: Path,
    z_index: int,
    *,
    overview_path: Path,
    roi_path: Path,
    remap: MesospimTiffRemapConfig,
) -> np.ndarray:
    """Load one remapped XY plane by stack Z index without reading the full stack."""
    path = path.expanduser().resolve()
    shape_zyx = tiff_shape(path)
    nz = shape_zyx[0]
    remapped_z = int(np.clip(int(z_index), 0, nz - 1))
    raw_page = _remapped_z_to_raw_page(remapped_z, nz, reverse_z=remap.reverse_z)

    plane = _read_raw_plane_yx(path, raw_page)
    if plane.ndim != 2:
        msg = f"Expected 2D TIFF page, got shape {plane.shape} for {path}"
        raise ValueError(msg)

    slab = remap_tiff_array_zyx(
        plane[np.newaxis, ...],
        path,
        overview_path=overview_path,
        roi_path=roi_path,
        remap=remap,
    )
    slab = _apply_path_extra_flips(
        slab,
        path,
        overview_path=overview_path,
        roi_path=roi_path,
        remap=remap,
    )
    return np.asarray(slab[0], dtype=np.float32)


def read_tiff_xy_slice_at_physical_um(
    path: Path,
    *,
    meta: dict[str, float | int | str],
    geometry: MesospimGeometryConfig,
    overview_path: Path,
    roi_path: Path,
    remap: MesospimTiffRemapConfig,
    cx_um: float,
    cy_um: float,
    cz_um: float,
) -> np.ndarray:
    """Load one XY plane at a physical point without reading the full stack."""
    path = path.expanduser().resolve()
    shape_zyx = tiff_shape(path)
    nz = shape_zyx[0]

    img = empty_image_from_shape(shape_zyx)
    apply_image_geometry(img, meta, geometry)
    idx = img.TransformPhysicalPointToContinuousIndex((float(cx_um), float(cy_um), float(cz_um)))
    remapped_z = int(np.clip(round(float(idx[2])), 0, nz - 1))
    raw_page = _remapped_z_to_raw_page(remapped_z, nz, reverse_z=remap.reverse_z)

    plane = _read_raw_plane_yx(path, raw_page)
    if plane.ndim != 2:
        msg = f"Expected 2D TIFF page, got shape {plane.shape} for {path}"
        raise ValueError(msg)

    slab = remap_tiff_array_zyx(
        plane[np.newaxis, ...],
        path,
        overview_path=overview_path,
        roi_path=roi_path,
        remap=remap,
    )
    slab = _apply_path_extra_flips(
        slab,
        path,
        overview_path=overview_path,
        roi_path=roi_path,
        remap=remap,
    )
    return np.asarray(slab[0], dtype=np.float32)


def read_tiff_as_float(
    path: Path,
    *,
    overview_path: Path,
    roi_path: Path,
    remap: MesospimTiffRemapConfig,
    validate_shape: bool = True,
) -> sitk.Image:
    """Load 3D stack as float32 with configured axis remapping."""
    path = path.expanduser().resolve()
    arr = np.asarray(tifffile.imread(str(path)))
    arr = _normalize_tiff_array(arr, path)
    arr = remap_tiff_array_zyx(
        arr,
        path,
        overview_path=overview_path,
        roi_path=roi_path,
        remap=remap,
    )
    arr = _apply_path_extra_flips(
        arr,
        path,
        overview_path=overview_path,
        roi_path=roi_path,
        remap=remap,
    )

    if validate_shape:
        expected_zyx = tiff_shape(path)
        if (arr.shape[0], arr.shape[1], arr.shape[2]) != expected_zyx:
            msg = (
                f"Array shape ZYX {arr.shape} != tiff_shape {expected_zyx} for {path}. "
                "Fix tiff_remap settings for this layout."
            )
            raise ValueError(msg)

    image = sitk.GetImageFromArray(arr.astype(np.float32, copy=False))
    if validate_shape:
        ref = empty_image_from_shape(tiff_shape(path))
        if image.GetSize() != ref.GetSize():
            msg = (
                f"Internal check failed: loaded GetSize {image.GetSize()} != "
                f"placeholder {ref.GetSize()}"
            )
            raise RuntimeError(msg)
    return image


def _load_tiff_raw_array(path: Path) -> np.ndarray:
    """Load a TIFF volume, using memmap when possible and falling back to a full read."""
    path = path.expanduser().resolve()
    try:
        return np.asarray(tifffile.memmap(str(path)))
    except ValueError:
        return np.asarray(tifffile.imread(str(path)))


def load_tiff_zyx_volume(
    path: Path,
    *,
    overview_path: Path,
    roi_path: Path,
    remap: MesospimTiffRemapConfig,
) -> np.ndarray:
    """Load a mesoSPIM stack and apply configured axis remapping."""
    path = path.expanduser().resolve()
    vol = _load_tiff_raw_array(path)
    vol = _normalize_tiff_array(vol, path)
    vol = remap_tiff_array_zyx(
        vol,
        path,
        overview_path=overview_path,
        roi_path=roi_path,
        remap=remap,
    )
    return _apply_path_extra_flips(
        vol,
        path,
        overview_path=overview_path,
        roi_path=roi_path,
        remap=remap,
    )


def load_tiff_zyx_memmap(
    path: Path,
    *,
    overview_path: Path,
    roi_path: Path,
    remap: MesospimTiffRemapConfig,
) -> np.ndarray:
    """Load a mesoSPIM stack (memmap when possible) with axis remapping."""
    return load_tiff_zyx_volume(
        path,
        overview_path=overview_path,
        roi_path=roi_path,
        remap=remap,
    )


def load_registered_canvas_zyx(path: Path) -> np.ndarray:
    """Load a registered full-overview canvas TIFF as ZYX float32."""
    vol = _load_tiff_raw_array(path)
    vol = _normalize_tiff_array(vol, path)
    return np.asarray(vol, dtype=np.float32)


def write_sitk_hyperstack_tiff(path: Path, image: sitk.Image) -> None:
    """Write Z,Y,X float volume as ImageJ-style hyperstack."""
    path = path.expanduser()
    path.parent.mkdir(parents=True, exist_ok=True)
    arr = sitk.GetArrayFromImage(image)
    sx, sy, sz = (float(s) for s in image.GetSpacing())
    tifffile.imwrite(
        path,
        arr.astype(np.float32, copy=False),
        imagej=True,
        resolution=(1.0 / sx, 1.0 / sy),
        metadata={"spacing": sz, "unit": "um"},
        compression="zlib",
    )


def sitk_to_itk(image: sitk.Image):
    """Convert SimpleITK image to ITK float image, preserving physical space."""
    import itk

    itk_img = itk.GetImageFromArray(sitk.GetArrayFromImage(image).astype(np.float32))
    itk_img.SetSpacing(tuple(float(s) for s in image.GetSpacing()))
    itk_img.SetOrigin(tuple(float(o) for o in image.GetOrigin()))
    itk_img.SetDirection(np.array(image.GetDirection(), dtype=float).reshape(3, 3))
    return itk_img
