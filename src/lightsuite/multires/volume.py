"""Volume I/O for multiresolution pair manifests."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import SimpleITK as sitk
import tifffile

from lightsuite.multires.models import ManifestVolumeSpec


def tiff_shape(path: Path) -> tuple[int, int, int]:
    """Read array shape from TIFF header only (no full load)."""
    with tifffile.TiffFile(path) as tf:
        shape = tf.series[0].shape
    if len(shape) == 2:
        return (1, shape[0], shape[1])
    if len(shape) == 3:
        return shape
    if len(shape) == 4:
        return shape[1:] if shape[0] == 1 else shape[:3]
    msg = f"Unexpected TIFF shape {shape} for {path}"
    raise ValueError(msg)


def discover_volume_shape(volume_path: Path) -> tuple[int, int, int]:
    """Infer ZYX shape from a hyperstack TIFF or plane-per-file folder."""
    volume_path = volume_path.expanduser().resolve()
    if volume_path.is_file():
        return tiff_shape(volume_path)
    if volume_path.is_dir():
        planes = _sorted_plane_files(volume_path)
        if not planes:
            msg = f"No TIFF planes found in {volume_path}"
            raise FileNotFoundError(msg)
        ny, nx = tiff_shape(planes[0])[1:]
        return len(planes), ny, nx
    msg = f"Volume path is neither file nor directory: {volume_path}"
    raise FileNotFoundError(msg)


def empty_image_from_shape(shape_zyx: tuple[int, int, int]) -> sitk.Image:
    """Placeholder volume for geometry-only checks (no pixel I/O)."""
    z, y, x = shape_zyx
    arr = np.zeros((z, y, x), dtype=np.float32)
    return sitk.GetImageFromArray(arr)


def _sorted_plane_files(folder: Path) -> list[Path]:
    paths: list[Path] = []
    for pattern in ("*.tif", "*.tiff", "*.TIF", "*.TIFF"):
        paths.extend(folder.glob(pattern))
    return sorted(paths, key=lambda p: p.name)


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


def _load_tiff_hyperstack(path: Path) -> np.ndarray:
    path = path.expanduser().resolve()
    try:
        return np.asarray(tifffile.memmap(str(path)))
    except ValueError:
        return np.asarray(tifffile.imread(str(path)))


def _load_plane_per_file_stack(folder: Path) -> np.ndarray:
    planes = _sorted_plane_files(folder)
    if not planes:
        msg = f"No TIFF planes in {folder}"
        raise FileNotFoundError(msg)
    first = _normalize_tiff_array(np.asarray(tifffile.imread(str(planes[0]))), planes[0])
    if first.ndim != 3 or first.shape[0] != 1:
        ny, nx = first.shape[-2:]
        stack = np.zeros((len(planes), ny, nx), dtype=np.float32)
        for index, plane_path in enumerate(planes):
            plane = _normalize_tiff_array(np.asarray(tifffile.imread(str(plane_path))), plane_path)
            stack[index] = plane[0] if plane.ndim == 3 else plane
        return stack
    nz, ny, nx = len(planes), first.shape[1], first.shape[2]
    stack = np.zeros((nz, ny, nx), dtype=np.float32)
    for index, plane_path in enumerate(planes):
        plane = _normalize_tiff_array(np.asarray(tifffile.imread(str(plane_path))), plane_path)
        stack[index] = plane[0]
    return stack


def load_volume_array(volume_path: Path) -> np.ndarray:
    """Load a 3D volume as ZYX float32 from a TIFF file or plane-per-file folder."""
    volume_path = volume_path.expanduser().resolve()
    if volume_path.is_file():
        return _normalize_tiff_array(_load_tiff_hyperstack(volume_path), volume_path).astype(
            np.float32, copy=False
        )
    if volume_path.is_dir():
        return _load_plane_per_file_stack(volume_path).astype(np.float32, copy=False)
    msg = f"Volume path not found: {volume_path}"
    raise FileNotFoundError(msg)


def apply_manifest_geometry(image: sitk.Image, spec: ManifestVolumeSpec) -> sitk.Image:
    """Apply spacing, origin, and direction from a manifest volume spec."""
    image.SetSpacing(tuple(float(v) for v in spec.spacing_um))
    image.SetOrigin(tuple(float(v) for v in spec.origin_um))
    image.SetDirection(tuple(float(v) for v in spec.direction))
    return image


def load_manifest_volume(
    spec: ManifestVolumeSpec,
    *,
    load_pixels: bool = True,
    manifest_dir: Path | None = None,
) -> sitk.Image:
    """Load a manifest volume with geometry applied (microscope-agnostic)."""
    volume_path = Path(spec.volume_path).expanduser()
    if not volume_path.is_absolute() and manifest_dir is not None:
        volume_path = (manifest_dir / volume_path).resolve()

    if load_pixels:
        arr = load_volume_array(volume_path)
        if tuple(arr.shape) != tuple(spec.shape_zyx):
            msg = (
                f"Loaded shape ZYX {arr.shape} != manifest shape_zyx {spec.shape_zyx} "
                f"for {volume_path}"
            )
            raise ValueError(msg)
        image = sitk.GetImageFromArray(arr.astype(np.float32, copy=False))
    else:
        image = empty_image_from_shape(tuple(spec.shape_zyx))

    return apply_manifest_geometry(image, spec)


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


def volume_spec_from_image(label: str, volume_path: Path, image: sitk.Image) -> ManifestVolumeSpec:
    """Build a manifest volume spec from a loaded SimpleITK image."""
    arr = sitk.GetArrayFromImage(image)
    return ManifestVolumeSpec(
        volume_path=str(volume_path),
        shape_zyx=[int(v) for v in arr.shape],
        spacing_um=[float(s) for s in image.GetSpacing()],
        origin_um=[float(o) for o in image.GetOrigin()],
        direction=[float(v) for v in image.GetDirection()],
    )


def manifest_geometry_report(
    label: str,
    spec: ManifestVolumeSpec,
    *,
    manifest_dir: Path | None = None,
) -> dict[str, object]:
    """Build a diagnostic geometry report from manifest fields."""
    from lightsuite.multires.spec_geometry import manifest_geometry_report_from_spec

    _ = manifest_dir
    return manifest_geometry_report_from_spec(label, spec)


def _resolve_volume_path(spec: ManifestVolumeSpec, manifest_dir: Path | None) -> Path:
    volume_path = Path(spec.volume_path).expanduser()
    if not volume_path.is_absolute() and manifest_dir is not None:
        volume_path = (manifest_dir / volume_path).resolve()
    return volume_path


def load_manifest_xy_slice(
    spec: ManifestVolumeSpec,
    *,
    physical_um: tuple[float, float, float],
    manifest_dir: Path | None = None,
) -> np.ndarray:
    """Load one XY plane at a physical point without reading the full stack."""
    from lightsuite.multires.spec_geometry import physical_to_continuous_index_xyz

    volume_path = _resolve_volume_path(spec, manifest_dir)
    index_xyz = physical_to_continuous_index_xyz(spec, physical_um)
    ix, iy, iz = (int(round(float(v))) for v in index_xyz)
    nz, ny, nx = (int(v) for v in spec.shape_zyx)
    ix = int(np.clip(ix, 0, nx - 1))
    iy = int(np.clip(iy, 0, ny - 1))
    iz = int(np.clip(iz, 0, nz - 1))

    if volume_path.is_file():
        with tifffile.TiffFile(volume_path) as tf:
            series = tf.series[0]
            if series.shape[0] == nz:
                plane = series.asarray(key=iz)
            else:
                plane = tf.pages[iz].asarray()
        plane = _normalize_tiff_array(np.asarray(plane), volume_path)
        if plane.ndim == 3:
            return np.asarray(plane[0], dtype=np.float32)
        return np.asarray(plane, dtype=np.float32)

    if volume_path.is_dir():
        planes = _sorted_plane_files(volume_path)
        if iz >= len(planes):
            msg = f"Z index {iz} out of range for {len(planes)} planes in {volume_path}"
            raise IndexError(msg)
        plane = _normalize_tiff_array(np.asarray(tifffile.imread(str(planes[iz]))), planes[iz])
        if plane.ndim == 3:
            return np.asarray(plane[0], dtype=np.float32)
        return np.asarray(plane, dtype=np.float32)

    msg = f"Volume path not found: {volume_path}"
    raise FileNotFoundError(msg)
