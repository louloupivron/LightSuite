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


def discover_volume_shape(
    volume_path: Path,
    *,
    expected_planes: int | None = None,
) -> tuple[int, int, int]:
    """Infer ZYX shape from a hyperstack TIFF or plane-per-file folder."""
    volume_path = volume_path.expanduser().resolve()
    if volume_path.is_file():
        return tiff_shape(volume_path)
    if volume_path.is_dir():
        planes = _sorted_plane_files(volume_path, expected_planes=expected_planes)
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


def _sorted_plane_files(folder: Path, *, expected_planes: int | None = None) -> list[Path]:
    """List plane TIFFs once per file (Windows glob is case-insensitive)."""
    paths: list[Path] = []
    for entry in folder.iterdir():
        if not entry.is_file():
            continue
        if entry.suffix.lower() not in {".tif", ".tiff"}:
            continue
        paths.append(entry)
    paths.sort(key=lambda p: p.name)
    return _resolve_smartspim_plane_files(paths, expected_planes=expected_planes)


def _resolve_smartspim_plane_files(
    paths: list[Path],
    *,
    expected_planes: int | None = None,
) -> list[Path]:
    """Keep one SmartSPIM laser channel when ``All_Channels`` holds interleaved planes."""
    if not paths:
        return paths

    from collections import defaultdict

    from lightsuite.io.smartspim_channels import parse_smartspim_channel_plane

    by_channel: dict[int, list[Path]] = defaultdict(list)
    for path in paths:
        parsed = parse_smartspim_channel_plane(path)
        if parsed is not None:
            by_channel[parsed[0]].append(path)

    if not by_channel:
        return paths

    for channel in sorted(by_channel):
        ch_paths = sorted(by_channel[channel], key=lambda p: p.name.lower())
        if expected_planes is not None and len(ch_paths) == expected_planes:
            return ch_paths

    if expected_planes is not None and len(paths) >= 2 * expected_planes and len(by_channel) >= 2:
        ch_paths = sorted(by_channel[min(by_channel)], key=lambda p: p.name.lower())
        if len(ch_paths) == expected_planes:
            return ch_paths

    if expected_planes is None and len(by_channel) >= 2:
        counts = {channel: len(by_channel[channel]) for channel in by_channel}
        if len(set(counts.values())) == 1:
            return sorted(by_channel[min(by_channel)], key=lambda p: p.name.lower())

    return paths


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
    from lightsuite.registration.volume import load_tiff_volume_zyx

    return np.asarray(load_tiff_volume_zyx(path))


_HYPERSTACK_MEMMAP: dict[str, np.ndarray] = {}


def _single_page_hyperstack_volume(volume_path: Path, shape_zyx: tuple[int, int, int]) -> np.ndarray:
    """Return a ZYX view for a one-page ImageJ hyperstack TIFF (mesoSPIM-style)."""
    key = str(volume_path.expanduser().resolve())
    cached = _HYPERSTACK_MEMMAP.get(key)
    if cached is not None:
        return cached
    try:
        vol = tifffile.memmap(key)
    except ValueError:
        vol = _load_tiff_hyperstack(volume_path)
    vol = _normalize_tiff_array(np.asarray(vol), volume_path)
    if tuple(vol.shape) != tuple(shape_zyx):
        msg = (
            f"Hyperstack shape ZYX {vol.shape} != manifest shape_zyx {shape_zyx} "
            f"for {volume_path}"
        )
        raise ValueError(msg)
    _HYPERSTACK_MEMMAP[key] = vol
    return vol


def _read_tiff_xy_plane(
    volume_path: Path,
    *,
    z_index: int,
    shape_zyx: tuple[int, int, int],
) -> np.ndarray:
    """Read one native Z plane from a hyperstack TIFF or page-per-slice stack."""
    nz, _ny, _nx = shape_zyx
    iz = int(np.clip(int(z_index), 0, nz - 1))
    with tifffile.TiffFile(volume_path) as tf:
        series = tf.series[0]
        if len(tf.pages) == 1 and series.ndim == 3 and int(series.shape[0]) == nz:
            return np.asarray(_single_page_hyperstack_volume(volume_path, shape_zyx)[iz])
        if series.ndim == 3 and int(series.shape[0]) == nz and len(tf.pages) > 1:
            plane = series.asarray(key=iz)
        else:
            if iz >= len(tf.pages):
                msg = f"Z index {iz} out of range for {len(tf.pages)} pages in {volume_path}"
                raise IndexError(msg)
            plane = tf.pages[iz].asarray()
    plane = _normalize_tiff_array(np.asarray(plane), volume_path)
    if plane.ndim == 3:
        return np.asarray(plane[0], dtype=np.float32)
    return np.asarray(plane, dtype=np.float32)


def _load_plane_per_file_stack(folder: Path, *, expected_planes: int | None = None) -> np.ndarray:
    planes = _sorted_plane_files(folder, expected_planes=expected_planes)
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


def load_volume_array(volume_path: Path, *, expected_planes: int | None = None) -> np.ndarray:
    """Load a 3D volume as ZYX float32 from a TIFF file or plane-per-file folder."""
    volume_path = volume_path.expanduser().resolve()
    if volume_path.is_file():
        return _normalize_tiff_array(_load_tiff_hyperstack(volume_path), volume_path).astype(
            np.float32, copy=False
        )
    if volume_path.is_dir():
        return _load_plane_per_file_stack(volume_path, expected_planes=expected_planes).astype(
            np.float32, copy=False
        )
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
        arr = load_volume_array(volume_path, expected_planes=int(spec.shape_zyx[0]))
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


def volume_spec_from_geometry(
    volume_path: Path,
    shape_zyx: tuple[int, int, int],
    spacing_um: tuple[float, float, float],
    origin_um: tuple[float, float, float],
    direction: list[float] | tuple[float, ...],
) -> ManifestVolumeSpec:
    """Build a manifest volume spec without allocating voxel data."""
    return ManifestVolumeSpec(
        volume_path=str(volume_path),
        shape_zyx=[int(v) for v in shape_zyx],
        spacing_um=[float(v) for v in spacing_um],
        origin_um=[float(v) for v in origin_um],
        direction=[float(v) for v in direction],
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
        return _read_tiff_xy_plane(volume_path, z_index=iz, shape_zyx=(nz, ny, nx))

    if volume_path.is_dir():
        planes = _sorted_plane_files(volume_path, expected_planes=nz)
        if iz >= len(planes):
            msg = f"Z index {iz} out of range for {len(planes)} planes in {volume_path}"
            raise IndexError(msg)
        plane = _normalize_tiff_array(np.asarray(tifffile.imread(str(planes[iz]))), planes[iz])
        if plane.ndim == 3:
            return np.asarray(plane[0], dtype=np.float32)
        return np.asarray(plane, dtype=np.float32)

    msg = f"Volume path not found: {volume_path}"
    raise FileNotFoundError(msg)


def load_manifest_xy_plane_at_z_index(
    spec: ManifestVolumeSpec,
    z_index: int,
    *,
    manifest_dir: Path | None = None,
) -> np.ndarray:
    """Load one full XY plane by Z index without reading the full stack."""
    _nz, ny, nx = (int(v) for v in spec.shape_zyx)
    return load_manifest_xy_crop(
        spec,
        z_index=z_index,
        start_xyz=[0, 0, 0],
        crop_size_xyz=[nx, ny, 1],
        manifest_dir=manifest_dir,
    )


def load_manifest_xy_crop(
    spec: ManifestVolumeSpec,
    *,
    z_index: int,
    start_xyz: list[int],
    crop_size_xyz: list[int],
    manifest_dir: Path | None = None,
) -> np.ndarray:
    """Load an XY crop from one Z plane without reading the full stack."""
    volume_path = _resolve_volume_path(spec, manifest_dir)
    ix0, iy0, _iz0 = start_xyz
    sx, sy, _sz = crop_size_xyz
    ix1 = ix0 + sx
    iy1 = iy0 + sy
    nz, ny, nx = (int(v) for v in spec.shape_zyx)
    iz = int(np.clip(z_index, 0, nz - 1))
    ix0 = int(np.clip(ix0, 0, nx - 1))
    iy0 = int(np.clip(iy0, 0, ny - 1))
    ix1 = int(np.clip(ix1, 0, nx))
    iy1 = int(np.clip(iy1, 0, ny))

    if volume_path.is_file():
        plane = _read_tiff_xy_plane(volume_path, z_index=iz, shape_zyx=(nz, ny, nx))
        return np.asarray(plane[iy0:iy1, ix0:ix1], dtype=np.float32)

    if volume_path.is_dir():
        planes = _sorted_plane_files(volume_path, expected_planes=nz)
        if iz >= len(planes):
            msg = f"Z index {iz} out of range for {len(planes)} planes in {volume_path}"
            raise IndexError(msg)
        plane = _normalize_tiff_array(np.asarray(tifffile.imread(str(planes[iz]))), planes[iz])
        if plane.ndim == 3:
            plane = plane[0]
        return np.asarray(plane[iy0:iy1, ix0:ix1], dtype=np.float32)

    msg = f"Volume path not found: {volume_path}"
    raise FileNotFoundError(msg)


def load_manifest_xyz_crop(
    spec: ManifestVolumeSpec,
    *,
    start_xyz: list[int],
    crop_size_xyz: list[int],
    manifest_dir: Path | None = None,
) -> sitk.Image:
    """Load a 3D XYZ crop plane-by-plane without materialising the full stack."""
    from lightsuite.multires.spec_geometry import index_xyz_to_physical

    ix0, iy0, iz0 = (int(v) for v in start_xyz)
    sx, sy, sz = (int(v) for v in crop_size_xyz)
    if min(sx, sy, sz) <= 0:
        msg = f"Invalid crop size {crop_size_xyz}"
        raise ValueError(msg)

    stack = np.empty((sz, sy, sx), dtype=np.float32)
    for dz in range(sz):
        stack[dz] = load_manifest_xy_crop(
            spec,
            z_index=iz0 + dz,
            start_xyz=[ix0, iy0, iz0],
            crop_size_xyz=[sx, sy, 1],
            manifest_dir=manifest_dir,
        )

    image = sitk.GetImageFromArray(stack)
    origin = index_xyz_to_physical(spec, (float(ix0), float(iy0), float(iz0)))
    image.SetSpacing(tuple(float(v) for v in spec.spacing_um))
    image.SetOrigin(tuple(float(v) for v in origin))
    image.SetDirection(tuple(float(v) for v in spec.direction))
    return image


def stream_resample_to_reference(
    moving_spec: ManifestVolumeSpec,
    reference: sitk.Image,
    *,
    manifest_dir: Path | None = None,
    reference_to_moving: np.ndarray | None = None,
    z_chunk: int | None = None,
    max_slab_bytes: int = 1_500_000_000,
) -> sitk.Image:
    """Resample a manifest volume onto ``reference`` without loading the full moving stack.

    For each Z chunk of the reference grid, only the moving voxels that can
    contribute to that chunk are loaded from disk, then SimpleITK resamples
    the slab onto the chunk. Chunk size is chosen so each moving slab stays
    near ``max_slab_bytes`` (default 1.5 GB).
    """
    from lightsuite.multires.geometry import physical_corners
    from lightsuite.multires.spec_geometry import crop_index_range_from_physical_box

    ref_size = reference.GetSize()  # x, y, z
    sx, sy, sz = (int(v) for v in ref_size)
    out = np.zeros((sz, sy, sx), dtype=np.float32)
    if reference_to_moving is None:
        transform: sitk.Transform = sitk.Transform(3, sitk.sitkIdentity)
        moving_from_ref = None
    else:
        moving_from_ref = np.asarray(reference_to_moving, dtype=float)
        transform = _sitk_affine_from_matrix(moving_from_ref)

    if z_chunk is None:
        # Estimate moving XY footprint from the full reference physical extent.
        ref_min = np.array(reference.TransformIndexToPhysicalPoint((0, 0, 0)), dtype=float)
        ref_max = np.array(
            reference.TransformIndexToPhysicalPoint((sx - 1, sy - 1, max(sz - 1, 0))),
            dtype=float,
        )
        phys_min = np.minimum(ref_min, ref_max)
        phys_max = np.maximum(ref_min, ref_max)
        if moving_from_ref is not None:
            corners = physical_corners(phys_min, phys_max)
            ones = np.ones((corners.shape[0], 1), dtype=float)
            moving_corners = (np.hstack([corners, ones]) @ moving_from_ref.T)[:, :3]
            phys_min = moving_corners.min(axis=0)
            phys_max = moving_corners.max(axis=0)
        _start, est_size = crop_index_range_from_physical_box(moving_spec, phys_min, phys_max)
        est_yx = max(1, int(est_size[0]) * int(est_size[1]))
        bytes_per_plane = est_yx * 4
        max_moving_planes = max(2, int(max_slab_bytes / max(bytes_per_plane, 1)))
        ref_dz = float(reference.GetSpacing()[2])
        mov_dz = float(moving_spec.spacing_um[2])
        planes_per_ref = max(1.0, ref_dz / max(mov_dz, 1e-6))
        z_chunk = max(1, int(max_moving_planes / planes_per_ref))

    chunk = max(1, int(z_chunk))
    for z0 in range(0, sz, chunk):
        z1 = min(sz, z0 + chunk)
        chunk_size = [sx, sy, z1 - z0]
        chunk_ref = sitk.RegionOfInterest(reference, chunk_size, [0, 0, z0])

        chunk_min = np.array(chunk_ref.TransformIndexToPhysicalPoint((0, 0, 0)), dtype=float)
        chunk_max = np.array(
            chunk_ref.TransformIndexToPhysicalPoint((sx - 1, sy - 1, z1 - z0 - 1)),
            dtype=float,
        )
        phys_min = np.minimum(chunk_min, chunk_max)
        phys_max = np.maximum(chunk_min, chunk_max)

        if moving_from_ref is not None:
            corners = physical_corners(phys_min, phys_max)
            ones = np.ones((corners.shape[0], 1), dtype=float)
            moving_corners = (np.hstack([corners, ones]) @ moving_from_ref.T)[:, :3]
            phys_min = moving_corners.min(axis=0)
            phys_max = moving_corners.max(axis=0)

        # Pad by one voxel so linear interpolation has neighbours at the edges.
        pad = np.asarray(moving_spec.spacing_um, dtype=float)
        phys_min = phys_min - pad
        phys_max = phys_max + pad
        start, crop_size = crop_index_range_from_physical_box(moving_spec, phys_min, phys_max)
        if min(crop_size) <= 0:
            continue

        moving_slab = load_manifest_xyz_crop(
            moving_spec,
            start_xyz=start,
            crop_size_xyz=crop_size,
            manifest_dir=manifest_dir,
        )
        resampled = sitk.Resample(
            moving_slab,
            chunk_ref,
            transform,
            sitk.sitkLinear,
            0.0,
            sitk.sitkFloat32,
        )
        out[z0:z1] = sitk.GetArrayFromImage(resampled)
        del moving_slab, resampled

    result = sitk.GetImageFromArray(out)
    result.CopyInformation(reference)
    return result


def _sitk_affine_from_matrix(matrix: np.ndarray) -> sitk.AffineTransform:
    transform = sitk.AffineTransform(3)
    transform.SetMatrix(matrix[:3, :3].reshape(-1).tolist())
    transform.SetTranslation(matrix[:3, 3].tolist())
    return transform


def write_embedded_crop_canvas(
    overview_spec: ManifestVolumeSpec,
    crop: sitk.Image,
    crop_start_index: list[int],
    output_path: Path,
    *,
    dtype=np.float32,
) -> None:
    """Write a full-overview canvas with ``crop`` pasted in, plane-by-plane.

    Avoids allocating a buffer the size of the overview (tens of GB). Pass
    ``dtype=np.uint8`` for masks to keep the canvas a quarter of the size.
    """
    output_path = output_path.expanduser()
    output_path.parent.mkdir(parents=True, exist_ok=True)

    nz, ny, nx = (int(v) for v in overview_spec.shape_zyx)
    cx, cy, cz = (int(v) for v in crop.GetSize())
    ix0, iy0, iz0 = (int(v) for v in crop_start_index)
    crop_arr = np.asarray(sitk.GetArrayViewFromImage(crop), dtype=dtype)
    sx, sy, sz = (float(v) for v in overview_spec.spacing_um)

    # Multi-page BigTIFF; ImageJ hyperstack metadata is set on the first page.
    with tifffile.TiffWriter(output_path, bigtiff=True) as tif:
        for iz in range(nz):
            plane = np.zeros((ny, nx), dtype=dtype)
            if iz0 <= iz < iz0 + cz:
                y1 = min(iy0 + cy, ny)
                x1 = min(ix0 + cx, nx)
                yy0 = max(iy0, 0)
                xx0 = max(ix0, 0)
                src_y0 = yy0 - iy0
                src_x0 = xx0 - ix0
                plane[yy0:y1, xx0:x1] = crop_arr[
                    iz - iz0,
                    src_y0 : src_y0 + (y1 - yy0),
                    src_x0 : src_x0 + (x1 - xx0),
                ]
            metadata = None
            if iz == 0:
                metadata = {
                    "spacing": sz,
                    "unit": "um",
                }
            tif.write(
                plane,
                compression="zlib",
                photometric="minisblack",
                metadata=metadata,
                resolution=(1.0 / sx, 1.0 / sy),
            )
