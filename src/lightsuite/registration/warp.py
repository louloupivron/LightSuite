"""Volume warping helpers matching MATLAB imwarp conventions."""

from __future__ import annotations

import numpy as np
from scipy.ndimage import affine_transform, map_coordinates


def _transform_points_internal(points: np.ndarray, matrix: np.ndarray) -> np.ndarray:
    if points.shape[0] == 0:
        return points
    hom = np.column_stack([points, np.ones(points.shape[0])])
    return (hom @ np.asarray(matrix, dtype=float).T)[:, :3]


def _array_indices_to_xyz(indices: np.ndarray) -> np.ndarray:
    """Convert array indices (Y, X, Z) to point cloud XYZ = (col, row, depth)."""
    return np.stack([indices[..., 1], indices[..., 0], indices[..., 2]], axis=-1)


def _xyz_to_array_coordinates(xyz: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Convert XYZ points back to ``map_coordinates`` input (axis0=Y, axis1=X, axis2=Z)."""
    return xyz[..., 1], xyz[..., 0], xyz[..., 2]

_SWAP_XY = np.array(
    [
        [0.0, 1.0, 0.0, 0.0],
        [1.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 1.0],
    ],
    dtype=float,
)


def affinetform_rows_to_internal(matrix_4x4: np.ndarray) -> np.ndarray:
    """Convert ``points @ A + t`` affinetform to internal ``_transform_points`` layout."""
    matrix = np.asarray(matrix_4x4, dtype=float)
    internal = np.eye(4, dtype=float)
    internal[:3, :3] = matrix[:3, :3].T
    internal[:3, 3] = matrix[:3, 3]
    return internal


def transform_points_affinetform(points: np.ndarray, affinetform: np.ndarray) -> np.ndarray:
    """Apply affinetform3d-style ``points @ A + t`` on 0-based XYZ points."""
    matrix = np.asarray(affinetform, dtype=float)
    return np.asarray(points, dtype=float) @ matrix[:3, :3] + matrix[:3, 3]


def matlab_voxel_affine_from_icp(matrix: np.ndarray) -> np.ndarray:
    """Convert Open3D ICP transform (0-based XYZ) to MATLAB 1-based affinetform3d."""
    shift_minus = np.eye(4)
    shift_minus[:3, 3] = -1.0
    shift_plus = np.eye(4)
    shift_plus[:3, 3] = 1.0
    return shift_minus @ np.asarray(matrix, dtype=float) @ shift_plus


def matlab_voxel_affine_to_icp(matrix: np.ndarray) -> np.ndarray:
    """Inverse of :func:`matlab_voxel_affine_from_icp` for 0-based point clouds."""
    shift_minus = np.eye(4)
    shift_minus[:3, 3] = -1.0
    shift_plus = np.eye(4)
    shift_plus[:3, 3] = 1.0
    return shift_plus @ np.asarray(matrix, dtype=float) @ shift_minus


def pixel_affine_for_volume(matrix_4x4: np.ndarray) -> np.ndarray:
    """Convert 1-based voxel affine to 0-based scipy indexing."""
    shift_minus = np.eye(4)
    shift_minus[:3, 3] = -1.0
    shift_plus = np.eye(4)
    shift_plus[:3, 3] = 1.0
    return shift_plus @ matrix_4x4 @ shift_minus


def swap_xy_transform(matrix_4x4: np.ndarray) -> np.ndarray:
    """Convert between XYZ point coords (col, row, depth) and array axes (Y, X, Z)."""
    matrix = np.asarray(matrix_4x4, dtype=float)
    return _SWAP_XY @ matrix @ _SWAP_XY


def imwarp_volume(
    volume: np.ndarray,
    tform: np.ndarray,
    output_shape: tuple[int, int, int],
    *,
    order: int = 1,
    cval: float = 0.0,
    point_coords: str = "array",
) -> np.ndarray:
    """Warp like MATLAB imwarp(volume, ref_in, tform, OutputView=ref_out).

    ``tform`` is the transform object passed to imwarp. imwarp applies inv(tform)
    when resampling. Use ``point_coords='xyz'`` for transforms fit on [X, Y, Z]
    point clouds (original_trans), and ``point_coords='array'`` for transforms
    already in volume axis order (affine control-point fits after [2,1,3] swap).
    """
    tform = np.asarray(tform, dtype=float)
    effective = np.linalg.inv(tform)
    if point_coords == "xyz":
        effective = swap_xy_transform(effective)
    elif point_coords != "array":
        msg = f"point_coords must be 'xyz' or 'array', got {point_coords!r}"
        raise ValueError(msg)

    matrix_0 = pixel_affine_for_volume(effective)
    return affine_transform(
        volume,
        matrix_0[:3, :3],
        offset=matrix_0[:3, 3],
        output_shape=output_shape,
        order=order,
        mode="constant",
        cval=cval,
    )


def warp_atlas_to_sample(
    atlas_volume: np.ndarray,
    sample_to_atlas: np.ndarray,
    sample_shape: tuple[int, int, int],
    *,
    order: int = 0,
) -> np.ndarray:
    """Warp atlas annotation/template onto the sample grid (initializeRegistration.m)."""
    internal = matlab_voxel_affine_to_icp(np.asarray(sample_to_atlas, dtype=float))
    grid_y, grid_x, grid_z = np.indices(sample_shape, dtype=float)
    sample_xyz = _array_indices_to_xyz(np.stack([grid_y, grid_x, grid_z], axis=-1))
    atlas_xyz = _transform_points_internal(
        sample_xyz.reshape(-1, 3),
        internal,
    ).reshape(sample_shape + (3,))
    coords = _xyz_to_array_coordinates(atlas_xyz)
    warped = map_coordinates(
        np.asarray(atlas_volume, dtype=float),
        coords,
        order=order,
        mode="constant",
        cval=0.0,
        prefilter=order > 1,
    )
    return warped.astype(atlas_volume.dtype, copy=False)


def warp_sample_to_atlas(
    sample_volume: np.ndarray,
    sample_to_atlas: np.ndarray,
    atlas_shape: tuple[int, int, int],
    *,
    order: int = 1,
) -> np.ndarray:
    """Warp sample volume onto the downsampled atlas grid (matchControlPoints.m)."""
    return imwarp_volume(
        sample_volume,
        sample_to_atlas,
        atlas_shape,
        order=order,
        point_coords="xyz",
    )


def warp_volume_affine(
    volume: np.ndarray,
    matrix_4x4: np.ndarray,
    output_shape: tuple[int, int, int],
    *,
    order: int = 1,
    cval: float = 0.0,
    point_coords: str = "array",
) -> np.ndarray:
    """Backward-compatible alias for atlas-to-sample warps used after affine registration."""
    return imwarp_volume(
        volume,
        matrix_4x4,
        output_shape,
        order=order,
        cval=cval,
        point_coords=point_coords,
    )
