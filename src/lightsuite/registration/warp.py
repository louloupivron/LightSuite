"""Volume warping helpers for registration previews."""

from __future__ import annotations

import numpy as np
from scipy.ndimage import affine_transform


def pixel_affine_for_volume(matrix_4x4: np.ndarray) -> np.ndarray:
    """Convert 1-based voxel affine to 0-based scipy indexing."""
    shift_minus = np.eye(4)
    shift_minus[:3, 3] = -1.0
    shift_plus = np.eye(4)
    shift_plus[:3, 3] = 1.0
    return shift_plus @ matrix_4x4 @ shift_minus


def warp_volume_affine(
    volume: np.ndarray,
    matrix_4x4: np.ndarray,
    output_shape: tuple[int, int, int],
    *,
    order: int = 1,
    cval: float = 0.0,
) -> np.ndarray:
    """Warp volume with a 4x4 mapping source voxels -> target voxels (1-based MATLAB style)."""
    matrix_0 = pixel_affine_for_volume(np.asarray(matrix_4x4, dtype=float))
    inverse = np.linalg.inv(matrix_0)
    return affine_transform(
        volume,
        inverse[:3, :3],
        offset=inverse[:3, 3],
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
    """Warp atlas data into the sample grid (initializeRegistration.m / imwarp invert)."""
    return warp_volume_affine(
        atlas_volume,
        np.linalg.inv(np.asarray(sample_to_atlas, dtype=float)),
        sample_shape,
        order=order,
    )
