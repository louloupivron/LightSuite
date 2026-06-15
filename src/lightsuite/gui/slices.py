"""Extract 2D slices from 3D volumes (volumeIdtoImage.m)."""

from __future__ import annotations

import numpy as np

from lightsuite.atlas.display import (
    SliceDisplayTransform,
    apply_slice_display_transform,
    canonical_view_slice,
    canonical_view_transform,
    map_display_pixels_to_slice,
    map_slice_pixels_to_display,
)

# Re-export transform types for tests and historical imports.
__all__ = [
    "SliceDisplayTransform",
    "apply_slice_display_transform",
    "canonical_view_slice",
    "canonical_view_transform",
    "layer_xy_from_slice_pixels",
    "map_display_pixels_to_slice",
    "map_slice_pixels_to_display",
    "prepare_display_slice",
    "slice_pixels_from_layer_xy",
    "volume_index_to_image",
]


def volume_index_to_image(volume: np.ndarray, chooserow: np.ndarray) -> np.ndarray:
    """Extract a 2D slice from a 3D volume given chooselist row."""
    slice_index = int(chooserow[0])
    axis = int(chooserow[1])  # 1-based MATLAB dim
    slices = [slice(None)] * 3
    slices[axis - 1] = slice_index - 1
    image = volume[tuple(slices)]
    return np.squeeze(image)


def prepare_display_slice(
    slice_2d: np.ndarray,
    cut_axis: int,
    atlas_provider: str,
) -> np.ndarray:
    """Map a native atlas-order slice to the canonical QC view for ``atlas_provider``."""
    return canonical_view_slice(slice_2d, atlas_provider=atlas_provider, cut_axis=cut_axis)


def layer_xy_from_slice_pixels(
    row: np.ndarray | float,
    col: np.ndarray | float,
    slice_shape: tuple[int, int],
    cut_axis: int,
    atlas_provider: str,
) -> np.ndarray:
    """Map slice (row, col) pixel coords to napari layer (x, y) after canonical display."""
    transform = canonical_view_transform(atlas_provider, cut_axis)
    if transform.rot90_k == 0 and not transform.flip_ud and not transform.flip_lr:
        return np.column_stack([np.asarray(col, dtype=float), np.asarray(row, dtype=float)])

    disp_row, disp_col = map_slice_pixels_to_display(row, col, slice_shape, transform)
    return np.column_stack([disp_col, disp_row])


def slice_pixels_from_layer_xy(
    xy: np.ndarray,
    slice_shape: tuple[int, int],
    cut_axis: int,
    atlas_provider: str,
) -> tuple[np.ndarray, np.ndarray]:
    """Inverse of :func:`layer_xy_from_slice_pixels` for click/drag events."""
    disp_col = np.asarray(xy[:, 0], dtype=float)
    disp_row = np.asarray(xy[:, 1], dtype=float)
    transform = canonical_view_transform(atlas_provider, cut_axis)
    if transform.rot90_k == 0 and not transform.flip_ud and not transform.flip_lr:
        return disp_row, disp_col

    return map_display_pixels_to_slice(disp_row, disp_col, slice_shape, transform)
