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
    "blank_image_alt",
    "canonical_view_slice",
    "canonical_view_transform",
    "layer_xy_from_slice_pixels",
    "map_display_pixels_to_slice",
    "map_slice_pixels_to_display",
    "match_points_atlas_slice",
    "match_points_sample_slice",
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


def blank_image_alt(slice_2d: np.ndarray, blank_flags: np.ndarray) -> np.ndarray:
    """Mask a slice half (``blankImage_alt.m``) so one region is annotated at a time.

    Flags outside ``{1, 2}`` mean "unknown" (for example a chooselist row recovered
    from MATLAB control points, where the flags are not stored) and leave the slice
    untouched rather than blanking an arbitrary half.
    """
    flags = np.asarray(blank_flags, dtype=int).ravel()
    if flags.size < 2 or not set(flags[:2].tolist()) <= {1, 2}:
        return slice_2d

    ny, nx = slice_2d.shape
    iy = np.arange(ny)
    ix = np.arange(nx)
    if flags[0] == 2:
        start = int(round(np.linspace(0, ny / 2, 2)[flags[0] - 1]))
        iy = np.arange(start, start + ny // 2)
    else:
        start = int(round(np.linspace(0, nx / 2, 2)[flags[1] - 1]))
        ix = np.arange(start, start + nx // 2)

    out = np.zeros_like(slice_2d)
    out[np.ix_(iy, ix)] = slice_2d[np.ix_(iy, ix)]
    return out


def match_points_sample_slice(
    volume: np.ndarray,
    chooserow: np.ndarray,
    atlas_provider: str,
) -> np.ndarray:
    """Sample panel slice: MATLAB half-mask applied before the canonical display remap."""
    raw = volume_index_to_image(volume, chooserow)
    masked = blank_image_alt(raw, np.asarray(chooserow, dtype=int).ravel()[2:4])
    return prepare_display_slice(masked, int(chooserow[1]), atlas_provider)


def match_points_atlas_slice(
    volume: np.ndarray,
    chooserow: np.ndarray,
    atlas_plane: int,
    atlas_provider: str,
) -> np.ndarray:
    """Atlas panel slice at ``atlas_plane`` along the chooselist cut axis."""
    row = np.asarray(chooserow, dtype=int).copy()
    row[0] = int(atlas_plane)
    return prepare_display_slice(
        volume_index_to_image(volume, row), int(row[1]), atlas_provider
    )


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
    """Map raw slice (row, col) to napari Points ``data`` as ``(row, col)``.

    Napari uses the same axis order as the image layer: axis 0 = rows, axis 1 = cols.
    """
    row_a = np.asarray(row, dtype=float)
    col_a = np.asarray(col, dtype=float)
    transform = canonical_view_transform(atlas_provider, cut_axis)
    if transform.rot90_k == 0 and not transform.flip_ud and not transform.flip_lr:
        return np.column_stack([row_a, col_a])

    disp_row, disp_col = map_slice_pixels_to_display(row_a, col_a, slice_shape, transform)
    return np.column_stack([disp_row, disp_col])


def slice_pixels_from_layer_xy(
    layer_rc: np.ndarray,
    slice_shape: tuple[int, int],
    cut_axis: int,
    atlas_provider: str,
) -> tuple[np.ndarray, np.ndarray]:
    """Inverse of :func:`layer_xy_from_slice_pixels` for napari click/drag events."""
    disp_row = np.asarray(layer_rc[:, 0], dtype=float)
    disp_col = np.asarray(layer_rc[:, 1], dtype=float)
    transform = canonical_view_transform(atlas_provider, cut_axis)
    if transform.rot90_k == 0 and not transform.flip_ud and not transform.flip_lr:
        return disp_row, disp_col

    return map_display_pixels_to_slice(disp_row, disp_col, slice_shape, transform)
