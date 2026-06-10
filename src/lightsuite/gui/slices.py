"""Extract 2D slices from 3D volumes (volumeIdtoImage.m)."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


def volume_index_to_image(volume: np.ndarray, chooserow: np.ndarray) -> np.ndarray:
    """Extract a 2D slice from a 3D volume given chooselist row."""
    slice_index = int(chooserow[0])
    axis = int(chooserow[1])  # 1-based MATLAB dim
    slices = [slice(None)] * 3
    slices[axis - 1] = slice_index - 1
    image = volume[tuple(slices)]
    return np.squeeze(image)


def needs_allen_display_reorientation(permvec: list[int] | None) -> bool:
    """Return True when 2D slices need Allen-style rotation for display.

    When the atlas Z axis is flipped in ``brain_orientation.txt`` / ``permvec`` (third
    entry negative), raw slices appear rotated relative to Allen QC conventions.
    """
    if permvec is None or len(permvec) != 3:
        return False
    return permvec[2] < 0


@dataclass(frozen=True)
class SliceDisplayTransform:
    """Display remap applied to a 2D slice (``np.rot90`` then optional flips)."""

    rot90_k: int = 0
    flip_ud: bool = False
    flip_lr: bool = False


def allen_display_transform(cut_axis: int) -> SliceDisplayTransform:
    """Per-view display transform when atlas Z is flipped (``permvec[2] < 0``).

    Verified visually against the Perens/Gubra-oriented sample volume in both the
    registered (``volumereg``) and atlas-grid (match-points) spaces:

    - cut axis 1 (Y): 90° CCW
    - cut axis 2 (X): 90° CCW
    - cut axis 3 (Z): 90° CW then horizontal flip
    """
    if cut_axis == 1:
        return SliceDisplayTransform(rot90_k=1)
    if cut_axis == 2:
        return SliceDisplayTransform(rot90_k=1)
    if cut_axis == 3:
        return SliceDisplayTransform(rot90_k=3, flip_lr=True)
    msg = f"cut_axis must be 1, 2, or 3, got {cut_axis}"
    raise ValueError(msg)


def allen_display_rot90_k(cut_axis: int) -> int:
    """Legacy helper: rotation k only (ignores post-rotation flips)."""
    return allen_display_transform(cut_axis).rot90_k


def _display_shape(orig_shape: tuple[int, int], rot90_k: int) -> tuple[int, int]:
    height, width = orig_shape
    return (width, height) if (rot90_k % 4) % 2 == 1 else (height, width)


def rotate_slice_pixel_coords(
    row: np.ndarray | float,
    col: np.ndarray | float,
    shape: tuple[int, int],
    k: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Map (row, col) in a slice through the same ``np.rot90`` applied to the image."""
    height, width = shape
    row_a = np.asarray(row, dtype=float)
    col_a = np.asarray(col, dtype=float)
    k = k % 4
    if k == 0:
        return row_a, col_a
    if k == 1:
        return col_a, height - 1 - row_a
    if k == 2:
        return height - 1 - row_a, width - 1 - col_a
    # k == 3, equivalent to k == -1
    return width - 1 - col_a, row_a


def map_slice_pixels_to_display(
    row: np.ndarray | float,
    col: np.ndarray | float,
    slice_shape: tuple[int, int],
    transform: SliceDisplayTransform,
) -> tuple[np.ndarray, np.ndarray]:
    """Map raw slice (row, col) to displayed image pixel coordinates."""
    disp_row, disp_col = rotate_slice_pixel_coords(row, col, slice_shape, transform.rot90_k)
    disp_height, disp_width = _display_shape(slice_shape, transform.rot90_k)
    if transform.flip_ud:
        disp_row = disp_height - 1 - disp_row
    if transform.flip_lr:
        disp_col = disp_width - 1 - disp_col
    return disp_row, disp_col


def map_display_pixels_to_slice(
    disp_row: np.ndarray | float,
    disp_col: np.ndarray | float,
    slice_shape: tuple[int, int],
    transform: SliceDisplayTransform,
) -> tuple[np.ndarray, np.ndarray]:
    """Inverse of :func:`map_slice_pixels_to_display`."""
    height, width = slice_shape
    disp_height, disp_width = _display_shape(slice_shape, transform.rot90_k)
    row_a = np.asarray(disp_row, dtype=float)
    col_a = np.asarray(disp_col, dtype=float)
    if transform.flip_lr:
        col_a = disp_width - 1 - col_a
    if transform.flip_ud:
        row_a = disp_height - 1 - row_a

    k = transform.rot90_k % 4
    if k == 0:
        return row_a, col_a
    if k == 1:
        orig_col = row_a
        orig_row = height - 1 - col_a
    elif k == 2:
        orig_row = height - 1 - row_a
        orig_col = width - 1 - col_a
    else:  # k == 3
        orig_row = col_a
        orig_col = width - 1 - row_a
    return orig_row, orig_col


def apply_slice_display_transform(
    slice_2d: np.ndarray,
    transform: SliceDisplayTransform,
) -> np.ndarray:
    """Apply rotation and optional flips for QC / GUI display."""
    out = np.rot90(slice_2d, k=transform.rot90_k) if transform.rot90_k else slice_2d
    if transform.flip_ud:
        out = np.flipud(out)
    if transform.flip_lr:
        out = np.fliplr(out)
    return out


def prepare_display_slice(
    slice_2d: np.ndarray,
    cut_axis: int,
    permvec: list[int] | None,
) -> np.ndarray:
    """Rotate a 2D slice for Allen-standard viewing when orientation flips atlas Z."""
    if not needs_allen_display_reorientation(permvec):
        return slice_2d
    return apply_slice_display_transform(slice_2d, allen_display_transform(cut_axis))


def layer_xy_from_slice_pixels(
    row: np.ndarray | float,
    col: np.ndarray | float,
    slice_shape: tuple[int, int],
    cut_axis: int,
    permvec: list[int] | None,
) -> np.ndarray:
    """Map slice (row, col) pixel coords to napari layer (x, y) after display rotation."""
    if not needs_allen_display_reorientation(permvec):
        stacked = np.column_stack([np.asarray(col, dtype=float), np.asarray(row, dtype=float)])
        return stacked

    transform = allen_display_transform(cut_axis)
    disp_row, disp_col = map_slice_pixels_to_display(row, col, slice_shape, transform)
    return np.column_stack([disp_col, disp_row])


def slice_pixels_from_layer_xy(
    xy: np.ndarray,
    slice_shape: tuple[int, int],
    cut_axis: int,
    permvec: list[int] | None,
) -> tuple[np.ndarray, np.ndarray]:
    """Inverse of :func:`layer_xy_from_slice_pixels` for click/drag events."""
    disp_col = np.asarray(xy[:, 0], dtype=float)
    disp_row = np.asarray(xy[:, 1], dtype=float)
    if not needs_allen_display_reorientation(permvec):
        return disp_row, disp_col

    transform = allen_display_transform(cut_axis)
    return map_display_pixels_to_slice(disp_row, disp_col, slice_shape, transform)
