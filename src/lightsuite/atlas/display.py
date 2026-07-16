"""Canonical anatomical display profiles for atlas QC plots and GUIs.

Display geometry is independent of sample ``permvec``: volumes shown in plots and
Napari tools are already in atlas axis order after registration permutation.

Plot panels use a fixed convention regardless of atlas native axis order:

- ``plot_dim`` 1 → coronal
- ``plot_dim`` 2 → sagittal
- ``plot_dim`` 3 → horizontal

Each atlas maps those panels to the volume ``cut_axis`` that actually contains
that anatomy, then applies a fixed in-plane transform for nose-up / nose-left QC.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

PLOT_VIEW_NAMES: dict[int, str] = {1: "coronal", 2: "sagittal", 3: "horizontal"}


@dataclass(frozen=True)
class SliceDisplayTransform:
    """Display remap applied to a 2D slice (``np.rot90`` then optional flips)."""

    rot90_k: int = 0
    flip_ud: bool = False
    flip_lr: bool = False


def _display_shape(orig_shape: tuple[int, int], rot90_k: int) -> tuple[int, int]:
    height, width = orig_shape
    return (width, height) if (rot90_k % 4) % 2 == 1 else (height, width)


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
    else:
        orig_row = col_a
        orig_col = width - 1 - row_a
    return orig_row, orig_col


# Transforms keyed by native volume cut_axis (calibrated on atlas templates).
# Allen 10 µm native cuts: 1=coronal, 2=horizontal, 3=sagittal.
_ALLEN_CUT_TRANSFORMS: dict[int, SliceDisplayTransform] = {
    1: SliceDisplayTransform(rot90_k=0, flip_ud=False, flip_lr=False),
    2: SliceDisplayTransform(rot90_k=0, flip_ud=True, flip_lr=False),
    3: SliceDisplayTransform(rot90_k=1, flip_ud=True, flip_lr=True),
}

# Perens / Gubra 20 µm native cuts: 1=sagittal, 2=coronal, 3=horizontal.
_PERENS_CUT_TRANSFORMS: dict[int, SliceDisplayTransform] = {
    1: SliceDisplayTransform(rot90_k=1, flip_ud=False, flip_lr=True),
    2: SliceDisplayTransform(rot90_k=1, flip_ud=False, flip_lr=True),
    3: SliceDisplayTransform(rot90_k=1, flip_ud=False, flip_lr=False),
}

# BrainGlobe ``perens_lsfm_mouse_20um`` (ASR, shape 621×323×461): 1=coronal, 2=horizontal,
# 3=sagittal. Calibrated against local Gubra NIfTIs so dim{1,2,3} QC panels match.
_PERENS_BRAINGLOBE_CUT_TRANSFORMS: dict[int, SliceDisplayTransform] = {
    1: SliceDisplayTransform(rot90_k=0, flip_ud=False, flip_lr=False),
    2: SliceDisplayTransform(rot90_k=0, flip_ud=True, flip_lr=True),
    3: SliceDisplayTransform(rot90_k=1, flip_ud=True, flip_lr=True),
}

# BrainGlobe ``princeton_mouse_20um`` (ASR, shape 640×352×540): same panel→axis map as Allen
# with cut-specific in-plane transforms calibrated against downsampled Allen 10 µm.
_PRINCETON_BRAINGLOBE_CUT_TRANSFORMS: dict[int, SliceDisplayTransform] = {
    1: SliceDisplayTransform(rot90_k=0, flip_ud=False, flip_lr=False),
    2: SliceDisplayTransform(rot90_k=0, flip_ud=False, flip_lr=False),
    3: SliceDisplayTransform(rot90_k=1, flip_ud=True, flip_lr=True),
}

# Plot panel dim (1=coronal, 2=sagittal, 3=horizontal) → volume cut_axis.
_ALLEN_PLOT_DIM_TO_CUT_AXIS: dict[int, int] = {1: 1, 2: 3, 3: 2}
_PERENS_PLOT_DIM_TO_CUT_AXIS: dict[int, int] = {1: 2, 2: 1, 3: 3}
_PERENS_BRAINGLOBE_PLOT_DIM_TO_CUT_AXIS: dict[int, int] = {1: 1, 2: 3, 3: 2}
_PRINCETON_BRAINGLOBE_PLOT_DIM_TO_CUT_AXIS: dict[int, int] = {1: 1, 2: 3, 3: 2}

# Spinal cord (Fiederling): MATLAB matchControlPointsSpine uses raw volumeIdtoImage slices
# with imagesc YDir reverse only — no rot90 / flips in the canonical QC map.
_CORD_CUT_TRANSFORMS: dict[int, SliceDisplayTransform] = {
    1: SliceDisplayTransform(rot90_k=0, flip_ud=False, flip_lr=False),
    2: SliceDisplayTransform(rot90_k=0, flip_ud=False, flip_lr=False),
    3: SliceDisplayTransform(rot90_k=0, flip_ud=False, flip_lr=False),
}
_CORD_PLOT_DIM_TO_CUT_AXIS: dict[int, int] = {1: 1, 2: 2, 3: 3}


@dataclass(frozen=True)
class AtlasDisplayProfile:
    """Fixed QC viewing conventions for one brain atlas provider."""

    provider: str
    axis_names: tuple[str, str, str]
    plot_dim_to_cut_axis: dict[int, int]
    cut_transforms: dict[int, SliceDisplayTransform]

    def cut_axis_for_plot_dim(self, plot_dim: int) -> int:
        if plot_dim not in self.plot_dim_to_cut_axis:
            msg = f"plot_dim must be 1, 2, or 3, got {plot_dim}"
            raise ValueError(msg)
        return self.plot_dim_to_cut_axis[plot_dim]

    def transform_for_cut_axis(self, cut_axis: int) -> SliceDisplayTransform:
        if cut_axis not in self.cut_transforms:
            msg = f"cut_axis must be 1, 2, or 3, got {cut_axis}"
            raise ValueError(msg)
        return self.cut_transforms[cut_axis]


_DISPLAY_PROFILES: dict[str, AtlasDisplayProfile] = {
    "allen": AtlasDisplayProfile(
        provider="allen",
        axis_names=("AP", "DV", "LR"),
        plot_dim_to_cut_axis=_ALLEN_PLOT_DIM_TO_CUT_AXIS,
        cut_transforms=_ALLEN_CUT_TRANSFORMS,
    ),
    "perens": AtlasDisplayProfile(
        provider="perens",
        axis_names=("L", "P", "S"),
        plot_dim_to_cut_axis=_PERENS_PLOT_DIM_TO_CUT_AXIS,
        cut_transforms=_PERENS_CUT_TRANSFORMS,
    ),
    "perens_brainglobe": AtlasDisplayProfile(
        provider="perens_brainglobe",
        axis_names=("P", "S", "L"),
        plot_dim_to_cut_axis=_PERENS_BRAINGLOBE_PLOT_DIM_TO_CUT_AXIS,
        cut_transforms=_PERENS_BRAINGLOBE_CUT_TRANSFORMS,
    ),
    "princeton_brainglobe": AtlasDisplayProfile(
        provider="princeton_brainglobe",
        axis_names=("AP", "DV", "LR"),
        plot_dim_to_cut_axis=_PRINCETON_BRAINGLOBE_PLOT_DIM_TO_CUT_AXIS,
        cut_transforms=_PRINCETON_BRAINGLOBE_CUT_TRANSFORMS,
    ),
    # BrainGlobe ``whs_sd_rat_39um`` (ASR): same panel→axis map as Princeton/Allen BG packs.
    "rat_brainglobe": AtlasDisplayProfile(
        provider="rat_brainglobe",
        axis_names=("AP", "DV", "LR"),
        plot_dim_to_cut_axis=_PRINCETON_BRAINGLOBE_PLOT_DIM_TO_CUT_AXIS,
        cut_transforms=_PRINCETON_BRAINGLOBE_CUT_TRANSFORMS,
    ),
    "cord": AtlasDisplayProfile(
        provider="cord",
        axis_names=("Y", "X", "Z"),
        plot_dim_to_cut_axis=_CORD_PLOT_DIM_TO_CUT_AXIS,
        cut_transforms=_CORD_CUT_TRANSFORMS,
    ),
}


def get_display_profile(atlas_provider: str) -> AtlasDisplayProfile:
    """Return the display profile for a brain atlas provider id."""
    key = atlas_provider.lower().strip()
    if key not in _DISPLAY_PROFILES:
        msg = (
            f"Unknown atlas provider {atlas_provider!r} for display. "
            f"Expected: {', '.join(sorted(_DISPLAY_PROFILES))}"
        )
        raise ValueError(msg)
    return _DISPLAY_PROFILES[key]


def display_provider_for_atlas(
    brain_atlas: str,
    *,
    atlas_source: str = "files",
) -> str:
    """Resolve the QC display profile id for an atlas provider and source backend."""
    atlas_id = brain_atlas.lower().strip()
    source = atlas_source.lower().strip()
    if source == "brainglobe":
        bg_key = f"{atlas_id}_brainglobe"
        if bg_key in _DISPLAY_PROFILES:
            return bg_key
    return atlas_id


def cut_axis_for_plot_dim(atlas_provider: str, plot_dim: int) -> int:
    """Map a plot panel (1=coronal, 2=sagittal, 3=horizontal) to a volume cut axis."""
    return get_display_profile(atlas_provider).cut_axis_for_plot_dim(plot_dim)


def canonical_view_transform(atlas_provider: str, cut_axis: int) -> SliceDisplayTransform:
    """Return the fixed 2D remap for one native volume cut axis."""
    return get_display_profile(atlas_provider).transform_for_cut_axis(cut_axis)


def canonical_view_slice(
    slice_2d: np.ndarray,
    *,
    atlas_provider: str,
    cut_axis: int,
) -> np.ndarray:
    """Map a native atlas-order 2D slice to the canonical QC orientation."""
    transform = canonical_view_transform(atlas_provider, cut_axis)
    return apply_slice_display_transform(slice_2d, transform)


def canonical_view_name(plot_dim: int) -> str:
    """Human-readable anatomical view label for a plot panel."""
    return PLOT_VIEW_NAMES.get(plot_dim, f"dim{plot_dim}")
