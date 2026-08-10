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
        return width - 1 - col_a, row_a
    if k == 2:
        return height - 1 - row_a, width - 1 - col_a
    return col_a, height - 1 - row_a


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
        orig_row = col_a
        orig_col = width - 1 - row_a
    elif k == 2:
        orig_row = height - 1 - row_a
        orig_col = width - 1 - col_a
    else:
        orig_row = height - 1 - col_a
        orig_col = row_a
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


def cut_axis_for_permuted_volume_axis(
    permuted_axis_0based: int,
    permute_sample_to_atlas: list[int],
) -> int:
    """Map a 0-based axis on the permuted registration grid to a native atlas cut axis."""
    if permuted_axis_0based not in {0, 1, 2}:
        msg = f"permuted_axis_0based must be 0, 1, or 2, got {permuted_axis_0based}"
        raise ValueError(msg)
    if len(permute_sample_to_atlas) != 3:
        msg = f"permute_sample_to_atlas must have length 3, got {permute_sample_to_atlas}"
        raise ValueError(msg)
    native_axis = abs(int(permute_sample_to_atlas[permuted_axis_0based])) - 1
    if native_axis not in {0, 1, 2}:
        msg = f"Invalid permute entry {permute_sample_to_atlas!r}"
        raise ValueError(msg)
    return native_axis + 1


def coronal_napari_permutation(atlas_provider: str) -> tuple[int, tuple[int, int, int]]:
    """Return 1-based coronal cut axis and a YXZ permutation that moves it last (Napari Z)."""
    coronal_cut = cut_axis_for_plot_dim(atlas_provider, 1)
    coronal_axis = coronal_cut - 1
    in_plane = tuple(axis for axis in range(3) if axis != coronal_axis)
    return coronal_cut, in_plane + (coronal_axis,)


def coronal_napari_permutation_for_registration(
    atlas_provider: str,
    permute_sample_to_atlas: list[int],
) -> tuple[int, tuple[int, int, int]]:
    """Return coronal cut axis and YXZ permutation for a permuted registration grid."""
    coronal_cut = cut_axis_for_plot_dim(atlas_provider, 1)
    coronal_axis: int | None = None
    for axis in range(3):
        if cut_axis_for_permuted_volume_axis(axis, permute_sample_to_atlas) == coronal_cut:
            coronal_axis = axis
            break
    if coronal_axis is None:
        msg = (
            f"Could not map coronal cut axis {coronal_cut} through "
            f"permute_sample_to_atlas={permute_sample_to_atlas!r}"
        )
        raise ValueError(msg)
    in_plane = tuple(axis for axis in range(3) if axis != coronal_axis)
    return coronal_cut, in_plane + (coronal_axis,)


def _volume_yxz_to_napari_coronal_zyx(
    volume_yxz: np.ndarray,
    *,
    atlas_provider: str,
    coronal_cut: int,
    coronal_perm: tuple[int, int, int],
) -> np.ndarray:
    """Shared Napari reorientation: move coronal to Z and apply canonical in-plane QC."""
    vol = np.asarray(volume_yxz)
    if vol.ndim != 3:
        msg = f"Expected 3D YXZ volume, got shape {vol.shape}"
        raise ValueError(msg)

    oriented = np.transpose(vol, coronal_perm)
    napari_vol = np.transpose(oriented, (2, 0, 1))

    transform = canonical_view_transform(atlas_provider, coronal_cut)
    if (
        transform.rot90_k == 0
        and not transform.flip_ud
        and not transform.flip_lr
    ):
        return napari_vol

    sample_slice = apply_slice_display_transform(napari_vol[0], transform)
    out = np.empty(
        (napari_vol.shape[0],) + sample_slice.shape,
        dtype=napari_vol.dtype,
    )
    out[0] = sample_slice
    for z in range(1, napari_vol.shape[0]):
        out[z] = apply_slice_display_transform(napari_vol[z], transform)
    return out


def atlas_volume_yxz_to_napari_zyx(
    volume_yxz: np.ndarray,
    *,
    atlas_provider: str,
) -> np.ndarray:
    """Reorient a native atlas-order YXZ volume for Napari (coronal = scroll axis Z).

    Applies the same per-slice canonical QC transforms used in registration plots so
    Perens / Gubra volumes are not upside-down relative to Allen when inspected.
    """
    coronal_cut, perm = coronal_napari_permutation(atlas_provider)
    return _volume_yxz_to_napari_coronal_zyx(
        volume_yxz,
        atlas_provider=atlas_provider,
        coronal_cut=coronal_cut,
        coronal_perm=perm,
    )


def registration_volume_yxz_to_napari_zyx(
    volume_yxz: np.ndarray,
    *,
    atlas_provider: str,
    permute_sample_to_atlas: list[int],
) -> np.ndarray:
    """Reorient a permuted registration-grid volume for Napari coronal scrolling."""
    coronal_cut, perm = coronal_napari_permutation_for_registration(
        atlas_provider,
        permute_sample_to_atlas,
    )
    return _volume_yxz_to_napari_coronal_zyx(
        volume_yxz,
        atlas_provider=atlas_provider,
        coronal_cut=coronal_cut,
        coronal_perm=perm,
    )


def _points_xyz_to_napari_coronal_zyx(
    coords_xyz_1based: np.ndarray,
    *,
    atlas_provider: str,
    volume_shape_yxz: tuple[int, int, int],
    coronal_cut: int,
    coronal_perm: tuple[int, int, int],
) -> np.ndarray:
    """Map 1-based ``(x, y, z)`` points to Napari ``(z, y, x)`` with coronal QC layout."""
    pts = np.asarray(coords_xyz_1based, dtype=float)
    if pts.size == 0:
        return np.zeros((0, 3), dtype=np.float64)
    if pts.ndim != 2 or pts.shape[1] < 3:
        msg = f"Expected Nx3+ points, got shape {pts.shape}"
        raise ValueError(msg)

    native_yxz = np.column_stack(
        [pts[:, 1] - 1.0, pts[:, 0] - 1.0, pts[:, 2] - 1.0]
    )
    permuted = native_yxz[:, coronal_perm]
    permuted_shape = tuple(int(volume_shape_yxz[axis]) for axis in coronal_perm)
    slice_shape = (permuted_shape[0], permuted_shape[1])

    transform = canonical_view_transform(atlas_provider, coronal_cut)
    disp_row, disp_col = map_slice_pixels_to_display(
        permuted[:, 0],
        permuted[:, 1],
        slice_shape,
        transform,
    )
    return np.column_stack([permuted[:, 2], disp_row, disp_col])


def atlas_points_xyz_to_napari_zyx(
    coords_xyz_1based: np.ndarray,
    *,
    atlas_provider: str,
    volume_shape_yxz: tuple[int, int, int],
) -> np.ndarray:
    """Map 1-based atlas ``(x, y, z)`` points to Napari ``(z, y, x)`` with coronal QC layout."""
    coronal_cut, perm = coronal_napari_permutation(atlas_provider)
    return _points_xyz_to_napari_coronal_zyx(
        coords_xyz_1based,
        atlas_provider=atlas_provider,
        volume_shape_yxz=volume_shape_yxz,
        coronal_cut=coronal_cut,
        coronal_perm=perm,
    )


def registration_points_xyz_to_napari_zyx(
    coords_xyz_1based: np.ndarray,
    *,
    atlas_provider: str,
    volume_shape_yxz: tuple[int, int, int],
    permute_sample_to_atlas: list[int],
) -> np.ndarray:
    """Map 1-based registration-grid points to Napari ZYX with coronal QC layout."""
    coronal_cut, perm = coronal_napari_permutation_for_registration(
        atlas_provider,
        permute_sample_to_atlas,
    )
    return _points_xyz_to_napari_coronal_zyx(
        coords_xyz_1based,
        atlas_provider=atlas_provider,
        volume_shape_yxz=volume_shape_yxz,
        coronal_cut=coronal_cut,
        coronal_perm=perm,
    )
