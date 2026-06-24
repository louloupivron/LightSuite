"""Tests for canonical atlas display profiles."""

from __future__ import annotations

import numpy as np
import pytest

from lightsuite.atlas.display import (
    PLOT_VIEW_NAMES,
    canonical_view_slice,
    canonical_view_transform,
    cut_axis_for_plot_dim,
    display_provider_for_atlas,
    get_display_profile,
    map_display_pixels_to_slice,
    map_slice_pixels_to_display,
)
from lightsuite.gui.slices import layer_xy_from_slice_pixels, slice_pixels_from_layer_xy, volume_index_to_image


@pytest.mark.parametrize("provider", ["allen", "perens", "perens_brainglobe", "princeton_brainglobe"])
def test_plot_dim_maps_to_distinct_cut_axes(provider: str) -> None:
    profile = get_display_profile(provider)
    mapped = {profile.cut_axis_for_plot_dim(d) for d in (1, 2, 3)}
    assert mapped == {1, 2, 3}
    assert PLOT_VIEW_NAMES == {1: "coronal", 2: "sagittal", 3: "horizontal"}


def test_allen_plot_dim_to_cut_axis_mapping() -> None:
    assert cut_axis_for_plot_dim("allen", 1) == 1
    assert cut_axis_for_plot_dim("allen", 2) == 3
    assert cut_axis_for_plot_dim("allen", 3) == 2


def test_perens_plot_dim_to_cut_axis_mapping() -> None:
    assert cut_axis_for_plot_dim("perens", 1) == 2
    assert cut_axis_for_plot_dim("perens", 2) == 1
    assert cut_axis_for_plot_dim("perens", 3) == 3


def test_perens_brainglobe_plot_dim_to_cut_axis_mapping() -> None:
    assert cut_axis_for_plot_dim("perens_brainglobe", 1) == 1
    assert cut_axis_for_plot_dim("perens_brainglobe", 2) == 3
    assert cut_axis_for_plot_dim("perens_brainglobe", 3) == 2


def test_display_provider_for_atlas_brainglobe_perens() -> None:
    assert display_provider_for_atlas("perens", atlas_source="files") == "perens"
    assert display_provider_for_atlas("perens", atlas_source="brainglobe") == "perens_brainglobe"
    assert display_provider_for_atlas("princeton", atlas_source="brainglobe") == "princeton_brainglobe"
    assert display_provider_for_atlas("allen", atlas_source="brainglobe") == "allen"


def test_princeton_brainglobe_plot_dim_to_cut_axis_mapping() -> None:
    assert cut_axis_for_plot_dim("princeton_brainglobe", 1) == 1
    assert cut_axis_for_plot_dim("princeton_brainglobe", 2) == 3
    assert cut_axis_for_plot_dim("princeton_brainglobe", 3) == 2


def test_perens_brainglobe_coronal_cut_is_identity() -> None:
    transform = canonical_view_transform("perens_brainglobe", 1)
    assert transform.rot90_k == 0
    assert not transform.flip_ud
    assert not transform.flip_lr


def test_allen_coronal_cut_is_identity() -> None:
    """Allen native coronal slices are already upright; no in-plane remap needed."""
    transform = canonical_view_transform("allen", 1)
    assert transform.rot90_k == 0
    assert not transform.flip_ud
    assert not transform.flip_lr


def test_perens_brainglobe_matches_file_perens_views() -> None:
    """BrainGlobe Perens QC panels align with local Gubra NIfTI conventions."""
    from pathlib import Path

    file_path = Path(
        "atlas/perens_20um/LSFM-mouse-brain-atlas-master/LSFM_atlas_files/gubra_template_olf.nii.gz"
    )
    bg_path = Path.home() / ".brainglobe/perens_lsfm_mouse_20um_v1.2/reference.tiff"
    if not file_path.is_file() or not bg_path.is_file():
        pytest.skip("Perens file + BrainGlobe atlases not available locally")

    from lightsuite.atlas.io import load_atlas_volume

    file_tpl = load_atlas_volume(file_path).astype(np.float32)
    bg_tpl = load_atlas_volume(bg_path).astype(np.float32)

    def ncc(a: np.ndarray, b: np.ndarray) -> float:
        a = a - a.mean()
        b = b - b.mean()
        return float((a * b).sum() / (np.linalg.norm(a) * np.linalg.norm(b)))

    for plot_dim in (1, 2, 3):
        ca_file = cut_axis_for_plot_dim("perens", plot_dim)
        ca_bg = cut_axis_for_plot_dim("perens_brainglobe", plot_dim)
        mid_f = file_tpl.shape[ca_file - 1] // 2
        mid_b = bg_tpl.shape[ca_bg - 1] // 2
        ref = canonical_view_slice(
            volume_index_to_image(file_tpl, np.array([mid_f, ca_file])),
            atlas_provider="perens",
            cut_axis=ca_file,
        )
        cand = canonical_view_slice(
            volume_index_to_image(bg_tpl, np.array([mid_b, ca_bg])),
            atlas_provider="perens_brainglobe",
            cut_axis=ca_bg,
        )
        assert ref.shape == cand.shape
        assert ncc(ref, cand) > 0.97


def test_perens_sagittal_cut_uses_rotation() -> None:
    transform = canonical_view_transform("perens", 1)
    assert transform.rot90_k == 1
    assert transform.flip_lr


def test_canonical_view_independent_of_permvec() -> None:
    sl = np.arange(12, dtype=np.uint8).reshape(3, 4)
    a = canonical_view_slice(sl, atlas_provider="allen", cut_axis=1)
    b = canonical_view_slice(sl, atlas_provider="allen", cut_axis=1)
    assert np.array_equal(a, b)


@pytest.mark.parametrize("provider", ["allen", "perens", "perens_brainglobe", "princeton_brainglobe"])
@pytest.mark.parametrize("cut_axis", [1, 2, 3])
def test_display_coordinate_roundtrip(provider: str, cut_axis: int) -> None:
    slice_shape = (24, 32)
    transform = canonical_view_transform(provider, cut_axis)
    row, col = 7.0, 11.0
    disp_row, disp_col = map_slice_pixels_to_display(row, col, slice_shape, transform)
    back_row, back_col = map_display_pixels_to_slice(disp_row, disp_col, slice_shape, transform)
    assert np.allclose(back_row, row)
    assert np.allclose(back_col, col)

    layer_xy = layer_xy_from_slice_pixels(row, col, slice_shape, cut_axis, provider)
    back_row2, back_col2 = slice_pixels_from_layer_xy(layer_xy, slice_shape, cut_axis, provider)
    assert np.allclose(back_row2, row)
    assert np.allclose(back_col2, col)
