"""Tests for expanded registration warp canvas padding."""

from __future__ import annotations

import numpy as np

from lightsuite.config.models import RegistrationCanvasMode
from lightsuite.registration.canvas import (
    RegistrationCanvas,
    WarpCanvasPadding,
    apply_canvas_sample_crop,
    compute_registration_canvas,
    compute_vd_warp_canvas_padding,
    crop_from_warp_canvas,
    offset_volume_indices,
    pad_volume_for_warp_canvas,
    undo_canvas_sample_crop,
)
from lightsuite.registration.volume import permute_brain_volume, unpermute_brain_volume
from lightsuite.registration.warp import warp_volume_affine


def test_compute_vd_warp_canvas_padding_matches_template_aspect() -> None:
    sample_shape = (669, 794, 346)
    atlas_shape = (461, 621, 323)
    padding = compute_vd_warp_canvas_padding(sample_shape, atlas_shape)
    assert padding.pad_before[0] == 0
    assert padding.pad_after[0] == 0
    assert padding.pad_before[1] == 0
    assert padding.pad_after[1] == 0
    assert padding.pad_before[2] + padding.pad_after[2] == 67
    assert padding.padded_shape(sample_shape) == (669, 794, 413)


def test_compute_vd_warp_canvas_padding_zero_when_sample_is_taller() -> None:
    sample_shape = (100, 200, 150)
    atlas_shape = (461, 621, 323)
    padding = compute_vd_warp_canvas_padding(sample_shape, atlas_shape)
    assert padding.is_zero


def test_pad_and_crop_round_trip() -> None:
    volume = np.arange(24, dtype=np.float32).reshape(2, 3, 4)
    padding = compute_vd_warp_canvas_padding((2, 3, 4), (4, 6, 5))
    padded = pad_volume_for_warp_canvas(volume, padding, constant=-1.0)
    restored = crop_from_warp_canvas(padded, padding, volume.shape)
    np.testing.assert_array_equal(restored, volume)


def test_offset_volume_indices_shifts_sample_points() -> None:
    points = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
    shifted = offset_volume_indices(points, (0, 0, 10))
    np.testing.assert_array_equal(shifted, np.array([[1.0, 2.0, 13.0], [4.0, 5.0, 16.0]]))


def test_warp_canvas_padding_recovers_near_edge_atlas_labels() -> None:
    sample_shape = (40, 60, 30)
    ann = np.zeros((20, 30, 20), dtype=np.float32)
    ann[5:15, 8:22, 8:18] = 100.0

    tform = np.eye(4, dtype=float)
    tform[2, 3] = -15.0

    padding = WarpCanvasPadding((0, 0, 10), (0, 0, 10))
    warped_tight = warp_volume_affine(ann, tform, sample_shape, order=0)
    warped_padded = warp_volume_affine(
        ann,
        tform,
        padding.padded_shape(sample_shape),
        order=0,
        output_origin=padding.pad_before,
    )

    n_tight = int(np.count_nonzero(warped_tight > 1))
    n_padded = int(np.count_nonzero(warped_padded > 1))
    assert n_padded > n_tight


def test_imwarp_output_origin_embeds_original_grid() -> None:
    volume = np.zeros((10, 10, 10), dtype=np.float32)
    volume[5, 5, 5] = 1.0
    padding = WarpCanvasPadding((0, 0, 3), (0, 0, 4))
    warped = warp_volume_affine(
        volume,
        np.eye(4),
        padding.padded_shape(volume.shape),
        order=0,
        output_origin=padding.pad_before,
    )
    z0 = padding.pad_before[2]
    assert float(warped[5, 5, z0 + 5]) == 1.0
    assert float(warped[5, 5, 5]) == 0.0


def test_registration_canvas_pad_uses_registration_resolution_shapes() -> None:
    sample_shape = (794, 346, 669)
    atlas_native_shape = (1320, 792, 1050)
    atlas_reg_shape = (660, 396, 525)

    canvas_native = compute_registration_canvas(
        sample_shape,
        atlas_native_shape,
        RegistrationCanvasMode.PAD,
    )
    canvas_reg = compute_registration_canvas(
        sample_shape,
        atlas_reg_shape,
        RegistrationCanvasMode.PAD,
    )

    assert canvas_native.working_shape == atlas_native_shape
    assert canvas_reg.working_shape == (794, 396, 669)
    assert canvas_reg.pad_before == (0, 25, 0)
    assert canvas_reg.pad_after == (0, 25, 0)


def test_registration_canvas_pad_aligns_sample_and_warped_atlas() -> None:
    sample_shape = (40, 60, 30)
    atlas_reg_shape = (40, 70, 30)
    volume = np.ones(sample_shape, dtype=np.float32)
    ann = np.zeros((20, 30, 20), dtype=np.float32)
    ann[5:15, 8:22, 8:18] = 100.0

    reg_canvas = compute_registration_canvas(
        sample_shape,
        atlas_reg_shape,
        RegistrationCanvasMode.PAD,
    )
    volume_work = apply_canvas_sample_crop(volume, reg_canvas)
    warp_pad = WarpCanvasPadding(reg_canvas.pad_before, reg_canvas.pad_after)
    working_shape = reg_canvas.working_shape
    assert volume_work.shape == working_shape

    tform = np.eye(4, dtype=float)
    tform[2, 3] = -5.0
    avaffine = warp_volume_affine(
        ann,
        tform,
        working_shape,
        order=0,
        output_origin=warp_pad.pad_before,
    )

    sample_mask = volume_work > 0
    atlas_mask = avaffine > 1
    assert np.any(sample_mask & atlas_mask)


def test_undo_canvas_sample_crop_inverts_pad() -> None:
    volume = np.arange(24, dtype=np.float32).reshape(2, 3, 4)
    canvas = RegistrationCanvas(
        mode="pad",
        pad_before=(1, 0, 2),
        pad_after=(0, 1, 1),
        sample_crop_start=(0, 0, 0),
        working_shape=(3, 4, 7),
    )
    cropped = apply_canvas_sample_crop(volume, canvas)
    restored = undo_canvas_sample_crop(cropped, canvas, volume.shape)
    np.testing.assert_array_equal(restored, volume)


def test_undo_canvas_sample_crop_takes_working_grid_volume() -> None:
    sample_shape = (794, 346, 669)
    canvas = compute_registration_canvas(
        sample_shape,
        (660, 396, 525),
        RegistrationCanvasMode.PAD,
    )
    warped = np.zeros(canvas.working_shape, dtype=np.float32)

    restored = undo_canvas_sample_crop(warped, canvas, sample_shape)

    assert restored.shape == sample_shape


def test_unpermute_brain_volume_inverts_permute() -> None:
    volume = np.arange(24, dtype=np.int32).reshape(2, 3, 4)
    permvec = [-2, 3, -1]
    permuted = permute_brain_volume(volume, permvec)
    restored = unpermute_brain_volume(permuted, permvec)
    np.testing.assert_array_equal(restored, volume)
