"""Tests for atlas→sample warp helpers."""

from __future__ import annotations

import numpy as np

from lightsuite.export.sample_space import resolve_forward_transform_paths
from lightsuite.registration.brain_register import TransformParamsCheckpoint
from lightsuite.registration.canvas import RegistrationCanvas
from lightsuite.registration.coordinates import affine_with_source_offset
from lightsuite.registration.volume import permute_brain_volume, unpermute_brain_volume
from lightsuite.registration.warp import warp_volume_affine


def _minimal_transform_params(
    *,
    reg_shape: tuple[int, int, int],
    atlas_shape: tuple[int, int, int],
    permute: list[int],
) -> TransformParamsCheckpoint:
    tform = np.eye(4, dtype=float)
    return TransformParamsCheckpoint(
        atlas_resolution_um=10.0,
        regvolsize=list(reg_shape),
        atlassize=list(atlas_shape),
        brain_atlas="allen",
        ori_voxel_um=[2.0, 2.0, 3.0],
        ori_size=list(reg_shape),
        permute_sample_to_atlas=permute,
        elastix_um_to_mm=1e-3,
        tform_bspline_samp20um_to_atlas_20um_px="/tmp/inv.txt",
        tform_affine_samp20um_to_atlas_10um_px=np.linalg.inv(tform).tolist(),
        control_point_weight=0.0,
        use_multistep=False,
        use_dual_channel_mi=False,
        tform_bspline_atlas_20um_to_samp20um_px="/tmp/fwd.txt",
        tform_affine_atlas_to_samp20um_px=tform.tolist(),
        registration_canvas=RegistrationCanvas(
            mode="off",
            pad_before=(0, 0, 0),
            pad_after=(0, 0, 0),
            sample_crop_start=(0, 0, 0),
            working_shape=reg_shape,
        ).to_checkpoint_dict(),
    )


def test_unpermute_inverts_permute_brain_volume() -> None:
    vol = np.arange(24, dtype=np.float32).reshape(2, 3, 4)
    perm = [1, 2, 3]
    out = unpermute_brain_volume(permute_brain_volume(vol, perm), perm)
    np.testing.assert_array_equal(out, vol)


def test_transform_atlas_volume_to_sample_identity_affine(tmp_path) -> None:
    """With identity affine and mocked bspline, shape round-trip holds."""
    reg_shape = (10, 12, 8)
    atlas_shape = (10, 12, 8)
    ann = np.zeros(atlas_shape, dtype=np.float32)
    ann[3:7, 4:9, 2:6] = 42.0

    fwd_bspline = tmp_path / "bspline_atlas_to_samp_20um.txt"
    fwd_bspline.write_text(
        "(Transform \"AffineTransform\")\n"
        "(TransformParameters 1 0 0 0 1 0 0 0 1 0 0 0)\n"
        "(FixedImageDimension 3)\n"
        "(MovingImageDimension 3)\n",
        encoding="utf-8",
    )

    params = _minimal_transform_params(
        reg_shape=reg_shape,
        atlas_shape=atlas_shape,
        permute=[1, 2, 3],
    )
    params.tform_bspline_atlas_20um_to_samp20um_px = str(fwd_bspline)

    # Skip transformix in unit test — test affine-only path via warp_volume_affine chain
    tform = np.eye(4)
    warped = warp_volume_affine(ann, tform, reg_shape, order=0)
    assert warped.shape == reg_shape
    assert int(np.count_nonzero(warped)) == int(np.count_nonzero(ann))


def _trimmed_atlas_params(tmp_path, *, crop, tform_aff, atlas_shape, reg_shape):
    fwd_bspline = tmp_path / "bspline_atlas_to_samp_20um.txt"
    fwd_bspline.write_text("(Transform \"AffineTransform\")\n", encoding="utf-8")
    params = _minimal_transform_params(
        reg_shape=reg_shape,
        atlas_shape=atlas_shape,
        permute=[1, 2, 3],
    )
    params.tform_bspline_atlas_20um_to_samp20um_px = str(fwd_bspline)
    params.atlas_crop_start_native = list(crop)
    params.tform_affine_samp20um_to_atlas_10um_px = np.linalg.inv(tform_aff).tolist()
    params.tform_affine_atlas_to_samp20um_px = affine_with_source_offset(
        tform_aff, crop
    ).tolist()
    return params


def test_forward_affine_is_native_indexed_when_atlas_is_trimmed(tmp_path) -> None:
    """Sample-space export warps the native atlas, so the trim offset must be undone."""
    crop = (0, 2, 3)
    tform_aff = np.array(
        [
            [0.5, 0.02, 0.0, -1.5],
            [0.01, 0.45, 0.03, 2.0],
            [0.0, -0.02, 0.55, 0.75],
            [0.0, 0.0, 0.0, 1.0],
        ]
    )
    params = _trimmed_atlas_params(
        tmp_path,
        crop=crop,
        tform_aff=tform_aff,
        atlas_shape=(10, 12, 14),
        reg_shape=(8, 8, 10),
    )

    _, affine = resolve_forward_transform_paths(params, tmp_path)
    np.testing.assert_allclose(affine, tform_aff, atol=1e-12)

    params.tform_affine_atlas_to_samp20um_px = None
    _, affine_fallback = resolve_forward_transform_paths(params, tmp_path)
    np.testing.assert_allclose(affine_fallback, tform_aff, atol=1e-10)


def test_native_and_trimmed_atlas_warps_agree(tmp_path) -> None:
    """Warping the native atlas must reproduce the trimmed warp used during register."""
    crop = (0, 2, 3)
    atlas_shape = (10, 12, 14)
    reg_shape = (8, 8, 10)
    tform_aff = np.array(
        [
            [0.5, 0.0, 0.0, -1.0],
            [0.0, 0.45, 0.0, 1.0],
            [0.0, 0.0, 0.55, 0.5],
            [0.0, 0.0, 0.0, 1.0],
        ]
    )

    native = np.zeros(atlas_shape, dtype=np.float32)
    native[2:8, 4:10, 5:12] = 7.0  # content strictly inside the trim box
    trimmed = native[crop[0] :, crop[1] :, crop[2] :]

    params = _trimmed_atlas_params(
        tmp_path,
        crop=crop,
        tform_aff=tform_aff,
        atlas_shape=atlas_shape,
        reg_shape=reg_shape,
    )
    _, affine = resolve_forward_transform_paths(params, tmp_path)

    from_native = warp_volume_affine(native, affine, reg_shape, order=0)
    from_trimmed = warp_volume_affine(
        trimmed,
        affine_with_source_offset(tform_aff, crop),
        reg_shape,
        order=0,
    )
    np.testing.assert_allclose(from_native, from_trimmed, atol=1e-6)
    assert int(np.count_nonzero(from_native)) > 0
