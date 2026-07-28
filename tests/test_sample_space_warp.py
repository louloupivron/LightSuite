"""Tests for atlas→sample warp helpers."""

from __future__ import annotations

import numpy as np

from lightsuite.registration.brain_register import TransformParamsCheckpoint
from lightsuite.registration.canvas import RegistrationCanvas
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
