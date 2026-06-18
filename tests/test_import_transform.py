"""Tests for sample-to-atlas point coordinate transforms."""

from __future__ import annotations

import numpy as np
import pytest

from lightsuite.import_.transform import (
    native_xyz_to_unpermuted_registration_yxz,
    permute_registration_indices_yxz,
    sample_points_to_registration_voxels,
)
from lightsuite.registration.brain_register import TransformParamsCheckpoint
from lightsuite.registration.volume import permute_brain_volume


def _transform_params() -> TransformParamsCheckpoint:
    return TransformParamsCheckpoint(
        atlas_resolution_um=10.0,
        regvolsize=[100, 80, 60],
        atlassize=[200, 180, 160],
        brain_atlas="allen",
        ori_voxel_um=[5.0, 5.0, 5.0],
        ori_size=[200, 160, 120],
        permute_sample_to_atlas=[1, 2, 3],
        elastix_um_to_mm=1e-3,
        tform_bspline_samp20um_to_atlas_20um_px="/tmp/bspline.txt",
        tform_affine_samp20um_to_atlas_10um_px=np.eye(4).tolist(),
        control_point_weight=0.2,
        use_multistep=True,
        use_dual_channel_mi=False,
    )


def test_native_xyz_to_unpermuted_registration_yxz() -> None:
    pts = np.array([[1.0, 1.0, 1.0], [81.0, 101.0, 61.0]])
    reg = native_xyz_to_unpermuted_registration_yxz(pts, ori_voxel_um=[5.0, 5.0, 5.0], registres_um=20.0)
    assert reg.shape == (2, 3)
    assert np.allclose(reg[0], [0.0, 0.0, 0.0])
    assert np.allclose(reg[1], [25.0, 20.0, 15.0])


def test_sample_points_to_registration_identity_perm() -> None:
    pts = np.array([[1.0, 1.0, 1.0], [81.0, 101.0, 61.0]])
    reg = sample_points_to_registration_voxels(pts, _transform_params(), registres_um=20.0)
    assert reg.shape == (2, 3)
    assert np.allclose(reg[0], [0.0, 0.0, 0.0])
    assert np.allclose(reg[1], [25.0, 20.0, 15.0])


def test_sample_points_to_registration_with_flip() -> None:
    params = _transform_params()
    params.permute_sample_to_atlas = [1, -2, 3]
    pts = np.array([[1.0, 1.0, 1.0]])
    reg = sample_points_to_registration_voxels(pts, params, registres_um=20.0)
    assert reg.shape == (1, 3)
    assert np.isfinite(reg).all()


def test_permute_registration_indices_matches_volume() -> None:
    shape = (671, 671, 341)
    vol = np.zeros(shape, dtype=np.float32)
    vol[24, 326, 239] = 1.0
    permuted = permute_brain_volume(vol, [-1, 3, 2])
    idx = np.argwhere(permuted > 0)[0].astype(float)

    reg = permute_registration_indices_yxz(
        np.array([[24.0, 326.0, 239.0]]),
        shape,
        [-1, 3, 2],
    )[0]
    assert np.allclose(reg, idx, atol=0.5)
