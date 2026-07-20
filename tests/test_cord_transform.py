"""Tests for spinal cord annotation coordinate transforms."""

from __future__ import annotations

import numpy as np

from lightsuite.import_.cord_transform import (
    crop_cord_registration_points,
    filter_points_for_cord_registration,
    native_xyz_to_cord_registration_yxz,
    permute_cord_registration_indices_yxz,
    registration_atlas_yxz_to_native_cloud_xyz,
)
from lightsuite.preprocess.cord_checkpoint import CordTransformParamsCheckpoint


def _transform_params() -> CordTransformParamsCheckpoint:
    return CordTransformParamsCheckpoint(
        tform_bspline_samp20um_to_atlas_20um_px="/tmp/bspline.txt",
        tform_affine_samp20um_to_atlas_20um_px=np.eye(4).tolist(),
        control_point_weight=0.2,
        samp_ikeeplong=[10, 50],
        samp_ikeepx=[20, 80],
        samp_ikeepy=[30, 90],
        how_to_perm=[1, 2, 3],
        slicetforms_path="/tmp/slicetforms.npy",
        sampleres_um=[1.8, 1.8, 1.8],
        registrationres_um=[20.0, 20.0, 20.0],
        tofliprc=False,
        atlassize=[127, 161, 200],
    )


def test_native_xyz_to_cord_registration_yxz() -> None:
    pts = np.array([[1.0, 1.0, 1.0], [112.0, 112.0, 112.0]])
    reg = native_xyz_to_cord_registration_yxz(
        pts,
        sampleres_um=[1.8, 1.8, 1.8],
        registrationres_um=20.0,
    )
    assert reg.shape == (2, 3)
    assert np.allclose(reg[0], [0.0, 0.0, 0.0])
    assert np.allclose(reg[1], [9.99, 9.99, 9.99], atol=0.01)


def test_crop_cord_registration_points() -> None:
    pts = np.array(
        [
            [29.0, 19.0, 9.0],
            [30.0, 20.0, 10.0],
            [90.0, 80.0, 50.0],
            [91.0, 81.0, 51.0],
        ]
    )
    cropped, keep = crop_cord_registration_points(
        pts,
        yrange=[30, 90],
        xrange=[20, 80],
        zrange=[10, 50],
    )
    assert keep.tolist() == [True, True, False, False]
    assert cropped.shape == (2, 3)
    assert np.allclose(cropped[0], [0.0, 0.0, 0.0])
    assert np.allclose(cropped[1], [1.0, 1.0, 1.0])


def test_filter_points_for_cord_registration() -> None:
    params = _transform_params()
    pts = np.array(
        [
            [1.0, 1.0, 1.0],
            [1000.0, 1000.0, 1000.0],
            [500.0, 500.0, 500.0],
        ]
    )
    filtered, keep = filter_points_for_cord_registration(pts, transform_params=params)
    assert filtered.shape[0] == int(keep.sum())
    assert filtered.shape[0] <= pts.shape[0]


def test_permute_cord_registration_indices_yxz() -> None:
    pts = np.array([[1.0, 2.0, 3.0]])
    permuted = permute_cord_registration_indices_yxz(pts, [2, 3, 1])
    assert np.allclose(permuted, [[2.0, 3.0, 1.0]])


def test_registration_atlas_yxz_to_native_cloud_xyz() -> None:
    pts = np.array([[0.0, 0.0, 0.0], [126.0, 160.0, 199.0]])
    cloud = registration_atlas_yxz_to_native_cloud_xyz(
        pts,
        atlassize=(127, 161, 200),
        template_native_shape=(1567, 253, 322),
    )
    assert cloud.shape == (2, 3)
    assert np.allclose(cloud[0], [1.0, 1.0, 1.0])
    assert np.allclose(cloud[1], [322.0, 253.0, 1567.0])
