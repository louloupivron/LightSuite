"""Tests for affine fit quality metrics."""

from __future__ import annotations

import numpy as np

from lightsuite.gui.affine import (
    affine_point_errors,
    fit_affine_transform,
    summarize_point_errors,
    transform_points,
    transform_points_inverse,
)
from lightsuite.registration.brain_register import _build_affine_fit_diagnostics
from lightsuite.registration.points import cloud_xyz_to_volume_indices, volume_indices_to_cloud_xyz


def test_affine_point_errors_identity() -> None:
    pts = np.array([[0.0, 0.0, 0.0], [1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]])
    mse, err = affine_point_errors(pts, pts, np.eye(4))
    assert mse == 0.0
    assert np.allclose(err, 0.0)


def test_build_affine_fit_diagnostics_auto_only() -> None:
    atlas = np.array(
        [
            [10.0, 20.0, 30.0],
            [12.0, 22.0, 32.0],
            [14.0, 24.0, 34.0],
            [16.0, 26.0, 36.0],
        ]
    )
    sample = atlas + np.array([1.0, 2.0, 0.5])
    tform, _ = fit_affine_transform(atlas, sample)
    diag = _build_affine_fit_diagnostics(
        af_atlas=atlas,
        af_sample=sample,
        tform_aff=tform,
        n_manual=0,
        n_auto=4,
        autocpsample_kept=sample.copy(),
        original_trans_vol=np.eye(4),
        cpaffine=transform_points(atlas, tform),
        cptshistology=sample,
        cpwt=0.1,
    )
    assert diag.n_total == 4
    assert diag.median_error_vox < 1.0
    assert diag.median_landmark_vox is not None
    assert diag.median_landmark_vox < 1.0


def test_coarse_auto_error_uses_forward_sample_to_atlas_transform() -> None:
    # original_trans_vol maps native sample -> atlas space and is already in volume
    # (Y, X, Z) axes. Inverting it, or re-applying the involutive swap_xy_transform,
    # reports a coarse error two orders of magnitude too large.
    atlas = np.array(
        [
            [10.0, 20.0, 30.0],
            [12.0, 22.0, 32.0],
            [14.0, 24.0, 34.0],
            [16.0, 26.0, 36.0],
        ]
    )
    original_trans_vol = np.eye(4)
    original_trans_vol[:3, 3] = [5.0, 0.0, 0.0]  # Y-only shift: asymmetric under X/Y swap
    # Auto sample points are stored in native sample space; forward-mapping them must
    # land on their atlas partners.
    sample = transform_points_inverse(atlas, original_trans_vol)
    tform, _ = fit_affine_transform(atlas, sample)

    diag = _build_affine_fit_diagnostics(
        af_atlas=atlas,
        af_sample=sample,
        tform_aff=tform,
        n_manual=0,
        n_auto=4,
        autocpsample_kept=sample.copy(),
        original_trans_vol=original_trans_vol,
        cpaffine=transform_points(atlas, tform),
        cptshistology=sample,
        cpwt=0.1,
    )
    assert diag.median_coarse_auto_vox is not None
    assert diag.median_coarse_auto_vox < 1e-6


def test_cloud_xyz_volume_index_roundtrip() -> None:
    xyz = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
    vol = cloud_xyz_to_volume_indices(xyz)
    assert np.allclose(vol, np.array([[2.0, 1.0, 3.0], [5.0, 4.0, 6.0]]))
    assert np.allclose(volume_indices_to_cloud_xyz(vol), xyz)


def test_summarize_point_errors_empty() -> None:
    summary = summarize_point_errors(np.array([]))
    assert summary["median"] == 0.0
