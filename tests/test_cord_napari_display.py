"""Tests for spinal cord Napari display alignment."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from lightsuite.gui.cord_napari_display import (
    align_sample_space_for_atlas_qc,
    flip_registration_z_points_xyz,
    flip_registration_z_volume_yxz,
    load_cord_tofliprc,
)


def test_flip_registration_z_points_xyz_is_involution() -> None:
    pts = np.array([[1.0, 2.0, 1.0], [3.0, 4.0, 4.0]])
    flipped = flip_registration_z_points_xyz(pts, nz=4)
    restored = flip_registration_z_points_xyz(flipped, nz=4)
    np.testing.assert_allclose(restored, pts)


def test_align_sample_space_for_atlas_qc_flips_z_when_tofliprc() -> None:
    ann = np.arange(24, dtype=np.int32).reshape(2, 3, 4)
    points = {"cells": np.array([[1.0, 1.0, 1.0], [2.0, 2.0, 4.0]])}
    template, annotation, channels, point_layers, hemisphere = align_sample_space_for_atlas_qc(
        template=ann.astype(np.float32),
        annotation=ann,
        channels={1: ann.astype(np.float32)},
        point_layers=points,
        hemisphere=ann.astype(np.uint8),
        tofliprc=True,
    )
    assert np.array_equal(annotation, flip_registration_z_volume_yxz(ann))
    assert np.array_equal(point_layers["cells"], flip_registration_z_points_xyz(points["cells"], nz=4))
    assert hemisphere is not None
    assert np.array_equal(hemisphere, flip_registration_z_volume_yxz(ann.astype(np.uint8)))
    assert np.array_equal(template, flip_registration_z_volume_yxz(ann.astype(np.float32)))
    assert np.array_equal(channels[1], template)


def test_load_cord_tofliprc_reads_transform_params(tmp_path: Path) -> None:
    save = tmp_path / "registered"
    save.mkdir()
    payload = {
        "tform_bspline_samp20um_to_atlas_20um_px": "",
        "tform_affine_samp20um_to_atlas_20um_px": [
            [1, 0, 0, 0],
            [0, 1, 0, 0],
            [0, 0, 1, 0],
            [0, 0, 0, 1],
        ],
        "control_point_weight": 0.5,
        "samp_ikeeplong": [1, 10],
        "samp_ikeepx": [1, 10],
        "samp_ikeepy": [1, 10],
        "how_to_perm": [1, 2, 3],
        "slicetforms_path": "",
        "sampleres_um": [20.0, 20.0, 20.0],
        "registrationres_um": [20.0, 20.0, 20.0],
        "tofliprc": True,
        "atlassize": [10, 10, 10],
    }
    (save / "transform_params.json").write_text(json.dumps(payload), encoding="utf-8")
    assert load_cord_tofliprc(save) is True
