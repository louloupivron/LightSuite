"""Tests for foreground bounding boxes and registration canvas."""

from __future__ import annotations

import numpy as np

from lightsuite.config.models import RegistrationCanvasMode
from lightsuite.registration.canvas import compute_registration_canvas
from lightsuite.registration.content_bbox import (
    ContentBox,
    bbox_from_annotation,
    crop_volume_yxz,
    sample_foreground_bbox,
)
from lightsuite.registration.coordinates import affine_with_source_offset


def test_bbox_from_annotation_with_margin() -> None:
    ann = np.zeros((50, 60, 40), dtype=np.uint16)
    ann[10:30, 15:45, 5:25] = 1
    box = bbox_from_annotation(ann, margin_vox=2)
    assert box is not None
    assert box.y0 == 8
    assert box.x0 == 13
    assert box.z0 == 3


def test_sample_foreground_bbox_largest_component() -> None:
    vol = np.zeros((80, 80, 40), dtype=np.uint16)
    vol[20:60, 25:55, 10:30] = 1000
    vol[0:5, 0:5, 0:5] = 500
    box = sample_foreground_bbox(vol, margin_vox=1, trim_z=False)
    assert box is not None
    assert 18 <= box.y0 <= 20
    assert 59 <= box.y1 <= 61


def test_affine_with_source_offset_shifts_translation() -> None:
    tform = np.eye(4)
    tform[:3, 3] = [1.0, 2.0, 3.0]
    shifted = affine_with_source_offset(tform, (10, 0, 0))
    assert shifted[0, 3] == 11.0


def test_compute_registration_canvas_pad() -> None:
    canvas = compute_registration_canvas((10, 20, 30), (10, 25, 30), RegistrationCanvasMode.PAD)
    assert canvas.working_shape == (10, 25, 30)
    assert canvas.pad_before[1] + canvas.pad_after[1] == 5


def test_crop_volume_round_trip() -> None:
    vol = np.arange(24, dtype=np.float32).reshape(2, 3, 4)
    box = ContentBox(0, 1, 1, 2, 0, 3)
    cropped = crop_volume_yxz(vol, box)
    assert cropped.shape == box.size_yxz
