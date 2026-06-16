"""Tests for mesoSPIM geometry helpers."""

from __future__ import annotations

import numpy as np
import SimpleITK as sitk

from lightsuite.mesospim.config_models import MesospimGeometryConfig
from lightsuite.mesospim.geometry import (
    apply_image_geometry,
    overlap_physical_bounds,
    physical_bounds,
    prepare_registration_pair,
)


def _meta(*, x_pos: float, y_pos: float, z_start: float, z_end: float) -> dict[str, float]:
    return {
        "Pixelsize in um": 5.0,
        "x_pos": x_pos,
        "y_pos": y_pos,
        "z_start": z_start,
        "z_end": z_end,
        "z_stepsize": 2.0,
    }


def _volume(shape_zyx: tuple[int, int, int]) -> sitk.Image:
    arr = np.zeros(shape_zyx, dtype=np.float32)
    return sitk.GetImageFromArray(arr)


def test_overlap_and_prepare_pair() -> None:
    geometry = MesospimGeometryConfig(
        stage_xy_is_center=True,
        itk_lateral_dim0_motor="x",
        lateral_flip=[1, 1],
    )
    meta_a = _meta(x_pos=0.0, y_pos=0.0, z_start=0.0, z_end=10.0)
    meta_b = _meta(x_pos=20.0, y_pos=0.0, z_start=0.0, z_end=10.0)

    fixed = _volume((5, 20, 20))
    moving = _volume((5, 20, 20))
    apply_image_geometry(fixed, meta_a, geometry)
    apply_image_geometry(moving, meta_b, geometry)

    overlap_min, overlap_max = overlap_physical_bounds(fixed, moving, margin_um=0.0)
    assert overlap_min[0] < overlap_max[0]
    assert overlap_min[1] <= overlap_max[1]

    fixed_crop, moving_resampled, box = prepare_registration_pair(fixed, moving, margin_um=0.0)
    assert fixed_crop.GetSize() == moving_resampled.GetSize()
    assert box[0].shape == (3,)
    assert box[1].shape == (3,)

    pmin, pmax = physical_bounds(fixed_crop)
    assert np.all(pmin <= pmax)
