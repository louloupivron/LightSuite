"""Tests for Fiederling atlas resampling and native upsampling."""

from __future__ import annotations

import numpy as np
import pandas as pd

from lightsuite.atlas.fiederling import FiederlingAtlasVolumes, resize_fiederling_atlas, upsample_to_fiederling_native


def _make_volumes() -> FiederlingAtlasVolumes:
    """Synthetic atlas with native SC_P56 layout (length axis first, 20 um)."""
    rng = np.random.default_rng(0)
    template = rng.random((40, 12, 16)).astype(np.float32)
    annotation = np.zeros((40, 12, 16), dtype=np.uint16)
    annotation[:, 4:8, 6:10] = 3
    return FiederlingAtlasVolumes(
        template=template,
        annotation=annotation,
        atlas_res_um=(20.0, 10.0, 10.0),
        segments=pd.DataFrame(),
        regions=pd.DataFrame(),
    )


def test_resize_moves_length_axis_last() -> None:
    tv, av = resize_fiederling_atlas(_make_volumes(), 20.0)
    assert tv.shape == (6, 8, 40)
    assert av.shape == (6, 8, 40)


def test_resize_length_axis_fully_occupied() -> None:
    _, av = resize_fiederling_atlas(_make_volumes(), 20.0)
    mask = av > 0
    occupied = mask.any(axis=(0, 1))
    assert occupied.mean() == 1.0
    assert mask.any(axis=(1, 2)).mean() < 1.0
    assert mask.any(axis=(0, 2)).mean() < 1.0


def test_upsample_to_fiederling_native_transposes_length_axis() -> None:
    volume = np.zeros((12, 16, 40), dtype=np.uint16)
    volume[5, 7, :] = 100
    native = np.zeros((40, 24, 32), dtype=np.uint8)
    out = upsample_to_fiederling_native(volume, native)
    assert out.shape == native.shape
    assert np.any(out[10:30, 9, 13] > 0)
