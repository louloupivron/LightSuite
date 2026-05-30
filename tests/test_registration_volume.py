"""Tests for registration TIFF loading."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import tifffile

from lightsuite.preprocess.slice_ops import write_z_downsampled_volume
from lightsuite.registration.volume import load_registration_volume


def test_load_registration_volume_stacks_multipage_tiff(tmp_path: Path) -> None:
    """Regression: tifffile may expose each IFD as a separate 2D series."""
    vol = np.arange(24, dtype=np.uint16).reshape(2, 3, 4)
    path = tmp_path / "reg.tif"
    write_z_downsampled_volume(vol, path, scale_z=1.0)

    with tifffile.TiffFile(path) as tif:
        assert len(tif.pages) == 4
        # Writer produces one 2D series per page in current tifffile versions.
        assert len(tif.series) >= 4
        assert tif.asarray().ndim == 2

    loaded = load_registration_volume(path)
    assert loaded.shape == (2, 3, 4)
    assert np.allclose(loaded, vol.astype(np.float32))


def test_load_registration_volume_single_page_is_2d(tmp_path: Path) -> None:
    path = tmp_path / "one.tif"
    tifffile.imwrite(path, np.zeros((5, 6), dtype=np.uint16))
    loaded = load_registration_volume(path)
    assert loaded.shape == (5, 6)
    assert loaded.ndim == 2
