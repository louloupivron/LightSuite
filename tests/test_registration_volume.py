"""Tests for registration TIFF loading."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import tifffile

from lightsuite.preprocess.slice_ops import write_z_downsampled_volume
from lightsuite.registration.volume import (
    load_registration_volume,
    load_tiff_volume_zyx,
    normalize_registration_volume,
)


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


def test_load_tiff_volume_zyx_multipage_canvas_style(tmp_path: Path) -> None:
    """Plane-by-plane overview canvases are multi-page 2D stacks (not shaped ZYX series)."""
    vol = np.arange(210, dtype=np.float32).reshape(5, 6, 7)
    path = tmp_path / "canvas.tif"
    with tifffile.TiffWriter(path, bigtiff=True) as tif:
        for iz in range(5):
            metadata = {"axes": "ZYX", "spacing": 2.0} if iz == 0 else None
            tif.write(vol[iz], compression="zlib", metadata=metadata)

    loaded_zyx = load_tiff_volume_zyx(path)
    assert loaded_zyx.shape == (5, 6, 7)
    assert np.allclose(loaded_zyx, vol)

    loaded_yxz = load_registration_volume(path)
    assert loaded_yxz.shape == (6, 7, 5)
    assert np.allclose(loaded_yxz, np.moveaxis(vol, 0, -1))

    import contextlib
    import io

    buf = io.StringIO()
    with contextlib.redirect_stderr(buf):
        load_registration_volume(path)
    assert "shaped series" not in buf.getvalue()


def test_load_registration_volume_single_page_is_2d(tmp_path: Path) -> None:
    path = tmp_path / "one.tif"
    tifffile.imwrite(path, np.zeros((5, 6), dtype=np.uint16))
    loaded = load_registration_volume(path)
    assert loaded.shape == (5, 6)
    assert loaded.ndim == 2


def test_normalize_registration_volume_uses_center_diagonal() -> None:
    volume = np.zeros((11, 11, 11), dtype=np.float32)
    volume[5, 5, 5] = 100.0
    normalized = normalize_registration_volume(volume)
    assert np.isclose(normalized[5, 5, 5], 0.5, rtol=0.05)
