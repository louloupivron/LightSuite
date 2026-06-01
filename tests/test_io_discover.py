"""Tests for TIFF stack discovery."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
import tifffile

from lightsuite.config.models import TiffLayout
from lightsuite.io.discover import discover_tiff_stack
from lightsuite.io.readers.tiff_stack import TiffStackReader


def _write_multipage_tiff(path: Path, arrays: list[np.ndarray]) -> None:
    tifffile.imwrite(path, np.stack(arrays, axis=0), photometric="minisblack")


def test_discover_channel_per_file_multi_tiff(tmp_path: Path) -> None:
    folder = tmp_path / "sample"
    folder.mkdir()
    slice_a = np.ones((10, 12), dtype=np.uint16) * 100
    slice_b = np.ones((10, 12), dtype=np.uint16) * 200
    _write_multipage_tiff(folder / "ch1.tif", [slice_a, slice_b])
    _write_multipage_tiff(folder / "ch2.tif", [slice_a + 1, slice_b + 1])

    discovery = discover_tiff_stack(folder, tiff_type=TiffLayout.CHANNEL_PER_FILE)
    assert discovery.nchans == 2
    assert discovery.nz == 2
    assert discovery.ny == 10
    assert discovery.nx == 12
    assert discovery.multitiffs is True


def test_discover_imagej_hyperstack_metadata(tmp_path: Path) -> None:
    folder = tmp_path / "sample"
    folder.mkdir()
    vol = np.arange(60, dtype=np.uint16).reshape(3, 4, 5)
    tifffile.imwrite(folder / "1-3-1x.tif", vol, imagej=True)

    discovery = discover_tiff_stack(folder, tiff_type=TiffLayout.CHANNEL_PER_FILE)
    assert discovery.nz == 3
    assert discovery.ny == 4
    assert discovery.nx == 5

    reader = TiffStackReader.open(folder, tiff_type=TiffLayout.CHANNEL_PER_FILE)
    sl = reader.get_slice(2, channel=0)
    assert sl.shape == (4, 5)
    assert sl[0, 0] == 20


def test_discover_single_ifd_volumetric(tmp_path: Path) -> None:
    folder = tmp_path / "sample"
    folder.mkdir()
    vol = np.arange(60, dtype=np.uint16).reshape(3, 4, 5)
    tifffile.imwrite(folder / "stack.tif", vol)

    discovery = discover_tiff_stack(folder, tiff_type=TiffLayout.CHANNEL_PER_FILE)
    assert discovery.nz == 3
    assert discovery.stack_read_mode == "memmap_zyx"
    assert discovery.ny == 4
    assert discovery.nx == 5


def test_planeperfile_rejects_multipage_stack(tmp_path: Path) -> None:
    folder = tmp_path / "sample"
    folder.mkdir()
    _write_multipage_tiff(folder / "stack.tif", [np.zeros((4, 5), dtype=np.uint16) for _ in range(3)])

    with pytest.raises(ValueError, match="channelperfile"):
        discover_tiff_stack(folder, tiff_type=TiffLayout.PLANE_PER_FILE)


def test_tiff_reader_get_slice(tmp_path: Path) -> None:
    folder = tmp_path / "sample"
    folder.mkdir()
    arr = np.arange(60, dtype=np.uint16).reshape(2, 5, 6)
    _write_multipage_tiff(folder / "volume.tif", [arr[0], arr[1]])

    reader = TiffStackReader.open(folder, tiff_type=TiffLayout.CHANNEL_PER_FILE)
    sl = reader.get_slice(1, channel=0)
    assert sl.shape == (5, 6)
    assert sl[0, 0] == 0
