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


def test_discover_compressed_multipage_uses_pages(tmp_path: Path) -> None:
    folder = tmp_path / "sample"
    folder.mkdir()
    slices = [np.full((10, 12), i, dtype=np.uint16) for i in range(3)]
    tifffile.imwrite(
        folder / "stack.tif",
        np.stack(slices, axis=0),
        photometric="minisblack",
        compression="lzw",
        bigtiff=True,
    )

    discovery = discover_tiff_stack(folder, tiff_type=TiffLayout.CHANNEL_PER_FILE)
    assert discovery.nz == 3
    assert discovery.stack_read_mode == "pages"
    assert discovery.use_native_tiff_pages is True

    reader = TiffStackReader.open(folder, tiff_type=TiffLayout.CHANNEL_PER_FILE)
    sl = reader.get_slice(3, channel=0)
    assert sl.shape == (10, 12)
    assert sl[0, 0] == 2


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


def test_volume_layout_from_labeled_axes_zyx_rat_shape() -> None:
    from lightsuite.io.discover import _volume_layout_from_labeled_axes, _volume_layout_from_series

    ny, nx, nz, mode = _volume_layout_from_labeled_axes((2121, 2843, 2086), "ZYX")
    assert (ny, nx, nz, mode) == (2843, 2086, 2121, "memmap_zyx")

    # Without axes metadata, the old heuristic mis-labels Z as the shortest axis.
    ny2, nx2, nz2, mode2 = _volume_layout_from_series((2121, 2843, 2086), None)
    assert (ny2, nx2, nz2, mode2) == (2121, 2843, 2086, "memmap_yxz")


@pytest.mark.skipif(
    not Path("/media/gbm/NVME2/ALICe-pipelines-data/rat/CTR-stitched.tif").is_file(),
    reason="rat sample not mounted",
)
def test_discover_rat_stitched_zyx_stack() -> None:
    discovery = discover_tiff_stack(
        Path("/media/gbm/NVME2/ALICe-pipelines-data/rat"),
        tiff_type=TiffLayout.CHANNEL_PER_FILE,
    )
    assert discovery.stack_read_mode == "memmap_zyx"
    assert discovery.ny == 2843
    assert discovery.nx == 2086
    assert discovery.nz == 2121


def test_planeperfile_rejects_multipage_stack(tmp_path: Path) -> None:
    folder = tmp_path / "sample"
    folder.mkdir()
    _write_multipage_tiff(folder / "stack.tif", [np.zeros((4, 5), dtype=np.uint16) for _ in range(3)])

    with pytest.raises(ValueError, match="channelperfile"):
        discover_tiff_stack(folder, tiff_type=TiffLayout.PLANE_PER_FILE)


def test_discover_planeperfile_multi_channel(tmp_path: Path) -> None:
    ch1 = tmp_path / "channel_488"
    ch2 = tmp_path / "channel_561"
    ch1.mkdir()
    ch2.mkdir()
    for z in range(3):
        plane = (np.arange(12, dtype=np.uint16).reshape(3, 4) + z * 10)
        tifffile.imwrite(ch1 / f"plane_{z:03d}.tif", plane)
        tifffile.imwrite(ch2 / f"plane_{z:03d}.tif", plane + 100)

    discovery = discover_tiff_stack(
        ch1,
        tiff_type=TiffLayout.PLANE_PER_FILE,
        channel_folders=(ch1, ch2),
    )
    assert discovery.nchans == 2
    assert discovery.nz == 3
    assert discovery.channel_plane_files is not None
    assert len(discovery.channel_plane_files) == 2

    reader = TiffStackReader(discovery)
    sl0 = reader.get_slice(2, channel=0)
    sl1 = reader.get_slice(2, channel=1)
    assert sl0[0, 0] == 10
    assert sl1[0, 0] == 110
    reader.close()


def test_discover_planeperfile_multi_channel_mismatched_nz(tmp_path: Path) -> None:
    ch1 = tmp_path / "ch1"
    ch2 = tmp_path / "ch2"
    ch1.mkdir()
    ch2.mkdir()
    tifffile.imwrite(ch1 / "z0.tif", np.zeros((4, 5), dtype=np.uint16))
    tifffile.imwrite(ch1 / "z1.tif", np.zeros((4, 5), dtype=np.uint16))
    tifffile.imwrite(ch2 / "z0.tif", np.zeros((4, 5), dtype=np.uint16))
    tifffile.imwrite(ch2 / "z1.tif", np.zeros((4, 5), dtype=np.uint16))
    tifffile.imwrite(ch2 / "z2.tif", np.zeros((4, 5), dtype=np.uint16))

    with pytest.raises(ValueError, match="mismatched dimensions"):
        discover_tiff_stack(
            ch1,
            tiff_type=TiffLayout.PLANE_PER_FILE,
            channel_folders=(ch1, ch2),
        )


def test_tiff_reader_get_slice(tmp_path: Path) -> None:
    folder = tmp_path / "sample"
    folder.mkdir()
    arr = np.arange(60, dtype=np.uint16).reshape(2, 5, 6)
    _write_multipage_tiff(folder / "volume.tif", [arr[0], arr[1]])

    reader = TiffStackReader.open(folder, tiff_type=TiffLayout.CHANNEL_PER_FILE)
    sl = reader.get_slice(1, channel=0)
    assert sl.shape == (5, 6)
    assert sl[0, 0] == 0
