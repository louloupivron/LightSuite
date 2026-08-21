"""Tests for SmartSPIM All_Channels → ChN split."""

from __future__ import annotations

from pathlib import Path

import pytest

from lightsuite.io.smartspim_channels import (
    SplitMode,
    discover_smartspim_channel_planes,
    parse_smartspim_channel_plane,
    split_smartspim_all_channels,
)


def _touch(path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(b"tiff")
    return path


def test_parse_smartspim_channel_plane() -> None:
    assert parse_smartspim_channel_plane(Path("Z000000_Ch0.tif")) == (0, "Z000000")
    assert parse_smartspim_channel_plane(Path("415_463_000040_Ch1.tiff")) == (
        1,
        "415_463_000040",
    )
    assert parse_smartspim_channel_plane(Path("readme.txt")) is None
    assert parse_smartspim_channel_plane(Path("Z000000.tif")) is None


def test_discover_and_split_symlink(tmp_path: Path) -> None:
    src = tmp_path / "All_Channels"
    for z in (0, 10, 20):
        for ch in (0, 1):
            _touch(src / f"Z{z:06d}_Ch{ch}.tif")
    _touch(src / "ImageList.txt")  # ignored

    by_ch = discover_smartspim_channel_planes(src)
    assert list(by_ch.keys()) == [0, 1]
    assert len(by_ch[0]) == 3

    result = split_smartspim_all_channels(src, mode=SplitMode.SYMLINK)
    assert result.plane_counts == {0: 3, 1: 3}
    assert (result.channel_dirs[0] / "Z000000_Ch0.tif").is_symlink()
    assert (result.channel_dirs[1] / "Z000010_Ch1.tif").resolve() == (
        src / "Z000010_Ch1.tif"
    ).resolve()
    # Originals remain for symlink mode
    assert (src / "Z000000_Ch0.tif").is_file()


def test_split_move_and_unequal_counts(tmp_path: Path) -> None:
    src = tmp_path / "All_Channels"
    _touch(src / "Z000000_Ch0.tif")
    _touch(src / "Z000000_Ch1.tif")
    _touch(src / "Z000010_Ch0.tif")

    with pytest.raises(ValueError, match="Unequal plane counts"):
        split_smartspim_all_channels(src)

    # Allow unequal when explicitly requested
    result = split_smartspim_all_channels(
        src, mode="move", require_equal_plane_counts=False
    )
    assert result.plane_counts == {0: 2, 1: 1}
    assert not (src / "Z000000_Ch0.tif").exists()
    assert (result.channel_dirs[0] / "Z000000_Ch0.tif").is_file()


def test_reject_channel_above_max(tmp_path: Path) -> None:
    src = tmp_path / "All_Channels"
    _touch(src / "Z000000_Ch0.tif")
    _touch(src / "Z000000_Ch3.tif")
    with pytest.raises(ValueError, match="channels 0–2"):
        discover_smartspim_channel_planes(src, max_channels=3)


def test_split_to_external_output(tmp_path: Path) -> None:
    src = tmp_path / "All_Channels"
    out = tmp_path / "by_channel"
    for ch in (0, 1, 2):
        _touch(src / f"Z000000_Ch{ch}.tif")
    result = split_smartspim_all_channels(src, output_dir=out, mode=SplitMode.SYMLINK)
    assert set(result.channel_dirs) == {0, 1, 2}
    assert (out / "Ch2" / "Z000000_Ch2.tif").is_symlink()
    assert (out / "Ch2" / "Z000000_Ch2.tif").resolve() == (src / "Z000000_Ch2.tif").resolve()
    assert (src / "Z000000_Ch2.tif").is_file()


def test_split_rejects_copy_mode(tmp_path: Path) -> None:
    src = tmp_path / "All_Channels"
    _touch(src / "Z000000_Ch0.tif")
    with pytest.raises(ValueError, match="copy mode is not supported"):
        split_smartspim_all_channels(src, mode="copy")


def test_discover_rejects_interleaved_smartspim_planes(tmp_path: Path) -> None:
    from lightsuite.config.models import TiffLayout
    from lightsuite.io.discover import discover_tiff_stack

    src = tmp_path / "All_Channels"
    for z in (0, 10):
        for ch in (0, 1):
            _touch(src / f"Z{z:06d}_Ch{ch}.tif")
    with pytest.raises(ValueError, match="interleaved"):
        discover_tiff_stack(src, tiff_type=TiffLayout.PLANE_PER_FILE)
