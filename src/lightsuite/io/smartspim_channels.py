"""Split SmartSPIM flat ``All_Channels`` stacks into per-channel plane folders.

SmartSPIM often writes interleaved planes in one directory::

    All_Channels/Z000000_Ch0.tif
    All_Channels/Z000000_Ch1.tif
    All_Channels/Z000010_Ch0.tif
    ...

LightSuite ``planeperfile`` multi-channel configs expect one folder per channel
(``Ch0/``, ``Ch1/``, …). This helper creates those folders with symlinks (default),
hardlinks, copies, or moves.
"""

from __future__ import annotations

import re
import shutil
from collections import defaultdict
from dataclasses import dataclass
from enum import Enum
from pathlib import Path

# Cap at three channels (Ch0–Ch2) for typical SmartSPIM dual/triple acquisitions.
MAX_CHANNELS = 3

_CHANNEL_RE = re.compile(
    r"^(?P<stem>.+)_Ch(?P<channel>\d+)\.(?P<ext>tif{1,2})$",
    re.IGNORECASE,
)


class SplitMode(str, Enum):
    SYMLINK = "symlink"
    HARDLINK = "hardlink"
    COPY = "copy"
    MOVE = "move"


@dataclass(frozen=True)
class ChannelSplitResult:
    """Outcome of splitting a flat SmartSPIM ``All_Channels`` folder."""

    source_dir: Path
    output_dir: Path
    channel_dirs: dict[int, Path]
    plane_counts: dict[int, int]
    mode: SplitMode


def parse_smartspim_channel_plane(path: Path) -> tuple[int, str] | None:
    """Return ``(channel_index, basename_stem_without_ChN)`` or ``None`` if not matching."""
    match = _CHANNEL_RE.match(path.name)
    if match is None:
        return None
    channel = int(match.group("channel"))
    return channel, match.group("stem")


def discover_smartspim_channel_planes(
    source_dir: Path,
    *,
    max_channels: int = MAX_CHANNELS,
) -> dict[int, list[Path]]:
    """Group TIFF planes under ``source_dir`` by channel index (0-based).

    Only files matching ``*_ChN.tif`` / ``*_ChN.tiff`` are included. Channels
    ``>= max_channels`` raise ``ValueError``. Non-matching TIFFs are ignored.
    """
    source_dir = source_dir.expanduser().resolve()
    if not source_dir.is_dir():
        msg = f"SmartSPIM All_Channels folder not found: {source_dir}"
        raise FileNotFoundError(msg)

    by_channel: dict[int, list[Path]] = defaultdict(list)
    for path in sorted(source_dir.iterdir()):
        if not path.is_file():
            continue
        if path.suffix.lower() not in {".tif", ".tiff"}:
            continue
        parsed = parse_smartspim_channel_plane(path)
        if parsed is None:
            continue
        channel, _stem = parsed
        if channel < 0:
            msg = f"Negative channel index in {path.name}"
            raise ValueError(msg)
        if channel >= max_channels:
            msg = (
                f"Found channel index {channel} in {path.name}, but only "
                f"channels 0–{max_channels - 1} are supported (max_channels={max_channels})."
            )
            raise ValueError(msg)
        by_channel[channel].append(path)

    if not by_channel:
        msg = (
            f"No SmartSPIM plane TIFFs matching '*_ChN.tif' found in {source_dir}. "
            "Expected names like Z000000_Ch0.tif or …_000040_Ch1.tif."
        )
        raise FileNotFoundError(msg)

    # Stable Z order within each channel (filename sort matches Z##### / stage keys).
    return {ch: sorted(paths, key=lambda p: p.name.lower()) for ch, paths in sorted(by_channel.items())}


def _place_file(src: Path, dest: Path, mode: SplitMode) -> None:
    if dest.exists() or dest.is_symlink():
        dest.unlink()
    if mode == SplitMode.SYMLINK:
        dest.symlink_to(src.resolve())
    elif mode == SplitMode.HARDLINK:
        dest.hardlink_to(src)
    elif mode == SplitMode.COPY:
        shutil.copy2(src, dest)
    elif mode == SplitMode.MOVE:
        shutil.move(str(src), str(dest))
    else:
        msg = f"Unsupported split mode: {mode}"
        raise ValueError(msg)


def split_smartspim_all_channels(
    source_dir: Path,
    *,
    output_dir: Path | None = None,
    mode: SplitMode | str = SplitMode.SYMLINK,
    max_channels: int = MAX_CHANNELS,
    require_equal_plane_counts: bool = True,
) -> ChannelSplitResult:
    """Create ``Ch0`` / ``Ch1`` / ``Ch2`` folders for LightSuite ``planeperfile``.

    Parameters
    ----------
    source_dir:
        Flat SmartSPIM ``All_Channels`` directory.
    output_dir:
        Parent for ``ChN`` folders. Defaults to ``source_dir`` (in-place subfolders).
    mode:
        ``symlink`` (default), ``hardlink``, ``copy``, or ``move``.
    max_channels:
        Highest allowed channel index + 1 (default 3 → Ch0–Ch2).
    require_equal_plane_counts:
        If True, raise when channels have different plane counts.
    """
    if max_channels < 1 or max_channels > MAX_CHANNELS:
        msg = f"max_channels must be between 1 and {MAX_CHANNELS}, got {max_channels}"
        raise ValueError(msg)

    mode_enum = SplitMode(mode) if not isinstance(mode, SplitMode) else mode
    source_dir = source_dir.expanduser().resolve()
    out = (output_dir or source_dir).expanduser().resolve()
    out.mkdir(parents=True, exist_ok=True)

    by_channel = discover_smartspim_channel_planes(source_dir, max_channels=max_channels)
    counts = {ch: len(paths) for ch, paths in by_channel.items()}
    if require_equal_plane_counts and len(set(counts.values())) > 1:
        msg = f"Unequal plane counts across channels: {counts}"
        raise ValueError(msg)

    channel_dirs: dict[int, Path] = {}
    for channel, paths in by_channel.items():
        ch_dir = out / f"Ch{channel}"
        ch_dir.mkdir(parents=True, exist_ok=True)
        for src in paths:
            dest = ch_dir / src.name
            _place_file(src, dest, mode_enum)
        channel_dirs[channel] = ch_dir

    return ChannelSplitResult(
        source_dir=source_dir,
        output_dir=out,
        channel_dirs=channel_dirs,
        plane_counts=counts,
        mode=mode_enum,
    )
