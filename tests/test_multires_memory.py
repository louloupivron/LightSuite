"""Tests for multires overlap memory preflight warnings."""

from __future__ import annotations

from pathlib import Path

import pytest

from lightsuite.multires.memory import (
    estimate_overlap_crop_memory,
    system_memory_info,
    warn_if_overlap_memory_exceeds_system,
)


def test_estimate_overlap_crop_memory() -> None:
    estimate = estimate_overlap_crop_memory([100, 200, 50], max_slab_bytes=100_000_000)
    assert estimate.crop_voxels == 100 * 200 * 50
    assert estimate.one_volume_gb == pytest.approx(estimate.crop_voxels * 4 / 1e9)
    assert estimate.estimated_peak_gb == pytest.approx(2.5 * estimate.one_volume_gb + 0.1)


def test_system_memory_info_parses_meminfo(tmp_path: Path) -> None:
    meminfo = tmp_path / "meminfo"
    meminfo.write_text(
        "\n".join(
            [
                "MemTotal:       33554432 kB",
                "MemFree:         1048576 kB",
                "MemAvailable:   16777216 kB",
            ]
        ),
        encoding="utf-8",
    )
    info = system_memory_info(meminfo_path=meminfo)
    assert info is not None
    assert info.total_gb == pytest.approx(32.0)
    assert info.available_gb == pytest.approx(16.0)


def test_warn_if_overlap_memory_exceeds_available(tmp_path: Path) -> None:
    meminfo = tmp_path / "meminfo"
    meminfo.write_text(
        "\n".join(
            [
                "MemTotal:        8388608 kB",
                "MemAvailable:    2097152 kB",
            ]
        ),
        encoding="utf-8",
    )
    # ~8.6 GB one volume → ~21.6 GB peak, far above 2 GB available.
    with pytest.warns(UserWarning, match="exceeds available RAM"):
        estimate = warn_if_overlap_memory_exceeds_system(
            [1500, 1500, 1000],
            max_slab_bytes=100_000_000,
            meminfo_path=meminfo,
        )
    assert estimate.estimated_peak_gb > 2.0


def test_warn_if_overlap_memory_close_to_available(tmp_path: Path) -> None:
    meminfo = tmp_path / "meminfo"
    meminfo.write_text(
        "\n".join(
            [
                "MemTotal:       33554432 kB",
                "MemAvailable:    4194304 kB",
            ]
        ),
        encoding="utf-8",
    )
    # one volume ~1.34 GB → peak ~3.45 GB, which is >75% of 4 GB available.
    with pytest.warns(UserWarning, match="close to available RAM"):
        warn_if_overlap_memory_exceeds_system(
            [800, 800, 550],
            max_slab_bytes=100_000_000,
            meminfo_path=meminfo,
        )
