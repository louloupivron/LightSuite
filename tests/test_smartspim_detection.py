"""Tests for SmartSPIM / LCT detection points JSON conversion."""

from __future__ import annotations

import json
from pathlib import Path

from lightsuite.import_.smartspim_detection import (
    convert_smartspim_points_json_to_csv,
    smartspim_zyx_to_native_xyz,
)


def test_smartspim_zyx_to_native_xyz() -> None:
    x, y, z = smartspim_zyx_to_native_xyz(105.0, 3260.0, 3225.0)
    assert (x, y, z) == (3226.0, 3261.0, 106.0)


def test_convert_smartspim_points_json_to_csv(tmp_path: Path) -> None:
    source = tmp_path / "488_points_Endogenous.json"
    source.write_text(
        json.dumps([[105.0, 3260.0, 3225.0], [0, 1, 2]]),
        encoding="utf-8",
    )
    out = tmp_path / "points.csv"
    n = convert_smartspim_points_json_to_csv(source, out)
    assert n == 2
    lines = out.read_text(encoding="utf-8").splitlines()
    assert lines[0] == "x,y,z"
    assert lines[1] == "3226.0,3261.0,106.0"
    assert lines[2] == "3.0,2.0,1.0"
