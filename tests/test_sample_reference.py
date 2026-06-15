"""Tests for sample_reference.json publication."""

from __future__ import annotations

import json
from pathlib import Path

from lightsuite.import_.sample_reference import (
    SAMPLE_SPACE_FORMAT,
    sample_reference_path,
    write_sample_reference,
)


def test_write_sample_reference(tmp_path: Path) -> None:
    path = write_sample_reference(
        tmp_path,
        sample_name="mouse1",
        ny=100,
        nx=200,
        nz=50,
        voxel_um=[2.0, 2.5, 3.0],
    )
    assert path == sample_reference_path(tmp_path)
    raw = json.loads(path.read_text(encoding="utf-8"))
    assert raw["format"] == SAMPLE_SPACE_FORMAT
    assert raw["shape_yxz"] == [100, 200, 50]
    assert raw["voxel_um"] == [2.0, 2.5, 3.0]
    assert raw["index_base"] == 1
    assert raw["orientation_applied"] is False
