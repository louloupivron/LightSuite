"""Tests for FIJI point-tool Results.csv conversion."""

from __future__ import annotations

from pathlib import Path

from lightsuite.import_.fiji import convert_fiji_points_to_csv, fiji_to_native_xyz


def test_fiji_um_to_native_xyz() -> None:
    x, y, z = fiji_to_native_xyz(1739.4, 1887.6, 499.0, voxel_um=[2.6, 2.6, 3.0])
    assert x == 670
    assert y == 726
    assert z == 499.0


def test_fiji_pixels_to_native_xyz() -> None:
    x, y, z = fiji_to_native_xyz(100.5, 200.5, 10.0, voxel_um=None)
    assert x == 101.0
    assert y == 201.0
    assert z == 10.0


def test_convert_fiji_points_to_csv(tmp_path: Path) -> None:
    source = tmp_path / "Results.csv"
    source.write_text(
        " ,X,Y,Slice\n"
        "1,1739.400,1887.600,499\n"
        "2,2823.600,3432.000,967\n",
        encoding="utf-8",
    )
    out = tmp_path / "points.csv"
    n = convert_fiji_points_to_csv(source, out, voxel_um=[2.6, 2.6, 3.0])
    assert n == 2
    lines = out.read_text(encoding="utf-8").splitlines()
    assert lines[0] == "x,y,z"
    assert "670,726,499.0" in lines[1]
    assert "1087,1321,967.0" in lines[2]
