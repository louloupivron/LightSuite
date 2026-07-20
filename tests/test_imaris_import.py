"""Tests for Imaris spot CSV conversion."""

from __future__ import annotations

from pathlib import Path

from lightsuite.import_.imaris import (
    convert_imaris_spots_by_component,
    convert_imaris_spots_to_points_csv,
    list_imaris_component_names,
    read_imaris_spot_csv,
    warn_if_positions_look_like_voxel_indices,
)


def _write_sample_csv(path: Path) -> None:
    path.write_text(
        "\n".join(
            [
                "",
                "Spot",
                " ==================== ",
                "Position X [µm],Position Y [µm],Position Z [µm],Category,Time,ID,Component Name,Original FileName,",
                "1131.019531,2974.889160,473.723511,Spot,1,125094,TAyellow,test.ims,",
                "1209.436890,3028.530273,480.679810,Spot,1,125143,TAyellow,test.ims,",
                "1068.126343,8186.522949,1459.500000,Spot,1,166601,MG_cyan,test.ims,",
            ]
        ),
        encoding="utf-8",
    )


def test_read_imaris_spot_csv(tmp_path: Path) -> None:
    csv_path = tmp_path / "spots.csv"
    _write_sample_csv(csv_path)
    header, rows = read_imaris_spot_csv(csv_path)
    assert "Position X [µm]" in header
    assert len(rows) == 3


def test_convert_imaris_spots_to_points_csv(tmp_path: Path) -> None:
    csv_path = tmp_path / "spots.csv"
    _write_sample_csv(csv_path)
    out_path = tmp_path / "points.csv"
    n = convert_imaris_spots_to_points_csv(
        csv_path,
        out_path,
        voxel_um=[1.8, 1.8, 1.8],
        component_name="TAyellow",
    )
    assert n == 2
    text = out_path.read_text(encoding="utf-8")
    assert "x,y,z" in text.splitlines()[0]
    assert "1653" in text or "1654" in text


def test_warn_if_positions_look_like_voxel_indices(tmp_path: Path) -> None:
    csv_path = tmp_path / "spots.csv"
    _write_sample_csv(csv_path)
    # Positions ~1e3–8e3 look like indices on a ~8793×2004×1931 grid, not µm at 1.8
    warning = warn_if_positions_look_like_voxel_indices(
        csv_path,
        voxel_um=[1.8, 1.8, 1.8],
        shape_yxz=(8793, 2004, 1931),
    )
    assert warning is not None
    assert "1,1,1" in warning
    assert (
        warn_if_positions_look_like_voxel_indices(
            csv_path,
            voxel_um=[1.0, 1.0, 1.0],
            shape_yxz=(8793, 2004, 1931),
        )
        is None
    )


def test_convert_imaris_spots_by_component(tmp_path: Path) -> None:
    csv_path = tmp_path / "spots.csv"
    _write_sample_csv(csv_path)
    components = list_imaris_component_names(csv_path)
    assert components == ["MG_cyan", "TAyellow"]
    written = convert_imaris_spots_by_component(
        csv_path,
        tmp_path / "converted",
        voxel_um=[1.8, 1.8, 1.8],
    )
    assert set(written) == {"MG_cyan", "TAyellow"}
    assert all(path.is_file() for path in written.values())
