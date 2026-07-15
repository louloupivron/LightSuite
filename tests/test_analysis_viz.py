"""Tests for analysis visualization (matplotlib Agg backend)."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import numpy as np
import pandas as pd
import pytest

from lightsuite.analysis.viz.io import (
    aggregate_by_division,
    load_region_plot_table,
    select_top_regions,
)
from lightsuite.analysis.viz.plots import plot_division_bars, plot_lr_scatter, plot_top_region_bars
from lightsuite.export.parcellation import ParcellationResult, write_parcellation_csv


def _tidy_rows() -> pd.DataFrame:
    rows = []
    spec = [
        (1, "Isocortex", "right", 100.0),
        (1, "Isocortex", "left", 110.0),
        (2, "Thalamus", "right", 50.0),
        (2, "Thalamus", "left", 60.0),
    ]
    for pidx, div, hemi, val in spec:
        rows.append(
            {
                "sample": "m1",
                "channel": 1,
                "atlas": "allen",
                "parcellation_index": pidx,
                "acronym": f"R{pidx}",
                "name": f"Region {pidx}",
                "structure": "MO",
                "division": div,
                "hemisphere": hemi,
                "metric": "median_intensity",
                "value": val,
            }
        )
    return pd.DataFrame(rows)


def test_load_region_plot_table_from_tidy_import_label(tmp_path: Path) -> None:
    path = tmp_path / "region_stats.csv"
    df = _tidy_rows()
    df["channel"] = "imaris_488_cells"
    df["metric"] = "cell_count"
    df.to_csv(path, index=False)
    table = load_region_plot_table(path, channel="imaris_488_cells", metric="cell_count")
    assert len(table) == 2


def test_parse_plot_channel() -> None:
    from lightsuite.analysis.viz.io import parse_plot_channel

    assert parse_plot_channel("1") == 1
    assert parse_plot_channel("imaris_488_cells") == "imaris_488_cells"


def test_load_region_plot_table_from_tidy(tmp_path: Path) -> None:
    path = tmp_path / "region_stats.csv"
    _tidy_rows().to_csv(path, index=False)
    table = load_region_plot_table(path, channel=1, metric="median_intensity")
    assert list(table.columns) == ["parcellation_index", "name", "structure", "division", "left", "right"]
    assert len(table) == 2
    row1 = table.set_index("parcellation_index").loc[1]
    assert row1["left"] == 110.0
    assert row1["right"] == 100.0


def test_load_region_plot_table_from_wide(tmp_path: Path) -> None:
    result = ParcellationResult(
        area_ids=np.array([1, 2], dtype=np.int64),
        median_over_areas=pd.DataFrame([[10.0, 20.0], [30.0, 40.0]]).to_numpy(dtype="float32"),
        std_over_areas=pd.DataFrame([[0.0, 0.0], [0.0, 0.0]]).to_numpy(dtype="float32"),
        volume_over_areas=pd.DataFrame([[1.0, 1.0], [1.0, 1.0]]).to_numpy(dtype="float32"),
    )
    path = tmp_path / "intensities.csv"
    write_parcellation_csv(path, result)
    meta = pd.DataFrame(
        {
            "parcellation_index": [1, 2],
            "name": ["A", "B"],
            "structure": ["MO", "TH"],
            "division": ["Isocortex", "Thalamus"],
        }
    )
    wide = pd.read_csv(path).merge(meta, on="parcellation_index")
    wide.to_csv(path, index=False)

    table = load_region_plot_table(path, metric="median_intensity")
    assert table.loc[0, "left"] == 20.0
    assert table.loc[0, "right"] == 10.0


def test_aggregate_by_division_means_regions() -> None:
    df = pd.DataFrame(
        [
            {"division": "Isocortex", "left": 10.0, "right": 20.0},
            {"division": "Isocortex", "left": 30.0, "right": 40.0},
            {"division": "Thalamus", "left": 5.0, "right": 7.0},
        ]
    )
    agg = aggregate_by_division(df, how="mean")
    iso = agg[agg["division"] == "Isocortex"].iloc[0]
    assert iso["left"] == 20.0
    assert iso["right"] == 30.0
    assert iso["count"] == 2


def test_aggregate_by_division_sums_regions() -> None:
    from lightsuite.analysis.viz.io import default_division_aggregate

    df = pd.DataFrame(
        [
            {"division": "Medulla", "left": 6.0, "right": 115.0},
            {"division": "Medulla", "left": 6.0, "right": 15.0},
        ]
    )
    agg = aggregate_by_division(df, how="sum")
    row = agg.iloc[0]
    assert row["left"] == 12.0
    assert row["right"] == 130.0
    assert row["total"] == 142.0
    assert default_division_aggregate("cell_count") == "sum"
    assert default_division_aggregate("median_intensity") == "mean"


def test_plot_division_bars_writes_png(tmp_path: Path) -> None:
    path = tmp_path / "region_stats.csv"
    _tidy_rows().to_csv(path, index=False)
    table = load_region_plot_table(path, channel=1)
    out = tmp_path / "division_bars.png"
    plot_division_bars(table, output_path=out)
    assert out.is_file()
    assert out.with_suffix(".csv").is_file()


def test_plot_lr_scatter_writes_png(tmp_path: Path) -> None:
    path = tmp_path / "region_stats.csv"
    _tidy_rows().to_csv(path, index=False)
    table = load_region_plot_table(path, channel=1)
    out = tmp_path / "lr_scatter.png"
    plot_lr_scatter(table, output_path=out, keep_divisions=["Isocortex", "Thalamus"])
    assert out.is_file()


def test_select_top_regions_ranks_by_total() -> None:
    df = pd.DataFrame(
        [
            {"parcellation_index": 1, "name": "A", "structure": "MO", "division": "Isocortex", "left": 10.0, "right": 5.0},
            {"parcellation_index": 2, "name": "B", "structure": "TH", "division": "Thalamus", "left": 50.0, "right": 40.0},
            {"parcellation_index": 3, "name": "C", "structure": "HY", "division": "Hypothalamus", "left": 30.0, "right": 20.0},
        ]
    )
    top = select_top_regions(df, top_n=2)
    assert list(top["parcellation_index"]) == [3, 2]
    assert top.iloc[-1]["total"] == 90.0


def test_plot_top_region_bars_writes_png(tmp_path: Path) -> None:
    path = tmp_path / "region_stats.csv"
    _tidy_rows().to_csv(path, index=False)
    table = load_region_plot_table(path, channel=1)
    out = tmp_path / "top_regions.png"
    plot_top_region_bars(table, top_n=2, metric="median_intensity", output_path=out)
    assert out.is_file()
    assert out.with_suffix(".csv").is_file()


def test_plot_group_division_bars(tmp_path: Path) -> None:
    from lightsuite.analysis.viz.cohort_plots import plot_group_division_bars

    summary = pd.DataFrame(
        {
            "group": ["control", "control", "treat", "treat"],
            "division": ["Isocortex", "Thalamus", "Isocortex", "Thalamus"],
            "mean": [10.0, 5.0, 20.0, 8.0],
            "sem": [1.0, 0.5, 2.0, 1.0],
        }
    )
    out = tmp_path / "group_div.png"
    plot_group_division_bars(summary, output_path=out)
    assert out.is_file()


def test_load_tidy_missing_channel_raises(tmp_path: Path) -> None:
    path = tmp_path / "region_stats.csv"
    _tidy_rows().to_csv(path, index=False)
    with pytest.raises(ValueError, match="No rows"):
        load_region_plot_table(path, channel=99, metric="median_intensity")
