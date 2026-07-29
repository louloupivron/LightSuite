"""Tests for spinal cord region-stats visualization."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import pandas as pd
import pytest

from lightsuite.analysis.viz.cord_io import (
    filter_cord_stats,
    filter_cord_stats_multi,
    load_cord_stats_csv,
    segment_grouped_totals_table,
    structure_heatmap_matrix,
    top_regions_table,
)
from lightsuite.analysis.viz.cord_plots import (
    plot_cord_coloc_overlap,
    plot_cord_division_profile,
    plot_cord_segment_bars,
    plot_cord_segment_grouped_bars,
    plot_cord_structure_heatmap,
    plot_cord_structure_panel,
    plot_cord_top_regions,
)


def _cord_stats_rows() -> pd.DataFrame:
    rows = []
    structures = [
        (201, "Lamina_I", "Lamina I Combined", "C1", "structure", 10.0),
        (201, "Lamina_I", "Lamina I Combined", "C2", "structure", 20.0),
        (71, "GM", "Gray Matter", "C1", "division", 15.0),
        (71, "GM", "Gray Matter", "C2", "division", 25.0),
        (130, "WM", "White matter", "C1", "division", 5.0),
        (130, "WM", "White matter", "C2", "division", 8.0),
        (7, "5Sp", "Lamina 5", "C1", "region", 3.0),
        (7, "5Sp", "Lamina 5", "C2", "region", 4.0),
        (8, "5SpL", "Lamina 5 Lateral", "C1", "region", 2.0),
    ]
    for pidx, acr, name, seg, level, val in structures:
        rows.append(
            {
                "sample": "op87",
                "channel": 1,
                "atlas": "fiederling",
                "parcellation_index": pidx,
                "acronym": acr,
                "name": name,
                "structure": "DH",
                "division": "GM",
                "segment": seg,
                "rollup_level": level,
                "hemisphere": "whole",
                "metric": "median_intensity" if level != "region" else "cell_count",
                "value": val,
            }
        )
    return pd.DataFrame(rows)


def _segments() -> pd.DataFrame:
    return pd.DataFrame({"Segment": ["C1", "C2"], "Start": [1, 15], "End": [14, 42]})


def test_filter_cord_stats_structure_level() -> None:
    df = _cord_stats_rows()
    sub = filter_cord_stats(df, channel=1, metric="median_intensity", rollup_level="structure")
    assert len(sub) == 2
    assert set(sub["segment"]) == {"C1", "C2"}


def test_structure_heatmap_matrix_shape() -> None:
    sub = filter_cord_stats(_cord_stats_rows(), channel=1, metric="median_intensity", rollup_level="structure")
    matrix, row_labels, col_labels = structure_heatmap_matrix(sub, segment_order=["C1", "C2"])
    assert matrix.shape == (1, 2)
    assert row_labels == ["Lamina I Combined"]
    assert col_labels == ["C1", "C2"]
    assert matrix.isna().sum().sum() == 0


def test_structure_heatmap_matrix_keeps_nan_and_crops() -> None:
    rows = _cord_stats_rows()
    # Add an empty-leading segment with no structure rows, and a gap segment via reindex.
    sub = filter_cord_stats(rows, channel=1, metric="median_intensity", rollup_level="structure")
    matrix, _rows, cols = structure_heatmap_matrix(
        sub,
        segment_order=["T1", "C1", "C2", "S1"],
        crop_empty_segments=True,
    )
    assert cols == ["C1", "C2"]
    # Missing T1/S1 dropped; existing values preserved.
    assert float(matrix.loc["Lamina I Combined", "C1"]) == 10.0


def test_plot_cord_structure_heatmap_writes_png(tmp_path: Path) -> None:
    sub = filter_cord_stats(_cord_stats_rows(), channel=1, metric="median_intensity", rollup_level="structure")
    out = tmp_path / "structure.png"
    plot_cord_structure_heatmap(sub, segment_order=["C1", "C2"], output_path=out)
    assert out.is_file()
    assert out.with_suffix(".csv").is_file()


def test_plot_cord_division_profile_writes_png(tmp_path: Path) -> None:
    sub = filter_cord_stats(_cord_stats_rows(), channel=1, metric="median_intensity", rollup_level="division")
    out = tmp_path / "division.png"
    plot_cord_division_profile(sub, _segments(), output_path=out)
    assert out.is_file()


def test_plot_cord_segment_bars_writes_png(tmp_path: Path) -> None:
    sub = filter_cord_stats(_cord_stats_rows(), channel=1, metric="cell_count", rollup_level="region")
    out = tmp_path / "segments.png"
    plot_cord_segment_bars(sub, segment_order=["C1", "C2"], output_path=out)
    assert out.is_file()


def test_segment_grouped_totals_table_shape() -> None:
    rows = []
    for label, c1, c2 in [
        ("imaris_a", 3.0, 1.0),
        ("imaris_b", 2.0, 4.0),
    ]:
        for seg, val in [("C1", c1), ("C2", c2)]:
            rows.append(
                {
                    "sample": "op87",
                    "channel": label,
                    "atlas": "fiederling",
                    "parcellation_index": 7,
                    "acronym": "5Sp",
                    "name": "Lamina 5",
                    "structure": "DH",
                    "division": "GM",
                    "segment": seg,
                    "rollup_level": "region",
                    "hemisphere": "whole",
                    "metric": "cell_count",
                    "value": val,
                }
            )
    df = pd.DataFrame(rows)
    sub = filter_cord_stats_multi(df, channels=["imaris_a", "imaris_b"], metric="cell_count")
    pivot = segment_grouped_totals_table(sub, channels=["imaris_a", "imaris_b"], segment_order=["C1", "C2"])
    assert list(pivot.columns) == ["segment", "imaris_a", "imaris_b"]
    assert pivot.loc[pivot.segment == "C1", "imaris_a"].iloc[0] == 3.0


def test_plot_cord_segment_grouped_bars_writes_png(tmp_path: Path) -> None:
    rows = []
    for label, c1, c2 in [
        ("imaris_a", 3.0, 1.0),
        ("imaris_b", 2.0, 4.0),
    ]:
        for seg, val in [("C1", c1), ("C2", c2)]:
            rows.append(
                {
                    "sample": "op87",
                    "channel": label,
                    "atlas": "fiederling",
                    "parcellation_index": 7,
                    "acronym": "5Sp",
                    "name": "Lamina 5",
                    "structure": "DH",
                    "division": "GM",
                    "segment": seg,
                    "rollup_level": "region",
                    "hemisphere": "whole",
                    "metric": "cell_count",
                    "value": val,
                }
            )
    sub = filter_cord_stats_multi(pd.DataFrame(rows), channels=["imaris_a", "imaris_b"])
    out = tmp_path / "grouped_segments.png"
    plot_cord_segment_grouped_bars(sub, channels=["imaris_a", "imaris_b"], output_path=out)
    assert out.is_file()
    assert out.with_suffix(".csv").is_file()


def test_top_regions_table_orders_by_value() -> None:
    sub = filter_cord_stats(_cord_stats_rows(), channel=1, metric="cell_count", rollup_level="region")
    top = top_regions_table(sub, top_n=2)
    assert len(top) == 2
    assert top["value"].max() == sub["value"].max()


def test_plot_cord_top_regions_writes_png(tmp_path: Path) -> None:
    sub = filter_cord_stats(_cord_stats_rows(), channel=1, metric="cell_count", rollup_level="region")
    out = tmp_path / "top_regions.png"
    plot_cord_top_regions(sub, top_n=2, output_path=out)
    assert out.is_file()


def test_plot_cord_structure_panel_writes_png(tmp_path: Path) -> None:
    df = _cord_stats_rows()
    out = tmp_path / "panel.png"
    plot_cord_structure_panel(
        df,
        channels=[1],
        metric="median_intensity",
        rollup_level="structure",
        segment_order=["C1", "C2"],
        output_path=out,
    )
    assert out.is_file()


def test_plot_cord_coloc_overlap_writes_png(tmp_path: Path) -> None:
    summary = pd.DataFrame(
        [
            {
                "source": "a",
                "target": "b",
                "n_source": 10,
                "n_target": 8,
                "n_overlap_source_to_target": 4,
                "n_overlap_target_to_source": 4,
                "frac_of_source": 0.4,
                "frac_of_target": 0.5,
                "tolerance_voxels": 2.0,
                "comparison": "pairwise",
            }
        ]
    )
    out = tmp_path / "overlap.png"
    plot_cord_coloc_overlap(summary, output_path=out)
    assert out.is_file()


def test_load_cord_stats_csv_requires_segment_column(tmp_path: Path) -> None:
    path = tmp_path / "bad.csv"
    pd.DataFrame({"metric": ["cell_count"], "value": [1], "channel": [1]}).to_csv(path, index=False)
    with pytest.raises(ValueError, match="segment"):
        load_cord_stats_csv(path)
