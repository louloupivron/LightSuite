"""Tests for spinal cord paper-style heatmaps."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from lightsuite.analysis.cord_counts import CORD_TIDY_COLUMNS
from lightsuite.analysis.cord_heatmap import (
    NO_DATA_COLOR,
    PAPER_STRUCTURE_ACRONYMS,
    build_cord_heatmap_figure,
    cord_stats_to_matrix,
    default_paper_segments,
    discover_cord_plot_options,
    parse_segment_range,
    plot_cord_heatmap,
    resolve_hemisphere_options,
    structure_acronym_to_label,
)


def _segments_csv(atlas_dir: Path) -> None:
    atlas_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(
        {
            "Segment": ["C1", "C2", "T1", "Co2", "Co3"],
            "Ref_Section": [1, 2, 3, 4, 5],
            "Start": [1, 2, 3, 4, 5],
            "End": [1, 2, 3, 4, 5],
        }
    ).to_csv(atlas_dir / "Segments.csv", index=False)


def _structure_stats_frame() -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for segment, lamina_idx in [("C1", 0), ("C2", 1), ("T1", 2)]:
        acronym = PAPER_STRUCTURE_ACRONYMS[lamina_idx]
        rows.append(
            {
                "sample": "s1",
                "channel": 1,
                "atlas": "fiederling",
                "parcellation_index": 201 + lamina_idx,
                "acronym": acronym,
                "name": f"Region {acronym}",
                "structure": "DH",
                "division": "GM",
                "segment": segment,
                "rollup_level": "structure",
                "hemisphere": "whole",
                "metric": "median_intensity",
                "value": float((lamina_idx + 1) * 10),
            }
        )
        rows.append(
            {
                "sample": "s1",
                "channel": "imaris_TA",
                "atlas": "fiederling",
                "parcellation_index": 220,
                "acronym": "df",
                "name": "Dorsal funiculus",
                "structure": "WM",
                "division": "WM",
                "segment": segment,
                "rollup_level": "structure",
                "hemisphere": "whole",
                "metric": "cell_count",
                "value": float(segment == "C1"),
            }
        )
    return pd.DataFrame(rows).reindex(columns=CORD_TIDY_COLUMNS)


def test_default_paper_segments_through_co2(tmp_path: Path) -> None:
    atlas_dir = tmp_path / "atlas"
    _segments_csv(atlas_dir)
    assert default_paper_segments(atlas_dir) == ["C1", "C2", "T1", "Co2"]


def test_parse_segment_range_colon_syntax(tmp_path: Path) -> None:
    atlas_dir = tmp_path / "atlas"
    _segments_csv(atlas_dir)
    assert parse_segment_range("C2:Co2", atlas_dir=atlas_dir) == ["C2", "T1", "Co2"]


def test_structure_acronym_labels() -> None:
    assert structure_acronym_to_label("Lamina_III") == "III"
    assert structure_acronym_to_label("df") == "df"


def test_cord_stats_to_matrix_orders_rows_and_columns() -> None:
    df = _structure_stats_frame()
    segments = ["C1", "C2", "T1"]
    matrix = cord_stats_to_matrix(
        df,
        metric="median_intensity",
        channel=1,
        segments=segments,
    )
    assert list(matrix.index) == segments
    assert list(matrix.columns[:3]) == ["Lamina_I", "Lamina_II", "Lamina_III"]
    assert matrix.loc["C1", "Lamina_I"] == 10.0
    assert matrix.loc["C2", "Lamina_II"] == 20.0


def test_cord_stats_to_matrix_cell_counts_by_label() -> None:
    df = _structure_stats_frame()
    matrix = cord_stats_to_matrix(
        df,
        metric="cell_count",
        channel="imaris_TA",
        segments=["C1", "C2", "T1"],
    )
    assert matrix.loc["C1", "df"] == 1.0
    assert matrix.loc["C2", "df"] == 0.0


def test_cord_stats_to_matrix_requires_channel_when_ambiguous() -> None:
    df = _structure_stats_frame()
    extra = df.iloc[[0]].copy()
    extra["channel"] = 2
    df = pd.concat([df, extra], ignore_index=True)
    with pytest.raises(ValueError, match="Multiple channels"):
        cord_stats_to_matrix(
            df,
            metric="median_intensity",
            channel=None,
            segments=["C1"],
        )


def test_plot_cord_heatmap_writes_png(tmp_path: Path) -> None:
    df = _structure_stats_frame()
    matrix = cord_stats_to_matrix(
        df,
        metric="median_intensity",
        channel=1,
        segments=["C1", "C2", "T1"],
    )
    output = tmp_path / "heatmap.png"
    plot_cord_heatmap(matrix, output, title="test")
    assert output.is_file()
    assert output.stat().st_size > 0


def test_discover_cord_plot_options_lists_metrics_and_channels() -> None:
    df = _structure_stats_frame()
    options = discover_cord_plot_options(df)
    assert "median_intensity" in options["metrics"]
    assert "cell_count" in options["metrics"]
    assert options["channels_by_metric"]["median_intensity"] == ["1"]
    assert options["channels_by_metric"]["cell_count"] == ["imaris_TA"]


def test_build_cord_heatmap_figure_returns_axes() -> None:
    df = _structure_stats_frame()
    matrix = cord_stats_to_matrix(
        df,
        metric="median_intensity",
        channel=1,
        segments=["C1", "C2"],
    )
    fig = build_cord_heatmap_figure(matrix, title="test")
    assert len(fig.axes) == 2


def test_build_cord_heatmap_figure_marks_missing_data_grey() -> None:
    import numpy as np
    from matplotlib import colors as mcolors

    matrix = pd.DataFrame(
        [[10.0, np.nan], [20.0, 30.0]],
        index=["C1", "C2"],
        columns=["Lamina_I", "Lamina_II"],
    )
    fig = build_cord_heatmap_figure(matrix, title="test")
    image = fig.axes[0].images[0]
    assert image.get_cmap().get_bad()[:3] == pytest.approx(mcolors.to_rgb(NO_DATA_COLOR))


def _split_hemisphere_stats_frame() -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for hemisphere, intensity, volume in (
        ("left", 10.0, 1.0),
        ("right", 30.0, 3.0),
    ):
        rows.extend(
            [
                {
                    "sample": "s1",
                    "channel": 1,
                    "atlas": "fiederling",
                    "parcellation_index": 201,
                    "acronym": "Lamina_I",
                    "name": "Lamina I",
                    "structure": "DH",
                    "division": "GM",
                    "segment": "C1",
                    "rollup_level": "structure",
                    "hemisphere": hemisphere,
                    "metric": "median_intensity",
                    "value": intensity,
                },
                {
                    "sample": "s1",
                    "channel": 1,
                    "atlas": "fiederling",
                    "parcellation_index": 201,
                    "acronym": "Lamina_I",
                    "name": "Lamina I",
                    "structure": "DH",
                    "division": "GM",
                    "segment": "C1",
                    "rollup_level": "structure",
                    "hemisphere": hemisphere,
                    "metric": "volume_mm3",
                    "value": volume,
                },
            ]
        )
    return pd.DataFrame(rows).reindex(columns=CORD_TIDY_COLUMNS)


def test_resolve_hemisphere_options_adds_whole_for_split_stats() -> None:
    df = _split_hemisphere_stats_frame()
    assert resolve_hemisphere_options(df) == ["whole", "left", "right"]


def test_cord_stats_to_matrix_whole_merges_split_hemispheres() -> None:
    df = _split_hemisphere_stats_frame()
    matrix = cord_stats_to_matrix(
        df,
        metric="median_intensity",
        channel=1,
        hemisphere="whole",
        segments=["C1"],
    )
    assert matrix.loc["C1", "Lamina_I"] == pytest.approx(25.0)


def test_discover_cord_plot_options_includes_whole_hemisphere() -> None:
    df = _split_hemisphere_stats_frame()
    options = discover_cord_plot_options(df)
    assert options["hemispheres"][0] == "whole"
