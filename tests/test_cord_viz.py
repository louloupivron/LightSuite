"""Tests for spinal cord region-stats visualization."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import pandas as pd
import pytest

from lightsuite.analysis.viz.cord_io import (
    df_subregion_table,
    filter_cord_stats,
    filter_cord_stats_multi,
    laminae_level_table,
    laminae_pct_gm_table,
    load_cord_stats_csv,
    segment_grouped_totals_table,
    structure_heatmap_matrix,
    top_regions_table,
)
from lightsuite.analysis.viz.cord_plots import (
    plot_cord_coloc_overlap,
    plot_cord_df_subregion_heatmap,
    plot_cord_division_profile,
    plot_cord_laminae_level_bars,
    plot_cord_laminae_pct_gm_bars,
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


def _laminae_stats_rows() -> pd.DataFrame:
    rows = []
    laminae = [
        ("Lamina_I", 201, 8.0, 2.0, 1.0),
        ("Lamina_V", 205, 10.0, 3.0, 1.2),
        ("Lamina_IX", 209, 12.0, 60.0, 1.5),
    ]
    segments = ["C4", "C5", "T10", "L5"]
    for seg in segments:
        for acronym, pidx, intensity, cells, volume in laminae:
            base = {
                "sample": "op87",
                "atlas": "fiederling",
                "parcellation_index": pidx,
                "acronym": acronym,
                "name": acronym.replace("_", " "),
                "structure": "DH" if pidx < 209 else "VH",
                "division": "GM",
                "segment": seg,
                "rollup_level": "structure",
                "hemisphere": "whole",
            }
            rows.append({**base, "channel": 1, "metric": "median_intensity", "value": intensity})
            rows.append({**base, "channel": 1, "metric": "volume_mm3", "value": volume})
            rows.append({**base, "channel": "imaris_cells", "metric": "cell_count", "value": cells})
    return pd.DataFrame(rows)


def _df_subregion_rows() -> pd.DataFrame:
    rows = []
    regions = [
        ("dcs", 62, 20.0, "region"),
        ("cu", 61, 8.0, "region"),
        ("gr", 63, 2.0, "region"),
        ("psdc", 64, 6.0, "region"),
        ("df", 74, 15.0, "structure"),
    ]
    for seg in ["C4", "C5", "T10"]:
        for acronym, pidx, val, level in regions:
            rows.append(
                {
                    "sample": "op87",
                    "channel": 1,
                    "atlas": "fiederling",
                    "parcellation_index": pidx,
                    "acronym": acronym,
                    "name": acronym,
                    "structure": "WM",
                    "division": "WM",
                    "segment": seg,
                    "rollup_level": level,
                    "hemisphere": "whole",
                    "metric": "median_intensity",
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


def test_division_profile_crops_empty_segments() -> None:
    from lightsuite.analysis.viz.cord_io import division_profile_table

    rows = _cord_stats_rows()
    # Add empty leading/trailing segments with zero intensity.
    extras = []
    for seg, start, end, val in [("T1", 50, 60, 0.0), ("S1", 70, 80, 0.0)]:
        for pidx, acr, name in [(71, "GM", "Gray Matter"), (130, "WM", "White matter")]:
            extras.append(
                {
                    "sample": "op87",
                    "channel": 1,
                    "atlas": "fiederling",
                    "parcellation_index": pidx,
                    "acronym": acr,
                    "name": name,
                    "structure": "SC",
                    "division": acr,
                    "segment": seg,
                    "rollup_level": "division",
                    "hemisphere": "whole",
                    "metric": "median_intensity",
                    "value": val,
                }
            )
    df = pd.concat([rows, pd.DataFrame(extras)], ignore_index=True)
    sub = filter_cord_stats(df, channel=1, metric="median_intensity", rollup_level="division")
    segments = pd.DataFrame(
        {
            "Segment": ["T1", "C1", "C2", "S1"],
            "Start": [50, 1, 15, 70],
            "End": [60, 14, 42, 80],
        }
    )
    profile = division_profile_table(sub, segments, crop_empty_segments=True, drop_nonpositive=True)
    assert set(profile["segment"].astype(str)) == {"C1", "C2"}
    assert profile["value"].isna().sum() == 0


def test_plot_cord_segment_bars_writes_png(tmp_path: Path) -> None:
    from lightsuite.analysis.viz.cord_io import segment_level_class

    assert segment_level_class("L5") == "L"
    assert segment_level_class("Co2") == "Co"
    sub = filter_cord_stats(_cord_stats_rows(), channel=1, metric="cell_count", rollup_level="region")
    out = tmp_path / "segments.png"
    plot_cord_segment_bars(sub, segment_order=["C1", "C2"], output_path=out)
    assert out.is_file()
    csv = pd.read_csv(out.with_suffix(".csv"))
    assert "level" in csv.columns
    assert "pct" in csv.columns


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
    assert top["plot_label"].str.contains("—").all()


def test_top_regions_table_omits_segment_suffix_when_filtered() -> None:
    sub = filter_cord_stats(_cord_stats_rows(), channel=1, metric="cell_count", rollup_level="region")
    top = top_regions_table(sub, top_n=2, segment="C1")
    assert not top["plot_label"].str.contains("@").any()
    assert top["plot_label"].str.startswith("5Sp").any() or top["plot_label"].str.contains("5Sp").any()


def test_plot_cord_top_regions_writes_png(tmp_path: Path) -> None:
    sub = filter_cord_stats(_cord_stats_rows(), channel=1, metric="cell_count", rollup_level="region")
    out = tmp_path / "top_regions.png"
    plot_cord_top_regions(sub, top_n=2, output_path=out)
    assert out.is_file()


def test_plot_cord_structure_panel_writes_png(tmp_path: Path) -> None:
    from lightsuite.analysis.viz.cord_io import align_structure_heatmap_matrices

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

    m1 = pd.DataFrame({"C1": [1.0], "C2": [0.0]}, index=["A"])
    m2 = pd.DataFrame({"C1": [0.0], "C2": [2.0]}, index=["B"])
    aligned, rows, cols = align_structure_heatmap_matrices(
        {"a": m1, "b": m2},
        row_order=["A", "B"],
        drop_empty_rows=True,
    )
    assert rows == ["A", "B"]
    assert cols == ["C1", "C2"]
    assert float(aligned["a"].loc["B"].fillna(0).sum()) == 0.0


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


def test_laminae_pct_gm_table_normalizes_to_100() -> None:
    table = laminae_pct_gm_table(
        _laminae_stats_rows(),
        intensity_channel=1,
        cell_channel="imaris_cells",
    )
    assert pytest.approx(table["intensity_pct_gm"].sum(), rel=1e-6) == 100.0
    assert pytest.approx(table["cell_pct_gm"].sum(), rel=1e-6) == 100.0
    assert table.loc[table["acronym"] == "Lamina_IX", "cell_pct_gm"].iloc[0] > 50.0


def test_laminae_level_table_groups_by_cord_level() -> None:
    table = laminae_level_table(
        _laminae_stats_rows(),
        channel=1,
        metric="median_intensity",
        levels=("C", "T", "L"),
    )
    assert list(table.columns) == ["lamina", "C", "T", "L"]
    assert table.loc[table["lamina"] == "IX", "C"].iloc[0] == 12.0


def test_df_subregion_table_includes_parent_df() -> None:
    table = df_subregion_table(_df_subregion_rows(), channel=1, metric="median_intensity")
    assert set(table["acronym"]) == {"dcs", "cu", "gr", "psdc", "df"}


def test_plot_cord_laminae_pct_gm_bars_writes_png(tmp_path: Path) -> None:
    out = tmp_path / "laminae_pct_gm.png"
    plot_cord_laminae_pct_gm_bars(
        _laminae_stats_rows(),
        intensity_channel=1,
        cell_channel="imaris_cells",
        output_path=out,
    )
    assert out.is_file()
    assert out.with_suffix(".csv").is_file()


def test_plot_cord_laminae_level_bars_writes_png(tmp_path: Path) -> None:
    out = tmp_path / "laminae_level.png"
    plot_cord_laminae_level_bars(
        _laminae_stats_rows(),
        channel=1,
        metric="median_intensity",
        output_path=out,
    )
    assert out.is_file()
    assert out.with_suffix(".csv").is_file()


def test_plot_cord_df_subregion_heatmap_writes_png(tmp_path: Path) -> None:
    out = tmp_path / "df_subregion.png"
    plot_cord_df_subregion_heatmap(
        _df_subregion_rows(),
        channel=1,
        metric="median_intensity",
        segment_order=["C4", "C5", "T10"],
        output_path=out,
    )
    assert out.is_file()
    assert out.with_suffix(".csv").is_file()


def _tiny_annotation_volume() -> np.ndarray:
    import numpy as np

    volume = np.zeros((32, 32, 3), dtype=np.uint16)
    volume[8:24, 8:24, 0] = 7
    volume[8:24, 8:24, 1] = 8
    volume[8:24, 8:24, 2] = 9
    return volume


def _segment_anatomy_stats_rows() -> pd.DataFrame:
    rows = []
    for segment, pidx, acr, val in [
        ("C1", 7, "5Sp", 10.0),
        ("C2", 8, "5SpL", 20.0),
        ("C3", 9, "5SpM", 30.0),
    ]:
        rows.append(
            {
                "sample": "op87",
                "channel": 1,
                "atlas": "fiederling",
                "parcellation_index": pidx,
                "acronym": acr,
                "name": acr,
                "structure": "DH",
                "division": "GM",
                "segment": segment,
                "rollup_level": "region",
                "hemisphere": "whole",
                "metric": "median_intensity",
                "value": val,
            }
        )
    return pd.DataFrame(rows)


def _tiny_hemisphere_volume() -> np.ndarray:
    import numpy as np

    volume = np.zeros((32, 32, 3), dtype=np.uint8)
    volume[:, :16, :] = 0
    volume[:, 16:, :] = 255
    return volume


def _bilateral_segment_stats_rows() -> pd.DataFrame:
    rows = []
    for segment, pidx, acr, left_val, right_val in [
        ("C1", 7, "5Sp", 10.0, 30.0),
        ("C2", 8, "5SpL", 20.0, 40.0),
    ]:
        for hemisphere, value in (("left", left_val), ("right", right_val)):
            rows.append(
                {
                    "sample": "op87",
                    "channel": 1,
                    "atlas": "fiederling",
                    "parcellation_index": pidx,
                    "acronym": acr,
                    "name": acr,
                    "structure": "DH",
                    "division": "GM",
                    "segment": segment,
                    "rollup_level": "region",
                    "hemisphere": hemisphere,
                    "metric": "median_intensity",
                    "value": value,
                }
            )
    return pd.DataFrame(rows)


def test_plot_cord_segment_anatomy_slice_writes_png(tmp_path: Path) -> None:
    import tifffile

    from lightsuite.analysis.viz.cord_segment_anatomy import plot_cord_segment_anatomy_slice

    ann_path = tmp_path / "annotation_registered.tiff"
    tifffile.imwrite(ann_path, _tiny_annotation_volume())
    segments_csv = tmp_path / "Segments.csv"
    segments_csv.write_text("Segment,Ref_Section,Start,End\nC1,1,0,0\nC2,2,1,1\nC3,3,2,2\n")

    out = tmp_path / "segment_slice.png"
    plot_cord_segment_anatomy_slice(
        _segment_anatomy_stats_rows(),
        segments=["C1", "C2", "C3"],
        annotation_path=ann_path,
        segments_csv=segments_csv,
        channel=1,
        metric="median_intensity",
        output_path=out,
        draw_outlines=False,
    )
    assert out.is_file()


def test_plot_cord_segment_anatomy_hemisphere_panel_writes_png(tmp_path: Path) -> None:
    import tifffile

    from lightsuite.analysis.viz.cord_segment_anatomy import (
        plot_cord_segment_anatomy_slice_hemisphere_panel,
    )

    ann_path = tmp_path / "annotation_registered.tiff"
    hem_path = tmp_path / "hemisphere_registered.tiff"
    tifffile.imwrite(ann_path, _tiny_annotation_volume())
    tifffile.imwrite(hem_path, _tiny_hemisphere_volume())
    segments_csv = tmp_path / "Segments.csv"
    segments_csv.write_text("Segment,Ref_Section,Start,End\nC1,1,0,0\nC2,2,1,1\n")

    out = tmp_path / "segment_lr_panel.png"
    plot_cord_segment_anatomy_slice_hemisphere_panel(
        _bilateral_segment_stats_rows(),
        segments=["C1", "C2"],
        annotation_path=ann_path,
        segments_csv=segments_csv,
        hemisphere_volume_path=hem_path,
        channel=1,
        metric="median_intensity",
        output_path=out,
        draw_outlines=False,
    )
    assert out.is_file()


def test_plot_cord_segment_anatomy_hemisphere_composite_writes_png(tmp_path: Path) -> None:
    import tifffile

    from lightsuite.analysis.viz.cord_segment_anatomy import (
        plot_cord_segment_anatomy_slice_hemisphere_composite,
    )

    ann_path = tmp_path / "annotation_registered.tiff"
    hem_path = tmp_path / "hemisphere_registered.tiff"
    tifffile.imwrite(ann_path, _tiny_annotation_volume())
    tifffile.imwrite(hem_path, _tiny_hemisphere_volume())
    segments_csv = tmp_path / "Segments.csv"
    segments_csv.write_text("Segment,Ref_Section,Start,End\nC1,1,0,0\nC2,2,1,1\n")

    out = tmp_path / "segment_lr_composite.png"
    plot_cord_segment_anatomy_slice_hemisphere_composite(
        _bilateral_segment_stats_rows(),
        segments=["C1", "C2"],
        annotation_path=ann_path,
        segments_csv=segments_csv,
        hemisphere_volume_path=hem_path,
        channel=1,
        metric="median_intensity",
        output_path=out,
        draw_outlines=False,
    )
    assert out.is_file()


@pytest.mark.skipif(
    __import__("importlib").util.find_spec("brainglobe_heatmap") is None,
    reason="brainglobe-heatmap not installed",
)
def test_plot_cord_segment_anatomy_bgh_writes_png(tmp_path: Path) -> None:
    from lightsuite.analysis.viz.cord_segment_anatomy import plot_cord_segment_anatomy_bgh

    segments_csv = tmp_path / "Segments.csv"
    segments_csv.write_text("Segment,Ref_Section,Start,End\nC1,1,0,0\n")

    out = tmp_path / "segment_bgh.png"
    plot_cord_segment_anatomy_bgh(
        _segment_anatomy_stats_rows(),
        segments=["C1"],
        segments_csv=segments_csv,
        channel=1,
        metric="median_intensity",
        output_path=out,
    )
    assert out.is_file()
