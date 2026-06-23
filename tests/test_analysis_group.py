"""Tests for cross-subject group analysis."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from lightsuite.analysis.cohort_models import CohortConfig, CohortSampleEntry, GroupAnalysisConfig
from lightsuite.analysis.group import (
    benjamini_hochberg,
    compare_groups_pair,
    filter_cohort_table,
    load_cohort_long,
    run_group_analysis,
    summarize_by_group,
)
from lightsuite.analysis.region_stats import TIDY_COLUMNS


def _tidy_row(
    *,
    sample: str,
    group: str,
    channel: int = 1,
    parcellation_index: int = 1,
    hemisphere: str = "right",
    metric: str = "median_intensity",
    value: float,
    division: str = "Isocortex",
) -> dict:
    return {
        "sample": sample,
        "channel": channel,
        "atlas": "allen",
        "parcellation_index": parcellation_index,
        "acronym": "MOp1",
        "name": "Primary motor area, Layer 1",
        "structure": "MO",
        "division": division,
        "hemisphere": hemisphere,
        "metric": metric,
        "value": value,
        "group": group,
    }


def _write_stats(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(path, index=False)


def test_load_cohort_long_concatenates(tmp_path: Path) -> None:
    p1 = tmp_path / "s1" / "volume_registered" / "region_stats.csv"
    p2 = tmp_path / "s2" / "volume_registered" / "region_stats.csv"
    _write_stats(p1, [_tidy_row(sample="s1", group="control", value=10.0)])
    _write_stats(p2, [_tidy_row(sample="s2", group="treatment", value=20.0)])

    from lightsuite.analysis.group import ResolvedCohortSample

    long = load_cohort_long(
        [
            ResolvedCohortSample("s1", "control", p1),
            ResolvedCohortSample("s2", "treatment", p2),
        ]
    )
    assert len(long) == 2
    assert set(long["group"]) == {"control", "treatment"}


def test_summarize_by_group_mean_and_sem() -> None:
    df = pd.DataFrame(
        [
            _tidy_row(sample="a", group="g1", value=10.0),
            _tidy_row(sample="b", group="g1", value=20.0),
            _tidy_row(sample="c", group="g2", value=30.0),
        ]
    )
    summary = summarize_by_group(df)
    g1 = summary[summary["group"] == "g1"].iloc[0]
    assert g1["n_subjects"] == 2
    assert g1["mean"] == 15.0
    assert g1["sem"] == pytest.approx(g1["std"] / np.sqrt(2))


def test_compare_groups_detects_difference() -> None:
    rows = []
    for sample, val in zip(["a", "b", "c", "d"], [1.0, 2.0, 1.5, 2.5]):
        rows.append(_tidy_row(sample=sample, group="control", value=val))
    for sample, val in zip(["e", "f", "g", "h"], [10.0, 11.0, 10.5, 11.5]):
        rows.append(_tidy_row(sample=sample, group="treatment", value=val))
    df = pd.DataFrame(rows)
    comp = compare_groups_pair(df, "control", "treatment", fdr_alpha=0.05)
    assert len(comp) == 1
    assert comp["p_value"].iloc[0] < 0.05
    assert comp["median_diff"].iloc[0] < 0


def test_division_rollup_summary() -> None:
    df = pd.DataFrame(
        [
            _tidy_row(sample="a", group="g1", parcellation_index=1, value=10.0, division="Iso"),
            _tidy_row(sample="a", group="g1", parcellation_index=2, value=20.0, division="Iso"),
            _tidy_row(sample="b", group="g1", parcellation_index=1, value=12.0, division="Iso"),
            _tidy_row(sample="b", group="g1", parcellation_index=2, value=18.0, division="Iso"),
        ]
    )
    # patch parcellation_index=2 row metadata
    df.loc[1, "acronym"] = "MOp2"
    df.loc[3, "acronym"] = "MOp2"
    summary = summarize_by_group(df, level="division")
    assert len(summary) == 1
    assert summary["division"].iloc[0] == "Iso"
    assert summary["mean"].iloc[0] == 15.0  # mean of per-subject means (15, 15)


def test_benjamini_hochberg_orders_q_values() -> None:
    p = np.array([0.01, 0.04, 0.03, 0.20])
    q, reject = benjamini_hochberg(p, alpha=0.05)
    assert reject[0]
    assert np.all(q[np.isfinite(q)] <= 1.0)


def test_filter_cohort_table_channels() -> None:
    df = pd.DataFrame(
        [
            _tidy_row(sample="a", group="g", channel=1, value=1.0),
            _tidy_row(sample="b", group="g", channel=2, value=2.0),
        ]
    )
    opts = GroupAnalysisConfig(channels=[1])
    filtered = filter_cohort_table(df, opts)
    assert len(filtered) == 1
    assert int(filtered["channel"].iloc[0]) == 1


def test_run_group_analysis_writes_summary(tmp_path: Path) -> None:
    paths = []
    for sid, grp, val in [("a", "control", 5.0), ("b", "control", 6.0), ("c", "treat", 50.0), ("d", "treat", 55.0)]:
        p = tmp_path / sid / "volume_registered" / "region_stats.csv"
        _write_stats(p, [_tidy_row(sample=sid, group=grp, value=val)])
        paths.append(p)

    cfg = CohortConfig(
        name="test",
        output_dir=tmp_path / "out",
        samples=[
            CohortSampleEntry(group="control", id="a", region_stats=paths[0]),
            CohortSampleEntry(group="control", id="b", region_stats=paths[1]),
            CohortSampleEntry(group="treat", id="c", region_stats=paths[2]),
            CohortSampleEntry(group="treat", id="d", region_stats=paths[3]),
        ],
        group_analysis=GroupAnalysisConfig(
            comparisons=[["control", "treat"]],
            rollups=["division"],
        ),
    )
    result = run_group_analysis(cfg)
    assert len(result.summary_by_region) == 2
    assert result.comparisons_by_region is not None
    assert len(result.comparisons_by_region) == 1
    assert result.summary_by_division is not None


def test_cohort_runner_writes_files(tmp_path: Path) -> None:
    from lightsuite.analysis.cohort_runner import run_cohort_group_analysis

    paths = []
    for sid, grp, val in [("a", "control", 5.0), ("b", "control", 6.0), ("c", "treat", 50.0), ("d", "treat", 55.0)]:
        p = tmp_path / sid / "volume_registered" / "region_stats.csv"
        _write_stats(p, [_tidy_row(sample=sid, group=grp, value=val)])
        paths.append(p)

    cfg = CohortConfig(
        name="test",
        output_dir=tmp_path / "out",
        samples=[
            CohortSampleEntry(group="control", id="a", region_stats=paths[0]),
            CohortSampleEntry(group="control", id="b", region_stats=paths[1]),
            CohortSampleEntry(group="treat", id="c", region_stats=paths[2]),
            CohortSampleEntry(group="treat", id="d", region_stats=paths[3]),
        ],
        group_analysis=GroupAnalysisConfig(comparisons=[["control", "treat"]]),
    )
    result = run_cohort_group_analysis(cfg)
    assert result.cohort_long_path is not None
    assert result.summary_region_path is not None
    assert result.comparisons_region_path is not None
    assert len(result.written_paths) >= 3


def test_cohort_sample_requires_source() -> None:
    with pytest.raises(ValueError, match="region_stats or config"):
        CohortSampleEntry(group="control")
