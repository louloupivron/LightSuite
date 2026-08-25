"""Tests for top-N region ranking from tidy region stats."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from lightsuite.analysis.region_stats import TIDY_COLUMNS
from lightsuite.analysis.top_regions import (
    maybe_write_top_n_regions_csv,
    resolve_rank_metric,
    top_n_regions_wide,
)


def _row(
    *,
    channel: int | str,
    parcellation_index: int,
    metric: str,
    value: float,
    hemisphere: str = "right",
    acronym: str = "R",
    name: str = "Region",
    division: str = "Cortex",
    segment: str | None = None,
    rollup_level: str | None = None,
) -> dict[str, object]:
    row: dict[str, object] = {
        "sample": "mouse",
        "channel": channel,
        "atlas": "allen",
        "parcellation_index": parcellation_index,
        "acronym": acronym,
        "name": name,
        "structure": "S",
        "division": division,
        "hemisphere": hemisphere,
        "metric": metric,
        "value": value,
    }
    if segment is not None:
        row["segment"] = segment
    if rollup_level is not None:
        row["rollup_level"] = rollup_level
    return row


def _tidy(rows: list[dict[str, object]]) -> pd.DataFrame:
    df = pd.DataFrame(rows)
    extra = [col for col in df.columns if col not in TIDY_COLUMNS]
    return df.reindex(columns=TIDY_COLUMNS + extra)


def test_ranks_by_summed_cell_count_and_keeps_metrics() -> None:
    tidy = _tidy(
        [
            _row(channel=1, parcellation_index=1, metric="cell_count", value=10, hemisphere="right", acronym="A"),
            _row(channel=1, parcellation_index=1, metric="cell_count", value=5, hemisphere="left", acronym="A"),
            _row(channel=1, parcellation_index=1, metric="median_intensity", value=1.0, hemisphere="right", acronym="A"),
            _row(channel=1, parcellation_index=1, metric="median_intensity", value=3.0, hemisphere="left", acronym="A"),
            _row(channel=1, parcellation_index=2, metric="cell_count", value=100, hemisphere="right", acronym="B"),
            _row(channel=1, parcellation_index=2, metric="cell_count", value=1, hemisphere="left", acronym="B"),
            _row(channel=1, parcellation_index=2, metric="median_intensity", value=9.0, hemisphere="right", acronym="B"),
            _row(channel=1, parcellation_index=0, metric="cell_count", value=999, hemisphere="right", acronym="bg"),
        ]
    )
    wide = top_n_regions_wide(tidy, n=1)
    assert len(wide) == 1
    assert wide.loc[0, "parcellation_index"] == 2
    assert wide.loc[0, "rank"] == 1
    assert wide.loc[0, "cell_count_right"] == 100
    assert wide.loc[0, "cell_count_left"] == 1
    assert wide.loc[0, "median_intensity_right"] == 9.0


def test_ranks_each_channel_separately() -> None:
    tidy = _tidy(
        [
            _row(channel="cells", parcellation_index=1, metric="cell_count", value=3, acronym="A"),
            _row(channel="cells", parcellation_index=2, metric="cell_count", value=8, acronym="B"),
            _row(channel=1, parcellation_index=1, metric="median_intensity", value=50, acronym="A"),
            _row(channel=1, parcellation_index=2, metric="median_intensity", value=10, acronym="B"),
        ]
    )
    wide = top_n_regions_wide(tidy, n=1)
    assert set(wide["channel"]) == {1, "cells"}
    by_channel = wide.set_index("channel")
    assert by_channel.loc["cells", "parcellation_index"] == 2
    assert by_channel.loc[1, "parcellation_index"] == 1


def test_mean_intensity_ranking_and_zero_disables_write(tmp_path: Path) -> None:
    tidy = _tidy(
        [
            _row(channel=1, parcellation_index=1, metric="median_intensity", value=2.0, hemisphere="right"),
            _row(channel=1, parcellation_index=1, metric="median_intensity", value=8.0, hemisphere="left"),
            _row(channel=1, parcellation_index=3, metric="median_intensity", value=9.0, hemisphere="right"),
            _row(channel=1, parcellation_index=3, metric="median_intensity", value=9.0, hemisphere="left"),
        ]
    )
    assert resolve_rank_metric(tidy, None) == "median_intensity"
    wide = top_n_regions_wide(tidy, n=1, rank_by="median_intensity")
    assert wide.loc[0, "parcellation_index"] == 3
    assert maybe_write_top_n_regions_csv(tmp_path, tidy, n=0) is None
    written = maybe_write_top_n_regions_csv(tmp_path, tidy, n=2, rank_by="median_intensity")
    assert written == tmp_path / "region_stats_top2.csv"
    assert written.is_file()


def test_cord_keeps_segment_and_drops_rollups() -> None:
    tidy = _tidy(
        [
            _row(
                channel=1,
                parcellation_index=7,
                metric="median_intensity",
                value=12.0,
                hemisphere="whole",
                segment="C1",
                rollup_level="region",
                acronym="L1",
            ),
            _row(
                channel=1,
                parcellation_index=8,
                metric="median_intensity",
                value=99.0,
                hemisphere="whole",
                segment="C1",
                rollup_level="division",
                acronym="GM",
            ),
            _row(
                channel=1,
                parcellation_index=9,
                metric="median_intensity",
                value=4.0,
                hemisphere="whole",
                segment="C2",
                rollup_level="region",
                acronym="L2",
            ),
        ]
    )
    wide = top_n_regions_wide(tidy, n=10)
    assert set(wide["parcellation_index"]) == {7, 9}
    assert "segment" in wide.columns
    assert "median_intensity" in wide.columns
    assert wide.loc[wide["parcellation_index"] == 7, "rank"].iloc[0] == 1
