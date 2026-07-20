"""Tests for Fiederling spinal cord cell counting."""

from __future__ import annotations

import numpy as np
import pandas as pd

from lightsuite.analysis.cord_counts import (
    assign_segment_names,
    count_points_in_cord_regions,
)
from lightsuite.analysis.cord_ontology import CordRegionTable


def _segments() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "Segment": ["C1", "C2"],
            "Start": [1, 3],
            "End": [2, 4],
        }
    )


def _annotation() -> np.ndarray:
    av = np.zeros((4, 4, 4), dtype=np.int32)
    av[:, :, 0] = 7
    av[:, :, 1] = 7
    av[:, :, 2] = 9
    av[:, :, 3] = 9
    return av


def _region_table() -> CordRegionTable:
    df = pd.DataFrame(
        {
            "parcellation_index": [7, 9],
            "acronym": ["a7", "a9"],
            "name": ["Region7", "Region9"],
            "structure": ["DH", "VH"],
            "division": ["GM", "GM"],
        }
    )
    return CordRegionTable(atlas="fiederling", df=df, source_csv=__import__("pathlib").Path("x.csv"))


def test_assign_segment_names() -> None:
    names = assign_segment_names(np.array([1, 2, 3, 4]), _segments())
    assert names.tolist() == ["C1", "C1", "C2", "C2"]


def test_count_points_in_cord_regions() -> None:
    points = np.array(
        [
            [1.0, 1.0, 1.0],
            [2.0, 2.0, 2.0],
            [3.0, 3.0, 3.0],
            [4.0, 4.0, 4.0],
        ]
    )
    tidy = count_points_in_cord_regions(
        points,
        _annotation(),
        _segments(),
        sample="op93",
        channel="imaris_test",
        region_table=_region_table(),
    )

    counts = tidy[tidy["metric"] == "cell_count"]
    r7_c1 = counts[(counts["parcellation_index"] == 7) & (counts["segment"] == "C1")]
    r9_c2 = counts[(counts["parcellation_index"] == 9) & (counts["segment"] == "C2")]
    assert r7_c1["value"].iloc[0] == 2.0
    assert r9_c2["value"].iloc[0] == 2.0
    assert set(counts["hemisphere"]) == {"whole"}

    density = tidy[
        (tidy["metric"] == "cell_density")
        & (tidy["parcellation_index"] == 7)
        & (tidy["segment"] == "C1")
    ]
    assert density["value"].iloc[0] > 0


def test_count_points_empty_returns_empty_frame() -> None:
    tidy = count_points_in_cord_regions(
        np.zeros((0, 3)),
        _annotation(),
        _segments(),
        sample="op93",
        channel="empty",
    )
    assert len(tidy) == 0
