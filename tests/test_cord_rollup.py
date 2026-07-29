"""Tests for spinal cord division/structure rollups."""

from __future__ import annotations

import pandas as pd

from lightsuite.analysis.cord_counts import CORD_TIDY_COLUMNS
from lightsuite.analysis.cord_rollup import (
    apply_cord_rollups,
    get_descendants,
    resolve_rollup_targets,
    rollup_cord_tidy,
)


def _mini_regions() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "id": [1, 2, 71, 90, 130, 201],
            "name": ["Lamina1", "Lamina2", "Gray Matter", "Dorsal horn", "White matter", "Lamina I Combined"],
            "acronym": ["1Sp", "2Ssp", "GM", "DH", "WM", "Lamina_I"],
            "parent_ID": [90, 90, 250, 71, 250, 90],
            "parent_acronym": ["DH", "DH", "SC", "GM", "SC", "DH"],
            "children_IDs": ["", "", "1,2", "", "", "1,2"],
        }
    )


def _finest_frame() -> pd.DataFrame:
    rows = [
        dict(
            sample="s1",
            channel=1,
            atlas="fiederling",
            parcellation_index=1,
            acronym="1Sp",
            name="Lamina1",
            structure="DH",
            division="GM",
            segment="C1",
            rollup_level="region",
            hemisphere="whole",
            metric="median_intensity",
            value=100.0,
        ),
        dict(
            sample="s1",
            channel=1,
            atlas="fiederling",
            parcellation_index=1,
            acronym="1Sp",
            name="Lamina1",
            structure="DH",
            division="GM",
            segment="C1",
            rollup_level="region",
            hemisphere="whole",
            metric="volume_mm3",
            value=10.0,
        ),
        dict(
            sample="s1",
            channel=1,
            atlas="fiederling",
            parcellation_index=2,
            acronym="2Ssp",
            name="Lamina2",
            structure="DH",
            division="GM",
            segment="C1",
            rollup_level="region",
            hemisphere="whole",
            metric="median_intensity",
            value=300.0,
        ),
        dict(
            sample="s1",
            channel=1,
            atlas="fiederling",
            parcellation_index=2,
            acronym="2Ssp",
            name="Lamina2",
            structure="DH",
            division="GM",
            segment="C1",
            rollup_level="region",
            hemisphere="whole",
            metric="volume_mm3",
            value=30.0,
        ),
    ]
    return pd.DataFrame(rows).reindex(columns=CORD_TIDY_COLUMNS)


def test_get_descendants_includes_children_ids() -> None:
    regions = _mini_regions()
    descendants = get_descendants(90, regions)
    assert 1 in descendants
    assert 2 in descendants


def test_resolve_division_targets() -> None:
    targets = resolve_rollup_targets("division", _mini_regions())
    acronyms = {acr for _id, acr, _name in targets}
    assert "GM" in acronyms
    assert "WM" in acronyms


def test_rollup_division_volume_weighted_intensity() -> None:
    rolled = rollup_cord_tidy(_finest_frame(), _mini_regions(), "division")
    gm = rolled[
        (rolled["rollup_level"] == "division")
        & (rolled["acronym"] == "GM")
        & (rolled["metric"] == "median_intensity")
    ]
    assert len(gm) == 1
    # weighted: (100*10 + 300*30) / 40 = 250
    assert gm["value"].iloc[0] == 250.0


def test_apply_cord_rollups_appends_levels() -> None:
    combined = apply_cord_rollups(_finest_frame(), _mini_regions(), ["division"])
    assert set(combined["rollup_level"]) == {"region", "division"}
    assert len(combined) > len(_finest_frame())


def _point_count_frame() -> pd.DataFrame:
    rows = [
        dict(
            sample="s1",
            channel="spots_a",
            atlas="fiederling",
            parcellation_index=1,
            acronym="1Sp",
            name="Lamina1",
            structure="DH",
            division="GM",
            segment="C1",
            rollup_level="region",
            hemisphere="whole",
            metric="cell_count",
            value=10.0,
        ),
        dict(
            sample="s1",
            channel="spots_a",
            atlas="fiederling",
            parcellation_index=1,
            acronym="1Sp",
            name="Lamina1",
            structure="DH",
            division="GM",
            segment="C1",
            rollup_level="region",
            hemisphere="whole",
            metric="volume_mm3",
            value=1.0,
        ),
        dict(
            sample="s1",
            channel="spots_a",
            atlas="fiederling",
            parcellation_index=1,
            acronym="1Sp",
            name="Lamina1",
            structure="DH",
            division="GM",
            segment="C1",
            rollup_level="region",
            hemisphere="whole",
            metric="cell_density",
            value=10.0,
        ),
    ]
    return pd.DataFrame(rows).reindex(columns=CORD_TIDY_COLUMNS)


def test_rollup_division_cell_density_from_counts_and_volume() -> None:
    rolled = rollup_cord_tidy(_point_count_frame(), _mini_regions(), "division")
    gm_density = rolled[
        (rolled["rollup_level"] == "division")
        & (rolled["acronym"] == "GM")
        & (rolled["metric"] == "cell_density")
    ]
    assert len(gm_density) == 1
    assert gm_density["value"].iloc[0] == 10.0
