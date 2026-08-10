"""Tests for spinal cord division/structure rollups."""

from __future__ import annotations

import pandas as pd

from lightsuite.analysis.cord_counts import CORD_TIDY_COLUMNS
from lightsuite.analysis.cord_rollup import (
    apply_cord_rollups,
    get_descendants,
    resolve_horn_acronym,
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


def _mini_regions_with_vh() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "id": [1, 18, 71, 90, 110, 130],
            "name": ["Lamina1", "Lamina7", "Gray Matter", "Dorsal horn", "Ventral horn", "White matter"],
            "acronym": ["1Sp", "7Sp", "GM", "DH", "VH", "WM"],
            "parent_ID": [90, 110, 250, 71, 71, 250],
            "parent_acronym": ["DH", "VH", "SC", "GM", "GM", "SC"],
            "children_IDs": ["", "", "1,18", "1", "18", ""],
        }
    )


def test_rollup_horn_splits_dorsal_ventral() -> None:
    rows = [
        dict(
            sample="s1", channel=1, atlas="fiederling",
            parcellation_index=1, acronym="1Sp", name="Lamina1",
            structure="DH", division="GM", segment="C1",
            rollup_level="region", hemisphere="whole",
            metric="cell_count", value=5.0,
        ),
        dict(
            sample="s1", channel=1, atlas="fiederling",
            parcellation_index=1, acronym="1Sp", name="Lamina1",
            structure="DH", division="GM", segment="C1",
            rollup_level="region", hemisphere="whole",
            metric="volume_mm3", value=2.0,
        ),
        dict(
            sample="s1", channel=1, atlas="fiederling",
            parcellation_index=18, acronym="7Sp", name="Lamina7",
            structure="VH", division="GM", segment="C1",
            rollup_level="region", hemisphere="whole",
            metric="cell_count", value=20.0,
        ),
        dict(
            sample="s1", channel=1, atlas="fiederling",
            parcellation_index=18, acronym="7Sp", name="Lamina7",
            structure="VH", division="GM", segment="C1",
            rollup_level="region", hemisphere="whole",
            metric="volume_mm3", value=4.0,
        ),
    ]
    df = pd.DataFrame(rows).reindex(columns=CORD_TIDY_COLUMNS)
    rolled = rollup_cord_tidy(df, _mini_regions_with_vh(), "horn")
    dh = rolled[(rolled["rollup_level"] == "horn") & (rolled["acronym"] == "DH") & (rolled["metric"] == "cell_count")]
    vh = rolled[(rolled["rollup_level"] == "horn") & (rolled["acronym"] == "VH") & (rolled["metric"] == "cell_count")]
    assert len(dh) == 1
    assert dh["value"].iloc[0] == 5.0
    assert len(vh) == 1
    assert vh["value"].iloc[0] == 20.0


def test_resolve_horn_targets() -> None:
    regions = _mini_regions_with_vh()
    targets = resolve_rollup_targets("horn", regions)
    acronyms = {acr for _id, acr, _name in targets}
    assert "DH" in acronyms
    assert "VH" in acronyms


def test_resolve_horn_acronym_walks_parent_chain() -> None:
    regions = _mini_regions_with_vh()
    assert resolve_horn_acronym(1, regions) == "DH"
    assert resolve_horn_acronym(18, regions) == "VH"


def _mini_regions_with_overlapping_horn_descendants() -> pd.DataFrame:
    """DH combined node lists a VH child in children_IDs (Fiederling-style overlap)."""
    return pd.DataFrame(
        {
            "id": [1, 18, 71, 90, 110, 207],
            "name": ["Lamina1", "Lamina7", "Gray Matter", "Dorsal horn", "Ventral horn", "Lamina VII Combined"],
            "acronym": ["1Sp", "7Sp", "GM", "DH", "VH", "Lamina_VII"],
            "parent_ID": [90, 110, 250, 71, 71, 90],
            "parent_acronym": ["DH", "VH", "SC", "GM", "GM", "DH"],
            "children_IDs": ["", "", "1,18", "1,207", "18", "18"],
        }
    )


def test_rollup_horn_uses_parent_walk_not_descendant_overlap() -> None:
    regions = _mini_regions_with_overlapping_horn_descendants()
    assert 18 in get_descendants(90, regions)
    rows = [
        dict(
            sample="s1",
            channel="spots",
            atlas="fiederling",
            parcellation_index=1,
            acronym="1Sp",
            name="Lamina1",
            structure="DH",
            division="GM",
            segment="L4",
            rollup_level="region",
            hemisphere="right",
            metric="cell_count",
            value=5.0,
        ),
        dict(
            sample="s1",
            channel="spots",
            atlas="fiederling",
            parcellation_index=18,
            acronym="7Sp",
            name="Lamina7",
            structure="VH",
            division="GM",
            segment="L4",
            rollup_level="region",
            hemisphere="right",
            metric="cell_count",
            value=20.0,
        ),
    ]
    df = pd.DataFrame(rows).reindex(columns=CORD_TIDY_COLUMNS)
    rolled = rollup_cord_tidy(df, regions, "horn")
    dh = rolled[(rolled["rollup_level"] == "horn") & (rolled["acronym"] == "DH") & (rolled["metric"] == "cell_count")]
    vh = rolled[(rolled["rollup_level"] == "horn") & (rolled["acronym"] == "VH") & (rolled["metric"] == "cell_count")]
    assert dh["value"].iloc[0] == 5.0
    assert vh["value"].iloc[0] == 20.0


def test_rollup_division_cell_density_from_counts_and_volume() -> None:
    rolled = rollup_cord_tidy(_point_count_frame(), _mini_regions(), "division")
    gm_density = rolled[
        (rolled["rollup_level"] == "division")
        & (rolled["acronym"] == "GM")
        & (rolled["metric"] == "cell_density")
    ]
    assert len(gm_density) == 1
    assert gm_density["value"].iloc[0] == 10.0
