"""Tests for the shared analysis region ontology."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from lightsuite.analysis.ontology import (
    build_ccf_to_division,
    build_ccf_to_structure,
    load_allen_region_table,
    load_perens_region_table,
    resolve_ccf_by_hierarchy,
    structures_parent_map,
)


def _write_allen_membership(path: Path) -> Path:
    """Minimal ABC long-form membership: 2 indices × 3 term sets we use."""
    rows = []
    spec = {
        1: {
            "division": ("Isocortex", "Isocortex", 315),
            "structure": ("Primary motor area", "MO", 985),
            "substructure": ("Primary motor area, Layer 1", "MOp1", 320),
        },
        2: {
            "division": ("Thalamus", "TH", 549),
            "structure": ("Ventral group", "VENT", 637),
            "substructure": ("Ventral posterolateral nucleus", "VPL", 718),
        },
    }
    for pidx, levels in spec.items():
        for term_set, (name, acronym, ccf) in levels.items():
            rows.append(
                {
                    "parcellation_index": pidx,
                    "parcellation_term_set_name": term_set,
                    "parcellation_term_name": name,
                    "parcellation_term_acronym": acronym,
                    "parcellation_term_label": f"AllenCCF-Ontology-2017-{ccf}",
                }
            )
    pd.DataFrame(rows).to_csv(path, index=False)
    return path


def test_load_allen_region_table(tmp_path: Path) -> None:
    csv = _write_allen_membership(tmp_path / "membership.csv")
    table = load_allen_region_table(csv)

    assert table.atlas == "allen"
    row = table.df.set_index("parcellation_index").loc[1]
    assert row["acronym"] == "MOp1"
    assert row["name"] == "Primary motor area, Layer 1"
    assert row["structure"] == "MO"
    assert row["division"] == "Isocortex"
    assert row["division_acronym"] == "Isocortex"
    assert int(row["ccf_id"]) == 320


def test_build_ccf_to_division_covers_all_levels(tmp_path: Path) -> None:
    csv = _write_allen_membership(tmp_path / "membership.csv")
    mapping = build_ccf_to_division(csv)
    # leaf, structure and division ccf ids all resolve to the division name
    assert mapping[320] == "Isocortex"  # substructure leaf
    assert mapping[985] == "Isocortex"  # structure
    assert mapping[315] == "Isocortex"  # division node
    assert mapping[718] == "Thalamus"


def test_load_perens_region_table_translates_division(tmp_path: Path) -> None:
    allen_csv = _write_allen_membership(tmp_path / "membership.csv")
    perens_csv = tmp_path / "ARA2.csv"
    pd.DataFrame(
        {
            "id": [320, 718, 99999, 5001],
            "parent_id": [985, 637, 985, 320],
            "name": ["MOp1", "VPL", "Made up", "Child layer"],
            "acronym": ["MOp1", "VPL", "MU", "CL1"],
        }
    ).to_csv(perens_csv, index=False)

    table = load_perens_region_table(perens_csv, allen_membership_csv=allen_csv)
    indexed = table.df.set_index("parcellation_index")
    assert table.atlas == "perens"
    assert indexed.loc[320, "division"] == "Isocortex"
    assert indexed.loc[320, "structure"] == "MO"
    assert indexed.loc[718, "division"] == "Thalamus"
    assert indexed.loc[718, "structure"] == "VENT"
    # parent 985 is a mapped structure-level CCF id
    assert indexed.loc[99999, "division"] == "Isocortex"
    assert indexed.loc[99999, "structure"] == "MO"
    # parent 320 is a mapped leaf id
    assert indexed.loc[5001, "division"] == "Isocortex"
    assert indexed.loc[5001, "structure"] == "MO"


def test_resolve_ccf_by_hierarchy_walks_parents() -> None:
    parent_map = {5002: 5003, 5003: 985, 985: 315, 315: None}
    value_map = {320: "Isocortex", 985: "Isocortex", 315: "Isocortex"}
    assert resolve_ccf_by_hierarchy(5002, value_map, parent_map) == "Isocortex"
    assert resolve_ccf_by_hierarchy(9999, value_map, parent_map) is None


def test_build_ccf_to_structure_maps_leaves_and_parents(tmp_path: Path) -> None:
    csv = _write_allen_membership(tmp_path / "membership.csv")
    mapping = build_ccf_to_structure(csv)
    assert mapping[320] == "MO"
    assert mapping[985] == "MO"
    assert mapping[718] == "VENT"


def test_allen_membership_missing_columns_raises(tmp_path: Path) -> None:
    bad = tmp_path / "bad.csv"
    pd.DataFrame({"parcellation_index": [1]}).to_csv(bad, index=False)
    with pytest.raises(ValueError, match="missing columns"):
        load_allen_region_table(bad)
