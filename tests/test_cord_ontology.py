"""Tests for Fiederling region table loading."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from lightsuite.analysis.cord_ontology import load_cord_region_table
from lightsuite.config.models import CordAtlasConfig


def test_load_cord_region_table(tmp_path: Path) -> None:
    atlas_dir = tmp_path / "atlas"
    atlas_dir.mkdir()
    pd.DataFrame(
        {
            "id": [1, 90],
            "name": ["Lamina 1", "Dorsal horn"],
            "acronym": ["1Sp", "DH"],
            "parent_ID": [90, 71],
            "parent_acronym": ["DH", "GM"],
        }
    ).to_csv(atlas_dir / "Atlas_Regions.csv", index=False)
    (atlas_dir / "Template.tif").write_bytes(b"")
    (atlas_dir / "Annotation.tif").write_bytes(b"")
    (atlas_dir / "Segments.csv").write_text("Segment,Start,End\nC1,1,2\n", encoding="utf-8")

    table = load_cord_region_table(CordAtlasConfig(atlas_dir=atlas_dir))
    row = table.df.loc[table.df["parcellation_index"] == 1].iloc[0]
    assert row["acronym"] == "1Sp"
    assert row["structure"] == "DH"
    assert row["division"] == "GM"
