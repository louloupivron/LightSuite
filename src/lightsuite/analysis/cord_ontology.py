"""Fiederling spinal cord region metadata for tidy cell-count tables."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from lightsuite.atlas.fiederling import resolve_fiederling_paths
from lightsuite.config.models import CordAtlasConfig


@dataclass(frozen=True)
class CordRegionTable:
    """Region metadata keyed on Fiederling ``Atlas_Regions.csv`` ids."""

    atlas: str
    df: pd.DataFrame
    source_csv: Path


def load_cord_region_table(atlas: CordAtlasConfig, *, atlas_id: str = "fiederling") -> CordRegionTable:
    """Load and normalize ``Atlas_Regions.csv`` for joins on ``parcellation_index``."""
    paths = resolve_fiederling_paths(atlas.atlas_dir)
    raw = pd.read_csv(paths.regions_csv)
    if "id" not in raw.columns:
        msg = f"Atlas_Regions.csv must contain an 'id' column: {paths.regions_csv}"
        raise ValueError(msg)

    regions = raw.copy()
    regions["parcellation_index"] = regions["id"].astype("int64")
    regions["acronym"] = regions["acronym"].astype(str) if "acronym" in regions.columns else ""
    regions["name"] = regions["name"].astype(str) if "name" in regions.columns else ""

    parent_acronym = (
        regions["parent_acronym"].astype(str)
        if "parent_acronym" in regions.columns
        else pd.Series("", index=regions.index)
    )
    regions["structure"] = parent_acronym

    id_to_parent = (
        regions.set_index("id")["parent_ID"].to_dict()
        if "parent_ID" in regions.columns
        else {}
    )
    id_to_parent_acronym = (
        regions.set_index("id")["parent_acronym"].to_dict()
        if "parent_acronym" in regions.columns
        else {}
    )

    def _division_for(region_id: int) -> str:
        parent_id = id_to_parent.get(region_id)
        if parent_id is None or pd.isna(parent_id):
            return ""
        parent_id = int(parent_id)
        parent_acr = str(id_to_parent_acronym.get(parent_id, ""))
        if parent_acr in {"GM", "WM", "SC", "CNS"}:
            return parent_acr
        grandparent = id_to_parent.get(parent_id)
        if grandparent is None or pd.isna(grandparent):
            return parent_acr
        return str(id_to_parent_acronym.get(int(grandparent), parent_acr))

    regions["division"] = [_division_for(int(rid)) for rid in regions["parcellation_index"]]

    keep = ["parcellation_index", "acronym", "name", "structure", "division"]
    return CordRegionTable(
        atlas=atlas_id,
        df=regions[keep].drop_duplicates(subset="parcellation_index"),
        source_csv=paths.regions_csv,
    )
