"""Shared region metadata table keyed on the Allen ontology.

Both supported atlases live in the Allen ontology id space:

- **Allen ABC** ``annotation_10.nii.gz`` stores a dense ``parcellation_index`` per
  voxel; ``parcellation_to_parcellation_term_membership.csv`` maps each index to
  acronyms/names at the *substructure / structure / division* levels and carries
  the underlying AllenCCF ontology id (``ccf_id``).
- **Perens (Gubra LSFM)** ``gubra_ano_olf.nii.gz`` stores raw AllenCCF structure
  ids directly; ``ARA2_annotation_info.csv`` provides ``id`` (== ``ccf_id``),
  ``acronym`` and ``name``.

``load_region_table`` returns a uniform table so downstream analysis (tidy region
stats, cell counts, heatmaps, cross-subject grouping) can join on
``parcellation_index`` and group by a common Allen ``division`` regardless of the
atlas a sample was registered to.
"""

from __future__ import annotations

import os
import re
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from lightsuite.atlas.registry import AtlasPaths

ALLEN_MEMBERSHIP_FILENAME = "parcellation_to_parcellation_term_membership.csv"

#: Columns guaranteed to exist on :attr:`RegionTable.df`.
REGION_TABLE_COLUMNS = (
    "parcellation_index",
    "acronym",
    "name",
    "structure",
    "division",
    "division_acronym",
    "ccf_id",
)


@dataclass(frozen=True)
class RegionTable:
    """Per-region metadata for one atlas, in a common Allen ontology vocabulary."""

    atlas: str
    df: pd.DataFrame
    source_csv: Path | None = None

    def metadata_by_index(self) -> pd.DataFrame:
        """Return a copy indexed by ``parcellation_index`` for fast joins."""
        return self.df.set_index("parcellation_index")


def _parse_ccf_id(term_label: object) -> int | None:
    """Extract the trailing integer ontology id from an AllenCCF term label."""
    match = re.search(r"(\d+)\s*$", str(term_label))
    return int(match.group(1)) if match else None


def _resolve_allen_membership_csv(atlas_dir: Path | None) -> Path | None:
    """Locate the Allen ABC membership CSV (atlas dir → env → cwd glob)."""
    if atlas_dir is not None:
        candidate = atlas_dir / ALLEN_MEMBERSHIP_FILENAME
        if candidate.is_file():
            return candidate
    env = os.environ.get("LIGHTSUITE_ALLEN_PARCELLATION_CSV", "").strip()
    if env:
        env_path = Path(env).expanduser()
        if env_path.is_file():
            return env_path
    for base in (atlas_dir, Path.cwd()):
        if base is None:
            continue
        for found in base.glob(f"**/{ALLEN_MEMBERSHIP_FILENAME}"):
            return found
    return None


def _pivot_allen_membership(df: pd.DataFrame) -> pd.DataFrame:
    """Pivot the long-form ABC membership table to one row per parcellation_index."""
    required = {
        "parcellation_index",
        "parcellation_term_set_name",
        "parcellation_term_name",
        "parcellation_term_acronym",
    }
    missing = required - set(df.columns)
    if missing:
        msg = f"Allen membership CSV missing columns: {sorted(missing)}"
        raise ValueError(msg)

    has_label = "parcellation_term_label" in df.columns

    def _level(name: str) -> pd.DataFrame:
        sel = df[df["parcellation_term_set_name"] == name]
        return sel.drop_duplicates(subset=["parcellation_index"]).set_index("parcellation_index")

    substructure = _level("substructure")
    structure = _level("structure")
    division = _level("division")

    out = pd.DataFrame(index=substructure.index.union(structure.index).union(division.index))
    out.index.name = "parcellation_index"

    out["acronym"] = substructure["parcellation_term_acronym"]
    out["name"] = substructure["parcellation_term_name"]
    out["structure"] = structure["parcellation_term_acronym"]
    out["division"] = division["parcellation_term_name"]
    out["division_acronym"] = division["parcellation_term_acronym"]
    if has_label:
        out["ccf_id"] = substructure["parcellation_term_label"].map(_parse_ccf_id)
    else:
        out["ccf_id"] = pd.NA

    return out.reset_index()


def load_allen_region_table(csv_path: Path) -> RegionTable:
    """Build a :class:`RegionTable` from the Allen ABC membership CSV."""
    df = pd.read_csv(csv_path)
    table = _pivot_allen_membership(df)
    table["parcellation_index"] = table["parcellation_index"].astype("int64")
    return RegionTable(atlas="allen", df=table, source_csv=csv_path)


def build_ccf_to_division(allen_csv: Path) -> dict[int, str]:
    """Map every AllenCCF ontology id (leaf/structure/division) to a division name.

    Used to translate other atlases (e.g. Perens, which stores raw CCF ids) into
    the common Allen division grouping.
    """
    df = pd.read_csv(allen_csv)
    if "parcellation_term_label" not in df.columns:
        return {}
    div_name = (
        df[df["parcellation_term_set_name"] == "division"]
        .drop_duplicates(subset=["parcellation_index"])
        .set_index("parcellation_index")["parcellation_term_name"]
    )
    mapping: dict[int, str] = {}
    for level in ("substructure", "structure", "division"):
        rows = df[df["parcellation_term_set_name"] == level]
        for pidx, label in zip(rows["parcellation_index"], rows["parcellation_term_label"]):
            ccf = _parse_ccf_id(label)
            if ccf is None:
                continue
            name = div_name.get(pidx)
            if isinstance(name, str) and name and name != "unassigned":
                mapping.setdefault(ccf, name)
    return mapping


def load_perens_region_table(
    csv_path: Path,
    *,
    allen_membership_csv: Path | None = None,
    atlas_id: str = "perens",
) -> RegionTable:
    """Build a :class:`RegionTable` from the Perens ``ARA2_annotation_info`` CSV.

    Perens ``id`` values are AllenCCF ontology ids, so ``parcellation_index`` and
    ``ccf_id`` are identical here. When an Allen membership CSV is available the
    common Allen ``division`` is filled in via :func:`build_ccf_to_division`.
    """
    df = pd.read_csv(csv_path)
    required = {"id", "acronym", "name"}
    missing = required - set(df.columns)
    if missing:
        msg = f"Perens annotation CSV missing columns: {sorted(missing)}"
        raise ValueError(msg)

    table = pd.DataFrame(
        {
            "parcellation_index": df["id"].astype("int64"),
            "acronym": df["acronym"].astype(str),
            "name": df["name"].astype(str),
            "structure": pd.NA,
            "division": pd.NA,
            "division_acronym": pd.NA,
            "ccf_id": df["id"].astype("int64"),
        }
    )

    if allen_membership_csv is not None and allen_membership_csv.is_file():
        ccf_to_div = build_ccf_to_division(allen_membership_csv)
        if ccf_to_div:
            table["division"] = table["ccf_id"].map(ccf_to_div)

    return RegionTable(atlas=atlas_id, df=table, source_csv=csv_path)


def load_region_table(
    atlas: AtlasPaths,
    *,
    allen_membership_csv: Path | None = None,
) -> RegionTable:
    """Load the region metadata table for a resolved atlas.

    Parameters
    ----------
    atlas:
        Resolved atlas paths (see :func:`lightsuite.atlas.registry.resolve_brain_atlas`).
    allen_membership_csv:
        Optional explicit path to the Allen ABC membership CSV. Used as the Allen
        source and, for the Perens atlas, as the CCF→division translation source.
        When omitted it is auto-resolved from the atlas dir / env / cwd.
    """
    if atlas.brain_atlas == "allen":
        csv_path = allen_membership_csv or _resolve_allen_membership_csv(atlas.atlas_dir)
        if csv_path is None or not csv_path.is_file():
            msg = (
                "Allen membership CSV not found. Set LIGHTSUITE_ALLEN_PARCELLATION_CSV "
                f"or place {ALLEN_MEMBERSHIP_FILENAME} on the atlas path."
            )
            raise FileNotFoundError(msg)
        return load_allen_region_table(csv_path)

    if atlas.brain_atlas in ("perens", "princeton") or atlas.atlas_source == "brainglobe":
        if atlas.structures_csv_path is None:
            msg = "Atlas structures CSV not found (ARA2 or BrainGlobe structures.csv)."
            raise FileNotFoundError(msg)
        allen_csv = allen_membership_csv or _resolve_allen_membership_csv(None)
        return load_perens_region_table(
            atlas.structures_csv_path,
            allen_membership_csv=allen_csv,
            atlas_id=atlas.brain_atlas,
        )

    msg = f"No region table loader for atlas '{atlas.brain_atlas}'."
    raise ValueError(msg)
