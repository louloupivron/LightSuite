"""Division-level atlas label volumes for interactive QC and masking.

Remaps fine atlas annotation IDs to coarse **division** integers (0 = unassigned)
using the Allen ABC membership table or the Perens → Allen ontology mapping.
Outputs are cached beside the atlas NIfTIs and reused across samples.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile

from lightsuite.analysis.ontology import (
    ALLEN_MEMBERSHIP_FILENAME,
    RegionTable,
    _resolve_allen_membership_csv,
    build_ccf_to_division,
    load_region_table,
    resolve_ccf_by_hierarchy,
    structures_parent_map,
)
from lightsuite.atlas.io import load_atlas_volume
from lightsuite.atlas.registry import (
    AtlasPaths,
    atlas_resolution_um_for_cache,
    uses_ccf_id_parcellation,
)

LEGEND_COLUMNS = ("division_id", "division_acronym", "division_name")


@dataclass(frozen=True)
class DivisionMapPaths:
    labels_tiff: Path
    legend_csv: Path


@dataclass(frozen=True)
class DivisionMapResult:
    labels: np.ndarray
    legend: pd.DataFrame
    paths: DivisionMapPaths


def division_cache_paths(atlas: AtlasPaths, *, fallback_resolution_um: float = 20.0) -> DivisionMapPaths:
    """Atlas-side cache paths for the division label volume and legend."""
    res = int(round(atlas_resolution_um_for_cache(atlas, fallback_resolution_um)))
    return DivisionMapPaths(
        labels_tiff=atlas.atlas_dir / f"division_labels_{res}um.tif",
        legend_csv=atlas.atlas_dir / f"division_id_legend_{res}um.csv",
    )


def ensure_division_map(atlas: AtlasPaths, *, force: bool = False) -> DivisionMapResult:
    """Load or build the cached division label volume for an atlas."""
    paths = division_cache_paths(atlas)
    if not force and paths.labels_tiff.is_file() and paths.legend_csv.is_file():
        labels = np.asarray(tifffile.imread(paths.labels_tiff), dtype=np.int32)
        legend = pd.read_csv(paths.legend_csv)
        return DivisionMapResult(labels=labels, legend=legend, paths=paths)

    labels, legend = build_division_labels(atlas)
    write_division_map(paths, labels, legend)
    return DivisionMapResult(labels=labels, legend=legend, paths=paths)


def build_division_labels(atlas: AtlasPaths) -> tuple[np.ndarray, pd.DataFrame]:
    """Build a division label volume (Y, X, Z) and legend table."""
    annotation = load_atlas_volume(atlas.annotation_path)
    name_to_acronym: dict[str, str] = {}
    if atlas.brain_atlas == "allen" and atlas.atlas_source == "files":
        membership = _resolve_allen_membership_csv(atlas.atlas_dir)
        if membership is None:
            msg = f"Allen membership CSV not found for division map ({ALLEN_MEMBERSHIP_FILENAME})."
            raise FileNotFoundError(msg)
        index_to_division = _allen_index_to_division(membership)
        name_to_acronym = _allen_division_name_to_acronym(membership)
    else:
        region_table = load_region_table(atlas)
        if uses_ccf_id_parcellation(atlas):
            index_to_division = _ccf_id_to_division_map(atlas, region_table)
        else:
            index_to_division = _perens_ccf_to_division(region_table)
        allen_csv = _resolve_allen_membership_csv(None)
        if allen_csv is not None:
            name_to_acronym = _allen_division_name_to_acronym(allen_csv)

    return _remap_annotation_to_divisions(
        annotation, index_to_division, name_to_acronym=name_to_acronym
    )


def write_division_map(
    paths: DivisionMapPaths,
    labels: np.ndarray,
    legend: pd.DataFrame,
) -> None:
    paths.labels_tiff.parent.mkdir(parents=True, exist_ok=True)
    tifffile.imwrite(paths.labels_tiff, labels.astype(np.int32))
    legend.to_csv(paths.legend_csv, index=False)


def _allen_index_to_division(membership_csv: Path) -> dict[int, str]:
    df = pd.read_csv(membership_csv)
    required = {"parcellation_index", "parcellation_term_set_name", "parcellation_term_name"}
    missing = required - set(df.columns)
    if missing:
        msg = f"Allen membership CSV missing columns: {sorted(missing)}"
        raise ValueError(msg)
    div = df[df["parcellation_term_set_name"] == "division"].drop_duplicates(
        subset=["parcellation_index"]
    )
    return {
        int(row["parcellation_index"]): str(row["parcellation_term_name"])
        for _, row in div.iterrows()
    }


def _perens_ccf_to_division(region_table: RegionTable) -> dict[int, str]:
    df = region_table.df.dropna(subset=["division"])
    return {
        int(row["parcellation_index"]): str(row["division"])
        for _, row in df.iterrows()
    }


def _ccf_id_to_division_map(
    atlas: AtlasPaths,
    region_table: RegionTable,
) -> dict[int, str]:
    """Map annotation voxel ids (Allen CCF ontology) to division names."""
    df = region_table.df.dropna(subset=["division", "ccf_id"])
    mapping: dict[int, str] = {
        int(row["ccf_id"]): str(row["division"])
        for _, row in df.iterrows()
    }

    allen_csv = _resolve_allen_membership_csv(atlas.atlas_dir)
    if allen_csv is None or not allen_csv.is_file():
        allen_csv = _resolve_allen_membership_csv(None)
    if allen_csv is None or not allen_csv.is_file():
        return mapping

    ccf_to_div = build_ccf_to_division(allen_csv)
    for ccf_id, division in ccf_to_div.items():
        mapping.setdefault(int(ccf_id), str(division))

    if atlas.structures_csv_path is None or not atlas.structures_csv_path.is_file():
        return mapping

    structures_df = pd.read_csv(atlas.structures_csv_path)
    parent_map = structures_parent_map(structures_df)
    for ccf_id in structures_df["id"].astype(int):
        if ccf_id in mapping:
            continue
        division = resolve_ccf_by_hierarchy(ccf_id, ccf_to_div, parent_map)
        if division is not None:
            mapping[int(ccf_id)] = str(division)
    return mapping


def _allen_division_name_to_acronym(membership_csv: Path) -> dict[str, str]:
    df = pd.read_csv(membership_csv)
    div = df[df["parcellation_term_set_name"] == "division"].drop_duplicates(
        subset=["parcellation_term_name"]
    )
    return {
        str(row["parcellation_term_name"]): str(row["parcellation_term_acronym"])
        for _, row in div.iterrows()
    }


def _remap_annotation_to_divisions(
    annotation: np.ndarray,
    index_to_division: dict[int, str],
    *,
    name_to_acronym: dict[str, str] | None = None,
) -> tuple[np.ndarray, pd.DataFrame]:
    """Remap fine annotation labels to integer division IDs."""
    av = np.asanyarray(annotation)
    unique_labels = np.unique(av)

    division_names: list[str] = []
    for label in unique_labels:
        label_int = int(label)
        if label_int == 0:
            division_names.append("unassigned")
        else:
            division_names.append(index_to_division.get(label_int, "unassigned"))

    unique_div_names = sorted(set(division_names))
    name_to_id: dict[str, int] = {}
    next_id = 1
    for name in unique_div_names:
        if name == "unassigned":
            name_to_id[name] = 0
        else:
            name_to_id[name] = next_id
            next_id += 1

    label_to_div_id = {
        int(lab): name_to_id[name]
        for lab, name in zip(unique_labels, division_names, strict=True)
    }

    lookup = np.zeros(unique_labels.shape, dtype=np.int64)
    for i, lab in enumerate(unique_labels):
        lookup[i] = label_to_div_id[int(lab)]

    inverse = np.searchsorted(unique_labels, av)
    division_labels = lookup[inverse].astype(np.int32)

    acronym_lookup = name_to_acronym or {}
    legend_rows = []
    for name in sorted(name_to_id, key=lambda n: name_to_id[n]):
        did = name_to_id[name]
        legend_rows.append(
            {
                "division_id": did,
                "division_acronym": acronym_lookup.get(name, name),
                "division_name": name,
            }
        )
    legend = pd.DataFrame(legend_rows, columns=list(LEGEND_COLUMNS))
    return division_labels, legend
