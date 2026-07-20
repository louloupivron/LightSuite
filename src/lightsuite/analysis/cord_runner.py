"""Assemble spinal cord region-stats tables from imported atlas-space points."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile
from rich.console import Console

from lightsuite.analysis.cord_counts import (
    CORD_TIDY_COLUMNS,
    count_points_in_cord_regions,
    load_atlas_points,
    write_cord_region_stats_csv,
)
from lightsuite.analysis.cord_ontology import load_cord_region_table
from lightsuite.atlas.fiederling import resolve_fiederling_paths
from lightsuite.config.models import SpinalCordPipelineConfig
from lightsuite.export.cord_registered import (
    REGISTERED_ANNOTATION_FILENAME,
    compute_registered_annotation_volume,
    volume_registered_dir,
)

console = Console()


@dataclass
class CordRegionStatsRunResult:
    combined_path: Path | None = None
    count_labels: list[str] = field(default_factory=list)
    n_rows: int = 0


def run_cord_region_stats(
    config: SpinalCordPipelineConfig,
    *,
    count_points: bool | None = None,
) -> CordRegionStatsRunResult:
    """Bin imported atlas-space points into Fiederling regions and segments."""
    register_path = volume_registered_dir(config)
    if not register_path.is_dir():
        msg = f"Missing {register_path}. Run 'lightsuite spinal export' or import-annotations first."
        raise FileNotFoundError(msg)

    do_counts = config.analysis.count_points if count_points is None else count_points
    if not do_counts:
        msg = "Point counting disabled (analysis.count_points: false)."
        raise ValueError(msg)

    region_table = load_cord_region_table(config.atlas)
    atlas_paths = resolve_fiederling_paths(config.atlas.atlas_dir)
    segments = pd.read_csv(atlas_paths.segments_csv)

    annotation_path = register_path / REGISTERED_ANNOTATION_FILENAME
    if annotation_path.is_file():
        annotation = tifffile.imread(annotation_path).astype(np.int32, copy=False)
    else:
        console.print(
            "[yellow]annotation_registered.tiff missing — computing from native atlas[/yellow]"
        )
        annotation = compute_registered_annotation_volume(config).astype(np.int32, copy=False)

    npz_paths = sorted(register_path.glob("*_atlas_coords.npz"))
    if not npz_paths:
        msg = (
            f"No *_atlas_coords.npz files in {register_path}. "
            "Run 'lightsuite spinal import-annotations' first."
        )
        raise FileNotFoundError(msg)

    wanted = config.analysis.point_labels
    frames: list[pd.DataFrame] = []
    count_labels: list[str] = []

    for npz_path in npz_paths:
        label = npz_path.stem.replace("_atlas_coords", "")
        if wanted is not None and label not in wanted:
            continue
        points = load_atlas_points(npz_path)
        tidy = count_points_in_cord_regions(
            points,
            annotation,
            segments,
            sample=config.sample.name,
            channel=label,
            region_table=region_table,
        )
        if len(tidy):
            per_label_path = register_path / f"{label}_region_counts.csv"
            write_cord_region_stats_csv(per_label_path, tidy)
            frames.append(tidy)
            count_labels.append(label)

    combined = (
        pd.concat(frames, ignore_index=True).reindex(columns=CORD_TIDY_COLUMNS)
        if frames
        else pd.DataFrame(columns=CORD_TIDY_COLUMNS)
    )
    combined_path: Path | None = None
    if len(combined):
        combined_path = register_path / "region_stats.csv"
        write_cord_region_stats_csv(combined_path, combined.reindex(columns=frames[0].columns))

    console.print(
        f"[green]Cord region stats:[/green] {len(combined)} rows "
        f"({len(count_labels)} point source(s))"
    )
    return CordRegionStatsRunResult(
        combined_path=combined_path,
        count_labels=count_labels,
        n_rows=int(len(combined)),
    )
