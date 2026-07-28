"""Assemble spinal cord region-stats tables from imported points."""

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
from lightsuite.analysis.counts import SAMPLE_POINTS_KEY
from lightsuite.analysis.cord_ontology import load_cord_region_table
from lightsuite.atlas.fiederling import resolve_fiederling_paths
from lightsuite.config.models import SpinalCordPipelineConfig
from lightsuite.export.cord_registered import (
    REGISTERED_ANNOTATION_FILENAME,
    compute_registered_annotation_volume,
    volume_registered_dir,
)
from lightsuite.export.cord_sample_space import (
    ANNOTATION_IN_SAMPLE,
    SEGMENTS_IN_SAMPLE,
    sample_space_dir,
)
from lightsuite.registration.volume import load_registration_volume

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
    stats_spaces: list[str] | None = None,
) -> CordRegionStatsRunResult:
    """Bin imported points into Fiederling regions and segments."""
    register_path = volume_registered_dir(config)
    if not register_path.is_dir():
        msg = f"Missing {register_path}. Run 'lightsuite spinal export' or import-annotations first."
        raise FileNotFoundError(msg)

    do_counts = config.analysis.count_points if count_points is None else count_points
    if not do_counts:
        msg = "Point counting disabled (analysis.count_points: false)."
        raise ValueError(msg)

    spaces = stats_spaces if stats_spaces is not None else config.analysis.stats_spaces
    spaces_set = {str(s).strip().lower() for s in spaces}

    region_table = load_cord_region_table(config.atlas)
    atlas_paths = resolve_fiederling_paths(config.atlas.atlas_dir)
    segments = pd.read_csv(atlas_paths.segments_csv)

    wanted = config.analysis.point_labels
    atlas_frames: list[pd.DataFrame] = []
    sample_frames: list[pd.DataFrame] = []
    count_labels: list[str] = []

    if "atlas" in spaces_set:
        annotation_path = register_path / REGISTERED_ANNOTATION_FILENAME
        if annotation_path.is_file():
            annotation = tifffile.imread(annotation_path).astype(np.int32, copy=False)
        else:
            console.print(
                "[yellow]annotation_registered.tiff missing — computing from native atlas[/yellow]"
            )
            annotation = compute_registered_annotation_volume(config).astype(np.int32, copy=False)

        for npz_path in sorted(register_path.glob("*_atlas_coords.npz")):
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
                atlas_frames.append(tidy)
                count_labels.append(label)

    if "sample" in spaces_set:
        sample_dir = sample_space_dir(config.sample.save_path.expanduser())
        ann_path = sample_dir / ANNOTATION_IN_SAMPLE
        seg_path = sample_dir / SEGMENTS_IN_SAMPLE
        if ann_path.is_file():
            annotation_sample = load_registration_volume(ann_path).astype(np.int32)
            segments_vol = None
            if seg_path.is_file():
                segments_vol = load_registration_volume(seg_path).astype(np.int32)
            registres_um = float(config.registration.resolution_um)
            voxel_um = (registres_um, registres_um, registres_um)
            for npz_path in sorted(register_path.glob("*_sample_coords.npz")):
                label = npz_path.stem.replace("_sample_coords", "")
                if wanted is not None and label not in wanted:
                    continue
                points = load_atlas_points(npz_path, key=SAMPLE_POINTS_KEY)
                tidy = count_points_in_cord_regions(
                    points,
                    annotation_sample,
                    segments,
                    sample=config.sample.name,
                    channel=label,
                    region_table=region_table,
                    voxel_um_yxz=voxel_um,
                    segments_volume=segments_vol,
                )
                if len(tidy):
                    sample_frames.append(tidy)
                    count_labels.append(f"{label}@sample")
        else:
            console.print(
                f"[yellow]Sample-space cord counts skipped:[/yellow] missing {ann_path}"
            )

    combined_path: Path | None = None
    if atlas_frames:
        combined = pd.concat(atlas_frames, ignore_index=True).reindex(columns=CORD_TIDY_COLUMNS)
        combined_path = register_path / "region_stats.csv"
        write_cord_region_stats_csv(combined_path, combined)
    else:
        combined = pd.DataFrame(columns=CORD_TIDY_COLUMNS)

    if sample_frames:
        sample_combined = pd.concat(sample_frames, ignore_index=True).reindex(columns=CORD_TIDY_COLUMNS)
        write_cord_region_stats_csv(
            sample_space_dir(config.sample.save_path.expanduser()) / "region_stats_sample.csv",
            sample_combined,
        )

    console.print(
        f"[green]Cord region stats:[/green] {len(combined)} atlas rows "
        f"({len(count_labels)} point source(s))"
    )
    return CordRegionStatsRunResult(
        combined_path=combined_path,
        count_labels=count_labels,
        n_rows=int(len(combined)),
    )
