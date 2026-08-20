"""Assemble spinal cord region-stats tables from intensities and imported points."""

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
from lightsuite.analysis.cord_parcellation import parcellate_cord_intensities
from lightsuite.analysis.intensity_metrics import filter_intensity_metric_rows
from lightsuite.analysis.cord_rollup import apply_cord_rollups
from lightsuite.analysis.cord_hemisphere import (
    cord_hemisphere_side_volume,
    sample_space_hemisphere_flip,
)
from lightsuite.analysis.counts import SAMPLE_POINTS_KEY
from lightsuite.analysis.cord_ontology import load_cord_region_table
from lightsuite.analysis.top_regions import maybe_write_top_n_regions_csv
from lightsuite.atlas.fiederling import resolve_fiederling_paths
from lightsuite.config.models import SpinalCordPipelineConfig
from lightsuite.export.cord_registered import (
    REGISTERED_ANNOTATION_FILENAME,
    compute_registered_annotation_volume,
    discover_registered_cord_paths,
    load_registered_hemisphere_volume,
    load_registered_stack,
    volume_registered_dir,
)
from lightsuite.export.cord_sample_space import (
    ANNOTATION_IN_SAMPLE,
    HEMISPHERE_IN_SAMPLE,
    SEGMENTS_IN_SAMPLE,
    sample_space_dir,
)
from lightsuite.registration.volume import load_registration_volume

console = Console()


def cord_analysis_requested(config: SpinalCordPipelineConfig) -> bool:
    """True when export/import should run region-stats assembly."""
    return bool(
        config.analysis.parcellate_intensities or config.analysis.count_points
    )


def resolve_cord_stats_spaces(
    export_spaces: list[str],
    config: SpinalCordPipelineConfig,
) -> list[str]:
    """Limit stats to spaces that were exported and configured in analysis.stats_spaces."""
    export_set = {str(space).strip().lower() for space in export_spaces}
    configured = [str(space).strip().lower() for space in config.analysis.stats_spaces]
    stats_spaces = [space for space in configured if space in export_set]
    if not stats_spaces:
        stats_spaces = sorted(export_set & {"atlas", "sample"})
    return stats_spaces


def maybe_run_cord_region_stats(
    config: SpinalCordPipelineConfig,
    *,
    export_spaces: list[str],
    run_region_stats: bool | None = None,
) -> CordRegionStatsRunResult | None:
    """Run region stats when analysis is enabled; skip quietly on missing prerequisites."""
    if run_region_stats is False:
        return None
    if not cord_analysis_requested(config):
        return None
    stats_spaces = resolve_cord_stats_spaces(export_spaces, config)
    if not stats_spaces:
        return None
    try:
        return run_cord_region_stats(config, stats_spaces=stats_spaces)
    except (ValueError, FileNotFoundError) as exc:
        console.print(f"[yellow]Region stats skipped:[/yellow] {exc}")
        return None


@dataclass
class CordRegionStatsRunResult:
    combined_path: Path | None = None
    top_n_path: Path | None = None
    rollup_paths: dict[str, Path] = field(default_factory=dict)
    intensity_channels: list[int] = field(default_factory=list)
    count_labels: list[str] = field(default_factory=list)
    n_rows: int = 0


def _load_atlas_annotation(config: SpinalCordPipelineConfig, register_path: Path) -> np.ndarray:
    annotation_path = register_path / REGISTERED_ANNOTATION_FILENAME
    if annotation_path.is_file():
        return tifffile.imread(annotation_path).astype(np.int32, copy=False)
    console.print("[yellow]annotation_registered.tiff missing — computing from native atlas[/yellow]")
    return compute_registered_annotation_volume(config).astype(np.int32, copy=False)


def _resolve_intensity_channels(
    config: SpinalCordPipelineConfig,
    *,
    available: dict[int, Path],
) -> dict[int, Path]:
    wanted = config.analysis.intensity_channels
    if wanted is None:
        return available
    selected = {ich: available[ich] for ich in wanted if ich in available}
    missing = [ich for ich in wanted if ich not in available]
    if missing:
        console.print(
            f"[yellow]Intensity channels not found in export:[/yellow] {missing} "
            f"(available: {sorted(available)})"
        )
    return selected


def run_cord_region_stats(
    config: SpinalCordPipelineConfig,
    *,
    count_points: bool | None = None,
    parcellate_intensities: bool | None = None,
    stats_spaces: list[str] | None = None,
) -> CordRegionStatsRunResult:
    """Assemble tidy region stats from registered intensities and/or imported points."""
    register_path = volume_registered_dir(config)
    if not register_path.is_dir():
        msg = f"Missing {register_path}. Run 'lightsuite spinal export' first."
        raise FileNotFoundError(msg)

    do_counts = config.analysis.count_points if count_points is None else count_points
    do_intensities = (
        config.analysis.parcellate_intensities
        if parcellate_intensities is None
        else parcellate_intensities
    )
    if not do_counts and not do_intensities:
        msg = "Both point counting and intensity parcellation are disabled."
        raise ValueError(msg)

    spaces = stats_spaces if stats_spaces is not None else config.analysis.stats_spaces
    spaces_set = {str(s).strip().lower() for s in spaces}

    region_table = load_cord_region_table(config.atlas)
    atlas_paths = resolve_fiederling_paths(config.atlas.atlas_dir)
    segments = pd.read_csv(atlas_paths.segments_csv)

    wanted_labels = config.analysis.point_labels
    atlas_frames: list[pd.DataFrame] = []
    sample_frames: list[pd.DataFrame] = []
    count_labels: list[str] = []
    intensity_channels: list[int] = []
    rollup_paths: dict[str, Path] = {}

    if "atlas" in spaces_set:
        annotation: np.ndarray | None = None
        hemisphere_side: np.ndarray | None = None
        split_hemispheres = bool(config.analysis.split_hemispheres)
        keep_whole = bool(config.analysis.hemisphere_keep_whole)

        def _ensure_annotation() -> np.ndarray:
            nonlocal annotation, hemisphere_side
            if annotation is None:
                annotation = _load_atlas_annotation(config, register_path)
            if split_hemispheres and hemisphere_side is None:
                hemisphere_mask = load_registered_hemisphere_volume(config, register_path)
                if hemisphere_mask.shape != annotation.shape:
                    msg = (
                        f"Hemisphere mask shape {hemisphere_mask.shape} != annotation {annotation.shape}. "
                        "Re-run 'lightsuite spinal export'."
                    )
                    raise ValueError(msg)
                hemisphere_side = cord_hemisphere_side_volume(
                    hemisphere_mask,
                    annotation,
                    flip=bool(config.analysis.hemisphere_flip),
                )
            return annotation

        if do_intensities:
            try:
                cord_paths = discover_registered_cord_paths(config)
            except FileNotFoundError as exc:
                if do_counts:
                    console.print(f"[yellow]Intensity parcellation skipped:[/yellow] {exc}")
                else:
                    raise
            else:
                channel_paths = _resolve_intensity_channels(config, available=cord_paths.registered_channels)
                if channel_paths:
                    ann = _ensure_annotation()
                    relative_to = config.analysis.relative_intensity_to
                    for ichan, ch_path in sorted(channel_paths.items()):
                        volume = load_registered_stack(ch_path)
                        if volume.shape != ann.shape:
                            msg = (
                                f"{ch_path.name} shape {volume.shape} != annotation {ann.shape}. "
                                "Re-run 'lightsuite spinal export'."
                            )
                            raise ValueError(msg)
                        tidy = parcellate_cord_intensities(
                            volume,
                            ann,
                            segments,
                            sample=config.sample.name,
                            channel=ichan,
                            region_table=region_table,
                            relative_to=relative_to,
                            hemisphere_side=hemisphere_side,
                            split_hemispheres=split_hemispheres,
                            keep_whole=keep_whole,
                            intensity_metrics=config.analysis.intensity_metrics,
                        )
                        tidy = filter_intensity_metric_rows(tidy, config.analysis.intensity_metrics)
                        if len(tidy):
                            per_chan_path = register_path / f"chan{ichan:02d}_region_stats.csv"
                            write_cord_region_stats_csv(per_chan_path, tidy)
                            atlas_frames.append(tidy)
                            intensity_channels.append(ichan)

        if do_counts:
            ann = _ensure_annotation()
            for npz_path in sorted(register_path.glob("*_atlas_coords.npz")):
                label = npz_path.stem.replace("_atlas_coords", "")
                if wanted_labels is not None and label not in wanted_labels:
                    continue
                points = load_atlas_points(npz_path)
                tidy = count_points_in_cord_regions(
                    points,
                    ann,
                    segments,
                    sample=config.sample.name,
                    channel=label,
                    region_table=region_table,
                    hemisphere_side=hemisphere_side,
                    split_hemispheres=split_hemispheres,
                    keep_whole=keep_whole,
                )
                if len(tidy):
                    per_label_path = register_path / f"{label}_region_counts.csv"
                    write_cord_region_stats_csv(per_label_path, tidy)
                    atlas_frames.append(tidy)
                    count_labels.append(label)

    if "sample" in spaces_set and do_counts:
        sample_dir = sample_space_dir(config.sample.save_path.expanduser())
        ann_path = sample_dir / ANNOTATION_IN_SAMPLE
        seg_path = sample_dir / SEGMENTS_IN_SAMPLE
        hem_path = sample_dir / HEMISPHERE_IN_SAMPLE
        split_hemispheres = bool(config.analysis.split_hemispheres)
        keep_whole = bool(config.analysis.hemisphere_keep_whole)
        if ann_path.is_file():
            annotation_sample = load_registration_volume(ann_path).astype(np.int32)
            segments_vol = None
            if seg_path.is_file():
                segments_vol = load_registration_volume(seg_path).astype(np.int32)
            hemisphere_side_sample: np.ndarray | None = None
            if split_hemispheres and hem_path.is_file():
                hemisphere_mask_sample = load_registration_volume(hem_path).astype(np.uint8)
                if hemisphere_mask_sample.shape != annotation_sample.shape:
                    msg = (
                        f"Sample hemisphere shape {hemisphere_mask_sample.shape} != "
                        f"annotation {annotation_sample.shape}. "
                        "Re-run 'lightsuite spinal export --space sample'."
                    )
                    raise ValueError(msg)
                hemisphere_side_sample = cord_hemisphere_side_volume(
                    hemisphere_mask_sample,
                    annotation_sample,
                    flip=sample_space_hemisphere_flip(
                        atlas_hemisphere_flip=bool(config.analysis.hemisphere_flip),
                    ),
                )
            elif split_hemispheres:
                console.print(
                    f"[yellow]Sample-space hemisphere split skipped:[/yellow] missing {hem_path}. "
                    "Re-run 'lightsuite spinal export --space sample'."
                )
            registres_um = float(config.registration.resolution_um)
            voxel_um = (registres_um, registres_um, registres_um)
            for npz_path in sorted(register_path.glob("*_sample_coords.npz")):
                label = npz_path.stem.replace("_sample_coords", "")
                if wanted_labels is not None and label not in wanted_labels:
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
                    hemisphere_side=hemisphere_side_sample,
                    split_hemispheres=split_hemispheres and hemisphere_side_sample is not None,
                    keep_whole=keep_whole,
                )
                if len(tidy):
                    sample_frames.append(tidy)
                    count_labels.append(f"{label}@sample")
        else:
            console.print(
                f"[yellow]Sample-space cord counts skipped:[/yellow] missing {ann_path}"
            )

    combined_path: Path | None = None
    top_n_path: Path | None = None
    if atlas_frames:
        combined = pd.concat(atlas_frames, ignore_index=True).reindex(columns=CORD_TIDY_COLUMNS)
        regions_df = pd.read_csv(atlas_paths.regions_csv)
        combined = apply_cord_rollups(combined, regions_df, config.analysis.rollups)
        combined_path = register_path / "region_stats.csv"
        write_cord_region_stats_csv(combined_path, combined)
        top_n_path = maybe_write_top_n_regions_csv(
            register_path,
            combined,
            n=config.analysis.top_n_regions,
            rank_by=config.analysis.top_n_rank_by,
        )
        for level in config.analysis.rollups:
            normalized = str(level).strip().lower()
            if normalized not in ("division", "structure", "horn"):
                continue
            subset = combined[combined["rollup_level"] == normalized]
            if len(subset):
                rollup_path = register_path / f"region_stats_{normalized}.csv"
                write_cord_region_stats_csv(rollup_path, subset)
                rollup_paths[normalized] = rollup_path
    else:
        combined = pd.DataFrame(columns=CORD_TIDY_COLUMNS)

    if sample_frames:
        sample_combined = pd.concat(sample_frames, ignore_index=True).reindex(columns=CORD_TIDY_COLUMNS)
        sample_dir = sample_space_dir(config.sample.save_path.expanduser())
        write_cord_region_stats_csv(sample_dir / "region_stats_sample.csv", sample_combined)
        maybe_write_top_n_regions_csv(
            sample_dir,
            sample_combined,
            n=config.analysis.top_n_regions,
            rank_by=config.analysis.top_n_rank_by,
            sample_space=True,
        )

    console.print(
        f"[green]Cord region stats:[/green] {len(combined)} atlas rows "
        f"({len(intensity_channels)} intensity channel(s), {len(count_labels)} point source(s)"
        f"{f', rollups: {sorted(rollup_paths)}' if rollup_paths else ''})"
    )
    if top_n_path is not None:
        console.print(
            f"[green]Top-{config.analysis.top_n_regions} regions:[/green] {top_n_path}"
        )
    return CordRegionStatsRunResult(
        combined_path=combined_path,
        top_n_path=top_n_path,
        rollup_paths=rollup_paths,
        intensity_channels=intensity_channels,
        count_labels=count_labels,
        n_rows=int(len(combined)),
    )
