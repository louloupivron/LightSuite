"""Orchestrate per-sample region-stats assembly (intensity + cell counts).

Reads the atlas-space artifacts already written by ``brain export`` (intensity
JSON/CSV) and ``brain import-annotations`` (``*_atlas_coords.npz``) and produces a
single tidy ``region_stats.csv`` with region names, plus per-source count tables.
No transformix run is required — this only consumes existing outputs.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from rich.console import Console

from lightsuite.analysis.counts import (
    SAMPLE_POINTS_KEY,
    count_points_in_regions,
    load_atlas_points,
)
from lightsuite.analysis.ontology import RegionTable, load_region_table
from lightsuite.analysis.region_stats import (
    concat_tidy,
    parcellation_result_to_tidy,
    write_region_stats_csv,
)
from lightsuite.atlas.io import load_atlas_volume
from lightsuite.atlas.registry import resolve_brain_atlas_from_config
from lightsuite.config.models import BrainPipelineConfig
from lightsuite.export.brain_export import _load_transform_params
from lightsuite.export.brain_sample_space import ANNOTATION_IN_SAMPLE, sample_space_dir
from lightsuite.export.parcellation import ParcellationResult
from lightsuite.preprocess.checkpoint import RegOptsCheckpoint
from lightsuite.registration.volume import load_registration_volume

console = Console()


@dataclass
class RegionStatsRunResult:
    combined_path: Path | None = None
    intensity_channels: list[int] = field(default_factory=list)
    count_labels: list[str] = field(default_factory=list)
    n_rows: int = 0


def _parcellation_result_from_json(path: Path) -> ParcellationResult:
    data = json.loads(path.read_text(encoding="utf-8"))
    return ParcellationResult(
        area_ids=np.asarray(data["areaidx"], dtype=np.int64),
        median_over_areas=np.asarray(data["medianoverareas"], dtype=np.float32),
        std_over_areas=np.asarray(data["stdoverareas"], dtype=np.float32),
        volume_over_areas=np.asarray(data["volumeoverareas"], dtype=np.float32),
    )


def run_region_stats(
    config: BrainPipelineConfig,
    *,
    count_points: bool | None = None,
    stats_spaces: list[str] | None = None,
) -> RegionStatsRunResult:
    """Assemble a tidy region-stats table for one sample from existing outputs."""
    save_path = config.sample.save_path.expanduser()
    register_path = save_path / "volume_registered"
    if not register_path.is_dir():
        msg = f"Missing {register_path}. Run 'lightsuite brain export' first."
        raise FileNotFoundError(msg)

    transform_params = _load_transform_params(save_path)
    atlas = resolve_brain_atlas_from_config(config.atlas)
    regopts_path = save_path / "regopts.json"
    checkpoint = RegOptsCheckpoint.load(regopts_path) if regopts_path.is_file() else None

    spaces = stats_spaces if stats_spaces is not None else config.analysis.stats_spaces
    spaces_set = {str(s).strip().lower() for s in spaces}

    region_table: RegionTable | None = None
    try:
        region_table = load_region_table(atlas)
    except FileNotFoundError as exc:
        console.print(f"[yellow]Region names unavailable:[/yellow] {exc}")

    atlas_frames: list = []
    sample_frames: list = []
    intensity_channels: list[int] = []
    count_labels: list[str] = []

    if "atlas" in spaces_set:
        for json_path in sorted(register_path.glob("chan*_intensities.json")):
            stem = json_path.stem
            try:
                channel = int(stem.replace("chan", "").split("_")[0])
            except ValueError:
                channel = stem
            result = _parcellation_result_from_json(json_path)
            tidy = parcellation_result_to_tidy(
                result,
                region_table,
                sample=config.sample.name,
                channel=channel,
                atlas=atlas.brain_atlas,
            )
            atlas_frames.append(tidy)
            if isinstance(channel, int):
                intensity_channels.append(channel)

    do_counts = config.analysis.count_points if count_points is None else count_points
    if do_counts and "atlas" in spaces_set:
        npz_paths = sorted(register_path.glob("*_atlas_coords.npz"))
        wanted = config.analysis.point_labels
        annotation = None
        for npz_path in npz_paths:
            label = npz_path.stem.replace("_atlas_coords", "")
            if wanted is not None and label not in wanted:
                continue
            if annotation is None:
                annotation = load_atlas_volume(atlas.annotation_path)
            points = load_atlas_points(npz_path)
            tidy = count_points_in_regions(
                points,
                annotation,
                atlas_id=atlas.brain_atlas,
                atlas_resolution_um=transform_params.atlas_resolution_um,
                sample=config.sample.name,
                channel=label,
                region_table=region_table,
                brainglobe_name=atlas.brainglobe_name,
            )
            if len(tidy):
                atlas_frames.append(tidy)
                count_labels.append(label)

    if do_counts and "sample" in spaces_set:
        sample_dir = sample_space_dir(save_path)
        ann_path = sample_dir / ANNOTATION_IN_SAMPLE
        if ann_path.is_file() and checkpoint is not None:
            from lightsuite.export.brain_export import _load_transform_params
            from lightsuite.export.brain_sample_space import (
                atlas_volumes_permuted_on_disk,
                load_sample_space_atlas_volume,
            )

            transform_params = _load_transform_params(save_path)
            permute = transform_params.permute_sample_to_atlas or [1, 2, 3]
            annotation_sample = load_sample_space_atlas_volume(
                ann_path,
                permute,
                permuted_on_disk=atlas_volumes_permuted_on_disk(sample_dir),
            ).astype(np.int32, copy=False)
            registres_um = float(checkpoint.registres_um)
            for npz_path in sorted(register_path.glob("*_sample_coords.npz")):
                label = npz_path.stem.replace("_sample_coords", "")
                if config.analysis.point_labels is not None and label not in config.analysis.point_labels:
                    continue
                points = load_atlas_points(npz_path, key=SAMPLE_POINTS_KEY)
                tidy = count_points_in_regions(
                    points,
                    annotation_sample,
                    atlas_id=atlas.brain_atlas,
                    atlas_resolution_um=registres_um,
                    sample=config.sample.name,
                    channel=label,
                    region_table=region_table,
                    brainglobe_name=atlas.brainglobe_name,
                )
                if len(tidy):
                    sample_frames.append(tidy)
                    count_labels.append(f"{label}@sample")
        else:
            console.print(
                f"[yellow]Sample-space counts skipped:[/yellow] missing {ann_path} "
                "or regopts.json (run export with sample space)."
            )

    combined_path: Path | None = None
    if atlas_frames:
        combined = concat_tidy(atlas_frames)
        combined_path = register_path / "region_stats.csv"
        write_region_stats_csv(combined_path, combined)
    else:
        combined = concat_tidy([])

    if sample_frames:
        sample_combined = concat_tidy(sample_frames)
        sample_stats_path = sample_space_dir(save_path) / "region_stats_sample.csv"
        write_region_stats_csv(sample_stats_path, sample_combined)

    console.print(
        f"[green]Region stats:[/green] {len(combined)} atlas rows "
        f"({len(intensity_channels)} intensity channel(s), {len(count_labels)} point source(s))"
    )
    return RegionStatsRunResult(
        combined_path=combined_path,
        intensity_channels=intensity_channels,
        count_labels=count_labels,
        n_rows=int(len(combined)),
    )
