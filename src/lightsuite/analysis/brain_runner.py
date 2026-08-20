"""Assemble brain region-stats tables from intensities and imported point clouds."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
from rich.console import Console

from lightsuite.analysis.counts import SAMPLE_POINTS_KEY, count_points_in_regions, load_atlas_points
from lightsuite.analysis.intensity_metrics import filter_intensity_metric_rows
from lightsuite.analysis.ontology import RegionTable
from lightsuite.analysis.region_stats import concat_tidy, write_region_stats_csv
from lightsuite.atlas.io import load_atlas_volume
from lightsuite.atlas.registry import AtlasPaths
from lightsuite.config.models import BrainPipelineConfig
from lightsuite.registration.brain_paths import (
    brain_stats_dir,
    iter_brain_import_paths,
    iter_brain_stats_paths,
)

console = Console()


@dataclass
class BrainRegionStatsRunResult:
    combined_path: Path | None = None
    count_labels: list[str] = field(default_factory=list)
    n_rows: int = 0


def _wanted_point_labels(config: BrainPipelineConfig) -> set[str] | None:
    labels = config.analysis.point_labels
    if labels is None:
        return None
    return {str(label).strip() for label in labels if str(label).strip()}


def _load_brain_annotation(atlas: AtlasPaths, shape: tuple[int, ...]) -> np.ndarray:
    annotation = np.asarray(load_atlas_volume(atlas.annotation_path)).astype(np.int32, copy=False)
    if tuple(annotation.shape) != tuple(shape):
        msg = (
            f"Atlas annotation shape {annotation.shape} != registered volume shape {shape}. "
            "Re-run brain export."
        )
        raise ValueError(msg)
    return annotation


def collect_brain_point_count_frames(
    config: BrainPipelineConfig,
    stats_path: Path,
    *,
    annotation: np.ndarray,
    region_table: RegionTable | None,
    atlas: AtlasPaths,
    transform_params: object,
    save_path: Path | None = None,
) -> list[pd.DataFrame]:
    """Return tidy cell_count / cell_density rows for atlas-space imports."""
    if not config.analysis.count_points:
        return []
    wanted = _wanted_point_labels(config)
    frames: list[pd.DataFrame] = []
    ml_axis = 2
    brainglobe_name = atlas.brainglobe_name
    npz_paths = (
        iter_brain_import_paths(save_path, "*_atlas_coords.npz")
        if save_path is not None
        else sorted(stats_path.glob("*_atlas_coords.npz"))
    )
    for npz_path in npz_paths:
        label = npz_path.stem.replace("_atlas_coords", "")
        if wanted is not None and label not in wanted:
            continue
        points = load_atlas_points(npz_path)
        tidy = count_points_in_regions(
            points,
            annotation,
            atlas_id=atlas.brain_atlas,
            atlas_resolution_um=float(transform_params.atlas_resolution_um),
            sample=config.sample.name,
            channel=label,
            region_table=region_table,
            ml_axis=ml_axis,
            brainglobe_name=brainglobe_name,
        )
        if len(tidy):
            per_label_path = stats_path / f"{label}_region_counts.csv"
            write_region_stats_csv(per_label_path, tidy)
            frames.append(tidy)
    return frames


def collect_brain_sample_point_count_frames(
    config: BrainPipelineConfig,
    save_path: Path,
    *,
    annotation_sample: np.ndarray,
    region_table: RegionTable | None,
    atlas: AtlasPaths,
    registres_um: float,
    transform_params: object,
) -> list[pd.DataFrame]:
    """Return tidy count rows for registration-grid sample-space imports."""
    if not config.analysis.count_points:
        return []
    if "sample" not in {str(s).strip().lower() for s in config.analysis.stats_spaces}:
        return []
    wanted = _wanted_point_labels(config)
    frames: list[pd.DataFrame] = []
    for npz_path in iter_brain_import_paths(save_path, "*_sample_coords.npz"):
        label = npz_path.stem.replace("_sample_coords", "")
        if wanted is not None and label not in wanted:
            continue
        points = load_atlas_points(npz_path, key=SAMPLE_POINTS_KEY)
        tidy = count_points_in_regions(
            points,
            annotation_sample,
            atlas_id=atlas.brain_atlas,
            atlas_resolution_um=float(registres_um),
            sample=config.sample.name,
            channel=f"{label}@sample",
            region_table=region_table,
            ml_axis=2,
            brainglobe_name=atlas.brainglobe_name,
        )
        if len(tidy):
            frames.append(tidy)
    return frames


def finalize_brain_region_stats(
    config: BrainPipelineConfig,
    stats_path: Path,
    tidy_frames: list[pd.DataFrame],
    *,
    annotation: np.ndarray | None = None,
    region_table: RegionTable | None = None,
    atlas: AtlasPaths | None = None,
    transform_params: object | None = None,
    save_path: Path | None = None,
) -> BrainRegionStatsRunResult:
    """Filter intensity metrics, append point counts, and write combined region_stats.csv."""
    metrics = config.analysis.intensity_metrics
    filtered = [filter_intensity_metric_rows(frame, metrics) for frame in tidy_frames]
    count_labels: list[str] = []

    if (
        config.analysis.count_points
        and annotation is not None
        and atlas is not None
        and transform_params is not None
    ):
        count_frames = collect_brain_point_count_frames(
            config,
            stats_path,
            annotation=annotation,
            region_table=region_table,
            atlas=atlas,
            transform_params=transform_params,
            save_path=save_path,
        )
        filtered.extend(count_frames)
        count_labels.extend(
            str(frame["channel"].iloc[0]) for frame in count_frames if len(frame) > 0
        )

    usable = [frame for frame in filtered if frame is not None and len(frame) > 0]
    combined_path: Path | None = None
    n_rows = 0
    if usable:
        combined = concat_tidy(usable)
        combined_path = stats_path / "region_stats.csv"
        write_region_stats_csv(combined_path, combined)
        n_rows = len(combined)
        console.print(
            f"[green]Brain region stats:[/green] {n_rows} rows "
            f"({len(tidy_frames)} intensity source(s), {len(count_labels)} point source(s))"
        )
    return BrainRegionStatsRunResult(
        combined_path=combined_path,
        count_labels=count_labels,
        n_rows=n_rows,
    )


def maybe_refresh_brain_region_stats(
    config: BrainPipelineConfig,
    *,
    export_spaces: list[str],
) -> BrainRegionStatsRunResult | None:
    """Rebuild region_stats.csv from existing per-channel tidy exports and point imports."""
    if not config.analysis.count_points:
        return None
    if "atlas" not in {str(s).strip().lower() for s in export_spaces}:
        return None
    save_path = config.sample.save_path.expanduser()
    stats_path = brain_stats_dir(save_path)

    from lightsuite.atlas.registry import resolve_brain_atlas_from_config
    from lightsuite.export.brain_export import _load_transform_params

    atlas = resolve_brain_atlas_from_config(config.atlas)
    transform_params = _load_transform_params(save_path)
    tidy_frames: list[pd.DataFrame] = []
    for path in iter_brain_stats_paths(save_path, "chan*_region_stats.csv"):
        frame = pd.read_csv(path)
        if len(frame):
            tidy_frames.append(frame)

    if not tidy_frames and not iter_brain_import_paths(save_path, "*_atlas_coords.npz"):
        return None

    try:
        shape = tuple(int(v) for v in transform_params.atlassize)
        annotation = _load_brain_annotation(atlas, shape)
    except (ValueError, OSError) as exc:
        console.print(f"[yellow]Brain region stats refresh skipped:[/yellow] {exc}")
        return None

    region_table = None
    try:
        from lightsuite.analysis.ontology import load_region_table

        region_table = load_region_table(atlas)
    except FileNotFoundError:
        pass

    return finalize_brain_region_stats(
        config,
        stats_path,
        tidy_frames,
        annotation=annotation,
        region_table=region_table,
        atlas=atlas,
        transform_params=transform_params,
        save_path=save_path,
    )


def maybe_write_brain_sample_region_stats(
    config: BrainPipelineConfig,
    *,
    tidy_frames: list[pd.DataFrame],
    annotation_sample: np.ndarray | None,
    region_table: RegionTable | None,
    atlas: AtlasPaths,
    registres_um: float,
    transform_params: object,
    stats_path: Path,
    save_path: Path,
) -> Path | None:
    """Write sample-space region_stats_sample.csv with optional point counts."""
    metrics = config.analysis.intensity_metrics
    filtered = [filter_intensity_metric_rows(frame, metrics) for frame in tidy_frames]
    if annotation_sample is not None:
        filtered.extend(
            collect_brain_sample_point_count_frames(
                config,
                save_path,
                annotation_sample=annotation_sample,
                region_table=region_table,
                atlas=atlas,
                registres_um=registres_um,
                transform_params=transform_params,
            )
        )
    usable = [frame for frame in filtered if frame is not None and len(frame) > 0]
    if not usable:
        return None
    combined_path = stats_path / "region_stats_sample.csv"
    write_region_stats_csv(combined_path, concat_tidy(usable))
    return combined_path
