"""Orchestrate naive unassigned-voxel registration QC for one sample."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile
from rich.console import Console

from lightsuite.analysis.division_map import ensure_division_map
from lightsuite.analysis.registration_qc import (
    RegistrationQcResult,
    compute_unassigned_registration_score,
    score_to_dataframe,
    threshold_sweep,
    write_threshold_sweep_plot,
)
from lightsuite.atlas.registry import resolve_brain_atlas_with_config
from lightsuite.config.models import BrainPipelineConfig
from lightsuite.export.brain_export import _load_transform_params
from lightsuite.gui.view_divisions_brain import discover_division_viewer_paths
from lightsuite.registration.volume import load_registration_volume

console = Console()


def _qc_output_dir(config: BrainPipelineConfig) -> Path:
    return config.sample.save_path.expanduser() / "registration_qc"


def _resolve_channel(paths, channel: int) -> tuple[str, Path]:
    key = f"channel {channel}"
    if key not in paths.channel_paths:
        available = ", ".join(sorted(paths.channel_paths))
        msg = f"Channel {channel} not found. Available: {available}"
        raise KeyError(msg)
    return key, paths.channel_paths[key]


def load_registration_qc_volumes(
    config: BrainPipelineConfig,
    *,
    channel: int = 1,
    space: str = "atlas",
) -> tuple[np.ndarray, np.ndarray, Path, Path]:
    """Load registered channel volume and division labels for QC."""
    space_key = space.strip().lower()
    if space_key == "sample":
        from lightsuite.export.brain_sample_space import sample_space_dir

        save_path = config.sample.save_path.expanduser()
        sample_dir = sample_space_dir(save_path)
        manifest_path = sample_dir / "sample_space_manifest.json"
        if not manifest_path.is_file():
            msg = f"Missing {manifest_path}. Run export with sample space."
            raise FileNotFoundError(msg)
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        ch_paths = manifest.get("channel_paths", {})
        vol_key = str(channel)
        if vol_key not in ch_paths:
            msg = f"Channel {channel} not in sample_space manifest."
            raise KeyError(msg)
        vol_path = Path(ch_paths[vol_key])
        labels_path = sample_dir / "division_labels_in_sample_20um.tif"
        if not labels_path.is_file():
            msg = f"Missing {labels_path}"
            raise FileNotFoundError(msg)
        volume = load_registration_volume(vol_path).astype(np.float32, copy=False)
        labels = np.asarray(tifffile.imread(labels_path), dtype=np.int32)
        if volume.shape != labels.shape:
            msg = f"Channel volume {volume.shape} != sample labels {labels.shape}"
            raise ValueError(msg)
        return volume, labels, vol_path, labels_path

    paths = discover_division_viewer_paths(config)
    _, vol_path = _resolve_channel(paths, channel)
    volume = load_registration_volume(vol_path).astype(np.float32, copy=False)
    labels = np.asarray(tifffile.imread(paths.division_labels), dtype=np.int32)
    if volume.shape != labels.shape:
        msg = (
            f"Channel volume {volume.shape} != division labels {labels.shape}. "
            "Rebuild division map or re-export volumes."
        )
        raise ValueError(msg)
    return volume, labels, vol_path, paths.division_labels


def suggest_threshold(volume: np.ndarray, *, percentile: float = 5.0) -> float:
    """Heuristic background threshold when the user does not pass one explicitly."""
    positive = volume[volume > 0]
    if positive.size == 0:
        return 1.0
    return float(max(1.0, np.percentile(positive, percentile) * 2.0))


def run_registration_qc(
    config: BrainPipelineConfig,
    *,
    channel: int = 1,
    threshold: float | None = None,
    inspect: bool = False,
    sweep: bool = False,
    sweep_points: int = 9,
    headless: bool = False,
    force_division_rebuild: bool = False,
    space: str = "atlas",
) -> RegistrationQcResult:
    """Compute (and optionally inspect / sweep) the unassigned division score."""
    space_key = space.strip().lower()
    if force_division_rebuild and space_key == "atlas":
        transform_params = _load_transform_params(config.sample.save_path.expanduser())
        atlas = resolve_brain_atlas_with_config(transform_params.brain_atlas, config.atlas)
        ensure_division_map(atlas, force=True)

    volume, labels, vol_path, _ = load_registration_qc_volumes(
        config, channel=channel, space=space_key
    )
    thr = float(threshold) if threshold is not None else suggest_threshold(volume)

    if inspect and not headless:
        from lightsuite.gui.registration_qc_brain import open_registration_qc_inspector

        open_registration_qc_inspector(volume, labels, threshold=thr)
    elif inspect and headless:
        console.print("[yellow]Skipping Napari inspect in --headless mode.[/yellow]")

    score = compute_unassigned_registration_score(volume, labels, thr)
    sweep_df: pd.DataFrame | None = None
    if sweep:
        sweep_df = threshold_sweep(volume, labels, thr, n_points=sweep_points)

    out_dir = _qc_output_dir(config)
    out_dir.mkdir(parents=True, exist_ok=True)
    stem = f"chan{channel:02d}" + ("_sample" if space_key == "sample" else "")

    meta = {
        "sample": config.sample.name,
        "channel": channel,
        "volume_path": str(vol_path),
    }
    score_csv = out_dir / f"{stem}_unassigned_score.csv"
    score_to_dataframe(score, extra=meta).to_csv(score_csv, index=False)

    sweep_csv: Path | None = None
    sweep_plot: Path | None = None
    if sweep_df is not None:
        sweep_csv = out_dir / f"{stem}_threshold_sweep.csv"
        sweep_df.to_csv(sweep_csv, index=False)
        sweep_plot = write_threshold_sweep_plot(
            sweep_df, out_dir / f"{stem}_threshold_sweep.png"
        )

    summary = {
        "sample": config.sample.name,
        "channel": channel,
        "threshold": thr,
        "score": score,
        "score_csv": str(score_csv),
        "sweep_csv": str(sweep_csv) if sweep_csv else None,
        "sweep_plot": str(sweep_plot) if sweep_plot else None,
    }
    (out_dir / f"{stem}_unassigned_score.json").write_text(
        json.dumps(summary, indent=2),
        encoding="utf-8",
    )

    console.print(
        f"[green]Registration QC[/green] channel {channel}: "
        f"{score['naive_unassigned_percent']:.2f}% unassigned above threshold "
        f"({score['unassigned_image_voxels']}/{score['image_voxels']} voxels @ {thr:.1f})"
    )
    console.print(f"  Saved: {score_csv}")

    return RegistrationQcResult(
        channel=channel,
        threshold=thr,
        score=score,
        sweep=sweep_df,
        score_csv=score_csv,
        sweep_csv=sweep_csv,
        sweep_plot=sweep_plot,
    )
