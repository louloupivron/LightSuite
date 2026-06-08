"""Refine automatic control-point pairs using AP slice correspondence."""

from __future__ import annotations

from pathlib import Path

import nibabel as nib
import numpy as np
from rich.console import Console

from lightsuite.atlas.registry import resolve_brain_atlas
from lightsuite.config.models import BrainPipelineConfig
from lightsuite.gui.brain_data import load_slice_correspondence, prepare_brain_align_slices_session
from lightsuite.gui.slice_correspondence import SliceCorrespondence, default_correspondence_path
from lightsuite.preprocess.checkpoint import RegOptsCheckpoint
from lightsuite.registration.ap_correspondence import ApFilterStats, filter_pairs_by_ap_correspondence
from lightsuite.registration.volume import resize_atlas_volume

console = Console()


def _atlas_reg_shape(config: BrainPipelineConfig, checkpoint: RegOptsCheckpoint) -> tuple[int, int, int]:
    atlas = resolve_brain_atlas(config.atlas.provider.value, config.atlas.atlas_dir)
    tv = np.asanyarray(nib.load(atlas.template_path).dataobj)
    downfac = float(
        checkpoint.downfac_reg or (config.atlas.resolution_um / checkpoint.registres_um)
    )
    tvreg = resize_atlas_volume(tv.astype(np.float32), downfac, nearest=False)
    return tuple(int(s) for s in tvreg.shape)


def _correspondence_matches_checkpoint(
    correspondence: SliceCorrespondence,
    checkpoint: RegOptsCheckpoint,
) -> bool:
    if checkpoint.original_trans is None:
        return False
    stored = np.asarray(correspondence.original_trans, dtype=float)
    current = np.asarray(checkpoint.original_trans, dtype=float)
    return bool(np.allclose(stored, current, atol=1e-3, rtol=1e-3))


def refine_brain_auto_points(
    config: BrainPipelineConfig,
    *,
    bootstrap_correspondence: bool = False,
    tolerance_vox: float | None = None,
    min_pairs_kept: int | None = None,
    force: bool = False,
) -> Path:
    """Filter init-registration auto pairs using ``slice_correspondence.json``."""
    save_path = config.sample.save_path.expanduser()
    regopts_path = save_path / "regopts.json"
    if not regopts_path.is_file():
        msg = f"Missing {regopts_path}. Run preprocess and init-registration first."
        raise FileNotFoundError(msg)

    checkpoint = RegOptsCheckpoint.load(regopts_path)
    if checkpoint.original_trans is None:
        msg = "regopts.json missing original_trans from init-registration."
        raise RuntimeError(msg)
    if not checkpoint.autocpsample or not checkpoint.autocpatlas:
        msg = "regopts.json has no automatic control-point pairs to refine."
        raise RuntimeError(msg)

    if checkpoint.auto_points_refined and not force:
        console.print(
            "[dim]Auto points already refined "
            f"(mode={checkpoint.auto_points_mode!r}). Use --force to re-run.[/dim]"
        )
        return regopts_path

    correspondence = load_slice_correspondence(save_path)
    if correspondence is None and bootstrap_correspondence:
        console.print("No slice_correspondence.json — bootstrapping via headless align-slices...")
        prepare_brain_align_slices_session(config)
        correspondence = load_slice_correspondence(save_path)

    if correspondence is None or not correspondence.confirmed_anchors():
        msg = (
            "slice_correspondence.json missing or has no confirmed anchors. "
            "Run align-slices or pass --bootstrap-correspondence."
        )
        raise FileNotFoundError(msg)

    if not _correspondence_matches_checkpoint(correspondence, checkpoint):
        console.print(
            "[yellow]Warning:[/yellow] slice_correspondence original_trans differs from "
            "regopts.json — re-run align-slices after init-registration."
        )

    reg = config.registration
    tol = float(tolerance_vox if tolerance_vox is not None else reg.ap_pair_tolerance_vox)
    min_kept = int(min_pairs_kept if min_pairs_kept is not None else reg.ap_pair_min_kept)
    atlas_shape = _atlas_reg_shape(config, checkpoint)

    sample_pts = np.asarray(checkpoint.autocpsample, dtype=float)
    atlas_pts = np.asarray(checkpoint.autocpatlas, dtype=float)
    original_trans = np.asarray(checkpoint.original_trans, dtype=float)

    filtered_sample, filtered_atlas, stats = filter_pairs_by_ap_correspondence(
        sample_pts,
        atlas_pts,
        correspondence,
        original_trans,
        atlas_shape,
        tolerance_vox=tol,
        min_pairs_kept=min_kept,
    )

    checkpoint.autocpsample = filtered_sample.tolist()
    checkpoint.autocpatlas = filtered_atlas.tolist()
    checkpoint.auto_points_source = "global_triage"
    checkpoint.auto_points_refined = True
    checkpoint.auto_points_mode = "ap_filter"
    checkpoint.auto_points_correspondence_path = default_correspondence_path(save_path).name
    checkpoint.save(regopts_path)

    stats.save(save_path / "auto_points_refine_stats.json")
    console.print(
        f"[green]Refined auto points:[/green] {stats.pairs_after}/{stats.pairs_before} pairs kept "
        f"({stats.pairs_removed_ap} removed, tol={stats.tolerance_vox:g} vox, "
        f"median AP residual {stats.median_ap_residual_before_vox:.1f} → "
        f"{stats.median_ap_residual_after_vox:.1f} vox)."
        if stats.median_ap_residual_before_vox is not None
        and stats.median_ap_residual_after_vox is not None
        else f"[green]Refined auto points:[/green] {stats.pairs_after}/{stats.pairs_before} pairs kept."
    )
    return regopts_path
