"""Initialize brain registration (initializeRegistration.m port)."""

from __future__ import annotations

import time
from pathlib import Path

import nibabel as nib
import numpy as np
from rich.console import Console

from lightsuite.atlas.registry import resolve_brain_atlas
from lightsuite.config.models import BrainPipelineConfig
from lightsuite.preprocess.checkpoint import RegOptsCheckpoint
from lightsuite.registration.align import (
    downsample_point_cloud,
    estimate_similarity_transform,
    similarity_scale,
    triage_and_match_clouds,
)
from lightsuite.registration.orientation import orientation_path, resolve_orientation, save_orientation
from lightsuite.registration.plots import save_initial_registration_previews
from lightsuite.registration.points import extract_atlas_points_gradient, extract_sample_points
from lightsuite.registration.volume import (
    load_registration_volume,
    normalize_registration_volume,
    permute_brain_volume,
    resize_atlas_volume,
)

console = Console()


def initialize_brain_registration(config: BrainPipelineConfig) -> RegOptsCheckpoint:
    """Coarse-align sample to atlas and update regopts checkpoint."""
    save_path = config.sample.save_path.expanduser()
    regopts_path = save_path / "regopts.json"
    if not regopts_path.is_file():
        msg = f"Missing checkpoint {regopts_path}. Run 'lightsuite brain preprocess' first."
        raise FileNotFoundError(msg)

    checkpoint = RegOptsCheckpoint.load(regopts_path)
    backvol = load_registration_volume(Path(checkpoint.regvolpath))
    downfac = config.atlas.resolution_um / checkpoint.registres_um

    atlas = resolve_brain_atlas(config.atlas.provider.value, config.atlas.atlas_dir)
    tv = np.asanyarray(nib.load(atlas.template_path).dataobj)
    av = np.asanyarray(nib.load(atlas.annotation_path).dataobj)
    tvreg = resize_atlas_volume(tv.astype(np.float32), downfac, nearest=False)
    avreg = resize_atlas_volume(av.astype(np.float32), downfac, nearest=True)

    if atlas.boundary_path is not None:
        boundary_full = np.asanyarray(nib.load(atlas.boundary_path).dataobj)
        console.print(f"Using atlas boundary volume [bold]{atlas.boundary_path.name}[/bold]")
    else:
        console.print(
            "[yellow]annotation_boundary_10.nii.gz not found in atlas_dir;[/yellow] "
            "deriving boundaries from annotation labels."
        )
        from lightsuite.registration.plots import boundary_volume_from_annotation

        boundary_full = boundary_volume_from_annotation(av)
    boundary_reg = resize_atlas_volume(boundary_full.astype(np.float32), downfac, nearest=True)
    boundary_reg = ((boundary_reg > 0).astype(np.uint8)) * 255

    permvec = resolve_orientation(config, save_path)
    orient_file = orientation_path(save_path)
    if config.registration.orientation is not None:
        save_orientation(save_path, permvec)
    elif orient_file.is_file():
        console.print(f"Loaded brain orientation {permvec} from {orient_file}")
    else:
        console.print(
            "[yellow]No brain_orientation.txt found;[/yellow] using default [1, 2, 3]. "
            "Run 'lightsuite brain check-orientation' to verify axes interactively."
        )

    newvol = normalize_registration_volume(backvol)
    console.print("Creating cloud for sample volume...", end=" ")
    t0 = time.perf_counter()
    volumereg = permute_brain_volume(newvol, permvec)
    ls_cloud = extract_sample_points(volumereg, config.registration.cloud_threshold)
    console.print(f"Done in {time.perf_counter() - t0:.1f}s. Found {ls_cloud.shape[0]} points.")

    console.print("Loading atlas and generating atlas cloud...", end=" ")
    t0 = time.perf_counter()
    tv_for_points = tvreg.copy()
    tv_for_points[avreg == 0] = 0
    tv_cloud = extract_atlas_points_gradient(tv_for_points, avreg, sigma=20.0, threshold=5.0)
    console.print(f"Done in {time.perf_counter() - t0:.1f}s. Found {tv_cloud.shape[0]} points.")

    if ls_cloud.shape[0] < 10 or tv_cloud.shape[0] < 10:
        msg = "Too few points extracted for coarse registration."
        raise RuntimeError(msg)

    if tv_cloud.shape[0] > 100_000:
        before = tv_cloud.shape[0]
        tv_cloud = downsample_point_cloud(tv_cloud, 100_000)
        console.print(
            f"Downsampled atlas cloud from {before:,} to {tv_cloud.shape[0]:,} points for alignment."
        )

    console.print("Estimating initial similarity transform...", end=" ")
    t0 = time.perf_counter()
    transform_icp, transform_matlab = estimate_similarity_transform(tv_cloud, ls_cloud)
    console.print(
        f"Done in {time.perf_counter() - t0:.1f}s "
        f"(scale={similarity_scale(transform_icp):.3f})."
    )

    console.print("Identifying candidate corresponding points...", end=" ")
    t0 = time.perf_counter()
    cpsample, cpatlas = triage_and_match_clouds(ls_cloud, tv_cloud, transform_icp)
    console.print(f"Done in {time.perf_counter() - t0:.1f}s. Pairs: {cpsample.shape[0]}")

    console.print("Saving initial registration previews...", end=" ")
    t0 = time.perf_counter()
    warped_boundary = save_initial_registration_previews(
        save_path, volumereg, boundary_reg, transform_matlab
    )
    console.print(
        f"Done in {time.perf_counter() - t0:.1f}s "
        f"({warped_boundary:,} warped boundary voxels)."
    )
    if warped_boundary == 0:
        console.print(
            "[yellow]Warning:[/yellow] no atlas boundaries overlapped the sample after warping. "
            "Check coarse alignment and ensure annotation_boundary_10.nii.gz is in atlas_dir."
        )

    checkpoint.permute_sample_to_atlas = permvec
    checkpoint.original_trans = transform_matlab.tolist()
    checkpoint.downfac_reg = downfac
    checkpoint.autocpsample = cpsample.tolist()
    checkpoint.autocpatlas = cpatlas.tolist()
    checkpoint.brain_atlas = config.atlas.provider.value
    checkpoint.save(regopts_path)
    console.print(f"Updated checkpoint [bold]{regopts_path}[/bold]")
    return checkpoint
