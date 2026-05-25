"""Initialize brain registration (initializeRegistration.m port)."""

from __future__ import annotations

import time
from pathlib import Path

import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np
from rich.console import Console

from lightsuite.atlas.registry import resolve_brain_atlas
from lightsuite.config.models import BrainPipelineConfig
from lightsuite.preprocess.checkpoint import RegOptsCheckpoint
from lightsuite.registration.align import estimate_similarity_transform, triage_and_match_clouds
from lightsuite.registration.points import extract_atlas_points_gradient, extract_sample_points
from lightsuite.registration.volume import (
    load_registration_volume,
    normalize_registration_volume,
    permute_brain_volume,
    resize_atlas_volume,
)

console = Console()


def _orientation_path(save_path: Path) -> Path:
    return save_path / "brain_orientation.txt"


def _load_orientation(
    config: BrainPipelineConfig,
    save_path: Path,
) -> list[int]:
    if config.registration.orientation is not None:
        permvec = list(config.registration.orientation)
    else:
        orient_file = _orientation_path(save_path)
        if orient_file.is_file():
            values = np.loadtxt(orient_file, dtype=np.int64)
            permvec = values.flatten().astype(int).tolist()
            console.print(f"Loaded brain orientation {permvec} from {orient_file}")
        else:
            permvec = [1, 2, 3]
            console.print(
                "[yellow]No brain_orientation.txt found;[/yellow] using default [1, 2, 3]. "
                "Set registration.orientation in config or add brain_orientation.txt."
            )
    if len(set(abs(v) for v in permvec)) != 3:
        msg = f"Invalid orientation permutation: {permvec}"
        raise ValueError(msg)
    return permvec


def _save_orientation_preview(
    save_path: Path,
    sample: np.ndarray,
    permvec: list[int],
) -> None:
    """Save midslice comparison PNGs (dim{1,2,3}_initial_registration.png)."""
    for idim in range(3):
        fig, ax = plt.subplots(figsize=(6, 6))
        mid = int(sample.shape[idim] / 2)
        if idim == 0:
            sl_sample = sample[mid, :, :]
        elif idim == 1:
            sl_sample = sample[:, mid, :]
        else:
            sl_sample = sample[:, :, mid]
        ax.imshow(sl_sample, cmap="gray", aspect="auto")
        ax.set_title(f"Sample dim {idim + 1} (perm={permvec})")
        ax.axis("off")
        out = save_path / f"dim{idim + 1}_initial_registration.png"
        fig.savefig(out, dpi=120, bbox_inches="tight")
        plt.close(fig)


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

    permvec = _load_orientation(config, save_path)
    if config.registration.orientation is not None:
        np.savetxt(_orientation_path(save_path), np.array(permvec, dtype=int), fmt="%d")

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

    console.print("Estimating initial similarity transform...", end=" ")
    t0 = time.perf_counter()
    transform = estimate_similarity_transform(tv_cloud, ls_cloud)
    console.print(f"Done in {time.perf_counter() - t0:.1f}s.")

    console.print("Identifying candidate corresponding points...", end=" ")
    t0 = time.perf_counter()
    cpsample, cpatlas = triage_and_match_clouds(ls_cloud, tv_cloud, transform)
    console.print(f"Done in {time.perf_counter() - t0:.1f}s. Pairs: {cpsample.shape[0]}")

    _save_orientation_preview(save_path, volumereg, permvec)

    checkpoint.permute_sample_to_atlas = permvec
    checkpoint.original_trans = transform.tolist()
    checkpoint.downfac_reg = downfac
    checkpoint.autocpsample = cpsample.tolist()
    checkpoint.autocpatlas = cpatlas.tolist()
    checkpoint.brain_atlas = config.atlas.provider.value
    checkpoint.save(regopts_path)
    console.print(f"Updated checkpoint [bold]{regopts_path}[/bold]")
    return checkpoint
