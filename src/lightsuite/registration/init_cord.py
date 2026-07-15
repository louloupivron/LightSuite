"""Initialize spinal cord registration (initializeCordRegistration.m port)."""

from __future__ import annotations

import time

import numpy as np
import tifffile
from rich.console import Console
from scipy import ndimage

from lightsuite.config.models import SpinalCordPipelineConfig
from lightsuite.preprocess.cord_checkpoint import CordRegOptsCheckpoint, SpinalAlignmentCheckpoint
from lightsuite.registration.cord_affine import (
    fit_cord_affine_atlas_to_straightvol,
    warp_cord_atlas_to_straightvol,
)
from lightsuite.registration.cord_longitudinal import (
    CORD_LONGITUDINAL_AXIS,
    load_longitudinal_correspondence,
    resolve_cord_z_transinit,
)
from lightsuite.registration.cord_paths import (
    cord_affine_transform_write_path,
    cord_cache_dir,
    cord_qc_dir,
    cord_save_path,
    cord_work_dir,
)
from lightsuite.registration.cord_plots import save_cord_annotation_preview
from lightsuite.registration.elastix.affine import run_affine_registration
from lightsuite.registration.straightening import (
    compute_straightening_transforms,
    save_slicetforms,
    transform_cord_images_slices,
    transform_cord_points_slices,
)
from lightsuite.registration.warp import warp_volume_affine

console = Console()


def _reduce_points(pts: np.ndarray, n_target: int) -> np.ndarray:
    """Port of reducePoints local function."""
    nlevels = int(pts[:, 2].max())
    downfac = n_target / max(nlevels, 1)
    chunks: list[np.ndarray] = []
    for ii in range(1, nlevels + 1):
        icurr = pts[:, 2] == ii
        if not np.any(icurr):
            continue
        n_down = max(6, int(np.floor(np.count_nonzero(icurr) / downfac)))
        subset = pts[icurr]
        if subset.shape[0] > n_down:
            idx = np.linspace(0, subset.shape[0] - 1, n_down, dtype=int)
            subset = subset[idx]
        chunks.append(subset)
    if not chunks:
        return np.zeros((0, 3), dtype=np.float32)
    return np.vstack(chunks).astype(np.float32)


def _remove_outliers(
    pts: np.ndarray,
    centers: np.ndarray,
    thresuse: float = 2.0,
) -> np.ndarray:
    """Port of removeOutliers local function."""
    ikeep = np.ones(pts.shape[0], dtype=bool)
    allds = np.linalg.norm(pts[:, :2] - centers[pts[:, 2].astype(int) - 1], axis=1)
    nmax = int(pts[:, 2].max())
    dsperslice = np.array(
        [np.median(np.abs(d - np.median(d))) * 1.4826 if d.size else 1.0 for d in _split_by_slice(allds, pts[:, 2], nmax)]
    )
    medperslice = np.array([np.median(d) if d.size else 0.0 for d in _split_by_slice(allds, pts[:, 2], nmax)])
    medperslice = ndimage.median_filter(medperslice, size=50, mode="nearest")
    dsperslice = ndimage.median_filter(dsperslice, size=50, mode="nearest")
    z_idx = pts[:, 2].astype(int) - 1
    inoise = (allds - medperslice[z_idx]) / np.maximum(dsperslice[z_idx], 1e-6) > thresuse
    if np.any(inoise):
        console.print(f"Removed {int(inoise.sum())} outlier points around the cord")
        ikeep = ~inoise
    return ikeep


def _split_by_slice(values: np.ndarray, slices: np.ndarray, nmax: int) -> list[np.ndarray]:
    out: list[np.ndarray] = []
    for ii in range(1, nmax + 1):
        out.append(values[slices == ii])
    return out


def initialize_cord_registration(config: SpinalCordPipelineConfig) -> CordRegOptsCheckpoint:
    """Apply straightening, z-scale, and elastix affine registration."""
    save_path = cord_save_path(config)
    cache_dir = cord_cache_dir(config)
    qc_dir = cord_qc_dir(config)
    regopts_path = save_path / "regopts.json"
    align_path = save_path / "spinal_alignment_opt.json"
    if not regopts_path.is_file():
        msg = f"Missing {regopts_path}. Run 'lightsuite spinal preprocess' first."
        raise FileNotFoundError(msg)
    if not align_path.is_file():
        msg = f"Missing {align_path}. Run 'lightsuite spinal straighten' first."
        raise FileNotFoundError(msg)

    checkpoint = CordRegOptsCheckpoint.load(regopts_path)
    align = SpinalAlignmentCheckpoint.load(align_path)
    regvol = tifffile.imread(checkpoint.regvol_path).astype(np.uint16)
    tv = tifffile.imread(checkpoint.tv_path).astype(np.float32)
    av = tifffile.imread(checkpoint.av_path).astype(np.uint16)
    samppts = np.load(checkpoint.smpts_path)
    tvpts = np.load(checkpoint.tvpts_path)

    nslices = checkpoint.ikeeprange[1] - checkpoint.ikeeprange[0] + 1
    target_center = (0.5 * tv.shape[1], 0.5 * tv.shape[0])
    tforms = compute_straightening_transforms(
        np.array(align.fit_x),
        np.array(align.fit_y),
        np.array(align.fit_theta),
        target_center,
        config.registration.target_orientation_deg,
    )

    t0 = time.perf_counter()
    sizetv = (tv.shape[0], tv.shape[1])
    straightvol = transform_cord_images_slices(regvol, tforms, sizetv)
    console.print(f"Straightening transforms applied in {time.perf_counter() - t0:.2f}s")

    centers = np.column_stack([align.fit_x, align.fit_y])
    ikeepori = _remove_outliers(samppts, centers)
    samppts2 = transform_cord_points_slices(samppts[ikeepori], tforms)
    _ = _reduce_points(samppts2, 10_000)
    _ = _reduce_points(tvpts, 10_000)

    correspondence = load_longitudinal_correspondence(save_path)
    if correspondence is None or not correspondence.has_confirmed_anchors(CORD_LONGITUDINAL_AXIS):
        console.print(
            "[yellow]No confirmed longitudinal_correspondence.json — using centered z-init. "
            "Run 'lightsuite spinal align-longitudinal' for partial-cord samples.[/yellow]"
        )
    else:
        console.print(
            f"Using longitudinal correspondence "
            f"({len(correspondence.confirmed_anchors(CORD_LONGITUDINAL_AXIS))} confirmed anchors)"
        )

    transinit = resolve_cord_z_transinit(nslices, tv.shape[2], correspondence)

    tvtemp = ndimage.median_filter(tv, size=3)
    atlasuse = warp_volume_affine(
        tvtemp,
        transinit,
        straightvol.shape,
        order=1,
        point_coords="array",
    )
    spacing_mm = config.registration.resolution_um * 1e-3
    avsim = warp_cord_atlas_to_straightvol(
        av,
        transinit=transinit,
        elastix_affine_path=None,
        output_shape=straightvol.shape,
        spacing_mm=spacing_mm,
        work_dir=cord_work_dir(config, "transformix", "init", "similarity"),
        nearest=True,
    )
    console.print(
        f"QC volume shapes: straightvol {straightvol.shape}, avsim {avsim.shape} "
        f"(atlas {av.shape} z-scaled to sample grid)"
    )

    volmax = float(np.quantile(straightvol, 0.999))
    volplot = np.clip(straightvol.astype(np.float32) / max(volmax, 1.0) * 255.0, 0, 255).astype(np.uint8)
    save_cord_annotation_preview(
        volplot,
        avsim.astype(np.uint16),
        qc_dir / "registration_initial_similarity.png",
    )

    affine_result = run_affine_registration(
        fixed_volume=straightvol.astype(np.float32),
        moving_volume=atlasuse,
        save_path=save_path,
        spacing_mm=spacing_mm,
        work_dir=cord_work_dir(config, "elastix", "affine"),
        output_path=cord_affine_transform_write_path(config),
    )
    avshow = warp_cord_atlas_to_straightvol(
        av,
        transinit=transinit,
        elastix_affine_path=affine_result.copied_transform_path,
        output_shape=straightvol.shape,
        spacing_mm=spacing_mm,
        work_dir=cord_work_dir(config, "transformix", "init", "affine"),
        nearest=True,
    )
    transaff = fit_cord_affine_atlas_to_straightvol(
        tv.shape,
        transinit,
        affine_result.copied_transform_path,
        spacing_mm,
    )
    console.print(
        f"QC volume shapes after affine: avshow {avshow.shape} "
        f"(matches straightvol {straightvol.shape})"
    )
    save_cord_annotation_preview(
        volplot,
        avshow.astype(np.uint16),
        qc_dir / "registration_initial_affine.png",
    )

    straightvol_path = cache_dir / "straightvol.tif"
    slicetforms_path = cache_dir / "slicetforms.npy"
    tifffile.imwrite(straightvol_path, straightvol)
    save_slicetforms(slicetforms_path, tforms)

    checkpoint.straightvol_path = str(straightvol_path)
    checkpoint.slicetforms_path = str(slicetforms_path)
    checkpoint.affine_atlas_to_samp = transaff.tolist()
    checkpoint.save(regopts_path)
    console.print(f"Updated checkpoint: {regopts_path}")
    return checkpoint
