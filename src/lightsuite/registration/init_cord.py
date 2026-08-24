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
    describe_cord_z_transinit_source,
    format_cord_z_transinit,
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
from lightsuite.reporter import emit_pipeline_message, format_duration

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

    emit_pipeline_message("Init-registration: loading registration-grid sample and atlas…")
    t0 = time.perf_counter()
    regvol = tifffile.imread(checkpoint.regvol_path).astype(np.uint16)
    tv = tifffile.imread(checkpoint.tv_path).astype(np.float32)
    av = tifffile.imread(checkpoint.av_path).astype(np.uint16)
    samppts = np.load(checkpoint.smpts_path)
    tvpts = np.load(checkpoint.tvpts_path)
    load_elapsed = time.perf_counter() - t0
    emit_pipeline_message(
        "Init-registration: volumes loaded in "
        f"{format_duration(load_elapsed)} "
        f"(sample {regvol.shape}, atlas {tv.shape})"
    )

    nslices = checkpoint.ikeeprange[1] - checkpoint.ikeeprange[0] + 1
    emit_pipeline_message(
        f"Init-registration: computing straightening transforms for {nslices} slices…"
    )
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
    straighten_elapsed = time.perf_counter() - t0
    emit_pipeline_message(
        "Init-registration: straightening done in "
        f"{format_duration(straighten_elapsed)} "
        f"(output {straightvol.shape})"
    )

    centers = np.column_stack([align.fit_x, align.fit_y])
    ikeepori = _remove_outliers(samppts, centers)
    samppts2 = transform_cord_points_slices(samppts[ikeepori], tforms)
    _ = _reduce_points(samppts2, 10_000)
    _ = _reduce_points(tvpts, 10_000)

    correspondence = load_longitudinal_correspondence(save_path)
    transinit = resolve_cord_z_transinit(nslices, tv.shape[2], correspondence)
    if correspondence is None or not correspondence.has_confirmed_anchors(CORD_LONGITUDINAL_AXIS):
        emit_pipeline_message(
            "Init-registration: no confirmed longitudinal_correspondence.json — "
            "using centered z-init (run align-longitudinal for partial-cord samples)"
        )
    else:
        emit_pipeline_message(
            "Init-registration: longitudinal z-init from "
            f"{describe_cord_z_transinit_source(nslices, tv.shape[2], correspondence)} "
            f"({format_cord_z_transinit(transinit)})"
        )

    tvtemp = ndimage.median_filter(tv, size=3)
    emit_pipeline_message("Init-registration: warping atlas with z-scale similarity transform…")
    t0 = time.perf_counter()
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
    similarity_elapsed = time.perf_counter() - t0
    emit_pipeline_message(
        "Init-registration: similarity warp done in "
        f"{format_duration(similarity_elapsed)} "
        f"(atlas {av.shape} → straightvol {straightvol.shape})"
    )

    volmax = float(np.quantile(straightvol, 0.999))
    volplot = np.clip(straightvol.astype(np.float32) / max(volmax, 1.0) * 255.0, 0, 255).astype(np.uint8)
    save_cord_annotation_preview(
        volplot,
        avsim.astype(np.uint16),
        qc_dir / "registration_initial_similarity.png",
    )
    emit_pipeline_message(
        f"Init-registration: wrote QC preview {qc_dir / 'registration_initial_similarity.png'}"
    )

    emit_pipeline_message("Init-registration: running elastix affine registration…")
    t0 = time.perf_counter()
    affine_result = run_affine_registration(
        fixed_volume=straightvol.astype(np.float32),
        moving_volume=atlasuse,
        save_path=save_path,
        spacing_mm=spacing_mm,
        work_dir=cord_work_dir(config, "elastix", "affine"),
        output_path=cord_affine_transform_write_path(config),
    )
    affine_elapsed = time.perf_counter() - t0
    emit_pipeline_message(
        "Init-registration: affine elastix done in "
        f"{format_duration(affine_elapsed)} "
        f"({affine_result.copied_transform_path.name})"
    )

    emit_pipeline_message("Init-registration: applying affine warp to atlas annotation…")
    t0 = time.perf_counter()
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
    affine_warp_elapsed = time.perf_counter() - t0
    emit_pipeline_message(
        "Init-registration: affine annotation warp done in "
        f"{format_duration(affine_warp_elapsed)}"
    )
    save_cord_annotation_preview(
        volplot,
        avshow.astype(np.uint16),
        qc_dir / "registration_initial_affine.png",
    )
    emit_pipeline_message(
        f"Init-registration: wrote QC preview {qc_dir / 'registration_initial_affine.png'}"
    )

    straightvol_path = cache_dir / "straightvol.tif"
    slicetforms_path = cache_dir / "slicetforms.npy"
    tifffile.imwrite(straightvol_path, straightvol)
    save_slicetforms(slicetforms_path, tforms)

    checkpoint.straightvol_path = str(straightvol_path)
    checkpoint.slicetforms_path = str(slicetforms_path)
    checkpoint.affine_atlas_to_samp = transaff.tolist()
    checkpoint.save(regopts_path)
    emit_pipeline_message(
        "Init-registration: wrote straightvol, slicetforms, and affine transform; "
        f"updated {regopts_path}"
    )
    return checkpoint
