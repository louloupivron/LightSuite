"""Spinal cord preprocessing (prepareCordSampleForRegistration.m port)."""

from __future__ import annotations

import time
from dataclasses import dataclass

import numpy as np
import tifffile
from skimage.measure import label
from skimage.morphology import dilation

from lightsuite.atlas.fiederling import (
    extract_atlas_point_cloud,
    extract_sample_point_cloud,
    load_fiederling_atlas_volumes,
    resize_fiederling_atlas,
    resolve_fiederling_paths,
)
from lightsuite.config.models import SpinalCordPipelineConfig
from lightsuite.gui.orientation_cord import ensure_cord_orientation
from lightsuite.io.cord_registration_cache import load_or_cache_cord_registration
from lightsuite.io.cord_volume import normalize_res_um
from lightsuite.preprocess.cord_checkpoint import CordRegOptsCheckpoint
from lightsuite.registration.cord_orientation import tofliprc_from_direction
from lightsuite.registration.cord_paths import cord_cache_dir, cord_save_path
from lightsuite.reporter import check_stage_cancelled, emit_pipeline_message, format_duration


@dataclass(frozen=True)
class CordPreprocessResult:
    checkpoint: CordRegOptsCheckpoint


def _robust_std(values: np.ndarray) -> float:
    return float(1.4826 * np.median(np.abs(values - np.median(values))))


def _brain_trim_range(ihigh: np.ndarray, tofliprc: bool, nz: int) -> list[int]:
    """Determine rostrocaudal slice range after removing brain-stem tissue."""
    ikeep = [1, nz]
    if ihigh.mean() <= 0.01:
        return ikeep
    if tofliprc:
        high_idx = np.flatnonzero(ihigh)
        if high_idx.size:
            ilast = int(high_idx[0]) + 1
            ikeep = [1, min(ilast + 1, nz)]
    else:
        normal_idx = np.flatnonzero(~ihigh)
        if normal_idx.size:
            ifirst = int(normal_idx[0]) + 1
            ikeep = [max(ifirst - 1, 1), nz]
    ikeep[0] = max(ikeep[0], 1)
    ikeep[1] = min(ikeep[1], nz)
    if ikeep[1] < ikeep[0]:
        msg = f"Could not determine a valid rostrocaudal slice range (ikeep={ikeep})."
        raise RuntimeError(msg)
    return ikeep


def preprocess_spinal_cord_sample(
    config: SpinalCordPipelineConfig,
    *,
    headless: bool = False,
) -> CordPreprocessResult:
    """Prepare cord sample and atlas for straightening / registration."""
    t0 = time.perf_counter()
    save_path = cord_save_path(config)
    cache_dir = cord_cache_dir(config)
    save_path.mkdir(parents=True, exist_ok=True)

    emit_pipeline_message(f"Preprocess: resolving longitudinal orientation for {config.sample.name}…")
    direction = ensure_cord_orientation(config, headless=headless)
    tofliprc = tofliprc_from_direction(direction)
    emit_pipeline_message(
        f"Preprocess: longitudinal orientation {direction} "
        f"(tofliprc={tofliprc}; from cord_orientation.txt or YAML)"
    )

    check_stage_cancelled()
    emit_pipeline_message(
        f"Preprocess: loading sample from {config.sample.source.path} "
        f"on the {config.registration.resolution_um:g} µm registration grid…"
    )
    t_load = time.perf_counter()
    registration = load_or_cache_cord_registration(config)
    load_elapsed = time.perf_counter() - t_load
    if registration.from_cache:
        emit_pipeline_message(
            f"Preprocess: reused cached registration-grid sample in {format_duration(load_elapsed)}"
        )
    else:
        emit_pipeline_message(
            f"Preprocess: built registration-grid sample in {format_duration(load_elapsed)}"
        )

    finvol = registration.volume
    sampleres = normalize_res_um(config.sample.voxel_um)
    regres = normalize_res_um([config.registration.resolution_um] * 3)
    regvolpaths = registration.regvolpaths
    sample_n_channels = registration.n_channels
    sample_native_orisize = registration.native_orisize
    sample_layout = registration.layout

    check_stage_cancelled()
    emit_pipeline_message("Preprocess: loading and resizing Fiederling atlas…")
    t_atlas = time.perf_counter()
    atlas_volumes = load_fiederling_atlas_volumes(config.atlas)
    tv, av = resize_fiederling_atlas(atlas_volumes, config.registration.resolution_um)
    tvpts = extract_atlas_point_cloud(tv, av, atlas_volumes.regions)
    emit_pipeline_message(
        f"Preprocess: atlas ready on registration grid in {format_duration(time.perf_counter() - t_atlas)}"
    )

    regchan = config.registration.channel_primary
    if sample_n_channels == 1:
        regchan = 1
    if regchan > sample_n_channels:
        msg = f"Registration channel {regchan} out of range (Nchan={sample_n_channels})"
        raise ValueError(msg)

    check_stage_cancelled()
    emit_pipeline_message("Preprocess: segmenting cord and trimming rostrocaudal range…")
    t_seg = time.perf_counter()
    regvolsize = finvol.shape[:3]
    im = int(np.argmax(regvolsize))
    axperm = [i for i in range(3) if i != im] + [im]
    regvol = np.transpose(finvol[:, :, :, regchan - 1], axperm)
    emit_pipeline_message(f"Preprocess: long axis at index {im + 1}, permuted to last dimension")

    rng = np.random.default_rng(1)
    isamprand = rng.choice(regvol.size, size=min(20_000, regvol.size), replace=False)
    samps = regvol.ravel()[isamprand]
    backval = int(np.bincount(samps[samps > 0].astype(np.int64)).argmax()) if np.any(samps > 0) else 0
    regvol_filled = regvol.copy()
    regvol_filled[regvol_filled == 0] = backval
    binvol = regvol_filled > 0
    binvol = dilation(binvol, footprint=np.ones((3, 3, 3), dtype=bool))
    labeled = label(binvol, connectivity=2)
    if labeled.max() == 0:
        msg = "Cord segmentation failed — no foreground voxels."
        raise RuntimeError(msg)
    counts = np.bincount(labeled.ravel())
    counts[0] = 0
    imax = int(np.argmax(counts))
    newvol = labeled == imax
    voxset = np.column_stack(np.where(newvol))
    voxset = voxset[:, [1, 0, 2]]

    cordarea = newvol.sum(axis=(0, 1)).astype(float)
    if tofliprc:
        emit_pipeline_message(
            "Preprocess: caudorostral sample +Z — flipping atlas rostrocaudal axis "
            f"(orientation={direction})"
        )
        tv = np.flip(tv, axis=2)
        av = np.flip(av, axis=2)
        tvpts[:, 2] = tv.shape[2] - tvpts[:, 2] + 1
    else:
        emit_pipeline_message(
            f"Preprocess: rostrocaudal sample +Z — atlas left in native order "
            f"(orientation={direction})"
        )

    med = float(np.median(cordarea))
    rstd = _robust_std(cordarea)
    ihigh = cordarea > (med + 3 * rstd)
    emit_pipeline_message(
        f"Preprocess: {100 * ihigh.mean():.2f}% of slices exceed expected cord area"
    )
    nz = regvol.shape[2]
    ikeep = _brain_trim_range(ihigh, tofliprc, nz)
    n_kept = ikeep[1] - ikeep[0] + 1
    n_removed = nz - n_kept
    if n_removed <= 0:
        emit_pipeline_message(
            f"Preprocess: longitudinal Z-trim keeping all {nz} slices (range {ikeep[0]}–{ikeep[1]})"
        )
    else:
        trim_end = "high-Z (caudal) end" if tofliprc else "low-Z (rostral) end"
        emit_pipeline_message(
            f"Preprocess: longitudinal Z-trim removed {n_removed} slice(s) from the "
            f"{trim_end}; keeping {n_kept}/{nz} slices (range {ikeep[0]}–{ikeep[1]})"
        )

    y0 = max(int(voxset[:, 1].min()) - 1, 0)
    y1 = min(int(voxset[:, 1].max()) + 1, regvol.shape[0] - 1)
    x0 = max(int(voxset[:, 0].min()) - 1, 0)
    x1 = min(int(voxset[:, 0].max()) + 1, regvol.shape[1] - 1)
    yrange = [y0 + 1, y1 + 1]
    xrange = [x0 + 1, x1 + 1]
    revoluse = regvol[y0 : y1 + 1, x0 : x1 + 1, :]
    newvol_crop = newvol[y0 : y1 + 1, x0 : x1 + 1, :]
    regvol_use = revoluse[:, :, ikeep[0] - 1 : ikeep[1]]
    smpts = extract_sample_point_cloud(revoluse, newvol_crop, (ikeep[0], ikeep[1]))
    emit_pipeline_message(
        f"Preprocess: cord segmentation and trim done in {format_duration(time.perf_counter() - t_seg)}"
    )

    check_stage_cancelled()
    emit_pipeline_message("Preprocess: writing cache artifacts and regopts.json…")
    t_write = time.perf_counter()
    regvol_path = cache_dir / "regvol_cord.tif"
    tv_path = cache_dir / "atlas_template.tif"
    av_path = cache_dir / "atlas_annotation.tif"
    smpts_path = cache_dir / "sample_points.npy"
    tvpts_path = cache_dir / "atlas_points.npy"
    tifffile.imwrite(regvol_path, regvol_use.astype(np.uint16))
    tifffile.imwrite(tv_path, tv.astype(np.float32))
    tifffile.imwrite(av_path, av.astype(np.uint16))
    np.save(smpts_path, smpts)
    np.save(tvpts_path, tvpts)

    native_y, native_x, native_z = sample_native_orisize
    atlas_paths = resolve_fiederling_paths(config.atlas.atlas_dir)
    checkpoint = CordRegOptsCheckpoint(
        sample_name=config.sample.name,
        data_folder=str(config.sample.source.path),
        lsfolder=str(save_path),
        orisize=[native_y, native_x, native_z],
        nchans=sample_n_channels,
        sampleres_um=sampleres.tolist(),
        registrationres_um=regres.tolist(),
        reg_channel=regchan,
        sample_perm=[p + 1 for p in axperm],
        tofliprc=tofliprc,
        ikeeprange=ikeep,
        xrange=xrange,
        yrange=yrange,
        regvol_path=str(regvol_path),
        tv_path=str(tv_path),
        av_path=str(av_path),
        smpts_path=str(smpts_path),
        tvpts_path=str(tvpts_path),
        atlas_res_um=list(atlas_volumes.atlas_res_um),
        segments_path=str(atlas_paths.segments_csv),
        regions_path=str(atlas_paths.regions_csv),
        tiff_type=sample_layout.value,
        regvolpaths=regvolpaths,
    )
    checkpoint.save(save_path / "regopts.json")

    from lightsuite.import_.sample_reference import write_sample_reference

    ny, nx, nz = native_y, native_x, native_z
    ref_path = write_sample_reference(
        save_path,
        sample_name=config.sample.name,
        ny=ny,
        nx=nx,
        nz=nz,
        voxel_um=list(config.sample.voxel_um),
    )
    emit_pipeline_message(
        f"Preprocess: wrote regopts.json and {ref_path.name} in "
        f"{format_duration(time.perf_counter() - t_write)}"
    )
    emit_pipeline_message(f"Preprocess: complete in {format_duration(time.perf_counter() - t0)}")

    return CordPreprocessResult(checkpoint=checkpoint)
