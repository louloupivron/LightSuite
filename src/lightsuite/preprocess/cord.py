"""Spinal cord preprocessing (prepareCordSampleForRegistration.m port)."""

from __future__ import annotations

import time
from dataclasses import dataclass

import numpy as np
import tifffile
from rich.console import Console
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
from lightsuite.io.cord_reader import read_spinal_cord_sample
from lightsuite.io.cord_volume import (
    cord_downsample_volume,
    normalize_res_um,
    registration_volume_path,
)
from lightsuite.preprocess.cord_checkpoint import CordRegOptsCheckpoint
from lightsuite.registration.cord_paths import cord_cache_dir, cord_save_path

console = Console()


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
            ilast = int(high_idx[-1]) + 1
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


def preprocess_spinal_cord_sample(config: SpinalCordPipelineConfig) -> CordPreprocessResult:
    """Prepare cord sample and atlas for straightening / registration."""
    save_path = cord_save_path(config)
    cache_dir = cord_cache_dir(config)
    save_path.mkdir(parents=True, exist_ok=True)

    console.print(
        f"Loading sample [bold]{config.sample.name}[/bold] from {config.sample.source.path}..."
    )
    sample = read_spinal_cord_sample(config)
    sampleres = normalize_res_um(config.sample.voxel_um)
    regres = normalize_res_um([config.registration.resolution_um] * 3)

    t0 = time.perf_counter()
    finvol = cord_downsample_volume(
        sample.volume,
        sampleres,
        regres,
        native_orisize=sample.native_orisize,
    )
    if finvol.shape != sample.volume.shape:
        console.print(f"Resampled to registration grid in {time.perf_counter() - t0:.2f}s")

    regvolpaths: dict[str, str] = {}
    for ich in range(finvol.shape[3]):
        out_path = registration_volume_path(cache_dir, ich + 1, config.registration.resolution_um)
        if out_path.is_file():
            out_path.unlink()
        tifffile.imwrite(out_path, finvol[:, :, :, ich].astype(np.uint16))
        regvolpaths[str(ich + 1)] = str(out_path)
    console.print(f"Cached {len(regvolpaths)} channel registration TIFF(s)")

    atlas_volumes = load_fiederling_atlas_volumes(config.atlas)
    tv, av = resize_fiederling_atlas(atlas_volumes, config.registration.resolution_um)
    tvpts = extract_atlas_point_cloud(tv, av, atlas_volumes.regions)

    regchan = config.registration.channel_primary
    if sample.n_channels == 1:
        regchan = 1
    if regchan > sample.n_channels:
        msg = f"Registration channel {regchan} out of range (Nchan={sample.n_channels})"
        raise ValueError(msg)

    regvolsize = finvol.shape[:3]
    im = int(np.argmax(regvolsize))
    axperm = [i for i in range(3) if i != im] + [im]
    regvol = np.transpose(finvol[:, :, :, regchan - 1], axperm)
    console.print(f"Long axis at index {im + 1}, permuted to last dimension")

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
    ncheck = max(1, int(np.ceil(0.05 * cordarea.size)))
    frontsum = float(cordarea[:ncheck].mean())
    backsum = float(cordarea[-ncheck:].mean())
    tofliprc = frontsum < backsum
    if tofliprc:
        console.print("Caudorostral direction detected — flipping atlas rostrocaudal axis")
        tv = np.flip(tv, axis=2)
        av = np.flip(av, axis=2)
        tvpts[:, 2] = tv.shape[2] - tvpts[:, 2] + 1

    med = float(np.median(cordarea))
    rstd = _robust_std(cordarea)
    ihigh = cordarea > (med + 3 * rstd)
    console.print(f"{100 * ihigh.mean():.2f}% of slices exceed expected cord area")
    ikeep = _brain_trim_range(ihigh, tofliprc, regvol.shape[2])

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

    native_y, native_x, native_z = sample.native_orisize
    atlas_paths = resolve_fiederling_paths(config.atlas.atlas_dir)
    checkpoint = CordRegOptsCheckpoint(
        sample_name=config.sample.name,
        data_folder=str(config.sample.source.path),
        lsfolder=str(save_path),
        orisize=[native_y, native_x, native_z],
        nchans=sample.n_channels,
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
        tiff_type=sample.layout.value,
        regvolpaths=regvolpaths,
    )
    checkpoint.save(save_path / "regopts.json")
    console.print(f"Wrote checkpoint: {save_path / 'regopts.json'}")
    return CordPreprocessResult(checkpoint=checkpoint)
