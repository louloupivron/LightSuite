"""Parcellation intensity statistics (generateRegisteredBrainVolumes.m)."""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd

from lightsuite.atlas.registry import AtlasPaths


@dataclass
class ParcellationResult:
    area_ids: np.ndarray
    median_over_areas: np.ndarray  # (n_areas, 2)
    std_over_areas: np.ndarray
    volume_over_areas: np.ndarray


def _std_per_group(values: np.ndarray) -> float:
    if values.size < 2:
        return 0.0
    return float(np.std(values.astype(np.float64), ddof=0))


def _accumulate_side(
    labels: np.ndarray,
    values: np.ndarray,
    area_ids: np.ndarray,
    voxel_mm3: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    n_bins = int(max(labels.max(initial=0), area_ids.max(initial=0))) + 2
    medians = np.zeros(n_bins, dtype=np.float32)
    stds = np.zeros(n_bins, dtype=np.float32)
    volumes = np.zeros(n_bins, dtype=np.float32)

    positive = labels > 0
    if np.any(positive):
        lab = labels[positive].astype(np.int64)
        vals = values[positive].astype(np.float64)
        for uid in np.unique(lab):
            mask = lab == uid
            group = vals[mask]
            bin_idx = int(uid) + 1
            medians[bin_idx] = np.median(group)
            stds[bin_idx] = _std_per_group(group)
            volumes[bin_idx] = float(mask.sum()) * voxel_mm3

    background = labels <= 0
    if np.any(background):
        bg = values[background].astype(np.float64)
        medians[0] = np.median(bg)
        stds[0] = _std_per_group(bg)

    lookup = np.zeros(len(area_ids), dtype=np.float32)
    lookup_std = np.zeros(len(area_ids), dtype=np.float32)
    lookup_vol = np.zeros(len(area_ids), dtype=np.float32)
    for row, aid in enumerate(area_ids.astype(np.int64)):
        if aid == 0:
            lookup[row] = medians[0]
            lookup_std[row] = stds[0]
            lookup_vol[row] = volumes[0]
        else:
            bin_idx = int(aid) + 1
            lookup[row] = medians[bin_idx]
            lookup_std[row] = stds[bin_idx]
            lookup_vol[row] = volumes[bin_idx]
    return lookup, lookup_std, lookup_vol


def compute_allen_parcellation(
    registered_volume: np.ndarray,
    atlas: AtlasPaths,
    atlas_resolution_um: float,
    *,
    parcellation_csv: Path | None = None,
) -> ParcellationResult:
    csv_path = parcellation_csv or _resolve_allen_parcellation_csv()
    if csv_path is None or not csv_path.is_file():
        msg = (
            "Allen parcellation CSV not found. Set LIGHTSUITE_ALLEN_PARCELLATION_CSV "
            "or place parcellation_to_parcellation_term_membership.csv on the atlas path."
        )
        raise FileNotFoundError(msg)

    parcelinfo = pd.read_csv(csv_path)
    substr = parcelinfo["parcellation_term_set_name"] == "substructure"
    area_ids = parcelinfo.loc[substr, "parcellation_index"].drop_duplicates().to_numpy(dtype=np.int64)

    av = np.asanyarray(nib.load(atlas.annotation_path).dataobj)
    n_z_half = av.shape[2] // 2
    voxel_mm3 = (atlas_resolution_um * 1e-3) ** 3

    median_over = np.full((len(area_ids), 2), np.nan, dtype=np.float32)
    std_over = np.full((len(area_ids), 2), np.nan, dtype=np.float32)
    volume_over = np.full((len(area_ids), 2), np.nan, dtype=np.float32)

    for side, z_slice in enumerate([slice(0, n_z_half), slice(n_z_half, av.shape[2])]):
        labels = av[:, :, z_slice].reshape(-1)
        values = registered_volume[:, :, z_slice].reshape(-1)
        med, std, vol = _accumulate_side(labels, values, area_ids, voxel_mm3)
        median_over[:, side] = med
        std_over[:, side] = std
        volume_over[:, side] = vol

    return ParcellationResult(
        area_ids=area_ids,
        median_over_areas=median_over,
        std_over_areas=std_over,
        volume_over_areas=volume_over,
    )


def compute_perens_parcellation(
    registered_volume: np.ndarray,
    atlas: AtlasPaths,
    atlas_resolution_um: float,
    *,
    ml_axis: int = 2,
) -> ParcellationResult:
    if atlas.structures_csv_path is None:
        msg = "Perens structures CSV not found beside atlas NIfTIs."
        raise FileNotFoundError(msg)

    st_tab = pd.read_csv(atlas.structures_csv_path)
    area_ids = st_tab["id"].to_numpy(dtype=np.int64)
    av = np.asanyarray(nib.load(atlas.annotation_path).dataobj)
    if tuple(av.shape) != tuple(registered_volume.shape):
        msg = f"Annotation shape {av.shape} != registered volume {registered_volume.shape}"
        raise ValueError(msg)

    axis = ml_axis - 1
    coords = np.indices(av.shape)[axis].astype(np.float32)
    brain = av > 0
    split_plane = round(float(coords[brain].mean()))
    side_lower = (coords <= split_plane) & brain
    side_upper = (coords > split_plane) & brain

    voxel_mm3 = (atlas_resolution_um * 1e-3) ** 3
    median_over = np.full((len(area_ids), 2), np.nan, dtype=np.float32)
    std_over = np.full((len(area_ids), 2), np.nan, dtype=np.float32)
    volume_over = np.full((len(area_ids), 2), np.nan, dtype=np.float32)

    for side, side_mask in enumerate([side_lower, side_upper]):
        labels = av[side_mask]
        values = registered_volume[side_mask]
        positive = labels > 0
        labels = labels[positive]
        values = values[positive]
        if labels.size == 0:
            continue
        for uid in np.unique(labels):
            row = np.where(area_ids == uid)[0]
            if row.size == 0:
                continue
            mask = labels == uid
            median_over[row[0], side] = np.median(values[mask])
            std_over[row[0], side] = _std_per_group(values[mask])
            volume_over[row[0], side] = float(mask.sum()) * voxel_mm3

    return ParcellationResult(
        area_ids=area_ids,
        median_over_areas=median_over,
        std_over_areas=std_over,
        volume_over_areas=volume_over,
    )


def write_parcellation_csv(path: Path, result: ParcellationResult) -> None:
    df = pd.DataFrame(
        {
            "parcellation_index": result.area_ids,
            "RightSideIntensity": result.median_over_areas[:, 0],
            "LeftSideIntensity": result.median_over_areas[:, 1],
            "RightSideIntensityStd": result.std_over_areas[:, 0],
            "LeftSideIntensityStd": result.std_over_areas[:, 1],
            "RightSideVolume[mm3]": result.volume_over_areas[:, 0],
            "LeftSideVolume[mm3]": result.volume_over_areas[:, 1],
        }
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)


def _resolve_allen_parcellation_csv() -> Path | None:
    env = os.environ.get("LIGHTSUITE_ALLEN_PARCELLATION_CSV", "").strip()
    if env:
        return Path(env).expanduser()
    for candidate in Path.cwd().glob("**/parcellation_to_parcellation_term_membership.csv"):
        return candidate
    return None
