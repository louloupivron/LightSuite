"""Fiederling et al. 2021 spinal cord atlas provider."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile
from scipy import ndimage
from skimage.transform import resize

from lightsuite.config.models import CordAtlasConfig


@dataclass(frozen=True)
class FiederlingAtlasPaths:
    atlas_dir: Path
    template_path: Path
    annotation_path: Path
    hemisphere_path: Path
    segments_csv: Path
    regions_csv: Path


@dataclass
class FiederlingAtlasVolumes:
    template: np.ndarray
    annotation: np.ndarray
    atlas_res_um: tuple[float, float, float]
    segments: pd.DataFrame
    regions: pd.DataFrame


def resolve_fiederling_paths(atlas_dir: Path) -> FiederlingAtlasPaths:
    """Resolve required Fiederling atlas files under ``atlas_dir``."""
    root = atlas_dir.expanduser().resolve()
    paths = FiederlingAtlasPaths(
        atlas_dir=root,
        template_path=root / "Template.tif",
        annotation_path=root / "Annotation.tif",
        hemisphere_path=root / "Hemisphere_Annotation.tif",
        segments_csv=root / "Segments.csv",
        regions_csv=root / "Atlas_Regions.csv",
    )
    missing = [
        name
        for name, path in (
            ("Template.tif", paths.template_path),
            ("Annotation.tif", paths.annotation_path),
            ("Segments.csv", paths.segments_csv),
            ("Atlas_Regions.csv", paths.regions_csv),
        )
        if not path.is_file()
    ]
    if missing:
        msg = f"Missing Fiederling atlas files in {root}: {', '.join(missing)}"
        raise FileNotFoundError(msg)
    return paths


def load_fiederling_atlas_volumes(atlas: CordAtlasConfig) -> FiederlingAtlasVolumes:
    """Load template + annotation TIFF stacks and CSV tables."""
    paths = resolve_fiederling_paths(atlas.atlas_dir)
    template = tifffile.imread(paths.template_path)
    annotation = tifffile.imread(paths.annotation_path)
    if template.ndim != 3 or annotation.ndim != 3:
        msg = "Fiederling atlas volumes must be 3D TIFF stacks."
        raise ValueError(msg)
    segments = pd.read_csv(paths.segments_csv)
    regions = pd.read_csv(paths.regions_csv)
    return FiederlingAtlasVolumes(
        template=template.astype(np.float32),
        annotation=annotation.astype(np.uint16),
        # Native SC_P56 atlas is (rostrocaudal, transverse, transverse) = (1567, 253, 322).
        # The rostrocaudal axis (dim0) is the coarse 20 µm slice axis; the transverse
        # plane (dims 1-2) is imaged at 10 µm. Resolution must match this axis order.
        atlas_res_um=(20.0, 10.0, 10.0),
        segments=segments,
        regions=regions,
    )


def upsample_to_fiederling_native(
    volume: np.ndarray,
    template_native: np.ndarray,
) -> np.ndarray:
    """Upsample a registered volume from the 20 µm grid to native 10×10×20 µm template shape.

    Registration volumes use the transverse plane in the first two dimensions and the
    rostrocaudal length axis last (Y, X, Z). The native Fiederling template stores the
    length axis first (Z, Y, X).
    """
    if volume.ndim != 3:
        msg = f"Expected 3D registration volume, got shape {volume.shape}"
        raise ValueError(msg)
    target = template_native.shape
    vol_native_axes = np.transpose(volume, (2, 0, 1))
    if vol_native_axes.shape == target:
        return vol_native_axes.astype(np.uint16)
    return resize(
        vol_native_axes,
        target,
        order=1,
        preserve_range=True,
        anti_aliasing=True,
    ).astype(np.uint16)


def resize_fiederling_atlas(
    volumes: FiederlingAtlasVolumes,
    output_res_um: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Resample atlas to registration resolution and orient length axis last.

    Port of loadSpinalCordAtlasAndPoints.m. The native atlas stores the
    rostrocaudal (length) axis first; MATLAB regopts.tv keeps the transverse
    plane in the first two dimensions and the length axis last (e.g. 127x161x1567).
    The straightening and registration code assumes this convention, so the
    resampled atlas must be permuted to match.
    """
    native = np.array(volumes.atlas_res_um, dtype=float)
    target = np.array([output_res_um, output_res_um, output_res_um], dtype=float)
    scale = native / target

    def _resample(volume: np.ndarray, order: int, anti_aliasing: bool) -> np.ndarray:
        new_shape = [int(round(s * sc)) for s, sc in zip(volume.shape, scale, strict=True)]
        return resize(
            volume,
            new_shape,
            order=order,
            preserve_range=True,
            anti_aliasing=anti_aliasing,
        )

    tv = _resample(volumes.template, 1, True).astype(np.float32)
    av = _resample(volumes.annotation, 0, False).astype(np.uint16)

    # Move the rostrocaudal (longest) axis to the last dimension so it becomes the
    # slice/length axis, matching the sample volume convention and MATLAB regopts.tv.
    long_axis = int(np.argmax(tv.shape))
    if long_axis != 2:
        perm = [a for a in range(3) if a != long_axis] + [long_axis]
        tv = np.ascontiguousarray(np.transpose(tv, perm))
        av = np.ascontiguousarray(np.transpose(av, perm))
    return tv, av


def _white_matter_ids(regions: pd.DataFrame) -> np.ndarray:
    name_col = regions["name"] if "name" in regions.columns else regions.iloc[:, 1]
    children_col = regions["children_IDs"] if "children_IDs" in regions.columns else None
    if children_col is None:
        return np.array([], dtype=int)
    lower = name_col.astype(str).str.lower()
    iwm = lower.str.contains("fasciculus") | lower.str.contains("funiculus")
    wmchildren: list[str] = []
    for raw in children_col[iwm]:
        if pd.isna(raw):
            continue
        wmchildren.extend(str(raw).split(","))
    ids = pd.to_numeric(wmchildren, errors="coerce")
    return ids[~np.isnan(ids)].astype(int)


def extract_atlas_point_cloud(
    template: np.ndarray,
    annotation: np.ndarray,
    regions: pd.DataFrame,
    *,
    downsample_grid: float = 6.0,
) -> np.ndarray:
    """Port of loadSpinalCordAtlasAndPoints.m atlas point extraction."""
    wm_ids = _white_matter_ids(regions)
    ptlist: list[np.ndarray] = []
    for islice in range(annotation.shape[2]):
        avcurr = annotation[:, :, islice]
        iswm = np.isin(avcurr, wm_ids)
        iswmtrans = ndimage.sobel(iswm.astype(float)) > 0
        isextpt = ndimage.sobel((avcurr > 0).astype(float)) > 0
        ipts = np.flatnonzero(iswmtrans | isextpt)
        if ipts.size == 0:
            continue
        rr, cc = np.unravel_index(ipts, avcurr.shape)
        pts = np.column_stack([cc + 1.0, rr + 1.0, np.full(rr.shape, islice + 1.0)])
        ptlist.append(_downsample_points_per_slice(pts, downsample_grid))
    if not ptlist:
        return np.zeros((0, 3), dtype=np.float32)
    return np.concatenate(ptlist, axis=0).astype(np.float32)


def extract_sample_point_cloud(
    regvol: np.ndarray,
    volgood: np.ndarray,
    ikeep: tuple[int, int],
    *,
    downsample_grid: float = 12.0,
) -> np.ndarray:
    """Port of spinalCordPointCloud.m."""
    ptlist: list[np.ndarray] = []
    thresuse = 1.0
    for islice in range(regvol.shape[2]):
        scurr = ndimage.median_filter(regvol[:, :, islice].astype(np.float32), size=3)
        igood = volgood[:, :, islice]
        gcurr = ndimage.sobel(scurr)
        gout = gcurr / np.maximum(scurr, 1e-6)
        ipts = np.flatnonzero((gout > thresuse) & igood)
        if ipts.size == 0:
            continue
        rr, cc = np.unravel_index(ipts, scurr.shape)
        pts = np.column_stack([cc + 1.0, rr + 1.0, np.full(rr.shape, islice + 1.0)])
        ptlist.append(_downsample_points_per_slice(pts, downsample_grid))
    if not ptlist:
        return np.zeros((0, 3), dtype=np.float32)
    ptvol = np.concatenate(ptlist, axis=0).astype(np.float32)
    z0, z1 = ikeep
    keep = (ptvol[:, 2] >= z0) & (ptvol[:, 2] <= z1)
    return ptvol[keep]


def _downsample_points_per_slice(points: np.ndarray, grid: float) -> np.ndarray:
    """Approximate MATLAB pcdownsample(..., 'nonuniformGridSample', grid)."""
    if points.shape[0] <= 6:
        return points
    keys = np.floor(points[:, :2] / grid).astype(int)
    _, idx = np.unique(keys, axis=0, return_index=True)
    return points[np.sort(idx)]
