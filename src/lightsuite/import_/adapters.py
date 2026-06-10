"""Readers for LCT and Arivis annotation exports."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
import tifffile
import zarr

from lightsuite.config.models import AnnotationImportConfig, AnnotationRole
from lightsuite.import_.models import ImportedMask, ImportedPoints
from lightsuite.import_.normalize import (
    reorder_coordinates,
    scale_coordinates_between_voxels,
    to_lightsuite_sample_indices,
)


def _read_lct_zarr_metadata(
    zarr_path: Path,
    level: str,
) -> tuple[list[float], tuple[int, int, int]]:
    attrs_path = zarr_path / ".zattrs"
    if not attrs_path.is_file():
        msg = f"Missing OME-NGFF metadata: {attrs_path}"
        raise FileNotFoundError(msg)
    attrs = json.loads(attrs_path.read_text(encoding="utf-8"))
    multiscales = attrs.get("multiscales", [])
    if not multiscales:
        msg = f"No multiscales metadata in {attrs_path}"
        raise ValueError(msg)
    datasets = multiscales[0].get("datasets", [])
    scale_um: list[float] | None = None
    for entry in datasets:
        if entry.get("path") == level:
            transforms = entry.get("coordinateTransformations", [])
            for tr in transforms:
                if tr.get("type") == "scale":
                    scale_um = [float(v) for v in tr["scale"]]
                    break
            break
    if scale_um is None:
        msg = f"Level {level!r} not found in {zarr_path}"
        raise ValueError(msg)

    arr = zarr.open_array(str(zarr_path / level), mode="r")
    # OME-NGFF declares z,y,x; LightSuite volumes use y,x,z.
    nz, ny, nx = (int(v) for v in arr.shape)
    voxel_um = [float(scale_um[2]), float(scale_um[1]), float(scale_um[0])]
    return voxel_um, (ny, nx, nz)


_MAX_MASK_BYTES = 4_000_000_000


def load_lct_zarr_mask(spec: AnnotationImportConfig) -> ImportedMask:
    zarr_path = spec.path.expanduser()
    if not zarr_path.is_dir():
        msg = f"LCT zarr path not found: {zarr_path}"
        raise FileNotFoundError(msg)
    voxel_um, shape_yxz = _read_lct_zarr_metadata(zarr_path, spec.level)
    if spec.voxel_um is not None:
        voxel_um = [float(v) for v in spec.voxel_um]

    arr = zarr.open_array(str(zarr_path / spec.level), mode="r")
    if arr.nbytes > _MAX_MASK_BYTES:
        msg = (
            f"LCT zarr level {spec.level!r} is {arr.nbytes / 1e9:.1f} GB — too large to load. "
            "Use a coarser pyramid level (e.g. level_05) in import.annotations."
        )
        raise ValueError(msg)
    # z,y,x -> y,x,z
    volume = np.asarray(arr[:], dtype=np.uint8)
    volume = np.transpose(volume, (1, 2, 0))
    label = spec.label or zarr_path.stem
    return ImportedMask(
        label=label,
        volume=volume,
        voxel_um=voxel_um,
        source_path=zarr_path,
        metadata={"level": spec.level, "shape_yxz": shape_yxz},
    )


def _load_mask_volume(path: Path) -> np.ndarray:
    """Load a 3D mask TIFF as (Y, X, Z), supporting page stacks and single 3D pages."""
    path = path.expanduser()
    with tifffile.TiffFile(path) as tif:
        if len(tif.pages) == 0:
            msg = f"No TIFF pages in {path}"
            raise ValueError(msg)
        if len(tif.pages) == 1:
            page = np.asarray(tif.pages[0].asarray(), dtype=np.float32)
            if page.ndim == 3:
                # tifffile.imwrite(3d_array) stores Z,Y,X in a single page.
                return np.transpose(page, (1, 2, 0))
            if page.ndim == 2:
                msg = f"2D TIFF mask not supported (expected Z-stack): {path}"
                raise ValueError(msg)
        planes = [np.asarray(page.asarray(), dtype=np.float32) for page in tif.pages]
        if any(plane.ndim != 2 for plane in planes):
            shapes = [plane.shape for plane in planes[:3]]
            msg = f"Expected 2D TIFF pages in {path}, got shapes {shapes}"
            raise ValueError(msg)
        stack_zyx = np.stack(planes, axis=0)
        return np.transpose(stack_zyx, (1, 2, 0))


def load_tiff_mask(spec: AnnotationImportConfig) -> ImportedMask:
    """Load a single- or multi-page TIFF mask as (Y, X, Z) uint8."""
    tiff_path = spec.path.expanduser()
    if not tiff_path.is_file():
        msg = f"TIFF mask not found: {tiff_path}"
        raise FileNotFoundError(msg)

    volume = _load_mask_volume(tiff_path)
    if volume.ndim != 3:
        msg = f"TIFF mask must be a 3D stack (Y, X, Z), got shape {volume.shape} in {tiff_path}"
        raise ValueError(msg)

    nbytes = int(np.prod(volume.shape))
    if nbytes > _MAX_MASK_BYTES:
        msg = (
            f"TIFF mask is {nbytes / 1e9:.1f} GB — too large to load in memory. "
            "Downsample the mask or use lct_zarr with a coarse pyramid level."
        )
        raise ValueError(msg)

    mask = (volume > 0).astype(np.uint8)
    label = spec.label or tiff_path.stem
    if spec.voxel_um is not None:
        voxel_um = [float(v) for v in spec.voxel_um]
    else:
        # Resolved from sample.voxel_um in brain_import when unset.
        voxel_um = [1.0, 1.0, 1.0]
    return ImportedMask(
        label=label,
        volume=mask,
        voxel_um=voxel_um,
        source_path=tiff_path,
        metadata={
            "format": "tiff_mask",
            "shape_yxz": tuple(int(v) for v in mask.shape),
            "voxel_um_from_sample": spec.voxel_um is None,
        },
    )


def load_lct_json_coords(spec: AnnotationImportConfig) -> ImportedPoints:
    json_path = spec.path.expanduser()
    raw = json.loads(json_path.read_text(encoding="utf-8"))
    if not isinstance(raw, list):
        msg = f"Expected list of coordinates in {json_path}"
        raise ValueError(msg)
    coords = np.asarray(raw, dtype=np.float64)
    if coords.ndim != 2 or coords.shape[1] < 3:
        msg = f"Expected Nx3 coordinate list in {json_path}, got shape {coords.shape}"
        raise ValueError(msg)

    xyz = reorder_coordinates(coords, axis_order=spec.axis_order)
    features = coords[:, 3:] if coords.shape[1] > 3 else None
    label = spec.label or json_path.stem
    return ImportedPoints(
        label=label,
        coordinates=xyz,
        features=features,
        source_path=json_path,
        metadata={"format": "lct_json_coords"},
    )


def load_arivis_csv(spec: AnnotationImportConfig) -> ImportedPoints:
    csv_path = spec.path.expanduser()
    with csv_path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
    if not rows:
        msg = f"No rows in Arivis CSV: {csv_path}"
        raise ValueError(msg)

    x_key = "X (px), Center of Mass (Intensities) #1"
    y_key = "Y (px), Center of Mass (Intensities) #1"
    z_key = "Z (px), Center of Mass (Intensities) #1"
    for key in (x_key, y_key, z_key):
        if key not in rows[0]:
            msg = f"Missing Arivis column {key!r} in {csv_path}"
            raise KeyError(msg)

    coords = np.column_stack(
        [
            [float(row[x_key]) for row in rows],
            [float(row[y_key]) for row in rows],
            [float(row[z_key]) for row in rows],
        ]
    )
    voxel_key = "VoxelCount, Volume"
    volume_key = "Volume, Volume (µm³)"
    features = None
    if voxel_key in rows[0]:
        features = np.column_stack(
            [
                [float(row[voxel_key]) for row in rows],
                [float(row.get(volume_key, 0.0)) for row in rows],
            ]
        )

    label = spec.label or csv_path.stem
    return ImportedPoints(
        label=label,
        coordinates=coords,
        features=features,
        source_path=csv_path,
        metadata={"format": "arivis_csv", "n_segments": len(rows)},
    )


def load_annotation(spec: AnnotationImportConfig) -> ImportedPoints | ImportedMask:
    from lightsuite.config.models import AnnotationFormat

    if spec.format == AnnotationFormat.LCT_ZARR:
        if spec.role != AnnotationRole.MASK:
            msg = "lct_zarr imports must use role: mask"
            raise ValueError(msg)
        return load_lct_zarr_mask(spec)
    if spec.format == AnnotationFormat.TIFF_MASK:
        if spec.role != AnnotationRole.MASK:
            msg = "tiff_mask imports must use role: mask"
            raise ValueError(msg)
        return load_tiff_mask(spec)
    if spec.format == AnnotationFormat.LCT_JSON_COORDS:
        return load_lct_json_coords(spec)
    if spec.format == AnnotationFormat.ARIVIS_CSV:
        return load_arivis_csv(spec)
    msg = f"Unsupported annotation format: {spec.format}"
    raise ValueError(msg)


def prepare_points_for_sample(
    points: ImportedPoints,
    spec: AnnotationImportConfig,
    *,
    target_voxel_um: list[float],
    target_size_yxz: tuple[int, int, int],
) -> ImportedPoints:
    """Normalize imported coordinates into LightSuite 1-based sample indices."""
    xyz = reorder_coordinates(points.coordinates, axis_order=spec.axis_order)
    if spec.voxel_um is not None:
        xyz = scale_coordinates_between_voxels(xyz, spec.voxel_um, target_voxel_um)
    xyz, in_bounds = to_lightsuite_sample_indices(
        xyz,
        index_base=spec.index_base,
        target_size_yxz=target_size_yxz,
    )
    return ImportedPoints(
        label=points.label,
        coordinates=xyz[in_bounds],
        features=points.features[in_bounds] if points.features is not None else None,
        source_path=points.source_path,
        metadata={
            **points.metadata,
            "n_input": int(points.coordinates.shape[0]),
            "n_in_bounds": int(in_bounds.sum()),
            "n_dropped": int(points.coordinates.shape[0] - in_bounds.sum()),
        },
    )
