"""Readers for LightSuite Sample Space v1 annotation exports."""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
import tifffile

from lightsuite.config.models import AnnotationFormat, AnnotationImportConfig
from lightsuite.import_.models import ImportedMask, ImportedPoints
from lightsuite.import_.normalize import filter_in_bounds_points
from lightsuite.import_.sample_reference import SampleReference

_MAX_MASK_BYTES = 4_000_000_000


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


def load_points_csv(spec: AnnotationImportConfig) -> ImportedPoints:
    """Load 1-based [x, y, z] voxel indices from a CSV with header columns x,y,z."""
    csv_path = spec.path.expanduser()
    if not csv_path.is_file():
        msg = f"Points CSV not found: {csv_path}"
        raise FileNotFoundError(msg)

    with csv_path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            msg = f"Points CSV must include a header row with x,y,z columns: {csv_path}"
            raise ValueError(msg)
        field_map = {name.strip().lower(): name for name in reader.fieldnames}
        for required in ("x", "y", "z"):
            if required not in field_map:
                msg = (
                    f"Points CSV missing required column {required!r} in {csv_path}. "
                    "Expected header: x,y,z"
                )
                raise KeyError(msg)
        rows = list(reader)

    if not rows:
        msg = f"No data rows in points CSV: {csv_path}"
        raise ValueError(msg)

    x_key, y_key, z_key = field_map["x"], field_map["y"], field_map["z"]
    coords = np.column_stack(
        [
            [float(row[x_key]) for row in rows],
            [float(row[y_key]) for row in rows],
            [float(row[z_key]) for row in rows],
        ]
    )
    extra_keys = [field_map[k] for k in field_map if k not in {"x", "y", "z"}]
    features = None
    if extra_keys:
        features = np.column_stack([[float(row[key]) for row in rows] for key in extra_keys])

    label = spec.label or csv_path.stem
    return ImportedPoints(
        label=label,
        coordinates=coords,
        features=features,
        source_path=csv_path,
        metadata={"format": "points_csv", "n_points": len(rows)},
    )


def load_mask_tiff(spec: AnnotationImportConfig) -> ImportedMask:
    """Load a native-resolution binary mask TIFF as (Y, X, Z) uint8."""
    tiff_path = spec.path.expanduser()
    if not tiff_path.is_file():
        msg = f"Mask TIFF not found: {tiff_path}"
        raise FileNotFoundError(msg)

    volume = _load_mask_volume(tiff_path)
    if volume.ndim != 3:
        msg = f"Mask TIFF must be a 3D stack (Y, X, Z), got shape {volume.shape} in {tiff_path}"
        raise ValueError(msg)

    nbytes = int(np.prod(volume.shape))
    if nbytes > _MAX_MASK_BYTES:
        msg = (
            f"Mask TIFF is {nbytes / 1e9:.1f} GB — too large to load in memory. "
            "Split the mask or process on a machine with more RAM."
        )
        raise ValueError(msg)

    mask = (volume > 0).astype(np.uint8)
    label = spec.label or tiff_path.stem
    return ImportedMask(
        label=label,
        volume=mask,
        voxel_um=[1.0, 1.0, 1.0],
        source_path=tiff_path,
        metadata={"format": "mask_tiff", "shape_yxz": tuple(int(v) for v in mask.shape)},
    )


def load_annotation(spec: AnnotationImportConfig) -> ImportedPoints | ImportedMask:
    if spec.format == AnnotationFormat.MASK_TIFF:
        return load_mask_tiff(spec)
    if spec.format == AnnotationFormat.POINTS_CSV:
        return load_points_csv(spec)
    msg = f"Unsupported annotation format: {spec.format}"
    raise ValueError(msg)


def prepare_points_for_sample(
    points: ImportedPoints,
    *,
    reference: SampleReference,
) -> ImportedPoints:
    """Validate 1-based native sample coordinates and drop out-of-bounds points."""
    ny, nx, nz = reference.shape_tuple
    xyz, in_bounds = filter_in_bounds_points(
        points.coordinates,
        target_size_yxz=(ny, nx, nz),
    )
    features = points.features[in_bounds] if points.features is not None else None
    return ImportedPoints(
        label=points.label,
        coordinates=xyz[in_bounds],
        features=features,
        source_path=points.source_path,
        metadata={
            **points.metadata,
            "n_input": int(points.coordinates.shape[0]),
            "n_in_bounds": int(in_bounds.sum()),
            "n_dropped": int(points.coordinates.shape[0] - in_bounds.sum()),
        },
    )
