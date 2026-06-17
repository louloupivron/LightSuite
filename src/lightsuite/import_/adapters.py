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

# Sanity cap on mask voxel count (~64 GiB as uint8). Raise if needed for very large samples.
_MAX_MASK_VOXELS = 64 * 1024**3


def _mask_tiff_shape_yxz(path: Path) -> tuple[int, int, int]:
    """Read (Y, X, Z) shape from a mask TIFF without loading the full volume."""
    with tifffile.TiffFile(path) as tif:
        if len(tif.pages) == 0:
            msg = f"No TIFF pages in {path}"
            raise ValueError(msg)
        if len(tif.pages) == 1:
            page = tif.pages[0]
            arr = page.asarray()
            if arr.ndim == 3:
                return int(arr.shape[1]), int(arr.shape[2]), int(arr.shape[0])
            if arr.ndim == 2:
                msg = f"2D TIFF mask not supported (expected Z-stack): {path}"
                raise ValueError(msg)
            msg = f"Unsupported mask page shape {arr.shape} in {path}"
            raise ValueError(msg)
        ny, nx = (int(v) for v in tif.pages[0].shape[:2])
        return ny, nx, len(tif.pages)


def _load_mask_volume(path: Path) -> np.ndarray:
    """Load a 3D mask TIFF as (Y, X, Z) uint8, supporting page stacks and single 3D pages."""
    path = path.expanduser()
    ny, nx, nz = _mask_tiff_shape_yxz(path)
    if ny * nx * nz > _MAX_MASK_VOXELS:
        msg = (
            f"Mask TIFF has {ny * nx * nz / 1e9:.1f}B voxels — exceeds the "
            f"{_MAX_MASK_VOXELS / 1e9:.0f}B voxel safety limit."
        )
        raise ValueError(msg)

    volume = np.zeros((ny, nx, nz), dtype=np.uint8)
    with tifffile.TiffFile(path) as tif:
        if len(tif.pages) == 1:
            page = np.asarray(tif.pages[0].asarray())
            if page.ndim == 3:
                stack_zyx = page
                for z in range(nz):
                    volume[:, :, z] = (stack_zyx[z] > 0)
                return volume
        for z, page in enumerate(tif.pages):
            plane = np.asarray(page.asarray())
            if plane.ndim != 2:
                shapes = plane.shape
                msg = f"Expected 2D TIFF pages in {path}, got shape {shapes} at page {z}"
                raise ValueError(msg)
            volume[:, :, z] = (plane > 0)
    return volume


def _column_is_numeric(rows: list[dict[str, str]], key: str) -> bool:
    """Return True if every non-empty cell in ``key`` parses as a float."""
    for row in rows:
        raw = row.get(key, "")
        if raw is None or str(raw).strip() == "":
            continue
        try:
            float(raw)
        except (TypeError, ValueError):
            return False
    return True


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
    numeric_keys = [key for key in extra_keys if _column_is_numeric(rows, key)]
    skipped_keys = [key for key in extra_keys if key not in numeric_keys]
    features = None
    if numeric_keys:
        features = np.column_stack(
            [
                [float(row[key]) if str(row.get(key, "")).strip() else np.nan for row in rows]
                for key in numeric_keys
            ]
        )

    label = spec.label or csv_path.stem
    metadata: dict = {"format": "points_csv", "n_points": len(rows)}
    if numeric_keys:
        metadata["feature_columns"] = numeric_keys
    if skipped_keys:
        metadata["skipped_non_numeric_columns"] = skipped_keys

    return ImportedPoints(
        label=label,
        coordinates=coords,
        features=features,
        source_path=csv_path,
        metadata=metadata,
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

    label = spec.label or tiff_path.stem
    return ImportedMask(
        label=label,
        volume=volume,
        voxel_um=[1.0, 1.0, 1.0],
        source_path=tiff_path,
        metadata={"format": "mask_tiff", "shape_yxz": tuple(int(v) for v in volume.shape)},
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
