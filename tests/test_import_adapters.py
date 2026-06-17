"""Tests for LightSuite Sample Space v1 annotation adapters."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
import tifffile

from lightsuite.config.models import AnnotationFormat, AnnotationImportConfig
from lightsuite.import_.adapters import load_mask_tiff, load_points_csv, prepare_points_for_sample
from lightsuite.import_.normalize import filter_in_bounds_points
from lightsuite.import_.sample_reference import (
    SampleReference,
    validate_mask_against_reference,
)


def _reference(*, ny: int = 10, nx: int = 10, nz: int = 10) -> SampleReference:
    return SampleReference.from_checkpoint(
        sample_name="test",
        ny=ny,
        nx=nx,
        nz=nz,
        voxel_um=[2.0, 2.0, 3.0],
    )


def test_load_points_csv(tmp_path: Path) -> None:
    csv_path = tmp_path / "cells.csv"
    csv_path.write_text(
        "x,y,z,intensity\n"
        "1.0,2.0,3.0,100\n"
        "4.5,5.5,6.5,200\n",
        encoding="utf-8",
    )
    spec = AnnotationImportConfig(format=AnnotationFormat.POINTS_CSV, path=csv_path)
    points = load_points_csv(spec)
    assert points.coordinates.shape == (2, 3)
    assert np.allclose(points.coordinates[0], [1.0, 2.0, 3.0])
    assert points.features is not None
    assert points.features.shape == (2, 1)


def test_load_points_csv_skips_non_numeric_extra_columns(tmp_path: Path) -> None:
    csv_path = tmp_path / "cells.csv"
    csv_path.write_text(
        "x,y,z,name,volume_um3\n"
        "1.0,2.0,3.0,Segment #001,27\n"
        "4.5,5.5,6.5,Segment #002,42\n",
        encoding="utf-8",
    )
    spec = AnnotationImportConfig(format=AnnotationFormat.POINTS_CSV, path=csv_path)
    points = load_points_csv(spec)
    assert points.coordinates.shape == (2, 3)
    assert points.features is not None
    assert points.features.shape == (2, 1)
    assert points.features[0, 0] == 27.0
    assert points.metadata["skipped_non_numeric_columns"] == ["name"]
    assert points.metadata["feature_columns"] == ["volume_um3"]


def test_load_points_csv_requires_xyz_header(tmp_path: Path) -> None:
    csv_path = tmp_path / "bad.csv"
    csv_path.write_text("a,b,c\n1,2,3\n", encoding="utf-8")
    spec = AnnotationImportConfig(format=AnnotationFormat.POINTS_CSV, path=csv_path)
    with pytest.raises(KeyError, match="x"):
        load_points_csv(spec)


def test_load_mask_tiff_stack(tmp_path: Path) -> None:
    stack = np.zeros((4, 6, 3), dtype=np.uint16)
    stack[1, 2, 1] = 255
    stack[3, 4, 2] = 128
    tiff_path = tmp_path / "mask.tif"
    tifffile.imwrite(tiff_path, stack, photometric="minisblack")

    spec = AnnotationImportConfig(format=AnnotationFormat.MASK_TIFF, path=tiff_path)
    mask = load_mask_tiff(spec)
    assert mask.volume.shape == (6, 3, 4)
    assert mask.volume[2, 1, 1] == 1
    assert mask.volume[4, 2, 3] == 1


def test_filter_in_bounds_points() -> None:
    pts = np.array([[1.0, 1.0, 1.0], [10.0, 10.0, 10.0], [0.5, 1.0, 1.0]])
    _, mask = filter_in_bounds_points(pts, target_size_yxz=(10, 10, 10))
    assert mask.tolist() == [True, True, False]


def test_prepare_points_drops_out_of_bounds() -> None:
    from lightsuite.import_.models import ImportedPoints

    points = ImportedPoints(
        label="t",
        coordinates=np.array([[1.0, 1.0, 1.0], [100.0, 100.0, 100.0]]),
    )
    prepared = prepare_points_for_sample(points, reference=_reference())
    assert prepared.coordinates.shape[0] == 1
    assert prepared.metadata["n_dropped"] == 1


def test_validate_mask_against_reference_ok() -> None:
    ref = _reference(ny=6, nx=4, nz=3)
    validate_mask_against_reference((6, 4, 3), [2.0, 2.0, 3.0], ref)


def test_validate_mask_against_reference_shape_mismatch() -> None:
    ref = _reference(ny=6, nx=4, nz=3)
    with pytest.raises(ValueError, match="shape"):
        validate_mask_against_reference((5, 4, 3), [2.0, 2.0, 3.0], ref)


def test_sample_reference_roundtrip(tmp_path: Path) -> None:
    ref = _reference(ny=100, nx=200, nz=50)
    path = ref.save(tmp_path / "sample_reference.json")
    loaded = SampleReference.load(path)
    assert loaded.shape_yxz == [100, 200, 50]
    assert loaded.format == "lightsuite_sample_space_v1"
    assert loaded.index_base == 1


def test_registration_shape_native_yxz_matches_preprocess_grid() -> None:
    from lightsuite.import_.brain_import import _registration_shape_native_yxz
    from lightsuite.registration.volume import permute_brain_volume

    shape = _registration_shape_native_yxz((2048, 2048, 1361), [6.55, 6.55, 5.0], 20.0)
    assert shape == (671, 671, 341)
    permuted = permute_brain_volume(np.zeros(shape, dtype=np.uint8), [-1, 3, 2])
    assert permuted.shape == (671, 341, 671)
