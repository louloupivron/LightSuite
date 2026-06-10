"""Tests for external annotation format adapters."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest
import tifffile

from lightsuite.config.models import AnnotationFormat, AnnotationImportConfig, AnnotationRole
from lightsuite.import_.adapters import (
    _read_lct_zarr_metadata,
    load_arivis_csv,
    load_lct_json_coords,
    load_lct_zarr_mask,
    load_tiff_mask,
    prepare_points_for_sample,
)
from lightsuite.import_.normalize import reorder_coordinates, to_lightsuite_sample_indices

IMPORT_ROOT = Path("/media/gbm/NVME2/ALICe-pipelines-data/registration_imports")


@pytest.mark.skipif(
    not (IMPORT_ROOT / "fromArivis/1_2ch_1X-features (2).csv").is_file(),
    reason="Arivis sample CSV not available",
)
def test_load_arivis_csv_sample() -> None:
    spec = AnnotationImportConfig(
        format=AnnotationFormat.ARIVIS_CSV,
        path=IMPORT_ROOT / "fromArivis/1_2ch_1X-features (2).csv",
        axis_order="xyz",
        index_base=0,
    )
    points = load_arivis_csv(spec)
    assert points.coordinates.shape[0] == 64
    assert points.coordinates.shape[1] == 3
    assert points.coordinates[0, 0] > 800


@pytest.mark.skipif(
    not (IMPORT_ROOT / "fromLCT/120_MCIrest_spinalcord_sensitive.json").is_file(),
    reason="LCT JSON sample not available",
)
def test_load_lct_json_coords_sample() -> None:
    spec = AnnotationImportConfig(
        format=AnnotationFormat.LCT_JSON_COORDS,
        path=IMPORT_ROOT / "fromLCT/120_MCIrest_spinalcord_sensitive.json",
        axis_order="zyx",
        index_base=0,
    )
    points = load_lct_json_coords(spec)
    assert points.coordinates.shape[0] == 256_423
    assert np.allclose(points.coordinates[0], [1863.0, 102.0, 2.0])


@pytest.mark.skipif(
    not (IMPORT_ROOT / "fromLCT/561_CFosSensitive_cells.zarr/level_01").is_dir(),
    reason="LCT zarr sample not available",
)
def test_load_lct_zarr_mask_metadata() -> None:
    zarr_path = IMPORT_ROOT / "fromLCT/561_CFosSensitive_cells.zarr"
    voxel_um, shape = _read_lct_zarr_metadata(zarr_path, "level_01")
    assert shape == (11692, 7441, 2057)
    assert voxel_um == pytest.approx([1.8, 1.8, 4.0])


def test_load_lct_zarr_mask_coarse_level(tmp_path: Path) -> None:
    import zarr

    root = tmp_path / "mask.zarr"
    store = zarr.storage.LocalStore(str(root))
    zarr.array(
        np.ones((4, 8, 6), dtype=np.uint8),
        store=store,
        path="level_01",
    )
    (root / ".zattrs").write_text(
        json.dumps(
            {
                "multiscales": [
                    {
                        "datasets": [
                            {
                                "path": "level_01",
                                "coordinateTransformations": [
                                    {"type": "scale", "scale": [4.0, 1.8, 1.8]}
                                ],
                            }
                        ]
                    }
                ]
            }
        ),
        encoding="utf-8",
    )
    spec = AnnotationImportConfig(
        format=AnnotationFormat.LCT_ZARR,
        path=root,
        role=AnnotationRole.MASK,
        level="level_01",
    )
    mask = load_lct_zarr_mask(spec)
    assert mask.volume.shape == (8, 6, 4)


def test_load_tiff_mask_stack(tmp_path: Path) -> None:
    stack = np.zeros((4, 6, 3), dtype=np.uint16)
    stack[1, 2, 1] = 255
    stack[3, 4, 2] = 128
    tiff_path = tmp_path / "mask.tif"
    tifffile.imwrite(tiff_path, stack, photometric="minisblack")

    spec = AnnotationImportConfig(
        format=AnnotationFormat.TIFF_MASK,
        path=tiff_path,
        role=AnnotationRole.MASK,
        voxel_um=[2.0, 2.0, 3.0],
    )
    mask = load_tiff_mask(spec)
    assert mask.volume.shape == (6, 3, 4)
    assert mask.volume[2, 1, 1] == 1
    assert mask.volume[4, 2, 3] == 1
    assert mask.voxel_um == pytest.approx([2.0, 2.0, 3.0])
    assert mask.metadata["voxel_um_from_sample"] is False


def test_load_tiff_mask_defaults_voxel_um_flag(tmp_path: Path) -> None:
    stack = np.ones((2, 2, 2), dtype=np.uint8)
    tiff_path = tmp_path / "mask.tif"
    tifffile.imwrite(tiff_path, stack)
    spec = AnnotationImportConfig(
        format=AnnotationFormat.TIFF_MASK,
        path=tiff_path,
        role=AnnotationRole.MASK,
    )
    mask = load_tiff_mask(spec)
    assert mask.metadata["voxel_um_from_sample"] is True


def test_reorder_coordinates_zyx_to_xyz() -> None:
    pts = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
    out = reorder_coordinates(pts, axis_order="zyx")
    assert np.allclose(out, [[3.0, 2.0, 1.0], [6.0, 5.0, 4.0]])


def test_to_lightsuite_sample_indices_zero_based() -> None:
    pts = np.array([[0.0, 0.0, 0.0], [9.0, 9.0, 9.0]])
    converted, mask = to_lightsuite_sample_indices(
        pts,
        index_base=0,
        target_size_yxz=(10, 10, 10),
    )
    assert np.allclose(converted[0], [1, 1, 1])
    assert mask.all()


def test_prepare_points_drops_out_of_bounds() -> None:
    from lightsuite.import_.models import ImportedPoints

    spec = AnnotationImportConfig(
        format=AnnotationFormat.ARIVIS_CSV,
        path=Path("dummy.csv"),
        axis_order="xyz",
        index_base=0,
    )
    points = ImportedPoints(
        label="t",
        coordinates=np.array([[0.0, 0.0, 0.0], [100.0, 100.0, 100.0]]),
    )
    prepared = prepare_points_for_sample(
        points,
        spec,
        target_voxel_um=[5.0, 5.0, 5.0],
        target_size_yxz=(10, 10, 10),
    )
    assert prepared.coordinates.shape[0] == 1
    assert prepared.metadata["n_dropped"] == 1
