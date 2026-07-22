"""Tests for SmartSPIM manifest geometry."""

from __future__ import annotations

from pathlib import Path

import pytest

from lightsuite.multires.spec_geometry import physical_center_from_spec
from lightsuite.multires.vendor.smartspim import (
    SmartspimGeometryConfig,
    build_smartspim_pair_manifest,
    parse_smartspim_metadata,
)


def test_stage_y_anchor_moves_roi_center_y() -> None:
    root = Path("/media/gbm/NVME2/ALICe-pipelines-data/Multi_RES_SCANs")
    if not root.is_dir():
        pytest.skip("SmartSPIM sample data not available")

    overview_tiles = [
        (463160, 430080, 1489),
        (535160, 430080, 1489),
        (463160, 487680, 1489),
        (535160, 487680, 1489),
    ]
    cortex_tiles = [
        (490520, 439090, 4970),
        (503300, 439090, 4970),
        (516080, 439090, 4970),
        (528860, 439090, 4970),
        (490520, 451870, 4970),
        (503300, 451870, 4970),
        (516080, 451870, 4970),
        (528860, 451870, 4970),
    ]
    overview_geometry = SmartspimGeometryConfig(lateral_flip=(1, -1))
    roi_geometry = SmartspimGeometryConfig(
        lateral_flip=(1, 1),
        anchor_roi_y_to_overview_stage=True,
    )
    manifest = build_smartspim_pair_manifest(
        sample_name="Multi_RES_SCANs",
        pair_label="cortex_9x_561",
        overview_path=root / "1_6X/All_Channels",
        roi_path=root / "9X/20260710_11_05_27_9X_4x2_cortex_cc/All_Channels",
        overview_meta_path=root / "1_6X/metadata.txt",
        roi_meta_path=root / "9X/20260710_11_05_27_9X_4x2_cortex_cc/metadata.txt",
        overview_tile_centers_stage=overview_tiles,
        roi_tile_centers_stage=cortex_tiles,
        overview_geometry=overview_geometry,
        roi_geometry=roi_geometry,
    )

    roi_center = physical_center_from_spec(manifest.roi)
    assert roi_center[0] == pytest.approx(5096.9, abs=50)
    assert roi_center[1] == pytest.approx(8768, abs=300)


def test_cerebellum_mosaic_geometry_with_full_overview_tiles() -> None:
    root = Path("/media/gbm/NVME2/ALICe-pipelines-data/Multi_RES_SCANs")
    if not root.is_dir():
        pytest.skip("SmartSPIM sample data not available")

    overview_tiles = parse_smartspim_metadata(root / "1_6X/metadata.txt").tile_centers_stage
    cerebellum_tiles = parse_smartspim_metadata(
        root / "9X/20260710_10_04_38_9X_cerebellum/metadata.txt"
    ).tile_centers_stage
    manifest = build_smartspim_pair_manifest(
        sample_name="Multi_RES_SCANs",
        pair_label="cerebellum_9x_561",
        overview_path=root / "1_6X/All_Channels",
        roi_path=root / "9X/20260710_10_04_38_9X_cerebellum/All_Channels",
        overview_meta_path=root / "1_6X/metadata.txt",
        roi_meta_path=root / "9X/20260710_10_04_38_9X_cerebellum/metadata.txt",
        overview_tile_centers_stage=overview_tiles,
        roi_tile_centers_stage=cerebellum_tiles,
        overview_geometry=SmartspimGeometryConfig(lateral_flip=(1, -1)),
        roi_geometry=SmartspimGeometryConfig(anchor_roi_y_to_overview_stage=True),
    )

    roi_center = physical_center_from_spec(manifest.roi)
    assert roi_center[0] == pytest.approx(4994.7, abs=50)
    assert roi_center[1] == pytest.approx(-1950.1, abs=300)


def _single_fov_overview_tiles() -> list[tuple[float, float, float]]:
    return [
        (463160, 430080, 1489),
        (535160, 430080, 1489),
        (463160, 487680, 1489),
        (535160, 487680, 1489),
    ]


def test_single_fov_corner_stage_interpretation() -> None:
    root = Path("/media/gbm/NVME2/ALICe-pipelines-data/Multi_RES_SCANs")
    if not root.is_dir():
        pytest.skip("SmartSPIM sample data not available")

    manifest = build_smartspim_pair_manifest(
        sample_name="Multi_RES_SCANs",
        pair_label="single_fov_561",
        overview_path=root / "1_6X/All_Channels",
        roi_path=root / "9X/20260710_11_03_38_9X_single__FOV/Ex_561_Em_561F",
        overview_meta_path=root / "1_6X/metadata.txt",
        roi_meta_path=root / "9X/20260710_11_03_38_9X_single__FOV/metadata.txt",
        overview_tile_centers_stage=_single_fov_overview_tiles(),
        overview_geometry=SmartspimGeometryConfig(lateral_flip=(1, -1)),
        roi_geometry=SmartspimGeometryConfig(
            stage_xy_is_center=False,
            stage_origin_offset_um=(2000.0, 587.5, 0.0),
        ),
    )

    roi_center = physical_center_from_spec(manifest.roi)
    assert roi_center[0] == pytest.approx(7994.5, abs=50)
    assert roi_center[1] == pytest.approx(6000.0, abs=50)
    assert manifest.overview.origin_um[2] == pytest.approx(1489.0, abs=0.1)
    assert manifest.roi.origin_um[2] == pytest.approx(4670.0, abs=0.1)
