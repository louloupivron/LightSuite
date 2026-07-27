"""Tests for SmartSPIM manifest geometry."""

from __future__ import annotations

from pathlib import Path

import pytest

from lightsuite.multires.spec_geometry import (
    physical_bounds_from_spec,
    physical_center_from_spec,
)
from lightsuite.multires.vendor.smartspim import (
    SmartspimGeometryConfig,
    _mosaic_pitch_correction_um,
    _overview_stage_y_anchors_bracketing_roi,
    _physical_y_at_overview_stage_row,
    build_smartspim_pair_manifest,
    parse_smartspim_metadata,
)
from lightsuite.multires.volume import discover_volume_shape

MULTI_RES_ROOT = Path("/media/gbm/NVME2/ALICe-pipelines-data/Multi_RES_SCANs")
OVERVIEW_GEOMETRY = SmartspimGeometryConfig(lateral_flip=(1, -1))
CORTEX_OVERVIEW_GEOMETRY = SmartspimGeometryConfig(
    lateral_flip=(1, -1),
    apply_mosaic_pitch_correction=True,
)
ROI_GEOMETRY = SmartspimGeometryConfig(lateral_flip=(1, -1))
CORTEX_ROI_GEOMETRY = ROI_GEOMETRY
CEREBELLUM_GEOMETRY = SmartspimGeometryConfig(lateral_flip=(1, -1), stage_z_is_center=True)

CORTEX_TILES = [
    (490520, 439090, 4970),
    (503300, 439090, 4970),
    (516080, 439090, 4970),
    (528860, 439090, 4970),
    (490520, 451870, 4970),
    (503300, 451870, 4970),
    (516080, 451870, 4970),
    (528860, 451870, 4970),
]


def _require_multi_res_root() -> Path:
    if not MULTI_RES_ROOT.is_dir():
        pytest.skip("SmartSPIM sample data not available")
    return MULTI_RES_ROOT


def _assert_roi_inside_overview(manifest) -> None:
    ov_min, ov_max = physical_bounds_from_spec(manifest.overview)
    roi_min, roi_max = physical_bounds_from_spec(manifest.roi)
    assert roi_min[0] >= ov_min[0] - 1e-6
    assert roi_min[1] >= ov_min[1] - 1e-6
    assert roi_min[2] >= ov_min[2] - 1e-6
    assert roi_max[0] <= ov_max[0] + 1e-6
    assert roi_max[1] <= ov_max[1] + 1e-6
    assert roi_max[2] <= ov_max[2] + 1e-6


def test_cortex_mosaic_geometry_with_unified_flip() -> None:
    root = _require_multi_res_root()
    manifest = build_smartspim_pair_manifest(
        sample_name="Multi_RES_SCANs",
        pair_label="cortex_9x_561",
        overview_path=root / "1_6X/All_Channels",
        roi_path=root / "9X/20260710_11_05_27_9X_4x2_cortex_cc/All_Channels",
        overview_meta_path=root / "1_6X/metadata.txt",
        roi_meta_path=root / "9X/20260710_11_05_27_9X_4x2_cortex_cc/metadata.txt",
        roi_tile_centers_stage=CORTEX_TILES,
        overview_geometry=CORTEX_OVERVIEW_GEOMETRY,
        roi_geometry=CORTEX_ROI_GEOMETRY,
    )

    roi_center = physical_center_from_spec(manifest.roi)
    assert roi_center[0] == pytest.approx(5096.9, abs=50)
    assert roi_center[1] == pytest.approx(4454.8, abs=50)
    assert roi_center[2] == pytest.approx(5788.5, abs=1.0)
    assert manifest.overview.origin_um[0] == pytest.approx(-9900.4, abs=1.0)
    assert manifest.overview.origin_um[1] == pytest.approx(13858.8, abs=0.1)
    assert manifest.overview.origin_um[2] == pytest.approx(1489.0, abs=0.1)
    assert manifest.roi.origin_um[2] == pytest.approx(4970.0, abs=0.1)


def test_overview_y_anchor_brackets_roi_stage_rows() -> None:
    root = _require_multi_res_root()
    overview_tiles = parse_smartspim_metadata(root / "1_6X/metadata.txt").tile_centers_stage
    cortex_tiles = parse_smartspim_metadata(
        root / "9X/20260710_11_05_27_9X_4x2_cortex_cc/metadata.txt"
    ).tile_centers_stage
    cerebellum_tiles = parse_smartspim_metadata(
        root / "9X/20260710_10_04_38_9X_cerebellum/metadata.txt"
    ).tile_centers_stage

    assert _overview_stage_y_anchors_bracketing_roi(overview_tiles, cortex_tiles) == (
        430080.0,
        487680.0,
    )
    assert _overview_stage_y_anchors_bracketing_roi(overview_tiles, cerebellum_tiles) == (
        487680.0,
        545280.0,
    )


def test_physical_y_at_overview_stage_row_uses_tile_row_pixel_index() -> None:
    root = _require_multi_res_root()
    overview_meta = parse_smartspim_metadata(root / "1_6X/metadata.txt")
    _, overview_ny, overview_nx = discover_volume_shape(root / "1_6X/All_Channels")
    mean_stage_y = 445480.0
    physical_y = _physical_y_at_overview_stage_row(
        mean_stage_y,
        overview_tiles=overview_meta.tile_centers_stage,
        shape_yx=(overview_ny, overview_nx),
        um_per_pix=overview_meta.um_per_pix,
        vres=overview_meta.vres,
        geometry=OVERVIEW_GEOMETRY,
    )
    assert physical_y == pytest.approx(-428.09, abs=1.0)


def test_cerebellum_mosaic_geometry_with_unified_flip() -> None:
    root = _require_multi_res_root()
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
        roi_tile_centers_stage=cerebellum_tiles,
        overview_geometry=CEREBELLUM_GEOMETRY,
        roi_geometry=CEREBELLUM_GEOMETRY,
    )

    roi_center = physical_center_from_spec(manifest.roi)
    assert roi_center[0] == pytest.approx(4995.0, abs=50)
    assert roi_center[1] == pytest.approx(5314.6, abs=50)
    assert manifest.overview.origin_um[2] == pytest.approx(-2411.0, abs=0.1)
    assert manifest.roi.origin_um[2] == pytest.approx(1336.5, abs=0.1)
    _assert_roi_inside_overview(manifest)


def test_cortex_upper_left_tile_corner_placement() -> None:
    root = _require_multi_res_root()
    manifest = build_smartspim_pair_manifest(
        sample_name="Multi_RES_SCANs",
        pair_label="cortex_9x_561_ul_corner",
        overview_path=root / "1_6X/All_Channels",
        roi_path=root / "9X/20260710_11_05_27_9X_4x2_cortex_cc/All_Channels",
        overview_meta_path=root / "1_6X/metadata.txt",
        roi_meta_path=root / "9X/20260710_11_05_27_9X_4x2_cortex_cc/metadata.txt",
        roi_tile_centers_stage=CORTEX_TILES,
        overview_geometry=OVERVIEW_GEOMETRY,
        roi_geometry=SmartspimGeometryConfig(
            lateral_flip=(1, -1),
            xy_placement="upper_left_tile_corner",
        ),
    )

    roi_center = physical_center_from_spec(manifest.roi)
    assert manifest.roi.origin_um[0] == pytest.approx(4195.2, abs=1.0)
    assert manifest.roi.origin_um[1] == pytest.approx(5228.7, abs=1.0)
    assert manifest.roi.origin_um[2] == pytest.approx(4970.0, abs=0.1)
    assert roi_center[0] == pytest.approx(6823.6, abs=50)
    assert roi_center[1] == pytest.approx(3878.3, abs=50)


def test_single_fov_stage_center_interpretation() -> None:
    root = _require_multi_res_root()
    manifest = build_smartspim_pair_manifest(
        sample_name="Multi_RES_SCANs",
        pair_label="single_fov_561",
        overview_path=root / "1_6X/All_Channels",
        roi_path=root / "9X/20260710_11_03_38_9X_single__FOV/Ex_561_Em_561F",
        overview_meta_path=root / "1_6X/metadata.txt",
        roi_meta_path=root / "9X/20260710_11_03_38_9X_single__FOV/metadata.txt",
        overview_geometry=OVERVIEW_GEOMETRY,
        roi_geometry=ROI_GEOMETRY,
    )

    roi_center = physical_center_from_spec(manifest.roi)
    assert roi_center[0] == pytest.approx(5284.9, abs=50)
    assert roi_center[1] == pytest.approx(4702.9, abs=50)
    assert manifest.overview.origin_um[2] == pytest.approx(1489.0, abs=0.1)
    assert manifest.roi.origin_um[2] == pytest.approx(4670.0, abs=0.1)
    _assert_roi_inside_overview(manifest)


def test_mosaic_pair_does_not_auto_enable_anchor() -> None:
    root = _require_multi_res_root()
    manifest = build_smartspim_pair_manifest(
        sample_name="Multi_RES_SCANs",
        pair_label="cortex_9x_561",
        overview_path=root / "1_6X/All_Channels",
        roi_path=root / "9X/20260710_11_05_27_9X_4x2_cortex_cc/All_Channels",
        overview_meta_path=root / "1_6X/metadata.txt",
        roi_meta_path=root / "9X/20260710_11_05_27_9X_4x2_cortex_cc/metadata.txt",
        roi_tile_centers_stage=CORTEX_TILES,
        overview_geometry=CORTEX_OVERVIEW_GEOMETRY,
        roi_geometry=CORTEX_ROI_GEOMETRY,
    )
    assert manifest.provenance["roi_y_anchored_to_overview_stage"] == "False"
    assert manifest.provenance["roi_x_anchored_to_overview_stage"] == "False"


def test_overview_mosaic_pitch_correction_matches_slice_qc_sweep() -> None:
    root = _require_multi_res_root()
    overview_meta = parse_smartspim_metadata(root / "1_6X/metadata.txt")
    geometry = SmartspimGeometryConfig(lateral_flip=(1, -1))
    dx, dy = _mosaic_pitch_correction_um(
        overview_meta.tile_centers_stage,
        hres=overview_meta.hres,
        vres=overview_meta.vres,
        um_per_pix=overview_meta.um_per_pix,
        geometry=geometry,
    )
    assert dx == pytest.approx(7280.0, abs=1.0)
    assert dy == pytest.approx(-11648.0, abs=1.0)  # computed but not applied yet (Y TBD)


def test_roi_y_origin_respects_lateral_flip_sign() -> None:
    root = _require_multi_res_root()
    manifest = build_smartspim_pair_manifest(
        sample_name="Multi_RES_SCANs",
        pair_label="cortex_9x_561_anchor",
        overview_path=root / "1_6X/All_Channels",
        roi_path=root / "9X/20260710_11_05_27_9X_4x2_cortex_cc/All_Channels",
        overview_meta_path=root / "1_6X/metadata.txt",
        roi_meta_path=root / "9X/20260710_11_05_27_9X_4x2_cortex_cc/metadata.txt",
        roi_tile_centers_stage=CORTEX_TILES,
        overview_geometry=OVERVIEW_GEOMETRY,
        roi_geometry=SmartspimGeometryConfig(
            lateral_flip=(1, -1),
            anchor_roi_y_to_overview_stage=True,
        ),
    )
    roi_center = physical_center_from_spec(manifest.roi)
    assert roi_center[1] == pytest.approx(-428.09, abs=1.0)


def test_stage_z_is_center_adjusts_origin() -> None:
    root = _require_multi_res_root()
    geometry = SmartspimGeometryConfig(lateral_flip=(1, -1), stage_z_is_center=True)
    manifest = build_smartspim_pair_manifest(
        sample_name="Multi_RES_SCANs",
        pair_label="cortex_z_center",
        overview_path=root / "1_6X/All_Channels",
        roi_path=root / "9X/20260710_11_05_27_9X_4x2_cortex_cc/All_Channels",
        overview_meta_path=root / "1_6X/metadata.txt",
        roi_meta_path=root / "9X/20260710_11_05_27_9X_4x2_cortex_cc/metadata.txt",
        roi_tile_centers_stage=CORTEX_TILES,
        overview_geometry=geometry,
        roi_geometry=geometry,
    )
    assert manifest.overview.origin_um[2] == pytest.approx(-2411.0, abs=0.1)
    assert manifest.roi.origin_um[2] == pytest.approx(4151.5, abs=0.1)


def test_roi_z_bounds_inside_overview() -> None:
    root = _require_multi_res_root()
    pairs = [
        (
            "cortex_9x_561",
            root / "9X/20260710_11_05_27_9X_4x2_cortex_cc/All_Channels",
            root / "9X/20260710_11_05_27_9X_4x2_cortex_cc/metadata.txt",
            CORTEX_TILES,
        ),
        (
            "cerebellum_9x_561",
            root / "9X/20260710_10_04_38_9X_cerebellum/All_Channels",
            root / "9X/20260710_10_04_38_9X_cerebellum/metadata.txt",
            parse_smartspim_metadata(
                root / "9X/20260710_10_04_38_9X_cerebellum/metadata.txt"
            ).tile_centers_stage,
        ),
        (
            "single_fov_561",
            root / "9X/20260710_11_03_38_9X_single__FOV/Ex_561_Em_561F",
            root / "9X/20260710_11_03_38_9X_single__FOV/metadata.txt",
            None,
        ),
    ]
    for pair_label, roi_path, roi_meta_path, roi_tiles in pairs:
        manifest = build_smartspim_pair_manifest(
            sample_name="Multi_RES_SCANs",
            pair_label=pair_label,
            overview_path=root / "1_6X/All_Channels",
            roi_path=roi_path,
            overview_meta_path=root / "1_6X/metadata.txt",
            roi_meta_path=roi_meta_path,
            roi_tile_centers_stage=roi_tiles,
            overview_geometry=(
                CORTEX_OVERVIEW_GEOMETRY
                if pair_label == "cortex_9x_561"
                else CEREBELLUM_GEOMETRY
                if pair_label == "cerebellum_9x_561"
                else OVERVIEW_GEOMETRY
            ),
            roi_geometry=(
                CORTEX_ROI_GEOMETRY
                if pair_label == "cortex_9x_561"
                else CEREBELLUM_GEOMETRY
                if pair_label == "cerebellum_9x_561"
                else ROI_GEOMETRY
            ),
        )
        if pair_label != "cortex_9x_561":
            _assert_roi_inside_overview(manifest)
