"""Tests for SmartSPIM manifest geometry.

Expected ROI positions are the image-content positions measured by normalised
cross-correlation of a mid-stack 9× plane against the 1.6× overview (see
``context.txt``). Metadata placement is expected to reproduce them to within the
cross-resolution stage drift, which is under ~100 µm on this dataset.
"""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import pytest

from lightsuite.multires.spec_geometry import (
    physical_bounds_from_spec,
    physical_center_from_spec,
)
from lightsuite.multires.vendor.smartspim import (
    NOMINAL_TILE_OVERLAP,
    SmartspimGeometryConfig,
    build_smartspim_pair_manifest,
    check_stage_scale,
    parse_smartspim_metadata,
    stage_pitch_overlap_fractions,
)

MULTI_RES_ROOT = Path("/media/gbm/NVME2/ALICe-pipelines-data/Multi_RES_SCANs")

OVERVIEW_PATH = MULTI_RES_ROOT / "1_6X/All_Channels"
OVERVIEW_META = MULTI_RES_ROOT / "1_6X/metadata.txt"

# pair_label -> (roi stack, roi metadata, empirical ROI centre XY in µm, ROI origin XYZ in µm)
PAIRS = {
    "cortex_9x_561": (
        MULTI_RES_ROOT / "9X/20260710_11_05_27_9X_4x2_cortex_cc/All_Channels",
        MULTI_RES_ROOT / "9X/20260710_11_05_27_9X_4x2_cortex_cc/metadata.txt",
        (51012.0, 44614.0),
        (48342.0, 43199.0, 4970.0),
    ),
    "cerebellum_9x_561": (
        MULTI_RES_ROOT / "9X/20260710_10_04_38_9X_cerebellum/All_Channels",
        MULTI_RES_ROOT / "9X/20260710_10_04_38_9X_cerebellum/metadata.txt",
        (49982.0, 53236.0),
        (46042.0, 50519.0, 2771.0),
    ),
    "single_fov_561": (
        MULTI_RES_ROOT / "9X/20260710_11_03_38_9X_single__FOV/Ex_561_Em_561F",
        MULTI_RES_ROOT / "9X/20260710_11_03_38_9X_single__FOV/metadata.txt",
        (52898.0, 47102.0),
        (52139.0, 46319.0, 4670.0),
    ),
}

# Overview stitched pixel [0, 0] = upper-left tile corner in stage µm.
OVERVIEW_ORIGIN_UM = (42316.0, 39808.0, 1489.0)

DRIFT_TOLERANCE_UM = 150.0

# Scans with more than one tile per axis, i.e. those with a stage pitch to check.
MOSAIC_META = [OVERVIEW_META, PAIRS["cortex_9x_561"][1], PAIRS["cerebellum_9x_561"][1]]


def _require_multi_res_root() -> Path:
    if not MULTI_RES_ROOT.is_dir():
        pytest.skip("SmartSPIM sample data not available")
    return MULTI_RES_ROOT


def _build(pair_label: str, **kwargs) -> object:
    roi_path, roi_meta, _centre, _origin = PAIRS[pair_label]
    return build_smartspim_pair_manifest(
        sample_name="Multi_RES_SCANs",
        pair_label=pair_label,
        overview_path=OVERVIEW_PATH,
        roi_path=roi_path,
        overview_meta_path=OVERVIEW_META,
        roi_meta_path=roi_meta,
        **kwargs,
    )


@pytest.mark.parametrize("meta_path", MOSAIC_META, ids=lambda p: p.parent.name)
def test_default_stage_scale_reproduces_hard_coded_tile_overlap(meta_path: Path) -> None:
    _require_multi_res_root()
    meta = parse_smartspim_metadata(meta_path)
    fractions = stage_pitch_overlap_fractions(
        meta.tile_centers_stage,
        meta=meta,
        geometry=SmartspimGeometryConfig(),
    )
    assert fractions, "expected at least one mosaic axis"
    for value in fractions.values():
        assert value == pytest.approx(NOMINAL_TILE_OVERLAP, abs=1e-4)


def test_wrong_stage_scale_is_rejected() -> None:
    _require_multi_res_root()
    meta = parse_smartspim_metadata(OVERVIEW_META)
    with pytest.raises(ValueError, match="hard-codes 10%"):
        check_stage_scale(
            meta.tile_centers_stage,
            meta=meta,
            geometry=SmartspimGeometryConfig(stage_coord_scale_um=0.01),
        )


def test_overview_origin_is_upper_left_tile_corner() -> None:
    _require_multi_res_root()
    manifest = _build("cortex_9x_561")
    assert manifest.overview.origin_um == pytest.approx(OVERVIEW_ORIGIN_UM, abs=0.1)
    assert manifest.overview.spacing_um == pytest.approx([4.0, 4.0, 10.0])
    assert manifest.overview.direction == pytest.approx([1, 0, 0, 0, 1, 0, 0, 0, 1])


def test_overview_placement_is_independent_of_the_roi() -> None:
    _require_multi_res_root()
    origins = {label: _build(label).overview.origin_um for label in PAIRS}
    for origin in origins.values():
        assert origin == pytest.approx(OVERVIEW_ORIGIN_UM, abs=0.1)


@pytest.mark.parametrize("pair_label", list(PAIRS))
def test_roi_origin_matches_upper_left_tile_corner(pair_label: str) -> None:
    _require_multi_res_root()
    _path, _meta, _centre, expected_origin = PAIRS[pair_label]
    manifest = _build(pair_label)
    assert manifest.roi.origin_um == pytest.approx(expected_origin, abs=0.1)


@pytest.mark.parametrize("pair_label", list(PAIRS))
def test_roi_lands_where_the_image_content_says_it_should(pair_label: str) -> None:
    _require_multi_res_root()
    _path, _meta, expected_centre, _origin = PAIRS[pair_label]
    manifest = _build(pair_label)
    centre = physical_center_from_spec(manifest.roi)
    assert centre[0] == pytest.approx(expected_centre[0], abs=DRIFT_TOLERANCE_UM)
    assert centre[1] == pytest.approx(expected_centre[1], abs=DRIFT_TOLERANCE_UM)


@pytest.mark.parametrize("pair_label", list(PAIRS))
def test_roi_is_fully_inside_the_overview(pair_label: str) -> None:
    _require_multi_res_root()
    manifest = _build(pair_label)
    ov_min, ov_max = physical_bounds_from_spec(manifest.overview)
    roi_min, roi_max = physical_bounds_from_spec(manifest.roi)
    assert (roi_min >= ov_min - 1e-6).all()
    assert (roi_max <= ov_max + 1e-6).all()


@pytest.mark.parametrize("pair_label", list(PAIRS))
def test_z_origin_is_the_stage_table_z_column(pair_label: str) -> None:
    _require_multi_res_root()
    _path, roi_meta, _centre, expected_origin = PAIRS[pair_label]
    manifest = _build(pair_label)
    meta = parse_smartspim_metadata(roi_meta)
    assert manifest.roi.origin_um[2] == pytest.approx(meta.tile_centers_stage[0][2], abs=0.1)
    assert manifest.roi.origin_um[2] == pytest.approx(expected_origin[2], abs=0.1)


def test_stage_z_is_center_shifts_origin_by_half_the_stack() -> None:
    _require_multi_res_root()
    manifest = _build("cortex_9x_561", geometry=SmartspimGeometryConfig(stage_z_is_center=True))
    # 1638 planes at 1 µm, table Z 4970 µm
    assert manifest.roi.origin_um[2] == pytest.approx(4970.0 - 1637 / 2.0, abs=0.1)


def test_lateral_flip_mirrors_origin_and_direction() -> None:
    _require_multi_res_root()
    manifest = _build("cortex_9x_561", geometry=SmartspimGeometryConfig(lateral_flip=(1, -1)))
    meta = parse_smartspim_metadata(OVERVIEW_META)
    max_stage_y = max(coord[1] for coord in meta.tile_centers_stage)
    expected = max_stage_y * 0.1 + meta.tile_height_um / 2.0
    assert manifest.overview.origin_um[1] == pytest.approx(expected, abs=0.1)
    assert manifest.overview.direction[4] == pytest.approx(-1.0)


def test_geometry_config_defaults_document_the_measured_convention() -> None:
    geometry = SmartspimGeometryConfig()
    assert geometry.stage_coord_scale_um == 0.1
    assert geometry.stage_xy_is_center is True
    assert geometry.stage_z_is_center is False
    assert geometry.lateral_flip == (1, 1)
    assert replace(geometry, lateral_flip=(1, -1)).lateral_flip == (1, -1)
