"""Tests for the tidy region-stats schema and hemisphere helper."""

from __future__ import annotations

import numpy as np
import pandas as pd

from lightsuite.analysis.hemisphere import hemisphere_side_masks, hemisphere_side_volume
from lightsuite.analysis.ontology import RegionTable
from lightsuite.analysis.region_stats import (
    TIDY_COLUMNS,
    parcellation_result_to_tidy,
    tidy_to_wide,
)
from lightsuite.export.parcellation import ParcellationResult


def _region_table() -> RegionTable:
    df = pd.DataFrame(
        {
            "parcellation_index": [1, 2],
            "acronym": ["MOp1", "VPL"],
            "name": ["Primary motor area, Layer 1", "Ventral posterolateral nucleus"],
            "structure": ["MO", "VENT"],
            "division": ["Isocortex", "Thalamus"],
            "division_acronym": ["Isocortex", "TH"],
            "ccf_id": [320, 718],
        }
    )
    return RegionTable(atlas="allen", df=df)


def _result() -> ParcellationResult:
    return ParcellationResult(
        area_ids=np.array([1, 2], dtype=np.int64),
        median_over_areas=np.array([[10.0, 20.0], [30.0, 40.0]], dtype=np.float32),
        mean_over_areas=np.array([[11.0, 21.0], [31.0, 41.0]], dtype=np.float32),
        std_over_areas=np.array([[1.0, 2.0], [3.0, 4.0]], dtype=np.float32),
        variance_over_areas=np.array([[1.0, 4.0], [9.0, 16.0]], dtype=np.float32),
        volume_over_areas=np.array([[0.1, 0.2], [0.3, 0.4]], dtype=np.float32),
    )


def test_parcellation_result_to_tidy_shape_and_metadata() -> None:
    tidy = parcellation_result_to_tidy(
        _result(), _region_table(), sample="m1", channel=2, atlas="allen"
    )
    assert list(tidy.columns) == TIDY_COLUMNS
    # 2 regions × 2 hemispheres × 3 default intensity metrics
    assert len(tidy) == 12

    mop_right_median = tidy[
        (tidy["parcellation_index"] == 1)
        & (tidy["hemisphere"] == "right")
        & (tidy["metric"] == "median_intensity")
    ]
    assert mop_right_median["value"].iloc[0] == 10.0
    assert mop_right_median["division"].iloc[0] == "Isocortex"
    assert (tidy["sample"] == "m1").all()
    assert (tidy["channel"] == 2).all()


def test_tidy_drops_nan_values() -> None:
    result = _result()
    result.median_over_areas[0, 0] = np.nan
    tidy = parcellation_result_to_tidy(result, None, sample="m1", channel=1, atlas="allen")
    rows = tidy[(tidy["parcellation_index"] == 1) & (tidy["metric"] == "median_intensity")]
    # only the left side survives for region 1's median
    assert set(rows["hemisphere"]) == {"left"}


def test_tidy_to_wide_roundtrip_columns() -> None:
    tidy = parcellation_result_to_tidy(
        _result(), _region_table(), sample="m1", channel=1, atlas="allen"
    )
    wide = tidy_to_wide(tidy)
    assert "RightSideIntensity" in wide.columns
    assert "LeftSideIntensity" in wide.columns
    assert "RightSideVolume[mm3]" in wide.columns
    assert "name" in wide.columns
    region2 = wide.set_index("parcellation_index").loc[2]
    assert region2["RightSideIntensity"] == 30.0
    assert region2["LeftSideIntensity"] == 40.0
    assert region2["division"] == "Thalamus"


def test_hemisphere_allen_split_is_z_midpoint() -> None:
    av = np.ones((2, 2, 4), dtype=np.int32)
    side0, side1 = hemisphere_side_masks(av, "allen")
    assert side0[:, :, :2].all() and not side0[:, :, 2:].any()
    assert side1[:, :, 2:].all() and not side1[:, :, :2].any()


def test_hemisphere_side_volume_labels() -> None:
    av = np.ones((2, 2, 4), dtype=np.int32)
    side = hemisphere_side_volume(av, "allen")
    assert set(np.unique(side)) == {0, 1}
    assert (side[:, :, :2] == 0).all()
    assert (side[:, :, 2:] == 1).all()


def test_hemisphere_perens_mean_plane_split() -> None:
    av = np.zeros((2, 2, 4), dtype=np.int32)
    av[..., 1:3] = 5  # brain occupies z=1,2 → mean plane = round(1.5) = 2
    side0, side1 = hemisphere_side_masks(av, "perens", ml_axis=3)
    # side0 = coord <= 2 within brain; side1 = coord > 2
    assert side0[..., 1].all() and side0[..., 2].all()
    assert not side1[..., 1].any() and not side1[..., 2].any()
    assert not side0[..., 0].any()  # background excluded
