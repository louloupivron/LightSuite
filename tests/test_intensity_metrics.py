"""Tests for configurable intensity metrics."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from lightsuite.analysis.intensity_metrics import (
    filter_intensity_metric_rows,
    normalize_intensity_metrics,
)
from lightsuite.analysis.region_stats import parcellation_result_to_tidy
from lightsuite.export.parcellation import ParcellationResult


def test_normalize_intensity_metrics_defaults_when_empty() -> None:
    assert normalize_intensity_metrics([]) == [
        "median_intensity",
        "std",
        "volume_mm3",
    ]


def test_normalize_intensity_metrics_rejects_unknown() -> None:
    with pytest.raises(ValueError, match="Unknown intensity metric"):
        normalize_intensity_metrics(["not_a_metric"])


def test_parcellation_result_to_tidy_respects_metric_filter() -> None:
    result = ParcellationResult(
        area_ids=np.array([1], dtype=np.int64),
        median_over_areas=np.array([[10.0, 20.0]], dtype=np.float32),
        mean_over_areas=np.array([[11.0, 21.0]], dtype=np.float32),
        std_over_areas=np.array([[1.0, 2.0]], dtype=np.float32),
        variance_over_areas=np.array([[1.0, 4.0]], dtype=np.float32),
        volume_over_areas=np.array([[0.1, 0.2]], dtype=np.float32),
    )
    tidy = parcellation_result_to_tidy(
        result,
        None,
        sample="m1",
        channel=1,
        atlas="allen",
        intensity_metrics=["mean_intensity", "variance"],
    )
    assert set(tidy["metric"]) == {"mean_intensity", "variance"}


def test_filter_intensity_metric_rows_keeps_point_metrics() -> None:
    df = pd.DataFrame(
        {
            "metric": ["median_intensity", "cell_count"],
            "value": [1.0, 2.0],
        }
    )
    filtered = filter_intensity_metric_rows(df, ["volume_mm3"])
    assert filtered["metric"].tolist() == ["cell_count"]
