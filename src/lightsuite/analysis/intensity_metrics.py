"""Configurable per-region intensity statistics for export."""

from __future__ import annotations

from collections.abc import Sequence

import pandas as pd

DEFAULT_INTENSITY_METRICS: tuple[str, ...] = (
    "median_intensity",
    "std",
    "volume_mm3",
)

OPTIONAL_INTENSITY_METRICS: tuple[str, ...] = (
    "mean_intensity",
    "variance",
)

INTENSITY_METRIC_NAMES: frozenset[str] = frozenset(
    (*DEFAULT_INTENSITY_METRICS, *OPTIONAL_INTENSITY_METRICS, "relative_median_intensity")
)

# Maps tidy metric name → ParcellationResult attribute suffix (…_over_areas).
PARCELLATION_METRIC_FIELDS: dict[str, str] = {
    "median_intensity": "median_over_areas",
    "mean_intensity": "mean_over_areas",
    "std": "std_over_areas",
    "variance": "variance_over_areas",
    "volume_mm3": "volume_over_areas",
}


def normalize_intensity_metrics(metrics: Sequence[str] | None) -> list[str]:
    """Return a de-duplicated metric list, preserving order; default when empty."""
    if not metrics:
        return list(DEFAULT_INTENSITY_METRICS)
    out: list[str] = []
    seen: set[str] = set()
    for item in metrics:
        key = str(item).strip()
        if not key or key in seen:
            continue
        if key not in INTENSITY_METRIC_NAMES:
            allowed = ", ".join(sorted(INTENSITY_METRIC_NAMES))
            msg = f"Unknown intensity metric {key!r}; allowed: {allowed}"
            raise ValueError(msg)
        seen.add(key)
        out.append(key)
    return out or list(DEFAULT_INTENSITY_METRICS)


def filter_intensity_metric_rows(df: pd.DataFrame, metrics: Sequence[str] | None) -> pd.DataFrame:
    """Keep non-intensity rows; drop intensity metrics not in the allowed set."""
    if df.empty or "metric" not in df.columns:
        return df
    allowed = set(normalize_intensity_metrics(metrics))
    is_intensity = df["metric"].isin(INTENSITY_METRIC_NAMES)
    return df.loc[~is_intensity | df["metric"].isin(allowed)].copy()
