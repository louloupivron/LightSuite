"""Cross-subject group summaries and pairwise comparisons on tidy region stats."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

from lightsuite.analysis.cohort_models import CohortConfig, CohortSampleEntry, GroupAnalysisConfig
from lightsuite.analysis.region_stats import METRICS, TIDY_COLUMNS

#: Region-level grouping keys (one row per group × region × measurement context).
_REGION_GROUP_KEYS = [
    "group",
    "channel",
    "hemisphere",
    "metric",
    "parcellation_index",
    "acronym",
    "name",
    "structure",
    "division",
]

#: Keys for pairwise tests (region metadata carried through).
_COMPARE_REGION_KEYS = [
    "channel",
    "hemisphere",
    "metric",
    "parcellation_index",
    "acronym",
    "name",
    "structure",
    "division",
]

_ROLLUP_COLUMN = {
    "division": "division",
    "structure": "structure",
}


@dataclass(frozen=True)
class ResolvedCohortSample:
    """A cohort member with a resolved region-stats path."""

    sample_id: str
    group: str
    region_stats_path: Path


@dataclass
class GroupAnalysisResult:
    """Outputs from a cross-subject group analysis run."""

    cohort_long: pd.DataFrame
    summary_by_region: pd.DataFrame
    comparisons_by_region: pd.DataFrame | None = None
    summary_by_division: pd.DataFrame | None = None
    summary_by_structure: pd.DataFrame | None = None
    comparisons_by_division: pd.DataFrame | None = None
    comparisons_by_structure: pd.DataFrame | None = None


def resolve_cohort_sample(entry: CohortSampleEntry) -> ResolvedCohortSample:
    """Resolve ``region_stats`` path and sample id for one cohort entry."""
    from lightsuite.config.loader import load_config

    if entry.region_stats is not None:
        stats_path = entry.region_stats.expanduser().resolve()
        sample_id = entry.id
        if sample_id is None:
            sample_id = stats_path.parent.parent.name
    else:
        assert entry.config is not None
        brain_cfg = load_config(entry.config)
        stats_path = (
            brain_cfg.sample.save_path.expanduser().resolve()
            / "volume_registered"
            / "region_stats.csv"
        )
        sample_id = entry.id or brain_cfg.sample.name

    if not stats_path.is_file():
        msg = f"region_stats not found for sample {sample_id!r}: {stats_path}"
        raise FileNotFoundError(msg)

    return ResolvedCohortSample(
        sample_id=sample_id,
        group=entry.group,
        region_stats_path=stats_path,
    )


def load_cohort_long(samples: list[ResolvedCohortSample]) -> pd.DataFrame:
    """Load and concatenate per-subject tidy tables with a ``group`` column."""
    frames: list[pd.DataFrame] = []
    for item in samples:
        df = pd.read_csv(item.region_stats_path)
        missing = set(TIDY_COLUMNS) - set(df.columns)
        if missing:
            msg = f"{item.region_stats_path} missing columns: {sorted(missing)}"
            raise ValueError(msg)
        df = df[TIDY_COLUMNS].copy()
        df["sample"] = item.sample_id
        df["group"] = item.group
        frames.append(df)
    if not frames:
        return pd.DataFrame(columns=[*TIDY_COLUMNS, "group"])
    return pd.concat(frames, ignore_index=True)


def filter_cohort_table(df: pd.DataFrame, options: GroupAnalysisConfig) -> pd.DataFrame:
    """Apply channel / metric / hemisphere filters from the cohort config."""
    out = df.copy()
    if options.channels is not None:
        allowed = {_normalize_channel(c) for c in options.channels}
        out = out[out["channel"].map(_normalize_channel).isin(allowed)]
    if options.metrics is not None:
        unknown = set(options.metrics) - set(METRICS)
        if unknown:
            msg = f"Unknown metrics: {sorted(unknown)}. Expected subset of {METRICS}."
            raise ValueError(msg)
        out = out[out["metric"].isin(options.metrics)]
    if options.hemispheres is not None:
        hemi = {h.lower().strip() for h in options.hemispheres}
        out = out[out["hemisphere"].astype(str).str.lower().isin(hemi)]
    return out


def summarize_by_group(
    df: pd.DataFrame,
    *,
    level: str = "region",
) -> pd.DataFrame:
    """Per-group mean, std, SEM and subject count at region or rollup level."""
    if level == "region":
        work = df
        group_keys = _REGION_GROUP_KEYS
    else:
        rollup_col = _ROLLUP_COLUMN.get(level)
        if rollup_col is None:
            msg = f"Unknown rollup level {level!r}. Use division or structure."
            raise ValueError(msg)
        work = _subject_mean_rollup(df, rollup_col)
        group_keys = ["group", "channel", "hemisphere", "metric", rollup_col]

    if work.empty:
        return pd.DataFrame()

    grouped = (
        work.groupby(group_keys, dropna=False)["value"]
        .agg(n_subjects="count", mean="mean", std="std")
        .reset_index()
    )
    grouped["sem"] = grouped["std"] / np.sqrt(grouped["n_subjects"].clip(lower=1))
    grouped["level"] = level
    return grouped


def compare_groups_pair(
    df: pd.DataFrame,
    group_a: str,
    group_b: str,
    *,
    level: str = "region",
    test: str = "mannwhitney",
    fdr_alpha: float = 0.05,
) -> pd.DataFrame:
    """Pairwise comparison between two groups at each region (or rollup level)."""
    if test != "mannwhitney":
        msg = f"Unsupported test {test!r}."
        raise ValueError(msg)

    if level == "region":
        region_keys = _COMPARE_REGION_KEYS
        work = df
    else:
        rollup_col = _ROLLUP_COLUMN.get(level)
        if rollup_col is None:
            msg = f"Unknown rollup level {level!r}."
            raise ValueError(msg)
        work = _subject_mean_rollup(df, rollup_col)
        region_keys = ["channel", "hemisphere", "metric", rollup_col]

    sub_a = work[work["group"] == group_a]
    sub_b = work[work["group"] == group_b]
    if sub_a.empty or sub_b.empty:
        return pd.DataFrame()

    rows: list[dict] = []
    for key_vals, chunk_a in sub_a.groupby(region_keys, dropna=False):
        if not isinstance(key_vals, tuple):
            key_vals = (key_vals,)
        chunk_b = sub_b
        for col, val in zip(region_keys, key_vals):
            chunk_b = chunk_b[chunk_b[col] == val]
        if chunk_b.empty:
            continue

        vals_a = chunk_a["value"].to_numpy(dtype=float)
        vals_b = chunk_b["value"].to_numpy(dtype=float)
        vals_a = vals_a[np.isfinite(vals_a)]
        vals_b = vals_b[np.isfinite(vals_b)]

        row = dict(zip(region_keys, key_vals))
        row["group_a"] = group_a
        row["group_b"] = group_b
        row["n_a"] = len(vals_a)
        row["n_b"] = len(vals_b)
        row["mean_a"] = float(np.mean(vals_a)) if len(vals_a) else np.nan
        row["mean_b"] = float(np.mean(vals_b)) if len(vals_b) else np.nan
        row["median_diff"] = (
            float(np.median(vals_a) - np.median(vals_b))
            if len(vals_a) and len(vals_b)
            else np.nan
        )

        if len(vals_a) >= 1 and len(vals_b) >= 1:
            try:
                stat = mannwhitneyu(vals_a, vals_b, alternative="two-sided")
                row["statistic"] = float(stat.statistic)
                row["p_value"] = float(stat.pvalue)
            except ValueError:
                row["statistic"] = np.nan
                row["p_value"] = np.nan
        else:
            row["statistic"] = np.nan
            row["p_value"] = np.nan
        rows.append(row)

    if not rows:
        return pd.DataFrame()

    result = pd.DataFrame(rows)
    q, reject = benjamini_hochberg(result["p_value"].to_numpy(dtype=float), alpha=fdr_alpha)
    result["q_value"] = q
    result["significant"] = reject
    result["level"] = level
    result["test"] = test
    return result


def run_group_analysis(config: CohortConfig) -> GroupAnalysisResult:
    """Run the full cross-subject pipeline for a cohort configuration."""
    resolved = [resolve_cohort_sample(s) for s in config.samples]
    long = load_cohort_long(resolved)
    filtered = filter_cohort_table(long, config.group_analysis)

    summary_region = summarize_by_group(filtered, level="region")

    comparisons_frames: list[pd.DataFrame] = []
    for pair in config.group_analysis.comparisons:
        if len(pair) != 2:
            continue
        comparisons_frames.append(
            compare_groups_pair(
                filtered,
                pair[0],
                pair[1],
                level="region",
                test=config.group_analysis.test.value,
                fdr_alpha=config.group_analysis.fdr_alpha,
            )
        )
    comparisons_region = (
        pd.concat(comparisons_frames, ignore_index=True) if comparisons_frames else None
    )

    summary_division = summary_structure = None
    comparisons_division = comparisons_structure = None

    for rollup in config.group_analysis.rollups:
        level = rollup.value
        summary = summarize_by_group(filtered, level=level)
        comp_frames: list[pd.DataFrame] = []
        for pair in config.group_analysis.comparisons:
            if len(pair) != 2:
                continue
            comp_frames.append(
                compare_groups_pair(
                    filtered,
                    pair[0],
                    pair[1],
                    level=level,
                    test=config.group_analysis.test.value,
                    fdr_alpha=config.group_analysis.fdr_alpha,
                )
            )
        comparisons = pd.concat(comp_frames, ignore_index=True) if comp_frames else None

        if level == "division":
            summary_division = summary
            comparisons_division = comparisons
        elif level == "structure":
            summary_structure = summary
            comparisons_structure = comparisons

    return GroupAnalysisResult(
        cohort_long=filtered,
        summary_by_region=summary_region,
        comparisons_by_region=comparisons_region,
        summary_by_division=summary_division,
        summary_by_structure=summary_structure,
        comparisons_by_division=comparisons_division,
        comparisons_by_structure=comparisons_structure,
    )


def benjamini_hochberg(
    p_values: np.ndarray,
    *,
    alpha: float = 0.05,
) -> tuple[np.ndarray, np.ndarray]:
    """Benjamini–Hochberg FDR correction. Returns (q_values, reject_mask)."""
    p = np.asarray(p_values, dtype=float)
    n = len(p)
    if n == 0:
        return p, np.zeros(0, dtype=bool)

    q_out = np.full(n, np.nan, dtype=float)
    finite = np.isfinite(p)
    if not finite.any():
        return q_out, np.zeros(n, dtype=bool)

    p_f = p[finite]
    m = len(p_f)
    order = np.argsort(p_f)
    ranked = p_f[order]
    q = ranked * m / (np.arange(1, m + 1))
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0.0, 1.0)

    q_finite = np.empty(m, dtype=float)
    q_finite[order] = q
    q_out[finite] = q_finite
    reject = np.zeros(n, dtype=bool)
    reject[finite] = q_finite <= alpha
    return q_out, reject


def _subject_mean_rollup(df: pd.DataFrame, rollup_col: str) -> pd.DataFrame:
    """Per-subject mean within a rollup column (division or structure)."""
    keys = ["sample", "group", "channel", "hemisphere", "metric", rollup_col]
    return (
        df.groupby(keys, dropna=False)["value"]
        .mean()
        .reset_index()
    )


def _normalize_channel(value: object) -> str:
    """Normalize channel for filter matching (int 1 == str \"1\")."""
    if isinstance(value, (int, np.integer)):
        return str(int(value))
    if isinstance(value, float) and value == int(value):
        return str(int(value))
    return str(value).strip()
