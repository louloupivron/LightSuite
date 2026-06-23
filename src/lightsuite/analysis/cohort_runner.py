"""Orchestrate cross-subject group analysis and write cohort outputs."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd
from rich.console import Console

from lightsuite.analysis.cohort_models import CohortConfig
from lightsuite.analysis.group import GroupAnalysisResult, run_group_analysis

console = Console()


@dataclass
class CohortRunResult:
    output_dir: Path
    cohort_long_path: Path | None = None
    summary_region_path: Path | None = None
    comparisons_region_path: Path | None = None
    summary_division_path: Path | None = None
    comparisons_division_path: Path | None = None
    summary_structure_path: Path | None = None
    comparisons_structure_path: Path | None = None
    n_subjects: int = 0
    n_groups: int = 0
    written_paths: list[Path] = field(default_factory=list)


def _write_csv(path: Path, df: pd.DataFrame) -> Path | None:
    if df is None or df.empty:
        return None
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)
    return path


def run_cohort_group_analysis(config: CohortConfig) -> CohortRunResult:
    """Run cross-subject summaries and optional pairwise comparisons."""
    result = run_group_analysis(config)
    out_dir = config.output_dir.expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    written: list[Path] = []
    paths: dict[str, Path | None] = {}

    spec = [
        ("cohort_long", result.cohort_long, "cohort_region_stats_long.csv"),
        ("summary_region", result.summary_by_region, "group_summary_by_region.csv"),
        ("comparisons_region", result.comparisons_by_region, "group_comparisons_by_region.csv"),
        ("summary_division", result.summary_by_division, "group_summary_by_division.csv"),
        ("comparisons_division", result.comparisons_by_division, "group_comparisons_by_division.csv"),
        ("summary_structure", result.summary_by_structure, "group_summary_by_structure.csv"),
        ("comparisons_structure", result.comparisons_by_structure, "group_comparisons_by_structure.csv"),
    ]

    for key, frame, filename in spec:
        path = _write_csv(out_dir / filename, frame)
        paths[key] = path
        if path is not None:
            written.append(path)

    n_subjects = result.cohort_long["sample"].nunique() if len(result.cohort_long) else 0
    n_groups = result.cohort_long["group"].nunique() if len(result.cohort_long) else 0

    console.print(
        f"[green]Cohort analysis '{config.name}':[/green] "
        f"{n_subjects} subject(s), {n_groups} group(s), {len(written)} file(s) → {out_dir}"
    )
    if result.comparisons_by_region is not None and len(result.comparisons_by_region):
        n_sig = int(result.comparisons_by_region["significant"].sum())
        console.print(f"  Region comparisons: {len(result.comparisons_by_region)} tests, {n_sig} significant (FDR)")

    return CohortRunResult(
        output_dir=out_dir,
        cohort_long_path=paths["cohort_long"],
        summary_region_path=paths["summary_region"],
        comparisons_region_path=paths["comparisons_region"],
        summary_division_path=paths["summary_division"],
        comparisons_division_path=paths["comparisons_division"],
        summary_structure_path=paths["summary_structure"],
        comparisons_structure_path=paths["comparisons_structure"],
        n_subjects=n_subjects,
        n_groups=n_groups,
        written_paths=written,
    )
