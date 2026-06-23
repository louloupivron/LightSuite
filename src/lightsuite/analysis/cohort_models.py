"""Pydantic models for cross-subject cohort analysis."""

from __future__ import annotations

from enum import Enum
from pathlib import Path
from typing import Annotated

from pydantic import BaseModel, Field, field_validator, model_validator


class GroupTest(str, Enum):
    MANNWHITNEY = "mannwhitney"


class RollupLevel(str, Enum):
    DIVISION = "division"
    STRUCTURE = "structure"


class CohortSampleEntry(BaseModel):
    """One subject in a cohort, with an experimental group label."""

    group: str = Field(min_length=1, description="Experimental group (e.g. control, treatment).")
    id: str | None = Field(
        default=None,
        description="Subject id; defaults to sample name from brain config or CSV.",
    )
    region_stats: Path | None = Field(
        default=None,
        description="Path to volume_registered/region_stats.csv for this subject.",
    )
    config: Path | None = Field(
        default=None,
        description="Brain pipeline YAML; region_stats is resolved from sample.save_path.",
    )

    @field_validator("region_stats", "config")
    @classmethod
    def expand_paths(cls, value: Path | None) -> Path | None:
        if value is None:
            return None
        return value.expanduser()

    @model_validator(mode="after")
    def require_stats_source(self) -> CohortSampleEntry:
        if self.region_stats is None and self.config is None:
            msg = "Each cohort sample needs region_stats or config."
            raise ValueError(msg)
        return self


class GroupAnalysisConfig(BaseModel):
    """Filters and statistics for cross-subject group comparison."""

    channels: list[int | str] | None = Field(
        default=None,
        description="Restrict to these channel ids (int) or point labels (str); null = all.",
    )
    metrics: list[str] | None = Field(
        default=None,
        description="Restrict to these metrics; null = all present in the data.",
    )
    hemispheres: list[str] | None = Field(
        default=None,
        description='Restrict to "right" / "left"; null = both.',
    )
    rollups: list[RollupLevel] = Field(
        default_factory=list,
        description="Also emit summaries rolled up to division or structure.",
    )
    comparisons: list[Annotated[list[str], Field(min_length=2, max_length=2)]] = Field(
        default_factory=list,
        description='Pairwise group pairs, e.g. [["control", "treatment"]].',
    )
    test: GroupTest = GroupTest.MANNWHITNEY
    fdr_alpha: float = Field(default=0.05, gt=0, le=1)


class CohortConfig(BaseModel):
    """Multi-subject cohort for cross-group region statistics."""

    name: str = Field(min_length=1)
    output_dir: Path
    samples: Annotated[list[CohortSampleEntry], Field(min_length=2)]
    group_analysis: GroupAnalysisConfig = Field(default_factory=GroupAnalysisConfig)

    @field_validator("output_dir")
    @classmethod
    def expand_output_dir(cls, value: Path) -> Path:
        return value.expanduser()
