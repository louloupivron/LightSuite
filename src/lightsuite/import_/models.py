"""Data models for external annotation import."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np


@dataclass
class ImportedPoints:
    """Cell coordinates in LightSuite native sample space."""

    label: str
    coordinates: np.ndarray
    features: np.ndarray | None = None
    source_path: Path | None = None
    metadata: dict = field(default_factory=dict)


@dataclass
class ImportedMask:
    """Binary or label mask on a known voxel grid (Y, X, Z)."""

    label: str
    volume: np.ndarray
    voxel_um: list[float]
    source_path: Path | None = None
    metadata: dict = field(default_factory=dict)


@dataclass
class AnnotationImportResult:
    """Outputs written under volume_registered/."""

    label: str
    kind: str
    atlas_points_path: Path | None = None
    atlas_mask_path: Path | None = None
    atlas_csv_path: Path | None = None
    sample_points_path: Path | None = None
    sample_mask_path: Path | None = None
    n_input: int = 0
    n_atlas: int = 0
    n_sample: int = 0
