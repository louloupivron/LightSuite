"""Naive registration QC from high-signal voxels in unassigned atlas divisions.

Port of the ``unassigned_registration_score`` notebook: if registration is imperfect,
tissue signal can appear in division ``id == 0``. The fraction of above-threshold
voxels landing in unassigned divisions is a simple alignment warning metric.

**Limitations** (see notebook): global threshold sensitivity, boundary biology, artifacts.
This is a proxy — not a geometric registration error measure.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

UNASSIGNED_DIVISION_ID = 0

SCORE_COLUMNS = (
    "threshold",
    "total_voxels",
    "image_voxels",
    "unassigned_voxels",
    "unassigned_image_voxels",
    "naive_unassigned_fraction",
    "naive_unassigned_percent",
)


@dataclass(frozen=True)
class RegistrationQcResult:
    channel: int
    threshold: float
    score: dict[str, float | int]
    sweep: pd.DataFrame | None = None
    score_csv: Path | None = None
    sweep_csv: Path | None = None
    sweep_plot: Path | None = None


def compute_unassigned_registration_score(
    volume: np.ndarray,
    division_labels: np.ndarray,
    threshold: float,
    *,
    unassigned_id: int = UNASSIGNED_DIVISION_ID,
) -> dict[str, float | int]:
    """Fraction of above-threshold voxels that fall in unassigned division labels."""
    vol = np.asarray(volume, dtype=np.float32)
    labels = np.asarray(division_labels, dtype=np.int32)
    if vol.shape != labels.shape:
        msg = f"Volume shape {vol.shape} != division labels {labels.shape}"
        raise ValueError(msg)

    image_mask = vol >= float(threshold)
    unassigned_mask = labels == int(unassigned_id)

    image_voxels = int(image_mask.sum())
    unassigned_image_voxels = int((image_mask & unassigned_mask).sum())

    if image_voxels == 0:
        fraction = np.nan
    else:
        fraction = float(unassigned_image_voxels / image_voxels)

    return {
        "threshold": float(threshold),
        "total_voxels": int(vol.size),
        "image_voxels": image_voxels,
        "unassigned_voxels": int(unassigned_mask.sum()),
        "unassigned_image_voxels": unassigned_image_voxels,
        "naive_unassigned_fraction": float(fraction) if np.isfinite(fraction) else np.nan,
        "naive_unassigned_percent": float(fraction * 100.0) if np.isfinite(fraction) else np.nan,
    }


def threshold_sweep(
    volume: np.ndarray,
    division_labels: np.ndarray,
    threshold: float,
    *,
    n_points: int = 9,
    span: tuple[float, float] = (0.5, 1.5),
    unassigned_id: int = UNASSIGNED_DIVISION_ID,
) -> pd.DataFrame:
    """Evaluate the score at ``n_points`` thresholds around the reference value."""
    if n_points < 2:
        msg = "n_points must be >= 2 for a sweep."
        raise ValueError(msg)
    lo_scale, hi_scale = span
    thresholds = np.linspace(
        max(0.0, float(threshold) * lo_scale),
        float(threshold) * hi_scale,
        int(n_points),
    )
    rows = [
        compute_unassigned_registration_score(
            volume, division_labels, float(t), unassigned_id=unassigned_id
        )
        for t in thresholds
    ]
    return pd.DataFrame(rows)


def write_threshold_sweep_plot(sweep: pd.DataFrame, path: Path) -> Path:
    """Save a simple threshold-sensitivity line plot."""
    import matplotlib.pyplot as plt

    path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(
        sweep["threshold"],
        sweep["naive_unassigned_percent"],
        marker="o",
    )
    ax.set_xlabel("threshold")
    ax.set_ylabel("naive unassigned %")
    ax.set_title("Threshold sensitivity of unassigned registration score")
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return path


def score_to_dataframe(score: dict[str, float | int], *, extra: dict | None = None) -> pd.DataFrame:
    """Single-row score table with optional metadata columns prepended."""
    row = dict(extra or {})
    row.update(score)
    return pd.DataFrame([row])
