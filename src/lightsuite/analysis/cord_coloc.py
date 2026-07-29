"""Pairwise colocalization overlap between imported spinal cord spot labels."""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

from lightsuite.analysis.counts import ATLAS_POINTS_KEY, load_atlas_points
from lightsuite.config.models import SpinalCordPipelineConfig
from lightsuite.export.cord_registered import volume_registered_dir


@dataclass(frozen=True)
class CordColocOverlapResult:
    summary: pd.DataFrame
    summary_path: Path | None = None
    plot_path: Path | None = None


def _slug_label(label: str) -> str:
    return str(label).strip()


def load_cord_coloc_points(
    register_path: Path,
    labels: list[str],
) -> dict[str, np.ndarray]:
    """Load atlas-space point clouds for the requested import labels."""
    register_path = register_path.expanduser().resolve()
    if not register_path.is_dir():
        msg = f"Missing {register_path}. Run 'lightsuite spinal import-annotations' first."
        raise FileNotFoundError(msg)

    points_by_label: dict[str, np.ndarray] = {}
    missing: list[str] = []
    for label in labels:
        slug = _slug_label(label)
        npz_path = register_path / f"{slug}_atlas_coords.npz"
        if not npz_path.is_file():
            missing.append(slug)
            continue
        points_by_label[slug] = load_atlas_points(npz_path, key=ATLAS_POINTS_KEY)

    if missing:
        msg = f"Missing atlas coord files for label(s): {missing}"
        raise FileNotFoundError(msg)
    return points_by_label


def match_points_within_tolerance(
    source: np.ndarray,
    target: np.ndarray,
    *,
    tolerance_voxels: float,
) -> tuple[int, np.ndarray]:
    """Count source points with a target neighbor within *tolerance_voxels*."""
    src = np.asarray(source, dtype=float)
    if src.ndim != 2 or src.shape[1] < 3:
        msg = f"Expected Nx3+ points, got shape {src.shape}"
        raise ValueError(msg)
    if src.size == 0:
        return 0, np.zeros(0, dtype=bool)
    tgt = np.asarray(target, dtype=float)
    if tgt.size == 0:
        return 0, np.zeros(len(src), dtype=bool)

    tree = cKDTree(tgt[:, :3])
    dists, _ = tree.query(src[:, :3], distance_upper_bound=float(tolerance_voxels))
    matched = np.isfinite(dists) & (dists <= float(tolerance_voxels))
    return int(matched.sum()), matched


def compute_pairwise_coloc_summary(
    points_by_label: dict[str, np.ndarray],
    *,
    tolerance_voxels: float = 2.0,
) -> pd.DataFrame:
    """Pairwise overlap counts and fractions between import labels."""
    labels = list(points_by_label)
    records: list[dict[str, object]] = []

    for source_label, target_label in combinations(labels, 2):
        source_pts = points_by_label[source_label]
        target_pts = points_by_label[target_label]
        n_overlap_st, mask_st = match_points_within_tolerance(
            source_pts, target_pts, tolerance_voxels=tolerance_voxels
        )
        n_overlap_ts, mask_ts = match_points_within_tolerance(
            target_pts, source_pts, tolerance_voxels=tolerance_voxels
        )
        records.append(
            {
                "source": source_label,
                "target": target_label,
                "n_source": int(len(source_pts)),
                "n_target": int(len(target_pts)),
                "n_overlap_source_to_target": n_overlap_st,
                "n_overlap_target_to_source": n_overlap_ts,
                "frac_of_source": n_overlap_st / len(source_pts) if len(source_pts) else np.nan,
                "frac_of_target": n_overlap_ts / len(target_pts) if len(target_pts) else np.nan,
                "tolerance_voxels": float(tolerance_voxels),
                "comparison": "pairwise",
            }
        )

    if len(labels) == 3:
        a, b, c = labels
        n_triple = _triple_overlap_count(
            points_by_label[a],
            points_by_label[b],
            points_by_label[c],
            tolerance_voxels=tolerance_voxels,
        )
        records.append(
            {
                "source": a,
                "target": f"{b} & {c}",
                "n_source": int(len(points_by_label[a])),
                "n_target": int(len(points_by_label[b]) + len(points_by_label[c])),
                "n_overlap_source_to_target": n_triple,
                "n_overlap_target_to_source": n_triple,
                "frac_of_source": n_triple / len(points_by_label[a]) if len(points_by_label[a]) else np.nan,
                "frac_of_target": np.nan,
                "tolerance_voxels": float(tolerance_voxels),
                "comparison": "triple",
            }
        )

    frame = pd.DataFrame.from_records(records)
    return frame


def _triple_overlap_count(
    a: np.ndarray,
    b: np.ndarray,
    c: np.ndarray,
    *,
    tolerance_voxels: float,
) -> int:
    """Count points in *a* with neighbors in both *b* and *c*."""
    if len(a) == 0:
        return 0
    _, mask_b = match_points_within_tolerance(a, b, tolerance_voxels=tolerance_voxels)
    _, mask_c = match_points_within_tolerance(a, c, tolerance_voxels=tolerance_voxels)
    return int(np.sum(mask_b & mask_c))


def run_cord_coloc_overlap(
    config: SpinalCordPipelineConfig,
    *,
    labels: list[str] | None = None,
    tolerance_voxels: float = 2.0,
    output_csv: Path | None = None,
) -> CordColocOverlapResult:
    """Compute colocalization overlap table for imported spot labels."""
    register_path = volume_registered_dir(config)
    label_list = labels if labels is not None else config.analysis.point_labels
    if not label_list:
        msg = "No import labels provided and analysis.point_labels is empty."
        raise ValueError(msg)

    points = load_cord_coloc_points(register_path, [str(label) for label in label_list])
    summary = compute_pairwise_coloc_summary(points, tolerance_voxels=tolerance_voxels)

    summary_path: Path | None = None
    if output_csv is not None:
        summary_path = output_csv.expanduser()
        summary_path.parent.mkdir(parents=True, exist_ok=True)
        summary.to_csv(summary_path, index=False)

    return CordColocOverlapResult(summary=summary, summary_path=summary_path)


__all__ = [
    "CordColocOverlapResult",
    "compute_pairwise_coloc_summary",
    "load_cord_coloc_points",
    "match_points_within_tolerance",
    "run_cord_coloc_overlap",
]
