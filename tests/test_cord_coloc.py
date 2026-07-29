"""Tests for spinal cord spot colocalization overlap."""

from __future__ import annotations

import numpy as np

from lightsuite.analysis.cord_coloc import (
    compute_pairwise_coloc_summary,
    match_points_within_tolerance,
)


def test_match_points_within_tolerance_finds_neighbor() -> None:
    source = np.array([[1.0, 1.0, 1.0], [10.0, 10.0, 10.0]])
    target = np.array([[1.2, 1.1, 1.0]])
    n_match, mask = match_points_within_tolerance(source, target, tolerance_voxels=2.0)
    assert n_match == 1
    assert mask.tolist() == [True, False]


def test_compute_pairwise_coloc_summary_fractions() -> None:
    points = {
        "a": np.array([[1.0, 1.0, 1.0], [20.0, 20.0, 20.0]]),
        "b": np.array([[1.1, 1.0, 1.0], [9.0, 9.0, 9.0]]),
    }
    summary = compute_pairwise_coloc_summary(points, tolerance_voxels=2.0)
    assert len(summary) == 1
    row = summary.iloc[0]
    assert row["n_overlap_source_to_target"] == 1
    assert row["frac_of_source"] == 0.5


def test_compute_pairwise_coloc_summary_triple_row() -> None:
    points = {
        "a": np.array([[1.0, 1.0, 1.0], [5.0, 5.0, 5.0]]),
        "b": np.array([[1.1, 1.0, 1.0]]),
        "c": np.array([[1.0, 1.2, 1.0]]),
    }
    summary = compute_pairwise_coloc_summary(points, tolerance_voxels=2.0)
    triple = summary[summary["comparison"] == "triple"]
    assert len(triple) == 1
    assert triple.iloc[0]["n_overlap_source_to_target"] == 1
