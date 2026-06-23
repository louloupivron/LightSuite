"""Tests for naive unassigned-voxel registration QC."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import numpy as np
import pandas as pd
import pytest

from lightsuite.analysis.registration_qc import (
    compute_unassigned_registration_score,
    threshold_sweep,
    write_threshold_sweep_plot,
)


def test_compute_score_all_unassigned_signal() -> None:
    volume = np.array([0, 100, 200, 50], dtype=np.float32)
    labels = np.array([1, 0, 0, 1], dtype=np.int32)
    score = compute_unassigned_registration_score(volume, labels, threshold=80.0)
    assert score["image_voxels"] == 2  # 100, 200
    assert score["unassigned_image_voxels"] == 2
    assert score["naive_unassigned_fraction"] == 1.0
    assert score["naive_unassigned_percent"] == 100.0


def test_compute_score_no_unassigned_signal() -> None:
    volume = np.array([0, 100, 200], dtype=np.float32)
    labels = np.array([0, 1, 1], dtype=np.int32)
    score = compute_unassigned_registration_score(volume, labels, threshold=50.0)
    assert score["unassigned_image_voxels"] == 0
    assert score["naive_unassigned_fraction"] == 0.0


def test_compute_score_mixed() -> None:
    volume = np.array([10, 100, 100, 100], dtype=np.float32)
    labels = np.array([0, 0, 1, 1], dtype=np.int32)
    score = compute_unassigned_registration_score(volume, labels, threshold=50.0)
    assert score["image_voxels"] == 3
    assert score["unassigned_image_voxels"] == 1
    assert score["naive_unassigned_fraction"] == pytest.approx(1 / 3)


def test_compute_score_shape_mismatch_raises() -> None:
    with pytest.raises(ValueError, match="shape"):
        compute_unassigned_registration_score(
            np.zeros((2, 2), dtype=np.float32),
            np.zeros((3, 3), dtype=np.int32),
            threshold=1.0,
        )


def test_threshold_sweep_n_points() -> None:
    volume = np.ones(4, dtype=np.float32) * 100
    labels = np.array([0, 0, 1, 1], dtype=np.int32)
    sweep = threshold_sweep(volume, labels, threshold=100.0, n_points=5)
    assert len(sweep) == 5
    assert "naive_unassigned_percent" in sweep.columns


def test_write_threshold_sweep_plot(tmp_path: Path) -> None:
    sweep = pd.DataFrame(
        {
            "threshold": [50.0, 100.0, 150.0],
            "naive_unassigned_percent": [1.0, 5.0, 10.0],
        }
    )
    out = tmp_path / "sweep.png"
    write_threshold_sweep_plot(sweep, out)
    assert out.is_file()
