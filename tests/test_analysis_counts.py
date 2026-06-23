"""Tests for per-region cell counting from atlas-space points."""

from __future__ import annotations

import numpy as np

from lightsuite.analysis.counts import atlas_points_to_voxel_indices, count_points_in_regions


def _annotation() -> np.ndarray:
    av = np.zeros((4, 4, 4), dtype=np.int32)
    av[:, :, :2] = 7  # right side (z < 2)
    av[:, :, 2:] = 9  # left side  (z >= 2)
    return av


def test_atlas_points_to_voxel_indices_is_1based_and_bounded() -> None:
    pts = np.array([[1.0, 1.0, 1.0], [4.0, 4.0, 4.0], [100.0, 1.0, 1.0]])
    idx = atlas_points_to_voxel_indices(pts, (4, 4, 4))
    # third point is out of bounds and dropped
    assert idx.shape == (2, 3)
    np.testing.assert_array_equal(idx[0], [0, 0, 0])
    np.testing.assert_array_equal(idx[1], [3, 3, 3])


def test_count_points_in_regions_counts_and_density() -> None:
    av = _annotation()
    # one point in the right half (label 7), one in the left half (label 9)
    points = np.array([[1.0, 1.0, 1.0], [4.0, 4.0, 4.0]])
    tidy = count_points_in_regions(
        points,
        av,
        atlas_id="allen",
        atlas_resolution_um=10.0,
        sample="m1",
        channel="cells",
    )

    counts = tidy[tidy["metric"] == "cell_count"]
    right7 = counts[(counts["parcellation_index"] == 7) & (counts["hemisphere"] == "right")]
    left9 = counts[(counts["parcellation_index"] == 9) & (counts["hemisphere"] == "left")]
    assert right7["value"].iloc[0] == 1.0
    assert left9["value"].iloc[0] == 1.0

    # region 7 fills the entire right side: 32 voxels × (10µm)³ = 32e-6 mm³
    density = tidy[
        (tidy["metric"] == "cell_density")
        & (tidy["parcellation_index"] == 7)
        & (tidy["hemisphere"] == "right")
    ]
    assert density["value"].iloc[0] == 1.0 / (32 * (10e-3) ** 3)


def test_count_points_empty_returns_empty_frame() -> None:
    tidy = count_points_in_regions(
        np.zeros((0, 3)),
        _annotation(),
        atlas_id="allen",
        atlas_resolution_um=10.0,
        sample="m1",
        channel="cells",
    )
    assert len(tidy) == 0
