"""Tests for Fiederling spinal cord hemisphere splitting."""

from __future__ import annotations

import numpy as np
import pandas as pd

from lightsuite.analysis.cord_counts import count_points_in_cord_regions
from lightsuite.analysis.cord_hemisphere import (
    cord_hemisphere_side_volume,
    hemisphere_label_from_side,
    sample_space_hemisphere_flip,
)
from lightsuite.analysis.cord_parcellation import parcellate_cord_intensities
from lightsuite.analysis.counts import atlas_points_to_voxel_indices


def _segments() -> pd.DataFrame:
    return pd.DataFrame({"Segment": ["C1", "C2"], "Start": [1, 3], "End": [2, 4]})


def _annotation_and_hemisphere() -> tuple[np.ndarray, np.ndarray]:
    ann = np.zeros((4, 4, 4), dtype=np.int32)
    ann[:2, :, 1] = 7
    ann[2:, :, 1] = 7
    ann[:, :, 2] = 9
    ann[:, :, 3] = 9
    hem = np.zeros((4, 4, 4), dtype=np.uint8)
    hem[:2, :, :] = 255
    hem[2:, :, :] = 0
    return ann, hem


def test_cord_hemisphere_side_volume_labels() -> None:
    ann, hem = _annotation_and_hemisphere()
    side = cord_hemisphere_side_volume(hem, ann, flip=False)
    assert side[0, 0, 1] == 0
    assert side[3, 0, 1] == 1
    flipped = cord_hemisphere_side_volume(hem, ann, flip=True)
    assert flipped[0, 0, 1] == 1
    assert flipped[3, 0, 1] == 0


def test_parcellate_split_hemispheres() -> None:
    ann, hem = _annotation_and_hemisphere()
    intensity = np.ones((4, 4, 4), dtype=np.float32)
    intensity[:2, :, 1] = 10.0
    intensity[2:, :, 1] = 20.0
    tidy = parcellate_cord_intensities(
        intensity,
        ann,
        _segments(),
        sample="op87",
        channel=1,
        hemisphere_side=cord_hemisphere_side_volume(hem, ann),
        split_hemispheres=True,
    )
    medians = tidy[tidy["metric"] == "median_intensity"]
    assert set(medians["hemisphere"]) == {"left", "right"}
    right = medians[(medians["hemisphere"] == "right") & (medians["parcellation_index"] == 7)]
    left = medians[(medians["hemisphere"] == "left") & (medians["parcellation_index"] == 7)]
    assert right["value"].iloc[0] == 10.0
    assert left["value"].iloc[0] == 20.0


def test_count_points_split_hemispheres() -> None:
    ann, hem = _annotation_and_hemisphere()
    points = np.array(
        [
            [1.0, 1.0, 2.0],
            [1.0, 1.0, 2.0],
            [3.0, 3.0, 4.0],
            [4.0, 4.0, 4.0],
        ]
    )
    tidy = count_points_in_cord_regions(
        points,
        ann,
        _segments(),
        sample="op87",
        channel="cells",
        hemisphere_side=cord_hemisphere_side_volume(hem, ann),
        split_hemispheres=True,
    )
    counts = tidy[tidy["metric"] == "cell_count"]
    right = counts[(counts["hemisphere"] == "right") & (counts["parcellation_index"] == 7)]
    left = counts[(counts["hemisphere"] == "left") & (counts["parcellation_index"] == 9)]
    assert right["value"].iloc[0] == 2.0
    assert left["value"].iloc[0] == 2.0


def test_sample_space_hemisphere_flip_complements_atlas_flip() -> None:
    assert sample_space_hemisphere_flip(atlas_hemisphere_flip=False) is True
    assert sample_space_hemisphere_flip(atlas_hemisphere_flip=True) is False


def test_sample_space_hemisphere_matches_atlas_after_lateral_mirror() -> None:
    """Warped sample masks are mirrored on Y; complementary flip restores atlas labels."""
    ann = np.ones((4, 4, 1), dtype=np.int32) * 7
    hem_atlas = np.zeros((4, 4, 1), dtype=np.uint8)
    hem_atlas[:2, :, 0] = 255
    hem_sample = np.flip(hem_atlas, axis=0)

    side_atlas = cord_hemisphere_side_volume(hem_atlas, ann, flip=False)
    side_sample = cord_hemisphere_side_volume(
        hem_sample,
        ann,
        flip=sample_space_hemisphere_flip(atlas_hemisphere_flip=False),
    )

    points = np.array([[1.0, 1.0, 1.0], [1.0, 3.0, 1.0], [1.0, 4.0, 1.0]])

    def labels(pts: np.ndarray, side: np.ndarray) -> list[str | None]:
        idx = atlas_points_to_voxel_indices(pts, ann.shape)
        return [hemisphere_label_from_side(int(side[y, x, z])) for y, x, z in idx]

    assert labels(points, side_atlas) == labels(points, side_sample)
