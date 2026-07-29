"""Tests for Fiederling spinal cord hemisphere splitting."""

from __future__ import annotations

import numpy as np
import pandas as pd

from lightsuite.analysis.cord_counts import count_points_in_cord_regions
from lightsuite.analysis.cord_hemisphere import cord_hemisphere_side_volume
from lightsuite.analysis.cord_parcellation import parcellate_cord_intensities


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
