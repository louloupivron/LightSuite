"""Tests for spinal cord longitudinal alignment helpers."""

from __future__ import annotations

import numpy as np
import pytest

from lightsuite.gui.chooselist import generate_cord_longitudinal_list
from lightsuite.gui.slice_correspondence import SliceAnchor, SliceCorrespondence
from lightsuite.registration.cord_affine import build_cord_z_transinit
from lightsuite.registration.cord_longitudinal import (
    CORD_LONGITUDINAL_AXIS,
    build_cord_z_transinit_from_correspondence,
    build_longitudinal_anchors,
    estimate_cord_atlas_plane,
    resolve_cord_z_transinit,
)


def test_generate_cord_longitudinal_list_is_ordered() -> None:
    chooselist = generate_cord_longitudinal_list(922, n_slices=20)
    assert chooselist.shape == (20, 4)
    assert np.all(chooselist[:, 1] == 3)
    assert np.all(np.diff(chooselist[:, 0]) >= 0)


def test_build_cord_z_transinit_from_correspondence_fits_line() -> None:
    corr = SliceCorrespondence.single_axis(
        CORD_LONGITUDINAL_AXIS,
        np.eye(4).tolist(),
        [
            SliceAnchor(sample_index=10, atlas_plane=100, confirmed=True),
            SliceAnchor(sample_index=510, atlas_plane=600, confirmed=True),
        ],
    )
    transinit = build_cord_z_transinit_from_correspondence(corr, atlas_depth=1000)
    assert transinit[2, 2] == pytest.approx(1.0, abs=0.05)
    assert transinit[2, 3] == pytest.approx(-90.0, abs=2.0)
    assert transinit[2, 2] * 99 + transinit[2, 3] == pytest.approx(9.0, abs=2.0)


def test_resolve_cord_z_transinit_uses_correspondence_when_confirmed() -> None:
    centered = build_cord_z_transinit(500, 1000)
    corr = SliceCorrespondence.single_axis(
        CORD_LONGITUDINAL_AXIS,
        np.eye(4).tolist(),
        [
            SliceAnchor(sample_index=1, atlas_plane=200, confirmed=True),
            SliceAnchor(sample_index=500, atlas_plane=900, confirmed=True),
        ],
    )
    fitted = resolve_cord_z_transinit(500, 1000, corr)
    assert not np.allclose(fitted, centered)
    assert fitted[2, 2] == pytest.approx(0.75, abs=0.05)


def test_resolve_cord_z_transinit_falls_back_without_confirmed() -> None:
    corr = SliceCorrespondence.single_axis(
        CORD_LONGITUDINAL_AXIS,
        np.eye(4).tolist(),
        [SliceAnchor(sample_index=10, atlas_plane=20, confirmed=False)],
    )
    centered = build_cord_z_transinit(400, 800)
    assert np.allclose(resolve_cord_z_transinit(400, 800, corr), centered)
    assert np.allclose(resolve_cord_z_transinit(400, 800, None), centered)


def test_estimate_cord_atlas_plane_matches_centered_init() -> None:
    nslices, atlas_depth = 922, 1567
    transinit = build_cord_z_transinit(nslices, atlas_depth)
    sample_index = 100
    atlas_plane = estimate_cord_atlas_plane(
        sample_index,
        nslices=nslices,
        atlas_depth=atlas_depth,
    )
    expected = (sample_index - 1 - transinit[2, 3]) / transinit[2, 2] + 1
    assert atlas_plane == int(np.round(expected))


def test_build_longitudinal_anchors_seeds_centered_planes() -> None:
    chooselist = generate_cord_longitudinal_list(100, n_slices=5)
    anchors = build_longitudinal_anchors(
        chooselist,
        nslices=100,
        atlas_depth=200,
        confirmed=False,
    )
    assert len(anchors) == 5
    assert all(anchor.atlas_plane > 0 for anchor in anchors)
    assert all(not anchor.confirmed for anchor in anchors)
