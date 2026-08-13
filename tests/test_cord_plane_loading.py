"""Tests for spinal cord plane-per-file parallel loading helpers."""

from __future__ import annotations

from lightsuite.io.cord_volume import _effective_cord_plane_workers


def test_effective_cord_plane_workers_respects_bounds() -> None:
    assert _effective_cord_plane_workers(4, 100) == 4
    assert _effective_cord_plane_workers(16, 100) == 8
    assert _effective_cord_plane_workers(4, 2) == 2
    assert _effective_cord_plane_workers(4, 1) == 1
    assert _effective_cord_plane_workers(0, 100) == 1
