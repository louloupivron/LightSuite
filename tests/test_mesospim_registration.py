"""Tests for mesoSPIM registration helpers."""

from __future__ import annotations

from lightsuite.mesospim.registration import sanitize_experiment_name


def test_sanitize_experiment_name() -> None:
    assert sanitize_experiment_name("baseline geom k2") == "baseline_geom_k2"
    assert sanitize_experiment_name("  run-01  ") == "run-01"
