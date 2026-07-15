"""Tests for spinal cord longitudinal orientation helpers."""

from __future__ import annotations

from pathlib import Path

import pytest

from lightsuite.registration.cord_orientation import (
    CAUDOROSTRAL,
    ROSTROCAUDAL,
    cord_orientation_path,
    direction_from_tofliprc,
    load_cord_orientation,
    normalize_direction,
    resolve_cord_orientation,
    save_cord_orientation,
    tofliprc_from_direction,
)


def test_tofliprc_from_direction() -> None:
    assert tofliprc_from_direction(ROSTROCAUDAL) is False
    assert tofliprc_from_direction(CAUDOROSTRAL) is True
    assert direction_from_tofliprc(True) == CAUDOROSTRAL
    assert direction_from_tofliprc(False) == ROSTROCAUDAL


def test_normalize_direction_is_case_insensitive() -> None:
    assert normalize_direction("  CaudoRostral  ") == CAUDOROSTRAL
    with pytest.raises(ValueError, match="Longitudinal direction"):
        normalize_direction("sideways")


def test_save_and_load_roundtrip(tmp_path: Path) -> None:
    path = save_cord_orientation(tmp_path, CAUDOROSTRAL, source="manual")
    assert path == cord_orientation_path(tmp_path)
    text = path.read_text(encoding="utf-8")
    assert "direction: caudorostral" in text
    assert "tofliprc: true" in text
    assert load_cord_orientation(tmp_path) == CAUDOROSTRAL


def test_load_returns_none_when_absent(tmp_path: Path) -> None:
    assert load_cord_orientation(tmp_path) is None


def test_load_falls_back_to_tofliprc_field(tmp_path: Path) -> None:
    cord_orientation_path(tmp_path).write_text("tofliprc: true\n", encoding="utf-8")
    assert load_cord_orientation(tmp_path) == CAUDOROSTRAL


def test_resolve_prefers_config_over_file(tmp_path: Path) -> None:
    save_cord_orientation(tmp_path, ROSTROCAUDAL)
    resolved = resolve_cord_orientation(
        tmp_path,
        config_direction=CAUDOROSTRAL,
        require=True,
    )
    assert resolved == CAUDOROSTRAL


def test_resolve_reads_file_when_no_config(tmp_path: Path) -> None:
    save_cord_orientation(tmp_path, CAUDOROSTRAL)
    assert resolve_cord_orientation(tmp_path, require=True) == CAUDOROSTRAL


def test_resolve_requires_orientation(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError, match="cord_orientation.txt"):
        resolve_cord_orientation(tmp_path, require=True)


def test_resolve_defaults_when_not_required(tmp_path: Path) -> None:
    assert resolve_cord_orientation(tmp_path, require=False) == ROSTROCAUDAL
