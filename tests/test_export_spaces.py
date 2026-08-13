"""Tests for export space CLI/GUI helpers."""

from __future__ import annotations

import pytest

from lightsuite.cli.spaces import (
    default_export_space_checks,
    export_spaces_from_checks,
    format_export_spaces,
    parse_spaces_option,
)


def test_export_spaces_from_checks_both() -> None:
    assert export_spaces_from_checks(atlas=True, sample=True) == ["atlas", "sample"]


def test_export_spaces_from_checks_requires_one() -> None:
    with pytest.raises(ValueError, match="at least one"):
        export_spaces_from_checks(atlas=False, sample=False)


def test_default_export_space_checks() -> None:
    assert default_export_space_checks(["atlas"]) == (True, False)
    assert default_export_space_checks(["sample"]) == (False, True)
    assert default_export_space_checks(["atlas", "sample"]) == (True, True)
    assert default_export_space_checks(None) == (True, False)


def test_format_export_spaces() -> None:
    assert format_export_spaces(["atlas", "sample"]) == "atlas + sample"
    assert format_export_spaces(["sample"]) == "sample"


def test_parse_spaces_option_both() -> None:
    assert parse_spaces_option("both") == ["atlas", "sample"]
