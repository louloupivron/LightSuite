"""Tests for doctor checks."""

from __future__ import annotations

from lightsuite.cli.doctor import _check_python, run_doctor


def test_python_version_check() -> None:
    result = _check_python()
    assert result.ok is True


def test_doctor_runs_without_config() -> None:
    report = run_doctor(config=None, strict=False)
    assert any(r.name == "Python" for r in report.results)
    assert any(r.name == "elastix binary" for r in report.results)
