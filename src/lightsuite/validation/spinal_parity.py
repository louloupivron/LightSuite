"""MATLAB vs Python parity checks for spinal cord MVP stages."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from lightsuite.registration.straightening_optimizer import run_straightening_optimizer


@dataclass
class ParityReport:
    passed: bool
    messages: list[str]


def compare_alignment_optimizer(reference_path: Path, *, rtol: float = 1e-4) -> ParityReport:
    """Compare optimizer output to a reference JSON fixture."""
    ref = json.loads(reference_path.read_text(encoding="utf-8"))
    user_cen = np.asarray(ref["user_cen"], dtype=float)
    user_ant = np.asarray(ref["user_ant"], dtype=float)
    user_pos = np.asarray(ref["user_pos"], dtype=float)
    fit = run_straightening_optimizer(user_cen, user_ant, user_pos)
    messages: list[str] = []
    passed = True
    for key in ("fit_x", "fit_y", "fit_theta"):
        if not np.allclose(fit[key], np.asarray(ref[key]), rtol=rtol, equal_nan=True):
            passed = False
            messages.append(f"{key} mismatch vs reference")
    if passed:
        messages.append("Straightening optimizer matches reference fixture.")
    return ParityReport(passed=passed, messages=messages)


def run_mvp_parity_checks(fixture_root: Path) -> ParityReport:
    """Run all automated parity checks available without MATLAB."""
    parity_dir = fixture_root / "parity"
    ref = parity_dir / "align_out_reference.json"
    if not ref.is_file():
        return ParityReport(passed=False, messages=[f"Missing reference fixture: {ref}"])
    return compare_alignment_optimizer(ref)
