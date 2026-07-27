"""Tests for multiresolution QA plots."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from lightsuite.multires.plots import save_registration_overlay_qc_plot


def test_save_registration_overlay_qc_plot(tmp_path: Path) -> None:
    overview = np.zeros((4, 8, 8), dtype=np.float32)
    roi = np.zeros((4, 8, 8), dtype=np.float32)
    overview[2, 3:6, 3:6] = 100.0
    roi[2, 3:6, 3:6] = 80.0

    out = tmp_path / "registration_overlay_qc.png"
    ncc = save_registration_overlay_qc_plot(
        overview_crop=overview,
        registered_roi=roi,
        output_path=out,
        pair_label="cortex_9x_561",
        z_index=2,
    )
    assert out.is_file()
    assert 0.0 < ncc <= 1.0


def test_registration_overlay_qc_requires_matching_shapes() -> None:
    overview = np.zeros((4, 8, 8), dtype=np.float32)
    roi = np.zeros((4, 10, 10), dtype=np.float32)
    try:
        save_registration_overlay_qc_plot(
            overview_crop=overview,
            registered_roi=roi,
            output_path=Path("/tmp/should_not_write.png"),
        )
    except ValueError as exc:
        assert "shapes must match" in str(exc)
    else:
        raise AssertionError("expected ValueError for mismatched shapes")
