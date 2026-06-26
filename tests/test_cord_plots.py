"""Tests for spinal cord QC plot helpers."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from lightsuite.registration.cord_plots import save_cord_annotation_preview


def test_save_cord_annotation_preview_requires_matching_shapes(tmp_path: Path) -> None:
    volume = np.zeros((10, 8, 12), dtype=np.uint8)
    annotation = np.zeros((10, 8, 11), dtype=np.uint16)
    with pytest.raises(ValueError, match="shapes must match"):
        save_cord_annotation_preview(volume, annotation, tmp_path / "bad.png")


def test_save_cord_annotation_preview_writes_png(tmp_path: Path) -> None:
    volume = np.zeros((40, 20, 30), dtype=np.uint8)
    annotation = np.zeros((40, 20, 30), dtype=np.uint16)
    volume[15:25, 8:12, 10:20] = 200
    annotation[16:24, 7:13, 12:18] = 5
    out = save_cord_annotation_preview(volume, annotation, tmp_path / "qc.png")
    assert out.is_file()
    assert out.stat().st_size > 0
