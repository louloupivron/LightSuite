"""Tests for spinal cord QC plot helpers."""

from __future__ import annotations

import threading
import warnings
from pathlib import Path

import matplotlib
import numpy as np
import pytest

from lightsuite.registration.cord_plots import _select_axial_slices, save_cord_annotation_preview


def test_save_cord_annotation_preview_requires_matching_shapes(tmp_path: Path) -> None:
    volume = np.zeros((10, 8, 12), dtype=np.uint8)
    annotation = np.zeros((10, 8, 11), dtype=np.uint16)
    with pytest.raises(ValueError, match="shapes must match"):
        save_cord_annotation_preview(volume, annotation, tmp_path / "bad.png")


def test_select_axial_slices_spans_tissue_extent() -> None:
    volume = np.zeros((40, 20, 921), dtype=np.uint8)
    volume[15:25, 8:12, 300:622] = 200
    # Sparse tail tissue beyond the main cord block (no annotation there in real QC).
    volume[18:22, 9:11, 900:921] = 50

    slices = _select_axial_slices(volume, n_show=12, bg_threshold=2)

    assert slices.shape == (12,)
    assert int(slices[0]) >= 300
    assert int(slices[-1]) < 700
    assert np.all(np.diff(slices) >= 0)


def test_select_axial_slices_independent_of_annotation() -> None:
    volume = np.zeros((40, 20, 921), dtype=np.uint8)
    volume[15:25, 8:12, 300:622] = 200

    slices_a = _select_axial_slices(volume, n_show=12, bg_threshold=2)
    slices_b = _select_axial_slices(volume, n_show=12, bg_threshold=2)

    assert np.array_equal(slices_a, slices_b)


def test_save_cord_annotation_preview_writes_png(tmp_path: Path) -> None:
    volume = np.zeros((40, 20, 30), dtype=np.uint8)
    annotation = np.zeros((40, 20, 30), dtype=np.uint16)
    volume[15:25, 8:12, 10:20] = 200
    annotation[16:24, 7:13, 12:18] = 5
    out = save_cord_annotation_preview(volume, annotation, tmp_path / "qc.png")
    assert out.is_file()
    assert out.stat().st_size > 0


def test_save_cord_annotation_preview_safe_on_background_thread(tmp_path: Path) -> None:
    matplotlib.use("QtAgg", force=True)
    volume = np.zeros((40, 20, 30), dtype=np.uint8)
    annotation = np.zeros((40, 20, 30), dtype=np.uint16)
    volume[15:25, 8:12, 10:20] = 200
    annotation[16:24, 7:13, 12:18] = 5
    errors: list[BaseException] = []

    def _run() -> None:
        try:
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                save_cord_annotation_preview(volume, annotation, tmp_path / "thread.png")
            gui_warnings = [
                w
                for w in caught
                if issubclass(w.category, UserWarning)
                and "outside of the main thread" in str(w.message)
            ]
            if gui_warnings:
                errors.append(RuntimeError(str(gui_warnings[0].message)))
        except BaseException as exc:
            errors.append(exc)

    thread = threading.Thread(target=_run)
    thread.start()
    thread.join()
    assert not errors, errors
    assert (tmp_path / "thread.png").is_file()
