"""Tests for mesoSPIM Napari inspect path resolution."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import tifffile
import yaml

from lightsuite.gui.inspect_mesospim import resolve_mesospim_inspect_paths
from lightsuite.mesospim.io import load_registered_canvas_zyx
from lightsuite.mesospim.meta import meta_path_for_tiff


def _write_stack(path: Path, shape_zyx: tuple[int, int, int]) -> None:
    arr = np.zeros(shape_zyx, dtype=np.uint16)
    tifffile.imwrite(path, arr, imagej=True)


def _write_meta(path: Path) -> None:
    meta_path_for_tiff(path).write_text(
        "\n".join(
            [
                "[Pixelsize in um] 5.0",
                "[x_pos] 0.0",
                "[y_pos] 0.0",
                "[z_start] 0.0",
                "[z_end] 8.0",
                "[z_stepsize] 2.0",
                "[z_planes] 5",
                "[x_pixels] 16",
                "[y_pixels] 16",
            ]
        ),
        encoding="utf-8",
    )


def test_resolve_mesospim_inspect_paths(tmp_path: Path) -> None:
    overview = tmp_path / "overview.tif"
    roi = tmp_path / "roi.tif"
    shape = (5, 16, 16)
    _write_stack(overview, shape)
    _write_stack(roi, shape)
    _write_meta(overview)
    _write_meta(roi)

    slug = "test_exp"
    reg_dir = tmp_path / "out" / "elastix_roi_to_overview" / slug
    reg_dir.mkdir(parents=True)
    registered = reg_dir / f"{slug}_roi_registered_to_overview_in_full_overview.tif"
    _write_stack(registered, shape)

    config_path = tmp_path / "mesospim.yaml"
    config_path.write_text(
        yaml.safe_dump(
            {
                "sample": {"name": "sample", "save_path": str(tmp_path / "out")},
                "mesospim": {
                    "overview": {"path": str(overview)},
                    "roi": {"path": str(roi)},
                    "registration": {
                        "experiment_name": slug,
                        "write_full_overview_canvas": True,
                    },
                },
            }
        ),
        encoding="utf-8",
    )

    from lightsuite.config.loader import load_mesospim_config

    cfg = load_mesospim_config(config_path)
    paths = resolve_mesospim_inspect_paths(cfg)
    assert paths.overview_path == overview.resolve()
    assert paths.registered_full_overview_path == registered.resolve()


def test_load_registered_canvas_zyx_compressed(tmp_path: Path) -> None:
    path = tmp_path / "registered.tif"
    arr = np.arange(2 * 4 * 4, dtype=np.float32).reshape(2, 4, 4)
    tifffile.imwrite(path, arr, imagej=True, compression="zlib")

    loaded = load_registered_canvas_zyx(path)
    assert loaded.shape == (2, 4, 4)
    assert loaded.dtype == np.float32
    assert loaded[0, 0, 0] == arr[0, 0, 0]
