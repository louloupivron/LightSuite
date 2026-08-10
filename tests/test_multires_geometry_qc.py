"""Tests for multiresolution geometry QC helpers."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import tifffile
import yaml

from lightsuite.config.loader import (
    load_multires_config,
    save_mesospim_lateral_flip_to_multires_config,
)
from lightsuite.gui.inspect_geometry_multires import run_multires_inspect_geometry
from lightsuite.mesospim.meta import meta_path_for_tiff
from lightsuite.multires.geometry_qc import (
    build_reference_specs_with_geometry,
    compute_geometry_qc_slice,
    lateral_flip_from_config,
    lateral_flip_tuple,
    mesospim_geometry_yaml_snippet,
)


def _write_stack(path: Path, shape_zyx: tuple[int, int, int], *, seed: int = 0) -> None:
    rng = np.random.default_rng(seed)
    arr = rng.integers(100, 200, size=shape_zyx, dtype=np.uint16)
    arr[:, 4:12, 4:12] = 4000
    tifffile.imwrite(path, arr, imagej=True)


def _write_meta(path: Path, *, x_pos: float = 100.0, y_pos: float = 500.0) -> None:
    path.write_text(
        "\n".join(
            [
                "[Pixelsize in um] 5.0",
                f"[x_pos] {x_pos}",
                f"[y_pos] {y_pos}",
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


def _write_channels_config(tmp_path: Path, overview: Path, roi: Path) -> Path:
    config_path = tmp_path / "multires.yaml"
    config_path.write_text(
        yaml.safe_dump(
            {
                "sample": {"name": "sample", "save_path": str(tmp_path / "out")},
                "multires": {
                    "pair_label": "test",
                    "channels": {
                        "488": {
                            "overview": str(overview),
                            "roi": str(roi),
                        }
                    },
                    "mesospim_geometry": {
                        "overview": {"lateral_flip": [-1, 1]},
                        "roi": {"lateral_flip": [-1, 1]},
                    },
                    "registration": {"reference_channel": "488"},
                },
            }
        ),
        encoding="utf-8",
    )
    return config_path


def test_lateral_flip_helpers() -> None:
    assert lateral_flip_tuple(flip_x=True, flip_y=False) == (-1, 1)
    assert mesospim_geometry_yaml_snippet((-1, 1)).strip().startswith("mesospim_geometry:")


def test_lateral_flip_from_config(tmp_path: Path) -> None:
    overview = tmp_path / "overview.tif"
    roi = tmp_path / "roi.tif"
    _write_stack(overview, (5, 16, 16))
    _write_stack(roi, (5, 16, 16))
    _write_meta(meta_path_for_tiff(overview))
    _write_meta(meta_path_for_tiff(roi))
    cfg = load_multires_config(_write_channels_config(tmp_path, overview, roi))
    assert lateral_flip_from_config(cfg) == (-1, 1)


def test_compute_geometry_qc_slice_matching_geometry(tmp_path: Path) -> None:
    overview = tmp_path / "overview.tif"
    roi = tmp_path / "roi.tif"
    _write_stack(overview, (5, 16, 16), seed=1)
    _write_stack(roi, (5, 16, 16), seed=1)
    _write_meta(meta_path_for_tiff(overview))
    _write_meta(meta_path_for_tiff(roi))

    cfg = load_multires_config(_write_channels_config(tmp_path, overview, roi))
    overview_spec, roi_spec, manifest_dir = build_reference_specs_with_geometry(cfg, (-1, 1))
    good = compute_geometry_qc_slice(
        overview_spec,
        roi_spec,
        manifest_dir=manifest_dir,
        lateral_flip=(-1, 1),
    )
    bad = compute_geometry_qc_slice(
        overview_spec,
        roi_spec,
        manifest_dir=manifest_dir,
        lateral_flip=(1, -1),
    )
    assert good.has_overlap
    assert good.physical_ncc > 0.5
    assert good.physical_ncc >= bad.physical_ncc


def test_run_multires_inspect_geometry_headless(tmp_path: Path) -> None:
    overview = tmp_path / "overview.tif"
    roi = tmp_path / "roi.tif"
    _write_stack(overview, (5, 16, 16))
    _write_stack(roi, (5, 16, 16))
    _write_meta(meta_path_for_tiff(overview))
    _write_meta(meta_path_for_tiff(roi))
    cfg = load_multires_config(_write_channels_config(tmp_path, overview, roi))
    result = run_multires_inspect_geometry(cfg, headless=True)
    assert result.has_overlap
    assert result.lateral_flip == (-1, 1)


def test_save_mesospim_lateral_flip_updates_existing_block(tmp_path: Path) -> None:
    overview = tmp_path / "overview.tif"
    roi = tmp_path / "roi.tif"
    _write_stack(overview, (5, 16, 16))
    _write_stack(roi, (5, 16, 16))
    _write_meta(meta_path_for_tiff(overview))
    _write_meta(meta_path_for_tiff(roi))
    config_path = _write_channels_config(tmp_path, overview, roi)

    saved = save_mesospim_lateral_flip_to_multires_config(config_path, (1, -1))
    text = saved.read_text(encoding="utf-8")
    assert "lateral_flip: [1, -1]" in text
    cfg = load_multires_config(saved)
    assert lateral_flip_from_config(cfg) == (1, -1)


def test_save_mesospim_lateral_flip_inserts_missing_block(tmp_path: Path) -> None:
    overview = tmp_path / "overview.tif"
    roi = tmp_path / "roi.tif"
    _write_stack(overview, (5, 16, 16))
    _write_stack(roi, (5, 16, 16))
    _write_meta(meta_path_for_tiff(overview))
    _write_meta(meta_path_for_tiff(roi))
    config_path = tmp_path / "multires.yaml"
    config_path.write_text(
        yaml.safe_dump(
            {
                "sample": {"name": "sample", "save_path": str(tmp_path / "out")},
                "multires": {
                    "pair_label": "test",
                    "channels": {
                        "488": {
                            "overview": str(overview),
                            "roi": str(roi),
                        }
                    },
                    "registration": {"reference_channel": "488"},
                },
            }
        ),
        encoding="utf-8",
    )

    saved = save_mesospim_lateral_flip_to_multires_config(config_path, (-1, 1))
    cfg = load_multires_config(saved)
    assert lateral_flip_from_config(cfg) == (-1, 1)


def test_headless_write_config(tmp_path: Path) -> None:
    overview = tmp_path / "overview.tif"
    roi = tmp_path / "roi.tif"
    _write_stack(overview, (5, 16, 16))
    _write_stack(roi, (5, 16, 16))
    _write_meta(meta_path_for_tiff(overview))
    _write_meta(meta_path_for_tiff(roi))
    config_path = _write_channels_config(tmp_path, overview, roi)
    cfg = load_multires_config(config_path)

    run_multires_inspect_geometry(
        cfg,
        config_path=config_path,
        headless=True,
        write_config=True,
    )
    cfg = load_multires_config(config_path)
    assert lateral_flip_from_config(cfg) == (-1, 1)
