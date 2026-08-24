"""Tests for spinal cord registration-grid cache reuse."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

import yaml

from lightsuite.config.loader import load_spinal_config
from lightsuite.io.cord_registration_cache import (
    REGISTER_CACHE_MANIFEST,
    load_or_cache_cord_registration,
    register_cache_valid,
)
from lightsuite.registration.cord_orientation import CAUDOROSTRAL, save_cord_orientation


def _ensure_fixtures(root: Path) -> None:
    if not (root / "atlas" / "Template.tif").is_file():
        from tests.fixtures.spinal_cord.build_fixtures import build_fixtures

        build_fixtures(root)


def _write_config(tmp_path: Path, fixture_root: Path) -> Path:
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    save = tmp_path / "results"
    save.mkdir()
    config_data = {
        "sample": {
            "name": "fixture_cord",
            "source": {
                "path": str(fixture_root / "sample"),
                "tiff_type": "channelperfile",
            },
            "scratch": str(scratch),
            "save_path": str(save),
            "voxel_um": [20.0, 20.0, 20.0],
        },
        "atlas": {"atlas_dir": str(fixture_root / "atlas")},
        "registration": {"channel_primary": 1, "resolution_um": 20},
    }
    config_path = tmp_path / "spinal.yaml"
    config_path.write_text(yaml.dump(config_data), encoding="utf-8")
    return config_path


def test_register_cache_reused_between_calls(tmp_path: Path) -> None:
    fixture_root = Path(__file__).resolve().parent / "fixtures" / "spinal_cord"
    _ensure_fixtures(fixture_root)
    config_path = _write_config(tmp_path, fixture_root)
    cfg = load_spinal_config(config_path)

    first = load_or_cache_cord_registration(cfg)
    assert first.from_cache is False
    cache_dir = Path(cfg.sample.save_path) / "cache"
    assert (cache_dir / REGISTER_CACHE_MANIFEST).is_file()
    assert register_cache_valid(cache_dir, cfg)

    with patch(
        "lightsuite.io.cord_registration_cache.read_spinal_cord_sample",
        side_effect=AssertionError("should not reload raw TIFFs"),
    ), patch(
        "lightsuite.io.cord_registration_cache.probe_cord_source",
        side_effect=AssertionError("should not re-probe raw TIFFs on cache hit"),
    ):
        second = load_or_cache_cord_registration(cfg)
    assert second.from_cache is True
    assert second.volume.shape == first.volume.shape


def test_register_cache_valid_without_probe(tmp_path: Path) -> None:
    fixture_root = Path(__file__).resolve().parent / "fixtures" / "spinal_cord"
    _ensure_fixtures(fixture_root)
    config_path = _write_config(tmp_path, fixture_root)
    cfg = load_spinal_config(config_path)
    load_or_cache_cord_registration(cfg)
    cache_dir = Path(cfg.sample.save_path) / "cache"

    with patch(
        "lightsuite.io.cord_registration_cache.probe_cord_source",
        side_effect=AssertionError("register_cache_valid should not probe raw TIFFs"),
    ):
        assert register_cache_valid(cache_dir, cfg) is True


def test_orphan_registration_tiff_without_manifest(tmp_path: Path) -> None:
    fixture_root = Path(__file__).resolve().parent / "fixtures" / "spinal_cord"
    _ensure_fixtures(fixture_root)
    config_path = _write_config(tmp_path, fixture_root)
    cfg = load_spinal_config(config_path)

    first = load_or_cache_cord_registration(cfg)
    assert first.from_cache is False
    cache_dir = Path(cfg.sample.save_path) / "cache"
    manifest = cache_dir / REGISTER_CACHE_MANIFEST
    assert manifest.is_file()
    manifest.unlink()

    with patch(
        "lightsuite.io.cord_registration_cache.read_spinal_cord_sample",
        side_effect=AssertionError("orphan TIFF should be reused without raw reload"),
    ):
        second = load_or_cache_cord_registration(cfg)

    assert second.from_cache is True
    assert second.volume.shape == first.volume.shape
    assert manifest.is_file()
    assert register_cache_valid(cache_dir, cfg)


def test_preprocess_reuses_cache_after_orientation_load(tmp_path: Path) -> None:
    fixture_root = Path(__file__).resolve().parent / "fixtures" / "spinal_cord"
    _ensure_fixtures(fixture_root)
    config_path = _write_config(tmp_path, fixture_root)
    cfg = load_spinal_config(config_path)
    save_cord_orientation(cfg.sample.save_path, CAUDOROSTRAL)

    load_or_cache_cord_registration(cfg)

    from lightsuite.preprocess.cord import preprocess_spinal_cord_sample

    with patch(
        "lightsuite.io.cord_registration_cache.read_spinal_cord_sample",
        side_effect=AssertionError("preprocess should reuse cached TIFFs"),
    ), patch(
        "lightsuite.io.cord_registration_cache.probe_cord_source",
        side_effect=AssertionError("preprocess should not re-probe raw TIFFs"),
    ):
        result = preprocess_spinal_cord_sample(cfg)

    assert result.checkpoint.regvolpaths
    assert Path(result.checkpoint.regvolpaths["1"]).is_file()


def test_orientation_preview_reused_without_volume_reload(tmp_path: Path) -> None:
    fixture_root = Path(__file__).resolve().parent / "fixtures" / "spinal_cord"
    _ensure_fixtures(fixture_root)
    config_path = _write_config(tmp_path, fixture_root)
    cfg = load_spinal_config(config_path)

    from lightsuite.gui.orientation_cord import load_cord_orientation_check_data

    first = load_cord_orientation_check_data(cfg)
    assert first.sample_longitudinal.ndim == 2
    assert first.atlas_longitudinal.ndim == 2

    with patch(
        "lightsuite.gui.orientation_cord.load_or_cache_cord_registration",
        side_effect=AssertionError("orientation preview should not reload volumes"),
    ):
        second = load_cord_orientation_check_data(cfg)
    assert second.sample_longitudinal.shape == first.sample_longitudinal.shape
    assert second.atlas_longitudinal.shape == first.atlas_longitudinal.shape
