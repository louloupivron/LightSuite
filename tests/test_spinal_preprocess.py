"""Integration tests for spinal cord preprocess MVP."""

from __future__ import annotations

from pathlib import Path

import pytest
import tifffile
import yaml

from lightsuite.config.loader import load_spinal_config
from lightsuite.gui.straighten_cord import run_spinal_straighten
from lightsuite.preprocess.cord import preprocess_spinal_cord_sample
from lightsuite.preprocess.cord_checkpoint import CordRegOptsCheckpoint, SpinalAlignmentCheckpoint
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


def test_spinal_preprocess_and_straighten_headless(tmp_path: Path) -> None:
    fixture_root = Path(__file__).resolve().parent / "fixtures" / "spinal_cord"
    _ensure_fixtures(fixture_root)
    config_path = _write_config(tmp_path, fixture_root)
    cfg = load_spinal_config(config_path)
    save_cord_orientation(cfg.sample.save_path, CAUDOROSTRAL)

    result = preprocess_spinal_cord_sample(cfg)
    assert (Path(result.checkpoint.lsfolder) / "regopts.json").is_file()
    assert (Path(result.checkpoint.lsfolder) / "sample_reference.json").is_file()
    assert result.checkpoint.tofliprc is True
    regvol = tifffile.imread(result.checkpoint.regvol_path)
    assert regvol.ndim == 3
    assert regvol.size > 0

    align_path = run_spinal_straighten(cfg, headless=True)
    assert align_path.is_file()
    align = SpinalAlignmentCheckpoint.load(align_path)
    assert len(align.fit_x) == regvol.shape[2]

    checkpoint = CordRegOptsCheckpoint.load(Path(result.checkpoint.lsfolder) / "regopts.json")
    assert checkpoint.regvol_path
    assert checkpoint.regvolpaths
    assert "1" in checkpoint.regvolpaths
    assert Path(checkpoint.regvolpaths["1"]).is_file()


def test_spinal_preprocess_requires_orientation(tmp_path: Path) -> None:
    fixture_root = Path(__file__).resolve().parent / "fixtures" / "spinal_cord"
    _ensure_fixtures(fixture_root)
    config_path = _write_config(tmp_path, fixture_root)
    cfg = load_spinal_config(config_path)
    with pytest.raises(FileNotFoundError, match="cord_orientation.txt"):
        preprocess_spinal_cord_sample(cfg)
