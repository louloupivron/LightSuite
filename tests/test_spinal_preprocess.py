"""Integration tests for spinal cord preprocess MVP."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import tifffile
import yaml

from lightsuite.config.loader import load_spinal_config
from lightsuite.gui.straighten_cord import run_spinal_straighten
from lightsuite.preprocess.cord import detect_tofliprc, preprocess_spinal_cord_sample
from lightsuite.preprocess.cord_checkpoint import CordRegOptsCheckpoint, SpinalAlignmentCheckpoint


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

    result = preprocess_spinal_cord_sample(cfg)
    assert (Path(result.checkpoint.lsfolder) / "regopts.json").is_file()
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


def test_detect_tofliprc_when_back_cross_section_is_larger() -> None:
    """Caudorostral samples have thicker cord at the high-Z end (tofliprc=True)."""
    regvol = np.zeros((20, 20, 100), dtype=np.uint16)
    regvol[5:15, 5:15, 30:] = 1000
    regvol[7:13, 7:13, :30] = 1000
    tv = np.zeros((10, 10, 200), dtype=np.float32)
    assert detect_tofliprc(regvol, tv) is True
    assert detect_tofliprc(np.flip(regvol, axis=2), tv) is False
