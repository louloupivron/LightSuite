"""Tests for pipeline stage graphs and orchestration helpers."""

from __future__ import annotations

import json
from pathlib import Path

import pytest
import yaml

from lightsuite.cli.stages import (
    StageState,
    brain_stage_specs,
    brain_stage_statuses,
    slice_stages,
)
from lightsuite.config.loader import load_config
from lightsuite.preprocess.checkpoint import RegOptsCheckpoint


def _write_brain_config(path: Path, save_path: Path) -> None:
    data = {
        "sample": {
            "name": "test_mouse",
            "source": {
                "format": "tiff_stack",
                "path": str(save_path / "source"),
                "tiff_type": "channelperfile",
            },
            "scratch": str(save_path / "scratch"),
            "save_path": str(save_path / "results"),
            "voxel_um": [5.0, 5.0, 5.0],
        },
        "atlas": {
            "provider": "allen",
            "resolution_um": 10,
            "atlas_dir": str(save_path / "atlas"),
        },
        "registration": {"resolution_um": 20, "channel_primary": 1},
        "detection": {"enabled": False},
        "compute": {"workers": 1},
        "export": {"registered_volume_format": "tiff"},
    }
    path.write_text(yaml.safe_dump(data), encoding="utf-8")
    (save_path / "source").mkdir(parents=True, exist_ok=True)
    (save_path / "scratch").mkdir(parents=True, exist_ok=True)
    (save_path / "results").mkdir(parents=True, exist_ok=True)
    (save_path / "atlas").mkdir(parents=True, exist_ok=True)


def test_brain_stage_specs_include_slice_stages_by_default(tmp_path: Path) -> None:
    cfg_path = tmp_path / "brain.yaml"
    _write_brain_config(cfg_path, tmp_path)
    cfg = load_config(cfg_path)
    ids = [spec.id for spec in brain_stage_specs(cfg)]
    assert "align-slices" in ids
    assert "refine-auto-points" not in ids


def test_brain_preprocess_status_done_with_regopts(tmp_path: Path) -> None:
    cfg_path = tmp_path / "brain.yaml"
    _write_brain_config(cfg_path, tmp_path)
    cfg = load_config(cfg_path)
    regvol = tmp_path / "results" / "regvol.tif"
    regvol.write_bytes(b"x")
    checkpoint = RegOptsCheckpoint(
        sample_name="test_mouse",
        ny=10,
        nx=10,
        nz=10,
        nchans=1,
        voxel_um=[5.0, 5.0, 5.0],
        registres_um=20.0,
        regvolpath=str(regvol),
        regvolpath_secondary=None,
        regvolpaths={"1": str(regvol)},
        tiff_type="channelperfile",
        channel_primary=1,
        channel_secondary=None,
    )
    checkpoint.save(tmp_path / "results" / "regopts.json")

    statuses = brain_stage_statuses(cfg)
    preprocess = next(item for item in statuses if item.stage.id == "preprocess")
    assert preprocess.state == StageState.DONE


def test_slice_stages_through(tmp_path: Path) -> None:
    cfg_path = tmp_path / "brain.yaml"
    _write_brain_config(cfg_path, tmp_path)
    cfg = load_config(cfg_path)
    specs = brain_stage_specs(cfg)
    selected = slice_stages(specs, from_stage=None, through_stage="init-registration")
    assert selected[-1].id == "init-registration"
    assert all(spec.id != "match-points" for spec in selected)


def test_load_config_raises_actionable_error(tmp_path: Path) -> None:
    from lightsuite.exceptions import LightsuiteConfigError

    bad = tmp_path / "bad.yaml"
    bad.write_text("sample:\n  name: test\n", encoding="utf-8")
    with pytest.raises(LightsuiteConfigError) as exc:
        load_config(bad)
    assert "config explain" in str(exc.value)
