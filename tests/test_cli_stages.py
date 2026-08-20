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
    spinal_stage_statuses,
)
from lightsuite.config.loader import load_config, load_spinal_config
from lightsuite.preprocess.checkpoint import RegOptsCheckpoint
from lightsuite.preprocess.cord_checkpoint import CordRegOptsCheckpoint


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


def test_brain_stage_specs_include_align_slices_when_correspondence_disabled(
    tmp_path: Path,
) -> None:
    cfg_path = tmp_path / "brain.yaml"
    _write_brain_config(cfg_path, tmp_path)
    raw = yaml.safe_load(cfg_path.read_text(encoding="utf-8"))
    raw["registration"]["use_slice_correspondence_affine"] = False
    raw["registration"]["use_slice_correspondence_landmarks"] = False
    cfg_path.write_text(yaml.safe_dump(raw), encoding="utf-8")
    cfg = load_config(cfg_path)
    ids = [spec.id for spec in brain_stage_specs(cfg)]
    align = next(spec for spec in brain_stage_specs(cfg) if spec.id == "align-slices")
    assert "align-slices" in ids
    assert align.optional is True


def test_brain_match_points_is_optional(tmp_path: Path) -> None:
    cfg_path = tmp_path / "brain.yaml"
    _write_brain_config(cfg_path, tmp_path)
    cfg = load_config(cfg_path)
    match_points = next(spec for spec in brain_stage_specs(cfg) if spec.id == "match-points")
    assert match_points.optional is True
    assert match_points.manual is True


def test_brain_import_segmentation_listed_when_configured(tmp_path: Path) -> None:
    cfg_path = tmp_path / "brain.yaml"
    _write_brain_config(cfg_path, tmp_path)
    cfg = load_config(cfg_path)
    ids = [spec.id for spec in brain_stage_specs(cfg)]
    assert "import-segmentation" not in ids
    assert "convert-annotations" not in ids
    assert "import-annotations" not in ids

    raw = yaml.safe_load(cfg_path.read_text(encoding="utf-8"))
    raw["import"] = {
        "converter": {
            "suite": "smartspim",
            "source": str(tmp_path / "points.json"),
            "label": "cells",
        }
    }
    cfg_path.write_text(yaml.safe_dump(raw), encoding="utf-8")
    (tmp_path / "points.json").write_text("[]", encoding="utf-8")
    cfg = load_config(cfg_path)
    specs = brain_stage_specs(cfg)
    ids = [spec.id for spec in specs]
    assert "import-segmentation" in ids
    assert "convert-annotations" not in ids
    assert "import-annotations" not in ids
    assert ids.index("export") < ids.index("import-segmentation") < ids.index("view-registration")
    stage = next(spec for spec in specs if spec.id == "import-segmentation")
    assert stage.optional is True

    statuses = brain_stage_statuses(cfg)
    status = next(item for item in statuses if item.stage.id == "import-segmentation")
    assert status.state == StageState.OPTIONAL
    assert "Config" in status.detail or "register" in status.detail

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


def _write_spinal_config(path: Path, save_path: Path) -> None:
    atlas = save_path.parent / "atlas"
    atlas.mkdir(parents=True, exist_ok=True)
    (atlas / "Template.tif").write_bytes(b"x")
    (atlas / "Annotation.tif").write_bytes(b"x")
    (atlas / "Segments.csv").write_text("Segment,Start,End\nC1,0,1\n", encoding="utf-8")
    (atlas / "Atlas_Regions.csv").write_text("id,name,children_IDs\n1,gm,1\n", encoding="utf-8")
    sample = save_path.parent / "sample"
    sample.mkdir(parents=True, exist_ok=True)
    save_path.mkdir(parents=True, exist_ok=True)
    data = {
        "sample": {
            "name": "test_cord",
            "source": {"path": str(sample), "tiff_type": "channelperfile"},
            "scratch": str(save_path.parent / "scratch"),
            "save_path": str(save_path),
            "voxel_um": [20.0, 20.0, 20.0],
        },
        "atlas": {"atlas_dir": str(atlas)},
        "registration": {"longitudinal_direction": "caudorostral"},
    }
    path.write_text(yaml.safe_dump(data), encoding="utf-8")


def test_spinal_preprocess_status_done_with_cord_regopts(tmp_path: Path) -> None:
    cfg_path = tmp_path / "spinal.yaml"
    save_path = tmp_path / "results"
    _write_spinal_config(cfg_path, save_path)
    cfg = load_spinal_config(cfg_path)
    regvol = save_path / "regvol.tif"
    regvol.write_bytes(b"x")
    checkpoint = CordRegOptsCheckpoint(
        sample_name="test_cord",
        data_folder=str(tmp_path / "sample"),
        lsfolder=str(save_path),
        orisize=[10, 10, 10],
        nchans=1,
        sampleres_um=[20.0, 20.0, 20.0],
        registrationres_um=[20.0, 20.0, 20.0],
        reg_channel=0,
        sample_perm=[1, 2, 3],
        tofliprc=False,
        ikeeprange=[0, 9],
        xrange=[0, 9],
        yrange=[0, 9],
        regvol_path=str(regvol),
        tv_path=str(save_path / "tv.tif"),
        av_path=str(save_path / "av.tif"),
        smpts_path=str(save_path / "smpts.npy"),
        tvpts_path=str(save_path / "tvpts.npy"),
        atlas_res_um=[10.0, 10.0, 20.0],
        segments_path=str(save_path / "segments.csv"),
        regions_path=str(save_path / "regions.csv"),
        tiff_type="channelperfile",
    )
    checkpoint.save(save_path / "regopts.json")

    statuses = spinal_stage_statuses(cfg)
    preprocess = next(item for item in statuses if item.stage.id == "preprocess")
    assert preprocess.state == StageState.DONE
