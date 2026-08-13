"""Tests for stage registry and reporter abstractions."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock

import pytest
import yaml

from lightsuite.cli.stage_registry import (
    StageContext,
    StageKind,
    WORKFLOWS,
    get_workflow,
    run_stage,
    runner_stage_ids,
    stage_kind,
    validate_registry,
)
from lightsuite.cli.stages import brain_stage_specs, spinal_stage_specs
from lightsuite.config.loader import load_config, load_spinal_config
from lightsuite.reporter import CallbackReporter, ConsoleReporter, NullReporter


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


def test_validate_registry_passes() -> None:
    validate_registry()


def test_workflows_cover_brain_spinal_multires() -> None:
    assert set(WORKFLOWS) == {"brain", "spinal", "multires"}


def test_get_workflow_unknown_raises() -> None:
    with pytest.raises(KeyError, match="Unknown workflow"):
        get_workflow("widefield")


def test_brain_runner_ids_cover_stage_specs(tmp_path: Path) -> None:
    cfg_path = tmp_path / "brain.yaml"
    _write_brain_config(cfg_path, tmp_path)
    cfg = load_config(cfg_path)
    spec_ids = {spec.id for spec in brain_stage_specs(cfg)}
    runner_ids = runner_stage_ids("brain")
    assert spec_ids <= runner_ids


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


def test_spinal_runner_ids_cover_stage_specs(tmp_path: Path) -> None:
    cfg_path = tmp_path / "spinal.yaml"
    save_path = tmp_path / "results"
    _write_spinal_config(cfg_path, save_path)
    cfg = load_spinal_config(cfg_path)
    spec_ids = {spec.id for spec in spinal_stage_specs(cfg)}
    runner_ids = runner_stage_ids("spinal")
    assert spec_ids <= runner_ids
    assert "view-registration" in spec_ids


def test_stage_kind_reflects_manual_flag(tmp_path: Path) -> None:
    cfg_path = tmp_path / "brain.yaml"
    _write_brain_config(cfg_path, tmp_path)
    cfg = load_config(cfg_path)
    specs = brain_stage_specs(cfg)
    match_points = next(spec for spec in specs if spec.id == "match-points")
    preprocess = next(spec for spec in specs if spec.id == "preprocess")
    assert stage_kind(match_points) == StageKind.INTERACTIVE
    assert stage_kind(preprocess) == StageKind.AUTO


def test_run_stage_unknown_raises(tmp_path: Path) -> None:
    cfg_path = tmp_path / "brain.yaml"
    _write_brain_config(cfg_path, tmp_path)
    cfg = load_config(cfg_path)
    ctx = StageContext(config_path=cfg_path)
    with pytest.raises(ValueError, match="Unsupported brain stage"):
        run_stage("brain", "not-a-stage", cfg, ctx)


def test_null_reporter_is_silent() -> None:
    reporter = NullReporter()
    reporter.message("hello")
    reporter.stage_start("Preprocess", manual=True, checkpoint_hint="regopts.json")
    reporter.stage_skip("Preprocess", "done", "ok")
    reporter.stage_failed("preprocess", RuntimeError("boom"))
    reporter.rerun_hint("brain", "preprocess", Path("cfg.yaml"))
    reporter.pipeline_complete(3)
    reporter.pipeline_no_stages()


def test_callback_reporter_forwards_events() -> None:
    on_message = MagicMock()
    on_stage_start = MagicMock()
    on_pipeline_complete = MagicMock()
    reporter = CallbackReporter(
        on_message=on_message,
        on_stage_start=on_stage_start,
        on_pipeline_complete=on_pipeline_complete,
    )
    reporter.message("status")
    reporter.stage_start("Register", manual=False, checkpoint_hint="transform_params.json")
    reporter.pipeline_complete(2)
    on_message.assert_called_once_with("status")
    on_stage_start.assert_called_once_with(
        "Register",
        manual=False,
        checkpoint_hint="transform_params.json",
    )
    on_pipeline_complete.assert_called_once_with(2)


def test_console_reporter_exposes_console() -> None:
    reporter = ConsoleReporter()
    assert reporter.console is not None
