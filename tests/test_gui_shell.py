"""Tests for workflow detection and GUI shell helpers."""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml

from lightsuite.config.workflow import detect_workflow, load_project
from lightsuite.gui.shell import (
    _strip_rich_markup,
    filter_stage_statuses,
    interactive_loading_message,
)
from lightsuite.cli.stages import StageSpec, StageState, StageStatus


def test_detect_workflow_brain() -> None:
    raw = {"sample": {"name": "x"}, "atlas": {"provider": "allen"}}
    assert detect_workflow(raw) == "brain"


def test_detect_workflow_spinal() -> None:
    raw = {"sample": {"name": "x"}, "atlas": {"atlas_dir": "/path/to/atlas"}}
    assert detect_workflow(raw) == "spinal"


def test_detect_workflow_multires() -> None:
    raw = {"sample": {"name": "x"}, "multires": {"pair_label": "pair1", "pair_manifest": "/p.json"}}
    assert detect_workflow(raw) == "multires"


def test_detect_workflow_brain_with_multires_link() -> None:
    raw = {
        "sample": {"name": "x"},
        "atlas": {"provider": "allen"},
        "multires": {"config": "examples/config/multiresolution/gilda_tg14_multires.yaml"},
    }
    assert detect_workflow(raw) == "brain"


def test_load_project_brain(tmp_path: Path) -> None:
    cfg = tmp_path / "brain.yaml"
    data = {
        "sample": {
            "name": "mouse",
            "source": {
                "format": "tiff_stack",
                "path": str(tmp_path / "src"),
                "tiff_type": "channelperfile",
            },
            "scratch": str(tmp_path / "scratch"),
            "save_path": str(tmp_path / "results"),
            "voxel_um": [5.0, 5.0, 5.0],
        },
        "atlas": {"provider": "allen", "resolution_um": 10, "atlas_dir": str(tmp_path / "atlas")},
        "registration": {"resolution_um": 20, "channel_primary": 1},
        "detection": {"enabled": False},
        "compute": {"workers": 1},
        "export": {"registered_volume_format": "tiff"},
    }
    cfg.write_text(yaml.safe_dump(data), encoding="utf-8")
    (tmp_path / "src").mkdir()
    (tmp_path / "scratch").mkdir()
    (tmp_path / "results").mkdir()
    (tmp_path / "atlas").mkdir()

    workflow, loaded = load_project(cfg)
    assert workflow == "brain"
    assert loaded.sample.name == "mouse"


def test_load_project_missing_file(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError):
        load_project(tmp_path / "missing.yaml")


def test_strip_rich_markup() -> None:
    assert _strip_rich_markup("[bold green]Done[/bold green]") == "Done"


def test_interactive_loading_message() -> None:
    assert interactive_loading_message("Match control points (GUI)", "match-points") == (
        "Loading Match control points (GUI)… (reading volumes and control-point session)"
    )
    assert interactive_loading_message("Custom stage", "unknown") == "Loading Custom stage…"


def test_filter_stage_statuses_hides_optional() -> None:
    optional = StageStatus(
        StageSpec("align-slices", "Align slices", "slice_correspondence.json", optional=True),
        StageState.OPTIONAL,
        "not done",
    )
    required = StageStatus(
        StageSpec("register", "Register", "transform_params.json"),
        StageState.PENDING,
        "missing",
    )
    statuses = [optional, required]
    hidden = filter_stage_statuses(statuses, show_optional=False)
    assert [item.stage.id for item in hidden] == ["register"]
    assert [item.stage.id for item in filter_stage_statuses(statuses, show_optional=True)] == [
        "align-slices",
        "register",
    ]
