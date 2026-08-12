"""Tests for napari stage attach infrastructure."""

from __future__ import annotations

from unittest.mock import MagicMock, patch

import pytest

from lightsuite.cli.stage_registry import StageContext
from lightsuite.cli.stages import brain_stage_specs, multires_stage_specs
from lightsuite.config.loader import load_config
from lightsuite.gui.stage_attach import STAGE_ATTACH, attach_stage, get_stage_attach
from lightsuite.gui.stage_controller import (
    DockStageController,
    clear_viewer_layers_safely,
    close_stage_or_viewer,
    mount_dock_widgets,
    remove_dock_widget,
    run_attached_stage,
)
from lightsuite.reporter import CallbackReporter, NullReporter
from tests.test_stage_registry import _write_brain_config


def test_stage_attach_covers_manual_brain_stages(tmp_path) -> None:
    cfg_path = tmp_path / "brain.yaml"
    _write_brain_config(cfg_path, tmp_path)
    cfg = load_config(cfg_path)
    manual_ids = {spec.id for spec in brain_stage_specs(cfg) if spec.manual}
    missing = manual_ids - {stage_id for wf, stage_id in STAGE_ATTACH if wf == "brain"}
    assert not missing, f"Missing attach functions for: {missing}"


def test_stage_attach_covers_multires_inspect_geometry() -> None:
    assert get_stage_attach("multires", "inspect-geometry") is not None


def test_stage_attach_covers_multires_inspect_registration() -> None:
    assert get_stage_attach("multires", "inspect-registration") is not None


def test_multires_stage_specs_include_inspect_geometry(tmp_path) -> None:
    from lightsuite.multires.config_models import MultiresPipelineConfig

    manifest = tmp_path / "pair.json"
    manifest.write_text("{}", encoding="utf-8")
    raw = {
        "sample": {
            "name": "test",
            "save_path": str(tmp_path / "results"),
            "scratch": str(tmp_path / "scratch"),
        },
        "multires": {
            "pair_label": "pair1",
            "pair_manifest": str(manifest),
        },
    }
    (tmp_path / "results").mkdir()
    (tmp_path / "scratch").mkdir()
    cfg = MultiresPipelineConfig.model_validate(raw)
    ids = [spec.id for spec in multires_stage_specs(cfg)]
    assert "inspect-geometry" in ids
    assert ids.index("inspect-geometry") < ids.index("check-geometry")
    assert "inspect-registration" in ids
    assert ids.index("register") < ids.index("inspect-registration")


def test_multires_manual_stages_have_attach(tmp_path) -> None:
    from lightsuite.multires.config_models import MultiresPipelineConfig

    manifest = tmp_path / "pair.json"
    manifest.write_text("{}", encoding="utf-8")
    raw = {
        "sample": {
            "name": "test",
            "save_path": str(tmp_path / "results"),
            "scratch": str(tmp_path / "scratch"),
        },
        "multires": {
            "pair_label": "pair1",
            "pair_manifest": str(manifest),
        },
    }
    (tmp_path / "results").mkdir()
    cfg = MultiresPipelineConfig.model_validate(raw)
    manual_ids = {spec.id for spec in multires_stage_specs(cfg) if spec.manual}
    missing = manual_ids - {stage_id for wf, stage_id in STAGE_ATTACH if wf == "multires"}
    assert not missing, f"Missing attach functions for: {missing}"


def test_get_stage_attach_unknown_returns_none() -> None:
    assert get_stage_attach("brain", "preprocess") is None
    assert get_stage_attach("brain", "not-a-stage") is None


def test_attach_stage_unknown_raises() -> None:
    viewer = MagicMock()
    ctx = StageContext(config_path="cfg.yaml")
    with pytest.raises(ValueError, match="No attach function"):
        attach_stage("brain", "preprocess", viewer, MagicMock(), ctx)


def test_clear_viewer_layers_safely_hides_text_before_clear() -> None:
    viewer = MagicMock()
    layer = MagicMock()
    layer.text = MagicMock()
    layers = MagicMock()
    layers.__iter__.return_value = iter([layer])
    viewer.layers = layers
    clear_viewer_layers_safely(viewer)
    assert layer.text.visible is False
    assert layer.visible is False
    layers.clear.assert_called_once()


def test_dock_stage_controller_mount_calls_refresh() -> None:
    viewer = MagicMock()
    refresh = MagicMock()
    widget = MagicMock()
    controller = DockStageController(
        dock_widgets=[(widget, "Test")],
        _refresh_fn=refresh,
    )
    controller.mount(viewer)
    viewer.window.add_dock_widget.assert_called_once_with(widget, area="right", name="Test")
    refresh.assert_called_once()


def test_mount_dock_widgets_adds_all() -> None:
    viewer = MagicMock()
    w1, w2 = MagicMock(), MagicMock()
    handles = mount_dock_widgets(viewer, [(w1, "One"), (w2, "Two")])
    assert viewer.window.add_dock_widget.call_count == 2
    assert len(handles) == 2


def test_remove_dock_widget_falls_back_to_native() -> None:
    viewer = MagicMock()
    widget = MagicMock()
    native = MagicMock()
    widget.native = native
    viewer.window.remove_dock_widget.side_effect = [LookupError("missing"), None]
    remove_dock_widget(viewer, widget)
    assert viewer.window.remove_dock_widget.call_args_list == [
        ((widget,),),
        ((native,),),
    ]


@patch("lightsuite.gui.stage_controller.require_napari")
def test_run_attached_stage_creates_viewer_and_runs(mock_require_napari) -> None:
    mock_napari = MagicMock()
    mock_viewer = MagicMock()
    mock_napari.Viewer.return_value = mock_viewer
    mock_require_napari.return_value = mock_napari
    controller = DockStageController()
    attach_fn = MagicMock(return_value=controller)

    run_attached_stage("Test title", attach_fn)

    mock_napari.Viewer.assert_called_once_with(title="Test title")
    attach_fn.assert_called_once_with(mock_viewer)
    mock_napari.run.assert_called_once()


def test_callback_reporter_pipeline_events() -> None:
    on_complete = MagicMock()
    reporter = CallbackReporter(on_pipeline_complete=on_complete)
    reporter.pipeline_complete(5)
    on_complete.assert_called_once_with(5)


def test_null_reporter_used_in_tests() -> None:
    reporter = NullReporter()
    reporter.stage_start("Register", checkpoint_hint="transform_params.json")
    reporter.pipeline_complete(1)


def test_dock_stage_controller_teardown_removes_docks() -> None:
    viewer = MagicMock()
    widget = MagicMock()
    dock_handle = MagicMock()
    viewer.window.add_dock_widget.return_value = dock_handle
    controller = DockStageController(dock_widgets=[(widget, "Test")])
    controller.mount(viewer)
    controller.teardown(viewer)
    viewer.window.remove_dock_widget.assert_called_once_with(dock_handle)


def test_close_stage_or_viewer_closes_standalone_viewer() -> None:
    viewer = MagicMock()
    viewer.window._lightsuite_shell = None
    close_stage_or_viewer(viewer)
    viewer.close.assert_called_once()


def test_close_stage_or_viewer_finishes_shell_stage() -> None:
    viewer = MagicMock()
    shell = MagicMock()
    viewer.window._lightsuite_shell = shell
    close_stage_or_viewer(viewer)
    shell.finish_interactive_stage.assert_called_once_with(refresh=True)
    viewer.close.assert_not_called()


def test_shell_auto_stage_does_not_import_napari_qt() -> None:
    import inspect

    from lightsuite.gui.shell import LightsuiteShell

    source = inspect.getsource(LightsuiteShell._run_auto_stage)
    assert "napari.qt" not in source
