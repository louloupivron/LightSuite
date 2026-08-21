"""Tests for napari stage attach infrastructure."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from lightsuite.cli.stage_registry import StageContext
from lightsuite.cli.stages import brain_stage_specs, multires_stage_specs, spinal_stage_specs
from lightsuite.config.loader import load_config, load_spinal_config
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
from tests.test_stage_registry import _write_brain_config, _write_spinal_config


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


def test_multires_stage_specs_include_inspect_geometry_for_mesospim(tmp_path) -> None:
    from lightsuite.multires.config_models import MultiresPipelineConfig

    overview = tmp_path / "overview"
    roi = tmp_path / "roi"
    overview.mkdir()
    roi.mkdir()
    raw = {
        "sample": {
            "name": "test",
            "save_path": str(tmp_path / "results"),
            "scratch": str(tmp_path / "scratch"),
        },
        "multires": {
            "vendor": {"suite": "mesospim"},
            "pair_label": "pair1",
            "channels": {
                "488": {
                    "overview": str(overview),
                    "roi": str(roi),
                }
            },
            "registration": {"reference_channel": "488"},
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


def test_multires_stage_specs_omit_match_points_for_metadata_mode(tmp_path) -> None:
    from lightsuite.multires.config_models import MultiresPipelineConfig

    overview = tmp_path / "overview"
    roi = tmp_path / "roi"
    overview.mkdir()
    roi.mkdir()
    raw = {
        "sample": {
            "name": "test",
            "save_path": str(tmp_path / "results"),
            "scratch": str(tmp_path / "scratch"),
        },
        "multires": {
            "vendor": {"suite": "smartspim"},
            "pair_label": "pair1",
            "geometry_mode": "metadata",
            "channels": {
                "488": {
                    "overview": str(overview),
                    "roi": str(roi),
                }
            },
            "registration": {"reference_channel": "488"},
        },
    }
    (tmp_path / "results").mkdir()
    (tmp_path / "scratch").mkdir()
    cfg = MultiresPipelineConfig.model_validate(raw)
    ids = [spec.id for spec in multires_stage_specs(cfg)]
    assert "match-points" not in ids
    assert "check-geometry" in ids


def test_multires_stage_specs_include_match_points_for_hybrid_mode(tmp_path) -> None:
    from lightsuite.multires.config_models import MultiresPipelineConfig

    overview = tmp_path / "overview"
    roi = tmp_path / "roi"
    overview.mkdir()
    roi.mkdir()
    raw = {
        "sample": {
            "name": "test",
            "save_path": str(tmp_path / "results"),
            "scratch": str(tmp_path / "scratch"),
        },
        "multires": {
            "vendor": {"suite": "smartspim"},
            "pair_label": "pair1",
            "geometry_mode": "hybrid",
            "channels": {
                "488": {
                    "overview": str(overview),
                    "roi": str(roi),
                }
            },
            "registration": {"reference_channel": "488"},
        },
    }
    (tmp_path / "results").mkdir()
    (tmp_path / "scratch").mkdir()
    cfg = MultiresPipelineConfig.model_validate(raw)
    ids = [spec.id for spec in multires_stage_specs(cfg)]
    assert "match-points" in ids
    assert ids.index("match-points") < ids.index("check-geometry")


def test_multires_stage_specs_omit_inspect_geometry_for_smartspim(tmp_path) -> None:
    from lightsuite.multires.config_models import MultiresPipelineConfig

    raw = {
        "sample": {
            "name": "test",
            "save_path": str(tmp_path / "results"),
            "scratch": str(tmp_path / "scratch"),
        },
        "multires": {
            "vendor": {"suite": "smartspim"},
            "pair_label": "pair1",
            "channels": {
                "488": {
                    "overview": str(tmp_path / "overview"),
                    "roi": str(tmp_path / "roi"),
                }
            },
            "registration": {"reference_channel": "488"},
        },
    }
    (tmp_path / "results").mkdir()
    (tmp_path / "scratch").mkdir()
    (tmp_path / "overview").mkdir()
    (tmp_path / "roi").mkdir()
    cfg = MultiresPipelineConfig.model_validate(raw)
    ids = [spec.id for spec in multires_stage_specs(cfg)]
    assert "inspect-geometry" not in ids
    assert "check-geometry" in ids


def test_multires_manual_stages_have_attach(tmp_path) -> None:
    from lightsuite.multires.config_models import MultiresPipelineConfig

    overview = tmp_path / "overview"
    roi = tmp_path / "roi"
    overview.mkdir()
    roi.mkdir()
    raw = {
        "sample": {
            "name": "test",
            "save_path": str(tmp_path / "results"),
            "scratch": str(tmp_path / "scratch"),
        },
        "multires": {
            "vendor": {"suite": "mesospim"},
            "pair_label": "pair1",
            "channels": {
                "488": {
                    "overview": str(overview),
                    "roi": str(roi),
                }
            },
            "registration": {"reference_channel": "488"},
        },
    }
    (tmp_path / "results").mkdir()
    cfg = MultiresPipelineConfig.model_validate(raw)
    manual_ids = {spec.id for spec in multires_stage_specs(cfg) if spec.manual}
    missing = manual_ids - {stage_id for wf, stage_id in STAGE_ATTACH if wf == "multires"}
    assert not missing, f"Missing attach functions for: {missing}"


def test_stage_attach_covers_manual_spinal_stages(tmp_path: Path) -> None:
    cfg_path = tmp_path / "spinal.yaml"
    save_path = tmp_path / "results"
    _write_spinal_config(cfg_path, save_path)
    cfg = load_spinal_config(cfg_path)
    manual_ids = {spec.id for spec in spinal_stage_specs(cfg) if spec.manual}
    missing = manual_ids - {stage_id for wf, stage_id in STAGE_ATTACH if wf == "spinal"}
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


def test_require_napari_raises_on_incomplete_vispy(monkeypatch) -> None:
    from lightsuite.gui import stage_controller

    monkeypatch.setattr(stage_controller, "_validate_pint_install", lambda: None)
    monkeypatch.setattr(stage_controller, "_validate_napari_install", lambda: None)
    monkeypatch.setattr(
        stage_controller,
        "_validate_vispy_install",
        lambda: (_ for _ in ()).throw(
            RuntimeError("missing spatial-filters.npy"),
        ),
    )
    with patch.dict("sys.modules", {"napari": MagicMock()}):
        with pytest.raises(RuntimeError, match="spatial-filters"):
            stage_controller.require_napari()


def test_validate_vispy_install_accepts_complete_package() -> None:
    from lightsuite.gui.stage_controller import _validate_vispy_install

    try:
        import vispy
    except ImportError:
        pytest.skip("vispy not installed")
    spatial_filters = (
        Path(vispy.__file__).resolve().parent / "io" / "_data" / "spatial-filters.npy"
    )
    if not spatial_filters.is_file():
        pytest.skip("vispy installed without spatial-filters data")
    _validate_vispy_install()


def test_validate_vispy_install_raises_when_data_missing(tmp_path, monkeypatch) -> None:
    from lightsuite.gui.stage_controller import _validate_vispy_install

    vispy_root = tmp_path / "vispy_pkg"
    vispy_root.mkdir()
    (vispy_root / "io" / "_data").mkdir(parents=True)
    fake_vispy = MagicMock()
    fake_vispy.__file__ = str(vispy_root / "__init__.py")
    monkeypatch.setitem(__import__("sys").modules, "vispy", fake_vispy)
    with pytest.raises(RuntimeError, match="reinstall-package vispy"):
        _validate_vispy_install()


def test_validate_napari_install_accepts_complete_package() -> None:
    from lightsuite.gui.stage_controller import _validate_napari_install

    try:
        import napari
    except ImportError:
        pytest.skip("napari not installed")
    logo = (
        Path(napari.__file__).resolve().parent
        / "resources"
        / "logos"
        / "gradient-plain-dark.svg"
    )
    if not logo.is_file():
        pytest.skip("napari installed without bundled logo resources")
    _validate_napari_install()


def test_validate_napari_install_raises_when_resources_missing(tmp_path, monkeypatch) -> None:
    from lightsuite.gui.stage_controller import _validate_napari_install

    napari_root = tmp_path / "napari_pkg"
    napari_root.mkdir()
    (napari_root / "resources" / "logos").mkdir(parents=True)
    fake_napari = MagicMock()
    fake_napari.__file__ = str(napari_root / "__init__.py")
    monkeypatch.setitem(__import__("sys").modules, "napari", fake_napari)
    with pytest.raises(RuntimeError, match="uv sync --extra gui"):
        _validate_napari_install()


def test_validate_pint_install_raises_when_package_missing(monkeypatch) -> None:
    from lightsuite.gui.stage_controller import _validate_pint_install

    monkeypatch.delitem(__import__("sys").modules, "pint", raising=False)
    with patch.dict("sys.modules", {"pint": None}):
        with pytest.raises(RuntimeError, match="uv sync --extra gui"):
            _validate_pint_install()


def test_validate_pint_install_raises_when_definitions_missing(tmp_path, monkeypatch) -> None:
    from lightsuite.gui.stage_controller import _validate_pint_install

    pint_root = tmp_path / "pint_pkg"
    pint_root.mkdir()
    fake_pint = MagicMock()
    fake_pint.__file__ = str(pint_root / "__init__.py")
    monkeypatch.setitem(__import__("sys").modules, "pint", fake_pint)
    with pytest.raises(RuntimeError, match="default_en.txt"):
        _validate_pint_install()
