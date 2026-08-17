"""Shared helpers for napari stages attachable to a unified LightSuite GUI."""

from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass, field
from typing import Any, Protocol, runtime_checkable


def require_napari() -> Any:
    """Import napari or raise with install instructions."""
    try:
        import napari
    except ImportError as exc:
        msg = "Napari GUI requires: uv sync --extra gui"
        raise RuntimeError(msg) from exc
    _validate_pint_install()
    _validate_vispy_install()
    _validate_napari_install()
    return napari


def _validate_napari_install() -> None:
    """Catch incomplete napari / PyOpenGL installs before Qt tries to create a GL context."""
    from pathlib import Path

    import napari

    logo = Path(napari.__file__).resolve().parent / "resources" / "logos" / "gradient-plain-dark.svg"
    if not logo.is_file():
        msg = (
            "The napari package in this environment is incomplete "
            "(missing resources/logos/gradient-plain-dark.svg). "
            "Repair with: uv sync --extra gui"
        )
        raise RuntimeError(msg)

    try:
        from OpenGL.arrays import numpymodule  # noqa: F401
    except ImportError as exc:
        msg = (
            "PyOpenGL is missing or incomplete (required for Napari). "
            "Repair with: uv sync --extra gui"
        )
        raise RuntimeError(msg) from exc


def _validate_vispy_install() -> None:
    """Napari image layers need vispy spatial-filter data shipped with the wheel."""
    try:
        import vispy
    except ImportError:
        return
    from pathlib import Path

    spatial_filters = (
        Path(vispy.__file__).resolve().parent / "io" / "_data" / "spatial-filters.npy"
    )
    if spatial_filters.is_file():
        return
    msg = (
        "The vispy package in this environment is incomplete "
        "(missing io/_data/spatial-filters.npy). "
        "Repair with: uv sync --extra gui --reinstall-package vispy "
        "(or run `uv sync --extra gui` if other GUI packages are missing)"
    )
    raise RuntimeError(msg)


def _validate_pint_install() -> None:
    """Napari layer scales need pint unit definitions shipped with the wheel."""
    try:
        import pint
    except ImportError as exc:
        msg = (
            "Pint is required for Napari physical scales but is not installed. "
            "Repair with: uv sync --extra gui"
        )
        raise RuntimeError(msg) from exc
    from pathlib import Path

    definitions = Path(pint.__file__).resolve().parent / "default_en.txt"
    if definitions.is_file():
        return
    msg = (
        "The pint package in this environment is incomplete (missing default_en.txt). "
        "Repair with: uv sync --extra gui --reinstall-package pint"
    )
    raise RuntimeError(msg)


def validate_gui_dependencies() -> None:
    """Re-check GUI runtime packages before opening interactive stages."""
    _validate_pint_install()
    _validate_vispy_install()


def require_magicgui() -> Any:
    """Import magicgui or raise with install instructions."""
    try:
        from magicgui import magicgui
    except ImportError as exc:
        msg = "Napari GUI requires: uv sync --extra gui"
        raise RuntimeError(msg) from exc
    return magicgui


def mount_dock_widgets(
    viewer: Any,
    widgets: Sequence[tuple[Any, str]],
    *,
    area: str = "right",
) -> list[Any]:
    """Add multiple dock widgets to a napari viewer.

    Returns the napari dock handles so they can be removed reliably (magicgui
    ``Label`` widgets are not always discoverable by inner-widget identity).
    """
    handles: list[Any] = []
    for widget, name in widgets:
        handles.append(viewer.window.add_dock_widget(widget, area=area, name=name))
    return handles


def remove_dock_widget(viewer: Any, widget: Any, *, dock_handle: Any = None) -> None:
    """Remove a docked widget, tolerating magicgui wrapper/native mismatches."""
    window = viewer.window
    if dock_handle is not None:
        try:
            window.remove_dock_widget(dock_handle)
            return
        except LookupError:
            pass
        except (AttributeError, RuntimeError, TypeError, ValueError):
            return

    for candidate in (widget, getattr(widget, "native", None)):
        if candidate is None:
            continue
        try:
            window.remove_dock_widget(candidate)
            return
        except LookupError:
            continue
        except (AttributeError, RuntimeError, TypeError, ValueError):
            return


def clear_viewer_layers_safely(viewer: Any) -> None:
    """Remove napari layers without tripping vispy Text draw races on teardown."""
    layers = list(getattr(viewer, "layers", ()))
    for layer in layers:
        text = getattr(layer, "text", None)
        if text is not None:
            try:
                text.visible = False
            except (AttributeError, RuntimeError, TypeError, ValueError):
                pass
        try:
            layer.visible = False
        except (AttributeError, RuntimeError, TypeError, ValueError):
            pass
    try:
        viewer.layers.clear()
    except (AttributeError, RuntimeError, TypeError, ValueError):
        for layer in layers:
            try:
                viewer.layers.remove(layer)
            except (LookupError, AttributeError, RuntimeError, TypeError, ValueError):
                pass


def defer_clear_viewer_layers(viewer: Any) -> None:
    """Defer layer removal to the next Qt event-loop tick."""
    try:
        from qtpy.QtCore import QTimer
    except ImportError:
        clear_viewer_layers_safely(viewer)
        return
    QTimer.singleShot(0, lambda: clear_viewer_layers_safely(viewer))


def close_stage_or_viewer(viewer: Any, *, refresh: bool = True) -> None:
    """Close a standalone stage viewer, or detach the stage inside the unified shell."""
    window = getattr(viewer, "window", None)
    shell = getattr(window, "_lightsuite_shell", None) if window is not None else None
    if shell is not None:
        shell.finish_interactive_stage(refresh=refresh)
        return
    viewer.close()


@runtime_checkable
class StageController(Protocol):
    """Interactive stage logic bound to an existing napari viewer."""

    def mount(self, viewer: Any) -> None:
        """Attach dock widgets and finalize the stage UI."""

    def refresh(self) -> None:
        """Refresh layers and status from the current stage state."""

    def teardown(self, viewer: Any) -> None:
        """Remove UI artifacts before another stage is attached."""


@dataclass
class DockStageController:
    """Default controller that mounts a list of dock widgets."""

    dock_widgets: list[tuple[Any, str]] = field(default_factory=list)
    _dock_handles: list[Any] = field(default_factory=list, repr=False)
    _refresh_fn: Callable[[], None] | None = None
    _teardown_fn: Callable[[], None] | None = None
    open_log_message: str | None = None
    result: Any = None

    def mount(self, viewer: Any) -> None:
        self._dock_handles = mount_dock_widgets(viewer, self.dock_widgets)
        self.refresh()

    def refresh(self) -> None:
        if self._refresh_fn is not None:
            self._refresh_fn()

    def teardown(self, viewer: Any) -> None:
        """Remove dock widgets and run optional cleanup before switching stages."""
        if self._dock_handles:
            pairs = list(zip(self.dock_widgets, self._dock_handles, strict=False))
            for (widget, _name), handle in reversed(pairs):
                remove_dock_widget(viewer, widget, dock_handle=handle)
        else:
            for widget, _name in reversed(self.dock_widgets):
                remove_dock_widget(viewer, widget)
        self._dock_handles.clear()
        if self._teardown_fn is not None:
            self._teardown_fn()


def run_attached_stage(
    title: str,
    attach_fn: Callable[[Any], StageController],
    *,
    before_run: Callable[[], None] | None = None,
) -> Any:
    """Create a viewer, attach a stage, and block until the viewer closes."""
    napari = require_napari()
    viewer = napari.Viewer(title=title)
    controller = attach_fn(viewer)
    controller.mount(viewer)
    if before_run is not None:
        before_run()
    napari.run()
    return getattr(controller, "result", None)
