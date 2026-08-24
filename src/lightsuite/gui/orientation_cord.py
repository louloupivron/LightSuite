"""Napari GUI to manually set spinal cord rostrocaudal orientation."""

from __future__ import annotations

import json
import threading
from collections.abc import Callable
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import tifffile
from rich.console import Console
from skimage.transform import resize

from lightsuite.atlas.fiederling import load_fiederling_atlas_volumes, resize_fiederling_atlas
from lightsuite.config.models import SpinalCordPipelineConfig
from lightsuite.gui.stage_controller import (
    DockStageController,
    close_stage_or_viewer,
    require_magicgui,
    run_attached_stage,
)
from lightsuite.reporter import (
    CallbackReporter,
    capture_pipeline_output,
    emit_pipeline_message,
    stage_cancellation,
)
from lightsuite.exceptions import StageCancelledError
from lightsuite.io.cord_registration_cache import (
    load_or_cache_cord_registration,
    manifest_register_fingerprint,
)
from lightsuite.preprocess.cord_checkpoint import CordRegOptsCheckpoint
from lightsuite.registration.cord_orientation import (
    CAUDOROSTRAL,
    ROSTROCAUDAL,
    cord_orientation_path,
    load_cord_orientation,
    resolve_cord_orientation,
    save_cord_orientation,
)
from lightsuite.registration.cord_paths import cord_cache_dir, cord_save_path

console = Console()

PANEL_GAP_X = 24
_VALID_DIRECTIONS = (ROSTROCAUDAL, CAUDOROSTRAL)
ORIENTATION_PREVIEW_FILENAME = "orientation_preview.npz"


@dataclass
class CordOrientationData:
    sample_longitudinal: np.ndarray  # (Z, Y) max projection, rostral (low Z) at top
    atlas_longitudinal: np.ndarray  # (Z, Y) atlas reference, rostral at top
    direction: str


def _normalize_display(regvol: np.ndarray) -> np.ndarray:
    data = regvol.astype(np.float32)
    positive = data[data > 0]
    if positive.size == 0:
        return data
    vmax = float(np.quantile(positive, 0.999))
    vmin = float(np.quantile(positive, 0.001))
    scaled = (data - vmin) / max(vmax - vmin, 1e-6)
    return np.clip(scaled, 0, 1)


def _longitudinal_max_projection(volume_yxz: np.ndarray) -> np.ndarray:
    """(Y, X, Z) volume → (Z, Y) max along X with rostral (low Z) at the top."""
    return _normalize_display(volume_yxz).max(axis=1).T


def _match_transverse_width(projection: np.ndarray, target_width: int) -> np.ndarray:
    """Resample projection columns so sample and atlas panels share the same width."""
    if projection.shape[1] == target_width:
        return projection
    return resize(
        projection,
        (projection.shape[0], target_width),
        order=1,
        preserve_range=True,
        anti_aliasing=True,
    ).astype(np.float32)


def _primary_channel_index(config: SpinalCordPipelineConfig, n_channels: int) -> int:
    regchan = config.registration.channel_primary
    if n_channels == 1 or regchan > n_channels:
        return 1
    return int(regchan)


def _longitudinal_last(volume_yxz: np.ndarray) -> np.ndarray:
    """Permute so the longest (rostrocaudal) axis is last, matching preprocess."""
    im = int(np.argmax(volume_yxz.shape))
    axperm = [i for i in range(3) if i != im] + [im]
    return np.transpose(volume_yxz, axperm)


def _load_atlas_template(config: SpinalCordPipelineConfig) -> np.ndarray:
    """Return atlas template on the registration grid (Y, X, Z)."""
    save_path = cord_save_path(config)
    cached = cord_cache_dir(config) / "atlas_template.tif"
    if cached.is_file():
        return tifffile.imread(cached).astype(np.float32)
    regopts_path = save_path / "regopts.json"
    if regopts_path.is_file():
        checkpoint = CordRegOptsCheckpoint.load(regopts_path)
        if checkpoint.tv_path and Path(checkpoint.tv_path).is_file():
            return tifffile.imread(checkpoint.tv_path).astype(np.float32)
    volumes = load_fiederling_atlas_volumes(config.atlas)
    tv, _ = resize_fiederling_atlas(volumes, config.registration.resolution_um)
    try:
        cached.parent.mkdir(parents=True, exist_ok=True)
        tifffile.imwrite(cached, tv.astype(np.float32))
    except OSError:
        pass
    return tv


def _orientation_preview_path(config: SpinalCordPipelineConfig) -> Path:
    return cord_cache_dir(config) / ORIENTATION_PREVIEW_FILENAME


def _orientation_preview_key(config: SpinalCordPipelineConfig) -> dict[str, Any]:
    return {
        "register_fingerprint": manifest_register_fingerprint(config),
        "atlas_dir": str(config.atlas.atlas_dir.expanduser().resolve()),
        "channel_primary": int(config.registration.channel_primary),
        "resolution_um": float(config.registration.resolution_um),
    }


def _resolve_stored_direction(config: SpinalCordPipelineConfig) -> str:
    save_path = cord_save_path(config)
    stored = load_cord_orientation(save_path)
    if stored is None and config.registration.longitudinal_direction is not None:
        stored = config.registration.longitudinal_direction
    return stored if stored in _VALID_DIRECTIONS else ROSTROCAUDAL


def _load_orientation_preview(config: SpinalCordPipelineConfig) -> CordOrientationData | None:
    """Load cached longitudinal projections when register cache and atlas are unchanged."""
    path = _orientation_preview_path(config)
    if not path.is_file():
        return None
    if manifest_register_fingerprint(config) is None:
        return None
    try:
        with np.load(path, allow_pickle=False) as archive:
            stored_key = json.loads(str(archive["fingerprint_json"]))
            if stored_key != _orientation_preview_key(config):
                return None
            return CordOrientationData(
                sample_longitudinal=archive["sample_longitudinal"].astype(np.float32),
                atlas_longitudinal=archive["atlas_longitudinal"].astype(np.float32),
                direction=str(archive["direction"]),
            )
    except (OSError, ValueError, KeyError, json.JSONDecodeError):
        return None


def _save_orientation_preview(config: SpinalCordPipelineConfig, data: CordOrientationData) -> None:
    path = _orientation_preview_path(config)
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        path,
        fingerprint_json=json.dumps(_orientation_preview_key(config), sort_keys=True),
        sample_longitudinal=data.sample_longitudinal.astype(np.float32),
        atlas_longitudinal=data.atlas_longitudinal.astype(np.float32),
        direction=np.array(data.direction),
    )


def _build_orientation_data_from_volumes(config: SpinalCordPipelineConfig) -> CordOrientationData:
    registration = load_or_cache_cord_registration(config)
    regchan = _primary_channel_index(config, registration.n_channels)
    regvol = _longitudinal_last(registration.volume[:, :, :, regchan - 1].astype(np.float32))

    atlas_vol = _load_atlas_template(config)
    atlas_longitudinal = _longitudinal_max_projection(atlas_vol)
    sample_longitudinal = _match_transverse_width(
        _longitudinal_max_projection(regvol),
        atlas_longitudinal.shape[1],
    )
    return CordOrientationData(
        sample_longitudinal=sample_longitudinal,
        atlas_longitudinal=atlas_longitudinal,
        direction=_resolve_stored_direction(config),
    )


def cord_orientation_missing(config: SpinalCordPipelineConfig) -> bool:
    """True when preprocess must collect orientation (no YAML value or saved file)."""
    if config.registration.longitudinal_direction is not None:
        return False
    return load_cord_orientation(cord_save_path(config)) is None


def ensure_cord_orientation(
    config: SpinalCordPipelineConfig,
    *,
    headless: bool = False,
) -> str:
    """Resolve cord orientation, opening Napari when it has not been set yet."""
    save_path = cord_save_path(config)
    if not headless and cord_orientation_missing(config):
        console.print(
            "[bold]Cord orientation required[/bold] — opening Napari to set the "
            "rostrocaudal direction (writes cord_orientation.txt)."
        )
        run_spinal_orientation(config, headless=False)
    return resolve_cord_orientation(
        save_path,
        config_direction=config.registration.longitudinal_direction,
        require=True,
    )


def load_cord_orientation_check_data(config: SpinalCordPipelineConfig) -> CordOrientationData:
    """Load sample + atlas longitudinal max projections for orientation picking."""
    preview = _load_orientation_preview(config)
    if preview is not None:
        preview.direction = _resolve_stored_direction(config)
        emit_pipeline_message(
            "Using cached orientation preview "
            f"({ORIENTATION_PREVIEW_FILENAME}; inputs unchanged)."
        )
        return preview

    emit_pipeline_message(
        "Check-orientation: building longitudinal max projections from sample and atlas…"
    )
    data = _build_orientation_data_from_volumes(config)
    _save_orientation_preview(config, data)
    return data


def _build_spinal_orientation_controller(
    viewer: Any,
    config: SpinalCordPipelineConfig,
    data: CordOrientationData,
) -> DockStageController:
    """Attach cord orientation controls to an existing napari viewer."""
    from napari.utils.notifications import show_info
    from qtpy.QtCore import QTimer

    magicgui = require_magicgui()
    save_path = cord_save_path(config)

    def _sample_projection(direction: str) -> np.ndarray:
        if direction == CAUDOROSTRAL:
            return np.flip(data.sample_longitudinal, axis=0)
        return data.sample_longitudinal

    viewer.dims.ndisplay = 2
    atlas_layer = viewer.add_image(
        data.atlas_longitudinal,
        name="atlas template (rostral at top)",
        colormap="magma",
    )
    sample_layer = viewer.add_image(
        _sample_projection(data.direction),
        name="sample",
        colormap="gray",
    )

    def _layout_panels() -> None:
        _h, w = atlas_layer.data.shape
        atlas_layer.translate = (0.0, 0.0)
        sample_layer.translate = (0.0, float(w + PANEL_GAP_X))

    def _set_direction(direction: str) -> None:
        data.direction = direction
        sample_layer.data = _sample_projection(direction)
        sample_layer.name = f"sample ({direction}; rostral at top)"
        _layout_panels()
        viewer.status = (
            f"direction={direction} | tofliprc={'true' if direction == CAUDOROSTRAL else 'false'} "
            "| align rostral ends at the top of both panels"
        )

    @magicgui(call_button="Rostrocaudal (rostral = top)")
    def rostrocaudal_button() -> None:
        _set_direction(ROSTROCAUDAL)

    @magicgui(call_button="Caudorostral (flip sample)")
    def caudorostral_button() -> None:
        _set_direction(CAUDOROSTRAL)

    @magicgui(call_button="Save orientation && close")
    def save_controls() -> None:
        path = save_cord_orientation(save_path, data.direction, source="manual")
        show_info(f"Saved {path} (direction={data.direction})")
        QTimer.singleShot(0, lambda: close_stage_or_viewer(viewer))

    _set_direction(data.direction)

    return DockStageController(
        dock_widgets=[
            (rostrocaudal_button, "Rostrocaudal"),
            (caudorostral_button, "Caudorostral"),
            (save_controls, "Save"),
        ],
        _refresh_fn=lambda: _set_direction(data.direction),
        result=save_path,
    )


@dataclass
class _AsyncCordOrientationController:
    """Load orientation projections off the Qt main thread, then mount controls."""

    viewer: Any
    config: SpinalCordPipelineConfig
    on_log: Callable[[str], None] | None = None
    cancel_event: threading.Event | None = None
    _inner: DockStageController | None = field(default=None, init=False, repr=False)
    _workers: list[Any] = field(default_factory=list, init=False, repr=False)
    _cancelled: bool = field(default=False, init=False, repr=False)
    _on_busy_finished: Callable[[], None] | None = field(default=None, init=False, repr=False)
    open_log_message: str | None = field(
        default="Loading orientation preview in background (see log for progress)…",
        init=False,
    )
    result: Path | None = field(default=None, init=False)

    def __post_init__(self) -> None:
        self.result = cord_save_path(self.config)

    def register_busy_finished(self, callback: Callable[[], None]) -> None:
        """Notify the GUI shell when background loading finishes."""
        self._on_busy_finished = callback

    def _finish_busy(self) -> None:
        callback = self._on_busy_finished
        self._on_busy_finished = None
        if callback is not None:
            callback()

    def mount(self, viewer: Any) -> None:
        if self._inner is not None:
            self._inner.mount(viewer)
            return
        viewer.status = "Loading orientation preview… (see log for progress)"
        from lightsuite.gui.qt_workers import start_background_task

        reporter = CallbackReporter(on_message=self.on_log) if self.on_log else None

        def _work() -> CordOrientationData:
            with stage_cancellation(self.cancel_event):
                with capture_pipeline_output(reporter):
                    return load_cord_orientation_check_data(self.config)

        def _on_success(data: CordOrientationData) -> None:
            if self._cancelled:
                return
            self._inner = _build_spinal_orientation_controller(viewer, self.config, data)
            self._inner.mount(viewer)
            if self.on_log is not None:
                self.on_log("Orientation preview ready.")
            self._finish_busy()

        def _on_failure(exc: BaseException) -> None:
            if self._cancelled:
                return
            if isinstance(exc, StageCancelledError):
                if self.on_log is not None:
                    self.on_log("Stage cancelled.")
            else:
                message = f"Failed to load orientation preview: {exc}"
                if self.on_log is not None:
                    self.on_log(message)
                try:
                    from napari.utils.notifications import show_warning

                    show_warning(message)
                except ImportError:
                    pass
            self._finish_busy()

        self._workers.append(
            start_background_task(_work, on_success=_on_success, on_failure=_on_failure)
        )

    def refresh(self) -> None:
        if self._inner is not None:
            self._inner.refresh()

    def teardown(self, viewer: Any) -> None:
        self._cancelled = True
        if self.cancel_event is not None:
            self.cancel_event.set()
        from lightsuite.gui.qt_workers import cancel_background_workers

        cancel_background_workers(self._workers)
        self._workers.clear()
        if self._inner is not None:
            self._inner.teardown(viewer)
        self._finish_busy()


def attach_spinal_orientation(
    viewer: Any,
    config: SpinalCordPipelineConfig,
    *,
    ctx: Any | None = None,
    data: CordOrientationData | None = None,
) -> DockStageController | _AsyncCordOrientationController:
    """Attach cord orientation controls to an existing napari viewer."""
    if data is not None:
        return _build_spinal_orientation_controller(viewer, config, data)
    on_log = getattr(ctx, "on_log", None) if ctx is not None else None
    cancel_event = getattr(ctx, "cancel_event", None) if ctx is not None else None
    return _AsyncCordOrientationController(
        viewer=viewer,
        config=config,
        on_log=on_log,
        cancel_event=cancel_event,
    )


def run_spinal_orientation(
    config: SpinalCordPipelineConfig,
    *,
    headless: bool = False,
) -> Path:
    """Open Napari to pick the rostral end; write cord_orientation.txt."""
    save_path = cord_save_path(config)
    if headless:
        stored = load_cord_orientation(save_path)
        if stored is None and config.registration.longitudinal_direction is not None:
            stored = config.registration.longitudinal_direction
        direction = stored if stored in _VALID_DIRECTIONS else ROSTROCAUDAL
        return save_cord_orientation(save_path, direction, source="headless")

    data = load_cord_orientation_check_data(config)
    title = f"LightSuite cord orientation — {config.sample.name}"

    def _attach(viewer: Any) -> DockStageController:
        return attach_spinal_orientation(viewer, config, data=data)

    run_attached_stage(
        title,
        _attach,
        before_run=lambda: console.print(
            "[bold]Cord orientation[/bold] — atlas template (left) and sample (right) longitudinal "
            "max projections, aligned side by side. Click Rostrocaudal or Caudorostral to flip the "
            "sample until rostral anatomy matches the atlas, then Save. "
            f"Writes cord_orientation.txt in {save_path}."
        ),
    )
    stored = load_cord_orientation(save_path)
    if stored is None:
        msg = (
            f"Orientation not saved. Use Save in the Napari window to write "
            f"{cord_orientation_path(save_path)}."
        )
        raise RuntimeError(msg)
    return cord_orientation_path(save_path)
