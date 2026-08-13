"""Napari GUI to manually set spinal cord rostrocaudal orientation."""

from __future__ import annotations

from dataclasses import dataclass
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
from lightsuite.io.cord_registration_cache import load_or_cache_cord_registration
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
    return tv


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
    registration = load_or_cache_cord_registration(config)
    regchan = _primary_channel_index(config, registration.n_channels)
    regvol = _longitudinal_last(registration.volume[:, :, :, regchan - 1].astype(np.float32))

    atlas_vol = _load_atlas_template(config)
    atlas_longitudinal = _longitudinal_max_projection(atlas_vol)
    sample_longitudinal = _match_transverse_width(
        _longitudinal_max_projection(regvol),
        atlas_longitudinal.shape[1],
    )

    save_path = cord_save_path(config)
    stored = load_cord_orientation(save_path)
    if stored is None and config.registration.longitudinal_direction is not None:
        stored = config.registration.longitudinal_direction
    direction = stored if stored in _VALID_DIRECTIONS else ROSTROCAUDAL

    return CordOrientationData(
        sample_longitudinal=sample_longitudinal,
        atlas_longitudinal=atlas_longitudinal,
        direction=direction,
    )


def attach_spinal_orientation(
    viewer: Any,
    config: SpinalCordPipelineConfig,
    *,
    data: CordOrientationData | None = None,
) -> DockStageController:
    """Attach cord orientation controls to an existing napari viewer."""
    from napari.utils.notifications import show_info
    from qtpy.QtCore import QTimer

    magicgui = require_magicgui()
    save_path = cord_save_path(config)
    if data is None:
        data = load_cord_orientation_check_data(config)

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

    return DockStageController(
        dock_widgets=[
            (rostrocaudal_button, "Rostrocaudal"),
            (caudorostral_button, "Caudorostral"),
            (save_controls, "Save"),
        ],
        _refresh_fn=lambda: _set_direction(data.direction),
        result=save_path,
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
