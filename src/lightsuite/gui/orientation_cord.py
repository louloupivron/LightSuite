"""Napari GUI to manually set spinal cord rostrocaudal orientation."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import tifffile
from rich.console import Console
from skimage.transform import resize

from lightsuite.atlas.fiederling import load_fiederling_atlas_volumes, resize_fiederling_atlas
from lightsuite.config.models import SpinalCordPipelineConfig
from lightsuite.io.cord_reader import read_spinal_cord_sample
from lightsuite.preprocess.cord_checkpoint import CordRegOptsCheckpoint
from lightsuite.registration.cord_orientation import (
    CAUDOROSTRAL,
    ROSTROCAUDAL,
    load_cord_orientation,
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


def load_cord_orientation_check_data(config: SpinalCordPipelineConfig) -> CordOrientationData:
    """Load sample + atlas longitudinal max projections for orientation picking."""
    save_path = cord_save_path(config)
    regopts_path = save_path / "regopts.json"

    regvol: np.ndarray | None = None
    if regopts_path.is_file():
        checkpoint = CordRegOptsCheckpoint.load(regopts_path)
        regchan = _primary_channel_index(config, checkpoint.nchans)
        cached = (checkpoint.regvolpaths or {}).get(str(regchan))
        if cached and Path(cached).is_file():
            regvol = _longitudinal_last(tifffile.imread(cached).astype(np.float32))

    if regvol is None:
        sample = read_spinal_cord_sample(config)
        regchan = _primary_channel_index(config, sample.n_channels)
        regvol = _longitudinal_last(sample.volume[:, :, :, regchan - 1].astype(np.float32))

    atlas_vol = _load_atlas_template(config)
    atlas_longitudinal = _longitudinal_max_projection(atlas_vol)
    sample_longitudinal = _match_transverse_width(
        _longitudinal_max_projection(regvol),
        atlas_longitudinal.shape[1],
    )

    stored = load_cord_orientation(save_path)
    if stored is None and config.registration.longitudinal_direction is not None:
        stored = config.registration.longitudinal_direction
    direction = stored if stored in _VALID_DIRECTIONS else ROSTROCAUDAL

    return CordOrientationData(
        sample_longitudinal=sample_longitudinal,
        atlas_longitudinal=atlas_longitudinal,
        direction=direction,
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

    try:
        import napari
        from magicgui import magicgui
        from napari.utils.notifications import show_info
        from qtpy.QtCore import QTimer
    except ImportError as exc:
        msg = "Napari GUI requires: uv sync --extra gui"
        raise RuntimeError(msg) from exc

    state = {"direction": data.direction}

    def _sample_projection(direction: str) -> np.ndarray:
        if direction == CAUDOROSTRAL:
            return np.flip(data.sample_longitudinal, axis=0)
        return data.sample_longitudinal

    def _layout_panels() -> None:
        _h, w = atlas_layer.data.shape
        atlas_layer.translate = (0.0, 0.0)
        sample_layer.translate = (0.0, float(w + PANEL_GAP_X))

    viewer = napari.Viewer(title=f"LightSuite cord orientation — {config.sample.name}")
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
    _layout_panels()

    def _set_direction(direction: str) -> None:
        state["direction"] = direction
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
        path = save_cord_orientation(save_path, state["direction"], source="manual")
        show_info(f"Saved {path} (direction={state['direction']})")
        QTimer.singleShot(0, viewer.close)

    viewer.window.add_dock_widget(rostrocaudal_button, area="right", name="Rostrocaudal")
    viewer.window.add_dock_widget(caudorostral_button, area="right", name="Caudorostral")
    viewer.window.add_dock_widget(save_controls, area="right", name="Save")
    _set_direction(data.direction)

    console.print(
        "[bold]Cord orientation[/bold] — atlas template (left) and sample (right) longitudinal "
        "max projections, aligned side by side. Click Rostrocaudal or Caudorostral to flip the "
        "sample until rostral anatomy matches the atlas, then Save. "
        f"Writes cord_orientation.txt in {save_path}."
    )
    napari.run()
    return save_cord_orientation(save_path, state["direction"], source="manual")
