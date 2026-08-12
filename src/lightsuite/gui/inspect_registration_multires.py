"""Napari viewer comparing the registered ROI against the overview after multires register."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
from rich.console import Console

from lightsuite.gui.stage_controller import DockStageController, run_attached_stage
from lightsuite.multires.checkpoint import MultiresRegOptsCheckpoint, multires_checkpoint_path
from lightsuite.multires.config_models import MultiresPipelineConfig

console = Console()

_GUI_HINT = "uv sync --extra gui --extra registration"


@dataclass(frozen=True)
class MultiresRegistrationInspectPaths:
    """Resolved layers for registered ROI vs overview QC."""

    checkpoint_path: Path
    pair_label: str
    overview_path: Path
    registered_roi_paths: dict[str, Path] = field(default_factory=dict)
    full_overview: bool = False


def discover_multires_registration_inspect_paths(
    config: MultiresPipelineConfig,
    *,
    full_overview: bool = False,
) -> MultiresRegistrationInspectPaths:
    """Find the overview and registered ROI volumes written by ``multires register``.

    Defaults to the overlap crop, where the ROI and overview share one small grid.
    ``full_overview`` switches to the full-overview canvas instead.
    """
    save_path = config.sample.save_path.expanduser()
    checkpoint_path = multires_checkpoint_path(save_path)
    if not checkpoint_path.is_file():
        msg = f"Missing {checkpoint_path}. Run 'lightsuite multires register' first."
        raise FileNotFoundError(msg)

    checkpoint = MultiresRegOptsCheckpoint.load(checkpoint_path)

    if full_overview:
        overview_path = Path(checkpoint.overview_volume_path).expanduser()
        roi_path_str = checkpoint.registered_roi_full_overview_path
        if not roi_path_str:
            msg = (
                "No full-overview canvas in multires_regopts.json. Re-run register with "
                "multires.registration.write_full_overview_canvas: true, or drop --full-overview."
            )
            raise FileNotFoundError(msg)
    else:
        if not checkpoint.cropped_overview_path:
            msg = "multires_regopts.json missing cropped_overview_path. Re-run 'multires register'."
            raise FileNotFoundError(msg)
        overview_path = Path(checkpoint.cropped_overview_path).expanduser()
        roi_path_str = checkpoint.registered_roi_path

    if not roi_path_str:
        msg = "multires_regopts.json missing registered_roi_path. Re-run 'multires register'."
        raise FileNotFoundError(msg)
    if not overview_path.is_file():
        msg = f"Missing overview volume {overview_path}. Re-run 'multires register'."
        raise FileNotFoundError(msg)

    reference_channel = checkpoint.reference_channel or "reference"
    registered: dict[str, Path] = {}
    roi_path = Path(roi_path_str).expanduser()
    if roi_path.is_file():
        registered[reference_channel] = roi_path

    # Extra channels only exist as overlap crops; skip them in full-overview mode.
    if not full_overview:
        for channel, path_str in (checkpoint.additional_channel_paths or {}).items():
            path = Path(path_str).expanduser()
            if path.is_file():
                registered[str(channel)] = path

    if not registered:
        msg = (
            f"No registered ROI volumes found for {checkpoint.pair_label}. "
            "Re-run 'lightsuite multires register'."
        )
        raise FileNotFoundError(msg)

    return MultiresRegistrationInspectPaths(
        checkpoint_path=checkpoint_path.resolve(),
        pair_label=checkpoint.pair_label,
        overview_path=overview_path.resolve(),
        registered_roi_paths={k: v.resolve() for k, v in registered.items()},
        full_overview=full_overview,
    )


def load_inspect_volume(path: Path) -> np.ndarray:
    """Load a registration TIFF lazily where possible (canvases can be tens of GB)."""
    from lightsuite.registration.volume import load_tiff_volume_zyx

    return load_tiff_volume_zyx(path)


def contrast_limits(volume: np.ndarray, *, max_sample: int = 4_000_000) -> tuple[float, float]:
    """Percentile contrast from a strided sample so memmapped canvases stay lazy."""
    flat = np.asarray(volume).reshape(-1)
    step = max(1, flat.size // max_sample)
    sample = np.asarray(flat[::step], dtype=np.float32)
    positive = sample[sample > 0]
    if positive.size:
        lo, hi = np.percentile(positive, (1.0, 99.5))
    else:
        lo, hi = float(sample.min()), float(sample.max())
    if hi <= lo:
        hi = lo + 1.0
    return float(lo), float(hi)


def attach_multires_inspect_registration(
    viewer: Any,
    config: MultiresPipelineConfig,
    paths: MultiresRegistrationInspectPaths,
) -> DockStageController:
    """Attach registered ROI vs overview layers to an existing napari viewer."""
    from napari.utils.notifications import show_info

    grid = "full overview" if paths.full_overview else "overlap crop"

    overview = load_inspect_volume(paths.overview_path)
    viewer.add_image(
        overview,
        name="overview",
        colormap="gray",
        blending="additive",
        contrast_limits=contrast_limits(overview),
    )

    channel_cmaps = ["magenta", "cyan", "yellow", "green"]
    for idx, (channel, path) in enumerate(sorted(paths.registered_roi_paths.items())):
        vol = load_inspect_volume(path)
        if vol.shape != overview.shape:
            console.print(
                f"[yellow]Channel {channel} shape {vol.shape} != overview {overview.shape}; "
                "layers may not align.[/yellow]"
            )
        viewer.add_image(
            vol,
            name=f"ROI registered — channel {channel}",
            colormap=channel_cmaps[idx % len(channel_cmaps)],
            blending="additive",
            opacity=0.6,
            contrast_limits=contrast_limits(vol),
        )

    def _notify() -> None:
        show_info(
            f"Loaded overview + {len(paths.registered_roi_paths)} registered ROI channel(s) "
            f"on the {grid} grid."
        )

    return DockStageController(
        dock_widgets=[],
        _refresh_fn=_notify,
        result=paths,
    )


def run_multires_inspect_registration(
    config: MultiresPipelineConfig,
    *,
    full_overview: bool = False,
    headless: bool = False,
) -> MultiresRegistrationInspectPaths:
    """Open Napari with the overview and every registered ROI channel overlaid."""
    paths = discover_multires_registration_inspect_paths(config, full_overview=full_overview)
    if headless:
        load_inspect_volume(paths.overview_path)
        for path in paths.registered_roi_paths.values():
            load_inspect_volume(path)
        return paths

    grid = "full overview" if paths.full_overview else "overlap crop"
    title = (
        f"LightSuite multires registration QC — "
        f"{config.sample.name} / {paths.pair_label} ({grid})"
    )

    def _attach(viewer: Any) -> DockStageController:
        return attach_multires_inspect_registration(viewer, config, paths)

    final = run_attached_stage(title, _attach)
    return final if final is not None else paths


__all__ = [
    "MultiresRegistrationInspectPaths",
    "discover_multires_registration_inspect_paths",
    "run_multires_inspect_registration",
]
