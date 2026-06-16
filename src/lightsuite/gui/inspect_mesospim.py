"""Napari viewer for mesoSPIM registration QA (overview vs full-canvas ROI)."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from lightsuite.mesospim.checkpoint import MesospimRegOptsCheckpoint, mesospim_checkpoint_path
from lightsuite.mesospim.config_models import MesospimPipelineConfig
from lightsuite.mesospim.io import load_registered_canvas_zyx, load_tiff_zyx_volume
from lightsuite.mesospim.registration import sanitize_experiment_name


@dataclass(frozen=True)
class MesospimInspectPaths:
    overview_path: Path
    registered_full_overview_path: Path


def _default_registered_full_overview_path(cfg: MesospimPipelineConfig) -> Path:
    meso = cfg.mesospim
    slug = sanitize_experiment_name(meso.registration.experiment_name)
    overview_stem = meso.overview.path.stem
    if overview_stem.lower().endswith(".tif"):
        overview_stem = overview_stem[:-4]
    roi_stem = meso.roi.path.stem
    if roi_stem.lower().endswith(".tif"):
        roi_stem = roi_stem[:-4]
    return (
        cfg.sample.save_path
        / "elastix_roi_to_overview"
        / slug
        / f"{slug}_{roi_stem}_registered_to_{overview_stem}_in_full_overview.tif"
    )


def resolve_mesospim_inspect_paths(cfg: MesospimPipelineConfig) -> MesospimInspectPaths:
    """Resolve overview and full-canvas registration TIFF paths."""
    overview_path = cfg.mesospim.overview.path.expanduser().resolve()
    if not overview_path.is_file():
        msg = f"Overview TIFF not found: {overview_path}"
        raise FileNotFoundError(msg)

    checkpoint_path = mesospim_checkpoint_path(cfg.sample.save_path)
    registered_path: Path | None = None
    if checkpoint_path.is_file():
        checkpoint = MesospimRegOptsCheckpoint.load(checkpoint_path)
        if checkpoint.registered_roi_full_overview_path:
            registered_path = Path(checkpoint.registered_roi_full_overview_path).expanduser()

    if registered_path is None or not registered_path.is_file():
        registered_path = _default_registered_full_overview_path(cfg)

    if not registered_path.is_file():
        msg = (
            f"Registered full-overview canvas not found: {registered_path}\n"
            "Run 'lightsuite mesospim register' with write_full_overview_canvas: true first."
        )
        raise FileNotFoundError(msg)

    return MesospimInspectPaths(
        overview_path=overview_path,
        registered_full_overview_path=registered_path.resolve(),
    )


def _contrast_limits(volume: np.ndarray) -> tuple[float, float]:
    positive = volume[volume > 0]
    if positive.size:
        lo, hi = np.percentile(positive, (1.0, 99.5))
    else:
        lo, hi = float(np.min(volume)), float(np.max(volume))
    if hi <= lo:
        hi = lo + 1.0
    return float(lo), float(hi)


def load_mesospim_inspect_volumes(cfg: MesospimPipelineConfig) -> tuple[np.ndarray, np.ndarray]:
    """Load overview and registered full-canvas volumes as ZYX float32 arrays."""
    paths = resolve_mesospim_inspect_paths(cfg)
    meso = cfg.mesospim
    overview = load_tiff_zyx_volume(
        paths.overview_path,
        overview_path=meso.overview.path,
        roi_path=meso.roi.path,
        remap=meso.tiff_remap,
    ).astype(np.float32, copy=False)

    registered = load_registered_canvas_zyx(paths.registered_full_overview_path)

    if overview.shape != registered.shape:
        msg = (
            f"Shape mismatch: overview ZYX {overview.shape} vs "
            f"registered canvas {registered.shape}"
        )
        raise ValueError(msg)

    return overview, registered


def run_mesospim_inspect(cfg: MesospimPipelineConfig, *, headless: bool = False) -> MesospimInspectPaths:
    """Open Napari to compare the 1× overview and embedded registration canvas."""
    paths = resolve_mesospim_inspect_paths(cfg)
    if headless:
        return paths

    try:
        import napari
        from napari.utils.notifications import show_info
    except ImportError as exc:
        msg = "Napari GUI requires: uv sync --extra gui"
        raise RuntimeError(msg) from exc

    overview, registered = load_mesospim_inspect_volumes(cfg)
    overview_limits = _contrast_limits(overview)
    registered_limits = _contrast_limits(registered)

    viewer = napari.Viewer(title=f"LightSuite mesospim — {cfg.sample.name}")
    viewer.add_image(
        overview,
        name="overview 1×",
        colormap="gray",
        blending="opaque",
        contrast_limits=overview_limits,
        scale=(1.0, 1.0, 1.0),
    )
    viewer.add_image(
        registered,
        name="registered ROI (full 1× canvas)",
        colormap="green",
        blending="additive",
        opacity=0.55,
        contrast_limits=registered_limits,
        scale=(1.0, 1.0, 1.0),
    )

    show_info(
        "Toggle layer visibility or opacity. Both volumes share the 1× overview grid; "
        "non-overlap voxels in the registered layer are zero."
    )
    napari.run()
    return paths
