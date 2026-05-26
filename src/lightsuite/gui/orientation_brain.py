"""Napari GUI for brain orientation checking (getBrainOrientation.m / brainreg)."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import nibabel as nib
import numpy as np
from rich.console import Console

from lightsuite.atlas.registry import resolve_brain_atlas
from lightsuite.config.models import BrainPipelineConfig
from lightsuite.gui.orientation_views import build_dual_panel_projections, mean_projections
from lightsuite.preprocess.checkpoint import RegOptsCheckpoint
from lightsuite.registration.orientation import (
    DEFAULT_PERMVEC,
    PERMUTATION_OPTIONS,
    indices_from_permvec,
    load_orientation_file,
    orientation_path,
    permute_for_atlas,
    permvec_from_indices,
    save_orientation,
    validate_permvec,
)
from lightsuite.registration.volume import load_registration_volume, resize_atlas_volume

console = Console()

PANEL_GAP_X = 24


@dataclass
class OrientationCheckData:
    sample_volume: np.ndarray
    atlas_volume: np.ndarray
    permvec: list[int]
    orientation_file: Path


def load_orientation_check_data(config: BrainPipelineConfig) -> OrientationCheckData:
    save_path = config.sample.save_path.expanduser()
    regopts_path = save_path / "regopts.json"
    if not regopts_path.is_file():
        msg = f"Missing {regopts_path}. Run 'lightsuite brain preprocess' first."
        raise FileNotFoundError(msg)

    checkpoint = RegOptsCheckpoint.load(regopts_path)
    sample = load_registration_volume(Path(checkpoint.regvolpath)).astype(np.float32)

    atlas = resolve_brain_atlas(config.atlas.provider.value, config.atlas.atlas_dir)
    template = np.asanyarray(nib.load(atlas.template_path).dataobj).astype(np.float32)
    downfac = config.atlas.resolution_um / checkpoint.registres_um
    atlas_reg = resize_atlas_volume(template, downfac, nearest=False)

    orient_file = orientation_path(save_path)
    if config.registration.orientation is not None:
        permvec = list(config.registration.orientation)
    elif orient_file.is_file():
        permvec = load_orientation_file(orient_file)
    else:
        permvec = DEFAULT_PERMVEC.copy()
    validate_permvec(permvec)

    return OrientationCheckData(
        sample_volume=sample,
        atlas_volume=atlas_reg,
        permvec=permvec,
        orientation_file=orient_file,
    )


def prepare_orientation_session(config: BrainPipelineConfig) -> Path:
    """Ensure brain_orientation.txt exists without opening Napari (for tests)."""
    data = load_orientation_check_data(config)
    return save_orientation(config.sample.save_path, data.permvec)


def run_brain_orientation_check(config: BrainPipelineConfig, *, headless: bool = False) -> Path:
    """Open Napari orientation checker or write orientation file in headless mode."""
    data = load_orientation_check_data(config)
    if headless:
        return save_orientation(config.sample.save_path, data.permvec)

    try:
        import napari
        from magicgui import magicgui
        from napari.utils.notifications import show_info
    except ImportError as exc:
        msg = "Napari GUI requires: uv sync --extra gui"
        raise RuntimeError(msg) from exc

    option_labels = [label for _value, label in PERMUTATION_OPTIONS]
    idx0, idx1, idx2 = indices_from_permvec(data.permvec)

    viewer = napari.Viewer(title=f"LightSuite orientation — {config.sample.name}")

    atlas_layer = viewer.add_image(
        np.zeros((10, 10), dtype=np.float32),
        name="atlas (reference)",
        colormap="gray",
        blending="opaque",
    )
    sample_layer = viewer.add_image(
        np.zeros((10, 10), dtype=np.float32),
        name="sample (input)",
        colormap="gray",
        blending="opaque",
    )

    def _show_projections(permvec: list[int]) -> None:
        sample_oriented = permute_for_atlas(data.sample_volume, permvec)
        atlas_proj = mean_projections(data.atlas_volume)
        sample_proj = mean_projections(sample_oriented)
        atlas_panel, sample_panel, sample_x = build_dual_panel_projections(
            atlas_proj,
            sample_proj,
            gap_x=PANEL_GAP_X,
        )

        atlas_layer.data = atlas_panel
        sample_layer.data = sample_panel
        sample_layer.translate = (sample_x, 0.0)
        sample_layer.scale = (1.0, 1.0)

        viewer.status = (
            f"permvec={permvec} | atlas {atlas_panel.shape} px, "
            f"sample {sample_panel.shape} px (native, 3 projections stacked)"
        )

    @magicgui(
        atlas_dim_1={"choices": option_labels, "label": "Map to atlas dim 1 (Y)"},
        atlas_dim_2={"choices": option_labels, "label": "Map to atlas dim 2 (X)"},
        atlas_dim_3={"choices": option_labels, "label": "Map to atlas dim 3 (Z)"},
        call_button="Update preview",
    )
    def controls(
        atlas_dim_1: str = option_labels[idx0],
        atlas_dim_2: str = option_labels[idx1],
        atlas_dim_3: str = option_labels[idx2],
    ) -> None:
        indices = (
            option_labels.index(atlas_dim_1),
            option_labels.index(atlas_dim_2),
            option_labels.index(atlas_dim_3),
        )
        try:
            permvec = permvec_from_indices(indices)
        except ValueError as exc:
            show_info(str(exc))
            return
        data.permvec = permvec
        _show_projections(permvec)

    @magicgui(call_button="Save orientation && close")
    def save_controls() -> None:
        try:
            validate_permvec(data.permvec)
        except ValueError as exc:
            show_info(str(exc))
            return
        path = save_orientation(config.sample.save_path, data.permvec)
        show_info(f"Saved {path}")
        viewer.close()

    viewer.window.add_dock_widget(controls, area="right", name="Orientation")
    viewer.window.add_dock_widget(save_controls, area="right", name="Save")
    controls.atlas_dim_1.value = option_labels[idx0]
    controls.atlas_dim_2.value = option_labels[idx1]
    controls.atlas_dim_3.value = option_labels[idx2]
    _show_projections(data.permvec)

    console.print(
        "[bold]Orientation checker[/bold] — atlas on the left, permuted sample on the right. "
        "Three axis projections are stacked in each panel. Adjust dropdowns and Update preview."
    )
    napari.run()
    return data.orientation_file
