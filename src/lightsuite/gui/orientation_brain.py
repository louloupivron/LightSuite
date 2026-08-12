"""Napari GUI for brain orientation checking (getBrainOrientation.m / brainreg)."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
from rich.console import Console

from lightsuite.atlas.io import load_atlas_volume
from lightsuite.atlas.registry import resolve_brain_atlas_content
from lightsuite.config.loader import save_orientation_to_config
from lightsuite.config.models import BrainPipelineConfig
from lightsuite.gui.brain_data import _normalize_display
from lightsuite.preprocess.checkpoint import RegOptsCheckpoint
from lightsuite.registration.orientation import (
    DEFAULT_PERMVEC,
    PERMUTATION_OPTIONS,
    indices_from_permvec,
    load_orientation_file,
    orientation_path,
    permute_for_atlas,
    permvec_from_indices,
    validate_permvec,
)
from lightsuite.registration.volume import load_registration_volume, resize_atlas_volume
from skimage.transform import resize

from lightsuite.gui.stage_controller import (
    DockStageController,
    close_stage_or_viewer,
    require_magicgui,
    require_napari,
    run_attached_stage,
)

console = Console()

PANEL_GAP_X = 24
# Napari holds multiple float32 copies of each layer; cap preview volumes to avoid OOM.
ORIENTATION_PREVIEW_MAX_BYTES = 256 * 1024 * 1024


@dataclass
class OrientationCheckData:
    sample_volume: np.ndarray
    atlas_volume: np.ndarray
    permvec: list[int]


def _downsample_for_orientation_preview(
    volume: np.ndarray,
    *,
    max_bytes: int = ORIENTATION_PREVIEW_MAX_BYTES,
) -> np.ndarray:
    """Shrink a volume for Napari when the full registration grid is too large."""
    vol = np.asarray(volume, dtype=np.float32)
    if vol.nbytes <= max_bytes:
        return vol
    scale = (max_bytes / vol.nbytes) ** (1.0 / 3.0)
    new_shape = tuple(max(1, int(dim * scale)) for dim in vol.shape)
    return resize(
        vol,
        new_shape,
        order=1,
        preserve_range=True,
        anti_aliasing=True,
    ).astype(np.float32)


def _atlas_for_orientation_check(template: np.ndarray, *, downfac: float) -> np.ndarray:
    """Atlas volume for the orientation GUI (preview only).

  When the atlas is coarser than the registration grid (``downfac > 1``), upscaling
  to registration voxels is required later in init-registration but allocates far too
  much memory for Napari (e.g. Waxholm 39 µm on a 20 µm grid ≈ 8 GiB float32).
  Native atlas voxels are sufficient to verify axis permutations.
    """
    if downfac > 1.0 or np.isclose(downfac, 1.0):
        return template.astype(np.float32, copy=False)
    return resize_atlas_volume(template.astype(np.float32), downfac, nearest=False)


def _atlas_for_display(atlas: np.ndarray) -> np.ndarray:
    """Return the atlas in its native voxel order for display.

    The atlas is the fixed reference: init-registration and match-points consume it
    in native NIfTI axis order, and the saved ``permvec`` maps the *sample* onto that
    order. Resizing the atlas to the permuted sample shape (as a previous version did)
    stretched it to a wrong aspect ratio and broke the orientation comparison badly for
    atlases whose native shape differs from the sample (e.g. Perens/gubra ``L,P,S``).
    """
    return atlas.astype(np.float32)


def load_orientation_check_data(config: BrainPipelineConfig) -> OrientationCheckData:
    save_path = config.sample.save_path.expanduser()
    regopts_path = save_path / "regopts.json"
    if not regopts_path.is_file():
        msg = f"Missing {regopts_path}. Run 'lightsuite brain preprocess' first."
        raise FileNotFoundError(msg)

    checkpoint = RegOptsCheckpoint.load(regopts_path)
    sample_full = load_registration_volume(Path(checkpoint.regvolpath)).astype(np.float32)

    atlas = resolve_brain_atlas_content(
        config.atlas,
        scratch=config.sample.scratch,
    ).paths
    template = load_atlas_volume(atlas.template_path).astype(np.float32)
    downfac = config.atlas.resolution_um / checkpoint.registres_um
    atlas_for_check = _atlas_for_orientation_check(template, downfac=downfac)
    sample = _downsample_for_orientation_preview(sample_full)
    atlas_reg = _downsample_for_orientation_preview(atlas_for_check)

    if sample.shape != sample_full.shape or atlas_reg.shape != atlas_for_check.shape:
        console.print(
            "[yellow]Orientation preview downsampled for Napari "
            f"(sample {sample_full.shape}→{sample.shape}, "
            f"atlas {atlas_for_check.shape}→{atlas_reg.shape}).[/yellow]"
        )
    if downfac > 1.0:
        console.print(
            "[dim]Atlas shown at native "
            f"{template.shape} @ {config.atlas.resolution_um:g} µm "
            f"(not upscaled ×{downfac:.2f} onto the {checkpoint.registres_um:g} µm grid "
            "for this preview).[/dim]"
        )

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
    )


def prepare_orientation_session(config: BrainPipelineConfig, config_path: Path) -> Path:
    """Write orientation to the pipeline YAML without opening Napari (for tests)."""
    data = load_orientation_check_data(config)
    return save_orientation_to_config(config_path, data.permvec)


def attach_brain_orientation_check(
    viewer: Any,
    config: BrainPipelineConfig,
    config_path: Path,
    *,
    data: OrientationCheckData | None = None,
) -> DockStageController:
    """Attach brain orientation controls to an existing napari viewer."""
    from napari.utils.notifications import show_info
    from qtpy.QtCore import QTimer

    magicgui = require_magicgui()
    if data is None:
        data = load_orientation_check_data(config)
    config_path = config_path.expanduser().resolve()

    option_labels = [label for _value, label in PERMUTATION_OPTIONS]
    idx0, idx1, idx2 = indices_from_permvec(data.permvec)

    atlas_layer = viewer.add_image(
        np.zeros((10, 10, 10), dtype=np.float32),
        name="atlas (reference)",
        colormap="gray",
        blending="opaque",
    )
    sample_layer = viewer.add_image(
        np.zeros((10, 10, 10), dtype=np.float32),
        name="sample (input)",
        colormap="gray",
        blending="opaque",
    )

    atlas_display = _normalize_display(_atlas_for_display(data.atlas_volume))
    atlas_layer.data = atlas_display

    def _update_preview(permvec: list[int]) -> None:
        sample_oriented = permute_for_atlas(data.sample_volume, permvec)
        sample_display = _normalize_display(sample_oriented)

        sample_layer.data = sample_display
        sample_layer.translate = (0.0, float(atlas_display.shape[1] + PANEL_GAP_X), 0.0)
        sample_layer.scale = (1.0, 1.0, 1.0)

        viewer.status = (
            f"permvec={permvec} | atlas {atlas_display.shape}, sample {sample_display.shape} "
            "(use Napari dimension slider and Ctrl+E to inspect axes)"
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
        _update_preview(permvec)

    @magicgui(call_button="Save orientation && close")
    def save_controls() -> None:
        try:
            validate_permvec(data.permvec)
        except ValueError as exc:
            show_info(str(exc))
            return
        path = save_orientation_to_config(config_path, data.permvec)
        show_info(f"Saved orientation to {path}")
        QTimer.singleShot(0, lambda: close_stage_or_viewer(viewer))

    controls.atlas_dim_1.value = option_labels[idx0]
    controls.atlas_dim_2.value = option_labels[idx1]
    controls.atlas_dim_3.value = option_labels[idx2]

    return DockStageController(
        dock_widgets=[
            (controls, "Orientation"),
            (save_controls, "Save"),
        ],
        _refresh_fn=lambda: _update_preview(data.permvec),
        result=config_path,
    )


def run_brain_orientation_check(
    config: BrainPipelineConfig,
    config_path: Path,
    *,
    headless: bool = False,
) -> Path:
    """Open Napari orientation checker or update config YAML in headless mode."""
    data = load_orientation_check_data(config)
    config_path = config_path.expanduser().resolve()
    if headless:
        return save_orientation_to_config(config_path, data.permvec)

    title = f"LightSuite orientation — {config.sample.name}"

    def _attach(viewer: Any) -> DockStageController:
        return attach_brain_orientation_check(viewer, config, config_path, data=data)

    run_attached_stage(
        title,
        _attach,
        before_run=lambda: console.print(
            "[bold]Orientation checker[/bold] — atlas on the left, permuted sample on the right. "
            "Scroll through slices and use Napari's axis-order control (Ctrl+E) to inspect projections. "
            "Adjust dropdowns and click Update preview."
        ),
    )
    return config_path
