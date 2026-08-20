"""Apply registration transforms to sample volumes (atlasSpaceFromVolumeParams.m)."""

from __future__ import annotations

import shutil
import tempfile
from pathlib import Path

import numpy as np

from lightsuite.registration.brain_register import TransformParamsCheckpoint
from lightsuite.registration.canvas import (
    RegistrationCanvas,
    WarpCanvasPadding,
    apply_canvas_sample_crop,
    pad_volume_for_warp_canvas,
)
from lightsuite.registration.coordinates import affine_with_source_offset
from lightsuite.registration.elastix.runner import run_transformix
from lightsuite.registration.volume import permute_brain_volume
from lightsuite.registration.warp import warp_volume_affine


def transform_volume_to_atlas(
    volume: np.ndarray,
    transform_params: TransformParamsCheckpoint,
    *,
    permute: list[int],
    spacing_mm: float,
    temp_dir: Path | None = None,
    nearest: bool = False,
) -> np.ndarray:
    """Warp a registration-resolution channel volume into atlas voxel space."""
    vol = permute_brain_volume(volume.astype(np.float32), permute)
    canvas = RegistrationCanvas.from_checkpoint_dict(transform_params.registration_canvas)
    warp_pad = WarpCanvasPadding((0, 0, 0), (0, 0, 0))
    if canvas is not None:
        # apply_canvas_sample_crop already pads to canvas.working_shape.
        vol = apply_canvas_sample_crop(vol, canvas)
        warp_pad = WarpCanvasPadding(canvas.pad_before, canvas.pad_after)
    elif transform_params.warp_canvas_pad_before is not None:
        # Legacy runs recorded only pad_before; the VD padding they applied was symmetric.
        pad = tuple(int(v) for v in transform_params.warp_canvas_pad_before)
        warp_pad = WarpCanvasPadding(pad, pad)
        vol = pad_volume_for_warp_canvas(vol, warp_pad)
    bspline_path = Path(transform_params.tform_bspline_samp20um_to_atlas_20um_px)
    if not bspline_path.is_file():
        msg = f"Missing B-spline transform file: {bspline_path}"
        raise FileNotFoundError(msg)

    owned = temp_dir is None
    if temp_dir is None:
        temp_dir = Path(tempfile.mkdtemp(prefix="lightsuite_transformix_"))
    else:
        temp_dir.mkdir(parents=True, exist_ok=True)

    try:
        volumereg = run_transformix(
            moving_volume=vol,
            transform_path=bspline_path,
            output_dir=temp_dir / "bspline",
            spacing_mm=spacing_mm,
            nearest=nearest,
        )
    finally:
        if owned:
            shutil.rmtree(temp_dir, ignore_errors=True)

    volumereg = np.abs(volumereg)
    affine = np.asarray(transform_params.tform_affine_samp20um_to_atlas_10um_px, dtype=float)
    if not warp_pad.is_zero:
        # volumereg is indexed on the padded working grid, the affine expects sample indices.
        affine = affine_with_source_offset(affine, tuple(-b for b in warp_pad.pad_before))
    atlas_shape = tuple(int(v) for v in transform_params.atlassize)
    interp_order = 0 if nearest else 1
    registered = warp_volume_affine(volumereg, affine, atlas_shape, order=interp_order)
    if nearest:
        return (registered > 0).astype(np.uint8)
    return np.clip(registered, 0, np.iinfo(np.uint16).max).astype(np.uint16)


def transform_mask_to_atlas(
    mask: np.ndarray,
    transform_params: TransformParamsCheckpoint,
    *,
    permute: list[int],
    spacing_mm: float,
    temp_dir: Path | None = None,
) -> np.ndarray:
    """Warp a registration-resolution binary mask into atlas space (nearest-neighbor)."""
    return transform_volume_to_atlas(
        (np.asarray(mask) > 0).astype(np.float32),
        transform_params,
        permute=permute,
        spacing_mm=spacing_mm,
        temp_dir=temp_dir,
        nearest=True,
    )
