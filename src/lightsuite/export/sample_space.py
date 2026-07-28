"""Apply registration transforms to atlas volumes (sample-space export)."""

from __future__ import annotations

import shutil
import tempfile
from pathlib import Path

import numpy as np

from lightsuite.registration.brain_register import TransformParamsCheckpoint
from lightsuite.registration.canvas import (
    RegistrationCanvas,
    WarpCanvasPadding,
    crop_from_warp_canvas,
    undo_canvas_sample_crop,
)
from lightsuite.registration.coordinates import affine_with_source_offset
from lightsuite.registration.elastix.runner import run_transformix
from lightsuite.registration.volume import unpermute_brain_volume
from lightsuite.registration.warp import warp_volume_affine


def resolve_forward_transform_paths(
    transform_params: TransformParamsCheckpoint,
    save_path: Path,
) -> tuple[Path, np.ndarray]:
    """Return B-spline (atlas→sample) path and atlas→sample affine matrix."""
    bspline_path: Path | None = None
    if transform_params.tform_bspline_atlas_20um_to_samp20um_px:
        candidate = Path(transform_params.tform_bspline_atlas_20um_to_samp20um_px)
        if candidate.is_file():
            bspline_path = candidate
    if bspline_path is None:
        fallback = save_path / "bspline_atlas_to_samp_20um.txt"
        if fallback.is_file():
            bspline_path = fallback
    if bspline_path is None:
        msg = (
            "Missing forward B-spline transform (bspline_atlas_to_samp_20um.txt). "
            "Re-run 'lightsuite brain register'."
        )
        raise FileNotFoundError(msg)

    if transform_params.tform_affine_atlas_to_samp20um_px is not None:
        affine = np.asarray(transform_params.tform_affine_atlas_to_samp20um_px, dtype=float)
    else:
        affine_inv = np.asarray(
            transform_params.tform_affine_samp20um_to_atlas_10um_px,
            dtype=float,
        )
        affine = np.linalg.inv(affine_inv)
        if transform_params.atlas_crop_start_native is not None:
            crop = tuple(int(v) for v in transform_params.atlas_crop_start_native)
            if any(crop):
                affine = affine_with_source_offset(affine, crop)

    return bspline_path, affine


def _registration_canvas(transform_params: TransformParamsCheckpoint) -> RegistrationCanvas | None:
    canvas = RegistrationCanvas.from_checkpoint_dict(transform_params.registration_canvas)
    if canvas is not None:
        return canvas
    if transform_params.warp_canvas_pad_before is not None:
        pad = tuple(int(v) for v in transform_params.warp_canvas_pad_before)
        if any(pad):
            regvol = tuple(int(v) for v in transform_params.regvolsize)
            working = tuple(s + pad[i] for i, s in enumerate(regvol))
            return RegistrationCanvas(
                mode="legacy_vd_pad",
                pad_before=pad,
                pad_after=(0, 0, 0),
                sample_crop_start=(0, 0, 0),
                working_shape=working,
            )
    regvol = tuple(int(v) for v in transform_params.regvolsize)
    return RegistrationCanvas(
        mode="off",
        pad_before=(0, 0, 0),
        pad_after=(0, 0, 0),
        sample_crop_start=(0, 0, 0),
        working_shape=regvol,
    )


def transform_atlas_volume_to_sample(
    atlas_volume: np.ndarray,
    transform_params: TransformParamsCheckpoint,
    *,
    save_path: Path,
    spacing_mm: float,
    temp_dir: Path | None = None,
    nearest: bool = True,
) -> np.ndarray:
    """Warp an atlas-grid volume onto the permuted registration sample grid."""
    canvas = _registration_canvas(transform_params)
    bspline_path, affine = resolve_forward_transform_paths(transform_params, save_path)
    warp_pad = WarpCanvasPadding(canvas.pad_before, canvas.pad_after)
    working_shape = warp_pad.padded_shape(canvas.working_shape)

    vol = atlas_volume.astype(np.float32)
    warped_affine = warp_volume_affine(
        vol,
        affine,
        working_shape,
        order=0 if nearest else 1,
    )

    owned = temp_dir is None
    if temp_dir is None:
        temp_dir = Path(tempfile.mkdtemp(prefix="lightsuite_transformix_sample_"))
    else:
        temp_dir.mkdir(parents=True, exist_ok=True)

    try:
        warped_bspline = run_transformix(
            moving_volume=warped_affine,
            transform_path=bspline_path,
            output_dir=temp_dir / "bspline",
            spacing_mm=spacing_mm,
            nearest=nearest,
        )
    finally:
        if owned:
            shutil.rmtree(temp_dir, ignore_errors=True)

    inner_shape = tuple(
        canvas.working_shape[i] - canvas.pad_before[i] - canvas.pad_after[i]
        for i in range(3)
    )
    cropped = crop_from_warp_canvas(warped_bspline, warp_pad, inner_shape)
    full_shape = tuple(int(v) for v in transform_params.regvolsize)
    embedded = undo_canvas_sample_crop(cropped, canvas, full_shape)
    return unpermute_brain_volume(embedded, transform_params.permute_sample_to_atlas)
