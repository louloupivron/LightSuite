"""Elastix registration for manifest-driven overview / ROI pairs."""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

import SimpleITK as sitk

from lightsuite.multires.geometry import voxel_count_gb
from lightsuite.multires.prepare import MultiresPreparedPair
from lightsuite.multires.volume import (
    sitk_to_itk,
    write_embedded_crop_canvas,
    write_sitk_hyperstack_tiff,
)


def sanitize_experiment_name(name: str) -> str:
    slug = re.sub(r"[^\w\-]+", "_", name.strip())
    if not slug:
        msg = "experiment_name must not be empty after sanitization"
        raise ValueError(msg)
    return slug


def build_elastix_parameter_object(stages: tuple[str, ...] | list[str]):
    """Build itk-elastix parameter object for the requested stages."""
    import itk

    parameter_object = itk.ParameterObject.New()
    stage_list = list(stages)
    for stage in stage_list:
        pmap = itk.ParameterObject.GetDefaultParameterMap(stage, 3)
        pmap["UseDirectionCosines"] = ["true"]
        pmap["WriteResultImage"] = ["false"]
        pmap["WriteIterationInfo"] = ["false"]
        parameter_object.AddParameterMap(pmap)
    last = parameter_object.GetNumberOfParameterMaps() - 1
    parameter_object.SetParameter(last, "WriteResultImage", "true")
    return parameter_object


@dataclass
class MultiresChannelRegistrationResult:
    channel: str
    registered_roi_path: Path
    registered_roi_full_overview_path: Path | None = None
    registration_overlay_qc_path: Path | None = None
    registration_slice_ncc: float | None = None


@dataclass
class MultiresRegistrationResult:
    output_dir: Path
    transform_paths: list[Path]
    cropped_overview_path: Path
    registered_roi_path: Path
    registered_roi_full_overview_path: Path | None
    overlap_box_um: tuple[list[float], list[float]]
    crop_start_index: list[int]
    fixed_voxel_gb: float
    moving_voxel_gb: float
    roi_to_overview_tform: list[list[float]] | None = None
    registration_overlay_qc_path: Path | None = None
    registration_slice_ncc: float | None = None
    reference_channel: str | None = None
    additional_channels: list[MultiresChannelRegistrationResult] = field(default_factory=list)


def apply_elastix_transforms(
    moving: sitk.Image,
    transform_paths: list[Path],
    *,
    reference: sitk.Image,
) -> sitk.Image:
    """Apply saved elastix transform parameter files to *moving* on *reference* grid."""
    import itk
    import numpy as np

    parameter_object = itk.ParameterObject.New()
    for path in transform_paths:
        parameter_object.AddParameterFile(str(path))

    # Point the final transform's output grid at the full-resolution reference.
    # Transforms estimated at registration_bin>1 store the binned Size/Spacing;
    # without this update transformix writes a coarse grid whose geometry is then
    # lost by array round-trips, producing an all-zero resample.
    last = parameter_object.GetNumberOfParameterMaps() - 1
    parameter_object.SetParameter(last, "Size", [str(int(v)) for v in reference.GetSize()])
    parameter_object.SetParameter(last, "Spacing", [str(float(v)) for v in reference.GetSpacing()])
    parameter_object.SetParameter(last, "Origin", [str(float(v)) for v in reference.GetOrigin()])
    parameter_object.SetParameter(
        last,
        "Direction",
        [str(float(v)) for v in reference.GetDirection()],
    )

    result_itk = itk.transformix_filter(
        sitk_to_itk(moving),
        transform_parameter_object=parameter_object,
    )
    result_sitk = sitk.GetImageFromArray(itk.GetArrayFromImage(result_itk).astype(np.float32))
    result_sitk.CopyInformation(reference)
    return result_sitk


def register_roi_to_overview(
    *,
    prepared: MultiresPreparedPair,
    output_dir: Path,
    experiment_slug: str,
    overview_stem: str,
    roi_stem: str,
    registration_bin: int,
    elastix_stages: list[str],
    write_full_overview_canvas: bool = True,
    pair_label: str | None = None,
    channel: str | None = None,
) -> MultiresRegistrationResult:
    """Run itk-elastix on a prepared overview / ROI pair."""
    import itk

    output_dir = output_dir.expanduser()
    output_dir.mkdir(parents=True, exist_ok=True)

    fixed_cropped = prepared.fixed_cropped
    moving = prepared.moving
    overlap_box = prepared.overlap_box
    crop_start_index = prepared.crop_start_index

    fixed_for_elastix = fixed_cropped
    moving_for_elastix = moving
    if registration_bin > 1:
        shrink = (registration_bin, registration_bin, registration_bin)
        fixed_for_elastix = sitk.Shrink(fixed_cropped, shrink)
        moving_for_elastix = sitk.Shrink(moving, shrink)

    fixed_gb = voxel_count_gb(fixed_for_elastix)
    moving_gb = voxel_count_gb(moving_for_elastix)

    parameter_object = build_elastix_parameter_object(elastix_stages)
    result_itk, _ = itk.elastix_registration_method(
        sitk_to_itk(fixed_for_elastix),
        sitk_to_itk(moving_for_elastix),
        parameter_object=parameter_object,
        output_directory=str(output_dir),
        log_to_console=False,
    )

    result_sitk = sitk.GetImageFromArray(itk.GetArrayFromImage(result_itk))
    result_sitk.CopyInformation(fixed_for_elastix)
    transform_paths = sorted(output_dir.glob("TransformParameters.*.txt"))
    if not transform_paths:
        msg = f"No elastix transform files written under {output_dir}"
        raise RuntimeError(msg)
    if registration_bin > 1:
        # Prefer re-applying transforms on the full-res moving image; if that
        # fails, upsample the binned elastix result in physical space.
        result_sitk = apply_elastix_transforms(
            moving,
            transform_paths,
            reference=fixed_cropped,
        )

    channel_prefix = f"{channel}_" if channel else ""
    cropped_overview_path = (
        output_dir / f"{experiment_slug}_{channel_prefix}{overview_stem}_cropped_overlap.tif"
    )
    registered_roi_path = (
        output_dir
        / f"{experiment_slug}_{channel_prefix}{roi_stem}_registered_to_{overview_stem}.tif"
    )
    sitk.WriteImage(fixed_cropped, str(cropped_overview_path), useCompression=True)
    sitk.WriteImage(result_sitk, str(registered_roi_path), useCompression=True)

    registration_overlay_qc_path = output_dir / (
        f"registration_overlay_qc_{channel}.png" if channel else "registration_overlay_qc.png"
    )
    registration_slice_ncc: float | None = None
    try:
        from lightsuite.multires.plots import save_registration_overlay_qc_plot

        registration_slice_ncc = save_registration_overlay_qc_plot(
            overview_crop=fixed_cropped,
            registered_roi=result_sitk,
            output_path=registration_overlay_qc_path,
            pair_label=pair_label,
        )
    except (ImportError, ValueError) as exc:
        registration_overlay_qc_path = None
        import warnings

        warnings.warn(f"Registration overlay QC plot skipped: {exc}", stacklevel=1)

    registered_roi_full_overview_path: Path | None = None
    if write_full_overview_canvas:
        registered_roi_full_overview_path = (
            output_dir
            / (
                f"{experiment_slug}_{channel_prefix}{roi_stem}_registered_to_"
                f"{overview_stem}_in_full_overview.tif"
            )
        )
        write_embedded_crop_canvas(
            prepared.overview_spec,
            result_sitk,
            crop_start_index,
            registered_roi_full_overview_path,
        )

    overlap_min, overlap_max = overlap_box
    roi_tform = None
    if prepared.landmark_fit is not None:
        roi_tform = prepared.landmark_fit.roi_to_overview_tform.tolist()

    return MultiresRegistrationResult(
        output_dir=output_dir,
        transform_paths=transform_paths,
        cropped_overview_path=cropped_overview_path,
        registered_roi_path=registered_roi_path,
        registered_roi_full_overview_path=registered_roi_full_overview_path,
        overlap_box_um=(overlap_min.tolist(), overlap_max.tolist()),
        crop_start_index=crop_start_index,
        fixed_voxel_gb=fixed_gb,
        moving_voxel_gb=moving_gb,
        roi_to_overview_tform=roi_tform,
        registration_overlay_qc_path=registration_overlay_qc_path,
        registration_slice_ncc=registration_slice_ncc,
        reference_channel=channel,
    )


def apply_registration_to_channel(
    *,
    prepared: MultiresPreparedPair,
    transform_paths: list[Path],
    output_dir: Path,
    experiment_slug: str,
    channel: str,
    overview_stem: str,
    roi_stem: str,
    write_full_overview_canvas: bool = True,
    pair_label: str | None = None,
) -> MultiresChannelRegistrationResult:
    """Apply saved elastix transforms to another channel without re-running elastix."""
    output_dir = output_dir.expanduser()
    output_dir.mkdir(parents=True, exist_ok=True)

    result_sitk = apply_elastix_transforms(
        prepared.moving,
        transform_paths,
        reference=prepared.fixed_cropped,
    )

    registered_roi_path = (
        output_dir
        / f"{experiment_slug}_{channel}_{roi_stem}_registered_to_{overview_stem}.tif"
    )
    sitk.WriteImage(result_sitk, str(registered_roi_path), useCompression=True)

    registration_overlay_qc_path = output_dir / f"registration_overlay_qc_{channel}.png"
    registration_slice_ncc: float | None = None
    try:
        from lightsuite.multires.plots import save_registration_overlay_qc_plot

        registration_slice_ncc = save_registration_overlay_qc_plot(
            overview_crop=prepared.fixed_cropped,
            registered_roi=result_sitk,
            output_path=registration_overlay_qc_path,
            pair_label=pair_label or channel,
        )
    except (ImportError, ValueError) as exc:
        registration_overlay_qc_path = None
        import warnings

        warnings.warn(f"Registration overlay QC plot skipped for {channel}: {exc}", stacklevel=1)

    registered_roi_full_overview_path: Path | None = None
    if write_full_overview_canvas:
        registered_roi_full_overview_path = (
            output_dir
            / (
                f"{experiment_slug}_{channel}_{roi_stem}_registered_to_"
                f"{overview_stem}_in_full_overview.tif"
            )
        )
        write_embedded_crop_canvas(
            prepared.overview_spec,
            result_sitk,
            prepared.crop_start_index,
            registered_roi_full_overview_path,
        )

    return MultiresChannelRegistrationResult(
        channel=channel,
        registered_roi_path=registered_roi_path,
        registered_roi_full_overview_path=registered_roi_full_overview_path,
        registration_overlay_qc_path=registration_overlay_qc_path,
        registration_slice_ncc=registration_slice_ncc,
    )


__all__ = [
    "MultiresChannelRegistrationResult",
    "MultiresRegistrationResult",
    "apply_elastix_transforms",
    "apply_registration_to_channel",
    "build_elastix_parameter_object",
    "register_roi_to_overview",
    "sanitize_experiment_name",
    "write_sitk_hyperstack_tiff",
]
