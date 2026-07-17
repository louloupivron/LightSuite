"""Elastix registration for manifest-driven overview / ROI pairs."""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

import SimpleITK as sitk

from lightsuite.multires.geometry import embed_crop_in_full_overview, resample_to_reference_grid, voxel_count_gb
from lightsuite.multires.prepare import MultiresPreparedPair
from lightsuite.multires.volume import sitk_to_itk, write_sitk_hyperstack_tiff


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
) -> MultiresRegistrationResult:
    """Run itk-elastix on a prepared overview / ROI pair."""
    import itk

    output_dir = output_dir.expanduser()
    output_dir.mkdir(parents=True, exist_ok=True)

    fixed = prepared.overview
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
    if registration_bin > 1:
        result_sitk = resample_to_reference_grid(result_sitk, fixed_cropped)

    cropped_overview_path = output_dir / f"{experiment_slug}_{overview_stem}_cropped_overlap.tif"
    registered_roi_path = (
        output_dir / f"{experiment_slug}_{roi_stem}_registered_to_{overview_stem}.tif"
    )
    sitk.WriteImage(fixed_cropped, str(cropped_overview_path), useCompression=True)
    sitk.WriteImage(result_sitk, str(registered_roi_path), useCompression=True)

    registered_roi_full_overview_path: Path | None = None
    if write_full_overview_canvas:
        full_canvas = embed_crop_in_full_overview(fixed, result_sitk, crop_start_index)
        registered_roi_full_overview_path = (
            output_dir
            / f"{experiment_slug}_{roi_stem}_registered_to_{overview_stem}_in_full_overview.tif"
        )
        sitk.WriteImage(full_canvas, str(registered_roi_full_overview_path), useCompression=True)

    transform_paths = sorted(output_dir.glob("TransformParameters.*.txt"))
    if not transform_paths:
        msg = f"No elastix transform files written under {output_dir}"
        raise RuntimeError(msg)

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
    )


__all__ = [
    "MultiresRegistrationResult",
    "build_elastix_parameter_object",
    "register_roi_to_overview",
    "sanitize_experiment_name",
    "write_sitk_hyperstack_tiff",
]
