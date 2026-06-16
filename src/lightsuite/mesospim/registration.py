"""Elastix registration for mesoSPIM overview / ROI pairs."""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

import SimpleITK as sitk

from lightsuite.mesospim.geometry import prepare_registration_pair, voxel_count_gb
from lightsuite.mesospim.io import read_tiff_as_float, sitk_to_itk


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
class MesospimRegistrationResult:
    output_dir: Path
    transform_paths: list[Path]
    cropped_overview_path: Path
    registered_roi_path: Path
    overlap_box_um: tuple[list[float], list[float]]
    fixed_voxel_gb: float
    moving_voxel_gb: float


def register_roi_to_overview(
    *,
    overview_path: Path,
    roi_path: Path,
    overview_meta: dict,
    roi_meta: dict,
    geometry,
    tiff_remap,
    output_dir: Path,
    experiment_slug: str,
    overview_stem: str,
    roi_stem: str,
    overlap_margin_um: float,
    registration_bin: int,
    elastix_stages: list[str],
) -> MesospimRegistrationResult:
    """Load stacks, prepare pair, run itk-elastix, and write outputs."""
    import itk

    overview_path = overview_path.expanduser().resolve()
    roi_path = roi_path.expanduser().resolve()
    output_dir = output_dir.expanduser()
    output_dir.mkdir(parents=True, exist_ok=True)

    fixed = read_tiff_as_float(
        overview_path,
        overview_path=overview_path,
        roi_path=roi_path,
        remap=tiff_remap,
    )
    roi_full = read_tiff_as_float(
        roi_path,
        overview_path=overview_path,
        roi_path=roi_path,
        remap=tiff_remap,
    )

    from lightsuite.mesospim.geometry import apply_image_geometry

    apply_image_geometry(fixed, overview_meta, geometry)
    apply_image_geometry(roi_full, roi_meta, geometry)

    fixed, moving, overlap_box = prepare_registration_pair(
        fixed,
        roi_full,
        margin_um=overlap_margin_um,
    )

    if registration_bin > 1:
        shrink = (registration_bin, registration_bin, registration_bin)
        fixed = sitk.Shrink(fixed, shrink)
        moving = sitk.Shrink(moving, shrink)

    fixed_gb = voxel_count_gb(fixed)
    moving_gb = voxel_count_gb(moving)

    parameter_object = build_elastix_parameter_object(elastix_stages)
    result_itk, _ = itk.elastix_registration_method(
        sitk_to_itk(fixed),
        sitk_to_itk(moving),
        parameter_object=parameter_object,
        output_directory=str(output_dir),
        log_to_console=False,
    )

    result_sitk = sitk.GetImageFromArray(itk.GetArrayFromImage(result_itk))
    result_sitk.CopyInformation(fixed)

    cropped_overview_path = output_dir / f"{experiment_slug}_{overview_stem}_cropped_overlap.tif"
    registered_roi_path = (
        output_dir / f"{experiment_slug}_{roi_stem}_registered_to_{overview_stem}.tif"
    )
    sitk.WriteImage(fixed, str(cropped_overview_path), useCompression=True)
    sitk.WriteImage(result_sitk, str(registered_roi_path), useCompression=True)

    transform_paths = sorted(output_dir.glob("TransformParameters.*.txt"))
    if not transform_paths:
        msg = f"No elastix transform files written under {output_dir}"
        raise RuntimeError(msg)

    overlap_min, overlap_max = overlap_box
    return MesospimRegistrationResult(
        output_dir=output_dir,
        transform_paths=transform_paths,
        cropped_overview_path=cropped_overview_path,
        registered_roi_path=registered_roi_path,
        overlap_box_um=(overlap_min.tolist(), overlap_max.tolist()),
        fixed_voxel_gb=fixed_gb,
        moving_voxel_gb=moving_gb,
    )
