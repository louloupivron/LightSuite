"""Orchestration for mesoSPIM overview / ROI registration stages."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from lightsuite.mesospim.checkpoint import MesospimRegOptsCheckpoint, mesospim_checkpoint_path
from lightsuite.mesospim.config_models import MesospimPipelineConfig
from lightsuite.mesospim.geometry import (
    apply_image_geometry,
    crop_to_physical_box,
    geometry_report,
    overlap_physical_bounds,
)
from lightsuite.mesospim.io import (
    empty_image_from_shape,
    read_tiff_as_float,
    tiff_shape,
    write_sitk_hyperstack_tiff,
)
from lightsuite.mesospim.meta import parse_mesospim_meta
from lightsuite.mesospim.plots import save_fov_overlap_plot, save_overlap_slice_plot
from lightsuite.mesospim.registration import register_roi_to_overview, sanitize_experiment_name


def _status(message: str) -> None:
    print(message, flush=True)


def _volume_stem(path: Path) -> str:
    name = path.name
    if name.lower().endswith(".tif"):
        return name[:-4]
    if name.lower().endswith(".tiff"):
        return name[:-5]
    return path.stem


def check_mesospim_geometry(cfg: MesospimPipelineConfig) -> MesospimRegOptsCheckpoint:
    """Validate FOV overlap and write geometry QA artifacts."""
    meso = cfg.mesospim
    overview_path = meso.overview.path
    roi_path = meso.roi.path
    overview_meta_path = meso.overview.resolved_meta_path()
    roi_meta_path = meso.roi.resolved_meta_path()

    meta_overview = parse_mesospim_meta(overview_meta_path)
    meta_roi = parse_mesospim_meta(roi_meta_path)

    shape_overview = tiff_shape(overview_path)
    shape_roi = tiff_shape(roi_path)
    rep_overview = geometry_report("overview", meta_overview, shape_overview, meso.geometry)
    rep_roi = geometry_report("roi", meta_roi, shape_roi, meso.geometry)

    img_overview = empty_image_from_shape(shape_overview)
    img_roi = empty_image_from_shape(shape_roi)
    apply_image_geometry(img_overview, meta_overview, meso.geometry)
    apply_image_geometry(img_roi, meta_roi, meso.geometry)
    overlap_min, overlap_max = overlap_physical_bounds(
        img_overview,
        img_roi,
        margin_um=meso.registration.overlap_margin_um,
    )

    geometry_dir = cfg.sample.save_path / "geometry"
    geometry_dir.mkdir(parents=True, exist_ok=True)

    report_path = geometry_dir / "geometry_report.json"
    serializable = {
        "overview": _serialize_report(rep_overview),
        "roi": _serialize_report(rep_roi),
        "overlap_box_um": [overlap_min.tolist(), overlap_max.tolist()],
    }
    report_path.write_text(json.dumps(serializable, indent=2), encoding="utf-8")

    overview_stem = _volume_stem(overview_path)
    cropped_path = geometry_dir / f"{overview_stem}_cropped_overlap.tif"
    overview_pixels = read_tiff_as_float(
        overview_path,
        overview_path=overview_path,
        roi_path=roi_path,
        remap=meso.tiff_remap,
    )
    apply_image_geometry(overview_pixels, meta_overview, meso.geometry)
    cropped_pixels, _crop_start = crop_to_physical_box(overview_pixels, overlap_min, overlap_max)
    write_sitk_hyperstack_tiff(cropped_path, cropped_pixels)

    title = f"{overview_stem} vs {_volume_stem(roi_path)} — stage frame (µm)"
    fov_plot_path = geometry_dir / "fov_overlap.png"
    save_fov_overlap_plot(
        rep_overview=rep_overview,
        rep_roi=rep_roi,
        overlap_min=overlap_min,
        overlap_max=overlap_max,
        output_path=fov_plot_path,
        title=title,
    )

    slice_plot_path = geometry_dir / "xy_slice_overlap.png"
    save_overlap_slice_plot(
        overview_path=overview_path,
        roi_path=roi_path,
        overview_meta=meta_overview,
        roi_meta=meta_roi,
        geometry=meso.geometry,
        tiff_remap=meso.tiff_remap,
        overlap_min=overlap_min,
        overlap_max=overlap_max,
        output_path=slice_plot_path,
    )

    experiment_slug = sanitize_experiment_name(meso.registration.experiment_name)
    checkpoint = MesospimRegOptsCheckpoint(
        sample_name=cfg.sample.name,
        overview_path=str(overview_path),
        roi_path=str(roi_path),
        overview_meta_path=str(overview_meta_path),
        roi_meta_path=str(roi_meta_path),
        experiment_slug=experiment_slug,
        overlap_box_um=[overlap_min.tolist(), overlap_max.tolist()],
        geometry_report_paths={
            "geometry_report_json": str(report_path),
            "fov_overlap_png": str(fov_plot_path),
            "xy_slice_overlap_png": str(slice_plot_path),
            "cropped_overview_preview": str(cropped_path),
        },
    )
    checkpoint.save(mesospim_checkpoint_path(cfg.sample.save_path))
    _status(f"Geometry QA: {geometry_dir}")
    return checkpoint


def run_mesospim_registration(cfg: MesospimPipelineConfig) -> MesospimRegOptsCheckpoint:
    """Register ROI stack to overview using metadata geometry and elastix."""
    meso = cfg.mesospim
    overview_path = meso.overview.path
    roi_path = meso.roi.path
    overview_meta_path = meso.overview.resolved_meta_path()
    roi_meta_path = meso.roi.resolved_meta_path()

    meta_overview = parse_mesospim_meta(overview_meta_path)
    meta_roi = parse_mesospim_meta(roi_meta_path)

    experiment_slug = sanitize_experiment_name(meso.registration.experiment_name)
    output_dir = cfg.sample.save_path / "elastix_roi_to_overview" / experiment_slug

    result = register_roi_to_overview(
        overview_path=overview_path,
        roi_path=roi_path,
        overview_meta=meta_overview,
        roi_meta=meta_roi,
        geometry=meso.geometry,
        tiff_remap=meso.tiff_remap,
        output_dir=output_dir,
        experiment_slug=experiment_slug,
        overview_stem=_volume_stem(overview_path),
        roi_stem=_volume_stem(roi_path),
        overlap_margin_um=meso.registration.overlap_margin_um,
        registration_bin=meso.registration.registration_bin,
        elastix_stages=meso.registration.elastix_stages,
        write_full_overview_canvas=meso.registration.write_full_overview_canvas,
    )

    checkpoint = MesospimRegOptsCheckpoint(
        sample_name=cfg.sample.name,
        overview_path=str(overview_path),
        roi_path=str(roi_path),
        overview_meta_path=str(overview_meta_path),
        roi_meta_path=str(roi_meta_path),
        experiment_slug=experiment_slug,
        overlap_box_um=list(result.overlap_box_um),
        elastix_output_dir=str(result.output_dir),
        transform_paths=[str(p) for p in result.transform_paths],
        cropped_overview_path=str(result.cropped_overview_path),
        registered_roi_path=str(result.registered_roi_path),
        registered_roi_full_overview_path=(
            str(result.registered_roi_full_overview_path)
            if result.registered_roi_full_overview_path is not None
            else None
        ),
        crop_start_index=result.crop_start_index,
    )
    checkpoint.save(mesospim_checkpoint_path(cfg.sample.save_path))
    if result.registered_roi_full_overview_path is not None:
        _status(f"Full overview canvas: {result.registered_roi_full_overview_path}")
    _status(f"Registered ROI: {result.registered_roi_path}")
    return checkpoint


def _serialize_report(report: dict[str, object]) -> dict[str, object]:
    out: dict[str, object] = {}
    for key, value in report.items():
        if isinstance(value, np.ndarray):
            out[key] = value.tolist()
        else:
            out[key] = value
    return out
