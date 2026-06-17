"""Orchestration for mesoSPIM overview / ROI registration stages."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from lightsuite.mesospim.checkpoint import MesospimRegOptsCheckpoint, mesospim_checkpoint_path
from lightsuite.mesospim.config_models import MesospimGeometryMode, MesospimPipelineConfig
from lightsuite.mesospim.geometry import (
    geometry_report,
    overlap_physical_bounds,
    transformed_bounds_in_target_space,
    voxel_geometry_report,
)
from lightsuite.mesospim.io import tiff_shape, write_sitk_hyperstack_tiff
from lightsuite.mesospim.meta import parse_mesospim_meta
from lightsuite.mesospim.plots import save_fov_overlap_plot, save_geometry_overlap_qc_plot
from lightsuite.mesospim.prepare import (
    prepare_mesospim_registration_pair,
    resolve_meta_paths,
)
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


def _serialize_report(report: dict[str, object]) -> dict[str, object]:
    out: dict[str, object] = {}
    for key, value in report.items():
        if isinstance(value, np.ndarray):
            out[key] = value.tolist()
        else:
            out[key] = value
    return out


def _landmark_roi_report(
    roi,
    fit,
) -> dict[str, object]:
    roi_min, roi_max = transformed_bounds_in_target_space(roi, fit.roi_to_overview_tform)
    center = 0.5 * (roi_min + roi_max)
    return {
        "label": "roi_in_overview_space",
        "phys_min": roi_min,
        "phys_max": roi_max,
        "phys_center": center,
        "landmark_rms_error_um": fit.rms_error_um,
        "landmark_fit_stats": fit.fit_stats,
    }


def check_mesospim_geometry(cfg: MesospimPipelineConfig) -> MesospimRegOptsCheckpoint:
    """Validate FOV overlap and write geometry QA artifacts."""
    meso = cfg.mesospim
    mode = meso.geometry_mode
    overview_path = meso.overview.path
    roi_path = meso.roi.path
    overview_meta_path, roi_meta_path = resolve_meta_paths(cfg)

    if mode == MesospimGeometryMode.METADATA:
        return _check_metadata_geometry(cfg)

    prepared = prepare_mesospim_registration_pair(cfg)
    overlap_min, overlap_max = prepared.overlap_box
    shape_overview = tiff_shape(overview_path)

    if mode == MesospimGeometryMode.LANDMARKS:
        assert meso.overview.voxel_um is not None and meso.roi.voxel_um is not None
        rep_overview = voxel_geometry_report("overview", shape_overview, meso.overview.voxel_um)
        rep_roi = _landmark_roi_report(prepared.roi, prepared.landmark_fit)
    else:
        meta_overview = parse_mesospim_meta(meso.overview.resolved_meta_path())
        rep_overview = geometry_report("overview", meta_overview, shape_overview, meso.geometry)
        rep_roi = _landmark_roi_report(prepared.roi, prepared.landmark_fit)

    geometry_dir = cfg.sample.save_path / "geometry"
    geometry_dir.mkdir(parents=True, exist_ok=True)

    report_path = geometry_dir / "geometry_report.json"
    serializable = {
        "geometry_mode": mode.value,
        "overview": _serialize_report(rep_overview),
        "roi": _serialize_report(rep_roi),
        "overlap_box_um": [overlap_min.tolist(), overlap_max.tolist()],
    }
    if prepared.landmark_fit is not None:
        serializable["landmark_fit"] = {
            "rms_error_um": prepared.landmark_fit.rms_error_um,
            "fit_stats": prepared.landmark_fit.fit_stats,
            "roi_to_overview_tform": prepared.landmark_fit.roi_to_overview_tform.tolist(),
        }
    report_path.write_text(json.dumps(serializable, indent=2), encoding="utf-8")

    overview_stem = _volume_stem(overview_path)
    cropped_path = geometry_dir / f"{overview_stem}_cropped_overlap.tif"

    title = f"{overview_stem} vs {_volume_stem(roi_path)} — landmark placement (µm)"
    fov_plot_path = geometry_dir / "fov_overlap.png"
    save_fov_overlap_plot(
        rep_overview=rep_overview,
        rep_roi=rep_roi,
        overlap_min=overlap_min,
        overlap_max=overlap_max,
        output_path=fov_plot_path,
        title=title,
    )

    qc_plot_path = geometry_dir / "geometry_overlap_qc.png"
    meta_overlap_min, meta_overlap_max = overlap_physical_bounds(
        prepared.overview,
        prepared.roi,
        margin_um=meso.registration.overlap_margin_um,
    )
    meta_overview_dict = parse_mesospim_meta(meso.overview.resolved_meta_path())
    meta_roi_dict = parse_mesospim_meta(meso.roi.resolved_meta_path())
    save_geometry_overlap_qc_plot(
        overview_path=overview_path,
        roi_path=roi_path,
        overview_meta=meta_overview_dict,
        roi_meta=meta_roi_dict,
        geometry=meso.geometry,
        tiff_remap=meso.tiff_remap,
        metadata_overlap_min=meta_overlap_min,
        metadata_overlap_max=meta_overlap_max,
        hybrid_overlap_min=overlap_min,
        hybrid_overlap_max=overlap_max,
        roi_to_overview=prepared.landmark_fit.roi_to_overview_tform,
        output_path=qc_plot_path,
        geometry_mode=mode.value,
    )

    write_sitk_hyperstack_tiff(cropped_path, prepared.fixed_cropped)

    experiment_slug = sanitize_experiment_name(meso.registration.experiment_name)
    landmark_session_path = str(meso.resolved_landmark_session_path(cfg.sample.save_path))
    checkpoint = MesospimRegOptsCheckpoint(
        sample_name=cfg.sample.name,
        overview_path=str(overview_path),
        roi_path=str(roi_path),
        overview_meta_path=overview_meta_path or "",
        roi_meta_path=roi_meta_path or "",
        experiment_slug=experiment_slug,
        geometry_mode=mode.value,
        landmark_session_path=landmark_session_path,
        overlap_box_um=[overlap_min.tolist(), overlap_max.tolist()],
        crop_start_index=prepared.crop_start_index,
        roi_to_overview_tform=prepared.landmark_fit.roi_to_overview_tform.tolist()
        if prepared.landmark_fit
        else None,
        landmark_rms_error_um=prepared.landmark_fit.rms_error_um if prepared.landmark_fit else None,
        geometry_report_paths={
            "geometry_report_json": str(report_path),
            "fov_overlap_png": str(fov_plot_path),
            "geometry_overlap_qc_png": str(qc_plot_path),
            "cropped_overview_preview": str(cropped_path),
        },
    )
    checkpoint.save(mesospim_checkpoint_path(cfg.sample.save_path))
    _status(f"Geometry QA ({mode.value}): {geometry_dir}")
    _status(f"Overlap QC plot: {qc_plot_path}")
    if prepared.landmark_fit is not None:
        _status(f"Landmark RMS error: {prepared.landmark_fit.rms_error_um:.2f} µm")
    return checkpoint


def _check_metadata_geometry(cfg: MesospimPipelineConfig) -> MesospimRegOptsCheckpoint:
    meso = cfg.mesospim
    overview_path = meso.overview.path
    roi_path = meso.roi.path
    overview_meta_path = meso.overview.resolved_meta_path()
    roi_meta_path = meso.roi.resolved_meta_path()

    meta_overview = parse_mesospim_meta(overview_meta_path)
    meta_roi = parse_mesospim_meta(roi_meta_path)

    prepared = prepare_mesospim_registration_pair(cfg)
    overlap_min, overlap_max = prepared.overlap_box

    shape_overview = tiff_shape(overview_path)
    shape_roi = tiff_shape(roi_path)
    rep_overview = geometry_report("overview", meta_overview, shape_overview, meso.geometry)
    rep_roi = geometry_report("roi", meta_roi, shape_roi, meso.geometry)

    geometry_dir = cfg.sample.save_path / "geometry"
    geometry_dir.mkdir(parents=True, exist_ok=True)

    report_path = geometry_dir / "geometry_report.json"
    serializable = {
        "geometry_mode": MesospimGeometryMode.METADATA.value,
        "overview": _serialize_report(rep_overview),
        "roi": _serialize_report(rep_roi),
        "overlap_box_um": [overlap_min.tolist(), overlap_max.tolist()],
    }
    report_path.write_text(json.dumps(serializable, indent=2), encoding="utf-8")

    overview_stem = _volume_stem(overview_path)
    cropped_path = geometry_dir / f"{overview_stem}_cropped_overlap.tif"

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

    qc_plot_path = geometry_dir / "geometry_overlap_qc.png"
    save_geometry_overlap_qc_plot(
        overview_path=overview_path,
        roi_path=roi_path,
        overview_meta=meta_overview,
        roi_meta=meta_roi,
        geometry=meso.geometry,
        tiff_remap=meso.tiff_remap,
        metadata_overlap_min=overlap_min,
        metadata_overlap_max=overlap_max,
        output_path=qc_plot_path,
        geometry_mode=MesospimGeometryMode.METADATA.value,
    )

    write_sitk_hyperstack_tiff(cropped_path, prepared.fixed_cropped)

    experiment_slug = sanitize_experiment_name(meso.registration.experiment_name)
    checkpoint = MesospimRegOptsCheckpoint(
        sample_name=cfg.sample.name,
        overview_path=str(overview_path),
        roi_path=str(roi_path),
        overview_meta_path=str(overview_meta_path),
        roi_meta_path=str(roi_meta_path),
        experiment_slug=experiment_slug,
        geometry_mode=MesospimGeometryMode.METADATA.value,
        overlap_box_um=[overlap_min.tolist(), overlap_max.tolist()],
        crop_start_index=prepared.crop_start_index,
        geometry_report_paths={
            "geometry_report_json": str(report_path),
            "fov_overlap_png": str(fov_plot_path),
            "geometry_overlap_qc_png": str(qc_plot_path),
            "cropped_overview_preview": str(cropped_path),
        },
    )
    checkpoint.save(mesospim_checkpoint_path(cfg.sample.save_path))
    _status(f"Geometry QA: {geometry_dir}")
    _status(f"Overlap QC plot: {qc_plot_path}")
    return checkpoint


def run_mesospim_registration(cfg: MesospimPipelineConfig) -> MesospimRegOptsCheckpoint:
    """Register ROI stack to overview using configured geometry and elastix."""
    meso = cfg.mesospim
    overview_path = meso.overview.path
    roi_path = meso.roi.path
    overview_meta_path, roi_meta_path = resolve_meta_paths(cfg)

    experiment_slug = sanitize_experiment_name(meso.registration.experiment_name)
    output_dir = cfg.sample.save_path / "elastix_roi_to_overview" / experiment_slug

    prepared = prepare_mesospim_registration_pair(cfg)
    result = register_roi_to_overview(
        cfg=cfg,
        prepared=prepared,
        output_dir=output_dir,
        experiment_slug=experiment_slug,
        overview_stem=_volume_stem(overview_path),
        roi_stem=_volume_stem(roi_path),
        registration_bin=meso.registration.registration_bin,
        elastix_stages=meso.registration.elastix_stages,
        write_full_overview_canvas=meso.registration.write_full_overview_canvas,
    )

    landmark_session_path = None
    if meso.geometry_mode != MesospimGeometryMode.METADATA:
        landmark_session_path = str(meso.resolved_landmark_session_path(cfg.sample.save_path))

    checkpoint = MesospimRegOptsCheckpoint(
        sample_name=cfg.sample.name,
        overview_path=str(overview_path),
        roi_path=str(roi_path),
        overview_meta_path=overview_meta_path or "",
        roi_meta_path=roi_meta_path or "",
        experiment_slug=experiment_slug,
        geometry_mode=meso.geometry_mode.value,
        landmark_session_path=landmark_session_path,
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
        roi_to_overview_tform=result.roi_to_overview_tform,
        landmark_rms_error_um=(
            prepared.landmark_fit.rms_error_um if prepared.landmark_fit is not None else None
        ),
    )
    checkpoint.save(mesospim_checkpoint_path(cfg.sample.save_path))
    if result.registered_roi_full_overview_path is not None:
        _status(f"Full overview canvas: {result.registered_roi_full_overview_path}")
    _status(f"Registered ROI: {result.registered_roi_path}")
    return checkpoint
