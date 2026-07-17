"""Orchestration for manifest-driven multiresolution registration."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from lightsuite.multires.checkpoint import MultiresRegOptsCheckpoint, multires_checkpoint_path
from lightsuite.multires.config_models import MultiresGeometryMode, MultiresPipelineConfig
from lightsuite.multires.geometry import physical_bounds, transformed_bounds_in_target_space
from lightsuite.multires.manifest import load_pair_manifest
from lightsuite.multires.models import MultiresPairManifest, serialize_report
from lightsuite.multires.plots import save_fov_overlap_plot, save_geometry_overlap_qc_plot
from lightsuite.multires.prepare import prepare_multires_registration_pair
from lightsuite.multires.registration import register_roi_to_overview, sanitize_experiment_name
from lightsuite.multires.volume import manifest_geometry_report, write_sitk_hyperstack_tiff


def _status(message: str) -> None:
    print(message, flush=True)


def _volume_stem(path: Path) -> str:
    name = path.name
    if name.lower().endswith(".tif"):
        return name[:-4]
    if name.lower().endswith(".tiff"):
        return name[:-5]
    return path.stem


def _landmark_roi_report(roi, fit) -> dict[str, object]:
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


def _manifest_geometry_report(label: str, spec, *, manifest_dir: Path) -> dict[str, object]:
    return manifest_geometry_report(label, spec, manifest_dir=manifest_dir)


def check_multires_geometry(cfg: MultiresPipelineConfig) -> MultiresRegOptsCheckpoint:
    """Validate FOV overlap and write geometry QA artifacts."""
    manifest_path = cfg.multires.pair_manifest
    manifest = load_pair_manifest(manifest_path)
    manifest_dir = manifest_path.parent
    mode = cfg.multires.geometry_mode

    prepared = prepare_multires_registration_pair(cfg, manifest=manifest)
    overlap_min, overlap_max = prepared.overlap_box

    rep_overview = _manifest_geometry_report("overview", manifest.overview, manifest_dir=manifest_dir)
    if prepared.landmark_fit is not None:
        rep_roi = _landmark_roi_report(prepared.roi, prepared.landmark_fit)
    else:
        rep_roi = _manifest_geometry_report("roi", manifest.roi, manifest_dir=manifest_dir)

    geometry_dir = cfg.sample.save_path / "geometry"
    geometry_dir.mkdir(parents=True, exist_ok=True)

    report_path = geometry_dir / "geometry_report.json"
    serializable = {
        "geometry_mode": mode.value,
        "pair_label": manifest.pair_label,
        "pair_manifest_path": str(manifest_path),
        "overview": serialize_report(rep_overview),
        "roi": serialize_report(rep_roi),
        "overlap_box_um": [overlap_min.tolist(), overlap_max.tolist()],
    }
    if prepared.landmark_fit is not None:
        serializable["landmark_fit"] = {
            "rms_error_um": prepared.landmark_fit.rms_error_um,
            "fit_stats": prepared.landmark_fit.fit_stats,
            "roi_to_overview_tform": prepared.landmark_fit.roi_to_overview_tform.tolist(),
        }
    report_path.write_text(json.dumps(serializable, indent=2), encoding="utf-8")

    overview_path = Path(manifest.overview.volume_path)
    roi_path = Path(manifest.roi.volume_path)
    overview_stem = _volume_stem(overview_path)
    cropped_path = geometry_dir / f"{overview_stem}_cropped_overlap.tif"

    title = f"{manifest.pair_label} — {overview_stem} vs {_volume_stem(roi_path)} ({mode.value})"
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
    roi_tform = prepared.landmark_fit.roi_to_overview_tform if prepared.landmark_fit else None
    save_geometry_overlap_qc_plot(
        overview=prepared.overview,
        roi=prepared.roi,
        overlap_min=overlap_min,
        overlap_max=overlap_max,
        output_path=qc_plot_path,
        geometry_mode=mode.value,
        roi_to_overview=roi_tform,
    )

    write_sitk_hyperstack_tiff(cropped_path, prepared.fixed_cropped)

    experiment_slug = sanitize_experiment_name(cfg.multires.registration.experiment_name)
    landmark_session_path = None
    if mode != MultiresGeometryMode.METADATA:
        landmark_session_path = str(
            cfg.multires.resolved_landmark_session_path(cfg.sample.save_path, manifest)
        )

    checkpoint = MultiresRegOptsCheckpoint(
        sample_name=cfg.sample.name,
        pair_label=manifest.pair_label,
        pair_manifest_path=str(manifest_path),
        experiment_slug=experiment_slug,
        overview_volume_path=manifest.overview.volume_path,
        roi_volume_path=manifest.roi.volume_path,
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
    checkpoint.save(multires_checkpoint_path(cfg.sample.save_path))
    _status(f"Geometry QA ({mode.value}): {geometry_dir}")
    _status(f"Overlap QC plot: {qc_plot_path}")
    if prepared.landmark_fit is not None:
        _status(f"Landmark RMS error: {prepared.landmark_fit.rms_error_um:.2f} µm")
    return checkpoint


def run_multires_registration(cfg: MultiresPipelineConfig) -> MultiresRegOptsCheckpoint:
    """Register ROI stack to overview using a pair manifest and elastix."""
    manifest_path = cfg.multires.pair_manifest
    manifest = load_pair_manifest(manifest_path)
    meso = cfg.multires

    experiment_slug = sanitize_experiment_name(meso.registration.experiment_name)
    output_dir = cfg.sample.save_path / "elastix_roi_to_overview" / experiment_slug

    prepared = prepare_multires_registration_pair(cfg, manifest=manifest)
    overview_stem = _volume_stem(Path(manifest.overview.volume_path))
    roi_stem = _volume_stem(Path(manifest.roi.volume_path))

    result = register_roi_to_overview(
        prepared=prepared,
        output_dir=output_dir,
        experiment_slug=experiment_slug,
        overview_stem=overview_stem,
        roi_stem=roi_stem,
        registration_bin=meso.registration.registration_bin,
        elastix_stages=meso.registration.elastix_stages,
        write_full_overview_canvas=meso.registration.write_full_overview_canvas,
    )

    landmark_session_path = None
    if meso.geometry_mode != MultiresGeometryMode.METADATA:
        landmark_session_path = str(
            meso.resolved_landmark_session_path(cfg.sample.save_path, manifest)
        )

    checkpoint = MultiresRegOptsCheckpoint(
        sample_name=cfg.sample.name,
        pair_label=manifest.pair_label,
        pair_manifest_path=str(manifest_path),
        experiment_slug=experiment_slug,
        overview_volume_path=manifest.overview.volume_path,
        roi_volume_path=manifest.roi.volume_path,
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
    checkpoint.save(multires_checkpoint_path(cfg.sample.save_path))
    if result.registered_roi_full_overview_path is not None:
        _status(f"Full overview canvas: {result.registered_roi_full_overview_path}")
    _status(f"Registered ROI: {result.registered_roi_path}")
    return checkpoint
