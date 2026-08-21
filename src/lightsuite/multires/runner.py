"""Orchestration for manifest-driven multiresolution registration."""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np

from lightsuite.multires.checkpoint import MultiresRegOptsCheckpoint, multires_checkpoint_path
from lightsuite.multires.config_models import (
    MultiresGeometryCheckLevel,
    MultiresGeometryMode,
    MultiresPipelineConfig,
)
from lightsuite.multires.geometry import (
    physical_corners,
    transform_physical_points,
)
from lightsuite.multires.landmarks import fit_landmark_transform
from lightsuite.multires.memory import warn_if_overlap_memory_exceeds_system
from lightsuite.multires.models import MultiresPairManifest, serialize_report
from lightsuite.multires.plots import (
    normalized_cross_correlation,
    save_fov_overlap_plot,
    save_geometry_overlap_qc_plot,
    save_geometry_slice_qc_plot,
)
from lightsuite.multires.prepare import load_landmark_session, prepare_multires_registration_pair
from lightsuite.multires.registration import (
    apply_registration_to_channel,
    register_roi_to_overview,
    sanitize_experiment_name,
)
from lightsuite.multires.resolve import resolve_pair_manifest
from lightsuite.multires.spec_geometry import (
    alignment_metrics_from_specs,
    crop_index_range_from_physical_box,
    manifest_geometry_report_from_spec,
    overlap_box_from_landmark_specs,
    overlap_physical_bounds_from_specs,
    physical_to_continuous_index_xyz,
    sitk_geometry_from_spec,
    transformed_bounds_from_spec,
)
from lightsuite.multires.volume import load_manifest_xy_crop, write_sitk_hyperstack_tiff
from lightsuite.reporter import check_stage_cancelled, emit_pipeline_message, format_duration


def _status(message: str) -> None:
    emit_pipeline_message(message)


def _volume_stem(path: Path) -> str:
    name = path.name
    if name.lower().endswith(".tif"):
        return name[:-4]
    if name.lower().endswith(".tiff"):
        return name[:-5]
    return path.stem


def _landmark_roi_report_from_spec(roi_spec, fit) -> dict[str, object]:
    roi_min, roi_max = transformed_bounds_from_spec(roi_spec, fit.roi_to_overview_tform)
    center = 0.5 * (roi_min + roi_max)
    return {
        "label": "roi_in_overview_space",
        "phys_min": roi_min,
        "phys_max": roi_max,
        "phys_center": center,
        "landmark_rms_error_um": fit.rms_error_um,
        "landmark_fit_stats": fit.fit_stats,
    }


def _lightweight_geometry_context(
    cfg: MultiresPipelineConfig,
    manifest: MultiresPairManifest,
):
    """Compute overlap and optional landmark fit without loading full stacks."""
    mode = cfg.multires.geometry_mode
    margin_um = cfg.multires.registration.overlap_margin_um
    overview_spec = manifest.overview
    roi_spec = manifest.roi
    landmark_fit = None

    if mode == MultiresGeometryMode.METADATA:
        overlap_min, overlap_max = overlap_physical_bounds_from_specs(
            overview_spec,
            roi_spec,
            margin_um=margin_um,
        )
    else:
        session_path = cfg.multires.resolved_landmark_session_path(cfg.sample.save_path, manifest)
        session = load_landmark_session(session_path)
        overview_geo = sitk_geometry_from_spec(overview_spec)
        roi_geo = sitk_geometry_from_spec(roi_spec)
        landmark_fit = fit_landmark_transform(
            overview=overview_geo,
            roi=roi_geo,
            session=session,
            fit_mode=cfg.multires.landmarks.fit_mode,
            min_pairs=cfg.multires.landmarks.min_pairs,
        )
        overlap_min, overlap_max = overlap_box_from_landmark_specs(
            overview_spec,
            roi_spec,
            landmark_fit.roi_to_overview_tform,
            margin_um=margin_um,
        )

    crop_start_index, crop_size = crop_index_range_from_physical_box(
        overview_spec,
        overlap_min,
        overlap_max,
    )
    warn_if_overlap_memory_exceeds_system(
        crop_size,
        max_slab_bytes=cfg.multires.registration.max_slab_bytes,
    )
    return overlap_min, overlap_max, crop_start_index, landmark_fit


def _write_geometry_artifacts(
    *,
    cfg: MultiresPipelineConfig,
    manifest: MultiresPairManifest,
    manifest_dir: Path,
    overlap_min: np.ndarray,
    overlap_max: np.ndarray,
    crop_start_index: list[int],
    landmark_fit,
    level: MultiresGeometryCheckLevel,
    prepared=None,
) -> MultiresRegOptsCheckpoint:
    check_stage_cancelled()
    mode = cfg.multires.geometry_mode
    rep_overview = manifest_geometry_report_from_spec("overview", manifest.overview)
    if landmark_fit is not None:
        rep_roi = _landmark_roi_report_from_spec(manifest.roi, landmark_fit)
    else:
        rep_roi = manifest_geometry_report_from_spec("roi", manifest.roi)

    geometry_dir = cfg.sample.save_path / "geometry" / manifest.pair_label
    geometry_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = cfg.multires.resolved_pair_manifest_path(
        cfg.sample.save_path,
        cfg.sample.name,
    )

    report_path = geometry_dir / "geometry_report.json"
    alignment_metrics = alignment_metrics_from_specs(
        manifest.overview,
        manifest.roi,
        overlap_min=overlap_min,
        overlap_max=overlap_max,
    )
    serializable = {
        "geometry_mode": mode.value,
        "geometry_check_level": level.value,
        "pair_label": manifest.pair_label,
        "pair_manifest_path": str(manifest_path),
        "overview": serialize_report(rep_overview),
        "roi": serialize_report(rep_roi),
        "overlap_box_um": [overlap_min.tolist(), overlap_max.tolist()],
        "alignment_metrics": serialize_report(alignment_metrics),
    }
    if landmark_fit is not None:
        serializable["landmark_fit"] = {
            "rms_error_um": landmark_fit.rms_error_um,
            "fit_stats": landmark_fit.fit_stats,
            "roi_to_overview_tform": landmark_fit.roi_to_overview_tform.tolist(),
        }

    overview_path = Path(manifest.overview.volume_path)
    roi_path = Path(manifest.roi.volume_path)
    overview_stem = _volume_stem(overview_path)
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
    cropped_path = geometry_dir / f"{overview_stem}_overview_crop.tif"
    geometry_report_paths = {
        "geometry_report_json": str(report_path),
        "fov_overlap_png": str(fov_plot_path),
        "geometry_overlap_qc_png": str(qc_plot_path),
    }

    roi_tform = landmark_fit.roi_to_overview_tform if landmark_fit is not None else None
    if level == MultiresGeometryCheckLevel.SLICE_QC:
        check_stage_cancelled()
        emit_pipeline_message("Check-geometry: loading mid-plane overlap crops for slice QC…")
        center_um = tuple(float(v) for v in 0.5 * (overlap_min + overlap_max))
        overview_start, overview_size = crop_index_range_from_physical_box(
            manifest.overview,
            overlap_min,
            overlap_max,
        )
        if roi_tform is None:
            roi_start, roi_size = crop_index_range_from_physical_box(
                manifest.roi,
                overlap_min,
                overlap_max,
            )
            roi_physical_um = center_um
        else:
            inv = np.linalg.inv(roi_tform)
            roi_corners = transform_physical_points(physical_corners(overlap_min, overlap_max), inv)
            roi_phys_min = roi_corners.min(axis=0)
            roi_phys_max = roi_corners.max(axis=0)
            roi_start, roi_size = crop_index_range_from_physical_box(
                manifest.roi,
                roi_phys_min,
                roi_phys_max,
            )
            roi_physical_um = tuple(
                float(v)
                for v in (inv @ np.array([*center_um, 1.0], dtype=float))[:3]
            )
        overview_z = int(
            round(physical_to_continuous_index_xyz(manifest.overview, center_um)[2])
        )
        roi_z = int(
            round(physical_to_continuous_index_xyz(manifest.roi, roi_physical_um)[2])
        )
        sl_overview = load_manifest_xy_crop(
            manifest.overview,
            z_index=overview_z,
            start_xyz=overview_start,
            crop_size_xyz=overview_size,
            manifest_dir=manifest_dir,
        )
        sl_roi = load_manifest_xy_crop(
            manifest.roi,
            z_index=roi_z,
            start_xyz=roi_start,
            crop_size_xyz=roi_size,
            manifest_dir=manifest_dir,
        )
        try:
            from lightsuite.multires.plots import _normalize_panel, _resample_to_shape

            roi_panel = _resample_to_shape(
                _normalize_panel(sl_roi),
                _normalize_panel(sl_overview).shape,
            )
            alignment_metrics["slice_ncc"] = normalized_cross_correlation(
                roi_panel,
                _normalize_panel(sl_overview),
            )
        except ValueError:
            alignment_metrics["slice_ncc"] = None
        serializable["alignment_metrics"] = serialize_report(alignment_metrics)
        save_geometry_slice_qc_plot(
            sl_overview=sl_overview,
            sl_roi=sl_roi,
            center_um=center_um,
            output_path=qc_plot_path,
            geometry_mode=mode.value,
            alignment_metrics=alignment_metrics,
        )
        if alignment_metrics.get("slice_ncc") is not None:
            emit_pipeline_message(
                f"Check-geometry: slice NCC {alignment_metrics['slice_ncc']:.3f}"
            )
    elif level == MultiresGeometryCheckLevel.FULL:
        assert prepared is not None
        emit_pipeline_message("Check-geometry: writing full overlap QC plot and preview crop…")
        save_geometry_overlap_qc_plot(
            overview=prepared.fixed_cropped,
            roi=prepared.moving,
            overlap_min=overlap_min,
            overlap_max=overlap_max,
            output_path=qc_plot_path,
            geometry_mode=mode.value,
            roi_to_overview=roi_tform,
            overview_crop=prepared.fixed_cropped,
            roi_crop=prepared.moving,
        )
        write_sitk_hyperstack_tiff(cropped_path, prepared.fixed_cropped)
        geometry_report_paths["cropped_overview_preview"] = str(cropped_path)

    report_path.write_text(json.dumps(serializable, indent=2), encoding="utf-8")
    emit_pipeline_message(f"Check-geometry: wrote geometry report → {report_path.name}")

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
        crop_start_index=crop_start_index,
        roi_to_overview_tform=landmark_fit.roi_to_overview_tform.tolist()
        if landmark_fit is not None
        else None,
        landmark_rms_error_um=landmark_fit.rms_error_um if landmark_fit is not None else None,
        geometry_report_paths=geometry_report_paths,
    )
    checkpoint.save(multires_checkpoint_path(cfg.sample.save_path))
    _status(f"Geometry QA ({mode.value}, {level.value}): {geometry_dir}")
    _status(f"FOV overlap plot: {fov_plot_path}")
    if level in (MultiresGeometryCheckLevel.SLICE_QC, MultiresGeometryCheckLevel.FULL):
        _status(f"Overlap QC plot: {qc_plot_path}")
    if landmark_fit is not None:
        _status(f"Landmark RMS error: {landmark_fit.rms_error_um:.2f} µm")
    return checkpoint


def check_multires_geometry(
    cfg: MultiresPipelineConfig,
    *,
    level: MultiresGeometryCheckLevel = MultiresGeometryCheckLevel.FULL,
) -> MultiresRegOptsCheckpoint:
    """Validate FOV overlap and write geometry QA artifacts."""
    check_stage_cancelled()
    t0 = time.perf_counter()
    mode = cfg.multires.geometry_mode
    emit_pipeline_message(
        f"Check-geometry: starting ({level.value}, geometry mode {mode.value})…"
    )

    emit_pipeline_message("Check-geometry: resolving pair manifest…")
    check_stage_cancelled()
    manifest, manifest_path = resolve_pair_manifest(cfg)
    manifest_dir = manifest_path.parent
    overview_stem = _volume_stem(Path(manifest.overview.volume_path))
    roi_stem = _volume_stem(Path(manifest.roi.volume_path))
    emit_pipeline_message(
        f"Check-geometry: pair {manifest.pair_label} — {overview_stem} vs {roi_stem} "
        f"({format_duration(time.perf_counter() - t0)})"
    )

    if level in (MultiresGeometryCheckLevel.METADATA_ONLY, MultiresGeometryCheckLevel.SLICE_QC):
        if mode != MultiresGeometryMode.METADATA:
            emit_pipeline_message(
                "Check-geometry: loading landmarks and fitting coarse transform…"
            )
        emit_pipeline_message("Check-geometry: computing overlap bounds…")
        check_stage_cancelled()
        t_overlap = time.perf_counter()
        overlap_min, overlap_max, crop_start_index, landmark_fit = _lightweight_geometry_context(
            cfg,
            manifest,
        )
        overlap_msg = (
            f"crop start {crop_start_index}, "
            f"overlap size "
            f"{[float(overlap_max[i] - overlap_min[i]) for i in range(3)]} µm"
        )
        if landmark_fit is not None:
            overlap_msg += f", landmark RMS {landmark_fit.rms_error_um:.2f} µm"
        emit_pipeline_message(
            f"Check-geometry: overlap bounds ready in {format_duration(time.perf_counter() - t_overlap)} "
            f"({overlap_msg})"
        )
        emit_pipeline_message(f"Check-geometry: writing QA artifacts ({level.value})…")
        check_stage_cancelled()
        checkpoint = _write_geometry_artifacts(
            cfg=cfg,
            manifest=manifest,
            manifest_dir=manifest_dir,
            overlap_min=overlap_min,
            overlap_max=overlap_max,
            crop_start_index=crop_start_index,
            landmark_fit=landmark_fit,
            level=level,
        )
        emit_pipeline_message(
            f"Check-geometry: complete in {format_duration(time.perf_counter() - t0)}"
        )
        return checkpoint

    emit_pipeline_message(
        "Check-geometry: loading overview overlap crop and resampling ROI…"
    )
    check_stage_cancelled()
    t_prep = time.perf_counter()
    prepared = prepare_multires_registration_pair(cfg, manifest=manifest)
    overlap_min, overlap_max = prepared.overlap_box
    fixed_size = prepared.fixed_cropped.GetSize()
    moving_size = prepared.moving.GetSize()
    prep_msg = (
        f"overview crop {tuple(fixed_size)}, ROI resampled {tuple(moving_size)}"
    )
    if prepared.landmark_fit is not None:
        prep_msg += f", landmark RMS {prepared.landmark_fit.rms_error_um:.2f} µm"
    emit_pipeline_message(
        f"Check-geometry: overlap volumes ready in {format_duration(time.perf_counter() - t_prep)} "
        f"({prep_msg})"
    )
    emit_pipeline_message("Check-geometry: writing QA artifacts (full)…")
    check_stage_cancelled()
    checkpoint = _write_geometry_artifacts(
        cfg=cfg,
        manifest=manifest,
        manifest_dir=manifest_dir,
        overlap_min=overlap_min,
        overlap_max=overlap_max,
        crop_start_index=prepared.crop_start_index,
        landmark_fit=prepared.landmark_fit,
        level=MultiresGeometryCheckLevel.FULL,
        prepared=prepared,
    )
    emit_pipeline_message(f"Check-geometry: complete in {format_duration(time.perf_counter() - t0)}")
    return checkpoint


def run_multires_registration(cfg: MultiresPipelineConfig) -> MultiresRegOptsCheckpoint:
    """Register ROI stack to overview using a pair manifest and elastix."""
    check_stage_cancelled()
    manifest, manifest_path = resolve_pair_manifest(cfg)
    meso = cfg.multires

    experiment_slug = sanitize_experiment_name(meso.registration.experiment_name)
    output_dir = cfg.sample.save_path / "elastix_roi_to_overview" / experiment_slug

    reference_channel = (
        meso.registration.reference_channel
        or manifest.resolved_reference_channel()
    )
    if manifest.channels and reference_channel:
        ref_specs = manifest.channel_specs(reference_channel)
        overview_spec = ref_specs.overview
        roi_spec = ref_specs.roi
    else:
        overview_spec = manifest.overview
        roi_spec = manifest.roi

    prepared = prepare_multires_registration_pair(
        cfg,
        manifest=manifest,
        overview_spec=overview_spec,
        roi_spec=roi_spec,
    )
    check_stage_cancelled()
    overview_stem = _volume_stem(Path(overview_spec.volume_path))
    roi_stem = _volume_stem(Path(roi_spec.volume_path))

    result = register_roi_to_overview(
        prepared=prepared,
        output_dir=output_dir,
        experiment_slug=experiment_slug,
        overview_stem=overview_stem,
        roi_stem=roi_stem,
        registration_bin=meso.registration.registration_bin,
        elastix_stages=meso.registration.elastix_stages,
        write_full_overview_canvas=meso.registration.write_full_overview_canvas,
        pair_label=manifest.pair_label,
        channel=reference_channel,
    )

    additional_channels = manifest.non_reference_channels(
        reference_channel=reference_channel,
        apply_transform_to=meso.registration.apply_transform_to,
    )
    for channel_name in additional_channels:
        check_stage_cancelled()
        channel_specs = manifest.channel_specs(channel_name)
        channel_prepared = prepare_multires_registration_pair(
            cfg,
            manifest=manifest,
            overview_spec=channel_specs.overview,
            roi_spec=channel_specs.roi,
        )
        channel_result = apply_registration_to_channel(
            prepared=channel_prepared,
            transform_paths=result.transform_paths,
            output_dir=output_dir,
            experiment_slug=experiment_slug,
            channel=channel_name,
            overview_stem=_volume_stem(Path(channel_specs.overview.volume_path)),
            roi_stem=_volume_stem(Path(channel_specs.roi.volume_path)),
            write_full_overview_canvas=meso.registration.write_full_overview_canvas,
            pair_label=manifest.pair_label,
        )
        result.additional_channels.append(channel_result)
        _status(f"Applied transform to channel {channel_name}: {channel_result.registered_roi_path}")

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
        overview_volume_path=overview_spec.volume_path,
        roi_volume_path=roi_spec.volume_path,
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
        registration_overlay_qc_path=(
            str(result.registration_overlay_qc_path)
            if result.registration_overlay_qc_path is not None
            else None
        ),
        registration_slice_ncc=result.registration_slice_ncc,
        reference_channel=reference_channel,
        additional_channel_paths={
            channel_result.channel: str(channel_result.registered_roi_path)
            for channel_result in result.additional_channels
        }
        if result.additional_channels
        else None,
    )
    checkpoint.save(multires_checkpoint_path(cfg.sample.save_path))
    if result.registration_overlay_qc_path is not None:
        _status(f"Registration overlay QC: {result.registration_overlay_qc_path}")
        if result.registration_slice_ncc is not None:
            _status(f"Mid-plane NCC: {result.registration_slice_ncc:.3f}")
    if result.registered_roi_full_overview_path is not None:
        _status(f"Full overview canvas: {result.registered_roi_full_overview_path}")
    _status(f"Registered ROI: {result.registered_roi_path}")
    return checkpoint
