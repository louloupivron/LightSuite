"""Brain registration stage (multiobjRegistration.m port)."""

from __future__ import annotations

import json
import shutil
import time
from dataclasses import asdict, dataclass
from pathlib import Path

import nibabel as nib
import numpy as np
from rich.console import Console
from scipy.spatial.distance import cdist

from lightsuite.atlas.registry import resolve_brain_atlas
from lightsuite.config.models import BrainPipelineConfig
from lightsuite.gui.affine import (
    affine_point_errors,
    fit_affine_transform,
    summarize_point_errors,
    transform_points,
    transform_points_inverse,
)
from lightsuite.registration.points import cloud_xyz_to_volume_indices
from lightsuite.registration.warp import swap_xy_transform, warp_sample_to_atlas
from lightsuite.gui.control_points import (
    ControlPointSession,
    default_session_path,
    load_registration_control_point_session,
)
from lightsuite.preprocess.checkpoint import RegOptsCheckpoint
from lightsuite.registration.elastix.invert import (
    invert_elastix_transform,
    write_inverted_transform_copy,
)
from lightsuite.registration.elastix.points import volume_indices_to_elastix_physical
from lightsuite.registration.elastix.runner import (
    read_elastix_landmark_metric_mm,
    run_bspline_registration,
    run_transformix,
)
from lightsuite.gui.brain_data import load_slice_correspondence
from lightsuite.registration.correspondence_affine import (
    append_correspondence_bspline_landmarks,
    apply_slice_correspondence_affine,
)
from lightsuite.registration.plots import save_registration_stage_previews
from lightsuite.registration.points_utils import thin_point_list
from lightsuite.registration.register_diagnostics import (
    RegistrationDiagnostics,
    classify_registration_status,
    landmark_mm_to_vox,
)
from lightsuite.registration.volume import (
    load_registration_volume,
    permute_brain_volume,
    resize_atlas_volume,
)
from lightsuite.registration.canvas import (
    compute_vd_warp_canvas_padding,
    crop_from_warp_canvas,
    offset_volume_indices,
    pad_volume_for_warp_canvas,
)
from lightsuite.registration.warp import warp_volume_affine

console = Console()
REGISTRATION_DIAGNOSTICS_FILENAME = "registration_diagnostics.json"


@dataclass
class AffineFitDiagnostics:
    """Affine control-point fit quality (voxel units at registration resolution)."""

    n_manual: int
    n_auto: int
    n_total: int
    mse: float
    median_error_vox: float
    p95_error_vox: float
    max_error_vox: float
    mean_error_vox: float
    median_error_manual_vox: float | None = None
    median_error_auto_vox: float | None = None
    median_coarse_auto_vox: float | None = None
    median_landmark_vox: float | None = None
    control_point_weight: float = 0.0

    def to_dict(self) -> dict:
        return asdict(self)

    def save(self, path: Path) -> None:
        path = path.expanduser()
        path.write_text(json.dumps(self.to_dict(), indent=2), encoding="utf-8")

@dataclass
class TransformParamsCheckpoint:
    """Registration outputs (transform_params.mat equivalent)."""

    atlas_resolution_um: float
    regvolsize: list[int]
    atlassize: list[int]
    brain_atlas: str
    ori_voxel_um: list[float]
    ori_size: list[int]
    permute_sample_to_atlas: list[int]
    elastix_um_to_mm: float
    tform_bspline_samp20um_to_atlas_20um_px: str
    tform_affine_samp20um_to_atlas_10um_px: list[list[float]]
    control_point_weight: float
    use_multistep: bool
    use_dual_channel_mi: bool
    dual_channel_mi_weight_autofluor: float | None = None
    dual_channel_mi_weight_signal: float | None = None
    channel_secondary: int | None = None
    warp_canvas_pad_before: list[int] | None = None

    def save(self, path: Path) -> None:
        path = path.expanduser()
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(asdict(self), indent=2), encoding="utf-8")

    @classmethod
    def load(cls, path: Path) -> TransformParamsCheckpoint:
        raw = json.loads(path.expanduser().read_text(encoding="utf-8"))
        return cls(**raw)


def validate_registration_inputs(
    config: BrainPipelineConfig,
) -> tuple[RegOptsCheckpoint, ControlPointSession]:
    save_path = config.sample.save_path.expanduser()
    regopts_path = save_path / "regopts.json"
    cp_path = default_session_path(save_path)

    if not regopts_path.is_file():
        msg = f"Missing {regopts_path}. Run preprocess and init-registration first."
        raise FileNotFoundError(msg)

    checkpoint = RegOptsCheckpoint.load(regopts_path)

    if checkpoint.original_trans is None:
        msg = "regopts.json missing original_trans from init-registration."
        raise RuntimeError(msg)

    if not cp_path.is_file():
        console.print(
            "[dim]No atlas2histology_tform.json — using auto control points from "
            "init-registration only (match-points is optional).[/dim]"
        )
    session = load_registration_control_point_session(
        save_path,
        original_trans=checkpoint.original_trans,
    )

    atlas_pts, sample_pts = session.paired_points_xyz()
    n_auto = len(checkpoint.autocpatlas or [])
    if atlas_pts.shape[0] < 4 and n_auto < 4:
        msg = (
            "Need at least 4 control points from init-registration (autocpatlas) "
            "or manual match-points."
        )
        raise RuntimeError(msg)
    if atlas_pts.shape[0] < 4:
        console.print(
            "[yellow]Warning:[/yellow] no manual control points — "
            "affine and B-spline use init-registration auto pairs only."
        )
    elif sample_pts.shape[0] < 4:
        console.print(
            "[yellow]Warning:[/yellow] fewer than 4 paired manual control points — "
            "registration quality may be poor."
        )

    if shutil.which("elastix") is None or shutil.which("transformix") is None:
        msg = "elastix and transformix must be on PATH for brain register."
        raise RuntimeError(msg)

    return checkpoint, session


def _build_affine_fit_diagnostics(
    *,
    af_atlas: np.ndarray,
    af_sample: np.ndarray,
    tform_aff: np.ndarray,
    n_manual: int,
    n_auto: int,
    autocpsample_kept: np.ndarray,
    original_trans: np.ndarray,
    cpaffine: np.ndarray,
    cptshistology: np.ndarray,
    cpwt: float,
) -> AffineFitDiagnostics:
    mse, errors = affine_point_errors(af_atlas, af_sample, tform_aff)
    summary = summarize_point_errors(errors)
    diag = AffineFitDiagnostics(
        n_manual=n_manual,
        n_auto=n_auto,
        n_total=int(af_atlas.shape[0]),
        mse=mse,
        median_error_vox=summary["median"],
        p95_error_vox=summary["p95"],
        max_error_vox=summary["max"],
        mean_error_vox=summary["mean"],
        control_point_weight=cpwt,
    )
    if n_manual > 0:
        manual_err = errors[:n_manual]
        diag.median_error_manual_vox = float(np.median(manual_err))
    if n_auto > 0:
        auto_err = errors[n_manual:]
        diag.median_error_auto_vox = float(np.median(auto_err))
        if autocpsample_kept.shape[0] == n_auto:
            coarse_atlas = transform_points_inverse(
                autocpsample_kept,
                swap_xy_transform(original_trans),
            )
            diag.median_coarse_auto_vox = float(
                np.median(np.linalg.norm(coarse_atlas - af_atlas[n_manual:], axis=1))
            )
    if cpaffine.shape[0] > 0 and cptshistology.shape[0] == cpaffine.shape[0]:
        diag.median_landmark_vox = float(
            np.median(np.linalg.norm(cpaffine - cptshistology, axis=1))
        )
    return diag


def _prepare_control_points(
    checkpoint: RegOptsCheckpoint,
    session: ControlPointSession,
    config: BrainPipelineConfig,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, float, AffineFitDiagnostics]:
    downfac = float(
        checkpoint.downfac_reg or (config.atlas.resolution_um / checkpoint.registres_um)
    )
    atlas_res = config.atlas.resolution_um
    distancethin = 1000.0 / atlas_res
    if not config.registration.augment_points:
        distancethin *= 3.0

    original_trans_vol = swap_xy_transform(
        np.asarray(
            session.ori_trans if session.ori_trans is not None else checkpoint.original_trans,
            dtype=float,
        )
    )

    autocpsample = cloud_xyz_to_volume_indices(
        np.asarray(checkpoint.autocpsample or [], dtype=float)
    )
    autocpatlas = cloud_xyz_to_volume_indices(
        np.asarray(checkpoint.autocpatlas or [], dtype=float)
    ) / downfac

    cptsatlas, cptshistology = session.paired_points_xyz()
    cptshistology = transform_points_inverse(cptshistology, original_trans_vol)
    cptsatlas = cptsatlas / downfac

    if cptsatlas.shape[0] > 0 and autocpatlas.shape[0] > 0:
        distances = cdist(cptsatlas, autocpatlas)
        keep_auto = np.all(distances > distancethin, axis=0)
        autocpsample = autocpsample[keep_auto]
        autocpatlas = autocpatlas[keep_auto]

    if autocpatlas.shape[0]:
        keep_idx = thin_point_list(autocpatlas, distancethin)
    else:
        keep_idx = np.zeros(0, dtype=bool)
    n_auto = int(np.count_nonzero(keep_idx))

    if keep_idx.size:
        af_atlas = np.vstack([cptsatlas, autocpatlas[keep_idx]])
        af_sample = np.vstack([cptshistology, autocpsample[keep_idx]])
    else:
        af_atlas = cptsatlas
        af_sample = cptshistology

    cpwt = config.registration.control_point_weight
    if af_atlas.shape[0] < 4:
        msg = "Need at least 4 control points (user or auto) for affine registration."
        raise RuntimeError(msg)

    tform_aff, _ = fit_affine_transform(af_atlas, af_sample)
    n_manual = int(cptsatlas.shape[0])
    if n_manual > 0:
        cpaffine = transform_points(cptsatlas, tform_aff)
        autocpsample_kept = autocpsample[keep_idx] if keep_idx.size else np.zeros((0, 3))
    else:
        # MATLAB multiobjRegistration.m "revert to automated mode": no manual points,
        # so drive the B-spline with the auto landmarks at half the control-point weight
        # while keeping the same image-dominant multistep schedule.
        cpwt = cpwt / 2.0
        autocpsample_kept = autocpsample[keep_idx]
        cptshistology = autocpsample_kept
        cpaffine = transform_points(autocpatlas[keep_idx], tform_aff)

    diagnostics = _build_affine_fit_diagnostics(
        af_atlas=af_atlas,
        af_sample=af_sample,
        tform_aff=tform_aff,
        n_manual=n_manual,
        n_auto=n_auto,
        autocpsample_kept=autocpsample_kept,
        original_trans=original_trans_vol,
        cpaffine=cpaffine,
        cptshistology=cptshistology,
        cpwt=cpwt,
    )

    return tform_aff, cpaffine, cptshistology, cpwt, diagnostics


def run_brain_registration(config: BrainPipelineConfig, *, use_multistep: bool = True) -> Path:
    """Run elastix B-spline registration and write transform_params.json."""
    checkpoint, session = validate_registration_inputs(config)
    save_path = config.sample.save_path.expanduser()
    spacing_mm = checkpoint.registres_um * 1e-3
    perm = checkpoint.permute_sample_to_atlas or [1, 2, 3]

    t0 = time.perf_counter()
    volume = load_registration_volume(Path(checkpoint.regvolpath))
    volume = permute_brain_volume(volume.astype(np.float32), perm)
    regvolsize = list(volume.shape)

    volume_secondary = None
    if checkpoint.regvolpath_secondary:
        volume_secondary = load_registration_volume(Path(checkpoint.regvolpath_secondary))
        volume_secondary = permute_brain_volume(
            volume_secondary.astype(np.float32),
            checkpoint.permute_sample_to_atlas or [1, 2, 3],
        )
        if volume_secondary.shape != volume.shape:
            msg = (
                f"Secondary volume shape {volume_secondary.shape} must match "
                f"primary {volume.shape}"
            )
            raise ValueError(msg)
    load_elapsed = time.perf_counter() - t0

    tform_aff, cpaffine, cptshistology, cpwt, affine_diag = _prepare_control_points(
        checkpoint, session, config
    )
    affine_diag.save(save_path / "affine_fit_stats.json")

    atlas = resolve_brain_atlas(config.atlas.provider.value, config.atlas.atlas_dir)
    tv = np.asanyarray(nib.load(atlas.template_path).dataobj).astype(np.float32)
    av = np.asanyarray(nib.load(atlas.annotation_path).dataobj).astype(np.float32)

    downfac = float(
        checkpoint.downfac_reg or (config.atlas.resolution_um / checkpoint.registres_um)
    )
    tvreg_shape = resize_atlas_volume(tv, downfac, nearest=False).shape
    original_trans_xyz = np.asarray(
        session.ori_trans if session.ori_trans is not None else checkpoint.original_trans,
        dtype=float,
    )
    sample_warped = warp_sample_to_atlas(volume, original_trans_xyz, tvreg_shape, order=1)
    correspondence = load_slice_correspondence(save_path)
    if (
        correspondence is not None
        and correspondence.has_confirmed_anchors()
        and checkpoint.original_trans is not None
    ):
        stored = np.asarray(correspondence.original_trans, dtype=float)
        current = np.asarray(checkpoint.original_trans, dtype=float)
        if not np.allclose(stored, current, rtol=0, atol=1e-3):
            console.print(
                "[yellow]Warning:[/yellow] slice_correspondence original_trans differs from "
                "regopts.json — correspondence affine uses current regopts transform."
            )
    tform_aff_before_corr = tform_aff.copy()
    tform_aff, corr_affine_stats = apply_slice_correspondence_affine(
        tform_aff,
        correspondence,
        sample_warped=sample_warped,
        original_trans=original_trans_xyz,
        downfac=downfac,
        enabled=config.registration.use_slice_correspondence_affine,
    )
    corr_affine_stats.save(save_path / "correspondence_affine_stats.json")
    if corr_affine_stats.applied:
        landmark_atlas = transform_points_inverse(cpaffine, tform_aff_before_corr)
        cpaffine = transform_points(landmark_atlas, tform_aff)
        console.print(
            "[green]Correspondence affine:[/green] "
            f"{corr_affine_stats.n_anchor_pairs} anchor pairs across "
            f"{corr_affine_stats.n_axes} axes · "
            f"median residual {corr_affine_stats.median_residual_before_vox:.2f} → "
            f"{corr_affine_stats.median_residual_after_vox:.2f} vox"
        )
    elif config.registration.use_slice_correspondence_affine and correspondence is not None:
        console.print(
            f"[dim]Slice correspondence affine not applied"
            f" ({corr_affine_stats.skip_reason}).[/dim]"
        )

    landmark_thin_vox = 1000.0 / config.atlas.resolution_um
    cpaffine, cptshistology, corr_landmark_stats = append_correspondence_bspline_landmarks(
        cpaffine,
        cptshistology,
        correspondence,
        sample_warped=sample_warped,
        original_trans=original_trans_xyz,
        downfac=downfac,
        tform_aff=tform_aff,
        min_distance_vox=landmark_thin_vox,
        max_landmarks=config.registration.correspondence_landmark_max_count,
        enabled=config.registration.use_slice_correspondence_landmarks,
    )
    corr_landmark_stats.save(save_path / "correspondence_landmark_stats.json")
    if corr_landmark_stats.applied:
        cpwt = max(cpwt, config.registration.correspondence_landmark_weight)
        console.print(
            "[green]Correspondence landmarks:[/green] "
            f"+{corr_landmark_stats.n_added_pairs} B-spline pairs "
            f"({corr_landmark_stats.n_merged_pairs} total, cpwt={cpwt:g})"
        )
    elif (
        config.registration.use_slice_correspondence_landmarks
        and correspondence is not None
    ):
        console.print(
            f"[dim]Slice correspondence landmarks not added"
            f" ({corr_landmark_stats.skip_reason}).[/dim]"
        )

    # Match multiobjRegistration.m: affine is fit in full-atlas index space (points / downfac_reg)
    # and imwarp uses full-resolution moving volumes, not the downsampled tvreg grid.
    warp_canvas = compute_vd_warp_canvas_padding(volume.shape, tv.shape)
    warp_shape = warp_canvas.padded_shape(volume.shape)
    if not warp_canvas.is_zero:
        console.print(
            "[dim]Warp canvas:[/dim] VD padding "
            f"{warp_canvas.pad_before[2]}+{warp_canvas.pad_after[2]} vox "
            f"({volume.shape} → {warp_shape})"
        )

    t0 = time.perf_counter()
    tvaffine = warp_volume_affine(
        tv,
        tform_aff,
        warp_shape,
        order=1,
        output_origin=warp_canvas.pad_before,
    )
    avaffine = warp_volume_affine(
        av,
        tform_aff,
        warp_shape,
        order=0,
        output_origin=warp_canvas.pad_before,
    )
    affine_warp_elapsed = time.perf_counter() - t0

    volume_padded = pad_volume_for_warp_canvas(volume, warp_canvas)
    volume_secondary_padded = (
        pad_volume_for_warp_canvas(volume_secondary, warp_canvas)
        if volume_secondary is not None
        else None
    )
    hi = float(np.quantile(volume, 0.999))
    voltoshow = np.clip(volume / max(hi, 1e-6) * 255.0, 0, 255).astype(np.uint8)
    voltoshow_padded = pad_volume_for_warp_canvas(voltoshow, warp_canvas)
    save_registration_stage_previews(
        save_path,
        config.sample.name,
        voltoshow_padded,
        avaffine,
        "affine_registration",
        atlas_provider=config.atlas.provider.value,
    )

    elastix_temp = save_path / "elastix_temp"
    cptshistology_padded = offset_volume_indices(cptshistology, warp_canvas.pad_before)
    moving_pts_mm = volume_indices_to_elastix_physical(cpaffine, spacing_mm)
    fixed_pts_mm = volume_indices_to_elastix_physical(cptshistology_padded, spacing_mm)

    t0 = time.perf_counter()
    bspline_result = run_bspline_registration(
        fixed_volume=volume_padded,
        moving_volume=tvaffine,
        fixed_secondary=volume_secondary_padded,
        moving_points_mm=moving_pts_mm,
        fixed_points_mm=fixed_pts_mm,
        output_dir=elastix_temp,
        save_path=save_path,
        spacing_mm=spacing_mm,
        control_point_weight=cpwt,
        n_histogram_bins=48,
        bspline_spatial_scale_mm=config.registration.bspline_spatial_scale_mm,
        use_multistep=use_multistep,
        dual_weight_autofluor=config.registration.dual_channel_mi_weight_autofluor,
        dual_weight_signal=config.registration.dual_channel_mi_weight_signal,
    )
    bspline_elapsed = time.perf_counter() - t0

    t0 = time.perf_counter()
    avreg_padded = run_transformix(
        moving_volume=np.rint(avaffine).astype(np.int32),
        transform_path=bspline_result.transform_path,
        output_dir=save_path / "transformix_annotation_temp",
        spacing_mm=spacing_mm,
        nearest=True,
    )
    transformix_elapsed = time.perf_counter() - t0
    avreg = crop_from_warp_canvas(avreg_padded, warp_canvas, volume.shape)
    annotation_label_voxels = int(np.sum(np.rint(avreg) > 1))
    final_landmark_mm = read_elastix_landmark_metric_mm(bspline_result.output_dir)
    final_landmark_vox = (
        landmark_mm_to_vox(final_landmark_mm, checkpoint.registres_um)
        if final_landmark_mm is not None
        else None
    )

    save_registration_stage_previews(
        save_path,
        config.sample.name,
        voltoshow_padded,
        avreg_padded,
        "bspline_registration",
        atlas_provider=config.atlas.provider.value,
    )

    inverse_dir = save_path / "elastix_inverse_temp"
    inverted = invert_elastix_transform(elastix_temp, inverse_dir)
    samp_to_atlas_path = save_path / "bspline_samp_to_atlas_20um.txt"
    write_inverted_transform_copy(inverted, samp_to_atlas_path)

    affine_inv = np.linalg.inv(tform_aff)
    reg = config.registration
    use_dual = volume_secondary is not None
    status, status_message, warnings = classify_registration_status(
        affine_median_error_vox=affine_diag.median_error_vox,
        affine_p95_error_vox=affine_diag.p95_error_vox,
        n_manual=affine_diag.n_manual,
        n_landmark_pairs=int(cptshistology.shape[0]),
        bspline_landmark_metric_vox=final_landmark_vox,
        annotation_label_voxels=annotation_label_voxels,
        use_multistep=use_multistep,
    )
    diagnostics = RegistrationDiagnostics(
        sample_shape=regvolsize,
        atlas_shape=list(tv.shape),
        orientation=list(perm),
        registration_resolution_um=float(checkpoint.registres_um),
        n_manual_pairs=affine_diag.n_manual,
        n_auto_pairs=affine_diag.n_auto,
        n_landmark_pairs=int(cptshistology.shape[0]),
        control_point_weight=cpwt,
        use_multistep=use_multistep,
        use_dual_channel_mi=use_dual,
        bspline_spatial_scale_mm=float(reg.bspline_spatial_scale_mm),
        dual_channel_mi_weight_autofluor=reg.dual_channel_mi_weight_autofluor if use_dual else None,
        dual_channel_mi_weight_signal=reg.dual_channel_mi_weight_signal if use_dual else None,
        affine_median_error_vox=affine_diag.median_error_vox,
        affine_p95_error_vox=affine_diag.p95_error_vox,
        affine_max_error_vox=affine_diag.max_error_vox,
        affine_median_manual_vox=affine_diag.median_error_manual_vox,
        affine_median_auto_vox=affine_diag.median_error_auto_vox,
        affine_median_coarse_auto_vox=affine_diag.median_coarse_auto_vox,
        bspline_landmark_metric_mm=final_landmark_mm,
        bspline_landmark_metric_vox=final_landmark_vox,
        annotation_label_voxels=annotation_label_voxels,
        load_elapsed_s=load_elapsed,
        affine_warp_elapsed_s=affine_warp_elapsed,
        bspline_elapsed_s=bspline_elapsed,
        transformix_elapsed_s=transformix_elapsed,
        status=status,
        status_message=status_message,
        warnings=warnings,
    )
    diag_path = save_path / REGISTRATION_DIAGNOSTICS_FILENAME
    diagnostics.save(diag_path)
    diagnostics.print_summary(console=console)

    transform_params = TransformParamsCheckpoint(
        atlas_resolution_um=config.atlas.resolution_um,
        regvolsize=regvolsize,
        atlassize=list(tv.shape),
        brain_atlas=config.atlas.provider.value,
        ori_voxel_um=checkpoint.voxel_um,
        ori_size=[checkpoint.ny, checkpoint.nx, checkpoint.nz],
        permute_sample_to_atlas=perm,
        elastix_um_to_mm=1e-3,
        tform_bspline_samp20um_to_atlas_20um_px=str(samp_to_atlas_path),
        tform_affine_samp20um_to_atlas_10um_px=affine_inv.tolist(),
        control_point_weight=cpwt,
        use_multistep=use_multistep,
        use_dual_channel_mi=use_dual,
        dual_channel_mi_weight_autofluor=reg.dual_channel_mi_weight_autofluor
        if use_dual
        else None,
        dual_channel_mi_weight_signal=reg.dual_channel_mi_weight_signal
        if use_dual
        else None,
        channel_secondary=checkpoint.channel_secondary,
        warp_canvas_pad_before=list(warp_canvas.pad_before) if not warp_canvas.is_zero else None,
    )
    out_json = save_path / "transform_params.json"
    transform_params.save(out_json)
    console.print(
        f"[dim]Checkpoint {out_json.name} · diagnostics {diag_path.name} · "
        f"affine_fit_stats.json · load {load_elapsed:.1f}s[/dim]"
    )
    return out_json
