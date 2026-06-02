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
from lightsuite.registration.warp import swap_xy_transform
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
from lightsuite.registration.plots import save_registration_stage_previews
from lightsuite.registration.points_utils import thin_point_list
from lightsuite.registration.volume import load_registration_volume, permute_brain_volume
from lightsuite.registration.warp import warp_volume_affine

console = Console()


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

    def print_summary(self) -> None:
        console.print(
            f"[bold]Affine fit[/bold] ({self.n_total} pairs: "
            f"{self.n_manual} manual, {self.n_auto} auto)"
        )
        console.print(
            f"  Residual: median {self.median_error_vox:.1f} vox, "
            f"p95 {self.p95_error_vox:.1f}, max {self.max_error_vox:.1f} "
            f"(MSE {self.mse:.1f})"
        )
        if self.median_error_manual_vox is not None:
            console.print(f"  Manual pairs median: {self.median_error_manual_vox:.1f} vox")
        if self.median_error_auto_vox is not None:
            console.print(f"  Auto pairs median: {self.median_error_auto_vox:.1f} vox")
        if self.median_coarse_auto_vox is not None:
            console.print(
                f"  Auto pairs vs coarse only: {self.median_coarse_auto_vox:.1f} vox "
                "(before affine)"
            )
        if self.median_landmark_vox is not None:
            console.print(
                f"  B-spline landmarks median: {self.median_landmark_vox:.1f} vox "
                f"(cpwt={self.control_point_weight:g})"
            )
        if self.median_error_vox > 25.0:
            console.print(
                "[yellow]Warning:[/yellow] high affine residual — check orientation, "
                "add manual match-points, or set augment_points: true."
            )
        elif self.median_error_vox > 15.0:
            console.print(
                "[dim]Affine residual is moderate; manual landmarks often improve B-spline.[/dim]"
            )


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
    console.print(f"Augmentation with {n_auto} auto-extracted control points.")

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

    console.print("Loading registration volume...", end=" ")
    t0 = time.perf_counter()
    volume = load_registration_volume(Path(checkpoint.regvolpath))
    perm = checkpoint.permute_sample_to_atlas or [1, 2, 3]
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
    console.print(f"Done in {time.perf_counter() - t0:.1f}s.")

    tform_aff, cpaffine, cptshistology, cpwt, affine_diag = _prepare_control_points(
        checkpoint, session, config
    )
    affine_diag.print_summary()
    affine_diag.save(save_path / "affine_fit_stats.json")
    console.print(f"Using {cptshistology.shape[0]} landmark pairs for B-spline.")

    atlas = resolve_brain_atlas(config.atlas.provider.value, config.atlas.atlas_dir)
    tv = np.asanyarray(nib.load(atlas.template_path).dataobj).astype(np.float32)
    av = np.asanyarray(nib.load(atlas.annotation_path).dataobj).astype(np.float32)

    # Match multiobjRegistration.m: affine is fit in full-atlas index space (points / downfac_reg)
    # and imwarp uses full-resolution moving volumes, not the downsampled tvreg grid.
    console.print("Applying affine pre-alignment to atlas...", end=" ")
    t0 = time.perf_counter()
    tvaffine = warp_volume_affine(tv, tform_aff, volume.shape, order=1)
    avaffine = warp_volume_affine(av, tform_aff, volume.shape, order=0)
    console.print(f"Done in {time.perf_counter() - t0:.1f}s.")

    hi = float(np.quantile(volume, 0.999))
    voltoshow = np.clip(volume / max(hi, 1e-6) * 255.0, 0, 255).astype(np.uint8)
    save_registration_stage_previews(
        save_path, config.sample.name, voltoshow, avaffine, "affine_registration"
    )

    elastix_temp = save_path / "elastix_temp"
    moving_pts_mm = volume_indices_to_elastix_physical(cpaffine, spacing_mm)
    fixed_pts_mm = volume_indices_to_elastix_physical(cptshistology, spacing_mm)

    bspline_result = run_bspline_registration(
        fixed_volume=volume,
        moving_volume=tvaffine,
        fixed_secondary=volume_secondary,
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

    console.print("Warping annotation with B-spline transform...", end=" ")
    t0 = time.perf_counter()
    avreg = run_transformix(
        moving_volume=np.rint(avaffine).astype(np.int32),
        transform_path=bspline_result.transform_path,
        output_dir=save_path / "transformix_annotation_temp",
        spacing_mm=spacing_mm,
        nearest=True,
    )
    final_landmark_mm = read_elastix_landmark_metric_mm(bspline_result.output_dir)
    landmark_note = (
        f", final landmark metric {final_landmark_mm:.3f} mm"
        if final_landmark_mm is not None
        else ""
    )
    console.print(
        f"Done in {time.perf_counter() - t0:.1f}s "
        f"({int(np.sum(np.rint(avreg) > 1)):,} label voxels after B-spline{landmark_note})."
    )
    save_registration_stage_previews(
        save_path, config.sample.name, voltoshow, avreg, "bspline_registration"
    )

    inverse_dir = save_path / "elastix_inverse_temp"
    inverted = invert_elastix_transform(elastix_temp, inverse_dir)
    samp_to_atlas_path = save_path / "bspline_samp_to_atlas_20um.txt"
    write_inverted_transform_copy(inverted, samp_to_atlas_path)

    affine_inv = np.linalg.inv(tform_aff)
    reg = config.registration
    transform_params = TransformParamsCheckpoint(
        atlas_resolution_um=config.atlas.resolution_um,
        regvolsize=regvolsize,
        atlassize=list(tv.shape),
        brain_atlas=config.atlas.provider.value,
        ori_voxel_um=checkpoint.voxel_um,
        ori_size=[checkpoint.ny, checkpoint.nx, checkpoint.nz],
        permute_sample_to_atlas=checkpoint.permute_sample_to_atlas or [1, 2, 3],
        elastix_um_to_mm=1e-3,
        tform_bspline_samp20um_to_atlas_20um_px=str(samp_to_atlas_path),
        tform_affine_samp20um_to_atlas_10um_px=affine_inv.tolist(),
        control_point_weight=cpwt,
        use_multistep=use_multistep,
        use_dual_channel_mi=volume_secondary is not None,
        dual_channel_mi_weight_autofluor=reg.dual_channel_mi_weight_autofluor
        if volume_secondary is not None
        else None,
        dual_channel_mi_weight_signal=reg.dual_channel_mi_weight_signal
        if volume_secondary is not None
        else None,
        channel_secondary=checkpoint.channel_secondary,
    )
    out_json = save_path / "transform_params.json"
    transform_params.save(out_json)
    console.print(f"[green]Registration complete.[/green] Transform params: {out_json}")
    return out_json
