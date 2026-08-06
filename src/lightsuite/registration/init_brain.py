"""Initialize brain registration (initializeRegistration.m port)."""

from __future__ import annotations

import time
from pathlib import Path

from lightsuite.atlas.io import load_atlas_volume
from lightsuite.atlas.trim import save_atlas_manifest_copy
import numpy as np
from rich.console import Console

from lightsuite.atlas.registry import atlas_display_provider_from_config, resolve_brain_atlas_content
from lightsuite.config.models import BrainPipelineConfig
from lightsuite.preprocess.checkpoint import RegOptsCheckpoint
from lightsuite.registration.align import (
    coarse_alignment_metrics,
    downsample_point_cloud,
    estimate_similarity_transform,
    similarity_scale,
    triage_and_match_clouds,
)
from lightsuite.registration.bcpd import find_bcpd_executable
from lightsuite.registration.init_diagnostics import (
    InitRegistrationDiagnostics,
    classify_init_registration_status,
)
from lightsuite.registration.orientation import orientation_path, resolve_orientation, save_orientation
from lightsuite.registration.plots import save_initial_registration_previews
from lightsuite.registration.points import (
    extract_atlas_points_gradient,
    extract_sample_points,
    extract_sample_points_stages,
)
from lightsuite.registration.volume import (
    load_registration_volume,
    normalize_registration_volume,
    permute_brain_volume,
    resize_atlas_volume,
)

console = Console()
INIT_DIAGNOSTICS_FILENAME = "init_registration_diagnostics.json"


def initialize_brain_registration(config: BrainPipelineConfig) -> RegOptsCheckpoint:
    """Coarse-align sample to atlas and update regopts checkpoint."""
    save_path = config.sample.save_path.expanduser()
    regopts_path = save_path / "regopts.json"
    if not regopts_path.is_file():
        msg = f"Missing checkpoint {regopts_path}. Run 'lightsuite brain preprocess' first."
        raise FileNotFoundError(msg)

    checkpoint = RegOptsCheckpoint.load(regopts_path)
    backvol = load_registration_volume(Path(checkpoint.regvolpath))
    downfac = config.atlas.resolution_um / checkpoint.registres_um

    atlas_content = resolve_brain_atlas_content(
        config.atlas,
        scratch=config.sample.scratch,
    )
    atlas = atlas_content.paths
    tv = load_atlas_volume(atlas.template_path)
    av = load_atlas_volume(atlas.annotation_path)
    tvreg = resize_atlas_volume(tv.astype(np.float32), downfac, nearest=False)
    avreg = resize_atlas_volume(av.astype(np.float32), downfac, nearest=True)

    if atlas.boundary_path is not None:
        boundary_full = load_atlas_volume(atlas.boundary_path)
    else:
        from lightsuite.registration.plots import boundary_volume_from_annotation

        boundary_full = boundary_volume_from_annotation(av)
    boundary_reg = resize_atlas_volume(boundary_full.astype(np.float32), downfac, nearest=True)
    boundary_reg = ((boundary_reg > 0).astype(np.uint8)) * 255

    permvec = resolve_orientation(config, save_path)
    orient_file = orientation_path(save_path)
    if config.registration.orientation is not None:
        save_orientation(save_path, permvec)
    elif not orient_file.is_file():
        console.print(
            "[dim]No brain_orientation.txt — using default [1, 2, 3]. "
            "Run check-orientation to verify axes.[/dim]"
        )

    newvol = normalize_registration_volume(backvol)
    t0 = time.perf_counter()
    volumereg = permute_brain_volume(newvol, permvec)
    ls_cloud, sample_stages = extract_sample_points_stages(
        volumereg,
        config.registration.cloud_threshold,
        subsample_fraction=config.registration.sample_cloud_subsample,
    )
    sample_cloud_elapsed = time.perf_counter() - t0

    t0 = time.perf_counter()
    tv_for_points = tvreg.copy()
    tv_for_points[avreg == 0] = 0
    tv_cloud = extract_atlas_points_gradient(tv_for_points, avreg, sigma=20.0, threshold=5.0)
    atlas_cloud_elapsed = time.perf_counter() - t0

    if ls_cloud.shape[0] < 10 or tv_cloud.shape[0] < 10:
        msg = "Too few points extracted for coarse registration."
        raise RuntimeError(msg)

    bcpd_path = find_bcpd_executable(config.registration.bcpd_path)
    if bcpd_path is None and tv_cloud.shape[0] > 100_000:
        tv_cloud = downsample_point_cloud(tv_cloud, 100_000)

    t0 = time.perf_counter()
    transform_icp, transform_matlab, backend = estimate_similarity_transform(
        tv_cloud,
        ls_cloud,
        bcpd_path=bcpd_path,
    )
    alignment_elapsed = time.perf_counter() - t0
    scale = similarity_scale(transform_icp)
    metrics = coarse_alignment_metrics(ls_cloud, tv_cloud, transform_icp)

    t0 = time.perf_counter()
    cpsample, cpatlas = triage_and_match_clouds(
        ls_cloud,
        tv_cloud,
        transform_icp,
        bcpd_path=bcpd_path,
    )
    triage_elapsed = time.perf_counter() - t0

    t0 = time.perf_counter()
    warped_boundary = save_initial_registration_previews(
        save_path,
        volumereg,
        avreg,
        transform_matlab,
        boundary_atlas=boundary_reg,
        atlas_provider=atlas_display_provider_from_config(config.atlas),
    )
    preview_elapsed = time.perf_counter() - t0

    status, status_message, warnings = classify_init_registration_status(
        median_error_vox=metrics["median_vox"],
        auto_pairs=int(cpsample.shape[0]),
        inlier_fraction=metrics["inlier_fraction"],
        similarity_scale=scale,
        warped_boundary_voxels=int(warped_boundary),
        sample_cloud_points=int(ls_cloud.shape[0]),
        atlas_cloud_points=int(tv_cloud.shape[0]),
        alignment_backend=backend,
    )
    if atlas.boundary_path is None:
        warnings.append(
            "annotation_boundary_10.nii.gz not found in atlas_dir — "
            "boundaries derived from annotation labels."
        )

    diagnostics = InitRegistrationDiagnostics(
        sample_shape=[int(v) for v in volumereg.shape],
        atlas_shape=[int(v) for v in tvreg.shape],
        orientation=list(permvec),
        registration_resolution_um=float(checkpoint.registres_um),
        cloud_threshold=float(config.registration.cloud_threshold),
        sample_cloud_subsample=float(config.registration.sample_cloud_subsample),
        sample_cloud_points=int(ls_cloud.shape[0]),
        sample_mask_points=sample_stages.mask_points,
        sample_trim_points=sample_stages.trim_points,
        sample_downsample_points=sample_stages.downsample_points,
        sample_denoise_points=sample_stages.denoise_points,
        atlas_cloud_points=int(tv_cloud.shape[0]),
        alignment_backend=backend,
        similarity_scale=scale,
        median_error_sample_to_atlas_vox=metrics["median_sample_to_atlas_vox"],
        median_error_atlas_to_sample_vox=metrics["median_atlas_to_sample_vox"],
        median_error_vox=metrics["median_vox"],
        p95_error_sample_to_atlas_vox=metrics["p95_sample_to_atlas_vox"],
        p95_error_atlas_to_sample_vox=metrics["p95_atlas_to_sample_vox"],
        inlier_fraction=metrics["inlier_fraction"],
        inlier_threshold_vox=metrics["inlier_threshold_vox"],
        auto_pairs=int(cpsample.shape[0]),
        warped_boundary_voxels=int(warped_boundary),
        alignment_elapsed_s=alignment_elapsed,
        triage_elapsed_s=triage_elapsed,
        preview_elapsed_s=preview_elapsed,
        status=status,
        status_message=status_message,
        warnings=warnings,
    )
    diag_path = save_path / INIT_DIAGNOSTICS_FILENAME
    diagnostics.save(diag_path)
    diagnostics.print_summary(console=console)
    console.print(
        f"[dim]Checkpoint {regopts_path.name} · diagnostics {diag_path.name} · "
        f"cloud extraction {sample_cloud_elapsed + atlas_cloud_elapsed:.1f}s[/dim]"
    )

    checkpoint.permute_sample_to_atlas = permvec
    checkpoint.original_trans = transform_matlab.tolist()
    checkpoint.downfac_reg = downfac
    checkpoint.autocpsample = cpsample.tolist()
    checkpoint.autocpatlas = cpatlas.tolist()
    checkpoint.auto_points_source = "global_triage"
    checkpoint.auto_points_refined = False
    checkpoint.auto_points_mode = None
    checkpoint.auto_points_correspondence_path = None
    checkpoint.brain_atlas = config.atlas.provider.value
    if atlas_content.is_trimmed:
        checkpoint.atlas_crop_start_native = list(atlas_content.crop_start_yxz)
        checkpoint.atlas_native_shape = list(atlas_content.native_shape)
        if atlas_content.manifest is not None:
            save_atlas_manifest_copy(atlas_content.manifest, save_path)
    checkpoint.save(regopts_path)
    return checkpoint
