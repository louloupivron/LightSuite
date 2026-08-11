"""Inject MATLAB clouds into the Python init/affine path (decisive parity test)."""

from __future__ import annotations

import json
import time
from dataclasses import replace
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

from lightsuite.config.loader import load_config
from lightsuite.config.models import BrainPipelineConfig
from lightsuite.gui.affine import fit_affine_transform, transform_points
from lightsuite.gui.control_points import ControlPointSession
from lightsuite.preprocess.checkpoint import RegOptsCheckpoint
from lightsuite.registration.align import estimate_similarity_transform, triage_and_match_clouds
from lightsuite.registration.bcpd import find_bcpd_executable
from lightsuite.registration.brain_register import _prepare_control_points
from lightsuite.registration.pc_downsample import pcdenoise, pcdownsample_random
from lightsuite.registration.points import extract_sample_points
from lightsuite.registration.points_utils import thin_point_list
from lightsuite.registration.volume import load_registration_volume, permute_brain_volume
from lightsuite.registration.warp import swap_xy_transform

np.set_printoptions(precision=3, suppress=True)

CONFIG = Path("examples/config/mesoSPIM/marianna.yaml")
PY_SAVE = Path(
    "/media/gbm/NVME2/ALICe-pipelines-data/Marianna_CMU/registered_allen_MATLAB_conterpart"
)
ML_CLOUD = Path(
    "/run/user/1002/gvfs/smb-share:server=h1data1.wysscenter.ch,share=computingdata"
    "/Alice/0003_CBT_WYSS_LIGHTSHEET/DATA/MESOSPIM/Marianna CMU/elastix_temp/cloud_debug_matlab"
)
ML_PAIRS = Path(
    "/run/user/1002/gvfs/smb-share:server=h1data1.wysscenter.ch,share=computingdata"
    "/Alice/0003_CBT_WYSS_LIGHTSHEET/DATA/MESOSPIM/Marianna CMU/elastix_temp"
)

MATLAB_SAMP_TO_ATLAS = np.array(
    [
        [2.0360, 0.3755, 0.0471, -146.6393],
        [-0.4070, 1.7069, -0.1748, 92.9538],
        [-0.0463, 0.1459, 1.8934, -283.9107],
        [0.0, 0.0, 0.0, 1.0],
    ]
)


def _load_cloud(path: Path, *, one_based: bool = True) -> np.ndarray:
    pts = np.loadtxt(path, dtype=float)
    if pts.ndim == 1:
        pts = pts.reshape(1, -1)
    if one_based:
        pts = pts - 1.0
    return pts


def _affine_samp_to_atlas_xyz(tform_aff_yxz: np.ndarray) -> np.ndarray:
    """Stored matrix: sample (cloud XYZ @ 20 µm) -> atlas (cloud XYZ @ 10 µm)."""
    return swap_xy_transform(np.linalg.inv(tform_aff_yxz))


def _compare_affine(label: str, tform_aff_yxz: np.ndarray) -> None:
    got = _affine_samp_to_atlas_xyz(tform_aff_yxz)
    t_diff = got[:3, 3] - MATLAB_SAMP_TO_ATLAS[:3, 3]
    centre_yxz = np.array([[421.5, 175.5, 421.5]])
    centre_disp = np.linalg.norm(
        transform_points(centre_yxz, np.linalg.inv(got))
        - transform_points(centre_yxz, np.linalg.inv(MATLAB_SAMP_TO_ATLAS))
    )
    print(f"\n{label}")
    print(f"  translation diff (vox @ 20µm): {t_diff}  |diff|={np.linalg.norm(t_diff):.2f}")
    print(f"  centre displacement (vox):     {centre_disp:.2f}")
    print(f"  max |linear diff|:             {np.abs(got[:3, :3] - MATLAB_SAMP_TO_ATLAS[:3, :3]).max():.4f}")


def _run_affine_from_checkpoint(
    checkpoint: RegOptsCheckpoint,
    config: BrainPipelineConfig,
    original_trans: np.ndarray,
    *,
    label: str,
) -> np.ndarray:
    session = ControlPointSession.empty(original_trans, n_slices=3)
    tform_aff, _, _, _, diag = _prepare_control_points(checkpoint, session, config)
    print(
        f"{label}: n_auto={diag.n_auto}  median_err={diag.median_error_vox:.2f} vox  "
        f"coarse_auto={diag.median_coarse_auto_vox:.2f} vox"
        if diag.median_coarse_auto_vox is not None
        else f"{label}: n_auto={diag.n_auto}  median_err={diag.median_error_vox:.2f} vox"
    )
    _compare_affine(label, tform_aff)
    return tform_aff


def _pair_overlap(
    ref_sample: np.ndarray,
    ref_atlas: np.ndarray,
    query_sample: np.ndarray,
    query_atlas: np.ndarray,
    tol: float,
) -> tuple[int, int]:
    tree_s = cKDTree(ref_sample)
    ds, is_ = tree_s.query(query_sample, k=1)
    da, ia = cKDTree(ref_atlas).query(query_atlas, k=1)
    matched = (is_ == ia) & (ds <= tol) & (da <= tol)
    return int(matched.sum()), int(query_sample.shape[0])


def extract_sample_stages(volume: np.ndarray, threshold: float, *, subsample: float) -> dict[str, int]:
    """Mirror extractSamplePoints stages for count diagnostics."""
    from lightsuite.registration.points import (
        _batch_intensity_threshold,
        _imgradient3_magnitude,
        _trim_border_points,
        _volume_mode_all,
    )

    rng = np.random.default_rng(1)
    grad = _imgradient3_magnitude(volume)
    overall_mode = _volume_mode_all(volume)
    sizevol = volume.shape
    batch = max(1, int(np.ceil(max(sizevol) / 3)))
    nb = [max(1, int(np.ceil(s / batch))) for s in sizevol]
    chunks: list[np.ndarray] = []

    for iby in range(nb[0]):
        ys = slice(iby * batch, min((iby + 1) * batch, sizevol[0]))
        for ibx in range(nb[1]):
            xs = slice(ibx * batch, min((ibx + 1) * batch, sizevol[1]))
            for ibz in range(nb[2]):
                zs = slice(ibz * batch, min((ibz + 1) * batch, sizevol[2]))
                from scipy import ndimage

                vol_curr = ndimage.median_filter(volume[ys, xs, zs], size=1)
                grad_curr = grad[ys, xs, zs].copy()
                thres_init = _batch_intensity_threshold(vol_curr, overall_mode, rng)
                grad_curr[vol_curr < thres_init] = 0
                with np.errstate(divide="ignore", invalid="ignore"):
                    ratio = grad_curr / vol_curr
                mask = (ratio > threshold) & (vol_curr > thres_init)
                if not np.any(mask):
                    continue
                rr, cc, dd = np.where(mask)
                y_idx = np.arange(ys.start, ys.stop)[rr]
                x_idx = np.arange(xs.start, xs.stop)[cc]
                z_idx = np.arange(zs.start, zs.stop)[dd]
                chunks.append(np.column_stack([x_idx, y_idx, z_idx]))

    if not chunks:
        return {"mask": 0, "trimmed": 0, "downsampled": 0, "denoised": 0}

    pts = np.vstack(chunks)
    n_mask = pts.shape[0]
    pts = _trim_border_points(pts, sizevol)
    n_trim = pts.shape[0]
    pts_ds = pcdownsample_random(pts, subsample, preserve_structure=True, seed=1)
    pts_dn = pcdenoise(pts_ds)
    return {
        "mask": n_mask,
        "trimmed": n_trim,
        "downsampled": pts_ds.shape[0],
        "denoised": pts_dn.shape[0],
    }


def main() -> None:
    config = load_config(CONFIG)
    base = RegOptsCheckpoint.load(PY_SAVE / "regopts.json")
    downfac = float(base.downfac_reg or 0.5)
    bcpd_path = find_bcpd_executable(config.registration.bcpd_path)

    ml_ls = _load_cloud(ML_CLOUD / "ls_cloud.txt")
    ml_tv = _load_cloud(ML_CLOUD / "tv_cloud.txt")
    ml_auto_s = np.loadtxt(ML_PAIRS / "matlab_autocpsample.txt", delimiter=",") - 1.0
    ml_auto_a = np.loadtxt(ML_PAIRS / "matlab_autocpatlas.txt", delimiter=",") - 1.0

    print("=" * 72)
    print("0. CLOUD COUNTS")
    print("=" * 72)
    print(f"MATLAB ls_cloud (0-based): {ml_ls.shape[0]:,}")
    print(f"MATLAB tv_cloud (0-based): {ml_tv.shape[0]:,}")
    py_ls = extract_sample_points(
        permute_brain_volume(
            load_registration_volume(Path(base.regvolpath)).astype(np.float32),
            base.permute_sample_to_atlas or [1, 2, 3],
        ),
        config.registration.cloud_threshold,
        subsample_fraction=config.registration.sample_cloud_subsample,
    )
    print(f"Python extract_sample_points (fixed): {py_ls.shape[0]:,}")
    tree = cKDTree(ml_ls)
    ds, _ = tree.query(py_ls, k=1)
    print(
        f"Python ⊂ MATLAB spatially (≤5 vox): {(ds <= 5).sum():,}/{py_ls.shape[0]:,} "
        f"({100 * (ds <= 5).mean():.1f}%)"
    )

    stages = extract_sample_stages(
        permute_brain_volume(
            load_registration_volume(Path(base.regvolpath)).astype(np.float32),
            base.permute_sample_to_atlas or [1, 2, 3],
        ),
        config.registration.cloud_threshold,
        subsample=config.registration.sample_cloud_subsample,
    )
    print("\nPython extraction stages:")
    for key, val in stages.items():
        pct = 100 * val / stages["mask"] if stages["mask"] else 0
        print(f"  {key:12s}: {val:>9,}  ({pct:.1f}% of mask)")

    print("\n" + "=" * 72)
    print("1. INJECT MATLAB CLOUDS → Python BCPD + triage")
    print("=" * 72)
    t0 = time.perf_counter()
    transform_icp, transform_matlab, backend = estimate_similarity_transform(
        ml_tv, ml_ls, bcpd_path=bcpd_path
    )
    print(f"BCPD backend={backend}  elapsed={time.perf_counter() - t0:.1f}s")
    print(f"original_trans translation: {transform_matlab[:3, 3]}")
    py_orig = np.asarray(base.original_trans, dtype=float)
    print(f"Python stored translation:  {py_orig[:3, 3]}")
    print(f"|Δt| vs Python stored:       {np.linalg.norm(transform_matlab[:3, 3] - py_orig[:3, 3]):.2f} vox")

    t0 = time.perf_counter()
    inj_sample, inj_atlas = triage_and_match_clouds(ml_ls, ml_tv, transform_icp, bcpd_path=bcpd_path)
    print(
        f"triage elapsed={time.perf_counter() - t0:.1f}s  "
        f"pairs={inj_sample.shape[0]:,}  (MATLAB export {ml_auto_s.shape[0]:,})"
    )
    for tol in (2.0, 5.0):
        matched, n = _pair_overlap(ml_auto_s, ml_auto_a, inj_sample, inj_atlas, tol)
        print(f"  overlap with MATLAB auto pairs @ {tol} vox: {matched:,}/{n:,} ({100*matched/n:.1f}%)")

    inj_ckpt = replace(
        base,
        original_trans=transform_matlab.tolist(),
        autocpsample=inj_sample.tolist(),
        autocpatlas=inj_atlas.tolist(),
    )
    _run_affine_from_checkpoint(
        inj_ckpt,
        config,
        transform_matlab,
        label="A) MATLAB clouds → Python init pairs → affine",
    )

    print("\n" + "=" * 72)
    print("2. INJECT MATLAB AUTO PAIRS (bypass triage)")
    print("=" * 72)
    ml_ckpt = replace(
        base,
        original_trans=transform_matlab.tolist(),
        autocpsample=ml_auto_s.tolist(),
        autocpatlas=ml_auto_a.tolist(),
    )
    _run_affine_from_checkpoint(
        ml_ckpt,
        config,
        transform_matlab,
        label="B) MATLAB auto pairs + MATLAB-cloud BCPD original_trans",
    )

    print("\n" + "=" * 72)
    print("3. PYTHON STORED PAIRS (baseline)")
    print("=" * 72)
    _run_affine_from_checkpoint(
        base,
        config,
        py_orig,
        label="C) Python regopts pairs (current pipeline)",
    )

    print("\n" + "=" * 72)
    print("4. THINNED LANDMARK COUNTS")
    print("=" * 72)
    distancethin = 1000.0 / config.atlas.resolution_um
    from lightsuite.registration.points import cloud_xyz_to_volume_indices

    for label, sample_xyz, atlas_xyz in (
        ("MATLAB export", ml_auto_s, ml_auto_a),
        ("Python inject", inj_sample, inj_atlas),
        ("Python stored", np.asarray(base.autocpsample), np.asarray(base.autocpatlas)),
    ):
        atlas_native = cloud_xyz_to_volume_indices(atlas_xyz) / downfac
        n_thin = int(thin_point_list(atlas_native, distancethin).sum())
        print(f"{label:16s}: {sample_xyz.shape[0]:,} pairs → {n_thin} thinned landmarks")

    stored = json.loads((PY_SAVE / "transform_params.json").read_text())
    py_stored = np.asarray(stored["tform_affine_samp20um_to_atlas_10um_px"])
    print("\nStored Python matrix centre displacement vs MATLAB:")
    centre_yxz = np.array([[421.5, 175.5, 421.5]])
    centre_disp = np.linalg.norm(
        transform_points(centre_yxz, np.linalg.inv(py_stored))
        - transform_points(centre_yxz, np.linalg.inv(MATLAB_SAMP_TO_ATLAS))
    )
    print(f"  {centre_disp:.2f} vox  (translation |diff|={np.linalg.norm(py_stored[:3,3]-MATLAB_SAMP_TO_ATLAS[:3,3]):.2f})")


if __name__ == "__main__":
    main()
