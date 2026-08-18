"""Diagnose sample-cloud mask-stage gap vs MATLAB extractSamplePoints."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy import ndimage
from scipy.spatial import cKDTree

from lightsuite.registration.points import (
    _batch_intensity_threshold,
    _imgradient3_magnitude,
    _trim_border_points,
    _volume_mode_all,
)
from lightsuite.registration.volume import load_registration_volume, permute_brain_volume

np.set_printoptions(precision=4, suppress=True)

PY_SAVE = Path(
    "/media/gbm/NVME2/ALICe-pipelines-data/Marianna_CMU/registered_allen_MATLAB_conterpart"
)
ML_CLOUD = Path(
    "/run/user/1002/gvfs/smb-share:server=h1data1.wysscenter.ch,share=computingdata"
    "/Alice/0003_CBT_WYSS_LIGHTSHEET/DATA/MESOSPIM/Marianna CMU/elastix_temp/cloud_debug_matlab"
)
THRESH = 5.0


def build_mask_points(
    volume: np.ndarray,
    threshold: float,
    *,
    per_batch_threshold: bool = True,
    global_thres_init: float | None = None,
) -> tuple[np.ndarray, dict[str, int]]:
    """Return pre-trim mask voxel coords (cloud XYZ) and stage counts."""
    rng = np.random.default_rng(1)
    grad = _imgradient3_magnitude(volume)
    overall_mode = _volume_mode_all(volume)
    if global_thres_init is None:
        flat = volume.ravel(order="F")
        sample_n = min(10_000, flat.size)
        idx = rng.choice(flat.size, size=sample_n, replace=False)
        global_thres_init = max(float(np.quantile(flat[idx], 0.05)) * 2.0, overall_mode)

    sizevol = volume.shape
    batch = max(1, int(np.ceil(max(sizevol) / 3)))
    nb = [max(1, int(np.ceil(s / batch))) for s in sizevol]
    chunks: list[np.ndarray] = []
    n_zeroed_by_thres = 0
    n_ratio_pass = 0
    n_vol_pass = 0

    for iby in range(nb[0]):
        ys = slice(iby * batch, min((iby + 1) * batch, sizevol[0]))
        for ibx in range(nb[1]):
            xs = slice(ibx * batch, min((ibx + 1) * batch, sizevol[1]))
            for ibz in range(nb[2]):
                zs = slice(ibz * batch, min((ibz + 1) * batch, sizevol[2]))
                vol_curr = ndimage.median_filter(volume[ys, xs, zs], size=1)
                grad_curr = grad[ys, xs, zs].copy()
                thres_init = (
                    _batch_intensity_threshold(vol_curr, overall_mode, rng)
                    if per_batch_threshold
                    else global_thres_init
                )
                low_vol = vol_curr < thres_init
                n_zeroed_by_thres += int(np.count_nonzero(low_vol))
                grad_curr[low_vol] = 0
                with np.errstate(divide="ignore", invalid="ignore"):
                    ratio = grad_curr / vol_curr
                vol_ok = vol_curr > thres_init
                ratio_ok = ratio > threshold
                n_vol_pass += int(np.count_nonzero(vol_ok))
                n_ratio_pass += int(np.count_nonzero(ratio_ok))
                mask = ratio_ok & vol_ok
                if not np.any(mask):
                    continue
                rr, cc, dd = np.where(mask)
                y_idx = np.arange(ys.start, ys.stop)[rr]
                x_idx = np.arange(xs.start, xs.stop)[cc]
                z_idx = np.arange(zs.start, zs.stop)[dd]
                chunks.append(np.column_stack([x_idx, y_idx, z_idx]))

    if not chunks:
        pts = np.zeros((0, 3), dtype=np.float64)
    else:
        pts = np.vstack(chunks)
    trimmed = _trim_border_points(pts, sizevol)
    stats = {
        "mask": pts.shape[0],
        "trimmed": trimmed.shape[0],
        "n_voxels": int(volume.size),
        "overall_mode": int(overall_mode),
        "global_thres_init": float(global_thres_init),
    }
    return pts, trimmed, stats


def voxel_values(
    volume: np.ndarray,
    grad: np.ndarray,
    points_xyz: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Sample vol and grad at integer cloud-XYZ voxel indices (0-based)."""
    pts = np.rint(points_xyz).astype(int)
    y, x, z = pts[:, 1], pts[:, 0], pts[:, 2]
    return volume[y, x, z], grad[y, x, z]


def batch_thres_init_for_points(
    volume: np.ndarray,
    points_xyz: np.ndarray,
    *,
    per_batch: bool = True,
) -> np.ndarray:
    """Per-point thres_init using the same batch tiling as extractSamplePoints."""
    rng = np.random.default_rng(1)
    overall_mode = _volume_mode_all(volume)
    sizevol = volume.shape
    batch = max(1, int(np.ceil(max(sizevol) / 3)))
    nb = [max(1, int(np.ceil(s / batch))) for s in sizevol]

    # Precompute per-batch thresholds in visit order.
    batch_thres: list[float] = []
    for iby in range(nb[0]):
        ys = slice(iby * batch, min((iby + 1) * batch, sizevol[0]))
        for ibx in range(nb[1]):
            xs = slice(ibx * batch, min((ibx + 1) * batch, sizevol[1]))
            for ibz in range(nb[2]):
                zs = slice(ibz * batch, min((ibz + 1) * batch, sizevol[2]))
                vol_curr = ndimage.median_filter(volume[ys, xs, zs], size=1)
                if per_batch:
                    t = _batch_intensity_threshold(vol_curr, overall_mode, rng)
                else:
                    flat = volume.ravel(order="F")
                    sample_n = min(10_000, flat.size)
                    idx = np.random.default_rng(1).choice(flat.size, size=sample_n, replace=False)
                    t = max(float(np.quantile(flat[idx], 0.05)) * 2.0, overall_mode)
                batch_thres.append(t)

    pts = np.rint(points_xyz).astype(int)
    out = np.empty(pts.shape[0], dtype=float)
    bi = 0
    for iby in range(nb[0]):
        y0, y1 = iby * batch, min((iby + 1) * batch, sizevol[0])
        for ibx in range(nb[1]):
            x0, x1 = ibx * batch, min((ibx + 1) * batch, sizevol[1])
            for ibz in range(nb[2]):
                z0, z1 = ibz * batch, min((ibz + 1) * batch, sizevol[2])
                in_batch = (
                    (pts[:, 1] >= y0)
                    & (pts[:, 1] < y1)
                    & (pts[:, 0] >= x0)
                    & (pts[:, 0] < x1)
                    & (pts[:, 2] >= z0)
                    & (pts[:, 2] < z1)
                )
                out[in_batch] = batch_thres[bi]
                bi += 1
    return out


def diagnose_matlab_points(
    volume: np.ndarray,
    grad: np.ndarray,
    ml_points: np.ndarray,
) -> None:
    vol, g = voxel_values(volume, grad, ml_points)
    thres = batch_thres_init_for_points(volume, ml_points, per_batch=True)
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = g / vol
    vol_ok = vol > thres
    ratio_ok = ratio > THRESH
    both_ok = vol_ok & ratio_ok
    print("\nMATLAB final-cloud points vs Python mask criteria:")
    print(f"  n points: {ml_points.shape[0]:,}")
    print(f"  vol > thres_init:          {vol_ok.sum():,} ({100*vol_ok.mean():.1f}%)")
    print(f"  grad/vol > {THRESH:g}:           {ratio_ok.sum():,} ({100*ratio_ok.mean():.1f}%)")
    print(f"  both (Python mask):        {both_ok.sum():,} ({100*both_ok.mean():.1f}%)")
    fail_vol = ~vol_ok & ratio_ok
    fail_ratio = vol_ok & ~ratio_ok
    fail_both = ~vol_ok & ~ratio_ok
    print(f"  fail vol only:             {fail_vol.sum():,}")
    print(f"  fail ratio only:           {fail_ratio.sum():,}")
    print(f"  fail both:                 {fail_both.sum():,}")
    print(f"  ratio  median pass={np.median(ratio[both_ok]):.2f}  "
          f"fail_ratio={np.median(ratio[vol_ok & ~ratio_ok]) if fail_ratio.any() else 0:.2f}")
    border = _trim_border_points(ml_points, volume.shape)
    in_border = ml_points.shape[0] - border.shape[0]
    print(f"  border-trim rejected:      {in_border:,} ({100*in_border/ml_points.shape[0]:.1f}%)")


def main() -> None:
    regopts = json.loads((PY_SAVE / "regopts.json").read_text())
    volume = permute_brain_volume(
        load_registration_volume(Path(regopts["regvolpath"])).astype(np.float32),
        regopts.get("permute_sample_to_atlas") or [1, 2, 3],
    )
    grad = _imgradient3_magnitude(volume)
    ml_final = np.loadtxt(ML_CLOUD / "ls_cloud.txt") - 1.0

    print("=" * 72)
    print("1. MASK COUNTS (pre-downsample)")
    print("=" * 72)
    mask_pts, trim_pts, stats = build_mask_points(volume, THRESH, per_batch_threshold=True)
    print(f"volume shape: {volume.shape}  voxels={stats['n_voxels']:,}")
    print(f"overall mode: {stats['overall_mode']}  global thres_init≈{stats['global_thres_init']:.2f}")
    print(f"per-batch mask:  {stats['mask']:,}")
    print(f"after trim:      {stats['trimmed']:,}")
    print(f"MATLAB implied pre-downsample (~524k/0.1): ~5,239,930")

    mask_global, trim_global, stats_g = build_mask_points(
        volume, THRESH, per_batch_threshold=False, global_thres_init=stats["global_thres_init"]
    )
    print(f"\nwith GLOBAL thres_init only: mask={stats_g['mask']:,}  trim={stats_g['trimmed']:,}")

    print("\n" + "=" * 72)
    print("2. MATLAB FINAL CLOUD vs PYTHON MASK")
    print("=" * 72)
    tree = cKDTree(mask_pts)
    ds, _ = tree.query(ml_final, k=1)
    print(f"MATLAB final points within 1 vox of Python mask: {(ds<=1).sum():,}/{ml_final.shape[0]:,} "
          f"({100*(ds<=1).mean():.1f}%)")
    print(f"within 5 vox: {(ds<=5).sum():,} ({100*(ds<=5).mean():.1f}%)")
    diagnose_matlab_points(volume, grad, ml_final)

    # Points MATLAB has that Python mask misses entirely (>1 vox)
    miss = ml_final[ds > 1]
    if miss.shape[0]:
        diagnose_matlab_points(volume, grad, miss)

    print("\n" + "=" * 72)
    print("3. GRADIENT IMPLEMENTATION SENSITIVITY")
    print("=" * 72)
    # Central differences magnitude (alternative)
    gx = np.zeros_like(volume)
    gy = np.zeros_like(volume)
    gz = np.zeros_like(volume)
    gx[:, 1:-1, :] = (volume[:, 2:, :] - volume[:, :-2, :]) / 2.0
    gy[1:-1, :, :] = (volume[2:, :, :] - volume[:-2, :, :]) / 2.0
    gz[:, :, 1:-1] = (volume[:, :, 2:] - volume[:, :, :-2]) / 2.0
    grad_cd = np.sqrt(gx**2 + gy**2 + gz**2)
    for name, g in (("sobel (current)", grad), ("central-diff", grad_cd)):
        vol_s, g_s = voxel_values(volume, g, ml_final)
        thres = batch_thres_init_for_points(volume, ml_final, per_batch=True)
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = g_s / vol_s
        both = (vol_s > thres) & (ratio > THRESH)
        print(f"  {name:18s}: {both.sum():,}/{ml_final.shape[0]:,} pass ({100*both.mean():.1f}%)")

    print("\n" + "=" * 72)
    print("4. RATIO THRESHOLD SWEEP ON MATLAB POINTS")
    print("=" * 72)
    vol_s, g_s = voxel_values(volume, grad, ml_final)
    thres = batch_thres_init_for_points(volume, ml_final, per_batch=True)
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = g_s / vol_s
    vol_ok = vol_s > thres
    for t in (3.0, 4.0, 4.5, 5.0, 5.5, 6.0):
        n = int((vol_ok & (ratio > t)).sum())
        print(f"  threshold={t:g}: {n:,} ({100*n/ml_final.shape[0]:.1f}%)")

    print("\n" + "=" * 72)
    print("5. PYTHON-ONLY MASK POINTS (not near MATLAB final cloud)")
    print("=" * 72)
    tree_ml = cKDTree(ml_final)
    ds_py, _ = tree_ml.query(mask_pts, k=1)
    py_only = mask_pts[ds_py > 5]
    print(f"Python mask points >5 vox from any MATLAB final point: {py_only.shape[0]:,} "
          f"({100*py_only.shape[0]/mask_pts.shape[0]:.1f}% of mask)")
    if py_only.shape[0]:
        vol_p, g_p = voxel_values(volume, grad, py_only)
        th_p = batch_thres_init_for_points(volume, py_only, per_batch=True)
        with np.errstate(divide="ignore", invalid="ignore"):
            r_p = g_p / vol_p
        print(f"  median ratio={np.median(r_p):.2f}  median vol={np.median(vol_p):.2f}  "
              f"median thres={np.median(th_p):.2f}")


if __name__ == "__main__":
    main()
