"""Compare MATLAB vs Python init-registration auto control point pairs."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

from lightsuite.gui.affine import fit_affine_transform, transform_points
from lightsuite.registration.points import cloud_xyz_to_volume_indices
from lightsuite.registration.points_utils import thin_point_list
from lightsuite.registration.warp import swap_xy_transform

np.set_printoptions(precision=3, suppress=True)

ML_DIR = Path(
    "/run/user/1002/gvfs/smb-share:server=h1data1.wysscenter.ch,share=computingdata"
    "/Alice/0003_CBT_WYSS_LIGHTSHEET/DATA/MESOSPIM/Marianna CMU/elastix_temp"
)
PY_DIR = Path(
    "/media/gbm/NVME2/ALICe-pipelines-data/Marianna_CMU/registered_allen_MATLAB_conterpart"
)

ml_sample = np.loadtxt(ML_DIR / "matlab_autocpsample.txt", delimiter=",")
ml_atlas = np.loadtxt(ML_DIR / "matlab_autocpatlas.txt", delimiter=",")
ro = json.loads((PY_DIR / "regopts.json").read_text())
py_sample = np.asarray(ro["autocpsample"], dtype=float)
py_atlas = np.asarray(ro["autocpatlas"], dtype=float)
downfac = float(ro["downfac_reg"])

print("=" * 72)
print("1. PAIR COUNTS")
print("=" * 72)
print(f"MATLAB export:  sample {ml_sample.shape[0]:,}  atlas {ml_atlas.shape[0]:,}")
print(f"Python regopts: sample {py_sample.shape[0]:,}  atlas {py_atlas.shape[0]:,}")
print(f"Gap: {ml_sample.shape[0] - py_sample.shape[0]:,} fewer in Python ({100*(ml_sample.shape[0]-py_sample.shape[0])/ml_sample.shape[0]:.1f}%)")

print("\n" + "=" * 72)
print("2. COORDINATE RANGES (cloud X,Y,Z @ 20 µm reg grid)")
print("=" * 72)
for name, pts in (("MATLAB sample", ml_sample), ("Python sample", py_sample),
                  ("MATLAB atlas", ml_atlas), ("Python atlas", py_atlas)):
    lo = pts.min(axis=0)
    hi = pts.max(axis=0)
    print(f"{name:14s}  X [{lo[0]:6.0f},{hi[0]:6.0f}]  Y [{lo[1]:6.0f},{hi[1]:6.0f}]  Z [{lo[2]:6.0f},{hi[2]:6.0f}]")


def match_pairs(
    ref_sample: np.ndarray,
    ref_atlas: np.ndarray,
    query_sample: np.ndarray,
    query_atlas: np.ndarray,
    tol: float,
) -> dict[str, float | int]:
    """For each query pair, find nearest ref pair by sample+atlas combined distance."""
    tree_s = cKDTree(ref_sample)
    tree_a = cKDTree(ref_atlas)
    ds, is_ = tree_s.query(query_sample, k=1)
    da, ia = tree_a.query(query_atlas, k=1)
    # same-index correspondence in ref
    same_idx = is_ == ia
    close_s = ds <= tol
    close_a = da <= tol
    matched = same_idx & close_s & close_a
  # also count "shifted": close in sample but atlas index differs or atlas far
    shifted = close_s & (~same_idx | ~close_a)
    missing = ~close_s
    return {
        "matched": int(matched.sum()),
        "shifted_sample_close": int(shifted.sum()),
        "missing": int(missing.sum()),
        "n_query": int(query_sample.shape[0]),
    }


print("\n" + "=" * 72)
print("3. PAIR OVERLAP (Python pairs vs MATLAB pool)")
print("=" * 72)
for tol in (1.0, 2.0, 5.0, 10.0):
    m = match_pairs(ml_sample, ml_atlas, py_sample, py_atlas, tol)
    print(
        f"tol={tol:4.1f} vox:  exact matches {m['matched']:5,}/{m['n_query']:,} "
        f"({100*m['matched']/m['n_query']:.1f}%)  "
        f"shifted {m['shifted_sample_close']:,}  missing {m['missing']:,}"
    )

print("\n" + "=" * 72)
print("4. REVERSE: MATLAB pairs not in Python pool")
print("=" * 72)
for tol in (2.0, 5.0):
    m = match_pairs(py_sample, py_atlas, ml_sample, ml_atlas, tol)
    print(
        f"tol={tol:4.1f} vox:  MATLAB pairs with Python counterpart "
        f"{m['matched']:5,}/{m['n_query']:,} ({100*m['matched']/m['n_query']:.1f}%)  "
        f"MATLAB-only ~{m['missing']:,}"
    )

print("\n" + "=" * 72)
print("5. SAMPLE-ONLY NN (ignoring atlas index coupling)")
print("=" * 72)
tree = cKDTree(ml_sample)
ds, _ = tree.query(py_sample, k=1)
for tol in (1.0, 2.0, 5.0):
    print(f"  Python sample points with MATLAB sample within {tol} vox: {(ds<=tol).sum():,}/{len(ds):,} ({100*(ds<=tol).mean():.1f}%)")
print(f"  median NN sample distance: {np.median(ds):.2f} vox  p95: {np.percentile(ds,95):.2f}")

da_tree = cKDTree(ml_atlas)
da, _ = da_tree.query(py_atlas, k=1)
print(f"  median NN atlas distance:  {np.median(da):.2f} vox  p95: {np.percentile(da,95):.2f}")

print("\n" + "=" * 72)
print("6. THINNED SUBSETS (distancethin=100, affine-fit landmarks)")
print("=" * 72)
distancethin = 100.0
ml_a_native = ml_atlas / downfac
py_a_native = cloud_xyz_to_volume_indices(py_atlas) / downfac
py_s_yxz = cloud_xyz_to_volume_indices(py_sample)

keep_ml = thin_point_list(ml_a_native, distancethin)
keep_py = thin_point_list(py_a_native, distancethin)
print(f"MATLAB thinned (XYZ atlas native): {keep_ml.sum():,}")
print(f"Python thinned (YXZ atlas native): {keep_py.sum():,}")

ml_thin_s = ml_sample[keep_ml]
ml_thin_a = ml_a_native[keep_ml]
py_thin_s = py_s_yxz[keep_py]
py_thin_a = py_a_native[keep_py]

# compare thinned sample locations
ds_thin, _ = cKDTree(ml_thin_s).query(py_thin_s, k=1)
print(f"Thinned sample NN: median {np.median(ds_thin):.1f} vox  p95 {np.percentile(ds_thin,95):.1f}")
print(f"Thinned points within 20 vox: {(ds_thin<=20).sum()}/{len(ds_thin)} ({100*(ds_thin<=20).mean():.0f}%)")

print("\n" + "=" * 72)
print("7. AFFINE SENSITIVITY: fit on MATLAB pairs vs Python pairs")
print("=" * 72)
# MATLAB convention: cloud XYZ, atlas native
t_ml, _ = fit_affine_transform(ml_a_native[keep_ml], ml_sample[keep_ml])
# Python convention: YXZ
t_py, _ = fit_affine_transform(py_a_native[keep_py], py_s_yxz[keep_py])

c = np.array([[421.5, 175.5, 421.5]])  # sample centre YXZ
# convert MATLAB fit to YXZ for comparison
t_ml_yxz = swap_xy_transform(np.linalg.inv(swap_xy_transform(np.linalg.inv(t_ml))))
# simpler: map centre through each forward atlas->sample
centre_disp = np.linalg.norm(transform_points(c, t_py) - transform_points(c, t_ml_yxz))
print(f"Centre displacement (20µm vox) MATLAB-pairs fit vs Python-pairs fit: {centre_disp:.1f}")

# stored transforms
tp = json.loads((PY_DIR / "transform_params.json").read_text())
stored_fwd = np.array(tp["tform_affine_atlas_to_samp20um_px"])
matlab_s2a = np.array(
    [[2.036, 0.3755, 0.0471, -146.6393], [-0.407, 1.7069, -0.1748, 92.9538],
     [-0.0463, 0.1459, 1.8934, -283.9107], [0, 0, 0, 1]]
)
matlab_fwd = np.linalg.inv(swap_xy_transform(matlab_s2a))

for name, t in (("Python stored", stored_fwd), ("Python regopts fit", t_py),
                ("MATLAB pairs refit", t_ml_yxz), ("MATLAB stored", matlab_fwd)):
    print(f"  {name:22s} centre -> {transform_points(c, t)}")
