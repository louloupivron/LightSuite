"""Diagnose Python vs MATLAB affine-step divergence (scratch)."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from lightsuite.gui.affine import fit_affine_transform, transform_points
from lightsuite.registration.points import cloud_xyz_to_volume_indices
from lightsuite.registration.warp import swap_xy_transform

np.set_printoptions(precision=4, suppress=True)

SAVE = Path(
    "/media/gbm/NVME2/ALICe-pipelines-data/Marianna_CMU/registered_allen_MATLAB_conterpart"
)

MATLAB_SAMP_TO_ATLAS = np.array(
    [
        [2.0360, 0.3755, 0.0471, -146.6393],
        [-0.4070, 1.7069, -0.1748, 92.9538],
        [-0.0463, 0.1459, 1.8934, -283.9107],
        [0.0, 0.0, 0.0, 1.0],
    ]
)

regopts = json.loads((SAVE / "regopts.json").read_text())
tparams = json.loads((SAVE / "transform_params.json").read_text())
python_samp_to_atlas = np.array(tparams["tform_affine_samp20um_to_atlas_10um_px"])

downfac = float(regopts["downfac_reg"])
autocpsample_xyz = np.asarray(regopts["autocpsample"], dtype=float)
autocpatlas_xyz = np.asarray(regopts["autocpatlas"], dtype=float)
print(f"auto pairs in regopts: {autocpsample_xyz.shape[0]:,}   downfac={downfac}")


def decompose(name: str, M: np.ndarray) -> None:
    L = M[:3, :3]
    det = np.linalg.det(L)
    u, s, vt = np.linalg.svd(L)
    print(f"{name}")
    print(f"   det={det:.4f}  singular values={s}  aniso={s.max() / s.min():.4f}")
    print(f"   translation={M[:3, 3]}")


print("\n" + "=" * 72)
print("1. STORED MATRICES AS-IS")
print("=" * 72)
decompose("MATLAB tform_affine_samp20um_to_atlas_10um_px", MATLAB_SAMP_TO_ATLAS)
decompose("Python tform_affine_samp20um_to_atlas_10um_px", python_samp_to_atlas)
print(f"   max |linear diff| = {np.abs(MATLAB_SAMP_TO_ATLAS[:3, :3] - python_samp_to_atlas[:3, :3]).max():.4f}")

print("\n" + "=" * 72)
print("2. AFTER SWAPPING MATLAB XYZ -> ARRAY (Y,X,Z)")
print("=" * 72)
matlab_swapped = swap_xy_transform(MATLAB_SAMP_TO_ATLAS)
print("MATLAB swapped to array order:\n", matlab_swapped)
print("\nPython stored:\n", python_samp_to_atlas)
print(f"\n   max |linear diff| = {np.abs(matlab_swapped[:3, :3] - python_samp_to_atlas[:3, :3]).max():.4f}")
print(f"   translation diff  = {python_samp_to_atlas[:3, 3] - matlab_swapped[:3, 3]}")
print(f"   translation |diff| = {np.linalg.norm(python_samp_to_atlas[:3, 3] - matlab_swapped[:3, 3]):.2f}")

print("\n" + "=" * 72)
print("3. REFIT FROM PYTHON CONTROL POINTS IN BOTH AXIS ORDERS")
print("=" * 72)
# Python convention: array (Y,X,Z)
sample_yxz = cloud_xyz_to_volume_indices(autocpsample_xyz)
atlas_yxz = cloud_xyz_to_volume_indices(autocpatlas_xyz) / downfac
# MATLAB convention: cloud (X,Y,Z), atlas divided by downfac only
sample_xyz = autocpsample_xyz
atlas_xyz = autocpatlas_xyz / downfac

fit_yxz, mse_yxz = fit_affine_transform(atlas_yxz, sample_yxz)
fit_xyz, mse_xyz = fit_affine_transform(atlas_xyz, sample_xyz)
inv_yxz = np.linalg.inv(fit_yxz)
inv_xyz = np.linalg.inv(fit_xyz)

print(f"fit in array (Y,X,Z)  mse={mse_yxz:.2f}")
print(f"fit in cloud (X,Y,Z)  mse={mse_xyz:.2f}")
print("\ninv(fit) in cloud XYZ order (directly comparable to MATLAB stored):")
print(inv_xyz)
print("\nMATLAB stored:")
print(MATLAB_SAMP_TO_ATLAS)
print(f"\n   max |linear diff| = {np.abs(inv_xyz[:3, :3] - MATLAB_SAMP_TO_ATLAS[:3, :3]).max():.4f}")
print(f"   translation diff   = {inv_xyz[:3, 3] - MATLAB_SAMP_TO_ATLAS[:3, 3]}")

print("\n   consistency check: swap(fit_yxz) vs fit_xyz")
print(f"   max diff = {np.abs(swap_xy_transform(inv_yxz) - inv_xyz).max():.6f}")

print("\n" + "=" * 72)
print("4. WHERE THE 241 THINNED POINTS LAND")
print("=" * 72)
for label, atlas_pts, sample_pts, M in (
    ("python YXZ fit", atlas_yxz, sample_yxz, fit_yxz),
    ("cloud XYZ fit", atlas_xyz, sample_xyz, fit_xyz),
):
    pred = transform_points(atlas_pts, M)
    err = np.linalg.norm(pred - sample_pts, axis=1)
    print(f"{label}: median={np.median(err):.2f} p95={np.percentile(err, 95):.2f} vox (all {len(err):,} pairs)")

print("\n" + "=" * 72)
print("5. DIAGNOSTIC BUG CHECK (median_coarse_auto_vox)")
print("=" * 72)
original_trans_vol = swap_xy_transform(np.asarray(regopts["original_trans"], dtype=float))
coarse = transform_points(sample_yxz, original_trans_vol)
print(f"coarse (reg-grid 20um atlas) vs af_atlas (native 10um):")
print(f"   median = {np.median(np.linalg.norm(coarse - atlas_yxz, axis=1)):.2f} vox   <- current stat")
print(f"   with /downfac applied to coarse:")
print(f"   median = {np.median(np.linalg.norm(coarse / downfac - atlas_yxz, axis=1)):.2f} vox")
