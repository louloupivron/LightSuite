"""Run Python BCPD on MATLAB's exact exported clouds and compare transforms."""

from __future__ import annotations

from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

from lightsuite.registration.bcpd import find_bcpd_executable, register_bcpd

np.set_printoptions(precision=4, suppress=True)

DEBUG = Path(
    "/run/user/1002/gvfs/smb-share:server=h1data1.wysscenter.ch,share=computingdata"
    "/Alice/0003_CBT_WYSS_LIGHTSHEET/DATA/MESOSPIM/Marianna CMU/elastix_temp/cloud_debug_matlab"
)

T_MATLAB = np.array(
    [
        [0.9271, 0.1870, 0.0514, -71.9302],
        [-0.1825, 0.9261, -0.0782, 4.1789],
        [-0.0657, 0.0667, 0.9425, -130.3184],
        [0.0, 0.0, 0.0, 1.0],
    ]
)

lspoints = np.loadtxt(DEBUG / "lspoints.txt")
tvpoints = np.loadtxt(DEBUG / "tvpoints.txt")
print(f"MATLAB lspoints {lspoints.shape}  tvpoints {tvpoints.shape}")

bcpd = find_bcpd_executable(None)
print(f"bcpd: {bcpd}")

kw = dict(
    bcpd_path=bcpd,
    outlier_ratio=0.01,
    gamma=1.0,
    beta=2.0,
    lambda_=50.0,
    convergence_tolerance=1e-8,
    normalize_common=False,
)

registered, _ = register_bcpd(tvpoints, lspoints, "similarity", **kw)
d_to_sample, _ = cKDTree(lspoints).query(registered, k=1)
d_to_reg, _ = cKDTree(registered).query(lspoints, k=1)
keep_atlas = d_to_sample <= 25.0
keep_sample = d_to_reg <= 25.0
print(f"pass-2 kept atlas {keep_atlas.sum():,}/{len(keep_atlas):,}  "
      f"sample {keep_sample.sum():,}/{len(keep_sample):,}")

_, aff = register_bcpd(tvpoints[keep_atlas], lspoints[keep_sample], "similarity", **kw)
# aff is the hybrid row-convention matrix used inside align.py:
#   forward (atlas -> sample) is  p @ aff[:3,:3] + aff[:3,3]
L = aff[:3, :3]
t = aff[:3, 3]

# --- current align.py behaviour -------------------------------------------
cur = np.linalg.inv(aff)
current = np.eye(4)
current[:3, :3] = cur[:3, :3].T
current[:3, 3] = cur[:3, 3]

# --- MATLAB behaviour: build the proper premultiply matrix, then invert ----
premul = np.eye(4)
premul[:3, :3] = L.T
premul[:3, 3] = t
fixed = np.linalg.inv(premul)

for name, T in (("current align.py", current), ("MATLAB-convention", fixed)):
    dt = T[:3, 3] - T_MATLAB[:3, 3]
    dl = np.abs(T[:3, :3] - T_MATLAB[:3, :3]).max()
    print(f"\n{name}")
    print(T)
    print(f"   translation {T[:3, 3]}  vs MATLAB {T_MATLAB[:3, 3]}")
    print(f"   |dt| = {np.linalg.norm(dt):.2f} vox   max|dLinear| = {dl:.4f}")
