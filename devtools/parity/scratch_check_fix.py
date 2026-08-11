"""Diagnose transform chain after parity fixes."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from lightsuite.registration.align import _estimate_similarity_bcpd, downsample_for_bcpd_similarity
from lightsuite.registration.bcpd import atlas_to_sample_affinetform, find_bcpd_executable, register_bcpd
from lightsuite.registration.warp import affinetform_rows_to_internal, matlab_voxel_affine_from_icp

np.set_printoptions(precision=4, suppress=True)

T_MATLAB = np.array(
    [
        [0.9271, 0.1870, 0.0514, -71.9302],
        [-0.1825, 0.9261, -0.0782, 4.1789],
        [-0.0657, 0.0667, 0.9425, -130.3184],
        [0.0, 0.0, 0.0, 1.0],
    ]
)

DEBUG = Path(
    "/run/user/1002/gvfs/smb-share:server=h1data1.wysscenter.ch,share=computingdata"
    "/Alice/0003_CBT_WYSS_LIGHTSHEET/DATA/MESOSPIM/Marianna CMU/elastix_temp/cloud_debug_matlab"
)

bcpd = find_bcpd_executable(None)
print("bcpd:", bcpd)

# Full chain via align helper on MATLAB clouds (1-based coords)
lspoints = np.loadtxt(DEBUG / "lspoints.txt")
tvpoints = np.loadtxt(DEBUG / "tvpoints.txt")
ls_cloud = np.loadtxt(DEBUG / "ls_cloud.txt")
tv_cloud = np.loadtxt(DEBUG / "tv_cloud.txt")

# MATLAB clouds are 1-based; convert to 0-based for Python pipeline
ls0 = ls_cloud - 1.0
tv0 = tv_cloud - 1.0

for label, sample, atlas in (
    ("MATLAB clouds 1-based via _estimate_similarity_bcpd", ls_cloud, tv_cloud),
    ("MATLAB clouds 0-based via _estimate_similarity_bcpd", ls0, tv0),
):
    t_icp, t_mat = _estimate_similarity_bcpd(atlas, sample, bcpd)
    dt = t_mat[:3, 3] - T_MATLAB[:3, 3]
    print(f"\n{label}")
    print("translation:", t_mat[:3, 3], " |dt|=", np.linalg.norm(dt))

# Downsample counts on Python-scale clouds
for name, cloud, div in (("sample", ls0, 10_000), ("atlas", tv0, 50_000)):
    sub = downsample_for_bcpd_similarity(cloud, div)
    print(f"{name}: full={cloud.shape[0]:,} bcpd={sub.shape[0]:,}")
