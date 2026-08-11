"""Compare Python downsampling against MATLAB nonuniformGridSample on the same clouds."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from lightsuite.registration.align import _matlab_cloud_subset

np.set_printoptions(precision=2, suppress=True)

DEBUG = Path(
    "/run/user/1002/gvfs/smb-share:server=h1data1.wysscenter.ch,share=computingdata"
    "/Alice/0003_CBT_WYSS_LIGHTSHEET/DATA/MESOSPIM/Marianna CMU/elastix_temp/cloud_debug_matlab"
)


def kdtree_leaf_sample(points: np.ndarray, max_num_points: int, seed: int = 1) -> np.ndarray:
    """Median-split kd-tree to 2^k leaves, one representative point per leaf."""
    n = points.shape[0]
    levels = max(0, int(np.ceil(np.log2(max(1.0, n / max_num_points)))))
    rng = np.random.default_rng(seed)
    groups = [np.arange(n)]
    for _ in range(levels):
        nxt = []
        for idx in groups:
            if idx.size <= 1:
                nxt.append(idx)
                continue
            pts = points[idx]
            axis = int(np.argmax(pts.max(axis=0) - pts.min(axis=0)))
            order = np.argsort(pts[:, axis], kind="stable")
            half = idx.size // 2
            nxt.append(idx[order[:half]])
            nxt.append(idx[order[half:]])
        groups = nxt
    out = [points[idx[rng.integers(idx.size)]] for idx in groups if idx.size]
    return np.vstack(out)


for label, full_file, sub_file, divisor in (
    ("SAMPLE", "ls_cloud.txt", "lspoints.txt", 10_000),
    ("ATLAS", "tv_cloud.txt", "tvpoints.txt", 50_000),
):
    full = np.loadtxt(DEBUG / full_file)
    matlab_sub = np.loadtxt(DEBUG / sub_file)
    max_num_points = max(6, int(round(full.shape[0] / divisor)))

    py_sub = _matlab_cloud_subset(full, divisor)
    kd_sub = kdtree_leaf_sample(full, max_num_points)

    print("=" * 72)
    print(f"{label}   full cloud = {full.shape[0]:,} points   maxNumPoints = {max_num_points}")
    print("=" * 72)
    for name, arr in (
        ("MATLAB nonuniformGridSample", matlab_sub),
        ("Python _matlab_cloud_subset", py_sub),
        ("kd-tree leaf sample (proposed)", kd_sub),
    ):
        c = arr.mean(axis=0)
        print(f"  {name:32s} n={arr.shape[0]:>7,}  centroid={c}")
    print(f"  centroid shift  Python vs MATLAB : "
          f"{np.linalg.norm(py_sub.mean(axis=0) - matlab_sub.mean(axis=0)):.2f} vox")
    print(f"  centroid shift  kd-tree vs MATLAB: "
          f"{np.linalg.norm(kd_sub.mean(axis=0) - matlab_sub.mean(axis=0)):.2f} vox")
    print()
