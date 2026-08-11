"""Diagnose Python vs MATLAB coarse-alignment parity (scratch, not part of the package)."""

from __future__ import annotations

import numpy as np

np.set_printoptions(precision=4, suppress=True)

# MATLAB original_trans (affinetform3d.A, premultiply), latest MATLAB run
T_MATLAB = np.array(
    [
        [0.9271, 0.1870, 0.0514, -71.9302],
        [-0.1825, 0.9261, -0.0782, 4.1789],
        [-0.0657, 0.0667, 0.9425, -130.3184],
        [0.0, 0.0, 0.0, 1.0],
    ]
)

# Python original_trans from the 12:55 run (regopts.json)
T_PYTHON = np.array(
    [
        [0.9325700266355924, 0.18972203405646193, 0.05799145099551907, -48.88133797205464],
        [-0.18484604742066837, 0.9321659879351032, -0.0770897066029254, -41.16350889850832],
        [-0.07203747183267195, 0.06415944127443594, 0.9485454559068348, -136.21130293721527],
        [0.0, 0.0, 0.0, 1.0],
    ]
)


def report_linear(name: str, T: np.ndarray) -> None:
    M = T[:3, :3]
    scale = float(np.cbrt(np.linalg.det(M)))
    R = M / scale
    angle = np.degrees(np.arccos(np.clip((np.trace(R) - 1) / 2, -1, 1)))
    print(f"{name}: scale={scale:.4f}  rot_angle={angle:.2f} deg  t={T[:3, 3]}")


print("=" * 70)
print("1. LINEAR PART COMPARISON")
print("=" * 70)
report_linear("MATLAB", T_MATLAB)
report_linear("Python", T_PYTHON)
print(f"linear diff (max abs) : {np.abs(T_MATLAB[:3, :3] - T_PYTHON[:3, :3]).max():.4f}")
print(f"translation diff      : {T_PYTHON[:3, 3] - T_MATLAB[:3, 3]}")
print(f"translation diff norm : {np.linalg.norm(T_PYTHON[:3, 3] - T_MATLAB[:3, 3]):.2f} vox")

print()
print("=" * 70)
print("2. INVERSION-CONVENTION HYPOTHESIS")
print("=" * 70)
# align.py does: aff = [[L, t],[0,1]] with forward p' = p @ L + t   (row convention)
#                sample_to_atlas = np.linalg.inv(aff)
#                transform_icp   = transpose of [:3,:3]
# => Python translation = -L^-1 t
# MATLAB does simtform3d.invert on the PROPER premultiply matrix [[L.T, t],[0,1]]
# => MATLAB translation = -(L^-1).T t
# Relation: t_python = M.T @ M^-1 @ t_matlab  = (R^2).T @ t_matlab
M = T_MATLAB[:3, :3]
predicted_t = M.T @ np.linalg.inv(M) @ T_MATLAB[:3, 3]
print(f"MATLAB translation             : {T_MATLAB[:3, 3]}")
print(f"predicted Python (bug model)   : {predicted_t}")
print(f"observed Python translation    : {T_PYTHON[:3, 3]}")
print(f"residual (predicted - observed): {predicted_t - T_PYTHON[:3, 3]}")
print(f"residual norm                  : {np.linalg.norm(predicted_t - T_PYTHON[:3, 3]):.2f} vox")
print(f"(compare with raw gap of       : {np.linalg.norm(T_PYTHON[:3, 3] - T_MATLAB[:3, 3]):.2f} vox)")

print()
print("=" * 70)
print("3. MATLAB nonuniformGridSample OUTPUT-SIZE MODEL")
print("=" * 70)
# MATLAB exported counts
cases = [
    ("sample (ls_cloud -> lspoints)", 523_993, 10_000, 52, 16_384),
    ("atlas  (tv_cloud -> tvpoints)", 10_691_773, 50_000, 214, 65_536),
]
for name, count, divisor, grid_arg_expected, out_expected in cases:
    grid_arg = max(6, int(round(count / divisor)))
    n_boxes_needed = count / grid_arg
    predicted_out = 2 ** int(np.ceil(np.log2(n_boxes_needed)))
    print(f"{name}")
    print(f"   Count={count:,}  maxNumPoints=max(6,round(Count/{divisor:.0e}))={grid_arg}"
          f"  (MATLAB reported {grid_arg_expected})")
    print(f"   kd-tree leaves = 2^ceil(log2({n_boxes_needed:.0f})) = {predicted_out:,}"
          f"  (MATLAB exported {out_expected:,})  "
          f"{'MATCH' if predicted_out == out_expected else 'MISMATCH'}")

print()
print("=" * 70)
print("4. PYTHON _matlab_cloud_subset ON THE SAME COUNTS")
print("=" * 70)
for name, count, divisor, _g, out_expected in cases:
    target = int(np.ceil(count / divisor))
    voxel = max(1.0, (count / target) ** (1 / 3))
    print(f"{name}")
    print(f"   Python target = ceil({count:,}/{divisor}) = {target}  -> voxel_down_sample(voxel={voxel:.1f})")
    print(f"   MATLAB feeds BCPD {out_expected:,} points; Python aims at ~{target} "
          f"(actual output set by voxel size, still orders of magnitude off)")
