"""Decisive voluse injection test: MATLAB volumereg -> Python mask extraction.

Copy the folder produced on the MATLAB machine (``voluse_export/``) into the
Python save path, then run::

    uv run python scratch_inject_matlab_voluse.py

Interpretation
--------------
- MATLAB volumereg -> Python mask ~5,248,642: extraction port is correct;
  fix upstream volume construction (normalize / permute / regvol).
- MATLAB volumereg -> Python mask still ~3,625,046: bug is inside mask logic;
  diff per-batch thresinit next.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from lightsuite.registration.points import (
    _extract_sample_mask_points,
    _volume_mode_all,
    extract_sample_points_stages,
)
from lightsuite.registration.volume import (
    load_registration_volume,
    normalize_registration_volume,
    permute_brain_volume,
)

CONFIG = Path("examples/config/mesoSPIM/marianna.yaml")
PY_SAVE = Path(
    "/media/gbm/NVME2/ALICe-pipelines-data/Marianna_CMU/registered_allen_MATLAB_conterpart"
)
EXPORT_DIR = PY_SAVE / "voluse_export"
THRESH = 5.0
MATLAB_TRIM_TARGET = 5_248_642
MATLAB_DENOISE_TARGET = 523_993


def load_matlab_volumereg(export_dir: Path) -> tuple[np.ndarray, dict]:
    meta_path = export_dir / "meta.json"
    vol_path = export_dir / "volumereg.f32"
    if not meta_path.is_file():
        raise FileNotFoundError(f"Missing {meta_path}")
    if not vol_path.is_file():
        raise FileNotFoundError(f"Missing {vol_path}")

    meta = json.loads(meta_path.read_text())
    shape = tuple(int(v) for v in meta["size_yxz"])
    flat = np.fromfile(vol_path, dtype=np.float32)
    expected = int(np.prod(shape))
    if flat.size != expected:
        msg = f"{vol_path.name}: expected {expected} values, got {flat.size}"
        raise ValueError(msg)
    # MATLAB column-major flattening == NumPy Fortran order.
    volume = flat.reshape(shape, order="F")
    return volume.astype(np.float32, copy=False), meta


def volume_stats(volume: np.ndarray, label: str) -> None:
    flat = volume.ravel(order="F")
    rng = np.random.default_rng(1)
    sample_n = min(10_000, flat.size)
    idx = rng.choice(flat.size, size=sample_n, replace=False)
    q05 = float(np.quantile(flat[idx], 0.05)) * 2.0
    print(f"\n{label}")
    print(f"  shape (Y,X,Z): {tuple(volume.shape)}")
    print(f"  dtype: {volume.dtype}")
    print(f"  min/max/mean: {volume.min():.6g} / {volume.max():.6g} / {volume.mean():.6g}")
    print(f"  mode(all): {_volume_mode_all(volume):.6g}")
    print(f"  global 5%*2 (rng=1): {q05:.6g}")
    print(f"  max(thres, mode): {max(q05, _volume_mode_all(volume)):.6g}")


def mask_counts(volume: np.ndarray, threshold: float) -> dict[str, int]:
    mask_pts, trim_pts = _extract_sample_mask_points(volume, threshold)
    _, stages = extract_sample_points_stages(volume, threshold, subsample_fraction=0.1)
    return {
        "mask_before_trim": int(mask_pts.shape[0]),
        "mask_after_trim": int(trim_pts.shape[0]),
        "downsample": stages.downsample_points,
        "denoise": stages.denoise_points,
    }


def load_permvec(save_path: Path) -> list[int]:
    orient_file = save_path / "brain_orientation.txt"
    if orient_file.is_file():
        vals = [int(v) for v in orient_file.read_text().split()]
        if len(vals) == 3:
            return vals
    return [1, 2, 3]


def build_python_volumereg(save_path: Path) -> np.ndarray:
    regopts = json.loads((save_path / "regopts.json").read_text())
    backvol = load_registration_volume(Path(regopts["regvolpath"]))
    permvec = load_permvec(save_path)
    newvol = normalize_registration_volume(backvol)
    return permute_brain_volume(newvol, permvec)


def compare_volumes(ml_vol: np.ndarray, py_vol: np.ndarray) -> None:
    if ml_vol.shape != py_vol.shape:
        print(f"\nShape mismatch: MATLAB {ml_vol.shape} vs Python {py_vol.shape}")
        return
    diff = ml_vol.astype(np.float64) - py_vol.astype(np.float64)
    corr = float(np.corrcoef(ml_vol.ravel(order="F"), py_vol.ravel(order="F"))[0, 1])
    print("\nMATLAB vs Python volumereg elementwise")
    print(f"  Pearson r: {corr:.8f}")
    print(f"  abs diff: max={np.max(np.abs(diff)):.6g}  mean={np.mean(np.abs(diff)):.6g}")
    print(f"  rel diff (|d|/|ml|): median={np.median(np.abs(diff)/(np.abs(ml_vol)+1e-12)):.6g}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--export-dir",
        type=Path,
        default=EXPORT_DIR,
        help="Folder with volumereg.f32 and meta.json from MATLAB",
    )
    args = parser.parse_args()
    export_dir = args.export_dir.expanduser()

    print("=" * 72)
    print("VOLUSE INJECTION TEST")
    print("=" * 72)
    print(f"Export dir: {export_dir}")

    ml_vol, meta = load_matlab_volumereg(export_dir)
    print("\nMATLAB export meta:")
    for key in sorted(meta):
        print(f"  {key}: {meta[key]}")

    py_vol = build_python_volumereg(PY_SAVE)

    volume_stats(ml_vol, "MATLAB volumereg (exported)")
    volume_stats(py_vol, "Python volumereg (current pipeline)")
    compare_volumes(ml_vol, py_vol)

    print("\n" + "=" * 72)
    print("MASK STAGE (threshold = 5)")
    print("=" * 72)

    ml_counts = mask_counts(ml_vol, THRESH)
    py_counts = mask_counts(py_vol, THRESH)

    def row(name: str, ml_n: int, py_n: int) -> None:
        print(f"  {name:22s}  MATLAB-via-Py {ml_n:>10,}   Python {py_n:>10,}")

    row("before trim", ml_counts["mask_before_trim"], py_counts["mask_before_trim"])
    row("after trim", ml_counts["mask_after_trim"], py_counts["mask_after_trim"])
    row("after 10% downsample", ml_counts["downsample"], py_counts["downsample"])
    row("after denoise", ml_counts["denoise"], py_counts["denoise"])

    print(f"\nMATLAB measured targets: trim ~{MATLAB_TRIM_TARGET:,}, denoise ~{MATLAB_DENOISE_TARGET:,}")

    ml_trim = ml_counts["mask_after_trim"]
    if abs(ml_trim - MATLAB_TRIM_TARGET) / MATLAB_TRIM_TARGET < 0.02:
        print("\nVERDICT: MATLAB volumereg reproduces MATLAB mask count in Python.")
        print("         -> Fix upstream volume construction (normalize / permute / regvol).")
    else:
        print(
            f"\nVERDICT: MATLAB volumereg -> Python trim = {ml_trim:,} "
            f"(expected ~{MATLAB_TRIM_TARGET:,})."
        )
        print("         -> Bug is inside mask extraction logic; diff per-batch thresinit.")


if __name__ == "__main__":
    main()
