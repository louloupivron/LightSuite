#!/usr/bin/env python3
"""MATLAB vs Python register parity for Marianna Yosi (register_Yosi_test)."""

from __future__ import annotations

import argparse
import json
import shutil
import sys
from pathlib import Path

import numpy as np
from scipy.io import loadmat

from lightsuite.config.loader import load_config
from lightsuite.gui.control_points import ControlPointSession, load_registration_control_point_session
from lightsuite.import_.matlab_control_points import (
    export_matlab_control_points,
    load_control_point_session_from_mat,
)
from lightsuite.import_.matlab_regopts import import_matlab_regopts, load_matlab_regopts_struct
from lightsuite.preprocess.checkpoint import RegOptsCheckpoint
from lightsuite.registration.brain_register import run_brain_registration
from lightsuite.registration.init_brain import initialize_brain_registration
from lightsuite.registration.orientation import save_orientation

MATLAB_SAVE = Path(
    "/run/user/1002/gvfs/smb-share:server=h1data1.wysscenter.ch,share=computingdata/Alice/0003_CBT_WYSS_LIGHTSHEET/DATA/MESOSPIM/Marianna CMU/MATLAB_python_comparison/register_Yosi_test"
)
CONFIG = Path("examples/config/mesoSPIM/marianna_yosi_parity.yaml")


def _link_or_copy(src: Path, dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists() or dst.is_symlink():
        return
    try:
        dst.symlink_to(src)
    except OSError:
        shutil.copyfile(src, dst)


def setup_python_run() -> Path:
    cfg = load_config(CONFIG)
    py_save = cfg.sample.save_path.expanduser()
    py_save.mkdir(parents=True, exist_ok=True)

    for name in (
        "chan_1_sample_register_20um.tif",
        "chan_2_sample_register_20um.tif",
        "atlas2histology_tform.mat",
        "brain_orientation.txt",
        "regopts.mat",
        "transform_params.mat",
        "regopts_python_export.mat",
    ):
        src = MATLAB_SAVE / name
        if src.is_file():
            _link_or_copy(src, py_save / name)

    orient = np.loadtxt(MATLAB_SAVE / "brain_orientation.txt", delimiter=",").astype(int).ravel().tolist()
    save_orientation(py_save, orient)

    import_matlab_regopts(MATLAB_SAVE, python_save_path=py_save)

    export_mat = py_save / "regopts_python_export.mat"
    checkpoint = RegOptsCheckpoint.load(py_save / "regopts.json")
    if export_mat.is_file():
        raw = loadmat(export_mat, squeeze_me=True)
        checkpoint.original_trans = np.asarray(raw["original_trans_matrix"], dtype=float).tolist()
        checkpoint.save(py_save / "regopts.json")
        print("Using original_trans from regopts_python_export.mat")
    else:
        print("No regopts_python_export.mat — running Python init-registration for original_trans")
        initialize_brain_registration(cfg)
        checkpoint = RegOptsCheckpoint.load(py_save / "regopts.json")
        opts = load_matlab_regopts_struct(MATLAB_SAVE / "regopts.mat")
        checkpoint.autocpsample = np.asarray(opts.autocpsample, dtype=float).tolist()
        checkpoint.autocpatlas = np.asarray(opts.autocpatlas, dtype=float).tolist()
        checkpoint.save(py_save / "regopts.json")
        print(
            "WARNING: original_trans from Python init (not MATLAB). "
            "Run export_regopts_for_python.m on the MATLAB machine for exact parity."
        )

    session = load_control_point_session_from_mat(
        py_save / "atlas2histology_tform.mat",
        original_trans=np.asarray(checkpoint.original_trans, dtype=float),
    )
    session.save(py_save / "atlas2histology_tform.json")
    matched, _, _ = session.point_counts()
    print(f"Python save_path: {py_save}")
    print(f"Manual control-point pairs: {matched}")
    print(f"Auto pairs in regopts: {len(checkpoint.autocpsample or [])}")
    return py_save


def export_cp_for_matlab() -> Path:
    """Write atlas2histology_tform.mat from Python Napari session for MATLAB register."""
    cfg = load_config(CONFIG)
    py_save = cfg.sample.save_path.expanduser()
    json_path = py_save / "atlas2histology_tform.json"
    if not json_path.is_file():
        msg = (
            f"Missing {json_path}. Run:\n"
            f"  uv run lightsuite brain match-points -c {CONFIG}\n"
            "Place pairs in Napari, click Save, then re-run export-cp-for-matlab."
        )
        raise FileNotFoundError(msg)

    out_mat = py_save / "atlas2histology_tform_from_python.mat"
    export_matlab_control_points(json_path, output_mat=out_mat)
    matched, total_sample, total_atlas = ControlPointSession.load(json_path).point_counts()
    print(f"Wrote {out_mat}")
    print(f"Matched pairs: {matched} ({total_sample} sample / {total_atlas} atlas)")
    print("\nMATLAB (new test folder recommended):")
    print("  1. Copy registration volumes + regopts.mat into e.g. register_Yosi_test_from_python/")
    print(f"  2. Copy {out_mat.name} as atlas2histology_tform.mat")
    print("  3. import_python_control_points_for_matlab('.../register_Yosi_test_from_python')")
    print("  4. Run register (multiobjRegistration) with the same cpwt / augment_points flags")
    return out_mat


    cfg = load_config(CONFIG)
    out = run_brain_registration(cfg, use_multistep=True)
    print(f"Wrote {out}")
    return out


def compare_outputs() -> int:
    cfg = load_config(CONFIG)
    py_save = cfg.sample.save_path.expanduser()
    py_tp = json.loads((py_save / "transform_params.json").read_text(encoding="utf-8"))
    py_affine = np.asarray(py_tp["tform_affine_samp20um_to_atlas_10um_px"], dtype=float)
    py_diag = json.loads((py_save / "registration_diagnostics.json").read_text(encoding="utf-8"))
    py_affine_diag = json.loads((py_save / "affine_fit_stats.json").read_text(encoding="utf-8"))

    checkpoint, session = (
        RegOptsCheckpoint.load(py_save / "regopts.json"),
        load_registration_control_point_session(
            py_save,
            original_trans=RegOptsCheckpoint.load(py_save / "regopts.json").original_trans,
            prefer_matlab=True,
        ),
    )
    atlas_pts, sample_pts = session.paired_points_xyz()

    print("=" * 72)
    print("CONTROL POINTS")
    print("=" * 72)
    print(f"Manual pairs loaded: {atlas_pts.shape[0]}")
    print(f"Auto pairs in regopts: {len(checkpoint.autocpsample or [])}")
    print(f"Landmarks sent to elastix: {py_diag['n_landmark_pairs']}")
    print(f"cpwt: Python={py_diag['control_point_weight']}")

    mat_tp_path = MATLAB_SAVE / "transform_params.mat"
    mat_export = MATLAB_SAVE / "regopts_python_export.mat"
    print("\n" + "=" * 72)
    print("AFFINE / REGISTER OUTPUT")
    print("=" * 72)
    print(
        f"Python affine median error: {py_affine_diag['median_error_vox']:.2f} vox "
        f"(manual={py_affine_diag.get('median_error_manual_vox')}, "
        f"auto={py_affine_diag.get('median_error_auto_vox')})"
    )
    if py_diag.get("bspline_landmark_metric_vox") is not None:
        print(f"Python B-spline landmark metric: {py_diag['bspline_landmark_metric_vox']:.2f} vox")

    if mat_export.is_file():
        raw = loadmat(mat_export, squeeze_me=True)
        mat_affine = np.asarray(raw["tform_affine_samp20um_to_atlas_10um_px_matrix"], dtype=float)
        lin_diff = np.abs(mat_affine[:3, :3] - py_affine[:3, :3]).max()
        trans_diff = np.linalg.norm(mat_affine[:3, 3] - py_affine[:3, 3])
        print(f"MATLAB vs Python affine (sample->atlas):")
        print(f"  max |linear diff| = {lin_diff:.4f}")
        print(f"  |translation diff| = {trans_diff:.2f} (20 µm voxels)")
    else:
        tp = loadmat(mat_tp_path, squeeze_me=True, struct_as_record=False)
        print(f"MATLAB transform_params.contol_pt_weight = {float(tp['contol_pt_weight'])}")
        print(
            "MATLAB affine matrix not readable (affinetform3d). "
            "Run export_regopts_for_python.m for numeric comparison."
        )

    mat_fix = MATLAB_SAVE / "elastix_temp" / "fixed.txt"
    py_fix = py_save / "elastix_temp" / "fixed.txt"
    if mat_fix.is_file() and py_fix.is_file():
        mat_pts = _read_elastix_points(mat_fix)
        py_pts = _read_elastix_points(py_fix)
        n = min(mat_pts.shape[0], py_pts.shape[0])
        if n:
            err = np.linalg.norm(mat_pts[:n] - py_pts[:n], axis=1)
            print("\n" + "=" * 72)
            print("ELASTIX FIXED LANDMARKS (sample side, mm)")
            print("=" * 72)
            print(f"pairs: MATLAB={mat_pts.shape[0]} Python={py_pts.shape[0]}")
            print(f"median |diff| = {np.median(err):.4f} mm  p95 = {np.percentile(err, 95):.4f} mm")

    return 0


def _read_elastix_points(path: Path) -> np.ndarray:
    lines = path.read_text(encoding="utf-8").splitlines()
    n = int(lines[1].strip())
    pts = []
    for line in lines[2 : 2 + n]:
        pts.append([float(x) for x in line.split()])
    return np.asarray(pts, dtype=float)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "command",
        choices=["setup", "register", "compare", "export-cp-for-matlab", "all"],
    )
    args = parser.parse_args(argv)

    if args.command == "export-cp-for-matlab":
        export_cp_for_matlab()
        return 0
    if args.command in ("setup", "all"):
        setup_python_run()
    if args.command in ("register", "all"):
        run_register()
    if args.command in ("compare", "all"):
        return compare_outputs()
    return 0


if __name__ == "__main__":
    sys.exit(main())
