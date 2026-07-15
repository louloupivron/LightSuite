"""Spinal cord B-spline registration (multiobjCordRegistration.m port)."""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import numpy as np
import tifffile
from rich.console import Console
from scipy import ndimage

from lightsuite.config.models import SpinalCordPipelineConfig
from lightsuite.gui.affine import fit_affine_transform, transform_points
from lightsuite.gui.control_points import ControlPointSession
from lightsuite.gui.cord_data import default_cord_session_path
from lightsuite.preprocess.cord_checkpoint import (
    CordRegOptsCheckpoint,
    CordTransformParamsCheckpoint,
)
from lightsuite.registration.cord_affine import warp_cord_atlas_to_straightvol
from lightsuite.registration.cord_longitudinal import (
    load_longitudinal_correspondence,
    resolve_cord_z_transinit,
)
from lightsuite.registration.cord_paths import (
    cord_affine_transform_path,
    cord_bspline_transform_write_path,
    cord_qc_dir,
    cord_save_path,
    cord_work_dir,
)
from lightsuite.registration.cord_plots import save_cord_annotation_preview
from lightsuite.registration.elastix.cord_bspline import build_cord_bspline_params
from lightsuite.registration.elastix.invert import invert_elastix_transform
from lightsuite.registration.elastix.mhd import scale_volume_for_elastix_mi, write_mhd
from lightsuite.registration.elastix.params import write_parameter_file
from lightsuite.registration.elastix.points import write_landmark_file
from lightsuite.registration.elastix.runner import clear_elastix_workspace, run_transformix

console = Console()


def run_spinal_registration(config: SpinalCordPipelineConfig) -> Path:
    """Run cord B-spline elastix registration and write transform_params.json."""
    if shutil.which("elastix") is None or shutil.which("transformix") is None:
        msg = "elastix and transformix must be on PATH for spinal register."
        raise RuntimeError(msg)

    save_path = cord_save_path(config)
    qc_dir = cord_qc_dir(config)
    regopts_path = save_path / "regopts.json"
    checkpoint = CordRegOptsCheckpoint.load(regopts_path)
    if checkpoint.straightvol_path is None or checkpoint.affine_atlas_to_samp is None:
        msg = "Missing init-registration outputs in regopts.json."
        raise RuntimeError(msg)

    straightvol = tifffile.imread(checkpoint.straightvol_path).astype(np.float32)
    tv = tifffile.imread(checkpoint.tv_path).astype(np.float32)
    av = tifffile.imread(checkpoint.av_path).astype(np.uint16)
    transaff = np.asarray(checkpoint.affine_atlas_to_samp, dtype=float)
    nslices = checkpoint.ikeeprange[1] - checkpoint.ikeeprange[0] + 1
    correspondence = load_longitudinal_correspondence(save_path)
    transinit = resolve_cord_z_transinit(nslices, tv.shape[2], correspondence)
    elastix_affine_path = cord_affine_transform_path(config)
    cpwt = config.registration.control_point_weight

    cptshistology = np.zeros((0, 3))
    cptsatlas = np.zeros((0, 3))
    cpaffine = np.zeros((0, 3))
    use_point_aff = False
    cp_path = default_cord_session_path(save_path)
    if cp_path.is_file() and cpwt > 0:
        session = ControlPointSession.load(cp_path)
        cptsatlas, cptshistology = session.paired_points_xyz()
        if cptshistology.shape[0] > 0:
            cpaffine = cptsatlas.copy()
            console.print(f"Found {cptshistology.shape[0]} user-defined control points.")
            if cptshistology.shape[0] > 16:
                atlas_ori = transform_points(
                    cptsatlas,
                    np.linalg.inv(transaff),
                )
                transaff, _ = fit_affine_transform(atlas_ori, cptshistology)
                cpaffine = transform_points(atlas_ori, transaff)
                use_point_aff = True

    spacing_mm = config.registration.resolution_um * 1e-3
    tvtemp = ndimage.median_filter(tv, size=3)
    if use_point_aff:
        from lightsuite.registration.warp import warp_volume_affine

        tvaffine = warp_volume_affine(tvtemp, transaff, straightvol.shape, order=1, point_coords="array")
        avaffine = warp_volume_affine(av.astype(np.float32), transaff, straightvol.shape, order=0, point_coords="array")
    else:
        tvaffine = warp_cord_atlas_to_straightvol(
            tvtemp,
            transinit=transinit,
            elastix_affine_path=elastix_affine_path,
            output_shape=straightvol.shape,
            spacing_mm=spacing_mm,
            work_dir=cord_work_dir(config, "transformix", "register", "tv_affine"),
            nearest=False,
        )
        avaffine = warp_cord_atlas_to_straightvol(
            av,
            transinit=transinit,
            elastix_affine_path=elastix_affine_path,
            output_shape=straightvol.shape,
            spacing_mm=spacing_mm,
            work_dir=cord_work_dir(config, "transformix", "register", "av_affine"),
            nearest=True,
        )

    volmax = float(np.quantile(straightvol, 0.999))
    volplot = np.clip(straightvol / max(volmax, 1.0) * 255.0, 0, 255).astype(np.uint8)
    if use_point_aff:
        save_cord_annotation_preview(
            volplot,
            avaffine.astype(np.uint16),
            qc_dir / "registration_point_affine.png",
        )

    elastix_temp = cord_work_dir(config, "elastix", "bspline")
    clear_elastix_workspace(elastix_temp)

    fixpath = elastix_temp / "fixed.txt"
    movpath = elastix_temp / "moving.txt"
    if cpaffine.size and cptshistology.size:
        volscale = spacing_mm
        write_landmark_file(movpath, (cpaffine - 1.0) * volscale)
        write_landmark_file(fixpath, (cptshistology - 1.0) * volscale)

    has_control_points = cpaffine.size > 0 and cptshistology.size > 0
    if cpwt > 0 and not has_control_points:
        console.print(
            "No control points found; running B-spline with mutual information only "
            f"(control_point_weight={cpwt} ignored until match-points is run)."
        )

    param_path = elastix_temp / "cord_bspline_parameters.txt"
    write_parameter_file(
        param_path,
        build_cord_bspline_params(
            control_point_weight=cpwt,
            fixed_shape=straightvol.shape,
            use_control_points=has_control_points,
        ),
    )

    fixed_mhd = elastix_temp / "fixed"
    moving_mhd = elastix_temp / "moving"
    write_mhd(scale_volume_for_elastix_mi(straightvol), fixed_mhd, [spacing_mm] * 3)
    write_mhd(scale_volume_for_elastix_mi(tvaffine), moving_mhd, [spacing_mm] * 3)

    cmd = [
        "elastix",
        "-f",
        str(fixed_mhd.with_suffix(".mhd")),
        "-m",
        str(moving_mhd.with_suffix(".mhd")),
        "-out",
        str(elastix_temp),
        "-p",
        str(param_path),
    ]
    if has_control_points:
        cmd.extend(["-fp", str(fixpath), "-mp", str(movpath)])

    proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if proc.returncode != 0:
        msg = f"elastix cord B-spline failed:\n{proc.stdout}\n{proc.stderr}"
        raise RuntimeError(msg)

    transforms = sorted(elastix_temp.glob("TransformParameters.*.txt"))
    if not transforms:
        msg = "No B-spline transform written."
        raise RuntimeError(msg)

    avreg = run_transformix(
        moving_volume=avaffine,
        transform_path=transforms[0],
        output_dir=cord_work_dir(config, "transformix", "register", "avreg"),
        spacing_mm=spacing_mm,
        nearest=True,
    )
    save_cord_annotation_preview(
        volplot,
        avreg.astype(np.uint16),
        qc_dir / "registration_bspline.png",
    )

    inv_path = invert_elastix_transform(
        elastix_temp,
        cord_work_dir(config, "elastix", "inverse"),
    )
    bspline_out = cord_bspline_transform_write_path(config)
    bspline_out.write_text(inv_path.read_text(encoding="utf-8"), encoding="utf-8")

    transaff_inv = np.linalg.inv(transaff)
    params = CordTransformParamsCheckpoint(
        tform_bspline_samp20um_to_atlas_20um_px=str(bspline_out),
        tform_affine_samp20um_to_atlas_20um_px=transaff_inv.tolist(),
        control_point_weight=cpwt,
        samp_ikeeplong=checkpoint.ikeeprange,
        samp_ikeepx=checkpoint.xrange,
        samp_ikeepy=checkpoint.yrange,
        how_to_perm=checkpoint.sample_perm,
        slicetforms_path=str(checkpoint.slicetforms_path or ""),
        sampleres_um=checkpoint.sampleres_um,
        registrationres_um=checkpoint.registrationres_um,
        tofliprc=checkpoint.tofliprc,
        atlassize=list(tv.shape),
    )
    out = save_path / "transform_params.json"
    params.save(out)
    console.print(f"Wrote {out}")
    return out
