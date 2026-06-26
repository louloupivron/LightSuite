"""Spinal cord z-scale + elastix affine warping (initializeCordRegistration.m port)."""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np

from lightsuite.gui.affine import fit_affine_transform, transform_points
from lightsuite.registration.elastix.runner import run_transformix
from lightsuite.registration.warp import warp_volume_affine


def build_cord_z_transinit(nslices: int, atlas_depth: int) -> np.ndarray:
    """Initial z-stretch mapping atlas depth onto the straightened sample grid."""
    z_scale = nslices * 0.98 / atlas_depth
    z_trans = nslices / 2 - z_scale * atlas_depth / 2
    transinit = np.eye(4, dtype=float)
    transinit[2, 2] = z_scale
    transinit[2, 3] = z_trans
    return transinit


def _parse_elastix_moving_to_fixed_matrix(path: Path) -> np.ndarray:
    """Return 4x4 moving→fixed transform in ITK physical coordinates (0-based origin)."""
    raw = path.expanduser().read_text(encoding="utf-8")
    pat_params = r"\(\s*TransformParameters\s+([^)]+)\)"
    pat_center = r"\(\s*CenterOfRotationPoint\s+([^)]+)\)"
    match = re.search(pat_params, raw)
    if match is None:
        msg = f"TransformParameters missing in {path}"
        raise ValueError(msg)
    params = np.fromstring(match.group(1), sep=" ")
    center = np.fromstring(re.search(pat_center, raw).group(1), sep=" ") if re.search(pat_center, raw) else np.zeros(3)
    r_mat = params[:9].reshape(3, 3).T
    t_vec = params[9:12]
    t_eff = t_vec + center - (r_mat @ center)
    m_elastix = np.eye(4, dtype=float)
    m_elastix[:3, :3] = r_mat
    m_elastix[:3, 3] = t_eff
    return np.linalg.inv(m_elastix)


def _volume_indices_to_itk_physical(points_yxz: np.ndarray, spacing_mm: float) -> np.ndarray:
    """Map array indices (Y, X, Z) to ITK physical coordinates."""
    sp = float(spacing_mm)
    return np.column_stack([points_yxz[:, 1] * sp, points_yxz[:, 0] * sp, points_yxz[:, 2] * sp])


def _itk_physical_to_volume_indices(points_xyz: np.ndarray, spacing_mm: float) -> np.ndarray:
    sp = float(spacing_mm)
    return np.column_stack([points_xyz[:, 1] / sp, points_xyz[:, 0] / sp, points_xyz[:, 2] / sp])


def apply_elastix_affine_volume_indices(
    points_yxz: np.ndarray,
    transform_path: Path,
    spacing_mm: float,
) -> np.ndarray:
    """Apply elastix affine (moving→fixed) to 0-based (Y, X, Z) voxel indices."""
    if points_yxz.size == 0:
        return points_yxz.reshape(0, 3)
    matrix = _parse_elastix_moving_to_fixed_matrix(transform_path)
    pts = np.asarray(points_yxz, dtype=float)
    hom = np.column_stack([_volume_indices_to_itk_physical(pts, spacing_mm), np.ones(pts.shape[0])])
    mapped = (hom @ matrix.T)[:, :3]
    return _itk_physical_to_volume_indices(mapped, spacing_mm)


def fit_cord_affine_atlas_to_straightvol(
    atlas_shape: tuple[int, int, int],
    transinit: np.ndarray,
    elastix_affine_path: Path,
    spacing_mm: float,
) -> np.ndarray:
    """Approximate composed atlas→straightvol affinetform for control-point tooling."""
    y_idx = np.linspace(0, atlas_shape[0] - 1, 8)
    x_idx = np.linspace(0, atlas_shape[1] - 1, 6)
    z_idx = np.linspace(0, atlas_shape[2] - 1, 8)
    yy, xx, zz = np.meshgrid(y_idx, x_idx, z_idx, indexing="ij")
    atlas_pts = np.column_stack([yy.ravel(), xx.ravel(), zz.ravel()])
    after_z = transform_points(atlas_pts, transinit)
    after_affine = apply_elastix_affine_volume_indices(after_z, elastix_affine_path, spacing_mm)
    matrix, _ = fit_affine_transform(atlas_pts, after_affine)
    return matrix


def write_inverse_elastix_affine(forward_path: Path, out_path: Path) -> Path:
    """Write the analytic inverse of an elastix ``AffineTransform`` parameter file.

    Elastix stores ``T(x) = A·(x − c) + t + c``. The inverse about the same centre ``c``
    is ``A⁻¹`` with translation ``−A⁻¹·t``; only ``TransformParameters`` changes, so the
    grid (``Size``/``Spacing``/``Origin``) stays identical and transformix resamples a
    sample-space volume back onto the (z-scaled) atlas-affine grid.
    """
    forward_path = Path(forward_path).expanduser()
    out_path = Path(out_path).expanduser()
    raw = forward_path.read_text(encoding="utf-8")
    match = re.search(r"\(\s*TransformParameters\s+([^)]+)\)", raw)
    if match is None:
        msg = f"TransformParameters missing in {forward_path}"
        raise ValueError(msg)
    params = np.fromstring(match.group(1), sep=" ")
    if params.size < 12:
        msg = f"Expected 12 affine parameters in {forward_path}, got {params.size}"
        raise ValueError(msg)
    a_mat = params[:9].reshape(3, 3)
    t_vec = params[9:12]
    a_inv = np.linalg.inv(a_mat)
    t_inv = -a_inv @ t_vec
    inv_params = np.concatenate([a_inv.reshape(-1), t_inv])
    inv_str = " ".join(f"{value:.16g}" for value in inv_params)
    text = re.sub(
        r"\(\s*TransformParameters\s+[^)]+\)",
        f"(TransformParameters {inv_str})",
        raw,
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(text, encoding="utf-8")
    return out_path


def warp_cord_straightvol_to_atlas(
    volume: np.ndarray,
    *,
    transinit: np.ndarray,
    elastix_affine_path: Path,
    atlas_shape: tuple[int, int, int],
    spacing_mm: float,
    work_dir: Path,
    nearest: bool,
) -> np.ndarray:
    """Warp a straightened-sample-grid volume into atlas space (inverse of
    :func:`warp_cord_atlas_to_straightvol`).

    The forward atlas→sample warp is ``atlas → transinit (z-scale) → elastix affine``.
    Inverting it requires applying the *inverse* elastix affine via transformix (the only
    method that matches elastix' XY affine) followed by the inverse z-scale matrix. Using a
    single composed matrix cannot reproduce the elastix affine and produces a sheared result.
    """
    order = 0 if nearest else 1
    work_dir = Path(work_dir).expanduser()
    work_dir.mkdir(parents=True, exist_ok=True)
    inverse_affine = write_inverse_elastix_affine(
        Path(elastix_affine_path),
        work_dir / "inverse_affine.txt",
    )
    vol_affine = run_transformix(
        moving_volume=np.asarray(volume, dtype=np.float32),
        transform_path=inverse_affine,
        output_dir=work_dir,
        spacing_mm=spacing_mm,
        nearest=nearest,
    )
    registered = warp_volume_affine(
        np.asarray(vol_affine, dtype=np.float32),
        np.linalg.inv(transinit),
        atlas_shape,
        order=order,
        point_coords="array",
    )
    if tuple(registered.shape) != tuple(atlas_shape):
        msg = (
            f"Warped sample shape {registered.shape} does not match atlas shape {atlas_shape}"
        )
        raise ValueError(msg)
    return registered


def warp_cord_atlas_to_straightvol(
    volume: np.ndarray,
    *,
    transinit: np.ndarray,
    elastix_affine_path: Path | None,
    output_shape: tuple[int, int, int],
    spacing_mm: float,
    work_dir: Path,
    nearest: bool,
) -> np.ndarray:
    """Warp atlas annotation/template onto the straightened sample grid."""
    order = 0 if nearest else 1
    vol_z = warp_volume_affine(
        np.asarray(volume, dtype=np.float32),
        transinit,
        output_shape,
        order=order,
        point_coords="array",
    )
    if elastix_affine_path is None or not Path(elastix_affine_path).expanduser().is_file():
        out = vol_z.astype(volume.dtype, copy=False)
    else:
        warped = run_transformix(
            moving_volume=vol_z,
            transform_path=Path(elastix_affine_path),
            output_dir=work_dir,
            spacing_mm=spacing_mm,
            nearest=nearest,
        )
        if nearest:
            out = np.rint(warped).astype(volume.dtype, copy=False)
        else:
            out = warped.astype(np.float32, copy=False)

    if tuple(out.shape) != tuple(output_shape):
        msg = (
            f"Warped cord atlas shape {out.shape} does not match requested output "
            f"{output_shape}"
        )
        raise ValueError(msg)
    return out
