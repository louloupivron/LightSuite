"""Per-slice straightening transforms (computeStraighteningTransforms.m port)."""

from __future__ import annotations

import math

import numpy as np
from scipy import ndimage


def circ_dist(a: float, b: float) -> float:
    """Signed circular distance a - b wrapped to (-pi, pi]."""
    return math.atan2(math.sin(a - b), math.cos(a - b))


def compute_straightening_transforms(
    fit_x: np.ndarray,
    fit_y: np.ndarray,
    fit_theta: np.ndarray,
    target_center: tuple[float, float],
    target_orientation_deg: float = 90.0,
) -> list[np.ndarray]:
    """Return per-slice 3x3 homogeneous rigid transforms (MATLAB pre-multiply)."""
    centers = np.column_stack([fit_x, fit_y])
    theta = fit_theta + math.pi
    if np.any(np.isnan(centers)):
        centers = _fill_nearest_2d(centers)
        theta = _fill_nearest_1d(theta)
    theta_target = math.radians(target_orientation_deg)
    tforms: list[np.ndarray] = []
    tx, ty = target_center
    for i in range(centers.shape[0]):
        cx, cy = centers[i]
        theta_rot = circ_dist(theta_target, float(theta[i]))
        c, s = math.cos(theta_rot), math.sin(theta_rot)
        t_to_origin = np.array([[1, 0, -cx], [0, 1, -cy], [0, 0, 1]], dtype=float)
        rot = np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]], dtype=float)
        t_to_target = np.array([[1, 0, tx], [0, 1, ty], [0, 0, 1]], dtype=float)
        tforms.append(t_to_target @ rot @ t_to_origin)
    return tforms


def transform_cord_images_slices(
    cordvol: np.ndarray,
    tforms: list[np.ndarray],
    output_size: tuple[int, int],
    *,
    fill_value: float | None = None,
) -> np.ndarray:
    """Apply per-slice 2D transforms (tranformCordImagesSlices.m)."""
    if fill_value is None:
        sample = cordvol.ravel()
        pos = sample[sample > 0]
        fill_value = float(np.bincount(pos.astype(np.int64)).argmax()) if pos.size else 0.0
    out_h, out_w = output_size
    volout = np.zeros((out_h, out_w, len(tforms)), dtype=cordvol.dtype)
    for ii, tform in enumerate(tforms):
        volout[:, :, ii] = _warp_slice_2d(cordvol[:, :, ii], tform, (out_h, out_w), fill_value)
    return volout


def transform_cord_points_slices(points: np.ndarray, tforms: list[np.ndarray]) -> np.ndarray:
    """Transform XY coordinates per slice (tranformCordPointsSlices.m)."""
    pts = np.asarray(points, dtype=float).copy()
    slice_ids = pts[:, 2].astype(int)
    for islice, tform in enumerate(tforms, start=1):
        mask = slice_ids == islice
        if not np.any(mask):
            continue
        xy = pts[mask, :2]
        hom = np.column_stack([xy, np.ones(xy.shape[0])])
        warped = (tform @ hom.T).T
        pts[mask, :2] = warped[:, :2]
    return pts


def save_slicetforms(path, tforms: list[np.ndarray]) -> None:
    """Persist 3x3 transforms as .npy stack."""
    from pathlib import Path

    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    np.save(p, np.stack(tforms, axis=0))


def load_slicetforms(path) -> list[np.ndarray]:
    from pathlib import Path

    arr = np.load(Path(path))
    return [arr[i] for i in range(arr.shape[0])]


def _warp_slice_2d(
    image: np.ndarray,
    tform: np.ndarray,
    output_shape: tuple[int, int],
    fill_value: float,
) -> np.ndarray:
    """Apply per-slice rigid transform matching MATLAB ``imwarp`` + ``imref2d``.

    ``tform`` maps input (x, y) = (column, row) to output (x, y). ``affine_transform``
    indexes (row, column), so swap axes when building the inverse resampling matrix.
    """
    inv = np.linalg.inv(tform)
    m_xy = inv[:2, :2]
    b_xy = inv[:2, 2]
    swap_rc_xy = np.array([[0.0, 1.0], [1.0, 0.0]])
    matrix = swap_rc_xy @ m_xy @ swap_rc_xy
    offset = swap_rc_xy @ b_xy
    return ndimage.affine_transform(
        image,
        matrix,
        offset=offset,
        output_shape=output_shape,
        order=1,
        mode="constant",
        cval=fill_value,
    )


def _fill_nearest_1d(values: np.ndarray) -> np.ndarray:
    out = values.astype(float).copy()
    valid = ~np.isnan(out)
    if not np.any(valid):
        return out
    idx = np.flatnonzero(valid)
    all_idx = np.arange(out.size)
    out[~valid] = np.interp(all_idx[~valid], idx, out[valid])
    return out


def _fill_nearest_2d(values: np.ndarray) -> np.ndarray:
    out = values.astype(float).copy()
    for col in range(out.shape[1]):
        out[:, col] = _fill_nearest_1d(out[:, col])
    return out
