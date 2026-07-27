"""QA plots for manifest-driven multiresolution geometry checks."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import SimpleITK as sitk

from lightsuite.multires.geometry import transform_physical_points


def xy_slice_at_physical_um(
    image_geo: sitk.Image,
    cx_um: float,
    cy_um: float,
    cz_um: float,
) -> np.ndarray:
    """One XY plane as numpy 2D at physical (cx, cy, cz)."""
    point = (float(cx_um), float(cy_um), float(cz_um))
    idx = image_geo.TransformPhysicalPointToContinuousIndex(point)
    _i, _j, k = (int(round(float(t))) for t in idx)
    size = image_geo.GetSize()
    k = int(np.clip(k, 0, size[2] - 1))
    slice_img = sitk.RegionOfInterest(image_geo, [size[0], size[1], 1], [0, 0, k])
    return np.asarray(sitk.GetArrayFromImage(slice_img)[0], dtype=np.float32)


def _normalize_panel(image: np.ndarray) -> np.ndarray:
    data = image.astype(np.float32, copy=False)
    if data.size == 0 or data.max() <= 0:
        return data
    positive = data[data > 0]
    hi = float(np.quantile(positive, 0.995)) if positive.size else float(data.max())
    return np.clip(data / max(hi, 1e-6), 0, 1)


def _resample_to_shape(image: np.ndarray, target_shape: tuple[int, int]) -> np.ndarray:
    """Resize a 2D panel to ``target_shape`` (rows, cols) for overlay comparison."""
    if image.shape == target_shape:
        return image.astype(np.float32, copy=False)
    from scipy.ndimage import zoom

    zoom_y = target_shape[0] / image.shape[0]
    zoom_x = target_shape[1] / image.shape[1]
    return np.asarray(zoom(image.astype(np.float32), (zoom_y, zoom_x), order=1), dtype=np.float32)


def normalized_cross_correlation(image_a: np.ndarray, image_b: np.ndarray) -> float:
    """Pearson correlation between two equally sized 2D panels."""
    if image_a.shape != image_b.shape:
        msg = f"NCC requires matching shapes, got {image_a.shape} and {image_b.shape}"
        raise ValueError(msg)
    a = image_a.astype(np.float64, copy=False).ravel()
    b = image_b.astype(np.float64, copy=False).ravel()
    a = a - a.mean()
    b = b - b.mean()
    denom = float(np.sqrt((a * a).sum() * (b * b).sum()))
    if denom < 1e-12:
        return 0.0
    return float((a * b).sum() / denom)


def _overlap_center_from_sitk(
    overview: sitk.Image,
    roi: sitk.Image,
    overlap_min: np.ndarray,
    overlap_max: np.ndarray,
    *,
    roi_to_overview: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray, tuple[float, float, float]]:
    center_um = 0.5 * (overlap_min + overlap_max)
    cx, cy, cz = (float(center_um[0]), float(center_um[1]), float(center_um[2]))
    sl_overview = xy_slice_at_physical_um(overview, cx, cy, cz)
    if roi_to_overview is None:
        sl_roi = xy_slice_at_physical_um(roi, cx, cy, cz)
    else:
        roi_center = transform_physical_points(
            np.array([[cx, cy, cz]], dtype=float),
            np.linalg.inv(roi_to_overview),
        )[0]
        sl_roi = xy_slice_at_physical_um(
            roi,
            float(roi_center[0]),
            float(roi_center[1]),
            float(roi_center[2]),
        )
    return sl_overview, sl_roi, (cx, cy, cz)


def save_fov_overlap_plot(
    *,
    rep_overview: dict[str, object],
    rep_roi: dict[str, object],
    overlap_min: np.ndarray,
    overlap_max: np.ndarray,
    output_path: Path,
    title: str,
) -> None:
    """Draw axis-aligned FOV boxes (XY + XZ) with volume centers marked."""
    output_path = output_path.expanduser()
    output_path.parent.mkdir(parents=True, exist_ok=True)

    def _as_vec(report: dict[str, object], key: str) -> np.ndarray:
        return np.asarray(report[key], dtype=float)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    def rect_xy(ax, pmin: np.ndarray, pmax: np.ndarray, **kwargs) -> None:
        w, h = float(pmax[0] - pmin[0]), float(pmax[1] - pmin[1])
        ax.add_patch(plt.Rectangle((float(pmin[0]), float(pmin[1])), w, h, fill=False, **kwargs))

    def rect_xz(ax, pmin: np.ndarray, pmax: np.ndarray, **kwargs) -> None:
        w, h = float(pmax[0] - pmin[0]), float(pmax[2] - pmin[2])
        ax.add_patch(plt.Rectangle((float(pmin[0]), float(pmin[2])), w, h, fill=False, **kwargs))

    o_min, o_max = _as_vec(rep_overview, "phys_min"), _as_vec(rep_overview, "phys_max")
    r_min, r_max = _as_vec(rep_roi, "phys_min"), _as_vec(rep_roi, "phys_max")
    o_center = _as_vec(rep_overview, "phys_center")
    r_center = _as_vec(rep_roi, "phys_center")

    ax = axes[0]
    rect_xy(ax, o_min, o_max, color="C0", lw=2, label=str(rep_overview.get("label", "overview")))
    rect_xy(ax, r_min, r_max, color="C1", lw=2, label=str(rep_roi.get("label", "roi")))
    rect_xy(ax, overlap_min, overlap_max, color="C2", lw=2, ls="--", label="overlap")
    ax.scatter(float(o_center[0]), float(o_center[1]), c="C0", s=40, zorder=5)
    ax.scatter(float(r_center[0]), float(r_center[1]), c="C1", s=40, zorder=5)
    ax.set_xlabel("x (µm)")
    ax.set_ylabel("y (µm)")
    ax.set_title("Axial (XY)")
    ax.invert_yaxis()
    ax.set_aspect("equal")
    ax.legend(loc="best")

    ax = axes[1]
    rect_xz(ax, o_min, o_max, color="C0", lw=2, label=str(rep_overview.get("label", "overview")))
    rect_xz(ax, r_min, r_max, color="C1", lw=2, label=str(rep_roi.get("label", "roi")))
    rect_xz(ax, overlap_min, overlap_max, color="C2", lw=2, ls="--", label="overlap")
    ax.scatter(float(o_center[0]), float(o_center[2]), c="C0", s=40, zorder=5)
    ax.scatter(float(r_center[0]), float(r_center[2]), c="C1", s=40, zorder=5)
    ax.set_xlabel("x (µm)")
    ax.set_ylabel("z (µm)")
    ax.set_title("Coronal (XZ)")
    ax.set_aspect("equal")
    ax.legend(loc="best")

    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)


def save_geometry_overlap_qc_plot(
    *,
    overview: sitk.Image,
    roi: sitk.Image,
    overlap_min: np.ndarray,
    overlap_max: np.ndarray,
    output_path: Path,
    geometry_mode: str,
    roi_to_overview: np.ndarray | None = None,
    overview_crop: sitk.Image | None = None,
    roi_crop: sitk.Image | None = None,
) -> None:
    """Side-by-side overlap crops at the shared overlap center."""
    center_um = tuple(float(v) for v in 0.5 * (overlap_min + overlap_max))
    if overview_crop is not None and roi_crop is not None:
        overview_arr = sitk.GetArrayViewFromImage(overview_crop)
        roi_arr = sitk.GetArrayViewFromImage(roi_crop)
        sl_overview = np.asarray(overview_arr[overview_arr.shape[0] // 2], dtype=np.float32)
        sl_roi = np.asarray(roi_arr[roi_arr.shape[0] // 2], dtype=np.float32)
    else:
        sl_overview, sl_roi, center_um = _overlap_center_from_sitk(
            overview,
            roi,
            overlap_min,
            overlap_max,
            roi_to_overview=roi_to_overview,
        )
    save_geometry_slice_qc_plot(
        sl_overview=sl_overview,
        sl_roi=sl_roi,
        center_um=center_um,
        output_path=output_path,
        geometry_mode=geometry_mode,
    )


def save_geometry_slice_qc_plot(
    *,
    sl_overview: np.ndarray,
    sl_roi: np.ndarray,
    center_um: tuple[float, float, float],
    output_path: Path,
    geometry_mode: str,
    alignment_metrics: dict[str, object] | None = None,
) -> None:
    """Side-by-side ROI and overview slices at the shared overlap center."""
    _ = (center_um, geometry_mode, alignment_metrics)
    output_path = output_path.expanduser()
    output_path.parent.mkdir(parents=True, exist_ok=True)

    roi_norm = _normalize_panel(sl_roi)
    overview_norm = _normalize_panel(sl_overview)

    fig, axes = plt.subplots(1, 2, figsize=(10, 5))
    axes[0].imshow(roi_norm, cmap="gray")
    axes[0].set_title("ROI")
    axes[0].axis("off")
    axes[1].imshow(overview_norm, cmap="gray")
    axes[1].set_title("Overview (overlap crop)")
    axes[1].axis("off")

    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close(fig)


def save_registration_overlay_qc_plot(
    *,
    overview_crop: sitk.Image | np.ndarray,
    registered_roi: sitk.Image | np.ndarray,
    output_path: Path,
    pair_label: str | None = None,
    z_index: int | None = None,
) -> float:
    """Three-panel post-registration QC: overview crop, registered ROI, R/G overlay.

    Returns mid-plane normalized cross-correlation (overview vs registered ROI).
    """
    output_path = output_path.expanduser()
    output_path.parent.mkdir(parents=True, exist_ok=True)

    ov_arr = (
        sitk.GetArrayFromImage(overview_crop)
        if isinstance(overview_crop, sitk.Image)
        else np.asarray(overview_crop, dtype=np.float32)
    )
    roi_arr = (
        sitk.GetArrayFromImage(registered_roi)
        if isinstance(registered_roi, sitk.Image)
        else np.asarray(registered_roi, dtype=np.float32)
    )
    if ov_arr.shape != roi_arr.shape:
        msg = (
            f"Overview crop and registered ROI shapes must match for overlay QC: "
            f"{ov_arr.shape} vs {roi_arr.shape}"
        )
        raise ValueError(msg)

    z = int(z_index if z_index is not None else ov_arr.shape[0] // 2)
    z = int(np.clip(z, 0, ov_arr.shape[0] - 1))
    sl_overview = np.asarray(ov_arr[z], dtype=np.float32)
    sl_roi = _resample_to_shape(np.asarray(roi_arr[z], dtype=np.float32), sl_overview.shape)

    overview_norm = _normalize_panel(sl_overview)
    roi_norm = _normalize_panel(sl_roi)
    ncc = normalized_cross_correlation(overview_norm, roi_norm)

    overlay = np.zeros((*overview_norm.shape, 3), dtype=np.float32)
    overlay[:, :, 0] = overview_norm
    overlay[:, :, 1] = roi_norm

    fig, axes = plt.subplots(1, 3, figsize=(14, 5))
    axes[0].imshow(overview_norm, cmap="gray")
    axes[0].set_title("Overview crop")
    axes[0].axis("off")
    axes[1].imshow(roi_norm, cmap="gray")
    axes[1].set_title("Registered ROI")
    axes[1].axis("off")
    axes[2].imshow(overlay)
    axes[2].set_title("Overlay R=ov G=roi")
    axes[2].axis("off")

    label = pair_label or "pair"
    fig.suptitle(f"{label} · mid Z · NCC={ncc:.3f}")
    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)
    return ncc
