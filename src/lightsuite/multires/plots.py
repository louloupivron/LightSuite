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
    ax.set_title("Lateral (XY)")
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
    ax.set_title("Sagittal (XZ)")
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
) -> None:
    """Side-by-side XY slices at overlap center."""
    sl_overview, sl_roi, center = _overlap_center_from_sitk(
        overview,
        roi,
        overlap_min,
        overlap_max,
        roi_to_overview=roi_to_overview,
    )
    save_geometry_slice_qc_plot(
        sl_overview=sl_overview,
        sl_roi=sl_roi,
        center_um=center,
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
) -> None:
    """Side-by-side XY slices from pre-loaded 2D panels."""
    output_path = output_path.expanduser()
    output_path.parent.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(1, 2, figsize=(10, 5))
    axes[0].imshow(_normalize_panel(sl_overview), cmap="gray")
    axes[0].set_title("Overview")
    axes[0].axis("off")
    axes[1].imshow(_normalize_panel(sl_roi), cmap="gray")
    axes[1].set_title("ROI")
    axes[1].axis("off")
    cx, cy, cz = center_um
    fig.suptitle(f"{geometry_mode} overlap @ ({cx:.0f}, {cy:.0f}, {cz:.0f}) µm")
    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close(fig)
