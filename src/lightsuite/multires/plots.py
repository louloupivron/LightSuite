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
    """Draw axis-aligned FOV boxes for overview, ROI, and overlap."""
    output_path = output_path.expanduser()
    output_path.parent.mkdir(parents=True, exist_ok=True)

    def _xy_box(report: dict[str, object], color: str, label: str) -> None:
        pmin = np.asarray(report["phys_min"], dtype=float)
        pmax = np.asarray(report["phys_max"], dtype=float)
        xs = [pmin[0], pmax[0], pmax[0], pmin[0], pmin[0]]
        ys = [pmin[1], pmin[1], pmax[1], pmax[1], pmin[1]]
        plt.plot(xs, ys, color=color, label=label)

    plt.figure(figsize=(8, 8))
    _xy_box(rep_overview, "C0", str(rep_overview.get("label", "overview")))
    _xy_box(rep_roi, "C1", str(rep_roi.get("label", "roi")))
    ox0, oy0 = float(overlap_min[0]), float(overlap_min[1])
    ox1, oy1 = float(overlap_max[0]), float(overlap_max[1])
    plt.plot([ox0, ox1, ox1, ox0, ox0], [oy0, oy0, oy1, oy1, oy0], "k--", label="overlap")
    plt.gca().set_aspect("equal")
    plt.xlabel("x (µm)")
    plt.ylabel("y (µm)")
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()


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
    output_path = output_path.expanduser()
    output_path.parent.mkdir(parents=True, exist_ok=True)

    sl_overview, sl_roi, center = _overlap_center_from_sitk(
        overview,
        roi,
        overlap_min,
        overlap_max,
        roi_to_overview=roi_to_overview,
    )
    fig, axes = plt.subplots(1, 2, figsize=(10, 5))
    axes[0].imshow(_normalize_panel(sl_overview), cmap="gray")
    axes[0].set_title("Overview")
    axes[0].axis("off")
    axes[1].imshow(_normalize_panel(sl_roi), cmap="gray")
    axes[1].set_title("ROI")
    axes[1].axis("off")
    cx, cy, cz = center
    fig.suptitle(f"{geometry_mode} overlap @ ({cx:.0f}, {cy:.0f}, {cz:.0f}) µm")
    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close(fig)
