"""QA plots for mesoSPIM geometry checks."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import SimpleITK as sitk

from lightsuite.mesospim.config_models import MesospimGeometryConfig, MesospimTiffRemapConfig
from lightsuite.mesospim.io import read_tiff_xy_slice_at_physical_um


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


def _overlap_center_from_tiff(
    *,
    overview_path: Path,
    roi_path: Path,
    overview_meta: dict,
    roi_meta: dict,
    geometry: MesospimGeometryConfig,
    tiff_remap: MesospimTiffRemapConfig,
    overlap_min: np.ndarray,
    overlap_max: np.ndarray,
    roi_to_overview: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray, tuple[float, float, float]]:
    """Load one overview / ROI XY plane per overlap box from disk (no full-volume read)."""
    from lightsuite.mesospim.geometry import transform_physical_points

    overview_path = overview_path.expanduser().resolve()
    roi_path = roi_path.expanduser().resolve()
    center_um = 0.5 * (overlap_min + overlap_max)
    cx, cy, cz = (float(center_um[0]), float(center_um[1]), float(center_um[2]))
    read_kwargs = {
        "geometry": geometry,
        "overview_path": overview_path,
        "roi_path": roi_path,
        "remap": tiff_remap,
    }
    sl_overview = read_tiff_xy_slice_at_physical_um(
        overview_path,
        meta=overview_meta,
        cx_um=cx,
        cy_um=cy,
        cz_um=cz,
        **read_kwargs,
    )
    if roi_to_overview is None:
        rcx, rcy, rcz = cx, cy, cz
    else:
        roi_center = transform_physical_points(
            np.array([[cx, cy, cz]], dtype=float),
            np.linalg.inv(roi_to_overview),
        )[0]
        rcx, rcy, rcz = float(roi_center[0]), float(roi_center[1]), float(roi_center[2])
    sl_roi = read_tiff_xy_slice_at_physical_um(
        roi_path,
        meta=roi_meta,
        cx_um=rcx,
        cy_um=rcy,
        cz_um=rcz,
        **read_kwargs,
    )
    return sl_overview, sl_roi, (cx, cy, cz)


def _overlap_center_from_sitk(
    overview: sitk.Image,
    roi: sitk.Image,
    overlap_min: np.ndarray,
    overlap_max: np.ndarray,
    *,
    roi_to_overview: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray, tuple[float, float, float]]:
    """Extract overview / ROI XY slices from in-memory volumes (tests / small stacks)."""
    from lightsuite.mesospim.geometry import transform_physical_points

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


def save_geometry_overlap_qc_plot(
    *,
    metadata_overlap_min: np.ndarray,
    metadata_overlap_max: np.ndarray,
    output_path: Path,
    geometry_mode: str,
    hybrid_overlap_min: np.ndarray | None = None,
    hybrid_overlap_max: np.ndarray | None = None,
    roi_to_overview: np.ndarray | None = None,
    overview_path: Path | None = None,
    roi_path: Path | None = None,
    overview_meta: dict | None = None,
    roi_meta: dict | None = None,
    geometry: MesospimGeometryConfig | None = None,
    tiff_remap: MesospimTiffRemapConfig | None = None,
    overview: sitk.Image | None = None,
    roi: sitk.Image | None = None,
) -> None:
    """Compare XY slices at overlap centers for metadata vs hybrid / landmark placement."""
    output_path = output_path.expanduser()
    output_path.parent.mkdir(parents=True, exist_ok=True)

    use_tiff = (
        overview_path is not None
        and roi_path is not None
        and overview_meta is not None
        and roi_meta is not None
        and geometry is not None
        and tiff_remap is not None
    )
    if use_tiff:
        extract = lambda o_min, o_max, tform=None: _overlap_center_from_tiff(
            overview_path=overview_path,
            roi_path=roi_path,
            overview_meta=overview_meta,
            roi_meta=roi_meta,
            geometry=geometry,
            tiff_remap=tiff_remap,
            overlap_min=o_min,
            overlap_max=o_max,
            roi_to_overview=tform,
        )
    elif overview is not None and roi is not None:
        extract = lambda o_min, o_max, tform=None: _overlap_center_from_sitk(
            overview,
            roi,
            o_min,
            o_max,
            roi_to_overview=tform,
        )
    else:
        msg = "save_geometry_overlap_qc_plot requires TIFF paths/meta or in-memory overview/roi"
        raise ValueError(msg)

    meta_overview, meta_roi, meta_center = extract(metadata_overlap_min, metadata_overlap_max)

    show_hybrid = hybrid_overlap_min is not None and hybrid_overlap_max is not None
    if show_hybrid:
        assert roi_to_overview is not None
        hybrid_overview, hybrid_roi, hybrid_center = extract(
            hybrid_overlap_min,
            hybrid_overlap_max,
            roi_to_overview,
        )

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    meta_mode_label = "Metadata (stage frame)"
    axes[0, 0].imshow(_normalize_panel(meta_overview), cmap="gray")
    axes[0, 0].set_title(
        f"Overview — {meta_mode_label}\nZ≈{meta_center[2]:.0f} µm",
    )
    axes[1, 0].imshow(_normalize_panel(meta_roi), cmap="gray")
    axes[1, 0].set_title(
        f"ROI — {meta_mode_label}\nZ≈{meta_center[2]:.0f} µm",
    )

    if show_hybrid:
        hybrid_mode_label = "Hybrid (landmarks)" if geometry_mode == "hybrid" else "Landmarks"
        axes[0, 1].imshow(_normalize_panel(hybrid_overview), cmap="gray")
        axes[0, 1].set_title(
            f"Overview — {hybrid_mode_label}\nZ≈{hybrid_center[2]:.0f} µm",
        )
        axes[1, 1].imshow(_normalize_panel(hybrid_roi), cmap="gray")
        axes[1, 1].set_title(
            f"ROI — {hybrid_mode_label}\nZ≈{hybrid_center[2]:.0f} µm",
        )
        suptitle = (
            "Geometry overlap QC — compare metadata stage placement with "
            f"{geometry_mode} landmark placement"
        )
    else:
        for ax in (axes[0, 1], axes[1, 1]):
            ax.axis("off")
            ax.text(
                0.5,
                0.5,
                "Hybrid / landmark comparison\nnot computed\n(geometry_mode=metadata)",
                ha="center",
                va="center",
                transform=ax.transAxes,
            )
        suptitle = "Geometry overlap QC — metadata (stage frame) placement"

    for ax in axes.ravel():
        if ax.images or ax.texts:
            ax.set_aspect("equal")
            ax.set_xticks([])
            ax.set_yticks([])

    fig.suptitle(suptitle)
    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)


def save_fov_overlap_plot(
    *,
    rep_overview: dict[str, object],
    rep_roi: dict[str, object],
    overlap_min: np.ndarray,
    overlap_max: np.ndarray,
    output_path: Path,
    title: str,
) -> None:
    output_path = output_path.expanduser()
    output_path.parent.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    def rect_xy(ax, pmin, pmax, **kwargs):
        w, h = pmax[0] - pmin[0], pmax[1] - pmin[1]
        ax.add_patch(plt.Rectangle((pmin[0], pmin[1]), w, h, fill=False, **kwargs))

    def rect_xz(ax, pmin, pmax, **kwargs):
        w, h = pmax[0] - pmin[0], pmax[2] - pmin[2]
        ax.add_patch(plt.Rectangle((pmin[0], pmin[2]), w, h, fill=False, **kwargs))

    ax = axes[0]
    rect_xy(
        ax,
        rep_overview["phys_min"],
        rep_overview["phys_max"],
        color="C0",
        lw=2,
        label="overview",
    )
    rect_xy(ax, rep_roi["phys_min"], rep_roi["phys_max"], color="C1", lw=2, label="ROI")
    rect_xy(ax, overlap_min, overlap_max, color="C2", lw=2, ls="--", label="overlap")
    ax.scatter(*rep_overview["phys_center"][:2], c="C0", s=40, zorder=5)
    ax.scatter(*rep_roi["phys_center"][:2], c="C1", s=40, zorder=5)
    ax.set_xlabel("physical X (µm)")
    ax.set_ylabel("physical Y (µm)")
    ax.set_title("Lateral (XY)")
    ax.set_aspect("equal")
    ax.legend(loc="best")

    ax = axes[1]
    rect_xz(
        ax,
        rep_overview["phys_min"],
        rep_overview["phys_max"],
        color="C0",
        lw=2,
        label="overview",
    )
    rect_xz(ax, rep_roi["phys_min"], rep_roi["phys_max"], color="C1", lw=2, label="ROI")
    rect_xz(ax, overlap_min, overlap_max, color="C2", lw=2, ls="--", label="overlap")
    ax.scatter(
        rep_overview["phys_center"][0],
        rep_overview["phys_center"][2],
        c="C0",
        s=40,
        zorder=5,
    )
    ax.scatter(rep_roi["phys_center"][0], rep_roi["phys_center"][2], c="C1", s=40, zorder=5)
    ax.set_xlabel("physical X (µm)")
    ax.set_ylabel("physical Z (µm)")
    ax.set_title("Sagittal (XZ)")
    ax.legend(loc="best")

    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)
