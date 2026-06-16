"""QA plots for mesoSPIM geometry checks."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import SimpleITK as sitk

from lightsuite.mesospim.config_models import MesospimGeometryConfig
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
    i, j, k = (int(round(float(t))) for t in idx)
    size = image_geo.GetSize()
    i = int(np.clip(i, 0, size[0] - 1))
    j = int(np.clip(j, 0, size[1] - 1))
    k = int(np.clip(k, 0, size[2] - 1))
    return np.asarray(sitk.GetArrayFromImage(image_geo)[k, :, :])


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


def save_overlap_slice_plot(
    *,
    overview_path: Path,
    roi_path: Path,
    overview_meta: dict,
    roi_meta: dict,
    geometry: MesospimGeometryConfig,
    tiff_remap,
    overlap_min: np.ndarray,
    overlap_max: np.ndarray,
    output_path: Path,
) -> None:
    output_path = output_path.expanduser()
    output_path.parent.mkdir(parents=True, exist_ok=True)

    overview_path = overview_path.expanduser().resolve()
    roi_path = roi_path.expanduser().resolve()

    center_um = 0.5 * (overlap_min + overlap_max)
    cx, cy, cz = float(center_um[0]), float(center_um[1]), float(center_um[2])
    sl_overview = read_tiff_xy_slice_at_physical_um(
        overview_path,
        meta=overview_meta,
        geometry=geometry,
        overview_path=overview_path,
        roi_path=roi_path,
        remap=tiff_remap,
        cx_um=cx,
        cy_um=cy,
        cz_um=cz,
    )
    sl_roi = read_tiff_xy_slice_at_physical_um(
        roi_path,
        meta=roi_meta,
        geometry=geometry,
        overview_path=overview_path,
        roi_path=roi_path,
        remap=tiff_remap,
        cx_um=cx,
        cy_um=cy,
        cz_um=cz,
    )

    fig, axes = plt.subplots(1, 2, figsize=(11, 5))
    axes[0].imshow(sl_overview, cmap="gray")
    axes[0].set_title("Overview XY @ overlap center")
    axes[1].imshow(sl_roi, cmap="gray")
    axes[1].set_title("ROI XY @ overlap center")
    for ax in axes:
        ax.set_aspect("equal")
    fig.suptitle("Tune mesospim.tiff_remap until both panels match")
    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)
