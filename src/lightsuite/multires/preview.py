"""Export small TIFF crops for visual alignment checks."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Literal

import numpy as np
import tifffile

from lightsuite.multires.config_models import MultiresPipelineConfig
from lightsuite.multires.manifest import load_pair_manifest
from lightsuite.multires.spec_geometry import (
    crop_index_range_from_physical_box,
    overlap_physical_bounds_from_specs,
    physical_to_continuous_index_xyz,
)
from lightsuite.multires.volume import _normalize_tiff_array, _resolve_volume_path, _sorted_plane_files


def _load_xy_crop_from_plane(
    volume_path: Path,
    *,
    z_index: int,
    start_xyz: list[int],
    crop_size_xyz: list[int],
) -> np.ndarray:
    ix0, iy0, _iz0 = start_xyz
    sx, sy, _sz = crop_size_xyz
    ix1 = ix0 + sx
    iy1 = iy0 + sy

    if volume_path.is_dir():
        planes = _sorted_plane_files(volume_path)
        plane = _normalize_tiff_array(
            np.asarray(tifffile.imread(str(planes[z_index]))),
            planes[z_index],
        )
        if plane.ndim == 3:
            plane = plane[0]
    else:
        with tifffile.TiffFile(volume_path) as tf:
            series = tf.series[0]
            if len(series.shape) == 3 and series.shape[0] > 1:
                plane = np.asarray(series.asarray(key=z_index))
            else:
                plane = np.asarray(tf.pages[z_index].asarray())
        plane = _normalize_tiff_array(plane, volume_path)
        if plane.ndim == 3:
            plane = plane[0]

    return np.asarray(plane[iy0:iy1, ix0:ix1], dtype=np.float32)


def _expand_crop_box(
    start_xyz: list[int],
    crop_size_xyz: list[int],
    *,
    max_size_xyz: tuple[int, int, int],
    margin_xy_vox: int,
) -> tuple[list[int], list[int]]:
    ix0, iy0, iz0 = start_xyz
    sx, sy, sz = crop_size_xyz
    nx, ny, nz = max_size_xyz
    ix0 = max(0, ix0 - margin_xy_vox)
    iy0 = max(0, iy0 - margin_xy_vox)
    ix1 = min(nx, ix0 + sx + 2 * margin_xy_vox)
    iy1 = min(ny, iy0 + sy + 2 * margin_xy_vox)
    return [ix0, iy0, iz0], [ix1 - ix0, iy1 - iy0, sz]


def _z_indices_around_center(nz: int, center_z: int, n_slices: int) -> list[int]:
    half = max(0, n_slices // 2)
    indices = list(range(center_z - half, center_z - half + n_slices))
    return [int(np.clip(index, 0, nz - 1)) for index in indices]


def _write_preview_stack(
    path: Path,
    stack_zyx: np.ndarray,
    *,
    spacing_xy_um: float,
    spacing_z_um: float,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tifffile.imwrite(
        path,
        stack_zyx.astype(np.float32, copy=False),
        imagej=True,
        resolution=(1.0 / spacing_xy_um, 1.0 / spacing_xy_um),
        metadata={"spacing": spacing_z_um, "unit": "um"},
        compression="zlib",
    )


def _write_preview_image(
    path: Path,
    image_yx: np.ndarray,
    *,
    spacing_xy_um: float,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tifffile.imwrite(
        path,
        image_yx.astype(np.float32, copy=False),
        imagej=True,
        resolution=(1.0 / spacing_xy_um, 1.0 / spacing_xy_um),
        metadata={"unit": "um"},
        compression="zlib",
    )


def _load_xy_max_projection(
    volume_path: Path,
    *,
    start_xyz: list[int],
    crop_size_xyz: list[int],
) -> np.ndarray:
    """Z-max projection over the overlap crop box."""
    _ix0, _iy0, iz0 = start_xyz
    _sx, _sy, sz = crop_size_xyz
    max_proj: np.ndarray | None = None
    for z_index in range(iz0, iz0 + sz):
        plane = _load_xy_crop_from_plane(
            volume_path,
            z_index=z_index,
            start_xyz=start_xyz,
            crop_size_xyz=crop_size_xyz,
        )
        if max_proj is None:
            max_proj = plane.copy()
        else:
            np.maximum(max_proj, plane, out=max_proj)
    if max_proj is None:
        msg = "Empty Z range for max projection"
        raise ValueError(msg)
    return max_proj


PreviewProjection = Literal["slices", "max"]


def export_alignment_preview_crops(
    cfg: MultiresPipelineConfig,
    *,
    output_dir: Path | None = None,
    n_slices: int = 5,
    margin_um: float = 100.0,
    projection: PreviewProjection = "slices",
) -> Path:
    """Export overview / ROI overlap crops for visual QC.

    With ``projection="slices"`` (default), writes small Z stacks around the overlap center.
    With ``projection="max"``, writes a single XY max-intensity projection over the full
    overlap Z extent instead of subvolumes.
    """
    manifest_path = cfg.multires.pair_manifest
    manifest = load_pair_manifest(manifest_path)
    manifest_dir = manifest_path.parent
    pair_label = manifest.pair_label

    overlap_min, overlap_max = overlap_physical_bounds_from_specs(
        manifest.overview,
        manifest.roi,
        margin_um=cfg.multires.registration.overlap_margin_um,
    )
    overview_start, overview_size = crop_index_range_from_physical_box(
        manifest.overview,
        overlap_min,
        overlap_max,
    )
    roi_start, roi_size = crop_index_range_from_physical_box(
        manifest.roi,
        overlap_min,
        overlap_max,
    )

    margin_xy_vox_overview = int(round(margin_um / manifest.overview.spacing_um[0]))
    margin_xy_vox_roi = int(round(margin_um / manifest.roi.spacing_um[0]))
    _, _, onz = (int(v) for v in manifest.overview.shape_zyx)
    onx, ony, _ = (
        int(manifest.overview.shape_zyx[2]),
        int(manifest.overview.shape_zyx[1]),
        int(manifest.overview.shape_zyx[0]),
    )
    _, _, rnz = (int(v) for v in manifest.roi.shape_zyx)
    rnx, rny, _ = (
        int(manifest.roi.shape_zyx[2]),
        int(manifest.roi.shape_zyx[1]),
        int(manifest.roi.shape_zyx[0]),
    )

    overview_start, overview_size = _expand_crop_box(
        overview_start,
        overview_size,
        max_size_xyz=(onx, ony, onz),
        margin_xy_vox=margin_xy_vox_overview,
    )
    roi_start, roi_size = _expand_crop_box(
        roi_start,
        roi_size,
        max_size_xyz=(rnx, rny, rnz),
        margin_xy_vox=margin_xy_vox_roi,
    )

    center_um = 0.5 * (overlap_min + overlap_max)
    overview_center_z = int(
        round(
            physical_to_continuous_index_xyz(manifest.overview, center_um)[2]
        )
    )
    roi_center_z = int(
        round(physical_to_continuous_index_xyz(manifest.roi, center_um)[2])
    )
    overview_path = _resolve_volume_path(manifest.overview, manifest_dir)
    roi_path = _resolve_volume_path(manifest.roi, manifest_dir)

    out_dir = output_dir or (cfg.sample.save_path / "geometry" / "alignment_preview" / pair_label)
    out_dir = out_dir.expanduser()
    out_dir.mkdir(parents=True, exist_ok=True)

    if projection == "max":
        overview_image = _load_xy_max_projection(
            overview_path,
            start_xyz=overview_start,
            crop_size_xyz=overview_size,
        )
        roi_image = _load_xy_max_projection(
            roi_path,
            start_xyz=roi_start,
            crop_size_xyz=roi_size,
        )
        overview_tiff = out_dir / f"{pair_label}_overview_overlap_maxproj.tif"
        roi_tiff = out_dir / f"{pair_label}_roi_overlap_maxproj.tif"
        _write_preview_image(
            overview_tiff,
            overview_image,
            spacing_xy_um=float(manifest.overview.spacing_um[0]),
        )
        _write_preview_image(
            roi_tiff,
            roi_image,
            spacing_xy_um=float(manifest.roi.spacing_um[0]),
        )
        overview_z_report: list[int] | dict[str, int] = {
            "start": overview_start[2],
            "end": overview_start[2] + overview_size[2] - 1,
            "count": overview_size[2],
        }
        roi_z_report: list[int] | dict[str, int] = {
            "start": roi_start[2],
            "end": roi_start[2] + roi_size[2] - 1,
            "count": roi_size[2],
        }
        overview_kind = "max_projection"
        roi_kind = "max_projection"
    else:
        overview_z = _z_indices_around_center(onz, overview_center_z, n_slices)
        roi_z = _z_indices_around_center(rnz, roi_center_z, n_slices)
        overview_z_report = overview_z
        roi_z_report = roi_z
        overview_stack = np.zeros(
            (len(overview_z), overview_size[1], overview_size[0]),
            dtype=np.float32,
        )
        roi_stack = np.zeros((len(roi_z), roi_size[1], roi_size[0]), dtype=np.float32)
        for out_index, z_index in enumerate(overview_z):
            overview_stack[out_index] = _load_xy_crop_from_plane(
                overview_path,
                z_index=z_index,
                start_xyz=overview_start,
                crop_size_xyz=overview_size,
            )
        for out_index, z_index in enumerate(roi_z):
            roi_stack[out_index] = _load_xy_crop_from_plane(
                roi_path,
                z_index=z_index,
                start_xyz=roi_start,
                crop_size_xyz=roi_size,
            )
        overview_tiff = out_dir / f"{pair_label}_overview_overlap_crop.tif"
        roi_tiff = out_dir / f"{pair_label}_roi_overlap_crop.tif"
        _write_preview_stack(
            overview_tiff,
            overview_stack,
            spacing_xy_um=float(manifest.overview.spacing_um[0]),
            spacing_z_um=float(manifest.overview.spacing_um[2]),
        )
        _write_preview_stack(
            roi_tiff,
            roi_stack,
            spacing_xy_um=float(manifest.roi.spacing_um[0]),
            spacing_z_um=float(manifest.roi.spacing_um[2]),
        )
        overview_kind = "slice_stack"
        roi_kind = "slice_stack"

    from lightsuite.multires.plots import save_geometry_slice_qc_plot

    qc_plot_path = out_dir / f"{pair_label}_overlap_{projection}_qc.png"
    save_geometry_slice_qc_plot(
        sl_overview=overview_image if projection == "max" else overview_stack[len(overview_stack) // 2],
        sl_roi=roi_image if projection == "max" else roi_stack[len(roi_stack) // 2],
        center_um=tuple(float(v) for v in center_um),
        output_path=qc_plot_path,
        geometry_mode=f"metadata ({projection})",
    )

    report = {
        "pair_label": pair_label,
        "pair_manifest_path": str(manifest_path),
        "projection": projection,
        "overlap_box_um": [overlap_min.tolist(), overlap_max.tolist()],
        "overlap_center_um": center_um.tolist(),
        "overview": {
            "volume_path": str(overview_path),
            "crop_start_xyz": overview_start,
            "crop_size_xyz": overview_size,
            "z_indices": overview_z_report,
            "kind": overview_kind,
            "output_tiff": str(overview_tiff),
        },
        "roi": {
            "volume_path": str(roi_path),
            "crop_start_xyz": roi_start,
            "crop_size_xyz": roi_size,
            "z_indices": roi_z_report,
            "kind": roi_kind,
            "output_tiff": str(roi_tiff),
        },
        "overlap_qc_png": str(qc_plot_path),
    }
    report_path = out_dir / f"{pair_label}_alignment_preview.json"
    report_path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    label = "max projection" if projection == "max" else "crop stack"
    print(f"Overview {label}: {overview_tiff}", flush=True)
    print(f"ROI {label}: {roi_tiff}", flush=True)
    print(f"Overlap QC plot: {qc_plot_path}", flush=True)
    print(f"Preview metadata: {report_path}", flush=True)
    return out_dir


__all__ = ["PreviewProjection", "export_alignment_preview_crops"]
