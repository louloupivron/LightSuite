#!/usr/bin/env python3
"""Warp Arivis segmentation from 2.5x native ROI space to 0.8x overview native space."""

from __future__ import annotations

import json
import re
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import SimpleITK as sitk
import tifffile

from lightsuite.multires.manifest import load_pair_manifest
from lightsuite.multires.spec_geometry import (
    index_xyz_to_physical,
    physical_to_continuous_index_xyz,
)
from lightsuite.multires.volume import apply_manifest_geometry
from lightsuite.registration.elastix.points import write_landmark_file

REGOPTS = Path(
    "/media/gbm/NVME2/ALICe-pipelines-data/Marianna_CMU/multiresolution_results/multires_regopts.json"
)
POINTS_IN = Path(
    "/media/gbm/NVME2/ALICe-pipelines-data/Marianna_CMU/segmentation/converted/arivis_2p5x_points.csv"
)
MASK_IN = Path(
    "/media/gbm/NVME2/ALICe-pipelines-data/Marianna_CMU/segmentation/converted/arivis_2p5x_mask.tif"
)
OUT_DIR = POINTS_IN.parent
POINTS_OUT = OUT_DIR / "arivis_2p5x_points_0p8x.csv"
MASK_OUT = OUT_DIR / "arivis_2p5x_mask_0p8x.tif"
WARP_TMP = OUT_DIR / "_warp_tmp"


def crop_reference_image(
    overview_spec,
    crop_start_index: list[int],
    crop_size_xyz: list[int],
) -> sitk.Image:
    """Geometry-only reference grid for elastix (matches registration prepare step)."""
    sx, sy, sz = (int(v) for v in crop_size_xyz)
    image = sitk.Image([sx, sy, sz], sitk.sitkUInt8)
    origin = index_xyz_to_physical(overview_spec, tuple(float(v) for v in crop_start_index))
    image.SetSpacing(tuple(float(v) for v in overview_spec.spacing_um))
    image.SetOrigin(tuple(float(v) for v in origin))
    image.SetDirection(tuple(float(v) for v in overview_spec.direction))
    return image


def write_embedded_crop_canvas_uint8(
    overview_spec,
    crop: sitk.Image,
    crop_start_index: list[int],
    output_path: Path,
) -> None:
    """Embed a uint8 crop into a full-overview BigTIFF (plane-by-plane)."""
    output_path = output_path.expanduser()
    output_path.parent.mkdir(parents=True, exist_ok=True)

    nz, ny, nx = (int(v) for v in overview_spec.shape_zyx)
    cx, cy, cz = (int(v) for v in crop.GetSize())
    ix0, iy0, iz0 = (int(v) for v in crop_start_index)
    crop_arr = np.asarray(sitk.GetArrayViewFromImage(crop), dtype=np.uint8)
    sx, sy, sz = (float(v) for v in overview_spec.spacing_um)

    with tifffile.TiffWriter(output_path, bigtiff=True) as tif:
        for iz in range(nz):
            plane = np.zeros((ny, nx), dtype=np.uint8)
            if iz0 <= iz < iz0 + cz:
                y1 = min(iy0 + cy, ny)
                x1 = min(ix0 + cx, nx)
                yy0 = max(iy0, 0)
                xx0 = max(ix0, 0)
                src_y0 = yy0 - iy0
                src_x0 = xx0 - ix0
                plane[yy0:y1, xx0:x1] = crop_arr[
                    iz - iz0,
                    src_y0 : src_y0 + (y1 - yy0),
                    src_x0 : src_x0 + (x1 - xx0),
                ]
            metadata = None
            if iz == 0:
                metadata = {
                    "axes": "ZYX",
                    "spacing": sz,
                    "unit": "um",
                    "loop": False,
                }
            tif.write(
                plane,
                compression="zlib",
                photometric="minisblack",
                metadata=metadata,
                resolution=(1.0 / sx, 1.0 / sy),
            )


def _parse_transformix_output(path: Path) -> np.ndarray:
    mapped: list[list[float]] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        match = re.search(r"OutputPoint\s*=\s*\[([^\]]+)\]", line)
        if match is None:
            continue
        mapped.append([float(v) for v in match.group(1).split()])
    return np.asarray(mapped, dtype=float)


def warp_points_physical(
    points_xyz_1based: np.ndarray,
    *,
    roi_spec,
    overview_spec,
    transform_path: Path,
    tmp_dir: Path,
) -> np.ndarray:
    """Map LightSuite 1-based XYZ indices through elastix (ROI -> overview physical)."""
    idx0 = np.asarray(points_xyz_1based, dtype=float) - 1.0
    phys_roi = np.vstack([index_xyz_to_physical(roi_spec, row) for row in idx0])

    tmp_dir.mkdir(parents=True, exist_ok=True)
    for pattern in ("inputPoints.txt", "outputpoints.txt", "transformix.log"):
        p = tmp_dir / pattern
        if p.is_file():
            p.unlink()

    input_path = tmp_dir / "inputPoints.txt"
    write_landmark_file(input_path, phys_roi)

    if shutil.which("transformix") is None:
        msg = "transformix not found on PATH"
        raise RuntimeError(msg)

    cmd = [
        "transformix",
        "-def",
        str(input_path),
        "-out",
        str(tmp_dir),
        "-tp",
        str(transform_path),
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if proc.returncode != 0:
        msg = f"transformix failed:\n{proc.stdout}\n{proc.stderr}"
        raise RuntimeError(msg)

    out_phys = _parse_transformix_output(tmp_dir / "outputpoints.txt")
    if out_phys.shape[0] != points_xyz_1based.shape[0]:
        msg = f"Expected {points_xyz_1based.shape[0]} points, got {out_phys.shape[0]}"
        raise RuntimeError(msg)

    out_idx0 = np.vstack(
        [physical_to_continuous_index_xyz(overview_spec, row) for row in out_phys]
    )
    return out_idx0 + 1.0


def apply_elastix_nearest(
    moving: sitk.Image,
    transform_paths: list[Path],
    *,
    reference: sitk.Image,
) -> sitk.Image:
    """Apply saved elastix transforms with nearest-neighbour interpolation."""
    import itk

    from lightsuite.multires.volume import sitk_to_itk

    parameter_object = itk.ParameterObject.New()
    for path in transform_paths:
        parameter_object.AddParameterFile(str(path))

    last = parameter_object.GetNumberOfParameterMaps() - 1
    parameter_object.SetParameter(last, "Size", [str(int(v)) for v in reference.GetSize()])
    parameter_object.SetParameter(last, "Spacing", [str(float(v)) for v in reference.GetSpacing()])
    parameter_object.SetParameter(last, "Origin", [str(float(v)) for v in reference.GetOrigin()])
    parameter_object.SetParameter(
        last,
        "Direction",
        [str(float(v)) for v in reference.GetDirection()],
    )
    parameter_object.SetParameter(last, "ResampleInterpolator", "FinalNearestNeighborInterpolator")

    result_itk = itk.transformix_filter(
        sitk_to_itk(moving),
        transform_parameter_object=parameter_object,
    )
    result_sitk = sitk.GetImageFromArray(
        itk.GetArrayFromImage(result_itk).astype(np.uint8, copy=False)
    )
    result_sitk.CopyInformation(reference)
    return result_sitk


def warp_mask(
    mask_path: Path,
    *,
    roi_spec,
    overview_spec,
    transform_paths: list[Path],
    crop_size_xyz: list[int],
    crop_start_index: list[int],
    output_path: Path,
) -> dict[str, int]:
    print(f"Loading mask {mask_path} (memmap)...")
    arr = tifffile.memmap(mask_path)
    if tuple(arr.shape) != tuple(roi_spec.shape_zyx):
        msg = f"Mask shape {arr.shape} != ROI manifest {roi_spec.shape_zyx}"
        raise ValueError(msg)

    moving = sitk.GetImageFromArray(np.asarray(arr > 0, dtype=np.uint8))
    moving = apply_manifest_geometry(moving, roi_spec)

    reference = crop_reference_image(overview_spec, crop_start_index, crop_size_xyz)
    print("Resampling mask onto overview crop (nearest neighbour)...")
    warped_crop = apply_elastix_nearest(moving, transform_paths, reference=reference)

    print(f"Embedding crop into full 0.8x overview -> {output_path}")
    write_embedded_crop_canvas_uint8(overview_spec, warped_crop, crop_start_index, output_path)

    crop_nz = int(np.count_nonzero(sitk.GetArrayFromImage(warped_crop)))
    return {
        "overview_shape_zyx": tuple(int(v) for v in overview_spec.shape_zyx),
        "labels_in_crop": crop_nz,
    }


def main() -> None:
    regopts = json.loads(REGOPTS.read_text(encoding="utf-8"))
    manifest = load_pair_manifest(regopts["pair_manifest_path"])
    roi_spec = manifest.roi
    overview_spec = manifest.overview
    transform_paths = [Path(p) for p in regopts["transform_paths"]]
    transform_final = transform_paths[-1]
    crop_start = [int(v) for v in regopts["crop_start_index"]]
    crop_ref = Path(regopts["cropped_overview_path"])
    crop_size_xyz = list(sitk.ReadImage(str(crop_ref)).GetSize())

    print("=== Warp points 2.5x -> 0.8x ===")
    df = pd.read_csv(POINTS_IN)
    pts_in = df[["x", "y", "z"]].to_numpy(dtype=float)
    pts_out = warp_points_physical(
        pts_in,
        roi_spec=roi_spec,
        overview_spec=overview_spec,
        transform_path=transform_final,
        tmp_dir=WARP_TMP / "points",
    )

    nx, ny, nz = (
        int(overview_spec.shape_zyx[2]),
        int(overview_spec.shape_zyx[1]),
        int(overview_spec.shape_zyx[0]),
    )
    in_bounds = (
        (pts_out[:, 0] >= 1)
        & (pts_out[:, 0] <= nx)
        & (pts_out[:, 1] >= 1)
        & (pts_out[:, 1] <= ny)
        & (pts_out[:, 2] >= 1)
        & (pts_out[:, 2] <= nz)
    )

    out_df = df.copy()
    out_df["x"] = pts_out[:, 0]
    out_df["y"] = pts_out[:, 1]
    out_df["z"] = pts_out[:, 2]
    out_df.to_csv(POINTS_OUT, index=False)
    print(f"Wrote {len(out_df)} points -> {POINTS_OUT}")
    print(f"In-bounds on 0.8x grid ({nx}x{ny}x{nz} XYZ): {in_bounds.sum()}/{len(in_bounds)}")

    print("\n=== Warp mask 2.5x -> 0.8x ===")
    mask_stats = warp_mask(
        MASK_IN,
        roi_spec=roi_spec,
        overview_spec=overview_spec,
        transform_paths=transform_paths,
        crop_size_xyz=crop_size_xyz,
        crop_start_index=crop_start,
        output_path=MASK_OUT,
    )
    print(f"Mask stats: {mask_stats}")
    print("Done.")


if __name__ == "__main__":
    main()
