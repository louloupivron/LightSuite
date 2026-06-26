"""Export registered spinal cord volumes (generateRegisteredCordVolume.m port)."""

from __future__ import annotations

import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import tifffile
from rich.console import Console

from lightsuite.atlas.fiederling import load_fiederling_atlas_volumes, upsample_to_fiederling_native
from lightsuite.config.models import SpinalCordPipelineConfig
from lightsuite.export.cord_registered import (
    REGISTERED_ANNOTATION_FILENAME,
    REGISTERED_TEMPLATE_FILENAME,
    CordRegisteredInspectPaths,
    ensure_registered_annotation_volume,
    export_layout_from_native,
    load_native_template_export_layout,
    warp_output_to_uint16,
)
from lightsuite.io.cord_volume import load_registration_volumes
from lightsuite.preprocess.cord_checkpoint import (
    CordRegOptsCheckpoint,
    CordTransformParamsCheckpoint,
)
from lightsuite.registration.cord_affine import (
    build_cord_z_transinit,
    warp_cord_straightvol_to_atlas,
)
from lightsuite.registration.cord_paths import (
    cord_affine_transform_path,
    cord_save_path,
    cord_work_dir,
)
from lightsuite.registration.elastix.runner import run_transformix
from lightsuite.registration.straightening import load_slicetforms, transform_cord_images_slices

console = Console()


@dataclass(frozen=True)
class CordExportResult:
    output_dir: Path
    channel_paths: list[Path]


def _load_transform_params(save_path: Path) -> CordTransformParamsCheckpoint:
    json_path = save_path / "transform_params.json"
    if not json_path.is_file():
        msg = f"Missing {json_path}. Run 'lightsuite spinal register' first."
        raise FileNotFoundError(msg)
    return CordTransformParamsCheckpoint.load(json_path)


def export_registered_cord_volumes(config: SpinalCordPipelineConfig) -> CordExportResult:
    """Register all channels and save slice-style TIFF outputs."""
    save_path = cord_save_path(config)
    checkpoint = CordRegOptsCheckpoint.load(save_path / "regopts.json")
    transform_params = _load_transform_params(save_path)

    console.print("Loading cached registration-grid volumes...")
    finvol = load_registration_volumes(checkpoint)

    perm = [p - 1 for p in transform_params.how_to_perm]
    finvol = np.transpose(finvol, perm + [3])
    yrange = transform_params.samp_ikeepy
    xrange = transform_params.samp_ikeepx
    zrange = transform_params.samp_ikeeplong
    finvol = finvol[yrange[0] - 1 : yrange[1], xrange[0] - 1 : xrange[1], zrange[0] - 1 : zrange[1], :]

    tforms = load_slicetforms(transform_params.slicetforms_path)
    sizetv = tuple(transform_params.atlassize[:2])
    raout = sizetv
    atlas_shape = tuple(transform_params.atlassize)
    bspline_path = Path(transform_params.tform_bspline_samp20um_to_atlas_20um_px)
    elastix_affine_path = cord_affine_transform_path(config)
    spacing_mm = config.registration.resolution_um * 1e-3
    nslices = checkpoint.ikeeprange[1] - checkpoint.ikeeprange[0] + 1
    transinit = build_cord_z_transinit(nslices, atlas_shape[2])

    atlas_volumes = load_fiederling_atlas_volumes(config.atlas)
    native_template = atlas_volumes.template

    output_dir = save_path / "volume_registered"
    output_dir.mkdir(parents=True, exist_ok=True)
    channel_paths: list[Path] = []
    t0 = time.perf_counter()

    for ich in range(finvol.shape[3]):
        currvol = transform_cord_images_slices(finvol[:, :, :, ich], tforms, raout)
        volumereg = run_transformix(
            moving_volume=currvol.astype(np.float32),
            transform_path=bspline_path,
            output_dir=cord_work_dir(config, "transformix", "export", f"ch{ich + 1}", "bspline"),
            spacing_mm=spacing_mm,
            nearest=False,
        )
        registered = warp_output_to_uint16(
            warp_cord_straightvol_to_atlas(
                volumereg,
                transinit=transinit,
                elastix_affine_path=elastix_affine_path,
                atlas_shape=atlas_shape,
                spacing_mm=spacing_mm,
                work_dir=cord_work_dir(config, "transformix", "export", f"ch{ich + 1}", "affine"),
                nearest=False,
            )
        )
        if transform_params.tofliprc:
            registered = np.flip(registered, axis=2)
        registered_native = upsample_to_fiederling_native(registered, native_template)
        out_path = output_dir / f"chan{ich + 1:02d}_channel{ich + 1}.tiff"
        export_stack = export_layout_from_native(registered_native)
        tifffile.imwrite(out_path, export_stack, imagej=True)
        channel_paths.append(out_path)
        console.print(f"Channel {ich + 1}/{finvol.shape[3]} exported in {time.perf_counter() - t0:.2f}s")

    annotation_paths = CordRegisteredInspectPaths(
        volume_registered_dir=output_dir,
        annotation_path=output_dir / REGISTERED_ANNOTATION_FILENAME,
        registered_channels={idx + 1: path for idx, path in enumerate(channel_paths)},
    )
    ensure_registered_annotation_volume(config, paths=annotation_paths, force_recompute=True)

    template_path = output_dir / REGISTERED_TEMPLATE_FILENAME
    template_stack = load_native_template_export_layout(config)
    tifffile.imwrite(template_path, template_stack.astype(np.uint16), imagej=True)

    return CordExportResult(output_dir=output_dir, channel_paths=channel_paths)
