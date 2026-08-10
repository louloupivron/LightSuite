"""Spinal cord sample-space export (atlas warped onto straightened registration grid)."""

from __future__ import annotations

import json
import re
import time
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import tifffile
from rich.console import Console
from scipy import ndimage

from lightsuite.analysis.cord_hemisphere import HEMISPHERE_IN_SAMPLE, load_fiederling_hemisphere_native
from lightsuite.atlas.fiederling import load_fiederling_atlas_volumes, resolve_fiederling_paths
from lightsuite.config.models import SpinalCordPipelineConfig
from lightsuite.export.cord_registered import warp_output_to_uint16
from lightsuite.io.cord_volume import load_registration_volumes
from lightsuite.io.tiff_write import save_registration_volume
from lightsuite.preprocess.cord_checkpoint import (
    CordRegOptsCheckpoint,
    CordTransformParamsCheckpoint,
)
from lightsuite.registration.cord_affine import warp_cord_atlas_to_straightvol
from lightsuite.registration.cord_longitudinal import (
    load_longitudinal_correspondence,
    resolve_cord_z_transinit,
)
from lightsuite.registration.cord_paths import (
    cord_affine_transform_path,
    cord_bspline_forward_transform_path,
    cord_save_path,
    cord_work_dir,
)
from lightsuite.registration.cord_plots import save_cord_annotation_preview
from lightsuite.registration.elastix.runner import run_transformix
from lightsuite.registration.straightening import load_slicetforms, transform_cord_images_slices
from lightsuite.registration.volume import load_registration_volume

console = Console()

SAMPLE_SPACE_SUBDIR = "sample_space"
ANNOTATION_IN_SAMPLE = "annotation_in_sample_20um.tif"
TEMPLATE_IN_SAMPLE = "template_in_sample_20um.tif"
SEGMENTS_IN_SAMPLE = "segments_in_sample_20um.tif"
MANIFEST_NAME = "sample_space_manifest.json"

__all__ = [
    "ANNOTATION_IN_SAMPLE",
    "CordSampleSpaceInspectPaths",
    "CordSampleSpaceInspectVolumes",
    "HEMISPHERE_IN_SAMPLE",
    "MANIFEST_NAME",
    "SAMPLE_SPACE_SUBDIR",
    "SEGMENTS_IN_SAMPLE",
    "TEMPLATE_IN_SAMPLE",
    "discover_cord_sample_space_paths",
    "export_cord_sample_space",
    "load_cord_sample_space_volumes",
    "sample_space_dir",
]
_SAMPLE_CHANNEL_PATTERN = re.compile(r"chan_(\d+)_sample_straight_20um\.tif$", re.IGNORECASE)


@dataclass(frozen=True)
class CordSampleSpaceInspectPaths:
    sample_space_dir: Path
    manifest_path: Path
    annotation_path: Path
    template_path: Path
    hemisphere_path: Path | None = None
    channel_paths: dict[int, Path] = field(default_factory=dict)
    point_npz_paths: dict[str, Path] = field(default_factory=dict)


@dataclass
class CordSampleSpaceInspectVolumes:
    annotation: np.ndarray
    template: np.ndarray
    channels: dict[int, np.ndarray]
    hemisphere: np.ndarray | None = None
    point_layers: dict[str, np.ndarray] = field(default_factory=dict)


def sample_space_dir(save_path: Path) -> Path:
    return save_path / "volume_registered" / SAMPLE_SPACE_SUBDIR


def _label_from_stem(stem: str, suffix: str) -> str:
    if stem.endswith(suffix):
        return stem[: -len(suffix)]
    return stem


def discover_cord_sample_space_paths(config: SpinalCordPipelineConfig) -> CordSampleSpaceInspectPaths:
    """Discover straightened sample-space export TIFFs under volume_registered/sample_space/."""
    save_path = config.sample.save_path.expanduser()
    out_dir = sample_space_dir(save_path)
    if not out_dir.is_dir():
        msg = (
            f"Missing {out_dir}. Run 'lightsuite spinal export --space sample' "
            "(or --space both) first."
        )
        raise FileNotFoundError(msg)

    manifest_path = out_dir / MANIFEST_NAME
    channel_paths: dict[int, Path] = {}
    annotation_path = out_dir / ANNOTATION_IN_SAMPLE
    template_path = out_dir / TEMPLATE_IN_SAMPLE
    hemisphere_path: Path | None = out_dir / HEMISPHERE_IN_SAMPLE

    if manifest_path.is_file():
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        for key, path_str in manifest.get("channel_paths", {}).items():
            channel_paths[int(key)] = Path(path_str).resolve()
        if manifest.get("annotation_path"):
            annotation_path = Path(manifest["annotation_path"]).resolve()
        if manifest.get("template_path"):
            template_path = Path(manifest["template_path"]).resolve()
        if manifest.get("hemisphere_path"):
            hemisphere_path = Path(manifest["hemisphere_path"]).resolve()
        elif not (out_dir / HEMISPHERE_IN_SAMPLE).is_file():
            hemisphere_path = None
    else:
        for path in sorted(out_dir.glob("chan_*_sample_straight_20um.tif")):
            match = _SAMPLE_CHANNEL_PATTERN.fullmatch(path.name)
            if match is None:
                continue
            channel_paths[int(match.group(1))] = path.resolve()

    if not channel_paths:
        msg = (
            f"No sample-space channel TIFFs found in {out_dir}. "
            "Expected files like chan_01_sample_straight_20um.tif."
        )
        raise FileNotFoundError(msg)
    if not annotation_path.is_file():
        msg = f"Missing warped annotation volume: {annotation_path}"
        raise FileNotFoundError(msg)
    if not template_path.is_file():
        msg = f"Missing warped template volume: {template_path}"
        raise FileNotFoundError(msg)

    point_npz_paths: dict[str, Path] = {}
    volume_registered = save_path / "volume_registered"
    for path in sorted(volume_registered.glob("*_sample_coords.npz")):
        label = _label_from_stem(path.stem, "_sample_coords")
        point_npz_paths[label] = path.resolve()

    if hemisphere_path is not None and not hemisphere_path.is_file():
        hemisphere_path = None

    return CordSampleSpaceInspectPaths(
        sample_space_dir=out_dir.resolve(),
        manifest_path=manifest_path.resolve(),
        annotation_path=annotation_path.resolve(),
        template_path=template_path.resolve(),
        hemisphere_path=hemisphere_path.resolve() if hemisphere_path is not None else None,
        channel_paths=channel_paths,
        point_npz_paths=point_npz_paths,
    )


def load_cord_sample_space_volumes(
    config: SpinalCordPipelineConfig,
    *,
    paths: CordSampleSpaceInspectPaths | None = None,
) -> CordSampleSpaceInspectVolumes:
    """Load straightened sample channels and warped atlas labels for Napari."""
    from lightsuite.analysis.counts import SAMPLE_POINTS_KEY, load_atlas_points

    paths = paths or discover_cord_sample_space_paths(config)

    channels: dict[int, np.ndarray] = {}
    expected_shape: tuple[int, ...] | None = None
    for ichan, path in paths.channel_paths.items():
        vol = load_registration_volume(path).astype(np.float32, copy=False)
        if expected_shape is None:
            expected_shape = vol.shape
        elif vol.shape != expected_shape:
            msg = f"{path.name} shape {vol.shape} != expected {expected_shape}"
            raise ValueError(msg)
        channels[ichan] = vol

    annotation = load_registration_volume(paths.annotation_path).astype(np.int32, copy=False)
    template = load_registration_volume(paths.template_path).astype(np.float32, copy=False)
    if expected_shape is not None:
        if annotation.shape != expected_shape:
            msg = (
                f"Annotation shape {annotation.shape} != channel shape {expected_shape}. "
                "Re-run 'lightsuite spinal export --space sample'."
            )
            raise ValueError(msg)
        if template.shape != expected_shape:
            msg = (
                f"Template shape {template.shape} != channel shape {expected_shape}. "
                "Re-run 'lightsuite spinal export --space sample'."
            )
            raise ValueError(msg)

    point_layers: dict[str, np.ndarray] = {}
    for label, npz_path in paths.point_npz_paths.items():
        point_layers[label] = load_atlas_points(npz_path, key=SAMPLE_POINTS_KEY)

    hemisphere: np.ndarray | None = None
    if paths.hemisphere_path is not None and paths.hemisphere_path.is_file():
        hemisphere = load_registration_volume(paths.hemisphere_path).astype(np.uint8, copy=False)
        if expected_shape is not None and hemisphere.shape != expected_shape:
            msg = (
                f"Hemisphere shape {hemisphere.shape} != channel shape {expected_shape}. "
                "Re-run 'lightsuite spinal export --space sample'."
            )
            raise ValueError(msg)

    return CordSampleSpaceInspectVolumes(
        annotation=annotation,
        template=template,
        channels=channels,
        hemisphere=hemisphere,
        point_layers=point_layers,
    )


def _build_native_segment_volume(
    native_shape: tuple[int, int, int],
    segments_df,
) -> np.ndarray:
    """Build (Z, Y, X) volume with segment row index per rostrocaudal slice."""
    nz, ny, nx = native_shape
    vol = np.zeros((nz, ny, nx), dtype=np.float32)
    for idx, row in enumerate(segments_df.itertuples(index=False), start=1):
        start = int(row.Start)
        end = int(row.End)
        vol[start - 1 : end, :, :] = float(idx)
    return vol


def export_cord_sample_space(
    config: SpinalCordPipelineConfig,
    *,
    transform_params: CordTransformParamsCheckpoint,
    checkpoint: CordRegOptsCheckpoint,
    save_volume: bool = True,
) -> Path:
    """Warp atlas labels onto the straightened registration grid."""
    save_path = cord_save_path(config)
    out_dir = sample_space_dir(save_path)
    if save_volume:
        out_dir.mkdir(parents=True, exist_ok=True)

    spacing_mm = config.registration.resolution_um * 1e-3
    elastix_affine_path = cord_affine_transform_path(config)
    bspline_fwd_path = cord_bspline_forward_transform_path(config)
    if not bspline_fwd_path.is_file():
        fwd = transform_params.tform_bspline_atlas20um_to_samp_20um_px
        if fwd and Path(fwd).is_file():
            bspline_fwd_path = Path(fwd)
        else:
            msg = (
                "Missing forward B-spline (bspline_atlas_to_samp_20um.txt). "
                "Re-run 'lightsuite spinal register'."
            )
            raise FileNotFoundError(msg)

    straightvol = tifffile.imread(checkpoint.straightvol_path).astype(np.float32)
    tv = tifffile.imread(checkpoint.tv_path).astype(np.float32)
    av = tifffile.imread(checkpoint.av_path).astype(np.uint16)
    nslices = checkpoint.ikeeprange[1] - checkpoint.ikeeprange[0] + 1
    correspondence = load_longitudinal_correspondence(save_path)
    transinit = resolve_cord_z_transinit(nslices, tv.shape[2], correspondence)
    target_shape = tuple(straightvol.shape)

    console.print("Warping Fiederling atlas to straightened sample grid...")
    t0 = time.perf_counter()
    tvtemp = ndimage.median_filter(tv, size=3)

    avaffine = warp_cord_atlas_to_straightvol(
        av,
        transinit=transinit,
        elastix_affine_path=elastix_affine_path,
        output_shape=target_shape,
        spacing_mm=spacing_mm,
        work_dir=cord_work_dir(config, "transformix", "sample_export", "av_affine"),
        nearest=True,
    )
    annotation_sample = run_transformix(
        moving_volume=avaffine.astype(np.float32),
        transform_path=bspline_fwd_path,
        output_dir=cord_work_dir(config, "transformix", "sample_export", "av_bspline"),
        spacing_mm=spacing_mm,
        nearest=True,
    )

    tvaffine = warp_cord_atlas_to_straightvol(
        tvtemp,
        transinit=transinit,
        elastix_affine_path=elastix_affine_path,
        output_shape=target_shape,
        spacing_mm=spacing_mm,
        work_dir=cord_work_dir(config, "transformix", "sample_export", "tv_affine"),
        nearest=False,
    )
    template_sample = run_transformix(
        moving_volume=tvaffine,
        transform_path=bspline_fwd_path,
        output_dir=cord_work_dir(config, "transformix", "sample_export", "tv_bspline"),
        spacing_mm=spacing_mm,
        nearest=False,
    )

    atlas_volumes = load_fiederling_atlas_volumes(config.atlas)
    paths = resolve_fiederling_paths(config.atlas.atlas_dir)
    seg_native = _build_native_segment_volume(
        atlas_volumes.template.shape,
        atlas_volumes.segments,
    )
    seg_affine = warp_cord_atlas_to_straightvol(
        seg_native,
        transinit=transinit,
        elastix_affine_path=elastix_affine_path,
        output_shape=target_shape,
        spacing_mm=spacing_mm,
        work_dir=cord_work_dir(config, "transformix", "sample_export", "seg_affine"),
        nearest=True,
    )
    segments_sample = run_transformix(
        moving_volume=seg_affine,
        transform_path=bspline_fwd_path,
        output_dir=cord_work_dir(config, "transformix", "sample_export", "seg_bspline"),
        spacing_mm=spacing_mm,
        nearest=True,
    )

    hem_native = load_fiederling_hemisphere_native(config.atlas.atlas_dir).astype(np.float32)
    hem_affine = warp_cord_atlas_to_straightvol(
        hem_native,
        transinit=transinit,
        elastix_affine_path=elastix_affine_path,
        output_shape=target_shape,
        spacing_mm=spacing_mm,
        work_dir=cord_work_dir(config, "transformix", "sample_export", "hem_affine"),
        nearest=True,
    )
    hemisphere_sample = run_transformix(
        moving_volume=hem_affine.astype(np.float32),
        transform_path=bspline_fwd_path,
        output_dir=cord_work_dir(config, "transformix", "sample_export", "hem_bspline"),
        spacing_mm=spacing_mm,
        nearest=True,
    )

    finvol = load_registration_volumes(checkpoint)
    perm = [p - 1 for p in transform_params.how_to_perm]
    finvol = np.transpose(finvol, perm + [3])
    yrange = transform_params.samp_ikeepy
    xrange = transform_params.samp_ikeepx
    zrange = transform_params.samp_ikeeplong
    finvol = finvol[
        yrange[0] - 1 : yrange[1],
        xrange[0] - 1 : xrange[1],
        zrange[0] - 1 : zrange[1],
        :,
    ]
    tforms = load_slicetforms(transform_params.slicetforms_path)
    sizetv = tuple(transform_params.atlassize[:2])
    channel_paths: dict[str, str] = {}

    if save_volume:
        save_registration_volume(
            np.rint(annotation_sample).astype(np.uint16),
            out_dir / ANNOTATION_IN_SAMPLE,
        )
        save_registration_volume(
            warp_output_to_uint16(template_sample),
            out_dir / TEMPLATE_IN_SAMPLE,
        )
        save_registration_volume(
            np.rint(segments_sample).astype(np.uint16),
            out_dir / SEGMENTS_IN_SAMPLE,
        )
        save_registration_volume(
            np.rint(hemisphere_sample).astype(np.uint8),
            out_dir / HEMISPHERE_IN_SAMPLE,
        )

        for ich in range(finvol.shape[3]):
            currvol = transform_cord_images_slices(finvol[:, :, :, ich], tforms, sizetv)
            ch_path = out_dir / f"chan_{ich + 1:02d}_sample_straight_20um.tif"
            save_registration_volume(
                warp_output_to_uint16(currvol),
                ch_path,
            )
            channel_paths[str(ich + 1)] = str(ch_path)

        volmax = float(np.quantile(straightvol, 0.999))
        volplot = np.clip(straightvol / max(volmax, 1.0) * 255.0, 0, 255).astype(np.uint8)
        save_cord_annotation_preview(
            volplot,
            np.rint(annotation_sample).astype(np.uint16),
            out_dir / "registration_export_sample.png",
            title="export sample",
        )

    registres_um = float(transform_params.registrationres_um[0])
    manifest = {
        "space": "sample",
        "grid": "straightened_registration",
        "shape_yxz": list(target_shape),
        "voxel_um": [registres_um, registres_um, registres_um],
        "straightened": True,
        "channel_paths": channel_paths,
        "annotation_path": str(out_dir / ANNOTATION_IN_SAMPLE),
        "template_path": str(out_dir / TEMPLATE_IN_SAMPLE),
        "segments_path": str(out_dir / SEGMENTS_IN_SAMPLE),
        "hemisphere_path": str(out_dir / HEMISPHERE_IN_SAMPLE),
        "segments_csv": str(paths.segments_csv),
    }
    manifest_path = out_dir / MANIFEST_NAME
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    console.print(
        f"Cord sample-space export done in {time.perf_counter() - t0:.1f}s under {out_dir}"
    )
    return manifest_path
