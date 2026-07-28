"""Build multiresolution pair manifests from mesoSPIM exports."""

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path

import SimpleITK as sitk

from lightsuite.mesospim.config_models import MesospimGeometryConfig, MesospimTiffRemapConfig
from lightsuite.mesospim.geometry import (
    apply_image_geometry,
    mesospim_geometry_fields,
    stitched_mosaic_geometry_fields,
)
from lightsuite.mesospim.io import empty_image_from_shape, read_tiff_as_float, tiff_shape
from lightsuite.mesospim.meta import meta_path_for_tiff, parse_mesospim_meta
from lightsuite.multires.models import (
    MANIFEST_FORMAT,
    ManifestVolumeSpec,
    MultiresChannelSpecs,
    MultiresPairManifest,
)
from lightsuite.multires.volume import (
    discover_volume_shape,
    volume_spec_from_geometry,
    volume_spec_from_image,
)


def mesospim_volume_shape(path: Path) -> tuple[int, int, int]:
    """Return ZYX shape for a mesoSPIM hyperstack TIFF or stitched plane-per-file folder."""
    path = path.expanduser().resolve()
    if path.is_dir():
        return discover_volume_shape(path)
    return tiff_shape(path)


def apply_mesospim_geometry_to_shape(
    shape_zyx: tuple[int, int, int],
    meta: dict[str, float | int | str],
    geometry: MesospimGeometryConfig,
) -> sitk.Image:
    """Attach mesoSPIM stage geometry to an arbitrary ZYX shape (no pixel I/O)."""
    image = empty_image_from_shape(shape_zyx)
    return apply_image_geometry(image, meta, geometry)


def apply_stitched_mosaic_geometry(
    stitched_shape_zyx: tuple[int, int, int],
    anchor_meta: dict[str, float | int | str],
    geometry: MesospimGeometryConfig,
) -> sitk.Image:
    """Attach geometry for a TeraStitcher-style mesoSPIM mosaic stitched along Y."""
    spacing, origin, direction = stitched_mosaic_geometry_fields(
        stitched_shape_zyx,
        anchor_meta,
        geometry,
    )
    image = sitk.Image(1, 1, 1, sitk.sitkFloat32)
    image.SetSpacing(spacing)
    image.SetOrigin(origin)
    image.SetDirection(direction)
    return image


def _resolve_overview_meta_path(
    overview_path: Path,
    overview_meta_path: Path | None,
) -> Path:
    if overview_meta_path is not None:
        return overview_meta_path.expanduser().resolve()
    if overview_path.is_file():
        return meta_path_for_tiff(overview_path).expanduser().resolve()
    msg = (
        "overview_meta_path is required when overview_path is a stitched folder "
        f"({overview_path})"
    )
    raise ValueError(msg)


def _overview_volume_spec(
    overview_path: Path,
    overview_meta: dict[str, float | int | str],
    geometry: MesospimGeometryConfig,
) -> ManifestVolumeSpec:
    shape_zyx = mesospim_volume_shape(overview_path)
    if overview_path.is_dir():
        spacing, origin, direction = stitched_mosaic_geometry_fields(
            shape_zyx,
            overview_meta,
            geometry,
        )
        return volume_spec_from_geometry(overview_path, shape_zyx, spacing, origin, direction)
    spacing, origin, direction = mesospim_geometry_fields(shape_zyx, overview_meta, geometry)
    return volume_spec_from_geometry(overview_path, shape_zyx, spacing, origin, direction)


def _roi_volume_spec(
    roi_path: Path,
    roi_meta: dict[str, float | int | str],
    geometry: MesospimGeometryConfig,
) -> ManifestVolumeSpec:
    if not roi_path.is_file():
        msg = f"ROI path must be a mesoSPIM hyperstack TIFF: {roi_path}"
        raise ValueError(msg)
    shape_zyx = tiff_shape(roi_path)
    spacing, origin, direction = mesospim_geometry_fields(shape_zyx, roi_meta, geometry)
    return volume_spec_from_geometry(roi_path, shape_zyx, spacing, origin, direction)


def build_mesospim_pair_manifest(
    *,
    sample_name: str,
    pair_label: str,
    overview_path: Path,
    roi_path: Path,
    geometry: MesospimGeometryConfig | None = None,
    tiff_remap: MesospimTiffRemapConfig | None = None,
    overview_meta_path: Path | None = None,
    roi_meta_path: Path | None = None,
    landmarks_path: Path | None = None,
    output_manifest_path: Path | None = None,
) -> MultiresPairManifest:
    """Convert mesoSPIM TIFF + meta sidecars into a LightSuite multires pair manifest.

    ``overview_path`` may be a single hyperstack TIFF or a stitched plane-per-file
    folder (for example TeraStitcher ``RES(...)`` output). Stitched folders require
    ``overview_meta_path`` pointing at an anchor tile meta file (typically the
    northern tile).
    """
    _ = tiff_remap
    geometry = geometry or MesospimGeometryConfig()
    overview_path = overview_path.expanduser().resolve()
    roi_path = roi_path.expanduser().resolve()

    overview_meta_path = _resolve_overview_meta_path(overview_path, overview_meta_path)
    roi_meta_path = (roi_meta_path or meta_path_for_tiff(roi_path)).expanduser().resolve()
    overview_meta = parse_mesospim_meta(overview_meta_path)
    roi_meta = parse_mesospim_meta(roi_meta_path)

    overview_spec = _overview_volume_spec(overview_path, overview_meta, geometry)
    roi_spec = _roi_volume_spec(roi_path, roi_meta, geometry)

    manifest = MultiresPairManifest(
        format=MANIFEST_FORMAT,
        sample_name=sample_name,
        pair_label=pair_label,
        overview=overview_spec,
        roi=roi_spec,
        provenance={
            "microscope": "mesospim",
            "source_overview": str(overview_path),
            "source_roi": str(roi_path),
            "overview_meta": str(overview_meta_path),
            "roi_meta": str(roi_meta_path),
            "conversion": "lightsuite.multires.vendor.mesospim",
            "conversion_timestamp": datetime.now(timezone.utc).isoformat(),
            "overview_layout": "stitched_folder" if overview_path.is_dir() else "hyperstack_tiff",
        },
        landmarks_path=str(landmarks_path) if landmarks_path is not None else None,
    )

    if output_manifest_path is not None:
        from lightsuite.multires.manifest import save_pair_manifest

        save_pair_manifest(manifest, output_manifest_path)

    return manifest


def build_mesospim_multichannel_pair_manifest(
    *,
    sample_name: str,
    pair_label: str,
    channels: dict[str, dict[str, Path]],
    reference_channel: str,
    geometry: MesospimGeometryConfig | None = None,
    overview_meta_path: Path | None = None,
    roi_meta_by_channel: dict[str, Path] | None = None,
    landmarks_path: Path | None = None,
    output_manifest_path: Path | None = None,
) -> MultiresPairManifest:
    """Build a multichannel manifest sharing geometry across lasers.

    Each channel entry must provide ``overview`` and ``roi`` paths. Overview paths
    may be stitched folders; ROI paths must be mesoSPIM hyperstack TIFFs.
    """
    if reference_channel not in channels:
        msg = f"reference_channel {reference_channel!r} missing from channels"
        raise KeyError(msg)

    geometry = geometry or MesospimGeometryConfig()
    roi_meta_by_channel = roi_meta_by_channel or {}
    channel_specs: dict[str, MultiresChannelSpecs] = {}

    for channel_name, paths in channels.items():
        overview_path = paths["overview"].expanduser().resolve()
        roi_path = paths["roi"].expanduser().resolve()
        overview_meta = _resolve_overview_meta_path(overview_path, overview_meta_path)
        roi_meta_path = (
            roi_meta_by_channel.get(channel_name) or meta_path_for_tiff(roi_path)
        ).expanduser().resolve()

        channel_specs[channel_name] = MultiresChannelSpecs(
            overview=_overview_volume_spec(
                overview_path,
                parse_mesospim_meta(overview_meta),
                geometry,
            ),
            roi=_roi_volume_spec(roi_path, parse_mesospim_meta(roi_meta_path), geometry),
        )

    ref_specs = channel_specs[reference_channel]
    manifest = MultiresPairManifest(
        format=MANIFEST_FORMAT,
        sample_name=sample_name,
        pair_label=pair_label,
        overview=ref_specs.overview,
        roi=ref_specs.roi,
        reference_channel=reference_channel,
        channels=channel_specs,
        provenance={
            "microscope": "mesospim",
            "conversion": "lightsuite.multires.vendor.mesospim",
            "conversion_timestamp": datetime.now(timezone.utc).isoformat(),
            "channels": ",".join(sorted(channels)),
            "reference_channel": reference_channel,
        },
        landmarks_path=str(landmarks_path) if landmarks_path is not None else None,
    )

    if output_manifest_path is not None:
        from lightsuite.multires.manifest import save_pair_manifest

        save_pair_manifest(manifest, output_manifest_path)

    return manifest


def load_mesospim_volumes_from_manifest(
    manifest: MultiresPairManifest,
    *,
    manifest_path: Path,
    geometry: MesospimGeometryConfig | None = None,
    tiff_remap: MesospimTiffRemapConfig | None = None,
):
    """Load mesoSPIM volumes using manifest paths plus axis remapping (notebook QA helper)."""
    geometry = geometry or MesospimGeometryConfig()
    tiff_remap = tiff_remap or MesospimTiffRemapConfig()
    overview_path = (manifest_path.parent / manifest.overview.volume_path).resolve()
    roi_path = (manifest_path.parent / manifest.roi.volume_path).resolve()
    if overview_path.is_file() is False:
        overview_path = Path(manifest.overview.volume_path).expanduser().resolve()
    if roi_path.is_file() is False:
        roi_path = Path(manifest.roi.volume_path).expanduser().resolve()

    overview_meta_path = Path(manifest.provenance.get("overview_meta", ""))
    if overview_meta_path.is_file():
        overview_meta = parse_mesospim_meta(overview_meta_path)
    elif overview_path.is_file():
        overview_meta = parse_mesospim_meta(meta_path_for_tiff(overview_path))
    else:
        msg = "Cannot resolve overview metadata for stitched folder manifest"
        raise ValueError(msg)
    roi_meta = parse_mesospim_meta(meta_path_for_tiff(roi_path))

    if overview_path.is_dir():
        from lightsuite.multires.volume import load_volume_array

        overview_arr = load_volume_array(overview_path)
        overview = sitk.GetImageFromArray(overview_arr.astype("float32", copy=False))
        spacing, origin, direction = stitched_mosaic_geometry_fields(
            mesospim_volume_shape(overview_path),
            overview_meta,
            geometry,
        )
        overview.SetSpacing(spacing)
        overview.SetOrigin(origin)
        overview.SetDirection(direction)
    else:
        overview = read_tiff_as_float(
            overview_path,
            overview_path=overview_path,
            roi_path=roi_path,
            remap=tiff_remap,
        )
        apply_image_geometry(overview, overview_meta, geometry)

    roi = read_tiff_as_float(
        roi_path,
        overview_path=overview_path,
        roi_path=roi_path,
        remap=tiff_remap,
    )
    apply_image_geometry(roi, roi_meta, geometry)
    return overview, roi
