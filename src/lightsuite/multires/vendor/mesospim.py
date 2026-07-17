"""Build multiresolution pair manifests from mesoSPIM exports."""

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path

from lightsuite.mesospim.config_models import MesospimGeometryConfig, MesospimTiffRemapConfig
from lightsuite.mesospim.geometry import apply_image_geometry
from lightsuite.mesospim.io import empty_image_from_shape, read_tiff_as_float, tiff_shape
from lightsuite.mesospim.meta import meta_path_for_tiff, parse_mesospim_meta
from lightsuite.multires.models import MANIFEST_FORMAT, ManifestVolumeSpec, MultiresPairManifest
from lightsuite.multires.volume import volume_spec_from_image


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
    """Convert mesoSPIM TIFF + meta sidecars into a LightSuite multires pair manifest."""
    geometry = geometry or MesospimGeometryConfig()
    tiff_remap = tiff_remap or MesospimTiffRemapConfig()
    overview_path = overview_path.expanduser().resolve()
    roi_path = roi_path.expanduser().resolve()

    overview_meta_path = (overview_meta_path or meta_path_for_tiff(overview_path)).expanduser()
    roi_meta_path = (roi_meta_path or meta_path_for_tiff(roi_path)).expanduser()
    overview_meta = parse_mesospim_meta(overview_meta_path)
    roi_meta = parse_mesospim_meta(roi_meta_path)

    overview_img = empty_image_from_shape(tiff_shape(overview_path))
    roi_img = empty_image_from_shape(tiff_shape(roi_path))
    apply_image_geometry(overview_img, overview_meta, geometry)
    apply_image_geometry(roi_img, roi_meta, geometry)

    overview_spec = volume_spec_from_image("overview", overview_path, overview_img)
    roi_spec = volume_spec_from_image("roi", roi_path, roi_img)

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

    overview = read_tiff_as_float(
        overview_path,
        overview_path=overview_path,
        roi_path=roi_path,
        remap=tiff_remap,
    )
    roi = read_tiff_as_float(
        roi_path,
        overview_path=overview_path,
        roi_path=roi_path,
        remap=tiff_remap,
    )
    overview_meta = parse_mesospim_meta(meta_path_for_tiff(overview_path))
    roi_meta = parse_mesospim_meta(meta_path_for_tiff(roi_path))
    apply_image_geometry(overview, overview_meta, geometry)
    apply_image_geometry(roi, roi_meta, geometry)
    return overview, roi
