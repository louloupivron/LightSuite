"""Example custom pair-manifest converter for ``multires.vendor.suite: custom``.

Contract — the module must define:

    build_pair_manifest(cfg, output) -> MultiresPairManifest | dict

``cfg`` is a :class:`~lightsuite.multires.config_models.MultiresPipelineConfig`.
Write the manifest JSON to ``output`` (or return a manifest object / dict and
LightSuite will save it). The returned manifest must use format
``lightsuite_multires_pair_v1``.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from lightsuite.multires.models import MANIFEST_FORMAT, ManifestVolumeSpec, MultiresPairManifest


def build_pair_manifest(cfg: Any, output: Path) -> MultiresPairManifest:
    """Minimal stub: copy reference-channel paths from ``multires.channels``."""
    multires = cfg.multires
    channels = multires.channels or {}
    ref = multires.registration.reference_channel or next(iter(channels))
    ref_paths = channels[ref]
    placeholder = ManifestVolumeSpec(
        volume_path=str(ref_paths.overview),
        shape_zyx=[1, 1, 1],
        spacing_um=[1.0, 1.0, 1.0],
        origin_um=[0.0, 0.0, 0.0],
    )
    manifest = MultiresPairManifest(
        format=MANIFEST_FORMAT,
        sample_name=cfg.sample.name,
        pair_label=multires.resolved_pair_label(cfg.sample.name),
        overview=placeholder,
        roi=ManifestVolumeSpec(
            volume_path=str(ref_paths.roi),
            shape_zyx=[1, 1, 1],
            spacing_um=[1.0, 1.0, 1.0],
            origin_um=[0.0, 0.0, 0.0],
        ),
        reference_channel=ref,
        provenance={
            "microscope": "custom",
            "conversion": str(Path(__file__).resolve()),
        },
    )
    _ = output
    return manifest
