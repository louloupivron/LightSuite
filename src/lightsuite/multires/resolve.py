"""Resolve or build a multires pair manifest from pipeline config."""

from __future__ import annotations

from pathlib import Path

from lightsuite.multires.config_models import MultiresPipelineConfig
from lightsuite.multires.manifest import load_pair_manifest, save_pair_manifest
from lightsuite.multires.models import MultiresPairManifest


def resolve_pair_manifest(
    cfg: MultiresPipelineConfig,
    *,
    rebuild: bool = True,
) -> tuple[MultiresPairManifest, Path]:
    """Return the pair manifest used by multires commands.

    If ``multires.channels`` is set, builds (or refreshes) a mesoSPIM pair
    manifest from those paths and writes it to ``pair_manifest`` (or a default
    under ``save_path/converted/``). Otherwise loads an existing pair manifest.
    """
    meso = cfg.multires
    manifest_path = meso.resolved_pair_manifest_path(cfg.sample.save_path, cfg.sample.name)

    if meso.channels:
        if not rebuild and manifest_path.is_file():
            return load_pair_manifest(manifest_path), manifest_path

        from lightsuite.multires.vendor.mesospim import build_mesospim_multichannel_pair_manifest

        reference_channel = meso.registration.reference_channel
        if reference_channel is None:
            msg = "registration.reference_channel is required when multires.channels is set"
            raise ValueError(msg)

        channels = {
            name: {
                "overview": channel.overview,
                "roi": channel.roi,
            }
            for name, channel in meso.channels.items()
        }
        overview_meta_path = meso.overview_meta_path
        if overview_meta_path is None:
            ref_paths = meso.channels[reference_channel]
            overview_meta_path = ref_paths.overview_meta_path

        roi_meta_by_channel = {
            name: channel.roi_meta_path
            for name, channel in meso.channels.items()
            if channel.roi_meta_path is not None
        }

        manifest = build_mesospim_multichannel_pair_manifest(
            sample_name=cfg.sample.name,
            pair_label=meso.resolved_pair_label(cfg.sample.name),
            channels=channels,
            reference_channel=reference_channel,
            overview_meta_path=overview_meta_path,
            roi_meta_by_channel=roi_meta_by_channel or None,
            output_manifest_path=manifest_path,
        )
        return manifest, manifest_path

    if not manifest_path.is_file():
        msg = f"Pair manifest does not exist: {manifest_path}"
        raise FileNotFoundError(msg)
    return load_pair_manifest(manifest_path), manifest_path


def ensure_pair_manifest(cfg: MultiresPipelineConfig) -> Path:
    """Build/refresh the pair manifest if needed and return its path."""
    _manifest, path = resolve_pair_manifest(cfg, rebuild=True)
    return path


__all__ = [
    "ensure_pair_manifest",
    "resolve_pair_manifest",
    "save_pair_manifest",
]
