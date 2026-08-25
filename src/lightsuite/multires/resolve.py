"""Resolve or build a multires pair manifest from pipeline config."""

from __future__ import annotations

from pathlib import Path

from lightsuite.mesospim.config_models import MesospimGeometryConfig
from lightsuite.multires.config_models import (
    MesospimGeometryOverride,
    MultiresPipelineConfig,
    MultiresVendorSuite,
)
from lightsuite.multires.manifest import load_pair_manifest, save_pair_manifest
from lightsuite.multires.models import MultiresPairManifest


def _merge_mesospim_geometry(override: MesospimGeometryOverride | None) -> MesospimGeometryConfig | None:
    if override is None:
        return None
    base = MesospimGeometryConfig().model_dump()
    for key, value in override.model_dump(exclude_none=True).items():
        base[key] = value
    return MesospimGeometryConfig(**base)


def _build_mesospim_manifest_from_channels(
    cfg: MultiresPipelineConfig,
    manifest_path: Path,
) -> MultiresPairManifest:
    from lightsuite.multires.vendor.mesospim import build_mesospim_multichannel_pair_manifest

    meso = cfg.multires
    reference_channel = meso.registration.reference_channel
    if reference_channel is None or meso.channels is None:
        msg = "registration.reference_channel and multires.channels are required for mesospim"
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

    overview_meta_by_channel = {
        name: channel.overview_meta_path
        for name, channel in meso.channels.items()
        if channel.overview_meta_path is not None
    }
    roi_meta_by_channel = {
        name: channel.roi_meta_path
        for name, channel in meso.channels.items()
        if channel.roi_meta_path is not None
    }

    mesospim_geometry = meso.mesospim_geometry
    overview_geometry = _merge_mesospim_geometry(
        mesospim_geometry.overview if mesospim_geometry is not None else None
    )
    roi_geometry = _merge_mesospim_geometry(
        mesospim_geometry.roi if mesospim_geometry is not None else None
    )
    return build_mesospim_multichannel_pair_manifest(
        sample_name=cfg.sample.name,
        pair_label=meso.resolved_pair_label(cfg.sample.name),
        channels=channels,
        reference_channel=reference_channel,
        overview_geometry=overview_geometry,
        roi_geometry=roi_geometry,
        overview_meta_path=overview_meta_path,
        overview_meta_by_channel=overview_meta_by_channel or None,
        roi_meta_by_channel=roi_meta_by_channel or None,
        output_manifest_path=manifest_path,
    )


def _build_smartspim_manifest_from_channels(
    cfg: MultiresPipelineConfig,
    manifest_path: Path,
) -> MultiresPairManifest:
    from lightsuite.multires.vendor.smartspim import build_smartspim_multichannel_pair_manifest

    meso = cfg.multires
    if not meso.channels:
        msg = "multires.channels is required for smartspim"
        raise ValueError(msg)

    reference_channel = meso.registration.reference_channel or next(iter(meso.channels))
    if reference_channel not in meso.channels:
        msg = f"registration.reference_channel {reference_channel!r} missing from multires.channels"
        raise KeyError(msg)

    channels = {
        name: {
            "overview": channel.overview,
            "roi": channel.roi,
        }
        for name, channel in meso.channels.items()
    }
    overview_meta_by_channel = {
        name: channel.overview_meta_path
        for name, channel in meso.channels.items()
        if channel.overview_meta_path is not None
    }
    roi_meta_by_channel = {
        name: channel.roi_meta_path
        for name, channel in meso.channels.items()
        if channel.roi_meta_path is not None
    }

    return build_smartspim_multichannel_pair_manifest(
        sample_name=cfg.sample.name,
        pair_label=meso.resolved_pair_label(cfg.sample.name),
        channels=channels,
        reference_channel=reference_channel,
        overview_meta_path=meso.overview_meta_path,
        overview_meta_by_channel=overview_meta_by_channel or None,
        roi_meta_by_channel=roi_meta_by_channel or None,
        output_manifest_path=manifest_path,
    )


def resolve_pair_manifest(
    cfg: MultiresPipelineConfig,
    *,
    rebuild: bool = True,
) -> tuple[MultiresPairManifest, Path]:
    """Return the pair manifest used by multires commands.

    Dispatches on ``multires.vendor.suite``:

    - ``mesospim`` / ``smartspim``: build from ``multires.channels``
    - ``manifest``: load an existing pair manifest JSON
    - ``custom``: run ``multires.vendor.custom_entry`` (``build_pair_manifest``)
    """
    meso = cfg.multires
    manifest_path = meso.resolved_pair_manifest_path(cfg.sample.save_path, cfg.sample.name)
    suite = meso.vendor.suite

    if suite == MultiresVendorSuite.MANIFEST:
        if not manifest_path.is_file():
            msg = f"Pair manifest does not exist: {manifest_path}"
            raise FileNotFoundError(msg)
        return load_pair_manifest(manifest_path), manifest_path

    if suite == MultiresVendorSuite.CUSTOM:
        if not rebuild and manifest_path.is_file():
            return load_pair_manifest(manifest_path), manifest_path
        if meso.vendor.custom_entry is None:
            msg = "multires.vendor.custom_entry is required when vendor.suite is custom"
            raise ValueError(msg)
        from lightsuite.multires.custom import run_custom_manifest_converter

        manifest = run_custom_manifest_converter(
            meso.vendor.custom_entry,
            cfg=cfg,
            output=manifest_path,
        )
        return manifest, manifest_path

    if suite == MultiresVendorSuite.SMARTSPIM:
        if not rebuild and manifest_path.is_file():
            return load_pair_manifest(manifest_path), manifest_path
        manifest = _build_smartspim_manifest_from_channels(cfg, manifest_path)
        return manifest, manifest_path

    # mesospim (default) or legacy configs with channels but no vendor block
    if meso.channels:
        if not rebuild and manifest_path.is_file():
            return load_pair_manifest(manifest_path), manifest_path
        manifest = _build_mesospim_manifest_from_channels(cfg, manifest_path)
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
