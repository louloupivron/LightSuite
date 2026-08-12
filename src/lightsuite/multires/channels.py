"""Helpers for multichannel mesoSPIM co-registration."""

from __future__ import annotations

from lightsuite.multires.config_models import MultiresPipelineConfig


def multires_channel_names(cfg: MultiresPipelineConfig) -> list[str]:
    """Sorted channel slugs declared in ``multires.channels``."""
    channels = cfg.multires.channels
    if not channels:
        return []
    return sorted(channels)


def resolved_reference_channel(cfg: MultiresPipelineConfig) -> str | None:
    return cfg.multires.registration.reference_channel


def default_apply_transform_to(
    cfg: MultiresPipelineConfig,
    reference_channel: str,
) -> list[str]:
    """All declared channels except the reference (Elastix apply-transform targets)."""
    return [name for name in multires_channel_names(cfg) if name != reference_channel]


def resolved_apply_transform_to(cfg: MultiresPipelineConfig) -> list[str]:
    ref = resolved_reference_channel(cfg)
    if ref is None:
        return []
    explicit = cfg.multires.registration.apply_transform_to
    if explicit is not None:
        return [name for name in explicit if name != ref]
    return default_apply_transform_to(cfg, ref)


__all__ = [
    "default_apply_transform_to",
    "multires_channel_names",
    "resolved_apply_transform_to",
    "resolved_reference_channel",
]
