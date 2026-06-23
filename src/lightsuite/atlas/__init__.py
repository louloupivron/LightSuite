"""Brain atlas providers."""

from lightsuite.atlas.io import load_atlas_volume
from lightsuite.atlas.registry import (
    AtlasPaths,
    atlas_display_provider,
    atlas_display_provider_from_config,
    resolve_brain_atlas,
    resolve_brain_atlas_from_config,
    resolve_brain_atlas_with_config,
)

__all__ = [
    "AtlasPaths",
    "atlas_display_provider",
    "atlas_display_provider_from_config",
    "load_atlas_volume",
    "resolve_brain_atlas",
    "resolve_brain_atlas_from_config",
    "resolve_brain_atlas_with_config",
]
