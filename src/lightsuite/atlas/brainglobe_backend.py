"""BrainGlobe Atlas API backend (optional dependency: brainglobe-atlasapi)."""

from __future__ import annotations

import re
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

import numpy as np

# BrainGlobe hemisphere labels: 1 = left, 2 = right (see brainglobe-atlasapi packagers).
_BRAINGLOBE_LEFT = 1
_BRAINGLOBE_RIGHT = 2

# Default API atlas names when ``atlas.brainglobe_name`` is omitted.
_DEFAULT_NAMES: dict[str, dict[int, str]] = {
    "allen": {10: "allen_mouse_10um"},
    "perens": {20: "perens_lsfm_mouse_20um"},
    "princeton": {20: "princeton_mouse_20um"},
    "rat": {39: "whs_sd_rat_39um"},
}

# BrainGlobe atlases exposed in the GUI (base registry names, without version suffix).
_LIGHTSUITE_BRAINGLOBE_ATLASES: tuple[str, ...] = (
    "allen_mouse_10um",
    "princeton_mouse_20um",
    "whs_sd_rat_39um",
    "perens_lsfm_mouse_20um",
)

_ATLAS_ENTRY_LABELS: dict[str, str] = {
    "allen": "Allen CCF (mouse)",
    "perens": "Perens / Gubra LSFM (mouse)",
    "princeton": "Princeton (mouse)",
    "rat": "Waxholm SD rat",
}


@dataclass(frozen=True)
class BrainGlobeAtlasEntry:
    """One BrainGlobe registry atlas supported by LightSuite."""

    name: str
    provider: str
    resolution_um: float
    label: str
    downloaded: bool = False


def infer_brainglobe_provider(atlas_name: str) -> str | None:
    """Map a BrainGlobe registry name to a LightSuite ``atlas.provider`` id."""
    name = atlas_name.lower().strip()
    if name.startswith("allen_mouse"):
        return "allen"
    if name.startswith("perens_"):
        return "perens"
    if name.startswith("princeton_mouse"):
        return "princeton"
    if name.startswith("whs_") or "_rat_" in name:
        return "rat"
    return None


def infer_brainglobe_resolution_um(atlas_name: str) -> float | None:
    """Parse isotropic resolution in µm from a BrainGlobe registry name."""
    match = re.search(r"(\d+)um", atlas_name.lower())
    if match is None:
        return None
    return float(match.group(1))


def brainglobe_atlas_entry_label(provider: str, resolution_um: float, *, downloaded: bool) -> str:
    """Human-readable dropdown label for one BrainGlobe atlas."""
    provider_label = _ATLAS_ENTRY_LABELS.get(provider, provider)
    label = f"{provider_label} — {resolution_um:g} µm"
    if downloaded:
        label += " (downloaded)"
    return label


def is_supported_lightsuite_brainglobe_atlas(atlas_name: str) -> bool:
    """True when ``atlas_name`` is one of the BrainGlobe atlases LightSuite exposes."""
    return brainglobe_name_key(atlas_name) in _LIGHTSUITE_BRAINGLOBE_ATLASES


def _brainglobe_catalog_sort_key(entry: BrainGlobeAtlasEntry) -> int:
    try:
        return _LIGHTSUITE_BRAINGLOBE_ATLASES.index(brainglobe_name_key(entry.name))
    except ValueError:
        return len(_LIGHTSUITE_BRAINGLOBE_ATLASES)


def _filter_lightsuite_brainglobe_catalog(
    entries: list[BrainGlobeAtlasEntry],
) -> list[BrainGlobeAtlasEntry]:
    filtered = [entry for entry in entries if is_supported_lightsuite_brainglobe_atlas(entry.name)]
    return sorted(filtered, key=_brainglobe_catalog_sort_key)


def _fallback_brainglobe_catalog() -> list[BrainGlobeAtlasEntry]:
    entries: list[BrainGlobeAtlasEntry] = []
    for base_name in _LIGHTSUITE_BRAINGLOBE_ATLASES:
        provider = infer_brainglobe_provider(base_name)
        resolution_um = infer_brainglobe_resolution_um(base_name)
        if provider is None or resolution_um is None:
            continue
        entries.append(
            BrainGlobeAtlasEntry(
                name=base_name,
                provider=provider,
                resolution_um=resolution_um,
                label=brainglobe_atlas_entry_label(provider, resolution_um, downloaded=False),
                downloaded=False,
            )
        )
    return entries


def list_lightsuite_brainglobe_atlases() -> list[BrainGlobeAtlasEntry]:
    """List BrainGlobe atlases LightSuite can register to.

    Uses the BrainGlobe ``last_versions.conf`` registry when
    ``brainglobe-atlasapi`` is installed; otherwise falls back to the built-in
    default name table.
    """
    try:
        require_brainglobe()
        from brainglobe_atlasapi.list_atlases import (
            get_all_atlases_lastversions,
            get_downloaded_atlases,
        )

        available = get_all_atlases_lastversions()
        downloaded = set(get_downloaded_atlases())
    except (ImportError, OSError, ValueError, KeyError):
        return _fallback_brainglobe_catalog()

    entries: list[BrainGlobeAtlasEntry] = []
    for name in sorted(available):
        if not is_supported_lightsuite_brainglobe_atlas(name):
            continue
        provider = infer_brainglobe_provider(name)
        resolution_um = infer_brainglobe_resolution_um(name)
        if provider is None or resolution_um is None:
            continue
        entries.append(
            BrainGlobeAtlasEntry(
                name=name,
                provider=provider,
                resolution_um=resolution_um,
                label=brainglobe_atlas_entry_label(
                    provider,
                    resolution_um,
                    downloaded=name in downloaded,
                ),
                downloaded=name in downloaded,
            )
        )
    if not entries:
        return _fallback_brainglobe_catalog()
    return _filter_lightsuite_brainglobe_catalog(entries)


def brainglobe_name_for_yaml(atlas_name: str) -> str:
    """Strip BrainGlobe version suffix for portable YAML (e.g. ``_v1.2``)."""
    return re.sub(r"_v\d.*$", "", atlas_name.strip())


def brainglobe_name_key(atlas_name: str) -> str:
    """Normalize a registry name for catalog lookup."""
    return brainglobe_name_for_yaml(atlas_name).lower()


def find_brainglobe_catalog_entry(
    entries: list[BrainGlobeAtlasEntry],
    *,
    brainglobe_name: str = "",
    provider: str = "",
    resolution_um: float | None = None,
) -> BrainGlobeAtlasEntry | None:
    """Match a catalog entry from a stored name or provider + resolution."""
    if brainglobe_name.strip():
        key = brainglobe_name_key(brainglobe_name)
        for entry in entries:
            if brainglobe_name_key(entry.name) == key:
                return entry
    provider_key = provider.lower().strip()
    if provider_key and resolution_um is not None:
        for entry in entries:
            if entry.provider == provider_key and abs(entry.resolution_um - resolution_um) < 0.5:
                return entry
    if provider_key:
        for entry in entries:
            if entry.provider == provider_key:
                return entry
    return entries[0] if entries else None


def resolve_brainglobe_name_for_form(
    *,
    brainglobe_name: str | None,
    provider: str,
    resolution_um: float,
) -> str:
    """Resolve the registry name stored in YAML / shown in the GUI."""
    if brainglobe_name and brainglobe_name.strip():
        return brainglobe_name.strip()
    return default_brainglobe_name(provider, resolution_um)


def require_brainglobe() -> None:
    """Raise ImportError with install hint when brainglobe-atlasapi is missing."""
    try:
        import brainglobe_atlasapi  # noqa: F401
    except ImportError as exc:
        msg = (
            "BrainGlobe atlas backend requires brainglobe-atlasapi. "
            "Install with: uv sync --extra atlas"
        )
        raise ImportError(msg) from exc


def default_brainglobe_name(provider: str, resolution_um: float) -> str:
    """Infer a BrainGlobe atlas registry name from provider and resolution."""
    provider = provider.lower().strip()
    table = _DEFAULT_NAMES.get(provider)
    if table is None:
        msg = f"No default BrainGlobe atlas for provider {provider!r}."
        raise ValueError(msg)
    res_key = int(round(resolution_um))
    if res_key in table:
        return table[res_key]
    if provider == "allen":
        return "allen_mouse_10um"
    if provider == "perens":
        return "perens_lsfm_mouse_20um"
    if provider == "princeton":
        return "princeton_mouse_20um"
    if provider == "rat":
        return "whs_sd_rat_39um"
    msg = f"Cannot infer BrainGlobe atlas for {provider} @ {resolution_um} µm."
    raise ValueError(msg)


@lru_cache(maxsize=8)
def _load_bg_atlas(atlas_name: str):
    require_brainglobe()
    from brainglobe_atlasapi.bg_atlas import BrainGlobeAtlas

    return BrainGlobeAtlas(atlas_name)


def resolve_brainglobe_paths(
    atlas_name: str,
    *,
    provider: str,
) -> tuple[Path, Path, Path | None, Path, float, Path | None]:
    """Return paths and metadata for a BrainGlobe atlas.

    Returns
    -------
    template_path, annotation_path, boundary_path, structures_csv, resolution_um, root_dir
    """
    bg = _load_bg_atlas(atlas_name)
    root = Path(bg.root_dir)
    template_path = root / "reference.tiff"
    annotation_path = root / "annotation.tiff"
    structures_csv = root / "structures.csv"
    if not template_path.is_file() or not annotation_path.is_file():
        msg = (
            f"BrainGlobe atlas {atlas_name!r} missing reference.tiff or annotation.tiff under {root}"
        )
        raise FileNotFoundError(msg)
    resolution = float(bg.resolution[0]) if getattr(bg, "resolution", None) is not None else 25.0
    boundary_path: Path | None = None
    return template_path, annotation_path, boundary_path, structures_csv, resolution, root


def load_brainglobe_hemispheres(atlas_name: str) -> np.ndarray:
    """Load the BrainGlobe hemisphere label volume (1=left, 2=right)."""
    bg = _load_bg_atlas(atlas_name)
    hem = getattr(bg, "hemispheres", None)
    if hem is None:
        msg = f"BrainGlobe atlas {atlas_name!r} has no hemisphere volume."
        raise FileNotFoundError(msg)
    return np.asanyarray(hem)


def hemisphere_side_masks_from_brainglobe(
    hemisphere_volume: np.ndarray,
    annotation: np.ndarray,
) -> list[np.ndarray]:
    """Build right/left masks from a BrainGlobe hemisphere volume."""
    hem = np.asanyarray(hemisphere_volume)
    ann = np.asanyarray(annotation)
    if hem.shape != ann.shape:
        msg = f"Hemisphere shape {hem.shape} != annotation {ann.shape}"
        raise ValueError(msg)
    brain = ann > 0
    side0 = (hem == _BRAINGLOBE_RIGHT) & brain
    side1 = (hem == _BRAINGLOBE_LEFT) & brain
    return [side0, side1]
