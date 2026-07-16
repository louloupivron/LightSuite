"""BrainGlobe Atlas API backend (optional dependency: brainglobe-atlasapi)."""

from __future__ import annotations

from functools import lru_cache
from pathlib import Path

import numpy as np

# BrainGlobe hemisphere labels: 1 = left, 2 = right (see brainglobe-atlasapi packagers).
_BRAINGLOBE_LEFT = 1
_BRAINGLOBE_RIGHT = 2

# Default API atlas names when ``atlas.brainglobe_name`` is omitted.
_DEFAULT_NAMES: dict[str, dict[int, str]] = {
    "allen": {10: "allen_mouse_10um", 25: "allen_mouse_25um", 50: "allen_mouse_50um", 100: "allen_mouse_100um"},
    "perens": {20: "perens_lsfm_mouse_20um", 25: "perens_multimodal_lsfm_25um"},
    "princeton": {20: "princeton_mouse_20um"},
    "rat": {39: "whs_sd_rat_39um"},
}


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
        return "allen_mouse_25um"
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
