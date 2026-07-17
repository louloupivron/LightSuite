"""Atlas path resolution (Python port of resolveBrainAtlasConfig.m)."""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

ATLAS_FILES: dict[str, dict[str, str]] = {
    "allen": {
        "template": "average_template_10.nii.gz",
        "annotation": "annotation_10.nii.gz",
        "boundary": "annotation_boundary_10.nii.gz",
    },
    "perens": {
        "template": "gubra_template_olf.nii.gz",
        "annotation": "gubra_ano_olf.nii.gz",
    },
}


@dataclass(frozen=True)
class AtlasPaths:
    brain_atlas: str
    atlas_dir: Path
    template_path: Path
    annotation_path: Path
    boundary_path: Path | None
    structures_csv_path: Path | None
    supports_parcellation: bool
    atlas_source: str = "files"
    brainglobe_name: str | None = None


if TYPE_CHECKING:
    from lightsuite.config.models import AtlasConfig


def _search_dirs(explicit: Path | None) -> list[Path]:
    dirs: list[Path] = []
    if explicit is not None:
        dirs.append(explicit.expanduser().resolve())
    env = os.environ.get("LIGHTSUITE_ATLAS_PATH", "")
    for part in env.split(os.pathsep):
        if part.strip():
            dirs.append(Path(part.strip()).expanduser().resolve())
    dirs.append(Path.cwd())
    for base in list(dirs):
        dirs.append(base / "LSFM_atlas_files")
        dirs.append(base / "LSFM_atlas_files" / "perens")
    return dirs


def _resolve_files_atlas(
    atlas_id: str,
    atlas_dir: Path | None,
) -> AtlasPaths:
    files = ATLAS_FILES[atlas_id]
    template_name = files["template"]
    annotation_name = files["annotation"]
    boundary_name = files.get("boundary")

    if atlas_dir is not None:
        base = atlas_dir.expanduser().resolve()
        if not base.is_dir():
            msg = f"atlas_dir is not a directory: {base}"
            raise FileNotFoundError(msg)
        template_path = base / template_name
        annotation_path = base / annotation_name
        if not template_path.is_file() or not annotation_path.is_file():
            msg = (
                f"Atlas files not found under {base}. "
                f"Expected {template_name} and {annotation_name}."
            )
            raise FileNotFoundError(msg)
        resolved_dir = base
    else:
        template_path = annotation_path = None
        resolved_dir = None
        for directory in _search_dirs(None):
            tpl = directory / template_name
            ann = directory / annotation_name
            if tpl.is_file() and ann.is_file():
                template_path = tpl
                annotation_path = ann
                resolved_dir = directory
                break
        if template_path is None or annotation_path is None or resolved_dir is None:
            msg = (
                f'Brain atlas "{atlas_id}" not found. Set atlas.atlas_dir in config or '
                f"add a folder containing {template_name} to LIGHTSUITE_ATLAS_PATH."
            )
            raise FileNotFoundError(msg)

    boundary_path: Path | None = None
    if boundary_name is not None:
        candidate = resolved_dir / boundary_name
        if candidate.is_file():
            boundary_path = candidate

    structures_csv: Path | None = None
    if atlas_id == "perens":
        csv_candidate = resolved_dir / "ARA2_annotation_info_avail_regions.csv"
        if csv_candidate.is_file():
            structures_csv = csv_candidate

    supports_parcellation = atlas_id == "allen" or structures_csv is not None

    return AtlasPaths(
        brain_atlas=atlas_id,
        atlas_dir=resolved_dir,
        template_path=template_path,
        annotation_path=annotation_path,
        boundary_path=boundary_path,
        structures_csv_path=structures_csv,
        supports_parcellation=supports_parcellation,
        atlas_source="files",
        brainglobe_name=None,
    )


def resolve_brain_atlas(
    brain_atlas: str = "allen",
    atlas_dir: Path | None = None,
    *,
    source: str = "files",
    brainglobe_name: str | None = None,
    resolution_um: float | None = None,
) -> AtlasPaths:
    """Resolve template and annotation paths for a brain atlas."""
    atlas_id = brain_atlas.lower().strip()
    if atlas_id not in ATLAS_FILES and source == "files":
        brainglobe_only = ("princeton", "rat")
        hint = ""
        if atlas_id in brainglobe_only:
            hint = f" Atlas '{atlas_id}' is BrainGlobe-only; set atlas.source: brainglobe."
        msg = (
            f"Unknown brain atlas '{brain_atlas}'. Expected: {', '.join(ATLAS_FILES)}."
            f"{hint}"
        )
        raise ValueError(msg)

    source = source.lower().strip()
    if source == "brainglobe":
        from lightsuite.atlas.brainglobe_backend import (
            default_brainglobe_name,
            resolve_brainglobe_paths,
        )

        bg_name = brainglobe_name or default_brainglobe_name(
            atlas_id,
            resolution_um if resolution_um is not None else 20.0,
        )
        template_path, annotation_path, boundary_path, structures_csv, _res, root = (
            resolve_brainglobe_paths(bg_name, provider=atlas_id)
        )
        structures_path = structures_csv if structures_csv.is_file() else None
        supports = structures_path is not None
        return AtlasPaths(
            brain_atlas=atlas_id,
            atlas_dir=root,
            template_path=template_path,
            annotation_path=annotation_path,
            boundary_path=boundary_path,
            structures_csv_path=structures_path,
            supports_parcellation=supports,
            atlas_source="brainglobe",
            brainglobe_name=bg_name,
        )

    if source != "files":
        msg = f"Unknown atlas source {source!r}. Expected 'files' or 'brainglobe'."
        raise ValueError(msg)

    return _resolve_files_atlas(atlas_id, atlas_dir)


def resolve_brain_atlas_from_config(cfg: AtlasConfig) -> AtlasPaths:
    """Resolve atlas paths from pipeline :class:`AtlasConfig`."""
    return resolve_brain_atlas(
        cfg.provider.value,
        cfg.atlas_dir,
        source=cfg.source.value,
        brainglobe_name=cfg.brainglobe_name,
        resolution_um=cfg.resolution_um,
    )


def resolve_brain_atlas_content(
    cfg: AtlasConfig,
    *,
    scratch: Path,
) -> "ResolvedAtlasContent":
    """Resolve atlas paths and optional trimmed working copies (Option A)."""
    from lightsuite.atlas.trim import ResolvedAtlasContent, trim_atlas_to_cache

    base = resolve_brain_atlas_from_config(cfg)
    return trim_atlas_to_cache(base, cfg, scratch=scratch)


def resolve_brain_atlas_with_config(brain_atlas: str, cfg: AtlasConfig) -> AtlasPaths:
    """Resolve atlas using a stored ``brain_atlas`` id and current atlas config."""
    return resolve_brain_atlas(
        brain_atlas,
        cfg.atlas_dir,
        source=cfg.source.value,
        brainglobe_name=cfg.brainglobe_name,
        resolution_um=cfg.resolution_um,
    )


def uses_ccf_id_parcellation(atlas: AtlasPaths) -> bool:
    """True when annotation voxels store Allen CCF structure ids (not ABC indices)."""
    if atlas.atlas_source == "brainglobe":
        return True
    return atlas.brain_atlas in ("perens", "princeton")


def atlas_resolution_um_for_cache(atlas: AtlasPaths, fallback: float) -> float:
    """Best-effort atlas voxel size in µm for cached division-map filenames."""
    if atlas.brainglobe_name:
        for token in ("10um", "15um", "20um", "25um", "39um", "50um"):
            if token in atlas.brainglobe_name:
                return float(token.replace("um", ""))
    if atlas.brain_atlas == "allen":
        return 10.0
    if atlas.brain_atlas == "rat":
        return 39.0
    if atlas.brain_atlas in ("perens", "princeton"):
        return 20.0
    return fallback


def atlas_display_provider(atlas: AtlasPaths) -> str:
    """QC display profile id for plots and GUIs (accounts for BrainGlobe axis order)."""
    from lightsuite.atlas.display import display_provider_for_atlas

    return display_provider_for_atlas(atlas.brain_atlas, atlas_source=atlas.atlas_source)


def atlas_display_provider_from_config(cfg: AtlasConfig) -> str:
    """QC display profile id from pipeline :class:`AtlasConfig`."""
    from lightsuite.atlas.display import display_provider_for_atlas

    return display_provider_for_atlas(cfg.provider.value, atlas_source=cfg.source.value)
