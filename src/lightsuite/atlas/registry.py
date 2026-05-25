"""Atlas path resolution (Python port of resolveBrainAtlasConfig.m)."""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path

ATLAS_FILES: dict[str, dict[str, str]] = {
    "allen": {
        "template": "average_template_10.nii.gz",
        "annotation": "annotation_10.nii.gz",
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
    structures_csv_path: Path | None
    supports_parcellation: bool


def _search_dirs(explicit: Path | None) -> list[Path]:
    dirs: list[Path] = []
    if explicit is not None:
        dirs.append(explicit.expanduser().resolve())
    env = os.environ.get("LIGHTSUITE_ATLAS_PATH", "")
    for part in env.split(os.pathsep):
        if part.strip():
            dirs.append(Path(part.strip()).expanduser().resolve())
    dirs.append(Path.cwd())
    # Common split Perens layout: template under perens/, annotation in parent
    for base in list(dirs):
        dirs.append(base / "LSFM_atlas_files")
        dirs.append(base / "LSFM_atlas_files" / "perens")
    return dirs


def resolve_brain_atlas(
    brain_atlas: str = "allen",
    atlas_dir: Path | None = None,
) -> AtlasPaths:
    """Resolve template and annotation NIfTI paths for a brain atlas."""
    atlas_id = brain_atlas.lower().strip()
    if atlas_id not in ATLAS_FILES:
        msg = f"Unknown brain atlas '{brain_atlas}'. Expected: {', '.join(ATLAS_FILES)}"
        raise ValueError(msg)

    files = ATLAS_FILES[atlas_id]
    template_name = files["template"]
    annotation_name = files["annotation"]

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
        structures_csv_path=structures_csv,
        supports_parcellation=supports_parcellation,
    )
