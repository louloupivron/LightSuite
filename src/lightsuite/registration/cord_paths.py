"""Path layout for spinal cord pipeline artifacts."""

from __future__ import annotations

from pathlib import Path

from lightsuite.config.models import SpinalCordPipelineConfig

AFFINE_ATLAS_TO_SAMP_FILENAME = "affine_atlas_to_samp_20um.txt"
BSPLINE_SAMP_TO_ATLAS_FILENAME = "bspline_samp_to_atlas_20um.txt"


def cord_save_path(config: SpinalCordPipelineConfig) -> Path:
    return config.sample.save_path.expanduser()


def cord_scratch_root(config: SpinalCordPipelineConfig) -> Path:
    """Per-sample scratch root for ephemeral elastix / transformix workspaces."""
    root = config.sample.scratch.expanduser() / config.sample.name
    root.mkdir(parents=True, exist_ok=True)
    return root


def cord_work_dir(config: SpinalCordPipelineConfig, *parts: str) -> Path:
    """Create (if needed) and return a scratch subdirectory."""
    path = cord_scratch_root(config).joinpath(*parts)
    path.mkdir(parents=True, exist_ok=True)
    return path


def cord_cache_dir(config: SpinalCordPipelineConfig) -> Path:
    path = cord_save_path(config) / "cache"
    path.mkdir(parents=True, exist_ok=True)
    return path


def cord_transforms_dir(config: SpinalCordPipelineConfig) -> Path:
    path = cord_save_path(config) / "transforms"
    path.mkdir(parents=True, exist_ok=True)
    return path


def cord_qc_dir(config: SpinalCordPipelineConfig) -> Path:
    path = cord_save_path(config) / "qc"
    path.mkdir(parents=True, exist_ok=True)
    return path


def resolve_cord_artifact(save_path: Path, subdir: str, filename: str) -> Path:
    """Resolve a persistent artifact, preferring ``subdir/`` over legacy flat layout."""
    save_path = save_path.expanduser()
    nested = save_path / subdir / filename
    if nested.is_file():
        return nested
    legacy = save_path / filename
    if legacy.is_file():
        return legacy
    return nested


def cord_affine_transform_path(config: SpinalCordPipelineConfig) -> Path:
    return resolve_cord_artifact(
        cord_save_path(config),
        "transforms",
        AFFINE_ATLAS_TO_SAMP_FILENAME,
    )


def cord_bspline_transform_path(config: SpinalCordPipelineConfig) -> Path:
    return resolve_cord_artifact(
        cord_save_path(config),
        "transforms",
        BSPLINE_SAMP_TO_ATLAS_FILENAME,
    )


def cord_affine_transform_write_path(config: SpinalCordPipelineConfig) -> Path:
    return cord_transforms_dir(config) / AFFINE_ATLAS_TO_SAMP_FILENAME


def cord_bspline_transform_write_path(config: SpinalCordPipelineConfig) -> Path:
    return cord_transforms_dir(config) / BSPLINE_SAMP_TO_ATLAS_FILENAME
