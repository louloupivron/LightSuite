"""Path layout for brain pipeline artifacts (qc / _work / durable root)."""

from __future__ import annotations

import shutil
from pathlib import Path

WORK_DIRNAME = "_work"
QC_DIRNAME = "qc"
PREVIEWS_DIRNAME = "previews"

# Ephemeral workspaces under ``save_path/_work/``.
ELASTIX_TEMP = "elastix"
ELASTIX_INVERSE_TEMP = "elastix_inverse"
TRANSFORMIX_ANNOTATION_TEMP = "transformix_annotation"
TRANSFORMIX_EXPORT_TEMP = "transformix_export"
TRANSFORMIX_SAMPLE_EXPORT_TEMP = "transformix_sample_export"
TRANSFORMIX_VIEW_TEMP = "transformix_view_registration"
IMPORT_ANNOTATIONS_TEMP = "import_annotations"

# Durable QC filenames under ``save_path/qc/``.
REGISTRATION_DIAGNOSTICS_FILENAME = "registration_diagnostics.json"
INIT_DIAGNOSTICS_FILENAME = "init_registration_diagnostics.json"
AFFINE_FIT_STATS_FILENAME = "affine_fit_stats.json"
CORRESPONDENCE_AFFINE_STATS_FILENAME = "correspondence_affine_stats.json"
CORRESPONDENCE_LANDMARK_STATS_FILENAME = "correspondence_landmark_stats.json"

# Legacy flat-root temp directory names (pre-_work layout).
_LEGACY_WORK_DIRS = (
    "elastix_temp",
    "elastix_inverse_temp",
    "transformix_annotation_temp",
    "transformix_export_temp",
    "transformix_sample_export_temp",
    "transformix_view_registration_temp",
    "import_annotations_temp",
)


def brain_work_root(save_path: Path) -> Path:
    """Return ``save_path/_work`` (created)."""
    root = Path(save_path).expanduser() / WORK_DIRNAME
    root.mkdir(parents=True, exist_ok=True)
    return root


def brain_work_dir(save_path: Path, *parts: str) -> Path:
    """Create and return a subdirectory under ``save_path/_work``."""
    path = brain_work_root(save_path).joinpath(*parts)
    path.mkdir(parents=True, exist_ok=True)
    return path


def brain_qc_dir(save_path: Path) -> Path:
    path = Path(save_path).expanduser() / QC_DIRNAME
    path.mkdir(parents=True, exist_ok=True)
    return path


def brain_qc_previews_dir(save_path: Path) -> Path:
    path = brain_qc_dir(save_path) / PREVIEWS_DIRNAME
    path.mkdir(parents=True, exist_ok=True)
    return path


def resolve_brain_artifact(save_path: Path, *relative: str) -> Path:
    """Resolve a durable artifact, preferring nested layout over legacy flat root.

    ``relative`` is the preferred path under ``save_path`` (e.g. ``("qc", "foo.json")``).
    Falls back to ``save_path / relative[-1]`` when the nested file is missing.
    """
    save_path = Path(save_path).expanduser()
    nested = save_path.joinpath(*relative)
    if nested.is_file():
        return nested
    legacy = save_path / relative[-1]
    if legacy.is_file():
        return legacy
    return nested


def brain_qc_file(save_path: Path, filename: str) -> Path:
    """Write-path for a QC JSON under ``qc/`` (directory created)."""
    return brain_qc_dir(save_path) / filename


def resolve_brain_qc_file(save_path: Path, filename: str) -> Path:
    """Read-path for a QC JSON: ``qc/`` first, then legacy flat root."""
    return resolve_brain_artifact(save_path, QC_DIRNAME, filename)


def remove_path(path: Path) -> None:
    """Remove a file or directory tree; ignore missing paths."""
    path = Path(path).expanduser()
    if path.is_dir():
        shutil.rmtree(path, ignore_errors=True)
    elif path.is_file() or path.is_symlink():
        path.unlink(missing_ok=True)


def cleanup_brain_work(save_path: Path, *parts: str) -> None:
    """Delete ``save_path/_work/<parts...>`` after a successful stage."""
    root = Path(save_path).expanduser()
    if not parts:
        remove_path(root / WORK_DIRNAME)
        return
    target = root / WORK_DIRNAME / Path(*parts)
    remove_path(target)
    work_root = root / WORK_DIRNAME
    if work_root.is_dir() and not any(work_root.iterdir()):
        remove_path(work_root)


def cleanup_legacy_brain_work(save_path: Path) -> None:
    """Remove pre-layout flat temp directories left at the save_path root."""
    root = Path(save_path).expanduser()
    for name in _LEGACY_WORK_DIRS:
        remove_path(root / name)
    work_root = root / WORK_DIRNAME
    if work_root.is_dir() and not any(work_root.iterdir()):
        remove_path(work_root)
