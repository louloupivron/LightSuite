"""Path layout for brain pipeline artifacts (qc / _work / durable root)."""

from __future__ import annotations

import shutil
from pathlib import Path

WORK_DIRNAME = "_work"
QC_DIRNAME = "qc"
STATS_DIRNAME = "stats"
IMPORTS_DIRNAME = "imports"
VOLUME_REGISTERED_DIRNAME = "volume_registered"
PREVIEWS_DIRNAME = "previews"  # legacy nested layout (pre-2026)
SAMPLE_SPACE_SUBDIR = "sample_space"

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
IMPORT_SUMMARY_FILENAME = "import_annotations_summary.json"

_QC_AUDIT_JSON_FILES = (
    REGISTRATION_DIAGNOSTICS_FILENAME,
    INIT_DIAGNOSTICS_FILENAME,
    AFFINE_FIT_STATS_FILENAME,
    CORRESPONDENCE_AFFINE_STATS_FILENAME,
    CORRESPONDENCE_LANDMARK_STATS_FILENAME,
)

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
    """Return ``save_path/qc`` for registration preview PNGs."""
    return brain_qc_dir(save_path)


def resolve_brain_qc_preview(save_path: Path, filename: str) -> Path:
    """Read-path for a preview PNG: ``qc/`` first, then legacy ``qc/previews/``."""
    save_path = Path(save_path).expanduser()
    direct = save_path / QC_DIRNAME / filename
    if direct.is_file():
        return direct
    legacy = save_path / QC_DIRNAME / PREVIEWS_DIRNAME / filename
    if legacy.is_file():
        return legacy
    return direct


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


def brain_stats_dir(save_path: Path) -> Path:
    """Return ``save_path/stats`` (created)."""
    path = Path(save_path).expanduser() / STATS_DIRNAME
    path.mkdir(parents=True, exist_ok=True)
    return path


def brain_imports_dir(save_path: Path) -> Path:
    """Return ``save_path/imports`` (created)."""
    path = Path(save_path).expanduser() / IMPORTS_DIRNAME
    path.mkdir(parents=True, exist_ok=True)
    return path


def brain_volume_registered_dir(save_path: Path) -> Path:
    """Return ``save_path/volume_registered`` (created)."""
    path = Path(save_path).expanduser() / VOLUME_REGISTERED_DIRNAME
    path.mkdir(parents=True, exist_ok=True)
    return path


def resolve_brain_stats_file(save_path: Path, filename: str) -> Path:
    """Read-path for a stats CSV/JSON: ``stats/`` then legacy ``volume_registered/``."""
    save_path = Path(save_path).expanduser()
    nested = save_path / STATS_DIRNAME / filename
    if nested.is_file():
        return nested
    legacy_vr = save_path / VOLUME_REGISTERED_DIRNAME / filename
    if legacy_vr.is_file():
        return legacy_vr
    legacy_ss = save_path / VOLUME_REGISTERED_DIRNAME / SAMPLE_SPACE_SUBDIR / filename
    if legacy_ss.is_file():
        return legacy_ss
    return nested


def resolve_brain_imports_file(save_path: Path, filename: str) -> Path:
    """Read-path for an import artifact: ``imports/`` then legacy ``volume_registered/``."""
    save_path = Path(save_path).expanduser()
    nested = save_path / IMPORTS_DIRNAME / filename
    if nested.is_file():
        return nested
    legacy = save_path / VOLUME_REGISTERED_DIRNAME / filename
    if legacy.is_file():
        return legacy
    return nested


def iter_brain_import_paths(save_path: Path, pattern: str) -> list[Path]:
    """Glob import artifacts in ``imports/`` and legacy ``volume_registered/``."""
    root = Path(save_path).expanduser()
    seen: set[Path] = set()
    ordered: list[Path] = []
    for base in (root / IMPORTS_DIRNAME, root / VOLUME_REGISTERED_DIRNAME):
        if not base.is_dir():
            continue
        for path in sorted(base.glob(pattern)):
            resolved = path.resolve()
            if resolved not in seen:
                seen.add(resolved)
                ordered.append(resolved)
    return ordered


def iter_brain_stats_paths(save_path: Path, pattern: str) -> list[Path]:
    """Glob stats tables in ``stats/`` and legacy export locations."""
    root = Path(save_path).expanduser()
    seen: set[Path] = set()
    ordered: list[Path] = []
    bases = [
        root / STATS_DIRNAME,
        root / VOLUME_REGISTERED_DIRNAME,
        root / VOLUME_REGISTERED_DIRNAME / SAMPLE_SPACE_SUBDIR,
    ]
    for base in bases:
        if not base.is_dir():
            continue
        for path in sorted(base.glob(pattern)):
            resolved = path.resolve()
            if resolved not in seen:
                seen.add(resolved)
                ordered.append(resolved)
    return ordered


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


def cleanup_brain_qc_audit_json(save_path: Path) -> None:
    """Remove QC diagnostic JSON checkpoints (preview PNGs under ``qc/`` are kept)."""
    root = Path(save_path).expanduser()
    for name in _QC_AUDIT_JSON_FILES:
        remove_path(root / QC_DIRNAME / name)
        remove_path(root / name)
