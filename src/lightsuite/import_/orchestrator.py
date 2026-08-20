"""Shared driver for annotation import across the brain, cord, and multires pipelines.

Each pipeline supplies an :class:`AnnotationImporter` that owns the coordinate
transform and the output filenames; everything around it — resolving specs,
loading, native-bounds validation, the summary JSON — lives here.
"""

from __future__ import annotations

import json
import shutil
import time
from pathlib import Path
from typing import Protocol, runtime_checkable

from rich.console import Console

from lightsuite.config.models import (
    AnnotationFormat,
    AnnotationImportConfig,
    ImportConfig,
)
from lightsuite.import_.adapters import load_annotation, prepare_points_for_sample
from lightsuite.import_.models import AnnotationImportResult, ImportedMask, ImportedPoints
from lightsuite.import_.sample_reference import SampleReference

console = Console()

SUMMARY_NAME = "import_annotations_summary.json"


def slug_for_label(label: str) -> str:
    """Filesystem-safe stem for annotation outputs."""
    cleaned = "".join(ch if ch.isalnum() or ch in "-_" else "_" for ch in label.strip())
    return cleaned or "annotation"


@runtime_checkable
class AnnotationImporter(Protocol):
    """Pipeline-specific transform + output writer for one annotation source."""

    @property
    def reference(self) -> SampleReference:
        """Native grid used to validate incoming coordinates and mask shapes."""

    def import_points(self, points: ImportedPoints, *, slug: str) -> AnnotationImportResult:
        """Transform in-bounds native points and write this pipeline's outputs."""

    def import_mask(self, mask: ImportedMask, *, slug: str) -> AnnotationImportResult:
        """Transform a native mask volume and write this pipeline's outputs."""


def resolve_annotation_specs(
    import_config: ImportConfig | None,
    annotations: list[AnnotationImportConfig] | None,
    *,
    save_path: Path | None = None,
) -> list[AnnotationImportConfig]:
    """Pick explicit annotations over the config block, with a shared error message.

    When ``import.annotations`` is empty after a vendor ``convert-annotations`` run,
    fall back to ``<save_path>/converted/convert_annotations_summary.json``.
    """
    if annotations is not None:
        return annotations
    if import_config is not None and import_config.annotations:
        return import_config.annotations

    if save_path is not None:
        from lightsuite.import_.convert import CONVERT_SUMMARY_NAME

        summary_path = Path(save_path).expanduser() / "converted" / CONVERT_SUMMARY_NAME
        if summary_path.is_file():
            payload = json.loads(summary_path.read_text(encoding="utf-8"))
            rows = payload.get("annotations") or []
            specs: list[AnnotationImportConfig] = []
            for row in rows:
                if not isinstance(row, dict) or not row.get("path"):
                    continue
                fmt_raw = str(row.get("format") or "points_csv")
                try:
                    fmt = AnnotationFormat(fmt_raw)
                except ValueError:
                    continue
                specs.append(
                    AnnotationImportConfig(
                        format=fmt,
                        path=Path(str(row["path"])),
                        label=str(row.get("label") or ""),
                    )
                )
            if specs:
                return specs

    msg = (
        "No import.annotations configured. Enable Import segmentation in Config, "
        "run Convert annotations for a vendor suite, or add native points_csv / "
        "mask_tiff layers (suite=Native)."
    )
    raise ValueError(msg)


def resolve_write_csv(import_config: ImportConfig | None, write_csv: bool | None) -> bool:
    if write_csv is not None:
        return write_csv
    return import_config.write_csv if import_config is not None else True


def require_transformix() -> None:
    if shutil.which("transformix") is None:
        msg = "transformix must be on PATH for annotation import."
        raise RuntimeError(msg)


def _summary_entry(result: AnnotationImportResult) -> dict[str, object]:
    return {
        "label": result.label,
        "kind": result.kind,
        "n_input": result.n_input,
        "n_atlas": result.n_atlas,
        "n_sample": result.n_sample,
        "atlas_points_path": str(result.atlas_points_path) if result.atlas_points_path else None,
        "atlas_mask_path": str(result.atlas_mask_path) if result.atlas_mask_path else None,
        "atlas_csv_path": str(result.atlas_csv_path) if result.atlas_csv_path else None,
        "sample_points_path": str(result.sample_points_path) if result.sample_points_path else None,
        "sample_mask_path": str(result.sample_mask_path) if result.sample_mask_path else None,
    }


def write_import_summary(output_dir: Path, results: list[AnnotationImportResult]) -> Path:
    summary_path = output_dir / SUMMARY_NAME
    summary_path.write_text(
        json.dumps([_summary_entry(r) for r in results], indent=2),
        encoding="utf-8",
    )
    return summary_path


def run_annotation_import(
    specs: list[AnnotationImportConfig],
    *,
    importer: AnnotationImporter,
    output_dir: Path,
    supports_masks: bool = True,
    mask_unsupported_message: str | None = None,
) -> list[AnnotationImportResult]:
    """Load, validate, and transform every annotation source, then write a summary."""
    reference = importer.reference
    output_dir.mkdir(parents=True, exist_ok=True)

    console.print(
        f"Importing {len(specs)} annotation source(s) "
        f"(native grid {reference.shape_yxz}, voxel_um={reference.voxel_um})..."
    )
    t0 = time.perf_counter()

    results: list[AnnotationImportResult] = []
    for spec in specs:
        if spec.format == AnnotationFormat.MASK_TIFF:
            if not supports_masks:
                raise NotImplementedError(
                    mask_unsupported_message or "Mask import is not supported by this pipeline."
                )
            results.append(_run_mask(spec, importer=importer))
        else:
            results.append(_run_points(spec, importer=importer))

    summary_path = write_import_summary(output_dir, results)
    console.print(
        f"[green]Import complete[/green] in {time.perf_counter() - t0:.1f}s. "
        f"Summary: {summary_path}"
    )
    return results


def _run_points(
    spec: AnnotationImportConfig,
    *,
    importer: AnnotationImporter,
) -> AnnotationImportResult:
    loaded = load_annotation(spec)
    if not isinstance(loaded, ImportedPoints):
        msg = f"Expected point annotation for {spec.path}"
        raise TypeError(msg)

    prepared = prepare_points_for_sample(loaded, reference=importer.reference)
    if prepared.coordinates.size == 0:
        console.print(f"[yellow]No in-bounds points for {prepared.label}[/yellow]")
        return AnnotationImportResult(
            label=prepared.label,
            kind="points",
            n_input=int(loaded.coordinates.shape[0]),
            n_atlas=0,
        )

    result = importer.import_points(prepared, slug=slug_for_label(prepared.label))
    if not result.n_input:
        result.n_input = int(loaded.coordinates.shape[0])
    return result


def _run_mask(
    spec: AnnotationImportConfig,
    *,
    importer: AnnotationImporter,
) -> AnnotationImportResult:
    loaded = load_annotation(spec)
    if not isinstance(loaded, ImportedMask):
        msg = f"Expected mask annotation for {spec.path}"
        raise TypeError(msg)
    return importer.import_mask(loaded, slug=slug_for_label(loaded.label))


__all__ = [
    "AnnotationImporter",
    "require_transformix",
    "resolve_annotation_specs",
    "resolve_write_csv",
    "run_annotation_import",
    "slug_for_label",
    "write_import_summary",
]
