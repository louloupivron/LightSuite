"""Validate LightSuite Sample Space v1 annotation files against sample_reference."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import tifffile

from lightsuite.config.models import AnnotationFormat, AnnotationImportConfig
from lightsuite.import_.adapters import load_points_csv
from lightsuite.import_.normalize import filter_in_bounds_points
from lightsuite.import_.sample_reference import SampleReference


@dataclass
class AnnotationValidationResult:
    """Outcome of validating one native annotation file."""

    label: str
    path: Path
    format: AnnotationFormat
    ok: bool
    n_points: int = 0
    n_in_bounds: int = 0
    n_dropped: int = 0
    shape_yxz: tuple[int, int, int] | None = None
    messages: list[str] = field(default_factory=list)


def validate_points_csv(
    path: Path,
    reference: SampleReference,
    *,
    label: str = "",
    min_in_bounds_fraction: float = 0.95,
) -> AnnotationValidationResult:
    """Check a points CSV against the native sample grid."""
    from lightsuite.config.models import AnnotationImportConfig as Spec

    path = path.expanduser()
    stem = label or path.stem
    result = AnnotationValidationResult(
        label=stem,
        path=path,
        format=AnnotationFormat.POINTS_CSV,
        ok=False,
    )
    if not path.is_file():
        result.messages.append(f"File not found: {path}")
        return result

    loaded = load_points_csv(
        Spec(format=AnnotationFormat.POINTS_CSV, path=path, label=stem)
    )
    ny, nx, nz = reference.shape_tuple
    _xyz, in_bounds = filter_in_bounds_points(
        loaded.coordinates,
        target_size_yxz=(ny, nx, nz),
    )
    n = int(loaded.coordinates.shape[0])
    n_in = int(in_bounds.sum())
    result.n_points = n
    result.n_in_bounds = n_in
    result.n_dropped = n - n_in
    if n == 0:
        result.messages.append("No points in CSV")
        return result
    frac = n_in / n
    if frac < min_in_bounds_fraction:
        result.messages.append(
            f"Only {100 * frac:.1f}% of points in bounds "
            f"(need ≥{100 * min_in_bounds_fraction:.0f}%); check axis order / index base"
        )
        return result
    if result.n_dropped:
        result.messages.append(f"Dropped {result.n_dropped} out-of-bounds point(s)")
    result.ok = True
    result.messages.append(f"{n_in}/{n} points in bounds on {reference.shape_tuple}")
    return result


def validate_mask_tiff(
    path: Path,
    reference: SampleReference,
    *,
    label: str = "",
) -> AnnotationValidationResult:
    """Check a 3D mask TIFF shape against ``sample_reference.shape_yxz``."""
    path = path.expanduser()
    stem = label or path.stem
    result = AnnotationValidationResult(
        label=stem,
        path=path,
        format=AnnotationFormat.MASK_TIFF,
        ok=False,
    )
    if not path.is_file():
        result.messages.append(f"File not found: {path}")
        return result

    with tifffile.TiffFile(path) as tif:
        if len(tif.pages) == 0:
            result.messages.append("Empty TIFF")
            return result
        first = tif.pages[0].asarray()
        if first.ndim == 2:
            ny, nx = first.shape
            nz = len(tif.pages)
            shape_yxz = (ny, nx, nz)
        elif first.ndim == 3:
            # single 3D page: assume ZYX
            shape_yxz = (int(first.shape[1]), int(first.shape[2]), int(first.shape[0]))
        else:
            result.messages.append(f"Unexpected TIFF ndim={first.ndim}")
            return result

    result.shape_yxz = shape_yxz
    expected = reference.shape_tuple
    if shape_yxz != expected:
        result.messages.append(f"Mask shape {shape_yxz} != sample_reference {expected}")
        return result
    result.ok = True
    result.messages.append(f"Mask shape matches {expected}")
    return result


def validate_annotation_specs(
    specs: list[AnnotationImportConfig],
    reference: SampleReference,
) -> list[AnnotationValidationResult]:
    """Validate every configured annotation layer."""
    results: list[AnnotationValidationResult] = []
    for spec in specs:
        label = spec.label or spec.path.stem
        if spec.format == AnnotationFormat.POINTS_CSV:
            results.append(validate_points_csv(spec.path, reference, label=label))
        elif spec.format == AnnotationFormat.MASK_TIFF:
            results.append(validate_mask_tiff(spec.path, reference, label=label))
        else:
            results.append(
                AnnotationValidationResult(
                    label=label,
                    path=spec.path,
                    format=spec.format,
                    ok=False,
                    messages=[f"Unsupported format: {spec.format}"],
                )
            )
    return results


def require_all_valid(results: list[AnnotationValidationResult]) -> None:
    """Raise ``ValueError`` if any validation failed."""
    failed = [r for r in results if not r.ok]
    if not failed:
        return
    parts = [f"{r.label} ({r.path}): {'; '.join(r.messages)}" for r in failed]
    msg = "Annotation validation failed:\n  - " + "\n  - ".join(parts)
    raise ValueError(msg)
