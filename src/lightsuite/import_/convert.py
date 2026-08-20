"""Convert vendor segmentation exports to LightSuite Sample Space (prepare stage)."""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from rich.console import Console

from lightsuite.config.models import (
    AnnotationConverterConfig,
    AnnotationFormat,
    AnnotationImportConfig,
    ImportConfig,
    SegmentationSuite,
)
from lightsuite.import_.arivis import convert_arivis_features_to_points_csv
from lightsuite.import_.custom import run_custom_converter
from lightsuite.import_.fiji import convert_fiji_points_to_csv
from lightsuite.import_.imaris import (
    convert_imaris_spots_to_points_csv,
    list_imaris_component_names,
)
from lightsuite.import_.sample_reference import SampleReference, load_sample_reference
from lightsuite.import_.smartspim_detection import convert_smartspim_points_json_to_csv
from lightsuite.import_.validate import (
    AnnotationValidationResult,
    require_all_valid,
    validate_annotation_specs,
)

console = Console()

CONVERT_SUMMARY_NAME = "convert_annotations_summary.json"


@dataclass
class ConvertAnnotationsResult:
    """Outcome of the convert-annotations prepare stage."""

    suite: str
    output_path: Path | None
    annotations: list[AnnotationImportConfig]
    validations: list[AnnotationValidationResult] = field(default_factory=list)
    n_converted: int = 0
    summary_path: Path | None = None


def _label_for(converter: AnnotationConverterConfig) -> str:
    if converter.label.strip():
        return converter.label.strip()
    if converter.source is not None:
        return converter.source.stem
    return "annotation"


def default_converter_output(save_path: Path, label: str) -> Path:
    return save_path.expanduser() / "converted" / f"{label}_points.csv"


def resolve_converter_output(
    converter: AnnotationConverterConfig,
    save_path: Path,
) -> Path:
    if converter.output is not None:
        return converter.output.expanduser()
    return default_converter_output(save_path, _label_for(converter))


def _voxel_um_or_reference(
    converter: AnnotationConverterConfig,
    reference: SampleReference,
) -> list[float]:
    if converter.voxel_um is not None:
        return [float(v) for v in converter.voxel_um]
    return [float(v) for v in reference.voxel_um]


def _slug_component(name: str) -> str:
    return "".join(ch if ch.isalnum() or ch in "-_" else "_" for ch in name)


def _run_vendor_convert(
    converter: AnnotationConverterConfig,
    *,
    output: Path,
    reference: SampleReference,
) -> tuple[list[AnnotationImportConfig], int]:
    """Run vendor conversion; return annotation specs written and total item count."""
    suite = converter.suite
    source = converter.source
    if source is None:
        msg = "converter.source is required"
        raise ValueError(msg)

    label = _label_for(converter)

    if suite == SegmentationSuite.SMARTSPIM:
        n = convert_smartspim_points_json_to_csv(source, output)
        return [AnnotationImportConfig(format=AnnotationFormat.POINTS_CSV, path=output, label=label)], n
    if suite == SegmentationSuite.FIJI:
        n = convert_fiji_points_to_csv(
            source,
            output,
            voxel_um=_voxel_um_or_reference(converter, reference),
        )
        return [AnnotationImportConfig(format=AnnotationFormat.POINTS_CSV, path=output, label=label)], n
    if suite == SegmentationSuite.IMARIS:
        voxel_um = _voxel_um_or_reference(converter, reference)
        components = list_imaris_component_names(source)
        if len(components) > 1:
            # Multi-spot Imaris: one Sample Space CSV per Component Name (spinal multi-channel).
            prefix = label or "imaris"
            output.parent.mkdir(parents=True, exist_ok=True)
            annotations: list[AnnotationImportConfig] = []
            n_total = 0
            for component in components:
                slug = _slug_component(component)
                layer_label = f"{prefix}_{slug}"
                out_path = output.parent / f"{layer_label}_points.csv"
                n_comp = convert_imaris_spots_to_points_csv(
                    source,
                    out_path,
                    voxel_um=voxel_um,
                    component_name=component,
                )
                if n_comp <= 0:
                    continue
                n_total += n_comp
                annotations.append(
                    AnnotationImportConfig(
                        format=AnnotationFormat.POINTS_CSV,
                        path=out_path,
                        label=layer_label,
                    )
                )
            return annotations, n_total
        n = convert_imaris_spots_to_points_csv(source, output, voxel_um=voxel_um)
        return [AnnotationImportConfig(format=AnnotationFormat.POINTS_CSV, path=output, label=label)], n
    if suite == SegmentationSuite.ARIVIS:
        n = convert_arivis_features_to_points_csv(source, output)
        return [AnnotationImportConfig(format=AnnotationFormat.POINTS_CSV, path=output, label=label)], n
    if suite == SegmentationSuite.CUSTOM:
        if converter.custom_entry is None:
            msg = "suite=custom requires converter.custom_entry"
            raise ValueError(msg)
        if converter.source is None:
            msg = "suite=custom requires converter.source"
            raise ValueError(msg)
        meta = run_custom_converter(
            converter.custom_entry,
            source=source,
            output=output,
            reference=reference,
        )
        fmt_raw = str(meta.get("format") or "points_csv")
        try:
            fmt = AnnotationFormat(fmt_raw)
        except ValueError as exc:
            msg = f"custom converter returned unknown format {fmt_raw!r}"
            raise ValueError(msg) from exc
        n = int(meta.get("n_points") or meta.get("n") or 0)
        return [AnnotationImportConfig(format=fmt, path=output, label=label)], n

    msg = f"Unsupported suite for conversion: {suite}"
    raise ValueError(msg)


def _write_summary(path: Path, payload: dict[str, Any]) -> Path:
    path = path.expanduser()
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return path


def run_convert_annotations(
    *,
    import_config: ImportConfig | None,
    save_path: Path,
    sample_voxel_um: list[float] | None = None,
) -> ConvertAnnotationsResult:
    """Convert (if configured) then validate annotations against ``sample_reference``.

    ``suite=native`` (or no converter) only validates ``import.annotations``.
    ``suite=custom`` always validates after the user script runs.
    """
    save_path = save_path.expanduser()
    reference = load_sample_reference(save_path)
    if sample_voxel_um is not None:
        # Prefer sample YAML voxel size when provided for Imaris/FIJI defaults
        pass

    if import_config is None:
        msg = (
            "No import section in config. Add import.annotations and/or import.converter "
            "in the YAML / Config panel."
        )
        raise ValueError(msg)

    converter = import_config.converter
    suite = converter.suite if converter is not None else SegmentationSuite.NATIVE

    annotations = list(import_config.annotations)
    n_converted = 0
    output_path: Path | None = None

    if converter is not None and suite != SegmentationSuite.NATIVE:
        output_path = resolve_converter_output(converter, save_path)
        annotations, n_converted = _run_vendor_convert(
            converter, output=output_path, reference=reference
        )
        console.print(
            f"[green]Converted[/green] {suite.value} → {len(annotations)} layer(s) "
            f"({n_converted} item(s))"
        )
        for spec in annotations:
            console.print(f"  • {spec.label or spec.path.stem}: {spec.path}")
    elif not annotations:
        msg = (
            "Nothing to convert or validate. Set import.converter.suite to a vendor "
            "(smartspim / fiji / imaris / arivis / custom) with source, or add "
            "import.annotations pointing at existing points_csv / mask_tiff."
        )
        raise ValueError(msg)
    else:
        console.print(
            f"[cyan]suite=native[/cyan] — validating {len(annotations)} annotation layer(s)"
        )

    validations = validate_annotation_specs(annotations, reference)
    for item in validations:
        status = "OK" if item.ok else "FAIL"
        console.print(f"  [{status}] {item.label}: {'; '.join(item.messages)}")
    require_all_valid(validations)

    if suite == SegmentationSuite.CUSTOM and converter is not None:
        console.print("[green]Custom converter output passed Sample Space validation[/green]")

    summary = {
        "suite": suite.value,
        "output": str(output_path) if output_path else None,
        "n_converted": n_converted,
        "annotations": [
            {"format": a.format.value, "path": str(a.path), "label": a.label or a.path.stem}
            for a in annotations
        ],
        "validations": [
            {
                "label": v.label,
                "ok": v.ok,
                "path": str(v.path),
                "format": v.format.value,
                "n_points": v.n_points,
                "n_in_bounds": v.n_in_bounds,
                "n_dropped": v.n_dropped,
                "shape_yxz": list(v.shape_yxz) if v.shape_yxz else None,
                "messages": v.messages,
            }
            for v in validations
        ],
    }
    summary_path = _write_summary(save_path / "converted" / CONVERT_SUMMARY_NAME, summary)

    return ConvertAnnotationsResult(
        suite=suite.value,
        output_path=output_path,
        annotations=annotations,
        validations=validations,
        n_converted=n_converted,
        summary_path=summary_path,
    )


def run_convert_annotations_for_pipeline(config: Any) -> ConvertAnnotationsResult:
    """Brain / spinal / multires entrypoint using ``config.sample.save_path``."""
    save_path = config.sample.save_path.expanduser()
    voxel = getattr(config.sample, "voxel_um", None)
    return run_convert_annotations(
        import_config=getattr(config, "import_config", None),
        save_path=save_path,
        sample_voxel_um=list(voxel) if voxel is not None else None,
    )
