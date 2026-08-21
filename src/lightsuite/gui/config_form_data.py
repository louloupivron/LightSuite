"""Form state and YAML round-trip helpers for the GUI config editor."""

from __future__ import annotations

import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml

from lightsuite.analysis.intensity_metrics import DEFAULT_INTENSITY_METRICS, normalize_intensity_metrics
from lightsuite.atlas.brainglobe_backend import (
    brainglobe_name_for_yaml,
    find_brainglobe_catalog_entry,
    infer_brainglobe_resolution_um,
    list_lightsuite_brainglobe_atlases,
)
from lightsuite.cli.spaces import default_export_space_checks, export_spaces_from_checks
from lightsuite.config.workflow import detect_workflow, is_brain_multires_link, load_project
from lightsuite.exceptions import LightsuiteConfigError

_REPO_ROOT = Path(__file__).resolve().parents[3]

_CONFIG_TEMPLATES: dict[str, Path] = {
    "brain": _REPO_ROOT / "examples/brain_lightsheet.yaml",
    "spinal": _REPO_ROOT / "examples/spinal_cord.yaml",
    "multires": _REPO_ROOT / "examples/config/multiresolution/gilda_tg14_multires.yaml",
}


@dataclass
class ChannelPaths:
    name: str
    overview: str = ""
    roi: str = ""
    overview_meta: str = ""
    roi_meta: str = ""


@dataclass
class AnnotationImportRow:
    format: str = "points_csv"
    path: str = ""
    label: str = ""


@dataclass
class AnnotationConverterState:
    suite: str = "native"
    source: str = ""
    label: str = ""
    # Optional Imaris/FIJI calibration; None → use sample / sample_reference voxel size.
    voxel_um: tuple[float, float, float] | None = None
    custom_entry: str = ""


def import_segmentation_enabled_from_raw(raw: dict[str, Any]) -> bool:
    """True when the YAML has vendor conversion and/or native annotation layers."""
    import_block = raw.get("import") or {}
    if not isinstance(import_block, dict):
        return False
    annotations = import_block.get("annotations")
    if isinstance(annotations, list) and any(
        isinstance(item, dict) and str(item.get("path") or "").strip() for item in annotations
    ):
        return True
    conv = import_block.get("converter") or {}
    if not isinstance(conv, dict):
        return False
    suite = str(conv.get("suite") or "native").strip().lower()
    if suite != "native":
        return True
    return bool(str(conv.get("source") or "").strip() or str(conv.get("custom_entry") or "").strip())


def import_converter_from_raw(raw: dict[str, Any]) -> AnnotationConverterState:
    import_block = raw.get("import") or {}
    conv = import_block.get("converter") or {}
    if not isinstance(conv, dict):
        return AnnotationConverterState()
    voxel = conv.get("voxel_um")
    voxel_tuple: tuple[float, float, float] | None = None
    if isinstance(voxel, list) and len(voxel) == 3:
        try:
            voxel_tuple = (float(voxel[0]), float(voxel[1]), float(voxel[2]))
        except (TypeError, ValueError):
            voxel_tuple = None
    return AnnotationConverterState(
        suite=str(conv.get("suite") or "native"),
        source=_path_str(conv.get("source")),
        label=str(conv.get("label") or ""),
        voxel_um=voxel_tuple,
        custom_entry=_path_str(conv.get("custom_entry")),
    )


def import_converter_to_raw(
    state: AnnotationConverterState,
    raw: dict[str, Any],
) -> dict[str, Any]:
    """Merge converter fields into the raw YAML mapping (preserves annotations)."""
    out = dict(raw)
    import_block = dict(out.get("import") or {})
    suite = (state.suite or "native").strip().lower()
    empty = (
        suite == "native"
        and not state.source.strip()
        and not state.custom_entry.strip()
        and not state.label.strip()
        and state.voxel_um is None
    )
    if empty:
        import_block.pop("converter", None)
        if import_block:
            out["import"] = import_block
        else:
            out.pop("import", None)
        return out

    converter: dict[str, Any] = {"suite": suite}
    if state.source.strip():
        converter["source"] = state.source.strip()
    if state.label.strip():
        converter["label"] = state.label.strip()
    if state.custom_entry.strip():
        converter["custom_entry"] = state.custom_entry.strip()
    if state.voxel_um is not None:
        converter["voxel_um"] = [
            float(state.voxel_um[0]),
            float(state.voxel_um[1]),
            float(state.voxel_um[2]),
        ]
    # Drop legacy converter.output if present — destination is always under save_path/converted/.
    import_block["converter"] = converter
    import_block.setdefault("write_csv", True)
    out["import"] = import_block
    return out


def _parse_orientation(raw: Any) -> tuple[int, int, int] | None:
    if not isinstance(raw, list) or len(raw) != 3:
        return None
    try:
        return (int(raw[0]), int(raw[1]), int(raw[2]))
    except (TypeError, ValueError):
        return None


def _format_orientation(orientation: tuple[int, int, int] | None) -> str:
    if orientation is None:
        return ""
    return f"{orientation[0]}, {orientation[1]}, {orientation[2]}"


def parse_orientation_text(text: str) -> tuple[int, int, int] | None:
    """Parse ``1, 2, 3`` or ``[1, -3, 2]`` into an orientation triple."""
    cleaned = text.strip().replace("[", "").replace("]", "")
    if not cleaned:
        return None
    parts = [part.strip() for part in cleaned.split(",")]
    if len(parts) != 3:
        return None
    try:
        return (int(parts[0]), int(parts[1]), int(parts[2]))
    except ValueError:
        return None


def import_annotations_from_raw(raw: dict[str, Any]) -> list[AnnotationImportRow]:
    import_block = raw.get("import") or {}
    annotations = import_block.get("annotations")
    rows: list[AnnotationImportRow] = []
    if isinstance(annotations, list):
        for item in annotations:
            if not isinstance(item, dict):
                continue
            rows.append(
                AnnotationImportRow(
                    format=str(item.get("format") or "points_csv"),
                    path=_path_str(item.get("path")),
                    label=str(item.get("label") or ""),
                )
            )
    return rows


def import_annotations_to_raw(
    rows: list[AnnotationImportRow],
    raw: dict[str, Any],
) -> dict[str, Any]:
    out = dict(raw)
    import_block = dict(out.get("import") or {})
    annotations: list[dict[str, str]] = []
    for row in rows:
        path = row.path.strip()
        if not path:
            continue
        entry: dict[str, str] = {
            "format": row.format.strip() or "points_csv",
            "path": path,
        }
        if row.label.strip():
            entry["label"] = row.label.strip()
        annotations.append(entry)
    if annotations:
        import_block["annotations"] = annotations
        import_block.setdefault("write_csv", True)
        out["import"] = import_block
    else:
        import_block.pop("annotations", None)
        if import_block:
            out["import"] = import_block
        else:
            out.pop("import", None)
    return out


def _apply_import_to_raw(
    state_converter: AnnotationConverterState,
    state_annotations: list[AnnotationImportRow],
    raw: dict[str, Any],
    *,
    enabled: bool,
) -> dict[str, Any]:
    """Write import block when segmentation import is enabled; otherwise remove it."""
    if not enabled:
        out = dict(raw)
        out.pop("import", None)
        return out
    out = import_annotations_to_raw(state_annotations, raw)
    return import_converter_to_raw(state_converter, out)


_BRAIN_ATLAS_PROVIDERS_FILES = ("allen", "perens")
_BRAIN_ATLAS_PROVIDERS_BRAINGLOBE = ("allen", "perens", "princeton", "rat")

_BRAIN_ATLAS_PROVIDER_LABELS: dict[str, str] = {
    "allen": "Allen CCF (mouse)",
    "perens": "Perens / Gubra LSFM (mouse)",
    "princeton": "Princeton (mouse)",
    "rat": "Waxholm SD rat",
}

_BRAIN_ATLAS_SOURCE_LABELS: dict[str, str] = {
    "files": "Local NIfTI files",
    "brainglobe": "BrainGlobe (auto-download)",
}


def brain_atlas_provider_label(provider: str) -> str:
    """Human-readable label for a brain atlas provider id."""
    key = provider.lower().strip()
    return _BRAIN_ATLAS_PROVIDER_LABELS.get(key, key)


def brain_atlas_provider_options(source: str) -> list[tuple[str, str]]:
    """Return ``(provider_id, label)`` pairs available for an atlas source."""
    source_key = source.lower().strip()
    providers = (
        _BRAIN_ATLAS_PROVIDERS_BRAINGLOBE
        if source_key == "brainglobe"
        else _BRAIN_ATLAS_PROVIDERS_FILES
    )
    return [(provider, brain_atlas_provider_label(provider)) for provider in providers]


def brain_atlas_help_text(source: str, provider: str) -> str:
    """Short GUI help for the selected atlas source and provider."""
    source_key = source.lower().strip()
    provider_key = provider.lower().strip()

    if source_key == "brainglobe":
        return (
            "Choose a BrainGlobe atlas — resolution is taken from the registry "
            "(needs uv sync --extra atlas). Re-run check-orientation when switching "
            "from local NIfTIs."
        )

    if provider_key == "allen":
        return (
            "Point atlas directory at a folder with average_template_10.nii.gz and "
            "annotation_10.nii.gz. Optional: annotation_boundary_10.nii.gz and "
            "parcellation_to_parcellation_term_membership.csv."
        )
    if provider_key == "perens":
        return (
            "Point atlas directory at Gubra LSFM_atlas_files with gubra_template_olf.nii.gz "
            "and gubra_ano_olf.nii.gz. Optional: ARA2_annotation_info_avail_regions.csv."
        )
    return "Set atlas directory to a folder containing the template and annotation NIfTIs."


def default_local_atlas_resolution_um(provider: str) -> float:
    """Default native atlas voxel size for local NIfTI layouts."""
    key = provider.lower().strip()
    if key == "perens":
        return 20.0
    if key == "princeton":
        return 20.0
    if key == "rat":
        return 39.0
    return 10.0


def brain_brainglobe_catalog() -> list[Any]:
    """BrainGlobe atlases LightSuite supports (from registry or built-in defaults)."""
    return list_lightsuite_brainglobe_atlases()


def resolve_brain_brainglobe_form_fields(
    *,
    brainglobe_name: str = "",
    provider: str = "allen",
    resolution_um: float = 10.0,
    catalog: list[Any] | None = None,
) -> tuple[str, str, float]:
    """Return ``(brainglobe_name, provider, resolution_um)`` for the config form."""
    entries = catalog if catalog is not None else brain_brainglobe_catalog()
    entry = find_brainglobe_catalog_entry(
        entries,
        brainglobe_name=brainglobe_name,
        provider=provider,
        resolution_um=resolution_um,
    )
    if entry is not None:
        return entry.name, entry.provider, float(entry.resolution_um)
    inferred = infer_brainglobe_resolution_um(brainglobe_name) if brainglobe_name.strip() else None
    return (
        brainglobe_name,
        provider,
        float(inferred if inferred is not None else resolution_um),
    )


@dataclass
class BrainFormState:
    sample_name: str = ""
    source_path: str = ""
    channel_paths: list[str] = field(default_factory=list)
    use_channel_list: bool = False
    tiff_type: str = "channelperfile"
    scratch: str = ""
    save_path: str = ""
    voxel_um: tuple[float, float, float] = (5.0, 5.0, 5.0)
    atlas_source: str = "files"
    atlas_provider: str = "allen"
    brainglobe_name: str = ""
    atlas_resolution_um: float = 10.0
    atlas_dir: str = ""
    channel_primary: int = 1
    channel_secondary: int | None = None
    registration_resolution_um: float = 20.0
    bspline_spatial_scale_mm: float = 0.64
    control_point_weight: float = 0.2
    augment_points: bool = False
    dual_channel_mi_weight_primary: float = 0.4
    dual_channel_mi_weight_secondary: float = 0.4
    orientation: tuple[int, int, int] | None = None
    canvas_mode: str = "off"
    import_segmentation: bool = False
    import_converter: AnnotationConverterState = field(default_factory=AnnotationConverterState)
    import_annotations: list[AnnotationImportRow] = field(default_factory=list)
    intensity_metrics: list[str] = field(default_factory=lambda: list(DEFAULT_INTENSITY_METRICS))
    stats_spaces: list[str] = field(default_factory=lambda: ["atlas"])
    workers: int = 4
    detection_enabled: bool = False


@dataclass
class SpinalFormState:
    sample_name: str = ""
    source_path: str = ""
    channel_paths: list[str] = field(default_factory=list)
    use_channel_list: bool = True
    tiff_type: str = "planeperfile"
    scratch: str = ""
    save_path: str = ""
    voxel_um: tuple[float, float, float] = (1.8, 1.8, 1.8)
    atlas_dir: str = ""
    channel_primary: int = 1
    registration_resolution_um: float = 20.0
    control_point_weight: float = 0.2
    import_segmentation: bool = False
    import_converter: AnnotationConverterState = field(default_factory=AnnotationConverterState)
    import_annotations: list[AnnotationImportRow] = field(default_factory=list)
    intensity_metrics: list[str] = field(default_factory=lambda: list(DEFAULT_INTENSITY_METRICS))
    stats_spaces: list[str] = field(default_factory=lambda: ["atlas"])
    parcellate_intensities: bool = True
    workers: int = 4


@dataclass
class MultiresFormState:
    sample_name: str = ""
    save_path: str = ""
    scratch: str = ""
    pair_label: str = ""
    vendor_suite: str = "mesospim"
    vendor_custom_entry: str = ""
    pair_manifest: str = ""
    reference_channel: str = ""
    experiment_name: str = "default"
    channels: list[ChannelPaths] = field(default_factory=list)
    geometry_mode: str = "metadata"
    landmark_fit_mode: str = "similarity"
    overlap_margin_um: float = 0.0
    write_full_overview_canvas: bool = True
    geometry_check_level: str = "full"
    import_segmentation: bool = False
    import_converter: AnnotationConverterState = field(default_factory=AnnotationConverterState)
    import_annotations: list[AnnotationImportRow] = field(default_factory=list)


def _infer_multires_vendor_suite(multires: dict[str, Any]) -> str:
    vendor = multires.get("vendor") or {}
    if isinstance(vendor, dict):
        suite = str(vendor.get("suite") or "").strip().lower()
        if suite:
            return suite
        if str(vendor.get("custom_entry") or "").strip():
            return "custom"
    channels_raw = multires.get("channels") or {}
    has_channels = isinstance(channels_raw, dict) and bool(channels_raw)
    has_manifest = bool(str(multires.get("pair_manifest") or "").strip())
    if has_manifest and not has_channels:
        return "manifest"
    return "mesospim"


def _multires_vendor_from_raw(multires: dict[str, Any]) -> tuple[str, str]:
    vendor = multires.get("vendor") or {}
    suite = _infer_multires_vendor_suite(multires)
    custom_entry = ""
    if isinstance(vendor, dict):
        custom_entry = _path_str(vendor.get("custom_entry"))
    return suite, custom_entry


def _apply_multires_vendor_to_raw(
    *,
    suite: str,
    custom_entry: str,
    multires: dict[str, Any],
) -> dict[str, Any]:
    suite_key = (suite or "mesospim").strip().lower()
    vendor: dict[str, Any] = {"suite": suite_key}
    if suite_key == "custom" and custom_entry.strip():
        vendor["custom_entry"] = custom_entry.strip()
    multires["vendor"] = vendor
    return multires


def load_template_raw(workflow: str) -> tuple[str, dict[str, Any]]:
    """Load a starter YAML template as ``(workflow, raw_dict)``."""
    key = workflow.strip().lower()
    template_path = _CONFIG_TEMPLATES.get(key)
    if template_path is None or not template_path.is_file():
        available = ", ".join(sorted(_CONFIG_TEMPLATES))
        msg = f"Unknown workflow {workflow!r}; choose from {available}"
        raise ValueError(msg)
    raw = yaml.safe_load(template_path.read_text(encoding="utf-8")) or {}
    if not isinstance(raw, dict):
        msg = f"Template {template_path} is not a YAML mapping."
        raise ValueError(msg)
    return key, raw


def load_raw_config(path: Path, *, workflow_hint: str | None = None) -> tuple[str, dict[str, Any]]:
    """Read a config file and return ``(workflow, raw_dict)``."""
    resolved = path.expanduser().resolve()
    raw = yaml.safe_load(resolved.read_text(encoding="utf-8")) or {}
    if not isinstance(raw, dict) or not raw:
        msg = f"Config must be a non-empty YAML mapping: {resolved}"
        raise LightsuiteConfigError(msg)
    workflow = detect_workflow(raw)
    if workflow_hint in {"brain", "spinal", "multires"} and workflow != workflow_hint:
        multires = raw.get("multires")
        if (
            workflow_hint == "multires"
            and isinstance(multires, dict)
            and multires
            and not is_brain_multires_link(multires)
        ):
            workflow = "multires"
    return workflow, raw


def dump_config_dict(data: dict[str, Any]) -> str:
    """Serialize a config mapping to YAML text."""
    text = yaml.safe_dump(data, sort_keys=False, default_flow_style=False)
    return text if text.endswith("\n") else text + "\n"


def _path_str(value: Any) -> str:
    if value is None:
        return ""
    return str(value)


def _voxel_tuple(raw: Any, default: tuple[float, float, float]) -> tuple[float, float, float]:
    if not isinstance(raw, list) or len(raw) != 3:
        return default
    try:
        return (float(raw[0]), float(raw[1]), float(raw[2]))
    except (TypeError, ValueError):
        return default


def _analysis_intensity_metrics_from_raw(raw: dict[str, Any]) -> list[str]:
    analysis = raw.get("analysis") or {}
    metrics = analysis.get("intensity_metrics")
    if metrics is None:
        return list(DEFAULT_INTENSITY_METRICS)
    return normalize_intensity_metrics(metrics)


def normalize_stats_spaces(spaces: list[str] | None) -> list[str]:
    """Normalize analysis.stats_spaces to atlas and/or sample (default atlas)."""
    atlas, sample = default_export_space_checks(spaces)
    if not atlas and not sample:
        return ["atlas"]
    return export_spaces_from_checks(atlas=atlas, sample=sample)


def _analysis_stats_spaces_from_raw(raw: dict[str, Any]) -> list[str]:
    analysis = raw.get("analysis") or {}
    spaces = analysis.get("stats_spaces")
    if spaces is None:
        return ["atlas"]
    if isinstance(spaces, str):
        spaces = [part.strip() for part in spaces.split(",") if part.strip()]
    if not isinstance(spaces, list):
        return ["atlas"]
    return normalize_stats_spaces([str(item) for item in spaces])


def _apply_analysis_to_raw(
    raw: dict[str, Any],
    *,
    intensity_metrics: list[str],
    stats_spaces: list[str] | None = None,
    parcellate_intensities: bool | None = None,
) -> dict[str, Any]:
    out = dict(raw)
    analysis = dict(out.get("analysis") or {})
    analysis["intensity_metrics"] = normalize_intensity_metrics(intensity_metrics)
    analysis["count_points"] = True
    if stats_spaces is not None:
        analysis["stats_spaces"] = normalize_stats_spaces(stats_spaces)
    if parcellate_intensities is not None:
        analysis["parcellate_intensities"] = parcellate_intensities
    out["analysis"] = analysis
    return out


def _optional_registration_float(
    registration: dict[str, Any],
    key: str,
    *legacy_keys: str,
    default: float,
) -> float:
    if key in registration and registration[key] is not None:
        return float(registration[key])
    for legacy in legacy_keys:
        if legacy in registration and registration[legacy] is not None:
            return float(registration[legacy])
    return default


def brain_form_from_raw(raw: dict[str, Any]) -> BrainFormState:
    sample = raw.get("sample") or {}
    source = sample.get("source") or {}
    atlas = raw.get("atlas") or {}
    registration = raw.get("registration") or {}
    compute = raw.get("compute") or {}
    detection = raw.get("detection") or {}
    analysis = raw.get("analysis") or {}
    channels = source.get("channels")
    channel_paths = [_path_str(item) for item in channels] if isinstance(channels, list) else []
    secondary = registration.get("channel_secondary")
    atlas_source = str(atlas.get("source") or "files")
    atlas_provider = str(atlas.get("provider") or "allen")
    atlas_resolution_um = float(atlas.get("resolution_um") or 10.0)
    brainglobe_name = _path_str(atlas.get("brainglobe_name"))
    if atlas_source == "brainglobe":
        brainglobe_name, atlas_provider, atlas_resolution_um = resolve_brain_brainglobe_form_fields(
            brainglobe_name=brainglobe_name,
            provider=atlas_provider,
            resolution_um=atlas_resolution_um,
        )
    return BrainFormState(
        sample_name=str(sample.get("name") or ""),
        source_path=_path_str(source.get("path")),
        channel_paths=channel_paths,
        use_channel_list=bool(channel_paths),
        tiff_type=str(source.get("tiff_type") or "channelperfile"),
        scratch=_path_str(sample.get("scratch")),
        save_path=_path_str(sample.get("save_path")),
        voxel_um=_voxel_tuple(sample.get("voxel_um"), (5.0, 5.0, 5.0)),
        atlas_source=atlas_source,
        atlas_provider=atlas_provider,
        brainglobe_name=brainglobe_name,
        atlas_resolution_um=atlas_resolution_um,
        atlas_dir=_path_str(atlas.get("atlas_dir")),
        channel_primary=int(registration.get("channel_primary") or 1),
        channel_secondary=int(secondary) if secondary is not None else None,
        registration_resolution_um=float(registration.get("resolution_um") or 20.0),
        bspline_spatial_scale_mm=float(registration.get("bspline_spatial_scale_mm") or 0.64),
        control_point_weight=float(registration.get("control_point_weight") or 0.2),
        augment_points=bool(registration.get("augment_points", False)),
        dual_channel_mi_weight_primary=_optional_registration_float(
            registration,
            "dual_channel_mi_weight_primary",
            "dual_channel_mi_weight_autofluor",
            default=0.4,
        ),
        dual_channel_mi_weight_secondary=_optional_registration_float(
            registration,
            "dual_channel_mi_weight_secondary",
            "dual_channel_mi_weight_signal",
            default=0.4,
        ),
        orientation=_parse_orientation(registration.get("orientation")),
        canvas_mode=str(registration.get("canvas_mode") or "off").lower(),
        import_segmentation=import_segmentation_enabled_from_raw(raw),
        import_converter=import_converter_from_raw(raw),
        import_annotations=import_annotations_from_raw(raw),
        intensity_metrics=_analysis_intensity_metrics_from_raw(raw),
        stats_spaces=_analysis_stats_spaces_from_raw(raw),
        workers=int(compute.get("workers") or 4),
        detection_enabled=bool(detection.get("enabled", False)),
    )


def brain_form_to_raw(state: BrainFormState, raw: dict[str, Any]) -> dict[str, Any]:
    out = dict(raw)
    sample = dict(out.get("sample") or {})
    source = dict(sample.get("source") or {})
    source["format"] = source.get("format") or "tiff_stack"
    source["tiff_type"] = state.tiff_type
    if state.use_channel_list:
        source["channels"] = [path for path in state.channel_paths if path.strip()]
        if source["channels"]:
            source["path"] = source["channels"][0]
        else:
            source.pop("path", None)
    else:
        source.pop("channels", None)
        if state.source_path.strip():
            source["path"] = state.source_path.strip()
        else:
            source.pop("path", None)
    sample["source"] = source
    sample["name"] = state.sample_name.strip() or "sample"
    sample["scratch"] = state.scratch.strip()
    sample["save_path"] = state.save_path.strip()
    sample["voxel_um"] = [state.voxel_um[0], state.voxel_um[1], state.voxel_um[2]]
    out["sample"] = sample

    atlas = dict(out.get("atlas") or {})
    atlas["provider"] = state.atlas_provider
    if state.atlas_source.strip().lower() == "brainglobe":
        atlas["source"] = "brainglobe"
        bg_name, provider, resolution_um = resolve_brain_brainglobe_form_fields(
            brainglobe_name=state.brainglobe_name,
            provider=state.atlas_provider,
            resolution_um=state.atlas_resolution_um,
        )
        atlas["provider"] = provider
        atlas["resolution_um"] = resolution_um
        yaml_name = brainglobe_name_for_yaml(bg_name)
        if yaml_name:
            atlas["brainglobe_name"] = yaml_name
        else:
            atlas.pop("brainglobe_name", None)
    else:
        atlas.pop("source", None)
        atlas.pop("brainglobe_name", None)
        atlas["resolution_um"] = state.atlas_resolution_um
    if state.atlas_dir.strip():
        atlas["atlas_dir"] = state.atlas_dir.strip()
    else:
        atlas.pop("atlas_dir", None)
    out["atlas"] = atlas

    registration = dict(out.get("registration") or {})
    registration["channel_primary"] = state.channel_primary
    registration["resolution_um"] = state.registration_resolution_um
    registration["bspline_spatial_scale_mm"] = state.bspline_spatial_scale_mm
    registration["control_point_weight"] = state.control_point_weight
    registration["augment_points"] = state.augment_points
    registration["dual_channel_mi_weight_primary"] = state.dual_channel_mi_weight_primary
    registration["dual_channel_mi_weight_secondary"] = state.dual_channel_mi_weight_secondary
    registration.pop("dual_channel_mi_weight_autofluor", None)
    registration.pop("dual_channel_mi_weight_signal", None)
    if state.orientation is None:
        registration.pop("orientation", None)
    else:
        registration["orientation"] = list(state.orientation)
    if state.canvas_mode and state.canvas_mode != "off":
        registration["canvas_mode"] = state.canvas_mode
    else:
        registration.pop("canvas_mode", None)
    if state.channel_secondary is None:
        registration["channel_secondary"] = None
    else:
        registration["channel_secondary"] = state.channel_secondary
    out["registration"] = registration

    compute = dict(out.get("compute") or {})
    compute["workers"] = state.workers
    out["compute"] = compute

    detection = dict(out.get("detection") or {})
    detection["enabled"] = state.detection_enabled
    out["detection"] = detection
    out = _apply_analysis_to_raw(
        out,
        intensity_metrics=state.intensity_metrics,
        stats_spaces=state.stats_spaces,
    )
    return _apply_import_to_raw(
        state.import_converter,
        state.import_annotations,
        out,
        enabled=state.import_segmentation,
    )


def spinal_form_from_raw(raw: dict[str, Any]) -> SpinalFormState:
    sample = raw.get("sample") or {}
    source = sample.get("source") or {}
    atlas = raw.get("atlas") or {}
    registration = raw.get("registration") or {}
    compute = raw.get("compute") or {}
    analysis = raw.get("analysis") or {}
    channels = source.get("channels")
    channel_paths = [_path_str(item) for item in channels] if isinstance(channels, list) else []
    return SpinalFormState(
        sample_name=str(sample.get("name") or ""),
        source_path=_path_str(source.get("path")),
        channel_paths=channel_paths,
        use_channel_list=bool(channel_paths),
        tiff_type=str(source.get("tiff_type") or "planeperfile"),
        scratch=_path_str(sample.get("scratch")),
        save_path=_path_str(sample.get("save_path")),
        voxel_um=_voxel_tuple(sample.get("voxel_um"), (1.8, 1.8, 1.8)),
        atlas_dir=_path_str(atlas.get("atlas_dir")),
        channel_primary=int(registration.get("channel_primary") or 1),
        registration_resolution_um=float(registration.get("resolution_um") or 20.0),
        control_point_weight=float(registration.get("control_point_weight") or 0.2),
        import_segmentation=import_segmentation_enabled_from_raw(raw),
        import_converter=import_converter_from_raw(raw),
        import_annotations=import_annotations_from_raw(raw),
        intensity_metrics=_analysis_intensity_metrics_from_raw(raw),
        stats_spaces=_analysis_stats_spaces_from_raw(raw),
        parcellate_intensities=bool(analysis.get("parcellate_intensities", True)),
        workers=int(compute.get("workers") or 4),
    )


def spinal_form_to_raw(state: SpinalFormState, raw: dict[str, Any]) -> dict[str, Any]:
    out = dict(raw)
    sample = dict(out.get("sample") or {})
    source = dict(sample.get("source") or {})
    source["format"] = source.get("format") or "tiff_stack"
    source["tiff_type"] = state.tiff_type
    if state.use_channel_list:
        source["channels"] = [path for path in state.channel_paths if path.strip()]
        if source["channels"]:
            source["path"] = source["channels"][0]
        else:
            source.pop("path", None)
    else:
        source.pop("channels", None)
        if state.source_path.strip():
            source["path"] = state.source_path.strip()
        else:
            source.pop("path", None)
    sample["source"] = source
    sample["name"] = state.sample_name.strip() or "sample"
    sample["scratch"] = state.scratch.strip()
    sample["save_path"] = state.save_path.strip()
    sample["voxel_um"] = [state.voxel_um[0], state.voxel_um[1], state.voxel_um[2]]
    out["sample"] = sample

    atlas = dict(out.get("atlas") or {})
    if state.atlas_dir.strip():
        atlas["atlas_dir"] = state.atlas_dir.strip()
    out["atlas"] = atlas

    registration = dict(out.get("registration") or {})
    registration["channel_primary"] = state.channel_primary
    registration["resolution_um"] = state.registration_resolution_um
    registration["control_point_weight"] = state.control_point_weight
    out["registration"] = registration

    compute = dict(out.get("compute") or {})
    compute["workers"] = state.workers
    out["compute"] = compute
    out = _apply_analysis_to_raw(
        out,
        intensity_metrics=state.intensity_metrics,
        stats_spaces=state.stats_spaces,
        parcellate_intensities=state.parcellate_intensities,
    )
    return _apply_import_to_raw(
        state.import_converter,
        state.import_annotations,
        out,
        enabled=state.import_segmentation,
    )


def multires_form_from_raw(raw: dict[str, Any]) -> MultiresFormState:
    sample = raw.get("sample") or {}
    multires = raw.get("multires") or {}
    registration = multires.get("registration") or {}
    channels_raw = multires.get("channels") or {}
    channels: list[ChannelPaths] = []
    if isinstance(channels_raw, dict):
        for name, item in channels_raw.items():
            if not isinstance(item, dict):
                continue
            channels.append(
                ChannelPaths(
                    name=str(name),
                    overview=_path_str(item.get("overview")),
                    roi=_path_str(item.get("roi")),
                    overview_meta=_path_str(item.get("overview_meta_path")),
                    roi_meta=_path_str(item.get("roi_meta_path")),
                )
            )
    landmarks = multires.get("landmarks") or {}
    vendor_suite, vendor_custom_entry = _multires_vendor_from_raw(multires)
    return MultiresFormState(
        sample_name=str(sample.get("name") or ""),
        save_path=_path_str(sample.get("save_path")),
        scratch=_path_str(sample.get("scratch")),
        pair_label=str(multires.get("pair_label") or ""),
        vendor_suite=vendor_suite,
        vendor_custom_entry=vendor_custom_entry,
        pair_manifest=_path_str(multires.get("pair_manifest")),
        reference_channel=str(registration.get("reference_channel") or ""),
        experiment_name=str(registration.get("experiment_name") or "default"),
        channels=channels,
        geometry_mode=str(multires.get("geometry_mode") or "metadata"),
        landmark_fit_mode=str(landmarks.get("fit_mode") or "similarity"),
        overlap_margin_um=float(registration.get("overlap_margin_um") or 0.0),
        write_full_overview_canvas=bool(registration.get("write_full_overview_canvas", True)),
        geometry_check_level=str(registration.get("geometry_check_level") or "full"),
        import_segmentation=import_segmentation_enabled_from_raw(raw),
        import_converter=import_converter_from_raw(raw),
        import_annotations=import_annotations_from_raw(raw),
    )


def multires_form_to_raw(state: MultiresFormState, raw: dict[str, Any]) -> dict[str, Any]:
    out = dict(raw)
    sample = dict(out.get("sample") or {})
    sample["name"] = state.sample_name.strip() or "sample"
    sample["save_path"] = state.save_path.strip()
    if state.scratch.strip():
        sample["scratch"] = state.scratch.strip()
    else:
        sample.pop("scratch", None)
    out["sample"] = sample

    multires = dict(out.get("multires") or {})
    if state.pair_label.strip():
        multires["pair_label"] = state.pair_label.strip()
    else:
        multires.pop("pair_label", None)
    if state.pair_manifest.strip() and state.vendor_suite == "manifest":
        multires["pair_manifest"] = state.pair_manifest.strip()
    else:
        # Built-in / custom vendors write under save_path/converted/; drop template leftovers.
        multires.pop("pair_manifest", None)

    multires = _apply_multires_vendor_to_raw(
        suite=state.vendor_suite,
        custom_entry=state.vendor_custom_entry,
        multires=multires,
    )
    if state.vendor_suite.strip().lower() != "mesospim":
        multires.pop("mesospim_geometry", None)

    multires["geometry_mode"] = state.geometry_mode
    landmarks = dict(multires.get("landmarks") or {})
    landmarks["fit_mode"] = state.landmark_fit_mode
    multires["landmarks"] = landmarks

    channels: dict[str, dict[str, str]] = {}
    for item in state.channels:
        name = item.name.strip()
        if not name:
            continue
        entry: dict[str, str] = {}
        if item.overview.strip():
            entry["overview"] = item.overview.strip()
        if item.roi.strip():
            entry["roi"] = item.roi.strip()
        if item.overview_meta.strip():
            entry["overview_meta_path"] = item.overview_meta.strip()
        if item.roi_meta.strip():
            entry["roi_meta_path"] = item.roi_meta.strip()
        if entry:
            channels[name] = entry
    if channels:
        multires["channels"] = channels
    else:
        multires.pop("channels", None)

    registration = dict(multires.get("registration") or {})
    if state.reference_channel.strip():
        registration["reference_channel"] = state.reference_channel.strip()
    else:
        registration.pop("reference_channel", None)
    registration["overlap_margin_um"] = state.overlap_margin_um
    registration["write_full_overview_canvas"] = state.write_full_overview_canvas
    level = (state.geometry_check_level or "full").strip().lower()
    if level and level != "full":
        registration["geometry_check_level"] = level
    else:
        registration.pop("geometry_check_level", None)
    experiment_name = state.experiment_name.strip() or "default"
    if experiment_name != "default":
        registration["experiment_name"] = experiment_name
    else:
        registration.pop("experiment_name", None)
    # Template leftovers (e.g. Gilda apply_transform_to: [555, 647]) must not survive
    # when the form only keeps a subset of channels.
    apply_raw = registration.get("apply_transform_to")
    if isinstance(apply_raw, list):
        kept = [str(name) for name in apply_raw if str(name) in channels]
        if kept:
            registration["apply_transform_to"] = kept
        else:
            registration.pop("apply_transform_to", None)
    multires["registration"] = registration
    out["multires"] = multires
    return _apply_import_to_raw(
        state.import_converter,
        state.import_annotations,
        out,
        enabled=state.import_segmentation,
    )


def try_validate_config_dict(
    data: dict[str, Any],
    *,
    workflow: str | None = None,
) -> tuple[str, Any] | LightsuiteConfigError:
    """Validate a config mapping via the normal loaders."""
    text = dump_config_dict(data)
    with tempfile.NamedTemporaryFile(
        mode="w",
        suffix=".yaml",
        encoding="utf-8",
        delete=False,
    ) as handle:
        handle.write(text)
        temp_path = Path(handle.name)
    try:
        detected = detect_workflow(data)
        use_workflow = workflow if workflow in {"brain", "spinal", "multires"} else detected
        multires = data.get("multires")
        if (
            workflow == "multires"
            and isinstance(multires, dict)
            and multires
            and not is_brain_multires_link(multires)
        ):
            use_workflow = "multires"
        if use_workflow == detected:
            return load_project(temp_path)
        from lightsuite.config.loader import load_config, load_multires_config, load_spinal_config

        if use_workflow == "multires":
            return "multires", load_multires_config(temp_path)
        if use_workflow == "spinal":
            return "spinal", load_spinal_config(temp_path)
        return "brain", load_config(temp_path)
    except (LightsuiteConfigError, FileNotFoundError, OSError, ValueError) as exc:
        if isinstance(exc, LightsuiteConfigError):
            return exc
        return LightsuiteConfigError(str(exc))
    finally:
        temp_path.unlink(missing_ok=True)
