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
from lightsuite.config.workflow import detect_workflow, load_project
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


@dataclass
class AnnotationImportRow:
    format: str = "points_csv"
    path: str = ""
    label: str = ""


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
    suite = (state_converter.suite or "native").strip().lower()
    if suite != "native":
        out = dict(raw)
        import_block = dict(out.get("import") or {})
        import_block.pop("annotations", None)
        out["import"] = import_block
        return import_converter_to_raw(state_converter, out)
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
    import_annotations: list[AnnotationImportRow] = field(default_factory=list)
    intensity_metrics: list[str] = field(default_factory=lambda: list(DEFAULT_INTENSITY_METRICS))
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
    import_annotations: list[AnnotationImportRow] = field(default_factory=list)
    intensity_metrics: list[str] = field(default_factory=lambda: list(DEFAULT_INTENSITY_METRICS))
    parcellate_intensities: bool = True
    workers: int = 4


@dataclass
class MultiresFormState:
    sample_name: str = ""
    save_path: str = ""
    scratch: str = ""
    pair_label: str = ""
    pair_manifest: str = ""
    reference_channel: str = ""
    channels: list[ChannelPaths] = field(default_factory=list)
    geometry_mode: str = "metadata"
    landmark_fit_mode: str = "similarity"
    overlap_margin_um: float = 0.0
    write_full_overview_canvas: bool = True
    import_annotations: list[AnnotationImportRow] = field(default_factory=list)


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


def load_raw_config(path: Path) -> tuple[str, dict[str, Any]]:
    """Read a config file and return ``(workflow, raw_dict)``."""
    resolved = path.expanduser().resolve()
    raw = yaml.safe_load(resolved.read_text(encoding="utf-8")) or {}
    if not isinstance(raw, dict) or not raw:
        msg = f"Config must be a non-empty YAML mapping: {resolved}"
        raise LightsuiteConfigError(msg)
    return detect_workflow(raw), raw


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


def _apply_analysis_to_raw(
    raw: dict[str, Any],
    *,
    intensity_metrics: list[str],
    parcellate_intensities: bool | None = None,
) -> dict[str, Any]:
    out = dict(raw)
    analysis = dict(out.get("analysis") or {})
    analysis["intensity_metrics"] = normalize_intensity_metrics(intensity_metrics)
    analysis["count_points"] = True
    if parcellate_intensities is not None:
        analysis["parcellate_intensities"] = parcellate_intensities
    out["analysis"] = analysis
    return out


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
        dual_channel_mi_weight_primary=float(
            registration.get("dual_channel_mi_weight_primary")
            or registration.get("dual_channel_mi_weight_autofluor")
            or 0.4
        ),
        dual_channel_mi_weight_secondary=float(
            registration.get("dual_channel_mi_weight_secondary")
            or registration.get("dual_channel_mi_weight_signal")
            or 0.4
        ),
        orientation=_parse_orientation(registration.get("orientation")),
        canvas_mode=str(registration.get("canvas_mode") or "off").lower(),
        import_annotations=import_annotations_from_raw(raw),
        intensity_metrics=_analysis_intensity_metrics_from_raw(raw),
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
    )
    return import_annotations_to_raw(state.import_annotations, out)


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
        import_annotations=import_annotations_from_raw(raw),
        intensity_metrics=_analysis_intensity_metrics_from_raw(raw),
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
        parcellate_intensities=state.parcellate_intensities,
    )
    return import_annotations_to_raw(state.import_annotations, out)


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
                )
            )
    landmarks = multires.get("landmarks") or {}
    return MultiresFormState(
        sample_name=str(sample.get("name") or ""),
        save_path=_path_str(sample.get("save_path")),
        scratch=_path_str(sample.get("scratch")),
        pair_label=str(multires.get("pair_label") or ""),
        pair_manifest=_path_str(multires.get("pair_manifest")),
        reference_channel=str(registration.get("reference_channel") or ""),
        channels=channels,
        geometry_mode=str(multires.get("geometry_mode") or "metadata"),
        landmark_fit_mode=str(landmarks.get("fit_mode") or "similarity"),
        overlap_margin_um=float(registration.get("overlap_margin_um") or 0.0),
        write_full_overview_canvas=bool(registration.get("write_full_overview_canvas", True)),
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
    if state.pair_manifest.strip():
        multires["pair_manifest"] = state.pair_manifest.strip()
    else:
        multires.pop("pair_manifest", None)

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
    multires["registration"] = registration
    out["multires"] = multires
    return import_annotations_to_raw(state.import_annotations, out)


def try_validate_config_dict(data: dict[str, Any]) -> tuple[str, Any] | LightsuiteConfigError:
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
        return load_project(temp_path)
    except (LightsuiteConfigError, FileNotFoundError, OSError, ValueError) as exc:
        if isinstance(exc, LightsuiteConfigError):
            return exc
        return LightsuiteConfigError(str(exc))
    finally:
        temp_path.unlink(missing_ok=True)
