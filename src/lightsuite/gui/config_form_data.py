"""Form state and YAML round-trip helpers for the GUI config editor."""

from __future__ import annotations

import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml

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


def brain_form_from_raw(raw: dict[str, Any]) -> BrainFormState:
    sample = raw.get("sample") or {}
    source = sample.get("source") or {}
    atlas = raw.get("atlas") or {}
    registration = raw.get("registration") or {}
    compute = raw.get("compute") or {}
    detection = raw.get("detection") or {}
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
    return out


def spinal_form_from_raw(raw: dict[str, Any]) -> SpinalFormState:
    sample = raw.get("sample") or {}
    source = sample.get("source") or {}
    atlas = raw.get("atlas") or {}
    registration = raw.get("registration") or {}
    compute = raw.get("compute") or {}
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
    out["registration"] = registration

    compute = dict(out.get("compute") or {})
    compute["workers"] = state.workers
    out["compute"] = compute
    return out


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
    return MultiresFormState(
        sample_name=str(sample.get("name") or ""),
        save_path=_path_str(sample.get("save_path")),
        scratch=_path_str(sample.get("scratch")),
        pair_label=str(multires.get("pair_label") or ""),
        pair_manifest=_path_str(multires.get("pair_manifest")),
        reference_channel=str(registration.get("reference_channel") or ""),
        channels=channels,
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
    multires["registration"] = registration
    out["multires"] = multires
    return out


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
