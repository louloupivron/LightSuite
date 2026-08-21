"""Detect pipeline workflow type from YAML configs."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml

from lightsuite.config.loader import load_config, load_multires_config, load_spinal_config
from lightsuite.exceptions import LightsuiteConfigError


_MULTIRES_PIPELINE_KEYS = frozenset(
    {
        "pair_manifest",
        "channels",
        "pair_label",
        "vendor",
        "geometry_mode",
        "landmarks",
        "registration",
        "mesospim_geometry",
    }
)


def is_brain_multires_link(multires: dict[str, Any]) -> bool:
    """True when ``multires`` only points at an external multires YAML (brain combined workflow)."""
    if not isinstance(multires, dict) or not multires:
        return False
    config_ref = multires.get("config")
    if not isinstance(config_ref, str) or not config_ref.strip():
        return False
    return not any(key in multires for key in _MULTIRES_PIPELINE_KEYS)


def is_multires_pipeline_config(raw: dict[str, Any]) -> bool:
    """True when ``multires`` is a full multires workflow block (not a brain link)."""
    multires = raw.get("multires")
    if not isinstance(multires, dict) or not multires:
        return False
    return not is_brain_multires_link(multires)


def detect_workflow(raw: dict[str, Any]) -> str:
    """Infer workflow name from parsed YAML (brain, spinal, or multires)."""
    if is_multires_pipeline_config(raw):
        return "multires"
    atlas = raw.get("atlas") or {}
    if isinstance(atlas, dict) and atlas.get("atlas_dir") and "provider" not in atlas:
        return "spinal"
    return "brain"


def load_project(config_path: str | Path) -> tuple[str, Any]:
    """Load a config file and return ``(workflow, validated_config)``."""
    path = Path(config_path).expanduser().resolve()
    if not path.is_file():
        msg = f"Config file not found: {path}"
        raise FileNotFoundError(msg)

    with path.open(encoding="utf-8") as handle:
        raw: dict[str, Any] = yaml.safe_load(handle) or {}

    workflow = detect_workflow(raw)
    try:
        if workflow == "multires":
            return workflow, load_multires_config(path)
        if workflow == "spinal":
            return workflow, load_spinal_config(path)
        return workflow, load_config(path)
    except LightsuiteConfigError:
        raise
    except Exception as exc:
        msg = f"Failed to load {workflow} config from {path}: {exc}"
        raise LightsuiteConfigError(msg) from exc
