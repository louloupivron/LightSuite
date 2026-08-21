"""Load a user-provided pair-manifest conversion script for vendor.suite=custom."""

from __future__ import annotations

import importlib.util
from pathlib import Path
from types import ModuleType

from lightsuite.multires.config_models import MultiresPipelineConfig
from lightsuite.multires.manifest import save_pair_manifest
from lightsuite.multires.models import MANIFEST_FORMAT, MultiresPairManifest

CUSTOM_ENTRY_FN = "build_pair_manifest"


def load_custom_manifest_converter_module(path: Path) -> ModuleType:
    """Import ``path`` as a module; must define ``build_pair_manifest``."""
    path = path.expanduser().resolve()
    if not path.is_file():
        msg = f"Custom pair manifest converter not found: {path}"
        raise FileNotFoundError(msg)
    if path.suffix.lower() != ".py":
        msg = f"Custom pair manifest converter must be a .py file: {path}"
        raise ValueError(msg)

    spec = importlib.util.spec_from_file_location("lightsuite_custom_pair_manifest_converter", path)
    if spec is None or spec.loader is None:
        msg = f"Could not load custom pair manifest converter: {path}"
        raise ImportError(msg)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    if not hasattr(module, CUSTOM_ENTRY_FN):
        msg = (
            f"{path} must define {CUSTOM_ENTRY_FN}(cfg, output) -> MultiresPairManifest | dict. "
            "See examples/multires/custom_manifest_converter_example.py."
        )
        raise AttributeError(msg)
    return module


def _coerce_manifest(result: object) -> MultiresPairManifest:
    if isinstance(result, MultiresPairManifest):
        return result
    if isinstance(result, dict):
        loaded = MultiresPairManifest.from_dict(result)
        if loaded.format != MANIFEST_FORMAT:
            msg = f"Custom converter manifest format must be {MANIFEST_FORMAT!r}, got {loaded.format!r}"
            raise ValueError(msg)
        return loaded
    msg = f"{CUSTOM_ENTRY_FN} must return MultiresPairManifest or dict, got {type(result).__name__}"
    raise TypeError(msg)


def run_custom_manifest_converter(
    entry: Path,
    *,
    cfg: MultiresPipelineConfig,
    output: Path,
) -> MultiresPairManifest:
    """Call ``build_pair_manifest`` and persist the manifest to ``output``."""
    module = load_custom_manifest_converter_module(entry)
    fn = getattr(module, CUSTOM_ENTRY_FN)
    result = fn(cfg, Path(output))
    manifest = _coerce_manifest(result)
    save_pair_manifest(manifest, output)
    return manifest
