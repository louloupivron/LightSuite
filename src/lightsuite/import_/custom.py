"""Load a user-provided conversion script for suite=custom."""

from __future__ import annotations

import importlib.util
from pathlib import Path
from types import ModuleType
from typing import Any

from lightsuite.import_.sample_reference import SampleReference

CUSTOM_ENTRY_FN = "convert_to_lightsuite"


def load_custom_converter_module(path: Path) -> ModuleType:
    """Import ``path`` as a module; must define ``convert_to_lightsuite``."""
    path = path.expanduser().resolve()
    if not path.is_file():
        msg = f"Custom converter not found: {path}"
        raise FileNotFoundError(msg)
    if path.suffix.lower() != ".py":
        msg = f"Custom converter must be a .py file: {path}"
        raise ValueError(msg)

    spec = importlib.util.spec_from_file_location("lightsuite_custom_annotation_converter", path)
    if spec is None or spec.loader is None:
        msg = f"Could not load custom converter: {path}"
        raise ImportError(msg)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    if not hasattr(module, CUSTOM_ENTRY_FN):
        msg = (
            f"{path} must define {CUSTOM_ENTRY_FN}(source, output, *, reference) -> dict. "
            "See docs/annotation_import.md (custom converters)."
        )
        raise AttributeError(msg)
    return module


def run_custom_converter(
    entry: Path,
    *,
    source: Path,
    output: Path,
    reference: SampleReference,
) -> dict[str, Any]:
    """Call ``convert_to_lightsuite`` and return its result dict."""
    module = load_custom_converter_module(entry)
    fn = getattr(module, CUSTOM_ENTRY_FN)
    result = fn(Path(source), Path(output), reference=reference)
    if not isinstance(result, dict):
        msg = f"{CUSTOM_ENTRY_FN} must return a dict, got {type(result).__name__}"
        raise TypeError(msg)
    return result
