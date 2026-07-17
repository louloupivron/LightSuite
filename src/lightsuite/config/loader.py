"""Load pipeline configuration from YAML files."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any

import yaml

from lightsuite.config.models import BrainPipelineConfig, SpinalCordPipelineConfig
from lightsuite.registration.orientation import validate_permvec


def load_mesospim_config(path: str | Path):
    """Load and validate a mesoSPIM overview / ROI YAML config."""
    from lightsuite.mesospim.config_models import MesospimPipelineConfig

    config_path = Path(path).expanduser().resolve()
    if not config_path.is_file():
        msg = f"Config file not found: {config_path}"
        raise FileNotFoundError(msg)

    with config_path.open(encoding="utf-8") as handle:
        raw: dict[str, Any] = yaml.safe_load(handle) or {}

    return MesospimPipelineConfig.model_validate(raw)


def load_multires_config(path: str | Path):
    """Load and validate a manifest-driven multiresolution YAML config."""
    from lightsuite.multires.config_models import MultiresPipelineConfig

    config_path = Path(path).expanduser().resolve()
    if not config_path.is_file():
        msg = f"Config file not found: {config_path}"
        raise FileNotFoundError(msg)

    with config_path.open(encoding="utf-8") as handle:
        raw: dict[str, Any] = yaml.safe_load(handle) or {}

    return MultiresPipelineConfig.model_validate(raw)


def load_config(path: str | Path) -> BrainPipelineConfig:
    """Load and validate a brain pipeline YAML config."""
    config_path = Path(path).expanduser().resolve()
    if not config_path.is_file():
        msg = f"Config file not found: {config_path}"
        raise FileNotFoundError(msg)

    with config_path.open(encoding="utf-8") as handle:
        raw: dict[str, Any] = yaml.safe_load(handle) or {}

    return BrainPipelineConfig.model_validate(raw)


def load_spinal_config(path: str | Path) -> SpinalCordPipelineConfig:
    """Load and validate a spinal cord pipeline YAML config."""
    config_path = Path(path).expanduser().resolve()
    if not config_path.is_file():
        msg = f"Config file not found: {config_path}"
        raise FileNotFoundError(msg)

    with config_path.open(encoding="utf-8") as handle:
        raw: dict[str, Any] = yaml.safe_load(handle) or {}

    return SpinalCordPipelineConfig.model_validate(raw)


def _format_orientation_line(permvec: list[int], indent: str) -> str:
    inner = ", ".join(str(v) for v in permvec)
    return f"{indent}orientation: [{inner}]"


def save_orientation_to_config(config_path: str | Path, permvec: list[int]) -> Path:
    """Update ``registration.orientation`` in a pipeline YAML config in place."""
    validate_permvec(permvec)
    path = Path(config_path).expanduser().resolve()
    if not path.is_file():
        msg = f"Config file not found: {path}"
        raise FileNotFoundError(msg)

    text = path.read_text(encoding="utf-8")
    lines = text.splitlines(keepends=True)
    orient_re = re.compile(r"^(\s*)orientation:\s*.+$")

    for index, line in enumerate(lines):
        match = orient_re.match(line.rstrip("\n\r"))
        if match is None:
            continue
        newline = "\n" if line.endswith("\n") else ""
        lines[index] = _format_orientation_line(permvec, match.group(1)) + newline
        path.write_text("".join(lines), encoding="utf-8")
        load_config(path)
        return path

    reg_start: int | None = None
    for index, line in enumerate(lines):
        if re.match(r"^registration:\s*$", line.rstrip("\n\r")):
            reg_start = index
            break
    if reg_start is None:
        msg = f"No registration: section in {path}"
        raise ValueError(msg)

    insert_at = reg_start + 1
    indent = "  "
    while insert_at < len(lines):
        stripped = lines[insert_at].rstrip("\n\r")
        if stripped == "" or stripped.lstrip().startswith("#"):
            insert_at += 1
            continue
        if re.match(r"^[^\s#]", stripped):
            break
        key_match = re.match(r"^(\s+)\S", stripped)
        if key_match is None:
            break
        indent = key_match.group(1)
        insert_at += 1

    lines.insert(insert_at, _format_orientation_line(permvec, indent) + "\n")
    path.write_text("".join(lines), encoding="utf-8")
    load_config(path)
    return path


def _format_yaml_list(values: list[int]) -> str:
    inner = ", ".join(str(v) for v in values)
    return f"[{inner}]"


def _update_or_insert_yaml_field(
    lines: list[str],
    *,
    section: str,
    key: str,
    value_line: str,
) -> bool:
    section_re = re.compile(rf"^{re.escape(section)}:\s*$")
    field_re = re.compile(rf"^(\s*){re.escape(key)}:\s*.+$")

    section_start: int | None = None
    section_indent = ""
    for index, line in enumerate(lines):
        if section_re.match(line.rstrip("\n\r")):
            section_start = index
            break
    if section_start is None:
        return False

    insert_at = section_start + 1
    while insert_at < len(lines):
        stripped = lines[insert_at].rstrip("\n\r")
        if stripped == "" or stripped.lstrip().startswith("#"):
            insert_at += 1
            continue
        if re.match(r"^[^\s#]", stripped):
            break
        key_match = re.match(r"^(\s+)\S", stripped)
        if key_match is None:
            break
        section_indent = key_match.group(1)
        if field_re.match(stripped):
            newline = "\n" if lines[insert_at].endswith("\n") else ""
            lines[insert_at] = f"{section_indent}{key}: {value_line}" + newline
            return True
        insert_at += 1

    lines.insert(insert_at, f"{section_indent}{key}: {value_line}\n")
    return True


def save_content_box_to_config(
    config_path: str | Path,
    *,
    target: str,
    box: list[int],
) -> Path:
    """Set manual content crop fields in a pipeline YAML config."""
    from lightsuite.registration.content_bbox import ContentBox
    from lightsuite.registration.content_probe import ContentProbeTarget

    probe_target = ContentProbeTarget(target)
    ContentBox.from_manual_box(box)
    path = Path(config_path).expanduser().resolve()
    if not path.is_file():
        msg = f"Config file not found: {path}"
        raise FileNotFoundError(msg)

    if probe_target == ContentProbeTarget.ATLAS:
        section = "atlas"
        mode_key = "content_trim"
        box_key = "content_box"
    else:
        section = "registration"
        mode_key = "sample_content_crop"
        box_key = "sample_content_box"

    text = path.read_text(encoding="utf-8")
    lines = text.splitlines(keepends=True)
    _update_or_insert_yaml_field(lines, section=section, key=mode_key, value_line="manual")
    _update_or_insert_yaml_field(
        lines,
        section=section,
        key=box_key,
        value_line=_format_yaml_list(box),
    )
    path.write_text("".join(lines), encoding="utf-8")
    load_config(path)
    return path


def load_cohort_config(path: str | Path):
    """Load and validate a cross-subject cohort YAML config."""
    from lightsuite.analysis.cohort_models import CohortConfig

    config_path = Path(path).expanduser().resolve()
    if not config_path.is_file():
        msg = f"Cohort config file not found: {config_path}"
        raise FileNotFoundError(msg)

    with config_path.open(encoding="utf-8") as handle:
        raw: dict[str, Any] = yaml.safe_load(handle) or {}

    return CohortConfig.model_validate(raw)
