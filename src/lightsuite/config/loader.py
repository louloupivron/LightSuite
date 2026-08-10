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


def _line_indent(line: str) -> int:
    match = re.match(r"^(\s*)", line.rstrip("\n\r"))
    return len(match.group(1)) if match else 0


def _find_section_line(lines: list[str], section: str) -> int | None:
    section_re = re.compile(rf"^{re.escape(section)}:\s*$")
    for index, line in enumerate(lines):
        if section_re.match(line.rstrip("\n\r")):
            return index
    return None


def _section_content_end(lines: list[str], section_idx: int) -> int:
    base_indent = _line_indent(lines[section_idx])
    index = section_idx + 1
    while index < len(lines):
        stripped = lines[index].rstrip("\n\r")
        if stripped == "" or stripped.lstrip().startswith("#"):
            index += 1
            continue
        if _line_indent(lines[index]) <= base_indent:
            break
        index += 1
    return index


def _find_child_section(
    lines: list[str],
    start: int,
    end: int,
    name: str,
    *,
    indent: int | None = None,
) -> int | None:
    pattern = re.compile(rf"^(\s+){re.escape(name)}:\s*$")
    for index in range(start, end):
        match = pattern.match(lines[index].rstrip("\n\r"))
        if match is None:
            continue
        if indent is not None and len(match.group(1)) != indent:
            continue
        return index
    return None


def _set_lateral_flip_in_range(
    lines: list[str],
    start: int,
    end: int,
    lateral_flip: tuple[int, int],
    *,
    default_indent: int,
) -> None:
    flip_value = _format_yaml_list([int(lateral_flip[0]), int(lateral_flip[1])])
    lateral_re = re.compile(r"^(\s*)lateral_flip:\s*.+$")
    for index in range(start, end):
        match = lateral_re.match(lines[index].rstrip("\n\r"))
        if match is None:
            continue
        newline = "\n" if lines[index].endswith("\n") else ""
        lines[index] = f"{match.group(1)}lateral_flip: {flip_value}" + newline
        return
    lines.insert(start, f"{' ' * default_indent}lateral_flip: {flip_value}\n")


def save_mesospim_lateral_flip_to_multires_config(
    config_path: str | Path,
    lateral_flip: tuple[int, int],
) -> Path:
    """Update ``multires.mesospim_geometry`` lateral_flip for overview and ROI."""
    fx, fy = (int(lateral_flip[0]), int(lateral_flip[1]))
    if fx not in (-1, 1) or fy not in (-1, 1):
        msg = "lateral_flip values must be +1 or -1"
        raise ValueError(msg)

    path = Path(config_path).expanduser().resolve()
    if not path.is_file():
        msg = f"Config file not found: {path}"
        raise FileNotFoundError(msg)

    lines = path.read_text(encoding="utf-8").splitlines(keepends=True)
    multires_idx = _find_section_line(lines, "multires")
    if multires_idx is None:
        msg = f"No multires: section in {path}"
        raise ValueError(msg)

    multires_end = _section_content_end(lines, multires_idx)
    multires_indent = _line_indent(lines[multires_idx])
    child_indent = multires_indent + 2
    flip_value = _format_yaml_list([fx, fy])

    geo_idx = _find_child_section(
        lines,
        multires_idx + 1,
        multires_end,
        "mesospim_geometry",
        indent=child_indent,
    )
    if geo_idx is None:
        insert_at = multires_end
        reg_idx = _find_child_section(
            lines,
            multires_idx + 1,
            multires_end,
            "registration",
            indent=child_indent,
        )
        if reg_idx is not None:
            insert_at = reg_idx
        block = [
            f"{' ' * child_indent}mesospim_geometry:\n",
            f"{' ' * (child_indent + 2)}overview:\n",
            f"{' ' * (child_indent + 4)}lateral_flip: {flip_value}\n",
            f"{' ' * (child_indent + 2)}roi:\n",
            f"{' ' * (child_indent + 4)}lateral_flip: {flip_value}\n",
        ]
        lines[insert_at:insert_at] = block
    else:
        geo_indent = _line_indent(lines[geo_idx])
        vol_indent = geo_indent + 2
        flip_indent = geo_indent + 4
        for volume in ("overview", "roi"):
            geo_end = _section_content_end(lines, geo_idx)
            vol_idx = _find_child_section(
                lines,
                geo_idx + 1,
                geo_end,
                volume,
                indent=vol_indent,
            )
            if vol_idx is None:
                lines.insert(
                    geo_end,
                    f"{' ' * vol_indent}{volume}:\n"
                    f"{' ' * flip_indent}lateral_flip: {flip_value}\n",
                )
                continue
            vol_end = _section_content_end(lines, vol_idx)
            _set_lateral_flip_in_range(
                lines,
                vol_idx + 1,
                vol_end,
                (fx, fy),
                default_indent=flip_indent,
            )

    path.write_text("".join(lines), encoding="utf-8")
    load_multires_config(path)
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
