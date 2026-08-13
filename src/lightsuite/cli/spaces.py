"""CLI helpers for registration output space selection."""

from __future__ import annotations

import typer

_VALID_SPACES = frozenset({"atlas", "sample"})


def export_spaces_from_checks(*, atlas: bool, sample: bool) -> list[str]:
    """Build an export space list from GUI checkboxes."""
    spaces: list[str] = []
    if atlas:
        spaces.append("atlas")
    if sample:
        spaces.append("sample")
    if not spaces:
        msg = "Select at least one export space: atlas and/or sample."
        raise ValueError(msg)
    return spaces


def default_export_space_checks(spaces: list[str] | None) -> tuple[bool, bool]:
    """Return (atlas, sample) checkbox defaults from YAML ``export.spaces``."""
    raw = spaces if spaces is not None else ["atlas"]
    normalized = {str(item).strip().lower() for item in raw}
    return "atlas" in normalized, "sample" in normalized


def format_export_spaces(spaces: list[str]) -> str:
    """Human-readable label for log messages."""
    normalized = [str(item).strip().lower() for item in spaces]
    if normalized == ["atlas", "sample"] or normalized == ["sample", "atlas"]:
        return "atlas + sample"
    return ", ".join(normalized)


def parse_spaces_option(space: str | None) -> list[str] | None:
    """Parse ``--space atlas|sample|both`` into a list or None (use YAML defaults)."""
    if space is None:
        return None
    key = space.strip().lower()
    if key == "both":
        return ["atlas", "sample"]
    if key in _VALID_SPACES:
        return [key]
    msg = f"Invalid --space {space!r}; use atlas, sample, or both."
    raise typer.BadParameter(msg)


def parse_view_space_option(space: str = "atlas") -> str:
    """Parse a single view/QC space (``atlas`` or ``sample``)."""
    key = space.strip().lower()
    if key in _VALID_SPACES:
        return key
    msg = f"Invalid --space {space!r}; use atlas or sample."
    raise typer.BadParameter(msg)
