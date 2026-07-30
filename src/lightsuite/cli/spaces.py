"""CLI helpers for registration output space selection."""

from __future__ import annotations

import typer

_VALID_SPACES = frozenset({"atlas", "sample"})


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
