"""LightSuite-specific errors with actionable hints."""

from __future__ import annotations

from pydantic import ValidationError


class LightsuiteError(Exception):
    """Base class for user-facing LightSuite errors."""


class LightsuiteConfigError(LightsuiteError):
    """Invalid or incomplete pipeline YAML."""


class LightsuiteCheckpointError(LightsuiteError):
    """Missing or stale checkpoint between pipeline stages."""


class LightsuiteExternalToolError(LightsuiteError):
    """Required external binary or library is missing."""


class StageCancelledError(LightsuiteError):
    """Pipeline stage stopped by the user."""


def format_validation_error(exc: ValidationError, *, config_path: str | None = None) -> str:
    """Turn Pydantic validation errors into a short, actionable message."""
    lines: list[str] = []
    if config_path:
        lines.append(f"Config validation failed: {config_path}")
    else:
        lines.append("Config validation failed.")
    for error in exc.errors():
        loc = ".".join(str(part) for part in error.get("loc", ()))
        msg = error.get("msg", "invalid value")
        lines.append(f"  • {loc}: {msg}")
    lines.append("Edit the YAML paths/fields above, or run: lightsuite config explain -c <config>")
    return "\n".join(lines)
