"""Load pipeline configuration from YAML files."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml

from lightsuite.config.models import BrainPipelineConfig


def load_config(path: str | Path) -> BrainPipelineConfig:
    """Load and validate a brain pipeline YAML config."""
    config_path = Path(path).expanduser().resolve()
    if not config_path.is_file():
        msg = f"Config file not found: {config_path}"
        raise FileNotFoundError(msg)

    with config_path.open(encoding="utf-8") as handle:
        raw: dict[str, Any] = yaml.safe_load(handle) or {}

    return BrainPipelineConfig.model_validate(raw)
