"""Parse mesoSPIM TIFF sidecar metadata."""

from __future__ import annotations

import re
from pathlib import Path

_META_LINE = re.compile(r"^\[([^\]]+)\]\s*(.*)$")


def meta_path_for_tiff(tiff_path: Path) -> Path:
    """Return the default mesoSPIM meta sidecar path for a stack TIFF."""
    tiff_path = tiff_path.expanduser()
    return tiff_path.parent / f"{tiff_path.name}_meta.txt"


def parse_mesospim_meta(path: Path) -> dict[str, float | int | str]:
    """Parse mesoSPIM export meta file into a flat dict (numeric values coerced)."""
    out: dict[str, float | int | str] = {}
    text = path.expanduser().read_text(encoding="utf-8", errors="replace")
    for raw in text.splitlines():
        line = raw.strip()
        match = _META_LINE.match(line)
        if not match:
            continue
        key, val = match.group(1).strip(), match.group(2).strip()
        if not val:
            continue
        try:
            if "." in val or "e" in val.lower():
                out[key] = float(val)
            else:
                out[key] = int(val)
        except ValueError:
            out[key] = val
    return out
