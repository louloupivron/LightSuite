"""Load and save multiresolution pair manifests."""

from __future__ import annotations

import json
from pathlib import Path

from lightsuite.multires.models import MANIFEST_FORMAT, MultiresPairManifest

__all__ = [
    "MANIFEST_FORMAT",
    "MultiresPairManifest",
    "load_pair_manifest",
    "save_pair_manifest",
]


def load_pair_manifest(path: str | Path) -> MultiresPairManifest:
    """Load and validate a pair manifest JSON file."""
    path = Path(path).expanduser().resolve()
    if not path.is_file():
        msg = f"Pair manifest not found: {path}"
        raise FileNotFoundError(msg)
    raw = json.loads(path.read_text(encoding="utf-8"))
    manifest = MultiresPairManifest.from_dict(raw)
    _validate_manifest_paths(manifest, manifest_path=path)
    return manifest


def save_pair_manifest(manifest: MultiresPairManifest, path: str | Path) -> Path:
    """Write a pair manifest JSON file."""
    path = Path(path).expanduser()
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(manifest.to_dict(), indent=2), encoding="utf-8")
    return path


def _validate_manifest_paths(manifest: MultiresPairManifest, *, manifest_path: Path) -> None:
    base = manifest_path.parent
    for label, spec in (("overview", manifest.overview), ("roi", manifest.roi)):
        volume_path = Path(spec.volume_path).expanduser()
        if not volume_path.is_absolute():
            volume_path = (base / volume_path).resolve()
        if not volume_path.exists():
            msg = f"Manifest {label} volume_path does not exist: {volume_path}"
            raise FileNotFoundError(msg)
