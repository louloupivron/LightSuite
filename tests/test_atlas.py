"""Tests for atlas path resolution and volume I/O."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
import tifffile

from lightsuite.atlas.io import load_atlas_volume
from lightsuite.atlas.registry import resolve_brain_atlas


def test_resolve_allen_atlas_explicit_dir(tmp_path: Path) -> None:
    atlas_dir = tmp_path / "allen"
    atlas_dir.mkdir()
    (atlas_dir / "average_template_10.nii.gz").write_bytes(b"")
    (atlas_dir / "annotation_10.nii.gz").write_bytes(b"")

    resolved = resolve_brain_atlas("allen", atlas_dir=atlas_dir)
    assert resolved.brain_atlas == "allen"
    assert resolved.template_path.name == "average_template_10.nii.gz"
    assert resolved.boundary_path is None
    assert resolved.supports_parcellation is True


def test_resolve_allen_atlas_with_boundary(tmp_path: Path) -> None:
    atlas_dir = tmp_path / "allen"
    atlas_dir.mkdir()
    (atlas_dir / "average_template_10.nii.gz").write_bytes(b"")
    (atlas_dir / "annotation_10.nii.gz").write_bytes(b"")
    (atlas_dir / "annotation_boundary_10.nii.gz").write_bytes(b"")

    resolved = resolve_brain_atlas("allen", atlas_dir=atlas_dir)
    assert resolved.boundary_path is not None
    assert resolved.boundary_path.name == "annotation_boundary_10.nii.gz"


def test_resolve_atlas_missing_raises(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError):
        resolve_brain_atlas("allen", atlas_dir=tmp_path)


def test_load_atlas_volume_tiff(tmp_path: Path) -> None:
    vol = np.arange(24, dtype=np.uint16).reshape(2, 3, 4)
    path = tmp_path / "reference.tiff"
    tifffile.imwrite(path, vol)
    loaded = load_atlas_volume(path)
    np.testing.assert_array_equal(loaded, vol)

