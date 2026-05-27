"""Tests for atlas path resolution."""

from __future__ import annotations

from pathlib import Path

import pytest

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
