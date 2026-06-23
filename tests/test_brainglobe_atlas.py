"""Tests for BrainGlobe atlas backend integration."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np
import pytest

from lightsuite.analysis.hemisphere import hemisphere_side_masks
from lightsuite.atlas.brainglobe_backend import (
    default_brainglobe_name,
    hemisphere_side_masks_from_brainglobe,
)
from lightsuite.atlas.registry import resolve_brain_atlas, uses_ccf_id_parcellation


def test_default_brainglobe_name_mapping() -> None:
    assert default_brainglobe_name("allen", 10.0) == "allen_mouse_10um"
    assert default_brainglobe_name("perens", 20.0) == "perens_lsfm_mouse_20um"
    assert default_brainglobe_name("perens", 25.0) == "perens_multimodal_lsfm_25um"


def test_uses_ccf_id_parcellation_flags() -> None:
    from lightsuite.atlas.registry import AtlasPaths

    files_perens = AtlasPaths(
        brain_atlas="perens",
        atlas_dir=Path("/tmp"),
        template_path=Path("/tmp/tpl.tif"),
        annotation_path=Path("/tmp/ann.tif"),
        boundary_path=None,
        structures_csv_path=Path("/tmp/s.csv"),
        supports_parcellation=True,
    )
    bg = AtlasPaths(
        brain_atlas="allen",
        atlas_dir=Path("/tmp"),
        template_path=Path("/tmp/ref.tiff"),
        annotation_path=Path("/tmp/ann.tiff"),
        boundary_path=None,
        structures_csv_path=Path("/tmp/structures.csv"),
        supports_parcellation=True,
        atlas_source="brainglobe",
        brainglobe_name="allen_mouse_25um",
    )
    assert uses_ccf_id_parcellation(files_perens) is True
    assert uses_ccf_id_parcellation(bg) is True


def test_hemisphere_side_masks_from_brainglobe_volume() -> None:
    hem = np.zeros((4, 4, 4), dtype=np.int8)
    ann = np.zeros((4, 4, 4), dtype=np.int16)
    hem[:, :, :2] = 2
    hem[:, :, 2:] = 1
    ann[:, :, :] = 100
    side0, side1 = hemisphere_side_masks_from_brainglobe(hem, ann)
    assert side0.sum() == ann[:, :, :2].size
    assert side1.sum() == ann[:, :, 2:].size


@patch("lightsuite.atlas.brainglobe_backend._load_bg_atlas")
def test_resolve_brain_atlas_brainglobe(mock_load: MagicMock, tmp_path: Path) -> None:
    root = tmp_path / "perens_lsfm_mouse_20um_v1.2"
    root.mkdir()
    (root / "reference.tiff").write_bytes(b"")
    (root / "annotation.tiff").write_bytes(b"")
    (root / "structures.csv").write_text("id,name,acronym,structure_id_path,parent_structure_id\n")

    mock_bg = MagicMock()
    mock_bg.root_dir = str(root)
    mock_bg.resolution = (20.0, 20.0, 20.0)
    mock_load.return_value = mock_bg

    resolved = resolve_brain_atlas(
        "perens",
        source="brainglobe",
        brainglobe_name="perens_lsfm_mouse_20um",
    )
    assert resolved.atlas_source == "brainglobe"
    assert resolved.brainglobe_name == "perens_lsfm_mouse_20um"
    assert resolved.supports_parcellation is True
    assert resolved.template_path.name == "reference.tiff"


def test_hemisphere_masks_use_brainglobe_when_name_set() -> None:
    ann = np.ones((2, 2, 2), dtype=np.int16)
    hem = np.array([[[2, 1], [2, 1]], [[2, 1], [2, 1]]], dtype=np.int8)
    with patch(
        "lightsuite.atlas.brainglobe_backend.load_brainglobe_hemispheres",
        return_value=hem,
    ):
        masks = hemisphere_side_masks(ann, "perens", brainglobe_name="perens_lsfm_mouse_20um")
    assert len(masks) == 2
    assert masks[0].sum() + masks[1].sum() == ann.size


def test_brainglobe_missing_dependency_raises() -> None:
    with patch(
        "lightsuite.atlas.brainglobe_backend.require_brainglobe",
        side_effect=ImportError("brainglobe-atlasapi"),
    ):
        with pytest.raises(ImportError, match="brainglobe"):
            resolve_brain_atlas("perens", source="brainglobe", brainglobe_name="x")
