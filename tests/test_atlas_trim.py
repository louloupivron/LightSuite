"""Tests for atlas content trim cache."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import tifffile

from lightsuite.atlas.registry import AtlasPaths
from lightsuite.atlas.trim import trim_atlas_to_cache
from lightsuite.config.models import AtlasConfig, AtlasSource, BrainAtlasId, ContentTrimMode


def test_trim_atlas_auto_caches_cropped_volumes(tmp_path: Path) -> None:
    atlas_dir = tmp_path / "atlas"
    atlas_dir.mkdir()
    shape = (40, 50, 30)
    tpl = np.zeros(shape, dtype=np.float32)
    ann = np.zeros(shape, dtype=np.uint16)
    tpl[8:32, 10:40, 5:25] = 1.0
    ann[8:32, 10:40, 5:25] = 7
    tpl_path = atlas_dir / "template.tif"
    ann_path = atlas_dir / "annotation.tif"
    tifffile.imwrite(tpl_path, tpl)
    tifffile.imwrite(ann_path, ann)

    paths = AtlasPaths(
        brain_atlas="perens",
        atlas_dir=atlas_dir,
        template_path=tpl_path,
        annotation_path=ann_path,
        boundary_path=None,
        structures_csv_path=None,
        supports_parcellation=False,
    )
    cfg = AtlasConfig(
        provider=BrainAtlasId.PERENS,
        resolution_um=20.0,
        atlas_dir=atlas_dir,
        source=AtlasSource.FILES,
        content_trim=ContentTrimMode.AUTO,
        content_margin_vox=0,
    )
    scratch = tmp_path / "scratch"
    resolved = trim_atlas_to_cache(paths, cfg, scratch=scratch)
    assert resolved.is_trimmed
    assert resolved.paths.template_path.is_file()
    from lightsuite.atlas.io import load_atlas_volume

    trimmed = load_atlas_volume(resolved.paths.template_path)
    assert trimmed.shape[0] < shape[0]
    assert resolved.crop_start_yxz[0] == 8
