"""Tests for spinal cord artifact path layout."""

from __future__ import annotations

from pathlib import Path

from lightsuite.config.models import (
    CordAtlasConfig,
    CordRegistrationConfig,
    CordSampleConfig,
    CordSampleSourceConfig,
    SpinalCordPipelineConfig,
)
from lightsuite.registration.cord_paths import (
    AFFINE_ATLAS_TO_SAMP_FILENAME,
    cord_affine_transform_path,
    cord_bspline_transform_path,
    cord_cache_dir,
    cord_qc_dir,
    cord_scratch_root,
    cord_transforms_dir,
    cord_work_dir,
    resolve_cord_artifact,
)


def _cfg(tmp_path: Path) -> SpinalCordPipelineConfig:
    return SpinalCordPipelineConfig(
        sample=CordSampleConfig(
            name="test_cord",
            source=CordSampleSourceConfig(format="tiff_stack", path=tmp_path),
            scratch=tmp_path / "scratch",
            save_path=tmp_path / "save",
            voxel_um=[20.0, 20.0, 20.0],
        ),
        atlas=CordAtlasConfig(atlas_dir=tmp_path),
        registration=CordRegistrationConfig(),
    )


def test_cord_work_dir_under_sample_scratch(tmp_path: Path) -> None:
    cfg = _cfg(tmp_path)
    work = cord_work_dir(cfg, "transformix", "init", "affine")
    assert work == tmp_path / "scratch" / "test_cord" / "transformix" / "init" / "affine"
    assert work.is_dir()


def test_cord_save_subdirs(tmp_path: Path) -> None:
    cfg = _cfg(tmp_path)
    assert cord_cache_dir(cfg) == tmp_path / "save" / "cache"
    assert cord_transforms_dir(cfg) == tmp_path / "save" / "transforms"
    assert cord_qc_dir(cfg) == tmp_path / "save" / "qc"


def test_resolve_cord_artifact_prefers_nested_layout(tmp_path: Path) -> None:
    save_path = tmp_path / "save"
    nested = save_path / "transforms" / AFFINE_ATLAS_TO_SAMP_FILENAME
    nested.parent.mkdir(parents=True)
    nested.write_text("nested", encoding="utf-8")
    legacy = save_path / AFFINE_ATLAS_TO_SAMP_FILENAME
    legacy.write_text("legacy", encoding="utf-8")
    assert resolve_cord_artifact(save_path, "transforms", AFFINE_ATLAS_TO_SAMP_FILENAME) == nested


def test_resolve_cord_artifact_falls_back_to_legacy(tmp_path: Path) -> None:
    save_path = tmp_path / "save"
    save_path.mkdir()
    legacy = save_path / AFFINE_ATLAS_TO_SAMP_FILENAME
    legacy.write_text("legacy", encoding="utf-8")
    resolved = resolve_cord_artifact(save_path, "transforms", AFFINE_ATLAS_TO_SAMP_FILENAME)
    assert resolved == legacy


def test_cord_transform_helpers(tmp_path: Path) -> None:
    cfg = _cfg(tmp_path)
    save_path = tmp_path / "save"
    nested = save_path / "transforms" / AFFINE_ATLAS_TO_SAMP_FILENAME
    nested.parent.mkdir(parents=True)
    nested.write_text("affine", encoding="utf-8")
    assert cord_affine_transform_path(cfg).read_text(encoding="utf-8") == "affine"
    assert cord_scratch_root(cfg) == tmp_path / "scratch" / "test_cord"
    assert not cord_bspline_transform_path(cfg).is_file()
