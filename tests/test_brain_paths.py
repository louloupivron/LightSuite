"""Tests for brain pipeline artifact path layout."""

from __future__ import annotations

from pathlib import Path

from lightsuite.registration.brain_paths import (
    AFFINE_FIT_STATS_FILENAME,
    ELASTIX_TEMP,
    IMPORTS_DIRNAME,
    INIT_DIAGNOSTICS_FILENAME,
    QC_DIRNAME,
    REGISTRATION_DIAGNOSTICS_FILENAME,
    STATS_DIRNAME,
    WORK_DIRNAME,
    brain_imports_dir,
    brain_qc_dir,
    brain_qc_file,
    brain_qc_previews_dir,
    brain_stats_dir,
    brain_work_dir,
    brain_work_root,
    cleanup_brain_qc_audit_json,
    cleanup_brain_work,
    cleanup_legacy_brain_work,
    iter_brain_import_paths,
    iter_brain_stats_paths,
    resolve_brain_artifact,
    resolve_brain_imports_file,
    resolve_brain_qc_file,
    resolve_brain_stats_file,
)


def test_brain_work_and_qc_dirs(tmp_path: Path) -> None:
    save = tmp_path / "save"
    work = brain_work_dir(save, ELASTIX_TEMP)
    assert work == save / WORK_DIRNAME / ELASTIX_TEMP
    assert work.is_dir()
    assert brain_work_root(save) == save / WORK_DIRNAME
    assert brain_qc_dir(save) == save / QC_DIRNAME
    assert brain_qc_previews_dir(save) == save / QC_DIRNAME
    assert brain_qc_file(save, AFFINE_FIT_STATS_FILENAME) == (
        save / QC_DIRNAME / AFFINE_FIT_STATS_FILENAME
    )


def test_resolve_brain_qc_prefers_nested(tmp_path: Path) -> None:
    save = tmp_path / "save"
    nested = brain_qc_file(save, REGISTRATION_DIAGNOSTICS_FILENAME)
    nested.write_text("nested", encoding="utf-8")
    legacy = save / REGISTRATION_DIAGNOSTICS_FILENAME
    legacy.write_text("legacy", encoding="utf-8")
    assert resolve_brain_qc_file(save, REGISTRATION_DIAGNOSTICS_FILENAME) == nested


def test_resolve_brain_qc_falls_back_to_legacy(tmp_path: Path) -> None:
    save = tmp_path / "save"
    save.mkdir()
    legacy = save / REGISTRATION_DIAGNOSTICS_FILENAME
    legacy.write_text("legacy", encoding="utf-8")
    assert resolve_brain_qc_file(save, REGISTRATION_DIAGNOSTICS_FILENAME) == legacy


def test_resolve_brain_artifact_generic(tmp_path: Path) -> None:
    save = tmp_path / "save"
    nested = save / "qc" / "dim1_initial_registration.png"
    nested.parent.mkdir(parents=True)
    nested.write_text("png", encoding="utf-8")
    assert resolve_brain_artifact(save, "qc", "dim1_initial_registration.png") == nested


def test_resolve_brain_qc_preview_legacy_subfolder(tmp_path: Path) -> None:
    from lightsuite.registration.brain_paths import resolve_brain_qc_preview

    save = tmp_path / "save"
    legacy = save / "qc" / "previews" / "dim1_initial_registration.png"
    legacy.parent.mkdir(parents=True)
    legacy.write_text("png", encoding="utf-8")
    assert resolve_brain_qc_preview(save, "dim1_initial_registration.png") == legacy


def test_cleanup_brain_work_and_legacy(tmp_path: Path) -> None:
    save = tmp_path / "save"
    work = brain_work_dir(save, ELASTIX_TEMP)
    (work / "junk.raw").write_bytes(b"x" * 10)
    legacy = save / "elastix_temp"
    legacy.mkdir(parents=True)
    (legacy / "junk.raw").write_bytes(b"y" * 10)
    cleanup_brain_work(save, ELASTIX_TEMP)
    cleanup_legacy_brain_work(save)
    assert not work.exists()
    assert not legacy.exists()
    assert (save / WORK_DIRNAME).is_dir() or not (save / WORK_DIRNAME).exists()


def test_cleanup_brain_qc_audit_json(tmp_path: Path) -> None:
    save = tmp_path / "save"
    previews = brain_qc_previews_dir(save)
    preview_png = previews / "dim1_initial_registration.png"
    preview_png.write_text("png", encoding="utf-8")
    nested = brain_qc_file(save, REGISTRATION_DIAGNOSTICS_FILENAME)
    nested.write_text("{}", encoding="utf-8")
    legacy = save / INIT_DIAGNOSTICS_FILENAME
    legacy.write_text("{}", encoding="utf-8")
    cleanup_brain_qc_audit_json(save)
    assert not nested.exists()
    assert not legacy.exists()
    assert preview_png.is_file()
    assert brain_qc_dir(save).is_dir()


def test_brain_stats_and_imports_dirs(tmp_path: Path) -> None:
    save = tmp_path / "save"
    assert brain_stats_dir(save) == save / STATS_DIRNAME
    assert brain_imports_dir(save) == save / IMPORTS_DIRNAME
    assert brain_stats_dir(save).is_dir()
    assert brain_imports_dir(save).is_dir()


def test_resolve_brain_stats_and_imports_with_legacy(tmp_path: Path) -> None:
    save = tmp_path / "save"
    nested_stats = brain_stats_dir(save) / "region_stats.csv"
    nested_stats.write_text("nested", encoding="utf-8")
    legacy_stats = save / "volume_registered" / "chan01_region_stats.csv"
    legacy_stats.parent.mkdir(parents=True)
    legacy_stats.write_text("legacy", encoding="utf-8")
    assert resolve_brain_stats_file(save, "region_stats.csv") == nested_stats
    assert resolve_brain_stats_file(save, "chan01_region_stats.csv") == legacy_stats

    nested_import = brain_imports_dir(save) / "cells_atlas_coords.npz"
    nested_import.write_bytes(b"not npz")
    assert resolve_brain_imports_file(save, "cells_atlas_coords.npz") == nested_import


def test_iter_brain_import_and_stats_paths(tmp_path: Path) -> None:
    save = tmp_path / "save"
    imports = brain_imports_dir(save)
    (imports / "a_atlas_coords.npz").write_bytes(b"x")
    vr = save / "volume_registered"
    vr.mkdir(parents=True)
    (vr / "b_atlas_coords.npz").write_bytes(b"y")
    stats = brain_stats_dir(save)
    (stats / "region_stats.csv").write_text("stats", encoding="utf-8")
    (vr / "region_stats_sample.csv").write_text("legacy", encoding="utf-8")

    import_paths = [p.name for p in iter_brain_import_paths(save, "*_atlas_coords.npz")]
    assert import_paths == ["a_atlas_coords.npz", "b_atlas_coords.npz"]

    stats_paths = [p.name for p in iter_brain_stats_paths(save, "region_stats*.csv")]
    assert stats_paths == ["region_stats.csv", "region_stats_sample.csv"]
