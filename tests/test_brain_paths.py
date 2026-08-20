"""Tests for brain pipeline artifact path layout."""

from __future__ import annotations

from pathlib import Path

from lightsuite.registration.brain_paths import (
    AFFINE_FIT_STATS_FILENAME,
    ELASTIX_TEMP,
    QC_DIRNAME,
    REGISTRATION_DIAGNOSTICS_FILENAME,
    WORK_DIRNAME,
    brain_qc_dir,
    brain_qc_file,
    brain_qc_previews_dir,
    brain_work_dir,
    brain_work_root,
    cleanup_brain_work,
    cleanup_legacy_brain_work,
    resolve_brain_artifact,
    resolve_brain_qc_file,
)


def test_brain_work_and_qc_dirs(tmp_path: Path) -> None:
    save = tmp_path / "save"
    work = brain_work_dir(save, ELASTIX_TEMP)
    assert work == save / WORK_DIRNAME / ELASTIX_TEMP
    assert work.is_dir()
    assert brain_work_root(save) == save / WORK_DIRNAME
    assert brain_qc_dir(save) == save / QC_DIRNAME
    assert brain_qc_previews_dir(save) == save / QC_DIRNAME / "previews"
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
    nested = save / "qc" / "previews" / "dim1_initial_registration.png"
    nested.parent.mkdir(parents=True)
    nested.write_text("png", encoding="utf-8")
    assert resolve_brain_artifact(save, "qc", "previews", "dim1_initial_registration.png") == nested


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
