"""Tests for convert-annotations prepare stage."""

from __future__ import annotations

import json
from pathlib import Path

import pytest
import yaml

from lightsuite.config.loader import load_config
from lightsuite.import_.convert import run_convert_annotations
from lightsuite.import_.sample_reference import SampleReference
from lightsuite.import_.validate import validate_points_csv


def _write_brain_config(path: Path, save_path: Path, *, import_block: dict) -> None:
    data = {
        "sample": {
            "name": "test_mouse",
            "source": {
                "format": "tiff_stack",
                "path": str(save_path / "source"),
                "tiff_type": "channelperfile",
            },
            "scratch": str(save_path / "scratch"),
            "save_path": str(save_path / "results"),
            "voxel_um": [5.0, 5.0, 5.0],
        },
        "atlas": {
            "provider": "allen",
            "resolution_um": 10,
            "atlas_dir": str(save_path / "atlas"),
        },
        "registration": {"resolution_um": 20, "channel_primary": 1},
        "detection": {"enabled": False},
        "compute": {"workers": 1},
        "export": {"registered_volume_format": "tiff"},
        "import": import_block,
    }
    path.write_text(yaml.safe_dump(data), encoding="utf-8")
    (save_path / "source").mkdir(parents=True, exist_ok=True)
    (save_path / "scratch").mkdir(parents=True, exist_ok=True)
    (save_path / "results").mkdir(parents=True, exist_ok=True)
    (save_path / "atlas").mkdir(parents=True, exist_ok=True)


def _write_reference(save_path: Path) -> SampleReference:
    ref = SampleReference.from_checkpoint(
        sample_name="test_mouse",
        ny=100,
        nx=100,
        nz=50,
        voxel_um=[5.0, 5.0, 5.0],
    )
    ref.save(save_path / "results" / "sample_reference.json")
    return ref


def test_convert_smartspim_and_validate(tmp_path: Path) -> None:
    cfg_path = tmp_path / "brain.yaml"
    results = tmp_path / "results"
    source_json = tmp_path / "points.json"
    source_json.write_text(json.dumps([[10, 20, 30], [11, 21, 31]]), encoding="utf-8")
    out_csv = tmp_path / "out_points.csv"
    _write_brain_config(
        cfg_path,
        tmp_path,
        import_block={
            "converter": {
                "suite": "smartspim",
                "source": str(source_json),
                "output": str(out_csv),
                "label": "cells",
            }
        },
    )
    _write_reference(tmp_path)
    cfg = load_config(cfg_path)
    result = run_convert_annotations(
        import_config=cfg.import_config,
        save_path=results,
    )
    assert result.n_converted == 2
    assert out_csv.is_file()
    assert result.summary_path is not None and result.summary_path.is_file()
    assert all(v.ok for v in result.validations)


def test_custom_converter_requires_validation(tmp_path: Path) -> None:
    cfg_path = tmp_path / "brain.yaml"
    results = tmp_path / "results"
    entry = Path("examples/annotation_sample/custom_converter_example.py").resolve()
    source = tmp_path / "native.csv"
    source.write_text("x,y,z\n10,20,30\n", encoding="utf-8")
    out_csv = tmp_path / "converted.csv"
    _write_brain_config(
        cfg_path,
        tmp_path,
        import_block={
            "converter": {
                "suite": "custom",
                "source": str(source),
                "output": str(out_csv),
                "custom_entry": str(entry),
                "label": "custom",
            }
        },
    )
    _write_reference(tmp_path)
    cfg = load_config(cfg_path)
    result = run_convert_annotations(import_config=cfg.import_config, save_path=results)
    assert result.n_converted == 1
    assert all(v.ok for v in result.validations)


def test_validate_points_rejects_out_of_bounds(tmp_path: Path) -> None:
    ref = SampleReference.from_checkpoint(
        sample_name="t", ny=10, nx=10, nz=10, voxel_um=[1, 1, 1]
    )
    csv_path = tmp_path / "bad.csv"
    csv_path.write_text("x,y,z\n1000,1000,1000\n", encoding="utf-8")
    result = validate_points_csv(csv_path, ref, label="bad")
    assert result.ok is False


def test_converter_custom_missing_entry_fails_load(tmp_path: Path) -> None:
    cfg_path = tmp_path / "brain.yaml"
    with pytest.raises(Exception):
        _write_brain_config(
            cfg_path,
            tmp_path,
            import_block={"converter": {"suite": "custom", "source": str(tmp_path / "x.json")}},
        )
        load_config(cfg_path)
