"""Tests for cooperative spinal stage cancellation."""

from __future__ import annotations

import threading
import time
from pathlib import Path
from unittest.mock import patch

import pytest
import yaml

from lightsuite.config.loader import load_spinal_config
from lightsuite.exceptions import StageCancelledError
from lightsuite.io.cord_reader import read_spinal_cord_sample
from lightsuite.preprocess.cord import preprocess_spinal_cord_sample
from lightsuite.reporter import iter_cancellable_process_map, stage_cancellation


def _ensure_fixtures(root: Path) -> None:
    if not (root / "atlas" / "Template.tif").is_file():
        from tests.fixtures.spinal_cord.build_fixtures import build_fixtures

        build_fixtures(root)


def _write_config(tmp_path: Path, fixture_root: Path) -> Path:
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    save = tmp_path / "results"
    save.mkdir()
    config_data = {
        "sample": {
            "name": "fixture_cord",
            "source": {
                "path": str(fixture_root / "sample"),
                "tiff_type": "channelperfile",
            },
            "scratch": str(scratch),
            "save_path": str(save),
            "voxel_um": [20.0, 20.0, 20.0],
        },
        "atlas": {"atlas_dir": str(fixture_root / "atlas")},
        "registration": {"channel_primary": 1, "resolution_um": 20},
    }
    config_path = tmp_path / "spinal.yaml"
    config_path.write_text(yaml.dump(config_data), encoding="utf-8")
    return config_path


def _slow_identity(value: int) -> int:
    time.sleep(0.05)
    return value


def test_iter_cancellable_process_map_aborts_promptly() -> None:
    cancel = threading.Event()
    seen: list[int] = []
    t0 = time.perf_counter()
    with stage_cancellation(cancel):
        with pytest.raises(StageCancelledError):
            for value in iter_cancellable_process_map(
                _slow_identity,
                list(range(40)),
                max_workers=2,
                chunksize=1,
            ):
                seen.append(value)
                if len(seen) >= 2:
                    cancel.set()
    elapsed = time.perf_counter() - t0
    assert len(seen) < 20
    # Waiting for all 40 jobs would take ~1s+; abort should be much faster.
    assert elapsed < 1.5


def test_read_spinal_cord_sample_respects_cancel_event(tmp_path: Path) -> None:
    fixture_root = Path(__file__).resolve().parent / "fixtures" / "spinal_cord"
    _ensure_fixtures(fixture_root)
    cfg = load_spinal_config(_write_config(tmp_path, fixture_root))

    cancel = threading.Event()
    cancel.set()
    with stage_cancellation(cancel):
        with pytest.raises(StageCancelledError):
            read_spinal_cord_sample(cfg)


def test_preprocess_respects_cancel_event(tmp_path: Path) -> None:
    fixture_root = Path(__file__).resolve().parent / "fixtures" / "spinal_cord"
    _ensure_fixtures(fixture_root)
    cfg = load_spinal_config(_write_config(tmp_path, fixture_root))

    from lightsuite.registration.cord_orientation import save_cord_orientation

    save_cord_orientation(cfg.sample.save_path, "rostrocaudal")

    cancel = threading.Event()

    def _cancel_after_first_check() -> None:
        cancel.set()

    with patch(
        "lightsuite.preprocess.cord.check_stage_cancelled",
        side_effect=_cancel_after_first_check,
    ), stage_cancellation(cancel):
        with pytest.raises(StageCancelledError):
            preprocess_spinal_cord_sample(cfg, headless=True)
