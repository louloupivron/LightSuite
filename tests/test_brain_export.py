"""Tests for brain export stage."""

from __future__ import annotations

import json
import shutil
from pathlib import Path

import numpy as np
import pytest

from lightsuite.export.brain_export import export_registered_brain_volumes
from lightsuite.export.parcellation import _accumulate_side, write_parcellation_csv
from lightsuite.registration.brain_register import TransformParamsCheckpoint


def test_accumulate_side_median() -> None:
    labels = np.array([0, 1, 1, 2, 2, 2], dtype=np.int32)
    values = np.array([10, 20, 40, 5, 5, 5], dtype=np.float32)
    area_ids = np.array([0, 1, 2], dtype=np.int64)
    med, std, vol = _accumulate_side(labels, values, area_ids, voxel_mm3=0.001)
    assert med[0] == 10
    assert med[1] == 30
    assert med[2] == 5
    assert std[1] == pytest.approx(10.0)
    assert vol[2] == pytest.approx(0.003)


def test_write_parcellation_csv(tmp_path: Path) -> None:
    from lightsuite.export.parcellation import ParcellationResult

    result = ParcellationResult(
        area_ids=np.array([1, 2]),
        median_over_areas=np.array([[1.0, 2.0], [3.0, 4.0]], dtype=np.float32),
        std_over_areas=np.zeros((2, 2), dtype=np.float32),
        volume_over_areas=np.ones((2, 2), dtype=np.float32),
    )
    path = tmp_path / "intensities.csv"
    write_parcellation_csv(path, result)
    text = path.read_text(encoding="utf-8")
    assert "RightSideIntensity" in text
    assert "parcellation_index" in text


def test_export_requires_transformix(tmp_path: Path) -> None:
    if shutil.which("transformix") is not None:
        pytest.skip("transformix is installed; validation-only test")

    from lightsuite.config.loader import load_config
    import yaml

    save = tmp_path / "results"
    save.mkdir()
    (save / "regopts.json").write_text(
        json.dumps(
            {
                "sample_name": "t",
                "ny": 1,
                "nx": 1,
                "nz": 1,
                "nchans": 1,
                "voxel_um": [10, 10, 10],
                "registres_um": 20,
                "regvolpath": str(save / "chan_1_sample_register_20um.tif"),
                "regvolpath_secondary": None,
                "regvolpaths": {"1": str(save / "chan_1_sample_register_20um.tif")},
                "tiff_type": "channelperfile",
                "channel_primary": 1,
                "channel_secondary": None,
            }
        ),
        encoding="utf-8",
    )
    TransformParamsCheckpoint(
        atlas_resolution_um=10,
        regvolsize=[1, 1, 1],
        atlassize=[2, 2, 2],
        brain_atlas="allen",
        ori_voxel_um=[10, 10, 10],
        ori_size=[1, 1, 1],
        permute_sample_to_atlas=[1, 2, 3],
        elastix_um_to_mm=1e-3,
        tform_bspline_samp20um_to_atlas_20um_px=str(save / "bspline.txt"),
        tform_affine_samp20um_to_atlas_10um_px=np.eye(4).tolist(),
        control_point_weight=0.2,
        use_multistep=True,
        use_dual_channel_mi=False,
    ).save(save / "transform_params.json")

    config_path = tmp_path / "config.yaml"
    config_path.write_text(
        yaml.dump(
            {
                "sample": {
                    "name": "t",
                    "source": {"path": str(tmp_path), "tiff_type": "channelperfile"},
                    "scratch": str(tmp_path / "scratch"),
                    "save_path": str(save),
                },
                "detection": {"enabled": False},
            }
        ),
        encoding="utf-8",
    )
    cfg = load_config(config_path)
    with pytest.raises(RuntimeError, match="transformix"):
        export_registered_brain_volumes(cfg)
