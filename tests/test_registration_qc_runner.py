"""Integration tests for registration QC runner (headless)."""

from __future__ import annotations

import json
from pathlib import Path

import nibabel as nib
import numpy as np
import yaml

from lightsuite.analysis.registration_qc_runner import run_registration_qc
from lightsuite.io.tiff_write import save_registration_volume
from lightsuite.registration.brain_register import TransformParamsCheckpoint


def _write_membership(atlas_dir: Path) -> None:
    import pandas as pd

    pd.DataFrame(
        {
            "parcellation_index": [1],
            "parcellation_term_set_name": ["division"],
            "parcellation_term_name": ["Isocortex"],
            "parcellation_term_acronym": ["Isocortex"],
        }
    ).to_csv(atlas_dir / "parcellation_to_parcellation_term_membership.csv", index=False)


def _write_transform_params(save_path: Path, *, shape: tuple[int, int, int]) -> None:
    save_path.mkdir(parents=True, exist_ok=True)
    bspline = save_path / "bspline.txt"
    bspline.write_text("(Transform AffineTransform)\n", encoding="utf-8")
    TransformParamsCheckpoint(
        atlas_resolution_um=10.0,
        regvolsize=list(shape),
        atlassize=list(shape),
        brain_atlas="allen",
        ori_voxel_um=[1.0, 1.0, 1.0],
        ori_size=[20, 20, 20],
        permute_sample_to_atlas=[1, 2, 3],
        elastix_um_to_mm=0.001,
        tform_bspline_samp20um_to_atlas_20um_px=str(bspline),
        tform_affine_samp20um_to_atlas_10um_px=np.eye(4).tolist(),
        control_point_weight=0.1,
        use_multistep=True,
        use_dual_channel_mi=False,
    ).save(save_path / "transform_params.json")


def test_run_registration_qc_headless(tmp_path: Path) -> None:
    atlas_dir = tmp_path / "atlas"
    atlas_dir.mkdir()
    shape = (4, 4, 4)
    annotation = np.zeros(shape, dtype=np.uint16)
    annotation[0, 0, 0] = 1
    nib.save(nib.Nifti1Image(annotation, np.eye(4)), str(atlas_dir / "annotation_10.nii.gz"))
    nib.save(nib.Nifti1Image(np.zeros(shape), np.eye(4)), str(atlas_dir / "average_template_10.nii.gz"))
    _write_membership(atlas_dir)

    (tmp_path / "data").mkdir()
    save_path = tmp_path / "results"
    vr = save_path / "volume_registered"
    vr.mkdir(parents=True)
    _write_transform_params(save_path, shape=shape)

    volume = np.zeros(shape, dtype=np.uint16)
    volume[0, 0, 0] = 200  # unassigned division 0 after map? need check - annotation 0 -> unassigned
    volume[1, 1, 1] = 200  # assigned
    save_registration_volume(volume, vr / "chan_01_registered_atlas.tif")

    config_path = tmp_path / "brain.yaml"
    config_path.write_text(
        yaml.safe_dump(
            {
                "sample": {
                    "name": "test",
                    "source": {"path": str(tmp_path / "data")},
                    "scratch": str(tmp_path / "scratch"),
                    "save_path": str(save_path),
                    "voxel_um": [1.0, 1.0, 1.0],
                },
                "atlas": {"provider": "allen", "atlas_dir": str(atlas_dir)},
            }
        ),
        encoding="utf-8",
    )

    from lightsuite.config.loader import load_config

    cfg = load_config(config_path)
    result = run_registration_qc(
        cfg,
        channel=1,
        threshold=100.0,
        sweep=True,
        headless=True,
    )
    assert result.score_csv is not None
    assert result.score_csv.is_file()
    assert result.sweep_csv is not None
    assert result.sweep_plot is not None
    assert result.score["image_voxels"] >= 1

    json_path = save_path / "registration_qc" / "chan01_unassigned_score.json"
    assert json_path.is_file()
    summary = json.loads(json_path.read_text(encoding="utf-8"))
    assert summary["channel"] == 1
