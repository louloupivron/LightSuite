"""Tests for Napari division viewer (headless path resolution)."""

from __future__ import annotations

import json
from pathlib import Path

import nibabel as nib
import numpy as np
import yaml

from lightsuite.gui.view_divisions_brain import (
    discover_division_viewer_paths,
    load_division_viewer_volumes,
    run_brain_division_viewer,
)
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


def test_division_viewer_headless(tmp_path: Path) -> None:
    atlas_dir = tmp_path / "atlas"
    atlas_dir.mkdir()
    shape = (4, 5, 6)
    annotation = np.zeros(shape, dtype=np.uint16)
    annotation[0, 0, 0] = 1
    nib.save(nib.Nifti1Image(annotation, np.eye(4)), str(atlas_dir / "annotation_10.nii.gz"))
    nib.save(nib.Nifti1Image(np.zeros(shape), np.eye(4)), str(atlas_dir / "average_template_10.nii.gz"))
    _write_membership(atlas_dir)

    save_path = tmp_path / "results"
    (tmp_path / "data").mkdir()
    vr = save_path / "volume_registered"
    vr.mkdir(parents=True)
    _write_transform_params(save_path, shape=shape)
    chan = np.ones(shape, dtype=np.uint16) * 50
    save_registration_volume(chan, vr / "chan_01_registered_atlas.tif")

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
    paths = run_brain_division_viewer(cfg, headless=True)
    assert paths.division_labels.is_file()
    assert paths.division_legend.is_file()
    assert "channel 1" in paths.channel_paths

    channels, labels, contrast = load_division_viewer_volumes(cfg, paths=paths, stride=2)
    assert len(channels) == 1
    assert labels.shape == (2, 3, 3)
    assert "channel 1" in contrast

    discovered = discover_division_viewer_paths(cfg)
    assert discovered.division_labels == paths.division_labels
