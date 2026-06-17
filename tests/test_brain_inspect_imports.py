"""Tests for brain import Napari inspect path resolution."""

from __future__ import annotations

import json
from pathlib import Path

import nibabel as nib
import numpy as np
import yaml

from lightsuite.gui.inspect_brain_imports import (
    atlas_points_to_napari_zyx,
    discover_brain_import_inspect_paths,
    load_brain_import_inspect_volumes,
)
from lightsuite.io.tiff_write import save_registration_volume
from lightsuite.registration.brain_register import TransformParamsCheckpoint


def _write_transform_params(save_path: Path, *, shape_yxz: tuple[int, int, int]) -> None:
    save_path.mkdir(parents=True, exist_ok=True)
    bspline = save_path / "bspline.txt"
    bspline.write_text("(Transform AffineTransform)\n", encoding="utf-8")
    checkpoint = TransformParamsCheckpoint(
        atlas_resolution_um=10.0,
        regvolsize=list(shape_yxz),
        atlassize=list(shape_yxz),
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
    )
    checkpoint.save(save_path / "transform_params.json")


def test_discover_brain_import_inspect_paths(tmp_path: Path) -> None:
    atlas_dir = tmp_path / "atlas"
    atlas_dir.mkdir()
    shape = (4, 5, 6)
    template = np.zeros(shape, dtype=np.uint16)
    annotation = np.zeros(shape, dtype=np.uint16)
    nib.save(nib.Nifti1Image(template, np.eye(4)), str(atlas_dir / "average_template_10.nii.gz"))
    nib.save(nib.Nifti1Image(annotation, np.eye(4)), str(atlas_dir / "annotation_10.nii.gz"))

    save_path = tmp_path / "results"
    vr = save_path / "volume_registered"
    vr.mkdir(parents=True)
    _write_transform_params(save_path, shape_yxz=shape)

    chan = np.ones(shape, dtype=np.uint16) * 100
    save_registration_volume(chan, vr / "chan_01_registered_atlas.tif")
    mask = np.zeros(shape, dtype=np.uint8)
    mask[1, 2, 3] = 1
    save_registration_volume(mask.astype(np.uint16), vr / "cells_registered_atlas.tif")
    np.savez_compressed(
        vr / "cells_atlas_coords.npz",
        atlasptcoords=np.array([[2.0, 3.0, 4.0]], dtype=np.float32),
        sampleptcoords=np.array([[10.0, 11.0, 12.0]], dtype=np.float32),
    )

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
    (tmp_path / "data").mkdir()

    from lightsuite.config.loader import load_config

    cfg = load_config(config_path)
    paths = discover_brain_import_inspect_paths(cfg)
    assert paths.registered_channels[1].name == "chan_01_registered_atlas.tif"
    assert "cells" in paths.point_npz_paths
    assert "cells" in paths.mask_paths

    volumes = load_brain_import_inspect_volumes(cfg, paths=paths)
    assert volumes.template.shape == shape
    assert volumes.registered_channels[1].shape == shape
    assert volumes.mask_layers["cells"].shape == shape
    assert volumes.point_layers["cells"].shape == (1, 3)


def test_atlas_points_to_napari_zyx() -> None:
    coords = np.array([[1.0, 2.0, 3.0], [10.0, 20.0, 30.0]])
    napari_pts = atlas_points_to_napari_zyx(coords)
    assert napari_pts[0].tolist() == [2.0, 1.0, 0.0]
    assert napari_pts[1].tolist() == [29.0, 19.0, 9.0]
