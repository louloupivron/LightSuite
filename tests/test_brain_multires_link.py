"""Tests for multires → brain view-registration ROI channel linking."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import yaml

from lightsuite.gui.brain_multires_link import (
    discover_multires_registered_roi_paths,
    embed_crop_in_overview_yxz,
    load_multires_roi_channels_on_registration_grid,
    resample_crop_to_registration_grid,
    resample_overview_volume_to_registration_grid,
    resolve_multires_checkpoint_path,
)
from lightsuite.import_.sample_reference import write_sample_reference
from lightsuite.io.tiff_write import save_registration_volume
from lightsuite.multires.checkpoint import MultiresRegOptsCheckpoint, multires_checkpoint_path
from lightsuite.preprocess.checkpoint import RegOptsCheckpoint
from lightsuite.registration.brain_register import TransformParamsCheckpoint


def _write_brain_checkpoints(save_path: Path, *, shape_yxz: tuple[int, int, int]) -> None:
    save_path.mkdir(parents=True, exist_ok=True)
    chan_path = save_path / "chan_1_sample_register_20um.tif"
    save_registration_volume(np.zeros(shape_yxz, dtype=np.uint16), chan_path)
    RegOptsCheckpoint(
        sample_name="test",
        ny=shape_yxz[0],
        nx=shape_yxz[1],
        nz=shape_yxz[2],
        nchans=1,
        voxel_um=[20.0, 20.0, 20.0],
        registres_um=20.0,
        regvolpath=str(chan_path),
        regvolpath_secondary=None,
        regvolpaths={"1": str(chan_path)},
        tiff_type="channelperfile",
        channel_primary=1,
        channel_secondary=None,
    ).save(save_path / "regopts.json")
    bspline = save_path / "bspline.txt"
    bspline.write_text("(Transform AffineTransform)\n", encoding="utf-8")
    TransformParamsCheckpoint(
        atlas_resolution_um=10.0,
        regvolsize=list(shape_yxz),
        atlassize=list(shape_yxz),
        brain_atlas="allen",
        ori_voxel_um=[20.0, 20.0, 20.0],
        ori_size=list(shape_yxz),
        permute_sample_to_atlas=[1, 2, 3],
        elastix_um_to_mm=0.001,
        tform_bspline_samp20um_to_atlas_20um_px=str(bspline),
        tform_affine_samp20um_to_atlas_10um_px=np.eye(4).tolist(),
        control_point_weight=0.1,
        use_multistep=True,
        use_dual_channel_mi=False,
    ).save(save_path / "transform_params.json")
    write_sample_reference(
        save_path,
        sample_name="test",
        ny=shape_yxz[0],
        nx=shape_yxz[1],
        nz=shape_yxz[2],
        voxel_um=[20.0, 20.0, 20.0],
    )


def test_resolve_multires_checkpoint_from_config(tmp_path: Path) -> None:
    multires_save = tmp_path / "multires_results"
    multires_save.mkdir()
    checkpoint_path = multires_checkpoint_path(multires_save)
    checkpoint_path.write_text("{}", encoding="utf-8")

    from lightsuite.config.models import BrainMultiresLinkConfig

    assert resolve_multires_checkpoint_path(
        BrainMultiresLinkConfig(checkpoint=checkpoint_path),
    ) == checkpoint_path


def test_discover_multires_registered_roi_paths_prefers_full_overview(tmp_path: Path) -> None:
    full = tmp_path / "roi_in_full_overview.tif"
    crop = tmp_path / "roi_crop.tif"
    full.touch()
    crop.touch()
    checkpoint = MultiresRegOptsCheckpoint(
        sample_name="t",
        pair_label="pair",
        pair_manifest_path=str(tmp_path / "pair.json"),
        experiment_slug="exp",
        overview_volume_path=str(tmp_path / "overview.tif"),
        roi_volume_path=str(tmp_path / "roi.tif"),
        registered_roi_path=str(crop),
        registered_roi_full_overview_path=str(full),
        reference_channel="488",
    )
    paths = discover_multires_registered_roi_paths(checkpoint, use_full_overview=True)
    assert paths["488"] == full.resolve()


def test_embed_crop_and_resample_to_registration_grid() -> None:
    overview_shape = (8, 10, 6)
    crop = np.zeros((4, 5, 3), dtype=np.float32)
    crop[1, 2, 1] = 500.0
    embedded = embed_crop_in_overview_yxz(crop, [2, 3, 1], overview_shape)
    assert embedded.shape == overview_shape
    assert float(embedded[4, 4, 2]) == 500.0

    down = resample_overview_volume_to_registration_grid(
        embedded,
        voxel_um=[20.0, 20.0, 20.0],
        registres_um=20.0,
        permute=[1, 2, 3],
        expected_shape=overview_shape,
    )
    assert down.shape == overview_shape
    assert float(down[4, 4, 2]) == 500.0

    crop_only = resample_crop_to_registration_grid(
        crop,
        [2, 3, 1],
        overview_shape,
        voxel_um=[20.0, 20.0, 20.0],
        registres_um=20.0,
        permute=[1, 2, 3],
        expected_shape=overview_shape,
    )
    assert crop_only.shape == overview_shape
    assert float(crop_only[4, 4, 2]) == 500.0


def test_load_multires_roi_channels_on_registration_grid(tmp_path: Path) -> None:
    overview_shape = (8, 10, 6)
    brain_save = tmp_path / "brain_results"
    _write_brain_checkpoints(brain_save, shape_yxz=overview_shape)

    multires_save = tmp_path / "multires_results"
    multires_save.mkdir()
    crop = np.zeros((4, 5, 3), dtype=np.uint16)
    crop[1, 2, 1] = 1234
    crop_path = multires_save / "crop.tif"
    save_registration_volume(crop, crop_path)

    MultiresRegOptsCheckpoint(
        sample_name="t",
        pair_label="pair",
        pair_manifest_path=str(tmp_path / "pair.json"),
        experiment_slug="exp",
        overview_volume_path=str(tmp_path / "overview.tif"),
        roi_volume_path=str(tmp_path / "roi.tif"),
        registered_roi_path=str(crop_path),
        registered_roi_full_overview_path=None,
        reference_channel="488",
        crop_start_index=[2, 3, 1],
        additional_channel_paths={"555": str(crop_path)},
    ).save(multires_checkpoint_path(multires_save))

    brain_yaml = tmp_path / "brain.yaml"
    brain_yaml.write_text(
        yaml.safe_dump(
            {
                "sample": {
                    "name": "test",
                    "source": {"path": str(tmp_path / "data")},
                    "scratch": str(tmp_path / "scratch"),
                    "save_path": str(brain_save),
                    "voxel_um": [20.0, 20.0, 20.0],
                },
                "atlas": {"provider": "allen", "atlas_dir": str(tmp_path / "atlas")},
                "multires": {"checkpoint": str(multires_checkpoint_path(multires_save))},
            }
        ),
        encoding="utf-8",
    )
    (tmp_path / "data").mkdir()

    from lightsuite.config.loader import load_config
    from lightsuite.export.brain_export import _load_transform_params
    from lightsuite.preprocess.checkpoint import RegOptsCheckpoint

    cfg = load_config(brain_yaml)
    checkpoint = RegOptsCheckpoint.load(brain_save / "regopts.json")
    transform_params = _load_transform_params(brain_save)
    channels = load_multires_roi_channels_on_registration_grid(
        cfg,
        checkpoint=checkpoint,
        transform_params=transform_params,
        expected_shape=overview_shape,
    )
    assert "488" in channels
    assert channels["488"].shape == overview_shape
    assert float(channels["488"][4, 4, 2]) == 1234.0
