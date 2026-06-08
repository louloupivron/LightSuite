"""Tests for AP-constrained auto control-point refinement."""

from __future__ import annotations

from pathlib import Path

import nibabel as nib
import numpy as np
import tifffile
import yaml

from lightsuite.config.loader import load_config
from lightsuite.gui.brain_data import prepare_brain_align_slices_session
from lightsuite.gui.slice_correspondence import SliceAnchor, SliceCorrespondence
from lightsuite.preprocess.brain import preprocess_lightsheet_volume
from lightsuite.preprocess.checkpoint import RegOptsCheckpoint
from lightsuite.registration.ap_correspondence import (
    ap_residuals_vox,
    filter_pairs_by_ap_correspondence,
)
from lightsuite.registration.init_brain import initialize_brain_registration
from lightsuite.registration.refine_auto_points import refine_brain_auto_points


def _correspondence_x_axis() -> SliceCorrespondence:
    return SliceCorrespondence.single_axis(
        2,
        np.eye(4).tolist(),
        [
            SliceAnchor(sample_index=10, atlas_plane=20, confirmed=True),
            SliceAnchor(sample_index=30, atlas_plane=50, confirmed=True),
        ],
    )


def test_ap_filter_keeps_consistent_pairs() -> None:
    corr = _correspondence_x_axis()
    atlas_shape = (40, 60, 30)
    # sample x=19 (1-based 20) expects atlas plane ~35 along X
    cpsample = np.array([[19.0, 5.0, 5.0], [19.0, 8.0, 8.0]], dtype=float)
    cpatlas = np.array([[34.0, 5.0, 5.0], [5.0, 8.0, 8.0]], dtype=float)
    filtered_s, filtered_a, stats = filter_pairs_by_ap_correspondence(
        cpsample,
        cpatlas,
        corr,
        np.eye(4),
        atlas_shape,
        tolerance_vox=12.0,
        min_pairs_kept=1,
    )
    assert filtered_s.shape[0] == 1
    assert np.allclose(filtered_s[0], cpsample[0])
    assert np.allclose(filtered_a[0], cpatlas[0])
    assert stats.pairs_removed == 1


def test_ap_filter_relaxes_tolerance_when_needed() -> None:
    corr = _correspondence_x_axis()
    atlas_shape = (40, 60, 30)
    cpsample = np.array([[19.0, 5.0, 5.0], [19.0, 6.0, 6.0]], dtype=float)
    cpatlas = np.array([[34.0, 5.0, 5.0], [30.0, 6.0, 6.0]], dtype=float)
    _, _, stats = filter_pairs_by_ap_correspondence(
        cpsample,
        cpatlas,
        corr,
        np.eye(4),
        atlas_shape,
        tolerance_vox=5.0,
        min_pairs_kept=2,
    )
    assert stats.pairs_after == 2
    assert stats.tolerance_vox == 5.0


def test_multi_axis_filter_uses_max_residual() -> None:
    atlas_shape = (40, 60, 30)
    corr = SliceCorrespondence(
        original_trans=np.eye(4).tolist(),
        axes={
            2: [
                SliceAnchor(sample_index=10, atlas_plane=20, confirmed=True),
                SliceAnchor(sample_index=30, atlas_plane=50, confirmed=True),
            ],
            3: [
                SliceAnchor(sample_index=10, atlas_plane=20, confirmed=True),
                SliceAnchor(sample_index=30, atlas_plane=50, confirmed=True),
            ],
        },
    )
    cpsample = np.array([[19.0, 5.0, 19.0], [19.0, 5.0, 5.0]], dtype=float)
    cpatlas = np.array([[34.0, 5.0, 19.0], [34.0, 5.0, 5.0]], dtype=float)
    _, _, stats = filter_pairs_by_ap_correspondence(
        cpsample,
        cpatlas,
        corr,
        np.eye(4),
        atlas_shape,
        tolerance_vox=12.0,
        min_pairs_kept=1,
    )
    assert stats.active_axes == [2, 3]
    assert stats.pairs_after == 1


def test_ap_residuals_identity_transform() -> None:
    corr = _correspondence_x_axis()
    cpsample = np.array([[19.0, 0.0, 0.0]], dtype=float)
    cpatlas = np.array([[34.0, 0.0, 0.0]], dtype=float)
    residuals = ap_residuals_vox(cpsample, cpatlas, corr, np.eye(4), (40, 60, 30))
    assert residuals.shape == (1,)
    assert residuals[0] == 0.0


def _write_channel_stack(path: Path, slices: list[np.ndarray]) -> None:
    tifffile.imwrite(path, np.stack(slices, axis=0), photometric="minisblack")


def test_refine_auto_points_integration(tmp_path: Path) -> None:
    data_dir = tmp_path / "data"
    data_dir.mkdir()
    scratch = tmp_path / "scratch"
    save = tmp_path / "results"
    save.mkdir()
    atlas_dir = tmp_path / "atlas"
    atlas_dir.mkdir()
    shape = (20, 20, 20)
    yy, xx, zz = np.mgrid[0 : shape[0], 0 : shape[1], 0 : shape[2]]
    template = (255 * np.exp(-((yy - 10) ** 2 + (xx - 10) ** 2 + (zz - 10) ** 2) / 30)).astype(
        np.float32
    )
    ann = np.zeros(shape, dtype=np.float32)
    ann[5:15, 5:15, 5:15] = 1
    nib.save(nib.Nifti1Image(template, np.eye(4)), atlas_dir / "average_template_10.nii.gz")
    nib.save(nib.Nifti1Image(ann, np.eye(4)), atlas_dir / "annotation_10.nii.gz")

    yy2, xx2 = np.mgrid[0:32, 0:32]
    ring = (((yy2 - 16) ** 2 + (xx2 - 16) ** 2) > 36) & (((yy2 - 16) ** 2 + (xx2 - 16) ** 2) < 100)
    base = (150 + 100 * ring).astype(np.uint16)
    _write_channel_stack(data_dir / "ch1.tif", [base + z * 5 for z in range(6)])

    config_data = {
        "sample": {
            "name": "test",
            "source": {"path": str(data_dir), "tiff_type": "channelperfile"},
            "scratch": str(scratch),
            "save_path": str(save),
            "voxel_um": [10.0, 10.0, 10.0],
        },
        "atlas": {"provider": "allen", "resolution_um": 10, "atlas_dir": str(atlas_dir)},
        "registration": {"resolution_um": 20, "channel_primary": 1, "orientation": [1, 2, 3]},
        "detection": {"enabled": False},
    }
    config_path = tmp_path / "config.yaml"
    config_path.write_text(yaml.dump(config_data), encoding="utf-8")
    cfg = load_config(config_path)
    preprocess_lightsheet_volume(cfg)
    initialize_brain_registration(cfg)

    checkpoint = RegOptsCheckpoint.load(save / "regopts.json")
    n_before = len(checkpoint.autocpsample or [])
    assert n_before > 0

    prepare_brain_align_slices_session(cfg)
    refine_brain_auto_points(cfg, force=True)

    checkpoint = RegOptsCheckpoint.load(save / "regopts.json")
    assert checkpoint.auto_points_refined is True
    assert checkpoint.auto_points_mode == "multi_axis_filter"
    assert len(checkpoint.autocpsample or []) <= n_before
    assert (save / "auto_points_refine_stats.json").is_file()

    refine_brain_auto_points(cfg, force=False)
    checkpoint2 = RegOptsCheckpoint.load(save / "regopts.json")
    assert len(checkpoint2.autocpsample or []) == len(checkpoint.autocpsample or [])
