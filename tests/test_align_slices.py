"""Tests for AP slice correspondence (align-slices stage)."""

from __future__ import annotations

from pathlib import Path

import nibabel as nib
import numpy as np
import tifffile
import yaml

from lightsuite.config.loader import load_config
from lightsuite.gui.brain_data import (
    apply_slice_correspondence_to_session,
    load_brain_match_points_data,
)
from lightsuite.gui.chooselist import (
    default_ap_cut_axis,
    generate_ap_alignment_list,
)
from lightsuite.gui.control_points import ControlPointSession
from lightsuite.gui.slice_correspondence import SliceAnchor, SliceCorrespondence
from lightsuite.preprocess.brain import preprocess_lightsheet_volume
from lightsuite.registration.init_brain import initialize_brain_registration
from lightsuite.gui.brain_data import prepare_brain_align_slices_session


def test_default_ap_cut_axis_is_longest() -> None:
    assert default_ap_cut_axis((80, 120, 60)) == 2


def test_generate_ap_alignment_list_shape() -> None:
    chooselist = generate_ap_alignment_list((100, 120, 80), cut_axis=2, n_slices=15)
    assert chooselist.shape == (15, 4)
    assert np.all(chooselist[:, 1] == 2)


def test_slice_correspondence_interpolation() -> None:
    corr = SliceCorrespondence(
        cut_axis=2,
        original_trans=np.eye(4).tolist(),
        anchors=[
            SliceAnchor(sample_index=10, atlas_plane=20, confirmed=True),
            SliceAnchor(sample_index=30, atlas_plane=50, confirmed=True),
        ],
    )
    assert corr.interpolate_atlas_plane(10, 2, 100) == 20
    assert corr.interpolate_atlas_plane(20, 2, 100) == 35
    assert corr.interpolate_atlas_plane(40, 2, 100) == 50
    assert corr.interpolate_atlas_plane(20, 3, 100) is None


def test_apply_slice_correspondence_to_session() -> None:
    chooselist = generate_ap_alignment_list((40, 60, 30), cut_axis=2, n_slices=5)
    session = ControlPointSession.empty(np.eye(4), chooselist.shape[0])
    corr = SliceCorrespondence(
        cut_axis=2,
        original_trans=np.eye(4).tolist(),
        anchors=[
            SliceAnchor(sample_index=int(row[0]), atlas_plane=10 + i * 5, confirmed=True)
            for i, row in enumerate(chooselist)
        ],
    )
    apply_slice_correspondence_to_session(session, chooselist, corr, (40, 60, 30))
    assert session.atlas_slice_indices is not None
    assert session.atlas_slice_indices[0] == 10
    assert session.atlas_slice_indices[-1] == 10 + 4 * 5


def _write_channel_stack(path: Path, slices: list[np.ndarray]) -> None:
    tifffile.imwrite(path, np.stack(slices, axis=0), photometric="minisblack")


def test_prepare_align_slices_and_match_points_integration(tmp_path: Path) -> None:
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

    corr_path = prepare_brain_align_slices_session(cfg)
    assert corr_path.is_file()
    corr = SliceCorrespondence.load(corr_path)
    assert len(corr.anchors) == 20
    assert all(anchor.confirmed for anchor in corr.anchors)

    match_data = load_brain_match_points_data(cfg)
    assert match_data.slice_correspondence is not None
    assert match_data.session.atlas_slice_indices is not None
    assert any(int(v) > 0 for v in match_data.session.atlas_slice_indices)
