"""Tests for brain init-registration."""

from __future__ import annotations

from pathlib import Path

import nibabel as nib
import numpy as np
import tifffile
import yaml

from lightsuite.config.loader import load_config
from lightsuite.preprocess.brain import preprocess_lightsheet_volume
from lightsuite.registration.init_brain import initialize_brain_registration


def _write_channel_stack(path: Path, slices: list[np.ndarray]) -> None:
    tifffile.imwrite(path, np.stack(slices, axis=0), photometric="minisblack")


def _write_mini_atlas(atlas_dir: Path, shape: tuple[int, int, int] = (20, 20, 20)) -> None:
    atlas_dir.mkdir(parents=True, exist_ok=True)
    yy, xx, zz = np.mgrid[0 : shape[0], 0 : shape[1], 0 : shape[2]]
    template = (255 * np.exp(-((yy - 10) ** 2 + (xx - 10) ** 2 + (zz - 10) ** 2) / 50)).astype(
        np.float32
    )
    annotation = np.zeros(shape, dtype=np.float32)
    annotation[5:15, 5:15, 5:15] = 1
    nib.save(nib.Nifti1Image(template, np.eye(4)), atlas_dir / "average_template_10.nii.gz")
    nib.save(nib.Nifti1Image(annotation, np.eye(4)), atlas_dir / "annotation_10.nii.gz")


def test_initialize_brain_registration(tmp_path: Path) -> None:
    data_dir = tmp_path / "data"
    data_dir.mkdir()
    scratch = tmp_path / "scratch"
    save = tmp_path / "results"
    save.mkdir()
    atlas_dir = tmp_path / "atlas"
    _write_mini_atlas(atlas_dir)

    # Sample with strong edges for gradient-based point extraction
    yy, xx = np.mgrid[0:32, 0:32]
    ring = (((yy - 16) ** 2 + (xx - 16) ** 2) > 36) & (((yy - 16) ** 2 + (xx - 16) ** 2) < 100)
    base = (150 + 100 * ring).astype(np.uint16)
    slices = [base + z * 5 for z in range(6)]
    _write_channel_stack(data_dir / "ch1.tif", slices)

    config_data = {
        "sample": {
            "name": "test",
            "source": {"path": str(data_dir), "tiff_type": "channelperfile"},
            "scratch": str(scratch),
            "save_path": str(save),
            "voxel_um": [10.0, 10.0, 10.0],
        },
        "atlas": {"provider": "allen", "resolution_um": 10, "atlas_dir": str(atlas_dir)},
        "registration": {
            "resolution_um": 20,
            "channel_primary": 1,
            "orientation": [1, 2, 3],
            "cloud_threshold": 0.5,
        },
        "detection": {"enabled": False},
    }
    config_path = tmp_path / "config.yaml"
    config_path.write_text(yaml.dump(config_data), encoding="utf-8")
    cfg = load_config(config_path)

    preprocess_lightsheet_volume(cfg)
    result = initialize_brain_registration(cfg)

    assert result.permute_sample_to_atlas == [1, 2, 3]
    assert result.original_trans is not None
    assert len(result.original_trans) == 4
    assert result.downfac_reg == 0.5
    assert result.autocpsample is not None
    assert (save / "dim1_initial_registration.png").is_file()
    assert (save / "dim1_initial_registration.png").stat().st_size > 10_000
