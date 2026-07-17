"""Tests for multiresolution pair manifests."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import tifffile
import yaml

from lightsuite.config.loader import load_multires_config
from lightsuite.multires.manifest import load_pair_manifest, save_pair_manifest
from lightsuite.multires.models import MANIFEST_FORMAT, ManifestVolumeSpec, MultiresPairManifest
from lightsuite.multires.prepare import prepare_multires_registration_pair
from lightsuite.multires.vendor.mesospim import build_mesospim_pair_manifest
from lightsuite.mesospim.meta import meta_path_for_tiff


def _write_stack(path: Path, shape_zyx: tuple[int, int, int]) -> None:
    arr = np.zeros(shape_zyx, dtype=np.uint16)
    tifffile.imwrite(path, arr, imagej=True)


def _write_meta(path: Path) -> None:
    path.write_text(
        "\n".join(
            [
                "[Pixelsize in um] 5.0",
                "[x_pos] 0.0",
                "[y_pos] 0.0",
                "[z_start] 0.0",
                "[z_end] 8.0",
                "[z_stepsize] 2.0",
                "[z_planes] 5",
                "[x_pixels] 16",
                "[y_pixels] 16",
            ]
        ),
        encoding="utf-8",
    )


def test_build_mesospim_pair_manifest(tmp_path: Path) -> None:
    overview = tmp_path / "overview.tif"
    roi = tmp_path / "roi.tif"
    _write_stack(overview, (5, 16, 16))
    _write_stack(roi, (5, 16, 16))
    _write_meta(meta_path_for_tiff(overview))
    _write_meta(meta_path_for_tiff(roi))

    manifest_path = tmp_path / "pair.json"
    manifest = build_mesospim_pair_manifest(
        sample_name="sample_a",
        pair_label="test_pair",
        overview_path=overview,
        roi_path=roi,
        output_manifest_path=manifest_path,
    )

    assert manifest.format == MANIFEST_FORMAT
    assert manifest.overview.shape_zyx == [5, 16, 16]
    assert manifest.provenance["microscope"] == "mesospim"
    loaded = load_pair_manifest(manifest_path)
    assert loaded.pair_label == "test_pair"


def test_load_multires_config(tmp_path: Path) -> None:
    overview = tmp_path / "overview.tif"
    roi = tmp_path / "roi.tif"
    _write_stack(overview, (5, 16, 16))
    _write_stack(roi, (5, 16, 16))
    _write_meta(meta_path_for_tiff(overview))
    _write_meta(meta_path_for_tiff(roi))

    manifest_path = tmp_path / "pair.json"
    build_mesospim_pair_manifest(
        sample_name="sample_a",
        pair_label="test_pair",
        overview_path=overview,
        roi_path=roi,
        output_manifest_path=manifest_path,
    )

    config_path = tmp_path / "multires.yaml"
    config_path.write_text(
        yaml.safe_dump(
            {
                "sample": {"name": "sample_a", "save_path": str(tmp_path / "out")},
                "multires": {"pair_manifest": str(manifest_path)},
            }
        ),
        encoding="utf-8",
    )

    cfg = load_multires_config(config_path)
    assert cfg.sample.name == "sample_a"
    assert cfg.multires.registration.elastix_stages == ["translation", "rigid"]


def test_prepare_multires_registration_pair(tmp_path: Path) -> None:
    overview = tmp_path / "overview.tif"
    roi = tmp_path / "roi.tif"
    _write_stack(overview, (5, 16, 16))
    _write_stack(roi, (5, 16, 16))
    _write_meta(meta_path_for_tiff(overview))
    _write_meta(meta_path_for_tiff(roi))

    manifest_path = tmp_path / "pair.json"
    build_mesospim_pair_manifest(
        sample_name="sample_a",
        pair_label="overlap_pair",
        overview_path=overview,
        roi_path=roi,
        output_manifest_path=manifest_path,
    )

    config_path = tmp_path / "multires.yaml"
    config_path.write_text(
        yaml.safe_dump(
            {
                "sample": {"name": "sample_a", "save_path": str(tmp_path / "out")},
                "multires": {"pair_manifest": str(manifest_path)},
            }
        ),
        encoding="utf-8",
    )
    cfg = load_multires_config(config_path)
    prepared = prepare_multires_registration_pair(cfg)
    assert prepared.fixed_cropped.GetSize() == prepared.moving.GetSize()
    assert prepared.crop_start_index == [0, 0, 0]


def test_save_and_load_manifest_roundtrip(tmp_path: Path) -> None:
    manifest = MultiresPairManifest(
        format=MANIFEST_FORMAT,
        sample_name="s",
        pair_label="p",
        overview=ManifestVolumeSpec(
            volume_path="overview.tif",
            shape_zyx=[2, 4, 4],
            spacing_um=[1.0, 1.0, 1.0],
            origin_um=[0.0, 0.0, 0.0],
        ),
        roi=ManifestVolumeSpec(
            volume_path="roi.tif",
            shape_zyx=[2, 4, 4],
            spacing_um=[1.0, 1.0, 1.0],
            origin_um=[0.0, 0.0, 0.0],
        ),
    )
    overview = tmp_path / "overview.tif"
    roi = tmp_path / "roi.tif"
    _write_stack(overview, (2, 4, 4))
    _write_stack(roi, (2, 4, 4))
    path = tmp_path / "pair.json"
    save_pair_manifest(manifest, path)
    loaded = load_pair_manifest(path)
    assert loaded.overview.shape_zyx == [2, 4, 4]
