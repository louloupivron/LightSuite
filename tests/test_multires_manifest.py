"""Tests for multiresolution pair manifests."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import tifffile
import yaml

from lightsuite.config.loader import load_multires_config
from lightsuite.mesospim.config_models import MesospimGeometryConfig
from lightsuite.mesospim.meta import meta_path_for_tiff, parse_mesospim_meta
from lightsuite.multires.manifest import load_pair_manifest, save_pair_manifest
from lightsuite.multires.models import MANIFEST_FORMAT, ManifestVolumeSpec, MultiresPairManifest
from lightsuite.multires.prepare import prepare_multires_registration_pair
from lightsuite.multires.vendor.mesospim import (
    build_mesospim_multichannel_pair_manifest,
    build_mesospim_pair_manifest,
    mesospim_volume_shape,
)
from lightsuite.mesospim.geometry import stitched_mosaic_geometry_fields


def _write_stack(path: Path, shape_zyx: tuple[int, int, int]) -> None:
    arr = np.zeros(shape_zyx, dtype=np.uint16)
    tifffile.imwrite(path, arr, imagej=True)


def _write_meta(path: Path, *, y_pixels: int = 16, x_pixels: int = 16, z_planes: int = 5) -> None:
    path.write_text(
        "\n".join(
            [
                "[Pixelsize in um] 5.0",
                "[x_pos] 100.0",
                "[y_pos] 500.0",
                "[z_start] 0.0",
                "[z_end] 8.0",
                "[z_stepsize] 2.0",
                f"[z_planes] {z_planes}",
                f"[x_pixels] {x_pixels}",
                f"[y_pixels] {y_pixels}",
            ]
        ),
        encoding="utf-8",
    )


def _write_stitched_folder(path: Path, shape_zyx: tuple[int, int, int]) -> None:
    path.mkdir(parents=True, exist_ok=True)
    nz, ny, nx = shape_zyx
    for iz in range(nz):
        plane = np.zeros((ny, nx), dtype=np.uint16)
        tifffile.imwrite(path / f"plane_{iz:04d}.tif", plane)


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


def test_build_mesospim_stitched_overview_manifest(tmp_path: Path) -> None:
    stitched = tmp_path / "RES_stitched"
    roi = tmp_path / "roi.tif"
    anchor_meta = tmp_path / "tile0_meta.txt"
    _write_stitched_folder(stitched, (5, 32, 12))
    _write_stack(roi, (5, 16, 16))
    _write_meta(anchor_meta, y_pixels=16, x_pixels=16, z_planes=5)
    _write_meta(meta_path_for_tiff(roi))

    manifest = build_mesospim_pair_manifest(
        sample_name="sample_a",
        pair_label="stitched_pair",
        overview_path=stitched,
        roi_path=roi,
        overview_meta_path=anchor_meta,
    )

    assert manifest.overview.shape_zyx == [5, 32, 12]
    assert manifest.provenance["overview_layout"] == "stitched_folder"
    assert mesospim_volume_shape(stitched) == (5, 32, 12)

    anchor = parse_mesospim_meta(anchor_meta)
    spacing, origin, direction = stitched_mosaic_geometry_fields(
        (5, 32, 12),
        anchor,
        MesospimGeometryConfig(),
    )
    assert spacing == (5.0, 5.0, 2.0)
    assert len(direction) == 9


def test_build_mesospim_multichannel_manifest(tmp_path: Path) -> None:
    overview_a = tmp_path / "overview_a.tif"
    overview_b = tmp_path / "overview_b.tif"
    roi_a = tmp_path / "roi_a.tif"
    roi_b = tmp_path / "roi_b.tif"
    for path, shape in (
        (overview_a, (5, 16, 16)),
        (overview_b, (5, 16, 16)),
        (roi_a, (5, 16, 16)),
        (roi_b, (5, 16, 16)),
    ):
        _write_stack(path, shape)
        _write_meta(meta_path_for_tiff(path))

    manifest_path = tmp_path / "pair.json"
    manifest = build_mesospim_multichannel_pair_manifest(
        sample_name="sample_a",
        pair_label="dual_channel",
        reference_channel="488",
        channels={
            "488": {"overview": overview_a, "roi": roi_a},
            "561": {"overview": overview_b, "roi": roi_b},
        },
        output_manifest_path=manifest_path,
    )

    assert manifest.reference_channel == "488"
    assert set(manifest.channel_names()) == {"488", "561"}
    assert manifest.non_reference_channels(reference_channel="488") == ["561"]
    loaded = load_pair_manifest(manifest_path)
    assert loaded.overview.volume_path.endswith("overview_a.tif")
