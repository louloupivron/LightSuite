"""Tests for manifest-only multiresolution geometry."""

from __future__ import annotations

import json
import threading
from pathlib import Path

import numpy as np
import pytest
import SimpleITK as sitk
import tifffile
import yaml

from lightsuite.config.loader import load_multires_config
from lightsuite.exceptions import StageCancelledError
from lightsuite.multires.config_models import MultiresGeometryCheckLevel
from lightsuite.multires.manifest import save_pair_manifest
from lightsuite.multires.models import MANIFEST_FORMAT, ManifestVolumeSpec, MultiresPairManifest
from lightsuite.multires.runner import check_multires_geometry
from lightsuite.multires.spec_geometry import (
    alignment_metrics_from_specs,
    manifest_geometry_report_from_spec,
    overlap_physical_bounds_from_specs,
    physical_bounds_from_spec,
)
from lightsuite.multires.volume import load_manifest_xy_slice, manifest_geometry_report
from lightsuite.reporter import stage_cancellation


def _write_plane_stack(folder: Path, shape_zyx: tuple[int, int, int]) -> None:
    folder.mkdir(parents=True, exist_ok=True)
    nz, ny, nx = shape_zyx
    for iz in range(nz):
        arr = np.full((ny, nx), iz + 1, dtype=np.uint16)
        tifffile.imwrite(folder / f"Z{iz:04d}.tif", arr)


def _overview_spec(path: Path) -> ManifestVolumeSpec:
    return ManifestVolumeSpec(
        volume_path=str(path),
        shape_zyx=[4, 8, 8],
        spacing_um=[2.0, 2.0, 2.0],
        origin_um=[0.0, 0.0, 0.0],
    )


def _roi_spec(path: Path) -> ManifestVolumeSpec:
    return ManifestVolumeSpec(
        volume_path=str(path),
        shape_zyx=[4, 4, 4],
        spacing_um=[2.0, 2.0, 2.0],
        origin_um=[4.0, 4.0, 0.0],
    )


def test_physical_bounds_from_spec_matches_sitk() -> None:
    spec = _overview_spec(Path("overview"))
    image = sitk.GetImageFromArray(np.zeros((4, 8, 8), dtype=np.float32))
    image.SetSpacing((2.0, 2.0, 2.0))
    image.SetOrigin((0.0, 0.0, 0.0))

    from lightsuite.multires.geometry import physical_bounds

    sitk_min, sitk_max = physical_bounds(image)
    spec_min, spec_max = physical_bounds_from_spec(spec)
    np.testing.assert_allclose(sitk_min, spec_min)
    np.testing.assert_allclose(sitk_max, spec_max)


def test_manifest_geometry_report_uses_spec_only() -> None:
    spec = _roi_spec(Path("roi"))
    report = manifest_geometry_report("roi", spec)
    assert report["shape_zyx"] == (4, 4, 4)
    assert report["phys_min"][0] == pytest.approx(4.0)


def test_spec_crop_matches_sitk_crop(tmp_path: Path) -> None:
    from lightsuite.multires.geometry import crop_to_physical_box
    from lightsuite.multires.spec_geometry import crop_index_range_from_physical_box

    overview_dir = tmp_path / "overview"
    roi_dir = tmp_path / "roi"
    _write_plane_stack(overview_dir, (4, 8, 8))
    _write_plane_stack(roi_dir, (4, 4, 4))

    manifest = MultiresPairManifest(
        format=MANIFEST_FORMAT,
        sample_name="sample",
        pair_label="pair",
        overview=_overview_spec(overview_dir),
        roi=_roi_spec(roi_dir),
    )
    overlap_min, overlap_max = overlap_physical_bounds_from_specs(manifest.overview, manifest.roi)
    spec_start, spec_size = crop_index_range_from_physical_box(
        manifest.overview,
        overlap_min,
        overlap_max,
    )

    image = sitk.GetImageFromArray(np.zeros((4, 8, 8), dtype=np.float32))
    image.SetSpacing((2.0, 2.0, 2.0))
    image.SetOrigin((0.0, 0.0, 0.0))
    _, sitk_start = crop_to_physical_box(image, overlap_min, overlap_max)
    assert spec_start == sitk_start
    assert spec_size == [4, 4, 4]


def test_alignment_metrics_from_specs() -> None:
    overview = _overview_spec(Path("overview"))
    roi = _roi_spec(Path("roi"))
    overlap_min, overlap_max = overlap_physical_bounds_from_specs(overview, roi)
    metrics = alignment_metrics_from_specs(overview, roi, overlap_min=overlap_min, overlap_max=overlap_max)
    assert metrics["center_offset_um"] == [0.0, 0.0, 0.0]
    assert metrics["center_offset_norm_um"] == pytest.approx(0.0)
    assert metrics["roi_overlap_fraction"] == pytest.approx(1.0)
    assert metrics["overview_overlap_fraction"] == pytest.approx(0.1836734693877551)


def test_check_multires_geometry_metadata_only(tmp_path: Path) -> None:
    overview_dir = tmp_path / "overview"
    roi_dir = tmp_path / "roi"
    _write_plane_stack(overview_dir, (4, 8, 8))
    _write_plane_stack(roi_dir, (4, 4, 4))

    manifest = MultiresPairManifest(
        format=MANIFEST_FORMAT,
        sample_name="sample",
        pair_label="pair",
        overview=_overview_spec(overview_dir),
        roi=_roi_spec(roi_dir),
    )
    manifest_path = tmp_path / "pair.json"
    save_pair_manifest(manifest, manifest_path)

    config_path = tmp_path / "multires.yaml"
    config_path.write_text(
        yaml.safe_dump(
            {
                "sample": {"name": "sample", "save_path": str(tmp_path / "out")},
                "multires": {
                    "vendor": {"suite": "manifest"},
                    "pair_manifest": str(manifest_path),
                },
            }
        ),
        encoding="utf-8",
    )
    cfg = load_multires_config(config_path)
    checkpoint = check_multires_geometry(cfg, level=MultiresGeometryCheckLevel.METADATA_ONLY)

    geometry_dir = tmp_path / "out" / "geometry" / "pair"
    assert (geometry_dir / "geometry_report.json").is_file()
    assert (geometry_dir / "fov_overlap.png").is_file()
    assert checkpoint.overlap_box_um is not None

    report = json.loads((geometry_dir / "geometry_report.json").read_text(encoding="utf-8"))
    assert report["roi"]["phys_min"][0] == pytest.approx(4.0)
    assert report["roi"]["phys_center"][0] == pytest.approx(7.0)
    assert report["alignment_metrics"]["center_offset_norm_um"] == pytest.approx(0.0)


def test_check_multires_geometry_respects_cancel_event(tmp_path: Path) -> None:
    overview_dir = tmp_path / "overview"
    roi_dir = tmp_path / "roi"
    _write_plane_stack(overview_dir, (4, 8, 8))
    _write_plane_stack(roi_dir, (4, 4, 4))

    manifest = MultiresPairManifest(
        format=MANIFEST_FORMAT,
        sample_name="sample",
        pair_label="pair",
        overview=_overview_spec(overview_dir),
        roi=_roi_spec(roi_dir),
    )
    manifest_path = tmp_path / "pair.json"
    save_pair_manifest(manifest, manifest_path)

    config_path = tmp_path / "multires.yaml"
    config_path.write_text(
        yaml.safe_dump(
            {
                "sample": {"name": "sample", "save_path": str(tmp_path / "out")},
                "multires": {
                    "vendor": {"suite": "manifest"},
                    "pair_manifest": str(manifest_path),
                },
            }
        ),
        encoding="utf-8",
    )
    cfg = load_multires_config(config_path)
    cancel = threading.Event()
    cancel.set()
    with stage_cancellation(cancel):
        with pytest.raises(StageCancelledError):
            check_multires_geometry(cfg, level=MultiresGeometryCheckLevel.METADATA_ONLY)


def test_load_manifest_xy_crop_matches_plane_per_file(tmp_path: Path) -> None:
    folder = tmp_path / "stack"
    _write_plane_stack(folder, (3, 5, 6))
    spec = ManifestVolumeSpec(
        volume_path=str(folder),
        shape_zyx=[3, 5, 6],
        spacing_um=[1.0, 1.0, 1.0],
        origin_um=[0.0, 0.0, 0.0],
    )
    from lightsuite.multires.volume import load_manifest_xy_crop, load_manifest_xy_slice

    full_plane = load_manifest_xy_slice(spec, physical_um=(1.0, 1.0, 1.0))
    crop = load_manifest_xy_crop(
        spec,
        z_index=1,
        start_xyz=[1, 1, 1],
        crop_size_xyz=[3, 2, 1],
    )
    assert crop.shape == (2, 3)
    np.testing.assert_array_equal(crop, full_plane[1:3, 1:4])


def test_load_manifest_xy_slice_plane_per_file(tmp_path: Path) -> None:
    folder = tmp_path / "stack"
    _write_plane_stack(folder, (3, 5, 6))
    spec = ManifestVolumeSpec(
        volume_path=str(folder),
        shape_zyx=[3, 5, 6],
        spacing_um=[1.0, 1.0, 1.0],
        origin_um=[0.0, 0.0, 0.0],
    )
    sl = load_manifest_xy_slice(spec, physical_um=(0.0, 0.0, 1.0))
    assert sl.shape == (5, 6)
    assert sl[0, 0] == pytest.approx(2.0)
