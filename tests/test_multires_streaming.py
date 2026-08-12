"""Tests for streamed multires overlap preparation (no full-volume load)."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import SimpleITK as sitk
import tifffile
import yaml

from lightsuite.config.loader import load_multires_config
from lightsuite.multires.geometry import prepare_registration_pair
from lightsuite.multires.manifest import save_pair_manifest
from lightsuite.multires.models import MANIFEST_FORMAT, ManifestVolumeSpec, MultiresPairManifest
from lightsuite.multires.prepare import prepare_multires_registration_pair
from lightsuite.multires.volume import (
    load_manifest_volume,
    load_manifest_xy_plane_at_z_index,
    load_manifest_xyz_crop,
    stream_resample_to_reference,
    write_embedded_crop_canvas,
)


def _write_plane_stack(folder: Path, arr: np.ndarray) -> None:
    folder.mkdir(parents=True, exist_ok=True)
    for z, plane in enumerate(arr):
        tifffile.imwrite(folder / f"Z{z:06d}.tif", plane.astype(np.float32))


def _spec(
    path: Path,
    arr: np.ndarray,
    *,
    spacing=(1.0, 1.0, 1.0),
    origin=(0.0, 0.0, 0.0),
) -> ManifestVolumeSpec:
    nz, ny, nx = arr.shape
    return ManifestVolumeSpec(
        volume_path=str(path),
        shape_zyx=[nz, ny, nx],
        spacing_um=list(spacing),
        origin_um=list(origin),
    )


def test_load_manifest_single_page_hyperstack(tmp_path: Path) -> None:
    arr = np.arange(4 * 6 * 8, dtype=np.uint16).reshape(4, 6, 8)
    path = tmp_path / "stack.tif"
    tifffile.imwrite(path, arr, imagej=True)

    spec = _spec(path, arr)
    for iz in (0, 2, 3):
        plane = load_manifest_xy_plane_at_z_index(spec, iz)
        np.testing.assert_array_equal(plane, arr[iz].astype(np.float32))

    crop = load_manifest_xyz_crop(spec, start_xyz=[1, 1, 1], crop_size_xyz=[3, 2, 2])
    np.testing.assert_array_equal(sitk.GetArrayFromImage(crop), arr[1:3, 1:3, 1:4])


def test_load_manifest_xyz_crop_matches_full_load(tmp_path: Path) -> None:
    arr = np.arange(5 * 8 * 10, dtype=np.float32).reshape(5, 8, 10)
    folder = tmp_path / "vol"
    _write_plane_stack(folder, arr)
    spec = _spec(folder, arr, origin=(10.0, 20.0, 30.0))

    full = load_manifest_volume(spec)
    crop = load_manifest_xyz_crop(spec, start_xyz=[2, 1, 1], crop_size_xyz=[4, 3, 2])
    expected = sitk.RegionOfInterest(full, [4, 3, 2], [2, 1, 1])
    np.testing.assert_allclose(
        sitk.GetArrayFromImage(crop),
        sitk.GetArrayFromImage(expected),
    )
    assert crop.GetOrigin() == expected.GetOrigin()
    assert crop.GetSpacing() == expected.GetSpacing()


def test_stream_resample_matches_full_resample(tmp_path: Path) -> None:
    overview = np.zeros((6, 12, 12), dtype=np.float32)
    overview[:, 2:10, 2:10] = 1.0
    roi = np.zeros((6, 12, 12), dtype=np.float32)
    # Bright blob shifted relative to overview content but same physical frame.
    zz, yy, xx = np.ogrid[:6, :12, :12]
    roi[((zz - 3) ** 2 + (yy - 5) ** 2 + (xx - 5) ** 2) < 6] = 5.0

    ov_dir = tmp_path / "overview"
    roi_dir = tmp_path / "roi"
    _write_plane_stack(ov_dir, overview)
    _write_plane_stack(roi_dir, roi)
    ov_spec = _spec(ov_dir, overview)
    roi_spec = _spec(roi_dir, roi)

    full_ov = load_manifest_volume(ov_spec)
    full_roi = load_manifest_volume(roi_spec)
    fixed_full, moving_full, _box, start = prepare_registration_pair(full_ov, full_roi)

    fixed_stream = load_manifest_xyz_crop(
        ov_spec,
        start_xyz=start,
        crop_size_xyz=list(fixed_full.GetSize()),
    )
    moving_stream = stream_resample_to_reference(roi_spec, fixed_stream, z_chunk=2)

    np.testing.assert_allclose(
        sitk.GetArrayFromImage(fixed_stream),
        sitk.GetArrayFromImage(fixed_full),
        atol=1e-5,
    )
    np.testing.assert_allclose(
        sitk.GetArrayFromImage(moving_stream),
        sitk.GetArrayFromImage(moving_full),
        atol=1e-4,
    )


def test_prepare_streams_without_full_volumes(tmp_path: Path) -> None:
    overview = np.random.default_rng(0).random((4, 16, 16), dtype=np.float32)
    roi = np.random.default_rng(1).random((4, 16, 16), dtype=np.float32)
    ov_dir = tmp_path / "overview"
    roi_dir = tmp_path / "roi"
    _write_plane_stack(ov_dir, overview)
    _write_plane_stack(roi_dir, roi)

    manifest = MultiresPairManifest(
        format=MANIFEST_FORMAT,
        sample_name="s",
        pair_label="p",
        overview=_spec(ov_dir, overview),
        roi=_spec(roi_dir, roi),
    )
    manifest_path = tmp_path / "pair.json"
    save_pair_manifest(manifest, manifest_path)

    config_path = tmp_path / "cfg.yaml"
    config_path.write_text(
        yaml.safe_dump(
            {
                "sample": {"name": "s", "save_path": str(tmp_path / "out")},
                "multires": {"pair_manifest": str(manifest_path)},
            }
        ),
        encoding="utf-8",
    )
    cfg = load_multires_config(config_path)
    prepared = prepare_multires_registration_pair(cfg)
    assert prepared.fixed_cropped.GetSize() == prepared.moving.GetSize()
    assert prepared.overview_spec.shape_zyx == [4, 16, 16]
    # Peak memory of prepared pair should be the two crops, not 2× full stacks.
    crop_voxels = int(np.prod(prepared.fixed_cropped.GetSize()))
    full_voxels = 4 * 16 * 16
    assert crop_voxels <= full_voxels


def test_write_embedded_crop_canvas_plane_by_plane(tmp_path: Path) -> None:
    overview_spec = ManifestVolumeSpec(
        volume_path="unused",
        shape_zyx=[5, 6, 7],
        spacing_um=[1.0, 1.0, 2.0],
        origin_um=[0.0, 0.0, 0.0],
    )
    crop_arr = np.arange(2 * 3 * 4, dtype=np.float32).reshape(2, 3, 4)
    crop = sitk.GetImageFromArray(crop_arr)
    crop.SetSpacing((1.0, 1.0, 2.0))
    out = tmp_path / "canvas.tif"
    write_embedded_crop_canvas(overview_spec, crop, [1, 2, 1], out)

    with tifffile.TiffFile(out) as tf:
        assert len(tf.pages) == 5
        pages = [page.asarray() for page in tf.pages]
    canvas = np.stack(pages, axis=0)
    assert canvas.shape == (5, 6, 7)
    assert canvas[0].sum() == 0
    np.testing.assert_array_equal(canvas[1:3, 2:5, 1:5], crop_arr)

    from lightsuite.gui.brain_multires_link import load_overview_native_volume_yxz
    from lightsuite.registration.volume import load_tiff_volume_zyx

    np.testing.assert_array_equal(load_tiff_volume_zyx(out), canvas.astype(np.float32))
    np.testing.assert_array_equal(
        load_overview_native_volume_yxz(out),
        np.moveaxis(canvas, 0, -1).astype(np.float32),
    )
