"""Tests for multiresolution landmark match-points helpers."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import tifffile
import yaml

from lightsuite.config.loader import load_multires_config
from lightsuite.gui.multires_data import (
    MultiresSliceSource,
    _default_z_indices,
    load_multires_match_points_data,
    prepare_multires_match_points_session,
)
from lightsuite.multires.landmark_session import MultiresLandmarkSession, default_landmark_session_path
from lightsuite.multires.manifest import save_pair_manifest
from lightsuite.multires.models import MANIFEST_FORMAT, ManifestVolumeSpec, MultiresPairManifest
from lightsuite.multires.volume import load_manifest_xy_plane_at_z_index


def _write_plane_stack(folder: Path, arr: np.ndarray) -> None:
    folder.mkdir(parents=True, exist_ok=True)
    for z, plane in enumerate(arr):
        tifffile.imwrite(folder / f"Z{z:06d}.tif", plane.astype(np.uint16))


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


def _write_pair_config(
    tmp_path: Path,
    overview: np.ndarray,
    roi: np.ndarray,
    *,
    overview_origin=(0.0, 0.0, 0.0),
    roi_origin=(2.0, 2.0, 1.0),
) -> Path:
    ov_dir = tmp_path / "overview"
    roi_dir = tmp_path / "roi"
    _write_plane_stack(ov_dir, overview)
    _write_plane_stack(roi_dir, roi)
    manifest = MultiresPairManifest(
        format=MANIFEST_FORMAT,
        sample_name="s",
        pair_label="p",
        overview=_spec(ov_dir, overview, origin=overview_origin),
        roi=_spec(roi_dir, roi, origin=roi_origin),
    )
    manifest_path = tmp_path / "pair.json"
    save_pair_manifest(manifest, manifest_path)
    config_path = tmp_path / "cfg.yaml"
    config_path.write_text(
        yaml.safe_dump(
            {
                "sample": {"name": "s", "save_path": str(tmp_path / "out")},
                "multires": {
                    "pair_manifest": str(manifest_path),
                    "geometry_mode": "hybrid",
                    "landmarks": {"fit_mode": "similarity", "min_pairs": 3},
                },
            }
        ),
        encoding="utf-8",
    )
    return config_path


def test_load_manifest_xy_plane_at_z_index(tmp_path: Path) -> None:
    arr = np.arange(4 * 6 * 8, dtype=np.uint16).reshape(4, 6, 8)
    folder = tmp_path / "vol"
    _write_plane_stack(folder, arr)
    spec = _spec(folder, arr)
    plane = load_manifest_xy_plane_at_z_index(spec, 2)
    assert plane.shape == (6, 8)
    np.testing.assert_array_equal(plane, arr[2].astype(np.float32))


def test_prepare_multires_match_points_session_headless(tmp_path: Path) -> None:
    overview = np.zeros((6, 12, 12), dtype=np.uint16)
    roi = np.zeros((4, 8, 8), dtype=np.uint16)
    config_path = _write_pair_config(tmp_path, overview, roi)
    cfg = load_multires_config(config_path)
    path = prepare_multires_match_points_session(cfg)
    expected = default_landmark_session_path(tmp_path / "out", "p")
    assert path == expected
    assert path.is_file()
    assert path == prepare_multires_match_points_session(cfg)


def test_multires_slice_source_lazy_read_and_cache(tmp_path: Path) -> None:
    arr = np.arange(4 * 8 * 8, dtype=np.uint16).reshape(4, 8, 8)
    folder = tmp_path / "vol"
    _write_plane_stack(folder, arr)
    source = MultiresSliceSource.from_spec(_spec(folder, arr), manifest_dir=tmp_path)
    sl0 = source.read_display_slice(0)
    sl0_cached = source.read_display_slice(0)
    assert sl0.shape == (8, 8)
    assert sl0 is sl0_cached
    sl2 = source.read_display_slice(2)
    assert sl2.shape == (8, 8)
    assert not np.allclose(sl0, sl2)


def test_default_z_indices_use_overlap_center(tmp_path: Path) -> None:
    overview = _spec(
        tmp_path / "ov",
        np.zeros((10, 20, 20), dtype=np.uint16),
        origin=(0.0, 0.0, 0.0),
    )
    roi = _spec(
        tmp_path / "roi",
        np.zeros((5, 8, 8), dtype=np.uint16),
        origin=(4.0, 4.0, 3.0),
    )
    oz, rz = _default_z_indices(overview, roi)
    assert oz == 5
    assert rz == 2


def test_match_points_resumes_at_latest_annotated_slices(tmp_path: Path) -> None:
    overview = np.zeros((10, 20, 20), dtype=np.uint16)
    roi = np.zeros((8, 12, 12), dtype=np.uint16)
    config_path = _write_pair_config(tmp_path, overview, roi)
    cfg = load_multires_config(config_path)
    from lightsuite.multires.resolve import resolve_pair_manifest

    manifest, _ = resolve_pair_manifest(cfg)
    session_path = cfg.multires.resolved_landmark_session_path(cfg.sample.save_path, manifest)
    session_path.parent.mkdir(parents=True, exist_ok=True)
    MultiresLandmarkSession(
        overview_points_zyx=[[7.0, 10.0, 11.0], [8.0, 12.0, 13.0]],
        roi_points_zyx=[[5.0, 4.0, 5.0], [6.0, 6.0, 7.0]],
    ).save(session_path)

    data = load_multires_match_points_data(cfg)
    default_oz, default_rz = _default_z_indices(data.overview_spec, data.roi_spec)
    assert (default_oz, default_rz) != (8, 6)
    assert data.initial_overview_z == 8
    assert data.initial_roi_z == 6


def test_match_points_hybrid_uses_metadata_overlap_crop(tmp_path: Path) -> None:
    # Overview spacing like SmartSPIM (~4 µm): 200 µm margin ≈ 50 voxels, not full FOV.
    overview = np.arange(6 * 400 * 400, dtype=np.uint16).reshape(6, 400, 400)
    roi = np.arange(4 * 40 * 40, dtype=np.uint16).reshape(4, 40, 40)
    ov_dir = tmp_path / "overview"
    roi_dir = tmp_path / "roi"
    _write_plane_stack(ov_dir, overview)
    _write_plane_stack(roi_dir, roi)
    manifest = MultiresPairManifest(
        format=MANIFEST_FORMAT,
        sample_name="s",
        pair_label="p",
        overview=_spec(ov_dir, overview, spacing=(4.0, 4.0, 4.0), origin=(0.0, 0.0, 0.0)),
        roi=_spec(roi_dir, roi, spacing=(1.0, 1.0, 1.0), origin=(200.0, 200.0, 4.0)),
    )
    manifest_path = tmp_path / "pair.json"
    save_pair_manifest(manifest, manifest_path)
    config_path = tmp_path / "cfg.yaml"
    config_path.write_text(
        yaml.safe_dump(
            {
                "sample": {"name": "s", "save_path": str(tmp_path / "out")},
                "multires": {
                    "pair_manifest": str(manifest_path),
                    "geometry_mode": "hybrid",
                },
            }
        ),
        encoding="utf-8",
    )
    cfg = load_multires_config(config_path)
    data = load_multires_match_points_data(cfg)
    assert data.crop_mode is True
    assert data.overview.is_xy_cropped is True
    assert data.overview.crop_size_yx[0] < 400
    assert data.overview.crop_size_yx[1] < 400

    plane = data.overview.read_display_slice(data.initial_overview_z)
    assert plane.shape == data.overview.crop_size_yx

    iy0, ix0 = data.overview.xy_origin_yx
    volume_pts = [[float(data.initial_overview_z), float(iy0 + 3), float(ix0 + 4)]]
    display = data.overview.display_xy_from_volume_zyx(volume_pts, data.initial_overview_z)
    assert np.allclose(display[0], [3.0, 4.0])
    restored = data.overview.volume_zyx_from_display_xy(
        display,
        data.initial_overview_z,
        [],
    )
    assert restored == volume_pts


def test_match_points_metadata_mode_uses_full_volumes(tmp_path: Path) -> None:
    overview = np.arange(6 * 400 * 400, dtype=np.uint16).reshape(6, 400, 400)
    roi = np.arange(4 * 40 * 40, dtype=np.uint16).reshape(4, 40, 40)
    ov_dir = tmp_path / "overview"
    roi_dir = tmp_path / "roi"
    _write_plane_stack(ov_dir, overview)
    _write_plane_stack(roi_dir, roi)
    manifest = MultiresPairManifest(
        format=MANIFEST_FORMAT,
        sample_name="s",
        pair_label="p",
        overview=_spec(ov_dir, overview, spacing=(4.0, 4.0, 4.0), origin=(0.0, 0.0, 0.0)),
        roi=_spec(roi_dir, roi, spacing=(1.0, 1.0, 1.0), origin=(200.0, 200.0, 4.0)),
    )
    manifest_path = tmp_path / "pair.json"
    save_pair_manifest(manifest, manifest_path)
    config_path = tmp_path / "cfg.yaml"
    config_path.write_text(
        yaml.safe_dump(
            {
                "sample": {"name": "s", "save_path": str(tmp_path / "out")},
                "multires": {
                    "pair_manifest": str(manifest_path),
                    "geometry_mode": "metadata",
                },
            }
        ),
        encoding="utf-8",
    )
    cfg = load_multires_config(config_path)
    data = load_multires_match_points_data(cfg)
    assert data.crop_mode is False
    assert data.overview.is_xy_cropped is False
    assert data.overview.crop_size_yx == (400, 400)


def test_legacy_landmarks_geometry_mode_maps_to_hybrid(tmp_path: Path) -> None:
    overview = np.zeros((6, 12, 12), dtype=np.uint16)
    roi = np.zeros((4, 8, 8), dtype=np.uint16)
    config_path = _write_pair_config(tmp_path, overview, roi)
    raw = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    raw["multires"]["geometry_mode"] = "landmarks"
    config_path.write_text(yaml.safe_dump(raw), encoding="utf-8")
    cfg = load_multires_config(config_path)
    from lightsuite.multires.config_models import MultiresGeometryMode

    assert cfg.multires.geometry_mode == MultiresGeometryMode.HYBRID


def test_match_points_applies_config_lateral_flip(tmp_path: Path) -> None:
    """Match-points must rebuild pair geometry from YAML lateral_flip, not a stale JSON."""
    import tifffile
    from lightsuite.mesospim.meta import meta_path_for_tiff

    overview = tmp_path / "overview.tif"
    roi = tmp_path / "roi.tif"
    tifffile.imwrite(overview, np.zeros((5, 16, 16), dtype=np.uint16), imagej=True)
    tifffile.imwrite(roi, np.zeros((5, 16, 16), dtype=np.uint16), imagej=True)
    meta_body = "\n".join(
        [
            "[Pixelsize in um] 5.0",
            "[x_pos] 100.0",
            "[y_pos] 500.0",
            "[z_start] 0.0",
            "[z_end] 8.0",
            "[z_stepsize] 2.0",
            "[z_planes] 5",
            "[x_pixels] 16",
            "[y_pixels] 16",
        ]
    )
    meta_path_for_tiff(overview).write_text(meta_body, encoding="utf-8")
    meta_path_for_tiff(roi).write_text(meta_body, encoding="utf-8")

    def _cfg(flip: list[int]):
        config_path = tmp_path / f"cfg_{flip[0]}_{flip[1]}.yaml"
        config_path.write_text(
            yaml.safe_dump(
                {
                    "sample": {"name": "s", "save_path": str(tmp_path / "out")},
                    "multires": {
                        "vendor": {"suite": "mesospim"},
                        "pair_label": "p",
                        "channels": {"488": {"overview": str(overview), "roi": str(roi)}},
                        "mesospim_geometry": {
                            "overview": {"lateral_flip": flip},
                            "roi": {"lateral_flip": flip},
                        },
                        "registration": {"reference_channel": "488"},
                    },
                }
            ),
            encoding="utf-8",
        )
        return load_multires_config(config_path)

    default = load_multires_match_points_data(_cfg([1, -1]))
    mirrored = load_multires_match_points_data(_cfg([-1, -1]))
    assert default.overview_spec.direction != mirrored.overview_spec.direction
    assert default.roi_spec.direction != mirrored.roi_spec.direction
