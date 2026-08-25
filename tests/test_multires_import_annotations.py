"""Tests for multires ROI→overview annotation import and registration inspect."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest
import tifffile
import yaml

from lightsuite.config.loader import load_multires_config
from lightsuite.gui.inspect_registration_multires import (
    contrast_limits,
    discover_multires_registration_inspect_paths,
    run_multires_inspect_registration,
)
from lightsuite.import_.models import ImportedPoints
from lightsuite.multires.checkpoint import MultiresRegOptsCheckpoint, multires_checkpoint_path
from lightsuite.multires.import_annotations import (
    MultiresAnnotationImporter,
    sample_reference_from_spec,
)
from lightsuite.multires.manifest import save_pair_manifest
from lightsuite.multires.models import MANIFEST_FORMAT, ManifestVolumeSpec, MultiresPairManifest


def _spec(path: Path, shape_zyx, *, spacing=(1.0, 1.0, 1.0), origin=(0.0, 0.0, 0.0)):
    nz, ny, nx = shape_zyx
    return ManifestVolumeSpec(
        volume_path=str(path),
        shape_zyx=[nz, ny, nx],
        spacing_um=list(spacing),
        origin_um=list(origin),
    )


def _write_config(tmp_path: Path, manifest_path: Path) -> Path:
    config_path = tmp_path / "cfg.yaml"
    config_path.write_text(
        yaml.safe_dump(
            {
                "sample": {"name": "s", "save_path": str(tmp_path / "out")},
                "multires": {
                    "vendor": {"suite": "manifest"},
                    "pair_manifest": str(manifest_path),
                },
                "import": {
                    "annotations": [
                        {
                            "format": "points_csv",
                            "path": str(tmp_path / "pts.csv"),
                            "label": "cells",
                        }
                    ]
                },
            }
        ),
        encoding="utf-8",
    )
    return config_path


def _build_pair(tmp_path: Path) -> Path:
    ov = tmp_path / "overview.tif"
    roi = tmp_path / "roi.tif"
    tifffile.imwrite(ov, np.zeros((6, 8, 10), dtype=np.uint16))
    tifffile.imwrite(roi, np.zeros((4, 5, 6), dtype=np.uint16))
    manifest = MultiresPairManifest(
        format=MANIFEST_FORMAT,
        sample_name="s",
        pair_label="p",
        overview=_spec(ov, (6, 8, 10)),
        roi=_spec(roi, (4, 5, 6), spacing=(0.5, 0.5, 0.5), origin=(1.0, 1.0, 1.0)),
    )
    manifest_path = tmp_path / "pair.json"
    save_pair_manifest(manifest, manifest_path)
    return manifest_path


def test_multires_config_accepts_import_block(tmp_path: Path) -> None:
    cfg = load_multires_config(_write_config(tmp_path, _build_pair(tmp_path)))
    assert cfg.import_config is not None
    assert cfg.import_config.annotations[0].label == "cells"


def test_sample_reference_from_spec_uses_roi_grid(tmp_path: Path) -> None:
    spec = _spec(tmp_path / "roi.tif", (4, 5, 6), spacing=(0.5, 0.5, 2.0))
    reference = sample_reference_from_spec(spec, sample_name="s")
    # ManifestVolumeSpec is (Z, Y, X); SampleReference is (Y, X, Z).
    assert reference.shape_yxz == [5, 6, 4]
    assert reference.voxel_um == [0.5, 0.5, 2.0]
    assert reference.index_base == 1


def test_importer_reference_bounds_reject_out_of_roi_points(tmp_path: Path) -> None:
    from lightsuite.import_.adapters import prepare_points_for_sample

    importer = MultiresAnnotationImporter(
        roi_spec=_spec(tmp_path / "roi.tif", (4, 5, 6)),
        overview_spec=_spec(tmp_path / "overview.tif", (6, 8, 10)),
        transform_paths=[tmp_path / "T.txt"],
        crop_start_index=[0, 0, 0],
        crop_size_xyz=[4, 4, 4],
        output_dir=tmp_path / "out",
        temp_dir=tmp_path / "tmp",
    )
    points = ImportedPoints(
        label="cells",
        coordinates=np.array([[1.0, 1.0, 1.0], [6.0, 5.0, 4.0], [99.0, 1.0, 1.0]]),
    )
    prepared = prepare_points_for_sample(points, reference=importer.reference)
    assert prepared.coordinates.shape[0] == 2


def test_points_warp_maps_through_manifest_geometry(tmp_path: Path, monkeypatch) -> None:
    """Identity transform in physical space still maps ROI indices onto overview indices."""
    from lightsuite.multires import import_annotations as mod

    roi_spec = _spec(
        tmp_path / "roi.tif",
        (4, 5, 6),
        spacing=(0.5, 0.5, 0.5),
        origin=(1.0, 1.0, 1.0),
    )
    overview_spec = _spec(tmp_path / "overview.tif", (6, 8, 10), spacing=(1.0, 1.0, 1.0))

    def fake_transformix(*, points_xyz, transform_path, output_dir):
        return np.asarray(points_xyz, dtype=float)

    monkeypatch.setattr(
        "lightsuite.registration.elastix.runner.run_transformix_physical_points",
        fake_transformix,
    )

    # ROI 1-based (1,1,1) → 0-based (0,0,0) → physical (1,1,1) → overview 0-based (1,1,1).
    out = mod.warp_points_roi_to_overview(
        np.array([[1.0, 1.0, 1.0]]),
        roi_spec=roi_spec,
        overview_spec=overview_spec,
        transform_path=tmp_path / "T.txt",
        temp_dir=tmp_path / "tmp",
    )
    assert np.allclose(out[0], [2.0, 2.0, 2.0])


def _write_checkpoint(tmp_path: Path, manifest_path: Path, *, full_canvas: bool) -> Path:
    save_path = tmp_path / "out"
    save_path.mkdir(parents=True, exist_ok=True)
    reg_dir = save_path / "elastix"
    reg_dir.mkdir(parents=True, exist_ok=True)

    crop = reg_dir / "p_overview_crop.tif"
    roi_reg = reg_dir / "p_roi_registered.tif"
    tifffile.imwrite(crop, np.ones((3, 4, 5), dtype=np.float32))
    tifffile.imwrite(roi_reg, np.ones((3, 4, 5), dtype=np.float32) * 2)

    full_path = None
    if full_canvas:
        full_path = reg_dir / "p_roi_in_full_overview.tif"
        tifffile.imwrite(full_path, np.zeros((6, 8, 10), dtype=np.float32))

    checkpoint = MultiresRegOptsCheckpoint(
        sample_name="s",
        pair_label="p",
        pair_manifest_path=str(manifest_path),
        experiment_slug="p",
        overview_volume_path=str(tmp_path / "overview.tif"),
        roi_volume_path=str(tmp_path / "roi.tif"),
        transform_paths=[str(reg_dir / "TransformParameters.0.txt")],
        cropped_overview_path=str(crop),
        registered_roi_path=str(roi_reg),
        registered_roi_full_overview_path=str(full_path) if full_path else None,
        crop_start_index=[0, 0, 0],
        reference_channel="488",
    )
    checkpoint.save(multires_checkpoint_path(save_path))
    return save_path


def test_discover_registration_inspect_defaults_to_crop(tmp_path: Path) -> None:
    manifest_path = _build_pair(tmp_path)
    config_path = _write_config(tmp_path, manifest_path)
    _write_checkpoint(tmp_path, manifest_path, full_canvas=True)

    cfg = load_multires_config(config_path)
    paths = discover_multires_registration_inspect_paths(cfg)
    assert not paths.full_overview
    assert paths.overview_path.name == "p_overview_crop.tif"
    assert "488" in paths.registered_roi_paths


def test_discover_registration_inspect_full_overview(tmp_path: Path) -> None:
    manifest_path = _build_pair(tmp_path)
    config_path = _write_config(tmp_path, manifest_path)
    _write_checkpoint(tmp_path, manifest_path, full_canvas=True)

    cfg = load_multires_config(config_path)
    paths = discover_multires_registration_inspect_paths(cfg, full_overview=True)
    assert paths.full_overview
    assert paths.overview_path.name == "overview.tif"
    assert paths.registered_roi_paths["488"].name == "p_roi_in_full_overview.tif"


def test_discover_registration_inspect_full_overview_missing_canvas(tmp_path: Path) -> None:
    manifest_path = _build_pair(tmp_path)
    config_path = _write_config(tmp_path, manifest_path)
    _write_checkpoint(tmp_path, manifest_path, full_canvas=False)

    cfg = load_multires_config(config_path)
    with pytest.raises(FileNotFoundError, match="write_full_overview_canvas"):
        discover_multires_registration_inspect_paths(cfg, full_overview=True)


def test_inspect_registration_headless_loads_layers(tmp_path: Path) -> None:
    manifest_path = _build_pair(tmp_path)
    config_path = _write_config(tmp_path, manifest_path)
    _write_checkpoint(tmp_path, manifest_path, full_canvas=True)

    cfg = load_multires_config(config_path)
    paths = run_multires_inspect_registration(cfg, headless=True)
    assert paths.overview_path.is_file()


def test_inspect_registration_requires_checkpoint(tmp_path: Path) -> None:
    manifest_path = _build_pair(tmp_path)
    cfg = load_multires_config(_write_config(tmp_path, manifest_path))
    with pytest.raises(FileNotFoundError, match="multires register"):
        discover_multires_registration_inspect_paths(cfg)


def test_contrast_limits_handles_flat_volume() -> None:
    lo, hi = contrast_limits(np.zeros((4, 4, 4), dtype=np.float32))
    assert hi > lo


def test_import_summary_written(tmp_path: Path, monkeypatch) -> None:
    """The shared orchestrator writes one summary JSON next to the outputs."""
    from lightsuite.import_.models import AnnotationImportResult
    from lightsuite.import_.orchestrator import write_import_summary

    out = tmp_path / "out"
    out.mkdir()
    summary = write_import_summary(
        out,
        [AnnotationImportResult(label="cells", kind="points", n_input=3, n_atlas=2)],
    )
    payload = json.loads(summary.read_text(encoding="utf-8"))
    assert payload[0]["label"] == "cells"
    assert payload[0]["n_atlas"] == 2
