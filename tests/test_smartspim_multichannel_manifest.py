"""Tests for SmartSPIM multichannel pair manifest building."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import tifffile
import yaml

from lightsuite.config.loader import load_multires_config
from lightsuite.multires.manifest import load_pair_manifest
from lightsuite.multires.models import MANIFEST_FORMAT
from lightsuite.multires.vendor.smartspim import build_smartspim_multichannel_pair_manifest


def _write_smartspim_meta_txt(path: Path, *, num_images: int = 5) -> None:
    path.write_text(
        "Obj\tH_Res\tV_Res\tum/pix\tZ step (um)\n"
        "LCT 3.6x\t2000\t1600\t1.8\t1.8\n"
        f"503010\t592920\t4832\t561\t0\t2\t1\t561F\t{num_images}\n",
        encoding="utf-8",
    )


def _write_smartspim_stack_folder(path: Path, *, num_planes: int = 5) -> None:
    path.mkdir(parents=True, exist_ok=True)
    for iz in range(num_planes):
        plane = np.zeros((1600, 2000), dtype=np.uint16)
        tifffile.imwrite(path / f"plane_{iz:04d}.tif", plane)


def test_build_smartspim_multichannel_pair_manifest(tmp_path: Path) -> None:
    overview_a = tmp_path / "overview_488"
    roi_a = tmp_path / "roi_488"
    overview_b = tmp_path / "overview_561"
    roi_b = tmp_path / "roi_561"
    _write_smartspim_stack_folder(overview_a)
    _write_smartspim_stack_folder(roi_a)
    _write_smartspim_stack_folder(overview_b)
    _write_smartspim_stack_folder(roi_b)
    ov_meta_a = tmp_path / "overview_488_meta.txt"
    roi_meta_a = tmp_path / "roi_488_meta.txt"
    ov_meta_b = tmp_path / "overview_561_meta.txt"
    roi_meta_b = tmp_path / "roi_561_meta.txt"
    _write_smartspim_meta_txt(ov_meta_a)
    _write_smartspim_meta_txt(roi_meta_a)
    _write_smartspim_meta_txt(ov_meta_b)
    _write_smartspim_meta_txt(roi_meta_b)

    manifest_path = tmp_path / "pair.json"
    manifest = build_smartspim_multichannel_pair_manifest(
        sample_name="sample_a",
        pair_label="smartspim_dual",
        channels={
            "488": {"overview": overview_a, "roi": roi_a},
            "561": {"overview": overview_b, "roi": roi_b},
        },
        reference_channel="488",
        overview_meta_by_channel={
            "488": ov_meta_a,
            "561": ov_meta_b,
        },
        roi_meta_by_channel={
            "488": roi_meta_a,
            "561": roi_meta_b,
        },
        output_manifest_path=manifest_path,
    )

    assert manifest.format == MANIFEST_FORMAT
    assert manifest.reference_channel == "488"
    assert set(manifest.channel_names()) == {"488", "561"}
    assert manifest.provenance["microscope"] == "smartspim"
    loaded = load_pair_manifest(manifest_path)
    assert loaded.channels is not None
    assert loaded.channels["561"].roi.volume_path.endswith("roi_561")


def test_resolve_smartspim_multichannel_from_config(tmp_path: Path) -> None:
    overview_a = tmp_path / "overview_488"
    roi_a = tmp_path / "roi_488"
    overview_b = tmp_path / "overview_561"
    roi_b = tmp_path / "roi_561"
    _write_smartspim_stack_folder(overview_a)
    _write_smartspim_stack_folder(roi_a)
    _write_smartspim_stack_folder(overview_b)
    _write_smartspim_stack_folder(roi_b)
    ov_meta = tmp_path / "ov_meta.txt"
    roi_meta_a = tmp_path / "roi_488_meta.txt"
    roi_meta_b = tmp_path / "roi_561_meta.txt"
    _write_smartspim_meta_txt(ov_meta)
    _write_smartspim_meta_txt(roi_meta_a)
    _write_smartspim_meta_txt(roi_meta_b)

    config_path = tmp_path / "config.yaml"
    config_path.write_text(
        yaml.safe_dump(
            {
                "sample": {"name": "sample_a", "save_path": str(tmp_path / "out")},
                "multires": {
                    "vendor": {"suite": "smartspim"},
                    "pair_label": "smartspim_dual",
                    "overview_meta_path": str(ov_meta),
                    "channels": {
                        "488": {
                            "overview": str(overview_a),
                            "roi": str(roi_a),
                            "roi_meta_path": str(roi_meta_a),
                        },
                        "561": {
                            "overview": str(overview_b),
                            "roi": str(roi_b),
                            "roi_meta_path": str(roi_meta_b),
                        },
                    },
                    "registration": {"reference_channel": "488"},
                },
            }
        ),
        encoding="utf-8",
    )

    cfg = load_multires_config(config_path)
    from lightsuite.multires.resolve import resolve_pair_manifest

    manifest, manifest_path = resolve_pair_manifest(cfg)
    assert manifest_path.is_file()
    assert manifest.reference_channel == "488"
    assert set(manifest.channel_names()) == {"488", "561"}
