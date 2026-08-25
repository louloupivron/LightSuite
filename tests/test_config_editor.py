"""Tests for GUI config form helpers."""

from __future__ import annotations

from pathlib import Path

import pytest

from lightsuite.config.loader import parse_config_yaml, write_config_yaml
from lightsuite.exceptions import LightsuiteConfigError
from lightsuite.gui.config_form_data import (
    brain_atlas_help_text,
    brain_atlas_provider_options,
    brain_form_from_raw,
    brain_form_to_raw,
    default_local_atlas_resolution_um,
    import_annotations_from_raw,
    import_annotations_to_raw,
    load_template_raw,
    multires_form_from_raw,
    multires_form_to_raw,
    parse_orientation_text,
    resolve_brain_brainglobe_form_fields,
    spinal_form_from_raw,
    spinal_form_to_raw,
    try_validate_config_dict,
)
from lightsuite.atlas.brainglobe_backend import BrainGlobeAtlasEntry


def test_parse_config_yaml_rejects_empty() -> None:
    with pytest.raises(LightsuiteConfigError, match="non-empty"):
        parse_config_yaml("")


def test_parse_config_yaml_rejects_invalid_syntax() -> None:
    with pytest.raises(LightsuiteConfigError, match="Invalid YAML"):
        parse_config_yaml("sample: [")


def test_write_config_yaml_roundtrip(tmp_path: Path) -> None:
    text = "sample:\n  name: mouse\n"
    out = write_config_yaml(tmp_path / "cfg.yaml", text)
    assert out.read_text(encoding="utf-8") == text
    assert parse_config_yaml(text)["sample"]["name"] == "mouse"


def test_try_validate_config_dict_brain(tmp_path: Path) -> None:
    (tmp_path / "src").mkdir()
    (tmp_path / "scratch").mkdir()
    (tmp_path / "results").mkdir()
    (tmp_path / "atlas").mkdir()
    data = {
        "sample": {
            "name": "mouse",
            "source": {
                "format": "tiff_stack",
                "path": str(tmp_path / "src"),
                "tiff_type": "channelperfile",
            },
            "scratch": str(tmp_path / "scratch"),
            "save_path": str(tmp_path / "results"),
            "voxel_um": [5.0, 5.0, 5.0],
        },
        "atlas": {"provider": "allen", "resolution_um": 10, "atlas_dir": str(tmp_path / "atlas")},
        "registration": {"resolution_um": 20, "channel_primary": 1},
        "detection": {"enabled": False},
        "compute": {"workers": 1},
        "export": {"registered_volume_format": "tiff"},
    }
    result = try_validate_config_dict(data)
    assert not isinstance(result, LightsuiteConfigError)
    workflow, config = result
    assert workflow == "brain"
    assert config.sample.name == "mouse"


def test_try_validate_config_dict_reports_missing_paths() -> None:
    _, raw = load_template_raw("brain")
    result = try_validate_config_dict(raw)
    assert isinstance(result, LightsuiteConfigError)


def test_load_template_raw_unknown_workflow() -> None:
    with pytest.raises(ValueError, match="Unknown workflow"):
        load_template_raw("not_a_workflow")


def test_load_template_raw_multires_is_blank_starter() -> None:
    workflow, raw = load_template_raw("multires")
    assert workflow == "multires"
    assert raw["sample"]["name"] == "example_multires"
    assert raw["multires"]["vendor"]["suite"] == "mesospim"
    assert "pair_manifest" not in (raw["multires"] or {})
    assert "488" in raw["multires"]["channels"]
    assert raw["multires"]["geometry_mode"] == "metadata"
    state = multires_form_from_raw(raw)
    assert state.vendor_suite == "mesospim"
    assert state.reference_channel == "488"
    assert len(state.channels) == 2


def test_tooltips_for_workflow_brain_and_spinal_differ_on_atlas() -> None:
    from lightsuite.gui.config_form_tooltips import tooltips_for_workflow

    brain = tooltips_for_workflow("brain")
    spinal = tooltips_for_workflow("spinal")
    multires = tooltips_for_workflow("multires")
    assert "atlas_dir" in brain
    assert "Fiederling" in spinal["atlas_dir"]
    assert "pair_manifest" in multires
    assert brain["atlas_dir"] != spinal["atlas_dir"]


def test_brain_form_roundtrip_preserves_extra_keys() -> None:
    raw = {
        "sample": {
            "name": "mouse",
            "source": {"format": "tiff_stack", "path": "/data", "tiff_type": "channelperfile"},
            "scratch": "/scratch",
            "save_path": "/out",
            "voxel_um": [1.0, 1.0, 1.0],
        },
        "atlas": {"provider": "allen", "resolution_um": 10, "atlas_dir": "/atlas"},
        "registration": {"resolution_um": 20, "channel_primary": 1},
        "export": {"write_pyramid": True},
    }
    state = brain_form_from_raw(raw)
    state.sample_name = "renamed"
    updated = brain_form_to_raw(state, raw)
    assert updated["sample"]["name"] == "renamed"
    assert updated["export"]["write_pyramid"] is True


def test_brain_atlas_provider_options_local_vs_brainglobe() -> None:
    local_ids = [provider for provider, _label in brain_atlas_provider_options("files")]
    bg_ids = [provider for provider, _label in brain_atlas_provider_options("brainglobe")]
    assert local_ids == ["allen", "perens"]
    assert bg_ids == ["allen", "perens", "princeton", "rat"]


def test_brain_atlas_help_text_mentions_orientation_switch() -> None:
    perens_bg = brain_atlas_help_text("brainglobe", "perens")
    assert "check-orientation" in perens_bg
    allen_local = brain_atlas_help_text("files", "allen")
    assert "average_template_10.nii.gz" in allen_local


def test_brain_form_roundtrip_brainglobe_source() -> None:
    raw = {
        "sample": {
            "name": "mouse",
            "source": {"format": "tiff_stack", "path": "/data", "tiff_type": "channelperfile"},
            "scratch": "/scratch",
            "save_path": "/out",
            "voxel_um": [1.0, 1.0, 1.0],
        },
        "atlas": {
            "provider": "princeton",
            "source": "brainglobe",
            "brainglobe_name": "princeton_mouse_20um",
            "resolution_um": 20,
        },
        "registration": {"resolution_um": 20, "channel_primary": 1},
    }
    state = brain_form_from_raw(raw)
    assert state.atlas_source == "brainglobe"
    assert state.atlas_provider == "princeton"
    assert state.atlas_resolution_um == 20.0
    updated = brain_form_to_raw(state, raw)
    assert updated["atlas"]["source"] == "brainglobe"
    assert updated["atlas"]["brainglobe_name"] == "princeton_mouse_20um"
    assert updated["atlas"]["resolution_um"] == 20.0
    assert "atlas_dir" not in updated["atlas"]


def test_resolve_brain_brainglobe_form_fields_uses_catalog_resolution() -> None:
    catalog = [
        BrainGlobeAtlasEntry(
            name="allen_mouse_10um_v4",
            provider="allen",
            resolution_um=10.0,
            label="Allen CCF (mouse) — 10 µm",
        )
    ]
    name, provider, resolution_um = resolve_brain_brainglobe_form_fields(
        brainglobe_name="allen_mouse_10um",
        provider="allen",
        resolution_um=25.0,
        catalog=catalog,
    )
    assert name == "allen_mouse_10um_v4"
    assert provider == "allen"
    assert resolution_um == 10.0


def test_default_local_atlas_resolution_um() -> None:
    assert default_local_atlas_resolution_um("allen") == 10.0
    assert default_local_atlas_resolution_um("perens") == 20.0


def test_brain_form_files_source_clears_brainglobe_name() -> None:
    raw = {
        "sample": {
            "name": "mouse",
            "source": {"format": "tiff_stack", "path": "/data", "tiff_type": "channelperfile"},
            "scratch": "/scratch",
            "save_path": "/out",
            "voxel_um": [1.0, 1.0, 1.0],
        },
        "atlas": {
            "provider": "allen",
            "source": "brainglobe",
            "brainglobe_name": "allen_mouse_10um",
            "resolution_um": 10,
            "atlas_dir": "/atlas",
        },
        "registration": {"resolution_um": 20, "channel_primary": 1},
    }
    state = brain_form_from_raw(raw)
    state.atlas_source = "files"
    updated = brain_form_to_raw(state, raw)
    assert "source" not in updated["atlas"]
    assert "brainglobe_name" not in updated["atlas"]
    assert updated["atlas"]["atlas_dir"] == "/atlas"


def test_spinal_form_channel_list() -> None:
    raw = {
        "sample": {
            "name": "cord",
            "source": {
                "format": "tiff_stack",
                "tiff_type": "planeperfile",
                "channels": ["/ch0", "/ch1"],
            },
            "scratch": "/scratch",
            "save_path": "/out",
            "voxel_um": [1.8, 1.8, 1.8],
        },
        "atlas": {"atlas_dir": "/atlas"},
        "registration": {"resolution_um": 20, "channel_primary": 1},
        "compute": {"workers": 2},
    }
    state = spinal_form_from_raw(raw)
    assert state.use_channel_list is True
    assert state.channel_paths == ["/ch0", "/ch1"]
    updated = spinal_form_to_raw(state, raw)
    assert updated["registration"]["channel_primary"] == 1


def test_spinal_form_preserves_longitudinal_direction() -> None:
    raw = {
        "sample": {
            "name": "cord",
            "source": {"format": "tiff_stack", "path": "/data", "tiff_type": "planeperfile"},
            "scratch": "/scratch",
            "save_path": "/out",
            "voxel_um": [1.8, 1.8, 1.8],
        },
        "atlas": {"atlas_dir": "/atlas"},
        "registration": {
            "resolution_um": 20,
            "channel_primary": 1,
            "longitudinal_direction": "caudorostral",
        },
        "compute": {"workers": 2},
    }
    state = spinal_form_from_raw(raw)
    state.workers = 8
    updated = spinal_form_to_raw(state, raw)
    assert updated["registration"]["longitudinal_direction"] == "caudorostral"
    assert updated["compute"]["workers"] == 8


def test_parse_orientation_text() -> None:
    assert parse_orientation_text("1, 2, 3") == (1, 2, 3)
    assert parse_orientation_text("[1, -3, 2]") == (1, -3, 2)
    assert parse_orientation_text("") is None
    assert parse_orientation_text("bad") is None


def test_brain_form_roundtrip_registration_advanced_fields() -> None:
    raw = {
        "sample": {
            "name": "mouse",
            "source": {"format": "tiff_stack", "path": "/data", "tiff_type": "channelperfile"},
            "scratch": "/scratch",
            "save_path": "/out",
            "voxel_um": [1.0, 1.0, 1.0],
        },
        "atlas": {"provider": "allen", "resolution_um": 10, "atlas_dir": "/atlas"},
        "registration": {
            "resolution_um": 20,
            "channel_primary": 1,
            "channel_secondary": 2,
            "bspline_spatial_scale_mm": 0.5,
            "control_point_weight": 0.3,
            "augment_points": True,
            "dual_channel_mi_weight_primary": 0.5,
            "dual_channel_mi_weight_secondary": 0.7,
            "orientation": [1, -3, 2],
            "canvas_mode": "pad",
        },
        "import": {
            "write_csv": True,
            "annotations": [
                {"format": "points_csv", "path": "/pts.csv", "label": "cells"},
            ],
        },
    }
    state = brain_form_from_raw(raw)
    assert state.bspline_spatial_scale_mm == 0.5
    assert state.augment_points is True
    assert state.orientation == (1, -3, 2)
    assert state.canvas_mode == "pad"
    assert len(state.import_annotations) == 1
    updated = brain_form_to_raw(state, raw)
    assert updated["registration"]["dual_channel_mi_weight_primary"] == 0.5
    assert updated["registration"]["dual_channel_mi_weight_secondary"] == 0.7
    assert "dual_channel_mi_weight_autofluor" not in updated["registration"]
    assert updated["registration"]["orientation"] == [1, -3, 2]
    assert updated["import"]["annotations"][0]["label"] == "cells"


def test_brain_form_reads_legacy_dual_channel_mi_weight_keys() -> None:
    raw = {
        "sample": {
            "name": "mouse",
            "source": {"format": "tiff_stack", "path": "/data", "tiff_type": "channelperfile"},
            "scratch": "/scratch",
            "save_path": "/out",
            "voxel_um": [1.0, 1.0, 1.0],
        },
        "atlas": {"provider": "allen", "resolution_um": 10, "atlas_dir": "/atlas"},
        "registration": {
            "resolution_um": 20,
            "channel_primary": 1,
            "dual_channel_mi_weight_autofluor": 0.6,
            "dual_channel_mi_weight_signal": 0.8,
        },
    }
    state = brain_form_from_raw(raw)
    assert state.dual_channel_mi_weight_primary == 0.6
    assert state.dual_channel_mi_weight_secondary == 0.8


def test_brain_form_preserves_zero_dual_channel_mi_weights() -> None:
    raw = {
        "sample": {
            "name": "mouse",
            "source": {"format": "tiff_stack", "path": "/data", "tiff_type": "channelperfile"},
            "scratch": "/scratch",
            "save_path": "/out",
            "voxel_um": [1.0, 1.0, 1.0],
        },
        "atlas": {"provider": "allen", "resolution_um": 10, "atlas_dir": "/atlas"},
        "registration": {
            "resolution_um": 20,
            "channel_primary": 1,
            "dual_channel_mi_weight_primary": 0.0,
            "dual_channel_mi_weight_secondary": 0.0,
        },
    }
    state = brain_form_from_raw(raw)
    assert state.dual_channel_mi_weight_primary == 0.0
    assert state.dual_channel_mi_weight_secondary == 0.0


def test_brain_form_roundtrip_analysis_metrics() -> None:
    raw = {
        "sample": {
            "name": "mouse",
            "source": {"format": "tiff_stack", "path": "/data", "tiff_type": "channelperfile"},
            "scratch": "/scratch",
            "save_path": "/out",
            "voxel_um": [1.0, 1.0, 1.0],
        },
        "atlas": {"provider": "allen", "resolution_um": 10, "atlas_dir": "/atlas"},
        "registration": {"resolution_um": 20, "channel_primary": 1},
        "analysis": {
            "intensity_metrics": ["mean_intensity", "std"],
            "stats_spaces": ["sample"],
        },
    }
    state = brain_form_from_raw(raw)
    assert state.intensity_metrics == ["mean_intensity", "std"]
    assert state.stats_spaces == ["sample"]
    updated = brain_form_to_raw(state, raw)
    assert updated["analysis"]["intensity_metrics"] == ["mean_intensity", "std"]
    assert updated["analysis"]["stats_spaces"] == ["sample"]
    assert updated["analysis"]["count_points"] is True


def test_brain_form_stats_spaces_defaults_to_atlas() -> None:
    raw = {
        "sample": {
            "name": "mouse",
            "source": {"format": "tiff_stack", "path": "/data", "tiff_type": "channelperfile"},
            "scratch": "/scratch",
            "save_path": "/out",
            "voxel_um": [1.0, 1.0, 1.0],
        },
        "atlas": {"provider": "allen", "resolution_um": 10, "atlas_dir": "/atlas"},
        "registration": {"resolution_um": 20, "channel_primary": 1},
    }
    state = brain_form_from_raw(raw)
    assert state.stats_spaces == ["atlas"]
    state.stats_spaces = ["atlas", "sample"]
    updated = brain_form_to_raw(state, raw)
    assert updated["analysis"]["stats_spaces"] == ["atlas", "sample"]


def test_spinal_form_roundtrip_control_point_and_import() -> None:
    raw = {
        "sample": {
            "name": "cord",
            "source": {"format": "tiff_stack", "path": "/data", "tiff_type": "planeperfile"},
            "scratch": "/scratch",
            "save_path": "/out",
            "voxel_um": [1.8, 1.8, 1.8],
        },
        "atlas": {"atlas_dir": "/atlas"},
        "registration": {"resolution_um": 20, "channel_primary": 1, "control_point_weight": 0.5},
        "import": {
            "annotations": [{"format": "mask_tiff", "path": "/mask.tif", "label": "region"}],
        },
    }
    state = spinal_form_from_raw(raw)
    assert state.control_point_weight == 0.5
    assert state.import_annotations[0].format == "mask_tiff"
    updated = spinal_form_to_raw(state, raw)
    assert updated["registration"]["control_point_weight"] == 0.5
    assert updated["import"]["annotations"][0]["path"] == "/mask.tif"


def test_spinal_form_always_enables_parcellate_intensities() -> None:
    raw = {
        "sample": {
            "name": "cord",
            "source": {"format": "tiff_stack", "path": "/data", "tiff_type": "planeperfile"},
            "scratch": "/scratch",
            "save_path": "/out",
            "voxel_um": [1.8, 1.8, 1.8],
        },
        "atlas": {"atlas_dir": "/atlas"},
        "analysis": {"parcellate_intensities": False},
    }
    updated = spinal_form_to_raw(spinal_form_from_raw(raw), raw)
    assert updated["analysis"]["parcellate_intensities"] is True


def test_multires_form_roundtrip_geometry_and_registration() -> None:
    raw = {
        "sample": {"name": "tg14", "save_path": "/out", "scratch": "/scratch"},
        "multires": {
            "pair_label": "pair",
            "geometry_mode": "hybrid",
            "landmarks": {"fit_mode": "affine"},
            "registration": {
                "reference_channel": "488",
                "experiment_name": "tg14_hybrid",
                "overlap_margin_um": -10.0,
                "write_full_overview_canvas": False,
            },
            "channels": {"488": {"overview": "/ov.tif", "roi": "/roi.tif"}},
        },
        "import": {
            "annotations": [{"format": "points_csv", "path": "/a.csv"}],
        },
    }
    state = multires_form_from_raw(raw)
    assert state.geometry_mode == "hybrid"
    assert state.landmark_fit_mode == "affine"
    assert state.overlap_margin_um == -10.0
    assert state.experiment_name == "tg14_hybrid"
    assert state.write_full_overview_canvas is False
    assert state.vendor_suite == "mesospim"
    updated = multires_form_to_raw(state, raw)
    assert updated["multires"]["geometry_mode"] == "hybrid"
    assert updated["multires"]["landmarks"]["fit_mode"] == "affine"
    assert updated["multires"]["registration"]["overlap_margin_um"] == -10.0
    assert updated["multires"]["registration"]["experiment_name"] == "tg14_hybrid"
    assert updated["multires"]["vendor"]["suite"] == "mesospim"
    assert updated["import"]["annotations"][0]["path"] == "/a.csv"


def test_multires_form_drops_mesospim_geometry_for_smartspim_vendor() -> None:
    raw = {
        "sample": {"name": "s", "save_path": "/out"},
        "multires": {
            "vendor": {"suite": "smartspim"},
            "geometry_mode": "metadata",
            "mesospim_geometry": {
                "overview": {"lateral_flip": [1, -1]},
                "roi": {"lateral_flip": [1, -1]},
            },
            "channels": {"640": {"overview": "/ov", "roi": "/roi"}},
        },
    }
    state = multires_form_from_raw(raw)
    state.vendor_suite = "smartspim"
    updated = multires_form_to_raw(state, raw)
    assert updated["multires"]["vendor"]["suite"] == "smartspim"
    assert "mesospim_geometry" not in updated["multires"]


def test_multires_form_roundtrip_vendor_manifest_and_custom() -> None:
    manifest_raw = {
        "sample": {"name": "s", "save_path": "/out"},
        "multires": {
            "vendor": {"suite": "manifest"},
            "pair_manifest": "/data/pair.json",
        },
    }
    manifest_state = multires_form_from_raw(manifest_raw)
    assert manifest_state.vendor_suite == "manifest"
    assert manifest_state.pair_manifest == "/data/pair.json"
    manifest_updated = multires_form_to_raw(manifest_state, manifest_raw)
    assert manifest_updated["multires"]["vendor"]["suite"] == "manifest"

    custom_raw = {
        "sample": {"name": "s", "save_path": "/out"},
        "multires": {
            "vendor": {"suite": "custom", "custom_entry": "/scripts/build.py"},
            "channels": {"488": {"overview": "/ov.tif", "roi": "/roi.tif"}},
        },
    }
    custom_state = multires_form_from_raw(custom_raw)
    assert custom_state.vendor_suite == "custom"
    assert custom_state.vendor_custom_entry == "/scripts/build.py"
    custom_updated = multires_form_to_raw(custom_state, custom_raw)
    assert custom_updated["multires"]["vendor"]["custom_entry"] == "/scripts/build.py"

    meta_raw = {
        "sample": {"name": "s", "save_path": "/out"},
        "multires": {
            "vendor": {"suite": "smartspim"},
            "channels": {
                "488": {
                    "overview": "/ov.tif",
                    "roi": "/roi.tif",
                    "overview_meta_path": "/ov_meta.json",
                    "roi_meta_path": "/roi_meta.txt",
                },
            },
        },
    }
    meta_state = multires_form_from_raw(meta_raw)
    assert meta_state.channels[0].overview_meta == "/ov_meta.json"
    assert meta_state.channels[0].roi_meta == "/roi_meta.txt"
    meta_updated = multires_form_to_raw(meta_state, meta_raw)
    ch = meta_updated["multires"]["channels"]["488"]
    assert ch["overview_meta_path"] == "/ov_meta.json"
    assert ch["roi_meta_path"] == "/roi_meta.txt"


def test_multires_form_prunes_stale_apply_transform_to() -> None:
    """Template leftovers must not block save when channels are reduced."""
    raw = {
        "sample": {"name": "s", "save_path": "/out"},
        "multires": {
            "pair_manifest": "/legacy/gilda_pair.json",
            "vendor": {"suite": "smartspim"},
            "channels": {"488": {"overview": "/ov.tif", "roi": "/roi.tif"}},
            "registration": {
                "reference_channel": "488",
                "apply_transform_to": ["555", "647"],
            },
        },
    }
    state = multires_form_from_raw(raw)
    updated = multires_form_to_raw(state, raw)
    assert "apply_transform_to" not in updated["multires"]["registration"]
    assert "pair_manifest" not in updated["multires"]
    assert updated["multires"]["vendor"]["suite"] == "smartspim"


def test_multires_form_roundtrips_geometry_check_level() -> None:
    raw = {
        "sample": {"name": "s", "save_path": "/out"},
        "multires": {
            "vendor": {"suite": "smartspim"},
            "channels": {"488": {"overview": "/ov", "roi": "/roi"}},
            "registration": {
                "reference_channel": "488",
                "geometry_check_level": "slice-qc",
            },
        },
    }
    state = multires_form_from_raw(raw)
    assert state.geometry_check_level == "slice-qc"
    updated = multires_form_to_raw(state, raw)
    assert updated["multires"]["registration"]["geometry_check_level"] == "slice-qc"

    state.geometry_check_level = "full"
    updated = multires_form_to_raw(state, raw)
    assert "geometry_check_level" not in updated["multires"]["registration"]


def test_import_annotations_to_raw_clears_empty_rows() -> None:
    raw = {"import": {"write_csv": True, "annotations": [{"format": "points_csv", "path": "/a.csv"}]}}
    rows = import_annotations_from_raw(raw)
    updated = import_annotations_to_raw(
        [rows[0], rows[0].__class__(format="points_csv", path="", label="")],
        raw,
    )
    assert len(updated["import"]["annotations"]) == 1
