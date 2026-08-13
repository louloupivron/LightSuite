"""Tests for spinal cord region-stats orchestration."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import tifffile

from lightsuite.analysis.cord_runner import (
    cord_analysis_requested,
    resolve_cord_stats_spaces,
    run_cord_region_stats,
)
from lightsuite.config.models import CordAtlasConfig, CordRegistrationConfig, CordSampleConfig, SpinalCordPipelineConfig
from lightsuite.export.cord_registered import REGISTERED_ANNOTATION_FILENAME


def _minimal_config(tmp_path: Path) -> SpinalCordPipelineConfig:
    save_path = tmp_path / "registered"
    vr = save_path / "volume_registered"
    vr.mkdir(parents=True)

    annotation = np.zeros((4, 4, 2), dtype=np.uint16)
    annotation[:, :, 0] = 7
    annotation[:, :, 1] = 7
    tifffile.imwrite(vr / REGISTERED_ANNOTATION_FILENAME, annotation, imagej=True)

    intensity = np.ones((4, 4, 2), dtype=np.uint16) * 120
    tifffile.imwrite(vr / "chan01_channel1.tiff", intensity, imagej=True)

    return SpinalCordPipelineConfig(
        sample=CordSampleConfig(
            name="test_cord",
            source={"format": "tiff_stack", "path": str(tmp_path)},
            scratch=tmp_path / "scratch",
            save_path=save_path,
            voxel_um=[1.8, 1.8, 1.8],
        ),
        atlas=CordAtlasConfig(atlas_dir=tmp_path / "atlas"),
        registration=CordRegistrationConfig(),
    )


def test_run_cord_region_stats_intensity_only(tmp_path: Path) -> None:
    atlas_dir = tmp_path / "atlas"
    atlas_dir.mkdir()
    (atlas_dir / "Template.tif").write_bytes(b"")
    (atlas_dir / "Annotation.tif").write_bytes(b"")
    pd.DataFrame({"Segment": ["C1"], "Start": [1], "End": [2]}).to_csv(
        atlas_dir / "Segments.csv", index=False
    )
    pd.DataFrame(
        {
            "id": [7],
            "acronym": ["a7"],
            "name": ["Region7"],
            "parent_ID": [1],
            "parent_acronym": ["GM"],
        }
    ).to_csv(atlas_dir / "Atlas_Regions.csv", index=False)

    config = _minimal_config(tmp_path)
    result = run_cord_region_stats(
        config,
        count_points=False,
        parcellate_intensities=True,
    )

    assert result.combined_path is not None
    assert result.combined_path.is_file()
    assert result.intensity_channels == [1]
    assert result.n_rows > 0

    df = pd.read_csv(result.combined_path)
    assert "median_intensity" in df["metric"].values
    assert (df["segment"] == "C1").any()


def test_resolve_cord_stats_spaces_intersects_export(tmp_path: Path) -> None:
    config = _minimal_config(tmp_path)
    config.analysis.stats_spaces = ["atlas", "sample"]
    assert resolve_cord_stats_spaces(["atlas"], config) == ["atlas"]
    assert resolve_cord_stats_spaces(["sample"], config) == ["sample"]


def test_cord_analysis_requested_respects_flags(tmp_path: Path) -> None:
    config = _minimal_config(tmp_path)
    assert cord_analysis_requested(config)
    config.analysis.parcellate_intensities = False
    config.analysis.count_points = False
    assert not cord_analysis_requested(config)
