"""Tests for spinal stats/plots GUI helpers."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock

import pandas as pd
import pytest

from lightsuite.analysis.cord_counts import CORD_TIDY_COLUMNS
from lightsuite.analysis.cord_heatmap import PAPER_STRUCTURE_ACRONYMS
from lightsuite.config.models import CordAtlasConfig, CordRegistrationConfig, CordSampleConfig, SpinalCordPipelineConfig
from lightsuite.gui.spinal_stats_plots import CordStatsPlotsPanel, attach_spinal_stats_plots


def _config_with_stats(tmp_path: Path, *, split_hemispheres: bool = False) -> SpinalCordPipelineConfig:
    save_path = tmp_path / "registered"
    vr = save_path / "volume_registered"
    vr.mkdir(parents=True)
    hemispheres = (
        (("left", 10.0), ("right", 30.0), ("whole", 25.0))
        if split_hemispheres
        else (("whole", 42.0),)
    )
    rows = []
    for hemisphere, value in hemispheres:
        rows.append(
            {
                "sample": "s1",
                "channel": 1,
                "atlas": "fiederling",
                "parcellation_index": 201,
                "acronym": PAPER_STRUCTURE_ACRONYMS[0],
                "name": "Lamina I",
                "structure": "DH",
                "division": "GM",
                "segment": "C1",
                "rollup_level": "structure",
                "hemisphere": hemisphere,
                "metric": "median_intensity",
                "value": value,
            }
        )
    pd.DataFrame(rows).reindex(columns=CORD_TIDY_COLUMNS).to_csv(
        vr / "region_stats.csv",
        index=False,
    )
    atlas_dir = tmp_path / "atlas"
    atlas_dir.mkdir()
    pd.DataFrame(
        {
            "Segment": ["C1", "Co2"],
            "Ref_Section": [1, 2],
            "Start": [1, 2],
            "End": [1, 2],
        }
    ).to_csv(atlas_dir / "Segments.csv", index=False)
    return SpinalCordPipelineConfig(
        sample=CordSampleConfig(
            name="test_cord",
            source={"format": "tiff_stack", "path": str(tmp_path)},
            scratch=tmp_path / "scratch",
            save_path=save_path,
            voxel_um=[1.8, 1.8, 1.8],
        ),
        atlas=CordAtlasConfig(atlas_dir=atlas_dir),
        registration=CordRegistrationConfig(),
    )


@pytest.mark.filterwarnings("ignore::DeprecationWarning")
def test_attach_spinal_stats_plots_returns_controller(tmp_path: Path) -> None:
    pytest.importorskip("qtpy")
    from qtpy.QtWidgets import QApplication

    app = QApplication.instance()
    if app is None:
        app = QApplication([])

    config = _config_with_stats(tmp_path)
    viewer = MagicMock()
    controller = attach_spinal_stats_plots(viewer, config)
    assert controller.open_log_message is not None
    controller.mount(viewer)
    viewer.window.add_dock_widget.assert_called()
    controller.teardown(viewer)
    app.processEvents()


@pytest.mark.filterwarnings("ignore::DeprecationWarning")
def test_stats_plots_panel_segment_range_and_csv(tmp_path: Path) -> None:
    pytest.importorskip("qtpy")
    from qtpy.QtWidgets import QApplication

    app = QApplication.instance()
    if app is None:
        app = QApplication([])

    config = _config_with_stats(tmp_path)
    stats_path = config.sample.save_path / "volume_registered" / "region_stats.csv"
    panel = CordStatsPlotsPanel(config, stats_df=pd.read_csv(stats_path), stats_path=stats_path)
    assert panel._segment_start.currentText() == "C1"
    assert panel._segment_end.currentText() == "Co2"
    assert panel._current_segments() == ["C1", "Co2"]
    assert panel._export_kind == "matrix"
    csv_path = tmp_path / "matrix.csv"
    from lightsuite.analysis.cord_heatmap import write_cord_heatmap_matrix_csv

    saved = write_cord_heatmap_matrix_csv(next(iter(panel._export_matrices.values())), csv_path)
    loaded = pd.read_csv(saved, index_col="segment")
    assert loaded.loc["C1", "Lamina_I"] == 42.0
    panel.teardown()
    app.processEvents()


@pytest.mark.filterwarnings("ignore::DeprecationWarning")
def test_stats_plots_panel_compare_view(tmp_path: Path) -> None:
    pytest.importorskip("qtpy")
    from qtpy.QtWidgets import QApplication

    app = QApplication.instance()
    if app is None:
        app = QApplication([])

    config = _config_with_stats(tmp_path, split_hemispheres=True)
    stats_path = config.sample.save_path / "volume_registered" / "region_stats.csv"
    panel = CordStatsPlotsPanel(config, stats_df=pd.read_csv(stats_path), stats_path=stats_path)
    compare_index = panel._view_combo.findData("compare")
    assert compare_index >= 0
    panel._view_combo.setCurrentIndex(compare_index)
    panel._update_plot()
    assert panel._export_kind == "compare"
    assert set(panel._export_matrices) == {"left", "right", "whole"}
    panel.teardown()
    app.processEvents()
