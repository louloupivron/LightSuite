"""Interactive stats / plots panel for spinal cord region tables."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg

from lightsuite.analysis.cord_heatmap import (
    MetricName,
    build_cord_heatmap_figure,
    cord_stats_to_matrix,
    default_heatmap_output_path,
    default_paper_segments,
    discover_cord_plot_options,
    metric_label,
    plot_cord_heatmap,
    resolve_region_stats_path,
)
from lightsuite.config.models import SpinalCordPipelineConfig
from lightsuite.gui.stage_controller import DockStageController, run_attached_stage


def _parse_channel_value(text: str) -> int | str:
    if text.isdigit():
        return int(text)
    return text


def _default_cmap(metric: str) -> str:
    if "intensity" in metric or metric == "std":
        return "hot"
    return "viridis"


class CordStatsPlotsPanel:
    """Qt panel to explore segment × lamina/WM heatmaps from region_stats.csv."""

    def __init__(
        self,
        config: SpinalCordPipelineConfig,
        *,
        stats_df: pd.DataFrame,
        stats_path: Path,
    ) -> None:
        from qtpy.QtWidgets import (
            QCheckBox,
            QComboBox,
            QFileDialog,
            QFormLayout,
            QHBoxLayout,
            QLabel,
            QMessageBox,
            QPushButton,
            QScrollArea,
            QVBoxLayout,
            QWidget,
        )

        self._config = config
        self._stats_df = stats_df
        self._stats_path = stats_path
        self._segments = default_paper_segments(config.atlas.atlas_dir)
        self._options = discover_cord_plot_options(stats_df)
        self._figure = None
        self._canvas: FigureCanvasQTAgg | None = None

        self.widget = QWidget()
        root = QVBoxLayout(self.widget)

        source_label = QLabel(f"Source: {stats_path}")
        source_label.setWordWrap(True)
        root.addWidget(source_label)

        form = QFormLayout()
        self._metric_combo = QComboBox()
        for metric in self._options["metrics"]:
            self._metric_combo.addItem(metric_label(str(metric)), str(metric))
        self._metric_combo.currentIndexChanged.connect(self._on_metric_changed)
        form.addRow("Plot type:", self._metric_combo)

        self._channel_combo = QComboBox()
        form.addRow("Channel / label:", self._channel_combo)

        self._hemisphere_combo = QComboBox()
        for hemisphere in self._options["hemispheres"]:
            self._hemisphere_combo.addItem(str(hemisphere))
        whole_index = self._hemisphere_combo.findText("whole")
        if whole_index >= 0:
            self._hemisphere_combo.setCurrentIndex(whole_index)
        form.addRow("Hemisphere:", self._hemisphere_combo)

        self._rollup_combo = QComboBox()
        rollup_levels = self._options["rollup_levels"]
        preferred = "structure" if "structure" in rollup_levels else rollup_levels[0]
        for level in rollup_levels:
            self._rollup_combo.addItem(str(level))
        rollup_index = self._rollup_combo.findText(preferred)
        if rollup_index >= 0:
            self._rollup_combo.setCurrentIndex(rollup_index)
        form.addRow("Rollup:", self._rollup_combo)

        self._normalize_combo = QComboBox()
        for label, value in (
            ("None", "none"),
            ("Per segment (row)", "row"),
            ("Per region (column)", "column"),
        ):
            self._normalize_combo.addItem(label, value)
        form.addRow("Normalize:", self._normalize_combo)

        self._cmap_combo = QComboBox()
        for cmap in ("hot", "viridis", "magma", "plasma", "inferno", "cividis"):
            self._cmap_combo.addItem(cmap)
        form.addRow("Colormap:", self._cmap_combo)

        self._log_scale_check = QCheckBox("Log color scale")
        form.addRow("", self._log_scale_check)

        root.addLayout(form)

        button_row = QHBoxLayout()
        self._update_button = QPushButton("Update plot")
        self._update_button.clicked.connect(self._update_plot)
        self._save_button = QPushButton("Save PNG…")
        self._save_button.clicked.connect(self._save_png)
        button_row.addWidget(self._update_button)
        button_row.addWidget(self._save_button)
        root.addLayout(button_row)

        self._status_label = QLabel("")
        self._status_label.setWordWrap(True)
        root.addWidget(self._status_label)

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        self._canvas_host = QWidget()
        self._canvas_layout = QVBoxLayout(self._canvas_host)
        scroll.setWidget(self._canvas_host)
        root.addWidget(scroll, stretch=1)

        self._on_metric_changed()
        self._update_plot()

    def _current_metric(self) -> str:
        return str(self._metric_combo.currentData())

    def _current_channel(self) -> int | str:
        return _parse_channel_value(str(self._channel_combo.currentData()))

    def _on_metric_changed(self) -> None:
        metric = self._current_metric()
        channels = self._options["channels_by_metric"].get(metric, [])
        self._channel_combo.blockSignals(True)
        self._channel_combo.clear()
        for channel in channels:
            label = f"Channel {channel}" if channel.isdigit() else channel
            self._channel_combo.addItem(label, channel)
        self._channel_combo.blockSignals(False)

        cmap = _default_cmap(metric)
        cmap_index = self._cmap_combo.findText(cmap)
        if cmap_index >= 0:
            self._cmap_combo.setCurrentIndex(cmap_index)

    def _build_matrix(self) -> pd.DataFrame:
        return cord_stats_to_matrix(
            self._stats_df,
            metric=self._current_metric(),  # type: ignore[arg-type]
            channel=self._current_channel(),
            rollup_level=str(self._rollup_combo.currentText()),
            hemisphere=str(self._hemisphere_combo.currentText()),
            segments=self._segments,
        )

    def _plot_title(self) -> str:
        metric = self._current_metric()
        channel = self._current_channel()
        title = f"{self._config.sample.name} — {metric_label(metric)} (ch {channel})"
        hemisphere = str(self._hemisphere_combo.currentText())
        if hemisphere != "whole":
            title = f"{title}, {hemisphere}"
        return title

    def _update_plot(self) -> None:
        from qtpy.QtWidgets import QMessageBox

        try:
            matrix = self._build_matrix()
            fig = build_cord_heatmap_figure(
                matrix,
                title=self._plot_title(),
                cmap=str(self._cmap_combo.currentText()),
                log_scale=self._log_scale_check.isChecked(),
                normalize=str(self._normalize_combo.currentData()),
            )
        except (ValueError, KeyError) as exc:
            self._status_label.setText(str(exc))
            QMessageBox.warning(self.widget, "Could not plot", str(exc))
            return

        if self._figure is not None:
            plt.close(self._figure)
        if self._canvas is not None:
            self._canvas_layout.removeWidget(self._canvas)
            self._canvas.deleteLater()
            self._canvas = None

        self._figure = fig
        self._canvas = FigureCanvasQTAgg(fig)
        self._canvas_layout.addWidget(self._canvas)
        self._status_label.setText(
            f"{matrix.shape[0]} segments × {matrix.shape[1]} regions"
        )

    def _save_png(self) -> None:
        from qtpy.QtWidgets import QFileDialog, QMessageBox

        metric = self._current_metric()
        channel = self._current_channel()
        default_path = default_heatmap_output_path(
            self._config,
            metric=metric,
            channel=channel,
        )
        path, _ = QFileDialog.getSaveFileName(
            self.widget,
            "Save heatmap",
            str(default_path),
            "PNG images (*.png)",
        )
        if not path:
            return
        try:
            matrix = self._build_matrix()
            saved = plot_cord_heatmap(
                matrix,
                Path(path),
                title=self._plot_title(),
                cmap=str(self._cmap_combo.currentText()),
                log_scale=self._log_scale_check.isChecked(),
                normalize=str(self._normalize_combo.currentData()),
            )
        except (ValueError, OSError) as exc:
            QMessageBox.warning(self.widget, "Could not save", str(exc))
            return
        self._status_label.setText(f"Saved {saved}")

    def teardown(self) -> None:
        if self._figure is not None:
            plt.close(self._figure)
            self._figure = None
        self._canvas = None


def attach_spinal_stats_plots(
    viewer: Any,
    config: SpinalCordPipelineConfig,
    *,
    stats_path: Path | None = None,
) -> DockStageController:
    """Attach the stats/plots dock widgets to a napari viewer."""
    resolved = resolve_region_stats_path(config, input_path=stats_path)
    stats_df = pd.read_csv(resolved)
    panel = CordStatsPlotsPanel(config, stats_df=stats_df, stats_path=resolved)

    def _teardown() -> None:
        panel.teardown()

    return DockStageController(
        dock_widgets=[(panel.widget, "Stats / plots")],
        _teardown_fn=_teardown,
        open_log_message=(
            f"Stats / plots — {config.sample.name} "
            f"({len(stats_df)} rows from {resolved.name})"
        ),
    )


def run_spinal_stats_plots(
    config: SpinalCordPipelineConfig,
    *,
    headless: bool = False,
    stats_path: Path | None = None,
) -> Any:
    """Open the stats/plots panel in napari or the unified shell."""
    if headless:
        msg = "Stats / plots requires an interactive GUI session."
        raise RuntimeError(msg)

    title = f"LightSuite spinal — {config.sample.name} stats"
    return run_attached_stage(
        title,
        lambda viewer: attach_spinal_stats_plots(viewer, config, stats_path=stats_path),
    )
