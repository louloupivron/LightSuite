"""CLI entry for the unified Napari GUI."""

from __future__ import annotations

from pathlib import Path

import typer


def register_gui_command(app: typer.Typer) -> None:
    @app.command("gui")
    def gui(
        config: str | None = typer.Option(
            None,
            "--config",
            "-c",
            help="Pipeline YAML config to load on startup.",
        ),
    ) -> None:
        """Open the unified Napari pipeline shell."""
        from lightsuite.gui.shell import launch_lightsuite_gui

        config_path = Path(config).expanduser().resolve() if config is not None else None
        launch_lightsuite_gui(config_path)
