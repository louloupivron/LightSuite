"""Register run / stages commands on workflow Typer apps."""

from __future__ import annotations

from pathlib import Path

import typer
from rich.console import Console
from rich.table import Table

from lightsuite.cli.runner import (
    run_brain_pipeline,
    run_multires_pipeline,
    run_spinal_pipeline,
)
from lightsuite.cli.stages import (
    StageState,
    StageStatus,
    brain_stage_statuses,
    multires_stage_statuses,
    spinal_stage_statuses,
)
from lightsuite.config.loader import load_config, load_multires_config, load_spinal_config

console = Console()

_STATE_ICONS = {
    StageState.DONE: "[green]✓[/green]",
    StageState.PENDING: "[red]✗[/red]",
    StageState.OPTIONAL: "[yellow]○[/yellow]",
    StageState.SKIPPED: "[dim]—[/dim]",
}


def _print_stage_table(statuses: list[StageStatus], *, save_path: Path) -> None:
    table = Table(title=f"Pipeline stages ({save_path})")
    table.add_column("Stage", style="cyan")
    table.add_column("Status")
    table.add_column("Checkpoint / detail")
    for item in statuses:
        spec = item.stage
        label = spec.title
        if spec.manual:
            label += " [dim](GUI)[/dim]"
        if spec.optional:
            label += " [dim](optional)[/dim]"
        table.add_row(
            label,
            _STATE_ICONS.get(item.state, "?"),
            item.detail or spec.checkpoint_hint,
        )
    console.print(table)


def _common_run_options():
    """Shared Typer options for run commands (documented in each command)."""
    return {}


def register_brain_commands(app: typer.Typer) -> None:
    @app.command("stages")
    def brain_stages(
        config: str = typer.Option(..., "--config", "-c", help="Pipeline YAML config."),
    ) -> None:
        """Show pipeline stage checklist and checkpoint status."""
        cfg = load_config(config)
        statuses = brain_stage_statuses(cfg)
        _print_stage_table(statuses, save_path=cfg.sample.save_path.expanduser())

    @app.command("run")
    def brain_run(
        config: str = typer.Option(..., "--config", "-c", help="Pipeline YAML config."),
        from_stage: str | None = typer.Option(
            None,
            "--from",
            help="First stage to run (default: first in pipeline).",
        ),
        through_stage: str | None = typer.Option(
            None,
            "--through",
            help="Last stage to run (default: last in pipeline).",
        ),
        resume: bool = typer.Option(
            False,
            "--resume",
            help="Skip stages that already have checkpoint artifacts.",
        ),
        include_optional: bool = typer.Option(
            False,
            "--include-optional",
            help="Run optional stages (align-slices, import-annotations).",
        ),
        headless: bool = typer.Option(
            False,
            "--headless",
            help="Run GUI stages in headless/test mode (HPC/CI only).",
        ),
        force_preprocess: bool = typer.Option(
            False,
            "--force-preprocess",
            help="Re-run preprocess even when cached TIFFs match the config fingerprint.",
        ),
    ) -> None:
        """Run brain pipeline stages in order (orchestrates individual stage commands)."""
        run_brain_pipeline(
            config,
            from_stage=from_stage,
            through_stage=through_stage,
            resume=resume,
            include_optional=include_optional,
            headless=headless,
            force_preprocess=force_preprocess,
        )


def register_spinal_commands(app: typer.Typer) -> None:
    @app.command("stages")
    def spinal_stages(
        config: str = typer.Option(..., "--config", "-c", help="Pipeline YAML config."),
    ) -> None:
        """Show spinal pipeline stage checklist and checkpoint status."""
        cfg = load_spinal_config(config)
        statuses = spinal_stage_statuses(cfg)
        _print_stage_table(statuses, save_path=cfg.sample.save_path.expanduser())

    @app.command("run")
    def spinal_run(
        config: str = typer.Option(..., "--config", "-c", help="Pipeline YAML config."),
        from_stage: str | None = typer.Option(None, "--from", help="First stage to run."),
        through_stage: str | None = typer.Option(None, "--through", help="Last stage to run."),
        resume: bool = typer.Option(False, "--resume", help="Skip completed stages."),
        include_optional: bool = typer.Option(
            True,
            "--include-optional/--no-include-optional",
            help="Run optional import-annotations and region-stats when configured.",
        ),
        headless: bool = typer.Option(
            False,
            "--headless",
            help="Run GUI stages in headless mode (HPC/CI only).",
        ),
    ) -> None:
        """Run spinal pipeline stages in order."""
        run_spinal_pipeline(
            config,
            from_stage=from_stage,
            through_stage=through_stage,
            resume=resume,
            include_optional=include_optional,
            headless=headless,
        )


def register_multires_commands(app: typer.Typer) -> None:
    @app.command("stages")
    def multires_stages(
        config: str = typer.Option(..., "--config", "-c", help="Multires pipeline YAML config."),
    ) -> None:
        """Show multires pipeline stage checklist and checkpoint status."""
        cfg = load_multires_config(config)
        statuses = multires_stage_statuses(cfg)
        _print_stage_table(statuses, save_path=cfg.sample.save_path.expanduser())

    @app.command("run")
    def multires_run(
        config: str = typer.Option(..., "--config", "-c", help="Multires pipeline YAML config."),
        from_stage: str | None = typer.Option(None, "--from", help="First stage to run."),
        through_stage: str | None = typer.Option(None, "--through", help="Last stage to run."),
        resume: bool = typer.Option(False, "--resume", help="Skip completed stages."),
        include_optional: bool = typer.Option(
            True,
            "--include-optional/--no-include-optional",
            help="Run optional match-points and import-annotations when configured.",
        ),
        headless: bool = typer.Option(
            False,
            "--headless",
            help="Run GUI stages in headless mode (HPC/CI only).",
        ),
    ) -> None:
        """Run multires pipeline stages in order."""
        run_multires_pipeline(
            config,
            from_stage=from_stage,
            through_stage=through_stage,
            resume=resume,
            include_optional=include_optional,
            headless=headless,
        )
