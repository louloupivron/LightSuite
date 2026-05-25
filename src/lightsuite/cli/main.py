"""LightSuite Typer CLI entry point."""

from __future__ import annotations

import typer

from lightsuite.cli.doctor import doctor_command

app = typer.Typer(
    name="lightsuite",
    help="LightSuite — mouse lightsheet and histology atlas registration.",
    no_args_is_help=True,
)
brain_app = typer.Typer(help="Brain lightsheet pipeline stages.")
app.add_typer(brain_app, name="brain")


@app.command("doctor")
def doctor(
    config: str | None = typer.Option(
        None,
        "--config",
        "-c",
        help="Optional YAML config to validate paths (scratch, atlas_dir, etc.).",
    ),
    strict: bool = typer.Option(
        False,
        "--strict",
        help="Treat optional atlas warnings as errors.",
    ),
) -> None:
    """Verify Python environment, Elastix, atlases, and optional GPU."""
    doctor_command(config_path=config, strict=strict)


@brain_app.command("validate-config")
def validate_config(
    config: str = typer.Option(..., "--config", "-c", help="Pipeline YAML config."),
) -> None:
    """Load and validate a brain pipeline configuration file."""
    from lightsuite.config.loader import load_config

    cfg = load_config(config)
    typer.echo(f"Config valid for sample '{cfg.sample.name}' ({cfg.sample.source.format.value}).")


def run() -> None:
    app()


if __name__ == "__main__":
    run()
