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


@brain_app.command("export")
def brain_export(
    config: str = typer.Option(..., "--config", "-c", help="Pipeline YAML config."),
    write_csv: bool = typer.Option(
        None,
        "--write-csv/--no-write-csv",
        help="Write parcellation intensity CSVs (default: export.write_cells_csv).",
    ),
    save_volume: bool = typer.Option(
        None,
        "--save-volume/--no-save-volume",
        help="Save registered atlas-space volumes (default: export.save_registered_volume).",
    ),
) -> None:
    """Apply transforms and export registered volumes (generateRegisteredBrainVolumes.m)."""
    from lightsuite.config.loader import load_config
    from lightsuite.export.brain_export import export_registered_brain_volumes

    cfg = load_config(config)
    result = export_registered_brain_volumes(
        cfg,
        write_csv=write_csv,
        save_registered_volume=save_volume,
    )
    if result.registered_volumes:
        typer.echo(f"Registered volumes: {len(result.registered_volumes)} channel(s)")
    if result.parcellation_paths:
        typer.echo(f"Parcellation CSVs: {len(result.parcellation_paths)} channel(s)")


@brain_app.command("register")
def brain_register(
    config: str = typer.Option(..., "--config", "-c", help="Pipeline YAML config."),
    single_step: bool = typer.Option(
        False,
        "--single-step",
        help="Use a single-resolution B-spline schedule (faster, lower quality).",
    ),
) -> None:
    """Run elastix registration (multiobjRegistration.m)."""
    from lightsuite.config.loader import load_config
    from lightsuite.registration.brain_register import run_brain_registration

    cfg = load_config(config)
    path = run_brain_registration(cfg, use_multistep=not single_step)
    typer.echo(f"Transform parameters: {path}")


@brain_app.command("match-points")
def brain_match_points(
    config: str = typer.Option(..., "--config", "-c", help="Pipeline YAML config."),
    headless: bool = typer.Option(
        False,
        "--headless",
        help="Load data and write an empty session without opening napari (for tests).",
    ),
) -> None:
    """Interactive control-point matching (matchControlPoints_unified.m)."""
    from lightsuite.config.loader import load_config
    from lightsuite.gui.match_points_brain import run_brain_match_points

    cfg = load_config(config)
    path = run_brain_match_points(cfg, headless=headless)
    typer.echo(f"Control points session: {path}")


@brain_app.command("init-registration")
def brain_init_registration(
    config: str = typer.Option(..., "--config", "-c", help="Pipeline YAML config."),
) -> None:
    """Coarse-align sample to atlas (initializeRegistration.m)."""
    from lightsuite.config.loader import load_config
    from lightsuite.registration.init_brain import initialize_brain_registration

    cfg = load_config(config)
    checkpoint = initialize_brain_registration(cfg)
    n_pairs = len(checkpoint.autocpsample or [])
    typer.echo(f"Initial registration complete. Auto control point pairs: {n_pairs}")


@brain_app.command("preprocess")
def brain_preprocess(
    config: str = typer.Option(..., "--config", "-c", help="Pipeline YAML config."),
) -> None:
    """Downsample sample volumes for registration (preprocessLightSheetVolume.m)."""
    from lightsuite.config.loader import load_config
    from lightsuite.preprocess.brain import preprocess_lightsheet_volume

    cfg = load_config(config)
    result = preprocess_lightsheet_volume(cfg)
    typer.echo(f"Primary registration volume: {result.checkpoint.regvolpath}")


def run() -> None:
    app()


if __name__ == "__main__":
    run()
