"""LightSuite Typer CLI entry point."""

from __future__ import annotations

from pathlib import Path

import typer

from lightsuite import __version__
from lightsuite.cli.doctor import doctor_command

app = typer.Typer(
    name="lightsuite",
    help="LightSuite — mouse lightsheet and histology atlas registration.",
    no_args_is_help=True,
)
brain_app = typer.Typer(help="Brain lightsheet pipeline stages.")
spinal_app = typer.Typer(help="Spinal cord lightsheet pipeline stages.")
mesospim_app = typer.Typer(help="mesoSPIM overview ↔ ROI registration.")
analysis_app = typer.Typer(help="Post-registration analysis (region stats, cell counts).")
app.add_typer(brain_app, name="brain")
app.add_typer(spinal_app, name="spinal")
app.add_typer(mesospim_app, name="mesospim")
app.add_typer(analysis_app, name="analysis")


def _version_callback(value: bool) -> None:
    if value:
        typer.echo(f"lightsuite {__version__}")
        raise typer.Exit()


@app.callback()
def main_callback(
    version: bool = typer.Option(
        False,
        "--version",
        "-V",
        help="Show package version and exit.",
        callback=_version_callback,
        is_eager=True,
    ),
) -> None:
    """LightSuite CLI."""


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


@brain_app.command("check-orientation")
def brain_check_orientation(
    config: str = typer.Option(..., "--config", "-c", help="Pipeline YAML config."),
    headless: bool = typer.Option(
        False,
        "--headless",
        help="Update registration.orientation in the config YAML without opening Napari.",
    ),
) -> None:
    """Interactive brain axis/orientation checker (getBrainOrientation.m)."""
    from lightsuite.config.loader import load_config
    from lightsuite.gui.orientation_brain import run_brain_orientation_check

    config_path = Path(config).expanduser().resolve()
    cfg = load_config(config_path)
    path = run_brain_orientation_check(cfg, config_path, headless=headless)
    typer.echo(f"Brain orientation saved to {path}")


@brain_app.command("import-annotations")
def brain_import_annotations(
    config: str = typer.Option(..., "--config", "-c", help="Pipeline YAML config."),
    write_csv: bool = typer.Option(
        None,
        "--write-csv/--no-write-csv",
        help="Write atlas coordinate CSVs for point imports.",
    ),
) -> None:
    """Register native sample-space annotations into atlas space."""
    from lightsuite.config.loader import load_config
    from lightsuite.import_.brain_import import run_brain_import_annotations

    cfg = load_config(config)
    results = run_brain_import_annotations(cfg, write_csv=write_csv)
    for item in results:
        typer.echo(f"{item.label}: {item.kind} — {item.n_atlas} atlas features")


@brain_app.command("inspect-imports")
def brain_inspect_imports(
    config: str = typer.Option(..., "--config", "-c", help="Pipeline YAML config."),
    headless: bool = typer.Option(
        False,
        "--headless",
        help="Validate inspect inputs without opening Napari.",
    ),
) -> None:
    """Napari QC: registered channels, atlas, and imported points/masks in atlas space."""
    from lightsuite.config.loader import load_config
    from lightsuite.gui.inspect_brain_imports import run_brain_inspect_imports

    cfg = load_config(config)
    paths = run_brain_inspect_imports(cfg, headless=headless)
    typer.echo(f"volume_registered: {paths.volume_registered_dir}")
    if paths.registered_channels:
        typer.echo(f"  channels: {sorted(paths.registered_channels)}")
    if paths.point_npz_paths:
        typer.echo(f"  point layers: {list(paths.point_npz_paths)}")
    if paths.mask_paths:
        typer.echo(f"  mask layers: {list(paths.mask_paths)}")


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
    run_brain_registration(cfg, use_multistep=not single_step)


@brain_app.command("refine-auto-points")
def brain_refine_auto_points(
    config: str = typer.Option(..., "--config", "-c", help="Pipeline YAML config."),
    bootstrap_correspondence: bool = typer.Option(
        False,
        "--bootstrap-correspondence",
        help="Run headless align-slices first if slice_correspondence.json is missing.",
    ),
    tolerance_vox: float | None = typer.Option(
        None,
        "--tolerance-vox",
        help="Override registration.ap_pair_tolerance_vox.",
    ),
    force: bool = typer.Option(
        False,
        "--force",
        help="Re-run even if auto points were already refined.",
    ),
) -> None:
    """Filter init-registration auto pairs using AP slice correspondence."""
    from lightsuite.config.loader import load_config
    from lightsuite.registration.refine_auto_points import refine_brain_auto_points

    cfg = load_config(config)
    path = refine_brain_auto_points(
        cfg,
        bootstrap_correspondence=bootstrap_correspondence,
        tolerance_vox=tolerance_vox,
        force=force,
    )
    typer.echo(f"Updated checkpoint: {path}")


@brain_app.command("align-slices")
def brain_align_slices(
    config: str = typer.Option(..., "--config", "-c", help="Pipeline YAML config."),
    headless: bool = typer.Option(
        False,
        "--headless",
        help="Auto-estimate slice correspondence without opening Napari (for tests).",
    ),
) -> None:
    """Interactive sample-to-atlas slice alignment before control-point matching."""
    from lightsuite.config.loader import load_config
    from lightsuite.gui.align_slices_brain import run_brain_align_slices

    cfg = load_config(config)
    path = run_brain_align_slices(cfg, headless=headless)
    typer.echo(f"Slice correspondence: {path}")


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
    initialize_brain_registration(cfg)


@brain_app.command("preprocess")
def brain_preprocess(
    config: str = typer.Option(..., "--config", "-c", help="Pipeline YAML config."),
    force: bool = typer.Option(
        False,
        "--force",
        help="Re-run downsampling even when cached registration TIFFs match the current config.",
    ),
) -> None:
    """Downsample sample volumes for registration (preprocessLightSheetVolume.m)."""
    from lightsuite.config.loader import load_config
    from lightsuite.preprocess.brain import preprocess_lightsheet_volume

    cfg = load_config(config)
    result = preprocess_lightsheet_volume(cfg, force=force)
    typer.echo(f"Primary registration volume: {result.checkpoint.regvolpath}")


@analysis_app.command("region-stats")
def analysis_region_stats(
    config: str = typer.Option(..., "--config", "-c", help="Pipeline YAML config."),
    count_points: bool = typer.Option(
        None,
        "--count-points/--no-count-points",
        help="Bin imported atlas-space points into cell counts (default: analysis.count_points).",
    ),
) -> None:
    """Assemble a tidy region_stats.csv (intensity + cell counts) with region names."""
    from lightsuite.analysis.runner import run_region_stats
    from lightsuite.config.loader import load_config

    cfg = load_config(config)
    result = run_region_stats(cfg, count_points=count_points)
    if result.combined_path is not None:
        typer.echo(f"Region stats: {result.combined_path} ({result.n_rows} rows)")
    else:
        typer.echo("No region stats produced (run 'lightsuite brain export' first).")


@analysis_app.command("validate-cohort")
def analysis_validate_cohort(
    config: str = typer.Option(..., "--config", "-c", help="Cohort YAML config."),
) -> None:
    """Load and validate a cross-subject cohort configuration file."""
    from lightsuite.config.loader import load_cohort_config

    cfg = load_cohort_config(config)
    typer.echo(
        f"Cohort valid: '{cfg.name}' with {len(cfg.samples)} sample(s), "
        f"output → {cfg.output_dir}"
    )


@analysis_app.command("group-stats")
def analysis_group_stats(
    config: str = typer.Option(..., "--config", "-c", help="Cohort YAML config."),
) -> None:
    """Cross-subject group summaries and optional pairwise comparisons."""
    from lightsuite.analysis.cohort_runner import run_cohort_group_analysis
    from lightsuite.config.loader import load_cohort_config

    cfg = load_cohort_config(config)
    result = run_cohort_group_analysis(cfg)
    for path in result.written_paths:
        typer.echo(f"  {path.name}")
    typer.echo(f"Output directory: {result.output_dir}")


def _resolve_plot_input(
    *,
    input_path: str | None,
    brain_config: str | None,
) -> "Path":
    from pathlib import Path

    from lightsuite.analysis.viz.io import resolve_region_stats_from_config

    if input_path and brain_config:
        raise typer.BadParameter("Use only one of --input or --config.")
    if brain_config:
        return resolve_region_stats_from_config(brain_config)
    if input_path:
        return Path(input_path).expanduser().resolve()
    raise typer.BadParameter("Provide --input or --config.")


@analysis_app.command("plot-division-bars")
def analysis_plot_division_bars(
    input_path: str | None = typer.Option(None, "--input", "-i", help="region_stats or intensities CSV."),
    brain_config: str | None = typer.Option(
        None, "--config", "-c", help="Brain YAML; uses volume_registered/region_stats.csv."
    ),
    output: str = typer.Option(..., "--output", "-o", help="Output PNG path."),
    channel: int | None = typer.Option(1, "--channel", help="Channel id (tidy tables only)."),
    metric: str = typer.Option("median_intensity", "--metric", help="Metric to plot."),
    title: str | None = typer.Option(None, "--title", help="Figure title."),
    exclude_division: list[str] | None = typer.Option(
        None, "--exclude-division", help="Division names to omit (repeatable)."
    ),
    dpi: int = typer.Option(200, "--dpi", help="Figure DPI."),
) -> None:
    """Bar plot: left vs right mean per Allen division."""
    from lightsuite.analysis.viz.io import load_region_plot_table
    from lightsuite.analysis.viz.plots import plot_division_bars

    csv_path = _resolve_plot_input(input_path=input_path, brain_config=brain_config)
    table = load_region_plot_table(csv_path, channel=channel, metric=metric)
    out = Path(output).expanduser()
    plot_division_bars(
        table,
        title=title,
        output_path=out,
        dpi=dpi,
        exclude_divisions=exclude_division,
    )
    typer.echo(f"Saved: {out.resolve()}")


@analysis_app.command("plot-lr-scatter")
def analysis_plot_lr_scatter(
    input_path: str | None = typer.Option(None, "--input", "-i", help="region_stats or intensities CSV."),
    brain_config: str | None = typer.Option(
        None, "--config", "-c", help="Brain YAML; uses volume_registered/region_stats.csv."
    ),
    output: str = typer.Option(..., "--output", "-o", help="Output PNG path."),
    channel: int | None = typer.Option(1, "--channel", help="Channel id (tidy tables only)."),
    metric: str = typer.Option("median_intensity", "--metric", help="Metric to plot."),
    title: str | None = typer.Option(None, "--title", help="Figure title."),
    keep_division: list[str] | None = typer.Option(
        None, "--keep-division", help="Only these divisions (repeatable)."
    ),
    exclude_division: list[str] | None = typer.Option(
        None, "--exclude-division", help="Division names to omit (repeatable)."
    ),
    axis_min: float | None = typer.Option(None, "--axis-min", help="Fixed scatter axis minimum."),
    axis_max: float | None = typer.Option(None, "--axis-max", help="Fixed scatter axis maximum."),
    dpi: int = typer.Option(200, "--dpi", help="Figure DPI."),
) -> None:
    """Scatter plot: left vs right per region, coloured by division."""
    from lightsuite.analysis.viz.io import load_region_plot_table
    from lightsuite.analysis.viz.plots import plot_lr_scatter

    csv_path = _resolve_plot_input(input_path=input_path, brain_config=brain_config)
    table = load_region_plot_table(csv_path, channel=channel, metric=metric)
    out = Path(output).expanduser()
    plot_lr_scatter(
        table,
        title=title,
        output_path=out,
        dpi=dpi,
        keep_divisions=keep_division,
        exclude_divisions=exclude_division,
        axis_min=axis_min,
        axis_max=axis_max,
    )
    typer.echo(f"Saved: {out.resolve()}")


@analysis_app.command("plot-group-division")
def analysis_plot_group_division(
    input_path: str = typer.Option(
        ...,
        "--input",
        "-i",
        help="group_summary_by_division.csv from cohort group-stats.",
    ),
    output: str = typer.Option(..., "--output", "-o", help="Output PNG path."),
    title: str | None = typer.Option(None, "--title", help="Figure title."),
    dpi: int = typer.Option(200, "--dpi", help="Figure DPI."),
) -> None:
    """Bar plot: group means by division (cohort output)."""
    import pandas as pd

    from lightsuite.analysis.viz.cohort_plots import plot_group_division_bars

    summary = pd.read_csv(input_path)
    out = Path(output).expanduser()
    plot_group_division_bars(summary, title=title, output_path=out, dpi=dpi)
    typer.echo(f"Saved: {out.resolve()}")


@analysis_app.command("build-division-map")
def analysis_build_division_map(
    config: str = typer.Option(..., "--config", "-c", help="Brain pipeline YAML config."),
    force: bool = typer.Option(False, "--force", help="Rebuild cached division labels."),
) -> None:
    """Build or refresh atlas-side division label volume + legend CSV."""
    from lightsuite.analysis.division_map import ensure_division_map
    from lightsuite.atlas.registry import resolve_brain_atlas_with_config
    from lightsuite.config.loader import load_config
    from lightsuite.export.brain_export import _load_transform_params

    cfg = load_config(config)
    transform_params = _load_transform_params(cfg.sample.save_path.expanduser())
    atlas = resolve_brain_atlas_with_config(transform_params.brain_atlas, cfg.atlas)
    result = ensure_division_map(atlas, force=force)
    typer.echo(f"Division labels: {result.paths.labels_tiff}")
    typer.echo(f"Legend: {result.paths.legend_csv} ({len(result.legend)} divisions)")


@analysis_app.command("view-divisions")
def analysis_view_divisions(
    config: str = typer.Option(..., "--config", "-c", help="Brain pipeline YAML config."),
    headless: bool = typer.Option(
        False,
        "--headless",
        help="Validate inputs and build division map without opening Napari.",
    ),
    stride: int = typer.Option(
        1,
        "--stride",
        min=1,
        help="Load every Nth voxel along each axis (faster interaction on large volumes).",
    ),
    force: bool = typer.Option(
        False,
        "--force",
        help="Rebuild cached atlas division labels before opening the viewer.",
    ),
) -> None:
    """Napari: toggle atlas divisions on registered channel volumes."""
    from lightsuite.config.loader import load_config
    from lightsuite.gui.view_divisions_brain import run_brain_division_viewer

    cfg = load_config(config)
    paths = run_brain_division_viewer(
        cfg,
        headless=headless,
        stride=stride,
        force_division_rebuild=force,
    )
    typer.echo(f"Division labels: {paths.division_labels}")
    typer.echo(f"Channels: {sorted(paths.channel_paths)}")


@analysis_app.command("registration-qc")
def analysis_registration_qc(
    config: str = typer.Option(..., "--config", "-c", help="Brain pipeline YAML config."),
    channel: int = typer.Option(1, "--channel", help="Registered channel to score."),
    threshold: float | None = typer.Option(
        None,
        "--threshold",
        "-t",
        help="Intensity threshold (above = signal). Auto-estimated from volume if omitted.",
    ),
    inspect: bool = typer.Option(
        False,
        "--inspect",
        help="Open Napari first to visually validate the threshold (requires gui extra).",
    ),
    sweep: bool = typer.Option(
        False,
        "--sweep",
        help="Also run a threshold sensitivity sweep and save CSV + plot.",
    ),
    sweep_points: int = typer.Option(9, "--sweep-points", min=2, help="Points in sweep."),
    headless: bool = typer.Option(
        False,
        "--headless",
        help="Skip Napari inspect; validate inputs and write score only.",
    ),
    force: bool = typer.Option(
        False,
        "--force",
        help="Rebuild cached atlas division labels before scoring.",
    ),
) -> None:
    """Naive registration QC: fraction of signal voxels in unassigned divisions."""
    from lightsuite.analysis.registration_qc_runner import run_registration_qc
    from lightsuite.config.loader import load_config

    cfg = load_config(config)
    result = run_registration_qc(
        cfg,
        channel=channel,
        threshold=threshold,
        inspect=inspect,
        sweep=sweep,
        sweep_points=sweep_points,
        headless=headless,
        force_division_rebuild=force,
    )
    if result.score_csv is not None:
        typer.echo(
            f"Unassigned fraction: {result.score['naive_unassigned_percent']:.2f}% "
            f"(threshold={result.threshold:g})"
        )
    if result.sweep_plot is not None:
        typer.echo(f"Sweep plot: {result.sweep_plot}")


@mesospim_app.command("validate-config")
def mesospim_validate_config(
    config: str = typer.Option(..., "--config", "-c", help="mesoSPIM pipeline YAML config."),
) -> None:
    """Load and validate a mesoSPIM overview / ROI YAML config."""
    from lightsuite.config.loader import load_mesospim_config

    cfg = load_mesospim_config(config)
    typer.echo(
        f"Config valid: {cfg.sample.name} "
        f"({cfg.mesospim.overview.path.name} → {cfg.mesospim.roi.path.name})"
    )


@mesospim_app.command("check-geometry")
def mesospim_check_geometry(
    config: str = typer.Option(..., "--config", "-c", help="mesoSPIM pipeline YAML config."),
) -> None:
    """Validate FOV overlap and write geometry QA artifacts."""
    from lightsuite.config.loader import load_mesospim_config
    from lightsuite.mesospim.runner import check_mesospim_geometry

    check_mesospim_geometry(load_mesospim_config(config))


@mesospim_app.command("register")
def mesospim_register(
    config: str = typer.Option(..., "--config", "-c", help="mesoSPIM pipeline YAML config."),
) -> None:
    """Register ROI stack to overview using metadata geometry and elastix."""
    from lightsuite.config.loader import load_mesospim_config
    from lightsuite.mesospim.runner import run_mesospim_registration

    run_mesospim_registration(load_mesospim_config(config))


@mesospim_app.command("match-points")
def mesospim_match_points(
    config: str = typer.Option(..., "--config", "-c", help="mesoSPIM pipeline YAML config."),
    headless: bool = typer.Option(
        False,
        "--headless",
        help="Create an empty landmark session without opening Napari (for tests).",
    ),
) -> None:
    """Interactive overview / ROI landmark placement (Napari)."""
    from lightsuite.config.loader import load_mesospim_config
    from lightsuite.gui.match_points_mesospim import run_mesospim_match_points

    cfg = load_mesospim_config(config)
    path = run_mesospim_match_points(cfg, headless=headless)
    typer.echo(f"Landmark session: {path}")


@mesospim_app.command("inspect")
def mesospim_inspect(
    config: str = typer.Option(..., "--config", "-c", help="mesoSPIM pipeline YAML config."),
    headless: bool = typer.Option(
        False,
        "--headless",
        help="Validate inspect inputs without opening Napari.",
    ),
) -> None:
    """Compare the 1× overview and registered ROI embedded in the full overview canvas."""
    from lightsuite.config.loader import load_mesospim_config
    from lightsuite.gui.inspect_mesospim import run_mesospim_inspect

    cfg = load_mesospim_config(config)
    paths =     run_mesospim_inspect(cfg, headless=headless)
    typer.echo(f"Overview: {paths.overview_path}")
    typer.echo(f"Registered canvas: {paths.registered_full_overview_path}")


@spinal_app.command("validate-config")
def spinal_validate_config(
    config: str = typer.Option(..., "--config", "-c", help="Spinal cord pipeline YAML config."),
) -> None:
    """Load and validate a spinal cord pipeline configuration file."""
    from lightsuite.config.loader import load_spinal_config

    cfg = load_spinal_config(config)
    typer.echo(
        f"Config valid for sample '{cfg.sample.name}' "
        f"({cfg.sample.source.tiff_type.value})."
    )


@spinal_app.command("preprocess")
def spinal_preprocess(
    config: str = typer.Option(..., "--config", "-c", help="Spinal cord pipeline YAML config."),
) -> None:
    """Prepare cord sample and atlas (prepareCordSampleForRegistration.m)."""
    from lightsuite.config.loader import load_spinal_config
    from lightsuite.preprocess.cord import preprocess_spinal_cord_sample

    cfg = load_spinal_config(config)
    result = preprocess_spinal_cord_sample(cfg)
    typer.echo(f"Wrote checkpoint: {result.checkpoint.lsfolder}/regopts.json")


@spinal_app.command("straighten")
def spinal_straighten(
    config: str = typer.Option(..., "--config", "-c", help="Spinal cord pipeline YAML config."),
    headless: bool = typer.Option(
        False,
        "--headless",
        help="Write default alignment without opening Napari (for tests).",
    ),
) -> None:
    """Interactive cord straightening GUI (spinal_cord_aligner.m)."""
    from lightsuite.config.loader import load_spinal_config
    from lightsuite.gui.straighten_cord import run_spinal_straighten

    cfg = load_spinal_config(config)
    path = run_spinal_straighten(cfg, headless=headless)
    typer.echo(f"Alignment checkpoint: {path}")


@spinal_app.command("init-registration")
def spinal_init_registration(
    config: str = typer.Option(..., "--config", "-c", help="Spinal cord pipeline YAML config."),
) -> None:
    """Apply straightening and coarse affine registration (initializeCordRegistration.m)."""
    from lightsuite.config.loader import load_spinal_config
    from lightsuite.registration.init_cord import initialize_cord_registration

    cfg = load_spinal_config(config)
    initialize_cord_registration(cfg)


@spinal_app.command("match-points")
def spinal_match_points(
    config: str = typer.Option(..., "--config", "-c", help="Spinal cord pipeline YAML config."),
    headless: bool = typer.Option(
        False,
        "--headless",
        help="Write an empty control-point session without Napari.",
    ),
) -> None:
    """Interactive control-point matching (matchControlPointsSpine.m)."""
    from lightsuite.config.loader import load_spinal_config
    from lightsuite.gui.match_points_cord import run_spinal_match_points

    cfg = load_spinal_config(config)
    path = run_spinal_match_points(cfg, headless=headless)
    typer.echo(f"Control points session: {path}")


@spinal_app.command("register")
def spinal_register(
    config: str = typer.Option(..., "--config", "-c", help="Spinal cord pipeline YAML config."),
) -> None:
    """Run B-spline registration (multiobjCordRegistration.m)."""
    from lightsuite.config.loader import load_spinal_config
    from lightsuite.registration.cord_register import run_spinal_registration

    cfg = load_spinal_config(config)
    path = run_spinal_registration(cfg)
    typer.echo(f"Transform params: {path}")


@spinal_app.command("export")
def spinal_export(
    config: str = typer.Option(..., "--config", "-c", help="Spinal cord pipeline YAML config."),
) -> None:
    """Export registered cord volumes (generateRegisteredCordVolume.m)."""
    from lightsuite.config.loader import load_spinal_config
    from lightsuite.export.cord_export import export_registered_cord_volumes

    cfg = load_spinal_config(config)
    result = export_registered_cord_volumes(cfg)
    typer.echo(f"Registered volumes in {result.output_dir} ({len(result.channel_paths)} channels)")


@spinal_app.command("view")
def spinal_view(
    config: str = typer.Option(..., "--config", "-c", help="Spinal cord pipeline YAML config."),
    headless: bool = typer.Option(
        False,
        "--headless",
        help="Validate view inputs without opening Napari.",
    ),
    recompute_annotation: bool = typer.Option(
        False,
        "--recompute-annotation",
        help="Rebuild annotation_registered.tiff from transform_params.json.",
    ),
) -> None:
    """Open Napari with registered sample channel(s) and warped atlas annotation."""
    from lightsuite.config.loader import load_spinal_config
    from lightsuite.gui.view_registered_cord import run_spinal_registered_view

    cfg = load_spinal_config(config)
    paths = run_spinal_registered_view(
        cfg,
        headless=headless,
        recompute_annotation=recompute_annotation,
    )
    typer.echo(f"Registered data: {paths.volume_registered_dir}")
    if paths.annotation_path.is_file():
        typer.echo(f"Annotation: {paths.annotation_path.name}")


@spinal_app.command("validate-parity")
def spinal_validate_parity(
    fixture_root: str = typer.Option(
        "tests/fixtures/spinal_cord",
        "--fixture-root",
        help="Directory containing parity reference fixtures.",
    ),
) -> None:
    """Run automated MATLAB-parity checks for MVP stages (optimizer reference)."""
    from lightsuite.validation.spinal_parity import run_mvp_parity_checks

    report = run_mvp_parity_checks(Path(fixture_root).expanduser().resolve())
    for message in report.messages:
        typer.echo(message)
    if not report.passed:
        raise typer.Exit(code=1)


def run() -> None:
    app()


if __name__ == "__main__":
    run()
