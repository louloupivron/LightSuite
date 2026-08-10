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
multires_app = typer.Typer(help="Manifest-driven overview ↔ ROI multiresolution registration.")
analysis_app = typer.Typer(help="Post-registration analysis (region stats, cell counts).")
app.add_typer(brain_app, name="brain")
app.add_typer(spinal_app, name="spinal")
app.add_typer(mesospim_app, name="mesospim")
app.add_typer(multires_app, name="multires")
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


@brain_app.command("probe-content-box")
def brain_probe_content_box(
    config: str = typer.Option(..., "--config", "-c", help="Pipeline YAML config."),
    target: str = typer.Option(
        "atlas",
        "--target",
        "-t",
        help="Crop target: atlas (native BrainGlobe grid) or sample (registration TIFF).",
    ),
    headless: bool = typer.Option(
        False,
        "--headless",
        help="Auto-detect and print the crop box without opening Napari.",
    ),
    write_config: bool = typer.Option(
        False,
        "--write-config",
        help="Write manual crop mode and box into the YAML (works with --headless).",
    ),
) -> None:
    """Probe or interactively pick atlas/sample content crop boxes."""
    from lightsuite.config.loader import load_config, save_content_box_to_config
    from lightsuite.gui.content_box_brain import run_content_box_picker
    from lightsuite.registration.content_probe import (
        ContentProbeTarget,
        format_content_box_report,
        load_content_probe_data,
    )

    config_path = Path(config).expanduser().resolve()
    cfg = load_config(config_path)
    probe_target = ContentProbeTarget(target)
    if headless:
        data = load_content_probe_data(cfg, target=probe_target)
        typer.echo(format_content_box_report(data))
        if write_config:
            path = save_content_box_to_config(
                config_path,
                target=probe_target.value,
                box=data.box.to_manual_list(),
            )
            typer.echo(f"Updated config: {path}")
        return

    box = run_content_box_picker(
        cfg,
        config_path,
        target=probe_target,
        write_config=write_config,
    )
    typer.echo(f"Selected box: {box.to_manual_list()}")


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


@brain_app.command("convert-fiji-points")
def brain_convert_fiji_points(
    source: str = typer.Option(..., "--source", "-s", help="FIJI point-tool Results.csv export."),
    output: str = typer.Option(..., "--output", "-o", help="Output points.csv path."),
    voxel_um: str | None = typer.Option(
        None,
        "--voxel-um",
        help=(
            "ImageJ XY calibration as comma-separated x,y,z in the same units as the "
            "Results table (typically µm). Omit when X/Y are already in pixels."
        ),
    ),
) -> None:
    """Convert FIJI point-tool Results.csv to LightSuite points.csv."""
    from pathlib import Path

    from lightsuite.import_.fiji import convert_fiji_points_to_csv

    source_path = Path(source).expanduser()
    if not source_path.is_file():
        msg = f"Source CSV not found: {source_path}"
        raise typer.BadParameter(msg)

    voxel: list[float] | None = None
    if voxel_um is not None:
        parts = [float(v.strip()) for v in voxel_um.split(",")]
        if len(parts) != 3:
            msg = "--voxel-um must have exactly three values: x,y,z"
            raise typer.BadParameter(msg)
        voxel = parts

    output_path = Path(output).expanduser()
    n = convert_fiji_points_to_csv(source_path, output_path, voxel_um=voxel)
    typer.echo(f"Wrote {n} points → {output_path}")


@brain_app.command("import-matlab-control-points")
def brain_import_matlab_control_points(
    mat: str = typer.Option(
        ...,
        "--mat",
        "-m",
        help="MATLAB atlas2histology_tform.mat (or other *tform.mat with control points).",
    ),
    output: str | None = typer.Option(
        None,
        "--output",
        "-o",
        help="Output JSON path (default: same folder as --mat, atlas2histology_tform.json).",
    ),
    save_path: str | None = typer.Option(
        None,
        "--save-path",
        help="Sample save_path from config; writes atlas2histology_tform.json there.",
    ),
) -> None:
    """Import MATLAB match-points session for Python register (parity with multiobjRegistration.m)."""
    from lightsuite.gui.control_points import ControlPointSession
    from lightsuite.import_.matlab_control_points import import_matlab_control_points

    mat_path = Path(mat).expanduser()
    if not mat_path.is_file():
        msg = f"MATLAB file not found: {mat_path}"
        raise typer.BadParameter(msg)

    if save_path is not None:
        out_path = Path(save_path).expanduser() / "atlas2histology_tform.json"
    elif output is not None:
        out_path = Path(output).expanduser()
    else:
        out_path = mat_path.with_name("atlas2histology_tform.json")

    written = import_matlab_control_points(mat_path, output_json=out_path)
    matched, total_sample, total_atlas = ControlPointSession.load(written).point_counts()
    typer.echo(
        f"Wrote {written} — {matched} matched pairs "
        f"({total_sample} sample / {total_atlas} atlas points across slices)."
    )


@brain_app.command("migrate-control-point-session")
def brain_migrate_control_point_session(
    json_path: str | None = typer.Option(
        None,
        "--json",
        "-j",
        help="atlas2histology_tform.json to migrate.",
    ),
    save_path: str | None = typer.Option(
        None,
        "--save-path",
        help="Sample save_path; migrates atlas2histology_tform.json there.",
    ),
    dry_run: bool = typer.Option(
        False,
        "--dry-run",
        help="Report how many points would be updated without writing.",
    ),
) -> None:
    """Fix in-plane axis swap in Napari-authored control-point sessions (pre-v2 schema).

    Do not run on sessions imported from MATLAB (``point_coord_source: matlab``).
    Re-import those from ``.mat`` if needed.
    """
    from lightsuite.gui.control_points import (
        COORD_SCHEMA_VERSION,
        POINT_COORD_SOURCE_MATLAB,
        ControlPointSession,
        default_session_path,
        migrate_napari_transposed_control_points,
        session_needs_napari_transpose_migration,
    )

    if save_path is not None:
        session_path = default_session_path(Path(save_path).expanduser())
    elif json_path is not None:
        session_path = Path(json_path).expanduser()
    else:
        msg = "Provide --json or --save-path."
        raise typer.BadParameter(msg)

    if not session_path.is_file():
        msg = f"Session not found: {session_path}"
        raise typer.BadParameter(msg)

    session = ControlPointSession.load(session_path)
    if session.point_coord_source == POINT_COORD_SOURCE_MATLAB:
        typer.echo(
            f"Skipping {session_path}: imported from MATLAB "
            f"(point_coord_source={POINT_COORD_SOURCE_MATLAB!r})."
        )
        raise typer.Exit(0)

    if not session_needs_napari_transpose_migration(session):
        typer.echo(
            f"No migration needed for {session_path} "
            f"(coord_schema_version={session.coord_schema_version})."
        )
        raise typer.Exit(0)

    if session.chooselist is not None:
        chooselist = session.chooselist
    else:
        typer.echo(
            "Warning: session has no chooselist; cannot migrate in-plane axes. "
            "Open match-points once to populate chooselist, then re-run."
        )
        raise typer.Exit(1)

    if dry_run:
        n_hist = sum(len(s) for s in session.histology_control_points)
        n_atlas = sum(len(s) for s in session.atlas_control_points)
        typer.echo(
            f"Would migrate {n_hist + n_atlas} point rows in {session_path} "
            f"(schema v{session.coord_schema_version} -> v{COORD_SCHEMA_VERSION})."
        )
        raise typer.Exit(0)

    updated = migrate_napari_transposed_control_points(session, chooselist)
    session.save(session_path)
    typer.echo(
        f"Migrated {updated} point rows in {session_path} "
        f"(coord_schema_version={COORD_SCHEMA_VERSION})."
    )


@brain_app.command("export-matlab-control-points")
def brain_export_matlab_control_points(
    json_path: str | None = typer.Option(
        None,
        "--json",
        "-j",
        help="atlas2histology_tform.json from Napari match-points.",
    ),
    save_path: str | None = typer.Option(
        None,
        "--save-path",
        help="Sample save_path; reads atlas2histology_tform.json there.",
    ),
    output: str | None = typer.Option(
        None,
        "--output",
        "-o",
        help="Output .mat path (default: atlas2histology_tform.mat beside the JSON).",
    ),
) -> None:
    """Export Napari control points for MATLAB register (multiobjRegistration.m)."""
    from lightsuite.gui.control_points import ControlPointSession, default_session_path
    from lightsuite.import_.matlab_control_points import export_matlab_control_points

    if save_path is not None:
        session_path = default_session_path(Path(save_path).expanduser())
    elif json_path is not None:
        session_path = Path(json_path).expanduser()
    else:
        msg = "Provide --json or --save-path."
        raise typer.BadParameter(msg)

    out_path = Path(output).expanduser() if output is not None else None
    written = export_matlab_control_points(session_path, output_mat=out_path)
    matched, total_sample, total_atlas = ControlPointSession.load(session_path).point_counts()
    typer.echo(
        f"Wrote {written} — {matched} matched pairs "
        f"({total_sample} sample / {total_atlas} atlas points across slices)."
    )


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
    space: str | None = typer.Option(
        None,
        "--space",
        help="Output space: atlas, sample, or both (default: export.spaces in YAML).",
    ),
) -> None:
    """Apply transforms and export registered volumes (generateRegisteredBrainVolumes.m)."""
    from lightsuite.cli.spaces import parse_spaces_option
    from lightsuite.config.loader import load_config
    from lightsuite.export.brain_export import export_registered_brain_volumes

    cfg = load_config(config)
    result = export_registered_brain_volumes(
        cfg,
        write_csv=write_csv,
        save_registered_volume=save_volume,
        spaces=parse_spaces_option(space),
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
    space: str | None = typer.Option(
        None,
        "--space",
        help="Stats space: atlas, sample, or both (default: analysis.stats_spaces).",
    ),
) -> None:
    """Assemble a tidy region_stats.csv (intensity + cell counts) with region names."""
    from lightsuite.analysis.runner import run_region_stats
    from lightsuite.cli.spaces import parse_spaces_option
    from lightsuite.config.loader import load_config

    cfg = load_config(config)
    result = run_region_stats(
        cfg,
        count_points=count_points,
        stats_spaces=parse_spaces_option(space),
    )
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
    channel: str = typer.Option(
        "1",
        "--channel",
        help="Imaging channel (1, 2, …) or import label such as imaris_488_cells (tidy tables only).",
    ),
    metric: str = typer.Option("median_intensity", "--metric", help="Metric to plot."),
    aggregate: str = typer.Option(
        "auto",
        "--aggregate",
        help="Division rollup: auto (sum for cell_count, mean otherwise), sum, or mean.",
    ),
    title: str | None = typer.Option(None, "--title", help="Figure title."),
    exclude_division: list[str] | None = typer.Option(
        None, "--exclude-division", help="Division names to omit (repeatable)."
    ),
    dpi: int = typer.Option(200, "--dpi", help="Figure DPI."),
) -> None:
    """Bar plot: left vs right per Allen division."""
    from lightsuite.analysis.viz.io import (
        default_division_aggregate,
        load_region_plot_table,
        parse_plot_channel,
    )
    from lightsuite.analysis.viz.plots import plot_division_bars

    agg_mode = aggregate.lower().strip()
    if agg_mode == "auto":
        agg_mode = default_division_aggregate(metric)
    elif agg_mode not in ("sum", "mean"):
        raise typer.BadParameter("--aggregate must be auto, sum, or mean.")
    csv_path = _resolve_plot_input(input_path=input_path, brain_config=brain_config)
    table = load_region_plot_table(
        csv_path, channel=parse_plot_channel(channel), metric=metric
    )
    out = Path(output).expanduser()
    plot_division_bars(
        table,
        title=title,
        output_path=out,
        dpi=dpi,
        exclude_divisions=exclude_division,
        aggregate=agg_mode,  # type: ignore[arg-type]
    )
    typer.echo(f"Saved: {out.resolve()} (aggregate={agg_mode})")


@analysis_app.command("plot-lr-scatter")
def analysis_plot_lr_scatter(
    input_path: str | None = typer.Option(None, "--input", "-i", help="region_stats or intensities CSV."),
    brain_config: str | None = typer.Option(
        None, "--config", "-c", help="Brain YAML; uses volume_registered/region_stats.csv."
    ),
    output: str = typer.Option(..., "--output", "-o", help="Output PNG path."),
    channel: str = typer.Option(
        "1",
        "--channel",
        help="Imaging channel (1, 2, …) or import label such as imaris_488_cells (tidy tables only).",
    ),
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
    from lightsuite.analysis.viz.io import load_region_plot_table, parse_plot_channel
    from lightsuite.analysis.viz.plots import plot_lr_scatter

    csv_path = _resolve_plot_input(input_path=input_path, brain_config=brain_config)
    table = load_region_plot_table(
        csv_path, channel=parse_plot_channel(channel), metric=metric
    )
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


@analysis_app.command("plot-top-regions")
def analysis_plot_top_regions(
    input_path: str | None = typer.Option(None, "--input", "-i", help="region_stats or intensities CSV."),
    brain_config: str | None = typer.Option(
        None, "--config", "-c", help="Brain YAML; uses volume_registered/region_stats.csv."
    ),
    output: str = typer.Option(..., "--output", "-o", help="Output PNG path."),
    top_n: int = typer.Option(10, "--top-n", "-n", help="Number of top regions to show."),
    channel: str = typer.Option(
        "1",
        "--channel",
        help="Imaging channel (1, 2, …) or import label such as imaris_488_cells (tidy tables only).",
    ),
    metric: str = typer.Option("cell_count", "--metric", help="Metric to plot."),
    title: str | None = typer.Option(None, "--title", help="Figure title."),
    keep_division: list[str] | None = typer.Option(
        None, "--keep-division", help="Only these divisions (repeatable)."
    ),
    exclude_division: list[str] | None = typer.Option(
        None, "--exclude-division", help="Division names to omit (repeatable)."
    ),
    dpi: int = typer.Option(200, "--dpi", help="Figure DPI."),
) -> None:
    """Bar plot: left vs right for the top N regions by total metric value."""
    from lightsuite.analysis.viz.io import load_region_plot_table, parse_plot_channel
    from lightsuite.analysis.viz.plots import plot_top_region_bars

    csv_path = _resolve_plot_input(input_path=input_path, brain_config=brain_config)
    table = load_region_plot_table(
        csv_path, channel=parse_plot_channel(channel), metric=metric
    )
    out = Path(output).expanduser()
    plot_top_region_bars(
        table,
        top_n=top_n,
        metric=metric,
        title=title,
        output_path=out,
        dpi=dpi,
        keep_divisions=keep_division,
        exclude_divisions=exclude_division,
    )
    typer.echo(f"Saved: {out.resolve()} (top_n={top_n})")


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


def _resolve_cord_plot_context(
    *,
    input_path: str | None,
    spinal_config: str | None,
) -> tuple["Path", Path | None, list[str]]:
    from pathlib import Path

    from lightsuite.analysis.viz.cord_io import load_segment_order, resolve_cord_region_stats_from_config
    from lightsuite.atlas.fiederling import resolve_fiederling_paths

    if input_path and spinal_config:
        raise typer.BadParameter("Use only one of --input or --config.")
    if spinal_config:
        from lightsuite.config.loader import load_spinal_config

        cfg = load_spinal_config(spinal_config)
        stats_path = resolve_cord_region_stats_from_config(spinal_config)
        segments_csv = resolve_fiederling_paths(cfg.atlas.atlas_dir).segments_csv
        return stats_path, segments_csv, load_segment_order(segments_csv)
    if input_path:
        return Path(input_path).expanduser().resolve(), None, []
    raise typer.BadParameter("Provide --input or --config (spinal cord YAML).")


def _resolve_cord_channel_list(
    channels: str | None,
    spinal_config: str | None,
) -> list[str]:
    from lightsuite.analysis.viz.cord_io import parse_plot_channels

    channel_list = parse_plot_channels(channels)
    if channel_list:
        return channel_list
    if spinal_config is None:
        raise typer.BadParameter("Provide --channels or --config with analysis.point_labels.")
    from lightsuite.config.loader import load_spinal_config

    cfg = load_spinal_config(spinal_config)
    labels = cfg.analysis.point_labels
    if not labels:
        raise typer.BadParameter("No --channels given and analysis.point_labels is empty in config.")
    return [str(label) for label in labels]


@analysis_app.command("plot-cord-structure")
def analysis_plot_cord_structure(
    output: str = typer.Option(..., "--output", "-o", help="Output PNG path."),
    input_path: str | None = typer.Option(None, "--input", "-i", help="Cord region_stats.csv."),
    spinal_config: str | None = typer.Option(
        None, "--config", "-c", help="Spinal YAML; uses volume_registered/region_stats.csv."
    ),
    channel: str = typer.Option("1", "--channel", help="Imaging channel or import label."),
    metric: str = typer.Option("median_intensity", "--metric", help="Metric to plot."),
    hemisphere: str | None = typer.Option(
        None,
        "--hemisphere",
        help="Filter to left or right hemisegment (requires split_hemispheres in region-stats).",
    ),
    title: str | None = typer.Option(None, "--title", help="Figure title."),
    dpi: int = typer.Option(200, "--dpi", help="Figure DPI."),
    crop_empty_segments: bool = typer.Option(
        True,
        "--crop-empty-segments/--no-crop-empty-segments",
        help="Drop leading/trailing segments with no data (default: on).",
    ),
) -> None:
    """Heatmap: structure (rows) × rostrocaudal segment (columns)."""
    from lightsuite.analysis.viz.cord_io import filter_cord_stats, load_cord_stats_csv, parse_plot_channel
    from lightsuite.analysis.viz.cord_plots import plot_cord_structure_heatmap

    stats_path, _segments_csv, segment_order = _resolve_cord_plot_context(
        input_path=input_path, spinal_config=spinal_config
    )
    table = filter_cord_stats(
        load_cord_stats_csv(stats_path),
        channel=parse_plot_channel(channel),
        metric=metric,
        rollup_level="structure",
        hemisphere=hemisphere,
    )
    out = Path(output).expanduser()
    plot_cord_structure_heatmap(
        table,
        segment_order=segment_order or None,
        title=title,
        metric=metric,
        output_path=out,
        dpi=dpi,
        crop_empty_segments=crop_empty_segments,
    )
    typer.echo(f"Saved: {out.resolve()}")


@analysis_app.command("plot-cord-structure-hemisphere-panel")
def analysis_plot_cord_structure_hemisphere_panel(
    output: str = typer.Option(..., "--output", "-o", help="Output PNG path."),
    input_path: str | None = typer.Option(None, "--input", "-i", help="Cord region_stats.csv."),
    spinal_config: str | None = typer.Option(
        None, "--config", "-c", help="Spinal YAML; uses volume_registered/region_stats.csv."
    ),
    channel: str = typer.Option("1", "--channel", help="Imaging channel or import label."),
    metric: str = typer.Option("median_intensity", "--metric", help="Metric to plot."),
    rollup_level: str = typer.Option("structure", "--rollup-level", help="Rollup level."),
    title: str | None = typer.Option(None, "--title", help="Figure title."),
    dpi: int = typer.Option(200, "--dpi", help="Figure DPI."),
    crop_empty_segments: bool = typer.Option(
        True,
        "--crop-empty-segments/--no-crop-empty-segments",
        help="Drop leading/trailing segments with no data (default: on).",
    ),
) -> None:
    """Side-by-side structure heatmaps for left and right hemisegments."""
    from lightsuite.analysis.viz.cord_io import load_cord_stats_csv, parse_plot_channel
    from lightsuite.analysis.viz.cord_plots import plot_cord_structure_hemisphere_panel

    stats_path, _segments_csv, segment_order = _resolve_cord_plot_context(
        input_path=input_path, spinal_config=spinal_config
    )
    out = Path(output).expanduser()
    plot_cord_structure_hemisphere_panel(
        load_cord_stats_csv(stats_path),
        channel=parse_plot_channel(channel),
        metric=metric,
        rollup_level=rollup_level,
        segment_order=segment_order or None,
        title=title,
        output_path=out,
        dpi=dpi,
        crop_empty_segments=crop_empty_segments,
    )
    typer.echo(f"Saved: {out.resolve()}")


@analysis_app.command("plot-cord-division-profile")
def analysis_plot_cord_division_profile(
    output: str = typer.Option(..., "--output", "-o", help="Output PNG path."),
    input_path: str | None = typer.Option(None, "--input", "-i", help="Cord region_stats.csv."),
    spinal_config: str | None = typer.Option(
        None, "--config", "-c", help="Spinal YAML; uses volume_registered/region_stats.csv."
    ),
    channel: str = typer.Option("1", "--channel", help="Imaging channel or import label."),
    metric: str = typer.Option("median_intensity", "--metric", help="Metric to plot."),
    title: str | None = typer.Option(None, "--title", help="Figure title."),
    z_voxel_um: float = typer.Option(20.0, "--z-voxel-um", help="Atlas Z voxel size in µm for mm axis."),
    dpi: int = typer.Option(200, "--dpi", help="Figure DPI."),
    crop_empty_segments: bool = typer.Option(
        True,
        "--crop-empty-segments/--no-crop-empty-segments",
        help="Drop leading/trailing segments with no data (default: on).",
    ),
) -> None:
    """Line plot: GM/WM signal vs rostrocaudal position."""
    import pandas as pd

    from lightsuite.analysis.viz.cord_io import filter_cord_stats, load_cord_stats_csv, parse_plot_channel
    from lightsuite.analysis.viz.cord_plots import plot_cord_division_profile

    stats_path, segments_csv, _segment_order = _resolve_cord_plot_context(
        input_path=input_path, spinal_config=spinal_config
    )
    if segments_csv is None or not segments_csv.is_file():
        raise typer.BadParameter("Provide --config (spinal YAML) so Segments.csv can be resolved.")
    table = filter_cord_stats(
        load_cord_stats_csv(stats_path),
        channel=parse_plot_channel(channel),
        metric=metric,
        rollup_level="division",
    )
    out = Path(output).expanduser()
    plot_cord_division_profile(
        table,
        pd.read_csv(segments_csv),
        title=title,
        metric=metric,
        z_voxel_um=z_voxel_um,
        output_path=out,
        dpi=dpi,
        crop_empty_segments=crop_empty_segments,
    )
    typer.echo(f"Saved: {out.resolve()}")


@analysis_app.command("plot-cord-segment-bars")
def analysis_plot_cord_segment_bars(
    output: str = typer.Option(..., "--output", "-o", help="Output PNG path."),
    input_path: str | None = typer.Option(None, "--input", "-i", help="Cord region_stats.csv."),
    spinal_config: str | None = typer.Option(
        None, "--config", "-c", help="Spinal YAML; uses volume_registered/region_stats.csv."
    ),
    channel: str = typer.Option("1", "--channel", help="Imaging channel or import label."),
    metric: str = typer.Option("cell_count", "--metric", help="Metric to plot."),
    title: str | None = typer.Option(None, "--title", help="Figure title."),
    min_total: float = typer.Option(
        0.0,
        "--min-total",
        help="Drop segments whose total is below this threshold.",
    ),
    dpi: int = typer.Option(200, "--dpi", help="Figure DPI."),
) -> None:
    """Bar chart: total metric per segment (summed over finest regions)."""
    from lightsuite.analysis.viz.cord_io import filter_cord_stats, load_cord_stats_csv, parse_plot_channel
    from lightsuite.analysis.viz.cord_plots import plot_cord_segment_bars

    stats_path, _segments_csv, segment_order = _resolve_cord_plot_context(
        input_path=input_path, spinal_config=spinal_config
    )
    table = filter_cord_stats(
        load_cord_stats_csv(stats_path),
        channel=parse_plot_channel(channel),
        metric=metric,
        rollup_level="region",
    )
    out = Path(output).expanduser()
    plot_cord_segment_bars(
        table,
        segment_order=segment_order or None,
        title=title,
        metric=metric,
        min_total=min_total,
        output_path=out,
        dpi=dpi,
    )
    typer.echo(f"Saved: {out.resolve()}")


@analysis_app.command("plot-cord-segment-grouped-bars")
def analysis_plot_cord_segment_grouped_bars(
    output: str = typer.Option(..., "--output", "-o", help="Output PNG path."),
    input_path: str | None = typer.Option(None, "--input", "-i", help="Cord region_stats.csv."),
    spinal_config: str | None = typer.Option(
        None, "--config", "-c", help="Spinal YAML; uses volume_registered/region_stats.csv."
    ),
    channels: str | None = typer.Option(
        None,
        "--channels",
        help="Comma-separated import labels (default: analysis.point_labels from --config).",
    ),
    metric: str = typer.Option("cell_count", "--metric", help="Metric to plot."),
    title: str | None = typer.Option(None, "--title", help="Figure title."),
    min_total: float = typer.Option(
        0.0,
        "--min-total",
        help="Drop segments whose summed metric across labels is below this threshold.",
    ),
    show_composition: bool = typer.Option(
        True,
        "--composition/--no-composition",
        help="Add a 100% stacked composition panel under the counts (default: on).",
    ),
    dpi: int = typer.Option(200, "--dpi", help="Figure DPI."),
) -> None:
    """Grouped bar chart: compare coloc import labels per rostrocaudal segment."""
    from lightsuite.analysis.viz.cord_io import (
        filter_cord_stats_multi,
        load_cord_stats_csv,
        parse_plot_channels,
    )
    from lightsuite.analysis.viz.cord_plots import plot_cord_segment_grouped_bars

    stats_path, _segments_csv, segment_order = _resolve_cord_plot_context(
        input_path=input_path, spinal_config=spinal_config
    )
    channel_list = _resolve_cord_channel_list(channels, spinal_config)

    table = filter_cord_stats_multi(
        load_cord_stats_csv(stats_path),
        channels=channel_list,
        metric=metric,
        rollup_level="region",
    )
    out = Path(output).expanduser()
    plot_cord_segment_grouped_bars(
        table,
        channels=channel_list,
        segment_order=segment_order or None,
        title=title,
        metric=metric,
        min_total=min_total,
        show_composition=show_composition,
        output_path=out,
        dpi=dpi,
    )
    typer.echo(f"Saved: {out.resolve()}")


@analysis_app.command("plot-cord-structure-panel")
def analysis_plot_cord_structure_panel(
    output: str = typer.Option(..., "--output", "-o", help="Output PNG path."),
    input_path: str | None = typer.Option(None, "--input", "-i", help="Cord region_stats.csv."),
    spinal_config: str | None = typer.Option(
        None, "--config", "-c", help="Spinal YAML; uses volume_registered/region_stats.csv."
    ),
    channels: str | None = typer.Option(
        None,
        "--channels",
        help="Comma-separated import labels (default: analysis.point_labels from --config).",
    ),
    metric: str = typer.Option("cell_count", "--metric", help="Metric to plot."),
    rollup_level: str = typer.Option("structure", "--rollup-level", help="Rollup level."),
    title: str | None = typer.Option(None, "--title", help="Figure title."),
    dpi: int = typer.Option(200, "--dpi", help="Figure DPI."),
    crop_empty_segments: bool = typer.Option(
        True,
        "--crop-empty-segments/--no-crop-empty-segments",
        help="Drop leading/trailing segments with no data (default: on).",
    ),
) -> None:
    """Side-by-side structure heatmaps for several import labels."""
    from lightsuite.analysis.viz.cord_io import load_cord_stats_csv
    from lightsuite.analysis.viz.cord_plots import plot_cord_structure_panel

    stats_path, _segments_csv, segment_order = _resolve_cord_plot_context(
        input_path=input_path, spinal_config=spinal_config
    )
    channel_list = _resolve_cord_channel_list(channels, spinal_config)
    out = Path(output).expanduser()
    plot_cord_structure_panel(
        load_cord_stats_csv(stats_path),
        channels=channel_list,
        metric=metric,
        rollup_level=rollup_level,
        segment_order=segment_order or None,
        title=title,
        output_path=out,
        dpi=dpi,
        crop_empty_segments=crop_empty_segments,
    )
    typer.echo(f"Saved: {out.resolve()}")


@analysis_app.command("plot-cord-top-regions")
def analysis_plot_cord_top_regions(
    output: str = typer.Option(..., "--output", "-o", help="Output PNG path."),
    input_path: str | None = typer.Option(None, "--input", "-i", help="Cord region_stats.csv."),
    spinal_config: str | None = typer.Option(
        None, "--config", "-c", help="Spinal YAML; uses volume_registered/region_stats.csv."
    ),
    channel: str = typer.Option(..., "--channel", help="Imaging channel or import label."),
    metric: str = typer.Option("cell_count", "--metric", help="Metric to plot."),
    top_n: int = typer.Option(15, "--top-n", help="Number of top region × segment rows."),
    segment: str | None = typer.Option(None, "--segment", help="Restrict to one segment (e.g. L5)."),
    title: str | None = typer.Option(None, "--title", help="Figure title."),
    dpi: int = typer.Option(200, "--dpi", help="Figure DPI."),
) -> None:
    """Horizontal bar chart of top region × segment combinations."""
    from lightsuite.analysis.viz.cord_io import filter_cord_stats, load_cord_stats_csv, parse_plot_channel
    from lightsuite.analysis.viz.cord_plots import plot_cord_top_regions

    stats_path, _segments_csv, _segment_order = _resolve_cord_plot_context(
        input_path=input_path, spinal_config=spinal_config
    )
    table = filter_cord_stats(
        load_cord_stats_csv(stats_path),
        channel=parse_plot_channel(channel),
        metric=metric,
        rollup_level="region",
    )
    out = Path(output).expanduser()
    plot_cord_top_regions(
        table,
        top_n=top_n,
        segment=segment,
        title=title,
        metric=metric,
        output_path=out,
        dpi=dpi,
    )
    typer.echo(f"Saved: {out.resolve()}")


def _parse_segment_list(value: str | None) -> list[str] | None:
    if value is None:
        return None
    segments = [part.strip() for part in str(value).split(",") if part.strip()]
    return segments or None


@analysis_app.command("plot-cord-laminae-pct-gm")
def analysis_plot_cord_laminae_pct_gm(
    output: str = typer.Option(..., "--output", "-o", help="Output PNG path."),
    input_path: str | None = typer.Option(None, "--input", "-i", help="Cord region_stats.csv."),
    spinal_config: str | None = typer.Option(
        None, "--config", "-c", help="Spinal YAML; uses volume_registered/region_stats.csv."
    ),
    intensity_channel: str = typer.Option("1", "--intensity-channel", help="Imaging channel for signal intensity."),
    cell_channel: str = typer.Option(..., "--cell-channel", help="Import label for cell counts."),
    intensity_metric: str = typer.Option(
        "median_intensity",
        "--intensity-metric",
        help="Intensity metric for % GM bars.",
    ),
    segments: str | None = typer.Option(
        None,
        "--segments",
        help="Comma-separated segment filter (e.g. C4,C5,C6,C7).",
    ),
    title: str | None = typer.Option(None, "--title", help="Figure title."),
    dpi: int = typer.Option(200, "--dpi", help="Figure DPI."),
) -> None:
    """Bar chart: % GM occupied by signal intensity vs. cell density across Rexed laminae."""
    from lightsuite.analysis.viz.cord_io import load_cord_stats_csv, parse_plot_channel
    from lightsuite.analysis.viz.cord_plots import plot_cord_laminae_pct_gm_bars

    stats_path, _segments_csv, _segment_order = _resolve_cord_plot_context(
        input_path=input_path, spinal_config=spinal_config
    )
    out = Path(output).expanduser()
    plot_cord_laminae_pct_gm_bars(
        load_cord_stats_csv(stats_path),
        intensity_channel=parse_plot_channel(intensity_channel),
        cell_channel=parse_plot_channel(cell_channel),
        segments=_parse_segment_list(segments),
        intensity_metric=intensity_metric,
        title=title,
        output_path=out,
        dpi=dpi,
    )
    typer.echo(f"Saved: {out.resolve()}")


@analysis_app.command("plot-cord-laminae-level-bars")
def analysis_plot_cord_laminae_level_bars(
    output: str = typer.Option(..., "--output", "-o", help="Output PNG path."),
    input_path: str | None = typer.Option(None, "--input", "-i", help="Cord region_stats.csv."),
    spinal_config: str | None = typer.Option(
        None, "--config", "-c", help="Spinal YAML; uses volume_registered/region_stats.csv."
    ),
    channel: str = typer.Option("1", "--channel", help="Imaging channel or import label."),
    metric: str = typer.Option("median_intensity", "--metric", help="Metric to plot."),
    levels: str = typer.Option("C,T,L", "--levels", help="Comma-separated cord levels (C,T,L,S)."),
    segments: str | None = typer.Option(
        None,
        "--segments",
        help="Comma-separated segment filter.",
    ),
    title: str | None = typer.Option(None, "--title", help="Figure title."),
    dpi: int = typer.Option(200, "--dpi", help="Figure DPI."),
) -> None:
    """Grouped bar chart: metric per Rexed lamina averaged across cervical/thoracic/lumbar levels."""
    from lightsuite.analysis.viz.cord_io import load_cord_stats_csv, parse_plot_channel
    from lightsuite.analysis.viz.cord_plots import plot_cord_laminae_level_bars

    stats_path, _segments_csv, _segment_order = _resolve_cord_plot_context(
        input_path=input_path, spinal_config=spinal_config
    )
    level_tuple = tuple(part.strip() for part in levels.split(",") if part.strip())
    out = Path(output).expanduser()
    plot_cord_laminae_level_bars(
        load_cord_stats_csv(stats_path),
        channel=parse_plot_channel(channel),
        metric=metric,
        segments=_parse_segment_list(segments),
        levels=level_tuple,  # type: ignore[arg-type]
        title=title,
        output_path=out,
        dpi=dpi,
    )
    typer.echo(f"Saved: {out.resolve()}")


@analysis_app.command("plot-cord-df-subregion-heatmap")
def analysis_plot_cord_df_subregion_heatmap(
    output: str = typer.Option(..., "--output", "-o", help="Output PNG path."),
    input_path: str | None = typer.Option(None, "--input", "-i", help="Cord region_stats.csv."),
    spinal_config: str | None = typer.Option(
        None, "--config", "-c", help="Spinal YAML; uses volume_registered/region_stats.csv."
    ),
    channel: str = typer.Option("1", "--channel", help="Imaging channel or import label."),
    metric: str = typer.Option("median_intensity", "--metric", help="Metric to plot."),
    segments: str | None = typer.Option(
        None,
        "--segments",
        help="Comma-separated segment filter.",
    ),
    include_parent_df: bool = typer.Option(
        True,
        "--include-parent-df/--no-include-parent-df",
        help="Include combined dorsal funiculus (df) row from structure rollup.",
    ),
    hemisphere: str | None = typer.Option(
        None,
        "--hemisphere",
        help="Filter to left or right hemisegment.",
    ),
    title: str | None = typer.Option(None, "--title", help="Figure title."),
    dpi: int = typer.Option(200, "--dpi", help="Figure DPI."),
    crop_empty_segments: bool = typer.Option(
        True,
        "--crop-empty-segments/--no-crop-empty-segments",
        help="Drop leading/trailing segments with no data (default: on).",
    ),
) -> None:
    """Heatmap: dorsal funiculus subregions (dcs, cu, gr, psdc, df) × rostrocaudal segment."""
    from lightsuite.analysis.viz.cord_io import load_cord_stats_csv, parse_plot_channel
    from lightsuite.analysis.viz.cord_plots import plot_cord_df_subregion_heatmap

    stats_path, _segments_csv, segment_order = _resolve_cord_plot_context(
        input_path=input_path, spinal_config=spinal_config
    )
    out = Path(output).expanduser()
    plot_cord_df_subregion_heatmap(
        load_cord_stats_csv(stats_path),
        channel=parse_plot_channel(channel),
        metric=metric,
        segments=_parse_segment_list(segments),
        segment_order=segment_order or None,
        include_parent_df=include_parent_df,
        hemisphere=hemisphere,
        title=title,
        output_path=out,
        dpi=dpi,
        crop_empty_segments=crop_empty_segments,
    )
    typer.echo(f"Saved: {out.resolve()}")


@analysis_app.command("plot-cord-horn-heatmap")
def analysis_plot_cord_horn_heatmap(
    output: str = typer.Option(..., "--output", "-o", help="Output PNG path."),
    input_path: str | None = typer.Option(None, "--input", "-i", help="Cord region_stats.csv."),
    spinal_config: str | None = typer.Option(
        None, "--config", "-c", help="Spinal YAML; uses volume_registered/region_stats.csv."
    ),
    channel: str = typer.Option("1", "--channel", help="Imaging channel or import label."),
    metric: str = typer.Option("median_intensity", "--metric", help="Metric to plot."),
    hemisphere: str | None = typer.Option(None, "--hemisphere", help="Filter to left or right."),
    title: str | None = typer.Option(None, "--title", help="Figure title."),
    dpi: int = typer.Option(200, "--dpi", help="Figure DPI."),
    crop_empty_segments: bool = typer.Option(True, "--crop-empty-segments/--no-crop-empty-segments"),
) -> None:
    """Heatmap: dorsal/ventral/central horn (rows) × rostrocaudal segment (columns)."""
    from lightsuite.analysis.viz.cord_io import filter_cord_stats, load_cord_stats_csv, parse_plot_channel
    from lightsuite.analysis.viz.cord_plots import plot_cord_horn_heatmap

    stats_path, _segments_csv, segment_order = _resolve_cord_plot_context(
        input_path=input_path, spinal_config=spinal_config
    )
    table = filter_cord_stats(
        load_cord_stats_csv(stats_path),
        channel=parse_plot_channel(channel),
        metric=metric,
        rollup_level="horn",
        hemisphere=hemisphere,
    )
    out = Path(output).expanduser()
    plot_cord_horn_heatmap(
        table,
        segment_order=segment_order or None,
        title=title,
        metric=metric,
        output_path=out,
        dpi=dpi,
        crop_empty_segments=crop_empty_segments,
    )
    typer.echo(f"Saved: {out.resolve()}")


@analysis_app.command("plot-cord-segment-anatomy-slice")
def analysis_plot_cord_segment_anatomy_slice(
    output: str = typer.Option(..., "--output", "-o", help="Output PNG path."),
    segments: str = typer.Option(
        ...,
        "--segments",
        help="Comma-separated segment list (e.g. T10,T12,T13,L2,S1).",
    ),
    spinal_config: str = typer.Option(
        ..., "--config", "-c", help="Spinal YAML (registered volumes + region_stats)."
    ),
    channel: str = typer.Option("1", "--channel", help="Imaging channel or import label."),
    metric: str = typer.Option("median_intensity", "--metric", help="Metric to plot."),
    rollup_level: str = typer.Option(
        "region",
        "--rollup-level",
        help="Rollup level for region values (region, structure, division, horn).",
    ),
    hemisphere: str | None = typer.Option(
        None,
        "--hemisphere",
        help="Optional left/right hemisegment mask (single-side mode only).",
    ),
    hemisphere_panel: bool = typer.Option(
        False,
        "--hemisphere-panel",
        help="Left/right columns per segment (requires hemisphere_registered.tiff).",
    ),
    hemisphere_composite: bool = typer.Option(
        False,
        "--hemisphere-composite",
        help="Both hemisegments on one cross-section per segment.",
    ),
    title: str | None = typer.Option(None, "--title", help="Figure title."),
    dpi: int = typer.Option(200, "--dpi", help="Figure DPI."),
    cmap: str = typer.Option("inferno", "--cmap", help="Matplotlib colormap."),
    vmin: float | None = typer.Option(None, "--vmin", help="Shared color scale minimum."),
    vmax: float | None = typer.Option(None, "--vmax", help="Shared color scale maximum."),
    ncol_max: int = typer.Option(4, "--ncol-max", help="Maximum columns in the panel grid."),
) -> None:
    """Per-segment anatomical heatmaps from registered annotation slices."""
    from lightsuite.analysis.viz.cord_io import load_cord_stats_csv, parse_plot_channel
    from lightsuite.analysis.viz.cord_segment_anatomy import (
        plot_cord_segment_anatomy_slice,
        plot_cord_segment_anatomy_slice_hemisphere_composite,
        plot_cord_segment_anatomy_slice_hemisphere_panel,
    )
    from lightsuite.config.loader import load_spinal_config
    from lightsuite.export.cord_registered import discover_registered_cord_paths

    if hemisphere_panel and hemisphere_composite:
        raise typer.BadParameter("Use only one of --hemisphere-panel or --hemisphere-composite.")
    if hemisphere and (hemisphere_panel or hemisphere_composite):
        raise typer.BadParameter(
            "--hemisphere is for single-side mode; use --hemisphere-panel or --hemisphere-composite instead."
        )

    segment_list = _parse_segment_list(segments)
    if not segment_list:
        raise typer.BadParameter("--segments must list at least one segment.")

    stats_path, segments_csv, _segment_order = _resolve_cord_plot_context(
        input_path=None, spinal_config=spinal_config
    )
    if segments_csv is None:
        raise typer.BadParameter("Segments.csv not found from config atlas_dir.")

    cfg = load_spinal_config(spinal_config)
    reg_paths = discover_registered_cord_paths(cfg)
    hem_path = reg_paths.volume_registered_dir / "hemisphere_registered.tiff"
    if (hemisphere_panel or hemisphere_composite) and not hem_path.is_file():
        raise typer.BadParameter(
            f"Missing {hem_path}. Run region-stats with analysis.split_hemispheres: true."
        )

    out = Path(output).expanduser()
    stats = load_cord_stats_csv(stats_path)
    channel_parsed = parse_plot_channel(channel)
    hemisphere_flip = bool(getattr(cfg.analysis, "hemisphere_flip", False))

    if hemisphere_panel:
        plot_cord_segment_anatomy_slice_hemisphere_panel(
            stats,
            segments=segment_list,
            annotation_path=reg_paths.annotation_path,
            segments_csv=segments_csv,
            hemisphere_volume_path=hem_path,
            channel=channel_parsed,
            metric=metric,
            rollup_level=rollup_level,
            hemisphere_flip=hemisphere_flip,
            title=title,
            output_path=out,
            dpi=dpi,
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
        )
    elif hemisphere_composite:
        plot_cord_segment_anatomy_slice_hemisphere_composite(
            stats,
            segments=segment_list,
            annotation_path=reg_paths.annotation_path,
            segments_csv=segments_csv,
            hemisphere_volume_path=hem_path,
            channel=channel_parsed,
            metric=metric,
            rollup_level=rollup_level,
            hemisphere_flip=hemisphere_flip,
            title=title,
            output_path=out,
            dpi=dpi,
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
            ncol_max=ncol_max,
        )
    else:
        hemisphere_path = hem_path if hem_path.is_file() else None
        plot_cord_segment_anatomy_slice(
            stats,
            segments=segment_list,
            annotation_path=reg_paths.annotation_path,
            segments_csv=segments_csv,
            hemisphere_volume_path=hemisphere_path,
            channel=channel_parsed,
            metric=metric,
            rollup_level=rollup_level,
            hemisphere=hemisphere,
            title=title,
            output_path=out,
            dpi=dpi,
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
            ncol_max=ncol_max,
        )
    typer.echo(f"Saved: {out.resolve()}")


@analysis_app.command("plot-cord-segment-anatomy-bgh")
def analysis_plot_cord_segment_anatomy_bgh(
    output: str = typer.Option(..., "--output", "-o", help="Output PNG path."),
    segments: str = typer.Option(
        ...,
        "--segments",
        help="Comma-separated segment list (e.g. T10,T12,T13,L2,S1).",
    ),
    spinal_config: str = typer.Option(
        ..., "--config", "-c", help="Spinal YAML (uses volume_registered/region_stats.csv)."
    ),
    channel: str = typer.Option("1", "--channel", help="Imaging channel or import label."),
    metric: str = typer.Option("median_intensity", "--metric", help="Metric to plot."),
    rollup_level: str = typer.Option(
        "region",
        "--rollup-level",
        help="Rollup level for region values (region, structure, division, horn).",
    ),
    hemisphere: str | None = typer.Option(
        None,
        "--hemisphere",
        help="Optional left or right hemisegment (brainrender mesh clipping).",
    ),
    title: str | None = typer.Option(None, "--title", help="Figure title."),
    dpi: int = typer.Option(200, "--dpi", help="Figure DPI."),
    cmap: str = typer.Option("inferno", "--cmap", help="Matplotlib colormap."),
    vmin: float | None = typer.Option(None, "--vmin", help="Shared color scale minimum."),
    vmax: float | None = typer.Option(None, "--vmax", help="Shared color scale maximum."),
    ncol_max: int = typer.Option(3, "--ncol-max", help="Maximum columns in the panel grid."),
    thickness_um: float = typer.Option(
        400.0,
        "--thickness-um",
        help="Slice thickness passed to brainglobe-heatmap.",
    ),
) -> None:
    """Per-segment anatomical heatmaps via brainglobe-heatmap (optional viz extra)."""
    from lightsuite.analysis.viz.cord_io import load_cord_stats_csv, parse_plot_channel
    from lightsuite.analysis.viz.cord_segment_anatomy import plot_cord_segment_anatomy_bgh

    segment_list = _parse_segment_list(segments)
    if not segment_list:
        raise typer.BadParameter("--segments must list at least one segment.")

    stats_path, segments_csv, _segment_order = _resolve_cord_plot_context(
        input_path=None, spinal_config=spinal_config
    )
    if segments_csv is None:
        raise typer.BadParameter("Segments.csv not found from config atlas_dir.")

    out = Path(output).expanduser()
    plot_cord_segment_anatomy_bgh(
        load_cord_stats_csv(stats_path),
        segments=segment_list,
        segments_csv=segments_csv,
        channel=parse_plot_channel(channel),
        metric=metric,
        rollup_level=rollup_level,
        hemisphere=hemisphere,
        title=title,
        output_path=out,
        dpi=dpi,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        ncol_max=ncol_max,
        thickness_um=thickness_um,
    )
    typer.echo(f"Saved: {out.resolve()}")


@analysis_app.command("cord-coloc-overlap")
def analysis_cord_coloc_overlap(
    output: str = typer.Option(..., "--output", "-o", help="Output PNG path."),
    spinal_config: str = typer.Option(..., "--config", "-c", help="Spinal cord YAML config."),
    channels: str | None = typer.Option(
        None,
        "--channels",
        help="Comma-separated import labels (default: analysis.point_labels).",
    ),
    tolerance_voxels: float = typer.Option(
        2.0,
        "--tolerance-voxels",
        help="Max atlas-voxel distance to call two spots colocalized.",
    ),
    title: str | None = typer.Option(None, "--title", help="Figure title."),
    dpi: int = typer.Option(200, "--dpi", help="Figure DPI."),
) -> None:
    """Compute and plot pairwise colocalization overlap between imported spot labels."""
    from lightsuite.analysis.cord_coloc import run_cord_coloc_overlap
    from lightsuite.analysis.viz.cord_plots import plot_cord_coloc_overlap
    from lightsuite.config.loader import load_spinal_config

    cfg = load_spinal_config(spinal_config)
    channel_list = _resolve_cord_channel_list(channels, spinal_config)
    out = Path(output).expanduser()
    csv_path = out.with_suffix(".csv")
    result = run_cord_coloc_overlap(
        cfg,
        labels=channel_list,
        tolerance_voxels=tolerance_voxels,
        output_csv=csv_path,
    )
    plot_cord_coloc_overlap(
        result.summary,
        title=title,
        output_path=out,
        dpi=dpi,
        save_csv=False,
    )
    typer.echo(f"Saved: {out.resolve()}")
    if result.summary_path is not None:
        typer.echo(f"Summary: {result.summary_path.resolve()}")


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
    space: str = typer.Option(
        "atlas",
        "--space",
        help="QC space: atlas (registered_atlas.tif) or sample (registration grid).",
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
        space=space,
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


@multires_app.command("validate-config")
def multires_validate_config(
    config: str = typer.Option(..., "--config", "-c", help="Multires pipeline YAML config."),
) -> None:
    """Load and validate a manifest-driven multiresolution YAML config."""
    from lightsuite.config.loader import load_multires_config
    from lightsuite.multires.resolve import resolve_pair_manifest

    cfg = load_multires_config(config)
    manifest, manifest_path = resolve_pair_manifest(cfg)
    channel_note = ""
    if manifest.channels:
        ref = cfg.multires.registration.reference_channel or manifest.resolved_reference_channel()
        extra = manifest.non_reference_channels(
            reference_channel=ref,
            apply_transform_to=cfg.multires.registration.apply_transform_to,
        )
        channel_note = f" channels={','.join(manifest.channel_names())} ref={ref}"
        if extra:
            channel_note += f" apply_to={','.join(extra)}"
    typer.echo(
        f"Config valid: {cfg.sample.name} / {manifest.pair_label} "
        f"({Path(manifest.overview.volume_path).name} → {Path(manifest.roi.volume_path).name})"
        f"{channel_note}"
    )
    if cfg.multires.channels:
        typer.echo(f"Pair manifest: {manifest_path}")


@multires_app.command("inspect-geometry")
def multires_inspect_geometry(
    config: str = typer.Option(..., "--config", "-c", help="Multires pipeline YAML config."),
    headless: bool = typer.Option(
        False,
        "--headless",
        help="Score the configured lateral_flip without opening Napari.",
    ),
    write_config: bool = typer.Option(
        False,
        "--write-config",
        help="Patch multires.mesospim_geometry in the YAML (uses current toggles in GUI, or configured flip in headless mode).",
    ),
) -> None:
    """Interactive mesoSPIM geometry QC (lateral_flip toggles + physical NCC)."""
    from lightsuite.config.loader import load_multires_config
    from lightsuite.gui.inspect_geometry_multires import run_multires_inspect_geometry

    cfg = load_multires_config(config)
    result = run_multires_inspect_geometry(
        cfg,
        config_path=config,
        headless=headless,
        write_config=write_config,
    )
    if headless:
        typer.echo(
            f"lateral_flip={list(result.lateral_flip)}  "
            f"physical_ncc={result.physical_ncc:.3f}  "
            f"has_overlap={result.has_overlap}"
        )


@multires_app.command("match-points")
def multires_match_points(
    config: str = typer.Option(..., "--config", "-c", help="Multires pipeline YAML config."),
    headless: bool = typer.Option(
        False,
        "--headless",
        help="Create an empty landmark session without opening Napari (for tests).",
    ),
) -> None:
    """Interactive overview / ROI landmark placement (Napari)."""
    from lightsuite.config.loader import load_multires_config
    from lightsuite.gui.match_points_multires import run_multires_match_points

    cfg = load_multires_config(config)
    path = run_multires_match_points(cfg, headless=headless)
    typer.echo(f"Landmark session: {path}")


@multires_app.command("check-geometry")
def multires_check_geometry(
    config: str = typer.Option(..., "--config", "-c", help="Multires pipeline YAML config."),
    level: str = typer.Option(
        "full",
        "--level",
        help="Geometry QA depth: metadata-only, slice-qc, or full (loads overlap crop).",
    ),
) -> None:
    """Validate FOV overlap and write geometry QA artifacts."""
    from lightsuite.config.loader import load_multires_config
    from lightsuite.multires.config_models import MultiresGeometryCheckLevel
    from lightsuite.multires.runner import check_multires_geometry

    try:
        check_level = MultiresGeometryCheckLevel(level)
    except ValueError as exc:
        allowed = ", ".join(item.value for item in MultiresGeometryCheckLevel)
        raise typer.BadParameter(f"level must be one of: {allowed}") from exc

    check_multires_geometry(load_multires_config(config), level=check_level)


@multires_app.command("export-preview")
def multires_export_preview(
    config: str = typer.Option(..., "--config", "-c", help="Multires pipeline YAML config."),
    output_dir: str | None = typer.Option(
        None,
        "--output-dir",
        help="Directory for preview TIFF crops (default: save_path/geometry/alignment_preview/<pair>).",
    ),
    n_slices: int = typer.Option(5, "--n-slices", min=1, help="Number of Z slices to export per volume."),
    margin_um: float = typer.Option(
        100.0,
        "--margin-um",
        help="Extra XY margin around the overlap crop (µm).",
    ),
    projection: str = typer.Option(
        "slices",
        "--projection",
        help="Export mode: 'slices' (small Z stacks) or 'max' (XY max projection over overlap Z).",
    ),
) -> None:
    """Export small overview / ROI TIFF crops around the overlap for visual QC."""
    from lightsuite.config.loader import load_multires_config
    from lightsuite.multires.preview import export_alignment_preview_crops

    if projection not in ("slices", "max"):
        msg = f"Unsupported projection {projection!r}; use 'slices' or 'max'"
        raise typer.BadParameter(msg)
    out = Path(output_dir).expanduser() if output_dir is not None else None
    export_alignment_preview_crops(
        load_multires_config(config),
        output_dir=out,
        n_slices=n_slices,
        margin_um=margin_um,
        projection=projection,  # type: ignore[arg-type]
    )


@multires_app.command("register")
def multires_register(
    config: str = typer.Option(..., "--config", "-c", help="Multires pipeline YAML config."),
) -> None:
    """Register ROI stack to overview using a pair manifest and elastix."""
    from lightsuite.config.loader import load_multires_config
    from lightsuite.multires.runner import run_multires_registration

    run_multires_registration(load_multires_config(config))


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


@spinal_app.command("check-orientation")
def spinal_check_orientation(
    config: str = typer.Option(..., "--config", "-c", help="Spinal cord pipeline YAML config."),
    headless: bool = typer.Option(
        False,
        "--headless",
        help="Write cord_orientation.txt without opening Napari (for tests).",
    ),
) -> None:
    """Manually set the cord rostrocaudal direction (writes cord_orientation.txt)."""
    from lightsuite.config.loader import load_spinal_config
    from lightsuite.gui.orientation_cord import run_spinal_orientation

    cfg = load_spinal_config(config)
    path = run_spinal_orientation(cfg, headless=headless)
    typer.echo(f"Cord orientation saved to {path}")


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


@spinal_app.command("align-longitudinal")
def spinal_align_longitudinal(
    config: str = typer.Option(..., "--config", "-c", help="Spinal cord pipeline YAML config."),
    headless: bool = typer.Option(
        False,
        "--headless",
        help="Auto-confirm centered z-init anchors without Napari (for tests).",
    ),
) -> None:
    """Interactive rostrocaudal sample-to-atlas alignment before init-registration."""
    from lightsuite.config.loader import load_spinal_config
    from lightsuite.gui.align_longitudinal_cord import run_spinal_align_longitudinal

    cfg = load_spinal_config(config)
    path = run_spinal_align_longitudinal(cfg, headless=headless)
    typer.echo(f"Longitudinal correspondence: {path}")


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
    space: str | None = typer.Option(
        None,
        "--space",
        help="Output space: atlas, sample, or both (default: export.spaces in YAML).",
    ),
) -> None:
    """Export registered cord volumes (generateRegisteredCordVolume.m)."""
    from lightsuite.cli.spaces import parse_spaces_option
    from lightsuite.config.loader import load_spinal_config
    from lightsuite.export.cord_export import export_registered_cord_volumes

    cfg = load_spinal_config(config)
    result = export_registered_cord_volumes(cfg, spaces=parse_spaces_option(space))
    typer.echo(f"Registered volumes in {result.output_dir} ({len(result.channel_paths)} channels)")


@spinal_app.command("import-annotations")
def spinal_import_annotations(
    config: str = typer.Option(..., "--config", "-c", help="Spinal cord pipeline YAML config."),
    write_csv: bool = typer.Option(
        None,
        "--write-csv/--no-write-csv",
        help="Write atlas coordinate CSVs for point imports.",
    ),
) -> None:
    """Register native sample-space annotations into Fiederling atlas space."""
    from lightsuite.config.loader import load_spinal_config
    from lightsuite.import_.cord_import import run_cord_import_annotations

    cfg = load_spinal_config(config)
    results = run_cord_import_annotations(cfg, write_csv=write_csv)
    for item in results:
        typer.echo(f"{item.label}: {item.kind} — {item.n_atlas} atlas features")


@spinal_app.command("convert-imaris-spots")
def spinal_convert_imaris_spots(
    source: str = typer.Option(..., "--source", "-s", help="Imaris Spot_*_Detailed.csv export."),
    output_dir: str = typer.Option(..., "--output-dir", "-o", help="Directory for points.csv files."),
    voxel_um: str = typer.Option(
        ...,
        "--voxel-um",
        help=(
            "Size of one LightSuite native voxel in the same units as Imaris Position columns "
            "(Image Properties → Voxel Size), as comma-separated x,y,z. "
            "Use 1,1,1 when the .ims is 1 µm isotropic on the same grid; "
            "use hybrid values such as 1,1,1.8 when XY matches LightSuite indices but Z "
            "plane counts differ; use the microscope size (e.g. 1.8,1.8,1.8) only when "
            "Imaris is calibrated to that size."
        ),
    ),
    label_prefix: str = typer.Option("imaris", help="Filename prefix for converted CSVs."),
    shape_yxz: str | None = typer.Option(
        None,
        "--shape-yxz",
        help="Optional native grid Y,X,Z (e.g. from sample_reference.json) to warn on unit mismatch.",
    ),
) -> None:
    """Convert Imaris multi-component spot CSV to LightSuite points.csv files."""
    from pathlib import Path

    from lightsuite.import_.imaris import (
        convert_imaris_spots_by_component,
        list_imaris_component_names,
        warn_if_positions_look_like_voxel_indices,
    )

    parts = [float(v.strip()) for v in voxel_um.split(",")]
    if len(parts) != 3:
        msg = "--voxel-um must have exactly three values: x,y,z"
        raise typer.BadParameter(msg)

    source_path = Path(source).expanduser()
    if not source_path.is_file():
        msg = f"Source CSV not found: {source_path}"
        raise typer.BadParameter(msg)

    shape: tuple[int, int, int] | None = None
    if shape_yxz is not None:
        shape_parts = [int(v.strip()) for v in shape_yxz.split(",")]
        if len(shape_parts) != 3:
            msg = "--shape-yxz must have exactly three integers: Y,X,Z"
            raise typer.BadParameter(msg)
        shape = (shape_parts[0], shape_parts[1], shape_parts[2])

    warning = warn_if_positions_look_like_voxel_indices(
        source_path, voxel_um=parts, shape_yxz=shape
    )
    if warning:
        typer.echo(f"WARNING: {warning}", err=True)

    components = list_imaris_component_names(source_path)
    if not components:
        msg = f"No component labels found in {source_path}"
        raise typer.BadParameter(msg)

    written = convert_imaris_spots_by_component(
        source_path,
        Path(output_dir),
        voxel_um=parts,
        label_prefix=label_prefix,
    )
    for component, path in written.items():
        typer.echo(f"{component}: {path}")
    typer.echo(f"Wrote {len(written)} component CSV(s) to {output_dir}")


@spinal_app.command("region-stats")
def spinal_region_stats(
    config: str = typer.Option(..., "--config", "-c", help="Spinal cord pipeline YAML config."),
    count_points: bool = typer.Option(
        None,
        "--count-points/--no-count-points",
        help="Bin imported atlas-space points into cell counts (default: analysis.count_points).",
    ),
    parcellate_intensities: bool = typer.Option(
        None,
        "--parcellate-intensities/--no-parcellate-intensities",
        help="Median intensity per region × segment from exported channels (default: analysis.parcellate_intensities).",
    ),
    space: str | None = typer.Option(
        None,
        "--space",
        help="Stats space: atlas, sample, or both (default: analysis.stats_spaces).",
    ),
) -> None:
    """Assemble region_stats.csv with per-region, per-segment intensities and cell counts."""
    from lightsuite.analysis.cord_runner import run_cord_region_stats
    from lightsuite.cli.spaces import parse_spaces_option
    from lightsuite.config.loader import load_spinal_config

    cfg = load_spinal_config(config)
    result = run_cord_region_stats(
        cfg,
        count_points=count_points,
        parcellate_intensities=parcellate_intensities,
        stats_spaces=parse_spaces_option(space),
    )
    if result.combined_path is not None:
        typer.echo(
            f"Region stats: {result.combined_path} ({result.n_rows} rows, "
            f"{len(result.intensity_channels)} intensity channel(s), "
            f"{len(result.count_labels)} point source(s))"
        )
    else:
        typer.echo(
            "No region stats produced. Run 'lightsuite spinal export' for intensities "
            "and/or 'lightsuite spinal import-annotations' for point counts."
        )


@spinal_app.command("inspect-imports")
def spinal_inspect_imports(
    config: str = typer.Option(..., "--config", "-c", help="Spinal cord pipeline YAML config."),
    space: str = typer.Option(
        "atlas",
        "--space",
        help="Inspect space: atlas (Fiederling export grid) or sample (straightened 20 µm grid).",
    ),
    headless: bool = typer.Option(
        False,
        "--headless",
        help="Validate inspect inputs without opening Napari.",
    ),
    recompute_annotation: bool = typer.Option(
        False,
        "--recompute-annotation",
        help="Rebuild annotation_registered.tiff from transform_params.json (atlas space only).",
    ),
) -> None:
    """Napari QC: registered channels, atlas, and imported points."""
    from lightsuite.cli.spaces import parse_view_space_option
    from lightsuite.config.loader import load_spinal_config
    from lightsuite.gui.inspect_cord_imports import run_cord_inspect_imports

    cfg = load_spinal_config(config)
    inspect_space = parse_view_space_option(space)
    paths = run_cord_inspect_imports(
        cfg,
        space=inspect_space,  # type: ignore[arg-type]
        headless=headless,
        recompute_annotation=recompute_annotation,
    )
    if paths.sample_space_dir is not None:
        typer.echo(f"sample_space: {paths.sample_space_dir}")
    else:
        typer.echo(f"volume_registered: {paths.volume_registered_dir}")
    if paths.registered_channels:
        typer.echo(f"  channels: {sorted(paths.registered_channels)}")
    if paths.point_npz_paths:
        typer.echo(f"  point layers: {list(paths.point_npz_paths)}")


@spinal_app.command("view")
def spinal_view(
    config: str = typer.Option(..., "--config", "-c", help="Spinal cord pipeline YAML config."),
    space: str = typer.Option(
        "atlas",
        "--space",
        help="View space: atlas (Fiederling export grid) or sample (straightened 20 µm grid).",
    ),
    headless: bool = typer.Option(
        False,
        "--headless",
        help="Validate view inputs without opening Napari.",
    ),
    recompute_annotation: bool = typer.Option(
        False,
        "--recompute-annotation",
        help="Rebuild annotation_registered.tiff from transform_params.json (atlas space only).",
    ),
) -> None:
    """Open Napari with registered sample channel(s) and warped atlas annotation."""
    from lightsuite.cli.spaces import parse_view_space_option
    from lightsuite.config.loader import load_spinal_config
    from lightsuite.export.cord_sample_space import CordSampleSpaceInspectPaths
    from lightsuite.gui.view_registered_cord import run_spinal_registered_view

    cfg = load_spinal_config(config)
    view_space = parse_view_space_option(space)
    paths = run_spinal_registered_view(
        cfg,
        space=view_space,  # type: ignore[arg-type]
        headless=headless,
        recompute_annotation=recompute_annotation,
    )
    if isinstance(paths, CordSampleSpaceInspectPaths):
        typer.echo(f"Sample-space data: {paths.sample_space_dir}")
        typer.echo(f"Annotation: {paths.annotation_path.name}")
        if paths.channel_paths:
            typer.echo(f"Channels: {sorted(paths.channel_paths)}")
        if paths.point_npz_paths:
            typer.echo(f"Point layers: {list(paths.point_npz_paths)}")
    else:
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
