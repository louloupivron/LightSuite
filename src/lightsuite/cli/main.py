"""LightSuite Typer CLI entry point."""

from __future__ import annotations

from pathlib import Path

import typer

from lightsuite import __version__
from lightsuite.cli.doctor import doctor_command
from lightsuite.cli.config_cmd import config_app
from lightsuite.cli.manifest_build import build_manifest as multires_build_manifest_cmd
from lightsuite.cli.pipeline_commands import (
    register_brain_commands,
    register_multires_commands,
    register_spinal_commands,
)
from lightsuite.cli.workflows import workflow_app
from lightsuite.exceptions import LightsuiteConfigError

app = typer.Typer(
    name="lightsuite",
    help="LightSuite — mouse lightsheet and histology atlas registration.",
    no_args_is_help=True,
)
brain_app = typer.Typer(help="Brain lightsheet pipeline stages.")
spinal_app = typer.Typer(help="Spinal cord lightsheet pipeline stages.")
multires_app = typer.Typer(help="Manifest-driven overview ↔ ROI multiresolution registration.")
app.add_typer(brain_app, name="brain")
app.add_typer(spinal_app, name="spinal")
app.add_typer(multires_app, name="multires")
app.add_typer(config_app, name="config")
app.add_typer(workflow_app, name="workflow")

register_brain_commands(brain_app)
register_spinal_commands(spinal_app)
register_multires_commands(multires_app)


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

    try:
        cfg = load_config(config)
    except LightsuiteConfigError as exc:
        typer.secho(str(exc), fg=typer.colors.RED, err=True)
        raise typer.Exit(code=1) from exc
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
    space: str = typer.Option(
        "atlas",
        "--space",
        help="Inspect space: atlas (Perens/Allen export grid) or sample (20 µm registration grid).",
    ),
    headless: bool = typer.Option(
        False,
        "--headless",
        help="Validate inspect inputs without opening Napari.",
    ),
) -> None:
    """Napari QC: registered channels, atlas, and imported points/masks."""
    from lightsuite.cli.spaces import parse_view_space_option
    from lightsuite.config.loader import load_config
    from lightsuite.gui.inspect_brain_imports import run_brain_inspect_imports

    cfg = load_config(config)
    inspect_space = parse_view_space_option(space)
    paths = run_brain_inspect_imports(
        cfg,
        space=inspect_space,  # type: ignore[arg-type]
        headless=headless,
    )
    if paths.sample_space_dir is not None:
        typer.echo(f"sample_space: {paths.sample_space_dir}")
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


@multires_app.command("inspect-registration")
def multires_inspect_registration(
    config: str = typer.Option(..., "--config", "-c", help="Multires pipeline YAML config."),
    full_overview: bool = typer.Option(
        False,
        "--full-overview",
        help="Load the full-overview canvas instead of the overlap crop (much larger).",
    ),
    headless: bool = typer.Option(
        False,
        "--headless",
        help="Resolve and load layers without opening Napari (for tests).",
    ),
) -> None:
    """Compare the registered ROI against the overview in Napari."""
    from lightsuite.config.loader import load_multires_config
    from lightsuite.gui.inspect_registration_multires import run_multires_inspect_registration

    paths = run_multires_inspect_registration(
        load_multires_config(config),
        full_overview=full_overview,
        headless=headless,
    )
    typer.echo(
        f"overview: {paths.overview_path}\n"
        f"registered ROI channels: {', '.join(sorted(paths.registered_roi_paths))}"
    )


@multires_app.command("build-manifest")
def multires_build_manifest(
    vendor: str = typer.Option(
        ...,
        "--vendor",
        help="Acquisition vendor: mesospim or smartspim.",
    ),
    sample_name: str = typer.Option(..., "--sample-name", help="Sample identifier."),
    pair_label: str = typer.Option(..., "--pair-label", help="Unique pair label."),
    overview: Path | None = typer.Option(
        None,
        "--overview",
        help="Overview volume path (single-channel mesospim/smartspim).",
    ),
    roi: Path | None = typer.Option(None, "--roi", help="ROI volume path (single-channel)."),
    overview_meta: Path | None = typer.Option(
        None,
        "--overview-meta",
        help="mesoSPIM overview meta sidecar (required for stitched overview folders).",
    ),
    roi_meta: Path | None = typer.Option(None, "--roi-meta", help="mesoSPIM ROI meta sidecar."),
    channels_json: str | None = typer.Option(
        None,
        "--channels-json",
        help='Multichannel mesoSPIM paths as JSON, e.g. \'{"488":{"overview":"...","roi":"..."}}\'',
    ),
    reference_channel: str | None = typer.Option(
        None,
        "--reference-channel",
        help="Reference channel for multichannel manifests.",
    ),
    output: Path = typer.Option(
        ...,
        "--output",
        "-o",
        help="Output pair manifest JSON path.",
    ),
    lateral_flip_overview: str | None = typer.Option(
        None,
        "--lateral-flip-overview",
        help="mesoSPIM overview lateral_flip as '1,-1' (optional).",
    ),
    lateral_flip_roi: str | None = typer.Option(
        None,
        "--lateral-flip-roi",
        help="mesoSPIM ROI lateral_flip as '1,-1' (optional).",
    ),
) -> None:
    """Build a multires pair manifest JSON from mesoSPIM or SmartSPIM paths."""
    multires_build_manifest_cmd(
        vendor=vendor,
        sample_name=sample_name,
        pair_label=pair_label,
        overview=overview,
        roi=roi,
        overview_meta=overview_meta,
        roi_meta=roi_meta,
        channels_json=channels_json,
        reference_channel=reference_channel,
        output=output,
        lateral_flip_overview=lateral_flip_overview,
        lateral_flip_roi=lateral_flip_roi,
    )


@multires_app.command("import-annotations")
def multires_import_annotations(
    config: str = typer.Option(..., "--config", "-c", help="Multires pipeline YAML config."),
    no_csv: bool = typer.Option(False, "--no-csv", help="Skip writing the overview points CSV."),
    full_overview_canvas: bool | None = typer.Option(
        None,
        "--full-overview-canvas/--crop-only",
        help="Write warped masks on the full overview grid, or only the overlap crop.",
    ),
) -> None:
    """Warp ROI-native segmentation into overview-native space."""
    from lightsuite.config.loader import load_multires_config
    from lightsuite.multires.import_annotations import run_multires_import_annotations

    run_multires_import_annotations(
        load_multires_config(config),
        write_csv=False if no_csv else None,
        full_overview_canvas=full_overview_canvas,
    )


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
