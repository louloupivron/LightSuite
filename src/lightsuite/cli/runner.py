"""Run multiple pipeline stages in order with resume support."""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path
from typing import Any

import typer
from rich.console import Console

from lightsuite.cli.stages import (
    StageSpec,
    StageState,
    StageStatus,
    brain_stage_specs,
    brain_stage_statuses,
    multires_stage_specs,
    multires_stage_statuses,
    slice_stages,
    spinal_stage_specs,
    spinal_stage_statuses,
)
from lightsuite.config.loader import load_config, load_multires_config, load_spinal_config

console = Console()


def _stage_is_complete(status: StageStatus) -> bool:
    return status.state in {StageState.DONE, StageState.SKIPPED}


def _should_run_stage(
    spec: StageSpec,
    status: StageStatus,
    *,
    resume: bool,
    include_optional: bool,
) -> bool:
    if resume and _stage_is_complete(status):
        return False
    if spec.optional and not include_optional and status.state == StageState.OPTIONAL:
        return False
    return True


def _find_status(statuses: list[StageStatus], stage_id: str) -> StageStatus:
    for item in statuses:
        if item.stage.id == stage_id:
            return item
    msg = f"Stage {stage_id!r} not found in status list"
    raise KeyError(msg)


def _echo_stage_start(spec: StageSpec) -> None:
    tag = " (manual GUI)" if spec.manual else ""
    console.print(f"\n[bold cyan]→ {spec.title}[/bold cyan]{tag}  [dim]{spec.checkpoint_hint}[/dim]")


def _run_brain_stage(
    spec: StageSpec,
    config_path: Path,
    *,
    headless: bool,
    force_preprocess: bool,
) -> None:
    cfg = load_config(config_path)
    if spec.id == "preprocess":
        from lightsuite.preprocess.brain import preprocess_lightsheet_volume

        preprocess_lightsheet_volume(cfg, force=force_preprocess)
        return
    if spec.id == "check-orientation":
        from lightsuite.gui.orientation_brain import run_brain_orientation_check

        run_brain_orientation_check(cfg, config_path, headless=headless)
        return
    if spec.id == "align-slices":
        from lightsuite.gui.align_slices_brain import run_brain_align_slices

        run_brain_align_slices(cfg, headless=headless)
        return
    if spec.id == "init-registration":
        from lightsuite.registration.init_brain import initialize_brain_registration

        initialize_brain_registration(cfg)
        return
    if spec.id == "match-points":
        from lightsuite.gui.match_points_brain import run_brain_match_points

        run_brain_match_points(cfg, headless=headless)
        return
    if spec.id == "register":
        from lightsuite.registration.brain_register import run_brain_registration

        run_brain_registration(cfg, use_multistep=True)
        return
    if spec.id == "export":
        from lightsuite.export.brain_export import export_registered_brain_volumes

        export_registered_brain_volumes(cfg)
        return
    if spec.id == "import-annotations":
        from lightsuite.import_.brain_import import run_brain_import_annotations

        run_brain_import_annotations(cfg)
        return
    msg = f"Unsupported brain stage: {spec.id}"
    raise ValueError(msg)


def _run_spinal_stage(
    spec: StageSpec,
    config_path: Path,
    *,
    headless: bool,
) -> None:
    cfg = load_spinal_config(config_path)
    if spec.id == "check-orientation":
        from lightsuite.gui.orientation_cord import run_spinal_orientation

        run_spinal_orientation(cfg, headless=headless)
        return
    if spec.id == "preprocess":
        from lightsuite.preprocess.cord import preprocess_spinal_cord_sample

        preprocess_spinal_cord_sample(cfg)
        return
    if spec.id == "straighten":
        from lightsuite.gui.straighten_cord import run_spinal_straighten

        run_spinal_straighten(cfg, headless=headless)
        return
    if spec.id == "align-longitudinal":
        from lightsuite.gui.align_longitudinal_cord import run_spinal_align_longitudinal

        run_spinal_align_longitudinal(cfg, headless=headless)
        return
    if spec.id == "init-registration":
        from lightsuite.registration.init_cord import initialize_cord_registration

        initialize_cord_registration(cfg)
        return
    if spec.id == "match-points":
        from lightsuite.gui.match_points_cord import run_spinal_match_points

        run_spinal_match_points(cfg, headless=headless)
        return
    if spec.id == "register":
        from lightsuite.registration.cord_register import run_spinal_registration

        run_spinal_registration(cfg)
        return
    if spec.id == "export":
        from lightsuite.export.cord_export import export_registered_cord_volumes

        export_registered_cord_volumes(cfg)
        return
    if spec.id == "import-annotations":
        from lightsuite.import_.cord_import import run_cord_import_annotations

        run_cord_import_annotations(cfg)
        return
    if spec.id == "region-stats":
        from lightsuite.analysis.cord_runner import run_cord_region_stats

        run_cord_region_stats(cfg)
        return
    msg = f"Unsupported spinal stage: {spec.id}"
    raise ValueError(msg)


def _run_multires_stage(
    spec: StageSpec,
    config_path: Path,
    *,
    headless: bool,
) -> None:
    cfg = load_multires_config(config_path)
    if spec.id == "match-points":
        from lightsuite.gui.match_points_multires import run_multires_match_points

        run_multires_match_points(cfg, headless=headless)
        return
    if spec.id == "check-geometry":
        from lightsuite.multires.runner import check_multires_geometry

        check_multires_geometry(cfg)
        return
    if spec.id == "register":
        from lightsuite.multires.runner import run_multires_registration

        run_multires_registration(cfg)
        return
    if spec.id == "import-annotations":
        from lightsuite.multires.import_annotations import run_multires_import_annotations

        run_multires_import_annotations(cfg)
        return
    msg = f"Unsupported multires stage: {spec.id}"
    raise ValueError(msg)


def _execute_pipeline(
    *,
    workflow: str,
    config_path: Path,
    specs: list[StageSpec],
    statuses: list[StageStatus],
    run_stage: Callable[..., None],
    from_stage: str | None,
    through_stage: str | None,
    resume: bool,
    include_optional: bool,
    headless: bool,
    force_preprocess: bool,
) -> None:
    selected = slice_stages(specs, from_stage=from_stage, through_stage=through_stage)
    status_by_id = {item.stage.id: item for item in statuses}
    ran = 0
    for spec in selected:
        status = status_by_id[spec.id]
        if not _should_run_stage(
            spec,
            status,
            resume=resume,
            include_optional=include_optional,
        ):
            console.print(
                f"[dim]Skipping {spec.title} ({status.state.value}) — {status.detail}[/dim]"
            )
            continue
        _echo_stage_start(spec)
        try:
            if workflow == "brain":
                run_stage(spec, config_path, headless=headless, force_preprocess=force_preprocess)
            else:
                run_stage(spec, config_path, headless=headless)
        except Exception as exc:
            console.print(f"[bold red]Stage failed:[/bold red] {spec.id} — {exc}")
            console.print(
                f"Re-run manually: lightsuite {workflow} {spec.id} -c {config_path}"
            )
            raise typer.Exit(code=1) from exc
        ran += 1
    if ran == 0:
        console.print("[yellow]No stages executed (all complete or optional).[/yellow]")
    else:
        console.print(f"\n[green]Completed {ran} stage(s).[/green]")


def run_brain_pipeline(
    config_path: str | Path,
    *,
    from_stage: str | None = None,
    through_stage: str | None = None,
    resume: bool = False,
    include_optional: bool = False,
    headless: bool = False,
    force_preprocess: bool = False,
) -> None:
    path = Path(config_path).expanduser().resolve()
    cfg = load_config(path)
    specs = brain_stage_specs(cfg)
    statuses = brain_stage_statuses(cfg)
    _execute_pipeline(
        workflow="brain",
        config_path=path,
        specs=specs,
        statuses=statuses,
        run_stage=_run_brain_stage,
        from_stage=from_stage,
        through_stage=through_stage,
        resume=resume,
        include_optional=include_optional,
        headless=headless,
        force_preprocess=force_preprocess,
    )


def run_spinal_pipeline(
    config_path: str | Path,
    *,
    from_stage: str | None = None,
    through_stage: str | None = None,
    resume: bool = False,
    include_optional: bool = True,
    headless: bool = False,
) -> None:
    path = Path(config_path).expanduser().resolve()
    cfg = load_spinal_config(path)
    specs = spinal_stage_specs(cfg)
    statuses = spinal_stage_statuses(cfg)
    _execute_pipeline(
        workflow="spinal",
        config_path=path,
        specs=specs,
        statuses=statuses,
        run_stage=_run_spinal_stage,
        from_stage=from_stage,
        through_stage=through_stage,
        resume=resume,
        include_optional=include_optional,
        headless=headless,
        force_preprocess=False,
    )


def run_multires_pipeline(
    config_path: str | Path,
    *,
    from_stage: str | None = None,
    through_stage: str | None = None,
    resume: bool = False,
    include_optional: bool = True,
    headless: bool = False,
) -> None:
    path = Path(config_path).expanduser().resolve()
    cfg = load_multires_config(path)
    specs = multires_stage_specs(cfg)
    statuses = multires_stage_statuses(cfg)
    _execute_pipeline(
        workflow="multires",
        config_path=path,
        specs=specs,
        statuses=statuses,
        run_stage=_run_multires_stage,
        from_stage=from_stage,
        through_stage=through_stage,
        resume=resume,
        include_optional=include_optional,
        headless=headless,
        force_preprocess=False,
    )
