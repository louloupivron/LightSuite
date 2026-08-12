"""Run multiple pipeline stages in order with resume support."""

from __future__ import annotations

from pathlib import Path

import typer

from lightsuite.cli.stage_registry import StageContext, get_workflow, run_stage
from lightsuite.cli.stages import (
    StageSpec,
    StageState,
    StageStatus,
    slice_stages,
)
from lightsuite.reporter import ConsoleReporter, Reporter


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


def _execute_pipeline(
    *,
    workflow: str,
    config_path: Path,
    specs: list[StageSpec],
    statuses: list[StageStatus],
    config: object,
    from_stage: str | None,
    through_stage: str | None,
    resume: bool,
    include_optional: bool,
    headless: bool,
    force_preprocess: bool,
    reporter: Reporter | None = None,
) -> None:
    sink = reporter or ConsoleReporter()
    selected = slice_stages(specs, from_stage=from_stage, through_stage=through_stage)
    status_by_id = {item.stage.id: item for item in statuses}
    ctx = StageContext(
        config_path=config_path,
        headless=headless,
        force_preprocess=force_preprocess,
    )
    ran = 0
    for spec in selected:
        status = status_by_id[spec.id]
        if not _should_run_stage(
            spec,
            status,
            resume=resume,
            include_optional=include_optional,
        ):
            sink.stage_skip(spec.title, status.state.value, status.detail)
            continue
        sink.stage_start(spec.title, manual=spec.manual, checkpoint_hint=spec.checkpoint_hint)
        try:
            run_stage(workflow, spec.id, config, ctx)
        except Exception as exc:
            sink.stage_failed(spec.id, exc)
            sink.rerun_hint(workflow, spec.id, config_path)
            raise typer.Exit(code=1) from exc
        ran += 1
    if ran == 0:
        sink.pipeline_no_stages()
    else:
        sink.pipeline_complete(ran)


def run_brain_pipeline(
    config_path: str | Path,
    *,
    from_stage: str | None = None,
    through_stage: str | None = None,
    resume: bool = False,
    include_optional: bool = False,
    headless: bool = False,
    force_preprocess: bool = False,
    reporter: Reporter | None = None,
) -> None:
    path = Path(config_path).expanduser().resolve()
    workflow = get_workflow("brain")
    cfg = workflow.load_config(path)
    _execute_pipeline(
        workflow="brain",
        config_path=path,
        specs=workflow.stage_specs(cfg),
        statuses=workflow.stage_statuses(cfg),
        config=cfg,
        from_stage=from_stage,
        through_stage=through_stage,
        resume=resume,
        include_optional=include_optional,
        headless=headless,
        force_preprocess=force_preprocess,
        reporter=reporter,
    )


def run_spinal_pipeline(
    config_path: str | Path,
    *,
    from_stage: str | None = None,
    through_stage: str | None = None,
    resume: bool = False,
    include_optional: bool = True,
    headless: bool = False,
    reporter: Reporter | None = None,
) -> None:
    path = Path(config_path).expanduser().resolve()
    workflow = get_workflow("spinal")
    cfg = workflow.load_config(path)
    _execute_pipeline(
        workflow="spinal",
        config_path=path,
        specs=workflow.stage_specs(cfg),
        statuses=workflow.stage_statuses(cfg),
        config=cfg,
        from_stage=from_stage,
        through_stage=through_stage,
        resume=resume,
        include_optional=include_optional,
        headless=headless,
        force_preprocess=False,
        reporter=reporter,
    )


def run_multires_pipeline(
    config_path: str | Path,
    *,
    from_stage: str | None = None,
    through_stage: str | None = None,
    resume: bool = False,
    include_optional: bool = True,
    headless: bool = False,
    reporter: Reporter | None = None,
) -> None:
    path = Path(config_path).expanduser().resolve()
    workflow = get_workflow("multires")
    cfg = workflow.load_config(path)
    _execute_pipeline(
        workflow="multires",
        config_path=path,
        specs=workflow.stage_specs(cfg),
        statuses=workflow.stage_statuses(cfg),
        config=cfg,
        from_stage=from_stage,
        through_stage=through_stage,
        resume=resume,
        include_optional=include_optional,
        headless=headless,
        force_preprocess=False,
        reporter=reporter,
    )
