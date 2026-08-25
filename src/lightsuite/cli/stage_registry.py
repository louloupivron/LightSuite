"""Unified stage registry for brain, spinal, and multires workflows."""

from __future__ import annotations

import threading
from collections.abc import Callable
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Any

from lightsuite.cli.stages import (
    StageSpec,
    StageStatus,
    brain_stage_specs,
    brain_stage_statuses,
    multires_stage_specs,
    multires_stage_statuses,
    spinal_stage_specs,
    spinal_stage_statuses,
)
from lightsuite.config.loader import load_config, load_multires_config, load_spinal_config


class StageKind(str, Enum):
    """Whether a stage runs headlessly or opens an interactive viewer."""

    AUTO = "auto"
    INTERACTIVE = "interactive"


@dataclass(frozen=True)
class StageContext:
    """Runtime options passed to every stage runner."""

    config_path: Path
    headless: bool = False
    force_preprocess: bool = False
    export_spaces: list[str] | None = None
    on_log: Callable[[str], None] | None = None
    cancel_event: threading.Event | None = None
    on_config_file_changed: Callable[[], None] | None = None


StageRunner = Callable[[Any, StageContext], Any]
ConfigLoader = Callable[[str | Path], Any]
StageSpecsFn = Callable[[Any], list[StageSpec]]
StageStatusesFn = Callable[[Any], list[StageStatus]]


@dataclass(frozen=True)
class WorkflowSpec:
    """Metadata and helpers for one pipeline workflow."""

    name: str
    load_config: ConfigLoader
    stage_specs: StageSpecsFn
    stage_statuses: StageStatusesFn
    runners: dict[str, StageRunner]


def _brain_preprocess(config: Any, ctx: StageContext) -> Any:
    from lightsuite.preprocess.brain import preprocess_lightsheet_volume

    return preprocess_lightsheet_volume(config, force=ctx.force_preprocess)


def _brain_check_orientation(config: Any, ctx: StageContext) -> Any:
    from lightsuite.gui.orientation_brain import run_brain_orientation_check

    return run_brain_orientation_check(config, ctx.config_path, headless=ctx.headless)


def _brain_align_slices(config: Any, ctx: StageContext) -> Any:
    from lightsuite.gui.align_slices_brain import run_brain_align_slices

    return run_brain_align_slices(config, headless=ctx.headless)


def _brain_init_registration(config: Any, ctx: StageContext) -> Any:
    from lightsuite.registration.init_brain import initialize_brain_registration

    return initialize_brain_registration(config)


def _brain_match_points(config: Any, ctx: StageContext) -> Any:
    from lightsuite.gui.match_points_brain import run_brain_match_points

    return run_brain_match_points(config, headless=ctx.headless)


def _brain_register(config: Any, ctx: StageContext) -> Any:
    from lightsuite.registration.brain_register import run_brain_registration

    return run_brain_registration(config, use_multistep=True)


def _brain_export(config: Any, ctx: StageContext) -> Any:
    from lightsuite.export.brain_export import export_registered_brain_volumes

    return export_registered_brain_volumes(config, spaces=ctx.export_spaces)


def _brain_import_annotations(config: Any, ctx: StageContext) -> Any:
    from lightsuite.import_.brain_import import run_brain_import_annotations

    return run_brain_import_annotations(config)


def _brain_view_registration(config: Any, ctx: StageContext) -> Any:
    from lightsuite.gui.view_registration_brain import run_brain_view_registration

    return run_brain_view_registration(config, space="sample", headless=ctx.headless)


def _spinal_check_orientation(config: Any, ctx: StageContext) -> Any:
    from lightsuite.gui.orientation_cord import run_spinal_orientation

    return run_spinal_orientation(config, headless=ctx.headless)


def _spinal_preprocess(config: Any, ctx: StageContext) -> Any:
    from lightsuite.preprocess.cord import preprocess_spinal_cord_sample

    return preprocess_spinal_cord_sample(config, headless=ctx.headless)


def _spinal_straighten(config: Any, ctx: StageContext) -> Any:
    from lightsuite.gui.straighten_cord import run_spinal_straighten

    return run_spinal_straighten(config, headless=ctx.headless)


def _spinal_align_longitudinal(config: Any, ctx: StageContext) -> Any:
    from lightsuite.gui.align_longitudinal_cord import run_spinal_align_longitudinal

    return run_spinal_align_longitudinal(config, headless=ctx.headless)


def _spinal_init_registration(config: Any, ctx: StageContext) -> Any:
    from lightsuite.registration.init_cord import initialize_cord_registration

    return initialize_cord_registration(config)


def _spinal_match_points(config: Any, ctx: StageContext) -> Any:
    from lightsuite.gui.match_points_cord import run_spinal_match_points

    return run_spinal_match_points(config, headless=ctx.headless)


def _spinal_register(config: Any, ctx: StageContext) -> Any:
    from lightsuite.registration.cord_register import run_spinal_registration

    return run_spinal_registration(config)


def _spinal_export(config: Any, ctx: StageContext) -> Any:
    from lightsuite.export.cord_export import export_registered_cord_volumes

    return export_registered_cord_volumes(config, spaces=ctx.export_spaces)


def _spinal_import_annotations(config: Any, ctx: StageContext) -> Any:
    from lightsuite.import_.cord_import import run_cord_import_annotations

    return run_cord_import_annotations(config)


def _spinal_view_registration(config: Any, ctx: StageContext) -> Any:
    from lightsuite.gui.view_registered_cord import run_spinal_view_registration

    return run_spinal_view_registration(config, space="sample", headless=ctx.headless)


def _spinal_region_stats(config: Any, ctx: StageContext) -> Any:
    from lightsuite.analysis.cord_runner import run_cord_region_stats

    return run_cord_region_stats(config)


def _spinal_plot_stats(config: Any, ctx: StageContext) -> Any:
    from lightsuite.gui.spinal_stats_plots import run_spinal_stats_plots

    return run_spinal_stats_plots(config, headless=ctx.headless)


def _multires_match_points(config: Any, ctx: StageContext) -> Any:
    from lightsuite.gui.match_points_multires import run_multires_match_points

    return run_multires_match_points(config, headless=ctx.headless)


def _multires_inspect_geometry(config: Any, ctx: StageContext) -> Any:
    from lightsuite.gui.inspect_geometry_multires import run_multires_inspect_geometry

    return run_multires_inspect_geometry(
        config,
        config_path=ctx.config_path,
        headless=ctx.headless,
    )


def _multires_check_geometry(config: Any, ctx: StageContext) -> Any:
    from lightsuite.multires.runner import check_multires_geometry

    return check_multires_geometry(config)


def _multires_register(config: Any, ctx: StageContext) -> Any:
    from lightsuite.multires.runner import run_multires_registration

    return run_multires_registration(config)


def _multires_inspect_registration(config: Any, ctx: StageContext) -> Any:
    from lightsuite.gui.inspect_registration_multires import run_multires_inspect_registration

    return run_multires_inspect_registration(config, headless=ctx.headless)


def _multires_import_annotations(config: Any, ctx: StageContext) -> Any:
    from lightsuite.multires.import_annotations import run_multires_import_annotations

    return run_multires_import_annotations(config)


_BRAIN_RUNNERS: dict[str, StageRunner] = {
    "preprocess": _brain_preprocess,
    "check-orientation": _brain_check_orientation,
    "align-slices": _brain_align_slices,
    "init-registration": _brain_init_registration,
    "match-points": _brain_match_points,
    "register": _brain_register,
    "export": _brain_export,
    "view-registration": _brain_view_registration,
    "import-annotations": _brain_import_annotations,
}

_SPINAL_RUNNERS: dict[str, StageRunner] = {
    "check-orientation": _spinal_check_orientation,
    "preprocess": _spinal_preprocess,
    "straighten": _spinal_straighten,
    "align-longitudinal": _spinal_align_longitudinal,
    "init-registration": _spinal_init_registration,
    "match-points": _spinal_match_points,
    "register": _spinal_register,
    "export": _spinal_export,
    "view-registration": _spinal_view_registration,
    "import-annotations": _spinal_import_annotations,
    "region-stats": _spinal_region_stats,
    "plot-stats": _spinal_plot_stats,
}

_MULTIRES_RUNNERS: dict[str, StageRunner] = {
    "inspect-geometry": _multires_inspect_geometry,
    "match-points": _multires_match_points,
    "check-geometry": _multires_check_geometry,
    "register": _multires_register,
    "inspect-registration": _multires_inspect_registration,
    "import-annotations": _multires_import_annotations,
}

WORKFLOWS: dict[str, WorkflowSpec] = {
    "brain": WorkflowSpec(
        name="brain",
        load_config=load_config,
        stage_specs=brain_stage_specs,
        stage_statuses=brain_stage_statuses,
        runners=_BRAIN_RUNNERS,
    ),
    "spinal": WorkflowSpec(
        name="spinal",
        load_config=load_spinal_config,
        stage_specs=spinal_stage_specs,
        stage_statuses=spinal_stage_statuses,
        runners=_SPINAL_RUNNERS,
    ),
    "multires": WorkflowSpec(
        name="multires",
        load_config=load_multires_config,
        stage_specs=multires_stage_specs,
        stage_statuses=multires_stage_statuses,
        runners=_MULTIRES_RUNNERS,
    ),
}


def get_workflow(workflow: str) -> WorkflowSpec:
    """Return workflow metadata; raises KeyError if unknown."""
    try:
        return WORKFLOWS[workflow]
    except KeyError as exc:
        known = ", ".join(sorted(WORKFLOWS))
        msg = f"Unknown workflow {workflow!r}; choose from: {known}"
        raise KeyError(msg) from exc


def stage_kind(spec: StageSpec) -> StageKind:
    """Classify a stage as automated or interactive."""
    return StageKind.INTERACTIVE if spec.manual else StageKind.AUTO


def run_stage(workflow: str, stage_id: str, config: Any, ctx: StageContext) -> Any:
    """Execute one registered stage for the given workflow."""
    spec = get_workflow(workflow)
    try:
        runner = spec.runners[stage_id]
    except KeyError as exc:
        known = ", ".join(sorted(spec.runners))
        msg = f"Unsupported {workflow} stage: {stage_id!r}; choose from: {known}"
        raise ValueError(msg) from exc
    return runner(config, ctx)


def validate_registry() -> None:
    """Ensure runner tables are internally consistent (dev/test helper)."""
    for workflow_name, workflow in WORKFLOWS.items():
        if not workflow.runners:
            msg = f"{workflow_name}: empty runner table"
            raise RuntimeError(msg)
        for stage_id, runner in workflow.runners.items():
            if not callable(runner):
                msg = f"{workflow_name}: runner for {stage_id!r} is not callable"
                raise RuntimeError(msg)


def runner_stage_ids(workflow: str) -> frozenset[str]:
    """Return all stage ids with registered runners for a workflow."""
    return frozenset(get_workflow(workflow).runners)
