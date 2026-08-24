"""Map pipeline stages to napari attach functions for the unified GUI."""

from __future__ import annotations

from collections.abc import Callable
from typing import Any

from lightsuite.cli.stage_registry import StageContext
from lightsuite.gui.stage_controller import StageController

AttachFactory = Callable[[Any, Any, StageContext], StageController]


def _brain_orientation(viewer: Any, config: Any, ctx: StageContext) -> StageController:
    from lightsuite.gui.orientation_brain import attach_brain_orientation_check

    return attach_brain_orientation_check(viewer, config, ctx.config_path)


def _brain_align_slices(viewer: Any, config: Any, ctx: StageContext) -> StageController:
    from lightsuite.gui.align_slices_brain import attach_brain_align_slices

    return attach_brain_align_slices(viewer, config)


def _brain_match_points(viewer: Any, config: Any, ctx: StageContext) -> StageController:
    from lightsuite.gui.match_points_brain import attach_brain_match_points

    return attach_brain_match_points(viewer, config)


def _brain_view_registration(viewer: Any, config: Any, ctx: StageContext) -> StageController:
    from lightsuite.gui.view_registration_brain import attach_brain_view_registration

    return attach_brain_view_registration(viewer, config, space="sample")


def _spinal_view_registration(viewer: Any, config: Any, ctx: StageContext) -> StageController:
    from lightsuite.gui.view_registered_cord import attach_spinal_view_registration

    return attach_spinal_view_registration(viewer, config, space="sample")


def _spinal_orientation(viewer: Any, config: Any, ctx: StageContext) -> StageController:
    from lightsuite.gui.orientation_cord import attach_spinal_orientation

    return attach_spinal_orientation(viewer, config, ctx=ctx)


def _spinal_straighten(viewer: Any, config: Any, ctx: StageContext) -> StageController:
    from lightsuite.gui.straighten_cord import attach_spinal_straighten

    return attach_spinal_straighten(viewer, config)


def _spinal_align_longitudinal(viewer: Any, config: Any, ctx: StageContext) -> StageController:
    from lightsuite.gui.align_longitudinal_cord import attach_spinal_align_longitudinal

    return attach_spinal_align_longitudinal(viewer, config)


def _spinal_match_points(viewer: Any, config: Any, ctx: StageContext) -> StageController:
    from lightsuite.gui.match_points_cord import attach_spinal_match_points

    return attach_spinal_match_points(viewer, config)


def _spinal_stats_plots(viewer: Any, config: Any, ctx: StageContext) -> StageController:
    from lightsuite.gui.spinal_stats_plots import attach_spinal_stats_plots

    return attach_spinal_stats_plots(viewer, config)


def _multires_match_points(viewer: Any, config: Any, ctx: StageContext) -> StageController:
    from lightsuite.gui.match_points_multires import attach_multires_match_points

    return attach_multires_match_points(viewer, config)


def _multires_inspect_geometry_attach(viewer: Any, config: Any, ctx: StageContext) -> StageController:
    from lightsuite.gui.inspect_geometry_multires import attach_multires_inspect_geometry

    return attach_multires_inspect_geometry(viewer, config, config_path=ctx.config_path)


def _multires_inspect_registration_attach(viewer: Any, config: Any, ctx: StageContext) -> StageController:
    from lightsuite.gui.inspect_registration_multires import (
        attach_multires_inspect_registration,
        discover_multires_registration_inspect_paths,
    )

    paths = discover_multires_registration_inspect_paths(config)
    return attach_multires_inspect_registration(viewer, config, paths)


STAGE_ATTACH: dict[tuple[str, str], AttachFactory] = {
    ("brain", "check-orientation"): _brain_orientation,
    ("brain", "align-slices"): _brain_align_slices,
    ("brain", "match-points"): _brain_match_points,
    ("brain", "view-registration"): _brain_view_registration,
    ("spinal", "check-orientation"): _spinal_orientation,
    ("spinal", "view-registration"): _spinal_view_registration,
    ("spinal", "straighten"): _spinal_straighten,
    ("spinal", "align-longitudinal"): _spinal_align_longitudinal,
    ("spinal", "match-points"): _spinal_match_points,
    ("spinal", "plot-stats"): _spinal_stats_plots,
    ("multires", "match-points"): _multires_match_points,
    ("multires", "inspect-geometry"): _multires_inspect_geometry_attach,
    ("multires", "inspect-registration"): _multires_inspect_registration_attach,
}


def get_stage_attach(workflow: str, stage_id: str) -> AttachFactory | None:
    """Return an attach factory for an interactive stage, or None if headless-only."""
    return STAGE_ATTACH.get((workflow, stage_id))


def attach_stage(
    workflow: str,
    stage_id: str,
    viewer: Any,
    config: Any,
    ctx: StageContext,
) -> StageController:
    """Attach an interactive stage to a napari viewer."""
    factory = get_stage_attach(workflow, stage_id)
    if factory is None:
        msg = f"No attach function for {workflow}/{stage_id!r}"
        raise ValueError(msg)
    return factory(viewer, config, ctx)
