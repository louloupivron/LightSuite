"""Pipeline stage graphs and checkpoint status for orchestration."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Any, Protocol

from lightsuite.config.models import BrainPipelineConfig, SpinalCordPipelineConfig
from lightsuite.gui.control_points import ControlPointSession, default_session_path
from lightsuite.gui.slice_correspondence import default_correspondence_path
from lightsuite.multires.checkpoint import MultiresRegOptsCheckpoint, multires_checkpoint_path
from lightsuite.preprocess.checkpoint import RegOptsCheckpoint
from lightsuite.preprocess.cord_checkpoint import CordRegOptsCheckpoint
from lightsuite.registration.cord_longitudinal import default_longitudinal_correspondence_path
from lightsuite.registration.cord_orientation import CORD_ORIENTATION_FILENAME
from lightsuite.registration.orientation import orientation_path


class StageState(str, Enum):
    PENDING = "pending"
    DONE = "done"
    SKIPPED = "skipped"
    OPTIONAL = "optional"


@dataclass(frozen=True)
class StageSpec:
    """One runnable pipeline step."""

    id: str
    title: str
    checkpoint_hint: str
    optional: bool = False
    manual: bool = False


@dataclass(frozen=True)
class StageStatus:
    stage: StageSpec
    state: StageState
    detail: str = ""


class _Config(Protocol):
    sample: Any


def _load_regopts(save_path: Path) -> RegOptsCheckpoint | None:
    path = save_path / "regopts.json"
    if not path.is_file():
        return None
    return RegOptsCheckpoint.load(path)


def _load_cord_regopts(save_path: Path) -> CordRegOptsCheckpoint | None:
    path = save_path / "regopts.json"
    if not path.is_file():
        return None
    return CordRegOptsCheckpoint.load(path)


def _load_multires_checkpoint(save_path: Path) -> MultiresRegOptsCheckpoint | None:
    path = multires_checkpoint_path(save_path)
    if not path.is_file():
        return None
    return MultiresRegOptsCheckpoint.load(path)


def _has_import_annotations(config: _Config) -> bool:
    import_cfg = getattr(config, "import_", None)
    if import_cfg is None:
        return False
    annotations = getattr(import_cfg, "annotations", None)
    return bool(annotations)


@dataclass
class _PreviewPipelineConfig:
    """Minimal config stand-in for GUI stage-list preview before save."""

    import_: Any = None
    multires: Any = None


_PREVIEW_CONFIG = _PreviewPipelineConfig()
_PREVIEW_STAGE_DETAIL = "Save a valid config file to run this stage"


def preview_stage_statuses(workflow: str) -> list[StageStatus]:
    """Return the default stage checklist for a workflow (all pending/optional)."""
    workflow_key = workflow.strip().lower()
    if workflow_key == "brain":
        specs = brain_stage_specs(_PREVIEW_CONFIG)  # type: ignore[arg-type]
    elif workflow_key == "spinal":
        specs = spinal_stage_specs(_PREVIEW_CONFIG)  # type: ignore[arg-type]
    elif workflow_key == "multires":
        specs = multires_stage_specs(_PREVIEW_CONFIG)
    else:
        msg = f"Unknown workflow {workflow!r} for stage preview."
        raise ValueError(msg)

    statuses: list[StageStatus] = []
    for spec in specs:
        state = StageState.OPTIONAL if spec.optional else StageState.PENDING
        statuses.append(StageStatus(spec, state, _PREVIEW_STAGE_DETAIL))
    return statuses


def brain_stage_specs(config: BrainPipelineConfig) -> list[StageSpec]:
    stages = [
        StageSpec("preprocess", "Preprocess", "regopts.json → regvolpath"),
        StageSpec(
            "check-orientation",
            "Check orientation",
            "brain_orientation.txt or registration.orientation",
            manual=True,
        ),
        StageSpec(
            "align-slices",
            "Align slices",
            "slice_correspondence.json",
            optional=True,
            manual=True,
        ),
        StageSpec(
            "init-registration",
            "Init registration",
            "regopts.json → original_trans",
        ),
        StageSpec(
            "match-points",
            "Match control points",
            "atlas2histology_tform.json",
            optional=True,
            manual=True,
        ),
        StageSpec("register", "Register", "transform_params.json"),
        StageSpec("export", "Export", "volume_registered/"),
        StageSpec(
            "view-registration",
            "View registration",
            "registration review (Napari)",
            optional=True,
            manual=True,
        ),
    ]
    if _has_import_annotations(config):
        stages.append(
            StageSpec(
                "import-annotations",
                "Import annotations",
                "volume_registered/*_atlas_coords.csv",
                optional=True,
            )
        )
    return stages


def spinal_stage_specs(config: SpinalCordPipelineConfig) -> list[StageSpec]:
    stages: list[StageSpec] = [
        StageSpec("preprocess", "Preprocess", "regopts.json"),
        StageSpec(
            "straighten",
            "Straighten",
            "spinal_alignment_opt.json",
            manual=True,
        ),
        StageSpec(
            "align-longitudinal",
            "Align longitudinal",
            "longitudinal_correspondence.json",
            manual=True,
        ),
        StageSpec("init-registration", "Init registration", "regopts.json → slicetforms"),
        StageSpec(
            "match-points",
            "Match control points",
            "atlas2histology_tform.json",
            optional=True,
            manual=True,
        ),
        StageSpec("register", "Register", "transform_params.json"),
        StageSpec("export", "Export", "volume_registered/ (+ region stats)"),
        StageSpec(
            "view-registration",
            "View registration",
            "channels + annotation + imports (Napari)",
            optional=True,
            manual=True,
        ),
    ]
    if _has_import_annotations(config):
        stages.append(
            StageSpec(
                "import-annotations",
                "Import annotations",
                "volume_registered/import_annotations_summary.json",
                optional=True,
            )
        )
    stages.append(
        StageSpec(
            "plot-stats",
            "Stats / plots",
            "plots/cord_heatmap_*.png",
            optional=True,
            manual=True,
        )
    )
    return stages


def multires_stage_specs(config: _Config) -> list[StageSpec]:
    stages: list[StageSpec] = []
    multires = getattr(config, "multires", None)
    landmarks = getattr(multires, "landmarks", None) if multires else None
    if landmarks is not None and getattr(landmarks, "session_path", None) is not False:
        stages.append(
            StageSpec(
                "match-points",
                "Match points",
                "multires landmark session",
                optional=True,
                manual=True,
            ),
        )
    stages.append(
        StageSpec(
            "inspect-geometry",
            "Inspect geometry",
            "multires.mesospim_geometry lateral_flip",
            optional=True,
            manual=True,
        ),
    )
    stages.extend(
        [
            StageSpec("check-geometry", "Check geometry", "geometry/ overlap QA"),
            StageSpec("register", "Register", "multires_regopts.json → transform_paths"),
            StageSpec(
                "inspect-registration",
                "Inspect registration",
                "registered ROI vs overview (Napari)",
                optional=True,
                manual=True,
            ),
        ]
    )
    if _has_import_annotations(config):
        stages.append(
            StageSpec(
                "import-annotations",
                "Import annotations",
                "annotations_in_overview/",
                optional=True,
            )
        )
    return stages


def _brain_stage_done(stage_id: str, save_path: Path, config: BrainPipelineConfig) -> tuple[bool, str]:
    regopts = _load_regopts(save_path)
    if stage_id == "preprocess":
        if regopts is None or not regopts.regvolpath:
            return False, "missing regopts.json or regvolpath"
        regvol = Path(regopts.regvolpath)
        return regvol.is_file(), str(regvol)
    if stage_id == "check-orientation":
        if config.registration.orientation is not None:
            return True, "registration.orientation in YAML"
        orient = orientation_path(save_path)
        return orient.is_file(), str(orient)
    if stage_id == "align-slices":
        path = default_correspondence_path(save_path)
        if not path.is_file():
            return False, str(path)
        from lightsuite.gui.brain_data import load_slice_correspondence

        corr = load_slice_correspondence(save_path)
        if corr is None or not corr.has_confirmed_anchors():
            return False, "no confirmed anchors"
        return True, str(path)
    if stage_id == "init-registration":
        if regopts is None or regopts.original_trans is None:
            return False, "missing original_trans"
        return True, "original_trans present"
    if stage_id == "match-points":
        session_path = default_session_path(save_path)
        if not session_path.is_file():
            return False, str(session_path)
        session = ControlPointSession.load(session_path)
        has_pairs = any(session.histology_control_points) and any(session.atlas_control_points)
        return has_pairs, str(session_path)
    if stage_id == "register":
        tpath = save_path / "transform_params.json"
        return tpath.is_file(), str(tpath)
    if stage_id == "export":
        vr = save_path / "volume_registered"
        if not vr.is_dir():
            return False, str(vr)
        tiffs = list(vr.glob("*.tif")) + list(vr.glob("*.tiff"))
        return bool(tiffs), f"{len(tiffs)} TIFF(s) in volume_registered/"
    if stage_id == "view-registration":
        vr = save_path / "volume_registered"
        if not vr.is_dir():
            return False, str(vr)
        tiffs = list(vr.glob("*.tif")) + list(vr.glob("*.tiff"))
        return bool(tiffs), "open Napari after export"
    if stage_id == "import-annotations":
        summary = save_path / "volume_registered" / "import_annotations_summary.json"
        if summary.is_file():
            return True, str(summary)
        csvs = list((save_path / "volume_registered").glob("*_atlas_coords.csv"))
        return bool(csvs), f"{len(csvs)} atlas-space CSV layer(s)"
    return False, "unknown stage"


def _spinal_stage_done(
    stage_id: str,
    save_path: Path,
    config: SpinalCordPipelineConfig,
) -> tuple[bool, str]:
    regopts = _load_cord_regopts(save_path)
    if stage_id == "check-orientation":
        if config.registration.longitudinal_direction is not None:
            return True, "registration.longitudinal_direction in YAML"
        path = save_path / CORD_ORIENTATION_FILENAME
        return path.is_file(), str(path)
    if stage_id == "preprocess":
        return regopts is not None, "regopts.json"
    if stage_id == "straighten":
        path = save_path / "spinal_alignment_opt.json"
        return path.is_file(), str(path)
    if stage_id == "align-longitudinal":
        path = default_longitudinal_correspondence_path(save_path)
        if not path.is_file():
            return False, str(path)
        from lightsuite.gui.slice_correspondence import SliceCorrespondence

        corr = SliceCorrespondence.load(path)
        return corr.has_confirmed_anchors(), str(path)
    if stage_id == "init-registration":
        tpath = save_path / "transform_params.json"
        if tpath.is_file():
            return True, "transform_params.json (post-init)"
        if regopts is None:
            return False, "missing regopts.json"
        return (
            regopts.affine_atlas_to_samp is not None,
            "affine_atlas_to_samp in regopts.json",
        )
    if stage_id == "match-points":
        from lightsuite.gui.cord_data import default_cord_session_path

        session_path = default_cord_session_path(save_path)
        if not session_path.is_file():
            return False, str(session_path)
        session = ControlPointSession.load(session_path)
        has_pairs = any(session.histology_control_points) and any(session.atlas_control_points)
        return has_pairs, str(session_path)
    if stage_id == "register":
        tpath = save_path / "transform_params.json"
        return tpath.is_file(), str(tpath)
    if stage_id == "export":
        vr = save_path / "volume_registered"
        if not vr.is_dir() or not any(vr.iterdir()):
            return False, str(vr)
        if config.analysis.parcellate_intensities:
            stats = list(vr.glob("region_stats*.csv")) + list(vr.glob("chan*_region_stats.csv"))
            if not stats:
                return False, "volumes exported; intensity region stats missing"
        return True, str(vr)
    if stage_id == "view-registration":
        vr = save_path / "volume_registered"
        if not vr.is_dir():
            return False, str(vr)
        tiffs = list(vr.glob("*.tif")) + list(vr.glob("*.tiff"))
        npz = list(vr.glob("*_atlas_coords.npz")) + list(vr.glob("*_sample_coords.npz"))
        if tiffs or npz:
            detail = (
                f"{len(tiffs)} TIFF(s), {len(npz)} import layer(s)"
                if npz
                else f"{len(tiffs)} TIFF(s)"
            )
            return True, detail
        return False, "run export and/or import-annotations first"
    if stage_id == "import-annotations":
        summary = save_path / "volume_registered" / "import_annotations_summary.json"
        return summary.is_file(), str(summary)
    if stage_id == "plot-stats":
        plots = save_path / "plots"
        stats = save_path / "volume_registered" / "region_stats.csv"
        if plots.is_dir() and any(plots.glob("cord_heatmap_*.png")):
            return True, str(plots)
        return stats.is_file(), "open Stats / plots after export"
    return False, "unknown stage"


def _multires_geometry_configured(config: Any) -> tuple[bool, str]:
    multires = config.multires
    geometry = multires.mesospim_geometry
    if geometry is None:
        return False, "multires.mesospim_geometry not set"
    overview = geometry.overview
    roi = geometry.roi
    if overview is None or roi is None:
        return False, "mesospim_geometry overview/roi incomplete"
    if overview.lateral_flip is None or roi.lateral_flip is None:
        return False, "lateral_flip not configured"
    return True, f"lateral_flip overview={list(overview.lateral_flip)} roi={list(roi.lateral_flip)}"


def _multires_stage_done(stage_id: str, save_path: Path, config: Any) -> tuple[bool, str]:
    checkpoint = _load_multires_checkpoint(save_path)
    if stage_id == "inspect-geometry":
        return _multires_geometry_configured(config)
    if stage_id == "match-points":
        multires = config.multires
        if multires.landmarks and multires.landmarks.session_path:
            path = Path(multires.landmarks.session_path).expanduser()
        else:
            path = save_path / f"multires_landmarks_{multires.pair_label}.json"
        return path.is_file(), str(path)
    if stage_id == "check-geometry":
        if checkpoint and checkpoint.geometry_report_paths:
            return True, "geometry reports in checkpoint"
        geom = save_path / "geometry"
        return geom.is_dir() and any(geom.rglob("*")), str(geom)
    if stage_id == "register":
        if checkpoint is None:
            return False, "missing multires_regopts.json"
        if not checkpoint.transform_paths:
            return False, "no transform_paths"
        return True, str(multires_checkpoint_path(save_path))
    if stage_id == "inspect-registration":
        if checkpoint is None or not checkpoint.registered_roi_path:
            return False, "run multires register first"
        path = Path(checkpoint.registered_roi_path).expanduser()
        return path.is_file(), str(path)
    if stage_id == "import-annotations":
        out = save_path / "annotations_in_overview"
        if not out.is_dir():
            return False, str(out)
        artifacts = list(out.glob("*.csv")) + list(out.glob("*.tif")) + list(out.glob("*.tiff"))
        return bool(artifacts), f"{len(artifacts)} file(s)"
    return False, "unknown stage"


def evaluate_stage_statuses(
    specs: list[StageSpec],
    *,
    done_fn,
    save_path: Path,
    config: Any,
) -> list[StageStatus]:
    statuses: list[StageStatus] = []
    for spec in specs:
        if spec.optional and spec.id in {
            "align-slices",
            "match-points",
            "inspect-geometry",
            "inspect-registration",
            "view-registration",
        }:
            done, detail = done_fn(spec.id, save_path, config)
            state = StageState.DONE if done else StageState.OPTIONAL
            statuses.append(StageStatus(spec, state, detail))
            continue
        done, detail = done_fn(spec.id, save_path, config)
        state = StageState.DONE if done else StageState.PENDING
        statuses.append(StageStatus(spec, state, detail))
    return statuses


def brain_stage_statuses(config: BrainPipelineConfig) -> list[StageStatus]:
    save_path = config.sample.save_path.expanduser()
    return evaluate_stage_statuses(
        brain_stage_specs(config),
        done_fn=_brain_stage_done,
        save_path=save_path,
        config=config,
    )


def spinal_stage_statuses(config: SpinalCordPipelineConfig) -> list[StageStatus]:
    save_path = config.sample.save_path.expanduser()
    return evaluate_stage_statuses(
        spinal_stage_specs(config),
        done_fn=_spinal_stage_done,
        save_path=save_path,
        config=config,
    )


def multires_stage_statuses(config: Any) -> list[StageStatus]:
    save_path = config.sample.save_path.expanduser()
    return evaluate_stage_statuses(
        multires_stage_specs(config),
        done_fn=_multires_stage_done,
        save_path=save_path,
        config=config,
    )


def slice_stages(
    specs: list[StageSpec],
    *,
    from_stage: str | None,
    through_stage: str | None,
) -> list[StageSpec]:
    ids = [spec.id for spec in specs]
    start = 0
    end = len(specs)
    if from_stage is not None:
        if from_stage not in ids:
            msg = f"Unknown stage {from_stage!r}; choose from: {', '.join(ids)}"
            raise ValueError(msg)
        start = ids.index(from_stage)
    if through_stage is not None:
        if through_stage not in ids:
            msg = f"Unknown stage {through_stage!r}; choose from: {', '.join(ids)}"
            raise ValueError(msg)
        end = ids.index(through_stage) + 1
    if start >= end:
        msg = f"Invalid stage range: from={from_stage!r} through={through_stage!r}"
        raise ValueError(msg)
    return specs[start:end]
