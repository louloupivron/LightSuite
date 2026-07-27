"""Microscope-agnostic overview ↔ ROI multiresolution registration."""

from lightsuite.multires.checkpoint import MultiresRegOptsCheckpoint, multires_checkpoint_path
from lightsuite.multires.config_models import MultiresGeometryMode, MultiresPipelineConfig
from lightsuite.multires.landmark_session import (
    LandmarkFitMode,
    MultiresLandmarkSession,
    default_landmark_session_path,
)
from lightsuite.multires.manifest import (
    MANIFEST_FORMAT,
    MultiresPairManifest,
    load_pair_manifest,
    save_pair_manifest,
)

__all__ = [
    "MANIFEST_FORMAT",
    "LandmarkFitMode",
    "MultiresGeometryMode",
    "MultiresLandmarkSession",
    "MultiresPairManifest",
    "MultiresPipelineConfig",
    "MultiresRegOptsCheckpoint",
    "check_multires_geometry",
    "default_landmark_session_path",
    "load_pair_manifest",
    "multires_checkpoint_path",
    "run_multires_registration",
    "save_pair_manifest",
]


def check_multires_geometry(*args, **kwargs):
    from lightsuite.multires.runner import check_multires_geometry as _impl

    return _impl(*args, **kwargs)


def run_multires_registration(*args, **kwargs):
    from lightsuite.multires.runner import run_multires_registration as _impl

    return _impl(*args, **kwargs)
