"""mesoSPIM overview ↔ ROI multiresolution registration."""

from lightsuite.mesospim.checkpoint import MesospimRegOptsCheckpoint, mesospim_checkpoint_path
from lightsuite.mesospim.config_models import MesospimGeometryMode, MesospimPipelineConfig
from lightsuite.mesospim.landmark_session import (
    MesospimLandmarkSession,
    default_landmark_session_path,
)

__all__ = [
    "MesospimGeometryMode",
    "MesospimLandmarkSession",
    "MesospimPipelineConfig",
    "MesospimRegOptsCheckpoint",
    "check_mesospim_geometry",
    "default_landmark_session_path",
    "mesospim_checkpoint_path",
    "run_mesospim_registration",
]


def check_mesospim_geometry(*args, **kwargs):
    from lightsuite.mesospim.runner import check_mesospim_geometry as _impl

    return _impl(*args, **kwargs)


def run_mesospim_registration(*args, **kwargs):
    from lightsuite.mesospim.runner import run_mesospim_registration as _impl

    return _impl(*args, **kwargs)
