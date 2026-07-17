"""Landmark session persistence for mesoSPIM overview / ROI placement."""

from __future__ import annotations

from pathlib import Path

from lightsuite.multires.landmark_session import LandmarkFitMode, MultiresLandmarkSession

MesospimLandmarkSession = MultiresLandmarkSession


def default_landmark_session_path(save_path: Path) -> Path:
    return save_path.expanduser() / "mesospim_landmarks.json"


__all__ = [
    "LandmarkFitMode",
    "MesospimLandmarkSession",
    "default_landmark_session_path",
]
