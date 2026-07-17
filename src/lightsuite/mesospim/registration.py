"""Elastix registration for mesoSPIM overview / ROI pairs."""

from __future__ import annotations

from pathlib import Path

from lightsuite.mesospim.config_models import MesospimPipelineConfig
from lightsuite.mesospim.prepare import MesospimPreparedPair, prepare_mesospim_registration_pair
from lightsuite.multires.registration import (
    MultiresRegistrationResult,
    build_elastix_parameter_object,
    register_roi_to_overview as _register_roi_to_overview,
    sanitize_experiment_name,
)

MesospimRegistrationResult = MultiresRegistrationResult


def register_roi_to_overview(
    *,
    cfg: MesospimPipelineConfig,
    prepared: MesospimPreparedPair | None = None,
    output_dir: Path,
    experiment_slug: str,
    overview_stem: str,
    roi_stem: str,
    registration_bin: int,
    elastix_stages: list[str],
    write_full_overview_canvas: bool = True,
) -> MesospimRegistrationResult:
    """Run itk-elastix on a prepared mesoSPIM overview / ROI pair."""
    if prepared is None:
        prepared = prepare_mesospim_registration_pair(cfg)
    return _register_roi_to_overview(
        prepared=prepared,
        output_dir=output_dir,
        experiment_slug=experiment_slug,
        overview_stem=overview_stem,
        roi_stem=roi_stem,
        registration_bin=registration_bin,
        elastix_stages=elastix_stages,
        write_full_overview_canvas=write_full_overview_canvas,
    )


__all__ = [
    "MesospimRegistrationResult",
    "build_elastix_parameter_object",
    "register_roi_to_overview",
    "sanitize_experiment_name",
]
