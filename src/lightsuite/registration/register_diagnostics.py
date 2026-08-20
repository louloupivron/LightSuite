"""Register-step quality metrics and checkpoint summary."""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Literal

from rich.console import Console

RegistrationStatus = Literal["good", "moderate", "poor", "failed"]

_STATUS_STYLE: dict[RegistrationStatus, str] = {
    "good": "green",
    "moderate": "yellow",
    "poor": "red",
    "failed": "bold red",
}

_STATUS_LABEL: dict[RegistrationStatus, str] = {
    "good": "GOOD",
    "moderate": "MODERATE",
    "poor": "POOR",
    "failed": "FAILED",
}


@dataclass
class RegistrationDiagnostics:
    """Affine + B-spline registration checkpoint (voxel units at registration resolution)."""

    sample_shape: list[int]
    atlas_shape: list[int]
    orientation: list[int]
    registration_resolution_um: float
    n_manual_pairs: int
    n_auto_pairs: int
    n_landmark_pairs: int
    control_point_weight: float
    use_multistep: bool
    use_dual_channel_mi: bool
    bspline_spatial_scale_mm: float
    bspline_bending_weight: float = 0.0
    dual_channel_mi_weight_primary: float | None = None
    dual_channel_mi_weight_secondary: float | None = None
    affine_median_error_vox: float = 0.0
    affine_p95_error_vox: float = 0.0
    affine_max_error_vox: float = 0.0
    affine_median_manual_vox: float | None = None
    affine_median_auto_vox: float | None = None
    affine_median_coarse_auto_vox: float | None = None
    bspline_landmark_metric_mm: float | None = None
    bspline_landmark_metric_vox: float | None = None
    annotation_label_voxels: int = 0
    load_elapsed_s: float = 0.0
    affine_warp_elapsed_s: float = 0.0
    bspline_elapsed_s: float = 0.0
    transformix_elapsed_s: float = 0.0
    status: RegistrationStatus = "moderate"
    status_message: str = ""
    warnings: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        return asdict(self)

    def save(self, path: Path) -> None:
        path = path.expanduser()
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(self.to_dict(), indent=2), encoding="utf-8")

    @classmethod
    def load(cls, path: Path) -> RegistrationDiagnostics:
        raw = json.loads(path.expanduser().read_text(encoding="utf-8"))
        return cls(**raw)

    def print_summary(self, *, console: Console | None = None) -> None:
        """Log GOOD/MODERATE/POOR/FAILED plus the status comment (metrics stay in JSON)."""
        out = console or Console()
        style = _STATUS_STYLE[self.status]
        label = _STATUS_LABEL[self.status]
        out.print(f"[{style}]{label}[/{style}] — {self.status_message}")
        for warning in self.warnings:
            out.print(f"[yellow]Warning:[/yellow] {warning}")


def landmark_mm_to_vox(metric_mm: float, registration_resolution_um: float) -> float:
    spacing_mm = registration_resolution_um * 1e-3
    return metric_mm / spacing_mm if spacing_mm > 0 else metric_mm


def classify_registration_status(
    *,
    affine_median_error_vox: float,
    affine_p95_error_vox: float,
    n_manual: int,
    n_landmark_pairs: int,
    bspline_landmark_metric_vox: float | None,
    annotation_label_voxels: int,
    use_multistep: bool,
) -> tuple[RegistrationStatus, str, list[str]]:
    """Return (status, message, warnings) for the register stage."""
    warnings: list[str] = []

    if n_manual == 0:
        warnings.append(
            "Auto-only landmarks — consider match-points for difficult samples."
        )
    if not use_multistep:
        warnings.append("Single-step B-spline schedule — lower quality than multistep.")
    if bspline_landmark_metric_vox is None:
        warnings.append("B-spline landmark metric unavailable from elastix IterationInfo.")
    if affine_p95_error_vox > 25 and affine_median_error_vox <= 15:
        warnings.append(
            "Affine p95 is high while median is low — localized landmark mismatch."
        )
    if annotation_label_voxels == 0:
        warnings.append(
            "No annotated brain regions overlapped the sample after B-spline."
        )

    if n_landmark_pairs < 4:
        return (
            "failed",
            "Too few landmark pairs for registration.",
            warnings,
        )

    landmark_bad = (
        bspline_landmark_metric_vox is not None and bspline_landmark_metric_vox > 25.0
    )
    landmark_moderate = (
        bspline_landmark_metric_vox is not None and bspline_landmark_metric_vox > 12.0
    )

    if affine_median_error_vox > 25 or landmark_bad or annotation_label_voxels == 0:
        return (
            "poor",
            "Registration quality is low — review previews and add manual landmarks.",
            warnings,
        )
    if (
        affine_median_error_vox > 15
        or landmark_moderate
        or (n_manual == 0 and affine_median_error_vox > 8)
    ):
        return (
            "moderate",
            "Registration usable; manual landmarks or align-slices may improve export.",
            warnings,
        )
    return (
        "good",
        "Registration ready for export.",
        warnings,
    )
