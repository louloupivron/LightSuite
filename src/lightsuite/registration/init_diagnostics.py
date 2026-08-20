"""Init-registration quality metrics and checkpoint summary."""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Literal

from rich.console import Console

InitRegistrationStatus = Literal["good", "moderate", "poor", "failed"]


@dataclass
class InitRegistrationDiagnostics:
    """Coarse init-registration checkpoint (voxel units at registration resolution)."""

    sample_shape: list[int]
    atlas_shape: list[int]
    orientation: list[int]
    registration_resolution_um: float
    cloud_threshold: float
    sample_cloud_subsample: float
    sample_cloud_points: int
    atlas_cloud_points: int
    alignment_backend: str
    similarity_scale: float
    median_error_sample_to_atlas_vox: float
    median_error_atlas_to_sample_vox: float
    median_error_vox: float
    p95_error_sample_to_atlas_vox: float
    p95_error_atlas_to_sample_vox: float
    inlier_fraction: float
    inlier_threshold_vox: float
    auto_pairs: int
    warped_boundary_voxels: int
    alignment_elapsed_s: float
    triage_elapsed_s: float
    preview_elapsed_s: float
    status: InitRegistrationStatus
    status_message: str
    warnings: list[str] = field(default_factory=list)
    sample_mask_points: int | None = None
    sample_trim_points: int | None = None
    sample_downsample_points: int | None = None
    sample_denoise_points: int | None = None

    def to_dict(self) -> dict:
        return asdict(self)

    def save(self, path: Path) -> None:
        path = path.expanduser()
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(self.to_dict(), indent=2), encoding="utf-8")

    @classmethod
    def load(cls, path: Path) -> InitRegistrationDiagnostics:
        raw = json.loads(path.expanduser().read_text(encoding="utf-8"))
        return cls(**raw)

    def print_summary(self, *, console: Console | None = None) -> None:
        """Log warnings only (status/metrics stay in the diagnostics JSON)."""
        if not self.warnings:
            return
        out = console or Console()
        for warning in self.warnings:
            out.print(f"[yellow]Warning:[/yellow] {warning}")


def classify_init_registration_status(
    *,
    median_error_vox: float,
    auto_pairs: int,
    inlier_fraction: float,
    similarity_scale: float,
    warped_boundary_voxels: int,
    sample_cloud_points: int,
    atlas_cloud_points: int,
    alignment_backend: str,
) -> tuple[InitRegistrationStatus, str, list[str]]:
    """Return (status, message, warnings) for the coarse init stage."""
    warnings: list[str] = []

    if alignment_backend == "icp":
        warnings.append(
            "BCPD not available — using Open3D ICP fallback (install bcpd for best results)."
        )
    if sample_cloud_points < 100:
        warnings.append(
            "Sparse sample cloud — try lowering registration.cloud_threshold or "
            "raising sample_cloud_subsample."
        )
    if atlas_cloud_points < 1_000:
        warnings.append("Few atlas features extracted — verify atlas_dir contents.")
    if not 0.75 <= similarity_scale <= 1.35:
        warnings.append(
            f"Similarity scale {similarity_scale:.3f} outside expected 0.75–1.35 range."
        )
    if warped_boundary_voxels == 0:
        warnings.append(
            "No atlas boundaries overlapped the sample after warping — check orientation "
            "and qc/previews/dim{1,2,3}_initial_registration.png."
        )
    if inlier_fraction < 0.35:
        warnings.append(
            f"Only {inlier_fraction * 100:.0f}% of sample points within "
            f"{25:g} vox of atlas features."
        )

    if auto_pairs < 4:
        return (
            "failed",
            "Too few auto control point pairs for downstream registration.",
            warnings,
        )
    if median_error_vox > 25 or auto_pairs < 20:
        return (
            "poor",
            "Coarse alignment likely wrong — verify orientation and preview PNGs.",
            warnings,
        )
    if median_error_vox > 15 or auto_pairs < 50 or inlier_fraction < 0.5:
        return (
            "moderate",
            "Usable starting point; align-slices or manual match-points recommended.",
            warnings,
        )
    return (
        "good",
        "Coarse alignment ready for match-points or register.",
        warnings,
    )
