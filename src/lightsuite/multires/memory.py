"""Memory estimates and preflight warnings for multires overlap crops."""

from __future__ import annotations

import warnings
from dataclasses import dataclass
from pathlib import Path


# Fixed overview crop + resampled moving crop are both held as float32, plus
# temporary SimpleITK / NumPy working copies during streaming and Elastix prep.
_PEAK_VOLUME_FACTOR = 2.5


@dataclass(frozen=True)
class SystemMemoryInfo:
    total_gb: float
    available_gb: float


@dataclass(frozen=True)
class OverlapMemoryEstimate:
    crop_size_xyz: tuple[int, int, int]
    crop_voxels: int
    one_volume_gb: float
    estimated_peak_gb: float
    max_slab_gb: float


def system_memory_info(meminfo_path: Path | None = None) -> SystemMemoryInfo | None:
    """Read total / available RAM from ``/proc/meminfo`` (Linux)."""
    path = Path("/proc/meminfo") if meminfo_path is None else Path(meminfo_path)
    if not path.is_file():
        return None

    values: dict[str, float] = {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith(("MemTotal:", "MemAvailable:")):
            continue
        key, raw, *_rest = line.split()
        values[key.rstrip(":")] = float(raw) / (1024.0 * 1024.0)

    total = values.get("MemTotal")
    available = values.get("MemAvailable")
    if total is None or available is None:
        return None
    return SystemMemoryInfo(total_gb=float(total), available_gb=float(available))


def estimate_overlap_crop_memory(
    crop_size_xyz: list[int] | tuple[int, int, int],
    *,
    max_slab_bytes: int = 500_000_000,
) -> OverlapMemoryEstimate:
    """Approximate peak RAM for materialising one overview/ROI overlap crop pair."""
    sx, sy, sz = (max(0, int(v)) for v in crop_size_xyz)
    voxels = sx * sy * sz
    one_volume_gb = voxels * 4.0 / 1e9
    max_slab_gb = max(0, int(max_slab_bytes)) / 1e9
    estimated_peak_gb = _PEAK_VOLUME_FACTOR * one_volume_gb + max_slab_gb
    return OverlapMemoryEstimate(
        crop_size_xyz=(sx, sy, sz),
        crop_voxels=voxels,
        one_volume_gb=one_volume_gb,
        estimated_peak_gb=estimated_peak_gb,
        max_slab_gb=max_slab_gb,
    )


def warn_if_overlap_memory_exceeds_system(
    crop_size_xyz: list[int] | tuple[int, int, int],
    *,
    max_slab_bytes: int = 500_000_000,
    meminfo_path: Path | None = None,
) -> OverlapMemoryEstimate:
    """Warn when the estimated overlap crop peak exceeds available or total RAM."""
    estimate = estimate_overlap_crop_memory(crop_size_xyz, max_slab_bytes=max_slab_bytes)
    system = system_memory_info(meminfo_path=meminfo_path)
    if system is None:
        return estimate

    sx, sy, sz = estimate.crop_size_xyz
    crop_desc = (
        f"overlap crop ~{sx}×{sy}×{sz} voxels "
        f"({estimate.one_volume_gb:.1f} GB float32 each; "
        f"estimated peak ~{estimate.estimated_peak_gb:.1f} GB including buffers)"
    )
    machine_desc = (
        f"machine has {system.available_gb:.1f} GB available "
        f"/ {system.total_gb:.1f} GB total"
    )

    if estimate.estimated_peak_gb > system.available_gb:
        warnings.warn(
            (
                f"Multires registration {crop_desc} exceeds available RAM "
                f"({machine_desc}). The process is likely to be killed (OOM). "
                "Reduce the overlap with a negative multires.registration.overlap_margin_um, "
                "increase registration_bin after a smaller crop, or free memory before continuing."
            ),
            UserWarning,
            stacklevel=2,
        )
    elif estimate.estimated_peak_gb > 0.75 * system.available_gb:
        warnings.warn(
            (
                f"Multires registration {crop_desc} is close to available RAM "
                f"({machine_desc}). Consider shrinking the overlap crop if the run is unstable."
            ),
            UserWarning,
            stacklevel=2,
        )
    return estimate
