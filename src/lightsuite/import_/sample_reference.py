"""Native sample-space reference frame (LightSuite Sample Space v1)."""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

SAMPLE_SPACE_FORMAT = "lightsuite_sample_space_v1"


@dataclass
class SampleReference:
    """Reference grid for external segmentation exports (native resolution only)."""

    format: str
    sample_name: str
    shape_yxz: list[int]
    voxel_um: list[float]
    index_base: int
    axis_order: str
    coordinate_units: str
    orientation_applied: bool

    @classmethod
    def from_checkpoint(
        cls,
        *,
        sample_name: str,
        ny: int,
        nx: int,
        nz: int,
        voxel_um: list[float],
    ) -> SampleReference:
        return cls(
            format=SAMPLE_SPACE_FORMAT,
            sample_name=sample_name,
            shape_yxz=[int(ny), int(nx), int(nz)],
            voxel_um=[float(v) for v in voxel_um],
            index_base=1,
            axis_order="xyz",
            coordinate_units="voxel_indices",
            orientation_applied=False,
        )

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    def save(self, path: Path) -> Path:
        path = path.expanduser()
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(self.to_dict(), indent=2), encoding="utf-8")
        return path

    @classmethod
    def load(cls, path: Path) -> SampleReference:
        raw = json.loads(path.expanduser().read_text(encoding="utf-8"))
        return cls(**raw)

    @property
    def shape_tuple(self) -> tuple[int, int, int]:
        ny, nx, nz = self.shape_yxz
        return int(ny), int(nx), int(nz)


def sample_reference_path(save_path: Path) -> Path:
    return save_path.expanduser() / "sample_reference.json"


def write_sample_reference(
    save_path: Path,
    *,
    sample_name: str,
    ny: int,
    nx: int,
    nz: int,
    voxel_um: list[float],
) -> Path:
    ref = SampleReference.from_checkpoint(
        sample_name=sample_name,
        ny=ny,
        nx=nx,
        nz=nz,
        voxel_um=voxel_um,
    )
    return ref.save(sample_reference_path(save_path))


def load_sample_reference(save_path: Path) -> SampleReference:
    path = sample_reference_path(save_path)
    if not path.is_file():
        msg = (
            f"Missing {path}. Run 'lightsuite brain preprocess' first to publish "
            "the native sample-space reference."
        )
        raise FileNotFoundError(msg)
    return SampleReference.load(path)


def validate_mask_against_reference(
    volume_shape_yxz: tuple[int, int, int],
    voxel_um: list[float],
    reference: SampleReference,
) -> None:
    expected = reference.shape_tuple
    if tuple(int(v) for v in volume_shape_yxz) != expected:
        msg = (
            f"Mask shape (Y, X, Z)={volume_shape_yxz} does not match native sample "
            f"reference {expected} from sample_reference.json. "
            "Resample the mask externally to native resolution before import."
        )
        raise ValueError(msg)
    ref_um = [float(v) for v in reference.voxel_um]
    if not all(abs(a - b) < 1e-6 for a, b in zip(voxel_um, ref_um, strict=True)):
        msg = (
            f"Mask voxel_um={voxel_um} does not match sample reference {ref_um}. "
            "Segment at native sample voxel size only."
        )
        raise ValueError(msg)
