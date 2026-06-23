"""Pydantic models mirroring MATLAB opts / regopts structs."""

from __future__ import annotations

from enum import Enum
from pathlib import Path
from typing import Annotated

from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator


class SourceFormat(str, Enum):
    AUTO = "auto"
    TIFF_STACK = "tiff_stack"
    CZI = "czi"
    OME_ZARR = "ome_zarr"
    IMARIS = "imaris"


class TiffLayout(str, Enum):
    PLANE_PER_FILE = "planeperfile"
    CHANNEL_PER_FILE = "channelperfile"


class BrainAtlasId(str, Enum):
    ALLEN = "allen"
    PERENS = "perens"


class AtlasSource(str, Enum):
    FILES = "files"
    BRAINGLOBE = "brainglobe"


class SampleSourceConfig(BaseModel):
    format: SourceFormat = SourceFormat.AUTO
    path: Path | None = None
    tiff_type: TiffLayout = TiffLayout.CHANNEL_PER_FILE
    channels: Annotated[list[Path], Field(min_length=1)] | None = Field(
        default=None,
        description=(
            "Optional list of planeperfile roots (one folder per channel, Terastitcher-style). "
            "Channel index follows list order. Requires tiff_type: planeperfile."
        ),
    )

    @field_validator("path")
    @classmethod
    def path_must_exist(cls, value: Path | None) -> Path | None:
        if value is None:
            return None
        if not value.expanduser().exists():
            msg = f"Sample source path does not exist: {value}"
            raise ValueError(msg)
        return value.expanduser().resolve()

    @field_validator("channels")
    @classmethod
    def channels_must_exist(cls, value: list[Path] | None) -> list[Path] | None:
        if value is None:
            return None
        resolved: list[Path] = []
        for item in value:
            expanded = item.expanduser().resolve()
            if not expanded.is_dir():
                msg = f"Channel folder not found: {expanded}"
                raise ValueError(msg)
            resolved.append(expanded)
        return resolved

    @model_validator(mode="after")
    def validate_source_paths(self) -> SampleSourceConfig:
        if self.channels is not None:
            if self.tiff_type != TiffLayout.PLANE_PER_FILE:
                msg = "source.channels is only supported with tiff_type: planeperfile"
                raise ValueError(msg)
            if self.path is None:
                self.path = self.channels[0]
            return self
        if self.path is None:
            msg = "sample.source.path is required when source.channels is not set"
            raise ValueError(msg)
        return self

    @property
    def channel_roots(self) -> tuple[Path, ...] | None:
        if self.channels is None:
            return None
        return tuple(self.channels)


class SampleConfig(BaseModel):
    name: str = Field(min_length=1)
    source: SampleSourceConfig
    scratch: Path
    save_path: Path
    voxel_um: Annotated[list[float], Field(min_length=3, max_length=3)] | None = None

    @field_validator("scratch", "save_path")
    @classmethod
    def expand_paths(cls, value: Path) -> Path:
        return value.expanduser()


class AtlasConfig(BaseModel):
    provider: BrainAtlasId = BrainAtlasId.ALLEN
    resolution_um: float = Field(default=10.0, gt=0)
    atlas_dir: Path | None = None
    source: AtlasSource = Field(
        default=AtlasSource.FILES,
        description="Atlas provider: local NIfTI files (files) or BrainGlobe Atlas API (brainglobe).",
    )
    brainglobe_name: str | None = Field(
        default=None,
        description=(
            "BrainGlobe registry name when source=brainglobe (e.g. perens_lsfm_mouse_20um). "
            "Inferred from provider and resolution_um when omitted."
        ),
    )

    @field_validator("atlas_dir")
    @classmethod
    def expand_atlas_dir(cls, value: Path | None) -> Path | None:
        if value is None:
            return None
        return value.expanduser()


class RegistrationConfig(BaseModel):
    resolution_um: float = Field(default=20.0, gt=0)
    channel_primary: int = Field(default=1, ge=1)
    channel_secondary: int | None = Field(default=None, ge=1)
    bspline_spatial_scale_mm: float = Field(default=0.64, gt=0)
    control_point_weight: float = Field(default=0.2, ge=0, le=1)
    augment_points: bool = False
    dual_channel_mi_weight_autofluor: float = Field(default=0.4, ge=0, le=1)
    dual_channel_mi_weight_signal: float = Field(default=0.4, ge=0, le=1)
    orientation: Annotated[list[int], Field(min_length=3, max_length=3)] | None = Field(
        default=None,
        description="Axis permutation e.g. [1, 2, 3]. Loaded from brain_orientation.txt if unset.",
    )
    cloud_threshold: float = Field(default=5.0, gt=0)
    sample_cloud_subsample: float = Field(
        default=0.1,
        gt=0,
        le=1.0,
        description="Random fraction of gradient sample points kept (MATLAB pcdownsample=0.1).",
    )
    outlier_ratio: float = Field(default=0.01, ge=0, le=1)
    bcpd_path: str | None = Field(
        default=None,
        description="Optional path to bcpd / bcpd.exe. Searched on PATH when unset.",
    )
    ap_pair_tolerance_vox: float = Field(
        default=12.0,
        gt=0,
        description="AP residual tolerance when refining auto control points (registration voxels).",
    )
    ap_pair_min_kept: int = Field(
        default=24,
        ge=4,
        description="Minimum auto pairs to keep after AP filtering (relaxes tolerance if needed).",
    )
    use_slice_correspondence_affine: bool = Field(
        default=True,
        description=(
            "When slice_correspondence.json has confirmed anchors, compose a "
            "correspondence-informed affine correction before B-spline registration."
        ),
    )
    use_slice_correspondence_landmarks: bool = Field(
        default=True,
        description=(
            "Add align-slices anchor pairs as extra B-spline landmarks in register."
        ),
    )
    correspondence_landmark_weight: float = Field(
        default=0.2,
        ge=0,
        le=1,
        description=(
            "Minimum landmark metric weight when correspondence B-spline landmarks are added."
        ),
    )
    correspondence_landmark_max_count: int = Field(
        default=96,
        ge=4,
        description="Maximum total B-spline landmark pairs after merging correspondence anchors.",
    )


class DetectionBackend(str, Enum):
    CLASSICAL = "classical"
    CELLPOSE = "cellpose"
    STARDIST = "stardist"


class DetectionConfig(BaseModel):
    enabled: bool = True
    backend: DetectionBackend = DetectionBackend.CLASSICAL
    cell_diameter_um: float = Field(default=14.0, gt=0)
    thresholds: Annotated[list[float], Field(min_length=2, max_length=2)] = [0.5, 0.4]
    channel: int | None = None
    debug: bool = False
    save_cell_images: bool = False
    write_to_csv: bool = False

    @field_validator("thresholds")
    @classmethod
    def thresholds_ordered(cls, value: list[float]) -> list[float]:
        if value[0] < value[1]:
            msg = "First detection threshold should be >= second threshold."
            raise ValueError(msg)
        return value


class ComputeConfig(BaseModel):
    use_gpu: bool = True
    workers: int = Field(default=4, ge=1)
    max_in_memory_scratch_gb: float = Field(
        default=24.0,
        gt=0,
        description=(
            "Keep the XY-downsampled scratch volume in RAM when it fits below this size; "
            "larger stacks use a disk memmap on sample.scratch."
        ),
    )


class ExportConfig(BaseModel):
    registered_volume_format: str = "ome_zarr"
    write_pyramid: bool = True
    write_cells_csv: bool = True
    save_registered_volume: bool = False


class AnalysisConfig(BaseModel):
    """Post-registration analysis (region stats, cell counts, taxonomy)."""

    write_tidy_csv: bool = Field(
        default=True,
        description="Emit long-form chanXX_region_stats.csv (with region names) during export.",
    )
    count_points: bool = Field(
        default=True,
        description="Bin imported atlas-space point clouds into per-region cell counts/densities.",
    )
    point_labels: list[str] | None = Field(
        default=None,
        description="Import labels to count (matches *_atlas_coords.npz stems); None = all found.",
    )


class AnnotationFormat(str, Enum):
    """LightSuite Sample Space v1 — native-resolution exports only."""

    POINTS_CSV = "points_csv"
    MASK_TIFF = "mask_tiff"


class AnnotationImportConfig(BaseModel):
    """Native sample-space annotation to register after brain register."""

    format: AnnotationFormat
    path: Path
    label: str = Field(
        default="",
        description="Output filename stem; defaults to the input file stem.",
    )

    @field_validator("path")
    @classmethod
    def expand_import_path(cls, value: Path) -> Path:
        return value.expanduser()


class ImportConfig(BaseModel):
    annotations: list[AnnotationImportConfig] = Field(default_factory=list)
    write_csv: bool = True


class BrainPipelineConfig(BaseModel):
    """Top-level brain lightsheet pipeline configuration."""

    model_config = ConfigDict(populate_by_name=True)

    sample: SampleConfig
    atlas: AtlasConfig = Field(default_factory=AtlasConfig)
    registration: RegistrationConfig = Field(default_factory=RegistrationConfig)
    detection: DetectionConfig = Field(default_factory=DetectionConfig)
    compute: ComputeConfig = Field(default_factory=ComputeConfig)
    export: ExportConfig = Field(default_factory=ExportConfig)
    analysis: AnalysisConfig = Field(default_factory=AnalysisConfig)
    import_config: ImportConfig | None = Field(default=None, alias="import")

    @model_validator(mode="after")
    def perens_atlas_resolution(self) -> BrainPipelineConfig:
        if self.atlas.provider == BrainAtlasId.PERENS and self.atlas.resolution_um != 20.0:
            # Perens LSFM atlas is 20 µm isotropic; warn via validation note in loader if needed.
            pass
        return self

    @property
    def data_folder(self) -> Path:
        """MATLAB-compatible alias for sample source path."""
        return self.sample.source.path

    @property
    def fproc(self) -> Path:
        """MATLAB-compatible alias for scratch directory."""
        return self.sample.scratch
