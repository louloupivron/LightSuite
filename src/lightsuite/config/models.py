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
    PRINCETON = "princeton"
    RAT = "rat"


class AtlasSource(str, Enum):
    FILES = "files"
    BRAINGLOBE = "brainglobe"


class ContentTrimMode(str, Enum):
    OFF = "off"
    AUTO = "auto"
    MANUAL = "manual"


class SampleContentCropMode(str, Enum):
    OFF = "off"
    AUTO = "auto"
    MANUAL = "manual"


class RegistrationCanvasMode(str, Enum):
    OFF = "off"
    PAD = "pad"
    CROP = "crop"
    UNION = "union"


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
    content_trim: ContentTrimMode = Field(
        default=ContentTrimMode.OFF,
        description="Trim atlas black padding for registration (export stays native atlas).",
    )
    content_margin_vox: int = Field(
        default=8,
        ge=0,
        description="Margin added around auto-detected atlas foreground bbox.",
    )
    content_box: Annotated[list[int], Field(min_length=6, max_length=6)] | None = Field(
        default=None,
        description="Manual atlas crop [y0, y1, x0, x1, z0, z1] inclusive when content_trim=manual.",
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
    sample_content_crop: SampleContentCropMode = Field(
        default=SampleContentCropMode.OFF,
        description="Crop registration TIFFs to sample foreground after preprocess.",
    )
    sample_content_margin_vox: int = Field(
        default=8,
        ge=0,
        description="Margin around auto-detected sample foreground crop.",
    )
    sample_content_box: Annotated[list[int], Field(min_length=6, max_length=6)] | None = Field(
        default=None,
        description="Manual sample crop [y0, y1, x0, x1, z0, z1] when sample_content_crop=manual.",
    )
    sample_content_trim_z: bool = Field(
        default=True,
        description="Remove sparse high-Z slices during auto sample crop (cord-style trim).",
    )
    canvas_mode: RegistrationCanvasMode = Field(
        default=RegistrationCanvasMode.OFF,
        description=(
            "Reconcile sample/atlas working grids before elastix: off=MATLAB parity "
            "(atlas warped to sample shape), pad/crop/union adjust the working canvas."
        ),
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
    spaces: list[str] = Field(
        default_factory=lambda: ["atlas"],
        description="Output coordinate spaces: atlas (sample→atlas) and/or sample (atlas→sample).",
    )
    save_sample_space_volume: bool = Field(
        default=True,
        description="Write warped atlas labels/template under volume_registered/sample_space/.",
    )


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
    parcellate_intensities: bool = Field(
        default=True,
        description="Compute per-region median intensity from registered channel volumes (spinal cord).",
    )
    intensity_channels: list[int] | None = Field(
        default=None,
        description="Registered channel indices to parcellate; None = all exported channels.",
    )
    relative_intensity_to: str = Field(
        default="none",
        description='Intensity normalization: "none" or "background" (relative to annotation id 0 per segment).',
    )
    point_labels: list[str] | None = Field(
        default=None,
        description="Import labels to count (matches *_atlas_coords.npz stems); None = all found.",
    )
    stats_spaces: list[str] = Field(
        default_factory=lambda: ["atlas"],
        description="Coordinate spaces for region_stats tables: atlas and/or sample.",
    )
    rollups: list[str] = Field(
        default_factory=list,
        description='Spinal cord rollups to append: "division" (GM/WM), "structure" (laminas/funiculi), and/or "horn" (DH/VH/C dorsal–ventral split).',
    )
    split_hemispheres: bool = Field(
        default=False,
        description="Split spinal cord stats into left/right using Hemisphere_Annotation.tif.",
    )
    hemisphere_flip: bool = Field(
        default=False,
        description="Swap left/right assignment for the Fiederling hemisphere mask (0/255).",
    )
    hemisphere_keep_whole: bool = Field(
        default=False,
        description="When split_hemispheres is true, also emit whole-cord summary rows.",
    )
    plots_dir: Path | None = Field(
        default=None,
        description="Directory for matplotlib analysis outputs; default is sample.save_path/plots.",
    )

    @field_validator("plots_dir")
    @classmethod
    def expand_plots_dir(cls, value: Path | None) -> Path | None:
        if value is None:
            return None
        return value.expanduser()


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


class CordTiffLayout(str, Enum):
    """Spinal cord TIFF layout (readSpinalCordSample.m)."""

    AUTO = "auto"
    PLANE_PER_FILE = "planeperfile"
    CHANNEL_PER_FILE = "channelperfile"
    MULTICHANNEL_SINGLE = "multichannel_single"


class CordAtlasConfig(BaseModel):
    """Fiederling et al. 2021 spinal cord atlas (external TIFF + CSV)."""

    atlas_dir: Path = Field(
        description="Directory containing Template.tif, Annotation.tif, Segments.csv, Atlas_Regions.csv.",
    )

    @field_validator("atlas_dir")
    @classmethod
    def expand_cord_atlas_dir(cls, value: Path) -> Path:
        return value.expanduser()


class CordRegistrationConfig(BaseModel):
    """Spinal cord registration settings."""

    resolution_um: float = Field(default=20.0, gt=0)
    channel_primary: int = Field(default=1, ge=1, description="Registration channel (MATLAB regchan).")
    control_point_weight: float = Field(default=0.2, ge=0, le=1)
    straightening_lambda_pos: float = Field(default=5000.0, gt=0)
    straightening_lambda_ang: float = Field(default=5000.0, gt=0)
    target_orientation_deg: float = Field(default=90.0)
    longitudinal_direction: str | None = Field(
        default=None,
        description=(
            "Optional override for sample +Z anatomy direction: "
            "'rostrocaudal' or 'caudorostral'. When unset, preprocess reads "
            "cord_orientation.txt from check-orientation."
        ),
    )


class CordSampleSourceConfig(BaseModel):
    format: SourceFormat = SourceFormat.TIFF_STACK
    path: Path | None = None
    tiff_type: CordTiffLayout = CordTiffLayout.AUTO
    channels: Annotated[list[Path], Field(min_length=1)] | None = Field(
        default=None,
        description=(
            "Optional list of planeperfile roots (one folder per channel, Terastitcher-style). "
            "Channel index follows list order. Requires tiff_type: planeperfile or auto."
        ),
    )
    skip_corrupt_slices: bool = Field(
        default=False,
        description="Skip unreadable plane-per-file slice TIFFs instead of failing (use sparingly).",
    )

    @field_validator("path")
    @classmethod
    def cord_path_must_exist(cls, value: Path | None) -> Path | None:
        if value is None:
            return None
        if not value.expanduser().exists():
            msg = f"Sample source path does not exist: {value}"
            raise ValueError(msg)
        return value.expanduser().resolve()

    @field_validator("channels")
    @classmethod
    def cord_channels_must_exist(cls, value: list[Path] | None) -> list[Path] | None:
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
    def validate_cord_source_paths(self) -> CordSampleSourceConfig:
        if self.channels is not None:
            if self.tiff_type not in (CordTiffLayout.PLANE_PER_FILE, CordTiffLayout.AUTO):
                msg = (
                    "source.channels is only supported with tiff_type: planeperfile "
                    "(or auto, which resolves to planeperfile)"
                )
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


class CordSampleConfig(BaseModel):
    name: str = Field(min_length=1)
    source: CordSampleSourceConfig
    scratch: Path
    save_path: Path
    voxel_um: Annotated[list[float], Field(min_length=3, max_length=3)]

    @field_validator("scratch", "save_path")
    @classmethod
    def expand_cord_paths(cls, value: Path) -> Path:
        return value.expanduser()


class SpinalCordPipelineConfig(BaseModel):
    """Top-level spinal cord lightsheet pipeline configuration."""

    model_config = ConfigDict(populate_by_name=True)

    sample: CordSampleConfig
    atlas: CordAtlasConfig
    registration: CordRegistrationConfig = Field(default_factory=CordRegistrationConfig)
    compute: ComputeConfig = Field(default_factory=ComputeConfig)
    export: ExportConfig = Field(default_factory=ExportConfig)
    analysis: AnalysisConfig = Field(default_factory=AnalysisConfig)
    import_config: ImportConfig | None = Field(default=None, alias="import")

    @property
    def data_folder(self) -> Path:
        path = self.sample.source.path
        if path is None:
            msg = "sample.source.path is unset"
            raise RuntimeError(msg)
        return path

    @property
    def lsfolder(self) -> Path:
        return self.sample.save_path.expanduser()


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
