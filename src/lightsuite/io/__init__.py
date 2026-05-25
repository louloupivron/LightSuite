"""Volume I/O and reader registry."""

from lightsuite.io.discover import TiffStackDiscovery, discover_tiff_stack
from lightsuite.io.readers.registry import open_volume
from lightsuite.io.readers.tiff_stack import TiffStackReader

__all__ = [
    "TiffStackDiscovery",
    "TiffStackReader",
    "discover_tiff_stack",
    "ensure_output_dirs",
    "open_sample_reader",
    "open_volume",
]


def open_sample_reader(config) -> TiffStackReader:
    """Open a volume reader for a BrainPipelineConfig sample."""
    from lightsuite.config.models import BrainPipelineConfig

    if not isinstance(config, BrainPipelineConfig):
        msg = "config must be BrainPipelineConfig"
        raise TypeError(msg)
    voxel = tuple(config.sample.voxel_um) if config.sample.voxel_um else None
    return TiffStackReader.open(
        config.sample.source.path,
        tiff_type=config.sample.source.tiff_type,
        voxel_um=voxel,
    )


def ensure_output_dirs(config) -> None:
    """Create scratch and save directories if missing."""
    config.sample.scratch.mkdir(parents=True, exist_ok=True)
    config.sample.save_path.mkdir(parents=True, exist_ok=True)
