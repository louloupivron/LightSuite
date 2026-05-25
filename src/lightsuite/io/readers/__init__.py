"""Volume reader protocol and registry."""

from lightsuite.io.readers.base import VolumeMetadata, VolumeReader
from lightsuite.io.readers.registry import open_volume, register_reader

__all__ = ["VolumeMetadata", "VolumeReader", "open_volume", "register_reader"]
