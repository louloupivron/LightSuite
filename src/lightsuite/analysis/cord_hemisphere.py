"""Left/right hemisegment split for the Fiederling spinal cord atlas."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import tifffile

from lightsuite.analysis.hemisphere import SIDE_LABELS

HEMISPHERE_ANNOTATION_FILENAME = "Hemisphere_Annotation.tif"
REGISTERED_HEMISPHERE_FILENAME = "hemisphere_registered.tiff"

# Fiederling ``Hemisphere_Annotation.tif`` uses 255 for one hemisegment and 0 for the other.
HEMI_ACTIVE_VALUE = 255


def resolve_hemisphere_annotation_path(atlas_dir: Path) -> Path:
    path = Path(atlas_dir).expanduser().resolve() / HEMISPHERE_ANNOTATION_FILENAME
    if not path.is_file():
        msg = (
            f"Missing {path}. Download the Fiederling atlas package "
            f"({HEMISPHERE_ANNOTATION_FILENAME}) to enable split_hemispheres."
        )
        raise FileNotFoundError(msg)
    return path


def load_fiederling_hemisphere_native(atlas_dir: Path) -> np.ndarray:
    """Load the native (Z, Y, X) hemisphere mask from the atlas directory."""
    return tifffile.imread(resolve_hemisphere_annotation_path(atlas_dir))


def cord_hemisphere_side_volume(
    hemisphere_mask: np.ndarray,
    annotation: np.ndarray,
    *,
    flip: bool = False,
) -> np.ndarray:
    """Map each in-cord voxel to side id ``0`` (right) or ``1`` (left).

    The Fiederling mask is binary (``0`` / ``255``). By default ``255`` maps to
    ``right`` and ``0`` to ``left``; set ``flip=True`` to swap. Voxels outside
    the annotation or with an unassigned mask value are ``-1``.
    """
    hem = np.asarray(hemisphere_mask)
    ann = np.asarray(annotation)
    if hem.shape != ann.shape:
        msg = f"Hemisphere mask shape {hem.shape} != annotation {ann.shape}"
        raise ValueError(msg)

    right_mask = hem == HEMI_ACTIVE_VALUE
    left_mask = hem == 0
    tissue = ann > 0
    side = np.full(ann.shape, -1, dtype=np.int8)
    if flip:
        side[right_mask & tissue] = 1
        side[left_mask & tissue] = 0
    else:
        side[right_mask & tissue] = 0
        side[left_mask & tissue] = 1
    return side


def hemisphere_label_from_side(side_id: int) -> str | None:
    if side_id < 0 or side_id >= len(SIDE_LABELS):
        return None
    return SIDE_LABELS[side_id]


__all__ = [
    "HEMI_ACTIVE_VALUE",
    "HEMISPHERE_ANNOTATION_FILENAME",
    "REGISTERED_HEMISPHERE_FILENAME",
    "cord_hemisphere_side_volume",
    "hemisphere_label_from_side",
    "load_fiederling_hemisphere_native",
    "resolve_hemisphere_annotation_path",
]
