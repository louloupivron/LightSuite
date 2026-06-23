"""Single source of truth for the left/right voxel split in atlas space.

Both region intensity statistics and cell counts must split the atlas volume into
the same two hemispheres so their ``right`` / ``left`` columns line up. Side index
``0`` maps to ``RightSideIntensity`` / right-hemisphere counts and side index ``1``
to the left, matching :func:`lightsuite.export.parcellation.write_parcellation_csv`.

Conventions:

- **BrainGlobe** (when ``brainglobe_name`` is set): use the packaged hemisphere
  volume (labels 1=left, 2=right).
- **Allen (files)**: split along the last array axis (Z) at its midpoint.
- **Perens / generic (files)**: split along ``ml_axis`` at the mean coordinate of
  labelled brain voxels.
"""

from __future__ import annotations

import numpy as np

#: Side index → hemisphere label. Index 0 is the "Right" column downstream.
SIDE_LABELS: tuple[str, str] = ("right", "left")


def hemisphere_side_masks(
    annotation: np.ndarray,
    atlas_id: str,
    *,
    ml_axis: int = 2,
    brainglobe_name: str | None = None,
) -> list[np.ndarray]:
    """Return ``[side0_mask, side1_mask]`` boolean volumes for an annotation array."""
    if brainglobe_name is not None:
        from lightsuite.atlas.brainglobe_backend import (
            hemisphere_side_masks_from_brainglobe,
            load_brainglobe_hemispheres,
        )

        hem = load_brainglobe_hemispheres(brainglobe_name)
        return hemisphere_side_masks_from_brainglobe(hem, annotation)

    av = np.asanyarray(annotation)
    atlas_id = atlas_id.lower().strip()

    if atlas_id == "allen":
        half = av.shape[2] // 2
        side0 = np.zeros(av.shape, dtype=bool)
        side1 = np.zeros(av.shape, dtype=bool)
        side0[:, :, :half] = True
        side1[:, :, half:] = True
        return [side0, side1]

    axis = ml_axis - 1
    if not 0 <= axis < av.ndim:
        msg = f"ml_axis {ml_axis} out of range for annotation with {av.ndim} dims"
        raise ValueError(msg)
    brain = av > 0
    if not brain.any():
        empty = np.zeros(av.shape, dtype=bool)
        return [empty, empty.copy()]
    coords = np.indices(av.shape)[axis].astype(np.float32)
    split_plane = round(float(coords[brain].mean()))
    side0 = (coords <= split_plane) & brain
    side1 = (coords > split_plane) & brain
    return [side0, side1]


def hemisphere_side_volume(
    annotation: np.ndarray,
    atlas_id: str,
    *,
    ml_axis: int = 2,
    brainglobe_name: str | None = None,
) -> np.ndarray:
    """Per-voxel side id volume: ``0`` (right), ``1`` (left), ``-1`` (neither)."""
    masks = hemisphere_side_masks(
        annotation,
        atlas_id,
        ml_axis=ml_axis,
        brainglobe_name=brainglobe_name,
    )
    side = np.full(np.asanyarray(annotation).shape, -1, dtype=np.int8)
    for index, mask in enumerate(masks):
        side[mask] = index
    return side
