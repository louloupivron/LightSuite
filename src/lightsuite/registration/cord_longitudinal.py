"""Longitudinal sample ↔ atlas correspondence for spinal cord z-init."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from lightsuite.gui.slice_correspondence import SliceAnchor, SliceCorrespondence
from lightsuite.registration.cord_affine import build_cord_z_transinit

CORD_LONGITUDINAL_AXIS = 3
LONGITUDINAL_CORRESPONDENCE_JSON = "longitudinal_correspondence.json"


def default_longitudinal_correspondence_path(save_path: Path) -> Path:
    return save_path.expanduser() / LONGITUDINAL_CORRESPONDENCE_JSON


def load_longitudinal_correspondence(save_path: Path) -> SliceCorrespondence | None:
    path = default_longitudinal_correspondence_path(save_path)
    if not path.is_file():
        return None
    return SliceCorrespondence.load(path)


def estimate_cord_atlas_plane(
    sample_index: int,
    *,
    nslices: int,
    atlas_depth: int,
) -> int:
    """Map a 1-based straightened sample z index to an atlas plane using centered z-init."""
    transinit = build_cord_z_transinit(nslices, atlas_depth)
    z_scale = float(transinit[2, 2])
    z_trans = float(transinit[2, 3])
    atlas_z0 = (float(sample_index) - 1.0 - z_trans) / z_scale
    return int(np.clip(int(np.round(atlas_z0)) + 1, 1, atlas_depth))


def build_cord_z_transinit_from_correspondence(
    correspondence: SliceCorrespondence,
    *,
    atlas_depth: int,
) -> np.ndarray:
    """Fit z-scale + offset from confirmed longitudinal anchors (axis 3)."""
    anchors = correspondence.confirmed_anchors(CORD_LONGITUDINAL_AXIS)
    if len(anchors) < 2:
        msg = "Need at least two confirmed longitudinal anchors to fit z-init."
        raise ValueError(msg)

    atlas_z = np.array([anchor.atlas_plane - 1 for anchor in anchors], dtype=float)
    sample_z = np.array([anchor.sample_index - 1 for anchor in anchors], dtype=float)
    design = np.column_stack([atlas_z, np.ones_like(atlas_z)])
    z_scale, z_trans = np.linalg.lstsq(design, sample_z, rcond=None)[0]

    transinit = np.eye(4, dtype=float)
    transinit[2, 2] = float(z_scale)
    transinit[2, 3] = float(z_trans)
    return transinit


def resolve_cord_z_transinit(
    nslices: int,
    atlas_depth: int,
    correspondence: SliceCorrespondence | None,
) -> np.ndarray:
    """Return z-init from longitudinal correspondence when confirmed, else centered default."""
    if (
        correspondence is not None
        and correspondence.has_confirmed_anchors(CORD_LONGITUDINAL_AXIS)
    ):
        return build_cord_z_transinit_from_correspondence(
            correspondence,
            atlas_depth=atlas_depth,
        )
    return build_cord_z_transinit(nslices, atlas_depth)


def build_longitudinal_anchors(
    chooselist: np.ndarray,
    *,
    nslices: int,
    atlas_depth: int,
    confirmed: bool = False,
) -> list[SliceAnchor]:
    """Seed axis-3 anchors from centered z-init estimates."""
    anchors: list[SliceAnchor] = []
    for row in chooselist:
        sample_index = int(row[0])
        plane = estimate_cord_atlas_plane(
            sample_index,
            nslices=nslices,
            atlas_depth=atlas_depth,
        )
        anchors.append(
            SliceAnchor(
                sample_index=sample_index,
                atlas_plane=plane,
                confirmed=confirmed,
            )
        )
    return anchors
