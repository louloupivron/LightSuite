"""Load straightened sample and atlas data for spinal cord longitudinal alignment."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import tifffile
from scipy import ndimage

from lightsuite.config.models import SpinalCordPipelineConfig
from lightsuite.gui.chooselist import generate_cord_longitudinal_list
from lightsuite.gui.cord_display import CORD_ATLAS_PROVIDER, normalize_cord_display
from lightsuite.gui.slices import prepare_display_slice, volume_index_to_image
from lightsuite.gui.slice_correspondence import SliceCorrespondence
from lightsuite.preprocess.cord_checkpoint import CordRegOptsCheckpoint, SpinalAlignmentCheckpoint
from lightsuite.registration.cord_longitudinal import (
    CORD_LONGITUDINAL_AXIS,
    build_longitudinal_anchors,
    default_longitudinal_correspondence_path,
    estimate_cord_atlas_plane,
    load_longitudinal_correspondence,
)
from lightsuite.registration.cord_paths import cord_save_path
from lightsuite.registration.straightening import (
    compute_straightening_transforms,
    transform_cord_images_slices,
)


@dataclass
class CordAlignLongitudinalData:
    sample_volume: np.ndarray
    atlas_template: np.ndarray
    chooselist: np.ndarray
    correspondence_path: Path
    correspondence: SliceCorrespondence
    nslices: int
    atlas_depth: int


def _load_straightened_sample(config: SpinalCordPipelineConfig) -> tuple[np.ndarray, int]:
    save_path = cord_save_path(config)
    regopts_path = save_path / "regopts.json"
    align_path = save_path / "spinal_alignment_opt.json"
    if not regopts_path.is_file():
        msg = f"Missing {regopts_path}. Run 'lightsuite spinal preprocess' first."
        raise FileNotFoundError(msg)
    if not align_path.is_file():
        msg = f"Missing {align_path}. Run 'lightsuite spinal straighten' first."
        raise FileNotFoundError(msg)

    checkpoint = CordRegOptsCheckpoint.load(regopts_path)
    align = SpinalAlignmentCheckpoint.load(align_path)
    regvol = tifffile.imread(checkpoint.regvol_path).astype(np.uint16)
    tv = tifffile.imread(checkpoint.tv_path).astype(np.float32)

    target_center = (0.5 * tv.shape[1], 0.5 * tv.shape[0])
    tforms = compute_straightening_transforms(
        np.array(align.fit_x),
        np.array(align.fit_y),
        np.array(align.fit_theta),
        target_center,
        config.registration.target_orientation_deg,
    )
    straightvol = transform_cord_images_slices(regvol, tforms, (tv.shape[0], tv.shape[1]))
    return straightvol.astype(np.float32), int(straightvol.shape[2])


def _ensure_longitudinal_correspondence(
    correspondence: SliceCorrespondence,
    *,
    chooselist: np.ndarray,
    nslices: int,
    atlas_depth: int,
) -> None:
    anchors = correspondence.anchors_for_axis(CORD_LONGITUDINAL_AXIS)
    if len(anchors) != chooselist.shape[0]:
        correspondence.axes[CORD_LONGITUDINAL_AXIS] = build_longitudinal_anchors(
            chooselist,
            nslices=nslices,
            atlas_depth=atlas_depth,
            confirmed=False,
        )


def load_cord_align_longitudinal_data(
    config: SpinalCordPipelineConfig,
) -> CordAlignLongitudinalData:
    """Load straightened sample, atlas template, and longitudinal correspondence anchors."""
    save_path = cord_save_path(config)
    straightvol, nslices = _load_straightened_sample(config)
    checkpoint = CordRegOptsCheckpoint.load(save_path / "regopts.json")
    tvtemp = ndimage.median_filter(
        tifffile.imread(checkpoint.tv_path).astype(np.float32),
        size=3,
    )
    atlas_depth = int(tvtemp.shape[2])
    chooselist = generate_cord_longitudinal_list(nslices)
    correspondence_path = default_longitudinal_correspondence_path(save_path)

    if correspondence_path.is_file():
        correspondence = SliceCorrespondence.load(correspondence_path)
    else:
        correspondence = SliceCorrespondence(
            original_trans=np.eye(4).tolist(),
            axes={},
            source="auto",
        )

    correspondence.original_trans = np.eye(4).tolist()
    _ensure_longitudinal_correspondence(
        correspondence,
        chooselist=chooselist,
        nslices=nslices,
        atlas_depth=atlas_depth,
    )

    return CordAlignLongitudinalData(
        sample_volume=straightvol,
        atlas_template=tvtemp,
        chooselist=chooselist,
        correspondence_path=correspondence_path,
        correspondence=correspondence,
        nslices=nslices,
        atlas_depth=atlas_depth,
    )


def align_longitudinal_pair(
    data: CordAlignLongitudinalData,
    slice_idx: int,
    *,
    atlas_plane: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Return sample and atlas transverse slices for one anchor."""
    row = np.asarray(data.chooselist[slice_idx - 1], dtype=int)
    cut_axis = int(row[1])
    sample = normalize_cord_display(
        prepare_display_slice(
            volume_index_to_image(data.sample_volume, row),
            cut_axis,
            CORD_ATLAS_PROVIDER,
        )
    )
    atlas_row = row.copy()
    atlas_row[0] = int(atlas_plane)
    atlas = normalize_cord_display(
        prepare_display_slice(
            volume_index_to_image(data.atlas_template, atlas_row),
            cut_axis,
            CORD_ATLAS_PROVIDER,
        )
    )
    return sample, atlas


def prepare_cord_align_longitudinal_session(config: SpinalCordPipelineConfig) -> Path:
    """Auto-confirm centered z-init anchors without opening Napari (tests / headless)."""
    data = load_cord_align_longitudinal_data(config)
    for anchor in data.correspondence.anchors_for_axis(CORD_LONGITUDINAL_AXIS):
        anchor.confirmed = True
    data.correspondence.source = "auto"
    data.correspondence.save(data.correspondence_path)
    return data.correspondence_path


def apply_longitudinal_correspondence_to_session(
    session,
    chooselist: np.ndarray,
    correspondence: SliceCorrespondence,
    atlas_shape: tuple[int, int, int],
) -> None:
    """Pre-fill atlas_slice_indices from saved longitudinal correspondence."""
    from lightsuite.gui.cord_display import atlas_cut_axis_size

    indices: list[int] = []
    for row in chooselist:
        chooserow = np.asarray(row, dtype=int)
        atlas_size = atlas_cut_axis_size(atlas_shape, chooserow)
        plane = correspondence.interpolate_atlas_plane(
            int(chooserow[0]),
            CORD_LONGITUDINAL_AXIS,
            atlas_size,
        )
        indices.append(int(plane or 0))
    session.atlas_slice_indices = indices
