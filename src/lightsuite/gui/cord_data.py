"""Load straightened sample and warped atlas data for spinal cord match-points."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import tifffile
from scipy import ndimage

from lightsuite.config.models import SpinalCordPipelineConfig
from lightsuite.gui.cord_align_data import apply_longitudinal_correspondence_to_session
from lightsuite.gui.cord_display import CORD_ATLAS_PROVIDER, normalize_cord_display
from lightsuite.gui.chooselist import generate_cord_control_point_list
from lightsuite.gui.control_points import ControlPointSession
from lightsuite.gui.slice_correspondence import SliceCorrespondence
from lightsuite.gui.slices import prepare_display_slice, volume_index_to_image
from lightsuite.preprocess.cord_checkpoint import CordRegOptsCheckpoint
from lightsuite.registration.cord_affine import warp_cord_atlas_to_straightvol
from lightsuite.registration.cord_longitudinal import (
    CORD_LONGITUDINAL_AXIS,
    load_longitudinal_correspondence,
    resolve_cord_z_transinit,
)
from lightsuite.registration.cord_paths import (
    cord_affine_transform_path,
    cord_save_path,
    cord_work_dir,
)

CORRESPONDING_POINTS_JSON = "corresponding_points.json"


def _normalize_display(image: np.ndarray) -> np.ndarray:
    return normalize_cord_display(image)


@dataclass
class CordMatchPointsData:
    sample_volume: np.ndarray
    atlas_template: np.ndarray
    atlas_annotation: np.ndarray
    chooselist: np.ndarray
    session: ControlPointSession
    session_path: Path
    affine_atlas_to_samp: np.ndarray
    longitudinal_correspondence: SliceCorrespondence | None = None


def default_cord_session_path(save_path: Path) -> Path:
    return save_path.expanduser() / CORRESPONDING_POINTS_JSON


def load_cord_match_points_data(config: SpinalCordPipelineConfig) -> CordMatchPointsData:
    """Load straightvol and atlas volumes warped into the same straightened grid."""
    save_path = cord_save_path(config)
    regopts_path = save_path / "regopts.json"
    if not regopts_path.is_file():
        msg = f"Missing {regopts_path}. Run init-registration first."
        raise FileNotFoundError(msg)

    checkpoint = CordRegOptsCheckpoint.load(regopts_path)
    if checkpoint.straightvol_path is None or checkpoint.affine_atlas_to_samp is None:
        msg = "regopts.json missing straightvol / affine_atlas_to_samp from init-registration."
        raise RuntimeError(msg)

    straightvol = tifffile.imread(checkpoint.straightvol_path).astype(np.float32)
    tv = tifffile.imread(checkpoint.tv_path).astype(np.float32)
    av = tifffile.imread(checkpoint.av_path).astype(np.uint16)
    affine_atlas_to_samp = np.asarray(checkpoint.affine_atlas_to_samp, dtype=float)
    nslices = checkpoint.ikeeprange[1] - checkpoint.ikeeprange[0] + 1
    correspondence = load_longitudinal_correspondence(save_path)

    transinit = resolve_cord_z_transinit(nslices, tv.shape[2], correspondence)
    elastix_affine_path = cord_affine_transform_path(config)
    spacing_mm = config.registration.resolution_um * 1e-3
    output_shape = straightvol.shape

    tvtemp = ndimage.median_filter(tv, size=3)
    atlas_template = warp_cord_atlas_to_straightvol(
        tvtemp,
        transinit=transinit,
        elastix_affine_path=elastix_affine_path,
        output_shape=output_shape,
        spacing_mm=spacing_mm,
        work_dir=cord_work_dir(config, "transformix", "match", "tv"),
        nearest=False,
    )
    atlas_annotation = warp_cord_atlas_to_straightvol(
        av,
        transinit=transinit,
        elastix_affine_path=elastix_affine_path,
        output_shape=output_shape,
        spacing_mm=spacing_mm,
        work_dir=cord_work_dir(config, "transformix", "match", "av"),
        nearest=True,
    )

    chooselist = generate_cord_control_point_list(output_shape[2])
    n_slices = int(chooselist.shape[0])
    session_path = default_cord_session_path(save_path)
    if session_path.is_file():
        session = ControlPointSession.load(session_path)
    else:
        session = ControlPointSession.empty(np.eye(4), n_slices)
        session.ori_trans = affine_atlas_to_samp.tolist()

    session.chooselist = chooselist.tolist()
    if len(session.histology_control_points) != n_slices:
        session.histology_control_points = [
            session.histology_control_points[i] if i < len(session.histology_control_points) else []
            for i in range(n_slices)
        ]
        session.atlas_control_points = [
            session.atlas_control_points[i] if i < len(session.atlas_control_points) else []
            for i in range(n_slices)
        ]

    if correspondence is not None and correspondence.has_confirmed_anchors(CORD_LONGITUDINAL_AXIS):
        has_manual_planes = (
            session.atlas_slice_indices is not None
            and any(int(v) > 0 for v in session.atlas_slice_indices)
        )
        if not has_manual_planes:
            apply_longitudinal_correspondence_to_session(
                session,
                chooselist,
                correspondence,
                tv.shape,
            )

    return CordMatchPointsData(
        sample_volume=straightvol,
        atlas_template=atlas_template,
        atlas_annotation=atlas_annotation,
        chooselist=chooselist,
        session=session,
        session_path=session_path,
        affine_atlas_to_samp=affine_atlas_to_samp,
        longitudinal_correspondence=correspondence,
    )


def atlas_cut_axis_size(atlas_shape: tuple[int, int, int], chooserow: np.ndarray) -> int:
    """Re-export for backward compatibility."""
    from lightsuite.gui.cord_display import atlas_cut_axis_size as _size

    return _size(atlas_shape, chooserow)


def chooserow_with_atlas_plane(chooserow: np.ndarray, atlas_plane: int) -> np.ndarray:
    row = np.asarray(chooserow, dtype=int).copy()
    row[0] = int(atlas_plane)
    return row


def resolve_atlas_plane_index(data: CordMatchPointsData, slice_idx: int) -> int:
    """Atlas plane for one chooselist entry (manual override, points, or auto-estimate)."""
    from lightsuite.gui.brain_data import estimate_atlas_plane_index

    chooserow = np.asarray(data.chooselist[slice_idx - 1], dtype=int)
    atlas_points = data.session.atlas_control_points[slice_idx - 1]
    if atlas_points:
        cut_axis = int(chooserow[1]) - 1
        return int(
            np.clip(
                int(np.round(float(np.median([p[cut_axis] for p in atlas_points])))),
                1,
                atlas_cut_axis_size(data.atlas_template.shape, chooserow),
            )
        )

    stored = data.session.atlas_slice_indices
    if stored is not None and len(stored) >= slice_idx and int(stored[slice_idx - 1]) > 0:
        return int(
            np.clip(
                int(stored[slice_idx - 1]),
                1,
                atlas_cut_axis_size(data.atlas_template.shape, chooserow),
            )
        )

    if data.longitudinal_correspondence is not None:
        plane = data.longitudinal_correspondence.interpolate_atlas_plane(
            int(chooserow[0]),
            CORD_LONGITUDINAL_AXIS,
            atlas_cut_axis_size(data.atlas_template.shape, chooserow),
        )
        if plane is not None:
            return int(plane)

    matrix = np.asarray(data.session.atlas2histology_tform, dtype=float)
    return estimate_atlas_plane_index(
        data.sample_volume,
        chooserow,
        matrix,
        data.atlas_template.shape,
    )


def set_atlas_plane_index(session: ControlPointSession, slice_idx: int, plane: int) -> None:
    n_slices = len(session.histology_control_points)
    if session.atlas_slice_indices is None or len(session.atlas_slice_indices) != n_slices:
        session.atlas_slice_indices = [0] * n_slices
    session.atlas_slice_indices[slice_idx - 1] = int(plane)


def slice_pair(
    data: CordMatchPointsData,
    slice_idx: int,
    *,
    atlas_plane: int | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    row = np.asarray(data.chooselist[slice_idx - 1], dtype=int)
    cut_axis = int(row[1])
    sample = _normalize_display(
        prepare_display_slice(
            volume_index_to_image(data.sample_volume, row),
            cut_axis,
            CORD_ATLAS_PROVIDER,
        )
    )
    plane = atlas_plane if atlas_plane is not None else resolve_atlas_plane_index(data, slice_idx)
    atlas_row = chooserow_with_atlas_plane(row, plane)
    atlas = _normalize_display(
        prepare_display_slice(
            volume_index_to_image(data.atlas_template, atlas_row),
            cut_axis,
            CORD_ATLAS_PROVIDER,
        )
    )
    return sample, atlas


def prepare_cord_match_points_session(config: SpinalCordPipelineConfig) -> Path:
    """Ensure corresponding_points.json exists (headless helper)."""
    data = load_cord_match_points_data(config)
    data.session.save(data.session_path)
    return data.session_path
