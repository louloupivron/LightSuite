"""Load volumes and atlas data for the brain control-point GUI."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import nibabel as nib
import numpy as np
from skimage.exposure import equalize_adapthist

from lightsuite.atlas.registry import resolve_brain_atlas
from lightsuite.config.models import BrainPipelineConfig
from lightsuite.gui.affine import transform_points_inverse
from lightsuite.gui.chooselist import (
    default_ap_cut_axis,
    generate_ap_alignment_list,
    generate_control_point_list,
)
from lightsuite.gui.control_points import ControlPointSession, default_session_path
from lightsuite.gui.slice_correspondence import (
    VOLUME_AXES,
    SliceAnchor,
    SliceCorrespondence,
    default_correspondence_path,
)
from lightsuite.gui.slices import volume_index_to_image
from lightsuite.preprocess.checkpoint import RegOptsCheckpoint
from lightsuite.registration.volume import (
    load_registration_volume,
    permute_brain_volume,
    resize_atlas_volume,
)


@dataclass
class BrainMatchPointsData:
    sample_volume: np.ndarray
    atlas_template: np.ndarray
    atlas_annotation: np.ndarray
    chooselist: np.ndarray
    session: ControlPointSession
    session_path: Path
    original_trans: np.ndarray
    auto_alignment: np.ndarray
    slice_correspondence: SliceCorrespondence | None = None


@dataclass
class BrainAlignSlicesData:
    sample_volume: np.ndarray
    atlas_template: np.ndarray
    chooselist_by_axis: dict[int, np.ndarray]
    original_trans: np.ndarray
    correspondence_path: Path
    correspondence: SliceCorrespondence
    axis_order: tuple[int, ...] = VOLUME_AXES

    def chooselist_for_axis(self, cut_axis: int) -> np.ndarray:
        return self.chooselist_by_axis[int(cut_axis)]


def _normalize_display(image: np.ndarray) -> np.ndarray:
    data = image.astype(np.float32)
    if data.max() <= 0:
        return data
    hi = float(np.quantile(data, 0.999))
    data = np.clip(data / max(hi, 1e-6), 0, 1)
    if data.ndim == 2:
        try:
            data = equalize_adapthist(data, clip_limit=0.01)
        except ValueError:
            pass
    return data


def _warp_sample_to_atlas_grid(
    volume: np.ndarray,
    original_trans: np.ndarray,
    target_shape: tuple[int, int, int],
) -> np.ndarray:
    """MATLAB imwarp(volload, Rvolume, original_trans, OutputView=Rmoving)."""
    from lightsuite.registration.warp import warp_sample_to_atlas

    return warp_sample_to_atlas(volume, original_trans, target_shape, order=1)


def _load_brain_registration_volumes(
    config: BrainPipelineConfig,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, float]:
    save_path = config.sample.save_path.expanduser()
    regopts_path = save_path / "regopts.json"
    checkpoint = RegOptsCheckpoint.load(regopts_path)
    if checkpoint.original_trans is None or checkpoint.permute_sample_to_atlas is None:
        msg = "Run init-registration before slice alignment or match-points."
        raise RuntimeError(msg)

    original_trans = np.asarray(checkpoint.original_trans, dtype=float)
    permvec = checkpoint.permute_sample_to_atlas
    downfac = float(checkpoint.downfac_reg or (config.atlas.resolution_um / checkpoint.registres_um))

    regvol = load_registration_volume(Path(checkpoint.regvolpath))
    regvol = permute_brain_volume(regvol.astype(np.float32), permvec)

    atlas = resolve_brain_atlas(config.atlas.provider.value, config.atlas.atlas_dir)
    tv = np.asanyarray(nib.load(atlas.template_path).dataobj).astype(np.float32)
    av = np.asanyarray(nib.load(atlas.annotation_path).dataobj).astype(np.float32)
    tvreg = resize_atlas_volume(tv, downfac, nearest=False)
    avreg = resize_atlas_volume(av, downfac, nearest=True)

    sample_warped = _warp_sample_to_atlas_grid(regvol, original_trans, tvreg.shape)
    return sample_warped, tvreg, avreg, original_trans, downfac


def _auto_atlas_plane_for_row(
    sample_volume: np.ndarray,
    chooserow: np.ndarray,
    original_trans: np.ndarray,
    atlas_shape: tuple[int, int, int],
) -> int:
    return estimate_atlas_plane_index(
        sample_volume,
        chooserow,
        original_trans,
        atlas_shape,
    )


def _build_correspondence_anchors(
    sample_volume: np.ndarray,
    chooselist: np.ndarray,
    original_trans: np.ndarray,
    atlas_shape: tuple[int, int, int],
    *,
    confirmed: bool,
) -> list[SliceAnchor]:
    anchors: list[SliceAnchor] = []
    for row in chooselist:
        chooserow = np.asarray(row, dtype=int)
        plane = _auto_atlas_plane_for_row(
            sample_volume,
            chooserow,
            original_trans,
            atlas_shape,
        )
        anchors.append(
            SliceAnchor(
                sample_index=int(chooserow[0]),
                atlas_plane=plane,
                confirmed=confirmed,
            )
        )
    return anchors


def load_slice_correspondence(save_path: Path) -> SliceCorrespondence | None:
    path = default_correspondence_path(save_path)
    if not path.is_file():
        return None
    return SliceCorrespondence.load(path)


def apply_slice_correspondence_to_session(
    session: ControlPointSession,
    chooselist: np.ndarray,
    correspondence: SliceCorrespondence,
    atlas_shape: tuple[int, int, int],
) -> None:
    """Pre-fill atlas_slice_indices from a saved correspondence curve."""
    n_slices = int(chooselist.shape[0])
    indices: list[int] = []
    for row in chooselist:
        chooserow = np.asarray(row, dtype=int)
        cut_axis = int(chooserow[1])
        atlas_size = atlas_cut_axis_size(atlas_shape, chooserow)
        plane = correspondence.interpolate_atlas_plane(
            int(chooserow[0]),
            cut_axis,
            atlas_size,
        )
        if plane is None:
            plane = 0
        indices.append(int(plane))
    session.atlas_slice_indices = indices


def _ensure_axis_correspondence(
    correspondence: SliceCorrespondence,
    *,
    cut_axis: int,
    sample_volume: np.ndarray,
    original_trans: np.ndarray,
    atlas_shape: tuple[int, int, int],
    chooselist: np.ndarray,
) -> None:
    anchors = correspondence.anchors_for_axis(cut_axis)
    if len(anchors) != chooselist.shape[0]:
        correspondence.axes[cut_axis] = _build_correspondence_anchors(
            sample_volume,
            chooselist,
            original_trans,
            atlas_shape,
            confirmed=False,
        )


def load_brain_align_slices_data(config: BrainPipelineConfig) -> BrainAlignSlicesData:
    sample_warped, tvreg, _avreg, original_trans, _downfac = _load_brain_registration_volumes(
        config
    )
    chooselist_by_axis = {
        axis: generate_ap_alignment_list(sample_warped.shape, cut_axis=axis)
        for axis in VOLUME_AXES
    }
    correspondence_path = default_correspondence_path(config.sample.save_path)

    if correspondence_path.is_file():
        correspondence = SliceCorrespondence.load(correspondence_path)
    else:
        correspondence = SliceCorrespondence(
            original_trans=original_trans.tolist(),
            axes={},
            source="auto",
        )

    correspondence.original_trans = original_trans.tolist()
    for axis, chooselist in chooselist_by_axis.items():
        _ensure_axis_correspondence(
            correspondence,
            cut_axis=axis,
            sample_volume=sample_warped,
            original_trans=original_trans,
            atlas_shape=tvreg.shape,
            chooselist=chooselist,
        )

    return BrainAlignSlicesData(
        sample_volume=sample_warped,
        atlas_template=tvreg,
        chooselist_by_axis=chooselist_by_axis,
        original_trans=original_trans,
        correspondence_path=correspondence_path,
        correspondence=correspondence,
    )


def prepare_brain_align_slices_session(config: BrainPipelineConfig) -> Path:
    """Auto-estimate and save slice correspondence on all axes without opening Napari."""
    data = load_brain_align_slices_data(config)
    for axis in data.axis_order:
        for anchor in data.correspondence.anchors_for_axis(axis):
            anchor.confirmed = True
    data.correspondence.source = "auto"
    data.correspondence.save(data.correspondence_path)
    return data.correspondence_path


def load_brain_match_points_data(config: BrainPipelineConfig) -> BrainMatchPointsData:
    save_path = config.sample.save_path.expanduser()
    sample_warped, tvreg, avreg, original_trans, _downfac = _load_brain_registration_volumes(
        config
    )
    chooselist = generate_control_point_list(sample_warped.shape)
    slice_correspondence = load_slice_correspondence(save_path)

    session_path = default_session_path(save_path)
    if session_path.is_file():
        session = ControlPointSession.load(session_path)
    else:
        session = ControlPointSession.empty(original_trans, chooselist.shape[0])
        session.chooselist = chooselist.tolist()

    if slice_correspondence is not None and slice_correspondence.has_confirmed_anchors():
        has_manual_planes = (
            session.atlas_slice_indices is not None
            and any(int(v) > 0 for v in session.atlas_slice_indices)
        )
        if not has_manual_planes:
            apply_slice_correspondence_to_session(
                session,
                chooselist,
                slice_correspondence,
                tvreg.shape,
            )

    auto_alignment = np.asarray(session.atlas2histology_tform, dtype=float)
    if np.allclose(auto_alignment, np.eye(4)):
        auto_alignment = original_trans.copy()

    return BrainMatchPointsData(
        sample_volume=sample_warped,
        atlas_template=tvreg,
        atlas_annotation=avreg,
        chooselist=chooselist,
        session=session,
        session_path=session_path,
        original_trans=original_trans,
        auto_alignment=auto_alignment,
        slice_correspondence=slice_correspondence,
    )


def atlas_cut_axis_size(atlas_shape: tuple[int, int, int], chooserow: np.ndarray) -> int:
    """Number of atlas planes along the chooselist cut axis (1-based indexing)."""
    cut_axis = int(chooserow[1]) - 1
    return int(atlas_shape[cut_axis])


def chooserow_with_atlas_plane(chooserow: np.ndarray, atlas_plane: int) -> np.ndarray:
    """Copy chooserow with atlas plane index along the cut axis (MATLAB atlas_slice)."""
    row = np.asarray(chooserow, dtype=int).copy()
    row[0] = int(atlas_plane)
    return row


def estimate_atlas_plane_index(
    sample_volume: np.ndarray,
    chooserow: np.ndarray,
    atlas_to_sample: np.ndarray,
    atlas_shape: tuple[int, int, int],
) -> int:
    """Estimate atlas plane along cut axis from inverse-transformed sample tissue (MATLAB update_slice)."""
    chooserow = np.asarray(chooserow, dtype=int)
    sample_slice = volume_index_to_image(sample_volume, chooserow)
    cut_axis = int(chooserow[1]) - 1
    plot_axes = [d for d in range(3) if d != cut_axis]
    h, w = sample_slice.shape
    tissue = sample_slice > float(np.quantile(sample_slice[sample_slice > 0], 0.05)) if np.any(
        sample_slice > 0
    ) else np.ones_like(sample_slice, dtype=bool)
    yy, xx = np.where(tissue)
    if yy.size == 0:
        yy, xx = np.mgrid[0:h, 0:w]
        yy, xx = yy.ravel(), xx.ravel()

    pts = np.zeros((yy.size, 3), dtype=float)
    pts[:, plot_axes[1]] = yy.astype(float) + 1.0
    pts[:, plot_axes[0]] = xx.astype(float) + 1.0
    pts[:, cut_axis] = float(chooserow[0])
    pts_fit = pts[:, [1, 0, 2]]
    atlas_pts = transform_points_inverse(pts_fit, np.asarray(atlas_to_sample, dtype=float))
    atlas_coords = atlas_pts[:, [1, 0, 2]]
    sluse = int(np.round(float(np.median(atlas_coords[:, cut_axis]))))
    return int(np.clip(sluse, 1, atlas_cut_axis_size(atlas_shape, chooserow)))


def resolve_atlas_plane_index(
    data: BrainMatchPointsData,
    slice_idx: int,
) -> int:
    """Atlas plane for one chooselist entry (manual, from points, or auto-estimated)."""
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

    if data.slice_correspondence is not None:
        atlas_size = atlas_cut_axis_size(data.atlas_template.shape, chooserow)
        plane = data.slice_correspondence.interpolate_atlas_plane(
            int(chooserow[0]),
            int(chooserow[1]),
            atlas_size,
        )
        if plane is not None:
            return plane

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
    data: BrainMatchPointsData,
    slice_idx: int,
    *,
    atlas_plane: int | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    row = np.asarray(data.chooselist[slice_idx - 1], dtype=int)
    sample = _normalize_display(volume_index_to_image(data.sample_volume, row))
    plane = atlas_plane if atlas_plane is not None else resolve_atlas_plane_index(data, slice_idx)
    atlas_row = chooserow_with_atlas_plane(row, plane)
    atlas = _normalize_display(volume_index_to_image(data.atlas_template, atlas_row))
    return sample, atlas


def prepare_brain_match_points_session(config: BrainPipelineConfig) -> Path:
    """Load match-point data and ensure session JSON exists (no GUI)."""
    data = load_brain_match_points_data(config)
    data.session.save(data.session_path)
    return data.session_path
