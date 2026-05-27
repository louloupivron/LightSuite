"""Load volumes and atlas data for the brain control-point GUI."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import nibabel as nib
import numpy as np
from skimage.exposure import equalize_adapthist

from lightsuite.atlas.registry import resolve_brain_atlas
from lightsuite.config.models import BrainPipelineConfig
from lightsuite.gui.chooselist import generate_control_point_list
from lightsuite.gui.control_points import ControlPointSession, default_session_path
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


def load_brain_match_points_data(config: BrainPipelineConfig) -> BrainMatchPointsData:
    save_path = config.sample.save_path.expanduser()
    regopts_path = save_path / "regopts.json"
    checkpoint = RegOptsCheckpoint.load(regopts_path)
    if checkpoint.original_trans is None or checkpoint.permute_sample_to_atlas is None:
        msg = "Run init-registration before match-points."
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
    chooselist = generate_control_point_list(sample_warped.shape)

    session_path = default_session_path(save_path)
    if session_path.is_file():
        session = ControlPointSession.load(session_path)
    else:
        session = ControlPointSession.empty(original_trans, chooselist.shape[0])
        session.chooselist = chooselist.tolist()

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
    )


def slice_pair(data: BrainMatchPointsData, slice_idx: int) -> tuple[np.ndarray, np.ndarray]:
    row = data.chooselist[slice_idx - 1]
    sample = _normalize_display(volume_index_to_image(data.sample_volume, row))
    atlas = _normalize_display(volume_index_to_image(data.atlas_template, row))
    return sample, atlas


def prepare_brain_match_points_session(config: BrainPipelineConfig) -> Path:
    """Load match-point data and ensure session JSON exists (no GUI)."""
    data = load_brain_match_points_data(config)
    data.session.save(data.session_path)
    return data.session_path
