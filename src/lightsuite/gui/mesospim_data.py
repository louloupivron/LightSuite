"""Load mesoSPIM volumes for the landmark match-points GUI."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from lightsuite.mesospim.config_models import MesospimPipelineConfig, MesospimTiffRemapConfig
from lightsuite.mesospim.io import read_tiff_xy_slice_at_z_index, tiff_shape
from lightsuite.mesospim.landmark_session import (
    MesospimLandmarkSession,
)

_SLICE_CACHE_MAX = 32


def _compute_display_hi(image: np.ndarray) -> float:
    data = image.astype(np.float32, copy=False)
    if data.max() <= 0:
        return 1.0
    positive = data[data > 0]
    if positive.size == 0:
        return float(data.max())
    return float(np.quantile(positive, 0.995))


def _normalize_display_fast(image: np.ndarray, hi: float) -> np.ndarray:
    data = image.astype(np.float32, copy=False)
    return np.clip(data / max(hi, 1e-6), 0, 1)


@dataclass
class MesospimSliceSource:
    """Lazy per-plane reader with a small display cache for GUI navigation."""

    path: Path
    overview_path: Path
    roi_path: Path
    remap: MesospimTiffRemapConfig
    shape_zyx: tuple[int, int, int]
    _display_hi: float = 1.0
    _cache: dict[int, np.ndarray] = field(default_factory=dict, repr=False)

    @classmethod
    def from_path(
        cls,
        path: Path,
        *,
        overview_path: Path,
        roi_path: Path,
        remap: MesospimTiffRemapConfig,
    ) -> MesospimSliceSource:
        resolved = path.expanduser().resolve()
        shape_zyx = tiff_shape(resolved)
        source = cls(
            path=resolved,
            overview_path=overview_path.expanduser().resolve(),
            roi_path=roi_path.expanduser().resolve(),
            remap=remap,
            shape_zyx=shape_zyx,
        )
        mid_z = shape_zyx[0] // 2
        source._display_hi = _compute_display_hi(source._read_remapped_plane(mid_z))
        return source

    def _read_remapped_plane(self, z_index: int) -> np.ndarray:
        return read_tiff_xy_slice_at_z_index(
            self.path,
            z_index,
            overview_path=self.overview_path,
            roi_path=self.roi_path,
            remap=self.remap,
        )

    def read_display_slice(self, z_index: int) -> np.ndarray:
        """Return one normalized YX plane for Napari display."""
        z = int(np.clip(z_index, 0, self.shape_zyx[0] - 1))
        cached = self._cache.get(z)
        if cached is not None:
            return cached
        display = _normalize_display_fast(self._read_remapped_plane(z), self._display_hi)
        if len(self._cache) >= _SLICE_CACHE_MAX:
            self._cache.pop(next(iter(self._cache)))
        self._cache[z] = display
        return display


@dataclass
class MesospimMatchPointsData:
    overview: MesospimSliceSource
    roi: MesospimSliceSource
    session: MesospimLandmarkSession
    session_path: Path
    fit_mode: str
    min_pairs: int

    @property
    def overview_shape_zyx(self) -> tuple[int, int, int]:
        return self.overview.shape_zyx

    @property
    def roi_shape_zyx(self) -> tuple[int, int, int]:
        return self.roi.shape_zyx


def load_mesospim_match_points_data(cfg: MesospimPipelineConfig) -> MesospimMatchPointsData:
    """Load overview / ROI slice readers and the landmark session."""
    meso = cfg.mesospim
    overview_path = meso.overview.path
    roi_path = meso.roi.path

    session_path = meso.resolved_landmark_session_path(cfg.sample.save_path)
    if session_path.is_file():
        session = MesospimLandmarkSession.load(session_path)
    else:
        session = MesospimLandmarkSession(fit_mode=meso.landmarks.fit_mode)

    session.fit_mode = meso.landmarks.fit_mode
    return MesospimMatchPointsData(
        overview=MesospimSliceSource.from_path(
            overview_path,
            overview_path=overview_path,
            roi_path=roi_path,
            remap=meso.tiff_remap,
        ),
        roi=MesospimSliceSource.from_path(
            roi_path,
            overview_path=overview_path,
            roi_path=roi_path,
            remap=meso.tiff_remap,
        ),
        session=session,
        session_path=session_path,
        fit_mode=meso.landmarks.fit_mode,
        min_pairs=meso.landmarks.min_pairs,
    )


def prepare_mesospim_match_points_session(cfg: MesospimPipelineConfig) -> Path:
    """Ensure landmark session exists (headless helper)."""
    data = load_mesospim_match_points_data(cfg)
    data.session.save(data.session_path)
    return data.session_path
