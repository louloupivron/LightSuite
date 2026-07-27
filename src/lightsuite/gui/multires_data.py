"""Load multiresolution manifest volumes for the landmark match-points GUI."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from lightsuite.multires.config_models import MultiresGeometryMode, MultiresPipelineConfig
from lightsuite.multires.landmark_session import MultiresLandmarkSession
from lightsuite.multires.manifest import load_pair_manifest
from lightsuite.multires.models import ManifestVolumeSpec, MultiresPairManifest
from lightsuite.multires.spec_geometry import (
    crop_index_range_from_physical_box,
    index_xyz_to_physical,
    overlap_physical_bounds_from_specs,
    physical_to_continuous_index_xyz,
)
from lightsuite.multires.volume import load_manifest_xy_crop

_SLICE_CACHE_MAX = 32

# Extra margin around metadata overlap for landmark placement (µm).
DEFAULT_MATCH_POINTS_MARGIN_UM = 200.0


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


def _expand_crop(
    start_xyz: list[int],
    crop_size_xyz: list[int],
    *,
    volume_shape_zyx: tuple[int, int, int],
    margin_xy_vox: int,
    margin_z_vox: int,
) -> tuple[list[int], list[int]]:
    """Expand an XYZ crop; clamp to volume bounds."""
    nz, ny, nx = volume_shape_zyx
    ix0, iy0, iz0 = start_xyz
    sx, sy, sz = crop_size_xyz
    ix0 = max(0, ix0 - margin_xy_vox)
    iy0 = max(0, iy0 - margin_xy_vox)
    iz0 = max(0, iz0 - margin_z_vox)
    ix1 = min(nx, ix0 + sx + 2 * margin_xy_vox)
    iy1 = min(ny, iy0 + sy + 2 * margin_xy_vox)
    iz1 = min(nz, iz0 + sz + 2 * margin_z_vox)
    return [ix0, iy0, iz0], [ix1 - ix0, iy1 - iy0, iz1 - iz0]


@dataclass
class MultiresSliceSource:
    """Lazy XY(-crop) plane reader with a small display cache.

    Landmark points remain in full-volume ZYX indices. Display coordinates are
    local to the XY crop window (``xy_origin_yx``).
    """

    spec: ManifestVolumeSpec
    manifest_dir: Path
    volume_shape_zyx: tuple[int, int, int]
    xy_origin_yx: tuple[int, int]
    crop_size_yx: tuple[int, int]
    z_min: int
    z_max: int
    _display_hi: float = 1.0
    _cache: dict[int, np.ndarray] = field(default_factory=dict, repr=False)

    @classmethod
    def from_crop(
        cls,
        spec: ManifestVolumeSpec,
        *,
        manifest_dir: Path,
        start_xyz: list[int],
        crop_size_xyz: list[int],
    ) -> MultiresSliceSource:
        nz, ny, nx = (int(v) for v in spec.shape_zyx)
        ix0, iy0, iz0 = (int(v) for v in start_xyz)
        sx, sy, sz = (int(v) for v in crop_size_xyz)
        z_min = int(np.clip(iz0, 0, max(nz - 1, 0)))
        z_max = int(np.clip(iz0 + max(sz, 1) - 1, z_min, max(nz - 1, 0)))
        source = cls(
            spec=spec,
            manifest_dir=manifest_dir,
            volume_shape_zyx=(nz, ny, nx),
            xy_origin_yx=(iy0, ix0),
            crop_size_yx=(sy, sx),
            z_min=z_min,
            z_max=z_max,
        )
        mid_z = int(np.clip((z_min + z_max) // 2, z_min, z_max))
        source._display_hi = _compute_display_hi(source._read_plane(mid_z))
        return source

    @classmethod
    def from_spec(
        cls,
        spec: ManifestVolumeSpec,
        *,
        manifest_dir: Path,
    ) -> MultiresSliceSource:
        """Full-volume fallback (no XY crop)."""
        nz, ny, nx = (int(v) for v in spec.shape_zyx)
        return cls.from_crop(
            spec,
            manifest_dir=manifest_dir,
            start_xyz=[0, 0, 0],
            crop_size_xyz=[nx, ny, nz],
        )

    @property
    def shape_zyx(self) -> tuple[int, int, int]:
        """Display shape: Z span of the crop window × crop YX."""
        return (self.z_max - self.z_min + 1, self.crop_size_yx[0], self.crop_size_yx[1])

    @property
    def is_xy_cropped(self) -> bool:
        _nz, ny, nx = self.volume_shape_zyx
        return self.crop_size_yx != (ny, nx) or self.xy_origin_yx != (0, 0)

    def clip_z(self, z_index: int) -> int:
        return int(np.clip(z_index, self.z_min, self.z_max))

    def _read_plane(self, z_index: int) -> np.ndarray:
        iy0, ix0 = self.xy_origin_yx
        sy, sx = self.crop_size_yx
        return load_manifest_xy_crop(
            self.spec,
            z_index=z_index,
            start_xyz=[ix0, iy0, 0],
            crop_size_xyz=[sx, sy, 1],
            manifest_dir=self.manifest_dir,
        )

    def read_display_slice(self, z_index: int) -> np.ndarray:
        """Return one normalized YX crop for Napari display."""
        z = self.clip_z(z_index)
        cached = self._cache.get(z)
        if cached is not None:
            return cached
        display = _normalize_display_fast(self._read_plane(z), self._display_hi)
        if len(self._cache) >= _SLICE_CACHE_MAX:
            self._cache.pop(next(iter(self._cache)))
        self._cache[z] = display
        return display

    def display_xy_from_volume_zyx(
        self,
        points_zyx: list[list[float]],
        z_index: int,
    ) -> np.ndarray:
        """Map full-volume ZYX points on this Z to crop-local (Y, X)."""
        if not points_zyx:
            return np.zeros((0, 2))
        iy0, ix0 = self.xy_origin_yx
        rows: list[list[float]] = []
        for pt in points_zyx:
            z, y, x = (float(v) for v in pt[:3])
            if int(round(z)) != int(z_index):
                continue
            rows.append([y - iy0, x - ix0])
        if not rows:
            return np.zeros((0, 2))
        return np.asarray(rows, dtype=float)

    def volume_zyx_from_display_xy(
        self,
        layer_xy: np.ndarray,
        z_index: int,
        existing_points: list[list[float]],
    ) -> list[list[float]]:
        """Replace points on this Z; store full-volume ZYX indices."""
        iy0, ix0 = self.xy_origin_yx
        kept = [pt for pt in existing_points if int(round(float(pt[0]))) != int(z_index)]
        new_points = [
            [float(z_index), float(row) + iy0, float(col) + ix0]
            for row, col in np.asarray(layer_xy, dtype=float)
        ]
        return kept + new_points

    def physical_z_um(self, z_index: int) -> float:
        z = self.clip_z(z_index)
        iy0, ix0 = self.xy_origin_yx
        return float(index_xyz_to_physical(self.spec, (float(ix0), float(iy0), float(z)))[2])

    def z_index_from_physical_z(self, z_um: float) -> int:
        iy0, ix0 = self.xy_origin_yx
        xy_phys = index_xyz_to_physical(self.spec, (float(ix0), float(iy0), float(self.z_min)))
        point = (float(xy_phys[0]), float(xy_phys[1]), float(z_um))
        z = int(round(physical_to_continuous_index_xyz(self.spec, point)[2]))
        return self.clip_z(z)


@dataclass
class MultiresMatchPointsData:
    overview: MultiresSliceSource
    roi: MultiresSliceSource
    overview_spec: ManifestVolumeSpec
    roi_spec: ManifestVolumeSpec
    session: MultiresLandmarkSession
    session_path: Path
    fit_mode: str
    min_pairs: int
    initial_overview_z: int
    initial_roi_z: int
    crop_mode: bool
    margin_um: float

    @property
    def overview_shape_zyx(self) -> tuple[int, int, int]:
        return self.overview.volume_shape_zyx

    @property
    def roi_shape_zyx(self) -> tuple[int, int, int]:
        return self.roi.volume_shape_zyx


def _default_z_indices(
    overview_spec: ManifestVolumeSpec,
    roi_spec: ManifestVolumeSpec,
) -> tuple[int, int]:
    """Prefer overlap-center Z when metadata FOVs intersect; else mid-volume."""
    overview_mid = int(overview_spec.shape_zyx[0]) // 2
    roi_mid = int(roi_spec.shape_zyx[0]) // 2
    try:
        overlap_min, overlap_max = overlap_physical_bounds_from_specs(
            overview_spec,
            roi_spec,
            margin_um=0.0,
        )
    except ValueError:
        return overview_mid, roi_mid

    center = 0.5 * (np.asarray(overlap_min, dtype=float) + np.asarray(overlap_max, dtype=float))
    oz = int(round(physical_to_continuous_index_xyz(overview_spec, tuple(center))[2]))
    rz = int(round(physical_to_continuous_index_xyz(roi_spec, tuple(center))[2]))
    oz = int(np.clip(oz, 0, overview_spec.shape_zyx[0] - 1))
    rz = int(np.clip(rz, 0, roi_spec.shape_zyx[0] - 1))
    return oz, rz


def _match_points_margin_um(cfg: MultiresPipelineConfig) -> float:
    return max(float(cfg.multires.registration.overlap_margin_um), DEFAULT_MATCH_POINTS_MARGIN_UM)


def _overlap_crops_from_metadata(
    overview_spec: ManifestVolumeSpec,
    roi_spec: ManifestVolumeSpec,
    *,
    margin_um: float,
) -> tuple[list[int], list[int], list[int], list[int]] | None:
    """Return overview/ROI XYZ crops around metadata overlap, or None if no overlap."""
    try:
        # Tight FOV intersection; XY margin is applied once via voxel expansion below.
        overlap_min, overlap_max = overlap_physical_bounds_from_specs(
            overview_spec,
            roi_spec,
            margin_um=0.0,
        )
    except ValueError:
        return None

    overview_start, overview_size = crop_index_range_from_physical_box(
        overview_spec,
        overlap_min,
        overlap_max,
    )
    roi_start, roi_size = crop_index_range_from_physical_box(
        roi_spec,
        overlap_min,
        overlap_max,
    )
    if min(overview_size) <= 0 or min(roi_size) <= 0:
        return None

    onz, ony, onx = (int(v) for v in overview_spec.shape_zyx)
    rnz, rny, rnx = (int(v) for v in roi_spec.shape_zyx)
    margin_xy_ov = max(1, int(round(margin_um / max(float(overview_spec.spacing_um[0]), 1e-6))))
    margin_xy_roi = max(1, int(round(margin_um / max(float(roi_spec.spacing_um[0]), 1e-6))))
    margin_z_ov = max(1, int(round(margin_um / max(float(overview_spec.spacing_um[2]), 1e-6))))
    margin_z_roi = max(1, int(round(margin_um / max(float(roi_spec.spacing_um[2]), 1e-6))))
    overview_start, overview_size = _expand_crop(
        overview_start,
        overview_size,
        volume_shape_zyx=(onz, ony, onx),
        margin_xy_vox=margin_xy_ov,
        margin_z_vox=margin_z_ov,
    )
    roi_start, roi_size = _expand_crop(
        roi_start,
        roi_size,
        volume_shape_zyx=(rnz, rny, rnx),
        margin_xy_vox=margin_xy_roi,
        margin_z_vox=margin_z_roi,
    )
    return overview_start, overview_size, roi_start, roi_size


def load_multires_match_points_data(
    cfg: MultiresPipelineConfig,
    *,
    manifest: MultiresPairManifest | None = None,
) -> MultiresMatchPointsData:
    """Load overview / ROI crop readers and the landmark session.

    Uses metadata FOV overlap (hybrid placement) to load only the overview crop
    around the predicted ROI location, plus the matching ROI crop — not the full
    overview mosaic.
    """
    if manifest is None:
        manifest = load_pair_manifest(cfg.multires.pair_manifest)
    manifest_dir = cfg.multires.pair_manifest.parent

    session_path = cfg.multires.resolved_landmark_session_path(cfg.sample.save_path, manifest)
    if session_path.is_file():
        session = MultiresLandmarkSession.load(session_path)
    else:
        session = MultiresLandmarkSession(fit_mode=cfg.multires.landmarks.fit_mode)

    session.fit_mode = cfg.multires.landmarks.fit_mode
    margin_um = _match_points_margin_um(cfg)
    crops = None
    if cfg.multires.geometry_mode == MultiresGeometryMode.HYBRID:
        crops = _overlap_crops_from_metadata(
            manifest.overview,
            manifest.roi,
            margin_um=margin_um,
        )

    if crops is not None:
        overview_start, overview_size, roi_start, roi_size = crops
        overview = MultiresSliceSource.from_crop(
            manifest.overview,
            manifest_dir=manifest_dir,
            start_xyz=overview_start,
            crop_size_xyz=overview_size,
        )
        roi = MultiresSliceSource.from_crop(
            manifest.roi,
            manifest_dir=manifest_dir,
            start_xyz=roi_start,
            crop_size_xyz=roi_size,
        )
        crop_mode = True
        oz, rz = _default_z_indices(manifest.overview, manifest.roi)
        overview_z = overview.clip_z(oz)
        roi_z = roi.clip_z(rz)
    else:
        overview = MultiresSliceSource.from_spec(manifest.overview, manifest_dir=manifest_dir)
        roi = MultiresSliceSource.from_spec(manifest.roi, manifest_dir=manifest_dir)
        crop_mode = False
        overview_z, roi_z = _default_z_indices(manifest.overview, manifest.roi)

    return MultiresMatchPointsData(
        overview=overview,
        roi=roi,
        overview_spec=manifest.overview,
        roi_spec=manifest.roi,
        session=session,
        session_path=session_path,
        fit_mode=cfg.multires.landmarks.fit_mode,
        min_pairs=cfg.multires.landmarks.min_pairs,
        initial_overview_z=overview_z,
        initial_roi_z=roi_z,
        crop_mode=crop_mode,
        margin_um=margin_um,
    )


def prepare_multires_match_points_session(cfg: MultiresPipelineConfig) -> Path:
    """Ensure landmark session exists (headless helper)."""
    data = load_multires_match_points_data(cfg)
    data.session.save(data.session_path)
    return data.session_path
