"""Spinal cord control-point matching GUI (matchControlPointsSpine.m port)."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import tifffile

from lightsuite.config.models import SpinalCordPipelineConfig
from lightsuite.preprocess.cord_checkpoint import CordRegOptsCheckpoint
from lightsuite.registration.warp import warp_volume_affine

CORRESPONDING_POINTS_JSON = "corresponding_points.json"


@dataclass
class CordControlPointSession:
    histology_control_points: list[list[list[float]]]
    atlas_control_points: list[list[list[float]]]

    def save(self, path: Path) -> None:
        import json

        payload = {
            "histology_control_points": self.histology_control_points,
            "atlas_control_points": self.atlas_control_points,
        }
        path.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    @classmethod
    def load(cls, path: Path) -> CordControlPointSession:
        import json

        raw = json.loads(path.read_text(encoding="utf-8"))
        return cls(
            histology_control_points=raw.get("histology_control_points", []),
            atlas_control_points=raw.get("atlas_control_points", []),
        )


def default_session_path(save_path: Path) -> Path:
    return save_path / CORRESPONDING_POINTS_JSON


def load_cord_match_data(config: SpinalCordPipelineConfig) -> tuple[CordRegOptsCheckpoint, np.ndarray, np.ndarray, np.ndarray]:
    save_path = config.sample.save_path.expanduser()
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
    aff = np.asarray(checkpoint.affine_atlas_to_samp, dtype=float)
    tv_warp = warp_volume_affine(tv, aff, straightvol.shape, order=1, point_coords="array")
    av_warp = warp_volume_affine(av.astype(np.float32), aff, straightvol.shape, order=0, point_coords="array")
    return checkpoint, straightvol, tv_warp, av_warp


def run_spinal_match_points(
    config: SpinalCordPipelineConfig,
    *,
    headless: bool = False,
) -> Path:
    """Interactive control-point matching or empty session in headless mode."""
    save_path = config.sample.save_path.expanduser()
    out_path = default_session_path(save_path)
    n_slices = 100
    session = CordControlPointSession(
        histology_control_points=[[] for _ in range(n_slices)],
        atlas_control_points=[[] for _ in range(n_slices)],
    )
    if headless:
        session.save(out_path)
        return out_path

    try:
        import napari
    except ImportError as exc:
        msg = "Napari GUI requires: uv sync --extra gui"
        raise RuntimeError(msg) from exc

    _, straightvol, tv_warp, av_warp = load_cord_match_data(config)
    vol_disp = np.clip(straightvol / max(float(np.quantile(straightvol, 0.999)), 1.0) * 255, 0, 255).astype(np.uint8)
    atlas_disp = np.clip(tv_warp / max(float(np.quantile(tv_warp, 0.999)), 1.0) * 255, 0, 255).astype(np.uint8)

    viewer = napari.Viewer(title="Spinal cord match points")
    viewer.add_image(vol_disp, name="sample", colormap="gray")
    viewer.add_image(atlas_disp, name="atlas", colormap="magma", opacity=0.5)
    hist_layer = viewer.add_points(name="histology", face_color="cyan", size=10)
    atlas_layer = viewer.add_points(name="atlas", face_color="yellow", size=10)

    state = {"hist": [], "atlas": [], "pair_idx": 0}

    def save_session() -> None:
        session.histology_control_points = [state["hist"]]
        session.atlas_control_points = [state["atlas"]]
        session.save(out_path)

    @viewer.bind_key("s")
    def _save(_v) -> None:
        save_session()

    napari.run()
    if not out_path.is_file():
        save_session()
    return out_path
