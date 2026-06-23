"""Napari threshold inspector for unassigned registration QC."""

from __future__ import annotations

import numpy as np

from lightsuite.gui.inspect_brain_imports import _contrast_limits, volume_yxz_to_napari_zyx


def open_registration_qc_inspector(
    volume: np.ndarray,
    division_labels: np.ndarray,
    *,
    threshold: float,
) -> None:
    """Open Napari showing above-threshold signal in unassigned (id 0) divisions."""
    try:
        import napari
    except ImportError as exc:
        msg = "Registration QC inspect requires: uv sync --extra gui"
        raise RuntimeError(msg) from exc

    vol = np.asarray(volume, dtype=np.float32)
    labels = np.asarray(division_labels, dtype=np.int32)
    lims = _contrast_limits(vol)

    viewer = napari.Viewer(title="Registration QC — threshold inspection")
    viewer.add_image(
        volume_yxz_to_napari_zyx(vol),
        name="registered volume",
        colormap="gray",
        contrast_limits=lims,
        blending="translucent",
        depiction="volume",
    )
    viewer.add_labels(
        volume_yxz_to_napari_zyx(labels),
        name="division labels",
        opacity=0.25,
        visible=False,
    )

    unassigned = labels == 0
    outside_signal = np.where(unassigned & (vol >= float(threshold)), vol, np.nan).astype(
        np.float32
    )
    viewer.add_image(
        volume_yxz_to_napari_zyx(outside_signal),
        name=f"id0 signal >= {threshold:g}",
        colormap="magenta",
        blending="additive",
        contrast_limits=lims,
        visible=True,
    )
    napari.run()
