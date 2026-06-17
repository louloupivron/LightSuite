"""Tests for mesoSPIM TIFF I/O helpers."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import tifffile

from lightsuite.mesospim.config_models import MesospimGeometryConfig, MesospimTiffRemapConfig
from lightsuite.mesospim.io import (
    _remapped_z_to_raw_page,
    read_tiff_xy_slice_at_physical_um,
    read_tiff_xy_slice_at_z_index,
    tiff_shape,
)
from lightsuite.mesospim.meta import meta_path_for_tiff


def _write_stack(path: Path, shape_zyx: tuple[int, int, int]) -> None:
    z, y, x = shape_zyx
    arr = np.arange(z * y * x, dtype=np.uint16).reshape(shape_zyx)
    tifffile.imwrite(path, arr, imagej=True)


def _write_meta(path: Path) -> None:
    path.write_text(
        "\n".join(
            [
                "[Pixelsize in um] 5.0",
                "[x_pos] 0.0",
                "[y_pos] 0.0",
                "[z_start] 0.0",
                "[z_end] 8.0",
                "[z_stepsize] 2.0",
                "[z_planes] 5",
                "[x_pixels] 16",
                "[y_pixels] 16",
            ]
        ),
        encoding="utf-8",
    )


def test_remapped_z_to_raw_page() -> None:
    assert _remapped_z_to_raw_page(2, 5, reverse_z=False) == 2
    assert _remapped_z_to_raw_page(2, 5, reverse_z=True) == 2


def test_read_tiff_xy_slice_at_physical_um(tmp_path: Path) -> None:
    overview = tmp_path / "overview.tif"
    roi = tmp_path / "roi.tif"
    _write_stack(overview, (5, 16, 16))
    _write_stack(roi, (5, 16, 16))
    _write_meta(meta_path_for_tiff(overview))
    _write_meta(meta_path_for_tiff(roi))

    geometry = MesospimGeometryConfig()
    remap = MesospimTiffRemapConfig()
    meta = {
        "Pixelsize in um": 5.0,
        "x_pos": 0.0,
        "y_pos": 0.0,
        "z_start": 0.0,
        "z_end": 8.0,
        "z_stepsize": 2.0,
    }

    sl = read_tiff_xy_slice_at_physical_um(
        overview,
        meta=meta,
        geometry=geometry,
        overview_path=overview,
        roi_path=roi,
        remap=remap,
        cx_um=40.0,
        cy_um=40.0,
        cz_um=4.0,
    )
    assert sl.shape == tiff_shape(overview)[1:]
    assert sl.dtype == np.float32


def test_read_tiff_xy_slice_at_z_index(tmp_path: Path) -> None:
    overview = tmp_path / "overview.tif"
    roi = tmp_path / "roi.tif"
    shape = (5, 16, 16)
    _write_stack(overview, shape)
    _write_stack(roi, shape)

    remap = MesospimTiffRemapConfig()
    sl = read_tiff_xy_slice_at_z_index(
        overview,
        2,
        overview_path=overview,
        roi_path=roi,
        remap=remap,
    )
    assert sl.shape == shape[1:]
    assert sl.dtype == np.float32
    assert sl[0, 0] == float(2 * 16 * 16)


def test_read_raw_plane_yx_volumetric_memmap(tmp_path: Path, monkeypatch) -> None:
    path = tmp_path / "vol.tif"
    path.touch()
    vol = np.arange(3 * 4 * 4, dtype=np.uint16).reshape(3, 4, 4)

    class FakePage:
        ndim = 3

    class FakeTiffFile:
        pages = [FakePage()]

        def __init__(self, *_args, **_kwargs) -> None:
            pass

        def __enter__(self):
            return self

        def __exit__(self, *_args) -> None:
            return None

    monkeypatch.setattr(tifffile, "TiffFile", FakeTiffFile)
    monkeypatch.setattr(tifffile, "memmap", lambda _key: vol)

    from lightsuite.mesospim.io import _MEMMAP_CACHE, _read_raw_plane_yx

    _MEMMAP_CACHE.clear()
    plane = _read_raw_plane_yx(path, 1)
    assert plane.shape == (4, 4)
    assert plane[0, 0] == vol[1, 0, 0]


def test_parse_mesospim_meta(tmp_path: Path) -> None:
    meta_path = tmp_path / "stack.tif_meta.txt"
    meta_path.write_text(
        "\n".join(
            [
                "[Pixelsize in um] 6.55",
                "[x_pos] 17167.48",
                "[y_pos] 63592.24",
                "[z_start] 24015.62",
                "[z_end] 17215.53",
                "[z_stepsize] 5.0",
                "[z_planes] 1361",
                "[x_pixels] 2048",
            ]
        ),
        encoding="utf-8",
    )

    from lightsuite.mesospim.meta import parse_mesospim_meta

    meta = parse_mesospim_meta(meta_path)
    assert meta["Pixelsize in um"] == 6.55
    assert meta["x_pos"] == 17167.48
    assert meta["z_planes"] == 1361


def test_meta_path_for_tiff() -> None:
    path = Path("/data/1X/1-561-1x.tif")
    assert meta_path_for_tiff(path) == Path("/data/1X/1-561-1x.tif_meta.txt")
