"""Tests for spinal cord affine warping helpers."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from lightsuite.registration.cord_affine import (
    build_cord_z_transinit,
    warp_cord_atlas_to_straightvol,
    write_inverse_elastix_affine,
)
from lightsuite.registration.elastix.affine import parse_elastix_affine
from lightsuite.registration.warp import warp_volume_affine


def test_build_cord_z_transinit_scales_depth() -> None:
    transinit = build_cord_z_transinit(1000, 200)
    assert transinit[2, 2] == pytest.approx(1000 * 0.98 / 200)
    assert transinit[2, 3] == pytest.approx(500 - transinit[2, 2] * 100)


def test_warp_cord_z_only_matches_imwarp() -> None:
    atlas = np.zeros((20, 10, 5), dtype=np.uint16)
    atlas[8:12, 4:6, 2] = 7
    transinit = build_cord_z_transinit(40, 5)
    out_shape = (20, 10, 40)
    a = warp_cord_atlas_to_straightvol(
        atlas,
        transinit=transinit,
        elastix_affine_path=None,
        output_shape=out_shape,
        spacing_mm=0.02,
        work_dir=Path("."),
        nearest=True,
    )
    b = warp_volume_affine(atlas.astype(np.float32), transinit, out_shape, order=0, point_coords="array")
    assert np.array_equal(a, b.astype(np.uint16))


def test_write_inverse_elastix_affine_inverts_matrix_and_translation(tmp_path: Path) -> None:
    a_mat = np.array(
        [[1.05, -0.02, 0.39], [-0.03, 1.01, 0.13], [0.002, 0.003, 1.28]],
        dtype=float,
    )
    t_vec = np.array([-0.02, -0.10, 0.40])
    params = " ".join(f"{v:.16g}" for v in np.concatenate([a_mat.reshape(-1), t_vec]))
    forward = tmp_path / "fwd.txt"
    forward.write_text(
        '(Transform "AffineTransform")\n'
        "(NumberOfParameters 12)\n"
        f"(TransformParameters {params})\n"
        "(CenterOfRotationPoint 1.6 1.25 10.51)\n"
        "(Size 161 126 1052)\n",
        encoding="utf-8",
    )
    out = write_inverse_elastix_affine(forward, tmp_path / "inv.txt")

    import re

    inv_params = np.fromstring(
        re.search(r"\(\s*TransformParameters\s+([^)]+)\)", out.read_text()).group(1),
        sep=" ",
    )
    a_inv = inv_params[:9].reshape(3, 3)
    t_inv = inv_params[9:12]
    assert np.allclose(a_inv, np.linalg.inv(a_mat))
    assert np.allclose(t_inv, -np.linalg.inv(a_mat) @ t_vec)
    # grid metadata is preserved so transformix resamples onto the same physical extent
    assert "(Size 161 126 1052)" in out.read_text()


def test_parse_elastix_affine_has_proper_homogeneous_row() -> None:
    fixture = Path(__file__).resolve().parent / "fixtures" / "spinal_cord" / "parity"
    candidates = list(fixture.glob("**/affine_atlas_to_samp*.txt"))
    if not candidates:
        pytest.skip("no elastix affine fixture")
    matrix = parse_elastix_affine(candidates[0])
    assert np.allclose(matrix[3, :3], 0.0)
    assert matrix[3, 3] == pytest.approx(1.0)
