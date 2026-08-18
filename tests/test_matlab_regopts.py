"""Tests for MATLAB regopts.mat import."""

from __future__ import annotations

from pathlib import Path

import numpy as np
from scipy.io import savemat

from lightsuite.import_.matlab_regopts import import_matlab_regopts


def _write_regopts_mat(path: Path) -> None:
    opts = np.dtype(
        [
            ("mousename", "O"),
            ("pxsize", "O"),
            ("atlasres", "O"),
            ("registres", "O"),
            ("brain_atlas", "O"),
            ("channelforregister", "O"),
            ("Ny", "O"),
            ("Nx", "O"),
            ("Nz", "O"),
            ("Nchans", "O"),
            ("tifftype", "O"),
            ("permute_sample_to_atlas", "O"),
            ("downfac_reg", "O"),
            ("autocpsample", "O"),
            ("autocpatlas", "O"),
        ]
    )
    row = np.zeros(1, dtype=opts)
    row[0]["mousename"] = "test"
    row[0]["pxsize"] = np.array([8.0, 8.0, 5.0])
    row[0]["atlasres"] = 20
    row[0]["registres"] = 20
    row[0]["brain_atlas"] = "perens"
    row[0]["channelforregister"] = 1
    row[0]["Ny"] = 100
    row[0]["Nx"] = 80
    row[0]["Nz"] = 60
    row[0]["Nchans"] = 1
    row[0]["tifftype"] = "channelperfile"
    row[0]["permute_sample_to_atlas"] = np.array([2, 1, 3])
    row[0]["downfac_reg"] = 1.0
    row[0]["autocpsample"] = np.array([[1.0, 2.0, 3.0]])
    row[0]["autocpatlas"] = np.array([[4.0, 5.0, 6.0]])
    savemat(path, {"opts": row}, do_compression=False)


def test_import_matlab_regopts(tmp_path: Path) -> None:
    mat_path = tmp_path / "regopts.mat"
    _write_regopts_mat(mat_path)
    (tmp_path / "chan_1_sample_register_20um.tif").write_bytes(b"")
    checkpoint = import_matlab_regopts(tmp_path)
    assert checkpoint.sample_name == "test"
    assert checkpoint.permute_sample_to_atlas == [2, 1, 3]
    assert len(checkpoint.autocpsample or []) == 1
