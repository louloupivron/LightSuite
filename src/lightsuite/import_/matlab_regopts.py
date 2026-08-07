"""Import MATLAB ``regopts.mat`` into a Python ``regopts.json`` checkpoint."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
from scipy.io import loadmat

from lightsuite.preprocess.checkpoint import RegOptsCheckpoint


def _mat_struct_fields(obj: Any) -> list[str]:
    if hasattr(obj, "_fieldnames"):
        return list(obj._fieldnames)
    return []


def _get_field(obj: Any, name: str, default: Any = None) -> Any:
    if hasattr(obj, name):
        return getattr(obj, name)
    return default


def _as_list_int(values: Any) -> list[int]:
    arr = np.asarray(values, dtype=int).ravel()
    return [int(v) for v in arr.tolist()]


def _as_list_float(values: Any) -> list[float]:
    arr = np.asarray(values, dtype=float).ravel()
    return [float(v) for v in arr.tolist()]


def _load_numeric_export(save_path: Path) -> dict[str, np.ndarray]:
    """Load optional ``regopts_python_export.mat`` (see ``demos/export_regopts_for_python.m``)."""
    export_path = save_path / "regopts_python_export.mat"
    if not export_path.is_file():
        return {}
    raw = loadmat(export_path, squeeze_me=True, struct_as_record=False)
    out: dict[str, np.ndarray] = {}
    for key in ("original_trans_matrix", "tform_affine_samp20um_to_atlas_10um_px_matrix"):
        if key in raw:
            matrix = np.asarray(raw[key], dtype=float)
            if matrix.shape == (4, 4):
                out[key] = matrix
    return out


def load_matlab_regopts_struct(mat_path: Path) -> Any:
    mat_path = mat_path.expanduser()
    if not mat_path.is_file():
        msg = f"MATLAB regopts not found: {mat_path}"
        raise FileNotFoundError(msg)
    raw = loadmat(mat_path, squeeze_me=True, struct_as_record=False)
    if "opts" not in raw:
        msg = f"{mat_path.name} does not contain an 'opts' struct."
        raise ValueError(msg)
    return raw["opts"]


def matlab_regopts_to_checkpoint(
    opts: Any,
    *,
    save_path: Path,
    channel_paths: dict[int, str] | None = None,
    original_trans: np.ndarray | None = None,
) -> RegOptsCheckpoint:
    """Build :class:`RegOptsCheckpoint` from a loaded MATLAB ``opts`` struct."""
    save_path = save_path.expanduser()
    ch_primary = int(_get_field(opts, "channelforregister", 1))
    ch_secondary = _get_field(opts, "channelforregister_secondary", None)
    if ch_secondary is not None:
        ch_secondary = int(ch_secondary)
        if ch_secondary <= 0:
            ch_secondary = None

    regvolpaths: dict[str, str] = {}
    if channel_paths:
        regvolpaths = {str(k): v for k, v in channel_paths.items()}
    else:
        for ich in range(1, int(_get_field(opts, "Nchans", 1)) + 1):
            candidate = save_path / f"chan_{ich}_sample_register_{int(_get_field(opts, 'registres', 20))}um.tif"
            if candidate.is_file():
                regvolpaths[str(ich)] = str(candidate)

    regvolpath = regvolpaths.get(str(ch_primary), "")
    if not regvolpath:
        legacy = _get_field(opts, "regvolpath", "")
        if legacy:
            name = Path(str(legacy)).name
            regvolpath = str(save_path / name)

    regvolpath_secondary = None
    if ch_secondary is not None:
        regvolpath_secondary = regvolpaths.get(str(ch_secondary))
        if regvolpath_secondary is None:
            legacy2 = _get_field(opts, "regvolpath_secondary", None)
            if legacy2:
                regvolpath_secondary = str(save_path / Path(str(legacy2)).name)

    autocpsample = np.asarray(_get_field(opts, "autocpsample", []), dtype=float)
    autocpatlas = np.asarray(_get_field(opts, "autocpatlas", []), dtype=float)
    if autocpsample.ndim == 1 and autocpsample.size == 3:
        autocpsample = autocpsample.reshape(1, 3)
    if autocpatlas.ndim == 1 and autocpatlas.size == 3:
        autocpatlas = autocpatlas.reshape(1, 3)

    checkpoint = RegOptsCheckpoint(
        sample_name=str(_get_field(opts, "mousename", "sample")),
        ny=int(_get_field(opts, "Ny", 0)),
        nx=int(_get_field(opts, "Nx", 0)),
        nz=int(_get_field(opts, "Nz", 0)),
        nchans=int(_get_field(opts, "Nchans", 1)),
        voxel_um=_as_list_float(_get_field(opts, "pxsize", [1.0, 1.0, 1.0])),
        registres_um=float(_get_field(opts, "registres", 20)),
        regvolpath=regvolpath,
        regvolpath_secondary=regvolpath_secondary,
        regvolpaths=regvolpaths,
        tiff_type=str(_get_field(opts, "tifftype", "channelperfile")),
        channel_primary=ch_primary,
        channel_secondary=ch_secondary,
        permute_sample_to_atlas=_as_list_int(_get_field(opts, "permute_sample_to_atlas", [1, 2, 3])),
        downfac_reg=float(_get_field(opts, "downfac_reg", 1.0)),
        autocpsample=autocpsample.tolist() if autocpsample.size else [],
        autocpatlas=autocpatlas.tolist() if autocpatlas.size else [],
        brain_atlas=str(_get_field(opts, "brain_atlas", "allen")).lower(),
    )
    if original_trans is not None:
        checkpoint.original_trans = np.asarray(original_trans, dtype=float).tolist()
    return checkpoint


def import_matlab_regopts(
    matlab_save_path: Path,
    *,
    python_save_path: Path | None = None,
    original_trans: np.ndarray | None = None,
) -> RegOptsCheckpoint:
    """Write ``regopts.json`` under ``python_save_path`` from MATLAB ``regopts.mat``."""
    matlab_save_path = matlab_save_path.expanduser()
    python_save_path = (python_save_path or matlab_save_path).expanduser()
    python_save_path.mkdir(parents=True, exist_ok=True)

    opts = load_matlab_regopts_struct(matlab_save_path / "regopts.mat")
    numeric = _load_numeric_export(matlab_save_path)
    if "original_trans_matrix" in numeric:
        matrix = np.asarray(numeric["original_trans_matrix"], dtype=float)
        # Legacy exports stored MATLAB ``tform.T`` (translation in row 4). New exports
        # store the transpose (translation in column 4) for Python ``transform_points``.
        if abs(matrix[3, 3] - 1.0) < 1e-6 and np.linalg.norm(matrix[:3, 3]) < 1e-6:
            if np.linalg.norm(matrix[3, :3]) > 1e-6:
                matrix = matrix.T
        original_trans = matrix

    checkpoint = matlab_regopts_to_checkpoint(
        opts,
        save_path=python_save_path,
        original_trans=original_trans,
    )
    checkpoint.save(python_save_path / "regopts.json")
    return checkpoint
