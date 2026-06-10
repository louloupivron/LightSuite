"""Orchestrate external annotation import after brain registration."""

from __future__ import annotations

import json
import shutil
import time
from pathlib import Path

import numpy as np
import pandas as pd
from rich.console import Console
from skimage.transform import resize as sk_resize

from lightsuite.config.models import (
    AnnotationFormat,
    AnnotationImportConfig,
    AnnotationRole,
    BrainPipelineConfig,
)
from lightsuite.export.atlas_space import transform_mask_to_atlas
from lightsuite.export.brain_export import _load_transform_params
from lightsuite.import_.adapters import load_annotation, prepare_points_for_sample
from lightsuite.import_.models import AnnotationImportResult, ImportedMask, ImportedPoints
from lightsuite.import_.transform import transform_points_to_atlas
from lightsuite.preprocess.checkpoint import RegOptsCheckpoint

console = Console()


def _mask_voxel_um(
    mask: ImportedMask,
    spec: AnnotationImportConfig,
    checkpoint: RegOptsCheckpoint,
) -> ImportedMask:
    """Resolve mask voxel size from import spec or the preprocessed sample."""
    if spec.voxel_um is not None:
        voxel_um = [float(v) for v in spec.voxel_um]
    elif mask.metadata.get("voxel_um_from_sample"):
        voxel_um = [float(v) for v in checkpoint.voxel_um]
    else:
        voxel_um = list(mask.voxel_um)
    return ImportedMask(
        label=mask.label,
        volume=mask.volume,
        voxel_um=voxel_um,
        source_path=mask.source_path,
        metadata=mask.metadata,
    )


def _slug(label: str) -> str:
    cleaned = "".join(ch if ch.isalnum() or ch in "-_" else "_" for ch in label.strip())
    return cleaned or "annotation"


def _resample_mask_to_registration(
    mask: ImportedMask,
    *,
    target_shape_yxz: tuple[int, int, int],
    native_voxel_um: list[float],
    registres_um: float,
) -> np.ndarray:
    """Resample an external mask onto the registration-resolution Y,X,Z grid."""
    vol = np.asarray(mask.volume)
    src_um = np.asarray(mask.voxel_um, dtype=np.float64)
    native_um = np.asarray(native_voxel_um, dtype=np.float64)
    reg_um = np.full(3, float(registres_um), dtype=np.float64)

    native_shape = tuple(
        max(1, int(round(vol.shape[i] * src_um[i] / native_um[i]))) for i in range(3)
    )
    if native_shape != vol.shape:
        vol = sk_resize(
            vol,
            native_shape,
            order=0,
            preserve_range=True,
            anti_aliasing=False,
        )

    reg_shape = tuple(
        max(1, int(round(native_shape[i] * native_um[i] / reg_um[i]))) for i in range(3)
    )
    if reg_shape != vol.shape:
        vol = sk_resize(
            vol,
            reg_shape,
            order=0,
            preserve_range=True,
            anti_aliasing=False,
        )

    if reg_shape != target_shape_yxz:
        vol = sk_resize(
            vol,
            target_shape_yxz,
            order=0,
            preserve_range=True,
            anti_aliasing=False,
        )
    return (vol > 0).astype(np.uint8)


def _import_points(
    spec: AnnotationImportConfig,
    config: BrainPipelineConfig,
    *,
    transform_params,
    checkpoint: RegOptsCheckpoint,
    output_dir: Path,
    temp_dir: Path,
    write_csv: bool,
) -> AnnotationImportResult:
    loaded = load_annotation(spec)
    if not isinstance(loaded, ImportedPoints):
        msg = f"Expected point annotation for {spec.path}"
        raise TypeError(msg)

    target_size = (checkpoint.ny, checkpoint.nx, checkpoint.nz)
    prepared = prepare_points_for_sample(
        loaded,
        spec,
        target_voxel_um=checkpoint.voxel_um,
        target_size_yxz=target_size,
    )
    if prepared.coordinates.size == 0:
        console.print(f"[yellow]No in-bounds points for {prepared.label}[/yellow]")
        return AnnotationImportResult(
            label=prepared.label,
            kind="points",
            n_input=int(loaded.coordinates.shape[0]),
            n_atlas=0,
        )

    atlas_coords = transform_points_to_atlas(
        prepared.coordinates,
        transform_params,
        registres_um=checkpoint.registres_um,
        temp_dir=temp_dir / _slug(prepared.label),
    )

    slug = _slug(prepared.label)
    npz_path = output_dir / f"{slug}_atlas_coords.npz"
    np.savez_compressed(
        npz_path,
        atlasptcoords=atlas_coords.astype(np.float32),
        sampleptcoords=prepared.coordinates.astype(np.float32),
    )

    csv_path: Path | None = None
    if write_csv:
        csv_path = output_dir / f"{slug}_atlas_coords.csv"
        frame = pd.DataFrame(
            atlas_coords[:, :3],
            columns=["atlas_x", "atlas_y", "atlas_z"],
        )
        if atlas_coords.shape[1] > 3:
            for idx in range(3, atlas_coords.shape[1]):
                frame[f"feature_{idx - 2}"] = atlas_coords[:, idx]
        frame.to_csv(csv_path, index=False)

    console.print(
        f"[green]{prepared.label}[/green]: {prepared.coordinates.shape[0]} sample points "
        f"→ {atlas_coords.shape[0]} atlas points"
    )
    return AnnotationImportResult(
        label=prepared.label,
        kind="points",
        atlas_points_path=npz_path,
        atlas_csv_path=csv_path,
        n_input=int(loaded.coordinates.shape[0]),
        n_atlas=int(atlas_coords.shape[0]),
    )


def _import_mask(
    spec: AnnotationImportConfig,
    config: BrainPipelineConfig,
    *,
    transform_params,
    checkpoint: RegOptsCheckpoint,
    output_dir: Path,
    temp_dir: Path,
) -> AnnotationImportResult:
    loaded = load_annotation(spec)
    if not isinstance(loaded, ImportedMask):
        msg = f"Expected mask annotation for {spec.path}"
        raise TypeError(msg)
    loaded = _mask_voxel_um(loaded, spec, checkpoint)

    reg_shape = tuple(int(v) for v in transform_params.regvolsize)
    mask_reg = _resample_mask_to_registration(
        loaded,
        target_shape_yxz=reg_shape,
        native_voxel_um=checkpoint.voxel_um,
        registres_um=checkpoint.registres_um,
    )
    atlas_mask = transform_mask_to_atlas(
        mask_reg,
        transform_params,
        permute=transform_params.permute_sample_to_atlas,
        spacing_mm=checkpoint.registres_um * 1e-3,
        temp_dir=temp_dir / _slug(loaded.label),
    )

    slug = _slug(loaded.label)
    out_path = output_dir / f"{slug}_registered_atlas.tif"
    from lightsuite.io.tiff_write import save_registration_volume

    save_registration_volume(atlas_mask, out_path)

    console.print(
        f"[green]{loaded.label}[/green]: mask warped to atlas "
        f"({int(np.count_nonzero(atlas_mask))} foreground voxels)"
    )
    return AnnotationImportResult(
        label=loaded.label,
        kind="mask",
        atlas_mask_path=out_path,
        n_input=int(np.count_nonzero(loaded.volume)),
        n_atlas=int(np.count_nonzero(atlas_mask)),
    )


def run_brain_import_annotations(
    config: BrainPipelineConfig,
    *,
    annotations: list[AnnotationImportConfig] | None = None,
    write_csv: bool | None = None,
) -> list[AnnotationImportResult]:
    """Import external LCT / Arivis annotations into atlas space."""
    if shutil.which("transformix") is None:
        msg = "transformix must be on PATH for annotation import."
        raise RuntimeError(msg)

    specs = annotations
    if specs is None:
        if config.import_config is None or not config.import_config.annotations:
            msg = (
                "No import.annotations configured. Add an 'import' section to the YAML "
                "or pass annotations explicitly."
            )
            raise ValueError(msg)
        specs = config.import_config.annotations

    if write_csv is None:
        write_csv = config.import_config.write_csv if config.import_config is not None else True

    save_path = config.sample.save_path.expanduser()
    regopts_path = save_path / "regopts.json"
    if not regopts_path.is_file():
        msg = f"Missing {regopts_path}. Run preprocess and register first."
        raise FileNotFoundError(msg)

    checkpoint = RegOptsCheckpoint.load(regopts_path)
    transform_params = _load_transform_params(save_path)
    output_dir = save_path / "volume_registered"
    output_dir.mkdir(parents=True, exist_ok=True)
    temp_root = save_path / "import_annotations_temp"
    temp_root.mkdir(parents=True, exist_ok=True)

    console.print(f"Importing {len(specs)} annotation source(s)...")
    t0 = time.perf_counter()
    results: list[AnnotationImportResult] = []
    for spec in specs:
        if (
            spec.role == AnnotationRole.MASK
            or spec.format in (AnnotationFormat.LCT_ZARR, AnnotationFormat.TIFF_MASK)
        ):
            results.append(
                _import_mask(
                    spec,
                    config,
                    transform_params=transform_params,
                    checkpoint=checkpoint,
                    output_dir=output_dir,
                    temp_dir=temp_root,
                )
            )
        else:
            results.append(
                _import_points(
                    spec,
                    config,
                    transform_params=transform_params,
                    checkpoint=checkpoint,
                    output_dir=output_dir,
                    temp_dir=temp_root,
                    write_csv=write_csv,
                )
            )

    summary_path = output_dir / "import_annotations_summary.json"
    summary_path.write_text(
        json.dumps(
            [
                {
                    "label": r.label,
                    "kind": r.kind,
                    "n_input": r.n_input,
                    "n_atlas": r.n_atlas,
                    "atlas_points_path": str(r.atlas_points_path) if r.atlas_points_path else None,
                    "atlas_mask_path": str(r.atlas_mask_path) if r.atlas_mask_path else None,
                    "atlas_csv_path": str(r.atlas_csv_path) if r.atlas_csv_path else None,
                }
                for r in results
            ],
            indent=2,
        ),
        encoding="utf-8",
    )
    console.print(
        f"[green]Import complete[/green] in {time.perf_counter() - t0:.1f}s. "
        f"Summary: {summary_path}"
    )
    return results
