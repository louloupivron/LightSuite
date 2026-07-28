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
    BrainPipelineConfig,
)
from lightsuite.export.atlas_space import transform_mask_to_atlas
from lightsuite.export.brain_export import _load_transform_params
from lightsuite.import_.adapters import load_annotation, prepare_points_for_sample
from lightsuite.import_.models import AnnotationImportResult, ImportedMask, ImportedPoints
from lightsuite.import_.sample_reference import (
    load_sample_reference,
    validate_mask_against_reference,
)
from lightsuite.import_.transform import (
    sample_points_to_registration_voxels,
    transform_points_to_atlas,
)
from lightsuite.preprocess.checkpoint import RegOptsCheckpoint
from lightsuite.preprocess.slice_ops import output_xy_shape, output_z_count
from lightsuite.registration.points import volume_indices_to_cloud_xyz

console = Console()


def _slug(label: str) -> str:
    cleaned = "".join(ch if ch.isalnum() or ch in "-_" else "_" for ch in label.strip())
    return cleaned or "annotation"


def _registration_shape_native_yxz(
    shape_yxz: tuple[int, int, int],
    voxel_um: list[float],
    registres_um: float,
) -> tuple[int, int, int]:
    """Downsampled native grid (Y, X, Z) before orientation permute — matches preprocess TIFFs."""
    ny, nx, nz = shape_yxz
    vx, vy, vz = (float(v) for v in voxel_um)
    scale_xy = vx / float(registres_um)
    scale_z = vz / float(registres_um)
    out_h, out_w = output_xy_shape(ny, nx, scale_xy)
    out_z = output_z_count(nz, scale_z)
    return out_h, out_w, out_z


def _resample_mask_native_to_registration(
    mask: np.ndarray,
    *,
    target_shape_yxz: tuple[int, int, int],
    native_voxel_um: list[float],
    registres_um: float,
) -> np.ndarray:
    """Downsample a native-resolution mask onto the unpermuted registration grid."""
    vol = np.asarray(mask)
    reg_shape = _registration_shape_native_yxz(vol.shape, native_voxel_um, registres_um)
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
    *,
    transform_params,
    checkpoint: RegOptsCheckpoint,
    reference,
    output_dir: Path,
    temp_dir: Path,
    write_csv: bool,
) -> AnnotationImportResult:
    loaded = load_annotation(spec)
    if not isinstance(loaded, ImportedPoints):
        msg = f"Expected point annotation for {spec.path}"
        raise TypeError(msg)

    prepared = prepare_points_for_sample(loaded, reference=reference)
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

    reg_yxz = sample_points_to_registration_voxels(
        prepared.coordinates,
        transform_params,
        registres_um=checkpoint.registres_um,
        content_crop_start=checkpoint.content_crop_start,
    )
    reg_coords = volume_indices_to_cloud_xyz(reg_yxz)

    slug = _slug(prepared.label)
    npz_path = output_dir / f"{slug}_atlas_coords.npz"
    sample_npz_path = output_dir / f"{slug}_sample_coords.npz"
    np.savez_compressed(
        npz_path,
        atlasptcoords=atlas_coords.astype(np.float32),
        sampleptcoords=prepared.coordinates.astype(np.float32),
    )
    np.savez_compressed(
        sample_npz_path,
        regptcoords=reg_coords.astype(np.float32),
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
        sample_points_path=sample_npz_path,
        n_input=int(loaded.coordinates.shape[0]),
        n_atlas=int(atlas_coords.shape[0]),
        n_sample=int(reg_coords.shape[0]),
    )


def _import_mask(
    spec: AnnotationImportConfig,
    *,
    transform_params,
    checkpoint: RegOptsCheckpoint,
    reference,
    output_dir: Path,
    temp_dir: Path,
) -> AnnotationImportResult:
    loaded = load_annotation(spec)
    if not isinstance(loaded, ImportedMask):
        msg = f"Expected mask annotation for {spec.path}"
        raise TypeError(msg)

    voxel_um = [float(v) for v in checkpoint.voxel_um]
    validate_mask_against_reference(loaded.volume.shape, voxel_um, reference)
    loaded = ImportedMask(
        label=loaded.label,
        volume=loaded.volume,
        voxel_um=voxel_um,
        source_path=loaded.source_path,
        metadata=loaded.metadata,
    )

    reg_shape = _registration_shape_native_yxz(
        loaded.volume.shape,
        voxel_um,
        checkpoint.registres_um,
    )
    mask_reg = _resample_mask_native_to_registration(
        loaded.volume,
        target_shape_yxz=reg_shape,
        native_voxel_um=voxel_um,
        registres_um=checkpoint.registres_um,
    )
    permuted_shape = tuple(int(v) for v in transform_params.regvolsize)
    from lightsuite.registration.volume import permute_brain_volume

    mask_sample = permute_brain_volume(mask_reg, transform_params.permute_sample_to_atlas)
    if mask_sample.shape != permuted_shape:
        mask_sample = _resample_mask_native_to_registration(
            mask_sample,
            target_shape_yxz=permuted_shape,
            native_voxel_um=[checkpoint.registres_um, checkpoint.registres_um, checkpoint.registres_um],
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
    sample_mask_path = output_dir / f"{slug}_in_sample_20um.tif"
    save_registration_volume(mask_sample.astype(np.uint16), sample_mask_path)

    console.print(
        f"[green]{loaded.label}[/green]: mask warped to atlas "
        f"({int(np.count_nonzero(atlas_mask))} foreground voxels) "
        f"and sample grid ({int(np.count_nonzero(mask_sample))} voxels)"
    )
    return AnnotationImportResult(
        label=loaded.label,
        kind="mask",
        atlas_mask_path=out_path,
        sample_mask_path=sample_mask_path,
        n_input=int(np.count_nonzero(loaded.volume)),
        n_atlas=int(np.count_nonzero(atlas_mask)),
        n_sample=int(np.count_nonzero(mask_sample)),
    )


def run_brain_import_annotations(
    config: BrainPipelineConfig,
    *,
    annotations: list[AnnotationImportConfig] | None = None,
    write_csv: bool | None = None,
) -> list[AnnotationImportResult]:
    """Import native sample-space annotations into atlas space."""
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
    reference = load_sample_reference(save_path)
    transform_params = _load_transform_params(save_path)
    output_dir = save_path / "volume_registered"
    output_dir.mkdir(parents=True, exist_ok=True)
    temp_root = save_path / "import_annotations_temp"
    temp_root.mkdir(parents=True, exist_ok=True)

    console.print(
        f"Importing {len(specs)} annotation source(s) "
        f"(native grid {reference.shape_yxz}, voxel_um={reference.voxel_um})..."
    )
    t0 = time.perf_counter()
    results: list[AnnotationImportResult] = []
    for spec in specs:
        if spec.format == AnnotationFormat.MASK_TIFF:
            results.append(
                _import_mask(
                    spec,
                    transform_params=transform_params,
                    checkpoint=checkpoint,
                    reference=reference,
                    output_dir=output_dir,
                    temp_dir=temp_root,
                )
            )
        else:
            results.append(
                _import_points(
                    spec,
                    transform_params=transform_params,
                    checkpoint=checkpoint,
                    reference=reference,
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
                    "sample_points_path": str(r.sample_points_path) if r.sample_points_path else None,
                    "sample_mask_path": str(r.sample_mask_path) if r.sample_mask_path else None,
                    "n_sample": r.n_sample,
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
