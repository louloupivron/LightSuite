"""Orchestrate external annotation import after spinal cord registration."""

from __future__ import annotations

import json
import shutil
import time
from pathlib import Path

import numpy as np
import pandas as pd
from rich.console import Console

from lightsuite.atlas.fiederling import load_fiederling_atlas_volumes
from lightsuite.config.models import (
    AnnotationFormat,
    AnnotationImportConfig,
    SpinalCordPipelineConfig,
)
from lightsuite.import_.adapters import load_annotation, prepare_points_for_sample
from lightsuite.import_.cord_transform import filter_points_for_cord_registration, transform_points_to_cord_atlas
from lightsuite.import_.models import AnnotationImportResult, ImportedPoints
from lightsuite.import_.sample_reference import load_sample_reference
from lightsuite.preprocess.cord_checkpoint import CordRegOptsCheckpoint, CordTransformParamsCheckpoint
from lightsuite.registration.cord_longitudinal import load_longitudinal_correspondence, resolve_cord_z_transinit
from lightsuite.registration.cord_paths import cord_affine_transform_path, cord_save_path
from lightsuite.registration.straightening import load_slicetforms

console = Console()


def _slug(label: str) -> str:
    cleaned = "".join(ch if ch.isalnum() or ch in "-_" else "_" for ch in label.strip())
    return cleaned or "annotation"


def _load_transform_params(save_path: Path) -> CordTransformParamsCheckpoint:
    json_path = save_path / "transform_params.json"
    if not json_path.is_file():
        msg = f"Missing {json_path}. Run 'lightsuite spinal register' first."
        raise FileNotFoundError(msg)
    return CordTransformParamsCheckpoint.load(json_path)


def _import_points(
    spec: AnnotationImportConfig,
    *,
    transform_params: CordTransformParamsCheckpoint,
    reference,
    output_dir: Path,
    temp_dir: Path,
    elastix_affine_path: Path,
    transinit: np.ndarray,
    tforms: list[np.ndarray],
    spacing_mm: float,
    template_native_shape: tuple[int, int, int],
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

    reg_coords, _ = filter_points_for_cord_registration(
        prepared.coordinates,
        transform_params=transform_params,
    )
    if reg_coords.size == 0:
        console.print(
            f"[yellow]No points inside registration crop for {prepared.label}[/yellow]"
        )
        return AnnotationImportResult(
            label=prepared.label,
            kind="points",
            n_input=int(loaded.coordinates.shape[0]),
            n_atlas=0,
        )

    atlas_coords = transform_points_to_cord_atlas(
        reg_coords,
        transform_params=transform_params,
        elastix_affine_path=elastix_affine_path,
        transinit=transinit,
        tforms=tforms,
        spacing_mm=spacing_mm,
        template_native_shape=template_native_shape,
        temp_dir=temp_dir / _slug(prepared.label),
    )

    slug = _slug(prepared.label)
    npz_path = output_dir / f"{slug}_atlas_coords.npz"
    np.savez_compressed(
        npz_path,
        atlasptcoords=atlas_coords.astype(np.float32),
        sampleptcoords=reg_coords.astype(np.float32),
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
        f"[green]{prepared.label}[/green]: {reg_coords.shape[0]} sample points "
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


def run_cord_import_annotations(
    config: SpinalCordPipelineConfig,
    *,
    annotations: list[AnnotationImportConfig] | None = None,
    write_csv: bool | None = None,
) -> list[AnnotationImportResult]:
    """Import native sample-space annotations into Fiederling atlas space."""
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

    save_path = cord_save_path(config)
    regopts_path = save_path / "regopts.json"
    if not regopts_path.is_file():
        msg = f"Missing {regopts_path}. Run preprocess and register first."
        raise FileNotFoundError(msg)

    reference = load_sample_reference(save_path)
    transform_params = _load_transform_params(save_path)
    if not transform_params.slicetforms_path:
        msg = "Missing slicetforms_path in transform_params.json. Run init-registration first."
        raise RuntimeError(msg)

    tforms = load_slicetforms(transform_params.slicetforms_path)
    elastix_affine_path = cord_affine_transform_path(config)
    if not elastix_affine_path.is_file():
        msg = f"Missing elastix affine transform: {elastix_affine_path}"
        raise FileNotFoundError(msg)

    regopts = CordRegOptsCheckpoint.load(regopts_path)
    nslices = regopts.ikeeprange[1] - regopts.ikeeprange[0] + 1
    correspondence = load_longitudinal_correspondence(save_path)
    transinit = resolve_cord_z_transinit(
        nslices,
        int(transform_params.atlassize[2]),
        correspondence,
    )
    spacing_mm = float(config.registration.resolution_um) * 1e-3

    atlas_volumes = load_fiederling_atlas_volumes(config.atlas)
    template_native_shape = tuple(int(v) for v in atlas_volumes.template.shape)

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
            msg = (
                "Mask import is not yet supported for the spinal cord pipeline. "
                "Convert segmentation to points_csv or use the brain pipeline for masks."
            )
            raise NotImplementedError(msg)
        results.append(
            _import_points(
                spec,
                transform_params=transform_params,
                reference=reference,
                output_dir=output_dir,
                temp_dir=temp_root,
                elastix_affine_path=elastix_affine_path,
                transinit=transinit,
                tforms=tforms,
                spacing_mm=spacing_mm,
                template_native_shape=template_native_shape,
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
