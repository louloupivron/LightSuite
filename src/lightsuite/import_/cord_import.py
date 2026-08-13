"""Orchestrate external annotation import after spinal cord registration."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
from rich.console import Console

from lightsuite.atlas.fiederling import load_fiederling_atlas_volumes
from lightsuite.config.models import (
    AnnotationImportConfig,
    SpinalCordPipelineConfig,
)
from lightsuite.import_.cord_transform import (
    filter_points_for_cord_registration,
    sample_points_to_straightened_grid,
    transform_points_to_cord_atlas,
)
from lightsuite.import_.models import AnnotationImportResult, ImportedMask, ImportedPoints
from lightsuite.import_.orchestrator import (
    require_transformix,
    resolve_annotation_specs,
    resolve_write_csv,
    run_annotation_import,
    slug_for_label,
)
from lightsuite.import_.sample_reference import SampleReference, load_sample_reference
from lightsuite.preprocess.cord_checkpoint import (
    CordRegOptsCheckpoint,
    CordTransformParamsCheckpoint,
)
from lightsuite.registration.cord_longitudinal import (
    load_longitudinal_correspondence,
    resolve_cord_z_transinit,
)
from lightsuite.registration.cord_paths import cord_affine_transform_path, cord_save_path
from lightsuite.registration.points import volume_indices_to_cloud_xyz
from lightsuite.registration.straightening import load_slicetforms

console = Console()

MASK_UNSUPPORTED = (
    "Mask import is not yet supported for the spinal cord pipeline. "
    "Convert segmentation to points_csv or use the brain pipeline for masks."
)


def _slug(label: str) -> str:
    return slug_for_label(label)


def _load_transform_params(save_path: Path) -> CordTransformParamsCheckpoint:
    json_path = save_path / "transform_params.json"
    if not json_path.is_file():
        msg = f"Missing {json_path}. Run 'lightsuite spinal register' first."
        raise FileNotFoundError(msg)
    return CordTransformParamsCheckpoint.load(json_path)


@dataclass
class CordAnnotationImporter:
    """Warp native sample points into Fiederling cord atlas space."""

    transform_params: CordTransformParamsCheckpoint
    sample_reference: SampleReference
    output_dir: Path
    temp_dir: Path
    elastix_affine_path: Path
    transinit: np.ndarray
    tforms: list[np.ndarray] = field(default_factory=list)
    spacing_mm: float = 0.02
    template_native_shape: tuple[int, int, int] = (0, 0, 0)
    write_csv: bool = True

    @property
    def reference(self) -> SampleReference:
        return self.sample_reference

    def import_points(self, points: ImportedPoints, *, slug: str) -> AnnotationImportResult:
        reg_coords, _ = filter_points_for_cord_registration(
            points.coordinates,
            transform_params=self.transform_params,
        )
        if reg_coords.size == 0:
            console.print(
                f"[yellow]No points inside registration crop for {points.label}[/yellow]"
            )
            return AnnotationImportResult(
                label=points.label,
                kind="points",
                n_input=int(points.coordinates.shape[0]),
                n_atlas=0,
            )

        atlas_coords = transform_points_to_cord_atlas(
            reg_coords,
            transform_params=self.transform_params,
            elastix_affine_path=self.elastix_affine_path,
            transinit=self.transinit,
            tforms=self.tforms,
            spacing_mm=self.spacing_mm,
            template_native_shape=self.template_native_shape,
            temp_dir=self.temp_dir / slug,
        )

        straight_yxz = sample_points_to_straightened_grid(
            reg_coords,
            transform_params=self.transform_params,
            tforms=self.tforms,
        )
        sample_reg_coords = volume_indices_to_cloud_xyz(straight_yxz)

        npz_path = self.output_dir / f"{slug}_atlas_coords.npz"
        sample_npz_path = self.output_dir / f"{slug}_sample_coords.npz"
        np.savez_compressed(
            npz_path,
            atlasptcoords=atlas_coords.astype(np.float32),
            sampleptcoords=reg_coords.astype(np.float32),
        )
        np.savez_compressed(
            sample_npz_path,
            regptcoords=sample_reg_coords.astype(np.float32),
            sampleptcoords=reg_coords.astype(np.float32),
        )

        csv_path: Path | None = None
        if self.write_csv:
            csv_path = self.output_dir / f"{slug}_atlas_coords.csv"
            frame = pd.DataFrame(
                atlas_coords[:, :3],
                columns=["atlas_x", "atlas_y", "atlas_z"],
            )
            if atlas_coords.shape[1] > 3:
                for idx in range(3, atlas_coords.shape[1]):
                    frame[f"feature_{idx - 2}"] = atlas_coords[:, idx]
            frame.to_csv(csv_path, index=False)

        console.print(
            f"[green]{points.label}[/green]: {reg_coords.shape[0]} sample points "
            f"→ {atlas_coords.shape[0]} atlas points"
        )
        return AnnotationImportResult(
            label=points.label,
            kind="points",
            atlas_points_path=npz_path,
            atlas_csv_path=csv_path,
            sample_points_path=sample_npz_path,
            n_input=int(points.coordinates.shape[0]),
            n_atlas=int(atlas_coords.shape[0]),
            n_sample=int(sample_reg_coords.shape[0]),
        )

    def import_mask(self, mask: ImportedMask, *, slug: str) -> AnnotationImportResult:
        raise NotImplementedError(MASK_UNSUPPORTED)


def run_cord_import_annotations(
    config: SpinalCordPipelineConfig,
    *,
    annotations: list[AnnotationImportConfig] | None = None,
    write_csv: bool | None = None,
) -> list[AnnotationImportResult]:
    """Import native sample-space annotations into Fiederling atlas space."""
    require_transformix()
    specs = resolve_annotation_specs(config.import_config, annotations)
    write_csv_resolved = resolve_write_csv(config.import_config, write_csv)

    save_path = cord_save_path(config)
    regopts_path = save_path / "regopts.json"
    if not regopts_path.is_file():
        msg = f"Missing {regopts_path}. Run preprocess and register first."
        raise FileNotFoundError(msg)

    transform_params = _load_transform_params(save_path)
    if not transform_params.slicetforms_path:
        msg = "Missing slicetforms_path in transform_params.json. Run init-registration first."
        raise RuntimeError(msg)

    elastix_affine_path = cord_affine_transform_path(config)
    if not elastix_affine_path.is_file():
        msg = f"Missing elastix affine transform: {elastix_affine_path}"
        raise FileNotFoundError(msg)

    regopts = CordRegOptsCheckpoint.load(regopts_path)
    nslices = regopts.ikeeprange[1] - regopts.ikeeprange[0] + 1
    transinit = resolve_cord_z_transinit(
        nslices,
        int(transform_params.atlassize[2]),
        load_longitudinal_correspondence(save_path),
    )
    atlas_volumes = load_fiederling_atlas_volumes(config.atlas)

    output_dir = save_path / "volume_registered"
    temp_root = save_path / "import_annotations_temp"
    temp_root.mkdir(parents=True, exist_ok=True)

    importer = CordAnnotationImporter(
        transform_params=transform_params,
        sample_reference=load_sample_reference(save_path),
        output_dir=output_dir,
        temp_dir=temp_root,
        elastix_affine_path=elastix_affine_path,
        transinit=transinit,
        tforms=load_slicetforms(transform_params.slicetforms_path),
        spacing_mm=float(config.registration.resolution_um) * 1e-3,
        template_native_shape=tuple(int(v) for v in atlas_volumes.template.shape),
        write_csv=write_csv_resolved,
    )
    results = run_annotation_import(
        specs,
        importer=importer,
        output_dir=output_dir,
        supports_masks=False,
        mask_unsupported_message=MASK_UNSUPPORTED,
    )
    if config.analysis.count_points:
        from lightsuite.analysis.cord_runner import maybe_run_cord_region_stats

        maybe_run_cord_region_stats(
            config,
            export_spaces=list(config.export.spaces),
        )
    return results
