"""Orchestrate external annotation import after brain registration."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from rich.console import Console
from skimage.transform import resize as sk_resize

from lightsuite.config.models import (
    AnnotationImportConfig,
    BrainPipelineConfig,
)
from lightsuite.export.atlas_space import transform_mask_to_atlas
from lightsuite.export.brain_export import _load_transform_params
from lightsuite.import_.models import AnnotationImportResult, ImportedMask, ImportedPoints
from lightsuite.import_.orchestrator import (
    require_transformix,
    resolve_annotation_specs,
    resolve_write_csv,
    run_annotation_import,
    slug_for_label,
)
from lightsuite.import_.sample_reference import (
    SampleReference,
    load_sample_reference,
    validate_mask_against_reference,
)
from lightsuite.import_.transform import (
    sample_points_to_registration_voxels,
    transform_points_to_atlas,
)
from lightsuite.preprocess.checkpoint import RegOptsCheckpoint
from lightsuite.preprocess.slice_ops import output_xy_shape, output_z_count
from lightsuite.registration.brain_paths import (
    IMPORT_ANNOTATIONS_TEMP,
    brain_imports_dir,
    brain_work_dir,
    cleanup_brain_work,
)
from lightsuite.registration.points import volume_indices_to_cloud_xyz

console = Console()


def _slug(label: str) -> str:
    return slug_for_label(label)


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


@dataclass
class BrainAnnotationImporter:
    """Warp native sample annotations into Perens/Allen atlas space."""

    transform_params: object
    checkpoint: RegOptsCheckpoint
    sample_reference: SampleReference
    output_dir: Path
    temp_dir: Path
    write_csv: bool = True

    @property
    def reference(self) -> SampleReference:
        return self.sample_reference

    def import_points(self, points: ImportedPoints, *, slug: str) -> AnnotationImportResult:
        atlas_coords = transform_points_to_atlas(
            points.coordinates,
            self.transform_params,
            registres_um=self.checkpoint.registres_um,
            temp_dir=self.temp_dir / slug,
            content_crop_start=self.checkpoint.content_crop_start,
        )

        reg_yxz = sample_points_to_registration_voxels(
            points.coordinates,
            self.transform_params,
            registres_um=self.checkpoint.registres_um,
            content_crop_start=self.checkpoint.content_crop_start,
        )
        # Sample-space overlays / Napari expect 1-based registration-grid xyz.
        reg_coords = volume_indices_to_cloud_xyz(reg_yxz) + 1.0

        npz_path = self.output_dir / f"{slug}_atlas_coords.npz"
        sample_npz_path = self.output_dir / f"{slug}_sample_coords.npz"
        np.savez_compressed(
            npz_path,
            atlasptcoords=atlas_coords.astype(np.float32),
            sampleptcoords=points.coordinates.astype(np.float32),
        )
        np.savez_compressed(
            sample_npz_path,
            regptcoords=reg_coords.astype(np.float32),
            sampleptcoords=points.coordinates.astype(np.float32),
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
            f"[green]{points.label}[/green]: {points.coordinates.shape[0]} sample points "
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
            n_sample=int(reg_coords.shape[0]),
        )

    def import_mask(self, mask: ImportedMask, *, slug: str) -> AnnotationImportResult:
        from lightsuite.io.tiff_write import save_registration_volume
        from lightsuite.registration.volume import permute_brain_volume

        voxel_um = [float(v) for v in self.checkpoint.voxel_um]
        validate_mask_against_reference(mask.volume.shape, voxel_um, self.sample_reference)

        registres_um = self.checkpoint.registres_um
        reg_shape = _registration_shape_native_yxz(mask.volume.shape, voxel_um, registres_um)
        mask_reg = _resample_mask_native_to_registration(
            mask.volume,
            target_shape_yxz=reg_shape,
            native_voxel_um=voxel_um,
            registres_um=registres_um,
        )

        permuted_shape = tuple(int(v) for v in self.transform_params.regvolsize)
        mask_sample = permute_brain_volume(
            mask_reg,
            self.transform_params.permute_sample_to_atlas,
        )
        if mask_sample.shape != permuted_shape:
            mask_sample = _resample_mask_native_to_registration(
                mask_sample,
                target_shape_yxz=permuted_shape,
                native_voxel_um=[registres_um, registres_um, registres_um],
                registres_um=registres_um,
            )

        atlas_mask = transform_mask_to_atlas(
            mask_reg,
            self.transform_params,
            permute=self.transform_params.permute_sample_to_atlas,
            spacing_mm=registres_um * 1e-3,
            temp_dir=self.temp_dir / slug,
        )

        out_path = self.output_dir / f"{slug}_registered_atlas.tif"
        save_registration_volume(atlas_mask, out_path)
        sample_mask_path = self.output_dir / f"{slug}_in_sample_20um.tif"
        save_registration_volume(mask_sample.astype(np.uint16), sample_mask_path)

        console.print(
            f"[green]{mask.label}[/green]: mask warped to atlas "
            f"({int(np.count_nonzero(atlas_mask))} foreground voxels) "
            f"and sample grid ({int(np.count_nonzero(mask_sample))} voxels)"
        )
        return AnnotationImportResult(
            label=mask.label,
            kind="mask",
            atlas_mask_path=out_path,
            sample_mask_path=sample_mask_path,
            n_input=int(np.count_nonzero(mask.volume)),
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
    require_transformix()
    save_path = config.sample.save_path.expanduser()
    specs = resolve_annotation_specs(
        config.import_config,
        annotations,
        save_path=save_path,
    )
    write_csv_resolved = resolve_write_csv(config.import_config, write_csv)

    regopts_path = save_path / "regopts.json"
    if not regopts_path.is_file():
        msg = f"Missing {regopts_path}. Run preprocess and register first."
        raise FileNotFoundError(msg)

    output_dir = brain_imports_dir(save_path)
    temp_root = brain_work_dir(save_path, IMPORT_ANNOTATIONS_TEMP)

    importer = BrainAnnotationImporter(
        transform_params=_load_transform_params(save_path),
        checkpoint=RegOptsCheckpoint.load(regopts_path),
        sample_reference=load_sample_reference(save_path),
        output_dir=output_dir,
        temp_dir=temp_root,
        write_csv=write_csv_resolved,
    )
    try:
        results = run_annotation_import(specs, importer=importer, output_dir=output_dir)
    finally:
        cleanup_brain_work(save_path, IMPORT_ANNOTATIONS_TEMP)
    if config.analysis.count_points:
        from lightsuite.analysis.brain_runner import maybe_refresh_brain_region_stats

        maybe_refresh_brain_region_stats(config, export_spaces=list(config.export.spaces))
    return results
