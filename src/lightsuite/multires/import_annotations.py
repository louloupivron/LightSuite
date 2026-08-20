"""Warp ROI-native segmentation into overview-native space after multires register.

Outputs stay in LightSuite Sample Space v1 on the *overview* grid, so they can be
handed straight to ``lightsuite brain import-annotations`` as ``points_csv`` /
``mask_tiff`` without any intermediate format.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from rich.console import Console

from lightsuite.config.models import AnnotationImportConfig
from lightsuite.import_.models import AnnotationImportResult, ImportedMask, ImportedPoints
from lightsuite.import_.orchestrator import (
    require_transformix,
    resolve_annotation_specs,
    resolve_write_csv,
    run_annotation_import,
)
from lightsuite.import_.sample_reference import SampleReference
from lightsuite.multires.checkpoint import MultiresRegOptsCheckpoint, multires_checkpoint_path
from lightsuite.multires.config_models import MultiresPipelineConfig
from lightsuite.multires.models import ManifestVolumeSpec
from lightsuite.multires.spec_geometry import (
    index_xyz_to_physical,
    physical_to_continuous_index_xyz,
)

console = Console()

OUTPUT_DIR_NAME = "annotations_in_overview"


def sample_reference_from_spec(spec: ManifestVolumeSpec, *, sample_name: str) -> SampleReference:
    """Build a Sample Space v1 reference from a manifest volume spec (no JSON needed)."""
    nz, ny, nx = (int(v) for v in spec.shape_zyx)
    return SampleReference.from_checkpoint(
        sample_name=sample_name,
        ny=ny,
        nx=nx,
        nz=nz,
        voxel_um=[float(v) for v in spec.spacing_um],
    )


def warp_points_roi_to_overview(
    points_xyz_1based: np.ndarray,
    *,
    roi_spec: ManifestVolumeSpec,
    overview_spec: ManifestVolumeSpec,
    transform_path: Path,
    temp_dir: Path,
) -> np.ndarray:
    """Map 1-based ROI voxel indices to 1-based overview voxel indices via elastix."""
    from lightsuite.registration.elastix.runner import run_transformix_physical_points

    pts = np.asarray(points_xyz_1based, dtype=float)
    if pts.size == 0:
        return np.zeros((0, 3), dtype=float)

    idx0 = pts[:, :3] - 1.0
    phys_roi = np.vstack([index_xyz_to_physical(roi_spec, row) for row in idx0])
    out_phys = run_transformix_physical_points(
        points_xyz=phys_roi,
        transform_path=transform_path,
        output_dir=temp_dir,
    )
    out_idx0 = np.vstack(
        [physical_to_continuous_index_xyz(overview_spec, row) for row in out_phys]
    )
    return out_idx0 + 1.0


def _points_in_bounds(points_xyz_1based: np.ndarray, spec: ManifestVolumeSpec) -> np.ndarray:
    nz, ny, nx = (int(v) for v in spec.shape_zyx)
    pts = np.asarray(points_xyz_1based, dtype=float)
    return (
        (pts[:, 0] >= 1)
        & (pts[:, 0] <= nx)
        & (pts[:, 1] >= 1)
        & (pts[:, 1] <= ny)
        & (pts[:, 2] >= 1)
        & (pts[:, 2] <= nz)
    )


def crop_reference_image(
    overview_spec: ManifestVolumeSpec,
    crop_start_index: list[int],
    crop_size_xyz: list[int],
):
    """Geometry-only reference grid matching the registration overlap crop."""
    import SimpleITK as sitk

    sx, sy, sz = (int(v) for v in crop_size_xyz)
    image = sitk.Image([sx, sy, sz], sitk.sitkUInt8)
    origin = index_xyz_to_physical(overview_spec, tuple(float(v) for v in crop_start_index))
    image.SetSpacing(tuple(float(v) for v in overview_spec.spacing_um))
    image.SetOrigin(tuple(float(v) for v in origin))
    image.SetDirection(tuple(float(v) for v in overview_spec.direction))
    return image


@dataclass
class MultiresAnnotationImporter:
    """Warp ROI-native annotations onto the overview native grid."""

    roi_spec: ManifestVolumeSpec
    overview_spec: ManifestVolumeSpec
    transform_paths: list[Path]
    crop_start_index: list[int]
    crop_size_xyz: list[int]
    output_dir: Path
    temp_dir: Path
    sample_name: str = "sample"
    write_csv: bool = True
    write_full_overview_canvas: bool = True

    @property
    def reference(self) -> SampleReference:
        return sample_reference_from_spec(self.roi_spec, sample_name=self.sample_name)

    def import_points(self, points: ImportedPoints, *, slug: str) -> AnnotationImportResult:
        overview_xyz = warp_points_roi_to_overview(
            points.coordinates,
            roi_spec=self.roi_spec,
            overview_spec=self.overview_spec,
            transform_path=self.transform_paths[-1],
            temp_dir=self.temp_dir / slug,
        )
        in_bounds = _points_in_bounds(overview_xyz, self.overview_spec)

        frame = pd.DataFrame(overview_xyz[:, :3], columns=["x", "y", "z"])
        if points.features is not None and points.features.shape[0] == overview_xyz.shape[0]:
            for idx in range(points.features.shape[1]):
                frame[f"feature_{idx + 1}"] = points.features[:, idx]
        csv_path = self.output_dir / f"{slug}_in_overview.csv"
        frame.to_csv(csv_path, index=False)

        npz_path = self.output_dir / f"{slug}_overview_coords.npz"
        np.savez_compressed(
            npz_path,
            overviewptcoords=overview_xyz.astype(np.float32),
            roiptcoords=points.coordinates.astype(np.float32),
        )

        nz, ny, nx = (int(v) for v in self.overview_spec.shape_zyx)
        console.print(
            f"[green]{points.label}[/green]: {points.coordinates.shape[0]} ROI points "
            f"→ {int(in_bounds.sum())}/{overview_xyz.shape[0]} inside the overview grid "
            f"({nx}×{ny}×{nz} XYZ)"
        )
        return AnnotationImportResult(
            label=points.label,
            kind="points",
            atlas_csv_path=csv_path,
            atlas_points_path=npz_path,
            n_input=int(points.coordinates.shape[0]),
            n_atlas=int(in_bounds.sum()),
            n_sample=int(overview_xyz.shape[0]),
        )

    def import_mask(self, mask: ImportedMask, *, slug: str) -> AnnotationImportResult:
        import SimpleITK as sitk

        from lightsuite.multires.registration import apply_elastix_transforms
        from lightsuite.multires.volume import apply_manifest_geometry, write_embedded_crop_canvas

        expected_zyx = tuple(int(v) for v in self.roi_spec.shape_zyx)
        # ImportedMask is (Y, X, Z); manifest specs are (Z, Y, X).
        mask_zyx = np.transpose(mask.volume, (2, 0, 1))
        if tuple(mask_zyx.shape) != expected_zyx:
            msg = (
                f"Mask shape (Z, Y, X)={tuple(mask_zyx.shape)} does not match the ROI "
                f"manifest {expected_zyx}. Segment on the ROI native grid before import."
            )
            raise ValueError(msg)

        moving = sitk.GetImageFromArray(np.asarray(mask_zyx > 0, dtype=np.uint8))
        moving = apply_manifest_geometry(moving, self.roi_spec)
        reference = crop_reference_image(
            self.overview_spec,
            self.crop_start_index,
            self.crop_size_xyz,
        )
        warped_crop = apply_elastix_transforms(
            moving,
            self.transform_paths,
            reference=reference,
            nearest=True,
            output_dtype=np.uint8,
        )
        crop_voxels = int(np.count_nonzero(sitk.GetArrayViewFromImage(warped_crop)))

        crop_path = self.output_dir / f"{slug}_in_overview_crop.tif"
        sitk.WriteImage(warped_crop, str(crop_path), useCompression=True)

        full_path: Path | None = None
        if self.write_full_overview_canvas:
            full_path = self.output_dir / f"{slug}_in_overview.tif"
            write_embedded_crop_canvas(
                self.overview_spec,
                warped_crop,
                self.crop_start_index,
                full_path,
                dtype=np.uint8,
            )

        console.print(
            f"[green]{mask.label}[/green]: mask warped to the overview overlap crop "
            f"({crop_voxels} foreground voxels)"
            + (f"; full canvas → {full_path.name}" if full_path is not None else "")
        )
        return AnnotationImportResult(
            label=mask.label,
            kind="mask",
            atlas_mask_path=full_path or crop_path,
            sample_mask_path=crop_path,
            n_input=int(np.count_nonzero(mask.volume)),
            n_atlas=crop_voxels,
            n_sample=crop_voxels,
        )


def _crop_size_from_checkpoint(checkpoint: MultiresRegOptsCheckpoint) -> list[int]:
    import SimpleITK as sitk

    if not checkpoint.cropped_overview_path:
        msg = "multires_regopts.json missing cropped_overview_path. Re-run 'multires register'."
        raise RuntimeError(msg)
    crop_ref = Path(checkpoint.cropped_overview_path).expanduser()
    if not crop_ref.is_file():
        msg = f"Missing overview crop reference {crop_ref}. Re-run 'multires register'."
        raise FileNotFoundError(msg)
    return [int(v) for v in sitk.ReadImage(str(crop_ref)).GetSize()]


def run_multires_import_annotations(
    config: MultiresPipelineConfig,
    *,
    annotations: list[AnnotationImportConfig] | None = None,
    write_csv: bool | None = None,
    full_overview_canvas: bool | None = None,
) -> list[AnnotationImportResult]:
    """Warp ROI-native annotations into overview-native space."""
    require_transformix()
    save_path = config.sample.save_path.expanduser()
    specs = resolve_annotation_specs(
        config.import_config,
        annotations,
        save_path=save_path,
    )
    write_csv_resolved = resolve_write_csv(config.import_config, write_csv)

    checkpoint_path = multires_checkpoint_path(save_path)
    if not checkpoint_path.is_file():
        msg = f"Missing {checkpoint_path}. Run 'lightsuite multires register' first."
        raise FileNotFoundError(msg)

    checkpoint = MultiresRegOptsCheckpoint.load(checkpoint_path)
    if not checkpoint.transform_paths:
        msg = "multires_regopts.json has no transform_paths. Run 'multires register' first."
        raise RuntimeError(msg)
    if not checkpoint.crop_start_index:
        msg = "multires_regopts.json has no crop_start_index. Run 'multires register' first."
        raise RuntimeError(msg)

    from lightsuite.multires.manifest import load_pair_manifest

    manifest = load_pair_manifest(checkpoint.pair_manifest_path)
    output_dir = save_path / OUTPUT_DIR_NAME
    temp_root = save_path / "import_annotations_temp"
    temp_root.mkdir(parents=True, exist_ok=True)

    if full_overview_canvas is None:
        full_overview_canvas = config.multires.registration.write_full_overview_canvas

    importer = MultiresAnnotationImporter(
        roi_spec=manifest.roi,
        overview_spec=manifest.overview,
        transform_paths=[Path(p) for p in checkpoint.transform_paths],
        crop_start_index=[int(v) for v in checkpoint.crop_start_index],
        crop_size_xyz=_crop_size_from_checkpoint(checkpoint),
        output_dir=output_dir,
        temp_dir=temp_root,
        sample_name=config.sample.name,
        write_csv=write_csv_resolved,
        write_full_overview_canvas=full_overview_canvas,
    )
    return run_annotation_import(specs, importer=importer, output_dir=output_dir)


__all__ = [
    "MultiresAnnotationImporter",
    "run_multires_import_annotations",
    "sample_reference_from_spec",
    "warp_points_roi_to_overview",
]
