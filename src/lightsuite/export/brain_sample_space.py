"""Brain pipeline sample-space export (atlas warped onto registration grid)."""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
from rich.console import Console

from lightsuite.analysis.division_map import ensure_division_map
from lightsuite.analysis.ontology import RegionTable, load_region_table
from lightsuite.analysis.region_stats import (
    concat_tidy,
    parcellation_result_to_tidy,
    write_region_stats_csv,
)
from lightsuite.atlas.io import load_atlas_volume
from lightsuite.atlas.registry import (
    atlas_display_provider_from_config,
    resolve_brain_atlas_from_config,
    uses_ccf_id_parcellation,
)
from lightsuite.config.models import BrainPipelineConfig
from lightsuite.export.parcellation import (
    compute_allen_parcellation,
    compute_perens_parcellation,
)
from lightsuite.export.sample_space import transform_atlas_volume_to_sample
from lightsuite.io.tiff_write import save_registration_volume
from lightsuite.preprocess.checkpoint import RegOptsCheckpoint
from lightsuite.registration.brain_register import TransformParamsCheckpoint
from lightsuite.registration.plots import (
    boundary_volume_from_annotation,
    save_registration_stage_previews,
)
from lightsuite.registration.volume import load_registration_volume

console = Console()

SAMPLE_SPACE_SUBDIR = "sample_space"
ANNOTATION_IN_SAMPLE = "annotation_in_sample_20um.tif"
TEMPLATE_IN_SAMPLE = "template_in_sample_20um.tif"
BOUNDARY_IN_SAMPLE = "annotation_boundary_in_sample_20um.tif"
DIVISION_IN_SAMPLE = "division_labels_in_sample_20um.tif"
MANIFEST_NAME = "sample_space_manifest.json"


def sample_space_dir(save_path: Path) -> Path:
    return save_path / "volume_registered" / SAMPLE_SPACE_SUBDIR


def export_brain_sample_space(
    config: BrainPipelineConfig,
    *,
    transform_params: TransformParamsCheckpoint,
    checkpoint: RegOptsCheckpoint,
    write_csv: bool,
    save_volume: bool,
    region_table: RegionTable | None,
) -> tuple[Path, dict[int, Path] | None]:
    """Warp atlas labels onto the registration grid and write sample-space stats."""
    save_path = config.sample.save_path.expanduser()
    out_dir = sample_space_dir(save_path)
    if save_volume or write_csv:
        out_dir.mkdir(parents=True, exist_ok=True)

    atlas = resolve_brain_atlas_from_config(config.atlas)
    spacing_mm = checkpoint.registres_um * 1e-3
    registres_um = float(checkpoint.registres_um)
    transformix_root = save_path / "transformix_sample_export_temp"
    transformix_root.mkdir(parents=True, exist_ok=True)

    console.print("Warping atlas volumes to sample registration grid...")
    t0 = time.perf_counter()

    av_native = load_atlas_volume(atlas.annotation_path)
    annotation_sample = transform_atlas_volume_to_sample(
        av_native,
        transform_params,
        save_path=save_path,
        spacing_mm=spacing_mm,
        temp_dir=transformix_root / "annotation",
        nearest=True,
    )

    tv_native = load_atlas_volume(atlas.template_path)
    template_sample = transform_atlas_volume_to_sample(
        tv_native,
        transform_params,
        save_path=save_path,
        spacing_mm=spacing_mm,
        temp_dir=transformix_root / "template",
        nearest=False,
    )

    if atlas.boundary_path is not None and atlas.boundary_path.is_file():
        boundary_native = load_atlas_volume(atlas.boundary_path)
        boundary_sample = transform_atlas_volume_to_sample(
            boundary_native,
            transform_params,
            save_path=save_path,
            spacing_mm=spacing_mm,
            temp_dir=transformix_root / "boundary",
            nearest=True,
        )
    else:
        boundary_sample = boundary_volume_from_annotation(annotation_sample)

    division_result = ensure_division_map(atlas)
    division_sample = transform_atlas_volume_to_sample(
        division_result.labels.astype(np.float32),
        transform_params,
        save_path=save_path,
        spacing_mm=spacing_mm,
        temp_dir=transformix_root / "division",
        nearest=True,
    )

    if save_volume:
        save_registration_volume(
            np.rint(annotation_sample).astype(np.uint16),
            out_dir / ANNOTATION_IN_SAMPLE,
        )
        save_registration_volume(
            np.clip(template_sample, 0, np.iinfo(np.uint16).max).astype(np.uint16),
            out_dir / TEMPLATE_IN_SAMPLE,
        )
        save_registration_volume(
            (boundary_sample > 0).astype(np.uint16),
            out_dir / BOUNDARY_IN_SAMPLE,
        )
        save_registration_volume(
            np.rint(division_sample).astype(np.uint16),
            out_dir / DIVISION_IN_SAMPLE,
        )

    channel_paths = {int(k): Path(v) for k, v in (checkpoint.regvolpaths or {}).items()}
    manifest = {
        "space": "sample",
        "grid": "registration",
        "shape_yxz": list(transform_params.regvolsize),
        "voxel_um": [registres_um, registres_um, registres_um],
        "straightened": False,
        "permute_sample_to_atlas": transform_params.permute_sample_to_atlas,
        "registration_canvas": transform_params.registration_canvas,
        "channel_paths": {str(k): str(v) for k, v in channel_paths.items()},
        "annotation_path": str(out_dir / ANNOTATION_IN_SAMPLE),
        "template_path": str(out_dir / TEMPLATE_IN_SAMPLE),
        "boundary_path": str(out_dir / BOUNDARY_IN_SAMPLE),
        "division_labels_path": str(out_dir / DIVISION_IN_SAMPLE),
    }
    manifest_path = out_dir / MANIFEST_NAME
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    region_stats_paths: dict[int, Path] | None = None
    tidy_frames: list[pd.DataFrame] = []
    if write_csv and atlas.supports_parcellation:
        ann_int = np.rint(annotation_sample).astype(np.int32)
        for ichan, volpath in sorted(channel_paths.items()):
            vol = load_registration_volume(volpath)
            if uses_ccf_id_parcellation(atlas):
                result = compute_perens_parcellation(
                    vol,
                    atlas,
                    transform_params.atlas_resolution_um,
                    annotation=ann_int,
                    voxel_um=registres_um,
                )
            else:
                result = compute_allen_parcellation(
                    vol,
                    atlas,
                    transform_params.atlas_resolution_um,
                    annotation=ann_int,
                    voxel_um=registres_um,
                )
            tidy = parcellation_result_to_tidy(
                result,
                region_table,
                sample=config.sample.name,
                channel=ichan,
                atlas=atlas.brain_atlas,
            )
            tidy_path = out_dir / f"chan{ichan:02d}_region_stats_sample.csv"
            write_region_stats_csv(tidy_path, tidy)
            if region_stats_paths is None:
                region_stats_paths = {}
            region_stats_paths[ichan] = tidy_path
            tidy_frames.append(tidy)

        if tidy_frames:
            combined = out_dir / "region_stats_sample.csv"
            write_region_stats_csv(combined, concat_tidy(tidy_frames))

    if save_volume and channel_paths:
        primary = min(channel_paths)
        vol = load_registration_volume(channel_paths[primary])
        hi = float(np.quantile(vol, 0.999))
        vol_u8 = np.clip(vol / max(hi, 1e-6) * 255.0, 0, 255).astype(np.uint8)
        save_registration_stage_previews(
            out_dir,
            config.sample.name,
            vol_u8,
            np.rint(annotation_sample),
            "export_sample",
            atlas_provider=atlas_display_provider_from_config(config.atlas),
        )

    console.print(
        f"Sample-space export done in {time.perf_counter() - t0:.1f}s under {out_dir}"
    )
    return manifest_path, region_stats_paths
