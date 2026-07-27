"""Rebuild Multi_RES_SCANs pair manifests with validated SmartSPIM geometry."""

from __future__ import annotations

from pathlib import Path

from lightsuite.multires.vendor.smartspim import (
    SmartspimGeometryConfig,
    build_smartspim_pair_manifest,
    parse_smartspim_metadata,
)

ROOT = Path("/media/gbm/NVME2/ALICe-pipelines-data/Multi_RES_SCANs")
OUT = ROOT / "registration_results/converted"


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)

    cortex_tiles = parse_smartspim_metadata(
        ROOT / "9X/20260710_11_05_27_9X_4x2_cortex_cc/metadata.txt"
    ).tile_centers_stage
    build_smartspim_pair_manifest(
        sample_name="Multi_RES_SCANs",
        pair_label="cortex_9x_561",
        overview_path=ROOT / "1_6X/All_Channels",
        roi_path=ROOT / "9X/20260710_11_05_27_9X_4x2_cortex_cc/All_Channels",
        overview_meta_path=ROOT / "1_6X/metadata.txt",
        roi_meta_path=ROOT / "9X/20260710_11_05_27_9X_4x2_cortex_cc/metadata.txt",
        roi_tile_centers_stage=cortex_tiles,
        overview_geometry=SmartspimGeometryConfig(
            lateral_flip=(1, -1),
            apply_mosaic_pitch_correction=True,
        ),
        roi_geometry=SmartspimGeometryConfig(lateral_flip=(1, -1)),
        output_manifest_path=OUT / "cortex_9x_561_pair.json",
    )
    print("wrote", OUT / "cortex_9x_561_pair.json")

    cerebellum_meta = ROOT / "9X/20260710_10_04_38_9X_cerebellum/metadata.txt"
    cerebellum_tiles = parse_smartspim_metadata(cerebellum_meta).tile_centers_stage
    z_center = SmartspimGeometryConfig(lateral_flip=(1, -1), stage_z_is_center=True)
    build_smartspim_pair_manifest(
        sample_name="Multi_RES_SCANs",
        pair_label="cerebellum_9x_561",
        overview_path=ROOT / "1_6X/All_Channels",
        roi_path=ROOT / "9X/20260710_10_04_38_9X_cerebellum/All_Channels",
        overview_meta_path=ROOT / "1_6X/metadata.txt",
        roi_meta_path=cerebellum_meta,
        roi_tile_centers_stage=cerebellum_tiles,
        overview_geometry=z_center,
        roi_geometry=z_center,
        output_manifest_path=OUT / "cerebellum_9x_561_pair.json",
    )
    print("wrote", OUT / "cerebellum_9x_561_pair.json")

    build_smartspim_pair_manifest(
        sample_name="Multi_RES_SCANs",
        pair_label="single_fov_561",
        overview_path=ROOT / "1_6X/All_Channels",
        roi_path=ROOT / "9X/20260710_11_03_38_9X_single__FOV/Ex_561_Em_561F",
        overview_meta_path=ROOT / "1_6X/metadata.txt",
        roi_meta_path=ROOT / "9X/20260710_11_03_38_9X_single__FOV/metadata.txt",
        overview_geometry=SmartspimGeometryConfig(lateral_flip=(1, -1)),
        roi_geometry=SmartspimGeometryConfig(
            lateral_flip=(1, -1),
            stage_z_is_center=True,
            stage_origin_offset_um=(2000.0, 587.5, 0.0),
        ),
        output_manifest_path=OUT / "single_fov_561_pair.json",
    )
    print("wrote", OUT / "single_fov_561_pair.json")


if __name__ == "__main__":
    main()
