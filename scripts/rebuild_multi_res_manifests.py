"""Rebuild Multi_RES_SCANs pair manifests from SmartSPIM stage metadata."""

from __future__ import annotations

from pathlib import Path

from lightsuite.multires.vendor.smartspim import build_smartspim_pair_manifest

ROOT = Path("/media/gbm/NVME2/ALICe-pipelines-data/Multi_RES_SCANs")
OUT = ROOT / "registration_results/converted"

OVERVIEW = ROOT / "1_6X/All_Channels"
OVERVIEW_META = ROOT / "1_6X/metadata.txt"

PAIRS = {
    "cortex_9x_561": (
        ROOT / "9X/20260710_11_05_27_9X_4x2_cortex_cc/All_Channels",
        ROOT / "9X/20260710_11_05_27_9X_4x2_cortex_cc/metadata.txt",
    ),
    "cerebellum_9x_561": (
        ROOT / "9X/20260710_10_04_38_9X_cerebellum/All_Channels",
        ROOT / "9X/20260710_10_04_38_9X_cerebellum/metadata.txt",
    ),
    "single_fov_561": (
        ROOT / "9X/20260710_11_03_38_9X_single__FOV/Ex_561_Em_561F",
        ROOT / "9X/20260710_11_03_38_9X_single__FOV/metadata.txt",
    ),
}


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    for pair_label, (roi_path, roi_meta_path) in PAIRS.items():
        output_path = OUT / f"{pair_label}_pair.json"
        build_smartspim_pair_manifest(
            sample_name="Multi_RES_SCANs",
            pair_label=pair_label,
            overview_path=OVERVIEW,
            roi_path=roi_path,
            overview_meta_path=OVERVIEW_META,
            roi_meta_path=roi_meta_path,
            output_manifest_path=output_path,
        )
        print("wrote", output_path)


if __name__ == "__main__":
    main()
