#!/usr/bin/env python3
"""Publish the OP39M2 mesoSPIM multichannel multires pair manifest."""

from __future__ import annotations

from pathlib import Path

from lightsuite.multires.vendor.mesospim import build_mesospim_multichannel_pair_manifest

DATA_ROOT = Path("/media/gbm/NVME2/ALICe-pipelines-data/spinal_cord/OP39M2")
SAVE_PATH = DATA_ROOT / "multiresolution_results"
CONVERTED_DIR = SAVE_PATH / "converted"
MANIFEST_PATH = CONVERTED_DIR / "OP39M2_spinal_cord_pair.json"

OVERVIEW_META = (
    DATA_ROOT / "2.5X/20250331_MES_OP39-M2_Mag2.5x_ch488_Tile0.tiff_meta.txt"
)


def main() -> None:
    CONVERTED_DIR.mkdir(parents=True, exist_ok=True)
    manifest = build_mesospim_multichannel_pair_manifest(
        sample_name="OP39M2",
        pair_label="spinal_cord_488_561",
        reference_channel="488",
        channels={
            "488": {
                "overview": DATA_ROOT / "2.5X/output/channel_488/RES(12486x2162x1163)",
                "roi": DATA_ROOT
                / "1.25X/20250331_MES_OP39-M2_Mag1.25x_Tile0_Ch488_Sh0_Rot0.tiff",
            },
            "561": {
                "overview": DATA_ROOT / "2.5X/output/channel_561/RES(12486x2162x1163)",
                "roi": DATA_ROOT
                / "1.25X/20250331_MES_OP39-M2_Mag1.25x_Tile0_Ch561_Sh0_Rot0.tiff",
            },
        },
        overview_meta_path=OVERVIEW_META,
        output_manifest_path=MANIFEST_PATH,
    )
    print(f"Wrote {MANIFEST_PATH}")
    print(f"Reference channel: {manifest.reference_channel}")
    print(f"Channels: {', '.join(manifest.channel_names())}")
    print(f"Overview shape (488): {manifest.overview.shape_zyx}")


if __name__ == "__main__":
    main()
