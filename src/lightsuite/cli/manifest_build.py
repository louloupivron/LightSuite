"""CLI for building multires pair manifests from vendor exports."""

from __future__ import annotations

import json
from pathlib import Path

import typer


def _parse_channels_json(raw: str) -> dict[str, dict[str, Path]]:
    data = json.loads(raw)
    if not isinstance(data, dict):
        msg = "channels JSON must be an object mapping channel name to {overview, roi} paths"
        raise typer.BadParameter(msg)
    out: dict[str, dict[str, Path]] = {}
    for name, paths in data.items():
        if not isinstance(paths, dict):
            raise typer.BadParameter(f"Channel {name!r} must be an object with overview/roi")
        out[str(name)] = {
            "overview": Path(paths["overview"]).expanduser(),
            "roi": Path(paths["roi"]).expanduser(),
        }
    return out


def build_manifest(
    vendor: str = typer.Option(
        ...,
        "--vendor",
        help="Acquisition vendor: mesospim or smartspim.",
    ),
    sample_name: str = typer.Option(..., "--sample-name", help="Sample identifier."),
    pair_label: str = typer.Option(..., "--pair-label", help="Unique pair label."),
    overview: Path | None = typer.Option(
        None,
        "--overview",
        help="Overview volume path (single-channel mesospim/smartspim).",
    ),
    roi: Path | None = typer.Option(None, "--roi", help="ROI volume path (single-channel)."),
    overview_meta: Path | None = typer.Option(
        None,
        "--overview-meta",
        help="mesoSPIM overview meta sidecar (required for stitched overview folders).",
    ),
    roi_meta: Path | None = typer.Option(None, "--roi-meta", help="mesoSPIM ROI meta sidecar."),
    channels_json: str | None = typer.Option(
        None,
        "--channels-json",
        help='Multichannel mesoSPIM paths as JSON, e.g. \'{"488":{"overview":"...","roi":"..."}}\'',
    ),
    reference_channel: str | None = typer.Option(
        None,
        "--reference-channel",
        help="Reference channel for multichannel manifests.",
    ),
    output: Path = typer.Option(
        ...,
        "--output",
        "-o",
        help="Output pair manifest JSON path.",
    ),
    lateral_flip_overview: str | None = typer.Option(
        None,
        "--lateral-flip-overview",
        help="mesoSPIM overview lateral_flip as '1,-1' (optional).",
    ),
    lateral_flip_roi: str | None = typer.Option(
        None,
        "--lateral-flip-roi",
        help="mesoSPIM ROI lateral_flip as '1,-1' (optional).",
    ),
) -> None:
    """Build a multires pair manifest JSON from mesoSPIM or SmartSPIM paths."""
    vendor_key = vendor.strip().lower()
    output = output.expanduser().resolve()

    if vendor_key == "mesospim":
        from lightsuite.mesospim.config_models import MesospimGeometryConfig
        from lightsuite.multires.vendor.mesospim import (
            build_mesospim_multichannel_pair_manifest,
            build_mesospim_pair_manifest,
        )

        overview_geometry = None
        roi_geometry = None
        if lateral_flip_overview:
            parts = [int(p.strip()) for p in lateral_flip_overview.split(",")]
            overview_geometry = MesospimGeometryConfig(lateral_flip=(parts[0], parts[1]))
        if lateral_flip_roi:
            parts = [int(p.strip()) for p in lateral_flip_roi.split(",")]
            roi_geometry = MesospimGeometryConfig(lateral_flip=(parts[0], parts[1]))

        if channels_json:
            channels = _parse_channels_json(channels_json)
            ref = reference_channel or next(iter(channels))
            manifest = build_mesospim_multichannel_pair_manifest(
                sample_name=sample_name,
                pair_label=pair_label,
                channels=channels,
                reference_channel=ref,
                overview_meta_path=overview_meta,
                overview_geometry=overview_geometry,
                roi_geometry=roi_geometry,
                output_manifest_path=output,
            )
        else:
            if overview is None or roi is None:
                raise typer.BadParameter("Single-channel mesospim requires --overview and --roi")
            manifest = build_mesospim_pair_manifest(
                sample_name=sample_name,
                pair_label=pair_label,
                overview_path=overview,
                roi_path=roi,
                overview_meta_path=overview_meta,
                roi_meta_path=roi_meta,
                overview_geometry=overview_geometry,
                roi_geometry=roi_geometry,
                output_manifest_path=output,
            )
        typer.echo(f"Wrote {output}")
        typer.echo(f"Reference channel: {manifest.reference_channel or 'n/a'}")
        if manifest.channels:
            typer.echo(f"Channels: {', '.join(manifest.channel_names())}")
        return

    if vendor_key == "smartspim":
        from lightsuite.multires.vendor.smartspim import build_smartspim_pair_manifest

        if overview is None or roi is None or overview_meta is None or roi_meta is None:
            raise typer.BadParameter(
                "smartspim requires --overview, --roi, --overview-meta, and --roi-meta"
            )
        build_smartspim_pair_manifest(
            sample_name=sample_name,
            pair_label=pair_label,
            overview_path=overview,
            roi_path=roi,
            overview_meta_path=overview_meta,
            roi_meta_path=roi_meta,
            output_manifest_path=output,
        )
        typer.echo(f"Wrote {output}")
        return

    raise typer.BadParameter("vendor must be mesospim or smartspim")
