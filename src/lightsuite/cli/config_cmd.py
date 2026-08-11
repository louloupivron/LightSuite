"""Config wizard, explain, and JSON Schema export."""

from __future__ import annotations

import json
from pathlib import Path

import typer
import yaml

from lightsuite.config.loader import (
    load_config,
    load_multires_config,
    load_spinal_config,
)

config_app = typer.Typer(help="Create and inspect pipeline YAML configs.")


def _prompt_path(label: str, default: str = "") -> str:
    value = typer.prompt(label, default=default or None)
    return str(value).strip()


def _prompt_floats(label: str, default: str) -> list[float]:
    raw = typer.prompt(label, default=default)
    parts = [p.strip() for p in raw.split(",")]
    return [float(p) for p in parts]


def _write_yaml(path: Path, data: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        yaml.safe_dump(data, sort_keys=False, default_flow_style=False),
        encoding="utf-8",
    )


@config_app.command("init")
def config_init(
    workflow: str = typer.Argument(
        ...,
        help="Workflow template: brain, spinal, or multires.",
    ),
    output: str = typer.Option(
        "my_sample.yaml",
        "--output",
        "-o",
        help="Output YAML path.",
    ),
) -> None:
    """Interactive wizard that writes a starter YAML from prompts."""
    out = Path(output).expanduser().resolve()
    if out.is_file() and not typer.confirm(f"{out} exists. Overwrite?", default=False):
        raise typer.Exit(code=1)

    workflow_key = workflow.strip().lower()
    if workflow_key == "brain":
        data = {
            "sample": {
                "name": _prompt_path("Sample name", "my_mouse"),
                "source": {
                    "format": "tiff_stack",
                    "path": _prompt_path("Stitched TIFF folder"),
                    "tiff_type": "channelperfile",
                },
                "scratch": _prompt_path("Scratch directory"),
                "save_path": _prompt_path("Results directory"),
                "voxel_um": _prompt_floats("Voxel size µm [x,y,z]", "5.26, 5.26, 5.0"),
            },
            "atlas": {
                "provider": "allen",
                "resolution_um": 10,
                "atlas_dir": _prompt_path("Allen atlas directory"),
            },
            "registration": {
                "resolution_um": 20,
                "channel_primary": 1,
                "channel_secondary": 2,
            },
            "detection": {"enabled": False},
            "compute": {"workers": 4},
            "export": {"registered_volume_format": "tiff", "write_pyramid": False},
        }
    elif workflow_key == "spinal":
        data = {
            "sample": {
                "name": _prompt_path("Sample name", "my_cord"),
                "source": {
                    "format": "tiff_stack",
                    "tiff_type": "planeperfile",
                    "channels": [_prompt_path("Channel 1 folder")],
                },
                "scratch": _prompt_path("Scratch directory"),
                "save_path": _prompt_path("Results directory"),
                "voxel_um": _prompt_floats("Voxel size µm [x,y,z]", "1.8, 1.8, 1.8"),
            },
            "atlas": {
                "atlas_dir": _prompt_path("Fiederling atlas directory (Segments.csv)"),
            },
            "registration": {
                "resolution_um": 20,
                "channel_primary": 1,
            },
            "compute": {"workers": 4},
        }
    elif workflow_key == "multires":
        data = {
            "sample": {
                "name": _prompt_path("Sample name", "my_multires"),
                "save_path": _prompt_path("Results directory"),
                "scratch": _prompt_path("Scratch directory"),
            },
            "multires": {
                "pair_label": _prompt_path("Pair label", "overview_roi"),
                "channels": {
                    "488": {
                        "overview": _prompt_path("Overview path (488)"),
                        "roi": _prompt_path("ROI path (488)"),
                    },
                },
                "registration": {
                    "reference_channel": "488",
                },
            },
        }
    else:
        raise typer.BadParameter("workflow must be brain, spinal, or multires")

    _write_yaml(out, data)
    typer.echo(f"Wrote {out}")
    typer.echo(f"Next: lightsuite doctor -c {out}")


@config_app.command("explain")
def config_explain(
    config: str = typer.Option(..., "--config", "-c", help="Pipeline YAML config."),
) -> None:
    """Print a human-readable summary of resolved paths and key settings."""
    path = Path(config).expanduser().resolve()
    raw = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    if "multires" in raw:
        cfg = load_multires_config(path)
        typer.echo(f"Workflow: multires  sample={cfg.sample.name!r}")
        typer.echo(f"  save_path: {cfg.sample.save_path}")
        typer.echo(f"  scratch:   {cfg.sample.scratch}")
        typer.echo(f"  pair_label: {cfg.multires.pair_label}")
        if cfg.multires.channels:
            typer.echo(f"  channels: {', '.join(cfg.multires.channels)}")
        if cfg.multires.pair_manifest:
            typer.echo(f"  pair_manifest: {cfg.multires.pair_manifest}")
        if _has_import_annotations(cfg):
            typer.echo(f"  import annotations: {len(cfg.import_.annotations)} layer(s)")
        return

    if "atlas" in raw and raw.get("atlas", {}).get("atlas_dir") and "provider" not in raw.get(
        "atlas", {}
    ):
        cfg = load_spinal_config(path)
        typer.echo(f"Workflow: spinal  sample={cfg.sample.name!r}")
        typer.echo(f"  save_path: {cfg.sample.save_path}")
        typer.echo(f"  scratch:   {cfg.sample.scratch}")
        typer.echo(f"  voxel_um:  {list(cfg.sample.voxel_um)}")
        typer.echo(f"  atlas_dir: {cfg.atlas.atlas_dir}")
        typer.echo(f"  tiff_type: {cfg.sample.source.tiff_type.value}")
        if cfg.import_ and cfg.import_.annotations:
            typer.echo(f"  import annotations: {len(cfg.import_.annotations)} layer(s)")
        return

    cfg = load_config(path)
    typer.echo(f"Workflow: brain  sample={cfg.sample.name!r}")
    typer.echo(f"  save_path: {cfg.sample.save_path}")
    typer.echo(f"  scratch:   {cfg.sample.scratch}")
    typer.echo(f"  voxel_um:  {list(cfg.sample.voxel_um)}")
    typer.echo(f"  atlas:     {cfg.atlas.provider.value} @ {cfg.atlas.resolution_um} µm")
    if cfg.atlas.atlas_dir:
        typer.echo(f"  atlas_dir: {cfg.atlas.atlas_dir}")
    typer.echo(f"  registres: {cfg.registration.resolution_um} µm")
    if cfg.import_ and cfg.import_.annotations:
        typer.echo(f"  import annotations: {len(cfg.import_.annotations)} layer(s)")


def _has_import_annotations(config) -> bool:
    import_cfg = getattr(config, "import_", None)
    return bool(import_cfg and import_cfg.annotations)


@config_app.command("schema")
def config_schema(
    workflow: str = typer.Argument(..., help="brain, spinal, or multires."),
    output: str | None = typer.Option(
        None,
        "--output",
        "-o",
        help="Write JSON Schema to this path (default: stdout).",
    ),
) -> None:
    """Export Pydantic JSON Schema for IDE autocomplete."""
    workflow_key = workflow.strip().lower()
    if workflow_key == "brain":
        from lightsuite.config.models import BrainPipelineConfig

        schema = BrainPipelineConfig.model_json_schema()
    elif workflow_key == "spinal":
        from lightsuite.config.models import SpinalCordPipelineConfig

        schema = SpinalCordPipelineConfig.model_json_schema()
    elif workflow_key == "multires":
        from lightsuite.multires.config_models import MultiresPipelineConfig

        schema = MultiresPipelineConfig.model_json_schema()
    else:
        raise typer.BadParameter("workflow must be brain, spinal, or multires")

    text = json.dumps(schema, indent=2)
    if output is None:
        typer.echo(text)
        return
    out = Path(output).expanduser()
    out.write_text(text, encoding="utf-8")
    typer.echo(f"Wrote {out.resolve()}")
