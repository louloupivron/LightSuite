"""Workflow decision tree for the CLI."""

from __future__ import annotations

import typer

WORKFLOW_GUIDE = """
LightSuite Python workflows
===========================

brain
  Whole-brain lightsheet → Allen or Perens atlas registration.
  Template: examples/brain_lightsheet.yaml
  Run:      lightsuite brain run -c my.yaml --through match-points

spinal
  Spinal cord lightsheet → Fiederling atlas + region stats.
  Template: examples/spinal_cord.yaml
  Run:      lightsuite spinal run -c my.yaml

multires
  Overview ↔ high-resolution ROI (mesoSPIM, SmartSPIM) without atlas.
  Template: declare channels in YAML or build a manifest first.
  Run:      lightsuite multires run -c my.yaml

Stay on MATLAB for
  • Widefield coronal slices (demos/ls_analyze_slice_volume.m)
  • Built-in cell detection
  • CZI native reader
  • Spinal cohort NNMF analysis

Combined multires → brain
  See docs/usage_combined_multires_brain.md
"""


def print_workflow_list() -> None:
    typer.echo(WORKFLOW_GUIDE.strip())


workflow_app = typer.Typer(help="Choose a LightSuite workflow.")


@workflow_app.command("list")
def workflow_list() -> None:
    """Print which workflow to use (brain vs spinal vs multires vs MATLAB)."""
    print_workflow_list()
