"""Convert SmartSPIM / LCT cell-detection JSON to LightSuite Sample Space v1 points CSV."""

from __future__ import annotations

import csv
import json
from pathlib import Path


def smartspim_zyx_to_native_xyz(
    z: float,
    y: float,
    x: float,
    *,
    index_base_in: int = 0,
    index_base_out: int = 1,
) -> tuple[float, float, float]:
    """Convert SmartSPIM detection ``[z, y, x]`` indices to LightSuite ``x, y, z``.

    Detection JSON from the SmartSPIM / LCT cell-detection stack stores **0-based**
    voxel indices as ``[z, y, x]`` on the same grid as the stitched ``All_Channels``
    TIFF stack. LightSuite expects **1-based** ``x, y, z``.
    """
    shift = float(index_base_out - index_base_in)
    return float(x) + shift, float(y) + shift, float(z) + shift


def load_smartspim_points_json(source_json: Path) -> list[tuple[float, float, float]]:
    """Load ``[[z, y, x], …]`` from a SmartSPIM detection points JSON file."""
    source_json = source_json.expanduser()
    if not source_json.is_file():
        msg = f"SmartSPIM points JSON not found: {source_json}"
        raise FileNotFoundError(msg)

    data = json.loads(source_json.read_text(encoding="utf-8"))
    if not isinstance(data, list):
        msg = f"Expected a JSON list of [z, y, x] triples in {source_json}"
        raise ValueError(msg)

    points: list[tuple[float, float, float]] = []
    for i, item in enumerate(data):
        if not isinstance(item, (list, tuple)) or len(item) != 3:
            msg = f"Point {i} in {source_json} is not a length-3 [z, y, x] list: {item!r}"
            raise ValueError(msg)
        z, y, x = (float(item[0]), float(item[1]), float(item[2]))
        points.append((z, y, x))
    return points


def convert_smartspim_points_json_to_csv(
    source_json: Path,
    output_csv: Path,
    *,
    index_base_in: int = 0,
    index_base_out: int = 1,
) -> int:
    """Write LightSuite ``points.csv`` (``x,y,z``) from SmartSPIM ``[[z,y,x],…]`` JSON.

    Returns the number of points written.
    """
    output_csv = output_csv.expanduser()
    output_csv.parent.mkdir(parents=True, exist_ok=True)

    points_zyx = load_smartspim_points_json(source_json)
    with output_csv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["x", "y", "z"])
        writer.writeheader()
        for z, y, x in points_zyx:
            xo, yo, zo = smartspim_zyx_to_native_xyz(
                z, y, x, index_base_in=index_base_in, index_base_out=index_base_out
            )
            writer.writerow({"x": xo, "y": yo, "z": zo})
    return len(points_zyx)
