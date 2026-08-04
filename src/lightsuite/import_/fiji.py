"""Convert FIJI / ImageJ point-tool Results.csv to LightSuite Sample Space v1 points CSV."""

from __future__ import annotations

import csv
from pathlib import Path


def _resolve_columns(fieldnames: list[str]) -> tuple[str, str, str]:
    field_map = {name.strip().lower(): name for name in fieldnames}
    missing = []
    x_key = field_map.get("x")
    y_key = field_map.get("y")
    z_key = field_map.get("slice") or field_map.get("z")
    if x_key is None:
        missing.append("x")
    if y_key is None:
        missing.append("y")
    if z_key is None:
        missing.append("slice or z")
    if missing:
        msg = f"FIJI Results CSV missing columns {missing}: {fieldnames}"
        raise KeyError(msg)
    return x_key, y_key, z_key


def fiji_to_native_xyz(
    x_val: float,
    y_val: float,
    z_val: float,
    *,
    voxel_um: list[float] | None,
) -> tuple[float, float, float]:
    """Convert FIJI point coordinates to 1-based native ``x, y, z`` voxel indices.

  When ``voxel_um`` is set, ``x_val`` and ``y_val`` are treated as calibrated units
  (typically µm from ImageJ spatial calibration). ``z_val`` is the 1-based slice index
  from the FIJI Results table.

  When ``voxel_um`` is ``None``, ``x_val`` and ``y_val`` are ImageJ pixel coordinates
  (fractional centers, e.g. 100.5 for pixel index 100).
    """
    if voxel_um is not None:
        vx, vy, _vz = (float(v) for v in voxel_um)
        x = int(x_val / vx) + 1
        y = int(y_val / vy) + 1
    else:
        x = float(x_val) + 0.5
        y = float(y_val) + 0.5
    z = float(z_val)
    return x, y, z


def convert_fiji_points_to_csv(
    source_csv: Path,
    output_csv: Path,
    *,
    voxel_um: list[float] | None = None,
) -> int:
    """Write a LightSuite ``points.csv`` from a FIJI point-tool Results export.

    Returns the number of points written.
    """
    source_csv = source_csv.expanduser()
    output_csv = output_csv.expanduser()
    output_csv.parent.mkdir(parents=True, exist_ok=True)

    with source_csv.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            msg = f"FIJI Results CSV must include a header row: {source_csv}"
            raise ValueError(msg)
        x_key, y_key, z_key = _resolve_columns(list(reader.fieldnames))
        rows = list(reader)

    n_written = 0
    with output_csv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["x", "y", "z"])
        writer.writeheader()
        for row in rows:
            try:
                x_val = float(row[x_key])
                y_val = float(row[y_key])
                z_val = float(row[z_key])
            except (TypeError, ValueError, KeyError):
                continue
            x, y, z = fiji_to_native_xyz(x_val, y_val, z_val, voxel_um=voxel_um)
            writer.writerow({"x": x, "y": y, "z": z})
            n_written += 1
    return n_written
