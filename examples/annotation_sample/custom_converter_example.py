"""Example custom annotation converter for ``import.converter.suite: custom``.

Contract — the module must define:

    convert_to_lightsuite(source, output, *, reference) -> dict

``reference`` is a :class:`~lightsuite.import_.sample_reference.SampleReference`.
Return at least ``{"format": "points_csv", "n_points": N}`` or
``{"format": "mask_tiff"}``. The convert-annotations stage always validates the
written file against ``sample_reference.json``.
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Any


def convert_to_lightsuite(
    source: Path,
    output: Path,
    *,
    reference: Any,
) -> dict[str, Any]:
    """Copy a native ``x,y,z`` CSV (or adapt your tool's columns here)."""
    source = Path(source)
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)

    with source.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            msg = f"CSV must have a header row: {source}"
            raise ValueError(msg)
        field_map = {name.strip().lower(): name for name in reader.fieldnames}
        for axis in ("x", "y", "z"):
            if axis not in field_map:
                msg = f"Missing column {axis!r} in {source}"
                raise KeyError(msg)
        rows = list(reader)

    n = 0
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["x", "y", "z"])
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    "x": float(row[field_map["x"]]),
                    "y": float(row[field_map["y"]]),
                    "z": float(row[field_map["z"]]),
                }
            )
            n += 1

    _ = reference  # available for bounds checks / voxel size if needed
    return {"format": "points_csv", "n_points": n}
