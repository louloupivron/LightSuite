"""Convert Imaris spot statistics exports to LightSuite Sample Space v1 points CSV."""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np


def read_imaris_spot_csv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    """Parse Imaris spot statistics CSV (skip preamble until ``Position X`` header)."""
    text = path.expanduser().read_text(encoding="utf-8-sig")
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    header_idx = next(i for i, line in enumerate(lines) if line.startswith("Position X"))
    header = [col.strip() for col in lines[header_idx].split(",") if col.strip()]
    rows: list[dict[str, str]] = []
    for line in lines[header_idx + 1 :]:
        parts = [part.strip() for part in line.split(",")]
        if len(parts) < len(header):
            parts.extend([""] * (len(header) - len(parts)))
        rows.append(dict(zip(header, parts, strict=False)))
    return header, rows


def _resolve_position_columns(header: list[str]) -> tuple[str, str, str]:
    field_map = {name.strip().lower(): name for name in header}
    for axis in ("x", "y", "z"):
        key = f"position {axis}"
        if key not in field_map and f"position {axis} [µm]" in field_map:
            field_map[key] = field_map[f"position {axis} [µm]"]
    missing = [axis for axis in ("x", "y", "z") if f"position {axis}" not in field_map]
    if missing:
        msg = f"Imaris CSV missing Position columns for {missing}: {header}"
        raise KeyError(msg)
    return field_map["position x"], field_map["position y"], field_map["position z"]


def imaris_um_to_native_xyz(
    x_um: float,
    y_um: float,
    z_um: float,
    voxel_um: list[float],
) -> tuple[float, float, float]:
    """Convert Imaris position values to 1-based native voxel indices.

    ``voxel_um`` is the size of one LightSuite native voxel in the same units as the
    Imaris Position columns. When the ``.ims`` was calibrated at 1 µm (or positions are
    already voxel indices labeled as µm), use ``[1, 1, 1]``.
    """
    vx, vy, vz = (float(v) for v in voxel_um)
    return (
        int(x_um / vx) + 1,
        int(y_um / vy) + 1,
        int(z_um / vz) + 1,
    )


def warn_if_positions_look_like_voxel_indices(
    source_csv: Path,
    *,
    voxel_um: list[float],
    shape_yxz: tuple[int, int, int] | None = None,
) -> str | None:
    """Return a warning when Position values look like voxel indices, not µm at ``voxel_um``."""
    header, rows = read_imaris_spot_csv(source_csv)
    if not rows:
        return None
    x_key, y_key, z_key = _resolve_position_columns(header)
    ums = []
    for row in rows[:5000]:
        try:
            ums.append((float(row[x_key]), float(row[y_key]), float(row[z_key])))
        except (TypeError, ValueError):
            continue
    if not ums:
        return None
    arr = np.asarray(ums, dtype=float)
    vmax = arr.max(axis=0)
    vx, vy, vz = (float(v) for v in voxel_um)
    # Expected physical extent if positions are true µm on the LightSuite grid
    if shape_yxz is not None:
        ny, nx, nz = (int(v) for v in shape_yxz)
        expected_um = np.array([nx * vx, ny * vy, nz * vz], dtype=float)
        # Indices fit the grid much better than µm at the stated voxel size
        as_index_fit = float(np.mean(np.minimum(vmax / np.array([nx, ny, nz], dtype=float), 1.0)))
        as_um_fit = float(np.mean(np.minimum(vmax / expected_um, 1.0)))
        if as_index_fit > 0.7 and as_um_fit < 0.55 and min(vx, vy, vz) > 1.01:
            return (
                f"Imaris Position maxima {vmax.round(1).tolist()} fit voxel indices on "
                f"grid {(ny, nx, nz)} better than µm at voxel_um={list(voxel_um)}. "
                "Re-run with --voxel-um matching Imaris Image Properties "
                "(often 1,1,1 when the .ims is 1 µm/voxel, or hybrid 1,1,<sample_z_um> "
                "when Imaris Z plane count differs from LightSuite nz)."
            )
    return None


def convert_imaris_spots_to_points_csv(
    source_csv: Path,
    output_csv: Path,
    *,
    voxel_um: list[float],
    component_name: str | None = None,
) -> int:
    """Write a LightSuite ``points.csv`` from an Imaris spot export.

    Returns the number of points written.
    """
    header, rows = read_imaris_spot_csv(source_csv)
    x_key, y_key, z_key = _resolve_position_columns(header)
    component_key = None
    for candidate in ("Component Name", "Collection", "Category"):
        if candidate in header:
            component_key = candidate
            break

    output_csv = output_csv.expanduser()
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    n_written = 0
    with output_csv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["x", "y", "z"])
        writer.writeheader()
        for row in rows:
            if component_name is not None and component_key is not None:
                if row.get(component_key, "").strip() != component_name:
                    continue
            try:
                x_um = float(row[x_key])
                y_um = float(row[y_key])
                z_um = float(row[z_key])
            except (TypeError, ValueError):
                continue
            x, y, z = imaris_um_to_native_xyz(x_um, y_um, z_um, voxel_um)
            writer.writerow({"x": x, "y": y, "z": z})
            n_written += 1
    return n_written


def list_imaris_component_names(path: Path) -> list[str]:
    """Return sorted unique Imaris component / collection labels in a spot CSV."""
    header, rows = read_imaris_spot_csv(path)
    for candidate in ("Component Name", "Collection", "Category"):
        if candidate in header:
            values = {row.get(candidate, "").strip() for row in rows}
            return sorted(v for v in values if v)
    return []


def convert_imaris_spots_by_component(
    source_csv: Path,
    output_dir: Path,
    *,
    voxel_um: list[float],
    label_prefix: str = "imaris",
) -> dict[str, Path]:
    """Write one ``points.csv`` per Imaris component name."""
    output_dir = output_dir.expanduser()
    output_dir.mkdir(parents=True, exist_ok=True)
    written: dict[str, Path] = {}
    for component in list_imaris_component_names(source_csv):
        slug = "".join(ch if ch.isalnum() or ch in "-_" else "_" for ch in component)
        out_path = output_dir / f"{label_prefix}_{slug}_points.csv"
        n = convert_imaris_spots_to_points_csv(
            source_csv,
            out_path,
            voxel_um=voxel_um,
            component_name=component,
        )
        if n > 0:
            written[component] = out_path
    return written
