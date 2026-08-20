"""Convert Arivis Blob Finder feature tables to LightSuite points.csv."""

from __future__ import annotations

import csv
from pathlib import Path


def _find_com_columns(header: list[str]) -> tuple[int, int, int]:
    lowered = [h.lower() for h in header]

    def pick(axis: str) -> int:
        for i, h in enumerate(lowered):
            if axis in h and "center of mass" in h:
                return i
        for i, h in enumerate(lowered):
            if h.startswith(f"{axis} (px)"):
                return i
        msg = f"Could not find {axis} COM column in header: {header}"
        raise KeyError(msg)

    return pick("x"), pick("y"), pick("z")


def _read_table(path: Path) -> tuple[list[str], list[list[str]]]:
    path = path.expanduser()
    suffix = path.suffix.lower()
    if suffix in {".xlsx", ".xlsm"}:
        try:
            import openpyxl
        except ImportError as exc:
            msg = "Reading Arivis .xlsx requires openpyxl (uv add openpyxl)"
            raise ImportError(msg) from exc
        wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
        ws = wb.active
        rows_iter = ws.iter_rows(values_only=True)
        header_row = next(rows_iter, None)
        if header_row is None:
            msg = f"Empty Arivis spreadsheet: {path}"
            raise ValueError(msg)
        header = [str(c) if c is not None else "" for c in header_row]
        data: list[list[str]] = []
        for row in rows_iter:
            data.append(["" if c is None else str(c) for c in row])
        return header, data

    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle)
        rows = list(reader)
    if not rows:
        msg = f"Empty Arivis CSV: {path}"
        raise ValueError(msg)
    return rows[0], rows[1:]


def convert_arivis_features_to_points_csv(
    source: Path,
    output_csv: Path,
    *,
    index_base_in: int = 0,
    index_base_out: int = 1,
) -> int:
    """Write LightSuite ``points.csv`` from Arivis ``*-features.csv`` / ``.xlsx``.

    COM pixel columns are treated as ``index_base_in`` (default 0-based) and shifted
    to ``index_base_out`` (LightSuite 1-based).
    """
    header, data_rows = _read_table(source)
    ix, iy, iz = _find_com_columns(header)
    type_col = header.index("Type") if "Type" in header else None
    shift = float(index_base_out - index_base_in)

    output_csv = output_csv.expanduser()
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    n_written = 0
    with output_csv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["x", "y", "z"])
        writer.writeheader()
        for row in data_rows:
            if type_col is not None and type_col < len(row) and row[type_col] != "Segment":
                continue
            try:
                x = float(row[ix]) + shift
                y = float(row[iy]) + shift
                z = float(row[iz]) + shift
            except (TypeError, ValueError, IndexError):
                continue
            writer.writerow({"x": x, "y": y, "z": z})
            n_written += 1
    return n_written
