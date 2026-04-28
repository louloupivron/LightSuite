#!/usr/bin/env python3
"""
Inspect Gubra / Perens LSFM atlas NIfTIs and ARA2 CSV; save figures and a text report.

Run from repo root (or this folder) with uv:
  cd scripts/lsfm_atlas_inspect && uv sync && uv run python inspect_lsfm_atlas.py --atlas-root /path/to/LSFM-mouse-brain-atlas-master

If gubra_template_olf.nii.gz lives in LSFM_atlas_files/perens/ but gubra_ano_olf.nii.gz
is in LSFM_atlas_files/, pass --atlas-root to the repo root; the script searches common layouts.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np
import pandas as pd

# Documented isotropic resolution (um) for Brainglobe / Perens packaging; NIfTI headers
# in this distribution often have placeholder zooms (1.0).
PERENS_ISOTROPIC_UM = 20.0


def find_atlas_files(atlas_root: Path) -> tuple[Path, Path, Path | None]:
    """
    Return (template_path, annotation_path, structures_csv_path).

    Tries:
      LSFM_atlas_files/{gubra_template_olf.nii.gz, gubra_ano_olf.nii.gz}
      LSFM_atlas_files/perens/gubra_template_olf.nii.gz with anno in parent
    """
    root = atlas_root.expanduser().resolve()
    candidates = [
        root / "LSFM_atlas_files",
        root,
    ]
    sub = root / "LSFM_atlas_files" / "perens"

    tpl_name = "gubra_template_olf.nii.gz"
    ann_name = "gubra_ano_olf.nii.gz"
    csv_name = "ARA2_annotation_info_avail_regions.csv"

    tpl = ann = None
    for base in candidates:
        a = base / ann_name
        t = base / tpl_name
        if a.is_file():
            ann = a
        if t.is_file():
            tpl = t
    if ann is None:
        raise FileNotFoundError(f"Could not find {ann_name} under {root}")
    if tpl is None and sub.is_dir():
        t = sub / tpl_name
        if t.is_file():
            tpl = t
    if tpl is None:
        raise FileNotFoundError(
            f"Could not find {tpl_name}. Expected next to {ann_name} or under LSFM_atlas_files/perens/."
        )

    csv_path = None
    for base in (ann.parent, ann.parent / "perens", root / "LSFM_atlas_files"):
        c = base / csv_name
        if c.is_file():
            csv_path = c
            break

    return tpl, ann, csv_path


def allen_lr_slab_hypothesis(ann: np.ndarray) -> dict:
    """
    Allen parcellation TIFF/NIfTI exports sometimes concatenate L/R hemispheres along one axis.
    Test each axis with even length: high label-set overlap between halves suggests mirrored slabs.
    """
    best = {"axis": None, "jaccard": -1.0}
    per_axis = []
    for axis in range(3):
        n = ann.shape[axis]
        if n % 2 != 0:
            per_axis.append({"axis": axis, "skipped": True, "reason": "odd length"})
            continue
        h = n // 2
        if axis == 0:
            left, right = ann[:h, :, :], ann[h:, :, :]
        elif axis == 1:
            left, right = ann[:, :h, :], ann[:, h:, :]
        else:
            left, right = ann[:, :, :h], ann[:, :, h:]
        ul, ur = np.unique(left), np.unique(right)
        inter = np.intersect1d(ul, ur)
        uni = np.union1d(ul, ur)
        jacc = float(len(inter) / max(len(uni), 1))
        per_axis.append(
            {
                "axis": axis,
                "jaccard_unique_labels": jacc,
                "n_unique_half_a": int(len(ul)),
                "n_unique_half_b": int(len(ur)),
            }
        )
        if jacc > best["jaccard"]:
            best = {"axis": axis, "jaccard": jacc}

    plausible = best["jaccard"] > 0.95 and best["axis"] is not None
    note = (
        "LightSuite Allen code uses size(av,3)/2; high jaccard on one axis would match that style."
    )
    if best["axis"] is None:
        note += " All axes have odd length — cannot split halves evenly; treat as single-volume atlas."
    return {
        "allen_lr_slab_plausible": plausible,
        "best_half_split_axis": best["axis"],
        "best_jaccard_unique_labels": best["jaccard"],
        "per_axis": per_axis,
        "note": note,
    }


def midpoint_report(ann: np.ndarray, tissue: np.ndarray | None) -> dict:
    """Rough COM along X for hemisphere-related workflows (ML axis often dim 1 in RAS-ish volumes)."""
    mask = ann > 0 if tissue is None else tissue > 0
    if not np.any(mask):
        return {}
    coords = np.argwhere(mask)
    com = coords.mean(axis=0)
    return {
        "brain_center_of_mass_ijk": tuple(float(x) for x in com),
        "volume_shape_ijk": tuple(int(x) for x in ann.shape),
    }


def plot_mosaic_slices(
    tpl: np.ndarray,
    ann: np.ndarray,
    out_png: Path,
    fixed_um: float,
) -> None:
    """Three rows: template / annotation / overlay; columns = cuts along dim0,1,2 mid."""
    fig, axes = plt.subplots(3, 3, figsize=(11, 10))
    ann_masked = np.ma.masked_where(ann == 0, ann).astype(np.float64)
    vmax_t = float(np.percentile(tpl[tpl > 0], 99.5)) if np.any(tpl > 0) else 1.0
    vmax_t = max(vmax_t, 1.0)

    for col, dim_name in enumerate(("axis0", "axis1", "axis2")):
        mid = [s // 2 for s in tpl.shape]
        if col == 0:
            sl_tpl, sl_ann = tpl[mid[0], :, :], ann[mid[0], :, :]
        elif col == 1:
            sl_tpl, sl_ann = tpl[:, mid[1], :], ann[:, mid[1], :]
        else:
            sl_tpl, sl_ann = tpl[:, :, mid[2]], ann[:, :, mid[2]]

        axes[0, col].imshow(sl_tpl, cmap="gray", vmin=0, vmax=vmax_t)
        axes[0, col].set_title(f"template {dim_name} mid")
        axes[0, col].axis("off")

        im = axes[1, col].imshow(sl_ann, cmap="nipy_spectral", interpolation="nearest")
        axes[1, col].set_title(f"labels {dim_name} mid")
        axes[1, col].axis("off")

        axes[2, col].imshow(sl_tpl, cmap="gray", vmin=0, vmax=vmax_t)
        axes[2, col].imshow(
            np.ma.masked_where(sl_ann == 0, sl_ann),
            cmap="nipy_spectral",
            alpha=0.45,
            interpolation="nearest",
        )
        axes[2, col].set_title("overlay")
        axes[2, col].axis("off")

    fig.suptitle(
        f"LSFM atlas (assume {fixed_um:g} µm isotropic if header zooms are placeholders)",
        fontsize=11,
    )
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    plt.close(fig)


def plot_label_counts(ids: np.ndarray, counts: np.ndarray, out_png: Path, top_n: int = 40) -> None:
    o = np.argsort(-counts)[:top_n]
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.barh(range(top_n), counts[o][::-1])
    ax.set_yticks(range(top_n))
    ax.set_yticklabels([str(int(ids[i])) for i in o[::-1]], fontsize=7)
    ax.set_xlabel("voxel count")
    ax.set_title(f"Top {top_n} label IDs by voxel count")
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    plt.close(fig)


def plot_csv_coverage(df: pd.DataFrame, in_vol: set[int], out_png: Path) -> None:
    ids_csv = set(df["id"].astype(np.int64).tolist())
    missing = sorted(ids_csv - in_vol)
    extra_in_vol = sorted(in_vol - ids_csv)
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.bar(
        ["CSV regions", "labels present in volume", "CSV ∩ volume", "CSV only (no voxels)", "Vol only (not in CSV)"],
        [
            len(ids_csv),
            len(in_vol),
            len(ids_csv & in_vol),
            len(missing),
            len(extra_in_vol),
        ],
        color=["#89c", "#8c9", "#4a4", "#c66", "#cc6"],
    )
    ax.set_ylabel("count")
    ax.set_title("ARA2 CSV vs annotation volume label IDs")
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    plt.close(fig)


def main() -> None:
    p = argparse.ArgumentParser(description="Inspect LSFM Gubra atlas files")
    p.add_argument(
        "--atlas-root",
        type=Path,
        default=None,
        help="Path to LSFM-mouse-brain-atlas-master (folder containing LSFM_atlas_files)",
    )
    p.add_argument(
        "--out-dir",
        type=Path,
        default=None,
        help="Output directory for PNG + report (default: ./lsfm_atlas_inspect_out under cwd)",
    )
    p.add_argument(
        "--resolution-um",
        type=float,
        default=PERENS_ISOTROPIC_UM,
        help="Physical isotropic resolution in µm (for reporting; NIfTI zooms may be wrong)",
    )
    args = p.parse_args()
    if args.atlas_root is None:
        raise SystemExit("Pass --atlas-root to your LSFM-mouse-brain-atlas-master folder.")

    out_dir = args.out_dir or (Path.cwd() / "lsfm_atlas_inspect_out")
    out_dir.mkdir(parents=True, exist_ok=True)

    tpl_path, ann_path, csv_path = find_atlas_files(args.atlas_root)
    print("Template:", tpl_path)
    print("Annotation:", ann_path)
    print("Structures CSV:", csv_path)

    tpl_img = nib.load(str(tpl_path))
    ann_img = nib.load(str(ann_path))
    tpl = np.asarray(tpl_img.dataobj)
    ann = np.asarray(ann_img.dataobj)

    zooms = tpl_img.header.get_zooms()[:3]
    report_lines = []
    report_lines.append("=== LSFM atlas inspection ===\n")
    report_lines.append(f"Template path: {tpl_path}")
    report_lines.append(f"Annotation path: {ann_path}")
    report_lines.append(f"Shape (both): {tpl.shape} {ann.shape}")
    report_lines.append(f"Match shapes: {tpl.shape == ann.shape}")
    report_lines.append(f"Template dtype: {tpl.dtype}, annotation dtype: {ann.dtype}")
    report_lines.append(f"NIfTI zooms (first 3): {zooms}")
    report_lines.append(
        f"Interpretation: if zooms are ~1.0, prefer documented isotropic spacing "
        f"({args.resolution_um:g} µm) from Brainglobe packaging / publication.\n"
    )

    u, c = np.unique(ann, return_counts=True)
    in_vol = set(int(x) for x in u.tolist())
    nz = int(np.sum(ann > 0))
    report_lines.append(f"Unique label values (incl. 0): {len(u)}, nonzero voxels: {nz}")
    report_lines.append(f"Label id min/max (nonzero): {ann[ann>0].min()} … {ann[ann>0].max()}")

    lr = allen_lr_slab_hypothesis(ann)
    report_lines.append("\n=== Allen-style half-volume (possible L|R slab) hypothesis ===")
    for k, v in lr.items():
        if k != "per_axis":
            report_lines.append(f"  {k}: {v}")
    report_lines.append("  per_axis:")
    for row in lr.get("per_axis", []):
        report_lines.append(f"    {row}")

    com = midpoint_report(ann, None)
    report_lines.append("\n=== Spatial summary ===")
    for k, v in com.items():
        report_lines.append(f"  {k}: {v}")

    plot_mosaic_slices(tpl, ann, out_dir / "fig_geometry_mosaic.png", args.resolution_um)
    plot_label_counts(u[u > 0], c[u > 0], out_dir / "fig_label_counts_top.png")

    if csv_path is not None:
        df = pd.read_csv(csv_path)
        ids_vol_nonzero = in_vol - {0}
        plot_csv_coverage(df, ids_vol_nonzero, out_dir / "fig_csv_vs_volume_ids.png")
        ids_csv = set(df["id"].astype(np.int64).tolist())
        report_lines.append("\n=== ARA2 CSV vs volume ===")
        report_lines.append(f"  Rows in CSV: {len(df)}, unique ids in CSV: {len(ids_csv)}")
        report_lines.append(f"  Label IDs in volume (nonzero): {len(ids_vol_nonzero)}")
        report_lines.append(f"  Intersection: {len(ids_csv & ids_vol_nonzero)}")
        report_lines.append(f"  CSV ids with zero voxels in this volume: {len(ids_csv - ids_vol_nonzero)}")
        report_lines.append(f"  Volume ids not listed in CSV: {len(ids_vol_nonzero - ids_csv)}")

    report_lines.append("\n=== Recommended LightSuite / Brainglobe strategy ===\n")
    report_lines.append(
        "1) Use gubra_template_olf.nii.gz + gubra_ano_olf.nii.gz from the same release; "
        "if your download splits files, put BOTH in one folder (or symlink) and set opts.atlas_dir.\n"
    )
    report_lines.append(
        "2) Do NOT reuse Allen heatmap parcellation that splits annotation along an axis into two slabs "
        "(e.g. size(av,3)/2). This LSFM volume is a single 3D grid; see half-slab hypothesis metrics "
        "above.\n"
    )
    report_lines.append(
        "3) Regional intensities for Brainglobe: aggregate by annotation ID with "
        "accumarray / pandas groupby on registered samples; map IDs to names via "
        "ARA2_annotation_info_avail_regions.csv (or Brainglobe structures.json from "
        "perens_lsfm_mouse packaging).\n"
    )
    report_lines.append(
        "4) Hemisphere separation: use spatial masks (split along ML axis through COM) "
        "or ontology leaves that are side-specific—do not assume Allen slab layout.\n"
    )
    report_lines.append(
        f"5) Use opts.atlasres = {args.resolution_um:g} with LightSuite when matching this atlas.\n"
    )

    report_text = "\n".join(report_lines)
    print("\n" + report_text)
    (out_dir / "report.txt").write_text(report_text, encoding="utf-8")
    print(f"\nWrote figures and report under: {out_dir.resolve()}")


if __name__ == "__main__":
    main()
