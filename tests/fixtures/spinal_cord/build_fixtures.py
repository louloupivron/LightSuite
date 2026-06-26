"""Build synthetic spinal cord test fixtures."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile

FIXTURE_ROOT = Path(__file__).resolve().parent


def build_fixtures(root: Path | None = None) -> None:
    root = FIXTURE_ROOT if root is None else Path(root)
    atlas_dir = root / "atlas"
    sample_dir = root / "sample"
    parity_dir = root / "parity"
    atlas_dir.mkdir(parents=True, exist_ok=True)
    sample_dir.mkdir(parents=True, exist_ok=True)
    parity_dir.mkdir(parents=True, exist_ok=True)

    ny, nx, nz = 40, 30, 24
    template = np.zeros((ny, nx, nz), dtype=np.float32)
    annotation = np.zeros((ny, nx, nz), dtype=np.uint16)
    for z in range(nz):
        cy, cx = ny // 2 + z // 5, nx // 2
        rr, cc = np.ogrid[:ny, :nx]
        mask = (rr - cy) ** 2 + (cc - cx) ** 2 < 36
        template[mask, z] = 1000 + z * 10
        annotation[mask, z] = 1 + (z % 5)

    tifffile.imwrite(atlas_dir / "Template.tif", template)
    tifffile.imwrite(atlas_dir / "Annotation.tif", annotation)
    pd.DataFrame({"Segment": ["C1", "C2"], "Start": [0, 12], "End": [11, 23]}).to_csv(
        atlas_dir / "Segments.csv",
        index=False,
    )
    pd.DataFrame(
        {
            "id": [1, 2],
            "name": ["gray matter", "lateral funiculus"],
            "children_IDs": ["1", "2"],
        }
    ).to_csv(atlas_dir / "Atlas_Regions.csv", index=False)

    n_long, n_y, n_x, n_z = 48, 40, 30, 20
    vol = np.zeros((n_long, n_y, n_x), dtype=np.uint16)
    for z in range(n_z):
        cy, cx = n_y // 2, n_x // 2 + z // 3
        for x in range(n_long):
            rr, cc = np.ogrid[:n_y, :n_x]
            mask = (rr - cy) ** 2 + (cc - cx) ** 2 < 49
            vol[x, mask] = 800 + z * 20
    tifffile.imwrite(sample_dir / "chan01_signal.tif", vol)

    plane_dir = sample_dir / "planeperfile"
    plane_dir.mkdir(parents=True, exist_ok=True)
    n_slices = 16
    for iz in range(n_slices):
        plane = np.zeros((n_y, n_x), dtype=np.uint16)
        cy, cx = n_y // 2, n_x // 2 + iz // 2
        rr, cc = np.ogrid[:n_y, :n_x]
        mask = (rr - cy) ** 2 + (cc - cx) ** 2 < 36
        plane[mask] = 600 + iz * 15
        tifffile.imwrite(plane_dir / f"slice_{iz:04d}.tif", plane)

    n = 12
    user_cen = np.column_stack([np.full(n, 16.0), np.full(n, 21.0)])
    user_ant = np.column_stack([np.full(n, 16.0), np.full(n, 13.0)])
    user_pos = np.column_stack([np.full(n, 16.0), np.full(n, 29.0)])
    from lightsuite.registration.straightening_optimizer import run_straightening_optimizer

    fit = run_straightening_optimizer(user_cen, user_ant, user_pos)
    parity = {
        "user_cen": user_cen.tolist(),
        "user_ant": user_ant.tolist(),
        "user_pos": user_pos.tolist(),
        "fit_x": fit["fit_x"].tolist(),
        "fit_y": fit["fit_y"].tolist(),
        "fit_theta": fit["fit_theta"].tolist(),
    }
    (parity_dir / "align_out_reference.json").write_text(json.dumps(parity, indent=2), encoding="utf-8")


if __name__ == "__main__":
    build_fixtures()
