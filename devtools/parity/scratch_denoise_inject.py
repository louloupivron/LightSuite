"""Denoise-only injection test: MATLAB downsample cloud -> Python pcdenoise.

On the MATLAB machine, add to ``extractSamplePoints.m`` just before ``pcdenoise``::

    export_dir = fullfile(fileparts(mfilename('fullpath')), '..', '..', 'voluse_export');
    if ~exist(export_dir, 'dir'), mkdir(export_dir); end
    writematrix(ptcloud.Location, fullfile(export_dir, 'downsample_pts.txt'));
    meta_down = struct('count', ptcloud.Count);
    fid = fopen(fullfile(export_dir, 'downsample_meta.json'), 'w');
    fprintf(fid, '%s', jsonencode(meta_down));
    fclose(fid);

Then run::

    uv run python scratch_denoise_inject.py --export-dir /path/to/voluse_export
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from lightsuite.registration.pc_downsample import pcdenoise

DEFAULT_EXPORT = Path(
    "/run/user/1002/gvfs/smb-share:server=h1data1.wysscenter.ch,share=computingdata"
    "/Alice/0003_CBT_WYSS_LIGHTSHEET/DATA/MESOSPIM/Marianna CMU/voluse_export"
)
MATLAB_DENOISE_TARGET = 524_010


def load_downsample(export_dir: Path) -> np.ndarray:
    pts_path = export_dir / "downsample_pts.txt"
    if not pts_path.is_file():
        raise FileNotFoundError(
            f"Missing {pts_path}. Export ptcloud.Location before pcdenoise in MATLAB."
        )
    pts = np.loadtxt(pts_path, dtype=np.float64, delimiter=",")
    if pts.ndim == 1:
        pts = pts.reshape(1, 3)
    return pts


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--export-dir", type=Path, default=DEFAULT_EXPORT)
    args = parser.parse_args()
    export_dir = args.export_dir.expanduser()

    down = load_downsample(export_dir)
    denoised = pcdenoise(down)
    removed = down.shape[0] - denoised.shape[0]

    print("=" * 72)
    print("DENOISE INJECTION TEST")
    print("=" * 72)
    print(f"Export dir: {export_dir}")
    print(f"Downsample input: {down.shape[0]:,}")
    print(f"Python pcdenoise: {denoised.shape[0]:,}  (removed {removed:,}, {100*removed/down.shape[0]:.2f}%)")
    print(f"MATLAB target:    {MATLAB_DENOISE_TARGET:,}  (removed ~854, ~0.17%)")

    meta_path = export_dir / "downsample_meta.json"
    if meta_path.is_file():
        meta = json.loads(meta_path.read_text())
        if "denoise_count" in meta:
            print(f"MATLAB measured:  {int(meta['denoise_count']):,}")


if __name__ == "__main__":
    main()
