"""Throwaway: do the landmark pairs themselves demand expansion or contraction?"""

from pathlib import Path

import numpy as np

PY = Path("/media/gbm/NVME2/ALICe-pipelines-data/Marianna_CMU/registered_perens/elastix_temp")
MAT = Path(
    "/run/user/1002/gvfs/smb-share:server=h1data1.wysscenter.ch,share=computingdata/Alice/"
    "0003_CBT_WYSS_LIGHTSHEET/DATA/MESOSPIM/Marianna CMU/elastix_temp/elastix_temp"
)


def read_pts(path: Path) -> tuple[str, np.ndarray]:
    lines = [ln.strip() for ln in path.read_text().splitlines() if ln.strip()]
    kind = lines[0].lower()
    n = int(lines[1])
    pts = np.array([[float(v) for v in ln.split()] for ln in lines[2 : 2 + n]])
    return kind, pts


def report(label: str, d: Path) -> None:
    kf, f = read_pts(d / "fixed.txt")
    km, m = read_pts(d / "moving.txt")
    print(f"=== {label} ===")
    print(f"  fixed.txt  : {kf}  n={len(f)}")
    print(f"  moving.txt : {km}  n={len(m)}")
    n = min(len(f), len(m))
    f, m = f[:n], m[:n]

    # T must map fixed -> moving. Required displacement:
    disp = m - f
    print(f"  required |T(f)-f| vox : median {np.median(np.linalg.norm(disp,axis=1)):.2f}"
          f"  p95 {np.percentile(np.linalg.norm(disp,axis=1),95):.2f}")

    c = f.mean(axis=0)
    rf = np.linalg.norm(f - c, axis=1)
    rm = np.linalg.norm(m - c, axis=1)
    dr = rm - rf
    print(f"  radius fixed  : mean {rf.mean():.1f} vox")
    print(f"  radius moving : mean {rm.mean():.1f} vox")
    print(f"  radial delta (moving-fixed) : mean {dr.mean():+.2f}  median {np.median(dr):+.2f} vox")
    print(f"  scale implied by landmarks  : {rm.mean()/rf.mean():.4f}")
    inward = float((dr < 0).mean())
    print(f"  pairs pulling T inward (=> annotation EXPANDS): {inward:.1%}")
    print(f"  fixed bbox  lo {np.round(f.min(0),0).tolist()} hi {np.round(f.max(0),0).tolist()}")
    print(f"  moving bbox lo {np.round(m.min(0),0).tolist()} hi {np.round(m.max(0),0).tolist()}")
    print()


report("Python", PY)
report("MATLAB", MAT)
print("scale < 1.0 => landmarks ask T to contract => warped annotation looks BIGGER.")
