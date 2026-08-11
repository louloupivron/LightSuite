# MATLAB / Python parity scripts

Ad-hoc scripts used during the Python migration to compare intermediate results
against MATLAB checkpoints. **Not part of the installed `lightsuite` package.**

| Script | Purpose |
|--------|---------|
| `scratch_affine_matlab_parity.py` | Affine transform parity |
| `scratch_autopair_overlap.py` | Auto-pair overlap diagnostics |
| `scratch_check_fix.py` | Checkpoint fix utilities |
| `scratch_denoise_inject.py` | Denoise injection for parity |
| `scratch_inject_matlab_clouds.py` | Inject MATLAB point clouds |
| `scratch_inject_matlab_voluse.py` | Inject MATLAB volume usage |
| `scratch_mask_parity.py` | Mask parity checks |
| `scratch_parity_bcpd.py` | BCPD registration parity |
| `scratch_parity_diag.py` | General parity diagnostics |
| `scratch_parity_downsample.py` | Downsample parity |
| `scratch_pts.py` | Point cloud utilities |
| `scratch_warp_arivis_to_0p8x.py` | Arivis warp experiment |
| `scratch_yosi_matlab_parity.py` | Sample-specific MATLAB parity |

Run from the repository root, e.g.:

```bash
uv run python devtools/parity/scratch_parity_diag.py
```
