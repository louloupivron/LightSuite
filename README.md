# LightSuite

[![Documentation Status](https://readthedocs.org/projects/lightsuite/badge/?version=latest)](https://lightsuite.readthedocs.io/en/latest/)
[![Python](https://img.shields.io/badge/Python-3.11%2B-blue.svg)](https://www.python.org/)
[![MATLAB](https://img.shields.io/badge/MATLAB-R2022b%2B-blue.svg)](https://www.mathworks.com/products/matlab.html)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

**LightSuite** registers large microscopy volumes to standard anatomical atlases and exports atlas-space intensities and cell coordinates. It is designed for 100 GB+ datasets from lightsheet and widefield imaging.

**Documentation:** [lightsuite.readthedocs.io](https://lightsuite.readthedocs.io/en/latest/)

---

## Python pipeline (brain lightsheet)

The Python version runs on branch **`feature/python-migration`** as a CLI with optional Napari GUIs. No MATLAB license required.

```bash
# Install
uv sync --extra dev --extra gui

# Verify environment
uv run lightsuite doctor -c examples/brain_lightsheet.yaml

# Brain pipeline
uv run lightsuite brain preprocess           -c my_sample.yaml
uv run lightsuite brain check-orientation    -c my_sample.yaml
uv run lightsuite brain init-registration    -c my_sample.yaml
uv run lightsuite brain match-points         -c my_sample.yaml
uv run lightsuite brain register             -c my_sample.yaml
uv run lightsuite brain export               -c my_sample.yaml --save-volume --write-csv
```

| Stage | Status |
|-------|--------|
| Preprocess, init registration, match points, register, export | Implemented |
| Cell detection & atlas cell mapping | Planned |
| Spinal cord & slice modules | MATLAB only (for now) |

**Guides:**

- [Installation (Python)](https://lightsuite.readthedocs.io/en/latest/installation/)
- [Brain lightsheet usage](https://lightsuite.readthedocs.io/en/latest/usage_lightsheet_brain/)
- Example config: [`examples/brain_lightsheet.yaml`](examples/brain_lightsheet.yaml)

**Requirements:** Python 3.11+, [uv](https://docs.astral.sh/uv/), [Elastix 5.1.0](https://github.com/SuperElastix/elastix/releases/tag/5.1.0) on `PATH`, Allen or Perens atlas NIfTIs.

![Example bspline registration](./images/example_bspline.PNG)

---

## MATLAB pipeline (legacy + spinal cord / slices)

The original MATLAB workflows remain available for all three modalities:

1. **Lightsheet brain** — `demos/ls_analyze_lightsheet_volume.m`
2. **Lightsheet spinal cord** — `demos/ls_analyze_spinal_cord.m`
3. **Widefield coronal slices** — `demos/ls_analyze_slice_volume.m`

**MATLAB requirements:** R2022b+, Image Processing and related toolboxes, [matlab_elastix](https://github.com/dimokaramanlis/matlab_elastix), [yamlmatlab](https://github.com/raacampbell/yamlmatlab), Elastix 5.1.0 on `PATH`, atlas files on the MATLAB path. Run `check_lightsuite_installation.m` to verify.

![Example spinal cord registration](./images/example_spinal_cord.PNG)
![Example slice registration](./images/example_slice_registration.png)

---

## Key features

* **Atlas registration** — Allen CCF v3 (2020) for brain; [Fiederling et al. (2021)](https://data.mendeley.com/datasets/4rrggzv5d5/1) for spinal cord
* **Interactive refinement** — manual control-point matching (Napari in Python; MATLAB GUI in legacy pipeline)
* **Deformable registration** — Elastix B-spline with landmark constraints
* **Atlas-space outputs** — registered volumes and regional intensity tables per channel
* **Cell detection** — 3D band-pass + local maxima (MATLAB; Python port planned)

---

## Support

* **Documentation:** [lightsuite.readthedocs.io](https://lightsuite.readthedocs.io/en/latest/)
* **Issues:** [GitHub Issues](https://github.com/dimokaramanlis/LightSuite/issues) — include OS, Python or MATLAB version, and the full error trace

---

## License and citation

LightSuite is distributed under the **GPL-3.0 License**. See [`LICENSE`](LICENSE).

A preprint is in preparation. If you use this software in your research, please link back to this repository.
