# LightSuite

[![Documentation Status](https://readthedocs.org/projects/lightsuite/badge/?version=latest)](https://lightsuite.readthedocs.io/en/latest/)
[![Python](https://img.shields.io/badge/Python-3.11%2B-blue.svg)](https://www.python.org/)
[![MATLAB](https://img.shields.io/badge/MATLAB-R2022b%2B-blue.svg)](https://www.mathworks.com/products/matlab.html)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

**LightSuite** registers large microscopy volumes to standard anatomical atlases and exports atlas-space intensities and cell coordinates. It is designed for 100 GB+ datasets from lightsheet and widefield imaging.

**Documentation:** [lightsuite.readthedocs.io](https://lightsuite.readthedocs.io/en/latest/)

---

## Python pipeline

LightSuite Python is a CLI companion to the MATLAB software: brain lightsheet, spinal cord,
and multiresolution registration with optional Napari GUIs. **No MATLAB license required**
for those workflows. See [Python vs MATLAB](https://lightsuite.readthedocs.io/en/latest/python_vs_matlab/).

### Getting started

```bash
# From a clone (development)
git clone https://github.com/dimokaramanlis/LightSuite.git
cd LightSuite
uv sync --extra dev --extra gui --extra registration

# Or after PyPI release:
# pip install "lightsuite[gui,registration]"

# Copy a template and edit paths
cp examples/brain_lightsheet.yaml my_sample.yaml

# Verify environment (Elastix on PATH, atlas paths, scratch disk)
uv run lightsuite doctor -c my_sample.yaml

# Brain pipeline (one stage at a time; orchestration CLI planned)
uv run lightsuite brain preprocess           -c my_sample.yaml
uv run lightsuite brain check-orientation    -c my_sample.yaml
uv run lightsuite brain init-registration    -c my_sample.yaml
uv run lightsuite brain match-points         -c my_sample.yaml
uv run lightsuite brain register             -c my_sample.yaml
uv run lightsuite brain export               -c my_sample.yaml --save-volume --write-csv
```

Other workflows: [`examples/spinal_cord.yaml`](examples/spinal_cord.yaml) (`lightsuite spinal`),
multires configs under [`examples/config/multiresolution/`](examples/config/multiresolution/).
Index: [`examples/README.md`](examples/README.md).

| Stage | Status |
|-------|--------|
| Preprocess, init registration, match points, register, export | Implemented |
| Spinal cord registration & per-sample region-stats | Implemented (`lightsuite spinal`) |
| Cell detection, slice module, CZI reader, spinal cohort, GPU detection | **Not in Python** — see [docs](https://lightsuite.readthedocs.io/en/latest/python_vs_matlab/) |

**Guides:**

- [Installation (Python)](https://lightsuite.readthedocs.io/en/latest/installation/)
- [Brain lightsheet usage](https://lightsuite.readthedocs.io/en/latest/usage_lightsheet_brain/)
- Example config: [`examples/brain_lightsheet.yaml`](examples/brain_lightsheet.yaml)

**Requirements:** Python 3.11+, [Elastix 5.1.0](https://github.com/SuperElastix/elastix/releases/tag/5.1.0) on `PATH`, Allen or Perens atlas NIfTIs. [uv](https://docs.astral.sh/uv/) is recommended for development; `pip install` works for end users.

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
* **Atlas-space outputs** — registered volumes (TIFF) and regional intensity tables per channel
* **External annotations** — import native `points.csv` / mask TIFF after registration

For MATLAB-only features (built-in cell detection, slice registration, cohort analysis), see the [Python vs MATLAB guide](https://lightsuite.readthedocs.io/en/latest/python_vs_matlab/).

---

## Support

* **Documentation:** [lightsuite.readthedocs.io](https://lightsuite.readthedocs.io/en/latest/)
* **Issues:** [GitHub Issues](https://github.com/dimokaramanlis/LightSuite/issues) — include OS, Python or MATLAB version, and the full error trace

---

## License and citation

LightSuite is distributed under the **GPL-3.0 License**. See [`LICENSE`](LICENSE).

A preprint is in preparation. If you use this software in your research, please link back to this repository.
