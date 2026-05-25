# Installation

LightSuite Python is installed from the repository with [`uv`](https://docs.astral.sh/uv/). You also need **Elastix 5.1.0** on your system `PATH` and atlas NIfTI files for registration.

## Requirements

| Component | Version / notes |
|-----------|-----------------|
| Python | 3.11 – 3.13 |
| [uv](https://docs.astral.sh/uv/getting-started/installation/) | Package and environment manager |
| [Elastix](https://github.com/SuperElastix/elastix/releases/tag/5.1.0) | 5.1.0 — `elastix` and `transformix` on `PATH` |
| Atlas NIfTIs | Allen CCF (10 µm) or Perens LSFM atlas |
| Scratch disk | ≥50 GB free (500 GB+ recommended for large brains) |

MATLAB is **not** required for the Python brain pipeline.

---

## 1. Clone the repository

```bash
git clone https://github.com/dimokaramanlis/LightSuite.git
cd LightSuite
git checkout feature/python-migration   # Python pipeline branch
```

---

## 2. Install Python dependencies

From the repository root:

```bash
# Core CLI + tests
uv sync --extra dev

# Napari GUI for match-points (recommended)
uv sync --extra gui

# All optional extras
uv sync --extra dev --extra gui --extra formats --extra gpu
```

This creates a virtual environment in `.venv` and installs the `lightsuite` command.

Verify:

```bash
uv run lightsuite --help
uv run lightsuite brain --help
```

---

## 3. Install Elastix 5.1.0

LightSuite calls `elastix` and `transformix` as subprocesses (same as the MATLAB pipeline).

### Download

Get the [5.1.0 release](https://github.com/SuperElastix/elastix/releases/tag/5.1.0) for your platform and unpack to a permanent location (e.g. `/opt/elastix` or `C:\elastix`).

### Add to PATH

**macOS / Linux** — add to `~/.zshrc` or `~/.bashrc`:

```bash
export PATH="/path/to/elastix/bin:$PATH"
```

**Windows** — add the folder containing `elastix.exe` and `transformix.exe` to the system PATH.

### Verify

Open a **new** terminal:

```bash
elastix --help
transformix --help
```

---

## 4. Download atlas data

### Allen brain atlas (default)

Registration expects these files in one directory:

- `average_template_10.nii.gz`
- `annotation_10.nii.gz`

Download from the [Allen Brain Cell Atlas](https://alleninstitute.github.io/abc_atlas_access/descriptions/Allen-CCF-2020.html) ([volume files](https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com/index.html#image_volumes/Allen-CCF-2020/20230630/)).

**Optional (for parcellation CSV export):**

- `parcellation_to_parcellation_term_membership.csv` — from the Allen SDK / metadata bundle

Point LightSuite at the atlas in either way:

```yaml
# In your config YAML
atlas:
  provider: allen
  resolution_um: 10
  atlas_dir: /path/to/allen_ccf_10um
```

or:

```bash
export LIGHTSUITE_ATLAS_PATH=/path/to/allen_ccf_10um
```

For Allen parcellation tables:

```bash
export LIGHTSUITE_ALLEN_PARCELLATION_CSV=/path/to/parcellation_to_parcellation_term_membership.csv
```

### Perens LSFM atlas

Set `atlas.provider: perens` and `resolution_um: 20`. Place `gubra_template_olf.nii.gz`, `gubra_ano_olf.nii.gz`, and optionally `ARA2_annotation_info_avail_regions.csv` in `atlas_dir`.

---

## 5. Verify installation

```bash
uv run lightsuite doctor
```

With a config file (checks scratch path, atlas, etc.):

```bash
uv run lightsuite doctor -c examples/brain_lightsheet.yaml
```

Use `--strict` to fail on optional atlas warnings.

Expected checks:

| Check | Required |
|-------|----------|
| Python ≥3.11 | Yes |
| `elastix` / `transformix` on PATH | Yes |
| Atlas NIfTIs | Yes (for registration) |
| GPU (CuPy) | Optional |
| Scratch free space | Warning if low |

---

## Optional extras

| Extra | Install | Purpose |
|-------|---------|---------|
| `gui` | `uv sync --extra gui` | Napari match-points GUI |
| `dev` | `uv sync --extra dev` | pytest, ruff |
| `formats` | `uv sync --extra formats` | Future OME-Zarr / Imaris readers |
| `gpu` | `uv sync --extra gpu` | CuPy (future GPU detection) |

---

## Troubleshooting

### `elastix not found`

Elastix is not on `PATH` in the shell where you run `lightsuite`. Confirm with `which elastix` (macOS/Linux) or `where elastix` (Windows).

### Atlas not found

Set `atlas.atlas_dir` in YAML or `LIGHTSUITE_ATLAS_PATH`. Both template and annotation NIfTIs must exist in that folder.

### `uv: command not found`

Install uv: `curl -LsSf https://astral.sh/uv/install.sh | sh`, then restart the terminal.

### Napari fails to launch

Install GUI extras: `uv sync --extra gui`. On Linux you may need a working Qt platform plugin (`QT_QPA_PLATFORM`).

### Tests

```bash
uv sync --extra dev
uv run pytest -v
```

Integration tests that call Elastix are skipped automatically when `elastix` is not on `PATH`.
