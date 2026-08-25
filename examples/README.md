# LightSuite example configs

## Start here (copy-paste templates)

| Workflow | Canonical template | Documentation |
|----------|-------------------|---------------|
| Brain lightsheet | [`brain_lightsheet.yaml`](brain_lightsheet.yaml) | [Brain usage guide](../docs/usage_lightsheet_brain.md) |
| Spinal cord | [`spinal_cord.yaml`](spinal_cord.yaml) | [Spinal usage guide](../docs/usage_spinal_cord.md) |
| Multiresolution (overview ↔ ROI) | [`multiresolution.yaml`](multiresolution.yaml) | [Multires usage guide](../docs/usage_multiresolution.md) |

```bash
cp examples/brain_lightsheet.yaml my_mouse.yaml
# or: cp examples/spinal_cord.yaml my_spinal.yaml
# or: cp examples/multiresolution.yaml my_multires.yaml
# edit paths, then:
uv run lightsuite doctor -c my_mouse.yaml
uv run lightsuite brain run -c my_mouse.yaml --through match-points
```

## Real sample configs (`config/`)

The YAML files under [`config/`](config/) are **machine-specific samples** from lab datasets. They contain absolute paths valid only on the developer workstation. Use them as references for field names and workflow structure — **do not copy paths verbatim**.

| Directory | Contents |
|-----------|----------|
| [`config/mesoSPIM/`](config/mesoSPIM/) | Brain lightsheet mesoSPIM acquisitions |
| [`config/spinal_cord/`](config/spinal_cord/) | Spinal cord samples (OP87, OP93, …) |
| [`config/multiresolution/`](config/multiresolution/) | Overview ↔ ROI multires manifests |
| [`config/smartspim/`](config/smartspim/) | SmartSPIM acquisitions |
| [`config/colm/`](config/colm/) | COLM microscope |
| [`config/rat/`](config/rat/) | Rat brain (non-mouse atlas) |

### Combined multires → brain atlas

Some mesoSPIM projects run multires registration first, then brain atlas registration on the overview-aligned data:

1. [`config/multiresolution/marianna_multires.yaml`](config/multiresolution/marianna_multires.yaml) — overview ↔ ROI
2. [`config/mesoSPIM/marianna_perens.yaml`](config/mesoSPIM/marianna_perens.yaml) — brain atlas on converted paths

## Notebooks (`notebooks/`)

Jupyter notebooks for vendor-specific conversion (Arivis, mesoSPIM, SmartSPIM) live in [`notebooks/`](notebooks/). Prefer the CLI manifest builders once available (`lightsuite multires build-manifest`); until then, see [multires usage](../docs/usage_multiresolution.md).
