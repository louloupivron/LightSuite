# Spinal cord lightsheet analysis

The **Python pipeline does not yet include spinal cord registration.** Use the MATLAB workflow described below until this module is ported.

For the Python brain pipeline, see [Brain lightsheet analysis](usage_lightsheet_brain.md).

---

## MATLAB workflow (current)

Spinal cord volumes register against the [Fiederling et al. (2021)](https://www.sciencedirect.com/science/article/pii/S2667237521001260) atlas. The pipeline straightens the cord, defines anatomical axes, and runs a registration workflow similar to the brain module.

### Getting started

1. Install LightSuite for MATLAB ([Installation](installation.md) — see legacy MATLAB section in repository `README.md`)
2. Download the Fiederling atlas from [Mendeley Data](https://data.mendeley.com/datasets/4rrggzv5d5/1)
3. Run `demos/ls_analyze_spinal_cord.m`

### Key differences from brain

- Dedicated GUI for central canal and anterior/posterior axis definition
- Cord-specific atlas geometry and control-point tools
- Entry script: `ls_analyze_spinal_cord.m`

### Python status

Spinal cord support is planned after the brain pipeline stabilizes. Track progress on branch `feature/python-migration` or open a GitHub issue if you need this prioritized.
