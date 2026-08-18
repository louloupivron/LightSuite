# MATLAB ↔ Python interoperability

Optional helpers for mixed MATLAB/Python registration workflows. The main
repository no longer ships the legacy MATLAB toolbox; these scripts and CLI
commands let you exchange checkpoints when you still run MATLAB elsewhere.

## Control points (match-points session)

| Direction | Tool |
|-----------|------|
| MATLAB → Python | `lightsuite brain import-matlab-control-points --mat atlas2histology_tform.mat` |
| Python → MATLAB | `lightsuite brain export-matlab-control-points --save-path <save_path>` |

After exporting from Python, verify in MATLAB before `register`:

```matlab
addpath('devtools/matlab_interop');
import_python_control_points_for_matlab('<save_path>');
```

Implementation: `src/lightsuite/import_/matlab_control_points.py`

## Registration options (`regopts`)

MATLAB `affinetform3d` objects are not readable in Python. On the MATLAB side,
export numeric 4×4 matrices once per checkpoint:

```matlab
addpath('devtools/matlab_interop');
export_regopts_for_python('<save_path>');
```

Then in Python (library API or parity scripts under `devtools/parity/`):

```python
from pathlib import Path
from lightsuite.import_.matlab_regopts import import_matlab_regopts

import_matlab_regopts(Path("<save_path>"))
```

Writes `regopts.json` beside the MATLAB `regopts.mat`.

## See also

- [`devtools/parity/`](../parity/README.md) — ad-hoc MATLAB/Python parity scratch scripts
