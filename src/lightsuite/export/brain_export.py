"""Export registered brain volumes (generateRegisteredBrainVolumes.m port)."""

from __future__ import annotations

import json
import shutil
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from rich.console import Console

from lightsuite.atlas.registry import resolve_brain_atlas
from lightsuite.config.models import BrainPipelineConfig
from lightsuite.export.atlas_space import transform_volume_to_atlas
from lightsuite.export.parcellation import (
    compute_allen_parcellation,
    compute_perens_parcellation,
    write_parcellation_csv,
)
from lightsuite.io.tiff_write import save_registration_volume
from lightsuite.preprocess.checkpoint import RegOptsCheckpoint
from lightsuite.registration.brain_register import TransformParamsCheckpoint
from lightsuite.registration.volume import load_registration_volume

console = Console()


@dataclass
class BrainExportResult:
    registered_volumes: dict[int, Path]
    parcellation_paths: dict[int, Path]
    all_medians: np.ndarray | None


def _load_transform_params(save_path: Path) -> TransformParamsCheckpoint:
    json_path = save_path / "transform_params.json"
    mat_hint = save_path / "transform_params.mat"
    if json_path.is_file():
        return TransformParamsCheckpoint.load(json_path)
    msg = f"Missing {json_path}. Run 'lightsuite brain register' first."
    if mat_hint.is_file():
        msg += f" (Found legacy {mat_hint}; re-run register in Python to produce JSON.)"
    raise FileNotFoundError(msg)


def export_registered_brain_volumes(
    config: BrainPipelineConfig,
    *,
    write_csv: bool | None = None,
    save_registered_volume: bool | None = None,
) -> BrainExportResult:
    """Apply transforms and optionally write registered volumes and parcellation CSVs."""
    if shutil.which("transformix") is None:
        msg = "transformix must be on PATH for brain export."
        raise RuntimeError(msg)

    save_path = config.sample.save_path.expanduser()
    regopts_path = save_path / "regopts.json"
    if not regopts_path.is_file():
        msg = f"Missing {regopts_path}."
        raise FileNotFoundError(msg)

    checkpoint = RegOptsCheckpoint.load(regopts_path)
    transform_params = _load_transform_params(save_path)
    permute = transform_params.permute_sample_to_atlas
    spacing_mm = checkpoint.registres_um * 1e-3

    write_csv = config.export.write_cells_csv if write_csv is None else write_csv
    save_vol = config.export.save_registered_volume if save_registered_volume is None else save_registered_volume

    register_path = save_path / "volume_registered"
    if save_vol:
        register_path.mkdir(parents=True, exist_ok=True)

    channel_paths = {
        int(k): Path(v) for k, v in (checkpoint.regvolpaths or {}).items()
    }
    if not channel_paths:
        msg = "regopts.json missing regvolpaths."
        raise RuntimeError(msg)

    atlas = resolve_brain_atlas(transform_params.brain_atlas, config.atlas.atlas_dir)
    n_chans = len(channel_paths)
    atlas_shape = tuple(int(v) for v in transform_params.atlassize)
    straightvol = np.zeros((*atlas_shape, n_chans), dtype=np.uint16)

    console.print("Applying transforms to registration volumes...")
    t0 = time.perf_counter()
    transformix_root = save_path / "transformix_export_temp"
    transformix_root.mkdir(parents=True, exist_ok=True)

    registered_paths: dict[int, Path] = {}
    for ichan, volpath in sorted(channel_paths.items()):
        console.print(f"Registering channel {ichan}: {volpath.name}")
        volume = load_registration_volume(volpath)
        registered = transform_volume_to_atlas(
            volume,
            transform_params,
            permute=permute,
            spacing_mm=spacing_mm,
            temp_dir=transformix_root / f"chan_{ichan:02d}",
        )
        straightvol[:, :, :, ichan - 1] = registered
        if save_vol:
            out_path = register_path / f"chan_{ichan:02d}_registered_atlas.tif"
            save_registration_volume(registered, out_path)
            registered_paths[ichan] = out_path
        console.print(f"Channel {ichan}/{n_chans} done in {time.perf_counter() - t0:.1f}s.")

    parcellation_paths: dict[int, Path] = {}
    all_medians: np.ndarray | None = None

    if write_csv and atlas.supports_parcellation:
        console.print("Calculating parcellation intensities...")
        for ichan in sorted(channel_paths):
            vol = straightvol[:, :, :, ichan - 1]
            try:
                if atlas.brain_atlas == "allen":
                    result = compute_allen_parcellation(
                        vol,
                        atlas,
                        transform_params.atlas_resolution_um,
                    )
                else:
                    result = compute_perens_parcellation(
                        vol,
                        atlas,
                        transform_params.atlas_resolution_um,
                    )
            except FileNotFoundError as exc:
                console.print(f"[yellow]Skipping parcellation for channel {ichan}:[/yellow] {exc}")
                continue

            register_path.mkdir(parents=True, exist_ok=True)
            json_path = register_path / f"chan{ichan:02d}_intensities.json"
            json_path.write_text(
                json.dumps(
                    {
                        "medianoverareas": result.median_over_areas.tolist(),
                        "stdoverareas": result.std_over_areas.tolist(),
                        "areaidx": result.area_ids.tolist(),
                        "volumeoverareas": result.volume_over_areas.tolist(),
                    },
                    indent=2,
                ),
                encoding="utf-8",
            )
            csv_path = register_path / f"chan{ichan:02d}_intensities.csv"
            write_parcellation_csv(csv_path, result)
            parcellation_paths[ichan] = csv_path
            if all_medians is None:
                n_areas = result.median_over_areas.shape[0]
                all_medians = np.full((n_areas, 2, n_chans), np.nan, dtype=np.float32)
            all_medians[:, :, ichan - 1] = result.median_over_areas
    elif write_csv:
        console.print(
            "[yellow]Parcellation CSV skipped:[/yellow] atlas does not expose structure metadata."
        )

    console.print(f"[green]Export complete.[/green] Output under {save_path}")
    return BrainExportResult(
        registered_volumes=registered_paths,
        parcellation_paths=parcellation_paths,
        all_medians=all_medians,
    )
