"""Installation and environment checks (replaces check_lightsuite_installation.m)."""

from __future__ import annotations

import os
import platform
import shutil
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from rich.console import Console
from rich.table import Table

from lightsuite.config.loader import load_config
from lightsuite.config.models import BrainPipelineConfig, CordTiffLayout, SpinalCordPipelineConfig
from lightsuite.atlas.fiederling import resolve_fiederling_paths
from lightsuite.atlas.registry import resolve_brain_atlas, resolve_brain_atlas_from_config
from lightsuite.registration.bcpd import find_bcpd_executable

console = Console()

ELASTIX_EXPECTED_VERSION = "5.1.0"
MIN_PYTHON = (3, 11)
MIN_SCRATCH_GB = 50.0


@dataclass
class CheckResult:
    name: str
    ok: bool
    detail: str
    required: bool = True


@dataclass
class DoctorReport:
    results: list[CheckResult] = field(default_factory=list)

    @property
    def passed(self) -> bool:
        return all(r.ok for r in self.results if r.required)

    def add(self, result: CheckResult) -> None:
        self.results.append(result)


def _check_python() -> CheckResult:
    version = sys.version_info[:3]
    ok = version >= MIN_PYTHON
    detail = f"{version[0]}.{version[1]}.{version[2]}"
    if not ok:
        detail += f" (requires >={MIN_PYTHON[0]}.{MIN_PYTHON[1]})"
    return CheckResult("Python", ok, detail)


def _check_package_install() -> CheckResult:
    """Verify the installed CLI includes the current brain subcommands."""
    import lightsuite
    import lightsuite.cli.main as main_module

    main_path = Path(main_module.__file__).resolve()
    source = main_path.read_text(encoding="utf-8")
    required_commands = ("check-orientation", "export", "match-points")
    missing = [name for name in required_commands if f'"{name}"' not in source and f"'{name}'" not in source]
    ok = not missing
    detail = f"v{lightsuite.__version__} from {main_path}"
    if missing:
        detail += (
            f"; missing commands: {', '.join(missing)}. "
            "On branch feature/python-migration run: "
            "uv sync --extra dev --reinstall-package lightsuite"
        )
    return CheckResult("LightSuite package", ok, detail)


def _check_elastix_binary(name: str) -> CheckResult:
    path = shutil.which(name)
    if path is None:
        return CheckResult(
            f"{name} binary",
            False,
            "Not found on PATH. Install Elastix 5.1.0 and add bin/ to PATH.",
        )

    try:
        proc = subprocess.run(
            [name, "--help"],
            capture_output=True,
            text=True,
            check=False,
            timeout=30,
        )
    except OSError as exc:
        return CheckResult(f"{name} binary", False, str(exc))

    if proc.returncode != 0:
        return CheckResult(f"{name} binary", False, f"Exit code {proc.returncode}")

    first_line = (proc.stdout or proc.stderr).splitlines()[0] if (proc.stdout or proc.stderr) else ""
    ok = ELASTIX_EXPECTED_VERSION in first_line or "elastix" in first_line.lower()
    detail = f"{path} — {first_line.strip() or 'help OK'}"
    if ELASTIX_EXPECTED_VERSION not in first_line:
        detail += f" (expected version {ELASTIX_EXPECTED_VERSION})"
    return CheckResult(f"{name} binary", ok, detail, required=True)


def _check_gpu(request_gpu: bool) -> CheckResult:
    if request_gpu:
        return CheckResult(
            "GPU (optional)",
            True,
            "compute.use_gpu=true is ignored — GPU-accelerated cell detection is not implemented in Python.",
            required=False,
        )
    return CheckResult(
        "GPU (optional)",
        True,
        "GPU-accelerated cell detection is not implemented in Python (MATLAB only).",
        required=False,
    )


def _check_disk(path: Path | None, label: str, min_gb: float) -> CheckResult:
    if path is None:
        return CheckResult(label, True, "Not configured.", required=False)

    target = path.expanduser()
    if not target.exists():
        return CheckResult(label, False, f"Path does not exist: {target}", required=False)

    usage = shutil.disk_usage(target)
    free_gb = usage.free / (1024**3)
    ok = free_gb >= min_gb
    detail = f"{target}: {free_gb:.1f} GB free (recommend >= {min_gb:.0f} GB)"
    return CheckResult(label, ok, detail, required=False)


def _check_brain_atlas(cfg: BrainPipelineConfig | None, strict: bool) -> list[CheckResult]:
    results: list[CheckResult] = []
    if cfg is None:
        for atlas_id in ("allen",):
            try:
                resolve_brain_atlas(atlas_id, atlas_dir=None)
                results.append(CheckResult(f"Atlas ({atlas_id})", True, "Found via LIGHTSUITE_ATLAS_PATH or cwd."))
            except FileNotFoundError as exc:
                results.append(
                    CheckResult(
                        f"Atlas ({atlas_id})",
                        False,
                        str(exc),
                        required=strict,
                    )
                )
        return results

    try:
        resolved = resolve_brain_atlas_from_config(cfg.atlas)
        detail = f"{resolved.template_path}"
        if resolved.atlas_source == "brainglobe":
            detail = f"{resolved.brainglobe_name} @ {resolved.template_path}"
        results.append(
            CheckResult(
                f"Atlas ({resolved.brain_atlas})",
                True,
                detail,
            )
        )
        if resolved.boundary_path is not None:
            results.append(
                CheckResult(
                    "Atlas boundary volume (optional)",
                    True,
                    str(resolved.boundary_path),
                    required=False,
                )
            )
        else:
            results.append(
                CheckResult(
                    "Atlas boundary volume (optional)",
                    False,
                    "annotation_boundary_10.nii.gz not found; init-registration previews "
                    "derive boundaries from annotation labels.",
                    required=False,
                )
            )
    except FileNotFoundError as exc:
        results.append(
            CheckResult(
                f"Atlas ({cfg.atlas.provider.value})",
                False,
                str(exc),
                required=True,
            )
        )
    return results


def _check_spinal_cord_atlas() -> CheckResult:
    search_dirs = _atlas_search_dirs(None)
    for directory in search_dirs:
        candidate = directory / "Segments.csv"
        if candidate.is_file():
            return CheckResult("Spinal cord atlas (optional)", True, str(candidate), required=False)
    return CheckResult(
        "Spinal cord atlas (optional)",
        True,
        "Segments.csv not found (needed only for cord pipeline).",
        required=False,
    )


def _atlas_search_dirs(explicit: Path | None) -> list[Path]:
    dirs: list[Path] = []
    if explicit is not None:
        dirs.append(explicit.expanduser().resolve())
    env = os.environ.get("LIGHTSUITE_ATLAS_PATH", "")
    for part in env.split(os.pathsep):
        if part.strip():
            dirs.append(Path(part.strip()).expanduser().resolve())
    dirs.append(Path.cwd())
    return dirs


def _check_bcpd(
    config: BrainPipelineConfig | None,
    *,
    config_path: Path | str | None = None,
) -> CheckResult:
    explicit = config.registration.bcpd_path if config else None
    search_roots = [Path(config_path).expanduser().resolve().parent] if config_path else None
    found = find_bcpd_executable(explicit, search_roots=search_roots)
    if found is not None:
        return CheckResult("BCPD (coarse registration)", True, str(found), required=False)
    return CheckResult(
        "BCPD (coarse registration)",
        True,
        "Not found — init-registration falls back to Open3D ICP (lower quality).",
        required=False,
    )


def _check_fiederling_atlas(cfg: SpinalCordPipelineConfig | None) -> CheckResult:
    if cfg is None:
        return _check_spinal_cord_atlas()
    try:
        paths = resolve_fiederling_paths(cfg.atlas.atlas_dir)
        return CheckResult(
            "Spinal cord atlas (Fiederling)",
            True,
            str(paths.atlas_dir),
            required=True,
        )
    except FileNotFoundError as exc:
        return CheckResult("Spinal cord atlas (Fiederling)", False, str(exc), required=True)


def _check_spinal_sample_resolution(cfg: SpinalCordPipelineConfig) -> CheckResult:
    """Warn when plane-per-file stacks declare native voxel size equal to registration grid."""
    from lightsuite.io.cord_volume import (
        _sorted_tiff_files,
        normalize_res_um,
        read_plane_tiff,
        resolve_cord_tiff_layout,
    )

    folder = cfg.sample.source.path
    if folder is None:
        return CheckResult("Spinal sample resolution", False, "sample.source.path is unset", required=False)
    channel_folders = cfg.sample.source.channel_roots
    try:
        layout = resolve_cord_tiff_layout(
            folder,
            cfg.sample.source.tiff_type,
            channel_folders=channel_folders,
        )
    except (FileNotFoundError, ValueError) as exc:
        return CheckResult("Spinal sample resolution", False, str(exc), required=False)

    sampleres = normalize_res_um(cfg.sample.voxel_um)
    regres = normalize_res_um([cfg.registration.resolution_um] * 3)
    resfac = sampleres / regres
    if not np.allclose(resfac, 1.0):
        nchan_note = f", {len(channel_folders)} channels" if channel_folders and len(channel_folders) > 1 else ""
        return CheckResult(
            "Spinal sample resolution",
            True,
            f"Native voxel_um {sampleres.tolist()} → registration {regres.tolist()}{nchan_note}",
            required=False,
        )

    if layout != CordTiffLayout.PLANE_PER_FILE:
        return CheckResult(
            "Spinal sample resolution",
            True,
            "voxel_um matches registration grid (small channel-per-file stack).",
            required=False,
        )

    probe = Path(channel_folders[0]) if channel_folders else folder
    files = _sorted_tiff_files(probe)
    if not files:
        return CheckResult("Spinal sample resolution", False, f"No TIFFs in {probe}", required=False)
    ny, nx = read_plane_tiff(files[0]).shape
    native = (ny, nx, len(files))
    native_gb = float(np.prod(native)) * 2.0 / 1e9
    nchan = len(channel_folders) if channel_folders else 1
    chan_note = f", {nchan} channel(s)" if nchan > 1 else ""
    detail = (
        f"Plane-per-file stack ({native[0]}×{native[1]}×{native[2]} px{chan_note}, "
        f"~{native_gb * nchan:.1f} GB) "
        f"but sample.voxel_um equals registration.resolution_um — no downsampling on load. "
        "Set sample.voxel_um to your native microscope voxel size (e.g. [1.8, 1.8, 4])."
    )
    return CheckResult("Spinal sample resolution", False, detail, required=False)


def run_doctor(
    config: BrainPipelineConfig | None = None,
    strict: bool = False,
    spinal_config: SpinalCordPipelineConfig | None = None,
    multires_config: object | None = None,
    *,
    config_path: Path | str | None = None,
) -> DoctorReport:
    report = DoctorReport()
    report.add(_check_python())
    report.add(_check_package_install())
    report.add(_check_elastix_binary("elastix"))
    report.add(_check_elastix_binary("transformix"))
    report.add(_check_bcpd(config, config_path=config_path))
    if spinal_config is not None:
        report.add(_check_fiederling_atlas(spinal_config))
        report.add(_check_spinal_sample_resolution(spinal_config))
    elif multires_config is not None:
        pair_label = multires_config.multires.resolved_pair_label(multires_config.sample.name)
        report.add(
            CheckResult(
                "Multires config",
                True,
                f"pair_label={pair_label}, save_path={multires_config.sample.save_path}",
                required=False,
            )
        )
    else:
        report.results.extend(_check_brain_atlas(config, strict))
        report.add(_check_spinal_cord_atlas())

    active = config or spinal_config or multires_config
    compute = getattr(active, "compute", None) if active else None
    request_gpu = compute.use_gpu if compute is not None else True
    report.add(_check_gpu(request_gpu))

    scratch = active.sample.scratch if active else None
    report.add(_check_disk(scratch, "Scratch disk", MIN_SCRATCH_GB))

    report.add(
        CheckResult(
            "Platform",
            True,
            f"{platform.system()} {platform.machine()}",
            required=False,
        )
    )
    return report


def doctor_command(config_path: str | None = None, strict: bool = False) -> None:
    config: BrainPipelineConfig | None = None
    spinal_config: SpinalCordPipelineConfig | None = None
    multires_config: object | None = None
    if config_path:
        from lightsuite.config.loader import load_config, load_multires_config, load_spinal_config

        try:
            config = load_config(config_path)
            console.print(f"[bold]Validating brain config:[/bold] {config_path}")
        except Exception:
            try:
                spinal_config = load_spinal_config(config_path)
                console.print(f"[bold]Validating spinal cord config:[/bold] {config_path}")
            except Exception:
                multires_config = load_multires_config(config_path)
                console.print(f"[bold]Validating multires config:[/bold] {config_path}")

    report = run_doctor(
        config,
        strict=strict,
        spinal_config=spinal_config,
        multires_config=multires_config,
        config_path=config_path,
    )

    table = Table(title="LightSuite doctor")
    table.add_column("Check")
    table.add_column("Status")
    table.add_column("Detail")

    for result in report.results:
        status = "[green]OK[/green]" if result.ok else "[red]FAIL[/red]"
        if not result.required and not result.ok:
            status = "[yellow]WARN[/yellow]"
        table.add_row(result.name, status, result.detail)

    console.print(table)

    if report.passed:
        console.print("[green]All required checks passed.[/green]")
        raise SystemExit(0)

    console.print("[red]One or more required checks failed.[/red]")
    raise SystemExit(1)
