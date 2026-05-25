"""Installation and environment checks (replaces check_lightsuite_installation.m)."""

from __future__ import annotations

import os
import platform
import shutil
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path

from rich.console import Console
from rich.table import Table

from lightsuite.config.loader import load_config
from lightsuite.config.models import BrainPipelineConfig
from lightsuite.atlas.registry import resolve_brain_atlas

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
    if not request_gpu:
        return CheckResult("GPU (optional)", True, "Disabled in config (compute.use_gpu=false).", required=False)

    try:
        import cupy  # noqa: F401
    except ImportError:
        return CheckResult(
            "GPU (optional)",
            True,
            "CuPy not installed. Install with: uv sync --extra gpu",
            required=False,
        )

    try:
        import cupy as cp

        device = cp.cuda.Device(0)
        device.use()
        name = cp.cuda.runtime.getDeviceProperties(0)["name"].decode()
        mem_gb = cp.cuda.runtime.memGetInfo()[1] / (1024**3)
        return CheckResult(
            "GPU (optional)",
            True,
            f"{name}, {mem_gb:.1f} GB free",
            required=False,
        )
    except Exception as exc:  # noqa: BLE001
        return CheckResult("GPU (optional)", False, str(exc), required=False)


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
        resolved = resolve_brain_atlas(cfg.atlas.provider.value, cfg.atlas.atlas_dir)
        results.append(
            CheckResult(
                f"Atlas ({resolved.brain_atlas})",
                True,
                f"{resolved.template_path}",
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


def run_doctor(
    config: BrainPipelineConfig | None = None,
    strict: bool = False,
) -> DoctorReport:
    report = DoctorReport()
    report.add(_check_python())
    report.add(_check_elastix_binary("elastix"))
    report.add(_check_elastix_binary("transformix"))
    report.results.extend(_check_brain_atlas(config, strict))
    report.add(_check_spinal_cord_atlas())

    request_gpu = config.compute.use_gpu if config else True
    report.add(_check_gpu(request_gpu))

    scratch = config.sample.scratch if config else None
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
    if config_path:
        config = load_config(config_path)
        console.print(f"[bold]Validating config:[/bold] {config_path}")

    report = run_doctor(config, strict=strict)

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
