"""Progress and log reporting for pipeline stages (CLI, GUI, tests)."""

from __future__ import annotations

import contextlib
import contextvars
import io
from collections.abc import Callable, Iterator
from pathlib import Path
from typing import Protocol, runtime_checkable

from rich.console import Console

_current_reporter: contextvars.ContextVar[Reporter | None] = contextvars.ContextVar(
    "lightsuite_reporter",
    default=None,
)


def get_active_reporter() -> Reporter | None:
    """Return the reporter active in the current context, if any."""
    return _current_reporter.get()


@contextlib.contextmanager
def capture_pipeline_output(reporter: Reporter | None) -> Iterator[None]:
    """Route Rich console output and stdout/stderr to *reporter*."""
    if reporter is None:
        yield
        return

    token = _current_reporter.set(reporter)
    original_print = Console.print

    def forwarding_print(self: Console, *args: object, **kwargs: object) -> None:
        original_print(self, *args, **kwargs)
        capture_kwargs = {key: value for key, value in kwargs.items() if key != "file"}
        with self.capture() as captured:
            original_print(self, *args, **capture_kwargs)
        text = captured.get().strip()
        active = _current_reporter.get()
        if text and active is not None:
            active.message(text)

    class _StreamWriter(io.TextIOBase):
        def write(self, s: str) -> int:
            active = _current_reporter.get()
            if active is not None and s:
                for line in s.splitlines():
                    stripped = line.strip()
                    if stripped:
                        active.message(stripped)
            return len(s)

        def flush(self) -> None:
            return None

    Console.print = forwarding_print  # type: ignore[method-assign]
    try:
        with contextlib.redirect_stdout(_StreamWriter()), contextlib.redirect_stderr(_StreamWriter()):
            yield
    finally:
        Console.print = original_print  # type: ignore[method-assign]
        _current_reporter.reset(token)


@runtime_checkable
class Reporter(Protocol):
    """Sink for pipeline progress messages."""

    def message(self, text: str) -> None:
        """Emit a general status line (may include Rich markup)."""

    def stage_start(self, title: str, *, manual: bool = False, checkpoint_hint: str = "") -> None:
        """Announce that a stage is starting."""

    def stage_skip(self, title: str, state: str, detail: str) -> None:
        """Announce that a stage was skipped."""

    def stage_failed(self, stage_id: str, exc: Exception) -> None:
        """Announce that a stage failed."""

    def rerun_hint(self, workflow: str, stage_id: str, config_path: Path) -> None:
        """Suggest a CLI command to re-run a failed stage."""

    def pipeline_complete(self, ran: int) -> None:
        """Announce successful pipeline completion."""

    def pipeline_no_stages(self) -> None:
        """Announce that no stages were executed."""


class ConsoleReporter:
    """Default Rich console reporter used by the CLI."""

    def __init__(self, console: Console | None = None) -> None:
        self._console = console or Console()

    @property
    def console(self) -> Console:
        return self._console

    def message(self, text: str) -> None:
        self._console.print(text)

    def stage_start(self, title: str, *, manual: bool = False, checkpoint_hint: str = "") -> None:
        tag = " (manual GUI)" if manual else ""
        hint = f"  [dim]{checkpoint_hint}[/dim]" if checkpoint_hint else ""
        self._console.print(f"\n[bold cyan]→ {title}[/bold cyan]{tag}{hint}")

    def stage_skip(self, title: str, state: str, detail: str) -> None:
        self._console.print(f"[dim]Skipping {title} ({state}) — {detail}[/dim]")

    def stage_failed(self, stage_id: str, exc: Exception) -> None:
        self._console.print(f"[bold red]Stage failed:[/bold red] {stage_id} — {exc}")

    def rerun_hint(self, workflow: str, stage_id: str, config_path: Path) -> None:
        self._console.print(f"Re-run manually: lightsuite {workflow} {stage_id} -c {config_path}")

    def pipeline_complete(self, ran: int) -> None:
        self._console.print(f"\n[green]Completed {ran} stage(s).[/green]")

    def pipeline_no_stages(self) -> None:
        self._console.print("[yellow]No stages executed (all complete or optional).[/yellow]")


class NullReporter:
    """No-op reporter for tests and headless automation."""

    def message(self, text: str) -> None:
        return None

    def stage_start(self, title: str, *, manual: bool = False, checkpoint_hint: str = "") -> None:
        return None

    def stage_skip(self, title: str, state: str, detail: str) -> None:
        return None

    def stage_failed(self, stage_id: str, exc: Exception) -> None:
        return None

    def rerun_hint(self, workflow: str, stage_id: str, config_path: Path) -> None:
        return None

    def pipeline_complete(self, ran: int) -> None:
        return None

    def pipeline_no_stages(self) -> None:
        return None


class CallbackReporter:
    """Forward reporter events to callables (e.g. Qt signals in a future GUI)."""

    def __init__(
        self,
        *,
        on_message: Callable[..., None] | None = None,
        on_stage_start: Callable[..., None] | None = None,
        on_stage_skip: Callable[..., None] | None = None,
        on_stage_failed: Callable[..., None] | None = None,
        on_rerun_hint: Callable[..., None] | None = None,
        on_pipeline_complete: Callable[..., None] | None = None,
        on_pipeline_no_stages: Callable[..., None] | None = None,
    ) -> None:
        self._on_message = on_message
        self._on_stage_start = on_stage_start
        self._on_stage_skip = on_stage_skip
        self._on_stage_failed = on_stage_failed
        self._on_rerun_hint = on_rerun_hint
        self._on_pipeline_complete = on_pipeline_complete
        self._on_pipeline_no_stages = on_pipeline_no_stages

    def message(self, text: str) -> None:
        if self._on_message is not None:
            self._on_message(text)

    def stage_start(self, title: str, *, manual: bool = False, checkpoint_hint: str = "") -> None:
        if self._on_stage_start is not None:
            self._on_stage_start(title, manual=manual, checkpoint_hint=checkpoint_hint)

    def stage_skip(self, title: str, state: str, detail: str) -> None:
        if self._on_stage_skip is not None:
            self._on_stage_skip(title, state, detail)

    def stage_failed(self, stage_id: str, exc: Exception) -> None:
        if self._on_stage_failed is not None:
            self._on_stage_failed(stage_id, exc)

    def rerun_hint(self, workflow: str, stage_id: str, config_path: Path) -> None:
        if self._on_rerun_hint is not None:
            self._on_rerun_hint(workflow, stage_id, config_path)

    def pipeline_complete(self, ran: int) -> None:
        if self._on_pipeline_complete is not None:
            self._on_pipeline_complete(ran)

    def pipeline_no_stages(self) -> None:
        if self._on_pipeline_no_stages is not None:
            self._on_pipeline_no_stages()
