"""Background Qt workers for the unified GUI shell.

Avoids ``napari.qt`` imports so auto stages still run when napari's public Qt
shim is missing or broken in the active environment.
"""

from __future__ import annotations

from collections.abc import Callable
from typing import Any, TypeVar

T = TypeVar("T")

# Keep worker QObject alive until its QThread finishes (prevents premature GC).
_ACTIVE_WORKERS: list[Any] = []


def wait_background_workers(workers: list[Any], *, timeout_ms: int = 60_000) -> None:
    """Block until background workers finish (e.g. when tearing down a stage)."""
    for worker in workers:
        thread = getattr(worker, "_lightsuite_thread", None)
        if thread is None:
            continue
        if thread.isRunning():
            thread.quit()
            thread.wait(timeout_ms)


def start_background_task(
    func: Callable[[], T],
    *,
    on_success: Callable[[T], None],
    on_failure: Callable[[BaseException], None],
) -> Any:
    """Run ``func`` on a worker thread; invoke callbacks on the Qt main thread."""
    from qtpy.QtCore import QObject, QThread, Signal

    class _Worker(QObject):
        returned = Signal(object)
        errored = Signal(object)

        def __init__(self, fn: Callable[[], T]) -> None:
            super().__init__()
            self._fn = fn

        def run(self) -> None:
            try:
                self.returned.emit(self._fn())
            except BaseException as exc:
                self.errored.emit(exc)

    thread = QThread()
    worker = _Worker(func)
    worker.moveToThread(thread)
    worker._lightsuite_thread = thread

    def _finish() -> None:
        thread.quit()

    def _release_worker() -> None:
        try:
            _ACTIVE_WORKERS.remove(worker)
        except ValueError:
            pass

    worker.returned.connect(on_success)
    worker.errored.connect(on_failure)
    worker.returned.connect(_finish)
    worker.errored.connect(_finish)
    thread.finished.connect(worker.deleteLater)
    thread.finished.connect(thread.deleteLater)
    thread.finished.connect(_release_worker)

    thread.started.connect(worker.run)
    _ACTIVE_WORKERS.append(worker)
    thread.start()
    return worker
