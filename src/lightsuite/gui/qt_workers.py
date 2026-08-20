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


def _thread_is_running(thread: Any) -> bool:
    """Return False if the C++ QThread object has already been deleted."""
    try:
        return bool(thread.isRunning())
    except RuntimeError:
        return False


def cancel_background_workers(workers: list[Any], *, timeout_ms: int = 5_000) -> None:
    """Signal workers to stop and wait briefly; safe to call from the main thread.

    Used when a new space switch starts while the previous one is still loading.
    Workers check ``load_state["cancelled"]`` and bail early; we just wait for
    them to finish cleanly rather than force-killing the thread.
    """
    for worker in workers:
        thread = getattr(worker, "_lightsuite_thread", None)
        if thread is None:
            continue
        if _thread_is_running(thread):
            thread.quit()
            try:
                thread.wait(timeout_ms)
            except RuntimeError:
                pass


def wait_background_workers(workers: list[Any], *, timeout_ms: int = 60_000) -> None:
    """Block until background workers finish (e.g. when tearing down a stage)."""
    for worker in workers:
        thread = getattr(worker, "_lightsuite_thread", None)
        if thread is None:
            continue
        if _thread_is_running(thread):
            thread.quit()
            try:
                thread.wait(timeout_ms)
            except RuntimeError:
                pass


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
        _done: bool = False

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
    # Store a Python-side flag so callers can check completion without touching
    # the C++ object (which may have been deleted by deleteLater).
    worker._lightsuite_thread = thread
    worker._lightsuite_done = False

    def _finish() -> None:
        worker._lightsuite_done = True
        thread.quit()

    def _release_worker() -> None:
        worker._lightsuite_thread = None
        try:
            _ACTIVE_WORKERS.remove(worker)
        except ValueError:
            pass

    worker.returned.connect(on_success)
    worker.errored.connect(on_failure)
    worker.returned.connect(_finish)
    worker.errored.connect(_finish)
    # deleteLater schedules C++ deletion; _release_worker clears our reference
    # beforehand so callers see _lightsuite_thread = None and skip isRunning().
    thread.finished.connect(_release_worker)
    thread.finished.connect(worker.deleteLater)
    thread.finished.connect(thread.deleteLater)

    thread.started.connect(worker.run)
    _ACTIVE_WORKERS.append(worker)
    thread.start()
    return worker
