"""Background Qt workers for the unified GUI shell.

Avoids ``napari.qt`` imports so auto stages still run when napari's public Qt
shim is missing or broken in the active environment.
"""

from __future__ import annotations

from collections.abc import Callable
from typing import Any, TypeVar

T = TypeVar("T")


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

    def _finish() -> None:
        thread.quit()

    worker.returned.connect(on_success)
    worker.errored.connect(on_failure)
    worker.returned.connect(_finish)
    worker.errored.connect(_finish)
    thread.finished.connect(worker.deleteLater)
    thread.finished.connect(thread.deleteLater)

    thread.started.connect(worker.run)
    thread.start()
    worker._lightsuite_thread = thread  # keep thread alive until finished
    return worker
