"""Background worker thread.

Analysis, auto-sizing and the verification suite all run off the GUI thread so
the window stays responsive and the progress bar actually moves.  The worker
takes any callable that accepts a ``progress=callable(fraction, message)``
keyword and emits its return value when it finishes.
"""

from __future__ import annotations

import traceback

from qtpy.QtCore import QThread, Signal


class Worker(QThread):
    progressed = Signal(float, str)
    succeeded = Signal(object)
    failed = Signal(str, str)

    def __init__(self, func, *args, parent=None, **kwargs):
        super().__init__(parent)
        self._func = func
        self._args = args
        self._kwargs = kwargs

    def _progress(self, fraction, message=""):
        self.progressed.emit(float(fraction), str(message))

    def run(self):                                   # pragma: no cover (thread)
        try:
            result = self._func(*self._args, progress=self._progress,
                                **self._kwargs)
        except Exception as exc:
            self.failed.emit(f"{type(exc).__name__}: {exc}",
                             traceback.format_exc())
            return
        self.succeeded.emit(result)
