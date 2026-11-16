"""Background worker thread.

Analysis, auto-sizing and the verification suite all run off the GUI thread so
the window stays responsive and the progress bar actually moves.  The worker
takes any callable that accepts a ``progress=callable(fraction, message)``
keyword and emits its return value when it finishes.
"""

from __future__ import annotations

import re
import traceback

from qtpy.QtCore import QThread, Signal

_CLASS_PREFIX = re.compile(r"^[A-Za-z_][A-Za-z0-9_.]*(?:Error|Exception): ")


def plain_message(message):
    """The text of a ``failed`` signal for a message box.

    The worker sends ``"ClassName: text"``.  This drops the class name and
    says in plain words what a singular stiffness matrix means.  The full
    message and the traceback still go to the Log.
    """
    text = _CLASS_PREFIX.sub("", str(message), count=1).strip()
    low = text.lower()
    if "unstable structure" in low:
        return text                      # the solver's own plain sentence
    if "singular" in low:
        return ("The stiffness matrix is singular, so the structure is a "
                "mechanism with these supports: part of it can move "
                "without resistance. Check the supports.")
    if text.startswith("'") and text.endswith("'") and len(text) < 60:
        return f"A required item is missing: {text}."     # a bare KeyError
    return text


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
