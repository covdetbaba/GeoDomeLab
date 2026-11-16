"""
GUI test: the Generate dome button, clicked through the real main window.

Added 27 October 2026 after Generate dome failed silently when the new dome
had more strut families than the one on screen (for example 1V -> 6V).  The
window runs off screen, so no display is needed.

Run with pytest::

    pytest -q tests/test_gui_generate.py

or standalone::

    python tests/test_gui_generate.py
"""

from __future__ import annotations

import os
import sys
import traceback

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ.setdefault("QT_QPA_FONTDIR", "C:/Windows/Fonts")
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import matplotlib                                                   # noqa: E402
matplotlib.use("Agg")

# Frequencies in an order that goes up and down in strut families
# (1V has 1 family, 4V 6, 6V 9, 10V 27).
SEQUENCE = [(1, 6.0), (6, 15.0), (2, 6.0), (4, 15.0), (10, 40.0), (3, 9.0),
            (8, 25.0), (1, 6.0), (7, 20.0), (5, 12.0), (9, 30.0), (6, 30.0)]


def _window():
    from qtpy import QtWidgets
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    from geodome.gui.main_window import MainWindow
    return app, MainWindow()


class _SlotErrors(list):
    """Collect exceptions raised inside Qt slots.

    PyQt aborts the whole process on an uncaught exception in a slot, so the
    test installs its own hook and turns those exceptions into failures.
    """

    def __enter__(self):
        self._old = sys.excepthook
        sys.excepthook = lambda t, v, tb: self.append(
            "".join(traceback.format_exception(t, v, tb)))
        return self

    def __exit__(self, *exc):
        sys.excepthook = self._old
        return False


def test_generate_follows_the_inputs():
    """Every click builds the requested dome, whatever was on screen before."""
    app, win = _window()
    with _SlotErrors() as errors:
        gp = win.geometry_panel
        for f, r in SEQUENCE:
            gp.frequency.setValue(f)
            gp.radius.setValue(r)
            gp.fraction.setValue(0.5)
            gp.generate_btn.click()
            app.processEvents()
            g = win.model.geom
            assert not errors, errors[-1]
            assert g.frequency == f and abs(g.radius - r) < 1e-9, (f, r)
            n = g.n_strut_types()
            assert len(win.model.group_sections) == n
            assert len(win.section_panel._combos) == n
            assert win.status_label.text().startswith(f"Generated {f}V")
    win.close()


def test_section_choice_is_kept_on_regenerate():
    """A tube picked for family A survives a change of frequency."""
    app, win = _window()
    with _SlotErrors() as errors:
        gp = win.geometry_panel
        combo = win.section_panel._combos[0]
        combo.setCurrentIndex(3)
        chosen = combo.currentText()
        gp.frequency.setValue(6)
        gp.radius.setValue(15.0)
        gp.generate_btn.click()
        app.processEvents()
        assert not errors, errors[-1]
        assert win.model.group_sections[0] == chosen
    win.close()


def _run_standalone():
    fns = [(n, f) for n, f in sorted(globals().items())
           if n.startswith("test_") and callable(f)]
    npass = 0
    for name, fn in fns:
        try:
            fn()
        except Exception as exc:
            print(f"  FAIL  {name}: {type(exc).__name__}: {exc}")
        else:
            npass += 1
            print(f"  pass  {name}")
    print(f"\n{npass} of {len(fns)} tests passed.")
    return 0 if npass == len(fns) else 1


if __name__ == "__main__":
    raise SystemExit(_run_standalone())
