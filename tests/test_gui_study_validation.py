"""
GUI tests: the parameter study (geodome/gui/study_view.py) and the
Validation menu (geodome/gui/validation_view.py), driven through the real
main window off screen.

What is checked
---------------
* install() adds the menu entries, the "Study" view tab and the "Study" and
  "Validation" result tabs (or the window already has them);
* the study dialog rejects bad input with a short message;
* a two-frequency study runs in the worker, fills the table and the four
  charts, exports CSV, and "Open selected row as model" makes that dome the
  model (geometry and sections);
* Cancel stops a study before the next point;
* the hand-calculation validation runs both scripts' compare functions,
  every comparison passes, and nothing is written into validation/;
* with OpenSees missing (monkeypatched) the cross-check action is disabled
  with a tooltip and a click does not start a run; with OpenSees present the
  quick cross-check passes and writes nothing into validation/;
* Export OpenSeesPy script writes static, modal and nonlinear scripts;
* a missing validation folder disables the actions with a tooltip;
* no exception escapes a Qt slot.

Run with pytest::

    pytest -q tests/test_gui_study_validation.py

or standalone::

    python tests/test_gui_study_validation.py
"""

from __future__ import annotations

import builtins
import contextlib
import csv
import os
import sys
import tempfile
import time
import traceback
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ.setdefault("QT_QPA_FONTDIR", "C:/Windows/Fonts")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import matplotlib                                                   # noqa: E402
matplotlib.use("Agg")


# ================================================================ helpers ==
def _app():
    from qtpy import QtWidgets
    return QtWidgets.QApplication.instance() or QtWidgets.QApplication([])


def _window():
    """A main window with both extensions (installed once)."""
    app = _app()
    from geodome.gui.main_window import MainWindow
    from geodome.gui import study_view, validation_view
    win = MainWindow()
    win.resize(1600, 1000)
    study = getattr(win, "study_extension", None) or study_view.install(win)
    val = (getattr(win, "validation_extension", None)
           or validation_view.install(win))
    app.processEvents()
    return app, win, study, val


def _wait(app, win, timeout=600.0):
    """Process events until the window is no longer busy."""
    t0 = time.time()
    while win._busy:
        app.processEvents()
        time.sleep(0.02)
        if time.time() - t0 > timeout:
            raise AssertionError("the worker did not finish in time")
    for _ in range(5):
        app.processEvents()


class _SlotErrors(list):
    """Collect exceptions raised inside Qt slots (PyQt would abort)."""

    def __enter__(self):
        self._old = sys.excepthook
        sys.excepthook = lambda t, v, tb: self.append(
            "".join(traceback.format_exception(t, v, tb)))
        return self

    def __exit__(self, *exc):
        sys.excepthook = self._old
        return False


class _NoDialogs(list):
    """Replace every modal message box and file dialog; record the calls."""

    def __init__(self, save_path=""):
        super().__init__()
        self.save_path = save_path

    def __enter__(self):
        from qtpy import QtWidgets
        mb, fd = QtWidgets.QMessageBox, QtWidgets.QFileDialog
        self._saved = [(mb, n, getattr(mb, n)) for n in
                       ("information", "warning", "critical", "question")]
        self._saved += [(fd, "getSaveFileName", fd.getSaveFileName)]
        for n in ("information", "warning", "critical"):
            setattr(mb, n, staticmethod(
                lambda *a, _n=n, **k: self.append((_n, a[1:3])) or mb.Ok))
        mb.question = staticmethod(
            lambda *a, **k: self.append(("question", a[1:3])) or mb.Yes)
        fd.getSaveFileName = staticmethod(
            lambda *a, **k: (self.save_path, ""))
        return self

    def __exit__(self, *exc):
        for owner, name, fn in self._saved:
            setattr(owner, name, fn)
        return False


@contextlib.contextmanager
def _watch_writes(folder):
    """Record every file opened for writing inside ``folder``."""
    folder = os.path.normcase(os.path.abspath(folder))
    hits = []
    real_open = builtins.open

    def spy(file, mode="r", *a, **k):
        if isinstance(file, (str, bytes, os.PathLike)) and any(
                c in str(mode) for c in "wax+"):
            p = os.path.normcase(os.path.abspath(os.fsdecode(file)))
            if p.startswith(folder):
                hits.append(p)
        return real_open(file, mode, *a, **k)

    builtins.open = spy
    try:
        yield hits
    finally:
        builtins.open = real_open


def _menu_texts(menu):
    return [a.text().replace("&", "") for a in menu.actions()
            if not a.isSeparator()]


def _top_menu(win, title):
    for a in win.menuBar().actions():
        if a.menu() is not None and \
                a.text().replace("&", "").strip() == title:
            return a.menu()
    return None


# ================================================================== tests ==
def test_install_adds_menus_and_tabs():
    app, win, study, val = _window()
    with _SlotErrors() as errors:
        views = [win.view_tabs.tabText(i) for i in range(win.view_tabs.count())]
        results = [win.result_tabs.tabText(i)
                   for i in range(win.result_tabs.count())]
        assert views.count("Study") == 1, views
        assert results.count("Study") == 1 and \
            results.count("Validation") == 1, results
        analysis = _top_menu(win, "Analysis")
        assert analysis is not None
        assert "Parameter study..." in _menu_texts(analysis)
        titles = [a.text().replace("&", "") for a in win.menuBar().actions()]
        assert "Validation" in titles
        assert titles.index("Validation") < titles.index("Help"), titles
        texts = _menu_texts(_top_menu(win, "Validation"))
        for t in ("Run verification suite (75 checks)",
                  "Run hand-calculation validation",
                  "Run OpenSees cross-check",
                  "Export OpenSeesPy script...",
                  "Open validation report (Markdown) folder"):
            assert t in texts, (t, texts)
        assert study.study_rows is None
        assert hasattr(study, "on_model_changed") and \
            hasattr(study, "on_results")
        assert hasattr(val, "on_model_changed") and hasattr(val, "on_results")
        study.on_model_changed(win.model)
        study.on_results(None)
        val.on_model_changed(win.model)
        val.on_results(None)
        assert not errors, errors[-1]
    win.close()


def test_study_inputs_are_checked():
    from geodome.gui import study_view as sv
    assert sv.parse_values("2, 4;6  8", "f", int) == [2, 4, 6, 8]
    assert sv.parse_values("2,2,3", "f", int) == [2, 3]
    assert sv.check_study_inputs("2,4", "8 12", "0.5") == \
        ([2, 4], [8.0, 12.0], [0.5])
    for bad in (("", "12", "0.5"), ("2, x", "12", "0.5"),
                ("2.5", "12", "0.5"), ("0", "12", "0.5"),
                ("11", "12", "0.5"), ("2", "-3", "0.5"),
                ("2", "12", "1.5"), ("2", "12", "0"), ("2", "nan", "0.5")):
        try:
            sv.check_study_inputs(*bad)
        except ValueError as exc:
            assert str(exc) and len(str(exc)) < 90, exc
        else:
            raise AssertionError(f"accepted {bad}")
    app, win, study, val = _window()
    dlg = sv.StudyDialog(win.model, win)
    assert dlg.run_button.isEnabled()
    dlg.frequencies.setText("2, 4, banana")
    assert not dlg.run_button.isEnabled()
    assert "banana" in dlg.error_label.text()
    dlg.frequencies.setText("2, 3")
    dlg.radii.setText("8, 10")
    assert dlg.run_button.isEnabled()
    assert dlg.points_label.text().startswith("4 domes")
    v = dlg.values()
    assert v["frequencies"] == [2, 3] and v["radii"] == [8.0, 10.0]
    assert v["do_buckling"] and not v["auto_size"]
    dlg.close()
    win.close()


def test_study_two_frequencies():
    """Run, table, charts, CSV export, x axis, open a row as the model."""
    app, win, study, val = _window()
    tmp = tempfile.mkdtemp(prefix="geodome_study_")
    csv_path = os.path.join(tmp, "study.csv")
    with _SlotErrors() as errors, _NoDialogs(save_path=csv_path) as boxes:
        assert study.run_study([2, 3], radii=[8.0], fractions=[0.5])
        _wait(app, win)
        assert not errors, errors[-1]
        rows = study.study_rows
        assert rows is not None and len(rows) == 2
        assert [r["frequency"] for r in rows] == [2, 3]
        assert all(r["status"] != "ERROR" for r in rows)
        # table
        assert study.table.rowCount() == 2
        assert study.table.item(0, 0).text() == "2"
        assert "2 domes analysed" in study.summary_label.text()
        # four charts with data
        axes = study.canvas.figure.axes
        assert len(axes) == 4
        for ax in axes:
            assert any([float(v) for v in ln.get_xdata()] == [2.0, 3.0]
                       for ln in ax.lines), ax.get_title()
        assert [t.get_text() for t in axes[0].get_xticklabels()] == \
            ["2", "3"]
        assert win.view_tabs.currentWidget() is study.view_tab
        assert win.result_tabs.currentWidget() is study.result_tab
        # other x axis
        study.x_combo.setCurrentIndex(study.x_combo.findData("nodes"))
        app.processEvents()
        assert study.canvas.figure.axes[0].get_xlabel() == "Nodes"
        # CSV export (through the file dialog, then directly)
        assert study.export_csv() == csv_path
        with open(csv_path, encoding="utf-8") as fh:
            data = list(csv.DictReader(fh))
        assert len(data) == 2 and data[1]["frequency"] == "3"
        assert float(data[0]["mass_kg"]) == rows[0]["mass_kg"]
        # open the selected row as the model
        study.table.selectRow(1)
        app.processEvents()
        assert study.open_button.isEnabled()
        model = study.open_row_as_model()
        assert model is not None and win.model is model
        assert any(b[0] == "question" for b in boxes)
        assert win.model.geom.frequency == 3
        assert abs(win.model.geom.radius - 8.0) < 1e-12
        assert list(win.model.group_sections) == rows[1]["group_sections"]
        assert win.results is None
        assert win.geometry_panel.frequency.value() == 3
        assert not errors, errors[-1]
    win.close()


def test_study_cancel():
    app, win, study, val = _window()
    with _SlotErrors() as errors, _NoDialogs():
        assert study.run_study([2, 3, 4, 5], radii=[8.0])
        assert not study.cancel_button.isHidden()
        study.cancel()
        _wait(app, win)
        rows = study.study_rows
        assert rows is not None and len(rows) < 4, len(rows)
        assert any(m.startswith("Stopped by the user") for m in rows.messages)
        assert study.table.rowCount() == len(rows)
        assert len(study.canvas.figure.axes) == 4
        assert study.cancel_button.isHidden()
        assert not errors, errors[-1]
    win.close()


def test_study_refuses_custom_supports():
    import copy
    from geodome import geometry as geo
    app, win, study, val = _window()
    base = copy.deepcopy(win.model)
    base.support_type = geo.CUSTOM
    with _SlotErrors() as errors, _NoDialogs() as boxes:
        assert not study.run_study([2, 3], base_model=base)
        assert not win._busy
        assert boxes and boxes[-1][0] == "information"
        assert not errors, errors[-1]
    win.close()


def test_hand_validation_writes_nothing():
    from geodome.gui import validation_view as vv
    app, win, study, val = _window()
    if not val.act_hand.isEnabled():
        print("    (skipped: the validation folder is not here)")
        win.close()
        return
    with _SlotErrors() as errors, _NoDialogs() as boxes, \
            _watch_writes(vv.VALIDATION_DIR) as writes:
        assert val.run_hand_validation()
        _wait(app, win)
        assert not errors, errors[-1]
    assert not writes, writes
    assert not boxes, boxes                       # no failure warning
    res = val.hand_result
    assert res["n"] > 400 and res["n_pass"] == res["n"], \
        [r for r in res["rows"] if not r["pass"]][:5]
    assert res["ab"][1] == 80 and res["ab"][0] == 80
    assert not res["info"]["missing"]
    assert val.hand_table.rowCount() == res["n"]
    assert f"{res['n']} of {res['n']}" in val.hand_label.text()
    headers = [val.hand_table.horizontalHeaderItem(j).text()
               for j in range(val.hand_table.columnCount())]
    for h in ("Case", "Quantity", "Hand", "Program", "Difference", "Pass"):
        assert h in headers, headers
    assert win.result_tabs.currentWidget() is val.result_tab
    val.fail_only.setChecked(True)
    assert val.hand_table.rowCount() == 0
    win.close()


def test_opensees_missing_disables_the_action():
    import geodome.opensees_bridge as ob
    app, win, study, val = _window()
    saved = ob.opensees_status
    ob.opensees_status = lambda: {"mode": None, "python": None}
    try:
        with _SlotErrors() as errors, _NoDialogs() as boxes:
            val.menu.aboutToShow.emit()           # what opening the menu does
            assert not val.act_opensees.isEnabled()
            tip = val.act_opensees.toolTip()
            assert "OpenSees was not found" in tip and ob.ENV_VAR in tip
            assert not val.run_opensees_crosscheck()
            assert not win._busy
            assert boxes and boxes[-1][0] == "information"
            # the other actions do not depend on OpenSees
            assert val.act_export.isEnabled() and val.act_verify.isEnabled()
            assert not errors, errors[-1]
    finally:
        ob.opensees_status = saved
    val.refresh_actions()
    win.close()


def test_opensees_quick_crosscheck():
    import geodome.opensees_bridge as ob
    from geodome.gui import validation_view as vv
    if not ob.opensees_available() or \
            vv.missing_files(vv.VALIDATION_DIR, (vv.OPENSEES_SCRIPT,)):
        print("    (skipped: no OpenSees interpreter or no validation folder)")
        return
    app, win, study, val = _window()
    with _SlotErrors() as errors, _NoDialogs() as boxes, \
            _watch_writes(vv.VALIDATION_DIR) as writes:
        assert val.act_opensees.isEnabled()
        assert val.run_opensees_crosscheck()
        _wait(app, win, timeout=900)
        assert not errors, errors[-1]
    assert not writes, writes
    res = val.opensees_result
    assert res["n"] > 50 and res["n_pass"] == res["n"], \
        [r for r in res["rows"] if not r["pass"]][:5]
    assert not boxes, boxes
    assert {s["section"].split()[0] for s in res["summary"]} == \
        {"P", "S", "M", "N", "X"}
    assert val.os_summary_table.rowCount() == len(res["summary"])
    assert val.os_table.rowCount() == res["n"]
    val.os_summary_table.selectRow(0)                 # filter one section
    app.processEvents()
    n_p = sum(1 for r in res["rows"] if r["section"] == "P")
    assert val.os_table.rowCount() == n_p
    win.close()


def test_export_openseespy_script():
    app, win, study, val = _window()
    tmp = tempfile.mkdtemp(prefix="geodome_ops_")
    with _SlotErrors() as errors, _NoDialogs() as boxes:
        for name, analysis in (("static.py", None), ("modal.py", "modal"),
                               ("nl.py", "nonlinear")):
            path = os.path.join(tmp, name)
            out = val.export_script(path=path, analysis=analysis)
            assert out == path, boxes
            text = Path(path).read_text(encoding="utf-8")
            assert "ops." in text and "import" in text
            compile(text, path, "exec")                # valid Python
        assert not boxes, boxes
        assert not errors, errors[-1]
    from geodome.gui import validation_view as vv
    assert vv.script_name(win.model, "LC2  1.2D + 1.6S", "static") == \
        f"opensees_{win.model.geom.frequency}V_LC2.py"
    assert vv.script_name(win.model, None, "modal").endswith("_modal.py")
    win.close()


def test_missing_validation_folder():
    app, win, study, val = _window()
    keep = val.folder
    val.folder = Path(tempfile.mkdtemp(prefix="geodome_noval_")) / "nothing"
    try:
        with _SlotErrors() as errors, _NoDialogs() as boxes:
            state = val.refresh_actions()
            assert not state["hand"] and not state["opensees"] \
                and not state["folder"]
            for a in (val.act_hand, val.act_opensees, val.act_folder):
                assert not a.isEnabled()
                assert "validation folder" in a.toolTip(), a.toolTip()
            assert not val.run_hand_validation()
            assert not val.run_opensees_crosscheck()
            assert not val.open_folder()
            assert not win._busy
            assert len(boxes) == 3
            assert val.act_export.isEnabled()
            assert not errors, errors[-1]
    finally:
        val.folder = keep
        val.refresh_actions()
    win.close()


def test_open_folder_and_verify_action():
    from qtpy import QtGui
    app, win, study, val = _window()
    opened, calls = [], []
    saved = QtGui.QDesktopServices.openUrl
    QtGui.QDesktopServices.openUrl = staticmethod(
        lambda url: opened.append(url.toLocalFile()) or True)
    win.on_verify = lambda: calls.append(1)
    try:
        with _SlotErrors() as errors, _NoDialogs():
            if val.act_folder.isEnabled():
                val.act_folder.trigger()
                assert opened and opened[0].replace("\\", "/").rstrip(
                    "/").endswith("validation/output"), opened
            val.act_verify.trigger()
            assert calls == [1]
            assert not errors, errors[-1]
    finally:
        QtGui.QDesktopServices.openUrl = saved
    win.close()


def _run_standalone():
    fns = [(n, f) for n, f in sorted(globals().items())
           if n.startswith("test_") and callable(f)]
    npass = 0
    for name, fn in fns:
        t0 = time.time()
        try:
            fn()
        except Exception as exc:
            print(f"  FAIL  {name}: {type(exc).__name__}: {exc}")
            traceback.print_exc()
        else:
            npass += 1
            print(f"  pass  {name}  ({time.time() - t0:.1f} s)")
    print(f"\n{npass} of {len(fns)} tests passed.")
    return 0 if npass == len(fns) else 1


if __name__ == "__main__":
    raise SystemExit(_run_standalone())
