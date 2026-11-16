"""
GUI tests of the integrated main window (SPEC_GUI H2 item 9).

The window is built once per test with every add-on loaded by
``MainWindow.__init__`` (nonlinear, modal, fabrication, parameter study,
validation), exactly as a user gets it.  The tests click through the real
window off screen:

* the layout: menus, toolbar, tabs, no duplicate entries, Help > Theory;
* every example in the Examples menu loads;
* the linear analysis of the 6V case study, and what the add-ons do with it;
* the nonlinear path of the von Mises truss reaches the closed-form limit
  point, and the busy state turns the other long runs off meanwhile;
* the modal analysis of validation Model A against hand calculation F7;
* the fabrication tables of a 2V dome (chord factors 0.54653 / 0.61803);
* a two-frequency parameter study;
* the Validation menu when OpenSees is missing (monkeypatched);
* the HTML report with every extra section, from File > Export HTML report
  and from File > Export everything;
* the DXF export;
* Save / Open round trip;
* no exception inside a Qt slot (the ``_SlotErrors`` pattern).

Run with pytest::

    pytest -q tests/test_gui_analysis.py

or standalone::

    python tests/test_gui_analysis.py
"""

from __future__ import annotations

import json
import math
import os
import sys
import tempfile
import time
import traceback

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ.setdefault("QT_QPA_FONTDIR", "C:/Windows/Fonts")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import matplotlib                                                   # noqa: E402
matplotlib.use("Agg")

import numpy as np                                                  # noqa: E402

from geodome import benchmarks as bm                                # noqa: E402
from geodome import geometry as geo                                 # noqa: E402


# ================================================================ helpers ==
def _window():
    from qtpy import QtWidgets
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    from geodome.gui.main_window import MainWindow
    win = MainWindow()
    win.resize(1600, 1000)
    assert not win.extension_errors, win.extension_errors[-1]
    nl = win.nonlinear_extension
    nl.quiet = True                        # its messages go to the log
    app.processEvents()
    return app, win


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


class _Dialogs:
    """Replace the blocking message boxes and file dialogs.

    Message boxes are recorded in ``messages`` as (kind, title, text);
    questions answer No.  The file dialogs return the given paths.
    """

    def __init__(self, open_path="", save_path="", folder=""):
        self.open_path = open_path
        self.save_path = save_path
        self.folder = folder
        self.messages = []

    def __enter__(self):
        from qtpy import QtWidgets
        mb, fd = QtWidgets.QMessageBox, QtWidgets.QFileDialog
        self._saved = [(mb, n, getattr(mb, n)) for n in
                       ("warning", "critical", "information", "question",
                        "about")]
        self._saved += [(fd, n, getattr(fd, n)) for n in
                        ("getOpenFileName", "getSaveFileName",
                         "getExistingDirectory")]

        def rec(kind):
            def f(parent, title, text, *a, **k):
                self.messages.append((kind, str(title), str(text)))
                return mb.Ok
            return staticmethod(f)

        for n in ("warning", "critical", "information", "about"):
            setattr(mb, n, rec(n))
        mb.question = staticmethod(
            lambda parent, title, text, *a, **k:
            self.messages.append(("question", str(title), str(text)))
            or mb.No)
        fd.getOpenFileName = staticmethod(lambda *a, **k: (self.open_path, ""))
        fd.getSaveFileName = staticmethod(lambda *a, **k: (self.save_path, ""))
        fd.getExistingDirectory = staticmethod(lambda *a, **k: self.folder)
        return self

    def __exit__(self, *exc):
        for owner, name, fn in self._saved:
            setattr(owner, name, fn)
        return False


def _wait(app, win, timeout=600.0):
    """Process events until the window is no longer busy."""
    t0 = time.time()
    while win._busy:
        if win._worker is not None:
            win._worker.wait(20)
        app.processEvents()
        if time.time() - t0 > timeout:
            raise AssertionError("the worker did not finish in time")
    for _ in range(3):
        app.processEvents()


def _analyse(app, win):
    win.on_analyse()
    assert win._busy, "the analysis did not start"
    _wait(app, win)
    assert win.results is not None, win.log_view.toPlainText()[-2000:]
    return win.results


def _small_dome(app, win, freq=2, radius=6.0, fraction=0.5):
    gp = win.geometry_panel
    gp.frequency.setValue(freq)
    gp.radius.setValue(radius)
    gp.fraction.setValue(fraction)
    gp.generate_btn.click()
    app.processEvents()
    assert win.model.geom.frequency == freq


def _titles(tabs):
    return [tabs.tabText(i) for i in range(tabs.count())]


def _plain(text):
    return text.replace("&", "")


def _menu(win, title):
    for a in win.menuBar().actions():
        if a.menu() is not None and _plain(a.text()) == title:
            return a.menu()
    raise AssertionError(f"no menu {title}")


def _menu_texts(menu):
    return [_plain(a.text()) for a in menu.actions() if not a.isSeparator()]


def _close(app, win):
    """Run pending draw_idle calls, then close (no dialogs may block)."""
    with _Dialogs():
        app.processEvents()
        win.close()
        app.processEvents()


def _dxf_entities(path):
    """{entity type: [(layer, [(x, y, z), ...]), ...]} of a DXF file."""
    with open(path, encoding="ascii") as fh:
        lines = [ln.rstrip("\r\n") for ln in fh]
    pairs = [(lines[i].strip(), lines[i + 1].strip())
             for i in range(0, len(lines) - 1, 2)]
    out, inside, cur = {}, False, None
    for code, value in pairs:
        if code == "2" and value == "ENTITIES":
            inside = True
            continue
        if not inside:
            continue
        if code == "0":
            if value == "ENDSEC":
                break
            cur = {"type": value, "layer": None, "xyz": {}}
            out.setdefault(value, []).append(cur)
        elif cur is not None and code == "8":
            cur["layer"] = value
        elif cur is not None and code in ("10", "20", "30", "11", "21", "31",
                                          "12", "22", "32", "13", "23", "33"):
            cur["xyz"][code] = float(value)
    return out


# ================================================================== tests ==
def test_layout_menus_toolbar_tabs():
    """Every add-on is in, once, and in the agreed order."""
    from qtpy import QtWidgets
    from geodome.gui import (nonlinear_view, modal_view, fabrication_view,
                             study_view, validation_view)
    app, win = _window()
    with _SlotErrors() as errors:
        kinds = [type(e).__name__ for e in win.extensions]
        assert kinds == ["NonlinearExtension", "ModalExtension",
                         "FabricationExtension", "StudyExtension",
                         "ValidationExtension"], kinds
        # a second install returns the same object and adds nothing
        n_view, n_res = win.view_tabs.count(), win.result_tabs.count()
        for mod, attr in ((nonlinear_view, "nonlinear_extension"),
                          (modal_view, "modal_extension"),
                          (fabrication_view, "fabrication_extension"),
                          (study_view, "study_extension"),
                          (validation_view, "validation_extension")):
            assert mod.install(win) is getattr(win, attr)
        assert len(win.extensions) == 5
        assert (win.view_tabs.count(), win.result_tabs.count()) == \
            (n_view, n_res)

        assert _titles(win.view_tabs) == [
            "Model", "Loads", "Deflected shape", "Axial force", "Utilisation",
            "Buckling mode", "Charts", "Nonlinear", "Modes", "Study"]
        assert _titles(win.result_tabs) == [
            "Summary", "Members", "Nodes", "Reactions", "Auto-sizing",
            "Verification", "Nonlinear", "Imperfection", "Modes",
            "Fabrication", "Study", "Validation", "Log"]

        tops = [_plain(a.text()) for a in win.menuBar().actions()]
        assert tops == ["File", "Examples", "Analysis", "View", "Validation",
                        "Help"], tops
        analysis = _menu_texts(_menu(win, "Analysis"))
        assert analysis == ["Generate dome", "Run analysis",
                            "Auto-size strut families", "Modal analysis...",
                            "Nonlinear analysis...", "Imperfection study...",
                            "Parameter study..."], analysis
        # one "Run verification suite" entry in the whole menu bar, on F11
        every = [t for a in win.menuBar().actions()
                 for t in _menu_texts(a.menu())]
        assert len(every) == len(set(every)), every
        verify = [t for t in every if t.startswith("Run verification suite")]
        assert verify == ["Run verification suite (75 checks)"], verify
        assert win.validation_extension.act_verify.shortcut().toString() \
            == "F11"
        acts = _menu(win, "Analysis").actions()       # no double lines
        assert not acts[0].isSeparator() and not acts[-1].isSeparator()
        assert not any(a.isSeparator() and b.isSeparator()
                       for a, b in zip(acts, acts[1:]))
        file_menu = _menu_texts(_menu(win, "File"))
        assert file_menu[-1] == "Exit"
        assert "Export geometry to CAD (DXF)..." in file_menu

        # toolbar: the original buttons first, then Nonlinear and Modes
        tb = [a.text() for a in win.toolbar.actions() if not a.isSeparator()]
        assert tb == ["Generate dome", "Run analysis", "Auto-size", "Verify",
                      "Export report", "Nonlinear", "Modes"], tb
        assert len([b for b in win.findChildren(QtWidgets.QToolBar)
                    if b.parent() is win]) == 1

        # the View menu lists the add-on view tabs when it opens
        win.menu_view.aboutToShow.emit()
        names = [a.text() for a in win.menu_view.actions()]
        for t in ("Nonlinear", "Modes", "Study"):
            assert names.count(t) == 1, names

        # Help > Theory is the HTML page of gui/theory.py
        from geodome.gui.theory import theory_html
        dlg = win.theory_dialog()
        assert isinstance(dlg.browser, QtWidgets.QTextBrowser)
        text = dlg.browser.toPlainText()
        for part in ("Nonlinear analysis", "Modal analysis",
                     "Fabrication geometry", "References"):
            assert part in text, part
        assert len(theory_html()) > 5000
        dlg.close()

        # every extension answers both hooks
        for e in win.extensions:
            assert callable(e.on_model_changed) and callable(e.on_results)
        assert win.windowTitle().startswith("GeoDomeLab 2.0.0")
        assert not errors, errors[-1]
    _close(app, win)


def test_every_preset_loads_from_the_examples_menu():
    app, win = _window()
    with _SlotErrors() as errors, _Dialogs() as dlg:
        assert set(win.example_actions) == set(bm.preset_keys())
        for key in bm.preset_keys():
            title, _ = bm.preset_info(key)
            win.example_actions[key].trigger()
            app.processEvents()
            ref = bm.make_preset(key)
            g = win.model.geom
            assert g.n_nodes == ref.geom.n_nodes, key
            assert g.n_members == ref.geom.n_members, key
            assert np.allclose(g.nodes, ref.geom.nodes, atol=1e-12)
            assert title in win.windowTitle(), (key, win.windowTitle())
            assert win.results is None
            # the add-ons were told: stale results are gone, tables rebuilt
            assert win.nonlinear_extension.nonlinear_result is None
            assert win.modal_extension.modal_result is None
            fx = win.fabrication_extension
            if fx.fabrication_tables is not None:
                assert len(fx.fabrication_tables["members"]) == g.n_members
            assert f"Example loaded: {title}" in win.log_view.toPlainText()
        assert not dlg.messages, dlg.messages
        assert not errors, errors[-1]
    _close(app, win)


def test_analysis_of_the_6v_case_study():
    app, win = _window()
    with _SlotErrors() as errors, _Dialogs() as dlg:
        win.example_actions["case_study_6v"].trigger()
        app.processEvents()
        g = win.model.geom
        assert g.frequency == 6
        res = _analyse(app, win)
        # the results keep a snapshot of the model (review GUI-03)
        assert res.model is not win.model
        assert res.model.to_dict() == win.model.to_dict()
        assert win.member_table.rowCount() == g.n_members
        assert res.buckling is not None and len(res.buckling.load_factors)
        assert math.isfinite(res.governing_util) and res.governing_util > 0
        assert win.status_label.text().startswith("Analysis complete")
        assert win.view_tabs.currentWidget() is win.canvas_util
        # the add-ons heard about the results
        fx = win.fabrication_extension
        assert fx.fabrication_tables is not None
        assert len(fx.fabrication_tables["members"]) == g.n_members
        # nothing is busy; every long-run action is back on
        assert not win._busy
        for a in win.busy_actions:
            if a is win.validation_extension.act_opensees:
                continue
            assert a.isEnabled(), a.text()
        extras = win.extension_results()
        assert set(extras) == {"fabrication_tables"}
        assert not dlg.messages, dlg.messages
        assert not errors, errors[-1]
    _close(app, win)


def test_nonlinear_von_mises_limit_point_and_busy_state():
    from geodome.gui import nonlinear_view as nv
    app, win = _window()
    ext = win.nonlinear_extension
    with _SlotErrors() as errors, _Dialogs() as dlg:
        assert win.load_preset("von_mises")
        s = nv.NonlinearDialog(win, win.model).values()
        assert s["method"] == "arc-length" and s["control_axis"] == "uz"
        assert ext.start_nonlinear(s)
        # while it runs, no other long run can start
        assert win._busy and ext.is_running()
        md, sy, va = (win.modal_extension, win.study_extension,
                      win.validation_extension)
        for a in (md.act_modal, md.act_modes, md.run_btn, sy.act_study,
                  va.act_verify, va.act_hand, va.act_opensees, va.act_export,
                  win.act_analyse, win.act_analyse_menu,
                  ext.act_nonlinear, ext.act_toolbar):
            assert not a.isEnabled(), a.text()
        va.menu.aboutToShow.emit()            # opening the menu keeps it off
        assert not va.act_hand.isEnabled()
        assert not win.menu_examples.isEnabled()
        # an export (Ctrl+R) waits; it must not end the busy state
        n_box = len(dlg.messages)
        win.on_export_html()
        win.on_export_all()
        assert win._busy and not md.act_modes.isEnabled()
        assert len(dlg.messages) == n_box
        assert ext.wait(300) == "done"
        app.processEvents()
        assert not win._busy
        for a in (md.act_modal, md.act_modes, md.run_btn, sy.act_study,
                  va.act_verify, va.act_hand, win.act_analyse,
                  ext.act_nonlinear):
            assert a.isEnabled(), a.text()

        r = ext.nonlinear_result
        a_, h, EA, P = bm.von_mises_parameters(win.model)
        lp = bm.von_mises_limit_point(a_, h, EA)
        assert r.limit_found and r.limit_type == "limit"
        assert abs(r.lambda_limit - lp["P"] / P) <= 1e-6 * lp["P"] / P
        assert abs(-r.control_limit - lp["w"]) <= 1e-6 * lp["w"]
        # the plot is filled
        ax = ext.view.path_canvas.ax_path
        labels = [ln.get_label() for ln in ax.get_lines()]
        assert "Corotational (large displacement)" in labels
        cor = [ln for ln in ax.get_lines()
               if ln.get_label() == "Corotational (large displacement)"][0]
        assert len(cor.get_xdata()) == r.n_points > 10
        assert "Limit point" in " ".join(t.get_text() for t in ax.texts)
        assert ext.view.slider.value() == r.limit_step
        assert ext.path_table.rowCount() == r.n_points
        assert win.extension_results()["nonlinear_result"] is r
        # a new model clears it
        win.load_preset("star_dome")
        assert ext.nonlinear_result is None
        assert "nonlinear_result" not in win.extension_results()
        assert not ext.errors, ext.errors[-1]
        assert not [m for m in dlg.messages if m[0] != "question"], \
            dlg.messages
        assert not errors, errors[-1]
    _close(app, win)


def test_modal_model_a_matches_hand_calc_f7():
    sys.path.insert(0, os.path.join(ROOT, "validation"))
    import hand_calc_ext as hce
    hand = hce.results()
    fz = hand["F7.f_z_dead"]["value"]
    fxy = hand["F7.f_xy_dead"]["value"]
    app, win = _window()
    md = win.modal_extension
    with _SlotErrors() as errors, _Dialogs() as dlg:
        assert win.load_preset("validation_a")
        assert md.run_modal({"n_modes": 12, "include_dead": True})
        assert win._busy and not win.nonlinear_extension.act_toolbar \
            .isEnabled()
        _wait(app, win)
        res = md.modal_result
        assert res is not None and res.n_modes == 3
        f = sorted(float(x) for x in res.frequencies)
        assert abs(f[0] - fz) <= 1e-9 * fz, (f, fz)
        assert abs(f[1] - fxy) <= 1e-9 * fxy and abs(f[2] - fxy) <= 1e-9 * fxy
        assert md.table.rowCount() == 3
        assert md.mode_combo.count() == 3
        assert win.view_tabs.currentWidget() is md.view
        md.start_animation()
        t0 = time.time()
        while time.time() - t0 < 0.3:
            app.processEvents()
            time.sleep(0.01)
        assert md.is_animating()
        win.view_tabs.setCurrentWidget(win.canvas_model)   # stops it
        app.processEvents()
        assert not md.is_animating()
        assert win.extension_results()["modal_result"] is res
        assert not dlg.messages, dlg.messages
        assert not errors, errors[-1]
    _close(app, win)


def test_fabrication_tables_of_a_2v_dome():
    app, win = _window()
    fx = win.fabrication_extension
    with _SlotErrors() as errors, _Dialogs() as dlg:
        _small_dome(app, win, 2, 6.0)
        t = fx.fabrication_tables
        assert t is not None
        cf = sorted(r["chord_factor"] for r in t["struts"])
        assert abs(cf[0] - 0.54653) < 5e-6 and abs(cf[1] - 0.61803) < 5e-6
        assert sum(r["count"] for r in t["struts"]) == \
            win.model.geom.n_members
        assert fx.tables["struts"].rowCount() == len(t["struts"]) == 2
        assert fx.tables["hubs"].rowCount() == len(t["hubs"]) > 0
        assert fx.tables["panels"].rowCount() == len(t["panels"]) > 0
        # hub offset: cut length = L - 2 offset
        fx.hub_offset.setValue(50.0)
        tt = fx.refresh()
        for r in tt["struts"]:
            assert abs(r["cut_length_mm"] - (r["length_mm"] - 100.0)) < 1e-6
        # cut list export
        with tempfile.TemporaryDirectory() as tmp:
            files = fx.export_csv(tmp)
            assert files and all(os.path.isfile(p) for p in files)
        assert win.extension_results()["fabrication_tables"] is \
            fx.fabrication_tables
        assert not dlg.messages, dlg.messages
        assert not errors, errors[-1]
    _close(app, win)


def test_parameter_study_two_frequencies():
    app, win = _window()
    sy = win.study_extension
    with _SlotErrors() as errors, _Dialogs() as dlg:
        assert sy.run_study([2, 3], radii=[8.0], fractions=[0.5])
        assert win._busy and not win.modal_extension.act_modes.isEnabled()
        _wait(app, win)
        rows = sy.study_rows
        assert rows is not None and len(rows) == 2
        assert [r["frequency"] for r in rows] == [2, 3]
        assert all(r["status"] != "ERROR" for r in rows)
        assert sy.table.rowCount() == 2
        assert len(sy.canvas.figure.axes) == 4
        assert win.view_tabs.currentWidget() is sy.view_tab
        assert win.result_tabs.currentWidget() is sy.result_tab
        # the study survives a model change (it is not tied to one model)
        win.load_preset("validation_b")
        assert sy.study_rows is rows
        assert win.extension_results()["study_rows"] is rows
        assert not dlg.messages, dlg.messages
        assert not errors, errors[-1]
    _close(app, win)


def test_validation_menu_with_opensees_missing():
    from geodome import opensees_bridge as ob
    app, win = _window()
    va = win.validation_extension
    saved = ob.opensees_status
    ob.opensees_status = lambda: {"mode": None, "python": None}
    try:
        with _SlotErrors() as errors, _Dialogs() as dlg:
            va.menu.aboutToShow.emit()
            assert not va.act_opensees.isEnabled()
            assert "OpenSees was not found" in va.act_opensees.toolTip()
            # the action refuses politely; nothing runs
            assert va.run_opensees_crosscheck() is False
            assert not win._busy
            assert dlg.messages and "OpenSees was not found" in \
                dlg.messages[-1][2]
            # a run of something else does not switch it back on
            _analyse(app, win)
            assert not va.act_opensees.isEnabled()
            # the script export does not need OpenSees
            with tempfile.TemporaryDirectory() as tmp:
                path = va.export_script(path=os.path.join(tmp, "s.py"))
                assert path and os.path.isfile(path)
                with open(path, encoding="utf-8") as fh:
                    compile(fh.read(), path, "exec")
            # the hand-calculation validation runs in the worker
            assert va.run_hand_validation()
            assert win._busy and not va.act_hand.isEnabled()
            _wait(app, win)
            h = va.hand_result
            assert h is not None and h["n"] > 0 and h["n_pass"] == h["n"], \
                (h or {}).get("n_pass")
            assert va.hand_table.rowCount() == h["n"]
            assert win.result_tabs.currentWidget() is va.result_tab
            assert va.act_hand.isEnabled() and not va.act_opensees.isEnabled()
            # the verification entry runs the window's suite (F11)
            calls = []
            win.on_verify = lambda: calls.append(1)
            va.act_verify.trigger()
            assert calls == [1]
            assert not errors, errors[-1]
    finally:
        ob.opensees_status = saved
        va.refresh_actions()
    _close(app, win)


def test_export_script_menu_entry_writes_the_file():
    """Review GUI-01: Validation > Export OpenSeesPy script... clicked from
    the menu (QAction.triggered passes checked=False) must open the dialog
    and write the file."""
    from qtpy import QtWidgets
    from geodome.gui import validation_view as vv
    app, win = _window()
    va = win.validation_extension
    saved = vv.ExportScriptDialog.exec
    vv.ExportScriptDialog.exec = lambda self: 1          # press OK
    try:
        with _SlotErrors() as errors, tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "dome_opensees.py")
            with _Dialogs(save_path=path) as dlg:
                va.menu.aboutToShow.emit()
                assert va.act_export.isEnabled()
                va.act_export.trigger()
                app.processEvents()
                assert os.path.isfile(path), dlg.messages
                assert not dlg.messages, dlg.messages
                with open(path, encoding="utf-8") as fh:
                    compile(fh.read(), path, "exec")
            assert not errors, errors[-1]
    finally:
        vv.ExportScriptDialog.exec = saved
    _close(app, win)


def test_results_never_outlive_their_inputs():
    """Review GUI-03: after a change of section, material, supports or an
    auto-size, the old static, modal and nonlinear results are gone; the
    results keep a snapshot of the model they were computed for."""
    app, win = _window()
    md = win.modal_extension
    with _SlotErrors() as errors, _Dialogs() as dlg:
        assert win.load_preset("validation_a")
        res = _analyse(app, win)
        assert res.model is not win.model
        assert md.run_modal({"n_modes": 3, "include_dead": True})
        _wait(app, win)
        assert md.modal_result is not None
        # a new section on every family
        sp = win.section_panel
        sp.apply_combo.setCurrentText("CHS 48.3 x 3.2")
        sp._apply_all()
        app.processEvents()
        assert win.model.group_sections == ["CHS 48.3 x 3.2"]
        assert res.model.group_sections == ["CHS 88.9 x 4.0"]  # snapshot
        assert win.results is None and md.modal_result is None
        assert "modal_result" not in win.extension_results()
        assert win.status_label.text().startswith("Inputs changed")
        n_box = len(dlg.messages)
        win.on_export_html()
        assert len(dlg.messages) == n_box + 1
        assert "Run the analysis first" in dlg.messages[-1][2]
        res2 = _analyse(app, win)
        assert res2.model.group_sections == ["CHS 48.3 x 3.2"]
        # material
        sp.material.setCurrentText("EN S355")
        app.processEvents()
        assert win.results is None
        # supports
        res3 = _analyse(app, win)
        gp = win.geometry_panel
        gp.set_support_type(geo.ROLLER3)
        app.processEvents()
        assert win.results is None and res3.model.support_type == geo.PINNED
        # a Loads tab edit: told, results kept until F9
        res4 = _analyse(app, win)
        win.load_panel.snow.setValue(win.load_panel.snow.value() + 0.5)
        app.processEvents()
        assert win.results is res4
        assert win.status_label.text().startswith("Inputs changed")
        # auto-size drops the modal result of the old sections
        assert md.run_modal({"n_modes": 3, "include_dead": True})
        _wait(app, win)
        assert md.modal_result is not None
        win.on_auto_size()
        _wait(app, win)
        _wait(app, win)                    # the re-analysis after sizing
        assert md.modal_result is None
        assert win.results is not None
        assert win.results.model.group_sections == win.model.group_sections
        assert not errors, errors[-1]
    _close(app, win)


def test_html_report_with_every_extra_section():
    from geodome.gui import nonlinear_view as nv
    from geodome import report as rp
    app, win = _window()
    nl, md, sy = (win.nonlinear_extension, win.modal_extension,
                  win.study_extension)
    tmp = tempfile.mkdtemp(prefix="geodome_report_")
    with _SlotErrors() as errors:
        # a study first: it is kept when the model changes
        assert sy.run_study([1, 2], radii=[6.0], fractions=[0.5])
        _wait(app, win)
        _small_dome(app, win, 2, 6.0)
        _analyse(app, win)
        s = nv.NonlinearDialog(win, win.model).values()
        assert nl.start_nonlinear(s, synchronous=True)
        dlg_i = nv.ImperfectionDialog(win, win.model)
        dlg_i.amplitudes.setText("span/300")
        assert nl.start_imperfection(dlg_i.values(), synchronous=True)
        assert md.run_modal({"n_modes": 6})
        _wait(app, win)
        extras = win.extension_results()
        assert set(extras) == set(rp.EXTRA_KEYS), sorted(extras)

        # File > Export HTML report...
        path = os.path.join(tmp, "report.html")
        with _Dialogs(save_path=path) as dlg:
            win.on_export_html()
        assert ("question", "Report", "Report written. Open it in your "
                "browser?") in dlg.messages, dlg.messages
        with open(path, encoding="utf-8") as fh:
            text = fh.read()
        for sid in ("nonlinear", "imperfection", "modal", "fabrication",
                    "study"):
            assert f"id='{sid}'" in text, sid
        # sections 4 (auto-sizing) and 7 (benchmark suite) need those runs
        for num in (1, 2, 3, 5, 6, 8, 9, 10, 11, 12):
            assert f">{num} &nbsp;" in text, num
        assert text.count("data:image/png;base64") >= 8
        assert "could not be built" not in text
        assert text.rstrip().endswith("</html>")

        # File > Export everything to a folder...
        out = os.path.join(tmp, "all")
        os.makedirs(out)
        with _Dialogs(folder=out) as dlg:
            win.on_export_all()
        assert not [m for m in dlg.messages if m[0] == "critical"], \
            dlg.messages
        with open(os.path.join(out, "dome_report.html"),
                  encoding="utf-8") as fh:
            text2 = fh.read()
        for sid in ("nonlinear", "imperfection", "modal", "fabrication",
                    "study"):
            assert f"id='{sid}'" in text2, sid
        assert os.path.isfile(os.path.join(out, "dome_model.json"))
        assert not nl.errors, nl.errors[-1]
        assert not errors, errors[-1]
    _close(app, win)


def test_dxf_export_from_the_file_menu():
    app, win = _window()
    fx = win.fabrication_extension
    tmp = tempfile.mkdtemp(prefix="geodome_dxf_")
    with _SlotErrors() as errors:
        _small_dome(app, win, 2, 6.0)
        g = win.model.geom
        path = os.path.join(tmp, "dome.dxf")
        with _Dialogs(save_path=path) as dlg:
            fx.act_dxf.trigger()
        assert os.path.isfile(path), dlg.messages
        ent = _dxf_entities(path)
        assert len(ent["LINE"]) == g.n_members
        assert len(ent["3DFACE"]) == len(g.faces)
        assert len(ent["POINT"]) == g.n_nodes
        nodes = np.asarray(g.nodes)
        for e, (i, j) in zip(ent["LINE"], np.asarray(g.members)):
            p = [e["xyz"][c] for c in ("10", "20", "30")]
            q = [e["xyz"][c] for c in ("11", "21", "31")]
            assert np.allclose(p, nodes[i], atol=1e-9)
            assert np.allclose(q, nodes[j], atol=1e-9)
            assert e["layer"].startswith("STRUT_")
        # coloured by utilisation once there are results
        _analyse(app, win)
        p2 = fx.export_dxf(os.path.join(tmp, "util.dxf"), by_utilisation=True)
        assert p2 and len(_dxf_entities(p2)["LINE"]) == g.n_members
        assert not errors, errors[-1]
    _close(app, win)


def test_save_open_round_trip():
    app, win = _window()
    tmp = tempfile.mkdtemp(prefix="geodome_json_")
    with _SlotErrors() as errors:
        assert win.load_preset("three_bar")
        win.load_panel.delta_T.setValue(25.0)
        _analyse(app, win)
        path = os.path.join(tmp, "three_bar.json")
        with _Dialogs(save_path=path) as dlg:
            win.on_save()
        assert os.path.isfile(path), dlg.messages
        saved = win.model.to_dict()
        win.on_new()
        assert win.model.geom.frequency == 4 and not win.model.geom.is_custom
        with _Dialogs(open_path=path) as dlg:
            win.on_open()
        assert not dlg.messages, dlg.messages
        assert win.model.geom.is_custom
        assert json.dumps(win.model.to_dict(), sort_keys=True) == \
            json.dumps(saved, sort_keys=True)
        assert "three_bar.json" in win.windowTitle()
        assert win.results is None
        assert win.load_panel.delta_T.value() == 25.0
        # the opened model analyses like the saved one
        r = _analyse(app, win)
        assert max(c.equilibrium_residual for c in r.cases.values()) < 1e-9
        assert not errors, errors[-1]
    _close(app, win)


def test_no_slot_errors_through_the_menus():
    """Open every add-on dialog from its menu entry and cancel it; switch
    every view and result tab; save each view tab as PNG."""
    from qtpy import QtWidgets
    app, win = _window()
    tmp = tempfile.mkdtemp(prefix="geodome_views_")
    rejected = []

    def reject_later():
        for w in QtWidgets.QApplication.topLevelWidgets():
            if isinstance(w, QtWidgets.QDialog) and w.isVisible():
                rejected.append(w.windowTitle())
                w.reject()

    with _SlotErrors() as errors, _Dialogs() as dlg:
        win.load_preset("star_dome")
        acts = {_plain(a.text()): a for a in win.menu_analysis.actions()}
        timer = __import__("qtpy").QtCore.QTimer()
        timer.timeout.connect(reject_later)
        timer.start(50)
        try:
            for name in ("Modal analysis...", "Nonlinear analysis...",
                         "Imperfection study...", "Parameter study..."):
                acts[name].trigger()
                app.processEvents()
        finally:
            timer.stop()
        assert len(rejected) >= 3, rejected
        assert not win._busy
        for i in range(win.view_tabs.count()):
            win.view_tabs.setCurrentIndex(i)
            app.processEvents()
            w = win.view_tabs.currentWidget()
            if hasattr(w, "save"):
                p = os.path.join(tmp, f"view_{i}.png")
                w.save(p, dpi=40)
                assert os.path.getsize(p) > 1000
        for i in range(win.result_tabs.count()):
            win.result_tabs.setCurrentIndex(i)
            app.processEvents()
        assert not [m for m in dlg.messages if m[0] == "critical"], \
            dlg.messages
        assert not win.nonlinear_extension.errors
        assert not errors, errors[-1]
    _close(app, win)


def test_view_menu_works_and_a_closed_window_is_freed():
    """The View menu lists and opens the add-on tabs.  Rebuilding it must not
    keep the window alive: lambda slots holding the window did, and
    collecting such a window later crashed the test process."""
    import gc
    import weakref
    app, win = _window()
    with _SlotErrors() as errors:
        for _ in range(2):
            win.menu_view.aboutToShow.emit()
        acts = {a.text(): a for a in win.menu_view.actions()}
        for name in ("Modes", "Study", "Nonlinear", "Model"):
            acts[name].trigger()
            app.processEvents()
            assert win.view_tabs.tabText(win.view_tabs.currentIndex()) == \
                name
        assert not errors, errors[-1]
    _close(app, win)
    ref = weakref.ref(win)
    del win, acts
    gc.collect()
    app.processEvents()
    gc.collect()
    assert ref() is None, "the closed window is still alive"


# ================================================================ driver ==
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
