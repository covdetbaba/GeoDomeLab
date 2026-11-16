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
* the interactive 3D viewer (gui/view3d.py): the legend column stays left
  of the model (1V, 6V, 10V, two sizes, and the static export), picking by
  a click, Ctrl + click and a click on empty space, the wheel zoom at the
  cursor, the right-drag pan, the keys, the view bar and the Find box, the
  link with the Members and Nodes tables, the camera kept on a redraw, and
  the Modes and Nonlinear views with the viewer on; a click next to a node
  picks the node, the details follow the results (analysis, input edit, new
  model), a legend of hundreds of families is cut to fit, the z axis and a
  two-line title stay in the figure, a fast wheel spin, a toolbar zoom kept
  on a redraw, figure exports with the hover on, and the number placement.
  The events go through the real canvas (QTest and Qt wheel and mouse
  events);
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
    import qtpy
    if qtpy.API_NAME.lower().startswith("pyside"):
        # PySide does not collect a Python reference cycle that runs through
        # a Qt parent (each add-on is a child of the window and keeps a
        # reference to it), so a closed window lives until the process
        # ends.  PyQt5 and PyQt6 free it.
        return
    ref = weakref.ref(win)
    del win, acts
    gc.collect()
    app.processEvents()
    gc.collect()
    assert ref() is None, "the closed window is still alive"


# ====================================================== 3D viewer tests ==
def _shown_window(size=(1600, 1000)):
    """A main window on the (offscreen) screen, so the canvases have their
    real sizes and take mouse and key events."""
    app, win = _window()
    win.resize(*size)
    win.show()
    app.processEvents()
    return app, win


def _draw(app, canvas):
    """Run pending redraws, then draw the canvas now."""
    app.processEvents()
    canvas.canvas.draw()
    app.processEvents()


def _qpoint(cv, x, y):
    """Qt widget point of a matplotlib display point (y up, device px)."""
    from qtpy import QtCore
    r = float(getattr(cv, "device_pixel_ratio", 1.0) or 1.0)
    return QtCore.QPoint(int(round(x / r)),
                         int(round((cv.figure.bbox.height - y) / r)))


def _click(app, cv, x, y, ctrl=False):
    from qtpy.QtCore import Qt
    from qtpy.QtTest import QTest
    mod = (Qt.KeyboardModifier.ControlModifier if ctrl
           else Qt.KeyboardModifier.NoModifier)
    QTest.mouseClick(cv, Qt.MouseButton.LeftButton, mod, _qpoint(cv, x, y))
    app.processEvents()


def _wheel_at(app, cv, x, y, notches=1, angle=None, pixel=0):
    """One Qt wheel event: ``notches`` mouse notches (angle delta 120
    each), or the given angle and pixel deltas (a touchpad)."""
    from qtpy import QtCore, QtGui, QtWidgets
    pos = QtCore.QPointF(_qpoint(cv, x, y))
    gpos = QtCore.QPointF(cv.mapToGlobal(pos.toPoint()))
    angle = 120 * notches if angle is None else angle
    ev = QtGui.QWheelEvent(pos, gpos, QtCore.QPoint(0, int(pixel)),
                           QtCore.QPoint(0, int(angle)),
                           QtCore.Qt.MouseButton.NoButton,
                           QtCore.Qt.KeyboardModifier.NoModifier,
                           QtCore.Qt.ScrollPhase.NoScrollPhase, False)
    QtWidgets.QApplication.sendEvent(cv, ev)
    app.processEvents()


def _drag(app, cv, x0, y0, x1, y1, button, steps=6):
    """Press, move in steps with the button held, release."""
    from qtpy import QtCore, QtGui, QtWidgets
    from qtpy.QtCore import Qt
    from qtpy.QtTest import QTest
    QTest.mousePress(cv, button, Qt.KeyboardModifier.NoModifier,
                     _qpoint(cv, x0, y0))
    for s in range(1, steps + 1):
        p = _qpoint(cv, x0 + (x1 - x0) * s / steps,
                    y0 + (y1 - y0) * s / steps)
        ev = QtGui.QMouseEvent(QtCore.QEvent.Type.MouseMove,
                               QtCore.QPointF(p),
                               QtCore.QPointF(cv.mapToGlobal(p)),
                               Qt.MouseButton.NoButton, button,
                               Qt.KeyboardModifier.NoModifier)
        QtWidgets.QApplication.sendEvent(cv, ev)
    QTest.mouseRelease(cv, button, Qt.KeyboardModifier.NoModifier,
                       _qpoint(cv, x1, y1))
    app.processEvents()


def _move(app, cv, x, y):
    """Move the mouse with no button held (hover)."""
    from qtpy import QtCore, QtGui, QtWidgets
    from qtpy.QtCore import Qt
    p = _qpoint(cv, x, y)
    ev = QtGui.QMouseEvent(QtCore.QEvent.Type.MouseMove, QtCore.QPointF(p),
                           QtCore.QPointF(cv.mapToGlobal(p)),
                           Qt.MouseButton.NoButton, Qt.MouseButton.NoButton,
                           Qt.KeyboardModifier.NoModifier)
    QtWidgets.QApplication.sendEvent(cv, ev)
    app.processEvents()


def _key(app, cv, key, shift=False, ctrl=False):
    from qtpy.QtCore import Qt
    from qtpy.QtTest import QTest
    mod = Qt.KeyboardModifier.NoModifier
    if shift:
        mod = Qt.KeyboardModifier.ShiftModifier
    if ctrl:
        mod = Qt.KeyboardModifier.ControlModifier
    QTest.keyClick(cv, key, mod)
    app.processEvents()


def _isolated_member(ctl, skip=()):
    """(member, x, y): a front member whose midpoint is at least 10 px from
    every other member and 12 px from every node on screen."""
    from geodome.gui import view3d as v3
    sc = ctl.scene
    px, py, pz = v3.project(ctl.ax, sc.nodes)
    i, j = sc.members[:, 0], sc.members[:, 1]
    mx, my = 0.5 * (px[i] + px[j]), 0.5 * (py[i] + py[j])
    reg = ctl._region_px()
    for k in np.argsort(0.5 * (pz[i] + pz[j])):
        k = int(k)
        if k in skip or not reg.contains(mx[k], my[k]):
            continue
        d, _ = v3.segment_distance(mx[k], my[k], px[i], py[i], px[j], py[j])
        d[k] = np.inf
        if d.min() > 10.0 and np.hypot(px - mx[k], py - my[k]).min() > 12.0:
            return k, float(mx[k]), float(my[k])
    raise AssertionError("no isolated member on screen")


def _ranges(ctl):
    lo, hi = ctl._limits()
    return hi - lo


def _no_view_errors(win):
    for c in (win.canvas_model, win.canvas_loads, win.canvas_defl,
              win.canvas_axial, win.canvas_util, win.canvas_buck,
              win.modal_extension.view,
              win.nonlinear_extension.view.path_canvas):
        assert not c.view3d.errors, c.view3d.errors[-1]


def test_view3d_legend_column_stays_left_of_the_model():
    """The strut-family legend is a column left of the drawing: no overlap
    with the projected model for 1V, 6V and 10V (27 families) at two window
    sizes, every row inside the figure, and the same in the static export
    used by the report."""
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from geodome import plotting as plot
    from geodome import report as rp
    from geodome.gui import view3d as v3
    from geodome.model import default_model

    def check(fig, ax, model):
        leg = ax.get_legend()
        assert leg is not None and leg.get_visible()
        r = fig.canvas.get_renderer()
        lb = leg.get_window_extent(r)
        W, H = fig.bbox.width, fig.bbox.height
        n_fam = model.geom.n_strut_types()
        assert len(leg.get_texts()) == n_fam
        assert lb.x0 >= 0.0 and lb.y0 >= -0.5 and lb.y1 <= H + 0.5, \
            (lb.extents, W, H)
        px, py, _ = v3.project(ax, model.geom.nodes)
        assert lb.x1 + 2.0 < px.min(), (lb.x1, px.min())
        # the 3D axes itself starts right of the legend column
        assert ax.get_position().x0 * W >= lb.x1 - 0.5
        assert ax.get_axes_locator().region.x0 * W >= lb.x1
        return lb

    app, win = _shown_window()
    with _SlotErrors() as errors:
        c = win.canvas_model
        for size in ((1600, 1000), (1280, 600)):
            win.resize(*size)
            app.processEvents()
            for f in (1, 6, 10):
                win.set_model(default_model(f, 24.0, 0.5))
                _draw(app, c)
                check(c.figure, c.axes, win.model)
        # 10V at the small size needs the short labels or two columns
        assert c.axes.get_axes_locator().choice[:2] != (1, False)
        # L hides the legend and the drawing takes the width; L shows it
        ctl = c.view3d
        x_with = ctl.layout.region.x0
        _key(app, c.canvas, __import__("qtpy").QtCore.Qt.Key.Key_L)
        _draw(app, c)
        assert not ctl.legend_on and not c.axes.get_legend().get_visible()
        assert ctl.layout.region.x0 < x_with
        assert not ctl.bar.legend_act.isChecked()
        ctl.bar.legend_act.trigger()
        _draw(app, c)
        assert ctl.legend_on and c.axes.get_legend().get_visible()
        assert abs(ctl.layout.region.x0 - x_with) < 1e-9
        # the utilisation legend follows the same rule
        win.set_model(default_model(4, 12.0, 0.5))
        _analyse(app, win)
        cu = win.canvas_util
        _draw(app, cu)
        lb = cu.axes.get_legend().get_window_extent(
            cu.figure.canvas.get_renderer())
        px, _, _ = v3.project(cu.axes, win.model.geom.nodes)
        assert lb.x1 + 2.0 < px.min()
        assert not errors, errors[-1]
        _no_view_errors(win)
    _close(app, win)

    # the static figure of the report keeps the column too
    for f in (1, 6, 10):
        m = default_model(f, 24.0, 0.5)
        fig = rp._fig3d(plot.draw_model, m)
        FigureCanvasAgg(fig)
        fig.canvas.draw()
        check(fig, fig.axes[0], m)


def test_view3d_pick_select_and_the_result_tables():
    """A click on a member selects it (details name its length), Ctrl +
    click adds and removes, a click on empty space clears; the Members and
    Nodes tables follow both ways; the axial and deflected views pick on
    the coordinates they draw and show their own value first."""
    from qtpy.QtCore import Qt
    from geodome.gui import view3d as v3
    app, win = _shown_window()
    with _SlotErrors() as errors, _Dialogs():
        assert win.load_preset("case_study_6v")
        _analyse(app, win)
        win.view_tabs.setCurrentWidget(win.canvas_model)
        win.result_tabs.setCurrentWidget(win.summary_view)
        c = win.canvas_model
        ctl = c.view3d
        _draw(app, c)
        assert ctl.selection is win.selection
        k, x, y = _isolated_member(ctl)
        # hover: the member under the cursor is highlighted (blitted over
        # the cached drawing), and the highlight goes on empty space
        assert ctl._bg is not None                   # cached by the draw
        _move(app, c.canvas, x, y)
        assert ctl._hover_key == ("member", k)
        assert ctl._hover_artist.get_visible()
        assert v3.tooltip(ctl.scene, "member", k).startswith(f"Member {k} (")
        reg = ctl._region_px()
        _move(app, c.canvas, reg.x1 - 4.0, reg.y0 + 4.0)
        assert ctl._hover_key is None and not ctl._hover_artist.get_visible()
        _click(app, c.canvas, x, y)
        assert win.selection.members == {k} and not win.selection.nodes
        L = float(win.model.geom.strut_lengths[k])
        text = ctl.info_text()
        assert f"Member {k}" in text and f"Length: {L:.4f} m" in text, text
        assert f"{L:.4f} m" in ctl.info_panel.text.toPlainText()
        # the Members row follows; the result tab stays where it was
        rows = {i.row() for i in
                win.member_table.selectionModel().selectedRows()}
        assert rows == {k}, rows
        assert win.result_tabs.currentWidget() is win.summary_view
        # the highlight is drawn over the member
        hl = [a for a in ctl._sel_artists if a.get_linestyle() == "-"]
        xs, ys, zs = hl[0].get_data_3d()
        ends = win.model.geom.nodes[win.model.geom.members[k]]
        assert np.allclose(np.c_[xs, ys, zs][:2], ends)
        # Ctrl + click adds a second member, and removes it again
        k2, x2, y2 = _isolated_member(ctl, skip={k})
        _click(app, c.canvas, x2, y2, ctrl=True)
        assert win.selection.members == {k, k2}
        _click(app, c.canvas, x2, y2, ctrl=True)
        assert win.selection.members == {k}
        # a click on empty space clears the selection and the table rows
        _click(app, c.canvas, reg.x1 - 4.0, reg.y0 + 4.0)
        assert win.selection.is_empty()
        assert not win.member_table.selectionModel().selectedRows()
        # a node: click on it
        sc = ctl.scene
        px, py, pz = v3.project(ctl.ax, sc.nodes)
        n = next(int(i) for i in np.argsort(pz)
                 if np.sort(np.hypot(px - px[i], py - py[i]))[1] > 16.0
                 and reg.contains(px[i], py[i]))
        _click(app, c.canvas, px[n], py[n])
        assert win.selection.nodes == {n} and not win.selection.members
        assert f"Node {n}" in ctl.info_text()
        rows = {i.row() for i in
                win.node_table.selectionModel().selectedRows()}
        assert rows == {n}
        # a double-click zooms to the member under the cursor
        from qtpy.QtTest import QTest
        r_before = _ranges(ctl)
        QTest.mouseDClick(c.canvas, Qt.MouseButton.LeftButton,
                          Qt.KeyboardModifier.NoModifier,
                          _qpoint(c.canvas, x, y))
        app.processEvents()
        assert win.selection.members == {k}
        assert np.all(_ranges(ctl) < 0.5 * r_before)
        _key(app, c.canvas, Qt.Key.Key_Home)
        assert np.allclose(_ranges(ctl), r_before)
        # table -> 3D: a row picked in the Members table is highlighted
        win.member_table.selectRow(k2)
        app.processEvents()
        assert win.selection.members == {k2}
        assert ctl._sel_artists
        # the axial view picks the same member and shows N first
        win.view_tabs.setCurrentWidget(win.canvas_axial)
        ca = win.canvas_axial
        _draw(app, ca)
        cta = ca.view3d
        title, rows_a, _ = cta.selection_details()
        assert title == f"Member {k2}" and rows_a[0][0].startswith("N (")
        cta.selection.clear()
        ka, xa, ya = _isolated_member(cta)
        _click(app, ca.canvas, xa, ya)
        assert win.selection.members == {ka}
        # the deflected view picks on the deformed coordinates
        win.view_tabs.setCurrentWidget(win.canvas_defl)
        cd = win.canvas_defl
        _draw(app, cd)
        ctd = cd.view3d
        assert not np.allclose(ctd.scene.nodes, win.model.geom.nodes)
        kd, xd, yd = _isolated_member(ctd)
        _click(app, cd.canvas, xd, yd)
        assert win.selection.members == {kd}
        assert ctd.selection_details()[1][0][0].startswith("|u| at ends")
        # Esc clears
        _key(app, cd.canvas, Qt.Key.Key_Escape)
        assert win.selection.is_empty()
        # a new model clears the selection (the numbers change)
        _click(app, cd.canvas, xd, yd)
        assert win.selection.members
        win.set_model(win.model.__class__(geom=win.model.geom))
        assert win.selection.is_empty()
        assert not errors, errors[-1]
        _no_view_errors(win)
    _close(app, win)


def test_view3d_wheel_zoom_at_the_cursor_and_drag():
    """The wheel zooms at the cursor (1:1:1 kept, the point under the
    cursor stays), right and middle drag pan without zooming, left drag
    rotates, and nothing fights the toolbar's own pan mode."""
    from qtpy.QtCore import Qt
    from geodome.gui import view3d as v3
    from geodome.model import default_model
    app, win = _shown_window()
    with _SlotErrors() as errors:
        win.set_model(default_model(6, 24.0, 0.5))
        c = win.canvas_model
        cv = c.canvas
        ctl = c.view3d
        _draw(app, c)
        box = np.asarray(ctl.ax.get_box_aspect(), float)
        r0 = _ranges(ctl)
        assert np.allclose(r0 / r0[0], box / box[0])
        sc = ctl.scene
        px, py, pz = v3.project(ctl.ax, sc.nodes)
        reg = ctl._region_px()
        inside = [int(i) for i in np.argsort(pz) if reg.contains(px[i], py[i])]
        n = inside[len(inside) // 3]
        x, y = _qpoint(cv, px[n], py[n]).x(), _qpoint(cv, px[n], py[n]).y()
        x_mpl, y_mpl = px[n], py[n]
        _wheel_at(app, cv, x_mpl, y_mpl, +1)
        r1 = _ranges(ctl)
        assert np.allclose(r1, r0 / v3.WHEEL_FACTOR), (r1, r0)
        assert np.allclose(r1 / r1[0], box / box[0])        # still 1:1:1
        qx, qy, _ = v3.project(ctl.ax, sc.nodes[[n]])
        # the point under the cursor stays (the cursor sits on a whole
        # pixel of the widget, so allow the rounding)
        cx = _qpoint(cv, x_mpl, y_mpl)
        assert abs(qx[0] - x_mpl) < 1.5 and abs(qy[0] - y_mpl) < 1.5, \
            (qx, qy, x_mpl, y_mpl, x, y, cx)
        _wheel_at(app, cv, x_mpl, y_mpl, +2)
        assert np.allclose(_ranges(ctl), r0 / v3.WHEEL_FACTOR ** 3)
        _wheel_at(app, cv, x_mpl, y_mpl, -3)
        assert np.allclose(_ranges(ctl), r0)
        # right drag pans: same ranges and angles, other centre
        e0, a0 = ctl.ax.elev, ctl.ax.azim
        lo0, hi0 = ctl._limits()
        mx, my = reg.x0 + 0.5 * reg.width, reg.y0 + 0.5 * reg.height
        _drag(app, cv, mx, my, mx + 60.0, my - 20.0, Qt.MouseButton.RightButton)
        lo1, hi1 = ctl._limits()
        assert np.allclose(hi1 - lo1, hi0 - lo0)
        assert np.linalg.norm(lo1 - lo0) > 1e-3
        assert (ctl.ax.elev, ctl.ax.azim) == (e0, a0)
        assert win.selection.is_empty()            # a drag is not a click
        # middle drag pans too
        _drag(app, cv, mx, my, mx - 40.0, my, Qt.MouseButton.MiddleButton)
        lo2, hi2 = ctl._limits()
        assert np.allclose(hi2 - lo2, hi0 - lo0)
        assert np.linalg.norm(lo2 - lo1) > 1e-3
        # left drag rotates and keeps the limits
        _drag(app, cv, mx, my, mx + 80.0, my + 30.0, Qt.MouseButton.LeftButton)
        assert abs(ctl.ax.azim - a0) > 5.0 and abs(ctl.ax.elev - e0) > 1.0
        assert np.allclose(ctl._limits()[0], lo2)
        assert win.selection.is_empty()
        # the toolbar's pan mode: the viewer does not rotate or select
        c.toolbar.pan()
        e1, a1 = ctl.ax.elev, ctl.ax.azim
        _drag(app, cv, mx, my, mx + 50.0, my, Qt.MouseButton.LeftButton)
        _click(app, cv, mx, my)
        assert (ctl.ax.elev, ctl.ax.azim) == (e1, a1)
        assert win.selection.is_empty()
        c.toolbar.pan()
        assert not errors, errors[-1]
        _no_view_errors(win)
    _close(app, win)


def test_view3d_keys_view_bar_find_and_hide():
    """Arrow keys orbit, Home resets, 1 to 6 and the bar's preset buttons
    set the view, + and - zoom, Find selects and frames, Esc clears, the
    toggles (L, M, N, A, O) work, H hides, Shift + H shows, I isolates and
    Ctrl + C copies."""
    from qtpy import QtWidgets
    from qtpy.QtCore import Qt
    from geodome.gui import view3d as v3
    from geodome.model import default_model
    app, win = _shown_window()
    with _SlotErrors() as errors:
        win.set_model(default_model(6, 24.0, 0.5))
        c = win.canvas_model
        cv = c.canvas
        ctl = c.view3d
        bar = ctl.bar
        _draw(app, c)
        e0, a0 = ctl.ax.elev, ctl.ax.azim
        home = ctl.camera()
        _key(app, cv, Qt.Key.Key_Left)
        assert abs(ctl.ax.azim - (a0 - 5.0)) < 1e-9
        _key(app, cv, Qt.Key.Key_Up)
        assert abs(ctl.ax.elev - (e0 + 5.0)) < 1e-9
        _key(app, cv, Qt.Key.Key_Right, shift=True)
        assert abs(ctl.ax.azim - (a0 + 10.0)) < 1e-9
        _key(app, cv, Qt.Key.Key_E)
        assert abs(ctl._roll() - 5.0) < 1e-9
        r0 = _ranges(ctl)
        _key(app, cv, Qt.Key.Key_Plus)
        assert np.allclose(_ranges(ctl), r0 / v3.KEY_ZOOM)
        _key(app, cv, Qt.Key.Key_Minus)
        assert np.allclose(_ranges(ctl), r0)
        _key(app, cv, Qt.Key.Key_Left, ctrl=True)              # pan
        assert np.allclose(_ranges(ctl), r0)
        _key(app, cv, Qt.Key.Key_Home)
        cam = ctl.camera()
        for k in ("elev", "azim", "roll", "xlim", "ylim", "zlim"):
            assert np.allclose(cam[k], home[k]), k
        _key(app, cv, Qt.Key.Key_2)
        assert ctl.ax.elev == 90.0
        _key(app, cv, Qt.Key.Key_3)
        assert (ctl.ax.elev, ctl.ax.azim) == (0.0, -90.0)
        _key(app, cv, Qt.Key.Key_R)
        assert ctl.ax.elev == home["elev"]
        bar.preset_actions["Top"].trigger()
        assert ctl.ax.elev == 90.0
        bar.preset_actions["Side"].trigger()
        assert (ctl.ax.elev, ctl.ax.azim) == (0.0, 0.0)
        bar.reset_act.trigger()
        bar.zoom_in_act.trigger()
        assert np.allclose(_ranges(ctl), r0 / v3.KEY_ZOOM)
        bar.zoom_out_act.trigger()
        assert np.allclose(_ranges(ctl), r0)
        # Find: m 12 selects and frames member 12; n 5 a node
        bar.find_edit.setText("m 12")
        bar.find_edit.returnPressed.emit()
        assert win.selection.members == {12} and not win.selection.nodes
        assert np.all(_ranges(ctl) < r0)
        bar.find_edit.setText("n 5")
        bar.find_edit.returnPressed.emit()
        assert win.selection.nodes == {5} and not win.selection.members
        bar.find_kind.setCurrentIndex(1)                  # plain = node
        bar.find_edit.setText("7")
        bar.find_edit.returnPressed.emit()
        assert win.selection.nodes == {7}
        bar.find_edit.setText("m 99999")
        bar.find_edit.returnPressed.emit()
        assert win.selection.nodes == {7}
        assert "There is no member 99999" in bar.msg.text()
        _key(app, cv, Qt.Key.Key_Escape)
        assert win.selection.is_empty()
        # F with nothing selected fits everything
        _key(app, cv, Qt.Key.Key_F)
        assert np.all(_ranges(ctl) > 0.5 * r0)
        # toggles
        _key(app, cv, Qt.Key.Key_M)
        _key(app, cv, Qt.Key.Key_N)
        assert ctl.member_labels and ctl.node_labels
        assert bar.member_act.isChecked() and win.chk_labels.isChecked()
        assert ctl._label_artists
        assert 0 < len(ctl._label_artists) <= v3.LABEL_CAP
        _key(app, cv, Qt.Key.Key_A)
        assert not ctl.axes_on and not ctl.ax._axis3don
        _key(app, cv, Qt.Key.Key_O)
        assert ctl.ortho and ctl.ax._focal_length == np.inf
        # kept when the view is drawn again
        win.chk_nodes.setChecked(False)
        app.processEvents()
        assert ctl.ortho and not ctl.ax._axis3don and ctl._label_artists
        for key in (Qt.Key.Key_M, Qt.Key.Key_N, Qt.Key.Key_A, Qt.Key.Key_P):
            _key(app, cv, key)
        assert not (ctl.member_labels or ctl.node_labels or ctl.ortho)
        assert ctl.axes_on and not ctl._label_artists
        # H hides the selected member, Shift + H shows it, I isolates
        bar.find_edit.setText("m 3")
        bar.find_kind.setCurrentIndex(0)
        bar.find_edit.returnPressed.emit()
        _key(app, cv, Qt.Key.Key_H)
        assert ctl.hidden is not None and ctl.hidden[3]
        lines = [a for a in ctl.ax.collections
                 if a.get_gid() == "geodome-members"][0]
        assert lines.get_linewidths()[3] == 0.0
        assert 3 not in win.selection.members
        assert not ctl.find("m 3") or ctl.hidden is None   # found: shown
        _key(app, cv, Qt.Key.Key_H, shift=True)
        assert ctl.hidden is None
        win.selection.set([3, 4])
        _key(app, cv, Qt.Key.Key_I)
        assert int((~ctl.hidden).sum()) == 2
        _key(app, cv, Qt.Key.Key_H, shift=True)
        # Ctrl + C copies the details
        win.selection.set([3])
        _key(app, cv, Qt.Key.Key_C, ctrl=True)
        clip = QtWidgets.QApplication.clipboard().text()
        assert "Member 3" in clip and "Length" in clip
        # the help window lists the controls (it does not block)
        dlg = bar.show_help()
        assert "Mouse wheel" in dlg.browser.toPlainText()
        dlg.close()
        # the window shortcuts still work while the drawing has the keys
        win.activateWindow()
        cv.setFocus()
        app.processEvents()
        if win.isActiveWindow() and \
                QtWidgets.QApplication.focusWidget() is cv:
            hits = []
            win.act_generate_menu.triggered.connect(
                lambda *a: hits.append("F5"))
            _key(app, cv, Qt.Key.Key_F5)
            assert hits == ["F5"]
        assert not errors, errors[-1]
        _no_view_errors(win)
    _close(app, win)


def test_view3d_camera_kept_on_redraw_and_reset_on_a_new_model():
    """The camera survives a redraw of the same model (check box, case,
    magnification) and goes back to the first view for a new model."""
    from qtpy.QtCore import Qt
    from geodome.model import default_model
    app, win = _shown_window()
    with _SlotErrors() as errors, _Dialogs():
        win.set_model(default_model(4, 12.0, 0.5))
        c = win.canvas_model
        ctl = c.view3d
        _draw(app, c)
        home = ctl.camera()
        for key in (Qt.Key.Key_Left, Qt.Key.Key_Left, Qt.Key.Key_Up,
                    Qt.Key.Key_Plus):
            _key(app, c.canvas, key)
        cam = ctl.camera()
        assert cam["azim"] != home["azim"]
        win.chk_nodes.setChecked(False)              # redraw, same model
        app.processEvents()
        assert ctl.ax is c.axes
        for k in ("elev", "azim", "xlim", "ylim", "zlim"):
            assert np.allclose(ctl.camera()[k], cam[k]), k
        # results views: change of case and of magnification
        _analyse(app, win)
        win.view_tabs.setCurrentWidget(win.canvas_axial)
        ca = win.canvas_axial
        _draw(app, ca)
        _key(app, ca.canvas, Qt.Key.Key_Right)
        _key(app, ca.canvas, Qt.Key.Key_Minus)
        cam_a = ca.view3d.camera()
        win.axial_combo.setCurrentIndex((win.axial_combo.currentIndex() + 1)
                                        % win.axial_combo.count())
        app.processEvents()
        for k in ("elev", "azim", "xlim", "ylim", "zlim"):
            assert np.allclose(ca.view3d.camera()[k], cam_a[k]), k
        win.view_tabs.setCurrentWidget(win.canvas_defl)
        cd = win.canvas_defl
        _draw(app, cd)
        _key(app, cd.canvas, Qt.Key.Key_Plus)
        cam_d = cd.view3d.camera()
        win.defl_scale.setValue(win.defl_scale.value() * 2.0)
        app.processEvents()
        for k in ("elev", "azim", "xlim", "ylim", "zlim"):
            assert np.allclose(cd.view3d.camera()[k], cam_d[k]), k
        # re-analysis (same model) keeps it too
        _analyse(app, win)
        for k in ("elev", "azim", "xlim"):
            assert np.allclose(ca.view3d.camera()[k], cam_a[k]), k
        # a new model: back to the first view
        win.set_model(default_model(5, 15.0, 0.5))
        app.processEvents()
        _draw(app, c)
        new = ctl.camera()
        assert (new["elev"], new["azim"]) == (22.0, -58.0)
        assert new["xlim"][1] > home["xlim"][1]          # the bigger dome
        assert not ctl._user_view and not ctl._user_limits
        assert not errors, errors[-1]
        _no_view_errors(win)
    _close(app, win)


def test_view3d_modes_animation_and_nonlinear_with_the_viewer():
    """The Modes view keeps animating with the viewer on (wheel zoom and
    keys while it runs), and the Nonlinear 3D shape zooms with the wheel
    while the 2D path next to it does not."""
    from qtpy.QtCore import Qt
    from geodome.gui import nonlinear_view as nv
    from geodome.gui import view3d as v3
    app, win = _shown_window()
    md = win.modal_extension
    nl = win.nonlinear_extension
    with _SlotErrors() as errors, _Dialogs():
        _small_dome(app, win, 2, 6.0)
        assert md.run_modal({"n_modes": 4})
        _wait(app, win)
        mc = md.view
        ctl = mc.view3d
        assert not ctl.picking and ctl.bar is not None
        _draw(app, mc)
        assert md.start_animation()
        r0 = _ranges(ctl)
        reg = ctl._region_px()
        x, y = reg.x0 + 0.5 * reg.width, reg.y0 + 0.5 * reg.height
        f0 = md.frames
        _wheel_at(app, mc.canvas, x, y, +2)
        _key(app, mc.canvas, Qt.Key.Key_Left)
        # no numbers here: they would stay put while the struts move
        _key(app, mc.canvas, Qt.Key.Key_M)
        _key(app, mc.canvas, Qt.Key.Key_N)
        assert not (ctl.member_labels or ctl.node_labels)
        assert not ctl._label_artists
        t0 = time.time()
        while time.time() - t0 < 0.4:
            app.processEvents()
            time.sleep(0.01)
        assert md.is_animating() and md.frames > f0
        assert np.allclose(_ranges(ctl), r0 / v3.WHEEL_FACTOR ** 2)
        # another mode keeps the zoom and the angle
        cam = ctl.camera()
        md.mode_combo.setCurrentIndex(1)
        app.processEvents()
        assert np.allclose(ctl.camera()["xlim"], cam["xlim"])
        assert abs(ctl.camera()["azim"] - cam["azim"]) < 1e-9
        md.stop_animation()
        assert win.selection.is_empty()            # no picking here
        # nonlinear: wheel over the 3D shape zooms it, over the path not
        assert win.load_preset("star_dome")
        s = nv.NonlinearDialog(win, win.model).values()
        assert nl.start_nonlinear(s, synchronous=True)
        win.view_tabs.setCurrentWidget(nl.view)
        app.processEvents()
        pc = nl.view.path_canvas
        _draw(app, pc)
        ctn = pc.view3d
        assert ctn.ax is pc.ax_3d
        r0 = _ranges(ctn)
        b3 = pc.ax_3d.bbox
        _wheel_at(app, pc.canvas, b3.x0 + 0.5 * b3.width,
                  b3.y0 + 0.5 * b3.height, +1)
        assert np.allclose(_ranges(ctn), r0 / v3.WHEEL_FACTOR)
        bp = pc.ax_path.bbox
        _wheel_at(app, pc.canvas, bp.x0 + 0.5 * bp.width,
                  bp.y0 + 0.5 * bp.height, +1)
        assert np.allclose(_ranges(ctn), r0 / v3.WHEEL_FACTOR)
        # a step of the slider redraws and keeps the zoom
        nl.view.slider.setValue(max(0, nl.view.slider.value() - 1))
        for _ in range(20):
            app.processEvents()
            time.sleep(0.01)
        assert ctn.ax is pc.ax_3d
        assert np.allclose(_ranges(ctn), r0 / v3.WHEEL_FACTOR)
        assert not nl.errors, nl.errors[-1]
        assert not errors, errors[-1]
        _no_view_errors(win)
    _close(app, win)


def _snap_nodes(ctl, count=3, min_len=14.0, min_gap=20.0):
    """Front nodes on screen whose members are all at least ``min_len`` px
    long and whose nearest other node is at least ``min_gap`` px away."""
    from geodome.gui import view3d as v3
    sc = ctl.scene
    px, py, pz = v3.project(ctl.ax, sc.nodes)
    reg = ctl._region_px()
    i, j = sc.members[:, 0], sc.members[:, 1]
    L = np.hypot(px[j] - px[i], py[j] - py[i])
    out = []
    for n in np.argsort(pz):
        n = int(n)
        if not reg.contains(px[n] - 10.0, py[n] - 10.0) or \
                not reg.contains(px[n] + 10.0, py[n] + 10.0):
            continue
        at = (sc.members == n).any(axis=1)
        if not at.any() or L[at].min() < min_len:
            continue
        if np.sort(np.hypot(px - px[n], py - py[n]))[1] < min_gap:
            continue
        out.append((n, float(px[n]), float(py[n])))
        if len(out) == count:
            break
    assert out, "no free node on screen"
    return out


def test_view3d_node_click_snaps_and_uses_the_exact_position():
    """A click 2 to 3 px off a node selects the node (every member that
    meets there passes through it, so a node used to need a click on its
    exact centre), in the Model, Deflected shape and Buckling mode views.
    A click on the middle of a member still selects the member.  The click
    position is the exact one: matplotlib rounds it down to a whole pixel,
    which moved it by one pixel in the views whose figure height is not a
    whole number of pixels."""
    from qtpy import QtCore
    from qtpy.QtCore import Qt
    from qtpy.QtTest import QTest
    app, win = _shown_window()
    with _SlotErrors() as errors, _Dialogs():
        assert win.load_preset("case_study_6v")
        _analyse(app, win)
        for c in (win.canvas_model, win.canvas_defl, win.canvas_buck):
            win.view_tabs.setCurrentWidget(c)
            _draw(app, c)
            ctl = c.view3d
            for n, x, y in _snap_nodes(ctl):
                for dx, dy in ((2.5, 0.0), (0.0, -2.5), (-2.0, 2.0)):
                    win.selection.clear()
                    _click(app, c.canvas, x + dx, y + dy)
                    assert win.selection.nodes == {n} and \
                        not win.selection.members, \
                        (ctl.name, n, dx, dy, win.selection.members)
                assert ctl.pick(x + 2.0, y)[:2] == ("node", n)
            k, x, y = _isolated_member(ctl)
            _click(app, c.canvas, x, y)
            assert win.selection.members == {k}, (ctl.name, k)
            # the position the viewer reads is the exact widget position
            seen = []
            cid = c.canvas.mpl_connect(
                "button_press_event",
                lambda e, ctl=ctl: seen.append(ctl._xy(e)))
            q = QtCore.QPoint(int(x), int(c.figure.bbox.height - y))
            QTest.mousePress(c.canvas, Qt.MouseButton.LeftButton,
                             Qt.KeyboardModifier.NoModifier, q)
            QTest.mouseRelease(c.canvas, Qt.MouseButton.LeftButton,
                               Qt.KeyboardModifier.NoModifier, q)
            c.canvas.mpl_disconnect(cid)
            app.processEvents()
            r = float(c.canvas.device_pixel_ratio or 1.0)
            want = (q.x() * r, c.figure.bbox.height - q.y() * r)
            assert np.allclose(seen[0], want, atol=1e-6), (seen, want)
        assert not errors, errors[-1]
        _no_view_errors(win)
    _close(app, win)


def test_view3d_details_follow_the_results():
    """The Model view's details show the results after an analysis (KL/r,
    utilisation, N in every combination, the node displacement) and lose
    them when an input is edited.  An input edit keeps the selection.  The
    emptied result views describe nothing (not the old model), keep the
    viewer's mouse (a right drag does not zoom) and keep the camera for
    the next analysis.  A new model clears all of it."""
    from qtpy.QtCore import Qt
    from geodome.model import default_model
    app, win = _shown_window()
    with _SlotErrors() as errors, _Dialogs():
        assert win.load_preset("case_study_6v")
        cm = win.canvas_model.view3d
        _draw(app, win.canvas_model)
        win.selection.set([12])
        names = [a for a, _ in cm.selection_details()[1]]
        assert "Length" in names and "KL/r" not in names
        _analyse(app, win)
        res = win.results
        win.view_tabs.setCurrentWidget(win.canvas_model)
        _draw(app, win.canvas_model)
        win.selection.set([12])
        rows = dict(cm.selection_details()[1])
        assert "KL/r" in rows and "Utilisation" in rows, rows
        for name in res.order:
            assert rows[f"N  {name}"] == \
                f"{res.cases[name].N[12] / 1e3:+,.2f} kN"
        assert "KL/r" in cm.info_panel.text.toPlainText()
        win.selection.set((), [5])
        assert any(a.startswith("u (") for a, _ in cm.selection_details()[1])
        # move the camera of the Axial force view
        ca = win.canvas_axial
        win.view_tabs.setCurrentWidget(ca)
        _draw(app, ca)
        _key(app, ca.canvas, Qt.Key.Key_Right)
        _key(app, ca.canvas, Qt.Key.Key_Plus)
        cam_a = ca.view3d.camera()
        # an input edit: the results go, the selection stays
        win.member_table.selectRow(12)
        app.processEvents()
        assert win.selection.members == {12}
        win._inputs_changed()
        app.processEvents()
        assert win.results is None
        assert win.selection.members == {12}
        for c in (win.canvas_model, win.canvas_loads):
            names = [a for a, _ in c.view3d.selection_details()[1]]
            assert "Length" in names and "KL/r" not in names, names
            assert not any(a.startswith("N  ") for a in names), names
        for c in (win.canvas_defl, win.canvas_axial, win.canvas_util,
                  win.canvas_buck):
            ctl = c.view3d
            title, rows, summary = ctl.selection_details()
            assert (title, rows, summary) == ("", [], "Nothing drawn here.")
            assert ctl.ax is c.axes and ctl.scene is None
            assert "Nothing drawn here" in ctl.info_panel.text.toPlainText()
        # the emptied view keeps the viewer's mouse: a right drag does not
        # zoom (matplotlib's own 3D mouse zooms with it)
        win.view_tabs.setCurrentWidget(ca)
        _draw(app, ca)
        reg = ca.view3d._region_px()
        mx, my = reg.x0 + 0.5 * reg.width, reg.y0 + 0.5 * reg.height
        lim0 = (ca.axes.get_xlim(), ca.axes.get_ylim(), ca.axes.get_zlim())
        _drag(app, ca.canvas, mx, my, mx + 60.0, my - 40.0,
              Qt.MouseButton.RightButton)
        lim1 = (ca.axes.get_xlim(), ca.axes.get_ylim(), ca.axes.get_zlim())
        assert np.allclose(lim0, lim1), (lim0, lim1)
        # the next analysis of the same model brings the camera back
        _analyse(app, win)
        win.view_tabs.setCurrentWidget(ca)
        _draw(app, ca)
        for k in ("elev", "azim", "xlim", "ylim", "zlim"):
            assert np.allclose(ca.view3d.camera()[k], cam_a[k]), k
        assert win.selection.members == {12}
        # a new model: the selection goes and nothing old is described
        win.set_model(default_model(4, 12.0, 0.5))
        app.processEvents()
        assert win.selection.is_empty()
        win.selection.set([3])
        i, j = win.model.geom.members[3]
        assert dict(cm.selection_details()[1])["End nodes"] == f"{i} - {j}"
        for c in (win.canvas_defl, win.canvas_axial, win.canvas_util,
                  win.canvas_buck):
            assert c.view3d.selection_details()[2] == "Nothing drawn here."
        assert not errors, errors[-1]
        _no_view_errors(win)
    _close(app, win)


def test_view3d_legend_of_many_families_titles_and_the_z_axis():
    """A custom truss with about 200 strut families: the legend shows the
    rows that fit, ends with '+ N more', stays inside the figure and
    within its width limit, and left of the model, at three sizes.  The
    z axis labels of the report's model figure stay inside the figure, and
    a two-line title stays inside a short figure."""
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure
    from geodome import loads as ld
    from geodome import plotting as plot
    from geodome import report as rp
    from geodome.gui import view3d as v3
    from geodome.model import DomeModel, default_model

    def figure(w, h):
        fig = Figure(figsize=(w, h), dpi=100)
        FigureCanvasAgg(fig)
        ax = fig.add_subplot(111, projection="3d")
        plot.axes_layout(ax)
        return fig, ax

    g = geo.generate_dome(4, 12.0, 0.5)
    rng = np.random.default_rng(7)
    nodes = g.nodes + rng.normal(0.0, 0.05, g.nodes.shape)
    cg = geo.custom_geometry(nodes, g.members, faces=g.faces,
                             base_nodes=g.base_nodes)
    model = DomeModel(geom=cg)
    n_fam = cg.n_strut_types()
    assert n_fam > 150, n_fam
    for w, h in ((9.42, 4.78), (5.0, 2.6), (3.8, 2.0)):
        fig, ax = figure(w, h)
        plot.draw_model(ax, model)
        fig.canvas.draw()
        r = fig.canvas.get_renderer()
        leg = ax.get_legend()
        lb = leg.get_window_extent(r)
        W, H = fig.bbox.width, fig.bbox.height
        assert lb.y0 >= -0.5 and lb.y1 <= H + 0.5, (w, h, lb.extents)
        assert lb.x0 >= 0.0 and lb.width <= plot._LEGEND_MAX_WIDTH * W + 1.0
        lay = ax.get_axes_locator()
        texts = [t.get_text() for t in leg.get_texts()]
        assert 1 <= lay.n_shown < n_fam
        assert len(texts) == lay.n_shown + 1
        assert texts[-1] == f"+ {n_fam - lay.n_shown} more", texts[-1]
        px, _, _ = v3.project(ax, nodes)
        assert lb.x1 < px.min()
    # the static export of the report: the z tick labels and the z label
    # stay inside the figure (the legend column took the width)
    for f in (1, 6, 10):
        fig = rp._fig3d(plot.draw_model, default_model(f, 24.0, 0.5))
        FigureCanvasAgg(fig)
        fig.canvas.draw()
        r = fig.canvas.get_renderer()
        ax = fig.axes[0]
        assert ax.zaxis.label.get_text() == "z [m]"
        zb = ax.zaxis.get_tightbbox(r)
        assert zb.x1 <= fig.bbox.width, (f, zb.x1, fig.bbox.width)
        assert ax.zaxis.label.get_window_extent(r).x1 <= fig.bbox.width
    # a two-line title (a combination with temperature) fits
    m6 = default_model(6, 24.0, 0.5)
    F = ld.combine(ld.basic_load_vectors(m6), m6.combinations[0])
    for h in (4.8, 3.2):
        fig, ax = figure(9.42, h)
        plot.draw_loads(ax, m6, F, title="Applied nodal loads -- LC1\n"
                                         "(+ temperature / settlement, "
                                         "not drawn)")
        fig.canvas.draw()
        tb = ax.title.get_window_extent(fig.canvas.get_renderer())
        assert tb.y1 <= fig.bbox.height + 0.5, (h, tb.extents)


def test_view3d_wheel_toolbar_camera_exports_labels_and_hover_limit():
    """One wheel event of five notches zooms four notches (not almost
    nothing), a touchpad's small delta zooms a fraction of a notch; a zoom
    made with the matplotlib toolbar survives a redraw of the same model;
    saving the figure as PDF, SVG or PNG with the hover highlight on gives
    no error; with the numbers on, a settled view that shows the same
    numbers is not drawn again; above the hover size limit the bar says
    that hover is off."""
    import tempfile
    from qtpy.QtCore import Qt
    from geodome.gui import view3d as v3
    from geodome.model import default_model
    app, win = _shown_window()
    with _SlotErrors() as errors:
        win.set_model(default_model(6, 24.0, 0.5))
        c = win.canvas_model
        cv = c.canvas
        ctl = c.view3d
        _draw(app, c)
        reg = ctl._region_px()
        x, y = reg.x0 + 0.6 * reg.width, reg.y0 + 0.5 * reg.height
        r0 = _ranges(ctl)
        _wheel_at(app, cv, x, y, angle=600)
        assert np.allclose(_ranges(ctl), r0 / v3.WHEEL_FACTOR ** 4)
        # zoomed in past half the first view: a hint names the A key
        assert "Press A" in ctl.bar.msg.text(), ctl.bar.msg.text()
        _key(app, cv, Qt.Key.Key_Home)
        _wheel_at(app, cv, x, y, angle=30, pixel=4)
        assert np.allclose(_ranges(ctl), r0 / v3.WHEEL_FACTOR ** 0.25)
        _key(app, cv, Qt.Key.Key_Home)
        assert np.allclose(_ranges(ctl), r0)
        # the toolbar's zoom rectangle, then a redraw of the same model
        c.toolbar.zoom()
        _drag(app, cv, x - 90.0, y - 60.0, x + 90.0, y + 60.0,
              Qt.MouseButton.LeftButton)
        c.toolbar.zoom()
        app.processEvents()
        r_tb = _ranges(ctl)
        assert np.all(r_tb < 0.9 * r0), (r_tb, r0)
        win.chk_nodes.setChecked(False)
        app.processEvents()
        assert np.allclose(_ranges(ctl), r_tb)
        # the toolbar's Home goes back, and a redraw keeps that too
        c.toolbar.home()
        app.processEvents()
        _draw(app, c)
        win.chk_nodes.setChecked(True)
        app.processEvents()
        assert np.allclose(_ranges(ctl), r0)
        # figure exports with the hover highlight on
        _draw(app, c)
        k, mx, my = _isolated_member(ctl)
        _move(app, cv, mx, my)
        assert ctl._hover_artist.get_visible()
        with tempfile.TemporaryDirectory() as tmp:
            for ext in ("pdf", "svg", "png"):
                c.figure.savefig(os.path.join(tmp, f"model.{ext}"))
                assert os.path.getsize(os.path.join(tmp, f"model.{ext}"))
        assert not ctl.errors, ctl.errors[-1]
        # numbers on: a settled view with the same numbers is not redrawn
        _key(app, cv, Qt.Key.Key_M)
        _draw(app, c)
        draws = []
        cid = cv.mpl_connect("draw_event", lambda e: draws.append(1))
        ctl._on_settled()
        app.processEvents()
        assert not draws
        _wheel_at(app, cv, x, y, +2)              # other numbers show
        ctl._on_settled()
        app.processEvents()
        assert draws
        cv.mpl_disconnect(cid)
        _key(app, cv, Qt.Key.Key_M)
        # above the hover limit the bar says so (once per model)
        old = v3.HOVER_MAX_MEMBERS
        v3.HOVER_MAX_MEMBERS = 10
        try:
            win.set_model(default_model(2, 6.0, 0.5))
            assert ctl._hover_artist is None
            assert "Hover is off" in ctl.bar.msg.text()
        finally:
            v3.HOVER_MAX_MEMBERS = old
        win.redraw_model()
        assert ctl._hover_artist is not None
        assert not errors, errors[-1]
        _no_view_errors(win)
    _close(app, win)


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
