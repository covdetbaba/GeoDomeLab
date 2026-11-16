"""
GUI tests for the inputs, models and presets (GeoDomeLab 2.0, task H1).

Checks, through the real main window running off screen:

* the menu bar and menus stay light with dark text under a dark palette
  (Windows dark mode);
* the Examples menu lists every preset of ``benchmarks.PRESETS``; each one
  loads, refreshes the panels, runs, and still matches its reference answer
  after a round trip through the panels;
* custom trusses: read-only geodesic inputs, save / open without
  regeneration, Generate dome returns to a geodesic dome;
* elastic supports, temperature, settlements, nodal loads and the T column
  reach the model and the analysis exactly;
* bad input gives a clear message and runs nothing;
* File > Import truss from CSV;
* the About box, the design wall thickness factor, ``model_changed()`` and
  ``extensions``;
* no exception escapes a Qt slot (the ``_SlotErrors`` pattern of
  tests/test_gui_generate.py).

Run with pytest::

    pytest -q tests/test_gui_inputs.py

or standalone::

    python tests/test_gui_inputs.py
"""

from __future__ import annotations

import json
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

import geodome                                                      # noqa: E402
from geodome import benchmarks as bm                                # noqa: E402
from geodome import csv_import                                      # noqa: E402
from geodome import geometry as geo                                 # noqa: E402
from geodome.model import DomeModel, LoadCombination                # noqa: E402
from geodome.postprocess import run_analysis                        # noqa: E402


# ================================================================ helpers ==
def _app():
    from qtpy import QtWidgets
    return QtWidgets.QApplication.instance() or QtWidgets.QApplication([])


def _window():
    app = _app()
    from geodome.gui.main_window import MainWindow
    win = MainWindow()
    win.resize(1600, 1000)
    return app, win


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


class _Dialogs:
    """Replace the blocking QMessageBox / QFileDialog calls.

    Message boxes are recorded in ``self.messages`` as (kind, title, text).
    ``open_path`` / ``save_path`` are returned by the file dialogs.
    """

    def __init__(self, open_path="", save_path="", answer=None):
        self.open_path = open_path
        self.save_path = save_path
        self.answer = answer
        self.messages = []

    def __enter__(self):
        from qtpy import QtWidgets
        mb, fd = QtWidgets.QMessageBox, QtWidgets.QFileDialog
        self._saved = {n: getattr(mb, n) for n in
                       ("warning", "critical", "information", "question",
                        "about")}
        self._saved_fd = {n: getattr(fd, n) for n in
                          ("getOpenFileName", "getSaveFileName")}

        def rec(kind):
            def f(parent, title, text, *a, **k):
                self.messages.append((kind, str(title), str(text)))
                return self.answer if self.answer is not None \
                    else QtWidgets.QMessageBox.Ok
            return staticmethod(f)

        for n in ("warning", "critical", "information", "about"):
            setattr(mb, n, rec(n))
        no = QtWidgets.QMessageBox.No
        mb.question = staticmethod(
            lambda parent, title, text, *a, **k:
            self.messages.append(("question", title, text)) or no)
        fd.getOpenFileName = staticmethod(
            lambda *a, **k: (self.open_path, ""))
        fd.getSaveFileName = staticmethod(
            lambda *a, **k: (self.save_path, ""))
        return self

    def __exit__(self, *exc):
        from qtpy import QtWidgets
        for n, f in self._saved.items():
            setattr(QtWidgets.QMessageBox, n, f)
        for n, f in self._saved_fd.items():
            setattr(QtWidgets.QFileDialog, n, f)
        return False

    def texts(self, kind=None):
        return [t for k, _, t in self.messages if kind is None or k == kind]


def _run(win, app, timeout=180.0):
    """Click Run analysis and wait for the worker; return the results."""
    win.on_analyse()
    t0 = time.time()
    while win._busy and time.time() - t0 < timeout:
        if win._worker is not None:
            win._worker.wait(25)
        app.processEvents()
    app.processEvents()
    assert not win._busy, "analysis did not finish"
    return win.results


def _small_dome(win, app, freq=2, radius=6.0):
    gp = win.geometry_panel
    gp.frequency.setValue(freq)
    gp.radius.setValue(radius)
    gp.fraction.setValue(0.5)
    gp.generate_btn.click()
    app.processEvents()


def _set_cell(table, row, col, text):
    from qtpy import QtWidgets
    item = table.item(row, col)
    if item is None:
        item = QtWidgets.QTableWidgetItem()
        table.setItem(row, col, item)
    item.setText(str(text))


def _luminance(img):
    vals = []
    for x in range(0, img.width(), 2):
        for y in range(0, img.height(), 2):
            c = img.pixelColor(x, y)
            vals.append(0.299 * c.red() + 0.587 * c.green() + 0.114 * c.blue())
    return np.asarray(vals)


def _dark_palette():
    from qtpy import QtGui
    pal = QtGui.QPalette()
    role = getattr(QtGui.QPalette, "ColorRole", QtGui.QPalette)
    for name, col in (("Window", "#202020"), ("WindowText", "#f0f0f0"),
                      ("Base", "#1b1b1b"), ("Text", "#f0f0f0"),
                      ("Button", "#2d2d2d"), ("ButtonText", "#f0f0f0"),
                      ("ToolTipBase", "#2b2b2b"), ("ToolTipText", "#f0f0f0")):
        pal.setColor(getattr(role, name), QtGui.QColor(col))
    return pal


# ================================================================== tests ==
def test_menu_bar_is_light_in_dark_mode():
    """Item 1: light menu bar, menus and tool bar with dark text."""
    from qtpy import QtCore
    from geodome.gui import style as st
    for rule in ("QMenuBar {", "QMenuBar::item", "QMenuBar::item:selected",
                 "QMenu {", "QMenu::item:selected", "QToolBar {",
                 "QToolTip {", "QMessageBox"):
        assert rule in st.STYLESHEET, rule
    app = _app()
    app.setPalette(_dark_palette())         # OS in dark mode, before start
    try:
        with _SlotErrors() as errors:
            app, win = _window()
            win.show()
            app.setPalette(_dark_palette())  # ... and switched while running
            app.processEvents()
            mb = win.menuBar()
            lum = _luminance(mb.grab().toImage())
            assert lum.mean() > 225, lum.mean()          # light background
            assert lum.min() < 110, lum.min()            # dark text on it
            tb = _luminance(win.toolbar.grab().toImage())
            assert tb.mean() > 225 and tb.min() < 110
            for title in ("File", "Examples"):
                menu = [a.menu() for a in mb.actions()
                        if a.text().replace("&", "") == title][0]
                menu.popup(win.mapToGlobal(QtCore.QPoint(20, 30)))
                app.processEvents()
                ml = _luminance(menu.grab().toImage())
                menu.hide()
                assert ml.mean() > 215 and ml.min() < 110, (title, ml.mean())
            assert not errors, errors[-1]
            win.close()
    finally:
        st.apply_light_theme(app)


def test_examples_menu_lists_every_preset():
    """Item 2: one action per preset, title + description status tip."""
    with _SlotErrors() as errors:
        app, win = _window()
        menu = win.menu_examples
        assert menu.title().replace("&", "") == "Examples"
        acts = [a for a in menu.actions() if not a.isSeparator()]
        assert [a.data() for a in acts] == bm.preset_keys()
        for a, key in zip(acts, bm.preset_keys()):
            title, desc = bm.preset_info(key)
            assert a.text() == title
            assert a.statusTip() == desc
        assert not errors, errors[-1]
        win.close()


def test_every_preset_loads_runs_and_round_trips():
    """Items 2, 3, 5: every preset loads, refreshes the panels and the
    window title, reads back through the panels unchanged, runs, and still
    matches its independent reference answer."""
    with _SlotErrors() as errors, _Dialogs() as dlg:
        app, win = _window()
        calls = []
        win.extensions.append(type("Rec", (), {
            "on_model_changed": lambda self, m: calls.append(m)})())
        for key in bm.preset_keys():
            title, _ = bm.preset_info(key)
            win.example_actions[key].trigger()
            app.processEvents()
            ref = bm.make_preset(key)
            m = win.model
            g = m.geom
            assert title in win.windowTitle()
            assert calls and calls[-1] is m
            assert f"Example loaded: {title}" in win.log_view.toPlainText()
            gp, lp = win.geometry_panel, win.load_panel
            assert gp.is_custom == g.is_custom
            assert gp.custom_note.isHidden() == (not g.is_custom)
            for w in (gp.frequency, gp.radius, gp.fraction):
                assert w.isEnabled() == (not g.is_custom)
            assert gp.schedule.rowCount() == g.n_strut_types()
            assert len(win.section_panel._combos) == g.n_strut_types()
            assert lp.nodal_table.rowCount() == len(ref.loads.nodal_loads)
            assert lp.settle_table.rowCount() == len(ref.loads.settlements)
            assert lp.combos.rowCount() == len(ref.combinations)
            assert lp.no_faces_note.isHidden() == (g.n_faces > 0)
            assert win.loadcase_combo.count() == len(ref.combinations)
            # the panels read back exactly what the preset holds
            assert win._read_inputs(), dlg.texts()
            assert m.support_type == ref.support_type
            assert m.custom_supports == ref.custom_supports
            assert m.loads.nodal_loads == ref.loads.nodal_loads
            assert m.loads.settlements == ref.loads.settlements
            assert m.loads.temperature_change == ref.loads.temperature_change
            assert [c.factors() for c in m.combinations] == \
                [c.factors() for c in ref.combinations]
            assert [c.serviceability for c in m.combinations] == \
                [c.serviceability for c in ref.combinations]
            np.testing.assert_array_equal(g.nodes, ref.geom.nodes)
            # a run through the worker
            res = _run(win, app)
            # the results keep a snapshot of the model (review GUI-03)
            assert res is not None and res.model is not m
            assert res.model.to_dict() == m.to_dict()
            assert win.member_table.rowCount() == g.n_members
            assert not dlg.messages, dlg.messages
            # and the reference answers still hold for the panel-read model
            rows = bm.verify_preset(key, model=m)
            bad = [r for r in rows if not r["pass"]]
            assert not bad, (key, bad[:2])
            assert not errors, errors[-1]
        win.close()


def test_custom_truss_save_and_open_without_regeneration():
    """Item 6: JSON v2 round trip of a custom truss through File menu."""
    with _SlotErrors() as errors, tempfile.TemporaryDirectory() as d:
        app, win = _window()
        for key in ("three_bar", "star_dome"):
            win.load_preset(key)
            win.design_panel.thickness_factor.setValue(0.93)
            before = win.model
            path = os.path.join(d, f"{key}.json")
            with _Dialogs(save_path=path) as dlg:
                win.on_save()
            assert not dlg.messages, dlg.messages
            with open(path, encoding="utf-8") as fh:
                data = json.load(fh)
            assert data["version"] == 2
            assert data["geometry"]["kind"] == "custom"
            assert data["design"]["design_thickness_factor"] == 0.93
            win.on_new()
            app.processEvents()
            assert not win.model.geom.is_custom
            with _Dialogs(open_path=path) as dlg:
                win.on_open()
            app.processEvents()
            assert not dlg.messages, dlg.messages
            m = win.model
            assert m.geom.is_custom and m.geom.title == before.geom.title
            np.testing.assert_array_equal(m.geom.nodes, before.geom.nodes)
            np.testing.assert_array_equal(m.geom.members,
                                          before.geom.members)
            assert m.custom_supports == before.custom_supports
            assert m.support_type == before.support_type
            assert m.loads.nodal_loads == before.loads.nodal_loads
            assert m.loads.settlements == before.loads.settlements
            assert m.loads.temperature_change == \
                before.loads.temperature_change
            assert [c.T for c in m.combinations] == \
                [c.T for c in before.combinations]
            assert m.design_thickness_factor == 0.93
            assert win.design_panel.thickness_factor.value() == 0.93
            assert win.geometry_panel.is_custom
            assert f"[{key}.json]" in win.windowTitle()
            res = _run(win, app)
            assert res is not None
        assert not errors, errors[-1]
        win.close()


def test_open_bad_file_shows_a_message():
    with _SlotErrors() as errors, tempfile.TemporaryDirectory() as d:
        app, win = _window()
        bad = os.path.join(d, "bad.json")
        with open(bad, "w", encoding="utf-8") as fh:
            fh.write("{not json")
        before = win.model
        with _Dialogs(open_path=bad) as dlg:
            win.on_open()
        assert dlg.texts("critical") and "bad.json" in dlg.texts()[0]
        assert win.model is before
        assert not errors, errors[-1]
        win.close()


def test_generate_dome_returns_from_a_custom_truss():
    """Item 3: Generate dome builds a geodesic dome from the inputs."""
    with _SlotErrors() as errors:
        app, win = _window()
        gp = win.geometry_panel
        gp.frequency.setValue(3)
        gp.radius.setValue(9.0)
        win.load_preset("three_bar")
        assert gp.is_custom and not gp.frequency.isEnabled()
        assert gp.support.currentData() == geo.CUSTOM
        assert not gp.custom_box.isHidden()
        assert "2 free direction" in gp.custom_summary.text()
        assert gp.generate_btn.isEnabled()
        gp.generate_btn.click()
        app.processEvents()
        g = win.model.geom
        assert not g.is_custom and g.frequency == 3 and g.radius == 9.0
        assert not gp.is_custom and gp.frequency.isEnabled()
        assert gp.custom_note.isHidden()
        assert win.model.support_type == geo.PINNED
        # per-node data of the old numbering is dropped, with a log line
        assert win.model.loads.settlements == []
        assert win.model.loads.nodal_loads == []
        assert win.load_panel.settle_table.rowCount() == 0
        assert "node numbers changed" in win.log_view.toPlainText()
        assert win._read_inputs()
        assert _run(win, app) is not None
        assert not errors, errors[-1]
        win.close()


def test_nan_chord_factors_and_no_panels_are_shown_cleanly():
    """Item 3: NaN chord factors read 'n/a'; models without panels say so."""
    with _SlotErrors() as errors:
        app, win = _window()
        win.load_preset("von_mises")
        sched = win.geometry_panel.schedule
        col = [sched.horizontalHeaderItem(j).text()
               for j in range(sched.columnCount())].index("Chord factor")
        assert sched.item(0, col).text() == "n/a"
        assert "nan" not in win.section_panel.mass_label.text().lower()
        assert not win.load_panel.no_faces_note.isHidden()
        win.redraw_model()
        win.redraw_loads()
        win.chk_panels.setChecked(True)
        app.processEvents()
        res = _run(win, app)
        assert res is not None
        win.redraw_charts()
        assert not errors, errors[-1]
        win.close()


def test_elastic_supports_ui_and_run():
    """Item 4: spring inputs shown only for Elastic; they reach the run."""
    with _SlotErrors() as errors, _Dialogs() as dlg:
        app, win = _window()
        _small_dome(win, app)
        gp = win.geometry_panel
        assert gp.spring_box.isHidden() and gp.custom_box.isHidden()
        gp.set_support_type(geo.ELASTIC)
        assert not gp.spring_box.isHidden()
        gp.spring_kx.setValue(20000.0)          # kN/m
        gp.spring_ky.setValue(20000.0)
        gp.spring_kz.setValue(5000.0)
        app.processEvents()
        assert "kx, ky, kz = 20, 20, 5 kN/mm" in gp.info.toPlainText()
        res = _run(win, app)
        assert res is not None and not dlg.messages, dlg.messages
        m = win.model
        assert m.support_type == geo.ELASTIC
        assert m.spring_stiffness == (2.0e7, 2.0e7, 5.0e6)
        # the same model analysed directly gives the same answer
        ref = DomeModel.from_dict(json.loads(json.dumps(m.to_dict())))
        rr = run_analysis(ref, do_buckling=False)
        for name, r in res.cases.items():
            np.testing.assert_allclose(r.u, rr.cases[name].u, rtol=1e-12,
                                       atol=1e-15)
        # the springs carry the load: their forces balance the loads
        r = res.cases[res.service_case]
        _, springs = m.support_conditions()
        assert len(springs) == 3 * len(m.geom.base_nodes)
        assert np.abs(r.spring_forces).max() > 0.0
        # softer than the pinned ring
        gp.set_support_type(geo.PINNED)
        res_p = _run(win, app)
        uz_e = np.abs(r.u[2::3]).max()
        uz_p = np.abs(res_p.cases[res_p.service_case].u[2::3]).max()
        assert uz_e > uz_p
        assert not errors, errors[-1]
        win.close()


def test_loads_tab_t_column_temperature_settlement_nodal_loads():
    """Item 5: the new load inputs reach the model and the analysis."""
    with _SlotErrors() as errors, _Dialogs() as dlg:
        app, win = _window()
        _small_dome(win, app)
        lp = win.load_panel
        m = win.model
        base = int(m.geom.base_nodes[0])
        crown = int(m.geom.crown_node)
        lp.delta_T.setValue(25.0)
        lp.add_settlement(base, 0.0, 0.0, -3.0)          # mm
        lp.add_nodal_load(crown, 2.0, 0.0, -10.0, "L")   # kN
        lp.add_nodal_load(crown, 0.0, 1.5, 0.0, "W")
        # buttons add sensible default rows too
        n0 = lp.settle_table.rowCount()
        lp._add_settlement()
        assert lp.settle_table.rowCount() == n0 + 1
        assert int(lp.settle_table.item(n0, 0).text()) in \
            [int(b) for b in m.geom.base_nodes]
        lp._remove_from(lp.settle_table)
        assert lp.settle_table.rowCount() == n0
        # a T combination typed in the table
        lp._add_row()
        i = lp.combos.rowCount() - 1
        _set_cell(lp.combos, i, 0, "LC9  1.2D + 1.0L + 1.0T")
        for j, v in ((1, "1.2"), (2, "0"), (3, "1"), (4, "0"), (5, "1")):
            _set_cell(lp.combos, i, j, v)
        assert lp.combos.horizontalHeaderItem(5).text() == "T"
        assert win._read_inputs(), dlg.texts()
        ls = m.loads
        assert ls.temperature_change == 25.0
        assert ls.settlements == [{"node": base, "ux": 0.0, "uy": 0.0,
                                   "uz": -0.003}]
        assert ls.nodal_loads == [
            {"node": crown, "Fx": 2000.0, "Fy": 0.0, "Fz": -10000.0,
             "case": "L"},
            {"node": crown, "Fx": 0.0, "Fy": 1500.0, "Fz": 0.0, "case": "W"}]
        c = m.combinations[-1]
        assert (c.name, c.D, c.L, c.T) == ("LC9  1.2D + 1.0L + 1.0T",
                                           1.2, 1.0, 1.0)
        res = _run(win, app)
        assert res is not None and not dlg.messages, dlg.messages
        # same answer as a model built directly from the saved JSON
        ref = DomeModel.from_dict(json.loads(json.dumps(m.to_dict())))
        rr = run_analysis(ref, do_buckling=False)
        r9, q9 = res.cases[c.name], rr.cases[c.name]
        np.testing.assert_allclose(r9.N, q9.N, rtol=1e-12, atol=1e-9)
        # the settlement is imposed exactly and T changes the forces
        assert abs(r9.u[3 * base + 2] - (-0.003)) < 1e-12
        lp.delta_T.setValue(0.0)
        lp.clear_node_tables()
        assert win._read_inputs()
        res0 = _run(win, app)
        assert np.abs(res0.cases[c.name].N - r9.N).max() > 1.0   # N
        # sync_from_model shows the values again (mm, kN)
        m.loads = ref.loads.copy()
        lp.sync_from_model(m)
        assert lp.delta_T.value() == 25.0
        assert lp.settle_table.item(0, 3).text() == "-3"
        assert lp.nodal_table.item(0, 3).text() == "-10"
        assert lp.nodal_table.cellWidget(1, 4).currentText() == "W"
        assert not errors, errors[-1]
        win.close()


def test_values_from_a_model_read_back_exactly():
    """Values the boxes show rounded (springs, dT, forces) are not changed
    by a read of the panels unless the user edits them."""
    with _SlotErrors() as errors:
        app, win = _window()
        _small_dome(win, app)
        m = win.model
        crown = int(m.geom.crown_node)
        base = int(m.geom.base_nodes[0])
        m.support_type = geo.ELASTIC
        m.spring_stiffness = (1.23456789e7, 2.0e7, 3.3333333e6)
        m.loads.temperature_change = 12.3456789
        m.loads.settlements = [{"node": base, "ux": 0.0, "uy": 0.0,
                                "uz": -1.23456789e-3}]
        m.loads.nodal_loads = [{"node": crown, "Fx": 1234.56789, "Fy": 0.0,
                                "Fz": -9876.54321, "case": "S"}]
        m.combinations = [LoadCombination("X", D=1.234567, T=0.3333333)]
        m.design_thickness_factor = 0.93
        win.set_model(m)
        assert win._read_inputs()
        assert m.spring_stiffness == (1.23456789e7, 2.0e7, 3.3333333e6)
        assert m.loads.temperature_change == 12.3456789
        assert m.loads.settlements[0]["uz"] == -1.23456789e-3
        assert m.loads.nodal_loads[0]["Fx"] == 1234.56789
        assert m.loads.nodal_loads[0]["Fz"] == -9876.54321
        assert (m.combinations[0].D, m.combinations[0].T) == \
            (1.234567, 0.3333333)
        # an edit is taken as typed
        win.geometry_panel.spring_kx.setValue(15000.0)
        _set_cell(win.load_panel.nodal_table, 0, 1, "2.5")
        assert win._read_inputs()
        assert m.spring_stiffness[0] == 1.5e7
        assert m.loads.nodal_loads[0]["Fx"] == 2500.0
        assert not errors, errors[-1]
        win.close()


def test_bad_inputs_give_a_clear_message_and_run_nothing():
    """Item 5: node numbers, supported settlement nodes, numbers, names."""
    with _SlotErrors() as errors:
        app, win = _window()
        _small_dome(win, app)
        lp = win.load_panel
        m = win.model
        n = m.geom.n_nodes
        crown = int(m.geom.crown_node)
        base = int(m.geom.base_nodes[0])

        def attempt(expect):
            with _Dialogs() as dlg:
                win.on_analyse()
                app.processEvents()
            assert not win._busy and win._worker is None
            assert dlg.texts("warning"), "no message shown"
            text = dlg.texts("warning")[-1]
            assert expect in text, text
            return text

        lp.add_nodal_load(n + 5, 0, 0, -1)
        attempt(f"node {n + 5} does not exist. Nodes are 0 to {n - 1}")
        lp.clear_node_tables()
        lp.add_nodal_load(crown, 0, 0, -1)
        _set_cell(lp.nodal_table, 0, 0, "1.5")
        attempt("is not a whole number")
        _set_cell(lp.nodal_table, 0, 0, "")
        attempt("Nodal loads, row 1: the node number is missing")
        _set_cell(lp.nodal_table, 0, 0, str(crown))
        _set_cell(lp.nodal_table, 0, 3, "-1,0")
        attempt("'-1,0' is not a number")
        lp.clear_node_tables()
        lp.add_settlement(crown, 0, 0, -2)
        attempt(f"node {crown} is not a supported node")
        lp.clear_node_tables()
        # roller supports: ux of a base node without a lateral restraint
        win.geometry_panel.set_support_type(geo.ROLLER3)
        fixed, _ = m.support_conditions()
        free_x = [int(b) for b in m.geom.base_nodes
                  if 3 * int(b) not in set(fixed.tolist())][0]
        lp.add_settlement(free_x, 1.0, 0, 0)
        attempt(f"node {free_x} is free in x")
        lp.clear_node_tables()
        win.geometry_panel.set_support_type(geo.PINNED)
        _set_cell(lp.combos, 0, 2, "abc")
        attempt("'abc' is not a number")
        lp._reset_rows()
        _set_cell(lp.combos, 1, 0, lp.combos.item(0, 0).text())
        attempt("Give each a different name")
        lp._reset_rows()
        # Custom supports without any data
        win.geometry_panel.set_support_type(geo.CUSTOM)
        attempt("Custom supports need per-node data")
        win.geometry_panel.set_support_type(geo.PINNED)
        # all fixed again: it runs
        lp.add_settlement(base, 0, 0, -1)
        assert win._read_inputs()
        assert _run(win, app) is not None
        assert not errors, errors[-1]
        win.close()


def _write(path, text):
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(text)


def test_csv_import_through_the_file_menu():
    """Item 6: File > Import truss from CSV (prefix, header, comments,
    1-based ids, springs, cases) and a clear message on a bad file."""
    with _SlotErrors() as errors, tempfile.TemporaryDirectory() as d:
        app, win = _window()
        p = os.path.join(d, "tripod_")
        _write(p + "nodes.csv", "# a tripod\nid,x,y,z\n"
               "1, 0, 0, 3\n2, 2, 0, 0\n3, -1, 1.7320508, 0\n"
               "4, -1, -1.7320508, 0\n")
        _write(p + "members.csv", "i,j\n1,2\n1,3\n1,4\n")
        _write(p + "supports.csv", "node,cx,cy,cz\n2,1,1,1\n3,1,1,1\n"
               "4,fixed,1,25000\n")
        _write(p + "loads.csv", "node,Fx,Fy,Fz,case\n1,5,0,-30,L\n"
               "1,0,0,-10,D\n")
        with _Dialogs(open_path=p + "nodes.csv") as dlg:
            win.on_import_csv()
        app.processEvents()
        assert not dlg.messages, dlg.messages
        m = win.model
        g = m.geom
        assert g.is_custom and g.n_nodes == 4 and g.n_members == 3
        assert g.title == "Imported truss (tripod)"
        assert m.support_type == geo.CUSTOM
        assert m.custom_supports == {1: [True, True, True],
                                     2: [True, True, True],
                                     3: [True, True, 25000.0 * 1e3]}
        assert m.loads.nodal_loads == [
            {"node": 0, "Fx": 5000.0, "Fy": 0.0, "Fz": -30000.0, "case": "L"},
            {"node": 0, "Fx": 0.0, "Fy": 0.0, "Fz": -10000.0, "case": "D"}]
        assert win.geometry_panel.is_custom
        assert win.geometry_panel.support.currentData() == geo.CUSTOM
        assert "renumbered" in win.log_view.toPlainText()
        assert "tripod" in win.windowTitle()
        res = _run(win, app)
        assert res is not None
        # equilibrium of the crown under LC2 = 1.2D + 1.6L + 0.5S
        r = res.cases["LC2  1.2D + 1.6L + 0.5S"]
        assert np.abs(r.reaction_balance).max() < 1e-6

        # a member that refers to a missing node: message, model unchanged
        _write(p + "members.csv", "i,j\n1,2\n1,3\n1,9\n")
        before = win.model
        with _Dialogs(open_path=p + "nodes.csv") as dlg:
            win.on_import_csv()
        assert win.model is before
        text = dlg.texts("warning")[0]
        assert "members.csv, line 4" in text and "9" in text
        # no members file at all
        os.remove(p + "members.csv")
        with _Dialogs(open_path=p + "nodes.csv") as dlg:
            win.on_import_csv()
        assert "No tripod_members.csv" in dlg.texts("warning")[0]
        assert win.model is before

        # the example shipped with the program imports and runs
        ex = os.path.join(ROOT, "examples", "csv_truss", "nodes.csv")
        with _Dialogs(open_path=ex) as dlg:
            win.on_import_csv()
        assert not dlg.messages, dlg.messages
        assert win.model.geom.n_faces == 4
        assert _run(win, app) is not None
        assert not errors, errors[-1]
        win.close()


def test_csv_reader_details():
    """Semicolons with decimal commas, groups, and error messages."""
    with tempfile.TemporaryDirectory() as d:
        _write(os.path.join(d, "nodes.csv"),
               "id;x;y;z\n0;0;0;0\n1;2,5;0;0\n2;1,25;1;0,5\n")
        _write(os.path.join(d, "members.csv"),
               "0;1;7\n1;2;3\n0;2;3\n")
        m, notes = csv_import.import_truss_csv(d)
        assert m.geom.nodes[1, 0] == 2.5 and m.geom.nodes[2, 2] == 0.5
        assert m.geom.strut_group.tolist() == [1, 0, 0]
        assert m.support_type == geo.PINNED       # no supports file
        assert any("lowest" in n for n in notes)
        assert m.geom.title.startswith("Imported truss (")
        _write(os.path.join(d, "members.csv"), "0;1;7\n1;2\n")
        try:
            csv_import.import_truss_csv(d)
        except csv_import.CSVImportError as exc:
            assert "group on every member or on none" in str(exc)
        else:
            raise AssertionError("mixed groups accepted")
        _write(os.path.join(d, "members.csv"), "0;1\n")
        _write(os.path.join(d, "supports.csv"), "0;1;1;-5\n")
        try:
            csv_import.import_truss_csv(d)
        except csv_import.CSVImportError as exc:
            assert "supports.csv, line 1" in str(exc) and "negative" in str(exc)
        else:
            raise AssertionError("negative spring accepted")
        _write(os.path.join(d, "supports.csv"), "0;1;1;1\n1;0;1;1\n")
        _write(os.path.join(d, "loads.csv"), "2;0;0;-1;T\n")
        try:
            csv_import.import_truss_csv(d)
        except csv_import.CSVImportError as exc:
            assert "must be D, S, L or W" in str(exc)
        else:
            raise AssertionError("case T accepted")


def test_about_box():
    """Item 7: version, author, course, link and the licence note."""
    with _SlotErrors() as errors:
        app, win = _window()
        text = win.about_text()
        for part in (geodome.__version__, geodome.__author__,
                     geodome.__url__, "CE 7405", "LSU", "2026",
                     "All rights reserved. Free to use for personal, "
                     "educational and research calculations with credit to "
                     "the author. See LICENSE."):
            assert part in text, part
        assert geodome.__version__ == "2.0.0"
        with _Dialogs() as dlg:
            win.on_about()
        assert dlg.messages[0][0] == "about" and \
            geodome.__author__ in dlg.messages[0][2]
        assert not errors, errors[-1]
        win.close()


def test_design_wall_thickness_factor():
    """Addendum: the spin box reaches the model and the strength checks."""
    with _SlotErrors() as errors:
        app, win = _window()
        sp = win.design_panel.thickness_factor
        assert (sp.minimum(), sp.maximum(), sp.singleStep(), sp.value()) == \
            (0.50, 1.00, 0.01, 1.00)
        assert "AISC 360-16 B4.2" in sp.toolTip()
        _small_dome(win, app)
        r1 = _run(win, app)
        sp.setValue(0.93)
        r2 = _run(win, app)
        assert win.model.design_thickness_factor == 0.93
        assert r2.governing_util > r1.governing_util     # thinner wall
        win.load_preset("validation_a")
        assert sp.value() == 1.00                         # synced back
        assert not errors, errors[-1]
        win.close()


def test_model_changed_hook_and_extensions():
    """Addendum: model_changed() and extensions, errors logged not raised."""
    class Recorder:
        def __init__(self):
            self.models, self.results = [], []

        def on_model_changed(self, model):
            self.models.append(model)

        def on_results(self, results):
            self.results.append(results)

    class Broken:
        def on_model_changed(self, model):
            raise RuntimeError("broken extension")

        def on_results(self, results):
            raise RuntimeError("broken results")

    with _SlotErrors() as errors, tempfile.TemporaryDirectory() as d:
        app, win = _window()
        assert isinstance(win.extensions, list)
        rec = Recorder()
        win.extensions.extend([Broken(), rec, object()])
        win.model_changed()
        assert rec.models == [win.model]
        assert "ERROR in Broken.on_model_changed" in win.log_view.toPlainText()
        win.on_new()                                         # new
        _small_dome(win, app)                                # generate
        win.load_preset("validation_b")                      # preset
        path = os.path.join(d, "m.json")
        with _Dialogs(save_path=path, open_path=path):
            win.on_save()
            win.on_open()                                    # open
        ex = os.path.join(ROOT, "examples", "csv_truss", "nodes.csv")
        with _Dialogs(open_path=ex):
            win.on_import_csv()                              # CSV import
        assert len(rec.models) == 6
        assert rec.models[-1] is win.model
        res = _run(win, app)
        assert rec.results == [res]
        assert "ERROR in Broken.on_results" in win.log_view.toPlainText()
        assert not errors, errors[-1]
        win.close()


def test_public_menus_toolbar_and_view_menu():
    """add_menu() inserts before Help; the View menu lists added tabs."""
    from qtpy import QtWidgets
    with _SlotErrors() as errors:
        app, win = _window()
        titles = [a.text().replace("&", "") for a in win.menuBar().actions()]
        # add-on modules may add menus; ours keep their order
        mine = [t for t in titles
                if t in ("File", "Examples", "Analysis", "View", "Help")]
        assert mine == ["File", "Examples", "Analysis", "View", "Help"]
        assert titles[:2] == ["File", "Examples"] and titles[-1] == "Help"
        win.add_menu("&Extra")
        titles = [a.text().replace("&", "") for a in win.menuBar().actions()]
        assert titles[-2:] == ["Extra", "Help"]
        assert isinstance(win.toolbar, QtWidgets.QToolBar)
        assert win.act_size in win.toolbar.actions()
        win.view_tabs.addTab(QtWidgets.QWidget(), "Extra view")
        win.menu_view.aboutToShow.emit()
        names = [a.text() for a in win.menu_view.actions()]
        assert "Extra view" in names and "Model" in names
        win.menu_view.aboutToShow.emit()                    # no duplicates
        names = [a.text() for a in win.menu_view.actions()]
        assert names.count("Extra view") == 1
        acts = [a.text().replace("&", "") for a in win.menu_file.actions()]
        assert "Import truss from CSV..." in acts
        assert not errors, errors[-1]
        win.close()


def test_report_title_uses_the_model_title():
    """Item 3: the HTML report of a custom truss is titled by its name."""
    with _SlotErrors() as errors:
        app, win = _window()
        assert win._project_title().startswith("Geodesic Dome Analysis -- 4V")
        win.load_preset("three_bar")
        t = win._project_title()
        assert win.model.geom.title in t and "nan" not in t.lower()
        assert not errors, errors[-1]
        win.close()


# ======================================================= review fixes ==
def test_mechanism_supports_give_a_plain_failure():
    """Review GUI-02 / GUI-10: Elastic supports with every k = 0 are a
    mechanism.  F9 must fail with a plain message, not report OK."""
    with _SlotErrors() as errors, _Dialogs() as dlg:
        app, win = _window()
        _small_dome(win, app)
        gp = win.geometry_panel
        gp.set_support_type(geo.ELASTIC)
        for w in (gp.spring_kx, gp.spring_ky, gp.spring_kz):
            w.setValue(0.0)
        app.processEvents()
        assert _run(win, app) is None
        crit = dlg.texts("critical")
        assert crit and "mechanism" in crit[-1], dlg.messages
        assert "LinAlgError" not in crit[-1]
        assert win.status_label.text() == "Failed."
        assert "LinAlgError" in win.log_view.toPlainText()  # full text logged
        # kz only: still free to slide and turn
        gp.spring_kz.setValue(50000.0)
        app.processEvents()
        assert _run(win, app) is None
        assert "rigid-body" in dlg.texts("critical")[-1]
        assert not errors, errors[-1]
        win.close()


def test_plain_failure_messages():
    """Review GUI-10: class names are dropped, a singular matrix is
    explained."""
    from geodome.gui.worker import plain_message
    assert plain_message("ValueError: The radius is bad.") == \
        "The radius is bad."
    t = plain_message("FloatingPointError: singular tangent (Factor is "
                      "exactly singular)")
    assert "mechanism" in t and "Error" not in t
    t = plain_message("LinAlgError: Unstable structure: node 2 can move "
                      "freely in uy (zero stiffness at DOF(s) [7]).")
    assert t.startswith("Unstable structure: node 2")
    assert "missing" in plain_message("KeyError: 'geometry'")


def test_open_checks_the_geodesic_inputs():
    """Review GUI-09 / GUI-10: out-of-range or broken model files are
    refused with a plain sentence; the model on screen is kept."""
    with _SlotErrors() as errors, tempfile.TemporaryDirectory() as d:
        app, win = _window()
        good = default = win.model.to_dict()
        cases = []
        for key, val, word in (("radius", -5.0, "radius"),
                               ("radius", 0.0, "radius"),
                               ("radius", 500.0, "radius"),
                               ("frequency", 25, "frequency"),
                               ("frequency", 2.5, "frequency"),
                               ("dome_fraction", 3.0, "fraction"),
                               ("dome_fraction", 0.05, "fraction")):
            data = json.loads(json.dumps(default))
            data["geometry"][key] = val
            cases.append((data, word))
        cases.append(({}, '"geometry"'))
        cases.append(({"geometry": {"kind": "geodesic", "radius": 5.0}},
                      '"frequency"'))
        before = win.model
        for k, (data, word) in enumerate(cases):
            p = os.path.join(d, f"case{k}.json")
            with open(p, "w", encoding="utf-8") as fh:
                json.dump(data, fh)
            t0 = time.time()
            with _Dialogs(open_path=p) as dlg:
                win.on_open()
            assert time.time() - t0 < 10.0          # no 25V generation
            crit = dlg.texts("critical")
            assert crit and word in crit[-1], (data.get("geometry"), crit)
            assert "Error" not in crit[-1]
            assert win.model is before
        bad = os.path.join(d, "bad.json")
        with open(bad, "w", encoding="utf-8") as fh:
            fh.write("{not json")
        with _Dialogs(open_path=bad) as dlg:
            win.on_open()
        assert "not a valid model file" in dlg.texts("critical")[-1]
        # a good file still opens
        p = os.path.join(d, "good.json")
        with open(p, "w", encoding="utf-8") as fh:
            json.dump(good, fh)
        with _Dialogs(open_path=p) as dlg:
            win.on_open()
        assert not dlg.messages and win.model is not before
        # the core refuses a bad radius in plain words too
        for R in (0.0, -5.0, float("nan")):
            try:
                geo.generate_dome(4, R, 0.5)
            except ValueError as exc:
                assert "radius" in str(exc)
            else:
                raise AssertionError(f"radius {R} accepted")
        assert not errors, errors[-1]
        win.close()


def test_generate_from_an_example_uses_the_default_loads():
    """Review GUI-07: Generate dome after a custom example must not keep the
    example's loads (often only a nodal load that no longer exists)."""
    from geodome.model import LoadSpec, default_combinations
    with _SlotErrors() as errors, _Dialogs() as dlg:
        app, win = _window()
        win.load_preset("von_mises")
        assert win.model.combinations[0].name.startswith("P")
        gp = win.geometry_panel
        gp.frequency.setValue(4)
        gp.radius.setValue(12.0)
        gp.fraction.setValue(0.5)
        gp.generate_btn.click()
        app.processEvents()
        m = win.model
        assert not m.geom.is_custom
        assert [c.name for c in m.combinations] == \
            [c.name for c in default_combinations()]
        ref = LoadSpec()
        assert m.loads.include_self_weight and m.loads.wind_enabled
        assert m.loads.snow_pressure == ref.snow_pressure
        assert "replaced by the default dome loads" in \
            win.log_view.toPlainText()
        res = _run(win, app)
        assert res is not None and res.governing_util > 0.0
        assert "No load is applied" not in win.status_label.text()
        # a model with no load at all says so
        win.load_panel.self_weight.setChecked(False)
        for w in (win.load_panel.dead, win.load_panel.snow):
            w.setValue(0.0)
        win.load_panel.wind_box.setChecked(False)
        res = _run(win, app)
        assert "No load is applied" in win.status_label.text()
        assert not errors, errors[-1]
        win.close()


def test_all_service_combinations_are_refused():
    """Review GUI-08: with every row ticked Service there is no strength
    check; the Loads tab must say so instead of reporting OK."""
    from qtpy.QtCore import Qt
    with _SlotErrors() as errors, _Dialogs() as dlg:
        app, win = _window()
        lp = win.load_panel
        for i in range(lp.combos.rowCount()):
            lp.combos.item(i, 6).setCheckState(Qt.Checked)
        assert _run(win, app) is None
        w = dlg.texts("warning")
        assert w and "untick Service" in w[-1], dlg.messages
        sections = list(win.model.group_sections)
        win.on_auto_size()
        app.processEvents()
        assert win.model.group_sections == sections      # nothing sized
        assert not errors, errors[-1]
        win.close()


def test_file_actions_wait_for_a_busy_run():
    """Review GUI-13: New, Open, Save and Import are off during a run and
    say why when reached by a shortcut."""
    with _SlotErrors() as errors, _Dialogs() as dlg:
        app, win = _window()
        win._set_busy(True, "Working...")
        try:
            for a in (win.act_new, win.act_open, win.act_save,
                      win.act_import):
                assert not a.isEnabled(), a.text()
            for slot in (win.on_new, win.on_open, win.on_save,
                         win.on_import_csv):
                win.status("(marker)")
                slot()
                assert win.status_label.text().startswith("Busy"), slot
            assert not dlg.messages
        finally:
            win._set_busy(False)
        assert all(a.isEnabled() for a in win.file_actions)
        assert not errors, errors[-1]
        win.close()


def test_window_fits_a_small_screen():
    """Review GUI-05: the window must shrink to fit a 1366 x 768 screen
    (work area about 1366 x 690), also after a nonlinear run."""
    with _SlotErrors() as errors:
        app, win = _window()
        win.show()
        app.processEvents()
        win.resize(1280, 600)
        app.processEvents()
        assert win.height() <= 680, win.height()
        assert win.width() <= 1280, win.width()
        ext = win.nonlinear_extension
        win.load_preset("von_mises")
        from geodome.gui import nonlinear_view as nv
        s = nv.NonlinearDialog(win, win.model).values()
        assert ext.start_nonlinear(s, synchronous=True)
        app.processEvents()
        win.resize(1280, 600)
        app.processEvents()
        assert win.height() <= 680, win.height()
        assert not errors, errors[-1]
        win.close()


def test_wind_direction_and_live_load_tooltips():
    """Review F7 / V13 / F2: the inputs say which way the wind blows and
    where the live load acts."""
    app, win = _window()
    lp = win.load_panel
    assert "blows towards" in lp.wind_az.toolTip()
    assert "LC6" in lp.live.toolTip()
    assert "not the ASCE 7" in lp.snow_unbalanced.toolTip()
    # review V07: aluminium is flagged next to the material choice
    sp = win.section_panel
    assert "do not apply" not in sp.material_info.text()
    sp.material.setCurrentText("Aluminium 6061-T6")
    assert "do not apply" in sp.material_info.text()
    win.close()


def _wheel(widget, notches=-1):
    from qtpy import QtCore, QtGui
    pos = QtCore.QPointF(widget.width() / 2, widget.height() / 2)
    gpos = QtCore.QPointF(widget.mapToGlobal(pos.toPoint()))
    return QtGui.QWheelEvent(pos, gpos, QtCore.QPoint(0, 0),
                             QtCore.QPoint(0, 120 * notches),
                             QtCore.Qt.MouseButton.NoButton,
                             QtCore.Qt.KeyboardModifier.NoModifier,
                             QtCore.Qt.ScrollPhase.NoScrollPhase, False)


def test_mouse_wheel_does_not_edit_an_unfocused_box():
    """Review P1b: one wheel notch over the dead-load box changed 0.35 to
    0.30 kPa although the box had not been clicked (Qt's WheelFocus).
    Spin and combo boxes now take the focus by click or Tab only, and the
    wheel scrolls the panel unless the box has the focus."""
    from qtpy import QtCore, QtWidgets
    from geodome.gui import style as st
    app, win = _window()
    with _SlotErrors() as errors:
        win.show()
        app.processEvents()
        assert win.load_preset("case_study_6v")
        app.processEvents()
        boxes = (win.findChildren(QtWidgets.QAbstractSpinBox)
                 + win.findChildren(QtWidgets.QComboBox))
        assert len(boxes) > 20
        for b in boxes:
            b.ensurePolished()
        wheel_focus = [b.objectName() or type(b).__name__ for b in boxes
                       if b.focusPolicy() == QtCore.Qt.FocusPolicy.WheelFocus]
        assert not wheel_focus, wheel_focus
        lp = win.load_panel
        win.input_tabs.setCurrentWidget(
            win.input_tabs.widget(win.input_tabs.indexOf(lp))
            if win.input_tabs.indexOf(lp) >= 0 else win.input_tabs.widget(2))
        app.processEvents()
        box = lp.dead
        before = box.value()
        box.clearFocus()
        app.processEvents()
        assert not box.hasFocus()
        for notches in (-1, 1, -3):
            app.sendEvent(box, _wheel(box, notches))
            app.processEvents()
            assert box.value() == before, (notches, box.value(), before)
        # the same for a combo box
        combo = win.section_panel.material
        text = combo.currentText()
        combo.clearFocus()
        app.sendEvent(combo, _wheel(combo, -1))
        app.processEvents()
        assert combo.currentText() == text
        # a focused box still follows the wheel
        win.activateWindow()
        box.setFocus(QtCore.Qt.FocusReason.MouseFocusReason)
        app.processEvents()
        if box.hasFocus():              # needs an active window off screen
            app.sendEvent(box, _wheel(box, -1))
            app.processEvents()
            assert box.value() < before
            box.setValue(before)
        # one guard only, however many windows are built
        assert st.install_wheel_guard() is st.install_wheel_guard()
    assert not errors, errors[0]
    win.close()


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
