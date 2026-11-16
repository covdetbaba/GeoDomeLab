"""
GUI tests for the nonlinear view (geodome/gui/nonlinear_view.py).

The main window runs off screen.  The tests install the nonlinear extension,
load the examples through the window, run the nonlinear path and the
imperfection study through the extension (in the worker thread and
synchronously) and check the numbers against closed forms and the core
module.

Run with pytest::

    pytest -q tests/test_gui_nonlinear.py

or standalone::

    python tests/test_gui_nonlinear.py
"""

from __future__ import annotations

import csv
import os
import sys
import tempfile
import traceback

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ.setdefault("QT_QPA_FONTDIR", "C:/Windows/Fonts")
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import matplotlib                                                   # noqa: E402
matplotlib.use("Agg")

import numpy as np                                                  # noqa: E402


def _window():
    from qtpy import QtWidgets
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    from geodome.gui.main_window import MainWindow
    from geodome.gui import nonlinear_view as nv
    win = MainWindow()
    win.resize(1600, 1000)
    ext = nv.install(win)
    ext.quiet = True                       # messages go to the log, no boxes
    return app, win, ext


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


def _load(win, ext, key):
    """Load an example the way the user does (Examples menu), or by hand
    on a window without that menu."""
    from geodome import benchmarks as bm
    if hasattr(win, "load_preset"):
        assert win.load_preset(key), key
        return win.model
    m = bm.make_preset(key)
    win.model = m
    win.results = None
    idx = win.geometry_panel.support.findData(m.support_type)
    if idx >= 0:
        win.geometry_panel.support.setCurrentIndex(idx)
    win.geometry_panel.show_model(m)
    win.section_panel.rebuild(m)
    win.load_panel.sync_from_model(m)
    win.design_panel.sync_from_model(m)
    ext.on_model_changed(m)
    return m


def _settings(win, **over):
    from geodome.gui import nonlinear_view as nv
    s = nv.NonlinearDialog(win, win.model).values()
    s.update(over)
    return s


def _labels(ax):
    return [ln.get_label() for ln in ax.get_lines()]


def _line(ax, label):
    for ln in ax.get_lines():
        if ln.get_label() == label:
            return ln
    raise AssertionError(f"no line {label!r} in {_labels(ax)}")


# ------------------------------------------------------------------ tests --
def test_install_adds_tab_actions_and_is_idempotent():
    from geodome.gui import nonlinear_view as nv
    app, win, ext = _window()
    with _SlotErrors() as errors:
        vt = [win.view_tabs.tabText(i) for i in range(win.view_tabs.count())]
        rt = [win.result_tabs.tabText(i)
              for i in range(win.result_tabs.count())]
        assert vt.count("Nonlinear") == 1, vt
        assert "Nonlinear" in rt and "Imperfection" in rt, rt
        if "Log" in rt:
            assert rt.index("Imperfection") < rt.index("Log")
        menus = [a.text().replace("&", "")
                 for m in win.menuBar().actions() if m.menu()
                 for a in m.menu().actions()]
        assert "Nonlinear analysis..." in menus
        assert "Imperfection study..." in menus
        tb = nv._main_toolbar(win)
        assert tb is not None and tb.parent() is win
        assert "Nonlinear" in [a.text() for a in tb.actions()]
        # a second install returns the same extension and adds nothing
        assert nv.install(win) is ext
        vt2 = [win.view_tabs.tabText(i) for i in range(win.view_tabs.count())]
        assert vt2 == vt
        if isinstance(getattr(win, "extensions", None), list):
            assert win.extensions.count(ext) == 1
        # empty state: placeholders, controls off, nothing stored
        win.view_tabs.setCurrentWidget(ext.view)
        app.processEvents()
        assert ext.nonlinear_result is None
        assert ext.imperfection_result is None
        assert not ext.view.slider.isEnabled()
        assert not errors, errors[-1]
        assert not ext.errors, ext.errors[-1]
    win.close()


def test_von_mises_limit_point_through_worker():
    """The worker run of the von Mises truss reaches the closed-form limit
    point and fills the plot."""
    from geodome import benchmarks as bm
    app, win, ext = _window()
    with _SlotErrors() as errors:
        _load(win, ext, "von_mises")
        s = _settings(win)
        assert s["control_node"] == int(win.model.geom.crown_node) == 2
        assert s["control_axis"] == "uz" and s["method"] == "arc-length"
        assert ext.start_nonlinear(s)
        assert ext.is_running() and win._busy
        assert ext.wait(120) == "done"
        assert not win._busy and not ext.is_running()
        r = ext.nonlinear_result
        a, h, EA, P = bm.von_mises_parameters(win.model)
        lp = bm.von_mises_limit_point(a, h, EA)
        assert r.limit_found and r.limit_type == "limit"
        assert abs(r.lambda_limit - lp["P"] / P) / (lp["P"] / P) < 1e-6
        assert abs(-r.control_limit - lp["w"]) / lp["w"] < 1e-6
        lam_cr = bm.von_mises_linear(a, h, EA, P)["lambda_cr"]
        assert abs(r.lambda_cr - lam_cr) / lam_cr < 1e-9
        # snap-through traced past w = 2h onto the stable branch
        assert float(np.max(-r.control)) > 2.0 * h
        # whole path on the closed form P(w)
        w = -r.control
        assert np.allclose(r.lam * P, bm.von_mises_load(w, a, h, EA),
                           rtol=1e-7, atol=1e-7 * P)
        # plot
        ax = ext.view.path_canvas.ax_path
        lab = _labels(ax)
        assert "Linear (first order)" in lab
        assert "P-Delta (linearised second order)" in lab
        cor = _line(ax, "Corotational (large displacement)")
        assert len(cor.get_xdata()) == r.n_points
        assert np.allclose(cor.get_ydata(), r.lam)
        assert np.allclose(cor.get_xdata(), -r.control * 1e3)  # down = +
        assert any(t.startswith(r"$\lambda_{cr}$") for t in lab)
        texts = " ".join(t.get_text() for t in ax.texts)
        assert "Limit point" in texts
        assert ext.view.slider.maximum() == r.n_points - 1
        assert ext.view.slider.value() == r.limit_step
        summary = ext.view.summary.text()
        assert "1.024" in summary and "5.328" in summary
        assert ext.path_table.rowCount() == r.n_points
        assert not errors, errors[-1]
        assert not ext.errors, ext.errors[-1]
    win.close()


def test_scrub_play_colour_and_exports():
    app, win, ext = _window()
    with _SlotErrors() as errors:
        _load(win, ext, "von_mises")
        assert ext.start_nonlinear(_settings(win), synchronous=True)
        r = ext.nonlinear_result
        v = ext.view
        for k in (0, r.n_points // 2, r.n_points - 1):
            ext.show_step(k)
            x, y = ext._marker.get_data()
            assert abs(y[0] - r.lam[k]) < 1e-12
            assert "Step %d:" % k in v.path_canvas.ax_3d.get_title()
        v.colour.setCurrentIndex(1)                      # axial force
        ext.show_step(3)
        v.undeformed.setChecked(False)
        ext.show_step(4)
        v.colour.setCurrentIndex(0)
        v.scale.setValue(5.0)
        v.auto_btn.click()
        v.true_btn.click()
        assert v.scale.value() == 1.0
        v.limit_btn.click()
        assert v.slider.value() == r.limit_step
        v.play_btn.setChecked(True)                      # one tick, then stop
        ext._play_tick()
        assert v.slider.value() == r.limit_step + 1
        v.play_btn.setChecked(False)
        app.processEvents()
        with tempfile.TemporaryDirectory() as tmp:
            png = os.path.join(tmp, "path.png")
            v.save(png, dpi=60)
            assert os.path.getsize(png) > 5000
            path_csv = ext.export_path_csv(os.path.join(tmp, "path.csv"))
            with open(path_csv, newline="", encoding="utf-8") as fh:
                rows = list(csv.reader(fh))
            assert len(rows) == r.n_points + 1
            assert abs(float(rows[-1][1]) - r.lam[-1]) < 1e-12
        assert not errors, errors[-1]
        assert not ext.errors, ext.errors[-1]
    win.close()


def test_case_study_6v_path_and_imperfection_study():
    """6V case study: synchronous path, then a two-amplitude imperfection
    study in the worker, then one imperfect run that must repeat the
    study's + sign result."""
    from geodome import nonlinear as nl
    from geodome.gui import nonlinear_view as nv
    app, win, ext = _window()
    with _SlotErrors() as errors:
        model = _load(win, ext, "case_study_6v")
        # the dialog default is the combination run_nonlinear picks itself
        assert nv.default_combination(model).name == \
            nl._find_combination(model, None).name == "LC2  1.2D + 1.6S"
        s = _settings(win)
        assert s["combination"] == "LC2  1.2D + 1.6S"
        assert s["control_node"] == int(model.geom.crown_node)
        assert ext.start_nonlinear(s, synchronous=True)
        r = ext.nonlinear_result
        assert r.limit_found and r.combination == "LC2  1.2D + 1.6S"
        assert 0.0 < r.limit_ratio < 1.0
        assert r.lambda_limit < r.lambda_cr
        assert "P-Delta (linearised second order)" in \
            _labels(ext.view.path_canvas.ax_path)

        # imperfection study, two amplitudes, through the worker
        dlg = nv.ImperfectionDialog(win, win.model)
        assert dlg.amplitudes.text() == nv.DEFAULT_AMPLITUDES_TEXT
        dlg.amplitudes.setText("span/1000, span/300")
        si = dlg.values()
        span = float(model.geom.base_diameter)
        assert [round(a, 12) for a, _ in si["amplitudes"]] == \
            [round(span / 1000, 12), round(span / 300, 12)]
        assert ext.start_imperfection(si)
        assert ext.wait(300) == "done"
        rows = ext.imperfection_result
        assert len(rows) == 3 and rows[0]["amplitude_m"] == 0.0
        assert abs(rows[1]["amplitude_over_span"] - 1e-3) < 1e-12
        assert abs(rows[2]["amplitude_over_span"] - 1.0 / 300) < 1e-12
        assert all(row["found"] for row in rows)
        assert abs(rows[0]["lambda_limit"] - r.lambda_limit) < \
            1e-6 * r.lambda_limit
        for row in rows[1:]:
            assert row["ratio"] < rows[0]["ratio"]
        assert ext.imperfection_labels == ["perfect", "span/1000", "span/300"]
        assert ext.imp_table.rowCount() == 3
        # review GUI-12: two amplitudes asked, two counted
        assert "2 amplitudes plus the perfect dome" in ext.last_status, \
            ext.last_status
        gov = _line(ext.view.imp_ax_ratio, "Governing (lower) sign")
        assert np.allclose(gov.get_ydata(), [row["ratio"] for row in rows])
        assert np.allclose(gov.get_xdata(),
                           [row["amplitude_over_span"] for row in rows])
        assert ext.view.tabs.currentIndex() == 1
        with tempfile.TemporaryDirectory() as tmp:
            p = ext.export_imperfection_csv(os.path.join(tmp, "imp.csv"))
            with open(p, newline="", encoding="utf-8") as fh:
                assert len(list(csv.reader(fh))) == 4

        # one imperfect run (+ mode 1, span/1000) = the study's + sign
        s1 = _settings(win, imperfection_mode=1,
                       imperfection_mm=span / 1000 * 1e3, stop_after_limit=5)
        assert ext.start_nonlinear(s1, synchronous=True)
        r1 = ext.nonlinear_result
        assert r1.limit_found
        assert abs(r1.lambda_limit - rows[1]["lambda_plus"]) < \
            1e-6 * rows[1]["lambda_plus"]
        assert abs(r1.lambda_cr - rows[0]["lambda_cr"]) < 1e-9 * r1.lambda_cr
        off = np.asarray(ext.nonlinear_offset)
        assert abs(np.max(np.abs(off)) - span / 1000) < 1e-12
        assert "imperfection" in ext.view.summary.text()
        assert not errors, errors[-1]
        assert not ext.errors, ext.errors[-1]
    win.close()


def test_failed_imperfection_sign_is_shown():
    """Review GUI-06: a sign whose run found no limit point reads 'failed'
    in the Imperfection table and is named in the summary."""
    from geodome.gui import nonlinear_view as nv
    nan = float("nan")
    assert nv._sign_cell(nan, 0.1) == "failed"
    assert nv._sign_cell(nan, 0.0) == "-"          # perfect: one sign only
    assert nv._sign_cell(16.30884, 0.1) == "16.309"
    app, win, ext = _window()
    with _SlotErrors() as errors:
        ext.imperfection_result = [
            {"amplitude_m": 0.0, "amplitude_over_span": 0.0,
             "amplitude_over_strut": 0.0, "lambda_limit": 34.2,
             "ratio": 0.605, "sign": 1.0, "lambda_plus": 34.2,
             "lambda_minus": nan, "lambda_critical": 34.2,
             "ratio_critical": 0.605, "lambda_cr": 56.5, "found": True},
            {"amplitude_m": 0.15, "amplitude_over_span": 0.005,
             "amplitude_over_strut": 0.07, "lambda_limit": 33.4,
             "ratio": 0.59, "sign": -1.0, "lambda_plus": nan,
             "lambda_minus": 33.4, "lambda_critical": 33.4,
             "ratio_critical": 0.59, "lambda_cr": 56.5, "found": True}]
        ext.imperfection_labels = ["perfect", "span/200"]
        ext.imperfection_settings = {"combination": "LC2  1.2D + 1.6S",
                                     "mode": 1, "method": "arc-length"}
        ext._fill_imp_table()
        cells = [ext.imp_table.item(1, j).text()
                 for j in range(ext.imp_table.columnCount())]
        assert "failed" in cells, cells
        assert "-" in [ext.imp_table.item(0, j).text()
                       for j in range(ext.imp_table.columnCount())]
        html = ext._imperfection_html()
        assert "one sign found no limit point for span/200" in html, html
        assert not errors, errors[-1]
    win.close()


def test_dialog_validation_and_amplitude_parsing():
    from geodome.gui import nonlinear_view as nv
    app, win, ext = _window()
    with _SlotErrors() as errors:
        _load(win, ext, "von_mises")
        dlg = nv.NonlinearDialog(win, win.model)
        dlg.node.setValue(0)                         # a pinned support
        try:
            dlg.values()
        except ValueError as exc:
            assert "support" in str(exc)
        else:
            raise AssertionError("a support DOF was accepted")
        dlg.node.setValue(2)
        dlg.axis.setCurrentText("uy")                # apex uy is restrained
        try:
            dlg.values()
        except ValueError:
            pass
        else:
            raise AssertionError("a restrained DOF was accepted")
        dlg.axis.setCurrentText("uz")
        dlg.imp_box.setChecked(True)
        dlg.imp_amp.setValue(0.0)
        try:
            dlg.values()
        except ValueError:
            pass
        else:
            raise AssertionError("a zero amplitude was accepted")
        dlg.imp_amp.setValue(-5.0)
        assert dlg.values()["imperfection_mm"] == -5.0
        dlg.target_auto.setChecked(False)
        dlg.target.setValue(2.5)
        dlg.inc_auto.setChecked(False)
        dlg.increment.setValue(0.05)
        v = dlg.values()
        assert v["lambda_max"] == 2.5 and v["first_increment"] == 0.05

        amps = nv.parse_amplitudes("span/1000, 25, L/500 ; /1000", 30.0)
        assert [round(a, 12) for a, _ in amps] == [0.025, 0.03, 0.06]
        assert [lab for _, lab in amps] == ["25 mm", "span/1000", "span/500"]
        for bad in ("", "span/abc", "-5", "span/0", "0"):
            try:
                nv.parse_amplitudes(bad, 30.0)
            except ValueError:
                pass
            else:
                raise AssertionError(f"{bad!r} was accepted")
        try:
            nv.parse_amplitudes("span/1000", float("nan"))
        except ValueError as exc:
            assert "mm" in str(exc)
        else:
            raise AssertionError("span/n accepted without a span")

        # a bad setting reaches the worker: clear message, no crash
        bad = _settings(win, imperfection_mode=7, imperfection_mm=10.0)
        assert ext.start_nonlinear(bad)
        assert ext.wait(60) == "failed"
        assert not win._busy and not ext.is_running()
        assert ext.last_message[0] == "critical"
        assert "buckling modes" in ext.last_message[2]
        assert not errors, errors[-1]
        assert not ext.errors, ext.errors[-1]
    win.close()


def test_menu_action_stop_and_model_change():
    """The menu action runs the dialog path; Stop ends a run cleanly; a new
    model clears the stored results."""
    app, win, ext = _window()
    with _SlotErrors() as errors:
        _load(win, ext, "von_mises")
        ext._exec_dialog = lambda dlg: True          # accept with defaults
        ext.act_nonlinear.trigger()
        assert ext.wait(120) == "done"
        assert ext.nonlinear_result is not None and \
            ext.nonlinear_result.limit_found

        _load(win, ext, "case_study_6v")
        if hasattr(win, "load_preset"):              # model_changed() hook
            assert ext.nonlinear_result is None
        ext.on_model_changed(win.model)
        assert ext.nonlinear_result is None
        assert not ext.view.slider.isEnabled()

        ext._exec_dialog = lambda dlg: True
        ext.act_imperfection.trigger()               # default 4 amplitudes
        assert ext.is_running()
        assert ext.stop_btn.isVisibleTo(win)
        ext.stop()
        assert ext.wait(120) == "cancelled"
        assert ext.imperfection_result is None
        assert not win._busy and not ext.is_running()
        assert not ext.stop_btn.isVisibleTo(win)
        assert ext.act_nonlinear.isEnabled()

        # Stop also ends a nonlinear path run (checked at every step)
        assert ext.start_nonlinear(_settings(win, stop_after_limit=200))
        ext.stop()
        assert ext.wait(120) == "cancelled"
        assert ext.nonlinear_result is None and not win._busy
        assert not errors, errors[-1]
        assert not ext.errors, ext.errors[-1]
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
            traceback.print_exc()
        else:
            npass += 1
            print(f"  pass  {name}")
    print(f"\n{npass} of {len(fns)} tests passed.")
    return 0 if npass == len(fns) else 1


if __name__ == "__main__":
    raise SystemExit(_run_standalone())
