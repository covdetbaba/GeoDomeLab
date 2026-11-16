"""
GUI tests for the modal view (gui/modal_view.py) and the fabrication view
(gui/fabrication_view.py), driven through the real main window off screen.

* install() adds the menu entries, the toolbar button and the tabs;
* modal analysis of Validation Model A runs in the worker and reproduces the
  hand calculation F7 (validation/hand_calc_ext.py): crown frequencies with
  and without the dead load as mass;
* the mode selector, the frequency label and the result table follow the
  result; the animation draws frames and stops when the tab changes, when the
  model changes and when the window closes;
* the fabrication tables of a 2V dome match the published chord factors, the
  hub offset changes the cut lengths, NaN shows as "-" for a custom truss,
  a too large offset gives a message (no crash), and the CSV and DXF exports
  write their files;
* no exception escapes a Qt slot, and no dialog blocks the run.

Run with pytest::

    pytest -q tests/test_gui_modal_fab.py

or standalone::

    python tests/test_gui_modal_fab.py
"""

from __future__ import annotations

import atexit
import os
import shutil
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
from geodome import modal as md                                     # noqa: E402
from geodome.model import default_model                             # noqa: E402

_TMP = tempfile.mkdtemp(prefix="geodome_gui_w2_")
atexit.register(shutil.rmtree, _TMP, ignore_errors=True)


# ================================================================ helpers ===
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
    """Record message boxes instead of showing them (they would block)."""

    NAMES = ("warning", "information", "critical", "question")

    def __enter__(self):
        from qtpy import QtWidgets
        box = QtWidgets.QMessageBox
        self._box = box
        self._old = {n: getattr(box, n) for n in self.NAMES}

        def make(name):
            def fake(parent, title, text, *a, **k):
                self.append((name, title, text))
                return box.Yes if name == "question" else box.Ok
            return fake
        for n in self.NAMES:
            setattr(box, n, make(n))
        return self

    def __exit__(self, *exc):
        for n, f in self._old.items():
            setattr(self._box, n, f)
        return False


def _window():
    """Main window with the modal and fabrication extensions installed once."""
    from qtpy import QtWidgets
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    from geodome.gui.main_window import MainWindow
    from geodome.gui import modal_view, fabrication_view
    win = MainWindow()
    exts = list(getattr(win, "extensions", []) or [])
    modal = next((e for e in exts
                  if isinstance(e, modal_view.ModalExtension)), None)
    if modal is None:
        modal = modal_view.install(win)
    fabx = next((e for e in exts
                 if isinstance(e, fabrication_view.FabricationExtension)),
                None)
    if fabx is None:
        fabx = fabrication_view.install(win)
    return app, win, modal, fabx


def _set_model(win, model):
    """Replace the model the way File > Open / Examples do."""
    if isinstance(model, str):
        model = bm.make_preset(model)
    if hasattr(win, "set_model"):
        win.set_model(model)
        return model
    # older main window: refresh the panels and tell the extensions
    win.model = model
    win.results = None
    win.geometry_panel.show_model(model)
    win.section_panel.rebuild(model)
    win.load_panel.sync_from_model(model)
    win.design_panel.sync_from_model(model)
    for e in getattr(win, "extensions", []):
        e.on_model_changed(model)
    return model


def _wait(app, win, timeout=60.0):
    t0 = time.time()
    while win._busy and time.time() - t0 < timeout:
        app.processEvents()
        time.sleep(0.005)
    app.processEvents()
    assert not win._busy, "the worker did not finish"


def _pump(app, seconds):
    t0 = time.time()
    while time.time() - t0 < seconds:
        app.processEvents()
        time.sleep(0.005)


def _run_modal(app, win, modal, **options):
    assert modal.run_modal(options), "the modal run did not start"
    _wait(app, win)
    return modal.modal_result


def _hand_f7():
    sys.path.insert(0, os.path.join(ROOT, "validation"))
    import hand_calc_ext as hce
    return hce.results()


def _tab_titles(tabs):
    return [tabs.tabText(i) for i in range(tabs.count())]


def _plain(text):
    return text.replace("&", "")


def _close(win):
    """Close the window, first running any pending matplotlib draw_idle
    calls, so that none fires later on a canvas that has been deleted."""
    from qtpy import QtWidgets
    app = QtWidgets.QApplication.instance()
    with _NoDialogs():
        app.processEvents()
        win.close()
        app.processEvents()


# ================================================================== tests ===
def test_install_adds_menus_toolbar_and_tabs():
    from geodome.gui import ext_util as eu
    app, win, modal, fabx = _window()
    with _SlotErrors() as errors:
        analysis = eu.find_menu(win, "Analysis")
        assert "Modal analysis..." in [_plain(a.text())
                                       for a in analysis.actions()]
        file_menu = eu.find_menu(win, "File")
        texts = [_plain(a.text()) for a in file_menu.actions()]
        assert "Export geometry to CAD (DXF)..." in texts
        # the DXF entry sits above Exit, and Exit stays last
        assert texts[-1] == "Exit"
        assert texts.index("Export geometry to CAD (DXF)...") < len(texts) - 1
        tb = eu.main_toolbar(win)
        assert "Modes" in [a.text() for a in tb.actions()]
        assert _tab_titles(win.view_tabs).count("Modes") == 1
        assert _tab_titles(win.result_tabs).count("Modes") == 1
        assert _tab_titles(win.result_tabs).count("Fabrication") == 1
        # Log stays the last result tab
        assert _tab_titles(win.result_tabs)[-1] == "Log"
        assert sum(1 for e in win.extensions if e is modal) == 1
        assert sum(1 for e in win.extensions if e is fabx) == 1
        assert _tab_titles(fabx.subtabs)[:4] == ["Summary", "Struts", "Hubs",
                                                 "Panels"]
        # the default model already has its fabrication tables
        assert fabx.fabrication_tables is not None
        assert modal.modal_result is None
        # a second install returns the same extensions, no duplicate tabs
        from geodome.gui import modal_view, fabrication_view
        n_view, n_res = win.view_tabs.count(), win.result_tabs.count()
        assert modal_view.install(win) is modal
        assert fabrication_view.install(win) is fabx
        assert win.view_tabs.count() == n_view
        assert win.result_tabs.count() == n_res
        assert [_plain(a.text()) for a in analysis.actions()].count(
            "Modal analysis...") == 1
    _close(win)
    assert not errors, errors[-1]


def test_modal_validation_a_matches_hand_calc_f7():
    """Crown frequencies of Model A = hand calculation F7 (both mass cases)."""
    hand = _hand_f7()
    app, win, modal, fabx = _window()
    with _SlotErrors() as errors, _NoDialogs() as boxes:
        _set_model(win, "validation_a")
        for dead, tag in ((True, "dead"), (False, "sw")):
            res = _run_modal(app, win, modal, n_modes=12, include_dead=dead)
            assert res is not None, boxes
            f = np.sort(res.frequencies)
            fz = hand[f"F7.f_z_{tag}"]["value"]
            fxy = hand[f"F7.f_xy_{tag}"]["value"]
            assert res.n_modes == 3                  # only the crown moves
            assert abs(f[0] - fz) <= 1e-9 * fz, (f, fz)
            assert abs(f[1] - fxy) <= 1e-9 * fxy, (f, fxy)
            assert abs(f[2] - fxy) <= 1e-9 * fxy, (f, fxy)
            Tz = hand[f"F7.T_z_{tag}"]["value"]
            assert abs(res.periods[0] - Tz) <= 1e-9 * Tz
            # the GUI result is the same as a direct call of the core module
            direct = md.modal_analysis(bm.make_preset("validation_a"),
                                       n_modes=12, include_dead=dead)
            assert np.allclose(res.frequencies, direct.frequencies,
                               rtol=1e-12, atol=0)
            # table: 3 rows, f [Hz] column shows the hand value
            assert modal.table.rowCount() == 3
            from geodome.gui import ext_util as eu
            fcol = eu.cell_texts(modal.table, 1)
            assert fcol[0] == f"{fz:,.4f}", fcol
            assert eu.cell_texts(modal.table, 5)[0] == "100.00"   # Mz %
            # tiny ratios (1e-29) show as 0.00, not in e-notation
            assert eu.cell_texts(modal.table, 3) == ["0.00", "100.00",
                                                     "0.00"]
            assert eu.cell_texts(modal.table, 8)[-1] == "100.00"  # cum z
            assert eu.cell_texts(modal.table, 9) == ["z", "x", "y"]
            # selector and label
            assert modal.mode_combo.count() == 3
            assert f"{fz:.3f}" in modal.mode_combo.itemText(0)
            assert f"{fz:.4f}" in modal.freq_label.text()
            assert "T =" in modal.freq_label.text()
            assert win.view_tabs.currentWidget() is modal.view
            assert win.result_tabs.currentWidget() is modal.result_widget
        assert not boxes, boxes
    _close(win)
    assert not errors, errors[-1]


def test_mode_selector_redraws_and_keeps_options():
    app, win, modal, fabx = _window()
    with _SlotErrors() as errors, _NoDialogs():
        _set_model(win, "case_study_6v")
        res = _run_modal(app, win, modal, n_modes=8, include_dead=True,
                         snow_fraction=0.25, mass_type="consistent")
        assert res.n_modes == 8 and res.mass_type == "consistent"
        assert modal.modal_options["snow_fraction"] == 0.25
        direct = md.modal_analysis(win.model, n_modes=8, include_dead=True,
                                   snow_fraction=0.25,
                                   mass_type="consistent")
        assert np.allclose(res.frequencies, direct.frequencies, rtol=1e-9)
        for k in (1, 4, 7):
            modal.mode_combo.setCurrentIndex(k)
            app.processEvents()
            assert f"{res.frequencies[k]:.4f}" in modal.freq_label.text()
            assert modal.view.axes.get_title().startswith(f"Mode {k + 1}")
        # double-click on a table row selects that mode
        modal._on_table_double_click(2, 0)
        assert modal.mode_combo.currentIndex() == 2
        modal.amplitude.setValue(0.2)
        app.processEvents()
        # CSV export of the mode table
        path = modal.export_csv(os.path.join(_TMP, "modes.csv"))
        assert path and os.path.getsize(path) > 0
        with open(path, encoding="utf-8") as fh:
            assert len(fh.read().strip().splitlines()) == 1 + 8
    _close(win)
    assert not errors, errors[-1]


def test_animation_runs_and_stops_cleanly():
    app, win, modal, fabx = _window()
    with _SlotErrors() as errors, _NoDialogs():
        win.resize(1400, 900)
        win.show()
        _set_model(win, "validation_b")
        _run_modal(app, win, modal, n_modes=6)
        win.view_tabs.setCurrentWidget(modal.view)

        # 1. it animates: frames advance and the strut positions move
        assert modal.start_animation()
        assert modal.is_animating() and modal.animate_btn.isChecked()
        assert modal.timer.interval() == 40                  # 25 fps
        frames0 = modal.frames
        phases = set()
        # Wait for three frames; a slow (CI) machine gets up to 10 s, but
        # a fast one stops as soon as the frames have advanced.
        t0 = time.time()
        while time.time() - t0 < 0.6 or (
                modal.frames - frames0 < 3 and time.time() - t0 < 10.0):
            app.processEvents()
            phases.add(round(modal.phase, 6))
            time.sleep(0.005)
        assert modal.frames - frames0 >= 3, modal.frames - frames0
        assert len(phases) >= 3, phases               # the struts move
        # a frame is base + s * displacement with |s| <= 1, and it is what
        # the line collection holds
        for _ in range(5):
            modal._advance()
            seg = modal.segments
            assert abs(modal.phase) <= 1.0
            assert np.allclose(seg, modal._base + modal.phase * modal._disp)
            held = getattr(modal._lc, "_segments3d", None)
            if held is not None:                       # matplotlib internals
                assert np.allclose(np.asarray(held, dtype=float), seg)
        assert np.max(np.abs(modal._disp)) > 0.0
        # frames are blitted over the cached static scene
        assert modal._blit and modal._bg is not None
        assert modal._lc.get_animated()
        # changing the mode while animating keeps the animation running
        modal.mode_combo.setCurrentIndex(1)
        app.processEvents()
        assert modal.is_animating()

        # 2. switching the view tab stops it
        win.view_tabs.setCurrentIndex(0)
        app.processEvents()
        assert not modal.is_animating()
        assert not modal.animate_btn.isChecked()
        assert modal.animate_btn.text() == "Animate"
        assert not modal._lc.get_animated() and modal._bg is None
        # the static shape is back at full amplitude
        assert modal.phase == 1.0
        assert np.allclose(modal.segments, modal._base + modal._disp)

        # 3. a model change stops it and clears the result
        win.view_tabs.setCurrentWidget(modal.view)
        modal.animate_btn.click()                       # the real toggle
        assert modal.is_animating()
        _pump(app, 0.1)
        _set_model(win, "validation_a")
        app.processEvents()
        assert not modal.is_animating()
        assert modal.modal_result is None
        assert modal.table.rowCount() == 0 and modal.mode_combo.count() == 0
        # Animate without a result does nothing
        modal.animate_btn.click()
        assert not modal.is_animating() and not modal.animate_btn.isChecked()

        # 4. closing the window stops it
        _run_modal(app, win, modal, n_modes=3)
        win.view_tabs.setCurrentWidget(modal.view)
        assert modal.start_animation()
        _pump(app, 0.1)
        _close(win)
        app.processEvents()
        assert not modal.is_animating()
    assert not errors, errors[-1]


def test_animation_frame_time():
    """One animation frame of the 6V case study (555 struts) must stay
    quick enough to look smooth (< 0.25 s here, off screen, with a full
    redraw; the 40 ms timer simply skips frames on a slow machine)."""
    app, win, modal, fabx = _window()
    with _SlotErrors() as errors, _NoDialogs():
        _set_model(win, "case_study_6v")
        _run_modal(app, win, modal, n_modes=6)
        win.view_tabs.setCurrentWidget(modal.view)
        modal.view.canvas.draw()
        t0 = time.perf_counter()
        for _ in range(5):
            modal._advance()
            modal.view.canvas.draw()
        dt = (time.perf_counter() - t0) / 5
        # Shared CI runners draw off screen without a GPU and are several
        # times slower; there the limit only catches a pathological redraw.
        assert dt < (1.0 if os.environ.get("CI") else 0.25), dt
    _close(win)
    assert not errors, errors[-1]


def test_modal_result_discarded_when_model_changes_during_run():
    app, win, modal, fabx = _window()
    with _SlotErrors() as errors, _NoDialogs():
        _set_model(win, "validation_b")
        assert modal.run_modal({"n_modes": 4})
        modal.on_model_changed(win.model)     # as if the model was replaced
        _wait(app, win)
        assert modal.modal_result is None
        assert modal.act_modes.isEnabled() and modal.run_btn.isEnabled()
    _close(win)
    assert not errors, errors[-1]


def test_modal_custom_trusses_do_not_crash():
    app, win, modal, fabx = _window()
    with _SlotErrors() as errors, _NoDialogs() as boxes:
        for key in ("von_mises", "three_bar", "star_dome"):
            _set_model(win, key)
            res = _run_modal(app, win, modal, n_modes=12)
            assert res is not None, (key, boxes)
            assert modal.table.rowCount() == res.n_modes
            assert modal.start_animation()
            modal._advance()
            modal.stop_animation()
    _close(win)
    assert not errors, errors[-1]


def test_number_format():
    from geodome.gui import ext_util as eu
    assert eu.fmt(float("nan")) == "-"
    assert eu.fmt(None) == "-"
    assert eu.fmt(np.float64("nan")) == "-"
    assert eu.fmt(8.9e-29, 2) == "0.00"
    assert eu.fmt(-1e-9, 2) == "0.00"
    assert eu.fmt(2.54e-14, 4, small_sci=True) == "2.54e-14"
    assert eu.fmt(1234.5678, 1) == "1,234.6"
    assert eu.fmt(30) == "30"
    assert eu.fmt(np.int64(7)) == "7"
    assert eu.fmt(True) == "yes"


def test_options_dialog_values():
    from geodome.gui.modal_view import ModalOptionsDialog, DEFAULT_OPTIONS
    app, win, modal, fabx = _window()
    with _SlotErrors() as errors:
        dlg = ModalOptionsDialog(win)
        assert dlg.options() == DEFAULT_OPTIONS
        assert DEFAULT_OPTIONS["include_dead"] is True
        assert DEFAULT_OPTIONS["snow_fraction"] == 0.0
        dlg.n_modes.setValue(20)
        dlg.include_dead.setChecked(False)
        dlg.snow.setValue(0.5)
        dlg.mass_type.setCurrentIndex(dlg.mass_type.findData("consistent"))
        assert dlg.options() == {"n_modes": 20, "include_dead": False,
                                 "snow_fraction": 0.5,
                                 "mass_type": "consistent"}
        dlg2 = ModalOptionsDialog(win, dlg.options())
        assert dlg2.options() == dlg.options()
    _close(win)
    assert not errors, errors[-1]


def test_fabrication_tables_2v_dome():
    from geodome.gui import ext_util as eu
    app, win, modal, fabx = _window()
    with _SlotErrors() as errors, _NoDialogs() as boxes:
        _set_model(win, default_model(2, 5.0, 0.5))
        tabs = fabx.fabrication_tables
        assert tabs is not None
        t = fabx.tables
        assert t["struts"].rowCount() == 2
        assert t["hubs"].rowCount() == 3
        assert t["panels"].rowCount() == 2
        assert t["summary"].rowCount() == len(tabs["summary"])
        # published 2V chord factors (Domerama chart): 0.54653, 0.61803
        assert eu.cell_texts(t["struts"], 6) == ["0.54653", "0.61803"]
        assert eu.cell_texts(t["struts"], 0) == ["A", "B"]
        assert eu.cell_texts(t["struts"], 2) == ["30", "35"]
        assert eu.cell_texts(t["struts"], 1) == ["red", "blue"]
        # 2V tangent angles 15.86 and 18.00 deg
        assert eu.cell_texts(t["struts"], 8) == ["15.86", "18.00"]
        assert eu.cell_texts(t["panels"], 0) == ["AAB", "BBB"]
        # no hub offset: cut length = length
        L = eu.cell_texts(t["struts"], 4)
        assert eu.cell_texts(t["struts"], 5) == L
        # 50 mm hub offset: cut = L - 100 mm
        fabx.hub_offset.setValue(50.0)
        fabx.refresh()
        cuts = [r["cut_length_mm"] for r in fabx.fabrication_tables["struts"]]
        lens = [r["length_mm"] for r in fabx.fabrication_tables["struts"]]
        assert np.allclose(np.array(lens) - np.array(cuts), 100.0)
        assert eu.cell_texts(t["struts"], 5)[0] == f"{cuts[0]:,.1f}"
        assert "50.0" in eu.cell_texts(t["summary"], 1)
        # the debounced refresh fires by itself too
        fabx.hub_offset.setValue(20.0)
        _pump(app, 0.4)
        cuts = [r["cut_length_mm"] for r in fabx.fabrication_tables["struts"]]
        assert np.allclose(np.array(lens) - np.array(cuts), 40.0)
        assert "2V dome" in fabx.info.text()
        assert not boxes
    _close(win)
    assert not errors, errors[-1]


def test_fabrication_custom_truss_shows_dash_for_nan():
    from geodome.gui import ext_util as eu
    app, win, modal, fabx = _window()
    with _SlotErrors() as errors, _NoDialogs():
        _set_model(win, "three_bar")
        t = fabx.tables
        assert t["struts"].rowCount() == 3
        assert eu.cell_texts(t["struts"], 6) == ["-", "-", "-"]  # chord f.
        # no faces: one spanning row with a short message
        assert t["panels"].rowCount() == 1
        assert "no faces" in t["panels"].item(0, 0).text()
        assert "no panels" in fabx.info.text()
        for key in fabx.tables:
            for i in range(t[key].rowCount()):
                for j in range(t[key].columnCount()):
                    it = t[key].item(i, j)
                    assert it is None or it.text().lower() != "nan"
    _close(win)
    assert not errors, errors[-1]


def test_fabrication_offset_too_large_gives_message():
    app, win, modal, fabx = _window()
    with _SlotErrors() as errors, _NoDialogs() as boxes:
        _set_model(win, default_model(2, 5.0, 0.5))
        fabx.hub_offset.setValue(2000.0)          # beyond the 1000 mm max
        assert fabx.hub_offset.value() == 1000.0
        _set_model(win, default_model(2, 1.0, 0.5))   # struts ~0.55 m
        assert fabx.fabrication_tables is None
        assert "hub offset" in fabx.info.text().lower()
        assert fabx.tables["struts"].rowCount() == 0
        files = fabx.export_csv(folder=_TMP)
        assert files == [] and boxes and boxes[-1][0] == "warning"
        fabx.hub_offset.setValue(10.0)
        fabx.refresh()
        assert fabx.fabrication_tables is not None
        assert fabx.tables["struts"].rowCount() == 2
    _close(win)
    assert not errors, errors[-1]


def test_fabrication_follows_generate_and_sections():
    """Generate dome (through model_changed) and a section change refresh
    the tables."""
    app, win, modal, fabx = _window()
    with _SlotErrors() as errors, _NoDialogs():
        if hasattr(win, "model_changed"):
            gp = win.geometry_panel
            gp.frequency.setValue(3)
            gp.radius.setValue(6.0)
            gp.fraction.setValue(0.5)
            gp.generate_btn.click()
            app.processEvents()
            n = win.model.geom.n_strut_types()
            assert len(fabx.fabrication_tables["struts"]) == n == 3
        else:
            print("        (main window has no model_changed yet: generate "
                  "part skipped)")
        mass0 = fabx.fabrication_tables["struts"][0]["mass_each_kg"]
        combo = win.section_panel._combos[0]
        combo.setCurrentIndex((combo.currentIndex() + 5) % combo.count())
        _pump(app, 0.4)
        mass1 = fabx.fabrication_tables["struts"][0]["mass_each_kg"]
        assert mass1 != mass0
    _close(win)
    assert not errors, errors[-1]


def test_export_cut_list_csv():
    app, win, modal, fabx = _window()
    with _SlotErrors() as errors, _NoDialogs():
        _set_model(win, default_model(2, 5.0, 0.5))
        out = os.path.join(_TMP, "cutlist")
        files = fabx.export_csv(folder=out)
        names = sorted(os.path.basename(f) for f in files)
        assert "dome_2V_fab_struts.csv" in names
        assert "dome_2V_fab_hubs.csv" in names
        assert "dome_2V_fab_panels.csv" in names
        with open(os.path.join(out, "dome_2V_fab_struts.csv"),
                  encoding="utf-8") as fh:
            lines = fh.read().strip().splitlines()
        assert len(lines) == 3 and lines[0].startswith("family,")
        # a custom truss: NaN is an empty cell, prefix from the title
        _set_model(win, "three_bar")
        files = fabx.export_csv(folder=out)
        assert any("planar_three_bar" in os.path.basename(f) for f in files)
    _close(win)
    assert not errors, errors[-1]


def test_dxf_export_from_the_gui():
    from geodome.postprocess import run_analysis
    app, win, modal, fabx = _window()
    with _SlotErrors() as errors, _NoDialogs():
        _set_model(win, default_model(2, 5.0, 0.5))
        path = fabx.export_dxf(os.path.join(_TMP, "dome2v"))
        assert path.endswith(".dxf") and os.path.exists(path)
        with open(path, encoding="ascii") as fh:
            vals = [s.strip() for s in fh.read().splitlines()]
        pairs = list(zip(vals[0::2], vals[1::2]))
        assert pairs.count(("0", "LINE")) == 65
        assert pairs.count(("0", "3DFACE")) == 40
        assert pairs.count(("0", "POINT")) == 26
        ents = pairs[pairs.index(("2", "ENTITIES")):]
        assert not any(c == "62" for c, _ in ents)   # layer colours only
        # with results: utilisation colours on every LINE
        win.results = run_analysis(win.model, do_buckling=False)
        assert fabx._usable_results() is win.results
        path = fabx.export_dxf(os.path.join(_TMP, "dome2v_util.dxf"),
                               by_utilisation=True)
        with open(path, encoding="ascii") as fh:
            vals = [s.strip() for s in fh.read().splitlines()]
        pairs = list(zip(vals[0::2], vals[1::2]))
        ents = pairs[pairs.index(("2", "ENTITIES")):]
        n62 = sum(1 for c, _ in ents if c == "62")
        assert n62 == 65
        # results of another model are not used
        _set_model(win, default_model(3, 5.0, 0.5))
        assert fabx._usable_results() is None
        # the menu action exists and is wired (dialog not opened here)
        assert fabx.act_dxf.text().replace("&", "") == \
            "Export geometry to CAD (DXF)..."
    _close(win)
    assert not errors, errors[-1]


def test_window_grab_of_the_new_tabs():
    """The tabs render off screen (non-empty pixmaps)."""
    app, win, modal, fabx = _window()
    with _SlotErrors() as errors, _NoDialogs():
        win.resize(1400, 900)
        win.show()
        _set_model(win, "validation_b")
        _run_modal(app, win, modal, n_modes=6)
        app.processEvents()
        pm = win.grab()
        assert not pm.isNull() and pm.width() > 100
        win.result_tabs.setCurrentWidget(fabx.widget)
        app.processEvents()
        assert not fabx.widget.grab().isNull()
    _close(win)
    assert not errors, errors[-1]


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
