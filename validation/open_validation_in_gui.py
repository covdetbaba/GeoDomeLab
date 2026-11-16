"""
Open one of the hand-checked validation models in the GeoDomeLab GUI.

The model is loaded through the GUI's own File > Open code and analysed with
its own Run analysis button, so what you see is the program exactly as it is.
Nothing in the geodome package is changed by this script.

Usage (Spyder: set the arguments in Run > Configuration per file, or edit
DEFAULT_CASE below and press F5):

    python validation/open_validation_in_gui.py            # Model B, case B1
    python validation/open_validation_in_gui.py A2         # Model A, 1.2D+1.6S
    python validation/open_validation_in_gui.py B2         # Model B, 1.2D+1.6S
    python validation/open_validation_in_gui.py B1 --screenshot out.png

The model files are written by validation/run_validation.py.
"""

from __future__ import annotations

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)

DEFAULT_CASE = "B1"
FILES = {
    "A2": "model_A_LC2.json",
    "B1": "model_B_LC1.json",
    "B2": "model_B_LC2.json",
}
# Which result view to show first for each case.
FIRST_VIEW = {"A2": "Utilisation", "B1": "Axial force", "B2": "Deflected shape"}


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    shot = None
    if "--screenshot" in argv:
        i = argv.index("--screenshot")
        shot = argv[i + 1]
        del argv[i:i + 2]
    case = (argv[0] if argv else DEFAULT_CASE).upper()
    path = os.path.join(HERE, "output", FILES[case])
    if not os.path.exists(path):
        raise SystemExit(f"{path} not found - run validation/run_validation.py first")

    from qtpy import QtWidgets
    from geodome.gui.main_window import MainWindow
    from geodome.postprocess import run_analysis

    app = QtWidgets.QApplication.instance()
    created = app is None
    if created:
        app = QtWidgets.QApplication(sys.argv)
    app.setStyle("Fusion")

    win = MainWindow()
    win.resize(1600, 1000)

    # Use the GUI's own File > Open code: answer its file dialog with our path.
    real_dialog = QtWidgets.QFileDialog.getOpenFileName
    QtWidgets.QFileDialog.getOpenFileName = staticmethod(
        lambda *a, **k: (path, "JSON (*.json)"))
    try:
        win.on_open()
    finally:
        QtWidgets.QFileDialog.getOpenFileName = real_dialog
    win.setWindowTitle(win.windowTitle() + f"  --  validation case {case}")
    win.show()

    if shot is None:
        win.on_analyse()               # the normal Run analysis button (F9)
        _select_view(win, FIRST_VIEW[case], delay=True)
        if created:
            (getattr(app, "exec", None) or app.exec_)()
        return win

    # Screenshot mode: same display code, analysis run in this thread so the
    # picture is taken after the results are drawn.
    win._read_inputs()
    res = run_analysis(win.model,
                       n_buckling_modes=win.design_panel.n_modes.value(),
                       do_buckling=win.design_panel.do_buckling.isChecked())
    win._analysis_done(res)
    _select_view(win, FIRST_VIEW[case])
    for _ in range(8):
        app.processEvents()
    win.grab().save(shot)
    print(f"saved {shot}")
    return win


def _select_view(win, name, delay=False):
    def go():
        for i in range(win.view_tabs.count()):
            if win.view_tabs.tabText(i) == name:
                win.view_tabs.setCurrentIndex(i)
        win.result_tabs.setCurrentIndex(0)          # Summary
    if delay:
        # wait until the background analysis has drawn its views
        win._worker.succeeded.connect(lambda *_: go())
    else:
        go()


if __name__ == "__main__":
    main()
