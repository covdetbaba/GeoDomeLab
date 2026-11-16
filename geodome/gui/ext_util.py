"""Small helpers shared by the GUI extension modules (modal, fabrication).

An extension module adds menu actions, toolbar buttons and tabs to an existing
:class:`MainWindow` without editing ``main_window.py``.  These helpers find the
right menu or toolbar position by the visible text of existing actions, so
they keep working when other extensions add their own entries.
"""

from __future__ import annotations

import math

from qtpy import QtWidgets
from qtpy.QtCore import Qt
from qtpy.QtGui import QColor


# ------------------------------------------------------------------ menus --
def _plain(text):
    return text.replace("&", "").strip().lower()


def find_menu(win, title):
    """The top-level menu whose title is ``title`` (ignoring '&'), or None."""
    for a in win.menuBar().actions():
        m = a.menu()
        if m is not None and _plain(a.text()) == _plain(title):
            return m
    return None


def ensure_menu(win, title):
    m = find_menu(win, title)
    return m if m is not None else win.menuBar().addMenu("&" + title)


def insert_in_first_group(menu, action):
    """Insert ``action`` at the end of the menu's first group (before the
    first separator), or append it when the menu has no separator."""
    for a in menu.actions():
        if a.isSeparator():
            menu.insertAction(a, action)
            return action
    menu.addAction(action)
    return action


def insert_before_exit(menu, action):
    """Insert ``action`` at the end of the group just above "Exit"."""
    acts = menu.actions()
    for i, a in enumerate(acts):
        if _plain(a.text()) in ("exit", "quit"):
            before = acts[i - 1] if i > 0 and acts[i - 1].isSeparator() else a
            menu.insertAction(before, action)
            return action
    menu.addAction(action)
    return action


def main_toolbar(win):
    tb = getattr(win, "toolbar", None)
    if isinstance(tb, QtWidgets.QToolBar):
        return tb
    # the window's own tool bars only (the matplotlib navigation bars inside
    # the canvases are QToolBars too)
    bars = [b for b in win.findChildren(QtWidgets.QToolBar)
            if b.parent() is win]
    return bars[0] if bars else None


def insert_toolbar_action(win, action):
    """Put ``action`` after the analysis buttons (Generate, Run analysis,
    Auto-size and any buttons other extensions added there), i.e. before the
    first separator that follows the Auto-size button.  Appends otherwise."""
    tb = main_toolbar(win)
    if tb is None:
        return None
    acts = tb.actions()
    anchor = getattr(win, "act_size", None)
    start = acts.index(anchor) + 1 if anchor in acts else None
    if start is not None:
        for a in acts[start:]:
            if a.isSeparator():
                tb.insertAction(a, action)
                return tb
    tb.addAction(action)
    return tb


def register(win, ext):
    """Add ``ext`` to ``win.extensions`` once (created if missing)."""
    exts = getattr(win, "extensions", None)
    if exts is None:
        exts = []
        win.extensions = exts
    if not any(e is ext for e in exts):
        exts.append(ext)
    return ext


def add_result_tab(win, widget, title):
    """Add a result tab before the Log tab (so Log stays last)."""
    tabs = win.result_tabs
    log = getattr(win, "log_view", None)
    idx = tabs.indexOf(log) if log is not None else -1
    if idx >= 0:
        tabs.insertTab(idx, widget, title)
    else:
        tabs.addTab(widget, title)
    return widget


# ------------------------------------------------------------- formatting --
def is_missing(v):
    if v is None:
        return True
    try:
        return isinstance(v, float) and not math.isfinite(v)
    except TypeError:
        return False


def fmt(v, decimals=3, small_sci=False):
    """Number as text; NaN / None as "-"; ints without decimals.

    Floats get ``decimals`` fixed decimals.  A value too small to show
    becomes 0 (never "-0.00"), unless ``small_sci`` is True: then it is
    written as 2.54e-14, which suits check values such as residuals.
    """
    if isinstance(v, bool):
        return "yes" if v else "no"
    if is_missing(v):
        return "-"
    if isinstance(v, int):
        return f"{v:,d}"
    if isinstance(v, float):
        if v != 0.0 and abs(v) < 10.0 ** (-decimals) * 0.5:
            if small_sci:
                return f"{v:.2e}"
            v = 0.0
        return f"{v:,.{decimals}f}"
    try:                                    # numpy scalars
        import numpy as np
        if isinstance(v, np.integer):
            return f"{int(v):,d}"
        if isinstance(v, np.floating):
            return fmt(float(v), decimals, small_sci)
    except Exception:                       # pragma: no cover
        pass
    return str(v)


def fill(table, columns, rows, empty_text=""):
    """Fill a QTableWidget.

    ``columns`` is a list of ``(key, header, decimals)``; ``decimals`` None
    means "text".  Numbers are right-aligned, NaN shows as "-".  With no rows
    and an ``empty_text``, one spanning row shows that text.
    """
    table.clearSpans()
    table.clearContents()
    table.setColumnCount(len(columns))
    table.setHorizontalHeaderLabels([c[1] for c in columns])
    if not rows:
        if empty_text:
            table.setRowCount(1)
            item = QtWidgets.QTableWidgetItem(empty_text)
            item.setForeground(QColor("#64748b"))
            table.setItem(0, 0, item)
            if len(columns) > 1:
                table.setSpan(0, 0, 1, len(columns))
        else:
            table.setRowCount(0)
        return
    table.setRowCount(len(rows))
    for i, row in enumerate(rows):
        for j, (key, _, dec) in enumerate(columns):
            v = row.get(key, "")
            numeric = (not isinstance(v, (str, bool)) and v is not None)
            text = fmt(v, 3 if dec is None else dec) if numeric else \
                ("-" if v is None else str(v))
            item = QtWidgets.QTableWidgetItem(text)
            if numeric:
                item.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
            table.setItem(i, j, item)
    table.resizeColumnsToContents()


def cell_texts(table, column):
    """All texts of one column (used by the tests)."""
    out = []
    for i in range(table.rowCount()):
        it = table.item(i, column)
        out.append(it.text() if it is not None else "")
    return out
