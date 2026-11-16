"""Application stylesheet and small widget helpers."""

from __future__ import annotations

import atexit

from qtpy import QtCore, QtWidgets
from qtpy.QtCore import Qt
from qtpy.QtGui import QColor

ACCENT = "#1f5f8b"
ACCENT_DARK = "#123a5c"
OK_GREEN = "#15803d"
BAD_RED = "#b4161b"
WARN = "#b45309"

STYLESHEET = f"""
QWidget {{
    font-family: "Segoe UI", "Calibri", sans-serif;
    font-size: 12px;
    color: #1c2733;
}}
QMainWindow, QDialog {{ background: #eef1f4; }}

QGroupBox {{
    border: 1px solid #c8d2db;
    border-radius: 5px;
    margin-top: 14px;
    background: #ffffff;
    padding-top: 6px;
}}
QGroupBox::title {{
    subcontrol-origin: margin;
    left: 10px;
    padding: 0 5px;
    color: {ACCENT_DARK};
    font-weight: 600;
}}

QTabWidget::pane {{ border: 1px solid #c8d2db; background: #ffffff; top: -1px; }}
QTabBar::tab {{
    background: #dde4ea;
    border: 1px solid #c8d2db;
    border-bottom: none;
    padding: 6px 13px;
    margin-right: 1px;
    border-top-left-radius: 4px;
    border-top-right-radius: 4px;
    color: #43566a;
}}
QTabBar::tab:selected {{ background: #ffffff; color: {ACCENT_DARK};
                         font-weight: 600; }}
QTabBar::tab:hover:!selected {{ background: #e8eef3; }}

QPushButton {{
    background: #ffffff;
    border: 1px solid #a9b8c6;
    border-radius: 4px;
    padding: 6px 13px;
}}
QPushButton:hover {{ background: #f0f5f9; border-color: {ACCENT}; }}
QPushButton:pressed {{ background: #dfe9f1; }}
QPushButton:disabled {{ color: #9aa7b3; background: #f2f4f6;
                        border-color: #d3dae1; }}
QPushButton#primary {{
    background: {ACCENT}; color: #ffffff; border: 1px solid {ACCENT_DARK};
    font-weight: 600; padding: 7px 17px;
}}
QPushButton#primary:hover {{ background: {ACCENT_DARK}; }}
QPushButton#primary:disabled {{ background: #94a9bb; border-color: #94a9bb; }}

QLineEdit, QSpinBox, QDoubleSpinBox, QComboBox {{
    background: #ffffff;
    border: 1px solid #b7c3ce;
    border-radius: 3px;
    padding: 3px 5px;
    min-height: 19px;
    selection-background-color: {ACCENT};
}}
QLineEdit:focus, QSpinBox:focus, QDoubleSpinBox:focus, QComboBox:focus {{
    border: 1px solid {ACCENT};
}}
QLineEdit:disabled, QSpinBox:disabled, QDoubleSpinBox:disabled,
QComboBox:disabled {{
    color: #9aa7b3; background: #f2f4f6; border-color: #d3dae1;
}}
QComboBox::drop-down {{ border: none; width: 17px; }}
QComboBox#cellcombo {{
    border: none; border-radius: 0; padding: 1px 4px; min-height: 0;
    background: transparent;
}}
QComboBox QAbstractItemView {{
    background: #ffffff; border: 1px solid {ACCENT};
    selection-background-color: {ACCENT}; selection-color: #ffffff;
}}

QHeaderView::section {{
    background: #2b4c6f; color: #ffffff; padding: 5px 7px;
    border: none; border-right: 1px solid #3d6287; font-weight: 600;
}}
QTableWidget {{
    background: #ffffff; gridline-color: #e3e9ee;
    alternate-background-color: #f6f9fb;
    selection-background-color: #cfe2f1; selection-color: #1c2733;
    border: 1px solid #c8d2db;
}}
QTableWidget::item {{ padding: 2px 4px; }}

QPlainTextEdit, QTextEdit {{
    background: #fbfcfd; border: 1px solid #c8d2db;
    font-family: "Consolas", "Courier New", monospace; font-size: 11px;
}}

QMenuBar {{
    background: #ffffff; color: #1c2733;
    border-bottom: 1px solid #c8d2db; padding: 1px 4px;
}}
QMenuBar::item {{
    background: transparent; color: #1c2733;
    padding: 4px 10px; border-radius: 3px;
}}
QMenuBar::item:selected {{ background: #e1ebf3; color: {ACCENT_DARK}; }}
QMenuBar::item:pressed {{ background: #cfe2f1; color: {ACCENT_DARK}; }}
QMenuBar::item:disabled {{ color: #9aa7b3; }}

QMenu {{
    background: #ffffff; color: #1c2733;
    border: 1px solid #b7c3ce; padding: 4px 0;
}}
QMenu::item {{
    background: transparent; color: #1c2733;
    padding: 5px 28px 5px 20px;
}}
QMenu::item:selected {{ background: {ACCENT}; color: #ffffff; }}
QMenu::item:disabled {{ color: #9aa7b3; background: transparent; }}
QMenu::separator {{ height: 1px; background: #dde4ea; margin: 4px 8px; }}

QToolTip {{
    background: #fffdf0; color: #1c2733;
    border: 1px solid #a9b8c6; padding: 4px 6px;
}}
QMessageBox {{ background: #eef1f4; }}
QMessageBox QLabel {{ color: #1c2733; background: transparent; }}

QToolBar {{
    background: #ffffff; border-bottom: 1px solid #c8d2db;
    spacing: 4px; padding: 4px 6px;
}}
QToolBar QToolButton {{
    padding: 5px 11px; border-radius: 4px; color: #1c2733;
    background: transparent;
}}
QToolBar QToolButton:hover {{ background: #e8eef3; }}
QToolBar QToolButton:pressed {{ background: #dfe9f1; }}
QToolBar QToolButton:disabled {{ color: #9aa7b3; }}
QScrollArea {{ background: transparent; }}
QScrollArea > QWidget > QWidget {{ background: transparent; }}

QStatusBar {{ background: #ffffff; border-top: 1px solid #c8d2db; }}
QStatusBar::item {{ border: none; }}

QProgressBar {{
    border: 1px solid #b7c3ce; border-radius: 3px; background: #eef1f4;
    text-align: center; height: 15px; font-size: 10px;
}}
QProgressBar::chunk {{ background: {ACCENT}; border-radius: 2px; }}

QSlider::groove:horizontal {{
    height: 4px; background: #ccd6df; border-radius: 2px;
}}
QSlider::handle:horizontal {{
    background: {ACCENT}; width: 13px; margin: -5px 0; border-radius: 7px;
}}
QSlider::handle:horizontal:disabled {{ background: #b7c3ce; }}
QSplitter::handle {{ background: #d4dce3; }}
QSplitter::handle:horizontal {{ width: 4px; }}
QSplitter::handle:vertical {{ height: 4px; }}

QCheckBox {{ spacing: 6px; }}
QLabel#heading {{ font-size: 13px; font-weight: 600; color: {ACCENT_DARK}; }}
QLabel#hint {{ color: #64748b; font-size: 11px; }}
QLabel#metric {{ font-family: "Consolas", monospace; font-size: 11px; }}
QLabel#note {{
    background: #fff7e6; border: 1px solid #f0d9a8; border-radius: 4px;
    padding: 6px 8px; color: #5c3b00;
}}
"""


def _arrow_rules():
    """Up/down arrows of spin boxes and combo boxes.

    Once the boxes above have a styled border, Qt no longer draws their
    arrows, so the fields looked like plain text boxes.  These rules draw
    small arrow images from ``gui/icons`` (Qt picks the @2x files on high-DPI
    screens).  Without the icon folder nothing changes.
    """
    import os
    folder = os.path.join(os.path.dirname(os.path.abspath(__file__)), "icons")
    if not os.path.isfile(os.path.join(folder, "arrow_down.png")):
        return ""
    d = folder.replace("\\", "/")
    return f"""
QAbstractSpinBox {{ padding-right: 18px; }}
QAbstractSpinBox::up-button, QAbstractSpinBox::down-button {{
    subcontrol-origin: border; width: 17px; border: none;
    background: transparent;
}}
QAbstractSpinBox::up-button {{ subcontrol-position: top right;
                               border-top-right-radius: 3px; }}
QAbstractSpinBox::down-button {{ subcontrol-position: bottom right;
                                 border-bottom-right-radius: 3px; }}
QAbstractSpinBox::up-button:hover, QAbstractSpinBox::down-button:hover {{
    background: #e1ebf3;
}}
QAbstractSpinBox::up-button:pressed, QAbstractSpinBox::down-button:pressed {{
    background: #cfe2f1;
}}
QAbstractSpinBox::up-arrow {{ image: url("{d}/arrow_up.png");
                             width: 9px; height: 5px; }}
QAbstractSpinBox::down-arrow {{ image: url("{d}/arrow_down.png");
                               width: 9px; height: 5px; }}
QAbstractSpinBox::up-arrow:disabled, QAbstractSpinBox::up-arrow:off {{
    image: url("{d}/arrow_up_disabled.png");
}}
QAbstractSpinBox::down-arrow:disabled, QAbstractSpinBox::down-arrow:off {{
    image: url("{d}/arrow_down_disabled.png");
}}
QComboBox::down-arrow {{ image: url("{d}/arrow_down.png");
                        width: 9px; height: 5px; }}
QComboBox::down-arrow:disabled {{
    image: url("{d}/arrow_down_disabled.png");
}}
QComboBox#cellcombo::down-arrow {{ width: 7px; height: 4px; }}
"""


STYLESHEET += _arrow_rules()


# ------------------------------------------------------------ light theme --
# Colours of the light palette.  The stylesheet above colours most widgets,
# but some parts (scroll bars, check-box indicators, spin-box arrows, the
# matplotlib toolbar icons) follow the palette.  Windows in dark mode hands
# Qt a dark palette, which gave dark text on a dark menu bar.  Forcing this
# palette makes the window look the same in light and dark OS themes.
_PALETTE = {
    "Window": "#eef1f4", "WindowText": "#1c2733",
    "Base": "#ffffff", "AlternateBase": "#f6f9fb",
    "ToolTipBase": "#fffdf0", "ToolTipText": "#1c2733",
    "PlaceholderText": "#7b8794",
    "Text": "#1c2733", "Button": "#f4f6f8", "ButtonText": "#1c2733",
    "BrightText": "#ffffff",
    "Light": "#ffffff", "Midlight": "#e3e9ee", "Mid": "#b7c3ce",
    "Dark": "#8a99a8", "Shadow": "#5b6b7b",
    "Highlight": ACCENT, "HighlightedText": "#ffffff",
    "Link": ACCENT, "LinkVisited": ACCENT_DARK,
}
_DISABLED = {"WindowText": "#9aa7b3", "Text": "#9aa7b3",
             "ButtonText": "#9aa7b3", "Highlight": "#94a9bb"}


def light_palette():
    """The light QPalette used by the application (independent of the OS)."""
    from qtpy.QtGui import QPalette
    pal = QPalette()
    role_enum = getattr(QPalette, "ColorRole", QPalette)
    group_enum = getattr(QPalette, "ColorGroup", QPalette)
    groups = [getattr(group_enum, g) for g in ("Active", "Inactive",
                                               "Disabled")]
    for name, colour in _PALETTE.items():
        role = getattr(role_enum, name, None)
        if role is None:
            continue
        for grp in groups:
            pal.setColor(grp, role, QColor(colour))
    for name, colour in _DISABLED.items():
        role = getattr(role_enum, name, None)
        if role is not None:
            pal.setColor(group_enum.Disabled, role, QColor(colour))
    return pal


def apply_light_theme(app=None):
    """Make the application light, whatever the OS colour scheme.

    Asks Qt (6.8 and later) for the light colour scheme and installs
    :func:`light_palette` on the application.  Safe to call more than once
    and a no-op when there is no QApplication.
    """
    app = app or QtWidgets.QApplication.instance()
    if app is None:
        return False
    try:
        hints = app.styleHints()
        scheme = getattr(getattr(Qt, "ColorScheme", None), "Light", None)
        if scheme is not None and hasattr(hints, "setColorScheme"):
            hints.setColorScheme(scheme)
    except Exception:                                # pragma: no cover
        pass
    app.setPalette(light_palette())
    return True


# ------------------------------------------------------------ wheel guard --
class WheelGuard(QtCore.QObject):
    """Keep the mouse wheel from editing an input box the user did not click.

    By default Qt gives a spin box or a combo box the focus when the wheel
    turns over it, and the box then changes its value.  Scrolling down a
    panel could therefore change a load or a section without the user
    noticing.  This application-wide filter gives those boxes
    ``StrongFocus`` (a click or Tab focuses them, the wheel does not) and,
    while a box does not have the focus, passes the wheel event on to the
    parent widget, so the panel scrolls instead.  A box that has the focus
    still responds to the wheel as usual.
    """

    _BOXES = (QtWidgets.QAbstractSpinBox, QtWidgets.QComboBox)
    _POLICY = (QtCore.QEvent.Type.Polish, QtCore.QEvent.Type.Show,
               QtCore.QEvent.Type.Enter)

    def eventFilter(self, obj, event):               # noqa: N802 (Qt name)
        try:
            t = event.type()
            if t == QtCore.QEvent.Type.Wheel:
                if isinstance(obj, self._BOXES) and not obj.hasFocus():
                    event.ignore()                   # the parent scrolls
                    return True
            elif t in self._POLICY and isinstance(obj, self._BOXES):
                if obj.focusPolicy() == Qt.FocusPolicy.WheelFocus:
                    obj.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
        except Exception:                            # pragma: no cover
            pass
        return False


_wheel_guard = None


def remove_wheel_guard():
    """Take the guard off the application.

    Called when the application quits and when Python exits, so that Qt
    never calls the Python filter while the interpreter shuts down.
    """
    global _wheel_guard
    guard, _wheel_guard = _wheel_guard, None
    if guard is None:
        return
    try:
        app = QtWidgets.QApplication.instance()
        if app is not None:
            app.removeEventFilter(guard)
        guard.setParent(None)
    except Exception:                                # pragma: no cover
        pass


def install_wheel_guard(app=None):
    """Install :class:`WheelGuard` on the application (once).

    Returns the guard, or None when there is no QApplication.
    """
    global _wheel_guard
    app = app or QtWidgets.QApplication.instance()
    if app is None:
        return None
    if _wheel_guard is None or _wheel_guard.parent() is not app:
        remove_wheel_guard()
        _wheel_guard = WheelGuard(app)
        app.installEventFilter(_wheel_guard)
        app.aboutToQuit.connect(remove_wheel_guard)
    return _wheel_guard


atexit.register(remove_wheel_guard)


# ----------------------------------------------------------------- helpers --
def heading(text):
    lab = QtWidgets.QLabel(text)
    lab.setObjectName("heading")
    return lab


def hint(text):
    lab = QtWidgets.QLabel(text)
    lab.setObjectName("hint")
    lab.setWordWrap(True)
    return lab


def note(text):
    """A boxed, word-wrapped note (used for the custom-model notice)."""
    lab = QtWidgets.QLabel(text)
    lab.setObjectName("note")
    lab.setWordWrap(True)
    return lab


def input_table(columns, min_height=84, max_height=132):
    """An editable table for row input (settlements, nodal loads)."""
    t = QtWidgets.QTableWidget(0, len(columns))
    t.setHorizontalHeaderLabels(columns)
    t.verticalHeader().setVisible(False)
    t.verticalHeader().setDefaultSectionSize(26)
    t.setAlternatingRowColors(True)
    t.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
    t.setMinimumHeight(min_height)
    t.setMaximumHeight(max_height)
    hh = t.horizontalHeader()
    hh.setSectionResizeMode(QtWidgets.QHeaderView.Stretch)
    return t


def fmt_float(v, digits=3):
    """Format a number for a table cell; NaN and inf read 'n/a'."""
    try:
        f = float(v)
    except (TypeError, ValueError):
        return str(v)
    if f != f or f in (float("inf"), float("-inf")):
        return "n/a"
    return f"{f:,.{digits}f}"


def hline():
    f = QtWidgets.QFrame()
    f.setFrameShape(QtWidgets.QFrame.HLine)
    f.setFrameShadow(QtWidgets.QFrame.Sunken)
    return f


def spin(minimum, maximum, value, step=1, suffix=""):
    w = QtWidgets.QSpinBox()
    w.setFocusPolicy(Qt.FocusPolicy.StrongFocus)     # see WheelGuard
    w.setRange(int(minimum), int(maximum))
    w.setSingleStep(int(step))
    w.setValue(int(value))
    if suffix:
        w.setSuffix(suffix)
    return w


def dspin(minimum, maximum, value, step=0.1, decimals=3, suffix=""):
    w = QtWidgets.QDoubleSpinBox()
    w.setFocusPolicy(Qt.FocusPolicy.StrongFocus)     # see WheelGuard
    w.setDecimals(decimals)
    w.setRange(float(minimum), float(maximum))
    w.setSingleStep(float(step))
    w.setValue(float(value))
    if suffix:
        w.setSuffix(suffix)
    return w


def table(columns, stretch_last=True, rows=0):
    t = QtWidgets.QTableWidget(rows, len(columns))
    t.setHorizontalHeaderLabels(columns)
    t.setAlternatingRowColors(True)
    t.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
    t.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
    t.verticalHeader().setVisible(False)
    t.verticalHeader().setDefaultSectionSize(21)
    hh = t.horizontalHeader()
    hh.setStretchLastSection(stretch_last)
    hh.setSectionResizeMode(QtWidgets.QHeaderView.ResizeToContents)
    return t


def fill_table(widget, rows, highlight=None, numeric_digits=3):
    """Populate a QTableWidget from a list of dicts."""
    widget.clearContents()
    if not rows:
        widget.setRowCount(0)
        return
    cols = list(rows[0].keys())
    widget.setColumnCount(len(cols))
    widget.setHorizontalHeaderLabels(cols)
    widget.setRowCount(len(rows))
    bad = QColor("#ffe0e0")
    for i, row in enumerate(rows):
        flag = bool(highlight(row)) if highlight else False
        for j, c in enumerate(cols):
            v = row.get(c, "")
            if isinstance(v, float):
                text = fmt_float(v, numeric_digits)
            else:
                text = str(v)
            item = QtWidgets.QTableWidgetItem(text)
            if isinstance(v, (int, float)) and not isinstance(v, bool):
                item.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
            if flag:
                item.setBackground(bad)
            widget.setItem(i, j, item)
    widget.resizeColumnsToContents()
