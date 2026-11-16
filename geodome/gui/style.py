"""Application stylesheet and small widget helpers."""

from __future__ import annotations

from qtpy import QtWidgets
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
QComboBox::drop-down {{ border: none; width: 17px; }}
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

QToolBar {{
    background: #ffffff; border-bottom: 1px solid #c8d2db;
    spacing: 4px; padding: 4px 6px;
}}
QToolBar QToolButton {{ padding: 5px 11px; border-radius: 4px; }}
QToolBar QToolButton:hover {{ background: #e8eef3; }}
QToolBar QToolButton:disabled {{ color: #9aa7b3; }}

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
QSplitter::handle {{ background: #d4dce3; }}
QSplitter::handle:horizontal {{ width: 4px; }}
QSplitter::handle:vertical {{ height: 4px; }}

QCheckBox {{ spacing: 6px; }}
QLabel#heading {{ font-size: 13px; font-weight: 600; color: {ACCENT_DARK}; }}
QLabel#hint {{ color: #64748b; font-size: 11px; }}
QLabel#metric {{ font-family: "Consolas", monospace; font-size: 11px; }}
"""


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


def hline():
    f = QtWidgets.QFrame()
    f.setFrameShape(QtWidgets.QFrame.HLine)
    f.setFrameShadow(QtWidgets.QFrame.Sunken)
    return f


def spin(minimum, maximum, value, step=1, suffix=""):
    w = QtWidgets.QSpinBox()
    w.setRange(int(minimum), int(maximum))
    w.setSingleStep(int(step))
    w.setValue(int(value))
    if suffix:
        w.setSuffix(suffix)
    return w


def dspin(minimum, maximum, value, step=0.1, decimals=3, suffix=""):
    w = QtWidgets.QDoubleSpinBox()
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
                text = f"{v:,.{numeric_digits}f}"
            else:
                text = str(v)
            item = QtWidgets.QTableWidgetItem(text)
            if isinstance(v, (int, float)) and not isinstance(v, bool):
                item.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
            if flag:
                item.setBackground(bad)
            widget.setItem(i, j, item)
    widget.resizeColumnsToContents()
