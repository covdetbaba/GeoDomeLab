"""Embedded matplotlib canvases.

``pyplot`` is deliberately never imported: the figures are built directly from
``matplotlib.figure.Figure`` so that running inside Spyder (whose IPython
console installs its own inline backend) cannot steal or close them.
"""

from __future__ import annotations

import qtpy                                   # noqa: F401  (sets QT_API first)
from qtpy import QtWidgets
from qtpy.QtCore import Qt

from matplotlib.backends.backend_qtagg import (FigureCanvasQTAgg,
                                               NavigationToolbar2QT)
from matplotlib.figure import Figure

from ..plotting import axes_layout


def _qt_alive(obj):
    """False once the C++ side of a Qt object has been deleted."""
    try:
        obj.objectName()
        return True
    except RuntimeError:
        return False


class SafeCanvas(FigureCanvasQTAgg):
    """FigureCanvasQTAgg that ends a draw quietly when its widget is deleted
    during the draw (a closed window freed by the garbage collector while a
    late redraw of it runs)."""

    def draw(self):
        try:
            super().draw()
        except RuntimeError:
            if _qt_alive(self):
                raise


class MplCanvas(QtWidgets.QWidget):
    """A figure + navigation toolbar + optional extra control strip."""

    def __init__(self, parent=None, projection=None, toolbar=True):
        super().__init__(parent)
        self.figure = Figure(figsize=(7.0, 5.2), dpi=100)
        self.figure.patch.set_facecolor("#ffffff")
        self.canvas = SafeCanvas(self.figure)
        # small enough for the window to fit a 1366 x 768 screen
        self.canvas.setMinimumSize(380, 200)
        self._projection = projection
        self.axes = self.figure.add_subplot(111, projection=projection)
        self._colorbar = None

        lay = QtWidgets.QVBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.setSpacing(0)

        self.controls = QtWidgets.QWidget()
        self.controls_layout = QtWidgets.QHBoxLayout(self.controls)
        self.controls_layout.setContentsMargins(8, 4, 8, 4)
        self.controls_layout.setSpacing(8)
        self.controls.setVisible(False)
        lay.addWidget(self.controls)
        self._layout = lay

        # the drawing, with room on its right for a side panel (view3d)
        self.body = QtWidgets.QWidget()
        self.body_layout = QtWidgets.QHBoxLayout(self.body)
        self.body_layout.setContentsMargins(0, 0, 0, 0)
        self.body_layout.setSpacing(0)
        self.body_layout.addWidget(self.canvas, 1)
        lay.addWidget(self.body, 1)
        if toolbar:
            self.toolbar = NavigationToolbar2QT(self.canvas, self)
            self.toolbar.setIconSize(self.toolbar.iconSize() * 0.8)
            lay.addWidget(self.toolbar)
        else:
            self.toolbar = None
        # the interactive 3D viewer, when one is installed (gui/view3d.py)
        self.view3d = None

    # ------------------------------------------------------------- helpers --
    def add_control(self, widget, stretch=0):
        self.controls.setVisible(True)
        self.controls_layout.addWidget(widget, stretch)
        return widget

    def add_control_stretch(self):
        self.controls_layout.addStretch(1)

    def insert_bar(self, widget):
        """Put a row (the 3D view bar) between the controls and the plot."""
        self._layout.insertWidget(self._layout.indexOf(self.body), widget)
        return widget

    def add_side_panel(self, widget):
        """Put a panel on the right of the plot (the selection details)."""
        self.body_layout.addWidget(widget)
        return widget

    def clear(self):
        """Drop the old axes (and any colour bar) and start a fresh one."""
        viewer = getattr(self, "view3d", None)
        if viewer is not None:
            viewer.before_clear()             # keeps a toolbar zoom or pan
        self.figure.clear()
        self._colorbar = None
        self.axes = self.figure.add_subplot(111, projection=self._projection)
        if self.toolbar is not None:
            try:
                self.toolbar.update()     # the back/forward history is stale
            except Exception:                        # pragma: no cover
                pass
        return self.axes

    def add_colorbar(self, mappable, label=""):
        if mappable is None or not hasattr(mappable, "set_array"):
            return None
        self._colorbar = self.figure.colorbar(mappable, ax=self.axes,
                                              shrink=0.60, pad=0.10)
        self._colorbar.set_label(label, fontsize=9)
        self._colorbar.ax.tick_params(labelsize=8)
        return self._colorbar

    def refresh(self, tight=False):
        if tight:
            try:
                self.figure.tight_layout()
            except Exception:
                pass
        self.canvas.draw_idle()

    def save(self, path, dpi=200):
        self.figure.savefig(path, dpi=dpi, facecolor="white")


class Canvas3D(MplCanvas):
    """A 3D plot.  The plotting layout (plotting.AxesLayout) places the axes,
    a legend column and the title, so the figure margins are small."""

    def __init__(self, parent=None, toolbar=True):
        super().__init__(parent, projection="3d", toolbar=toolbar)
        self.figure.subplots_adjust(left=0.005, right=0.995, bottom=0.01,
                                    top=0.99)
        axes_layout(self.axes)

    def clear(self):
        ax = super().clear()
        axes_layout(ax)
        return ax


class Canvas2D(MplCanvas):
    def __init__(self, parent=None, toolbar=True):
        super().__init__(parent, projection=None, toolbar=toolbar)


class DualCanvas2D(QtWidgets.QWidget):
    """Two stacked 2D plots used for the chart tab."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.figure = Figure(figsize=(7.0, 7.4), dpi=100)
        self.figure.patch.set_facecolor("#ffffff")
        self.canvas = SafeCanvas(self.figure)
        self.ax_top = self.figure.add_subplot(211)
        self.ax_bottom = self.figure.add_subplot(212)
        lay = QtWidgets.QVBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.addWidget(self.canvas, 1)
        self.toolbar = NavigationToolbar2QT(self.canvas, self)
        self.toolbar.setIconSize(self.toolbar.iconSize() * 0.8)
        lay.addWidget(self.toolbar)
        # Lay the two plots out again when the pane is resized; otherwise a
        # shorter pane cut off the lower axis title.
        self.canvas.mpl_connect("resize_event", self._on_resize)

    def _on_resize(self, _event=None):
        try:
            self.figure.tight_layout(h_pad=2.4)
        except Exception:
            pass

    def refresh(self):
        try:
            self.figure.tight_layout(h_pad=2.4)
        except Exception:
            pass
        self.canvas.draw_idle()

    def save(self, path, dpi=200):
        self.figure.savefig(path, dpi=dpi, facecolor="white")
