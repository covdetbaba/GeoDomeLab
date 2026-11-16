"""Qt front end for GeoDomeLab.

The GUI is written against :mod:`qtpy`, so the identical source runs on PyQt5,
PyQt6 and PySide6 -- whichever binding the user's Spyder/Anaconda install
happens to ship.  PyQt6 is the binding of the CI runs.  The GUI tests also
pass with PyQt5 5.15 and PySide6 6.9 (local runs).  Under PySide6 a closed
main window is not freed until the process ends.
"""

from .main_window import MainWindow, launch

__all__ = ["MainWindow", "launch"]
