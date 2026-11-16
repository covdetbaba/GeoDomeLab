"""Qt front end for GeoDomeLab.

The GUI is written against :mod:`qtpy`, so the identical source runs on PyQt5,
PyQt6 and PySide6 -- whichever binding the user's Spyder/Anaconda install
happens to ship.
"""

from .main_window import MainWindow, launch

__all__ = ["MainWindow", "launch"]
