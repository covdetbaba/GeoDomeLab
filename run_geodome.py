#!/usr/bin/env python
"""
GeoDomeLab launcher -- run this file.

Spyder
------
Open this file and press F5.  If the window appears and immediately freezes,
set the IPython console's graphics backend to Qt:

    Tools -> Preferences -> IPython console -> Graphics -> Backend : Qt5
    (then restart the kernel)

"Automatic" also works.  With the Inline backend the script still runs, but
Spyder will block until the window is closed.

Command line
------------
    python run_geodome.py                 start the GUI
    python run_geodome.py --verify        run the verification suite only
    python run_geodome.py --demo          run the worked example, no GUI
    python run_geodome.py --batch 6 15 0.5 out/   analyse and export headless
    python run_geodome.py --version       print the version
    python run_geodome.py --help          list these options

The options are handled by ``geodome.__main__.main``, so this file,
``python -m geodome`` and the installed ``geodomelab`` command behave the
same way (an unknown option prints the usage and returns exit code 2).
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from geodome.__main__ import USAGE, batch  # noqa: E402,F401  (kept for scripts)


def main(argv=None):
    from geodome.__main__ import main as _main
    return _main(argv)


if __name__ == "__main__":
    raise SystemExit(main())
