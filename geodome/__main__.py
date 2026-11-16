"""
Command line entry point: ``python -m geodome`` and the ``geodomelab``
command (installed by ``pip install -e .``).

It behaves like ``run_geodome.py`` and takes the same options:

    python -m geodome                    start the GUI
    python -m geodome --verify           run the verification suite only
    python -m geodome --demo             run the worked example, no GUI
    python -m geodome --batch F R FRACTION OUT
                                         analyse a dome and export the
                                         results, no GUI
    python -m geodome --version          print the version
    python -m geodome --help             show this text

``--demo`` runs ``examples/worked_example.py`` of the repository, so it
needs the cloned repository (an editable install keeps it in place).
"""

from __future__ import annotations

import importlib.util
import os
import sys

from . import __version__

USAGE = f"""GeoDomeLab {__version__} -- geodesic dome space-truss analyzer

usage:
    python -m geodome                 start the GUI
    python -m geodome --verify        run the verification suite only
    python -m geodome --demo          run the worked example, no GUI
    python -m geodome --batch F R FRACTION OUT
                                      analyse a dome and export the
                                      results, no GUI (defaults 4 12 0.5
                                      output)
    python -m geodome --version       print the version
    python -m geodome --help          show this text

The command geodomelab (after pip install -e .) takes the same options.
"""

# The repository root: the folder that holds geodome/, examples/ and
# validation/ in a clone (also after an editable install).
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def main(argv=None):
    """Run the command line.  Returns the exit code."""
    argv = list(sys.argv[1:] if argv is None else argv)

    if "--help" in argv or "-h" in argv:
        print(USAGE)
        return 0

    if "--version" in argv:
        print(f"GeoDomeLab {__version__}")
        return 0

    if "--verify" in argv:
        from . import verification as vf
        rows = vf.run_all()
        print(vf.format_report(rows))
        return 0 if all(r["pass"] for r in rows) else 1

    if "--demo" in argv:
        return demo()

    if "--batch" in argv:
        i = argv.index("--batch")
        rest = argv[i + 1:]
        try:
            freq = int(rest[0]) if len(rest) > 0 else 4
            radius = float(rest[1]) if len(rest) > 1 else 12.0
            frac = float(rest[2]) if len(rest) > 2 else 0.5
        except ValueError:
            print("--batch needs: frequency (whole number), radius [m], "
                  "dome fraction, output folder.  Example: --batch 6 15 0.5 "
                  "out")
            return 2
        out = rest[3] if len(rest) > 3 else "output"
        return batch(freq, radius, frac, out)

    unknown = [a for a in argv if a.startswith("-")]
    if unknown:
        print(f"Unknown option: {' '.join(unknown)}\n")
        print(USAGE)
        return 2

    from .gui import launch
    launch()
    return 0


def demo():
    """Run ``examples/worked_example.py`` from the repository."""
    path = os.path.join(REPO_ROOT, "examples", "worked_example.py")
    if not os.path.isfile(path):
        print("The worked example is in the examples folder of the "
              "repository.  Clone https://github.com/covdetbaba/GeoDomeLab "
              "and run: python run_geodome.py --demo")
        return 1
    if REPO_ROOT not in sys.path:
        sys.path.insert(0, REPO_ROOT)
    spec = importlib.util.spec_from_file_location("geodome_worked_example",
                                                  path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.main()


def batch(frequency, radius, fraction, out_dir):
    """Headless analyse-and-export, the same as ``run_geodome.py --batch``."""
    import matplotlib
    matplotlib.use("Agg")
    from .model import default_model
    from .postprocess import run_analysis, summary_text
    from . import report as rp

    os.makedirs(out_dir, exist_ok=True)
    model = default_model(frequency, radius, fraction)
    res = run_analysis(model)
    print(summary_text(res))
    files = rp.export_csv(res, out_dir)
    files += rp.export_figures(res, out_dir)
    html = os.path.join(out_dir, "dome_report.html")
    rp.export_html(res, html)
    files.append(html)
    print(f"\n{len(files)} files written to {os.path.abspath(out_dir)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
