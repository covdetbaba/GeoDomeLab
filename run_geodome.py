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
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)

    # ---------------------------------------------------------- verify ------
    if "--verify" in argv:
        from geodome import verification as vf
        rows = vf.run_all()
        print(vf.format_report(rows))
        return 0 if all(r["pass"] for r in rows) else 1

    # ------------------------------------------------------------ demo ------
    if "--demo" in argv:
        from examples.worked_example import main as demo_main
        return demo_main()

    # ----------------------------------------------------------- batch ------
    if "--batch" in argv:
        i = argv.index("--batch")
        rest = argv[i + 1:]
        freq = int(rest[0]) if len(rest) > 0 else 4
        radius = float(rest[1]) if len(rest) > 1 else 12.0
        frac = float(rest[2]) if len(rest) > 2 else 0.5
        out = rest[3] if len(rest) > 3 else "output"
        return batch(freq, radius, frac, out)

    # ------------------------------------------------------------- GUI ------
    from geodome.gui import launch
    launch()
    return 0


def batch(frequency, radius, fraction, out_dir):
    """Headless analyse-and-export, useful for scripted parameter studies."""
    import matplotlib
    matplotlib.use("Agg")
    from geodome.model import default_model
    from geodome.postprocess import run_analysis, summary_text
    from geodome import report as rp

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
