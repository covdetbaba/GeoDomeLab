"""
Worked example -- the case study used in the written report and the live demo.

A 6V Class I geodesic dome, 15 m sphere radius (30 m span, 15 m rise), carrying
dead, snow and wind load, is analysed, auto-sized, re-analysed and exported.
Run it headless with::

    python examples/worked_example.py
    python run_geodome.py --demo

Everything it writes lands in ``output/worked_example/``.
"""

from __future__ import annotations

import os
import sys

import matplotlib
matplotlib.use("Agg")                       # headless: no window, no pyplot

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np                                               # noqa: E402

from geodome import geometry as geo                              # noqa: E402
from geodome import report as rp                                 # noqa: E402
from geodome.model import DomeModel, LoadSpec, LoadCombination   # noqa: E402
from geodome.optimize import auto_size                           # noqa: E402
from geodome.postprocess import run_analysis, summary_text       # noqa: E402

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                   "output", "worked_example")


def build_model():
    """The case-study dome: 30 m span exhibition hall, Baton Rouge loading."""
    geometry = geo.generate_dome(frequency=6, radius=15.0, dome_fraction=0.50)

    model = DomeModel(
        geom=geometry,
        material_name="ASTM A500 Gr. C (round)",
        support_type=geo.PINNED,
    )
    model.group_sections = ["CHS 114.3 x 5.0"] * geometry.n_strut_types()

    model.loads = LoadSpec(
        include_self_weight=True,
        dead_pressure=350.0,        # 0.35 kPa  insulated panels + services
        snow_pressure=600.0,        # 0.60 kPa  ground snow, projected
        snow_unbalanced=False,
        live_pressure=0.0,
        wind_enabled=True,
        wind_speed=58.0,            # 58 m/s ~ 130 mph 3 s gust (Gulf coast)
        wind_Kz=0.95, wind_Kzt=1.0, wind_Kd=0.95, wind_G=0.85,
        crown_point_load=0.0,
    )
    model.combinations = [
        LoadCombination("LC1  1.4D", D=1.4),
        LoadCombination("LC2  1.2D + 1.6S", D=1.2, S=1.6),
        LoadCombination("LC3  1.2D + 1.0W + 0.5S", D=1.2, W=1.0, S=0.5),
        LoadCombination("LC4  0.9D + 1.0W", D=0.9, W=1.0),
        LoadCombination("SLS  1.0D + 1.0S", D=1.0, S=1.0, serviceability=True),
    ]
    model.K_effective = 1.0
    model.deflection_limit = 240.0
    return model


def rule(title):
    print("\n" + "=" * 78)
    print("  " + title)
    print("=" * 78)


def main():
    os.makedirs(OUT, exist_ok=True)
    model = build_model()

    rule("STEP 1  geometry")
    for line in model.summary_lines():
        print("  " + line)
    print()
    print(f"  {'family':<8s}{'length [m]':>12s}{'chord factor':>15s}{'count':>8s}")
    for s in model.geom.strut_schedule():
        print(f"  {s['label']:<8s}{s['length']:>12.4f}"
              f"{s['chord_factor']:>15.6f}{s['count']:>8d}")

    rule("STEP 2  first-pass analysis with a trial section")
    res = run_analysis(model, n_buckling_modes=6)
    print(summary_text(res))

    rule("STEP 3  fully stressed auto-sizing")
    sizing = auto_size(model, target_util=0.95, max_iter=12,
                       enforce_deflection=True, slenderness_limit=200.0)
    print(f"  {'iter':<6s}{'mass [kg]':>12s}{'max util':>11s}"
          f"{'deflection':>14s}   sections")
    for row in sizing.table():
        print(f"  {row['iteration']:<6d}{row['mass_kg']:>12,.0f}"
              f"{row['max_utilisation']:>11.3f}{row['deflection']:>14s}"
              f"   {row['sections']}")
    direction = "lighter" if sizing.saving >= 0 else "heavier"
    print(f"\n  {sizing.start_mass:,.0f} kg -> {sizing.final_mass:,.0f} kg "
          f"({abs(sizing.saving) * 100:.1f} % {direction}), "
          f"{'converged' if sizing.converged else 'not converged'}")
    for m in sizing.messages:
        print("  * " + m)

    rule("STEP 4  final analysis of the sized dome")
    res = run_analysis(model, n_buckling_modes=6)
    print(summary_text(res))

    rule("STEP 5  export")
    files = rp.export_csv(res, OUT, prefix="case_study")
    files += rp.export_figures(res, OUT, prefix="case_study")
    html = os.path.join(OUT, "case_study_report.html")
    from geodome import verification as vf
    checks = vf.run_all()
    n_pass = sum(1 for c in checks if c["pass"])
    print(f"  verification suite: {n_pass} of {len(checks)} checks passed")
    rp.export_html(res, html, sizing=sizing, verification=checks,
                   project_title="Case study -- 30 m span 6V geodesic dome",
                   author="GeoDomeLab worked example")
    files.append(html)
    model.save(os.path.join(OUT, "case_study_model.json"))
    files.append(os.path.join(OUT, "case_study_model.json"))
    for f in files:
        print("  " + os.path.relpath(f))
    print(f"\n  {len(files)} files written to {OUT}")

    rule("STEP 6  parameter study -- dome fraction vs steel mass and stability")
    print(f"  {'fraction':>9s}{'rise [m]':>10s}{'span [m]':>10s}"
          f"{'mass [kg]':>12s}{'max util':>11s}{'lambda_cr':>12s}{'defl':>13s}")
    study = []
    for frac in (0.30, 0.40, 0.50, 0.60, 0.70):
        m2 = build_model()
        m2.geom = geo.generate_dome(6, 15.0, frac)
        m2.group_sections = ["CHS 114.3 x 5.0"] * m2.geom.n_strut_types()
        r2 = run_analysis(m2, n_buckling_modes=3)
        lam = (r2.buckling.load_factors[0]
               if r2.buckling is not None and len(r2.buckling.load_factors)
               else float("nan"))
        # Where does the critical mode actually live?  Weight the nodes by the
        # square of the mode amplitude and report the height of the peak.
        where = "-"
        if r2.buckling is not None and len(r2.buckling.load_factors):
            mag = np.linalg.norm(r2.buckling.mode(0), axis=1)
            z = m2.geom.nodes[:, 2]
            where = f"{z[int(np.argmax(mag))] / m2.geom.height:.2f}"
        study.append({"frac": m2.geom.dome_fraction, "rise": m2.geom.height,
                      "span": m2.geom.base_diameter, "mass": m2.steel_mass(),
                      "util": r2.governing_util, "lam": lam,
                      "defl": r2.deflection_ratio, "mode_z": where,
                      "ref": r2.buckling.reference_case
                      if r2.buckling is not None else "-"})
        print(f"  {m2.geom.dome_fraction:>9.3f}{m2.geom.height:>10.2f}"
              f"{m2.geom.base_diameter:>10.2f}{m2.steel_mass():>12,.0f}"
              f"{r2.governing_util:>11.3f}{lam:>12.1f}"
              f"   span/{r2.deflection_ratio:<7.0f}  mode at z/H = {where}")

    # State what the numbers actually show rather than a preconception.
    lo, hi = study[0], study[-1]
    print(f"\n  Steel mass rises {hi['mass'] / lo['mass']:.2f}x from the "
          f"{lo['frac']:.2f} cap to the {hi['frac']:.2f} dome because there is "
          f"\n  more surface to clad and carry.  The dome also becomes "
          f"*more* flexible, not\n  less (span/{lo['defl']:.0f} -> "
          f"span/{hi['defl']:.0f}): past the hemisphere the added surface "
          f"brings load\n  on faster than the added rise brings stiffness, and "
          f"the base ring stops\n  growing once the cut passes the equator.")
    print(f"\n  The buckling results split into two regimes.  Up to the "
          f"hemisphere the\n  critical mode sits at the crown (z/H ~ "
          f"{study[0]['mode_z']}) under the snow combination and\n  lambda_cr "
          f"is essentially constant at {study[0]['lam']:.0f}, because the "
          f"crown cap is geometrically\n  identical in all of them.  For the "
          f"deeper domes the governing case switches to\n  wind "
          f"({hi['ref']}) and the mode migrates down to the windward\n  flank "
          f"(z/H ~ {hi['mode_z']}) at lambda_cr = {hi['lam']:.0f}.  Truncation "
          f"height is therefore a cost\n  decision below the hemisphere and a "
          f"wind-stability decision above it.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
