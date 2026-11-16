"""
Parameter study: how frequency, radius and depth change a geodesic dome.

The question a designer asks before anything else is *which* dome to build:
a 2V dome has a few long struts, an 8V dome hundreds of short ones.  This
module answers it by regenerating the dome for every combination of

    frequency  nu   (Class I breakdown, 2V, 3V, ...)
    radius     R    [m]
    fraction        retained part of the sphere diameter (0.5 = hemisphere)

and running the complete GeoDomeLab analysis on each one: the same direct
stiffness solve, load combinations, AISC 360 member checks, deflection check
and linear buckling analysis as a single model (``postprocess.run_analysis``),
optionally preceded by the fully-stressed-design auto-sizing of the strut
families (``optimize.auto_size``).  Nothing is re-implemented here, so every
study point is exactly what the user would get by building that dome by hand.

What carries over from the base model
-------------------------------------
Everything that does not refer to a node number: material, support *type*
(pinned, roller, elastic springs with their stiffness), load data (pressures,
wind, crown point load, uniform temperature change), load combinations
including the T factor, K, phi factors and the deflection limit.  Inputs that
are tied to node numbers of the base geometry (support settlements, explicit
nodal loads, CUSTOM per-node supports) have no meaning on a regenerated mesh:
settlements and nodal loads are left out (and a message says so), and a
CUSTOM support table raises ``ValueError``.

Sections: with ``auto_size=False`` every strut family of every dome gets one
uniform section (``section`` if given, else the base model's section when it
uses one throughout, else the heaviest section used in the base model), so
the rows compare geometries with the same member.  With ``auto_size=True``
each dome starts from that section and is then sized family by family.

Output columns (one row per study point)
----------------------------------------
* counts: nodes, struts, panels, strut families, supported nodes.  For an
  even frequency cut at the equator the Class I hemisphere has

      V = (10 nu^2 + 5 nu + 2) / 2,   E = (30 nu^2 + 5 nu) / 2,   F = 10 nu^2

  (half of the full-sphere counts 10 nu^2 + 2, 30 nu^2, 20 nu^2 plus the
  5 nu nodes and struts of the equator, which both halves share; Euler's
  V - E + F = 1 for a disc).  ``tests/test_study.py`` checks this.
* steel mass  M = sum_k rho A_k L_k                               [kg]
* mass per plan area  M / A_plan, A_plan = sum_f |((b - a) x (c - a))_z| / 2
  over the cladding triangles (the area of the base polygon for a dome)
* max utilisation  max over strength combinations and members of
  |N| / (phi Pn), AISC 360-16 D2 (tension) and E3 (flexural buckling)
* governing combination (and member)
* lambda_cr: smallest root of det(K + lambda Kg(N_ref)) = 0, N_ref = member
  forces of the governing combination (linearised buckling, McGuire,
  Gallagher & Ziemian 2000, Ch. 9).  NaN when not requested or when the
  dome has no compressive buckling mode.
* span / deflection = D_base / max_i |u_i| under the serviceability
  combination, against the deflection limit
* status (OK, OVERSTRESSED, DEFLECTION, UNSTABLE, ERROR) and the two solver
  identities (equilibrium residual and strain-energy mismatch, both ~1e-12)
* wall-clock time of the point [s] (geometry + sizing + analysis + buckling)

Independent checks of a study (used in the tests)
-------------------------------------------------
Dimensional analysis gives exact scaling laws for a fixed section when every
load is a pressure (no self weight, no point load): the nodal loads scale as
R^2, the bar stiffnesses EA/L as 1/R, so

    u ~ R^3,   N ~ R^2,   Kg ~ N/L ~ R,   lambda_cr ~ K/Kg ~ R^-2,
    span/deflection ~ R/R^3 = R^-2,   mass ~ R,   mass/plan area ~ 1/R.

A study at R and 2R must reproduce these ratios to round-off.  The 2V mass
has a closed form from the chord factors 2 sin(arctan(2)/4) (30 struts) and
2 sin(18 deg) (35 struts) of the hemisphere (Kenner 1976; Clinton 1971).

Reading the results
-------------------
With one fixed section, a higher frequency means shorter, stockier struts
(utilisation falls quickly because phi Pn grows with 1/(KL/r)^2), but more
steel, and a *lower* lambda_cr: a pin-jointed single-layer lattice gets its
out-of-plane stiffness only from the facet geometry, and the local rise of a
node over its neighbours shrinks like h ~ L^2 / (2R).  A rough estimate from
the shallow two-bar (von Mises) truss, P_cr ~ EA (h/L)^3 with a nodal load
P ~ q L^2, gives lambda_cr ~ EA L / (q R^3), i.e. roughly 1/nu at fixed
section (an order-of-magnitude argument, not a design formula; self weight
and the global mode shape change the exponent).  That is why the study
reports lambda_cr next to the utilisation; instability is the main design
problem of single-layer reticulated shells (Gioncu 1995).

References
----------
* W. McGuire, R. H. Gallagher, R. D. Ziemian (2000). Matrix Structural
  Analysis, 2nd ed., Wiley. Ch. 3 (direct stiffness), Ch. 9 (geometric
  stiffness and elastic critical load).
* ANSI/AISC 360-16, Specification for Structural Steel Buildings, D2, E3.
* ASCE/SEI 7-16, Minimum Design Loads, Ch. 2 (load combinations).
* H. Kenner (1976). Geodesic Math and How to Use It. Univ. of California
  Press (chord factors and breakdown counts).
* J. D. Clinton (1971). Advanced Structural Geometry Studies, Part I.
  NASA CR-1734 (chord factor definition).
* V. Gioncu (1995). Buckling of reticulated shells: state-of-the-art.
  Int. J. Space Structures 10(1), 1-46.
"""

from __future__ import annotations

import copy
import csv
import math
import time
from dataclasses import fields

import numpy as np

from . import geometry as geo
from . import optimize as opt
from . import sections as sx
from .model import DomeModel
from .postprocess import run_analysis


# ------------------------------------------------------------- columns -----
# (key, short header, unit).  The order is the CSV column order.
STUDY_COLUMNS = [
    ("frequency", "nu", "V"),
    ("radius_m", "R", "m"),
    ("fraction_requested", "frac req", "-"),
    ("fraction", "frac", "-"),
    ("span_m", "span", "m"),
    ("rise_m", "rise", "m"),
    ("nodes", "nodes", "-"),
    ("struts", "struts", "-"),
    ("panels", "panels", "-"),
    ("families", "fam", "-"),
    ("supports", "supp", "-"),
    ("strut_min_m", "Lmin", "m"),
    ("strut_max_m", "Lmax", "m"),
    ("total_strut_length_m", "sum L", "m"),
    ("sections", "sections", ""),
    ("mass_kg", "mass", "kg"),
    ("mass_per_plan_kg_m2", "mass/plan", "kg/m2"),
    ("max_util", "max util", "-"),
    ("governing_combination", "governing", ""),
    ("governing_member", "member", "-"),
    ("members_over", "n>1", "-"),
    ("lambda_cr", "lambda_cr", "-"),
    ("buckling_case", "buckling case", ""),
    ("max_service_disp_mm", "max |u|", "mm"),
    ("crown_disp_mm", "crown uz", "mm"),
    ("span_over_deflection", "span/defl", "-"),
    ("deflection_limit", "limit", "-"),
    ("status", "status", ""),
    ("auto_sized", "sized", ""),
    ("sizing_converged", "converged", ""),
    ("sizing_iterations", "iter", "-"),
    ("equilibrium_residual", "equil", "-"),
    ("energy_error", "energy", "-"),
    ("time_s", "time", "s"),
    ("notes", "notes", ""),
]
COLUMN_KEYS = [c[0] for c in STUDY_COLUMNS]

# Long names for axis labels.
COLUMN_LABELS = {
    "frequency": "Frequency nu",
    "radius_m": "Radius R [m]",
    "fraction_requested": "Dome fraction (requested)",
    "fraction": "Dome fraction",
    "span_m": "Span [m]",
    "rise_m": "Rise [m]",
    "nodes": "Nodes",
    "struts": "Struts",
    "panels": "Panels",
    "families": "Strut families",
    "supports": "Supported nodes",
    "strut_min_m": "Shortest strut [m]",
    "strut_max_m": "Longest strut [m]",
    "total_strut_length_m": "Total strut length [m]",
    "mass_kg": "Steel mass [kg]",
    "mass_per_plan_kg_m2": "Mass / plan area [kg/m2]",
    "max_util": "Max utilisation",
    "members_over": "Members over 1.0",
    "lambda_cr": "Buckling factor lambda_cr",
    "max_service_disp_mm": "Max service displacement [mm]",
    "crown_disp_mm": "Crown displacement [mm]",
    "span_over_deflection": "Span / deflection",
    "time_s": "Time [s]",
}

SWEEP_KEYS = ("frequency", "radius_m", "fraction_requested")


class StudyRows(list):
    """The rows of a study (a plain ``list`` of dicts) plus two attributes.

    ``messages`` -- study-level notes (inputs left out, cancellation, errors);
    ``settings`` -- the options the study ran with.
    """

    def __init__(self, rows=(), messages=None, settings=None):
        super().__init__(rows)
        self.messages = list(messages or [])
        self.settings = dict(settings or {})


# ------------------------------------------------------------ helpers ------
def _as_list(values, name):
    if values is None:
        return None
    if np.isscalar(values):
        values = [values]
    out = list(values)
    if not out:
        raise ValueError(f"{name}: give at least one value.")
    return out


def _check_frequencies(frequencies):
    out = []
    for f in _as_list(frequencies, "frequencies"):
        try:
            v = float(f)
        except (TypeError, ValueError):
            raise ValueError(f"frequencies: {f!r} is not a number.") from None
        if not math.isfinite(v) or v != round(v) or v < 1:
            raise ValueError(f"frequencies: {f!r} is not a whole number >= 1 "
                             f"(the Class I frequency is an integer).")
        out.append(int(round(v)))
    return out


def _check_radii(radii):
    out = []
    for r in radii:
        v = float(r)
        if not (math.isfinite(v) and v > 0.0):
            raise ValueError(f"radii: {r!r} must be a positive length in m.")
        out.append(v)
    return out


def _check_fractions(fractions):
    out = []
    for f in fractions:
        v = float(f)
        if not (math.isfinite(v) and 0.0 < v <= 1.0):
            raise ValueError(f"fractions: {f!r} must lie in (0, 1] "
                             f"(0.5 = hemisphere).")
        out.append(v)
    return out


def _check_base(base_model):
    if base_model.support_type == geo.CUSTOM:
        raise ValueError(
            "The base model uses CUSTOM per-node supports. Node numbers change "
            "when the dome is regenerated, so they cannot be carried over to "
            "the study. Choose pinned, roller or elastic supports for the "
            "study.")
    if base_model.support_type not in geo.SUPPORT_TYPES:
        raise ValueError(f"Unknown support type {base_model.support_type!r}.")


def uniform_section(base_model, section=None):
    """The single section used for every family, and why.

    Returns ``(section_name, explanation)``.
    """
    if section is not None:
        if section not in sx.SECTIONS:
            raise ValueError(f"Unknown section {section!r}; choose one of the "
                             f"catalogue names (geodome.sections).")
        return section, f"every strut family uses {section} (given)"
    names = list(dict.fromkeys(base_model.group_sections))
    if len(names) == 1:
        return names[0], (f"every strut family uses the base model's "
                          f"section {names[0]}")
    heaviest = max(names, key=lambda n: sx.SECTIONS[n].A)
    return heaviest, (f"the base model uses {len(names)} different sections; "
                      f"every strut family uses the heaviest of them, "
                      f"{heaviest}")


def _node_specific_load_fields(loads):
    """LoadSpec list fields whose entries refer to node numbers."""
    out = {}
    for f in fields(loads):
        val = getattr(loads, f.name)
        if isinstance(val, list) and val and all(
                isinstance(e, dict) and "node" in e for e in val):
            out[f.name] = len(val)
    return out


def study_model(base_model, frequency, radius, fraction, section=None):
    """The :class:`DomeModel` of one study point.

    The geometry is regenerated with ``generate_dome(frequency, radius,
    fraction)``; every other setting is a deep copy of ``base_model`` (see the
    module docstring for what is left out).  Every strut family gets the
    uniform section of :func:`uniform_section`.
    """
    _check_base(base_model)
    sec, _ = uniform_section(base_model, section)
    geom = geo.generate_dome(int(frequency), float(radius), float(fraction))

    kw = {f.name: copy.deepcopy(getattr(base_model, f.name))
          for f in fields(base_model) if f.name != "geom"}
    loads = kw["loads"]
    for name in _node_specific_load_fields(loads):
        setattr(loads, name, [])
    kw.update(geom=geom, group_sections=[sec] * geom.n_strut_types(),
              loads=loads, custom_supports={})
    return DomeModel(**kw)


def model_from_row(base_model, row):
    """Rebuild the model behind a study row (geometry *and* sections)."""
    m = study_model(base_model, row["frequency"], row["radius_m"],
                    row["fraction_requested"])
    secs = row.get("group_sections")
    if secs:
        m.group_sections = list(secs)
    return m


def _sections_label(group_sections):
    names = list(dict.fromkeys(group_sections))
    if len(names) == 1:
        return f"all: {names[0]}"
    return "; ".join(f"{chr(ord('A') + g) if g < 26 else f'S{g}'}: {s}"
                     for g, s in enumerate(group_sections))


def _geometry_notes(geom):
    keep = ("snapped", "stepped", "support set")
    return [n for n in geom.notes if any(k in n for k in keep)]


def _status(row):
    if row.get("status") == "ERROR":
        return "ERROR"
    bad = []
    if not row["max_util"] <= 1.0:
        bad.append("OVERSTRESSED")
    if not row["span_over_deflection"] >= row["deflection_limit"]:
        bad.append("DEFLECTION")
    lam = row["lambda_cr"]
    if math.isfinite(lam) and lam < 1.0:
        bad.append("UNSTABLE")
    return "OK" if not bad else " + ".join(bad)


def _nan_row(frequency, radius, fraction):
    row = {k: float("nan") for k in COLUMN_KEYS}
    row.update(frequency=int(frequency), radius_m=float(radius),
               fraction_requested=float(fraction), sections="",
               governing_combination="", buckling_case="", status="ERROR",
               auto_sized=False, sizing_converged=None, sizing_iterations=None,
               notes="", group_sections=[])
    return row


def _same_geometry(a, b):
    return (a.n_nodes == b.n_nodes and a.n_members == b.n_members
            and np.array_equal(a.members, b.members)
            and np.array_equal(a.nodes, b.nodes))


# --------------------------------------------------------------- study -----
def parameter_study(base_model, frequencies, radii=None, fractions=None,
                    auto_size=False, do_buckling=True, progress=None,
                    section=None, target_util=0.95, max_iter=12,
                    enforce_deflection=True, cancel=None):
    """Analyse a grid of geodesic domes built from ``base_model``'s settings.

    Parameters
    ----------
    base_model : DomeModel
        Supplies material, supports, loads, combinations and design factors.
    frequencies : int or sequence of int
        Class I frequencies to study (whole numbers >= 1).
    radii : float or sequence, optional
        Sphere radii [m]; default: the base model's radius.
    fractions : float or sequence, optional
        Dome fractions in (0, 1]; default: the base model's fraction.  The
        generator snaps a fraction onto a ring of nodes, so the row keeps the
        requested and the actual value.  Requests that snap onto the same
        geometry are analysed once and the row is reused (noted).
    auto_size : bool
        Size every dome's strut families by fully-stressed design
        (``optimize.auto_size(target_util, max_iter, enforce_deflection)``)
        before the final analysis.
    do_buckling : bool
        Compute lambda_cr (one linear buckling eigensolve per point).
    progress : callable(fraction, message), optional
        Progress callback, as for ``run_analysis``.
    section : str, optional
        Uniform section for every family (see :func:`uniform_section`).
    cancel : callable() -> bool, optional
        Polled before every point; returning True stops the study and the
        rows done so far are returned (with a message).

    Returns
    -------
    StudyRows
        A list of row dicts (keys :data:`COLUMN_KEYS` plus ``group_sections``,
        the per-family section list), ordered radius -> fraction -> frequency,
        with ``.messages`` and ``.settings``.
    """
    _check_base(base_model)
    freqs = _check_frequencies(frequencies)
    g0 = base_model.geom
    if radii is None:
        if g0.is_custom or not (np.isfinite(g0.radius) and g0.radius > 0):
            raise ValueError("The base model is a custom truss with no sphere "
                             "radius; give the radii to study.")
        radii = [float(g0.radius)]
    radii = _check_radii(_as_list(radii, "radii"))
    if fractions is None:
        if g0.is_custom or not np.isfinite(g0.dome_fraction):
            raise ValueError("The base model is a custom truss with no dome "
                             "fraction; give the fractions to study.")
        fractions = [float(g0.dome_fraction)]
    fractions = _check_fractions(_as_list(fractions, "fractions"))
    sec, sec_note = uniform_section(base_model, section)
    if auto_size and not (0.0 < float(target_util) <= 1.0):
        raise ValueError("target_util must lie in (0, 1].")

    messages = []
    dropped = _node_specific_load_fields(base_model.loads)
    if dropped:
        what = ", ".join(f"{n} {k.replace('_', ' ')}"
                         for k, n in dropped.items())
        messages.append(
            f"Node-numbered input does not carry over to a regenerated dome, "
            f"so it was left out of every study point: {what}.")
    messages.append(("Sections: " + sec_note
                     + (", then auto-sized by fully-stressed design."
                        if auto_size else ".")))

    settings = {
        "frequencies": freqs, "radii": radii, "fractions": fractions,
        "auto_size": bool(auto_size), "do_buckling": bool(do_buckling),
        "section": sec, "target_util": float(target_util),
        "max_iter": int(max_iter), "enforce_deflection": bool(enforce_deflection),
        "material": base_model.material_name,
        "support_type": base_model.support_type,
        "combinations": [c.name for c in base_model.combinations],
    }

    def report(f, msg):
        if progress is not None:
            progress(float(min(max(f, 0.0), 1.0)), str(msg))

    points = [(nu, R, fr) for R in radii for fr in fractions for nu in freqs]
    n = len(points)
    rows = StudyRows(messages=messages, settings=settings)
    cache = []                                 # [(geom, row)]

    for i, (nu, R, fr) in enumerate(points):
        if cancel is not None and cancel():
            rows.messages.append(f"Stopped by the user after {i} of {n} "
                                 f"study points.")
            break
        lo, hi = i / n, (i + 1) / n
        tag = f"{nu}V, R = {R:g} m, fraction {fr:g} ({i + 1}/{n})"
        report(lo, f"Study point {tag}")
        t0 = time.perf_counter()
        try:
            model = study_model(base_model, nu, R, fr, sec)
            hit = next((r for g, r in cache if _same_geometry(g, model.geom)),
                       None)
            if hit is not None:
                row = copy.deepcopy(hit)
                row["fraction_requested"] = float(fr)
                row["notes"] = "; ".join(
                    s for s in (f"same geometry as the request "
                                f"fraction {hit['fraction_requested']:g} "
                                f"(snapped); results reused",
                                row["notes"]) if s)
                row["time_s"] = time.perf_counter() - t0
                rows.append(row)
                report(hi, f"Done {tag}")
                continue
            row = _analyse_point(
                model, nu, R, fr, auto_size=auto_size, do_buckling=do_buckling,
                target_util=target_util, max_iter=max_iter,
                enforce_deflection=enforce_deflection,
                report=lambda f, msg: report(lo + (hi - lo) * f,
                                             f"{tag}: {msg}"))
            row["time_s"] = time.perf_counter() - t0
            cache.append((model.geom, row))
        except Exception as exc:                     # keep the other points
            row = _nan_row(nu, R, fr)
            row["notes"] = f"{type(exc).__name__}: {exc}"
            row["time_s"] = time.perf_counter() - t0
            rows.messages.append(f"{tag} failed: {row['notes']}")
        rows.append(row)
        report(hi, f"Done {tag}")

    report(1.0, f"Parameter study complete ({len(rows)} of {n} points).")
    return rows


def _analyse_point(model, nu, R, fr, auto_size, do_buckling, target_util,
                   max_iter, enforce_deflection, report):
    g = model.geom
    sizing = None
    if auto_size:
        sizing = opt.auto_size(
            model, target_util=target_util, max_iter=max_iter,
            enforce_deflection=enforce_deflection,
            progress=lambda f, msg: report(0.7 * f, msg))
    res = run_analysis(model, n_buckling_modes=1, do_buckling=do_buckling,
                       progress=lambda f, msg: report(
                           (0.7 if auto_size else 0.0)
                           + (0.3 if auto_size else 1.0) * f, msg))

    L = g.strut_lengths
    lam = float("nan")
    bcase = ""
    notes = _geometry_notes(g)
    if res.buckling is not None:
        bcase = res.buckling.reference_case
        if len(res.buckling.load_factors):
            lam = float(res.buckling.load_factors[0])
        else:
            notes.append("no compressive buckling mode under the reference "
                         "combination")
    for msg in res.messages:
        if msg.startswith("Buckling analysis failed"):
            notes.append(msg)
    if sizing is not None:
        notes += sizing.messages

    row = {
        "frequency": int(nu),
        "radius_m": float(R),
        "fraction_requested": float(fr),
        "fraction": float(g.dome_fraction),
        "span_m": float(g.base_diameter),
        "rise_m": float(g.height),
        "nodes": int(g.n_nodes),
        "struts": int(g.n_members),
        "panels": int(g.n_faces),
        "families": int(g.n_strut_types()),
        "supports": int(len(model.support_nodes())),
        "strut_min_m": float(L.min()),
        "strut_max_m": float(L.max()),
        "total_strut_length_m": float(L.sum()),
        "sections": _sections_label(model.group_sections),
        "group_sections": list(model.group_sections),
        "mass_kg": float(model.steel_mass()),
        "mass_per_plan_kg_m2": float(model.mass_per_plan_area()),
        "max_util": float(res.governing_util),
        "governing_combination": res.governing_case,
        "governing_member": int(res.governing_member),
        "members_over": int(np.sum(res.util_env > 1.0)),
        "lambda_cr": lam,
        "buckling_case": bcase,
        "max_service_disp_mm": float(res.max_service_disp * 1e3),
        "crown_disp_mm": float(res.crown_disp * 1e3),
        "span_over_deflection": float(res.deflection_ratio),
        "deflection_limit": float(model.deflection_limit),
        "auto_sized": bool(auto_size),
        "sizing_converged": (bool(sizing.converged) if sizing is not None
                             else None),
        "sizing_iterations": (int(sizing.iterations) if sizing is not None
                              else None),
        "equilibrium_residual": float(max(r.equilibrium_residual
                                          for r in res.cases.values())),
        "energy_error": float(max(r.energy_error for r in res.cases.values())),
        "time_s": float("nan"),
        "notes": "; ".join(notes),
    }
    row["status"] = _status(row)
    return row


# --------------------------------------------------------------- output ----
def best_design(rows, objective="mass_kg", min_lambda_cr=None):
    """The row with the smallest ``objective`` among the acceptable designs.

    Acceptable = status "OK" (utilisation <= 1, deflection within the limit,
    lambda_cr >= 1 when computed) and, when ``min_lambda_cr`` is given,
    lambda_cr >= min_lambda_cr.  Returns ``None`` when no row qualifies.
    """
    ok = []
    for r in rows:
        if r.get("status") != "OK":
            continue
        val = r.get(objective)
        if val is None or not math.isfinite(float(val)):
            continue
        if min_lambda_cr is not None:
            lam = r.get("lambda_cr", float("nan"))
            if not (math.isfinite(lam) and lam >= float(min_lambda_cr)):
                continue
        ok.append(r)
    return min(ok, key=lambda r: float(r[objective])) if ok else None


def _csv_value(v):
    if v is None:
        return ""
    if isinstance(v, (bool, np.bool_)):
        return "yes" if v else "no"
    if isinstance(v, (float, np.floating)):
        return repr(float(v))                 # round-trips exactly
    if isinstance(v, (list, tuple)):
        return "; ".join(str(x) for x in v)
    return str(v)


def study_to_csv(rows, path):
    """Write the study rows to ``path`` (UTF-8 CSV, one row per point).

    Columns follow :data:`STUDY_COLUMNS` (units are in the key names).
    Floats are written with ``repr`` so they read back exactly; booleans as
    yes/no; ``None`` as an empty cell.
    """
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(COLUMN_KEYS)
        for r in rows:
            w.writerow([_csv_value(r.get(k)) for k in COLUMN_KEYS])
    return path


def _fmt(key, v):
    if v is None:
        return "-"
    if isinstance(v, str):
        return v
    if isinstance(v, (bool, np.bool_)):
        return "yes" if v else "no"
    v = float(v)
    if not math.isfinite(v):
        return "n/a" if math.isnan(v) else "inf"
    if key in ("frequency", "nodes", "struts", "panels", "families",
               "supports", "members_over", "governing_member",
               "sizing_iterations"):
        return f"{int(v):d}"
    if key in ("fraction", "fraction_requested", "max_util"):
        return f"{v:.3f}"
    if key in ("span_over_deflection", "mass_kg", "deflection_limit"):
        return f"{v:.0f}"
    if key in ("equilibrium_residual", "energy_error"):
        return f"{v:.1e}"
    if key == "lambda_cr":
        return f"{v:.2f}"
    if key == "time_s":
        return f"{v:.2f}"
    return f"{v:.2f}"


DEFAULT_TABLE_COLUMNS = ("frequency", "radius_m", "fraction", "nodes",
                         "struts", "families", "mass_kg",
                         "mass_per_plan_kg_m2", "max_util",
                         "governing_combination", "lambda_cr",
                         "span_over_deflection", "status", "time_s")


def study_table(rows, columns=DEFAULT_TABLE_COLUMNS):
    """Plain-text table of the study (list of lines, ASCII)."""
    head = {k: h for k, h, _ in STUDY_COLUMNS}
    unit = {k: u for k, _, u in STUDY_COLUMNS}
    cells = [[_fmt(k, r.get(k)) for k in columns] for r in rows]
    titles = [head.get(k, k) for k in columns]
    units = [(f"[{unit[k]}]" if unit.get(k) not in (None, "", "-") else "")
             for k in columns]
    width = [max([len(t), len(u)] + [len(c[j]) for c in cells])
             for j, (t, u) in enumerate(zip(titles, units))]
    line = lambda vals: "  ".join(v.rjust(w) for v, w in zip(vals, width))  # noqa: E731
    out = [line(titles), line(units), "  ".join("-" * w for w in width)]
    out += [line(c) for c in cells]
    return out


def _series_keys(rows, x):
    varying = []
    for k in SWEEP_KEYS:
        if k == x or (x == "fraction" and k == "fraction_requested"):
            continue
        if len({r.get(k) for r in rows}) > 1:
            varying.append(k)
    return varying


def _series_label(key_vals):
    parts = []
    for k, v in key_vals:
        if k == "frequency":
            parts.append(f"{int(v)}V")
        elif k == "radius_m":
            parts.append(f"R = {v:g} m")
        else:
            parts.append(f"fraction {v:g}")
    return ", ".join(parts)


def plot_study(fig, rows, x="frequency",
               ys=("mass_kg", "max_util", "lambda_cr", "span_over_deflection"),
               logy=()):
    """Plot study columns ``ys`` against ``x`` on ``fig`` (one panel each).

    One line per combination of the other swept parameters.  Filled markers
    are acceptable designs (status OK), hollow markers fail a check.  Dashed
    reference lines: utilisation 1.0, lambda_cr 1.0 and the span/deflection
    limit.  ``fig`` may be ``None`` (a new matplotlib Figure is made; pyplot
    is not used).  Returns ``(fig, axes)``.
    """
    from matplotlib.figure import Figure
    from matplotlib.ticker import MaxNLocator

    if isinstance(ys, str):
        ys = (ys,)
    ys = list(ys)
    if fig is None:
        fig = Figure(figsize=(7.0, 2.4 * max(len(ys), 1)))
    fig.clear()
    good = [r for r in rows if r.get("status") != "ERROR"]
    axes = list(fig.subplots(len(ys), 1, sharex=True, squeeze=False)[:, 0])
    skeys = _series_keys(good, x)
    groups = {}
    for r in good:
        groups.setdefault(tuple(r.get(k) for k in skeys), []).append(r)

    for ax, yk in zip(axes, ys):
        for gk, grp in groups.items():
            grp = sorted(grp, key=lambda r: float(r[x]))
            xs = np.array([float(r[x]) for r in grp])
            yv = np.array([float(r[yk]) for r in grp])
            ok = np.array([r.get("status") == "OK" for r in grp])
            label = _series_label(list(zip(skeys, gk))) if skeys else None
            (ln,) = ax.plot(xs, yv, "-", lw=1.4, label=label)
            col = ln.get_color()
            if ok.any():
                ax.plot(xs[ok], yv[ok], "o", color=col, ms=5)
            if (~ok).any():
                ax.plot(xs[~ok], yv[~ok], "o", mfc="none", mec=col, ms=6)
        if yk in ("max_util", "lambda_cr"):
            ax.axhline(1.0, ls="--", lw=0.9, color="0.5")
        if yk == "span_over_deflection":
            lims = {float(r["deflection_limit"]) for r in good
                    if math.isfinite(float(r.get("deflection_limit", "nan")))}
            for lim in lims:
                ax.axhline(lim, ls="--", lw=0.9, color="0.5")
        if yk in logy:
            ax.set_yscale("log")
        ax.set_ylabel(COLUMN_LABELS.get(yk, yk))
        ax.grid(True, alpha=0.3)
    axes[-1].set_xlabel(COLUMN_LABELS.get(x, x))
    if x in ("frequency", "nodes", "struts", "panels", "families"):
        axes[-1].xaxis.set_major_locator(MaxNLocator(integer=True))
    if skeys and len(groups) > 1:
        axes[0].legend(fontsize="small")
    title = "Parameter study"
    if any(r.get("auto_sized") for r in good):
        title += " (auto-sized)"
    axes[0].set_title(title + "  -  filled: design OK, hollow: fails a check",
                      fontsize="medium")
    try:
        fig.tight_layout()
    except Exception:                                  # pragma: no cover
        pass
    return fig, axes
