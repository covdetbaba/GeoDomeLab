"""
Run the full analysis suite over a set of load combinations and envelope it.

The :class:`AnalysisResults` object produced here is what the GUI tables, the
plots and the exported report all read from, so the solving logic stays in one
place and the presentation layers stay thin.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from . import loads as ld
from . import solver as slv
from . import buckling as bk


@dataclass
class AnalysisResults:
    model: object
    cases: dict = field(default_factory=dict)        # name -> StaticResult
    order: list = field(default_factory=list)        # combination names, in order

    caps: dict = field(default_factory=dict)         # member capacities
    util: dict = field(default_factory=dict)         # name -> utilisation array

    # envelopes over the strength combinations
    N_max: np.ndarray = None                         # most tensile
    N_min: np.ndarray = None                         # most compressive
    N_max_case: np.ndarray = None
    N_min_case: np.ndarray = None
    util_env: np.ndarray = None
    util_case: np.ndarray = None

    governing_case: str = ""
    governing_member: int = -1
    governing_util: float = 0.0

    service_case: str = ""
    max_service_disp: float = 0.0
    max_service_node: int = -1
    crown_disp: float = 0.0
    deflection_ratio: float = float("inf")
    deflection_ok: bool = True

    # Linear buckling.  ``buckling`` is the full result (all requested
    # modes) for the critical combination, the strength combination with
    # the lowest positive critical load factor.  ``buckling_cases`` holds the
    # result of every strength combination, so ``lambda_cr_by_case`` gives
    # the first load factor of each (nan when it has no compressive mode).
    buckling: object = None
    buckling_cases: dict = field(default_factory=dict)
    messages: list = field(default_factory=list)

    # ------------------------------------------------------------ helpers --
    def strength_cases(self):
        return [n for n in self.order if not self._combo(n).serviceability]

    def _combo(self, name):
        for c in self.model.combinations:
            if c.name == name:
                return c
        raise KeyError(name)

    def case(self, name):
        return self.cases[name]

    def worst_case_result(self):
        return self.cases.get(self.governing_case) or next(iter(self.cases.values()))

    @property
    def lambda_cr_by_case(self):
        """First critical load factor of every analysed combination."""
        out = {}
        for name, b in self.buckling_cases.items():
            out[name] = (float(b.load_factors[0]) if len(b.load_factors)
                         else float("nan"))
        return out

    @property
    def lambda_cr(self):
        """Lowest critical load factor over the strength combinations."""
        if self.buckling is not None and len(self.buckling.load_factors):
            return float(self.buckling.load_factors[0])
        return float("nan")

    @property
    def lambda_cr_case(self):
        """Combination that gives :attr:`lambda_cr`."""
        return self.buckling.reference_case if self.buckling is not None \
            else ""


def run_analysis(model, n_buckling_modes=6, do_buckling=True, progress=None):
    """Solve every load combination, envelope the results and check the design.

    ``progress`` is an optional ``callable(fraction, message)`` used by the GUI
    to drive its progress bar from the worker thread.
    """
    def report(f, msg):
        if progress is not None:
            progress(float(f), str(msg))

    report(0.02, "Assembling load cases...")
    basic = ld.basic_load_vectors(model)

    res = AnalysisResults(model=model)
    combos = list(model.combinations)
    if not combos:
        raise ValueError("No load combinations are defined.")

    nstep = max(len(combos), 1)
    for i, c in enumerate(combos):
        report(0.05 + 0.55 * i / nstep, f"Solving {c.name} ...")
        F = ld.combine(basic, c)
        r = slv.solve_model(model, F, name=c.name,
                            T_factor=float(getattr(c, "T", 0.0)))
        res.cases[c.name] = r
        res.order.append(c.name)

    report(0.65, "Checking member capacities...")
    caps = bk.member_capacities(model)
    res.caps = caps

    nm = model.geom.n_members
    N_max = np.full(nm, -np.inf)
    N_min = np.full(nm, np.inf)
    N_max_case = np.zeros(nm, dtype=int)
    N_min_case = np.zeros(nm, dtype=int)
    util_env = np.zeros(nm)
    util_case = np.zeros(nm, dtype=int)

    strength = [c for c in combos if not c.serviceability]
    for idx, c in enumerate(strength):
        r = res.cases[c.name]
        u, _ = bk.utilisation(model, r.N, caps)
        res.util[c.name] = u

        hi = r.N > N_max
        N_max[hi] = r.N[hi]
        N_max_case[hi] = idx
        lo = r.N < N_min
        N_min[lo] = r.N[lo]
        N_min_case[lo] = idx
        up = u > util_env
        util_env[up] = u[up]
        util_case[up] = idx

    for c in combos:                                  # also store service ratios
        if c.serviceability and c.name not in res.util:
            res.util[c.name], _ = bk.utilisation(model, res.cases[c.name].N, caps)

    res.N_max, res.N_min = N_max, N_min
    res.N_max_case, res.N_min_case = N_max_case, N_min_case
    res.util_env, res.util_case = util_env, util_case

    if not strength:
        res.messages.append(
            "No strength combination (every row is a service case), so no "
            "member was checked. Untick Service on at least one combination.")
    if len(strength):
        gm = int(np.argmax(util_env))
        res.governing_member = gm
        res.governing_util = float(util_env[gm])
        res.governing_case = strength[int(util_case[gm])].name

    # ------------------------------------------------------ serviceability --
    report(0.72, "Checking serviceability...")
    service = [c for c in combos if c.serviceability]
    if service:
        sc = service[0]
    else:
        sc = combos[0]
    res.service_case = sc.name
    rs = res.cases[sc.name]
    node, dmax = rs.max_displacement()
    res.max_service_node, res.max_service_disp = node, dmax
    res.crown_disp = float(abs(rs.u[3 * model.geom.crown_node + 2]))
    span = model.geom.base_diameter
    res.deflection_ratio = span / max(dmax, 1e-12)
    res.deflection_ok = res.deflection_ratio >= model.deflection_limit

    # ----------------------------------------------------------- buckling --
    # Every strength combination is a design load, so the critical one is
    # the combination with the lowest positive load factor; it need not be
    # the combination that governs the member utilisation.
    if do_buckling:
        names = [c.name for c in strength] or list(res.order)
        for k, name in enumerate(names):
            report(0.80 + 0.12 * k / max(len(names), 1),
                   f"Linear buckling eigenanalysis ({name})...")
            try:
                b = bk.linear_buckling(model, res.cases[name].N,
                                       n_buckling_modes)
            except Exception as exc:
                res.messages.append(
                    f"Buckling analysis failed for {name}: {exc}")
                continue
            b.reference_case = name
            res.buckling_cases[name] = b
        lam = {n: v for n, v in res.lambda_cr_by_case.items()
               if np.isfinite(v)}
        if lam:
            ref = min(lam, key=lam.get)
        else:
            ref = res.governing_case or res.order[0]
        res.buckling = res.buckling_cases.get(ref)
        if len(lam) > 1 and res.governing_case and \
                ref != res.governing_case:
            res.messages.append(
                f"The lowest critical load factor, {lam[ref]:.2f}, is under "
                f"{ref}, not under the combination that governs the member "
                f"checks ({res.governing_case}).")

    # ------------------------------------------------------------ messages --
    report(0.95, "Collecting diagnostics...")
    worst_equil = max(r.equilibrium_residual for r in res.cases.values())
    worst_energy = max(r.energy_error for r in res.cases.values())
    res.messages.append(
        f"Verification: max equilibrium residual {worst_equil:.2e}, "
        f"max strain-energy mismatch {worst_energy:.2e} (both should be ~1e-12).")

    T_cases = [c.name for c in combos if abs(getattr(c, "T", 0.0)) > 1e-12]
    if model.loads.has_self_straining:
        if T_cases:
            res.messages.append(
                f"Self-straining load T (temperature / settlement) applied in "
                f"{len(T_cases)} combination(s): {', '.join(T_cases)}.")
        else:
            res.messages.append(
                "A temperature change or support settlement is defined, but "
                "no combination has a T factor, so it was not applied.")

    # A basic load that no strength combination uses is silently dropped
    # by the envelope: say so (the live load is the usual case).
    load_names = {"D": "dead load", "S": "snow load",
                  "L": "live load (live pressure, crown load and L nodal "
                       "loads)",
                  "W": "wind load"}
    for key in ("D", "S", "L", "W"):
        if not np.any(basic[key] != 0.0):
            continue
        if any(abs(c.factors()[key]) > 1e-12 for c in strength):
            continue
        art = "an" if key in ("L", "S") else "a"
        res.messages.append(
            f"The {load_names[key]} is defined, but no strength combination "
            f"has {art} {key} factor, so the strength check does not include "
            f"it.")

    mat = model.material
    if not getattr(mat, "aisc_checks_apply", True):
        res.messages.append(
            f"Member checks use AISC 360 (steel). They are not valid for "
            f"{mat.name} ({mat.design_code}); use the results for stiffness "
            f"and forces only.")

    slender = caps["slenderness"]
    if np.max(slender) > 200.0:
        n_over = int(np.sum(slender > 200.0))
        res.messages.append(
            f"{n_over} strut(s) exceed the KL/r = 200 slenderness guidance "
            f"(max {np.max(slender):.0f}).")
    if np.any(~np.isfinite(caps["Q"])) or np.any(caps["Q"] < 1.0):
        n_sl = int(np.sum(caps["Q"] < 1.0))
        res.messages.append(
            f"{n_sl} strut(s) use a locally slender round HSS (D/t > 0.11E/Fy); "
            f"their compressive strength uses the AISC 360-16 E7 effective "
            f"area Ae.")
    if res.governing_util > 1.0:
        res.messages.append(
            f"OVERSTRESSED: member {res.governing_member} reaches "
            f"{res.governing_util:.2f} of capacity under {res.governing_case}.")
    if not res.deflection_ok:
        res.messages.append(
            f"Deflection span/{res.deflection_ratio:.0f} exceeds the "
            f"span/{model.deflection_limit:.0f} limit.")

    report(1.0, "Analysis complete.")
    return res


# ----------------------------------------------------------- summary text ---
def summary_text(res):
    """Multi-line plain-text summary used by the GUI and the HTML report."""
    m = res.model
    g = m.geom
    L = []
    L.append("=" * 74)
    L.append("  SPACE TRUSS -- ANALYSIS SUMMARY" if g.is_custom
             else "  GEODESIC DOME -- ANALYSIS SUMMARY")
    L.append("=" * 74)
    L += ["  " + s for s in m.summary_lines()]
    L.append("")
    L.append("-" * 74)
    L.append("  UNFACTORED LOAD RESULTANTS")
    L.append("-" * 74)
    L += ["  " + s for s in ld.describe_loads(m)]
    L.append("")
    L.append("-" * 74)
    L.append("  LOAD COMBINATION RESULTS")
    L.append("-" * 74)
    L.append(f"  {'Combination':<24s}{'max T [kN]':>12s}{'max C [kN]':>12s}"
             f"{'max |u| [mm]':>14s}{'max util':>10s}")
    n_service = 0
    for name in res.order:
        r = res.cases[name]
        _, dmax = r.max_displacement()
        u = res.util.get(name)
        umax = float(np.max(u)) if u is not None else float("nan")
        service = res._combo(name).serviceability
        n_service += service
        mark = " *" if service else ""
        L.append(f"  {name + mark:<24s}{r.N.max() / 1e3:12.1f}"
                 f"{r.N.min() / 1e3:12.1f}{dmax * 1e3:14.2f}{umax:10.3f}")
    if n_service:
        L.append("  * service combination: its utilisation is shown for "
                 "information only and")
        L.append("    is not part of the strength check.")
    L.append("")
    L.append("-" * 74)
    L.append("  DESIGN CHECK (AISC 360 LRFD)")
    L.append("-" * 74)
    if not getattr(m.material, "aisc_checks_apply", True):
        L.append(f"  Note: AISC 360 is a steel code. These checks are not "
                 f"valid for {m.material.name}.")
    if res.governing_member >= 0:
        i, j = g.members[res.governing_member]
        grp = int(g.strut_group[res.governing_member])
        sec = m.section_of_group(grp)
        cap = res.caps
        L.append(f"  Governing member       #{res.governing_member} "
                 f"(nodes {i}-{j}, family "
                 f"{chr(ord('A') + grp)}, {sec.name})")
        L.append(f"  Governing combination  {res.governing_case}")
        L.append(f"  Utilisation            {res.governing_util:.3f}  "
                 f"-> {'OK' if res.governing_util <= 1.0 else 'OVERSTRESSED'}")
        L.append(f"  Slenderness KL/r       "
                 f"{cap['slenderness'][res.governing_member]:.1f}")
        L.append(f"  phi_c Pn               "
                 f"{cap['phiPn_c'][res.governing_member] / 1e3:.1f} kN")
        L.append(f"  phi_t Pn               "
                 f"{cap['phiPn_t'][res.governing_member] / 1e3:.1f} kN")
        L.append(f"  Members over 1.0       "
                 f"{int(np.sum(res.util_env > 1.0))} of {g.n_members}")
        L.append(f"  Mean utilisation       {float(np.mean(res.util_env)):.3f}")
    L.append("")
    L.append(f"  Serviceability case    {res.service_case}")
    L.append(f"  Max nodal displacement {res.max_service_disp * 1e3:.2f} mm "
             f"at node {res.max_service_node}")
    L.append(f"  Crown vertical disp.   {res.crown_disp * 1e3:.2f} mm")
    L.append(f"  Deflection ratio       span/{res.deflection_ratio:.0f} "
             f"(limit span/{m.deflection_limit:.0f}) -> "
             f"{'OK' if res.deflection_ok else 'EXCEEDED'}")
    L.append("")
    if res.buckling is not None and len(res.buckling.load_factors):
        L.append("-" * 74)
        L.append("  LINEAR BUCKLING  (critical combination: "
                 f"{res.buckling.reference_case})")
        L.append("-" * 74)
        for k, lam in enumerate(res.buckling.load_factors):
            L.append(f"  Mode {k + 1:<2d} critical load factor  {lam:10.3f}")
        L.append(f"  {res.buckling.info}")
        by_case = res.lambda_cr_by_case
        if len(by_case) > 1:
            L.append("  First load factor of each strength combination:")
            for name, lam in by_case.items():
                txt = f"{lam:10.3f}" if np.isfinite(lam) else \
                    "  no compressive mode"
                L.append(f"    {name:<24s}{txt}")
    L.append("")
    if res.messages:
        L.append("-" * 74)
        L.append("  NOTES")
        L.append("-" * 74)
        for msg in res.messages:
            L.append(f"  * {msg}")
    if res.buckling is not None:
        for msg in res.buckling.warnings:
            L.append(f"  * {msg}")
    L.append("=" * 74)
    return "\n".join(L)


# ------------------------------------------------------------- tabulation ---
def _support_nodes(model):
    """Nodes with a restraint or a spring (base nodes for older models)."""
    try:
        return np.asarray(model.support_nodes(), dtype=int)
    except AttributeError:                           # pragma: no cover
        return np.asarray(model.geom.base_nodes, dtype=int)


def member_table(res):
    """Rows of per-member results for the GUI table and CSV export."""
    m = res.model
    g = m.geom
    caps = res.caps
    rows = []
    for k in range(g.n_members):
        i, j = g.members[k]
        grp = int(g.strut_group[k])
        # Envelope over the strength combinations: the largest (most
        # tensile) and the smallest (most compressive) axial force.  A strut
        # that is always in compression has two negative values.
        n_max = res.N_max[k]
        n_min = res.N_min[k]
        rows.append({
            "member": k,
            "node_i": int(i),
            "node_j": int(j),
            "family": chr(ord("A") + grp),
            "section": m.group_sections[grp],
            "length_m": float(g.strut_lengths[k]),
            "KL_r": float(caps["slenderness"][k]),
            "N_max_kN": float(n_max / 1e3),
            "N_min_kN": float(n_min / 1e3),
            "phiPn_t_kN": float(caps["phiPn_t"][k] / 1e3),
            "phiPn_c_kN": float(caps["phiPn_c"][k] / 1e3),
            "utilisation": float(res.util_env[k]),
            "status": "OK" if res.util_env[k] <= 1.0 else "OVER",
        })
    return rows


def node_table(res, case=None):
    """Nodal displacement rows for one load case (default: serviceability)."""
    m = res.model
    g = m.geom
    name = case or res.service_case
    r = res.cases[name]
    d = r.displacement_vectors()
    supported = set(_support_nodes(m).tolist())
    rows = []
    for k in range(g.n_nodes):
        rows.append({
            "node": k,
            "x_m": float(g.nodes[k, 0]),
            "y_m": float(g.nodes[k, 1]),
            "z_m": float(g.nodes[k, 2]),
            "ring": int(g.ring_index[k]),
            "ux_mm": float(d[k, 0] * 1e3),
            "uy_mm": float(d[k, 1] * 1e3),
            "uz_mm": float(d[k, 2] * 1e3),
            "|u|_mm": float(np.linalg.norm(d[k]) * 1e3),
            "support": "yes" if k in supported else "",
        })
    return rows


def reaction_table(res, case=None):
    """Support reaction rows for one load case (default: governing)."""
    m = res.model
    g = m.geom
    name = case or res.governing_case or res.order[0]
    r = res.cases[name]
    R = r.reaction_vectors()
    rows = []
    for k in _support_nodes(m):
        rx, ry, rz = R[k]
        rows.append({
            "node": int(k),
            "x_m": float(g.nodes[k, 0]),
            "y_m": float(g.nodes[k, 1]),
            "Rx_kN": float(rx / 1e3),
            "Ry_kN": float(ry / 1e3),
            "Rz_kN": float(rz / 1e3),
            "thrust_kN": float(np.hypot(rx, ry) / 1e3),
        })
    tot = R.sum(axis=0)
    rows.append({
        "node": "TOTAL", "x_m": "", "y_m": "",
        "Rx_kN": float(tot[0] / 1e3), "Ry_kN": float(tot[1] / 1e3),
        "Rz_kN": float(tot[2] / 1e3), "thrust_kN": "",
    })
    return rows
