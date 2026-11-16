"""Nodal load vector assembly and LRFD load combination.

Every routine returns a flat global load vector of length ``3 * n_nodes``
ordered ``[Fx0, Fy0, Fz0, Fx1, ...]``, in newtons, so the vectors can simply
be scaled and added when combinations are formed.

Area loads are lumped at the nodes using the tributary areas computed in
``geometry`` (one third of every adjacent triangular panel):

    * dead cladding  -> applied on the **faceted surface** area
    * snow and live  -> applied on the **horizontally projected** area
    * wind           -> pressure normal to the surface, on surface area

Explicit nodal loads (``LoadSpec.nodal_loads``: node, Fx, Fy, Fz, case) are
added to the vector of their case (D, S, L or W).

The self-straining load T (uniform temperature change + support settlement,
ASCE 7 notation) is *not* a nodal load vector: it enters the solver as member
initial strains alpha*dT and prescribed support displacements
(``solver.solve_model(..., T_factor)``).  ``combine`` therefore skips the T
factor; ``postprocess.run_analysis`` passes it to the solver instead.
"""

from __future__ import annotations

import numpy as np

from .units import G_ACCEL


# ------------------------------------------------------------ dead loads ----
def self_weight(model):
    """Member self weight, lumped half to each end node."""
    g = model.geom
    F = np.zeros(3 * g.n_nodes)
    rho = model.material.rho
    L = g.member_lengths()
    A = model.member_areas()
    w = A * L * rho * G_ACCEL                     # total weight per member [N]
    for m, (i, j) in enumerate(g.members):
        F[3 * i + 2] -= 0.5 * w[m]
        F[3 * j + 2] -= 0.5 * w[m]
    return F


def cladding_dead(model):
    """Superimposed dead pressure acting on the faceted surface, vertically."""
    g = model.geom
    F = np.zeros(3 * g.n_nodes)
    p = model.loads.dead_pressure
    F[2::3] -= p * g.trib_area
    return F


def dead_load(model):
    """Total dead load vector D."""
    F = cladding_dead(model)
    if model.loads.include_self_weight:
        F = F + self_weight(model)
    return F


# ------------------------------------------------------- snow / live load ---
def _projected_pressure_load(model, pressure, unbalanced=False):
    g = model.geom
    F = np.zeros(3 * g.n_nodes)
    if abs(pressure) < 1e-12:
        return F
    factor = np.ones(g.n_nodes)
    if unbalanced:
        # ASCE-style unbalanced roof load: 1.5 p on the leeward half, 0.5 p on
        # the windward half, split by the plane x = 0.
        factor = np.where(g.nodes[:, 0] >= 0.0, 1.5, 0.5)
    F[2::3] -= pressure * g.trib_area_proj * factor
    return F


def snow_load(model):
    return _projected_pressure_load(model, model.loads.snow_pressure,
                                    model.loads.snow_unbalanced)


def live_load(model):
    F = _projected_pressure_load(model, model.loads.live_pressure)
    if abs(model.loads.crown_point_load) > 1e-12:
        F[3 * model.geom.crown_node + 2] -= model.loads.crown_point_load
    return F


# -------------------------------------------------------------- wind load ---
# Simplified external pressure coefficient for a domed roof, tabulated along
# the meridian measured from the windward springing line.  theta = 0 deg is the
# windward base, 90 deg the crown, 180 deg the leeward base.  The shape follows
# the trend of ASCE 7 Figure 27.3-2 for a hemispherical dome and is documented
# in the report as a deliberate simplification: positive pressure over the
# lower windward face, strong suction over the crown, moderate suction leeward.
_CP_THETA = np.array([0.0, 30.0, 60.0, 90.0, 120.0, 150.0, 180.0])
_CP_VALUE = np.array([0.80, 0.20, -0.60, -1.00, -0.80, -0.50, -0.40])


def wind_pressure_coefficients(model):
    """External pressure coefficient Cp at every node (outward positive = suction)."""
    g = model.geom
    az = np.deg2rad(model.loads.wind_azimuth)
    # Unit vector the wind blows *towards*.
    d = np.array([np.cos(az), np.sin(az), 0.0])
    n = g.normals
    cos_t = np.clip(-np.dot(n, d), -1.0, 1.0)
    theta = np.degrees(np.arccos(cos_t))
    return np.interp(theta, _CP_THETA, _CP_VALUE), theta


def velocity_pressure(model):
    """ASCE 7 velocity pressure qz = 0.613 Kz Kzt Kd V^2  [Pa], V in m/s."""
    ls = model.loads
    return 0.613 * ls.wind_Kz * ls.wind_Kzt * ls.wind_Kd * ls.wind_speed ** 2


def wind_load(model):
    """Wind pressure applied normal to the dome surface at every node."""
    g = model.geom
    F = np.zeros(3 * g.n_nodes)
    if not model.loads.wind_enabled:
        return F
    q = velocity_pressure(model) * model.loads.wind_G
    Cp, _ = wind_pressure_coefficients(model)
    p = q * Cp                                        # [Pa]
    # Cp > 0 pushes inward (along -n); Cp < 0 (suction) pulls outward (+n).
    vec = (-p * g.trib_area)[:, None] * g.normals
    F[0::3] = vec[:, 0]
    F[1::3] = vec[:, 1]
    F[2::3] = vec[:, 2]
    return F


# ------------------------------------------------------------ combination ---
CASE_KEYS = ("D", "S", "L", "W")


def nodal_load_vectors(model):
    """Explicit nodal loads, one vector per case: {'D','S','L','W'} -> (3n,).

    Every entry of ``model.loads.nodal_loads`` is a dict
    ``{"node": int, "Fx": N, "Fy": N, "Fz": N, "case": "D"|"S"|"L"|"W"}``
    (missing force components are 0, a missing case is "L").
    """
    n = model.geom.n_nodes
    out = {k: np.zeros(3 * n) for k in CASE_KEYS}
    for k, p in enumerate(getattr(model.loads, "nodal_loads", None) or []):
        nd = int(p["node"])
        if nd < 0 or nd >= n:
            raise ValueError(f"Nodal load {k}: node {nd} is outside "
                             f"0..{n - 1}.")
        case = str(p.get("case", "L") or "L").strip().upper()
        if case not in out:
            raise ValueError(f"Nodal load {k}: case {case!r} is not one of "
                             f"{', '.join(CASE_KEYS)} (temperature and "
                             f"settlement belong to T, not to nodal loads).")
        for a, key in enumerate(("Fx", "Fy", "Fz")):
            v = float(p.get(key, 0.0) or 0.0)
            if not np.isfinite(v):
                raise ValueError(f"Nodal load {k}: {key} is not finite.")
            out[case][3 * nd + a] += v
    return out


def basic_load_vectors(model):
    """Return the dictionary of unfactored load vectors {'D','S','L','W'}.

    Explicit nodal loads are included in their case.  The self-straining
    load T is not a load vector and is not in this dictionary.
    """
    basic = {
        "D": dead_load(model),
        "S": snow_load(model),
        "L": live_load(model),
        "W": wind_load(model),
    }
    if getattr(model.loads, "nodal_loads", None):
        extra = nodal_load_vectors(model)
        for k in CASE_KEYS:
            if np.any(extra[k] != 0.0):
                basic[k] = basic[k] + extra[k]
    return basic


def combine(basic, combination):
    """Factored load vector for one :class:`LoadCombination`.

    Keys without a load vector in ``basic`` (the self-straining factor T)
    are skipped; T is applied by the solver (``solve_model(T_factor=...)``).
    """
    F = np.zeros_like(basic["D"])
    for key, factor in combination.factors().items():
        if key not in basic:
            continue
        if abs(factor) > 1e-12:
            F = F + factor * basic[key]
    return F


def load_totals(F):
    """Resultant of a global load vector: (sum Fx, sum Fy, sum Fz) in N."""
    return float(F[0::3].sum()), float(F[1::3].sum()), float(F[2::3].sum())


def describe_loads(model):
    """Human-readable breakdown of the unfactored load cases."""
    basic = basic_load_vectors(model)
    lines = []
    names = {"D": "Dead (self weight + cladding)",
             "S": "Snow (projected)",
             "L": "Live / point",
             "W": "Wind"}
    for k in ("D", "S", "L", "W"):
        fx, fy, fz = load_totals(basic[k])
        lines.append(f"{names[k]:<34s} Fx = {fx / 1e3:9.2f} kN   "
                     f"Fy = {fy / 1e3:9.2f} kN   Fz = {fz / 1e3:9.2f} kN")
    if model.loads.wind_enabled:
        q = velocity_pressure(model)
        lines.append(f"{'Velocity pressure qz':<34s} {q / 1e3:.3f} kPa "
                     f"(V = {model.loads.wind_speed:.1f} m/s)")
    nl = getattr(model.loads, "nodal_loads", None) or []
    if nl:
        lines.append(f"{'Explicit nodal loads':<34s} {len(nl)} entr"
                     f"{'y' if len(nl) == 1 else 'ies'} (included above)")
    lines += describe_self_straining(model)
    return lines


def describe_self_straining(model):
    """Text lines describing the self-straining load T."""
    ls = model.loads
    dT = float(getattr(ls, "temperature_change", 0.0) or 0.0)
    sets = getattr(ls, "settlements", None) or []
    name = "Self-straining T"
    if dT == 0.0 and not any(
            float(s.get(k, 0.0) or 0.0) != 0.0
            for s in sets for k in ("ux", "uy", "uz")):
        return [f"{name:<34s} none (no temperature change or settlement)"]
    lines = []
    if dT != 0.0:
        a = model.material.alpha
        lines.append(f"{name + ': temperature':<34s} dT = {dT:+.1f} C on every "
                     f"strut, alpha = {a * 1e6:.1f}e-6 /C, free strain "
                     f"alpha*dT = {a * dT:.3e}")
    if sets:
        mx = max(max(abs(float(s.get(k, 0.0) or 0.0))
                     for k in ("ux", "uy", "uz")) for s in sets)
        lines.append(f"{name + ': settlement':<34s} {len(sets)} support "
                     f"node(s), largest movement {mx * 1e3:.2f} mm")
    return lines
