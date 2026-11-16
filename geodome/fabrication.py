"""
Fabrication schedule: cut list, hub schedule, panel schedule, dihedral angles.

The analysis treats every strut as a line between two hub centres.  A
builder needs more than that: how long to cut each tube, at what angle it
meets its hub, how many different hubs and panels there are, where each one
goes, and at what angle two neighbouring panels meet.  This module derives
all of it from the node coordinates of a :class:`geometry.DomeGeometry`
(generated geodesic dome or custom truss) and, when a :class:`DomeModel` is
given, from its sections and material (tube masses).

Every angle is computed from the vectors of the actual mesh and, where a
closed form exists, checked against it in the same table, so each schedule
carries its own proof of correctness.

Notation
--------
Hub (node) positions ``x_i``; strut ``k`` from hub ``i`` to hub ``j`` with
nominal length ``L = |x_j - x_i|`` (hub centre to hub centre).  Sphere centre
``c`` and radius ``R``: for a generated dome the generating sphere, for a
custom truss the least-squares sphere of :func:`geometry.best_fit_sphere`
(``R`` = NaN when the nodes are not spherical).

Struts (cut list)
-----------------
Chord factor (Clinton 1971, Sec. 1.4; Kenner 1976)::

    cf = L / R

Cut length and tube mass for a hub offset ``o`` (distance from the hub
centre to the end of the tube, equal at both ends)::

    L_cut = L - 2 o ,          m = rho * A * L_cut

One schedule row per strut family; a family whose individual lengths differ
by more than the cut tolerance (default 0.5 mm; the analysis families use
1e-4 R, which at 8V-10V on a large sphere can lump lengths 1-2 mm apart) is
split into cut groups U.1, U.2, ... so the cut list is exact.

Axial angle -- "an angle formed by an element and a radius from the center
of the polyhedron meeting in a common point" (Clinton 1971, Sec. 1.4).  It is
computed from vectors at *both* ends of every strut::

    alpha_i = angle(x_j - x_i, -u_i) ,   u_i = (x_i - c) / |x_i - c|
            = atan2(|d x (-u_i)|, d . (-u_i)),   d = x_j - x_i

and checked against the closed form for hubs on the sphere::

    alpha = arccos(L / (2 R))

(with c = 0 and |x_i| = |x_j| = R:  d . (-x_i) = R^2 - x_i . x_j = L^2 / 2,
so cos alpha = (L^2/2) / (L R) = L / (2R); both ends give the same angle).
The tangent ("bend") angle is the complement (Domerama glossary: "tangential
angle ... namely 90 minus the axial angle")::

    tau = 90 deg - alpha = arcsin(L / (2R))

i.e. the angle between the strut and the hub's tangent plane, half the
central angle subtended by the strut.  It is the bend of a flattened-end
pipe strut and the socket tilt of a hub.  When the nodes are not spherical
the hub axis ``u_i`` is the area-weighted outward panel normal (with no
panels: the unit vector pointing away from the hub's struts) and the
closed-form columns are NaN.

Panels (cladding)
-----------------
::

    area     A = |(x_b - x_a) x (x_c - x_a)| / 2
    corner   cos gamma_a = (b^2 + c^2 - a^2) / (2 b c)    (opposite side a)

Identity (tested): the sum of all panel areas equals the faceted surface
area, i.e. the sum of the nodal tributary areas the pressure loads use.

Dihedral angles (interior angle between two panels across a strut, measured
inside the dome).  From the consistently oriented outward unit normals::

    theta = atan2(|n1 x n2|, n1 . n2)          (angle between the normals)
    delta = 180 deg - theta    convex fold (third corner of panel 2 lies
                               on the inner side of panel 1)
    delta = 180 deg + theta    re-entrant (valley) fold

``bend = 180 deg - delta`` is the angle a panel bracket must turn.
Independent closed form for panels inscribed in the sphere, which uses only
the side lengths and R::

    delta = sum_{f = 1, 2} atan2( sqrt(R^2 - rho_f^2),  rho_f cos gamma_f )

where ``rho_f = a b c / (4 A_f)`` is the panel circumradius and ``gamma_f``
the panel corner opposite the shared strut.  Derivation: the foot of the
perpendicular from c onto the panel plane is the panel circumcentre O_f
(c is equidistant from the three corners), at distance sqrt(R^2 - rho_f^2)
from c and at signed distance rho_f cos gamma_f from the strut midpoint M
(the distance from a circumcentre to a chord); in the plane through M normal
to the strut, the panel trace and the line M-c therefore enclose the angle
atan2(sqrt(R^2 - rho_f^2), rho_f cos gamma_f), and the two panels lie on
either side of M-c.  For acute panels this is the F8 hand-calculation form
asin(d_1/h) + asin(d_2/h), h = sqrt(R^2 - L^2/4).

Hubs
----
Valence = number of struts at the hub.  Hub type = (location, multiset of
strut families, cyclic order of the families counter-clockwise seen from
outside, canonical over rotation and mirror).  Location: ``base`` (a
supported node), ``edge`` (on a free edge of the panel mesh, not supported)
or ``interior``.  The in-plane gaps are the angles between consecutive struts
projected on the hub's tangent plane (what a hub plate or a printed hub needs
together with the tangent angle of each socket).  For base and edge hubs the
order starts after the largest gap, so the open side is the last gap.

Counts of a complete Class I sphere (T = nu^2; tested): nodes 10 nu^2 + 2,
struts 30 nu^2, panels 20 nu^2; twelve hubs of valence 5, all others 6;
V - E + F = 2 (Euler).

Stock (cutting) plan
--------------------
Optional.  A bar of stock length ``S`` cut with kerf ``k`` yields
``n = floor((S + k) / (L_cut + k))`` pieces of one family (n pieces need
n L_cut + (n - 1) k <= S), so a family of N struts needs ceil(N / n) bars.
Mixing families of the same section on one bar can save bars; that bin
packing problem is solved with the first-fit-decreasing heuristic (Johnson
1973) and compared with the lower bound
max(ceil(sum(L_cut + k) / (S + k)), number of pieces with L_cut + k > (S + k)/2)
(no two pieces longer than half a bar can share one).
The better of the two plans (single-family or mixed) is reported as cutting
patterns.

Published values reproduced (tests/test_fabrication.py)
-------------------------------------------------------
2V chord factors 0.54653 / 0.61803 (Domerama chord-factor chart; closed
forms 2 sin(15.8587 deg) and 2 sin 18 deg), 2V tangent angles 15.86 / 18 deg (Domerama), the 3V to 8V dome
chord factors and strut counts of the Domerama chart ("method 1",
icosahedral; 5 printed decimals to within one unit, 8 printed decimals to
within 1e-8), the full-sphere counts, the regular icosahedron and
octahedron dihedral angles, and every value of the F8 hand calculation
(validation/hand_calc_ext.py) to 1e-15 relative.

References
----------
* J. D. Clinton (1971). Advanced Structural Geometry Studies, Part I:
  Polyhedral Subdivision Concepts for Structural Applications. NASA CR-1734,
  Sec. 1.4 (chord factor, axial angle).
* H. Kenner (1976). Geodesic Math and How to Use It. University of
  California Press (2V chord factors; chord-factor tables).
* Domerama chord-factor charts, https://domerama.com/calculators/chord-charts
  (2V-6V "method 1" icosahedral chord factors and counts).
* D. S. Johnson (1973). Near-optimal bin packing algorithms. PhD thesis, MIT
  (first-fit decreasing).

Units: SI internally (m, kg); the tables show lengths in mm, angles in
degrees, areas in m^2 and masses in kg, with the unit in every column name.
"""

from __future__ import annotations

import csv
import math
import os

import numpy as np

from . import geometry as geo

__all__ = [
    "TABLE_ORDER", "TABLE_TITLES", "STRUT_COLOURS",
    "fabrication_tables", "export_fabrication_csv", "fabrication_text",
    "sphere_of", "family_labels", "stock_plan",
]

DEG = 180.0 / math.pi

# Display order of the tables (the GUI and the CSV export use it).
TABLE_ORDER = ("summary", "struts", "hubs", "panels", "dihedrals",
               "members", "nodes", "faces", "stock", "stock_patterns",
               "notes")

TABLE_TITLES = {
    "summary": "Fabrication summary",
    "struts": "Strut schedule (cut list)",
    "hubs": "Hub schedule",
    "panels": "Panel schedule",
    "dihedrals": "Dihedral angles between panels",
    "members": "Strut assembly list",
    "nodes": "Hub assembly list",
    "faces": "Panel assembly list",
    "stock": "Stock bars per strut family",
    "stock_patterns": "Cutting patterns",
    "notes": "Notes",
}

# Colour codes for strut families, the usual practice in dome kits.
STRUT_COLOURS = ("red", "blue", "green", "yellow", "orange", "purple",
                 "white", "black", "brown", "grey", "pink", "cyan")

_LOCATION_ORDER = {"interior": 0, "edge": 1, "base": 2}


# ================================================================ helpers ===
def family_labels(n_families):
    """Family labels A, B, ..., Z, S26, S27, ... (as geometry.strut_schedule)."""
    return [chr(ord("A") + g) if g < 26 else f"S{g}" for g in range(n_families)]


def _resolve(model):
    """Return ``(model_or_None, geometry)`` for a DomeModel or a DomeGeometry."""
    if isinstance(model, geo.DomeGeometry):
        return None, model
    g = getattr(model, "geom", None)
    if isinstance(g, geo.DomeGeometry):
        return model, g
    raise TypeError("fabrication_tables needs a DomeModel or a DomeGeometry, "
                    f"got {type(model).__name__}.")


def _angle_rows(a, b):
    """Row-wise angle [rad] between vectors a and b (atan2, accurate near 0)."""
    cr = np.linalg.norm(np.cross(a, b), axis=1)
    dt = np.einsum("ij,ij->i", a, b)
    return np.arctan2(cr, dt)


def _heron(a, b, c):
    """Triangle area from its sides (Kahan's stable form of Heron's rule)."""
    a, b, c = sorted((float(a), float(b), float(c)), reverse=True)
    q = (a + (b + c)) * (c - (a - b)) * (c + (a - b)) * (a + (b - c))
    return 0.25 * math.sqrt(max(q, 0.0))


def sphere_of(geom):
    """Sphere through the hubs: ``(centre (3,), R, max_dev)``.

    Generated geodesic dome: the generating sphere, recovered exactly from
    the radial node normals (``c = x_i - R n_i`` for every node).  Custom
    truss: the least-squares sphere (:func:`geometry.best_fit_sphere`).
    ``max_dev = max_i | |x_i - c| - R |`` (0 to round-off on a geodesic
    dome).  All NaN when the nodes do not define a sphere.
    """
    nan3 = np.full(3, np.nan)
    P = np.asarray(geom.nodes, dtype=float)
    R = float(geom.radius)
    if not (np.isfinite(R) and R > 0.0) or len(P) < 2:
        return nan3, float("nan"), float("nan")
    centre = None
    if getattr(geom, "kind", "geodesic") == "geodesic":
        N = np.asarray(geom.normals, dtype=float)
        if N.shape == P.shape and np.allclose(np.linalg.norm(N, axis=1), 1.0,
                                              atol=1e-9):
            centre = np.mean(P - R * N, axis=0)
    if centre is None:
        c, Rfit, _ = geo.best_fit_sphere(P)
        if not np.isfinite(Rfit):
            return nan3, float("nan"), float("nan")
        centre, R = np.asarray(c, dtype=float), float(Rfit)
    dev = float(np.max(np.abs(np.linalg.norm(P - centre, axis=1) - R)))
    return centre, R, dev


def _hub_axes(P, M, normals, centre):
    """Outward unit hub axis at every node (see the module docstring)."""
    n = len(P)
    if np.all(np.isfinite(centre)):
        u = P - centre
        nr = np.linalg.norm(u, axis=1)
        ok = nr > 0.0
        u[ok] /= nr[ok][:, None]
        if np.all(ok):
            return u
        bad = ~ok
    else:
        u = np.zeros((n, 3))
        bad = np.ones(n, dtype=bool)
    N = np.asarray(normals, dtype=float)
    if N.shape == P.shape:
        nn = np.linalg.norm(N, axis=1)
        use = bad & (nn > 1e-12)
        u[use] = N[use] / nn[use][:, None]
        bad &= ~use
    if np.any(bad):
        acc = np.zeros((n, 3))
        d = P[M[:, 1]] - P[M[:, 0]]
        d /= np.linalg.norm(d, axis=1)[:, None]
        np.add.at(acc, M[:, 0], -d)                  # away from the struts
        np.add.at(acc, M[:, 1], d)
        for i in np.where(bad)[0]:
            a = acc[i]
            na = np.linalg.norm(a)
            u[i] = a / na if na > 1e-9 else np.array([0.0, 0.0, 1.0])
    return u


def _orient_faces(P, F, ref):
    """Orient every panel outwards, consistently across shared struts.

    Panels are first made consistent with their neighbours (two panels that
    share a strut must run along it in opposite directions), then each
    connected patch is turned so that  sum_f (x_b - x_a) x (x_c - x_a) .
    (centroid_f - ref) > 0  -- for a closed surface that sum is 6 x volume
    whatever ``ref`` is; for an open cap ``ref`` is the sphere centre (or
    the node centroid).  A patch for which the sum vanishes is turned to
    face +z.  Returns ``(oriented_faces, unit_normals, areas, notes)``.
    """
    notes = []
    nf = len(F)
    if nf == 0:
        return np.zeros((0, 3), dtype=int), np.zeros((0, 3)), np.zeros(0), \
            notes
    F = [tuple(int(v) for v in t) for t in F]

    edge_faces = {}
    for f, (a, b, c) in enumerate(F):
        for i, j in ((a, b), (b, c), (c, a)):
            edge_faces.setdefault((min(i, j), max(i, j)), []).append(f)

    def directed(tri):
        a, b, c = tri
        return ((a, b), (b, c), (c, a))

    flip = np.full(nf, -1, dtype=int)
    comp = np.full(nf, -1, dtype=int)
    orientable = True
    ncomp = 0
    for s in range(nf):
        if comp[s] >= 0:
            continue
        comp[s], flip[s] = ncomp, 0
        stack = [s]
        while stack:
            f = stack.pop()
            tri = F[f][::-1] if flip[f] else F[f]
            for i, j in directed(tri):
                nb = edge_faces[(min(i, j), max(i, j))]
                if len(nb) != 2:
                    continue
                g = nb[0] if nb[1] == f else nb[1]
                need = 1 if (i, j) in directed(F[g]) else 0
                if comp[g] < 0:
                    comp[g], flip[g] = ncomp, need
                    stack.append(g)
                elif flip[g] != need:
                    orientable = False
        ncomp += 1

    out = np.array([F[f][::-1] if flip[f] else F[f] for f in range(nf)],
                   dtype=int)
    A, B, C = P[out[:, 0]], P[out[:, 1]], P[out[:, 2]]
    cross = np.cross(B - A, C - A)
    cen = (A + B + C) / 3.0
    scale = float(np.ptp(P, axis=0).max()) or 1.0
    mag = np.linalg.norm(cross, axis=1)

    if not orientable:
        notes.append("The panel mesh is not orientable (a Moebius-type "
                     "patch): each panel was oriented on its own, away from "
                     "the reference point.")
        s = np.einsum("ij,ij->i", cross, cen - ref)
        tol = 1e-9 * mag * scale
        turn = (s < -tol) | ((np.abs(s) <= tol) & (cross[:, 2] < 0.0))
    else:
        turn = np.zeros(nf, dtype=bool)
        for k in range(ncomp):
            sel = comp == k
            s = float(np.einsum("ij,ij->", cross[sel], cen[sel] - ref))
            tol = 1e-9 * float(mag[sel].sum()) * scale
            if abs(s) > tol:
                t = s < 0.0
            else:
                t = float(cross[sel, 2].sum()) < 0.0
            turn[sel] = t
    out[turn] = out[turn][:, ::-1]
    cross[turn] = -cross[turn]
    ok = mag > 0.0
    unit = np.zeros_like(cross)
    unit[ok] = cross[ok] / mag[ok][:, None]
    if not np.all(ok):
        notes.append(f"{int(np.sum(~ok))} panel(s) have zero area "
                     f"(collinear corners).")
    return out, unit, 0.5 * mag, notes


def _canonical_cycle(fids, az, open_fan=False):
    """Canonical (rotation and mirror) cyclic order of the struts at a hub.

    ``fids`` and ``az`` are family ids and azimuths [rad] sorted
    counter-clockwise.  Returns ``(fid_tuple, gaps_rad_tuple)`` of the
    lexicographically smallest of the 2 v rotations/reflections; gap k is
    the angle between struts k and k+1 (the last one closes the circle).
    With ``open_fan`` (base and edge hubs) only the rotations that start
    right after the largest gap are allowed, so the struts read from one
    side of the opening to the other and the open side is the last gap.
    """
    v = len(fids)
    if v == 0:
        return (), ()
    best = None
    for seq, a in ((list(fids), list(az)),
                   (list(fids)[::-1], [-x for x in list(az)[::-1]])):
        gaps = [((a[(k + 1) % v] - a[k]) % (2.0 * math.pi)) if v > 1
                else 2.0 * math.pi for k in range(v)]
        gmax = max(gaps)
        for r in range(v):
            if open_fan and gaps[r - 1] < gmax - 1e-9:
                continue
            s = tuple(seq[r:] + seq[:r])
            g = tuple(gaps[r:] + gaps[:r])
            key = (s, tuple(round(x, 9) for x in g))
            if best is None or key < best[0]:
                best = (key, s, g)
    return best[1], best[2]


def _hub_cycles(P, M, fam, axes):
    """Per node: family ids and azimuths of its struts, counter-clockwise."""
    n = len(P)
    nbr = [[] for _ in range(n)]
    for k, (i, j) in enumerate(M):
        nbr[i].append((int(j), int(fam[k])))
        nbr[j].append((int(i), int(fam[k])))
    out = []
    for i in range(n):
        e3 = axes[i]
        t = np.eye(3)[int(np.argmin(np.abs(e3)))]
        e1 = t - (t @ e3) * e3
        e1 /= np.linalg.norm(e1)
        e2 = np.cross(e3, e1)
        items = []
        for j, fid in nbr[i]:
            d = P[j] - P[i]
            items.append((math.atan2(float(d @ e2), float(d @ e1)), fid))
        items.sort()
        out.append(([f for _, f in items], [a for a, _ in items]))
    return out


def _multiset_label(fids, labels):
    if len(fids) == 0:
        return "(none)"
    ids, cnt = np.unique(np.asarray(fids, dtype=int), return_counts=True)
    return "+".join(f"{c}{labels[i] if i >= 0 else '?'}"
                    for i, c in zip(ids, cnt))


def _seq_label(fids, labels, sep="-"):
    return sep.join(labels[i] if i >= 0 else "?" for i in fids)


def _panel_label(fids, labels):
    names = [labels[i] if i >= 0 else "?" for i in fids]
    return ("" if all(len(s) == 1 for s in names) else "-").join(names)


# ========================================================== stock plan ======
def stock_plan(pieces, stock_length, kerf=0.0):
    """Cut pieces from stock bars.

    ``pieces``: list of ``(family_label, section, cut_length_m, count)``.
    Returns ``(family_rows, pattern_rows, section_rows, notes)``:

    * single-family plan: ``n = floor((S + k)/(L_cut + k))`` pieces per bar,
      ``ceil(N / n)`` bars per family;
    * mixed plan per section: first-fit decreasing, each piece occupying
      ``L_cut + k`` of a bar of capacity ``S + k``;
    * lower bound per section ``ceil(sum(L_cut + k) / (S + k))``.

    The better plan of each section (fewer bars; single-family on a tie)
    is listed as cutting patterns.
    """
    S = float(stock_length)
    k = float(kerf)
    if not (np.isfinite(S) and S > 0.0):
        raise ValueError(f"Stock length must be a positive length, got {S!r}.")
    if not (np.isfinite(k) and k >= 0.0):
        raise ValueError(f"Saw kerf must be >= 0, got {k!r}.")
    tol = 1e-9
    notes = []
    fam_rows = []
    sections = {}
    for label, section, c, N in pieces:
        c, N = float(c), int(N)
        n = int(math.floor((S + k) / (c + k) + tol)) if c > 0 else 0
        if c > S + tol:
            n = 0
        if n == 0:
            notes.append(f"Strut family {label}: cut length {c * 1e3:.1f} mm "
                         f"is longer than the {S * 1e3:.0f} mm stock; order "
                         f"it as a special length (not included in the bars).")
            fam_rows.append({
                "family": label, "section": section,
                "cut_length_mm": c * 1e3, "count": N,
                "pieces_per_bar": 0, "bars": float("nan"),
                "offcut_per_bar_mm": float("nan"),
                "waste_pct": float("nan")})
            continue
        bars = int(math.ceil(N / n))
        off = S - (n * c + (n - 1) * k)
        waste = 100.0 * (bars * S - N * c) / (bars * S)
        fam_rows.append({
            "family": label, "section": section, "cut_length_mm": c * 1e3,
            "count": N, "pieces_per_bar": n, "bars": bars,
            "offcut_per_bar_mm": off * 1e3, "waste_pct": waste})
        sections.setdefault(section, []).append((label, c, N, n, bars))

    pat_rows, sec_rows = [], []
    for section, fams in sections.items():
        # ---- single-family plan
        single_bars = sum(b for *_, b in fams)
        single_patterns = []
        for label, c, N, n, bars in fams:
            full, rest = divmod(N, n)
            if full:
                single_patterns.append(((((label, n),), c * n + k * (n - 1)),
                                        full))
            if rest:
                single_patterns.append(((((label, rest),),
                                         c * rest + k * (rest - 1)), 1))
        # ---- mixed plan: first-fit decreasing
        items = []
        for label, c, N, n, bars in fams:
            items += [(c, label)] * N
        items.sort(key=lambda t: (-t[0], t[1]))
        cap = S + k
        bins = []                                   # [remaining, {label: n}]
        for c, label in items:
            need = c + k
            for b in bins:
                if b[0] >= need - tol:
                    b[0] -= need
                    b[1][label] = b[1].get(label, 0) + 1
                    break
            else:
                bins.append([cap - need, {label: 1}])
        mixed_bars = len(bins)
        lens = {label: c for label, c, *_ in fams}
        agg = {}
        for rem, cnt in bins:
            key = tuple(sorted(cnt.items()))
            used = sum(lens[lb] * q for lb, q in key) \
                + k * (sum(q for _, q in key) - 1)
            agg[key] = (used, agg.get(key, (0.0, 0))[1] + 1)
        mixed_patterns = [((key, used), cnt) for key, (used, cnt)
                          in agg.items()]
        total_need = sum((c + k) * N for _, c, N, *_ in fams)
        # Lower bound: total length, and pieces longer than half a bar
        # (no two of them can share one bar).
        n_long = sum(N for _, c, N, *_ in fams if c + k > 0.5 * cap + tol)
        lower = max(int(math.ceil(total_need / cap - tol)), n_long)
        use_mixed = mixed_bars < single_bars
        chosen = mixed_patterns if use_mixed else single_patterns
        n_pieces = sum(N for _, _, N, *_ in fams)
        cut_total = sum(c * N for _, c, N, *_ in fams)
        best = mixed_bars if use_mixed else single_bars
        sec_rows.append({
            "section": section, "pieces": n_pieces,
            "bars_single_family": single_bars, "bars_mixed": mixed_bars,
            "lower_bound": lower,
            "plan": "mixed (first-fit decreasing)" if use_mixed
            else "single family per bar",
            "bars": best,
            "waste_pct": 100.0 * (best * S - cut_total) / (best * S)})
        chosen.sort(key=lambda t: (-t[1], t[0][0]))
        for (key, used), cnt in chosen:
            pat_rows.append({
                "section": section,
                "pattern": " + ".join(f"{q} x {lb}" for lb, q in key),
                "bars": cnt,
                "used_mm": used * 1e3,
                "offcut_mm": (S - used) * 1e3})
    return fam_rows, pat_rows, sec_rows, notes


def _strut_row(label, index, sel, L, o, ax_i, ax_j, ax_f, R, spherical,
               section, section_name, rho, length_tol, notes):
    """One row of the strut schedule for the struts ``sel`` (one cut group)."""
    Lg = L[sel]
    Lm = float(Lg.mean())
    cut = Lm - 2.0 * o
    alphas = np.concatenate([ax_i[sel], ax_j[sel]])
    a_mean = float(alphas.mean())
    if spherical:
        af = float(np.arccos(np.clip(Lm / (2.0 * R), -1.0, 1.0)))
        ref = np.concatenate([ax_f[sel], ax_f[sel]])
        chk = float(np.max(np.abs(alphas - ref)))
    else:
        af = chk = float("nan")
    mass_each = rho * section.A * cut if section is not None else float("nan")
    spread = float(np.ptp(Lg))
    if spread > length_tol:
        notes.append(f"Cut group {label}: individual lengths still differ by "
                     f"{spread * 1e3:.2f} mm (a chain of close lengths); cut "
                     f"these struts to their own lengths (members table).")
    return {
        "family": label,
        "colour": STRUT_COLOURS[index] if index < len(STRUT_COLOURS) else "",
        "count": int(len(sel)),
        "section": section_name,
        "length_mm": Lm * 1e3,
        "cut_length_mm": cut * 1e3,
        "chord_factor": Lm / R if spherical else float("nan"),
        "axial_angle_deg": a_mean * DEG,
        "tangent_angle_deg": 90.0 - a_mean * DEG,
        "mass_each_kg": mass_each,
        "mass_total_kg": mass_each * len(sel),
        "total_cut_length_m": cut * len(sel),
        "length_spread_mm": spread * 1e3,
        "axial_spread_deg": float(np.ptp(alphas) * DEG),
        "axial_formula_deg": af * DEG,
        "axial_check_deg": chk * DEG,
    }


# ===================================================== main entry point ====
def fabrication_tables(model, hub_offset=0.0, stock_length=None, kerf=0.0,
                       length_tol=5e-4):
    """Builder-facing schedules for a dome model (or a bare geometry).

    Parameters
    ----------
    model : DomeModel or DomeGeometry
        With a DomeModel the strut tables include the section and the tube
        mass, and ``base`` hubs are the model's supported nodes; with a bare
        geometry they are ``geom.base_nodes`` and the masses are NaN.
    hub_offset : float
        Distance [m] from a hub centre to the end of a tube, the same at both
        ends (hub radius plus connector).  ``cut length = L - 2 hub_offset``.
    stock_length, kerf : float, optional
        Stock bar length and saw kerf [m]; when ``stock_length`` is given the
        result also holds a stock (cutting) plan.
    length_tol : float
        Cut tolerance [m], default 0.5 mm.  The analysis groups struts into
        families with the geometry tolerance (1e-4 R on a geodesic dome, so
        2 mm on a 20 m sphere); a family whose sorted lengths have a gap
        larger than ``length_tol`` is split into cut groups ``U.1, U.2, ...``
        in the strut schedule (and the stock plan), with a note.  Hubs,
        panels, dihedrals and sections keep the analysis family; the members
        table gives both (``family`` and ``cut_group``).

    Returns
    -------
    dict
        ``{name: [row dict, ...]}`` for the names in :data:`TABLE_ORDER`
        (``stock`` and ``stock_patterns`` only with a stock length).  Every
        column name carries its unit (``_mm``, ``_m``, ``_m2``, ``_deg``,
        ``_kg``, ``_pct``); NaN means "not defined" (for example a chord
        factor of a non-spherical truss).
    """
    mdl, g = _resolve(model)
    P = np.asarray(g.nodes, dtype=float)
    M = np.asarray(g.members, dtype=int).reshape(-1, 2)
    F = np.asarray(g.faces, dtype=int).reshape(-1, 3)
    n, m, nf = len(P), len(M), len(F)
    fam = np.asarray(g.strut_group, dtype=int)
    nfam = int(fam.max()) + 1 if m else 0
    labels = family_labels(nfam)

    o = float(hub_offset)
    if not (np.isfinite(o) and o >= 0.0):
        raise ValueError(f"Hub offset must be a length >= 0, got {hub_offset!r}.")
    d = P[M[:, 1]] - P[M[:, 0]]
    L = np.linalg.norm(d, axis=1)
    if m and 2.0 * o >= float(L.min()):
        raise ValueError(
            f"Hub offset {o * 1e3:.1f} mm leaves no tube: the shortest strut "
            f"is {L.min() * 1e3:.1f} mm between hub centres, so the offset "
            f"must be below {L.min() * 1e3 / 2:.1f} mm.")
    length_tol = float(length_tol)
    if not (np.isfinite(length_tol) and length_tol >= 0.0):
        raise ValueError(f"Cut tolerance must be >= 0, got {length_tol!r}.")

    notes = []
    centre, R, dev = sphere_of(g)
    spherical = bool(np.isfinite(R))
    if not spherical:
        notes.append("The hubs do not lie on a sphere: chord factors and the "
                     "closed-form angle checks are not defined (NaN); axial "
                     "angles use the panel normals as hub axes (without "
                     "panels: the direction away from the hub's struts).")
    ref = centre if spherical else P.mean(axis=0)

    # ------------------------------------------------------------ angles --
    axes = _hub_axes(P, M, g.normals, centre)
    ax_i = _angle_rows(d, -axes[M[:, 0]])
    ax_j = _angle_rows(-d, -axes[M[:, 1]])
    if spherical:
        ax_f = np.arccos(np.clip(L / (2.0 * R), -1.0, 1.0))
    else:
        ax_f = np.full(m, np.nan)

    # ------------------------------------------------------------ panels --
    Fo, nrm, area, fnotes = _orient_faces(P, F, ref)
    notes += fnotes
    edge_member = {(min(int(i), int(j)), max(int(i), int(j))): k
                   for k, (i, j) in enumerate(M)}
    Fo = [tuple(int(v) for v in t) for t in Fo]
    edge_faces = {}
    for f, (a, b, c) in enumerate(Fo):
        for i, j in ((a, b), (b, c), (c, a)):
            edge_faces.setdefault((min(i, j), max(i, j)), []).append(f)

    def edge_fid(i, j):
        k = edge_member.get((min(i, j), max(i, j)))
        return int(fam[k]) if k is not None else -1

    n_loose = sum(1 for e in edge_faces if e not in edge_member)
    if n_loose:
        notes.append(f"{n_loose} panel edge(s) are not struts; their family "
                     f"is shown as '?'.")
    n_nonmanifold = sum(1 for v in edge_faces.values() if len(v) > 2)
    if n_nonmanifold:
        notes.append(f"{n_nonmanifold} strut(s) are shared by more than two "
                     f"panels; no dihedral angle is given for them.")

    # per face: sides sorted by (family, length), angle opposite each side
    face_rec = []
    for f, (a, b, c) in enumerate(Fo):
        sides = []
        for (i, j), opp in (((a, b), c), ((b, c), a), ((c, a), b)):
            fid = edge_fid(int(i), int(j))
            Ls = float(np.linalg.norm(P[j] - P[i]))
            u, v = P[i] - P[opp], P[j] - P[opp]
            ang = math.atan2(float(np.linalg.norm(np.cross(u, v))),
                             float(u @ v))
            sides.append((fid if fid >= 0 else 10 ** 9, Ls, ang, fid))
        cyc = [s[3] for s in sides]                 # ccw seen from outside
        sides.sort(key=lambda s: (s[0], s[1]))
        fids = [s[3] for s in sides]
        face_rec.append({"fids": tuple(fids), "cyc": cyc,
                         "sides": [s[1] for s in sides],
                         "angles": [s[2] for s in sides],
                         "label": _panel_label(fids, labels)})

    # ---- panel types
    ptypes = {}
    for f, rec in enumerate(face_rec):
        ptypes.setdefault(rec["fids"], []).append(f)
    panel_rows = []
    for key in sorted(ptypes, key=lambda t: tuple(x if x >= 0 else 10 ** 9
                                                  for x in t)):
        fl = ptypes[key]
        S = np.array([face_rec[f]["sides"] for f in fl])
        Ang = np.array([face_rec[f]["angles"] for f in fl])
        Ar = area[fl]
        hand = ""
        if len(set(key)) == 3 and min(key) >= 0:
            x, y, z = key
            pos = 0
            for f in fl:
                cyc = face_rec[f]["cyc"]
                r = cyc.index(x)
                if cyc[(r + 1) % 3] == y:
                    pos += 1
            lx, ly, lz = labels[x], labels[y], labels[z]
            hand = (f"{pos} {lx}{ly}{lz} / {len(fl) - pos} {lx}{lz}{ly} "
                    f"(ccw from outside)")
        panel_rows.append({
            "panel_type": face_rec[fl[0]]["label"],
            "count": len(fl),
            "side_1_mm": float(S[:, 0].mean() * 1e3),
            "side_2_mm": float(S[:, 1].mean() * 1e3),
            "side_3_mm": float(S[:, 2].mean() * 1e3),
            "angle_1_deg": float(Ang[:, 0].mean() * DEG),
            "angle_2_deg": float(Ang[:, 1].mean() * DEG),
            "angle_3_deg": float(Ang[:, 2].mean() * DEG),
            "area_m2": float(Ar.mean()),
            "total_area_m2": float(Ar.sum()),
            "handedness": hand,
            "size_spread_mm": float(np.max(np.abs(S - S.mean(axis=0))) * 1e3),
        })
    ptype_of_face = {}
    for f, rec in enumerate(face_rec):
        ptype_of_face[f] = rec["label"]

    # ---- dihedral angles
    edge_dihedral = {}
    edge_dihedral_formula = {}
    for e, fl in edge_faces.items():
        if len(fl) != 2:
            continue
        f1, f2 = fl
        n1, n2 = nrm[f1], nrm[f2]
        theta = math.atan2(float(np.linalg.norm(np.cross(n1, n2))),
                           float(n1 @ n2))
        q = [v for v in Fo[f2] if v not in e][0]
        s = float(n1 @ (P[q] - P[e[0]]))
        delta = math.pi - theta if s <= 0.0 else math.pi + theta
        edge_dihedral[e] = delta
        if spherical:
            # closed form from side lengths only (module docstring)
            Le = float(np.linalg.norm(P[e[1]] - P[e[0]]))
            tot = 0.0
            for f in (f1, f2):
                r = [v for v in Fo[f] if v not in e][0]     # corner opposite
                b_ = float(np.linalg.norm(P[r] - P[e[0]]))
                c_ = float(np.linalg.norm(P[r] - P[e[1]]))
                Af = _heron(Le, b_, c_)
                if Af <= 0.0:
                    tot = float("nan")
                    break
                rho = Le * b_ * c_ / (4.0 * Af)              # circumradius
                cg = (b_ * b_ + c_ * c_ - Le * Le) / (2.0 * b_ * c_)
                tot += math.atan2(math.sqrt(max(R * R - rho * rho, 0.0)),
                                  rho * cg)
            edge_dihedral_formula[e] = tot
    groups = {}
    for e, delta in edge_dihedral.items():
        fid = edge_fid(*e)
        pair = tuple(sorted((ptype_of_face[edge_faces[e][0]],
                             ptype_of_face[edge_faces[e][1]])))
        groups.setdefault((fid if fid >= 0 else 10 ** 9, pair), []).append(e)
    dihedral_rows = []
    for key in sorted(groups):
        es = groups[key]
        fid = key[0] if key[0] < 10 ** 9 else -1
        vals = np.array([edge_dihedral[e] for e in es])
        if spherical:
            fv = np.array([edge_dihedral_formula[e] for e in es])
            fmean = float(np.mean(fv) * DEG)
            chk = float(np.max(np.abs(vals - fv)) * DEG)
        else:
            fmean = chk = float("nan")
        dihedral_rows.append({
            "strut_family": labels[fid] if fid >= 0 else "?",
            "panels": " | ".join(key[1]),
            "count": len(es),
            "dihedral_deg": float(vals.mean() * DEG),
            "bend_deg": float(180.0 - vals.mean() * DEG),
            "dihedral_spread_deg": float(np.ptp(vals) * DEG),
            "dihedral_formula_deg": fmean,
            "dihedral_check_deg": chk,
        })

    # -------------------------------------------------------------- hubs --
    valence = np.bincount(M.ravel(), minlength=n) if m else np.zeros(n, int)
    if mdl is not None:
        try:
            supported = set(int(v) for v in mdl.support_nodes())
        except Exception:                             # e.g. empty CUSTOM table
            supported = set(int(v) for v in g.base_nodes)
    else:
        supported = set(int(v) for v in g.base_nodes)
    free_edge_nodes = set()
    for e, fl in edge_faces.items():
        if len(fl) == 1:
            free_edge_nodes.update(e)

    def location(i):
        if i in supported:
            return "base"
        if i in free_edge_nodes:
            return "edge"
        return "interior"

    cycles = _hub_cycles(P, M, fam, axes)
    hub_key = []
    hub_gaps = []
    for i in range(n):
        fids, az = cycles[i]
        loc = location(i)
        seq, gaps = _canonical_cycle(fids, az, open_fan=(loc != "interior"))
        hub_key.append((loc, _multiset_label(fids, labels),
                        _seq_label(seq, labels) if seq else "(none)",
                        int(valence[i])))
        hub_gaps.append(np.asarray(gaps, dtype=float))
    hub_groups = {}
    for i, key in enumerate(hub_key):
        hub_groups.setdefault(key, []).append(i)
    order = sorted(hub_groups, key=lambda k: (_LOCATION_ORDER[k[0]], -k[3],
                                              k[1], k[2]))
    hub_label = {}
    hub_rows = []
    for t, key in enumerate(order, start=1):
        nodes_t = hub_groups[key]
        lab = f"H{t}"
        for i in nodes_t:
            hub_label[i] = lab
        rep = nodes_t[0]
        g0 = np.sort(hub_gaps[rep])
        spread = 0.0
        if len(g0):
            for i in nodes_t[1:]:
                spread = max(spread, float(np.max(np.abs(
                    np.sort(hub_gaps[i]) - g0))))
        hub_rows.append({
            "hub_type": lab,
            "location": key[0],
            "valence": key[3],
            "struts": key[1],
            "arrangement": key[2],
            "count": len(nodes_t),
            "gaps_deg": " ".join(f"{x * DEG:.2f}" for x in hub_gaps[rep]),
            "gap_spread_deg": spread * DEG,
            "example_node": int(rep),
        })

    # ------------------------------------------------------------ struts --
    if mdl is not None:
        rho = float(mdl.material.rho)
        secs = [mdl.section_of_group(gid) for gid in range(nfam)]
        sec_names = list(mdl.group_sections)
    else:
        rho = float("nan")
        secs = [None] * nfam
        sec_names = [""] * nfam
    # Cut groups: the analysis groups struts into families with the geometry
    # tolerance (1e-4 R for a geodesic dome, i.e. 2 mm on a 20 m sphere),
    # which can lump struts that differ by more than a shop tolerance at
    # high frequencies.  Each family is split into cut groups whose sorted
    # lengths have no gap larger than ``length_tol`` (labels U.1, U.2, ...
    # when a family splits).  Hubs, panels and sections keep the family.
    strut_rows = []
    cut_group = [""] * m
    for gid in range(nfam):
        sel = np.where(fam == gid)[0]
        if len(sel) == 0:
            continue
        srt = sel[np.argsort(L[sel], kind="stable")]
        clusters, start = [], 0
        for t_ in range(1, len(srt) + 1):
            if t_ == len(srt) or L[srt[t_]] - L[srt[t_ - 1]] > length_tol:
                clusters.append(np.sort(srt[start:t_]))
                start = t_
        if len(clusters) > 1:
            span = float(np.ptp(L[sel]))
            notes.append(
                f"Strut family {labels[gid]} has {len(clusters)} distinct "
                f"lengths (spread {span * 1e3:.2f} mm, more than the "
                f"{length_tol * 1e3:.2f} mm cut tolerance): listed as cut "
                f"groups {labels[gid]}.1 to {labels[gid]}.{len(clusters)}.")
        for ci, cl in enumerate(clusters):
            lab = labels[gid] if len(clusters) == 1 else \
                f"{labels[gid]}.{ci + 1}"
            for k in cl:
                cut_group[int(k)] = lab
            strut_rows.append(_strut_row(
                lab, len(strut_rows), cl, L, o, ax_i, ax_j, ax_f, R,
                spherical, secs[gid], sec_names[gid], rho, length_tol, notes))

    # ------------------------------------------------- assembly lists -----
    member_rows = []
    for k, (i, j) in enumerate(M):
        e = (min(int(i), int(j)), max(int(i), int(j)))
        fl = edge_faces.get(e, [])
        member_rows.append({
            "member": k, "node_i": int(i), "node_j": int(j),
            "family": labels[int(fam[k])],
            "cut_group": cut_group[k],
            "length_mm": float(L[k] * 1e3),
            "cut_length_mm": float((L[k] - 2.0 * o) * 1e3),
            "hub_i": hub_label[int(i)], "hub_j": hub_label[int(j)],
            "axial_i_deg": float(ax_i[k] * DEG),
            "axial_j_deg": float(ax_j[k] * DEG),
            "panels": len(fl),
            "dihedral_deg": (float(edge_dihedral[e] * DEG)
                             if e in edge_dihedral else float("nan")),
        })
    node_rows = []
    for i in range(n):
        node_rows.append({
            "node": i, "x_m": float(P[i, 0]), "y_m": float(P[i, 1]),
            "z_m": float(P[i, 2]), "location": hub_key[i][0],
            "valence": int(valence[i]), "hub_type": hub_label[i],
            "struts": hub_key[i][1],
        })
    face_rows = []
    for f, (a, b, c) in enumerate(Fo):
        face_rows.append({
            "face": f, "node_1": int(a), "node_2": int(b), "node_3": int(c),
            "panel_type": ptype_of_face[f], "area_m2": float(area[f]),
        })

    # -------------------------------------------------------- stock plan --
    tables_stock = None
    if stock_length is not None:
        pieces = [(r["family"], r["section"] or "(section)",
                   r["cut_length_mm"] / 1e3, r["count"]) for r in strut_rows]
        fam_rows, pat_rows, sec_rows, snotes = stock_plan(
            pieces, stock_length, kerf)
        notes += snotes
        tables_stock = (fam_rows, pat_rows, sec_rows)

    # ----------------------------------------------------------- summary --
    n_panels_of = [len(edge_faces.get(e, ())) for e in edge_member]
    n_interior_edges = sum(1 for v in n_panels_of if v == 2)
    n_boundary = sum(1 for v in n_panels_of if v == 1)
    n_bare = sum(1 for v in n_panels_of if v == 0)
    clad = float(area.sum())
    total_L = float(L.sum())
    total_cut = float(total_L - 2.0 * o * m)
    mass_total = float(sum(r["mass_total_kg"] for r in strut_rows))
    ax_dev = float(np.max(np.abs(np.concatenate([ax_i - ax_f, ax_j - ax_f])))
                   ) * DEG if spherical and m else float("nan")
    dih_dev = (max((abs(edge_dihedral[e] - edge_dihedral_formula[e])
                    for e in edge_dihedral), default=0.0) * DEG
               if spherical else float("nan"))

    def row(q, v, u="", note=""):
        return {"quantity": q, "value": v, "unit": u, "note": note}

    summary = [
        row("structure", g.title or ("custom truss" if g.is_custom
                                     else "geodesic dome")),
        row("hubs (nodes)", n, "-"),
        row("struts (members)", m, "-"),
        row("panels (faces)", nf, "-"),
        row("Euler characteristic V - E + F", n - m + nf, "-",
            "2 for a closed sphere, 1 for a dome cap"),
        row("strut families", nfam, "-", "analysis families"),
        row("cut groups", len(strut_rows), "-",
            f"strut schedule rows (cut tolerance {length_tol * 1e3:g} mm)"),
        row("hub types", len(hub_rows), "-"),
        row("panel types", len(panel_rows), "-"),
        row("hub offset", o * 1e3, "mm", "hub centre to tube end"),
        row("total strut length (hub to hub)", total_L, "m"),
        row("total cut length", total_cut, "m"),
        row("total tube mass (cut lengths)", mass_total, "kg",
            "hubs and connectors not included"),
        row("cladding area (sum of panels)", clad, "m2"),
        row("faceted surface area (geometry)", float(g.surface_area), "m2",
            "sum of the nodal tributary areas"),
        row("struts with two panels / one panel / none",
            f"{n_interior_edges} / {n_boundary} / {n_bare}", "-"),
        row("sphere radius", R, "m", "best-fit sphere for a custom truss"),
        row("sphere centre z", float(centre[2]), "m"),
        row("max hub distance from the sphere", dev, "m"),
        row("max axial angle: vectors - arccos(L/2R)", ax_dev, "deg",
            "closed-form check"),
        row("max dihedral: vectors - closed form", dih_dev, "deg",
            "closed-form check"),
    ]
    if tables_stock is not None:
        summary.append(row("stock length", float(stock_length) * 1e3, "mm"))
        summary.append(row("saw kerf", float(kerf) * 1e3, "mm"))
        for s in tables_stock[2]:
            summary.append(row(
                f"bars of {s['section']}", s["bars"], "-",
                f"{s['plan']}; single family {s['bars_single_family']}, "
                f"mixed {s['bars_mixed']}, lower bound {s['lower_bound']}, "
                f"waste {s['waste_pct']:.1f} %"))

    out = {
        "summary": summary,
        "struts": strut_rows,
        "hubs": hub_rows,
        "panels": panel_rows,
        "dihedrals": dihedral_rows,
        "members": member_rows,
        "nodes": node_rows,
        "faces": face_rows,
    }
    if tables_stock is not None:
        out["stock"] = tables_stock[0]
        out["stock_patterns"] = tables_stock[1]
    out["notes"] = [{"note": s} for s in notes]
    return out


# ================================================================ export ====
def _fmt(key, value):
    """CSV cell: NaN -> empty, floats rounded by the unit in the column name."""
    if isinstance(value, (bool, np.bool_)):
        return "yes" if value else "no"
    if isinstance(value, (int, np.integer)):
        return int(value)
    if isinstance(value, (float, np.floating)):
        v = float(value)
        if not math.isfinite(v):
            return ""
        if key == "chord_factor":
            return f"{v:.8f}"
        for suffix, dp in (("_mm", 2), ("_m2", 4), ("_m", 4), ("_deg", 4),
                           ("_kg", 3), ("_pct", 2)):
            if key.endswith(suffix):
                return f"{v:.{dp}f}"
        return f"{v:.10g}"
    return value


def export_fabrication_csv(tables, folder, prefix="dome", **kwargs):
    """Write every fabrication table as ``<prefix>_fab_<table>.csv``.

    ``tables`` is the dict of :func:`fabrication_tables`, or a DomeModel /
    DomeGeometry (then ``kwargs`` go to :func:`fabrication_tables`).  Floats
    are rounded by unit (mm 2 dp, m 4 dp, deg 4 dp, kg 3 dp, chord factor
    8 dp); NaN is written as an empty cell.  Returns the written paths in
    :data:`TABLE_ORDER` order.  Empty tables are skipped.
    """
    if not isinstance(tables, dict):
        tables = fabrication_tables(tables, **kwargs)
    os.makedirs(folder, exist_ok=True)
    written = []
    for name in TABLE_ORDER:
        rows = tables.get(name)
        if not rows:
            continue
        path = os.path.join(folder, f"{prefix}_fab_{name}.csv")
        cols = list(rows[0].keys())
        with open(path, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            w.writerow(cols)
            for r in rows:
                w.writerow([_fmt(c, r.get(c, "")) for c in cols])
        written.append(path)
    return written


def fabrication_text(tables, max_rows=40):
    """Plain-text rendering of the summary, strut, hub and panel schedules."""
    def fmt_cell(k, v):
        s = _fmt(k, v)
        return str(s)

    out = []
    for name in ("summary", "struts", "hubs", "panels", "dihedrals",
                 "stock_patterns", "notes"):
        rows = tables.get(name) or []
        if not rows:
            continue
        cols = list(rows[0].keys())
        if name == "hubs":
            cols = [c for c in cols if c != "gaps_deg"]
        cells = [[fmt_cell(c, r.get(c, "")) for c in cols]
                 for r in rows[:max_rows]]
        width = [max(len(c), *(len(row[i]) for row in cells))
                 for i, c in enumerate(cols)]
        out.append(TABLE_TITLES.get(name, name))
        out.append("-" * len(out[-1]))
        out.append("  ".join(c.ljust(w) for c, w in zip(cols, width)))
        for row in cells:
            out.append("  ".join(s.ljust(w) for s, w in zip(row, width)))
        if len(rows) > max_rows:
            out.append(f"... {len(rows) - max_rows} more rows")
        out.append("")
    return "\n".join(out)
