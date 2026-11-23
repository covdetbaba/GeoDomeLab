"""
Class I (icosahedral, "alternate") geodesic sphere and dome generation.

Pipeline
--------
1.  Build the regular icosahedron in *pole orientation* (a vertex at each
    pole, two staggered rings of five vertices at z = +-1/sqrt(5)).
2.  Subdivide every one of the 20 spherical triangles into ``freq**2`` smaller
    triangles using barycentric interpolation, then radially project each new
    point onto the unit sphere.  This is the Class I / "alternate" breakdown.
3.  Merge duplicated vertices that are shared between adjacent faces.
4.  Scale to the requested radius, truncate at a horizontal plane that is
    snapped onto an existing ring of nodes, and translate so the base ring
    sits at z = 0.
5.  Extract the unique strut (edge) set and group struts by length, which
    yields the classic geodesic *chord factor* schedule.

Vertex count of a complete Class I sphere of frequency ``v`` is ``10*v**2 + 2``
and the strut count is ``30*v**2`` -- both are asserted in the test-suite.

General trusses
---------------
:func:`custom_geometry` wraps any pin-jointed space truss (nodes, members,
optional cladding faces) in the same :class:`DomeGeometry` container, so the
solver, loads, design checks and reports work unchanged.  Its ``radius`` is
the least-squares best-fit sphere radius (NaN for a non-spherical truss).

Supports
--------
:func:`support_conditions` returns rigid restraints and grounded springs for
the four support types PINNED, ROLLER3, ELASTIC (springs kx, ky, kz at every
base node) and CUSTOM (per-node fixed / free / spring table).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

try:                                   # SciPy is used only for fast merging
    from scipy.spatial import cKDTree
    _HAVE_KDTREE = True
except Exception:                      # pragma: no cover - SciPy always present
    _HAVE_KDTREE = False


# ============================================================== icosahedron ==
def icosahedron():
    """Return ``(vertices, faces)`` of a unit icosahedron in pole orientation.

    vertices : (12, 3) float array, every row a unit vector
    faces    : (20, 3) int array of vertex indices
    """
    z_ring = 1.0 / np.sqrt(5.0)
    r_ring = 2.0 / np.sqrt(5.0)

    verts = [np.array([0.0, 0.0, 1.0])]                     # 0  : north pole
    for i in range(5):                                      # 1-5: upper ring
        a = np.deg2rad(72.0 * i)
        verts.append(np.array([r_ring * np.cos(a), r_ring * np.sin(a), z_ring]))
    for i in range(5):                                      # 6-10: lower ring
        a = np.deg2rad(72.0 * i + 36.0)
        verts.append(np.array([r_ring * np.cos(a), r_ring * np.sin(a), -z_ring]))
    verts.append(np.array([0.0, 0.0, -1.0]))                # 11 : south pole

    U = lambda i: 1 + (i % 5)        # noqa: E731  upper-ring index helper
    L = lambda i: 6 + (i % 5)        # noqa: E731  lower-ring index helper

    faces = []
    for i in range(5):                                   # 5 polar caps (top)
        faces.append([0, U(i), U(i + 1)])
    for i in range(5):                                   # 10 equatorial band
        faces.append([U(i), L(i), U(i + 1)])
        faces.append([L(i), L(i + 1), U(i + 1)])
    for i in range(5):                                   # 5 polar caps (bottom)
        faces.append([11, L(i + 1), L(i)])

    return np.asarray(verts, dtype=float), np.asarray(faces, dtype=int)


# ================================================================ merging ====
def merge_vertices(points, tol=1e-7):
    """Collapse coincident points.

    Returns ``(unique_points, index_map)`` where ``index_map[i]`` is the row of
    ``unique_points`` that original point ``i`` was merged into.
    """
    points = np.asarray(points, dtype=float)
    n = len(points)
    index_map = np.full(n, -1, dtype=int)
    unique = []

    if _HAVE_KDTREE and n > 1:
        tree = cKDTree(points)
        groups = tree.query_ball_point(points, r=tol)
        for i in range(n):
            if index_map[i] != -1:
                continue
            new_id = len(unique)
            unique.append(points[i])
            for j in groups[i]:
                if index_map[j] == -1:
                    index_map[j] = new_id
    else:                                              # pragma: no cover
        bucket = {}
        decimals = max(0, int(round(-np.log10(tol))) - 1)
        for i, p in enumerate(points):
            key = tuple(np.round(p, decimals))
            if key not in bucket:
                bucket[key] = len(unique)
                unique.append(p)
            index_map[i] = bucket[key]

    return np.asarray(unique, dtype=float), index_map


# ============================================================ subdivision ====
def subdivide_class1(vertices, faces, freq):
    """Class I subdivision of a triangulated sphere, projected to unit radius.

    Every input face is split into ``freq**2`` triangles.  Node ``P(i, j)`` of
    the triangle (A, B, C) is the barycentric point

        P(i, j) = [ (v - i)*A + (i - j)*B + j*C ] / v ,
        0 <= j <= i <= v
    """
    freq = int(freq)
    if freq < 1:
        raise ValueError("Geodesic frequency must be >= 1")

    raw_pts, raw_tris = [], []
    for fa, fb, fc in faces:
        A, B, C = vertices[fa], vertices[fb], vertices[fc]
        # Row-major index of P(i, j) inside this face's local point list.
        local = {}
        for i in range(freq + 1):
            for j in range(i + 1):
                p = ((freq - i) * A + (i - j) * B + j * C) / float(freq)
                p = p / np.linalg.norm(p)          # radial projection
                local[(i, j)] = len(raw_pts)
                raw_pts.append(p)
        for i in range(freq):
            for j in range(i + 1):
                raw_tris.append([local[(i, j)], local[(i + 1, j)],
                                 local[(i + 1, j + 1)]])
                if j < i:
                    raw_tris.append([local[(i, j)], local[(i + 1, j + 1)],
                                     local[(i, j + 1)]])

    raw_pts = np.asarray(raw_pts, dtype=float)
    pts, imap = merge_vertices(raw_pts, tol=1e-6)
    tris = imap[np.asarray(raw_tris, dtype=int)]

    # Drop any degenerate triangle produced by the merge (should not happen).
    keep = np.array([len(set(t)) == 3 for t in tris], dtype=bool)
    tris = tris[keep]

    # De-duplicate faces that may be shared between neighbouring base faces.
    seen, unique_tris = set(), []
    for t in tris:
        key = tuple(sorted(t))
        if key not in seen:
            seen.add(key)
            unique_tris.append(t)
    return pts, np.asarray(unique_tris, dtype=int)


def edges_from_faces(faces):
    """Unique undirected edge list ``(m, 2)`` with ``i < j`` in every row."""
    s = set()
    for a, b, c in faces:
        for i, j in ((a, b), (b, c), (c, a)):
            s.add((min(i, j), max(i, j)))
    return np.asarray(sorted(s), dtype=int)


# ============================================================= containers ====
@dataclass
class DomeGeometry:
    """Everything purely geometric about one generated dome."""

    frequency: int
    radius: float
    dome_fraction: float

    nodes: np.ndarray                 # (n, 3) coordinates, base ring at z = 0
    members: np.ndarray               # (m, 2) node indices
    faces: np.ndarray                 # (f, 3) node indices (cladding panels)
    normals: np.ndarray               # (n, 3) outward unit surface normals
    base_nodes: np.ndarray            # (nb,) indices of the supported ring
    crown_node: int                   # index of the highest node
    ring_index: np.ndarray            # (n,) ring number, 0 = crown ring
    ring_heights: np.ndarray          # (nr,) z of each ring, descending
    strut_group: np.ndarray           # (m,) strut-type id 0 = shortest
    strut_lengths: np.ndarray         # (m,) chord lengths
    trib_area: np.ndarray             # (n,) surface tributary area  [m^2]
    trib_area_proj: np.ndarray        # (n,) projected tributary area [m^2]
    notes: list = field(default_factory=list)
    # "geodesic" for generate_dome(), "custom" for custom_geometry().  For a
    # custom truss frequency = 0, dome_fraction = NaN and radius is the
    # best-fit sphere radius (NaN when the points do not define a sphere).
    kind: str = "geodesic"
    title: str = ""

    # ------------------------------------------------------------- queries --
    @property
    def n_nodes(self):
        return len(self.nodes)

    @property
    def n_members(self):
        return len(self.members)

    @property
    def n_faces(self):
        return len(self.faces)

    @property
    def is_custom(self):
        return self.kind == "custom"

    @property
    def height(self):
        """Rise of the structure: crown elevation above the lowest support."""
        if self.is_custom:
            zb = (self.nodes[self.base_nodes, 2] if len(self.base_nodes)
                  else self.nodes[:, 2])
            return float(self.nodes[:, 2].max() - zb.min())
        return float(self.nodes[:, 2].max())

    @property
    def base_diameter(self):
        """Span used for the deflection limit.

        Geodesic dome: diameter of the base ring (centred on the z axis).
        Custom truss: largest horizontal distance between two supported
        nodes (the plan diameter of the support set), 0 if fewer than two.
        """
        b = self.nodes[self.base_nodes]
        if self.is_custom:
            if len(b) < 2:
                return 0.0
            d = b[:, None, :2] - b[None, :, :2]
            return float(np.sqrt((d ** 2).sum(axis=2)).max())
        return 2.0 * float(np.max(np.hypot(b[:, 0], b[:, 1])))

    @property
    def length_scale(self):
        """Characteristic size for plot scaling [m].

        The sphere radius for a geodesic dome (and a spherical custom truss);
        half the largest bounding-box side when the radius is NaN.
        """
        R = self.radius
        if np.isfinite(R) and R > 0.0:
            return R
        span = float(np.ptp(self.nodes, axis=0).max())
        return 0.5 * span if span > 0.0 else 1.0

    @property
    def surface_area(self):
        return float(self.trib_area.sum())

    @property
    def plan_area(self):
        return float(self.trib_area_proj.sum())

    def member_lengths(self):
        p = self.nodes[self.members[:, 0]]
        q = self.nodes[self.members[:, 1]]
        return np.linalg.norm(q - p, axis=1)

    def strut_schedule(self):
        """Per-strut-type summary: id, label, length, chord factor, count."""
        out = []
        for g in range(int(self.strut_group.max()) + 1):
            sel = self.strut_group == g
            L = float(np.mean(self.strut_lengths[sel]))
            R = self.radius
            out.append({
                "id": g,
                "label": chr(ord("A") + g) if g < 26 else f"S{g}",
                "length": L,
                # NaN for a custom truss with no defining sphere
                "chord_factor": (L / R if np.isfinite(R) and R > 0.0
                                 else float("nan")),
                "count": int(sel.sum()),
            })
        return out

    def n_strut_types(self):
        return int(self.strut_group.max()) + 1


# ============================================================== generation ===
def _cluster_1d(values, tol):
    """Sort values and return (cluster_id_per_value, cluster_representatives)."""
    order = np.argsort(values)
    ids = np.empty(len(values), dtype=int)
    reps = []
    current = -1
    last = None
    for k in order:
        v = values[k]
        if last is None or abs(v - last) > tol:
            current += 1
            reps.append(v)
            last = v
        else:
            # running mean keeps the representative centred
            reps[current] = 0.5 * (reps[current] + v)
            last = v
        ids[k] = current
    return ids, np.asarray(reps, dtype=float)


def _tributary_areas(nodes, faces, centre=None):
    """One third of each adjacent panel area, surface and horizontal-projected.

    The projected (plan) area of a panel is ``A max(n_z, 0)`` with ``n`` its
    outward unit normal, so a panel below the equator of a super-hemisphere,
    which faces down, takes no snow or live load and adds nothing to the plan
    area.  The outward side is the side away from ``centre`` (the sphere
    centre); without a centre the plan area is ``A |n_z|`` (custom trusses
    whose faces have no defined outward side).
    """
    n = len(nodes)
    surf = np.zeros(n)
    proj = np.zeros(n)
    if centre is not None:
        centre = np.asarray(centre, dtype=float)
        if not np.all(np.isfinite(centre)):
            centre = None
    for a, b, c in faces:
        u = nodes[b] - nodes[a]
        v = nodes[c] - nodes[a]
        cross = np.cross(u, v)
        area = 0.5 * np.linalg.norm(cross)
        if centre is None:
            area_p = 0.5 * abs(cross[2])
        else:
            cen = (nodes[a] + nodes[b] + nodes[c]) / 3.0
            if float(np.dot(cross, cen - centre)) < 0.0:
                cross = -cross                  # orient outward
            area_p = 0.5 * max(float(cross[2]), 0.0)
        for k in (a, b, c):
            surf[k] += area / 3.0
            proj[k] += area_p / 3.0
    return surf, proj


def boundary_nodes(faces):
    """Nodes on the free edge of a triangulated patch.

    An edge that belongs to exactly one triangle is a boundary edge; the nodes
    on those edges are the ones that have to be supported.  This is the only
    definition that stays correct for an arbitrary truncation: taking "every
    node at the cut elevation" instead leaves the real free edge of the mesh
    dangling whenever the cut does not land on a complete horizontal ring.
    """
    count = {}
    for a, b, c in faces:
        for i, j in ((a, b), (b, c), (c, a)):
            key = (min(i, j), max(i, j))
            count[key] = count.get(key, 0) + 1
    nodes = set()
    for (i, j), n in count.items():
        if n == 1:
            nodes.add(i)
            nodes.add(j)
    return np.asarray(sorted(nodes), dtype=int)


def _truncate_at(pts, tris, z_cut, ztol):
    """Keep the panels entirely above ``z_cut``; return (faces, boundary)."""
    keep = np.array([np.all(pts[t][:, 2] >= z_cut - ztol) for t in tris])
    kept = tris[keep]
    if len(kept) == 0:
        return kept, np.empty(0, dtype=int)
    return kept, boundary_nodes(kept)


def generate_dome(frequency=4, radius=10.0, dome_fraction=0.5,
                  strut_tol=1e-4, prefer_level_base=True):
    """Build a truncated Class I geodesic dome.

    Parameters
    ----------
    frequency : int
        Geodesic frequency ``v`` (1 = bare icosahedron, 2V, 3V, ... 10V).
    radius : float
        Sphere radius [m].
    dome_fraction : float
        Fraction of the sphere *diameter* retained, measured from the crown
        down.  ``0.5`` is a hemisphere, ``0.3`` a shallow cap, ``0.7`` a
        "super-hemisphere".  The requested cut is snapped onto the nearest
        existing ring of nodes so the base ring is always a closed polygon of
        real nodes.
    strut_tol : float
        Relative tolerance used when grouping struts into length families.

    Returns
    -------
    DomeGeometry
    """
    try:
        frequency = int(frequency)
        radius = float(radius)
        fraction = float(dome_fraction)
    except (TypeError, ValueError):
        raise ValueError("The frequency, radius and dome fraction must be "
                         "numbers.") from None
    if frequency < 1:
        raise ValueError(f"The frequency must be 1 or more (got {frequency}).")
    if not (math.isfinite(radius) and radius > 0.0):
        raise ValueError(f"The sphere radius must be above 0 m (got "
                         f"{radius:g} m).")
    if not (math.isfinite(fraction) and 0.0 < fraction <= 1.0):
        raise ValueError(f"The dome fraction must be above 0 and at most 1 "
                         f"(got {fraction:g}).")
    notes = []

    verts, faces0 = icosahedron()
    if frequency == 1:
        pts, tris = verts.copy(), faces0.copy()
    else:
        pts, tris = subdivide_class1(verts, faces0, frequency)

    pts = pts * radius

    # ---- ring detection ----------------------------------------------------
    ztol = max(1e-9, 1e-4 * radius)
    _, ring_z_asc = _cluster_1d(pts[:, 2], ztol)
    ring_z = ring_z_asc[::-1]                       # descending: crown first

    # ---- choose the cut level ----------------------------------------------
    # Every candidate level is actually truncated and its free edge inspected,
    # because only the resulting boundary tells us whether the base will be a
    # level ring or a stepped (zigzag) edge.
    z_target = radius * (1.0 - 2.0 * float(dome_fraction))
    candidates = []
    for k in range(1, len(ring_z)):
        kept, bnd = _truncate_at(pts, tris, ring_z[k], ztol)
        if len(kept) == 0 or len(bnd) < 3:
            continue
        zb = pts[bnd]
        step = float(zb[:, 2].max() - zb[:, 2].min())
        # Circumradius of the free edge: a cut so low that the dome is almost a
        # closed sphere leaves a tiny hole that cannot sensibly be supported.
        r_edge = float(np.max(np.hypot(zb[:, 0], zb[:, 1])))
        if r_edge < 0.25 * radius:
            continue
        candidates.append({"k": k, "z": float(ring_z[k]),
                           "n_boundary": int(len(bnd)),
                           "step": step, "level": step <= 10.0 * ztol,
                           "z_low": float(zb[:, 2].min())})
    if not candidates:
        raise ValueError("No valid truncation level; use frequency >= 2.")

    def miss(c):
        return abs(c["z_low"] - z_target)

    best = min(candidates, key=miss)
    if prefer_level_base:
        level_only = [c for c in candidates if c["level"]]
        if level_only:
            best_level = min(level_only, key=miss)
            # Take the level base only when it costs little: a flat springing
            # line is worth a small change of rise, not a different dome.
            if miss(best_level) <= miss(best) + 0.06 * radius:
                best = best_level

    z_cut = best["z"]
    z_base = best["z_low"]
    actual_fraction = 0.5 * (1.0 - z_base / radius)

    if abs(actual_fraction - dome_fraction) > 0.02:
        notes.append(
            f"Requested dome fraction {dome_fraction:.3f} snapped to "
            f"{actual_fraction:.3f}: a geodesic mesh can only be cut along its "
            f"own strut lines, and {len(candidates)} such level(s) exist at "
            f"{frequency}V.")
    if not best["level"]:
        notes.append(
            f"The base edge is stepped by {best['step']:.3f} m rather than "
            f"level; the supports follow the strut lines. Use an even "
            f"frequency with a 0.50 fraction for a flat springing line.")

    # ---- truncation --------------------------------------------------------
    tris, _ = _truncate_at(pts, tris, z_cut, ztol)

    used = np.unique(tris)
    remap = -np.ones(len(pts), dtype=int)
    remap[used] = np.arange(len(used))
    nodes = pts[used]
    faces = remap[tris]

    normals = nodes / np.linalg.norm(nodes, axis=1, keepdims=True)

    members = edges_from_faces(faces)

    # ---- supports = the free edge of the retained cap ----------------------
    base_nodes = boundary_nodes(faces)
    if len(base_nodes) < 3:
        raise ValueError("The truncated cap has no usable free edge; "
                         "choose a different dome fraction.")

    # Any node left with fewer than three struts is a mechanism unless it is
    # supported, so fold those into the support set and say so.
    degree = np.zeros(len(nodes), dtype=int)
    for i, j in members:
        degree[i] += 1
        degree[j] += 1
    weak = np.setdiff1d(np.where(degree < 3)[0], base_nodes)
    if len(weak):
        base_nodes = np.union1d(base_nodes, weak)
        notes.append(f"{len(weak)} node(s) with fewer than three struts were "
                     f"added to the support set to keep the model stable.")

    _, ring_rep_asc = _cluster_1d(nodes[:, 2], ztol)
    ring_rep = ring_rep_asc[::-1]
    ring_index = np.array([int(np.argmin(np.abs(ring_rep - z)))
                           for z in nodes[:, 2]], dtype=int)

    crown_node = int(np.argmax(nodes[:, 2]))

    # ---- strut families ----------------------------------------------------
    L = np.linalg.norm(nodes[members[:, 1]] - nodes[members[:, 0]], axis=1)
    strut_group, _ = _cluster_1d(L, strut_tol * radius)

    # ---- translate so the lowest support sits on the ground ----------------
    nodes = nodes.copy()
    nodes[:, 2] -= z_base
    ring_rep = ring_rep - z_base

    # sphere centre after the translation: (0, 0, -z_base)
    trib, trib_p = _tributary_areas(nodes, faces,
                                    centre=(0.0, 0.0, -float(z_base)))

    notes.append(f"{len(nodes)} nodes, {len(members)} struts, {len(faces)} panels.")
    notes.append(f"{int(strut_group.max()) + 1} distinct strut length(s).")

    return DomeGeometry(
        frequency=frequency, radius=radius, dome_fraction=actual_fraction,
        nodes=nodes, members=members, faces=faces, normals=normals,
        base_nodes=base_nodes, crown_node=crown_node,
        ring_index=ring_index, ring_heights=ring_rep,
        strut_group=strut_group, strut_lengths=L,
        trib_area=trib, trib_area_proj=trib_p, notes=notes,
        kind="geodesic", title=f"{frequency}V Class I geodesic dome")


# ======================================================= custom geometry ====
def best_fit_sphere(points, rel_tol=0.05):
    """Least-squares sphere through a point cloud.

    The algebraic (Kasa) fit writes the sphere |x - c|^2 = R^2 as the linear
    equation

        |x|^2 = 2 c . x + d ,      d = R^2 - |c|^2

    and solves it for (c, d) in the least-squares sense; R = sqrt(d + |c|^2).
    The points are centred first for conditioning.

    Returns ``(centre, R, rms)`` where ``rms`` is the RMS of |x - c| - R.  The
    radius (and the centre) is NaN when the points cannot define a sphere:
    fewer than four points, all points coplanar (or collinear), or a fit whose
    RMS distance error exceeds ``rel_tol * R`` (the structure is not
    spherical, so a "sphere radius" would be meaningless).
    """
    nan3 = np.full(3, np.nan)
    P = np.asarray(points, dtype=float)
    if P.ndim != 2 or P.shape[1] != 3 or len(P) < 4:
        return nan3, float("nan"), float("nan")
    mean = P.mean(axis=0)
    Q = P - mean
    size = float(np.max(np.linalg.norm(Q, axis=1)))
    if size <= 0.0:
        return nan3, float("nan"), float("nan")
    sv = np.linalg.svd(Q / size, compute_uv=False)
    if sv[-1] < 1e-9 * sv[0]:                       # coplanar / collinear
        return nan3, float("nan"), float("nan")
    A = np.column_stack([2.0 * Q, np.ones(len(Q))])
    b = np.sum(Q ** 2, axis=1)
    sol, *_ = np.linalg.lstsq(A, b, rcond=None)
    c = sol[:3]
    R2 = sol[3] + float(c @ c)
    if not np.isfinite(R2) or R2 <= 0.0:
        return nan3, float("nan"), float("nan")
    R = float(np.sqrt(R2))
    rms = float(np.sqrt(np.mean((np.linalg.norm(Q - c, axis=1) - R) ** 2)))
    if rms > rel_tol * R:
        return nan3, float("nan"), rms
    return c + mean, R, rms


def _face_normals_at_nodes(nodes, faces, ref_point):
    """Area-weighted outward unit normal at every node (zeros if no faces).

    Each triangle normal n_f = (b - a) x (c - a) (length = twice its area) is
    oriented to point away from ``ref_point`` (the best-fit sphere centre, or
    the centroid of the nodes); when the face is edge-on to that point the
    upward (+z) orientation is used.  The node normal is the normalised sum of
    the oriented normals of its adjacent faces.
    """
    n = len(nodes)
    acc = np.zeros((n, 3))
    if len(faces) == 0:
        return acc
    scale = float(np.ptp(nodes, axis=0).max()) or 1.0
    for a, b, c in faces:
        nf = np.cross(nodes[b] - nodes[a], nodes[c] - nodes[a])
        mag = np.linalg.norm(nf)
        if mag <= 0.0:
            continue
        cen = (nodes[a] + nodes[b] + nodes[c]) / 3.0
        s = float(np.dot(nf, cen - ref_point))
        if abs(s) > 1e-6 * mag * scale:
            if s < 0.0:
                nf = -nf
        elif nf[2] < 0.0:
            nf = -nf
        for k in (a, b, c):
            acc[k] += nf
    nrm = np.linalg.norm(acc, axis=1)
    ok = nrm > 0.0
    acc[ok] /= nrm[ok][:, None]
    return acc


def custom_geometry(nodes, members, faces=None, base_nodes=None, groups=None,
                    title="Custom truss", strut_tol=1e-4):
    """Wrap an arbitrary pin-jointed space truss in a :class:`DomeGeometry`.

    Parameters
    ----------
    nodes : (n, 3) array_like
        Node coordinates [m].  They are used exactly as given (no translation).
    members : (m, 2) array_like of int
        Node indices of every bar (zero based, no repeats, no self loops).
    faces : (f, 3) array_like of int, optional
        Cladding triangles.  When given, the surface normals and the surface /
        horizontally projected tributary areas (one third of each adjacent
        triangle) are computed from them, so dead, snow, live and wind
        pressures act on the structure exactly as on a generated dome.  When
        omitted, normals and tributary areas are zero: only self weight and
        explicit nodal loads (``LoadSpec.nodal_loads``) load the truss.
    base_nodes : array_like of int, optional
        Supported nodes.  Default: every node within 1e-6 x (model size) of
        the lowest elevation.
    groups : array_like of int, optional
        Strut family of every member (any integer labels; they are renumbered
        0, 1, 2, ... in ascending label order).  Default: families by length,
        exactly as for a geodesic dome (bars whose lengths differ by less than
        ``strut_tol`` x the longest bar share a family, 0 = shortest).
    title : str
        Name shown in summaries and reports.

    The returned geometry has ``kind = "custom"``, ``frequency = 0``,
    ``dome_fraction = NaN`` and ``radius`` = the best-fit sphere radius of the
    nodes (:func:`best_fit_sphere`), NaN when the nodes are not spherical; the
    chord factors are then NaN as well.  ``crown_node`` is the highest node.
    """
    notes = []
    P = np.array(nodes, dtype=float)
    if P.ndim != 2 or P.shape[1] != 3 or len(P) < 2:
        raise ValueError("nodes must be an (n, 3) array with n >= 2.")
    if not np.all(np.isfinite(P)):
        raise ValueError("Node coordinates must be finite numbers.")
    n = len(P)

    M = np.array(members, dtype=float)
    if M.ndim != 2 or M.shape[1] != 2 or len(M) < 1:
        raise ValueError("members must be an (m, 2) array with m >= 1.")
    if not np.all(M == np.round(M)):
        raise ValueError("Member node indices must be integers.")
    M = M.astype(int)
    if M.min() < 0 or M.max() >= n:
        raise ValueError(f"A member refers to a node outside 0..{n - 1}.")
    if np.any(M[:, 0] == M[:, 1]):
        k = int(np.where(M[:, 0] == M[:, 1])[0][0])
        raise ValueError(f"Member {k} connects node {M[k, 0]} to itself.")
    _, counts = np.unique(np.sort(M, axis=1), axis=0, return_counts=True)
    if np.any(counts > 1):
        raise ValueError(f"{int(np.sum(counts > 1))} member(s) are defined "
                         f"twice between the same pair of nodes.")

    if faces is None or len(faces) == 0:
        Fc = np.zeros((0, 3), dtype=int)
    else:
        Fc = np.array(faces, dtype=int)
        if Fc.ndim != 2 or Fc.shape[1] != 3:
            raise ValueError("faces must be an (f, 3) array of node indices.")
        if Fc.min() < 0 or Fc.max() >= n:
            raise ValueError(f"A face refers to a node outside 0..{n - 1}.")
        for k, t in enumerate(Fc):
            if len(set(t.tolist())) != 3:
                raise ValueError(f"Face {k} repeats a node.")

    L = np.linalg.norm(P[M[:, 1]] - P[M[:, 0]], axis=1)
    if np.any(L <= 0.0):
        k = int(np.argmin(L))
        raise ValueError(f"Member {k} has zero length (coincident nodes).")
    size = float(np.ptp(P, axis=0).max())

    # ---- supports ----------------------------------------------------------
    if base_nodes is None:
        ztol = max(1e-9, 1e-6 * size)
        zmin = float(P[:, 2].min())
        base = np.where(P[:, 2] <= zmin + ztol)[0].astype(int)
        notes.append(f"Supports default to the {len(base)} lowest node(s).")
    else:
        base = np.unique(np.asarray(base_nodes, dtype=int))
        if len(base) and (base.min() < 0 or base.max() >= n):
            raise ValueError(f"A base node is outside 0..{n - 1}.")

    # ---- sphere, normals, tributary areas ----------------------------------
    centre, R, rms = best_fit_sphere(P)
    if np.isfinite(R):
        notes.append(f"Best-fit sphere radius {R:.4f} m (RMS distance error "
                     f"{rms:.2e} m).")
    else:
        notes.append("The nodes do not define a sphere: radius and chord "
                     "factors are reported as NaN.")
    ref = centre if np.all(np.isfinite(centre)) else P.mean(axis=0)
    normals = _face_normals_at_nodes(P, Fc, ref)
    if len(Fc):
        trib, trib_p = _tributary_areas(P, Fc, centre=centre)
    else:
        trib, trib_p = np.zeros(n), np.zeros(n)
        notes.append("No faces given: area loads (cladding, snow, live, "
                     "wind) are zero; use self weight and nodal loads.")

    # ---- rings (elevation levels) -------------------------------------------
    ztol_r = max(1e-9, 1e-4 * size)
    _, ring_asc = _cluster_1d(P[:, 2], ztol_r)
    ring_rep = ring_asc[::-1]
    ring_index = np.array([int(np.argmin(np.abs(ring_rep - z)))
                           for z in P[:, 2]], dtype=int)
    crown_node = int(np.argmax(P[:, 2]))

    # ---- strut families -----------------------------------------------------
    if groups is None:
        strut_group, _ = _cluster_1d(L, strut_tol * float(L.max()))
    else:
        gr = np.asarray(groups)
        if gr.shape != (len(M),):
            raise ValueError("groups must give one family per member.")
        labels, strut_group = np.unique(gr.astype(int), return_inverse=True)
        strut_group = strut_group.reshape(-1).astype(int)
        if not np.array_equal(labels, np.arange(len(labels))):
            notes.append("Strut family labels renumbered 0.."
                         f"{len(labels) - 1}.")

    nfam = int(strut_group.max()) + 1
    notes.append(f"{n} nodes, {len(M)} struts, {len(Fc)} panels.")
    notes.append(f"{nfam} distinct strut famil{'y' if nfam == 1 else 'ies'}.")

    return DomeGeometry(
        frequency=0, radius=R, dome_fraction=float("nan"),
        nodes=P, members=M, faces=Fc, normals=normals,
        base_nodes=base, crown_node=crown_node,
        ring_index=ring_index, ring_heights=ring_rep,
        strut_group=np.asarray(strut_group, dtype=int), strut_lengths=L,
        trib_area=trib, trib_area_proj=trib_p, notes=notes,
        kind="custom", title=str(title))


# =============================================================== supports ====
PINNED = "pinned"
ROLLER3 = "roller3"
ELASTIC = "elastic"
CUSTOM = "custom"

SUPPORT_TYPES = {
    PINNED: "Pinned base ring (ux, uy, uz fixed at every base node)",
    ROLLER3: "Vertical rollers + three lateral restraints at ~120 deg",
    ELASTIC: "Elastic springs kx, ky, kz at every base node",
    CUSTOM: "Custom per-node supports (fixed / free / spring per direction)",
}

DEFAULT_SPRING_STIFFNESS = (5.0e7, 5.0e7, 5.0e7)     # N/m  (= 50 kN/mm)


def _condition(value, where):
    """Classify one support entry -> ('fixed' | 'free' | 'spring', k)."""
    if isinstance(value, (bool, np.bool_)):
        return ("fixed", 0.0) if bool(value) else ("free", 0.0)
    if value is None:
        return "free", 0.0
    if isinstance(value, str):
        v = value.strip().lower()
        if v in ("fixed", "pinned", "true", "rigid"):
            return "fixed", 0.0
        if v in ("free", "false", "none", ""):
            return "free", 0.0
        try:
            value = float(v)
        except ValueError:
            raise ValueError(f"{where}: unknown support entry "
                             f"{value!r}") from None
    k = float(value)
    if not (k >= 0.0):                    # catches negatives and NaN
        raise ValueError(f"{where}: spring stiffness must be >= 0 N/m, "
                         f"got {value!r}.")
    if k == 0.0:
        return "free", 0.0
    if math.isinf(k):
        return "fixed", 0.0
    return "spring", k


def support_conditions(geom, support_type=PINNED, spring_stiffness=None,
                       custom_supports=None):
    """Rigid restraints and elastic springs for a support type.

    Returns ``(fixed, springs)``: ``fixed`` is the sorted int array of
    restrained DOF numbers (``3*node + {0: ux, 1: uy, 2: uz}``) and
    ``springs`` a ``{dof: k}`` dict of grounded translational springs [N/m]
    (those DOFs stay free in the solver).

    * PINNED  -- ux, uy, uz fixed at every base node.
    * ROLLER3 -- uz fixed at every base node, ux and uy fixed at the three
      base nodes nearest the bearings 0, 120 and 240 deg.
    * ELASTIC -- springs ``spring_stiffness = (kx, ky, kz)`` [N/m] at every
      base node (default :data:`DEFAULT_SPRING_STIFFNESS`; 0 = free in that
      direction, inf = rigid).
    * CUSTOM  -- ``custom_supports = {node: [cx, cy, cz]}`` with each entry
      ``True`` (fixed), ``False`` (free) or a positive stiffness [N/m]
      (``inf`` = fixed).  Only the listed nodes are supported; base nodes
      that are not listed are free.
    """
    base = np.asarray(geom.base_nodes, dtype=int)
    springs = {}

    if support_type in (PINNED, ROLLER3):
        fixed = []
        if support_type == PINNED:
            for nd in base:
                fixed += [3 * nd, 3 * nd + 1, 3 * nd + 2]
        else:
            for nd in base:
                fixed.append(3 * nd + 2)                    # uz everywhere
            b = geom.nodes[base]
            ang = np.degrees(np.arctan2(b[:, 1], b[:, 0])) % 360.0
            for target in (0.0, 120.0, 240.0):
                d = np.abs((ang - target + 180.0) % 360.0 - 180.0)
                nd = int(base[int(np.argmin(d))])
                fixed += [3 * nd, 3 * nd + 1]
        return np.unique(np.asarray(fixed, dtype=int)), springs

    if support_type == ELASTIC:
        ks = DEFAULT_SPRING_STIFFNESS if spring_stiffness is None \
            else tuple(spring_stiffness)
        if len(ks) != 3:
            raise ValueError("spring_stiffness must be (kx, ky, kz) in N/m.")
        fixed = []
        for a, kval in enumerate(ks):
            kind, k = _condition(float(kval), f"spring_stiffness[{a}]")
            for nd in base:
                if kind == "fixed":
                    fixed.append(3 * int(nd) + a)
                elif kind == "spring":
                    springs[3 * int(nd) + a] = k
        return np.unique(np.asarray(fixed, dtype=int)), springs

    if support_type == CUSTOM:
        if not custom_supports:
            raise ValueError(
                "Support type 'custom' needs DomeModel.custom_supports = "
                "{node: [cx, cy, cz]} (True = fixed, False = free, "
                "positive number = spring in N/m); the table is empty.")
        n = len(geom.nodes)
        fixed = []
        for node, conds in custom_supports.items():
            nd = int(node)
            if nd < 0 or nd >= n:
                raise ValueError(f"custom_supports: node {node} is outside "
                                 f"0..{n - 1}.")
            if isinstance(conds, (str, bytes)) or len(conds) != 3:
                raise ValueError(f"custom_supports[{node}] must list three "
                                 f"conditions [cx, cy, cz].")
            for a, c in enumerate(conds):
                kind, k = _condition(c, f"custom_supports[{node}][{a}]")
                if kind == "fixed":
                    fixed.append(3 * nd + a)
                elif kind == "spring":
                    springs[3 * nd + a] = k
        return np.unique(np.asarray(fixed, dtype=int)), springs

    raise ValueError(f"Unknown support type {support_type!r}")


def support_dofs(geom, support_type=PINNED, spring_stiffness=None,
                 custom_supports=None):
    """Return the sorted array of rigidly restrained DOF numbers.

    DOF numbering is ``3*node + {0: ux, 1: uy, 2: uz}``.  Springs (ELASTIC,
    CUSTOM) are not included; see :func:`support_conditions`.
    """
    return support_conditions(geom, support_type, spring_stiffness,
                              custom_supports)[0]
