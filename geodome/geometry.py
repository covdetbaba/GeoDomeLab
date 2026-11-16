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
"""

from __future__ import annotations

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
    def height(self):
        return float(self.nodes[:, 2].max())

    @property
    def base_diameter(self):
        b = self.nodes[self.base_nodes]
        return 2.0 * float(np.max(np.hypot(b[:, 0], b[:, 1])))

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
            out.append({
                "id": g,
                "label": chr(ord("A") + g) if g < 26 else f"S{g}",
                "length": L,
                "chord_factor": L / self.radius,
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


def _tributary_areas(nodes, faces):
    """One third of each adjacent panel area, surface and horizontal-projected."""
    n = len(nodes)
    surf = np.zeros(n)
    proj = np.zeros(n)
    for a, b, c in faces:
        u = nodes[b] - nodes[a]
        v = nodes[c] - nodes[a]
        cross = np.cross(u, v)
        area = 0.5 * np.linalg.norm(cross)
        area_p = 0.5 * abs(cross[2])
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
    frequency = int(frequency)
    radius = float(radius)
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

    trib, trib_p = _tributary_areas(nodes, faces)

    notes.append(f"{len(nodes)} nodes, {len(members)} struts, {len(faces)} panels.")
    notes.append(f"{int(strut_group.max()) + 1} distinct strut length(s).")

    return DomeGeometry(
        frequency=frequency, radius=radius, dome_fraction=actual_fraction,
        nodes=nodes, members=members, faces=faces, normals=normals,
        base_nodes=base_nodes, crown_node=crown_node,
        ring_index=ring_index, ring_heights=ring_rep,
        strut_group=strut_group, strut_lengths=L,
        trib_area=trib, trib_area_proj=trib_p, notes=notes)


# =============================================================== supports ====
PINNED = "pinned"
ROLLER3 = "roller3"

SUPPORT_TYPES = {
    PINNED: "Pinned base ring (ux, uy, uz fixed at every base node)",
    ROLLER3: "Vertical rollers + three lateral restraints at ~120 deg",
}


def support_dofs(geom, support_type=PINNED):
    """Return the sorted array of globally restrained DOF numbers.

    DOF numbering is ``3*node + {0: ux, 1: uy, 2: uz}``.
    """
    fixed = []
    base = np.asarray(geom.base_nodes, dtype=int)

    if support_type == PINNED:
        for nd in base:
            fixed += [3 * nd, 3 * nd + 1, 3 * nd + 2]

    elif support_type == ROLLER3:
        for nd in base:
            fixed.append(3 * nd + 2)                    # uz everywhere
        b = geom.nodes[base]
        ang = np.degrees(np.arctan2(b[:, 1], b[:, 0])) % 360.0
        for target in (0.0, 120.0, 240.0):
            d = np.abs((ang - target + 180.0) % 360.0 - 180.0)
            nd = int(base[int(np.argmin(d))])
            fixed += [3 * nd, 3 * nd + 1]
    else:
        raise ValueError(f"Unknown support type {support_type!r}")

    return np.unique(np.asarray(fixed, dtype=int))
