"""
Three-dimensional space-truss direct stiffness solver.

Element formulation
-------------------
Each strut is a two-node pin-ended bar with three translational degrees of
freedom per node.  With direction cosines

    {c} = (x2 - x1)/L ,  (y2 - y1)/L ,  (z2 - z1)/L

the 6x6 element stiffness matrix in global coordinates is

          EA  |  [c c^T]   -[c c^T] |
    k  =  --  |                     |
          L   | -[c c^T]    [c c^T] |

which is the familiar ``k = T^T k_local T`` with the 1x6 transformation
``T = [-c^T  c^T]`` and ``k_local = EA/L``.

Solution
--------
The global matrix is assembled in sparse COO form, partitioned into free (f)
and supported (s) degrees of freedom, and the free partition solved:

    [Kff Kfs] {uf}   {Ff}
    [Ksf Kss] {us} = {Fs}      with  us = 0

    Kff uf = Ff            ->  uf
    R = Ksf uf - Fs            ->  support reactions

Member axial force follows from the element compatibility relation

    N = (EA/L) * {c}^T ({u2} - {u1})

with tension positive.  Three independent verification quantities are returned
with every solution: the global equilibrium residual, the reaction/applied
load balance, and the strain-energy identity

    1/2 {F}^T {u}  ==  1/2 sum_m  N_m^2 L_m / (E A)_m

Self-straining loads and elastic supports
-----------------------------------------
``solve_static`` also accepts three optional ingredients (all default off, in
which case every number is bit-identical to the plain solve above):

* **Initial (free) member strain** ``eps0`` -- e.g. ``alpha * dT`` for a
  uniform temperature change.  The member force becomes

      N = EA * (e/L - eps0),         e = {c}^T ({u2} - {u1})

  and the strain is carried to the nodes as the equivalent (fixed-end) load

      {F0} = sum_m  EA_m eps0_m [ -c ; +c ]

  (a heated bar pushes its end nodes apart).  ``F0`` is self-equilibrated, so
  it changes no resultant and no reaction of a statically determinate truss.

* **Prescribed support displacements** ``{dof: value}`` on restrained DOFs
  (support settlement).  With ``us`` the prescribed values (0 where missing)

      Kff uf = Ff + F0f - Kfs us ,      R_s = Ksf uf + Kss us - Fs - F0s

* **Grounded translational springs** ``{dof: k}``.  Those DOFs stay *free*;
  ``k`` is added to the diagonal of ``Kff``.  The spring force on the node,
  ``-k (u - delta)``, is reported in ``reactions`` (and ``spring_forces``) so
  that applied loads + reactions still sum to zero.  A ``prescribed`` value
  on a spring DOF is a movement ``delta`` of the spring's ground end (an
  elastic support that settles); it adds ``k delta`` to the right-hand side.

With these terms the energy identity generalises (Clapeyron's theorem with
self-straining and support work).  Writing ``K_s`` for the bar stiffness only,
equilibrium at every DOF is ``K_s u = F + F0 + R`` with R = support reactions
+ spring forces.  Pre-multiplying by u^T and using u^T K_s u = sum EA/L e^2
and u^T F0 = sum EA eps0 e gives

    1/2 (F + R)^T u  ==  1/2 sum_m N_m e_m
                     ==  1/2 sum_m N_m^2 L_m/(EA)_m + 1/2 sum_m N_m eps0_m L_m

(the left side is ``strain_energy_ext``, the right ``strain_energy_int``).
Equivalently, with the springs inside K and no settlement,
1/2 u^T K u = sum 1/2 (EA/L) e^2 + sum 1/2 k u^2.  Because R at the
restrained DOFs is defined from K u, the identity tests the free-DOF
equations exactly as the plain one does: Ue - Ui = -1/2 uf^T (residual)f.
For a pure self-straining case both sides can be zero while the individual
terms are not (fixed-fixed heated bar: N e = 0), so ``energy_error`` is
measured relative to the largest single energy term (``energy_scale``).
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

try:
    import scipy.sparse as sp
    import scipy.sparse.linalg as spla
    _HAVE_SCIPY = True
except Exception:                                    # pragma: no cover
    _HAVE_SCIPY = False


# ---------------------------------------------------------------- results ---
@dataclass
class StaticResult:
    """Everything produced by one static solve."""

    name: str
    u: np.ndarray                 # (3n,) nodal displacements         [m]
    N: np.ndarray                 # (m,)  member axial forces, + = tension [N]
    stress: np.ndarray            # (m,)  axial stress                [Pa]
    reactions: np.ndarray         # (3n,) reaction vector (0 at free dofs) [N]
    F: np.ndarray                 # (3n,) applied load vector         [N]
    lengths: np.ndarray           # (m,)  member lengths              [m]

    # verification / diagnostics
    equilibrium_residual: float = 0.0
    reaction_balance: np.ndarray = field(default_factory=lambda: np.zeros(3))
    strain_energy_ext: float = 0.0
    strain_energy_int: float = 0.0
    solve_info: str = ""

    # self-straining / elastic-support extras (zeros for a plain solve)
    spring_forces: np.ndarray = field(default_factory=lambda: np.zeros(0))
    eps0: np.ndarray = None       # (m,) initial strain used, or None
    prescribed: dict = field(default_factory=dict)   # {dof: value} applied
    strain_energy_springs: float = 0.0   # 1/2 sum k (u - delta)^2   [J]
    energy_scale: float = 0.0     # largest energy term (see module doc)

    # ------------------------------------------------------------ queries --
    def displacement_vectors(self):
        return self.u.reshape(-1, 3)

    def displacement_magnitude(self):
        return np.linalg.norm(self.displacement_vectors(), axis=1)

    def max_displacement(self):
        d = self.displacement_magnitude()
        i = int(np.argmax(d))
        return i, float(d[i])

    def reaction_vectors(self):
        return self.reactions.reshape(-1, 3)

    def total_reaction(self):
        R = self.reaction_vectors()
        return R.sum(axis=0)

    def max_tension(self):
        i = int(np.argmax(self.N))
        return i, float(self.N[i])

    def max_compression(self):
        i = int(np.argmin(self.N))
        return i, float(self.N[i])

    @property
    def energy_error(self):
        denom = max(abs(self.strain_energy_ext), self.energy_scale, 1e-30)
        return abs(self.strain_energy_ext - self.strain_energy_int) / denom


# --------------------------------------------------------------- geometry ---
def direction_cosines(nodes, members):
    """Return ``(L, C)`` with lengths (m,) and unit direction cosines (m, 3)."""
    p = nodes[members[:, 0]]
    q = nodes[members[:, 1]]
    v = q - p
    L = np.linalg.norm(v, axis=1)
    if np.any(L <= 0.0):
        bad = int(np.argmin(L))
        raise ValueError(f"Member {bad} has zero length.")
    return L, v / L[:, None]


def element_stiffness(EA, L, c):
    """6x6 global element stiffness matrix for one bar."""
    cc = np.outer(c, c)
    k = np.empty((6, 6))
    k[0:3, 0:3] = cc
    k[0:3, 3:6] = -cc
    k[3:6, 0:3] = -cc
    k[3:6, 3:6] = cc
    return (EA / L) * k


def element_dofs(members):
    """(m, 6) array of global DOF numbers for every member."""
    i = members[:, 0]
    j = members[:, 1]
    return np.column_stack([3 * i, 3 * i + 1, 3 * i + 2,
                            3 * j, 3 * j + 1, 3 * j + 2])


# --------------------------------------------------------------- assembly ---
def assemble_stiffness(nodes, members, EA, L=None, C=None):
    """Assemble the global stiffness matrix (sparse CSR if SciPy is present)."""
    if L is None or C is None:
        L, C = direction_cosines(nodes, members)
    ndof = 3 * len(nodes)
    dofs = element_dofs(members)

    nm = len(members)
    rows = np.empty(36 * nm, dtype=np.int64)
    cols = np.empty(36 * nm, dtype=np.int64)
    vals = np.empty(36 * nm, dtype=float)

    for m in range(nm):
        ke = element_stiffness(EA[m], L[m], C[m])
        d = dofs[m]
        s = 36 * m
        rows[s:s + 36] = np.repeat(d, 6)
        cols[s:s + 36] = np.tile(d, 6)
        vals[s:s + 36] = ke.ravel()

    if _HAVE_SCIPY:
        K = sp.coo_matrix((vals, (rows, cols)), shape=(ndof, ndof)).tocsr()
    else:                                            # pragma: no cover
        K = np.zeros((ndof, ndof))
        np.add.at(K, (rows, cols), vals)
    return K


def assemble_geometric_stiffness(nodes, members, N, L=None, C=None):
    """Global geometric (initial stress) stiffness matrix.

    For a bar carrying axial force ``N`` the consistent geometric stiffness is

              N  |  [I - c c^T]   -[I - c c^T] |
        kg =  -  |                             |
              L  | -[I - c c^T]    [I - c c^T] |

    which captures the transverse restoring/destabilising effect of the axial
    force and is what makes the linearised buckling eigenproblem possible.
    """
    if L is None or C is None:
        L, C = direction_cosines(nodes, members)
    ndof = 3 * len(nodes)
    dofs = element_dofs(members)
    I3 = np.eye(3)

    nm = len(members)
    rows = np.empty(36 * nm, dtype=np.int64)
    cols = np.empty(36 * nm, dtype=np.int64)
    vals = np.empty(36 * nm, dtype=float)

    for m in range(nm):
        c = C[m]
        M = I3 - np.outer(c, c)
        kg = np.empty((6, 6))
        kg[0:3, 0:3] = M
        kg[0:3, 3:6] = -M
        kg[3:6, 0:3] = -M
        kg[3:6, 3:6] = M
        kg *= (N[m] / L[m])
        d = dofs[m]
        s = 36 * m
        rows[s:s + 36] = np.repeat(d, 6)
        cols[s:s + 36] = np.tile(d, 6)
        vals[s:s + 36] = kg.ravel()

    if _HAVE_SCIPY:
        return sp.coo_matrix((vals, (rows, cols)), shape=(ndof, ndof)).tocsr()
    Kg = np.zeros((ndof, ndof))                      # pragma: no cover
    np.add.at(Kg, (rows, cols), vals)
    return Kg


# ------------------------------------------------------------------ solve ---
def free_dofs(ndof, fixed):
    mask = np.ones(ndof, dtype=bool)
    mask[np.asarray(fixed, dtype=int)] = False
    return np.where(mask)[0]


_AXIS = ("ux", "uy", "uz")


def _dof_name(d):
    return f"DOF {d} (node {d // 3}, {_AXIS[d % 3]})"


def initial_strain_load(nodes, members, EA, eps0, L=None, C=None):
    """Equivalent nodal load of member initial strains.

        {F0} = sum_m EA_m eps0_m [ -c_m ; +c_m ]

    ``eps0`` may be a scalar or an (m,) array.  The vector is self-equilibrated
    (every member contributes two equal and opposite forces).
    """
    nodes = np.asarray(nodes, dtype=float)
    members = np.asarray(members, dtype=int)
    if L is None or C is None:
        L, C = direction_cosines(nodes, members)
    e0 = np.broadcast_to(np.asarray(eps0, dtype=float), (len(members),))
    f = (np.asarray(EA, dtype=float) * e0)[:, None] * C      # (m, 3)
    F0 = np.zeros(3 * len(nodes))
    d = element_dofs(members)
    np.add.at(F0, d[:, 0:3], -f)
    np.add.at(F0, d[:, 3:6], f)
    return F0


def solve_static(nodes, members, EA, F, fixed, name="LC", prescribed=None,
                 eps0=None, springs=None):
    """Solve one static load case and return a :class:`StaticResult`.

    Parameters
    ----------
    nodes, members, EA, F, fixed, name
        As before: coordinates (n, 3), connectivity (m, 2), axial rigidity
        (m,), global load vector (3n,) and the restrained DOF numbers.
    prescribed : dict {dof: value}, optional
        Prescribed displacement [m] of restrained DOFs (support settlement);
        restrained DOFs not listed stay at zero.  A value on a spring DOF is
        a movement of the spring's ground end.  A non-zero value on a DOF that
        is neither restrained nor a spring raises ``ValueError``.
    eps0 : float or (m,) array, optional
        Initial (free) member strain, e.g. ``alpha * dT``; N = EA (e/L - eps0).
    springs : dict {dof: k}, optional
        Grounded translational springs [N/m] on free DOFs.  A spring on a
        restrained DOF is ignored (the restraint governs) and noted in
        ``solve_info``.

    Equations and the generalised energy identity: see the module docstring.
    With all three left at ``None`` the arithmetic is exactly the plain solve.
    """
    nodes = np.asarray(nodes, dtype=float)
    members = np.asarray(members, dtype=int)
    EA = np.asarray(EA, dtype=float)
    F = np.asarray(F, dtype=float)

    ndof = 3 * len(nodes)
    L, C = direction_cosines(nodes, members)
    K = assemble_stiffness(nodes, members, EA, L, C)

    fd = free_dofs(ndof, fixed)
    sd = np.asarray(sorted(set(range(ndof)) - set(fd.tolist())), dtype=int)
    info = []

    # ---- optional ingredients ---------------------------------------------
    fixed_set = set(sd.tolist())
    spr = {}
    for dof, k in (springs or {}).items():
        dof = int(dof)
        k = float(k)
        if dof < 0 or dof >= ndof:
            raise ValueError(f"Spring at DOF {dof} is outside 0..{ndof - 1}.")
        if not (np.isfinite(k) and k > 0.0):
            raise ValueError(f"Spring stiffness at {_dof_name(dof)} must be a "
                             f"positive finite number (got {k!r}).")
        if dof in fixed_set:
            info.append(f"Spring at restrained {_dof_name(dof)} ignored.")
            continue
        spr[dof] = spr.get(dof, 0.0) + k

    us_map, delta = {}, {}
    for dof, val in (prescribed or {}).items():
        dof = int(dof)
        val = float(val)
        if not np.isfinite(val):
            raise ValueError(f"Prescribed displacement at {_dof_name(dof)} "
                             f"is not finite.")
        if dof < 0 or dof >= ndof:
            raise ValueError(f"Prescribed DOF {dof} is outside 0..{ndof - 1}.")
        if val == 0.0:
            continue
        if dof in fixed_set:
            us_map[dof] = val
        elif dof in spr:
            delta[dof] = val
        else:
            raise ValueError(
                f"Prescribed displacement {val:g} m at {_dof_name(dof)}, which "
                f"is neither restrained nor a spring support.  Settlements "
                f"apply to supported directions only.")

    e0 = None
    if eps0 is not None:
        e0 = np.array(np.broadcast_to(np.asarray(eps0, dtype=float),
                                      (len(members),)))
        if not np.all(np.isfinite(e0)):
            raise ValueError("Initial strains must be finite.")
        if not np.any(e0 != 0.0):
            e0 = None

    F0 = initial_strain_load(nodes, members, EA, e0, L, C) \
        if e0 is not None else None

    u = np.zeros(ndof)
    for dof, val in us_map.items():
        u[dof] = val
    spr_dofs = np.asarray(sorted(spr), dtype=int)
    spr_k = np.asarray([spr[d] for d in spr_dofs], dtype=float)
    spr_delta = np.asarray([delta.get(int(d), 0.0) for d in spr_dofs])
    pos = {int(d): i for i, d in enumerate(fd)}

    # generalised load scale for the residual norm (0 for a plain solve)
    extra_load_norms = []
    if F0 is not None:
        extra_load_norms.append(float(np.linalg.norm(F0)))

    if len(fd) == 0:
        info.insert(0, "Every degree of freedom is restrained.")
    else:
        rhs = F[fd]
        if F0 is not None:
            rhs = rhs + F0[fd]
        if us_map:
            Ku_s = K @ u if _HAVE_SCIPY else K.dot(u)   # only us is non-zero
            rhs = rhs - Ku_s[fd]
            extra_load_norms.append(float(np.linalg.norm(Ku_s)))
        if delta:
            kd = np.zeros(ndof)
            for d, k, dl in zip(spr_dofs, spr_k, spr_delta):
                kd[d] = k * dl
            rhs = rhs + kd[fd]
            extra_load_norms.append(float(np.linalg.norm(kd)))
        if _HAVE_SCIPY:
            Kff = K[fd][:, fd].tocsc()
            if len(spr_dofs):
                ii = np.asarray([pos[int(d)] for d in spr_dofs], dtype=int)
                Kspr = sp.coo_matrix((spr_k, (ii, ii)),
                                     shape=Kff.shape).tocsc()
                Kff = (Kff + Kspr).tocsc()
            # Guard against an unstable mechanism: a zero diagonal means a node
            # has no stiffness in that direction at all.
            diag = Kff.diagonal()
            if np.any(np.abs(diag) < 1e-12 * max(np.abs(diag).max(), 1.0)):
                bad = fd[np.where(np.abs(diag) < 1e-12 * np.abs(diag).max())[0]]
                raise np.linalg.LinAlgError(
                    f"Unstable structure: zero stiffness at DOF(s) "
                    f"{bad[:8].tolist()} (node {bad[0] // 3}).")
            try:
                uf = spla.spsolve(Kff, rhs)
                if not np.all(np.isfinite(uf)):
                    raise RuntimeError("non-finite solution")
                info.insert(0, "Sparse LU (SuperLU).")
            except Exception as exc:                 # pragma: no cover
                info.insert(0, f"Sparse solve fell back to dense ({exc}).")
                uf = np.linalg.solve(Kff.toarray(), rhs)
        else:                                        # pragma: no cover
            Kff = K[np.ix_(fd, fd)].copy()
            for d, k in zip(spr_dofs, spr_k):
                Kff[pos[int(d)], pos[int(d)]] += k
            uf = np.linalg.solve(Kff, rhs)
            info.insert(0, "Dense LAPACK solve.")
        u[fd] = uf

    # ---- member forces ----------------------------------------------------
    d = element_dofs(members)
    du = u[d[:, 3:6]] - u[d[:, 0:3]]
    elong = np.einsum("ij,ij->i", C, du)
    N = EA / L * elong
    if e0 is not None:
        N = N - EA * e0

    # ---- reactions --------------------------------------------------------
    Ku = K @ u if _HAVE_SCIPY else K.dot(u)
    reactions = np.zeros(ndof)
    if F0 is None:
        reactions[sd] = Ku[sd] - F[sd]
    else:
        reactions[sd] = Ku[sd] - F[sd] - F0[sd]
    spring_forces = np.zeros(ndof)
    if len(spr_dofs):
        spring_forces[spr_dofs] = -spr_k * (u[spr_dofs] - spr_delta)
        reactions[spr_dofs] += spring_forces[spr_dofs]

    # ---- verification quantities -----------------------------------------
    residual = Ku - F - reactions
    if F0 is not None:
        residual = residual - F0
    rnorm = float(np.linalg.norm(residual)
                  / max(np.linalg.norm(F), *extra_load_norms, 1e-12))

    applied = F.reshape(-1, 3).sum(axis=0)
    reac = reactions.reshape(-1, 3).sum(axis=0)
    balance = applied + reac

    Ue = 0.5 * float(F @ u)
    Ui = 0.5 * float(np.sum(N ** 2 * L / EA))
    Us = 0.0
    scale = 0.0
    if e0 is not None or us_map or len(spr_dofs):
        UR = 0.5 * float(reactions @ u)
        Ue = Ue + UR
        Ui_th = 0.5 * float(np.sum(N * e0 * L)) if e0 is not None else 0.0
        if len(spr_dofs):
            Us = 0.5 * float(np.sum(spr_k * (u[spr_dofs] - spr_delta) ** 2))
        # Reference energies of the imposed actions, so that a case whose
        # true work is zero (a determinate truss turning on a settled
        # support, a heated bar between walls) is still measured against a
        # physical magnitude rather than round-off:
        #   fully restrained thermal strain   1/2 sum EA L eps0^2
        #   prescribed motion, rest held      1/2 us^T K us
        #   spring ground movement            1/2 sum k delta^2
        U_th = (0.5 * float(np.sum(EA * L * e0 ** 2))
                if e0 is not None else 0.0)
        U_us = 0.0
        if us_map:
            us_vec = np.zeros(ndof)
            for dof, val in us_map.items():
                us_vec[dof] = val
            U_us = 0.5 * abs(float(us_vec @ (K @ us_vec if _HAVE_SCIPY
                                             else K.dot(us_vec))))
        U_dl = 0.5 * float(np.sum(spr_k * spr_delta ** 2))
        scale = max(abs(0.5 * float(F @ u)), abs(UR), Ui, abs(Ui_th), Us,
                    U_th, U_us, U_dl)
        Ui = Ui + Ui_th

    return StaticResult(
        name=name, u=u, N=N, stress=np.zeros_like(N), reactions=reactions,
        F=F, lengths=L, equilibrium_residual=rnorm, reaction_balance=balance,
        strain_energy_ext=Ue, strain_energy_int=Ui, solve_info=" ".join(info),
        spring_forces=spring_forces, eps0=e0,
        prescribed={**us_map, **delta}, strain_energy_springs=Us,
        energy_scale=scale)


def model_self_straining(model, T_factor=1.0):
    """Initial strains and prescribed displacements of the T load.

    Returns ``(eps0, prescribed)`` for ``T_factor`` times the model's
    self-straining load: ``eps0 = T_factor * alpha * dT`` (an (m,) array, or
    ``None`` when zero) and ``prescribed = {dof: T_factor * settlement}``
    built from ``model.loads.settlements`` (``None`` when empty).
    """
    if T_factor == 0.0:
        return None, None
    ls = model.loads
    eps0 = None
    dT = float(getattr(ls, "temperature_change", 0.0) or 0.0)
    if dT != 0.0:
        eps0 = np.full(model.geom.n_members,
                       T_factor * model.material.alpha * dT)
    prescribed = {}
    n = model.geom.n_nodes
    for k, s in enumerate(getattr(ls, "settlements", None) or []):
        nd = int(s["node"])
        if nd < 0 or nd >= n:
            raise ValueError(f"Settlement {k}: node {nd} is outside "
                             f"0..{n - 1}.")
        for a, key in enumerate(("ux", "uy", "uz")):
            v = float(s.get(key, 0.0) or 0.0)
            if v != 0.0:
                dof = 3 * nd + a
                prescribed[dof] = prescribed.get(dof, 0.0) + T_factor * v
    return eps0, (prescribed or None)


def solve_model(model, F, name="LC", T_factor=0.0):
    """Solve ``model`` for a given global load vector.

    ``T_factor`` scales the model's self-straining load T (uniform
    temperature change on every strut and the support settlements in
    ``model.loads``).  Supports and springs come from
    ``model.support_conditions()``.
    """
    fixed, springs = model.support_conditions()
    eps0, prescribed = model_self_straining(model, float(T_factor))
    res = solve_static(model.geom.nodes, model.geom.members, model.EA(),
                       F, fixed, name=name, prescribed=prescribed,
                       eps0=eps0, springs=springs or None)
    res.stress = res.N / model.member_areas()
    return res
