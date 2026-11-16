"""
Natural frequencies and mode shapes of a pin-jointed space truss.

Undamped free vibration
-----------------------
For free vibration of the unrestrained ("free") degrees of freedom,
``u(t) = phi sin(omega t)``, the equation of motion ``M u'' + K u = 0`` gives
the generalised symmetric eigenproblem

    K phi_n = omega_n^2 M phi_n ,        n = 1, 2, ...                  (1)

* ``K`` is the bar stiffness of ``solver.assemble_stiffness`` (the same matrix
  as the static solve), with the grounded springs of ELASTIC / CUSTOM supports
  added to its diagonal.  Rigid restraints from ``model.support_conditions()``
  are removed (their rows and columns are dropped).
* ``f_n = omega_n / (2 pi)`` [Hz] and ``T_n = 1 / f_n`` [s].

Mass matrix
-----------
A pin-jointed node has three translational DOFs and no rotary inertia.  The
default **lumped** mass puts half of every strut's mass at each end node and
adds the non-structural masses at the nodes, the same in x, y and z:

    m_i = sum_{bars at i} rho A L / 2             struts (always present)
        + q_D A_trib,i / g                         cladding dead pressure
        + max(0, -F_z) / g  of "D" nodal loads     hung dead loads
        + s p_s A_proj,i / g                       fraction s of the snow
        + user nodal masses                                              (2)

    M = diag(m_0, m_0, m_0, m_1, m_1, m_1, ...)

``A_trib`` and ``A_proj`` are the surface and plan tributary areas used by
``loads`` (one third of each adjacent panel), so with self weight switched on
the mass is exactly the dead load divided by g:  ``sum m_i g = -sum F_z(D)``.
The strut mass is physical and is therefore always included, whatever
``LoadSpec.include_self_weight`` says about the static load.

The **consistent** bar mass (option ``mass_type="consistent"``) is

    m_e = (rho A L / 6) [ 2 I   I ]                                       (3)
                        [  I   2 I ]

in every direction (the same translation mass as (2); OpenSees ``-cMass 1``).
Lumped mass gives frequencies below, and consistent mass above, those of the
continuous bar; both converge as O(h^2) (``tests/test_modal.py`` shows it on a
fixed-free rod against the exact discrete and continuous solutions).

Solution of (1)
---------------
* Free DOFs with neither stiffness nor mass (nodes connected to nothing) are
  left out with a warning.
* Massless free DOFs (possible only with zero-density bars) are removed by
  exact static condensation (Guyan reduction is exact when those DOFs carry
  no mass):

      K* = K_mm - K_m0 K_00^-1 K_0m ,       phi_0 = -K_00^-1 K_0m phi_m

* Up to ``dense_limit`` (default 2400) dynamic DOFs: dense LAPACK
  generalised symmetric eigensolve (``scipy.linalg.eigh(K*, M*)``).
* Above it: shift-invert Lanczos (ARPACK ``eigsh`` with ``sigma = 0``, the
  modes nearest zero frequency), massless components recomputed from the
  condensation relation, then a Rayleigh-Ritz step on the computed subspace.
  A singular K (a mechanism) makes ``sigma = 0`` fail; a small negative shift
  is used instead.

Every mode is mass-normalised, ``phi_n^T M phi_n = 1``, so that
``phi_n^T K phi_n = omega_n^2`` (Rayleigh quotient) and ``Phi^T M Phi = I``.

Participation factors and effective modal mass  (Chopra 2023, Ch. 13)
---------------------------------------------------------------------
With ``iota_d`` the influence vector of a unit rigid translation in direction
d (1 at every free DOF in that direction, 0 elsewhere):

    L_nd     = phi_n^T M iota_d                      (excitation factor)
    Gamma_nd = L_nd / M_n = L_nd                     (M_n = 1)  [kg^0.5]
    M*_nd    = L_nd^2 / M_n = Gamma_nd^2                        [kg]
    sum over ALL modes of M*_nd = iota_d^T M iota_d                       (4)

because a complete set of M-orthonormal modes satisfies Phi Phi^T = M^-1.
The right side of (4) is the mass on the free DOFs in direction d
(``free_mass``); the effective-mass ratios are ``M*_nd / free_mass_d``, and
their running sum tells how many modes are needed to capture, say, the 90 %
of the mass that response-spectrum practice (ASCE 7, Ch. 12) asks for.  The
ground ends of spring supports move with the ground, so ``iota_d`` is a unit
translation for ELASTIC supports too.  With consistent mass the coupling of
free DOFs to the restrained ones (M_fs) is not part of ``iota_d^T M iota_d``.

Repeated frequencies
--------------------
A dome on a regular ring has n-fold symmetry about z, so many frequencies come
in exactly equal pairs.  Inside a repeated eigenvalue any rotation of the
modes is an equally valid basis, and the computed basis is arbitrary.  With
``align_degenerate=True`` each cluster of equal eigenvalues is rotated (QR of
its 3-column participation matrix) so that the first mode carries all of the
cluster's x participation, the next one none of the x and all of the remaining
y participation, and so on: a sway pair becomes one pure-x and one pure-y
mode.  The total effective mass of the cluster is unchanged.

Pre-stress (stress stiffening)
------------------------------
Optionally the geometric stiffness of a static state is added,

    [K + K_g(N)] phi = omega^2 M phi                                      (5)

with ``K_g`` from ``solver.assemble_geometric_stiffness`` (the matrix of the
linear buckling analysis).  Tension raises and compression lowers the
frequencies.  Since ``det(K + lambda_cr K_g) = 0`` defines the linear
buckling load, the lowest frequency of (5) is exactly zero at
``N = lambda_cr N_ref``, with the buckling mode as its mode shape (the
vibration criterion of stability); ``tests/test_modal.py`` checks this
identity against ``buckling.linear_buckling``.  ``frequency_load_curve``
traces omega^2 against the load factor of one combination.

Verification quantities returned with every solution
-----------------------------------------------------
``orthogonality_error = max |Phi^T M Phi - I|``,
``stiffness_error = max |Phi^T K Phi - diag(omega^2)| / max |omega^2|`` and
``residual = max_n ||K phi_n - omega_n^2 M phi_n|| /
((||K|| + |omega_n^2| ||M||) ||phi_n||)`` (normwise backward error).

References
----------
* A. K. Chopra (2023). Dynamics of Structures, 6th ed., Pearson.  Ch. 10
  (free vibration, orthogonality, normalisation), Ch. 13 (participation
  factor, effective modal mass).
* W. McGuire, R. H. Gallagher, R. D. Ziemian (2000). Matrix Structural
  Analysis, 2nd ed., Wiley.  Ch. 3 (direct stiffness assembly), Ch. 9
  (geometric stiffness).
* K. J. Bathe (2014). Finite Element Procedures, 2nd ed.  (lumped and
  consistent mass, static condensation, generalised eigensolution, cited at
  book level).
* R. R. Craig Jr., A. J. Kurdila (2006). Fundamentals of Structural Dynamics,
  2nd ed., Wiley.  (mass-normalised modes, Rayleigh-Ritz, axial vibration of
  a uniform rod; cited at book level).
"""

from __future__ import annotations

import csv
import math
from dataclasses import dataclass, field

import numpy as np

try:
    import scipy.linalg as sla
    import scipy.sparse as sp
    import scipy.sparse.linalg as spla
    _HAVE_SCIPY = True
except Exception:                                    # pragma: no cover
    _HAVE_SCIPY = False

from .solver import assemble_stiffness, assemble_geometric_stiffness, \
    direction_cosines, element_dofs, free_dofs
from .units import G_ACCEL

LUMPED = "lumped"
CONSISTENT = "consistent"
MASS_TYPES = (LUMPED, CONSISTENT)

_AXES = ("x", "y", "z")
_DEGENERATE_RTOL = 1e-8          # relative gap below which omega^2 are equal
_ZERO_RTOL = 1e-10               # |omega^2| below this x (max Kii/Mii) is 0
_DENSE_MAX = 12000               # never build a dense problem larger than this


def _require_scipy():
    if not _HAVE_SCIPY:                              # pragma: no cover
        raise ImportError("Modal analysis needs SciPy (scipy.linalg and "
                          "scipy.sparse.linalg).  Install scipy >= 1.7.")


# ===================================================================== mass ==
@dataclass
class MassModel:
    """Lumped nodal masses of a model, split by origin.

    ``nodal`` is the lumped mass of every node [kg], the same in x, y and z.
    ``parts`` maps each source ("struts", "cladding", "dead nodal loads",
    "snow", "added") to its (n,) contribution; the parts sum to ``nodal``.
    ``line_mass`` is rho*A of every strut [kg/m] (used for consistent mass).
    """

    nodal: np.ndarray
    parts: dict
    line_mass: np.ndarray
    notes: list = field(default_factory=list)

    @property
    def total(self):
        return float(self.nodal.sum())

    def totals(self):
        """{source: total mass [kg]} for every non-zero source."""
        return {k: float(v.sum()) for k, v in self.parts.items()
                if np.any(v != 0.0)}

    @property
    def non_structural(self):
        """(n,) nodal mass without the strut part [kg]."""
        out = np.zeros_like(self.nodal)
        for k, v in self.parts.items():
            if k != "struts":
                out = out + v
        return out


def _nodal_mass_input(nodal_masses, n):
    """Normalise the user ``nodal_masses`` argument to an (n,) array [kg].

    Accepts None, a dict {node: kg}, a list of {"node": i, "mass": kg} dicts,
    or an array of length n.
    """
    out = np.zeros(n)
    if nodal_masses is None:
        return out
    if isinstance(nodal_masses, dict):
        items = list(nodal_masses.items())
    elif (isinstance(nodal_masses, (list, tuple)) and nodal_masses
          and isinstance(nodal_masses[0], dict)):
        items = [(d["node"], d["mass"]) for d in nodal_masses]
    else:
        arr = np.asarray(nodal_masses, dtype=float).reshape(-1)
        if arr.shape != (n,):
            raise ValueError(f"nodal_masses must have one value per node "
                             f"({n}), got {arr.size}.")
        items = list(enumerate(arr))
    for node, mass in items:
        nd = int(node)
        if nd < 0 or nd >= n:
            raise ValueError(f"nodal_masses: node {node} is outside "
                             f"0..{n - 1}.")
        mv = float(mass)
        if not (np.isfinite(mv) and mv >= 0.0):
            raise ValueError(f"nodal_masses: the mass at node {nd} must be a "
                             f"finite number >= 0 kg (got {mass!r}).")
        out[nd] += mv
    return out


def lumped_masses(model, include_dead=True, nodal_masses=None,
                  snow_fraction=0.0):
    """Lumped nodal masses of ``model`` (equation (2) of the module doc).

    Parameters
    ----------
    include_dead : bool
        Add the superimposed dead load as mass: cladding pressure x surface
        tributary area / g, and the downward part of explicit nodal loads of
        case "D" / g.
    nodal_masses : dict {node: kg}, list of {"node", "mass"}, or (n,) array
        Extra point masses (equipment, hubs, lanterns).
    snow_fraction : float
        Fraction of the (balanced) snow load counted as mass, 0 by default.
        Seismic-weight rules (ASCE 7, Ch. 12) include part of the roof snow
        only for heavy snow; the choice is the engineer's.

    Returns a :class:`MassModel`.
    """
    from . import loads as _loads

    g = model.geom
    n = g.n_nodes
    notes = []
    rho = float(model.material.rho)
    A = model.member_areas()
    L = g.member_lengths()
    mu = rho * A                                     # kg/m of every strut
    struts = np.zeros(n)
    np.add.at(struts, g.members[:, 0], 0.5 * mu * L)
    np.add.at(struts, g.members[:, 1], 0.5 * mu * L)

    clad = np.zeros(n)
    dead_nodal = np.zeros(n)
    if include_dead:
        q = float(model.loads.dead_pressure)
        if q < 0.0:
            raise ValueError("A negative dead pressure cannot be a mass.")
        clad = q * np.asarray(g.trib_area, dtype=float) / G_ACCEL
        if getattr(model.loads, "nodal_loads", None):
            FD = _loads.nodal_load_vectors(model)["D"]
            Fz = FD[2::3]
            dead_nodal = np.maximum(-Fz, 0.0) / G_ACCEL
            if np.any(Fz > 0.0):
                notes.append("Upward 'D' nodal loads were not counted as "
                             "mass.")
            if np.any(FD[0::3] != 0.0) or np.any(FD[1::3] != 0.0):
                notes.append("Horizontal components of 'D' nodal loads are "
                             "not mass; only the downward part was used.")

    snow = np.zeros(n)
    s = float(snow_fraction)
    if not (np.isfinite(s) and 0.0 <= s <= 1.0):
        raise ValueError(f"snow_fraction must be between 0 and 1 (got "
                         f"{snow_fraction!r}).")
    if s > 0.0:
        ps = float(model.loads.snow_pressure)
        snow = s * max(ps, 0.0) * np.asarray(g.trib_area_proj,
                                             dtype=float) / G_ACCEL

    added = _nodal_mass_input(nodal_masses, n)
    parts = {"struts": struts, "cladding": clad,
             "dead nodal loads": dead_nodal, "snow": snow, "added": added}
    nodal = struts + clad + dead_nodal + snow + added
    return MassModel(nodal=nodal, parts=parts, line_mass=mu, notes=notes)


def _per_dof(masses, n):
    """(3n,) per-DOF mass from None, an (n,) or a (3n,) array."""
    if masses is None:
        return np.zeros(3 * n)
    m = np.asarray(masses, dtype=float).reshape(-1)
    if m.size == n:
        m = np.repeat(m, 3)
    elif m.size != 3 * n:
        raise ValueError(f"masses must have n = {n} (same in x, y, z) or "
                         f"3n = {3 * n} (per DOF) entries, got {m.size}.")
    if not np.all(np.isfinite(m)) or np.any(m < 0.0):
        raise ValueError("Nodal masses must be finite and >= 0.")
    return m


def mass_matrix(nodes, members, masses=None, line_mass=None,
                mass_type=LUMPED):
    """Global mass matrix (3n x 3n, sparse CSR).

    ``masses``: (n,) nodal masses (same in x, y, z) or (3n,) per-DOF masses,
    always lumped.  ``line_mass``: rho*A of every bar [kg/m] (scalar or (m,)),
    distributed lumped (rho A L / 2 per end and direction) or consistent
    (equation (3)).
    """
    _require_scipy()
    if mass_type not in MASS_TYPES:
        raise ValueError(f"mass_type must be one of {MASS_TYPES}, got "
                         f"{mass_type!r}.")
    nodes = np.asarray(nodes, dtype=float)
    members = np.asarray(members, dtype=int)
    n = len(nodes)
    ndof = 3 * n
    diag = _per_dof(masses, n)
    rows = [np.arange(ndof)]
    cols = [np.arange(ndof)]
    vals = [diag]
    if line_mass is not None and len(members):
        mu = np.array(np.broadcast_to(np.asarray(line_mass, dtype=float),
                                      (len(members),)))
        if not np.all(np.isfinite(mu)) or np.any(mu < 0.0):
            raise ValueError("line_mass (rho*A) must be finite and >= 0.")
        L, _ = direction_cosines(nodes, members)
        mb = mu * L                                  # mass of every bar
        d = element_dofs(members)
        for a in range(3):
            i, j = d[:, a], d[:, 3 + a]
            if mass_type == LUMPED:
                rows += [i, j]
                cols += [i, j]
                vals += [0.5 * mb, 0.5 * mb]
            else:
                c = mb / 6.0
                rows += [i, j, i, j]
                cols += [i, j, j, i]
                vals += [2.0 * c, 2.0 * c, c, c]
    M = sp.coo_matrix((np.concatenate(vals),
                       (np.concatenate(rows), np.concatenate(cols))),
                      shape=(ndof, ndof)).tocsr()
    return M


# ================================================================== result ==
@dataclass
class ModalResult:
    """Everything produced by one modal analysis.

    Arrays over modes have length k (the number of modes returned), in
    ascending order of omega^2.  ``modes`` holds full-DOF, mass-normalised
    vectors (zeros at restrained DOFs), DOF = 3*node + {0: x, 1: y, 2: z}.
    """

    omega2: np.ndarray            # (k,) eigenvalues omega^2      [rad^2/s^2]
    omega: np.ndarray             # (k,) circular frequency       [rad/s]
    frequencies: np.ndarray       # (k,) natural frequency        [Hz]
    periods: np.ndarray           # (k,) period (inf if f = 0)    [s]
    modes: np.ndarray             # (3n, k) phi^T M phi = 1    [kg^-0.5]
    participation: np.ndarray     # (k, 3) Gamma = L = phi^T M iota [kg^0.5]
    effective_mass: np.ndarray    # (k, 3) Gamma^2                [kg]
    effective_mass_ratio: np.ndarray   # (k, 3) M* / free_mass      [-]
    dof_mass: np.ndarray          # (3n,) lumped-equivalent mass per DOF [kg]
    total_mass: np.ndarray        # (3,) mass of the whole model, x/y/z [kg]
    free_mass: np.ndarray         # (3,) mass on the free DOFs, x/y/z [kg]
    free_dofs: np.ndarray         # free DOF numbers in the eigenproblem
    mass_type: str = LUMPED
    n_free: int = 0               # free DOFs in the eigenproblem
    n_dynamic: int = 0            # free DOFs with mass (= possible modes)
    method: str = ""
    info: str = ""
    warnings: list = field(default_factory=list)
    notes: list = field(default_factory=list)
    mass_parts: dict = field(default_factory=dict)   # {source: kg}
    prestress: str = ""           # description of K_g, "" if none
    orthogonality_error: float = 0.0
    stiffness_error: float = 0.0
    residual: float = 0.0

    # ------------------------------------------------------------ queries --
    @property
    def n_modes(self):
        return len(self.omega2)

    @property
    def nodal_mass(self):
        """(n, 3) lumped-equivalent mass of every node in x, y, z [kg]."""
        return self.dof_mass.reshape(-1, 3)

    @property
    def cumulative_ratio(self):
        """(k, 3) running sum of the effective-mass ratios."""
        return np.cumsum(self.effective_mass_ratio, axis=0)

    @property
    def total_participating(self):
        """(3,) effective mass captured by the returned modes [kg]."""
        return self.effective_mass.sum(axis=0)

    @property
    def participating_ratio(self):
        """(3,) fraction of the free mass captured by the returned modes."""
        return self.effective_mass_ratio.sum(axis=0)

    @property
    def fundamental_frequency(self):
        f = self.frequencies
        return float(f[0]) if len(f) else float("nan")

    @property
    def unstable(self):
        """(k,) True where omega^2 < 0 (pre-stress beyond stability)."""
        return ~np.isfinite(self.omega)

    def mode(self, k):
        """(n, 3) mass-normalised displacement of mode ``k`` (0-based)."""
        return self.modes[:, k].reshape(-1, 3)

    def mode_shape(self, k):
        """(n, 3) mode ``k`` scaled so that the largest component is 1."""
        v = self.modes[:, k]
        a = np.max(np.abs(v))
        return (v / (a if a > 0.0 else 1.0)).reshape(-1, 3)

    def dominant_direction(self, k, threshold=0.01):
        """'x', 'y' or 'z' with the largest effective-mass ratio of mode k,
        or '-' when no direction reaches ``threshold``."""
        r = self.effective_mass_ratio[k]
        a = int(np.argmax(r))
        return _AXES[a] if r[a] >= threshold else "-"

    def modes_for_ratio(self, target=0.90):
        """{'x','y','z': number of modes whose cumulative effective mass
        reaches ``target``} (None when the returned modes do not reach it,
        or when there is no free mass in that direction)."""
        cum = self.cumulative_ratio
        out = {}
        for a, ax in enumerate(_AXES):
            if not self.free_mass[a] > 0.0 or not len(cum):
                out[ax] = None
                continue
            hit = np.where(cum[:, a] >= target - 1e-12)[0]
            out[ax] = int(hit[0]) + 1 if len(hit) else None
        return out

    def table(self):
        """One dict per mode (for GUI tables, CSV and reports)."""
        rows = []
        cum = self.cumulative_ratio
        for k in range(self.n_modes):
            rows.append({
                "mode": k + 1,
                "frequency_Hz": float(self.frequencies[k]),
                "period_s": float(self.periods[k]),
                "omega_rad_s": float(self.omega[k]),
                "omega2": float(self.omega2[k]),
                "Gamma_x": float(self.participation[k, 0]),
                "Gamma_y": float(self.participation[k, 1]),
                "Gamma_z": float(self.participation[k, 2]),
                "Meff_x_kg": float(self.effective_mass[k, 0]),
                "Meff_y_kg": float(self.effective_mass[k, 1]),
                "Meff_z_kg": float(self.effective_mass[k, 2]),
                "ratio_x": float(self.effective_mass_ratio[k, 0]),
                "ratio_y": float(self.effective_mass_ratio[k, 1]),
                "ratio_z": float(self.effective_mass_ratio[k, 2]),
                "cum_x": float(cum[k, 0]),
                "cum_y": float(cum[k, 1]),
                "cum_z": float(cum[k, 2]),
                "dominant": self.dominant_direction(k),
            })
        return rows

    def summary_lines(self):
        tm = self.total_mass
        same = np.allclose(tm, tm[0], rtol=1e-12, atol=0.0)
        lines = [
            f"Modal analysis         {self.mass_type} mass, {self.n_modes} "
            f"mode(s), {self.n_free} free DOF ({self.method})",
            ("Total mass             " + (f"{tm[0]:,.1f} kg" if same else
             " / ".join(f"{v:,.1f}" for v in tm) + " kg (x/y/z)")),
        ]
        if self.mass_parts:
            parts = ", ".join(f"{k} {v:,.1f}" for k, v in
                              self.mass_parts.items())
            lines.append(f"Mass sources           {parts} kg")
        lines.append("Free (participating)   "
                     + " / ".join(f"{v:,.1f}" for v in self.free_mass)
                     + " kg (x/y/z)")
        if self.prestress:
            lines.append(f"Pre-stress (K_g)       {self.prestress}")
        pr = self.participating_ratio * 100.0
        lines.append(f"Captured by the modes  x {pr[0]:.1f} %, "
                     f"y {pr[1]:.1f} %, z {pr[2]:.1f} %")
        need = self.modes_for_ratio(0.90)
        more = f"> {self.n_modes}"
        lines.append("Modes for 90 % mass    " + ", ".join(
            f"{ax} {need[ax] if need[ax] is not None else more}"
            for ax in _AXES))
        lines.append(f"Checks                 |Phi'M Phi - I| = "
                     f"{self.orthogonality_error:.1e}, |Phi'K Phi - W| = "
                     f"{self.stiffness_error:.1e}, residual = "
                     f"{self.residual:.1e}")
        lines.append("")
        lines.append(f"{'mode':>4s} {'f [Hz]':>10s} {'T [s]':>9s} "
                     f"{'Mx %':>7s} {'My %':>7s} {'Mz %':>7s} "
                     f"{'sum x':>7s} {'sum y':>7s} {'sum z':>7s}  dir")
        for r in self.table():
            T = r["period_s"]
            Ts = "inf" if math.isinf(T) else f"{T:9.4f}"
            lines.append(
                f"{r['mode']:>4d} {r['frequency_Hz']:>10.4f} {Ts:>9s} "
                f"{100 * r['ratio_x']:>7.2f} {100 * r['ratio_y']:>7.2f} "
                f"{100 * r['ratio_z']:>7.2f} {100 * r['cum_x']:>7.2f} "
                f"{100 * r['cum_y']:>7.2f} {100 * r['cum_z']:>7.2f}  "
                f"{r['dominant']}")
        for w in self.warnings:
            lines.append("WARNING: " + w)
        for n_ in self.notes:
            lines.append("Note: " + n_)
        return lines

    def summary_text(self):
        return "\n".join(self.summary_lines())


def export_modal_csv(result, path):
    """Write ``result.table()`` to a CSV file and return the path."""
    rows = result.table()
    keys = list(rows[0].keys()) if rows else ["mode"]
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow(r)
    return path


# ================================================================== solver ==
def _clusters(w, scale):
    """Index ranges [s, e) of consecutive (nearly) equal eigenvalues."""
    out = []
    s = 0
    for i in range(1, len(w) + 1):
        if i == len(w) or abs(w[i] - w[i - 1]) > (
                _DEGENERATE_RTOL * max(abs(w[i]), abs(w[i - 1]))
                + 1e-12 * scale):
            out.append((s, i))
            s = i
    return out


def _condensation(Kd, ml, dyn):
    """-K_00^-1 K_0m for the massless DOFs ``ml`` (dense), or raise."""
    K00 = Kd[np.ix_(ml, ml)]
    K0m = Kd[np.ix_(ml, dyn)]
    lu, piv = sla.lu_factor(K00)
    dU = np.abs(np.diag(lu))
    if dU.size and (dU.min() <= 1e-13 * max(dU.max(), 1e-300)):
        raise ValueError("The massless degrees of freedom form a mechanism "
                         "(their stiffness block is singular); give those "
                         "nodes some mass or restrain them.")
    return -sla.lu_solve((lu, piv), K0m)


def solve_modes(nodes, members, EA, masses, fixed, line_mass=None,
                springs=None, n_modes=12, mass_type=LUMPED, N_prestress=None,
                dense_limit=2400, align_degenerate=True):
    """Natural frequencies and modes of a space truss (low-level API).

    Parameters
    ----------
    nodes, members, EA, fixed
        As for :func:`solver.solve_static`: coordinates (n, 3) [m],
        connectivity (m, 2), axial rigidity (m,) [N] and restrained DOFs.
    masses : None, (n,) or (3n,) array
        Lumped nodal masses [kg] (same in x, y, z, or per DOF).
    line_mass : None, float or (m,) array
        rho*A of the bars [kg/m], lumped or consistent per ``mass_type``.
    springs : dict {dof: k}, optional
        Grounded springs [N/m] on free DOFs (springs on restrained DOFs are
        ignored, as in the static solver).
    n_modes : int or None
        Number of modes to return (None = all).
    N_prestress : (m,) array, optional
        Member axial forces [N] whose geometric stiffness is added (5).
    dense_limit : int
        Largest number of dynamic DOFs solved densely.
    align_degenerate : bool
        Rotate repeated-frequency clusters to align with x, y, z.

    Returns a :class:`ModalResult`.  Equations: see the module docstring.
    """
    _require_scipy()
    if mass_type not in MASS_TYPES:
        raise ValueError(f"mass_type must be one of {MASS_TYPES}, got "
                         f"{mass_type!r}.")
    nodes = np.asarray(nodes, dtype=float)
    members = np.asarray(members, dtype=int)
    EA = np.asarray(EA, dtype=float)
    n = len(nodes)
    ndof = 3 * n
    if n_modes is not None:
        if isinstance(n_modes, str):
            if n_modes.strip().lower() != "all":
                raise ValueError("n_modes must be a positive integer, None "
                                 "or 'all'.")
            n_modes = None
        elif int(n_modes) < 1:
            raise ValueError("n_modes must be at least 1.")

    warnings, notes = [], []
    L, C = direction_cosines(nodes, members)
    K = assemble_stiffness(nodes, members, EA, L, C)
    pre = False
    if N_prestress is not None:
        Np = np.asarray(N_prestress, dtype=float).reshape(-1)
        if Np.shape != (len(members),):
            raise ValueError(f"N_prestress needs one force per member "
                             f"({len(members)}), got {Np.size}.")
        if not np.all(np.isfinite(Np)):
            raise ValueError("N_prestress must be finite.")
        if np.any(Np != 0.0):
            K = (K + assemble_geometric_stiffness(nodes, members, Np, L,
                                                  C)).tocsr()
            pre = True

    fixed = np.unique(np.asarray(fixed, dtype=int).reshape(-1))
    if fixed.size and (fixed.min() < 0 or fixed.max() >= ndof):
        raise ValueError(f"A restrained DOF is outside 0..{ndof - 1}.")
    fixed_set = set(fixed.tolist())
    sd, kv = [], []
    for dof, k in (springs or {}).items():
        dof, k = int(dof), float(k)
        if dof < 0 or dof >= ndof:
            raise ValueError(f"Spring at DOF {dof} is outside 0..{ndof - 1}.")
        if not (np.isfinite(k) and k > 0.0):
            raise ValueError(f"Spring stiffness at DOF {dof} must be a "
                             f"positive finite number (got {k!r}).")
        if dof in fixed_set:
            notes.append(f"Spring at restrained DOF {dof} ignored.")
            continue
        sd.append(dof)
        kv.append(k)
    if sd:
        K = (K + sp.coo_matrix((kv, (sd, sd)), shape=K.shape)).tocsr()

    M = mass_matrix(nodes, members, masses, line_mass, mass_type)
    dof_mass = np.asarray(M.sum(axis=1)).reshape(-1)
    total_mass = dof_mass.reshape(-1, 3).sum(axis=0)

    fd = free_dofs(ndof, fixed)
    Kf = K[fd][:, fd].tocsr()
    Mf = M[fd][:, fd].tocsr()

    # ---- DOFs with neither stiffness nor mass are not part of the problem --
    krow = np.asarray(abs(Kf).sum(axis=1)).reshape(-1)
    mrow = np.asarray(abs(Mf).sum(axis=1)).reshape(-1)
    loose = (krow == 0.0) & (mrow == 0.0)
    if np.any(loose):
        bad = fd[loose]
        warnings.append(
            f"{int(loose.sum())} free DOF(s) have neither stiffness nor mass "
            f"and were left out (first: DOF {int(bad[0])}, node "
            f"{int(bad[0]) // 3}); restrain them or connect those nodes.")
        keep = ~loose
        fd = fd[keep]
        Kf = Kf[keep][:, keep].tocsr()
        Mf = Mf[keep][:, keep].tocsr()
        mrow = mrow[keep]
    nf = len(fd)
    dyn = np.where(mrow > 0.0)[0]
    ml = np.where(mrow == 0.0)[0]
    nd = len(dyn)
    if nd == 0:
        raise ValueError("No free degree of freedom carries mass; there is "
                         "nothing to vibrate.")

    # influence vectors (rigid unit translation of the free DOFs) and M iota
    R = np.zeros((nf, 3))
    for a in range(3):
        iota = (fd % 3 == a).astype(float)
        R[:, a] = Mf @ iota
    free_mass = np.array([float(((fd % 3 == a).astype(float)) @ R[:, a])
                          for a in range(3)])

    kdiag = Kf.diagonal()[dyn]
    mdiag = Mf.diagonal()[dyn]
    lam_scale = float(np.max(np.abs(kdiag) / np.maximum(mdiag, 1e-300))) \
        if nd else 1.0
    lam_scale = max(lam_scale, 1e-300)

    want = nd if n_modes is None else min(int(n_modes), nd)
    extra = 4 if align_degenerate else 0
    use_dense = nd <= int(dense_limit) or want >= nd - 1
    if use_dense and nd > _DENSE_MAX:
        raise ValueError(f"{want} of {nd} modes would need a dense "
                         f"{nd} x {nd} eigensolve; ask for fewer modes (the "
                         f"sparse solver is used for n_modes < {nd - 1}) or "
                         f"raise dense_limit deliberately.")

    method = ""
    if not use_dense:
        kk = min(want + extra, nd - 1)
        try:
            V, w, method = _sparse_modes(Kf, Mf, dyn, ml, kk, lam_scale)
        except Exception as exc:                     # pragma: no cover
            if nd > _DENSE_MAX:
                raise
            warnings.append(f"Sparse eigensolver failed ({exc}); solved "
                            f"densely instead.")
            use_dense = True
    if use_dense:
        kk = min(want + extra, nd)
        Kd = Kf.toarray()
        Md = Mf.toarray()
        if len(ml):
            Tm = _condensation(Kd, ml, dyn)
            Ks = Kd[np.ix_(dyn, dyn)] + Kd[np.ix_(dyn, ml)] @ Tm
        else:
            Ks = Kd
        Ms = Md[np.ix_(dyn, dyn)]
        Ks = 0.5 * (Ks + Ks.T)
        Ms = 0.5 * (Ms + Ms.T)
        if kk < nd:
            w, Vd = sla.eigh(Ks, Ms, subset_by_index=[0, kk - 1])
        else:
            w, Vd = sla.eigh(Ks, Ms)
        V = np.zeros((nf, len(w)))
        V[dyn] = Vd
        if len(ml):
            V[ml] = Tm @ Vd
        method = (f"dense LAPACK generalised eigensolve, {nd} dynamic DOF"
                  + (f", {len(ml)} massless DOF condensed" if len(ml) else ""))

    order = np.argsort(w)
    w, V = np.asarray(w)[order], V[:, order]

    # ---- re-normalise (phi' M phi = 1) -------------------------------------
    nrm = np.sqrt(np.einsum("ij,ij->j", V, Mf @ V))
    V = V / np.where(nrm > 0.0, nrm, 1.0)

    # ---- align repeated-frequency clusters with x, y, z ---------------------
    sqrt_m = math.sqrt(max(float(free_mass.max()), 1e-300))
    if align_degenerate:
        for s, e in _clusters(w, lam_scale):
            if e - s < 2:
                continue
            Pc = V[:, s:e]
            Lc = Pc.T @ R                            # (p, 3)
            # directions in which the cluster participates at all (a column
            # of round-off would otherwise pick an arbitrary first vector)
            sig = [a for a in range(3)
                   if np.linalg.norm(Lc[:, a]) > 1e-8 * sqrt_m]
            if not sig:
                continue
            Q, _ = np.linalg.qr(Lc[:, sig], mode="complete")
            V[:, s:e] = Pc @ Q
            # Rayleigh quotients of the rotated vectors
            KV = Kf @ V[:, s:e]
            w[s:e] = np.einsum("ij,ij->j", V[:, s:e], KV)

    V, w = V[:, :want], w[:want]

    # ---- sign convention: dominant participation (else largest entry) > 0 --
    Lp = V.T @ R
    for i in range(want):
        a = int(np.argmax(np.abs(Lp[i])))
        if abs(Lp[i, a]) > 1e-8 * sqrt_m:
            sgn = np.sign(Lp[i, a])
        else:
            sgn = np.sign(V[int(np.argmax(np.abs(V[:, i]))), i])
        if sgn < 0:
            V[:, i] = -V[:, i]
            Lp[i] = -Lp[i]

    # ---- participation and effective mass ----------------------------------
    eff = Lp ** 2
    ratio = np.zeros_like(eff)
    for a in range(3):
        if free_mass[a] > 0.0:
            ratio[:, a] = eff[:, a] / free_mass[a]

    # ---- verification quantities -------------------------------------------
    KV = Kf @ V
    MV = Mf @ V
    G = V.T @ MV
    orth = float(np.max(np.abs(G - np.eye(want)))) if want else 0.0
    H = V.T @ KV
    # relative to the largest |omega^2| returned, but never to less than the
    # zero-frequency threshold (all-mechanism results would divide by ~0)
    wmax = max(float(np.max(np.abs(w))) if want else 0.0,
               _ZERO_RTOL * lam_scale, 1e-300)
    stiff = float(np.max(np.abs(H - np.diag(w)))) / wmax if want else 0.0
    kn = float(spla.norm(Kf, 1))
    mn = float(spla.norm(Mf, 1))
    res = 0.0
    for i in range(want):
        r = KV[:, i] - w[i] * MV[:, i]
        den = (kn + abs(w[i]) * mn) * float(np.linalg.norm(V[:, i]))
        res = max(res, float(np.linalg.norm(r)) / max(den, 1e-300))

    # ---- frequencies --------------------------------------------------------
    tol0 = _ZERO_RTOL * lam_scale
    zero = np.abs(w) <= tol0
    neg = w < -tol0
    omega = np.sqrt(np.clip(w, 0.0, None))
    omega[zero] = 0.0
    omega[neg] = np.nan
    freq = omega / (2.0 * math.pi)
    with np.errstate(divide="ignore", invalid="ignore"):
        period = np.where(freq > 0.0, 1.0 / np.where(freq > 0.0, freq, 1.0),
                          np.inf)
    period[neg] = np.nan

    if np.any(zero):
        idx = np.where(zero)[0] + 1
        warnings.append(
            f"{int(zero.sum())} mode(s) have zero frequency (mode "
            f"{', '.join(str(i) for i in idx[:6])}"
            f"{', ...' if len(idx) > 6 else ''}): "
            + ("the pre-stress is at the stability limit or the structure is "
               "a mechanism." if pre else
               "the structure is a mechanism (or can move as a rigid body) "
               "with these supports."))
    if np.any(neg):
        warnings.append(
            f"{int(neg.sum())} mode(s) have omega^2 < 0: the pre-stressed "
            f"structure is unstable (the pre-stress exceeds the buckling "
            f"load); their frequencies are reported as NaN.")
        if not use_dense:
            notes.append("The sparse solver returns the eigenvalues nearest "
                         "zero; more negative ones may exist.")

    cum = ratio.sum(axis=0)
    for a in range(3):
        if free_mass[a] > 0.0 and cum[a] < 0.90 and want < nd:
            notes.append(f"The {want} mode(s) capture {100 * cum[a]:.1f} % "
                         f"of the {_AXES[a]} mass; ask for more modes if a "
                         f"90 % modal mass is needed.")

    full = np.zeros((ndof, want))
    full[fd] = V
    info = (f"{method}; {want} of {nd} modes; {mass_type} mass"
            + ("; with geometric (pre-stress) stiffness" if pre else ""))
    return ModalResult(
        omega2=w, omega=omega, frequencies=freq, periods=period, modes=full,
        participation=Lp, effective_mass=eff, effective_mass_ratio=ratio,
        dof_mass=dof_mass, total_mass=total_mass, free_mass=free_mass,
        free_dofs=fd, mass_type=mass_type, n_free=nf, n_dynamic=nd,
        method=method, info=info, warnings=warnings, notes=notes,
        orthogonality_error=orth, stiffness_error=stiff, residual=res)


def _sparse_modes(Kf, Mf, dyn, ml, kk, lam_scale):
    """Shift-invert Lanczos + massless recovery + Rayleigh-Ritz."""
    Kc = Kf.tocsc()
    Mc = Mf.tocsc()
    how = "sigma = 0"
    try:
        w, V = spla.eigsh(Kc, k=kk, M=Mc, sigma=0.0, which="LM", tol=0.0)
        if not (np.all(np.isfinite(w)) and np.all(np.isfinite(V))):
            raise RuntimeError("non-finite eigenpairs")
    except Exception:
        shift = 1e-6 * lam_scale
        w, V = spla.eigsh(Kc, k=kk, M=Mc, sigma=-shift, which="LM", tol=0.0)
        how = f"sigma = -{shift:.3g} (K is singular)"
    if len(ml):
        K00 = Kc[ml][:, ml].tocsc()
        K0m = Kc[ml][:, dyn]
        rhs = K0m @ V[dyn]
        V[ml] = -np.asarray(spla.spsolve(K00, rhs)).reshape(len(ml), -1)
    Kr = V.T @ (Kc @ V)
    Mr = V.T @ (Mc @ V)
    Kr = 0.5 * (Kr + Kr.T)
    Mr = 0.5 * (Mr + Mr.T)
    w, Y = sla.eigh(Kr, Mr)
    V = V @ Y
    method = (f"sparse shift-invert Lanczos (ARPACK, {how}) + Rayleigh-Ritz, "
              f"{len(dyn)} dynamic DOF"
              + (f", {len(ml)} massless DOF recovered" if len(ml) else ""))
    return V, w, method


# =========================================================== model level ====
def _find_combination(model, name):
    combos = list(model.combinations)
    for c in combos:
        if c.name == name:
            return c
    key = str(name).strip().lower()
    hits = [c for c in combos if c.name.strip().lower() == key]
    if not hits and key:
        first = key.split()[0]
        hits = [c for c in combos
                if (c.name.strip().split() or [""])[0].lower() == first]
    if len(hits) == 1:
        return hits[0]
    raise ValueError(f"Load combination {name!r} not found (or not unique); "
                     f"available: {', '.join(c.name for c in combos)}.")


def prestress_forces(model, prestress, factor=1.0):
    """Member forces for the geometric stiffness of :func:`modal_analysis`.

    ``prestress`` is a combination name (or :class:`LoadCombination`) whose
    linear static member forces are used (self-straining T included with the
    combination's T factor), or an (m,) array of member forces [N].
    Returns ``(N, description)``; ``(None, "")`` when ``prestress`` is None.
    """
    if prestress is None:
        return None, ""
    f = float(factor)
    if not np.isfinite(f):
        raise ValueError("prestress_factor must be finite.")
    if isinstance(prestress, str) or hasattr(prestress, "factors"):
        from . import loads as _loads
        from .solver import solve_model
        combo = (_find_combination(model, prestress)
                 if isinstance(prestress, str) else prestress)
        F = _loads.combine(_loads.basic_load_vectors(model), combo)
        r = solve_model(model, F, name=combo.name,
                        T_factor=float(getattr(combo, "T", 0.0)))
        return f * r.N, (f"{f:g} x member forces of {combo.name.strip()}"
                         f" ({combo.label()})")
    N = np.asarray(prestress, dtype=float).reshape(-1)
    if N.shape != (model.geom.n_members,):
        raise ValueError(f"prestress needs one force per member "
                         f"({model.geom.n_members}), got {N.size}.")
    return f * N, f"{f:g} x given member forces"


def modal_analysis(model, n_modes=12, include_dead=True, nodal_masses=None,
                   snow_fraction=0.0, mass_type=LUMPED, prestress=None,
                   prestress_factor=1.0, dense_limit=2400,
                   align_degenerate=True):
    """Natural frequencies, mode shapes and modal mass of a :class:`DomeModel`.

    Parameters
    ----------
    n_modes : int or None
        Modes to return (lowest first); None or "all" = every mode.
    include_dead : bool
        Superimposed dead load as mass (cladding pressure and downward "D"
        nodal loads, divided by g).  The strut mass is always included.
    nodal_masses : dict {node: kg}, list of {"node", "mass"} or (n,) array
        Extra point masses.
    snow_fraction : float
        Fraction of the balanced snow load added as mass (default 0).
    mass_type : "lumped" (default) or "consistent"
        Distribution of the strut mass; the non-structural masses are lumped.
    prestress : None, combination name, LoadCombination or (m,) forces [N]
        Adds the geometric stiffness of that state (stress stiffening).
    prestress_factor : float
        Multiplies the pre-stress forces (e.g. a fraction of lambda_cr).
    dense_limit, align_degenerate
        See :func:`solve_modes`.

    Supports and springs come from ``model.support_conditions()``.
    Returns a :class:`ModalResult` (frequencies [Hz], periods [s], omega
    [rad/s], mass-normalised full-DOF modes, participation factors and
    effective modal mass ratios in x, y, z, total participating mass,
    verification quantities and messages).
    """
    mm = lumped_masses(model, include_dead=include_dead,
                       nodal_masses=nodal_masses,
                       snow_fraction=snow_fraction)
    fixed, springs = model.support_conditions()
    N_pre, desc = prestress_forces(model, prestress, prestress_factor)
    g = model.geom
    res = solve_modes(g.nodes, g.members, model.EA(), mm.non_structural,
                      fixed, line_mass=mm.line_mass,
                      springs=springs or None, n_modes=n_modes,
                      mass_type=mass_type, N_prestress=N_pre,
                      dense_limit=dense_limit,
                      align_degenerate=align_degenerate)
    res.mass_parts = mm.totals()
    res.notes = list(mm.notes) + res.notes
    res.prestress = desc
    return res


def frequency_load_curve(model, combination, load_factors, n_modes=3,
                         **kwargs):
    """Lowest frequencies as the load of one combination is scaled up.

    For every factor ``lam`` the geometric stiffness of ``lam`` x the linear
    member forces of ``combination`` is added (equation (5)).  The curve of
    omega_1^2 against ``lam`` reaches zero at the linear buckling factor
    lambda_cr (vibration criterion of stability); it is useful to show how
    much stiffness the dome loses under its design load.

    ``kwargs`` go to :func:`modal_analysis` (masses, mass_type, ...).
    Returns a dict with ``load_factors`` (q,), ``omega2`` (q, n_modes),
    ``frequencies`` (q, n_modes) [Hz, NaN where omega^2 < 0] and
    ``combination`` (its name).
    """
    if isinstance(combination, str):
        combination = _find_combination(model, combination)
    N_ref, _ = prestress_forces(model, combination, 1.0)
    lf = np.asarray(load_factors, dtype=float).reshape(-1)
    W = np.full((len(lf), int(n_modes)), np.nan)
    Fq = np.full_like(W, np.nan)
    for i, lam in enumerate(lf):
        r = modal_analysis(model, n_modes=n_modes, prestress=N_ref,
                           prestress_factor=float(lam), **kwargs)
        k = min(r.n_modes, W.shape[1])
        W[i, :k] = r.omega2[:k]
        Fq[i, :k] = r.frequencies[:k]
    return {"load_factors": lf, "omega2": W, "frequencies": Fq,
            "combination": combination.name}
