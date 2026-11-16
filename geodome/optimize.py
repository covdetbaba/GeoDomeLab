"""
Fully-stressed-design auto-sizing of the strut families.

A geodesic dome is fabricated from a handful of repeated strut lengths, so the
design variables here are the *families* (A, B, C, ...) rather than the
individual members -- typically six to ten variables instead of several
hundred.  That matches how the structure is actually built and keeps the search
space small enough for a simple, transparent, fully-stressed iteration:

    repeat
        analyse the current design
        for every family
            take the envelope tension and compression demands
            pick the lightest catalogue section whose LRFD capacity covers them
            at the target utilisation, while respecting KL/r <= 200
        if nothing changed and the deflection limit is satisfied: stop
        otherwise, if the deflection limit is violated, step every family up

The result is a local optimum, not a global one, but it is the classical
fully-stressed design and it converges in a handful of iterations because the
force distribution in a triangulated dome is almost independent of the member
sizing.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from . import buckling as bk
from . import sections as sx
from .postprocess import run_analysis


@dataclass
class SizingStep:
    iteration: int
    sections: list
    mass: float
    max_util: float
    deflection_ratio: float
    note: str = ""


@dataclass
class SizingResult:
    converged: bool
    iterations: int
    history: list = field(default_factory=list)
    start_mass: float = 0.0
    final_mass: float = 0.0
    messages: list = field(default_factory=list)

    @property
    def saving(self):
        if self.start_mass <= 0:
            return 0.0
        return 1.0 - self.final_mass / self.start_mass

    def table(self):
        return [{
            "iteration": s.iteration,
            "mass_kg": round(s.mass, 1),
            "max_utilisation": round(s.max_util, 3),
            "deflection": f"span/{s.deflection_ratio:.0f}",
            "sections": ", ".join(s.sections),
            "note": s.note,
        } for s in self.history]


def _select_section(model, L_design, N_t, N_c, target, min_index=0,
                    slenderness_limit=200.0):
    """Lightest catalogue section covering the demands of one family."""
    mat = model.material
    K = model.K_effective
    need_t = abs(N_t) / max(target, 1e-6)
    need_c = abs(N_c) / max(target, 1e-6)

    for i in range(min_index, len(sx.CATALOGUE)):
        s = sx.CATALOGUE[i]
        if K * L_design / s.r > slenderness_limit:
            continue
        if bk.tension_capacity(s, mat, model.phi_t) < need_t:
            continue
        cap_c, *_ = bk.compression_capacity(s, mat, L_design, K, model.phi_c)
        if cap_c < need_c:
            continue
        return s.name
    return sx.SECTION_NAMES[-1]


def auto_size(model, target_util=0.95, max_iter=12, enforce_deflection=True,
              slenderness_limit=200.0, progress=None):
    """Run the fully-stressed design loop, mutating ``model.group_sections``.

    Returns a :class:`SizingResult`.  ``model`` is left holding the final
    (accepted) section assignment.
    """
    def report(f, msg):
        if progress is not None:
            progress(float(f), str(msg))

    geom = model.geom
    n_groups = geom.n_strut_types()
    start_mass = model.steel_mass()

    history = []
    converged = False
    messages = []
    bump = 0                       # global catalogue offset for deflection

    for it in range(1, max_iter + 1):
        report(0.05 + 0.9 * (it - 1) / max_iter, f"Sizing iteration {it}...")
        res = run_analysis(model, do_buckling=False)

        new_sections = list(model.group_sections)
        for g in range(n_groups):
            sel = geom.strut_group == g
            L_design = float(np.max(geom.strut_lengths[sel]))
            N_t = float(np.max(res.N_max[sel]))
            N_c = float(np.min(res.N_min[sel]))
            new_sections[g] = _select_section(
                model, L_design, N_t, N_c, target_util,
                min_index=bump, slenderness_limit=slenderness_limit)

        note = ""
        changed = new_sections != list(model.group_sections)
        defl_ok = res.deflection_ok or not enforce_deflection

        history.append(SizingStep(
            iteration=it, sections=list(model.group_sections),
            mass=model.steel_mass(), max_util=float(res.governing_util),
            deflection_ratio=float(res.deflection_ratio),
            note=("deflection governs" if not res.deflection_ok else "")))

        if not changed and defl_ok:
            converged = True
            break

        if not changed and not defl_ok:
            # Sections are stress-adequate but the dome is too flexible: raise
            # the floor of the catalogue search and try again.
            bump += 1
            if bump >= len(sx.CATALOGUE):
                messages.append("Deflection limit cannot be met with the "
                                "available catalogue; relax the limit or "
                                "increase the dome frequency.")
                break
            note = "stepping every family up to satisfy the deflection limit"
            new_sections = [sx.next_larger(s) for s in model.group_sections]

        model.group_sections = new_sections
        if note:
            history[-1].note = note

    report(0.98, "Final check...")
    final = run_analysis(model, do_buckling=False)
    history.append(SizingStep(
        iteration=len(history) + 1, sections=list(model.group_sections),
        mass=model.steel_mass(), max_util=float(final.governing_util),
        deflection_ratio=float(final.deflection_ratio),
        note="final" + ("" if final.deflection_ok else " (deflection exceeded)")))

    if not converged:
        messages.append(f"Stopped after {max_iter} iterations without a fully "
                        f"stable section set; the last design is reported.")
    if converged and final.governing_util < 0.85 * target_util:
        # Explain the gap rather than leaving it looking like a failed search.
        gm = int(np.argmax(final.util_env))
        grp = int(geom.strut_group[gm])
        cur = model.group_sections[grp]
        i = sx.section_index(cur)
        reason = ("it is the smallest section in the catalogue"
                  if i == 0 else
                  f"the next size down ({sx.SECTION_NAMES[i - 1]}) fails the "
                  f"compression or slenderness check")
        messages.append(
            f"Final utilisation {final.governing_util:.2f} sits below the "
            f"{target_util:.2f} target because the catalogue is discrete: for "
            f"the governing family ({chr(ord('A') + grp)}, {cur}) {reason}. "
            f"Slender struts in a single-layer dome are usually governed by "
            f"buckling, not by stress.")
    if final.governing_util > 1.0:
        messages.append(f"Final design still shows utilisation "
                        f"{final.governing_util:.2f}; the heaviest catalogue "
                        f"section is insufficient.")

    report(1.0, "Auto-sizing complete.")
    return SizingResult(converged=converged, iterations=len(history),
                        history=history, start_mass=start_mass,
                        final_mass=model.steel_mass(), messages=messages)
