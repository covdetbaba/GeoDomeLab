# Acknowledgment of AI Assistance

**Project:** GeoDomeLab — Interactive Analysis and Design of Geodesic Domes
**Course:** CE 7405 Statically Indeterminate Structures, final project
**Student:** Jovdat Babazade

---

## Tool used

Anthropic **Claude** (via Claude Code), used as a coding assistant throughout
development, as permitted by the course AI policy.

---

## What the assistant was used for

**Code authoring.** The assistant wrote the bulk of the Python source across all
sixteen modules, working from my specification of the scope, the methods to use,
the module structure and the feature list. This includes the geometry generator,
the stiffness and geometric-stiffness assembly, the load generators, the AISC
check routines, the optimisation loop, the matplotlib renderers, the Qt
interface and the export layer.

**Debugging.** The assistant diagnosed and fixed the support-placement defect
described in §6.3 of the report, after a parameter study produced displacements
three orders of magnitude too large for every truncation except the hemisphere.

**Documentation.** Docstrings, this document, the user guide, and drafts of the
proposal, report and presentation outline.

**Numerical experiments.** Running the parameter studies and benchmarks whose
results appear in §5.5 and in the timing table.

---

## What was *not* taken on trust

Every quantitative claim in the report is backed by a reference the program did
not generate. Specifically:

| Claim | Independent reference |
|---|---|
| Geometry generator is correct | Published Class I chord-factor tables (2V, 4V); *n* = 10ν² + 2, *m* = 30ν², *f* = 20ν²; Euler's *n* − *m* + *f* = 2 |
| Solver is correct | Closed-form solutions for a single bar, a symmetric tripod and a determinate two-bar node, derived by hand |
| Member forces are statically admissible | Method-of-sections free bodies at three cut heights — a check that never touches the stiffness matrix |
| Force distribution is physically right | Membrane shell theory, *n*φ = *pR*/2, agreeing to 2.6 % |
| Model has no mechanisms | Direct eigenvalue computation of **K**ff, not merely "the solver returned an answer" |
| Buckling formulation is correct | Linearity of the eigenproblem: doubling the reference force must halve every load factor |

Two instances are worth stating plainly, because they are the reason the
verification suite is structured the way it is:

1. **A plausible-looking check was wrong and was removed.** An early "mesh
   convergence" check compared crown deflection across frequencies. That is not
   a convergence test — each frequency is a *different structure*, not a refined
   mesh of the same one — and it was replaced with the method-of-sections and
   membrane-theory checks, which are genuine.

2. **A textbook value did not match, and the discrepancy was investigated rather
   than tuned away.** The hoop-force reversal came out at 35° against the
   smooth-shell value of 45°. Rather than loosening the tolerance, the cause was
   identified (in a triangulated lattice the diagonals also carry hoop force, so
   the ring struts alone under-represent the resultant) and confirmed by showing
   the value converges on 35° at 6V, 8V and 10V. The suite now asserts the
   *reversal*, which is real, and reports the location with its explanation.

The engineering judgments in the report — the choice of idealisation, the
limitations in §6.2, the interpretation of the governing-combination switch in
§5.3, and the recommendation to apply an imperfection knock-down factor to λcr —
are mine.

---

## Hand-calculation validation (added 26 October 2026)

At the instructor's suggestion, the program was checked against hand
calculations of two small models (`validation/`). The step-by-step hand
calculation, its arithmetic, the comparison script and the validation report
were prepared with the AI assistant. To keep the check independent, the hand
calculation script does not import any GeoDomeLab code: it rebuilds the
geometry, element matrices, loads, AISC capacity and buckling factor from first
principles, and its totals are closed by separate equilibrium checks at every
node and for every free body. The program was run through its normal input
path with no code changes.

## Proposal, references and presentation (added 26 October 2026)

The written proposal (on the course template), the presentation and its
speech were drafted with the AI assistant. The reference list was compiled
with its help, and each entry was checked against a library catalogue,
publisher page or citing record before it was used; nothing is cited that
was not found.

## GUI bug fix (added 27 October 2026)

While using the program I found that *Generate dome* did nothing when I went
from a dome with few strut families to one with more (1V to 6V). The
assistant traced it to the order of two lines in `on_generate` (the old
section rows were read back before being rebuilt for the new geometry) and
swapped them. The analysis code was not touched, and the tests, the
verification suite and the hand-calc comparison were re-run with the same
results.

## Reproducibility

Nothing in the report has to be taken on faith:

```bash
python run_geodome.py --verify     # the 75 verification checks
python tests/test_geodome.py       # the 30 unit tests
python run_geodome.py --demo       # the complete case study of §5
```

All 75 checks and all 30 tests pass, and `--demo` regenerates every number in
§5 and §5.5 of the report.
