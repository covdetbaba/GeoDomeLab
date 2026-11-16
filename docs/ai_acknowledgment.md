# Acknowledgment of AI Assistance

**Project:** GeoDomeLab 2.0, analysis, stability and fabrication of single-layer geodesic domes
**Course:** CE 7405 Statically Indeterminate Structures, Louisiana State University, Fall 2026, final project
**Author:** Jovdat Babazade

## Tool used

Anthropic Claude, through Claude Code, was used as a coding and writing
assistant during the whole project, as the course AI policy allows.

## What the assistant was used for

* **Code.** The assistant wrote most of the Python source from the author's
  specifications. Each specification fixed the scope, the method and its
  equations, the module layout, the files a task could change, and the tests
  and tolerances the code had to pass. This covers the geometry generator, the
  stiffness and geometric stiffness assembly, the load generators, the
  AISC 360-16 member checks, the sizing loop, the nonlinear, modal and
  fabrication modules, the OpenSees bridge, the Qt interface and the exports.
* **Tests and validation scripts.** It wrote the unit tests, the
  hand-calculation scripts that print every step, and the scripts that compare
  the program with the hand values and with OpenSees.
* **Review and debugging.** Separate assistant sessions reviewed the code
  for mathematical, interface and validation errors. Every finding was
  reproduced before it was fixed, and each fix has a regression test.
* **Documentation.** It drafted the README, the user guide, the theory and
  validation notes and this page.

## What was not taken on trust

No result is accepted because the program produced it. Each feature is
checked against a reference that shares no code with the program
(details in [VALIDATION.md](VALIDATION.md)):

| Check | Independent reference |
|---|---|
| Hand calculations, Models A and B | Stiffness method by hand, in scripts that do not import the program |
| Hand calculations, cases F1 to F10 | Force method, closed forms and AISC 360-16 equations by hand (temperature, settlement, springs, wind, natural frequencies, snap-through, fabrication angles, member strength) |
| Built-in verification suite | Closed forms, published Class I chord factors, Euler's formula, method-of-sections free bodies, membrane theory |
| Nonlinear path | Closed-form von Mises truss; published star dome results (Greco and Venturini 2006; Iu and Bradford 2012) |
| Independent program | OpenSees: statics, modes, buckling and the corotational path, with a 1 ppm control run that the comparison detects |

Some checks exist because AI-generated work contained mistakes. Examples:

1. The slender round HSS strength first followed the Q-factor method of
   AISC 360-10. It now follows the effective-area method of AISC 360-16,
   Section E7.
2. An early "mesh convergence" check compared domes of different frequency.
   Each frequency is a different structure, so the check was removed and
   replaced by method-of-sections and membrane-theory checks.
3. The default load combinations first missed ASCE 7-16 combination
   1.2D + 1.6S + 0.5W and any live load combination. A review found it and
   both were added.
4. Linear buckling was first reported only for the combination that governs
   the member checks. It is now computed for every strength combination and
   the lowest load factor is reported.

## The author's part

The author defined the problem and the scope, made the engineering decisions
(the pin-jointed truss idealisation, the methods, the case-study loads, the
validation structures and their tolerances, OpenSees as the independent
program), reviewed the equations against the course text and the cited
references, and accepted or rejected every change.

## Reproducibility

```bash
python run_all_checks.py          # every test, verification and validation step
python run_geodome.py --verify    # the built-in verification suite
python run_geodome.py --demo      # the 6V case study
```
