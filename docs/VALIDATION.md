# How GeoDomeLab is verified and validated

GeoDomeLab is checked in five ways:

1. **Built-in checks.** A 75-check verification suite, plus equilibrium and
   energy checks that every solve reports about itself. Part of the suite
   compares with closed forms and published values; the rest are identities
   (equilibrium, energy, symmetry), which test consistency, not accuracy.
2. **Unit and GUI tests.** 321 tests in 17 files.
3. **Hand calculations.** Two small domes solved by hand (80 comparisons), and
   ten more cases F1 to F10 (422 comparisons) from a hand-calculation file
   that does not import the program.
4. **Published benchmarks.** The 24-member star dome, the von Mises two-bar
   truss and the published geodesic chord factors.
5. **OpenSees.** The same models solved again by an independent finite
   element program (1069 comparisons).

Every number on this page comes from the output files in
`validation/output/` and from a run of `python run_all_checks.py` on
16 November 2026 (Windows, Python 3.13, PyQt6 6.9, OpenSees 3.8.0 through a
separate interpreter). Two more runs on the same day, made before the
review fixes added 31 tests and cases F9 and F10:

- Windows, Python 3.10 with PyQt5 5.15: all 286 tests of that version pass,
  and the verification and hand-calculation steps pass.
- Linux, Python 3.12 with PyQt6 6.12, OpenSees in-process: every
  verification and validation step passes, and 283 of the 286 tests pass.
  Two GUI tests that time the mode animation (`test_animation_frame_time`,
  `test_animation_runs_and_stops_cleanly`) failed on that slower machine. A
  third timed test passed when it ran alone.

## Run everything

```bash
python run_all_checks.py
```

It runs every step below in its own process, with Qt offscreen, and prints
one table. It exits with code 0 only when every step passes. Steps that need
OpenSees are skipped, not failed, when OpenSees is missing. The OpenSees step
runs the quick set and writes its report to a temporary folder, so the full
report in `validation/output` stays as it is. GitHub Actions runs the same
command on every push (`.github/workflows/ci.yml`).

Options: `--jobs N` runs the test files N at a time. `--only TEXT` runs only
the matching steps. `--skip-opensees` leaves OpenSees out. `--opensees-full`
runs the full OpenSees matrix and rewrites its report. `--verbose` prints the
output of every step.

Result of the run on 16 November 2026:

| Step | Command | Result |
|---|---|---|
| Unit and GUI tests | `python tests/test_*.py` (each file) | 321 of 321 pass, 1 skipped |
| Built-in verification | `python run_geodome.py --verify` | 75 of 75 pass |
| Models A and B | `python validation/run_validation.py` | 80 of 80 pass |
| Hand file F1 to F10, internal checks | `python validation/hand_calc_ext.py` | 120 of 120 pass, 290 reference values |
| Program vs F1 to F10 | `python validation/run_validation_ext.py` | 422 of 422 pass, all 290 reference values used |
| OpenSees, quick set | `python validation/opensees_crosscheck.py --quick` | 486 of 486 pass |
| OpenSees, full set | `python validation/opensees_crosscheck.py` | 1069 of 1069 pass (the report in `validation/output`, made the same day) |

The skipped test compares the in-process and subprocess ways of running
OpenSees, and it needs both on the same computer.

## 1. Built-in checks

### 1.1 Verification suite (75 checks)

Run it with `python run_geodome.py --verify`, the Verify button (F11), or
Validation > Run verification suite. Source: `geodome/verification.py`.

| Group | Checks | What is compared | Worst difference |
|---|---|---|---|
| Closed-form truss solutions | 8 | single bar, tripod and two-bar node: displacements, bar forces, reactions | 2.4e-16 (relative) |
| Geodesic geometry | 33 | 2V and 4V chord factors vs the published tables; vertex and strut counts of 1V to 6V spheres and Euler's V − E + F = 2; Maxwell's count m − (3n − 6) = 0 for closed spheres; no zero-energy mode with supports; surface area and crown height of 4V and 8V hemispheres | chord factors 1.6e-6 (tolerance 2e-5); faceted surface area 1.9 % below the sphere (tolerance 3 %) |
| Section and code checks | 6 | CHS area, inertia and radius of gyration; AISC E3 elastic and inelastic branches and the Euler load, with E, Fy, D and t typed in (not taken from the program's material or section) | 0 |
| Global behaviour | 28 | equilibrium residual, strain-energy identity and load/reaction balance for every combination; five-fold symmetry; free bodies cut at three heights (method of sections); membrane theory n_φ = pR/2 at three heights; hoop force sign along the meridian; stiffening with frequency; buckling eigenvalue scaling | equilibrium residual 6.4e-15, energy identity 1.1e-15, symmetry 6.9e-15, load/reaction balance 3.5e-10 (tolerance 1e-6); membrane theory 2.6 % (tolerance 5 %, because the faceted truss is not a smooth shell) |

### 1.2 Checks inside every solve

Every static solve computes three checks of its own answer:

- the equilibrium residual ‖K u − F‖ / ‖F‖ at the free DOFs,
- the balance of applied loads and reactions, Σ F + Σ R = 0,
- the strain-energy identity ½ F^T u = Σ N² L / (2EA), extended for initial
  strains, settlements and springs (see THEORY.md).

The Summary tab shows the worst equilibrium residual and energy mismatch over
all combinations. The HTML report has a Solution verification section. Every
nonlinear step stores its residual. Every modal result stores
‖Φ^T M Φ − I‖, ‖Φ^T K Φ − Ω²‖ and the eigen residual.

## 2. Unit and GUI tests

Each file runs alone (`python tests/test_x.py`, no pytest needed) and also
under pytest. The GUI tests open the real window offscreen and click through
it. They also catch any exception raised inside a Qt slot.

| File | Tests | What it covers |
|---|---|---|
| `test_geodome.py` | 30 | geometry, sections, loads, solver, AISC checks, buckling |
| `test_foundation.py` | 44 | temperature, settlement, springs, custom supports, custom geometry, JSON v2, mechanism detection, default combinations, unused-load and aluminium messages, wind direction, unbalanced snow, snow on deep domes, lowest buckling factor over the strength combinations |
| `test_nonlinear.py` | 24 | tangent vs finite differences, von Mises closed form, star dome, 6V limit point, imperfections (including a crown snap close to zero load, a buckling-mode sign that is the same on every platform, and the local rise left at the hub), the double bifurcation of a three-bar tent |
| `test_modal.py` | 25 | closed-form frequencies, M-orthonormality, effective mass sum rule |
| `test_fabrication.py` | 20 | chord factors 2V to 8V, valences, angles, dihedrals, panel areas, cut lists |
| `test_study.py` | 20 | parameter studies, CSV round trip, auto-sizing with the B4.2 wall and deflection floors |
| `test_benchmarks.py` | 21 | every example model and its reference answers, the pass rule of run_validation.py |
| `test_opensees.py` | 20 + 1 skipped | OpenSees bridge, strain measure, thermal sign, statics, modes, paths |
| `test_validation_ext.py` | 15 | the F1 to F10 comparison |
| `test_cad.py` | 11 | DXF export read back: entity counts and coordinates |
| `test_report_extras.py` | 14 | HTML report sections for nonlinear, modal, fabrication and study results; no NaN text for a custom truss |
| `test_gui_generate.py` | 2 | Generate dome from 1V to 10V, up and down |
| `test_gui_inputs.py` | 26 | examples, custom models, supports, loads, CSV import, save and open, mechanism and file-check messages, busy file actions, small-screen size |
| `test_gui_analysis.py` | 15 | the integrated window: menus and tabs, examples, 6V analysis, nonlinear, modal, fabrication, study, validation without OpenSees, the script export menu entry, results cleared after input changes, HTML and DXF export, save and open |
| `test_gui_nonlinear.py` | 7 | nonlinear view, path scrubbing, imperfection study, failed imperfection signs |
| `test_gui_modal_fab.py` | 16 | Modes view and animation, Fabrication tab, DXF menu |
| `test_gui_study_validation.py` | 11 | Study view, Validation menu with and without OpenSees |
| **Total** | **321 + 1 skipped** | |

## 3. Hand calculations

### 3.1 Models A and B (80 comparisons)

`validation/hand_calc.py` solves two small domes by hand. It does not use the
program. `validation/run_validation.py` builds the same models with the
program's own generator and compares.

- **Model A**, 1V pentagonal cap of a 5 m sphere: 6 nodes, 10 struts, 3 free
  DOF. Solved by the direct stiffness method with every matrix written out.
  Case A1 is a crown load (+20 kN in x, −60 kN in z). Case A2 is 1.2D + 1.6S.
- **Model B**, 1V dome on its lower ring: 11 nodes, 25 struts, 18 free DOF.
  Solved by virtual work with five-fold symmetry (3 unknowns). Case B1 is a
  50 kN crown load. Case B2 is 1.2D + 1.6S. Model B reaches below the
  sphere's equator: 5 of its 10 band panels face down, so they carry no snow
  (snow acts on the 5 top panels and the 5 upward band panels).

Compared: geometry, nodal loads, stiffness matrix terms, displacements,
member forces, reactions, AISC capacities and utilisation, and λ_cr. A row passes when |program − hand| / |hand| ≤ 1e-12
(1e-10 %), or |program| ≤ 1e-9 in the row unit where the hand answer is 0.

Result: 80 of 80 pass. The largest relative difference is 1.2e-13 % (67
rows). Where the hand answer is zero (13 rows) the largest program value is
1.0e-11 kN/m. Table: `validation/output/comparison.md` and `.csv`.

Example values (they also appear in the Examples menu): Model A, A1: crown
ux = +0.2724 mm, uz = −1.0697 mm, spoke force −32.23 kN. Model A, A2:
utilisation 0.2726, λ_cr = 5506.1. Model B, B1: crown uz = −1.6451 mm.

### 3.2 Cases F1 to F10 (422 comparisons)

`validation/hand_calc_ext.py` writes step-by-step hand calculations to
`validation/output/hand_calc_ext.md`. It gives 290 reference values and checks
itself 120 times (120 of 120 within tolerance). It never imports the program,
and the comparison script checks that the reference values are computed
before the program is loaded.

| Case | Structure | What is checked |
|---|---|---|
| F1 | Heated fixed-fixed bars | N = −EAαΔT, reactions, a stepped bar, the ΔT that reaches the Euler load |
| F2 | Model A under a temperature change | free expansion of the crown, stress-free spokes, reactions; one heated spoke |
| F3 | Model A with a support movement | a 3 x 3 stiffness solve for a settling support and a spreading support |
| F4 | Planar three-bar truss | the force (flexibility) method for a load, a heated bar and a settlement, independent of the stiffness method |
| F5 | Bar and elastic spring | springs in series and in parallel, a heated bar on a spring |
| F6 | von Mises two-bar truss | the exact path P(w), both limit points, the snap-through target, linearised λ_cr and the sway root |
| F7 | Model A natural frequencies | crown frequencies and effective masses, with and without the dead load as mass |
| F8 | 2V geodesic geometry | chord factors, axial, bevel, dihedral and face angles, areas and counts; chord factors also vs the printed values of the Domerama chord-factor chart |
| F9 | Wind on Model B | qz = 0.613 Kz Kzt Kd V², Cp of every node along the meridian, the nodal forces −qz G Cp A n and their resultant, with the default wind data. This checks the implementation of the wind field; the 7-point Cp table itself is an illustrative simplification, not ASCE 7-16 Fig. 27.3-2 |
| F10 | AISC 360-16 strength of round HSS | the inelastic E3 branch (CHS 88.9 x 4.0, L = 2 m: φPn = 239.578 kN), the B4.2 design wall 0.93 t through a model (223.996 kN, tension 309.095 kN), the E7 effective area of a slender wall (D 323.9 x 4.0, L = 3 m: Ae/Ag = 0.93871, φPn = 1112.935 kN) and the 0.45 E/Fy scope limit (φPn = 0) |

`validation/run_validation_ext.py` runs each case through the normal program
API and writes `validation/output/comparison_ext.md` and `.csv`.

Comparison rule: d = |program − hand| / |hand| for a non-zero reference.
Where the exact answer is zero, d = |program| / s, with s the largest hand
value of the same unit in the same case. Tolerances:

| Kind | Tolerance |
|---|---|
| linear statics, buckling and geometry | 1e-9 |
| modal | 1e-9 |
| nonlinear path | 1e-6 |
| counts | exact |
| published 5-decimal chord factors | 5e-6 (absolute) |
| program identities | 1e-10 |

| Case | Rows | Pass | Worst d (linear) | Worst d (modal) | Worst d (nonlinear) | Worst identity |
|---|---|---|---|---|---|---|
| F1 | 21 | 21 | 1.6e-15 | - | - | 7.4e-17 |
| F2 | 52 | 52 | 6.6e-16 | - | - | 2.2e-16 |
| F3 | 66 | 66 | 1.2e-15 | - | - | 1.6e-16 |
| F4 | 78 | 78 | 7.2e-16 | - | - | 3.0e-16 |
| F5 | 21 | 21 | 2.2e-16 | - | - | 1.7e-16 |
| F6 | 32 | 32 | 1.2e-13 | - | 4.9e-10 | 0.0e+00 |
| F7 | 34 | 34 | - | 7.0e-16 | - | 4.4e-16 |
| F8 | 52 | 52 | 3.2e-15 | - | - | - |
| F9 | 52 | 52 | 2.6e-15 | - | - | 2.8e-16 |
| F10 | 14 | 14 | 0 | - | - | - |
| **All** | **422** | **422** | | | | |

The two printed chord factors (0.54653 and 0.61803) differ from the program by
3.1e-6 and 4.0e-6, which is the rounding of the published 5 decimals.

Some F quantities have no single program output, for example a single heated
bar (the input heats every strut). These are computed by a short script step
on the program's own objects. The `source` column of the CSV says which.

## 4. Published benchmarks

### 4.1 The 24-member star dome

The shallow star dome of Hangai and Kawamata is a classic snap-through test.
GeoDomeLab uses the geometry given by Greco & Venturini (2006, Sec. 5.3) and
Iu & Bradford (2012, Sec. 5.5): crown 8.216 cm above the supports, inner ring
of radius 25 cm at 6.216 cm, six pinned supports at radius 50 cm. The outer
nodes use the exact hexagon value 43.30127 cm. The papers print 43.3, which
changes the first limit load by 2.5e-6. Results are given as P/EA, so the
section does not matter. The preset uses a crown load of 1e-4 EA, so the load
factor reads 10⁴ P/EA. The original 1972 curves of Hangai and Kawamata were
not used, because Iu & Bradford report that they disagree with later
solutions.

| Quantity | Published | GeoDomeLab |
|---|---|---|
| Crown displacement at the first limit point | "near 0.75 cm" (Greco & Venturini 2006) | 0.76844 cm |
| First limit load | about 3 kN for EA = 10 796 kN, read from a figure (Greco & Venturini 2006) | P/EA = 3.15655e-4, that is 3.41 kN for EA = 10 796 kN |
| First zero-load state | crown about 19 mm, node 2 about −1.16 mm (Iu & Bradford 2012) | crown 1.8838 cm, node 2 moves up 0.1162 cm |
| Exact mirror states (closed form, N = 0) | crown at 2 x 2.000 cm and dome at 2 x 8.216 cm | P = 0 and N = 0 to 1e-12 |
| Linear buckling | | λ_cr P/EA = 1.52097e-3, so λ_limit/λ_cr = 0.2075 |

The published values are rounded or read from figures. GeoDomeLab matches
them within that precision. The 5- and 6-digit values in the right-hand column
are not published numbers; they are computed in this work (see 4.1.1). The first limit point, the zero-load state and
the load minimum (−2.76000e-4 EA at 3.0278 cm) are also found by an
independent Newton solver written inside `tests/test_benchmarks.py`.
`tests/test_nonlinear.py` checks the limit point with both arc-length and
displacement control, and the second maximum (8.8654e-3 EA near 10.54 cm).
OpenSees follows the same path (case N-star-dome, refined limit λ differs by
8.75e-14).

The linear buckling factor is 4.8 times the true first limit load. This is
why GeoDomeLab has a nonlinear analysis.

#### 4.1.1 Star dome reference values (computed in this work)

The precise star dome values used by the program's checks
(`geodome/benchmarks.py`, `tests/test_benchmarks.py`,
`tests/test_nonlinear.py`, the OpenSees cross-check) were computed in this
work for an engineering-strain corotational truss under crown displacement
control. They are reproduced by the independent Newton solver written inside
`tests/test_benchmarks.py`, which does not use the program's solver, and they
agree with Greco & Venturini (2006, Sec. 5.3 and Fig. 8) and Iu & Bradford
(2012, Sec. 5.5) to the precision those papers print. EA is any value; loads
are P/EA, displacements are in cm (the papers' unit).

| Quantity | Value |
|---|---|
| Initial crown stiffness | 0.090430 EA per m (9.04299e-4 EA per cm) |
| Linearised buckling load | λ_cr P/EA = 1.52097e-3 |
| First limit point | P/EA = 3.15655e-4 at a crown deflection of 0.76844 cm (node 2 −0.0491 cm) |
| First zero-load state | crown 1.8838 cm, node 2 moves up 0.1162 cm |
| Load minimum | P/EA = −2.76000e-4 at 3.0278 cm |
| Exact mirror state (every bar at its original length) | crown 4.0000 cm, node 2 at 0 |
| Second load maximum | P/EA = 8.8654e-3 at about 10.54 cm |

A grid sample of the flat first maximum gives 0.7685 cm; refining it by
golden-section search on the path gives 0.76844 cm, the value used.

### 4.2 The von Mises two-bar truss

Two bars of half span b = 2.5 m and rise h = 0.125 m meet at a loaded apex
(CHS 88.9 x 4.0, P_ref = 10 kN). With engineering strain the exact path is
known in closed form (the problem studied by Pecknold, Ghaboussi & Healey
1985):

- P(w) = 2EA y (1/l − 1/L), with y = h − w, l = √(b² + y²), L = √(b² + h²);
- limit point at l* = (b²L)^(1/3);
- linearised critical load λ_cr = 2EA h³ / (L b² P_ref).

GeoDomeLab traces the whole path with arc-length: through the limit point,
along the unstable branch and onto the inverted stable branch. Every traced
point is compared with P(w) (case F6, 32 rows, worst 4.9e-10 on the path).
The limit load is P* = 10.24 kN at w* = 52.86 mm, and λ_cr = 5.33, so the
linear buckling factor overestimates the snap-through load 5.2 times.
OpenSees `corotTruss` gives the same path (vs the closed form 2.6e-13) and the
same limit load (1.4e-15).

### 4.3 Geodesic chord factors

Chord factor = strut length / sphere radius (Kenner 1976 describes the Class I
breakdown).
The printed reference values come from the Domerama chord-factor chart
(<https://domerama.com/calculators/chord-charts>, "method 1" icosahedral,
read 16 November 2026), which gives five decimals.

- 2V: the closed forms 2 sin(18°) = 0.618034 and
  2 sin(arccos(1/√5)/4) = 0.546533 agree with the program to 1e-13. The chart
  values 0.54653 and 0.61803 agree to the printed 5 decimals.
- 2V and 4V: the verification suite checks 8 chord factors against the chart
  values (2V 0.54653 and 0.61803; 4V 0.25318, 0.29453, 0.29524, 0.29859,
  0.31287 and 0.32492), all within the rounding of the fifth decimal.
- 3V to 8V: `tests/test_fabrication.py` checks every chord factor and strut
  count against the Domerama chord-factor chart. 2V to 6V agree within the
  printed decimals. The 7V and 8V chart rounds or merges some lengths, so
  there each length is matched to the nearest printed value within 1e-5.
- Full spheres have 10ν² + 2 nodes, 30ν² struts and 20ν² faces, 12 hubs of
  5 struts and the rest of 6.

Kenner's book itself was not read for this project, and the Domerama page
does not name its source. The printed values were therefore checked against
the closed forms and the Domerama chart only. (2V B = 2 sin 18° = 1/φ is the edge of the icosidodecahedron
formed by the icosahedron's edge midpoints, not the icosahedron edge, whose
chord factor is 1.05146, the 1V value.)

## 5. OpenSees cross-check (1069 comparisons)

`validation/opensees_crosscheck.py` solves each structure in GeoDomeLab and
in OpenSees 3.8.0 (McKenna, Scott & Fenves 2010, through OpenSeesPy, Zhu,
McKenna & Scott 2018) and compares the results.

**Shared inputs.** The two programs share only the input data written by
`opensees_bridge.to_opensees_spec`: node coordinates, connectivity, E and A,
supports, the factored nodal load vector of the combination, the thermal
strain T·α·ΔT, the settlements T·δ and, for modes, the masses.

**Built by OpenSees itself.** `Truss` and `corotTruss` elements with an
`Elastic` material, `InitStrainMaterial(−αΔT)` for temperature, `sp`
constraints for settlement, `zeroLength` springs to fixed ground nodes,
element mass `-rho` (lumped, or consistent with `-cMass 1`) plus nodal
`mass`, `eigen`, and Newton iteration with `DisplacementControl` for the
nonlinear path.

**Two probes first.** OpenSees' strain measure and thermal sign were found
numerically, not assumed. A rotated and stretched bar shows that `corotTruss`
uses engineering strain on the original length, like GeoDomeLab (difference
1.5e-15). Green-Lagrange and logarithmic strain would differ by 6.6e-2 and
6.9e-2. A
heated bar between walls shows that `InitStrainMaterial` needs −αΔT. A third
probe raises E by 1 ppm in OpenSees only and sees a 1.0e-6 difference, so the
agreement below is a real agreement, not a program compared with itself.

| Section | Cases | Comparisons | Worst difference | Target |
|---|---|---|---|---|
| Formulation probes | 3 | 6 | 1.5e-15 | 1e-9 |
| Linear statics: 1V to 8V domes on pinned, roller and elastic supports, every default combination, temperature and settlement, a custom support mix, other dome fractions, Models A and B, the three-bar truss and the 6V case study | 34 | 914 | u 2.19e-13, N 1.72e-13, R 6.14e-14 | 1e-8 |
| Modal: the same domes, all modes of 1V to 3V, lumped and consistent mass, Model A, star dome, three-bar truss | 35 | 97 | f 2.72e-13, 1 − MAC 3.3e-15, repeated pairs 4.0e-15 | 1e-6 |
| Nonlinear paths: von Mises (shallow and deep), star dome, 2V, 3V, 4V and 6V domes, 3V on springs | 8 | 47 | λ at matched states 4.14e-11, displacements 3.86e-11 | 1e-5 |
| Exported OpenSeesPy scripts: static, modal, nonlinear | 3 | 5 | 6.3e-15 (static), 9.9e-13 (nonlinear residual) | 1e-8 |
| **Total** | **83** | **1069** | | |

Static and modal differences are normwise relative (largest error over the
largest value). Repeated-frequency mode pairs are compared as subspaces,
because any rotation of a repeated pair is an equally valid pair. For the
nonlinear paths OpenSees is driven to every GeoDomeLab state by displacement
control. Where the path snaps back, control switches to another DOF.
OpenSees then finds λ and every displacement by itself. OpenSees' own
ArcLength integrator turned back at the limit points of both snap-through
benchmarks, so it is used only as an extra check of equilibrium.

**What OpenSees does not check.** The load generation and the
non-structural masses are shared inputs, so an error in them would not show
here. Checked by hand elsewhere: the dead and snow tributary loads and the
crown load on Models A and B, the wind field on Model B (F9), and E, α and ρ
through the typed-in constants of the hand files. Not checked by a hand
calculation: the unbalanced snow pattern and the default combination factors
(unit tests only), and the Cp table itself (an illustrative simplification).
The AISC member checks (hand-checked in Model A and F10), the linear buckling
eigenvalue (OpenSees has no linear buckling analysis) and the fabrication
geometry are not part of OpenSees.

**Running it.** In-process when `import openseespy.opensees` works in the
current Python. Otherwise a separate Python named in the environment variable
`GEODOME_OPENSEES_PYTHON`, or a conda environment called `opensees`, is used
through a subprocess. Without either, the cross-check is skipped.

```bash
python validation/opensees_crosscheck.py --quick   # 486 comparisons, a few seconds
python validation/opensees_crosscheck.py           # 1069 comparisons, rewrites the report
```

Report: `validation/output/opensees_crosscheck.md` and `.csv`.

## 6. Rerun each check

| Check | Command | Output |
|---|---|---|
| Everything | `python run_all_checks.py` | the table on screen |
| Verification suite | `python run_geodome.py --verify` | on screen; also Validation menu and the HTML report |
| One test file | `python tests/test_geodome.py` (or `pytest`) | on screen |
| Hand calculation, Models A and B | `python validation/hand_calc.py` | the derivation on screen |
| Program vs Models A and B | `python validation/run_validation.py` | `validation/output/comparison.md`, `.csv`, sketches, model files |
| Hand calculation F1 to F10 | `python validation/hand_calc_ext.py` | `validation/output/hand_calc_ext.md` |
| Program vs F1 to F10 | `python validation/run_validation_ext.py` | `validation/output/comparison_ext.md`, `.csv` |
| OpenSees | `python validation/opensees_crosscheck.py [--quick] [--out DIR]` | `validation/output/opensees_crosscheck.md`, `.csv` |
| A validation model in the window | `python validation/open_validation_in_gui.py B1` | the GUI |
| Example self-checks | `python -m geodome.benchmarks` | on screen |

In the window, Validation > Run hand-calculation validation and Run OpenSees
cross-check show the same comparisons in the Validation tab without
rewriting the files.

## 7. Limits of the validation

- No published linear analysis of a full geodesic dome with numbers that
  could be reproduced was found. Linear statics of domes therefore rest on
  the hand calculations (Models A and B, F2, F3, F7), the built-in identities
  and the OpenSees cross-check.
- The wind pressure coefficient along the meridian is a simplified curve that
  follows the trend of the ASCE 7 domed-roof chart. It is not the full code
  procedure, and F9 checks only that the program applies it as documented.
  The vertex-normal lumping of the wind pressure gives 5.1 % less than the
  exact resultant at 2V, 1.4 % at 4V and 0.6 % at 6V.
- The unbalanced snow pattern (1.5 p on x ≥ 0, 0.5 p on x < 0) is a simple
  test pattern, not the ASCE 7-16 curved-roof rule.
- The member checks follow AISC 360 (steel). For aluminium they do not apply;
  the program says so in every run.
- The nonlinear analysis is elastic with pin-ended bars. It reports when a
  strut reaches its Euler load or design strength but does not reduce its
  stiffness. Joint stiffness, strut buckling and yielding along the path are
  not modelled.
- Only the first load maximum is refined. Later peaks are path samples.
- The self-straining load T is not part of the nonlinear path.
- Published benchmark values read from figures or rounded in the text are
  matched only to their printed precision.

## References

- Greco, M., Venturini, W. S. (2006). Stability analysis of three-dimensional
  trusses. *Latin American Journal of Solids and Structures* 3(3), 325-344.
- Iu, C. K., Bradford, M. A. (2012). Higher-order non-linear analysis of steel
  structures. Part I: Elastic second-order formulation. *Advanced Steel
  Construction* 8(2), 168-182. doi:10.18057/IJASC.2012.8.2.5.
- Kenner, H. (1976). *Geodesic Math and How to Use It*. University of
  California Press, Berkeley.
- McKenna, F., Scott, M. H., Fenves, G. L. (2010). Nonlinear finite-element
  analysis software architecture using object composition. *Journal of
  Computing in Civil Engineering* 24(1), 95-107.
  doi:10.1061/(ASCE)CP.1943-5487.0000002.
- Pecknold, D. A., Ghaboussi, J., Healey, T. J. (1985). Snap-through and
  bifurcation in a simple structure. *Journal of Engineering Mechanics*
  111(7), 909-922. doi:10.1061/(ASCE)0733-9399(1985)111:7(909).
- Zhu, M., McKenna, F., Scott, M. H. (2018). OpenSeesPy: Python library for
  the OpenSees finite element framework. *SoftwareX* 7, 6-11.
  doi:10.1016/j.softx.2017.10.009.
- Domerama, geodesic dome chord-factor charts,
  https://domerama.com/calculators/chord-charts (accessed November 2026).

More references for the methods are in [THEORY.md](THEORY.md).
