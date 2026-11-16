# GeoDomeLab — Interactive Analysis and Design of Geodesic Domes by the Direct Stiffness Method

**Matrix Structural Analysis · Final Project Report**
Source code: `GeoDomeLab/` · Launcher: `run_geodome.py`

---

## Abstract

GeoDomeLab is a Python application that generates, analyses, checks and sizes
single-layer geodesic domes. The dome is modelled as a three-dimensional
pin-jointed space truss and solved by the direct stiffness method; members are
checked to AISC 360 LRFD; system stability is assessed by a linear buckling
eigenanalysis using the geometric stiffness matrix; and the strut families can
be auto-sized by fully stressed design. A Qt graphical interface lets the user
define the structure and the loading and view every result interactively, and a
single self-contained HTML report can be exported.

The program is verified by 75 independent checks — closed-form truss solutions,
published geodesic chord-factor tables, Euler's and Maxwell's counting formulas,
method-of-sections free bodies, membrane shell theory, symmetry, and
strain-energy and equilibrium identities. All 75 pass; equilibrium residuals are
of order 10⁻¹⁵ and the computed meridional stress resultant agrees with membrane
shell theory to within 2.6 %.

A 30 m span, 15 m rise 6V dome is carried through as a case study. Auto-sizing
reduces the structural steel from 22 445 kg to 7 504 kg, a 66.6 % saving, and in
doing so changes the governing load combination from snow (1.2D + 1.6S) to wind
uplift (0.9D + 1.0W) — a result that would be easy to miss without enveloping
every combination.

---

## 1. Introduction

### 1.1 The structure

A geodesic dome is built by subdividing the faces of a regular icosahedron and
projecting the new vertices onto a sphere. The resulting triangulated lattice is
a space truss: struts carry axial force only, and the curvature of the surface is
what lets a very light structure span a long distance. The dome analysed in the
case study spans 30 m with 60.3 mm diameter tubes.

Three things make it a natural subject for a matrix-methods project:

1. **It cannot be done by hand.** A 6V hemisphere has 196 nodes, 555 struts and
   588 degrees of freedom. The method is the one from class; only the scale is
   different.
2. **It is three-dimensional.** The extension from the 2D truss element to the
   3D bar is a single extra direction cosine, but it changes the bookkeeping,
   the support conditions and the stability behaviour.
3. **It is checkable.** Published chord-factor tables fix the geometry, membrane
   shell theory predicts the force distribution, and Euler's and Maxwell's
   formulas fix the topology — so the program can be validated against things it
   did not compute.

### 1.2 Scope and limitations

The analysis is **first-order linear elastic** with **pin-jointed** connections.
This is the standard idealisation for a single-layer dome with bolted hub
connectors, and it is conservative for stiffness because real hubs carry some
moment. Second-order (P-Δ) effects are not traced incrementally; instead the
destabilising effect of the axial forces is captured by a *linearised* buckling
eigenanalysis, which gives the critical load factor but not the post-buckling
path. Section 7 states the consequences of this explicitly.

---

## 2. Theory and formulation

### 2.1 Geometry: the Class I breakdown

A regular icosahedron is built in pole orientation: a vertex at each pole and two
staggered rings of five vertices at *z* = ±1/√5. Each of its 20 faces (A, B, C)
is subdivided by barycentric interpolation,

> **P**(*i*, *j*) = [(ν − *i*)**A** + (*i* − *j*)**B** + *j* **C**] / ν,  0 ≤ *j* ≤ *i* ≤ ν

and each new point is projected radially onto the unit sphere, **p** ← **p**/|**p**|.
Points shared between adjacent faces are merged with a KD-tree at a tolerance of
10⁻⁶. The result for a complete sphere is

| quantity | value |
|----------|-------|
| vertices | *n* = 10ν² + 2 |
| struts   | *m* = 30ν² |
| panels   | *f* = 20ν² |

all three of which are asserted in the verification suite together with Euler's
*n* − *m* + *f* = 2.

**Chord factors.** Struts are grouped into families by length; the chord factor
is the family length divided by the sphere radius. These are the numbers a
fabricator works from, and they are published for the standard breakdowns, so
they are the first validation of the generator:

| | computed | published | error |
|---|---|---|---|
| 2V A | 0.5465330 | 0.546533 | 1.1 × 10⁻⁷ |
| 2V B | 0.6180340 | 0.618034 | 1.8 × 10⁻⁸ |
| 4V A | 0.2531850 | 0.253185 | 1.6 × 10⁻⁶ |
| 4V B | 0.2945310 | 0.294531 | 5.6 × 10⁻⁷ |
| 4V C | 0.2952420 | 0.295242 | 6.5 × 10⁻⁷ |
| 4V D | 0.2985880 | 0.298588 | 4.5 × 10⁻⁷ |
| 4V E | 0.3128690 | 0.312869 | 2.2 × 10⁻⁷ |
| 4V F | 0.3249200 | 0.324920 | 9.4 × 10⁻⁷ |

(2V B is 1/φ, the reciprocal golden ratio — the icosahedron edge itself.)

**Truncation.** The dome is cut on a horizontal plane snapped to an existing ring
of struts, because a geodesic mesh can only be cut along its own members. The
supports are then placed on the **topological free edge** of the retained cap —
the nodes lying on an edge that belongs to only one panel — rather than on "all
nodes at the cut elevation". This distinction matters: except at the equator of
an even-frequency dome, the cut elevation contains only a small sub-ring of
nodes, and supporting those alone leaves the real free edge of the mesh
unsupported. Section 6.3 shows what that error looked like before it was found.

### 2.2 The space-truss element

Each strut is a two-node pin-ended bar with three translational degrees of
freedom per node, numbered 3*k*, 3*k*+1, 3*k*+2 for node *k*. With

> **c** = [(*x*₂−*x*₁)/*L*, (*y*₂−*y*₁)/*L*, (*z*₂−*z*₁)/*L*]ᵀ

the axial deformation is *e* = **c**ᵀ(**u**₂ − **u**₁), so the 1 × 6
transformation is **T** = [−**c**ᵀ **c**ᵀ] and

> **k** = **T**ᵀ (*EA*/*L*) **T** = (*EA*/*L*) · ⎡ **cc**ᵀ −**cc**ᵀ ⎤
>                                                ⎣ −**cc**ᵀ  **cc**ᵀ ⎦

This is the same derivation as the plane truss, with a three-component direction
cosine vector instead of two.

### 2.3 Assembly and solution

The global matrix is assembled in sparse coordinate form (36 entries per member)
and converted to CSR. Degrees of freedom are partitioned into free (f) and
supported (s):

> ⎡ **K**ff **K**fs ⎤ ⎧**u**f⎫ = ⎧**F**f⎫
> ⎣ **K**sf **K**ss ⎦ ⎩**u**s⎭   ⎩**F**s⎭

With **u**s = **0**, the free partition is solved by sparse LU,
**K**ff **u**f = **F**f, and the reactions follow by back-substitution,
**R** = **K**sf **u**f − **F**s. Member forces come from the element
compatibility relation,

> *N* = (*EA*/*L*) **c**ᵀ(**u**₂ − **u**₁),  tension positive.

A zero diagonal in **K**ff is trapped before solving and reported as an unstable
structure rather than allowed to produce nonsense.

**Static determinacy.** For a *closed* Class I sphere, *m* = 30ν² and
3*n* − 6 = 3(10ν² + 2) − 6 = 30ν², so **the complete geodesic sphere is exactly
isostatic by Maxwell's count, at every frequency**. Truncating it destroys that
balance — a 3V cap has 120 struts against 3*n* − 6 = 132 — so an open cap
floating in space is a mechanism with 12 internal inextensional modes. The dome
depends on its supports. The verification suite therefore does not merely check
that the solver runs: it computes the eigenvalues of **K**ff directly and
confirms that the smallest is positive, i.e. that no zero-energy mode survives
once the free edge is restrained.

| model | nodes | struts | free DOF | redundancies | λmin/λmax of **K**ff |
|-------|-------|--------|----------|--------------|----------------------|
| 2V pinned | 26 | 65 | 48 | 17 | 2.26 × 10⁻² |
| 4V pinned | 91 | 250 | 213 | 37 | 4.84 × 10⁻³ |
| 6V pinned | 196 | 555 | 498 | 57 | 2.03 × 10⁻³ |
| 4V roller | 91 | 250 | 247 | 3 | 1.57 × 10⁻³ |
| 6V roller | 196 | 555 | 552 | 3 | 5.59 × 10⁻⁴ |

Condition numbers run from 44 to 1 800 — the systems are well conditioned, and
the vertical-roller support option is only three times redundant, which is a
useful thing to be able to show.

### 2.4 Loads

All area loads are lumped to the nodes using tributary areas computed as one
third of each adjacent triangular panel.

| load | applied on | direction |
|------|-----------|-----------|
| self weight | member mass, half to each end | vertical |
| cladding dead | faceted **surface** area | vertical |
| snow, live | horizontally **projected** area | vertical |
| wind | faceted surface area | **normal** to the surface |

Wind uses the ASCE 7 velocity pressure *q*z = 0.613 *K*z *K*zt *K*d *V*²
(SI, *V* in m/s) with a gust factor *G*, and a pressure coefficient that varies
along the meridian from the windward springing (θ = 0°) over the crown (θ = 90°)
to the leeward springing (θ = 180°):

| θ [°] | 0 | 30 | 60 | 90 | 120 | 150 | 180 |
|-------|---|----|----|----|-----|-----|-----|
| *C*p | +0.80 | +0.20 | −0.60 | −1.00 | −0.80 | −0.50 | −0.40 |

This follows the trend of ASCE 7 Figure 27.3-2 for a domed roof and is
deliberately a simplification: it is axisymmetric about the wind axis, whereas
the real field is not. It is adequate to show the structurally important fact —
that wind on a dome is dominantly **uplift** — and that fact is what changes the
governing load combination in Section 5.

Load combinations are user-editable; the defaults are the ASCE 7 strength
combinations 1.4D, 1.2D + 1.6S, 1.2D + 1.0W + 0.5S and 0.9D + 1.0W, plus an
unfactored serviceability case for the deflection check. Strength combinations
are **enveloped** member by member.

### 2.5 Member design checks (AISC 360 LRFD)

Tension, yielding on the gross section: φt*P*n = 0.90 *F*y *A*g.

Compression, flexural buckling of a doubly symmetric member (Chapter E3):

> *F*e = π²*E*/(*KL*/*r*)²
> *F*cr = 0.658^(*F*y/*F*e) · *F*y   if *KL*/*r* ≤ 4.71√(*E*/*F*y)
> *F*cr = 0.877 *F*e                 otherwise
> φc*P*n = 0.90 *F*cr *A*g

Round hollow sections are doubly symmetric, so there is one slenderness ratio per
strut and no weak axis to track. Sections with *D*/*t* > 0.11*E*/*F*y are
locally slender and receive the Chapter E7 reduction *Q* = 0.038*E*/(*F*y *D*/*t*)
+ 2/3, which is applied rather than merely warned about.

### 2.6 Linear buckling

A single-layer dome does not fail by individual strut buckling; it fails by
*system* instability — snap-through of a node or of a whole cap. That is captured
by the linearised buckling problem. For a bar carrying axial force *N*, the
consistent geometric stiffness is

> **k**g = (*N*/*L*) · ⎡  **I** − **cc**ᵀ   −(**I** − **cc**ᵀ) ⎤
>                      ⎣ −(**I** − **cc**ᵀ)  **I** − **cc**ᵀ  ⎦

which represents the transverse force developed when a member carrying axial load
moves sideways. Assembling **K**g from the member forces of a reference load case
gives

> [**K** + λ**K**g]{φ} = {0}

Rather than solving this unsymmetric-looking form directly, it is rearranged as

> (−**K**g){φ} = (1/λ) **K** {φ}

a symmetric generalised eigenproblem whose right-hand matrix **K**ff is positive
definite. The *largest* eigenvalue μ of that pair gives the *smallest* critical
load factor λcr = 1/μ, which is numerically the well-conditioned end of the
spectrum. Problems below 2 400 free DOF are solved densely with
`scipy.linalg.eigh`; larger ones use shift-invert Lanczos.

λcr is the multiple of the *applied* load at which the dome becomes unstable.

### 2.7 Fully stressed auto-sizing

Geodesic domes are fabricated from a handful of repeated strut lengths, so the
design variables are the **families** (typically six to nine) rather than the
hundreds of individual members. The loop is:

```
repeat
    analyse the current design, envelope all strength combinations
    for each family:
        take the envelope tension and compression demands
        choose the lightest catalogue section whose LRFD capacity
          covers them at the target utilisation, subject to KL/r <= 200
    if nothing changed and the deflection limit is met: stop
    if nothing changed but the deflection limit is violated:
        raise the floor of the catalogue search and repeat
```

This converges in two or three iterations because the force distribution in a
triangulated dome is almost independent of the member sizing — the structure is
close to statically determinate in its load path. The result is a local optimum
(classical fully stressed design), not a global one, and the report says so.

---

## 3. Software

### 3.1 Architecture

```
geodome/
  geometry.py     Class I generation, truncation, boundary detection, chord factors
  sections.py     CHS catalogue; A, I, r, S computed from D and t
  materials.py    steel and aluminium grades
  model.py        DomeModel container; JSON save/load
  loads.py        dead / snow / live / wind vectors; LRFD combination
  solver.py       direction cosines, K and Kg assembly, partitioned solve
  buckling.py     AISC capacities; linear buckling eigenanalysis
  optimize.py     fully stressed auto-sizing
  postprocess.py  multi-case envelope, summary, result tables
  plotting.py     matplotlib 3D renderers (no pyplot anywhere)
  report.py       CSV / PNG / self-contained HTML export
  verification.py 75-check benchmark suite
  gui/            canvas, panels, background worker, main window, stylesheet
```

The analysis layer has no knowledge of Qt and the GUI layer contains no
mechanics, so the whole program also runs headless
(`run_geodome.py --batch 6 15 0.5 out/`) for parameter studies.

Two deliberate choices are worth noting. First, the GUI is written against
`qtpy`, so the identical source runs on PyQt5, PyQt6 or PySide6 — whichever the
marker's Anaconda happens to provide. Second, `plotting.py` never imports
`pyplot`; figures are built from `matplotlib.figure.Figure` directly, so Spyder's
inline backend cannot capture or close the embedded canvases.

### 3.2 Interface

![GUI](../output/screenshots/gui_util.png)

The window is divided into an input dock (Geometry, Sections, Loads, Design), a
set of 3D views (Model, Loads, Deflected shape, Axial force, Utilisation,
Buckling mode, Charts) and a result area (Summary, Members, Nodes, Reactions,
Auto-sizing, Verification, Log). Analysis, sizing and verification all run on a
background thread so the window stays responsive and the progress bar moves.

### 3.3 Units

Internally everything is SI (m, N, Pa, kg). Conversion happens only at the
interface: pressures are entered in kPa, point loads in kN, and results are shown
in mm, kN and MPa.

---

## 4. Verification

The suite (`run_geodome.py --verify`, or F11 in the GUI) runs **75 checks** and
all 75 pass. The point of each is that the reference value does *not* come from
the program.

### 4.1 Closed-form solutions

| check | computed | reference | rel. error |
|---|---|---|---|
| Single bar, tip displacement | 3.750000 × 10⁻⁴ m | *PL*/*EA* | 0 |
| Single bar, axial force | 50 000 N | *P* | 0 |
| Tripod, leg axial force | −29 814.2 N | −*PL*/3*h* | 2.4 × 10⁻¹⁶ |
| Tripod, apex deflection | −3.726780 × 10⁻⁴ m | −*PL*³/3*EAh*² | 1.5 × 10⁻¹⁶ |
| Tripod, vertical reaction | 80 000 N | *P* | 0 |
| Two-bar node, bar forces | *P*x, −*P*z | determinate equilibrium | 0 |

### 4.2 Geometry and topology

Chord factors as tabulated in Section 2.1. Vertex, strut and panel counts and
Euler's formula at 1V, 2V, 3V, 4V and 6V, all exact. Maxwell's count
*m* = 3*n* − 6 exact for the closed sphere at 2V, 4V and 6V. The faceted
hemisphere area converges on 2π*R*² from below: 1.89 % low at 4V, 0.48 % low
at 8V.

### 4.3 Statics, independent of the stiffness solution

**Method of sections.** A horizontal plane is passed through the dome, the
resultant of the member forces crossing it is summed, and compared with the
resultant of the loads applied above it. This never touches the stiffness matrix
again — it is a pure free-body check on the computed member forces.

| cut height | residual (normalised) |
|---|---|
| *z* = 2.50 m | 2.8 × 10⁻¹⁶ |
| *z* = 5.00 m | 6.5 × 10⁻¹⁶ |
| *z* = 7.50 m | 1.0 × 10⁻¹⁵ |

### 4.4 Comparison with membrane shell theory

For a spherical dome carrying a uniform load *p* per unit **projected** area,
membrane theory gives a meridional stress resultant that is *constant over the
whole shell*: *n*φ = −*pR*/2. Extracting it from the 8V model by section cut
(*n*φ = |ΣN·c_z| / (2π ρ sin φ), with ρ = *R* sin φ):

| cut | computed *n*φ | theory *pR*/2 | error |
|---|---|---|---|
| *z*/*R* = 0.30 | 4 869 N/m | 5 000 N/m | 2.6 % |
| *z*/*R* = 0.50 | 5 105 N/m | 5 000 N/m | 2.1 % |
| *z*/*R* = 0.70 | 4 899 N/m | 5 000 N/m | 2.0 % |

The residual 2 % is the discretisation: a smooth pressure field is being lumped
onto discrete nodes.

**Hoop behaviour.** The classic dome result — hoop compression near the crown
reversing to hoop tension near the springing — is reproduced. The reversal is
asserted; its *location* is not, and the reason is worth recording. Smooth-shell
theory puts the reversal at φ = 45° for projected load, but the model gives
≈ 35°. That is not an error: in a *triangulated* lattice the diagonal struts also
carry hoop force, so taking only the struts that lie along a parallel circle
under-represents the hoop resultant and shifts the apparent reversal up the
meridian. Measured at 6V, 8V and 10V the apparent reversal converges on 35°, not
on 45°, confirming it is a property of the extraction rather than a mesh error.

### 4.5 Internal identities

For every load combination of a 4V dome:

| identity | worst value |
|---|---|
| ‖**Ku** − **F** − **R**‖ / ‖**F**‖ | 6.4 × 10⁻¹⁵ |
| ½**F**ᵀ**u** vs ½Σ*N*²*L*/*EA* | 1.1 × 10⁻¹⁵ |
| Σ**F** + Σ**R** | 3.5 × 10⁻¹⁰ N |
| five-fold symmetry of *u*z | 6.9 × 10⁻¹⁵ |
| buckling: λ(2**N**) vs λ(**N**)/2 | 0 |

### 4.6 Unit tests

A separate 30-test suite (`tests/test_geodome.py`) covers the catalogue ordering,
model JSON round-tripping, load linearity, the AISC branch switch, the
monotonicity of capacity with length, the fact that buckling mode shapes are zero
at supported DOF, that the utilisation envelope bounds every individual case, and
that **every truncation in the usable range leaves no unsupported node with fewer
than three struts**. All 30 pass.

---

## 5. Case study — 30 m span exhibition hall

**Structure.** 6V Class I hemisphere, *R* = 15 m, 30 m span, 15 m rise.
196 nodes, 555 struts, 360 panels, 9 strut families, 30 pinned base nodes.
ASTM A500 Gr. C (*E* = 200 GPa, *F*y = 345 MPa).
Faceted surface 1 401.7 m², footprint 701.7 m².

**Loading.** 0.35 kPa superimposed dead on the surface; 0.60 kPa snow on the
projection; 58 m/s basic wind (*K*z = 0.95, *K*d = 0.95, *G* = 0.85, giving
*q*z = 1.861 kPa). Unfactored resultants:

| case | *F*x [kN] | *F*z [kN] |
|---|---|---|
| Dead (self weight + cladding) | 0 | −710.7 |
| Snow | 0 | −421.0 |
| Wind | +263.9 | **+754.3** |

The wind resultant is a net **uplift of 754 kN** — larger than the entire snow
load and comparable with the dead load. This is the single most important thing
the load generator shows, and it is what drives the result below.

### 5.1 First pass — trial section CHS 114.3 × 5.0

| combination | max T [kN] | max C [kN] | max \|u\| [mm] | max util |
|---|---|---|---|---|
| 1.4D | 27.3 | −19.4 | 2.24 | 0.058 |
| 1.2D + 1.6S | 48.1 | −32.5 | 4.58 | **0.094** |
| 1.2D + 1.0W + 0.5S | 45.4 | −30.9 | 2.34 | 0.092 |
| 0.9D + 1.0W | 32.3 | −24.4 | 2.06 | 0.077 |
| 1.0D + 1.0S (SLS) | 35.0 | −23.3 | 3.26 | 0.069 |

Governing member #416 (family G), utilisation 0.094, *KL*/*r* = 79.8.
Crown deflection 3.26 mm = span/9192 against a span/240 limit.
Mass 22 445 kg (32.0 kg/m² of plan). λcr = 56.6.

The trial design is enormously over-strong — utilisation under 10 %. This is
typical of a dome: the surface is so efficient in membrane action that strength
is almost never the binding constraint.

### 5.2 Auto-sizing

| iteration | mass [kg] | max util | deflection | sections |
|---|---|---|---|---|
| 1 | 22 445 | 0.094 | span/9192 | CHS 114.3 × 5.0 (all) |
| 2 | 7 504 | 0.808 | span/3428 | CHS 60.3 × 3.2 (all) |
| 3 | 7 504 | 0.808 | span/3428 | converged |

**22 445 kg → 7 504 kg, a 66.6 % saving**, in two iterations. Mass per unit plan
area falls from 32.0 to 10.7 kg/m².

### 5.3 Final design

| combination | max T [kN] | max C [kN] | max \|u\| [mm] | max util |
|---|---|---|---|---|
| 1.4D | 21.8 | −15.5 | 5.28 | 0.420 |
| 1.2D + 1.6S | 43.4 | −29.6 | 12.49 | 0.751 |
| 1.2D + 1.0W + 0.5S | 40.8 | −27.8 | 6.41 | 0.724 |
| **0.9D + 1.0W** | 33.8 | −28.0 | 6.57 | **0.808** |
| 1.0D + 1.0S (SLS) | 31.0 | −20.8 | 8.75 | 0.544 |

Governing member #452 (family I, CHS 60.3 × 3.2), *KL*/*r* = 160.7,
φc*P*n = 34.6 kN against φt*P*n = 178.2 kN. Deflection span/3428. λcr = 41.9.

**The governing combination has moved from snow to wind uplift.** The mechanism
is clear once stated: 0.9D + 1.0W is the combination that *minimises* the
stabilising dead load while applying full uplift, and lightening the structure by
two thirds removed most of the dead load that had been holding the dome down.
Optimising for the original governing case changed which case governs. Enveloping
every combination at every iteration is what catches this; checking only the
snow case would have produced a design that is 8 % overstressed under wind.

Note also how much the capacity ratio has shifted: the governing strut's
compression capacity (34.6 kN) is now only 19 % of its tension capacity
(178.2 kN), because *KL*/*r* = 161 puts it deep in the elastic buckling branch.
The design is buckling-governed, not stress-governed.

### 5.4 Stability

λcr = 41.9 on the governing combination. The mode shape is a **local
snap-through of the crown cap** — 100 % of the modal strain energy lies above
0.8 of the dome height, with the peak at the crown node itself. It is not a
global collapse mode and it is not an individual strut buckling.

Eigenvalues come in near-degenerate pairs (64.961 twice, 80.559 twice in the 4V
model) which is the expected signature of the five-fold rotational symmetry of an
icosahedral dome, and is itself a check on the eigensolver.

λcr = 41.9 looks like a large margin, but single-layer domes are notoriously
**imperfection sensitive**: a linear eigenvalue overestimates the real capacity,
and a knock-down factor of 0.3–0.5 on λcr is normal practice, giving an effective
margin nearer 13–21. The program prints this caution whenever λcr < 2 and reports
it in the summary regardless.

### 5.5 Parameter study — truncation height

6V, *R* = 15 m, CHS 114.3 × 5.0 throughout:

| fraction | rise [m] | span [m] | mass [kg] | max util | λcr | deflection | mode at *z*/*H* |
|---|---|---|---|---|---|---|---|
| 0.317 | 9.50 | 27.91 | 12 499 | 0.109 | 56.6 | span/14645 | 1.00 |
| 0.406 | 12.19 | 29.47 | 15 405 | 0.087 | 56.5 | span/13232 | 1.00 |
| 0.500 | 15.00 | 30.00 | 22 445 | 0.094 | 56.6 | span/9192 | 1.00 |
| 0.594 | 17.81 | 29.55 | 26 210 | 0.119 | 107.6 | span/7462 | 0.32 |
| 0.683 | 20.50 | 28.38 | 29 850 | 0.138 | 90.0 | span/5897 | 0.41 |

Two findings:

1. **Deeper domes are heavier and, relative to span, more flexible.** Mass rises
   2.39× from the shallow cap to the 0.68 dome because there is more surface to
   clad and carry, while deflection worsens from span/14645 to span/5897. Past
   the hemisphere the base ring stops growing, so added surface brings load on
   faster than added rise brings stiffness.
2. **The buckling behaviour splits into two regimes.** Up to the hemisphere the
   critical mode is the crown snap-through and λcr is essentially constant at
   56.5–56.6 — because the crown cap is *geometrically identical* in all of
   them, so the same local instability governs regardless of how far down the
   dome is cut. Beyond the hemisphere the governing combination switches to wind
   and the critical mode migrates down to the windward flank (*z*/*H* ≈ 0.3–0.4).
   Truncation height is a cost decision below the hemisphere and a wind-stability
   decision above it.

---

## 6. Discussion

### 6.1 What the matrix method made visible

The axial force plot (Figure 04) shows membrane behaviour directly: every
meridional strut in compression, and a band of hoop **tension** in the lower
rings, shown in red. A dome pushes outwards at its base, and the lower rings are
the tie that resists it. This is the result that membrane shell theory predicts
analytically, arrived at here purely from assembling and solving **K**.

### 6.2 Where the idealisation matters

* **Pinned joints.** Real hubs have some rotational stiffness, so the model
  underestimates stiffness and, more importantly, under-predicts the local
  buckling resistance of the crown cap. The computed λcr is conservative in that
  respect and unconservative with respect to imperfections.
* **Linear buckling.** The eigenvalue gives the bifurcation load of the perfect
  structure. Snap-through of a shallow cap is a *limit point*, not a
  bifurcation, and a full nonlinear incremental analysis would give a lower
  value. This is the single largest modelling limitation, and it is why the
  knock-down factor is recommended rather than optional.
* **Wind.** The axisymmetric *C*p distribution is a simplification. It captures
  net uplift correctly, which is what governs, but it will not give a correct
  asymmetric load for a torsional check.

### 6.3 A bug worth reporting

An early version placed supports at "every node whose *z* equals the cut
elevation". For a hemisphere of even frequency that is correct — the equator *is*
a complete ring. For every other truncation it is wrong: the cut elevation
contains only the 5- or 10-node sub-ring inherited from the icosahedron's own
vertex latitudes, and the actual free edge of the retained mesh sits above it.
The symptom was a parameter study in which the hemisphere behaved perfectly and
every other dome fraction returned displacements three orders of magnitude too
large, utilisations of 86, and a failed eigensolve.

What made it diagnosable was not the displacement magnitude but a connectivity
count: nodes with fewer than three struts and a support set of 5 nodes out of
141. The fix was to define the support set topologically — nodes on an edge
belonging to exactly one panel — which is correct for any cut. The unit test
`test_every_truncation_is_stable` now checks 15 combinations of frequency and
truncation and asserts that no unsupported node has fewer than three struts.

The episode is worth recording because the failure was *silent in the
hemisphere case*, which is exactly the case anyone would test first.

---

## 7. Conclusions

1. The direct stiffness method extends to three-dimensional space trusses with
   one extra direction cosine, and handles a 555-member geodesic dome without
   difficulty. Equilibrium residuals are of order 10⁻¹⁵ and the strain-energy
   identity closes to the same order.
2. The implementation reproduces published geodesic chord factors to better than
   2 × 10⁻⁶ and membrane shell theory to within 2.6 %, so the geometry and the
   mechanics are both independently validated.
3. The closed geodesic sphere is exactly isostatic by Maxwell's count at every
   frequency; a truncated dome is not, and relies on its supports. Confirming
   that **K**ff is positive definite — not merely that the solver returns an
   answer — is a necessary part of validating a space-truss program.
4. Geodesic domes are strength-governed only at absurdly conservative sizing.
   Fully stressed design cut the case-study dome by 66.6 % and left it governed by
   **compression buckling of slender struts under wind uplift**, not by stress.
5. Optimisation can change which load combination governs. In the case study it
   moved from 1.2D + 1.6S to 0.9D + 1.0W, because removing dead load removes the
   restraint against uplift. Enveloping all combinations inside the optimisation
   loop is essential, not a refinement.
6. System stability, not member strength, is the real limit on a single-layer
   dome. The critical mode is a local crown snap-through, and its load factor is
   insensitive to truncation height below the hemisphere.

### Suggested further work

* Geometric nonlinear (incremental Newton–Raphson) analysis to trace the
  snap-through limit point rather than estimate it from an eigenvalue.
* Imperfection sensitivity: seed the model with a scaled first buckling mode and
  re-run, to quantify the knock-down factor instead of quoting a rule of thumb.
* Class II and III breakdowns, and elliptical/non-spherical domes.
* Rigid-jointed (space frame) elements, to bound the benefit of moment-resisting
  hubs.

---

## 8. References

1. AISC, *Specification for Structural Steel Buildings*, ANSI/AISC 360.
   Chapter D (tension), Chapter E3 (flexural buckling), Chapter E7 (slender
   elements).
2. ASCE/SEI 7, *Minimum Design Loads and Associated Criteria for Buildings and
   Other Structures*. Velocity pressure; Figure 27.3-2, domed roofs; Section 2.3
   LRFD load combinations.
3. Kenner, H., *Geodesic Math and How to Use It*, University of California Press
   — Class I breakdown and chord factor tables.
4. Timoshenko, S. and Woinowsky-Krieger, S., *Theory of Plates and Shells* —
   membrane theory of spherical shells.
5. Maxwell, J. C., "On the calculation of the equilibrium and stiffness of
   frames", *Philosophical Magazine*, 1864.
6. Course notes, Matrix Structural Analysis — direct stiffness method, geometric
   stiffness and linearised buckling.

---

## 9. Acknowledgment of AI assistance

See `docs/ai_acknowledgment.md` for the full disclosure required by the course AI
policy.

---

## Appendix A — Verification suite output

Reproduce with `python run_geodome.py --verify`. 75 of 75 checks pass.

## Appendix B — Source code

Listed in `README.md` §4. 5 872 lines of Python in total:

| component | files | lines |
|---|---|---|
| `geodome/` — analysis | 14 | 3 024 |
| `geodome/gui/` — interface | 6 | 1 792 |
| `tests/test_geodome.py` | 1 | 379 |
| `examples/worked_example.py` | 1 | 180 |
| `run_geodome.py` | 1 | 83 |
| *(verification suite, counted inside `geodome/`)* | 1 | 543 |

The case study is reproducible in full with `python run_geodome.py --demo`.
