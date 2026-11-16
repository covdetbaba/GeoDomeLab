# GeoDomeLab theory and method

This page summarises the methods in GeoDomeLab 2.0 and the equations behind
them. The module docstrings in `geodome/` give the same equations next to the
code. Help > Theory and method in the window shows a shorter version.

Books are cited by chapter. Section and equation numbers are given only where
they were checked against the source.

## Contents

1. [Conventions](#1-conventions)
2. [Geometry](#2-geometry)
3. [Direct stiffness method](#3-direct-stiffness-method)
4. [Loads and combinations](#4-loads-and-combinations)
5. [Member checks to AISC 360-16](#5-member-checks-to-aisc-360-16)
6. [Linear buckling](#6-linear-buckling)
7. [Geometrically nonlinear analysis](#7-geometrically-nonlinear-analysis)
8. [Imperfection sensitivity](#8-imperfection-sensitivity)
9. [Modal analysis](#9-modal-analysis)
10. [Fabrication geometry](#10-fabrication-geometry)
11. [Auto-sizing and parameter studies](#11-auto-sizing-and-parameter-studies)
12. [Independent check with OpenSees](#12-independent-check-with-opensees)
13. [References](#13-references)

## 1. Conventions

- SI units inside the program: m, N, Pa, kg, °C.
- Tension is positive. The z axis points up.
- Three translations per node. The DOF numbers of node $n$ are $3n$ ($x$),
  $3n+1$ ($y$) and $3n+2$ ($z$).
- Every strut is a straight pin-ended bar with axial force only.

## 2. Geometry

GeoDomeLab builds Class I icosahedral geodesic domes. Each of the 20 faces of
an icosahedron is split into $\nu^2$ triangles, where $\nu$ is the frequency.
Every new point is projected radially onto the sphere of radius $R$. A full
sphere has

$$
V = 10\nu^2 + 2, \qquad E = 30\nu^2, \qquad F = 20\nu^2, \qquad V - E + F = 2 .
$$

Twelve nodes have 5 struts and the rest have 6. The closed sphere is exactly
isostatic: $E = 3V - 6$ (Maxwell's count).

The dome is cut on a horizontal plane at a ring of existing nodes, so the
base is always a closed polygon. Supports go on the free edge of the cut mesh.
Struts of equal length (within a tolerance) form one strut family. The chord
factor of a family is

$$
\mathrm{cf} = \frac{L}{R} .
$$

For 2V the closed forms are $\mathrm{cf} = 2\sin 18^\circ = 0.618034$ and
$\mathrm{cf} = 2\sin\!\left(\tfrac14 \arccos\tfrac{1}{\sqrt5}\right) = 0.546533$
(Kenner 1976 prints 0.61803 and 0.54653, as reproduced by the Domerama
chord-factor chart). $2\sin 18^\circ = 1/\varphi$ is the edge of the
icosidodecahedron formed by the icosahedron's edge midpoints.

Any other space truss can be read from a model file or from CSV files. Its
coordinates are used as given. Its "radius" is the best-fit sphere radius
when the nodes lie on a sphere, and it is undefined otherwise.

## 3. Direct stiffness method

For a bar from node $i$ to node $j$ with length $L$, area $A$, modulus $E$ and
unit vector $\mathbf c = (\mathbf x_j - \mathbf x_i)/L$, the element
stiffness in global axes is

$$
\mathbf k = \frac{EA}{L}
\begin{bmatrix} \mathbf c\mathbf c^T & -\mathbf c\mathbf c^T \\
-\mathbf c\mathbf c^T & \mathbf c\mathbf c^T \end{bmatrix} .
$$

The element matrices are added into the global stiffness $\mathbf K$
(sparse). The DOFs are split into free ($f$) and supported ($s$) sets
(McGuire, Gallagher & Ziemian 2000, Ch. 3):

$$
\mathbf K_{ff}\,\mathbf u_f = \mathbf F_f + \mathbf F_{0f} - \mathbf K_{fs}\,\mathbf u_s ,
$$

$$
\mathbf R = \mathbf K_{sf}\,\mathbf u_f + \mathbf K_{ss}\,\mathbf u_s - \mathbf F_s - \mathbf F_{0s} .
$$

- $\mathbf u_s$ holds prescribed support movements (settlements).
- An initial (free) strain $\varepsilon_0$, for example $\alpha\,\Delta T$,
  gives the equivalent load
  $\mathbf F_0 = \sum EA\,\varepsilon_0 \,[-\mathbf c;\ +\mathbf c]$.
- The axial force is
  $N = EA\left(\dfrac{e}{L} - \varepsilon_0\right)$ with the elongation
  $e = \mathbf c^T(\mathbf u_j - \mathbf u_i)$.
- An elastic support is a free DOF with its spring stiffness $k$ added to
  the diagonal of $\mathbf K_{ff}$. Its spring force $-k\,u$ is reported as a
  reaction.

Every solve checks itself:

- the equilibrium residual $\lVert \mathbf K_{ff}\mathbf u_f - \mathbf F^{*}_f\rVert / \lVert \mathbf F^{*}_f\rVert$,
  where $\mathbf F^{*}_f$ is the right-hand side of the first equation above,
- the balance $\sum \mathbf F + \sum \mathbf R = \mathbf 0$,
- the strain-energy identity. Without initial strains and settlements it is
  $\tfrac12 \mathbf F^T \mathbf u = \sum \dfrac{N^2 L}{2EA} + \sum \tfrac12 k u^2$.
  The solver uses a general form that also holds with them.

## 4. Loads and combinations

**Dead load D.** Strut self weight $\rho A L g$, half to each end node, plus
a cladding pressure on the faceted surface area.

**Snow S and live load L.** A pressure on the horizontally projected area
of the panels that face up: a panel takes $p\,A\max(n_z, 0)$ with its
outward normal $\mathbf n$, so the panels below the equator of a deep dome
carry none. Unbalanced snow puts $1.5p$ on the half $x \ge 0$ and $0.5p$ on
the other half. This is a simple test pattern, not the ASCE 7-16 curved-roof
rule (Fig. 7.6-2). A crown point load can be added to L.

**Wind W.** A pressure $p = q_z\,G\,C_p$ normal to the surface, with the
velocity pressure

$$
q_z = 0.613\,K_z K_{zt} K_d V^2 \quad [\mathrm{Pa},\ V \text{ in m/s}] .
$$

$C_p$ is interpolated along the meridian angle $\theta$, measured from the
windward base ($0^\circ$) over the crown ($90^\circ$) to the leeward base
($180^\circ$):

| $\theta$ | 0° | 30° | 60° | 90° | 120° | 150° | 180° |
|---|---|---|---|---|---|---|---|
| $C_p$ | +0.80 | +0.20 | −0.60 | −1.00 | −0.80 | −0.50 | −0.40 |

This 7-point curve is an illustrative simplification that follows the
trend of the ASCE 7 domed-roof chart; it is not read from Fig. 27.3-2. The
wind direction is the direction the wind blows towards, measured from $+x$
(0 means wind from $-x$, so the windward base is at $x = -R$). Each node $i$
takes $-p_i A_i \mathbf n_i$ with its surface tributary area $A_i$ and its
radial normal $\mathbf n_i$. Under a uniform pressure this gives a little
less than the exact resultant of the faceted surface, $p$ times the plan
area: 5.1 % less at 2V, 1.4 % at 4V, 0.6 % at 6V and 0.35 % at 8V.

**Tributary areas.** Area loads go to the nodes by tributary area: one third
of each triangular panel that meets the node.

**Nodal loads** can be added to any of the cases D, S, L and W.

**Self-straining load T.** A uniform temperature change $\Delta T$ and support
settlements $\boldsymbol\delta$ (ASCE 7 notation). T is not a nodal load. A
combination with factor $\gamma_T$ uses

$$
\varepsilon_0 = \gamma_T\,\alpha\,\Delta T, \qquad \mathbf u_s = \gamma_T\,\boldsymbol\delta .
$$

A statically determinate truss gets no member force from T. A dome with
redundant members does. ASCE 7-16 Sec. 2.3.4 covers combinations with T and
asks for a load factor of at least 1.0 on T. The thermal coefficients are
$11.7\times10^{-6}$/°C for ASTM steels (the value conventionally used for
structural steel at ambient temperature), $12\times10^{-6}$/°C for EN steels
(EN 1993-1-1, 3.2.6) and $23\times10^{-6}$/°C for aluminium
(EN 1999-1-1, 3.2.5).

**Combinations.** The default strength combinations follow ASCE 7-16,
Sec. 2.3.1, with S standing for roof live, snow or rain:

| Name | Combination | ASCE 7-16 Sec. 2.3.1 |
|---|---|---|
| LC1 | 1.4D | 1 |
| LC2 | 1.2D + 1.6S | 3 |
| LC3 | 1.2D + 1.0W + 0.5S | 4 |
| LC4 | 0.9D + 1.0W | 5 |
| LC5 | 1.2D + 1.6S + 0.5W | 3, with 0.5W |
| LC6 | 1.2D + 1.6L + 0.5S | 2 |
| SLS | 1.0D + 1.0S (serviceability, for the deflection check) | |

When the live load L is not zero, ASCE 7 also adds 1.0L to combinations 3
and 4; add those rows by hand. A run names any load case (D, S, L or W) that
is defined but used by no strength combination. Every combination is
solved. The member checks use the envelope of the
strength combinations. The deflection check compares the span with the
largest nodal displacement under the serviceability combinations.

## 5. Member checks to AISC 360-16

Each strut is a round HSS (CHS) with gross area $A_g$, radius of gyration $r$,
effective length $KL$ and LRFD factors $\phi_t = \phi_c = 0.90$ (AISC 2016).

Tension, yielding (D2):

$$
\phi_t P_n = \phi_t F_y A_g .
$$

Compression, flexural buckling (E3):

$$
F_e = \frac{\pi^2 E}{(KL/r)^2}, \qquad
F_{cr} = \begin{cases}
0.658^{F_y/F_e}\,F_y & KL/r \le 4.71\sqrt{E/F_y} \\
0.877\,F_e & \text{otherwise}
\end{cases}
$$

Slender round HSS (E7), with the effective area $A_e$:

$$
A_e = \begin{cases}
A_g & D/t \le 0.11\,E/F_y \\
\left[\dfrac{0.038\,E}{F_y\,(D/t)} + \dfrac23\right] A_g & 0.11\,E/F_y < D/t < 0.45\,E/F_y
\end{cases}
\qquad \phi_c P_n = \phi_c F_{cr} A_e .
$$

Above $D/t = 0.45\,E/F_y$ the specification does not apply, and the program
sets the capacity to zero. The utilisation of a strut is

$$
U = \max\left(\frac{N_t}{\phi_t P_{n,t}},\ \frac{|N_c|}{\phi_c P_{n,c}}\right)
$$

over the strength envelope. The design wall thickness factor (B4.2, 0.93 for
ERW HSS such as ASTM A500) scales $t$ in the strength checks only, not in the
stiffness. Tension rupture at the connection is not checked, because the hub
detail is not modelled. AISC 360 is a steel specification: for aluminium the
program still gives forces and displacements, and every run says that the
member checks do not apply.

## 6. Linear buckling

The geometric stiffness of a bar with axial force $N$ is

$$
\mathbf k_g = \frac{N}{L}
\begin{bmatrix} \mathbf I - \mathbf c\mathbf c^T & -(\mathbf I - \mathbf c\mathbf c^T) \\
-(\mathbf I - \mathbf c\mathbf c^T) & \mathbf I - \mathbf c\mathbf c^T \end{bmatrix} .
$$

With $\mathbf K_g$ assembled from the member forces $N_{\mathrm{ref}}$ of the
governing combination, the critical load factor $\lambda_{cr}$ solves

$$
\left[\mathbf K + \lambda\,\mathbf K_g(N_{\mathrm{ref}})\right]\boldsymbol\phi = \mathbf 0
$$

(McGuire et al. 2000, Ch. 9). It is solved as the symmetric problem
$(-\mathbf K_g)\boldsymbol\phi = \mu\,\mathbf K\boldsymbol\phi$ with
$\lambda = 1/\mu$. The largest $\mu$ gives the smallest positive
$\lambda_{cr}$. Problems up to 2400 free DOFs use a dense solver, larger ones
a sparse shift-invert solver. Springs are part of $\mathbf K$.

$\lambda_{cr}$ is the bifurcation load factor of the perfect structure. A
shallow dome can snap through well below it, so GeoDomeLab also traces the
nonlinear path.

## 7. Geometrically nonlinear analysis

**Corotational bar.** Displacements can be large. Strains stay small and the
material stays elastic. With the current bar vector
$\mathbf d = \mathbf d_0 + \mathbf u_j - \mathbf u_i$, current length
$l = \lVert\mathbf d\rVert$ and unit vector $\mathbf c = \mathbf d/l$
(Crisfield 1991, Ch. 3):

$$
\varepsilon = \frac{l - L}{L}, \qquad N = EA\,\varepsilon, \qquad
\mathbf f = N\,[-\mathbf c;\ +\mathbf c] ,
$$

$$
\mathbf k_T = \begin{bmatrix} \mathbf k_3 & -\mathbf k_3 \\ -\mathbf k_3 & \mathbf k_3 \end{bmatrix},
\qquad
\mathbf k_3 = \frac{EA}{L}\,\mathbf c\mathbf c^T + \frac{N}{l}\left(\mathbf I - \mathbf c\mathbf c^T\right) .
$$

This is engineering strain on the original length, the same measure and
tangent as the OpenSees `corotTruss` element. The tests check the tangent
against finite differences.

**Equilibrium.** The load is $\lambda\,\mathbf F_{\mathrm{ref}}$, where
$\mathbf F_{\mathrm{ref}}$ is the D, S, L and W part of one combination. The
self-straining load T is left out. Each iteration forms the tangent again
(full Newton-Raphson) and solves

$$
\mathbf r = \lambda\,\mathbf F_{\mathrm{ref}} - \mathbf f(\mathbf u), \qquad
\mathbf K_T\,\mathbf v_r = \mathbf r, \qquad
\mathbf K_T\,\mathbf v_t = \mathbf F_{\mathrm{ref}}
$$

(McGuire et al. 2000, Ch. 12). Convergence needs
$\lVert\mathbf r\rVert / \lVert\lambda\mathbf F_{\mathrm{ref}}\rVert < 10^{-8}$
and a small displacement correction.

**Load control** fixes $\lambda$ in each step. It cannot pass a load maximum.

**Displacement control** fixes one DOF $c$ in each step and sets
$\delta\lambda = -v_{r,c}/v_{t,c}$. It passes a load maximum but fails where
the control displacement snaps back.

**Arc-length** (Riks 1979; Crisfield 1981; Crisfield 1991, Ch. 9). Each step
keeps the cylindrical constraint
$\Delta\mathbf u^T\Delta\mathbf u = \Delta l^2$. With
$\delta\mathbf u = \mathbf v_r + \delta\lambda\,\mathbf v_t$, the corrector
solves

$$
a\,\delta\lambda^2 + b\,\delta\lambda + c = 0, \qquad
a = \mathbf v_t\cdot\mathbf v_t, \quad
b = 2\,\mathbf v_t\cdot(\Delta\mathbf u + \mathbf v_r), \quad
c = \lVert\Delta\mathbf u + \mathbf v_r\rVert^2 - \Delta l^2 ,
$$

and keeps the root whose new increment is closest in direction to the step
so far. The predictor takes the sign of the previous increment. The arc
length adapts to the number of iterations, and a failed step is cut back.
This follows the path through load maxima and snap-backs.

**Stability along the path.** The program tracks the sign of
$\det\mathbf K_T$, the number of negative eigenvalues of $\mathbf K_T$ (up to
600 free DOFs) and the current stiffness parameter

$$
S_p = \frac{\mathbf F^T\mathbf K_0^{-1}\mathbf F}{\mathbf F^T\mathbf K_T^{-1}\mathbf F} .
$$

The first load maximum $\lambda_{\mathrm{limit}}$ is refined and classed as a
limit point or a bifurcation by the load projection of the critical mode.
Above 600 free DOFs only the sign of $\det\mathbf K_T$ is available. A pair
of equal eigenvalues crossing zero together (a double bifurcation, common in
symmetric domes) leaves that sign unchanged, so it can be missed; the run
says so, and $\lambda_{\mathrm{limit}}$ should then be compared with
$\lambda_{cr}$. A bracketed bifurcation is accurate to about $10^{-9}$
(relative). A first arc-length step that would go to a negative load is cut
back, and a path that dips below zero before its first maximum is not
reported as a limit point.

**Comparison paths.** The plot also shows the linear path and the linearised
second-order (P-Delta style) path

$$
\mathbf u = \lambda\,\mathbf K^{-1}\mathbf F_{\mathrm{ref}}, \qquad
\mathbf u = \lambda\left(\mathbf K + \lambda\,\mathbf K_{g,\mathrm{ref}}\right)^{-1}\mathbf F_{\mathrm{ref}} ,
$$

and reports $\lambda_{\mathrm{limit}}/\lambda_{cr}$. The program reports the
$\lambda$ at which the first strut reaches its Euler load and its design
strength. The path itself does not include strut buckling or yielding.

**Closed-form check.** For the von Mises two-bar truss with half span $b$,
rise $h$, $L = \sqrt{b^2 + h^2}$, apex deflection $w$, $y = h - w$ and
$l = \sqrt{b^2 + y^2}$, the exact engineering-strain path is

$$
P(w) = 2EA\,y\left(\frac1l - \frac1L\right),
$$

with the limit point at $l^{*} = (b^2 L)^{1/3}$ and the linearised critical
load $P_{cr} = 2EA\,h^3/(L\,b^2)$. These closed forms are derived here for
engineering strain; the problem itself is the one studied by Pecknold,
Ghaboussi & Healey (1985).

## 8. Imperfection sensitivity

Single-layer domes are sensitive to imperfections (Gioncu 1995). The
imperfection study moves the nodes to

$$
\mathbf X = \mathbf X_0 + a\,\boldsymbol\phi ,
$$

where $\boldsymbol\phi$ is a linear buckling mode scaled so that its largest
component is 1 and $a$ is the amplitude. Both signs of $a$ are run, and the
lower limit load is kept. The result is $\lambda_{\mathrm{limit}}/\lambda_{cr}$
against $a/\text{span}$. Bruno, Sassone & Venuti (2016) report amplitudes
from about span/10 000 to span/300 in use. They note that the Chinese code
JGJ 61-2003 limits construction imperfections to span/300, and that there is
no agreed rule for combining modes into the decisive shape (also Bulenda &
Knippers 2001). The dialog starts with span/1000, span/500, span/300 and
span/200.

**Local modes and the first step.** The first buckling mode of a geodesic
dome is often a dimple at a single five-valent hub, such as the crown. Then
the useful measure of $a$ is the local rise of that hub above its
neighbours, not the span. In the 6V, $R$ = 15 m case study the crown stands
0.198 m above its five neighbours, so span/300 (100 mm) and span/200
(150 mm) flatten most of it. The crown then snaps through at a very low
load, like a shallow von Mises truss. Such a snap can lie within a
millimetre of the start of the path. An imperfect run therefore starts with
a first arc length of $2\times10^{-5}$ times the length scale (0.3 mm for
$R$ = 15 m), and the step then grows to its usual size. A nominal first step
(about 0.3 m here) jumps over the snap and finds the later load maximum of
the inverted crown instead.

## 9. Modal analysis

Undamped free vibration (Chopra 2023, Ch. 10 and 13):

$$
\mathbf K\boldsymbol\phi_n = \omega_n^2\,\mathbf M\boldsymbol\phi_n, \qquad
f_n = \frac{\omega_n}{2\pi}, \qquad T_n = \frac{1}{f_n}, \qquad
\boldsymbol\phi_n^T\mathbf M\boldsymbol\phi_n = 1 .
$$

**Mass.** Lumped by default: half of each strut's mass $\rho A L$ at each
end, plus the cladding dead load $/g$, downward dead nodal loads $/g$, an
optional share of the snow load and any added masses. The same mass acts in
$x$, $y$ and $z$. The consistent bar mass is an option:

$$
\mathbf m = \frac{\rho A L}{6}\begin{bmatrix} 2\mathbf I & \mathbf I \\ \mathbf I & 2\mathbf I\end{bmatrix} .
$$

**Participation.** With $\boldsymbol\iota_d$ equal to 1 at every free DOF in
direction $d$ and 0 elsewhere:

$$
\Gamma_{nd} = \boldsymbol\phi_n^T\mathbf M\boldsymbol\iota_d, \qquad
M^{*}_{nd} = \Gamma_{nd}^2, \qquad
\sum_{n} M^{*}_{nd} = \text{free mass in direction } d .
$$

The running sum shows how many modes reach 90 % of the mass. With a
pre-stress the problem becomes
$(\mathbf K + \mathbf K_g)\boldsymbol\phi = \omega^2\mathbf M\boldsymbol\phi$:
compression lowers the frequencies, and the first one reaches zero at
$\lambda_{cr}$. Each result checks $\boldsymbol\Phi^T\mathbf M\boldsymbol\Phi = \mathbf I$,
$\boldsymbol\Phi^T\mathbf K\boldsymbol\Phi = \mathrm{diag}(\omega^2)$ and the
eigen residual.

## 10. Fabrication geometry

$$
\mathrm{cf} = \frac{L}{R}, \qquad
\alpha_{\mathrm{axial}} = \arccos\frac{L}{2R}, \qquad
\tau = 90^\circ - \alpha_{\mathrm{axial}} = \arcsin\frac{L}{2R}, \qquad
L_{\mathrm{cut}} = L - 2o .
$$

- The axial angle is the angle between a strut and the radius to its hub
  (Clinton 1971). The program computes it from the mesh vectors at both ends
  and checks it against the closed form. For 2V the axial angles are
  72.000° and 74.141°.
- $o$ is the hub offset, the distance from the hub centre to the tube end.
- The dihedral angle between two panels is $\delta = 180^\circ - \theta$,
  where $\theta$ is the angle between their outward normals (a convex fold).
  A re-entrant fold has $\delta = 180^\circ + \theta$.
- Hub types group the nodes by location (interior, edge, base) and by the
  order of the strut families around the hub. Panel types group the
  triangles by the families of their sides.

## 11. Auto-sizing and parameter studies

**Auto-sizing** is a fully stressed design on the strut families:

1. analyse the current design and envelope the strength combinations;
2. give each family the lightest catalogue section that carries its envelope
   tension and compression at the target utilisation, with $KL/r \le 200$,
   using the design wall (B4.2 factor) exactly as the final check does;
3. if the deflection limit fails, step every family up one size and keep
   that size as the family's lowest choice from then on;
4. repeat until nothing changes. The run counts as converged only when the
   final design meets the target utilisation and the deflection limit.

The force distribution of a triangulated dome hardly depends on the member
sizes, so the loop converges in a few iterations. The result is a good
design, not a proven optimum.

**A parameter study** rebuilds the dome for every combination of frequency,
radius and dome fraction, and runs the full analysis on each one (optionally
with auto-sizing). It reports steel mass, mass per plan area, largest
utilisation, $\lambda_{cr}$ and span/deflection.

## 12. Independent check with OpenSees

`geodome/opensees_bridge.py` writes the model as plain data (nodes, members,
E, A, supports, the factored load vector, thermal strains, settlements and
masses) and solves it again in OpenSees (McKenna, Scott & Fenves 2010)
through OpenSeesPy (Zhu, McKenna & Scott 2018): `Truss` and `corotTruss`
elements, `InitStrainMaterial` for temperature, `sp` constraints for
settlement, `zeroLength` springs, `eigen` for modes and displacement control
for the nonlinear path. OpenSees builds and solves everything else itself.
The results and the limits of this check are in
[VALIDATION.md](VALIDATION.md).

## 13. References

- AISC (2016). *Specification for Structural Steel Buildings*
  (ANSI/AISC 360-16). American Institute of Steel Construction, Chicago, IL.
- ASCE (2017). *Minimum Design Loads and Associated Criteria for Buildings and
  Other Structures* (ASCE/SEI 7-16). American Society of Civil Engineers,
  Reston, VA. doi:10.1061/9780784414248.
- Bruno, L., Sassone, M., Venuti, F. (2016). Effects of the equivalent
  geometric nodal imperfections on the stability of single layer grid shells.
  *Engineering Structures* 112, 184-199. doi:10.1016/j.engstruct.2016.01.017.
- Bulenda, T., Knippers, J. (2001). Stability of grid shells. *Computers &
  Structures* 79(12), 1161-1174. doi:10.1016/S0045-7949(01)00011-6.
- CEN (2005). EN 1993-1-1, *Eurocode 3: Design of steel structures, Part 1-1:
  General rules and rules for buildings*. Clause 3.2.6.
- CEN (2007). EN 1999-1-1, *Eurocode 9: Design of aluminium structures,
  Part 1-1: General structural rules*. Clause 3.2.5.
- Chopra, A. K. (2023). *Dynamics of Structures: Theory and Applications to
  Earthquake Engineering*, 6th ed. Pearson, Hoboken, NJ.
- Clinton, J. D. (1971). *Advanced Structural Geometry Studies, Part I:
  Polyhedral Subdivision Concepts for Structural Applications*. NASA CR-1734.
- Crisfield, M. A. (1981). A fast incremental/iterative solution procedure
  that handles "snap-through". *Computers & Structures* 13(1-3), 55-62.
  doi:10.1016/0045-7949(81)90108-5.
- Crisfield, M. A. (1991). *Non-linear Finite Element Analysis of Solids and
  Structures, Volume 1: Essentials*. Wiley, Chichester.
- Gioncu, V. (1995). Buckling of reticulated shells: state-of-the-art.
  *International Journal of Space Structures* 10(1), 1-46.
  doi:10.1177/026635119501000101.
- Greco, M., Venturini, W. S. (2006). Stability analysis of three-dimensional
  trusses. *Latin American Journal of Solids and Structures* 3(3), 325-344.
- Iu, C. K., Bradford, M. A. (2012). Higher-order non-linear analysis of steel
  structures. Part I: Elastic second-order formulation. *Advanced Steel
  Construction* 8(2), 168-182. doi:10.18057/IJASC.2012.8.2.5.
- Kenner, H. (1976). *Geodesic Math and How to Use It*. University of
  California Press, Berkeley.
- McGuire, W., Gallagher, R. H., Ziemian, R. D. (2000). *Matrix Structural
  Analysis*, 2nd ed. Wiley, New York.
- McKenna, F., Scott, M. H., Fenves, G. L. (2010). Nonlinear finite-element
  analysis software architecture using object composition. *Journal of
  Computing in Civil Engineering* 24(1), 95-107.
  doi:10.1061/(ASCE)CP.1943-5487.0000002.
- Pecknold, D. A., Ghaboussi, J., Healey, T. J. (1985). Snap-through and
  bifurcation in a simple structure. *Journal of Engineering Mechanics*
  111(7), 909-922. doi:10.1061/(ASCE)0733-9399(1985)111:7(909).
- Riks, E. (1979). An incremental approach to the solution of snapping and
  buckling problems. *International Journal of Solids and Structures* 15(7),
  529-551. doi:10.1016/0020-7683(79)90081-7.
- Zhu, M., McKenna, F., Scott, M. H. (2018). OpenSeesPy: Python library for
  the OpenSees finite element framework. *SoftwareX* 7, 6-11.
  doi:10.1016/j.softx.2017.10.009.
