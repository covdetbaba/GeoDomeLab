"""
Help > Theory: a short guide to the methods used by GeoDomeLab, as HTML.

``THEORY_HTML`` is one self-contained HTML page (Qt rich text and any web
browser can show it).  ``theory_html()`` returns it; the main window shows it
in a read-only text browser.  The equations are plain text in ``<pre>``
blocks.  Every citation was checked against the source as far as it was
seen (abstracts are marked); books are cited by chapter, never by an unverified section or equation number.

This module imports nothing, so it loads without Qt or matplotlib.
"""

THEORY_HTML = """<!doctype html>
<html><head><meta charset='utf-8'><title>GeoDomeLab theory and method</title>
<style>
body { font-family: 'Segoe UI', Calibri, Arial, sans-serif; font-size: 10pt;
       color: #1c2733; background-color: #ffffff; }
h1 { font-size: 16pt; color: #123a5c; margin-bottom: 2px; }
h2 { font-size: 12.5pt; color: #123a5c; margin-top: 18px; margin-bottom: 4px; }
p, li { margin-top: 3px; margin-bottom: 5px; }
pre { font-family: Consolas, 'Courier New', monospace; font-size: 9.5pt;
      background-color: #f3f6f9; padding: 6px; margin-top: 4px;
      margin-bottom: 8px; }
.sub { color: #5b6b7a; }
.ref { margin-bottom: 5px; }
</style></head><body>

<h1>GeoDomeLab 2.0: theory and method</h1>
<p class='sub'>A short guide to the methods in the program. Units are SI
inside the program (m, N, Pa, kg). Tension is positive. DOF number =
3 &times; node + (0 for x, 1 for y, 2 for z).</p>

<p><a href='#geometry'>1 Geometry</a> |
<a href='#dsm'>2 Direct stiffness method</a> |
<a href='#loads'>3 Loads</a> |
<a href='#checks'>4 Member checks</a> |
<a href='#buckling'>5 Linear buckling</a> |
<a href='#nonlinear'>6 Nonlinear analysis</a> |
<a href='#modal'>7 Modal analysis</a> |
<a href='#fabrication'>8 Fabrication geometry</a> |
<a href='#sizing'>9 Auto-sizing and studies</a> |
<a href='#vv'>10 Verification and validation</a> |
<a href='#refs'>References</a></p>

<h2 id='geometry'><a name='geometry'></a>1 Geometry</h2>
<p>The dome is a Class I icosahedral geodesic sphere of frequency &nu;.
Each of the 20 faces of an icosahedron is split into &nu;&sup2; triangles,
and every new point is projected radially onto the sphere. A full sphere has
10&nu;&sup2; + 2 nodes, 30&nu;&sup2; struts and 20&nu;&sup2; panels, and
V &minus; E + F = 2. The dome is cut on a horizontal plane at a ring of
nodes. Struts of the same length form one family. A custom truss (any nodes
and members) can also be opened from a model file.</p>

<h2 id='dsm'><a name='dsm'></a>2 Direct stiffness method (3D truss)</h2>
<p>Each strut is a pin-ended bar with three translations per node. With the
unit vector c from node 1 to node 2 and the length L:</p>
<pre>
        EA  |  c c^T   -c c^T |
  k  =  --  |                 |          6 x 6, global axes
        L   | -c c^T    c c^T |
</pre>
<p>The element matrices are added into the global stiffness K. The DOFs are
split into free (f) and supported (s) sets:</p>
<pre>
  Kff uf = Ff + F0f - Kfs us                displacements
  R = Ksf uf + Kss us - Fs - F0s            support reactions
  N = EA (e/L - eps0),  e = c^T (u2 - u1)   axial force, tension +
  F0 = sum EA eps0 [ -c ; +c ]              load of an initial strain
</pre>
<p>us holds prescribed support movements (settlement). eps0 is an initial
strain, for example &alpha; &Delta;T. Elastic supports stay free and add
their spring stiffness k to the diagonal of Kff. Every solve reports three
checks: the equilibrium residual, the balance of loads and reactions, and the
strain-energy identity 1/2 F^T u = sum N^2 L / (2 EA), extended for initial
strains and settlements. Reference: McGuire, Gallagher and Ziemian (2000),
Ch. 3.</p>

<h2 id='loads'><a name='loads'></a>3 Loads</h2>
<ul>
<li>Dead load D: strut self weight (half to each end node) and a cladding
pressure on the faceted surface area.</li>
<li>Snow S and live load L: a pressure on the horizontal projected area
of the panels that face up (panels below the equator of a deep dome face
down and carry none). Snow can be unbalanced: 1.5 p on the half x &ge; 0 and
0.5 p on x &lt; 0. That is a simple test pattern, not the ASCE 7-16
curved-roof rule (Fig. 7.6-2).</li>
<li>Wind W: p = qz G Cp, normal to the surface, with
qz = 0.613 Kz Kzt Kd V^2 [Pa, V in m/s]. Cp follows a 7-point curve along
the meridian: +0.8 at the windward base, &minus;1.0 at the crown and
&minus;0.4 at the leeward base. This table is an illustrative simplification
that follows the trend of the ASCE 7 domed-roof chart; it is not read from
Fig. 27.3-2. The wind direction is the direction the wind blows towards,
measured from +x (0 = wind from &minus;x). Each node takes
&minus;p A n with its tributary area A and its radial normal n; under a
uniform pressure this gives 5.1 % less than the exact resultant (p times the
plan area) at 2V, 1.4 % at 4V and 0.6 % at 6V.</li>
<li>Area loads go to the nodes by tributary area: one third of each panel
that meets the node.</li>
<li>Nodal loads can be added to any case (D, S, L or W).</li>
</ul>
<p>The default strength combinations follow ASCE 7-16 Sec. 2.3.1, with S
standing for roof live, snow or rain:</p>
<pre>
  LC1   1.4D                       combination 1
  LC2   1.2D + 1.6S                combination 3
  LC3   1.2D + 1.0W + 0.5S         combination 4
  LC4   0.9D + 1.0W                combination 5
  LC5   1.2D + 1.6S + 0.5W         combination 3 with 0.5W
  LC6   1.2D + 1.6L + 0.5S         combination 2 (live load L)
  SLS   1.0D + 1.0S                serviceability (deflection check)
</pre>
<p>When the live load L is not zero, ASCE 7 also adds 1.0L to combinations 3
and 4 (1.2D + 1.6S + 1.0L, 1.2D + 1.0W + 1.0L + 0.5S); add those rows. The
member checks use the envelope of the strength combinations. A run names
any load case that no strength combination uses.</p>
<p><b>Self-straining load T.</b> T is a uniform temperature change plus any
support settlements (ASCE 7 notation). It is not a nodal load. The
temperature gives every strut the initial strain eps0 = &alpha; &Delta;T;
settlements are prescribed support displacements. The T column of a
combination scales both. ASCE 7-16 Sec. 2.3.4 covers combinations with T;
use a load factor of at least 1.0 on T. A statically determinate truss gets
no member force from T. A dome with redundant members does.</p>

<h2 id='checks'><a name='checks'></a>4 Member checks (AISC 360-16, LRFD)</h2>
<pre>
  Tension, yielding (D2)       &phi;_t Pn = 0.90 Fy Ag
  Compression (E3)             Fe  = &pi;^2 E / (KL/r)^2
                               Fcr = 0.658^(Fy/Fe) Fy   if KL/r &lt;= 4.71 sqrt(E/Fy)
                               Fcr = 0.877 Fe           otherwise
                               &phi;_c Pn = 0.90 Fcr Ae
  Slender round HSS (E7)       Ae = Ag                              D/t &lt;= 0.11 E/Fy
                               Ae = [0.038 E / (Fy D/t) + 2/3] Ag   0.11 E/Fy &lt; D/t &lt; 0.45 E/Fy
  Utilisation                  U = |N| / (&phi; Pn), envelope of the strength combinations
</pre>
<p>K is the effective length factor (1.0 for pin-ended struts); the &phi;
factors are 0.90 by default. Above D/t = 0.45 E/Fy the specification does not
apply, and the program sets the capacity to zero. The design wall thickness
factor (Sec. B4.2) can be set to 0.93 for ERW HSS such as ASTM A500. It
changes the strength checks only, not the stiffness. Auto-sizing uses the
same design wall. Tension rupture at the connection is not checked, because
the hub detail is not modelled. AISC 360 is a steel code: for aluminium the
program still gives forces and displacements, but every run says that the
member checks do not apply.</p>

<h2 id='buckling'><a name='buckling'></a>5 Linear buckling</h2>
<pre>
  [K + &lambda; Kg(N_ref)] &phi; = 0

         N  |  (I - c c^T)   -(I - c c^T) |
  kg  =  -  |                             |     geometric stiffness of a bar
         L  | -(I - c c^T)    (I - c c^T) |
</pre>
<p>N_ref are the member forces of the governing combination. The problem is
solved as (&minus;Kg) &phi; = (1/&lambda;) K &phi;: the largest eigenvalue
&mu; gives the smallest critical load factor &lambda;_cr = 1/&mu;. Springs
add to K. &lambda;_cr is the bifurcation load factor of the perfect dome.
Real domes are sensitive to imperfections, so the true limit load can be much
lower (see Section 6). Reference: McGuire, Gallagher and Ziemian (2000),
Ch. 9.</p>

<h2 id='nonlinear'><a name='nonlinear'></a>6 Geometrically nonlinear analysis</h2>
<p>The bars are corotational with engineering strain (Crisfield 1991, Ch. 3).
Displacements can be large; strains stay small and the material stays
elastic.</p>
<pre>
  d = d0 + uj - ui,    l = |d|,    c = d / l
  eps = (l - L) / L,   N = EA eps,  f = N [ -c ; +c ]
  k_T = [ k3  -k3 ; -k3  k3 ],   k3 = (EA/L) c c^T + (N/l) (I - c c^T)
</pre>
<p>The load is &lambda; F_ref, where F_ref is the D, S, L and W part of one
combination (T is left out). Equilibrium needs
r = &lambda; F_ref &minus; f(u) = 0. Each iteration forms the tangent again
(full Newton-Raphson) and solves K_T dv_r = r and K_T dv_t = F_ref
(McGuire et al. 2000, Ch. 12).</p>
<ul>
<li>Load control: &lambda; is fixed in each step. It cannot pass a load
maximum.</li>
<li>Displacement control: one DOF (the crown uz by default) is fixed in each
step, and d&lambda; = &minus;v_r,c / v_t,c. It passes a load maximum but
fails at a snap-back.</li>
<li>Arc-length (Crisfield 1981; Crisfield 1991, Ch. 9): each step obeys the
cylindrical constraint Du^T Du = Dl^2. The corrector solves the quadratic
below and keeps the root closest in direction to the step so far. It follows
the path through load maxima and snap-backs. The step length adapts to the
number of iterations.</li>
</ul>
<pre>
  a d&lambda;^2 + b d&lambda; + c = 0
  a = v_t . v_t,   b = 2 v_t . (Du + v_r),   c = |Du + v_r|^2 - Dl^2
</pre>
<p>Stability is tracked with the sign of det K_T, the number of negative
eigenvalues (up to 600 free DOFs) and the stiffness parameter
S = (F^T K0^-1 F) / (F^T K_T^-1 F). The first load maximum &lambda;_limit is
refined, then called a limit point or a bifurcation. The program also draws
the linear path u = &lambda; K^-1 F and the P-Delta path
u = &lambda; (K + &lambda; Kg)^-1 F, and reports the ratio
&lambda;_limit / &lambda;_cr. It reports the &lambda; at which the first strut
reaches its Euler load and its design strength; the path itself does not
include strut buckling or yielding.</p>
<p><b>Imperfections.</b> The nodes are moved by a &times; &phi;, where &phi; is a
buckling mode scaled so that its largest component is 1. Both signs are run,
and the lower limit load is kept. Bruno et al. (2016) report amplitudes from
about span/10000 to span/300 in use. The Imperfection study dialog starts
with span/1000, span/500, span/300 and span/200.</p>

<h2 id='modal'><a name='modal'></a>7 Modal analysis</h2>
<pre>
  K &phi;_n = &omega;_n^2 M &phi;_n,    f_n = &omega;_n / (2 &pi;),    T_n = 1 / f_n
  &phi;_n^T M &phi;_n = 1               mass-normalised modes
  &Gamma;_nd = &phi;_n^T M &iota;_d            participation factor, direction d
  M*_nd = &Gamma;_nd^2                effective modal mass
  sum over all modes of M*_nd = free mass in direction d
</pre>
<p>&iota;_d is 1 at every free DOF in direction d and 0 elsewhere. The
default mass is lumped: half of each strut's mass at each end, plus the
cladding dead load / g, hung dead nodal loads / g, an optional share of the
snow and any added masses, the same in x, y and z. A consistent bar mass,
(&rho; A L / 6) [2I I ; I 2I], is an option. The running sum of the effective
mass shows how many modes capture a given share of the mass; the program
reports the count for 90 %. With a pre-stress the problem is
(K + Kg) &phi; = &omega;^2 M &phi;: compression lowers the frequencies, and the
first one reaches zero at &lambda;_cr. Checks: &Phi;^T M &Phi; = I,
&Phi;^T K &Phi; = diag(&omega;^2) and the eigen residual. Reference: Chopra
(2023), Ch. 10 and 13.</p>

<h2 id='fabrication'><a name='fabrication'></a>8 Fabrication geometry</h2>
<pre>
  chord factor      cf = L / R
  axial angle       &alpha; = arccos(L / 2R)        strut to hub radius
  tangent angle     &tau; = 90 deg - &alpha; = arcsin(L / 2R)
  cut length        L_cut = L - 2 o           o = hub offset
  dihedral angle    &delta; = 180 deg - &theta;           &theta; = angle between the outward
                                              normals of two panels (convex fold)
</pre>
<p>The axial angle is the angle between a strut and the radius to its hub
(Clinton 1971). The program computes it from the mesh vectors at both ends
and checks it against the closed form. The breakdown is Kenner's (1976)
Class I, and the chord factors are checked against the Domerama
chord-factor chart.
For a 2V dome the closed forms give cf = 2 sin(18 deg) = 0.618034 and
cf = 2 sin(arccos(1/sqrt5) / 4) = 0.546533, with axial angles 72.000 and
74.141 deg. Hub types group the nodes by location (interior, edge, base) and
by the order of the strut families around the hub. Panel types group the
triangles by the families of their sides. A re-entrant fold has
&delta; = 180 deg + &theta;.</p>

<h2 id='sizing'><a name='sizing'></a>9 Auto-sizing and parameter studies</h2>
<p>Auto-sizing is a fully stressed design on the strut families. It analyses
the dome, gives each family the lightest catalogue section that carries its
envelope force at the target utilisation with KL/r &le; 200, and repeats. If
the deflection limit fails, every family steps up one size and never goes
below that size again. The result is a good
design, not a proven optimum. A parameter study rebuilds the dome for each
frequency, radius and dome fraction, and runs the full analysis on each
one.</p>

<h2 id='vv'><a name='vv'></a>10 Verification and validation</h2>
<ul>
<li><b>Built-in checks</b> (menu "Run verification suite", or
python run_geodome.py --verify): 75 checks. Some compare with closed-form
solutions and published chord factors: single bars and tripods, Euler's
formula, sections cut through the dome, membrane theory (n_&phi; = p R / 2).
The others are identities (symmetry, equilibrium, strain energy): they test
consistency, not accuracy.</li>
<li><b>Hand calculations</b>: validation/run_validation.py compares two
hand-solved dome models (80 quantities). validation/run_validation_ext.py adds
temperature, settlement, springs, the force method, the von Mises path,
modal values, fabrication angles, the wind field (F9) and the AISC E3, E7 and
B4.2 strength (F10).</li>
<li><b>Benchmarks</b>: the von Mises two-bar truss in closed form (derived
here; the problem of Pecknold et al. 1985) and the 24-member star dome
(Greco and Venturini 2006; the precise values are computed in this work and
agree with their figures).</li>
<li><b>OpenSees</b> (McKenna et al. 2010; Zhu et al. 2018): the same model
runs with the corotTruss element, which uses the same strain measure and
tangent. Static, modal and nonlinear results are compared. OpenSees does not
check the load generation, the AISC checks or &lambda;_cr.</li>
</ul>
<p>Every static solve also reports its own equilibrium residual and energy
error, and every nonlinear point its residual.</p>

<h2 id='refs'><a name='refs'></a>References</h2>
<p class='ref'>AISC (2016). <i>Specification for Structural Steel Buildings</i>
(ANSI/AISC 360-16). American Institute of Steel Construction, Chicago, IL.</p>
<p class='ref'>ASCE (2017). <i>Minimum Design Loads and Associated Criteria
for Buildings and Other Structures</i> (ASCE/SEI 7-16). American Society of
Civil Engineers, Reston, VA. doi:10.1061/9780784414248.</p>
<p class='ref'>Bruno, L., Sassone, M., Venuti, F. (2016). Effects of the
equivalent geometric nodal imperfections on the stability of single layer grid
shells. <i>Engineering Structures</i> 112, 184-199.
doi:10.1016/j.engstruct.2016.01.017.</p>
<p class='ref'>Chopra, A. K. (2023). <i>Dynamics of Structures: Theory and
Applications to Earthquake Engineering</i>, 6th ed. Pearson, Hoboken, NJ.
Ch. 10 and 13.</p>
<p class='ref'>Clinton, J. D. (1971). <i>Advanced Structural Geometry
Studies, Part I: Polyhedral Subdivision Concepts for Structural
Applications</i>. NASA CR-1734.</p>
<p class='ref'>Crisfield, M. A. (1981). A fast incremental/iterative solution
procedure that handles "snap-through". <i>Computers &amp; Structures</i>
13(1-3), 55-62. doi:10.1016/0045-7949(81)90108-5.</p>
<p class='ref'>Crisfield, M. A. (1991). <i>Non-linear Finite Element Analysis
of Solids and Structures, Volume 1: Essentials</i>. Wiley, Chichester.
Ch. 3 and 9.</p>
<p class='ref'>Greco, M., Venturini, W. S. (2006). Stability analysis of
three-dimensional trusses. <i>Latin American Journal of Solids and
Structures</i> 3(3), 325-344.</p>
<p class='ref'>Kenner, H. (1976). <i>Geodesic Math and How to Use It</i>.
University of California Press, Berkeley.</p>
<p class='ref'>McGuire, W., Gallagher, R. H., Ziemian, R. D. (2000).
<i>Matrix Structural Analysis</i>, 2nd ed. Wiley, New York. Ch. 3, 9 and
12.</p>
<p class='ref'>McKenna, F., Scott, M. H., Fenves, G. L. (2010). Nonlinear
finite-element analysis software architecture using object composition.
<i>Journal of Computing in Civil Engineering</i> 24(1), 95-107.
doi:10.1061/(ASCE)CP.1943-5487.0000002.</p>
<p class='ref'>Pecknold, D. A., Ghaboussi, J., Healey, T. J. (1985).
Snap-through and bifurcation in a simple structure. <i>Journal of Engineering
Mechanics</i> 111(7), 909-922. doi:10.1061/(ASCE)0733-9399(1985)111:7(909).</p>
<p class='ref'>Zhu, M., McKenna, F., Scott, M. H. (2018). OpenSeesPy: Python
library for the OpenSees finite element framework. <i>SoftwareX</i> 7, 6-11.
doi:10.1016/j.softx.2017.10.009.</p>

</body></html>
"""


def theory_html():
    """The Theory help page as one HTML string."""
    return THEORY_HTML
