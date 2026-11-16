# GeoDomeLab extended hand calculations (validation cases F1 to F10)

This document is produced by validation/hand_calc_ext.py. The script does not import the GeoDomeLab package. Each case defines a structure, solves it with a hand method (equilibrium, compatibility, the force method or the stiffness method as taught in the course), prints every intermediate number and ends with an independent check. The reference values feed validation/run_validation_ext.py, which runs the same structures through GeoDomeLab and compares.

Conventions. z is vertical (up), tension is positive. A member from node i to node j has the unit vector c = (x_j - x_i)/L, the elongation e = c . (u_j - u_i) and the force N = EA (e/L - eps0) with the free strain eps0 = alpha dT. The member pushes or pulls its i end with +N c and its j end with -N c. A reaction is minus the sum of member forces and loads at its node. Values are kept in SI units and printed in kN, mm and MPa.

| parameter | value | source |
|---|---|---|
| E (A500 Gr. C) | 200 GPa | AISC 360 / ASTM A500 |
| Fy | 345 MPa | ASTM A500 Gr. C round |
| alpha | 1.17e-05 1/degC | AISC Manual 6.5e-6 /degF for carbon steel |
| rho | 7850 kg/m3 | steel |
| dT | +30 degC | uniform rise above erection |
| CHS 88.9 x 4.0 | A = 1066.88 mm2, I = 963398 mm4 | pi/4 (D^2 - d^2) |
| CHS 76.1 x 4.0 | A = 906.035 mm2 | pi/4 (D^2 - d^2) |
| Model A | R = 5 m, L = 5.25731 m, k = EA/L = 40586.7 kN/m | hand_calc.py |

## F1  Heated fixed-fixed bar

Purpose. The most basic self-straining case. A bar that is held at both ends cannot expand, so heating it produces a compressive force that does not depend on its length. This checks the sign and size of the thermal strain eps0 = alpha dT and the thermal reactions in GeoDomeLab.

### F1a  Inclined bar pinned at both ends

Data. CHS 88.9 x 4.0: A = 1066.88 mm2, I = 963398 mm4. E = 200 GPa, alpha = 1.17e-05 1/degC, dT = +30 degC. Node 1 = (0, 0, 0) m, node 2 = (2, 1, 2) m, both pinned (all three translations fixed).

```text
L = |x2 - x1| = sqrt(2^2 + 1^2 + 2^2) = 3 m
c = (x2 - x1)/L = (0.666667, 0.333333, 0.666667)
EA = 213376.97 kN

Compatibility: the ends cannot move, so the total elongation is 0
  e = N L/(EA) + alpha dT L = 0
  N = -EA alpha dT
  free thermal elongation alpha dT L = 1.053 mm
  N = -213376.97 x 1.17e-05 x 30 = -74.8953 kN  (compression)
  sigma = N/A = -E alpha dT = -70.2 MPa  (0.2035 Fy)

Reactions (force of the bar on node 1 is +N c, on node 2 is -N c):
  R1 = -N c = (49.9302, 24.9651, 49.9302) kN
  R2 = +N c = (-49.9302, -24.9651, -49.9302) kN
```

Independent check by the stiffness method. With every degree of freedom fixed, u = 0 and the reactions equal minus the thermal equivalent nodal loads F0 = EA alpha dT [-c ; +c].

```text
F0 at node 1 = -EA alpha dT c = (-49.9302, -24.9651, -49.9302) kN  ->  R1 = -F0 = (49.9302, 24.9651, 49.9302) kN
F0 at node 2 = +EA alpha dT c = (49.9302, 24.9651, 49.9302) kN  ->  R2 = -F0 = (-49.9302, -24.9651, -49.9302) kN
```

- Total elongation e = -2.17e-19 m (must be 0).
- Sum of reactions = (0, 0, 0) N, moment about the origin = (0, 0, 0) N m.
- Reactions by equilibrium and by -F0 differ by 0 N.
- Engineering note: the Euler load of this 3 m pin-ended bar is Pe = pi^2 EI/L^2 = 211.297 kN, so the bar would buckle at dT = Pe/(EA alpha) = 84.64 degC. The compressive thermal force does not depend on L, but the buckling temperature falls with L^2.

### F1b  Stepped bar between two fixed ends (one redundant)

Data. Collinear bar along x: node 1 (0,0,0) fixed, node 2 (2,0,0) free in x (y, z restrained), node 3 (3.5,0,0) fixed. Segment 1-2 is CHS 88.9 x 4.0 (EA1 = 213377 kN), segment 2-3 is CHS 76.1 x 4.0 (A = 906.035 mm2, EA2 = 181207 kN). Both segments are heated by dT = +30 degC. Force method with the redundant X = axial force N (equal in both segments by equilibrium of node 2): release node 3 in x, let the bar expand freely, then push it back.

```text
free expansion at the release  d_T = alpha dT (L1 + L2) = 1.2285 mm
flexibility  f = L1/EA1 + L2/EA2 = 1.76509e-08 m/N
compatibility  d_T + f N = 0
  N = -d_T / f = -69.5998 kN
node 2 moves by the elongation of segment 1
  u2 = (N/EA1 + alpha dT) L1 = 0.0496351 mm
reactions  R1x = -N = 69.5998 kN,  R3x = +N = -69.5998 kN
```

Independent check by the stiffness method (one free DOF): (EA1/L1 + EA2/L2) u2 = F0 with F0 = alpha dT (EA1 - EA2), the j end load of segment 1 plus the i end load of segment 2.

```text
u2 = alpha dT (EA1 - EA2)/(EA1/L1 + EA2/L2) = 0.0496351 mm
segment 2 elongation (N/EA2 + alpha dT) L2 = -0.0496351 mm = -u2 (node 3 fixed)
```

- u2 by the two methods differs by -8.13e-20 m.
- Compatibility of segment 2: e2 + u2 = -1.08e-19 m.

## F2  Model A (1V cap) under a uniform temperature change

Structure. Validation Model A of hand_calc.py: the 1V pentagonal cap of a sphere of radius R = 5 m. Node 1 is the crown (free), nodes 2 to 6 are the pinned base ring at 0, 72, 144, 216 and 288 deg. Members 1 to 5 are the spokes (ring node to crown), members 6 to 10 the ring (2-3, 3-4, 4-5, 5-6, 6-2). All ten struts have the 1V length L = 5.25731 m and the section CHS 88.9 x 4.0, EA = 213377 kN, k = EA/L = 40586.7 kN/m. Crown height H = R - R/sqrt5 = 2.76393 m, ring plan radius a = 2R/sqrt5 = 4.47214 m.

```text
spoke direction (ring node -> crown)  c_i = (-s cos th_i, -s sin th_i, t)
  s = a/L = 0.850651,  t = H/L = 0.525731
free block  Kff = k diag(5/2 s^2, 5/2 s^2, 5 t^2) = diag(73422.1, 73422.1, 56089.5) kN/m
```

### F2a  Every strut heated by dT = +30 degC

Hand method (compatibility first). The five ring struts join two pinned nodes, so each one is the fixed-fixed bar of F1 and carries N = -EA alpha dT. The crown is free. If it rises by w, every spoke lengthens by c_i . (0,0,w) = t w, the same for all five. Choosing t w = alpha dT L makes each spoke exactly as long as it wants to be, so the spokes carry no force. Then the crown is in equilibrium with zero spoke forces, every member is compatible and every force obeys the material law. A linear elastic solution that satisfies all three is the unique solution.

```text
alpha dT L = 1.84532 mm
crown rise  w = alpha dT L / t = alpha dT L^2 / H
  w = 1.84532 / 0.525731 = 3.51 mm
spokes  N = k (t w) - EA alpha dT = 0
ring    N = -EA alpha dT = -74.8953 kN

Reaction at a ring node (two ring bars in compression push the
node outwards, the unit vectors to the two neighbours add to
-2 sin36 e_r):
  R = -2 sin36 EA alpha dT e_r   (inward)
  R_radial = -88.0447 kN,  R_z = 0
```

| member | type | elongation [mm] | N [kN] |
|---|---|---|---|
| 1 | spoke | 1.84532 | 1.45519e-14 |
| 2 | spoke | 1.84532 | 1.45519e-14 |
| 3 | spoke | 1.84532 | 2.91038e-14 |
| 4 | spoke | 1.84532 | 1.45519e-14 |
| 5 | spoke | 1.84532 | 2.91038e-14 |
| 6 | ring | 0 | -74.8953 |
| 7 | ring | 0 | -74.8953 |
| 8 | ring | 0 | -74.8953 |
| 9 | ring | 0 | -74.8953 |
| 10 | ring | 0 | -74.8953 |

Check by the stiffness method. The thermal equivalent load at the crown (the j end of every spoke) is F0 = sum EA alpha dT c_i = EA alpha dT (0, 0, 5t). The ring loads act on supported DOFs only.

```text
F0 = (2.18279e-14, -7.27596e-15, 196.874) kN
u = Kff^-1 F0 = (2.97293e-16, -9.90977e-17, 3.51) mm
```

- Crown rise by compatibility and by Kff^-1 F0 differ by 4.34e-19 m.
- Crown equilibrium residual 5.41e-11 N, sum of reactions 6.27e-11 N, moment 8.98e-11 N m (a self-straining load has self-equilibrated reactions).
- Radial reaction from the joint sums equals the formula to 5.82e-11 N at all five supports (five-fold symmetry).
- Engineering note: the ring force -74.9 kN does not depend on the strut length. The AISC E3 design strength of this 5.257 m strut is phi Pn = 54.31 kN (hand_calc.py), so a 30 degC rise alone gives a utilisation of 1.38. A ring between pinned supports needs either a temperature check or a sliding bearing. The spokes are stress free because the crown can rise 3.51 mm, which is L^2/H^2 = 3.618 times the free rise alpha dT H of an unrestrained cap.

### F2b  Only spoke 1 heated by dT = +30 degC (sun on one side)

A non-uniform temperature field. Only member 1 (node 2 to the crown) gets eps0 = alpha dT. The crown now moves sideways and the five spokes share a self-equilibrated set of forces. Stiffness method with the thermal equivalent load of member 1 at the crown.

```text
F0 = EA alpha dT c_1 = 74.8953 x (-0.850651, 0, 0.525731) = (-63.7098, 0, 39.3748) kN
Kff u = F0   (Gaussian elimination, multipliers 0, 0, 0)
u = (-0.86772, 0, 0.702) mm
closed form  u = alpha dT L ( -1/(2.5 s), 0, 1/(5 t) )
            = (-0.86772, 0, 0.702) mm
N_1 = EA alpha dT (0.4 + 0.2 - 1) = -0.4 EA alpha dT
N_i = EA alpha dT (0.4 cos th_i + 0.2)   for i = 2..5
```

| member | N [kN] | N / (EA alpha dT) | closed form |
|---|---|---|---|
| 1 | -29.9581 | -0.4 | -0.4 |
| 2 | 24.2366 | 0.323607 | 0.323607 |
| 3 | -9.25757 | -0.123607 | -0.123607 |
| 4 | -9.25757 | -0.123607 | -0.123607 |
| 5 | 24.2366 | 0.323607 | 0.323607 |
| 6 | 0 | 0 | 0 |
| 7 | 0 | 0 | 0 |
| 8 | 0 | 0 | 0 |
| 9 | 0 | 0 | 0 |
| 10 | 0 | 0 | 0 |

- Spoke forces match the closed form to 2.22e-16 (relative to EA alpha dT).
- Crown equilibrium residual 1.51e-11 N, sum of reactions 1.51e-11 N.
- The five spoke forces add to zero, as vertical equilibrium of the unloaded crown requires (all spokes have the same rise t).

## F3  Model A with one support movement

Prescribed support displacements are the second self-straining load. The stiffness method with prescribed values us on the restrained DOFs solves Kff uf = Ff - Kfs us. The crown is linked to a ring node only through that node's spoke, so the coupling block is Kfs = -k c c^T for that spoke and the right-hand side is k c (c . us). Ring struts do not touch the crown.

### F3a  Hand node 2 settles vertically by 5 mm

```text
us (node 2) = (0, 0, -5) mm,   spoke 1 c = (-0.850651, 0, 0.525731)
Kfs us = -k c (c . us),  c . us = -2.62866 mm
right-hand side  -Kfs us = (90.7546, 0, -56.0895) kN

3 x 3 system (kN/m, mm -> kN):
  [  73422.05      0      0 ] {D1}   {  90.75465}
  [         0 73422.05      0 ] {D2} = {         0}
  [         0      0 56089.46 ] {D3}   { -56.08946}
  elimination multipliers (row, pivot, m): (2,1,0), (3,1,0), (3,2,0)
  back substitution  D3 = b3/U33, D2 = b2/U22, D1 = b1/U11
  D = (1.23607, 0, -1) mm
```

| member | type | elongation [mm] | N [kN] |
|---|---|---|---|
| 1 | spoke | 1.05146 | 42.6754 |
| 2 | spoke | -0.850651 | -34.5251 |
| 3 | spoke | 0.32492 | 13.1874 |
| 4 | spoke | 0.32492 | 13.1874 |
| 5 | spoke | -0.850651 | -34.5251 |
| 6 | ring | 0 | 0 |
| 7 | ring | 0 | 0 |
| 8 | ring | 0 | 0 |
| 9 | ring | 0 | 0 |
| 10 | ring | 0 | 0 |

| support | Rx [kN] | Ry [kN] | Rz [kN] |
|---|---|---|---|
| 2 | 36.3019 | 0 | -22.4358 |
| 3 | -9.07546 | -27.9314 | 18.1509 |
| 4 | -9.07546 | 6.59371 | -6.93304 |
| 5 | -9.07546 | -6.59371 | -6.93304 |
| 6 | -9.07546 | 27.9314 | 18.1509 |

Closed form. With c_1 = (-s, 0, t) and us = (0, 0, -delta) the right-hand side is k t delta (s, 0, -t), so

```text
D1 = t delta/(2.5 s),  D2 = 0,  D3 = -delta/5
   = (1.23607, 0, -1) mm
N_1 = 0.4 k t delta,   N_i = -k t delta (0.4 cos th_i + 0.2)
k t delta = 106.688 kN
```

Independent checks. Joint equilibrium of the crown, global equilibrium of the self-equilibrated reactions, and Clapeyron's theorem for a prescribed displacement: the strain energy stored in the struts equals the work of the support reaction on the settlement, U = sum N^2 L/(2EA) = 1/2 R2 . us.

- Crown displacement and spoke forces equal the closed form to 2.17e-19 m and 5.55e-17 (relative).
- Crown equilibrium 1.42e-11 N, sum of reactions 1.42e-11 N, moment 4.12e-11 N m.
- U = 56.0895 J, 1/2 R2 . us = 56.0895 J, difference -7.11e-15 J.
- The ring struts stay stress free because they are horizontal and a vertical settlement does not change their length (first-order theory).
- Engineering note: a 5 mm settlement puts -34.53 kN of compression into spokes 2 and 5. Settlement forces scale with EA, so stiffer struts attract more of them.

### F3b  Hand node 3 spreads radially outwards by 5 mm

```text
us (node 3) = (1.54508, 4.75528, 0) mm,   spoke 2 c = (-0.262866, -0.809017, 0.525731)
Kfs us = -k c (c . us),  c . us = -4.25325 mm
right-hand side  -Kfs us = (45.3773, 139.657, -90.7546) kN

3 x 3 system (kN/m, mm -> kN):
  [  73422.05      0      0 ] {D1}   {  45.37732}
  [         0 73422.05      0 ] {D2} = {   139.657}
  [         0      0 56089.46 ] {D3}   { -90.75465}
  elimination multipliers (row, pivot, m): (2,1,0), (3,1,0), (3,2,0)
  back substitution  D3 = b3/U33, D2 = b2/U22, D1 = b1/U11
  D = (0.618034, 1.90211, -1.61803) mm
```

| member | type | elongation [mm] | N [kN] |
|---|---|---|---|
| 1 | spoke | -1.37638 | -55.8628 |
| 2 | spoke | 1.7013 | 69.0502 |
| 3 | spoke | -1.37638 | -55.8628 |
| 4 | spoke | 0.525731 | 21.3377 |
| 5 | spoke | 0.525731 | 21.3377 |
| 6 | ring | 2.93893 | 119.281 |
| 7 | ring | 2.93893 | 119.281 |
| 8 | ring | 0 | 0 |
| 9 | ring | 0 | 0 |
| 10 | ring | 0 | 0 |

| support | Rx [kN] | Ry [kN] | Rz [kN] |
|---|---|---|---|
| 2 | 22.5921 | -96.5006 | 29.3688 |
| 3 | 61.4824 | 189.223 | -36.3019 |
| 4 | -74.999 | -64.7914 | 29.3688 |
| 5 | -14.6844 | -10.6688 | -11.2179 |
| 6 | 5.60895 | -17.2626 | -11.2179 |

Closed form. c_2 . us = -s Delta, so the crown moves by 0.4 Delta along e_r(72 deg) and drops by s Delta/(5 t). The two ring struts at node 3 are stretched by Delta sin36.

```text
D = 0.4 Delta e_r(72) + (0, 0, -s Delta/(5 t))
  = (0.618034, 1.90211, -1.61803) mm
N_2 = 0.4 k s Delta,  N_i = -k s Delta (0.4 cos(th_i - 72) + 0.2)
ring 6 and 7: N = k Delta sin36 = 119.281 kN
```

- Crown displacement, spoke and ring forces equal the closed form to 3.25e-19 m, 1.67e-16 (relative) and 1.46e-11 N.
- Crown equilibrium 6.2e-11 N, sum of reactions 6.08e-11 N, moment 1.65e-10 N m.
- U = 497.403 J equals 1/2 R3 . us = 497.403 J (difference 5.68e-14 J).
- F3a and F2b have the same pattern 0.4 / (0.4 cos + 0.2): heating one spoke and moving its support are both a change of that spoke's free length.

## F4  Planar three-bar truss by the force (flexibility) method

Purpose. A statically indeterminate truss solved by the force method, which shares nothing with the stiffness method used by GeoDomeLab: the unknowns are forces, compatibility is written with the principle of virtual forces, and displacements come last from the unit-load method. Three load types are treated: a point load, a temperature change in one bar and a support settlement.

Structure. Bars lie in the x-z plane and meet at the free node D = (0, 0, 3) m. Supports (pinned): A = (-4, 0, 0), B = (0, 0, 0), C = (3, 0, 0) m. Node D is restrained out of plane (uy = 0) so it has two free DOFs (ux, uz). E = 200 GPa, alpha = 1.17e-05 1/degC.

| bar | from-to | section | A [mm2] | L [m] | cx | cz | L/EA [m/MN] |
|---|---|---|---|---|---|---|---|
| 1 | A-D | CHS 88.9 x 4.0 | 1066.88 | 5 | 0.8 | 0.6 | 0.0234327 |
| 2 | B-D | CHS 76.1 x 4.0 | 906.035 | 3 | 0 | 1 | 0.0165556 |
| 3 | C-D | CHS 88.9 x 4.0 | 1066.88 | 4.24264 | -0.707107 | 0.707107 | 0.0198833 |

Loads. Case P: Px = 15 kN, Pz = -60 kN at D. Case T: bar 2 heated by dT = +30 degC. Case S: support C settles by (0, 0, -5) mm.

### Step 1  Degree of indeterminacy and the primary structure

```text
plane truss: i = m + r - 2j = 3 + 6 - 2*4 = 1
redundant X1 = N2 (force in bar 2), released by cutting bar 2
primary structure = bars 1 and 3 (a determinate two-bar truss)
joint D: forces on D are -N_i c_i, so sum N_i c_i = P
```

### Step 2  Primary structure under the load (n0) and under X1 = 1 (n1)

```text
[c1 c3] n = rhs,  det = c1x c3z - c3x c1z = 0.989949
n_1 = (rhs_x c3z - c3x rhs_z)/det,  n_3 = (c1x rhs_z - rhs_x c1z)/det
load:      rhs = P   = (15, -60) kN  ->  n0 = (-32.1429, 0, -57.5787) kN
X1 = 1:    rhs = -c2 = (0, -1)  ->  n1 = (-0.714286, 1, -0.808122)
unit x:    rhs = (1,0)  ->  n_x = (0.714286, 0, -0.606092)
unit z:    rhs = (0,1)  ->  n_z = (0.714286, 0, 0.808122)
```

### Step 3  Flexibility coefficient and compatibility terms

```text
virtual forces (X1 = 1 system) x real deformations:
  sum n1_i e_i - sum Rbar_k . delta_k = 0   (gap at the cut = 0)
  e_i = (n0_i + X1 n1_i) L_i/EA_i + alpha dT_i L_i
  f11 X1 + D1P + D1T + D1S = 0
  f11 = sum n1_i^2 L_i/EA_i
  D1P = sum n1_i n0_i L_i/EA_i
  D1T = sum n1_i alpha dT_i L_i
  D1S = -sum Rbar_k . delta_k,  Rbar_k = -n1_i c_i at the support of bar i
```

| bar | n0 [kN] | n1 [-] | n1^2 L/EA [m/MN] | n1 n0 L/EA [mm] | n1 alpha dT L [mm] |
|---|---|---|---|---|---|
| 1 | -32.1429 | -0.714286 | 0.0119555 | 0.537996 | 0 |
| 2 | 0 | 1 | 0.0165556 | 0 | 1.053 |
| 3 | -57.5787 | -0.808122 | 0.012985 | 0.925183 | 0 |
| sum |  |  | 0.0414961 | 1.46318 | 1.053 |

```text
Rbar_C = -n1_3 c3 = (-0.571429, 0.571429),  delta_C = (0, -5) mm
D1S = 2.85714 mm
f11 = 4.14961e-08 m/N
```

### Step 4  Redundant, bar forces, displacements and reactions

```text
X1 = -(D1P + D1T + D1S)/f11,   N_i = n0_i + X1 n1_i
u_D . e = sum_{1,3} n_e,i e_i - sum Rtilde_e,k . delta_k   (unit load e on the primary)
R_k = -N_i c_i at the support of bar i
```

| case | X1 = N2 [kN] | N1 [kN] | N2 [kN] | N3 [kN] | ux_D [mm] | uz_D [mm] |
|---|---|---|---|---|---|---|
| load | -35.2606 | -6.95671 | -35.2606 | -29.0838 | 0.234053 | -0.583762 |
| temperature | -25.3759 | 18.1256 | -25.3759 | 20.5068 | 0.0562506 | 0.632886 |
| settlement | -68.8532 | 49.1809 | -68.8532 | 55.6418 | 2.29548 | -1.13991 |
| combined | -129.49 | 60.3498 | -129.49 | 47.0648 | 2.58579 | -1.09079 |

| case | RA_x | RA_z | RB_x | RB_z | RC_x | RC_z |
|---|---|---|---|---|---|---|
| load | 5.56537 | 4.17403 | 0 | 35.2606 | -20.5654 | 20.5654 |
| temperature | -14.5005 | -10.8754 | 0 | 25.3759 | 14.5005 | -14.5005 |
| settlement | -39.3447 | -29.5085 | 0 | 68.8532 | 39.3447 | -39.3447 |
| combined | -48.2798 | -36.2099 | 0 | 129.49 | 33.2798 | -33.2798 |

Reactions in kN. The out-of-plane reaction at D and all y reactions are zero.

### Step 5  Independent checks

Compatibility of every bar (the elongation computed from the unit-load displacements equals N L/EA + alpha dT L, including the cut bar 2), equilibrium of joint D, global force and moment equilibrium, and a second solution by the displacement method (2 x 2 stiffness matrix of joint D solved by Cramer's rule).

| case | max compat. [m] | joint D [N] | sum R + P [N] | moment [N m] | \|dN\| vs stiffness [N] | \|du\| [m] |
|---|---|---|---|---|---|---|
| load | 2.2e-19 | 1.8e-12 | 1.8e-12 | 0 | 7.3e-12 | 3.4e-19 |
| temperature | 0 | 3.6e-12 | 3.6e-12 | 7.3e-12 | 1.1e-11 | 6.8e-20 |
| settlement | 2.2e-19 | 7.3e-12 | 7.3e-12 | 1.5e-11 | 1.5e-11 | 4.8e-19 |
| combined | 3.3e-19 | 1.5e-11 | 1.5e-11 | 1.5e-11 | 1.5e-11 | 4.3e-19 |

- Superposition: X1(P) + X1(T) + X1(S) = -129.49 kN, combined case -129.49 kN.
- Engineering note: heating bar 2 alone by 30 degC gives N2 = -25.38 kN, and a 5 mm settlement of support C gives N2 = -68.85 kN, compared with -35.26 kN from the 60 kN load. Self-straining loads are of the same order as gravity loads in a stiff indeterminate truss.

## F5  Bar and elastic spring support

Purpose. Elastic (spring) supports. A CHS 88.9 x 4.0 bar of length L = 3 m lies on the x axis from node 1 (0,0,0) to node 2 (3,0,0). Its axial stiffness is EA/L = 71125.7 kN/m. The spring has k = 50000 kN/m (GeoDomeLab's default support spring) and acts in x. The y and z DOFs of both nodes are restrained, so each case has one or two free DOFs. The spring force on its node is -k u, reported as the reaction.

### F5a  In series: spring at node 1, load at node 2

```text
node 2:  bar force N = P (statically determinate)
node 1:  the spring carries the same force, k u1 = N
flexibilities add:  u2 = P/k + P L/EA,  k_eq = 1/(1/k + L/EA)
  u1 = 1 mm,  u2 = 1.70298 mm,  N = 50 kN
  spring reaction at node 1 = -k u1 = -50 kN,  k_eq = 29360.3 kN/m
```

- Energy: 1/2 P u2 = 42.5745 J, 1/2 N^2 L/EA + 1/2 k u1^2 = 42.5745 J.
- Equilibrium: P + spring reaction = 0 N.

### F5b  In parallel: node 1 fixed, spring and load at node 2

```text
node 2:  (EA/L + k) u2 = P
  u2 = 0.412794 mm
  bar N = (EA/L) u2 = 29.3603 kN  (share EA/L/(EA/L + k) = 0.587206)
  spring reaction = -k u2 = -20.6397 kN,  fixed reaction at node 1 = -N = -29.3603 kN
```

- Equilibrium: P + R1 + spring = -3.64e-12 N.
- Energy: 1/2 P u2 = 1/2 (EA/L) u2^2 + 1/2 k u2^2 holds to 0 J.

### F5c  Heated bar on a spring: node 1 fixed, spring at node 2, dT = +30 degC

```text
thermal equivalent load at node 2:  F0 = EA alpha dT
(EA/L + k) u2 = EA alpha dT
  u2 = alpha dT L / (1 + k L/EA)
  N  = EA (u2/L - alpha dT) = -EA alpha dT k/(k + EA/L)
  u2 = 0.618327 mm   (free expansion 1.053 mm)
  N  = -30.9164 kN    (fully fixed -74.8953 kN)
  spring reaction -k u2 = -30.9164 kN = N
```

- Node 2 equilibrium: the compressed bar pushes node 2 by -N = 30.9164 kN and the spring pulls back by -30.9164 kN.
- Limits: k -> 0 gives the free expansion alpha dT L with N = 0, k -> infinity gives the fixed bar of F1 (N -> -EA alpha dT = -74.8953 kN).
- Engineering note: an elastic support with k comparable to EA/L relieves 58.7 % of the fully restrained thermal force.

## F6  von Mises two-bar truss: snap-through

Purpose. The classic benchmark for geometrically nonlinear (large displacement, small strain) truss analysis. Two CHS 88.9 x 4.0 bars run from pinned supports at (-a, 0, 0) and (+a, 0, 0) to the apex (0, 0, h), with a = 2.5 m and h = 0.125 m (h/a = 0.05, a shallow truss). The apex is restrained out of plane (uy = 0) and carries a downward load P. Strain is the engineering strain (l - L)/L, the measure GeoDomeLab's corotational bar uses. EA = 213377 kN.

### Equilibrium path (exact, closed form)

```text
L = sqrt(a^2 + h^2) = 2.503123049 m
apex deflection w (down), y = h - w, l = sqrt(a^2 + y^2)
eps = (l - L)/L,  N = EA eps  (both bars, by symmetry)
vertical equilibrium of the apex in the deformed position:
  P = -2 N y/l
  P(w) = 2 EA (L - l)/L * (h - w)/l
  (the same as F_z = 2 EA (l - L)/L (h - w)/l with F_z positive up)
  evaluated with L - l = (L^2 - l^2)/(L + l) = w (2h - w)/(L + l)
  (identical algebra, no loss of digits when l is close to L)
tangent  dP/dw = 2 EA (1/L - a^2/l^3)
         = 2 [ (EA/L)(y/l)^2 + (N/l)(a/l)^2 ]   (material + geometric)
at w = 0  dP/dw = 2 EA h^2/L^3 = 425.159 kN/m  (the linear stiffness)
```

### Limit point from dP/dw = 0

```text
dP/dw = 0  <=>  a^2/l^3 = 1/L  <=>  l* = (a^2 L)^(1/3)
written without cancellation, z = y^2/a^2, r = h/a:
  (1 + z)^3 = 1 + r^2   <=>   z^3 + 3 z^2 + 3 z - r^2 = 0
  z* = (1 + r^2)^(1/3) - 1 = expm1(log1p(r^2)/3)
  z* = 0.000832639851790493
  y* = a sqrt(z*) = 0.07213874876715412 m
  l* = a sqrt(1 + z*) = 2.501040583252037 m   ((a^2 L)^(1/3) = 2.501040583252037 m)
  w* = h - y* = 0.05286125123284588 m  = 0.4228900099 h
  P* = P(w*) = 10240.49959685 N = 10.2405 kN
  eps* = (l* - L)/L = -0.000831947,  N* = EA eps* = -177.518 kN,  sigma* = -166.4 MPa
Newton iteration on p(z) = z^3 + 3z^2 + 3z - r^2 from z = r^2/3:
  it  1: z = 0.00083263985227100912   p = +2.084e-06
  it  2: z = 0.00083263985179049299   p = +1.444e-12
  it  3: z = 0.00083263985179049299   p = +0.000e+00
  Newton y* - closed form y* = 0 m  (relative 0)
second derivative  d2P/dw2 = -6 EA a^2 y/l^5 < 0  (a maximum)
  d2P/dw2(w*) = -5.89853e+06 N/m^2
shallow-truss approximation  P* ~ 2 EA h^3/(3 sqrt3 a^3) = 10.2661 kN,  w* ~ h (1 - 1/sqrt3) = 0.0528312 m
```

### Snap-through target and the full path

```text
Under load control the apex jumps at P* from w* to the stable
inverted branch (w > 2h) where P(w) = P* again. Bisection on
P(w) - P* over [2h, w_hi]:
  w_snap = 0.269367614677 m  = 2.1549409 h  (49 bisections)
  bar force there N = 88.8146 kN (tension, the bars are
  longer than L in the inverted position)
second limit point  w = h + y* = 0.197139 m, P = -P*
zero-load points  w = 0, w = h (bars horizontal), w = 2h
```

| w/h | w [mm] | strain | N [kN] | P [kN] | P/P* |
|---|---|---|---|---|---|
| 0.1 | 12.5 | -0.000236936 | -50.5566 | 4.5455 | 0.443875 |
| 0.25 | 31.25 | -0.00054566 | -116.431 | 8.72621 | 0.852128 |
| 0.5 | 62.5 | -0.0009356 | -199.635 | 9.97865 | 0.97443 |
| 1 | 125 | -0.00124766 | -266.222 | 0 | 0 |
| 1.5 | 187.5 | -0.0009356 | -199.635 | -9.97865 | -0.97443 |
| 2 | 250 | 0 | 0 | 0 | 0 |
| 2.5 | 312.5 | 0.00155739 | 332.311 | 49.7071 | 4.85397 |

### Linearised buckling factor (2-DOF eigenproblem)

```text
reference load P_ref = 10 kN downward, linear member force
  N_ref = -P_ref L/(2h) = -100.125 kN
apex DOFs (u horizontal, w vertical), bars c = (+-a, h)/L:
  K  = (EA/L) sum c c^T        = diag(2EA a^2/L^3, 2EA h^2/L^3)
  Kg = (N_ref/L) sum (I - c c^T) = diag(2N_ref h^2/L^3, 2N_ref a^2/L^3)
  K  = diag(170063, 425.159) kN/m
  Kg = diag(-0.199501, -79.8005) kN/m
det(K + lambda Kg) = (Kuu + lambda Guu)(Kww + lambda Gww) = 0
  lambda_u = -Kuu/Guu = EA a^2/(-N_ref h^2)   (sway mode (1,0))
  lambda_w = -Kww/Gww = EA h^2/(-N_ref a^2)   (symmetric mode (0,1))
  lambda_u = 852443.0035,  lambda_w = 5.327768772
  lambda_cr = 5.327768772  ->  P_cr,lin = lambda_cr P_ref = 53.2777 kN
  closed form P_cr,lin = 2 EA h^3/(a^2 L) = 53.2777 kN
  det(K + lambda_cr Kg) = 0
```

- The exact limit load is P* = 10.2405 kN. The linearised eigenvalue gives P_cr,lin = 53.2777 kN, 5.2026 times larger. For a shallow truss the ratio tends to 3 sqrt3 = 5.1962. Linear buckling ignores the change of geometry before the limit point and is unconservative for shallow lattices, which is why GeoDomeLab 2.0 adds the nonlinear path.
- Sway stiffness at the limit point Kuu(w*) = 170347 kN/m > 0, so the apex reaches the symmetric limit point before any sway bifurcation.
- At the limit point the bar stress is -166.4 MPa < Fy = 345 MPa and |N*| = 177.5 kN is below the Euler load of one bar, 303.5 kN, so the elastic two-bar model is physically meaningful.

## F7  Model A crown natural frequencies (lumped mass)

Purpose. Free vibration. With pinned ring nodes only the crown moves, so Model A has three dynamic DOFs and the mass matrix is m I (the same lumped mass in x, y and z). Lumping rule: each strut puts half of its mass rho A L at each end node. Optionally the superimposed dead load (cladding pressure q times the surface tributary area, divided by g) is added as mass.

```text
rho A L = 7850 x 0.00106688 x 5.25731 = 44.0302 kg per strut
self weight at the crown  m_sw = 5 x 1/2 rho A L = 110.076 kg
crown surface tributary area = 5/3 x (sqrt3/4) L^2 = 19.947 m2
dead load mass  m_q = q A_trib / g = 300 x 19.947 / 9.80665 = 610.207 kg
stiffness  Kxx = Kyy = 5/2 k s^2,  Kzz = 5 k t^2
  Kxx = 73422.1 kN/m,  Kzz = 56089.5 kN/m
(K - omega^2 m I) phi = 0 with diagonal K:
  omega_x = omega_y = sqrt(Kxx/m),  omega_z = sqrt(Kzz/m)
  f = omega/(2 pi),  T = 1/f,  mass-normalised phi = e/sqrt(m)
  f_z/f_xy = sqrt(2) t/s, independent of the mass
```

| mass | m [kg] | f_x [Hz] | f_y [Hz] | f_z [Hz] | T_z [s] | d_st [mm] |
|---|---|---|---|---|---|---|
| self weight only | 110.076 | 129.983 | 129.983 | 113.61 | 0.00880207 | 0.0192456 |
| self weight + dead | 720.283 | 50.8138 | 50.8138 | 44.4129 | 0.022516 | 0.125934 |

- The x and y modes form a repeated pair: the five-fold symmetric cap is equally stiff in every horizontal direction (Kxx = Kyy, Kxy = 0).
- Check 1: the vertical frequency equals the single-DOF static-deflection formula f = sqrt(g/d_st)/(2 pi) with d_st = m g/Kzz.
- Check 2: the Rayleigh quotient of the static deflected shape under a 1 g horizontal load equals omega_x^2 (exact here because the shape is the mode).
- Check 3: each mode is a pure crown translation, so its effective modal mass is the whole free mass m in its own direction and 0 in the others (100 % participation). The modes are orthogonal with respect to M and K because both matrices are diagonal.

## F8  2V Class I geodesic geometry: chord factors and angles

Purpose. Fabrication data. The 2V Class I sphere splits each of the 20 icosahedron faces into four triangles through the edge midpoints, which are pushed out radially to the sphere. Every value below is computed twice, once from closed-form spherical trigonometry and once from the vertex vectors of an independently built 2V sphere, and the chord factors are compared with the 2V chord factors printed by the Domerama chord-factor chart (method 1, icosahedral): A = 0.54653, B = 0.61803. The Class I breakdown itself is described by Kenner (1976).

### Closed form

```text
icosahedron: central angle of an edge  cos th_ico = 1/sqrt5
  th_ico = 63.4349488229 deg,  1V chord = 2 sin(th_ico/2) = 1.051462224
strut A (5-valent vertex -> edge midpoint):  th_A = th_ico/2
  th_A = 31.7174744115 deg,  CF_A = 2 sin(th_A/2) = 0.546533057825
strut B (midpoint -> midpoint in one face):
  cos th_B = (v1+v2).(v1+v3)/|v1+v2|^2 = (1 + 3c)/(2 + 2c), c = 1/sqrt5
         = 0.809016994374948 = cos 36 deg = (1 + sqrt5)/4
  CF_B = 2 sin 18 deg = (sqrt5 - 1)/2 = 0.61803398875
axial angle (strut to the radius at its hub):
  (p_j - p_i).(-p_i) = R^2 - p_i.p_j = L^2/2  ->  cos alpha = L/(2R)
  alpha_A = acos(CF_A/2) = 74.14126279 deg,  alpha_B = acos(CF_B/2) = 72 deg
  (the end-cut or bevel angle to the tangent plane is 90 - alpha = th/2)
```

```text
panels: AAB (isosceles, 60 of them) and BBB (equilateral, 20)
  AAB apex angle (at the 5-valent hub) = acos(1 - B^2/(2A^2)) = 68.86197642 deg
  AAB base angles = 55.56901179 deg,  BBB = 60 deg
  areas (R = 1): AAB = (B/4) sqrt(4A^2 - B^2) = 0.1393,  BBB = (sqrt3/4) B^2 = 0.165396
dihedral angle across an edge (face circumradius rho = abc/(4 Area)):
  d_f = sqrt(R^2 - rho^2),  h_e = sqrt(R^2 - (L/2)^2)
  delta = asin(d_1/h_e) + asin(d_2/h_e)
  A edge (AAB | AAB): rho = 0.331309, d = 0.943522, h_e = 0.961938
    delta_A = 157.5410761 deg,  bend = 180 - delta = 22.45892392 deg
  B edge (AAB | BBB): rho = 0.331309, 0.356822, d = 0.943522, 0.934172, h_e = 0.951057
    delta_B = 161.9708979 deg,  bend = 18.02910211 deg
  anchor: the same formula for the icosahedron gives 138.189685104 deg,
  the exact value acos(-sqrt5/3) = 138.189685104 deg
```

### From the vertex vectors of an independently built 2V sphere

| quantity | vectors | closed form | difference |
|---|---|---|---|
| nodes V | 42 | 10 v^2 + 2 = 42 | 0 |
| struts E | 120 | 30 v^2 = 120 | 0 |
| panels F | 80 | 20 v^2 = 80 | 0 |
| Euler V - E + F | 2 | 2 | 0 |
| 5-valent hubs | 12 | 12 | 0 |
| 6-valent hubs | 30 | 30 | 0 |
| A struts | 60 | 60 | 0 |
| B struts | 60 | 60 | 0 |
| CF_A | 0.546533057825 | 0.546533057825 | 0 |
| CF_B | 0.61803398875 | 0.61803398875 | 2.22e-16 |
| axial A [deg] (max dev.) | 74.14126279 | 74.14126279 | 1.27e-14 |
| axial B [deg] (max dev.) | 72 | 72 | 2.54e-14 |
| dihedral A [deg] | 157.5410761 | 157.5410761 | 5.09e-14 |
| dihedral B [deg] | 161.9708979 | 161.9708979 | 2.54e-14 |
| panel area sum (R = 1) | 11.6659313917 | 11.6659313917 | 7.11e-15 |

Hub types: 30 x 2A+4B, 12 x 5A+0B. Panel types: 60 x AAB, 20 x BBB. Every strut has exactly two panels: True. All panels are oriented outwards: True. Largest distance of a node from the unit sphere: 2.22e-16. Spread of the dihedral angle inside each family: A 1.02e-13 deg, B 7.63e-14 deg. Faceted area / sphere area = 0.92834532.

- Domerama chart: A = 0.54653, B = 0.61803. Computed: A = 0.54653306, B = 0.61803399. Both agree to the five published decimals (differences 3.06e-06 and 3.99e-06).
- For a sphere of R = 5 m the struts are 2.73267 m (A) and 3.09017 m (B).
- The axial angle is the same at both ends of a strut because both hubs lie on the sphere. The dihedral angle is the angle between the two panels measured inside the dome. The bend angle (180 - dihedral) is the angle between the panel normals and is what a hub or panel bracket must turn.

## F9  Wind on validation Model B

Purpose. The wind load generator. Model B (1V dome on its lower ring, R = 5 m, 11 nodes, 15 equilateral panels) gets the program's default wind: an external pressure p = qz G Cp normal to the surface, lumped at the nodes on the surface tributary area A_i (one third of each adjacent panel) along the radial unit normal n_i. Cp > 0 pushes inward. The wind blows towards azimuth 0 deg (towards +x), so the windward springing is on the -x side.

```text
qz = 0.613 Kz Kzt Kd V^2            (ASCE 7-16 Eq. 26.10-1, SI)
   = 0.613 x 0.9 x 1 x 0.95 x 45^2 = 1061.332875 Pa
qz G = 1061.332875 x 0.85 = 902.1329437 Pa
theta_i = arccos(-n_i . d)   (0 = windward springing, 90 = crown, 180 = leeward)
Cp_i  = straight-line interpolation of the table 0:0.8, 30:0.2, 60:-0.6, 90:-1, 120:-0.8, 150:-0.5, 180:-0.4
F_i   = -qz G Cp_i A_i n_i
```

| node | A_trib [m2] | n | theta [deg] | Cp | F [kN] |
|---|---|---|---|---|---|
| 1 | 19.947 | (0, 0, 1) | 90 | -1 | (0, 0, 17.9948) |
| 2 | 19.947 | (0.894427, 0, 0.447214) | 153.43495 | -0.48855017 | (7.86324, 0, 3.93162) |
| 3 | 19.947 | (0.276393, 0.850651, 0.447214) | 106.04506 | -0.89303295 | (4.44163, 13.6699, 7.1867) |
| 4 | 19.947 | (-0.723607, 0.525731, 0.447214) | 43.646927 | -0.16391806 | (-2.1344, 1.55074, 1.31913) |
| 5 | 19.947 | (-0.723607, -0.525731, 0.447214) | 43.646927 | -0.16391806 | (-2.1344, -1.55074, 1.31913) |
| 6 | 19.947 | (0.276393, -0.850651, 0.447214) | 106.04506 | -0.89303295 | (4.44163, -13.6699, 7.1867) |
| 7 | 11.9682 | (0.723607, 0.525731, -0.447214) | 136.35307 | -0.63646927 | (4.97254, 3.61276, -3.0732) |
| 8 | 11.9682 | (-0.276393, 0.850651, -0.447214) | 73.954943 | -0.7860659 | (-2.34577, 7.21953, -3.79553) |
| 9 | 11.9682 | (-0.894427, 1.09536e-16, -0.447214) | 26.565051 | 0.26869898 | (2.59483, 0, 1.29742) |
| 10 | 11.9682 | (-0.276393, -0.850651, -0.447214) | 73.954943 | -0.7860659 | (-2.34577, -7.21953, -3.79553) |
| 11 | 11.9682 | (0.723607, -0.525731, -0.447214) | 136.35307 | -0.63646927 | (4.97254, -3.61276, -3.0732) |

- Velocity pressure qz = 1061.332875 Pa, qz G = 902.1329437 Pa.
- Crown: n = (0, 0, 1), theta = 90 deg, Cp = -1.00 (strong suction), F = qz G (5/3) A_face upward = 17.9948108 kN.
- Resultant (20.3261, 3.63798e-15, 26.4981) kN: drag towards +x and net uplift.
- Checks: the panel areas by cross product, the tributary areas sum to the surface, every node lies on the sphere, Fy of the resultant is zero and mirror nodes carry mirror forces.

## F10  AISC 360-16 compression strength of round HSS

Purpose. The member strength routines beyond the elastic branch that Model A reaches (its spoke has KL/r = 175): the inelastic E3 branch, the B4.2 design wall thickness 0.93 t for ERW HSS, the E7 effective area of a slender round HSS and the scope limit. ASTM A500 Gr. C: E = 200 GPa, Fy = 345 MPa, K = 1, phi_c = 0.9.

```text
KL/r,  Fe = pi^2 E / (KL/r)^2                       (E3-4)
Fcr = 0.658^(Fy/Fe) Fy   if KL/r <= 4.71 sqrt(E/Fy) (E3-2)
Fcr = 0.877 Fe           otherwise                  (E3-3)
D/t <= 0.11 E/Fy: Ae = Ag; up to 0.45 E/Fy:
   Ae/Ag = 0.038 E/(Fy D/t) + 2/3                    (E7-7)
phi_c Pn = 0.90 Fcr Ae                               (E3-1, E7-1)
4.71 sqrt(E/Fy) = 113.40347,  0.11 E/Fy = 63.768116,  0.45 E/Fy = 260.86957
```

| section and length | A [mm2] | r [mm] | KL/r | Fe [MPa] | Fcr [MPa] | D/t | Ae/Ag | phi Pn [kN] |
|---|---|---|---|---|---|---|---|---|
| CHS 88.9 x 4.0, L = 2 m | 1066.885 | 30.04998 | 66.55579 | 445.6133 | 249.5095 | 22.225 | 1 | 239.5781 |
| CHS 88.9 x 4.0 with t_des = 0.93 t, L = 2 m | 995.4752 | 30.14438 | 66.34735 | 448.4175 | 250.0156 | 23.8978 | 1 | 223.9959 |
| D 323.9 x 4.0, L = 3 m | 4019.982 | 113.1106 | 26.52272 | 2806.038 | 327.6952 | 80.975 | 0.9387134 | 1112.935 |
| D 323.9 x 1.2, L = 3 m | 1216.55 | 114.0925 | 26.29446 | 2854.968 | 327.9844 | 269.917 | 0 | 0 |

- Inelastic: KL/r = 66.5558 < 113.403, Fcr = 249.5095 MPa, phi Pn = 239.5781 kN.
- B4.2: the 0.93 t wall lowers A and r; phi Pn = 223.9959 kN, phi_t Pn = 309.0951 kN.
- E7: D/t = 80.975 > 0.11 E/Fy = 63.7681, Ae/Ag = 0.9387134, phi Pn = 1112.935 kN.
- Scope: D/t = 269.917 > 0.45 E/Fy = 260.87: outside the specification, reported as phi Pn = 0.

## Summary of independent checks

120 of 120 internal checks are within tolerance. A check compares two independent routes to the same number (equilibrium against compatibility, the force method against the displacement method, energy against work, closed form against vectors, or a published value).

| case | check | residual | tolerance | status |
|---|---|---|---|---|
| F1a | total elongation N L/EA + alpha dT L [m] | 2.17e-19 | 1e-15 | ok |
| F1a | R1 + R2 [N] | 0 | 1e-09 | ok |
| F1a | moment of reactions about the origin [N m] | 0 | 1e-09 | ok |
| F1a | R1 (equilibrium) - R1 (stiffness method) [N] | 0 | 1e-09 | ok |
| F1b | u2 force method - u2 stiffness method [m] | 8.13e-20 | 1e-15 | ok |
| F1b | segment 2 elongation + u2 [m] | 1.08e-19 | 1e-15 | ok |
| F2a | crown uz compatibility - stiffness [m] | 4.34e-19 | 1e-15 | ok |
| F2a | crown ux, uy by stiffness [m] | 3.13e-19 | 1e-18 | ok |
| F2a | crown equilibrium [N] | 5.41e-11 | 1e-08 | ok |
| F2a | sum of reactions [N] | 6.27e-11 | 1e-08 | ok |
| F2a | moment of reactions [N m] | 8.98e-11 | 1e-07 | ok |
| F2a | radial reaction joint sum - formula [N] | 5.82e-11 | 1e-08 | ok |
| F2b | spoke forces vs closed form [-] | 2.22e-16 | 1e-13 | ok |
| F2b | crown equilibrium [N] | 1.51e-11 | 1e-08 | ok |
| F2b | sum of reactions [N] | 1.51e-11 | 1e-08 | ok |
| F2b | sum of spoke forces (z equilibrium) [N] | 7.28e-12 | 1e-08 | ok |
| F3a | crown displacement vs closed form [m] | 2.17e-19 | 1e-17 | ok |
| F3a | spoke forces vs closed form [-] | 5.55e-17 | 1e-13 | ok |
| F3a | crown equilibrium [N] | 1.42e-11 | 1e-08 | ok |
| F3a | sum of reactions [N] | 1.42e-11 | 1e-08 | ok |
| F3a | moment of reactions [N m] | 4.12e-11 | 1e-07 | ok |
| F3a | energy: U - 1/2 R.us [J] | 7.11e-15 | 1e-09 | ok |
| F3b | crown displacement vs closed form [m] | 3.25e-19 | 1e-17 | ok |
| F3b | spoke forces vs closed form [-] | 1.67e-16 | 1e-13 | ok |
| F3b | ring forces vs closed form [N] | 1.46e-11 | 1e-08 | ok |
| F3b | crown equilibrium [N] | 6.2e-11 | 1e-08 | ok |
| F3b | sum of reactions [N] | 6.08e-11 | 1e-08 | ok |
| F3b | moment of reactions [N m] | 1.65e-10 | 1e-07 | ok |
| F3b | energy: U - 1/2 R.us [J] | 5.68e-14 | 1e-09 | ok |
| F4 load | compatibility of all bars [m] | 2.17e-19 | 1e-15 | ok |
| F4 load | equilibrium of joint D [N] | 1.82e-12 | 1e-08 | ok |
| F4 load | sum of reactions + load [N] | 1.82e-12 | 1e-08 | ok |
| F4 load | moment about B [N m] | 0 | 1e-08 | ok |
| F4 load | N force method - displacement method [N] | 7.28e-12 | 1e-08 | ok |
| F4 load | u force method - displacement method [m] | 3.35e-19 | 1e-16 | ok |
| F4 temperature | compatibility of all bars [m] | 0 | 1e-15 | ok |
| F4 temperature | equilibrium of joint D [N] | 3.64e-12 | 1e-08 | ok |
| F4 temperature | sum of reactions + load [N] | 3.64e-12 | 1e-08 | ok |
| F4 temperature | moment about B [N m] | 7.28e-12 | 1e-08 | ok |
| F4 temperature | N force method - displacement method [N] | 1.09e-11 | 1e-08 | ok |
| F4 temperature | u force method - displacement method [m] | 6.78e-20 | 1e-16 | ok |
| F4 settlement | compatibility of all bars [m] | 2.17e-19 | 1e-15 | ok |
| F4 settlement | equilibrium of joint D [N] | 7.28e-12 | 1e-08 | ok |
| F4 settlement | sum of reactions + load [N] | 7.28e-12 | 1e-08 | ok |
| F4 settlement | moment about B [N m] | 1.46e-11 | 1e-08 | ok |
| F4 settlement | N force method - displacement method [N] | 1.46e-11 | 1e-08 | ok |
| F4 settlement | u force method - displacement method [m] | 4.85e-19 | 1e-16 | ok |
| F4 combined | compatibility of all bars [m] | 3.25e-19 | 1e-15 | ok |
| F4 combined | equilibrium of joint D [N] | 1.46e-11 | 1e-08 | ok |
| F4 combined | sum of reactions + load [N] | 1.46e-11 | 1e-08 | ok |
| F4 combined | moment about B [N m] | 1.46e-11 | 1e-08 | ok |
| F4 combined | N force method - displacement method [N] | 1.46e-11 | 1e-08 | ok |
| F4 combined | u force method - displacement method [m] | 4.34e-19 | 1e-16 | ok |
| F4 | superposition of X1 [N] | 0 | 1e-08 | ok |
| F5a | Clapeyron 1/2 P u2 - (bar + spring energy) [J] | 0 | 1e-12 | ok |
| F5a | load + spring reaction [N] | 0 | 1e-09 | ok |
| F5b | load + reactions [N] | 3.64e-12 | 1e-09 | ok |
| F5b | Clapeyron [J] | 0 | 1e-12 | ok |
| F5c | N from u2 - closed form [N] | 7.28e-12 | 1e-08 | ok |
| F5c | node 2 equilibrium -N + spring [N] | 3.64e-12 | 1e-08 | ok |
| F5c | k -> infinity gives F1 [relative] | 7.22e-15 | 1e-13 | ok |
| F5c | k -> 0 gives free expansion [relative] | 1.4e-14 | 1e-13 | ok |
| F6 | dP/dw(0) analytic - finite difference [relative] | 3.15e-13 | 1e-06 | ok |
| F6 | dP/dw(0) - 2 EA h^2/L^3 [relative] | 2.22e-16 | 1e-12 | ok |
| F6 | Newton y* - closed form y* [relative] | 0 | 8.9e-16 | ok |
| F6 | l* = a sqrt(1+z*) vs (a^2 L)^(1/3) [relative] | 0 | 1e-14 | ok |
| F6 | dP/dw at the limit point [relative to K0] | 6.85e-17 | 1e-12 | ok |
| F6 | 2EA(1/L - a^2/l*^3) at the limit point [relative to K0] | 1.11e-13 | 1e-12 | ok |
| F6 | P* = -2 N* y*/l* [N] | 0 | 1e-07 | ok |
| F6 | P(w_snap) - P* [N] | 2.55e-11 | 1e-06 | ok |
| F6 | P(h) [N] | 0 | 1e-09 | ok |
| F6 | P(2h) [N] | 0 | 1e-06 | ok |
| F6 | P(h + y*) + P* [N] | 0 | 1e-06 | ok |
| F6 | lambda_cr P_ref - 2EA h^3/(a^2 L) [relative] | 0 | 1e-13 | ok |
| F6 | no sway bifurcation before the limit (Kuu(w*) > 0) | 0 | 0.5 | ok |
| F7 sw | f_z vs static deflection sqrt(g/d_st) [rel] | 0 | 1e-14 | ok |
| F7 sw | Rayleigh quotient vs omega_x^2 [rel] | 4.44e-16 | 1e-14 | ok |
| F7 sw | f_z/f_xy - sqrt2 t/s [-] | 2.22e-16 | 1e-13 | ok |
| F7 dead | f_z vs static deflection sqrt(g/d_st) [rel] | 0 | 1e-14 | ok |
| F7 dead | Rayleigh quotient vs omega_x^2 [rel] | 3.33e-16 | 1e-14 | ok |
| F7 dead | f_z/f_xy - sqrt2 t/s [-] | 2.22e-16 | 1e-13 | ok |
| F8 | CF_B - (sqrt5 - 1)/2 | 3.33e-16 | 1e-15 | ok |
| F8 | cos th_B - cos 36 | 1.11e-16 | 1e-15 | ok |
| F8 | dihedral formula on the icosahedron - acos(-sqrt5/3) [rad] | 4.44e-16 | 1e-13 | ok |
| F8 | AAB triangle acute (circumcentre inside) | 0 | 0.5 | ok |
| F8 | V, E, F counts | 0 | 0 | ok |
| F8 | unclassified struts | 0 | 0 | ok |
| F8 | valences 12 x 5 and 30 x 6 | 0 | 0 | ok |
| F8 | A and B strut counts 60/60 | 0 | 0 | ok |
| F8 | axial angle vectors - acos(L/2R) [rad] | 4.44e-16 | 1e-12 | ok |
| F8 | dihedral A vectors - formula [rad] | 8.88e-16 | 1e-12 | ok |
| F8 | dihedral B vectors - formula [rad] | 4.44e-16 | 1e-12 | ok |
| F8 | dihedral spread inside a family [rad] | 1.78e-15 | 1e-12 | ok |
| F8 | panel area sum vectors - formula | 7.11e-15 | 1e-12 | ok |
| F8 | CF_A vs Domerama 0.54653 (5 decimals) | 3.06e-06 | 5e-06 | ok |
| F8 | CF_B vs Domerama 0.61803 (5 decimals) | 3.99e-06 | 5e-06 | ok |
| F8 | every strut shared by two panels | 0 | 0.5 | ok |
| F9 | 15 panel areas by cross product - 15 sqrt3/4 L^2 [m2] | 0 | 1e-11 | ok |
| F9 | sum of tributary areas - surface [m2] | 0 | 1e-11 | ok |
| F9 | node 1 on the sphere \|x - c\| - R [m] | 0 | 1e-12 | ok |
| F9 | node 2 on the sphere \|x - c\| - R [m] | 0 | 1e-12 | ok |
| F9 | node 3 on the sphere \|x - c\| - R [m] | 0 | 1e-12 | ok |
| F9 | node 4 on the sphere \|x - c\| - R [m] | 0 | 1e-12 | ok |
| F9 | node 5 on the sphere \|x - c\| - R [m] | 0 | 1e-12 | ok |
| F9 | node 6 on the sphere \|x - c\| - R [m] | 0 | 1e-12 | ok |
| F9 | node 7 on the sphere \|x - c\| - R [m] | 8.88e-16 | 1e-12 | ok |
| F9 | node 8 on the sphere \|x - c\| - R [m] | 0 | 1e-12 | ok |
| F9 | node 9 on the sphere \|x - c\| - R [m] | 0 | 1e-12 | ok |
| F9 | node 10 on the sphere \|x - c\| - R [m] | 0 | 1e-12 | ok |
| F9 | node 11 on the sphere \|x - c\| - R [m] | 0 | 1e-12 | ok |
| F9 | crown Cp = -1.00 exactly | 0 | 0 | ok |
| F9 | resultant Fy = 0 (mirror symmetry about x-z) [N] | 3.64e-12 | 1e-09 | ok |
| F9 | mirror pair 3/6: Fx, Fz equal, Fy opposite [N] | 2.73e-12 | 1e-09 | ok |
| F9 | mirror pair 4/5: Fx, Fz equal, Fy opposite [N] | 4.55e-12 | 1e-09 | ok |
| F9 | mirror pair 7/11: Fx, Fz equal, Fy opposite [N] | 4.55e-13 | 1e-09 | ok |
| F9 | mirror pair 8/10: Fx, Fz equal, Fy opposite [N] | 1.36e-12 | 1e-09 | ok |
| F10 | inelastic branch: KL/r below 4.71 sqrt(E/Fy) | 0 | 0 | ok |
| F10 | E7 case: D/t between 0.11 and 0.45 E/Fy | 0 | 0 | ok |
| F10 | scope case: D/t above 0.45 E/Fy | 0 | 0 | ok |
| F10 | E7: 0.038/lambda + 2/3 - Q [-] | 0 | 1e-15 | ok |

## Reference values (results() dictionary)

All values in SI units. These are the numbers the program must reproduce.

| key | value | unit | formula | description |
|---|---|---|---|---|
| F1.L | 3 | m | L = \|x2 - x1\| | length of the inclined bar |
| F1.N | -74895.31753 | N | N = -E*A*alpha*dT | axial force in the heated fixed-fixed bar |
| F1.sigma | -70200000 | Pa | sigma = -E*alpha*dT | axial stress in the heated fixed-fixed bar |
| F1.R1_x | 49930.21169 | N | R1 = -N*c | reaction x at node 1 (0,0,0) |
| F1.R2_x | -49930.21169 | N | R2 = +N*c | reaction x at node 2 (2,1,2) |
| F1.R1_y | 24965.10584 | N | R1 = -N*c | reaction y at node 1 (0,0,0) |
| F1.R2_y | -24965.10584 | N | R2 = +N*c | reaction y at node 2 (2,1,2) |
| F1.R1_z | 49930.21169 | N | R1 = -N*c | reaction z at node 1 (0,0,0) |
| F1.R2_z | -49930.21169 | N | R2 = +N*c | reaction z at node 2 (2,1,2) |
| F1.free_elongation | 0.001053 | m | e_T = alpha*dT*L | elongation if the bar were free to expand |
| F1.dT_euler | 84.63689564 | degC | dT_cr = pi^2*I/(alpha*A*L^2) | temperature rise that buckles the 3 m bar (information) |
| F1.stepped_N | -69599.82832 | N | N = -alpha*dT*(L1+L2)/(L1/EA1 + L2/EA2) | axial force in both segments of the heated stepped bar |
| F1.stepped_u2 | 4.9635058e-05 | m | u2 = (N/EA1 + alpha*dT)*L1 | x displacement of the middle node of the stepped bar |
| F1.stepped_R1x | 69599.82832 | N | R1x = -N | x reaction at node 1 of the stepped bar |
| F2.crown_uz | 0.00351 | m | w = alpha*dT*L^2/H | crown rise of Model A, all struts heated |
| F2.crown_ux | 0 | m | symmetry | crown x movement (zero) |
| F2.N_spoke | 0 | N | N = k*t*w - EA*alpha*dT = 0 | force in each spoke (stress free) |
| F2.N_ring | -74895.31753 | N | N = -E*A*alpha*dT | force in each ring strut (fixed-fixed) |
| F2.R_radial | -88044.72622 | N | R_r = -2*sin(36deg)*E*A*alpha*dT | radial reaction at every support node (+ outward) |
| F2.R_z | 0 | N | spokes carry no force | vertical reaction at every support node |
| F2.R2_x | -88044.72622 | N | R = -2*sin36*EA*alpha*dT*e_r | reaction x at hand node 2, all struts heated |
| F2.R2_y | 0 | N | R = -2*sin36*EA*alpha*dT*e_r | reaction y at hand node 2, all struts heated |
| F2.R2_z | 0 | N | R = -2*sin36*EA*alpha*dT*e_r | reaction z at hand node 2, all struts heated |
| F2.R3_x | -27207.31667 | N | R = -2*sin36*EA*alpha*dT*e_r | reaction x at hand node 3, all struts heated |
| F2.R3_y | -83735.5106 | N | R = -2*sin36*EA*alpha*dT*e_r | reaction y at hand node 3, all struts heated |
| F2.R3_z | 0 | N | R = -2*sin36*EA*alpha*dT*e_r | reaction z at hand node 3, all struts heated |
| F2.R4_x | 71229.67978 | N | R = -2*sin36*EA*alpha*dT*e_r | reaction x at hand node 4, all struts heated |
| F2.R4_y | -51751.39162 | N | R = -2*sin36*EA*alpha*dT*e_r | reaction y at hand node 4, all struts heated |
| F2.R4_z | 0 | N | R = -2*sin36*EA*alpha*dT*e_r | reaction z at hand node 4, all struts heated |
| F2.R5_x | 71229.67978 | N | R = -2*sin36*EA*alpha*dT*e_r | reaction x at hand node 5, all struts heated |
| F2.R5_y | 51751.39162 | N | R = -2*sin36*EA*alpha*dT*e_r | reaction y at hand node 5, all struts heated |
| F2.R5_z | 0 | N | R = -2*sin36*EA*alpha*dT*e_r | reaction z at hand node 5, all struts heated |
| F2.R6_x | -27207.31667 | N | R = -2*sin36*EA*alpha*dT*e_r | reaction x at hand node 6, all struts heated |
| F2.R6_y | 83735.5106 | N | R = -2*sin36*EA*alpha*dT*e_r | reaction y at hand node 6, all struts heated |
| F2.R6_z | 0 | N | R = -2*sin36*EA*alpha*dT*e_r | reaction z at hand node 6, all struts heated |
| F2.oneSpoke_crown_ux | -0.0008677197202 | m | u = Kff^-1 * EA*alpha*dT*c_1 | crown x displacement, only spoke 1 heated |
| F2.oneSpoke_crown_uy | 0 | m | u = Kff^-1 * EA*alpha*dT*c_1 | crown y displacement, only spoke 1 heated |
| F2.oneSpoke_crown_uz | 0.000702 | m | u = Kff^-1 * EA*alpha*dT*c_1 | crown z displacement, only spoke 1 heated |
| F2.oneSpoke_N1 | -29958.12701 | N | N1 = -0.4*EA*alpha*dT, Ni = EA*alpha*dT*(0.4cos(th_i)+0.2) | force in spoke 1, only spoke 1 heated |
| F2.oneSpoke_N2 | 24236.63387 | N | N1 = -0.4*EA*alpha*dT, Ni = EA*alpha*dT*(0.4cos(th_i)+0.2) | force in spoke 2, only spoke 1 heated |
| F2.oneSpoke_N3 | -9257.570367 | N | N1 = -0.4*EA*alpha*dT, Ni = EA*alpha*dT*(0.4cos(th_i)+0.2) | force in spoke 3, only spoke 1 heated |
| F2.oneSpoke_N4 | -9257.570367 | N | N1 = -0.4*EA*alpha*dT, Ni = EA*alpha*dT*(0.4cos(th_i)+0.2) | force in spoke 4, only spoke 1 heated |
| F2.oneSpoke_N5 | 24236.63387 | N | N1 = -0.4*EA*alpha*dT, Ni = EA*alpha*dT*(0.4cos(th_i)+0.2) | force in spoke 5, only spoke 1 heated |
| F3.vert_crown_ux | 0.001236067977 | m | Kff*u = k*c1*(c1.us); u = (t*d/(2.5s), 0, -d/5) | crown x displacement, node 2 settles 5 mm |
| F3.vert_crown_uy | 0 | m | Kff*u = k*c1*(c1.us); u = (t*d/(2.5s), 0, -d/5) | crown y displacement, node 2 settles 5 mm |
| F3.vert_crown_uz | -0.001 | m | Kff*u = k*c1*(c1.us); u = (t*d/(2.5s), 0, -d/5) | crown z displacement, node 2 settles 5 mm |
| F3.vert_N1 | 42675.39461 | N | N1 = 0.4 k t d, Ni = -k t d (0.4 cos th_i + 0.2), ring 0 | force in member 1, node 2 settles |
| F3.vert_N2 | -34525.11948 | N | N1 = 0.4 k t d, Ni = -k t d (0.4 cos th_i + 0.2), ring 0 | force in member 2, node 2 settles |
| F3.vert_N3 | 13187.42218 | N | N1 = 0.4 k t d, Ni = -k t d (0.4 cos th_i + 0.2), ring 0 | force in member 3, node 2 settles |
| F3.vert_N4 | 13187.42218 | N | N1 = 0.4 k t d, Ni = -k t d (0.4 cos th_i + 0.2), ring 0 | force in member 4, node 2 settles |
| F3.vert_N5 | -34525.11948 | N | N1 = 0.4 k t d, Ni = -k t d (0.4 cos th_i + 0.2), ring 0 | force in member 5, node 2 settles |
| F3.vert_N6 | 0 | N | N1 = 0.4 k t d, Ni = -k t d (0.4 cos th_i + 0.2), ring 0 | force in member 6, node 2 settles |
| F3.vert_N7 | 0 | N | N1 = 0.4 k t d, Ni = -k t d (0.4 cos th_i + 0.2), ring 0 | force in member 7, node 2 settles |
| F3.vert_N8 | 0 | N | N1 = 0.4 k t d, Ni = -k t d (0.4 cos th_i + 0.2), ring 0 | force in member 8, node 2 settles |
| F3.vert_N9 | 0 | N | N1 = 0.4 k t d, Ni = -k t d (0.4 cos th_i + 0.2), ring 0 | force in member 9, node 2 settles |
| F3.vert_N10 | 0 | N | N1 = 0.4 k t d, Ni = -k t d (0.4 cos th_i + 0.2), ring 0 | force in member 10, node 2 settles |
| F3.vert_R2_x | 36301.85892 | N | R = -(sum of member forces on the node) | reaction x at hand node 2, node 2 settles |
| F3.vert_R2_y | 0 | N | R = -(sum of member forces on the node) | reaction y at hand node 2, node 2 settles |
| F3.vert_R2_z | -22435.78267 | N | R = -(sum of member forces on the node) | reaction z at hand node 2, node 2 settles |
| F3.vert_R3_x | -9075.46473 | N | R = -(sum of member forces on the node) | reaction x at hand node 3, node 2 settles |
| F3.vert_R3_y | -27931.40839 | N | R = -(sum of member forces on the node) | reaction y at hand node 3, node 2 settles |
| F3.vert_R3_z | 18150.92946 | N | R = -(sum of member forces on the node) | reaction z at hand node 3, node 2 settles |
| F3.vert_R4_x | -9075.46473 | N | R = -(sum of member forces on the node) | reaction x at hand node 4, node 2 settles |
| F3.vert_R4_y | 6593.711088 | N | R = -(sum of member forces on the node) | reaction y at hand node 4, node 2 settles |
| F3.vert_R4_z | -6933.038126 | N | R = -(sum of member forces on the node) | reaction z at hand node 4, node 2 settles |
| F3.vert_R5_x | -9075.46473 | N | R = -(sum of member forces on the node) | reaction x at hand node 5, node 2 settles |
| F3.vert_R5_y | -6593.711088 | N | R = -(sum of member forces on the node) | reaction y at hand node 5, node 2 settles |
| F3.vert_R5_z | -6933.038126 | N | R = -(sum of member forces on the node) | reaction z at hand node 5, node 2 settles |
| F3.vert_R6_x | -9075.46473 | N | R = -(sum of member forces on the node) | reaction x at hand node 6, node 2 settles |
| F3.vert_R6_y | 27931.40839 | N | R = -(sum of member forces on the node) | reaction y at hand node 6, node 2 settles |
| F3.vert_R6_z | 18150.92946 | N | R = -(sum of member forces on the node) | reaction z at hand node 6, node 2 settles |
| F3.vert_energy | 56.08945667 | J | U = sum N^2 L/(2EA) = 1/2 R.us | strain energy, node 2 settles |
| F3.spread_crown_ux | 0.0006180339887 | m | u = 0.4*D*e_r(72) + (0,0,-s*D/(5t)) | crown x displacement, node 3 spreads 5 mm |
| F3.spread_crown_uy | 0.001902113033 | m | u = 0.4*D*e_r(72) + (0,0,-s*D/(5t)) | crown y displacement, node 3 spreads 5 mm |
| F3.spread_crown_uz | -0.001618033989 | m | u = 0.4*D*e_r(72) + (0,0,-s*D/(5t)) | crown z displacement, node 3 spreads 5 mm |
| F3.spread_N1 | -55862.81678 | N | N2 = 0.4 k s D, ring 6,7 = k D sin36 | force in member 1, node 3 spreads |
| F3.spread_N2 | 69050.23896 | N | N2 = 0.4 k s D, ring 6,7 = k D sin36 | force in member 2, node 3 spreads |
| F3.spread_N3 | -55862.81678 | N | N2 = 0.4 k s D, ring 6,7 = k D sin36 | force in member 3, node 3 spreads |
| F3.spread_N4 | 21337.6973 | N | N2 = 0.4 k s D, ring 6,7 = k D sin36 | force in member 4, node 3 spreads |
| F3.spread_N5 | 21337.6973 | N | N2 = 0.4 k s D, ring 6,7 = k D sin36 | force in member 5, node 3 spreads |
| F3.spread_N6 | 119281.3541 | N | N2 = 0.4 k s D, ring 6,7 = k D sin36 | force in member 6, node 3 spreads |
| F3.spread_N7 | 119281.3541 | N | N2 = 0.4 k s D, ring 6,7 = k D sin36 | force in member 7, node 3 spreads |
| F3.spread_N8 | 0 | N | N2 = 0.4 k s D, ring 6,7 = k D sin36 | force in member 8, node 3 spreads |
| F3.spread_N9 | 0 | N | N2 = 0.4 k s D, ring 6,7 = k D sin36 | force in member 9, node 3 spreads |
| F3.spread_N10 | 0 | N | N2 = 0.4 k s D, ring 6,7 = k D sin36 | force in member 10, node 3 spreads |
| F3.spread_R2_x | 22592.07058 | N | R = -(sum of member forces on the node) | reaction x at hand node 2, node 3 spreads |
| F3.spread_R2_y | -96500.64261 | N | R = -(sum of member forces on the node) | reaction y at hand node 2, node 3 spreads |
| F3.spread_R2_z | 29368.82079 | N | R = -(sum of member forces on the node) | reaction z at hand node 2, node 3 spreads |
| F3.spread_R3_x | 61482.41775 | N | R = -(sum of member forces on the node) | reaction x at hand node 3, node 3 spreads |
| F3.spread_R3_y | 189223.4249 | N | R = -(sum of member forces on the node) | reaction y at hand node 3, node 3 spreads |
| F3.spread_R3_z | -36301.85892 | N | R = -(sum of member forces on the node) | reaction z at hand node 3, node 3 spreads |
| F3.spread_R4_x | -74999.0236 | N | R = -(sum of member forces on the node) | reaction x at hand node 4, node 3 spreads |
| F3.spread_R4_y | -64791.37393 | N | R = -(sum of member forces on the node) | reaction y at hand node 4, node 3 spreads |
| F3.spread_R4_z | 29368.82079 | N | R = -(sum of member forces on the node) | reaction z at hand node 4, node 3 spreads |
| F3.spread_R5_x | -14684.4104 | N | R = -(sum of member forces on the node) | reaction x at hand node 5, node 3 spreads |
| F3.spread_R5_y | -10668.84865 | N | R = -(sum of member forces on the node) | reaction y at hand node 5, node 3 spreads |
| F3.spread_R5_z | -11217.89133 | N | R = -(sum of member forces on the node) | reaction z at hand node 5, node 3 spreads |
| F3.spread_R6_x | 5608.945667 | N | R = -(sum of member forces on the node) | reaction x at hand node 6, node 3 spreads |
| F3.spread_R6_y | -17262.55974 | N | R = -(sum of member forces on the node) | reaction y at hand node 6, node 3 spreads |
| F3.spread_R6_z | -11217.89133 | N | R = -(sum of member forces on the node) | reaction z at hand node 6, node 3 spreads |
| F3.spread_energy | 497.4032081 | J | U = sum N^2 L/(2EA) = 1/2 R.us | strain energy, node 3 spreads |
| F4.f11 | 4.149612784e-08 | m/N | f11 = sum n1_i^2 L_i/(E A_i) | flexibility coefficient of the redundant N2 |
| F4.N1_load | -6956.712954 | N | N_i = n0_i + X1*n1_i, X1 = -(D1P+D1T+D1S)/f11 | force in bar 1 (A-D), case load |
| F4.N2_load | -35260.60186 | N | N_i = n0_i + X1*n1_i, X1 = -(D1P+D1T+D1S)/f11 | force in bar 2 (B-D), case load |
| F4.N3_load | -29083.82568 | N | N_i = n0_i + X1*n1_i, X1 = -(D1P+D1T+D1S)/f11 | force in bar 3 (C-D), case load |
| F4.ux_load | 0.0002340532636 | m | unit-load method on the primary structure | x displacement of node D, case load |
| F4.uz_load | -0.0005837620405 | m | unit-load method on the primary structure | z displacement of node D, case load |
| F4.RAx_load | 5565.370364 | N | R_k = -N_i*c_i | reaction x at support A, case load |
| F4.RAz_load | 4174.027773 | N | R_k = -N_i*c_i | reaction z at support A, case load |
| F4.RBx_load | 0 | N | R_k = -N_i*c_i | reaction x at support B, case load |
| F4.RBz_load | 35260.60186 | N | R_k = -N_i*c_i | reaction z at support B, case load |
| F4.RCx_load | -20565.37036 | N | R_k = -N_i*c_i | reaction x at support C, case load |
| F4.RCz_load | 20565.37036 | N | R_k = -N_i*c_i | reaction z at support C, case load |
| F4.N1_temperature | 18125.61548 | N | N_i = n0_i + X1*n1_i, X1 = -(D1P+D1T+D1S)/f11 | force in bar 1 (A-D), case temperature |
| F4.N2_temperature | -25375.86167 | N | N_i = n0_i + X1*n1_i, X1 = -(D1P+D1T+D1S)/f11 | force in bar 2 (B-D), case temperature |
| F4.N3_temperature | 20506.79299 | N | N_i = n0_i + X1*n1_i, X1 = -(D1P+D1T+D1S)/f11 | force in bar 3 (C-D), case temperature |
| F4.ux_temperature | 5.625062267e-05 | m | unit-load method on the primary structure | x displacement of node D, case temperature |
| F4.uz_temperature | 0.0006328862184 | m | unit-load method on the primary structure | z displacement of node D, case temperature |
| F4.RAx_temperature | -14500.49239 | N | R_k = -N_i*c_i | reaction x at support A, case temperature |
| F4.RAz_temperature | -10875.36929 | N | R_k = -N_i*c_i | reaction z at support A, case temperature |
| F4.RBx_temperature | 0 | N | R_k = -N_i*c_i | reaction x at support B, case temperature |
| F4.RBz_temperature | 25375.86167 | N | R_k = -N_i*c_i | reaction z at support B, case temperature |
| F4.RCx_temperature | 14500.49239 | N | R_k = -N_i*c_i | reaction x at support C, case temperature |
| F4.RCz_temperature | -14500.49239 | N | R_k = -N_i*c_i | reaction z at support C, case temperature |
| F4.N1_settlement | 49180.88585 | N | N_i = n0_i + X1*n1_i, X1 = -(D1P+D1T+D1S)/f11 | force in bar 1 (A-D), case settlement |
| F4.N2_settlement | -68853.2402 | N | N_i = n0_i + X1*n1_i, X1 = -(D1P+D1T+D1S)/f11 | force in bar 2 (B-D), case settlement |
| F4.N3_settlement | 55641.82063 | N | N_i = n0_i + X1*n1_i, X1 = -(D1P+D1T+D1S)/f11 | force in bar 3 (C-D), case settlement |
| F4.ux_settlement | 0.002295483985 | m | unit-load method on the primary structure | x displacement of node D, case settlement |
| F4.uz_settlement | -0.001139909867 | m | unit-load method on the primary structure | z displacement of node D, case settlement |
| F4.RAx_settlement | -39344.70868 | N | R_k = -N_i*c_i | reaction x at support A, case settlement |
| F4.RAz_settlement | -29508.53151 | N | R_k = -N_i*c_i | reaction z at support A, case settlement |
| F4.RBx_settlement | 0 | N | R_k = -N_i*c_i | reaction x at support B, case settlement |
| F4.RBz_settlement | 68853.2402 | N | R_k = -N_i*c_i | reaction z at support B, case settlement |
| F4.RCx_settlement | 39344.70868 | N | R_k = -N_i*c_i | reaction x at support C, case settlement |
| F4.RCz_settlement | -39344.70868 | N | R_k = -N_i*c_i | reaction z at support C, case settlement |
| F4.N1_combined | 60349.78838 | N | N_i = n0_i + X1*n1_i, X1 = -(D1P+D1T+D1S)/f11 | force in bar 1 (A-D), case combined |
| F4.N2_combined | -129489.7037 | N | N_i = n0_i + X1*n1_i, X1 = -(D1P+D1T+D1S)/f11 | force in bar 2 (B-D), case combined |
| F4.N3_combined | 47064.78794 | N | N_i = n0_i + X1*n1_i, X1 = -(D1P+D1T+D1S)/f11 | force in bar 3 (C-D), case combined |
| F4.ux_combined | 0.002585787871 | m | unit-load method on the primary structure | x displacement of node D, case combined |
| F4.uz_combined | -0.00109078569 | m | unit-load method on the primary structure | z displacement of node D, case combined |
| F4.RAx_combined | -48279.83071 | N | R_k = -N_i*c_i | reaction x at support A, case combined |
| F4.RAz_combined | -36209.87303 | N | R_k = -N_i*c_i | reaction z at support A, case combined |
| F4.RBx_combined | 0 | N | R_k = -N_i*c_i | reaction x at support B, case combined |
| F4.RBz_combined | 129489.7037 | N | R_k = -N_i*c_i | reaction z at support B, case combined |
| F4.RCx_combined | 33279.83071 | N | R_k = -N_i*c_i | reaction x at support C, case combined |
| F4.RCz_combined | -33279.83071 | N | R_k = -N_i*c_i | reaction z at support C, case combined |
| F5.series_u1 | 0.001 | m | u1 = P/k | spring node displacement |
| F5.series_u2 | 0.001702981197 | m | u2 = P/k + P*L/EA | loaded node displacement, bar and spring in series |
| F5.series_N | 50000 | N | N = P | bar force, series |
| F5.series_spring_reaction | -50000 | N | R = -k*u1 = -P | spring reaction at node 1, series |
| F5.parallel_u2 | 0.0004127944563 | m | u2 = P/(EA/L + k) | loaded node displacement, bar and spring in parallel |
| F5.parallel_N | 29360.27719 | N | N = (EA/L)*u2 | bar force, parallel |
| F5.parallel_spring_reaction | -20639.72281 | N | R = -k*u2 | spring reaction at node 2, parallel |
| F5.parallel_R1x | -29360.27719 | N | R1 = -N | fixed reaction at node 1, parallel |
| F5.thermal_u2 | 0.0006183274376 | m | u2 = alpha*dT*L/(1 + k*L/EA) | spring node displacement, heated bar |
| F5.thermal_N | -30916.37188 | N | N = -EA*alpha*dT*k/(k + EA/L) | bar force, heated bar on a spring |
| F5.thermal_spring_reaction | -30916.37188 | N | R = -k*u2 = N | spring reaction, heated bar on a spring |
| F6.P_at_w0.1h | 4545.497953 | N | P = 2EA(L-l)/L*(h-w)/l | load on the path at w = 0.1 h |
| F6.P_at_w0.25h | 8726.21398 | N | P = 2EA(L-l)/L*(h-w)/l | load on the path at w = 0.25 h |
| F6.P_at_w0.5h | 9978.654483 | N | P = 2EA(L-l)/L*(h-w)/l | load on the path at w = 0.5 h |
| F6.P_at_w1h | 0 | N | P = 2EA(L-l)/L*(h-w)/l | load on the path at w = 1 h |
| F6.P_at_w1.5h | -9978.654483 | N | P = 2EA(L-l)/L*(h-w)/l | load on the path at w = 1.5 h |
| F6.P_at_w2h | 0 | N | P = 2EA(L-l)/L*(h-w)/l | load on the path at w = 2 h |
| F6.P_at_w2.5h | 49707.09372 | N | P = 2EA(L-l)/L*(h-w)/l | load on the path at w = 2.5 h |
| F6.L | 2.503123049 | m | L = sqrt(a^2 + h^2) | initial bar length |
| F6.EA | 213376973 | N | EA | axial rigidity of each bar |
| F6.K0 | 425158.6052 | N/m | K0 = 2EA h^2/L^3 | initial vertical stiffness of the apex |
| F6.l_limit | 2.501040583 | m | l* = (a^2 L)^(1/3) | bar length at the limit point |
| F6.w_limit | 0.05286125123 | m | w* = h - sqrt(l*^2 - a^2) | apex deflection at the limit point |
| F6.P_limit | 10240.4996 | N | P* = 2EA(L-l*)/L*(h-w*)/l* | limit (snap-through) load |
| F6.N_limit | -177518.3623 | N | N* = EA(l*-L)/L | bar force at the limit point |
| F6.strain_limit | -0.0008319471394 | - | eps* = (l*-L)/L | engineering strain at the limit point |
| F6.w_snap | 0.2693676147 | m | P(w_snap) = P*, w_snap > 2h | apex deflection on the stable inverted branch at P* |
| F6.N_snap | 88814.58642 | N | N = EA(l-L)/L at w_snap | bar force after snap-through at P* |
| F6.w_second_limit | 0.1971387488 | m | w = h + y* | apex deflection at the second limit point (P = -P*) |
| F6.P_ref | 10000 | N | stated | reference load of the eigenproblem |
| F6.N_ref | -100124.922 | N | N_ref = -P_ref*L/(2h) | linear bar force under P_ref |
| F6.lambda_cr | 5.327768772 | - | lambda = EA h^2/(-N_ref a^2) | linearised buckling factor at P_ref (symmetric mode) |
| F6.lambda_sway | 852443.0035 | - | lambda = EA a^2/(-N_ref h^2) | linearised buckling factor of the sway mode |
| F6.P_cr_linear | 53277.68772 | N | P_cr = 2EA h^3/(a^2 L) | linearised critical load |
| F6.ratio_linear_to_limit | 5.202645361 | - | P_cr,lin / P* | overestimate of linear buckling (-> 3 sqrt3 when shallow) |
| F7.m_crown_sw | 110.0755587 | kg | m = 5*rho*A*L/2 | lumped crown mass, self weight only |
| F7.f_xy_sw | 129.9834146 | Hz | f = sqrt(5/2 k s^2/m)/(2 pi) | horizontal crown frequency (repeated pair), self weight only |
| F7.f_z_sw | 113.6096702 | Hz | f = sqrt(5 k t^2/m)/(2 pi) | vertical crown frequency, self weight only |
| F7.T_z_sw | 0.008802067629 | s | T = 1/f_z | vertical crown period, self weight only |
| F7.omega_xy_sw | 816.7098805 | rad/s | omega = sqrt(Kxx/m) | horizontal circular frequency, self weight only |
| F7.omega_z_sw | 713.8306102 | rad/s | omega = sqrt(Kzz/m) | vertical circular frequency, self weight only |
| F7.phi_sw | 0.09531352924 | 1/sqrt(kg) | phi = 1/sqrt(m) | mass-normalised crown mode amplitude, self weight only |
| F7.m_crown_dead | 720.2827512 | kg | m = 5*rho*A*L/2 + q*A_trib/g | lumped crown mass, self weight + dead |
| F7.f_xy_dead | 50.81382449 | Hz | f = sqrt(5/2 k s^2/m)/(2 pi) | horizontal crown frequency (repeated pair), self weight + dead |
| F7.f_z_dead | 44.41291113 | Hz | f = sqrt(5 k t^2/m)/(2 pi) | vertical crown frequency, self weight + dead |
| F7.T_z_dead | 0.02251597507 | s | T = 1/f_z | vertical crown period, self weight + dead |
| F7.omega_xy_dead | 319.2726754 | rad/s | omega = sqrt(Kxx/m) | horizontal circular frequency, self weight + dead |
| F7.omega_z_dead | 279.0545507 | rad/s | omega = sqrt(Kzz/m) | vertical circular frequency, self weight + dead |
| F7.phi_dead | 0.03726048406 | 1/sqrt(kg) | phi = 1/sqrt(m) | mass-normalised crown mode amplitude, self weight + dead |
| F7.meff_x_sw | 110.0755587 | kg | m_eff = (phi^T M r)^2 = m | effective modal mass of the x mode, self weight only |
| F7.meff_y_sw | 110.0755587 | kg | m_eff = (phi^T M r)^2 = m | effective modal mass of the y mode, self weight only |
| F7.meff_z_sw | 110.0755587 | kg | m_eff = (phi^T M r)^2 = m | effective modal mass of the z mode, self weight only |
| F7.meff_x_dead | 720.2827512 | kg | m_eff = (phi^T M r)^2 = m | effective modal mass of the x mode, self weight + dead |
| F7.meff_y_dead | 720.2827512 | kg | m_eff = (phi^T M r)^2 = m | effective modal mass of the y mode, self weight + dead |
| F7.meff_z_dead | 720.2827512 | kg | m_eff = (phi^T M r)^2 = m | effective modal mass of the z mode, self weight + dead |
| F8.CF_A | 0.5465330578 | - | CF_A = 2 sin(acos(1/sqrt5)/4) | 2V chord factor of strut A (Domerama 0.54653) |
| F8.CF_B | 0.6180339887 | - | CF_B = 2 sin 18 deg = (sqrt5-1)/2 | 2V chord factor of strut B (Domerama 0.61803) |
| F8.L_A_R5 | 2.732665289 | m | L = CF_A*R | strut A length for R = 5 m |
| F8.L_B_R5 | 3.090169944 | m | L = CF_B*R | strut B length for R = 5 m |
| F8.axial_A_deg | 74.14126279 | deg | alpha = acos(L/(2R)) | axial angle of strut A (strut to hub radius) |
| F8.axial_B_deg | 72 | deg | alpha = acos(L/(2R)) | axial angle of strut B (strut to hub radius) |
| F8.bevel_A_deg | 15.85873721 | deg | 90 - alpha = th_A/2 | strut A end angle to the tangent plane |
| F8.bevel_B_deg | 18 | deg | 90 - alpha = th_B/2 | strut B end angle to the tangent plane |
| F8.dihedral_A_deg | 157.5410761 | deg | asin(d1/h)+asin(d2/h), d = sqrt(R^2-rho^2) | interior dihedral angle across A struts (AAB \| AAB) |
| F8.dihedral_B_deg | 161.9708979 | deg | asin(d1/h)+asin(d2/h), d = sqrt(R^2-rho^2) | interior dihedral angle across B struts (AAB \| BBB) |
| F8.bend_A_deg | 22.45892392 | deg | 180 - dihedral | angle between panel normals across A struts |
| F8.bend_B_deg | 18.02910211 | deg | 180 - dihedral | angle between panel normals across B struts |
| F8.face_AAB_apex_deg | 68.86197642 | deg | acos(1 - B^2/(2A^2)) | AAB panel angle at the 5-valent hub |
| F8.face_AAB_base_deg | 55.56901179 | deg | (180 - apex)/2 | AAB panel base angles |
| F8.area_AAB_R1 | 0.139300145 | - | (B/4) sqrt(4A^2 - B^2) | AAB panel area on the unit sphere |
| F8.area_BBB_R1 | 0.1653961346 | - | (sqrt3/4) B^2 | BBB panel area on the unit sphere |
| F8.area_total_R1 | 11.66593139 | - | 60 AAB + 20 BBB | faceted area of the 2V sphere, R = 1 (x R^2 for any R) |
| F8.nodes | 42 | - | 10 v^2 + 2 | 2V Class I full sphere: nodes of the sphere |
| F8.struts | 120 | - | 30 v^2 | 2V Class I full sphere: struts of the sphere |
| F8.panels | 80 | - | 20 v^2 | 2V Class I full sphere: panels of the sphere |
| F8.struts_A | 60 | - | count | 2V Class I full sphere: A struts |
| F8.struts_B | 60 | - | count | 2V Class I full sphere: B struts |
| F8.hubs_valence5 | 12 | - | count | 2V Class I full sphere: 5-valent hubs (type 5A) |
| F8.hubs_valence6 | 30 | - | count | 2V Class I full sphere: 6-valent hubs (type 2A+4B) |
| F8.panels_AAB | 60 | - | count | 2V Class I full sphere: AAB panels |
| F8.panels_BBB | 20 | - | count | 2V Class I full sphere: BBB panels |
| F9.Cp_1 | -1 | - | Cp(theta), theta = acos(-n.d) | pressure coefficient at hand node 1 |
| F9.Fx_1 | 0 | N | F = -qz*G*Cp*A_trib*n | wind force x at hand node 1 |
| F9.Fy_1 | 0 | N | F = -qz*G*Cp*A_trib*n | wind force y at hand node 1 |
| F9.Fz_1 | 17994.81084 | N | F = -qz*G*Cp*A_trib*n | wind force z at hand node 1 |
| F9.Cp_2 | -0.4885501706 | - | Cp(theta), theta = acos(-n.d) | pressure coefficient at hand node 2 |
| F9.Fx_2 | 7863.238501 | N | F = -qz*G*Cp*A_trib*n | wind force x at hand node 2 |
| F9.Fy_2 | 0 | N | F = -qz*G*Cp*A_trib*n | wind force y at hand node 2 |
| F9.Fz_2 | 3931.61925 | N | F = -qz*G*Cp*A_trib*n | wind force z at hand node 2 |
| F9.Cp_3 | -0.8930329524 | - | Cp(theta), theta = acos(-n.d) | pressure coefficient at hand node 3 |
| F9.Fx_3 | 4441.627443 | N | F = -qz*G*Cp*A_trib*n | wind force x at hand node 3 |
| F9.Fy_3 | 13669.92366 | N | F = -qz*G*Cp*A_trib*n | wind force y at hand node 3 |
| F9.Fz_3 | 7186.704168 | N | F = -qz*G*Cp*A_trib*n | wind force z at hand node 3 |
| F9.Cp_4 | -0.1639180562 | - | Cp(theta), theta = acos(-n.d) | pressure coefficient at hand node 4 |
| F9.Fx_4 | -2134.404458 | N | F = -qz*G*Cp*A_trib*n | wind force x at hand node 4 |
| F9.Fy_4 | 1550.735611 | N | F = -qz*G*Cp*A_trib*n | wind force y at hand node 4 |
| F9.Fz_4 | 1319.134501 | N | F = -qz*G*Cp*A_trib*n | wind force z at hand node 4 |
| F9.Cp_5 | -0.1639180562 | - | Cp(theta), theta = acos(-n.d) | pressure coefficient at hand node 5 |
| F9.Fx_5 | -2134.404458 | N | F = -qz*G*Cp*A_trib*n | wind force x at hand node 5 |
| F9.Fy_5 | -1550.735611 | N | F = -qz*G*Cp*A_trib*n | wind force y at hand node 5 |
| F9.Fz_5 | 1319.134501 | N | F = -qz*G*Cp*A_trib*n | wind force z at hand node 5 |
| F9.Cp_6 | -0.8930329524 | - | Cp(theta), theta = acos(-n.d) | pressure coefficient at hand node 6 |
| F9.Fx_6 | 4441.627443 | N | F = -qz*G*Cp*A_trib*n | wind force x at hand node 6 |
| F9.Fy_6 | -13669.92366 | N | F = -qz*G*Cp*A_trib*n | wind force y at hand node 6 |
| F9.Fz_6 | 7186.704168 | N | F = -qz*G*Cp*A_trib*n | wind force z at hand node 6 |
| F9.Cp_7 | -0.6364692711 | - | Cp(theta), theta = acos(-n.d) | pressure coefficient at hand node 7 |
| F9.Fx_7 | 4972.543773 | N | F = -qz*G*Cp*A_trib*n | wind force x at hand node 7 |
| F9.Fy_7 | 3612.764523 | N | F = -qz*G*Cp*A_trib*n | wind force y at hand node 7 |
| F9.Fz_7 | -3073.201062 | N | F = -qz*G*Cp*A_trib*n | wind force z at hand node 7 |
| F9.Cp_8 | -0.7860659049 | - | Cp(theta), theta = acos(-n.d) | pressure coefficient at hand node 8 |
| F9.Fx_8 | -2345.766896 | N | F = -qz*G*Cp*A_trib*n | wind force x at hand node 8 |
| F9.Fy_8 | 7219.528158 | N | F = -qz*G*Cp*A_trib*n | wind force y at hand node 8 |
| F9.Fz_8 | -3795.530567 | N | F = -qz*G*Cp*A_trib*n | wind force z at hand node 8 |
| F9.Cp_9 | 0.2686989765 | - | Cp(theta), theta = acos(-n.d) | pressure coefficient at hand node 9 |
| F9.Fx_9 | 2594.833772 | N | F = -qz*G*Cp*A_trib*n | wind force x at hand node 9 |
| F9.Fy_9 | 0 | N | F = -qz*G*Cp*A_trib*n | wind force y at hand node 9 |
| F9.Fz_9 | 1297.416886 | N | F = -qz*G*Cp*A_trib*n | wind force z at hand node 9 |
| F9.Cp_10 | -0.7860659049 | - | Cp(theta), theta = acos(-n.d) | pressure coefficient at hand node 10 |
| F9.Fx_10 | -2345.766896 | N | F = -qz*G*Cp*A_trib*n | wind force x at hand node 10 |
| F9.Fy_10 | -7219.528158 | N | F = -qz*G*Cp*A_trib*n | wind force y at hand node 10 |
| F9.Fz_10 | -3795.530567 | N | F = -qz*G*Cp*A_trib*n | wind force z at hand node 10 |
| F9.Cp_11 | -0.6364692711 | - | Cp(theta), theta = acos(-n.d) | pressure coefficient at hand node 11 |
| F9.Fx_11 | 4972.543773 | N | F = -qz*G*Cp*A_trib*n | wind force x at hand node 11 |
| F9.Fy_11 | -3612.764523 | N | F = -qz*G*Cp*A_trib*n | wind force y at hand node 11 |
| F9.Fz_11 | -3073.201062 | N | F = -qz*G*Cp*A_trib*n | wind force z at hand node 11 |
| F9.qz | 1061.332875 | Pa | qz = 0.613*Kz*Kzt*Kd*V^2 | velocity pressure, default wind data |
| F9.qG | 902.1329437 | Pa | qz*G | velocity pressure times gust factor |
| F9.Rx | 20326.072 | N | sum of the nodal wind forces | wind resultant x |
| F9.Ry | 0 | N | sum of the nodal wind forces | wind resultant y |
| F9.Rz | 26498.06105 | N | sum of the nodal wind forces | wind resultant z |
| F10.phiPn_inel | 239578.0961 | N | phi*Fcr*Q*Ag (E3 + E7) | design compression strength, CHS 88.9 x 4.0, L = 2 m |
| F10.Fcr_inel | 249509476.9 | Pa | E3-2 | critical stress, CHS 88.9 x 4.0, L = 2 m |
| F10.Q_inel | 1 | - | 1, or 0.038E/(Fy D/t) + 2/3 (E7-7) | effective area ratio Ae/Ag, CHS 88.9 x 4.0, L = 2 m |
| F10.KLr_inel | 66.5557865 | - | K*L/r | slenderness, CHS 88.9 x 4.0, L = 2 m |
| F10.phiPn_b42 | 223995.9127 | N | phi*Fcr*Q*Ag (E3 + E7) | design compression strength, CHS 88.9 x 4.0 with t_des = 0.93 t, L = 2 m |
| F10.Fcr_b42 | 250015616.3 | Pa | E3-2 | critical stress, CHS 88.9 x 4.0 with t_des = 0.93 t, L = 2 m |
| F10.Q_b42 | 1 | - | 1, or 0.038E/(Fy D/t) + 2/3 (E7-7) | effective area ratio Ae/Ag, CHS 88.9 x 4.0 with t_des = 0.93 t, L = 2 m |
| F10.KLr_b42 | 66.34735148 | - | K*L/r | slenderness, CHS 88.9 x 4.0 with t_des = 0.93 t, L = 2 m |
| F10.phiPn_e7 | 1112934.907 | N | phi*Fcr*Q*Ag (E3 + E7) | design compression strength, D 323.9 x 4.0, L = 3 m |
| F10.Fcr_e7 | 327695234.8 | Pa | E3-2 | critical stress, D 323.9 x 4.0, L = 3 m |
| F10.Q_e7 | 0.9387134158 | - | 1, or 0.038E/(Fy D/t) + 2/3 (E7-7) | effective area ratio Ae/Ag, D 323.9 x 4.0, L = 3 m |
| F10.KLr_e7 | 26.52271999 | - | K*L/r | slenderness, D 323.9 x 4.0, L = 3 m |
| F10.phiPn_scope | 0 | N | phi*Fcr*Q*Ag (E3 + E7) | design compression strength, D 323.9 x 1.2, L = 3 m |
| F10.phiPn_t_b42 | 309095.0519 | N | phi_t*Fy*Ag(t_des) (D2-1) | design tension yield strength, CHS 88.9 x 4.0 with t_des = 0.93 t |

## References

- Kassimali, A. Matrix Analysis of Structures, 2nd ed. Cengage, 2012 (direct stiffness method for plane and space trusses, member forces from end displacements, support settlements and temperature changes by equivalent joint loads).
- Hibbeler, R. C. Structural Analysis, 10th ed. Pearson, 2018 (force method of analysis for trusses, principle of virtual work, unit-load method, temperature and fabrication errors, support settlement).
- Kenner, H. Geodesic Math and How to Use It. University of California Press, 1976 (Class I breakdown).
- Domerama, Chord factor charts, https://domerama.com/calculators/chord-charts, read 16 November 2026 (method 1 icosahedral chord factors to five decimals; 2V 0.54653 and 0.61803; closed forms 2 sin(15.8587 deg) and 2 sin 18 deg).
- ASCE/SEI 7-16 Minimum Design Loads, Eq. 26.10-1 (velocity pressure qz = 0.613 Kz Kzt Kd V^2 in SI). The 7-point Cp table of F9 is the program's illustrative simplification, not values read from Fig. 27.3-2; F9 checks the implementation, not the table.
- AISC 360-16 Specification for Structural Steel Buildings, Sec. B4.2 (design wall thickness 0.93 t for ERW HSS), Table B4.1a, Eqs. E3-2 to E3-4 and E7-7 (F10).
- Crisfield, M. A. Non-linear Finite Element Analysis of Solids and Structures, Vol. 1. Wiley, 1991 (shallow two-bar truss, limit point, snap-through).
- Bazant, Z. P., and Cedolin, L. Stability of Structures. Oxford University Press, 1991 (von Mises truss, linearised versus exact limit load).
- Chopra, A. K. Dynamics of Structures, 5th ed. Pearson, 2017 (free vibration, lumped mass, mass-normalised modes, effective modal mass).
- AISC 360-22 Specification for Structural Steel Buildings, Chapter E (compression strength used for the utilisation notes), and the AISC Steel Construction Manual, 16th ed. (coefficient of thermal expansion of structural steel, 6.5e-6 /degF = 11.7e-6 /degC).
