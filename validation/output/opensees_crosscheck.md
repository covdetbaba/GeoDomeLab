# GeoDomeLab vs OpenSees: independent cross-check

Generated 2026-11-16 by `validation/opensees_crosscheck.py`.  
GeoDomeLab 2.0.0 (Python 3.13.9); OpenSees 3.8.0 via OpenSeesPy, run subprocess (separate OpenSees interpreter).  Total time 13.4 s.

**Result: 835 of 835 comparisons pass.**

## What is compared

Each structure is solved by GeoDomeLab and by OpenSees, an independent finite element framework (McKenna, Scott & Fenves 2010; OpenSeesPy, Zhu, McKenna & Scott 2018).  The programs share only the input data written by `opensees_bridge.to_opensees_spec`: node coordinates, connectivity, E and A, supports, the factored nodal load vector of the combination, the thermal strain T*alpha*dT, the settlements T*delta and (modal) the masses.  OpenSees builds and solves everything else itself: `Truss` and `corotTruss` elements with an `Elastic` material, `InitStrainMaterial(-alpha dT)` for temperature, `sp` constraints in the load pattern for settlement, `zeroLength` springs to fixed ground nodes for elastic supports, element mass `-rho` (lumped, or consistent with `-cMass 1`) plus nodal `mass`, `eigen` (ARPACK or LAPACK), and Newton iteration with `DisplacementControl` for the nonlinear path.

Load generation (tributary areas, the wind pressure coefficients) and the non-structural masses are therefore *inputs* here; they are validated by hand calculation in `validation/run_validation.py` and `validation/run_validation_ext.py`.

Differences are normwise relative: for displacements max|u_GDL - u_OS| / max|u_GDL| (likewise N, with the fully restrained thermal force EA*alpha*dT as a floor of the scale, and R, with max|F| as a floor); natural frequencies pointwise |f_GDL - f_OS| / f_GDL; mode shapes by 1 - MAC, and inside clusters of exactly repeated frequencies (the symmetric pairs of a dome, whose basis is arbitrary) by 1 - subspace MAC; the nonlinear load factor at matched states by max|lambda_GDL - lambda_OS| / max|lambda| over the path.  Targets (spec 2.E): static 1e-08, frequencies 1e-06, nonlinear lambda 1e-05.

## Summary (worst case of every quantity)

| section | cases | comparisons | pass | quantity | worst difference | target | worst case |
|---|---:|---:|---:|---|---:|---:|---|
| Formulation probes | 3 | 6 | 6 | N_OS vs engineering strain EA(l-L)/L | 1.46e-15 | <= 1e-09 | P1 - |
|  |  |  |  | u2 vs closed form | 1.16e-16 | <= 1e-09 | P2 dT = 30 C |
|  |  |  |  | N1, N2 vs closed form | 1.55e-16 | <= 1e-09 | P2 dT = 30 C |
| Linear static analysis | 34 | 680 | 680 | displacements u | 2.19e-13 | <= 1e-08 | S-8V-elastic LC3  1.2D + 1.0W + 0.5S |
|  |  |  |  | member forces N | 1.72e-13 | <= 1e-08 | S-8V-elastic LC3  1.2D + 1.0W + 0.5S |
|  |  |  |  | reactions R (incl. springs) | 6.14e-14 | <= 1e-08 | S-8V-roller3 LC1  1.4D |
|  |  |  |  | OpenSees residual at free DOFs / \|F\| | 4.09e-14 | <= 1e-08 | S-7V-elastic SLS  1.0D + 1.0S |
|  |  |  |  | spring forces | 4.59e-14 | <= 1e-08 | S-8V-elastic LC3  1.2D + 1.0W + 0.5S |
| Modal analysis | 35 | 97 | 97 | natural frequencies f | 2.72e-13 | <= 1e-06 | M-8V-roller3 lumped mass |
|  |  |  |  | 1 - MAC (single modes) | 3.33e-15 | <= 1e-06 | M-8V-elastic lumped mass |
|  |  |  |  | 1 - subspace MAC (repeated pairs) | 4.00e-15 | <= 1e-06 | M-8V-elastic lumped mass |
| Geometrically nonlinear analysis | 8 | 47 | 47 | lambda at matched states / max\|lambda\| | 4.14e-11 | <= 1e-05 | N-case6V-LC2 LC2  1.2D + 1.6S |
|  |  |  |  | lambda pointwise (\|lambda\| > 5 % of max) | 4.38e-11 | <= 1e-05 | N-case6V-LC2 LC2  1.2D + 1.6S |
|  |  |  |  | displacements at matched states | 3.86e-11 | <= 1e-05 | N-case6V-LC2 LC2  1.2D + 1.6S |
|  |  |  |  | lambda at the refined first limit point | 3.16e-11 | <= 1e-05 | N-von-mises-deep P  1.0L  apex reference load |
|  |  |  |  | GDL equilibrium residual of OpenSees states | 3.37e-11 | <= 1e-08 | N-star-dome P  1.0L  crown reference load |
|  |  |  |  | GDL residual of OpenSees ArcLength points | 6.38e-12 | <= 1e-08 | N-star-dome P  1.0L  crown reference load |
|  |  |  |  | OpenSees P(w) vs exact closed form | 2.16e-13 | <= 1e-05 | N-von-mises-deep P  1.0L  apex reference load |
|  |  |  |  | OpenSees limit load vs closed form P* | 1.05e-12 | <= 1e-05 | N-von-mises-deep P  1.0L  apex reference load |
| Exported standalone OpenSeesPy scripts | 3 | 5 | 5 | exported static script: u | 6.32e-15 | <= 1e-08 | X1 T1  1.2D + 1.6S + 1.0T |
|  |  |  |  | exported static script: N | 5.79e-15 | <= 1e-08 | X1 T1  1.2D + 1.6S + 1.0T |
|  |  |  |  | exported static script: R | 2.00e-15 | <= 1e-08 | X1 T1  1.2D + 1.6S + 1.0T |
|  |  |  |  | exported modal script: f (10 modes) | 2.76e-15 | <= 1e-06 | X2 - |
|  |  |  |  | exported nonlinear script: GDL residual of its states | 9.86e-13 | <= 1e-08 | X3 LC2 |

## P. Formulation probes

**Strain measure of `corotTruss` (decided numerically).**  One bar of length L = 2 m was rotated by 44.6 degrees (prescribed transverse end displacement 0.8 L) and pulled axially (P = 0.1 EA) to l/L = 1.140328.  OpenSees returned N = 2.806554437e+07 N.  Predictions: engineering strain EA(l-L)/L = 2.806554437e+07 N (difference 1.5e-15), Green-Lagrange EA(l^2-L^2)/(2L^2) = 3.003473132e+07 N (6.6e-02), logarithmic EA ln(l/L) = 2.626313928e+07 N (6.9e-02).  OpenSees 3.8.0 `corotTruss` therefore uses **engineering strain on the original length**, the same measure as GeoDomeLab's corotational bar (and as the OpenSees source, CorotTruss.cpp).  The nonlinear paths can be compared at any strain; the deep von Mises case below reaches 33 % strain.

**Sign of `InitStrainMaterial`.**  Bar 1 (heated, free strain eps0 = 3.5100e-04) and bar 2 (twice the area) in series between two walls.  Closed form u2 = eps0 L/(1 + r), N = -EA eps0 r/(1 + r) with r = 2: u2 = 2.340000000e-04 m, N = -46800.000000 N.  OpenSees with epsInit = -eps0: u2 = 2.340000000e-04 m, N1 = -46800.000000 N, N2 = -46800.000000 N.  So `InitStrainMaterial(tag, E, -alpha dT)` gives N = EA(eps - alpha dT): heating compresses a restrained bar.

**Sensitivity control.**  The same 4V comparison with E raised by 1 ppm in the OpenSees input only must show a displacement difference of about 1e-6.  It does (row P3), so the round-off-level agreement elsewhere is a real agreement of two independent computations, not a program compared with itself.

| case | structure | quantity | difference | target | pass | note |
|---|---|---|---:|---:|---|---|
| P1 | single corotTruss bar | N_OS vs engineering strain EA(l-L)/L | 1.46e-15 | <= 1e-09 | yes | bar rotated 44.6 deg, stretched to l/L = 1.1403; N_OS = 2.806554e+07 N |
| P1 | single corotTruss bar | N_OS vs Green-Lagrange EA(l^2-L^2)/(2L^2) (must differ) | 6.56e-02 | >= 1e-03 | yes | different strain measure |
| P1 | single corotTruss bar | N_OS vs logarithmic EA ln(l/L) (must differ) | 6.86e-02 | >= 1e-03 | yes | different strain measure |
| P2 | heated bar + bar between walls | u2 vs closed form | 1.16e-16 | <= 1e-09 | yes | u2 = 2.340000e-04 m (exact 2.340000e-04), N = -4.680000e+04 N (exact -4.680000e+04) |
| P2 | heated bar + bar between walls | N1, N2 vs closed form | 1.55e-16 | <= 1e-09 | yes | InitStrainMaterial(-alpha dT) gives N = EA(eps - alpha dT) |
| P3 | 4V R=12 m f=0.5 | du with E x (1 + 1e-6) in OpenSees only (must be ~1e-6) | 1.00e-06 | in [5.0e-07, 2.0e-06] | yes | the measure resolves a 1 ppm input change, 100x below the static target |

## S. Linear static analysis

Geodesic domes: radius 12 m, dome fraction 0.5 requested (the generator cuts at the nearest node ring, so odd frequencies come out lower; the structure column shows the actual fraction f), default sections (CHS 88.9 x 4.0 unless stated), default loads and the ASCE 7 combinations LC1 1.4D, LC2 1.2D + 1.6S, LC3 1.2D + 1.0W + 0.5S, LC4 0.9D + 1.0W, SLS 1.0D + 1.0S.  Springs (elastic supports) 50 MN/m in x, y and z.  T cases: dT = +30 C on every strut and two base nodes settling (10 mm down, and 3 / -2 / -5 mm in x / y / z; vertical only on the roller supports), combinations 1.2D + 1.6S + 1.0T, 1.0T and 1.2D + 1.0W + 0.5S + 1.2T.  Custom: rigid, spring (20 and 30 MN/m; 100 and 50 MN/m) and free directions mixed at the base, one settling rigid support, one spring support whose ground end moves, and dT = -25 C.  The two extra fractions use unbalanced snow and wind from 30 degrees.

| case | structure | supports | combination | displacements u | member forces N | reactions R (incl. springs) | OpenSees residual at free DOFs / \|F\| | spring forces | pass | note |
|---|---|---|---|---:|---:|---:|---:|---:|---|---|
| S-1V-pinned | 1V R=12 m f=0.276 | pinned | LC1  1.4D | 0 | 1.84e-16 | 7.01e-17 | 2.10e-16 |  | yes | 3 free DOF, \|u\|max = 2.22 mm, \|N\|max = 19.74 kN, OS system FullGeneral |
| S-1V-pinned | 1V R=12 m f=0.276 | pinned | LC2  1.2D + 1.6S | 0 | 1.74e-16 | 1.32e-16 | 1.66e-16 |  | yes | 3 free DOF, \|u\|max = 9.404 mm, \|N\|max = 83.61 kN, OS system FullGeneral |
| S-1V-pinned | 1V R=12 m f=0.276 | pinned | LC3  1.2D + 1.0W + 0.5S | 0 | 1.36e-16 | 1.22e-17 | 3.35e-17 |  | yes | 3 free DOF, \|u\|max = 0.1882 mm, \|N\|max = 1.673 kN, OS system FullGeneral |
| S-1V-pinned | 1V R=12 m f=0.276 | pinned | LC4  0.9D + 1.0W | 0 | 4.08e-16 | 1.55e-16 | 2.59e-16 |  | yes | 3 free DOF, \|u\|max = 3.008 mm, \|N\|max = 26.74 kN, OS system FullGeneral |
| S-1V-pinned | 1V R=12 m f=0.276 | pinned | SLS  1.0D + 1.0S | 0 | 1.30e-16 | 9.92e-17 | 1.74e-16 |  | yes | 3 free DOF, \|u\|max = 6.274 mm, \|N\|max = 55.78 kN, OS system FullGeneral |
| S-1V-roller3 | 1V R=12 m f=0.276 | 3 rollers | LC1  1.4D | 1.54e-16 | 1.28e-16 | 1.05e-16 | 2.10e-16 |  | yes | 7 free DOF, \|u\|max = 2.821 mm, \|N\|max = 28.38 kN, OS system FullGeneral |
| S-1V-roller3 | 1V R=12 m f=0.276 | 3 rollers | LC2  1.2D + 1.6S | 1.45e-16 | 1.82e-16 | 1.66e-16 | 2.65e-16 |  | yes | 7 free DOF, \|u\|max = 11.95 mm, \|N\|max = 120.2 kN, OS system FullGeneral |
| S-1V-roller3 | 1V R=12 m f=0.276 | 3 rollers | LC3  1.2D + 1.0W + 0.5S | 2.60e-16 | 4.18e-16 | 2.33e-16 | 1.86e-16 |  | yes | 7 free DOF, \|u\|max = 1.666 mm, \|N\|max = 17.41 kN, OS system FullGeneral |
| S-1V-roller3 | 1V R=12 m f=0.276 | 3 rollers | LC4  0.9D + 1.0W | 4.11e-17 | 1.17e-16 | 1.04e-16 | 2.33e-16 |  | yes | 7 free DOF, \|u\|max = 2.637 mm, \|N\|max = 31.06 kN, OS system FullGeneral |
| S-1V-roller3 | 1V R=12 m f=0.276 | 3 rollers | SLS  1.0D + 1.0S | 2.18e-16 | 5.44e-16 | 9.92e-17 | 9.92e-17 |  | yes | 7 free DOF, \|u\|max = 7.971 mm, \|N\|max = 80.19 kN, OS system FullGeneral |
| S-1V-elastic | 1V R=12 m f=0.276 | elastic springs | LC1  1.4D | 1.34e-16 | 3.69e-16 | 2.80e-16 | 5.26e-16 | 2.80e-16 | yes | 18 free DOF, \|u\|max = 3.227 mm, \|N\|max = 19.74 kN, OS system FullGeneral |
| S-1V-elastic | 1V R=12 m f=0.276 | elastic springs | LC2  1.2D + 1.6S | 1.27e-16 | 3.48e-16 | 1.32e-16 | 4.30e-16 | 1.32e-16 | yes | 18 free DOF, \|u\|max = 13.62 mm, \|N\|max = 83.61 kN, OS system FullGeneral |
| S-1V-elastic | 1V R=12 m f=0.276 | elastic springs | LC3  1.2D + 1.0W + 0.5S | 2.70e-16 | 4.57e-16 | 2.93e-16 | 2.93e-16 | 2.93e-16 | yes | 18 free DOF, \|u\|max = 0.8019 mm, \|N\|max = 6.472 kN, OS system FullGeneral |
| S-1V-elastic | 1V R=12 m f=0.276 | elastic springs | LC4  0.9D + 1.0W | 2.64e-16 | 3.87e-16 | 1.29e-16 | 3.11e-16 | 1.29e-16 | yes | 18 free DOF, \|u\|max = 3.282 mm, \|N\|max = 28.22 kN, OS system FullGeneral |
| S-1V-elastic | 1V R=12 m f=0.276 | elastic springs | SLS  1.0D + 1.0S | 5.73e-16 | 1.04e-15 | 2.23e-16 | 6.20e-16 | 2.23e-16 | yes | 18 free DOF, \|u\|max = 9.09 mm, \|N\|max = 55.78 kN, OS system FullGeneral |
| S-2V-pinned | 2V R=12 m f=0.5 | pinned | LC1  1.4D | 1.21e-15 | 1.58e-15 | 5.39e-16 | 7.18e-16 |  | yes | 48 free DOF, \|u\|max = 2.877 mm, \|N\|max = 19.63 kN, OS system FullGeneral |
| S-2V-pinned | 2V R=12 m f=0.5 | pinned | LC2  1.2D + 1.6S | 1.10e-15 | 1.69e-15 | 7.84e-16 | 1.00e-15 |  | yes | 48 free DOF, \|u\|max = 12.65 mm, \|N\|max = 81.83 kN, OS system FullGeneral |
| S-2V-pinned | 2V R=12 m f=0.5 | pinned | LC3  1.2D + 1.0W + 0.5S | 1.59e-15 | 1.60e-15 | 7.32e-16 | 1.19e-15 |  | yes | 48 free DOF, \|u\|max = 2.999 mm, \|N\|max = 42.1 kN, OS system FullGeneral |
| S-2V-pinned | 2V R=12 m f=0.5 | pinned | LC4  0.9D + 1.0W | 1.66e-15 | 1.68e-15 | 8.44e-16 | 1.15e-15 |  | yes | 48 free DOF, \|u\|max = 2.358 mm, \|N\|max = 26.03 kN, OS system FullGeneral |
| S-2V-pinned | 2V R=12 m f=0.5 | pinned | SLS  1.0D + 1.0S | 7.21e-16 | 1.53e-15 | 7.30e-16 | 8.44e-16 |  | yes | 48 free DOF, \|u\|max = 8.419 mm, \|N\|max = 54.57 kN, OS system FullGeneral |
| S-2V-roller3 | 2V R=12 m f=0.5 | 3 rollers | LC1  1.4D | 1.10e-15 | 1.67e-15 | 5.87e-16 | 8.29e-16 |  | yes | 62 free DOF, \|u\|max = 3.161 mm, \|N\|max = 25.67 kN, OS system FullGeneral |
| S-2V-roller3 | 2V R=12 m f=0.5 | 3 rollers | LC2  1.2D + 1.6S | 8.94e-16 | 1.86e-15 | 7.99e-16 | 1.22e-15 |  | yes | 62 free DOF, \|u\|max = 13.58 mm, \|N\|max = 85.96 kN, OS system FullGeneral |
| S-2V-roller3 | 2V R=12 m f=0.5 | 3 rollers | LC3  1.2D + 1.0W + 0.5S | 1.38e-15 | 9.44e-16 | 8.43e-16 | 7.52e-16 |  | yes | 62 free DOF, \|u\|max = 3.445 mm, \|N\|max = 46.25 kN, OS system FullGeneral |
| S-2V-roller3 | 2V R=12 m f=0.5 | 3 rollers | LC4  0.9D + 1.0W | 1.10e-15 | 2.26e-15 | 1.00e-15 | 1.27e-15 |  | yes | 62 free DOF, \|u\|max = 2.759 mm, \|N\|max = 27.9 kN, OS system FullGeneral |
| S-2V-roller3 | 2V R=12 m f=0.5 | 3 rollers | SLS  1.0D + 1.0S | 8.62e-16 | 1.87e-15 | 5.82e-16 | 1.15e-15 |  | yes | 62 free DOF, \|u\|max = 9.053 mm, \|N\|max = 58.27 kN, OS system FullGeneral |
| S-2V-elastic | 2V R=12 m f=0.5 | elastic springs | LC1  1.4D | 1.04e-15 | 3.06e-15 | 8.98e-16 | 1.75e-15 | 8.98e-16 | yes | 78 free DOF, \|u\|max = 3.747 mm, \|N\|max = 19.63 kN, OS system FullGeneral |
| S-2V-elastic | 2V R=12 m f=0.5 | elastic springs | LC2  1.2D + 1.6S | 1.00e-15 | 2.58e-15 | 7.53e-16 | 2.35e-15 | 7.53e-16 | yes | 78 free DOF, \|u\|max = 15.17 mm, \|N\|max = 81.83 kN, OS system FullGeneral |
| S-2V-elastic | 2V R=12 m f=0.5 | elastic springs | LC3  1.2D + 1.0W + 0.5S | 1.71e-15 | 1.94e-15 | 7.31e-16 | 1.23e-15 | 7.31e-16 | yes | 78 free DOF, \|u\|max = 3.794 mm, \|N\|max = 42.11 kN, OS system FullGeneral |
| S-2V-elastic | 2V R=12 m f=0.5 | elastic springs | LC4  0.9D + 1.0W | 1.65e-15 | 1.54e-15 | 8.81e-16 | 1.93e-15 | 8.81e-16 | yes | 78 free DOF, \|u\|max = 2.636 mm, \|N\|max = 26.05 kN, OS system FullGeneral |
| S-2V-elastic | 2V R=12 m f=0.5 | elastic springs | SLS  1.0D + 1.0S | 1.03e-15 | 2.70e-15 | 1.10e-15 | 1.92e-15 | 1.10e-15 | yes | 78 free DOF, \|u\|max = 10.15 mm, \|N\|max = 54.57 kN, OS system FullGeneral |
| S-3V-pinned | 3V R=12 m f=0.414 | pinned | LC1  1.4D | 1.92e-15 | 2.74e-15 | 7.44e-16 | 9.29e-16 |  | yes | 93 free DOF, \|u\|max = 1.803 mm, \|N\|max = 14.58 kN, OS system FullGeneral |
| S-3V-pinned | 3V R=12 m f=0.414 | pinned | LC2  1.2D + 1.6S | 1.39e-15 | 3.08e-15 | 1.21e-15 | 1.45e-15 |  | yes | 93 free DOF, \|u\|max = 8.084 mm, \|N\|max = 59.01 kN, OS system FullGeneral |
| S-3V-pinned | 3V R=12 m f=0.414 | pinned | LC3  1.2D + 1.0W + 0.5S | 3.35e-15 | 2.61e-15 | 1.04e-15 | 2.87e-15 |  | yes | 93 free DOF, \|u\|max = 1.942 mm, \|N\|max = 29.99 kN, OS system FullGeneral |
| S-3V-pinned | 3V R=12 m f=0.414 | pinned | LC4  0.9D + 1.0W | 2.99e-15 | 3.51e-15 | 1.79e-15 | 3.67e-15 |  | yes | 93 free DOF, \|u\|max = 1.269 mm, \|N\|max = 16.57 kN, OS system FullGeneral |
| S-3V-pinned | 3V R=12 m f=0.414 | pinned | SLS  1.0D + 1.0S | 1.61e-15 | 2.47e-15 | 9.71e-16 | 1.06e-15 |  | yes | 93 free DOF, \|u\|max = 5.374 mm, \|N\|max = 39.48 kN, OS system FullGeneral |
| S-3V-roller3 | 3V R=12 m f=0.414 | 3 rollers | LC1  1.4D | 3.37e-15 | 3.54e-15 | 2.90e-15 | 1.47e-15 |  | yes | 117 free DOF, \|u\|max = 2.173 mm, \|N\|max = 23.54 kN, OS system FullGeneral |
| S-3V-roller3 | 3V R=12 m f=0.414 | 3 rollers | LC2  1.2D + 1.6S | 2.93e-15 | 3.50e-15 | 3.34e-15 | 1.76e-15 |  | yes | 117 free DOF, \|u\|max = 9.335 mm, \|N\|max = 80.58 kN, OS system FullGeneral |
| S-3V-roller3 | 3V R=12 m f=0.414 | 3 rollers | LC3  1.2D + 1.0W + 0.5S | 5.48e-15 | 3.70e-15 | 2.84e-15 | 1.42e-15 |  | yes | 117 free DOF, \|u\|max = 2.372 mm, \|N\|max = 32.98 kN, OS system FullGeneral |
| S-3V-roller3 | 3V R=12 m f=0.414 | 3 rollers | LC4  0.9D + 1.0W | 5.00e-15 | 4.94e-15 | 1.46e-15 | 1.90e-15 |  | yes | 117 free DOF, \|u\|max = 1.539 mm, \|N\|max = 19.16 kN, OS system FullGeneral |
| S-3V-roller3 | 3V R=12 m f=0.414 | 3 rollers | SLS  1.0D + 1.0S | 3.20e-15 | 4.32e-15 | 3.01e-15 | 1.81e-15 |  | yes | 117 free DOF, \|u\|max = 6.222 mm, \|N\|max = 54.56 kN, OS system FullGeneral |
| S-3V-elastic | 3V R=12 m f=0.414 | elastic springs | LC1  1.4D | 1.84e-15 | 5.47e-15 | 1.32e-15 | 3.27e-15 | 1.32e-15 | yes | 138 free DOF, \|u\|max = 2.352 mm, \|N\|max = 14.63 kN, OS system FullGeneral |
| S-3V-elastic | 3V R=12 m f=0.414 | elastic springs | LC2  1.2D + 1.6S | 3.01e-15 | 5.29e-15 | 2.31e-15 | 2.99e-15 | 2.31e-15 | yes | 138 free DOF, \|u\|max = 9.796 mm, \|N\|max = 59.14 kN, OS system FullGeneral |
| S-3V-elastic | 3V R=12 m f=0.414 | elastic springs | LC3  1.2D + 1.0W + 0.5S | 2.14e-15 | 3.28e-15 | 1.50e-15 | 2.72e-15 | 1.50e-15 | yes | 138 free DOF, \|u\|max = 2.429 mm, \|N\|max = 29.82 kN, OS system FullGeneral |
| S-3V-elastic | 3V R=12 m f=0.414 | elastic springs | LC4  0.9D + 1.0W | 3.67e-15 | 4.84e-15 | 2.30e-15 | 2.75e-15 | 2.30e-15 | yes | 138 free DOF, \|u\|max = 1.358 mm, \|N\|max = 16.53 kN, OS system FullGeneral |
| S-3V-elastic | 3V R=12 m f=0.414 | elastic springs | SLS  1.0D + 1.0S | 1.34e-15 | 5.70e-15 | 1.41e-15 | 3.53e-15 | 1.41e-15 | yes | 138 free DOF, \|u\|max = 6.542 mm, \|N\|max = 39.57 kN, OS system FullGeneral |
| S-4V-pinned | 4V R=12 m f=0.5 | pinned | LC1  1.4D | 5.13e-15 | 9.49e-15 | 1.94e-15 | 3.40e-15 |  | yes | 213 free DOF, \|u\|max = 1.922 mm, \|N\|max = 16.19 kN, OS system FullGeneral |
| S-4V-pinned | 4V R=12 m f=0.5 | pinned | LC2  1.2D + 1.6S | 5.15e-15 | 1.09e-14 | 2.68e-15 | 4.72e-15 |  | yes | 213 free DOF, \|u\|max = 7.792 mm, \|N\|max = 57.31 kN, OS system FullGeneral |
| S-4V-pinned | 4V R=12 m f=0.5 | pinned | LC3  1.2D + 1.0W + 0.5S | 9.45e-15 | 5.86e-15 | 2.44e-15 | 5.22e-15 |  | yes | 213 free DOF, \|u\|max = 2.248 mm, \|N\|max = 34.14 kN, OS system FullGeneral |
| S-4V-pinned | 4V R=12 m f=0.5 | pinned | LC4  0.9D + 1.0W | 8.33e-15 | 8.80e-15 | 3.36e-15 | 1.15e-14 |  | yes | 213 free DOF, \|u\|max = 1.278 mm, \|N\|max = 17.91 kN, OS system FullGeneral |
| S-4V-pinned | 4V R=12 m f=0.5 | pinned | SLS  1.0D + 1.0S | 5.03e-15 | 1.18e-14 | 2.43e-15 | 3.65e-15 |  | yes | 213 free DOF, \|u\|max = 5.213 mm, \|N\|max = 38.71 kN, OS system FullGeneral |
| S-4V-roller3 | 4V R=12 m f=0.5 | 3 rollers | LC1  1.4D | 1.22e-14 | 1.13e-14 | 1.08e-14 | 3.87e-15 |  | yes | 247 free DOF, \|u\|max = 2.031 mm, \|N\|max = 18.23 kN, OS system FullGeneral |
| S-4V-roller3 | 4V R=12 m f=0.5 | 3 rollers | LC2  1.2D + 1.6S | 8.91e-15 | 1.25e-14 | 9.80e-15 | 4.32e-15 |  | yes | 247 free DOF, \|u\|max = 8.106 mm, \|N\|max = 62.24 kN, OS system FullGeneral |
| S-4V-roller3 | 4V R=12 m f=0.5 | 3 rollers | LC3  1.2D + 1.0W + 0.5S | 1.72e-14 | 9.39e-15 | 9.88e-15 | 5.18e-15 |  | yes | 247 free DOF, \|u\|max = 2.421 mm, \|N\|max = 37.69 kN, OS system FullGeneral |
| S-4V-roller3 | 4V R=12 m f=0.5 | 3 rollers | LC4  0.9D + 1.0W | 1.07e-14 | 9.49e-15 | 5.23e-15 | 3.71e-15 |  | yes | 247 free DOF, \|u\|max = 2.087 mm, \|N\|max = 24.93 kN, OS system FullGeneral |
| S-4V-roller3 | 4V R=12 m f=0.5 | 3 rollers | SLS  1.0D + 1.0S | 6.87e-15 | 8.46e-15 | 6.91e-15 | 3.59e-15 |  | yes | 247 free DOF, \|u\|max = 5.429 mm, \|N\|max = 42.15 kN, OS system FullGeneral |
| S-4V-elastic | 4V R=12 m f=0.5 | elastic springs | LC1  1.4D | 7.56e-15 | 1.19e-14 | 2.00e-15 | 5.26e-15 | 2.00e-15 | yes | 273 free DOF, \|u\|max = 2.41 mm, \|N\|max = 16.32 kN, OS system FullGeneral |
| S-4V-elastic | 4V R=12 m f=0.5 | elastic springs | LC2  1.2D + 1.6S | 5.92e-15 | 1.04e-14 | 1.82e-15 | 6.73e-15 | 1.82e-15 | yes | 273 free DOF, \|u\|max = 9.09 mm, \|N\|max = 57.74 kN, OS system FullGeneral |
| S-4V-elastic | 4V R=12 m f=0.5 | elastic springs | LC3  1.2D + 1.0W + 0.5S | 8.44e-15 | 8.54e-15 | 3.01e-15 | 6.70e-15 | 3.01e-15 | yes | 273 free DOF, \|u\|max = 2.671 mm, \|N\|max = 34.07 kN, OS system FullGeneral |
| S-4V-elastic | 4V R=12 m f=0.5 | elastic springs | LC4  0.9D + 1.0W | 9.30e-15 | 1.07e-14 | 6.45e-15 | 1.19e-14 | 6.45e-15 | yes | 273 free DOF, \|u\|max = 1.469 mm, \|N\|max = 17.83 kN, OS system FullGeneral |
| S-4V-elastic | 4V R=12 m f=0.5 | elastic springs | SLS  1.0D + 1.0S | 5.96e-15 | 1.04e-14 | 1.98e-15 | 6.05e-15 | 1.98e-15 | yes | 273 free DOF, \|u\|max = 6.111 mm, \|N\|max = 39 kN, OS system FullGeneral |
| S-5V-pinned | 5V R=12 m f=0.448 | pinned | LC1  1.4D | 1.42e-14 | 2.16e-14 | 4.17e-15 | 7.22e-15 |  | yes | 303 free DOF, \|u\|max = 1.469 mm, \|N\|max = 12.02 kN, OS system FullGeneral |
| S-5V-pinned | 5V R=12 m f=0.448 | pinned | LC2  1.2D + 1.6S | 9.67e-15 | 2.22e-14 | 7.59e-15 | 6.09e-15 |  | yes | 303 free DOF, \|u\|max = 5.917 mm, \|N\|max = 43.68 kN, OS system FullGeneral |
| S-5V-pinned | 5V R=12 m f=0.448 | pinned | LC3  1.2D + 1.0W + 0.5S | 3.00e-14 | 1.83e-14 | 7.52e-15 | 8.15e-15 |  | yes | 303 free DOF, \|u\|max = 1.721 mm, \|N\|max = 26.39 kN, OS system FullGeneral |
| S-5V-pinned | 5V R=12 m f=0.448 | pinned | LC4  0.9D + 1.0W | 2.27e-14 | 1.68e-14 | 1.15e-14 | 1.00e-14 |  | yes | 303 free DOF, \|u\|max = 1.04 mm, \|N\|max = 13.72 kN, OS system FullGeneral |
| S-5V-pinned | 5V R=12 m f=0.448 | pinned | SLS  1.0D + 1.0S | 1.03e-14 | 1.42e-14 | 6.79e-15 | 6.49e-15 |  | yes | 303 free DOF, \|u\|max = 3.96 mm, \|N\|max = 29.45 kN, OS system FullGeneral |
| S-5V-roller3 | 5V R=12 m f=0.448 | 3 rollers | LC1  1.4D | 2.41e-14 | 2.07e-14 | 1.56e-14 | 3.96e-15 |  | yes | 347 free DOF, \|u\|max = 1.666 mm, \|N\|max = 17.43 kN, OS system FullGeneral |
| S-5V-roller3 | 5V R=12 m f=0.448 | 3 rollers | LC2  1.2D + 1.6S | 1.60e-14 | 1.63e-14 | 1.59e-14 | 7.15e-15 |  | yes | 347 free DOF, \|u\|max = 6.495 mm, \|N\|max = 58.36 kN, OS system FullGeneral |
| S-5V-roller3 | 5V R=12 m f=0.448 | 3 rollers | LC3  1.2D + 1.0W + 0.5S | 3.13e-14 | 1.77e-14 | 8.01e-15 | 5.55e-15 |  | yes | 347 free DOF, \|u\|max = 1.929 mm, \|N\|max = 31.17 kN, OS system FullGeneral |
| S-5V-roller3 | 5V R=12 m f=0.448 | 3 rollers | LC4  0.9D + 1.0W | 1.87e-14 | 1.84e-14 | 5.46e-15 | 5.26e-15 |  | yes | 347 free DOF, \|u\|max = 1.3 mm, \|N\|max = 18.16 kN, OS system FullGeneral |
| S-5V-roller3 | 5V R=12 m f=0.448 | 3 rollers | SLS  1.0D + 1.0S | 2.05e-14 | 1.97e-14 | 1.52e-14 | 5.63e-15 |  | yes | 347 free DOF, \|u\|max = 4.357 mm, \|N\|max = 39.52 kN, OS system FullGeneral |
| S-5V-elastic | 5V R=12 m f=0.448 | elastic springs | LC1  1.4D | 1.81e-14 | 2.96e-14 | 7.23e-15 | 8.49e-15 | 7.23e-15 | yes | 378 free DOF, \|u\|max = 1.847 mm, \|N\|max = 11.99 kN, OS system FullGeneral |
| S-5V-elastic | 5V R=12 m f=0.448 | elastic springs | LC2  1.2D + 1.6S | 1.83e-14 | 3.19e-14 | 5.65e-15 | 1.55e-14 | 5.65e-15 | yes | 378 free DOF, \|u\|max = 6.96 mm, \|N\|max = 43.61 kN, OS system FullGeneral |
| S-5V-elastic | 5V R=12 m f=0.448 | elastic springs | LC3  1.2D + 1.0W + 0.5S | 1.98e-14 | 1.92e-14 | 5.21e-15 | 1.26e-14 | 5.21e-15 | yes | 378 free DOF, \|u\|max = 2.01 mm, \|N\|max = 26.29 kN, OS system FullGeneral |
| S-5V-elastic | 5V R=12 m f=0.448 | elastic springs | LC4  0.9D + 1.0W | 1.71e-14 | 1.79e-14 | 9.93e-15 | 1.90e-14 | 9.93e-15 | yes | 378 free DOF, \|u\|max = 1.159 mm, \|N\|max = 13.66 kN, OS system FullGeneral |
| S-5V-elastic | 5V R=12 m f=0.448 | elastic springs | SLS  1.0D + 1.0S | 1.72e-14 | 2.62e-14 | 5.04e-15 | 1.02e-14 | 5.04e-15 | yes | 378 free DOF, \|u\|max = 4.68 mm, \|N\|max = 29.4 kN, OS system FullGeneral |
| S-6V-pinned | 6V R=12 m f=0.5 | pinned | LC1  1.4D | 2.86e-14 | 3.71e-14 | 1.22e-14 | 1.46e-14 |  | yes | 498 free DOF, \|u\|max = 1.534 mm, \|N\|max = 14.56 kN, OS system FullGeneral |
| S-6V-pinned | 6V R=12 m f=0.5 | pinned | LC2  1.2D + 1.6S | 2.47e-14 | 4.09e-14 | 1.75e-14 | 1.54e-14 |  | yes | 498 free DOF, \|u\|max = 5.704 mm, \|N\|max = 44.13 kN, OS system FullGeneral |
| S-6V-pinned | 6V R=12 m f=0.5 | pinned | LC3  1.2D + 1.0W + 0.5S | 3.31e-14 | 3.48e-14 | 1.84e-14 | 2.15e-14 |  | yes | 498 free DOF, \|u\|max = 1.773 mm, \|N\|max = 27.31 kN, OS system FullGeneral |
| S-6V-pinned | 6V R=12 m f=0.5 | pinned | LC4  0.9D + 1.0W | 3.77e-14 | 3.91e-14 | 3.01e-14 | 3.11e-14 |  | yes | 498 free DOF, \|u\|max = 0.9875 mm, \|N\|max = 14.65 kN, OS system FullGeneral |
| S-6V-pinned | 6V R=12 m f=0.5 | pinned | SLS  1.0D + 1.0S | 2.39e-14 | 3.41e-14 | 1.41e-14 | 1.32e-14 |  | yes | 498 free DOF, \|u\|max = 3.839 mm, \|N\|max = 30.18 kN, OS system FullGeneral |
| S-6V-roller3 | 6V R=12 m f=0.5 | 3 rollers | LC1  1.4D | 3.97e-14 | 3.91e-14 | 2.04e-14 | 1.73e-14 |  | yes | 552 free DOF, \|u\|max = 1.593 mm, \|N\|max = 16.64 kN, OS system FullGeneral |
| S-6V-roller3 | 6V R=12 m f=0.5 | 3 rollers | LC2  1.2D + 1.6S | 4.06e-14 | 5.00e-14 | 3.66e-14 | 1.34e-14 |  | yes | 552 free DOF, \|u\|max = 5.861 mm, \|N\|max = 49.64 kN, OS system FullGeneral |
| S-6V-roller3 | 6V R=12 m f=0.5 | 3 rollers | LC3  1.2D + 1.0W + 0.5S | 1.05e-13 | 4.76e-14 | 5.09e-14 | 1.02e-14 |  | yes | 552 free DOF, \|u\|max = 2.003 mm, \|N\|max = 31.03 kN, OS system FullGeneral |
| S-6V-roller3 | 6V R=12 m f=0.5 | 3 rollers | LC4  0.9D + 1.0W | 6.30e-14 | 4.45e-14 | 3.03e-14 | 8.40e-15 |  | yes | 552 free DOF, \|u\|max = 1.642 mm, \|N\|max = 25.96 kN, OS system FullGeneral |
| S-6V-roller3 | 6V R=12 m f=0.5 | 3 rollers | SLS  1.0D + 1.0S | 2.96e-14 | 3.67e-14 | 3.52e-14 | 1.08e-14 |  | yes | 552 free DOF, \|u\|max = 3.947 mm, \|N\|max = 34 kN, OS system FullGeneral |
| S-6V-elastic | 6V R=12 m f=0.5 | elastic springs | LC1  1.4D | 3.90e-14 | 8.38e-14 | 1.37e-14 | 2.38e-14 | 1.37e-14 | yes | 588 free DOF, \|u\|max = 1.891 mm, \|N\|max = 14.67 kN, OS system FullGeneral |
| S-6V-elastic | 6V R=12 m f=0.5 | elastic springs | LC2  1.2D + 1.6S | 4.19e-14 | 7.81e-14 | 2.17e-14 | 2.41e-14 | 2.17e-14 | yes | 588 free DOF, \|u\|max = 6.594 mm, \|N\|max = 44.48 kN, OS system FullGeneral |
| S-6V-elastic | 6V R=12 m f=0.5 | elastic springs | LC3  1.2D + 1.0W + 0.5S | 4.11e-14 | 3.54e-14 | 1.30e-14 | 1.87e-14 | 1.30e-14 | yes | 588 free DOF, \|u\|max = 2.078 mm, \|N\|max = 27.41 kN, OS system FullGeneral |
| S-6V-elastic | 6V R=12 m f=0.5 | elastic springs | LC4  0.9D + 1.0W | 4.26e-14 | 3.23e-14 | 1.28e-14 | 2.59e-14 | 1.28e-14 | yes | 588 free DOF, \|u\|max = 1.12 mm, \|N\|max = 14.6 kN, OS system FullGeneral |
| S-6V-elastic | 6V R=12 m f=0.5 | elastic springs | SLS  1.0D + 1.0S | 3.21e-14 | 6.42e-14 | 1.75e-14 | 2.23e-14 | 1.75e-14 | yes | 588 free DOF, \|u\|max = 4.459 mm, \|N\|max = 30.42 kN, OS system FullGeneral |
| S-7V-pinned | 7V R=12 m f=0.463 | pinned | LC1  1.4D | 4.06e-14 | 5.15e-14 | 1.40e-14 | 1.42e-14 |  | yes | 633 free DOF, \|u\|max = 1.279 mm, \|N\|max = 11.84 kN, OS system FullGeneral |
| S-7V-pinned | 7V R=12 m f=0.463 | pinned | LC2  1.2D + 1.6S | 3.43e-14 | 6.05e-14 | 2.14e-14 | 2.14e-14 |  | yes | 633 free DOF, \|u\|max = 4.701 mm, \|N\|max = 36.84 kN, OS system FullGeneral |
| S-7V-pinned | 7V R=12 m f=0.463 | pinned | LC3  1.2D + 1.0W + 0.5S | 6.84e-14 | 6.60e-14 | 2.16e-14 | 2.54e-14 |  | yes | 633 free DOF, \|u\|max = 1.44 mm, \|N\|max = 22.59 kN, OS system FullGeneral |
| S-7V-pinned | 7V R=12 m f=0.463 | pinned | LC4  0.9D + 1.0W | 6.44e-14 | 5.54e-14 | 3.95e-14 | 3.07e-14 |  | yes | 633 free DOF, \|u\|max = 0.8483 mm, \|N\|max = 11.96 kN, OS system FullGeneral |
| S-7V-pinned | 7V R=12 m f=0.463 | pinned | SLS  1.0D + 1.0S | 2.93e-14 | 4.95e-14 | 1.82e-14 | 1.81e-14 |  | yes | 633 free DOF, \|u\|max = 3.166 mm, \|N\|max = 25.14 kN, OS system FullGeneral |
| S-7V-roller3 | 7V R=12 m f=0.463 | 3 rollers | LC1  1.4D | 6.76e-14 | 5.23e-14 | 3.40e-14 | 1.61e-14 |  | yes | 697 free DOF, \|u\|max = 1.403 mm, \|N\|max = 16.1 kN, OS system FullGeneral |
| S-7V-roller3 | 7V R=12 m f=0.463 | 3 rollers | LC2  1.2D + 1.6S | 7.20e-14 | 9.24e-14 | 4.55e-14 | 1.76e-14 |  | yes | 697 free DOF, \|u\|max = 5.031 mm, \|N\|max = 47.96 kN, OS system FullGeneral |
| S-7V-roller3 | 7V R=12 m f=0.463 | 3 rollers | LC3  1.2D + 1.0W + 0.5S | 8.27e-14 | 6.79e-14 | 3.29e-14 | 2.18e-14 |  | yes | 697 free DOF, \|u\|max = 1.587 mm, \|N\|max = 27.64 kN, OS system FullGeneral |
| S-7V-roller3 | 7V R=12 m f=0.463 | 3 rollers | LC4  0.9D + 1.0W | 1.28e-13 | 7.00e-14 | 2.53e-14 | 1.73e-14 |  | yes | 697 free DOF, \|u\|max = 1.05 mm, \|N\|max = 19.32 kN, OS system FullGeneral |
| S-7V-roller3 | 7V R=12 m f=0.463 | 3 rollers | SLS  1.0D + 1.0S | 8.77e-14 | 7.66e-14 | 3.77e-14 | 1.77e-14 |  | yes | 697 free DOF, \|u\|max = 3.395 mm, \|N\|max = 32.85 kN, OS system FullGeneral |
| S-7V-elastic | 7V R=12 m f=0.463 | elastic springs | LC1  1.4D | 6.75e-14 | 9.56e-14 | 1.83e-14 | 3.02e-14 | 1.83e-14 | yes | 738 free DOF, \|u\|max = 1.581 mm, \|N\|max = 11.84 kN, OS system FullGeneral |
| S-7V-elastic | 7V R=12 m f=0.463 | elastic springs | LC2  1.2D + 1.6S | 5.27e-14 | 8.56e-14 | 2.52e-14 | 4.00e-14 | 2.52e-14 | yes | 738 free DOF, \|u\|max = 5.465 mm, \|N\|max = 36.87 kN, OS system FullGeneral |
| S-7V-elastic | 7V R=12 m f=0.463 | elastic springs | LC3  1.2D + 1.0W + 0.5S | 7.42e-14 | 5.72e-14 | 1.57e-14 | 2.73e-14 | 1.57e-14 | yes | 738 free DOF, \|u\|max = 1.696 mm, \|N\|max = 22.59 kN, OS system FullGeneral |
| S-7V-elastic | 7V R=12 m f=0.463 | elastic springs | LC4  0.9D + 1.0W | 6.27e-14 | 6.95e-14 | 2.54e-14 | 2.63e-14 | 2.54e-14 | yes | 738 free DOF, \|u\|max = 0.9654 mm, \|N\|max = 11.92 kN, OS system FullGeneral |
| S-7V-elastic | 7V R=12 m f=0.463 | elastic springs | SLS  1.0D + 1.0S | 6.05e-14 | 8.99e-14 | 1.88e-14 | 4.09e-14 | 1.88e-14 | yes | 738 free DOF, \|u\|max = 3.698 mm, \|N\|max = 25.16 kN, OS system FullGeneral |
| S-8V-pinned | 8V R=12 m f=0.5 | pinned | LC1  1.4D | 3.34e-14 | 4.04e-14 | 2.24e-14 | 3.55e-15 |  | yes | 903 free DOF, \|u\|max = 1.323 mm, \|N\|max = 13.11 kN, OS system UmfPack |
| S-8V-pinned | 8V R=12 m f=0.5 | pinned | LC2  1.2D + 1.6S | 3.19e-14 | 4.95e-14 | 2.94e-14 | 4.02e-15 |  | yes | 903 free DOF, \|u\|max = 4.56 mm, \|N\|max = 35.76 kN, OS system UmfPack |
| S-8V-pinned | 8V R=12 m f=0.5 | pinned | LC3  1.2D + 1.0W + 0.5S | 5.58e-14 | 4.48e-14 | 2.15e-14 | 4.56e-15 |  | yes | 903 free DOF, \|u\|max = 1.479 mm, \|N\|max = 22.72 kN, OS system UmfPack |
| S-8V-pinned | 8V R=12 m f=0.5 | pinned | LC4  0.9D + 1.0W | 9.68e-14 | 6.87e-14 | 2.96e-14 | 7.58e-15 |  | yes | 903 free DOF, \|u\|max = 0.8033 mm, \|N\|max = 12.57 kN, OS system UmfPack |
| S-8V-pinned | 8V R=12 m f=0.5 | pinned | SLS  1.0D + 1.0S | 2.90e-14 | 4.94e-14 | 1.59e-14 | 4.84e-15 |  | yes | 903 free DOF, \|u\|max = 3.086 mm, \|N\|max = 24.69 kN, OS system UmfPack |
| S-8V-roller3 | 8V R=12 m f=0.5 | 3 rollers | LC1  1.4D | 9.63e-14 | 5.31e-14 | 6.14e-14 | 3.17e-15 |  | yes | 977 free DOF, \|u\|max = 1.361 mm, \|N\|max = 14.82 kN, OS system UmfPack |
| S-8V-roller3 | 8V R=12 m f=0.5 | 3 rollers | LC2  1.2D + 1.6S | 4.10e-14 | 5.44e-14 | 3.21e-14 | 3.29e-15 |  | yes | 977 free DOF, \|u\|max = 4.654 mm, \|N\|max = 39.89 kN, OS system UmfPack |
| S-8V-roller3 | 8V R=12 m f=0.5 | 3 rollers | LC3  1.2D + 1.0W + 0.5S | 6.19e-14 | 6.36e-14 | 2.45e-14 | 1.98e-15 |  | yes | 977 free DOF, \|u\|max = 1.838 mm, \|N\|max = 28.58 kN, OS system UmfPack |
| S-8V-roller3 | 8V R=12 m f=0.5 | 3 rollers | LC4  0.9D + 1.0W | 1.01e-13 | 4.15e-14 | 3.37e-14 | 1.98e-15 |  | yes | 977 free DOF, \|u\|max = 1.549 mm, \|N\|max = 24.96 kN, OS system UmfPack |
| S-8V-roller3 | 8V R=12 m f=0.5 | 3 rollers | SLS  1.0D + 1.0S | 4.26e-14 | 7.56e-14 | 5.90e-14 | 4.32e-15 |  | yes | 977 free DOF, \|u\|max = 3.152 mm, \|N\|max = 27.58 kN, OS system UmfPack |
| S-8V-elastic | 8V R=12 m f=0.5 | elastic springs | LC1  1.4D | 5.01e-14 | 7.61e-14 | 7.54e-15 | 5.68e-15 | 7.54e-15 | yes | 1023 free DOF, \|u\|max = 1.616 mm, \|N\|max = 13.23 kN, OS system UmfPack |
| S-8V-elastic | 8V R=12 m f=0.5 | elastic springs | LC2  1.2D + 1.6S | 3.93e-14 | 8.19e-14 | 1.31e-14 | 7.82e-15 | 1.31e-14 | yes | 1023 free DOF, \|u\|max = 5.247 mm, \|N\|max = 36.11 kN, OS system UmfPack |
| S-8V-elastic | 8V R=12 m f=0.5 | elastic springs | LC3  1.2D + 1.0W + 0.5S | 2.19e-13 | 1.72e-13 | 4.59e-14 | 6.68e-15 | 4.59e-14 | yes | 1023 free DOF, \|u\|max = 1.728 mm, \|N\|max = 22.77 kN, OS system UmfPack |
| S-8V-elastic | 8V R=12 m f=0.5 | elastic springs | LC4  0.9D + 1.0W | 1.11e-13 | 8.83e-14 | 4.35e-14 | 1.14e-14 | 4.35e-14 | yes | 1023 free DOF, \|u\|max = 0.9236 mm, \|N\|max = 12.4 kN, OS system UmfPack |
| S-8V-elastic | 8V R=12 m f=0.5 | elastic springs | SLS  1.0D + 1.0S | 4.88e-14 | 6.08e-14 | 1.08e-14 | 5.56e-15 | 1.08e-14 | yes | 1023 free DOF, \|u\|max = 3.568 mm, \|N\|max = 24.93 kN, OS system UmfPack |
| S-validation_a | 1V R=5 m f=0.276 | pinned | A1  1.0L | 0 | 2.26e-16 | 6.06e-17 | 2.12e-16 |  | yes | 3 free DOF, \|u\|max = 1.07 mm, \|N\|max = 32.23 kN, OS system FullGeneral |
| S-validation_a | 1V R=5 m f=0.276 | pinned | A2  1.2D + 1.6S | 0 | 2.46e-16 | 9.35e-17 | 2.34e-16 |  | yes | 3 free DOF, \|u\|max = 0.6937 mm, \|N\|max = 14.8 kN, OS system FullGeneral |
| S-validation_b | 1V R=5 m f=0.724 | pinned | B1  1.0L | 1.97e-16 | 3.83e-16 | 1.82e-16 | 1.53e-16 |  | yes | 18 free DOF, \|u\|max = 1.645 mm, \|N\|max = 19.02 kN, OS system FullGeneral |
| S-validation_b | 1V R=5 m f=0.724 | pinned | B2  1.2D + 1.6S | 2.72e-16 | 5.67e-16 | 1.72e-16 | 2.15e-16 |  | yes | 18 free DOF, \|u\|max = 1.859 mm, \|N\|max = 19.25 kN, OS system FullGeneral |
| S-three_bar | Planar three-bar indeterminate truss | custom | P  1.0L | 0 | 2.58e-17 | 3.03e-17 | 1.82e-16 | 0 | yes | 2 free DOF, \|u\|max = 0.5838 mm, \|N\|max = 35.26 kN, OS system FullGeneral |
| S-three_bar | Planar three-bar indeterminate truss | custom | T  1.0T  temperature + settlement | 3.25e-17 | 6.80e-16 | 5.19e-16 | 1.25e-15 | 0 | yes | 2 free DOF, \|u\|max = 5 mm, \|N\|max = 35.02 kN, OS system FullGeneral |
| S-three_bar | Planar three-bar indeterminate truss | custom | P + T  1.0L + 1.0T | 2.98e-17 | 9.71e-17 | 2.07e-16 | 4.93e-16 | 0 | yes | 2 free DOF, \|u\|max = 5 mm, \|N\|max = 70.28 kN, OS system FullGeneral |
| S-case_study_6v | 6V R=15 m f=0.5 | pinned | LC1  1.4D | 2.49e-14 | 4.33e-14 | 9.91e-15 | 1.23e-14 |  | yes | 498 free DOF, \|u\|max = 2.238 mm, \|N\|max = 27.31 kN, OS system FullGeneral |
| S-case_study_6v | 6V R=15 m f=0.5 | pinned | LC2  1.2D + 1.6S | 3.26e-14 | 5.49e-14 | 1.41e-14 | 1.49e-14 |  | yes | 498 free DOF, \|u\|max = 4.582 mm, \|N\|max = 48.14 kN, OS system FullGeneral |
| S-case_study_6v | 6V R=15 m f=0.5 | pinned | LC3  1.2D + 1.0W + 0.5S | 4.76e-14 | 2.87e-14 | 2.17e-14 | 2.92e-14 |  | yes | 498 free DOF, \|u\|max = 2.332 mm, \|N\|max = 45.39 kN, OS system FullGeneral |
| S-case_study_6v | 6V R=15 m f=0.5 | pinned | LC4  0.9D + 1.0W | 4.93e-14 | 4.01e-14 | 2.07e-14 | 2.03e-14 |  | yes | 498 free DOF, \|u\|max = 2.043 mm, \|N\|max = 32.32 kN, OS system FullGeneral |
| S-case_study_6v | 6V R=15 m f=0.5 | pinned | SLS  1.0D + 1.0S | 3.26e-14 | 4.78e-14 | 1.93e-14 | 1.44e-14 |  | yes | 498 free DOF, \|u\|max = 3.264 mm, \|N\|max = 34.97 kN, OS system FullGeneral |
| S-4V-T-pinned | 4V R=12 m f=0.5 | pinned | T1  1.2D + 1.6S + 1.0T | 6.45e-15 | 3.21e-15 | 7.77e-16 | 1.98e-15 |  | yes | 213 free DOF, \|u\|max = 10.08 mm, \|N\|max = 278.4 kN, OS system FullGeneral |
| S-4V-T-pinned | 4V R=12 m f=0.5 | pinned | T2  1.0T | 4.94e-15 | 2.38e-15 | 6.87e-16 | 1.10e-15 |  | yes | 213 free DOF, \|u\|max = 10 mm, \|N\|max = 278.4 kN, OS system FullGeneral |
| S-4V-T-pinned | 4V R=12 m f=0.5 | pinned | T3  1.2D + 1.0W + 0.5S + 1.2T | 5.71e-15 | 3.00e-15 | 6.62e-16 | 7.72e-16 |  | yes | 213 free DOF, \|u\|max = 12 mm, \|N\|max = 334.1 kN, OS system FullGeneral |
| S-4V-T-roller3 | 4V R=12 m f=0.5 | 3 rollers | T1  1.2D + 1.6S + 1.0T | 1.10e-14 | 2.01e-14 | 1.69e-14 | 4.13e-15 |  | yes | 247 free DOF, \|u\|max = 27.37 mm, \|N\|max = 145.2 kN, OS system FullGeneral |
| S-4V-T-roller3 | 4V R=12 m f=0.5 | 3 rollers | T2  1.0T | 6.65e-15 | 2.92e-14 | 1.96e-14 | 4.30e-15 |  | yes | 247 free DOF, \|u\|max = 29.47 mm, \|N\|max = 99.29 kN, OS system FullGeneral |
| S-4V-T-roller3 | 4V R=12 m f=0.5 | 3 rollers | T3  1.2D + 1.0W + 0.5S + 1.2T | 5.52e-15 | 1.61e-14 | 8.60e-15 | 4.02e-15 |  | yes | 247 free DOF, \|u\|max = 33.6 mm, \|N\|max = 142 kN, OS system FullGeneral |
| S-4V-T-elastic | 4V R=12 m f=0.5 | elastic springs | T1  1.2D + 1.6S + 1.0T | 9.92e-15 | 9.34e-15 | 1.75e-15 | 4.74e-15 | 1.75e-15 | yes | 273 free DOF, \|u\|max = 8.835 mm, \|N\|max = 119.2 kN, OS system FullGeneral |
| S-4V-T-elastic | 4V R=12 m f=0.5 | elastic springs | T2  1.0T | 8.44e-15 | 5.70e-15 | 1.37e-15 | 1.87e-15 | 1.37e-15 | yes | 273 free DOF, \|u\|max = 6.163 mm, \|N\|max = 118.9 kN, OS system FullGeneral |
| S-4V-T-elastic | 4V R=12 m f=0.5 | elastic springs | T3  1.2D + 1.0W + 0.5S + 1.2T | 1.10e-14 | 7.55e-15 | 1.69e-15 | 2.59e-15 | 1.69e-15 | yes | 273 free DOF, \|u\|max = 8.654 mm, \|N\|max = 143.1 kN, OS system FullGeneral |
| S-3V-custom | 3V R=10 m f=0.414 | custom | LC1  1.4D | 2.62e-15 | 2.87e-15 | 1.20e-15 | 1.31e-15 | 4.51e-16 | yes | 118 free DOF, \|u\|max = 1.198 mm, \|N\|max = 13.3 kN, OS system FullGeneral |
| S-3V-custom | 3V R=10 m f=0.414 | custom | LC2  1.2D + 1.6S | 2.48e-15 | 3.78e-15 | 1.18e-15 | 1.47e-15 | 5.92e-16 | yes | 118 free DOF, \|u\|max = 5.081 mm, \|N\|max = 44.25 kN, OS system FullGeneral |
| S-3V-custom | 3V R=10 m f=0.414 | custom | LC3  1.2D + 1.0W + 0.5S | 3.08e-15 | 2.91e-15 | 8.72e-16 | 1.27e-15 | 4.36e-16 | yes | 118 free DOF, \|u\|max = 1.266 mm, \|N\|max = 22.2 kN, OS system FullGeneral |
| S-3V-custom | 3V R=10 m f=0.414 | custom | LC4  0.9D + 1.0W | 2.48e-15 | 3.27e-15 | 1.03e-15 | 2.03e-15 | 2.79e-16 | yes | 118 free DOF, \|u\|max = 0.8085 mm, \|N\|max = 12.26 kN, OS system FullGeneral |
| S-3V-custom | 3V R=10 m f=0.414 | custom | SLS  1.0D + 1.0S | 2.88e-15 | 3.39e-15 | 1.03e-15 | 1.48e-15 | 5.06e-16 | yes | 118 free DOF, \|u\|max = 3.39 mm, \|N\|max = 30.03 kN, OS system FullGeneral |
| S-3V-custom | 3V R=10 m f=0.414 | custom | T1  1.2D + 1.6S + 1.0T | 3.37e-15 | 4.05e-15 | 1.17e-15 | 1.83e-15 | 4.43e-16 | yes | 118 free DOF, \|u\|max = 9.272 mm, \|N\|max = 155.6 kN, OS system FullGeneral |
| S-5V-f0.625 | 5V R=15 m f=0.648 | pinned | LC1  1.4D | 2.91e-14 | 2.45e-14 | 7.39e-15 | 9.65e-15 |  | yes | 393 free DOF, \|u\|max = 3.731 mm, \|N\|max = 30.41 kN, OS system FullGeneral |
| S-5V-f0.625 | 5V R=15 m f=0.648 | pinned | LC2  1.2D + 1.6S | 1.47e-14 | 2.11e-14 | 8.77e-15 | 6.85e-15 |  | yes | 393 free DOF, \|u\|max = 17.2 mm, \|N\|max = 113.2 kN, OS system FullGeneral |
| S-5V-f0.625 | 5V R=15 m f=0.648 | pinned | LC3  1.2D + 1.0W + 0.5S | 1.94e-14 | 1.63e-14 | 8.03e-15 | 7.34e-15 |  | yes | 393 free DOF, \|u\|max = 5.501 mm, \|N\|max = 63.74 kN, OS system FullGeneral |
| S-5V-f0.625 | 5V R=15 m f=0.648 | pinned | LC4  0.9D + 1.0W | 3.89e-14 | 2.48e-14 | 8.43e-15 | 1.52e-14 |  | yes | 393 free DOF, \|u\|max = 2.421 mm, \|N\|max = 30.42 kN, OS system FullGeneral |
| S-5V-f0.625 | 5V R=15 m f=0.648 | pinned | SLS  1.0D + 1.0S | 1.36e-14 | 1.79e-14 | 1.08e-14 | 7.52e-15 |  | yes | 393 free DOF, \|u\|max = 11.42 mm, \|N\|max = 76.2 kN, OS system FullGeneral |
| S-3V-f0.375 | 3V R=8 m f=0.406 | pinned | LC1  1.4D | 1.80e-15 | 2.54e-15 | 4.16e-16 | 8.32e-16 |  | yes | 48 free DOF, \|u\|max = 0.3908 mm, \|N\|max = 5.545 kN, OS system FullGeneral |
| S-3V-f0.375 | 3V R=8 m f=0.406 | pinned | LC2  1.2D + 1.6S | 1.74e-15 | 2.09e-15 | 5.75e-16 | 8.63e-16 |  | yes | 48 free DOF, \|u\|max = 2.495 mm, \|N\|max = 33.01 kN, OS system FullGeneral |
| S-3V-f0.375 | 3V R=8 m f=0.406 | pinned | LC3  1.2D + 1.0W + 0.5S | 1.41e-15 | 1.86e-15 | 9.21e-16 | 9.90e-16 |  | yes | 48 free DOF, \|u\|max = 0.5785 mm, \|N\|max = 8.323 kN, OS system FullGeneral |
| S-3V-f0.375 | 3V R=8 m f=0.406 | pinned | LC4  0.9D + 1.0W | 2.28e-15 | 1.66e-15 | 1.28e-15 | 1.48e-15 |  | yes | 48 free DOF, \|u\|max = 0.2258 mm, \|N\|max = 6.196 kN, OS system FullGeneral |
| S-3V-f0.375 | 3V R=8 m f=0.406 | pinned | SLS  1.0D + 1.0S | 7.99e-16 | 1.18e-15 | 5.76e-16 | 1.12e-15 |  | yes | 48 free DOF, \|u\|max = 1.629 mm, \|N\|max = 21.58 kN, OS system FullGeneral |

## M. Modal analysis

Mass: strut mass rho*A*L/2 to each end (OpenSees `-rho`), or the consistent bar mass (`-cMass 1`), plus the cladding dead load / g at the nodes (OpenSees `mass`).  'all' compares every mode of the model (OpenSees `-fullGenLapack`).  OpenSees mode shapes are compared without any normalisation (MAC is scale-free).

| case | structure | supports | combination | natural frequencies f | 1 - MAC (single modes) | 1 - subspace MAC (repeated pairs) | pass | note |
|---|---|---|---|---:|---:|---:|---|---|
| M-1V-pinned | 1V R=12 m f=0.276 | pinned | lumped mass | 0 | 0 | 0 | yes | 3 modes (1 single, 2 in repeated clusters), f1 = 12.516074 Hz, OS -fullGenLapack |
| M-1V-roller3 | 1V R=12 m f=0.276 | 3 rollers | lumped mass | 3.44e-16 | 1.11e-16 | 0 | yes | 7 modes (5 single, 2 in repeated clusters), f1 = 10.100402 Hz, OS -fullGenLapack |
| M-1V-elastic | 1V R=12 m f=0.276 | elastic springs | lumped mass | 1.29e-15 | 2.22e-16 | 2.22e-16 | yes | 12 modes (2 single, 10 in repeated clusters), f1 = 10.996040 Hz, OS -genBandArpack |
| M-2V-pinned | 2V R=12 m f=0.5 | pinned | lumped mass | 1.34e-15 | 2.22e-16 | 9.99e-16 | yes | 12 modes (2 single, 10 in repeated clusters), f1 = 7.972006 Hz, OS -genBandArpack |
| M-2V-roller3 | 2V R=12 m f=0.5 | 3 rollers | lumped mass | 1.87e-15 | 4.44e-16 |  | yes | 12 modes (12 single, 0 in repeated clusters), f1 = 5.577929 Hz, OS -genBandArpack |
| M-2V-elastic | 2V R=12 m f=0.5 | elastic springs | lumped mass | 2.71e-15 | 5.55e-16 | 9.99e-16 | yes | 12 modes (2 single, 10 in repeated clusters), f1 = 7.222081 Hz, OS -genBandArpack |
| M-3V-pinned | 3V R=12 m f=0.414 | pinned | lumped mass | 7.45e-15 | 4.44e-16 | 1.55e-15 | yes | 12 modes (2 single, 10 in repeated clusters), f1 = 10.727312 Hz, OS -genBandArpack |
| M-3V-roller3 | 3V R=12 m f=0.414 | 3 rollers | lumped mass | 5.96e-15 | 6.66e-16 |  | yes | 12 modes (12 single, 0 in repeated clusters), f1 = 7.155026 Hz, OS -genBandArpack |
| M-3V-elastic | 3V R=12 m f=0.414 | elastic springs | lumped mass | 9.09e-15 | 5.55e-16 | 1.67e-15 | yes | 12 modes (2 single, 10 in repeated clusters), f1 = 9.774876 Hz, OS -genBandArpack |
| M-4V-pinned | 4V R=12 m f=0.5 | pinned | lumped mass | 1.68e-14 | 7.77e-16 | 1.78e-15 | yes | 12 modes (2 single, 10 in repeated clusters), f1 = 10.040369 Hz, OS -genBandArpack |
| M-4V-roller3 | 4V R=12 m f=0.5 | 3 rollers | lumped mass | 5.68e-14 | 6.66e-16 |  | yes | 12 modes (12 single, 0 in repeated clusters), f1 = 5.972223 Hz, OS -genBandArpack |
| M-4V-elastic | 4V R=12 m f=0.5 | elastic springs | lumped mass | 1.66e-14 | 1.67e-15 | 7.77e-16 | yes | 12 modes (2 single, 10 in repeated clusters), f1 = 8.964835 Hz, OS -genBandArpack |
| M-5V-pinned | 5V R=12 m f=0.448 | pinned | lumped mass | 1.23e-14 | 7.77e-16 | 1.22e-15 | yes | 12 modes (2 single, 10 in repeated clusters), f1 = 12.046212 Hz, OS -genBandArpack |
| M-5V-roller3 | 5V R=12 m f=0.448 | 3 rollers | lumped mass | 8.15e-14 | 8.88e-16 |  | yes | 12 modes (12 single, 0 in repeated clusters), f1 = 7.575058 Hz, OS -genBandArpack |
| M-5V-elastic | 5V R=12 m f=0.448 | elastic springs | lumped mass | 1.34e-14 | 8.88e-16 | 2.78e-15 | yes | 12 modes (2 single, 10 in repeated clusters), f1 = 10.800846 Hz, OS -genBandArpack |
| M-6V-pinned | 6V R=12 m f=0.5 | pinned | lumped mass | 3.41e-14 | 8.88e-16 | 1.89e-15 | yes | 12 modes (2 single, 10 in repeated clusters), f1 = 11.500382 Hz, OS -genBandArpack |
| M-6V-roller3 | 6V R=12 m f=0.5 | 3 rollers | lumped mass | 2.79e-14 | 1.22e-15 |  | yes | 12 modes (12 single, 0 in repeated clusters), f1 = 6.236573 Hz, OS -genBandArpack |
| M-6V-elastic | 6V R=12 m f=0.5 | elastic springs | lumped mass | 5.57e-14 | 7.77e-16 | 1.89e-15 | yes | 12 modes (2 single, 10 in repeated clusters), f1 = 10.204937 Hz, OS -genBandArpack |
| M-7V-pinned | 7V R=12 m f=0.463 | pinned | lumped mass | 2.89e-14 | 2.22e-16 | 1.33e-15 | yes | 12 modes (2 single, 10 in repeated clusters), f1 = 13.063831 Hz, OS -genBandArpack |
| M-7V-roller3 | 7V R=12 m f=0.463 | 3 rollers | lumped mass | 1.16e-13 | 9.99e-16 |  | yes | 12 modes (12 single, 0 in repeated clusters), f1 = 7.701355 Hz, OS -genBandArpack |
| M-7V-elastic | 7V R=12 m f=0.463 | elastic springs | lumped mass | 8.20e-14 | 5.55e-16 | 1.78e-15 | yes | 12 modes (2 single, 10 in repeated clusters), f1 = 11.632923 Hz, OS -genBandArpack |
| M-8V-pinned | 8V R=12 m f=0.5 | pinned | lumped mass | 1.26e-13 | 0 | 3.77e-15 | yes | 12 modes (2 single, 10 in repeated clusters), f1 = 12.592208 Hz, OS -genBandArpack |
| M-8V-roller3 | 8V R=12 m f=0.5 | 3 rollers | lumped mass | 2.72e-13 | 1.55e-15 |  | yes | 12 modes (12 single, 0 in repeated clusters), f1 = 6.793846 Hz, OS -genBandArpack |
| M-8V-elastic | 8V R=12 m f=0.5 | elastic springs | lumped mass | 4.39e-14 | 3.33e-15 | 4.00e-15 | yes | 12 modes (2 single, 10 in repeated clusters), f1 = 11.134832 Hz, OS -genBandArpack |
| M-1V-all-lumped | 1V R=12 m f=0.276 | pinned | lumped mass | 0 | 0 | 0 | yes | 3 modes (1 single, 2 in repeated clusters), f1 = 12.516074 Hz, OS -fullGenLapack |
| M-1V-all-consistent | 1V R=12 m f=0.276 | pinned | consistent mass | 1.40e-16 | 0 | 0 | yes | 3 modes (1 single, 2 in repeated clusters), f1 = 12.664503 Hz, OS -fullGenLapack |
| M-2V-all-lumped | 2V R=12 m f=0.5 | pinned | lumped mass | 2.34e-15 | 3.33e-16 | 1.11e-15 | yes | 48 modes (10 single, 38 in repeated clusters), f1 = 7.972006 Hz, OS -fullGenLapack |
| M-2V-all-consistent | 2V R=12 m f=0.5 | pinned | consistent mass | 3.10e-15 | 4.44e-16 | 7.77e-16 | yes | 48 modes (10 single, 38 in repeated clusters), f1 = 8.005464 Hz, OS -fullGenLapack |
| M-3V-all-lumped | 3V R=12 m f=0.414 | pinned | lumped mass | 1.08e-14 | 4.44e-16 | 1.33e-15 | yes | 93 modes (19 single, 74 in repeated clusters), f1 = 10.727312 Hz, OS -fullGenLapack |
| M-3V-all-consistent | 3V R=12 m f=0.414 | pinned | consistent mass | 1.62e-14 | 6.66e-16 | 9.99e-16 | yes | 93 modes (19 single, 74 in repeated clusters), f1 = 10.770052 Hz, OS -fullGenLapack |
| M-case6V-lumped | 6V R=15 m f=0.5 | pinned | lumped mass | 2.57e-14 | 9.99e-16 | 2.22e-15 | yes | 20 modes (4 single, 16 in repeated clusters), f1 = 9.525865 Hz, OS -genBandArpack |
| M-case6V-consistent | 6V R=15 m f=0.5 | pinned | consistent mass | 5.16e-14 | 2.22e-16 | 2.55e-15 | yes | 20 modes (4 single, 16 in repeated clusters), f1 = 9.542726 Hz, OS -genBandArpack |
| M-validation_a | 1V R=5 m f=0.276 | pinned | lumped mass | 1.40e-16 | 0 | 0 | yes | 3 modes (1 single, 2 in repeated clusters), f1 = 44.412911 Hz, OS -fullGenLapack |
| M-star_dome | 24-member shallow star dome | pinned | lumped mass | 1.25e-14 | 1.11e-16 | 4.44e-16 | yes | 21 modes (7 single, 14 in repeated clusters), f1 = 274.324562 Hz, OS -fullGenLapack |
| M-three_bar | Planar three-bar indeterminate truss | custom | lumped mass | 1.25e-16 | 0 |  | yes | 2 modes (2 single, 0 in repeated clusters), f1 = 163.357857 Hz, OS -fullGenLapack |

## N. Geometrically nonlinear analysis

GeoDomeLab traces the path (cylindrical arc-length, Crisfield 1981) and refines the first limit point.  OpenSees (`corotTruss`, Newton, unbalance tolerance 1e-10 |F_ref|) is then driven by `DisplacementControl` to **every GeoDomeLab state**, including the refined limit point: where the path snaps back in the control displacement, control switches to another DOF that is monotone over the next stretch (piecewise displacement control, `opensees_bridge.path_segments`).  At each state OpenSees finds lambda and all displacements by itself; the table compares them.  In addition every OpenSees state is put into GeoDomeLab's own corotational equilibrium equations (residual row).

OpenSees' own `ArcLength` integrator (fixed arc length 2 mm, alpha = 0, cylindrical) was also run on the benchmarks as an independent path tracer.  Its predictor takes the sign of the previous load increment (ArcLength.cpp), which is not a reliable direction test at a limit point: N-von-mises: it went as far as -198 mm in the loading direction (lambda = -1.0238 there; GeoDomeLab's first limit point is at -52.86 mm) and then turned back, ending at -88 mm with lambda = 0.71886 after 300 steps. N-star-dome: it went as far as -7.879 mm in the loading direction (lambda = 3.155 there; GeoDomeLab's first limit point is at -7.684 mm) and then turned back, ending at 170.7 mm with lambda = -2885.4 after 125 steps.  Every point it produced still satisfies GeoDomeLab's equilibrium equations (residual row), but because it can turn back at a limit point instead of continuing along the path, the matched comparison drives OpenSees by displacement control instead.

| case | structure | supports | combination | lambda at matched states / max\|lambda\| | lambda pointwise (\|lambda\| > 5 % of max) | displacements at matched states | lambda at the refined first limit point | GDL equilibrium residual of OpenSees states | GDL residual of OpenSees ArcLength points | OpenSees P(w) vs exact closed form | OpenSees limit load vs closed form P* | OpenSees P/EA at limit vs reference 3.15655e-4 (absolute, 6 digits) | pass | note |
|---|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---|
| N-von-mises | von Mises two-bar truss | custom | P  1.0L  apex reference load | 9.60e-14 | 4.48e-13 | 0 | 3.77e-14 | 3.73e-13 | 3.14e-13 | 9.60e-14 | 1.28e-14 |  | yes | 22/22 states, 1 control segment(s) on DOF(s) [8]; lambda in [-1.022, 3.977], control in [-303.9, 0] mm, max \|strain\| 1.31e-03; \|control displacement\| exceeded 0.3 m; lambda_limit GDL 1.02404996, OS 1.02404996 (limit); OpenSees states checked with GeoDomeLab's corotational equations; 301 independent OpenSees ArcLength points, control in [-198, 0] mm, max lambda 1.02383 (step limit reached); 22 points, w up to 303.9 mm (h = 125 mm); P* = 10.2405 kN at w* = 52.8613 mm |
| N-von-mises-deep | von Mises two-bar truss | custom | P  1.0L  apex reference load | 1.04e-11 | 3.16e-11 | 0 | 3.16e-11 | 7.38e-13 |  | 2.16e-13 | 1.05e-12 |  | yes | 13/13 states, 1 control segment(s) on DOF(s) [8]; lambda in [-38.51, 120.9], control in [-2600, 0] mm, max \|strain\| 3.34e-01; \|control displacement\| exceeded 2.4 m; lambda_limit GDL 39.9875436, OS 39.9875436 (limit); OpenSees states checked with GeoDomeLab's corotational equations; 13 points, w up to 2600 mm (h = 1000 mm); P* = 39987.544 kN at w* = 490.175 mm |
| N-star-dome | 24-member shallow star dome | pinned | P  1.0L  crown reference load | 8.93e-13 | 3.65e-12 | 2.66e-13 | 2.80e-14 | 3.37e-11 | 6.38e-12 |  |  | 4.04e-10 | yes | 284/284 states, 5 control segment(s) on DOF(s) [2, 3, 5]; lambda in [-88.65, 88.64], control in [-171.8, 0] mm, max \|strain\| 1.69e-02; \|control displacement\| exceeded 0.17 m; lambda_limit GDL 3.156545962, OS 3.156545962 (limit); OpenSees states checked with GeoDomeLab's corotational equations; 126 independent OpenSees ArcLength points, control in [-7.879, 170.7] mm, max lambda 3.15497 (control limit reached); P/EA = 3.156545962e-04; reference = research notes 5d (agrees with Greco & Venturini 2006, Fig. 8) |
| N-2V-LC2 | 2V R=12 m f=0.5 | pinned | LC2  1.2D + 1.6S | 1.24e-13 | 1.84e-13 | 6.85e-13 | 0 | 8.27e-13 |  |  |  |  | yes | 24/24 states, 1 control segment(s) on DOF(s) [2]; lambda in [0, 86.73], control in [-2637, 0] mm, max \|strain\| 9.22e-02; first limit point passed (10 steps traced beyond it); lambda_limit GDL 86.73435804, OS 86.73435804 (limit); OpenSees states checked with GeoDomeLab's corotational equations |
| N-3V-LC2 | 3V R=12 m f=0.414 | pinned | LC2  1.2D + 1.6S | 4.79e-13 | 9.26e-13 | 5.27e-13 | 6.27e-15 | 1.53e-12 |  |  |  |  | yes | 23/23 states, 1 control segment(s) on DOF(s) [2]; lambda in [0, 48.72], control in [-578.7, 0] mm, max \|strain\| 1.49e-02; first limit point passed (10 steps traced beyond it); lambda_limit GDL 48.71855229, OS 48.71855229 (limit); OpenSees states checked with GeoDomeLab's corotational equations |
| N-4V-LC2 | 4V R=12 m f=0.5 | pinned | LC2  1.2D + 1.6S | 1.45e-11 | 1.45e-11 | 1.40e-11 | 3.04e-14 | 1.57e-12 |  |  |  |  | yes | 22/22 states, 2 control segment(s) on DOF(s) [2, 3]; lambda in [0, 36.45], control in [-356.2, 0] mm, max \|strain\| 1.22e-02; first limit point passed (10 steps traced beyond it); lambda_limit GDL 36.4518813, OS 36.4518813 (limit); OpenSees states checked with GeoDomeLab's corotational equations |
| N-3V-elastic-LC2 | 3V R=12 m f=0.414 | elastic springs | LC2  1.2D + 1.6S | 2.86e-12 | 2.86e-12 | 2.50e-12 | 6.54e-15 | 7.71e-13 |  |  |  |  | yes | 21/21 states, 2 control segment(s) on DOF(s) [2, 3]; lambda in [0, 48.89], control in [-622.8, 0] mm, max \|strain\| 1.49e-02; first limit point passed (10 steps traced beyond it); lambda_limit GDL 48.90615315, OS 48.90615315 (limit); OpenSees states checked with GeoDomeLab's corotational equations |
| N-case6V-LC2 | 6V R=15 m f=0.5 | pinned | LC2  1.2D + 1.6S | 4.14e-11 | 4.38e-11 | 3.86e-11 | 2.57e-14 | 6.45e-12 |  |  |  |  | yes | 23/23 states, 1 control segment(s) on DOF(s) [2]; lambda in [0, 34.24], control in [-193.5, 0] mm, max \|strain\| 5.85e-03; first limit point passed (10 steps traced beyond it); lambda_limit GDL 34.23695771, OS 34.23695771 (limit); OpenSees states checked with GeoDomeLab's corotational equations |

Star dome, a sample of matched states (lambda = 10^4 P/EA):

| state | crown w [mm] | lambda GeoDomeLab | lambda OpenSees | difference |
|---|---:|---:|---:|---:|
| 1 | -0.4203 | 0.3671144888 | 0.3671144888 | 2.9e-14 |
| limit | -7.6844 | 3.1565459615 | 3.1565459615 | 8.8e-14 |
| 20 | -37.3156 | -1.3364927800 | -1.3364927800 | 1.3e-11 |
| 40 | -71.7513 | 44.0725174585 | 44.0725174585 | 1.1e-13 |
| 60 | -96.7417 | 84.2540257104 | 84.2540257104 | 5.5e-13 |
| 80 | -115.3256 | 80.5816363140 | 80.5816363140 | 1.9e-12 |
| 100 | -127.5603 | 30.7330392497 | 30.7330392497 | 1.4e-12 |
| 120 | -120.3049 | -46.8765845669 | -46.8765845669 | 2.6e-13 |
| 140 | -83.5504 | -2.5791775270 | -2.5791775270 | 3.5e-14 |
| 160 | -47.4489 | 47.3989289206 | 47.3989289206 | 4.1e-11 |
| 180 | -36.1430 | -24.9226007496 | -24.9226007496 | 5.2e-13 |
| 200 | -47.7372 | -78.1472574649 | -78.1472574649 | 3.8e-13 |
| 220 | -65.8412 | -85.7696782090 | -85.7696782090 | 3.8e-13 |
| 240 | -90.2384 | -48.4675505232 | -48.4675505232 | 2.5e-12 |
| 260 | -123.7397 | -0.3480338777 | -0.3480338777 | 3.3e-13 |
| 280 | -165.5019 | 1.1739767677 | 1.1739767677 | 9.7e-14 |
| 283 | -171.8360 | 11.3607106962 | 11.3607106962 | 5.8e-12 |

## X. Exported standalone scripts

`opensees_bridge.write_openseespy_script` writes a readable OpenSeesPy script (explicit `ops.node`, `ops.element`, ... lines) that anyone can run with `pip install openseespy`.  The scripts below were written for a 3V dome on elastic supports with temperature and settlement, run in the OpenSees interpreter and compared with GeoDomeLab.

| case | structure | quantity | difference | target | pass | note |
|---|---|---|---:|---:|---|---|
| X1 | 3V R=10 m f=0.414 | exported static script: u | 6.32e-15 | <= 1e-08 | yes | write_openseespy_script + run_script |
| X1 | 3V R=10 m f=0.414 | exported static script: N | 5.79e-15 | <= 1e-08 | yes |  |
| X1 | 3V R=10 m f=0.414 | exported static script: R | 2.00e-15 | <= 1e-08 | yes |  |
| X2 | 3V R=10 m f=0.414 | exported modal script: f (10 modes) | 2.76e-15 | <= 1e-06 | yes |  |
| X3 | 3V R=10 m f=0.414 | exported nonlinear script: GDL residual of its states | 9.86e-13 | <= 1e-08 | yes | 20 DisplacementControl states to 400 mm, lambda up to 60.633 |

## Notes and limits

* The differences are at round-off level (1e-16 to 1e-11), many orders of magnitude below the targets.  The largest ones occur on the biggest models and along the nonlinear paths, where GeoDomeLab's own convergence tolerance (1e-8 on the residual) is the limit, not the formulation.
* What OpenSees cannot check: the generation of the loads and masses (inputs here; see the hand-calculation validation), the AISC member checks, the linear buckling eigenvalue (OpenSees has no built-in linearised buckling analysis; GeoDomeLab's is checked against closed forms in the benchmark suite), and the fabrication geometry.
* The nonlinear path excludes the self-straining load T in both programs (GeoDomeLab's documented scope).
* Mode shapes inside a repeated-frequency cluster are compared as subspaces because any rotation of a repeated pair is an equally valid pair of modes.

References: F. McKenna, M. H. Scott, G. L. Fenves (2010), J. Comput. Civ. Eng. 24(1), 95-107; M. Zhu, F. McKenna, M. H. Scott (2018), SoftwareX 7, 6-11; OpenSees source files CorotTruss.cpp and InitStrainMaterial.cpp; M. A. Crisfield (1981), Comput. Struct. 13, 55-62 and (1991) Non-linear Finite Element Analysis of Solids and Structures, Vol. 1, Wiley; R. J. Allemang (2003), Sound and Vibration 37(8), 14-23 (MAC); Pecknold, Ghaboussi & Healey (1985), J. Eng. Mech. 111(7), 909-922 (von Mises truss); Greco & Venturini (2006) (star dome).

Full rows: `validation/output/opensees_crosscheck.csv` (835 rows).
