"""Unit conversion constants.

Everything inside GeoDomeLab is SI (m, N, Pa, kg).  These helpers exist only so
that the GUI and the printed reports can show engineering-friendly numbers.
"""

# ---------------------------------------------------------------- length ----
M_TO_MM = 1.0e3
MM_TO_M = 1.0e-3
M_TO_FT = 3.280839895
M_TO_IN = 39.37007874

# ----------------------------------------------------------------- force ----
N_TO_KN = 1.0e-3
KN_TO_N = 1.0e3
N_TO_KIP = 2.248089431e-4

# ---------------------------------------------------------------- stress ----
PA_TO_MPA = 1.0e-6
MPA_TO_PA = 1.0e6
PA_TO_KSI = 1.450377377e-7

# -------------------------------------------------------------- pressure ----
PA_TO_KPA = 1.0e-3
KPA_TO_PA = 1.0e3          # 1 kN/m^2 == 1 kPa
PSF_TO_PA = 47.88025898

# -------------------------------------------------------------- stiffness ----
# Spring supports are entered in N/m internally; bearing data sheets quote
# kN/mm (= 1e6 N/m) or kip/in (= 4448.2216 N / 0.0254 m).
KN_PER_MM_TO_N_PER_M = 1.0e6
N_PER_M_TO_KN_PER_MM = 1.0e-6
KIP_PER_IN_TO_N_PER_M = 175126.8352

# ------------------------------------------------------------ temperature ----
# Temperature *differences* only (the self-straining load uses dT, so no
# offset is involved).  A coefficient of expansion per deg F times 1.8 gives
# the value per deg C: 6.5e-6 /F (AISC) = 11.7e-6 /C.
DELTA_F_TO_DELTA_C = 5.0 / 9.0
DELTA_C_TO_DELTA_F = 9.0 / 5.0
ALPHA_PER_F_TO_PER_C = 9.0 / 5.0

# ------------------------------------------------------------------- misc ----
G_ACCEL = 9.80665          # m/s^2
KG_TO_TONNE = 1.0e-3
KG_TO_LB = 2.204622622
MPH_TO_MS = 0.44704


def fmt(value, unit="", digits=3):
    """Format a float compactly with a trailing unit label."""
    if value is None:
        return "-"
    av = abs(value)
    if av != 0.0 and (av < 1e-3 or av >= 1e6):
        s = f"{value:.{digits}e}"
    else:
        s = f"{value:,.{digits}f}"
    return f"{s} {unit}".strip()
