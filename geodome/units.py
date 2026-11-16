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
