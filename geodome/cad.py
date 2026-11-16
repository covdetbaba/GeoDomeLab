"""
CAD export: the dome geometry as an ASCII DXF file (AutoCAD Release 12).

A DXF file is a plain list of *group pairs*: a group code (an integer on its
own line) followed by a value on the next line.  The code says what the value
is (Autodesk, "DXF Reference"):

    0   start of an entity, table entry or section marker
    2   name (section, table, layer, linetype)
    8   layer of an entity
    9   name of a header variable
    6   linetype name of a layer
    62  colour number (AutoCAD Colour Index, ACI 1..255)
    70  integer flags / counts
    10, 20, 30     first point  x, y, z
    11, 21, 31     second point x, y, z   (and 12.., 13.. for 3DFACE)

The file written here has four sections::

    HEADER    $ACADVER AC1009 (R12), $INSBASE, $EXTMIN / $EXTMAX (the model
              box), $INSUNITS 6 (= metres), $PDMODE / $PDSIZE (hub symbol)
    TABLES    LTYPE (CONTINUOUS) and LAYER: "0", one STRUT_<family> layer per
              strut family, PANELS and HUBS, each with an ACI colour
    BLOCKS    empty
    ENTITIES  one LINE per strut       on layer STRUT_<family>
              one 3DFACE per panel     on layer PANELS
              one POINT per node (hub) on layer HUBS

A 3DFACE always has four corners; a triangle repeats its third corner as the
fourth (the DXF convention for a three-sided face).  Entities are written in
model order: LINE k is member k, 3DFACE k is face k and POINT k is node k, so
a CAD user can match an entity to the analysis tables.

Coordinates are in metres, the program's internal unit, written in fixed
notation with 12 decimals (1 pm), so they read back exactly to 1e-12 m.  The
file is pure ASCII with CRLF line ends and needs no third-party package.

Colours.  The strut layers use the colour codes of the fabrication cut list
(``fabrication.STRUT_COLOURS``: family A red, B blue, C green, ...), mapped to
the nearest ACI colour.  When analysis results are given, every LINE also
gets its own colour (group 62) from its utilisation band, using the bands
of the Utilisation view (up to 0.25, 0.50, 0.75, 0.90, 1.00, and over 1.00).
The struts stay on their family layers, so a CAD user can switch back to
family colours by setting the entity colour to ByLayer.

Reference: Autodesk, AutoCAD Release 12 DXF format ("DXF Reference",
group codes, HEADER, TABLES and ENTITIES sections; header variable INSUNITS,
value 6 = metres).
"""

from __future__ import annotations

import math

import numpy as np

from . import geometry as geo

__all__ = ["export_dxf", "dxf_text", "layer_table", "strut_layer",
           "utilisation_aci", "FAMILY_ACI", "UTIL_BOUNDS", "UTIL_ACI",
           "PANEL_LAYER", "HUB_LAYER", "PANEL_ACI", "HUB_ACI"]

PANEL_LAYER = "PANELS"
HUB_LAYER = "HUBS"

# ACI numbers for the fabrication colour codes red, blue, green, yellow,
# orange, purple, white, black, brown, grey, pink, cyan (in that order).
FAMILY_ACI = (1, 5, 3, 2, 30, 190, 7, 250, 32, 8, 221, 4)
# Further families (13 and up) cycle through these distinct hues.
_EXTRA_ACI = (10, 50, 90, 130, 170, 210, 20, 60, 100, 140, 180, 220, 240)
PANEL_ACI = 9          # light grey
HUB_ACI = 7            # white on a dark CAD background, black on a light one

# Utilisation bands (as plotting.draw_utilisation) and their ACI colours.
UTIL_BOUNDS = (0.0, 0.25, 0.50, 0.75, 0.90, 1.00)
UTIL_ACI = (4, 3, 2, 30, 1, 6)     # cyan, green, yellow, orange, red, magenta


# ================================================================ helpers ===
def _resolve(model):
    """``(geometry, model_or_None)`` for a DomeModel or a bare DomeGeometry."""
    if isinstance(model, geo.DomeGeometry):
        return model, None
    g = getattr(model, "geom", None)
    if isinstance(g, geo.DomeGeometry):
        return g, model
    raise TypeError("export_dxf needs a DomeModel or a DomeGeometry, "
                    f"got {type(model).__name__}.")


def _family_label(g):
    """Family label as in ``DomeGeometry.strut_schedule`` (A..Z, S26, ...)."""
    return chr(ord("A") + g) if g < 26 else f"S{g}"


def strut_layer(g):
    """Layer name of strut family ``g`` (0-based): STRUT_A, STRUT_B, ..."""
    return f"STRUT_{_family_label(int(g))}"


def _family_aci(g):
    g = int(g)
    if g < len(FAMILY_ACI):
        return FAMILY_ACI[g]
    return _EXTRA_ACI[(g - len(FAMILY_ACI)) % len(_EXTRA_ACI)]


def utilisation_aci(util):
    """ACI colour of each utilisation value (array in, int array out).

    Bands: u <= 0.25, 0.25 < u <= 0.50, 0.50 < u <= 0.75, 0.75 < u <= 0.90,
    0.90 < u <= 1.00 and u > 1.00 (overstressed).  A value on a bound goes
    to the lower band, so u = 1.00 exactly is still "OK", as in the member
    table.
    """
    u = np.asarray(util, dtype=float)
    idx = np.searchsorted(np.asarray(UTIL_BOUNDS[1:]), u, side="left")
    idx = np.clip(idx, 0, len(UTIL_ACI) - 1)
    return np.asarray(UTIL_ACI)[idx]


def _num(v):
    """A coordinate as fixed-point text with up to 12 decimals (no exponent)."""
    v = float(v)
    if not math.isfinite(v):
        raise ValueError(f"Cannot write a non-finite coordinate ({v}) to DXF.")
    s = f"{v:.12f}".rstrip("0")
    if s.endswith("."):
        s += "0"
    if s in ("-0.0", "-0"):
        s = "0.0"
    return s


def _utilisation(results, n_members):
    """Utilisation per member from an AnalysisResults or an array, or None."""
    if results is None:
        return None
    u = getattr(results, "util_env", results)
    if u is None:
        return None
    u = np.asarray(u, dtype=float).ravel()
    if len(u) != n_members:
        raise ValueError(
            f"The results have {len(u)} utilisation values but the model has "
            f"{n_members} struts. Run the analysis again for this model.")
    return u


def layer_table(model):
    """``[(layer name, ACI colour)]`` in the order they are written."""
    g, _ = _resolve(model)
    m = len(g.members)
    nfam = int(np.max(g.strut_group)) + 1 if m else 0
    layers = [("0", 7)]
    layers += [(strut_layer(k), _family_aci(k)) for k in range(nfam)]
    layers += [(PANEL_LAYER, PANEL_ACI), (HUB_LAYER, HUB_ACI)]
    return layers


# ================================================================= writer ===
def dxf_text(model, results=None):
    """The DXF file content as one string (CRLF line ends added on write).

    ``model`` is a DomeModel or a DomeGeometry.  ``results`` is optional: an
    AnalysisResults (its ``util_env`` is used) or an array with one
    utilisation per strut; then every LINE gets a colour by utilisation band.
    """
    g, _ = _resolve(model)
    P = np.asarray(g.nodes, dtype=float).reshape(-1, 3)
    M = np.asarray(g.members, dtype=int).reshape(-1, 2)
    F = np.asarray(g.faces, dtype=int).reshape(-1, 3) if len(g.faces) else \
        np.zeros((0, 3), dtype=int)
    fam = np.asarray(g.strut_group, dtype=int).ravel()
    if len(fam) != len(M):
        raise ValueError("strut_group must have one entry per strut.")
    if len(P) and not np.all(np.isfinite(P)):
        raise ValueError("The model has non-finite node coordinates.")
    util = _utilisation(results, len(M))
    colours = utilisation_aci(util) if util is not None else None

    out = []
    put = out.append

    def pair(code, value):
        put(f"{int(code):>3d}")
        put(str(value))

    def point(base, xyz):
        pair(10 + base, _num(xyz[0]))
        pair(20 + base, _num(xyz[1]))
        pair(30 + base, _num(xyz[2]))

    # ------------------------------------------------------------ HEADER --
    lo = P.min(axis=0) if len(P) else np.zeros(3)
    hi = P.max(axis=0) if len(P) else np.zeros(3)
    if len(M):
        Ls = np.linalg.norm(P[M[:, 1]] - P[M[:, 0]], axis=1)
        pdsize = 0.05 * float(Ls.min())
    else:
        pdsize = 0.05 * float(np.max(hi - lo)) if len(P) else 0.1
    if not (pdsize > 0.0):
        pdsize = 0.1
    pair(0, "SECTION")
    pair(2, "HEADER")
    pair(9, "$ACADVER")
    pair(1, "AC1009")
    pair(9, "$INSBASE")
    point(0, (0.0, 0.0, 0.0))
    pair(9, "$EXTMIN")
    point(0, lo)
    pair(9, "$EXTMAX")
    point(0, hi)
    pair(9, "$INSUNITS")
    pair(70, 6)
    pair(9, "$PDMODE")
    pair(70, 34)
    pair(9, "$PDSIZE")
    pair(40, _num(pdsize))
    pair(0, "ENDSEC")

    # ------------------------------------------------------------ TABLES --
    layers = layer_table(g)
    pair(0, "SECTION")
    pair(2, "TABLES")
    pair(0, "TABLE")
    pair(2, "LTYPE")
    pair(70, 1)
    pair(0, "LTYPE")
    pair(2, "CONTINUOUS")
    pair(70, 0)
    pair(3, "Solid line")
    pair(72, 65)
    pair(73, 0)
    pair(40, "0.0")
    pair(0, "ENDTAB")
    pair(0, "TABLE")
    pair(2, "LAYER")
    pair(70, len(layers))
    for name, aci in layers:
        pair(0, "LAYER")
        pair(2, name)
        pair(70, 0)
        pair(62, int(aci))
        pair(6, "CONTINUOUS")
    pair(0, "ENDTAB")
    pair(0, "ENDSEC")

    # ------------------------------------------------------------ BLOCKS --
    pair(0, "SECTION")
    pair(2, "BLOCKS")
    pair(0, "ENDSEC")

    # ---------------------------------------------------------- ENTITIES --
    pair(0, "SECTION")
    pair(2, "ENTITIES")
    for k, (i, j) in enumerate(M):
        pair(0, "LINE")
        pair(8, strut_layer(fam[k]))
        if colours is not None:
            pair(62, int(colours[k]))
        point(0, P[i])
        point(1, P[j])
    for a, b, c in F:
        pair(0, "3DFACE")
        pair(8, PANEL_LAYER)
        point(0, P[a])
        point(1, P[b])
        point(2, P[c])
        point(3, P[c])          # triangle: the fourth corner repeats the third
    for p in P:
        pair(0, "POINT")
        pair(8, HUB_LAYER)
        point(0, p)
    pair(0, "ENDSEC")
    pair(0, "EOF")
    return "\n".join(out) + "\n"


def export_dxf(model, path, results=None):
    """Write the geometry of ``model`` to an ASCII DXF R12 file at ``path``.

    One LINE per strut (layer ``STRUT_<family>``), one 3DFACE per panel
    (layer ``PANELS``), one POINT per node (layer ``HUBS``); units metres
    (``$INSUNITS`` 6).  With ``results`` (an AnalysisResults or one
    utilisation value per strut) each strut is coloured by its utilisation
    band.  A custom truss without faces simply has no 3DFACE entities.
    Returns ``path``.
    """
    text = dxf_text(model, results=results)
    with open(path, "w", encoding="ascii", newline="\r\n") as fh:
        fh.write(text)
    return path
