"""Scheme A (user-selected, 2026-09-24) -- the single colour source for the ED/SF beautify pass.

Values copied verbatim from 00_audit/palettes/palettes.tsv, column A_harmonised.
Import:  sys.path.insert(0, <beautify>/common); from palA import *
"""
import csv
from pathlib import Path

_T = Path("${WORK_DIR}/fix/beautify/00_audit/palettes/palettes.tsv")
A = {}
with open(_T) as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        A[(r["group"], r["key"])] = r["A_harmonised"]

# materials
FL, TN, TK, NS = A["material", "FL"], A["material", "TN"], A["material", "TK"], A["material", "NS"]
NIG, EOL, EG11, EO12 = A["material", "Nigerian"], A["material", "E_oleifera"], A["material", "EG11"], A["material", "EO12"]
# populations
POP = [A["population", f"Pop{i}"] for i in range(1, 5)]
# SV types
SV = {k: A["sv", k] for k in ("INS", "DEL", "DUP", "INV", "TRA", "COMPLEX", "SYN")}
# ASE classes
ASE = {k: A["ase", k] for k in ("NoDiff", "HapDom", "Sub", "NoASE")}
# direction (A-biased = warm, B-biased = cool; ED3d convention)
EXPAND, CONTRACT, SWITCH = A["direction", "Expansion"], A["direction", "Contraction"], A["direction", "Switch"]
# GWAS
SNP, SVG, THR = A["gwas", "SNP"], A["gwas", "SV_gwas"], A["gwas", "Threshold"]
# colour maps
SEQ = A["colormap", "sequential"].split()
SEQ_WARM = A["colormap", "sequential_warm"].split()
DIV = A["colormap", "diverging"].split()
FL_TN = A["colormap", "FL_TN"].split()

# softened WGCNA module colours (keep colour *names*, lower saturation; issues.tsv SF3 row)
WGCNA_SOFT = {
    "blue": "#3B5BA9", "magenta": "#C04FA0", "greenyellow": "#A8D46F", "turquoise": "#3FB6B0",
    "brown": "#9C6B3F", "yellow": "#E8C94A", "green": "#5DA85A", "red": "#D0504A", "black": "#3A3A3A",
    "pink": "#E79AB8", "purple": "#8E6BB8", "tan": "#CDB48C", "salmon": "#EE9A84", "cyan": "#63C5D8",
    "midnightblue": "#2C3E73", "lightcyan": "#CDEBEE", "grey60": "#999999", "lightgreen": "#A6D8A0",
    "lightyellow": "#F4EFC0", "royalblue": "#4A6FC9", "darkred": "#8E2F2F", "darkgreen": "#2F6B45",
    "darkturquoise": "#2A9FA6", "darkgrey": "#6E6E6E", "orange": "#E8963C", "darkorange": "#C8742E",
    "white": "#F2F2F2", "skyblue": "#86BEE0", "saddlebrown": "#8A5A34", "steelblue": "#4F7FA8",
    "paleturquoise": "#B5E2DE", "violet": "#C48ED8", "darkolivegreen": "#5B6B3A", "darkmagenta": "#8E3F86",
    "grey": "#BDBDBD",
}

LETTER_PT = 8          # panel letters: 8 pt Arial Bold lower case (Nature)
MIN_PT = 5             # body-text floor


def tint(hexc, f=0.4):
    """f = fraction of the colour kept (0.4 -> 40 % tint on white)."""
    h = hexc.lstrip("#")
    r, g, b = (int(h[i:i + 2], 16) for i in (0, 2, 4))
    return "#%02X%02X%02X" % tuple(round(255 - (255 - v) * f) for v in (r, g, b))


def sci_tex(p, digits=2):
    """'3.8 \\times 10^{-143}' for matplotlib mathtext (P-value style)."""
    m, e = f"{p:.{digits - 1}e}".split("e")
    return rf"{m} \times 10^{{{int(e)}}}"
