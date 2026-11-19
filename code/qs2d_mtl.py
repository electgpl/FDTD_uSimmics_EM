# -*- coding: utf-8 -*-
"""
qs2d_mtl.py  -  ELECTGPL
2D quasi-static cross-section solver (Z0, Z0e/Z0o, eps_eff e/o) + MTL S-matrix of the coupled
section, in seconds. Writes a .s4p that kicad_em.py overlays on the openEMS results.
Ports: [line1@start, line1@end, line2@start, line2@end] = [input, through, coupled, isolated].
Only the uniform coupled section is modelled (bends, chamfers and feed lines are not).
Edit the CONFIG block and run (Thonny, no command-line arguments).
"""

# ============================ CONFIG ======================================
OUT_S4P     = "examples/coupler/em_out/mtl_qs2d.s4p"
GRID_MM     = 0.05            # FD cell (0.05 mm: ~1 % on Z0, a few seconds)
WIDTH_MM    = 40.0            # cross-section width = side walls (PEC box on the board edge)
SIDE_WALLS  = True
H_MM        = 1.6             # substrate
ER, TAND    = 4.4, 0.02
AIR_MM      = 4.8             # air between top copper and lid
LID_MM      = 1.6             # lid dielectric (0 = none)
LID_ER, LID_TAND = 4.4, 0.02
LID_METAL   = True            # PEC on the lid's outer face (False = open top)
STRIPS_MM   = [(18.5, 21.5), (24.5, 27.5)]   # (y0, y1) of line 1 (main) and line 2 (coupled)
COUPLED_LEN_MM = 38.0         # straight coupled section length
F_START, F_STOP, N_POINTS = 0.3e9, 3.0e9, 541
BAND        = (902e6, 928e6)
VARIANTS    = True            # also print lid -> air and no-lid variants
# ==========================================================================

import os
import sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__)) if "__file__" in globals() else os.getcwd()
sys.path.insert(0, HERE)
from qemviz import qs2d, touchstone

mm = 1e-3


def cross_section(lid_mm=LID_MM, lid_er=LID_ER, lid_tand=LID_TAND, metal=LID_METAL):
    layers = [(0, H_MM * mm, ER, TAND)]
    top = (H_MM + AIR_MM + lid_mm) * mm
    if lid_mm > 0:
        layers.append(((H_MM + AIR_MM) * mm, top, lid_er, lid_tand))
    return dict(width=WIDTH_MM * mm, side_walls=SIDE_WALLS, layers=layers,
                top=top if metal else None, top_air=None if metal else 30 * mm,
                strips=[(a * mm, b * mm, H_MM * mm) for a, b in STRIPS_MM])


def report(tag, r, f, S):
    k = int(np.argmin(np.abs(f - 915e6)))
    C = -20 * np.log10(abs(S[k, 2, 0])); I = -20 * np.log10(abs(S[k, 3, 0]))
    kmax = int(np.argmax(np.abs(S[:, 2, 0])))
    print(f"{tag:22s} Z0={r['Z0']:6.2f}  Z0e={r['Z0e']:6.2f}  Z0o={r['Z0o']:6.2f}  "
          f"ee={r['ee']:.3f}  eo={r['eo']:.3f} | 915 MHz: C={C:5.2f} I={I:5.2f} D={I-C:5.2f} dB | "
          f"Cmax {-20*np.log10(abs(S[kmax,2,0])):.2f} dB @ {f[kmax]/1e9:.3f} GHz")


f = np.linspace(F_START, F_STOP, N_POINTS)
r = qs2d.solve(cross_section(), d=GRID_MM * mm)
S = qs2d.mtl_sparams(r["L"], r["C"], r["Gw"], COUPLED_LEN_MM * mm, f)
report("nominal (lid FR4)", r, f, S)
out = OUT_S4P if os.path.isabs(OUT_S4P) else os.path.join(HERE, OUT_S4P)
os.makedirs(os.path.dirname(out), exist_ok=True)
touchstone.write(out, f, S, comment="2D quasi-static + MTL, coupled section only")
print("Touchstone:", out)
if VARIANTS:
    for tag, kw in (("lid FR4 -> air", dict(lid_er=1.0, lid_tand=0.0)),
                    ("no lid (open top)", dict(lid_mm=0, metal=False))):
        rv = qs2d.solve(cross_section(**kw), d=GRID_MM * mm)
        report(tag, rv, f, qs2d.mtl_sparams(rv["L"], rv["C"], rv["Gw"], COUPLED_LEN_MM * mm, f))
