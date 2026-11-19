# -*- coding: utf-8 -*-
"""
kicad_em.py  -  ELECTGPL
KiCad 9 .kicad_pcb -> openEMS FDTD: N-port S-matrix, Touchstone, coupler metrics, |E| maps and
phase animations. Materials, ports and enclosure come from a TOML sidecar next to the board.
Edit the CONFIG block and run (Thonny, no command-line arguments).
"""

# ============================ CONFIG ======================================
KICAD_PCB       = "examples/coupler/DirectionalCoupler.kicad_pcb"
SIDECAR         = "examples/coupler/DirectionalCoupler.toml"   # TOML (stdlib tomllib)
OUT_DIR         = "examples/coupler/em_out"                     # relative to this script
OE_DLL_DIR      = None          # Windows: folder with CSXCAD.dll if it is not found automatically
SYMMETRY        = "mirror_x"    # None | "mirror_x" | "mirror_y"  (fewer runs; ports must match)
FORCE_RERUN     = False
MESH = dict(
    res_mm=0.5,                 # fine cell at copper edges (thirds rule)
    lambda_div=20,              # max cell = lambda/20 at f_stop in the densest dielectric
    max_board_mm=1.0,           # max cell in x/y (keeps >= 3 cells across gaps of coupled lines)
    n_sub=4,                    # cells across each substrate layer
    n_lid=2,                    # cells across the lid dielectric
    max_z_mm=1.0,
    end_db=-40.0,
    max_ts=1_000_000,
    n_points=541,
    verbose=1,
)
DUMP_Z_FRAC     = 0.25          # dump plane depth below the top copper, as a fraction of h
                                #   (0.25 = middle of the upper half of the substrate)
DUMP_FREQS      = [915e6, 1830e6]
DUMP_PORT       = 1             # excited port for the field maps
VTK_TIME_DUMP   = False         # time-domain VTK of the same plane (ParaView animation of the pulse)
MAKE_MP4        = True          # 1920x1080, 24 fps, H.264 crf 16 (needs ffmpeg or imageio-ffmpeg)
MAKE_GIF        = True
ANIM_FRAMES     = 72
DB_RANGE        = 50.0          # dynamic range of the |E| maps (coupled line ~ -30 dB)
PROFILES        = [("coupled line, y = 108", (121.0, 108.0), (159.0, 108.0))]  # |E| cuts, KiCad mm
OVERLAY_S4P     = ["examples/coupler/em_out/mtl_qs2d.s4p"]   # extra curves (MTL, VNA .sNp ...)
SHOW            = False         # open the figures at the end
# ==========================================================================

import os
import sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__)) if "__file__" in globals() else os.getcwd()
sys.path.insert(0, HERE)
import matplotlib
if not SHOW:
    matplotlib.use("Agg")
import matplotlib.pyplot as plt
from qemviz import sidecar, fe_kicad, plots, touchstone, fields, oe_backend


def P(p):
    return p if os.path.isabs(p) else os.path.join(HERE, p)


out = P(OUT_DIR)
os.makedirs(out, exist_ok=True)
cfg = sidecar.load(P(SIDECAR))
model = fe_kicad.load(P(KICAD_PCB), cfg)
print(model.summary())
plots.layout(model, os.path.join(out, "01_layout.png"))

names = [p.name for p in model.ports]
h = max(d.z1 - d.z0 for d in model.dielectrics)
zd = max(model.copper_z.values()) - DUMP_Z_FRAC * h
dumps = dict(port=DUMP_PORT, z=zd, freqs=DUMP_FREQS, vtk_time=VTK_TIME_DUMP)
band = cfg["sweep"].get("band")
mcfg = dict(MESH, n_points=cfg["sweep"].get("n_points", MESH["n_points"]),
            f_loss=0.5 * (band[0] + band[1]) if band else None)
f, S, key = oe_backend.run_nport(model, mcfg, out, dll_dir=OE_DLL_DIR, symmetry=SYMMETRY,
                                 dumps=dumps, force=FORCE_RERUN)
kb = (f >= band[0]) & (f <= band[1]) if band else slice(None)
print(f"reciprocity error in band: max |Sij - Sji| = {oe_backend.reciprocity_error(S[kb]):.2e}")

N = len(names)
ts = os.path.join(out, f"{model.name}.s{N}p")
touchstone.write(ts, f, S, comment=f"{model.name} openEMS FDTD, ports {names}")
print("Touchstone:", ts)
plots.sparams(f, S, names, os.path.join(out, "02_sparams.png"), band, f"{model.name}: |S| (openEMS)")

roles = cfg.get("roles", {})
if all(k in roles for k in ("input", "through", "coupled", "isolated")):
    r = {k: names.index(roles[k]) for k in ("input", "through", "coupled", "isolated")}
    overlays = []
    for pth in OVERLAY_S4P:
        if os.path.exists(P(pth)):
            fo, So, _ = touchstone.read(P(pth))
            overlays.append((os.path.splitext(os.path.basename(pth))[0], fo, So))
    plots.coupler_metrics(f, S, r, overlays, os.path.join(out, "03_coupler_metrics.png"), band,
                          f"{model.name}: coupling, isolation, directivity, RL, IL")

    def at(fq, SS, ff):
        k = int(np.argmin(np.abs(ff - fq)))
        i, t, c, iso = r["input"], r["through"], r["coupled"], r["isolated"]
        C = -plots.db(SS[k, c, i]); I = -plots.db(SS[k, iso, i])
        return C, I, I - C, -plots.db(SS[k, i, i]), -plots.db(SS[k, t, i])
    print("\n  f [MHz]   source        C [dB]  I [dB]  D [dB]  RL [dB]  IL [dB]")
    fpts = sorted(set([band[0], 0.5 * (band[0] + band[1]), band[1]] if band else []) | {915e6, 1830e6})
    for fq in fpts:
        for lbl, ff, SS in [("openEMS", f, S)] + overlays:
            print(f"  {fq/1e6:7.1f}   {lbl:12s} " + " ".join(f"{v:7.2f}" for v in at(fq, SS, ff)))

# ---------------- field maps and animations ----------------
sp = os.path.join(out, f"oe_{key}_exc{DUMP_PORT}", "E_plane.h5")
if os.path.exists(sp):
    for fq in DUMP_FREQS:
        x, y, z, E, fr = fields.read_fd(sp, fq)
        E2 = E[..., 0]
        tag = f"{fr/1e6:.0f}MHz"
        fields.top_view(model, x, y, E2, fr, os.path.join(out, f"04_E_top_{tag}.png"), DB_RANGE,
                        title=f"|E| at z = {zd:+.2f} mm, {fr/1e6:.1f} MHz, port {DUMP_PORT} excited (openEMS)")
        if PROFILES:
            fig, ax = plt.subplots(figsize=(9, 4.5))
            mag = np.sqrt(np.sum(np.abs(E2) ** 2, axis=0)); ref = mag.max()
            for lbl, (xa, ya), (xb, yb) in PROFILES:
                t = np.linspace(0, 1, 400)
                px, py = xa + (xb - xa) * t, model.ky(ya + (yb - ya) * t)
                ix = np.abs(x[:, None] - px).argmin(0); iy = np.abs(y[:, None] - py).argmin(0)
                d = np.hypot(px - xa, model.ky(py) - ya)
                ax.plot(xa + (xb - xa) * t if xb != xa else d, 20 * np.log10(mag[ix, iy] / ref + 1e-15), label=lbl)
            ax.set_xlabel("x [mm] (KiCad)"); ax.set_ylabel("|E| [dB rel. max]"); ax.grid(True, alpha=0.4)
            ax.legend(); ax.set_title(f"|E| profile at z = {zd:+.2f} mm, {fr/1e6:.1f} MHz", fontsize=10)
            fig.savefig(os.path.join(out, f"06_E_profile_{tag}.png"), dpi=150, bbox_inches="tight"); plt.close(fig)
        fields.phase_animation(model, x, y, E2, fr, os.path.join(out, f"05_E_phase_{tag}"),
                               frames=ANIM_FRAMES, db_range=DB_RANGE, gif=MAKE_GIF, mp4=MAKE_MP4)
else:
    print("no field dump found:", sp)
print("Files in:", out)
if SHOW:
    plt.show()
