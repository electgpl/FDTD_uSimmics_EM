# -*- coding: utf-8 -*-
"""
qucs_em_viz.py  -  ELECTGPL
UNIVERSAL visualizer for uSimmics / QucsStudio EM layout projects (.pcb + .dat).
Independent of the antenna type: patch arrays, MIFA, IFA, dipoles, filters, lines, etc.

Always (any design):
  01_layout.png      actual copper from the .pcb, ports and shorts to ground, overall size
  02_s11.png         |Sii| in dB, shaded ISM bands, -10 dB bandwidth (+ CSV export)
  03_smith.png       Sii on a Smith chart
If a radiation pattern is available (see PATTERN_SOURCE):
  04_cuts_polar.png  cuts at phi = 0 and phi = 90 deg
  05_pattern_3d.png  3D pattern over the PCB (interactive window: rotate with the mouse)
  06_field_map_db.png  |E|*r on a vertical plane
  wave.gif (+ live animated wave window); pattern_3d_rot.gif only if MAKE_ROT_GIF

Pattern sources (PATTERN_SOURCE = "auto" tries them in this order: dat -> openems -> patch):
  "dat"   : far-field computed by uSimmics (enable "Radiation" in the layout and simulate).
            The right option for ANY antenna. The field on the plane is synthesized from
            the far-field, so it is only valid outside the near-field region
            (r < max(lambda/2pi, antenna size) is masked).
  "openems": our own FDTD simulation with openEMS built from the .pcb primitives
            (uSimmics does not need to export anything). It also produces the REAL near field:
            07 |E| over the board, 08 current density |J| on the copper, 09 vertical cut,
            nearfield.gif, nearfield_side.gif. Results are cached (no re-simulation unless
            something changes).
  "patch" : analytical model of rectangular patches (cavity model + array factor), only if
            PATCH_REGEX finds patches. No mutual coupling, no feed network.
  "none"  : layout and S-parameters only.
Everything is set in the CONFIG block (meant for Thonny, no command-line arguments).
"""

# ============================ CONFIG ======================================
QUCS_DIR        = "examples"                 # project folder (absolute or relative to the working dir)
PROJECT         = "Radar24G_Array1x4_Layout" # base name (.pcb / .dat)
PATTERN_SOURCE  = "auto"        # "auto" | "dat" | "openems" | "patch" | "none"
F_PATTERN_HZ    = None          # None: .pcb Radiation frequency if inside the sweep; otherwise S11 minimum
CUT_PHI_DEG     = None          # field-map plane; None = auto (narrowest-beam plane)
PATCH_REGEX     = r"^(P\d+)_"   # "patch" mode only
POL_AXIS        = "auto"        # "patch" mode: "x", "y" or "auto"
ELEM_PHASE_DEG  = None          # "patch" mode: per-element phases (None = in phase)
ELEM_AMPL       = None          # "patch" mode: amplitudes (None = uniform)
DB_RANGE        = 30.0
FIELD_HALF_LAM  = 5.0           # field map: +-5 lambda0 (adapts to any frequency)
FIELD_RES_LAM   = 1 / 30        # field-map resolution in lambda0
MAKE_GIFS       = True
MAKE_ROT_GIF    = False         # rotating 3D GIF (the 3D window can be rotated with the mouse anyway)
N_FRAMES        = 24
OUT_DIR         = "qucs_viz"    # relative to QUCS_DIR
SHOW            = True          # False to run in batch without windows
SHOW_WAVE_LIVE  = True
CMAP            = "turbo"
# --- backend openEMS (FDTD) ---
OE_ENABLE       = True          # in "auto" mode, use openEMS if the .dat has no far-field
OE_DLL_DIR      = None          # Windows: openEMS folder if Python cannot find the DLLs (e.g. r"C:\openEMS")
OE_GROUND       = "auto"        # bottom-side ground: "auto" | "none" | "full" | "largest" | [(x0,y0,x1,y1) in m]
                                #   auto: "full" if there are patches; otherwise "none" (drawn copper only)
OE_EXTRA_TOP    = []            # extra top copper [(x0,y0,x1,y1) in m], e.g. the rest of the board ground
OE_BOARD_MARGIN_MM = 1.0        # substrate beyond the copper
OE_RES_MM       = None          # fine cell on the board; None = auto (half of the narrow traces, >= 2 cells per port)
OE_END_DB       = -40.0         # energy end criterion
OE_MAX_TS       = 400000
OE_AIR_LAM      = 0.25          # air around, in lambda at the lowest frequency
OE_AIR_UP_LAM   = 1.0           # air ABOVE the board (lambda): sets how much of the outgoing wave is visible
                                #   in 09 / nearfield_side.gif. Larger = nicer and slower
OE_NEAR_DZ_MM   = 0.5           # height of the |E| plane above the copper
OE_N_FDUMP      = 21            # frequencies of the field dumps inside the sweep
OE_FORCE_RERUN  = False
BANDS_HZ        = {             # bands inside the sweep are shaded
    "ISM 433": (433.05e6, 434.79e6), "SRD 868": (863e6, 870e6), "ISM 915": (902e6, 928e6),
    "ISM 2.4": (2.400e9, 2.4835e9), "UNII 5": (5.150e9, 5.850e9), "ISM 24": (24.00e9, 24.25e9),
}
# ==========================================================================

import os
import re
import struct
import numpy as np
import matplotlib
if not SHOW:
    matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import animation, cm, colors
from matplotlib.patches import Rectangle, Circle
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

C0 = 299792458.0
UNITS = {"": 1.0, "m": 1.0, "mm": 1e-3, "um": 1e-6, "µm": 1e-6, "mil": 25.4e-6,
         "Hz": 1.0, "kHz": 1e3, "MHz": 1e6, "GHz": 1e9}


def to_si(txt):
    m = re.match(r"\s*([-+0-9.eE]+)\s*(\S*)", txt)
    return float(m.group(1)) * UNITS.get(m.group(2), 1.0)


# ============================ LECTURA .pcb =================================
def parse_pcb(path):
    """Rectangle NAME x y dx|dy|w|h  -> copper in [x+dx, x+dx+w] x [y+dy, y+dy+h]
       Lumped  NAME x y num dir Z ... wx wy -> num > 0: port; num = 0: lumped element (e.g. short to ground)"""
    props, rects, lumped = {}, [], []
    sect = None
    with open(path, encoding="utf-8", errors="replace") as fh:
        for ln in fh:
            s = ln.strip()
            if s in ("<Properties>", "<Items>"):
                sect = s; continue
            if s.startswith("</"):
                sect = None; continue
            p = s.split()
            if sect == "<Properties>" and "=" in s:
                k, v = s.split("=", 1); props[k] = v
            elif sect == "<Items>" and p and p[0] == "Rectangle":
                x, y = float(p[2]), float(p[3])
                dx, dy, w, h = (float(v) for v in p[4].split("|"))
                rects.append((p[1], x + dx, y + dy, x + dx + w, y + dy + h))
            elif sect == "<Items>" and p and p[0] == "Lumped":
                lumped.append({"name": p[1], "x": float(p[2]), "y": float(p[3]),
                               "port": int(float(p[4])) if len(p) > 4 else 0,
                               "z": float(p[6]) if len(p) > 6 else 0.0,
                               # last two fields: port width along x and along y (one of them is 0)
                               "wx": float(p[8]) if len(p) > 8 else 0.0,
                               "wy": float(p[9]) if len(p) > 9 else 0.0})
                lumped[-1]["w"] = max(lumped[-1]["wx"], lumped[-1]["wy"])
            elif sect == "<Items>" and p and p[0] not in ("Rectangle", "Lumped"):
                print(f"WARNING: primitive '{p[0]}' not supported yet (ignored in the drawing)")
    sub = props.get("Substrate", "").split("|")
    substrate = {"er": float(sub[1]), "h": to_si(sub[2]), "tand": float(sub[3] or 0),
                 "t": to_si(sub[5])} if len(sub) > 5 else None
    fr = props.get("Frequency", "").split("|")
    sweep = (to_si(fr[0]), to_si(fr[1]), int(fr[2])) if len(fr) == 3 else None
    rad = props.get("Radiation", "").split("|")
    radiation = {"enabled": rad[0].strip() == "1", "f": to_si(rad[1])} if len(rad) > 1 else None
    return props, rects, lumped, substrate, sweep, radiation


# ============================ LECTURA .dat (QucsData binario) ==============
def parse_dat(path):
    """uSimmics/QucsStudio binary QucsData dataset.
    'QucsData' + 4 bytes + int32 header length; header: [int32 type][int32 N][name\0]
    [dependencies\0] (dependencies may carry matrix dimensions: 'frequency 1 1').
    Block sizes are NOT stored explicitly: the type code is used (3 = complex dependent,
    7 = real dependent, 5 = real independent) and, for unknown codes, the possible
    combinations (element count x real/complex) are tried until the sum matches the file size."""
    import itertools
    d = open(path, "rb").read()
    if not d.startswith(b"QucsData"):
        raise ValueError("Not a binary QucsData dataset")
    hlen = struct.unpack("<i", d[12:16])[0]
    p, end, ents = 16, 16 + hlen, []
    while p < end:
        typ, n = struct.unpack("<ii", d[p:p + 8]); p += 8
        strs = []
        while p < end and d[p] >= 0x20:
            q = d.index(b"\0", p); strs.append(d[p:q].decode("latin-1")); p = q + 1
        toks = " ".join(strs[1:]).split()
        deps = [t for t in toks if not re.fullmatch(r"\d+", t)]
        dims = [int(t) for t in toks if re.fullmatch(r"\d+", t)]
        ents.append({"type": typ, "n": n, "name": strs[0] if strs else "?", "deps": deps, "dims": dims})
    remaining = len(d) - end

    def header_table():
        return "\n".join(f"   type={e['type']:3d}  N={e['n']:7d}  {e['name']:<16s} deps={e['deps']} dims={e['dims']}"
                         for e in ents)

    nmap = {e["name"]: e["n"] for e in ents if not e["deps"]}
    cands = []
    for e in ents:
        md = int(np.prod(e["dims"])) if e["dims"] else 1
        cnt = {e["n"], e["n"] * md}
        if e["deps"] and all(x in nmap for x in e["deps"]):
            pd = int(np.prod([nmap[x] for x in e["deps"]]))
            cnt |= {pd, pd * md}
        # observed codes: 3 = complex dependent, 7 = real dependent, 5 = real independent
        sizes = {3: [16], 7: [8], 5: [8]}.get(e["type"], [(16 if e["deps"] else 8), (8 if e["deps"] else 16)])
        cands.append([(c, el) for el in sizes for c in sorted(cnt)])
    sol = None
    for combo in itertools.product(*cands):
        if sum(c * el for c, el in combo) == remaining:
            sol = combo; break
    if sol is None:
        raise ValueError(f"Could not resolve the .dat format ({remaining} data bytes). "
                         f"Header:\n{header_table()}")
    out, p = {}, end
    for e, (c, el) in zip(ents, sol):
        v = np.frombuffer(d[p:p + c * el], "<c16" if el == 16 else "<f8").copy()
        p += c * el
        out[e["name"]] = {"v": v, "deps": e["deps"], "dims": e["dims"]}
    for k, e in out.items():
        lens = [out[x]["v"].size for x in e["deps"] if x in out]
        md = int(np.prod(e["dims"])) if e["dims"] else 1
        try:
            if md > 1:     # matrix per point (order NOT verified for >1 port)
                e["v"] = e["v"].reshape([*lens, *e["dims"]])
            elif len(lens) > 1:
                e["v"] = e["v"].reshape(lens, order="F")
        except ValueError:
            print(f"WARNING .dat: could not reshape '{k}' ({e['v'].size} values, deps {lens})")
    print(".dat header (resolved blocks):\n" + "\n".join(
        f"   type={e['type']:3d}  N={e['n']:7d}  {e['name']:<16s} deps={e['deps']} dims={e['dims']}"
        f"  -> {c} x {'complex' if el == 16 else 'real'}" for e, (c, el) in zip(ents, sol)))
    return out


def get_sii(dat, i=1):
    """Returns (f, Sii) accepting 'S11', 'S[1,1]' or the 'S' matrix."""
    for k in (f"S{i}{i}", f"S[{i},{i}]"):
        if k in dat and dat[k]["deps"]:
            return dat[dat[k]["deps"][0]]["v"], dat[k]["v"].ravel()
    if "S" in dat and dat["S"]["deps"]:
        v = dat["S"]["v"]
        f = dat[dat["S"]["deps"][0]]["v"]
        if v.ndim == 1 and v.size == f.size:
            return f, v
        if v.ndim >= 3:
            return f, v[..., i - 1, i - 1].ravel()
    return None, None


# ============================ UTILIDADES ===================================
def eps_eff(er, h, w):
    return (er + 1) / 2 + (er - 1) / 2 / np.sqrt(1 + 12 * h / w)


def f_res(L, W, er, h):
    ee = eps_eff(er, h, W)
    dL = 0.412 * h * (ee + 0.3) * (W / h + 0.264) / ((ee - 0.258) * (W / h + 0.8))
    return C0 / (2 * (L + 2 * dL) * np.sqrt(ee))


def interp2(xg, yg, V, x, y):
    """Bilinear interpolation on a regular grid (no scipy). Out of range -> 0."""
    x = np.asarray(x, float); y = np.asarray(y, float)
    ix = np.clip(np.searchsorted(xg, x) - 1, 0, len(xg) - 2)
    iy = np.clip(np.searchsorted(yg, y) - 1, 0, len(yg) - 2)
    tx = (x - xg[ix]) / (xg[ix + 1] - xg[ix]); ty = (y - yg[iy]) / (yg[iy + 1] - yg[iy])
    out = (V[ix, iy] * (1 - tx) * (1 - ty) + V[ix + 1, iy] * tx * (1 - ty) +
           V[ix, iy + 1] * (1 - tx) * ty + V[ix + 1, iy + 1] * tx * ty)
    bad = (x < xg[0] - 1e-9) | (x > xg[-1] + 1e-9) | (y < yg[0] - 1e-9) | (y > yg[-1] + 1e-9)
    return np.where(bad, 0.0, out)


def ff_axes(deps):
    """Identifies the angular axes of a far-field variable (theta/phi or elevation/azimuth)."""
    th = next((x for x in deps if re.search(r"theta|elev", x, re.I)), None)
    ph = next((x for x in deps if re.search(r"phi|azim", x, re.I)), None)
    return th, ph


def hpbw(tt, g, circular=False):
    """-3 dB width of the main lobe (the maximum closest to broadside)."""
    n = len(g); gmax = g.max()
    cand = np.where(g >= gmax - 0.01)[0]
    i0 = int(cand[np.argmin(np.abs(tt[cand]))])
    l = r = 0
    while l < n and g[(i0 - l - 1) % n if circular else max(i0 - l - 1, 0)] >= gmax - 3 and \
            (circular or i0 - l - 1 >= 0):
        l += 1
    while r < n and g[(i0 + r + 1) % n if circular else min(i0 + r + 1, n - 1)] >= gmax - 3 and \
            (circular or i0 + r + 1 <= n - 1):
        r += 1
    return np.rad2deg((l + r) * (tt[1] - tt[0]))



# ============================ BACKEND openEMS ==============================
def oe_import():
    """On Windows, Python >= 3.8 does not search the PATH for DLLs: the openEMS folder
    (the one containing openEMS.exe, CSXCAD.dll, openEMS.dll) must be registered with os.add_dll_directory."""
    if os.name == "nt":
        cands = [OE_DLL_DIR, os.environ.get("OPENEMS_INSTALL_PATH"), r"C:\openEMS",
                 os.path.join(os.environ.get("ProgramFiles", r"C:\Program Files"), "openEMS"),
                 os.path.expanduser(r"~\openEMS")]
        for c in cands:
            if c and os.path.isfile(os.path.join(c, "CSXCAD.dll")):
                os.add_dll_directory(c)
                print("openEMS: DLLs from", c)
                break
        else:
            print("openEMS: CSXCAD.dll not found; set OE_DLL_DIR to the openEMS folder")
    from CSXCAD import ContinuousStructure
    from openEMS import openEMS
    from openEMS.physical_constants import EPS0
    return ContinuousStructure, openEMS, EPS0


def run_openems(rects, lumped, sub, f_lo, f_hi, sim_path, run, log=print):
    """Builds the structure from the .pcb primitives (drawing unit: mm) and, if
    run=True, runs the FDTD. Returns ports (S11, Pacc), NF2FF and geometry."""
    ContinuousStructure, openEMS, EPS0 = oe_import()
    U = 1e-3
    er, h, tand = sub["er"], sub["h"] / U, sub.get("tand", 0.0)
    f_c = 0.5 * (f_lo + f_hi)
    R = [(n, x0 / U, y0 / U, x1 / U, y1 / U) for n, x0, y0, x1, y1 in rects]
    extra = [tuple(np.asarray(b, float) / U) for b in OE_EXTRA_TOP]
    allr = [r[1:] for r in R] + extra
    exs = np.array([[r[0], r[2]] for r in allr]).ravel(); eys = np.array([[r[1], r[3]] for r in allr]).ravel()
    m = OE_BOARD_MARGIN_MM
    bx0, bx1, by0, by1 = exs.min() - m, exs.max() + m, eys.min() - m, eys.max() + m

    FDTD = openEMS(NrTS=OE_MAX_TS, EndCriteria=10 ** (OE_END_DB / 10))
    FDTD.SetGaussExcite(f_c, max(0.6 * (f_hi - f_lo), 0.05 * f_c))
    FDTD.SetBoundaryCond(["MUR"] * 6)
    CSX = ContinuousStructure(); FDTD.SetCSX(CSX)
    mesh = CSX.GetGrid(); mesh.SetDeltaUnit(U)
    CSX.AddMaterial("substrate", epsilon=er, kappa=2 * np.pi * f_c * EPS0 * er * tand).AddBox(
        [bx0, by0, -h], [bx1, by1, 0], priority=0)
    cu = CSX.AddMetal("copper")
    for r in allr:
        cu.AddBox([r[0], r[1], 0], [r[2], r[3], 0], priority=10)

    gm = OE_GROUND
    if gm == "auto":
        gm = "full" if any(re.match(PATCH_REGEX, r[0]) for r in R) else "none"
    gnd = []
    if gm == "full":
        gnd = [(bx0, by0, bx1, by1)]
    elif gm == "largest":
        gnd = [max(allr, key=lambda r: (r[2] - r[0]) * (r[3] - r[1]))]
    elif isinstance(gm, (list, tuple)):
        gnd = [tuple(np.asarray(b, float) / U) for b in gm]
    for b in gnd:
        cu.AddBox([b[0], b[1], -h], [b[2], b[3], -h], priority=10)
    log(f"openEMS: bottom ground = {gm!r}")

    def inside(b, x, y):
        return b[0] - 1e-6 <= x <= b[2] + 1e-6 and b[1] - 1e-6 <= y <= b[3] + 1e-6

    lx, ly, ports = [], [], []
    for L in lumped:
        x, y, w = L["x"] / U, L["y"] / U, max(L["w"] / U, 0.1)
        over = any(inside(b, x, y) for b in gnd)
        if L["port"] == 0:
            if over:
                cu.AddBox([x - w / 2, y - w / 2, -h], [x + w / 2, y + w / 2, 0], priority=10)
                lx += [x - w / 2, x + w / 2]; ly += [y - w / 2, y + w / 2]
                log(f"openEMS: '{L['name']}' -> via to ground")
            else:
                log(f"openEMS: '{L['name']}' ({L['z']:.0f} Ω) is not over a bottom ground: ignored")
            continue
        if over:                                          # vertical port: trace -> bottom ground
            if L["wx"] > 0:
                start, stop = [x - w / 2, y, -h], [x + w / 2, y, 0]
            else:
                start, stop = [x, y - w / 2, -h], [x, y + w / 2, 0]
            d = "z"
            log(f"openEMS: port {L['name']} vertical (trace -> bottom ground)")
        else:                                             # in-plane port (gap to the nearest copper)
            best = None
            for r in allr:
                if inside(r, x, y):
                    continue
                if r[0] - 1e-6 <= x <= r[2] + 1e-6:
                    for e in (r[1], r[3]):
                        if best is None or abs(e - y) < best[0]:
                            best = (abs(e - y), "y", e)
                if r[1] - 1e-6 <= y <= r[3] + 1e-6:
                    for e in (r[0], r[2]):
                        if best is None or abs(e - x) < best[0]:
                            best = (abs(e - x), "x", e)
            if best is None:
                raise RuntimeError(f"no reference found for port {L['name']}")
            if best[1] == "y":
                start, stop = [x - w / 2, min(y, best[2]), 0], [x + w / 2, max(y, best[2]), 0]
            else:
                start, stop = [min(x, best[2]), y - w / 2, 0], [max(x, best[2]), y + w / 2, 0]
            d = best[1]
            log(f"openEMS: port {L['name']} in-plane, gap {best[0]:.3f} mm along {d}")
        lx += [start[0], stop[0]]; ly += [start[1], stop[1]]
        ports.append((L["name"], FDTD.AddLumpedPort(L["port"], L["z"] or 50, start, stop, d,
                                                     excite=1.0, priority=5)))

    lam_min = C0 / f_hi / U
    lam_sub = lam_min / np.sqrt(er)
    widths = [min(r[2] - r[0], r[3] - r[1]) for r in allr if min(r[2] - r[0], r[3] - r[1]) > 0]
    res = OE_RES_MM or max(np.percentile(widths, 10) / 2, lam_sub / 200)   # half of the narrow traces
    pw = [max(L["w"] / U, 0.1) for L in lumped if L["port"] > 0]
    if pw and res > min(pw) / 2:      # at least 2 cells across every port, otherwise it may not connect
        log(f"openEMS: fine cell reduced {res:.3f} -> {min(pw)/2:.3f} mm (port width)")
        res = min(pw) / 2
    air = OE_AIR_LAM * C0 / f_lo / U
    air_up = max(air, OE_AIR_UP_LAM * C0 / f_lo / U)
    def merge(v, tol):
        """Merges nearly coincident lines (edges differing by microns create tiny
        cells and an absurd timestep)."""
        v = np.sort(np.asarray(v, float)); out = [v[0]]
        for a in v[1:]:
            if a - out[-1] > tol:
                out.append(a)
        return out

    tol = 0.25 * res
    mesh.AddLine("x", merge(list(exs) + lx + [bx0, bx1, bx0 - air, bx1 + air], tol))
    mesh.AddLine("y", merge(list(eys) + ly + [by0, by1, by0 - air, by1 + air], tol))
    nz = max(3, int(np.ceil(h / res)))
    mesh.AddLine("z", list(np.linspace(-h, 0, nz + 1)) + [-h - air, air_up])
    for a in "xyz":
        mesh.SmoothMeshLines(a, lam_min / 20, 1.4)
    for a, a0, a1 in (("x", bx0, bx1), ("y", by0, by1)):
        ln = np.array(mesh.GetLines(a)); ln = ln[(ln >= a0) & (ln <= a1)]
        new = []
        for u0, u1 in zip(ln[:-1], ln[1:]):
            k = int(np.ceil((u1 - u0) / res))
            new += list(np.linspace(u0, u1, k + 1)[1:-1]) if k > 1 else []
        if new:
            mesh.AddLine(a, new)
    # final cleanup of tiny cells, NEVER dropping copper/port edges (dropping one disconnects the port)
    keep = {"x": np.array(sorted(set(merge(list(exs) + lx + [bx0, bx1], tol)))),
            "y": np.array(sorted(set(merge(list(eys) + ly + [by0, by1], tol)))),
            "z": np.array([-h, 0.0])}
    for a in "xyz":
        tol_f = 0.2 * min(res, h / 3)
        out = []
        for v in np.sort(np.asarray(mesh.GetLines(a), float)):
            is_k = np.min(np.abs(keep[a] - v)) < 1e-9
            near_k = np.min(np.abs(keep[a] - v)) < tol_f
            if near_k and not is_k:
                continue                                  # auxiliary line glued to an edge: drop it
            if out and v - out[-1][0] < tol_f:
                if is_k and not out[-1][1]:
                    out[-1] = (v, True)
                continue
            out.append((v, is_k))
        mesh.ClearLines(a); mesh.AddLine(a, [v for v, _ in out])

    dmin = {a: np.diff(np.asarray(mesh.GetLines(a))).min() for a in "xyz"}
    log("openEMS: smallest cell " + ", ".join(f"{a}={v*1e3:.1f} um" for a, v in dmin.items()))
    nx, ny, nzz = (len(mesh.GetLines(a)) for a in "xyz")
    log(f"openEMS: fine cell {res:.3f} mm, mesh {nx} x {ny} x {nzz} = {nx*ny*nzz/1e6:.2f} Mcells")

    nf2ff = FDTD.CreateNF2FFBox()
    fd = list(np.linspace(f_lo, f_hi, OE_N_FDUMP))
    zt = OE_NEAR_DZ_MM
    cxm, cym = 0.5 * (bx0 + bx1), 0.5 * (by0 + by1)
    CSX.AddDump("E_top", dump_type=10, file_type=1, frequency=fd).AddBox(
        [bx0 - 0.5 * air, by0 - 0.5 * air, zt], [bx1 + 0.5 * air, by1 + 0.5 * air, zt])
    CSX.AddDump("J_cu", dump_type=13, file_type=1, frequency=fd).AddBox([bx0, by0, 0], [bx1, by1, 0])
    CSX.AddDump("E_xz", dump_type=10, file_type=1, frequency=fd).AddBox(
        [bx0 - air, cym, -h - air], [bx1 + air, cym, air_up])
    CSX.AddDump("E_yz", dump_type=10, file_type=1, frequency=fd).AddBox(
        [cxm, by0 - air, -h - air], [cxm, by1 + air, air_up])
    if run:
        os.makedirs(sim_path, exist_ok=True)
        FDTD.Run(sim_path, cleanup=True, verbose=1)
    f = np.linspace(f_lo, f_hi, 401)
    pres = {}
    for name, p in ports:
        p.CalcPort(sim_path, f)
        if np.all(np.abs(p.uf_ref / p.uf_inc) > 0.999):
            log(f"WARNING: port {name} shows |S11| = 1 at every frequency (open/short): the port is not "
                f"connected to the copper in the mesh. Results from this run are NOT valid.")
        pres[name] = {"f": f, "s11": p.uf_ref / p.uf_inc,
                      "pacc": np.real(0.5 * p.uf_tot * np.conj(p.if_tot))}
    return {"ports": pres, "nf2ff": nf2ff, "center": [cxm * U, cym * U, 0.0],
            "board": (bx0, by0, bx1, by1), "copper": allr, "h": h}


def read_fd(path, freq):
    """Reads an openEMS frequency-domain HDF5 dump: returns x, y, z [mm] and field (3, Nx, Ny, Nz)."""
    import h5py
    with h5py.File(path, "r") as hf:
        x, y, z = (np.asarray(hf["Mesh/" + a][()], float) for a in "xyz")
        if max(np.abs(x).max(), np.abs(y).max(), np.abs(z).max()) < 1.0:    # stored in meters
            x, y, z = x * 1e3, y * 1e3, z * 1e3
        g = hf["FieldData/FD"]
        ks = sorted(g.keys(), key=lambda k: int(k[1:]))
        fr = np.array([g[k].attrs["frequency"] for k in ks]).ravel()
        k = ks[int(np.argmin(np.abs(fr - freq)))]
        v = g[k][()]
        if v.dtype.names:
            v = v[v.dtype.names[0]] + 1j * v[v.dtype.names[1]]
        if v.ndim == 4 and v.shape[0] != 3 and v.shape[-1] == 3:
            v = np.moveaxis(v, -1, 0)
        return x, y, z, v, float(fr[int(np.argmin(np.abs(fr - freq)))])


# ============================ CARGA ========================================
base = os.path.join(QUCS_DIR, PROJECT)
out_dir = os.path.join(QUCS_DIR, OUT_DIR)
os.makedirs(out_dir, exist_ok=True)
props, rects, lumped, sub, sweep, rad = parse_pcb(base + ".pcb")
dat = {}
if os.path.exists(base + ".dat"):
    try:
        dat = parse_dat(base + ".dat")
    except Exception as ex:
        print("ERROR reading the .dat (continuing without it):", ex)
else:
    print("WARNING: no .dat (simulate in uSimmics first). Only the layout is drawn.")
ports = [l for l in lumped if l["port"] > 0]
shorts = [l for l in lumped if l["port"] == 0]
print(f"Layout: {len(rects)} rectangles | ports: {[p['name'] for p in ports]} | "
      f"lumped: {[s['name'] for s in shorts]}")
if sub:
    print(f"Substrate: er={sub['er']}, h={sub['h']*1e3:.3f} mm, t={sub['t']*1e6:.0f} um")
print("Variables in the .dat:", {k: (v["deps"], v["dims"], v["v"].shape) for k, v in dat.items()})

xs = np.array([[r[1], r[3]] for r in rects]).ravel()
ys = np.array([[r[2], r[4]] for r in rects]).ravel()
ant_size = max(np.ptp(xs), np.ptp(ys))
cx, cy = (xs.min() + xs.max()) / 2, (ys.min() + ys.max()) / 2

fS, s11 = get_sii(dat, 1)
if F_PATTERN_HZ:
    F0 = F_PATTERN_HZ
elif rad and sweep and sweep[0] <= rad["f"] <= sweep[1]:
    F0 = rad["f"]
elif s11 is not None:
    F0 = fS[int(np.argmin(np.abs(s11)))]
else:
    F0 = rad["f"] if rad else 2.45e9
LAM = C0 / F0
K = 2 * np.pi / LAM
print(f"Pattern frequency: {F0/1e9:.4f} GHz (lambda0 = {LAM*1e3:.2f} mm, antenna = {ant_size*1e3:.1f} mm "
      f"= {ant_size/LAM:.3f} lambda0)")

# ============================ openEMS (if applicable) ======================
ffk_dat = [k for k, v in dat.items() if len(v["deps"]) == 2 and all(ff_axes(v["deps"]))]
OE = None
if PATTERN_SOURCE == "openems" or (PATTERN_SOURCE == "auto" and OE_ENABLE and not ffk_dat):
    if not (sub and sweep):
        print("openEMS: Substrate or Frequency missing in the .pcb")
    else:
        import hashlib
        key = hashlib.md5(open(base + ".pcb", "rb").read() + repr((
            OE_GROUND, OE_EXTRA_TOP, OE_BOARD_MARGIN_MM, OE_RES_MM, OE_END_DB, OE_MAX_TS, OE_AIR_LAM,
            OE_NEAR_DZ_MM, OE_N_FDUMP, PATCH_REGEX, OE_AIR_UP_LAM, "backend-v3")).encode()).hexdigest()[:10]
        sim_path = os.path.join(out_dir, "openems_" + key)
        done = os.path.join(sim_path, "done.txt")
        need_run = OE_FORCE_RERUN or not os.path.exists(done)
        try:
            print("openEMS:", "simulating (may take several minutes)..." if need_run else "using cache " + sim_path)
            OE = run_openems(rects, lumped, sub, sweep[0], sweep[1], sim_path, need_run)
            if need_run:
                open(done, "w").write("ok")
        except ImportError as ex:
            print("openEMS not available in this Python:", ex)
        except Exception as ex:
            print("openEMS failed:", repr(ex))
if OE is not None and OE["ports"]:
    pO = next(iter(OE["ports"].values()))
    iO = int(np.argmin(np.abs(pO["s11"])))
    print(f"openEMS: S11 min = {20*np.log10(abs(pO['s11'][iO])):.2f} dB @ {pO['f'][iO]/1e9:.4f} GHz")
    if not F_PATTERN_HZ and s11 is None and not (rad and sweep and sweep[0] <= rad["f"] <= sweep[1]):
        F0 = pO["f"][iO]; LAM = C0 / F0; K = 2 * np.pi / LAM
        print(f"Pattern frequency (openEMS S11 minimum): {F0/1e9:.4f} GHz")

# ============================ 01 LAYOUT ====================================
fig, ax = plt.subplots(figsize=(8, 8))
for nm, x0, y0, x1, y1 in rects:
    ax.add_patch(Rectangle((x0 * 1e3, y0 * 1e3), (x1 - x0) * 1e3, (y1 - y0) * 1e3, fc="#c8961e", ec="none"))
for p_ in ports:
    ax.plot(p_["x"] * 1e3, p_["y"] * 1e3, "o", color="#e53935", ms=9, mec="w")
    ax.annotate(f"{p_['name']} ({p_['z']:.0f} Ω)", (p_["x"] * 1e3, p_["y"] * 1e3), xytext=(8, 8),
                textcoords="offset points", color="#ff8a80", weight="bold")
for s_ in shorts:
    ax.plot(s_["x"] * 1e3, s_["y"] * 1e3, "s", color="#90caf9", ms=8, mec="k")
    ax.annotate(f"{s_['name']} ({s_['z']:.0f} Ω)", (s_["x"] * 1e3, s_["y"] * 1e3), xytext=(8, -14),
                textcoords="offset points", color="#bbdefb")
ax.set_aspect("equal"); ax.autoscale(); ax.margins(0.08)
ax.set_facecolor("#1f5e2e"); ax.grid(True, color="w", lw=0.2)
ax.set_xlabel("x [mm]"); ax.set_ylabel("y [mm]")
t = f"{PROJECT}   {np.ptp(xs)*1e3:.2f} x {np.ptp(ys)*1e3:.2f} mm"
if sub:
    t += f"\nεr={sub['er']}, h={sub['h']*1e3:.3f} mm, t={sub['t']*1e6:.0f} µm"
ax.set_title(t, fontsize=10)
fig.savefig(os.path.join(out_dir, "01_layout.png"), dpi=150, bbox_inches="tight")

# ============================ 02/03 S11 ====================================
curves = []
if s11 is not None:
    curves.append(("uSimmics", fS, s11, "#1565c0"))
if OE is not None and OE["ports"]:
    curves.append(("openEMS", pO["f"], pO["s11"], "#e65100"))
for lbl, ff_, ss_, _c in curves:
    db = 20 * np.log10(np.abs(ss_) + 1e-15); i = int(np.argmin(db))
    zin = 50 * (1 + ss_[i]) / (1 - ss_[i])
    txt = f"[{lbl}] S11 min = {db[i]:.2f} dB @ {ff_[i]/1e9:.4f} GHz (Zin = {zin.real:.1f}{zin.imag:+.1f}j Ω)"
    below = db <= -10
    if below[i]:
        l = r = i
        while l > 0 and below[l - 1]: l -= 1
        while r < len(db) - 1 and below[r + 1]: r += 1
        txt += f"  BW(-10 dB) = {(ff_[r]-ff_[l])/1e6:.1f} MHz"
    print(txt)
    for bn, (b0, b1) in BANDS_HZ.items():
        mk = (ff_ >= b0) & (ff_ <= b1)
        if mk.any():
            print(f"   worst S11 in {bn}: {db[mk].max():.2f} dB")
for lbl, ff_, ss_, _c in curves:            # numeric export to compare with uSimmics
    zz_ = 50 * (1 + ss_) / (1 - ss_)
    np.savetxt(os.path.join(out_dir, f"02_s11_{lbl}.csv"),
               np.c_[ff_ / 1e9, ss_.real, ss_.imag, 20 * np.log10(np.abs(ss_) + 1e-15), zz_.real, zz_.imag],
               delimiter=";", header="f_GHz;re_S11;im_S11;S11_dB;re_Zin;im_Zin", comments="", fmt="%.6g")
if curves:
    fig, ax = plt.subplots(figsize=(8, 4.5))
    fmin = min(c[1].min() for c in curves); fmax = max(c[1].max() for c in curves)
    for bn, (b0, b1) in BANDS_HZ.items():
        if b1 > fmin and b0 < fmax:
            ax.axvspan(b0 / 1e9, b1 / 1e9, color="#4caf50", alpha=0.2, label=bn)
    for lbl, ff_, ss_, c_ in curves:
        db = 20 * np.log10(np.abs(ss_) + 1e-15); i = int(np.argmin(db))
        ax.plot(ff_ / 1e9, db, lw=2, color=c_, label=f"|S11| {lbl}")
        ax.plot(ff_[i] / 1e9, db[i], "o", color=c_)
        ax.annotate(f"{db[i]:.1f} dB @ {ff_[i]/1e9:.4f} GHz", (ff_[i] / 1e9, db[i]), xytext=(8, -4),
                    textcoords="offset points", color=c_, fontsize=8)
    ax.axhline(-10, color="gray", ls="--", lw=1)
    ax.set_xlabel("Frequency [GHz]"); ax.set_ylabel("|S11| [dB]"); ax.grid(True, alpha=0.4)
    ax.legend(loc="lower left", fontsize=8); ax.set_title(PROJECT, fontsize=10)
    fig.savefig(os.path.join(out_dir, "02_s11.png"), dpi=150, bbox_inches="tight")

    fig, ax = plt.subplots(figsize=(6.5, 6.5))
    tt = np.linspace(0, 2 * np.pi, 400)
    ax.plot(np.cos(tt), np.sin(tt), "k", lw=1)
    for r_ in (0.2, 0.5, 1, 2, 5):
        ax.plot(r_ / (1 + r_) + np.cos(tt) / (1 + r_), np.sin(tt) / (1 + r_), color="gray", lw=0.5)
    for x_ in (0.2, 0.5, 1, 2, 5):
        for sg in (1, -1):
            zz = np.logspace(-3, 3, 600) + 1j * sg * x_
            g = (zz - 1) / (zz + 1); ax.plot(g.real, g.imag, color="gray", lw=0.5)
    ax.plot([-1, 1], [0, 0], color="gray", lw=0.5)
    for lbl, ff_, ss_, c_ in curves:
        ax.plot(ss_.real, ss_.imag, lw=2, color=c_, label=lbl)
        for bn, (b0, b1) in BANDS_HZ.items():
            mk = (ff_ >= b0) & (ff_ <= b1)
            if mk.any():
                ax.plot(ss_.real[mk], ss_.imag[mk], lw=4, color="#4caf50")
    ax.set_aspect("equal"); ax.set_xlim(-1.05, 1.05); ax.set_ylim(-1.05, 1.05); ax.axis("off")
    ax.legend(loc="upper left", fontsize=8); ax.set_title("S11 - Smith chart (Z0 = 50 Ω); green = ISM band")
    fig.savefig(os.path.join(out_dir, "03_smith.png"), dpi=150, bbox_inches="tight")

# ============================ PATTERN SOURCE ===============================
pat = None              # function pat(theta, phi) -> normalized power
full_sphere = False
near_model = None       # function E(PX, PY, ZZ) for patch mode
src = "none"
cut_phi = np.deg2rad(CUT_PHI_DEG) if CUT_PHI_DEG is not None else 0.0

# --- (a) far-field from the .dat
for k in ("directivity_dBi", "directivity", "efficiency", "gain_dBi"):
    if k in dat and dat[k]["v"].size == 1:
        v = float(np.real(dat[k]["v"][0]))
        print(f"uSimmics: {k} = {v:.4g}" + (f"  ({100*v:.1f} %)" if k == "efficiency" and v <= 1 else ""))
if PATTERN_SOURCE in ("auto", "dat"):
    ffk = list(ffk_dat)
    if ffk:
        try:
            e0 = dat[ffk[0]]
            tn, pn = ff_axes(e0["deps"])
            tv = np.real(dat[tn]["v"]).astype(float); pv = np.real(dat[pn]["v"]).astype(float)
            if max(np.abs(tv).max(), np.abs(pv).max()) > 2 * np.pi + 1e-6:
                tv, pv = np.deg2rad(tv), np.deg2rad(pv)
            if re.search("elev", tn, re.I) and not (tv.min() >= -1e-6 and tv.max() > np.pi / 2 + 1e-3):
                print(f"Far-field: '{tn}' = elevation above the xy plane "
                      f"({np.rad2deg(tv.min()):.0f}..{np.rad2deg(tv.max()):.0f} deg) -> theta = 90 deg - elev")
                tv = np.pi / 2 - tv
            print(f"Far-field: theta {np.rad2deg(tv.min()):.0f}..{np.rad2deg(tv.max()):.0f} deg ({tv.size}), "
                  f"phi {np.rad2deg(pv.min()):.0f}..{np.rad2deg(pv.max()):.0f} deg ({pv.size})")

            def grid(k):
                A = dat[k]["v"].reshape(dat[e0["deps"][0]]["v"].size, dat[e0["deps"][1]]["v"].size)
                return A if e0["deps"].index(tn) == 0 else A.T          # -> [theta, phi]

            cplx = [k for k in ffk if np.iscomplexobj(dat[k]["v"])]
            if cplx:            # field components: power = sum of |E|^2
                use = cplx
                Pf = sum(np.abs(grid(k)) ** 2 for k in cplx)
            else:               # real magnitudes only: linear or dB
                use = ffk[:1]
                g = grid(use[0]).real
                isdb = "db" in use[0].lower() or g.min() < 0
                Pf = 10 ** (g / 10) if isdb else np.abs(g)
                print(f"Far-field: using '{use[0]}' as power ({'dB' if isdb else 'linear'})")
            ffk = use
            ot, op = np.argsort(tv), np.argsort(pv)
            tv, pv, Pf = tv[ot], pv[op], Pf[ot][:, op]
            Pf = Pf / Pf.max()
            if pv.max() - pv.min() < 2 * np.pi - 1e-6:          # close the phi periodicity
                pv = np.r_[pv, pv[0] + 2 * np.pi]; Pf = np.c_[Pf, Pf[:, :1]]
            p0 = pv.min()
            pat = lambda th, ph: interp2(tv, pv, Pf, th, np.mod(ph - p0, 2 * np.pi) + p0)
            full_sphere = tv.max() > np.pi / 2 + 1e-3
            src = "uSimmics far-field (" + ", ".join(ffk) + ")"
        except Exception as ex:
            print("Could not use the far-field from the .dat:", repr(ex))
    elif PATTERN_SOURCE == "dat":
        print("The .dat has no far-field. Enable 'Radiation' in uSimmics and re-simulate.")

# --- (a2) far-field from openEMS (NF2FF)
if pat is None and OE is not None and PATTERN_SOURCE in ("auto", "openems"):
    try:
        thd, phd = np.arange(0, 181, 2.0), np.arange(0, 361, 4.0)
        ffr = OE["nf2ff"].CalcNF2FF(sim_path, F0, thd, phd, center=OE["center"], read_cached=False)
        Pf_oe = np.abs(ffr.E_norm[0]) ** 2
        Pf_oe = Pf_oe / Pf_oe.max()
        tv_oe, pv_oe = np.deg2rad(thd), np.deg2rad(phd)
        pat = lambda th, ph: interp2(tv_oe, pv_oe, Pf_oe, th, np.mod(ph, 2 * np.pi))
        full_sphere = True
        pacc0 = np.interp(F0, pO["f"], pO["pacc"])
        eff = ffr.Prad[0] / pacc0 if pacc0 > 0 else np.nan
        src = "openEMS FDTD (NF2FF)"
        s0 = 20 * np.log10(abs(np.interp(F0, pO["f"], np.abs(pO["s11"]))))
        print(f"openEMS NF2FF @ {F0/1e9:.4f} GHz: Dmax = {10*np.log10(ffr.Dmax[0]):.2f} dBi, "
              f"radiation efficiency ~ {100*eff:.1f} % (mismatch loss excluded)")
        if s0 > -3:
            print(f"WARNING: S11 = {s0:.1f} dB at that frequency: the antenna is mismatched; the efficiency "
                  f"(Prad/Pacc) is numerically unreliable and the pattern is not the operating one.")
    except Exception as ex:
        print("NF2FF failed:", repr(ex))

# --- (b) analytical patch model
if pat is None and PATTERN_SOURCE in ("auto", "patch") and sub:
    groups = {}
    for nm, x0, y0, x1, y1 in rects:
        m = re.match(PATCH_REGEX, nm)
        if m:
            groups.setdefault(m.group(1), []).append((x0, y0, x1, y1))
    if groups:
        pn_ = sorted(groups, key=lambda s: int(re.sub(r"\D", "", s) or 0))
        N = len(pn_)
        bb = np.array([[min(r[0] for r in groups[n]), min(r[1] for r in groups[n]),
                        max(r[2] for r in groups[n]), max(r[3] for r in groups[n])] for n in pn_])
        cen = np.c_[(bb[:, 0] + bb[:, 2]) / 2, (bb[:, 1] + bb[:, 3]) / 2]
        dimx, dimy = (bb[:, 2] - bb[:, 0]).mean(), (bb[:, 3] - bb[:, 1]).mean()
        fx, fy = f_res(dimx, dimy, sub["er"], sub["h"]), f_res(dimy, dimx, sub["er"], sub["h"])
        pol = POL_AXIS if POL_AXIS in ("x", "y") else ("x" if abs(fx - F0) < abs(fy - F0) else "y")
        Lp, Wp = (dimx, dimy) if pol == "x" else (dimy, dimx)
        ROT = 0.0 if pol == "x" else np.pi / 2
        pos = cen - cen.mean(axis=0)
        A = np.ones(N) if ELEM_AMPL is None else np.asarray(ELEM_AMPL, float)
        PH = np.zeros(N) if ELEM_PHASE_DEG is None else np.deg2rad(ELEM_PHASE_DEG)
        print(f"Patch mode: {pn_}, {dimx*1e3:.2f}x{dimy*1e3:.2f} mm; cavity model fx={fx/1e9:.2f} GHz, "
              f"fy={fy/1e9:.2f} GHz -> E-plane along {pol}")
        if N > 1:
            u = cen[-1] - cen[0]; u /= np.linalg.norm(u)
            if CUT_PHI_DEG is None:
                cut_phi = np.arctan2(u[1], u[0]) % np.pi
            print(f"Mean pitch = {np.mean(np.linalg.norm(np.diff(cen, axis=0), axis=1))*1e3:.2f} mm")
        if ELEM_PHASE_DEG is None:
            print("WARNING: element phases assumed equal (feed network not solved).")

        def element(theta, phi):
            ph = phi - ROT
            st, ct = np.sin(theta), np.cos(theta)
            slot = np.sinc(K * Wp / 2 * st * np.sin(ph) / np.pi)
            pair = np.cos(K * Lp / 2 * st * np.cos(ph))
            polf = np.sqrt(np.cos(ph) ** 2 + (ct * np.sin(ph)) ** 2)
            return np.where(theta <= np.pi / 2, np.abs(slot * pair * polf), 0.0)

        def field_ff(theta, phi):
            s = np.zeros(np.broadcast(theta, phi).shape, complex)
            for n in range(N):
                s += A[n] * np.exp(1j * (K * np.sin(theta) * (pos[n, 0] * np.cos(phi) +
                                                              pos[n, 1] * np.sin(phi)) + PH[n]))
            return element(theta, phi) * np.abs(s)

        _th = np.linspace(0, np.pi / 2, 181)
        _P = field_ff(*np.meshgrid(_th, np.linspace(0, 2 * np.pi, 361), indexing="ij")) ** 2
        _pmax = _P.max()
        pat = lambda th, ph: field_ff(th, ph) ** 2 / _pmax
        full_sphere = False
        acx, acy = cen.mean(axis=0)

        def near_model(PX, PY, ZZ):
            E = np.zeros_like(ZZ, complex)
            for n in range(N):
                dx, dy = PX - pos[n, 0], PY - pos[n, 1]
                r = np.sqrt(dx ** 2 + dy ** 2 + ZZ ** 2)
                E += A[n] * np.exp(1j * PH[n]) * element(np.arccos(ZZ / r), np.arctan2(dy, dx)) * \
                    np.exp(-1j * K * r) / r
            return E
        src = "analytical (patch)"

if pat is None:
    print("\nNo radiation pattern. To get it for ANY antenna: in uSimmics, in the\n"
          "layout properties enable 'Radiation' at the frequency of interest and re-simulate.\n"
          "The script detects the far-field in the .dat automatically.")
    print("Files in:", os.path.abspath(out_dir))
    if SHOW:
        plt.show()
    raise SystemExit(0)

print("Pattern source:", src)

# ============================ METRICS AND CUTS =============================
tmax = np.pi if full_sphere else np.pi / 2
th = np.linspace(0, tmax, 361)
ph = np.linspace(0, 2 * np.pi, 721)
TH, PHI = np.meshgrid(th, ph, indexing="ij")
Pg = pat(TH, PHI)
D = 4 * np.pi * Pg.max() / (np.sum(Pg * np.sin(TH)) * (th[1] - th[0]) * (ph[1] - ph[0]))


def cut(phi0):
    tt = np.linspace(-tmax, tmax, 1801)
    e2 = pat(np.abs(tt), np.where(tt >= 0, phi0, phi0 + np.pi))
    return tt, 10 * np.log10(e2 / max(e2.max(), 1e-30) + 1e-12)


cuts = [(0.0, "φ = 0° (xz plane)"), (np.pi / 2, "φ = 90° (yz plane)")]
fig, ax = plt.subplots(subplot_kw={"projection": "polar"}, figsize=(6.5, 6.5))
bws = []
for phi0, lab in cuts:
    tt, g = cut(phi0)
    bw = hpbw(tt[:-1], g[:-1], circular=full_sphere) if full_sphere else hpbw(tt, g)
    bws.append(bw)
    txt = "omni" if bw > 300 else f"HPBW {bw:.1f}°"
    print(f"{lab}: {txt}")
    ax.plot(tt, np.clip(g, -DB_RANGE, 0), lw=2, label=f"{lab} - {txt}")
ax.set_theta_zero_location("N"); ax.set_theta_direction(-1)
if not full_sphere:
    ax.set_thetamin(-90); ax.set_thetamax(90)
ax.set_rlim(-DB_RANGE, 0)
ax.set_title(f"@ {F0/1e9:.4f} GHz   D ≈ {10*np.log10(D):.2f} dBi\n[{src}]", fontsize=9)
ax.legend(loc="lower center", bbox_to_anchor=(0.5, -0.18), fontsize=8)
fig.savefig(os.path.join(out_dir, "04_cuts_polar.png"), dpi=150, bbox_inches="tight")
print(f"Directivity ≈ {10*np.log10(D):.2f} dBi [{src}]")

# ============================ 3D ===========================================
th3 = np.linspace(0, tmax, 91 if not full_sphere else 121)
ph3 = np.linspace(0, 2 * np.pi, 181)
T3, P3 = np.meshgrid(th3, ph3, indexing="ij")
g3 = 10 * np.log10(pat(T3, P3) + 1e-12)
r3 = np.clip((g3 + DB_RANGE) / DB_RANGE, 0, 1)
SCALE = 1.4 * ant_size * 1e3
X3 = cx * 1e3 + SCALE * r3 * np.sin(T3) * np.cos(P3)
Y3 = cy * 1e3 + SCALE * r3 * np.sin(T3) * np.sin(P3)
Z3 = SCALE * r3 * np.cos(T3)
norm = colors.Normalize(-DB_RANGE, 0)
fc = plt.get_cmap(CMAP)(norm(np.clip(g3, -DB_RANGE, 0)))


def draw_3d(fig):
    ax = fig.add_subplot(projection="3d")
    mg = 0.15 * ant_size * 1e3
    b = (xs.min() * 1e3 - mg, xs.max() * 1e3 + mg, ys.min() * 1e3 - mg, ys.max() * 1e3 + mg)
    ax.add_collection3d(Poly3DCollection([[(b[0], b[2], 0), (b[1], b[2], 0), (b[1], b[3], 0), (b[0], b[3], 0)]],
                                         facecolor="#1f5e2e", alpha=0.9))
    cu = [[(r[1] * 1e3, r[2] * 1e3, 0.02), (r[3] * 1e3, r[2] * 1e3, 0.02),
           (r[3] * 1e3, r[4] * 1e3, 0.02), (r[1] * 1e3, r[4] * 1e3, 0.02)] for r in rects]
    ax.add_collection3d(Poly3DCollection(cu, facecolor="#d4a017", edgecolor="none"))
    ax.plot_surface(X3, Y3, Z3, facecolors=fc, rstride=1, cstride=1, linewidth=0,
                    antialiased=False, shade=False, alpha=0.8)
    fig.colorbar(cm.ScalarMappable(norm=norm, cmap=CMAP), ax=ax, shrink=0.6, label="Relative gain [dB]")
    L3 = SCALE
    ax.set_xlim(cx * 1e3 - L3, cx * 1e3 + L3); ax.set_ylim(cy * 1e3 - L3, cy * 1e3 + L3)
    ax.set_zlim(-L3 if full_sphere else 0, L3 if full_sphere else 2 * L3)
    ax.set_box_aspect((1, 1, 1 if full_sphere else 1))
    ax.set_xlabel("x [mm]"); ax.set_ylabel("y [mm]"); ax.set_zlabel("z [mm]")
    ax.set_title(f"3D pattern @ {F0/1e9:.4f} GHz [{src}]", fontsize=9)
    ax.view_init(elev=22, azim=-60)
    return ax


fig3 = plt.figure(figsize=(8, 7.5))
ax3 = draw_3d(fig3)          # static window: rotate with the mouse
fig3.savefig(os.path.join(out_dir, "05_pattern_3d.png"), dpi=150, bbox_inches="tight")

# ============================ FIELD ON A VERTICAL PLANE ===================
if CUT_PHI_DEG is None and near_model is None:      # plane of the narrowest beam (the most illustrative)
    cut_phi = cuts[int(np.argmin(bws))][0]
    print(f"Field-map plane: φ = {np.rad2deg(cut_phi):.0f}° (narrowest beam)")
u_c = np.array([np.cos(cut_phi), np.sin(cut_phi)])
half, res = FIELD_HALF_LAM * LAM, FIELD_RES_LAM * LAM
s = np.arange(-half, half + res, res)
z = np.arange(-half if full_sphere else res / 2, half + res, res)
SS, ZZ = np.meshgrid(s, z)
R = np.sqrt(SS ** 2 + ZZ ** 2)
if near_model is not None:
    E = near_model(SS * u_c[0], SS * u_c[1], ZZ) * R
    mask = np.zeros_like(R, bool)
    fld_lbl = "analytical model (superposition of sources)"
else:
    thp = np.arccos(np.clip(ZZ / np.maximum(R, 1e-12), -1, 1))
    php = np.where(SS >= 0, cut_phi, cut_phi + np.pi)
    E = np.sqrt(pat(thp, php)) * np.exp(-1j * K * R)          # |E|*r from the far-field
    r_nf = max(LAM / (2 * np.pi), ant_size)
    mask = R < r_nf
    fld_lbl = f"synthesized from the far-field (r < {r_nf*1e3:.0f} mm not modeled)"
ref = np.abs(E[(ZZ > LAM / 2) & ~mask]).max() if np.any((ZZ > LAM / 2) & ~mask) else np.abs(E).max()
Edb = 20 * np.log10(np.abs(E) / ref + 1e-12)
Edb[mask] = np.nan
ext = [s[0] * 1e3, s[-1] * 1e3, z[0] * 1e3, z[-1] * 1e3]
cmap_nan = plt.get_cmap(CMAP).copy(); cmap_nan.set_bad("#555555")
plane_txt = f"vertical plane at φ = {np.rad2deg(cut_phi):.0f}°"


def decorate(ax):
    ax.axhline(0, color="#1f5e2e", lw=3)
    if mask.any():
        ax.add_patch(Circle((0, 0), max(LAM / (2 * np.pi), ant_size) * 1e3, fill=False, ec="w", ls="--", lw=1))
    ax.set_xlabel(f"Horizontal distance in the {plane_txt} [mm]"); ax.set_ylabel("z [mm]")


fig2, ax2 = plt.subplots(figsize=(8, 6.5))
im = ax2.imshow(np.clip(Edb, -DB_RANGE, 0), origin="lower", extent=ext, cmap=cmap_nan,
                vmin=-DB_RANGE, vmax=0, aspect="equal")
decorate(ax2)
ax2.set_title(f"|E|·r - {fld_lbl}", fontsize=9)
fig2.colorbar(im, ax=ax2, label="dB rel.")
fig2.savefig(os.path.join(out_dir, "06_field_map_db.png"), dpi=150, bbox_inches="tight")

En = np.clip(np.abs(E) / ref, 0, 1) * np.exp(1j * np.angle(E))
w_amp = np.sqrt(np.abs(En))


def wave_frame(i):
    f = np.real(En * np.exp(1j * 2 * np.pi * i / N_FRAMES)) / (np.abs(En) + 1e-9) * w_amp
    f[mask] = np.nan
    return f


def draw_wave(fig):
    ax = fig.add_subplot()
    cmw = plt.get_cmap("RdBu_r").copy(); cmw.set_bad("#555555")
    im_ = ax.imshow(wave_frame(0), origin="lower", extent=ext, cmap=cmw, vmin=-1, vmax=1, aspect="equal")
    decorate(ax)
    ax.set_title(f"Re{{E·e^(jωt)}} @ {F0/1e9:.4f} GHz - {fld_lbl}", fontsize=9)
    return im_


if MAKE_GIFS:
    def save_gif(fig, func, name):
        path = os.path.join(out_dir, name)
        try:
            animation.FuncAnimation(fig, func, frames=N_FRAMES, blit=False).save(
                path, writer=animation.PillowWriter(fps=12))
            print("GIF written:", os.path.abspath(path))
        except Exception as ex:
            print(f"ERROR generating {name}: {ex!r}  (Pillow installed? pip install pillow)")

    fw = Figure(figsize=(7, 6)); FigureCanvasAgg(fw)
    imw = draw_wave(fw)
    save_gif(fw, lambda i: (imw.set_array(wave_frame(i)),), "wave.gif")
    if MAKE_ROT_GIF:
        fr = Figure(figsize=(8, 7.5)); FigureCanvasAgg(fr)
        axr = draw_3d(fr)
        save_gif(fr, lambda i: (axr.view_init(elev=22, azim=-60 + 360 * i / N_FRAMES),), "pattern_3d_rot.gif")


# ============================ REAL NEAR FIELD (openEMS) ===================
if OE is not None:
    import traceback
    bx0, by0, bx1, by1 = OE["board"]

    def outline(ax):
        for r in OE["copper"]:
            ax.add_patch(Rectangle((r[0], r[1]), r[2] - r[0], r[3] - r[1], fill=False, ec="w", lw=0.6))
        ax.add_patch(Rectangle((bx0, by0), bx1 - bx0, by1 - by0, fill=False, ec="#9e9e9e", lw=0.8, ls="--"))
        ax.set_aspect("equal"); ax.set_xlabel("x [mm]"); ax.set_ylabel("y [mm]")

    def heat(ax, x, y, V, title):
        im_ = ax.pcolormesh(x, y, np.clip(20 * np.log10(V / np.nanmax(V) + 1e-12), -DB_RANGE, 0).T,
                            cmap=CMAP, vmin=-DB_RANGE, vmax=0, shading="nearest")
        ax.set_title(title, fontsize=9)
        return im_

    def uniform(x, y, F, n=400):
        """Resamples a complex field (3, Nx, Ny) onto a uniform grid for imshow/animation."""
        if np.ptp(x) >= np.ptp(y):
            xu = np.linspace(x[0], x[-1], n); yu = np.linspace(y[0], y[-1], max(2, int(n * np.ptp(y) / np.ptp(x))))
        else:
            yu = np.linspace(y[0], y[-1], n); xu = np.linspace(x[0], x[-1], max(2, int(n * np.ptp(x) / np.ptp(y))))
        XU, YU = np.meshgrid(xu, yu, indexing="ij")
        G = np.array([interp2(x, y, F[c].real, XU, YU) + 1j * interp2(x, y, F[c].imag, XU, YU) for c in range(3)])
        return xu, yu, G

    def animate(xu, yu, G, title, name, deco):
        ref = np.sqrt(np.sum(np.abs(G) ** 2, axis=0)).max()

        def frame(i):
            inst = np.sqrt(np.sum(np.real(G * np.exp(1j * 2 * np.pi * i / N_FRAMES)) ** 2, axis=0))
            return np.clip(20 * np.log10(inst / ref + 1e-12), -DB_RANGE, 0).T

        fn = Figure(figsize=(8, 7)); FigureCanvasAgg(fn)
        axn = fn.add_subplot()
        imn = axn.imshow(frame(0), origin="lower", extent=[xu[0], xu[-1], yu[0], yu[-1]],
                         cmap=CMAP, vmin=-DB_RANGE, vmax=0)
        deco(axn); fn.colorbar(imn, ax=axn, label="dB rel.")
        axn.set_title(title, fontsize=9)
        save_gif(fn, lambda i: (imn.set_array(frame(i)),), name)

    def step(name, fn_):
        try:
            fn_()
        except Exception:
            print(f"ERROR in {name}:\n" + traceback.format_exc())

    def fig07():
        x, y, z, Ev, fE = read_fd(os.path.join(sim_path, "E_top.h5"), F0)
        Emag = np.sqrt(np.sum(np.abs(Ev[..., 0]) ** 2, axis=0))
        fig, ax = plt.subplots(figsize=(8, 7))
        im = heat(ax, x, y, Emag, f"|E| a {OE_NEAR_DZ_MM} mm above the copper @ {fE/1e9:.4f} GHz (openEMS)")
        outline(ax); fig.colorbar(im, ax=ax, label="dB rel.")
        fig.savefig(os.path.join(out_dir, "07_nearfield_E_top.png"), dpi=150, bbox_inches="tight")
        print("Figure: 07_nearfield_E_top.png")
        if MAKE_GIFS:
            xu, yu, G = uniform(x, y, Ev[..., 0])
            animate(xu, yu, G, f"|E(t)| over the board @ {fE/1e9:.4f} GHz (openEMS)", "nearfield.gif", outline)

    def fig08():
        x, y, z, Jv, fJ = read_fd(os.path.join(sim_path, "J_cu.h5"), F0)
        Jmag = np.sqrt(np.sum(np.abs(Jv[..., 0]) ** 2, axis=0))
        fig, ax = plt.subplots(figsize=(8, 7))
        im = heat(ax, x, y, Jmag, f"|J| on the copper @ {fJ/1e9:.4f} GHz (openEMS)")
        outline(ax); fig.colorbar(im, ax=ax, label="dB rel.")
        fig.savefig(os.path.join(out_dir, "08_surface_current_J.png"), dpi=150, bbox_inches="tight")
        print("Figure: 08_surface_current_J.png")

    def fig09():
        vert = "E_yz" if abs(np.sin(cut_phi)) > 0.7 else "E_xz"
        xv, yv, zv, Vv, fV = read_fd(os.path.join(sim_path, vert + ".h5"), F0)
        hor = yv if vert == "E_yz" else xv
        Vp = Vv[:, 0, :, :] if vert == "E_yz" else Vv[:, :, 0, :]
        k_ = 3                                     # crop: cells next to the MUR boundary show artifacts
        hor, zv, Vp = hor[k_:-k_], zv[k_:-k_], Vp[:, k_:-k_, k_:-k_]
        c0 = 0.5 * (by0 + by1) if vert == "E_yz" else 0.5 * (bx0 + bx1)
        hl = ("y" if vert == "E_yz" else "x") + " [mm]"
        lam_mm = C0 / fV * 1e3

        def deco(ax):
            ax.axhline(0, color="k", lw=0.6); ax.axhline(-OE["h"], color="k", lw=0.6)
            ax.set_aspect("equal"); ax.set_xlabel(hl); ax.set_ylabel("z [mm]")

        # |E|*r: compensates the 1/r spreading so the beam is visible, not only the reactive field
        HH, ZZv = np.meshgrid(hor - c0, zv, indexing="ij")
        Rv = np.sqrt(HH ** 2 + ZZv ** 2)
        Vmag = np.sqrt(np.sum(np.abs(Vp) ** 2, axis=0)) * Rv
        far = Rv > lam_mm / 4
        fig, ax = plt.subplots(figsize=(8, 6.5))
        im = ax.pcolormesh(hor, zv, np.clip(20 * np.log10(Vmag / Vmag[far].max() + 1e-12), -DB_RANGE, 0).T,
                           cmap=CMAP, vmin=-DB_RANGE, vmax=0, shading="nearest")
        ax.set_title(f"|E|·r on the {vert[2:]} plane @ {fV/1e9:.4f} GHz (openEMS)", fontsize=9)
        deco(ax); fig.colorbar(im, ax=ax, label="dB rel.")
        fig.savefig(os.path.join(out_dir, "09_nearfield_E_side.png"), dpi=150, bbox_inches="tight")
        print("Figure: 09_nearfield_E_side.png")
        if MAKE_GIFS:
            # real wave: dominant component away from the board, Re{E e^jwt}*r with compression
            comp = int(np.argmax([np.sum(np.abs(Vp[c][far]) ** 2) for c in range(3)]))
            xu, zu, G = uniform(hor, zv, np.array([Vp[comp] * Rv, Vp[comp] * 0, Vp[comp] * 0]))
            Gc = G[0]
            XU, ZU = np.meshgrid(xu - c0, zu, indexing="ij")
            ref = np.abs(Gc[np.sqrt(XU ** 2 + ZU ** 2) > lam_mm / 4]).max()

            def frame(i):
                v = np.real(Gc * np.exp(1j * 2 * np.pi * i / N_FRAMES)) / ref
                return (np.sign(v) * np.sqrt(np.clip(np.abs(v), 0, 1))).T

            fn = Figure(figsize=(8, 7)); FigureCanvasAgg(fn)
            axn = fn.add_subplot()
            imn = axn.imshow(frame(0), origin="lower", extent=[xu[0], xu[-1], zu[0], zu[-1]],
                             cmap="RdBu_r", vmin=-1, vmax=1)
            deco(axn)
            axn.set_title(f"Re{{E{'xyz'[comp]}·e^(jωt)}}·r on the {vert[2:]} plane @ {fV/1e9:.4f} GHz (openEMS)",
                          fontsize=9)
            save_gif(fn, lambda i: (imn.set_array(frame(i)),), "nearfield_side.gif")

    if not MAKE_GIFS:
        def save_gif(*a):
            pass
    step("07 / nearfield.gif", fig07)
    step("08", fig08)
    step("09 / nearfield_side.gif", fig09)

print("Files in:", os.path.abspath(out_dir))
if SHOW:
    if SHOW_WAVE_LIVE:
        figw = plt.figure(figsize=(7, 6))
        imlive = draw_wave(figw)
        wave_anim = animation.FuncAnimation(figw, lambda i: (imlive.set_array(wave_frame(i)),),
                                            frames=N_FRAMES, interval=80, blit=False)
    plt.show()
