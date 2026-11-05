#!/usr/bin/env python3
"""Full-page cover art: Wilson Hall built from quantum-computing schematics.

Blueprint idiom -- deep FNAL-blue plate, white/pale-blue linework, sparse
orange accents. Complexity rises with altitude under CONTINUOUS WEIGHTS:
every row glyph carries a Gaussian preference over floor number, so ground
hardware (junctions, SQUIDs, spins, trapped ions, cold atoms, 1-2 qubit
gates) fades gradually into NISQ objects, then integrated hardware and
register-level primitive-gate circuits with the SQMS SRF cavities, then
abstraction (tensor networks, plaquettes, braids), with the deterministic
set pieces above: the cryo plant, and the QEC story (surface code, 15:1
factory, lattice surgery, logical qubit), and the gauge-theory lattice
glowing at the apex.

BUILDING GEOMETRY, measured from the reference photo of the elevation
(OIP-1016144596.jpg, repo root): base half-width 0.37 H, top 0.25 H, outer
flare concentrated in the bottom third (P = 3.2). The atrium glass is a
third of the base width, tapers out by t = T1, continues as a narrow
window strip to t = T2, and ENDS there -- solid concrete above, with the
parapet nubs on the roof.

The page is full-bleed letter (8.5 x 11); the report title is typeset over
the sky by LaTeX (\\FNALCoverPage in fnal-report.cls), so the art stays
text-free for SVG editing. Same seed, same art.
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle, Ellipse, PathPatch, Rectangle
from matplotlib.path import Path

import cartoon as ct
import cover_glyphs as cg
import fnalfig as ff
from cover_piece_bubble import draw_bubble_chamber
from cover_piece_cosmology import draw_cosmology
from cover_piece_detector import draw_detector_slice
from cover_piece_event import draw_event_display
from cover_piece_latticeqcd import draw_lattice_qcd
from cover_piece_loopcloud import draw_loop_cloud
from cover_piece_nunucleus import draw_nu_nucleus

# ------------------------------------------------------------ palette -----
PLATE = "#0d2c56"
WHITE = "#f2f6fb"
PALE = "#9fc4e8"
ORANGE = "#e8833a"
GLOW = "#f4a259"

# ------------------------------------------------------------ canvas ------
PW, PH = 8.5, 11.0
FRAME1, FRAME2 = 0.10, 0.18
GRID = 0.50

GROUND = 1.95
POND0 = 0.38
H = 5.65
CX = 4.25

TITLE_BAND = 1.75        # sky depth under the top frame reserved for LaTeX

# ---------------------------------------------------- building curves -----
W0, W1, P = 2.51, 1.79, 3.2   # 0.445 H / 0.316 H, measured
S0 = 0.99                # glass half-width at ground: 0.175 H, measured
T1, QTAPER = 0.42, 1.0   # LINEAR taper (straight edges), ends t=0.42
T2, SMIN = 0.85, 0.48    # constant strip 0.085 H to the sill at t=0.85
NUB_WB, NUB_WT, NUB_H = 0.73, 0.57, 0.40  # base/top widths + height, measured
NFLOOR = 16
FLOOR_H = H / NFLOOR
CROSSOVERS = (4,)   # one visible crossover band, t~0.23

RNG_SEED = 20260806


def x_out(t):
    return W1 + (W0 - W1) * (1 - np.asarray(t)) ** P


def s_slot(t):
    """Gap half-width: eased glass taper to SMIN by T1, then the open strip
    continues to the roof. Glass ends at T2; above that the gap is sky."""
    t = np.asarray(t, dtype=float)
    taper = SMIN + (S0 - SMIN) * np.clip(1 - (t / T1) ** QTAPER, 0, None)
    return np.where(t < T1, taper, SMIN)


# ------------------------------------------------------------ layout ------
# Continuous altitude weights: WEIGHTS[name] = (mu, sigma, amplitude) over
# floor number; a row at floor f samples every glyph within 2.5 sigma.
WEIGHTS = {
    "jjunction":   (0.5, 2.0, 0.9),
    "transmon":    (1.0, 2.4, 1.0),
    "squid":       (1.6, 2.2, 0.8),
    "spin_qubit":  (1.3, 2.0, 0.7),
    "blade_trap":  (2.2, 2.2, 0.8),
    "cold_atoms":  (2.8, 2.4, 0.9),
    "gateH":       (1.8, 2.8, 0.9),
    "gateX":       (2.0, 2.8, 0.7),
    "gateT":       (2.4, 2.8, 0.7),
    "gateRz":      (3.0, 2.8, 0.6),
    "cnot":        (2.4, 3.0, 0.9),
    "cz":          (2.2, 2.8, 0.7),
    "swap":        (3.4, 2.8, 0.7),
    "meter":       (2.6, 3.0, 0.7),
    "bloch":       (4.6, 2.2, 0.9),
    "ion_chain":   (4.4, 2.2, 0.9),
    "resonator":   (4.2, 2.4, 0.8),
    "chip5":       (5.2, 2.0, 0.9),
    "toffoli":     (5.0, 2.2, 0.8),
    "circuit3":    (5.6, 2.4, 0.9),
    "laser":       (3.2, 2.0, 0.8),
    "atom":        (3.0, 2.2, 0.9),
    "pulse_seq":   (6.5, 2.5, 1.0),
    "circuit5":    (8.5, 2.5, 1.1),
    "hemt_amp":    (5.5, 2.2, 0.7),
    "wirebond_chip": (6.0, 2.0, 0.8),
    "circulator":  (6.5, 2.0, 0.7),
    "tweezer":     (7.0, 1.9, 0.9),
    "heavy_hex":   (7.2, 1.8, 1.0),
    "qpu_die":     (7.8, 1.9, 1.0),
    "srf_cavity":  (10.0, 2.0, 2.2),
    "group_circuit": (8.4, 1.9, 1.3),
    "tensor_net":  (9.8, 2.4, 0.9),
    "braid":       (10.0, 2.2, 0.8),
    "plaquette":   (10.2, 2.4, 0.9),
    "sho":         (4.0, 2.5, 0.9),
    "cryostat":    (10.6, 2.0, 0.65),
    "qldpc":       (11.6, 2.0, 1.0),
    "code713":     (9.6, 2.4, 0.95),
    "code_gross":  (11.8, 2.4, 1.1),
    "code_iceberg": (8.4, 2.4, 0.85),
    "code513":     (8.9, 2.4, 0.85),
    "code_shor9":  (9.9, 2.4, 0.9),
    "code_rm15":   (10.7, 2.4, 0.9),
    "code_golay":  (11.3, 2.4, 0.9),
    "code_qudit":  (12.0, 2.4, 0.9),
    "code_nkdq":   (12.3, 2.4, 0.95),
    "trotter_block": (10.6, 2.4, 1.0),
    "peps":        (10.2, 2.4, 0.9),
    "mpo":         (10.4, 2.4, 0.9),
    "nu_e":        (13.2, 2.0, 0.6),
    "nu_mu":       (13.8, 2.0, 0.6),
    "nu_tau":      (14.4, 2.0, 0.6),
    "tn_loop":     (10.8, 2.3, 0.95),
    "mera":        (10.9, 2.4, 0.9),
    "braket_bell": (13.4, 2.2, 0.8),
    "braket":      (12.6, 2.2, 0.9),
}

# The mosaic region (floors 0-10) is filled by dart-throwing rather than
# rows: positions are uniform in the tower cross-section, sizes follow the
# altitude with jitter, and a candidate is kept only if it clears every
# accepted glyph. Deterministic via the seed.
SCATTER_FLOORS = (0.0, 15.85)
SCATTER_N = 380
SIZE_LO, SIZE_HI = 0.155, 0.375
SIZE_JITTER = (0.86, 1.16)
GAP = 0.035               # clear space demanded between glyph bounding boxes
SAME_MIN = 0.85           # minimum distance between two copies of one glyph
# Second fill pass: small circuits and notation packed into the gaps the
# big-glyph scatter leaves on the upper floors.
FILLER = dict(floors=(10.3, 15.55), n=70, size=(0.19, 0.26),
              menu=["circuit3", "toffoli", "swap", "cnot", "gateT",
                    "braket", "nu_e", "nu_mu", "nu_tau", "code713",
                    "code513", "plaquette", "braid", "tn_loop",
                    "resonator", "cz"])
# Text glyphs are far more recognisable than pictorial ones, so their
# repeats must sit much further apart.
SAME_MIN_MAP = {n: 1.95 for n in
                ("code713", "code_gross", "code_iceberg", "code513",
                 "code_shor9", "code_rm15", "code_golay", "code_qudit",
                 "code_nkdq", "nu_osc", "nu_e", "nu_mu", "nu_tau", "u_nu",
                 "pmns_matrix",
                 "braket",
                 "braket_bell", "ghz", "shor_eq", "kubo", "qcd_ham")}
SIZE_MULT = {"qcd_ham": 1.55, "shor_eq": 1.30, "kubo": 1.5,
             "braket_bell": 1.15, "ghz": 1.2, "code713": 1.25,
             "code_gross": 1.3, "code_iceberg": 1.2, "code513": 1.25,
             "code_shor9": 1.25, "code_rm15": 1.25, "code_golay": 1.25,
             "code_qudit": 1.3, "code_nkdq": 1.25, "nu_osc": 1.35, "u_nu": 1.25, "pmns_matrix": 1.35,
             "nu_e": 1.1, "nu_mu": 1.1, "nu_tau": 1.1,
             "circuit5": 1.30, "pulse_seq": 1.22, "molecule": 1.10,
             "braket": 1.15, "atom": 1.05}

# Half-extents (width, height) of each glyph's bounding box in units of its
# size parameter. Collision is tested on these rectangles, not on circles,
# so wide glyphs (cavities, circuits, equations) cannot kiss their
# neighbors diagonally.
EXTENT = {
    "transmon": (0.45, 0.45), "jjunction": (0.50, 0.30),
    "gateH": (0.52, 0.32), "gateX": (0.52, 0.32), "gateT": (0.52, 0.32),
    "gateRz": (0.52, 0.32), "cnot": (0.50, 0.48), "cz": (0.44, 0.36),
    "swap": (0.44, 0.37), "toffoli": (0.46, 0.41), "meter": (0.36, 0.28),
    "bloch": (0.42, 0.42), "ion_chain": (0.50, 0.20),
    "chip5": (0.46, 0.46), "resonator": (0.47, 0.18),
    "circuit3": (0.54, 0.33), "heavy_hex": (0.38, 0.38),
    "srf_cavity": (0.64, 0.29), "tweezer": (0.42, 0.44),
    "qpu_die": (0.48, 0.48), "group_circuit": (0.53, 0.35),
    "squid": (0.40, 0.32), "spin_qubit": (0.32, 0.30),
    "blade_trap": (0.36, 0.40), "cold_atoms": (0.35, 0.35),
    "tensor_net": (0.40, 0.24), "plaquette": (0.37, 0.37),
    "braid": (0.36, 0.40), "wirebond_chip": (0.46, 0.46),
    "hemt_amp": (0.50, 0.28), "circulator": (0.46, 0.30),
    "laser": (0.48, 0.17), "atom": (0.48, 0.30), "molecule": (0.45, 0.35),
    "circuit5": (0.64, 0.40), "pulse_seq": (0.57, 0.42),
    "braket": (0.56, 0.25), "shor_eq": (0.50, 0.42),
    "qcd_ham": (0.72, 0.38), "kubo": (0.68, 0.44),
    "peps": (0.44, 0.34), "mera": (0.44, 0.42), "braket_bell": (0.42, 0.34),
    "mpo": (0.50, 0.24), "tn_loop": (0.46, 0.30),
    "nu_e": (0.30, 0.22), "nu_mu": (0.30, 0.22), "nu_tau": (0.30, 0.22),
    "nu_osc": (0.66, 0.24), "u_nu": (0.52, 0.22),
    "pmns_matrix": (0.64, 0.36),
    "code713": (0.52, 0.20), "code_gross": (0.72, 0.20),
    "code_iceberg": (0.46, 0.20), "code513": (0.52, 0.20),
    "code_shor9": (0.52, 0.20), "code_rm15": (0.56, 0.20),
    "code_golay": (0.56, 0.20), "code_qudit": (0.68, 0.24),
    "code_nkdq": (0.58, 0.24),
    "wilson_loop": (0.42, 0.34), "trotter_block": (0.60, 0.44),
    "ghz": (0.46, 0.34),
    "chiral_scurve": (0.50, 0.35),
    "sho": (0.42, 0.42), "qldpc": (0.50, 0.40),
    "dilfridge": (0.48, 0.66), "cryostat": (0.44, 0.50),
    "rack": (0.36, 0.54), "surface_code": (0.44, 0.44),
    "msd_factory": (0.56, 0.60), "lattice_surgery": (0.58, 0.28),
    "logical_ket": (0.42, 0.34), "color_code": (0.40, 0.40),
}

# INVERTED top: QEC substrate lowest, application glyphs above it, and
# the flagship equations as the crown -- the goal, not the machinery.
SET_PIECES = [
    dict(rows=[(10, 3)], size=0.62, pitch=1.05,
         menu=["dilfridge", "cryostat", "rack", "hemt_amp"]),
    # QEC substrate row
    dict(rows=[(12.9, 1)], size=0.40, pitch=1.12,
         menu=["surface_code", "msd_factory", "lattice_surgery",
               "logical_ket", "color_code", "code_qudit"]),
    # application row
    dict(rows=[(13.9, 1)], size=0.42, pitch=1.30,
         menu=["pmns_matrix", "chiral_scurve", "molecule",
               "wilson_loop", "u_nu", "ghz"]),
    # the crown: sized so two plaques fit INSIDE each tower's width
    dict(rows=[(14.8, 1.1)], size=0.43, pitch=1.42,
         menu=["qcd_ham", "kubo", "nu_osc", "shor_eq"]),
]

GLYPHS = {
    "transmon": cg.transmon,
    "jjunction": cg.jjunction,
    "gateH": lambda ax, xy, s, c, lw, accent=None: cg.gate1q(ax, xy, s, c, lw, accent, "H"),
    "gateX": lambda ax, xy, s, c, lw, accent=None: cg.gate1q(ax, xy, s, c, lw, None, "X"),
    "gateT": lambda ax, xy, s, c, lw, accent=None: cg.gate1q(ax, xy, s, c, lw, None, "T"),
    "gateRz": lambda ax, xy, s, c, lw, accent=None: cg.gate1q(ax, xy, s, c, lw, None, "Z"),
    "cnot": cg.cnot,
    "cz": cg.cz,
    "swap": cg.swapg,
    "toffoli": cg.toffoli,
    "meter": cg.meter,
    "bloch": cg.bloch,
    "ion_chain": cg.ion_chain,
    "chip5": cg.chip5,
    "resonator": cg.resonator,
    "circuit3": cg.circuit3,
    "heavy_hex": cg.heavy_hex,
    "srf_cavity": cg.srf_cavity,
    "tweezer": cg.tweezer_array,
    "qpu_die": cg.qpu_die,
    "group_circuit": cg.group_circuit,
    "atom": cg.atom,
    "molecule": cg.molecule,
    "circuit5": cg.circuit5,
    "pulse_seq": cg.pulse_seq,
    "chiral_scurve": cg.chiral_scurve,
    "sho": cg.sho,
    "qldpc": cg.qldpc,
    "kubo": cg.kubo,
    "code713": lambda ax, xy, s, c, lw, accent=None: cg.code_params(ax, xy, s, c, lw, accent, 0),
    "code_gross": lambda ax, xy, s, c, lw, accent=None: cg.code_params(ax, xy, s, c, lw, accent, 2),
    "code_iceberg": lambda ax, xy, s, c, lw, accent=None: cg.code_params(ax, xy, s, c, lw, accent, 3),
    "code513": lambda ax, xy, s, c, lw, accent=None: cg.code_params(ax, xy, s, c, lw, accent, 1),
    "code_shor9": lambda ax, xy, s, c, lw, accent=None: cg.code_params(ax, xy, s, c, lw, accent, 4),
    "code_rm15": lambda ax, xy, s, c, lw, accent=None: cg.code_params(ax, xy, s, c, lw, accent, 5),
    "code_golay": lambda ax, xy, s, c, lw, accent=None: cg.code_params(ax, xy, s, c, lw, accent, 6),
    "code_qudit": lambda ax, xy, s, c, lw, accent=None: cg.code_params(ax, xy, s, c, lw, accent, 7),
    "code_nkdq": lambda ax, xy, s, c, lw, accent=None: cg.code_params(ax, xy, s, c, lw, accent, 8),
    "wilson_loop": cg.wilson_loop,
    "trotter_block": cg.trotter_block,
    "ghz": cg.ghz,
    "peps": cg.peps,
    "mpo": cg.mpo,
    "nu_e": lambda ax, xy, s, c, lw, accent=None: cg.nu_ket(ax, xy, s, c, lw, accent, 0),
    "nu_mu": lambda ax, xy, s, c, lw, accent=None: cg.nu_ket(ax, xy, s, c, lw, accent, 1),
    "nu_tau": lambda ax, xy, s, c, lw, accent=None: cg.nu_ket(ax, xy, s, c, lw, accent, 2),
    "nu_osc": cg.nu_osc,
    "u_nu": cg.u_nu,
    "pmns_matrix": cg.pmns_matrix,
    "tn_loop": cg.tn_loop,
    "mera": cg.mera,
    "braket_bell": lambda ax, xy, s, c, lw, accent=None: cg.braket(ax, xy, s, c, lw, accent, seed=1),
    "braket": cg.braket,
    "shor_eq": cg.shor_eq,
    "qcd_ham": cg.qcd_ham,
    "wirebond_chip": cg.wirebond_chip,
    "hemt_amp": cg.hemt_amp,
    "circulator": cg.circulator,
    "laser": cg.laser,
    "color_code": cg.color_code,
    "squid": cg.squid,
    "spin_qubit": cg.spin_qubit,
    "blade_trap": cg.blade_trap,
    "cold_atoms": cg.cold_atoms,
    "tensor_net": cg.tensor_net,
    "plaquette": cg.plaquette,
    "braid": cg.braid,
    "dilfridge": cg.dilfridge,
    "cryostat": cg.cryostat,
    "rack": cg.rack,
    "surface_code": cg.surface_code,
    "msd_factory": cg.msd_factory,
    "lattice_surgery": cg.lattice_surgery,
    "logical_ket": cg.logical_ket,
}

ACCENTED = {"dilfridge", "lattice_surgery", "logical_ket", "msd_factory",
            "bloch", "spin_qubit", "blade_trap", "laser",
            "gateH", "group_circuit", "qcd_ham", "shor_eq", "kubo",
            "nu_osc", "pmns_matrix"}
ACCENT_EVERY = 3


# ------------------------------------------------------- construction -----

def building_paths():
    t = np.linspace(0, 1, 200)
    xo = x_out(t)
    y = GROUND + t * H
    outer = ([(CX - xo[i], y[i]) for i in range(len(t))]
             + [(CX + xo[i], y[i]) for i in range(len(t) - 1, -1, -1)])
    ts = np.linspace(0, 1, 140)
    ss = s_slot(ts)
    ys = GROUND + ts * H
    hole = ([(CX - ss[i], ys[i]) for i in range(len(ts))]
            + [(CX + ss[i], ys[i]) for i in range(len(ts) - 1, -1, -1)])
    verts = outer + [outer[0]] + hole[::-1] + [hole[::-1][0]]
    codes = ([Path.MOVETO] + [Path.LINETO] * (len(outer) - 1) + [Path.CLOSEPOLY]
             + [Path.MOVETO] + [Path.LINETO] * (len(hole) - 1) + [Path.CLOSEPOLY])
    return Path(verts, codes)


def draw_building(ax):
    t = np.linspace(0, 1, 200)
    y = GROUND + t * H
    xo = x_out(t)
    roof = GROUND + H
    ax.plot(CX - xo, y, color=WHITE, lw=1.5, zorder=5)
    ax.plot(CX + xo, y, color=WHITE, lw=1.5, zorder=5)
    # separate rooflines: the towers do NOT join across the gap
    for sgn in (-1, 1):
        ax.plot([CX + sgn * x_out(1), CX + sgn * SMIN], [roof, roof],
                color=WHITE, lw=1.5, zorder=5)
    # parapets, measured: inner edge VERTICAL at the gap, only the outer
    # side slopes (base wider than top), flat top
    from matplotlib.patches import Polygon as _Poly
    for sgn in (-1, 1):
        x_in = SMIN
        pts = [(CX + sgn * x_in, roof),
               (CX + sgn * (x_in + NUB_WB), roof),
               (CX + sgn * (x_in + NUB_WT), roof + NUB_H),
               (CX + sgn * x_in, roof + NUB_H)]
        ax.add_patch(_Poly(pts, closed=True, facecolor=PLATE,
                           edgecolor=WHITE, lw=1.2, zorder=6))
    # Glyphs set INTO the parapets, like the collage inside the towers: the
    # two ends of the report's story. Left nub holds the logical qubit (the
    # substrate every chapter runs on); right nub holds the SRF cavity (the
    # Fermilab hardware that delivers it).
    nub_cx = SMIN + NUB_WT / 2
    cg.logical_ket(ax, (CX - nub_cx, GROUND + H + 0.52 * NUB_H), 0.48,
                   WHITE, 0.7, accent=ORANGE)
    cg.srf_cavity(ax, (CX + nub_cx, GROUND + H + 0.50 * NUB_H), 0.40,
                  WHITE, 0.7, accent=ORANGE, cells=4)
    # inner edges run the full height: glass below T2, open gap above
    ts = np.linspace(0, 1, 140)
    ys = GROUND + ts * H
    ss = s_slot(ts)
    ax.plot(CX - ss, ys, color=WHITE, lw=1.2, zorder=5)
    ax.plot(CX + ss, ys, color=WHITE, lw=1.2, zorder=5)
    # top of the glass: a light sill line, not a structural cap
    ax.plot([CX - SMIN, CX + SMIN], [GROUND + T2 * H] * 2, color=PALE,
            lw=0.7, alpha=0.8, zorder=5)
    # glass: three mullions -> FOUR panes per row, persisting to the sill,
    # with a horizontal every floor plus half-floor lines in the wide wedge
    tg = np.linspace(0, T2, 120)
    yg = GROUND + tg * H
    sg = s_slot(tg)
    for f in (-0.5, 0.0, 0.5):
        ax.plot(CX + f * sg, yg, color=PALE, lw=0.45, alpha=0.6, zorder=4)
    for k in range(1, NFLOOR):
        tt = k / NFLOOR
        if tt >= T2:
            break
        yy = GROUND + k * FLOOR_H
        sk = float(s_slot(tt))
        ax.plot([CX - sk, CX + sk], [yy, yy], color=PALE, lw=0.4, alpha=0.55,
                zorder=4)
    for k in range(1, int(T1 * NFLOOR * 2) + 1, 2):
        tt = k / (2 * NFLOOR)
        yy = GROUND + tt * H
        sk = float(s_slot(tt))
        ax.plot([CX - sk, CX + sk], [yy, yy], color=PALE, lw=0.3, alpha=0.35,
                zorder=4)
    # floor lines
    for k in range(1, NFLOOR):
        tt = k / NFLOOR
        yy = GROUND + k * FLOOR_H
        xok, sk = float(x_out(tt)), float(s_slot(tt))
        if sk > 0.02:
            ax.plot([CX - xok, CX - sk], [yy, yy], color=PALE, lw=0.4,
                    alpha=0.40, zorder=3)
            ax.plot([CX + sk, CX + xok], [yy, yy], color=PALE, lw=0.4,
                    alpha=0.40, zorder=3)
        else:
            ax.plot([CX - xok, CX + xok], [yy, yy], color=PALE, lw=0.4,
                    alpha=0.40, zorder=3)
    # crossover bridges, inside the glass
    for k in CROSSOVERS:
        tt = k / NFLOOR
        yy = GROUND + k * FLOOR_H
        sk = float(s_slot(tt))
        ax.add_patch(Rectangle((CX - sk, yy), 2 * sk, FLOOR_H * 0.55,
                               facecolor=WHITE, alpha=0.10,
                               edgecolor=WHITE, lw=0.8, zorder=5))


def _menu_at(floor):
    names, ws = [], []
    for name, (mu, sig, amp) in WEIGHTS.items():
        z = abs(floor - mu) / sig
        if z <= 2.5:
            names.append(name)
            ws.append(amp * np.exp(-0.5 * z * z))
    w = np.array(ws)
    return names, w / w.sum()


def fill_tiers(ax, clip_patch):
    rng = np.random.default_rng(RNG_SEED)
    accent_count = 0
    lines0, patches0, texts0 = len(ax.lines), len(ax.patches), len(ax.texts)
    accepted = []

    def record(x, y, sz, name, accent=None):
        ew, eh = EXTENT.get(name, (0.5, 0.5))
        accepted.append((x, y, ew * sz, eh * sz, name))
        GLYPHS[name](ax, (x, y), sz, WHITE, 0.6, accent=accent)

    def try_place(x, y, sz, name, accent=None):
        ew, eh = EXTENT.get(name, (0.5, 0.5))
        hw, hh = ew * sz, eh * sz
        for (px, py, pw, ph, pn) in accepted:
            if abs(x - px) < hw + pw + GAP and abs(y - py) < hh + ph + GAP:
                return False
            if pn == name and np.hypot(x - px, y - py) < SAME_MIN_MAP.get(
                    name, SAME_MIN):
                return False
        record(x, y, sz, name, accent)
        return True

    # 1. deterministic set pieces, placed first so the darts avoid them
    for piece in SET_PIECES:
        size, pitch = piece["size"], piece["pitch"]
        menu = piece["menu"]
        ix = 0
        for (bf, span) in piece["rows"]:
            t_mid = (bf + span / 2) / NFLOOR
            y_mid = GROUND + (bf + span / 2) * FLOOR_H
            xok, sk = float(x_out(t_mid)), float(s_slot(t_mid))
            pad = 0.06 if sk > 0.02 else 0.10
            for (xa, xb) in ((CX - xok + 0.03, CX - sk - pad),
                             (CX + sk + pad, CX + xok - 0.03)):
                width = xb - xa
                step = size * pitch
                n = max(1, int(width / step + 0.30))
                margin = (width - n * step) / 2 + step / 2
                for i in range(n):
                    if ix >= len(menu):
                        break               # unique-only: never repeat
                    name = menu[ix]
                    ix += 1
                    record(xa + margin + i * step, y_mid + 0.02, size, name,
                           ORANGE if name in ACCENTED else None)

    # 2. guaranteed content, placed before the random darts
    for (fx, fy, fs, nm, acc) in (
            (-0.60, 8.30, 0.40, "srf_cavity", None),
            (+0.55, 8.80, 0.40, "srf_cavity", None),
            (+0.60, 7.55, 0.38, "group_circuit", ORANGE),
            (-0.55, 3.10, 0.26, "laser", None)):
        placed_seed = False
        for dy in (0.0, 0.18, -0.18, 0.36, -0.36, 0.55):
            for fxx in (fx, np.sign(fx) * min(0.92, abs(fx) + 0.25),
                        np.sign(fx) * max(0.10, abs(fx) - 0.25)):
                fy_ = fy + dy / FLOOR_H * 0  # dy applied in inches below
                t_ = fy_ / NFLOOR
                xok, sk = float(x_out(t_)), float(s_slot(t_))
                xr = CX + np.sign(fxx) * (sk + 0.10
                                          + abs(fxx) * (xok - sk - 0.30))
                if try_place(xr, GROUND + fy_ * FLOOR_H + dy, fs, nm, acc):
                    placed_seed = True
                    break
            if placed_seed:
                break
        if not placed_seed:
            print(f"WARNING: seeded glyph {nm!r} could not be placed")

    # 3. organic scatter over the whole mosaic, flowing around everything
    f_lo, f_hi = SCATTER_FLOORS
    tries = 0
    placed = len(accepted)
    while placed < SCATTER_N and tries < SCATTER_N * 140:
        tries += 1
        f = rng.uniform(f_lo, f_hi)
        t_ = f / NFLOOR
        y = GROUND + f * FLOOR_H
        frac = (f - f_lo) / (f_hi - f_lo)
        sz = (SIZE_LO + (SIZE_HI - SIZE_LO) * frac) * rng.uniform(*SIZE_JITTER)
        xok, sk = float(x_out(t_)), float(s_slot(t_))
        names, probs = _menu_at(f)
        name = names[int(rng.choice(len(names), p=probs))]
        sz *= SIZE_MULT.get(name, 1.0)
        ew = EXTENT.get(name, (0.5, 0.5))[0]
        xa = sk + 0.05 + ew * sz
        xb = xok - 0.03 - ew * sz
        if xb <= xa:
            continue
        x = CX + rng.choice([-1.0, 1.0]) * rng.uniform(xa, xb)
        accent = None
        if name in ACCENTED:
            accent_count += 1
            accent = ORANGE if accent_count % ACCENT_EVERY == 0 else None
        if try_place(x, y, sz, name, accent):
            placed += 1

    # fine fill: small circuits and notation into the remaining upper gaps
    f_lo2, f_hi2 = FILLER["floors"]
    menu2 = FILLER["menu"]
    placed2, tries2 = 0, 0
    while placed2 < FILLER["n"] and tries2 < FILLER["n"] * 110:
        tries2 += 1
        f = rng.uniform(f_lo2, f_hi2)
        t_ = f / NFLOOR
        y = GROUND + f * FLOOR_H
        sz = rng.uniform(*FILLER["size"])
        name = menu2[int(rng.integers(0, len(menu2)))]
        sz *= SIZE_MULT.get(name, 1.0)
        xok, sk = float(x_out(t_)), float(s_slot(t_))
        ew = EXTENT.get(name, (0.5, 0.5))[0]
        xa = sk + 0.05 + ew * sz
        xb = xok - 0.03 - ew * sz
        if xb <= xa:
            continue
        x = CX + rng.choice([-1.0, 1.0]) * rng.uniform(xa, xb)
        accent = None
        if name in ACCENTED:
            accent_count += 1
            accent = ORANGE if accent_count % ACCENT_EVERY == 0 else None
        if try_place(x, y, sz, name, accent):
            placed2 += 1

    for coll, n0 in ((ax.lines, lines0), (ax.patches, patches0),
                     (ax.texts, texts0)):
        for art in list(coll)[n0:]:
            art.set_clip_path(clip_patch)


# ------------------------------------------------------------ scenery -----

def shower(ax, xy, n_straight, n_curved, rmax, seed, ang0, spread, alpha):
    """Parton shower in the classic cascade idiom: gently curved flowing
    branches (photons wavy, gluons coiled), terminating in cluster blobs
    that hadronise into short-legged dots."""
    rng = np.random.default_rng(seed)
    x0, y0 = xy
    ax.add_patch(Circle((x0, y0), 0.045, facecolor=PALE, edgecolor="none",
                        alpha=alpha + 0.15, zorder=2))
    jet = np.deg2rad(ang0 + spread * rng.uniform(0.25, 0.45))

    def curved(px, py, a, L, lw, bow):
        """A bowed segment; returns its endpoint and exit angle."""
        ex, ey = px + L * np.cos(a), py + L * np.sin(a)
        nx, ny = -np.sin(a), np.cos(a)
        t = np.linspace(0, 1, 40)
        xs = px + t * (ex - px) + np.sin(np.pi * t) * bow * L * nx
        ys = py + t * (ey - py) + np.sin(np.pi * t) * bow * L * ny
        ax.plot(xs, ys, color=PALE, lw=lw, alpha=alpha, zorder=2,
                solid_capstyle="round")
        return ex, ey

    def cluster(px, py, a):
        """Terminal hadronisation: a blob, then dotted decay legs."""
        ax.add_patch(Ellipse((px, py), 0.075, 0.042,
                             angle=np.degrees(a), facecolor=PALE,
                             edgecolor="none", alpha=alpha + 0.10, zorder=2))
        for _ in range(int(rng.integers(2, 4))):
            aa = a + rng.normal(0, 0.55)
            LL = rng.uniform(0.07, 0.15)
            hx, hy = px + LL * np.cos(aa), py + LL * np.sin(aa)
            ax.plot([px, hx], [py, hy], color=PALE, lw=0.4,
                    alpha=alpha, zorder=2)
            ax.plot([hx], [hy], "o", ms=1.5, color=PALE,
                    alpha=alpha + 0.25, zorder=2)

    def track(px, py, a, L, lw, depth, bow_sign):
        u = rng.uniform()
        bow = bow_sign * rng.uniform(0.10, 0.24)
        if u < 0.62:
            ex, ey = curved(px, py, a, L, lw, bow)
        elif u < 0.80:
            ex, ey = px + L * np.cos(a), py + L * np.sin(a)
            ct.wavy(ax, (px, py), (ex, ey), n=max(3, int(L * 5)),
                    amp=0.020 * 0.65 ** depth, color=PALE, lw=lw * 0.85,
                    z=2, alpha=alpha)
        else:
            ex, ey = px + L * np.cos(a), py + L * np.sin(a)
            ct.coil(ax, (px, py), (ex, ey), n=max(3, int(L * 4.5)),
                    amp=0.017 * 0.65 ** depth, color=PALE, lw=lw * 0.8,
                    z=2, alpha=alpha)
        if depth < 2 and L > 0.45 and rng.uniform() < 0.6:
            nb = int(rng.integers(2, 4))
            for _ in range(nb):
                track(ex, ey, a + rng.normal(0, 0.28) + 0.25 * bow,
                      L * rng.uniform(0.35, 0.55), lw * 0.8, depth + 1,
                      bow_sign * rng.choice([1, 1, -1]))
        elif depth >= 1 and rng.uniform() < 0.55:
            cluster(ex, ey, a)

    for i in range(n_straight):
        if rng.uniform() < 0.7:
            a = jet + rng.normal(0, 0.18)
        else:
            a = np.deg2rad(ang0 + spread * rng.uniform(0, 1))
        L = rng.uniform(0.45, 1.0) * rmax * (1.15 if abs(a - jet) < 0.2 else 0.8)
        track(x0 + 0.07 * np.cos(a), y0 + 0.07 * np.sin(a), a, L, 0.55, 0,
              rng.choice([-1, 1]))

    for i in range(n_curved):
        a0 = np.deg2rad(ang0 + spread * rng.uniform(0.1, 0.9))
        L = rng.uniform(0.5, 0.85) * rmax
        curl = rng.choice([-1, 1]) * rng.uniform(1.2, 2.6)
        tt = np.linspace(0, 1, 60)
        aa = a0 + curl * tt ** 1.5
        rr = 0.07 + L * tt
        ax.plot(x0 + rr * np.cos(aa), y0 + rr * np.sin(aa), color=PALE,
                lw=0.55, alpha=alpha, zorder=2, solid_capstyle="round")


def vignette_nu_nucleus(ax, xy, s, alpha=0.45):
    x, y = xy
    # incoming neutrino, dashed, ending ON the cluster edge
    ax.plot([x - 0.85 * s, x - 0.13 * s], [y + 0.30 * s, y + 0.05 * s],
            color=PALE, lw=0.55, ls=(0, (2, 2)), alpha=alpha, zorder=2)
    # area-uniform blue-noise in a disc: reads as a round nucleus, not a
    # center-clumped smear
    rng = np.random.default_rng(6)
    pts, tries = [], 0
    while len(pts) < 11 and tries < 400:
        tries += 1
        a = rng.uniform(0, 2 * np.pi)
        r = 0.16 * s * np.sqrt(rng.uniform())
        px_, py_ = x + r * np.cos(a), y + r * np.sin(a)
        if all(np.hypot(px_ - q[0], py_ - q[1]) > 0.058 * s for q in pts):
            pts.append((px_, py_))
    for (px_, py_) in pts:
        ax.add_patch(Circle((px_, py_), 0.040 * s, facecolor=PALE,
                            edgecolor="none", alpha=alpha + 0.1, zorder=2))
    # outgoing lepton, leaving FROM the cluster edge
    ct.arrow(ax, (x + 0.11 * s, y + 0.09 * s),
             (x + 0.72 * s, y + 0.44 * s), color=PALE, lw=0.6, mut=6, z=2,
             alpha=alpha + 0.1)
    # boson radiated FROM the cluster edge, heading away
    ct.wavy(ax, (x + 0.06 * s, y - 0.13 * s), (x + 0.30 * s, y - 0.52 * s),
            n=3, amp=0.035 * s, color=PALE, lw=0.5, z=2)


def vignette_bubble(ax, xy, s, alpha=0.42):
    x, y = xy
    ax.add_patch(Circle((x, y), 0.30 * s, facecolor="none", edgecolor=PALE,
                        lw=0.6, alpha=alpha + 0.1, zorder=2))
    for a in np.deg2rad([30, 120, 210, 300]):
        ct.arrow(ax, (x + 0.32 * s * np.cos(a), y + 0.32 * s * np.sin(a)),
                 (x + 0.52 * s * np.cos(a), y + 0.52 * s * np.sin(a)),
                 color=PALE, lw=0.55, mut=5, z=2, alpha=alpha)
    rng = np.random.default_rng(11)
    for _ in range(9):
        a, r = rng.uniform(0, 2 * np.pi), rng.uniform(0.62, 0.95) * s
        ax.plot([x + r * np.cos(a)], [y + r * np.sin(a)], "o", ms=1.4,
                color=PALE, alpha=alpha, zorder=2)


def vignette_phase(ax, xy, s, alpha=0.42):
    """QCD phase diagram, drawn like a working sketch: crossover band that
    narrows into the critical endpoint, first-order line beyond it, chemical
    freeze-out curve with BES points, nuclear-matter point, CSC wedge."""
    x, y = xy
    # axes with arrowheads and ticks
    ax.annotate("", xy=(x + 1.02 * s, y), xytext=(x, y),
                arrowprops=dict(arrowstyle="-|>", color=PALE, lw=0.55,
                                mutation_scale=5, alpha=alpha), zorder=2)
    ax.annotate("", xy=(x, y + 0.78 * s), xytext=(x, y),
                arrowprops=dict(arrowstyle="-|>", color=PALE, lw=0.55,
                                mutation_scale=5, alpha=alpha), zorder=2)
    ax.text(x - 0.055 * s, y + 0.72 * s, r"$T$", fontsize=4.6, color=PALE,
            alpha=alpha + 0.15, ha="right", va="center", zorder=2)
    ax.text(x + 0.99 * s, y - 0.075 * s, r"$\mu_B$", fontsize=4.6, color=PALE,
            alpha=alpha + 0.15, ha="center", va="top", zorder=2)
    for tx in (0.25, 0.50, 0.75):
        ax.plot([x + tx * s] * 2, [y, y + 0.022 * s], color=PALE, lw=0.4,
                alpha=alpha, zorder=2)
    for ty in (0.26, 0.52):
        ax.plot([x, x + 0.022 * s], [y + ty * s] * 2, color=PALE, lw=0.4,
                alpha=alpha, zorder=2)
    # crossover band (two converging dashes) narrowing into the CEP
    tt = np.linspace(0, 1, 60)
    xc, yc = 0.46, 0.475                      # critical endpoint
    curve = 0.55 * s * np.sqrt(np.clip(1 - 0.72 * tt ** 2, 0, None))
    band = 0.045 * s * (1 - tt)
    mask = tt <= xc / 0.80
    for sgn in (-1, 1):
        ax.plot(x + tt[mask] * 0.80 * s, y + (curve + sgn * band)[mask],
                color=PALE, lw=0.5, ls=(0, (2.6, 1.6)), alpha=alpha + 0.05,
                zorder=2)
    # first-order line beyond the CEP, solid, heavier
    mask2 = tt >= xc / 0.80
    ax.plot(x + tt[mask2] * 0.80 * s, y + curve[mask2], color=PALE, lw=0.7,
            alpha=alpha + 0.15, zorder=2)
    # critical endpoint: orange dot with an uncertainty cross
    ax.plot([x + xc * s], [y + yc * s], "o", ms=2.4, color=GLOW,
            alpha=alpha + 0.30, zorder=3)
    ax.plot([x + (xc - 0.05) * s, x + (xc + 0.05) * s],
            [y + yc * s] * 2, color=GLOW, lw=0.45, alpha=alpha + 0.10, zorder=2)
    ax.plot([x + xc * s] * 2,
            [y + (yc - 0.035) * s, y + (yc + 0.035) * s], color=GLOW, lw=0.45,
            alpha=alpha + 0.10, zorder=2)
    # chemical freeze-out: dotted curve under the transition, BES points on it
    fo = 0.44 * s * np.sqrt(np.clip(1 - 0.92 * tt ** 2, 0, None))
    ax.plot(x + tt * 0.86 * s, y + fo, color=PALE, lw=0.45,
            ls=(0, (1.2, 1.8)), alpha=alpha, zorder=2)
    for tb in (0.10, 0.28, 0.46, 0.62, 0.76):
        fb = 0.44 * np.sqrt(max(1 - 0.92 * tb ** 2, 0))
        ax.plot([x + tb * 0.86 * s], [y + fb * s], "+", ms=2.6, mew=0.5,
                color=PALE, alpha=alpha + 0.10, zorder=2)
    # nuclear matter point on the mu axis + liquid-gas stub
    ax.plot([x + 0.68 * s], [y + 0.012 * s], "o", ms=1.6, color=PALE,
            mfc="none", mew=0.5, alpha=alpha + 0.10, zorder=2)
    ax.plot([x + 0.60 * s, x + 0.68 * s], [y + 0.055 * s, y + 0.012 * s],
            color=PALE, lw=0.4, alpha=alpha, zorder=2)
    # color-superconductor wedge, hatched
    for k in range(4):
        ax.plot([x + (0.80 + 0.045 * k) * s, x + (0.88 + 0.045 * k) * s],
                [y, y + 0.10 * s], color=PALE, lw=0.4, alpha=alpha * 0.9,
                zorder=2)
    # region labels
    ax.text(x + 0.16 * s, y + 0.62 * s, "QGP", fontsize=4.2, color=PALE,
            alpha=alpha + 0.12, ha="left", zorder=2)
    ax.text(x + 0.13 * s, y + 0.16 * s, "hadron gas", fontsize=3.8,
            color=PALE, alpha=alpha + 0.08, ha="left", zorder=2)
    ax.text(x + 0.955 * s, y + 0.135 * s, "CSC", fontsize=3.6, color=PALE,
            alpha=alpha + 0.05, ha="right", zorder=2)


def draw_sky(ax):
    # Four detailed HEP vignettes replace the corner showers and the two small
    # mid-sky doodles. The corner pieces overhang the frame deliberately --
    # only a rotated slice of each shows, like the showers they replace -- so
    # every artist they add is clipped to the inner frame.
    clip = Rectangle((FRAME2, FRAME2), PW - 2 * FRAME2, PH - 2 * FRAME2,
                     facecolor="none", edgecolor="none")
    ax.add_patch(clip)

    def clipped(drawfn, **kw):
        n = (len(ax.lines), len(ax.patches), len(ax.collections), len(ax.texts))
        drawfn(ax, **kw)
        for coll, n0 in zip((ax.lines, ax.patches, ax.collections, ax.texts), n):
            for art in list(coll)[n0:]:
                art.set_clip_path(clip)

    # Detector slice: design center pinned to the inner-frame corner, wedge
    # axis (design +90 deg) aimed at the towers, so the fan sweeps from the
    # corner toward Wilson Hall and only the inner quadrant shows.
    # Six mid-sky vignettes distributed around the lab: two per side strip,
    # two in the band above the towers, sized to fill their slots.
    aim = np.degrees(np.arctan2((GROUND + H - 0.6) - (PH - FRAME2),
                                CX - FRAME2))
    clipped(draw_detector_slice, center=(0.02, PH + 0.02),
            rot_deg=aim - 93, scale=2.30, alpha=0.62)
    clipped(draw_bubble_chamber, center=(8.02, 10.48), rot_deg=28,
            scale=1.30, alpha=0.62)
    # band above the towers, pushed up and outward from the parapets
    draw_nu_nucleus(ax, center=(3.08, 8.72), rot_deg=0,
                    scale=0.95, alpha=0.72)
    draw_loop_cloud(ax, center=(5.95, 8.62), rot_deg=12,
                    scale=0.95, alpha=0.72)
    # left strip
    draw_lattice_qcd(ax, center=(1.25, GROUND + 3.65), rot_deg=8,
                     scale=1.25, alpha=0.80)
    vignette_phase(ax, (0.35, GROUND + 0.90), 1.55)
    # right strip: event display raised, cosmology up into its old slot,
    # clear of the obelisk
    draw_event_display(ax, center=(PW - 1.25, GROUND + 5.35), rot_deg=-12,
                       scale=0.95, alpha=0.78)
    draw_cosmology(ax, center=(PW - 1.40, GROUND + 3.97), rot_deg=-8,
                   scale=1.08, alpha=0.80)
    ax.plot([FRAME2 + 0.15, PW - FRAME2 - 0.15],
            [GROUND + 4.95, GROUND + 5.45],
            color=PALE, lw=0.5, ls=(0, (1.5, 2.5)), alpha=0.26, zorder=1)


def draw_ground_and_pond(ax):
    ax.plot([FRAME2, PW - FRAME2], [GROUND, GROUND], color=WHITE, lw=1.0,
            zorder=5)
    depth = (GROUND - POND0) / H * 0.85
    for a0, a1, t0, t1 in ((0.30, 0.18, 0.0, 0.45), (0.18, 0.06, 0.45, 1.0)):
        t = np.linspace(t0, t1, 60)
        y_r = GROUND - t * H * depth
        xo = x_out(t)
        wob = 0.022 * np.sin(y_r * 22)
        for sgn in (-1, 1):
            ax.plot(CX + sgn * xo + wob, y_r, color=PALE, lw=0.75,
                    alpha=(a0 + a1) / 2, zorder=2)
    for k in (3, 7, 11):
        tt = k / NFLOOR
        yy = GROUND - tt * H * depth
        xok = float(x_out(tt))
        ax.plot([CX - xok, CX + xok], [yy, yy], color=PALE, lw=0.35,
                alpha=0.14, zorder=2)
    ax.add_patch(Ellipse((CX, GROUND - H * depth - 0.16), 0.95, 0.16,
                         facecolor=GLOW, alpha=0.06, edgecolor="none",
                         zorder=2))
    ax.add_patch(Ellipse((CX, GROUND - H * depth - 0.16), 0.52, 0.08,
                         facecolor=GLOW, alpha=0.09, edgecolor="none",
                         zorder=2))
    for (rx, ry, sc) in ((1.60, GROUND - 0.52, 1.3), (7.32, GROUND - 0.80, 1.05)):
        for k in range(3):
            w = sc * (0.16 + 0.14 * k)
            ax.add_patch(Ellipse((rx, ry), 2 * w, 0.42 * w, facecolor="none",
                                 edgecolor=PALE, lw=0.45,
                                 alpha=0.30 - 0.08 * k, zorder=3))
    rng = np.random.default_rng(5)
    for _ in range(12):
        rx = rng.uniform(1.0, PW - 1.0)
        ry = rng.uniform(POND0 + 0.10, GROUND - 0.14)
        hw = rng.uniform(0.25, 0.75)
        ax.plot([rx - hw, rx + hw], [ry, ry], color=PLATE, lw=2.6, zorder=3)
        ax.plot([rx - hw * 0.8, rx + hw * 0.8], [ry + 0.016, ry + 0.016],
                color=PALE, lw=0.4, alpha=0.30, zorder=3)


def draw_shore(ax):
    rng = np.random.default_rng(21)
    for cxg in (0.75, 1.25, 6.95, 7.72):
        n = int(rng.integers(5, 8))
        for _ in range(n):
            dx = rng.uniform(-0.09, 0.09)
            h = rng.uniform(0.08, 0.18)
            lean = rng.uniform(-0.04, 0.04)
            ax.plot([cxg + dx, cxg + dx + lean], [GROUND, GROUND + h],
                    color=PALE, lw=0.5, alpha=0.65, zorder=4)
    for cxg, hh in ((1.00, 0.24), (7.05, 0.21)):
        ax.plot([cxg, cxg], [GROUND, GROUND + hh], color=PALE, lw=0.55,
                alpha=0.7, zorder=4)
        ax.add_patch(Ellipse((cxg, GROUND + hh + 0.03), 0.032, 0.075,
                             facecolor=PALE, edgecolor="none", alpha=0.7,
                             zorder=4))


def draw_obelisk(ax):
    ox = 7.32
    base, top = GROUND - 0.45, GROUND + 3.05
    t = np.linspace(0, 1, 90)
    yy = base + t * (top - base)
    w = 0.112 - 0.056 * np.exp(-((t - 0.55) / 0.34) ** 2)
    for sgn in (-1, 1):
        ax.plot(ox + sgn * w, yy, color=WHITE, lw=1.0, zorder=6)
    ax.plot([ox - w[-1], ox + w[-1]], [top, top], color=WHITE, lw=1.0,
            zorder=6)
    ax.plot([ox, ox], [base + 0.08, top - 0.05], color=PALE, lw=0.35,
            alpha=0.5, zorder=6)
    tr = np.linspace(0, 1, 50)
    yr = base - tr * (top - base) * 0.28
    wr = 0.112 - 0.056 * np.exp(-((tr - 0.55) / 0.34) ** 2)
    wob = 0.016 * np.sin(yr * 26)
    for sgn in (-1, 1):
        ax.plot(ox + sgn * wr + wob, yr, color=PALE, lw=0.65, alpha=0.25,
                zorder=2)
    for k in range(2):
        ww = 0.22 + 0.17 * k
        ax.add_patch(Ellipse((ox, base + 0.02), 2 * ww, 0.34 * ww,
                             facecolor="none", edgecolor=PALE, lw=0.45,
                             alpha=0.30 - 0.10 * k, zorder=3))


def draw_flagpoles(ax):
    rng = np.random.default_rng(13)
    for k in range(11):
        fx = CX - 1.85 + k * 0.37
        hh = 0.26 + rng.uniform(-0.012, 0.012)
        ax.plot([fx, fx], [GROUND, GROUND + hh], color=PALE, lw=0.4,
                alpha=0.55, zorder=6)
        ax.plot([fx, fx + 0.055], [GROUND + hh, GROUND + hh - 0.018],
                color=PALE, lw=0.55, alpha=0.5, zorder=6)


def draw_bison(ax):
    cg.bison(ax, (7.58, GROUND + 0.13), 0.36, WHITE, alpha=0.88)
    cg.bison(ax, (8.16, GROUND + 0.10), 0.23, WHITE, alpha=0.74)


def draw_plate(ax):
    ax.add_patch(Rectangle((0, 0), PW, PH, facecolor=PLATE, edgecolor="none",
                           zorder=0))
    for inset, lw in ((FRAME1, 1.3), (FRAME2, 0.6)):
        ax.add_patch(Rectangle((inset, inset), PW - 2 * inset, PH - 2 * inset,
                               facecolor="none", edgecolor=WHITE, lw=lw,
                               zorder=8))
    for gx in np.arange(FRAME2, PW - FRAME2, GRID):
        ax.plot([gx, gx], [FRAME2, PH - FRAME2], color=WHITE, lw=0.3,
                alpha=0.045, zorder=1)
    for gy in np.arange(FRAME2, PH - FRAME2, GRID):
        ax.plot([FRAME2, PW - FRAME2], [gy, gy], color=WHITE, lw=0.3,
                alpha=0.045, zorder=1)


def main():
    ff.setup()
    fig = plt.figure(figsize=(PW, PH), dpi=300)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, PW)
    ax.set_ylim(0, PH)
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)

    draw_plate(ax)
    draw_sky(ax)
    draw_ground_and_pond(ax)

    clip = PathPatch(building_paths(), facecolor="none", edgecolor="none",
                     zorder=1)
    ax.add_patch(clip)
    ax.add_patch(PathPatch(building_paths(), facecolor=WHITE, alpha=0.030,
                           edgecolor="none", zorder=2))

    fill_tiers(ax, clip)
    draw_building(ax)
    draw_shore(ax)
    draw_obelisk(ax)
    draw_flagpoles(ax)
    draw_bison(ax)

    out = ff.FIGDIR / "cover_wilson_hall.pdf"
    fig.savefig(out, facecolor=PLATE)
    print(f"wrote {out}")
    out_svg = ff.FIGDIR / "cover_wilson_hall.svg"
    fig.savefig(out_svg, facecolor=PLATE)
    print(f"wrote {out_svg}")


if __name__ == "__main__":
    main()
