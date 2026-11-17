#!/usr/bin/env python3
"""Ch. 4 cartoon: where the quantum device sits inside a lattice-QCD workflow.

Read left to right. A gauge configuration is generated classically on HPC, as
it always has been (drawn as a stack of stored configurations). From its links
one builds the Dirac matrix M -- sparse, but of order 10^4 on a side (per
flavor), and it is the determinant and the trace of its inverse that the Monte
Carlo needs. Those two numbers, and only those, are handed to a quantum
coprocessor that block-encodes M, applies a QSVT polynomial, and estimates the
trace.

The point of the figure is the size of the last panel: the device holds one
linear-algebra primitive, not the gauge field, which is why this chapter's
register is the smallest in the report while its depth is not.

v2 layout: drawing units are inches on a 6.5in-wide canvas, so font sizes in
points map directly onto the page. The figure is set at 0.95 textwidth, so
every label is >= 7.4pt here (>= 7pt printed).

Usage: fig_ch04_hybrid_lqcd.py [OUTDIR]
"""

import argparse
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import fnalfig as ff

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import (FancyArrowPatch, FancyBboxPatch, PathPatch,
                                Polygon, Rectangle)
from matplotlib.path import Path as MPath

# Fermilab identity (fnal-report.cls) plus one restrained accent.
NAVY = ff.FNAL_BLUE        # 0,33,71
LOGO = "#0061A1"           # 0,97,161
ACCENT = ff.ORANGE         # the two traces, as in the original
CL_FILL = "#F2F3F5"        # classical zone tint
CL_EDGE = "#C3C9D0"
SLATE = "#55616E"          # classical header band
Q_FILL = "#EAF2F9"         # quantum zone tint
Q_EDGE = "#9DBBD6"
ACCENT_LINK = LOGO

W, H = 6.5, 2.42            # canvas in inches (old figure: 6.48 x 2.71 in)
S_MAIN, S_SUB, S_BIG = 8.0, 7.4, 9.0


def T(ax, xy, s, size=S_SUB, color=ff.DARKGREY, ha="center", va="center",
      z=20, **kw):
    return ax.text(*xy, s, fontsize=size, color=color, ha=ha, va=va,
                   zorder=z, **kw)


def rbox(ax, x, y, w, h, fc, ec="none", lw=0.8, r=0.08, z=1, **kw):
    p = FancyBboxPatch((x + r, y + r), w - 2 * r, h - 2 * r,
                       boxstyle=f"round,pad={r}", fc=fc, ec=ec, lw=lw,
                       zorder=z, **kw)
    ax.add_patch(p)
    return p


def arrow(ax, p0, p1, color=ff.INK, lw=1.0, mut=8, z=8, rad=0.0):
    ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle="-|>", mutation_scale=mut,
                                 color=color, lw=lw, zorder=z,
                                 shrinkA=0, shrinkB=0,
                                 connectionstyle=f"arc3,rad={rad}",
                                 joinstyle="miter", capstyle="butt"))


# ----------------------------------------------------------------- glyphs --

def _cube_pts(origin, size, n, skew):
    sx, sy = skew
    ox, oy = origin
    d = size / (n - 1)
    return lambda i, j, k: np.array([ox + i * d + sx * k * d,
                                     oy + j * d + sy * k * d])


def _faces(P, n):
    m = n - 1
    front = [P(0, 0, 0), P(m, 0, 0), P(m, m, 0), P(0, m, 0)]
    top = [P(0, m, 0), P(m, m, 0), P(m, m, m), P(0, m, m)]
    right = [P(m, 0, 0), P(m, 0, m), P(m, m, m), P(m, m, 0)]
    hull = [P(0, 0, 0), P(m, 0, 0), P(m, 0, m), P(m, m, m), P(0, m, m),
            P(0, m, 0)]
    return front, top, right, hull


def cube_card(ax, origin, size, n, skew, fc, ec, z):
    """A stored configuration further back in the stack: silhouette only."""
    P = _cube_pts(origin, size, n, skew)
    front, top, right, hull = _faces(P, n)
    for poly, f in ((top, fc[0]), (right, fc[1]), (front, fc[2])):
        ax.add_patch(Polygon(poly, closed=True, fc=f, ec=ec, lw=0.5,
                             zorder=z, joinstyle="round"))


def lattice_cube(ax, origin, size, n=4, skew=(0.42, 0.30), z=4):
    """Oblique 3D lattice: tinted opaque faces, depth-faded links, sites."""
    P = _cube_pts(origin, size, n, skew)
    front, top, right, hull = _faces(P, n)
    ax.add_patch(Polygon(hull, closed=True, fc="white", ec="none", zorder=z))
    for poly, fc in ((top, "#DDE8F3"), (right, "#CCDCEC"), (front, "#F5F8FC")):
        ax.add_patch(Polygon(poly, closed=True, fc=fc, ec="none", zorder=z))
    for k in range(n):
        depth = 1.0 - 0.6 * k / (n - 1)  # back planes lighter
        for i in range(n):
            for j in range(n):
                for (a, b) in (((i, j, k), (i + 1, j, k)),
                               ((i, j, k), (i, j + 1, k)),
                               ((i, j, k), (i, j, k + 1))):
                    if max(b) > n - 1:
                        continue
                    pa, pb = P(*a), P(*b)
                    ax.plot([pa[0], pb[0]], [pa[1], pb[1]], color="#7D8B9B",
                            lw=0.45, alpha=depth, zorder=z + 1,
                            solid_capstyle="round")
    # one highlighted link U_l on the front face
    pa, pb = P(1, 1, 0), P(2, 1, 0)
    ax.plot([pa[0], pb[0]], [pa[1], pb[1]], color=ACCENT_LINK, lw=1.8,
            zorder=z + 2, solid_capstyle="round")
    for k in range(n):
        depth = 1.0 - 0.55 * k / (n - 1)
        pts = np.array([P(i, j, k) for i in range(n) for j in range(n)])
        ax.scatter(pts[:, 0], pts[:, 1], s=5.0 * (0.55 + 0.45 * depth),
                   color=NAVY, alpha=depth, lw=0, zorder=z + 3)
    return P

def sparse_matrix(ax, x, y, w, n=28):
    """Square matrix with a 4D-stencil-like sparsity pattern in brackets."""
    ax.add_patch(Rectangle((x, y), w, w, fc="white", ec="none", zorder=3))
    c = w / n
    offs = {0: 1.0, 1: 0.55, 3: 0.40, 9: 0.30}
    for i in range(n):
        for o, a in offs.items():
            for j in {i - o, i + o}:
                if 0 <= j < n:
                    ax.add_patch(Rectangle((x + j * c + 0.08 * c,
                                            y + (n - 1 - i) * c + 0.08 * c),
                                           0.84 * c, 0.84 * c,
                                           fc=NAVY if o == 0 else LOGO,
                                           alpha=0.85 if o == 0 else a,
                                           ec="none", zorder=4))
    t = 0.10 * w
    pad = 0.05
    for xx, s in ((x - pad, 1), (x + w + pad, -1)):
        ax.plot([xx + s * t, xx, xx, xx + s * t],
                [y - pad, y - pad, y + w + pad, y + w + pad],
                color=ff.INK, lw=1.0, zorder=5, solid_joinstyle="miter")


def pill(ax, xc, yc, w, h, text, color, size=8.6):
    rbox(ax, xc - w / 2, yc - h / 2, w, h, fc="white", ec=color, lw=0.9,
         r=h / 2 - 0.005, z=9)
    T(ax, (xc, yc), text, size=size, color=color)


# ------------------------------------------------------------------- main --

def main(outdir):
    ff.setup()
    plt.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42})
    fig = plt.figure(figsize=(W, H))
    ax = fig.add_axes([0, 0, 1, 1])
    ff.blank(ax)
    ax.set_xlim(0, W)
    ax.set_ylim(0, H)

    ymid = 1.24      # pipeline axis
    y_lo1, y_lo2 = 0.46, 0.27   # two rows of sub-labels under the pipeline

    # ---- (a) classical zone: the stored gauge configurations ----
    cx0, cx1 = 0.02, 1.46
    hb = 0.30  # header band height, shared by both zones
    rbox(ax, cx0, 0.04, cx1 - cx0, H - 0.06, fc=CL_FILL, ec=CL_EDGE, lw=0.8,
         r=0.10, z=0)
    rbox(ax, cx0, H - 0.02 - hb, cx1 - cx0, hb, fc=SLATE, r=0.10, z=1)
    ax.add_patch(Rectangle((cx0, H - 0.02 - hb), cx1 - cx0, 0.12, fc=SLATE,
                           ec="none", zorder=1))
    rbox(ax, cx0, 0.04, cx1 - cx0, H - 0.06, fc="none", ec=CL_EDGE, lw=0.8,
         r=0.10, z=2)  # outline over the header band
    T(ax, ((cx0 + cx1) / 2, H - 0.02 - hb / 2), "classical HMC on HPC",
      size=S_MAIN, color="white")
    size, n, skew = 0.62, 4, (0.40, 0.30)
    ox, oy = 0.15, 0.78
    cards = ((("#D3DEEA", "#C4D3E3", "#E3EAF2"), (0.30, 0.22)),
             (("#DFE8F1", "#D0DDEA", "#EDF2F7"), (0.15, 0.11)))
    for m, (fc, (dx, dy)) in enumerate(cards):
        cube_card(ax, (ox + dx, oy + dy), size, n, skew, fc=fc,
                  ec="#A3B5C8", z=2 + m)
    lattice_cube(ax, (ox, oy), size, n=n, skew=skew, z=5)
    T(ax, ((cx0 + cx1) / 2, y_lo1), r"gauge links $U_\ell$", size=S_MAIN,
      color=ff.INK)
    T(ax, ((cx0 + cx1) / 2, y_lo2),
      r"$V=8^3\!\times\!16$, $a\!\sim\!0.12\,$fm", size=S_SUB)

    # ---- (b) the Dirac matrix built from those links ----
    ax0, ax1 = cx1 + 0.06, cx1 + 0.31
    arrow(ax, (ax0, ymid), (ax1, ymid), color=ff.INK, lw=1.1)
    T(ax, ((ax0 + ax1) / 2, ymid + 0.14), "build", size=S_SUB)

    mw = 0.78
    mx, my = ax1 + 0.15, ymid - mw / 2
    mcx = mx + mw / 2
    sparse_matrix(ax, mx, my, mw)
    T(ax, (mcx, H - 0.02 - hb / 2), r"$M$  (Dirac operator)", size=S_MAIN,
      color=ff.INK)
    T(ax, (mcx, my + mw + 0.20), r"sparse, $2.5\!\times\!10^{4}$ a side",
      size=S_SUB)
    T(ax, (mcx, y_lo1 - 0.03), r"$\kappa\!\sim\!10^{2}$", size=S_SUB)

    # what the Monte Carlo actually needs from it
    pw, ph = 0.74, 0.25
    px = mx + mw + 0.22 + pw / 2
    y_ld, y_tr = ymid + 0.27, ymid - 0.27
    pill(ax, px, y_ld, pw, ph, r"$\log\det M$", ACCENT)
    pill(ax, px, y_tr, pw, ph, r"$\mathrm{Tr}\,M^{-1}$", ACCENT)
    x_m1 = mx + mw + 0.09
    for y0, y1 in ((ymid + 0.06, y_ld), (ymid - 0.06, y_tr)):
        arrow(ax, (x_m1, y0), (px - pw / 2 - 0.03, y1), color=ACCENT,
              lw=0.9, mut=7)
    T(ax, (px, y_lo1 - 0.03), "the only two", size=S_SUB)
    T(ax, (px, y_lo2), "numbers handed off", size=S_SUB)

    # ---- (c) quantum coprocessor ----
    qx0, qx1 = px + pw / 2 + 0.27, W - 0.02
    qy0, qy1 = 0.04, H - 0.02
    rbox(ax, qx0, qy0, qx1 - qx0, qy1 - qy0, fc=Q_FILL, ec=Q_EDGE, lw=0.8,
         r=0.10, z=0)
    rbox(ax, qx0, qy1 - hb, qx1 - qx0, hb, fc=NAVY, r=0.10, z=1)
    ax.add_patch(Rectangle((qx0, qy1 - hb), qx1 - qx0, 0.12, fc=NAVY,
                           ec="none", zorder=1))
    rbox(ax, qx0, qy0, qx1 - qx0, qy1 - qy0, fc="none", ec=Q_EDGE, lw=0.8,
         r=0.10, z=2)  # outline over the header band
    qc = (qx0 + qx1) / 2
    T(ax, (qc, qy1 - hb / 2), "quantum coprocessor", size=S_MAIN + 0.5,
      color="white")

    # the hand-off: both traces merge into one arrow onto the device
    xr, xm = px + pw / 2 + 0.03, qx0 - 0.10
    g = xm - xr
    for y0 in (y_ld, y_tr):
        ax.add_patch(PathPatch(MPath([(xr, y0), (xr + 0.55 * g, y0),
                                     (xm - 0.55 * g, ymid), (xm, ymid)],
                                    [MPath.MOVETO, MPath.CURVE4, MPath.CURVE4,
                                     MPath.CURVE4]),
                               fc="none", ec=ACCENT, lw=0.9, zorder=8,
                               capstyle="round"))
    arrow(ax, (xm - 0.005, ymid), (qx0 - 0.01, ymid), color=ACCENT, lw=0.9,
          mut=7)

    T(ax, (qc, 1.86), "83 logical qubits", size=S_BIG, color=LOGO,
      weight="bold")
    T(ax, (qc, 1.645), "13 site + 4 column + 4 color", size=S_SUB)
    T(ax, (qc, 1.485), "+ 12 QPE + 50 ancilla", size=S_SUB)

    # circuit: BE on all three wires, QSVT on the bottom two, read on top
    wx0, wx1 = qx0 + 0.32, qx1 - 0.08
    dy = 0.165
    wy = [0.92, 0.92 + dy, 0.92 + 2 * dy]  # |0>, |0>, |psi>, bottom to top
    for yy in wy:
        ax.plot([wx0, wx1], [yy, yy], color=ff.GREY, lw=0.8, zorder=10)
    for yy, lab in zip(wy, (r"$|0\rangle$", r"$|0\rangle$", r"$|\psi\rangle$")):
        T(ax, (wx0 - 0.04, yy), lab, size=S_SUB, ha="right")
    span = wx1 - wx0
    gx = [wx0 + 0.24, wx0 + 0.86, wx1 - 0.40]  # spaced to fit the labels
    gw = 0.30
    rbox(ax, gx[0] - gw / 2, wy[0] - 0.10, gw, wy[2] - wy[0] + 0.20,
         fc=LOGO, ec=LOGO, lw=0.8, r=0.03, z=12)
    T(ax, (gx[0], wy[1]), "BE", size=S_SUB, color="white", z=13)
    gw2 = 0.48
    rbox(ax, gx[1] - gw2 / 2, wy[0] - 0.09, gw2, wy[1] - wy[0] + 0.18,
         fc=NAVY, ec=NAVY, lw=0.8, r=0.03, z=12)
    T(ax, (gx[1], (wy[0] + wy[1]) / 2), "QSVT", size=S_SUB, color="white",
      z=13)
    rbox(ax, gx[2] - 0.22, wy[2] - 0.085, 0.44, 0.17, fc="white", ec=ff.INK,
         lw=0.8, r=0.03, z=12)
    T(ax, (gx[2], wy[2]), "read", size=S_SUB, color=ff.INK, z=13)

    for xg, lab in zip(gx, (r"encode $M$", r"apply $1/x$", "estimate trace")):
        T(ax, (xg, 0.66), lab, size=S_SUB)
    ax.plot([qx0 + 0.25, qx1 - 0.25], [0.555, 0.555], color=Q_EDGE, lw=0.6,
            zorder=3)
    T(ax, (qc, y_lo1 - 0.03), "one linear-algebra primitive,", size=S_SUB,
      color=NAVY)
    T(ax, (qc, y_lo2), "not the gauge field", size=S_SUB, color=NAVY)

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    out = outdir / "ch04_hybrid_lqcd.pdf"
    fig.savefig(out)
    print(f"wrote {out}")
    return fig, ax, out


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("outdir", nargs="?", default=str(ff.FIGDIR),
                    help="output directory (default: Figures/)")
    a = ap.parse_args()
    main(a.outdir)
