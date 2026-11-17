#!/usr/bin/env python3
"""Ch. 6 cartoon: the light-cone problem, and the two ways around it.

Three panels in one row at text width (layout since 2026-10-12, author decision):
(a) the featured route (inject a colour-singlet pair, evolve, count by species),
(b) the fragmentation instance (struck quark, shower, D(z)), (c) the light-cone
obstruction that applies to (b) only.
"""

import numpy as np
from matplotlib.patches import Circle, Ellipse

import cartoon as ct
import fnalfig as ff
import matplotlib.pyplot as plt


def caption(ax, s):
    """Panel caption, two lines, set on a common baseline below all three panels."""
    ax.annotate(s, xy=(0.5, 0.01), xycoords=("axes fraction", "figure fraction"),
                ha="center", va="bottom", fontsize=7, color=ff.INK,
                annotation_clip=False, linespacing=1.1)


def lightcone(ax):
    ff.blank(ax)
    ax.set_aspect("auto")
    ax.set_xlim(-1.15, 1.55)
    ax.set_ylim(-0.30, 1.65)

    # light cone
    ax.plot([-1.05, 1.05], [1.05, 1.05], color="none")
    for s in (-1, 1):
        ax.plot([0, s * 1.05], [0, 1.05], color=ff.GREY, lw=0.7, ls=(0, (3, 2)))
    ct.txt(ax, (0.80, 1.16), "light cone", size=7, color=ff.GREY)

    # the light-cone separation (what is wanted)
    ct.arrow(ax, (0, 0), (0.82, 0.82), color=ff.ORANGE, lw=1.3, style="-|>")
    ct.txt(ax, (0.55, 0.40), "light-cone\nseparation", size=7, color=ff.ORANGE, ha="left",
           bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none"))

    # the equal-time surrogate in a boosted frame (what is computed)
    ct.arrow(ax, (0, 0), (0.92, 0.0), color=ff.BLUE, lw=1.3, style="-|>")
    ct.txt(ax, (0.46, -0.19), r"equal-time, boost $P$", size=7, color=ff.BLUE)

    ax.plot([0], [0], "o", ms=3.2, color=ff.INK, zorder=6)

    # matching, and where it fails
    # Set below the left cone line, which owns the upper-left corner.
    bb = dict(boxstyle="round,pad=0.12", fc="white", ec="none")
    ct.txt(ax, (-1.10, 0.52), "LaMET matching", size=7, color=ff.INK,
           ha="left", bbox=bb)
    ct.txt(ax, (-1.10, 0.34), r"needs large $zP$", size=7,
           color=ff.MIDGREY, ha="left", bbox=bb)
    ct.txt(ax, (-1.10, 0.17), r"breaks down as $z\to 0$", size=7,
           color=ff.ORANGE, ha="left", bbox=bb)
    ct.arrow(ax, (-0.44, 0.60), (-0.10, 0.86), color=ff.MIDGREY, lw=0.7, mut=6)

    caption(ax, "(c)  light-cone obstruction,\n(b) only")


def fragmentation(ax):
    ff.blank(ax)
    ax.set_aspect("auto")
    ax.set_xlim(0, 6.0)
    ax.set_ylim(-0.55, 2.75)

    # A boosted pion streaming in: Lorentz-contracted along the boost, with
    # motion lines behind it, in the same idiom as the beam in Fig. 2.1.
    for k, dy in enumerate((0.20, 0.0, -0.20)):
        ax.plot([0.02, 0.42 - 0.06 * abs(k - 1)], [1.42 + dy] * 2,
                color=ff.GREY, lw=0.6, alpha=0.8, zorder=2,
                solid_capstyle="round")
    ax.add_patch(Ellipse((0.78, 1.42), 0.34, 0.62, facecolor=ff.BLUE,
                         alpha=0.14, edgecolor=ff.BLUE, lw=0.9, zorder=3))
    rng0 = np.random.default_rng(3)
    for a in (0.7, 2.4):
        ax.add_patch(Circle((0.78 + 0.09 * np.cos(a), 1.42 + 0.17 * np.sin(a)),
                            0.062, facecolor=ff.BLUE, edgecolor="white",
                            lw=0.4, zorder=4))
    ct.txt(ax, (0.78, 0.92), r"boosted $\pi$", size=7, color=ff.BLUE)
    ct.txt(ax, (0.60, 2.10), r"$P\!\sim\!2\,$GeV", size=7, color=ff.MIDGREY)

    # the current comes in and terminates on the pion; the struck quark leaves
    # from that point, so the three objects form one connected line
    ct.wavy(ax, (1.46, 2.14), (0.94, 1.62), n=3, amp=0.065, color=ff.ORANGE,
            lw=0.9)
    ct.txt(ax, (1.62, 2.20), r"$J^\mu$", size=7.4, color=ff.ORANGE)
    ct.arrow(ax, (0.98, 1.44), (2.30, 1.44), color=ff.INK, lw=1.0)
    ct.txt(ax, (1.88, 1.70), "struck quark", size=7)

    # a gluon radiated off the quark line itself, not floating beside it
    ct.coil(ax, (1.84, 1.42), (2.44, 0.96), n=4, amp=0.065, color=ff.ORANGE,
            lw=0.75)

    # the shower: three partons leaving the end of the quark line and reaching
    # the hadrons they become, rather than arrows stopping in open space
    for tip in ((3.22, 1.92), (3.30, 1.18), (3.14, 0.62)):
        ct.arrow(ax, (2.32, 1.44), tip, color=ff.GREY, lw=0.6, mut=5, z=2,
                 shrinkB=4)

    # final-state hadrons, each a little blob, sorted by momentum fraction
    zs = [(3.52, 1.98, 0.20), (3.62, 1.10, 0.26), (3.98, 1.62, 0.16),
          (4.06, 0.78, 0.13), (3.46, 0.44, 0.11)]
    for n, (hx, hy, hr) in enumerate(zs):
        ct.hadron(ax, (hx, hy), r=hr, quarks=2, color=ff.BLUE, seed=n + 1,
                  z=4, alpha=0.16, lw=0.7)
    ct.txt(ax, (3.80, 2.44), "hadrons", size=7, color=ff.BLUE)

    # binned into the fragmentation function
    x0 = 4.55
    for i in range(5):
        h = 0.28 + 0.86 * np.exp(-((i - 1.1) ** 2) / 2.4)
        ct.box(ax, (x0 + i * 0.26, 0.16), 0.21, h, color=ff.BLUE, fc=ff.BLUE,
               lw=0, r=0.03, alpha=0.28 + 0.11 * i, z=3)
    ax.plot([x0 - 0.06, x0 + 5 * 0.26], [0.16, 0.16], color=ff.INK, lw=0.7)
    ct.txt(ax, (6.0, -0.10), r"momentum fraction $z$", size=7, ha="right")
    ct.txt(ax, (x0 + 0.62, 1.76), r"$D^h_q(z)$", size=7.6, color=ff.BLUE)
    ct.arrow(ax, (4.30, 1.20), (4.62, 1.20), color=ff.GREY, lw=0.7, mut=6)

    caption(ax, "(b)  fragmentation instance:\nstruck quark to $D^h_q(z)$")


def featured(ax):
    """The featured route: a colour-singlet q-qbar source injected into the vacuum, evolved in real time for
    the 30-50 a window, its hadrons tallied by species with N_h = U Pi_h U^dagger."""
    ff.blank(ax)
    ax.set_aspect("auto")
    ax.set_xlim(0, 6.0)
    ax.set_ylim(-0.55, 2.85)

    # the injected pair: quark and antiquark joined by one flux link, moving apart along the source axis
    qx, ax_y = (0.55, 1.45), 1.42
    ax.plot([0.62, 1.38], [ax_y, ax_y], color=ff.ORANGE, lw=2.2, alpha=0.75, zorder=3,
            solid_capstyle="round")
    for x, lab in ((0.58, r"$q$"), (1.42, r"$\bar q$")):
        ax.add_patch(Circle((x, ax_y), 0.085, facecolor=ff.BLUE, edgecolor="white", lw=0.5, zorder=5))
        ct.txt(ax, (x, ax_y + 0.30), lab, size=7.4, color=ff.BLUE)
    ct.arrow(ax, (0.50, ax_y - 0.26), (0.18, ax_y - 0.26), color=ff.MIDGREY, lw=0.7, mut=5)
    ct.arrow(ax, (1.50, ax_y - 0.26), (1.82, ax_y - 0.26), color=ff.MIDGREY, lw=0.7, mut=5)
    ct.txt(ax, (1.62, 0.62), "color-singlet source,\n" r"$P\!\sim\!1.5$--$3.1$ GeV", size=7, color=ff.MIDGREY)
    ct.txt(ax, (1.00, 2.20), "injected into\nthe vacuum", size=7, color=ff.ORANGE)

    # real-time evolution across the window
    ct.arrow(ax, (1.95, ax_y), (3.05, ax_y), color=ff.INK, lw=1.0)
    ct.txt(ax, (2.50, ax_y + 0.24), r"evolve $e^{-iHt}$", size=7)
    ct.txt(ax, (2.50, ax_y - 0.26), r"$t = 30$--$50\,a$", size=7, color=ff.MIDGREY)

    # the hadrons that come out: mesons (2 quarks), a baryon pair (3), a kaon pair marked s
    had = [((3.42, 2.02), 0.19, 2, "M"), ((3.46, 1.14), 0.20, 2, "M"), ((3.92, 1.66), 0.22, 3, "B"),
           ((3.90, 0.78), 0.22, 3, r"$\bar B$"), ((3.40, 0.44), 0.17, 2, "K")]
    for n, (xy, r, nq, lab) in enumerate(had):
        ct.hadron(ax, xy, r=r, quarks=nq, color=ff.BLUE, seed=n + 11, z=4, alpha=0.16, lw=0.7)
        ct.txt(ax, (xy[0] + r + 0.17, xy[1]), lab, size=7, color=ff.BLUE)
    ct.txt(ax, (3.62, -0.24), "hadrons", size=7, color=ff.BLUE)

    # the species tally: N_h = U Pi_h U^dagger, three bars
    x0 = 4.85
    ct.arrow(ax, (4.34, 1.20), (4.72, 1.20), color=ff.GREY, lw=0.7, mut=6)
    for i, (lab, h) in enumerate((("M", 1.10), ("B", 0.42), ("K", 0.30))):
        ct.box(ax, (x0 + i * 0.34, 0.16), 0.26, h, color=ff.BLUE, fc=ff.BLUE, lw=0, r=0.03,
               alpha=0.30 + 0.12 * i, z=3)
        ct.txt(ax, (x0 + i * 0.34 + 0.13, -0.08), lab, size=7)
    ax.plot([x0 - 0.06, x0 + 3 * 0.34], [0.16, 0.16], color=ff.INK, lw=0.7)
    ct.txt(ax, (6.1, 2.55), r"count by species, $\mathcal{N}_h$", size=7, ha="right", color=ff.BLUE)
    ct.txt(ax, (6.1, 2.22), r"$B/M$, strangeness", size=7, ha="right", color=ff.MIDGREY)

    caption(ax, "(a)  featured route: evolve a\ncolor-singlet pair, count species")


def main():
    # canvas = printed size (full textwidth), so every label prints at >= 7pt
    ff.setup(truetype=True, tight=False)
    W, H = ff.TEXTWIDTH_IN, 1.63
    fig = plt.figure(figsize=(W, H), dpi=300)
    # panel widths and gaps in inches: (a), (b), (c)
    x, widths, gap, bot = 0.03, (2.33, 2.30, 1.58), 0.13, 0.30
    axes = []
    for w in widths:
        axes.append(fig.add_axes([x / W, bot / H, w / W, (H - bot - 0.01) / H]))
        x += w + gap
    featured(axes[0])
    fragmentation(axes[1])
    lightcone(axes[2])
    ff.save(fig, "ch06_collider")


if __name__ == "__main__":
    main()
