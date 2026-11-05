#!/usr/bin/env python3
"""Ch. 9 cartoon: what chiral fermions are, and what they cost.

(a) The problem. On a naive lattice the left- and right-handed components of a
    Wilson fermion sit on the same site and the Wilson term couples them, so
    chirality is violated at O(a) -- exactly the channel the weak interaction
    lives in, and exactly the channel an anomaly-mediated observable measures.
(b) One cure, paid for in qubits: domain wall. Give the fermion a fifth
    dimension and the two handednesses bind to opposite walls, so their mixing
    falls exponentially in L5. The register is multiplied by L5.
(c) The other cure, paid for in depth: overlap. No extra dimension, but the
    matrix sign function of the Wilson-Dirac kernel must be synthesized, and
    QSP does that with a polynomial whose degree grows with the condition
    number set by the spectral gap.

In the 100-1000 LQ regime the register binds first, which is why the chapter
takes the overlap route.
"""

import numpy as np

import cartoon as ct
import fnalfig as ff
import matplotlib.pyplot as plt


def wilson(ax):
    ff.blank(ax)
    ax.set_aspect("auto")
    ax.set_xlim(-0.35, 4.05)
    ax.set_ylim(-0.95, 2.55)

    y = 1.10
    ax.plot([0.10, 3.70], [y, y], color=ff.GREY, lw=0.7, zorder=2)
    for i in range(5):
        x = 0.35 + i * 0.75
        # left- and right-handed components sharing one site
        ax.add_patch(plt.Circle((x - 0.085, y), 0.115, facecolor=ff.BLUE,
                                edgecolor="white", lw=0.5, zorder=5))
        ax.add_patch(plt.Circle((x + 0.085, y), 0.115, facecolor=ff.ORANGE,
                                edgecolor="white", lw=0.5, zorder=5))
        # the Wilson term, mixing them on the spot
        ct.wavy(ax, (x - 0.075, y + 0.16), (x + 0.075, y + 0.16), n=2,
                amp=0.05, color=ff.INK, lw=0.7, z=6)

    ct.txt(ax, (0.35, 0.60), r"$\psi_L$", size=6.5, color=ff.BLUE)
    ct.txt(ax, (0.95, 0.60), r"$\psi_R$", size=6.5, color=ff.ORANGE)
    ct.txt(ax, (1.90, 1.86), r"Wilson term mixes them at $O(a)$", size=6.2,
           color=ff.INK)
    ct.arrow(ax, (1.90, 1.72), (1.60, 1.36), color=ff.INK, lw=0.7, mut=6)
    ct.txt(ax, (1.90, 0.10), "chirality broken by the\ndiscretization itself",
           size=6, color=ff.MIDGREY)

    ax.annotate("(a)  a naive lattice fermion", xy=(0.5, -0.185),
                xycoords="axes fraction", ha="center", fontsize=6.8,
                color=ff.INK, annotation_clip=False)


def domainwall(ax):
    ff.blank(ax)
    ax.set_aspect("auto")
    ax.set_xlim(-0.55, 3.95)
    ax.set_ylim(-0.95, 2.55)

    n5 = 7
    for j in range(n5):
        yy = 0.28 + j * 0.27
        ct.lattice(ax, 0.62, yy, 5, 1, 0.52, color=ff.GREY, lw=0.4,
                   site_ms=1.5, alpha=0.5, z=2)
    for yy, col, lab in ((0.28, ff.BLUE, r"$\psi_L$"),
                         (0.28 + (n5 - 1) * 0.27, ff.ORANGE, r"$\psi_R$")):
        ct.lattice(ax, 0.62, yy, 5, 1, 0.52, color=col, lw=0.9, site_ms=3.0,
                   site_color=col, z=4)
        ct.txt(ax, (3.02, yy), lab, size=6.5, color=col, ha="left")

    ct.arrow(ax, (0.24, 0.28), (0.24, 0.28 + (n5 - 1) * 0.27), color=ff.INK,
             lw=0.7, style="<|-|>", mut=6)
    ct.txt(ax, (0.06, 1.09), r"$L_5$", size=7, ha="right")

    ct.txt(ax, (1.72, 2.24), "handedness separated", size=6.2, color=ff.INK)
    ct.txt(ax, (1.72, -0.24), r"mixing $\sim e^{-L_5}$;  register $\times\,L_5$",
           size=6, color=ff.MIDGREY)

    ax.annotate("(b)  domain wall: pay in qubits", xy=(0.5, -0.185),
                xycoords="axes fraction", ha="center", fontsize=6.8,
                color=ff.INK, annotation_clip=False)


def overlap(ax):
    ff.blank(ax)
    ax.set_aspect("auto")
    x = np.linspace(-1, 1, 900)

    ax.plot([-1.05, 1.05], [0, 0], color=ff.INK, lw=0.6, zorder=2)
    ax.plot([0, 0], [-1.30, 1.30], color=ff.INK, lw=0.6, zorder=2)
    ax.plot(x, np.sign(x), color=ff.GREY, lw=1.5, zorder=3, ls=(0, (3, 2)))

    for d, a in ((9, 0.30), (25, 0.60), (81, 1.0)):
        ax.plot(x, np.tanh(d * x * 0.55), color=ff.BLUE, lw=1.0, alpha=a,
                zorder=4)
    ct.txt(ax, (0.50, 0.46), "higher QSP\ndegree", size=5.8, color=ff.BLUE)
    ct.arrow(ax, (0.48, 0.68), (0.18, 0.94), color=ff.BLUE, lw=0.7, mut=6)

    ax.axvspan(-0.10, 0.10, color=ff.ORANGE, alpha=0.16, lw=0, zorder=1)
    ct.txt(ax, (0.0, -1.56), r"gap $\Delta$", size=6, color=ff.ORANGE)

    ax.set_xlim(-1.05, 1.05)
    ax.set_ylim(-1.38, 1.46)
    ct.txt(ax, (0.74, -1.16), r"spec$(H_W)$", size=6.2)
    ct.txt(ax, (-1.0, 1.26), r"$\mathrm{sgn}(H_W)$", size=7, ha="left",
           color=ff.MIDGREY)
    ct.txt(ax, (-1.0, -1.56), r"depth $\sim\kappa\log\varepsilon^{-1}$",
           size=6, color=ff.BLUE, ha="left")

    ax.annotate("(c)  overlap: pay in depth", xy=(0.5, -0.185),
                xycoords="axes fraction", ha="center", fontsize=6.8,
                color=ff.INK, annotation_clip=False)


def main():
    ff.setup()
    fig, axes = plt.subplots(1, 3, figsize=(6.5, 2.05), dpi=300,
                             gridspec_kw=dict(width_ratios=[1.02, 1.0, 0.92],
                                              wspace=0.16))
    wilson(axes[0])
    domainwall(axes[1])
    overlap(axes[2])
    fig.subplots_adjust(left=0.02, right=0.985, top=0.97, bottom=0.17)
    ff.save(fig, "ch09_chiral")


if __name__ == "__main__":
    main()
