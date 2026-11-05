#!/usr/bin/env python3
"""Ch. 11 cartoon: why a qudit is the natural register for a gauge link.

Left: an SRF cavity's Fock ladder supplies d = 32 levels directly, so one cavity
carries one group element of G. The qubit alternative spends ceil(log2 |G|)
qubits on the same link and then has to synthesise group multiplication and the
non-Abelian Fourier transform out of Clifford+T, and leaks into states that are
not group elements at all.

Right: the stack this chapter has to deliver above that hardware -- it is an
enabling-technology chapter, so the figure is an architecture, not a physical
process.
"""

import numpy as np

import cartoon as ct
import fnalfig as ff
import matplotlib.pyplot as plt


def encoding(ax):
    ff.blank(ax)
    ax.set_aspect("auto")
    ax.set_xlim(0, 5.0)
    ax.set_ylim(0, 3.5)

    # SRF cavity with its Fock ladder
    ct.box(ax, (0.22, 0.55), 1.20, 2.35, color=ff.BLUE, fc="#f2f8fd", lw=1.0,
           r=0.30)
    for j in range(7):
        y = 0.85 + j * 0.31
        lw = 0.9 if j < 6 else 0.6
        ax.plot([0.42, 1.22], [y, y], color=ff.BLUE, lw=lw,
                alpha=1.0 if j < 6 else 0.45, zorder=6)
    ct.txt(ax, (0.82, 3.06), "SRF cavity", size=6.4, color=ff.BLUE)
    ct.txt(ax, (0.82, 0.34), r"$d\!=\!32$ Fock levels", size=6, color=ff.MIDGREY)
    ct.txt(ax, (1.52, 1.78), r"$\vdots$", size=7, color=ff.BLUE, ha="left")

    ct.arrow(ax, (1.78, 1.72), (2.52, 1.72), color=ff.INK, lw=0.9)
    ct.txt(ax, (2.15, 1.92), "one link", size=5.8)

    # the gauge link it encodes
    ax.plot([2.85, 4.45], [1.72, 1.72], color=ff.ORANGE, lw=1.6, zorder=4)
    for xx in (2.85, 4.45):
        ax.plot([xx], [1.72], "o", ms=4.6, color=ff.INK, mec="white", mew=0.5,
                zorder=5)
    ct.txt(ax, (3.65, 2.02), r"$U_\ell \in G$", size=7.5, color=ff.ORANGE)
    ct.txt(ax, (3.65, 1.36), r"$|G|\leq 32$", size=6, color=ff.MIDGREY)

    ct.txt(ax, (2.55, 0.72), "qubit route: "
           r"$\lceil\log_2|G|\rceil$ qubits,", size=6, color=ff.MIDGREY,
           ha="left")
    ct.txt(ax, (2.55, 0.44), r"synthesized $U_{\rm mult}$, $U_F$, and leakage",
           size=6, color=ff.MIDGREY, ha="left")

    ax.annotate("(a)  one cavity, one group element", xy=(0.5, -0.045),
                xycoords="axes fraction", ha="center", fontsize=6.8,
                color=ff.INK, annotation_clip=False)


def stack(ax):
    ff.blank(ax)
    ax.set_aspect("auto")
    ax.set_xlim(0, 4.0)
    ax.set_ylim(0, 3.5)

    layers = [
        ("physics specification", ff.GREEN, 2.74),
        ("middleware: gate error " r"$\to$" " observable error", ff.GREEN, 2.20),
        (r"compiler: $U_{\rm mult}$, $U_{\rm inv}$, $U_F$", ff.ORANGE, 1.66),
        ("calibration: ML pulse synthesis", ff.ORANGE, 1.12),
        (r"SRF qudit hardware, $\epsilon_q\!\leq\!10^{-4}$", ff.BLUE, 0.58),
    ]
    for lab, col, y in layers:
        ct.box(ax, (0.18, y), 3.64, 0.42, color=col, fc="white", lw=0.9, r=0.07)
        ct.txt(ax, (2.0, y + 0.21), lab, size=6.1, color=col)

    # arrows sit fully in the gaps between stack boxes (TR-04 review fix)
    for y in (2.73, 2.19, 1.65, 1.11):
        ct.arrow(ax, (2.0, y), (2.0, y - 0.10), color=ff.GREY, lw=0.7, mut=5)

    ct.txt(ax, (2.0, 3.30), "the stack this chapter must deliver", size=6.4,
           color=ff.INK)

    ax.annotate("(b)  enabling-technology architecture", xy=(0.5, -0.045),
                xycoords="axes fraction", ha="center", fontsize=6.8,
                color=ff.INK, annotation_clip=False)


def main():
    ff.setup()
    fig, axes = plt.subplots(1, 2, figsize=(6.5, 2.45), dpi=300,
                             gridspec_kw=dict(width_ratios=[1.12, 0.88],
                                              wspace=0.10))
    encoding(axes[0])
    stack(axes[1])
    fig.subplots_adjust(left=0.01, right=0.99, top=0.98, bottom=0.075)
    ff.save(fig, "ch11_qudit")


if __name__ == "__main__":
    main()
