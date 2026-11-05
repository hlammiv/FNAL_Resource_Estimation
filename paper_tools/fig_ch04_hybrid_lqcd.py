#!/usr/bin/env python3
"""Ch. 4 cartoon: where the quantum device sits inside a lattice-QCD workflow.

Read left to right. A gauge configuration is generated classically on HPC, as
it always has been. From its links one builds the Dirac matrix M -- sparse, but
of order 10^4 on a side (per flavor), and it is the determinant and the trace of its inverse
that the Monte Carlo needs. Those two numbers, and only those, are handed to a
quantum coprocessor that block-encodes M, applies a QSVT polynomial, and
estimates the trace.

The point of the figure is the size of the last box: the device holds one
linear-algebra primitive, not the gauge field, which is why this chapter's
register is the smallest in the report while its depth is not.
"""

import cartoon as ct
import fnalfig as ff
import matplotlib.pyplot as plt


def main():
    ff.setup()
    fig, ax = plt.subplots(figsize=(6.5, 2.85), dpi=300)
    ff.blank(ax)
    ax.set_xlim(0, 13.0)
    ax.set_ylim(0, 5.4)

    # ---------- (a) the gauge configuration, generated classically ----------
    ct.cube(ax, (0.55, 1.55), 1.85, nx=4, color=ff.GREY, lw=0.45,
            site_ms=1.9, site_color=ff.BLUE, alpha=0.95, skew=(0.34, 0.28))
    ct.txt(ax, (1.85, 0.92), r"gauge links $U_\ell$", size=6.2, color=ff.INK)
    ct.txt(ax, (1.85, 0.60), r"$V=8^3\!\times\!16$, $a\!\sim\!0.12\,$fm",
           size=6, color=ff.MIDGREY)

    ct.txt(ax, (1.85, 4.30), "classical HMC on HPC", size=6.5, color=ff.MIDGREY)

    # ---------- (b) the Dirac matrix built from those links ----------
    ct.arrow(ax, (3.62, 2.72), (4.42, 2.72), color=ff.INK, lw=1.0)
    ct.txt(ax, (4.02, 2.94), "build", size=6, color=ff.MIDGREY)

    ct.matrix_box(ax, (4.60, 1.62), 2.05, 2.20, color=ff.INK, lw=1.0,
                  band=9, fill=ff.BLUE)
    ct.txt(ax, (5.62, 4.02), r"$M$  (Dirac operator)", size=6.8, color=ff.INK)
    ct.txt(ax, (5.62, 1.28), r"sparse, $2.5\!\times\!10^{4}$ a side",
           size=6, color=ff.MIDGREY)
    ct.txt(ax, (5.62, 0.94), r"$\kappa\!\sim\!10^{2}$", size=6,
           color=ff.MIDGREY)

    # what the Monte Carlo actually needs from it
    ct.txt(ax, (7.30, 3.28), r"$\log\det M$", size=7.5, color=ff.ORANGE)
    ct.txt(ax, (7.30, 2.28), r"$\mathrm{Tr}\,M^{-1}$", size=7.5, color=ff.ORANGE)
    ct.arrow(ax, (6.78, 2.98), (7.30, 3.10), color=ff.ORANGE, lw=0.9, mut=6)
    ct.arrow(ax, (6.78, 2.46), (7.30, 2.46), color=ff.ORANGE, lw=0.9, mut=6)
    ct.txt(ax, (7.55, 1.88), "the only two\nnumbers handed off", size=5.8,
           color=ff.MIDGREY)

    # ---------- (c) the quantum coprocessor ----------
    ct.box(ax, (8.42, 0.98), 4.32, 3.62, color=ff.BLUE, fc="#f5fafd", lw=1.1,
           r=0.14, z=2)
    ct.txt(ax, (10.58, 4.28), "quantum coprocessor", size=6.8, color=ff.BLUE)

    ct.arrow(ax, (8.02, 2.72), (8.38, 2.72), color=ff.ORANGE, lw=1.0)

    ct.circuit(ax, (9.02, 1.90), 3.36, nwires=3, color=ff.INK, lw=0.7, z=5,
               wire_color=ff.GREY,
               gates=[(0.11, (0, 2), "BE"), (0.44, (0, 1), "QSVT"),
                      (0.80, 2, "read")])
    ct.txt(ax, (8.88, 1.90), r"$|0\rangle$", size=6, color=ff.MIDGREY, ha="right")
    ct.txt(ax, (8.88, 2.16), r"$|0\rangle$", size=6, color=ff.MIDGREY, ha="right")
    ct.txt(ax, (8.88, 2.42), r"$|\psi\rangle$", size=6, color=ff.MIDGREY,
           ha="right")

    ct.txt(ax, (9.40, 1.34), r"encode $M$", size=5.6, color=ff.MIDGREY)
    ct.txt(ax, (10.53, 1.34), r"apply $1/x$", size=5.6, color=ff.MIDGREY)
    ct.txt(ax, (11.80, 1.34), "estimate trace", size=5.6, color=ff.MIDGREY)

    ct.txt(ax, (10.58, 3.52), r"$83$ logical qubits", size=7,
           color=ff.BLUE)
    ct.txt(ax, (10.58, 3.06),
           "13 site + 4 column + 4 color\n+ 12 QPE + 50 ancilla",
           size=5.6, color=ff.MIDGREY)

    ct.txt(ax, (10.58, 0.60), "one linear-algebra primitive, not the gauge field",
           size=6.2, color=ff.INK)

    fig.subplots_adjust(left=0.005, right=0.995, top=0.99, bottom=0.01)
    ff.save(fig, "ch04_hybrid_lqcd")


if __name__ == "__main__":
    main()
