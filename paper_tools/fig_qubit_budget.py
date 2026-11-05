#!/usr/bin/env python3
"""Rebuild of the overview qubit-budget figure, in house style.

Contours of constant logical-qubit count in the (sites per dimension, qubits
per site) plane, for 2D and 3D lattices. Where a theory's qubits-per-site line
crosses a contour is the largest lattice that fits that budget. Reproduces the
content of the previous Figures/lattice_qubits.png, which had no source script.

n_lq = q_site * N^D, so at fixed n_lq the contour is q = n_lq / N^D.
"""

import numpy as np

import cartoon as ct
import fnalfig as ff
import matplotlib.pyplot as plt

# qubits per site for each digitization, as used through the report
THEORIES = [
    (r"$\mathbb{Z}_2$", 1.0, ff.MIDGREY),
    (r"$\mathbb{Z}_3$", 2.0, ff.MIDGREY),
    ("SU(2)", 12.0, ff.MIDGREY),
    ("SU(3)", 27.0, ff.MIDGREY),
]


def main():
    ff.setup()
    fig, ax = plt.subplots(figsize=(5.4, 3.5), dpi=300)

    N = np.linspace(2, 10, 400)

    # budget contours
    for nlq, alpha in ((1e2, 1.0), (1e3, 0.45)):
        ax.plot(N, nlq / N ** 2, color=ff.BLUE, lw=1.5, alpha=alpha, zorder=4)
        ax.plot(N, nlq / N ** 3, color=ff.ORANGE, lw=1.5, alpha=alpha,
                ls=(0, (4, 2)), zorder=4)

    # the region a 2028-class machine opens, in 2D
    ax.fill_between(N, 1.0, np.minimum(1e2 / N ** 2, 40), color=ff.BLUE,
                    alpha=0.07, lw=0, zorder=1)

    # qubits-per-site lines, and the largest lattice each admits
    for lab, q, col in THEORIES:
        ax.plot([2, 10], [q, q], color=col, lw=0.6, ls=":", zorder=3)
        ct.txt(ax, (10.15, q), lab, size=6.2, color=ff.INK, ha="left")
        for nlq, D, mcol in ((1e2, 2, ff.BLUE), (1e2, 3, ff.ORANGE)):
            n = (nlq / q) ** (1 / D)
            if 2 <= n <= 10:
                ax.plot([n], [q], "o", ms=4.0, color=mcol, mec="white", mew=0.5,
                        zorder=6)

    ax.set_xscale("linear")
    ax.set_yscale("log")
    ax.set_xlim(2, 10)
    ax.set_ylim(1, 40)
    ax.set_xlabel("sites per dimension $N$")
    ax.set_ylabel("qubits per lattice site")
    ax.set_yticks([1, 2, 5, 10, 20, 40])
    ax.set_yticklabels(["1", "2", "5", "10", "20", "40"])
    ax.minorticks_off()

    # Curve labels tangent to the contours, house style: no boxed legend.
    fig.tight_layout(pad=0.35)
    fig.canvas.draw()
    ff.label_along(ax, 4.2, lambda n: 1e2 / n ** 2, r"$n_{lq}=10^{2}$, 2D",
                   ff.BLUE, side=1, pad=6)
    ff.label_along(ax, 8.2, lambda n: 1e3 / n ** 2, r"$10^{3}$, 2D",
                   ff.BLUE, side=1, pad=6)
    ff.label_along(ax, 3.1, lambda n: 1e2 / n ** 3, r"$10^{2}$, 3D",
                   ff.ORANGE, side=-1, pad=6)
    ff.label_along(ax, 6.3, lambda n: 1e3 / n ** 3, r"$10^{3}$, 3D",
                   ff.ORANGE, side=-1, pad=6)

    ct.txt(ax, (5.6, 1.28), "opens by 2028 (2D)", size=6.2, color=ff.BLUE)

    # what another decade of logical qubits buys
    ax.annotate("", xy=(9.13, 12.0), xytext=(3.05, 12.0),
                arrowprops=dict(arrowstyle="-|>", color=ff.DARKGREY, lw=0.9,
                                mutation_scale=7), zorder=7)
    ct.txt(ax, (5.05, 13.9), r"$\times 10$ logical qubits", size=6,
           color=ff.DARKGREY)

    ff.save(fig, "qubit_budget")


if __name__ == "__main__":
    main()
