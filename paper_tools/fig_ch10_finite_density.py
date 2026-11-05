#!/usr/bin/env python3
"""Ch. 10 cartoon: the region of the QCD phase diagram the sign problem closes.

Euclidean Monte Carlo is controlled in a wedge around the mu_B = 0 axis -- the
severity of the sign problem is set by mu_B/T, so the usable region opens
upward, not sideways. Outside it sit the conjectured critical point and the
cold dense matter of neutron-star interiors.

Each region carries a sketch of what the matter actually is: separate hadrons
at low T and low mu_B, a deconfined plasma at high T, and closely packed
nucleons at high mu_B. A Hamiltonian simulation has no sign problem by
construction; the cost reappears as circuit depth.
"""

import numpy as np
from matplotlib.patches import Circle, Ellipse, Polygon

import cartoon as ct
import fnalfig as ff
import matplotlib.pyplot as plt


# The axes are not equal-aspect (T and mu_B carry no common unit), so a data
# circle renders as an oval. AR is the x:y stretch that makes a blob look round.
AR = 2.15


def blob(ax, xy, rx, color, alpha, lw=0.6, z=5, ec=None):
    ax.add_patch(Ellipse(xy, 2 * rx, 2 * rx * AR, facecolor=color, alpha=alpha,
                         edgecolor=color if ec is None else ec, lw=lw, zorder=z))


def hadron_gas(ax, xy, n=5, seed=1):
    """Separate color-singlet blobs, well apart: the confined phase."""
    rng = np.random.default_rng(seed)
    x0, y0 = xy
    pts = []
    while len(pts) < n:
        p = (x0 + rng.uniform(-0.085, 0.085), y0 + rng.uniform(-0.085, 0.085) * AR)
        if all(np.hypot(p[0] - q[0], (p[1] - q[1]) / AR) > 0.046 for q in pts):
            pts.append(p)
    for k, p in enumerate(pts):
        blob(ax, p, 0.016, ff.BLUE, 0.20, lw=0.6)
        for m in range(3):
            a = 2 * np.pi * m / 3 + 0.5 * k
            blob(ax, (p[0] + 0.0062 * np.cos(a), p[1] + 0.0062 * AR * np.sin(a)),
                 0.0044, ff.BLUE, 0.95, lw=0.0, z=6)


def plasma(ax, xy, n=44, seed=2):
    """Deconfined quarks and gluons, no blobs left."""
    rng = np.random.default_rng(seed)
    x0, y0 = xy
    for _ in range(n):
        p = (x0 + rng.uniform(-0.15, 0.15), y0 + rng.uniform(-0.075, 0.075))
        ax.plot([p[0]], [p[1]], "o", ms=rng.uniform(1.1, 2.4),
                color=ff.ORANGE if rng.uniform() < 0.6 else ff.YELLOW,
                alpha=0.9, zorder=5)


def packed(ax, xy, seed=3):
    """Nucleons squeezed until they touch: neutron-star interior."""
    x0, y0 = xy
    for i in range(6):
        for j in range(3):
            px = x0 - 0.092 + i * 0.037 + (0.018 if j % 2 else 0)
            py = y0 - 0.037 * AR + j * 0.037 * AR
            blob(ax, (px, py), 0.019, ff.GREEN, 0.80, lw=0.35, z=5, ec="white")


def main():
    ff.setup()
    fig, ax = plt.subplots(figsize=(5.5, 3.15), dpi=300)
    ff.despine(ax, keep=("left", "bottom"))
    ax.set_aspect("auto")

    mu = np.linspace(0, 1.05, 400)
    Tc = 0.86 * np.sqrt(np.clip(1 - (mu / 1.16) ** 2, 0, None))
    cp = 0.46
    icp = int(np.argmin(np.abs(mu - cp)))

    # --- the wedge Euclidean Monte Carlo can actually reach ---------------
    # severity is set by mu_B/T, so the usable region is a wedge opening upward
    slope = 0.34
    ax.add_patch(Polygon([[0, 0], [slope * 1.24, 1.24], [0, 1.24]],
                         closed=True, facecolor=ff.BLUE, alpha=0.13, lw=0,
                         zorder=1))
    ax.plot([0, slope * 1.24], [0, 1.24], color=ff.BLUE, lw=1.0,
            ls=(0, (4, 2)), zorder=3)
    ct.txt(ax, (0.045, 1.10), "Euclidean MC\ncontrolled", size=6, color=ff.BLUE,
           ha="left")
    ct.txt(ax, (0.255, 0.155), r"$\mu_B/T$ grows $\to$ sign problem",
           size=6, color=ff.ORANGE, ha="left")

    # --- transition lines and the conjectured critical point --------------
    ax.plot(mu[:icp], Tc[:icp], color=ff.GREY, lw=1.4, ls=(0, (4, 2)), zorder=4)
    ax.plot(mu[icp:], Tc[icp:], color=ff.INK, lw=1.5, zorder=4)
    ax.plot([mu[icp]], [Tc[icp]], "o", ms=5.0, color=ff.ORANGE, mec="white",
            mew=0.6, zorder=6)
    ct.txt(ax, (0.63, 0.92), "conjectured\ncritical point", size=6,
           color=ff.ORANGE, ha="left")
    ct.arrow(ax, (0.625, 0.895), (0.495, 0.805), color=ff.ORANGE, lw=0.7, mut=6)

    # --- what the matter is, in each region -------------------------------
    plasma(ax, (0.56, 1.06))
    ct.txt(ax, (0.60, 1.20), "quark-gluon plasma", size=6, color=ff.MIDGREY)

    hadron_gas(ax, (0.155, 0.28))
    ct.txt(ax, (0.155, 0.075), "hadron gas", size=6, color=ff.MIDGREY)

    packed(ax, (0.90, 0.26))
    ct.txt(ax, (0.90, 0.055), "dense nuclear matter\n(neutron-star interior)",
           size=6, color=ff.MIDGREY)

    # --- the grid this chapter prices -------------------------------------
    gm = np.linspace(0.34, 0.80, 5)
    gt = np.linspace(0.44, 0.80, 5)
    G1, G2 = np.meshgrid(gm, gt)
    ax.plot(G1.ravel(), G2.ravel(), "+", ms=3.0, mew=0.7, color=ff.GREEN,
            zorder=5)
    ct.txt(ax, (0.60, 0.36), r"$25$-point $(T,\mu_B)$ scan", size=6,
           color=ff.GREEN)

    ax.set_xlim(0, 1.12)
    ax.set_ylim(0, 1.30)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlabel(r"baryon chemical potential $\mu_B$", fontsize=7)
    ax.set_ylabel(r"temperature $T$", fontsize=7)

    fig.subplots_adjust(left=0.055, right=0.985, top=0.97, bottom=0.095)
    ff.save(fig, "ch10_finite_density")


if __name__ == "__main__":
    main()
