#!/usr/bin/env python3
"""Ch. 5 cartoon: from a heavy-ion collision to the shear viscosity.

(a) What the experiment does: two nuclei collide, a droplet of quark-gluon
    plasma forms and expands anisotropically into the detector. How much of the
    initial spatial asymmetry survives into the final momentum asymmetry is
    what viscosity controls -- a nearly perfect fluid keeps it, a viscous one
    washes it out.
(b) What the device does: an unequal-time stress-tensor autocorrelator in a
    thermal state.
(c) The Kubo relation: eta is the zero-frequency slope of the spectral
    function. That limit is the analytic continuation that makes the Euclidean
    route ill-posed, and is why this observable needs real time.
"""

import numpy as np
from matplotlib.patches import Circle, Ellipse

import cartoon as ct
import fnalfig as ff
import matplotlib.pyplot as plt


def collision(ax):
    ff.blank(ax)
    ax.set_xlim(-2.35, 2.35)
    ax.set_ylim(-2.15, 2.35)

    # detector barrels, seen down the beam
    ct.detector_rings(ax, (0, 0), [1.66, 1.92, 2.16], color=ff.GREY, lw=0.6,
                      alpha=0.45)

    # the almond-shaped overlap: an expanding, anisotropic droplet
    for rr, a in ((1.24, 0.07), (0.96, 0.11), (0.68, 0.17), (0.40, 0.26)):
        ax.add_patch(Ellipse((0, 0), 2 * rr * 0.60, 2 * rr, facecolor=ff.ORANGE,
                             alpha=a, edgecolor="none", zorder=2))
    ct.txt(ax, (0, -0.04), "QGP", size=6.2, color="#7a3300", z=6)

    # spectator nuclei receding along the beam
    for s in (-1, 1):
        ax.add_patch(Circle((s * 1.10, 0), 0.28, facecolor=ff.BLUE, alpha=0.16,
                            edgecolor=ff.BLUE, lw=0.7, ls=(0, (2, 2)), zorder=1))

    # elliptic flow: more and longer tracks along the short axis of the almond
    rng = np.random.default_rng(4)
    for _ in range(90):
        a = rng.uniform(0, 2 * np.pi)
        w = 1.0 + 0.55 * np.cos(2 * a)
        if rng.uniform() > 0.30 + 0.50 * w / 1.55:
            continue
        r0 = 0.58
        L = (1.15 + 0.45 * np.cos(2 * a)) * rng.uniform(0.85, 1.15)
        ax.plot([r0 * np.cos(a), (r0 + L) * np.cos(a)],
                [r0 * np.sin(a), (r0 + L) * np.sin(a)],
                color=ff.ORANGE, lw=0.55, alpha=0.75, zorder=3,
                solid_capstyle="round")

    ct.txt(ax, (0, -1.86), "spatial asymmetry $\\to$ momentum asymmetry",
           size=5.9, color=ff.MIDGREY)
    ax.annotate("(a)  a heavy-ion collision", xy=(0.5, -0.075),
                xycoords="axes fraction", ha="center", fontsize=6.8,
                color=ff.INK, annotation_clip=False)


def correlator(ax):
    ff.blank(ax)
    ax.set_aspect("auto")
    t = np.linspace(0, 1, 500)
    c = np.exp(-3.1 * t) * np.cos(6.0 * t) * 0.86

    ax.plot([0, 1.02], [0, 0], color=ff.INK, lw=0.7, zorder=3)
    ax.plot(t, c, color=ff.ORANGE, lw=1.4, zorder=4)
    ax.fill_between(t, 0, c, color=ff.ORANGE, alpha=0.12, lw=0, zorder=2)

    ax.set_xlim(0, 1.05)
    ax.set_ylim(-0.85, 1.05)
    ct.txt(ax, (0.52, -0.80), r"time separation $t$", size=6.5)
    ct.txt(ax, (0.02, 0.97), r"$\langle T^{xy}(0)\,T^{xy}(t)\rangle_T$", size=7,
           ha="left")
    ct.txt(ax, (0.62, 0.52), "thermal state,\nreal-time evolution", size=6,
           color=ff.MIDGREY, ha="left")
    ax.annotate("(b)  what the device measures", xy=(0.5, -0.155),
                xycoords="axes fraction", ha="center", fontsize=6.8,
                color=ff.INK, annotation_clip=False)


def spectral(ax):
    ff.blank(ax)
    ax.set_aspect("auto")
    w = np.linspace(0, 1, 500)
    rho = 0.92 * np.tanh(2.6 * w) * np.exp(-0.55 * w)

    ax.plot([0, 1.02], [0, 0], color=ff.INK, lw=0.7, zorder=3)
    ax.plot([0, 0], [0, 1.02], color=ff.INK, lw=0.7, zorder=3)
    ax.plot(w, rho, color=ff.BLUE, lw=1.4, zorder=4)

    # the transport coefficient is the zero-frequency slope
    slope = 0.92 * 2.6
    xs = np.linspace(0, 0.40, 20)
    ax.plot(xs, slope * xs, color=ff.ORANGE, lw=1.0, ls=(0, (4, 2)), zorder=5)
    ct.txt(ax, (0.46, 0.82), r"slope $\propto \eta$", size=6.5, color=ff.ORANGE,
           ha="left")
    ct.arrow(ax, (0.44, 0.78), (0.26, 0.62), color=ff.ORANGE, lw=0.7, mut=6)

    ax.plot([0], [0], "o", ms=3.4, color=ff.ORANGE, zorder=6)
    ct.txt(ax, (0.02, -0.115), r"$\omega\to 0$", size=6, color=ff.ORANGE,
           ha="left")

    # What the slope buys: eta/s against the conjectured lower bound.
    ct.txt(ax, (0.99, 0.30), r"$\eta/s$ to 20%,", size=6, color=ff.BLUE,
           ha="right")
    ct.txt(ax, (0.99, 0.18), r"vs. the $1/4\pi$ bound", size=6, color=ff.MIDGREY,
           ha="right")

    ax.set_xlim(-0.03, 1.05)
    ax.set_ylim(-0.20, 1.12)
    ct.txt(ax, (0.52, -0.155), r"frequency $\omega$", size=6.5)
    ct.txt(ax, (0.06, 1.05), r"$\rho(\omega)$", size=7, ha="left")
    ax.annotate(r"(c)  Kubo relation $\to\ \eta$", xy=(0.5, -0.155),
                xycoords="axes fraction", ha="center", fontsize=6.8,
                color=ff.INK, annotation_clip=False)


def main():
    ff.setup()
    fig, axes = plt.subplots(1, 3, figsize=(6.5, 2.25), dpi=300,
                             gridspec_kw=dict(width_ratios=[0.90, 1.05, 1.0],
                                              wspace=0.17))
    collision(axes[0])
    correlator(axes[1])
    spectral(axes[2])
    fig.subplots_adjust(left=0.02, right=0.985, top=0.97, bottom=0.16)
    ff.save(fig, "ch05_qgp_transport")


if __name__ == "__main__":
    main()
