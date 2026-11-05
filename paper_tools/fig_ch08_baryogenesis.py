#!/usr/bin/env python3
"""Ch. 8 cartoon: the two calculations of a baryogenesis program.

(a) Nucleation. The field sits in a metastable false vacuum; the device measures
    the survival probability, whose exponential decay gives the rate directly,
    with no bounce action and no fluctuation determinant.
(b) Wall scattering. Particles and antiparticles reflect differently off the
    moving wall because the Yukawa coupling is complex, leaving a net charge
    ahead of it. That asymmetry is the seed the sphalerons reprocess.
"""

import numpy as np
from matplotlib.patches import Circle

import cartoon as ct
import fnalfig as ff
import matplotlib.pyplot as plt


def nucleation(ax):
    ff.blank(ax)
    ax.set_aspect("auto")
    ax.set_xlim(-1.15, 1.35)
    ax.set_ylim(-0.55, 1.55)

    # degenerate-well potential with a metastable minimum
    phi = np.linspace(-1.15, 1.25, 900)
    V = 0.72 * (phi ** 2 - 0.72) ** 2 - 0.30 * phi
    V = (V - V.min()) / (V.max() - V.min()) * 1.05 - 0.10
    ax.plot(phi, V, color=ff.INK, lw=1.1, zorder=4)

    # Locate the two minima on the rescaled curve rather than guessing: the
    # markers have to sit on the line, not near it.
    interior = np.where((np.diff(np.sign(np.diff(V))) > 0))[0] + 1
    left = interior[phi[interior] < 0][0]
    right = interior[phi[interior] > 0][-1]
    top = int(np.argmax(V[(phi > phi[left]) & (phi < phi[right])])
              + np.searchsorted(phi, phi[left]))

    ax.plot([phi[left]], [V[left]], "o", ms=5.0, color=ff.ORANGE, mec="white",
            mew=0.6, zorder=6)
    ct.txt(ax, (phi[left], V[left] + 0.30), "false\nvacuum", size=6,
           color=ff.ORANGE)
    ax.plot([phi[right]], [V[right]], "o", ms=4.2, color=ff.GREY, mec="white",
            mew=0.5, zorder=5)
    ct.txt(ax, (phi[right] + 0.06, V[right] + 0.20), "true", size=6,
           color=ff.MIDGREY)

    ct.arrow(ax, (phi[left] + 0.14, V[left] + 0.04),
             (phi[right] - 0.20, V[right] + 0.03), color=ff.ORANGE, lw=0.8,
             ls=(0, (3, 2)), mut=6, z=5)
    ct.txt(ax, (0.5 * (phi[left] + phi[right]), V[top] + 0.14), "tunnel",
           size=5.8, color=ff.ORANGE)

    ct.txt(ax, (0.10, -0.50), r"$\phi$", size=6.5)
    ct.txt(ax, (-1.10, 1.18), r"$V_{\rm eff}(\phi,T)$", size=7, ha="left")

    # inset: survival probability
    axi = ax.inset_axes([0.60, 0.62, 0.38, 0.34])
    t = np.linspace(0, 1, 200)
    axi.plot(t, np.exp(-2.6 * t), color=ff.BLUE, lw=1.1)
    axi.set_xticks([]); axi.set_yticks([])
    for s in axi.spines.values():
        s.set_linewidth(0.5); s.set_color(ff.GREY)
    axi.set_title(r"$P_{\rm FV}(t)$", fontsize=5.8, color=ff.BLUE, pad=1.5)

    ax.annotate("(a)  nucleation", xy=(0.5, -0.145),
                xycoords="axes fraction", ha="center", fontsize=6.8,
                color=ff.INK, annotation_clip=False)


def wall(ax):
    ff.blank(ax)
    ax.set_aspect("equal")
    ax.set_xlim(-2.5, 2.5)
    ax.set_ylim(-2.15, 2.35)

    rng = np.random.default_rng(8)

    # The false-vacuum plasma outside, still symmetric.
    for _ in range(120):
        x, y = rng.uniform(-2.45, 2.45), rng.uniform(-2.05, 2.05)
        if np.hypot(x, y) < 1.30:
            continue
        ax.plot([x], [y], "o", ms=rng.uniform(1.0, 2.2),
                color=ff.BLUE if rng.uniform() < 0.5 else ff.ORANGE,
                alpha=0.55, zorder=2)

    # An expanding bubble of true vacuum, with a wall of finite thickness.
    for rr, a in ((1.30, 0.10), (1.20, 0.14), (1.10, 0.20)):
        ax.add_patch(Circle((0, 0), rr, facecolor=ff.GREEN, alpha=a,
                            edgecolor="none", zorder=1))
    ax.add_patch(Circle((0, 0), 1.16, facecolor="none", edgecolor=ff.GREEN,
                        lw=1.3, zorder=4))
    ct.txt(ax, (0, 0.30), "broken phase", size=6, color=ff.GREEN, z=6)
    ct.txt(ax, (0, -0.02), "(true vacuum)", size=5.6, color=ff.GREEN, z=6)

    # wall velocity
    for a in np.deg2rad([35, 145, 215, 325]):
        p0 = (1.22 * np.cos(a), 1.22 * np.sin(a))
        p1 = (1.62 * np.cos(a), 1.62 * np.sin(a))
        ct.arrow(ax, p0, p1, color=ff.GREEN, lw=0.9, mut=6, z=5)
    ct.txt(ax, (1.42, 1.72), r"wall, $v_w$", size=6, color=ff.GREEN, z=6)

    # A particle reflects weakly off the wall; an antiparticle reflects
    # strongly, because the Yukawa phase is complex. That difference is Delta n.
    ct.arrow(ax, (-2.30, 0.86), (-1.28, 0.62), color=ff.BLUE, lw=1.0, mut=6)
    ct.arrow(ax, (-1.24, 0.50), (-1.78, 0.18), color=ff.BLUE, lw=0.9, mut=6)
    ct.txt(ax, (-2.32, 1.06), r"$p$", size=7, color=ff.BLUE, ha="left")

    ct.arrow(ax, (-2.30, -0.72), (-1.28, -0.52), color=ff.ORANGE, lw=1.0, mut=6)
    ct.arrow(ax, (-1.24, -0.64), (-2.24, -1.14), color=ff.ORANGE, lw=1.5, mut=8)
    ct.txt(ax, (-2.32, -0.52), r"$\bar p$", size=7, color=ff.ORANGE, ha="left")

    # the resulting excess piled up ahead of the wall
    for _ in range(26):
        a = rng.uniform(np.pi * 0.62, np.pi * 1.38)
        r = rng.uniform(1.34, 1.80)
        ax.plot([r * np.cos(a)], [r * np.sin(a)], "o", ms=rng.uniform(1.6, 2.8),
                color=ff.ORANGE, alpha=0.9, zorder=3)
    ct.txt(ax, (-1.52, 1.86), r"excess $\bar p$ ahead of the wall" "\n" r"$\Rightarrow\ \Delta n\neq0$",
           size=6, color=ff.ORANGE, ha="center", z=6)

    ct.txt(ax, (1.58, -1.72), "symmetric phase", size=6, color=ff.MIDGREY)

    ax.annotate("(b)  wall scattering", xy=(0.5, -0.115),
                xycoords="axes fraction", ha="center", fontsize=6.8,
                color=ff.INK, annotation_clip=False)


def main():
    ff.setup()
    fig, axes = plt.subplots(1, 2, figsize=(6.5, 2.35), dpi=300,
                             gridspec_kw=dict(width_ratios=[1.02, 0.98],
                                              wspace=0.13))
    nucleation(axes[0])
    wall(axes[1])
    fig.subplots_adjust(left=0.03, right=0.985, top=0.97, bottom=0.15)
    ff.save(fig, "ch08_baryogenesis")


if __name__ == "__main__":
    main()
