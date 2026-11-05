#!/usr/bin/env python3
"""Ch. 2 cartoon: what a neutrino-nucleus response calculation has to deliver.

Left: the inclusive process -- a DUNE-energy neutrino strikes an argon nucleus,
a charged lepton comes out, and everything hadronic is summed over. Right: why
that is hard on a quantum computer. The response W(omega) is a dense forest of
excited states, and what the device actually returns is not the forest but its
convolution with a finite-width window, whose width sets the cost.
"""

import numpy as np

import cartoon as ct
import fnalfig as ff
import matplotlib.pyplot as plt


def process(ax):
    ff.blank(ax)
    ax.set_xlim(0, 5.0)
    ax.set_ylim(0, 3.2)

    ct.nucleus(ax, (1.55, 1.5), r=0.72, n=13, color=ff.BLUE, alt=ff.SKY, seed=7)
    ct.txt(ax, (1.55, 0.50), r"$^{40}$Ar", size=7, color=ff.BLUE)

    # incoming neutrino
    ct.arrow(ax, (0.05, 2.62), (0.72, 2.26), color=ff.INK, lw=1.0)
    ct.txt(ax, (0.20, 2.80), r"$\nu_\mu$", size=7.5)

    # W exchange into the nucleus
    ct.wavy(ax, (0.72, 2.26), (1.52, 1.72), n=4, amp=0.075,
            color=ff.ORANGE, lw=0.95)
    ct.txt(ax, (1.20, 2.15), r"$W$", size=7, color=ff.ORANGE)

    # outgoing lepton -- the measured leg
    ct.arrow(ax, (0.72, 2.26), (1.92, 3.02), color=ff.INK, lw=1.0)
    ct.txt(ax, (2.06, 3.06), r"$\mu^-$", size=7.5)

    # hadronic final state, summed over
    rng = np.random.default_rng(11)
    for k in range(7):
        ang = np.deg2rad(rng.uniform(-42, 40))
        L = rng.uniform(0.85, 1.45)
        p0 = (1.55 + 0.62 * np.cos(ang), 1.5 + 0.62 * np.sin(ang))
        p1 = (p0[0] + L * np.cos(ang), p0[1] + L * np.sin(ang))
        ct.arrow(ax, p0, p1, color=ff.GREY, lw=0.7, mut=5, z=2)
    ct.txt(ax, (4.02, 1.62), "hadrons,\nsummed over", size=6, color=ff.MIDGREY,
           ha="center")

    ct.txt(ax, (2.5, 0.06), r"inclusive: $\nu_\mu\,^{40}$Ar $\to \mu^- X$",
           size=6.5, color=ff.INK)


def response(ax):
    ff.blank(ax)
    ax.set_aspect("auto")
    w = np.linspace(0, 1, 900)

    # A dense forest of discrete excitations: this is the true finite-volume
    # spectrum, and it is not what the device can resolve.
    rng = np.random.default_rng(4)
    levels = np.sort(rng.uniform(0.06, 0.95, 46))
    strength = np.exp(-((levels - 0.42) ** 2) / (2 * 0.19 ** 2)) * rng.uniform(0.35, 1.0, 46)
    for L, s in zip(levels, strength):
        ax.plot([L, L], [0, s], color=ff.GREY, lw=0.6, zorder=2)

    # What is actually returned: the same strength through a finite window.
    sigma = 0.055
    smeared = np.zeros_like(w)
    for L, s in zip(levels, strength):
        smeared += s * np.exp(-((w - L) ** 2) / (2 * sigma ** 2))
    smeared /= smeared.max() / 0.93
    ax.plot(w, smeared, color=ff.BLUE, lw=1.3, zorder=4)

    # the window itself, drawn at the peak
    wc = 0.42
    g = 0.30 * np.exp(-((w - wc) ** 2) / (2 * sigma ** 2))
    ax.fill_between(w, 0, g, color=ff.ORANGE, alpha=0.30, lw=0, zorder=3)
    ax.plot(w, g, color=ff.ORANGE, lw=0.9, zorder=3)
    ax.annotate("", xy=(wc - sigma, 0.075), xytext=(wc + sigma, 0.075),
                arrowprops=dict(arrowstyle="<|-|>", color=ff.ORANGE, lw=0.9,
                                mutation_scale=6), zorder=6)
    ct.txt(ax, (0.022, 0.70), r"resolution window" "\n" r"$\sigma\!\sim\!50\,$MeV",
           size=6, color=ff.ORANGE, ha="left")
    ct.arrow(ax, (0.175, 0.655), (0.352, 0.145), color=ff.ORANGE, lw=0.7,
             mut=5, z=6)

    ax.set_xlim(0, 1.0)
    ax.set_ylim(0, 1.16)
    ax.plot([0, 1.0], [0, 0], color=ff.INK, lw=0.7)
    ct.txt(ax, (0.5, -0.115), r"energy transfer $\omega$", size=6.5)
    ct.txt(ax, (0.035, 1.05), r"$W^{\mu\nu}(\omega)$", size=7, ha="left")
    ct.txt(ax, (0.985, 0.83), "window-smeared\nresponse (returned)", size=6,
           color=ff.BLUE, ha="right")
    ct.txt(ax, (0.985, 0.36), "true finite-volume\nlevels (unresolved)", size=6,
           color=ff.MIDGREY, ha="right")


def main():
    ff.setup()
    fig, axes = plt.subplots(1, 2, figsize=(6.5, 2.35), dpi=300,
                             gridspec_kw=dict(width_ratios=[1.05, 1], wspace=0.10))
    process(axes[0])
    response(axes[1])
    fig.subplots_adjust(left=0.005, right=0.995, top=0.99, bottom=0.085)
    ff.save(fig, "ch02_nu_nucleus")


if __name__ == "__main__":
    main()
