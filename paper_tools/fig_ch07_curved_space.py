#!/usr/bin/env python3
"""Ch. 7 cartoon: particle production by an expanding background.

A mode enters the de Sitter era as pure positive frequency and leaves it mixed.
The mixing coefficient beta_k is the whole observable: |beta_k|^2 is the
occupation number of the out-vacuum, and it is nonzero only because the metric
is time dependent. Nothing about this is accessible to a Euclidean calculation,
because the background itself carries the time dependence.
"""

import numpy as np

import cartoon as ct
import fnalfig as ff
import matplotlib.pyplot as plt


def expansion(ax):
    ff.blank(ax)
    ax.set_aspect("auto")
    ax.set_xlim(-1.72, 1.38)
    ax.set_ylim(-0.14, 2.66)

    tt = np.linspace(0, 1, 240)
    y = tt * 2.06
    half = 0.20 + 1.00 * tt ** 1.85

    # Occupation grows through the era and saturates: the same shape the right
    # panel plots against momentum, here against time.
    def occ(f):
        return np.exp(-((np.log(np.clip(f, 1e-3, None) / 0.42)) ** 2) / 0.62)

    # A chaotic, colored background inside the light cone: field fluctuation
    # amplitude growing with the expansion.
    rng = np.random.default_rng(12)
    ny, nx = 190, 150
    gy = np.linspace(0, 1, ny)[:, None]
    gx = np.linspace(-1, 1, nx)[None, :]
    field = np.zeros((ny, nx))
    for k in range(7):                      # a few Fourier modes = mottled field
        kx, ky = rng.uniform(2, 13), rng.uniform(2, 11)
        ph = rng.uniform(0, 2 * np.pi)
        field += np.sin(kx * gx * np.pi + ph) * np.cos(ky * gy * np.pi + ph)
    field *= (0.25 + 1.4 * gy ** 1.6)       # fluctuations grow as it expands
    hw = 0.20 + 1.00 * gy ** 1.85
    field = np.where(np.abs(gx) <= hw / 1.20, field, np.nan)
    ax.imshow(field, extent=(-1.20, 1.20, 0, 2.06), origin="lower",
              aspect="auto", cmap="PuOr", vmin=-3.2, vmax=3.2, alpha=0.42,
              zorder=1, interpolation="bilinear")

    for s_ in (-1, 1):
        ax.plot(s_ * half, y, color=ff.GREEN, lw=1.0, zorder=3)
    ct.txt(ax, (1.31, 1.90), r"$a(t)$", size=7, color=ff.GREEN, ha="right")

    # quanta popping out of the vacuum, denser as the occupation rises
    npart = 0
    for _ in range(420):
        f = rng.uniform(0.05, 1.0)
        if rng.uniform() > occ(f) * 0.95:
            continue
        hwf = 0.20 + 1.00 * f ** 1.85
        px = rng.uniform(-0.88, 0.88) * hwf
        ax.plot([px], [f * 2.06], "o", ms=rng.uniform(1.4, 2.8),
                color=ff.BLUE if rng.uniform() < 0.5 else ff.ORANGE,
                alpha=0.85, zorder=4)
        npart += 1

    # the mode itself: clean in, mixed out
    ty = np.linspace(0.08, 2.02, 700)
    wave = 0.15 * np.sin(12.0 * ty)
    wave = np.where(ty < 1.02, wave, wave + 0.55 * 0.15 * np.sin(18.6 * ty + 0.8))
    ax.plot(wave, ty, color=ff.INK, lw=0.9, zorder=5)

    ax.plot([-1.20, 1.30], [1.02, 1.02], color=ff.GREY, lw=0.6, ls=(0, (3, 2)),
            zorder=3)

    ct.txt(ax, (-1.48, 0.30), "in-vacuum:\npure positive\nfrequency", size=5.8,
           color=ff.MIDGREY, ha="left")
    ct.txt(ax, (0.02, 2.16), "out-state: mixed, populated", size=5.8,
           color=ff.BLUE, va="bottom")
    ct.arrow(ax, (-1.64, 0.06), (-1.64, 2.42), color=ff.INK, lw=0.7, mut=6)
    ct.txt(ax, (-1.52, 2.48), r"$t$", size=6.5)

    ax.annotate("(a)  quanta created by the expansion", xy=(0.5, -0.135),
                xycoords="axes fraction", ha="center", fontsize=6.8,
                color=ff.INK, annotation_clip=False)


def bogoliubov(ax):
    ff.blank(ax)
    ax.set_aspect("auto")
    k = np.linspace(0.06, 1.0, 400)
    n = 0.92 * np.exp(-((np.log(k / 0.24)) ** 2) / 0.95)

    ax.plot([0, 1.04], [0, 0], color=ff.INK, lw=0.7, zorder=3)
    ax.plot([0, 0], [0, 1.06], color=ff.INK, lw=0.7, zorder=3)
    ax.plot(k, n, color=ff.BLUE, lw=1.4, zorder=4)
    ax.fill_between(k, 0, n, color=ff.BLUE, alpha=0.13, lw=0, zorder=2)

    # the resolvable band -- what a finite lattice actually delivers
    ax.axvspan(0.14, 0.62, color=ff.ORANGE, alpha=0.13, lw=0, zorder=1)
    ct.txt(ax, (0.38, 1.00), "resolved\nmomentum bins", size=6, color=ff.ORANGE)

    ax.set_xlim(-0.02, 1.06)
    ax.set_ylim(-0.16, 1.20)
    ct.txt(ax, (0.52, -0.085), r"comoving momentum $k$", size=6.5)
    ct.txt(ax, (0.05, 1.12), r"$\langle n_k\rangle = |\beta_k|^2$", size=7,
           ha="left", color=ff.BLUE)

    ax.annotate(r"(b)  the observable", xy=(0.5, -0.135),
                xycoords="axes fraction", ha="center", fontsize=6.8,
                color=ff.INK, annotation_clip=False)


def main():
    ff.setup()
    fig, axes = plt.subplots(1, 2, figsize=(6.5, 2.25), dpi=300,
                             gridspec_kw=dict(width_ratios=[0.92, 1.08],
                                              wspace=0.14))
    expansion(axes[0])
    bogoliubov(axes[1])
    fig.subplots_adjust(left=0.03, right=0.985, top=0.97, bottom=0.15)
    ff.save(fig, "ch07_curved_space")


if __name__ == "__main__":
    main()
