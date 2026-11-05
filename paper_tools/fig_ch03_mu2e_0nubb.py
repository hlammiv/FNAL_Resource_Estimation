#!/usr/bin/env python3
"""Ch. 3 cartoon: the two lepton-number/flavor-violating matrix elements, and
why they cost differently on a quantum device.

Both are drawn as contact interactions: the BSM physics is integrated out into
a point-like vertex, and everything expensive is the nuclear state the vertex
is sandwiched between.

(a) Coherent mu->e conversion on 27Al. The muon converts on a nucleon and the
    nucleus is left in its ground state, so one register serves for both sides
    of the matrix element.
(b) Neutrinoless double-beta decay. Two neutrons convert to two protons and two
    electrons leave; parent and daughter are different nuclei. Both live in one
    Fock register: the parent is prepared, the operator applied, and the result
    projected onto the daughter ground state -- a second ground-state projection,
    the factor of two in state-preparation cost that separates the Al cell, which
    fits the 2033 envelope, from the 0nubb pairs, which do not (R14 option A).
"""

import numpy as np
from matplotlib.patches import Circle

import cartoon as ct
import fnalfig as ff
import matplotlib.pyplot as plt


def vertex(ax, xy, r=0.20, color=ff.ORANGE, z=6):
    """Fuzzy contact vertex: the BSM operator, with its scale integrated out."""
    for k, (rr, a) in enumerate(((r * 1.75, 0.10), (r * 1.35, 0.18), (r, 0.85))):
        ax.add_patch(Circle(xy, rr, facecolor=color, alpha=a,
                            edgecolor="none", zorder=z + k))


def mu2e(ax):
    ff.blank(ax)
    ax.set_xlim(0, 5.0)
    ax.set_ylim(0, 3.15)

    ct.nucleus(ax, (2.86, 1.24), r=0.76, n=17, color=ff.BLUE, alt=ff.SKY, seed=5)
    ct.txt(ax, (2.60, 0.08), r"$^{27}$Al, same state in and out",
           size=6.2, color=ff.BLUE)

    vx = (1.90, 1.86)
    vertex(ax, vx, r=0.17)

    # incoming muon and outgoing electron, both meeting the same vertex
    ct.arrow(ax, (0.10, 2.72), vx, color=ff.INK, lw=1.05, shrinkB=5)
    ct.txt(ax, (0.26, 2.90), r"$\mu^-$", size=7.5)
    ct.arrow(ax, vx, (4.10, 2.86), color=ff.INK, lw=1.05, shrinkA=5)
    ct.txt(ax, (4.24, 2.94), r"$e^-$", size=7.5)

    # The nucleon leg: the vertex is a four-fermion contact, two leptons
    # against a nucleon, so it needs a hadronic line as well as the two leptons.
    # symmetric V about the vertex-to-nucleus axis: out of the cloud on one
    # arm, back into it on the other, so the nucleus ends where it began
    ct.arrow(ax, (2.39, 1.30), vx, color=ff.BLUE, lw=1.0, shrinkB=6, mut=6)
    ct.arrow(ax, vx, (2.61, 1.64), color=ff.BLUE, lw=1.0, shrinkA=6, mut=6)
    ct.txt(ax, (2.20, 1.26), "N", size=6.5, color=ff.BLUE)

    ct.txt(ax, (1.52, 1.64), r"$\hat{\mathcal{O}}_{\rm BSM}$", size=7.5,
           color=ff.ORANGE)

    ct.txt(ax, (2.6, 3.02), "(a)  Mu2e: one register", size=7, color=ff.INK)


def zeronubb(ax):
    ff.blank(ax)
    ax.set_xlim(0, 5.4)
    ax.set_ylim(0, 3.15)

    ct.nucleus(ax, (1.02, 1.30), r=0.72, n=15, color=ff.BLUE, alt=ff.SKY, seed=2)
    ct.txt(ax, (1.02, 0.16), r"$^{76}$Ge", size=7, color=ff.BLUE)
    ct.nucleus(ax, (4.36, 1.30), r=0.72, n=15, color=ff.GREEN, alt="#7fd4bd",
               seed=9)
    ct.txt(ax, (4.36, 0.16), r"$^{76}$Se", size=7, color=ff.GREEN)

    vx = (2.69, 1.44)
    vertex(ax, vx, r=0.19)
    ct.txt(ax, (2.69, 1.06), r"$\hat{\mathcal{O}}^{0\nu}$", size=7.5,
           color=ff.ORANGE)

    # two neutrons in, two protons out -- the lepton-number-violating core
    for dy in (0.20, -0.20):
        ct.arrow(ax, (1.80, 1.44 + dy), vx, color=ff.BLUE, lw=0.9, mut=6,
                 shrinkB=6)
        ct.arrow(ax, vx, (3.58, 1.44 + dy), color=ff.GREEN, lw=0.9, mut=6,
                 shrinkA=6)
    ct.txt(ax, (1.95, 1.78), r"$2n$", size=6.2, color=ff.BLUE)
    ct.txt(ax, (3.45, 1.78), r"$2p$", size=6.2, color=ff.GREEN)

    # the two emitted electrons
    ct.arrow(ax, vx, (1.98, 2.66), color=ff.INK, lw=1.0, mut=6, shrinkA=6)
    ct.arrow(ax, vx, (3.42, 2.66), color=ff.INK, lw=1.0, mut=6, shrinkA=6)
    ct.txt(ax, (1.82, 2.76), r"$e^-$", size=7.5)
    ct.txt(ax, (3.58, 2.76), r"$e^-$", size=7.5)

    ct.txt(ax, (2.7, 0.06), "parent prepared, daughter projected",
           size=6.2, color=ff.MIDGREY)
    ct.txt(ax, (2.7, 3.00), r"(b)  $0\nu\beta\beta$: one register, two projections", size=7,
           color=ff.INK)


def main():
    ff.setup()
    fig, axes = plt.subplots(1, 2, figsize=(6.5, 2.35), dpi=300,
                             gridspec_kw=dict(width_ratios=[0.96, 1.04],
                                              wspace=0.05))
    mu2e(axes[0])
    zeronubb(axes[1])
    fig.subplots_adjust(left=0.005, right=0.995, top=0.99, bottom=0.01)
    ff.save(fig, "ch03_mu2e_0nubb")


if __name__ == "__main__":
    main()
