#!/usr/bin/env python3
"""Ch. 9: T per shot of the 2033 target across the counted variations, against the 10^9 T target.

Every value comes from the Ch. 9 model. The Yukawa band reruns the model at y rho_max = 0.3 and 2
(Assumption yukawa_weight_band).
"""
import dataclasses

import matplotlib.pyplot as plt

import estimates
from estimates import common as c
import fnalfig as ff


def data():
    m = estimates.load(9)
    a = m.Assumptions()
    i = m.model(a, "2033").intermediates
    w_lo, w_hi = a.yukawa_weight_band.lo, a.yukawa_weight_band.hi
    # weight = y rho_max: the band's ends are y = 0.3 at rho_max a = 1 and y = 1 at rho_max a = 2
    y_lo = dataclasses.replace(a, yukawa_coupling_y=c.Assumed(w_lo, "fig"))
    y_hi = dataclasses.replace(a, higgs_rho_max_lat=c.Assumed(w_hi, "fig"))
    band = tuple(m.model(x, "2033").intermediates["t_per_shot"] for x in (y_lo, y_hi))
    return i, band, (w_lo, w_hi), [
        ("Condensate (thermal state alone)", i["t_per_shot_condensate"]),
        (r"$5a$ shot", i["t_per_shot_short"]),
        (r"$10a$ shot ($\kappa=30$)", i["t_per_shot"]),
        (r"$10a$ evolution alone, QSP floor $2\lambda t$", i["t_per_shot_floor"]),
        (r"$10a$ shot, Trotter-step range", i["t_per_shot_range"]),
        (rf"$10a$ shot, $y\rho_{{\max}}={w_lo:g}$–${w_hi:g}$", band),
        (r"$10a$ shot, worst-case commutator bound", i["t_per_shot_bound"]),
        (r"$10a$ shot, $\kappa=10^{2}$", i["t_per_shot_kappa_r22"]),
        (r"$10a$ shot, $\kappa=10^{3}$–$10^{4}$", i["t_per_shot_kappa_physical"]),
    ]


def sci(x):
    e = len(str(int(x))) - 1
    return f"{x / 10 ** e:.1f}", e


def label(v):
    if isinstance(v, tuple):
        (a, ea), (b, eb) = sci(v[0]), sci(v[1])
        if ea == eb:
            return rf"${a}$–${b}\times 10^{{{eb}}}$"
        return rf"${a}\times 10^{{{ea}}}$–${b}\times 10^{{{eb}}}$"
    a, e = sci(v)
    return rf"${a}\times 10^{{{e}}}$"


def main():
    ff.setup()
    i, _, _, rows = data()
    fig, ax = plt.subplots(figsize=(6.5, 2.75))
    n = len(rows)
    ax.axvline(1e9, color=ff.ORANGE, ls="--", lw=0.9)
    ax.axvline(i["t_per_shot"], color=ff.LIGHTGREY, ls=":", lw=0.7)
    for k, (name, v) in enumerate(rows):
        y = n - 1 - k
        col = ff.FNAL_BLUE if name.startswith(r"$10a$ shot ($\kappa=30") else ff.BLUE
        if isinstance(v, tuple):
            ax.plot(v, [y, y], color=col, lw=5, solid_capstyle="butt")
            r = v[1]
        else:
            ax.plot([v], [y], "o", color=col, ms=6.5 if col == ff.FNAL_BLUE else 4.5)
            r = v
        ax.text(r * 1.28, y, label(v), va="center", ha="left", fontsize=7, color=ff.INK)
    ax.text(1e9 * 1.15, 0, "2033 target", color=ff.ORANGE, fontsize=7, ha="left", va="center")
    ax.set(xscale="log", xlim=(3e8, 3e14), ylim=(-0.6, n - 0.4), yticks=range(n),
           yticklabels=[r[0] for r in reversed(rows)], xlabel="T per shot")
    ax.tick_params(axis="y", length=0)
    ax.tick_params(axis="y", which="both", right=False)
    fig.tight_layout()
    return fig


if __name__ == "__main__":
    import sys
    fig = main()
    fig.savefig(sys.argv[1] + ".pdf"); fig.savefig(sys.argv[1] + ".png", dpi=130)
