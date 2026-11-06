#!/usr/bin/env python3
"""Ch. 10: (a) 2033 campaign years against the five-year window by quark mass and sampling; (b) the model's kappa_2.

Every value comes from the Ch. 10 model. Rows other than the reweighted one are shots times wall time per shot,
serial on one machine.
"""
import matplotlib.pyplot as plt

import estimates
import fnalfig as ff


def sig2(x):
    return f"{float(f'{x:.2g}'):g}"


def sci(x):
    e = len(str(int(x))) - 1
    return rf"${x / 10 ** e:.1f}\times 10^{{{e}}}$"


def main():
    ff.setup()
    m = estimates.load(10)
    i = m.model(m.Assumptions(), "2033").intermediates
    wall = i["wall_per_shot_s"]
    year = m.SECONDS_PER_YEAR
    window = i["campaign_horizon_yr"]
    rows = [(" 313 MeV, reweighted", i["shots_total"]),
            (" 313 MeV, every point", i["shots_total_direct"]),
            (" 400 MeV, every point", i["grid_shots_sc_by_m_q"][400.0]),
            (" 500 MeV, every point", i["grid_shots_sc_by_m_q"][500.0])]

    fig, (ax, bx) = plt.subplots(1, 2, figsize=(6.5, 2.35), width_ratios=[1.6, 1])

    labels = []
    n = len(rows) + 1
    for k, (name, shots) in enumerate(rows):
        y = n - 1 - (k + 1 if k else 0)
        lo, hi = (shots * w / year for w in wall)
        ax.barh(y, hi - lo, left=lo, height=0.5, color=ff.BLUE)
        ax.text(hi * 1.1, y, f"{sig2(lo)}–{sig2(hi)} yr", va="center", ha="left", fontsize=7, color=ff.INK)
        labels.append((y, f"$m_q=${name}\n{sci(shots)} shots"))
    d_lo, d_hi = i["floor_wall_campaign_yr"]
    y = n - 2
    ax.barh(y, d_hi - d_lo, left=d_lo, height=0.3, color=ff.SKY)
    ax.text(d_lo / 1.1, y, "depth limit", va="center", ha="right", fontsize=7, color=ff.INK)
    ax.axvline(window, color=ff.INK, ls="--", lw=0.8)
    ax.text(window * 1.07, n - 0.42, "five-year window", fontsize=7, color=ff.INK, va="center")
    ax.set(xscale="log", xlim=(0.15, 5000), ylim=(-0.5, n - 0.15), xlabel="serial campaign on one machine (yr)")
    ax.set_yticks([y for y, _ in labels], [t for _, t in labels], fontsize=7)
    ax.tick_params(axis="y", length=0)

    grid = i["kappa2_by_point"]
    temps = sorted({T for T, _ in grid})
    mus = sorted({mu for _, mu in grid})
    for mu in mus:
        strong = mu in (0, 600)
        bx.plot(temps, [grid[(T, mu)] for T in temps], color=ff.BLUE if mu == 0 else ff.ORANGE if mu == 600 else ff.LIGHTGREY,
                lw=1.2 if strong else 0.6, marker="o" if strong else None, ms=2.5)
    bx.set(yscale="log", xlim=(92, 228), ylim=(3e-3, 9), xticks=[100, 150, 200], xlabel=r"$T$ (MeV)", ylabel=r"$\kappa_2$")
    bx.text(103, grid[(100, 0)] * 0.92, r"$6.7\times 10^{-3}$", fontsize=7, color=ff.INK, va="center")
    bx.text(203, grid[(200, 0)] * 0.85, f"{grid[(200, 0)]:.2f}", fontsize=7, color=ff.INK, va="center")
    bx.text(203, grid[(200, 600)] * 1.15, f"{grid[(200, 600)]:.1f}", fontsize=7, color=ff.INK, va="center")
    bx.text(155, 0.06, r"$\mu_B=0$", fontsize=7, color=ff.INK, ha="left")
    bx.text(98, 3.3, r"$\mu_B=600$ MeV", fontsize=7, color=ff.INK, ha="left")
    ax.text(0.015, 0.97, "(a)", transform=ax.transAxes, va="top", fontsize=8)
    bx.text(0.97, 0.96, "(b)", transform=bx.transAxes, va="top", ha="right", fontsize=8)
    fig.tight_layout(w_pad=1.0)
    return fig


if __name__ == "__main__":
    import sys
    fig = main()
    fig.savefig(sys.argv[1] + ".pdf"); fig.savefig(sys.argv[1] + ".png", dpi=130)
