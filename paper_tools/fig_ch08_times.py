#!/usr/bin/env python3
"""Ch. 8: 2033 wall times on one machine against the five-year window, from the Ch. 8 model.

The first three rows are lower bounds. The last adds the Gibbs preparation at 3-20 mixing sweeps to each nucleation
shot (bath_mixing_sensitivity); those sweeps are an estimate, not a measurement.
"""
import matplotlib.pyplot as plt

import estimates
import fnalfig as ff


def main():
    ff.setup()
    m = estimates.load(8)
    i = m.model(m.Assumptions(), "2033").intermediates
    yr = i["wall_campaign_s"] / i["campaign_serial_yr"]
    first, plus, camp = i["first_result_yr"], i["scan_plus_first_yr"], i["wall_campaign_s"] / yr
    (lo, hi), (tlo, thi) = i["nucleation_scan_yr_with_prep_band"], i["nucleation_shot_t_with_prep_band"]
    null_days = (i["first_result_null_yr"] - first) * 365.25
    # the printed values
    assert round(i["first_result_days"]) == 38 and round(plus, 1) == 3.3 and round(camp, 1) == 3.8
    assert round(lo) == 78 and float(f"{hi:.1g}") == 500 and f"{tlo:.1e}"[:3] == "4.2" and f"{thi:.1e}"[:3] == "2.7"
    assert round(null_days) == 6 and round(i["scan_plus_first_null_yr"], 1) == 3.3
    assert i["bath_mixing_sweeps"] == (3.0, 20.0)
    window = i["campaign_serial_yr"] / i["campaign_over_horizon"]

    fig, ax = plt.subplots(figsize=(6.5, 1.75))
    rows = [("First result", first), ("Nucleation scan and first result", plus), ("Campaign", camp)]
    n = len(rows) + 1
    for k, (name, v) in enumerate(rows):
        y = n - 1 - k
        ax.annotate("", xy=(v * 1.9, y), xytext=(v, y), arrowprops=dict(arrowstyle="-|>", color=ff.BLUE, lw=1.1, mutation_scale=7, shrinkA=0, shrinkB=0))
        ax.plot(v, y, "o", color=ff.BLUE, ms=4)
    ax.text(first * 2.1, n - 1, r"$\gtrsim 38$ d; null circuit $+6$ d", va="center", fontsize=7)
    ax.text(plus / 1.2, n - 2, r"$\gtrsim 3.3$ yr", va="center", ha="right", fontsize=7)
    ax.text(camp / 1.2, n - 3, r"$\gtrsim 3.8$ yr", va="center", ha="right", fontsize=7)
    ax.barh(0, hi - lo, left=lo, height=0.42, color=ff.ORANGE)
    ax.text(hi * 1.12, 0, "$78$\u2013$500$ yr", va="center", fontsize=7)
    ax.text(window / 1.15, 0, "$0.4$\u2013$2.7\\times 10^{9}$ T per shot", va="center", ha="right", fontsize=7)
    ax.axvline(window, color=ff.INK, ls="--", lw=0.8)
    ax.text(window * 1.1, n - 0.72, "five-year window", fontsize=7, color=ff.INK, va="center")
    ax.set(xscale="log", xlim=(0.03, 4000), ylim=(-0.55, n - 0.4), yticks=range(n),
           yticklabels=["Nucleation scan with Gibbs\npreparation, $3$\u2013$20$ sweeps"] + [r[0] for r in reversed(rows)],
           xlabel="wall time on one machine (yr)")
    ax.tick_params(axis="y", length=0)
    fig.tight_layout()
    return fig


if __name__ == "__main__":
    import sys
    fig = main()
    fig.savefig(sys.argv[1] + ".pdf"); fig.savefig(sys.argv[1] + ".png", dpi=200)
