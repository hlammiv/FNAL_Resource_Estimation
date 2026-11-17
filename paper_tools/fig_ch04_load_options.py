#!/usr/bin/env python3
"""Ch. 4: 2033 campaign time for each way of loading the link-variable oracle, against the five-year window.

Campaign walls come from the Ch. 4 model. The multiplexed form has F* = 1, ten times the ten-factory baseline;
the unary-iteration form F* = 4; three parallel subtrees reach the baseline rate.
"""
import matplotlib.pyplot as plt

import estimates
import fnalfig as ff


def main():
    ff.setup(truetype=True)
    m = estimates.load(4)
    i = m.model(m.Assumptions(), "2033").intermediates
    base = i["wall_campaign_baseline_yr"]
    rows = [("Multiplexed rotation ($F^*=1$)", tuple(10 * x for x in base)),
            ("Unary-iteration load ($F^*=4$)", i["wall_campaign_unary_iteration_yr"]),
            (r"Three parallel subtrees ($F^*\approx 12$)", base)]
    fig, ax = plt.subplots(figsize=(6.5, 1.5))
    for y, (name, (lo, hi)) in enumerate(reversed(rows)):
        ax.barh(y, hi - lo, left=lo, height=0.5, color=ff.BLUE if hi <= i["campaign_horizon_yr"] else ff.ORANGE)
        left = lo < i["campaign_horizon_yr"] < hi * 1.8
        ax.text(lo / 1.08 if left else hi * 1.08, y, f"{lo:.2g}--{hi:.2g} yr".replace("--", "–"), va="center", ha="right" if left else "left", fontsize=7)
    ax.axvline(i["campaign_horizon_yr"], color=ff.INK, ls="--", lw=0.8)
    ax.text(i["campaign_horizon_yr"] * 1.04, 2.45, "five-year window", fontsize=7, color=ff.INK)
    ax.set(xscale="log", xlim=(0.3, 60), ylim=(-0.5, 2.8), yticks=range(3), yticklabels=[r[0] for r in reversed(rows)], xlabel="campaign on one machine (yr)")
    fig.tight_layout()
    return fig


if __name__ == "__main__":
    import sys
    fig = main()
    fig.savefig(sys.argv[1] + ".pdf"); fig.savefig(sys.argv[1] + ".png", dpi=130)
