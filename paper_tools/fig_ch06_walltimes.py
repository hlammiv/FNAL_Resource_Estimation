#!/usr/bin/env python3
"""Ch. 6: serial wall time on one machine per deliverable, against the five-year window.

Every value comes from the Ch. 6 model. An arrow marks a value that is a lower bound: the species counts take the
electric term at its stated lower bound and omit the source insertion and the final Fourier layer, and the
fragmentation count also omits the hadron-momentum transform.
"""
import matplotlib.pyplot as plt

import estimates
import fnalfig as ff


def data():
    m = estimates.load(6)
    i = m.model(m.Assumptions(), "2033").intermediates
    yr = m.SECONDS_PER_YEAR
    d_per_yr = yr / 86400
    first = tuple(x / d_per_yr for x in i["wall_first_2033_days"])
    return i, [
        (r"$B/M$ first result", first, True, r"$13$–$103$ d", ""),
        (r"$B/M$ campaign (trend + $L_\parallel=4$ check)", i["wall_campaign_2033_yr"], True, r"$3.3$–$26$ yr", r"high end $5.2\times$ window"),
        (r"$L_\parallel=4$ check", i["wall_box_check_one_mode_yr"], True, r"$0.46$–$3.7$ yr", ""),
        (r"Static dipole Im$\,V$, three temperatures", tuple(x / yr for x in i["dipole_static"]["wall_scan_s"]), False, r"$0.47$–$0.80$ yr", ""),
        (r"Light-like dipole $\hat q$, three temperatures", tuple(x / yr for x in i["dipole_lightlike"]["wall_scan_s"]), False, r"$3.7$–$6.8$ d", ""),
        (r"Fragmentation function (FCC era)", i["wall_ff_yr"], True, r"$\geq 12$–$1400$ yr", r"$2.3$–$280\times$ window"),
    ]


def main():
    ff.setup()
    i, rows = data()
    h = i["campaign_horizon_yr"]
    n = len(rows)
    fig, ax = plt.subplots(figsize=(6.5, 2.1))
    ax.axvline(h, color=ff.ORANGE, ls="--", lw=0.9)
    for k, (name, (lo, hi), bound, txt, extra) in enumerate(rows):
        y = n - 1 - k
        ax.plot([lo, hi], [y, y], color=ff.BLUE, lw=5, solid_capstyle="butt")
        r = hi
        if bound:
            ax.annotate("", xy=(hi * 2.1, y), xytext=(hi, y), arrowprops=dict(arrowstyle="-|>", color=ff.BLUE, lw=0.9, mutation_scale=7, shrinkA=0, shrinkB=0))
            r = hi * 2.1
        s = txt + (f"; {extra}" if extra else "")
        if lo < h and r * 30 > h and hi < h:                  # the label would cross the window line: put it left of the bar
            ax.text(lo / 1.2, y, s, va="center", ha="right", fontsize=7, color=ff.INK)
        else:
            ax.text(r * 1.2, y, s, va="center", ha="left", fontsize=7, color=ff.INK)
    ax.text(h * 1.15, n - 0.35, "five-year window", color=ff.ORANGE, fontsize=7, ha="left", va="center")
    ax.set(xscale="log", xlim=(1e-3, 3e6), ylim=(-0.6, n), yticks=range(n), yticklabels=[r[0] for r in reversed(rows)],
           xlabel="serial wall time on one machine (yr)")
    ax.tick_params(axis="y", which="both", length=0)
    fig.tight_layout()
    return fig


if __name__ == "__main__":
    import sys
    fig = main()
    fig.savefig(sys.argv[1] + ".pdf"); fig.savefig(sys.argv[1] + ".png", dpi=130)
