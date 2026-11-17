#!/usr/bin/env python3
"""Ch. 8: the nucleation measurement, from the Ch. 8 model (markup board C8R-2, 2026-11-09).

-ln P_FV(t) in the exponential regime is the decrement Gamma V t, linear in t. Each point measures it at two times
with binomial errors sqrt(delta / N). The model does not compute the energy spread Delta E, so the Zeno region
before the exponential regime is shaded at the chapter's stated crossover (1/m_phi), not computed.

Usage: fig_ch08_nucleation.py [OUTDIR]   writes OUTDIR/ch08_nucleation.pdf (default: Figures/).
The canvas is exactly 6.5 in wide (no tight bbox), so at \\textwidth every font prints at its nominal size.
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
sys.path[:0] = [str(REPO / "paper_tools"), str(REPO)]

import numpy as np
import matplotlib.pyplot as plt

import estimates
import fnalfig as ff


def main(outdir=None):
    ff.setup()
    plt.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42, "savefig.bbox": None})
    m = estimates.load(8)
    a = m.Assumptions()
    i28 = m.model(a, "2028").intermediates
    i33 = m.model(a, "2033").intermediates

    # 2028: vacuum decay on 6^2, two times at 10 and 20 of 20 steps (1 and 2 / m_phi)
    t28 = i28["t_window_over_mphi"]
    d28 = i28["slope_decrement_late"]
    r28 = i28["slope_time_ratio"]
    n28 = i28["slope_shots_per_time"]
    # 2033: thermal decay on the 15 x 8 cylinder, two times at 33 and 100 of 100 steps, at T_c and 0.9 T_c
    t33 = 10.0
    r33 = i33["nucleation_early_steps"] / 100
    d33 = i33["nucleation_decrement_late"]
    n33 = i33["nucleation_shots_per_time"]
    d09 = d33 * i33["decrement_ratio_0p9tc"]
    n09 = i33["nucleation_shots_per_time_0p9tc"]
    # the values the chapter prints
    assert t28 == 2.0 and r28 == 0.5 and round(n28, -3) == 18000 and abs(d28 - 8.33e-3) < 1e-5
    assert round(r33, 2) == 0.33 and abs(d33 - 2.5e-3) < 1e-12 and round(n33, -3) == 30000
    assert round(i33["decrement_ratio_0p9tc"], 2) == 0.5 and round(n09, -3) == 59000

    fig = plt.figure(figsize=ff.WIDE)
    ax = fig.add_axes([0.085, 0.155, 0.905, 0.825])
    ax.axvspan(0.4, 1.0, color=ff.PALEGREY, alpha=0.5, lw=0)
    ax.text(0.63, 1.3e-2, "Zeno regime\n(not computed)", ha="center", va="center", fontsize=7.5, color=ff.DARKGREY)
    # label: text, x where it sits along the curve, side (+1 above, -1 below)
    rows = [(r"2028: vacuum, $6^2$", t28, d28, r28, n28, ff.BLUE, 3.6, -1),
            (r"2033: $T_c$, $15\times 8$", t33, d33, r33, n33, ff.ORANGE, 1.75, 1),
            (r"2033: $0.9\,T_c$", t33, d09, r33, n09, ff.GREEN, 5.5, -1)]
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(0.4, 12); ax.set_ylim(1e-4, 3e-2)
    for name, tl, dl, r, n, c, xl, side in rows:
        g = dl / tl                                   # decrement per unit 1/m_phi time
        tt = np.geomspace(1.0, 12.0, 50)
        ax.plot(tt, g * tt, color=c, lw=0.9)
        ts = np.array([r * tl, tl]); ds = g * ts
        ax.errorbar(ts, ds, yerr=np.sqrt(ds / n), fmt="o", color=c, ms=3.5, capsize=2, lw=0.8)
        ff.label_along(ax, xl, lambda x, g=g: g * x, name, c, side=side, pad=8.0, fontsize=8)
    ax.set_xlabel(r"$m_\phi t$", labelpad=1); ax.set_ylabel(r"$-\ln P_{\rm FV}(t)$", labelpad=2)
    ff.despine(ax)
    outdir = Path(outdir) if outdir else ff.FIGDIR; outdir.mkdir(parents=True, exist_ok=True)
    out = outdir / "ch08_nucleation.pdf"
    fig.savefig(out)
    print(f"wrote {out}")
    return fig


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else None)
