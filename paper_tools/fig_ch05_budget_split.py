#!/usr/bin/env python3
"""Ch. 5 per-shot T-budget at the 2033 target, split between the quench thermal-state preparation and the fit
(real-time) evolution, against the RFI targets.

Every bar length is read from the Ch. 5 model (estimates.load(5)); nothing is hand-set. Quench ruling (author,
2026-10-05): each shot is ONE circuit, a coupling ramp (one a_s) plus a thermalization of c / T at the chapter's step,
then the fit evolution to the shot's own fit time (r25 R4). Temperature ruling (author, 2026-10-05): T = 1.5 T_c^lat
with T_c^lat = 150-200 MeV, so the band has a third coordinate. Solid bars at the cheap corner (T = 300 MeV, tau = 1/k,
c = 2), light bars to the expensive corner (T = 225 MeV, tau = 1/(2 pi T), c = 2 pi). Rows: the longest circuit at a_s = 0.2 fm (the box), the second
spacing's late circuit (T3, own step), and the reference circuit evolved to 4 fm/c after the quench (record). The E-rho-OQ
and Gibbs rows of the previous version are retired with those routes.
"""

import numpy as np

import cartoon as ct
import fnalfig as ff
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.transforms import ScaledTranslation

import estimates
from resources import ENVELOPES


def routes():
    """(label, prep_lo, prep_hi, evo_lo, evo_hi, total_lo, total_hi) per row, from the Ch. 5 model."""
    m = estimates.load(5)
    a = m.Assumptions()
    i = m.model(a, "2033").intermediates
    prep_lo, prep_hi = i["quench_prep_t"]                           # 5.93e9, 2.25e10 (152 / 572 steps)
    fit_lo, fit_hi = i["t_fit_deepest"]                             # 7.41e8, 1.26e9 (19 / 32 fit steps)
    tot_lo, tot_hi = i["t_deepest"]                                 # 6.67e9, 2.38e10: the box
    assert np.isclose(tot_lo, prep_lo + fit_lo) and np.isclose(tot_hi, prep_hi + fit_hi)
    p2_lo, p2_hi = i["second_spacing_quench_prep_t"]                # 196 / 755 steps at dt = 0.0075 fm/c
    f2_lo, f2_hi = i["second_spacing_t_fit_deepest"]                # 19 / 48 fit steps
    t2_lo, t2_hi = i["second_spacing_own_step_t_deepest"]           # 8.40e9, 3.17e10
    assert np.isclose(t2_hi, i["t_deepest_campaign"][1])
    ref = i["t_per_shot_2033_mu0_evolution"]                        # 1.570e10: the 4 fm/c reference evolution (record)
    return [
        ("quench + Trotter, longest circuit\n" r"($a_s{=}0.2$ fm)", prep_lo, prep_hi, fit_lo, fit_hi, tot_lo, tot_hi),
        ("same, second spacing\n" r"($a_s{=}0.15$ fm) late circuit", p2_lo, p2_hi, f2_lo, f2_hi, t2_lo, t2_hi),
        (r"quench + Trotter to 4 fm/$c$" "\n(reference)", prep_lo, prep_hi, ref, ref, prep_lo + ref, prep_hi + ref),
    ]


TOT_X = 10.95         # log10 x of the per-shot total column (right of the 3.8e10 reference bar)
INCLUDE = 0.96        # \includegraphics width in app03_qgp_transport.tex
FS = 7.0              # smallest label: canvas is the printed size, so this prints at 7pt


def sci(x):
    e = int(np.floor(np.log10(x)))
    return rf"${x / 10**e:.1f}{{\times}}10^{{{e}}}$"


def sci_range(lo, hi):
    """One value, or 'a--b x 10^e' when both ends share the exponent, else two numbers."""
    if np.isclose(lo, hi, rtol=0.005):
        return sci(hi)
    e_lo, e_hi = int(np.floor(np.log10(lo))), int(np.floor(np.log10(hi)))
    if e_lo == e_hi:
        return rf"${lo / 10**e_lo:.1f}$–${hi / 10**e_hi:.1f}{{\times}}10^{{{e_hi}}}$"
    return sci(lo) + "–" + sci(hi)


def ratio_range(lo, hi):
    return f"{hi:.2g}" if np.isclose(lo, hi, rtol=0.005) else f"{lo:.2g}–{hi:.2g}"


def main():
    ff.setup(truetype=True, tight=False)
    rows = routes()
    ceil28 = ENVELOPES["2028"]["t"][1]
    ceil33 = ENVELOPES["2033"]["t"][1]
    fig, ax = plt.subplots(figsize=(ff.TEXTWIDTH_IN * INCLUDE, 1.79), dpi=300)

    h = 0.30
    for k, (lab, p_lo, p_hi, evo, evo_hi, tot, tot_hi) in enumerate(rows):
        y = len(rows) - 1 - k
        # prep: solid to the fast corner, lighter to the slow corner
        ax.barh(y + 0.17, np.log10(p_lo), height=h, left=0, color=ff.BLUE,
                alpha=0.85, zorder=4, lw=0)
        if p_hi > p_lo:
            ax.barh(y + 0.17, np.log10(p_hi / p_lo), height=h, left=np.log10(p_lo),
                    color=ff.BLUE, alpha=0.35, zorder=4, lw=0)
        ax.barh(y - 0.17, np.log10(evo), height=h, left=0, color=ff.ORANGE,
                alpha=0.85, zorder=4, lw=0)
        if evo_hi > evo:
            ax.barh(y - 0.17, np.log10(evo_hi / evo), height=h, left=np.log10(evo),
                    color=ff.ORANGE, alpha=0.35, zorder=4, lw=0)
        ct.txt(ax, (-0.25, y), lab, size=FS, ha="right", color=ff.INK, linespacing=1.0)
        # totals in one column right of every bar, clear of the target lines
        ax.annotate(sci_range(tot, tot_hi) + "\n" + rf"({ratio_range(tot / ceil33, tot_hi / ceil33)}$\times$)",
                    xy=(TOT_X, y), fontsize=FS,
                    color=ff.DARKGREY, ha="left", va="center", linespacing=1.0)

    for x, lab, col in ((np.log10(ceil28), r"2028 target ($10^{5}$)", ff.GREY),
                        (np.log10(ceil33), r"2033 target ($10^{9}$)", ff.INK)):
        # stop above the inline key row (y = -0.60), which the lines used to cross
        y0, y1 = -0.72, len(rows) - 0.18          # the ylim set below
        ax.axvline(x, ymin=(-0.40 - y0) / (y1 - y0), color=col, lw=0.9,
                   ls=(0, (4, 2)), zorder=3)
        # label to the left of its line, so the 2033 one clears the totals column
        ax.annotate(lab, xy=(x - 0.08, len(rows) - 0.34), fontsize=FS, color=col,
                    ha="right", va="bottom")

    ax.annotate("total per shot\n(" + r"$\times$ 2033 target)", xy=(TOT_X, len(rows) - 0.34),
                fontsize=FS, color=ff.DARKGREY, ha="left", va="bottom", linespacing=1.0)
    ax.set_xlim(0, 13.3)
    ax.set_ylim(-0.72, len(rows) - 0.18)
    ax.set_yticks([])
    ax.set_xticks(range(0, 11, 2))
    ax.set_xticklabels([rf"$10^{{{k}}}$" for k in range(0, 11, 2)])
    ax.set_xlabel("T-gates per shot")
    ff.despine(ax, keep=("bottom",))
    ax.tick_params(left=False, right=False, top=False)
    ax.tick_params(axis="x", direction="out")   # inward ticks cut the inline key

    # inline key, no boxed legend: [swatch] text  [swatch] text along the header row, from the canvas's left
    # edge; the corner definitions sit on the key row under the bars
    LEFT = 0.285
    fig.subplots_adjust(left=LEFT, right=0.985, top=0.86, bottom=0.235)
    ky = len(rows) - 0.34
    head = ax.get_yaxis_transform()             # x: axes fraction, y: data
    inch = fig.dpi_scale_trans + ScaledTranslation(0, ky, head)
    x = 3.0 - LEFT * fig.get_figwidth() * 72     # points, measured from the axes' left edge
    rend = fig.canvas.get_renderer()
    for lab, col in (("quench preparation", ff.BLUE), ("fit evolution", ff.ORANGE)):
        ax.add_patch(Rectangle((x / 72, 1.2 / 72), 5.5 / 72, 5.5 / 72, transform=inch, color=col,
                               alpha=0.85, lw=0, clip_on=False))
        t = ax.annotate(lab, xy=(0, ky), xycoords=head, xytext=(x + 8.0, 0), textcoords="offset points",
                        fontsize=FS, color=col, ha="left", va="bottom", annotation_clip=False)
        x += 8.0 + t.get_window_extent(rend).width * 72 / fig.dpi + 9.0
    ax.annotate(r"(solid: $T{=}300$ MeV, $\tau{=}1/k$, $t_{\rm th}{=}2/T$; light: to $225$ MeV, $1/2\pi T$, $2\pi/T$)",
                xy=(0.05, -0.60), fontsize=FS, color=ff.DARKGREY,
                ha="left", va="center")

    ff.save(fig, "ch05_budget_split")


if __name__ == "__main__":
    for r in routes():
        print(f"{r[0]:48s} prep {r[1]:.3g}-{r[2]:.3g}  evo {r[3]:.4g}-{r[4]:.4g}  total {r[5]:.4g}-{r[6]:.4g}")
    main()
