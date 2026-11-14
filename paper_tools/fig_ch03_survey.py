#!/usr/bin/env python3
"""Ch. 3: per-shot state-preparation depth across the valence spaces, and shots across the truncation ladder.

(a) T per shot against active modes, one band per pair, with the 2033 target. (b) Projection shots for 30% on
|M| against the truncation, from lambda/|M| and the reference overlap. All values come from the Ch. 3 model.
"""
import matplotlib.pyplot as plt

import estimates
import fnalfig as ff


def data():
    m = estimates.load(3)
    a = m.Assumptions()
    i = m.model(a, "2033").intermediates
    rows = [("$^{27}$Al", 24, i["prep_al_sd"]), ("$^{48}$Ca/Ti", 40, i["prep_pf"]), ("$^{76}$Ge/Se", 44, i["prep_jj44"]),
            ("$^{100}$Mo/Ru", 54, i["prep_mo_ru"]), ("$jj55$ pairs", 64, i["prep_jj55"])]
    # shots = (lambda/|M| / (2 x 0.3))^2 / p_ref, the model's projection count; full pf is its lower bound at p_ref = 1
    v = lambda t: t.lo
    mid = (v(a.lam_over_m_f7p3p1) / 0.6) ** 2 / v(a.p_ref_f7p3p1_record)
    return rows, i["simplified_0nubb"], mid


def main(out=None):
    ff.setup()
    rows, sim, mid = data()
    fig, (l, r) = plt.subplots(1, 2, figsize=(6.5, 2.6), gridspec_kw={"width_ratios": [1.5, 1]})
    for name, n, (lo, hi) in rows:
        fits = hi <= 1e9
        l.plot([n, n], [lo, hi], color=ff.GREEN if fits else ff.BLUE, lw=5, solid_capstyle="butt")
        l.annotate(name, (n, hi), xytext=(0, 3), textcoords="offset points", ha="center", fontsize=7, color=ff.INK)
    l.axhline(1e9, color=ff.ORANGE, ls="--", lw=0.9)
    l.text(66.5, 0.75e9, "2033 target", color=ff.ORANGE, fontsize=7, ha="right")
    # 2026-11-13: one rotation per Hamiltonian Pauli string per Trotter step; the string count is not a power of
    # the mode count (jj44 carries fewer strings than pf), so the fourth-power guide is gone
    l.set(yscale="log", xlim=(20, 70), ylim=(5e8, 5e11), xlabel="active modes", ylabel="T per shot, state preparation")
    b, apf = sim["f7p3_projection"], sim["pf_projection_0plus"]
    pts = [(r"$f_{7/2}p_{3/2}$", b["shots"][0]), (r"$+\,p_{1/2}$", mid), (r"full $pf$", apf["shots_lower_bound"][0])]
    x = range(len(pts))
    r.bar(x, [p[1] for p in pts], color=ff.BLUE, width=0.55)
    for k, p in zip(x, pts):
        r.text(k, p[1] * 1.25, f"{p[1]:.1e}".replace("e+0", r"$\times 10^{") + "}$", ha="center", fontsize=7)
    r.set(yscale="log", ylim=(1e4, 3e7), xticks=list(x), xticklabels=[p[0] for p in pts], ylabel=r"projection shots, $30\%$ on $|M|$")
    l.text(0.02, 0.92, "(a)", transform=l.transAxes); r.text(0.04, 0.92, "(b)", transform=r.transAxes)
    fig.tight_layout()
    return fig


if __name__ == "__main__":
    import sys
    fig = main()
    ff.setup()
    fig.savefig(sys.argv[1] + ".pdf"); fig.savefig(sys.argv[1] + ".png", dpi=130)
