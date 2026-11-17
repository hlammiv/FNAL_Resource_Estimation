#!/usr/bin/env python3
"""Ch. 8: the comparator's reflection asymmetry A_R(k), its thermal weight, and S_int across the wall profiles.

Data: applications/ch08_baryogenesis/sint/sint_results.json (scatter.py). The chapter's comparator is the run with
|y| phi = 0.5, width 2 a (= 1/m_phi), Delta theta_C = 0.1. f_eq is Fermi-Dirac at a T = 0.5 (T_c); 0.9 T_c is a T = 0.45.

Usage: fig_ch08_asymmetry.py [OUTDIR]   writes OUTDIR/ch08_asymmetry.pdf (default: Figures/).
The canvas is exactly 6.5 x 2.0 in (no tight bbox), so at \\textwidth every font prints at its nominal size.
"""
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
sys.path[:0] = [str(REPO / "paper_tools"), str(REPO)]

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

import estimates
import fnalfig as ff

DATA = REPO / "applications/ch08_baryogenesis/sint/sint_results.json"
DATA15 = REPO / "applications/ch08_baryogenesis/sint/sint15_results.json"


def s_int(rows, aT):
    f = np.array([1 / (np.exp(r["E"] / aT) + 1) for r in rows])
    return float((f * np.array([r["A_R"] for r in rows])).sum() / f.sum())


def grid(rows, key):
    kz = sorted({r["kz"] for r in rows})
    ky = sorted({((r["ky"] + np.pi) % (2 * np.pi)) - np.pi for r in rows})
    z = np.full((len(ky), len(kz)), np.nan)
    for r in rows:
        y = ((r["ky"] + np.pi) % (2 * np.pi)) - np.pi
        z[min(range(len(ky)), key=lambda i: abs(ky[i] - y)), kz.index(r["kz"])] = r[key]
    dz, dy = kz[1] - kz[0], ky[1] - ky[0]
    return np.r_[kz, kz[-1] + dz] - dz / 2, np.r_[ky, ky[-1] + dy] - dy / 2, z


def main(outdir=None):
    ff.setup()
    plt.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42, "savefig.bbox": None})
    runs = json.loads(DATA.read_text())
    wanted = [r for r in runs if (r["wall"], r["width"], r["dtheta"]) == (0.5, 2.0, 0.1)]
    assert len(wanted) == 1
    run = wanted[0]
    r15 = next(r for r in json.loads(DATA15.read_text()) if (r["wall"], r["width"], r["dtheta"], r["sigma"]) == (0.5, 2.0, 0.1, 1.5))
    m = estimates.load(8)
    i = m.model(m.Assumptions(), "2033").intermediates
    a = m.Assumptions()
    p0, prior = i["flag_contrast_p0"], float(a.dn_branch_prior.lo)
    amax = run["max_abs_A"]

    # the printed values
    assert f"{s_int(run['rows'], 0.5):.1e}" == "2.1e-02" and f"{s_int(run['rows'], 0.45):.2e}"[:4] in ("2.35", "2.36", "2.34")
    assert round(i["dn_estimate_headroom"], 1) == 1.3
    assert f"{r15['S_int']:.1e}" == "1.6e-02" and f"{r15['S_int_09']:.1e}" == "1.7e-02" and i["s_int_taken"] == 1.6e-2   # SR-ch8-1
    assert round(amax, 2) == 0.06 and round(p0 * amax, 3) == 0.037 and round(i["dn_estimate"], 3) == 0.056
    band = [r for r in runs if r["dtheta"] == 0.1]
    vals = [s_int(r["rows"], t) for r in band for t in (0.5, 0.45)]
    lo, hi = i["s_int_band"]
    assert lo - 5e-4 <= min(vals) and max(vals) <= hi + 5e-4

    W, H = 6.5 * 72, 2.0 * 72                        # canvas in points (old PDF printed 145 pt tall); axes placed in points

    def rect(x0, w):
        return [x0 / W, 31 / H, w / W, (H - 46) / H]

    fs = 7.5                                         # in-panel annotation and colorbar size (prints at 7.5 pt)
    fig = plt.figure(figsize=(W / 72, H / 72))
    axa = fig.add_axes(rect(37, 98))
    cax = fig.add_axes(rect(160, 6))
    axb = fig.add_axes(rect(240, 74))
    cbb = fig.add_axes(rect(318, 6))
    axc = fig.add_axes(rect(366, 99))
    for ax, title in ((axa, r"(a) $A_R(k)$"), (axb, r"(b) $f_{\rm eq}$ at $T_c$"), (axc, r"(c) $S_{\rm int}$")):
        ax.set_title(title, loc="left", fontsize=8, pad=3)

    cmap = LinearSegmentedColormap.from_list("div", [ff.BLUE, "#ffffff", ff.ORANGE])
    edges = prior / p0
    x, y, z = grid(run["rows"], "A_R")
    pc = axa.pcolormesh(x, y, z, cmap=cmap, vmin=-edges, vmax=edges, edgecolors="#ffffff", linewidth=0.3)
    axa.set(xlabel=r"$ak_z$", ylabel=r"$ak_y$", yticks=[-np.pi, -np.pi / 2, 0, np.pi / 2], yticklabels=[r"$-\pi$", r"$-\pi/2$", "0", r"$\pi/2$"])
    axa.xaxis.labelpad = 1; axa.yaxis.labelpad = 1
    cb = fig.colorbar(pc, cax=cax, ticks=[-0.05, 0, 0.05])
    cax.yaxis.set_ticks_position("left"); cax.yaxis.set_label_position("left")
    cb.ax.set_yticklabels(["-0.05", "0", "0.05"], fontsize=fs)
    cb.outline.set_linewidth(0.4); cax.tick_params(length=2, left=True, right=False, pad=2)
    for v in (edges, -edges):
        cax.axhline(v, color=ff.INK, lw=1.4)
    for v in (amax, -amax):
        cax.axhline(v, color=ff.INK, lw=0.7, ls=":")
    tr = cax.get_yaxis_transform()
    # the prior label sits level with the colorbar's top edge, the largest-bin label below its dotted line
    for v, label, dy, va, rp in ((edges, "depth prior\n" + r"$|p_0A_R|<0.048$", 3, "center", (0, 0.5)),
                                 (amax, "largest bin\n" + r"$p_0A_R=0.037$," + "\n" + r"$1.3$ times below", -12, "top", (0, 0.9))):
        cax.annotate(label, (1, v), xycoords=tr, xytext=(7, dy), textcoords="offset points", fontsize=fs, va=va, ha="left",
                     linespacing=1.15, annotation_clip=False,
                     arrowprops=dict(arrowstyle="-", lw=0.4, color=ff.MIDGREY, shrinkA=1, shrinkB=0, relpos=rp))

    x, y, z = grid(run["rows"], "f_eq")
    pb = axb.pcolormesh(x, y, z, cmap="Greys", vmin=0, vmax=0.2, edgecolors="#ffffff", linewidth=0.3)
    axb.set(xlabel=r"$ak_z$", yticks=[-np.pi, -np.pi / 2, 0, np.pi / 2], yticklabels=[])
    axb.xaxis.labelpad = 1
    cb2 = fig.colorbar(pb, cax=cbb, ticks=[0, 0.1, 0.2])
    cb2.ax.set_yticklabels(["0", "0.1", "0.2"], fontsize=fs)
    cb2.outline.set_linewidth(0.4); cbb.tick_params(length=2, right=True, left=False, pad=2)

    order = [(0.3, 2.0), (0.5, 2.0), (1.0, 2.0), (0.5, 1.0)]
    labels = ["0.3,\n2", "0.5,\n2", "1,\n2", "0.5,\n1"]
    axc.axhspan(lo * 100, hi * 100, color=ff.PALEGREY, alpha=0.5, lw=0)
    for j, key in enumerate(order):
        r = next(q for q in band if (q["wall"], q["width"]) == key)
        axc.plot(j, 100 * s_int(r["rows"], 0.5), "o", color=ff.BLUE, ms=3.5)
        axc.plot(j, 100 * s_int(r["rows"], 0.45), "s", color=ff.ORANGE, ms=3.2)
    axc.plot(1, 100 * s_int(run["rows"], 0.5), "o", mfc="none", mec=ff.INK, ms=7, mew=0.7)
    axc.plot(1, 100 * s_int(run["rows"], 0.45), "o", mfc="none", mec=ff.INK, ms=7, mew=0.7)
    axc.plot(1, 100 * r15["S_int"], "D", color=ff.INK, ms=4)
    axc.text(1.35, 100 * r15["S_int"], r"$1.6$, 15-site" + "\n" + r"axis at $t_{\max}$", fontsize=fs, va="center", color=ff.INK,
             linespacing=1.1)
    axc.set(xlim=(-0.6, 3.6), ylim=(1.3, 3.5), xticks=range(4), yticks=[1.5, 2, 2.5, 3], ylabel=r"$S_{\rm int}\;(10^{-2})$")
    axc.set_xticklabels(labels, fontsize=fs, linespacing=1.0)
    axc.set_xlabel(r"$|y|\phi$, width$/a$", labelpad=1)
    axc.yaxis.labelpad = 1
    axc.tick_params(axis="x", pad=2)
    # key above the band: marker, then the comparator's value at that temperature
    for yk, mk, ms, c, txt in ((3.32, "o", 3.5, ff.BLUE, r"$2.1$ at $T_c$"), (3.06, "s", 3.2, ff.ORANGE, r"$2.4$ at $0.9\,T_c$")):
        axc.plot(-0.35, yk, mk, color=c, ms=ms, clip_on=False)
        axc.text(-0.12, yk, txt, fontsize=fs, va="center", color=c)
    axc.text(3.5, 2.74, "band $1.9$–$2.7$", fontsize=fs, ha="right", va="bottom", color=ff.DARKGREY)

    outdir = Path(outdir) if outdir else ff.FIGDIR; outdir.mkdir(parents=True, exist_ok=True)
    out = outdir / "ch08_asymmetry.pdf"
    fig.savefig(out)
    print(f"wrote {out}")
    return fig


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else None)
