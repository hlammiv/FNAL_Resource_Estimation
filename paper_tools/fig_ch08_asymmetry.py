#!/usr/bin/env python3
"""Ch. 8: the comparator's reflection asymmetry A_R(k), its thermal weight, and S_int across the wall profiles.

Data: applications/ch08_baryogenesis/sint/sint_results.json (scatter.py). The chapter's comparator is the run with
|y| phi = 0.5, width 2 a (= 1/m_phi), Delta theta_C = 0.1. f_eq is Fermi-Dirac at a T = 0.5 (T_c); 0.9 T_c is a T = 0.45.
"""
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

import estimates
import fnalfig as ff

DATA = Path(__file__).resolve().parent.parent / "applications/ch08_baryogenesis/sint/sint_results.json"


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


def main():
    ff.setup()
    runs = json.loads(DATA.read_text())
    wanted = [r for r in runs if (r["wall"], r["width"], r["dtheta"]) == (0.5, 2.0, 0.1)]
    assert len(wanted) == 1
    run = wanted[0]
    m = estimates.load(8)
    i = m.model(m.Assumptions(), "2033").intermediates
    a = m.Assumptions()
    p0, prior = i["flag_contrast_p0"], float(a.dn_branch_prior.lo)
    amax = run["max_abs_A"]

    # the printed values
    assert f"{s_int(run['rows'], 0.5):.1e}" == "2.1e-02" and f"{s_int(run['rows'], 0.45):.2e}"[:4] in ("2.35", "2.36", "2.34")
    assert round(i["dn_estimate_headroom"], 1) == 1.3
    assert round(amax, 2) == 0.06 and round(p0 * amax, 3) == 0.037 and round(i["dn_estimate"], 3) == 0.056
    band = [r for r in runs if r["dtheta"] == 0.1]
    vals = [s_int(r["rows"], t) for r in band for t in (0.5, 0.45)]
    lo, hi = i["s_int_band"]
    assert lo - 5e-4 <= min(vals) and max(vals) <= hi + 5e-4

    fig = plt.figure(figsize=(6.5, 2.35))
    y0, h = 0.19, 0.70
    axa = fig.add_axes([0.065, y0, 0.235, h])
    cax = fig.add_axes([0.307, y0, 0.012, h])
    axb = fig.add_axes([0.505, y0, 0.175, h])
    cbb = fig.add_axes([0.685, y0, 0.012, h])
    axc = fig.add_axes([0.795, y0, 0.195, h])

    cmap = LinearSegmentedColormap.from_list("div", [ff.BLUE, "#ffffff", ff.ORANGE])
    edges = prior / p0
    x, y, z = grid(run["rows"], "A_R")
    pc = axa.pcolormesh(x, y, z, cmap=cmap, vmin=-edges, vmax=edges, edgecolors="#ffffff", linewidth=0.3)
    axa.set(xlabel=r"$ak_z$", ylabel=r"$ak_y$", yticks=[-np.pi, -np.pi / 2, 0, np.pi / 2], yticklabels=[r"$-\pi$", r"$-\pi/2$", "0", r"$\pi/2$"])
    axa.text(0.97, 0.965, r"(a) $A_R(k)$", transform=axa.transAxes, ha="right", va="top", fontsize=8)
    cb = fig.colorbar(pc, cax=cax, ticks=[-0.05, 0, 0.05])
    cb.ax.set_yticklabels(["-0.05", "0", "0.05"], fontsize=7)
    cb.outline.set_linewidth(0.4); cax.tick_params(length=2, right=True, left=False)
    for v in (edges, -edges):
        cax.axhline(v, color=ff.INK, lw=1.4)
    for v in (amax, -amax):
        cax.axhline(v, color=ff.INK, lw=0.7, ls=":")
    tr = cax.get_yaxis_transform()
    for v, label, yt in ((edges, "depth prior\n" + r"$|p_0A_R|<0.048$", 0.070), (amax, "largest bin\n" + r"$p_0A_R=0.037$," + "\n" + r"$1.3$ times below", 0.026)):
        cax.annotate(label, (1, v), xycoords=tr, xytext=(5.0, yt), textcoords=tr, fontsize=6.3, va="center", ha="left",
                     annotation_clip=False, arrowprops=dict(arrowstyle="-", lw=0.4, color=ff.MIDGREY, shrinkA=1, shrinkB=0))

    x, y, z = grid(run["rows"], "f_eq")
    pb = axb.pcolormesh(x, y, z, cmap="Greys", vmin=0, vmax=0.2, edgecolors="#ffffff", linewidth=0.3)
    axb.set(xlabel=r"$ak_z$", yticks=[-np.pi, -np.pi / 2, 0, np.pi / 2], yticklabels=[])
    axb.text(0.97, 0.965, r"(b) $f_{\rm eq}$ at $T_c$", transform=axb.transAxes, ha="right", va="top", fontsize=8, color=ff.INK)
    cb2 = fig.colorbar(pb, cax=cbb, ticks=[0, 0.1, 0.2])
    cb2.ax.set_yticklabels(["0", "0.1", "0.2"], fontsize=7)
    cb2.outline.set_linewidth(0.4); cbb.tick_params(length=2, right=True, left=False)

    order = [(0.3, 2.0), (0.5, 2.0), (1.0, 2.0), (0.5, 1.0)]
    labels = ["0.3,\n2", "0.5,\n2", "1,\n2", "0.5,\n1"]
    axc.axhspan(lo * 100, hi * 100, color=ff.PALEGREY, alpha=0.5, lw=0)
    for j, key in enumerate(order):
        r = next(q for q in band if (q["wall"], q["width"]) == key)
        axc.plot(j, 100 * s_int(r["rows"], 0.5), "o", color=ff.BLUE, ms=3.5)
        axc.plot(j, 100 * s_int(r["rows"], 0.45), "s", color=ff.ORANGE, ms=3.2)
    axc.plot(1, 100 * s_int(run["rows"], 0.5), "o", mfc="none", mec=ff.INK, ms=7, mew=0.7)
    axc.plot(1, 100 * s_int(run["rows"], 0.45), "o", mfc="none", mec=ff.INK, ms=7, mew=0.7)
    axc.set(xlim=(-0.6, 3.6), ylim=(1.55, 3.3), xticks=range(4), yticks=[2, 2.5, 3], ylabel=r"$S_{\rm int}\;(10^{-2})$")
    axc.set_xticklabels(labels, fontsize=6.5)
    axc.set_xlabel(r"$|y|\phi$, width$/a$", labelpad=1)
    axc.yaxis.labelpad = 1
    axc.tick_params(axis="x", pad=2)
    axc.text(0.03, 0.965, r"(c) $2.1$ at $T_c$", transform=axc.transAxes, va="top", fontsize=6.8, color=ff.BLUE)
    axc.text(0.03, 0.875, r"$2.4$ at $0.9\,T_c$", transform=axc.transAxes, va="top", fontsize=6.8, color=ff.ORANGE)
    axc.text(3.5, 1.62, "band $1.9$\u2013$2.7$", fontsize=6.5, ha="right", va="bottom", color=ff.DARKGREY)
    return fig


if __name__ == "__main__":
    import sys
    fig = main()
    fig.savefig(sys.argv[1] + ".pdf"); fig.savefig(sys.argv[1] + ".png", dpi=200)
