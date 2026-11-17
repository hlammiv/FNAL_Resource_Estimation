#!/usr/bin/env python3
"""Ch. 7: step error of the 2033 proxy against m_phi Dt, with the exact reduced check and the 5^3 scaling.

Each curve runs verlet_proxy for one (lambda, g) set at Phi0 = 0.31 m_phi, m_psi = 0.1 m_phi, against m_phi Dt = 1e-3:
the largest relative error of n_k over shells with n_k >= 0.1 (all fermion shells). The sweep takes a few minutes, so
it is cached in <stem>.json beside the output. The exact marker and the scaling ratio are the Ch. 7 model inputs
exact_step_err_reduced_2033 and exact_over_proxy_step_2033 (from exact_reduced_step, about half an hour per coupling set at T = 1000).
"""
import json
import os
from multiprocessing import Pool
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from estimates import ch07_curved_space as ch7
import fnalfig as ff

SETS = [(1.0, 1.0), (3.0, 0.5), (1.0, 2.0)]
GRID = [1.0, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01]
REF = 1e-3


def _run(args):
    lam, g, mdt = args
    s, f = ch7.verlet_proxy(mdt, lam, 0.31, g)
    return (s.tolist(), f.tolist())


def sweep(cache):
    if cache.exists():
        return json.loads(cache.read_text())
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    jobs = [(lam, g, mdt) for lam, g in SETS for mdt in [REF] + GRID]
    with Pool(8) as p:
        res = dict(zip(jobs, p.map(_run, jobs)))
    out = []
    for lam, g in SETS:
        ref_s, ref_f = map(np.array, res[(lam, g, REF)])
        assert ref_s.max() < 0.1          # every shell above 0.1 is a fermion shell
        keep = ref_f >= 0.1
        err = [float(np.max(np.abs(np.array(res[(lam, g, m)][1]) / ref_f - 1)[keep])) for m in GRID]
        out.append({"lam": lam, "g": g, "err": err, "n_scalar_max": float(ref_s.max())})
    data = {"mdt": GRID, "ref": REF, "sets": out}
    cache.write_text(json.dumps(data, indent=1))
    return data


def main(stem):
    ff.setup(truetype=True, tight=False)   # canvas = printed size (0.72 textwidth)
    d = sweep(Path(stem + ".json"))
    a = ch7.Assumptions()
    x = np.array(d["mdt"])
    errs = np.array([s["err"] for s in d["sets"]])
    # the printed points: 20-46% at m_phi Dt = 1, at most 1.6% at 0.05
    top = errs[:, 0]
    assert (round(100 * top.min()), round(100 * top.max())) == (20, 46), top
    i05 = d["mdt"].index(a.accuracy_omega_dt.value)
    assert round(100 * errs[:, i05].max(), 1) == 100 * a.proxy_step_err_2033.value, errs[:, i05]

    fig, ax = plt.subplots(figsize=(0.72 * ff.TEXTWIDTH_IN, 2.83))
    colors = [ff.BLUE, ff.ORANGE, ff.GREEN]
    for s, c in zip(d["sets"], colors):
        ax.plot(x, s["err"], "-o", color=c, ms=2.8, lw=1.0)
    ax.text(0.0125, 2.0e-4, r"$(3,\,0.5)$", color=ff.ORANGE, fontsize=7, ha="left", va="top")
    ax.text(0.0088, 3.6e-3, r"$(\lambda,g)=(1,\,2)$", color=ff.GREEN, fontsize=7, ha="left", va="bottom")
    ax.text(0.145, 0.058, r"$(1,\,1)$", color=ff.BLUE, fontsize=7, ha="left", va="bottom")
    # slope-2 guide, offset below the data
    g = np.array([0.03, 0.09])
    ax.plot(g, 0.105 * g ** 2, color=ff.GREY, ls=":", lw=0.9)
    ax.text(0.1, 7e-4, "slope 2", color=ff.MIDGREY, fontsize=7, ha="left", va="center")
    # 5% criterion
    ax.axhline(0.05, color=ff.INK, ls="--", lw=0.7)
    ax.text(0.0088, 0.056, r"$5\%$", color=ff.INK, fontsize=7, va="bottom")
    # chosen step and the two checks
    dt = a.accuracy_omega_dt.value
    ax.axvline(dt, color=ff.LIGHTGREY, lw=0.7, zorder=0)
    ax.text(dt * 0.94, 0.8, r"$m_\phi\Delta t=0.05$", fontsize=7, ha="right", va="top", color=ff.DARKGREY)
    proxy = a.proxy_step_err_2033.value
    ax.plot(dt, proxy, "o", mfc="none", mec=ff.INK, ms=7, mew=0.9, zorder=5)
    ax.annotate(r"proxy, at most $1.6\%$", (dt, proxy), xytext=(0.0088, 0.0165), fontsize=7, va="center",
                arrowprops=dict(arrowstyle="-", lw=0.4, color=ff.MIDGREY, shrinkB=3.5))
    scaled = proxy * a.exact_over_proxy_step_2033.value
    ax.plot(dt * 1.0, scaled, "s", mfc="none", mec=ff.PURPLE, ms=4.5, mew=0.9, zorder=6)
    ax.annotate(r"scaled to $5^3$ ($\times 1.2$), about $2\%$", (dt, scaled), xytext=(0.0086, 0.032), fontsize=7, color=ff.PURPLE, va="center",
                arrowprops=dict(arrowstyle="-", lw=0.4, color=ff.PURPLE, shrinkB=2.5))
    exact = a.exact_step_err_reduced_2033.value
    ax.plot(dt, exact, "D", color=ff.PURPLE, ms=4, zorder=6)
    ax.annotate("exact, reduced instance,\n" + r"at most $0.3\%$", (dt, exact), xytext=(0.13, 2.2e-3), fontsize=7, color=ff.PURPLE, va="center",
                arrowprops=dict(arrowstyle="-", lw=0.4, color=ff.PURPLE, shrinkB=2.5))
    ax.set(xscale="log", yscale="log", xlim=(0.008, 1.3), ylim=(1e-4, 1), xlabel=r"$m_\phi\Delta t$",
           ylabel=r"largest relative step error of $\langle n_k\rangle$")
    fig.tight_layout(pad=0.2)
    return fig


if __name__ == "__main__":
    import sys
    fig = main(sys.argv[1])
    fig.savefig(sys.argv[1] + ".pdf"); fig.savefig(sys.argv[1] + ".png", dpi=200)
