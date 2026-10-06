#!/usr/bin/env python3
"""LP-minimal odd QET polynomial for the inverse (Ch. 4, round R11, ruling ch04-dinv-500-1e4 option B').

Odd-parity copy of apply_log/ch04_roundB_lp_degree.py. Not part of the `estimates` package: it needs
numpy and scipy. It is the provenance of `Assumptions.d_inv_2033 = 500` and `d_inv_codesign = 5e3`.

Problem. Over odd polynomials p of degree d, written in the Chebyshev basis T_1, T_3, ..., T_d,

    minimize   max over x in [1/kappa, 1] of | p(x) - f(x) |,   f(x) = 1/(2 kappa x),
    subject to | p(x) | <= 1 on [-1, 1]      (odd, so imposed on [0, 1])

f is 1/x normalized so that its peak, at x = 1/kappa, is 1/2 (the room the QET bound leaves).
The relative error eps_rel is the uniform error divided by that peak, err / 0.5.

Finding (R11 triage, reproduced here): err(d) = 0.5 exp(-d/kappa), i.e.

    d_inv(LP) = kappa ln(1/eps_rel),    eps_rel measured against the peak of 1/x.

No kappa appears inside the log. At kappa = 1e2 and eps_rel = e^-5 = 0.674% the degree is 500;
at kappa = 1e3 the same eps_rel gives 5e3 (extrapolated from the law: the LP runs out of memory
near d ~ 5e3 on the default grids).

    python3 ch04_r11_lp_odd.py 30 121 [more degrees ...]
    python3 ch04_r11_lp_odd.py 100 401 501 601
"""
import math
import sys

import numpy as np
from scipy.optimize import linprog


def grids(kappa, n_in, n_all):
    lam = 1.0 / kappa
    k = np.arange(n_in)
    x_in = 0.5 * (1 + lam) - 0.5 * (1 - lam) * np.cos(np.pi * k / (n_in - 1))   # [1/kappa, 1]
    k = np.arange(n_all)
    x_cheb = np.cos(np.pi * k / (2 * (n_all - 1)))                               # [0, 1], dense at 1
    x_lin = np.linspace(0, 1, n_all)
    x_gap = np.linspace(0, lam, n_all // 2)                                      # below 1/kappa
    return x_in, np.unique(np.concatenate([x_cheb, x_lin, x_gap, x_in]))


def cheb_odd(x, d):
    theta = np.arccos(np.clip(x, -1, 1))
    return np.cos(np.outer(theta, np.arange(1, d + 1, 2)))


def min_err(kappa, d, n_in=None, n_all=None):
    """Smallest achievable uniform error at odd degree d."""
    n_in = n_in or max(2000, 8 * d)
    n_all = n_all or max(3000, 8 * d)
    x_in, x_all = grids(kappa, n_in, n_all)
    f = 1.0 / (2 * kappa * x_in)
    a_in, a_all = cheb_odd(x_in, d), cheb_odd(x_all, d)
    n = a_in.shape[1]
    cost = np.zeros(n + 1)
    cost[-1] = 1
    ones = np.ones((len(x_in), 1))
    zeros = np.zeros((len(x_all), 1))
    a_ub = np.block([[a_in, -ones], [-a_in, -ones], [a_all, zeros], [-a_all, zeros]])
    b_ub = np.concatenate([f, -f, np.ones(len(x_all)), np.ones(len(x_all))])
    res = linprog(cost, A_ub=a_ub, b_ub=b_ub, bounds=[(None, None)] * n + [(0, None)], method="highs")
    if res.status != 0:
        return float("inf")
    return res.x[-1]


if __name__ == "__main__":
    kappa = float(sys.argv[1])
    for d in (int(s) for s in sys.argv[2:]):
        e = min_err(kappa, d)
        law = 0.5 * math.exp(-d / kappa)
        print(f"kappa={kappa:g} d={d}: err={e:.4g}  eps_rel={e / 0.5:.4g}  law 0.5 e^(-d/kappa)={law:.4g}  "
              f"kappa ln(1/eps_rel)={kappa * math.log(0.5 / e):.1f}")
