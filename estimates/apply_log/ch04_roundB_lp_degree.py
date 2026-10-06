#!/usr/bin/env python3
"""LP-minimal QET degree for the normalized logarithm (Ch. 4, rulings R7 and 'SHIFT + 1%', 2026-09-28).

Not part of the `estimates` package: it needs numpy and scipy, and the package is
standard-library only. It is the provenance of `Assumptions.d_log_2028 = 84` and of the
other LP degrees carried in estimates/ch04_hybrid_lqcd.py.

Problem. Over even polynomials p of degree d, written in the Chebyshev basis T_0, T_2, ..., T_d,

    minimize   max over x in [lam, 1] of | p(x) - f(x) |,   f(x) = log|x| / log(1/lam) + shift,
    subject to | p(x) | <= 1 on [-1, 1]

(the second condition is what a QET polynomial must obey). `lam` is the smallest eigenvalue the
QET sees, lambda_min / s. The minimal degree for a tolerance is the smallest even d whose optimum
is at or below it. Both conditions are imposed on grids; the defaults are converged to four
digits in the error at the degrees quoted below (checked with 2.5x finer grids).

The chapter USES the shift of arXiv:2407.13080 Sec. IV, shift = 1/2 ('A shift of 1/2 works well.
This known shift can always be subtracted at the end of the calculation'): the target then sits in
the middle of [-1, 1] and the polynomial has room below it for |x| < lam. The tolerance is the one
that guarantees 1% on log det W at V = 4^4, 0.01 x 150.5528 / (256 ln 27) = 1.7844e-3.

    python3 ch04_roundB_lp_degree.py 0.037037037 1.7844e-3 0.5   -> d = 84    (the 2028 box)
    python3 ch04_roundB_lp_degree.py 0.037037037 1e-2 0.5        -> d = 48
    python3 ch04_roundB_lp_degree.py 0.037037037 1.7844e-3       -> d = 314   (no shift)
    python3 ch04_roundB_lp_degree.py 0.037037037 1e-2            -> d = 112   (no shift)
    python3 ch04_roundB_lp_degree.py 0.16 1e-2                   -> d = 36    (no shift)
    python3 ch04_roundB_lp_degree.py 0.037037037 1e-3            -> d = 510   (no shift)

Results of this pass (lam = 0.16/4.32):
    shift 1/2   d    46        48        82        84        86
                err  1.104e-2  9.939e-3  1.866e-3  1.699e-3  1.548e-3
    no shift    d    36       64       110      112      312       314       508       510
                err  4.36e-2  2.68e-2  1.03e-2  9.91e-3  1.795e-3  1.770e-3  1.006e-3  9.985e-4
The unshifted error is not geometric in d, so interpolating between two degrees is not reliable.
The shifted optimum does not depend on the size of the shift between 0.3 and 0.6.
"""
import math
import sys

import numpy as np
from scipy.optimize import linprog


def grids(lam, n_in, n_all):
    k = np.arange(n_in)
    x_in = 0.5 * (1 + lam) - 0.5 * (1 - lam) * np.cos(np.pi * k / (n_in - 1))   # [lam, 1]
    k = np.arange(n_all)
    x_cheb = np.cos(np.pi * k / (2 * (n_all - 1)))                               # [0, 1], dense at 1
    x_lin = np.linspace(0, 1, n_all)
    x_gap = np.linspace(0, lam, n_all // 2)                                      # below lam
    return x_in, np.unique(np.concatenate([x_cheb, x_lin, x_gap, x_in]))


def cheb_even(x, d):
    theta = np.arccos(np.clip(x, -1, 1))
    return np.cos(np.outer(theta, np.arange(0, d + 1, 2)))


def min_err(lam, d, n_in=None, n_all=None, shift=0.0):
    """Smallest achievable uniform error at even degree d; (error, Chebyshev coefficients)."""
    n_in = n_in or max(2000, 12 * d)
    n_all = n_all or max(3000, 12 * d)
    x_in, x_all = grids(lam, n_in, n_all)
    f = np.log(x_in) / math.log(1 / lam) + shift
    a_in, a_all = cheb_even(x_in, d), cheb_even(x_all, d)
    n = a_in.shape[1]
    cost = np.zeros(n + 1)
    cost[-1] = 1                                   # variables: coefficients, then the error t
    ones = np.ones((len(x_in), 1))
    zeros = np.zeros((len(x_all), 1))
    a_ub = np.block([[a_in, -ones], [-a_in, -ones], [a_all, zeros], [-a_all, zeros]])
    b_ub = np.concatenate([f, -f, np.ones(len(x_all)), np.ones(len(x_all))])
    res = linprog(cost, A_ub=a_ub, b_ub=b_ub, bounds=[(None, None)] * n + [(0, None)], method="highs")
    if res.status != 0:
        return float("inf"), None
    return res.x[-1], res.x[:-1]


def min_degree(lam, tol, shift=0.0, d_max=2048):
    """Smallest even degree whose optimum is <= tol (the optimum is non-increasing in d)."""
    lo, d = 0, 2
    while min_err(lam, d, shift=shift)[0] > tol:
        lo, d = d, 2 * d
        if d > d_max:
            raise ValueError(f"no degree <= {d_max} reaches {tol}")
    hi = d
    while hi - lo > 2:
        mid = 2 * ((lo + hi) // 4)
        if mid <= lo:
            mid = lo + 2
        if min_err(lam, mid, shift=shift)[0] <= tol:
            hi = mid
        else:
            lo = mid
    return hi


if __name__ == "__main__":
    lam, tol = float(sys.argv[1]), float(sys.argv[2])
    shift = float(sys.argv[3]) if len(sys.argv) > 3 else 0.0
    d = min_degree(lam, tol, shift)
    print(f"lam={lam:.6g} tol={tol:g} shift={shift:g}: d_min={d}  "
          f"err(d)={min_err(lam, d, shift=shift)[0]:.5g}  err(d-2)={min_err(lam, d - 2, shift=shift)[0]:.5g}")
