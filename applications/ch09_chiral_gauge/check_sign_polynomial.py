#!/usr/bin/env python3
"""Numerical check of the Ch. 9 sign-function polynomial (ERROR_LEDGER_2033["sign_function"]).

The model prices d_sgn = ceil(kappa ln 1/delta_sgn) = 208 at kappa = alpha/Delta = 30 and delta_sgn = 1e-3.
QSVT needs a real odd polynomial p with |p(x)| <= 1 on [-1, 1] and |p(x) - sign(x)| <= delta on
1/kappa <= |x| <= 1 (the spectrum of the normalized kernel, gap Delta, norm alpha). This script:
  1. solves, by linear program in the odd Chebyshev basis, the minimax problem
       min_c max_{x in [1/kappa, 1]} |p(x) - 1|   subject to |p(x)| <= 1 on [0, 1],
     on a Chebyshev-clustered grid, at the priced degree;
  2. re-checks the result on a dense independent grid (error, boundedness, parity by construction);
     then divides p by its dense-grid maximum when that exceeds 1, so |p| <= 1 holds on the dense grid;
  3. bisects for the smallest odd degree that reaches delta_sgn under the same constraints.
d_sgn = 208 is even; an odd polynomial has odd degree, so the check runs at 209 and also at 207.
What it does not do: compute QSP phase factors, or check that the kernel's spectrum lies inside 1/kappa <= |x| <= 1.
Usage: python applications/ch09_chiral_gauge/check_sign_polynomial.py [--kappa 30] [--delta 1e-3] [--out f.json]
Runtime: a few seconds."""
import argparse, json, math
import numpy as np
from numpy.polynomial import chebyshev as C
from scipy.optimize import linprog


def grid(kappa, n):
    a = 1.0 / kappa
    t = np.cos(np.pi * (np.arange(n) + 0.5) / n)                  # Chebyshev nodes on [-1, 1]
    band = a + (1 - a) * (t + 1) / 2                                # approximation band [a, 1]
    gap = a * (t + 1) / 2                                           # [0, a]: boundedness only
    return np.sort(band), np.sort(gap)


def basis(x, deg):
    return np.stack([C.chebval(x, np.eye(deg + 1)[k]) for k in range(1, deg + 1, 2)], axis=1)


def minimax(kappa, deg, n_grid=None):
    n_grid = n_grid or max(1500, 8 * deg)
    band, gap = grid(kappa, n_grid)
    A, G = basis(band, deg), basis(gap, deg)
    m = A.shape[1]
    # variables: c (m), e ; minimize e
    cost = np.r_[np.zeros(m), 1.0]
    ones_b, ones_g = np.ones((len(band), 1)), np.ones((len(gap), 1))
    A_ub = np.vstack([np.hstack([A, -ones_b]),                      # p - 1 <= e
                      np.hstack([-A, -ones_b]),                     # 1 - p <= e
                      np.hstack([A, 0 * ones_b]),                   # p <= 1 on band
                      np.hstack([G, 0 * ones_g]),                   # |p| <= 1 on gap
                      np.hstack([-G, 0 * ones_g])])
    b_ub = np.r_[np.ones(len(band)), -np.ones(len(band)), np.ones(len(band)), np.ones(len(gap)), np.ones(len(gap))]
    r = linprog(cost, A_ub=A_ub, b_ub=b_ub, bounds=[(None, None)] * m + [(0, None)], method="highs")
    if r.status != 0:
        raise RuntimeError(r.message)
    coef = np.zeros(deg + 1)
    coef[1::2] = r.x[:m]
    return coef, r.x[-1]


def dense_check(coef, kappa, n=400_001):
    """Rescale to max |p| <= 1 on a dense grid of [0, 1] (odd p, so [-1, 0] follows), then measure the band error."""
    x = np.linspace(0, 1, n)
    p = C.chebval(x, coef)
    peak = float(np.max(np.abs(p)))
    scale = max(peak, 1.0)
    p = p / scale
    band = x >= 1 / kappa
    return dict(err_band=float(np.max(np.abs(p[band] - 1))), max_abs=float(np.max(np.abs(p))), lp_peak=peak)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--kappa", type=float, default=30.0)
    ap.add_argument("--delta", type=float, default=1e-3)
    ap.add_argument("--out")
    a = ap.parse_args()
    d_priced = math.ceil(a.kappa * math.log(1 / a.delta))
    d_odd = d_priced if d_priced % 2 else d_priced + 1                # an odd polynomial has odd degree
    coef, e_lp = minimax(a.kappa, d_odd)
    chk = dense_check(coef, a.kappa)
    chk_lo = dense_check(minimax(a.kappa, d_odd - 2)[0], a.kappa)
    lo, hi = 3, d_odd
    while hi - lo > 2:                                               # smallest odd degree reaching delta
        mid = (lo + hi) // 2
        mid += 1 - mid % 2
        if dense_check(minimax(a.kappa, mid)[0], a.kappa)["err_band"] <= a.delta:
            hi = mid
        else:
            lo = mid
    out = dict(kappa=a.kappa, delta=a.delta, d_priced=d_priced, d_checked=d_odd, lp_error=e_lp,
               lp_peak_before_rescale=chk["lp_peak"], dense_error=chk["err_band"], dense_max_abs=chk["max_abs"],
               dense_error_at_d_minus_2=chk_lo["err_band"],
               passes=chk["err_band"] <= a.delta and chk_lo["err_band"] <= a.delta and chk["max_abs"] <= 1.0,
               smallest_odd_degree_found=hi)
    for k, v in out.items():
        print(f"{k:28s} {v}")
    if a.out:
        with open(a.out, "w") as f:
            json.dump(out, f, indent=1)
    return 0 if out["passes"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
