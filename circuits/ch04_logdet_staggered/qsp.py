"""QSP phase factors for a real polynomial of definite parity (numpy and scipy only).

Two conventions appear, and both are defined here by their 2x2 products so that tests can check them.

Wx convention (d + 1 phases):
    U(x) = e^{i phi_0 Z} W(x) e^{i phi_1 Z} W(x) ... W(x) e^{i phi_d Z},   W(x) = e^{i arccos(x) X}
    target: Re <0|U(x)|0> = f(x), with the phases symmetric, phi_j = phi_{d-j}.

Reflection convention (d phases), the one a QSVT circuit applies (Gilyen, Su, Low, Wiebe 2019,
Corollary 8 and Definition 15):
    U_R(x) = prod_{j=1..d} [ e^{i phi_j Z} R(x) ],   R(x) = [[x, sqrt(1-x^2)], [sqrt(1-x^2), -x]]
    <0|U_R(x)|0> = P(x), Re P = f.
    Leftmost factor (j = 1) is applied last in a circuit.

The solver is the symmetric-phase least-squares problem of Dong, Meng, Whaley and Lin (2021), started
from (pi/4, 0, ..., 0, pi/4) and solved here with scipy.optimize.least_squares on the reduced phases.
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import least_squares


def wx_response(phases, x):
    """<0|U(x)|0> in the Wx convention, vectorized over x. Returns a complex array."""
    x = np.atleast_1d(np.asarray(x, dtype=float))
    s = np.sqrt(np.clip(1.0 - x * x, 0.0, None))
    # row vector <0| e^{i phi_0 Z} = (e^{i phi_0}, 0), then multiply on the right by W e^{i phi_j Z}
    a = np.exp(1j * phases[0]) * np.ones_like(x, dtype=complex)
    b = np.zeros_like(a)
    for phi in phases[1:]:
        a, b = a * x + b * (1j * s), a * (1j * s) + b * x
        a, b = a * np.exp(1j * phi), b * np.exp(-1j * phi)
    return a


def reflection_response(phases, x):
    """<0| prod_j e^{i phi_j Z} R(x) |0> in the reflection convention, vectorized over x."""
    x = np.atleast_1d(np.asarray(x, dtype=float))
    s = np.sqrt(np.clip(1.0 - x * x, 0.0, None))
    a = np.ones_like(x, dtype=complex)
    b = np.zeros_like(a)
    for phi in phases:
        a, b = a * np.exp(1j * phi), b * np.exp(-1j * phi)
        a, b = a * x + b * s, a * s - b * x
    return a


def wx_to_reflection(phases_wx):
    """d + 1 Wx phases -> d reflection phases with the same <0|.|0> element.

    R(x) = -i e^{i pi/4 Z} W(x) e^{i pi/4 Z}, so
    phi_1 = phi^W_0 + phi^W_d + (d - 1) pi/2 and phi_j = phi^W_{j-1} - pi/2 for j = 2..d.
    Phases are reduced to (-pi, pi]."""
    p = np.asarray(phases_wx, dtype=float)
    d = len(p) - 1
    out = np.empty(d)
    out[0] = p[0] + p[d] + (d - 1) * np.pi / 2.0
    out[1:] = p[1:d] - np.pi / 2.0
    return -((-out + np.pi) % (2.0 * np.pi) - np.pi)


def _full_from_reduced(red, d):
    """Symmetric phase vector of length d + 1 from its first ceil((d+1)/2) entries."""
    n = d + 1
    full = np.empty(n)
    half = (n + 1) // 2
    full[:half] = red
    full[n - half:] = red[::-1]
    return full


def solve_symmetric_phases(cheb_coeffs, tol=1e-12, max_nfev=400):
    """Symmetric Wx phases with Re <0|U(x)|0> = f(x) = sum_k c_k T_k(x).

    `cheb_coeffs` is the full Chebyshev coefficient vector c_0..c_d of a polynomial of definite parity
    (the coefficients of the other parity must vanish) with max |f| < 1 on [-1, 1].
    Returns (phases, residual) where residual is the largest |Re <0|U|0> - f| on a fine grid."""
    c = np.asarray(cheb_coeffs, dtype=float)
    d = len(c) - 1
    parity = d % 2
    if np.abs(c[1 - parity::2]).max(initial=0.0) > 1e-14:
        raise ValueError("polynomial does not have the parity of its degree")
    n_red = (d + 2) // 2
    k = np.arange(1, n_red + 1)
    nodes = np.cos((2 * k - 1) * np.pi / (4.0 * n_red))            # positive Chebyshev nodes
    target = np.polynomial.chebyshev.chebval(nodes, c)

    def resid(red):
        return wx_response(_full_from_reduced(red, d), nodes).real - target

    red0 = np.zeros(n_red)
    red0[0] = np.pi / 4.0
    if d == 0:
        red0[0] = np.pi / 4.0
    sol = least_squares(resid, red0, xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=max_nfev, method="lm")
    phases = _full_from_reduced(sol.x, d)
    grid = np.cos(np.linspace(0.0, np.pi, 4001))
    err = np.abs(wx_response(phases, grid).real - np.polynomial.chebyshev.chebval(grid, c)).max()
    if err > max(tol, 1e-9):
        raise RuntimeError(f"QSP phase solve did not converge: residual {err:.3e}")
    return phases, float(err)
