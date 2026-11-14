"""Checks of su2rt: N_CS normalization (winding-1 vacuum), Gauss law, energy conservation, equipartition."""
import numpy as np, math, sys, time, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from su2rt import *
rng = np.random.default_rng(0)
for N in (8, 12, 16):
    c = (N - 1) / 2.0
    x = np.indices((N, N, N)).astype(float) - c
    r = np.sqrt((x * x).sum(0)); R = 0.45 * N
    F = np.where(r < R, math.pi * r / R, math.pi)
    v = 2 * F * x / np.where(r > 0, r, 1)
    g = qexp(v)[:, None]
    U = np.zeros((4, 1, 3, N, N, N))
    for i in range(3):
        U[:, :, i] = qmul(g, qdag(sh(g, i, 1)))
    U1 = np.zeros_like(U); U1[0] = 1
    print(N, 'dN along geodesic to the winding-1 vacuum', dNcs_path(U1, U, nsub=64), 'plaq energy', plaq_energy(U))
U, e, F = thermalize(2, 6, 6.0, rng, n_ref=20)
print('ekin*beta per link', (0.5 * (e * e).sum((0, 2, 3, 4, 5)) / (3 * 216) * 6.0), '(expect 1: 2 of 3 components)')
print('gauss', np.abs(gauss(U, e)).max())
E0 = energy(U, e); U2, e2, F2 = leapfrog(U, e, 0.2, 200, F); print('rel drift dtau=0.2', (energy(U2, e2) - E0) / E0, 'gauss', np.abs(gauss(U2, e2)).max())
U, e, F = thermalize(16, 4, 2.5, rng, n_ref=2)
t = time.time(); run(4, 2.5, B=16, t_max=5, n_ref=1); print('s per unit time, 4^3 x 16', (time.time() - t) / 5)
