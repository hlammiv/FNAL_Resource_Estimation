# Ch. 8: the integrated reflection asymmetry S_int from the free-fermion comparator (2026-10-12)

`scatter.py` scatters Gaussian packets of the symmetric-phase positive-energy eigenstates of the chapter's two
2-component Wilson flavours (r = 1, a m_f = 0.5) off the C-violating tanh wall
|y| phi(z) e^{i theta_C(z)} psi_L^dag beta psi_R + h.c., prof(z) = (1 + tanh((z - z0)/width))/2, theta_C = Delta theta_C x
prof, on a long open chain (N = 600-800 sites) per k_y sector of the L_y = 8 cylinder, and reads the reflected
weight R(p) on the incident side. A_R(k) = R_L(k) - R_R(k): the C' symmetry maps the flavour-L antiparticle onto the
flavour-R particle, which `checks.py` confirms directly by scattering the flavour-L hole (R(hole_L) = R_R to 1e-4).
S_int = sum_k f_eq(k) A_R(k) / sum_k f_eq(k) over the incident modes (k_z = 2 pi n/15, n = 1..7; k_y = 2 pi m/8),
f_eq Fermi-Dirac at a T = 0.5 (beta m_phi = 1, a m_phi = 0.5) and at 0.9 T_c.

| |y| phi | width (a) | Delta theta_C | S_int (T_c) | S_int (0.9 T_c) | max per-bin |A_R| |
|---|---|---|---|---|---|---|
| 0.5 | 1 | 0.1 | +1.92e-02 | +2.09e-02 | 0.059 |
| 0.5 | 2 | 0.1 | +2.14e-02 | +2.35e-02 | 0.056 |
| 0.3 | 2 | 0.1 | +2.42e-02 | +2.68e-02 | 0.055 |
| 1.0 | 2 | 0.1 | +2.39e-02 | +2.57e-02 | 0.050 |
| 0.5 | 2 | 0.2 | +4.28e-02 | +4.70e-02 | 0.111 |
| 0.5 | 2 | 0.0 | 0 | 0 | 0.000 |

The chapter takes the |y| phi = 0.5, width 2 a row (S_int = 2.14e-2) and prints the band 1.9-2.7e-2 over both
temperatures. Per-bin A_R changes sign across the axis (+0.047, 0.001, 0, 0, -0.045, -0.029, -0.001 at k_y = 0), so the
thermal average is a partial cancellation; the middle bins transmit (R < 0.01).

**Sign fix (verification pass, 2026-10-12).** The first run launched packets from the eigenvectors of s3 H_Bloch s3
(the sin terms of `bloch()` had the wrong sign relative to `chain_h`), i.e. s3 u instead of u, which mixed in a
negative-energy component of weight 1 - (M/E)^2 that stayed on the incident side and scaled every bin's A_R by
(M/E)^2. That run gave S_int = 1.33e-2 (kept in `sint_results_before_signfix.json`); the corrected values above are
1.6-1.7 times larger bin by bin, with the same sign pattern, linearity and zero at Delta theta_C = 0. Converged against N (600 vs 1000) and packet width (sigma 25 vs 40 sites) to four digits except
at the zone-edge bin (slow, small weight). The comparator is a static prescribed wall in the free theory; the
interacting, plasma-dressed wall is what the chapter's calculation determines. `sint_results.json` holds every bin.
