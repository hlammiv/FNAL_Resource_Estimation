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
| 0.5 | 1 | 0.1 | +1.16e-2 | +1.27e-2 | 0.055 |
| 0.5 | 2 | 0.1 | +1.33e-2 | +1.46e-2 | 0.051 |
| 0.3 | 2 | 0.1 | +1.50e-2 | +1.67e-2 | 0.045 |
| 1.0 | 2 | 0.1 | +1.41e-2 | +1.52e-2 | 0.039 |
| 0.5 | 2 | 0.2 | +2.65e-2 | +2.92e-2 | 0.103 |
| 0.5 | 2 | 0.0 | 0 | 0 | 0 |

The chapter takes the |y| phi = 0.5, width 2 a row (S_int = 1.33e-2) and prints the band 1.2-1.7e-2. Per-bin A_R
changes sign across the axis (+0.031, 0.001, 0, 0, -0.038, -0.027, -0.002 at k_y = 0), so the thermal average is a
partial cancellation. Converged against N (600 vs 1000) and packet width (sigma 25 vs 40 sites) to four digits except
at the zone-edge bin (slow, small weight). The comparator is a static prescribed wall in the free theory; the
interacting, plasma-dressed wall is what the chapter's calculation determines. `sint_results.json` holds every bin.
