# Ch. 8: the integrated reflection asymmetry S_int from the free-fermion comparator (2026-10-12)

**The 15-site axis (2026-11-09, `scatter15.py`; referee SR-ch8-1).** The comparator below scatters on an 800-site open chain (packet at N/4, wall at N/2, each bin read at t = 2(z_wall - z_pk)/v_g), which the chapter's 15 x 8 instance cannot do: there the wall sits at z = 7 of 15, the packet starts 4.5 sites before it, and the transmitted part reaches the open far end and returns inside t_max = 10/m_phi = 20 a. `scatter15.py` reruns the same Hamiltonian, wall and packets on the 15-site axis and reads every bin on one time grid; the chapter and the model (`s_int`) take its value at t_max, the end of the priced evolution, which is what the first-result circuit reads. Section 'Results on the 15-site axis' below; `sint15_results.json` holds every bin at every time.

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

## Results on the 15-site axis (`scatter15.py`, 2026-11-09)

Same Hamiltonian, wall (|y| phi = 0.5, width 2 a, Delta theta_C = 0.1) and packet construction as `scatter.py`, on N = 15 sites with the wall at z0 = 7.0 (sites 0-6 incident), the packet at z = 2.5 with sigma_z = 1.5 sites (2.0 as a check), and the incident-side weight read at t = 0, 1, ..., 20 a. The 7 incident bins are k_z = 2 pi n/15, n = 1..7 (k_z = 0 has v_g = 0 and is not an incident packet; the chapter printed 8 bins until 2026-11-09), times 8 k_y sectors = 56 bins. A 1.5-site packet has sigma_k ~ 0.33, close to the bin spacing 0.42, so a 15-site bin is a momentum band, not a sharp k_z; the chapter's Eq. (deltan) says so.

S_int(t), f_eq-weighted over the 56 bins (sigma_z = 2.0 in brackets):

| t (a) | t (1/m_phi) | S_int at T_c | at 0.9 T_c |
|---|---|---|---|
| 10 | 5 | 5.6e-3 (6.5e-3) | |
| 12 | 6 | 8.5e-3 (9.5e-3) | |
| 14 | 7 | 1.14e-2 (1.24e-2) | |
| 16 | 8 | 1.40e-2 (1.49e-2) | |
| 18 | 9 | 1.55e-2 | |
| 20 = t_max | 10 | **1.57e-2 (1.59e-2)** | 1.69e-2 (1.71e-2) |

The model takes **S_int = 1.6e-2** (`s_int`), the t_max value, stable against the packet width; `s_int_readout_band` = (0.6e-2, 1.6e-2) records the sweep over readout times 10-20 a. By t_max the transmitted part of the fast bins has returned from the open end (free-chain R0 below), which is what the priced circuit would also see. The 800-site value 2.14e-2 (`s_int_long_chain`) is never reached on this axis. Over wall profiles at t_max (sigma_z = 1.5): 1.47e-2 at width 1 a, 0.82e-2 at |y| phi = 0.3, 1.92e-2 at 1.0, 0 at Delta theta_C = 0. Largest per-bin |A_R| at t_max 0.027 (0.028 at sigma_z = 2.0), p0 A_R = 0.018, 2.7 times below the depth prior 0.048 (the model's `dn_estimate` keeps the 800-site 0.056, the tighter check of the prior).

Per-bin readout windows at k_y = 0 (t_clear: first time the free packet has left the incident side, R0 < 0.2; t_last: last time before the returning part raises R0 again), with A_R at t_last and at t_max:

| k_z | v_g | f_eq (T_c) | window (a) | A_R at t_last | A_R at t_max | 800-site A_R |
|---|---|---|---|---|---|---|
| 0.419 | 0.85 | 0.194 | 10-17 | +0.0267 | +0.0257 | +0.047 |
| 0.838 | 1.00 | 0.097 | 6-13 | +0.0115 | +0.0132 | +0.001 |
| 1.257 | 0.94 | 0.045 | 6-15 | +0.0027 | +0.0056 | 0.000 |
| 1.676 | 0.79 | 0.022 | 7-17 | -0.0023 | -0.0005 | 0.000 |
| 2.094 | 0.60 | 0.013 | 10-20 | -0.0156 | -0.0156 | -0.045 |
| 2.513 | 0.37 | 0.008 | 19-20 | -0.0210 | -0.0210 | -0.029 |
| 2.932 | 0.13 | 0.007 | never inside t_max | | -0.0111 | -0.001 |

The campaign's per-bin estimate reads each bin inside its window (the two slowest bins, 1.5% of the thermal weight, have none inside t_max and are read at t_max); the first result reads the f_eq-weighted superposition once, at t_max. The 800-site chain remains the converged comparator of Fig. ch08_asymmetry; the figure marks the 15-site t_max value on panel (c).

## Readout study (`window_study.py`, 2026-11-11; referee follow-up F2)

The referee asked whether the first result's single readout at t_max, after the fast bins' transmitted parts have
returned from the open far end, still measures reflection, whether the first result and the campaign need separate
signal calibrations, and whether the flag efficiency p0 holds at the readout time. `window_study.py` answers all three
on the comparator above (same wall, packets and 56 bins); `window_study.json` holds every bin.

**Far-end returns are a 1% effect.** The same wall and packets on an axis padded with 45 sites at both ends (nothing
returns inside 20 a; 80 sites gives the same to five digits) separate the parts of the 15-site value at T_c:

| quantity | S_int |
|---|---|
| 15-site axis at t_max (what the first result reads; the model's `s_int`) | 1.568e-2 |
| padded axis at t_max (no returns) | 1.556e-2 |
| far-end return contribution | +1.2e-4 (0.8%) |
| padded axis at 60 a (scattering complete, these packets) | 1.647e-2 |
| 800-site chain, sharp k (`s_int_long_chain`) | 2.14e-2 |
| campaign estimator (each bin at its t_last, slowest bins at t_max) | 1.499e-2 |

So the first-result signal is reflection from the wall to within 1%. What separates it from the converged 2.14e-2 is
the 1.5-site packets' momentum band (23%) and incomplete scattering at t_max (5%), not the boundary. An absorbing
boundary or a longer axis is not needed for the first result. At 0.9 T_c the return part is 2e-6.

**Separate calibrations.** First result: 1.568e-2 at t_max, the value the model prices (9 / (p0^2 S_int^2 m_top)).
Campaign: each bin is read in its own window; the per-bin query count 1 / (p0^2 eps_C^2 m_top) is set by the absolute
precision eps_C = 3e-3, not by the signal, so the campaign's cost does not depend on the readout calibration. Its
windowed estimator 1.499e-2 is recorded, not priced.

**Flag efficiency at the readout time.** The model's a = 1/2 + p0 A_R/2 takes the vacuum charge fluctuation on the
incident half as independent of the packet. The exact many-body state here is the wall vacuum (every negative-energy
orbital of all 8 k_y sectors) plus the packet's positive-energy part (at least 97% of its weight), evolved to the
readout time; the incident-half charge distribution is its full counting statistics. The effective efficiency
p_eff = (a_part + a_anti - 1) / A_R is 0.669 for the first result at t_max and 0.662 for the campaign readout (per bin
0.640-0.675 where |A_R| > 5e-3), against the model's free-vacuum p0 = 0.658 and the wall-vacuum 0.644 the same code
gives with no packet. The model's p0 is thus within 2% and on the conservative side; nothing is repriced.
