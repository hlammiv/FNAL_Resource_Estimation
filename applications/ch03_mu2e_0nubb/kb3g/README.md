# Ch. 3: lambda/|M| and p_ref for the 48Ca -> 48Ti ladder with KB3G (ruling 2026-10-12)

The chapter's 0nubb validation ladder needs, for each truncated pf space, the ratio lambda/|M| of the Pauli-LCU
1-norm of the two-body operator to its ground-state-to-ground-state matrix element, and the overlap p_ref of the
closed-f7/2 reference with the 48Ca ground state. Until 2026-10-12 these came from a schematic-interaction toy
(a delta force fit to the 42Ca 0+ -> 2+ splitting). They are now computed with the same KB3G interaction that
prices the 0+ gaps, by exact M-scheme Lanczos in each kept space (no renormalization when orbits are dropped),
with GXPF1A as the interaction spread. The operator is unchanged: the toy's sigma.sigma tau+ tau+ with a 1/r12
radial form on harmonic-oscillator orbitals (b from hw = 41 A^-1/3), so lambda is the toy's (19.05 / 47.50 /
208.36) and only |M| and p_ref move.

| space | dims (48Ca, 48Ti, M = 0) | lambda/\|M\| KB3G | GXPF1A | toy (printed before) | p_ref KB3G | toy |
|---|---|---|---|---|---|---|
| f7/2 p3/2 | 57, 5296 | 197.04 | 145.69 | 69.1 | 0.9825 | 0.93 |
| f7/2 p3/2 p1/2 | 325, 24453 | 445.44 | 355.34 | 164 | 0.9709 | 0.90 |
| full pf | 12022, 634744 | 1403.9 | 1650.1 | 600-720 | 0.9022 | - |

The KB3G ground states sit closer to closed f7/2 and share less of O|48Ca> with the 48Ti ground state (ground-state
share of ||O|48Ca>||^2: 1.3% / 1.5% / 0.4% against the toy's ~10%).

Rabi-readout check (`rabi_kb3g.py`, the toy's protocol: evolve the reference under H_p - dE + detuning, H_n and
theta (O + O^dag) for a quarter period, scan the detuning, fit P_daughter = A sin^2(W t)): with KB3G in f7/2 p3/2 the
fit recovers |M| to 1 / 4 / 7 / 15 % at theta|M| = 0.1 / 0.2 / 0.3 / 0.5 MeV, the resonance shifted by 0.10 / 0.40 /
0.75 / 2.25 MeV. The toy's "6% at <= 0.5 MeV" had |M| 2.9 times larger, so the same theta|M| was a 2.9 times weaker
drive. In f7/2 p3/2 p1/2 (the printed row's space) the same check gives 4% and 10% at 0.3 and 0.5 MeV (resonance
shifted by 0.60 and 2.00 MeV; `rabi_kb3g_f7p3p1.log`, 2.7 h of which 2.6 h built the sparse Hamiltonian column by
column). The chapter keeps the 0.3-0.5 MeV band and prints both.

## Files

- `sm.py`: minimal pn-factorized M-scheme shell-model code reading KSHELL `.snt` files (the referee-M3 gap code).
- `toyop.py`: the simplobs toy (operator, delta-force Hamiltonian, Pauli-LCU lambda), reconstructed from the session
  record; `run(["f7","p3"])` reproduces lambda/|M| = 69.1 and p_ref = 0.929.
- `kb3g_lam.py <snt> <space> [k] [tol]`: the computation; writes `lam_<snt>_<space>.json`. The vectorized route
  applies the operator in the interaction file's mode numbering; the dict-based cross-check it prints is valid
  only when the toy's orbit order matches the file's (it does in f7p3 and f7p3p1, not in pf), see the log.
- `rabi_kb3g.py <snt> <space> <theta|M| list>`: the Rabi check (logs `rabi_kb3g_f7p3*.log`).
- `implications.py [snt]`: pushes the JSON results through `estimates/ch03_mu2e_0nubb.py` and prints the rows
  before and after.
- `fetch_snt.sh`: downloads `kb3g.snt` and `gxpf1a.snt` from the KSHELL distribution (not mirrored).

Full pf needs numpy >= 2 (`np.bitwise_count`) and about 1.5 GB; the 48Ti Lanczos took 31 min (KB3G) and 15 min
(GXPF1A) on one core.
