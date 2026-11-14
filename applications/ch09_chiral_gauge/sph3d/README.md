# Ch. 9 rate instance: classical real-time SU(2) Chern-Simons diffusion (2026-11-14)

Classical-statistical inputs for the 2033 rate instance in `estimates/ch09_chiral_gauge.py` (`sph_rate_point`,
`SPH3D_CLASSICAL`): pure SU(2), Kogut-Susskind Hamiltonian, canonical ensemble at beta_L = 4/(g^2 a T).

- `su2rt.py`: lattice code (leapfrog, Gauss-law projection, gradient-flow cooling, naive Kubo and cooled N_CS).
- `tests_su2.py`: N_CS normalization on a winding-1 vacuum (0.81, 0.92, 0.95 on 8^3, 12^3, 16^3), Gauss law, energy drift.
- `drive.py N beta B t_max dt_meas seed [depth] [tag]`: one run; `jobs.txt` lists the runs (on a 32-core host, about 1 h).
- `analyze.py` -> `kappa_table.json`: kappa = Gamma/(alpha T)^4, tree and one-loop a, jackknife errors, per (N, beta_L).
- `summary.py` -> `summary.json`: finite-volume ratios, topological share (`zfactor.py`), comparison with
  Moore-Rummukainen hep-ph/9906259 Tables 1-2 at beta_L = 8.6.
- `classical_inputs.py N beta`: the S(t) and Re G(chi, t) that `SPH3D_CLASSICAL` holds (4^3, beta_L = 3).
- `quench.py`: equilibration from U = 1 (product state), the thermalization assumption.

Run data (`runs/`, about 70 MB) are not committed; `jobs.txt` regenerates them.
