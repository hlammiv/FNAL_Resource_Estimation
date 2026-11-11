# Ch. 7: the 2026-10-04 step calibration (exact_over_proxy_step_2033 = 1.23), recovered 2026-11-11

The chapter's 2033 step error, '~2% (scaled from a reduced instance)', is the 5^3 mean-field proxy error
(`proxy_step_err_2033` = 0.016) times `exact_over_proxy_step_2033` = 1.23. The 2026-11-11 referee follow-up (F3)
found that the scripts behind 1.23 were not in the repository, and that `../full_step_check.py` uses its own proxy
definition and read grid, so it does not regenerate 1.23. These are the original scripts, rebuilt from the session
record of the 2026-10-04 run (each file's write and its later in-place patch, replayed in order); only paths changed.

| file | what it does |
|---|---|
| `exact_step.py` | exact many-body evolution under the chapter's splitting on a 3-site ring (K = 16, one Wilson fermion, Hilbert dim 81,920), `fixed` and `chi` (comoving) registers, at m_phi dt = 0.05, 0.025, 0.0125, read at t = 12.5, 25, 50, 100, 200, 400, 700, 1000 |
| `proxy_red.py` | the mean-field proxy (the chapter's `verlet_proxy` construction) on the same ring, same splitting and readout |
| `proxy_full.py` | `verlet_proxy` on the 5^3 instance at the same three steps and a dt = 0.001 reference |
| `analyze.py` | Richardson extrapolation per read time; writes `summary.json` |
| `ratio.py` | the calibration as computed on 2026-10-04: per coupling set, the largest exact (chi) step error over read times divided by the largest ring-proxy step error |
| `run.sh` | runs all of it (about an hour on 12 cores) |

Definition. For each (lambda, g) in (1, 1), (3, 0.5), (1, 2): exact = max over read times of
|n(dt) - n_R| / n_R for the fermion shell k = 2 pi/(3b), with n_R the Richardson value from dt, dt/2, dt/4;
proxy = the same on the ring proxy. `exact_over_proxy_step_2033` is the largest ratio over the three sets.

Rerun on 2026-11-11 (`ratio.json`; numpy 2.5, scipy 1.18), identical to the 2026-10-04 record:

| (lambda, g) | exact | ring proxy | ratio | 5^3 proxy | per-set extrapolation |
|---|---|---|---|---|---|
| (1, 1) | 1.38e-3 | 1.42e-3 | 0.973 | 0.0146 | 0.0142 |
| (3, 0.5) | 3.05e-3 | 2.48e-3 | 1.230 | 0.0033 | 0.0040 |
| (1, 2) | 1.14e-3 | 2.20e-3 | 0.518 | 0.0163 | 0.0085 |

The model's 2% is the largest 5^3 proxy error times the largest ratio, which come from different coupling sets;
the largest per-set product is 1.4%. So 2% is an upper envelope over the three sets, not a single-set value.
What the calibration does not establish: that the ratio measured on a 3-site ring carries over to 5^3. That
extrapolation is the stated assumption behind 'scaled from a reduced instance'.

An earlier copy of `exact_step.py` without the 2026-10-04 20:08 patch (the kinetic term of the ground-state
Hamiltonian summed over sites) gives ratios 0.871, 1.019, 0.440; the patch is part of the recorded run.
