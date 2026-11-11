# Ch. 3: gap, readout and 4He data-size inputs (recovered 2026-11-11)

The Ch. 3 model marks these fields COMPUTED: the 0+ gaps of the truncated KB3G Hamiltonians for 48Ca and 48Ti
(`gap_0plus_*`), the 27Al USDB gaps (`gap_al27_usdb_*`), the 27Al readout inputs (`sd_spin_expectation`,
`sd_element_expectation`, `direct_readout_variance`, `projection_success`) and the 4He data size D (`he4_d_*`).
They were computed on 2026-10-05 and 2026-10-09 in scratch folders that no longer existed when the
2026-11-11 referee follow-up (F1) found the README pointing at absent files. The drivers here were rebuilt from the
session record of those runs: each file's original write plus every later patch, replayed in order. Nothing was
rewritten by hand except the import paths, which now point at this repository.

| file | what it computes |
|---|---|
| `sm.py` | M-scheme shell model for KSHELL `.snt` interactions, proton-neutron factorized, Lanczos with a J^2 shift |
| `run_gaps.py` | lowest states at fixed M with a (J^2 - J0(J0+1))^2 penalty; writes `gap_<snt>_<space>_<nuc>_J<J0>.json` |
| `run_gaps_proj.py` | the same with a J-projected start vector (used to check the 48Ti full-pf run; not needed to reproduce) |
| `al27_sd.py` | 27Al USDB ground state (J = M = 5/2): <S_p>, <O_SD>, Var(O_SD), reference overlaps; writes `al27_sd_readout.json` |
| `he4D.py` | 4He in the 40-mode space with the Minnesota NN force and a Lawson term; D at fidelity and dipole targets |
| `reproduce.sh` | runs all of the above (fetch the interactions first with `../kb3g/fetch_snt.sh`) |
| `compare.py` | checks every output against the value the model holds; exits nonzero on any difference |

Regeneration on 2026-11-11 (numpy 2.5, scipy 1.18, sympy for Clebsch-Gordan; `reproduce.log`): all 20 checked values
match the model, the gaps to 1e-6 MeV, the readout inputs to 1e-5, and D exactly. The KB3G 48Ti full-pf case
(M-dimension 634,744) takes about an hour on four threads; everything else takes minutes. `sm.py` needs numpy >= 2
(`np.bitwise_count`).

`../kb3g/sm.py` is an earlier copy of the same code used by the lambda/|M| ladder; it is left unchanged.
