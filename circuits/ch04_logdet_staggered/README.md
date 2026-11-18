# Ch. 4 2028 benchmark as an explicit circuit: log det W, free staggered, V = 4^4

This directory holds a gate-level circuit for the 2028 benchmark of Chapter 4 (hybrid lattice QCD), with tests against the exact matrix and counts obtained by walking the circuit. The model it is compared with is `estimates/ch04_hybrid_lqcd.py` (`_model_2028`). Nothing outside this directory was changed.

## 1. What the circuit is

One shot of the benchmark, 18 qubits, built by `build.benchmark_circuit(variant)`:

1. **Block encoding** `U = D_s O_c O_A D_s` of `W/s` (arXiv:2407.13080 Sec. II.2, Figs. 1-3), with `W = M^dag M` the free staggered operator on a periodic 4^4 lattice, `m0 = 0.4`, `K = 1`, `s = 4.32`.
   - `D_s`: Hadamards on the 4-qubit column register.
   - `O_c` (shift oracle): controlled shifts by two sites on the coordinate registers, as Toffoli ladders with clean ancillas.
   - `O_A` (entry oracle): y-rotations on the flag qubit controlled on the top column bit, written as `rz` between Cliffords plus CNOTs.
2. **QET** of degree `d = 84`: 84 queries alternating `U` and `U^dag`, and after each query a projector-controlled phase `exp(i phi (2 Pi - 1))`, `Pi` the projector on flag and column register all zero. The AND of those 5 qubits is written onto the signal qubit by 4 Toffolis on 3 workspace qubits, the phase is an `rz` on the signal qubit, and 4 more Toffolis undo the AND.
3. **Hadamard test**: only the phases are controlled on the Hadamard-test ancilla (2 `rz` and 2 CNOT each). With the phases off the queries cancel in pairs, so no query needs a control. `<Z>` on the ancilla is `<n| f(W/s) |n>` for the site basis state `|n>` of the shot, with `f` the degree-84 polynomial.

The polynomial is the one the model assumes: the LP-minimal even polynomial for `log|x|/log(s/lambda_min) + 1/2` on `[lambda_min/s, 1]`, computed here by the repository's own script `estimates/apply_log/ch04_roundB_lp_degree.py` (`min_err`, loaded read-only). Its 84 phases are solved in `qsp.py` and stored in `polynomial_d84.json`.

Gate set of every circuit: `h, x, s, sdg, cx, ccx, rz`. Every arbitrary rotation is an explicit `rz` with its angle. No rotation in the benchmark is a multiple of pi/4, so `N_T_fixed = 0`.

Three variants are built. The first two are the two ends of the model's band. The third is not in the model.

| variant | `O_c` per query | `O_A` per query | phases |
|---|---|---|---|
| `literal` | 8 shifts, each a 4-control X (5 Toffolis), as Fig. 1 draws it | two singly-controlled `R_y` as Fig. 2 draws it | 84 projector phases |
| `merged` | 4 shifts, each a 3-control X (3 Toffolis). On L = 4 the +2 and -2 shifts are the same | one multiplexed rotation | 84 projector phases |
| `compact` | 6 Toffolis, controls shared between directions | one multiplexed rotation | 83 projector phases and one `rz` on the Hadamard-test ancilla |

## 2. How it maps to the chapter

| chapter (app12, 2028 benchmark) | here |
|---|---|
| free staggered, V = 4^4 periodic, m0 = 0.4, K = 1 | `build.BENCHMARK`, `build.Instance` |
| block encoding of W/s, s = 4.32 | `append_block_encoding`. `Instance.s_min` gives 4.32 and theta_0 = 0, theta_1 = 5.5086, as the model's `be_subnormalization` and `oa_angles` |
| "shift oracle O_c ... into Toffoli ladders" | `append_shift_oracle`, `mcx_ladder` |
| "entry oracle O_A and the projector-controlled phase shifts into synthesized rotations" | `append_entry_oracle`, `append_projector_phase`. The rotations stay as `rz` |
| QET of the normalized log shifted by 1/2, d = 84 | `append_qet`, `polynomial_d84.json` |
| trace by the Hadamard test | `hadamard_test_circuit` |
| 18 logical qubits: 8 site, 4 column, 1 flag, 1 signal, 3 workspace, 1 Hadamard | `registers`, the same six registers in that order |

Register order in the circuit and in the QASM files: `site[8]`, `col[4]`, `flag`, `sig`, `work[3]`, `had`. Coordinate mu holds `site[2 mu]` (low bit) and `site[2 mu + 1]`. The column register is direction bits `col[0..1]`, sign bit `col[2]`, top bit `col[3]` (0 for a hop, 1 for the diagonal), which is the paper's l = 0..7 for hops and 8..15 for the diagonal.

## 3. Files

| file | content |
|---|---|
| `build.py` | circuit generator. `Instance(L, D, m0, K, s, variant)`, the block builders, `hadamard_test_circuit`, `benchmark_circuit` |
| `qsp.py` | QSP phase solver for a real polynomial of definite parity, and the two phase conventions |
| `polynomial_d84.json` | LP coefficients and phases of the benchmark polynomial (`python build.py --regenerate` recomputes them) |
| `count.py`, `counts.json` | counts walked from the circuit, per block and total, next to the model's line items |
| `test_circuit.py` | pytest suite, 86 tests |
| `sim.py` | statevector helper (qiskit-aer, 4 threads) |
| `synthesis_sensitivity.py`, `synthesis_sensitivity.json` | response of the result to rotation-angle errors (Section 7) |
| `logdet_4x4x4x4_d84_{literal,merged,compact}.qasm` | OpenQASM 3 of the full benchmark shot, site 0, with the final measurement. 173 kB, 96 kB, 81 kB |

Run from the repository root:

```
OMP_NUM_THREADS=4 python -m pytest circuits/ch04_logdet_staggered/test_circuit.py -q
python circuits/ch04_logdet_staggered/count.py --qasm     # tables, counts.json, the three QASM files
```

Run this test file by itself. The sibling directories under `circuits/` use the same module names (`build`, `count`, `test_circuit`), and pytest stops at collection when two of them are given in one session.

Other sizes: `build.hadamard_test_circuit(build.Instance(L=8, D=4, variant="literal"), phases)` with any even-length list of reflection-convention phases (`qsp.wx_to_reflection(qsp.solve_symmetric_phases(cheb)[0])`). `L` is a power of two from 4 up and `D` is a power of two. `merged` and `compact` need L = 4, and `compact` needs D = 2 or 4.

## 4. What was verified, and at what size

All references are numpy objects built from the definitions in the test file, not from `build.py`: the matrix `M` of Eq. 16 on the periodic lattice, `W = M^T M`, the index map `c(n, l)` of Eqs. 4-8 with add2 and sub2, and `f(W/s)` by eigendecomposition of the dense `W`. Tolerance on every compared amplitude: `1e-9`. The worst deviation observed on the full circuit is `1e-13`.

| what | reference | size |
|---|---|---|
| **full benchmark circuit**, all three variants, sites 0, 37, 173 and 222 | on the Hadamard ancilla = 1 branch, the real part of the 256 amplitudes with every other ancilla zero equals the column `f(W/s)|n>`. On the ancilla = 0 branch the state is exactly the input. `<Z>` equals `Tr f(W/s)/V` | 18 qubits, d = 84, the full instance |
| log det from the noiseless `<Z>` | exact dense `log det W = 150.5528` (also the model's `logdet_W`) | circuit gives 150.1774, 0.25% low, inside the chapter's 1% polynomial guarantee |
| `O_c` | permutation `c(n, l)` on a random state of site and column registers, ancillas back in 0 | L=4 D=4 (three variants), L=4 D=2, L=8 D=2, L=16 D=1, L=8 D=4 (22 qubits) |
| `O_A` | flag amplitudes from the matrix entries, per column value | same instances, and s = 6.5 where theta_0 is not 0 |
| block encoding `U` and `U^dag` | `<0|U|0> = W/s` applied to 4 basis and 3 random site vectors. `U^dag U = 1` on a random state of all qubits | same instances |
| projector phase, controlled and not | diagonal `exp(i phi (2 Pi - 1))` on a random state, ancillas back in 0 | 18, 12 and 11 qubits |
| QET and Hadamard test with random even polynomials, d = 2 and 8 | `f(W/s)` column | L=4 D=2, L=4 D=4 at s = 6.5, L=8 D=2, L=16 D=1 |
| polynomial | within the model's tolerance 1.784e-3 of the shifted normalized log on a 200001-point grid (error 1.70e-3), bounded by 0.674, reproduced by the stored phases to 1e-12 | d = 84 |
| ladders and AND tree | exact permutation matrices, 2n - 3 Toffolis as `common.mcx_toffoli_count` | up to 6 controls |
| QASM export | parses back to the same gate sequence, angles to 1e-12 | all three files |

Tests that fail if a block is wrong. Each builds a faulty circuit and asserts that the same check that passes on the correct circuit sees the fault:

- hopping entry with the wrong sign (theta_1 replaced by 2 pi - theta_1): block check error above 1e-2.
- one hopping term dropped (one Toffoli onto the site register removed), all three variants: `O_c` and block checks fail.
- full d = 84 circuit with one QET phase sign flipped, and with one Toffoli of one `O_c` removed: column check error above 1e-3.

Result of the last run on this machine: `86 passed`, 30 to 60 s with 4 threads, depending on the load of the machine.

## 5. Counts

Counted by `count.py` from the circuit object (`walk_counts`, `t_depth`). Model values come from a live call of `estimates.ch04_hybrid_lqcd.model(Assumptions(), '2028')`: `.lq`, `.hard_ops`, and the intermediates `be_queries`, `oc_toffolis_per_query`, `oa_rotations_per_query`, `projector_phases`, `projector_toffolis_each`, `projector_rotations_each`, `rotations_per_shot`, `t_per_rotation`, `t_depth_per_shot`, `eps_syn`.

The totals of each variant (4032, 1680 and 1168 Toffolis, 336, 336 and 335 rotations, no exact T) are also asserted by formula in `test_gate_set_and_counts_are_walked_from_the_circuit`, with no Toffoli allowed outside `O_c` and `proj_and`. A regenerated `counts.json` thus cannot hide a change of the circuit.

T convention: `common.toffoli_t(n, "textbook")` (7 T), `common.t_per_rotation(common.eps_rot_for(N_rot, eps_syn))` with `N_rot` from the circuit. The chapter model departs from the report default: `eps_syn = 0.1 x` the polynomial tolerance `= 1.784e-4` per shot. That value is used for the comparison. The total at the default `EPS_SYN = 1e-2` is listed too.

| line item | circuit, literal | model, literal | ratio | circuit, merged | model, merged | ratio |
|---|---|---|---|---|---|---|
| logical qubits (register) | 18 | 18 | 1.000 | 18 | 18 | 1.000 |
| qubits touched by a gate | 14 | not a model item | | 14 | | |
| `O_c` Toffolis | 3360 | 3360 | 1.000 | 1008 | 1008 | 1.000 |
| `O_A` rotations | 168 | 168 | 1.000 | 168 | 168 | 1.000 |
| projector Toffolis | 672 | 680 | 0.988 | 672 | 680 | 0.988 |
| projector rotations | 168 | 170 | 0.988 | 168 | 170 | 0.988 |
| N_Toff | 4032 | 4040 | 0.998 | 1680 | 1688 | 0.995 |
| N_rot | 336 | 338 | 0.994 | 336 | 338 | 0.994 |
| N_T_fixed | 0 | 0 | | 0 | 0 | |
| T per rotation | 21.186 | 21.191 | 1.000 | 21.186 | 21.191 | 1.000 |
| **T per shot** (eps_syn = 1.784e-4) | **35342.4** | 35442.4 | 0.997 | **18878.4** | 18978.4 | 0.995 |
| T per shot (eps_syn = 1e-2) | 34220.2 | not used by the model | | 17756.2 | | |
| T-depth, commuting order | 7455.4 | 7417.4 | 1.005 | 7455.4 | 7417.4 | 1.005 |
| T-depth, as written | 10899.4 | | 1.469 | 7455.4 | | 1.005 |

Per block, per shot (literal / merged / compact):

| block | Toffoli | arbitrary rz | CNOT | 1-qubit Clifford |
|---|---|---|---|---|
| `D_s` | 0 | 0 | 0 | 840 |
| `O_A` | 0 | 168 | 168 | 672 |
| `O_c` | 3360 / 1008 / 504 | 0 | 0 / 0 / 168 | 1176 / 504 / 168 |
| `proj_and` | 672 / 672 / 664 | 0 | 0 | 840 / 840 / 830 |
| `proj_rot` | 0 | 168 / 168 / 166 | 168 / 168 / 166 | 0 |
| `had` | 0 | 0 / 0 / 1 | 0 | 2 |

The `compact` variant has 1168 Toffolis, 335 rotations, 15272.4 T per shot and T-depth 7407.2. That is 0.80 of the model's lower end. It is a different compile, not a correction to the model.

**T-depth convention.** The model's: one layer per Toffoli, `t_per_rotation` layers per synthesized rotation, and gates that commute may overlap, so `O_c` runs under the `O_A` rotations. `count.t_depth(qc, t_rot, "commuting")` is the weighted longest path under that rule. Two gates are taken to commute when on every shared qubit both are diagonal or both are X-type. `"as_written"` orders every pair of gates that share a qubit.

**Verdict on the comparison.** No T or rotation line item differs by more than 1.2%. The model's count is confirmed for both ends of its band. The circuit is 100 T below the model at each end because QSVT of degree d has d projector phases, not d + 1 (Section 6, item 1).

One line item differs by more than 20%: **qubits touched by a gate, 14 against the 18-qubit register**. Both numbers are right and they answer different questions. The register is 18, as the model says, and the circuit and QASM files allocate 18. But free `W` only couples sites two hops apart, so the low bit of each of the four coordinates is never the target or control of any gate. At V = 4^4 the matrix is 16 identical blocks, each `(m0^2 + 2K^2) - (K^2/2) sum_mu X_mu` on the four high coordinate bits (test `test_low_coordinate_bits_are_never_acted_on`). The benchmark as specified can run on 14 logical qubits with the same T-count. This holds for free `W` at any even L, and stops holding once links are loaded.

## 6. Choices the circuit forced that the model leaves open

1. **Number of phases.** The model counts d + 1 = 85 projector phases (the QSP convention with a signal rotation at both ends). A QSVT circuit of even degree d applies d = 84. In a Hadamard test the last of them acts on a state that is projected onto `Pi` by the overlap, so it can be one `rz` on the ancilla: 83 projector phases. `literal` and `merged` use 84. `compact` uses 83.
2. **Which convention the phases are in.** The solver works in the Wx convention (85 symmetric phases). The circuit needs the reflection convention (84 phases). `qsp.wx_to_reflection` converts, and the full-circuit test checks the chain.
3. **Controlling `U_Phi`.** Only the phases are controlled, each as `rz(phi) cx rz(-phi) cx` on the signal qubit. This is the model's assumption (2 rotations per phase) and it works because d is even and the queries alternate `U`, `U^dag`.
4. **Uncomputing the AND.** The circuit uncomputes with 4 Toffolis, so it is unitary and can be simulated. The model's count also charges 8 Toffolis per phase, but its T-depth (`and_tree_depth_2028 = 3`) assumes a measurement-based uncompute of depth 0. In the circuit 2 compute layers and 2 uncompute layers per phase lie on the critical path, which is why the depth is 7455 and not 7417. A measurement-based uncompute would remove 4 of the 8 Toffolis per phase from the T-count, which the model does not do.
5. **Polarity of the top column bit.** `O_A` and `O_c` both test the top column bit for 0. With X gates around each negative control the two oracles do not commute gate by gate, and a schedule that keeps the gate order runs `O_c` after `O_A`. The "as written" row shows the size of that effect, 10899 for the literal variant. The circuit complements the top bit once per query, next to the Hadamards, and both oracles use it as a positive control. Then they commute and the model's "`O_c` runs underneath" holds. `O_c` literal is 40 Toffoli layers and the two `O_A` rotations are 42.4 layers, so the margin is thin at this synthesis tolerance.
6. **`O_A` at theta_0 = 0.** At `s = 4.32` the diagonal rotation is the identity and is not emitted. The two remaining rotations are `R_y(+-theta_1/2)`, `theta_1/2 = 2.7543`. At generic `s` the literal form emits 4 rotations and the multiplexed form 2 (tested at s = 6.5), which is the model's `oa_rotations_2028 = (2, 4)`.
7. **`R_y` as `rz`.** `R_y(theta) = S H R_z(theta) H S^dag`.
8. **Adder for the shift.** On L = 4, add2 is one multi-controlled X on the high coordinate bit. For larger L the builder uses the cascade of multi-controlled X that the model's `shift2_toffolis` assumes, and the Toffoli counts agree with that function at L = 8 and 16 (test `test_shift_oracle_toffolis_follow_the_models_formula`).
9. **Workspace.** 3 clean qubits. The projector AND uses all 3, the literal shifts use 2 of them.
10. **Input state.** A site basis state per shot, drawn at random classically to estimate the trace. For the free operator every site gives the same value, so the exported QASM uses site 0.
11. **Polynomial normalization.** The LP polynomial is bounded by 0.674 on [-1, 1], so the phases exist without rescaling and the solver reaches a residual of 6e-15.
12. **Column register for other D.** The builder needs D a power of two, so that the column register splits into direction bits, a sign bit and a top bit.

## 7. What was not verified

- **Rotation synthesis.** The `rz` gates are not synthesized into Clifford+T. `1.15 log2(1/eps) + 9.2` is the model's formula applied to the circuit's `N_rot`. No gridsynth or repeat-until-success circuit was built or counted.
- **The synthesis ancilla.** The formula above is the repeat-until-success fit, and `estimates/common.py` records one ancilla held while such a rotation is in flight (`RUS_ANCILLA = 1`). The circuit keeps `rz` gates and allocates none, and the model's 18 does not include it either. With that synthesis the register is 19 while a rotation runs.
- **Phases close to a Clifford angle.** In the reflection convention most QET phases sit near -pi/2. Of the 84, 12 are within `eps_rot = 7.29e-4` rad of a multiple of pi/4 and 24 are within `2 eps_rot`, the angle shift used in the sensitivity study below. The closest is 2.7e-4 rad away. All of them are counted as arbitrary rotations here. A synthesizer run at this tolerance should be checked on these angles.
- **The synthesis error budget**, only partly. `synthesis_sensitivity.py` shifts every `rz` angle of the `merged` benchmark circuit by `+-2 eps_rot` with random signs (an operator-norm error of `eps_rot = 7.29e-4` each), 24 sign patterns. For one fixed pattern `<Z>` moves by 4.6e-3 rms, which is 2.6 times the polynomial tolerance. Averaged over patterns the shift is -3.9e-5 (standard error 0.5e-5), 0.22 of the model's `eps_syn = 1.78e-4`. So the chapter's 0.1% synthesis error holds if the rotation errors are redrawn from shot to shot (randomized synthesis, as the chapter states). It does not hold for one fixed deterministic synthesis: 4.6e-3 on the normalized trace is 2.6% of log det W. This is an angle-error model with independent signs, not a simulation of synthesized gate sequences. A fixed synthesis repeats the same error on the 84 copies of each `O_A` angle, which this model does not capture.
- **Minimality of d = 84.** The LP was rerun at d = 84 only (error 1.6993e-3, the model's `d_log_err`). The repository test `test_2028_lp_degree_is_minimal` covers d = 82.
- **Sampling.** No shots were drawn. The tests compare exact expectation values. The 3.1e5 shots, the wall time and the logical error rate of the model are not touched here.
- **Layout, routing, factories.** The circuit is logical-level with all-to-all connectivity. The T-depth is a longest-path count under the stated weights, not a schedule.
- **The full QET at L > 4.** `O_c`, `O_A` and the block encoding are checked at L = 8 (D = 2 and D = 4) and L = 16 (D = 1). The QET is checked at L = 8, D = 2 and L = 16, D = 1 with d = 2 and 8 only.
- **The gauged circuit.** Nothing here applies to the 2033 or post-2033 rows, which have no compiled link oracle.

Tools that would help and are not installed: gridsynth (to replace each `rz` by a Clifford+T sequence and count T directly), Qualtran or pyLIQTR (an independent count of the same block encoding). pyzx and stim are installed and were not needed.
