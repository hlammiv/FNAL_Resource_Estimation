# Ch. 5, 2028 benchmark: one Trotter link-step on the binary tetrahedral group, as an explicit circuit

This directory holds a gate-level circuit for the 2028 benchmark of Chapter 5 (transport in the quark-gluon plasma): 2+1D pure-gauge SU(2) digitized on the binary tetrahedral group 2T (24 elements, 5 qubits per link), Kogut-Susskind Hamiltonian, 4x4 periodic lattice (32 links, 16 plaquettes). It builds the group primitives, one electric term per link and one magnetic term per plaquette, and from them one Trotter step. The counts are read off the circuit object by code and set next to the line items of the chapter model, `estimates/ch05_qgp_transport.py`.

Headline: the circuit needs 1.80e4 T for the one-step first result and 1.12e5 T for the six-step full benchmark, against 9.47e4 and 5.94e5 in the model. The ratio is 0.19 in both. The reasons are listed in section 6. The circuits here are derived from the group law. They are not transcriptions of the published 2T circuits that the model quotes.

## 1. Files

| file | content |
|---|---|
| `build.py` | the generator. Group tables, the primitives, `plaquette`, `electric_link`, `trotter_circuit`, `single_plaquette_step`, `lattice_step(L, ...)`, `benchmark_circuit(tier)`, and the walkers `count_gates`, `count_by_block`, `weighted_depth`. Run as a script it writes the three QASM files. |
| `test_circuit.py` | pytest checks against numpy references (section 4). |
| `verify_independent.py` | a second checker that imports nothing from `build.py`. It builds 2T from quaternion arithmetic, reads the QASM text with its own parser and runs it with its own simulators (section 4). Added by the independent review. Run as a script it checks the three QASM files and recounts them. |
| `count.py` | walks the benchmark circuits, imports the model, writes `counts.json`, prints the tables. |
| `counts.json` | every count below, per block and total, with the model values and ratios. |
| `bt_4x4_first_result_1step.qasm` | OpenQASM 3, the 4x4 lattice, one step (176 qubits, 10,080 gates, 0.16 MB). |
| `bt_4x4_full_benchmark_6steps.qasm` | OpenQASM 3, the 4x4 lattice, six steps (60,480 gates, 0.97 MB). |
| `bt_single_plaquette_1step.qasm` | OpenQASM 3, one open plaquette, one step (22 qubits, 1,048 gates). This is the instance the tests check against the exact unitary. |

Run from the repository root:

```
python circuits/ch05_bt_link_step/build.py        # writes the QASM files
python circuits/ch05_bt_link_step/count.py        # writes counts.json, prints the tables
python -m pytest circuits/ch05_bt_link_step/test_circuit.py -q -p no:cacheprovider
python circuits/ch05_bt_link_step/verify_independent.py   # the exported files, block by block, and a recount
```

The tests take 35 to 50 s on one machine with 4 threads. The largest statevector has 22 qubits.

Other sizes: `build.lattice_step(L, theta_e, theta_b, steps, order, n_workspace, uncompute)` for any even L, and `build.trotter_circuit(n_links, plaquettes, ...)` for any list of plaquettes.

## 2. How the circuit maps to the chapter

The chapter's Eq. (Hks) without fermions is H = (g^2/2a) sum_links E^2 - (1/g^2 a) sum_plaq (Tr U_p + h.c.). One step of the circuit is

  exp(-i theta_e sum_links E^2) exp(+i theta_b sum_plaq Re Tr U_p),  theta_e = g^2 dt / 2a,  theta_b = 2 dt / (g^2 a).

- Register: g = (-1)^m i^n j^o l^(p+2q) stored as |q p o n m>, the register of arXiv:2208.12309 Eq. (pres). A test checks the paper's worked example (the element stored as |23>).
- Lattice: link (x, y, mu) is register 2(x + L y) + mu. The plaquette at site s is U_x(s) U_y(s+x) U_x(s+y)^-1 U_y(s)^-1.
- Model instance: L = 4, d = 2, a_s = 0.2 fm, first result one step of 0.045 fm/c, full benchmark six steps of 0.025 fm/c, 16 workspace qubits. `count.py` asserts that the generator defaults equal these model values (`Assumptions.L_2028`, `dt_first_2028_fm`, `dt_over_at_2028`, `first_result_steps_2028`, `trotter_steps_2028`, `workspace_2028`).
- A model "step" has one electric layer and one magnetic layer (2 U_F and 1/2 U_Tr per link). The circuit's default `order=1` is that step.

Gate set: Clifford (x, z, s, sdg, h, cx, cz, swap), t and tdg, ccx, and rz with an angle that is not a multiple of pi/4. Rotations by multiples of pi/4 are written as Clifford or T gates. Every arbitrary rotation stays an explicit rz with its angle, so the circuit can be resynthesized. The walker raises an error on any other gate.

## 3. The blocks

All blocks are derived here. Write g = z l^c with z in Q8 and c in Z3. Conjugation by l acts on the exponents of z as (n, o) -> (o, n+o) with no sign, which is two CNOTs.

| block | action | construction | cost per call |
|---|---|---|---|
| `u_inv` | g -> g^-1 | sign of z^-1 (1 Toffoli), controlled power of the conjugation (3 Toffoli), swap of p and q | 4 Toffoli |
| `u_mul` | (g, h) -> (g, gh) | gh = z_g phi^c_g(z_h) l^(c_g + c_h): controlled conjugation (3 Toffoli), Q8 product (2 Toffoli), controlled mod-3 shift, which is affine on (p, q) (3 Toffoli) | 8 Toffoli, no ancilla |
| `u_tr` | phase exp(i theta Re Tr g) | Z-string (Walsh) expansion of Re Tr g with the 8 unused states left free. An exhaustive search over supports finds 8 strings as the minimum | 8 rotations |
| `u_fourier` | sends each isotypic component of the regular representation onto a set of basis states | H on m, a Clifford Bell-basis change on (n, o), one doubly controlled order-3 Clifford (exact in Clifford+T), H on n and o, a Z3 Fourier transform on (p, q) up to output phases | 2 Toffoli + 10 T + 2 rotations, 1 ancilla |
| `u_phase` | phase exp(-i theta_e f(irrep)) in that basis | one flag (2 Toffoli to compute, 2 to uncompute) and 5 Z-string rotations | 4 Toffoli + 5 rotations, 2 ancillas |
| `electric_link` | exp(-i theta_e E^2) | `u_fourier`, `u_phase`, `u_fourier` inverse | 8 Toffoli + 20 T + 9 rotations |
| `plaquette` | exp(i theta_b Re Tr U_p) | product built in place on one link, `u_tr`, uncompute. Gate sequence of arXiv:1903.08807 Fig. 1 | 6 `u_mul` + 6 `u_inv` + 1 `u_tr` = 72 Toffoli + 8 rotations |

The electric term depends on the irrep only. So the circuit does not need the group Fourier transform with a fixed basis inside every irrep block. It needs a unitary that sorts the irreps. `u_fourier` is that unitary and nothing more. It is not the transform of arXiv:2408.00075 and cannot replace it where matrix elements inside an irrep matter.

## 4. What was verified

Reference objects are built with numpy from the definitions: the 24 elements as 2x2 matrices, their products, inverses and traces, and the link electric term as the Cayley-graph Laplacian sum over gamma of (1 - L_gamma) built from matrix products and exponentiated with `scipy.linalg.expm`. None of them uses the gate constructions. The 2x2 matrices do come from the class `BinaryTetrahedral` in `build.py`, and the angles of `u_tr` and `u_phase` are fitted to tables from the same class. The last three rows of the table remove that shared piece. Tolerance: 1e-9 on every amplitude after removing one global phase. Observed errors are at or below 1e-14.

| check | size | reference | observed max error |
|---|---|---|---|
| `u_inv` on all 24 elements | 5 qubits, full unitary | matrix inverse | exact |
| `u_mul` and its inverse on all 576 pairs | 10 qubits, full unitary | matrix product | exact |
| `u_tr` at two angles | 5 qubits, full unitary | exp(i theta Re Tr g) | 7e-16 |
| `electric_link` at three angles | 5 + 2 ancilla, full unitary | expm of the Laplacian. The 24x24 block is unitary, so nothing leaves the register and both ancillas return to zero | 6e-15 |
| `u_fourier` | 5 + 1 ancilla | the Laplacian becomes diagonal | below 1e-9 |
| electric eigenvalues | 24x24 | spectrum of the Laplacian with multiplicities 1, 4, 9, 8, 2 | exact |
| `plaquette`, both uncompute variants | 20 qubits, statevector on a random state of the gauge register | exp(i theta_b Re Tr(g1 g2 g3^-1 g4^-1)) on all 24^4 configurations | 4e-18 (amplitudes of order 2e-3) |
| one-plaquette Trotter step, order 1 and order 2, and two steps | 22 qubits (20 + 2 ancilla), statevector | expm on each link, then the plaquette phase | 4e-17, norm outside the subspace 5e-16 |
| magnetic layer of the full 4x4 benchmark lattice, both variants, and of 2x2 | 160 qubits (40 for 2x2), bit-level simulation of 6 random group configurations each | sum over sites of Re Tr of the plaquette from the lattice definition. All registers restored | below 1e-9 |
| benchmark structure | 176 qubits | each link has one electric block, each of the 16 plaquettes appears once, each link sits in two plaquettes | exact |
| QASM | single plaquette | reload with `qiskit.qasm3.loads`, same gate names in order and same counts | exact |
| the 176-qubit benchmark circuit, one step and six steps, block by block (`verify_independent.py`) | a fresh QASM export is cut into contiguous blocks that exhaust the gate list: 48 blocks for one step, 288 for six. Each electric block is run as a dense unitary on 5 + 2 qubits. Each plaquette block is run on all 24^4 configurations of its four links | quaternion group table. exp(-i theta_e Laplacian) by eigendecomposition, exp(i theta_b Re Tr) of the torus plaquette, with theta_e and theta_b computed in the checker from H_KS. Every link has one electric block per step and every plaquette of the torus appears once per step | electric 6e-15, norm outside the register 6e-16, plaquette 7e-16 |
| single-plaquette step, same checker | 22 qubits, numpy statevector from the QASM text | same references | 2e-16 |
| files on disk | the three QASM files | text equal to a fresh export of the generator | exact |

Tests that fail for a wrong circuit: the opposite sign of theta_b (error 1.2e-2 against amplitudes of 2e-3), a plaquette with one inversion left out, `u_mul` with one Toffoli removed, and an electric phase with the three-dimensional irrep's eigenvalue set to zero. Each is asserted to miss the reference.

The block-by-block rows were added by the independent review. Before them, a generator fault that changed the electric angle only on links 4 and up of the lattice passed all 26 tests, because the lattice electric layer was checked by block count only. With them it fails. Since every block is checked on the whole physical domain of the qubits it touches, with ancillas in zero before and after, and the blocks are contiguous and exhaust the gate list, the 176-qubit circuit is the first-order product of the block operators on the physical subspace, up to a global phase. The same review ran eleven single-fault mutations of `build.py` (a flipped sign, a dropped Z-string, a changed angle, a wrong Toffoli control, a dropped CNOT, a wrong link in a plaquette, a skipped block, a dropped uncompute). Each fails at least one test.

Result of the last run: 29 passed.

## 5. What was not verified

- The full 4x4 step was not run as one 176-qubit statevector. It is checked as a product of blocks, each exact on its whole physical domain (section 4). That check covers the default circuits (order 1, Fig. 1 uncompute, 16 workspace). The other variants of section 6 on the lattice rest on the single-block tests, and for `uncompute="dagger"` on 6 random basis configurations of the magnetic layer.
- No periodic lattice was checked as a statevector. The smallest, 2x2, has 40 qubits.
- Rotation synthesis. The rz gates stay symbolic and c_T is the model formula `common.t_per_rotation`. No synthesizer (gridsynth) is installed.
- The T-depth is a schedule on the circuit's own qubits with no limit on magic states and no routing. It is not a compiled lattice-surgery depth.
- The Trotter error, the value of g^2, and the physical choice of the electric term (section 7, items 3 and 7).
- The Hadamard test, the initial state and the observable A_k. They are not part of the link-step and are not built.
- The published circuits of arXiv:2208.12309 and arXiv:2408.00075 were not rebuilt. The model's per-primitive numbers were compared with those papers' tables, which they match, not with a rerun of their figures.
- Minimality. The counts are upper bounds for these constructions.
- The action on the 8 unused states of each link register (p = q = 1) is not specified. The tests check that states of the gauge register stay in it.

## 6. Counts

From `python circuits/ch05_bt_link_step/count.py` (`counts.json`). Convention from `estimates/common.py`: 7 T per Toffoli (`toffoli_t`), rotations at c_T = `t_per_rotation(eps_rot_for(N_rot))` with eps_syn = 1e-2 per shot. The chapter model states no departure. T-depth: the model's convention (comment above `SCHED_2028`): a Toffoli is one layer, a T gate one layer, a rotation c_T layers. Model values are from `evolution_2028(a, steps)`, `step_depth_2028(a, t_rot)`, `model(a, "2028").intermediates` and `groups.GROUPS["2T"].primitives`. The model reports T counts, not Toffoli counts. The model N_Toffoli in the tables is derived in `count.py` as the model's calls per link-step times its Toffolis per primitive, 172 per link-step.

Primitives, per call, and calls per link-step:

| primitive | model Toffoli | model rotations | model ancilla | model calls | circuit Toffoli | circuit rotations | circuit T | circuit ancilla | circuit calls |
|---|---|---|---|---|---|---|---|---|---|
| U_mul | 22 | 0 | 1 | 6 | 8 | 0 | 0 | 0 | 3 |
| U_inv | 4 | 0 | 0 | 3 | 4 | 0 | 0 | 0 | 3 |
| U_Tr | 0 | 11 | 0 | 0.5 | 0 | 8 | 0 | 0 | 0.5 |
| U_F (model: U_FFT) | 14 | 42 | 2 | 2 | 2 | 2 | 10 | 1 | 2 |
| U_phase | not in the model | | | | 4 | 5 | 0 | 2 | 1 |

First result, one step on 32 links:

| item | circuit | model | circuit / model |
|---|---|---|---|
| logical qubits | 176 (160 + 16 workspace) | 177-178 (160 + 16 + 1-2 Hadamard-test ancillas) | 0.99 |
| N_Toffoli | 1,408 | 5,504 | 0.26 |
| N_rot | 416 | 2,864 | 0.15 |
| N_T_fixed | 640 | 0 | |
| eps_rot | 4.90e-3 | 1.87e-3 | |
| c_T | 18.02 | 19.62 | 0.92 |
| T total | 1.799e4 | 9.473e4 | 0.19 |
| T total at the model's c_T | 1.866e4 | 9.473e4 | 0.20 |
| T-depth | 936 | 1,626 (low) to 8,721 (high) | 0.58 to 0.11 |

Full benchmark, six steps:

| item | circuit | model | circuit / model |
|---|---|---|---|
| N_Toffoli | 8,448 | 33,024 | 0.26 |
| N_rot | 2,496 | 17,184 | 0.15 |
| N_T_fixed | 3,840 | 0 | |
| eps_rot | 2.00e-3 | 7.63e-4 | |
| c_T | 19.51 | 21.11 | 0.92 |
| T per link-step, magnetic | 330.0 | 1,124.1 | 0.29 |
| T per link-step, electric | 251.6 | 1,969.2 | 0.13 |
| T per link-step | 581.6 | 3,093.3 | 0.19 |
| T total | 1.117e5 | 5.939e5 | 0.19 |
| T-depth | 5,963 | 10,407 (low) to 55,714 (high) | 0.57 to 0.11 |

Per block, full benchmark (T at c_T = 19.51):

| block | calls | Toffoli | rotations | T gates | T total |
|---|---|---|---|---|---|
| U_mul | 576 | 4,608 | 0 | 0 | 32,256 |
| U_inv | 576 | 2,304 | 0 | 0 | 16,128 |
| U_Tr | 96 | 0 | 768 | 0 | 14,983 |
| U_F | 384 | 768 | 768 | 3,840 | 24,199 |
| U_phase | 192 | 768 | 960 | 0 | 24,105 |

Variants of the first result (same walker):

| variant | logical qubits | Toffoli | rotations | T gates | T total | T-depth |
|---|---|---|---|---|---|---|
| default (order 1, Fig. 1 uncompute, 16 workspace) | 176 | 1,408 | 416 | 640 | 1.799e4 | 936 |
| uncompute with U_mul run backwards (4 U_inv per plaquette) | 176 | 1,280 | 416 | 640 | 1.710e4 | 932 |
| second-order step (two magnetic half layers) | 176 | 2,560 | 544 | 640 | 2.849e4 | 1,301 |
| 2 workspace qubits | 162 | 1,408 | 416 | 640 | 1.799e4 | 5,028 |
| 64 workspace qubits | 224 | 1,408 | 416 | 640 | 1.799e4 | 498 |

With 16 workspace qubits the one-step circuit has N_T / D_T = 19, above the ten-factory baseline.

### Where the circuit and the model differ by more than 20 percent

Every T line item differs by more than 20 percent. Which side I believe, item by item:

1. **U_mul calls per link-step, 3 against 6.** I believe the circuit. The model takes 6(d-1) from the multiplicity table that the 2T, 2O and Sigma papers reprint from arXiv:2203.02823. That paper's own plaquette figure has six multiplications per plaquette, and so do Fig. 1 and Table I (6 N_P) of arXiv:1903.08807. With (d-1)/2 plaquettes per link, six per plaquette is 3(d-1) per link. The same table has 3(d-1) for U_inv and (d-1)/2 for U_Tr, which do agree with six inversions and one trace per plaquette. The U_mul row is twice the figure in both of its columns (the rectangle figure has 10 multiplications, the table 20(d-1)). I found no construction in those papers that needs twelve. The independent review read the gate labels of both figures of arXiv:2203.02823 from the PDF text: the plaquette figure has 6 U_x, 6 U_-1 and 1 U_Tr, the rectangle figure 10 U_x and 8 U_-1. The published closed forms do use the table value. For 2O, the 2863(d-1) of arXiv:2312.10285 is 6 x 392 + 3 x 112 + 350/2. If this reading is right it also applies to every model that uses `PRIMCOST`. In the Ch. 5 model it alone moves the magnetic term from 1,124 to 662 T per link-step.
2. **U_mul per call, 8 Toffoli against 22.** I believe both numbers as descriptions of different circuits. The model's 154 T is the table value of arXiv:2208.12309. As I read its figure, the published circuit conjugates the control's Q8 part by the target's Z3 exponent, which needs four triply controlled gates. Writing the product as z_g phi^c_g(z_h) l^(c_g + c_h) puts the control on the Z3 exponent of the control register and needs only Toffolis. The 8-Toffoli circuit is checked on all 576 pairs.
3. **U_Tr, 8 rotations against 11.** I believe the circuit. It is a property of the function Re Tr g on this register, found by exhaustive search and checked on the full unitary. The formula for H_Tr as typeset in arXiv:2208.12309 did not reproduce Re Tr g on the register of Eq. (pres) in my check, so I did not use it. I did not look for a relabeling that repairs it.
4. **Electric term, 252 against 1,969 T per link-step.** I believe the circuit for this Hamiltonian. The model counts two full fast Fourier transforms of arXiv:2408.00075 at 14 Toffoli and 42 rotations each. The 42 is consistent with three Z3 transforms at 14 transpiler rotations each. That is my inference from the paper's text, I did not transcribe its figure. The circuit here uses one Z3 transform per `u_fourier` with 2 rotations and 4 T, and replaces the rest by exact Clifford+T gates, because only the irrep label is needed. The count would not carry over to a term that needs the full transform.
5. **U_phase.** The model has no electric phase gate for 2T (`groups.GROUPS["2T"]` has no `U_phi`). The circuit needs one per link-step: 4 Toffoli and 5 rotations, 126 T at c_T = 19.5. Here the model is low, by 4 percent of its own link-step. The multiplicity tables of arXiv:1903.08807 and arXiv:2203.02823 both list one U_phase per link-step.
6. **Rotation count and tolerance.** With 13 rotations per link-step in place of 89.5, eps_rot is 2.6 times larger and c_T drops from 21.1 to 19.5. This is a consequence, not a separate difference.
7. **T-depth.** The circuit's 936 layers for one step lie below the model's low schedule (1,626). The model's schedule uses the published Toffoli and rotation counts, so this follows from items 1 to 4.

Taken together: the model's 2028 line items quote the published tables correctly as far as I could check, with the U_mul multiplicity as the one place where the published table and the published figures disagree. An explicit circuit built for this Hamiltonian is about five times smaller. The six-step full benchmark is then 1.12 times the 1e5 T reference count, not 5.9 times. The circuit was not adjusted toward the model at any point.

## 7. Choices the circuit forced

1. **Plaquette circuit.** The model gives multiplicities only. The circuit builds g3^-1 g4^-1 g1 g2 in place on link 2 and follows arXiv:1903.08807 Fig. 1 (forward multiplications only, 6 U_inv). `uncompute="dagger"` runs U_mul backwards and saves two inversions per plaquette.
2. **U_mul multiplicity.** Six per plaquette, see item 1 above.
3. **Electric term for a discrete group.** The model does not fix it for 2T. Following the chapter text, E^2 is the Cayley-graph Laplacian. The generating set is taken as the 8 elements with Re Tr = 1, the nearest neighbors of the identity, and the result is scaled so that the faithful irrep carries the SU(2) Casimir 3/4. Eigenvalues: 0 (trivial), 0.75 (2), 1.5 (3), 1.875 (the two complex doublets), 2.25 (the two complex singlets). `electric_link` accepts any other eigenvalues that are equal on conjugate irreps, at the same gate count.
4. **Fourier transform.** An irrep-sorting unitary, not the full transform (section 3).
5. **Trotter order.** The model counts a step as one electric and one magnetic layer and calls the steps second order. For six steps the merged second-order sequence has seven magnetic layers. For the one-step first result the second-order step built here has two magnetic half layers, 2.85e4 T (`order=2`). Splitting the electric layer instead, two electric half layers around one magnetic layer, has fewer Toffolis and was not built. The default follows the model's count.
6. **Workspace.** Each electric block uses 2 ancillas. With the model's 16 workspace qubits, 8 links run at once. The circuit holds no Hadamard-test ancilla, so it has 176 qubits where the model has 177-178. The magnetic term uses no ancilla.
7. **Angles.** The chapter does not fix g^2 for the 2028 benchmark and the model reads no anisotropy for the couplings. The exported QASM uses g^2 = 1 as a stand-in (`G2_PLACEHOLDER`) with a = a_s. Counts do not depend on the angles unless one is a multiple of pi/4.
8. **Uncomputation.** All ancilla uncomputation is done with explicit Toffolis. Measurement-based uncomputation would remove 4 of the 8 Toffolis of `electric_link`.
9. **Diagonal phases.** U_Tr and U_phase are Z-string rotations. A flag-based U_Tr (Toffolis in place of rotations) was not built.

## 8. Caveats

- gridsynth, Qualtran and pyLIQTR are not installed and were not used. Installing gridsynth would let the rz gates be synthesized and c_T checked against the formula.
- pyzx and stim are installed but were not used. No automated T-count optimization was run on the circuit.
- All counts are for one machine, serial, one shot.
