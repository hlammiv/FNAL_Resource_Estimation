# Ch. 2, 2028 benchmark: gate-level circuit of the four-nucleon charge correlator on 3^3

This folder holds an explicit logical-level circuit for the 2028 benchmark of Chapter 2 (neutrino-nucleus response),
a generator for other sizes, tests against an independent numpy reference, and counts obtained by walking the circuit.
The chapter model is `estimates/ch02_nu_nucleus.py`. The physics instance is defined in
`applications/app01_neutrino_nucleus.tex` of the paper ("Concrete instantiations").

Result in one line: the walked counts equal the model's three line items exactly at 2, 4 and 6 steps
(1,192 rotations, 2,240 Toffolis, 744 other T, 36 qubits, 38,948 T at 6 steps). The circuit also gives two things the
model does not have: a T-depth, and the number of qubits once the rotations are synthesized by repeat-until-success.

## Files

| file | content |
|---|---|
| `build.py` | `Instance`, `benchmark_instance(steps)`, `build_circuit(inst)`, and the block emitters. `python build.py` writes the QASM files. |
| `count.py` | walks the circuit, counts per block, compares with the model, computes depth. `python count.py` writes `counts.json`. |
| `counts.json` | all counts, per block and total, the model values and the ratios |
| `test_circuit.py` | pytest suite against the numpy reference. `python test_circuit.py` writes `verification.json`. |
| `verification.json` | the errors observed in the complete-circuit checks |
| `ch02_3x3x3_steps6_re_p0p0.qasm3` | OpenQASM 3 of the deepest benchmark circuit (0.27 MB, 15,671 gates) |
| `ch02_3x3x3_steps6_re_p0p0_merged.qasm3` | the same circuit with the ancilla phase gates merged (variant, see below) |
| `independent_check.py` | second check, added in the verification pass: its own QASM parser, recount and depth, an exact simulation of the full 36-qubit file, and its own reference. It imports neither `count.py` nor `test_circuit.py`. |
| `independent_check.json` | what that check found for the 6-step circuits (`independent_check_steps2.json` and `independent_check_steps4.json` for the shorter ones) |

```
python build.py                 # QASM 3 exports
python count.py                 # counts.json and the comparison table
python -m pytest test_circuit.py -q    # 39 tests, 70 to 90 s on one machine with 4 threads
python test_circuit.py          # verification.json
python independent_check.py --count    # recount of the QASM text, seconds
python independent_check.py --pairs    # full 36-qubit check of all four current pairs, 6 to 7 min, 1.5 GB
```

Other instances come from the same generator, for example

```python
from build import benchmark_instance, build_circuit
from dataclasses import replace
built = build_circuit(benchmark_instance(steps=4, current_pair=(0, 1), component="im"))
small = build_circuit(replace(benchmark_instance(6), dim=2, n_particles=3))   # 3 nucleons on 3^2, 20 qubits
```

`Instance` takes the dimension, the number of nucleons, the couplings, the preparation layers, the time, the number of
steps, the momentum, the reference energy, the current pair, and the component. Only L = 3 sites per axis is
implemented, because the kinetic block uses the three-site identity below.

## What the circuit is

One Hadamard test of

    C_ab(tau) = <psi| U_a(q)^dagger exp(-i (H - E0) tau) U_b(q) |psi>,   U_f(q) = exp(i q x_f),

for one nucleon of each of four species on a periodic 3^3 lattice with

    H = sum_f [6 t - t (hops to the two other sites of each axis)] + C0 sum_pairs delta + D0 sum_triples delta.

The chapter's correlator is the average of the four C_ab over the two protons (species 0 and 1), real and imaginary
parts, so one estimate uses eight circuits that differ only in which coordinate carries the current and in one S gate.
All eight have the same counts (`counts.json`, `all_eight_tests_same_count`).

Gate set: h, s, sdg, x, cx, t, tdg, ccx, rz, and one measurement of the ancilla. Every arbitrary rotation is an rz
with its angle. None of the 1,192 rz angles is a multiple of pi/4 (tested).

Register: 24 position qubits (two per three-valued coordinate, x = bit0 + 2 bit1), 11 clean workspace qubits, one
Hadamard ancilla. All 36 are used.

| chapter statement | block in `build.py` | one instance: rz / Toffoli / T |
|---|---|---|
| product of four uniform one-particle states | `prepare3` on 12 coordinates | 12 / 0 / 24 |
| contact terms as phases through reversible position comparisons | `potential_layer` (6 pairs, 4 triples) | 10 / 148 / 0 |
| the same, controlled on the Hadamard ancilla | `potential_layer(control)` | 30 / 148 / 0 |
| kinetic term in its eigenbasis: basis change, phase, inverse | `kinetic_layer` (12 coordinates) | 60 / 0 / 48 |
| the same, controlled | `kinetic_layer(control)` | 61 / 24 / 48 |
| charge insertion exp(i q x_p) | `current` | 6 / 0 / 0 |
| reference-energy phase | one rz on the ancilla | 1 / 0 / 0 |

The kinetic block uses K_coord = t I - 3t |s><s| on three sites, with |s> the uniform state. `prepare3` maps |00> to
|s> with one rz and a controlled Hadamard (two T). A pair comparison uses 6 comparison bits and a chain of 5
Toffolis, a triple comparison 12 bits and 11 Toffolis, each chain computed and uncomputed. The triple chain sets the
11 workspace qubits. The sequence is eight preparation layers exp(-i dt_V V) exp(-i dt_K K), then the controlled
second-order evolution K/2, (V, K) x (n-1), V, K/2.

## Count comparison

Counts are from `count.py`, which walks `circuit.data`. T totals use `estimates/common.py`: `toffoli_t` (7 T) and
`t_per_rotation(eps_rot_for(N_rot))` with eps_synth = 1e-2. Model values are `counts('2028', n)` and
`circuit(Assumptions(), '2028', n)` in `estimates/ch02_nu_nucleus.py`.

| steps | item | circuit | model | circuit / model |
|---|---|---|---|---|
| 6 | logical qubits | 36 | 36 (`model(a,'2028').lq`) | 1.000 |
| 6 | N_rot | 1,192 | 1,192 | 1.000 |
| 6 | N_Toff | 2,240 | 2,240 | 1.000 |
| 6 | N_T fixed | 744 | 744 | 1.000 |
| 6 | T total, c_T = 18.896 | 38,948.3 | 38,948.3 | 1.0000 |
| 4 | N_rot / N_Toff / N_T fixed | 1,010 / 1,896 / 648 | 1,010 / 1,896 / 648 | 1.000 |
| 4 | T total | 32,866.4 | 32,866.4 | 1.0000 |
| 2 | N_rot / N_Toff / N_T fixed | 828 / 1,552 / 552 | 828 / 1,552 / 552 | 1.000 |
| 2 | T total | 26,811.8 | 26,811.8 | 1.0000 |
| per step | N_rot / N_Toff / N_T fixed | 91 / 172 / 48 | 91 / 172 / 48 | 1.000 |

Per block family at 6 steps (T at the whole circuit's c_T):

| family | blocks | N_rot | N_Toff | N_T fixed | Clifford | T | share |
|---|---|---|---|---|---|---|---|
| uniform initial state | 1 | 12 | 0 | 24 | 180 | 250.8 | 0.6% |
| preparation, potential | 8 | 80 | 1,184 | 0 | 2,688 | 9,799.7 | 25.2% |
| preparation, kinetic | 8 | 480 | 0 | 384 | 3,456 | 9,454.2 | 24.3% |
| current in and out | 2 | 12 | 0 | 0 | 8 | 226.8 | 0.6% |
| controlled kinetic | 7 | 427 | 168 | 336 | 3,024 | 9,580.7 | 24.6% |
| controlled potential | 6 | 180 | 888 | 0 | 2,136 | 9,617.3 | 24.7% |
| energy phase, H, measurement | 3 | 1 | 0 | 0 | 2 | 18.9 | 0.0% |
| total | | 1,192 | 2,240 | 744 | 11,494 | 38,948.3 | |

The study's own recount file (`neutrino_feasibility/larger_circuit_counts_steps6.json`) lists the same numbers after
the initial state (12 / 0 / 24), after preparation (572 / 1,184 / 408) and for the complete circuit.

No line item differs from the model. That is expected. The model's literals were tallied from the study's
`larger_circuit.py`, and this circuit is the same construction written out gate by gate. The agreement confirms that
the tally rules (three rotations per controlled phase, one T per Ry(pi/4)) are what an explicit decomposition gives.
It is not an independent derivation of the construction. The circuit was not adjusted to reach it.

## T-depth

The chapter model defines no T-depth (`depth_status = serial_t_count_bound`). `count.py` computes the longest path
through the circuit counting non-Clifford gates only, under three weightings.

| convention | depth |
|---|---|
| non-Clifford layers: each T, Toffoli or arbitrary rz is one layer | 2,657 |
| T layers: rz = c_T sequential T, Toffoli = Qiskit's 7-T circuit (T-depth 4) | 15,757 |
| T layers: rz = c_T sequential T, Toffoli = one layer | 10,078 |

The second row is the one comparable to the other chapters, which count a rotation as c_T sequential T
(`estimates/ch04_hybrid_lqcd.py`). With it, `estimates.common.depth_exports` gives F* = N_T / D_T = 2.47, below the
10 factories the 1 us per T convention assumes, and a shortest campaign of 164,424 x 15,757 x 10 us = 25,909 s, or
7.2 h on one machine. The model's serial bound is 6,420 s, or 1.8 h. With the Toffoli as one layer the floor is
16,570 s, or 4.6 h. The model file calls its wall times lower bounds (module docstring). The paper says that the wall
times are serial T-counts with no T-depth derived, and it calls the 2033 ones lower bounds. It does not say this of the
2028 figure of 1.8 h. This puts a number on it for this layout: 2.6 to 4.0 times longer.

The depth is that of the exported (p0, p0) circuit. It depends slightly on which proton carries the incoming current:
15,757 for proton 0 and 15,831 for proton 1, with the same counts (`independent_check.json`).

The depth is high because every evolution gate is controlled on the one ancilla and every contact reuses the same
workspace, so little runs side by side. Only the uncontrolled kinetic layers of the preparation are parallel.

## What was verified

The reference is in `test_circuit.py` and never calls the circuit code. It holds the state as a tensor over the
coordinates, builds the 3 x 3 hopping matrix and the contact potential from their definitions, and applies the
chapter's sequence of exponentials. Tolerance: 1e-9 on the largest difference of any amplitude or matrix element.
Observed errors are below 2e-14 (`verification.json`).

Every block type at the full-instance size:

- uniform-state block and its inverse, 2 qubits, dense unitary
- kinetic block of one coordinate against expm of the hopping matrix, uncontrolled (2 qubits) and controlled
  (4 qubits, control off is the identity, workspace returns to zero), three values of dt including a negative one
- a kinetic layer on two coordinates with the exact global phase, uncontrolled and controlled
- current insertion, 3 active qubits, both signs
- pair contact, 18 qubits, uncontrolled and controlled, all 2^13 basis inputs including the unused pattern 11
- triple contact, 30 qubits, uncontrolled and controlled, all 2^19 basis inputs. These blocks contain only x, cx,
  ccx and rz, so they are a permutation times a diagonal and a basis-state propagation checks them completely. The
  propagator is itself checked against a dense unitary on a 6-qubit contact.
- a full potential layer (6 pairs, 4 triples, both couplings) for 4 nucleons on 3^2, every legal basis state

Complete Hadamard-test circuits, all eight preparation layers and six steps with the benchmark's angles, the full
final statevector compared with the reference including the global phase, and the ancilla expectation value compared
with Re or Im of C_ab:

| instance | qubits | gates | max amplitude error | error of <Z> |
|---|---|---|---|---|
| 4 nucleons on 3^1, all four pairs, Re and Im | 12 | 5,310 | 6.7e-15 | 2.2e-14 |
| 2 nucleons on 3^3 | 18 | 4,787 | 5.0e-15 | 4.2e-15 |
| 3 nucleons on 3^2 | 20 | 5,992 | 1.4e-14 | 1.2e-14 |
| 4 nucleons on 3^2, the largest that fits in 24 qubits | 24 | 10,491 | 4.4e-15 | 9.1e-15 |

The 24-qubit run takes about 50 s with 4 threads on one machine.

Tests that fail when something is wrong:

- a reference with the triple term dropped, the pair sign flipped, q flipped, t flipped, 5 steps instead of 6, the
  wrong current pair or a shifted E0 differs from the circuit by more than 1e-3 in 2-norm (the smallest is 2.9e-3,
  for the step count)
- a circuit with one Toffoli removed, one rz sign flipped or one T removed fails the same comparison
- the 30-qubit triple contact with the opposite angle fails the basis-state check

The reference on the full instance (27^4 = 531,441 coordinate states, numpy only) reproduces the chapter's numbers:
E0 = -52.67220 MeV by Lanczos, prepared energy -52.45933 MeV, overlap squared 0.997943, and the study's
product-formula correlators at 2, 4 and 6 steps to 1e-7. The reference's product formula converges at second order
to exp(-i H tau) from a matrix of explicit elements (4 nucleons on 3^1).

The exported QASM 3 file reloads to the same 15,671 gates, qubits and angles as the generator's circuit.

### The full 36-qubit circuit (`independent_check.py`, added in the verification pass)

The exported file is simulated exactly, gate for gate, from its QASM text. No 2^36 vector is needed, for two reasons
that the script checks on the gate list and does not assume. Between its two Hadamards the ancilla is touched only by
rz and by cx with the ancilla as control, so the run splits into an ancilla = 0 and an ancilla = 1 branch. In each
branch every workspace qubit is back in |0> at the end of every group of consecutive gates, so the state kept is that
of the 24 position qubits (2^24 amplitudes). A group is a run without h, evaluated on every basis state of its support,
or a run on at most four position qubits, multiplied out as a dense matrix.

The reference is a second one, separate from `test_circuit.py`: the kinetic term by FFT in momentum space, the ground
state by Lanczos on the nearest-neighbour stencil, the exact correlator by a Chebyshev expansion of exp(-i H tau).
It reads the eight preparation angles from the study file and types everything else from the chapter.

| check on the 6-step circuits | result |
|---|---|
| prepared state against the reference, all 2^24 amplitudes | max error 1.4e-15, weight outside the legal coordinates 3e-15 |
| energy and ground-state overlap of the state the circuit prepares | -52.45933 MeV and 0.997943 |
| ancilla = 0 branch against the prepared state | 2.6e-15 |
| ancilla = 1 branch against U_a^dagger S(tau) U_b psi | 3.1e-15 (the nearest wrong pair is off by 1.2e-2) |
| C_00 from the exported file | 0.653018430 - 0.608951249i, error below 4e-13 |
| merged-variant file | the same C_00, error below 4e-13 |
| average of the four pairs, from four gate-level circuits | 0.511346 - 0.326093i, equal to the study's value to 3e-13 |
| bias against the recomputed exact correlator 0.522251 - 0.316777i | 0.01434, 2.35% of abs(C) |

The other three pairs are written by `build.py` and read back as QASM text. The same check on the 2-step and 4-step
circuits (`--steps 2`, `--steps 4`, all generated) gives 0.637891 - 0.129481i and 0.584156 - 0.240493i with biases
0.00609 and 0.01098, the chapter's table. The reference also recomputes the preparation-only bias 0.00414 and the
static strength 0.65942. The recount in the same script reads the QASM text with its own parser and its own typed-out
7-T Toffoli, and gives the same counts and depths as `count.py`.

A mutation pass on a scratch copy of `build.py` changed one sign, angle, term or ordering at a time (22 cases).
Twenty were caught by the tests against the reference. A wrong tau in `benchmark_instance` was caught only by the
comparison with the stored QASM file, so it passed once the file was regenerated. `test_benchmark_instance_is_the_chapters`
now pins the benchmark's parameters to the chapter's values.

## What was not verified

- The pytest suite itself stops at 24 qubits. The full 36-qubit circuit is covered by `independent_check.py` only,
  which is exact but is not a 2^36 statevector run.
- The exact correlator in the pytest suite is the study's. It is recomputed in `independent_check.py`.
- Rotation synthesis. The 1,192 rz are not compiled to Clifford+T. The T total uses the fit in `common.py`, so it
  is an expected count, not a compiled one. The synthesis error budget of 1e-2 per shot was not simulated.
- Any layout, routing, lattice-surgery schedule, factory count or logical-noise effect.
- The shot counts and the estimator.
- Lattices other than L = 3 per axis. A larger L needs a different kinetic block.

## Choices the circuit forced

1. Ry(theta) is written S H Rz(theta) H Sdg. The controlled Hadamard in `prepare3` is Ry(pi/4) CX Ry(-pi/4), two T.
2. A controlled phase is three rz and two cx. The model's tally assumes this.
3. A phase gate P(theta) is emitted as rz(theta), and the circuit's global phase is tracked so each block equals its
   target exactly.
4. The controlled kinetic block computes the |s> flag into a workspace qubit with a Toffoli pair and puts a
   controlled phase on it. The chapter says only "a basis change, a phase, and the inverse basis change".
5. The identity part of the kinetic term is a separate rz on the ancilla in every controlled kinetic layer, and the
   reference-energy phase is one more. They are not merged in the baseline circuit.
6. Workspace is taken from the lowest index, and the ten contacts of a layer run one after another. Two disjoint
   pair contacts would fit side by side in the 11 workspace qubits. That changes the depth and not the count.
7. Pairs and triples are applied in lexicographic order. They commute.
8. The momentum is along axis 0, the protons are species 0 and 1, and the exported file is the (p0, p0) real part.
9. E0 is read from the study's eigensolve (-52.67220268742125 MeV) and the ancilla phase is (E0 - 24 t) tau, because
   the circuit evolves the hopping part without the constant 6 t per nucleon.
10. The preparation angles are read at full double precision from `unitary_prep_results.json`.
11. Only the ancilla is measured. The study's QASM 2 files measure all 36 qubits.
12. Toffolis stay as ccx and count 7 T each. The weighted depth uses Qiskit's 7-T circuit, which has T-depth 4.
13. eps_rot is set from this circuit's own N_rot: 2.90e-3 at 6 steps, c_T = 18.896.

## Where the circuit and the model differ

Counts: nowhere. The three line items and the T total agree exactly.

Qubits once rotations are synthesized. The model and the circuit both have 36. The report's rotation cost is the
repeat-until-success fit, which holds one clean ancilla per rotation in flight (`common.RUS_ANCILLA = 1`, itemized in
the registers of Chs. 6, 9 and 10). `count.py` finds 104 rotations that are applied while all 11 workspace qubits are
busy, inside the triple contacts: 80 on the workspace flag and 24 on the Hadamard ancilla. The 24 commute out of the
comparison. The 80 do not. So with repeat-until-success synthesis the register is 37, or the 80 flag rotations
need ancilla-free synthesis at a higher T cost. The circuit is the side to believe here. 36 is right for the
unsynthesized circuit, and the chapter should say 37 if it keeps the report's synthesis convention, or state that
those rotations use ancilla-free synthesis.

Wall time. See the T-depth section. The printed 1.8 h assumes ten factories are kept busy. With F* = 2.47 they are
not, and the campaign takes 7.2 h on one machine under the depth convention above. The depth-limited number is the
one to believe for this layout. The model file states that its walls are lower bounds. The paper prints 1.8 h for
2028 without that qualifier, so the chapter should either add it or quote the depth-limited time.

## Variant built here: merged ancilla phases

Between its two Hadamards the ancilla is only ever a control, so every phase gate on it commutes with everything and
the 156 of them merge into one rz (`Instance(merge_ancilla_phases=True)`). The variant is checked against the same
reference on 12 and 20 qubits.

| | N_rot | N_Toff | N_T fixed | c_T | T total | T / model | T-depth (report weighting) |
|---|---|---|---|---|---|---|---|
| baseline | 1,192 | 2,240 | 744 | 18.896 | 38,948.3 | 1.000 | 15,757 |
| merged | 1,037 | 2,240 | 744 | 18.781 | 35,899.6 | 0.922 | 15,635 |

This is 13% fewer rotations and 8% fewer T. It also removes the 24 ancilla rotations from the triple contacts. The
baseline is kept as the comparison circuit because it is what the model counts.

## Other recompilations to try (not built, no numbers claimed)

- Every Toffoli targets a clean workspace qubit and is later undone by the same Toffoli. Measurement-based
  uncomputation (arXiv:1709.06648) applies to all of these pairs.
- The controlled phase on a comparison flag could be one more Toffoli pair with the ancilla and a single rz. The
  triple contact would then need a twelfth workspace qubit.
- The last basis change of the final kinetic layer is outside the ancilla's backward light cone on 11 of the 12
  coordinates and can be dropped. The basis changes at the boundary between preparation and evolution cancel on the
  same 11 coordinates.
- Disjoint pair contacts in parallel, to lower the depth.

## Tools

Used: Python 3, numpy, scipy, qiskit 2.4.1, qiskit-aer 0.17.2. Nothing was installed. gridsynth would let the 1,192
rz be synthesized exactly and the c_T fit be checked against a compiled T count. It is not installed, so that check
was not done.
