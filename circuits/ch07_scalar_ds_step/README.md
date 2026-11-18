# Ch. 7, 2028 benchmark: one Stormer-Verlet step of a free scalar on de Sitter, as an explicit circuit

This directory holds a gate-level logical circuit for the evolution part of the Chapter 7 2028 benchmark
(`applications/app04_curved_space.tex` in the paper, `estimates/ch07_curved_space.py` here). The instance is a
free real scalar on a periodic $4^3$ lattice, the field amplitude on a $K=8$ grid per site (3 qubits per site,
192 qubits), one symmetric second-order step of $\Delta t = 0.5/H_{\rm inf}$ with the scale factor at the step
midpoint, and the per-step grid rescale and shear of the chapter. The circuit is built by code from the instance
parameters, checked against dense matrices and state vectors on reduced lattices, and counted by walking the
gate list. The counts sit next to the model's line items below.

The circuit is the step only. The in-vacuum preparation (30,000 T, stated and not derived in the model) and the
three readout settings are not built.

## Files

| file | content |
|---|---|
| `build.py` | `Instance`, `step_parameters` (the classical coefficients), the block builders, `build_step`, and the optional `merge_rz` pass |
| `test_circuit.py` | 36 pytest checks, references listed under "What was verified" |
| `count.py` | walks the circuit, writes `counts.json` and the OpenQASM 3 file |
| `counts.json` | per-block and total counts for five variants, the model line items with their source, ratios, block index ranges |
| `independent_check.py` | a second check by the reviewer that shares no code with `test_circuit.py` or `count.py` (see below) |
| `ch07_scalar_ds_step_4x4x4_K8.qasm` | OpenQASM 3 of the full benchmark step, 192 qubits, 13,056 instructions, 250 kB |

Other sizes come from the same generator, for example
`build.build_step(build.Instance(dims=(5, 5, 5), n_q=4, m=1.2, n_steps=2))`.
`dims` can have any number of dimensions. With `n_steps > 1` the kinetic half-steps of neighbouring steps merge,
as in the model.

## The step and how it maps to the chapter

In units $H_{\rm inf}=1$, with $a=e^{t}$ and the lattice Laplacian $\Delta$ at spacing $b$,

$$H(t)=\sum_x\Bigl[\frac{\pi_x^2}{2a^3}+\frac{a^3m^2}{2}\phi_x^2\Bigr]+\frac{a}{2}\,\phi^{T}\Delta\,\phi .$$

One step is $K(h/2)\,D(h)\,K(h/2)$ with the coefficients at $a_m=e^{h/2}$, where $K$ is the $\pi^2$ block and $D$
the diagonal block (`_ds2028_blocks` in the model). The register of each site holds $z'=Sz$, $z=(\phi,\pi)$,
with $S$ lower triangular (a dilation and a shear) refit to the continuum site covariance at each kinetic-block
boundary (`ds2028_grid_check`, `frame='adapted'`). In that frame a kinetic block is an $x^2$ phase, a $p^2$
phase and an $x^2$ phase, and $D$ is an $x$ phase with $Q/s^2$. On the grid $x=-(dx/2)\sum_q 2^q Z_q$ with
$dx=\sqrt{2\pi/K}$, so every $x^2$ or $x_xx_y$ phase is a set of $ZZ$ rotations, and $p=F^\dagger xF$ with the
centered DFT $F$.

The blocks of the circuit, in the order they act:

| block label | what it is | gates per site or link | chapter / model item |
|---|---|---|---|
| `shear_in` | $x^2$ phase of the first frame change | 3 $ZZ$ per site | "rescale and shear", not a model line item |
| `centering_in` | $R_z$ layer of the centered DFT | 1 Rz, 1 T, 1 Clifford per site | `centering_layer` |
| `k1_qft_inv` | swap-free inverse QFT | 2 controlled-S, 1 controlled-$R_3$ per site | `qft_pair_rot`, `qft_pair_cs` |
| `k1_pi2` | $p^2$ phase, bit-reversed weights | 3 $ZZ$ per site | `pi2_block` |
| `k1_qft` | swap-free QFT | as `k1_qft_inv` | `qft_pair_rot`, `qft_pair_cs` |
| `d1_site` | site part of $D$ plus the two inner shears | 3 $ZZ$ per site | `diag_block`, $V\,C(n_q,2)$ |
| `d1_link` | $\phi_x\phi_y$ on the 192 links | 9 $ZZ$ per link | `diag_block`, links $\times\,n_q^2$ |
| `k2_qft_inv`, `k2_pi2`, `k2_qft` | second kinetic block | as the first | as the first |
| `centering_out` | inverse $R_z$ layer | 1 Rz, 1 T, 1 Clifford per site | `centering_layer` |
| `shear_out` | $x^2$ phase of the last frame change | 3 $ZZ$ per site | not a model line item |

Gate set: H, CX, S, Sdg, T, Tdg and Rz. A $ZZ$ rotation is CX, Rz, CX. A controlled phase is two single-qubit
phases, CX, one phase, CX, so a controlled-S is three T gates and a controlled-$R_3$ is three Rz($\pm\pi/8$).
Phases at multiples of $\pi/4$ are emitted as exact S and T gates. Every Rz in the file is a rotation that needs
synthesis, and its angle is kept. There is no Toffoli and no ancilla.

## Counts

All counts come from `count.py`, which walks the instruction list of the circuit object. T totals use
`estimates.common`: `t_per_rotation(eps_rot_for(N_rot))` with `EPS_SYN` = $10^{-2}$ per shot, and 7 T per
Toffoli (none here). Model values come from `shot_counts(Assumptions(), '2028', 1)`, `register(a, '2028')` and
`model(a, '2028')` in `estimates/ch07_curved_space.py`. One machine, serial.

Headline variant (`full`: the exact one-step map with both end shears, textbook decompositions):

| line item | circuit | model | ratio | model value read from |
|---|---|---|---|---|
| logical qubits used by the circuit | 192 | 192 | 1.000 | `register()['system_qubits']` |
| logical qubits against the model total | 192 | 222 to 242 | 0.86 to 0.79 | `register()['lq']`, system + $N_{\rm anc}$ = 30 to 50 |
| rotations, diagonal block $D$ | 1920 | 1920 | 1.000 | `shot_counts()['rot']['diag']` |
| rotations, $\pi^2$ phases | 384 | 384 | 1.000 | `shot_counts()['rot']['pi2']` |
| rotations, QFT controlled-$R_3$ | 768 | 768 | 1.000 | `shot_counts()['rot']['qft']` |
| rotations, centering layers | 128 | 128 | 1.000 | `shot_counts()['rot']['centering']` |
| rotations, frame shear at the two ends | 384 | 0 | none | not a model line item |
| $N_{\rm rot}$ | 3584 | 3200 | 1.120 | `shot_counts()['n_rot']` |
| T, QFT controlled-S | 1536 | 1536 | 1.000 | `shot_counts()['t_exact']['qft_cs']` |
| T, centering layers | 128 | 128 | 1.000 | `shot_counts()['t_exact']['centering']` |
| $N_{T,\rm fixed}$ | 1664 | 1664 | 1.000 | sum of `shot_counts()['t_exact']` |
| $N_{\rm Toff}$ | 0 | 0 | none | the model has no Toffoli |
| $\epsilon_{\rm rot}$ | $1.670\times10^{-3}$ | $1.768\times10^{-3}$ | 0.945 | `shot_counts()['eps_rot']` |
| T per rotation | 19.809 | 19.715 | 1.005 | `shot_counts()['t_per_rot']` |
| T, one step | 72,661 | 64,753 | 1.122 | `shot_counts()['t_total']` |

Variants, same instance:

| variant | $N_{\rm rot}$ | $N_{T,\rm fixed}$ | T, one step | ratio to model | step + stated 30,000 T prep |
|---|---|---|---|---|---|
| `full` | 3584 | 1664 | 72,661 | 1.122 | 102,661 (above the $10^5$ target) |
| `in_shear_only` | 3392 | 1664 | 68,703 | 1.061 | 98,703 |
| `no_outer_shear` (the model's accounting) | 3200 | 1664 | 64,753 | 1.000 | 94,753 |
| `full_merged_rz` | 3456 | 1344 | 69,701 | 1.076 | 99,701 |
| `no_outer_shear_merged_rz` | 3072 | 1344 | 61,806 | 0.954 | 91,806 |

The last column adds the model's stated preparation (`t_prep_2028`), which is not built here, so it mixes a
circuit count with a stated number. `merge_rz` is a peephole pass that commutes Z phases through CX controls and
merges them. It removes 128 rotations and 320 T, in the QFTs and the centering layers. After the pass the
per-block attribution in `counts.json` is approximate, because a merged phase keeps the label of its last part.

### Depth

Depth is the longest path through the circuit that counts only non-Clifford gates. The model defines its
T-depth as rotation layers times T per rotation (2033 note in the model docstring), and `t_depth_report` uses
that rule: a rotation counts `t_per_rotation`, an exact T counts 1.

| quantity (`full`) | value |
|---|---|
| non-Clifford layers | 59 |
| rotation layers | 42 |
| T-depth, report convention, any number of rotations in flight | 847 |
| largest number of rotations in one layer | 128 |
| rotation layers with at most 50, 30, 10 rotations in flight | 74, 122, 360 |
| T-depth estimate at 50, 30, 10 in flight (layers $\times$ T per rotation) | 1466, 2417, 7131 |

The model's 2028 T-depth, 346 to 1308, is an input (`t_depth_2028`) from an audit of the retired $K=4$
four-step circuit and its source file is not in the repository, so it is not reproduced here. The circuit's 847
lies inside that range. It needs 128 rotations in flight at its widest layer. If each rotation in flight holds
one synthesis ancilla (`RUS_ANCILLA` = 1) and only $N_{\rm anc}$ = 30 to 50 are available, the depth is about
1.5 to $2.4\times10^3$, above the model's upper value. The width-capped numbers come from a greedy schedule in
circuit order and are not optimized. The order of gates inside each commuting layer follows an edge colouring
(6 link colours, 3 rounds per link), which sets the depth and not the counts.

## What was verified

Tolerance $10^{-9}$ everywhere: largest entry-wise deviation for unitaries, 2-norm of the state difference for
state vectors, after removing one global phase. Observed deviations are $10^{-14}$ or smaller. The references of
the dense-matrix and state-vector checks are numpy constructions from the operator definitions ($x$ = grid
values, $F$ from its formula, $p=F^\dagger xF$, `scipy.linalg.expm`) or the model's own functions, and do not use
the circuit code. Three rows compare the circuit with itself and say so in the reference column: the centering
merge, the per-site blocks of the full $4^3$ circuit, and the Rz merge pass. They check structure. The angles
they rely on are checked by the other rows.

| check | size | reference |
|---|---|---|
| site $x^2$ phase | $n_q$ = 2, 3, 4 | dense $\exp(-i\gamma x^2/2)$ |
| centered DFT (centering + swap-free QFT) | $n_q$ = 2, 3, 4 | $F$ from its definition |
| kinetic block | $n_q$ = 2, 3, 4 | `expm(-i u p^2/2)` |
| link block | 2 sites, $n_q$ = 2, 3 | dense $\exp(-i\kappa\,x\otimes x)$ |
| centering layers cancel between kinetic blocks | 3-site ring, 9 qubits | the unmerged circuit |
| step coefficients | rings, $2^3$, $3^2$, $4^3$, one and two steps | $S_{\rm out}\,(KDK)\,S_{\rm in}^{-1}$ from the Hamiltonian, real-space covariance |
| full step, every unitary entry | 3-site ring, $K=8$, 9 qubits | the model's map (`_ds2028_setup`, `_ds2028_blocks`, frame rule) |
| two merged steps | 3-site ring, $K=4$, 6 qubits | same |
| full step, random state | 4-site ring, $K=8$, 12 qubits | same |
| full step, random state | $3\times3$ torus, $K=4$, 18 qubits | reference map built in the test |
| full step, random state | $2^3$ lattice, $K=8$, 24 qubits | reference map built in the test |
| Gaussian in-state through the circuit | 4-site ring, $K=8$, $m$ = 1.0 and 1.4 | shell errors equal `ds2028_grid_check(m, 8, 'adapted', 1, 4, 0.5, merged=True)` to $10^{-8}$, largest 0.14% |
| full $4^3$ circuit, per-site blocks | 192 qubits | gate-for-gate equal to the single-site blocks above, on all 64 sites |
| full $4^3$ circuit, diagonal block | 192 qubits, 1920 rotations | the $ZZ$ angles decode to $x^{T}Qx/2$ of the reference map |
| Rz merge pass | 3-site ring | same unitary, global phase included |
| OpenQASM 3 export | 192 qubits | loads back with the same counts and angles |

Tests that fail when a block is wrong: a flipped sign of the gradient coupling (deviation 1.7), one dropped
rotation in each of six blocks (deviations 0.017 to 1.1), a $p^2$ phase without bit reversal, the QFT in the
wrong direction, a kinetic block without centering, and the step with one end shear removed.

### Independent check

`independent_check.py` was written in review and repeats the main checks by a different route. It reads the
gate list from the `.qasm` text and counts it with its own counter. It defines $p$ as the generator of the grid
shift, $P=i\log T/dx$ with $T|j\rangle=|j+1\rangle$ and an antiperiodic wrap, and not through the DFT formula.
It solves the frame coefficients from the Hamiltonian with its own mode sums.

| check | result |
|---|---|
| counts from the `.qasm` text | 192 qubits, all acted on, no gate wider than two qubits, 0 Toffoli, 3584 rotations, 1664 T, 72,661 T, depths 59, 42 and 847 |
| kinetic block against `expm` with the shift-generator $P$ | $3\times10^{-15}$ |
| full step, 3-site ring, $K=8$, product of dense `expm` of many-site operators | $9\times10^{-15}$ |
| full 192-qubit file as (site factors) $\cdot$ (link $ZZ$ rotations) $\cdot$ (site factors), each $8\times8$ site factor against a dense reference on all 64 sites | $7\times10^{-15}$ |
| the 1728 link $ZZ$ angles of the 192-qubit file against the Hamiltonian, own neighbour list | $3\times10^{-16}$ |
| Gaussian in-state through the circuit, 4-site ring, $K=8$, shell errors against the continuum map | 0.14% at $m=1.0$, 0.10% at $m=1.4$ |

The factorization is exact because every gate between the first and the last inter-site gate of a site is part
of an inter-site $ZZ$ rotation, which the script asserts. With the two site factors and the link angles each
equal to their reference, the 192-qubit unitary equals the reference map with the $4^3$ coefficients.

## What was not verified

- The 192-qubit unitary as a whole was not simulated. The full circuit is checked block by block as listed
  above, and as an exact factorization in `independent_check.py`. A 3D lattice with six distinct neighbours
  per site was not run as a state vector. The $2^3$ lattice has three distinct neighbours with doubled
  coupling, and the $3\times3$ torus has four.
- The in-vacuum preparation and the three readout settings are not built.
- The rotations are left as Rz. The T total applies the report formula to the rotation count. No synthesized
  Clifford+T sequences were generated and no synthesis error was measured.
- The accuracy of the $K=8$ grid on the $4^3$ lattice. The 0.14% figure is the model's four-site ring check,
  which the circuit reproduces. The site frame on $4^3$ uses the 3D site covariance.
- The model's T-depth input (346 to 1308).

## Choices the circuit forced that the model leaves open

1. **Frame shear at the two ends of the shot.** The chapter says the rescale and shear change only rotation
   angles. That holds for the shears between kinetic blocks, which merge into the site rotations of $D$. The
   first and the last shear of a shot have no neighbouring site rotation to merge into. Each is 192 rotations.
2. **Field mass.** $m=1.0\,H_{\rm inf}$ (the chapter gives 1 to 1.4). It sets angles only.
3. **Lattice spacing and in-state in 3D.** $b$ from $m^2+12/b^2=4$ and Bunch-Davies mode functions at the
   lattice momenta with the $a=1$ ground state for $k=0$, the model's ring definitions carried to $4^3$.
4. **Fourier transform.** $F$ is the centered DFT with the sign of `ds2028_grid_check`. It is built from the
   inverse QFT without SWAP gates, and the $p^2$ phase is applied with bit-reversed weights.
5. **Centering layers.** Merged to one at the start and one at the end of the shot, as the model counts. The
   unmerged form is `merge_centering=False`.
6. **Controlled phases.** The textbook form, three single-qubit phases and two CX each, as the model assumes.
7. **$ZZ$ rotations.** CX, Rz, CX on the two data qubits. No parity ancilla.
8. **Gate order inside commuting layers.** Edge colouring for depth.
9. **Ancilla.** None. The model's $N_{\rm anc}$ = 30 to 50 is an allowance for rotation synthesis, which this
   circuit leaves abstract.
10. **Synthesis tolerance.** $\epsilon_{\rm rot}$ is set from the rotations of the step alone, as in the model.
    A built preparation would add its rotations to $N_{\rm rot}$.

## Where the circuit and the model differ

No total differs by more than 20%. Two line items do.

**End shears: circuit 384 rotations, model 0.** Without them the circuit has exactly the model's counts
(`no_outer_shear`: 3200 rotations, 1664 T, 64,753 T). With them the step is 72,661 T, 12% above the model, and
with the stated preparation the shot is 102,661 T, above the $10^5$ target. I believe the circuit for what a
standalone step costs. For the benchmark shot the two layers differ:

- The last shear is a diagonal phase. It does not change a field-basis measurement (tested). The chapter's other
  two settings change the final $\pi^2$ half-step and, as described, do not need it either. That was read from
  the chapter and not built. Dropping it is safe, so the model is right to leave it out.
- The first shear acts on the prepared state. It can be folded into the preparation only if the Gaussian network
  prepares the state with that extra site phase. The preparation is stated and not built, so this is open. If
  the layer stays in the step, the step is 68,703 T (`in_shear_only`, 6% above the model) and the shot with the
  stated preparation is 98,703 T, inside the target.

**Logical qubits: circuit 192, model 222 to 242.** The circuit uses no ancilla because rotations are kept as Rz.
The model's extra 30 to 50 qubits are not derived from a circuit. The depth table shows what they allow if each
rotation in flight needs one.

Smaller points. The model's note that the busiest qubit sees 36 rotation layers counts 6 rotations per QFT pair
on one qubit. The circuit has 34 rotation strings on its busiest qubit in the model's accounting and a longest
rotation path of 37 layers. The merge pass shows that 128 rotations and 320 T of the textbook form are removable
without ancilla.

## Reproduce

```
cd FNAL_Resource_Estimation
python circuits/ch07_scalar_ds_step/count.py                                   # about 1 s, rewrites counts.json and the .qasm
OMP_NUM_THREADS=4 python -m pytest circuits/ch07_scalar_ds_step/test_circuit.py -q   # 36 tests, about 1 min
OMP_NUM_THREADS=4 python circuits/ch07_scalar_ds_step/independent_check.py           # about 3 s
```

Needs numpy, scipy, qiskit and qiskit-aer. The largest test is a 24-qubit state vector (peak memory 1.6 GB, at most
4 threads). A rotation synthesizer (gridsynth or a repeat-until-success compiler) is not installed here. With
one, the Rz gates of the exported file could be replaced by explicit Clifford+T sequences and the formula-based
T total checked.
