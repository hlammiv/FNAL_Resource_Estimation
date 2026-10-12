# Larger 2028 candidate: four nucleons on a 3×3×3 lattice

**Decision, 2026-10-04:** the user accepted this as the 2028 candidate ("ok, lets use it then") and requested a separate 2033 model closer to 1,000 logical qubits. Manuscript integration remains pending discussion of the paired 2033 proposal; see `../2033_model_study/PROPOSAL.md`.

**Recommendation for discussion:** use the compact-coordinate, deterministic-preparation candidate below in preference to the 2×2×2 toy. It improves the spatial model and admits a longer-time correlator, while its estimated complete-circuit cost remains near the 2028 reference budget. This is a feasibility study; the manuscript and its published resource tables have not been changed.

## What improved

| Quantity | First trial | Larger candidate |
|---|---:|---:|
| Spatial sites | 8 | 27 |
| Box side at a=1.4 fm | 2.8 fm | 4.2 fm |
| Smallest nonzero axial momentum | 443 MeV | 295 MeV |
| Algorithmic logical qubits | 13 | 36 |
| Preparation | Generic relative-state loader | Eight optimized kinetic/contact layers |
| Complete short-time circuit, leading T estimate | 77,524 | 54,482 |
| Extended-time complete circuit, leading T estimate | Not costed in first trial | 80,302 |

The larger circuit is cheaper at short times because its preparation is structured, not because more sites are intrinsically cheaper. The 36 qubits comprise 24 position bits, 11 reusable clean workspace bits, and one Hadamard ancilla. A 20-bit fully packed position register is possible in principle, but here each three-valued coordinate uses two bits, simplifying the circuit. This differs deliberately from the earlier rough 20-system-qubit comparison.

## Model and observable

We retain one particle of each spin–isospin species, with the SU(4)-symmetric pair and triple contact interaction from the lattice Hamiltonian and parameter table of [Roggero et al.](https://arxiv.org/abs/1911.06368). At a=1.4 fm the parameters are t=10.5794, C0=−98.2265511, D0=127.839693 MeV. No fit to the tiny box or physical helium has been performed.

The observable is C(τ)=⟨ψ0|J† exp[−i(H−E0)τ] J|ψ0⟩, where J=(exp(iq x_p1)+exp(iq x_p2))/2 and q=2π/(3a). Each of the four terms is a unitary Hadamard test with a known coefficient 1/4. Uniformly sampling the four pairs gives a bounded ±1 estimator without postselection or an assumed block-encoding normalization. Each accepted execution includes a fresh preparation; all executions are accepted in the ideal unitary design.

The physical full fixed-species space has 27^4=531,441 states. Translation symmetry reduces the classical calculation to 27^3=19,683 relative-coordinate states. Ground-state preparation stays in zero total momentum. The current changes the total momentum, and the response calculation explicitly uses the appropriate nonzero-momentum sector; it is not evolved in the ground-state sector by mistake.

This is still a finite-volume model benchmark, not a controlled prediction for helium scattering. The box is small and the axial q a=2π/3 gives substantial lattice dispersion artifacts. The momentum is also high for a controlled pionless interpretation. The ground energy changes from −95.02 MeV in the eight-site trial to −52.67 MeV here, which reinforces the need for a volume/discretization study. These energies are not predictions for the physical helium binding energy. The model can validate the quantum response workflow and resolve the effect of its interactions; it cannot support a claimed reduction in DUNE's uncertainty.

## Preparation actually tested

A classical energy optimization selects eight layers of exp(−i β K) exp(−i γ V), applied to the product of four spatially uniform one-particle states. This yields:

- Exact ground-state energy: **−52.67220 MeV**.
- Prepared-state energy: **−52.45933 MeV**.
- Ground-state overlap squared: **0.9979429**.
- No filters, heralded projections, or amplitude amplification.

Parameters and states are saved in `unitary_prep_results.json` and `unitary_state_8.npy`. Classical optimization and exact validation are intentional for this classically checkable example; no quantum advantage is claimed. A 2/4/6/8-layer sweep attained ground-state fidelities of approximately 76%, 93%, 98.87%, and 99.79%, respectively. These are numerical results for this Hamiltonian, not guarantees for larger nuclei.

A correlated Jastrow trial was also tested. It reached 99.53% fidelity, but the straightforward sequence of local filters succeeded with probability only 1.13×10^−5. That route is not the recommendation. Its data are preserved to distinguish a good variational wavefunction from an efficient preparation circuit.

## Complete circuit and resource accounting

For an L=3 coordinate, the kinetic matrix is t I−3t|s⟩⟨s| with |s⟩=(|0⟩+|1⟩+|2⟩)/√3. Each coordinate kinetic evolution is therefore an explicitly constructed basis change, phase, and inverse basis change. This identity avoids a generic Hamiltonian simulation of all site-to-site hops. Controlled evolution retains the identity phase on the Hadamard ancilla.

Contact phases compute equality of particle positions, apply a phase, and uncompute the workspace. Pair comparisons use six position bits, triple comparisons twelve. The triple comparison sets the peak workspace of eleven bits. All comparator Toffolis, current phases, preparation rotations, and controlled kinetic phases are included.

| Complete circuit | Evolution steps | Arbitrary rotations | Toffolis | Other exact T | Leading total T |
|---|---:|---:|---:|---:|---:|
| τ=0.01 MeV^−1, 1.97 fm/c | 2 | 828 | 1,552 | 552 | 54,482 |
| τ=0.03 MeV^−1, 5.92 fm/c | 6 | 1,192 | 2,240 | 744 | 80,302 |

The counts include eight preparation layers and a single measurement shot. They use seven T per Toffoli, and a coherent rotation-approximation allowance of 0.005 divided across the arbitrary rotations. The report's leading Ross–Selinger formula gives 3 log2(1/ε) per rotation. This is a synthesis estimate, **not an exact compiled T count or a rigorous upper bound**; [Ross–Selinger](https://arxiv.org/abs/1403.2975) includes a subleading term. An illustrative addition of ten T per rotation gives 62,762 and 92,222 T respectively. That check indicates some headroom but is not a bound.

The report's randomized repeat-until-success convention would give lower expected counts, approximately 26,812 and 38,948 T. It is not used to establish the headline fit because its runtime and approximation assumptions differ. Factories, physical qubits, routing, idle errors, and wall time are not inferred from algorithmic width or T count.

Explicit unsynthesized QASM circuits for both depths are saved. Decomposing their controlled phases and Toffolis independently with Qiskit reproduces the above arbitrary-rotation and exact-T counts. Import these QASM files with `qiskit.qasm2.LEGACY_CUSTOM_INSTRUCTIONS` because Qiskit's QASM2 exporter uses its legacy `p`/`cp` names.

## Numerical observable checks

| Time (fm/c) | Exact C | Prepared, product-formula C | Absolute complex bias |
|---|---|---|---:|
| 0 | 0.659422 | 0.658549 | 0.000873 |
| 1.97 | 0.640631 − 0.124040i | 0.637891 − 0.129481i | 0.006093 |
| 3.95 | 0.590847 − 0.231793i | 0.584156 − 0.240493i | 0.010975 |
| 5.92 | 0.522251 − 0.316777i | 0.511346 − 0.326093i | 0.014342 |

The longest-time bias is **2.35% of |C|**, including preparation bias and product-formula error, before gate synthesis, logical noise, and shot noise. At that time preparation alone contributes 0.00414 absolute and the product formula 0.01243; they should not be conflated. Doubling the evolution steps reduces the latter to 0.00309, consistent with second-order convergence, but the extra-depth calculation is a classical validation check, not the budgeted circuit.

The free ground-state static strength is 0.5, compared with 0.65942 for the interacting ground state, a 31.9% enhancement. At the last time the free response is 0.28997−0.40733i. Thus the chosen measurements resolve model interaction effects rather than only an imperceptible circuit difference. The free comparison changes the Hamiltonian and its ground state; it is not an isolated measurement of one correlation term.

A plausible first campaign is three nonzero times (1.97, 3.95, 5.92 fm/c), both real and imaginary parts, with ±0.02 statistical error per component. Hoeffding's bound and a union bound give **27,404 shots per estimate**, or **164,424 shots** for six estimates at at least 95% simultaneous coverage. A conservative bound using the deepest estimated circuit for every shot is 1.32×10^10 T, before an additional static-strength check. This is a pilot accuracy target, not a claim of 5% total uncertainty after all hardware errors. The bounded synthesis allowance can add at most 0.01 absolute bias to a Hadamard expectation; logical-noise contributions still need a hardware model.

## Comparison with occupation encoding

The same physical system uses 108 occupation qubits, or 109 with the Hadamard ancilla. An explicit comparison construction uses the same eight optimized kinetic/contact layers, diagonal onsite Pauli phases, and exact three-mode kinetic basis changes. Its leading complete-circuit estimates are **0.72 million T** at the short time and **1.14 million T** at the longer time. The report's more optimistic RUS convention still gives about 0.24 and 0.37 million T.

These are costs of the specified construction, not lower bounds or a proof against occupation encoding. Jordan–Wigner parity strings contribute Clifford operations in this comparison. An optimized occupation encoding might do better, but the available qubits alone do not make that construction fit. The compact-coordinate circuit currently supplies the stronger candidate.

## Validation and reproduction

1. `python3 editorial_review/neutrino_feasibility/larger_box.py`
2. `python3 editorial_review/neutrino_feasibility/unitary_prep.py`
3. `python3 editorial_review/neutrino_feasibility/larger_response_check.py`
4. `python3 editorial_review/neutrino_feasibility/larger_circuit.py --steps 2 --tau .01`
5. `python3 editorial_review/neutrino_feasibility/larger_circuit.py --steps 6 --tau .03`
6. `python3 editorial_review/neutrino_feasibility/encoding_comparison.py`

Dependencies: NumPy, SciPy, Qiskit, and the repository's existing synthesis conventions. Setting `OPENBLAS_NUM_THREADS=1` avoids unnecessary thread overhead on these small calculations. `larger_response_check.py` currently reruns the baseline study when importing its Hamiltonian construction.

Validation includes the exact eigensolver residual, analytic versus finite-difference variational gradients, independent sparse versus Fourier kinetic operators, a lifted full-coordinate check of the nonzero-momentum Hamiltonian, local gate-matrix checks including the control-off branch and clean workspace, exhaustive legal-coordinate triple-equality checks, and gate recounting after decomposition. The physical-sector dynamics are numerically simulated; a dense 36-qubit statevector is not. Exact Clifford+T synthesis and a physical-device schedule remain unperformed.

**Next manuscript decision:** adopt this larger example as a proposed 2028 short-time charge-response benchmark, with its model limitations and complete circuit scope explicit. It should replace the one-chiral-step claim rather than inherit that claim's experimental-validation language. The 2033 target still needs a separate audit of improved volume, interaction, and energy resolution.
