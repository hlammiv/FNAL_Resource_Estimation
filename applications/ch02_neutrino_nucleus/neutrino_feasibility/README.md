# Candidate 2028 four-nucleon response benchmark

**Discussion draft, not a manuscript revision.** This study constructs and checks a complete arbitrary-rotation circuit for a small interacting four-nucleon model. Its Clifford+T resource estimate suggests a fit near 10^5 T gates, but exact rotation synthesis remains to be done. It does not establish a physical helium response or quantum advantage.

**Follow-up:** the [3×3×3 study](LARGER_BOX.md) now provides a larger candidate with deterministic structured preparation. The original trial below is retained for comparison.

## Scientific question and model

Can a fault-tolerant processor reproduce the change in a short-time charge correlator caused by nuclear interactions in a classically checkable model?

Take one particle of each of the four spin–isospin species (two protons and two neutrons), on a periodic 2×2×2 cubic lattice. This is an alpha-particle-sector model, with SU(4)-symmetric pionless contact interactions. In this fixed-species sector each particle position uses three qubits. No exchange between identical particles of the same species is needed because there is only one of each; the encoding represents the corresponding fermionic occupation sector.

Use t=10.5794 MeV, C0=−98.2265511 MeV and D0=127.839693 MeV, at a=1.4 fm, from the lattice Hamiltonian and parameter table in [Roggero et al., arXiv:1911.06368](https://arxiv.org/abs/1911.06368). The potential is C0 times the number of coincident pairs plus D0 times the number of coincident triples. The kinetic term is −2t times the sum of the twelve position-bit X operators. The factor two keeps both directed periodic bonds when the two neighbours coincide on an L=2 lattice. An additive 24t is omitted from both H and E0, leaving H−E0 unchanged.

The box is only 2.8 fm across. Its nonzero axial momentum q=π/a is 443 MeV, at the lattice cutoff. This is explicitly a regulator-scale nuclear-model benchmark, not a controlled pionless-EFT prediction at physical scattering kinematics. The calculated energy, −95.02 MeV including the kinetic constant, is not a prediction for the physical helium binding energy. The source couplings are reused, not refitted to this tiny volume. These are substantial scientific limitations, not just resource caveats.

## Observable and complete circuit

Define J=ρ(q)/2=(Z0+Z3)/2, using the x-position bits of the two protons. Measure

C(τ)=⟨ψ0|J exp[−i(H−E0)τ] J|ψ0⟩.

The four terms Z_a U Z_b are unitary and can each be measured by a Hadamard test. Uniformly sample the proton pair (a,b); the ±1 outcome estimates C directly after the 1/4 average. There is no unspecified current normalization, amplitude amplification, or postselection. The unnormalized charge correlator is 4C. Real and imaginary parts use separate measurement settings.

The ground state is calculated classically, as appropriate for validation. Its preparation is included in every shot. Translation symmetry factors it into a positive real state on nine relative-coordinate qubits and a uniform three-qubit origin. A uniformly controlled Ry tree uses 511 elementary rotations, followed by nine CNOTs to restore absolute coordinates. The implementation explicitly generates and checks the loader; the method is based on [Möttönen et al.](https://arxiv.org/abs/quant-ph/0407010).

Time evolution uses the symmetric second-order product formula with exact commuting kinetic and potential layers. The potential has 210 nonidentity diagonal Pauli strings, obtained by a Walsh transform rather than a scaling guess. At two steps, merged kinetic half-steps give 36 Pauli rotations and the potential gives 420. Each controlled Pauli rotation is decomposed into two elementary rotations and Clifford gates. Two additional rotations retain the potential's identity phase and one applies the known E0 phase on the Hadamard ancilla.

Total: **13 algorithmic logical qubits; 1,426 elementary rotations**, including the 511 preparation rotations. No magic-state factory or error-correction overhead is included in the register. There is no need to fill the entire 150–250-qubit reference machine.

## Numerical result and accuracy

Exact sparse diagonalization and exact time evolution are compared with the product formula. At τ=0, the normalized strength is **0.631216**, compared with **0.500000** in the free, uncorrelated ground state: a 26.2% enhancement. At τ=0.01 MeV^−1 (1.97 fm/c),

- Exact interacting C: **0.595555 − 0.180132i**.
- Two-step circuit C: **0.593646 − 0.185290i**.
- Free ground-state C: **0.455895 − 0.205329i**.

The two-step complex error is **0.005500**, or 0.88% of |C|. Four and eight steps reduce it to 0.001346 and 0.000335, checking second-order convergence. At earlier times 0.0025 and 0.005 MeV^−1 the two-step errors are smaller. The proposed result is three short-time correlator points plus a static-strength check, not a reconstructed spectrum. The free comparison changes both the Hamiltonian and its ground state. A same-interacting-H uncorrelated diagnostic is also recorded in the JSON, but is not called a ground-state response.

## Gate estimate and remaining qualification

Using a total coherent rotation-approximation allowance of 0.005 gives a per-rotation tolerance of 0.005/1426. The existing report's Ross–Selinger leading estimate, 3 log2(1/ε) T per rotation, gives **about 77,500 T** for the entire shot. The unitary approximation allowance implies at most 0.01 absolute bias in a bounded Hadamard observable, separately from Trotter error and logical faults.

This is **not an exact compiled T count or a rigorous T-count upper bound**: the leading synthesis formula omits its subleading contribution. Adding an illustrative 10 T per rotation would give about 91,800 T; that sensitivity check is not a bound either. [Ross–Selinger](https://arxiv.org/abs/1403.2975) explicitly gives a subleading term. The report's randomized repeat-until-success convention gives approximately 27,200 expected T, but has different error/runtime assumptions and is not the conservative headline used here.

For six estimates (real and imaginary parts at three times), Hoeffding's bound gives 27,404 shots per estimate for ±0.02 statistical error with at least 95% simultaneous coverage, using a union bound. That is 164,424 shots and about 1.27×10^10 T at the leading per-shot estimate. Preparation is repeated and included in all of them. A static-strength sample set is additional; it needs only preparation and computational-basis measurement. Logical noise, Clifford scheduling, factory throughput, and wall time require a hardware model; no wall-time claim is made here.

## Reproduction and next decision

Run `python3 editorial_review/neutrino_feasibility/study.py`, then `python3 editorial_review/neutrino_feasibility/circuit_check.py` from the manuscript directory. The first requires NumPy/SciPy and uses the repository's synthesis conventions. The second requires Qiskit. The scripts save numerical results, loader amplitudes/angles, and an explicit unsynthesized QASM circuit. `response_check.png` visualizes the numerical comparison.

Checks include the eigenvector residual, exact reconstruction by the relative-coordinate loader, the diagonal Hamiltonian's Pauli decomposition, all four Hadamard terms, and agreement between the assembled 13-qubit statevector circuit and the independently implemented product formula. Rotation synthesis itself is not validated by those checks.

**Recommendation:** retain this as a concrete low-cost validation candidate, subject to agreement that its tiny box and simplified interaction provide enough scientific content for the 2028 example. It cannot inherit the current chapter's physical-data comparison or precision-response claims. A 2033 example should increase physical fidelity and control volume/momentum artifacts, not merely run this same toy longer. No existing chapter, resource model, summary table, or figure has been changed by this study.
