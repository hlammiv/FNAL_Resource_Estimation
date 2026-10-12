# Chapter 2: a near-1000-qubit 2033 alternative for discussion

The user accepted the 36-qubit, four-nucleon 2028 candidate on 2026-10-04. This note investigates a separate 2033 model. The manuscript and its resource tables have not been changed.

## Candidate and decision

A two-channel, two-nucleon weak-response benchmark on a periodic 6×6×6 lattice has a construction-level estimate of **887 logical qubits and 5.21×10^8 T per deepest circuit**, including deterministic source preparation. Adding an illustrative ten T per arbitrary rotation gives 5.79×10^8. These are leading synthesis estimates, not exact compiled counts or certified hardware fits.

The scientific tradeoff is material: this is a larger box and much longer evolution, but fewer nucleons than the proposed 2028 calculation. Moreover, 864 of these qubits represent the full four-mode occupation lattice. The particular charged-current transition studied here can be restricted to two active modes per site; a compact position encoding is smaller still. Thus **887 is an implementation choice, not a demonstrated minimum or a reason this problem needs a 1000-qubit machine**. If the chapter needs a many-body progression, a four-nucleon response on the same lattice is the better next investigation, but its preparation has not been validated or costed here. The two-body loader and error results cannot be reused as a four-body estimate.

## Physics definition

Use nearest-neighbor nonrelativistic hopping t=10.5794 MeV at a=1.4 fm, and independent on-site contacts in the spin-triplet/isospin-singlet and spin-singlet/isospin-triplet channels. The box is 8.4 fm on each side. All six antisymmetric local two-nucleon channels are included in the four-mode Hamiltonian. There is no three-body force, Coulomb interaction, effective-range correction, tensor force, or two-body axial current in this benchmark.

Fit the triplet contact to infinite-volume binding B=2.224575 MeV and the singlet contact to a_s=−23.76 fm using the lattice Green function. With relative dispersion ε(k)=4t Σ_i(1−cos k_i), normalized Brillouin-zone average G(E)=⟨1/(E−ε)⟩, the conventions are

    1/C_t = G(−B),       1/C_s = G(0) + a/(8π t a_s).

Numerical quadrature gives C_t=−93.62604 MeV and C_s≈−82.19588 MeV. The singlet contact changes by about 0.00035 MeV between quadrature orders 100 and 160; the quoted calculation retains the finite-quadrature value, not an exact continuum integral. These are newly fitted benchmark couplings, not the four-body 2028 couplings.

Prepare the zero-momentum triplet ground state with spin projection +1. The normalized, zero-momentum-transfer one-body charged axial insertion converts its n↑ into p↓, leaving a spatially symmetric pp singlet source. In this fixed sector it has unit norm and can be implemented by a known fermionic flavor permutation, or absorbed into the output labels of the source loader. The response is

    C(τ) = ⟨ψ_t| exp[−i(H_s−E_t)τ] |ψ_t⟩.

Known weak/spin normalization factors are stripped off; this is a normalized channel correlator. The zero-momentum-transfer insertion is the allowed-transition limit, not a full neutrino-kinematics calculation. Without Coulomb and higher-order operators, it is not a prediction of the physical ν_e+d cross section. Published neutrino-deuteron studies provide the scientific motivation and demonstrate why more complete currents matter: [Acharya and Bacca, PRC 101, 015505](https://journals.aps.org/prc/abstract/10.1103/PhysRevC.101.015505).

## Observable and numerical checks

Compute real and imaginary parts at τ=0.05, 0.10, 0.25, and 0.50 MeV⁻¹ (9.87–98.66 fm/c). This extends the longest 2028 time by a factor of 16.7. It can check interaction-dependent spectral interference, but four time samples alone do not reconstruct a spectrum. No energy-resolution claim or reconstructed cross section is included.

For L=6 the initial energy is −4.14722 MeV and the lowest final energy is −2.01216 MeV. The first two transition energies/strengths are 2.13506 MeV/0.92895 and 14.32518 MeV/0.06371. The lowest final state is below zero even though the infinite-volume singlet is unbound. This is a finite-volume level, not evidence for a bound physical diproton. At L=4 the initial energy is −7.40905 MeV, confirming substantial remaining volume effects. A volume extrapolation is not claimed.

Exact spectra use the 216-dimensional zero-total-momentum relative-coordinate space. Product-formula validation evolves all 216²=46,656 two-position amplitudes, so intermediate even/odd bond matchings are not incorrectly forced to preserve one-site translation symmetry. Six kinetic bond matchings are traversed forward and backward around the contact evolution in each symmetric step.

| τ (MeV⁻¹) | Steps at Δτ=0.0025 | Absolute complex correlator error |
|---:|---:|---:|
| 0.05 | 20 | 0.00178 |
| 0.10 | 40 | 0.00312 |
| 0.25 | 100 | 0.00458 |
| 0.50 | 200 | 0.00997 |

Halving the step gives 0.00250 error at the last time, consistent with second-order convergence. These are numerical errors for the specified source and observable, not operator-norm error bounds. Gate-synthesis and hardware errors are separate.

## Complete-circuit construction estimate

The source is classically computed deliberately: this is an independently checkable benchmark, not a quantum-advantage claim. Embed the full two-position wavefunction into two eight-bit registers. A generic real-state loader on 16 bits uses at most 2^16−1=65,535 arbitrary rotations with Clifford multiplexing. Convert the two positions to two one-hot occupation registers with 432 eight-controlled NOTs. Six clean scratch bits suffice for a 13-Toffoli construction of each; erase the binary labels with occupation-controlled CNOTs. The source is prepared deterministically, with no success probability or overlap amplification hidden in the estimate.

| Resource | Count |
|---|---:|
| Full four-mode system register | 864 |
| Binary source-loader registers | 16 |
| Clean source-loader scratch | 6 |
| Hadamard ancilla | 1 |
| Conservative allocated logical qubits | **887** |
| Controlled hopping rotations per symmetric step | 20,736 |
| Controlled contact rotations per step | 7,776 |
| Complete deepest-circuit arbitrary rotations | 5,767,936 |
| Loader Toffolis | 5,616 |
| Leading complete-circuit T estimate | **520,944,149** |
| Illustrative extra ten T per rotation | **578,623,509** |

There are 3×216 bonds per flavor, four flavors, two Pauli rotations per hopping gate, two kinetic half-sweeps, and two elementary rotations per controlled Pauli rotation. The on-site interaction has 18 nonidentity Pauli terms. Its diagonal sum commutes with its off-diagonal sum, and the eight off-diagonal Pauli terms mutually commute; their ordered exponentials therefore implement the contact exactly. The script verifies these commutators. Jordan–Wigner strings and basis changes add Clifford gates; no physical routing schedule is inferred. Scalar kinetic/contact/reference-energy phases are combined into one controlled phase for the fixed two-particle sector.

Counts charge seven T per Toffoli and 3 log₂(N_rot/0.005) T per rotation. The latter is only the leading [Ross–Selinger](https://arxiv.org/abs/1403.2975) estimate. The 0.005 total coherent synthesis allowance can add at most 0.01 to a bounded Hadamard expectation, separately from the product-formula error. No full 887-qubit QASM or exact Clifford+T synthesis has been generated. There is ample estimated count headroom, but physical failure rates, Clifford depth, factories, and routing remain unpriced.

For illustration, eight bounded estimates (real/imaginary at four times), each to ±0.05 statistical error with 95% simultaneous coverage, require ceil[2 log(2×8/0.05)/0.05²]=4,615 shots each by Hoeffding plus a union bound: 36,920 shots. Applying the deepest cost to every shot gives a conservative leading gate-volume estimate of 1.92×10^13 T. This is distinct from the per-circuit 10^9 budget and is not a hardware runtime estimate. This modest first campaign does not precisely determine the small excited-state strength.

## Reproduction

    OPENBLAS_NUM_THREADS=1 python3 editorial_review/2033_model_study/two_body.py
    OPENBLAS_NUM_THREADS=1 python3 editorial_review/2033_model_study/resources.py

The JSON files preserve spectra, quadrature checks, tested errors, and all resource arithmetic. No manuscript replacement is proposed without discussing the fewer-particles / larger-box tradeoff and the choice of occupation encoding.
