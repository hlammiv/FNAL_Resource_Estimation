# Fully mobile triton response: proposed 2033 benchmark

Prepared 4 October 2026 following the user's instruction to investigate triton. **Study only: the manuscript and its resource tables remain unchanged.**

## Result for discussion

A specified three-mobile-nucleon calculation on a 6×6×6 lattice uses a conservatively allocated **908 logical qubits** and an estimated **7.23×10^8 T** for its deepest complete correlator circuit. An illustrative extra ten T per arbitrary rotation raises this to **8.02×10^8 T**. Preparation, current insertion, controlled evolution, and readout are included. Exact Clifford+T synthesis and a device schedule have not been performed; these are constructive resource estimates, not certified hardware fits.

The calculation has nontrivial inelastic strength: **29.35%** of the unit-normalized proton-charge response lies outside the recoiling ground-state pole. At the longest proposed time, **39.47 fm/c**, the tested product-formula bias is **0.00731 absolute**, or **1.13% of the full correlator magnitude**. Relative to the elastic-subtracted correlator magnitude it is **3.30%**. Gate synthesis, statistical errors, and logical noise are separate.

The main physics limitation is substantial finite-volume dependence. This is a **finite-volume triton-inspired benchmark**, not a converged physical breakup cross section. The coupling calibration and observable below make that distinction explicit.

## Published foundation and the extension

[Roggero et al., PRD 101, 074038 (2020)](https://arxiv.org/abs/1911.06368) gives a lattice nuclear Hamiltonian, quantum resource analysis, and a triton-inspired demonstration with one static nucleon on a tiny two-dimensional lattice. [Baroni et al., PRD 105, 074503 (2022)](https://doi.org/10.1103/PhysRevD.105.074503) computes real-time response correlators for that simplified model.

This study makes all three particles mobile in three dimensions and distinguishes the singlet and triplet two-body contacts. Its resource numbers are new calculations, not numbers taken from those papers. The leading-order three-body interaction has a standard EFT motivation; see [Bedaque, Hammer and van Kolck](https://arxiv.org/abs/nucl-th/9906032).

## Precisely specified model

- One proton and two identical neutrons; total spin projection +1/2. The calculated ground state is in the total-spin-1/2 sector, verified numerically. Neutron exchange is antisymmetric.
- 216 spatial sites, periodic boundary conditions, lattice spacing 1.4 fm, box side 8.4 fm.
- Nearest-neighbor hopping t=10.5794 MeV. Each particle's kinetic diagonal is 6t.
- Two on-site pair channels: C_t=−93.62604 MeV and C_s=−82.19588 MeV. These are the preceding two-body study's lattice-Green-function fits to infinite-volume deuteron binding 2.224575 MeV and isospin-symmetric singlet scattering length −23.76 fm. The singlet quadrature uncertainty was checked separately in that study.
- Three-body interaction D Σ_x Σ_{α<β<γ} n_{xα}n_{xβ}n_{xγ}, with **D=126.95409 MeV fitted to E_0=−8.4818 MeV in the L=6 box**. This is explicitly a finite-box renormalization condition, not an infinite-volume fit.
- No effective-range correction, tensor force, or two-body electroweak current. The singlet coupling is isospin symmetric; nucleon mass differences are omitted.

The spatial calculation uses two relative vectors, r_n1−r_p and r_n2−r_p, and three spin components labeled by the particle with down spin. Its working dimension is 3×216²=139,968 before symmetry reduction. The center of mass is translation invariant initially. After the probe the total momentum is nonzero; evolution is calculated in that momentum sector, with proton recoil included.

All four occupation-mode types participate in the spin-dependent model. In the computed ground state, the p↓n↑n↑ spin component carries about 0.806% probability; the other two labeled spin components each carry about 49.597%. Thus the fourth mode is not simply an unused register added to enlarge the count. Nevertheless, occupation encoding is not qubit minimal: a first-quantized implementation could use much less than 908 qubits. No claim that this problem intrinsically requires 1,000 qubits or gives quantum advantage is made.

## Observable and useful physics

The probe is the proton charge density at q=(2π/(6a),0,0), approximately 148 MeV. With exactly one proton, its insertion J(q)=exp(iq x_p) is unitary. The measured observable is

    C(τ,q) = ⟨ψ_0|J† exp[−i(H−E_0)τ] J|ψ_0⟩.

One Hadamard test measures either its real or imaginary part with bounded ±1 outcomes and no postselection. The proposed first time grid is τ=0.05, 0.10, 0.20 MeV⁻¹, corresponding to 9.87, 19.73, 39.47 fm/c. The largest time is 6.67 times the proposed 2028 endpoint.

The moving ground-state energy is −4.31648 MeV, so the elastic recoil energy is 4.16532 MeV. Its response weight is 0.706545. Subtracting the pole gives inelastic strength 0.293455. In this checkable benchmark, the elastic energy and weight are independently calculated classically; the resource estimate is for the quantum correlator plus classical pole subtraction. A separate quantum determination of the elastic form factor is not included or claimed.

The finite-volume deuteron-plus-neutron threshold diagnostic at this total momentum is about +1.30069 MeV in absolute energy, above the moving trimer. The response is sensitive to interacting excited states, with substantial weight beyond the trimer pole. This is not an extraction of separate d+n and p+n+n yields or their angular distributions. Finite-volume levels and a few correlator points cannot simply be called an infinite-volume breakup spectrum.

The first nonzero momentum is near the pionless breakdown scale; therefore the model should not be presented as a controlled high-accuracy EFT prediction at this momentum. It remains a specified nuclear-model test of preparation, charge insertion, recoil, spin-channel mixing, and three-body dynamics.

## Numerical checks

| τ (MeV⁻¹) | Steps, Δτ=0.002 | Exact full C | Absolute product-formula error |
|---:|---:|---|---:|
| 0.05 | 25 | 0.877022 − 0.337544i | 0.002046 |
| 0.10 | 50 | 0.666322 − 0.530120i | 0.003949 |
| 0.20 | 100 | 0.264612 − 0.592591i | 0.007313 |

Evolution uses symmetric V/2–K–V/2 steps, with **exact kinetic evolution**, not nearest-neighbor even/odd Trotter splitting. The circuit therefore uses fermionic Fourier basis changes, matching the numerical FFT calculation. Contacts are also exponentiated exactly within each potential layer.

Halving Δτ reduces the endpoint error from 0.007313 to 0.001827, consistent with second-order convergence. The ground-state residual is below 10^−10 MeV; the moving-state residual is below 10^−10 MeV. Independent real-space hopping and momentum-space kinetic actions agree to approximately 1.5×10^−13. Increasing the reference Krylov space from 100 to 140 vectors changes the tested correlators by less than 2×10^−15. These are numerical convergence checks, not rigorous operator-norm error bounds.

An independent 2³-box calculation builds the full 736-dimensional occupation-basis Hamiltonian with explicit fermionic hopping signs and local spin/isospin interactions. It agrees with the relative-coordinate implementation's ground energy to 1.5×10^−14 MeV and its charge correlators to 1.3×10^−14. This checks the spin-channel mapping, identical-neutron treatment, current, and recoil together; it does not simulate a dense 908-qubit statevector.

See [the discussion figure](triton_check.png) for the full and elastic-subtracted correlators and the volume dependence.

The same couplings give:

| Lattice | Box side | Ground energy |
|---|---:|---:|
| 4³ | 5.6 fm | −15.12348 MeV |
| 6³ | 8.4 fm | −8.48180 MeV, fitted |
| 8³, classical check | 11.2 fm | −6.44750 MeV |

This visible drift is why the matching energy must not be advertised as a prediction of physical triton binding. Larger-volume matching and continuum/volume analysis would be further physics work.

## Deterministic preparation construction

The exact classical reference wavefunction is embedded in a 20-bit register: six three-bit relative coordinates and two spin-label bits. A generic real-amplitude loader needs at most 2²⁰−1 arbitrary rotations. Three uniform six-valued center coordinates add at most 21 rotations. Saved conditional rotation angles reconstruct the signed amplitudes with norm error below 4×10^−16.

Reversible modulo-six translations form the three absolute positions. Spin labels are decoded; the two neutron position-spin keys are reversibly sorted before conversion to occupations. Because the input is antisymmetric, the sorting flag factors into the known |−⟩ state and can be reset with Clifford gates. Identical-neutron-mode amplitudes vanish. This avoids incorrectly mapping two distinguishable labeled particles into one fermionic occupation state.

Occupation writes are controlled on each particle's position and spin. Binary labels are then erased from the occupied modes: proton labels by CNOTs; sorted neutron labels by a scan using the parity of preceding neutron occupations to distinguish first and second. The scan parity returns to zero for two neutrons. There is no ground-state filtering, preparation-overlap assumption, or rejected-shot multiplier.

These reversible operations have an explicit conservative allowance of 60,066 Toffolis. The full large loader and 908-qubit circuit have not been exported or independently compiled; the checked objects are the loader amplitudes, exchange symmetry, local Hamiltonian algebra, Fourier decomposition, and construction-level operation counts.

## Register and gate accounting

| Allocated register | Qubits |
|---|---:|
| Four occupation modes on each of 216 sites | 864 |
| Three absolute position registers | 27 |
| Spin labels and decoding | 5 |
| Clean arithmetic/comparator scratch | 10 |
| Neutron sorting flag | 1 |
| Hadamard ancilla | 1 |
| **Total** | **908** |

For each kinetic layer, a separable 3D Fourier transform uses 108 six-mode blocks per flavor. A conservative generic decomposition of a block uses 15 complex Givens rotations at at most four arbitrary rotations each, plus six phase rotations. Basis change and inverse across four flavors cost 57,024 rotations. Controlled momentum phases add 1,728. The basis changes themselves need not be controlled: their control-off product cancels.

The complete on-site interaction has 22 nonidentity Pauli terms, including the three-body contact. Its diagonal sum commutes with its off-diagonal sum, and the off-diagonal Pauli strings mutually commute; these facts are checked in the script. A controlled potential layer costs 9,504 rotations. Across n symmetric steps, adjacent potential half-layers merge, leaving n+1 potential layers.

| Complete circuit | Arbitrary rotations | Leading T | With illustrative +10 T/rotation |
|---|---:|---:|---:|
| τ=0.05, 25 steps | 2,764,933 | 2.41×10^8 | 2.69×10^8 |
| τ=0.10, 50 steps | 4,471,333 | 3.99×10^8 | 4.44×10^8 |
| τ=0.20, 100 steps | 7,884,133 | **7.23×10^8** | **8.02×10^8** |

All rows include the same fresh preparation, charge insertion, and 60,066-Toffoli allowance. The count uses seven T per Toffoli and the leading synthesis estimate 3 log₂(N_rot/0.005) per rotation. As in the preceding studies, the [Ross–Selinger](https://arxiv.org/abs/1403.2975) subleading synthesis term is not rigorously bounded by the illustrative +10 allowance. A total coherent synthesis allowance of 0.005 permits up to 0.01 additional bias in a bounded Hadamard expectation. The 1.13% numerical evolution error is not a total experimental error budget.

The 200-step refinement used to check numerical convergence would cost about 1.39×10^9 leading T with this construction, above the 2033 target. It is a classical verification run, not the proposed budgeted circuit. The generic Fourier and loader constructions leave possible optimization opportunities, but no savings from those are assumed.

## Repetitions and limits of the fit

A first result at τ=0.10, estimating both components to ±0.05 statistical error with 95% simultaneous coverage, needs 3,506 shots per component using Hoeffding plus a union bound: 7,012 complete executions. Its leading gate volume is about 2.80×10^12 T. This is distinct from the per-circuit 10^9 budget.

A three-time campaign with ±0.02 per component and 95% simultaneous coverage uses 27,404 shots per component, or 164,424 executions. Charging the appropriate circuit for each time gives approximately 7.47×10^13 T. No physical runtime follows without Clifford depth, routing, magic-state supply, feed-forward, and logical-error assumptions. The construction fits the estimated algorithmic width and per-circuit T target; it does not establish a device-level campaign schedule.

## Reproduction and proposed manuscript decision

    OPENBLAS_NUM_THREADS=1 python3 editorial_review/2033_model_study/triton.py
    OPENBLAS_NUM_THREADS=1 python3 editorial_review/2033_model_study/triton_resources.py
    OPENBLAS_NUM_THREADS=1 python3 editorial_review/2033_model_study/triton_crosscheck.py
    MPLCONFIGDIR=/tmp/triton-matplotlib OPENBLAS_NUM_THREADS=1 python3 editorial_review/2033_model_study/triton_plot.py

Inputs include the independently saved two-body contact-fit results. Outputs preserve model parameters, fitting history, states, conditional loader angles, error checks, and resource arithmetic.

**Recommendation:** discuss this as the 2033 companion to the accepted 2028 benchmark, under a title such as “Fully mobile three-nucleon charge response on a 6³ lattice.” Do not label the first deliverable an exclusive triton-breakup cross section or a volume-converged neutrino prediction. No manuscript change has been made pending that discussion.
