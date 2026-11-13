# Independent check of the first-quantized scoping study

Date 2026-11-12. Read-only check of `README.md`, `scoping.py`, `spagnoli_lemmas.py`, `spin_contact_check.py`,
`results.json`. Nothing in the study, the model or the paper was changed. `scoping.py` was not imported. The verification scripts
(`v_table5.py`, `v_counts.py`, `v_spin_contact.py`, `v_gqsp40.py`, each under a minute on two threads) were run outside the repository and are not shipped.
Source paper: arXiv:2507.22814, fetched and read as text.

## Verdicts

| claim | verdict | what to change |
|---|---|---|
| (a) triton, 48 LQ, 2.6e7 T, 49 d | confirmed | It is the 8^3 box, not "the same triton" as the 6^3 box. 46 to 48 LQ depending on how the equality ladders are held. |
| (b) 12C, 160 LQ, 3.6e8 T, 4.2 yr, 66 d | arithmetic confirmed to 0.1%; not supported as a feasibility statement | The assumed quantity is the overlap squared, 0.1 (amplitude 0.32). An overlap of 0.1 gives 28 yr. The centre-of-mass motion of the trial state is not treated; see section 7. |
| (c) 40Ar exceeds 1e9 T | confirmed, conditional | Only the ground-state filter exceeds 1e9. The evolution is 1.5e8 to 2.1e8 T. Floor 3.0e9 T at the assumed 1.46 MeV gap. |
| (d) commuting form of the spin-dependent contact | confirmed | Holds on Pauli-allowed states only. The chapter model satisfies the condition. |
| (e) Table 5 to about 1% for the qubitization rows | confirmed | QSP and GQSP within 1.2%. The product-formula rows are 5% to 44% high, as the README says. |

## 1. Registers

Triton, species-labeled (p, n, n), 8^3.
System: 3 x (9 coordinate + 1 spin) = 30. The loader's 20 qubits (six 3-bit relative coordinates, 2 spin-label bits),
the 9 centre bits and 1 spin-decoding bit become these 30; no qubit is left over. No sorting flag is needed because
the loaded wavefunction is already antisymmetric in the two neutron registers.
Workspace: 3 pair flags, a 7-qubit AND ladder for the 9-bit equality test, 1 triple flag, 3 for the AND of
(control, exclusive flag, two spin bits), 1 synthesis, 1 Hadamard = 16. Total 46.
The study has 48 (8 ladder, 5 AND chain). Agreement within 2.

One inconsistency. The study's 62 Toffolis per potential layer count 8 Toffolis to compute and 8 to uncompute each
flag. That needs the ladder of each flag kept alive, so three ladders (21 qubits) with three stored flags, which gives
60 LQ. With one shared ladder each flag costs 15 + 15 Toffolis and the layer costs 110. So: 46 to 48 LQ at about 110
Toffolis per layer, or 60 LQ at 68. The per-shot T count is the same to three digits (2.612e7 or 2.615e7).

12C, species-labeled, 8^3.
System: 12 x 10 = 120. Workspace: 24 (2-bit occupancy counter per particle), 8 (ladder including the flag, kept
alive while the flag is used, consistent with 16 Toffolis per equality), 5 (AND of control, flag, two counter bits, two
spin bits), 3 (synthesis, Hadamard, filter flag) = 40. Total 160. Agrees.
16O 208 and 40Ar 496 follow the same formula, 12 A + 16.

The labeled register is valid. The contact acts on spin only, so an isospin label is never changed. Checked
numerically below: the labeled spectrum equals the (Z, N) sector of the occupation-mode Hamiltonian.

## 2. T counts

Own derivation of each layer, then `estimates/common.py` for c_T (RUS, eps_syn = 1e-2) and 7 T per Toffoli.

Kinetic layer, per coordinate: QFT_8 and inverse contain 4 controlled-S (3 T each) and 2 controlled-T (3 synthesized
rotations each). cos(2 pi p/8) has 4 Walsh terms (verified), each controlled by the Hadamard qubit, 2 rotations each.
Total 14 rotations and 12 T. Agrees.

Potential layer, A >= 12, per pair: three sweeps of an equality test at 8 + 8 Toffolis, 2 + 2 Toffolis for the two
counter increments and their inverse, 5 + 5 for the six-input AND, one rotation. Per particle: three controlled
Walsh phases on the 2-bit counter, 6 rotations. Total P (6 x 8 + 14) Toffolis and P + 6 A rotations. Agrees.
The scalar part of exp(-i th SWAP) folds into g(1), so it adds nothing.

Output of `v_counts.py` (ratio = mine / study):

| component | 3H 8^3 mine | ratio | 12C 8^3 mine | ratio |
|---|---|---|---|---|
| kinetic layer | 3,199 T | 1.000 | 13,444 T | 1.000 |
| potential layer | 770 T (68 Toffoli) or 1,064 T (110) | 1.13 or 1.57 | 32,207 T | 1.000 |
| one step | 3,969 or 4,263 T | 1.02 or 1.10 | 45,651 T | 1.000 |
| loader / Slater preparation | 2.572e7 | 1.000 | 9.55e5 (not rederived, see below) | 1.000 |
| filter per attempt | none | | 3.450e8 | 0.999 |
| evolution to 0.20 | 3.98e5 or 4.27e5 | 1.02 or 1.10 | 9.16e6 | 1.000 |
| deepest circuit | 2.612e7 or 2.615e7 | 1.000, 1.001 | 3.551e8 | 0.999 |
| expected T per accepted shot | same | | 8.123e8 | 0.999 |
| first result | 2.10 d | 1.000 | 65.6 d | 0.999 |
| campaign | 49.4 d | 1.000 | 4.21 yr | 0.999 |

No factor-of-two disagreement. The largest is the triton potential layer (1.13 to 1.57, cause above); it moves the
per-shot count by 0.1%. The 0.1% on the 12C filter is the mean of |t| for a Gaussian, 0.798 sigma against the study's 0.8.
16O: 5.05e8 T, 5.96 yr. 40Ar: 1.77e10 T, 211 yr. Both agree.

Slater preparation. I did not rederive the count of arXiv:2509.07279. A floor of one 511-rotation orbital load per
particle is 1.6e5 T; the study has 9.5e5 T. Either is under 0.3% of the 12C count.

Cross-check against the paper's own lemmas converted to 7 T per Toffoli and the report's c_T (SU(4), uncontrolled):

| A | Lemma 3 potential step | study-type potential step | Lemma 1 kinetic step (Eq. (19), b = 25) | hopping / QFT_8 kinetic step |
|---|---|---|---|---|
| 3 | 560 T | 770 T | 8,442 T | 3,199 T |
| 12 | 42,495 T | 32,207 T | 33,768 T | 13,443 T |

The study's kinetic layer is 2.5 times below the converted Lemma 1. The cause is a different circuit for a different
dispersion: the chapter's hopping term through an exact QFT_8 and four Walsh phases, against squaring plus phase
kickback for q^2. This is legitimate for the chapter Hamiltonian. With Lemma 1 the 12C step would be 6.6e4 T, 1.44
times the study's, and every 12C number would scale by that factor.

## 3. Claim (d)

Analytic. For n nucleons on one site with total spin S, the sum of pair SWAPs is S(S+1) + n(n-4)/4. It vanishes for
n = 3 only if S = 1/2 and for n = 4 only if S = 0. With four spin-isospin states per site and antisymmetry, the
on-site state is spatially symmetric, so n = 3 has S = 1/2 and n = 4 has S = 0, and n >= 5 cannot occur. The same
holds when antisymmetry is imposed within each species only. Then V = sum_i g(n_i) + dC sum e_ij SWAP_ij.
Two exclusive-pair terms either share a particle, and then have disjoint position support, or act on different spins.
They commute with each other and with the diagonal part on the whole space.

Condition, stated precisely: the identity holds on states in which every site with three or more nucleons has the
minimal total spin. The Pauli principle with exactly four internal states per site guarantees this. It fails on the
unsymmetrized space, by 3 dC for a spin-3/2 triple and 6 dC for a spin-2 quadruple. No hard-core or occupancy
constraint beyond Pauli is needed. The chapter model (four modes per site, C_t, C_s, D) satisfies it. The layer is
exact only if the prepared state is antisymmetric within each species; the loader and the Slater step must supply that.

Numerical (`v_spin_contact.py`; second-quantized reference built by me with Jordan-Wigner and channel projectors):

| case | (V_full - V_excl) on allowed states | same on the full space | exclusive-term commutators | plain pair commutators | spectrum vs occupation model | layer as product of gadgets vs exp(-i V dt) |
|---|---|---|---|---|---|---|
| 2 nucleons, 3 sites | 0 | 0 | 0 | 0 | | 1e-15 |
| 3 nucleons, 2 sites | 1.3e-14 | 17.15 = 3 dC | 0 | 32.7 | 2e-13 | 1.7e-15 (plain pair product: 1.4e-2) |
| 3 nucleons, 3 sites | 1.9e-14 | 17.15 | 0 | 32.7 | | 2.9e-15 (plain: 1.7e-2) |
| 4 nucleons, 2 sites | 6.8e-14 | 34.29 = 6 dC | 0 | 216 (random vectors) | 4e-13 | |
| labeled p n n, 2 sites | 2e-16 | | | | 3e-14 | |
| labeled p p n n, 2 sites | 6e-14 | | | | 3e-13 | |
| labeled p p n n n, 2 sites | 9e-14 | | | | 3e-13 | |

V_excl maps the antisymmetric space into itself (leakage 5e-13). Levels found: C_t, C_s, 3 Cbar + D = -136.779,
6 Cbar + 4 D = -19.649 MeV.

The count P (6 (dm-1) + 14) Toffolis and P + 6 A rotations follows from the construction (section 2). It is a
construction-level count. No circuit was compiled, by the study or by me.

## 4. Shots and wall times

- Shot counts follow the chapter's Eq. (Nshot): 27,404 per estimate (m = 6, eta = 0.02, alpha = 0.05), 164,424 for the
  campaign, 2 x 3,506 = 7,012 for the first result (m = 2, eta = 0.05). Reproduced.
- These are accepted shots. Attempts per accepted shot are 1/overlap^2 = 10. A rejected attempt is charged the Slater
  step plus 2.02 of 14 rodeo cycles (rejection at the first failed cycle, per-cycle survival 0.506). Expected T per
  accepted 12C shot: 9 x 5.07e7 + 3.55e8 = 8.1e8. Reproduced.
- Arithmetic: 12C campaign 2 x 27,404 x (sum over three times of 8.0e8 to 8.1e8) = 1.33e14 T = 1.33e8 s = 4.21 yr.
  First result 7,012 x 8.08e8 = 5.66e12 T = 65.6 d. Triton: 4.27e12 T = 49.4 d and 1.82e11 T = 2.1 d. Serial, one machine.
- For 12C the insertion is an average over Z^2 = 36 unitary pairs, sampled uniformly as in the 2028 box. The bound on the
  estimator is unchanged, so the shot count is the chapter's. The study counts one insertion per shot; a pair needs two.
  The difference is 9 gates.

## 5. Assumptions

All five are tagged Assumed in the README table. Whether the values are favorable:

- Overlap. The assumed quantity is the overlap squared, 0.1. The proposed sentence says "overlap of 0.1". Read as an
  amplitude that is 0.01 squared, and the campaign is 28 yr, not 4.2. The value has no source. See section 7 for why
  0.1 is favorable for a localized determinant in a periodic box.
- Gaps. 4.44 (12C, 2+), 6.05 (16O, 0+), 1.46 MeV (40Ar, 2+) are the correct physical first excited states. They are
  proxies; the model's spectrum is not computed. For open-shell 12C and 40Ar the 2+ is the right state to use, since a
  single determinant is not a J = 0 state. They are not the gap that the filter sees; see section 7.
- Trotter step. dt = 0.002 sqrt(3/A) keeps the triton's measured prefactor and uses the A scaling of the paper's
  bound (alpha_2 is linear in A). The rigorous bound needs 54 times more steps at every A. The triton prefactor was
  measured on 6^3 at q = 148 MeV for one observable. Three times more steps puts 12C at 1.1e9 T and 12.7 yr.
- Filter. Rodeo with Gaussian times, sigma = 3/gap, 14 cycles, total time 33.5/gap, residual excited weight 1e-3.
  The cycle factor [1 + exp(-dE^2 sigma^2/2)]/2 is the published one. The constant is not favorable: the Lin-Tong
  filter at the same leakage corresponds to 7.4/gap. The hidden input is E_0. An error of 0.1 gap (0.44 MeV for 12C)
  cuts the ground-state acceptance to 0.73 over 14 cycles; 0.2 gap gives 0.30. The filter also projects onto the
  eigenstate of the Trotterized Hamiltonian, so E_0 must be that operator's. The model has a sign problem for
  C_t != C_s with a repulsive D, so a classical E_0 at the 0.1% level is not assured. A quantum energy estimate is not counted.

## 6. Claim (e)

`v_table5.py`, written from the printed Eqs. (89), (90), (107)-(110), (130)-(133), (137), (138), (143)-(145), (217),
(222)-(225). Ratio to the published Table 5:

| method | A | T at t_cross | qubits | T at t_dw | qubits |
|---|---|---|---|---|---|
| PF2 | 16 | 1.443 | 207 / 196 | 1.383 | 201 / 190 |
| PF2 | 40 | 1.070 | 427 / 416 | 1.055 | 423 / 410 |
| QSP | 16 | 0.997 | 234 / 235 | 0.988 | 233 / 234 |
| QSP | 40 | 1.002 | 499 / 500 | 1.000 | 497 / 498 |
| GQSP | 16 | 0.998 | 235 / 235 | 0.991 | 234 / 234 |
| GQSP | 40 | 1.000 | 500 / 500 | 0.997 | 498 / 498 |

Same as the study to the digits shown. The potential term alone gives 0.96 to 0.99 of the published product-formula
totals, so the published kinetic part is smaller than Eq. (19) gives. The claim is true for the qubitization rows only.

## 7. What a referee would attack and the README does not disclose

1. Centre-of-mass motion (main point). A Slater determinant of localized orbitals is not a state of zero total
   momentum. The Hamiltonian's lowest excitation above the zero-momentum ground state is the nucleus moving with
   P = 2 pi/L = 111 MeV: P^2/(2 A m) = 0.54 MeV (12C), 0.41 MeV (16O), 0.16 MeV (40Ar). A filter of width gap/3 = 1.5 MeV
   does not remove these. Two consequences.
   - If the target is the chapter's zero-momentum ground state, the gap is 0.54 MeV, and the zero-momentum weight of a
     localized determinant is of order one over the number of distinct centre positions, 1e-2 or smaller. At a
     0.54 MeV gap the 12C count is 2.9e9 T per shot and 35 yr at overlap squared 0.1, and 234 yr at 0.01.
   - If the moving copies are kept, the measured correlator is an average over total momenta. The recoil phase
     differs by P q tau / (A m), about 0.2 to 0.4 rad at tau = 0.20, against an evolution error target of 0.007.
   The triton box avoids this by loading the zero-momentum state exactly. The README does not mention it.
2. "Overlap" against "overlap squared", as above.
3. "The same triton". The 48 LQ and 2.6e7 T are for 8^3 (E_0 = -6.45 MeV, q = 111 MeV); the box in the chapter is 6^3
   (E_0 = -8.48 MeV, q = 148 MeV, D fitted there). On 6^3 the study has 48 LQ and 2.7e7 T with a bound, not a circuit,
   for the mod-6 Fourier step. "Nothing assumed" in the README table is too strong for 8^3: only E_0 was computed
   there (`triton_results.json`), not the correlator or the bias at 100 steps. The count is insensitive, since the
   evolution is 1.5% of it.
4. Signal size. With J = (1/Z) sum_p e^{i q x_p} at q = 111 MeV the normalized 12C correlator is mostly the elastic
   pole; the inelastic strength is of order (1 - F^2)/Z, roughly 0.1 (my estimate, not computed). The chapter's
   tolerance of +-0.02 is 9% of the triton's elastic-subtracted signal and would be 20% or more here. There is no
   classical elastic subtraction for 12C.
5. The model. The couplings give -136.8 MeV for three nucleons on a site and -19.6 MeV for four. Whether 12C is bound
   against tritons or alphas on 8^3 is unknown. The README says this; it should travel with any 12C number.
6. Verifiability. No exact reference exists for A >= 12, so a 12C row is not a benchmark in the chapter's sense
   (disclosed in the README).
7. The 40Ar count through GQSP rests on a spin-SWAP SELECT that is a sketch. My floor with the paper's SU(4) SELECT
   alone and no spin term is 3.0e9 T for the filter at overlap squared 0.1 (1.7e9 at 0.9); with the sketch 5.5e9,
   against the study's 5.4e9. The filter falls below 1e9 T only for a gap above 4.4 MeV (floor) or 8 MeV (sketch).
8. Minor: the triton's 62 Toffolis against 48 qubits (section 1); one insertion counted instead of two; 13 extra
   potential layers in the filter from the per-cycle n + 1 rule. All below 0.1%.

## 8. Wording I would allow

(a) "On an $8^3$ lattice with the same couplings, a first-quantized layout of this calculation needs about 48 LQ and
$2.6\times10^7$ T per shot, nearly all in the $2^{20}$-rotation loader, and the same campaign takes 49 d."
Not "the same triton", unless the 6^3 numbers (48 LQ, 2.7e7 T, Fourier step bounded and not constructed) are quoted.

(b) "A scoping count for $^{12}$C on $8^3$ gives 160 LQ, $3.6\times10^8$ T per shot, a first result in 66 d and a
4.2 yr campaign. It assumes a squared mean-field overlap of 0.1, a 4.44 MeV gap, and a Trotter step scaled from the
triton's as $\sqrt{3/A}$. None of these is computed for this model, and the count does not treat the center-of-mass
motion of the trial state." If the last clause cannot be carried, I would not print the 4.2 yr.

(c) "For $^{40}$Ar the evolution is about $2\times10^8$ T, but filtering to the ground state at a 1.46 MeV gap
exceeds $10^9$ T per shot ($5.6\times10^9$ to $1.8\times10^{10}$)."

(d) "On antisymmetric states the spin-dependent contact is a sum of mutually commuting terms, so a potential layer is
exact." The Toffoli and rotation counts are a construction-level estimate and should be called that.

(e) "Our implementation of the published formulas reproduces the QSP and GQSP rows of Table 5 to about 1%." Not the
product-formula rows.
