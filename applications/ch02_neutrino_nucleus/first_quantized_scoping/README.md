# First-quantized nuclei with A >= 12 as the Ch. 2 2033 benchmark: scoping study

Prepared 2026-11-12 for the accepted card on the 2033 instance of `app01_neutrino_nucleus.tex`.
Nothing in the paper, in `estimates/` or in the tests was changed. Everything here is new and lives in this folder.

Question. Should the 2033 benchmark be a first-quantized nucleus with A >= 12 (12C, 16O or 40Ar on an 8^3 lattice)
instead of the occupation-mode triton on 6^3?

## 0. Result

Short answer: not as a replacement. A first-quantized 12C on 8^3 fits the register and the per-shot T target, and its
campaign is within or near the five-year window, but only under three inputs that nobody has computed. 40Ar does not fit.
The first-quantized triton is far smaller than the present one and leaves the 1000-qubit class entirely.

| instance, method | LQ | T per shot (deepest circuit) | first result | campaign | fits 1000 LQ | fits 1e9 T | fits 5 yr | the verdict hinges on |
|---|---|---|---|---|---|---|---|---|
| 3H 6^3, occupation modes (current box) | 908 | 2.07e8 | 9.4 d | 250 d | yes | yes | yes | nothing assumed; 864 of the 908 qubits are an encoding choice |
| 3H 8^3, first quantized, PF2 | 48 | 2.61e7 | 2.1 d | 49.4 d | yes | yes | yes | nothing assumed; the 2^20-rotation loader is 98% of the count |
| 3H 6^3, first quantized, PF2 | 48 | 2.68e7 | 2.1 d | 50.2 d | yes | yes | yes | same; the mod-6 Fourier step is a bound, not a circuit |
| 12C 8^3, PF2, calibrated step | 160 | 3.56e8 | 65.6 d | 4.2 yr | yes | yes | yes | the Trotter step (triton-calibrated; the rigorous bound gives 2.01e10 T). Then overlap^2 >= 0.08 |
| 12C 8^3, GQSP | 223 | 1.85e8 | 134.6 d | 8.6 yr | yes | yes | near | overlap^2 >= 0.16 at a 4.44 MeV gap |
| 16O 8^3, PF2, calibrated step | 208 | 5.05e8 | 92.8 d | 6.0 yr | yes | yes | near | the Trotter step, then overlap^2 >= 0.13 |
| 16O 8^3, GQSP | 275 | 2.41e8 | 169.3 d | 10.9 yr | yes | yes | near | overlap^2 >= 0.20 at a 6.05 MeV gap |
| 40Ar 8^3, PF2, calibrated step | 496 | 1.77e10 | 9.0 yr | 211.4 yr | yes | no | no | the 1.46 MeV gap; no overlap rescues it |
| 40Ar 8^3, GQSP | 589 | 5.60e9 | 12.0 yr | 281.5 yr | yes | no | no | the 1.46 MeV gap; needs gap >= 8 MeV and overlap^2 >= 0.9 together |

Baseline for the A >= 12 rows: chapter Hamiltonian with C_t != C_s and the three-body contact, species-labeled registers,
mean-field overlap squared 0.1, gap equal to the physical first excited state, charge correlator at the chapter's three
times and tolerances, serial shots on one machine at 1 us per T. "T per shot" is the deepest complete circuit: one
successful preparation attempt, insertion and controlled evolution to tau = 0.20 MeV^-1. Wall times also count the
rejected preparation attempts. "near" means within 3x of the target, as in the report.

The three assumptions that matter most:

1. **Mean-field overlap (assumed 0.1, plausible 1e-3 to 0.5).** The campaign wall scales roughly as 1/overlap^2. 12C
   needs overlap^2 >= 0.08 (product formula) or 0.16 (GQSP) for five years. At 0.01 every A >= 12 campaign takes
   decades. No value for this lattice contact model exists.
2. **Product-formula step (assumed, calibrated on the triton).** The step is the chapter's measured triton step scaled
   by sqrt(3/A). The rigorous bound of the source paper needs 54 times more steps, which puts 12C at 2.0e10 T per
   shot. The GQSP rows do not use this assumption.
3. **Gap of the finite-volume model (assumed equal to the physical first excited state).** The filter length scales as
   1/gap. Leading-order pionless EFT finds no 16O that is stable against breakup into four alphas (arXiv:1701.06516),
   so the gap in this model may be set by cluster thresholds in an 11.2 fm box and could be near 1 to 2 MeV. At 1 MeV
   the 12C campaign takes 19 yr (product formula) or 38 yr (GQSP).

What could not be settled: whether the contact model binds 12C or 16O at all with the triton couplings; the overlap;
the Trotter error for A >= 12; a compiled circuit for any of the new pieces; T-depth. There is also no exact classical
reference for A >= 12 on 8^3, so the benchmark would lose the check against "the exact finite-volume model" that both
present boxes rely on.

Recommendation: **present both.** Keep the triton as the validated 2033 box. State next to the 908 that a
first-quantized layout of the same calculation needs 48 logical qubits and 2.6e7 T on 8^3 (or on 6^3), so the register
is visibly a choice. Give 12C on 8^3 as a scoped extension in the post-2033 paragraph, with its assumptions named.
A 12C box would carry: 160 LQ (120 system, 40 workspace); 3.6e8 T per shot, of which the filter is 3.5e8 and the
evolution 9.2e6; 8.1e8 T expected per accepted shot; first result 7,012 accepted shots, 5.7e12 T, 66 d; campaign
164,424 accepted shots, 1.3e14 T, 4.2 yr; eps_l <= 2.8e-10. The GQSP variant: 223 LQ, 1.9e8 T per shot, 1.7e9 T
expected per accepted shot, 135 d, 8.6 yr, eps_l <= 5.4e-10.

A side finding for line 221 of the chapter. The "about 500 qubits and 3e7 T" for argon is Table 5 of arXiv:2507.22814:
4 T per Toffoli, t = 0.062 MeV^-1, error 0.1, SU(4) contacts. In the report convention the same point is 5.8e7 T. At
the chapter's tau = 0.20 with spin-dependent contacts the evolution alone is 2.1e8 T on 589 LQ.

## 1. Register

System register. Eq. (5) of arXiv:2507.22814 gives (3m+2) A qubits: 132, 176 and 440 for 12C, 16O and 40Ar at m = 3.
This study uses a species-labeled register instead: protons sit in registers 0..Z-1 and neutrons in Z..A-1, each species
is antisymmetrized on its own, and the isospin qubit is dropped. That gives (3m+1) A = 120, 160 and 400. The chapter
Hamiltonian never changes a particle's isospin label in this form (section 3(i)), and `spin_contact_check.py` confirms
that the labeled spectrum equals the one-proton sector of the fully antisymmetrized one (24 states, difference 1e-13).
The reason to do this is the observable, not the 2 A qubits (section 3(iv)).

| instance | method | system | workspace | total | share of 1000 |
|---|---|---|---|---|---|
| 3H 8^3 | PF2 | 30 | 18: 3 pair flags, 8 equality ladder, 5 AND chain, 1 synthesis, 1 Hadamard | 48 | 5% |
| 3H 8^3 | GQSP | 30 | 74: 65 block encoding and phase gradients (Eq. (132)), 9 second selection | 104 | 10% |
| 12C 8^3 | PF2 | 120 | 40: 24 occupancy counters, 8 ladder, 5 AND chain, 3 synthesis, Hadamard, filter | 160 | 16% |
| 12C 8^3 | GQSP | 120 | 103: 85 block encoding and phase gradients, 18 second selection register and flags | 223 | 22% |
| 16O 8^3 | PF2 | 160 | 48 | 208 | 21% |
| 16O 8^3 | GQSP | 160 | 115 | 275 | 28% |
| 40Ar 8^3 | PF2 | 400 | 96 | 496 | 50% |
| 40Ar 8^3 | GQSP | 400 | 189 | 589 | 59% |
| 12C 6^3 | PF2 | 120 | 40 | 160 | 16% |

The Slater-determinant step needs Z-1 and N-1 ancillas (arXiv:2509.07279); these reuse the counter qubits. With the
Eq. (5) register the totals are 235, 291 and 629 (GQSP). No instance comes near 1000 logical qubits. 12C on 6^3 has the
same register as on 8^3, because 6 sites also need 3 bits. It is not informative beyond that: the Fourier step mod 6
has no power-of-two circuit and is counted here only as a bound (63% more T per step).

## 2. Evolution cost from the published lemmas

### 2a. Check: Table 5 of arXiv:2507.22814

`spagnoli_lemmas.py` implements Lemmas 1 and 3, Theorems 2 to 5, Appendix A and the norms of Appendix C in the paper's
own convention. Ratios to the published value are in parentheses.

| method | A | T at t_cross: ours / paper | LQ: ours / paper | T at t_dw: ours / paper | LQ: ours / paper |
|---|---|---|---|---|---|
| PF2 | 16 | 1.88e+09 / 1.3e+09 (1.44) | 207 / 196 | 1.13e+08 / 8.2e+07 (1.38) | 201 / 190 |
| PF2 | 40 | 3.53e+10 / 3.3e+10 (1.07) | 427 / 416 | 2.22e+09 / 2.1e+09 (1.05) | 423 / 410 |
| QSP | 16 | 6.98e+07 / 7e+07 (1.00) | 234 / 235 | 1.11e+07 / 1.12e+07 (0.99) | 233 / 234 |
| QSP | 40 | 4.01e+08 / 4e+08 (1.00) | 499 / 500 | 6.4e+07 / 6.4e+07 (1.00) | 497 / 498 |
| GQSP | 16 | 3.74e+07 / 3.75e+07 (1.00) | 235 / 235 | 5.93e+06 / 5.99e+06 (0.99) | 234 / 234 |
| GQSP | 40 | 2.11e+08 / 2.11e+08 (1.00) | 500 / 500 | 3.37e+07 / 3.38e+07 (1.00) | 498 / 498 |

- QSP and GQSP: T counts agree to 1%. GQSP qubit counts agree exactly. QSP is one qubit low because the paper prints
  the GQSP count for both.
- Product formula: the potential term alone gives 0.96 and 0.98 of the published totals at
  t_dw (A = 16, 40), so the published kinetic contribution must be small. With T_DIAG = 8 b (2m-1)/2 as printed in
  Eq. (19) the totals come out 5% to 44% high. With T_DIAG = 4 b w_H and w_H = 1, the smallest value Lemma 6 allows,
  the ratios are 1.07 and 1.00 at t_dw and 1.10 and 1.01 at t_cross. The
  paper does not say which was used. Our product-formula qubit count is 11 to 13 above the published one and we could
  not find the cause.
- The published qubit counts use 3m A system qubits, without the spin and isospin qubits of Eq. (5). That is the
  410 against 440 puzzle noted in `ch2_M1_encoding_lit.md`.

The kinetic discrepancy does not enter the results below. This study does not use the phase-kickback kinetic step.

### 2b. Published SU(4) lemmas in the report convention

Evolution only, uncontrolled, SU(4) couplings of the paper, tau = 0.20 MeV^-1, error 0.00731. Toffoli-derived T counts
are multiplied by 7/4 and each rotation costs c_T = 1.15 log2(1/eps_rot) + 9.2 with eps_rot = sqrt(1e-2/N_rot). The
product-formula step count is the paper's bound, Eq. (54) with alpha_2 of Eq. (217). Register with the Eq. (5) system.

| A | PF2 steps (bound) | PF2 T/step | PF2 T | PF2 LQ | GQSP queries | GQSP T/query | GQSP T | GQSP LQ |
|---|---|---|---|---|---|---|---|---|
| 12 | 23132 | 8.51e+04 | 1.97e+09 | 199 | 3461 | 5714 | 1.98e+07 | 216 |
| 16 | 26711 | 1.61e+05 | 4.29e+09 | 245 | 4607 | 7228 | 3.33e+07 | 268 |
| 40 | 42233 | 1.80e+06 | 7.61e+10 | 513 | 11483 | 1.64e+04 | 1.89e+08 | 581 |
  paper point A=16 t=0.0618 eps=0.1: PF2 1.85e+08 (paper conv. 1.13e+08), GQSP 1.02e+07 (paper conv. 5.93e+06)
  paper point A=40 t=0.0624 eps=0.1: PF2 3.49e+09 (paper conv. 2.22e+09), GQSP 5.84e+07 (paper conv. 3.37e+07)

The bound-based product formula exceeds 1e9 T already for 12C. GQSP is 2e7 to 2e8 T.

## 3. What the paper leaves out and the chapter needs

### (i) Spin-dependent two-body contact, C_t != C_s

Appendix E of the paper writes the pair term as delta_ij (C_0 + C_1 sigma_i.sigma_j), notes that these terms no longer
commute, and leaves the count to future work. With Cbar = (C_t+C_s)/2 = -87.911 MeV and dC = (C_t-C_s)/2 = -5.715 MeV
the pair term is delta_ij (Cbar + dC SWAP_ij), where SWAP exchanges the two spin qubits.

Construction. On an allowed state a site holds at most four nucleons. For three or four on one site the pair SWAPs sum
to zero (total spin 1/2 or 0). So the interaction equals

    V = sum_i g(n_i) + dC sum_(i<j) e_ij SWAP_ij,   g(n) = Cbar n/2 + D n(n-1)/6,

where n_i counts the other particles on particle i's site and e_ij = delta_ij [n_i = 1] selects pairs that are alone
on their site. All terms of this form commute, so a potential layer is exact, as in the chapter's occupation circuit.
`spin_contact_check.py` verifies: equality with the non-commuting form on antisymmetric states (1e-14); zero
commutators between the exclusive terms on the full space (the plain pair terms give 0.41); the on-site energies
C_t (x3), C_s (x3), 3 Cbar + D (x4), 6 Cbar + 4 D, identical to the chapter's 22 Pauli terms in
`2033_model_study/triton_resources.json`; and the two-qubit identity
exp(-i th SWAP) = e^(-i th) B^dag diag(1,1,1,e^(2 i th)) B with B = (H x 1) CNOT.

Circuit for one controlled layer, A particles, P = A(A-1)/2 pairs, dm = 9 coordinate bits:

1. For each pair, test equality of the coordinates (dm CNOTs and one C^dm X) and add the flag into the 2-bit counters
   n_i and n_j. Uncompute the flag.
2. Three controlled Walsh phases on each counter apply g(n_i): 6 A rotations.
3. For each pair, recompute the flag, map the spin singlet to |11> with B, take the AND of (flag, n_i = 1, both spin
   bits, Hadamard control) and rotate it: one rotation per pair.
4. Uncompute the counters with a third sweep.

Count: P (6 (dm-1) + 14) Toffolis and P + 6 A rotations per layer, with compute and uncompute both counted at 7 T per
Toffoli as in Lemma 3. This is O(A^2) where Lemma 3 is O(A^3), and it includes the three-body contact. Per layer:
4,092 Toffolis for 12C, 7,440 for 16O, 48,360 for 40Ar. Measurement-based uncomputation would roughly halve the
Toffolis and is not assumed. Status: construction-level count, verified algebra, no compiled circuit.

Trotter error. With exact layers the only error is between T and V. We extend alpha_2 of Eq. (217) to this
interaction (C -> Cbar plus |dC| terms for doubly and triply occupied sites) and use the hopping dispersion's largest
single-particle kinetic energy, 12 t. This is our extension and has not been reviewed.

Block encoding (GQSP route). The paper block-encodes the diagonal part with C = C_0 - C_1 and leaves the SWAP part
open. We use an LCU over pairs: select a second particle with a second unary register and swap network, compare
coordinates, apply the spin SWAP when they match and a sign-averaged identity otherwise. Bounds: one-norm
|dC| A(A-1)/2 (377, 686, 4458 MeV for 12C, 16O, 40Ar), added to A (6 t + (3|Cbar| + 4 D)/2) = 449.3 A MeV; extra
SELECT cost 2(A-1)(dm+2) + 2(dm-1) + 5 Toffolis plus a second PREPARE. The total norm is 5768, 7874 and 22428 MeV. We
did not construct this SELECT gate by gate or verify that it stays self-inverse. The count is a bound by construction
sketch.

### (ii) Three-body contact as the chapter uses it

The chapter's term D sum over on-site triples is the paper's V_3 with G = D. In the product-formula route it is the
D n(n-1)/6 part of g(n) above and adds no gates. In the GQSP route it is inside the paper's potential block encoding
(Lemma 5), norm included above. D = 126.954 MeV was fitted to the triton in the 6^3 box. It is not refitted for
A >= 12 or for 8^3, and nothing here says what binding energies it gives there.

### (iii) Kinetic term

The chapter uses nearest-neighbor hopping, evolved exactly in momentum space. The paper uses a q^2 dispersion by
squaring and phase kickback. We keep the chapter's Hamiltonian. On 8 sites cos(2 pi p/8) has four Walsh terms, so one
controlled kinetic layer costs, per coordinate, two exact QFT_8 (each two controlled-S, 6 T, and one controlled-T,
3 rotations) and 8 rotations: 14 rotations and 12 T. The Fourier transforms are not controlled.

### (iv) Ground-state preparation

The chapter's triton is loaded from a classical wavefunction with no filter and no overlap assumption. That is not
possible for A >= 12. We use:

- Mean-field state. One Slater determinant per species, antisymmetrized with the recursive algorithm of Rule et al.
  (arXiv:2509.07279, Quantum 10, 2056): n(n+1)/2 uses of the orbital unitaries and n(n-1)/2 of their inverses,
  n(n-1)/2 controlled swaps with a C^k X, and n^2 - 4n + 3 rotations, for n = Z and n = N. Each orbital unitary is
  counted as a real-amplitude loader on the 9 spatial qubits, 511 rotations. Cost: 9e5 T (12C), 1.6e6 (16O),
  1.1e7 (40Ar). This is under 1% of every per-shot count. The antisymmetrization of composite states in
  arXiv:2512.16138 was not needed.
- Filter to the correlated ground state. Assumed overlap squared gamma^2 = 0.1 and gap = 4.44, 6.05, 1.46 MeV (the
  physical first excited states of 12C, 16O, 40Ar, used as proxies). Assumed residual excited weight 1e-3. Assumed
  that E_0 is known to much better than the gap from classical lattice Monte Carlo; a quantum energy estimate is not
  counted.
  - Product-formula route: rodeo filter (arXiv:2009.04092) with Gaussian times of rms 3/gap. Each cycle halves the
    excited weight, 14 cycles, total controlled evolution 33.6/gap (7.6 MeV^-1 for 12C, 38 times the response
    evolution). A rejected attempt stops at its first failed cycle, after 2.0 cycles on average. The cycle factor is
    the published one; the schedule is ours and is not optimized.
  - GQSP route: eigenstate filter of Lin and Tong (arXiv:1910.14596, Lemma 2), degree 2 l with
    l = sqrt(2) lambda ln(2/eps_f)/gap. That is 19,276 walk steps for 12C against 2,330 for the evolution. It has one
    flag at the end, so a rejected attempt costs the full filter.
- Attempts per accepted shot: 1/gamma^2 = 10. Amplitude amplification and a staged filter with early rejection would
  lower this and are not counted.

The filter is 97% of the 12C per-shot count in the product-formula route and 89% in the GQSP route. Preparation,
not evolution, decides every A >= 12 verdict.

### (v) Response measurement and shots

Same observable, times and tolerances as the chapter: Re and Im of the normalized charge correlator at tau = 0.05,
0.10, 0.20 MeV^-1, Hadamard test, first result at tau = 0.10 to +-0.05 (7,012 shots), campaign of six estimates at
+-0.02 (164,424 shots). The momentum is q = 2 pi/(8a) = 110.7 MeV. On 8 sites e^(i q x) is a T, an S and a Z gate.

The register layout matters here. With species-labeled registers J = (1/Z) sum_p e^(i q x_p) is an average of Z
unitaries, sampled classically as in the 2028 benchmark, and the shot counts are the chapter's. With the fully
antisymmetrized register of Eq. (5), J needs the projector (1 + tau_z)/2 on every particle. Its LCU norm is A/Z per
insertion, so the same absolute tolerance needs (A/Z)^4 times the shots: 16 for 12C and 16O, 24 for 40Ar. The 12C GQSP
campaign would then take 139 yr. We could not remove this factor inside the Eq. (5) layout.

Shots counted here are accepted shots. The chapter subtracts a classically computed elastic pole. No such reference
exists for A >= 12, so the deliverable would be the full correlator only.

### (vi) Logical error rate

eps_l <= 0.1/N_T at the deepest circuit: 3.8e-9 for the first-quantized triton, 2.8e-10 (12C, product formula),
5.4e-10 (12C, GQSP), 2.0e-10 and 4.2e-10 for 16O, 5.6e-12 and 1.8e-11 for 40Ar. The current triton needs 4.8e-10.

## 4. Per-shot counts and wall times

Report time model: 1 us per T plus 0.1 ms per attempt, one machine, serial shots. No T-depth is derived, as in the
chapter, so every wall time is a lower bound. The qubitized circuits are sequential and may not be able to use ten
factories; the report quotes F* of 1 to 50 for such circuits.

"cal." is the calibrated step, "bound" the rigorous one. "eq5" rows use the fully antisymmetrized register and show
the (A/Z)^4 shot factor. T/step is one second-order step (PF2) or one walk step (GQSP).

| instance | method | encoding | LQ | T/step | T evol (0.20) | T filter/attempt | T deepest circuit | E[T]/accepted shot | campaign shots | first result | campaign | eps_l | fits 1e9 T | fits 5 yr |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 3H 8^3 | PF2 cal. | labeled | 48 | 3878 | 3.89e+05 | 0.00e+00 | 2.61e+07 | 2.61e+07 | 1.64e+05 | 2.1 d | 49.4 d | 3.8e-09 | yes | yes |
| 3H 6^3 | PF2 cal. | labeled | 48 | 1.05e+04 | 1.05e+06 | 0.00e+00 | 2.68e+07 | 2.68e+07 | 1.64e+05 | 2.1 d | 50.2 d | 3.7e-09 | yes | yes |
| 3H 8^3 | PF2 bound | labeled | 48 | 3937 | 2.13e+07 | 0.00e+00 | 4.75e+07 | 4.75e+07 | 1.64e+05 | 3.0 d | 73.5 d | 2.1e-09 | yes | yes |
| 3H 8^3 | GQSP | labeled | 104 | 3539 | 2.01e+06 | 0.00e+00 | 2.77e+07 | 2.77e+07 | 1.64e+05 | 2.2 d | 51.2 d | 3.6e-09 | yes | yes |
| 12C 6^3 | PF2 cal. | labeled | 160 | 7.46e+04 | 1.50e+07 | 5.65e+08 | 5.80e+08 | 1.32e+09 | 1.64e+05 | 106.8 d | 6.9 yr | 1.7e-10 | yes | near |
| 12C 8^3 | PF2 cal. | labeled | 160 | 4.57e+04 | 9.16e+06 | 3.46e+08 | 3.56e+08 | 8.13e+08 | 1.64e+05 | 65.6 d | 4.2 yr | 2.8e-10 | yes | yes |
| 12C 8^3 | PF2 bound | labeled | 160 | 4.78e+04 | 5.17e+08 | 1.96e+10 | 2.01e+10 | 4.56e+10 | 1.64e+05 | 10.1 yr | 236.3 yr | 5.0e-12 | no | no |
| 12C 8^3 | GQSP | labeled | 223 | 8482 | 1.98e+07 | 1.64e+08 | 1.85e+08 | 1.67e+09 | 1.64e+05 | 134.6 d | 8.6 yr | 5.4e-10 | yes | near |
| 16O 8^3 | PF2 cal. | labeled | 208 | 7.57e+04 | 1.75e+07 | 4.86e+08 | 5.05e+08 | 1.15e+09 | 1.64e+05 | 92.8 d | 6.0 yr | 2.0e-10 | yes | near |
| 16O 8^3 | PF2 bound | labeled | 208 | 7.86e+04 | 9.84e+08 | 2.73e+10 | 2.83e+10 | 6.38e+10 | 1.64e+05 | 14.1 yr | 330.4 yr | 3.5e-12 | no | no |
| 16O 8^3 | GQSP | labeled | 275 | 1.06e+04 | 3.37e+07 | 2.05e+08 | 2.41e+08 | 2.10e+09 | 1.64e+05 | 169.3 d | 10.9 yr | 4.2e-10 | yes | near |
| 40Ar 8^3 | PF2 cal. | labeled | 496 | 4.17e+05 | 1.53e+08 | 1.75e+10 | 1.77e+10 | 4.06e+10 | 1.64e+05 | 9.0 yr | 211.4 yr | 5.6e-12 | no | no |
| 40Ar 8^3 | PF2 bound | labeled | 496 | 4.26e+05 | 8.42e+09 | 9.68e+11 | 9.77e+11 | 2.24e+12 | 1.64e+05 | 495.9 yr | 11631.9 yr | 1.0e-13 | no | no |
| 40Ar 8^3 | GQSP | labeled | 589 | 2.36e+04 | 2.12e+08 | 5.38e+09 | 5.60e+09 | 5.41e+10 | 1.64e+05 | 12.0 yr | 281.5 yr | 1.8e-11 | no | no |
| 12C 8^3 | GQSP | eq5 | 235 | 8483 | 1.98e+07 | 1.64e+08 | 1.86e+08 | 1.68e+09 | 2.63e+06 | 5.9 yr | 139.1 yr | 5.4e-10 | yes | no |
| 16O 8^3 | GQSP | eq5 | 291 | 1.06e+04 | 3.37e+07 | 2.05e+08 | 2.42e+08 | 2.12e+09 | 2.63e+06 | 7.5 yr | 175.5 yr | 4.1e-10 | yes | no |
| 40Ar 8^3 | GQSP | eq5 | 629 | 2.36e+04 | 2.12e+08 | 5.38e+09 | 5.61e+09 | 5.42e+10 | 4.01e+06 | 293.3 yr | 6879.2 yr | 1.8e-11 | no | no |
| 3H 6^3 occupation (current) | PF2 | charge | 908 | 1.7e6 | 1.8e8 | none | 2.07e8 | 2.07e8 | 1.64e+05 | 9.4 d | 250.0 d | 4.8e-10 | yes | yes |

For the first-quantized triton the loader is 2.57e7 of the 2.61e7 T. One evolution step costs 3.9e3 T (136 rotations,
62 Toffolis, 108 T) against 1.7e6 T in the occupation circuit. Within 1e9 T the first-quantized triton could run
2.5e5 steps, tau = 500 MeV^-1.

## 5. Sensitivity of the verdict

Campaign wall in years, with the per-shot T of the deepest circuit in parentheses. Rows are the assumed overlap
squared, columns the assumed gap in MeV. The third column is the baseline gap.

12C PF2 (calibrated steps)
| overlap^2 | gap 1.0 | gap 2.0 | gap 4.44 | gap 8.0 |
|---|---|---|---|---|
| 0.9 | 4.2 (7.84e+08) | 2.1 (3.94e+08) | 0.95 (1.81e+08) | 0.54 (1.05e+08) |
| 0.5 | 7.6 (1.23e+09) | 3.8 (6.17e+08) | 1.7 (2.81e+08) | 0.96 (1.59e+08) |
| 0.2 | 12 (1.46e+09) | 6.1 (7.28e+08) | 2.8 (3.31e+08) | 1.5 (1.87e+08) |
| 0.1 | 19 (1.57e+09) | 9.4 (7.84e+08) | 4.2 (3.56e+08) | 2.4 (2.01e+08) |
| 0.01 | 1.3e+02 (1.91e+09) | 63 (9.52e+08) | 29 (4.31e+08) | 16 (2.42e+08) |
| 0.001 | 1.2e+03 (2.36e+09) | 6e+02 (1.18e+09) | 2.7e+02 (5.31e+08) | 1.5e+02 (2.97e+08) |

12C GQSP
| overlap^2 | gap 1.0 | gap 2.0 | gap 4.44 | gap 8.0 |
|---|---|---|---|---|
| 0.9 | 2.5 (4.44e+08) | 1.3 (2.32e+08) | 0.62 (1.16e+08) | 0.37 (7.34e+07) |
| 0.5 | 6.1 (5.97e+08) | 3.1 (3.09e+08) | 1.4 (1.50e+08) | 0.82 (9.25e+07) |
| 0.2 | 18 (6.94e+08) | 8.8 (3.57e+08) | 4 (1.72e+08) | 2.3 (1.05e+08) |
| 0.1 | 38 (7.50e+08) | 19 (3.85e+08) | 8.6 (1.85e+08) | 4.8 (1.12e+08) |
| 0.01 | 4.7e+02 (9.17e+08) | 2.3e+02 (4.69e+08) | 1.1e+02 (2.22e+08) | 59 (1.32e+08) |
| 0.001 | 5.5e+03 (1.08e+09) | 2.8e+03 (5.49e+08) | 1.2e+03 (2.58e+08) | 6.9e+02 (1.52e+08) |

16O PF2 (calibrated steps)
| overlap^2 | gap 1.0 | gap 2.0 | gap 6.05 | gap 8.0 |
|---|---|---|---|---|
| 0.9 | 8 (1.50e+09) | 4 (7.57e+08) | 1.4 (2.61e+08) | 1 (2.01e+08) |
| 0.5 | 15 (2.36e+09) | 7.3 (1.18e+09) | 2.4 (4.00e+08) | 1.8 (3.06e+08) |
| 0.2 | 24 (2.79e+09) | 12 (1.40e+09) | 3.9 (4.70e+08) | 3 (3.59e+08) |
| 0.1 | 36 (3.01e+09) | 18 (1.50e+09) | 6 (5.05e+08) | 4.5 (3.86e+08) |
| 0.01 | 2.4e+02 (3.66e+09) | 1.2e+02 (1.83e+09) | 40 (6.10e+08) | 31 (4.65e+08) |
| 0.001 | 2.3e+03 (4.52e+09) | 1.1e+03 (2.26e+09) | 3.8e+02 (7.51e+08) | 2.9e+02 (5.71e+08) |

16O GQSP
| overlap^2 | gap 1.0 | gap 2.0 | gap 6.05 | gap 8.0 |
|---|---|---|---|---|
| 0.9 | 4.3 (7.58e+08) | 2.2 (3.97e+08) | 0.8 (1.54e+08) | 0.63 (1.25e+08) |
| 0.5 | 10 (1.02e+09) | 5.2 (5.27e+08) | 1.8 (1.98e+08) | 1.4 (1.58e+08) |
| 0.2 | 30 (1.18e+09) | 15 (6.09e+08) | 5.1 (2.25e+08) | 3.9 (1.78e+08) |
| 0.1 | 65 (1.28e+09) | 33 (6.57e+08) | 11 (2.41e+08) | 8.3 (1.90e+08) |
| 0.01 | 8e+02 (1.57e+09) | 4e+02 (8.00e+08) | 1.3e+02 (2.88e+08) | 1e+02 (2.26e+08) |
| 0.001 | 9.4e+03 (1.84e+09) | 4.7e+03 (9.37e+08) | 1.6e+03 (3.33e+08) | 1.2e+03 (2.60e+08) |

40Ar PF2 (calibrated steps)
| overlap^2 | gap 1.0 | gap 2.0 | gap 1.46 | gap 8.0 |
|---|---|---|---|---|
| 0.9 | 69 (1.30e+10) | 35 (6.54e+09) | 48 (8.91e+09) | 9 (1.74e+09) |
| 0.5 | 1.2e+02 (2.03e+10) | 62 (1.02e+10) | 86 (1.39e+10) | 16 (2.65e+09) |
| 0.2 | 2e+02 (2.40e+10) | 1e+02 (1.20e+10) | 1.4e+02 (1.65e+10) | 26 (3.11e+09) |
| 0.1 | 3.1e+02 (2.58e+10) | 1.5e+02 (1.30e+10) | 2.1e+02 (1.77e+10) | 39 (3.34e+09) |
| 0.01 | 2.1e+03 (3.14e+10) | 1e+03 (1.57e+10) | 1.4e+03 (2.15e+10) | 2.6e+02 (4.02e+09) |
| 0.001 | 2e+04 (3.88e+10) | 9.8e+03 (1.94e+10) | 1.3e+04 (2.65e+10) | 2.5e+03 (4.94e+09) |

40Ar GQSP
| overlap^2 | gap 1.0 | gap 2.0 | gap 1.46 | gap 8.0 |
|---|---|---|---|---|
| 0.9 | 27 (4.79e+09) | 14 (2.50e+09) | 19 (3.35e+09) | 4 (7.92e+08) |
| 0.5 | 65 (6.43e+09) | 33 (3.33e+09) | 45 (4.48e+09) | 8.8 (9.98e+08) |
| 0.2 | 1.9e+02 (7.47e+09) | 95 (3.85e+09) | 1.3e+02 (5.19e+09) | 24 (1.13e+09) |
| 0.1 | 4.1e+02 (8.08e+09) | 2.1e+02 (4.15e+09) | 2.8e+02 (5.60e+09) | 52 (1.20e+09) |
| 0.01 | 5e+03 (9.88e+09) | 2.5e+03 (5.05e+09) | 3.5e+03 (6.83e+09) | 6.3e+02 (1.43e+09) |
| 0.001 | 5.9e+04 (1.16e+10) | 3e+04 (5.91e+09) | 4.1e+04 (8.02e+09) | 7.5e+03 (1.64e+09) |

Overlap needed for a five-year campaign at the baseline gap: 12C 0.08 (PF2, calibrated) or 0.16 (GQSP); 16O 0.13 or
0.20; 40Ar not reachable with PF2, 0.999 with GQSP. With the rigorous step bound no A >= 12 product-formula case fits
at any overlap.

Is a first-quantized triton a better 2033 benchmark than the occupation-mode one? It is the same physics at 48 instead
of 908 logical qubits and 2.6e7 instead of 2.07e8 T, with a 49 d campaign instead of 250 d. It is fully checkable and
assumes nothing. But it is a 50-qubit calculation, so it no longer represents a 1000-qubit machine; it needs the 2033
machine only for its T count. The 8^3 box is the natural one (power-of-two Fourier transform, and the chapter already
has E_0 = -6.4475 MeV there). Its cost is set by the loader, which grows as the square of the volume: a 16^3 box would
need 2^26 rotations, about 1.8e9 T. So the first-quantized triton cannot grow in volume by loading. It can grow in
evolution time.

## 6. Inputs and provenance

| input | value | tag | source | note, range |
|---|---|---|---|---|
| t_hop, C_t, C_s, D | 10.5794, -93.62604, -82.19588, 126.95409 MeV | Cited | app01_neutrino_nucleus.tex 2033 model; TRITON.md | Chapter couplings, fitted for A<=3 in the 6^3 box. Not refitted for A>=12 or for 8^3. |
| lattice spacing | 1.4 fm | Cited | app01; arXiv:1911.06368 |  |
| time grid | 0.05, 0.10, 0.20 MeV^-1 | Cited | app01 2033 box | Kept. q = 2 pi/(8a) = 110.7 MeV on 8^3. |
| shots | 7012 first result; 27404 x 6 campaign | Cited | app01 Eq. (Nshot) | Hoeffding plus union bound; +-0.05 (2 estimates) and +-0.02 (6 estimates), 95% simultaneous. |
| system register | (3m+2) A qubits, m=3 | Cited | arXiv:2507.22814 Eq. (5) | Labeled triton: (3m+1) A = 30. |
| Lemma 1, Lemma 3, Thm 2-5, App. A, Eqs. (90), (108), (217) | formulas | Cited | arXiv:2507.22814 | Implemented in spagnoli_lemmas.py; Table 5 reproduced (README section 2). |
| 7 T/Toffoli, RUS synthesis, eps_syn=1e-2 | common.py | Cited | main-overview.tex conventions; estimates/common.py |  |
| 1 us per T, 0.1 ms per shot, 0.1 faults per shot | report | Cited | main-overview.tex | Serial T count on one machine; no T-depth derived. |
| spin-dependent contact, exclusive-pair form | V_excl | Computed | spin_contact_check.py | Equal to the chapter's on-site interaction on antisymmetric states; all terms commute. Circuit is a construction-level count, not compiled. |
| potential layer counts | P(6(dm-1)+14) Toffoli, P+6A rotations | Computed | scoping.py potential_layer | Three equality sweeps; compute and uncompute both counted at 7 T per Toffoli. |
| kinetic layer counts (8^3) | 14 rotations + 12 T per coordinate | Computed | scoping.py kinetic_layer; Walsh count in spin_contact_check.py | Exact QFT_8 (2 CS + 1 CT each way) and 4 controlled Walsh phases for the hopping dispersion. |
| kinetic layer counts (6^3) | <= 40 rotations + 16 Toffoli per coordinate | Computed (bound) | scoping.py kinetic_layer | CRT relabeling, generic two-qubit QFT_3; no circuit constructed. |
| generalized alpha_2 for C_t != C_s | Eq. (217) with C -> Cbar and +|dC| terms | Computed | scoping.py alpha2_chapter | Our extension of the paper's semi-norm bound; hopping dispersion uses max single-particle kinetic energy 12 t. |
| block-encoding norm | lambda = A(6t + (3|Cbar|+4D)/2) + |dC| A(A-1)/2 | Computed (bound) | scoping.py lam_chapter | Paper Eq. (90) form with hopping kinetic LCU and a pair-SWAP LCU. The paper leaves the SWAP block encoding to future work. |
| spin SELECT cost | 2(A-1)(dm+2)+2(dm-1)+5 Toffoli plus a second PREPARE | Computed (bound) | scoping.py walk_step | Second swap network and comparator; construction-level, self-inverse property not verified. |
| Slater-determinant preparation | A^2 orbital unitaries, A(A-1)/2 swap+C^kX, A^2-4A+3 rotations | Cited | arXiv:2509.07279 Sec. 2-3 | Each orbital unitary costed here as a 511-rotation real loader on 9 spatial qubits (Computed). |
| eigenstate filter degree | 2 l, l = ln(2/eps_f) lambda sqrt(2)/gap | Cited | arXiv:1910.14596 Lemma 2 | Needs E_0 known to much better than the gap. That knowledge is Assumed. |
| evolution error target | 0.00731 | Assumed | this study | Half the +-0.02 tolerance; chapter measured 0.0073. T scales as eps^-1/2 (PF) or log (QSP). |
| PF step count, calibrated | 100 sqrt(A/3) (tau/0.2)^1.5 sqrt(0.00731/eps) | Assumed | chapter triton bias, extensive scaling | No A>=12 numerics exist. Range: this value up to the rigorous bound (60-110 times more steps). |
| PF step count, bound | sqrt(4 tau^3 alpha_2/eps) | Cited | arXiv:2507.22814 Lemma 7 | Worst-case semi-norm; 100 times above the triton's measured need. |
| ground-state gap | 4.44 (12C), 6.05 (16O), 1.46 (40Ar) MeV | Assumed | physical first excited states as proxies | The model's finite-volume gap is unknown. Range 1-8 MeV; LO pionless nuclei may sit near cluster thresholds. |
| mean-field overlap squared | 0.1 | Assumed | no source | Range 1e-3 to 0.5. Sets attempts per accepted shot. |
| filter infidelity | 0.001 | Assumed | this study | Enters logarithmically. |
| rodeo filter time | N_cyc x 0.8 x 3/gap | Assumed | arXiv:2009.04092 cycle factor [1+exp(-gap^2 sigma^2/2)]/2, sigma=3/gap | Product-formula route only. Constants from the published cycle factor; schedule not optimized. |
| charge-operator LCU norm | (A/Z)^2 on the correlator, (A/Z)^4 on shots | Computed | scoping.py shot_factor | Antisymmetrized registers carry no proton label. Isoscalar density has factor 1. |
| qubitized-route extra qubits | A + 6 | Assumed | this study | Second unary register and flags; not laid out. |
| current triton numbers | 908 LQ, 2.07e8 T, 9.4 d, 250 d | Cited | app01 2033 box |  |

## 7. Not settled

- No circuit in this study is compiled. The exclusive-pair layer has verified algebra and a gate-level description. The
  spin SELECT of the GQSP route and the mod-6 Fourier step are bounds from sketches.
- Product-formula step for A >= 12. No numerics exist; the calibrated and rigorous step counts differ by 54.
- Whether the triton-fitted contact model binds 12C, 16O or 40Ar on 8^3, and with what gap and mean-field overlap.
  The couplings are not refitted.
- E_0 for the filter is assumed known.
- There is no exact classical reference for A >= 12, and no elastic subtraction.
- T-depth, factory schedule, and whether the qubitized circuits can feed ten factories.
- The kinetic part of the paper's product-formula rows in Table 5 (section 2a).
- Rodeo schedule constants, amplitude amplification, staged filters, Hamming-weight phasing and measurement-based
  uncomputation were not optimized. They could lower the counts by small integer factors.

## 8. Files and reproduction

| file | content |
|---|---|
| `spagnoli_lemmas.py` | arXiv:2507.22814 formulas in the paper's convention; reproduces Table 5 |
| `spin_contact_check.py`, `spin_contact_check.json` | numerical checks of section 3(i) and of the labeled register |
| `scoping.py` | the model in the report convention; imports `estimates/common.py` |
| `results.json` | every number in this README, with the inputs and their tags |
| `run_output.txt` | printed tables of the last run |

    python3 spagnoli_lemmas.py
    OMP_NUM_THREADS=1 python3 spin_contact_check.py
    python3 scoping.py > run_output.txt

Each runs in a few seconds on one thread. The remote host was not needed. Sources read: arXiv:2507.22814 (full text),
2509.07279, 2512.16138, 1910.14596 (Lemma 2), 1711.10460, and the abstracts of 2009.04092, 1701.06516, 1712.10246.
The rodeo cycle factor [1 + exp(-gap^2 sigma^2/2)]/2 is quoted from memory of arXiv:2009.04092 and should be checked
against the paper before it is cited.
