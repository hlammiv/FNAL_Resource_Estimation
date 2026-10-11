"""Ch. 8 — Baryogenesis (bubble nucleation and C-violating wall scattering).
Reproduces the LQ and hard-op numbers of applications/app05_baryogenesis.tex.

The chapter's cost model is Eqs. (Nq_baryo)-(Nshot_baryo), app05:81-83:

    N_q      = V [ log2 N_phi + 2 N_f ] + N_anc
    N_gate   = N_Trot * V * c_T(N_phi, N_f)          per shot
    N_shot   = max( N_rare, C_est eps_C^-2 m_max^-1 N_bin )   (depth-capped MLAE; Heisenberg floor eps_C^-1 N_bin)

ESTIMATOR CONSTANT (R16 author ruling; R15 derivation kept). C_est = 1, every circuit at m_max = 33. The authors
adopt |A_R| < sin(pi/66) = 0.048 as the branch prior for the Delta n arm (R16, 'that value sounds fine'); the
~1e-3 signal estimate sits 48x inside it. C_est is Stated(1, author prior). The R15 analysis below is why a prior
is needed (no bound is derived from the chapter's inputs); Suzuki's LIS constant 1.457 is kept as a sensitivity
(intermediates mlae_lis_*), not priced.
R15 analysis. The Delta n readout encodes the asymmetry as a = (1 + Delta)/2: a
control qubit in |+> selects the particle or antiparticle packet and the good outcome is (particle, reflected)
or (antiparticle, transmitted), so a = [P_R(p) + 1 - P_R(pbar)]/2. Suzuki's Fisher information (1904.10246
Manuscript_v2.tex Eq. Fisher_final, :238) per circuit of m = 2k+1 applications is m^2/(a(1-a)) for a, i.e.
m^2/(1 - Delta^2) for Delta. All N circuits at one depth m: queries N m = (1 - Delta^2)/(eps_C^2 m), constant 1
(the 1 - Delta^2 = 1 - 1e-6 is dropped). That requires a prior that keeps m theta in one monotone branch of
sin^2(m theta), theta = arcsin sqrt(a) = pi/4 + arcsin(Delta)/2. For odd m, m pi/4 is the center of a branch,
so the condition is |Delta| < sin(pi/(2 m)) = 0.0476 at m = 33. No state-independent bound gives this: the
~1e-3 is an order-of-magnitude estimate of the very quantity being measured, and the 'suppressed below
Delta theta_C' statement is itself a scaling estimate (|y| and the wall-profile integral are not fixed); even
taken as a bound, |Delta| < 0.1 is about twice sin(pi/66) = 0.048 (R15b). Without the R16 prior the 2033 schedule
would be Suzuki's LIS (:268, :271): m = 1, 3, ..., m_max, equal shots N per depth. Queries N sum m; Fisher N sum m^2;
at equal variance the query count is C_LIS/(eps_C^2 m_max), C_LIS = m_max sum m / sum m^2 = 33 x 289 / 6545
= 1.45714 (EIS 1, 3, 5, 9, 17, 33: 1.5020). Circuits rise 2.83x (r14 quotes circuit ratios for Ch. 3).

c_T IS DERIVED HERE (R12, item 14; app05:85). It was a stated undecomposed total (5e2 / 1e3 / +5e2).
Per site per Trotter step, n_q = log2 N_phi qubits on the centered field-amplitude grid, phi linear in Z
(arxiv_2210_07985 Eq. phiqub, main.tex:341-346; arxiv_2108_10793 half-integer sampling, bosons.tex:2318):

  scalar (term-level counts from the papers; rotations only, no T given there)
    on-site polynomial a phi + b phi^2 + c phi^3 + d phi^4 (mass, phi^4, the bias of the biased well, the
      gradient's phi^2 share): every Z-string of weight 1..min(4, n_q), one R_Z each
      = C(n_q,1)+C(n_q,2)+C(n_q,3)+C(n_q,4)                  arxiv_2407_13819 Eq. rzTrotter (:1756-1762);
                                                              arxiv_2210_07985 Eqs. trdispl/trphi2/trphi4 (:374-435)
      -> 15 at n_q=4, 3 at n_q=2
    gradient phi_i phi_j: n_q^2 ZZ per link, links_per_site = d = 2 (periodic)  arxiv_2407_13819 E_D (log2k+1)^2 (:1760);
                                                              arxiv_2210_07985 Eq. trphijphik (:409-416)
      -> 32 at n_q=4, 8 at n_q=2
    pi^2 = F (sum ZZ) F^dag on the centered grid: C(n_q,2) ZZ  arxiv_2210_07985 Eq. trpi2 (:397-407)
      -> 6 at n_q=4, 1 at n_q=2. The centered F = D QFT D (Eq. fftqft, :350-356); the inner D's commute
      with the diagonal ZZ and cancel, the outer D's commute with every phi-diagonal term and cancel between
      steps, the first/last merge into the on-site weight-1 rotations (no extra rotations; derived here).
    two exact n_q-qubit QFTs per site per step (derived here; the paper's AQFT formula, arxiv_2407_13819
      Eq. approxQFT_ampTrot :1769, gives ~5e2 T per 4-qubit QFT and is not used). Controlled phase pi/2^k
      between qubits at distance k, (n_q - k) of them:
        k=1  controlled-S = T T CNOT T^dag CNOT        3 T exact
        k=2  controlled-T = Toffoli into ancilla + T  7 + 1 = 8 T exact (measurement uncompute)
        k>=3 Toffoli into ancilla + R_Z(pi/2^k)       7 T + 1 synthesized rotation
      -> per QFT at n_q=4: 3x3 + 2x8 + 1x7 = 32 T + 1 rotation; at n_q=2: 3 T.
    TOTAL scalar per site-step: n_q=4 (N_phi=16): 55 rotations + 64 T;  n_q=2 (N_phi=4): 12 rotations + 6 T.

  two 2-component Wilson flavors (L, R) with Yukawa, r = 1, under Jordan-Wigner (derived here; no paper
  gives a term count for this Hamiltonian). H = sum psi^dag beta (m+2r) psi
      - 1/2 sum_{x,j} [psi_x^dag (r beta - i alpha_j) psi_{x+j} + h.c.]
      + |y| phi [e^{i theta_C} psi_L^dag beta psi_R + h.c.],  beta = sigma_3, alpha_j = sigma_j.
    hopping: at r = 1, A_j = (sigma_3 - i sigma_j)/2 has rank 1, A_j = |v_j><u_j| with u_j orthogonal to v_j.
      Rotating each site's 2 modes to (v_j, u_j) (a pi/4 Givens with Clifford phases = two pi/8 Pauli
      rotations = 2 T; in and out = 4 T per site per flavor per direction) leaves one c^dag c + h.c. per
      link per flavor = (XX + YY)/2 x Z-string: 2 R_Z. Per site: links 2 x flavors 2 x 2 = 8 R_Z;
      d 2 x flavors 2 x 4 T = 16 T. JW strings are parity ladders: CNOTs only.
    mass + Wilson term: (m+2r)(n_1 - n_2) per flavor: 2 R_Z per flavor -> 4.
    Yukawa: per spinor component a, e^{i theta} c_La^dag c_Ra + h.c. = e^{-i theta n_La} (c^dag c + h.c.)
      e^{i theta n_La}; phi (x) (XX+YY)/2 = n_q Z-terms x 2 strings -> 2 n_q R_Z per component, plus the
      theta_C(z) phase frame in and out, 2 R_Z per component. Per site: 2 x 2 n_q + 2 x 2
      -> 16 + 4 = 20 at n_q=4, 8 + 4 = 12 at n_q=2.
    TOTAL fermion per site-step: 32 R_Z + 16 T at N_phi=16; 24 R_Z + 16 T at N_phi=4.
    4-component spinor: A_j rank 2, 4 mass terms, 4 Yukawa bilinears: every fermion count x2 (components/2).

  pricing (report rules): T per rotation = 1.15 log2(1/eps_rot) + 9.2 (common.t_per_rotation "rus"), with
  eps_rot = sqrt(1e-2 / N_rot) and N_rot the synthesized rotations in ONE circuit of that instance
  (common.eps_rot_for, R-TOL); Toffoli = 7 T (R5); c_T = n_rot t_rot + n_T. So c_T differs between
  instances of the same Hamiltonian through t_rot.

  2033 deepest-circuit budget (R13 author ruling): the depth-capped MLAE schedule runs circuits of up to
  m_max base-circuit applications, each with its own eps_syn = 1e-2. The base circuit's rotations are
  synthesized to the deepest circuit's budget, N_rot = m_max x N_rot(base). m_max = floor(1e9 / T_base) and
  T_base depends on m_max through eps_rot, so the pair is iterated to its fixed point from m = 1 (the base
  circuit's own budget): m = 1 -> T 2.6555e7 -> m 37 -> T 2.9683e7 -> m 33 -> T 2.9583e7 -> m 33 (fixed).
  Every start m in [28, 39] maps to 33 and every m <= 27 maps to 34, which maps to 33: the fixed point is unique.

  Instances (site-steps, N_rot, eps_rot, t_rot, c_T, raw T):
    2028 6^2 N_phi=16 scalar      432,   23760, 6.49e-4, 21.379, 1239.8, 5.356e5; /3.3 = 1.623e5
    8^2 N_phi=4 family            768,    9216, 1.04e-3, 20.593,  253.1, 1.944e5; /3.3 = 5.89e4
    8^2 N_phi=16 family           768,   42240, 4.87e-4, 21.856, 1266.1, 9.723e5; /3.3 = 2.95e5
    6^2 N_phi=16 + fermions       432,   37584, 5.16e-4, 21.759, 1973.0, 8.523e5; /3.3 = 2.58e5
    interim 4^2 N_phi=4 + ferm.   192,    6912, 1.20e-3, 20.354,  754.8, 1.449e5; /3.3 = 4.39e4
    2033 base-circuit budget    12000, 1044000, 9.79e-5, 24.517, 2212.9, 2.656e7 (superseded R13; recorded)
    2033 deepest-circuit budget 12000, 33 x 1044000 = 34452000, 1.70e-5, 27.417, 2465.3, 2.958e7 (raw; not banked)

INPUTS AND SOURCES
  V_2028 = 36 (6^2), N_phi = 16 -> 4 qubits/site, N_f = 0, N_anc = 30      app05:122-127
  N_Trot(2028) = 12 steps of dt = 0.1/m_phi to t = 1.2/m_phi (box and prose) app05:82, 91, 129
  banked c_T reduction: ~3.3x printed in the box (R4 convention A); assumed development target, not derived or cited (R17 c)
    (R11); applies to every sector (R11 ch08-2e5-with-fermions A); not banked in 2033 (R11 A)
  V_2033 = 120 (15x8), N_f = 2, 2-component spinor, N_anc = 40 (30 bath + 10 workspace)  app05:156-164, 95
  N_Trot(2033) = t_max/dt = (10/m_phi)/(0.1/m_phi) = 100                     app05:82
  eps_C = 1e-3, N_bin = 100, m_max = floor(1e9/T_base) = 33 (fixed point), C_est/(eps_C^2 m_max) = 3.0e4/bin, C_est = 1 (R16)  app05:83, 88-89, 166-168
  2033 hard-op envelope 1e9, 2028 cap 1e5 (DOE RFI)
  t_gate_s = 1 us per T; shot_overhead_s = 0.1 ms (register init, final readout, decode; R17 g_rate)  app05:130, 151
    E27 (r23): one machine everywhere; no machine count is formed. R25 tiers and walls: see the wall-time section below.
  fault budget 0.1 expected faults per shot (R3)                             app05:131, 167
  interim variant: 4^2, N_phi=4, N_f=2 -> 96 system + 17 workspace = 113, printed ~115  app05:176-177
  The 2028 run drops the fermion arm for resolution (N_phi=4 cannot resolve the biased well; the fermion
  arm enters at N_phi=16, 288 system qubits at 6^2, and waits for 2033); the T numbers are facts, not
  the reason (R11b gap-logic; R12 light touch, app05:91, 181).
WHAT IS NOT DERIVED HERE
  The 3-4x reduction (development target); the box banks 3.3 exactly, printed results rounded once (R4).
  R17 (c): the chapter states it as an assumed development target, neither derived nor cited; 'derived here'
  in the box applies to c_T only.
R17 (e) CHECK: the Fermion_Primitives (FP) hop prices gauge-covariant hops; Ch.8 has no gauge field, so no
  color squish or diagonalizer enters and groups.hop_link_cost does not apply. Mass + Wilson term: 4 R_Z per
  site-step = FP section_mass.tex 2N at d = 2, N = N_c N_f = 2 (groups.wilson_mass_rotations(1, 2, 2)).
  Spinor basis change: FP's 8N T per link (W+ and W- at the two ends) would give 32 T per site-step; Ch.8 uses
  16 T because at r = 1 the hop is rank one with the same vector at both ends, so one pi/4 Givens per site per
  direction (in and out) serves every link of that direction. At FP's count c_T(2033) would rise by 16 T
  (0.6%). Nothing moves.
  Gibbs-state prep (2033), adiabatic ramp (2028), the Lindblad dissipator and the basin projector: carried
  at 0 T, UNSOURCED. Boundary conditions: c_T counts d = 2 links per site (periodic); open boundaries
  remove at most 1/L of the links (e.g. 8 of 240 on the 15x8 cylinder), not taken.
  MLAE queries/bin and m_max cited to arxiv_1904_10246, arxiv_2012_03348; C_est = 1 rests on the author-adopted
  branch prior (R16); the LIS sensitivity constant 1.457 derived here from Suzuki's Fisher information (R15). Since R13 the 2033 rotations are
  synthesized to the deepest MLAE circuit's budget (fixed point above); the 2028 box and the interim
  instance are single circuits and keep their own-circuit budget. Since r25 the 2033 nucleation circuits are
  single circuits too and are synthesized to their own budget (2.656e7 T at t, 8.45e6 T at t/3).
  The f_eq-weighted integrated asymmetry S_int is NOT computed: the first-result tier takes S_int = 1e-3 (the
  chapter's per-k estimate) and its Delta-n wall scales as (1e-3/S_int)^2 (NEEDS_AUTHOR).
REFEREE REVISION (2026-10-04; editorial_review/REFEREE_REPORT.md B2, B3, G2; responses/ch08.md)
  G2: the 2033 headline (Result.hard_ops, PUBLISHED, Table 1.1) is the deepest executed circuit, m_max base
  circuits = 33 x 2.958e7 = 9.76e8 T, printed '>~9.8e8' because state preparation is unpriced. The base circuit
  (2.958e7) stays the per-query unit of every wall time (intermediates t_per_shot, base_circuit_s). Result.shots is
  the number of depth-33 circuits (9.18e4), so shots x hard_ops = queries x base circuit. The landscape row is
  conditional (preparation unpriced).
  B2 (derived here): the Delta-n base circuit is the unitary A = U_evo U_pk U_vac. U_vac is the vacuum of the wall
  background (Eq. deltan as written: f_eq carried classically, as a bin or a label amplitude); the Lindblad bath with
  mid-circuit resets supplies no A^dagger, and a purified Gibbs state (thermofield double) would double the
  960-qubit system register. Each Grover iterate = A^dagger, A, S_chi, S_0. Bounds derived here (not added to the
  count): S_chi = incident-side charge of <= n_modes = 480 occupations counted into a 9-bit register by controlled
  increments (<= 10 Toffolis each, ~20 ancillas; the arxiv_1902_10673 adder tree needs ~476 ancillas the register
  lacks), computed and uncomputed; S_0 = controlled phase on the ~1000-qubit input, lq - 2 Toffolis; together
  7.42e4 T per iterate (0.25% of the base circuit). S_chi is a threshold flag, not a projector onto a reflected
  particle: the Slater-determinant vacuum has no sharp half-system charge, and the C-odd vacuum polarization of
  the wall is a background; contrast and background are not quantified (NEEDS_AUTHOR). Packet: one fermion mode, <= n_modes - 1 Givens per charge, 2 charges,
  2 rotations per Givens, x2 for the selector control: 1.05e5 T. Free-fermion Slater determinant: <= n(n-1)/2
  Givens of two rotations (as Ch. 7): 6.30e6 T (21% of the base circuit). First-result label-controlled packets:
  <= 2^7 x packet = 1.34e7 T. The interacting dressing of U_vac is not priced. r = preparation / evolution.
  Depth cap with preparation: A costs (1 + r) x evolution; m_max = floor(1e9 / T_A) at the rotation budget of
  round(m (1 + r)) base circuits, iterated to its fixed point (prep_sensitivity). Delta-n walls grow ~(1 + r)^2:
  r = 1 -> m_max 16, odd depth 15, 12.5 yr per campaign Delta-n point (2.84 at r = 0).
  The 30 bath ancillas are not used by any Delta-n circuit (the bath serves the nucleation circuits only), so the
  first result's rotation workspace A = 33 needs no borrowing ruling.
OPEN ITEMS B2 / B5 / G4 (2026-10-05; responses/ch08.md 'Open items'; scratch scripts chopen/ch08/)
  B2 preparation, now priced per application of A (the r26 'bounds, not added' are superseded):
    scalar product state, N_phi - 1 = 15 rotations per site (binary-tree amplitude loading about phi_cl(z)): 1800 rot;
    fermion vacuum of the z-dependent wall, block-diagonal in k_y: 8 sectors x (N - eta) eta = 30 x 30 Givens
      (arxiv_1711_05395 Eq. n_gates), 2 rotations each: 14400 rot (was the dense 480*479/2 bound, 6.3e6 T);
    inverse 8-point fermionic FFT on 60 chains, rotation-free (arxiv_1902_10673 :1102): F_2 = 2 T, twiddles <= 1 T,
      32 T per chain, 1920 T;
    selector-controlled packet in one k_y sector, <= 59 Givens per charge: 472 rot.
    Total 4.59e5 T = 1.55% of the evolution (r_free). The rotations join the synthesis budget. A depth-m circuit also
    executes (m - 1)/2 reflection pairs (7.42e4 T each), so the headline is 33 x 3.0056e7 + 16 x 7.42e4 = 9.93e8 T.
    Depth 33 leaves 0.71 Trotter steps of headroom per application: any dressing ramp lowers the depth. The ramp
    (gradient, anharmonicity and Yukawa of the fluctuations) is r x the evolution, unpriced: r = 0.2 -> depth 27,
    9.7 yr per campaign point; r = 1 -> depth 15, 29 yr. First-result circuits carry 2^7 label-controlled packets
    (<= 1.66e6 T), so their A is 3.17e7 and they run at depth 31.
  B2 flag contrast (flag_contrast): particle branch [Q_inc >= 1], antiparticle [Q_inc >= 0] give
    a = 1/2 + delta/2 + p0 A_R/2, p0 = P(vacuum fluctuation q = 0) from the free-vacuum full counting statistics on the
    incident half: 0.571, 0.658, 0.767 at a m_f = 0.25, 0.5, 1 (one-flavor variances 0.32, 0.23, 0.14 reproduce the
    verifier). Queries x 1/p0^2 = 2.31 at a m_f = 0.5 (assumed = a m_phi). delta = 0: C' (c -> W c^dag, W = sigma_1 in
    flavor x sigma_1 in spinor) commutes with H term by term at any real phi (c_prime_residual = 0 on a C-violating
    wall), so the vacuum Q_inc distribution is symmetric. C' maps a flavor-L particle onto a flavor-R antiparticle
    with the same density at all times, so the flavor-summed Delta n vanishes exactly; the readout uses flavor-L
    packets in both branches (A_R = R(pL) - R(pbarL) = R(pL) - R(pR)). Walls: campaign Delta-n point 6.67 yr (was 2.84),
    first-result point 246 d (was 93), first result 1.41 yr (0.58), campaign 70.0 yr (31.7), 14.0x the horizon (6.3),
    scan + first result 4.67 yr (3.84).
  B5 (filtered_jump_t): arxiv_2311_09207 Thm. L_cost: ~beta of controlled Hamiltonian simulation per unit Lindblad
    time on a patch of radius v_LR beta; ||sum A+ A|| <= 1 with adjoints included. beta = 1/m_phi = 10 dt; patch 5x5;
    c_T at the nucleation circuit's own budget (2212.9): 5.53e5 T per unit time. A sweep over the >= 2V = 240 scalar
    jumps: 1.33e8 T = 5.0 nucleation circuits; the 1e9 envelope holds 7.3 sweeps. Lower-end estimate (uncontrolled,
    polylogs dropped); the sweeps to mix are unknown, so Gibbs prep and dissipator stay at 0 T (UNSOURCED).
  G4: Poisson faults, mean lam; p_f = 1 - e^-lam = 0.095 at the R3 budget. a_obs = (1 - p_f) a + p_f q_f, so the A_R
    estimate is scaled by 1 - p_f and shifted by 2 p_f (q_f - 1/2)/(m p0) <= p_f/(m p0) = 4.4e-3. Below the first
    result's sigma = S_int/3 = 3.3e-4 at m = 31: p_f <= 6.8e-3, eps_l <= 6.9e-12 per T; per bin (1e-3, m = 33):
    p_f <= 0.022, eps_l <= 2.2e-11. Idle locations: 1000 LQ x 993 s / 10 us = 9.93e10 qubit-cycles; at 1e-10 they add
    ~10 faults; 0.1 faults over T + idle needs 1.0e-12 per location (6.8e-14 for the first-result bias). No fault is
    heralded; the final occupation readout's total-charge check rejects X/Y faults on fermion qubits only. Null
    control: flavor-R antiparticle branch, a = 1/2 exactly by C'; equal-split difference quadruples the queries.
WALL TIME AND FACTORIES
  t_gate_s = 1e-6 s per T, shot_overhead_s = 1e-4 s per shot, ONE machine, serial (E27).
  Rulings R1/R2 (H. Lamm 2026-10-02) and the shot audit: Gamma/V is the slope of ln P_FV, so every nucleation
  point measures P_FV at two times, t/3 and t, with 3x the printed single-time count at each (2028: 1.5e4 at 4 and
  12 steps; 2033: 3e4 at 33 and 100 steps). Each circuit is priced at its own depth.
  T-depth per shot from factory.json (ch08): 2028 banked 5.4e3-7.0e3 (A = 29-30 RUS ancillas in flight;
  if the 3.3x cuts T-count but not depth, 1.78e4-2.31e4 raw and F* = 7-9, not taken); 2033 per base circuit
  7.2e5 (A = 40, the 30 bath ancillas borrowed during the closed scattering segments) to 2.87e6 (A = 10 workspace).
  F* = N_T / D_T: 2028 23-30, 2033 10.3-41. Both >= 10, so the 10-factory baseline wall stands.
  First-result Delta-n circuits carry the ~7-qubit packet label, which leaves A = 3 workspace qubits:
  t_depth_2033_first 8.7e5 (A = 33, bath borrowed; F* 34, as the chapter states) to 9.6e6 (A = 3; F* 3.1, the
  first result depth-bound at ~1.7 yr instead of ~0.58 yr).
  Mixed-depth runs enter depth_exports as full-circuit equivalents, shots x T_i / T_headline (depth scales with
  T under the ancilla-bound schedule); the walls themselves are summed circuit by circuit.
  2028 (one tier): one parameter point, slope to 20% -> wall_campaign_s ~54 min; wall_first_result_s = None.
  2033 first result: 3-sigma detection of Delta n at T_c and 0.9 T_c (9/(S_int^2 m_max) = 2.73e5 queries per
  point at S_int = 1e-3) plus the slope at both temperatures -> wall_first_result_s ~0.58 yr.
  2033 campaign: 10 Delta-n points at +-1e-3 per bin over 100 bins (3.03e6 queries each) plus the 100-point
  slope scan -> wall_campaign_s ~32 yr, 6.3x the horizon; the scan plus the first result fit (3.8 yr).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, fields, replace

from estimates.common import (RFI_SRC, RFI_LQ_2028, RFI_T_2028, RFI_LQ_2033, RFI_T_2033, RFI_EPS_L,
                              SECONDS_PER_YEAR as _SECONDS_PER_YEAR)   # shared anchors (2026-10-09)
from estimates.common import (Assumed, Cited, Stated, Uncited, Tagged, Primitive,
                              CircuitStatus, Result, Published, ERAS,
                              EPS_SYN, EPS_SYN_SRC, SYNTHESIS_SRC, eps_rot_for, t_per_rotation, toffoli_t,
                              depth_exports, REACTION_TIME_S)

TEX = "app05"          # applications/app05_baryogenesis.tex
SEC_PER_YR = 365.25 * 86400.0
SEC_PER_DAY = 86400.0
DERIVED = "derived here (R12 item 14)"
SRC_SCALAR = "arxiv_2407_13819 Eq. rzTrotter (:1756-1762); arxiv_2210_07985 Eqs. trdispl-trphi4 (:374-435)"
SRC_KIN = "arxiv_2210_07985 Eqs. fftqft, trpi2 (:350-356, :397-407)"

# inputs that must be strictly positive for the model to mean anything
_POSITIVE = ("v_2028", "n_phi", "n_phi_small", "spinor_components", "lz_sites", "lperp_sites",
             "a_mphi", "v_family", "v_interim", "links_per_site", "lattice_dims", "eps_syn",
             "ct_reduction", "ct_reduction_banked", "n_trot_2028",
             "t_max_2033", "dt", "fault_budget_per_shot",
             "shots_2028", "shots_nucleation_2033", "eps_c", "n_bin", "envelope_2033", "cap_2028",
             "accuracy", "t_gate_s", "n_scan_points", "n_dn_points",
             "slope_times", "slope_early_fraction_2028", "slope_early_fraction_2033", "decrement_ref_time_2028",
             "detection_sigma", "n_points_first",
             "s_int", "t_depth_2028", "t_depth_2028_ref_steps", "t_depth_2033", "t_depth_2033_first",
             "sigma_over_mphi2", "bath_mixing_sweeps", "bath_mixing_sweeps_central", "n_trot_family",
             "horizon_yr", "t_2028_box", "am_f", "beta_mphi", "lr_patch_radius", "idle_cycle_s")


@dataclass(frozen=True)
class Assumptions:
    # ---- lattice and registers ------------------------------------------
    v_2028: Tagged = Stated(36, f"{TEX}:49,122", "6^2 (2+1)D lattice, 36 sites")
    n_phi: Tagged = Stated(16, f"{TEX}:122,157", "N_phi=16 (4 qubits/site); resolves both minima and the barrier")
    n_phi_small: Tagged = Stated(4, f"{TEX}:134,176", "N_phi=4 (2 qubits/site) in the family row and the interim variant")
    n_f_2033: Tagged = Stated(2, f"{TEX}:156,163", "2 Wilson-fermion flavors; 2N_f=4 fermion qubits/site")
    spinor_components: Tagged = Stated(2, f"{TEX}:95", "2^floor((D+1)/2)=2-component Dirac spinor in 2+1D")
    n_anc_2028: Tagged = Stated(30, f"{TEX}:93,126-127", "~30 ancilla: Trotter control, basin projector")
    n_anc_2033: Tagged = Stated(40, f"{TEX}:95,164", "40 ancilla (30 bath + 10 workspace)")
    n_bath_anc_2033: Tagged = Stated(30, f"{TEX}:159,164", "Lindblad bath, 30 shared ancillas (reset each step)")
    n_anc_generic: Tagged = Stated((10, 40), f"{TEX}:86", "'N_anc ~ 10-40' in the scaling section")
    lz_sites: Tagged = Stated(15, f"{TEX}:114,156", "15-site scattering axis")
    lperp_sites: Tagged = Stated(8, f"{TEX}:114,156", "8-site transverse direction")
    a_mphi: Tagged = Stated(0.5, f"{TEX}:157", "a m_phi = 0.5")
    v_family: Tagged = Stated(64, f"{TEX}:134", "8^2 family rows")
    v_interim: Tagged = Stated(16, f"{TEX}:176", "4^2 interim scalar+Wilson-fermion variant")
    workspace_interim_prior: Tagged = Stated(30, f"{TEX}:177", "'~30-qubit amplitude-estimation workspace' of a phase-register estimator, "
                                             "of which ~13 qubits are the phase register; MLAE needs no phase register")
    workspace_interim: Tagged = Stated(17, f"{TEX}:177", "'~17-qubit amplitude-estimation workspace ... a workspace budget, "
                                       "not an itemized count' = 30 - 13 (MLAE needs no phase register)")
    lq_interim_quoted: Tagged = Stated(115, f"{TEX}:177", "'The total is ~115 LQ'; 96+17=113")
    lq_two_register_quoted: Tagged = Stated(190, f"{TEX}:177", "'two-register protocol ... ~190': 96+96=192")

    # ---- c_T derivation inputs ---------------------------------
    lattice_dims: Tagged = Stated(2, f"{TEX}:49,156", "2+1D: d = 2 hopping directions")
    links_per_site: Tagged = Assumed(2, "d = 2 links per site, as for periodic boundaries (arxiv_2407_13819 E_D); periodic on the "
                                     "square 2028 lattices; the 2033 15x8 cylinder is open along its 15-site axis, "
                                     "8 of 240 links fewer (~1% of c_T), so the 2033 count is slightly conservative", "derived here")
    wilson_r: Tagged = Assumed(1, "Wilson parameter r = 1: the hopping matrix (r beta - i alpha_j)/2 has rank 1", "derived here")
    eps_syn: Tagged = Cited(EPS_SYN, EPS_SYN_SRC, "report-wide: total synthesis error per circuit 1e-2; eps_rot = sqrt(eps_syn/N_rot)")
    rot_synthesis: Tagged = Cited("rus", SYNTHESIS_SRC, "report-wide: 1.15 log2(1/eps_rot) + 9.2 T per rotation")
    toffoli_convention: Tagged = Cited("textbook", "main-overview:47 (report-wide)", "7 T per Toffoli everywhere")
    # printed c_T values (echoes of the derived numbers; checked against the model in the tests)
    c_t_2028_quoted: Tagged = Stated(1.26e3, "r12 record", "c_T 1263.1 T/site/step at 20 steps (the rotations of a longer shot "
                                     "are synthesized tighter; 1239.8 at 12 steps); derived here, no longer printed since the Oct 6 cut")
    c_t_2033_quoted: Tagged = Stated(2.5e3, f"{TEX}:91,165", "'c_T ~ 2.5e3 (1.6e3 scalar, 8.9e2 fermion)' at the deepest circuit's budget")
    c_t_range_quoted: Tagged = Stated((2.5e2, 2.5e3), f"{TEX}:85", "'c_T runs from 2.5e2 to 2.5e3 across the instances below'")
    ct_reduction: Tagged = Assumed((3, 4), "model-internal sanity span for the banked ~3.3x (not printed in the chapter); "
                                   "assumed development target, not derived or cited", "")
    ct_reduction_banked: Tagged = Assumed(3.3, "assumed development target, not derived or cited: "
                                          "'divided by an assumed ~3.3x reduction that is a development target, not derived "
                                          "or cited'; carried exactly, printed results rounded once",
                                          f"{TEX}:91,126,162,163,179")
    t_2028_box: Tagged = Stated(2.8e5, f"{TEX}:94", "'2.8e5 (9.1e5 before the assumed 3.3x step reduction)': 9.095e5/3.3 "
                                "= 2.756e5 at 20 steps (ruling 2026-10-11; was 1.6e5 at 12 steps)")
    t_with_fermions_6x6_raw_quoted: Tagged = Stated(8.5e5, f"{TEX}:91", "'8.5e5 T raw': 1973.0 x 36 x 12 = 8.523e5")
    t_with_fermions_6x6_quoted: Tagged = Stated(2.6e5, f"{TEX}:91", "'~2.6e5 T with the same reduction applied to every sector'")
    c_t_with_fermions_6x6_quoted: Tagged = Stated(2.0e3, f"{TEX}:91", "'c_T ~ 2.0e3'")
    system_qubits_with_fermions_6x6_quoted: Tagged = Stated(288, f"{TEX}:91", "'288 system qubits' = 36 x (4 + 4)")

    # ---- Trotter --------------------------------------------------------------
    n_trot_2028: Tagged = Stated(20, f"{TEX}:59,97", "'20 steps at the same delta t', to t = 2/m_phi: both slope times "
                                 "(10 and 20 steps) inside the chapter's 1-2/m_phi exponential window (ruling "
                                 "2026-10-11, was 12 steps with the early time 0.4/m_phi in the Zeno regime)")
    n_trot_family: Tagged = Stated(12, "r16 record", "the 8^2 family rows and the interim 4^2 variant (no longer printed) "
                                   "were set at 12 steps to t = 1.2/m_phi; kept there as records")
    t_max_2033: Tagged = Stated(10.0, f"{TEX}:82,100", "t_max = 10/m_phi (units of 1/m_phi)")
    dt: Tagged = Stated(0.1, f"{TEX}:82,91", "'delta t = 0.1/m_phi'")
    t_window_2028: Tagged = Stated((1.0, 2.0), f"{TEX}:91,101,142", "'t ~ 1-2/m_phi', Zeno-to-exponential crossover")

    # ---- shots and readout ----------------------------------------------------
    shots_2028: Tagged = Stated(5e3, "r23 record", "single-time 20% count at t = 1.2/m_phi (r23 box '5e3'); it fixes the "
                                "ln P_FV decrement there, delta = 1/(N eps^2) = 5e-3, from which the slope shots follow")
    decrement_ref_time_2028: Tagged = Assumed(1.2, "time (1/m_phi) at which shots_2028 set the decrement; in the exponential "
                                              "regime the decrement grows linearly with t", "derived here 2026-10-11")
    n_rare_range: Tagged = Stated((1e3, 1e4), f"{TEX}:83", "N_rare ~ 1e3-1e4 in Eq. Nshot_baryo")
    shots_nucleation_2033: Tagged = Stated(1e4, f"{TEX}:169", "nucleation arm '1e4 shots' per grid point")
    eps_c: Tagged = Stated(1e-3, f"{TEX}:88,160", "reflection asymmetry ~1e-3 at Delta theta_C = 0.1")
    n_bin: Tagged = Stated(8, f"{TEX}:63,115", "'the 8 bins the 15-site axis resolves': floor(15/2) + 1 distinct |k_z| "
                           "(ruling 2026-10-11, batch 4 item 3b; was 'N_bin ~ 1e2', which the lattice cannot resolve)")
    mlae_queries_per_bin_quoted: Tagged = Stated(7.0e4, f"{TEX}:83", "'~7.0e4 base-circuit queries per kinematic bin' (C_est = 1, "
                                                 "x 1/p0^2 flag contrast)")
    c_est: Tagged = Stated(1.0, f"{TEX}:88", "author-adopted branch prior |A_R| < sin(pi/66) = 0.048: "
                           "every circuit at m_max, Cramer-Rao constant one (Suzuki Fisher_final)")
    dn_branch_prior: Tagged = Stated(0.048, f"{TEX}:88", "'We adopt this as a prior on |A_R|': |A_R| < sin(pi/66) (author choice)")
    mlae_schedule: Tagged = Cited("single", "arxiv_1904_10246 Manuscript_v2.tex:238",
                                  "all circuits at m_max; LIS m = 1, 3, ..., m_max (:268, :271) kept as sensitivity")
    dn_estimate: Tagged = Stated(1e-3, f"{TEX}:88", "'|R_p|^2 - |R_pbar|^2 ~ ... ~ 1e-3 at Delta theta_C = 0.1' (an estimate, not a bound)")
    dn_state_independent_bound: Tagged = Assumed(1.0, "no derived bound on |Delta| follows from the chapter's inputs (the 'below Delta theta_C' "
                                                 "scaling, taken as a bound, gives 0.1 > 0.048); recorded only: the "
                                                 "single depth rests on the author prior (c_est), not on a derived bound", "derived here")
    m_max_quoted: Tagged = Stated(33, f"{TEX}:88,166", "'m_max = floor(1e9/3.0e7) = 33' (fixed point)")
    envelope_2033: Tagged = Cited(RFI_T_2033, "DOE RFI 2026", f"1e9 hard-op envelope, {TEX}:88,91,166")
    cap_2028: Tagged = Cited(RFI_T_2028, "DOE RFI 2026", f"1e5 T cap, {TEX}:128,181")
    accuracy: Tagged = Stated(0.2, f"{TEX}:49,160,201", "20% relative accuracy on Gamma/V")
    fault_budget_per_shot: Tagged = Assumed(0.1, "0.1 expected faults per shot, report-wide",
                                            "main:272 (report-wide)")
    eps_l_2028_quoted: Tagged = Stated(4e-7, f"{TEX}:96", "'Required epsilon_l <~ 4e-7' (0.1 expected faults per 2.7e5-T shot; "
                                       "was 6e-7 at 12 steps)")
    eps_l_2033_deepest_quoted: Tagged = Stated(1e-10, f"{TEX}:167", "'<~ 1e-10 for the deepest MLAE circuit (33 x 3.0e7 T)'")
    eps_l_2033_base_quoted: Tagged = Stated(3e-9, f"{TEX}:167", "'<~ 3e-9 for the base circuit'")

    # ---- wall time and campaign -------------------------------------------------
    t_gate_s: Tagged = Stated(1e-6, f"{TEX}:89,91", "'1 us per T-gate' (seconds per T)")
    shot_overhead_s: Tagged = Assumed(1e-4, "per-shot overhead t0 ~ 0.1 ms: register initialization (~1 logical cycle, 10 us), final transversal readout and decode (63 us decoder latency, Google_QEC_below_threshold); the wall is gate-limited, with no separate shot rate", f"{TEX}:130")
    n_scan_points: Tagged = Stated(100, f"{TEX}:89,151,176", "'~100-point (T, Delta theta_C) scan'")
    n_dn_points: Tagged = Stated(10, f"{TEX}:89", "Delta n arm 'at ~10 representative points'")
    horizon_yr: Tagged = Stated(5, f"{TEX}:89,151,175", "'5-year campaign horizon', one machine, serial (no machine count is formed)")

    # ---- two-tier schedule (shot audit and factory analysis) ----------------
    slope_times: Tagged = Assumed(2, "Gamma/V is the slope of ln P_FV past the Zeno time, so two times (t/3 and t), "
                                  "not one", "shot_audit.json ch8 (1)")
    slope_early_fraction_2028: Tagged = Assumed(0.5, "2028 early time t/2: 10 of 20 steps (1.0 and 2.0/m_phi), both in the "
                                                "chapter's exponential window; shots per time follow from slope_shots "
                                                "(was t/3, 4 of 12 steps, with 0.4/m_phi in the Zeno regime)",
                                                "ruling 2026-10-11")
    slope_early_fraction_2033: Tagged = Assumed(1 / 3, "2033 early time t/3: 33 of 100 steps; a long lever arm is ~5x "
                                                "cheaper than t/1.5", "shot_audit.json ch8 lever (5)")
    detection_sigma: Tagged = Assumed(3, "the 2033 first result is a 3-sigma detection of Delta n", "H. Lamm 2026-10-02")
    n_points_first: Tagged = Stated(2, f"{TEX}:143", "first result at T = T_c and 0.9 T_c (the box's two temperatures)")
    s_int: Tagged = Assumed(1e-3, "f_eq-weighted integrated asymmetry; NOT computed; taken at the per-k estimate "
                            "(app05:88); the first-result Delta-n wall scales as (1e-3/S_int)^2", "NEEDS_AUTHOR")
    label_register_first: Tagged = Assumed(7, "packet-label register of the f_eq-weighted packet superposition, "
                                           "taken from the 10-qubit workspace", "shot_audit.json ch8 lever (1)")
    t_depth_2028: Tagged = Assumed((5.4e3, 7.0e3), "T-depth per 12-step shot, banked: ceil(1980/29) x 21.38 + ... per step "
                                   "(A = 29-30 RUS ancillas) to the serial structural schedule, x 12 / 3.3; the model "
                                   "scales it by n_trot_2028 / t_depth_2028_ref_steps", "factory.json ch08 run 1")
    t_depth_2028_ref_steps: Tagged = Assumed(12, "step count of the factory.json run behind t_depth_2028; depth is "
                                             "proportional to the step count", "factory.json ch08 run 1")
    t_depth_2033: Tagged = Assumed((7.2e5, 2.87e6), "T-depth per base circuit: A = 40 (bath ancillas borrowed during "
                                   "the closed scattering segments) to A = 10 (workspace only)", "factory.json ch08 runs 2-5")
    t_depth_2033_first: Tagged = Assumed((8.7e5, 9.6e6), "T-depth per first-result Delta-n base circuit: the 7-qubit "
                                         "label leaves A = 33 (bath borrowed) to A = 3 (workspace only). Scaled from "
                                         "t_depth_2033 by the rotation-layer count per step, ceil(10440/A): "
                                         "7.2e5 x 317/261 = 8.74e5 and 2.87e6 x 3480/1044 = 9.57e6",
                                         "factory.json ch08 runs 2-5, rescaled")

    # ---- state preparation, flag contrast, bath and fault model --------------------
    am_f: Tagged = Assumed(0.5, "fermion Wilson mass a m_f in the symmetric phase: not fixed by the chapter, taken equal to "
                           "a m_phi; flag contrast p0 = 0.57, 0.66, 0.77 at a m_f = 0.25, 0.5, 1 (flag_contrast)",
                           "NEEDS_AUTHOR (flag contrast)")
    slater_givens: Tagged = Cited("(N - eta) eta", "arxiv_1711_05395 Eq. n_gates (fhm_arXiv4.tex:346-350)",
                                  "Givens rotations for one Slater determinant of eta particles in N modes")
    fft8_rotation_free: Tagged = Cited(True, "arxiv_1902_10673 main.tex:1102",
                                       "FFFT of side 8 needs no arbitrary rotations; F_2 = two pi/8 Pauli rotations "
                                       "= 2 T, twiddles powers of T (<= 1 T each), derived here")
    beta_mphi: Tagged = Stated(1.0, f"{TEX}:99", "'plasma correlation time ~1/T ~ 1/m_phi' at T ~ T_c: beta = 1/m_phi")
    lr_patch_radius: Tagged = Assumed(2, "Lieb-Robinson radius v_LR beta in sites: beta = 1/m_phi = 2a at a m_phi = 0.5, "
                                      "v_LR ~ 1 site per unit a; patch (2R+1)^2 = 25 sites", "derived here")
    gibbs_cost_model: Tagged = Cited("t beta", "arxiv_2311_09207 Thm. L_cost and lattice estimate (main.tex:398-436)",
                                     "Hamiltonian simulation time ~beta per unit Lindblad time on a patch of radius "
                                     "v_LR beta; jumps normalized ||sum_a A^a+ A^a|| <= 1 (main.tex:346), adjoints included")
    idle_cycle_s: Tagged = Assumed(1e-5, "one logical cycle ~10 us (as in shot_overhead_s); idle qubit-cycles = "
                                      "LQ x circuit time / cycle", "derived here")
    bath_mixing_sweeps: Tagged = Assumed((3, 20), "sweeps of the filtered dissipator to mix the Gibbs state, a dynamic-"
                                         "exponent heuristic (not measured): relaxation time tau = xi^z = 4-4.5 sweeps at "
                                         "xi = 1/m_phi = 2a and z = 2-2.17 (model A); warm start from a mean-field thermal "
                                         "state 2-3 tau = 3-14 sweeps, cold start to 1% local accuracy 7-21; printed as "
                                         "words only, the nucleation walls stay unpriced lower bounds (ruling 2026-10-11)",
                                         "estimate 2026-10-10 (scratch ch8_heur)")
    bath_mixing_sweeps_central: Tagged = Assumed(8, "central value of bath_mixing_sweeps (record; nothing printed uses it)",
                                                 "estimate 2026-10-10")

    # ---- finite-volume saddle behind the nucleation decrement (ruling 2026-10-11) --------------
    sigma_over_mphi2: Tagged = Assumed(1.32, "wall tension sigma / m_phi^2 of the 2033 scalar potential, chosen so that the "
                                       "wall-pair (slab) saddle on the 15x8 cylinder reproduces the decrement 2.5e-3 at "
                                       "t = 10/m_phi behind the printed 1e4 single-time shots; thin-wall phi^4: "
                                       "sigma = m^3/(3 lambda), so lambda/m_phi = 0.25. The chapter never fixed sigma or "
                                       "lambda; this is the stated assumption", "derived here 2026-10-10 (scratch ch8_heur)")
    latent_heat_ell: Tagged = Assumed(1.0, "free-energy difference Delta F = ell m_phi^3 (1 - T/T_c) between the wells "
                                      "(thin wall, T_c = m_phi from beta_mphi); lowers the slab action by ell (1 - T/T_c) "
                                      "L_perp over a strip of width 1/m_phi", "assumed 2026-10-10")

    # ---- interim variant (printed echoes) ------------------------------------------
    t_interim_quoted: Tagged = Stated(1.4e5, f"{TEX}:178,181", "'~1.4e5 T raw (c_T ~ 7.5e2: 2.5e2 + 5.0e2, x 16 x 12)'")
    gap_interim_quoted: Tagged = Stated(1.4, f"{TEX}:178,181", "'~1.4x above the 2028 hard-op cap'")
    t_interim_banked_quoted: Tagged = Stated(4.4e4, f"{TEX}:181", "'brings it to ~4.4e4 if applied to every sector'")
    t_interim_scalar_share_quoted: Tagged = Stated(1.1e5, f"{TEX}:181", "'~1.1e5 if applied to the scalar share alone, "
                                                   "below the benchmark's 1.6e5'")

    # ---- utility box: no dollar figure (referee G5, 2026-10-04; paper-wide ruling 2026-10-09) ----------------

    def __post_init__(self):
        for f in fields(self):
            v = getattr(self, f.name)
            if not isinstance(v, Tagged):
                raise TypeError(f"{f.name} must be Tagged")
            if isinstance(v.value, str):
                continue
            if v.is_range and not (v.lo <= v.hi):
                raise ValueError(f"{f.name}: range {v.value} is inverted")
            if v.lo < 0:
                raise ValueError(f"{f.name}: negative")
            if f.name in _POSITIVE and v.lo <= 0:
                raise ValueError(f"{f.name}: must be positive")
        for name in ("n_phi", "n_phi_small"):
            n = int(getattr(self, name).lo)
            if n < 2 or n & (n - 1):
                raise ValueError(f"{name}={n} must be a power of two >= 2; log2 N_phi is the qubits per site")
        if self.spinor_components.lo not in (2, 4):
            raise ValueError("spinor must be 2- or 4-component")
        if self.wilson_r.lo != 1:
            raise ValueError("the hopping count uses the rank-1 Wilson projector, r = 1")
        if self.n_bath_anc_2033.lo > self.n_anc_2033.lo:
            raise ValueError("bath ancillas exceed the total 2033 ancilla budget (negative workspace)")
        if not (self.accuracy.lo < 1 and self.eps_c.lo < 1):
            raise ValueError("accuracy and eps_c are fractions < 1")
        if not (self.ct_reduction.lo <= self.ct_reduction_banked.lo <= self.ct_reduction.hi):
            raise ValueError(f"ct_reduction_banked={self.ct_reduction_banked.lo} is outside the "
                             f"{self.ct_reduction.value} span the box quotes")
        if not (0 < self.fault_budget_per_shot.lo < 1):
            raise ValueError("fault_budget_per_shot is an expected fault count per shot, in (0, 1)")
        for name in ("slope_early_fraction_2028", "slope_early_fraction_2033"):
            if not (0 < getattr(self, name).lo < 1):
                raise ValueError(f"{name} is the early time over the late time, in (0, 1)")
        if int(self.slope_times.lo) != 2:
            raise ValueError("the slope is priced from exactly two times")
        if int(self.n_bin.lo) > int(self.lz_sites.lo) // 2 + 1:
            raise ValueError(f"n_bin={int(self.n_bin.lo)} exceeds the {int(self.lz_sites.lo) // 2 + 1} distinct |k_z| "
                             "the scattering axis resolves")
        if not (self.bath_mixing_sweeps.lo <= self.bath_mixing_sweeps_central.lo <= self.bath_mixing_sweeps.hi):
            raise ValueError("bath_mixing_sweeps_central is outside bath_mixing_sweeps")
        if self.label_register_first.lo > self.n_anc_2033.lo - self.n_bath_anc_2033.lo:
            raise ValueError("the packet-label register does not fit in the 2033 workspace")


# --------------------------------------------------------------------------- #
# Gate-level counts per site per Trotter step (derived here; see module docstring)
# --------------------------------------------------------------------------- #

def qubits_per_site(n_phi: int, n_f: int, components: int = 2) -> int:
    """log2 N_phi + N_f * components  (the chapter's 2 N_f is the 2-component case)."""
    return int(math.log2(n_phi)) + n_f * components


def qft_cost(n_q: int, tof_t: float) -> tuple[int, float]:
    """(synthesized rotations, exact T) of one exact n_q-qubit QFT, derived here.

    Controlled phase pi/2^k between qubits at distance k (n_q - k of them): k=1 controlled-S, 3 T;
    k=2 controlled-T via a Toffoli into an ancilla plus T, tof_t + 1 T; k>=3 Toffoli plus one synthesized
    R_Z(pi/2^k). The AND ancilla is uncomputed by measurement (0 T)."""
    rot, t = 0, 0.0
    for k in range(1, n_q):
        n = n_q - k
        if k == 1:
            t += 3 * n
        elif k == 2:
            t += (tof_t + 1) * n
        else:
            t += tof_t * n
            rot += n
    return rot, t


def scalar_counts(n_q: int, links: int, tof_t: float) -> dict:
    """Per site per step: synthesized rotations by term and exact T of the N_phi = 2^n_q scalar."""
    q_rot, q_t = qft_cost(n_q, tof_t)
    rot = {
        "onsite_polynomial": sum(math.comb(n_q, w) for w in range(1, min(4, n_q) + 1)),
        "gradient_links": links * n_q * n_q,
        "kinetic_pi2_zz": math.comb(n_q, 2),
        "qft_controlled_phase": 2 * q_rot,
    }
    return {"rot": rot, "t_exact": {"qft_exact": 2 * q_t}}


def fermion_counts(n_q: int, n_f: int, components: int, links: int, dims: int) -> dict:
    """Per site per step: rotations and exact T of n_f Wilson flavors (r = 1) plus the Yukawa to the
    n_q-qubit scalar, under Jordan-Wigner. Counts scale with components/2 (rank of the Wilson projector,
    number of mass and Yukawa bilinears). The Yukawa needs L/R pairs: n_f // 2 of them."""
    h = components // 2
    pairs = n_f // 2
    rot = {
        "wilson_hopping": links * n_f * 2 * h,          # (XX+YY)/2 x Z-string per link per flavor, per rank
        "mass_wilson_onsite": n_f * components,          # (m+2r)(n_1 - n_2) per flavor
        "yukawa_phi_bilinear": pairs * components * 2 * n_q,   # phi (n_q Z) x (XX+YY)/2 per component
        "yukawa_phase_frame": pairs * components * 2,    # e^{i theta_C(z) n_L} in and out per component
    }
    t_exact = {"hopping_basis_change": dims * n_f * 4 * h}   # pi/4 Givens (2 T) in and out per direction
    return {"rot": rot, "t_exact": t_exact}


def price_instance(a: Assumptions, v: int, n_steps: int, n_phi: int, n_f: int, components: int,
                   reduction: float = 1.0, budget_circuits: int = 1, extra_rot: int = 0) -> dict:
    """Price one instance: counts per site-step, N_rot per shot, eps_rot (R-TOL), c_T, raw T, and
    Primitives (each t_each divided by `reduction`, the banked development-target factor, if any).

    `budget_circuits` = number of base-circuit applications in the circuit whose eps_syn budget the
    rotations must meet (1 for a single circuit; m_max for the deepest MLAE circuit, R13):
    eps_rot = sqrt(eps_syn / (budget_circuits x N_rot(base))).
    `extra_rot` = synthesized rotations per application outside the Trotter steps (the priced B2 preparation),
    counted into the budget: N_rot = budget_circuits x (N_rot(base) + extra_rot)."""
    n_q = int(math.log2(n_phi))
    links = int(a.links_per_site.lo)
    dims = int(a.lattice_dims.lo)
    tof_t = toffoli_t(1, a.toffoli_convention.value)
    sc = scalar_counts(n_q, links, tof_t)
    fc = fermion_counts(n_q, n_f, components, links, dims) if n_f else {"rot": {}, "t_exact": {}}
    site_steps = v * n_steps
    rot_scalar = sum(sc["rot"].values())
    rot_ferm = sum(fc["rot"].values())
    t_scalar = sum(sc["t_exact"].values())
    t_ferm = sum(fc["t_exact"].values())
    n_rot_site_step = rot_scalar + rot_ferm
    n_rot_shot = n_rot_site_step * site_steps
    if budget_circuits < 1:
        raise ValueError("budget_circuits must be >= 1")
    n_rot_budget = (n_rot_shot + extra_rot) * budget_circuits
    eps_rot = eps_rot_for(n_rot_budget, float(a.eps_syn.lo))
    t_rot = t_per_rotation(eps_rot, a.rot_synthesis.value)
    ct_scalar = rot_scalar * t_rot + t_scalar
    ct_ferm = rot_ferm * t_rot + t_ferm
    c_t = ct_scalar + ct_ferm
    prims = []
    for sector, counts in (("scalar", sc), ("fermion", fc)):
        src = SRC_SCALAR if sector == "scalar" else DERIVED
        for name, n in counts["rot"].items():
            s = SRC_KIN if name == "kinetic_pi2_zz" else (DERIVED if name == "qft_controlled_phase" else src)
            prims.append(Primitive(f"{sector}_{name}", n * site_steps, t_rot / reduction, CircuitStatus.COMPILED, s,
                                   f"{n} R_Z/site/step at eps_rot={eps_rot:.3g} ({t_rot:.3f} T each)"
                                   + (f", synthesized to the {budget_circuits}-application circuit's budget"
                                      if budget_circuits != 1 else "")
                                   + (f", / banked {reduction:g}x (uncited development target)" if reduction != 1 else "")))
        for name, t in counts["t_exact"].items():
            prims.append(Primitive(f"{sector}_{name}", site_steps, t / reduction, CircuitStatus.COMPILED, DERIVED,
                                   f"{t:g} exact T/site/step (Toffoli = {tof_t:g} T)"
                                   + (f", / banked {reduction:g}x" if reduction != 1 else "")))
    return {
        "n_q": n_q, "site_steps": site_steps,
        "rot_scalar": rot_scalar, "rot_fermion": rot_ferm, "t_exact_scalar": t_scalar, "t_exact_fermion": t_ferm,
        "rot_by_term": {**{f"scalar_{k}": x for k, x in sc["rot"].items()}, **{f"fermion_{k}": x for k, x in fc["rot"].items()}},
        "n_rot_site_step": n_rot_site_step, "n_rot_shot": n_rot_shot, "n_rot_budget": n_rot_budget,
        "budget_circuits": budget_circuits, "eps_rot": eps_rot, "t_rot": t_rot,
        "c_t_scalar": ct_scalar, "c_t_fermion": ct_ferm, "c_t": c_t, "t_raw": c_t * site_steps,
        "primitives": tuple(prims),
    }


def slope_shots(eps: float, delta_late: float, r: float) -> float:
    """Shots at EACH of two times t1 = r t2 and t2 for a relative accuracy eps on the slope of ln P_FV.

    Binomial statistics at small decrement: var(ln P_i) = delta_i / N with delta_1 = r delta_2 (the decrement grows
    linearly with t in the exponential regime). The slope estimate (ln P_1 - ln P_2) / (t_2 - t_1) then has
    relative variance (1 + r) / (N delta_2 (1 - r)^2), so N = (1 + r) / (eps^2 delta_2 (1 - r)^2).
    At r = 1/3 this is 3x the single-time count 1 / (eps^2 delta_2) (the r25 'slope_shot_factor = 3');
    at r = 1/2 it is 6x. Derived here 2026-10-11."""
    if not (0 < r < 1) or not (0 < delta_late < 1) or eps <= 0:
        raise ValueError("need 0 < r < 1, 0 < delta_late < 1, eps > 0")
    return (1.0 + r) / (eps ** 2 * delta_late * (1.0 - r) ** 2)


def slab_saddle(a: "Assumptions", t_over_tc: float, t_mphi: float, sigma: float | None = None) -> dict:
    """Finite-volume nucleation on the 15 x 8 cylinder from a wall-pair (slab) saddle (heuristic, derived here
    2026-10-10). A thin-wall bubble of tension s = sigma/m^2 against Delta F = ell m^3 x, x = 1 - T/T_c, has
    S_2/T >= 4 pi s^2 / (ell T) (at x = 1/2) and a radius s/(ell x) that must fit in L_perp/2 = 2/m_phi; on this
    lattice it never does, and the saddle is a pair of walls wrapping the L_perp = 8a direction with barrier
    2 s L_perp, lowered by ell x L_perp over a strip of width 1/m. Prefactor: translational zero mode along L_z and
    the Gaussian factor sqrt(S/2 pi) at scale m_phi. Everything in units of m_phi (T_c = m_phi, beta_mphi = 1).
    Returns the action over T, the decrement Gamma t over the window t_mphi (in 1/m_phi), and the bubble bound."""
    s = float(a.sigma_over_mphi2.lo) if sigma is None else float(sigma)
    ell = float(a.latent_heat_ell.lo)
    am = float(a.a_mphi.lo)
    lperp = int(a.lperp_sites.lo) * am                        # 4.0 / m_phi
    lz = int(a.lz_sites.lo) * am                              # 7.5 / m_phi
    x = 1.0 - t_over_tc
    action_over_t = (2.0 * s * lperp - ell * x * lperp) / t_over_tc
    rate = lz * math.sqrt(action_over_t / (2.0 * math.pi)) * math.exp(-action_over_t)   # per unit 1/m_phi time
    bubble_min_action_over_t = 4.0 * math.pi * s ** 2 / (ell * t_over_tc)
    return {"action_over_t": action_over_t, "decrement": rate * t_mphi,
            "bubble_min_action_over_t": bubble_min_action_over_t,
            "slab_dominates": action_over_t < bubble_min_action_over_t,
            "lambda_over_mphi_thin_wall": 1.0 / (3.0 * s)}


def bounce_action_implied(gamma_over_area: float) -> float:
    """Vacuum (2+1D, three Euclidean dimensions) S_E with (S/2 pi)^{3/2} e^{-S} = Gamma/A in units of m_phi;
    bisection on [1, 60] (record for the 2028 benchmark; derived here 2026-10-10)."""
    f = lambda x: (x / (2 * math.pi)) ** 1.5 * math.exp(-x) - gamma_over_area
    lo, hi = 1.0, 60.0
    if f(lo) < 0 or f(hi) > 0:
        raise ValueError("Gamma/A outside the bracket")
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if f(mid) > 0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def _banked(raw: float, red: Tagged) -> tuple[float, float]:
    """The 3-4x span the 2028 box quotes, applied to a raw count: (raw/4, raw/3)."""
    return (raw / red.hi, raw / red.lo)


def eps_l_required(t_ops: float, fault_budget: float) -> float:
    """R3: required logical error rate for `fault_budget` expected faults over `t_ops` hard ops."""
    return fault_budget / t_ops


def _within(x: float, rng: Tagged) -> bool:
    return rng.lo <= x <= rng.hi


# --------------------------------------------------------------------------- #
# 2028
# --------------------------------------------------------------------------- #

def _model_2028(a: Assumptions) -> Result:
    v = int(a.v_2028.lo)
    n_phi = int(a.n_phi.lo)
    qps = qubits_per_site(n_phi, 0)
    system = v * qps
    lq = system + int(a.n_anc_2028.lo)
    n_trot = int(a.n_trot_2028.lo)
    red = float(a.ct_reduction_banked.lo)                   # 3.3, the one factor the box banks
    p = price_instance(a, v, n_trot, n_phi, 0, 2, reduction=red)
    c_t = p["c_t"]                                          # 1263.1 at 20 steps (1239.8 at 12)
    site_steps = p["site_steps"]
    raw = p["t_raw"]                                        # 9.095e5 (5.356e5 at 12 steps)
    t_point = raw / red                                     # 1.623e5 -> '~1.6e5'
    banked = _banked(raw, a.ct_reduction)                   # (1.34e5, 1.79e5)
    cap = float(a.cap_2028.lo)
    dt = float(a.dt.lo)
    t_window = n_trot * dt
    steps_from_window = (a.t_window_2028.lo / dt, a.t_window_2028.hi / dt)
    box_t = float(a.t_2028_box.lo)

    shots = float(a.shots_2028.lo)
    gt = float(a.t_gate_s.lo)
    t0 = float(a.shot_overhead_s.lo)
    gates_only_s = t_point * gt * shots                     # 811.5 s, T gates alone
    gate_time_raw_s = raw * gt * shots                      # 2678 s, T gates alone
    wall_s = shots * (t_point * gt + t0)                    # 812.0 s = 13.5 min -> '~14 min'
    wall_raw_s = shots * (raw * gt + t0)                    # 2678.5 s = 44.6 min -> '~45 min at the unreduced c_T'

    # R25 (R2 + shot audit), re-timed 2026-10-11: the slope of ln P_FV from two times t/2 and t = 2/m_phi (10 and
    # 20 steps), both inside the chapter's 1-2/m_phi exponential window (the r25 schedule, 0.4 and 1.2/m_phi, put the
    # early time in the Zeno regime). Shots per time from slope_shots at the decrement the r23 single-time count
    # fixed: delta(1.2/m_phi) = 1/(5e3 x 0.2^2) = 5e-3, growing linearly with t. Each circuit at its own depth and
    # its own synthesis budget.
    eps = float(a.accuracy.lo)
    n_early = int(round(n_trot * float(a.slope_early_fraction_2028.lo)))     # 10 steps
    r_times = n_early / n_trot                                              # 0.5
    delta_ref = 1.0 / (shots * eps ** 2)                                    # 5e-3 at 1.2/m_phi
    delta_late = delta_ref * t_window / float(a.decrement_ref_time_2028.lo)  # 8.33e-3 at 2.0/m_phi
    p_early = price_instance(a, v, n_early, n_phi, 0, 2, reduction=red)
    t_early = p_early["t_raw"] / red                                        # 1.34e5
    shots_per_time = slope_shots(eps, delta_late, r_times)                  # 1.8e4
    slope_shot_factor = shots_per_time * eps ** 2 * delta_late              # 6 (was 3 at r = 1/3)
    slope_wall_s = shots_per_time * ((t_point + t_early) * gt + 2 * t0)      # 7.29e3 s = 2.0 h
    slope_wall_raw_s = shots_per_time * ((raw + p_early["t_raw"]) * gt + 2 * t0)   # 2.4e4 s = 6.7 h
    shots_equiv = shots_per_time * (1 + t_early / t_point)                  # full-circuit equivalents for depth_exports
    depth_scale = n_trot / float(a.t_depth_2028_ref_steps.lo)               # 20 / 12: depth grows with the steps
    t_depth = (float(a.t_depth_2028.lo) * depth_scale, float(a.t_depth_2028.hi) * depth_scale)   # 9.0e3-1.17e4
    dx = depth_exports(t_point, t_depth, shots_equiv, a.t_gate_s, a.shot_overhead_s)
    # record: the vacuum bounce action the decrement implies on 6^2 (A = 36 a^2 = 9/m^2): Gamma/A = delta/(A t)
    area_2028 = v * float(a.a_mphi.lo) ** 2
    bounce_action = bounce_action_implied(delta_late / (area_2028 * t_window))   # ~8.0

    # family rows at the r16 record's 12 steps (n_trot_family; no longer printed)
    vf = int(a.v_family.lo)
    n_small = int(a.n_phi_small.lo)
    n_fam = int(a.n_trot_family.lo)
    fam4 = price_instance(a, vf, n_fam, n_small, 0, 2)
    fam16 = price_instance(a, vf, n_fam, n_phi, 0, 2)
    fam4_lq = vf * qubits_per_site(n_small, 0) + int(a.n_anc_2028.lo)
    fam16_lq = vf * qps + int(a.n_anc_2028.lo)

    # the fermion arm at N_phi=16 on 6^2 (app05:91): a fact, not the reason it waits
    n_f = int(a.n_f_2033.lo)
    comp = int(a.spinor_components.lo)
    wf = price_instance(a, v, n_fam, n_phi, n_f, comp)
    with_f_raw = wf["t_raw"]                                # 8.523e5 (12-step record)
    with_f_banked_all_point = with_f_raw / red              # 2.583e5 -> '~2.6e5'
    with_f_system = v * qubits_per_site(n_phi, n_f, comp)   # 288

    # interim scalar-fermion variant
    vi = int(a.v_interim.lo)
    interim_qps = qubits_per_site(n_small, n_f, comp)
    interim_system = vi * interim_qps                        # 96
    phase_register_bits = math.ceil(math.log2(1.0 / (float(a.accuracy.lo) * float(a.eps_c.lo))))  # 13, dropped in R11b
    workspace_from_prior = int(a.workspace_interim_prior.lo) - phase_register_bits   # 17
    interim_lq = interim_system + int(a.workspace_interim.lo)  # 113 -> '~115'
    interim_lq_with_trotter_anc = interim_lq + int(a.n_anc_2028.lo)
    two_register_candidates = (interim_system + interim_system, interim_system + 2 * int(a.workspace_interim.lo))
    ip = price_instance(a, vi, n_fam, n_small, n_f, comp)
    interim_t = ip["t_raw"]                                  # 1.449e5 -> '~1.4e5' (12-step record)
    interim_t_banked_all = interim_t / red                   # 4.39e4 -> '~4.4e4'
    interim_t_banked_scalar_only = (ip["c_t_scalar"] / red + ip["c_t_fermion"]) * ip["site_steps"]   # 1.114e5
    gap = interim_t / cap                                    # 1.449 -> '~1.4'

    breakdown = p["primitives"] + (
        Primitive("adiabatic_state_prep", 1, 0.0, CircuitStatus.UNSOURCED, f"{TEX}:65,107,123",
                  "adiabatic ramp to the false vacuum; not priced in the chapter"),
        Primitive("basin_projector_and_readout", 1, 0.0, CircuitStatus.UNSOURCED, f"{TEX}:127",
                  "metastable-basin projector; direct-sampling readout; not priced"),
    )
    eps_l = eps_l_required(t_point, float(a.fault_budget_per_shot.lo))   # 6.16e-7
    inter = {
        "qubits_per_site": qps, "system_qubits": system, "n_anc": int(a.n_anc_2028.lo), "lq": lq,
        "rot_per_site_step": p["n_rot_site_step"],            # 55
        "rot_by_term": p["rot_by_term"],
        "t_exact_per_site_step": p["t_exact_scalar"],         # 64
        "n_rot_shot": p["n_rot_shot"],                        # 23760
        "eps_rot": p["eps_rot"], "t_per_rotation": p["t_rot"],
        "c_t": c_t, "c_t_quoted": float(a.c_t_2028_quoted.lo),
        "site_steps": site_steps, "t_raw_per_shot": raw,
        "ct_reduction_banked": red, "t_banked_range": banked,
        "t_representative": t_point, "t_box_quoted": box_t,
        "t_box_within_banked_range": banked[0] <= box_t <= banked[1],
        "t_over_cap": t_point / cap, "t_over_cap_range": (banked[0] / cap, banked[1] / cap),
        "epsilon_l_required": eps_l, "epsilon_l_quoted": float(a.eps_l_2028_quoted.lo),
        "n_trot": n_trot, "t_window_over_mphi": t_window,
        "t_window_within_prose": _within(t_window, a.t_window_2028),
        "n_trot_from_prose_window": steps_from_window,
        "shots": shots, "shots_within_n_rare_range": _within(shots, a.n_rare_range),
        "wall_min": wall_s / 60.0, "wall_min_raw": wall_raw_s / 60.0, "shot_overhead_s": t0,
        "shot_overhead_fraction": t0 / (t_point * gt + t0),
        "gate_time_min": gates_only_s / 60.0, "gate_time_min_raw": gate_time_raw_s / 60.0,
        "family_8x8_nphi4_lq": fam4_lq, "family_8x8_nphi4_c_t": fam4["c_t"],
        "family_8x8_nphi4_t_raw": fam4["t_raw"], "family_8x8_nphi4_t_banked": fam4["t_raw"] / red,
        "family_8x8_nphi16_lq": fam16_lq, "family_8x8_nphi16_c_t": fam16["c_t"],
        "family_8x8_nphi16_t_raw": fam16["t_raw"], "family_8x8_nphi16_t_banked": fam16["t_raw"] / red,
        "with_fermions_6x6_c_t": wf["c_t"], "with_fermions_6x6_rot_per_site_step": wf["n_rot_site_step"],
        "with_fermions_6x6_t_raw": with_f_raw, "with_fermions_6x6_t_banked_all_sectors_point": with_f_banked_all_point,
        "with_fermions_6x6_system_qubits": with_f_system,
        "with_fermions_6x6_t_raw_quoted": float(a.t_with_fermions_6x6_raw_quoted.lo),
        "with_fermions_6x6_t_quoted": float(a.t_with_fermions_6x6_quoted.lo),
        "interim_qubits_per_site": interim_qps, "interim_system_qubits": interim_system,
        "interim_lq": interim_lq, "interim_lq_quoted": int(a.lq_interim_quoted.lo),
        "interim_lq_with_2028_ancilla": interim_lq_with_trotter_anc,
        "interim_two_register_lq_quoted": int(a.lq_two_register_quoted.lo),
        "interim_two_register_lq_candidates": two_register_candidates,
        "interim_phase_register_bits": phase_register_bits,
        "interim_workspace": int(a.workspace_interim.lo), "interim_workspace_from_prior": workspace_from_prior,
        "interim_site_steps": ip["site_steps"], "interim_rot_per_site_step": ip["n_rot_site_step"],
        "interim_c_t": ip["c_t"], "interim_c_t_scalar": ip["c_t_scalar"], "interim_c_t_fermion": ip["c_t_fermion"],
        "interim_t_from_inputs": interim_t, "interim_t_quoted": float(a.t_interim_quoted.lo),
        "interim_t_banked_all_sectors": interim_t_banked_all,
        "interim_t_banked_quoted": float(a.t_interim_banked_quoted.lo),
        "interim_fits_cap_with_banked_reduction": interim_t_banked_all <= cap,
        "interim_t_banked_scalar_only": interim_t_banked_scalar_only,
        "interim_t_scalar_share_quoted": float(a.t_interim_scalar_share_quoted.lo),
        "interim_fits_cap_scalar_share_only": interim_t_banked_scalar_only <= cap,
        "interim_scalar_share_le_benchmark": interim_t_banked_scalar_only <= t_point,
        "gap_factor_vs_cap": gap, "gap_factor_vs_cap_quoted": float(a.gap_interim_quoted.lo),
        # R25 slope run (single tier) and R9 exports
        "wall_min_single_time_r23": wall_s / 60.0, "wall_min_raw_single_time_r23": wall_raw_s / 60.0,
        "slope_early_steps": n_early, "slope_t_early": t_early, "slope_t_early_over_late": t_early / t_point,
        "slope_time_ratio": r_times, "slope_decrement_ref": delta_ref, "slope_decrement_late": delta_late,
        "slope_shot_factor": slope_shot_factor,
        "slope_shots_per_time": shots_per_time, "slope_shots_total": 2 * shots_per_time,
        "slope_both_times_in_window": a.t_window_2028.lo <= n_early * dt and t_window <= a.t_window_2028.hi,
        "n_trot_family_record": n_fam, "t_depth_scaled": t_depth,
        "bounce_action_implied_2028": bounce_action, "area_2028_over_mphi2": area_2028,
        "slope_shots_full_equiv": shots_equiv,
        "slope_wall_min": slope_wall_s / 60.0, "slope_wall_min_raw": slope_wall_raw_s / 60.0,
        "slope_wall_over_single_time": slope_wall_s / wall_s,
        **{k: dx[k] for k in ("t_per_shot", "t_depth_per_shot", "f_star", "floor_wall_s", "factories_for_1yr")},
        "baseline_ok": dx["baseline_ok"],
        "shot_s": t_point * gt + t0,                         # 'Maximum circuit time 0.16 s'
        "wall_first_result_s": None,                         # one tier: the benchmark point is the first result
        "wall_campaign_s": slope_wall_s,
    }
    notes = (
        f"c_T derived here: 55 R_Z + 64 T per site-step for the N_phi=16 scalar; N_rot = 55 x 432 = "
        f"{p['n_rot_shot']}; eps_rot = sqrt(1e-2/N_rot) = {p['eps_rot']:.3e}; {p['t_rot']:.3f} T/rotation; "
        f"c_T = {c_t:.1f}; raw {raw:.4g}; / banked 3.3 = {t_point:.4g} printed '2.8e5' (20 steps; was 1.6e5 at 12).",
        "The banked 3.3x is the uncited development target; each primitive's t_each is divided by it.",
        f"0.1 expected faults per shot gives epsilon_l = {eps_l:.2e}; the box prints '<~4e-7'.",
        "State prep (adiabatic ramp) and the basin projector are not priced; carried at 0 T UNSOURCED.",
        "The 2028 run drops the fermion arm for resolution: N_phi=4 cannot resolve the biased well and the fermion "
        "arm enters at N_phi=16 (288 system qubits, c_T 1973, 8.5e5 raw, 2.6e5 with 3.3x), which waits for 2033.",
    )
    notes = notes + (
        f"Slope: P_FV at {n_early} and {n_trot} steps (both in the 1-2/m_phi window), {shots_per_time:.3g} shots each "
        f"({t_early:.4g} and {t_point:.4g} T; decrement {delta_late:.3g} at the late time, {slope_shot_factor:.1f}x the "
        f"single-time count); {slope_wall_s / 3600:.1f} h per point ({slope_wall_raw_s / 3600:.1f} h raw). "
        f"T-depth {t_depth[0]:.3g}-{t_depth[1]:.3g} (factory.json x 20/12), F* {dx['f_star'][0]:.1f}-{dx['f_star'][1]:.1f}; if the 3.3x "
        "cut T-count but not depth, F* = 7-9 and the depth floor would be ~1.1-1.4x the quoted wall (not taken).",
    )
    return Result(era="2028", lq=(lq, lq), hard_ops=(t_point, t_point), breakdown=breakdown,
                  intermediates=inter, shots=(2 * shots_per_time, 2 * shots_per_time),
                  wall_time_s=(slope_wall_s, slope_wall_s),
                  epsilon_l=(eps_l, eps_l), notes=notes)


# --------------------------------------------------------------------------- #
# MLAE estimator constant (R15, derived here; see module docstring)
# --------------------------------------------------------------------------- #

def mlae_depths(m_max: int) -> list[int]:
    """Suzuki LIS depths m = 2k+1 <= m_max (base-circuit applications per circuit)."""
    if m_max < 1:
        raise ValueError("m_max must be >= 1")
    return list(range(1, m_max + 1, 2))


def balanced_branch_bound(m: int) -> float:
    """Largest |Delta| that keeps m theta in one monotone branch of sin^2(m theta) for the balanced encoding
    a = (1 + Delta)/2, theta = pi/4 + arcsin(Delta)/2 (m odd: m pi/4 is a branch center, half-width pi/4)."""
    if m % 2 == 0:
        raise ValueError("odd depth required")
    return math.sin(math.pi / (2 * m))


def lis_query_constant(m_max: int) -> float:
    """Queries of LIS over queries of all-at-m_max at equal Cramer-Rao variance: m_max sum m / sum m^2."""
    ms = mlae_depths(m_max)
    return ms[-1] * sum(ms) / sum(x * x for x in ms)


def lis_circuit_constant(m_max: int) -> float:
    """Circuits of LIS over circuits of all-at-m_max at equal variance: (K+1) m_max^2 / sum m^2."""
    ms = mlae_depths(m_max)
    return len(ms) * ms[-1] ** 2 / sum(x * x for x in ms)


def lis_branch_margin(m_max: int, n_per_depth: float) -> float:
    """Smallest ratio, over the LIS ladder, of the next depth's branch half-width pi/(4 m_next) to the theta
    standard deviation sqrt(1/(4 N sum_{m <= current} m^2)) from the depths already run."""
    ms = mlae_depths(m_max)
    out = math.inf
    for k in range(len(ms) - 1):
        sig = 1.0 / (2.0 * math.sqrt(n_per_depth * sum(x * x for x in ms[:k + 1])))
        out = min(out, (math.pi / (4 * ms[k + 1])) / sig)
    return out


COUNTER_TOF_PER_INPUT = 10   # controlled increment of a ceil(log2(n+1)) = 9-bit counter, <= 10 Toffolis
GIVENS_ROTATIONS = 2         # a complex Givens rotation = two synthesized rotations (the Ch. 7 convention)


def mlae_reflection_t(n_modes: int, lq: int, tof_t: float) -> float:
    """T per Grover iterate for the two reflections (B2, derived here; executed in the deepest circuit).
    S_chi: incident-side charge counted into a 9-bit register by n_modes controlled increments of <= 10 Toffolis
    each (~20 ancillas, reused between applications of A), computed and uncomputed; the threshold compare is
    negligible. The adder tree of arxiv_1902_10673 App. A would need n - 1 ~ 476 ancillas, which the 40-ancilla
    register does not have (verifier, 2026-10-04). S_0: controlled phase on the lq-qubit input, lq - 2 Toffolis."""
    return (2 * n_modes * COUNTER_TOF_PER_INPUT + (lq - 2)) * tof_t


def packet_rotations(n_sector_modes: int) -> int:
    """One packet per charge (particle or antiparticle) inside one transverse-momentum sector: <= n - 1 Givens
    rotations of two synthesized rotations, each controlled by the selector qubit (x2). Derived here (B2)."""
    return 2 * (n_sector_modes - 1) * GIVENS_ROTATIONS * 2


def packet_t(n_modes: int, t_rot: float) -> float:
    return packet_rotations(n_modes) * t_rot


def slater_t_max(n_modes: int, t_rot: float) -> float:
    """Superseded dense bound (referee revision 2026-10-04): <= n(n-1)/2 Givens on all n modes. Recorded only."""
    return n_modes * (n_modes - 1) / 2 * GIVENS_ROTATIONS * t_rot


def slater_rotations_ky(lz: int, lperp: int, n_f: int, comp: int) -> int:
    """Free-fermion vacuum of a z-dependent wall (B2, 2026-10-05). The periodic transverse axis makes the one-body
    Hamiltonian block-diagonal in k_y: lperp sectors of N = lz n_f comp modes at half filling eta = N/2, each a
    Slater determinant of (N - eta) eta Givens rotations (arxiv_1711_05395 Eq. n_gates)."""
    n = lz * n_f * comp
    eta = n // 2
    return lperp * (n - eta) * eta * GIVENS_ROTATIONS


def fft_transverse_t(lperp: int, chains: int) -> float:
    """Inverse fermionic FFT along the transverse axis, once per chain (lz n_f comp chains). (n/2) log2 n F_2 gates,
    each two pi/8 Pauli rotations = 2 T, and <= (n/2)(log2 n - 1) twiddle phases, powers of T (<= 1 T each).
    Rotation-free at side 8 (arxiv_1902_10673 main.tex:1102); T count derived here."""
    k = int(math.log2(lperp))
    if 2 ** k != lperp:
        raise ValueError("the transverse FFFT needs a power-of-two side")
    return chains * (lperp // 2 * k * 2 + lperp // 2 * (k - 1))


def scalar_product_rotations(v: int, n_phi: int) -> int:
    """Scalar product state: an on-site wave function about phi_cl(z) on the N_phi-point grid, binary-tree amplitude
    loading, N_phi - 1 rotations per site. The gradient coupling, the anharmonicity and the Yukawa coupling of the
    fluctuations are switched on by the dressing ramp r (not priced). Derived here (B2)."""
    return v * (n_phi - 1)


def _dims_2033(a: "Assumptions") -> dict:
    lz, lperp = int(a.lz_sites.lo), int(a.lperp_sites.lo)
    n_f, comp, n_phi = int(a.n_f_2033.lo), int(a.spinor_components.lo), int(a.n_phi.lo)
    v = lz * lperp
    lq = v * qubits_per_site(n_phi, n_f, comp) + int(a.n_anc_2033.lo)
    return dict(lz=lz, lperp=lperp, n_f=n_f, comp=comp, n_phi=n_phi, v=v, lq=lq, n_modes=v * n_f * comp,
                sector_modes=lz * n_f * comp, n_trot=int(round(float(a.t_max_2033.lo) / float(a.dt.lo))))


def free_prep(a: "Assumptions", label_bits: int = 0) -> dict:
    """Priced part of U_pk U_vac per application of A (B2, 2026-10-05): scalar product state, k_y-sector Slater
    determinant, transverse FFFT, selector-controlled packet (x 2^label_bits for the first result's label-controlled
    packets, a bound). Returns synthesized rotations and exact T."""
    d = _dims_2033(a)
    rot_scalar = scalar_product_rotations(d["v"], d["n_phi"])
    rot_slater = slater_rotations_ky(d["lz"], d["lperp"], d["n_f"], d["comp"])
    rot_packet = packet_rotations(d["sector_modes"]) * 2 ** label_bits
    t_fft = fft_transverse_t(d["lperp"], d["lz"] * d["n_f"] * d["comp"])
    return dict(rot_scalar=rot_scalar, rot_slater=rot_slater, rot_packet=rot_packet,
                rot=rot_scalar + rot_slater + rot_packet, t_exact=t_fft)


def dn_schedule(a: "Assumptions", prep: dict, r_ramp: float = 0.0) -> dict:
    """Depth-capped MLAE schedule of the Delta-n circuits (G2/B2/B3). One application of A = evolution x (1 + r_ramp)
    (r_ramp: the unpriced dressing ramp in units of the evolution, same Trotter step) + the priced free preparation.
    A depth-m circuit = m applications + (m - 1)/2 reflection pairs, its rotations synthesized to its own budget
    (R13, budget_circuits = m). The circuit cost f(m) = m t_A(m) + (m - 1)/2 refl rises with m (t_A(m) is
    nondecreasing), so the fitting depths are an interval [1, m_max]. Downward scan from the bound
    m0 = floor((env + refl/2) / (t_A(1) + refl/2)) >= m_max: m_max = the first m that fits, m_top = the largest odd
    one (m_max or m_max - 1), priced at its own budget. (Replaces a fixed-point iteration that cycled, e.g.
    33 -> 32 -> 33, for some r; verifier 2026-10-05.)"""
    if r_ramp < 0:
        raise ValueError("r is a cost ratio, >= 0")
    d = _dims_2033(a)
    tof_t = toffoli_t(1, a.toffoli_convention.value)
    refl = mlae_reflection_t(d["n_modes"], d["lq"], tof_t)
    env = float(a.envelope_2033.lo)
    n_rot_evo = price_instance(a, d["v"], d["n_trot"], d["n_phi"], d["n_f"], d["comp"])["n_rot_shot"]
    extra = prep["rot"] + int(round(r_ramp * n_rot_evo))

    def at(m: int) -> tuple:
        p = price_instance(a, d["v"], d["n_trot"], d["n_phi"], d["n_f"], d["comp"], budget_circuits=m, extra_rot=extra)
        t_evo = p["t_raw"]
        t_prep = prep["rot"] * p["t_rot"] + prep["t_exact"]
        t_a = (1 + r_ramp) * t_evo + t_prep
        return p, t_evo, t_prep, t_a, m * t_a + (m - 1) / 2 * refl

    t_a1 = at(1)[3]
    m0 = math.floor((env + refl / 2) / (t_a1 + refl / 2))
    if m0 < 1:
        raise ValueError(f"one application of A ({t_a1:.3g} T) exceeds the {env:.3g} envelope; "
                         "no depth-capped MLAE schedule fits")
    iterations, m_max = [], None
    for m in range(m0, 0, -1):
        cost = at(m)[4]
        iterations.append((m, cost, cost <= env))
        if cost <= env:
            m_max = m
            break
    m_top = m_max if m_max % 2 else m_max - 1
    p, t_evo, t_prep, t_a, deepest = at(m_top)
    return dict(m_max=m_max, m_top=m_top, m0=m0, p=p, t_evo=t_evo, t_prep=t_prep, t_a=t_a, refl=refl,
                deepest=deepest, iterations=tuple(iterations), r_ramp=r_ramp)


def circuit_s(a: "Assumptions", s: dict) -> float:
    """Serial time of one depth-m_top circuit: T x 1 us + t0."""
    return s["deepest"] * float(a.t_gate_s.lo) + float(a.shot_overhead_s.lo)


def _one_thread():
    """BLAS thread limit for the small dense algebra below (OpenBLAS oversubscription stalls a busy host)."""
    try:
        from threadpoolctl import threadpool_limits
        return threadpool_limits(1)
    except ImportError:            # pragma: no cover
        import contextlib
        return contextlib.nullcontext()


def _wilson_h_ky(a: "Assumptions", am: float, ky: float):
    """One-body Hamiltonian of the chapter's two Wilson flavors (app05:35; model docstring) at transverse momentum
    k_y on the cylinder open along z (L_z sites x n_f x 2 components), no wall: the symmetric phase."""
    import numpy as np
    d = _dims_2033(a)
    lz, nf = d["lz"], d["n_f"]
    s1 = np.array([[0, 1], [1, 0]], complex)
    s2 = np.array([[0, -1j], [1j, 0]])
    s3 = np.diag([1.0, -1.0]).astype(complex)
    rw = float(a.wilson_r.lo)
    hop = lambda alpha: -0.5 * (rw * s3 - 1j * alpha)
    onsite = s3 * (am + 2 * rw) + hop(s2) * np.exp(1j * ky) + hop(s2).conj().T * np.exp(-1j * ky)
    h = np.zeros((lz * nf * 2,) * 2, complex)
    for z in range(lz):
        for f in range(nf):
            i = (z * nf + f) * 2
            h[i:i + 2, i:i + 2] += onsite
            if z + 1 < lz:
                j = ((z + 1) * nf + f) * 2
                h[i:i + 2, j:j + 2] += hop(s1)
                h[j:j + 2, i:i + 2] += hop(s1).conj().T
    return h


def flag_contrast(a: "Assumptions", am_f: float | None = None, wall: float = 0.0, dtheta: float = 0.0) -> dict:
    """B2 flag contrast (derived here, 2026-10-05). The good flag thresholds the incident-side charge Q_inc = packet
    + vacuum fluctuation q: particle branch [Q_inc >= 1] (reflected), antiparticle branch [Q_inc >= 0] (transmitted).
    Then a = 1/2 + delta/2 + p0 A_R / 2, p0 = P(q = 0), delta = P(q >= 1) - P(q <= -1) (= 0 by C', c_prime_residual).
    q follows the full counting statistics of the free vacuum: independent Bernoulli(nu_i), nu_i the eigenvalues of
    the correlation matrix restricted to the incident half (z < L_z/2), computed sector by sector in k_y (the half is
    y-uniform). Two flavors, 2-component Wilson, r = 1, 15 x 8 cylinder open along z. Queries scale as 1/p0^2.
    wall > 0: the vacuum of the wall background U_vac actually prepares (|y| phi(z) = wall x tanh profile, theta_C
    = dtheta x profile; real-space _wilson_h, incident half z < L_z/2); wall = 0 reproduces the free value."""
    import numpy as np
    am = float(a.am_f.lo) if am_f is None else am_f
    d = _dims_2033(a)
    n_half = (d["lz"] // 2) * d["n_f"] * d["comp"]
    nus = []
    with _one_thread():
        if wall:
            e, v = np.linalg.eigh(_wilson_h(a, am, wall, dtheta))
            occ = v[:n_half * d["lperp"], e < 0]
            nus.extend(np.clip(np.linalg.eigvalsh(occ.conj() @ occ.T), 0.0, 1.0))
        else:
            for n in range(d["lperp"]):
                e, v = np.linalg.eigh(_wilson_h_ky(a, am, 2 * math.pi * n / d["lperp"]))
                occ = v[:n_half, e < 0]
                nus.extend(np.clip(np.linalg.eigvalsh(occ.conj() @ occ.T), 0.0, 1.0))
    dist = np.array([1.0])
    for x in nus:
        dist = np.convolve(dist, [1.0 - x, x])
    q = np.arange(len(dist)) - len(nus) // 2
    mean = float((dist * q).sum())
    return dict(am_f=am, p0=float(dist[q == 0][0]), var=float((dist * q * q).sum() - mean ** 2), mean=mean,
                delta=float(dist[q >= 1].sum() - dist[q <= -1].sum()))


def _wilson_h(a: "Assumptions", am: float, wall: float = 0.0, dtheta: float = 0.0, width: float = 1.0):
    """Real-space one-body Hamiltonian, 15 x 8 cylinder open along z, periodic along y, with a tanh wall
    |y| phi(z) = wall x prof(z) and theta_C(z) = dtheta x prof(z) in |y| phi [e^{i theta} psi_L^dag beta psi_R + h.c.]."""
    import numpy as np
    d = _dims_2033(a)
    lz, ly, nf = d["lz"], d["lperp"], d["n_f"]
    s1 = np.array([[0, 1], [1, 0]], complex)
    s2 = np.array([[0, -1j], [1j, 0]])
    s3 = np.diag([1.0, -1.0]).astype(complex)
    rw = float(a.wilson_r.lo)
    n = lz * ly * nf * 2
    h = np.zeros((n, n), complex)
    idx = lambda z, y, f: ((z * ly + y) * nf + f) * 2
    prof = lambda z: 0.5 * (1 + math.tanh((z - (lz - 1) / 2) / width))
    for z in range(lz):
        for y in range(ly):
            for f in range(nf):
                i = idx(z, y, f)
                h[i:i + 2, i:i + 2] += s3 * (am + 2 * rw)
                for dz, dy, alpha in ((1, 0, s1), (0, 1, s2)):
                    z2, y2 = z + dz, (y + dy) % ly
                    if z2 >= lz:
                        continue
                    j = idx(z2, y2, f)
                    hop = -0.5 * (rw * s3 - 1j * alpha)
                    h[i:i + 2, j:j + 2] += hop
                    h[j:j + 2, i:i + 2] += hop.conj().T
            if nf == 2 and wall:
                i, j = idx(z, y, 0), idx(z, y, 1)
                blk = wall * prof(z) * np.exp(1j * dtheta * prof(z)) * s3
                h[i:i + 2, j:j + 2] += blk
                h[j:j + 2, i:i + 2] += blk.conj().T
    return h


def c_prime_residual(a: "Assumptions", wall: float = 0.5, dtheta: float = 0.1) -> float:
    """max |W h^* W^dag + h| for W = sigma_1 (flavor) x sigma_1 (spinor) on every site, on a C-violating tanh wall.
    Zero means the Fock-space unitary C' (c -> W c^dagger) commutes with H: it maps a flavor-L particle onto a
    flavor-R antiparticle with the same position density at all times, so the flavor-summed Delta n vanishes and the
    flag background delta = 0. It holds term by term (mass, each hop, Yukawa at any real phi), so for the interacting
    H and every Trotter factor (derived here). W is a site-local permutation: (flavor, spin) -> (1-flavor, 1-spin)."""
    import numpy as np
    h = _wilson_h(a, float(a.am_f.lo), wall, dtheta)
    perm = np.arange(h.shape[0]) ^ 3      # flip the flavor and spin bits of the local index ((site*2+f)*2+s)
    return float(np.abs(h.conj()[np.ix_(perm, perm)] + h).max())


def filtered_jump_t(a: "Assumptions", c_t: float) -> float:
    """B5 (derived here from the arxiv_2311_09207 cost model): one unit of Lindblad time = controlled evolution of a
    patch of radius v_LR beta for a time ~beta. T >= patch sites x (beta/dt) steps x c_T (uncontrolled, polylogs dropped:
    a lower-end estimate)."""
    rr = int(a.lr_patch_radius.lo)
    patch = min((2 * rr + 1), int(a.lz_sites.lo)) * min((2 * rr + 1), int(a.lperp_sites.lo))
    steps = float(a.beta_mphi.lo) / float(a.dt.lo)
    return patch * steps * c_t


def prep_sensitivity(a: "Assumptions", r: float) -> dict:
    """B2/B3: the Delta-n schedule and walls when the unpriced dressing ramp costs r x the evolution (on top of the
    priced free preparation). Queries per bin 1/(p0^2 eps_C^2 m_top) (C_est = 1 holds a fortiori at shallower odd
    depth); walls serial, circuit by circuit."""
    p0 = flag_contrast(a)["p0"]
    s = dn_schedule(a, free_prep(a), r)
    s1 = dn_schedule(a, free_prep(a, int(a.label_register_first.lo)), r)
    eps = float(a.eps_c.lo)
    q_bin = 1.0 / (p0 ** 2 * eps ** 2 * s["m_top"])
    circ = q_bin * int(a.n_bin.lo) / s["m_top"]
    first_q = float(a.detection_sigma.lo) ** 2 / (p0 ** 2 * float(a.s_int.lo) ** 2 * s1["m_top"])
    return {"r": r, "m_max": s["m_max"], "m_top": s["m_top"], "t_A": s["t_a"], "deepest_t": s["deepest"],
            "queries_per_bin": q_bin, "dn_point_yr": circ * circuit_s(a, s) / SEC_PER_YR,
            "first_m_top": s1["m_top"], "first_dn_point_days": first_q / s1["m_top"] * circuit_s(a, s1) / SEC_PER_DAY}


# 2033
# --------------------------------------------------------------------------- #

def _model_2033(a: Assumptions) -> Result:
    v = int(a.lz_sites.lo) * int(a.lperp_sites.lo)
    n_f = int(a.n_f_2033.lo)
    comp = int(a.spinor_components.lo)
    n_phi = int(a.n_phi.lo)
    qps = qubits_per_site(n_phi, n_f, comp)                  # 8
    system = v * qps                                          # 960
    n_anc = int(a.n_anc_2033.lo)
    lq = system + n_anc                                       # 1000
    workspace = n_anc - int(a.n_bath_anc_2033.lo)             # 10
    four_comp_system = v * qubits_per_site(n_phi, n_f, 4)     # 1440

    n_trot = int(round(float(a.t_max_2033.lo) / float(a.dt.lo)))  # 100
    envelope = float(a.envelope_2033.lo)
    # R13: rotations synthesized to the deepest MLAE circuit's budget; m_max by a downward scan (dn_schedule).
    # B2 (2026-10-05): one application of A = evolution + the priced free preparation (scalar product state,
    # k_y-sector Slater determinant, transverse FFFT, selector-controlled packet); a depth-m circuit also carries
    # (m - 1)/2 reflection pairs. The interacting dressing ramp r is not priced (prep_sensitivity).
    p_base_budget = price_instance(a, v, n_trot, n_phi, n_f, comp)   # m = 1: 2.656e7 (the nucleation circuit)
    prep = free_prep(a)
    sched = dn_schedule(a, prep)
    iterations = sched["iterations"]
    m_max = sched["m_max"]                                     # 33
    p = sched["p"]
    c_t = p["c_t"]                                             # 2466
    t_per_step = c_t * v                                       # 2.959e5
    t_per_shot = p["t_raw"]                                    # evolution 2.959e7
    t_a = sched["t_a"]                                         # one application of A, 3.0e7
    contrast = flag_contrast(a)                                # p0 0.658 at a m_f = 0.5
    p0 = contrast["p0"]
    # verifier 2026-10-05: U_vac prepares the wall-background vacuum, where p0 is lower; the free p0 makes every
    # Delta-n query count a lower bound. |y| phi_wall = 0.3, 0.5, 1.0 (<= 2 m_f) at the chapter's Delta theta_C = 0.1.
    p0_wall = {w: flag_contrast(a, wall=w, dtheta=0.1)["p0"] for w in (0.3, 0.5, 1.0)}   # 0.654, 0.646, 0.603
    eps = float(a.eps_c.lo)
    per_bin_shot_noise = 1.0 / eps ** 2
    per_bin_heisenberg = 1.0 / eps
    # R16: estimator constant. Single depth (constant 1) when a derived bound or the author-adopted prior keeps
    # m_max theta in one branch (R16: prior |A_R| < sin(pi/66)); otherwise Suzuki LIS 1, 3, ..., m_max.
    m_top = m_max if m_max % 2 else m_max - 1                  # odd deepest depth (33)
    branch_bound = balanced_branch_bound(m_top)                # 0.0476
    derived_bound_ok = float(a.dn_state_independent_bound.lo) < branch_bound   # False (R15)
    prior_ok = (a.mlae_schedule.value == "single" and float(a.c_est.lo) == 1.0
                and round(float(a.dn_branch_prior.lo), 3) == round(branch_bound, 3)
                and float(a.dn_estimate.lo) < branch_bound)  # True (R16)
    single_depth = derived_bound_ok or prior_ok
    if a.mlae_schedule.value not in ("single", "LIS"):
        raise ValueError("only the single-depth and Suzuki LIS schedules are priced")
    c_lis = lis_query_constant(m_top)                          # 1.45714 (sensitivity)
    c_est = 1.0 if single_depth else c_lis                     # 1 (R16)
    per_bin_single_depth = 1.0 / (eps ** 2 * m_max)            # 3.03e4 (constant 1, contrast 1: r25 record)
    per_bin_ideal = c_est / (eps ** 2 * m_top)                 # 3.03e4 at unit flag contrast (r26 record)
    per_bin_mlae = c_est / (p0 ** 2 * eps ** 2 * m_top)        # 7.0e4: B2 flag contrast, queries x 1/p0^2
    per_bin_lis = c_lis / (p0 ** 2 * eps ** 2 * m_top)         # (sensitivity)
    depths = [m_top] if single_depth else mlae_depths(m_top)
    n_per_depth = per_bin_mlae / sum(depths)                   # 918.3 circuits per depth per bin
    circuits_per_bin = n_per_depth * len(depths)               # 918.3
    lis_n_per_depth = per_bin_lis / sum(mlae_depths(m_top))    # 152.8 (sensitivity)
    n_bin = int(a.n_bin.lo)                                    # 8 distinct |k_z| (ruling 2026-10-11; was 100)
    n_bin_distinct_kz = int(a.lz_sites.lo) // 2 + 1            # 8
    queries = per_bin_mlae * n_bin                             # 5.6e5 (was 7.0e6 at 100 bins)
    queries_heisenberg = per_bin_heisenberg * n_bin
    overhead = per_bin_mlae / per_bin_heisenberg               # 30.3

    # E27 (r23): every wall time is serial on one machine, shots x per-shot time at the report convention
    # (1 us per T-gate plus t0 = 0.1 ms per circuit). No machine count is formed anywhere.
    t0 = float(a.shot_overhead_s.lo)
    horizon = float(a.horizon_yr.lo)
    base_circuit_s = t_a * float(a.t_gate_s.lo)            # 30.0 s, one application of A
    n_circuits = circuits_per_bin * n_bin                      # 2.1e5 circuits of depth m_top
    deepest_t_pre = sched["deepest"]
    dn_arm_s = n_circuits * circuit_s(a, sched)                # 2.1e8 s = 6.6 yr per campaign Delta-n point
    lis_circuits = lis_n_per_depth * len(mlae_depths(m_top)) * n_bin
    lis_dn_arm_s = per_bin_lis * n_bin * base_circuit_s + lis_circuits * t0   # (sensitivity; reflections dropped)
    nucl_single_time_s = float(a.shots_nucleation_2033.lo) * (t_per_shot * float(a.t_gate_s.lo) + t0)   # r23 record
    # R25 (R2 + shot audit): the nucleation slope from t/3 and t, 3x the single-time count at each time. These are
    # single circuits (no MLAE), so each is synthesized to its own budget (R-TOL per circuit) at its own depth.
    gt = float(a.t_gate_s.lo)
    n_early = int(round(n_trot * float(a.slope_early_fraction_2033.lo)))      # 33 steps
    t_nucl_late = p_base_budget["t_raw"]                                       # 2.6555e7 (own budget)
    t_nucl_early = price_instance(a, v, n_early, n_phi, n_f, comp)["t_raw"]    # 8.446e6
    acc = float(a.accuracy.lo)
    delta_2033 = 1.0 / (float(a.shots_nucleation_2033.lo) * acc ** 2)         # 2.5e-3 at t = 10/m_phi (box 1e4 shots)
    nucl_shots_per_time = slope_shots(acc, delta_2033, n_early / n_trot)      # 2.96e4 -> '3e4'
    nucl_s = nucl_shots_per_time * ((t_nucl_late + t_nucl_early) * gt + 2 * t0)   # 1.04e6 s = 12.0 d
    # Finite-volume saddle (ruling 2026-10-11, words in the chapter): on 15 x 8 no critical bubble fits and the
    # saddle is a wall pair wrapping L_perp; the decrement 2.5e-3 corresponds to sigma = 1.32 m_phi^2, Gamma/V is finite
    # at T_c, and at 0.9 T_c the action over T grows so the decrement roughly halves (2x the shots at that point).
    slab_tc = slab_saddle(a, 1.0, float(a.t_max_2033.lo))
    slab_09 = slab_saddle(a, 0.9, float(a.t_max_2033.lo))
    decrement_ratio_09 = slab_09["decrement"] / slab_tc["decrement"]          # 0.50
    nucl_shots_per_time_09 = nucl_shots_per_time / decrement_ratio_09         # 5.9e4
    nucl_s_09 = nucl_s / decrement_ratio_09                                   # 24 d
    n_dn = int(a.n_dn_points.lo)
    dn_point_yr = dn_arm_s / SEC_PER_YR                        # 2.84
    dn_campaign_yr = n_dn * dn_point_yr                        # 28.4 yr serial
    n_scan = float(a.n_scan_points.lo)
    nucl_scan_yr = n_scan * nucl_s / SEC_PER_YR                # 3.33 yr serial
    campaign_serial_yr = dn_campaign_yr + nucl_scan_yr         # 31.7 yr serial on one machine
    campaign_over_horizon = campaign_serial_yr / horizon       # 6.35: cost reduction needed to fit the horizon

    # R25 first result (R1/R2): 3-sigma detection of the f_eq-weighted Delta n at the two box temperatures, plus
    # the nucleation slope at both. Queries (1 - S^2)/(sigma^2 m) with sigma = S/3; S_int not computed (1e-3 taken).
    s_int = float(a.s_int.lo)
    n_sig = float(a.detection_sigma.lo)
    n_first = int(a.n_points_first.lo)
    # B2: the label-controlled packets (2^7 x one packet, a bound) lengthen A, so these circuits run shallower.
    sched_first = dn_schedule(a, free_prep(a, int(a.label_register_first.lo)))
    m_first = sched_first["m_top"]                             # 31
    first_queries = c_est * n_sig ** 2 / (p0 ** 2 * s_int ** 2 * m_first)   # 6.7e5 per point
    first_circuits = first_queries / m_first                   # 2.2e4 per point
    first_dn_s = first_circuits * circuit_s(a, sched_first)    # 2.2e7 s = 2.5e2 d
    first_s = n_first * first_dn_s + nucl_s + nucl_s_09        # 4.6e7 s = 528 d = 1.45 yr (T_c and 0.9 T_c slopes)
    # what fits in the horizon on one machine: the nucleation scan (which contains the two first-result
    # temperatures) plus the first-result Delta-n points, then whole campaign Delta-n points
    scan_plus_first_yr = nucl_scan_yr + n_first * first_dn_s / SEC_PER_YR   # 3.84 yr
    dn_points_in_horizon = min(n_dn, math.floor((horizon - nucl_scan_yr) / dn_point_yr))   # 0
    fits_in_horizon_yr = scan_plus_first_yr
    lis_campaign_yr = n_dn * lis_dn_arm_s / SEC_PER_YR + nucl_scan_yr   # 44.7 yr serial (sensitivity)

    # R9 exports. Campaign: Delta-n base circuits plus nucleation circuits as full-circuit equivalents.
    # (equivalents in units of the evolution, whose T-depth factory.json gives; the preparation and reflections add
    # 1.6% and 0.25% and are carried in the T total, at the same T per layer)
    campaign_equiv = (n_dn * n_circuits * deepest_t_pre + n_scan * nucl_shots_per_time * (t_nucl_late + t_nucl_early)) / t_per_shot
    first_equiv = (n_first * first_circuits * sched_first["deepest"]
                   + (nucl_shots_per_time + nucl_shots_per_time_09) * (t_nucl_late + t_nucl_early)) / t_per_shot
    dx = depth_exports(t_per_shot, a.t_depth_2033, campaign_equiv, a.t_gate_s, a.shot_overhead_s)
    # The first-result Delta-n circuits carry the 7-qubit label, so they use t_depth_2033_first (A = 33 borrowed to
    # A = 3 not). The nucleation circuits in first_equiv carry no label; pricing them at this band too is conservative.
    dx_first = depth_exports(t_per_shot, a.t_depth_2033_first, first_equiv, a.t_gate_s, a.shot_overhead_s)
    # Without the borrowing (A = 3, F* ~ 3 < 10) each first-result base circuit is depth-bound at D x t_r.
    first_base_unborrowed_s = max(t_per_shot * gt, float(a.t_depth_2033_first.hi) * REACTION_TIME_S)   # 95.7 s per evolution
    first_dn_unborrowed_s = (first_circuits * sched_first["deepest"] / t_per_shot * first_base_unborrowed_s
                             + first_circuits * t0)   # 2.61e7 s = 302 d
    first_unborrowed_yr = (n_first * first_dn_unborrowed_s + nucl_s + nucl_s_09) / SEC_PER_YR   # 4.4 yr
    dx_dn_point = depth_exports(t_per_shot, a.t_depth_2033, n_circuits * deepest_t_pre / t_per_shot, a.t_gate_s,
                                a.shot_overhead_s)

    t_if_banked_point = t_per_shot / float(a.ct_reduction_banked.lo)   # not banked (R11 A); record only

    fb = float(a.fault_budget_per_shot.lo)
    deepest_t = sched["deepest"]                                         # 9.9e8 (33 x A + 16 reflection pairs)
    eps_l_base = eps_l_required(t_per_shot, fb)                          # 3.38e-9
    eps_l_deepest = eps_l_required(deepest_t, fb)                        # 1.02e-10

    # B2 preparation, priced per application of A (derived here, 2026-10-05; was 'bounded, not added')
    n_modes = v * n_f * comp                                   # 480
    sector_modes = int(a.lz_sites.lo) * n_f * comp             # 60 per k_y sector
    t_rot = p["t_rot"]
    refl_t = sched["refl"]                                     # 7.42e4 per Grover iterate
    pk_t = prep["rot_packet"] * t_rot                          # 1.3e4 (one k_y sector)
    slater_t = prep["rot_slater"] * t_rot                      # 3.9e5 (8 sectors x 900 Givens)
    scalar_prep_t = prep["rot_scalar"] * t_rot                 # 4.9e4
    fft_t = prep["t_exact"]                                    # 1.9e3
    slater_dense_t = slater_t_max(n_modes, t_rot)              # 6.3e6, superseded r26 bound (record)
    label_pk_t = 2 ** int(a.label_register_first.lo) * pk_t    # 1.7e6, first-result circuits only
    r_free = sched["t_prep"] / t_per_shot                      # 0.0155
    sens_r1 = prep_sensitivity(a, 1.0)
    sens_r02 = prep_sensitivity(a, 0.2)
    ramp_steps_at_depth33 = (envelope - deepest_t) / m_top / (t_per_shot / n_trot)   # < 1: any ramp lowers the depth
    c_prime = c_prime_residual(a)
    # G4 (derived here): fault model at the deepest circuit. Undetected logical faults Poisson, mean lam; a faulty
    # circuit returns its flag with an unknown probability q_f. Then a_obs = (1 - p_f) a + p_f q_f, so the estimate of
    # A_R is scaled by (1 - p_f) and shifted by 2 p_f (q_f - 1/2) / (m p0) <= p_f / (m p0) (worst case q_f in {0, 1}).
    p_fault = 1.0 - math.exp(-fb)                              # 0.095 at the R3 budget
    bias_worst = p_fault / (m_top * p0)                        # 4.4e-3 in A_R units
    sigma_bin = eps                                            # per-bin target
    sigma_first = s_int / n_sig                                # 3.3e-4
    p_fault_bin = m_top * p0 * sigma_bin                       # 0.022
    p_fault_first = m_first * p0 * sigma_first                 # 6.8e-3
    eps_l_bias_bin = -math.log(1 - p_fault_bin) / deepest_t
    eps_l_bias_first = -math.log(1 - p_fault_first) / sched_first["deepest"]
    qubit_cycles = lq * deepest_t * gt / float(a.idle_cycle_s.lo)   # 9.9e10 idle/Clifford locations
    lam_idle_at_box = qubit_cycles * eps_l_required(deepest_t, fb)      # ~10 faults if idles fail at the T rate
    eps_l_with_idle = fb / (deepest_t + qubit_cycles)          # 1.0e-12 per location for 0.1 faults
    qubit_cycles_first = lq * sched_first["deepest"] * gt / float(a.idle_cycle_s.lo)   # of the depth-31 circuit
    eps_l_with_idle_first = -math.log(1 - p_fault_first) / (sched_first["deepest"] + qubit_cycles_first)   # 7e-14
    null_control_factor = 4.0                                  # equal-split difference of two equal-precision estimates
    # verifier 2026-10-05: the null route costed. It removes only the bias common to the signal and null circuits
    # (different circuits). Delta-n queries x 4 at each first-result point; nucleation unchanged.
    first_null_yr = (n_first * null_control_factor * first_dn_s + nucl_s + nucl_s_09) / SEC_PER_YR   # ~5.5 yr
    scan_plus_first_null_yr = nucl_scan_yr + n_first * null_control_factor * first_dn_s / SEC_PER_YR   # ~8.7 yr
    # B5 (derived here from the arxiv_2311_09207 cost model): one unit of Lindblad time = one filtered jump on a
    # patch; ||sum A+ A|| <= 1 with adjoints included means a sweep over >= 2V scalar jumps (a_x and a_x^dag) is 2V
    # units. Nucleation circuits are single circuits at their own budget (c_T at m = 1).
    jump_t = filtered_jump_t(a, p_base_budget["c_t"])          # 5.5e5
    n_jumps_min = 2 * v                                        # 240
    sweep_t = n_jumps_min * jump_t                             # 1.3e8
    sweeps_in_envelope = (envelope - t_nucl_late) / sweep_t    # 7.3
    # Mixing-time band (words only in the chapter; ruling 2026-10-11): if every nucleation shot paid n sweeps of
    # Gibbs preparation (each shot ends in a projective measurement), the shot, the point, the scan and the first
    # result would be:
    def _with_prep(n_sweeps: float) -> dict:
        t_prep = n_sweeps * sweep_t
        shot_late = t_nucl_late + t_prep
        point_s = nucl_shots_per_time * ((shot_late + t_nucl_early + t_prep) * gt + 2 * t0)
        return {"sweeps": n_sweeps, "prep_t": t_prep, "nucleation_shot_t": shot_late,
                "fits_envelope": shot_late <= envelope,
                "nucleation_point_days": point_s / SEC_PER_DAY,
                "nucleation_scan_yr": n_scan * point_s / SEC_PER_YR,
                "first_result_yr": (n_first * first_dn_s + point_s * (1 + 1 / decrement_ratio_09)) / SEC_PER_YR,
                "campaign_yr": dn_campaign_yr + n_scan * point_s / SEC_PER_YR}
    mixing = {k: _with_prep(n) for k, n in (("lo", float(a.bath_mixing_sweeps.lo)),
                                            ("central", float(a.bath_mixing_sweeps_central.lo)),
                                            ("hi", float(a.bath_mixing_sweeps.hi)))}
    # G2: the headline is the deepest executed circuit: m_top applications of A plus (m_top - 1)/2 reflection pairs.
    deep_prims = tuple(replace(q, count=q.count * m_top,
                               note=q.note + f"; x {m_top} applications of A (deepest MLAE circuit, G2)")
                       for q in p["primitives"])
    breakdown = deep_prims + (
        Primitive("vacuum_prep_scalar_product", prep["rot_scalar"] * m_top, t_rot, CircuitStatus.COMPILED, DERIVED,
                  f"on-site amplitude loading, {int(a.n_phi.lo) - 1} rotations per site, per application (B2)"),
        Primitive("vacuum_prep_slater_ky", prep["rot_slater"] * m_top, t_rot, CircuitStatus.COMPILED,
                  a.slater_givens.src, f"{int(a.lperp_sites.lo)} k_y sectors x (N - eta) eta = 900 Givens, N = "
                  f"{sector_modes}, 2 rotations each, per application (B2)"),
        Primitive("vacuum_prep_fft_transverse", m_top, fft_t, CircuitStatus.COMPILED, a.fft8_rotation_free.src,
                  "inverse 8-point fermionic FFT on 60 chains, rotation-free, per application (B2)"),
        Primitive("packet_prep", prep["rot_packet"] * m_top, t_rot, CircuitStatus.COMPILED, DERIVED,
                  f"selector-controlled packet in one k_y sector, <= {sector_modes - 1} Givens per charge (B2)"),
        Primitive("mlae_reflections", (m_top - 1) // 2, refl_t, CircuitStatus.COMPILED, DERIVED,
                  "S_chi (controlled-increment counter) + S_0 per Grover iterate (B2)"),
        Primitive("vacuum_prep_dressing", m_top, 0.0, CircuitStatus.UNSOURCED, f"{TEX}:83",
                  "adiabatic ramp switching on the scalar gradient, anharmonicity and Yukawa coupling of the "
                  "fluctuations about the wall; r x the evolution per application; not priced (B2)"),
        Primitive("lindblad_dissipator_step", n_trot, 0.0, CircuitStatus.UNSOURCED, f"{TEX}:55,99",
                  "nucleation circuits only: filtered dissipator each step; unpriced (B5 estimate in intermediates)"),
        Primitive("gibbs_state_prep", 1, 0.0, CircuitStatus.UNSOURCED, f"{TEX}:91,106",
                  f"nucleation circuits only: Gibbs preparation, >= {sweep_t:.2g} T per sweep (B5); basin projection; unpriced"),
    )
    inter = {
        "volume": v, "qubits_per_site": qps, "system_qubits": system, "n_anc": n_anc,
        "n_anc_within_generic_range": _within(n_anc, a.n_anc_generic),
        "workspace_anc_implied": workspace, "lq": lq, "four_component_system_qubits": four_comp_system,
        "lz_over_mphi": int(a.lz_sites.lo) * float(a.a_mphi.lo),
        "n_trot": n_trot,
        "rot_per_site_step": p["n_rot_site_step"], "rot_by_term": p["rot_by_term"],
        "t_exact_per_site_step": p["t_exact_scalar"] + p["t_exact_fermion"],
        "n_rot_shot": p["n_rot_shot"], "n_rot_budget": p["n_rot_budget"],
        "eps_rot": p["eps_rot"], "t_per_rotation": p["t_rot"],
        "m_max_iterations": tuple(iterations),
        "base_budget_eps_rot": p_base_budget["eps_rot"], "base_budget_c_t": p_base_budget["c_t"],
        "base_budget_t_per_shot": p_base_budget["t_raw"],
        "base_budget_m_max": math.floor(envelope / p_base_budget["t_raw"]),
        "c_t": c_t, "c_t_scalar": p["c_t_scalar"], "c_t_fermion": p["c_t_fermion"],
        "c_t_quoted": float(a.c_t_2033_quoted.lo),
        "t_per_step": t_per_step, "fraction_of_envelope": t_per_shot / envelope,
        "m_max": m_max, "m_max_quoted": int(a.m_max_quoted.lo),
        "per_bin_shot_noise": per_bin_shot_noise, "per_bin_heisenberg": per_bin_heisenberg,
        "per_bin_mlae": per_bin_mlae, "per_bin_mlae_quoted": float(a.mlae_queries_per_bin_quoted.lo),
        "mlae_single_depth": single_depth, "mlae_estimator_constant": c_est,
        "mlae_branch_bound_needed": branch_bound, "dn_estimate": float(a.dn_estimate.lo),
        "dn_estimate_headroom": branch_bound / float(a.dn_estimate.lo),
        "mlae_depths": tuple(depths), "mlae_circuits_per_depth_per_bin": n_per_depth,
        "mlae_circuits_per_bin": circuits_per_bin,
        "mlae_lis_circuit_constant": lis_circuit_constant(m_top),
        "mlae_lis_branch_margin_sigma": lis_branch_margin(m_top, lis_n_per_depth),
        "mlae_lis_query_constant": c_lis, "mlae_lis_per_bin": per_bin_lis, "mlae_lis_queries_total": per_bin_lis * n_bin,
        "mlae_lis_circuits_per_depth_per_bin": lis_n_per_depth,
        "mlae_lis_dn_arm_yr": lis_dn_arm_s / SEC_PER_YR, "mlae_lis_campaign_serial_yr": lis_campaign_yr,
        "mlae_lis_campaign_over_horizon": lis_campaign_yr / horizon,
        "mlae_derived_bound_ok": derived_bound_ok, "mlae_author_prior_ok": prior_ok,
        "n_bin": n_bin, "n_bin_distinct_kz": n_bin_distinct_kz,
        "per_bin_single_depth_r13": per_bin_single_depth, "queries_total_single_depth_r13": per_bin_single_depth * n_bin,
        "mlae_overhead_vs_heisenberg": overhead, "queries_heisenberg_total": queries_heisenberg,
        "queries_total": queries, "base_circuit_s": base_circuit_s,
        "dn_arm_yr": dn_arm_s / SEC_PER_YR, "nucleation_point_days": nucl_s / SEC_PER_DAY,
        "mlae_circuits_total": n_circuits, "shot_overhead_s": t0,
        "dn_campaign_serial_yr": dn_campaign_yr, "nucleation_scan_serial_yr": nucl_scan_yr,
        "campaign_serial_yr": campaign_serial_yr, "horizon_yr": horizon,
        "campaign_over_horizon": campaign_over_horizon,
        "dn_points_in_horizon": dn_points_in_horizon, "fits_in_horizon_yr": fits_in_horizon_yr,
        # R25 slope, tiers and R9 exports
        "nucleation_point_days_single_time_r23": nucl_single_time_s / SEC_PER_DAY,
        "nucleation_shots_per_time": nucl_shots_per_time, "nucleation_early_steps": n_early,
        "nucleation_decrement_late": delta_2033, "nucleation_slope_shot_factor": nucl_shots_per_time * acc ** 2 * delta_2033,
        # finite-volume saddle (slab_saddle), ruling 2026-10-11
        "slab_action_over_t_tc": slab_tc["action_over_t"], "slab_decrement_tc": slab_tc["decrement"],
        "slab_reproduces_assumed_decrement": math.isclose(slab_tc["decrement"], delta_2033, rel_tol=0.02),
        "bubble_min_action_over_t_tc": slab_tc["bubble_min_action_over_t"], "slab_dominates": slab_tc["slab_dominates"],
        "lambda_over_mphi_thin_wall": slab_tc["lambda_over_mphi_thin_wall"],
        "slab_action_over_t_0p9tc": slab_09["action_over_t"], "decrement_ratio_0p9tc": decrement_ratio_09,
        "nucleation_shots_per_time_0p9tc": nucl_shots_per_time_09, "nucleation_point_days_0p9tc": nucl_s_09 / SEC_PER_DAY,
        # bath mixing band (words only)
        "bath_mixing_sweeps": (float(a.bath_mixing_sweeps.lo), float(a.bath_mixing_sweeps.hi)),
        "bath_mixing_sensitivity": mixing,
        "bath_prep_dominates_nucleation_shot": mixing["lo"]["prep_t"] > t_nucl_late,
        "nucleation_shot_t_with_prep_band": (mixing["lo"]["nucleation_shot_t"], mixing["hi"]["nucleation_shot_t"]),
        "nucleation_scan_yr_with_prep_band": (mixing["lo"]["nucleation_scan_yr"], mixing["hi"]["nucleation_scan_yr"]),
        "first_result_yr_with_prep_band": (mixing["lo"]["first_result_yr"], mixing["hi"]["first_result_yr"]),
        "nucleation_t_late_own_budget": t_nucl_late, "nucleation_t_early": t_nucl_early,
        "s_int_taken": s_int, "first_detection_sigma": n_sig, "first_points": n_first,
        "first_queries_per_point": first_queries, "first_circuits_per_point": first_circuits,
        "first_dn_point_days": first_dn_s / SEC_PER_DAY, "first_result_days": first_s / SEC_PER_DAY,
        "first_result_yr": first_s / SEC_PER_YR, "scan_plus_first_yr": scan_plus_first_yr,
        "first_label_register": int(a.label_register_first.lo),
        "first_workspace_left": workspace - int(a.label_register_first.lo),
        "campaign_shots_full_equiv": campaign_equiv, "first_shots_full_equiv": first_equiv,
        "campaign_serial_yr_check": dx["serial_yr"],
        **{k: dx[k] for k in ("t_per_shot", "t_depth_per_shot", "f_star", "floor_wall_s", "factories_for_1yr")},
        "baseline_ok": dx["baseline_ok"],
        "first_t_depth_per_shot": dx_first["t_depth_per_shot"], "first_f_star": dx_first["f_star"],
        "first_baseline_ok": dx_first["baseline_ok"],
        "first_dn_point_days_unborrowed": first_dn_unborrowed_s / SEC_PER_DAY,
        "first_result_yr_unborrowed": first_unborrowed_yr,
        "first_floor_wall_s": dx_first["floor_wall_s"], "first_factories_for_1yr": dx_first["factories_for_1yr"],
        "dn_point_floor_wall_s": dx_dn_point["floor_wall_s"],
        "dn_point_factories_for_1yr": dx_dn_point["factories_for_1yr"],
        "wall_first_result_s": first_s,
        "wall_campaign_s": campaign_serial_yr * SEC_PER_YR,
        # retired (E27, r23): machines_over_horizon 5.87 and machines_needed 6 (ceil); the report assumes one machine
        "t_per_shot_if_reduction_banked": _banked(t_per_shot, a.ct_reduction),
        "t_per_shot_if_reduction_banked_point": t_if_banked_point,
        "deepest_mlae_circuit_t": deepest_t, "headline_t": deepest_t, "base_circuit_t": t_per_shot,
        "deepest_circuit_s": deepest_t * gt + t0,
        # B2 construction, priced (2026-10-05)
        "n_fermion_modes": n_modes, "sector_modes": sector_modes, "mlae_reflection_t_per_iterate": refl_t,
        "mlae_reflection_fraction": refl_t / t_per_shot,
        "packet_t": pk_t, "slater_t": slater_t, "slater_fraction": slater_t / t_per_shot,
        "scalar_prep_t": scalar_prep_t, "fft_t": fft_t, "free_prep_t": sched["t_prep"], "r_free": r_free,
        "slater_dense_bound_t_r26": slater_dense_t, "t_A": t_a,
        "label_packets_t": label_pk_t, "label_packets_fraction": label_pk_t / t_per_shot,
        "first_m_top": m_first, "first_t_A": sched_first["t_a"], "first_deepest_t": sched_first["deepest"],
        "flag_contrast_p0": p0, "flag_contrast_var": contrast["var"], "flag_vacuum_background": contrast["delta"],
        "query_factor_contrast": 1.0 / p0 ** 2, "per_bin_unit_contrast_r26": per_bin_ideal,
        "c_prime_residual": c_prime,
        "prep_r1_m_max": sens_r1["m_max"], "prep_r1_m_top": sens_r1["m_top"],
        "prep_r1_dn_point_yr": sens_r1["dn_point_yr"], "prep_r1_first_dn_point_days": sens_r1["first_dn_point_days"],
        "prep_r1_dn_wall_ratio": sens_r1["dn_point_yr"] / (dn_arm_s / SEC_PER_YR),
        "prep_r02_m_top": sens_r02["m_top"], "prep_r02_dn_point_yr": sens_r02["dn_point_yr"],
        "prep_r02_first_dn_point_days": sens_r02["first_dn_point_days"],
        # B3/G2: headroom left at depth 33 for the unpriced ramp, in Trotter steps per application (< 1)
        "ramp_steps_left_at_depth33": ramp_steps_at_depth33,
        "prep_r_depth33_exceeds_envelope": float(a.envelope_2033.lo) / deepest_t - 1.0,
        # G4 fault / bias model
        "fault_prob_per_circuit": p_fault, "fault_bias_worst_A_R": bias_worst, "fault_scale_bias": p_fault,
        "fault_prob_for_bias_below_bin_sigma": p_fault_bin, "fault_prob_for_bias_below_first_sigma": p_fault_first,
        "epsilon_l_bias_bin": eps_l_bias_bin, "epsilon_l_bias_first": eps_l_bias_first,
        "idle_qubit_cycles_deepest": qubit_cycles, "idle_faults_at_box_epsilon_l": lam_idle_at_box,
        "epsilon_l_with_idle": eps_l_with_idle, "epsilon_l_with_idle_first_bias": eps_l_with_idle_first,
        "null_control_query_factor": null_control_factor, "first_result_null_yr": first_null_yr,
        "scan_plus_first_null_yr": scan_plus_first_null_yr,
        "flag_contrast_p0_wall": p0_wall,
        "query_factor_contrast_wall_max": 1.0 / min(p0_wall.values()) ** 2,
        # B5 filtered-bath cost model
        "filtered_jump_t": jump_t, "bath_jumps_min": n_jumps_min, "bath_sweep_t": sweep_t,
        "bath_sweep_over_nucleation_circuit": sweep_t / t_nucl_late, "bath_sweeps_in_envelope": sweeps_in_envelope,
        "epsilon_l_required_base_circuit": eps_l_base,
        "epsilon_l_required_deepest_mlae_circuit": eps_l_deepest,
        "epsilon_l_quoted_base_circuit": float(a.eps_l_2033_base_quoted.lo),
        "epsilon_l_quoted_deepest_mlae_circuit": float(a.eps_l_2033_deepest_quoted.lo),
    }
    notes = (
        f"c_T derived here: 87 R_Z + 80 T per site-step; N_rot(base) = 87 x 12000 = {p['n_rot_shot']}; "
        f"synthesized to the deepest MLAE circuit's budget, N_rot = {m_max} x ({p['n_rot_shot']} + {prep['rot_scalar']} + {prep['rot_slater']} + {prep['rot_packet']}) = {p['n_rot_budget']}; "
        f"eps_rot = {p['eps_rot']:.3e}; {p['t_rot']:.3f} T/rotation; c_T = {c_t:.1f} "
        f"(scalar {p['c_t_scalar']:.1f}, fermion {p['c_t_fermion']:.1f}); raw, no reduction banked.",
        f"m_max downward scan from m0 = {sched['m0']} (m: m T_A(m) + (m-1)/2 refl vs 1e9): "
        + "; ".join(f"{mi}: {ci:.5g} {'fits' if ok else 'over'}" for mi, ci, ok in iterations) + ".",
        f"m_max = {m_max}; estimator constant {c_est:.5g} ({'single depth' if single_depth else 'Suzuki LIS 1..' + str(m_top)}; "
        f"|Delta| < {branch_bound:.4f} needed for a single depth: author-adopted prior; no derived bound; "
        f"LIS sensitivity {c_lis:.5g}); "
        f"{per_bin_mlae:.4g}/bin; {queries:.4g} total over {n_bin} bins (the distinct |k_z| on the 15-site axis; "
        f"ruling 2026-10-11, was 100); cited to arxiv_1904_10246 and arxiv_2012_03348.",
        "Gibbs prep and the Lindblad dissipator (nucleation circuits) and the Delta-n dressing ramp are not priced; "
        f"headline = deepest circuit {deepest_t:.4g} T, a lower bound (printed '>~9.9e8').",
        f"A = U_evo U_pk U_vac; priced free preparation {sched['t_prep']:.3g} T per application (r_free "
        f"{r_free:.4f}: Slater {slater_t:.3g}, scalar {scalar_prep_t:.3g}, FFFT {fft_t:.3g}, packet {pk_t:.3g}); "
        f"reflections {refl_t:.3g} per iterate; flag contrast p0 {p0:.3f} (queries x {1 / p0 ** 2:.3g}); C' residual "
        f"{c_prime:.1e}; at r = 1, m_max {sens_r1['m_max']} (odd {sens_r1['m_top']}), "
        f"{sens_r1['dn_point_yr']:.3g} yr per campaign Delta-n point.",
        f"Bath: filtered jump {jump_t:.3g} T per unit Lindblad time; sweep over {n_jumps_min} jumps {sweep_t:.3g} T "
        f"({sweep_t / t_nucl_late:.2g} nucleation circuits); {sweeps_in_envelope:.2g} sweeps fit 1e9. Mixing band "
        f"{a.bath_mixing_sweeps.value} sweeps (heuristic): a nucleation shot would be {mixing['lo']['nucleation_shot_t']:.2g}-"
        f"{mixing['hi']['nucleation_shot_t']:.2g} T, the scan {mixing['lo']['nucleation_scan_yr']:.0f}-"
        f"{mixing['hi']['nucleation_scan_yr']:.0f} yr, the first result {mixing['lo']['first_result_yr']:.1f}-"
        f"{mixing['hi']['first_result_yr']:.0f} yr; words only, not priced.",
        f"Saddle on 15x8 (slab_saddle): wall pair, S/T = {slab_tc['action_over_t']:.1f} at sigma = "
        f"{float(a.sigma_over_mphi2.lo)} m^2 reproduces the decrement {slab_tc['decrement']:.2e} (assumed {delta_2033:.1e}); "
        f"bubble minimum S/T {slab_tc['bubble_min_action_over_t']:.1f}; at 0.9 T_c the decrement is x{decrement_ratio_09:.2f}, "
        f"{nucl_s_09 / SEC_PER_DAY:.0f} d for that point.",
        f"Fault model: p_f {p_fault:.3f}; worst bias {bias_worst:.2g} in A_R; eps_l for bias < first-result sigma "
        f"{eps_l_bias_first:.2g} per T; idle qubit-cycles {qubit_cycles:.3g} -> {lam_idle_at_box:.3g} faults at the box "
        f"eps_l; {eps_l_with_idle:.2g} per location for 0.1 faults.",
        f"epsilon_l = {eps_l_base:.2e} base circuit, {eps_l_deepest:.2e} deepest MLAE circuit "
        f"({deepest_t:.3g} T); printed '<~3e-9' and '<~1e-10'.",
        f"Tiers: first result {first_s / SEC_PER_DAY:.0f} d (2 x 3-sigma Delta n at {first_dn_s / SEC_PER_DAY:.1f} d, "
        f"S_int = {s_int:g} taken, not computed; nucleation slopes {nucl_s / SEC_PER_DAY:.1f} d at T_c and "
        f"{nucl_s_09 / SEC_PER_DAY:.0f} d at 0.9 T_c); campaign "
        f"{campaign_serial_yr:.2f} yr, {campaign_over_horizon:.2f}x the horizon; scan + first result {scan_plus_first_yr:.2f} yr fits.",
        f"T-depth {a.t_depth_2033.value} per base circuit (factory.json), F* {dx['f_star'][0]:.1f}-{dx['f_star'][1]:.1f}: "
        "the 10-factory baseline is just fed at A = 10. The first result's ~7-qubit packet label comes from the "
        "10-qubit workspace, leaving A = 3; with the bath borrowed (A = 33) its T-depth is "
        f"{a.t_depth_2033_first.lo:.2g}, F* {dx_first['f_star'][1]:.1f}; without, F* {dx_first['f_star'][0]:.1f} and the "
        f"first result is depth-bound at {first_unborrowed_yr:.2f} yr (NEEDS_AUTHOR).",
    )
    return Result(era="2033", lq=(lq, lq), hard_ops=(deepest_t, deepest_t), breakdown=breakdown,
                  intermediates=inter, shots=(n_circuits, n_circuits),
                  wall_time_s=(dn_arm_s, dn_arm_s), epsilon_l=(eps_l_deepest, eps_l_base), notes=notes)


def model(a: Assumptions, era: str) -> Result:
    """Return the Result for '2028' | '2033'. The chapter has no co-design instance."""
    if era == "2028":
        return _model_2028(a)
    if era == "2033":
        return _model_2033(a)
    if era in ERAS:
        raise ValueError("Ch. 8 quotes no codesign instance; the interim scalar-fermion variant "
                         "is an INSTANCE_ROWS row under 2028")
    raise ValueError(f"unknown era {era!r}")


PUBLISHED = {
    "2028": Published(lq=(174, 174), hard_ops=(2.8e5, 2.8e5),
                      src=f"{TEX}:2028 box ('174', '2.8e5 (9.1e5 before the assumed 3.3x step reduction)'): c_T 1239.8 "
                          "x 36 x 20 = 9.095e5 raw / 3.3 = 2.756e5 (c_T 1263.1 at the 20-step budget; ruling 2026-10-11; was 1.623e5 at 12)",
                      rel_tol=0.05),
    "2033": Published(lq=(1000, 1000), hard_ops=(9.9e8, 9.9e8),
                      src=f"{TEX}:2033 box ('1000', 'Deepest circuit >~9.9e8 T at depth 33 (MLAE), dressing ramp unpriced'). "
                          "Exact 33 x 3.0056e7 + 16 x 7.42e4 = 9.930e8 (G2 2026-10-04; B2 free preparation and reflections "
                          "priced 2026-10-05, was 9.763e8); base circuit 3.0e7",
                      rel_tol=0.05),
}


def INSTANCE_ROWS(a, era, r):
    """Rows for resources.json: [(label, (lq_lo, lq_hi), (t_lo, t_hi), {extra}), ...]."""
    base = dict(t_bound=None, conditional=False, codesign=False, reconciled=False)
    if era == "2028":
        i = r.intermediates
        return [
            (r"scalar bubble $6^2$", r.lq, r.hard_ops,
             dict(base, reconciled=True,
                  note=f"{TEX}:94 prints 2.8e5 = 9.095e5/3.3 at 20 steps (c_T 1263.1 derived here, R12; "
                       "slope times 10 and 20 steps, ruling 2026-10-11)")),
            (r"family $8^2$, $N_\phi{=}4$", (i["family_8x8_nphi4_lq"],) * 2,
             (i["family_8x8_nphi4_t_banked"],) * 2, dict(base, note="r16 record (no longer printed): 158 LQ, 5.9e4 T (1.944e5/3.3) at 12 steps")),
            (r"family $8^2$, $N_\phi{=}16$", (i["family_8x8_nphi16_lq"],) * 2,
             (i["family_8x8_nphi16_t_banked"],) * 2, dict(base, note="r16 record (no longer printed): 286 LQ, 2.9e5 T (9.723e5/3.3) at 12 steps")),
            # the interim 4^2 scalar+Wilson variant was removed from the chapter text (2026-10-08);
            # its intermediates are kept below as a record and are no longer published
        ]
    if era == "2033":
        return [(r"wall scattering $15{\times}8$", r.lq, r.hard_ops,
                 dict(base, reconciled=True, conditional=True,
                      note=f"{TEX}:2033 box: deepest MLAE circuit 33 x 3.0056e7 + 16 x 7.42e4 = 9.93e8 T (G2; B2 free "
                           "preparation and reflections priced), a lower bound: dressing ramp unpriced (B2/B3)"))]
    return []


if __name__ == "__main__":
    a = Assumptions()
    for era in PUBLISHED:
        r = model(a, era)
        print(era, "LQ", r.lq, "T", r.hard_ops, "breakdown", r.breakdown_total())
        for k, v in r.intermediates.items():
            print(f"   {k:40s} {v}")
