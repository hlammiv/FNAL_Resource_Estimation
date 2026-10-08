"""Ch. 7 — QFT in curved spacetime. Reproduces the LQ and hard-op numbers of
applications/app04_curved_space.tex.

The chapter's cost model (app04:62-67) is
    N_q          = V (ceil(log2 K) + 4 N_f) + N_anc
    N_gate/shot  = (T_phys / dt) * Gamma_step
    N_shot       = O(eps^-2) * N_modes
Since R12 (2026-09-30) Gamma_step is DERIVED HERE gate by gate (app04:84), not
Stated. The digitization is the uniform field-amplitude grid of Macridin et al.
(arXiv:2108.10793) with the term circuits of Li-Macridin-Mrenna-Spentzouris
(arXiv:2210.07985 Sec. III.C, main.tex:338-458): phi_j is linear in the n_q site
Z's (Eq. phiqub, :338-341), so
    phi^2      : C(n_q,2) ZZ rotations per site          (Eq. trphi2, :382-395)
    phi^4      : C(n_q,4) Z^4 + C(n_q,2) ZZ per site     (Eq. trphi4, :418-440); the ZZ merge with phi^2
    pi^2       : C(n_q,2) ZZ per site between F and F^-1 (Eq. trpi2, :397-407)
    phi_j phi_l: n_q^2 ZZ per link                       (Eq. trphijphik, :409-416)
    F          : Rz-layer . QFT . Rz-layer               (Eq. fftqft, :351-353, App. :1082-1125)
The paper gives CNOT counts only (Table I, :442-455), so every T figure is ours.
Derived here, with no paper source: the QFT T price (controlled-S = 3 T; each
smaller controlled phase = 3 synthesized R_z), the Wilson hop (8 nonzero entries of
gamma^0(r - gamma_i), Dirac rep, r = 1; 2 JW Pauli strings each -> 16 per link),
the Wilson mass (4 single-Z per site) and the Yukawa g phi psibar psi (app04:84,
author-confirmed 2026-09-30; 4 n_q ZZ per site; tr gamma^0 = 0 leaves no phi-only string). The step is the
symmetric second-order (Stormer-Verlet) splitting the chapter's stability argument
uses (app04:81): with a hopping block, the diagonal block is outermost and merges
across steps and the hopping block is applied twice; with none (2028) the pi^2 block
is outermost. All rotations of one shot are synthesized under R-TOL
(eps_rot = sqrt(1e-2/N_rot), 1.15 log2(1/eps) + 9.2 T).

INPUTS AND SOURCES
  L, K, N_f (2028)     4, 8, 0        app04:112-113   4^3 free real scalar, K=8 -> 3 qubits/site (option (c), 2026-10-05)
  L, K, N_f (2033)     5, 16, 1       app04:137-138   5^3 phi^4 (K=16 -> 4 q/site) + 1 Wilson flavor
  Wilson q/site/flavor 4              app04:69        "standard Wilson-fermion encoding"
  N_anc                30-50          app04:70,120    ramp (2033 only) / out-mode / Trotter control
  N_Hadamard           0              app04:71        R6 (r25): direct computational-basis readout, no Hadamard test (was 1-2)
  dt (2028)            0.5/H_inf      app04:73,123    omega_max ~ 2 H_inf
  N_steps (2028)       1              app04:77        one macro-step, as many as the cap holds at K = 8 (option (c); was 4 at K = 4)
  T prep (2028)        3e4            app04:78        Gaussian-network in-vacuum, "~150 T per system qubit" (2e4 x 192/128); Stated
  hard-op cap (2028)   1e5            app04:122       "the 1e5 cap" (RFI 2028 floor)
  T_phys (2033)        10/H_inf       app04:81,139
  dt (2033)            5e-4/H_inf     app04 (U2)      m_phi dt = 0.05: step-convergence proxy (verlet_proxy); was 0.01/H_inf
  commutator lever     1 (none)       app04           no lever
  stability window     omega dt < 2   app04           Stormer-Verlet map on the pi^2 + m^2 phi^2 part; a validity check only
  hard-op ceiling 2033 1e9            DOE RFI 2026    reference budget; the 2033 shot is 14.9x it (U2 reprice)
  ramp window (2033)   5/H_inf        app04           adiabatic switch-on of the same H (R17 ruling; was 10/m_phi)
  prep steps (2033)    500            derived         window / dt_ramp, dt_ramp = 0.01/H_inf (m_phi dt = 1, phibar = 0)
  condensate boost     V n_q = 500    app04 (U1)      exp(i pibar sum_x phi_x) after the ramp: single-Z rotations, once
  per-step T           derived        app04:84        term counts above, priced at R-TOL
  readout settings     3 (2028), 2 (2033)  app04:79   R6 (r25): grouped equal-time two-point readout
  shots/setting (2028) 150            app04:79        5% on the worst IR shell (3 complex modes): >= 400/3 = 134
  shots/profile (2033) N_set ((n+1/2)/n)^2 / (m' eps^2), n = 1, m' = 3, x1-2 margin   app04:79 (shot audit r25)
  precision (2033)     30% first result, 10% campaign   R1 (r25)
  profiles (2033)      1 first result, 10-20 campaign   app04:88 (shot audit: 10-20 serve the EOS; was 50)
  t_gate_s             1 us           app04:88        1 us per T
  shot_overhead_s      0.1 ms         report rule (main-overview.tex:175)
  T-depth per shot     factory.json   2028 346-1308; 2033 27 / 92 rotation layers per step (r25) x T/rot x steps
  fault budget         0.1/shot       R3: 0.1 expected logical faults per shot, report-wide; T gates only (G4)
  eps_l (2033)         <~7e-12        app04           0.1 / 1.486e10 = 6.73e-12
  RFI eps_l            1e-8           app04           -> ~149 expected faults per shot; no split-circuit route (U4)
  campaign horizon     5 yr           app04           the 10-20-profile campaign takes 0.71-2.8 yr on one machine
  envelope steps       O(10)          app04:90        needs <~ 7e3 T/step at K = 8, a cut of about 9x (was 8e3, 2.3x at K = 4)
  de Sitter variant    20 steps       app04:103       T_phys = 10/H_inf at dt = 0.5/H_inf: ~1.2e7 T
  H_inf                1e-5 M_Pl      app04:139,159   enters no arithmetic
  utility              3% x $109M/yr  app04:164       DOE_HEP_FY25_CJ Cosmic Frontier line

WHAT IS NOT DERIVED HERE
  The 2028 Gaussian-network prep (3e4 T at K = 8, Stated; app04:78 says so). The 2033 ramp's
  starting free vacuum (Gaussian network for 500 scalar qubits plus the free-fermion
  vacuum) is not in the chapter's headline and is not priced (~7.5e4 T at 150 T/qubit,
  1e-4 of the shot). N_anc (Stated). The Yukawa form
  is Stated (app04:84, author-confirmed). The 2033 fermion out-basis change (FFFT + Bogoliubov) is bounded,
  not priced: <= N(N-1)/2 Givens rotations on N = 500 modes, ~7.4e6 T (5e-4 of the shot). The prep-only vacuum
  reference run (2.4% of a shot) and the Yukawa readout setting (no T) are not priced. The U2 step comes from
  verlet_proxy and is checked by exact_reduced_step (exact, 3-site ring, actual couplings): <= 0.3% there, <= 2%
  scaled to 5^3. The same check shows that a field grid fixed at a(t_in) fails (fermion n_k off by up to 30% at
  m_phi t = 50, a^3 ~ 3; fixed_grid_err_t50), so the register holds chi = a^{3/2} phi. Its squeeze gamma(chi p + p chi)/2
  costs nothing: p^2/2 + gamma(chi p + p chi)/2 = U (p^2/2) U^dag - gamma^2 chi^2/2 with U = exp(-i gamma chi^2/2),
  and the chi^2 phases are ZZ angle changes on the phi^2 rotations of the diagonal block (enc = "shear" in
  exact_reduced_step; matches the explicit squeeze to 1e-6 in n_k). The shot-count
  formula's band (n ~ 1), worst shell (m' = 3) and x2 margin come from the r25 shot audit (estimates).

WALL TIME AND FACTORIES
  t_gate_s = 1e-6 s per T, shot_overhead_s = 1e-4 s per shot, ONE machine, serial.
  T-depth per shot from factory.json (Ch. 7 runs): 2028 low 346 / high 1308 (K = 4 four-step circuit, parallel prep;
  kept for the K = 8 one-step circuit, which has 36 rotation layers on its busiest qubit against 57; factory.json
  not rerun; a fully sequential prep would be ~3e4 deep, F* ~ 3, immaterial to the T-bound 43 s); 2033 low 1.11e6 / high 3.79e6 (naive JW tail
  blocking 8.47e7, F* = 12). F* = N_T / D_T = 72-273 (2028), 266-908 (2033): >= 10, no fix needed. The
  10-factory baseline holds 10 RUS rotations in flight, 10 ancillas inside N_anc.
  First result (2033): 1 a(t) profile, 2 settings, 30% per shell in the n ~ 1 band, 17-34 shots
  -> wall_first_result_s 4.8-9.5 h. Campaign: 10-20 profiles at 10% per shell, 150-300 shots each
  -> wall_campaign_s 17.5-70 d. 2028: one tier (3 settings x 150 shots, 42.5 s); wall_first_result_s None.

HISTORY
  2026-09-28 (apply_log/ch07.md): R3 fault budget, R4 carry-exact, mechanical batch.
  2026-09-30 (apply_log/r11_ch07.md): prep in the 2033 headline (A); no x10 lever (a).
  2026-09-30 (apply_log/r12_ch07.md): per-step T derived here. 2028: 1e5 -> 1.886e4
  T/step, shot 1.2e5 -> 3.886e4 (fits; 4 macro-steps would). 2033: 1e7 -> 6.638e5
  T/step, shot 1.01e10 -> 6.707e8 (fits the 1e9 ceiling with no lever).
  2026-09-30 (apply_log/r13_ch07.md): Yukawa g phi psibar psi author-confirmed;
  yukawa_on Assumed -> Stated (app04:84). No number moves.
  2026-10-01 (apply_log/r17_ch07.md): author ruling (d). 2028 runs the four macro-steps the
  cap holds: shot 3.886e4 -> 9.437e4 T. 2033 ramp window 10/m_phi -> 5/H_inf (500 steps):
  shot 6.707e8 -> 1.008e9 T, 1% over the 1e9 ceiling (small overshoot accepted per ruling).
  2026-10-02 (apply_log/r23_ch07.md): author ruling E27, one machine everywhere. The input
  n_machines (10) and the intermediates n_machines, scan_machine_weeks (499.9),
  scan_machine_years (9.58), scan_wall_time_weeks (49.99), machines_for_campaign_horizon
  (1.92), machines_for_campaign_horizon_whole (2), machines_for_one_year (9.58) are retired
  (RETIRED_R23 keeps the old values). In their place: the serial scan time on one machine
  (9.58 yr), the whole profiles one machine completes in the 5-year horizon (26), and the
  cost reduction the full scan needs to fit (9.58 / 5 = 1.92x). Wall times now carry the
  report's 0.1 ms per-shot overhead (2028: 377.47 -> 377.87 s; 2033 profile: +0.6 s). No
  per-shot T or LQ moves.
  2026-10-04 (U2 exact check): exact_reduced_step runs the chapter's splitting exactly on a 3-site ring at the
  actual couplings, K = 16 and one Wilson fermion: step error on the fermion n_k <= 0.3% at m_phi dt = 0.05,
  <= 1.23x the proxy on the same lattice, so <= 2% scaled to 5^3; the step stands and nothing is repriced. The
  check also finds that the grid fixed at a(t_in) cannot hold the field (fermion n_k off by up to 30% at
  m_phi t = 50). The comoving register's squeeze is a conjugation of the pi^2 block by chi^2 phases that merge
  into the diagonal block (U2 verifier, enc = "shear"), so it adds no rotation and the count is not conditional.
  No number moves.
  2026-10-04 (editorial_review/responses/ch07.md): referee U1-U4, G4. U1: the condensate is a boost
  exp(i pibar sum phi) of the ramped vacuum (+500 rotations). U2: the 2033 step is set by accuracy, not
  stability: m_phi dt 1 -> 0.05 (verlet_proxy), 1.5e3 -> 2.05e4 steps; shot 1.008e9 -> 1.486e10 T (14.9x the
  1e9 budget), eps_l 1e-10 -> 7e-12, shot 17 min -> 4.1 h, first result 2.9-5.8 d, campaign 0.71-2.8 yr. U4:
  split-circuit intermediates removed. G4: idle qubit-cycles and the one-site-fault bias estimate added. 2028
  numbers unchanged.
  2026-10-02 (apply_log/r25_ch07.md): rulings R1, R6, R9 (H. Lamm). Grouped equal-time readout by direct
  computational-basis measurement: 2028 3 settings x 150 = 450 shots (was 400 x 10 = 4e3), 42.5 s (was
  6.3 min); no Hadamard qubits (LQ 159-180 -> 158-178, box 160-180 unchanged). 2033 two tiers: first result
  1 profile at 30%, 17-34 shots, 4.8-9.5 h; campaign 10-20 profiles at 10%, 1.5-6e3 shots, 17.5-70 d (was
  50 x 6e3 = 3e5 shots, 9.58 yr, 26 profiles in 5 yr). gate_time_s -> t_gate_s; R9 depth exports added.
  RETIRED_R25 keeps the old shot fields and scan values. No per-shot T moves.
  2026-10-05 smaller ruling (b) (Henry: "do ... the smaller things"; Claude's recommended default, the NEEDS_AUTHOR
  option (c)): the 2028 benchmark runs one macro-step at K = 8 (the priced four-step K = 4 run missed its 5% target,
  ds2028_grid_check). K 4 -> 8, steps 4 -> 1, prep 2e4 -> 3e4 (Stated 2e4 scaled by the system qubits, 192/128):
  LQ 158-178 -> 222-242 (box 160-180 -> 220-240), T 94,366 -> 94,753 (9.4e4 -> 9.5e4), every shell <= 0.15% at
  m ~ 1-1.5 H_inf. Shots (450) unchanged; wall 42.51 -> 42.68 s (43 s). RETIRED_2028_K4 keeps the old values.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, fields

from estimates.common import (Assumed, Cited, Stated, Uncited, Tagged, Primitive,
                              CircuitStatus, Result, Published, ERAS,
                              t_per_rotation, eps_rot_for, hwp_group, EPS_SYN, EPS_SYN_SRC,
                              depth_exports)

TEX = "app04"
SECONDS_PER_WEEK = 7 * 86400.0
SECONDS_PER_YEAR = 365.25 * 86400.0
SRC_LMMS = "arXiv:2210.07985 Sec. III.C (main.tex:338-458)"
DERIVED = "derived here (app04:84)"

# E27 (2026-10-02): the report assumes one quantum machine. These multi-machine quantities were
# printed through r22 and are no longer inputs or intermediates; the values are kept as a record only.
RETIRED_R23 = {
    "n_machines": 10,                               # "~10 would finish it in a year"
    "scan_machine_weeks": 499.9,                    # "~500 machine-weeks"
    "scan_machine_years": 9.58,                     # "9.6 machine-years"
    "scan_wall_time_weeks": 49.99,                  # scan / 10 machines
    "machines_for_campaign_horizon": 1.92,          # scan / 5 yr
    "machines_for_campaign_horizon_whole": 2,       # "Two shot-parallel machines"
    "machines_for_one_year": 9.58,                  # "~10"
}

# R1/R6 (2026-10-02, r25): per-mode shot pricing and the 50-profile scan are replaced by grouped two-point
# readout and two tiers. Old inputs and printed values, kept as a record only.
RETIRED_R25 = {
    "shots_per_mode_2028": 400, "n_modes_2028": 10, "shots_2028": 4000, "wall_2028_s": 377.866,
    "shots_per_mode_2033": 100, "n_modes_2033": 60, "shots_per_profile_2033": 6000,
    "wall_per_profile_weeks": 10.0, "n_profiles": 50, "scan_serial_years": 9.58,
    "profiles_in_campaign_horizon": 26, "scan_reduction_to_fit_horizon": 1.92,
    "subroutine_calls_in_campaign_horizon": 1.56e5, "n_hadamard": (1, 2),
}


# smaller ruling (b), 2026-10-05: the four-step K = 4 benchmark (R17 (d)) is replaced by option (c). Record only.
RETIRED_2028_K4 = {
    "K_2028": 4, "n_steps_2028": 4, "t_prep_2028": 2e4, "lq_2028": (158, 178), "box_lq_2028": (160, 180),
    "t_evolution_2028": 74366.42, "t_per_shot_2028": 94366.42, "t_per_step_2028": 18591.6, "wall_2028_s": 42.510,
    "grid_ir_err_end": (0.111, 0.239),
}


@dataclass(frozen=True)
class Assumptions:
    # ---- lattice and field content -----------------------------------------
    L_2028:   Tagged = Stated(4, "app04:112", "Volume: V = 4^3 (64 sites)")
    K_2028:   Tagged = Stated(8, "app04:113", "free real phi, (K = 8, 3 qubits/site); smaller ruling (b) 2026-10-05, "
                                         "NEEDS_AUTHOR option (c) (was K = 4)")
    N_f_2028: Tagged = Stated(0, "app04:113", "free real phi: no fermions")
    phi4_2028: Tagged = Stated(0, "app04:113", "free scalar: no phi^4 vertex")
    L_2033:   Tagged = Stated(5, "app04:137", "Volume: V = 5^3 (125 sites)")
    K_2033:   Tagged = Stated(16, "app04:138", "phi^4 (K = 16, 4 qubits/site)")
    N_f_2033: Tagged = Stated(1, "app04:138", "+1 Wilson fermion (4 qubits/site)")
    phi4_2033: Tagged = Stated(1, "app04:138", "phi^4 self-coupling")
    wilson_qubits_per_site: Tagged = Stated(
        4, "app04:69", "standard Wilson-fermion encoding (4 qubits/site per flavor); no cite")
    # ---- per-term circuit inputs (derived here, app04:84) ------------------
    links_per_site: Tagged = Assumed(3, "periodic cubic lattice: 3V links", src=DERIVED)
    wilson_hop_bilinears_per_link: Tagged = Assumed(
        8, "nonzero entries of gamma^0 (r - gamma_i), Dirac rep, r = 1, each i "
           "(checked by wilson_hop_nonzeros())", src=DERIVED)
    pauli_strings_per_bilinear: Tagged = Assumed(
        2, "JW: c_a^dag c_b + h.c. with a real or imaginary coefficient = 2 Pauli strings", src=DERIVED)
    wilson_mass_strings_per_flavor: Tagged = Assumed(
        4, "(m + 3r) psibar psi = sum_a (+-) n_a: 4 single-Z strings per site", src=DERIVED)
    yukawa_on: Tagged = Stated(
        1, "app04:84", "author-confirmed 2026-09-30: g phi psibar psi, the scalar->fermion transfer "
           "of app04:100; 4 n_q ZZ per site per flavor (tr gamma^0 = 0: no phi-only string)")
    t_per_controlled_s: Tagged = Assumed(
        3, "controlled-S = three pi/4 phases (T, T, T^dag) + 2 CNOT", src=DERIVED)
    rot_per_controlled_phase: Tagged = Assumed(
        3, "controlled-R_k, k >= 3 = three R_z(2 pi / 2^(k+1)) + 2 CNOT, each synthesized", src=DERIVED)
    trotter_order: Tagged = Stated(
        2, "app04:81,84", "symmetric second-order (Stormer-Verlet) splitting")
    eps_syn: Tagged = Cited(EPS_SYN, EPS_SYN_SRC, "R-TOL: 1e-2 synthesis error per shot")
    # ---- ancilla registers -------------------------------------------------
    n_anc: Tagged = Stated((30, 50), "app04:70",
                           "N_anc = 30-50 covers the ramp (2033 only), Trotter-control, out-mode registers")
    n_hadamard: Tagged = Stated(0, "app04:70",
                                "R6 (H. Lamm 2026-10-02, r25): every observable is an equal-time two-point function "
                                "read by direct measurement in the computational basis; no Hadamard test (was 1-2)")
    # ---- 2028 time stepping -------------------------------------------------
    dt_Hinf_2028: Tagged = Stated(0.5, "app04:73", "Delta t ~ 0.5/H_inf")
    omega_max_Hinf_2028: Tagged = Stated(2, "app04:73,123", "omega_max ~ 2 H_inf")
    n_steps_2028: Tagged = Stated(
        1, "app04:77", "one Trotter macro-step at K = 8, as many as the cap holds with the prep (R17 (d) rule; "
                       "smaller ruling (b) 2026-10-05, option (c); was four at K = 4)")
    t_prep_2028: Tagged = Stated(3e4, "app04:78",
                                 "Gaussian-network preparation ... ~150 T per system qubit: the Stated 2e4 at 128 "
                                 "system qubits scaled to 192 (K = 8); 'stated, not derived'")
    hard_op_cap_2028: Tagged = Stated(1e5, "app04:122", "'the 1e5 cap' (RFI 2028 floor)")
    # ---- 2033 time stepping -------------------------------------------------
    T_phys_Hinf_2033: Tagged = Stated(10, "app04:81", "T_phys H_inf = 10 (ten Hubble times)")
    dt_Hinf_2033: Tagged = Stated(
        5e-4, "app04 (U2)", "Delta t = 5e-4/H_inf, m_phi Delta t = 0.05: <= 1.6% step error on n_k in the n >~ 0.1 band "
                            "in all three proxy coupling sets at Phi0 = 0.31 m_phi (verlet_proxy; 0.1 gives up to 7%). "
                            "Was 0.01/H_inf")
    accuracy_omega_dt: Tagged = Assumed(
        0.05, "U2 accuracy criterion: m_phi Delta t at which the proxy step error on n_k is <= 5%",
        src="verlet_proxy (editorial_review/responses/ch07.md)")
    lattice_spacing_mphi_2033: Tagged = Stated(
        3.46, "app04 (U3)", "b ~ 3.5/m_phi at 2033 (5^3 box of side 17.3/m_phi)")
    condensate_amplitude_2033: Tagged = Stated(
        0.31, "app04 (U1, verifier 2026-10-04)", "Phi0 ~ 0.3 m_phi: m b^3 Phi0^2 / 2 ~ 2 quanta per site, the K = 16 "
                                                 "limit (coherent_fidelity); the proxy sets use this amplitude")
    proxy_couplings_2033: Tagged = Assumed(
        3, "(lambda, g) = (1, 1), (3, 0.5), (1, 2) at Phi0 = 0.31, m_psi = 0.1 m_phi (verlet_proxy)",
        src="editorial_review/responses/ch07.md")
    omega_lattice_max_over_mphi: Tagged = Stated(
        1.41, "app04 (U3)", "highest free lattice frequency sqrt(1 + 12/(b m)^2) at b = 3.46/m_phi")
    condensate_boost: Tagged = Stated(
        1, "app04 (U1)", "the condensate is the boost exp(i pibar sum_x phi_x) of the ramped vacuum: V n_q single-Z "
                         "rotations in the field basis, once")
    omega_max_Hinf_2033: Tagged = Stated(1e2, "app04:73", "omega_max ~ 1e2 H_inf (m_phi ~ 1e2 H_inf)")
    commutator_lever: Tagged = Stated(
        1, "app04", "no commutator lever: the step is set by the U2 accuracy criterion")
    stability_omega_dt: Tagged = Stated(
        2, "app04:81", "Stormer-Verlet map on the quadratic part is stable only for omega dt < 2")
    hard_op_ceiling_2033: Tagged = Cited(1e9, "DOE RFI 2026", "2033 reference budget; the U2 shot is 14.9x it")
    ramp_window_Hinf_2033: Tagged = Stated(
        5, "app04:81", "adiabatic switch-on window 5/H_inf (author ruling R17 (d), 2026-10-01; was 10/m_phi)")
    dt_ramp_Hinf_2033: Tagged = Stated(0.01, "app04 (U2)", "ramp at Delta t = 0.01/H_inf (m_phi dt = 1; phibar = 0 in the ramp, so omega dt <~ 1.4 < 2); "
                                       "a sudden switch to the evolution step leaves 5e-3 to 3e-2 quanta per mode")
    dt_sweep_steps_2033: Tagged = Stated(
        50, "app04 (U2 verifier)", "linear Delta t sweep 0.01 -> 5e-4 /H_inf over 50 extra steps at the end of the "
                                   "ramp: leaves < 2e-5 quanta per mode (switch_occupation); not priced")
    exact_step_err_reduced_2033: Tagged = Assumed(
        3.1e-3, "U2 exact check (2026-10-04): largest relative step error of the fermion n_k at m_phi dt = 0.05 in "
                "exact_reduced_step (3-site ring at b = 3.46, K = 16, one Wilson fermion, the chapter's a(t), Phi0, m_psi; "
                "comoving register), over the three coupling sets and 8 times t = 12.5-1000; Richardson from dt, dt/2, "
                "dt/4, ratio 4.0. Per set 1.4e-3, 3.1e-3, 1.1e-3", src="exact_reduced_step")
    proxy_step_err_2033: Tagged = Assumed(
        0.016, "largest 5^3 verlet_proxy step error at m_phi dt = 0.05 on fermion shells with n_k > 0.1 (set 3; "
               "against 1e-3 and against Richardson from 0.05/0.025/0.0125 alike)", src="verlet_proxy")
    exact_over_proxy_step_2033: Tagged = Assumed(
        1.23, "largest ratio of exact to proxy step error on the same reduced lattice (set 2; sets 1, 3: 0.97, 0.52). "
              "Scales the 5^3 proxy error to the actual couplings: <= 0.016 x 1.23 = 2%", src="exact_reduced_step")
    fixed_grid_err_t50: Tagged = Assumed(
        0.30, "largest relative error of the fermion n_k (shell k = 2 pi/(3b)) on the K = 16 grid fixed at a = 1, "
              "exact_reduced_step enc = 'fixed' vs 'chi' at m_phi t = 50 (a^3 = 3.06): 0.391 vs 0.300 (set 2), "
              "0.735 vs 0.900 (set 3, -18%), 0.869 vs 0.898 (set 1). Fermion fine at t = 25; scalar already wrong "
              "at t = 25 (0.20-0.30 vs 0.015-0.034). A fixed K = 32 grid fails by t = 100 (scalar 1.0-1.2 vs "
              "0.02-0.06)", src="exact_reduced_step")
    fixed_grid_gap_a3_4: Tagged = Assumed(
        0.932, "K = 16 site gap on the grid fixed at a = 1, at a^3 = 4 (m_phi t ~ 67): 7% low; 0.02 at a^3 = 256 "
               "(t = 1000). Not the headline: the exact run on that grid fails earlier (fixed_grid_err_t50). The "
               "comoving register chi = a^{3/2} phi keeps the grid matched; its squeeze gamma (chi p + p chi)/2 is "
               "applied as exp(-i gamma chi^2/2) exp(-i h p^2/2) exp(i gamma chi^2/2) with -gamma^2 chi^2/2 in the "
               "diagonal block, ZZ angle changes only, no added rotation", src="fixed_grid_gap")
    prep_fraction_printed: Tagged = Stated(0.025, "app04", "the ramp is '2.5%' of the evolution")
    # ---- 2028 field-grid check and G6 comparator (chapter-open pass, 2026-10-05) ----------
    grid_2028_fixed_err_step1: Tagged = Assumed(
        (0.40, 2.49), "U2 for 2028 (ds2028_grid_check): largest relative error of the IR-shell two-point functions "
                      "after the first step on a K = 4 grid fixed at t_in, 4-site ring at the 2028 spacing "
                      "(m^2 + 12/b^2 = 4 H^2), m = 0.5, 1.0, 1.4, 1.5 H_inf: 0.40, 1.02, 2.04, 2.49. The comoving "
                      "register of 2033 (chi = a^{3/2} phi, continuous shear) also fails (0.63-0.86 at m = 0.5, 1)",
        src="ds2028_grid_check")
    grid_2028_k4_two_steps_ir: Tagged = Assumed(
        (0.027, 0.038), "IR-shell error after two steps (T_phys = 1/H_inf) on a K = 4 grid refitted every step to "
                        "the reference site covariance (frame = 'adapted': rescale + shear), priced merged circuit "
                        "(merged = True), m = 1.0 and 1.4 H_inf (0.059 at 1.5; 0.44 at 0.5; k = 0 shell 0.07-0.16 and "
                        "top shell 0.03-0.05 at m = 1-1.4). Split circuit: 0.030, 0.035", src="ds2028_grid_check")
    grid_2028_k4_adapted_ir_end: Tagged = Assumed(
        (0.111, 0.239), "IR-shell error after the four priced steps, K = 4, adapted frame, merged circuit, m = 1.4 "
                        "and 1.5 H_inf (0.20 at m = 1, 0.66 at 0.5; split circuit 0.105, 0.237)",
        src="ds2028_grid_check")
    grid_2028_k8_one_step: Tagged = Assumed(
        (1.04e-3, 1.42e-3), "option (c), the priced 2028 circuit: largest error on any shell after one step, K = 8, "
                            "adapted frame, merged circuit, m = 1.0, 1.4, 1.5 H_inf: 1.42e-3, 1.05e-3, 1.04e-3 "
                            "(0.116 at m = 0.5)", src="ds2028_grid_check")
    map_vs_desitter_2028_one_step: Tagged = Assumed(
        (0.133, 0.284), "largest relative deviation of the continuum one-step map from exact de Sitter at "
                        "T_phys = 0.5/H_inf, m = 1.0 and 1.5 H_inf (0.246 at 1.4, 0.078 at 0.5)",
        src="ds2028_map_vs_exact")
    grid_2028_k8_adapted_end: Tagged = Assumed(
        1.4e-3, "largest error on any shell after four steps, K = 8, adapted frame, m = 1.4-1.5 H_inf (1.1e-3, 1.4e-3, "
                "merged and split agree to 1e-4; at m = 1: 0.26 merged, 0.30 split; K = 8 one step at m = H: 1.4e-3 "
                "on the k = 0 shell; two steps at m = H: 0.029). Merged vs split differ by <= 0.006 at m >= 1.4 H_inf, "
                "up to 0.04 at m = H", src="ds2028_grid_check")
    map_vs_desitter_2028: Tagged = Assumed(
        (0.19, 0.87), "largest relative deviation of the continuum four-step map (dt = 0.5/H_inf, midpoint a) from "
                      "exact de Sitter evolution on the ring's lattice modes at T_phys = 2/H_inf, m = 0.5-1.5 H_inf; "
                      "second-order convergence checked (0.188, 0.052, 0.013, 0.003 at m = 0.5)",
        src="ds2028_map_vs_exact")
    twopi_lo_fermion_err_late: Tagged = Assumed(
        0.053, "G6: largest relative error of the Gaussian truncation (Hartree + mean-field Yukawa with "
               "fermion back-reaction; twopi_lo) on the fermion n_k of the U2 reduced instance, t = 200-1000/m_phi, "
               "against the exact K = 16 runs at m_phi dt = 0.0125; per set 5.3%, 2.9%, 3.5% (t = 1000: 5.2, 2.9, "
               "1.2%)", src="twopi_lo")
    twopi_lo_fermion_err_early: Tagged = Assumed(
        0.18, "same, all read times from t = 12.5/m_phi (set 3 at t = 12.5: 0.514 vs 0.626)", src="twopi_lo")
    twopi_lo_scalar_err_late: Tagged = Assumed(
        0.66, "same, scalar shell n (k = 2 pi/(3b)), t >= 200: up to 66% low (two entries high: +6.5% set 1 and "
              "+2.5% set 3, both at t = 200)", src="twopi_lo")
    twopi_lo_scalar_err_early: Tagged = Assumed(
        0.96, "same, scalar shell, t = 12.5-100: up to 96% (set 3 at t = 25: 0.0014 vs 0.034; set 2 at t = 50: "
              "0.0054 vs 0.082)", src="twopi_lo")
    # ---- shots, timing, error rate ----------------------------------------
    # R6 (r25): grouped equal-time two-point readout, direct computational-basis measurement
    n_settings_2028: Tagged = Stated(
        3, "app04:79", "field basis, conjugate basis, sheared field basis: the full (phi, pi) covariance "
                       "the de Sitter mode functions fix; the cross term is O(1) at omega ~ H")
    precision_2028: Tagged = Stated(0.05, "app04:79", "~5% on each IR shell")
    worst_shell_modes_2028: Tagged = Assumed(
        3, "lowest nonzero shell of 4^3 holds 6 momenta = 3 complex modes (enumerated, shot audit r25)",
        src="shot_audit.json ch 7")
    shots_per_setting_2028: Tagged = Stated(
        150, "app04:79", "150 per setting >= 400/3 = 134 (5% on the 3-mode shell), shot audit r25")
    n_settings_2033: Tagged = Stated(
        2, "app04:79,88", "field and conjugate basis; the energy-form n_k needs no phi-pi cross term")
    fermion_shell_modes_2033: Tagged = Assumed(
        12, "smallest nonzero 5^3 shell: 6 momenta x 2 helicities of the fermion, each read in its out-basis",
        src="DERIVED")
    scalar_n_max_proxy: Tagged = Assumed(
        0.057, "largest scalar shell occupation in the three proxy sets at Phi0 = 0.31 (verlet_proxy)",
        src="verlet_proxy (editorial_review/responses/ch07.md)")
    occupation_band_2033: Tagged = Stated(
        1, "app04:79,88", "n_k ~ 1: the band beyond classical-statistical simulation (shot audit r25)")
    worst_shell_modes_2033: Tagged = Assumed(
        3, "smallest nonzero shell of 5^3 holds 6 momenta = 3 complex modes (enumerated, shot audit r25)",
        src="shot_audit.json ch 7")
    shot_margin_2033: Tagged = Assumed(
        (1, 2), "1 = strict worst shell at n = 1; 2 covers n ~ 0.5 (x1.8) or loss of shell averaging in a "
                "resonant non-Gaussian state (shot audit r25)", src="shot_audit.json ch 7")
    precision_first_2033: Tagged = Stated(0.3, "app04:79,88", "R1 (H. Lamm 2026-10-02): first result at 30% statistical error")
    precision_campaign_2033: Tagged = Stated(0.1, "app04:79,88", "10% per shell, the chapter's physics target")
    profiles_first_2033: Tagged = Stated(1, "app04:79,88", "first claim: one a(t) profile at one realistic coupling")
    profiles_campaign_2033: Tagged = Stated(
        (10, 20), "app04:79,88", "10-20 profiles serve a smooth equation-of-state fit (shot audit r25; was 50)")
    t_gate_s: Tagged = Stated(1e-6, "app04:88", "at 1 us per T (R9 name; was gate_time_s)")
    shot_overhead_s: Tagged = Assumed(
        1e-4, "per-shot overhead ~0.1 ms (register init, final readout, decode), report rule", src="main-overview.tex:175")
    # R9 (r25): T-depth per shot from factory.json, low / high
    t_depth_2028: Tagged = Assumed(
        (346, 1308), "factory.json ch 7 2028 (K = 4, four steps): low (parity ancilla per string) / high (13 disjoint "
                     "layers, parallel prep). Kept for the K = 8 one-step circuit (option (c)): it has 36 rotation "
                     "layers on its busiest qubit (D 20, two pi^2 halves 2 x 2, two QFT pairs 2 x 6) against 57 "
                     "(4 x 13 + 5); factory.json not rerun", src="factory.json ch 7")
    t_depth_layers_per_step_2033: Tagged = Assumed(
        (27, 92), "factory.json ch 7 2033: low 27 / high 92 rotation layers per step (r25; 1.11e6-3.79e6 at the "
                  "retired 1500 steps); naive JW tail blocking 8.47e7 then", src="factory.json ch 7")
    eps_l_2033: Tagged = Stated(7e-12, "app04",
                                "printed 'eps_l <~ 7e-12' (0.1 faults per shot); exact 6.73e-12")
    fault_budget_per_shot: Tagged = Assumed(
        0.1, "R3 (2026-09-28): 0.1 expected logical faults per shot in the T gates, report-wide; ~90% of shots "
             "fault-free (exp(-0.1) = 0.905), not a heralded acceptance (G4)", src="TRACKED_CHANGES.md R3")
    survival_target: Tagged = Stated(0.9, "app04", "fault-free fraction ~90%")
    eps_l_rfi: Tagged = Cited(1e-8, "DOE RFI 2026", "first-generation RFI figure")
    idle_cycle_s: Tagged = Stated(1e-5, "main-overview.tex:33", "10 us logical cycle of the reference machine (G4 idle count)")
    site_energy_max_quanta: Tagged = Assumed(
        19.4, "highest level of the K = 16 digitized site oscillator above its ground state (digitized_site(16))",
        src="digitized_site, arXiv:2108.10793 grid")
    fermion_fault_dn_max: Tagged = Assumed(
        0.119, "largest shell-averaged added occupation from one X fault on a fermion qubit (JW string) in the free "
               "Wilson vacuum of 5^3, b m_psi = 0.346 (fermion_fault_shift, worst j = 248)",
        src="fermion_fault_shift")
    campaign_years: Tagged = Stated(5, "app04:161", "Campaign horizon 5 years (one machine, serial; E27)")
    n_steps_envelope_2028: Tagged = Stated(
        10, "app04", "'an O(10)-step de Sitter run inside 1e5 hard ops'")
    H_inf_over_Mpl: Tagged = Stated(
        1e-5, "app04:139", "'H_inf = 1e-5 M_Pl setting the reheating scale'; enters no arithmetic")
    # ---- utility box -------------------------------------------------------
    cosmic_frontier_musd_per_yr: Tagged = Cited(109, "DOE_HEP_FY25_CJ",
                                                "$109M/yr DOE Cosmic Frontier Experimental Physics (app04:164)")
    utility_share: Tagged = Assumed(0.03, "3% share; 'a stated, rescalable convention'", src="app04:164")

    def __post_init__(self):
        for f in fields(self):
            v = getattr(self, f.name)
            if not isinstance(v, Tagged):
                raise TypeError(f"{f.name} must be Tagged")
            if v.is_range and v.lo > v.hi:
                raise ValueError(f"{f.name}: range lo > hi")
            if v.lo < 0:
                raise ValueError(f"{f.name}: negative")
        for L, K in ((self.L_2028, self.K_2028), (self.L_2033, self.K_2033)):
            if int(L.lo) < 2 or int(K.lo) < 2:
                raise ValueError("lattice side and field-grid size K must be >= 2")
            if int(K.lo) & (int(K.lo) - 1):
                raise ValueError("K must be a power of two: the uniform grid uses log2 K qubits per site")
        if self.trotter_order.lo != 2:
            raise ValueError("only the symmetric second-order splitting is modeled (app04:81,84)")
        # the chapter's stability rule: dt <~ 1/omega_max (app04:73)
        for dt, w in ((self.dt_Hinf_2028, self.omega_max_Hinf_2028),
                      (self.dt_Hinf_2033, self.omega_max_Hinf_2033)):
            if dt.lo * w.lo > 1.0 + 1e-9:
                raise ValueError("Delta t exceeds 1/omega_max")
        if self.commutator_lever.lo < 1:
            raise ValueError("commutator lever must be >= 1 (it only ever loosens)")
        # a lever lengthens the 2033 step; the product formula must stay inside the
        # Stormer-Verlet window omega dt < 2 on the m_phi mode (app04:81, R11)
        if (self.dt_Hinf_2033.lo * self.commutator_lever.lo * self.omega_max_Hinf_2033.lo
                >= self.stability_omega_dt.lo):
            raise ValueError("levered 2033 step leaves the omega dt < 2 stability window")
        if (self.dt_Hinf_2033.lo * self.commutator_lever.lo * self.omega_max_Hinf_2033.lo
                > self.accuracy_omega_dt.lo + 1e-12):
            raise ValueError("2033 step fails the U2 accuracy criterion (m_phi dt <= 0.05)")
        # the ramp Trotterizes the same H, m_phi^2 phi^2 included, so the same window applies (R12)
        if (self.dt_ramp_Hinf_2033.lo * self.omega_max_Hinf_2033.lo
                >= self.stability_omega_dt.lo):
            raise ValueError("2033 ramp step leaves the omega dt < 2 stability window")
        if self.ramp_window_Hinf_2033.lo <= 0:
            raise ValueError("ramp window must be positive")
        if not (0 < self.survival_target.lo < 1):
            raise ValueError("survival target must be in (0,1)")
        if self.fault_budget_per_shot.lo <= 0:
            raise ValueError("fault budget must be positive")


# --------------------------------------------------------------------------- #
# Register (Eq. Nq_curved)
# --------------------------------------------------------------------------- #

def _era_inputs(a: Assumptions, era: str):
    if era == "2028":
        return a.L_2028, a.K_2028, a.N_f_2028, a.phi4_2028
    if era == "2033":
        return a.L_2033, a.K_2033, a.N_f_2033, a.phi4_2033
    raise ValueError(f"Ch. 7 has no {era!r} box")


def register(a: Assumptions, era: str) -> dict:
    L, K, Nf, _ = _era_inputs(a, era)
    V = int(L.lo) ** 3
    q_scalar = math.ceil(math.log2(int(K.lo)))
    q_fermion = int(a.wilson_qubits_per_site.lo) * int(Nf.lo)
    system = V * (q_scalar + q_fermion)
    lo = system + a.n_anc.lo + a.n_hadamard.lo
    hi = system + a.n_anc.hi + a.n_hadamard.hi
    # app04:53: Hilbert dimension ~ K^V (or 2^{4 N_f V} for fermions); log2 of each
    return dict(V=V, qubits_per_site_scalar=q_scalar, qubits_per_site_fermion=q_fermion,
                system_qubits=system, n_anc=(a.n_anc.lo, a.n_anc.hi),
                n_hadamard=(a.n_hadamard.lo, a.n_hadamard.hi), lq=(lo, hi),
                log2_hilbert_dim_scalar=V * q_scalar, log2_hilbert_dim_fermion=V * q_fermion,
                log2_hilbert_dim=system)


# --------------------------------------------------------------------------- #
# Per-term circuit counts (derived here, app04:84)
# --------------------------------------------------------------------------- #

def wilson_hop_nonzeros(r: float = 1.0) -> list[int]:
    """Nonzero entries of gamma^0 (r - gamma_i), i = 1..3, Dirac representation (pure Python)."""
    s = [((0, 1), (1, 0)), ((0, -1j), (1j, 0)), ((1, 0), (0, -1))]
    g0 = [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, -1, 0], [0, 0, 0, -1]]
    out = []
    for si in s:
        gi = [[0] * 4 for _ in range(4)]
        for x in range(2):
            for y in range(2):
                gi[x][2 + y] = si[x][y]
                gi[2 + x][y] = -si[x][y]
        m = [[r * (x == y) - gi[x][y] for y in range(4)] for x in range(4)]
        prod = [[sum(g0[x][z] * m[z][y] for z in range(4)) for y in range(4)] for x in range(4)]
        out.append(sum(1 for x in range(4) for y in range(4) if abs(prod[x][y]) > 1e-12))
    return out


def centering_layer(n_q: int) -> tuple[int, int]:
    """(synthesized rotations, T gates) per site in one R_z layer of the centered transform
    F = Rz . QFT . Rz, angles 2^(n_q-1-q) delta, delta = pi (N-1)/N (arXiv:2210.07985 Eq. fftqft).
    A multiple of pi/2 is Clifford, an odd multiple of pi/4 is one T, anything else is synthesized."""
    N = 2 ** n_q
    delta = math.pi * (N - 1) / N
    rot = t = 0
    for q in range(n_q):
        m = 2 ** (n_q - 1 - q) * delta / (math.pi / 4)
        if abs(m - round(m)) < 1e-9:
            t += round(m) % 2
        else:
            rot += 1
    return rot, t


def step_terms(a: Assumptions, era: str) -> dict:
    """Rotation and T counts of ONE application of each block on the whole lattice."""
    L, K, Nf, phi4 = _era_inputs(a, era)
    V = int(L.lo) ** 3
    n = int(round(math.log2(int(K.lo))))
    nf = int(Nf.lo)
    links = int(a.links_per_site.lo) * V
    c2 = math.comb(n, 2)
    c4 = math.comb(n, 4) if int(phi4.lo) else 0
    hop = int(a.wilson_hop_bilinears_per_link.lo) * int(a.pauli_strings_per_bilinear.lo) * nf
    ferm_site = nf * (int(a.wilson_mass_strings_per_flavor.lo) + int(a.yukawa_on.lo) * 4 * n)
    cs_per_qft = n - 1                                    # controlled-R_2 = controlled-S
    cr_per_qft = c2 - cs_per_qft                          # controlled-R_k, k >= 3
    layer_rot, layer_t = centering_layer(n)
    return dict(
        V=V, n_q=n, links=links,
        scalar_site_strings=c2 + c4,                      # phi^2 (+gradient diagonal, mass) and phi^4
        gradient_link_strings=n * n,                      # phi_j phi_l
        fermion_site_strings=ferm_site,                   # Wilson mass + Yukawa
        pi2_site_strings=c2,                              # pi^2 in the conjugate basis
        hop_link_strings=hop,                             # Wilson hop
        cs_per_qft=cs_per_qft, cr_per_qft=cr_per_qft,
        diag_block=V * (c2 + c4) + links * n * n + V * ferm_site,       # "A"
        pi2_block=V * c2,                                                # "B" (rotations only)
        hop_block=links * hop,                                           # "C"
        qft_pair_rot=V * 2 * cr_per_qft * int(a.rot_per_controlled_phase.lo),   # QFT + QFT^-1
        qft_pair_cs=V * 2 * cs_per_qft,
        layer_rot=V * layer_rot, layer_t=V * layer_t,
    )


def shot_counts(a: Assumptions, era: str, n_steps: float, extra_rot: int = 0) -> dict:
    """Block applications over a shot of `n_steps` symmetric second-order steps.

    With a hopping block (2033): D/2 [C_sym . B] D/2, D merges across steps -> D x (n+1), B x n,
    C x 2n (every hop string twice; the middle group once in truth, a <1% overcount), QFT pairs x n.
    Without (2028): B/2 D B/2, B outermost -> D x n, B x (n+1), QFT pairs x (n+1).
    The two centering R_z layers of F cancel between consecutive pi^2 blocks (the blocks between
    them are diagonal or act on the fermion qubits), so a shot has two of them."""
    s = step_terms(a, era)
    if s["hop_block"]:
        nD, nB, nC = n_steps + 1, n_steps, 2 * n_steps
    else:
        nD, nB, nC = n_steps, n_steps + 1, 0
    rot = dict(diag=nD * s["diag_block"], pi2=nB * s["pi2_block"], hop=nC * s["hop_block"],
               qft=nB * s["qft_pair_rot"], centering=2 * s["layer_rot"], boost=extra_rot)
    n_rot = sum(rot.values())
    t_exact = dict(qft_cs=nB * s["qft_pair_cs"] * a.t_per_controlled_s.lo, centering=2 * s["layer_t"])
    eps = eps_rot_for(n_rot, a.eps_syn.lo)
    t_rot = t_per_rotation(eps)
    total = n_rot * t_rot + sum(t_exact.values())
    per_step_rot = s["diag_block"] + s["pi2_block"] + 2 * s["hop_block"] + s["qft_pair_rot"]
    per_step_t = per_step_rot * t_rot + s["qft_pair_cs"] * a.t_per_controlled_s.lo
    return dict(terms=s, n_blocks=(nD, nB, nC), rot=rot, n_rot=n_rot, t_exact=t_exact,
                eps_rot=eps, t_per_rot=t_rot, t_total=total,
                per_step_rot=per_step_rot, per_step_t=per_step_t)


def _primitives(a: Assumptions, sc: dict, era: str) -> tuple:
    s, rot, tr, eps = sc["terms"], sc["rot"], sc["t_per_rot"], sc["eps_rot"]
    nD, nB, nC = sc["n_blocks"]
    note = f"eps_rot = {eps:.2e} over N_rot = {sc['n_rot']:.0f} (R-TOL), {tr:.2f} T each"
    diag_note = (f"{nD:g} x [{s['V']} sites x {s['scalar_site_strings']} (phi^2+phi^4) + {s['links']} links x "
                 f"{s['gradient_link_strings']} (phi_j phi_l) + {s['V']} x {s['fermion_site_strings']} "
                 f"(Wilson mass + Yukawa)]; ")
    prims = [
        Primitive("diagonal_z_string_rotations", rot["diag"], tr, CircuitStatus.COMPILED,
                  f"{SRC_LMMS} Eqs. trphi2, trphi4, trphijphik; fermion strings {DERIVED}", diag_note + note),
        Primitive("pi2_conjugate_basis_rotations", rot["pi2"], tr, CircuitStatus.COMPILED,
                  f"{SRC_LMMS} Eq. trpi2", f"{nB:g} x {s['V']} sites x {s['pi2_site_strings']} ZZ; " + note),
        Primitive("qft_small_controlled_phase_rotations", rot["qft"], tr, CircuitStatus.COMPILED,
                  f"QFT per arXiv:2210.07985 Eq. fftqft; T pricing {DERIVED}",
                  f"{nB:g} QFT pairs x {s['V']} sites x 2 x {s['cr_per_qft']} controlled-R_k(k>=3) x "
                  f"{a.rot_per_controlled_phase.lo:g} R_z; " + note),
        Primitive("qft_controlled_s", nB * s["qft_pair_cs"], a.t_per_controlled_s.lo, CircuitStatus.COMPILED,
                  f"QFT per arXiv:2210.07985 Eq. fftqft; T pricing {DERIVED}",
                  f"{nB:g} QFT pairs x {s['V']} sites x 2 x {s['cs_per_qft']} controlled-S at 3 T (Clifford+T)"),
        Primitive("qft_centering_layers", 2 * s["layer_t"], 1.0, CircuitStatus.COMPILED,
                  f"arXiv:2210.07985 Eq. fftqft (App. main.tex:1082-1125); {DERIVED}",
                  "two R_z layers per shot; odd multiples of pi/4 cost one T"),
    ]
    if rot["centering"]:
        prims.append(Primitive("qft_centering_rotations", rot["centering"], tr, CircuitStatus.COMPILED,
                               f"arXiv:2210.07985 Eq. fftqft; {DERIVED}",
                               "two R_z layers per shot, non-Clifford+T angles; " + note))
    if rot.get("boost"):
        prims.append(Primitive("condensate_boost_rotations", rot["boost"], tr, CircuitStatus.COMPILED,
                               "derived here (U1)", f"exp(i pibar sum_x phi_x): {s['V']} sites x {s['n_q']} single-Z, "
                               "once after the ramp; " + note))
    if rot["hop"]:
        prims.append(Primitive("wilson_hop_pauli_strings", rot["hop"], tr, CircuitStatus.COMPILED, DERIVED,
                               f"{nC:g} x {s['links']} links x {s['hop_link_strings']} strings "
                               f"(8 bilinears of gamma^0(1-gamma_i) x 2 JW strings); " + note))
    return tuple(prims)


# --------------------------------------------------------------------------- #
# Model
# --------------------------------------------------------------------------- #

def model(a: Assumptions, era: str) -> Result:
    """Return the Result for '2028' | '2033'. The chapter has no codesign box."""
    if era not in ERAS:
        raise ValueError(era)
    if era == "codesign":
        raise ValueError("Ch. 7 quotes no codesign instance (app04 has two boxes only)")
    reg = register(a, era)
    gt = a.t_gate_s.lo
    if era == "2028":
        return _model_2028(a, reg, gt)
    return _model_2033(a, reg, gt)


def _model_2028(a: Assumptions, reg: dict, gt: float) -> Result:
    n_steps = a.n_steps_2028.lo
    dt = a.dt_Hinf_2028.lo
    T_phys = n_steps * dt                                   # 0.5 Hubble times (1 x 0.5/H_inf; option (c))
    sc = shot_counts(a, "2028", n_steps)
    t_evol = sc["t_total"]                                  # 64,753 (derived here; K = 4 four steps: 74,366)
    t_step = t_evol / n_steps                               # 64,753 per macro-step at the shot's eps_rot
    t_one_step = shot_counts(a, "2028", 1)["t_total"]       # 64,753: one step alone at its own eps_rot
    t_prep = a.t_prep_2028.lo                               # 3e4 (Stated, scaled to 192 system qubits)
    t_total = t_evol + t_prep                               # 94,753 (was 94,366)
    cap = a.hard_op_cap_2028.lo
    # how many macro-steps the cap holds with the prep (app04:77, :89)
    n_fit, shot_at = 0, {}
    while True:
        tn = shot_counts(a, "2028", n_fit + 1)["t_total"] + t_prep
        shot_at[n_fit + 1] = tn
        if tn > cap:
            break
        n_fit += 1
    # R6: three computational-basis settings; every mode is read from the same shots
    n_set = int(a.n_settings_2028.lo)
    samples = 1.0 / a.precision_2028.lo ** 2                # 400 samples for 5%
    shots_min_per_setting = math.ceil(samples / a.worst_shell_modes_2028.lo - 1e-9)   # 134
    if a.shots_per_setting_2028.lo < shots_min_per_setting:
        raise ValueError("2028 shots per setting below the 5% worst-shell minimum")
    shots = n_set * a.shots_per_setting_2028.lo             # 450
    t0 = a.shot_overhead_s.lo                               # 0.1 ms per shot (report rule)
    t_shot_s = t_total * gt + t0                            # 1 us per T + 0.1 ms
    wall_s = shots * t_shot_s                               # serial on one machine
    n_env = a.n_steps_envelope_2028.lo
    if n_steps > n_fit:
        raise ValueError(f"2028: {n_steps:g} macro-steps exceed the {n_fit} the cap holds")
    t_step_envelope = (cap - t_prep) / n_env                # 8e3 T/step
    shot_10 = shot_counts(a, "2028", n_env)["t_total"] + t_prep
    survival_rfi = math.exp(-t_total * a.eps_l_rfi.lo)
    faults_rfi = t_total * a.eps_l_rfi.lo
    eps_l_for_target = -math.log(a.survival_target.lo) / t_total
    eps_l_for_budget = a.fault_budget_per_shot.lo / t_total
    hwp = hwp_group(32, sc["eps_rot"], synthesis="rus")
    breakdown = _primitives(a, sc, "2028") + (
        Primitive("gaussian_network_prep", 1, t_prep, CircuitStatus.SCALING, a.t_prep_2028.src,
                  "normal-mode circuit for the Bunch-Davies vacuum; '~150 T per system qubit' (2e4/128 = 156); "
                  "Stated, not derived (app04:78)"),)
    s = sc["terms"]
    inter = {
        "V": reg["V"],                                              # 64
        "qubits_per_site": reg["qubits_per_site_scalar"],           # 3 (K = 8)
        "system_qubits": reg["system_qubits"],                      # 192
        "log2_hilbert_dim": reg["log2_hilbert_dim"],                # 192
        "n_anc": reg["n_anc"],
        "n_hadamard": reg["n_hadamard"],
        "lq_sum": reg["lq"],                                        # (222, 242); box says 220-240
        "dt_Hinf": dt,
        "omega_max_Hinf": a.omega_max_Hinf_2028.lo,
        "dt_times_omega_max": dt * a.omega_max_Hinf_2028.lo,
        "expansion_time_over_step": 1.0 / dt,
        "time_ordering_over_trotter_error": 1.0 / a.omega_max_Hinf_2028.lo,
        "n_trotter_steps": n_steps,                                 # 1 (option (c); was 4)
        "T_phys_Hinf": T_phys,
        # per-term counts (app04:84)
        "links": s["links"],                                        # 192
        "scalar_site_strings": s["scalar_site_strings"],            # 3 (C(3,2))
        "gradient_link_strings": s["gradient_link_strings"],        # 9
        "pi2_site_strings": s["pi2_site_strings"],                  # 3
        "diag_block_rotations": s["diag_block"],                    # 1920 = 64 x 3 + 192 x 9
        "pi2_block_rotations": s["pi2_block"],                      # 192
        "qft_pair_rotations": s["qft_pair_rot"],                    # 384 = 64 x 2 x 1 x 3
        "n_rot_per_shot": sc["n_rot"],                              # 3200 = 1920 + 2 x 192 + 2 x 384 + 128 centering
        "t_exact_per_shot": sum(sc["t_exact"].values()),            # 1664 = 1536 cS + 128 centering
        "t_controlled_s": sc["t_exact"]["qft_cs"],                  # 1536 (2 QFT pairs)
        "t_centering": sc["t_exact"]["centering"],                  # 128
        "eps_rot": sc["eps_rot"],                                   # 1.77e-3
        "t_per_rotation": sc["t_per_rot"],                          # 19.72
        "t_per_step": t_step,                                       # 64,753 (one step: the evolution)
        "t_single_step_alone": t_one_step,                          # 64,753
        "t_evolution": t_evol,                                      # 64,753 ("6.5e4")
        "t_prep": t_prep,                                           # 3e4
        "t_prep_per_system_qubit": t_prep / reg["system_qubits"],   # 156 ("~150")
        "t_per_shot": t_total,                                      # 94,753 ("9.5e4")
        "hard_op_cap": cap,
        "ratio_to_cap": t_total / cap,                              # 0.948 ("0.95 of the cap")
        "macro_steps_that_fit": n_fit,                              # 1 ("one"; a second step: 1.5e5)
        "shot_t_at_n_steps": shot_at,                               # {1: 94,753, 2: 1.47e5}
        "envelope_steps": n_env,
        "t_per_step_for_envelope": t_step_envelope,                 # 7e3
        "step_reduction_for_envelope": t_step / t_step_envelope,    # 9.25 ("about 9x")
        "shot_t_at_envelope_steps": shot_10,                        # 10 derived steps + prep
        "hwp32_t_over_plain": hwp["t"] / (32 * sc["t_per_rot"]),    # lever not taken (app04:84)
        "hwp32_ancilla": hwp["ancilla"],                            # 31 = 32 - w(32) (E26, r22; was 37)
        "eps_l_rfi": a.eps_l_rfi.lo,
        "expected_faults_at_rfi": faults_rfi,                       # 9.5e-4
        "shot_survival_at_rfi": survival_rfi,                       # 0.99905 ("99.91%")
        "fault_budget_per_shot": a.fault_budget_per_shot.lo,
        "eps_l_for_fault_budget": eps_l_for_budget,
        "eps_l_for_survival_target": eps_l_for_target,
        "n_settings": n_set,                                        # 3 (R6)
        "samples_for_precision": samples,                           # 400
        "worst_shell_modes": a.worst_shell_modes_2028.lo,           # 3
        "shots_min_per_setting": shots_min_per_setting,             # 134
        "shots_per_setting": a.shots_per_setting_2028.lo,           # 150
        "shots": shots,                                             # 450
        "stat_precision_worst_shell": 1 / math.sqrt(a.shots_per_setting_2028.lo * a.worst_shell_modes_2028.lo),
        "readout_extra_t": 0.0,                                     # settings reuse the last pi^2 half-step
        "t_gate_s": gt,
        "shot_overhead_s": t0,
        "t_shot_s": t_shot_s,                                       # 0.0949 (0.1 ms overhead included)
        "wall_time_s": wall_s,                                      # 42.68 s ("43 s"), serial on one machine
        "wall_time_min": wall_s / 60.0,
    }
    dx = depth_exports(t_total, a.t_depth_2028.value, shots, gt, t0)
    inter.update({k: dx[k] for k in ("t_per_shot", "t_depth_per_shot", "f_star", "floor_wall_s",
                                     "factories_for_1yr", "wall_serial_s", "baseline_ok", "fits_1yr")})
    inter["wall_first_result_s"] = None                             # one tier (a benchmark)
    inter["wall_campaign_s"] = wall_s
    if abs(dx["wall_serial_s"][0] - wall_s) > 1e-9 * wall_s:
        raise AssertionError("depth_exports serial wall disagrees with the model wall")
    notes = (
        f"LQ: {reg['system_qubits']} + (30-50) = {reg['lq'][0]}-{reg['lq'][1]}, no Hadamard qubits (R6); box 220-240 "
        "(app04:70, box). Option (c) (smaller ruling (b), 2026-10-05): one step at K = 8.",
        f"T/step derived here (app04:84): B/2 D B/2 on 4^3, n_q = {s['n_q']}: D = {s['diag_block']}, B = {s['pi2_block']}, "
        f"QFT pair {s['qft_pair_rot']} rotations; {n_steps:g} step(s): N_rot = {sc['n_rot']:.0f} "
        f"at eps {sc['eps_rot']:.2e} ({sc['t_per_rot']:.2f} T), + {sc['t_exact']['qft_cs']:.0f} T controlled-S "
        f"+ {sc['t_exact']['centering']:.0f} T centering = {t_evol:.0f} T ('6.5e4').",
        f"T/shot: {t_evol:.0f} + {t_prep:.0f} prep (Stated, scaled) = {t_total:.0f}, {t_total / cap:.2f} of the cap (app04:77).",
        f"Macro-steps that fit with prep: {n_fit} ({shot_at[n_fit]:.0f} T); {n_fit + 1} would be "
        f"{shot_at[n_fit + 1]:.0f} (app04:77,89).",
        f"Readout (R6): {n_set} computational-basis settings x {a.shots_per_setting_2028.lo:g} shots "
        f"(>= 400/3 = {shots_min_per_setting}); no extra T (the settings reuse the last pi^2 half-step).",
        f"Wall, one machine, serial: {shots:g} x ({t_total:.0f} x 1 us + 0.1 ms) = {wall_s:.2f} s ('43 s').",
        f"T-depth {dx['t_depth_per_shot'][0]:.0f}-{dx['t_depth_per_shot'][1]:.0f}, F* = "
        f"{dx['f_star'][0]:.0f}-{dx['f_star'][1]:.0f} (>= 10).",
        f"O(10)-step envelope (app04:90): <= {t_step_envelope:.0f} T/step, a {t_step / t_step_envelope:.2f}x cut; "
        f"10 derived steps + prep would be {shot_10:.3g} T.",
        f"HWP (32-rotation groups) would cost {hwp['t'] / (32 * sc['t_per_rot']):.2f} of plain synthesis at "
        f"{hwp['ancilla']} workspace qubits per group; not taken (app04:84).",
    )
    return Result(era="2028", lq=reg["lq"], hard_ops=(t_total, t_total), breakdown=breakdown,
                  intermediates=inter, shots=(shots, shots), wall_time_s=(wall_s, wall_s),
                  epsilon_l=None, notes=notes)


def _model_2033(a: Assumptions, reg: dict, gt: float) -> Result:
    T_phys, dt = a.T_phys_Hinf_2033.lo, a.dt_Hinf_2033.lo
    n_steps_raw = T_phys / dt                               # 2e4 (m_phi dt = 0.05, U2)
    lever = a.commutator_lever.lo                           # 1 (no lever)
    n_steps_net = n_steps_raw / lever                       # 2e4
    n_prep = round(a.ramp_window_Hinf_2033.lo / a.dt_ramp_Hinf_2033.lo)   # 500 switch-on steps at 0.01/H_inf
    n_total = n_steps_net + n_prep                          # 20,500 steps in one product
    st = step_terms(a, "2033")
    n_boost = int(a.condensate_boost.lo) * st["V"] * st["n_q"]            # 500 single-Z rotations, once (U1)
    sc = shot_counts(a, "2033", n_total, extra_rot=n_boost)
    t_step = sc["per_step_t"]                               # 7.246e5 at the shot's eps_rot
    t_evol = n_steps_net * t_step                           # 1.449e10
    t_raw = n_steps_raw * t_step
    t_prep = n_prep * t_step                                # 3.62e8
    t_boost = n_boost * sc["t_per_rot"]                     # 1.5e4
    t_total = sc["t_total"]                                 # 1.486e10
    t_closing = t_total - t_evol - t_prep - t_boost         # closing D block + centering layers
    ceiling = a.hard_op_ceiling_2033.lo
    t0 = a.shot_overhead_s.lo                               # 0.1 ms per shot (report rule)
    t_shot_s = t_total * gt + t0                            # 1.486e4 s (4.1 h)
    # R1/R6 (r25): two-point readout in N_set settings; per profile N = N_set ((n+1/2)/n)^2 / (m' eps^2)
    n_set = a.n_settings_2033.lo
    nb = a.occupation_band_2033.lo
    m_sh = a.worst_shell_modes_2033.lo
    def shots_per_profile(eps):
        base = n_set * ((nb + 0.5) / nb) ** 2 / (m_sh * eps ** 2)
        return base, (math.ceil(base * a.shot_margin_2033.lo - 1e-9), math.ceil(base * a.shot_margin_2033.hi - 1e-9))
    base_first, spp_first = shots_per_profile(a.precision_first_2033.lo)     # 16.67 -> (17, 34)
    base_camp, spp_camp = shots_per_profile(a.precision_campaign_2033.lo)    # 150 -> (150, 300)
    shots_first = (a.profiles_first_2033.lo * spp_first[0], a.profiles_first_2033.hi * spp_first[1])
    shots_camp = (a.profiles_campaign_2033.lo * spp_camp[0], a.profiles_campaign_2033.hi * spp_camp[1])
    wall_first = (shots_first[0] * t_shot_s, shots_first[1] * t_shot_s)     # 2.9-5.8 d
    wall_camp = (shots_camp[0] * t_shot_s, shots_camp[1] * t_shot_s)        # 0.71-2.8 yr
    # readout circuits (R6): scalar conjugate setting adds one per-site QFT; the fermion out-basis change is bounded
    lr, lt = centering_layer(st["n_q"])
    t_qft_readout = st["V"] * (st["cs_per_qft"] * a.t_per_controlled_s.lo
                               + (st["cr_per_qft"] * a.rot_per_controlled_phase.lo + lr) * sc["t_per_rot"] + lt)
    n_ferm_modes = st["V"] * reg["qubits_per_site_fermion"]                  # 500
    n_givens = n_ferm_modes * (n_ferm_modes - 1) // 2                       # 124,750
    t_ferm_readout = 2 * n_givens * sc["t_per_rot"]                         # ~7.4e6 (bound)
    # U3: the vacuum reference (prep-only run: the ramp without the boost or the evolution) costs this per shot
    t_vacuum_reference = (n_prep + 1) * t_step
    omega_dt = dt * lever * a.omega_max_Hinf_2033.lo        # 0.05
    omega_dt_ramp = a.dt_ramp_Hinf_2033.lo * a.omega_max_Hinf_2033.lo   # 1.0 (stable: phibar = 0 in the ramp)
    # U2: Stormer-Verlet frequency error of a quadratic mode, cos(w~ dt) = 1 - (w dt)^2 / 2
    def freq_err(x):
        return math.acos(1 - x * x / 2) / x - 1
    freq_err_old = freq_err(1.0)                            # 0.0472 at the retired m_phi dt = 1
    phase_err_old = freq_err_old * T_phys * a.omega_max_Hinf_2033.lo   # 47 rad over m_phi T = 1000
    freq_err_new = freq_err(omega_dt)                       # 1.0e-4
    # ramp -> evolution step switch: the coarse-step map vacuum is squeezed against the fine-step one (U2 verifier)
    w_hi = a.omega_lattice_max_over_mphi.lo
    n_switch = (switch_occupation(1.0, omega_dt_ramp, omega_dt),             # 5.2e-3 at the m_phi mode
                switch_occupation(w_hi, omega_dt_ramp, omega_dt))            # 3.0e-2 at the top mode
    n_sweep = int(a.dt_sweep_steps_2033.lo)
    n_switch_swept = max(switch_occupation(w, omega_dt_ramp, omega_dt, n_sweep) for w in (1.0, 1.2, w_hi))
    t_sweep = n_sweep * t_step                              # 3.6e7, not priced
    e_folds = (2.0 / 3.0) * math.log(1.0 + 1.5 * T_phys)   # 1.85
    n_steps_ds = T_phys / a.dt_Hinf_2028.lo                 # 20
    t_ds_variant = shot_counts(a, "2033", n_steps_ds)["t_total"]   # 1.19e7 (its own R-TOL eps)
    budget = a.fault_budget_per_shot.lo
    eps_l_for_budget = budget / t_total                     # 6.73e-12
    faults_at_printed = t_total * a.eps_l_2033.lo           # 0.104
    survival = math.exp(-budget)                            # fault-free fraction, not a heralded acceptance (G4)
    survival_at_printed = math.exp(-faults_at_printed)
    eps_l_for_target = -math.log(a.survival_target.lo) / t_total
    faults_at_rfi = t_total * a.eps_l_rfi.lo                # ~149 expected faults per shot at 1e-8 (U4)
    # G4: excluded fault locations. Idle logical qubit-cycles in one shot, at the 10 us cycle
    idle_qubit_cycles = reg["lq"][1] * t_shot_s / a.idle_cycle_s.lo   # 1.6e12
    idle_over_t = idle_qubit_cycles / t_total               # ~105
    # G4: observable-level bias from unflagged one-site faults (estimate). A fault confined to one site register
    # leaves at most E_max - E_0 quanta there; spread over V modes it adds <= that / V to each n_k, i.e. a relative
    # shift <= 2 (E_max - E_0) / V of the vacuum seed (1/2) that resonance amplifies. Faulted fraction 1 - e^{-0.1}.
    e_site_max = a.site_energy_max_quanta.lo
    dn_per_fault = e_site_max / reg["V"]                    # 0.155 per scalar mode, late fault
    f_faulted = 1 - math.exp(-budget)                       # 0.095
    dn_ferm = a.fermion_fault_dn_max.lo                     # 0.119 per fermion mode (JW string, worst case)
    bias_ferm = {n: f_faulted * dn_ferm / n for n in (0.1, 0.5, 1.0)}   # 11%, 2.3%, 1.1%
    n_ferm_bias_10pct = f_faulted * dn_ferm / 0.1           # 0.113: bias <= 10% above this occupation
    # U2 verifier: shots on the fermion modes, binomial n(1-n) per mode, m' = 12 in the smallest shell
    m_f = a.fermion_shell_modes_2033.lo
    n_ferm_min = 1.0 / (1.0 + base_camp * m_f * a.precision_campaign_2033.lo ** 2)   # 0.053 (same at 30%)
    n_ferm_min_corr = 1.0 / (1.0 + base_camp * a.precision_campaign_2033.lo ** 2)    # 0.4 if m' = 1
    ns = a.scalar_n_max_proxy.lo
    g_s = ((ns + 0.5) / ns) ** 2
    scalar_shots_10pct = n_set * g_s / (m_sh * a.precision_campaign_2033.lo ** 2)  # 6.4e3 per profile
    scalar_eps_at_priced = tuple(math.sqrt(n_set * g_s / (m_sh * N)) for N in spp_camp)   # 0.65, 0.46
    # E27: one machine, every shot serial. The campaign fits the 5-year horizon.
    horizon_s = a.campaign_years.lo * SECONDS_PER_YEAR
    layers = a.t_depth_layers_per_step_2033.value
    t_depth = tuple(l * sc["t_per_rot"] * n_total for l in layers)   # 1.6-5.6e7
    dx = depth_exports(t_total, t_depth, shots_camp, gt, t0)
    dx_first = depth_exports(t_total, t_depth, shots_first, gt, t0)
    hwp = hwp_group(32, sc["eps_rot"], synthesis="rus")
    s = sc["terms"]
    breakdown = _primitives(a, sc, "2033")
    inter = {
        "V": reg["V"],
        "qubits_per_site_scalar": reg["qubits_per_site_scalar"],
        "qubits_per_site_fermion": reg["qubits_per_site_fermion"],
        "system_qubits": reg["system_qubits"],
        "log2_hilbert_dim_scalar": reg["log2_hilbert_dim_scalar"],
        "log2_hilbert_dim_fermion": reg["log2_hilbert_dim_fermion"],
        "log2_hilbert_dim": reg["log2_hilbert_dim"],
        "H_inf_over_Mpl": a.H_inf_over_Mpl.lo,
        "n_anc": reg["n_anc"],
        "n_hadamard": reg["n_hadamard"],
        "lq_sum": reg["lq"],
        "T_phys_Hinf": T_phys,
        "e_folds_matter_like": e_folds,
        "dt_Hinf": dt,
        "omega_max_Hinf": a.omega_max_Hinf_2033.lo,
        "dt_times_omega_max": omega_dt,                             # 0.05 (U2)
        "expansion_time_over_step": 1.0 / dt,                       # 2e3
        "time_ordering_over_trotter_error": 1.0 / a.omega_max_Hinf_2033.lo,
        "stability_omega_dt": a.stability_omega_dt.lo,
        "omega_dt_taken": omega_dt,
        "omega_dt_ramp": omega_dt_ramp,                             # 1.0
        "verlet_freq_err_at_omega_dt_1": freq_err_old,              # 0.0472 (referee U2)
        "verlet_phase_err_at_omega_dt_1": phase_err_old,            # 47 rad
        "verlet_freq_err_taken": freq_err_new,                      # 1.0e-4
        "step_err_exact_reduced": a.exact_step_err_reduced_2033.lo,  # 3.1e-3 (exact_reduced_step, chi register)
        "step_err_extrapolated": a.proxy_step_err_2033.lo * a.exact_over_proxy_step_2033.lo,   # 0.016 x 1.23 = 0.020
        "switch_occupation": n_switch,                              # (5.2e-3, 3.0e-2) per mode, sudden
        "switch_occupation_swept": n_switch_swept,                  # < 2e-5 after a 50-step linear sweep
        "dt_sweep_steps": n_sweep,                                  # 50
        "t_dt_sweep": t_sweep,                                      # 3.6e7 (not priced)
        "dt_sweep_fraction": t_sweep / t_total,                     # 0.24%
        "condensate_amplitude": a.condensate_amplitude_2033.lo,     # 0.31 m_phi
        "condensate_quanta_per_site": a.condensate_amplitude_2033.lo ** 2 * a.lattice_spacing_mphi_2033.lo ** 3 / 2,
        "condensate_boost_rotations": n_boost,                      # 500 (U1)
        "t_condensate_boost": t_boost,                              # 1.5e4
        # per-term counts (app04:84)
        "links": s["links"],                                        # 375
        "scalar_site_strings": s["scalar_site_strings"],            # 7
        "gradient_link_strings": s["gradient_link_strings"],        # 16
        "fermion_site_strings": s["fermion_site_strings"],          # 20
        "pi2_site_strings": s["pi2_site_strings"],                  # 6
        "hop_link_strings": s["hop_link_strings"],                  # 16
        "wilson_hop_nonzeros": wilson_hop_nonzeros(),               # [8, 8, 8]
        "diag_block_rotations": s["diag_block"],                    # 9375
        "pi2_block_rotations": s["pi2_block"],                      # 750
        "hop_block_rotations": s["hop_block"],                      # 6000 (applied twice per step)
        "qft_pair_rotations": s["qft_pair_rot"],                    # 2250
        "qft_pair_controlled_s": s["qft_pair_cs"],                  # 750 -> 2250 T
        "per_step_rotations": sc["per_step_rot"],                   # 24,375
        "per_step_exact_t": s["qft_pair_cs"] * a.t_per_controlled_s.lo,   # 2250
        "n_rot_per_shot": sc["n_rot"],                              # 4.997e8
        "eps_rot": sc["eps_rot"],                                   # 4.47e-6
        "t_per_rotation": sc["t_per_rot"],                          # 29.64
        "n_trotter_steps_raw": n_steps_raw,                         # 2e4
        "t_per_step": t_step,                                       # 7.246e5
        "t_evolution_raw": t_raw,
        "commutator_lever": lever,
        "n_trotter_steps_net": n_steps_net,
        "t_evolution_net": t_evol,                                  # 1.449e10
        "n_prep_steps": n_prep,                                     # 500
        "dt_ramp_Hinf": a.dt_ramp_Hinf_2033.lo,
        "ramp_window_Hinf": n_prep * a.dt_ramp_Hinf_2033.lo,        # 5
        "ramp_window_over_m_phi_period": n_prep * a.dt_ramp_Hinf_2033.lo * a.omega_max_Hinf_2033.lo,  # 500/m_phi
        "t_prep": t_prep,                                           # 3.62e8
        "t_closing": t_closing,                                     # closing D block + centering
        "t_per_shot": t_total,                                      # 1.486e10
        "hard_op_ceiling": ceiling,
        "ratio_to_ceiling": t_total / ceiling,                      # 14.9
        "overshoot_fraction": t_total / ceiling - 1.0,
        "prep_fraction": t_prep / t_evol,                           # 0.025
        "t_vacuum_reference_per_shot": t_vacuum_reference,
        "vacuum_reference_fraction": t_vacuum_reference / t_total,  # 0.024
        "de_sitter_variant_steps": n_steps_ds,                      # 20
        "de_sitter_variant_t": t_ds_variant,                        # 1.19e7 ("~1.2e7")
        "de_sitter_variant_over_envelope": t_ds_variant / ceiling,
        "hwp32_t_over_plain": hwp["t"] / (32 * sc["t_per_rot"]),
        "hwp32_ancilla": hwp["ancilla"],                            # 31 (E26, r22; was 37)
        # R1/R6 shots (r25)
        "n_settings": n_set,                                        # 2
        "occupation_band": nb,                                      # n ~ 1
        "worst_shell_modes": m_sh,                                  # 3
        "shot_margin": (a.shot_margin_2033.lo, a.shot_margin_2033.hi),
        "precision_first": a.precision_first_2033.lo,               # 0.3
        "precision_campaign": a.precision_campaign_2033.lo,         # 0.1
        "shots_per_profile_first_strict": base_first,               # 16.67
        "shots_per_profile_first": spp_first,                       # (17, 34)
        "shots_per_profile_campaign_strict": base_camp,             # 150
        "shots_per_profile_campaign": spp_camp,                     # (150, 300)
        "profiles_first": a.profiles_first_2033.lo,                 # 1
        "profiles_campaign": (a.profiles_campaign_2033.lo, a.profiles_campaign_2033.hi),   # (10, 20)
        "shots_first": shots_first,                                 # (17, 34)
        "shots_campaign": shots_camp,                               # (1500, 6000)
        "shots": shots_camp,
        "t_qft_readout": t_qft_readout,                             # ~4e4: scalar conjugate setting
        "n_fermion_modes": n_ferm_modes,                            # 500
        "n_givens_bound": n_givens,                                 # 124,750
        "t_fermion_readout_bound": t_ferm_readout,                  # ~7.4e6 (not priced)
        "fermion_readout_fraction_bound": t_ferm_readout / t_total, # 5e-4
        "t_gate_s": gt,
        "shot_overhead_s": t0,
        "t_shot_s": t_shot_s,                                       # 1.486e4 s
        "t_shot_hours": t_shot_s / 3600.0,                          # 4.1 h
        "wall_first_result_d": (wall_first[0] / 86400, wall_first[1] / 86400),   # 2.9-5.8 d
        "wall_campaign_yr": (wall_camp[0] / SECONDS_PER_YEAR, wall_camp[1] / SECONDS_PER_YEAR),   # 0.71-2.8 yr
        "wall_time_s": wall_camp,
        "campaign_years": a.campaign_years.lo,
        "campaign_fits_horizon": wall_camp[1] <= horizon_s,         # True
        "subroutine_calls": shots_camp,
        # R9 exports (campaign run), first-result run alongside
        "t_per_shot": dx["t_per_shot"],
        "t_depth_per_shot": dx["t_depth_per_shot"],
        "f_star": dx["f_star"],                                     # (265, 903)
        "floor_wall_s": dx["floor_wall_s"],
        "factories_for_1yr": dx["factories_for_1yr"],
        "wall_serial_s": dx["wall_serial_s"],
        "baseline_ok": dx["baseline_ok"],
        "fits_1yr": dx["fits_1yr"],
        "floor_wall_first_s": dx_first["floor_wall_s"],
        "wall_first_result_s": wall_first,
        "wall_campaign_s": wall_camp,
        "fault_budget_per_shot": budget,
        "eps_l": eps_l_for_budget,                                  # 6.73e-12
        "eps_l_printed": a.eps_l_2033.lo,                           # 7e-12
        "eps_l_for_fault_budget": eps_l_for_budget,
        "expected_faults_at_printed_eps": faults_at_printed,        # 0.104
        "shot_survival": survival,                                  # fault-free fraction 0.905
        "shot_survival_at_printed_eps": survival_at_printed,
        "eps_l_for_survival_target": eps_l_for_target,
        "eps_l_rfi": a.eps_l_rfi.lo,
        "eps_l_rfi_over_requirement": a.eps_l_rfi.lo / eps_l_for_budget,   # 1486 ("three orders")
        "expected_faults_at_rfi": faults_at_rfi,                    # ~149 (U4: no split-circuit route)
        "idle_qubit_cycles": idle_qubit_cycles,                     # 1.6e12
        "idle_qubit_cycles_over_t": idle_over_t,                    # ~105
        "fault_bias_dn_per_fault": dn_per_fault,                    # 0.155 (scalar, late)
        "fault_bias_dn_fermion": dn_ferm,                           # 0.12 (fermion JW string, late)
        "fault_bias_fermion": bias_ferm,                            # {0.1: 0.113, 0.5: 0.023, 1: 0.011}
        "fault_bias_fermion_n_10pct": n_ferm_bias_10pct,            # 0.113
        "fermion_n_min_priced": n_ferm_min,                         # 0.053
        "fermion_n_min_priced_correlated": n_ferm_min_corr,         # 0.4
        "scalar_shots_10pct_per_profile": scalar_shots_10pct,       # 6.4e3
        "scalar_eps_at_priced_shots": scalar_eps_at_priced,         # (0.65, 0.46)
    }
    notes = (
        "LQ: 1000 + (30-50) = 1030-1050, no Hadamard qubits (R6); box 1050. The condensate boost needs no ancilla.",
        f"T/step derived here (app04): D = 125 x 7 + 375 x 16 + 125 x 20 = 9375, B = 125 x 6 = 750, "
        f"C = 375 x 16 = 6000 (x2), QFT pair 2250 rot + 750 controlled-S; {sc['per_step_rot']} rotations x "
        f"{sc['t_per_rot']:.2f} T (eps {sc['eps_rot']:.2e}, N_rot {sc['n_rot']:.0f}) + 2250 T = {t_step:.4g} T.",
        f"Step (U2): m_phi dt = {omega_dt:g} from the step-convergence proxy (verlet_proxy), checked exactly on a "
        f"reduced instance (exact_reduced_step: <= {100 * a.exact_step_err_reduced_2033.lo:.1f}%, "
        f"<= {100 * a.proxy_step_err_2033.lo * a.exact_over_proxy_step_2033.lo:.0f}% scaled to 5^3); the retired m_phi dt = 1 "
        f"had a {100 * freq_err_old:.2f}% Verlet frequency error, {phase_err_old:.0f} rad over m_phi T = 1000.",
        f"T/shot: {n_steps_net:.0f} x {t_step:.4g} = {t_evol:.4g} + {n_prep} ramp steps at 0.01/H_inf = {t_prep:.4g} "
        f"+ boost {n_boost} rot = {t_boost:.3g} + closing {t_closing:.3g} = {t_total:.4g}, "
        f"{t_total / ceiling:.1f}x the 1e9 reference budget (beyond the 1.5x small-overshoot rule).",
        f"Shots per profile (R6): {n_set:g} x ((n+1/2)/n)^2 / (m' eps^2) at n = {nb:g}, m' = {m_sh:g}: "
        f"{base_first:.2f} at 30% -> {spp_first}, {base_camp:.0f} at 10% -> {spp_camp} (x1-2 margin).",
        f"First result (R1), one machine: 1 profile, {shots_first[0]}-{shots_first[1]} shots x {t_shot_s:.4g} s = "
        f"{wall_first[0] / 86400:.2f}-{wall_first[1] / 86400:.2f} d.",
        f"Campaign: {a.profiles_campaign_2033.lo:g}-{a.profiles_campaign_2033.hi:g} profiles, "
        f"{shots_camp[0]}-{shots_camp[1]} shots = {wall_camp[0] / SECONDS_PER_YEAR:.2f}-"
        f"{wall_camp[1] / SECONDS_PER_YEAR:.2f} yr.",
        f"T-depth {dx['t_depth_per_shot'][0]:.3g}-{dx['t_depth_per_shot'][1]:.3g}, F* = "
        f"{dx['f_star'][0]:.0f}-{dx['f_star'][1]:.0f}; floor {dx['floor_wall_s'][0]:.3g}-{dx['floor_wall_s'][1]:.3g} s.",
        f"Readout: scalar conjugate setting +{t_qft_readout:.3g} T; fermion out-basis <= {n_givens} Givens = "
        f"{t_ferm_readout:.3g} T ({100 * t_ferm_readout / t_total:.2f}% of the shot), not priced; vacuum reference "
        f"run {100 * t_vacuum_reference / t_total:.1f}% per shot, not priced.",
        f"eps_l: 0.1 / {t_total:.4g} = {eps_l_for_budget:.3g}, printed '<~7e-12' ({faults_at_printed:.3f} faults); "
        f"{faults_at_rfi:.0f} expected faults per shot at the RFI 1e-8. No split-circuit route (U4).",
        f"G4: idle qubit-cycles {idle_qubit_cycles:.3g} = {idle_over_t:.0f}x N_T; late-fault bias <= "
        f"{100 * f_faulted:.1f}% x {dn_ferm:.2f} / n_k on a fermion shell ({100 * bias_ferm[0.5]:.1f}% at n = 0.5).",
        f"Fermion shots: binomial, m' = {m_f:g}; the priced shots reach target for n_k >= {n_ferm_min:.3f}; scalar "
        f"shells (n <= {ns:g}) read to {scalar_eps_at_priced[1]:.2f}-{scalar_eps_at_priced[0]:.2f}.",
        f"de Sitter variant: 20 steps = {t_ds_variant:.3g} T at its own R-TOL eps ('~1.2e7').",
    ) + ((f"WARNING: prep_fraction = {t_prep / t_evol:.3f}; the chapter prints "
          f"{a.prep_fraction_printed.lo:g}.",)
         if abs(t_prep / t_evol - a.prep_fraction_printed.lo) > 0.0005 else ())
    return Result(era="2033", lq=reg["lq"], hard_ops=(t_total, t_total), breakdown=breakdown,
                  intermediates=inter, shots=shots_camp, wall_time_s=wall_camp,
                  epsilon_l=(eps_l_for_budget, eps_l_for_budget), notes=notes)


def utility(a: Assumptions) -> dict:
    """The utility box (app04:158-164): 3% of the $109M/yr Cosmic Frontier line."""
    per_yr = a.utility_share.lo * a.cosmic_frontier_musd_per_yr.lo      # 3.27 -> "~$3M/yr"
    return {"utility_musd_per_yr": per_yr,
            "utility_musd_per_yr_rounded": round(per_yr),
            "utility_musd_campaign": round(per_yr) * a.campaign_years.lo}


# --------------------------------------------------------------------------- #
# Referee checks U1-U3 (2026-10-04). numpy/scipy are imported lazily: the cost
# model above stays pure Python; only these checks and their tests need them.
# --------------------------------------------------------------------------- #

def switch_occupation(w: float, wh1: float, wh2: float, n_sweep: int = 0) -> float:
    """U2 verifier: quanta left in a mode of frequency w (units m_phi) when the Stormer-Verlet step changes from
    m_phi dt = wh1 to wh2. The kick-drift-kick map conserves p^2 + w_eff^2 x^2, w_eff = w sqrt(1 - (w h)^2 / 4), so its
    vacuum depends on the step. n_sweep = 0 is a sudden switch, (r + 1/r)/4 - 1/2 with r = w_eff(h1) / w_eff(h2);
    otherwise the step is swept linearly from h1 to h2 over n_sweep steps, starting in the h1 map vacuum."""
    def weff(h):
        return w * math.sqrt(1 - (w * h) ** 2 / 4)
    if n_sweep == 0:
        r = weff(wh1) / weff(wh2)
        return (r + 1 / r) / 4 - 0.5
    we = weff(wh1)
    xx, pp, xp = 1 / (2 * we), we / 2, 0.0                  # covariance of the h1 map vacuum
    for i in range(n_sweep):
        h = wh1 + (wh2 - wh1) * i / (n_sweep - 1)
        for m in (((1, 0), (-w * w * h / 2, 1)), ((1, h), (0, 1)), ((1, 0), (-w * w * h / 2, 1))):
            (a11, a12), (a21, a22) = m
            xx, xp, pp = (a11 * a11 * xx + 2 * a11 * a12 * xp + a12 * a12 * pp,
                          a11 * a21 * xx + (a11 * a22 + a12 * a21) * xp + a12 * a22 * pp,
                          a21 * a21 * xx + 2 * a21 * a22 * xp + a22 * a22 * pp)
    we2 = weff(wh2)
    return (pp + we2 ** 2 * xx) / (2 * we2) - 0.5


def fermion_fault_shift(L: int = 5, bm: float = 0.346, r: float = 1.0, js=None) -> float:
    """G4 (verifier): worst shell-averaged occupation added by one X fault on a fermion qubit, in the free Wilson
    vacuum (Hamiltonian form, alpha_i / beta, lattice units, mass b m_psi). In Jordan-Wigner order X_j = P_{<j}
    (c_j + c_j^dag): conjugation flips c_i for i > j (up to a global sign) and swaps c_j <-> c_j^dag. Occupations are
    read in the free out-mode basis. js = qubit positions to try (default all)."""
    import itertools
    import numpy as np
    s1 = [np.array([[0, 1], [1, 0]]), np.array([[0, -1j], [1j, 0]]), np.diag([1, -1])]
    beta = np.diag([1, 1, -1, -1]).astype(complex)
    alpha = [np.block([[np.zeros((2, 2)), s], [s, np.zeros((2, 2))]]) for s in s1]
    V = L ** 3
    idx = lambda x, y, z: ((x % L) * L + (y % L)) * L + (z % L)
    h = np.zeros((4 * V, 4 * V), complex)
    for x, y, z in itertools.product(range(L), repeat=3):
        A = slice(4 * idx(x, y, z), 4 * idx(x, y, z) + 4)
        h[A, A] += beta * (bm + 3 * r)
        for d, (dx, dy, dz) in enumerate(((1, 0, 0), (0, 1, 0), (0, 0, 1))):
            B = slice(4 * idx(x + dx, y + dy, z + dz), 4 * idx(x + dx, y + dy, z + dz) + 4)
            h[A, B] += -0.5j * alpha[d] - 0.5 * r * beta
            h[B, A] += 0.5j * alpha[d] - 0.5 * r * beta
    E, U = np.linalg.eigh(h)
    neg, pos = U[:, E < 0], U[:, E > 0]
    G = neg.conj() @ neg.T                                   # <c_a^dag c_b>
    shells = np.round(E[E > 0], 6)
    labels = np.unique(shells)
    N = 4 * V
    worst = 0.0
    for j in (range(N) if js is None else js):
        s = np.where(np.arange(N) < j, 1.0, -1.0)
        s[j] = 1.0
        Gp = np.outer(s, s) * G
        Gp[j, :] = 0
        Gp[:, j] = 0
        Gp[j, j] = 1 - G[j, j]
        dn = np.real(np.einsum("am,ab,bm->m", pos, Gp, pos.conj()))
        worst = max(worst, max(dn[shells == e].mean() for e in labels))
    return worst


def digitized_site(K: int) -> dict:
    """One site oscillator (unit frequency and mass) on the K-point grid of arXiv:2108.10793
    (Eq. at bosons.tex:830, Delta_phi = sqrt(2 pi / K), conjugate grid by the centered DFT, :841).
    Returns the gap, vacuum <x^2> (continuum 1/2), and the highest level above the ground state."""
    import numpy as np
    from scipy.linalg import eigh
    dx = math.sqrt(2 * math.pi / K)
    j = np.arange(K) - (K - 1) / 2
    x = j * dx
    F = np.exp(-1j * np.outer(j, j) * 2 * np.pi / K) / np.sqrt(K)
    P2 = F.conj().T @ np.diag(x ** 2) @ F
    w, v = eigh(0.5 * P2 + np.diag(0.5 * x ** 2))
    return dict(gap=w[1] - w[0], x2=float(np.sum(abs(v[:, 0]) ** 2 * x ** 2)),
                e_max_above_ground=w[-1] - w[0], x=x, P2=P2, ground=v[:, 0])


def coherent_fidelity(K: int, alpha2: float, periods: float = 100, dt: float = 0.05) -> float:
    """U1: boost the digitized vacuum by exp(i pbar x), pbar = sqrt(2 alpha2) (diagonal in the field basis), evolve
    under the digitized H, and return the minimum fidelity with the continuum coherent state over `periods`."""
    import numpy as np
    from scipy.linalg import expm
    d = digitized_site(K)
    x = d["x"]
    U = expm(-1j * (0.5 * d["P2"] + np.diag(0.5 * x ** 2)) * dt)
    pbar = math.sqrt(2 * alpha2)
    st = np.exp(1j * pbar * x) * d["ground"]
    worst, t = 1.0, 0.0
    for n in range(int(2 * math.pi * periods / dt) + 1):
        if n % 63 == 0:
            ref = np.exp(-(x - pbar * math.sin(t)) ** 2 / 2 + 1j * pbar * math.cos(t) * x)
            worst = min(worst, abs(np.vdot(ref / np.linalg.norm(ref), st)) ** 2)
        st = U @ st
        t += dt
    return worst


def verlet_proxy(mdt: float, lam: float, Phi0: float, g: float, m_psi: float = 0.1,
                 ks=None, H_over_m: float = 0.01, T_m: float = 1000.0):
    """U2: mean-field / linear proxy of the 2033 Trotter step (units m_phi = 1).

    The symmetric splitting D/2 B D/2 acts on a quadratic Hamiltonian exactly as the classical kick-drift-kick map,
    so this IS the circuit's action on (i) the condensate (mean field of the boosted state), (ii) linear scalar
    fluctuations about it (mass^2 1 + 3 lam phibar^2) and (iii) free fermion modes with mass m_psi + g phibar and
    hop k/a. Coefficients use a at the step midpoint and phibar at the kick time; a(t) = (1 + 1.5 H t)^{2/3}.
    The condensate starts at a zero crossing with velocity Phi0 (the boost). Vacuum and readout are those of the
    step map itself (the ramp at the same step prepares the map's vacuum; a prep-only run calibrates the readout).
    Returns (n_scalar, n_fermion) per k at T. Default k: the 5^3 shells at b = 3.46/m_phi."""
    import numpy as np
    if ks is None:
        ks = np.sqrt(np.array([1, 2, 3, 4, 5, 6, 8, 9.])) * 2 * np.pi / 17.3
    ks = np.asarray(ks, float)
    N = int(round(T_m / mdt))
    h = T_m / N
    a_of = lambda t: (1 + 1.5 * H_over_m * t) ** (2 / 3)

    def ferm_vac(M, k, sign):
        e = np.exp(-1j * M * h / 2)
        c, s = np.cos(k * h), np.sin(k * h)
        out = np.empty((len(k), 2), complex)
        for i in range(len(k)):
            U = np.array([[e * c[i] * e, -1j * s[i]], [-1j * s[i], np.conj(e) * c[i] * np.conj(e)]])
            w, V = np.linalg.eig(U)
            q = -np.angle(w) / h
            out[i] = V[:, np.argmin(q) if sign < 0 else np.argmax(q)]
        return out / np.linalg.norm(out, axis=1)[:, None]

    ph, pp = 0.0, Phi0
    W0 = np.sqrt(ks ** 2 + 1.0)
    Om = W0 * np.sqrt(1 - (W0 * h) ** 2 / 4)
    u = 1 / np.sqrt(2 * Om) + 0j
    pu = -1j * Om * u
    v = ferm_vac(m_psi, ks, -1)
    t = 0.0
    for _ in range(N):
        am = a_of(t + h / 2)
        a3 = am ** 3
        for half in (0, 1):
            pp -= a3 * (ph + lam * ph ** 3) * h / 2
            pu = pu - (a3 * (1.0 + 3 * lam * ph ** 2) + am * ks ** 2) * u * h / 2
            e = np.exp(-1j * (m_psi + g * ph) * h / 2)
            v = v * np.array([e, np.conj(e)])[None, :]
            if half == 0:
                ph += pp / a3 * h
                u = u + pu / a3 * h
                c, s = np.cos(ks / am * h), np.sin(ks / am * h)
                v = np.stack([c * v[:, 0] - 1j * s * v[:, 1], -1j * s * v[:, 0] + c * v[:, 1]], 1)
        t += h
    aT = a_of(T_m)
    a3 = aT ** 3
    W = np.sqrt(1.0 + 3 * lam * ph ** 2 + ks ** 2 / aT ** 2)
    Om = W * np.sqrt(1 - (W * h) ** 2 / 4)
    ns = (abs(pu) ** 2 / (2 * a3) + a3 * Om ** 2 * abs(u) ** 2 / 2) / Om - 0.5
    M = m_psi + g * ph
    kk = ks / aT
    e = np.exp(-1j * M * h / 2)
    c, s = np.cos(kk * h), np.sin(kk * h)
    wv = np.empty((len(ks), 2), complex)
    for i in range(len(ks)):
        U = np.array([[e * c[i] * e, -1j * s[i]], [-1j * s[i], np.conj(e) * c[i] * np.conj(e)]])
        w, V = np.linalg.eig(U)
        wv[i] = V[:, np.argmax(-np.angle(w) / h)]
    nf = abs(np.sum(np.conj(wv) * v, 1)) ** 2
    return ns, nf



def fixed_grid_gap(K: int, a3: float) -> float:
    """U2 exact check: gap (units m_phi) of one site on the K-point grid set at a = 1 when the site Hamiltonian is
    p^2/(2 a^3) + a^3 x^2/2 (x = sqrt(m b^3) phi, fixed grid). The site vacuum narrows as a^{-3/2}; the grid does not."""
    import numpy as np
    from scipy.linalg import eigh
    d = digitized_site(K)
    w = eigh(0.5 * d["P2"] / a3 + np.diag(0.5 * a3 * d["x"] ** 2), eigvals_only=True)
    return float(w[1] - w[0])


def exact_reduced_step(lam: float, g: float, hs, enc: str = "chi", T_m: float = 1000.0, rec=(1000.0,),
                       L: int = 3, K: int = 16, b: float = 3.46, Phi0: float = 0.31, m_psi: float = 0.1,
                       H_over_m: float = 0.01, r: float = 1.0):
    """U2 exact step check on a reduced instance (units m_phi = 1). Exact many-body evolution of a 1D periodic ring
    of L sites at the chapter's b, with the site normalization v = b^3 of the 5^3 instance (so per-site quanta,
    digitization and couplings match), a real scalar on the K-point grid of arXiv:2108.10793 and one 2-component
    Wilson fermion (gamma0 = s3, gamma0 gamma1 = s1) at half filling (L = 3: dim 20 x 16^3 = 81920).
    Step = the chapter's symmetric splitting D(h/2) [Kin(h) Hop(h)] D(h/2), coefficients at a(t + h/2); D holds every
    term diagonal in the field/occupation basis (mass, phi^4, gradient, Wilson on-site, m_psi, Yukawa), Kin the per-site
    pi^2 (exact via the centered DFT), Hop the fermion hop (exact, it acts on the other register).
    enc = "fixed": phi on the grid set at a = 1. enc = "chi": the comoving field a^{3/2} phi, which adds
    gamma (x p + p x)/2 with gamma = (3/2) adot/a, applied as D/2 Dil/2 [Kin Hop] Dil/2 D/2 (explicit squeeze).
    enc = "shear": the same H written as p^2/2 + gamma(xp+px)/2 = U (p^2/2) U^dag - gamma^2 x^2/2, U = exp(-i gamma
    x^2/2): Kin -> U Kin U^dag and -gamma^2 x^2/2 joins D. The x^2 phases are ZZ angle changes on the phi^2
    rotations, so this is the priced circuit (U2 verifier).
    Initial state: ground state of the digitized H at a = 1, then the boost exp(i sqrt(v) Phi0 sum_x x_x).
    Readout at each time in rec: fermion particle number per |k| (positive-energy out-modes of h(k) at
    M = m_psi + g phibar) and the scalar quasiparticle n at k = 2 pi/(L b). Returns {h: [dict per rec time]}."""
    import numpy as np
    from scipy.sparse.linalg import eigsh, LinearOperator
    v = b ** 3
    a_of = lambda t: (1 + 1.5 * H_over_m * t) ** (2 / 3)
    dx = math.sqrt(2 * math.pi / K)
    jj = np.arange(K) - (K - 1) / 2
    X = jj * dx
    F = np.exp(-1j * np.outer(jj, jj) * 2 * np.pi / K) / np.sqrt(K)
    P2 = F.conj().T @ np.diag(X ** 2) @ F
    Pm = F.conj().T @ np.diag(X) @ F
    G = (np.diag(X) @ Pm + Pm @ np.diag(X)) / 2
    gw, gV = np.linalg.eigh((G + G.conj().T) / 2)
    NM = 2 * L
    states = [st for st in range(2 ** NM) if bin(st).count("1") == L]
    sidx = {st: i for i, st in enumerate(states)}
    NF = len(states)

    def cdag_c(i, j):
        M = np.zeros((NF, NF))
        for st in states:
            if not (st >> j) & 1:
                continue
            s1_ = st ^ (1 << j)
            sg = (-1) ** bin(st & ((1 << j) - 1)).count("1")
            if (s1_ >> i) & 1:
                continue
            sg *= (-1) ** bin(s1_ & ((1 << i) - 1)).count("1")
            M[sidx[s1_ | (1 << i)], sidx[st]] += sg
        return M
    CC = [[cdag_c(i, j) for j in range(NM)] for i in range(NM)]
    nmode = np.array([[(st >> i) & 1 for i in range(NM)] for st in states], float)
    cbc = np.stack([nmode[:, 2 * x] - nmode[:, 2 * x + 1] for x in range(L)], 1)
    s1 = np.array([[0, 1], [1, 0.]])
    s3 = np.array([[1, 0], [0, -1.]])
    T = np.roll(np.eye(L), 1, axis=1)
    h_hop1 = (-(r / 2) * np.kron(T + T.T, s3) + np.kron((T - T.T) / 2j, s1)) / b
    H_hop1 = sum(h_hop1[i, j] * CC[i][j] for i in range(NM) for j in range(NM))
    hw, hV = np.linalg.eigh((H_hop1 + H_hop1.conj().T) / 2)
    ks = 2 * np.pi * np.arange(L) / (L * b)
    ks[ks > np.pi / b] -= 2 * np.pi / b
    xs = np.stack(np.meshgrid(*([X] * L), indexing="ij"), 0).reshape(L, -1)
    grad = sum((xs[x] - xs[(x + 1) % L]) ** 2 for x in range(L)) / (2 * b ** 2)
    quad, quart = (xs ** 2).sum(0) / 2, (xs ** 4).sum(0) / (4 * v)
    wil, yuk = cbc.sum(1), cbc @ xs / math.sqrt(v)
    shape = (NF,) + (K,) * L

    def diag(a):
        if enc == "fixed":
            sc, yk = a ** 3 * (quad + lam * quart) + a * grad, yuk
        else:
            sc, yk = quad + lam * quart / a ** 3 + grad / a ** 2, yuk / a ** 1.5
        return sc[None, :] + ((m_psi + r / (a * b)) * wil)[:, None] + g * yk

    def site(psi, U, add=False):
        out = np.zeros_like(psi) if add else psi
        for ax in range(1, L + 1):
            y = np.moveaxis(np.tensordot(U, psi if add else out, axes=([1], [ax])), 0, ax)
            out = out + y if add else y
        return out
    eig = lambda w_, V_, th: V_ @ (np.exp(-1j * th * w_)[:, None] * V_.conj().T)
    kinU = lambda th: F.conj().T @ (np.exp(-1j * th * X ** 2)[:, None] * F)

    d0 = diag(1.0).ravel()
    op = LinearOperator((d0.size, d0.size), dtype=complex, matvec=lambda q: d0 * q + (
        site(q.reshape(shape), P2 / 2, add=True) + np.tensordot(H_hop1, q.reshape(shape), axes=([1], [0]))).ravel())
    psi0 = eigsh(op, k=1, which="SA", tol=1e-12)[1][:, 0].reshape(shape)
    psi0 = psi0 * np.exp(1j * math.sqrt(v) * Phi0 * xs.sum(0)).reshape((1,) + (K,) * L)

    def hk(k, a, M):
        return (M + r / (a * b) * (1 - np.cos(k * b))) * s3 + np.sin(k * b) / (a * b) * s1

    def observe(psi, t):
        a = a_of(t)
        ps = (abs(psi) ** 2).sum(0).reshape(-1)
        xm = xs @ ps
        Cx = (xs * ps) @ xs.T - np.outer(xm, xm)
        pq = (abs(site(psi, F)) ** 2).sum(0).reshape(-1)
        pm = xs @ pq
        Cp = (xs * pq) @ xs.T - np.outer(pm, pm)
        phibar = xm.mean() / math.sqrt(v) / (1.0 if enc == "fixed" else a ** 1.5)
        sx, sp = (a ** 3, a ** -3) if enc == "fixed" else (1.0, 1.0)
        k1 = 2 * math.pi / (L * b)
        ph = np.exp(-1j * k1 * b * np.arange(L))
        W = math.sqrt(1 + 3 * lam * phibar ** 2 + 4 * math.sin(k1 * b / 2) ** 2 / (a * b) ** 2)
        ns = (sp * np.real(ph.conj() @ Cp @ ph) / L + sx * W ** 2 * np.real(ph.conj() @ Cx @ ph) / L) / (2 * W) - 0.5
        Mf = psi.reshape(NF, -1)
        RHO = Mf @ Mf.conj().T
        rho = np.array([[np.trace(CC[i][j] @ RHO) for j in range(NM)] for i in range(NM)])
        nf = {}
        for k in ks:
            up = np.linalg.eigh(hk(k, a, m_psi + g * phibar))[1][:, 1]
            w_ = np.kron(np.exp(1j * k * b * np.arange(L)) / math.sqrt(L), up)
            nf.setdefault(round(abs(k), 8), []).append(float(np.real(w_ @ rho @ w_.conj())))
        return dict(t=t, phibar=float(phibar), n_scalar=float(ns),
                    n_ferm=[float(np.mean(q)) for _, q in sorted(nf.items())])

    out = {}
    for h in hs:
        psi, t = psi0.copy(), 0.0
        recn = {int(round(tr / h)): tr for tr in rec}
        res = []
        for n in range(1, int(round(T_m / h)) + 1):
            am = a_of(t + h / 2)
            ph = np.exp(-1j * (h / 2) * diag(am)).reshape(shape)
            gm = 1.5 * H_over_m / (1 + 1.5 * H_over_m * (t + h / 2))
            if enc == "shear":
                ph = ph * np.exp(1j * (h / 2) * gm ** 2 * (xs ** 2).sum(0) / 2).reshape((1,) + (K,) * L)
            psi = psi * ph
            if enc == "chi":
                Ud = eig(gw, gV, gm * h / 2)
                psi = site(site(psi, Ud), kinU(h / 2))
            elif enc == "shear":
                sh = np.exp(-1j * gm * X ** 2 / 2)
                psi = site(psi, sh[:, None] * kinU(h / 2) * np.conj(sh)[None, :])
            else:
                psi = site(psi, kinU(h / (2 * am ** 3)))
            psi = np.tensordot(eig(hw, hV, h / am), psi, axes=([1], [0]))
            if enc == "chi":
                psi = site(psi, Ud)
            psi = psi * ph
            t = n * h
            if n in recn:
                res.append(observe(psi, recn[n]))
        out[h] = res
    return out


def _ds2028_setup(m: float, L: int):
    """Ring of L sites at the 2028 spacing (3D omega_max = 2 H_inf at a = 1: m^2 + 12/b^2 = 4), units H_inf = 1.
    Returns b, the real orthonormal Fourier basis E (columns k = 0, 1, ...), lattice momenta, and the in-state
    covariance: Bunch-Davies (massive dS mode functions at the lattice momentum, eta_in = -1) for k != 0, the a = 1
    ground state for k = 0 (IR-singular for BD). Per-site x = sqrt(v) phi; v drops out of a free field."""
    import numpy as np
    from scipy.special import hankel1, h1vp
    b = math.sqrt(12.0 / (4.0 - m * m))
    n_ = np.arange(L)
    kt = 2 / b * np.abs(np.sin(np.pi * n_ / L))
    cols = []
    for n in range(L):
        if n == 0:
            cols.append(np.ones(L) / math.sqrt(L))
        elif 2 * n < L:
            cols.append(math.sqrt(2 / L) * np.cos(2 * np.pi * n * n_ / L))
        elif 2 * n == L:
            cols.append(np.cos(np.pi * n_) / math.sqrt(L))
        else:
            cols.append(math.sqrt(2 / L) * np.sin(2 * np.pi * n * n_ / L))
    E = np.stack(cols, 1)
    nu = math.sqrt(9 / 4 - m * m)
    ph = np.exp(1j * math.pi * (nu + 0.5) / 2)
    sxx, sxp, spp = np.zeros(L), np.zeros(L), np.zeros(L)
    for n in range(L):
        if n == 0:
            sxx[n], spp[n] = 1 / (2 * m), m / 2
            continue
        z = kt[n]
        q = math.sqrt(math.pi) / 2 * hankel1(nu, z) * ph
        p = math.sqrt(math.pi) / 2 * ph * (-1.5 * hankel1(nu, z) - z * h1vp(nu, z))   # a^3 dq/dt at t = 0
        sxx[n], spp[n], sxp[n] = abs(q) ** 2, abs(p) ** 2, float(np.real(q * np.conj(p)))
    C = np.block([[E @ np.diag(sxx) @ E.T, E @ np.diag(sxp) @ E.T], [E @ np.diag(sxp) @ E.T, E @ np.diag(spp) @ E.T]])
    T = np.roll(np.eye(L), 1, axis=1)
    lap = (2 * np.eye(L) - T - T.T) / b ** 2
    return b, E, C, lap


def _ds2028_blocks(m: float, lap, L: int, h: float, n: int):
    """The chapter's 2028 step K(h/2) D(h) K(h/2), coefficients at a(t + h/2), a = e^t (app04:82)."""
    import numpy as np
    seq = []
    for s in range(n):
        am = math.exp(s * h + h / 2)
        seq += [("K", h / 2 / am ** 3), ("D", h * (am ** 3 * m * m * np.eye(L) + am * lap)), ("K", h / 2 / am ** 3)]
    return seq


def _shell_err(Cd, Cr, E, L):
    """Max relative error of <phi_k phi_k>, <pi_k pi_k>, and of Re<phi_k pi_k> (in units of their geometric mean)."""
    out = []
    for n in range(L):
        e = E[:, n]
        d = (e @ Cd[:L, :L] @ e, e @ Cd[L:, L:] @ e, e @ Cd[:L, L:] @ e)
        r = (e @ Cr[:L, :L] @ e, e @ Cr[L:, L:] @ e, e @ Cr[:L, L:] @ e)
        out.append(max(abs(d[0] / r[0] - 1), abs(d[1] / r[1] - 1), abs(d[2] - r[2]) / math.sqrt(r[0] * r[1])))
    return out


def ds2028_grid_check(m: float, K: int, frame: str = "adapted", n_steps: int = 4, L: int = 4, h: float = 0.5,
                      merged: bool = False):
    """U2 for 2028 (2026-10-05): exact digitized evolution of the 2028 free scalar in de Sitter on a reduced instance,
    against the continuum Gaussian evolved by the same map (the chapter's reference).

    Ring of L = 4 sites at the 2028 spacing (its lowest shell and frequency are those of the 4^3 IR shell), K-point
    grid of arXiv:2108.10793 (spacing sqrt(2 pi/K), centered DFT), state = the in-state Gaussian sampled on the grid.
    Register frame z' = S z per site, S lower triangular (dilation + shear):
      'fixed'   : S balanced on the in-state site marginal, held for the whole run (the chapter's grid);
      'adapted' : S re-balanced on the reference site covariance at each kinetic-block boundary.
    Each kinetic half-step in the frame is S_out Up(u) S_in^-1 = Lx(c1) Up(u') Lx(c2): x^2 and p^2 phases on the
    existing phi^2 / pi^2 rotations, so no rotation is added; D in a lower-triangular frame is an x-phase with
    Q / s^2. Returns, after each step, the list of per-shell errors (index 1 = the IR shell).
    merged = True joins the two pi^2 half-steps between neighboring steps into one kinetic block (n + 1 QFT pairs
    for n steps, the priced circuit) and returns only the final-time errors (one-element list)."""
    import numpy as np
    b, E, C, lap = _ds2028_setup(m, L)
    seq = _ds2028_blocks(m, lap, L, h, n_steps)
    if merged:
        mseq = []
        for blk in seq:
            if mseq and blk[0] == "K" and mseq[-1][0] == "K":
                mseq[-1] = ("K", mseq[-1][1] + blk[1])
            else:
                mseq.append(blk)
        seq = mseq
    covs = [C]
    for typ, v in seq:
        M = np.eye(2 * L)
        if typ == "K":
            M[:L, L:] = v * np.eye(L)
        else:
            M[L:, :L] = -v
        covs.append(M @ covs[-1] @ M.T)

    def bal(Cm):
        sxx, sxp, spp = Cm[0, 0], Cm[0, L], Cm[L, L]
        sc = math.sqrt(math.sqrt(sxx * spp - sxp * sxp) / sxx)
        return np.array([[sc, 0.0], [-sxp / (sc * sxx), 1 / sc]])
    S0 = bal(C)
    frames = []
    for i, Cm in enumerate(covs):
        if i > 0 and seq[i - 1][0] == "D":
            frames.append(frames[-1])
        else:
            frames.append(S0 if frame == "fixed" else bal(Cm))
    big = lambda S: np.kron(S, np.eye(L))
    dx = math.sqrt(2 * math.pi / K)
    jj = np.arange(K) - (K - 1) / 2
    X = jj * dx
    F = np.exp(-1j * np.outer(jj, jj) * 2 * np.pi / K) / np.sqrt(K)
    Pm = F.conj().T @ np.diag(X) @ F
    xs = np.stack(np.meshgrid(*([X] * L), indexing="ij"), 0).reshape(L, -1)
    shape = (K,) * L
    Cr = big(S0) @ C @ big(S0).T
    M0 = ((Cr[L:, :L] + Cr[:L, L:].T) / 2 + 0.5j * np.eye(L)) @ np.linalg.inv(Cr[:L, :L])
    psi = np.exp(0.5j * np.einsum("in,ij,jn->n", xs, M0, xs))
    psi /= np.linalg.norm(psi)

    def site(ps, U, only=None):
        ps = ps.reshape(shape)
        for ax in (range(L) if only is None else (only,)):
            ps = np.moveaxis(np.tensordot(U, ps, axes=([1], [ax])), 0, ax)
        return ps.reshape(-1)
    xph = lambda ps, Q: ps * np.exp(-0.5j * np.einsum("in,ij,jn->n", xs, Q, xs))

    def cov(ps):
        pr = abs(ps) ** 2
        pq = abs(site(ps, F)) ** 2
        Cxp = np.array([[np.real(np.vdot(ps, xs[i] * site(ps, Pm, j))) for j in range(L)] for i in range(L)])
        return np.block([[(xs * pr) @ xs.T, Cxp], [Cxp.T, (xs * pq) @ xs.T]])
    out = []
    for i, (typ, v) in enumerate(seq):
        Sin, Sout = frames[i], frames[i + 1]
        if typ == "D":
            psi = xph(psi, v / Sin[0, 0] ** 2)
            continue
        A = Sout @ np.array([[1.0, v], [0.0, 1.0]]) @ np.linalg.inv(Sin)
        up = A[0, 1]
        psi = xph(psi, (1 - A[0, 0]) / up * np.eye(L))
        psi = site(psi, F.conj().T @ np.diag(np.exp(-0.5j * up * X ** 2)) @ F)
        psi = xph(psi, (1 - A[1, 1]) / up * np.eye(L))
        if (i == len(seq) - 1) if merged else (i % 3 == 2):
            Si = big(np.linalg.inv(Sout))
            out.append(_shell_err(Si @ cov(psi) @ Si.T, covs[i + 1], E, L))
    return out


def ds2028_map_vs_exact(m: float, h: float = 0.5, n_steps: int = 4, L: int = 4) -> float:
    """Largest relative deviation, over the ring's shells, of the continuum n-step map from exact de Sitter evolution
    (the full linear ODE) at T_phys = n h. The 2028 reference is the map itself (app04:74); this is its distance from
    de Sitter physics."""
    import numpy as np
    from scipy.integrate import solve_ivp
    b, E, C, lap = _ds2028_setup(m, L)
    I, Z = np.eye(L), np.zeros((L, L))

    def f(t, y):
        a = math.exp(t)
        A = np.block([[Z, I / a ** 3], [-(a ** 3 * m * m * I + a * lap), Z]])
        return (A @ y.reshape(2 * L, 2 * L)).ravel()
    sol = solve_ivp(f, (0, n_steps * h), np.eye(2 * L).ravel(), rtol=1e-12, atol=1e-13, method="DOP853")
    Mx = sol.y[:, -1].reshape(2 * L, 2 * L)
    Cm = C
    for typ, v in _ds2028_blocks(m, lap, L, h, n_steps):
        M = np.eye(2 * L)
        if typ == "K":
            M[:L, L:] = v * np.eye(L)
        else:
            M[L:, :L] = -v
        Cm = M @ Cm @ M.T
    return max(_shell_err(Cm, Mx @ C @ Mx.T, E, L))


def twopi_lo(lam: float, g: float, rec=(1000.0,), L: int = 3, b: float = 3.46, Phi0: float = 0.31,
             m_psi: float = 0.1, H_over_m: float = 0.01, r: float = 1.0, rtol: float = 1e-9):
    """G6 (2026-10-05): Gaussian truncation (below the two-loop 2PI) on the U2 reduced instance (same lattice, register, a(t), boost
    and readout as exact_reduced_step enc = 'chi', continuum scalar). Hartree for lam phi^4 (Gaussian scalar),
    mean-field Yukawa (Slater-determinant fermion in the field of xbar) with the fermion tadpole back-reacting on the
    condensate. Scalar-fermion scattering (the two-loop 2PI diagrams of arXiv:hep-ph/0212404) is not included.
    Initial state: self-consistent Gaussian/Hartree-Fock ground state at a = 1, then the boost pbar += sqrt(v) Phi0.
    Returns [dict(t, phibar, n_scalar, n_ferm per |k|)] at the times in rec."""
    import numpy as np
    from scipy.integrate import solve_ivp
    v = b ** 3
    a_of = lambda t: (1 + 1.5 * H_over_m * t) ** (2 / 3)
    gam_of = lambda t: 1.5 * H_over_m / (1 + 1.5 * H_over_m * t)
    s1 = np.array([[0, 1], [1, 0.]])
    s3 = np.array([[1, 0], [0, -1.]])
    T = np.roll(np.eye(L), 1, axis=1)
    h_hop1 = (-(r / 2) * np.kron(T + T.T, s3) + np.kron((T - T.T) / 2j, s1)) / b
    sgn = np.kron(np.eye(L), s3)
    ks = 2 * np.pi * np.arange(L) / (L * b)
    ks[ks > np.pi / b] -= 2 * np.pi / b
    w2g = 4 * np.sin(ks * b / 2) ** 2 / b ** 2
    NM = 2 * L
    hsp = lambda a, x: h_hop1 / a + (m_psi + r / (a * b) + g * x / (math.sqrt(v) * a ** 1.5)) * sgn
    pbp = lambda rho: float(np.real(np.trace(sgn @ rho))) / L
    xbar, D = 0.0, 0.5
    for _ in range(500):
        U = np.linalg.eigh(hsp(1.0, xbar))[1][:, :L]
        rho = U.conj() @ U.T                                   # <c_i^dag c_j>
        Dn = float(np.mean(1 / (2 * np.sqrt(1 + 3 * lam * (xbar ** 2 + D) / v + w2g))))
        xn = xbar
        for _ in range(50):
            xn -= (xn + lam * (xn ** 3 + 3 * xn * Dn) / v + g * pbp(rho) / math.sqrt(v)) / (1 + lam * (3 * xn ** 2 + 3 * Dn) / v)
        if abs(xn - xbar) < 1e-13 and abs(Dn - D) < 1e-13:
            break
        xbar, D = 0.5 * (xbar + xn), Dn
    w2 = 1 + 3 * lam * (xbar ** 2 + D) / v + w2g
    G = np.stack([np.diag([1 / (2 * np.sqrt(q)), np.sqrt(q) / 2]) for q in w2])
    y0 = np.concatenate([[xbar, math.sqrt(v) * Phi0], G.ravel(), rho.ravel().view(float)])
    unpack = lambda y: (y[0], y[1], y[2:2 + 4 * L].reshape(L, 2, 2),
                        np.ascontiguousarray(y[2 + 4 * L:]).view(complex).reshape(NM, NM))

    def f(t, y):
        a, gm = a_of(t), gam_of(t)
        x, p, G, rho = unpack(y)
        D = float(np.mean(G[:, 0, 0]))
        dp = -x - lam * (x ** 3 + 3 * x * D) / (v * a ** 3) - gm * p - g * pbp(rho) / (math.sqrt(v) * a ** 1.5)
        w2 = 1 + 3 * lam * (x ** 2 + D) / (v * a ** 3) + w2g / a ** 2
        dG = np.stack([(lambda A: A @ G[i] + G[i] @ A.T)(np.array([[gm, 1.0], [-w2[i], -gm]])) for i in range(L)])
        h = hsp(a, x)
        drho = 1j * (h.T @ rho - rho @ h.T)
        return np.concatenate([[p + gm * x, dp], dG.ravel(), drho.ravel().view(float)])
    sol = solve_ivp(f, (0, max(rec)), y0, t_eval=sorted(rec), rtol=rtol, atol=1e-11, method="DOP853")
    out = []
    for j, t in enumerate(sol.t):
        x, p, G, rho = unpack(sol.y[:, j])
        a = a_of(t)
        phibar = x / math.sqrt(v) / a ** 1.5
        W = math.sqrt(1 + 3 * lam * phibar ** 2 + w2g[1] / a ** 2)
        M = m_psi + g * phibar
        nf = {}
        for k in ks:
            hk = (M + r / (a * b) * (1 - np.cos(k * b))) * s3 + np.sin(k * b) / (a * b) * s1
            w_ = np.kron(np.exp(1j * k * b * np.arange(L)) / math.sqrt(L), np.linalg.eigh(hk)[1][:, 1])
            nf.setdefault(round(abs(k), 8), []).append(float(np.real(w_ @ rho @ w_.conj())))
        out.append(dict(t=float(t), phibar=float(phibar), n_scalar=float((G[1, 1, 1] + W ** 2 * G[1, 0, 0]) / (2 * W) - 0.5),
                        n_ferm=[float(np.mean(q)) for _, q in sorted(nf.items())]))
    return out


PUBLISHED = {
    "2028": Published(lq=(220, 240), hard_ops=(9.5e4, 9.5e4),
                      src="app04:2028 box ('220--240', '9.5e4 T-gates', wall '43 s'); requirements "
                          "summary. Smaller ruling (b) (2026-10-05), option (c): one macro-step at K = 8 "
                          "(6.475e4 T) + 3e4 prep = 9.475e4 (was four steps at K = 4, 9.437e4, 160-180 LQ).",
                      rel_tol=0.10),
    "2033": Published(lq=(1050, 1050), hard_ops=(1.5e10, 1.5e10),
                      src="app04:2033 box ('1050', '1.5e10 T-gates'); requirements summary. Referee U2 "
                          "(2026-10-04): m_phi dt = 0.05; 2e4 evolution + 500 ramp steps x 7.246e5 + boost + "
                          "closing block = 1.486e10, 14.9x the 1e9 reference budget (was 1.0e9 at m_phi dt = 1). "
                          "LQ 1000+30-50 = 1030-1050 (prose, no Hadamard qubits, R6).",
                      rel_tol=0.10),
}


def INSTANCE_ROWS(a: Assumptions, era: str, r: Result):
    """One costed instance per era (each box quotes one)."""
    if era == "2028":
        return [(r"free scalar dS $4^3$", r.lq, r.hard_ops,
                 {"shots": r.shots[0], "wall_time_s": r.wall_time_s[0],
                  "t_evolution": r.intermediates["t_evolution"], "t_prep": r.intermediates["t_prep"]})]
    if era == "2033":
        i = r.intermediates
        return [(r"$\phi^4$+Wilson preheating $5^3$", r.lq, r.hard_ops,
                 {"shots": r.shots[0], "wall_time_s": r.wall_time_s[0],
                  "shots_hi": r.shots[1], "wall_time_s_hi": r.wall_time_s[1],
                  "shots_first": i["shots_first"], "wall_first_result_s": i["wall_first_result_s"],
                  "t_evolution_raw": r.intermediates["t_evolution_raw"],
                  "t_evolution_net": r.intermediates["t_evolution_net"],
                  "t_prep": r.intermediates["t_prep"],
                  "conditional": False,     # U2: step checked exactly; the comoving-register squeeze is free (shear)
                  "t_bound": None})]
    return []
