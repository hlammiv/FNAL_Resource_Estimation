# HISTORICAL MODEL: superseded 2028/2033 instances; only codesign scaling remains active.
"""Ch. 2 — Neutrino–nucleus response (DUNE). Reproduces the LQ and hard-op numbers of
applications/app01_neutrino_nucleus.tex.

Three eras, three routes:
  2028      route (i)   second-quantized chiral EFT, JW, Trotter; DERIVED in prose
                        (N_orb^4/2 Pauli strings x T per R_Z x N_Trotter; T per R_Z from the
                        circuit's own tolerance under ruling R-TOL).
  2033      route (ii)  site-basis LO pionless EFT, qubitized walk + QPE (4He, 12C);
            route (iii) first-quantized QSP (40Ar). Route (ii) uses the lattice pionless
                        Hamiltonian and qubitization walk of arXiv:1911.06368 (lambda, bins and
                        the four registers are DERIVED from its cited t, C0, D0 and walk
                        ancilla count; ruling ch02-route-ii-iii-stated, R11b). Route (iii) takes
                        its register and evolution cost from arXiv:2507.22814 (evolution only):
                        Tab. 5 GQSP, eta=40, 8^3, eps=0.1, 3.38e7 T / 498 qubits at Delta omega =
                        100 MeV, REPRODUCED from its Thm. 5 (3.389e7) and REPRICED at 7 T per
                        Toffoli with R-TOL rotations (5.88e7; 5.90e8 at 10 MeV, our extrapolation)
                        (R12, apply_log/r12_ch02.md).
                        Every composed T number (insertion/amplification schedule) stays STATED:
                        the chapter says it is the authors' composition, not itemized.
  codesign  route (i)   converged-basis 12C, N_orb 150-200; DERIVED in prose with the
                        same per-step formula and N_Trotter = 200.

Revised after the adversarial verification pass (see NEEDS_AUTHOR.md): the 12C
lambda scaling now starts from the stated lambda = 19.1 GeV (not from the 29.3 MeV bin,
which made the W_q = 13 check circular; W_q = 13 is now stated by the chapter and the
bin is checked against it); tau_max is derived from sigma (the chapter's own
"tau_max ~ 1/sigma"); the fault budget is one report-wide input (0.1 expected faults
per shot, ruling R3 2026-09-28) and the old unit-fault figure is kept only for comparison;
the model raises when a stated peak register is smaller than the site register it must
hold; and every chapter number the verifier found missing (kinematic coverage, the
configuration count, the Toffoli-to-T equivalents, the campaign-aggregated ops, the
route-(i) 40Ar top-end ratio, the V^2 scaling checks, the 2028/2033 ancilla ratio, the
Stated-vs-derived diffs) is emitted in `intermediates`.

INPUTS AND SOURCES
  N_orb (2028)            10          app01:125,196  active 0s-0p-1s0d on A=4 (Stated)
  Pauli strings / step    N_orb^4/2   arxiv_2312_05344 via app01:106 (Cited)
  T per R_Z               RUS fit 1.15 log2(1/eps_rot) + 9.2   bocharovRoettelerSvore2015,campbell2017
                                      (the chapter's formula, common.t_per_rotation "rus"; Cited).
  eps_rot                 sqrt(eps_syn / N_rot) per circuit, eps_syn = 1e-2 per shot: ruling R-TOL
                          (H. Lamm, 2026-09-29; common.eps_rot_for). N_rot = synthesized rotations in
                          ONE shot of that circuit: 5e3 (2028 chiral) -> 1.414e-3, 20.086 T;
                          3e3 (pionless, implied 1e3/step x 3) -> 1.826e-3, 19.662 T;
                          5.0625e10-1.6e11 (codesign) -> 4.444e-7-2.5e-7, 33.467-34.421 T.
                          Replaces the fixed eps_gs = 1e-4 and the '~30 T' (RUS 24.5, rounded up).
                          The compiled SU(3) papers use the slope-only 1.15 log2(1/eps) (emitted
                          for comparison only).
  N_Trotter (2028)        1 chiral / 3 pionless     app01:126 (Stated)
  pionless per-step gain  ~5x         arxiv_2312_05344 via app01:126 (Cited)
  2028 ancilla split      80 QPE+window, 50 LCU/BE, 1-2 Hadamard   app01:203-204 (Stated; the box marks the
                          80 and 50 '(budgeted, not derived)', ruling ch02-ancilla-4x (a), R11 2026-09-30)
  route-(ii) walk         arxiv_1911_06368 appendix 'Gate cost of the qubiterate' (Cited, app01:108):
                          Gamma = 38M terms, N_A = ceil(log2 Gamma); per controlled step 2N_A-1 ancillas,
                          15*2^(N_A-1)+14N_A-37 T and 2^N_A-N_A-2 U(2), a lower bound (Clifford controlled
                          unitaries neglected, one prepare per step). One R_y per U(2) (real prepare
                          amplitudes; ruling ch02-no-compiled-walk (i)); 2^W_q - 1 steps per shot.
  eps (accuracy)          5%          app01:57,289 (Stated); N_comp = 5 app01:90 (Stated)
  wall-time convention    t_gate_s = 1 us per T-gate (hard operation), report-wide; ruling R8
                          (H. Lamm, 2026-09-29), stated in Ch. 1 (Assumed). Was 10 us per logical cycle
                          before R8. Plus shot_overhead_s = 0.1 ms per shot (init, readout, decode).
                          R9 (2026-10-02) renamed logical_cycle_s -> t_gate_s; values unchanged.
  one machine             ruling E27 (H. Lamm, 2026-10-02): every wall time is serial on ONE machine.
                          No machine count enters. Co-design: 1.72e4-1.12e5 yr serial (printed
                          '~1.7e4-1.1e5'); 5 years holds 28-93 whole shots; cost reduction to fit = serial
                          years / 5 = 3436-22338 (printed '~3e3-2e4x'), an algorithmic/co-design reduction,
                          not a machine count. The R8 multi-machine fields are retired (RETIRED_E27).
  shot framing (r25)      rulings R1/R2 (H. Lamm, 2026-10-02) and the shot audit (shot_audit.json, Ch. 2):
                          QPE drops each shot into one bin, so a bin of weight p reaches relative error eps
                          after (1-p)/(p eps^2) shots per (shell, component). 2033 4He: 5% only on bins with
                          p >= 0.1 (R2) -> 3600 per (shell, comp); 3 box shells (derived from the 3^3 lattice);
                          first result R_L only at 30% (R1) -> 100 per shell, 300 shots; campaign 5 components
                          -> 5.4e4 shots. Normalization ||J Psi0||^2 is common to every bin: <= 2% for the
                          campaign (combined 5.4%), <= 10% for the first result, from ground-state snapshots
                          N = v/eps^2 with v = 2/3 (audit estimate, 1.7e3 at 2%), one snapshot set per
                          component, each snapshot one ground-state load (<= 1.5% of 6.7e8 T at the block
                          rung, 28.5% of 9.6e8 at the CI ceiling). Replaces the Stated 5.7e4-3.1e5 band
                          (no derivation existed; the audit reconstructs it as 5% on bins down to p ~ 0.02).
                          2028: J|Psi0> is NOT amplified (no room under 1e5 T); the Hadamard test estimates
                          C/alpha^2, so relative eps costs (1-mu^2)/(eps mu)^2 per (comp, Re/Im) with
                          mu = ||J Psi0||^2/alpha^2 = 1.5/10^2 (alpha ~ 10 illustrative, HO-basis alpha not
                          computed; audit). First result Re C(tau) of R_L at one tau, 30%; full benchmark
                          L and T, Re and Im, 5%. Replaces '~2e3 shots, ~3 min' (the +/-1 bound was applied
                          to the unnormalized quantity).
  RFI envelopes           1e5 T / 150-250 LQ (2028); 1e9 T / 1000 LQ (2033)  DOE_RFI_2026 (Cited)
  fault budget            0.1 expected faults per shot, app01:100,137 (Stated; R3 report-wide).
                          The pre-R3 'order unity' figure is Assumed(1.0) for comparison only.
  lambda                  M x 706.76 MeV from arxiv_1911_06368 (Cited t = 10.5794, C0 = -98.2266, D0 = 127.8397 MeV
                          at a = 1.4 fm; lambda/M = 24t + 3|C0+D0| + (3/2)|C0+2D0| + D0): 19.08 / 45.23 / 88.35 GeV
                          on 3^3 / 4^3 / 5^3; printed '19.1' / '45.2' (app01:135,144) checked to 1%.
                          W_q = 12 app01:137 (Stated); W_q = 13 for 12C app01:152,234 (Stated)
  box sides               4.2 / 5.6 / 7.0 fm for 3^3 / 4^3 / 5^3   app01:137,150 (Stated);
                          hbar c = 197.327 MeV fm (PDG)
  N_orb (converged 12C)   150-200     app01:65,166 (Stated); N_orb (16O) 200-300 app01:94
  sigma (DUNE)            10 MeV      app01:90,145 (Stated) -> tau_max = 1/sigma = 100 GeV^-1
  N_Trotter (codesign)    200         app01:166 (Stated; step 0.5 GeV^-1 not derived)
  N_anc                   200         app01:94: 16O 800-1200 system -> 1000-1400 LQ (Stated)
  N_anc bracket           100-300     app01:90 (Stated)
  configurations          5 comp x 5-10 q x 3 nuclei ~ 1e2 app01:279 (Stated)
  tau grid (codesign)     Nyquist step pi/q0_max at |q0| <= 1 GeV (app01:166): 3.14 GeV^-1, ceil(100/3.14) = 32
                          points to tau_max = 100 GeV^-1 (ruling ch02-codesign-shots (c), R11b). N_bin = (5-10) x 32
                          = 160-320; shots 400 x 5 x N_bin = 3.2-6.4e5 (printed '~3-6e5'); per configuration
                          400 x 32 = 12800 (printed '~1.3e4').
  Toffoli -> T            chapter names no convention; both 7 (textbook) and 4 (jones) emitted
WHAT IS NOT DERIVED HERE
  Every route-(ii)/(iii) 2033 T figure: 4He 6.7-9.6e8 T; 12C 2.7-2.8e9 T; 12C moments 'a few
  x 1e7 T' (not separately costed; R11b); 40Ar 500-601 LQ (601 composed); the composed
  factors (0.15-9.3x, <=1.6x), the retained-amplitude fractions (68-70%, <=1.5%, 15.6-28.5%),
  the walk share (~35%), the +1-bit (1.34-1.38x) and +1-cube (5.4-5.6x) factors, the 5^3 box
  (10-11x; both candidate laws fall in it), the tau-grid companion (63-65%), the Euclidean
  variant (1.5-3.1x), and every 2033 eps_l value: all Stated (the 4He shot count is now derived, r25);
  app01:135 says the composition is the authors' own and is not itemized. The route-(ii)
  registers (141/130/292/280 LQ) and lambda are DERIVED from arxiv_1911_06368 (R11b). The walk alone now has a COMPILED lower bound per query from
  the 1911.06368 appendix (4He: 4095 x 6.887e4 = 2.82e8 T, 29-42% of the stated range,
  bracketing the stated ~35%); the rest of the composition is still not derived.
  The model checks what CAN be checked from them (eps_l from the fault budgets, wall time
  from the cycle, campaign years, the 29.3 MeV bin from lambda and W_q, the site-basis
  registers from 4 modes/site, the V^2 scaling the +1-cube factor implies, the kinematic
  coverage from the box side) and records the rest as-is. Because the 2033 T is Stated it
  does not respond to lattice size or W_q; the +1-bit response is emitted separately.
  Also Stated: N_Trotter = 200 (codesign), the route-(i) 40Ar extrapolation.
  R12: the 40Ar evolution (5.9e7-5.9e8 T, derived here from 2507.22814 Thm. 5 at 7 T per Toffoli), the
  8^3 comparison (Tab. 5 of 2507.22814: 10x T, 6.2x qubits at t_dw; 28x T at t_cross) and the
  first-quantized register form A ceil(log2 4 N_orb) are now Cited/derived; '~700 LQ first-quantized'
  and '500x, 2100 vs 500 LQ' are gone (no source supported them).
  No gauge group enters this chapter, so groups.py is not used.
  Not derived (r25): the 2028 current normalization alpha ~ 10 and ||J Psi0||^2 ~ 1.5 (audit, illustrative);
  the snapshot variance v = 2/3; the bin-weight floor p >= 0.1 is a ruling, the actual bin weights are not
  computed (audit estimate: peak 0.2-0.35, ~5 populated 29.3 MeV bins).
WALL TIME AND FACTORIES
  t_gate_s = 1e-6 s per T, shot_overhead_s = 1e-4 s per shot, ONE machine, serial (E27).
  T-depth per shot from factory.json (Ch. 2 runs), low/high:
    2028 (route i, 1.00428e5 T): commuting Pauli-string rotations on disjoint mode quadruples run at the same
      time, 4 N_orb / 4 = 10 per layer, up to 8 x 10 = 80 with the commuting strings of each double excitation;
      D_T = (5e3 / W) x 20.086 = 1.26e3 (W = 80) - 1.00e4 (W = 10). F* = 10-80: the baseline holds.
    2033 4He (route ii): QPE steps are sequential. PREPARE's level-k multiplexor is 2^k commuting rotations,
      spread over the walk's own 2N_A-1 = 21 ancillas: D_prep = sum_k ceil(2^k/21) x 26.24 x 2035/2047.
      SELECT is sequential over 2^11 indices: T-depth 1 per index (measurement-uncomputed AND) to 3 per
      7-T Toffoli (3 x 15477/7). Walk F* = 14.29 / 7.32; the un-itemized remainder is given the walk's
      parallelism. D_T = 4.69e7 (6.7e8 T, depth-1 SELECT) - 1.31e8 (9.6e8 T, depth-3 SELECT).
      F* (cross-paired, an artifact of pairing the depth-1 compile with t_lo and depth-3 with t_hi) 5.1-20.5;
      physical F* per compile 7.32-14.29 (f_star_physical). Rule 12: at the depth-3 end the shot is
      reaction-limited (1311 s, 1.37x), so the box walls run from the depth-1 compile at 6.7e8 T to the
      depth-3 compile at 9.6e8 T, never the 10-factory baseline at the depth-3 end.
  First result (2033): R_L, 3 shells, 30% on bins with p >= 0.1, 300 shots + 67 normalization snapshots
    (<= 10%, combined 31.6%) -> wall_first_result_s = 2.017e5-4.115e5 s = 2.33-4.76 d, 11.2-21.9 min per shot.
    Campaign: 5 components, 3 shells, 5% on p >= 0.1, 5.4e4 shots + normalization snapshots (1667 per
    component, x1-3 for the four grouped-basis components: 8.3e3-2.2e4) -> wall_campaign_s =
    3.626e7-7.672e7 s = 1.149-2.431 yr (fits the 5-year horizon).
  First result (2028): Re C(tau), R_L, one tau, 30%: 4.94e4 shots, 82.7 min. Full benchmark: L+T, Re+Im, 5%:
    7.11e6 shots, 8.27 d (alpha ~ 10, C ~ 1.5, mu = 0.015; shots scale as alpha^4).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, fields

from estimates.common import (Assumed, Cited, Stated, Tagged, Primitive, CircuitStatus,
                              Result, Published, ERAS, SYNTHESIS_SRC, T_PER_TOFFOLI,
                              EPS_SYN, EPS_SYN_SRC, eps_rot_for, t_per_rotation, toffoli_t,
                              depth_exports, REACTION_TIME_S)

TEX = "app01"
YEAR_S = 3600.0 * 24 * 365.25
MONTH_S = YEAR_S / 12
# The chapter's synthesis formula: the full RUS fit 1.15 log2(1/eps) + 9.2 (R-TOL keeps it; only eps moves)
SYNTH_MODEL = "rus"
# E27 (H. Lamm, 2026-10-02): "we don't want to anywhere assume we have multiple machines". The R8 multi-machine
# inputs and the strings they backed are retired. History only; nothing in the model reads this.
RETIRED_E27 = {
    "parallel_machines": (10, 100),
    "campaign_months_on_10": (1.5, 11),          # box '~1.5-11 months on 10 shot-parallel machines'
    "campaign_months_on_100": (0.15, 1.1),       # prose '0.15-1.1 months on 10^2 shot-parallel machines'
    "parallelism_for_horizon": (3e3, 2e4),       # box 'a 5-yr horizon needs ~3e3-2e4x parallelism'; the same ratio is
                                                 # now printed as the cost reduction needed on one machine
}
# r25 (rulings R1/R2, shot audit, 2026-10-02): the Stated 4He campaign band and the statements built on it are
# replaced by the derived first-result / campaign tiers. History only; nothing in the model reads this.
RETIRED_R25 = {
    "shots_4he": (5.7e4, 3.1e5),                 # '5.7e4-3.1e5 shots over the kinematic grid' (no derivation on disk)
    "campaign_years": (1.2, 9.4),                # '1.2-9.4 years on one machine'
    "campaign_over_horizon_4he": 1.9,            # 'upper end x1.9 over 5 years'
    "shots_in_horizon_4he": 1.6e5,               # '5 years on one machine holds 1.6e5 of its 3.1e5 shots'
    "shots_2028": 2e3,                           # '~2e3 (single tau point, full mu-nu tensor)'
    "total_wall_2028_min": 3.35,                 # '~3 min of dedicated time'
}


@dataclass(frozen=True)
class Assumptions:
    # ---- route (i): shared cost model -----------------------------------------------
    pauli_fraction: Tagged = Cited(0.5, "arxiv_2312_05344",
                                   "app01:106 'N_orb^4/2 Pauli strings per step, the factor 1/2 "
                                   "being the index symmetry of the two-body term'")
    # R-TOL (2026-09-29): total synthesis error per shot, report-wide; each circuit sets
    # eps_rot = sqrt(eps_syn / N_rot). Replaces eps_gs = 1e-4 and the fixed '~30 T per R_Z'.
    eps_syn: Tagged = Assumed(EPS_SYN, "ruling R-TOL: total synthesis error per shot 1e-2, report-wide; "
                                       "eps_rot = sqrt(eps_syn/N_rot) per circuit (randomized synthesis)",
                              EPS_SYN_SRC)
    # ---- 2028 --------------------------------------------------------------------------
    n_orb_2028: Tagged = Stated(10, f"{TEX}:125,196", "'N_orb=10 (active 0s-0p-1s0d basis on A=4)'")
    n_trotter_2028: Tagged = Stated(1, f"{TEX}:126", "'depth-limited to N_Trotter=1 low-order step'")
    pionless_gain: Tagged = Cited(5, "arxiv_2312_05344",
                                  "app01:126 'per-step T-count drops by ~5x' in pionless EFT")
    n_trotter_pionless: Tagged = Stated(3, f"{TEX}:126,207", "'we run N_Trotter=3'")
    # ruling ch02-ancilla-4x (a), R11: keep ~170 (T1); the box marks these two rows '(budgeted, not derived)'
    anc_qpe_window: Tagged = Stated(80, f"{TEX}:203", "'~80 QPE energy register + window ancilla (budgeted, not derived)'")
    anc_lcu_be: Tagged = Stated(50, f"{TEX}:204", "'~50 LCU/BE dirty (budgeted, not derived)'")
    anc_hadamard: Tagged = Stated((1, 2), f"{TEX}:204", "'1--2 Hadamard'")
    lq_2028_prose: Tagged = Stated(170, f"{TEX}:94", "'4He on route (i) fits in ~170 LQ (40 system + ~130 ancilla)'; box says '~170'")
    eps_target: Tagged = Stated(0.05, f"{TEX}:57,289", "'5% relative accuracy'")
    n_comp: Tagged = Stated(5, f"{TEX}:90", "'five for the inclusive weak response of an unpolarized target'")
    n_bin_2028: Tagged = Stated(1, f"{TEX}:208", "'single tau point, full mu-nu tensor'")
    t_gate_s: Tagged = Assumed(1e-6, "1 us per T-gate (hard operation), the report-wide wall-time "
                                     "convention of ruling R8 (H. Lamm, 2026-09-29; convention in "
                                     "Ch. 1); was 10e-6. R9 (2026-10-02) renamed it from logical_cycle_s. "
                                     "app01 2028 box '~0.10 s'; 2033 box '11-16 min'",
                               f"{TEX}:2028,2033 boxes")
    shot_overhead_s: Tagged = Assumed(1e-4, "per-shot overhead t0 ~ 0.1 ms (register initialization, final readout, "
                                            "decode), the report rule stated in Ch. 1; as Ch. 8 shot_overhead_s. "
                                            "0.1% of the 2028 shot, < 2e-7 of a 2033 shot: no printed number moves (E27)",
                                      "E27; ch08_baryogenesis.py")
    # r25 (R1, shot audit): the 2028 Hadamard test is unamplified; it estimates C/alpha^2
    alpha_2028: Tagged = Assumed(10, "block-encoding normalization of the current in the N_orb=10 oscillator basis, "
                                     "illustrative (shot audit; the HO-basis alpha is not computed; the lattice "
                                     "rho(q) has alpha ~ 27)", "shot_audit Ch. 2")
    jpsi_norm_sq_2028: Tagged = Assumed(1.5, "||J Psi0||^2 ~ Z S_L(q) scale used by the shot audit with alpha ~ 10 "
                                             "(mu = 0.015); illustrative", "shot_audit Ch. 2")
    n_comp_2028_full: Tagged = Assumed(2, "the (e,e') validation uses the L and T responses only (shot audit: "
                                          "'generous on N_comp')", "shot_audit Ch. 2")
    n_parts_2028: Tagged = Assumed(2, "Re and Im of C(tau) are separate Hadamard circuits", "shot_audit Ch. 2")
    eps_first: Tagged = Assumed(0.3, "ruling R1 (H. Lamm, 2026-10-02): minimal first result at 30% statistical error",
                                "R1")
    # r25 2028 T-depth (factory.json Ch. 2): rotations per layer, 4 N_orb/4 disjoint quadruples .. x8 commuting
    # strings of one double excitation (each concurrent rotation holds one parity ancilla; the 2028 machine has room)
    strings_per_excitation_2028: Tagged = Assumed(8, "the 8 Pauli strings of one double excitation commute pairwise "
                                                     "(factory.json Ch. 2, 2028 run)", "factory.json")
    sigma_2028_MeV: Tagged = Stated(50, f"{TEX}:128,197", "'coarse Gaussian window (sigma ~ 50 MeV)'")
    # ---- RFI envelopes --------------------------------------------------------------------
    rfi_2028_t: Tagged = Cited(1e5, "DOE_RFI_2026", "app01:126,192 '10^5 hard ops'")
    rfi_2028_lq: Tagged = Cited((150, 250), "DOE_RFI_2026", "app01:192")
    rfi_2033_t: Tagged = Cited(1e9, "DOE_RFI_2026", "app01:166,216 '~10^9 hard ops'")
    rfi_2033_lq: Tagged = Cited(1000, "DOE_RFI_2026", "app01:216 '>~ 10^3 LQ'")
    eps_l_2028_floor: Tagged = Cited(1e-8, "DOE_RFI_2026", "app01:185 '2028 RFI floor eps_l <= 1e-8'")
    eps_l_2028: Tagged = Stated(1e-6, f"{TEX}:2028 box", "'<~ 1e-6 (<= 0.1 expected fault at 1.004e5 T)'; "
                                "0.1/1.00428e5 = 9.96e-7 (R-TOL; was 7e-7 at 1.5e5 T)")
    # ---- fault budget: one report-wide convention (ruling R3, 2026-09-28) -----------
    fault_budget: Tagged = Stated(0.1, f"{TEX}:100,137",
                                  "'expected number of logical faults per circuit to <~ 0.1'; '<= 0.1-expected-fault budget per shot'")
    fault_budget_unit: Tagged = Assumed(1.0, "the pre-R3 'order unity' convention app01:100 used to print; no longer "
                                             "in the chapter, kept only so the comparison values can be emitted", "R3")
    # ---- 2033 primary: 4He, route (ii) ------------------------------------------------
    t_4he: Tagged = Stated((6.7e8, 9.6e8), f"{TEX}:137,224", "'composed per-shot cost ... 6.7-9.6e8 T'")
    lq_4he: Tagged = Cited(141, "arxiv_1911_06368", "app01:137,227 'peak register is 141 LQ'; derived 108 site + "
                           "W_q 12 + 2N_A-1 21 (R11b retag; the model raises if the components disagree)")
    lattice_4he: Tagged = Stated(3, f"{TEX}:137,221", "'3^3 lattice (box side 4.2 fm)'")
    box_4he_fm: Tagged = Stated(4.2, f"{TEX}:137", "'box side 4.2 fm'")
    hbar_c_MeV_fm: Tagged = Cited(197.327, "PDG", "hbar c, for |q| = 2 pi hbar c / L")
    q_coverage_3cube_GeV: Tagged = Stated((0.30, 0.51), f"{TEX}:246", "'|q| ~ 0.30-0.51 GeV (3^3)'")
    q_coverage_4cube_GeV: Tagged = Stated((0.22, 0.77), f"{TEX}:246", "'0.22-0.77 GeV (4^3)'")
    w_q: Tagged = Stated(12, f"{TEX}:137,222", "'W_q = 12'")
    # R11b ch02-route-ii-iii-stated: lambda derived from the 1911.06368 lattice pionless couplings (Table tab:hparams,
    # a = 1.4 fm); C0 is attractive (negative), stored as its magnitude because inputs must be positive
    hop_t_MeV: Tagged = Cited(10.5794, "arxiv_1911_06368", "Table tab:hparams: t at a = 1.4 fm (from Rokash2013)")
    c0_abs_MeV: Tagged = Cited(98.2265511, "arxiv_1911_06368", "Table tab:hparams: C0 = -98.2265511 MeV (sign applied in _lambda_per_site)")
    d0_MeV: Tagged = Cited(127.839693, "arxiv_1911_06368", "Table tab:hparams: D0 = 127.839693 MeV")
    lambda_GeV: Tagged = Cited(19.1, "arxiv_1911_06368", "app01:135,144 'lambda ~ 19.1 GeV on the 3^3 lattice'; derived "
                               "27 x 706.76 MeV = 19.083 GeV from the cited t, C0, D0 (R11b retag; checked to 1%)")
    bin_4he_MeV: Tagged = Stated(29.3, f"{TEX}:144,223", "'physical bin of 29.3 MeV'")
    sigma_dune_MeV: Tagged = Stated(10, f"{TEX}:90,145", "'sigma ~ 10 MeV' DUNE-relevant window; tau_max ~ 1/sigma. "
                                    "Compared with the bin at equal FWHM (ruling ch02-sigma-window-vs-register (a))")
    bin_target_MeV: Tagged = Stated(5, f"{TEX}:145", "'A 5 MeV-binned 4He response ... three more register bits'")
    eps_l_4he: Tagged = Stated(1e-10, f"{TEX}:137,226", "'eps_l ~ 1e-10, floor-class'")
    block_frac: Tagged = Stated((0.68, 0.70), f"{TEX}:137", "'68-70% at the zero-momentum-block rung'")
    gs_prep_frac_block: Tagged = Stated(0.015, f"{TEX}:139", "'<= 1.5% of the composed cost at and below the block rung'")
    gs_prep_frac_ci: Tagged = Stated((0.156, 0.285), f"{TEX}:139", "'15.6-28.5% at the unrestricted-CI ceiling'")
    walk_frac: Tagged = Stated(0.35, f"{TEX}:141", "'qubitization walk accounts for only approximately 35%'")
    # ruling ch02-no-compiled-walk (i), R11: the walk priced per query from the 1911.06368 appendix
    walk_terms_per_site: Tagged = Cited(38, "arxiv_1911_06368", "appendix 'Gate cost of the qubiterate': "
                                        "Gamma_K = 24M + Gamma_V = 14M = 38M for N_f = 4; app01:108,137")
    walk_rot_per_u2: Tagged = Assumed(1, "one R_y per prepare-tree node: the amplitudes sqrt(lambda_j) are real "
                                         "(ruling ch02-no-compiled-walk (i)); the paper counts 3 R_z per U(2) "
                                         "'to be conservative'", "arxiv_1911_06368")
    # r25 (R1, R2, shot audit): the 4He shot count is derived from the QPE bin statistics (was Stated 5.7e4-3.1e5)
    p_bin_min: Tagged = Assumed(0.1, "ruling R2 (H. Lamm, 2026-10-02): 5% only on bins holding >= 10% of the "
                                     "strength", "R2")
    n_comp_first: Tagged = Assumed(1, "ruling R1 / shot audit: the first result is R_L alone (rho(q) insertion)", "R1")
    eps_norm_campaign: Tagged = Assumed(0.02, "normalization ||J Psi0||^2 is common to every bin, so <= 2% "
                                              "(combined with 5%: 5.4%); shot audit", "shot_audit Ch. 2")
    eps_norm_first: Tagged = Assumed(0.1, "first-result normalization <= 10% (combined with 30%: 31.6%)", "r25")
    norm_grouped_factor: Tagged = Assumed((1, 3), "snapshot multiplier for the four non-L components: only "
                                                  "rho^dag rho is Z-diagonal; R_T and the axial sum rules need "
                                                  "grouped non-diagonal bases, 'a few times' R_L (shot audit)",
                                          "shot_audit Ch. 2")
    norm_snapshot_var: Tagged = Assumed(2 / 3, "relative variance per ground-state snapshot of the sum-rule "
                                               "estimator, N = v/eps^2 (shot audit: 1.7e3 at 2%)", "shot_audit Ch. 2")
    # r25 2033 T-depth (factory.json Ch. 2, 2033 run)
    prepare_ancilla_4he: Tagged = Assumed(21, "PREPARE's commuting multiplexor rotations spread over the walk's own "
                                              "2N_A-1 = 21 ancillas, idle during PREPARE (factory.json)", "factory.json")
    select_depth_per_index: Tagged = Assumed((1, 3), "SELECT unary iteration over 2^N_A indices: T-depth 1 per index "
                                                     "(measurement-uncomputed AND) to 3 per 7-T Toffoli "
                                                     "(factory.json bracket)", "factory.json")
    plus_bit_factor: Tagged = Stated((1.34, 1.38), f"{TEX}:141", "'Adding one register bit ... factor of 1.34-1.38'")
    block_plus_bit_frac: Tagged = Stated((0.91, 0.97), f"{TEX}:141", "'91-97% of the depth limit'")
    plus_cube_factor: Tagged = Stated((5.4, 5.6), f"{TEX}:141", "'one cube at fixed resolution ... 5.4-5.6'")
    tau_grid_frac: Tagged = Stated((0.63, 0.65), f"{TEX}:142", "'63-65% of the envelope ... on 130 LQ'")
    tau_grid_lq: Tagged = Cited(130, "arxiv_1911_06368", "app01:142 'on 130 LQ'; derived 108 + 1 Hadamard + 21 (R11b retag)")
    # E27 (H. Lamm, 2026-10-02): one machine only. The R8 fields parallel_machines, campaign_months and
    # campaign_months_100 are retired (values kept in RETIRED_E27); the campaign is serial on one machine.
    modes_per_site: Tagged = Assumed(4, "p/n x spin up/down per site; app01:148 'occupation of a fixed set "
                                        "of proton and neutron modes at each lattice site'; 4x27=108 and "
                                        "4x64=256 reproduce 141 and 292 with ~30 ancilla", f"{TEX}:148")
    spin_isospin_bits: Tagged = Cited(2, "arxiv_2507_22814",
                                      "Sec. 3 register Q = (dm+2) eta [:180]: 2 spin-isospin qubits per nucleon. Its Thm. 5 "
                                      "qubit formula (Eq. 132) counts only n_s = eta d m, since the SU(4) H never acts on them; "
                                      "a spin/isospin-changing current insertion does (R12 verifier)")
    # ---- 2033 near miss and on-ramp: 12C, route (ii) --------------------------------
    lattice_12c: Tagged = Stated(4, f"{TEX}:150,233", "'smallest realizable cube that holds it is 4^3'")
    box_12c_fm: Tagged = Stated(5.6, f"{TEX}:150", "'4^3, with box side 5.6 fm'")
    box_5cube_fm: Tagged = Stated(7.0, f"{TEX}:150", "'A 5^3 box, 7.0 fm'")
    t_12c: Tagged = Stated((2.7e9, 2.8e9), f"{TEX}:152,234", "'2.7-2.8e9 T at 34.7 MeV physical bins'")
    lq_12c: Tagged = Cited(292, "arxiv_1911_06368", "app01:150,235 'from 141 to 292 LQ'; derived 256 + W_q 13 + 23 (R11b retag)")
    bin_12c_MeV: Tagged = Stated(34.7, f"{TEX}:152,234", "'34.7 MeV physical bins'")
    w_q_12c: Tagged = Stated(13, f"{TEX}:152,234", "'(W_q = 13; ...)'")
    eps_l_12c: Tagged = Stated(3.6e-11, f"{TEX}:152,235", "'eps_l ~ 3.6e-11' (0.1 fault at 2.8e9 T)")
    # R11b ch02-5cube-factor fallback: 'x10-11 ... at comparable bins (W_q=14, 33.9 MeV)'; both candidate laws fall
    # in it (V^2: 10.30-10.68; one W_q bit x one LCU-index bit: 10.78-11.18); the author did not say which produced
    # the old 10.8-11.1
    five_cube_factor: Tagged = Stated((10, 11), f"{TEX}:150", "'5^3 box ... x10-11 the depth envelope at comparable bins "
                                      "(W_q=14, 33.9 MeV)'")
    w_q_5cube: Tagged = Stated(14, f"{TEX}:150", "'(W_q=14, 33.9 MeV)'")
    te_pe_gain: Tagged = Cited(8, "arxiv_1911_06368", "app01:152 'up to ~8x cheaper'; different depth convention, not folded in")
    # R11b ch02-3p4e7-coincidence fallback: the on-ramp was never separately costed; 'a few x 1e7 T'
    t_onramp: Tagged = Stated((2e7, 5e7), f"{TEX}:154,239", "'a few x 1e7 T (ground-state load at the block rung plus "
                              "one insertion and a short moment readout; not separately costed)'")
    lq_onramp: Tagged = Cited(280, "arxiv_1911_06368", "app01:154,239 '280 LQ'; derived 256 + 1 Hadamard + 23 (R11b retag)")
    t_onramp_eps_ref: Tagged = Stated(3e7, f"{TEX}:154,240,288", "'eps_l <~ 3e-9 at 3e7 T' (reference depth for the printed eps_l)")
    eps_l_onramp: Tagged = Stated(3e-9, f"{TEX}:154,240,288", "'eps_l <~ 3e-9 at 3e7 T' (0.1 fault / 3e7 = 3.33e-9, rounded once)")
    euclid_factor: Tagged = Stated((1.5, 3.1), f"{TEX}:154", "'Euclidean tau-grid variant ... x1.5-3.1'")
    # ---- 2033 stretch: 40Ar, route (iii) --------------------------------------------
    A_ar: Tagged = Stated(40, f"{TEX}:44", "'A = 40, Z = 18 for argon'")
    # R12 (2507.22814 on disk): route (iii) is the GQSP algorithm of arXiv:2507.22814 (Spagnoli, Lissoni,
    # Roggero; v3, Quantum 2026). Equation numbers are the published PDF's (main_final.tex lines in notes).
    lattice_ar: Tagged = Cited(8, "arxiv_2507_22814", "Tab. 5 caption 'All estimates use a 8x8x8 lattice' (m = 3); "
                               "main_final.tex:1422")
    fq_dim: Tagged = Cited(3, "arxiv_2507_22814", "d = 3 spatial dimensions (Sec. 4, main_final.tex:1362)")
    fq_eps: Tagged = Cited(0.1, "arxiv_2507_22814", "Tab. 5 caption 'target error eps = 0.1' (main_final.tex:1422)")
    fq_kin_MeV: Tagged = Cited(10.58, "arxiv_2507_22814", "Tab. 4 hbar^2/(2 mu a^2) = 10.58 MeV at a = 1.4 fm (main_final.tex:1356)")
    fq_c_abs_MeV: Tagged = Cited(98.23, "arxiv_2507_22814", "Tab. 4 C = -98.23 MeV (sign applied in _gqsp_lambda)")
    fq_g_MeV: Tagged = Cited(127.84, "arxiv_2507_22814", "Tab. 4 G = 127.84 MeV")
    fq_spacing_fm: Tagged = Cited(1.4, "arxiv_2507_22814", "Sec. 4 'lattice spacing a = 1.4 fm' (main_final.tex:1332)")
    fq_mass_MeV: Tagged = Cited(939, "arxiv_2507_22814", "Sec. 4 'nucleon mass mu = 939 MeV', for t_cross Eq. (137)")
    fq_cross_E_MeV: Tagged = Cited(10, "arxiv_2507_22814", "Eq. (137): crossing time of an E = 10 MeV nucleon")
    fq_heavy_H_per_nucleon_MeV: Tagged = Cited(18, "arxiv_2507_22814", "Eq. (138): Delta H = ||T|| + ||V|| + 18 eta "
                                               "(main_final.tex:1346)")
    fq_dw_MeV: Tagged = Cited(100, "arxiv_2507_22814", "Eq. (138) text 'we consider Delta omega = 100 MeV' (main_final.tex:1346)")
    fq_dw_extrap_MeV: Tagged = Stated(10, f"{TEX}:156,249", "'Delta omega ~ 10 MeV': our extrapolation with the source's "
                                      "formulas; the source costs only 100 MeV")
    fq_t_src_per_toffoli: Tagged = Cited(4, "arxiv_2507_22814", "Thm. A.4 proof and :443 'implement each Toffoli with 4 T "
                                         "gates' (Gidney temporary AND); repriced at T_PER_TOFFOLI['textbook'] = 7 (R5)")
    # Tab. 5 (main_final.tex:1404-1425), eta = 40 unless named; T-counts in the source's 4-T convention
    tab5_gqsp_t_dw: Tagged = Cited(3.38e7, "arxiv_2507_22814", "Tab. 5 GQSP eta=40 T-count [t_dw]")
    tab5_gqsp_q_dw: Tagged = Cited(498, "arxiv_2507_22814", "Tab. 5 GQSP eta=40 qubits [t_dw]")
    tab5_gqsp_t_cross: Tagged = Cited(2.11e8, "arxiv_2507_22814", "Tab. 5 GQSP eta=40 T-count [t_cross]")
    tab5_gqsp_q_cross: Tagged = Cited(500, "arxiv_2507_22814", "Tab. 5 GQSP eta=40 qubits [t_cross]")
    tab5_gqsp16_t_dw: Tagged = Cited(5.99e6, "arxiv_2507_22814", "Tab. 5 GQSP eta=16 T-count [t_dw]")
    tab5_gqsp16_t_cross: Tagged = Cited(3.75e7, "arxiv_2507_22814", "Tab. 5 GQSP eta=16 T-count [t_cross]")
    tab5_gqsp16_q: Tagged = Cited((234, 235), "arxiv_2507_22814", "Tab. 5 GQSP eta=16 qubits [t_dw], [t_cross]")
    tab5_trotter2q_t_dw: Tagged = Cited(3.47e8, "arxiv_2507_22814", "Tab. 5 second-order Trotter, 2nd quantization "
                                        "[97] = arxiv_2312_05344, eta=40, T-count [t_dw]")
    tab5_trotter2q_t_cross: Tagged = Cited(5.89e9, "arxiv_2507_22814", "Tab. 5 same row, T-count [t_cross]")
    tab5_trotter2q_q: Tagged = Cited(3072, "arxiv_2507_22814", "Tab. 5 same row, qubits (both times)")
    lq_ar: Tagged = Stated((580, 601), f"{TEX}:156,159,243", "'~580' (derived 578 / 580: Tab. 5 498 / 500 + 2 eta "
                                                            "spin-isospin, R12) ... 'worst-case 601-LQ' (Stated)")
    lq_ar_requirements: Tagged = Stated((580, 601), f"{TEX}:286", "'~580-601 LQ (40Ar, conditional)' in the requirements table")
    ar_amplified_factor: Tagged = Stated((0.15, 9.3), f"{TEX}:158,244", "'0.15-9.3x the depth envelope across the amplified bracket'")
    ar_postselected_factor: Tagged = Stated(1.6, f"{TEX}:159,244", "'<= 1.6x throughout'")
    ar_repeats: Tagged = Stated(1e3, f"{TEX}:159", "'up to 10^3 repeats on the campaign axis'")
    eps_l_ar: Tagged = Stated((1.1e-11, 6.7e-10), f"{TEX}:159,245,288", "'from 6.7e-10 (the 0.15x floor) to 1.1e-11 (the 9.3x ceiling)'")
    eps_l_ar_postselected: Tagged = Stated(6.3e-11, f"{TEX}:159", "'6.3e-11 post-selected' (0.1 fault at 1.6e9 T)")
    ar_window_MeV: Tagged = Stated((10, 100), f"{TEX}:156,247", "'Delta omega = 10-100 MeV'")
    # ---- co-design: converged 12C, route (i) ----------------------------------------
    n_orb_explore: Tagged = Stated(20, f"{TEX}:166", "'exploratory N_orb=20 active space'")
    n_orb_conv: Tagged = Stated((150, 200), f"{TEX}:65,166", "'N_orb ~ 150-200 needed for a converged A=12 response'")
    n_orb_16o: Tagged = Stated((200, 300), f"{TEX}:94", "'16O needs N_orb ~ 200-300'")
    n_orb_ar_route_i: Tagged = Stated((250, 500), f"{TEX}:96,119", "'A=40 implies N_orb ~ 250-500'")
    t_ar_route_i: Tagged = Stated((1e13, 5e14), f"{TEX}:119", "'~1e13-5e14 T-gates per shot' (R-TOL; was 1e13-4e14 "
                                  "at 30 T per R_Z); derived 1.37e13-4.75e14")
    tau_max_GeVinv: Tagged = Stated(100, f"{TEX}:90,166", "'tau_max ~ 1/sigma ~ 100 GeV^-1'; model derives it from sigma_dune_MeV")
    n_trotter_codesign: Tagged = Stated(200, f"{TEX}:166", "'in N_Trotter ~ 200 steps'; step size not derived")
    window_trunc_factor: Tagged = Stated(2, f"{TEX}:90,168", "'tau_max ~ 2/sigma ... doubling the Trotter depth'")
    n_anc: Tagged = Stated(200, f"{TEX}:94", "16O: '800-1200 system qubits, ~1000-1400 LQ' implies N_anc = 200")
    n_anc_bracket: Tagged = Stated((100, 300), f"{TEX}:90", "'ancilla budget N_anc ~ 100-300'")
    A_c: Tagged = Stated(12, f"{TEX}:58", "carbon-12")
    # R11b ch02-codesign-shots (c): Nyquist tau grid; the model derives 3.2-6.4e5 and checks the printed '~3-6e5'
    shots_codesign: Tagged = Stated((3e5, 6e5), f"{TEX}:172,253", "'~3-6e5 shots (32-point tau grid)'")
    q0_max_GeV: Tagged = Stated(1.0, f"{TEX}:166,172", "'|q_0| <= 1 GeV' sets the Nyquist tau step pi/q0_max")
    horizon_yr: Tagged = Stated(5, f"{TEX}:139,185,273", "'5-year campaign horizon'")
    # E27: single-machine statements for the co-design campaign (derived in _model_codesign, checked in the test)
    codesign_serial_yr: Tagged = Stated((1.7e4, 1.1e5), f"{TEX}:185", "'~1.7e4-1.1e5 years' serial on one machine "
                                        "(3.2e5 x 19.6 d - 6.4e5 x 63.7 d = 1.718e4-1.117e5 yr)")
    codesign_shots_in_horizon: Tagged = Stated((28, 93), f"{TEX}:185,249", "'one machine completes 28-93 shots in "
                                               "5 years' (whole shots: floor of 28.65 and 93.13)")
    codesign_reduction: Tagged = Stated((3e3, 2e4), f"{TEX}:185,249", "'~3e3-2e4x' cost reduction to fit 5 years on "
                                        "one machine (serial years / 5 = 3436-22338); an algorithmic/co-design "
                                        "reduction, not a machine count")
    qubitized_toffoli: Tagged = Cited((4e14, 7e15), "arxiv_2607_21563", "app01:170 'e_max=8 (660 orbitals) ... at 0.01 MeV'")
    qubitized_precision_MeV: Tagged = Cited(0.01, "arxiv_2607_21563", "app01:170")
    rescale_bin_MeV: Tagged = Stated(30, f"{TEX}:171", "'rescaled to the ~30 MeV bin width'")
    eps_l_codesign: Tagged = Stated((2e-14, 6e-14), f"{TEX}:185",
                                    "'<= 0.1 expected fault at 1.7-5.5e12 hard ops gives eps_l ~ 2-6e-14' (R3, R-TOL; was "
                                    "'2-7e-14' at 1.5-4.8e12 T); 0.1/5.507e12 = 1.82e-14, 0.1/1.694e12 = 5.90e-14")
    eps_l_codesign_shorthand: Tagged = Stated(1e-13, f"{TEX}:185", "'~1e12 hard ops at eps_l ~ 1e-13/op'; 0.1/1e12, the "
                                              "same shorthand as line 100's 1e-10 at 1e9 ops")
    # ---- requirements summary (app01:277-278) ----------------------------------------
    configs_momentum_transfers: Tagged = Stated((5, 10), f"{TEX}:277", "'5-10 momentum transfers'")
    n_nuclei: Tagged = Stated(3, f"{TEX}:277", "'4He/12C/40Ar'")
    configs_stated: Tagged = Stated(1e2, f"{TEX}:277", "'Number of distinct circuit configurations ~1e2'")
    shots_per_config_stated: Tagged = Stated(1.3e4, f"{TEX}:280", "'Shots per configuration ~1.3e4 on the tau-grid routes (32 tau points "
                                             "x eps^-2 = 400 at eps = 5%, Eq. Nshot); the 2033 primary reads energy "
                                             "bins and needs fewer (next row)' (R11b, was '~1e3')")

    def __post_init__(self):
        for f in fields(self):
            v = getattr(self, f.name)
            if not isinstance(v, Tagged):
                raise TypeError(f"{f.name} must be Tagged")
            if v.is_range and v.lo > v.hi:
                raise ValueError(f"{f.name}: lo > hi")
            if v.lo <= 0:
                raise ValueError(f"{f.name}: must be positive")
        if not (0 < self.eps_target.lo < 1):
            raise ValueError("eps_target must be a fraction")
        if self.pauli_fraction.lo > 1:
            raise ValueError("pauli_fraction is a fraction of N_orb^4")
        for name in ("block_frac", "gs_prep_frac_block", "gs_prep_frac_ci", "walk_frac", "tau_grid_frac"):
            if getattr(self, name).hi > 1:
                raise ValueError(f"{name} is a fraction")
        if self.gs_prep_frac_block.lo + self.walk_frac.lo > 1:
            raise ValueError("gs_prep + walk shares exceed the composed cost")


# --------------------------------------------------------------------------- #
# Route (i) cost model, shared by 2028 and codesign
# --------------------------------------------------------------------------- #

def _route_i(a: Assumptions, n_orb: float, n_trotter: float, gain: float = 1.0) -> dict:
    """One route-(i) circuit under R-TOL: N_orb^4/2 Pauli strings per step (divided by a per-step
    gain for the pionless variant), N_rot = strings x N_Trotter rotations per shot, eps_rot from
    that N_rot, T per rotation from the chapter's RUS fit. Returns every quantity by name."""
    strings = a.pauli_fraction.lo * n_orb ** 4 / gain
    n_rot = strings * n_trotter
    eps = eps_rot_for(n_rot, a.eps_syn.lo)
    t_rot = t_per_rotation(eps, SYNTH_MODEL)
    return {"strings": strings, "n_rot": n_rot, "eps_rot": eps, "t_rot": t_rot,
            "t_step": strings * t_rot, "t_shot": n_rot * t_rot}


def _shots(a: Assumptions, n_bin: float) -> float:
    """Eq. (Nshot): eps^-2 x N_bin x N_comp."""
    return a.eps_target.lo ** -2 * n_bin * a.n_comp.lo


def _site_register(a: Assumptions, L: int) -> int:
    return a.modes_per_site.lo * L ** 3


def _fq_register(A: int, M: int) -> int:
    return A * math.ceil(math.log2(M))


# --------------------------------------------------------------------------- #
# Route (iii): GQSP evolution of arXiv:2507.22814 (R12). Equation numbers are the
# published PDF's; main_final.tex line numbers in brackets.
# --------------------------------------------------------------------------- #

def _gqsp_lambda(a: Assumptions, eta: int, m: int) -> tuple[float, float]:
    """Eq. (90) [:1080]: lambda_H = eta((3|C|+4G)/2 + d K 2^{2(m-1)}), K = (hbar^2/2 mu a^2)(2 pi/2^m)^2 (Eq. kkin)."""
    d = a.fq_dim.lo
    K = a.fq_kin_MeV.lo * (2 * math.pi / 2 ** m) ** 2
    return eta * ((3 * a.fq_c_abs_MeV.lo + 4 * a.fq_g_MeV.lo) / 2 + d * K * 2 ** (2 * (m - 1))), K


def _gqsp_norm_v(a: Assumptions, eta: int) -> float:
    """Eqs. (203)-(204) [:2585-2603]: ||V2+V3|| <= max(a0, a1, a2) with C0 = C < 0, G > 0."""
    c, g = -a.fq_c_abs_MeV.lo, a.fq_g_MeV.lo
    a0 = abs(c) * (eta // 2)
    a1 = abs(3 * c + g) * (eta // 3) + (abs(c) if eta % 3 == 2 else 0)
    r = eta % 4
    a2 = abs(6 * c + 4 * g) * (eta // 4) + (abs(c) if r == 2 else max(abs(c), abs(3 * c + g)) if r == 3 else 0)
    return max(a0, a1, a2)


def _gqsp_t_dw(a: Assumptions, eta: int, m: int, dw_MeV: float) -> dict:
    """Eq. (138) [:1342-1346]: t = (ceil(Delta H / Delta omega) - 1) 2 pi / Delta H, with
    Delta H = ||T|| + ||V|| + 18 eta; ||T|| = d K eta 2^{2m-2} (Eq. 194 [:2524])."""
    _, K = _gqsp_lambda(a, eta, m)
    norm_t = a.fq_dim.lo * K * eta * 2 ** (2 * m - 2)
    norm_v = _gqsp_norm_v(a, eta)
    dh = norm_t + norm_v + a.fq_heavy_H_per_nucleon_MeV.lo * eta
    n = math.ceil(dh / dw_MeV) - 1
    return {"t": n * 2 * math.pi / dh, "delta_h": dh, "norm_t": norm_t, "norm_v": norm_v, "n_steps": n}


def _gqsp_t_cross(a: Assumptions, m: int) -> float:
    """Eq. (137) [:1335]: t_cross = (a L / hbar c) sqrt(mu / 2E), L = 2^m."""
    return (a.fq_spacing_fm.lo * 2 ** m / a.hbar_c_MeV_fm.lo
            * math.sqrt(a.fq_mass_MeV.lo / (2 * a.fq_cross_E_MeV.lo)))


def _t_rot_2507(eps: float) -> float:
    """Thm. A.3, Eq. (144) [:2180]: the source's rotation synthesis, 0.57 log2(1/eps) + 8.83."""
    return 0.57 * math.log2(1 / eps) + 8.83


def _t_qft_2507(n: int, eps: float) -> float:
    """Thm. A.2, Eq. (143) [:2176]: 7N - 11 + sum_{n=3}^{N-1} (8 min(ceil(log2(N/eps)), n) - 15)."""
    return 7 * n - 11 + sum(8 * min(math.ceil(math.log2(n / eps)), k) - 15 for k in range(3, n))


def _gqsp_2507(a: Assumptions, eta: int, m: int, t: float, eps: float) -> dict:
    """Thm. 5 (G-QSP for pionless EFT), Eqs. (130)-(132) [:1249-1270], with T^P, T^S from Thm. 3 Eq. (89)
    [:948-950] in the source's own conventions (4 T per Toffoli, T_ROT of Eq. 144):
      T = Q T^S(e') + 2(Q+1) T^P(e') + 3(Q+1) T_ROT(eps/(6(Q+1))) + 2Q T_MCX(eta+2m+10),
      Q = ceil(2 lambda_H t + 3 ln(6/(eps/2))), e' = eps/(4 lambda_H t),
      T^P = 12 n_eta + 16 b_r + 8m - 64 + 4 T_ROT(e'/32),   b_r = ceil(log2(18 pi^2/e')/2),
      T^S = 24(eta-1)dm + 4(d+2)(m-1) + 16 + 2d T_QFT(m, e'/4d) - 8 eta + 52,
      T_MCX(N) = 4(N-1)  (Thm. A.4).
    Qubits (Eq. 132): eta d m + eta + (2+d)m + 13 + b_r + b_QFT + max(2 b_QFT - 1, dm + 6, eta + 2m + 9),
      b_r = ceil(log2(72 pi^2 lambda_H t / eps)/2), b_QFT = min(m-1, ceil(log2(16 d lambda_H t m/eps))) + 1.
    The polynomial-order log is natural (it reproduces all eight Tab. 5 GQSP entries; log2 misses eta=16 at t_dw
    by 0.5%)."""
    d = a.fq_dim.lo
    lam, _ = _gqsp_lambda(a, eta, m)
    q = math.ceil(2 * lam * t + 3 * math.log(6 / (eps / 2)))
    e1 = eps / (4 * lam * t)
    n_eta = math.ceil(math.log2(eta))
    b_r = math.ceil(0.5 * math.log2(18 * math.pi ** 2 / e1))
    qft = _t_qft_2507(m, e1 / (4 * d))
    tof_p = 3 * n_eta + 4 * b_r + 2 * m - 16           # Lemma 4 proof: PREPARE Toffolis [:684]
    t_p = 4 * tof_p + 4 * _t_rot_2507(e1 / 32)          # = 12 n_eta + 16 b_r + 8m - 64 + 4 T_ROT
    t_s = 24 * (eta - 1) * d * m + 4 * (d + 2) * (m - 1) + 16 + 2 * d * qft - 8 * eta + 52
    tof_s = (t_s - 2 * d * qft) / 4                     # every non-QFT term of T^S is 4 x integer Toffolis
    n_a = eta + 2 * m + 10
    t_mcx = 4 * (n_a - 1)
    t_rot_q = _t_rot_2507(eps / (6 * (q + 1)))
    total = q * t_s + 2 * (q + 1) * t_p + 3 * (q + 1) * t_rot_q + 2 * q * t_mcx
    b_r2 = math.ceil(0.5 * math.log2(72 * math.pi ** 2 * lam * t / eps))
    b_qft = min(m - 1, math.ceil(math.log2(16 * d * lam * t * m / eps))) + 1
    qubits = eta * d * m + eta + (2 + d) * m + 13 + b_r2 + b_qft + max(2 * b_qft - 1, d * m + 6, eta + 2 * m + 9)
    return {"lambda_h": lam, "q": q, "eps_h": e1, "n_eta": n_eta, "b_r": b_r, "t_qft": qft,
            "toffoli_prep": tof_p, "t_prep": t_p, "t_select": t_s, "toffoli_select": tof_s,
            "toffoli_refl": n_a - 1, "t_mcx": t_mcx, "t_rot_qsp": t_rot_q, "t_total": total, "qubits": qubits,
            "rotations": 4 * 2 * (q + 1) + 3 * (q + 1)}


def _gqsp_repriced(a: Assumptions, eta: int, m: int, t: float, eps: float) -> dict:
    """R5 reprice of _gqsp_2507 (derived here): every Toffoli at 7 T (T_PER_TOFFOLI['textbook']); the QFT's
    2d T_QFT, built on 4-T phase-gradient adders (Nam et al.), scaled by 7/4 (an upper bound: at m = 3 it is
    60 T of 8272 per SELECT); the 11(Q+1) rotations (4 per PREPARE x 2(Q+1), 3 per GQSP SU(2) x (Q+1))
    synthesized under R-TOL at eps_rot = sqrt(1e-2 / N_rot) with the RUS fit. The source's eps still sets Q,
    b_r and the phase-gradient widths."""
    s = _gqsp_2507(a, eta, m, t, eps)
    tpt = T_PER_TOFFOLI["textbook"]
    d, q = a.fq_dim.lo, s["q"]
    n_rot = s["rotations"]
    eps_rot = eps_rot_for(n_rot, a.eps_syn.lo)
    t_rot = t_per_rotation(eps_rot, SYNTH_MODEL)
    sel_each = s["toffoli_select"] * tpt + 2 * d * s["t_qft"] * tpt / a.fq_t_src_per_toffoli.lo
    parts = (
        Primitive("gqsp_select", q, sel_each, CircuitStatus.COMPILED, "arxiv_2507_22814",
                  f"Thm. 3 T^S: {s['toffoli_select']:.0f} Toffoli x 7 T + QFT {2 * d * s['t_qft']:.0f} T x 7/4; "
                  "repriced at 7 T per Toffoli, derived here"),
        Primitive("gqsp_prepare", 2 * (q + 1), s["toffoli_prep"] * tpt, CircuitStatus.COMPILED, "arxiv_2507_22814",
                  f"Thm. 3 T^P: 3n_eta+4b_r+2m-16 = {s['toffoli_prep']} Toffoli x 7 T (rotations separate); derived here"),
        Primitive("gqsp_reflection", 2 * q, s["toffoli_refl"] * tpt, CircuitStatus.COMPILED, "arxiv_2507_22814",
                  f"Thm. 5: C^(eta+2m+10)X ladder, {s['toffoli_refl']} Toffoli x 7 T, measurement uncompute; derived here"),
        Primitive("rot_synth", n_rot, t_rot, CircuitStatus.COMPILED, "arxiv_2507_22814;" + SYNTHESIS_SRC,
                  f"11(Q+1) = {n_rot} R_z (4 per PREPARE, 3 per SU(2)); R-TOL eps_rot = {eps_rot:.4g}, {t_rot:.3f} T"),
    )
    return {**s, "n_rot": n_rot, "eps_rot": eps_rot, "t_rot": t_rot, "parts": parts,
            "t_total_repriced": sum(p.t_total for p in parts)}


def _q_coverage_GeV(a: Assumptions, L: int, box_fm: float) -> tuple[float, float]:
    """|q| reach of an L^3 periodic box: 2 pi hbar c / L_box times |n| in {1 .. sqrt(3) floor(L/2)}."""
    q_min = 2 * math.pi * a.hbar_c_MeV_fm.lo / box_fm       # MeV
    q_max = q_min * math.sqrt(3) * (L // 2)
    return q_min / 1e3, q_max / 1e3


def _walk_per_query(a: Assumptions, M: int) -> dict:
    """Controlled-qubiterate lower bound of arxiv_1911_06368 (appendix 'Gate cost of the qubiterate')
    at M sites: N_A = ceil(log2(38 M)); 2N_A-1 ancillas, 15*2^(N_A-1)+14N_A-37 T, 2^N_A-N_A-2 U(2)."""
    n_a = math.ceil(math.log2(a.walk_terms_per_site.lo * M))
    return {"n_a": n_a, "ancilla": 2 * n_a - 1,
            "t": 15 * 2 ** (n_a - 1) + 14 * n_a - 37,
            "u2": 2 ** n_a - n_a - 2}


def _replace_rot_per_u2(a: Assumptions, k: int) -> Assumptions:
    import dataclasses
    return dataclasses.replace(a, walk_rot_per_u2=Assumed(k, "comparison: 1911.06368's 3 R_z per U(2)",
                                                          "arxiv_1911_06368"))


def _walk_compiled(a: Assumptions, M: int, w_q: int) -> dict:
    """The walk of one shot: 2^W_q - 1 controlled steps; the U(2) are R_y, synthesized at the R-TOL
    tolerance of the walk's own rotation count (a lower bound: the rest of the shot's rotations are
    not counted in N_rot, which can only lower eps_rot and raise T per rotation slightly)."""
    q = _walk_per_query(a, M)
    n_query = 2 ** w_q - 1
    n_rot = n_query * q["u2"] * a.walk_rot_per_u2.lo
    eps = eps_rot_for(n_rot, a.eps_syn.lo)
    t_rot = t_per_rotation(eps, SYNTH_MODEL)
    per_query = q["t"] + q["u2"] * a.walk_rot_per_u2.lo * t_rot
    return {**q, "n_query": n_query, "n_rot": n_rot, "eps_rot": eps, "t_rot": t_rot,
            "t_per_query": per_query, "t_walk": n_query * per_query}


def _lambda_per_site_MeV(a: Assumptions) -> float:
    """1-norm per site of the arxiv_1911_06368 lattice pionless Hamiltonian in Pauli form (N_f = 4, 3D):
    kinetic 4 species x 6 neighbors x (XX+YY) x t/2 = 24t; Z: 4 x (3/4)|C0+D0|; ZZ: 6 pairs x |C0+2D0|/4;
    ZZZ: 4 triples x D0/4 (the identity term is dropped). 706.76 MeV at a = 1.4 fm."""
    t, c0, d0 = a.hop_t_MeV.lo, -a.c0_abs_MeV.lo, a.d0_MeV.lo
    return 24 * t + 3 * abs(c0 + d0) + 1.5 * abs(c0 + 2 * d0) + d0


def _n_shells(L: int) -> int:
    """Distinct nonzero |n|^2 reachable in an L^3 periodic box with |n_i| <= L//2: the momentum shells
    (3^3: |n|^2 = 1, 2, 3 -> 295, 417, 511 MeV at 4.2 fm)."""
    h = L // 2
    return len({x * x + y * y + z * z for x in range(-h, h + 1) for y in range(-h, h + 1)
                for z in range(-h, h + 1)} - {0})


def _bin_shots(p: float, eps: float) -> float:
    """QPE bin of weight p (multinomial): relative error eps after (1-p)/(p eps^2) shots."""
    return (1 - p) / (p * eps ** 2)


def _hadamard_shots(mu: float, eps: float) -> float:
    """Unamplified Hadamard test of a quantity mu in [-1, 1] (here C/alpha^2): (1-mu^2)/(eps mu)^2."""
    return (1 - mu ** 2) / (eps * mu) ** 2


def _depth_4he(a: Assumptions, w: dict, sel_depth: float) -> dict:
    """T-depth of one route-(ii) walk query (factory.json Ch. 2): PREPARE's level-k multiplexor is 2^k commuting
    rotations spread over K ancillas, ceil(2^k/K) layers of one synthesized R_y each (scaled by u2/(2^N_A-1),
    the model's 2^N_A-N_A-2 count); SELECT sequential over 2^N_A indices at sel_depth per index (1: measurement-
    uncomputed AND; 3: 7-T Toffoli, i.e. 3 x T_select/7); the controlled reflection 3 x (8N_A-9)/7."""
    n_a, K = w["n_a"], a.prepare_ancilla_4he.lo
    d_prep = sum(math.ceil(2 ** k / K) for k in range(n_a)) * w["t_rot"] * w["u2"] / (2 ** n_a - 1)
    d_sel = 2 ** n_a if sel_depth == 1 else sel_depth * w["t"] / 7
    d_refl = 3 * (8 * n_a - 9) / 7
    d_q = d_prep + d_sel + d_refl
    return {"d_prep": d_prep, "d_select": d_sel, "d_refl": d_refl, "d_query": d_q,
            "f_walk": w["t_per_query"] / d_q}


def _tau_max_from_sigma(sigma_MeV: float) -> float:
    """tau_max ~ 1/sigma in GeV^-1 (app01:90)."""
    return 1.0 / (sigma_MeV * 1e-3)


# --------------------------------------------------------------------------- #
# Eras
# --------------------------------------------------------------------------- #

def _model_2028(a: Assumptions) -> Result:
    n_orb = a.n_orb_2028.lo
    ch = _route_i(a, n_orb, a.n_trotter_2028.lo)              # 5e3 rotations, eps 1.414e-3, 20.086 T
    strings, t_step, t_chiral = ch["strings"], ch["t_step"], ch["t_shot"]   # 5e3, 1.00428e5, 1.00428e5
    # pionless: its own circuit, 1e3 strings/step (the Cited ~5x per-step gain read as a rotation
    # count, implied) x N_Trotter=3 = 3e3 rotations -> eps 1.826e-3, 19.662 T
    pl = _route_i(a, n_orb, a.n_trotter_pionless.lo, a.pionless_gain.lo)
    t_pionless_step, t_pionless = pl["t_step"], pl["t_shot"]  # 1.966e4, 5.899e4
    pl4 = _route_i(a, n_orb, 4, a.pionless_gain.lo)             # ceiling check: 7.960e4
    pl5 = _route_i(a, n_orb, 5, a.pionless_gain.lo)             # 1.00428e5, at the cap
    lq_sys = 4 * n_orb                                       # 40
    lq_lo = lq_sys + a.anc_qpe_window.lo + a.anc_lcu_be.lo + a.anc_hadamard.lo
    lq_hi = lq_sys + a.anc_qpe_window.hi + a.anc_lcu_be.hi + a.anc_hadamard.hi
    # r25 (R1, shot audit): J|Psi0> is not amplified (no room under 1e5 T); the Hadamard test estimates
    # C/alpha^2 = mu, so relative eps costs (1-mu^2)/(eps mu)^2 shots per (component, Re/Im part, tau point)
    mu = a.jpsi_norm_sq_2028.lo / a.alpha_2028.lo ** 2                       # 0.015
    n_first_each = _hadamard_shots(mu, a.eps_first.lo)                        # 4.937e4
    n_full_each = _hadamard_shots(mu, a.eps_target.lo)                        # 1.777e6
    shots_first = n_first_each * a.n_bin_2028.lo                              # Re C, R_L, one tau: 4.94e4
    shots = n_full_each * a.n_bin_2028.lo * a.n_comp_2028_full.lo * a.n_parts_2028.lo   # L+T, Re+Im: 7.11e6
    shots_bounded = _shots(a, a.n_bin_2028.lo)               # retired '~2e3' (+/-1 bound on the unnormalized C)
    # one machine, serial (E27): 1 us per T-gate (R8) plus t0 ~ 0.1 ms per shot
    wall = t_chiral * a.t_gate_s.lo + a.shot_overhead_s.lo   # 0.1005 s
    total = shots * wall                                     # 7.15e5 s = 8.27 d
    total_first = shots_first * wall                         # 4963 s = 82.7 min
    # T-depth (factory.json): W rotations per layer, W = 4 N_orb/4 = 10 .. 8 x 10 = 80
    w_lo = lq_sys / 4
    w_hi = w_lo * a.strings_per_excitation_2028.lo
    depth = (ch["n_rot"] / w_hi * ch["t_rot"], ch["n_rot"] / w_lo * ch["t_rot"])   # 1255.4, 1.0043e4
    dx = depth_exports(t_chiral, depth, shots, a.t_gate_s, a.shot_overhead_s)
    dx_first = depth_exports(t_chiral, depth, shots_first, a.t_gate_s, a.shot_overhead_s)
    eps_l_2028 = a.fault_budget.lo / t_chiral                # 9.96e-7, box '<~ 1e-6'
    rot = Primitive("rot_synth", ch["n_rot"], ch["t_rot"],
                    CircuitStatus.SCALING, "arxiv_2312_05344;" + SYNTHESIS_SRC,
                    f"N_orb^4/2 Pauli strings x N_Trotter={a.n_trotter_2028.lo}; R-TOL eps_rot="
                    f"sqrt({a.eps_syn.lo:g}/{ch['n_rot']:g})={ch['eps_rot']:.4g}, RUS fit {ch['t_rot']:.3f} T")
    inter = {
        "pauli_strings_per_step": strings,
        "rotations_per_step": strings,
        "rotations_per_shot_chiral": ch["n_rot"],
        "eps_syn": a.eps_syn.lo,
        "eps_rot_chiral": ch["eps_rot"],
        "t_per_rotation_chiral": ch["t_rot"],
        # literature: the compiled SU(3) circuits on disk price R_Z at 1.15 log2(1/eps), no offset
        "t_per_rotation_rus_slope_at_eps_rot": t_per_rotation(ch["eps_rot"], "rus-slope"),
        "t_per_step_chiral": t_step,
        "n_trotter_chiral": a.n_trotter_2028.lo,
        "t_per_shot_chiral": t_chiral,
        "over_cap_factor": t_chiral / a.rfi_2028_t.lo,
        # overshoot: 0.43% over the 1e5 reference; printed '~1e5, 0.4% over' (author, 2026-09-29)
        "over_cap_percent": 100 * (t_chiral / a.rfi_2028_t.lo - 1),
        "n_trotter_affordable_chiral": math.floor(a.rfi_2028_t.lo / t_step),
        "rotations_per_step_pionless": pl["strings"],
        "rotations_per_shot_pionless": pl["n_rot"],
        "eps_rot_pionless": pl["eps_rot"],
        "t_per_rotation_pionless": pl["t_rot"],
        "t_per_step_pionless": t_pionless_step,
        "n_trotter_pionless": a.n_trotter_pionless.lo,
        "t_per_shot_pionless": t_pionless,
        "pionless_fits_budget": t_pionless <= a.rfi_2028_t.lo,
        # step ceiling of the pionless circuit at the 1e5 reference (NOT the instance, which stays
        # at N_Trotter=3): N_T=4 fits outright (7.960e4); N_T=5 is 5e3 rotations, the chiral
        # circuit's count, so 1.00428e5, 0.43% over, printed 'at the cap' (app01:126)
        "t_per_shot_pionless_4": pl4["t_shot"],
        "t_per_shot_pionless_5": pl5["t_shot"],
        "pionless_4_fits_budget": pl4["t_shot"] <= a.rfi_2028_t.lo,
        "pionless_5_over_cap_percent": 100 * (pl5["t_shot"] / a.rfi_2028_t.lo - 1),
        "sigma_2028_MeV": a.sigma_2028_MeV.lo,
        "lq_system_jw": lq_sys,
        "lq_qpe_window_anc": a.anc_qpe_window.lo,
        "lq_lcu_be_anc": a.anc_lcu_be.lo,
        "lq_hadamard_anc": (a.anc_hadamard.lo, a.anc_hadamard.hi),
        "lq_ancilla_total": (lq_lo - lq_sys, lq_hi - lq_sys),
        "lq_total": (lq_lo, lq_hi),
        "lq_4he_prose": a.lq_2028_prose.lo,
        "lq_inside_rfi_2028_envelope": a.rfi_2028_lq.lo <= lq_lo and lq_hi <= a.rfi_2028_lq.hi,
        "inv_eps_sq": a.eps_target.lo ** -2,
        "shots_per_component": a.eps_target.lo ** -2,
        "shots": shots,
        "hadamard_mu": mu,
        "alpha_2028": a.alpha_2028.lo,
        "jpsi_amplified": False,
        "shots_per_part_first": n_first_each,
        "shots_per_part_full": n_full_each,
        "shots_first_result": shots_first,
        "shots_full_benchmark": shots,
        "shots_bounded_retired": shots_bounded,
        "wall_per_shot_s": wall,
        "shot_overhead_s": a.shot_overhead_s.lo,
        "shot_overhead_fraction": a.shot_overhead_s.lo / wall,
        "total_wall_h": total / 3600,
        "total_wall_min": total / 60,
        "total_wall_days": total / 86400,
        "first_result_wall_min": total_first / 60,
        # R9 exports (full benchmark run; first result alongside)
        "rot_layer_width": (w_lo, w_hi),
        "t_per_shot": dx["t_per_shot"],
        "t_depth_per_shot": dx["t_depth_per_shot"],
        "f_star": dx["f_star"],
        "floor_wall_s": dx["floor_wall_s"],
        "factories_for_1yr": dx["factories_for_1yr"],
        "baseline_ok": dx["baseline_ok"],
        "wall_serial_s": dx["wall_serial_s"],
        "wall_first_result_s": dx_first["wall_serial_s"],
        "wall_campaign_s": dx["wall_serial_s"],
        "floor_wall_first_result_s": dx_first["floor_wall_s"],
        "eps_l_from_fault_budget": eps_l_2028,
        "eps_l_2028_stated": a.eps_l_2028.lo,
        "eps_l_unit_fault": a.fault_budget_unit.lo / t_chiral,
        "expected_faults_at_rfi_floor": t_chiral * a.eps_l_2028_floor.lo,
        "rfi_2028_floor_is_ample": a.eps_l_2028_floor.lo <= eps_l_2028,
    }
    return Result(era="2028", lq=(lq_lo, lq_hi), hard_ops=(t_chiral, t_chiral),
                  breakdown=(rot,), intermediates=inter,
                  shots=(shots, shots), wall_time_s=(wall, wall), epsilon_l=(eps_l_2028, eps_l_2028),
                  notes=("R-TOL: chiral N_Trotter=1 is 5e3 rotations at eps_rot = sqrt(1e-2/5e3) = 1.414e-3, "
                         "20.086 T each: 1.00428e5 T, 0.43% over the 1e5 RFI cap (printed '~1e5, 0.4% over'). "
                         "Pionless (1e3 rotations/step implied, N_Trotter=3, eps 1.826e-3, 19.662 T): 5.90e4 T, fits; N_T=4 7.96e4 fits, N_T=5 1.004e5 at the cap.",
                         "2028 eps_l row: 0.1/1.00428e5 = 9.96e-7 (box '<~ 1e-6'); the RFI 1e-8 floor is ample.",
                         "Rotation count is a Pauli-string counting argument (SCALING), not a compiled circuit.",
                         "The compiled-convention slope (no offset) at the same eps_rot is emitted for comparison."))


def _model_2033(a: Assumptions) -> Result:
    cyc = a.t_gate_s.lo
    env = a.rfi_2033_t.lo
    t_lo, t_hi = a.t_4he.lo, a.t_4he.hi
    fb, fu = a.fault_budget.lo, a.fault_budget_unit.lo
    # --- what can be checked from the stated numbers -------------------------------
    eps_l = (fb / t_hi, fb / t_lo)
    t0 = a.shot_overhead_s.lo
    horizon_s = a.horizon_yr.lo * YEAR_S
    # one machine, serial (E27): shots x per-shot time; no machine count enters
    wall = (t_lo * cyc + t0, t_hi * cyc + t0)                # 670-960 s
    # r25 (R1, R2, shot audit): shots from the QPE bin statistics, 5% (campaign) / 30% (first result) on
    # every bin with p >= 0.1, on the box's momentum shells
    n_shells = _n_shells(a.lattice_4he.lo)                   # 3
    per_cell = _bin_shots(a.p_bin_min.lo, a.eps_target.lo)  # 3600 per (shell, component)
    per_cell_first = _bin_shots(a.p_bin_min.lo, a.eps_first.lo)   # 100
    shots_camp = n_shells * a.n_comp.lo * per_cell           # 5.4e4
    shots_first = n_shells * a.n_comp_first.lo * per_cell_first   # 300
    # normalization ||J Psi0||^2 from ground-state snapshots, one set per component; one snapshot = one
    # ground-state load: <= 1.5% of 6.7e8 (block rung) .. 28.5% of 9.6e8 (CI ceiling)
    snap_t = (a.gs_prep_frac_block.lo * t_lo, a.gs_prep_frac_ci.hi * t_hi)   # 1.005e7, 2.736e8
    snap_wall = (snap_t[0] * cyc + t0, snap_t[1] * cyc + t0)
    per_comp_snap = a.norm_snapshot_var.lo / a.eps_norm_campaign.lo ** 2                  # 1667
    # R_L at 1x; the other four components at norm_grouped_factor (1 .. 3)x: 8333 .. 2.17e4
    n_snap_camp_band = tuple(per_comp_snap * (1 + (a.n_comp.lo - 1) * g)
                             for g in (a.norm_grouped_factor.lo, a.norm_grouped_factor.hi))
    n_snap_camp = n_snap_camp_band[0]                                                    # 8333 (all Z-basis)
    n_snap_first = a.n_comp_first.lo * a.norm_snapshot_var.lo / a.eps_norm_first.lo ** 2  # 66.7
    norm_combined = math.hypot(a.eps_target.lo, a.eps_norm_campaign.lo)                  # 0.0539
    camp_s = tuple(shots_camp * w + n_snap_camp * sw for w, sw in zip(wall, snap_wall))     # 3.626e7, 5.412e7 (depth-1, Z-basis)
    first_s = tuple(shots_first * w + n_snap_first * sw for w, sw in zip(wall, snap_wall))  # 2.017e5, 3.062e5
    camp_yr = (camp_s[0] / YEAR_S, camp_s[1] / YEAR_S)       # 1.149-1.715 yr
    # T-depth (factory.json): walk F* at the depth-1 / depth-3 SELECT; the remainder at the walk's parallelism
    w3d = _walk_compiled(a, a.lattice_4he.lo ** 3, a.w_q.lo)
    dep1 = _depth_4he(a, w3d, a.select_depth_per_index.lo)  # F_walk 14.29
    dep3 = _depth_4he(a, w3d, a.select_depth_per_index.hi)  # F_walk 7.32
    d_shot = (t_lo / dep1["f_walk"], t_hi / dep3["f_walk"])  # 4.69e7, 1.311e8
    dx = depth_exports((t_lo, t_hi), d_shot, shots_camp, a.t_gate_s, a.shot_overhead_s)
    dx_first = depth_exports((t_lo, t_hi), d_shot, shots_first, a.t_gate_s, a.shot_overhead_s)
    # the box assumes the depth-1 SELECT (F* >= 10); with the depth-3 SELECT a shot is reaction-limited
    shot_d3 = (max(t_lo * cyc, t_lo / dep3["f_walk"] * REACTION_TIME_S) + t0,
               max(t_hi * cyc, t_hi / dep3["f_walk"] * REACTION_TIME_S) + t0)            # 915, 1311 s
    camp_d3 = tuple(shots_camp * w + n_snap_camp * sw for w, sw in zip(shot_d3, snap_wall))
    first_d3 = tuple(shots_first * w + n_snap_first * sw for w, sw in zip(shot_d3, snap_wall))
    # CONTRACT rule 12 (f_star < 10 at the depth-3 SELECT): the box quotes the corrected wall, never the
    # 10-factory baseline. Low end: 6.7e8 T, depth-1 SELECT (F* 14.3, baseline holds), Z-basis normalization.
    # High end: 9.6e8 T, depth-3 SELECT (F* 7.3, reaction-limited), grouped normalization at 3x.
    shot_box = (wall[0], shot_d3[1])                                                    # 670, 1311 s
    first_box = (first_s[0], first_d3[1])                                               # 2.017e5, 4.115e5 s
    camp_box = (camp_s[0], shots_camp * shot_d3[1] + n_snap_camp_band[1] * snap_wall[1])   # 3.626e7, 7.673e7 s
    camp_box_yr = (camp_box[0] / YEAR_S, camp_box[1] / YEAR_S)                          # 1.149, 2.431 yr
    # depth band per T end and per SELECT compile, not cross-paired across T ends (physical F* 7.32-14.29)
    dx_pairs = [depth_exports(t_, t_ / f_, shots_camp, a.t_gate_s, a.shot_overhead_s)
                for t_ in (t_lo, t_hi) for f_ in (dep1["f_walk"], dep3["f_walk"])]
    f_star_phys = (min(d["f_star"][0] for d in dx_pairs), max(d["f_star"][1] for d in dx_pairs))
    # R11b: lambda derived from the cited 1911.06368 couplings; the printed 19.1 GeV is checked to 1%
    lam_site = _lambda_per_site_MeV(a)                      # 706.76 MeV
    lam = lam_site * a.lattice_4he.lo ** 3                   # 19.083 GeV on 3^3
    if abs(a.lambda_GeV.lo * 1e3 / lam - 1) > 0.01:
        raise ValueError(f"printed lambda {a.lambda_GeV.lo} GeV disagrees with M x lambda/site = {lam / 1e3:.3f} GeV")
    two_wq = 2 ** a.w_q.lo
    bin_from_lambda = 2 * math.pi * lam / two_wq
    lam_implied = a.bin_4he_MeV.lo * two_wq / (2 * math.pi)  # what 29.3 MeV would need
    step_reg = lam / two_wq                                  # 4.66 MeV: lambda/2^Wq, the bin without its 2 pi
    # ruling ch02-sigma-window-vs-register (a): compare bin and Gaussian window at equal FWHM
    fwhm_per_sigma = 2 * math.sqrt(2 * math.log(2))          # 2.3548
    sigma_eq_fwhm = a.bin_4he_MeV.lo / fwhm_per_sigma         # 12.44 MeV: sigma a 29.3 MeV bin matches
    bin_plus_one_bit = bin_from_lambda / 2                    # 14.65 MeV
    extra_bits = math.ceil(math.log2(a.bin_4he_MeV.lo / a.bin_target_MeV.lo))
    block_t = (a.block_frac.lo * env, a.block_frac.hi * env)
    plus_bit_from_walk = (1 - a.walk_frac.lo) + 2 * a.walk_frac.lo
    block_plus_bit = (a.block_frac.lo * a.plus_bit_factor.lo, a.block_frac.hi * a.plus_bit_factor.hi)
    plus_bit_t = (t_lo * a.plus_bit_factor.lo, t_hi * a.plus_bit_factor.hi)
    plus_bit_t_walk = (t_lo * plus_bit_from_walk, t_hi * plus_bit_from_walk)
    plus_cube_t = (t_lo * a.plus_cube_factor.lo, t_hi * a.plus_cube_factor.hi)
    tau_grid_t = (a.tau_grid_frac.lo * env, a.tau_grid_frac.hi * env)
    gs_ci_t = (a.gs_prep_frac_ci.lo * t_hi, a.gs_prep_frac_ci.hi * t_hi)
    # lattice geometry and kinematic coverage
    L3, L4 = a.lattice_4he.lo, a.lattice_12c.lo
    a_fm = a.box_4he_fm.lo / L3                              # 1.4 fm
    box_12c_from_a = a_fm * L4                               # 5.6 fm
    box_5_from_a = a_fm * 5                                  # 7.0 fm
    q3 = _q_coverage_GeV(a, L3, a.box_4he_fm.lo)
    q4 = _q_coverage_GeV(a, L4, a.box_12c_fm.lo)
    vol_ratio_43 = (L4 ** 3) / (L3 ** 3)                     # 64/27
    vol_ratio_54 = (5 ** 3) / (L4 ** 3)                      # 125/64
    # site-basis registers; a stated peak register smaller than its site register is nonsense input
    reg3, reg4 = _site_register(a, L3), _site_register(a, L4)
    for label, lq_stated, reg in (("4He", a.lq_4he.lo, reg3), ("4He tau-grid", a.tau_grid_lq.lo, reg3),
                                  ("12C", a.lq_12c.lo, reg4), ("12C on-ramp", a.lq_onramp.lo, reg4)):
        if lq_stated < reg:
            raise ValueError(f"{label}: stated peak register {lq_stated} LQ is smaller than the "
                             f"{a.modes_per_site.lo} x L^3 = {reg} site register it must hold")
    # 12C: lambda scales with volume from the STATED lambda (not from the bin, which was circular)
    lam_12c = lam * vol_ratio_43
    bin_12c_at_wq = 2 * math.pi * lam_12c / two_wq
    wq_12c_implied = math.log2(2 * math.pi * lam_12c / a.bin_12c_MeV.lo)
    eps_12c = (fb / a.t_12c.hi, fb / a.t_12c.lo)             # 0.1 fault at 2.8e9 .. 2.7e9 T (R3)
    eps_12c_unit = (fu / a.t_12c.hi, fu / a.t_12c.lo)        # the pre-R3 unit-fault figure, comparison only
    bin_12c_at_stated_wq = 2 * math.pi * lam_12c / 2 ** a.w_q_12c.lo
    over_12c = (a.t_12c.lo / env, a.t_12c.hi / env)
    five_cube_t = (a.five_cube_factor.lo * env, a.five_cube_factor.hi * env)
    five_cube_if_v2 = (over_12c[0] * vol_ratio_54 ** 2, over_12c[1] * vol_ratio_54 ** 2)
    # second candidate law: one W_q bit (x2) times one LCU-index bit (38 x 125 = 4750 > 4096: N_A 12 -> 13,
    # the 1911.06368 per-query T goes 30851 -> 61585)
    lcu_bit_ratio = _walk_per_query(a, 125)["t"] / _walk_per_query(a, L4 ** 3)["t"]
    five_cube_if_bit_lcu = (over_12c[0] * 2 * lcu_bit_ratio, over_12c[1] * 2 * lcu_bit_ratio)
    lam_5 = lam_site * 125
    bin_5cube = 2 * math.pi * lam_5 / 2 ** a.w_q_5cube.lo   # 33.88 MeV at W_q = 14
    # on-ramp (R11b fallback): 'a few x 1e7 T', not separately costed; eps_l printed at the 3e7 reference
    eps_onramp = (fb / a.t_onramp.hi, fb / a.t_onramp.lo)    # 2e-9 .. 5e-9
    eps_onramp_ref = fb / a.t_onramp_eps_ref.lo              # 3.33e-9, printed '<~ 3e-9 at 3e7 T'
    euclid_t = (a.t_onramp.lo * a.euclid_factor.lo, a.t_onramp.hi * a.euclid_factor.hi)
    # r25 (shot audit): depth-capped amplitude estimation (MLAE) on the on-ramp at the 3e7 reference: Grover depth
    # 2k+1 up to the 1e9 envelope; the deepest circuit then needs eps_l ~ 1e-10
    mlae_kmax = math.floor((env / a.t_onramp_eps_ref.lo - 1) / 2)          # 16
    mlae_deepest = (2 * mlae_kmax + 1) * a.t_onramp_eps_ref.lo            # 9.9e8
    # 40Ar
    M_ar = a.lattice_ar.lo ** 3
    fq_ar = _fq_register(a.A_ar.lo, M_ar)
    site_ar = _site_register(a, a.lattice_ar.lo)
    # R12: route-(iii) evolution from arXiv:2507.22814 Thm. 5, reproduced then repriced (derived here)
    eta, m_ar, eps_fq = a.A_ar.lo, int(math.log2(a.lattice_ar.lo)), a.fq_eps.lo
    tdw = _gqsp_t_dw(a, eta, m_ar, a.fq_dw_MeV.lo)             # t = 154 x 2pi/15419.5 = 0.06275 MeV^-1
    tdw10 = _gqsp_t_dw(a, eta, m_ar, a.fq_dw_extrap_MeV.lo)    # 1541 x 2pi/15419.5 = 0.6279 MeV^-1
    tcr = _gqsp_t_cross(a, m_ar)                               # 0.3889 MeV^-1
    src100 = _gqsp_2507(a, eta, m_ar, tdw["t"], eps_fq)        # 3.389e7 T, 498 qubits (Tab. 5: 3.38e7, 498)
    src10 = _gqsp_2507(a, eta, m_ar, tdw10["t"], eps_fq)       # 3.408e8 T, 500 qubits (extrapolation)
    srccr = _gqsp_2507(a, eta, m_ar, tcr, eps_fq)              # 2.111e8 T, 500 qubits (Tab. 5: 2.11e8, 500)
    src16 = (_gqsp_2507(a, 16, m_ar, _gqsp_t_dw(a, 16, m_ar, a.fq_dw_MeV.lo)["t"], eps_fq),
             _gqsp_2507(a, 16, m_ar, tcr, eps_fq))             # 5.996e6 / 3.743e7 (Tab. 5: 5.99e6 / 3.75e7)
    rep100 = _gqsp_repriced(a, eta, m_ar, tdw["t"], eps_fq)    # 5.875e7 T
    rep10 = _gqsp_repriced(a, eta, m_ar, tdw10["t"], eps_fq)   # 5.905e8 T
    repcr = _gqsp_repriced(a, eta, m_ar, tcr, eps_fq)          # 3.657e8 T
    # the reproduction gate: all eight Tab. 5 GQSP entries within 5% (T) and exactly (qubits), else no reprice
    for lab, got, tab, q_got, q_tab in (
            ("eta=40 t_dw", src100["t_total"], a.tab5_gqsp_t_dw.lo, src100["qubits"], a.tab5_gqsp_q_dw.lo),
            ("eta=40 t_cross", srccr["t_total"], a.tab5_gqsp_t_cross.lo, srccr["qubits"], a.tab5_gqsp_q_cross.lo),
            ("eta=16 t_dw", src16[0]["t_total"], a.tab5_gqsp16_t_dw.lo, src16[0]["qubits"], a.tab5_gqsp16_q.lo),
            ("eta=16 t_cross", src16[1]["t_total"], a.tab5_gqsp16_t_cross.lo, src16[1]["qubits"], a.tab5_gqsp16_q.hi)):
        if abs(got / tab - 1) > 0.05 or q_got != q_tab:
            raise ValueError(f"2507.22814 Tab. 5 GQSP {lab} not reproduced by Thm. 5: {got:.4g} T / {q_got} vs {tab:.4g} / {q_tab}")
    t_ar_ev = (rep100["t_total_repriced"], rep10["t_total_repriced"])
    # R12: register with spin-isospin qubits (source Q = (dm+2) eta); Tab. 5 counts only eta d m
    lq_ar_ev = (src100["qubits"] + eta * a.spin_isospin_bits.lo, src10["qubits"] + eta * a.spin_isospin_bits.lo)
    if abs(lq_ar_ev[0] - a.lq_ar.lo) > 5 or lq_ar_ev[1] > a.lq_ar.hi:
        raise ValueError(f"40Ar register {lq_ar_ev} does not match the printed ~{a.lq_ar.lo}-{a.lq_ar.hi}")
    # the Stated composed factors were set on the source's 4-T evolution; the reprice shift per evolution
    ar_shift = (rep100["t_total_repriced"] - src100["t_total"], rep10["t_total_repriced"] - src10["t_total"])
    ar_amp = (a.ar_amplified_factor.lo * env, a.ar_amplified_factor.hi * env)
    ar_post = a.ar_postselected_factor.lo * env
    eps_ar_unit = (fu / ar_amp[1], fu / ar_amp[0])
    eps_ar_budget = (fb / ar_amp[1], fb / ar_amp[0])   # 1.1e-11 .. 6.7e-10 (R3)
    fq_bits_4he = _fq_register(4, L3 ** 3) + 4 * a.spin_isospin_bits.lo          # 28
    fq_bits_12c = _fq_register(a.A_c.lo, L4 ** 3) + a.A_c.lo * a.spin_isospin_bits.lo   # 96
    # cross-era ancilla
    anc_2028 = a.anc_qpe_window.lo + a.anc_lcu_be.lo + a.anc_hadamard.lo   # 131
    # route-(ii) registers from 1911.06368: site register + (W_q | 1 Hadamard) + (2N_A - 1)
    wk3, wk4 = _walk_per_query(a, L3 ** 3), _walk_per_query(a, L4 ** 3)   # N_A = 11, 12
    reg_4he = reg3 + a.w_q.lo + wk3["ancilla"]                # 108 + 12 + 21 = 141
    reg_tau = reg3 + 1 + wk3["ancilla"]                       # 108 + 1 + 21 = 130
    reg_12c = reg4 + a.w_q_12c.lo + wk4["ancilla"]            # 256 + 13 + 23 = 292
    reg_onramp = reg4 + 1 + wk4["ancilla"]                    # 256 + 1 + 23 = 280
    for label, printed, derived in (("4He", a.lq_4he.lo, reg_4he), ("4He tau-grid", a.tau_grid_lq.lo, reg_tau),
                                    ("12C", a.lq_12c.lo, reg_12c), ("12C on-ramp", a.lq_onramp.lo, reg_onramp)):
        if printed != derived:
            raise ValueError(f"{label}: printed {printed} LQ != 1911.06368 components {derived}")
    # compiled walk (lower bound) per shot: 4He W_q=12, 12C W_q=13
    w3 = _walk_compiled(a, L3 ** 3, a.w_q.lo)                 # 4095 x 6.887e4 = 2.820e8
    w4 = _walk_compiled(a, L4 ** 3, a.w_q_12c.lo)             # 8191 x 1.427e5 = 1.169e9
    # breakdown at the lo end of the stated composed range, split by the stated shares
    rep = t_lo
    gs = Primitive("gs_prep_load", 1, a.gs_prep_frac_block.lo * rep, CircuitStatus.UNSOURCED,
                   f"{TEX}:139", "'<= 1.5% of the composed cost at and below the block rung'; no circuit shown")
    walk_t = Primitive("walk_qubiterate_t", w3["n_query"], w3["t"], CircuitStatus.COMPILED, "arxiv_1911_06368",
                       f"lower bound per controlled qubiterate (appendix 'Gate cost of the qubiterate'): "
                       f"15*2^(N_A-1)+14N_A-37 = {w3['t']} T at N_A={w3['n_a']} (Gamma=38x27); "
                       f"2^W_q-1 = {w3['n_query']} steps at W_q={a.w_q.lo}; Clifford controlled unitaries neglected")
    walk_rot = Primitive("walk_prepare_ry", w3["n_rot"], w3["t_rot"], CircuitStatus.COMPILED,
                         "arxiv_1911_06368;" + SYNTHESIS_SRC,
                         f"2^N_A-N_A-2 = {w3['u2']} U(2) per step, one R_y each (real amplitudes, ruling "
                         f"ch02-no-compiled-walk (i)); R-TOL eps_rot=sqrt({a.eps_syn.lo:g}/{w3['n_rot']:g})="
                         f"{w3['eps_rot']:.4g}, RUS fit {w3['t_rot']:.3f} T")
    rest = Primitive("current_insertion_amplification_readout", 1,
                     (1 - a.gs_prep_frac_block.lo) * rep - w3["t_walk"], CircuitStatus.UNSOURCED,
                     f"{TEX}:141", "'the operator insertion and amplification schedule dominate'; remainder after "
                                   "gs prep and the compiled walk, no circuit shown")
    inter = {
        "envelope_fraction_4he": (t_lo / env, t_hi / env),
        "block_rung_t": block_t,
        "eps_l_from_fault_budget": eps_l,
        "eps_l_4he_stated": a.eps_l_4he.lo,
        "eps_l_4he_unit_fault": (fu / t_hi, fu / t_lo),
        "wall_per_shot_h": (wall[0] / 3600, wall[1] / 3600),
        "wall_per_shot_min": (wall[0] / 60, wall[1] / 60),
        "shot_overhead_s": t0,
        # r25 (R1, R2): first result and campaign tiers, single machine, serial (E27)
        "n_shells": n_shells,
        "bin_weight_floor": a.p_bin_min.lo,
        "shots_per_shell_comp_campaign": per_cell,
        "shots_per_shell_comp_first": per_cell_first,
        "shots_campaign": shots_camp,
        "shots_first_result": shots_first,
        "norm_snapshot_t": snap_t,
        "norm_snapshots_campaign": n_snap_camp,
        "norm_snapshots_first": n_snap_first,
        "norm_snapshots_per_component_2pct": a.norm_snapshot_var.lo / a.eps_norm_campaign.lo ** 2,   # 1667
        "norm_wall_campaign_s": (n_snap_camp * snap_wall[0], n_snap_camp * snap_wall[1]),
        "norm_wall_first_s": (n_snap_first * snap_wall[0], n_snap_first * snap_wall[1]),
        "error_combined_campaign": norm_combined,
        "first_result_days": (first_s[0] / 86400, first_s[1] / 86400),
        "campaign_serial_yr": camp_yr,
        "norm_snapshots_campaign_grouped": n_snap_camp_band,                              # 8333 .. 2.17e4
        "wall_per_shot_box_s": shot_box,
        "wall_per_shot_box_min": (shot_box[0] / 60, shot_box[1] / 60),                   # 11.2 .. 21.9
        "first_result_days_box": (first_box[0] / 86400, first_box[1] / 86400),           # 2.33 .. 4.76
        "campaign_yr_box": camp_box_yr,
        "campaign_box_fits_horizon": (camp_box[0] <= horizon_s, camp_box[1] <= horizon_s),
        "campaign_serial_months": (camp_s[0] / MONTH_S, camp_s[1] / MONTH_S),
        "horizon_yr": a.horizon_yr.lo,
        "campaign_fits_horizon": (camp_s[0] <= horizon_s, camp_s[1] <= horizon_s),      # (True, True)
        "campaign_over_horizon": (camp_s[0] / horizon_s, camp_s[1] / horizon_s),        # 0.230, 0.343
        # R9 exports (campaign run; first-result run alongside)
        "t_per_shot": dx["t_per_shot"],
        "t_depth_per_shot": dx["t_depth_per_shot"],
        "f_star": dx["f_star"],
        "floor_wall_s": dx["floor_wall_s"],
        "factories_for_1yr": dx["factories_for_1yr"],
        "baseline_ok": dx["baseline_ok"],
        "fits_1yr": dx["fits_1yr"],
        # rule 12: the box walls (corrected at the depth-3 end), not the 10-factory baseline
        "wall_first_result_s": first_box,
        "wall_campaign_s": camp_box,
        "wall_first_result_s_baseline": first_s,
        "wall_campaign_s_baseline": camp_s,
        # f_star above is depth_exports' cross-paired band (t_lo/d_hi, t_hi/d_lo); with d_shot pairing the
        # depth-1 compile at t_lo and the depth-3 compile at t_hi, its 5.1 and 20.5 are an artifact of mixing
        # the two T ends. The physical per-copy band is f_star_physical (factory.json '7.3-14.3 per copy').
        "f_star_physical": f_star_phys,
        "floor_wall_first_result_s": dx_first["floor_wall_s"],
        "factories_for_1yr_first_result": dx_first["factories_for_1yr"],
        "walk_f_star_select_depth1": dep1["f_walk"],
        "walk_f_star_select_depth3": dep3["f_walk"],
        "walk_depth_per_query": (dep1["d_query"], dep3["d_query"]),
        "prepare_depth_per_query": dep1["d_prep"],
        "t_depth_per_shot_paired": ((t_lo / dep1["f_walk"], t_lo / dep3["f_walk"]),
                                    (t_hi / dep1["f_walk"], t_hi / dep3["f_walk"])),
        "wall_per_shot_depth3_select_s": shot_d3,
        "wall_first_result_s_depth3_select": first_d3,
        "wall_campaign_s_depth3_select": camp_d3,
        "campaign_yr_depth3_select": (camp_d3[0] / YEAR_S, camp_d3[1] / YEAR_S),
        "depth3_stretch": (shot_d3[0] / wall[0], shot_d3[1] / wall[1]),
        "lambda_per_site_MeV": lam_site,
        "lambda_3cube_GeV": lam / 1e3,
        "lambda_4cube_GeV": lam_site * L4 ** 3 / 1e3,
        "lambda_5cube_GeV": lam_5 / 1e3,
        "lambda_printed_GeV": a.lambda_GeV.lo,
        "bin_from_lambda_wq_MeV": bin_from_lambda,
        "lambda_implied_by_bin_GeV": lam_implied / 1e3,
        "register_step_MeV": step_reg,                       # app01:144 'lambda/2^Wq ~ 4.7 MeV is the register step'
        "sigma_dune_fwhm_MeV": fwhm_per_sigma * a.sigma_dune_MeV.lo,   # 23.5
        "sigma_matching_bin_equal_fwhm_MeV": sigma_eq_fwhm,   # 12.4, printed '~12'
        "bin_std_uniform_MeV": a.bin_4he_MeV.lo / math.sqrt(12),   # 8.46 (equal-variance reading, not printed)
        "bin_plus_one_bit_MeV": bin_plus_one_bit,             # 14.6
        "sigma_matching_plus_one_bit_MeV": bin_plus_one_bit / fwhm_per_sigma,   # 6.2 <= 10: met
        "bin_meets_sigma_dune": sigma_eq_fwhm <= a.sigma_dune_MeV.lo,
        "plus_one_bit_meets_sigma_dune": bin_plus_one_bit / fwhm_per_sigma <= a.sigma_dune_MeV.lo,
        "bin_plus_three_bits_MeV": bin_from_lambda / 8,       # 3.66 <= 5
        "extra_bits_for_5MeV": extra_bits,
        "plus_bit_factor_from_walk_share": plus_bit_from_walk,
        "t_4he_plus_one_bit": plus_bit_t,
        "t_4he_plus_one_bit_from_walk_share": plus_bit_t_walk,
        "block_plus_bit_frac": block_plus_bit,
        "block_plus_bit_frac_stated": (a.block_plus_bit_frac.lo, a.block_plus_bit_frac.hi),
        "plus_cube_t": plus_cube_t,
        "plus_cube_exceeds_envelope": plus_cube_t[0] > env,
        "plus_cube_factor_if_v_squared": vol_ratio_43 ** 2,
        "gs_prep_t_ci_ceiling": gs_ci_t,
        "tau_grid_t": tau_grid_t,
        "lattice_spacing_fm": a_fm,
        "box_12c_fm_from_spacing": box_12c_from_a,
        "box_5cube_fm_from_spacing": box_5_from_a,
        "box_5cube_fm_stated": a.box_5cube_fm.lo,
        "q_coverage_3cube_GeV": q3,
        "q_coverage_3cube_stated_GeV": (a.q_coverage_3cube_GeV.lo, a.q_coverage_3cube_GeV.hi),
        "q_coverage_4cube_GeV": q4,
        "q_coverage_4cube_stated_GeV": (a.q_coverage_4cube_GeV.lo, a.q_coverage_4cube_GeV.hi),
        "lq_4he_site_register": reg3,
        "lq_4he_ancilla_implied": a.lq_4he.lo - reg3,
        "lq_4he_tau_grid_ancilla_implied": a.tau_grid_lq.lo - reg3,
        "lq_2028_ancilla_total": anc_2028,
        "walk_n_a_4he": wk3["n_a"],
        "walk_ancilla_4he": wk3["ancilla"],
        "walk_n_a_12c": wk4["n_a"],
        "walk_ancilla_12c": wk4["ancilla"],
        "lq_4he_from_components": reg_4he,
        "lq_4he_tau_grid_from_components": reg_tau,
        "lq_12c_from_components": reg_12c,
        "lq_12c_onramp_from_components": reg_onramp,
        # compiled walk lower bound (1911.06368 appendix), per query and per shot
        "walk_t_per_query_4he": w3["t"],
        "walk_u2_per_query_4he": w3["u2"],
        "walk_queries_4he": w3["n_query"],
        "walk_rotations_per_shot_4he": w3["n_rot"],
        "walk_eps_rot_4he": w3["eps_rot"],
        "walk_t_per_rotation_4he": w3["t_rot"],
        "walk_t_per_query_total_4he": w3["t_per_query"],
        "walk_t_per_shot_4he": w3["t_walk"],
        "walk_share_4he_compiled": (w3["t_walk"] / t_hi, w3["t_walk"] / t_lo),   # 0.29-0.42
        "walk_share_brackets_stated": w3["t_walk"] / t_hi <= a.walk_frac.lo <= w3["t_walk"] / t_lo,
        "walk_t_per_query_total_12c": w4["t_per_query"],
        "walk_t_per_shot_12c": w4["t_walk"],
        "walk_share_12c_compiled": (w4["t_walk"] / a.t_12c.hi, w4["t_walk"] / a.t_12c.lo),   # 0.42-0.43
        # option (ii) of the ruling, for comparison only: the paper's 3 R_z per U(2)
        "walk_t_per_shot_4he_3rz": _walk_compiled(
            _replace_rot_per_u2(a, 3), L3 ** 3, a.w_q.lo)["t_walk"],
        "ancilla_ratio_2028_over_4he_implied": anc_2028 / (a.lq_4he.lo - reg3),
        "lq_12c_site_register": reg4,
        "lq_12c_ancilla_implied": a.lq_12c.lo - reg4,
        "lq_12c_onramp_ancilla_implied": a.lq_onramp.lo - reg4,
        "lq_12c_5cube_site_register": _site_register(a, 5),
        "lambda_12c_volume_scaled_GeV": lam_12c / 1e3,
        "lambda_12c_from_bin_implied_GeV": lam_implied * vol_ratio_43 / 1e3,
        "bin_12c_at_wq12_MeV": bin_12c_at_wq,
        "wq_12c_implied_by_bin": wq_12c_implied,
        "wq_12c_stated": a.w_q_12c.lo,
        "bin_12c_at_stated_wq_MeV": bin_12c_at_stated_wq,
        "eps_l_12c_from_fault_budget": eps_12c,
        "eps_l_12c_unit_fault": eps_12c_unit,
        "eps_l_12c_stated": a.eps_l_12c.lo,
        "over_envelope_12c": over_12c,
        "five_cube_t": five_cube_t,
        "five_cube_factor_stated": (a.five_cube_factor.lo, a.five_cube_factor.hi),
        "five_cube_factor_if_v_squared": five_cube_if_v2,
        "five_cube_factor_if_bit_x_lcu_bit": five_cube_if_bit_lcu,
        "five_cube_lcu_bit_ratio": lcu_bit_ratio,
        "five_cube_bin_MeV": bin_5cube,
        "five_cube_wq": a.w_q_5cube.lo,
        "te_pe_lever_t": (a.t_12c.lo / a.te_pe_gain.lo, a.t_12c.hi / a.te_pe_gain.lo),
        "t_onramp": (a.t_onramp.lo, a.t_onramp.hi),
        "eps_l_onramp_from_fault_budget": eps_onramp,
        "eps_l_onramp_at_ref": eps_onramp_ref,
        "eps_l_onramp_stated": a.eps_l_onramp.lo,
        "euclid_variant_t": euclid_t,
        "onramp_mlae_kmax": mlae_kmax,
        "onramp_mlae_deepest_t": mlae_deepest,
        "onramp_mlae_faults_at_stated_eps_l": mlae_deepest * eps_onramp_ref,   # 3.3
        "onramp_mlae_eps_l_needed": fb / mlae_deepest,                         # 1.01e-10
        "crossover_A_log2M_4he": _fq_register(4, L3 ** 3),
        "crossover_M_4he": L3 ** 3,
        "crossover_4M_4he": reg3,
        "crossover_fq_with_spin_isospin_4he": fq_bits_4he,
        "crossover_A_log2M_12c": _fq_register(a.A_c.lo, L4 ** 3),
        "crossover_M_12c": L4 ** 3,
        "crossover_4M_12c": reg4,
        "crossover_fq_with_spin_isospin_12c": fq_bits_12c,
        "fq_wins_register_4he": fq_bits_4he < reg3,
        "lq_ar_fq_system": fq_ar,
        "lq_ar_fq_with_anc": (fq_ar + a.n_anc_bracket.lo, fq_ar + a.n_anc_bracket.hi),
        "lq_ar_site_basis_8cubed": site_ar,
        "lq_ar_evolution_with_spin_isospin": lq_ar_ev,          # (578, 580), printed '~580'
        "ar_composed_reprice_shift_per_evolution": ar_shift,    # (2.49e7, 2.50e8), printed '2.5e7 ... 2.5e8'
        "ar_composed_one_evolution_repriced": (ar_amp[0] + ar_shift[0], ar_amp[1] + ar_shift[1]),   # not printed
        # app01:114 (R12): the source's Tab. 5 comparison, second-quantized Trotter [2312.05344] vs GQSP
        "tab5_trotter2q_over_gqsp_t_dw": a.tab5_trotter2q_t_dw.lo / a.tab5_gqsp_t_dw.lo,        # 10.27 -> '10x'
        "tab5_trotter2q_over_gqsp_q_dw": a.tab5_trotter2q_q.lo / a.tab5_gqsp_q_dw.lo,           # 6.17 -> '6.2x'
        "tab5_trotter2q_over_gqsp_t_cross": a.tab5_trotter2q_t_cross.lo / a.tab5_gqsp_t_cross.lo,   # 27.9 -> '28x'
        "tab5_trotter2q_over_gqsp_repriced_t_dw": a.tab5_trotter2q_t_dw.lo / rep100["t_total_repriced"],   # 5.9x
        # R12 GQSP reproduction and reprice (derived here from Thm. 5 / Eq. 138)
        "fq_delta_h_MeV": tdw["delta_h"], "fq_norm_t_MeV": tdw["norm_t"], "fq_norm_v_MeV": tdw["norm_v"],
        "fq_lambda_h_MeV": src100["lambda_h"],
        "fq_t_dw100": tdw["t"], "fq_t_dw10": tdw10["t"], "fq_t_cross": tcr,
        "fq_t_ratio_10_to_100": tdw10["t"] / tdw["t"],
        "fq_q_dw100": src100["q"], "fq_q_dw10": src10["q"], "fq_q_cross": srccr["q"],
        "fq_src_t_dw100": src100["t_total"], "fq_src_t_dw10": src10["t_total"], "fq_src_t_cross": srccr["t_total"],
        "fq_src_qubits": (src100["qubits"], src10["qubits"], srccr["qubits"]),
        "fq_src_eta16": (src16[0]["t_total"], src16[1]["t_total"], src16[0]["qubits"], src16[1]["qubits"]),
        "fq_src_over_tab5_dw100": src100["t_total"] / a.tab5_gqsp_t_dw.lo,
        "fq_src_over_tab5_cross": srccr["t_total"] / a.tab5_gqsp_t_cross.lo,
        "fq_components_dw100": {k: src100[k] for k in ("t_select", "toffoli_select", "t_qft", "t_prep", "toffoli_prep",
                                                       "b_r", "n_eta", "toffoli_refl", "t_mcx", "t_rot_qsp", "rotations")},
        "fq_components_dw10": {k: src10[k] for k in ("toffoli_prep", "b_r", "t_rot_qsp", "rotations")},
        "fq_rep_t_dw100": rep100["t_total_repriced"], "fq_rep_t_dw10": rep10["t_total_repriced"],
        "fq_rep_t_cross": repcr["t_total_repriced"],
        "fq_rep_eps_rot": (rep100["eps_rot"], rep10["eps_rot"]), "fq_rep_t_rot": (rep100["t_rot"], rep10["t_rot"]),
        "fq_rep_over_src": rep100["t_total_repriced"] / src100["t_total"],
        "ar_evolution_breakdown_dw100": rep100["parts"],
        "ar_evolution_breakdown_dw10": rep10["parts"],
        "ar_evolution_t": t_ar_ev,
        "eps_l_ar_evolution_only": (fb / t_ar_ev[1], fb / t_ar_ev[0]),   # not printed (composed rows are)
        "lq_ar_requirements_summary": (a.lq_ar_requirements.lo, a.lq_ar_requirements.hi),
        "ar_window_ratio": a.ar_window_MeV.hi / a.ar_window_MeV.lo,
        "ar_evolution_ratio_10_to_100_MeV": t_ar_ev[1] / t_ar_ev[0],
        "ar_composed_amplified_t": ar_amp,
        "ar_composed_postselected_t_max": ar_post,
        "ar_postselected_campaign_t_max": ar_post * a.ar_repeats.lo,
        "ar_evolution_envelope_fraction": (t_ar_ev[0] / env, t_ar_ev[1] / env),
        "eps_l_ar_unit_fault": eps_ar_unit,
        "eps_l_ar_from_fault_budget": eps_ar_budget,
        "eps_l_ar_postselected_budget": fb / ar_post,
        "eps_l_ar_postselected_stated": a.eps_l_ar_postselected.lo,
        "eps_l_ar_stated": (a.eps_l_ar.lo, a.eps_l_ar.hi),
        "eps_l_ar_lo_implied_depth": fb / a.eps_l_ar.lo,
        "lq_ar_envelope_fraction": (a.lq_ar.lo / a.rfi_2033_lq.lo, a.lq_ar.hi / a.rfi_2033_lq.lo),
    }
    return Result(era="2033", lq=(reg_4he, reg_4he), hard_ops=(t_lo, t_hi),
                  breakdown=(gs, walk_t, walk_rot, rest), intermediates=inter,
                  shots=(shots_camp, shots_camp), wall_time_s=wall, epsilon_l=eps_l,
                  notes=("R11b (ch02-route-ii-iii-stated): route (ii) is the 1911.06368 lattice pionless Hamiltonian "
                         "and walk; lambda = M x 706.76 MeV (19.08 / 45.23 GeV) and the four registers are derived from "
                         "its cited inputs. Route (iii) register and evolution cost are from 2507.22814 (evolution only), "
                         "Tab. 5 reproduced from its Thm. 5 and repriced at 7 T per Toffoli (R12). "
                         "The composed T figures stay Stated: app01:135 says the insertion/amplification schedule is the "
                         "authors' composition and is not itemized. Ledger R-069.",
                         "5^3 box (R11b fallback): printed x10-11 at W_q=14 (33.9 MeV); V^2 gives 10.30-10.68, one W_q "
                         "bit x one LCU bit gives 10.78-11.18; both fall in it.",
                         "On-ramp (R11b fallback): 'a few x 1e7 T', not separately costed, carried as (2e7, 5e7); eps_l "
                         "printed '<~ 3e-9 at 3e7 T' (0.1/3e7 = 3.33e-9).",
                         "Walk (R11, ch02-no-compiled-walk (i)): COMPILED lower bound from the 1911.06368 appendix, "
                         "4095 x (15477 T + 2035 R_y x 26.24 T) = 2.82e8 T, 29-42% of 6.7-9.6e8, bracketing the "
                         "stated ~35%. 12C at W_q=13: 8191 x 1.427e5 = 1.17e9, 42-43% of 2.7-2.8e9 (not printed).",
                         "Registers (R11, ch02-ancilla-4x (a)): 141 = 108 + W_q 12 + 2N_A-1 21; 130 = 108+1+21; "
                         "292 = 256+13+23; 280 = 256+1+23. The 2028 ~130 ancilla stay budgeted, not derived.",
                         "Headline = primary 4He instance; 12C, on-ramp, 40Ar are in INSTANCE_ROWS.",
                         "Breakdown at the lo end of the stated range: gs prep at the stated <=1.5%, the compiled "
                         "walk (two COMPILED rows), and the remainder insertion+amplification+readout (UNSOURCED).",
                         "The 34.7 MeV 12C bin is 2 pi x 19.1 GeV x 64/27 / 2^13 = 34.73 MeV at the stated W_q=13.",
                         "Every eps_l row is 0.1 expected faults per shot (R3): 4He 1.0-1.5e-10 (printed ~1e-10), "
                         "12C 3.6e-11, on-ramp <~3e-9 at 3e7 T, 40Ar 1.1e-11-6.7e-10 amplified / 6.3e-11 post-selected, 2028 6.7e-7.",
                         "Stated T does not respond to lattice size or W_q; the +1-bit response is emitted "
                         "from the stated factor and from the walk share separately."))


def _model_codesign(a: Assumptions) -> Result:
    cyc = a.t_gate_s.lo
    anchor = a.rfi_2033_t.lo
    n_lo, n_hi = a.n_orb_conv.lo, a.n_orb_conv.hi
    nt = a.n_trotter_codesign.lo
    tau_max = _tau_max_from_sigma(a.sigma_dune_MeV.lo)       # 100 GeV^-1, from sigma
    c_lo, c_hi = _route_i(a, n_lo, nt), _route_i(a, n_hi, nt)   # each end its own eps (R-TOL)
    s_lo, ts_lo, t_lo = c_lo["strings"], c_lo["t_step"], c_lo["t_shot"]   # 2.53e8, 8.47e9, 1.694e12
    s_hi, ts_hi, t_hi = c_hi["strings"], c_hi["t_step"], c_hi["t_shot"]   # 8e8, 2.75e10, 5.507e12
    # the 2/sigma window doubles N_Trotter: a different circuit, so its own eps
    w = a.window_trunc_factor.lo
    t2_lo, t2_hi = _route_i(a, n_lo, nt * w)["t_shot"], _route_i(a, n_hi, nt * w)["t_shot"]   # 3.45e12, 1.12e13
    eps_l_cd = (a.fault_budget.lo / t_hi, a.fault_budget.lo / t_lo)   # 1.82e-14, 5.90e-14 -> '2-6e-14'
    cx = _route_i(a, a.n_orb_explore.lo, nt)                 # 1.6e7 rotations, eps 2.5e-5, 26.78 T
    s_x, ts_x, t_x = cx["strings"], cx["t_step"], cx["t_shot"]   # 8e4, 2.14e6, 4.28e8
    lq_sys = (4 * n_lo, 4 * n_hi)                            # 600-800
    lq = (lq_sys[0] + a.n_anc.lo, lq_sys[1] + a.n_anc.lo)    # 800-1000
    lq_16o = (4 * a.n_orb_16o.lo + a.n_anc.lo, 4 * a.n_orb_16o.hi + a.n_anc.lo)
    # Eq. (Nq) first-quantized line (R12): log2 of the single-particle state count per nucleon, 4 N_orb states
    # (2507.22814 :112, 'log2(Omega) qubits ... d log2 M position + 2 spin-isospin'), derived here
    fq_sys_12c = (_fq_register(a.A_c.lo, 4 * n_lo), _fq_register(a.A_c.lo, 4 * n_hi))          # 120, 120
    fq_sys_16o = (_fq_register(16, 4 * a.n_orb_16o.lo), _fq_register(16, 4 * a.n_orb_16o.hi))  # 160, 176
    fq_min = fq_sys_12c[0] + a.n_anc.lo
    t0 = a.shot_overhead_s.lo
    horizon_s = a.horizon_yr.lo * YEAR_S
    wall = (t_lo * cyc + t0, t_hi * cyc + t0)                # one machine, serial (E27)
    # R11b ch02-codesign-shots (c): Nyquist tau grid for |q0| <= q0_max out to tau_max
    tau_step = math.pi / a.q0_max_GeV.lo                     # 3.1416 GeV^-1
    tau_points = math.ceil(tau_max / tau_step)               # ceil(31.83) = 32
    n_bin = (a.configs_momentum_transfers.lo * tau_points, a.configs_momentum_transfers.hi * tau_points)   # 160, 320
    shots = (_shots(a, n_bin[0]), _shots(a, n_bin[1]))       # 3.2e5, 6.4e5
    # E27: one machine. Serial campaign = shots x per-shot time; the reduction needed to fit the horizon is
    # serial years / 5, an algorithmic/co-design reduction in shots x per-shot cost, never a machine count
    serial_yr = (shots[0] * wall[0] / YEAR_S, shots[1] * wall[1] / YEAR_S)             # 1.718e4, 1.117e5
    reduction = (serial_yr[0] / a.horizon_yr.lo, serial_yr[1] / a.horizon_yr.lo)       # 3436, 22338
    shots_in_horizon = (horizon_s / wall[1], horizon_s / wall[0])                      # 28.65, 93.13
    # r25 (factory.json): rotations on disjoint quadruples run at the same time, at most n_anc = 200 per layer
    # (one parity ancilla each), so factories beyond ~200 sit idle. Floor = shots x (T / 200) x 10 us.
    cd_depth = (t_lo / a.n_anc.lo, t_hi / a.n_anc.lo)
    cd_floor_yr = (shots[0] * cd_depth[0] * REACTION_TIME_S / YEAR_S,
                   shots[1] * cd_depth[1] * REACTION_TIME_S / YEAR_S)                  # 859, 5585
    shots_per_config = a.eps_target.lo ** -2 * tau_points    # 12800, printed '~1.3e4'
    q_scale = a.qubitized_precision_MeV.lo / a.rescale_bin_MeV.lo
    q_tof = (a.qubitized_toffoli.lo * q_scale, a.qubitized_toffoli.hi * q_scale)
    ar_i = (_route_i(a, a.n_orb_ar_route_i.lo, nt)["t_shot"], _route_i(a, a.n_orb_ar_route_i.hi, nt)["t_shot"])
    ar_i_2s = (_route_i(a, a.n_orb_ar_route_i.lo, nt * w)["t_shot"],
               _route_i(a, a.n_orb_ar_route_i.hi, nt * w)["t_shot"])
    # requirements summary (app01:277-278)
    configs = (a.n_comp.lo * a.configs_momentum_transfers.lo * a.n_nuclei.lo,
               a.n_comp.lo * a.configs_momentum_transfers.hi * a.n_nuclei.lo)
    rot = Primitive("rot_synth", c_lo["n_rot"], c_lo["t_rot"], CircuitStatus.SCALING,
                    "arxiv_2312_05344;" + SYNTHESIS_SRC,
                    f"N_orb={n_lo}: N_orb^4/2 strings x N_Trotter={nt}; R-TOL eps_rot={c_lo['eps_rot']:.4g}, "
                    f"RUS fit {c_lo['t_rot']:.3f} T; lo end of 150-200")
    inter = {
        "rotations_per_shot": (c_lo["n_rot"], c_hi["n_rot"]),
        "eps_rot": (c_lo["eps_rot"], c_hi["eps_rot"]),
        "t_per_rotation": (c_lo["t_rot"], c_hi["t_rot"]),
        "explore_rotations_per_shot": cx["n_rot"],
        "explore_eps_rot": cx["eps_rot"],
        "explore_t_per_rotation": cx["t_rot"],
        "explore_strings_per_step": s_x,
        "explore_t_per_step": ts_x,
        "explore_t_per_shot": t_x,
        "explore_inside_anchor": t_x <= anchor,
        "explore_lq": 4 * a.n_orb_explore.lo + a.n_anc.lo,
        "n_trotter": nt,
        "tau_max_from_sigma_GeVinv": tau_max,
        "tau_max_GeVinv_stated": a.tau_max_GeVinv.lo,
        "trotter_step_GeVinv": tau_max / nt,
        "strings_per_step": (s_lo, s_hi),
        "t_per_step": (ts_lo, ts_hi),
        "t_per_shot": (t_lo, t_hi),
        "t_per_shot_2sigma_window": (t2_lo, t2_hi),
        "over_anchor_factor": (t_lo / anchor, t_hi / anchor),
        "campaign_hard_ops_aggregated": (shots[0] * t_lo, shots[1] * t_hi),   # app01:172 '~5e17-4e18'
        "lq_system_jw": lq_sys,
        "lq_total": lq,
        "lq_total_full_anc_bracket": (lq_sys[0] + a.n_anc_bracket.lo, lq_sys[1] + a.n_anc_bracket.hi),
        "lq_16o_system": (4 * a.n_orb_16o.lo, 4 * a.n_orb_16o.hi),
        "lq_16o_total": lq_16o,
        "lq_fq_minimal": fq_min,
        "lq_fq_system_12c": fq_sys_12c,
        "lq_fq_system_16o": fq_sys_16o,
        "fq_compression_12c": (4 * n_lo / fq_sys_12c[0], 4 * n_hi / fq_sys_12c[1]),          # 5.0, 6.7
        "fq_compression_16o": (4 * a.n_orb_16o.lo / fq_sys_16o[0], 4 * a.n_orb_16o.hi / fq_sys_16o[1]),   # 5.0, 6.8
        "wall_per_shot_yr": (wall[0] / YEAR_S, wall[1] / YEAR_S),
        "wall_per_shot_days": (wall[0] / 86400, wall[1] / 86400),
        "shots": shots,
        "shots_stated": (a.shots_codesign.lo, a.shots_codesign.hi),
        "tau_step_GeVinv": tau_step,
        "tau_points": tau_points,
        "n_bin": n_bin,
        "shots_per_component": a.eps_target.lo ** -2,
        "shots_per_config_from_eq_nshot": shots_per_config,
        "shots_per_config_stated": a.shots_per_config_stated.lo,
        "configurations": configs,
        "configurations_stated": a.configs_stated.lo,
        "campaign_shots_from_requirements": (configs[0] * shots_per_config, configs[1] * shots_per_config),
        "campaign_shots_per_nucleus_from_requirements": (configs[0] / a.n_nuclei.lo * shots_per_config,
                                                         configs[1] / a.n_nuclei.lo * shots_per_config),
        "campaign_shots_stated_product": a.configs_stated.lo * a.shots_per_config_stated.lo,
        "shot_overhead_s": t0,
        "horizon_yr": a.horizon_yr.lo,
        "campaign_serial_yr": serial_yr,
        "campaign_serial_yr_stated": (a.codesign_serial_yr.lo, a.codesign_serial_yr.hi),
        "cost_reduction_for_horizon": reduction,
        "cost_reduction_for_horizon_stated": (a.codesign_reduction.lo, a.codesign_reduction.hi),
        "shots_in_horizon": shots_in_horizon,
        "shots_in_horizon_whole": (math.floor(shots_in_horizon[0]), math.floor(shots_in_horizon[1])),   # 28, 93
        "codesign_t_depth_per_shot": cd_depth,
        "codesign_f_star": (t_lo / cd_depth[0], t_hi / cd_depth[1]),
        "codesign_floor_wall_yr": cd_floor_yr,
        "codesign_floor_over_horizon": (cd_floor_yr[0] / a.horizon_yr.lo, cd_floor_yr[1] / a.horizon_yr.lo),
        "shots_in_horizon_stated": (a.codesign_shots_in_horizon.lo, a.codesign_shots_in_horizon.hi),
        "expected_faults_1e12_at_rfi_2028_floor": 1e12 * a.eps_l_2028_floor.lo,
        "eps_l_unit_fault": (a.fault_budget_unit.lo / t_hi, a.fault_budget_unit.lo / t_lo),
        # app01:185 prints '2-7e-14' at the R3 0.1-fault budget and '~1e12 ops at ~1e-13/op' as the round shorthand
        "eps_l_from_fault_budget": eps_l_cd,
        "eps_l_codesign_stated": (a.eps_l_codesign.lo, a.eps_l_codesign.hi),
        "eps_l_1e12_shorthand": a.fault_budget.lo / 1e12,
        "eps_l_1e12_shorthand_stated": a.eps_l_codesign_shorthand.lo,
        "qubitized_rescaled_toffoli": q_tof,
        # chapter names no Toffoli->T convention (app01:171 compares Toffoli to T directly)
        "qubitized_rescaled_t_textbook": (toffoli_t(q_tof[0], "textbook"), toffoli_t(q_tof[1], "textbook")),
        "qubitized_rescaled_t_jones": (toffoli_t(q_tof[0], "jones"), toffoli_t(q_tof[1], "jones")),
        "qubitized_over_trotter_toffoli_as_t": (q_tof[0] / t_lo, q_tof[1] / t_hi),
        "qubitized_over_trotter_textbook": (toffoli_t(q_tof[0], "textbook") / t_lo,
                                            toffoli_t(q_tof[1], "textbook") / t_hi),
        "ar_route_i_t_per_shot": ar_i,
        "ar_route_i_t_per_shot_2sigma": ar_i_2s,
        "ar_route_i_t_stated": (a.t_ar_route_i.lo, a.t_ar_route_i.hi),
        "ar_route_i_top_end_ratio_stated_over_derived": a.t_ar_route_i.hi / ar_i_2s[1],
        "ar_route_i_low_end_ratio_stated_over_derived": a.t_ar_route_i.lo / ar_i[0],
        "ar_route_i_lq_system": (4 * a.n_orb_ar_route_i.lo, 4 * a.n_orb_ar_route_i.hi),
    }
    return Result(era="codesign", lq=lq, hard_ops=(t_lo, t_hi),
                  breakdown=(rot,), intermediates=inter,
                  shots=shots, wall_time_s=wall,
                  epsilon_l=eps_l_cd,
                  notes=("Derived (R-TOL): (N_orb^4/2) x 200 steps = 5.06e10-1.6e11 rotations at N_orb=150-200, "
                         "eps_rot 4.44e-7-2.5e-7, 33.47-34.42 T each: 1.694e12-5.507e12 T.",
                         "tau_max = 1/sigma = 100 GeV^-1 is derived from sigma_dune; N_Trotter=200 is Stated "
                         "(step 0.5 GeV^-1, no error bound shown), so hard_ops scale with N_Trotter, not with tau_max.",
                         "LQ: 4 N_orb + N_anc with N_anc=200 (implied by the 16O row) gives 800-1000; the box, the "
                         "requirements row, app01:166 and Table 1.1 print '~800-1000' (R11, ch02-codesign-lq (i)).",
                         "Shots (R11b, ch02-codesign-shots (c)): Nyquist tau step pi/(1 GeV) = 3.14 GeV^-1, 32 points "
                         "to 100 GeV^-1; N_bin = (5-10) x 32 = 160-320; 400 x 5 x N_bin = 3.2-6.4e5 (printed '~3-6e5'). "
                         "Aggregate 5.42e17-3.52e18 ('~5e17-4e18'). One machine (E27): 1.72e4-1.12e5 yr serial, "
                         "28-93 whole shots in 5 yr, cost reduction to fit 3.4e3-2.2e4 ('~3e3-2e4x'). Requirements: 400 x 32 = 12800 per configuration ('~1.3e4').",
                         "Route-(i) 40Ar (R-TOL, each circuit its own eps): chapter's own formula gives "
                         "1.37e13-2.34e14 (2.79e13-4.75e14 at 2/sigma); the prose prints 1e13-5e14.",
                         "app01:185 eps_l: 0.1 fault at 1.7-5.5e12 T = 1.82e-14-5.90e-14, printed '2-6e-14' / 'of order "
                         "1e-13-1e-14'; the '~1e12 ops at ~1e-13/op' shorthand is 0.1/1e12 (R3 applied in the FIX pass; "
                         "the pre-FIX '1e-12-1e-13' was the unit-fault reading).",
                         f"Toffoli->T convention not named by the chapter; both textbook ({T_PER_TOFFOLI['textbook']}) "
                         f"and jones ({T_PER_TOFFOLI['jones']}) equivalents emitted."))


def model(a: Assumptions, era: str) -> Result:
    """Return the Result for '2028' | '2033' | 'codesign'."""
    if era == "2028":
        return _model_2028(a)
    if era == "2033":
        return _model_2033(a)
    if era == "codesign":
        return _model_codesign(a)
    raise ValueError(f"era must be one of {ERAS}, got {era!r}")


PUBLISHED = {
    "2028": Published(lq=(170, 170), hard_ops=(1.004e5, 1.004e5),
                      src=f"{TEX}:2028 box ('~170', '~1e5 T-gates ... 1.004e5, 0.4% over'); R-TOL exact product "
                          "5e3 x 20.086 = 1.00428e5 (was 5e3 x 30 = 1.5e5)",
                      rel_tol=0.05),
    "2033": Published(lq=(141, 141), hard_ops=(6.7e8, 9.6e8),
                      src=f"{TEX}:2033 box lines 225-226, primary 4He instance ('141 peak', '6.7-9.6e8 T'); "
                          "other instances via INSTANCE_ROWS", rel_tol=0.05),
    "codesign": Published(lq=(800, 1000), hard_ops=(1.7e12, 5.5e12),
                          src=f"{TEX}:co-design box ('~800-1000 LQ', '~1.7-5.5e12 T/shot'; R-TOL, was 1.5-4.8e12). "
                              "4 N_orb + N_anc = 600-800 + 200; printed as the range under R11 "
                              "ch02-codesign-lq (i) (was '~1000')", rel_tol=0.10),
}


def INSTANCE_ROWS(a: Assumptions, era: str, r: Result):
    """Rows for resources.json: [(label, (lq_lo, lq_hi), (t_lo, t_hi), {extra}), ...]."""
    i = r.intermediates
    if era == "2028":
        return [
            ("4He chiral NN+3N, N_Trotter=1", r.lq, r.hard_ops, {"route": "i", "nucleus": "4He"}),
            ("4He pionless, N_Trotter=3", r.lq,
             (i["t_per_shot_pionless"], i["t_per_shot_pionless"]), {"route": "i", "nucleus": "4He"}),
        ]
    if era == "2033":
        return [
            ("4He inclusive, 3^3 site basis, W_q=12 (primary)", (a.lq_4he.lo, a.lq_4he.lo),
             (a.t_4he.lo, a.t_4he.hi), {"route": "ii", "nucleus": "4He"}),
            ("4He tau-grid correlator (companion)", (a.tau_grid_lq.lo, a.tau_grid_lq.lo),
             i["tau_grid_t"], {"route": "ii", "nucleus": "4He"}),
            ("12C inclusive, 4^3 site basis (near miss)", (a.lq_12c.lo, a.lq_12c.lo),
             (a.t_12c.lo, a.t_12c.hi), {"route": "ii", "nucleus": "12C"}),
            ("12C moments / sum rules (on-ramp)", (a.lq_onramp.lo, a.lq_onramp.lo),
             (a.t_onramp.lo, a.t_onramp.hi), {"route": "ii", "nucleus": "12C", "t_note": "not separately costed"}),
            ("40Ar inclusive window, evolution only", i["lq_ar_evolution_with_spin_isospin"],
             i["ar_evolution_t"], {"route": "iii", "nucleus": "40Ar",
                                   "t_note": "arxiv_2507_22814 GQSP repriced at 7 T/Toffoli (derived here)"}),
            ("40Ar composed, amplified prep (conditional)", (i["lq_ar_evolution_with_spin_isospin"][0], a.lq_ar.hi),
             i["ar_composed_amplified_t"], {"route": "iii", "nucleus": "40Ar",
                                            "t_note": "Stated; set on the source's 4-T evolution (lower bound)"}),
            ("40Ar composed, post-selected prep (conditional)", (a.lq_ar.hi, a.lq_ar.hi),
             (i["ar_evolution_t"][0], i["ar_composed_postselected_t_max"]), {"route": "iii", "nucleus": "40Ar",
                 "t_note": "upper end Stated; set on the source's 4-T evolution (lower bound)"}),
        ]
    if era == "codesign":
        return [
            # 16O ('heavier still', 1000-1400 LQ) has no stated T count; it stays in
            # intermediates['lq_16o_total'] rather than as a row with an empty T.
            ("12C converged basis, N_orb 150-200", r.lq, r.hard_ops, {"route": "i", "nucleus": "12C"}),
        ]
    return []
