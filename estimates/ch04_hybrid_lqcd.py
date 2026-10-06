"""Ch. 4 — Hybrid lattice QCD: trace-log and trace-inverse on a quantum computer.
Reproduces the LQ and hard-op numbers of applications/app12_hybrid_lqcd.tex.

Two circuits are priced, and they are not the same circuit (ruling R7, H. Lamm, 2026-09-28).

2028: the FREE-FERMION block encoding of arXiv:2407.13080 Sec. II.2 at V = 4^4. W = M^dag M has
2D+1 = 9 nonzero entries per column, a 4-qubit column register, one rotation flag, and
B = D_s O_c O_A D_s with O_A two y-rotations that do not depend on the site (Table I: O(D log D log V)).
The paper gives the circuit structure and NO T-counts and no adder circuit, so every T figure here is
DERIVED: Toffoli ladders for the controlled shifts and the projector phases, synthesized rotations
through common.t_per_rotation. The breakdown is primitive-level. The QET target is the paper's
shifted one, log|x|/log(s/lambda_min) + 1/2 (Sec. IV), at the tolerance that guarantees 1% on
log det W (ruling 'SHIFT + 1%', H. Lamm, 2026-09-28).

2033 and the 1e4-LQ rung: the GAUGED block encoding of M (Tr M^-1), priced at the query level as D V log2 V
T-gates with a unit prefactor (ruling R3, 2026-10-02: M has 2D+1 = 9 entries per column, so the D^2 of the
33-entry M^dag M does not apply; Eq. Ngate_hybrid, app12:89). The gauged 4^4 log det comparison of the 2028
section block-encodes W = M^dag M and keeps D^2 V log2 V. No Toffoli or rotation count exists for either.

INPUTS AND SOURCES
  D = 4                          app12:48 (4D Euclidean lattice)                              Stated
  L = 4, V = 4^4 = 256           app12:110, ruling R7 (2^4 periodic is degenerate)            Assumed
  m0 = 0.4                       app12:110                                                     Stated
  K = 1                          arxiv_2407_13080 txt:266-267 "typically set to one"          Cited
  9 nonzeros per column          arxiv_2407_13080 txt:323-324                                 Cited
  column register 4 qubits       arxiv_2407_13080 txt:103-105, Fig. 1                         Cited
  1 block-encoding flag          arxiv_2407_13080 Figs. 2-3                                   Cited
  W/s, s >= 2(m0^2+2K^2), 4K^2   arxiv_2407_13080 Eqs. 25-26 (txt:327-340)                    Cited
  QME: copy of |n>, m QPE, 1 a   arxiv_2407_13080 Sec. III (txt:583-594), Eq. 42              Cited
  QET signal 1, workspace 3, Hadamard ancilla 1      ancilla decomposition, this work         Assumed
  O_c: 4-8 shifts, 3-4 controls, clean-ancilla ladders (common.mcx_toffoli_count)             Assumed
  O_A: 2-4 synthesized rotations at generic angles (Fig. 2: two singly-controlled R_y)        Assumed
       at s = 2(m0^2+2K^2), theta_0 = 0 (Eq. 25): 2 at both ends (Clifford+T-exact dropped, R-TOL)  derived
  projector phase: 8 Toffolis + 2 rotations, d+1 of them                                      Assumed
  shift of the target by 1/2     arxiv_2407_13080 Sec. IV, txt:697-703                        Cited
  1% on log det W                polynomial guarantee and sampling target, app12:116-117      Assumed
  d = 84                         LP-minimal degree, this work (not in the paper)              Assumed
  Toffoli = 7 T                  main-overview:47, ruling R5                                  Cited
  eps_syn = 0.1 x tol = 1.784e-4 ruling SYN-BUDGET (report default EPS_SYN = 1e-2 not used)  Stated
                                 eps_rot = sqrt(eps_syn/N_rot) (common.eps_rot_for) -> 7.266e-4
                                 at N_rot = 338 (both ends), full RUS fit -> 21.19 T
  eps = 1.784e-3 (derived from the 1%), rms           app12:117 -> shots = Var/eps^2, Var = 1    Stated
  Hadamard-test per-shot variance 1 - mu^2 <= 1; the bound 1 is used (ruling R16, rms convention of Chs. 3, 8)
  T cap 1e5 (2028), 1e9 (2033); eps_l 1e-8 (2028)   DOE_RFI_2026                              Cited
  1 us per T gate                app12:128 (Ch. overview convention)                          Cited
  shot overhead t0 ~ 0.1 ms      main-overview:48 (register init, readout, decode)            Assumed
  campaign horizon 5 yr          app12 Requirements summary                                    Stated
  ONE machine, serial            ruling E27 (H. Lamm, 2026-10-02): no wall time divides by a machine count
  p_f = 0.1 expected faults/shot app12:180; ruling R3                                         Stated
  LCU 3-10x (QSP-degree lever dropped, R11b: d at LP floor)   app12:137                     Stated
  even-odd sublattice 2x         arxiv_2407_13080 Sec. V (txt:755-760)                        Cited
  V_2033 = 8^3x16, SU(3), N_f = 2+1   app12:144-145                                           Stated
  kappa_2033 = 1e2, am = 0.04        app12 1000-LQ box; kappa = (4+am)/am (shot audit, R3)      Stated
  column index of M: 4 qubits (2D+1 = 9 entries), D V log2 V per query   ruling R3            Stated
  readout: block-encoding success probability P = f^2 am c, f = 1/2, c = (1/V_F) Re Tr M^-1     shot audit, R3
    c = 0.0275 (free 8^3x16, antiperiodic t, computed here) to 0.08 (interacting, audit estimate)   Assumed
  first result: 1 configuration at 30% on P (R1); campaign: 3 configurations at 20% (R3 '1-3')   Stated
  Banks-Casher demonstration: nu(Lambda)/V_F at a = 0.15 fm, Lambda = 0.10-0.15, d = 640/430, 5 cfg,
    Sigma^(1/3) = 270 MeV, 4 tastes (simplobs.json C2; ruling R7)                              Stated/Assumed
  d_inv = kappa ln(1/eps_rel), eps_rel = e^-5 = 0.674% of the peak of 1/x: 500 at kappa = 1e2,
    5e3 at kappa = 1e3. LP-minimal odd degree, this work (apply_log/ch04_r11_lp_odd.py)       Assumed
  n_QPE = 12, ancilla = 50            app12:152                                                Stated
  eps_2033 = 1e-2 rms per config     the retired Hadamard-test target; kept for the 1e4-LQ shot plan (R16)  Stated
  N_cfg = 1e2                         retired at 2033 (R3); informational                     Stated
  V_codesign = 24^3x48, kappa = 1e3, per-shot reduction >= 2.6x    app12 1e4-LQ paragraph (R3 reprice)  Stated
  synth_to_poly_2033 = 0.1: eps_syn = 0.1 x 0.5 e^-5 = 3.369e-4 at the 1000-LQ and 1e4-LQ rungs (E19 P1)  Stated
  V_production = 48^3x96              app12:33,180                                             Stated
  N_cfg_production = 1e3              app12:196                                                Stated
WHAT IS NOT DERIVED HERE
  d = 84 (and 48, 112, 314, 36, 510): outputs of a linear program that needs numpy/scipy and
    cannot be rerun inside this standard-library package. They are carried Assumed, tied to the
    lambda_min/s and the tolerance they were computed for (the model flags a mismatch), and
    reproduced by apply_log/ch04_roundB_lp_degree.py (pinned by a test that skips without scipy).
  The ancilla decomposition (signal, workspace, Hadamard ancilla) and the Toffoli ladders: this
    work's compilation of the structure the paper draws. The paper specifies no adder circuit.
  scaling_prefactor = 1 T per unit of D^2 V log2 V for the GAUGED circuit: a working prefactor,
    nowhere justified (Uncited); every gauged BE-query Primitive is SCALING.
  d_inv = 500 and 5e3: the LP law d = kappa ln(1/eps_rel) (R11); the LP itself needs scipy and is
    reproduced by apply_log/ch04_r11_lp_odd.py (pinned at kappa = 30 by a test that skips without
    scipy). At kappa = 1e3 the LP runs out of memory near d ~ 5e3, so 5e3 is the law extrapolated.
  n_QPE and ancilla counts of the gauged rungs (12 + 50): asserted (Stated).
  The interacting condensate c = 0.08 (shot audit estimate) and the Banks-Casher P (Sigma^(1/3) = 270 MeV,
    4 tastes, uncertain by ~2x in a 1.2 fm box): not derived; they set the shot bands.
  The 2033 and 1e4-LQ T-depths: structural inferences in factory.json from a unit-prefactor SCALING count,
    carried as the ratio F* = 1-4 (D_T = N_T / F*).
WALL TIME AND FACTORIES
  t_gate_s = 1e-6 s per T, shot_overhead_s = 1e-4 s per shot, ONE machine, serial (E27). Factory model of
  common.depth_exports: one T per 10 us per factory, one 10 us reaction per T-depth layer; the report's 1 us per T
  is the 10-factory baseline. Corrected shot time = max(N_T t_gate_s, D_T x 10 us) + t0 (R10: F* < 10 everywhere).
  2028 (one tier, factory.json 'scheduled'): D_T = d (2 t_rot) + (d+1)(3 + 2 t_rot) = 7.42e3 at both ends (the O_A and
    phase rotations are one sequential chain; O_c runs underneath). F* = 2.6-4.8 < 10: 74 ms per shot, 6.5 h
    (was 19-36 ms, 1.7-3.1 h). wall_first_result_s = None; wall_campaign_s = 6.5 h.
  2033 (factory.json '1000-LQ' runs, ratio carried to D V log2 V): D_T = N_T/4 (unary-iteration link load) to N_T
    (multiplexed rotation on the flag), F* = 1-4. Corrected shot 532-2130 s (serial baseline 213 s).
    First result: 1 configuration, 30% on P -> 1.39e4-4.04e4 shots -> wall_first_result_s 86 d - 2.7 yr.
    Campaign: 3 configurations, 20% on P -> 9.37e4-2.73e5 shots -> wall_campaign_s 1.6-18.4 yr (3.7x the horizon
    at the top). Workspace lever (factory.json, unpriced): 3 parallel unary-iteration subtrees on ~100 idle LQ give
    F* ~ 12 and restore the baseline 0.63-1.84 yr.
  1e4-LQ rung (prose, factory.json): F* = 1-4 after the 2.6x reduction; 2.5e5-9.9e5 s per shot; 1e7 calls
    7.8e4-3.1e5 yr. No exports (not a 2028/2033 era).
APPLIED 2026-09-28 (apply_log/ch04.md, apply_log/ch04_roundB.md)
  R3: p_f = 0.1 expected faults per shot; the 2028 box prints eps_l <~ 2.5e-6 (0.1 / 4.07e4).
  R4: carry exact, round the printed result once.
  R7: 2028 repriced as the free-fermion circuit at V = 4^4: 18 LQ (Hadamard route) / 34 (QME
    route). Rung label ~30-LQ -> ~20-LQ. The 1e4-LQ rung prints its 92-LQ register and has a
    PUBLISHED co-design entry. Prose adopts log2 V.
  SHIFT + 1%: the paper's 1/2 shift and the tolerance 1.784e-3: d 112 -> 84, 2.0e4-4.1e4 T per
    shot (was 2.7e4-5.4e4), 1.4e6 shots (was 4.6e4), 8-16 h. Synthesis error in the report's
    convention (T4, common.eps_per_rotation): N eps^2.
APPLIED 2026-09-29 (apply_log/rtol_ch04.md)
  R-TOL: the fixed 1e-4 per rotation is replaced by eps_rot = sqrt(1e-2 / N_rot), N_rot the synthesized
    rotations of ONE shot of that circuit. The two ends of the 2028 band are two circuits (multiplexed O_A,
    merged shifts: 338 rotations; literal: 506), so each carries its own eps: 5.439e-3 and 4.446e-3, 17.851
    and 18.185 T per rotation (full RUS fit, unchanged formula). Per shot 1.785e4-3.748e4 (was 2.009e4-
    4.067e4); the unshifted d = 314 comparison is its own circuit (1258 / 1886 rotations). The synthesis
    error per shot is now 1e-2 at both ends, 5.6x the polynomial tolerance 1.784e-3 (was 3.4-5.1e-6).
    2033 and codesign: SCALING BE queries with no rotation split, not affected; the informational
    qsp_phase_eps (Uncited 1e-4) is outside the per-shot T and is left as it was.
  R-TOL FIX (verifier, same day): R-TOL does not count Clifford+T-exact rotations. At the chapter's
    s = 2(m0^2+2K^2) = 4.32, Eq. 25 gives theta_0 = 2 arccos(1) = 0, so the controlled R_y(theta_0) of O_A is
    the identity: 0 rotations, 0 T. The literal end drops from 4 to 2 synthesized rotations per query
    (the multiplexed end stays 2, angles +-theta_1/2). Both ends: N_rot = 338, eps_rot = 5.439e-3, 17.851 T.
    Per shot 1.785e4-3.431e4 (was 1.785e4-3.748e4); margin 2.91-5.60x; eps_l 2.914e-6; unshifted d = 314
    6.78e4-1.29e5 (1258 rotations at both ends); 2.58e10-4.96e10 T, 7.17-13.79 h.
  SYN-BUDGET (H. Lamm, 2026-09-29; apply_log/ch04_synbudget.md): EPS_SYN = 1e-2 stays the report-wide
    default, but the 2028 circuit sets eps_syn = 0.1 x the polynomial tolerance 1.784e-3 = 1.784e-4 per shot
    (input synth_to_poly_2028 replaces eps_syn_2028). eps_rot = sqrt(1.784e-4/338) = 7.266e-4, 21.19 T per
    rotation. Per shot 1.898e4-3.544e4; margin 2.82-5.27x; eps_l 2.821e-6; unshifted d = 314 7.20e4-1.34e5
    (eps_rot 3.766e-4); 2.745e10-5.126e10 T, 7.62-14.24 h. Synthesis on log det W: 0.1506 = 0.1%; combined
    worst case 1% + 1% + 0.1% = 2.1% (printed ~2%).
APPLIED 2026-09-30 (apply_log/r11_ch04.md)
  ch04-dinv-500-1e4 (B'): d_inv = kappa ln(1/eps_poly), the LP minimum at eps_rel = e^-5 = 0.674% of the
    peak of 1/x: 500 at kappa = 1e2 (2033 box unchanged, 8.52e8 T), 5e3 at kappa = 1e3 (was 1e4). The
    co-design instance 24^3x48 is 5e3 x 2.0533e8 = 1.027e12 T (was 2.05e12), PUBLISHED (1e12, 1e12).
  ch04-48x96-reduction-factor (a, B' figures): the deliverable is a per-shot T reduction of ~10x
    (1.027e12 / 1e11 = 10.3x; was 'cut d_inv by 1e2x'). Levers: LCU 3-10x x even-odd 2x = 6-20x; the
    degree is at the LP floor. 48^3x96 is 5e3 x 3.9647e9 = 1.98e13 T and needs 198x (~2e2x).
  ch04-qme-vs-incoherent-1000lq: prose and Eq. Nshot now say Hadamard-test sampling, 1/eps^2, at every
    rung; QME is the upgrade path. The model already priced it; the QME count stays as an alternative.
APPLIED 2026-10-01 (apply_log/r16_ch4.md)
  R16 rms shot counting (H. Lamm): eps is an rms (one-standard-deviation) additive error, as in Chs. 3 and 8.
    Hadamard-test shots = Var/eps^2 with per-shot variance Var = 1 - mu^2 <= 1; the bound 1 is used. This
    replaces ln(1/delta)/eps^2 at delta = 1e-2 (a factor ln 100 = 4.605 fewer shots). QME: 1/eps calls
    (constant one, the Ch. 8 convention). delta_2028 and delta_2033 removed; input shot_variance added.
    2028: 1.4464e6 -> 3.1408e5 shots; 2.745e10-5.126e10 -> 5.961e9-1.113e10 T; 7.62-14.24 h -> 1.656-3.092 h.
    2033: 4.605e4 -> 1e4 shots/cfg; 4.605e6 -> 1e6 shots; 3.92e15 -> 8.52e14 T; 3.92e9 -> 8.52e8 s (27.0 yr).
    Codesign: 4.605e7 -> 1e7 subroutine calls. A 99% two-sided interval at the same eps needs z^2 = 6.63x.
APPLIED 2026-10-01 (apply_log/r17_ch04.md, ruling E19 a = P1)
  2028 SYN-BUDGET split (E9) kept. The same 0.1x rule now sets the 1000-LQ and 1e4-LQ synthesis budget:
    eps_syn = 0.1 x 0.5 e^-5 = 3.369e-4 on the normalized trace (0.067% of the peak of 1/x); the report default
    1e-2 would be 2.97x the polynomial error. +2.81 T per rotation, absorbed in the unit prefactor; no T moves.
  24^3x48 deliverable: per-shot reduction ~10x -> at least 10.3x (1.0266e12 / 1e11 = 10.27x; leaves 9.967e10,
    fits); LCU at least 5.2x with the even-odd 2x (5.13x needed). 48^3x96 after it: 1.925e12 (not printed).
APPLIED 2026-10-02 (apply_log/r23_ch04.md, ruling E27 "no multiple machines", H. Lamm)
  Every wall time is serial on ONE machine at shots x (T x 1 us + t0), t0 = 0.1 ms (input shot_overhead_s)
    (superseded by r25: depth-bound walls).
    The model never had a machine count; the retired statement is the prose 'wall-clock time falls linearly
    with device count and duty cycle' (app12, after the 1000-LQ box).
  2028: 19.1-35.5 ms per shot (printed '19--36 ms', was '19--35 ms' without t0); 1.664-3.101 h ('1.7--3.1 h').
  2033: 852 s per shot, 1e6 shots = 8.52e8 s = 27.0 yr serial = 5.40x the 5-year horizon. Five years hold
    1.852e5 shots = 18 whole configurations at 1e4 shots (ensemble mean 2.4e-3); one configuration is 98.6 days.
    The full 100 configurations need a 5.4x reduction in shots x per-shot T (algorithmic, not a machine count).
  1e4-LQ rung: after the 10.3x per-shot reduction a shot is 9.967e10 T = 9.97e4 s (27.7 h); the 1e7 calls
    are 9.97e11 s = 3.16e4 yr = 6.3e3x the horizon; five years hold 1.58e3 shots, 16% of one configuration.
    QME (1/eps calls) supplies 1e2x of that at 1e2x the coherent depth. No T or LQ moves.
APPLIED 2026-10-02 (apply_log/r25_ch04.md, rulings R1, R3, R7, R9, R10, H. Lamm)
  R3: the Tr M^-1 rungs block-encode M (9 entries, 4-qubit column index) at D V log2 V per query. 2033: 85 -> 83 LQ,
    8.52e8 -> 2.13e8 T. 1e4-LQ rung: 92 -> 90 LQ, 1.03e12 -> 2.57e11 T, reduction 10.3x -> 2.6x; 48^3x96
    1.98e13 -> 4.96e12 T, ~2e2x -> ~50x. Readout at 2033: success probability P = f^2 am c (am = 0.04) replaces the
    Hadamard test; the eps = 1e-2 / 1e2-configuration plan is retired (it did not resolve mu = 5.5e-4 - 1.6e-3).
  R1: first result 1 configuration at 30%; campaign 3 configurations at 20%. R7: Banks-Casher mode number at
    a = 0.15 fm as the 2033 demonstration row (1.8-2.7e8 T at D V log2 V).
  R9/R10: depth_exports on every 2028/2033 era; corrected (depth-bound) walls everywhere F* < 10.
"""

from __future__ import annotations

import itertools
import math
from dataclasses import dataclass

from estimates.common import (Assumed, Cited, Stated, Uncited, Tagged, Primitive,
                              CircuitStatus, Result, Published, ERAS, SYNTHESIS_SRC,
                              T_PER_TOFFOLI, mcx_toffoli_count, EPS_SYN, EPS_SYN_SRC, eps_rot_for,
                              t_per_rotation, eps_per_rotation, toffoli_t, depth_exports,
                              REACTION_TIME_S, FACTORY_BASELINE)

TEX = "app12"
SYNBUDGET_SRC = ("ruling SYN-BUDGET (H. Lamm, 2026-09-29): EPS_SYN = 1e-2 stays the report-wide default; a circuit "
                 "may set a tighter budget when its accuracy target needs one. Ch. 4 2028: eps_syn = 0.1 x the "
                 "polynomial tolerance 1.784e-3 (app12:125)")
PAPER = "arxiv_2407_13080"
FJ = "factory.json (T-depth audit, 2026-10-02) ch04"
SECONDS_PER_YEAR = 365.25 * 86400.0      # 3.15576e7 s (as ch07_curved_space.py)


@dataclass(frozen=True)
class Assumptions:
    # EVERY field is a Tagged value. Ranges are (lo, hi) tuples.
    # --- lattice / formulation -------------------------------------------
    D: Tagged = Stated(4, "app12:48", "D-dimensional Euclidean lattice; 4D (app12:200)")
    n_s: Tagged = Stated(1, "app12:48", "n_s=1 for staggered fermions")
    log_base: Tagged = Stated(2, "app12:55,89", "'we evaluate the gauged cost as D^2 V log_2 V' (R7); "
                                                "1.7e6 at V=8192 = 16*8192*13")
    scaling_prefactor: Tagged = Uncited(
        1.0, "GAUGED circuit only: app12:55 'a working prefactor of one T per unit' of D^2 V log2 V; "
             "no compiled circuit and no citation for the prefactor")
    toffoli_convention: Tagged = Cited("textbook", "main-overview:47 (ruling R5, H. Lamm 2026-09-28)",
                                       "7 T per Toffoli; app12:125 '7 T per Toffoli'")
    # --- 2028 tech demo: free staggered W at V = 4^4 (ruling R7) -----------
    L_2028: Tagged = Assumed(4, "linear extent, periodic; V = L^D = 4^4 = 256 (app12:110). Ruling R7: the old "
                                "V = 2^4 instance is degenerate (M = m0 x identity)", "app12:26,110")
    m0_2028: Tagged = Stated(0.4, "app12:110", "bare mass m_0 = 0.4")
    K_2028: Tagged = Cited(1.0, PAPER, "txt:266-267 'The K coupling is typically set to one'; app12:110 'K=1'")
    nnz_free_2028: Tagged = Cited(9, PAPER, "Sec. II.2, txt:323-324 'only nine nonzero entries in each "
                                            "row/column' = 2D+1; app12:111")
    column_qubits_2028: Tagged = Cited(4, PAPER, "Sec. II.1, txt:103-105 'l = 0-7 ... 8 <= l < 16', Fig. 1 draws "
                                                 "l_0..l_3; reused for free fermions (txt:317-319). The paper's "
                                                 "'six qubits' (txt:349) is the U(1)-gauged matrix, 33 nonzeros")
    be_flag_qubits_2028: Tagged = Cited(1, PAPER, "Figs. 2-3: the |0> qubit the O_A y-rotations act on")
    qet_signal_qubits_2028: Tagged = Assumed(1, "QET signal qubit carrying the projector-controlled phases; the "
                                                "standard construction, not drawn in the paper")
    mcx_workspace_2028: Tagged = Assumed(3, "clean workspace for the Toffoli ladders: the 5-control projector AND "
                                            "needs 5-2 = 3; the 4-control shifts of O_c need 2 of the same 3")
    hadamard_qubits_2028: Tagged = Assumed(1, "Hadamard-test ancilla. The Hadamard test is this chapter's choice; "
                                              "the paper's trace step is QME (Sec. III)")
    n_qpe_2028: Tagged = Assumed(8, "m = 8 QPE register on the QME upgrade path; the paper's example runs "
                                    "m = 6, 7, 8, 9 (txt:554)")
    qme_sign_qubits_2028: Tagged = Cited(1, PAPER, "txt:592-594 'QME ... requires one additional qubit, which we "
                                                   "shall label a'")
    poly_shift_2028: Tagged = Cited(0.5, PAPER, "Sec. IV, txt:697-703 'it helps to shift the entire polynomial away "
                                                "from y = -1 ... A shift of 1/2 works well. This known shift can "
                                                "always be subtracted at the end of the calculation'; app12:100,113")
    target_rel_2028: Tagged = Assumed(0.01, "relative error on log det W: guaranteed for the polynomial by the uniform "
                                            "tolerance, and resolved by the sampling as an rms error "
                                            "(app12:116-117; ruling 'SHIFT + 1%'; rms since R16)")
    d_log_2028: Tagged = Assumed(84, "LP-minimal degree, this work: smallest even d with uniform error <= 1.784e-3 on "
                                     "the shifted target log|x|/log(s/lambda_min) + 1/2 over [lambda_min/s, 1] and "
                                     "|p| <= 1 on [-1, 1], at lambda_min/s = 0.0370. Not in the paper")
    d_log_lam_norm: Tagged = Assumed(0.16 / 4.32, "the lambda_min/s the LP degrees here were computed for; the "
                                                  "model flags d as stale if the instance moves off it")
    d_log_err: Tagged = Assumed(1.6993e-3, "LP optimum at d = 84 with the shift: the uniform error reached")
    d_log_err_below: Tagged = Assumed(1.8660e-3, "LP optimum at d - 2 = 82 with the shift: above the tolerance, so 84 "
                                                 "is minimal")
    d_shifted_tol_1e2: Tagged = Assumed(48, "LP-minimal degree with the shift at tolerance 1e-2 (err(48) = 9.939e-3, "
                                            "err(46) = 1.104e-2); app12:100")
    d_unshifted_tol_1e2: Tagged = Assumed(112, "LP-minimal degree WITHOUT the shift at tolerance 1e-2 (err(112) = "
                                               "9.907e-3, err(110) = 1.026e-2): the round-B box before this ruling; "
                                               "app12:100")
    d_unshifted_2028: Tagged = Assumed(314, "LP-minimal degree WITHOUT the shift at the tolerance 1.784e-3 "
                                            "(err(314) = 1.7704e-3, err(312) = 1.7945e-3): what the shift saves; "
                                            "app12:100,136")
    d_lp_at_lam016: Tagged = Assumed(36, "record: LP-minimal degree at lambda_min = 0.16, no shift, tolerance 1e-2 "
                                         "(err 8.99e-3; err(34) = 1.011e-2): what the old 16+12 = 28 box was sized for")
    d_lp_at_tol_1e3: Tagged = Assumed(510, "record: LP-minimal degree at lambda_min/s = 0.0370, no shift, tolerance "
                                           "1e-3: err(510) = 9.985e-4, err(508) = 1.0061e-3")
    d_fig12_paper: Tagged = Cited((64, 70), PAPER, "Fig. 12 caption, txt:644-647: total degree d = 64, 66, 68, 70 "
                                                   "with d_r = 60, d = d_f + d_r (txt:686); the paper states NO "
                                                   "lambda_min for it")
    oc_shift_gates_2028: Tagged = Assumed((4, 8), "controlled shifts in O_c: 8 as Fig. 1 draws them (add2 for +mu, "
                                                  "sub2 for -mu); 4 when the +2mu and -2mu shifts, which coincide on "
                                                  "L = 4, are merged")
    oc_controls_2028: Tagged = Assumed((3, 4), "controls per shift: the 4 column qubits as drawn; 3 when the sign "
                                               "bit is dropped by the merge")
    oa_rotations_2028: Tagged = Assumed((2, 4), "O_A at GENERIC angles: Fig. 2 draws two singly-controlled R_y. 4 "
                                                "synthesized rotations if each controlled R_y is 2 rotations + 2 "
                                                "CNOT; 2 if the pair is compiled as one multiplexed rotation. The "
                                                "model drops rotations whose angle is Clifford+T exact (R-TOL): at "
                                                "s = 2(m0^2+2K^2) theta_0 = 0, so both ends are 2 (oa_synthesized_rotations)")
    projector_toffolis_2028: Tagged = Assumed(8, "one controlled projector phase: the AND of the 5 block-encoding "
                                                 "qubits (4 column + 1 flag) is computed onto the signal qubit by a "
                                                 "4-Toffoli ladder and uncomputed by 4 more")
    and_tree_depth_2028: Tagged = Cited(3, FJ + " 2028 run", "the 5-control AND onto the signal qubit is a 3-layer "
                                                        "tree on the 3 workspace qubits; measurement uncompute is "
                                                        "Clifford (T-depth 0)")
    projector_rotations_2028: Tagged = Assumed(2, "the phase rotation controlled on the Hadamard-test ancilla: 2 "
                                                  "synthesized rotations + 2 CNOT. For even d the block-encoding "
                                                  "queries need no control")
    synth_to_poly_2028: Tagged = Stated(0.1, SYNBUDGET_SRC,
                                        "synthesis budget of the 2028 circuit as a fraction of the polynomial "
                                        "tolerance: eps_syn = 0.1 x 1.784e-3 = 1.784e-4 per shot, tighter than the "
                                        "report-wide EPS_SYN = 1e-2 because the 1% guarantee on log det W rests on "
                                        "the polynomial tolerance; eps_rot = sqrt(eps_syn / N_rot) from the "
                                        "circuit's own rotation count (common.eps_rot_for). Replaces eps_syn = 1e-2")
    shot_variance: Tagged = Stated(1.0, "ruling R16 (H. Lamm, 2026-10-01); app12:92",
                                   "per-shot variance of the +-1 Hadamard-test outcome, 1 - mu^2 <= 1; the bound 1 "
                                   "is used. Shots = Var/eps^2 at rms additive eps (the convention of Chs. 3 and 8). "
                                   "Replaces confidence delta = 1e-2 and ln(1/delta)/eps^2")
    eps_l_2028: Tagged = Cited(1e-8, "DOE_RFI_2026", "app12:26: the RFI's eps_l=1e-8. The box prints the "
                                                     "requirement 0.1/T instead (R3)")
    t_cap_2028: Tagged = Cited(1e5, "DOE_RFI_2026", "app12:26,126: <=1e5 T per circuit")
    t_gate_s: Tagged = Cited(1e-6, "app12:128", "1 us per T-gate; convention in Ch. overview")
    shot_overhead_s: Tagged = Assumed(1e-4, "per-shot overhead t0 ~ 0.1 ms (register init, final readout, decode), "
                                            "report rule, main-overview:48; as Ch. 8 shot_overhead_s. Wall = shots x "
                                            "(T x 1 us + t0), serial on ONE machine (ruling E27, 2026-10-02)")
    campaign_horizon_yr: Tagged = Stated(5, "app12:Requirements", "'Campaign horizon 5 years'. Serial time on one "
                                            "machine is compared with it; no machine count enters (ruling E27)")
    # --- the gauged V=4^4 comparison (app12:28,94,136) ---------------------
    lcu_saving: Tagged = Stated((3, 10), "app12:137", "LCU-based BE compression '~3-10x saving'")
    # R11b (ch04-lever-ii, author option 2): the former qsp_degree_saving = Stated((2, 3)) lever
    # ('~2-3x lower polynomial degree') is DROPPED. d = 84 is already the LP minimum at the stated
    # uniform tolerance for the free instance, so no QSP family lowers it; the gauged instance carries
    # d = 84 as a working figure. The levers are LCU x sublattice only: 6-20x against a 27.5x gap.
    sublattice_saving: Tagged = Cited(2, PAPER, "Sec. V, txt:755-760 'only half the matrix need be block-encoded, "
                                                "since the even and odd sub-lattices are completely independent'")
    # --- 1000-LQ rung -----------------------------------------------------
    V_2033: Tagged = Stated(8192, "app12:144", "V=8^3x16")
    VF_box_2033: Tagged = Stated(2.5e4, "app12:144", "'V_F = N_c V ~= 2.5e4 per flavor' as printed (was 1.6e5; ch04-VF-8x8x8x16)")
    Nc_2033: Tagged = Stated(3, "app12:145", "SU(3) staggered")
    Nf_2033: Tagged = Stated(3, "app12:145", "N_f=2+1")
    kappa_2033: Tagged = Stated(1e2, "app12:146", "kappa ~ 1e2 (am ~ 0.04, m_pi ~ 600 MeV, near the strange mass; was m_pi ~ 300-400 MeV)")
    d_inv_eps_rel: Tagged = Assumed(math.exp(-5), "uniform error of the 1/x polynomial relative to the peak of 1/x "
                                                  "on [1/kappa, 1]: e^-5 = 0.674% (app12:147,169; ruling R11 "
                                                  "ch04-dinv-500-1e4 B')")
    d_inv_2033: Tagged = Assumed(500, "LP-minimal odd degree, this work: kappa ln(1/eps_rel) = 1e2 x 5 "
                                      "(apply_log/ch04_r11_lp_odd.py: err(501) = 3.334e-3 = 0.5 e^-5.01); "
                                      "app12:147 'd_inv = 500, the smallest degree ... (linear program, this work)'")
    neighbor_box_2033: Tagged = Stated(4, "app12 1000-LQ box", "'13 site + 4 column index + 4 color' (R3, 2026-10-02: "
                                                        "the box block-encodes M, 2D+1 = 9 entries; was 6, the "
                                                        "M^dag M count)")
    n_qpe_2033: Tagged = Stated(12, "app12:152", "12 QPE")
    n_anc_2033: Tagged = Stated(50, "app12:152", "~50 ancilla (BE workspace + phase loader + QME)")
    eps_2033: Tagged = Stated(1e-2, "app12 1e4-LQ paragraph", "eps=1e-2 rms per configuration (R16): the Hadamard-test "
                                                           "target the 1000-LQ box used before R3; kept as the 1e4-LQ "
                                                           "rung's inherited shot plan (1e4 per configuration)")
    n_cfg_2033: Tagged = Stated(100, "retired 2033 plan (R3)", "~1e2 configurations; informational since R3")
    # --- 1000-LQ rung after rulings R1, R3, R7, R10 (H. Lamm, 2026-10-02) -----------------------------
    be_D_power_M: Tagged = Stated(1, "ruling R3 (H. Lamm, 2026-10-02); shot audit ch04 needs_author 4",
                                  "the Tr M^-1 rungs block-encode M, 2D+1 = 9 entries per column: D V log2 V per "
                                  "query. D^2 V log2 V (33 entries of M^dag M) stays for the gauged 4^4 log det")
    L_s_2033: Tagged = Stated(8, "app12 1000-LQ box", "V = 8^3 x 16, spatial extent")
    L_t_2033: Tagged = Stated(16, "app12 1000-LQ box", "V = 8^3 x 16, temporal extent (antiperiodic for the free c)")
    am_2033: Tagged = Stated(0.04, "ruling R3; shot audit ch04 change (1)",
                             "bare mass am ~ 0.04: kappa = (4 + am)/am ~ 1e2 at subnormalization 4 + am (am = 4/99)")
    condensate_interacting_2033: Tagged = Assumed(0.08, "shot audit ch04: interacting estimate of c = (1/V_F) Re Tr M^-1 "
                                                        "at am = 0.04, 'up to 0.08'; the free value is computed here")
    peak_norm_2033: Tagged = Assumed(0.5, "QET polynomial normalized to 1/2 at sigma_min (the 0.5 e^-5 tolerance of "
                                          "the uncertainty row); P = f^2 am c (shot audit readout lever 1)")
    rel_first_2033: Tagged = Stated(0.30, "ruling R1 (H. Lamm, 2026-10-02)", "first result at 30% statistical error")
    n_cfg_first_2033: Tagged = Stated(1, "ruling R1; shot audit ch04 minimal_first_result", "one configuration")
    rel_campaign_2033: Tagged = Stated(0.20, "ruling R3", "campaign: 20% relative against the classical value")
    n_cfg_campaign_2033: Tagged = Stated(3, "ruling R3 '1-3 configurations'", "three priced, the top of 1-3")
    parallelism_2033: Tagged = Cited((1, 4), FJ + " 1000-LQ runs",
                                     "F* = N_T / D_T: 1 for a multiplexed rotation on the flag qubit (every T on the "
                                     "critical path), 4 for the unary-iteration (QROM) link load. Ratio-preserving: "
                                     "D_T = N_T / F*, carried to the D V log2 V count")
    subtree_lq: Tagged = Cited(34, FJ + " 1000-LQ campaign", "~2 log2 L LQ per parallel unary-iteration subtree, "
                                                            "F* ~ 4k for k subtrees; unpriced co-design lever")
    bc_a_fm: Tagged = Stated(0.15, "simplobs.json ch04 C2; ruling R7", "coarse a = 0.15 fm variant (MILC HISQ class)")
    bc_hbarc_gev_fm: Tagged = Cited(0.1973, "PDG", "hbar c = 0.1973 GeV fm")
    bc_sigma_third_gev: Tagged = Assumed(0.270, "simplobs.json ch04 C2: Sigma^(1/3) = 270 MeV; P uncertain by ~2x in "
                                                "the small box")
    bc_tastes: Tagged = Stated(4, "simplobs.json ch04 C2", "4 staggered tastes")
    bc_lambda_lo: Tagged = Stated(0.10, "simplobs.json ch04 C2", "Lambda = 0.10 (lattice units), d = 640")
    bc_lambda_hi: Tagged = Stated(0.15, "simplobs.json ch04 C2", "Lambda = 0.15 (lattice units, 197 MeV), d = 430")
    bc_d_lambda_lo: Tagged = Assumed(640, "simplobs.json ch04 C2: step-filter LP degree at Lambda = 0.10")
    bc_d_lambda_hi: Tagged = Assumed(430, "simplobs.json ch04 C2: step-filter LP degree at Lambda = 0.15")
    bc_n_cfg: Tagged = Stated(5, "simplobs.json ch04 C2", "5 configurations; gauge variance of nu ~ nu per config")
    eps_l_2033: Tagged = Stated(1e-10, "app12:154", "eps_l ~ 1e-10")
    t_cap_2033: Tagged = Cited(1e9, "DOE_RFI_2026", "app12:156: ~1e9 hard-op envelope")
    V_16x4: Tagged = Stated(65536, "app12:22", "'88 LQ at V=16^4'")
    V_128x4: Tagged = Stated(2 ** 28, "app12:22", "'100 LQ at V=128^4'")
    # --- 1e4-LQ (codesign) rung ------------------------------------------
    p_f: Tagged = Stated(0.1, "app12:180", "<=0.1 expected fault per shot (app12:94: p_f/eps_l); "
                                             "report-wide convention by ruling R3 (2026-09-28)")
    eps_l_codesign: Tagged = Stated(1e-12, "app12:180", "eps_l ~ 1e-12")
    V_codesign: Tagged = Stated(24 ** 3 * 48, "app12:180", "V ~ 24^3x48")
    kappa_codesign: Tagged = Stated(1e3, "app12:180", "kappa ~ 1e3 at physical mass")
    d_inv_codesign: Tagged = Assumed(5e3, "kappa ln(1/eps_rel) = 1e3 x 5 at the 2033 tolerance: the LP law, "
                                          "extrapolated (the LP runs out of memory near d ~ 5e3); app12:33,182 "
                                          "'d_inv ~= 5e3' (R11; was 1e4)")
    per_shot_reduction: Tagged = Stated(2.6, "app12 1e4-LQ paragraph", "'cut the per-shot T by at least 2.6x' (R3 reprice "
                                        "at D V log2 V, 2026-10-02: 2.567e11 / 1e11; was 10.3x at D^2 V log2 V, E19 P1)")
    parallelism_codesign: Tagged = Cited((1, 4), FJ + " 1e4-LQ run", "F* = 1-4 on the 92-LQ register; LCU keeps the "
                                                                 "PREP+SELECT unary-iteration form, so F* does not change")
    synth_to_poly_2033: Tagged = Stated(0.1, SYNBUDGET_SRC + " extended to the 1000-LQ and 1e4-LQ rungs (author, E19 P1)",
                                        "total synthesis error = 0.1 x the polynomial error 0.5 e^-5 = 3.369e-3 on the "
                                        "normalized trace (peak of 1/x = 1/2), i.e. 3.369e-4 = 0.067% of the peak; "
                                        "app12:155,162. Informational: the SCALING count has a unit prefactor, so no "
                                        "hard_ops change")
    V_production: Tagged = Stated(48 ** 3 * 96, "app12:180", "physical-mass 48^3x96")
    n_cfg_production: Tagged = Stated(1e3, "app12:196", "'x 1e3 production cfg'")
    subroutine_calls_box_1e4: Tagged = Stated(1e7, "app12:197",
                                              "'~1e7' (1e4 shots/cfg x 1e3 production cfg); R16 rms, was 4.6e7")
    hvp_site_variance: Tagged = Assumed(1.0, "referee H2 (2026-10-04): gauge variance per site of the local one-link "
                                             "current term <x|Gamma M^-1|x>, order one; not measured (NEEDS_AUTHOR, "
                                             "George). Sets the noise-matched disconnected-HVP shot count")
    n_qpe_codesign: Tagged = Stated(12, "app12:180", "'12 QPE', carried from the 1000-LQ box (app12:152)")
    n_anc_codesign: Tagged = Stated(50, "app12:180", "'~50 ancilla', carried from the 1000-LQ box (app12:152)")
    # --- informational only: QSP phase loader synthesis ------------------
    qsp_phase_eps: Tagged = Uncited(
        1e-4, "no synthesis tolerance stated for the phase loader (app12:100); "
              "1e-4 is the value this informational figure carried before R-TOL (there is no longer a fixed "
              "report-wide tolerance; R-TOL would set it from the loader's own rotation count); NOT in the "
              "chapter's per-shot T and not printed")

    def __post_init__(self):
        for name in ("D", "n_s", "L_2028", "V_2033", "Nc_2033", "Nf_2033",
                     "V_16x4", "V_128x4", "V_codesign", "V_production",
                     "n_qpe_2028", "n_qpe_2033", "n_anc_2033",
                     "n_qpe_codesign", "n_anc_codesign", "neighbor_box_2033",
                     "nnz_free_2028", "column_qubits_2028", "be_flag_qubits_2028",
                     "qet_signal_qubits_2028", "mcx_workspace_2028", "hadamard_qubits_2028",
                     "qme_sign_qubits_2028", "d_log_2028", "d_lp_at_lam016", "d_lp_at_tol_1e3",
                     "d_shifted_tol_1e2", "d_unshifted_tol_1e2", "d_unshifted_2028",
                     "projector_toffolis_2028", "projector_rotations_2028",
                     "sublattice_saving", "be_D_power_M", "L_s_2033", "L_t_2033", "n_cfg_first_2033",
                     "n_cfg_campaign_2033", "subtree_lq", "bc_tastes", "bc_d_lambda_lo", "bc_d_lambda_hi", "bc_n_cfg",
                     "and_tree_depth_2028"):
            v = getattr(self, name)
            if not v.is_range and (v.lo < 1 or int(v.lo) != v.lo):
                raise ValueError(f"{name} must be a positive integer, got {v.value}")
        L = int(self.L_2028.lo)
        if L < 2 or L & (L - 1):
            raise ValueError(f"L_2028 must be a power of two >= 2 (binary coordinate registers), got {L}")
        if int(self.d_log_2028.lo) % 2:
            raise ValueError("d_log_2028 must be even: the normalized log|x| is an even function")
        if 2 ** int(self.column_qubits_2028.lo) < int(self.nnz_free_2028.lo):
            raise ValueError("column register cannot index the nonzero entries")
        for name in ("shot_variance", "eps_2033", "p_f",
                     "eps_l_2028", "eps_l_2033", "eps_l_codesign", "target_rel_2028", "synth_to_poly_2028",
                     "synth_to_poly_2033", "am_2033", "condensate_interacting_2033", "peak_norm_2033",
                     "rel_first_2033", "rel_campaign_2033", "bc_lambda_lo", "bc_lambda_hi",
                     "d_log_lam_norm", "d_log_err", "d_log_err_below", "poly_shift_2028", "d_inv_eps_rel",
                     "hvp_site_variance"):
            v = getattr(self, name)
            if not (0 < v.lo <= 1):
                raise ValueError(f"{name} must lie in (0, 1], got {v.value}")
        for name in ("t_cap_2028", "t_cap_2033", "t_gate_s",
                     "m0_2028", "K_2028", "kappa_2033", "d_inv_2033", "kappa_codesign",
                     "d_inv_codesign", "per_shot_reduction", "scaling_prefactor", "n_cfg_2033",
                     "n_cfg_production", "lcu_saving", "campaign_horizon_yr",
                     "oc_shift_gates_2028", "oc_controls_2028", "oa_rotations_2028"):
            if getattr(self, name).lo <= 0:
                raise ValueError(f"{name} must be positive")
        if self.shot_overhead_s.lo < 0:
            raise ValueError("shot_overhead_s must be non-negative")
        for name in ("lcu_saving", "oc_shift_gates_2028", "oc_controls_2028",
                     "oa_rotations_2028", "d_fig12_paper"):
            v = getattr(self, name)
            if v.is_range and v.lo > v.hi:
                raise ValueError(f"{name} range must be (lo, hi)")
        for name in ("parallelism_2033", "parallelism_codesign"):
            v = getattr(self, name)
            if v.lo < 1 or (v.is_range and v.lo > v.hi):
                raise ValueError(f"{name} must be a (lo, hi) range of factor counts >= 1")
        if int(self.L_s_2033.lo) ** 3 * int(self.L_t_2033.lo) != int(self.V_2033.lo):
            raise ValueError("V_2033 must equal L_s^3 x L_t")
        if self.bc_lambda_lo.lo >= self.bc_lambda_hi.lo or self.bc_d_lambda_lo.lo < self.bc_d_lambda_hi.lo:
            raise ValueError("Banks-Casher: Lambda_lo < Lambda_hi, and the lower cut needs the higher degree")
        for name in ("bc_a_fm", "bc_hbarc_gev_fm", "bc_sigma_third_gev"):
            if getattr(self, name).lo <= 0:
                raise ValueError(f"{name} must be positive")
        if self.d_log_err.lo > self.d_log_err_below.lo:
            raise ValueError("the LP error cannot grow with the degree")
        if self.log_base.lo not in (2, math.e, 10):
            raise ValueError("log_base must be 2, e or 10")
        if self.toffoli_convention.value not in T_PER_TOFFOLI:
            raise ValueError(f"unknown Toffoli convention {self.toffoli_convention.value!r}")


# --------------------------------------------------------------------------- #
# Eq. (Nq_hybrid) and Eq. (Ngate_hybrid), app12:88-89
# --------------------------------------------------------------------------- #

def ceil_log2(n: int) -> int:
    """ceil(log2 n) exactly for integer n >= 1."""
    return (int(n) - 1).bit_length()


def site_qubits(V: int) -> int:
    return ceil_log2(V)


def neighbor_qubits(D: int) -> int:
    """Column (displacement) index of the GAUGED M^dagger M: ceil(log2(2D^2+1)).

    With links, M^dagger M connects displacements {0, +-2e_mu, +-e_mu +- e_nu (mu != nu)}:
    1 + 2D + 4 C(D,2) = 2D^2 + 1 = 33 at D = 4, six qubits (arXiv:2407.13080 txt:347-349)."""
    return ceil_log2(2 * D * D + 1)


def neighbor_qubits_free(D: int) -> int:
    """Column index of the FREE staggered W: ceil(log2(2D+1)) = 4 at D = 4.

    Without links the mixed mu != nu terms cancel (arXiv:2407.13080 Eqs. 23-24), leaving the
    diagonal and the 2D two-hop neighbours: 2D + 1 = 9 nonzero entries."""
    return ceil_log2(2 * D + 1)


def color_qubits(Nc: int) -> int:
    return ceil_log2(Nc * Nc)


def register(V: int, D: int, Nc: int, n_qpe: int, n_anc: int) -> int:
    return site_qubits(V) + neighbor_qubits(D) + color_qubits(Nc) + n_qpe + n_anc


def scaling_units(D: int, V: int, base: float) -> float:
    """D^2 V log V, the gauged block-encoding cost of Eq. Ngate_hybrid as the boxes evaluate it."""
    return D * D * V * math.log(V, base)


def column_qubits_M(D: int) -> int:
    """Column index of the gauged (or free) staggered M itself: 2D + 1 = 9 nonzero entries per column (the
    diagonal and the 2D nearest neighbours), ceil(log2 9) = 4 qubits. Ruling R3: the Tr M^-1 rungs block-encode M."""
    return ceil_log2(2 * D + 1)


def register_M(V: int, D: int, Nc: int, n_qpe: int, n_anc: int) -> int:
    """Eq. Nq with the column index of M (the Tr M^-1 rungs, R3)."""
    return site_qubits(V) + column_qubits_M(D) + color_qubits(Nc) + n_qpe + n_anc


def scaling_units_M(D: int, V: int, base: float, power: int = 1) -> float:
    """D^power V log V: the gauged block encoding of M at power 1 (ruling R3, shot audit: 9 entries, not the 33 of
    M^dag M); power 2 is scaling_units."""
    return D ** power * V * math.log(V, base)


def free_condensate(L: int, Lt: int, m: float) -> float:
    """c = (1/V_F) Re Tr M^-1 for free staggered fermions on L^3 x Lt, periodic in space, antiperiodic in time:
    (1/V) sum_p m / (m^2 + sum_mu sin^2 p_mu), K = 1 (the 2028 normalization). 0.0275 at 8^3 x 16, am = 0.04."""
    s = [math.sin(2.0 * math.pi * n / L) ** 2 for n in range(L)]
    st = [math.sin(math.pi * (2 * n + 1) / Lt) ** 2 for n in range(Lt)]
    tot = 0.0
    for x, y, z in itertools.product(s, repeat=3):
        q = x + y + z
        tot += sum(m / (m * m + q + w) for w in st)
    return tot / (L ** 3 * Lt)


def shots_success_probability(P: float, rel: float) -> float:
    """Shots for relative standard error rel on a success probability P: (1 - P)/(P rel^2) (binomial)."""
    return (1.0 - P) / (P * rel * rel)


def corrected_shot_s(t: float, depth: float, t_gate_s: float, t0: float) -> float:
    """One shot on ONE machine when the factories cannot all be fed (R10): max(N_T t_gate_s, D_T t_r) + t0."""
    return max(t * t_gate_s, depth * REACTION_TIME_S) + t0


def scaling_units_logDV(D: int, V: int, base: float) -> float:
    """The variant arXiv:2407.13080 Table I writes: D^2 V log(DV). Informational."""
    return D * D * V * math.log(D * V, base)


def shots_incoherent(eps: float, variance: float) -> float:
    """Var/eps^2: Hadamard-test shots for rms additive error eps (Eq. Nshot_hybrid; R16)."""
    return variance / eps ** 2


def shots_qme(eps: float) -> float:
    """1/eps: QET calls for rms additive error eps under QME, constant one (app12:53,72,92; R16)."""
    return 1.0 / eps


def interval_99_factor() -> float:
    """Shots for a two-sided 99% normal interval of half-width eps over the rms count: z_{0.995}^2 = 6.63."""
    from statistics import NormalDist
    return NormalDist().inv_cdf(0.995) ** 2


def d_inv_formula(kappa: float, eps: float, base: float) -> float:
    """kappa log(kappa/eps): the CKS-class degree (app12:52), kept for comparison only."""
    return kappa * math.log(kappa / eps, base)


def d_inv_lp(kappa: float, eps_rel: float) -> float:
    """kappa ln(1/eps_rel): the LP-minimal odd degree for 1/x on [1/kappa, 1] at uniform error eps_rel
    relative to the peak of 1/x (app12:52,89; apply_log/ch04_r11_lp_odd.py, R11)."""
    return kappa * math.log(1.0 / eps_rel)


def depth_budget(p_f: float, eps_l: float) -> float:
    """Hard ops a shot absorbs at logical-fault probability p_f: p_f/eps_l (app12:94,180)."""
    return p_f / eps_l


def _after_saving(t: float, lo: float, hi: float) -> tuple[float, float]:
    """(lo, hi) T after dividing by a saving factor in [lo, hi]; lo T = largest saving."""
    return (t / hi, t / lo)


# --------------------------------------------------------------------------- #
# Free staggered W = M^dagger M, arXiv:2407.13080 Sec. II.2
# --------------------------------------------------------------------------- #

def free_W_row(site: tuple, L: int, m0: float, K: float) -> dict:
    """One row of W by Eq. (24): (m0^2 + D K^2/2) delta - (K^2/4) sum_mu (delta_{a+2mu} + delta_{a-2mu}).

    D K^2/2 is the paper's 2K^2 at D = 4. Periodic. On L = 4 the +2mu and -2mu sites coincide,
    so the row has D + 1 = 5 entries, not 2D + 1 = 9; on L = 2 every shift returns to the site."""
    D = len(site)
    row = {site: m0 * m0 + D * K * K / 2.0}
    for mu in range(D):
        for hop in (2, -2):
            t = list(site)
            t[mu] = (t[mu] + hop) % L
            t = tuple(t)
            row[t] = row.get(t, 0.0) - K * K / 4.0
    return {k: v for k, v in row.items() if abs(v) > 1e-15}


def free_W_spectrum(L: int, D: int, m0: float, K: float) -> dict:
    """{eigenvalue: multiplicity} of the periodic free W: m0^2 + K^2 sum_mu sin^2(2 pi n_mu / L).

    Plane waves diagonalize Eq. (24): 2K^2 - (K^2/2) sum_mu cos(2 p_mu) = K^2 sum_mu sin^2 p_mu."""
    spec: dict = {}
    s2 = [math.sin(2.0 * math.pi * n / L) ** 2 for n in range(L)]
    for ns in itertools.product(range(L), repeat=D):
        lam = round(m0 * m0 + K * K * sum(s2[n] for n in ns), 12)
        spec[lam] = spec.get(lam, 0) + 1
    return dict(sorted(spec.items()))


def be_subnormalization(D: int, m0: float, K: float, column_qubits: int) -> float:
    """Smallest s for which the O_A angles of Eqs. (25)-(26) exist.

    Half of the 2^c column slots carry the diagonal, so cos(theta_0/2) = 2 (m0^2 + D K^2/2)/s;
    one slot carries each hop, so cos(theta_1/2) = -2^c (K^2/4)/s. At D = 4, c = 4:
    s >= 2(m0^2 + 2K^2) and s >= 4K^2."""
    return max(2.0 * (m0 * m0 + D * K * K / 2.0), 2 ** column_qubits * K * K / 4.0)


def oa_angles(D: int, m0: float, K: float, column_qubits: int, s: float) -> tuple[float, float]:
    """(theta_0, theta_1) of O_A, arXiv:2407.13080 Eqs. 25-26, in the normalization of be_subnormalization:
    cos(theta_0/2) = 2 (m0^2 + D K^2/2)/s, cos(theta_1/2) = -2^c (K^2/4)/s."""
    c0 = min(1.0, 2.0 * (m0 * m0 + D * K * K / 2.0) / s)
    c1 = max(-1.0, -(2 ** column_qubits) * (K * K / 4.0) / s)
    return 2.0 * math.acos(c0), 2.0 * math.acos(c1)


def _exact(phi: float) -> bool:
    """R_y(phi) is exact in Clifford+T when phi is a multiple of pi/4 (conjugate of a power of T)."""
    k = phi / (math.pi / 4.0)
    return math.isclose(k, round(k), abs_tol=1e-9)


def oa_synthesized_rotations(theta0: float, theta1: float) -> tuple[int, int]:
    """Synthesized rotations per O_A query, (multiplexed, literal). R-TOL: Clifford+T-exact rotations are not
    synthesized. Literal: each singly-controlled R_y(theta) = R_y(theta/2) CNOT R_y(-theta/2) CNOT.
    Multiplexed pair: R_y((theta0+theta1)/2) CNOT R_y((theta0-theta1)/2) CNOT."""
    literal = sum(0 if _exact(t / 2.0) else 2 for t in (theta0, theta1))
    multiplexed = sum(0 if _exact(a) else 1 for a in ((theta0 + theta1) / 2.0, (theta0 - theta1) / 2.0))
    return multiplexed, literal


def shift2_toffolis(L: int, n_controls: int) -> int:
    """Toffolis in one controlled 'add2' (or 'sub2') on a log2(L)-bit coordinate register.

    Adding 2 is an increment of the register with its lowest bit removed. On L = 4 that is a
    single X on the high bit, so the controlled shift is one C^n X. On L = 2 it is the identity.
    Larger L needs a controlled incrementer, priced as its cascade of multi-controlled X gates.
    Clean-ancilla ladders (common.mcx_toffoli_count, 'linear'). The paper specifies no adder."""
    bits = ceil_log2(L)
    return sum(mcx_toffoli_count(n_controls + j, "linear") for j in range(bits - 1))


# --------------------------------------------------------------------------- #
# Model
# --------------------------------------------------------------------------- #

def _model_2028(a: Assumptions) -> Result:
    """Free staggered W at V = L^D, the free-fermion block encoding of arXiv:2407.13080 Sec. II.2,
    QET of the normalized log at degree d, Hadamard-test trace."""
    D, L = int(a.D.lo), int(a.L_2028.lo)
    V = L ** D
    m0, K = a.m0_2028.lo, a.K_2028.lo
    conv = a.toffoli_convention.value
    d = int(a.d_log_2028.lo)

    # --- the matrix ---------------------------------------------------------
    spec = free_W_spectrum(L, D, m0, K)
    lam_min, lam_max = min(spec), max(spec)
    logdet = sum(n * math.log(lam) for lam, n in spec.items())             # 150.5528
    row = free_W_row((0,) * D, L, m0, K)
    degenerate = len(row) == 1                                             # L = 2: W = m0^2 x identity
    col = int(a.column_qubits_2028.lo)
    s = be_subnormalization(D, m0, K, col)                                 # 4.32
    lam_norm = lam_min / s                                                 # 0.0370
    log_norm = math.log(1.0 / lam_norm)                                    # 3.2958
    v_log_s = V * math.log(s)                                              # 374.59, added back classically

    # --- register -----------------------------------------------------------
    site = site_qubits(V)
    flag, sig = int(a.be_flag_qubits_2028.lo), int(a.qet_signal_qubits_2028.lo)
    work, had = int(a.mcx_workspace_2028.lo), int(a.hadamard_qubits_2028.lo)
    n_qpe, sign = int(a.n_qpe_2028.lo), int(a.qme_sign_qubits_2028.lo)
    overhead = flag + sig + work + had                                     # 6
    lq = site + col + overhead                                             # 18
    lq_qme = lq - had + site + n_qpe + sign                                # 34

    # --- T per block-encoding query and per shot -----------------------------
    # --- accuracy target (needed first: it sets the synthesis budget) ---------
    rel = a.target_rel_2028.lo                                             # 1%
    tol = rel * abs(logdet) / (V * log_norm)                               # 1.784e-3
    # Synthesis budget (ruling SYN-BUDGET, H. Lamm 2026-09-29): 0.1 x the polynomial tolerance, tighter
    # than the report-wide EPS_SYN = 1e-2, so synthesis stays negligible beside the 1% guarantee.
    eps_syn = a.synth_to_poly_2028.lo * tol                                # 1.784e-4 per shot
    merge = L == 4                                                         # +2mu = -2mu only on L = 4
    g_lo = a.oc_shift_gates_2028.lo if merge else a.oc_shift_gates_2028.hi
    c_lo = a.oc_controls_2028.lo if merge else a.oc_controls_2028.hi
    oc_tof = (g_lo * shift2_toffolis(L, int(c_lo)),
              a.oc_shift_gates_2028.hi * shift2_toffolis(L, int(a.oc_controls_2028.hi)))   # (12, 40)
    theta0, theta1 = oa_angles(D, m0, K, col, s)                          # (0, 5.5086): theta_0 = 0 at this s
    oa_generic = (a.oa_rotations_2028.lo, a.oa_rotations_2028.hi)         # (2, 4) at generic angles
    oa_rot = oa_synthesized_rotations(theta0, theta1)                     # (2, 2): R_y(theta_0) is the identity
    assert all(r <= g for r, g in zip(oa_rot, oa_generic))
    p_tof, p_rot = int(a.projector_toffolis_2028.lo), int(a.projector_rotations_2028.lo)
    n_phase = d + 1

    def circuit(deg: int, end: int) -> dict:
        """One end of the band is one circuit (0: merged shifts + multiplexed O_A; 1: literal). R-TOL:
        its tolerance comes from ITS rotation count, sqrt(eps_syn / N_rot); N_rot does not depend on eps."""
        r = oa_rot[end]
        n = deg * r + (deg + 1) * p_rot                                    # synthesized rotations per shot
        e = eps_rot_for(n, eps_syn)
        t = t_per_rotation(e)                                              # full RUS fit, unchanged
        q = toffoli_t(oc_tof[end], conv) + r * t
        ph = toffoli_t(p_tof, conv) + p_rot * t
        return {"n_rot": n, "eps_rot": e, "t_rot": t, "t_query": q, "t_phase": ph,
                "t_shot": deg * q + (deg + 1) * ph}

    def per_shot(deg: int) -> tuple[float, float]:
        return tuple(circuit(deg, end)["t_shot"] for end in (0, 1))

    ends = (circuit(d, 0), circuit(d, 1))
    n_rot = tuple(c_["n_rot"] for c_ in ends)                              # (338, 338)
    eps_rot = tuple(c_["eps_rot"] for c_ in ends)                          # (7.266e-4, 7.266e-4)
    t_rot = tuple(c_["t_rot"] for c_ in ends)                              # (21.191, 21.191)
    t_query = tuple(c_["t_query"] for c_ in ends)                          # (126.4, 322.4)
    t_phase = tuple(c_["t_phase"] for c_ in ends)                          # (98.38, 98.38)
    t_shot = per_shot(d)                                                   # (1.898e4, 3.544e4)
    cap = a.t_cap_2028.lo
    margin = (cap / t_shot[1], cap / t_shot[0])                            # (2.82, 5.27)
    # total synthesis error in the report's convention (T4): common.eps_per_rotation's randomized
    # model is eps_rot = sqrt(eps_total / N), i.e. eps_total = N eps_rot^2 = eps_syn by construction.
    synth = tuple(n * e ** 2 for n, e in zip(n_rot, eps_rot))              # (1.784e-4, 1.784e-4)
    assert all(math.isclose(eps_per_rotation(s_, n, "incoherent"), e) for s_, n, e in zip(synth, n_rot, eps_rot))
    assert all(math.isclose(s_, eps_syn) for s_ in synth)
    eps_l_req = (a.p_f.lo / t_shot[1], a.p_f.lo / t_shot[0])               # (2.82e-6, 5.27e-6)

    # --- accuracy -------------------------------------------------------------
    shift = a.poly_shift_2028.lo                                           # 1/2
    poly_abs = V * log_norm * tol                                          # 1.5055 = 1% of log det W
    d_ok = (math.isclose(lam_norm, a.d_log_lam_norm.lo, rel_tol=1e-6)
            and a.d_log_err.lo <= tol < a.d_log_err_below.lo)
    eps_s = tol                                                            # sampling resolves the same 1%
    samp_abs = V * log_norm * eps_s
    # synthesis error eps_syn is a shift of the normalized trace, the same units as tol: x V log(s/lambda_min)
    synth_abs = V * log_norm * eps_syn                                     # 0.1506 = 0.1% of log det W
    shift_term = V * log_norm * shift                                      # 421.87, subtracted classically
    d_u = int(a.d_unshifted_2028.lo)
    ends_u = (circuit(d_u, 0), circuit(d_u, 1))                            # its own N_rot: 1258, 1258
    t_shot_u = per_shot(d_u)                                               # (7.20e4, 1.34e5) without the shift

    # --- sampling ---------------------------------------------------------------
    shots = shots_incoherent(eps_s, a.shot_variance.lo)                    # 3.14e5 (R16 rms)
    t0 = a.shot_overhead_s.lo                                              # 0.1 ms per shot (E27: one machine)
    wall_shot = tuple(t * a.t_gate_s.lo + t0 for t in t_shot)              # (0.0191, 0.0355) s
    wall_total = tuple(shots * w for w in wall_shot)                       # (5.99e3, 1.12e4) s = 1.66-3.10 h
    t_total = tuple(shots * x for x in t_shot)                             # (5.96e9, 1.11e10)
    # --- T-depth and factories (R9, R10; factory.json 'scheduled'): the 2 O_A rotations of a query and the controlled
    # phase rotations form one sequential chain on the flag and signal qubits; O_c uses the column register only as
    # controls and runs underneath. Query depth 2 t_rot, phase depth (AND tree) + 2 t_rot, all d + (d+1) sequential.
    and_depth = a.and_tree_depth_2028.lo
    t_depth = d * oa_rot[0] * t_rot[0] + n_phase * (and_depth + p_rot * t_rot[0])   # 7.42e3 at both ends
    dx = depth_exports(t_shot, (t_depth, t_depth), shots, a.t_gate_s.lo, t0)
    wall_shot_corr = tuple(corrected_shot_s(t, t_depth, a.t_gate_s.lo, t0) for t in t_shot)   # (0.0743, 0.0743) s
    wall_corr = tuple(shots * w for w in wall_shot_corr)                   # 2.33e4 s = 6.48 h

    # --- the gauged V = L^D instance: a different circuit (app12:28,94,136) -------
    units = scaling_units(D, V, a.log_base.lo)                             # 32768
    t_gauged = d * a.scaling_prefactor.lo * units                          # 2.75e6
    factor = t_gauged / cap                                                # 27.5
    sub = a.sublattice_saving.lo
    lcu = a.lcu_saving
    lev = (lcu.lo * sub, lcu.hi * sub)                                     # (6, 20): (i) LCU x (iii) sublattice
    after_lcu = _after_saving(t_gauged, lcu.lo, lcu.hi)
    after_levers = _after_saving(t_gauged, lev[0], lev[1])                 # (1.376e5, 4.588e5)

    tof_t = T_PER_TOFFOLI[conv]
    breakdown = (
        Primitive("O_c_controlled_shift_toffolis", d * oc_tof[1], tof_t, CircuitStatus.COMPILED,
                  f"{PAPER} Fig. 1, Sec. II.2 (txt:317-319)",
                  f"literal end: {d} queries x {a.oc_shift_gates_2028.hi} controlled add2/sub2 as Fig. 1 draws them, "
                  f"{int(a.oc_controls_2028.hi)} controls each. ASSUMED here: on L = 4 add2 is one X on the high "
                  f"coordinate bit, each C^4X is a clean-ancilla ladder of 5 Toffolis, 7 T per Toffoli. The paper "
                  f"specifies no adder and no T-count. Merged end: {oc_tof[0]} Toffolis per query"),
        Primitive("O_A_entry_rotations", d * oa_rot[1], t_rot[1], CircuitStatus.COMPILED,
                  f"{PAPER} Fig. 2, Eqs. 25-26; {SYNTHESIS_SRC}",
                  f"literal end: {d} queries x 2 singly-controlled R_y as Fig. 2 draws them. ASSUMED here: each "
                  f"controlled R_y is 2 synthesized rotations + 2 CNOT, randomized synthesis; at s = {s:.4g} "
                  f"theta_0 = {theta0:.4g} (Eq. 25), so R_y(theta_0) is the identity and is not synthesized "
                  f"(R-TOL): {oa_rot[1]} rotations per query, theta_1 = {theta1:.4g}. SYN-BUDGET: eps_rot = "
                  f"sqrt({eps_syn:.4g}/{n_rot[1]}) = {eps_rot[1]:.4g}. Multiplexed end: {oa_rot[0]} rotations per "
                  f"query, {n_rot[0]} per shot, eps_rot = {eps_rot[0]:.4g}"),
        Primitive("projector_phase_toffolis", n_phase * p_tof, tof_t, CircuitStatus.SCALING,
                  "this work; QET construction of Martyn_QSVT_unified_2021, gilyen2018QSingValTransfArXiv",
                  f"{n_phase} controlled projector phases x {p_tof} Toffolis. Counted from the standard "
                  f"projector-controlled phase, which {PAPER} does not draw: compute and uncompute the AND of the "
                  f"{col + flag} block-encoding qubits on {work} clean workspace qubits"),
        Primitive("projector_phase_rotations", n_phase * p_rot, t_rot[1], CircuitStatus.SCALING,
                  f"this work; {SYNTHESIS_SRC}",
                  f"{n_phase} phases x {p_rot} synthesized rotations: each QET phase controlled on the "
                  f"Hadamard-test ancilla. Not drawn in {PAPER}, whose trace step is QME. Literal end, "
                  f"eps_rot={eps_rot[1]:.4g} (SYN-BUDGET)"),
    )
    inter = {
        # the matrix
        "V": V, "L": L, "nnz_per_column": int(a.nnz_free_2028.lo), "nnz_formula_2D_plus_1": 2 * D + 1,
        "nnz_distinct_sites_this_L": len(row), "degenerate_M_is_m0_identity": degenerate,
        "W_diagonal": m0 * m0 + D * K * K / 2.0, "W_hop": -K * K / 4.0,
        "W_spectrum": spec, "lambda_min_W": lam_min, "lambda_max_W": lam_max, "kappa_W": lam_max / lam_min,
        "logdet_W": logdet,
        "be_subnormalization_s": s, "lambda_min_over_s": lam_norm, "lambda_max_over_s": lam_max / s,
        "log_s_over_lambda_min": log_norm, "V_log_s": v_log_s, "sum_log_lambda_over_s": logdet - v_log_s,
        "poly_shift": shift, "shift_term_subtracted": shift_term,
        "normalized_trace_shifted": (logdet - v_log_s) / (V * log_norm) + shift,
        # register
        "site_qubits": site, "column_qubits": col, "column_qubits_formula": neighbor_qubits_free(D),
        "column_qubits_gauged": neighbor_qubits(D),
        "be_flag_qubits": flag, "qet_signal_qubits": sig, "mcx_workspace_qubits": work,
        "hadamard_qubits": had, "overhead_V_independent": overhead,
        "projector_controls": col + flag, "projector_workspace_needed": col + flag - 2,
        "lq_hadamard_route": lq, "n_qpe": n_qpe, "qme_sign_qubits": sign, "qme_purification_qubits": site,
        "lq_qme_route": lq_qme,
        # polynomial
        "d": d, "d_is_for_this_instance": d_ok,
        "target_rel_error": rel, "poly_tolerance": tol,
        "lp_error_at_d": a.d_log_err.lo, "lp_error_at_d_minus_2": a.d_log_err_below.lo,
        "d_shifted_at_tolerance_1e-2": int(a.d_shifted_tol_1e2.lo),
        "d_unshifted_at_tolerance_1e-2": int(a.d_unshifted_tol_1e2.lo),
        "d_unshifted": d_u, "t_per_shot_unshifted": t_shot_u,
        "rotations_per_shot_unshifted": tuple(c_["n_rot"] for c_ in ends_u),
        "eps_rot_unshifted": tuple(c_["eps_rot"] for c_ in ends_u),
        "d_lp_at_lambda_min_0.16": int(a.d_lp_at_lam016.lo),
        "d_lp_at_tolerance_1e-3": int(a.d_lp_at_tol_1e3.lo),
        "d_paper_fig12": a.d_fig12_paper.value,
        "poly_abs_error_bound": poly_abs, "poly_rel_error_bound": poly_abs / abs(logdet),
        "sampling_eps": eps_s, "sampling_abs_error": samp_abs, "sampling_rel_error": samp_abs / abs(logdet),
        "sampling_is_rms": True, "shot_variance": a.shot_variance.lo,
        "shots_99pct_interval_factor": interval_99_factor(),
        "combined_worst_case_rel_error": (poly_abs + samp_abs) / abs(logdet),
        "synthesis_abs_error": synth_abs, "synthesis_rel_error": synth_abs / abs(logdet),
        "combined_worst_case_with_synthesis_rel_error": (poly_abs + samp_abs + synth_abs) / abs(logdet),
        "synth_to_poly": a.synth_to_poly_2028.lo, "eps_syn_report_default": EPS_SYN,
        # T
        "eps_syn": eps_syn, "eps_rot": eps_rot, "t_per_rotation": t_rot, "t_per_toffoli": tof_t,
        "oc_toffolis_per_query": oc_tof, "oc_t_per_query": tuple(toffoli_t(x, conv) for x in oc_tof),
        "oa_theta0": theta0, "oa_theta1": theta1,
        "oa_rotations_generic_angles": oa_generic,
        "oa_rotations_per_query": oa_rot, "t_per_query": t_query,
        "projector_phases": n_phase, "projector_toffolis_each": p_tof, "projector_rotations_each": p_rot,
        "t_per_projector_phase": t_phase,
        "be_queries": d,
        "t_cap": cap, "margin_vs_cap": margin,
        "rotations_per_shot": n_rot,
        "synthesis_error_randomized": synth,
        "synthesis_error_over_poly_tolerance": eps_syn / tol,
        "synthesis_error_if_coherent": tuple(n * e for n, e in zip(n_rot, eps_rot)),
        "eps_l_required": eps_l_req,
        "depth_budget_pf_over_eps_l": depth_budget(a.p_f.lo, a.eps_l_2028.lo),   # 1e7
        # sampling
        "shots": shots, "shots_qme_alternative": shots_qme(eps_s),
        "wall_per_shot_s": wall_shot, "wall_total_s": wall_total, "shot_overhead_s": t0,
        "wall_total_h": tuple(w / 3600 for w in wall_total),
        "t_total": t_total,
        # T-depth and factories (R9, R10): the 1 us/T baseline above needs F* >= 10; here F* = 2.6-4.8
        **{k: dx[k] for k in ("t_per_shot", "t_depth_per_shot", "f_star", "floor_wall_s", "factories_for_1yr")},
        "wall_serial_baseline_s": dx["wall_serial_s"], "baseline_ok": dx["baseline_ok"], "fits_1yr": dx["fits_1yr"],
        "and_tree_depth": and_depth, "t_depth_per_query": oa_rot[0] * t_rot[0],
        "t_depth_per_phase": and_depth + p_rot * t_rot[0],
        "wall_per_shot_corrected_s": wall_shot_corr, "wall_total_corrected_s": wall_corr,
        "wall_total_corrected_h": tuple(w / 3600 for w in wall_corr),
        "corrected_over_baseline": tuple(c_ / b for c_, b in zip(wall_corr, wall_total)),
        "wall_first_result_s": None,                                       # one tier: the benchmark itself
        "wall_campaign_s": wall_corr,
        # gauged comparison
        "gauged_units_D2Vlog2V": units, "gauged_units_D2Vlog2DV": scaling_units_logDV(D, V, a.log_base.lo),
        "gauged_t_per_shot": t_gauged, "gauged_factor_over_cap": factor,
        "gauged_fits_fault_budget": t_gauged <= depth_budget(a.p_f.lo, a.eps_l_2028.lo),
        "gauged_after_lcu": after_lcu, "gauged_after_levers": after_levers,
        "saving_levers": lev,
        "lcu_alone_closes_gap": after_lcu[0] <= cap,
        "gauged_saving_needed": factor,
        "levers_close_gap_at_upper_end": after_levers[0] <= cap,
        "levers_close_gap_at_lower_end": after_levers[1] <= cap,
        "gauged_over_cap_at_lever_top": after_levers[0] / cap,                # 1.376
    }
    notes = (
        f"R7: priced as the FREE-FERMION block encoding of {PAPER} Sec. II.2 at V = {L}^{D}. The paper gives the "
        "circuit structure (Figs. 1-3, Table I: O(D log D log V)) and no T-counts; every T figure is derived here.",
        f"LQ: {site} site + {col} column + {flag} flag + {sig} signal + {work} workspace + {had} Hadamard = {lq}; "
        f"QME route {lq_qme} (copy of the site register, {n_qpe} QPE, {sign} sign qubit, no Hadamard ancilla).",
        f"T: {d} queries x ({oc_tof[0]}-{oc_tof[1]} Toffolis + {oa_rot[0]}-{oa_rot[1]} rotations) + {n_phase} phases x "
        f"({p_tof} Toffolis + {p_rot} rotations) = {t_shot[0]:.4g}-{t_shot[1]:.4g}. The breakdown is the literal end.",
        f"The circuit block-encodes W/s, s = {s:.4g}: the QET sees lambda_min/s = {lam_norm:.4f}, not {lam_min:.2f}; "
        f"log det W = sum log(lambda/s) + V log s, and V log s = {v_log_s:.2f} is added back classically "
        "(the paper does not state this).",
        f"Target: the paper's shifted one, log|x|/log(s/lambda_min) + {shift:g} (Sec. IV). Classically the shift "
        f"term V log(s/lambda_min) x {shift:g} = {shift_term:.2f} is subtracted and V log s is added.",
        f"Accuracy, two kinds of claim. Polynomial: the uniform tolerance {tol:.4g} GUARANTEES "
        f"V log(s/lambda_min) x tol = {poly_abs:.4f} = {100 * rel:g}% of {logdet:.2f}. Sampling: eps = {eps_s:.4g} on "
        f"the normalized trace, {100 * rel:g}% rms, {shots:.5g} shots at Var = {a.shot_variance.lo:g} (R16). Combined worst case "
        f"{100 * (poly_abs + samp_abs) / abs(logdet):g}%.",
        f"d = {d} is an LP output carried as an input (err({d}) = {a.d_log_err.lo:.4e} <= tol < err({d - 2}) = "
        f"{a.d_log_err_below.lo:.4e}). Without the shift the same tolerance needs d = {d_u}: "
        f"{t_shot_u[0]:.3g}-{t_shot_u[1]:.3g} T per shot.",
        f"Synthesis (SYN-BUDGET; report convention T4): eps_syn = {a.synth_to_poly_2028.lo:g} x tol = "
        f"{eps_syn:.4g} per shot (report-wide default {EPS_SYN:g}); each end is its own circuit, "
        f"{n_rot[0]} / {n_rot[1]} rotations at eps_rot = {eps_rot[0]:.4g} / {eps_rot[1]:.4g}, "
        f"{t_rot[0]:.4g} / {t_rot[1]:.4g} T each (full RUS fit). eps_total = N eps^2 = {eps_syn:.4g} on the "
        f"normalized trace = {synth_abs:.4f} on log det W ({100 * synth_abs / abs(logdet):.3g}%). Combined worst case "
        f"with synthesis {100 * (poly_abs + samp_abs + synth_abs) / abs(logdet):.4g}% (was 1e-2 = 5.6x tol under R-TOL).",
        f"Gauged V = {L}^{D}: {d} x D^2 V log2 V = {d} x {units:.0f} = {t_gauged:.3g} T at one T per unit, "
        f"{factor:.1f}x the cap. LCU with the sublattice halving gives {lev[0]:g}-{lev[1]:g}x; at the top "
        f"{after_levers[0]:.4g} T, {after_levers[0] / cap:.3g}x the cap, so it does not fit at any lever setting "
        "(R11b: the QSP-degree lever is dropped; d is at the LP floor).",
        f"T-depth (R9, factory.json): {d} x {oa_rot[0] * t_rot[0]:.4g} + {n_phase} x {and_depth + p_rot * t_rot[0]:.4g} = "
        f"{t_depth:.4g}; F* = {dx['f_star'][0]:.2f}-{dx['f_star'][1]:.2f} < 10, so the 1 us/T baseline cannot be fed "
        f"(R10). A shot takes {1e3 * wall_shot_corr[0]:.1f} ms (T-depth x 10 us + t0), the run {wall_corr[0] / 3600:.2f} h "
        f"(baseline {wall_total[0] / 3600:.2f}-{wall_total[1] / 3600:.2f} h). One factory suffices for a year.",
    )
    return Result(era="2028", lq=(lq, lq), hard_ops=t_shot, breakdown=breakdown,
                  intermediates=inter, shots=(shots, shots), wall_time_s=wall_corr,
                  epsilon_l=eps_l_req, notes=notes)


def _model_2033(a: Assumptions) -> Result:
    """1000-LQ rung after rulings R1, R3, R7, R10 (H. Lamm, 2026-10-02): Tr M^-1 on 8^3x16 at kappa ~ 1e2, the
    block encoding of M at D V log2 V per query, read out as the block-encoding success probability
    P = f^2 am (1/V_F) Re Tr M^-1. Two tiers: one configuration at 30% (first result), three at 20% (campaign).
    Walls are depth-bound (F* = 1-4 < 10). The Banks-Casher mode number is the demonstration row (R7)."""
    D, V, Nc, Nf = int(a.D.lo), int(a.V_2033.lo), int(a.Nc_2033.lo), int(a.Nf_2033.lo)
    base, pref = a.log_base.lo, a.scaling_prefactor.lo
    n_qpe, n_anc = int(a.n_qpe_2033.lo), int(a.n_anc_2033.lo)
    p_D = int(a.be_D_power_M.lo)
    tg, t0 = a.t_gate_s.lo, a.shot_overhead_s.lo
    yr = SECONDS_PER_YEAR

    # --- register (Eq. Nq with the column index of M, R3) ------------------------------------------------
    site, col_M, col = site_qubits(V), column_qubits_M(D), color_qubits(Nc)
    lq = site + col_M + col + n_qpe + n_anc                     # 13+4+4+12+50 = 83
    lq_box_sum = site + int(a.neighbor_box_2033.lo) + col + n_qpe + n_anc   # 83 (box addends)
    lq_incoherent = lq - n_qpe                                  # 71
    lq_16x4 = register_M(int(a.V_16x4.lo), D, Nc, n_qpe, n_anc)   # 86
    lq_128x4 = register_M(int(a.V_128x4.lo), D, Nc, n_qpe, n_anc) # 98
    lq_W_column = register(V, D, Nc, n_qpe, n_anc)              # 85: the M^dag M column index, before R3
    VF = int(a.n_s.lo) * Nc * Nf * V                            # 7.4e4 flavor-stacked (not printed)
    VF_single = int(a.n_s.lo) * Nc * V                          # 2.5e4 per flavor

    # --- per-shot T (Eq. Ngate, M at D V log2 V, R3) --------------------------------------------------------
    d_inv = a.d_inv_2033.lo                                     # 500 (LP, R11)
    d_inv_law = d_inv_lp(a.kappa_2033.lo, a.d_inv_eps_rel.lo)   # 500.0
    d_inv_ln = d_inv_formula(a.kappa_2033.lo, a.eps_2033.lo, math.e)   # 921
    d_inv_log10 = d_inv_formula(a.kappa_2033.lo, a.eps_2033.lo, 10)    # 400
    units = scaling_units_M(D, V, base, p_D)                    # 4.26e5 = 4 x 8192 x 13
    units_W = scaling_units(D, V, base)                         # 1.70e6, the M^dag M price before R3
    t_query = pref * units
    t_shot = d_inv * t_query                                    # 2.13e8
    t_shot_W = d_inv * pref * units_W                           # 8.52e8 (before R3)
    budget = depth_budget(a.p_f.lo, a.eps_l_2033.lo)           # 1e9
    eps_l_req = a.p_f.lo / t_shot                               # 4.69e-10 at 0.1 expected faults per shot

    # --- readout: block-encoding success probability (R3; shot audit readout lever 1) ------------------------
    am, f = a.am_2033.lo, a.peak_norm_2033.lo
    kappa_from_am = (4.0 + am) / am                             # 101
    c_free = free_condensate(int(a.L_s_2033.lo), int(a.L_t_2033.lo), am)   # 0.0275
    c_int = a.condensate_interacting_2033.lo                    # 0.08
    P = (f * f * am * c_free, f * f * am * c_int)               # (2.75e-4, 8.0e-4)
    mu_had = (0.5 * am * c_free, 0.5 * am * c_int)              # Hadamard-test mean (m/2) c: 5.5e-4 to 1.6e-3
    r1, r2 = a.rel_first_2033.lo, a.rel_campaign_2033.lo
    n1, n2 = a.n_cfg_first_2033.lo, a.n_cfg_campaign_2033.lo
    shots_cfg_first = (shots_success_probability(P[1], r1), shots_success_probability(P[0], r1))   # (1.39e4, 4.04e4)
    shots_cfg_camp = (shots_success_probability(P[1], r2), shots_success_probability(P[0], r2))    # (3.12e4, 9.09e4)
    shots_first = tuple(n1 * x for x in shots_cfg_first)
    shots_camp = tuple(n2 * x for x in shots_cfg_camp)          # (9.37e4, 2.73e5)
    shots_had_first = (1.0 / (r1 * r1 * mu_had[1] ** 2), 1.0 / (r1 * r1 * mu_had[0] ** 2))   # 30% by Hadamard test
    p_over_had = (shots_had_first[0] / shots_cfg_first[0], shots_had_first[1] / shots_cfg_first[1])
    old_eps_over_mu = (a.eps_2033.lo / mu_had[1], a.eps_2033.lo / mu_had[0])   # (6.3, 18): the retired plan

    # --- T-depth and factories (R9, R10; factory.json, ratio carried to the M count) --------------------------
    par = a.parallelism_2033
    t_depth = (t_shot / par.hi, t_shot / par.lo)                # (5.3e7, 2.13e8)
    dx1 = depth_exports(t_shot, t_depth, shots_first, tg, t0)
    dx2 = depth_exports(t_shot, t_depth, shots_camp, tg, t0)
    w_shot = (corrected_shot_s(t_shot, t_depth[0], tg, t0), corrected_shot_s(t_shot, t_depth[1], tg, t0))   # (533, 2130) s
    w_shot_base = t_shot * tg + t0                              # 213 s at the 10-factory baseline
    wall_first = (shots_first[0] * w_shot[0], shots_first[1] * w_shot[1])      # (7.39e6, 8.60e7) s
    wall_camp = (shots_camp[0] * w_shot[0], shots_camp[1] * w_shot[1])         # (4.99e7, 5.81e8) s
    wall_camp_qrom = (shots_camp[0] * w_shot[0], shots_camp[1] * w_shot[0])    # F* = 4 only: (1.58, 4.60) yr
    horizon_s = a.campaign_horizon_yr.lo * yr
    horizon_ratio = (wall_camp[0] / horizon_s, wall_camp[1] / horizon_s)       # (0.32, 3.68)
    k_needed = math.ceil(FACTORY_BASELINE / par.hi)             # 3 parallel subtrees -> F* ~ 12
    extra_lq_fix = k_needed * a.subtree_lq.lo                   # ~100 idle LQ (unpriced)

    # --- Banks-Casher mode number (R7; simplobs.json C2 at a = 0.15 fm) --------------------------------------
    a_inv = a.bc_hbarc_gev_fm.lo / a.bc_a_fm.lo                 # 1.315 GeV
    s_lat = a.bc_tastes.lo * (a.bc_sigma_third_gev.lo / a_inv) ** 3   # a^3 Sigma x tastes
    bc = []
    for lam, deg in ((a.bc_lambda_hi.lo, a.bc_d_lambda_hi.lo), (a.bc_lambda_lo.lo, a.bc_d_lambda_lo.lo)):
        P_bc = 2.0 * lam * s_lat / (math.pi * Nc)               # nu / V_F
        nu = P_bc * VF_single                                   # modes below Lambda per configuration
        var_rel = r1 * r1 - 1.0 / (a.bc_n_cfg.lo * nu)          # 30% total, less the gauge variance of nu
        n_bc = (1.0 - P_bc) / (P_bc * var_rel)
        t_bc = deg * t_query
        bc.append({"lambda": lam, "lambda_mev": 1e3 * lam * a_inv, "d": deg, "P": P_bc, "nu": nu,
                   "shots": n_bc, "t_per_shot": t_bc})
    bc_t = (bc[0]["t_per_shot"], bc[1]["t_per_shot"])           # (1.83e8, 2.73e8)
    bc_shots = (bc[0]["shots"], bc[1]["shots"])                 # (1.10e4, 1.73e4)
    bc_wall = (bc_shots[0] * corrected_shot_s(bc_t[0], bc_t[0] / par.hi, tg, t0),
               bc_shots[1] * corrected_shot_s(bc_t[1], bc_t[1] / par.lo, tg, t0))   # (5.0e6, 4.7e7) s
    bc_wall_base = tuple(n * (t * tg + t0) for n, t in zip(bc_shots, bc_t))

    n_phase = int(d_inv) + 1
    t_phase = n_phase * t_per_rotation(a.qsp_phase_eps.lo)

    breakdown = (
        Primitive("BE_query_M_SU3_staggered_8^3x16", d_inv, t_query, CircuitStatus.SCALING,
                  "app12 1000-LQ box; arxiv_2407_13080; ruling R3",
                  f"'d_inv . D V log2 V ~ 500 . 4.3e5'; D V log2 V = {units:.4g} for M (9 entries per column) "
                  f"with unit prefactor; no compiled circuit"),
    )
    inter = {
        "site_qubits": site, "neighbor_qubits": col_M, "color_qubits": col,
        "neighbor_qubits_box": int(a.neighbor_box_2033.lo),
        "neighbor_qubits_M_dag_M": neighbor_qubits(D),                                  # 6, before R3
        "n_qpe": n_qpe, "n_anc": n_anc, "overhead_V_independent": n_qpe + n_anc,   # 62
        "lq_eq": lq, "lq_box_sum": lq_box_sum, "lq_incoherent_route": lq_incoherent,
        "lq_with_W_column_index": lq_W_column,
        "lq_16x4": lq_16x4, "lq_128x4": lq_128x4, "site_qubits_128x4": site_qubits(int(a.V_128x4.lo)),
        "qubits_per_4096x_volume": ceil_log2(4096),                                  # 12
        "V_F": VF, "V_F_single_flavor": VF_single, "V_F_box": a.VF_box_2033.lo,
        "kappa": a.kappa_2033.lo, "am": am, "kappa_from_am": kappa_from_am,
        "d_inv": d_inv, "d_inv_lp_law": d_inv_law, "d_inv_eps_rel": a.d_inv_eps_rel.lo,
        "d_inv_formula_ln": d_inv_ln, "d_inv_formula_log10": d_inv_log10,
        "d_inv_implied_prefactor_ln": d_inv / d_inv_ln,                              # 0.54
        "scaling_units_DVlogV": units, "scaling_units_D2VlogV": units_W,
        "scaling_units_D2VlogDV": scaling_units_logDV(D, V, base),
        "t_per_query": t_query,
        "t_per_shot_before_R3": t_shot_W, "reprice_factor": t_shot_W / t_shot,       # 4 = D
        "t_per_shot_d_inv_formula_ln": d_inv_ln * t_query,
        "t_cap": a.t_cap_2033.lo,
        "depth_budget_pf_over_eps_l": budget, "fits_budget": t_shot <= budget,
        "eps_l_required": eps_l_req,
        # readout
        "peak_norm": f, "condensate_free": c_free, "condensate_interacting": c_int,
        "P_success": P, "mu_hadamard": mu_had,
        "rel_first": r1, "n_cfg_first": n1, "rel_campaign": r2, "n_cfg_campaign": n2,
        "shots_per_cfg_first": shots_cfg_first, "shots_per_cfg_campaign": shots_cfg_camp,
        "shots_first": shots_first, "shots_campaign": shots_camp,
        "shots_hadamard_first": shots_had_first, "p_readout_gain": p_over_had,
        "retired_eps_over_mu": old_eps_over_mu,
        "t_total_first": tuple(n * t_shot for n in shots_first),
        "t_total_campaign": tuple(n * t_shot for n in shots_camp),
        # T-depth, factories and walls (R9, R10)
        **{k: dx2[k] for k in ("t_per_shot", "t_depth_per_shot", "f_star", "floor_wall_s", "factories_for_1yr")},
        "factories_for_1yr_first": dx1["factories_for_1yr"], "floor_wall_first_s": dx1["floor_wall_s"],
        "baseline_ok": dx2["baseline_ok"], "fits_1yr": dx2["fits_1yr"], "fits_1yr_first": dx1["fits_1yr"],
        "wall_per_shot_s": w_shot, "wall_per_shot_baseline_s": w_shot_base,
        "wall_serial_baseline_first_s": dx1["wall_serial_s"], "wall_serial_baseline_campaign_s": dx2["wall_serial_s"],
        "wall_first_result_s": wall_first, "wall_campaign_s": wall_camp,
        "wall_first_result_days": tuple(w / 86400.0 for w in wall_first),
        "wall_first_result_yr": tuple(w / yr for w in wall_first),
        "wall_campaign_yr": tuple(w / yr for w in wall_camp),
        "wall_campaign_unary_iteration_yr": tuple(w / yr for w in wall_camp_qrom),
        "wall_campaign_baseline_yr": tuple(w / yr for w in dx2["wall_serial_s"]),
        "shot_overhead_s": t0, "campaign_horizon_yr": a.campaign_horizon_yr.lo,
        "horizon_ratio_campaign": horizon_ratio,
        "fits_horizon_one_machine": (wall_camp[0] <= horizon_s, wall_camp[1] <= horizon_s),
        "subtrees_for_baseline": k_needed, "subtree_f_star": k_needed * par.hi,
        "extra_lq_for_baseline": extra_lq_fix,
        # Banks-Casher demonstration (R7)
        "bc_a_inv_gev": a_inv, "bc_a3_sigma_tastes": s_lat, "bc": tuple(bc),
        "bc_t_per_shot": bc_t, "bc_shots": bc_shots, "bc_wall_s": bc_wall,
        "bc_wall_days": tuple(w / 86400.0 for w in bc_wall), "bc_wall_yr": tuple(w / yr for w in bc_wall),
        "bc_wall_baseline_days": tuple(w / 86400.0 for w in bc_wall_base),
        # QME and synthesis (unchanged rules)
        "t_per_shot_qme_coherent": t_shot / a.eps_2033.lo,                          # 2.1e10: ~1e2 QETs in one circuit
        "qsp_phase_rotations": n_phase, "qsp_phase_rotations_t": t_phase,
        "synth_to_poly": a.synth_to_poly_2033.lo,
        "poly_tolerance_normalized": 0.5 * a.d_inv_eps_rel.lo,                         # 3.369e-3
        "eps_syn_2033": a.synth_to_poly_2033.lo * 0.5 * a.d_inv_eps_rel.lo,            # 3.369e-4
        "eps_syn_2033_frac_of_peak": a.synth_to_poly_2033.lo * a.d_inv_eps_rel.lo,     # 6.74e-4 (0.067%)
        "synth_over_poly_at_report_default": EPS_SYN / (0.5 * a.d_inv_eps_rel.lo),   # 2.97
        "synth_extra_t_per_rotation": 1.15 * math.log2(math.sqrt(
            EPS_SYN / (a.synth_to_poly_2033.lo * 0.5 * a.d_inv_eps_rel.lo))),          # +2.81 T (absorbed)
    }
    notes = (
        f"R3: block-encodes M (2D+1 = 9 entries, {col_M}-qubit column index): {site} site + {col_M} column + {col} "
        f"color + {n_qpe} QPE + {n_anc} anc = {lq} (was {lq_W_column} with the M^dag M index).",
        f"T: d_inv = {d_inv:g} x D V log2 V = {units:.4g} -> {t_shot:.4g} (was {t_shot_W:.4g} at D^2 V log2 V).",
        f"Readout (R3): P = f^2 am c = {P[0]:.3g}-{P[1]:.3g} at f = {f:g}, am = {am:g}, c = {c_free:.4f} (free) to "
        f"{c_int:g}. The Hadamard test would need {shots_had_first[0]:.3g}-{shots_had_first[1]:.3g} shots for 30%; "
        f"the retired eps = 1e-2 was {old_eps_over_mu[0]:.1f}-{old_eps_over_mu[1]:.0f}x the signal.",
        f"First result: {n1:g} configuration at {100 * r1:g}%: {shots_first[0]:.4g}-{shots_first[1]:.4g} shots, "
        f"{wall_first[0] / 86400:.0f} d - {wall_first[1] / yr:.2f} yr. Campaign: {n2:g} at {100 * r2:g}%: "
        f"{shots_camp[0]:.4g}-{shots_camp[1]:.4g} shots, {wall_camp[0] / yr:.2f}-{wall_camp[1] / yr:.1f} yr.",
        f"F* = {dx2['f_star'][0]:.0f}-{dx2['f_star'][1]:.0f} (factory.json): a shot is {w_shot[0]:.0f}-{w_shot[1]:.0f} s, "
        f"not {w_shot_base:.0f} s. {k_needed} parallel subtrees (~{extra_lq_fix:g} idle LQ) restore the baseline: "
        f"{dx2['wall_serial_s'][0] / yr:.2f}-{dx2['wall_serial_s'][1] / yr:.2f} yr.",
        f"Banks-Casher (R7): {bc_t[0]:.3g}-{bc_t[1]:.3g} T, {bc_shots[0]:.3g}-{bc_shots[1]:.3g} shots, "
        f"{bc_wall[0] / 86400:.0f} d - {bc_wall[1] / yr:.2f} yr.",
    )
    return Result(era="2033", lq=(lq, lq), hard_ops=(t_shot, t_shot), breakdown=breakdown,
                  intermediates=inter, shots=shots_camp, wall_time_s=wall_camp,
                  epsilon_l=(eps_l_req, eps_l_req), notes=notes)


def _model_codesign(a: Assumptions) -> Result:
    """The 1e4-LQ rung, post-2033 (ruling R7). Tr M^-1 at 24^3x48, kappa ~ 1e3, block-encoding M at D V log2 V
    (R3): ~2.6e11 T on a 90-LQ register before the per-shot reduction."""
    D, Nc = int(a.D.lo), int(a.Nc_2033.lo)
    base, pref = a.log_base.lo, a.scaling_prefactor.lo
    n_qpe, n_anc = int(a.n_qpe_codesign.lo), int(a.n_anc_codesign.lo)
    V, Vp = int(a.V_codesign.lo), int(a.V_production.lo)
    p_D = int(a.be_D_power_M.lo)
    tg, t0 = a.t_gate_s.lo, a.shot_overhead_s.lo

    lq = register_M(V, D, Nc, n_qpe, n_anc)                     # 20+4+4+12+50 = 90
    lq_prod = register_M(Vp, D, Nc, n_qpe, n_anc)               # 94
    d_inv = a.d_inv_codesign.lo                                 # 5e3 (LP law, R11)
    d_inv_law = d_inv_lp(a.kappa_codesign.lo, a.d_inv_eps_rel.lo)   # 5000.0
    d_inv_ln = d_inv_formula(a.kappa_codesign.lo, a.eps_2033.lo, math.e)   # 1.15e4
    d_inv_2033_ln = d_inv_formula(a.kappa_2033.lo, a.eps_2033.lo, math.e)  # 921
    units = scaling_units_M(D, V, base, p_D)                    # 5.13e7
    t_query = pref * units
    t_shot = d_inv * t_query                                    # 2.567e11 (box ~2.6e11)
    t_shot_W = d_inv * pref * scaling_units(D, V, base)         # 1.027e12 before R3
    budget = depth_budget(a.p_f.lo, a.eps_l_codesign.lo)       # 1e11
    red = a.per_shot_reduction.lo                               # 2.6
    t_shot_reduced = t_shot / red                               # 9.87e10: fits
    t_shot_prod = d_inv * pref * scaling_units_M(D, Vp, base, p_D)   # 4.96e12
    t_shot_prod_reduced = t_shot_prod / red
    sub = a.sublattice_saving.lo
    levers = (a.lcu_saving.lo * sub, a.lcu_saving.hi * sub)     # (6, 20): LCU x even-odd; degree at the LP floor

    shots_cfg = shots_incoherent(a.eps_2033.lo, a.shot_variance.lo)   # 1e4 (R16 rms, inherited)
    n_cfg_prod = a.n_cfg_production.lo                          # 1e3
    calls = shots_cfg * n_cfg_prod                              # 1e7
    par = a.parallelism_codesign
    depth_reduced = (t_shot_reduced / par.hi, t_shot_reduced / par.lo)
    w_shot = (corrected_shot_s(t_shot_reduced, depth_reduced[0], tg, t0),
              corrected_shot_s(t_shot_reduced, depth_reduced[1], tg, t0))         # (2.47e5, 9.87e5) s
    w_shot_base = t_shot_reduced * tg + t0                      # 9.87e4 s at the 10-factory baseline
    wall = (calls * w_shot[0], calls * w_shot[1])               # (2.47e12, 9.87e12) s
    horizon_s = a.campaign_horizon_yr.lo * SECONDS_PER_YEAR
    horizon_ratio = (wall[0] / horizon_s, wall[1] / horizon_s)  # (1.6e4, 6.3e4)
    shots_in_horizon = (horizon_s / w_shot[1], horizon_s / w_shot[0])   # (160, 639)
    shots_cfg_qme = shots_qme(a.eps_2033.lo)                    # 100 QETs per config, coherent
    # referee H2: the disconnected-HVP loops L_mu(t) = Tr(Gamma_{mu,t} M^-1) are not positive, so the P readout does
    # not apply; on the Hadamard route a shot's mean is p0 am <x|Gamma M^-1|x>. Matching the shot noise on every loop
    # to its gauge fluctuation (variance hvp_site_variance per site) needs V_F / (var (p0 am)^2) shots per configuration:
    # the estimator V_F/(p0 am) mean(b) has shot variance V_F^2/((p0 am)^2 N); setting it to V_F var gives N. The count
    # scales inversely with var, so a site variance below one raises it. One Gamma_mu insertion per shot: the count
    # is per flavor and per current direction (x3 for the spatial HVP currents).
    am_prod = 4.0 / (a.kappa_codesign.lo - 1.0)                 # kappa = (4 + am)/am -> am = 4/(kappa - 1) = 0.004
    vf_prod = Nc * Vp                                           # 3.19e7 per flavor (n_s = 1)
    vf_cd = Nc * V                                              # 1.99e6 per flavor at 24^3x48
    p0 = a.peak_norm_2033.lo
    hvp_shots_cfg = vf_prod / (a.hvp_site_variance.lo * (p0 * am_prod) ** 2)    # 7.95e12 at var = 1
    hvp_shots_cfg_cd = vf_cd / (a.hvp_site_variance.lo * (p0 * am_prod) ** 2)   # 4.97e11 at var = 1 (24^3x48)
    hvp_current_dirs = 3                                        # spatial Gamma_mu, one insertion per shot
    # referee H2 mapping (2026-10-05; scratch chopen/ch04_h2_hvp/hvp_map.py). C_disc(t) ~ <(L_l-L_s)(t)(L_l-L_s)(0)>/9
    # from the loops L_k(t) = Tr(Gamma_{k,t} M^-1); staggered eps-symmetry makes them purely imaginary (Im Hadamard
    # test). Gamma_k = two signed shifts on loaded links, BE at alpha = 1, one query against d_inv queries of M:
    # ~V log2 V against d_inv D V log2 V, a share 1/(D d_inv). A shot reads one diagonal element, so its variance on
    # L_k(t) is (3 V_s)^2/(p0 am)^2, against <= 3 V_s ||Gamma M^-1||^2 <= 3 V_s/am^2 for a slice-diluted noise vector
    # (one sparse solve): ratio >= 3 V_s/p0^2. QME per slice at noise matching (delta^2 = 3 V_s var): sqrt(3 V_s var)
    # /(p0 am var) calls; classical vectors per slice <= 1/(am^2 var).
    l_t_cd = 48                                                 # 24^3x48 time extent (V = 24^3 x 48)
    v_s_cd = V // l_t_cd                                        # 13824
    hvp_insertion_share = 1.0 / (D * d_inv)                     # 5e-5 of the per-shot T
    hvp_var_ratio_lb_cd = Nc * v_s_cd / p0 ** 2                 # 1.66e5: Hadamard shot vs classical noise vector
    hvp_qme_calls_cfg_cd = hvp_current_dirs * l_t_cd * math.sqrt(Nc * v_s_cd * a.hvp_site_variance.lo) / (
        p0 * am_prod * a.hvp_site_variance.lo)                  # 1.46e7
    hvp_classical_vectors_cfg_ub_cd = hvp_current_dirs * l_t_cd / (am_prod ** 2 * a.hvp_site_variance.lo)  # 8.98e6

    breakdown = (
        Primitive("BE_query_M_SU3_staggered_24^3x48", d_inv, t_query, CircuitStatus.SCALING,
                  "app12 1e4-LQ paragraph; arxiv_2407_13080; ruling R3",
                  f"'d_inv ~= 5e3 costs ~2.6e11 T'; D V log2 V = {units:.4g} (M, 9 entries per column)"),
    )
    inter = {
        "site_qubits": site_qubits(V), "column_qubits": column_qubits_M(D), "color_qubits": color_qubits(Nc),
        "n_qpe": n_qpe, "n_anc": n_anc,
        "lq_eq": lq, "lq_with_W_column_index": register(V, D, Nc, n_qpe, n_anc),     # 92 before R3
        "site_qubits_production": site_qubits(Vp), "lq_production": lq_prod,
        "kappa": a.kappa_codesign.lo, "d_inv": d_inv, "d_inv_lp_law": d_inv_law,
        "d_inv_formula_ln": d_inv_ln,
        "d_inv_implied_prefactor_ln": d_inv / d_inv_ln,                                  # 0.43
        "d_inv_ratio_vs_2033": d_inv / a.d_inv_2033.lo,                                  # 10 (ratio of the kappas)
        "d_inv_ratio_vs_2033_formula_ln": d_inv_ln / d_inv_2033_ln,                      # 12.5
        "scaling_units_DVlogV": units, "t_per_shot": t_shot, "t_per_shot_before_R3": t_shot_W,
        "depth_budget_pf_over_eps_l": budget, "fits_budget": t_shot <= budget,
        "factor_over_budget": t_shot / budget,                                           # 2.57
        "factor_over_rfi_2033_cap": t_shot / a.t_cap_2033.lo,                            # 257
        "per_shot_reduction": red,
        "t_per_shot_after_reduction": t_shot_reduced,
        "fits_budget_after_reduction": t_shot_reduced <= budget,
        "lever_range": levers,
        "fits_budget_at_lever_top": t_shot / levers[1] <= budget,
        "fits_budget_at_lever_bottom": t_shot / levers[0] <= budget,
        "lcu_alone_closes": t_shot / a.lcu_saving.lo <= budget,                          # 3x: yes
        "even_odd_alone_closes": t_shot / sub <= budget,                                 # 2x: no
        "t_per_shot_production": t_shot_prod,
        "t_per_shot_production_after_reduction": t_shot_prod_reduced,
        "production_fits_after_reduction": t_shot_prod_reduced <= budget,
        "production_reduction_needed": t_shot_prod / budget,                             # 49.6 (~50)
        "production_beyond_levers": t_shot_prod / budget > levers[1],
        "shots_per_cfg": shots_cfg, "n_cfg_production": n_cfg_prod,
        "subroutine_calls_1e4_rung": calls,
        "subroutine_calls_1e4_rung_box": a.subroutine_calls_box_1e4.lo,
        "shot_overhead_s": t0, "campaign_horizon_yr": a.campaign_horizon_yr.lo,
        "f_star": (par.lo, par.hi), "t_depth_after_reduction": depth_reduced,
        "wall_per_shot_after_reduction_s": w_shot,                                       # (2.47e5, 9.87e5) s
        "wall_per_shot_after_reduction_days": tuple(w / 86400.0 for w in w_shot),        # (2.9, 11.4) d
        "wall_per_shot_baseline_s": w_shot_base,
        "wall_serial_after_reduction_s": wall,
        "wall_serial_after_reduction_yr": tuple(w / SECONDS_PER_YEAR for w in wall),     # (7.8e4, 3.1e5) yr
        "fits_horizon_one_machine": wall[0] <= horizon_s,                                # False
        "horizon_ratio_serial": horizon_ratio,
        "cost_reduction_to_fit_horizon": horizon_ratio,
        "shots_in_horizon": shots_in_horizon,                                            # (160, 639)
        "cfg_fraction_in_horizon": tuple(x / shots_cfg for x in shots_in_horizon),
        "qme_call_saving": shots_cfg / shots_cfg_qme,                                    # 1e2 (1/eps vs 1/eps^2)
        "am_production": am_prod, "v_f_production": vf_prod,
        "hvp_noise_matched_shots_per_cfg": hvp_shots_cfg,                                # ~8e12 (referee H2)
        "hvp_noise_matched_over_box": hvp_shots_cfg / shots_cfg,                         # ~8e8 x the 1e4 lower bound
        "v_f_codesign": vf_cd,
        "hvp_noise_matched_shots_per_cfg_24c48": hvp_shots_cfg_cd,                      # ~5e11 (same lattice as box)
        "hvp_noise_matched_over_box_24c48": hvp_shots_cfg_cd / shots_cfg,               # ~5e7
        "hvp_current_directions": hvp_current_dirs,
        "hvp_noise_matched_shots_per_cfg_all_dirs": hvp_current_dirs * hvp_shots_cfg,   # 2.4e13 at 48^3x96
        "hvp_noise_matched_shots_per_cfg_all_dirs_24c48": hvp_current_dirs * hvp_shots_cfg_cd,  # 1.5e12
        "hvp_noise_matched_over_box_all_dirs_24c48": hvp_current_dirs * hvp_shots_cfg_cd / shots_cfg,  # 1.5e8
        "hvp_insertion_share_of_shot": hvp_insertion_share,                              # 5e-5
        "hvp_shot_vs_noise_vector_variance_ratio_lb_24c48": hvp_var_ratio_lb_cd,         # 1.66e5
        "hvp_qme_calls_per_cfg_24c48": hvp_qme_calls_cfg_cd,                             # 1.46e7
        "hvp_classical_vectors_per_cfg_ub_24c48": hvp_classical_vectors_cfg_ub_cd,       # 8.98e6
        "hvp_noise_matched_wall_per_cfg_yr_24c48": tuple(hvp_current_dirs * hvp_shots_cfg_cd * w / SECONDS_PER_YEAR
                                                         for w in w_shot),               # (1.2e10, 4.7e10) yr
        "hvp_qme_coherent_t_per_circuit_24c48": hvp_qme_calls_cfg_cd / (hvp_current_dirs * l_t_cd) * t_shot_reduced,
        "hvp_qme_coherent_over_budget_24c48": hvp_qme_calls_cfg_cd / (hvp_current_dirs * l_t_cd) * t_shot_reduced
                                              / budget,                          # 1.0e5 (H2 verifier #1)
        "synth_to_poly": a.synth_to_poly_2033.lo,
        "poly_tolerance_normalized": 0.5 * a.d_inv_eps_rel.lo,                         # 3.369e-3
        "eps_syn_2033": a.synth_to_poly_2033.lo * 0.5 * a.d_inv_eps_rel.lo,            # 3.369e-4
        "eps_syn_2033_frac_of_peak": a.synth_to_poly_2033.lo * a.d_inv_eps_rel.lo,     # 6.74e-4 (0.067%)
        "synth_over_poly_at_report_default": EPS_SYN / (0.5 * a.d_inv_eps_rel.lo),   # 2.97
        "synth_extra_t_per_rotation": 1.15 * math.log2(math.sqrt(
            EPS_SYN / (a.synth_to_poly_2033.lo * 0.5 * a.d_inv_eps_rel.lo))),          # +2.81 T (absorbed)
    }
    notes = (
        f"R3 + R7: post-2033 rung, block-encodes M: {site_qubits(V)} site + {column_qubits_M(D)} column + "
        f"{color_qubits(Nc)} color + {n_qpe} QPE + {n_anc} anc = {lq} LQ (was 92); 48^3x96 gives {lq_prod}.",
        f"{t_shot:.4g} T (was {t_shot_W:.4g}) is {t_shot / budget:.2f}x the 1e11 budget at eps_l = 1e-12; the "
        f"{red:g}x reduction leaves {t_shot_reduced:.4g}. LCU alone (3-10x) closes it; even-odd alone (2x) does not.",
        f"48^3x96: {t_shot_prod:.3g} T needs ~{t_shot_prod / budget:.0f}x, beyond the {levers[1]:g}x of the levers.",
        f"One machine, F* = {par.lo:g}-{par.hi:g}: {w_shot[0]:.3g}-{w_shot[1]:.3g} s per shot; {calls:.2g} calls are "
        f"{wall[0] / SECONDS_PER_YEAR:.3g}-{wall[1] / SECONDS_PER_YEAR:.3g} yr. QME supplies 1e2x at 1e2x the depth. "
        "The 1e4 shots per configuration are a lower bound for disconnected HVP: noise-matching the current loops "
        f"needs ~{hvp_shots_cfg_cd:.2g} per configuration per flavor and current direction on 24^3x48 and "
        f"~{hvp_shots_cfg:.2g} on 48^3x96 at order-one site variance; a smaller variance raises it (referee H2). "
        f"Three currents: {hvp_current_dirs * hvp_shots_cfg_cd:.2g} per cfg at 24^3x48; a shot has >= "
        f"{hvp_var_ratio_lb_cd:.2g}x the variance of a classical noise vector; QME ~{hvp_qme_calls_cfg_cd:.2g} calls "
        f"against <= {hvp_classical_vectors_cfg_ub_cd:.2g} classical vectors per cfg.",
    )
    return Result(era="codesign", lq=(lq, lq), hard_ops=(t_shot, t_shot), breakdown=breakdown,
                  intermediates=inter, shots=None, wall_time_s=None,
                  epsilon_l=(a.eps_l_codesign.lo, a.eps_l_codesign.lo), notes=notes)


def model(a: Assumptions, era: str) -> Result:
    """Return the Result for '2028' | '2033' | 'codesign'."""
    if era == "2028":
        return _model_2028(a)
    if era == "2033":
        return _model_2033(a)
    if era == "codesign":
        return _model_codesign(a)
    raise ValueError(f"unknown era {era!r}; expected one of {ERAS}")


PUBLISHED = {
    "2028": Published(lq=(18, 18), hard_ops=(1.9e4, 3.5e4),
                      src="app12:2028 box lines 119-122: 'Logical qubits 18: 8 site + 4 column index + 1 "
                          "block-encoding flag + 1 QET signal + 3 multi-control workspace + 1 Hadamard-test "
                          "ancilla' (QME upgrade path 34); 'Per-shot 1.9--3.5e4 T-gates, derived here' (ruling "
                          "SYN-BUDGET 2026-09-29, eps_syn = 0.1 x 1.784e-3, eps_rot = sqrt(eps_syn/N_rot); was "
                          "1.8--3.4e4 under R-TOL at eps_syn = 1e-2, Clifford+T-exact R_y(theta_0 = 0) "
                          "not synthesized; 2.0--4.1e4 at the fixed 1e-4 under 'SHIFT + 1%', d = 84; 2.7--5.4e4 "
                          "at d = 112 under R7; '~30 LQ, ~3e4 T' at V=2^4 before it). Model 1.898e4-3.544e4: 0.1%, 1.3%",
                      rel_tol=0.02),
    "2033": Published(lq=(83, 83), hard_ops=(2.1e8, 2.1e8),
                      src="app12 1000-LQ box: 'Logical qubits 83 (13 site + 4 column index + 4 color + 12 QPE + 50 "
                          "ancilla)'; 'Per-shot T 2.1e8' = d_inv D V log2 V = 500 x 4.26e5 (ruling R3, 2026-10-02: M has "
                          "9 entries per column; was 85 LQ, '~9e8' at D^2 V log2 V). Model 2.130e8 vs 2.1e8: 1.4%",
                      rel_tol=0.02),
    "codesign": Published(lq=(90, 90), hard_ops=(2.6e11, 2.6e11),
                          src="app12 1e4-LQ paragraph: 'V ~ 24^3x48 ... d_inv ~= 5e3 costs 2.6e11 T by Eq. Ngate_hybrid "
                              "on a 90-LQ algorithmic register (20 site + 4 column index + 4 color + 12 QPE + ~50 "
                              "ancilla)' (R3 reprice at D V log2 V, 2026-10-02; was 92 LQ, ~1e12 T). Model 2.567e11: 1.3%",
                          rel_tol=0.02),
}


def INSTANCE_ROWS(a, era, r):
    """Rows for resources.json: [(label, (lq_lo, lq_hi), (t_lo, t_hi), {extra}), ...]."""
    if era == "2028":
        i = r.intermediates
        return [(f"log det W, free staggered V={i['L']}^{int(a.D.lo)} (Hadamard test)", r.lq, r.hard_ops,
                 {"volume": f"{i['L']}^{int(a.D.lo)}", "observable": "Tr log W",
                  "eps_l": r.epsilon_l[0], "shots": i["shots"],
                  "lq_qme_route": i["lq_qme_route"], "degree": i["d"]})]
    if era == "2033":
        i = r.intermediates
        return [("Tr M^-1, SU(3) staggered V=8^3x16, kappa~1e2 (success-probability readout)", r.lq, r.hard_ops,
                 {"volume": "8^3x16", "observable": "Tr M^-1", "eps_l": r.epsilon_l[0],
                  "shots_first": i["shots_first"], "shots_campaign": i["shots_campaign"],
                  "n_cfg_campaign": i["n_cfg_campaign"]}),
                ("Banks-Casher mode number, V=8^3x16, a=0.15 fm (demonstration)", r.lq, i["bc_t_per_shot"],
                 {"volume": "8^3x16", "observable": "nu(Lambda)/V_F", "eps_l": a.p_f.lo / i["bc_t_per_shot"][1],
                  "shots": i["bc_shots"]})]
    if era == "codesign":
        return [("Tr M^-1, V=24^3x48, kappa~1e3, d_inv~5e3 (before reduction)", r.lq, r.hard_ops,
                 {"volume": "24^3x48", "observable": "disconnected HVP", "eps_l": a.eps_l_codesign.lo})]
    return []
