"""Ch. 5 — QGP transport coefficients. Reproduces the LQ and hard-op numbers of
applications/app03_qgp_transport.tex.

Two boxes: a 2028 pure-gauge SU(2)->2T benchmark on a 4^2 (2+1D) lattice, and a
2033 Sigma(216x3) + 3-staggered-field target on 3^3 at H_I (s648 author ruling, 2026-10-05; was Sigma(72x3) at
H_KS; the rulings2 g^2 = 3.5 / a_s = 0.69 fm instance is retired and the box is back at a_s = 0.2 fm). Both count qubits with
Eq. (Nq) at app03:68,  N_q = 3 L^d q_G + N_c N_stag L^d + N_anc, with q_G the
register of the compiled primitive set (ruling R1, 2026-09-28: 2T 5; Sigma(216x3) 11 since s648, Sigma(72x3) 9 before),
and count T-gates as  (links) x (T per link per Trotter step) x (Trotter steps),
the per-link-step cost being the sum of the cited papers' magnetic and electric
terms at H_KS from the papers' tables (ruling R2, 2026-09-28; PAPER_COSTS.md), priced since E20 (r18, 2026-10-01)
at 7 T per Toffoli and the full fit 1.15 log2(1/eps) + 9.2 per rotation (n_rot = the paper's log coefficient / 1.15;
numbers below are at E20, r18 Ch. 5 apply), eps set per circuit
by ruling R-TOL (eps_rot = sqrt(1e-2 / N_rot), common.eps_rot_for; 2028 7.628e-4, 2033 reference shot 7.99e-6).
The 2033 fermion hop and mass (r17, 2026-10-01) are priced from the Fermion_Primitives draft's gate tables at
the full fit 1.15 log2(1/eps) + 9.2, at the same circuit eps. Since r19 (author ruling E21 (1), 2026-10-01) the
hop holds the colour squish and parity flags once per link (share="link") with the frame undo kept (undo=True);
the r17 per-application recompute (share="draft") is a conservative sensitivity.
Every primitive cost comes from groups.py (magnetic_per_link and the primitives
table); nothing is retyped here.

INPUTS AND SOURCES
  L_2028=4, dim=2, a_t=0.1 fm/c, a_s=0.2 fm            app03:111 (box)
  Tc = 1.12 sqrt(sigma), sqrt(sigma)=440 MeV, T=1.5Tc     app03:112; arxiv_0803_2128
  2T: |G|=24, 5 qubits/link                              GROUPS["2T"].link_qubits (arxiv_2208_12309)
  2T magnetic term per link-step (H_KS, d=2)             groups.magnetic_per_link: 6 U_mul + 3 U_inv + U_Tr/2
                                                          = 1.124e3 T at eps=7.63e-4 (arxiv_2208_12309 tab:tgatecost)
  2T electric term per link-step                         2 x U_FFT (arxiv_2408_00075, Murairi et al.) = 1.969e3 T
  6 second-order Trotter steps at dt = a_t                app03:74,122
  1-2 Hadamard ancilla + 16 workspace (R10, r25)          app03 2028 box; factory.json Ch. 5
  shots: 1e3/timeslice x 7 timeslices t = 0..6 a_t        app03 2028 box (r25: each slice at its own depth, R4;
                                                          first result t = 0 and one step, R1)
  t_gate_s = 1 us per T = 10 us logical cycle / 10 factories   app03:78
  per-shot overhead ~0.1 ms (init, readout, decode)      report rule (Ch. 1, R9); 0.2 s of the 2028 first
                                                          result's 94.9 s, 115 s of the 2033 campaign at tau = 1/k
  L_2033=3, a_s=0.2 fm, a_t=0.1 fm/c                     app03:135
  Sigma(72x3): |G|=216, 9 qubits/link                    GROUPS["S72x3"].link_qubits (arxiv_2511_17437 :150)
  Sigma(72x3) magnetic term per link-step (H_KS, d=3)    groups.magnetic_per_link: 12 U_mul + 6 U_inv + U_Tr
                                                          = 1.888e4 T at eps=8.68e-6 (arxiv_2511_17437 tab:tgatecost)
  Sigma(72x3) electric term per link-step                2 x U_FFT LOWER BOUND (532 + 448 rot, arxiv_2511_17437 :905)
                                                          + U_phi (256 rot) = 3.39e4 T at the full fit; no FFT exists, a floor
  N_stag=3, N_c=3                                        app03:72,138
  staggered hop per link-step (FP draft, r17; E21)       groups.hop_link_cost("S72x3", 3): 3,672 Toffoli + 360 T + 2,933
                                                          rotations at the full fit = 1.098e5 T at eps=8.68e-6 (unpublished
                                                          gate counts, FermionPrimitives_unpub; squish held per link, E21 (1);
                                                          frame undo ESTIMATED here: FP draws V_C^dag but does not count it)
  staggered mass per site-step (r17, option M-B)          groups.staggered_mass_site(3, 3): one HWP(9) = 4 rot + 7 Toffoli
  dt = a_t/10 = 0.01 fm/c; t_max = 4 fm/c is the reference shot (per-step prices); each priced shot runs only to
                                                          its own fit time (r25, R4)
  100-200 ancilla                                        app03:144
  THERMAL STATE BY QUENCH (author ruling 2026-10-05, replaces E-rho-OQ at mu_B = 0 and the KMS Gibbs sampler at
    mu_B > 0): electric vacuum (or a Fock state in the target N_B sector), energy tuned by a coupling ramp of one
    a_s (quench_ramp_a_s, priced as ceil(a_s / dt) steps of the same Trotter circuit, a floor), then evolved for
    t_th = c / T with c = 2 (small-lattice toys relaxed) to 2 pi (Hayata-Hidaka arxiv_2011_09814) at the a_t/10
    grid step, then the fit evolution. One circuit per shot (prep + fit) at its own R-TOL eps. ETH / typicality:
    arxiv_cond-mat_9403051, arxiv_0708_1324, arxiv_1509_06411, arxiv_0902_0927, arxiv_2303_14264, arxiv_2308_16202.
    2033 (temperature ruling 2026-10-05, T = 1.5 T_c^lat = 225-300 MeV): 1/T = 0.658-0.877 fm/c = 65.8-87.7 steps;
    ramp 20 + thermalization 132 (c = 2, 300 MeV) to 552 (c = 2 pi, 225 MeV) = 152-572 prep steps (5.56e9-2.12e10 T);
    [at the retired 160 MeV: 123.3 steps, 267-795 prep steps, 9.80e9-2.95e10 T];
    2028: 1/T = 0.2669 fm/c, ramp 8 + 22-68 = 30-76 prep steps (quench_*_2028 records; the 2028 box stays non-thermal).
    Same circuit at mu_B = 0 and mu_B > 0 (the Fock-state sector fixes N_B), so one thermal route, no shot penalty.
  BE trims 10x                                           app03 route (a), budget split (Stated; BE not synthesized)
  r25 shots: Fisher, two fit points at Gamma t = 0.5, 1.8 (shot audit design B), Cbar = 0.16 (Stated), 30% (first
    result, R1) / 20% (campaign) on Gamma(k_min); campaign 4 x 4 (T, mu_B) x k_min x 2 spacings; decay band
    tau = 1/k_min = 0.0955 fm/c (T1, 2026-10-05; was the hydrodynamic 0.049) to 1/(2 pi T) = 0.140 fm/c at 225 MeV
  r23 reference shots: 1e3/grid-pt x (4x4 (T,mu_B) x 3 k x 2 t), records only   app03 r23 box
  2 lattice spacings in the 5-yr campaign                app03:197
  ONE machine, serial (E27, H. Lamm 2026-10-02)          every wall time = shots x (T x 1 us + 0.1 ms) on a single
                                                          machine; no machine count anywhere. Over-horizon campaigns
                                                          state the algorithmic reduction needed (serial yr / 5)
                                                          and what fits in 5 yr (apply_log/r23_ch05.md)
  0.1 expected faults per shot -> eps_l = 0.1/N          TRACKED_CHANGES.md ruling R3 (2026-09-28); app03:149
  RFI 2033 ceiling ~1e9 hard ops, eps_l 1e-8             DOE_RFI_2026; app03:149
  T = 1.5 T_c^lat, T_c^lat = 150-200 MeV ASSUMED        app03:42,135 (author ruling 2026-10-05: the simulated theory
    (Tc_lat_2033_MeV, T_over_Tc_2033) -> T = 225-300 MeV   has its own crossover, measured in situ; T = 160 MeV was
                                                          physical QCD. Bracket: HotQCD T_c^0 = 132 MeV arxiv_1903_04801,
                                                          physical 156.5 MeV arxiv_1812_08235; heavier pions raise it,
                                                          extra tastes lower it). LT = 0.68-0.91 (< 1: a small hot box),
                                                          k_min = 6.9-9.2 T. The cheap corner of every band sits at T hi.
  T_c^lat diagnostic (tc_diag_*)                         10 ramp energies x 100 prep-only shots (condensate read out in
                                                          the computational basis), 0.18-0.67 yr, ~6% of the first result
  m_pi = 300-400 MeV, m_pi L >= 4 target, 12^3 extension app03:51,77
  Sigma(360x3): 12 qubits/link compiled, 11 dense          QC/discrete_gates/s1080/paper_s1080.tex:95; app03:162
  $47.454M/yr NP Heavy Ion research x 5% -> "$2.4M/yr", 5-yr campaign "$12M"   app03:186 (utility box; not a resource number)
WHAT IS NOT DERIVED HERE
  Synthesis tolerance per rotation: ruling R-TOL (H. Lamm 2026-09-29), total synthesis error 1e-2 per
    shot, each circuit's eps_rot = sqrt(1e-2 / N_rot) with N_rot its own synthesized-rotation count
    (rotations_per_link_step x links x steps; 2028 17,184 -> 7.628e-4, 21.11 T at the full fit; 2033 evolution
    1.5658e8 -> 7.99e-6, 28.7 T; conjectured BE N_rot/10 -> 2.746e-5, the N/10
    reading flagged for the author). Replaces the fixed 1e-4 (T4, INDEX R2, E3).
    The papers' own fiducial is eps_T = 1e-8 total, linearly split; over this circuit's 1.326e8 rotations that is
    eps = 7.54e-17 per rotation, at which the Sigma(72x3) per-link-step is 3.36e5 (magnetic 1.92e4, electric floor 8.26e4,
    FP hop + mass 2.34e5; E21 sharing) and the mu_B=0 shot 1.09e10 (r18 verifier; earlier rounds priced eps = 1e-8 per
    rotation; at share='draft' 4.10e5 and 1.33e10)
    (app03:78 states this; intermediates *_at_papers_fiducial_eps).
  The Sigma(72x3) electric term is the paper's LOWER BOUND for a fast transform that has not
    been constructed (CircuitStatus.CONJECTURE in groups.py). The compiled dense U_F is
    ~270x larger (490x the magnetic term) at the 2033 eps; both are intermediates.
  Trotter step counts: 6 at 2028 (the full benchmark), so its deepest shot lands 5.9x over the 1e5-T envelope and the
    box says so (app03 2028 box; gap bullet). Ruling (b) 2026-10-01: the benchmark needs a >=5.9x cut and the named
    levers give 3-10x, so only the upper part of that range closes it. Since r25 the one-step first result fits.
    At 2033 each shot runs only to its own fit time (3-52 steps; the slow-end late shot takes 52 rather than its 36
    grid steps under the step rule, 2026-10-04); the 400-step 4 fm/c shot is a reference for per-step prices only.
  Step rule (2026-10-04, Claude's decision on the G7 item): the state-dependent second-order estimate at eps = 0.1
    sets the step (trotter_check, trotter_check_2028); 2028 runs six steps of a_t/4 to 0.15 fm/c and a one-step first
    result of 0.045 fm/c, at the same T as before.
  Quench thermalization time c / T with c in [2, 2 pi] (quench_c, Assumed): ETH at the chapter's couplings on 3^3,
    the relaxation of the momentum density within t_th and the thermal variance of the prepared state are ASSUMED;
    checked only by exact evolution on toys (Z2 up to 5x3, 2T up to 2x2, Z2 + staggered fermion N = 16; scratchpad
    quench/). The ramp length (one a_s) is a floor. The preparation runs at the grid step, not the state-dependent
    fit rule (it needs energy conservation, O(dt^2), not a faithful unitary); a coarser prep step (2-5x) is an
    UNPRICED lever (NEEDS_AUTHOR). The E-rho-OQ penalties (G8, computed 2026-10-05) and the 1e9 Gibbs prep are
    retired from the model with the methods (author ruling 2026-10-05; the computed values stay in NEEDS_AUTHOR).
  10x BE trim: Stated; the Sigma(72x3) block encoding does not exist yet.
  mu_B grid {0,200,400,600} MeV: Stated (ch05-mub-split, ruled option A); 1 of 4 values at mu_B=0 (4 and 12 points
    per spacing, k_min only); since the quench ruling every point costs the same per shot.
  2028 shots (1e3 per timeslice): no target precision is set (NEEDS_AUTHOR). The prose prices the first result at
    R1's 30% with an assumed relative fall-off delta = 0.1 (first_result_falloff_2028): 8.46e4 shots, 1.1 h.
  r23 record (retired from the text): 1e3 shots per grid point glossed as eps^-2 = 25 x ~40 (ch05-shots-eps2);
    N_grid = 96 = 4x4x3x2 with 24 of 96 points on the E-rho-OQ rail; these survive only as *_r23 / record keys.
  eta/s = 0.15: Stated (ch05-eta-s-tau); the hydrodynamic tau = T/((eta/s) k^2) = 0.069-0.092 fm/c at 225-300 MeV is a
    record (the chapter says it has no basis at this k).
  The 2033 staggered hop (app03:85 'Trotter evolution of the electric, magnetic, and fermion-hopping
    terms') is priced from the unpublished gate counts of the Fermion_Primitives draft (FermionPrimitives_unpub,
    r17 ruling (e), H. Lamm 2026-10-01: "adopt the fermion primitive numbers ... estimate or approximate the
    missing pieces"), through the shared groups.hop_link_cost("S72x3", n_stag=3) with its defaults
    (share="link", undo=True, mcx="mbu", phasing="hwp"; E21 (1), 2026-10-01, "make the switch"): per link per
    step the colour squish (2 x 1,219 Toffoli at MBU, compute + uncompute) and parity flags (892) computed ONCE
    and held through V_g(x), V_g(y), hop, V_g^dag(x), V_g^dag(y) (3,330 Toffoli); 4 colour-frame applications
    x 3 fields x 184 colour rotations (2 N_angles, E21 (3); 2,208 rotations), the V_g^dag half ESTIMATED HERE
    (FP draws V_C^dag in big_fig_group_agnostic.pdf, section_hamiltonian.tex:86-91, but its 2 C^G is two sites
    forward only, resources.tex:15,37-41); the hop squish + flags once per link (326 Toffoli); 2 x 10 classes
    x 3 colours x 3 fields = 180 controlled diagonalizers (360 T + 720 rotations); one HWP(18) phasing group
    (16 Toffoli + 5 rotations). In all 3,672 Toffoli + 360 T + 2,933 rotations per link-step, 1.098e5 T at the
    2033 eps (15.7x the retired R-HOP 6.97e3 at the same eps). The undo is 1,104 rotations, 3.15e4 T.
    Sensitivities at the same eps: share="draft" (the r17 headline, squish + parity per frame application and
    hop squish per field: 14,314 Toffoli) 1.843e5, shot 7.683e9; no undo (the draft's 2 C^G
    rotation layers, squish per link) 7.83e4; no undo at share="draft" 1.061e5.
  The staggered mass (r17, option M-B): one R_Z per colour-flavour copy per site (FP section_mass.tex:3-10), the
    9 copies on a site as one HWP(9) group, 4 rotations + 7 Toffoli = 163 T per site-step, 4.4e3 T per step
    (0.02% of the step). mu_B N_B commutes with H and folds into the same angle at no cost.
  The 2033 breakdown is ONE shot at a_s = 0.2 fm at the deepest fit point of the slow end (tau = 1/(2 pi T), c =
    2 pi): quench prep (795 steps) + fit evolution (52 steps) of the same primitives; the second spacing's late shot
    (1054 + 80 steps) is an intermediate (t_deepest_campaign).
  Which early transient the fit excludes (design A / B / C, a 7x spread in shots) and Cbar are NOT derived: design B
    and Cbar = 0.16 are carried (NEEDS_AUTHOR). The 2028 shot count (1e3 per timeslice) has no stated precision.
  The decay time at k_min ~ 7-9 T: the hydrodynamic 0.049 fm/c has no basis there; the band runs from the collisionless
    dephasing time 1/k (T1, 2026-10-05) to 1/(2 pi T). The second spacing is priced at its own step (T3, 2026-10-05).
WALL TIME AND FACTORIES
  t_gate_s = 1e-6 s per T, shot_overhead_s = 1e-4 s per shot, ONE machine, serial. Each shot is charged
    max(N_T t_gate_s, D_T x 10 us) + 0.1 ms at the HIGH depth schedule, so a shot with F* < 10 pays its depth.
  T-depth per shot from factory.json (Ch. 5 runs; schedules re-evaluated here at each shot's own t_rot by
    step_depth_2028 / step_depth_2033, which reproduce factory.json's 628 / 3882 (2028, 64 workspace) and 4.68e4 /
    6.05e5 (2033) per step), low/high; the quench prep is the same circuit, so its depth scales with its steps.
    F* = N_T / D_T: 2028 10.7-57 with 16 workspace qubits (R10; ~1.7 without); 2033 ~35 (high, every shot, quench included).
  2028. First result: t = 0 and one step, 2e3 shots, 94.9 s -> wall_first_result_s. Full benchmark: 7 timeslices,
    7e3 shots, 2059 s (34.3 min) -> wall_campaign_s. Exports are the full benchmark, shot-averaged.
  2033 (s648 Sigma(216x3), H_I, a_s = 0.2 fm; quench prep 2026-10-05). Two corners: (tau = 1/k, c = 2) and
    (tau = 1/(2 pi T), c = 2 pi), at T = 300 and 225 MeV (temperature ruling 2026-10-05). First result: Gamma(k_min),
    T = 1.5 T_c^lat, mu_B = 0, 30%: 1.60e4-1.53e4 shots, 3.1-10.7 yr -> wall_first_result_s (the other six corners
    3.5-11.0 yr, first_result_corners). Campaign: 20%, 16 points x 2 spacings (0.2 and 0.15 fm, each at its own
    step), 1.15e6-1.14e6 shots -> wall_campaign_s, 254-938 yr, 51-188x the 5-yr horizon, a co-design gap; 5 yr do not
    cover one 20% grid point (7.0-24 yr). Exports: the campaign at the (300 MeV, 1/k, c = 2) corner over both
    spacings, shot-averaged. (At the retired 160 MeV: first result 5.3-15.3 yr, campaign 435-1325 yr; before the
    quench: deepest 6.8e8-3.0e9 with 1e8 Gibbs prep on the mu_B > 0 rail, first result 107-283 d, campaign 23.9-75.7 yr.)
"""
from __future__ import annotations

import math
from dataclasses import dataclass

from estimates.common import (Assumed, Cited, Stated, Uncited, Tagged, Primitive,
                              CircuitStatus, Result, Published, ERAS,
                              t_per_rotation, eps_per_rotation, toffoli_t,
                              EPS_SYN, EPS_SYN_SRC, eps_rot_for, hwp_synth_rotations, hwp_ancilla,
                              depth_exports, REACTION_TIME_S, FACTORY_BASELINE, FACTORY_SRC)
from estimates.groups import (GROUPS, PRIMCOST, MAGNETIC, magnetic_per_link, electric_per_link,
                              rhop_link, fermion_hop_counts, hop_link_cost, staggered_mass_site,
                              FP_KEY, FP_DIR, FERMION_DRAFT)

ETH_KEYS = ("arxiv_cond-mat_9403051,arxiv_0708_1324,arxiv_1509_06411,arxiv_0902_0927,arxiv_2303_14264,"
            "arxiv_2308_16202,arxiv_2011_09814")

TEX = "app03"
HBARC_MEV_FM = 197.327          # PDG; used for k = 2 pi / L and tau conversions
SECONDS_PER_YEAR = 3.156e7
R3 = "TRACKED_CHANGES.md 'RULINGS, NEEDS_AUTHOR round A' R3 (2026-09-28)"
R25 = "apply_log/r25_ch05.md: author rulings R1-R10 (H. Lamm 2026-10-02)"


@dataclass(frozen=True)
class Assumptions:
    # ---- shared physical constants -------------------------------------------
    hbarc: Tagged = Cited(HBARC_MEV_FM, "PDG", "hbar c in MeV fm")
    seconds_per_year: Tagged = Cited(SECONDS_PER_YEAR, "PDG", "Julian year")

    # ---- per-link-step pricing (R2 tables, E20 pricing)-------------------------------
    hamiltonian: Tagged = Stated("KS", f"{TEX}:26,78", "H_KS; tab:primcost multiplicities for H_KS (the 2028 2T box)")
    hamiltonian_2033: Tagged = Stated("I", f"{TEX}:eq:Hks; arxiv_2203_02823",
                                      "s648 author ruling (H. Lamm 2026-10-05): the 2033 box is priced at the improved "
                                      "Hamiltonian H_I (tab:primcost H_I row: 4 U_F, 2 U_Ph, 3(d-1)/2 U_Tr, 2 + 11(d-1) "
                                      "U_inv, 4 + 26(d-1) U_mul per link-step; Gustafson_S648_inprep). Was H_KS with "
                                      "Sigma(72x3). 'KS' restores the H_KS multiplicities for the new group (sensitivity)")
    eps_syn: Tagged = Cited(EPS_SYN, EPS_SYN_SRC,
                            "ruling R-TOL: total synthesis error per shot, report-wide; each circuit sets its own "
                            "per-rotation tolerance eps_rot = sqrt(eps_syn / N_rot) (common.eps_rot_for, randomized "
                            "synthesis), N_rot = synthesized rotations in one shot of that circuit. Replaces the fixed "
                            "1e-4 per rotation (R-TOL). The papers' own fiducial is eps_T = 1e-8 total with a linear split "
                            "(PAPER_COSTS.md sec. 3-4); app03:78 prices that as a sensitivity")
    papers_fiducial_eps_T: Tagged = Cited(1e-8, "arxiv_2511_17437,arxiv_2405_05973,arxiv_2208_12309",
                                          "papers' benchmark: eps_T = 1e-8 total synthesis error, split linearly over all "
                                          "rotations (PAPER_COSTS.md sec. 4); used only to price the sensitivity at app03:78")
    synthesis_convention: Tagged = Cited(
        "rus", "arxiv_2208_12309,arxiv_2511_17437,arxiv_2408_00075; apply_log/r18_core.md",
        "E20 ruling (2026-10-01, supersedes R2's papers' convention): the papers' tables with every rotation at the "
        "full fit 1.15 log2(1/eps) + 9.2 (n_rot = log coefficient / 1.15), 7 T/Toffoli. The papers' slope-only "
        "price is kept as the *_papers record")
    toffoli_convention: Tagged = Cited("textbook", "arxiv_2208_12309,arxiv_2511_17437", "7 T per Toffoli in every cited paper")

    # ---- 2028 box (app03:107-129) --------------------------------------------
    L_2028: Tagged = Stated(4, f"{TEX}:111", "Volume: V = 4^2 transverse (2+1D)")
    dim_2028: Tagged = Stated(2, f"{TEX}:111", "2+1D: two spatial dimensions")
    a_t_fm: Tagged = Stated(0.1, f"{TEX}:111,135", "a_t = 0.1 fm/c in both boxes")
    a_s_fm: Tagged = Stated(0.2, f"{TEX}:111,135", "a_s = 0.2 fm in both boxes")
    group_2028: Tagged = Stated("2T", f"{TEX}:114", "SU(2) -> 2T binary tetrahedral, |G|=24")
    sqrt_sigma_MeV: Tagged = Stated(440, f"{TEX}:112", "sqrt(sigma) = 440 MeV")
    tc_over_sqrt_sigma: Tagged = Cited(1.12, "arxiv_0803_2128", "T_c ~ 1.12 sqrt(sigma) for 2+1D SU(2)")
    T_over_Tc_2028: Tagged = Stated(1.5, f"{TEX}:112", "T = 1.5 T_c ~ 740 MeV")
    trotter_steps_2028: Tagged = Stated(6, f"{TEX}:122", "x 6 2nd-order Trotter steps")
    dt_over_at_2028: Tagged = Assumed(0.25, "step rule (2026-10-04, Claude's decision on the G7 item, recorded in "
                                            "NEEDS_AUTHOR): the state-dependent second-order estimate (trotter_check_2028) "
                                            "at eps = 0.1 sets the step; six steps of a_t/4 reach t_max = 0.15 fm/c at "
                                            "0.083 (was Delta t = a_t, 6 steps to 0.6 fm/c, free-field error 0.78). Same "
                                            "T per shot: the gate count of a step does not depend on Delta t")
    dt_first_2028_fm: Tagged = Assumed(0.045, "step rule: the first result is t = 0 and one step of 0.045 fm/c, inside "
                                              "the rule (0.081; the largest one-step Delta t at 0.1 is 0.048 fm/c); one "
                                              "step costs the same at any Delta t")
    dt_over_at_2033: Tagged = Stated(0.1, f"{TEX}:74,147", "the 2033 box takes Delta t = a_t/10: the grid of the fit "
                                     "times and the step of every shot that meets the step rule at it; the slow-end late "
                                     "shot takes more steps (trotter_rule_2033, own_depth_2033)")
    trotter_rule_2033: Tagged = Assumed("state", "STEP RULE (2026-10-04, Claude's decision on the G7 item, NEEDS_AUTHOR "
                                                 "closed): each fit shot runs max(t / (a_t/10), N_sd) steps, N_sd = ceil("
                                                 "sqrt(Lambda_sd t^3 / eps)) the fewest steps at which the state-dependent "
                                                 "second-order estimate (connected variance, Ch. 9) meets eps = 0.1. Only the "
                                                 "slow-end late shot binds: 36 -> 44 steps (0.147 -> 0.099). In a mixed state the "
                                                 "estimate is an upper estimate of the part that grows with t; at T a_s = 0.16 "
                                                 "it is almost all zero-point (0.03 of 32.8 grows), so 44 is conservative "
                                                 "(free field alone: 36). 'grid' restores "
                                                 "the a_t/10 step everywhere (record)")
    anc_2028: Tagged = Stated((1, 2), f"{TEX}:120", "+ 1--2 Hadamard ancilla")
    shots_per_timeslice: Tagged = Stated(1e3, f"{TEX}:121", "10^3/timeslice; the target precision on the early fall-off is "
                                                          "not stated (NEEDS_AUTHOR), so the count stays as printed")
    t_gate_s: Tagged = Stated(1e-6, f"{TEX}:78,192", "1 us per T gate, i.e. ~10 "
                                                       "factories at the 10 us logical cycle; R9 unified name (was gate_time_us)")
    shot_overhead_s: Tagged = Assumed(1e-4, "per-shot overhead ~0.1 ms (register init, final readout, decode), report "
                                            "rule; as Ch. 8 shot_overhead_s. E27: every wall time is serial on ONE "
                                            "machine, shots x (T x 1 us + t0); R9 unified name",
                                      "main-overview.tex 'A reference machine'; ch08_baryogenesis.py")
    rfi_2028_t: Tagged = Cited(1e5, "DOE_RFI_2026", "first-generation ~1e5 T/Toffoli hard ops (app03:74,78,123)")
    per_step_cut_levers: Tagged = Stated((3, 10), f"{TEX}:78,123,180",
                                         "the named per-step levers give 3-10x (ruling (b) 2026-10-01: the 2028 box needs "
                                         ">=5.9x at E20, so only the upper part of the range closes it)")

    # ---- 2028 tiers and workspace (r25: R1, R4, R10) ----------------------------
    first_result_steps_2028: Tagged = Assumed(1, "R1/R4: the 2028 first result is the fall-off from t = 0 to one step "
                                                 "(two timeslices), each shot evolved to its own time; one step at its own "
                                                 "R-TOL eps is 9.47e4 T, inside 1e5 (shot audit Ch. 5)", R25)
    first_result_falloff_2028: Tagged = Assumed(0.1, "NOT A BOX NUMBER. Relative fall-off delta of the normalized "
                                                    "correlator between t = 0 and one step (0.045 fm/c since the step rule; for an even "
                                                    "correlator delta ~ t^2, so 0.1 at a_t would be ~0.02 here), used only to price the "
                                                    "first result at 30% (R1) in the prose: shot-audit formula "
                                                    "N = 2 (1 - c0^2) / (c0^2 r^2 delta^2), c0 = cbar. delta is open "
                                                    "(NEEDS_AUTHOR, 2028 shot count)", "shot_audit.json Ch. 5")
    workspace_2028: Tagged = Assumed(16, "R10: primitive workspace on top of the 160 link + 1-2 Hadamard qubits, so 8 U_FFT "
                                         "(2 ancilla each) run at once and F* >= 10 (10.7-10.9 at the high schedule). "
                                         "factory.json Ch. 5 2028 row: +2 workspace gives F* = 1.5, +64 gives 25.5",
                                     f"{R25}; factory.json")

    # ---- 2033 box (app03:131-171) --------------------------------------------
    L_2033: Tagged = Stated(3, f"{TEX}:135", "Volume: V = 3^3")
    dim_2033: Tagged = Stated(3, f"{TEX}:135", "3+1D")
    group_2033: Tagged = Stated("S216x3", f"{TEX}:137", "s648 ruling (2026-10-05): Sigma(216x3) subset SU(3), |G| = 648, "
                                                         "11 qubits/link in the draft's qubit encoding (Gustafson_S648_inprep, "
                                                         "in preparation). Was Sigma(72x3), |G| = 216, 9 qubits/link")
    hop_group_2033: Tagged = Assumed("S72x3", "s648: no hop circuit exists for Sigma(216x3) (the Fermion_Primitives draft "
                                              "covers 2T, 2O, Sigma(36x3), Sigma(72x3)). The 2033 hop is priced at the "
                                              "Sigma(72x3) counts as a PLACEHOLDER, not a bound: Sigma(216x3) has 24 classes "
                                              "against 16 and a wider register, so its colour frame is likely larger. "
                                              "NEEDS_AUTHOR")
    # s648 (author ruling 2026-10-05, (2)): freezing from the literature, not Monte Carlo
    beta_f_S72x3: Tagged = Cited(3.18, "arxiv_2511_17437", "isotropic Wilson 3+1D freezing beta_f = 3.18(3); g^2_f = 6 / "
                                                            "beta_f = 1.89")
    beta_f_S216x3: Tagged = Cited(3.80, "Gustafson_S648_inprep", "isotropic Wilson 3+1D freezing beta_f = 3.80(5) (draft "
                                                                  "tab:subgroups :190); g^2_f = 1.58")
    nonfreezing_HI_assumed: Tagged = Assumed(True, "s648 ruling (2): a non-freezing H_I trajectory for Sigma(216x3) down to "
                                                    "the box spacings (0.2 fm and the second spacing 0.15 fm) is ASSUMED by "
                                                    "analogy with S(1080), whose Wilson action plus one extra single-plaquette "
                                                    "trace term reached a ~ 0.08 fm and reproduced T_c (arxiv_1906_11213; an "
                                                    "adjoint-term action failed for S(648), reviewed there). H_I (Symanzik "
                                                    "O(a^2) improvement) is a different action, so the analogy is weak; H_I "
                                                    "plus such a plaquette term (unpriced) is the fallback. Not established. "
                                                    "Every crystal-like SU(3) subgroup freezes before the scaling regime with "
                                                    "the Wilson / KS action (arxiv_2112_08482, arxiv_2511_17437, draft), and "
                                                    "the Hamiltonian (anisotropic) limit moves freezing to stronger coupling (Gustafson_S648_inprep, in preparation; author ruling 2026-10-05)",
                                             "arxiv_1906_11213,arxiv_2203_02330")
    n_stag: Tagged = Stated(3, f"{TEX}:138", "3 staggered fields (N_f = 4+2 tastes)")
    n_c: Tagged = Stated(3, f"{TEX}:72", "N_c = 3 colors per staggered field per site")
    hop_rule: Tagged = Stated("FP-draft", f"{TEX}:85; {FP_KEY} ({FP_DIR}); apply_log/r17_core.md",
                              "the 2033 staggered hop, per link per Trotter step, by the shared groups.hop_link_cost "
                              "from the unpublished Fermion_Primitives gate counts (r17 ruling (e), H. Lamm 2026-10-01), "
                              "defaults share='link' (E21 (1): squish + parity held per link), undo=True (V_g^dag ESTIMATED: drawn in "
                              "FP, not counted), mcx='mbu', phasing='hwp'; share='draft' is a sensitivity")
    fermion_synthesis: Tagged = Cited("rus", f"{FP_KEY} section_resources.tex:7-10; apply_log/r17_core.md rule 11",
                                      "fermion hop and mass rotations at the full fit 1.15 log2(1/eps) + 9.2, as the draft "
                                      "prices them; never the slope-only variant for new prices (r17)")
    mass_rule: Tagged = Stated("HWP-site", f"{FP_KEY} section_mass.tex:3-10; apply_log/r17_core.md rule 12",
                               "staggered mass: the N_c N_stag equal-angle R_Z per site as one HWP group (option M-B); "
                               "mu_B N_B folds into the same angle")
    anc_2033: Tagged = Stated((100, 200), f"{TEX}:144", "+ 100--200 ancilla")
    t_max_2033_fm: Tagged = Stated(4.0, f"{TEX}:53,141,147", "t_max = 4 fm/c")
    # ---- 2033 temperature (author ruling 2026-10-05, Ch. 5 2033 temperature) ------------------------------------
    # T = 160 MeV was borrowed from physical QCD. The simulated theory (Sigma(216x3)/H_I, 3 unrooted staggered fields =
    # 4+2 tastes, m_pi = 300-400 MeV, L = 0.6 fm) has its own crossover T_c^lat, unknown, to be measured in situ
    # (chiral condensate / susceptibility from prep-only diagnostic shots, tc_diag_*). The 2033 target is T = 1.5 T_c^lat,
    # the 2028 convention; for pricing T_c^lat = 150-200 MeV, i.e. T = 225-300 MeV, carried as a band through every
    # T-dependent number (1/T steps, t_th, tau = 1/(2 pi T), fit steps, eps, walls, k_min / T, LT).
    Tc_lat_2033_MeV: Tagged = Assumed((150.0, 200.0),
                                      "T_c^lat of the simulated lattice theory, NOT KNOWN; bracketed by the 2+1-flavour "
                                      "QCD crossover: HotQCD chiral-limit T_c^0 = 132(+3)(-6) MeV (arxiv_1903_04801, Ding "
                                      "et al. PRL 123, 062002) and physical-mass 156.5(1.5) MeV (arxiv_1812_08235, HotQCD "
                                      "PLB 795, 15); T_c rises roughly linearly with m_pi, so the heavier pions (300-400 "
                                      "MeV) push it up, while the extra light tastes (4+2) pull it down. Measured in situ "
                                      "(tc_diag_*) before the campaign", "arxiv_1903_04801,arxiv_1812_08235")
    T_over_Tc_2033: Tagged = Stated(1.5, f"{TEX}:42,135", "T = 1.5 T_c^lat, the same convention as the 2028 benchmark "
                                                         "(author ruling 2026-10-05; was T = 160 MeV, physical QCD)")
    tc_diag_points: Tagged = Assumed(10, "in-situ T_c^lat diagnostic: ramp energies scanned to locate the crossover "
                                         "(the staggered condensate drop / susceptibility peak); each point is the quench "
                                         "preparation alone, read out in the computational basis (the mass term is diagonal)")
    tc_diag_shots_per_point: Tagged = Assumed(100, "shots per diagnostic energy: one shot reads all 243 fermion modes, so "
                                                   "~1e2 shots fix the site-averaged condensate to a few percent and its "
                                                   "variance (the susceptibility) to ~20%; an order of magnitude, not derived")
    eta_s_assumed: Tagged = Stated(0.15, f"{TEX}:79", "eta/s ~ 0.15 in the hydrodynamic tau = T/((eta/s) k^2), a record "
                                                     "(the chapter says it has no basis at this k)")
    eta_s_hydro_range: Tagged = Cited((0.08, 0.2), "JETSCAPE_Bayesian_2020",
                                      "app03:7: Bayesian fits pin eta/s in 0.08--0.2")
    # ---- thermal state by quench (author ruling 2026-10-05; route (c)) --------------------------------------
    quench_c: Tagged = Assumed((2.0, 2 * math.pi),
                               "thermalization time t_th = c / T of the quench route: a simple gauge-invariant state "
                               "(electric vacuum, or a Fock state in the target N_B sector) is evolved under the same "
                               "Trotterized Hamiltonian until its correlators match the microcanonical ensemble (ETH, "
                               "dynamical typicality; canonical up to O(1/V)). c = 2: the exact toys (Z2 up to 5x3, 2T "
                               "up to 2x2, Z2 + staggered fermion N = 16) relax within 0.6-2.5 / T; c = 2 pi: the "
                               "3+1D SU(2) small-lattice study of Hayata-Hidaka. ETH at the chapter's couplings on 3^3, "
                               "and the relaxation of the momentum density (a conserved density, which needs the "
                               "average over random initial states) within t_th, are ASSUMED. Priced at the a_t/10 "
                               "grid step (energy conservation, not a faithful unitary; a coarser prep step is an "
                               "unpriced lever, NEEDS_AUTHOR)", ETH_KEYS)
    quench_ramp_a_s: Tagged = Assumed(1.0, "energy tuning of the initial state by a coupling ramp of one a_s (lattice "
                                           "time unit), priced as ceil(a_s / dt) steps of the same Trotter circuit, a "
                                           "FLOOR (a single magnetic layer is the alternative). On 2T 2x2 a ramp of a "
                                           "few a_s moves T_eff from 0.80 to 0.30 (scratchpad quench/res "
                                           "out_2tramp_2x2)", "scratchpad quench/run/g2t_ramp2.py")
    be_trim: Tagged = Stated(10, f"{TEX}:91,147", "qubitized BE trims 10x; Sigma(72x3) BE not yet synthesized")
    ceiling_2033_t: Tagged = Cited(1e9, "DOE_RFI_2026", "'the RFI's ~1e9-hard-op 2033 envelope' (app03 budget split (classical-sample rail),149); R3: the ceiling")
    rfi_eps_l: Tagged = Cited(1e-8, "DOE_RFI_2026", "the RFI's logical error rate; eps_l is printed as a requirement instead")
    faults_per_shot: Tagged = Cited(0.1, R3, "0.1 expected faults per shot, report-wide: eps_l = 0.1/N (app03:149)")
    shots_per_gridpt: Tagged = Stated(1e3, "app03 r23 box (retired r25)",
                                      "RECORD: the r23 10^3/grid-pt (eps^-2 = 25 x ~40) of the 4 fm/c reference campaign; "
                                      "r25 shots come from fisher_shots")
    grid_T: Tagged = Stated(4, f"{TEX}:44,144", "4 x 4 in (T, mu_B)")
    grid_muB: Tagged = Stated(4, f"{TEX}:44,144", "4 x 4 in (T, mu_B)")
    grid_k: Tagged = Stated(3, "app03 r23 box (retired r25)", "RECORD: x 3 k of the 4 fm/c reference campaign; r25: campaign_k")
    grid_t: Tagged = Stated(2, "app03 r23 box (retired r25)", "RECORD: x 2 t-windows; r25: the two fit points of fit_x")
    n_grid_prose: Tagged = Stated(96, "app03 r23 box (retired r25)", "RECORD: N_grid = 96 of the reference campaign")
    mu_b_zero_fraction: Tagged = Stated(0.25, f"{TEX}:78,144",
                                        "mu_B in {0,200,400,600} MeV: 1 of 4 at mu_B=0; r25 campaign '4 points at "
                                        "mu_B=0 and 12 at mu_B>0' (ch05-mub-split option A). Since the quench ruling "
                                        "(2026-10-05) every point costs the same per shot (one thermal route)")
    n_spacings: Tagged = Stated(2, f"{TEX}:78,187", "x 2 lattice spacings")
    second_spacing_a_s_fm: Tagged = Assumed(0.15, "referee T3 (verifier 2026-10-04): an example finer spacing for the "
                                                  "campaign's second spacing ('e.g. a_s = 0.15 fm'); the box prices that "
                                                  "spacing at the first spacing's per-shot cost (same 3^3 lattice, same "
                                                  "Delta t = 0.01 fm/c). NEEDS_AUTHOR: choose it.")
    # ---- referee G7 (2026-10-04): Trotter error at the chosen Delta t = a_t/10 (trotter_check) ----
    trotter_eps: Tagged = Assumed(0.1, "G7: target Trotter error of a shot, ||(U - e^{i theta} S^N) rho^{1/2}||_2 in the "
                                       "simulated state, the measure and value of Ch. 9 (app06 'at eps ~ 0.1')", "app06")
    trotter_g2: Tagged = Assumed(1.0, "G7: coupling g^2 for the worst-case bound only (not tuned for Sigma(216x3) / H_I at "
                                      "a_s = 0.2 fm; the bound is minimized over g^2 as a check, trotter_check)")
    trotter_dispersion_2033: Tagged = Assumed("I", "s648: the step rule's free-field model for H_I is the tree-level "
                                                   "Symanzik dispersion w^2 = sum_i (khat_i^2 + khat_i^4 / 12) (the same "
                                                   "reading as Ch. 6, free_gauge_frequencies_h); H_I's two-link electric "
                                                   "term, which lowers w at large k, is left out, so the estimate is on the "
                                                   "high side. In the free field the improved electric and magnetic terms "
                                                   "each still commute among themselves, so the two-group second-order "
                                                   "formula applies mode by mode. 'KS' restores the H_KS dispersion (record)")
    sigma72_emax_over_fund: Tagged = Cited(2.0, "arxiv_2511_17437", "tab:eham: the Sigma(72x3) electric eigenvalues "
                                           "f(rho) run to 72 against 36 for the faithful 3; with the 3 normalized to the "
                                           "SU(3) Casimir 4/3 the top of the electric spectrum is (g^2/2) 8/3 (G7). "
                                           "NORMALIZATION (2026-10-04, Claude's decision, NEEDS_AUTHOR closed; stated at "
                                           "app03 eq:Hks): E^a E^a -> (4/3) f(rho) / f(3). f is the Cayley-graph Laplacian "
                                           "f = |Gamma| - Re sum_Gamma chi/dim (arxiv_2511_17437 eq:electric, Gamma from the "
                                           "transfer matrix); Gamma is the 54 elements of trace 1 (checked: the group built "
                                           "from the paper's generators C, E, V, X has order 216, f(3) = 2|Gamma|/3 = 36, "
                                           "f(8) = |Gamma| - sum(|Tr|^2 - 1)/8 = 54; test_eham_normalization), far "
                                           "from the identity, so f is not proportional to the Casimir and f(8)/f(3) = 1.5 "
                                           "is the group's, not a misprint. Matching the fundamental keeps the strong-coupling "
                                           "flux energy (g^2/2a) C_F per link that the workflow's string-tension tuning fixes; "
                                           "matching the adjoint (C_max = 4) would raise the worst-case bound 1.5x")
    trotter_n_adj: Tagged = Assumed(8, "G7: weak-coupling (free-field) evaluation treats Sigma(72x3) as SU(3): 8 colour "
                                       "copies of each transverse gauge oscillator; a discrete group has no such expansion")
    campaign_horizon_yr: Tagged = Stated(5, f"{TEX}:196", "Campaign horizon 5 years")
    # ---- r25: each shot at its own depth (R4), two tiers (R1), Gibbs cut required (R5) ----------------------
    cbar: Tagged = Stated(0.16, f"{TEX}:78", "normalized correlator Cbar ~ 0.16 of the +/-1 Hadamard-test estimator "
                                             "(the old variance factor ~40 = 1/Cbar^2); not derived (NEEDS_AUTHOR)")
    fit_x: Tagged = Assumed((0.5, 1.8), "fit points at Gamma t = 0.5 and 1.8 (shot audit design B): the fit excludes "
                                        "the early transient Gamma t < 0.5, since the symmetrized C(t) is even in t; the "
                                        "points are both rounded UP onto the Delta t grid (Gamma t = 0.52, 1.88 at tau = 1/k = 0.0955; 0.50, 1.86 at 1/(2 pi T) = 0.140, 225 MeV). "
                                        "Designs A (from t = 0, 0.37x the shots) and C (Gamma t >= 1, 2.7x) are the "
                                        "sensitivities (NEEDS_AUTHOR: which transient the fit excludes)", "shot_audit.json Ch. 5")
    fit_frac_early: Tagged = Assumed(0.2, "share of the shots at the early fit point (design B)", "shot_audit.json Ch. 5")
    target_first: Tagged = Assumed(0.3, "R1: the first result, Gamma(k_min) at T = 1.5 T_c^lat, mu_B = 0, to 30% statistical", R25)
    target_campaign: Tagged = Stated(0.2, f"{TEX}:43", "the campaign: Gamma(k_min) to 20% statistical at every grid point")
    campaign_k: Tagged = Assumed(1, "R4 / shot audit: the campaign prices k_min only; a k -> 0 extrapolation has no "
                                    "meaning at L = 0.6 fm (k_min ~ 7-9 T, app03:74)", "shot_audit.json Ch. 5")
    tau_slow_2piT: Tagged = Assumed(1.0, "slow end of the decay band, tau = factor / (2 pi T) = 0.140 fm/c at 225 MeV "
                                         "(0.196 at the retired 160 MeV): at k ~ 7-9 T the hydrodynamic tau has no basis (app03:79), and "
                                         "1/(2 pi T) is the thermal decay scale of non-hydrodynamic modes (shot audit's 0.2)",
                                    "shot_audit.json Ch. 5")
    tau_fast_1k: Tagged = Assumed(1.0, "T1 (referee; 2026-10-05): fast end of the decay band = factor x 1/k_min, the "
                                       "collisionless (free-streaming) dephasing time of the transverse momentum density "
                                       "at k_min ~ 7-9 T (0.0955 fm/c on 3^3 at a_s = 0.2 fm). Replaces the hydrodynamic "
                                       "T/((eta/s) k^2) = 0.049 fm/c, which the chapter says has no basis at this k; that "
                                       "value stays a record (tau_decay_at_eta_s_assumed_fm)", "app03:79")
    second_spacing_own_step: Tagged = Assumed(True, "T3 (2026-10-05): the campaign's second spacing is priced at its own "
                                                    "step under the step rule (a_t scaled with a_s at fixed anisotropy, so "
                                                    "Delta t ||H|| is held), not at the first spacing's per-shot cost, which "
                                                    "let Delta t ||H|| grow by 4/3 and broke the chapter's own step rule. "
                                                    "False restores the old convention (record)", "app03:72,79")
    # G5 open item 4 (Claude's decision, 2026-10-04): RHIC operations end for EIC construction, so the base is the
    # NP Heavy Ion research line (FY 2026 CJ, Science/NP p. 220: FY24 and FY25 enacted $47,454K; FY26 request $37,004K)
    heavy_ion_research_musd_per_yr: Tagged = Cited(47.454, "DOE_NP_FY26_budget", "app03:186: '$47.5M/yr (FY24--FY25 enacted; $37.0M requested for FY26)'")
    utility_fraction: Tagged = Assumed(0.05, "app03:186: '5% of it gives ~$2.4M/yr ... The percentage is the stated, rescalable assumption'", f"{TEX}:186")
    m_pi_MeV: Tagged = Stated((300, 400), f"{TEX}:51,138", "m_pi = 300--400 MeV")
    m_pi_L_target: Tagged = Stated(4, f"{TEX}:77", "requires m_pi L >~ 4")
    L_ext: Tagged = Stated(12, f"{TEX}:77", "~12^3 lattice at a_s = 0.2 fm")
    L_ext_small: Tagged = Stated(4, f"{TEX}:79", "the 4^3-volume extension")
    L_ext_fm: Tagged = Stated(2.3, f"{TEX}:77", "L >~ 2.3 fm; k_min ~ 3.4 T is quoted at this L")
    sigma360_order: Tagged = Stated(1080, f"{TEX}:87", "Sigma(360x3): ceil(log2 1080) = 11 qubits/link dense")
    sigma360_qubits_compiled: Tagged = Cited(12, "QC/discrete_gates/s1080/paper_s1080.tex:95",
                                             "'This presentation requires ~12 qubits'; not in groups.py (no S360x3 entry)")
    lanl_lq: Tagged = Cited(36000, "Baertschi2024LANL", "36,000-LQ estimate for physical-mass QGP")

    def __post_init__(self):
        for name in ("L_2028", "L_2033", "trotter_steps_2028", "n_stag", "n_c"):
            if getattr(self, name).lo <= 0:
                raise ValueError(f"{name} must be positive")
        for name in ("anc_2028", "anc_2033", "m_pi_MeV", "eta_s_hydro_range", "quench_c", "Tc_lat_2033_MeV"):
            v = getattr(self, name)
            if not v.is_range or v.lo > v.hi:
                raise ValueError(f"{name} must be a (lo, hi) range")
        if self.Tc_lat_2033_MeV.lo <= 0 or self.T_over_Tc_2033.lo <= 0:
            raise ValueError("Tc_lat_2033_MeV and T_over_Tc_2033 must be positive")
        if self.tc_diag_points.lo < 0 or self.tc_diag_shots_per_point.lo < 0:
            raise ValueError("tc_diag_points and tc_diag_shots_per_point must be non-negative")
        if self.quench_c.lo <= 0 or self.quench_ramp_a_s.lo < 0:
            raise ValueError("quench_c must be positive and quench_ramp_a_s non-negative")
        if not (0.0 <= self.mu_b_zero_fraction.lo <= 1.0):
            raise ValueError("mu_b_zero_fraction must be in [0, 1]")
        for gf in ("group_2028", "group_2033"):
            gk = getattr(self, gf).value
            if gk not in GROUPS:
                raise ValueError("unknown gauge group")
            if not set(MAGNETIC) <= set(GROUPS[gk].primitives) or "U_FFT" not in GROUPS[gk].primitives:
                raise ValueError(f"{gf}={gk!r}: groups.py has no magnetic/FFT primitive table for it")
        if self.hamiltonian.value not in PRIMCOST or self.hamiltonian_2033.value not in PRIMCOST:
            raise ValueError("hamiltonian must be a tab:primcost key ('KS' | 'I')")
        if self.trotter_dispersion_2033.value not in ("KS", "I"):
            raise ValueError("trotter_dispersion_2033 must be 'KS' or 'I'")
        if self.fermion_synthesis.value not in ("rus", "watson"):
            raise ValueError("fermion_synthesis: full-fit models only (r17: never the slope-only variant)")
        if self.toffoli_convention.value not in ("textbook", "jones"):
            raise ValueError("toffoli_convention must be 'textbook' or 'jones'")
        if not (0.0 <= self.utility_fraction.lo <= 1.0):
            raise ValueError("utility_fraction must be in [0, 1]")
        if not (0.0 < self.faults_per_shot.lo <= 1.0):
            raise ValueError("faults_per_shot must be in (0, 1]")
        if self.shot_overhead_s.lo < 0:
            raise ValueError("shot_overhead_s must be non-negative")

    def T_2033(self) -> tuple[float, float]:
        """The 2033 box temperature band (MeV): T = T_over_Tc_2033 x T_c^lat over the T_c^lat range, (lo, hi) = (225, 300).
        The cheap corner of every band sits at T hi (shortest 1/T), the expensive one at T lo (author ruling 2026-10-05)."""
        return (self.T_over_Tc_2033.lo * self.Tc_lat_2033_MeV.lo, self.T_over_Tc_2033.lo * self.Tc_lat_2033_MeV.hi)


# --------------------------------------------------------------------------- #
# Lattice counting, Eq. (Nq) at app03:68
# --------------------------------------------------------------------------- #

def n_links(L: int, dim: int) -> int:
    return dim * L ** dim


def n_plaquettes(L: int, dim: int) -> int:
    """Periodic lattice: dim(dim-1)/2 plaquettes per site."""
    return dim * (dim - 1) // 2 * L ** dim


def gauge_qubits(L: int, dim: int, qubits_per_link: int) -> int:
    return n_links(L, dim) * qubits_per_link


def fermion_qubits(L: int, dim: int, n_c: int, n_stag: int) -> int:
    return n_c * n_stag * L ** dim


def lq_range(L: int, dim: int, qubits_per_link: int, n_c: int, n_stag: int,
             anc: Tagged) -> tuple[float, float]:
    base = gauge_qubits(L, dim, qubits_per_link) + fermion_qubits(L, dim, n_c, n_stag)
    return (base + anc.lo, base + anc.hi)


# --------------------------------------------------------------------------- #
# Per-link-per-step pricing from the papers' tables (groups.py), ruling R2; rotations at E20
# --------------------------------------------------------------------------- #

def electric_lower_bound_per_link(group: str, ham: str, d: int, eps: float, papers: bool = False) -> float:
    """nF x U_FFT (the COMPILED transform where one exists, else the paper's stated
    lower bound for one) + U_phi where the paper lists it. For 2T this equals
    groups.electric_per_link; for Sigma(72x3), whose U_FFT is CONJECTURE (no FFT
    constructed, arxiv_2511_17437 :905), groups.electric_per_link falls back to the
    dense U_F and this function prices the ruled floor instead. Rotations at the full fit (E20); papers=True
    gives the papers' slope-only price (legacy record)."""
    g = GROUPS[group]
    price = (lambda pc: pc.t_papers(eps)) if papers else (lambda pc: pc.t(eps))
    t = PRIMCOST[ham]["U_F"](d) * price(g.primitives["U_FFT"])
    if "U_phi" in g.primitives:
        t += phi_mult(group, ham) * price(g.primitives["U_phi"])     # s648: Sigma(216x3) H_I calls U_Ph twice
    return t


def phi_mult(group: str, ham: str) -> int:
    """U_phi calls per link per step: 1 unless the group's paper says otherwise (groups.Group.phi_mult; s648:
    Sigma(216x3) 1 at H_KS, 2 at H_I, Gustafson_S648_inprep tab:primcost)."""
    return GROUPS[group].phi_mult.get(ham, 1)


def ham_for(a: Assumptions, group: str) -> str:
    """The Hamiltonian a group is priced at: the 2033 group at hamiltonian_2033 (s648 ruling: H_I), every other
    group (the 2028 2T benchmark) at hamiltonian (H_KS)."""
    return a.hamiltonian_2033.value if group == a.group_2033.value else a.hamiltonian.value


def hop_group(a: Assumptions, group: str) -> str:
    """The Fermion_Primitives entry a gauge group's hop is priced from. Sigma(216x3) has none (s648): the 2033 hop
    is priced at hop_group_2033 (Sigma(72x3)) counts as a PLACEHOLDER (Assumed, NEEDS_AUTHOR)."""
    if group in FERMION_DRAFT:
        return group
    if group == a.group_2033.value:
        return a.hop_group_2033.value
    raise KeyError(f"no Fermion_Primitives counts for {group!r}")


def hop_link(a: Assumptions, group: str, eps: float, **kw) -> dict:
    """The staggered hop per link per Trotter step (r17): the shared groups.hop_link_cost with this chapter's
    N_stag at the circuit's tolerance `eps` (R-TOL), from the unpublished Fermion_Primitives gate counts, rotations
    at the full fit. Keyword arguments (share, undo, mcx, phasing) pass through for sensitivities."""
    return hop_link_cost(hop_group(a, group), int(a.n_stag.lo), eps=eps, synthesis=a.fermion_synthesis.value,
                         toffoli_convention=a.toffoli_convention.value, **kw)


def rhop_link_legacy(a: Assumptions, group: str, eps: float) -> dict:
    """The retired rule R-HOP (6 U_mul + 2 HWP(9)) in the chapter's old convention, for the comparison only."""
    return rhop_link(hop_group(a, group), int(a.n_stag.lo), int(a.n_c.lo), eps,
                     synthesis="rus-slope", toffoli_convention=a.toffoli_convention.value)   # the retired print


def mass_site(a: Assumptions, eps: float) -> dict:
    """Staggered mass per site per Trotter step (r17, option M-B): one HWP(N_c N_stag) group."""
    return staggered_mass_site(int(a.n_c.lo), int(a.n_stag.lo), eps=eps, synthesis=a.fermion_synthesis.value,
                               toffoli_convention=a.toffoli_convention.value)


def hop_primitives(a: Assumptions, group: str, links: int, steps: float, eps: float,
                   sites: int = 0) -> tuple[Primitive, ...]:
    """Per-shot Primitives of the staggered hop (one line per hop piece, each a per-link-step block) and, with
    sites > 0, the staggered mass (one HWP group per site-step), over `steps` Trotter steps at tolerance `eps`."""
    hl = hop_link(a, group, eps)
    out = []
    for it in hl["items"]:
        out.append(Primitive(f"hop_{it['name']}_{group}", links * steps, it["t"], it["status"], it["src"],
                             f"from unpublished gate counts (Fermion_Primitives draft), per link-step: "
                             f"{it['toffoli']} Toffoli + {it['t_direct']} T + {it['n_rot']} rotations at "
                             f"{hl['t_rot']:.4g} T (full fit, eps={eps:.4g}); {it['note']}"))
    if sites:
        ms = mass_site(a, eps)
        out.append(Primitive(f"mass_hwp_{group}", sites * steps, ms["t"], CircuitStatus.SCALING, ms["src"],
                             f"staggered mass, one HWP({ms['k']}) group per site-step: {ms['n_rot']} rotations + "
                             f"{ms['n_toffoli']} Toffoli (option M-B; mu_B N_B folded into the same angle)"))
    return tuple(out)


def rotations_per_link_step(a: Assumptions, group: str, d: int, hop: bool) -> float:
    """Synthesized (arbitrary-angle) rotations per link per Trotter step (ruling R-TOL's N_rot per
    link-step). Each paper primitive contributes its implied rotation count PrimitiveCost.rot
    (= t_log / 1.15) times its tab:primcost multiplicity; with hop=True the Fermion_Primitives hop adds its
    rotations (groups.fermion_hop_counts) and the staggered mass adds one HWP(N_c N_stag) group per site, i.e.
    its rotations / d per link (sites = links / d). The 7-T Toffolis and any exact T are not rotations.
    Independent of eps, so the tolerance follows from it without iteration."""
    g = GROUPS[group]
    ham = ham_for(a, group)
    m = PRIMCOST[ham]
    n = sum(m[p](d) * g.primitives[p].rot for p in MAGNETIC) + m["U_F"](d) * g.primitives["U_FFT"].rot
    if "U_phi" in g.primitives:
        n += phi_mult(group, ham) * g.primitives["U_phi"].rot
    if hop:
        n += fermion_hop_counts(hop_group(a, group), int(a.n_stag.lo))["n_rot"]
        n += hwp_synth_rotations(int(a.n_c.lo) * int(a.n_stag.lo)) / d
    return n


def link_step_primitives(a: Assumptions, group: str, d: int, links: int, steps: float, eps: float
                         ) -> tuple[tuple[Primitive, ...], dict]:
    """The per-shot Primitives of the Trotter evolution: each paper primitive with its
    tab:primcost multiplicity x links x steps, t_each = the paper's a + (b / 1.15)(1.15 log2(1/eps) + 9.2) (E20).
    Returns (primitives, per-link-step pieces)."""
    g = GROUPS[group]
    ham = ham_for(a, group)
    m = PRIMCOST[ham]
    prims = []
    for p in MAGNETIC:
        pc = g.primitives[p]
        prims.append(Primitive(f"{p}_{group}", m[p](d) * links * steps, pc.t(eps), pc.status, pc.src,
                               f"magnetic term: {m[p](d):g}/link-step (tab:primcost H_{ham}); "
                               f"{pc.t_const:g} + {pc.n_rot:g} x (1.15 log2(1/eps) + 9.2) T at eps={eps:.4g} (R-TOL, E20)"))
    fft = g.primitives["U_FFT"]
    prims.append(Primitive(f"U_FFT_{group}", m["U_F"](d) * links * steps, fft.t(eps), fft.status, fft.src,
                           (f"electric term: {m['U_F'](d):g} transforms/link-step; "
                            + ("COMPILED fast transform" if g.has_fft else
                               "NOT CONSTRUCTED: the paper's stated lower bound, a floor (ruling R2)"))))
    if "U_phi" in g.primitives:
        ph = g.primitives["U_phi"]
        pm = phi_mult(group, ham)
        prims.append(Primitive(f"U_phi_{group}", pm * links * steps, ph.t(eps), ph.status, ph.src,
                               f"electric phase, {pm}/link-step at H_{ham} (paper lists it separately)"))
    mag = magnetic_per_link(group, ham, d, eps)
    elec = electric_lower_bound_per_link(group, ham, d, eps)
    elec_dense = electric_per_link(group, ham, d, eps, fft=False)
    # the same pieces at the papers' slope-only price (legacy record, R2) and in the report's convention (T4, E20,
    # equal to mag/elec at 7 T per Toffoli), so a "full fit over papers" sentence is checked, not asserted
    conv = a.toffoli_convention.value
    mag_report = sum(m[p](d) * g.primitives[p].t_report(eps, conv) for p in MAGNETIC)
    elec_report = m["U_F"](d) * fft.t_report(eps, conv) + (phi_mult(group, ham) * g.primitives["U_phi"].t_report(eps, conv)
                                                          if "U_phi" in g.primitives else 0.0)
    pieces = dict(magnetic_per_link=mag, electric_per_link=elec, per_link=mag + elec,
                  electric_dense_per_link=elec_dense, magnetic_report=mag_report, electric_report=elec_report,
                  magnetic_papers=magnetic_per_link(group, ham, d, eps, papers=True),
                  electric_papers=electric_lower_bound_per_link(group, ham, d, eps, papers=True),
                  fft_compiled=g.has_fft)
    return tuple(prims), pieces


# --------------------------------------------------------------------------- #
# r25 (R4): one circuit per evolution depth, each with its own R-TOL tolerance
# --------------------------------------------------------------------------- #

def evolution_2028(a: Assumptions, steps: int) -> dict:
    """The 2028 pure-gauge evolution of `steps` Trotter steps at its own R-TOL tolerance (one shot)."""
    gk = a.group_2028.value
    L, d = int(a.L_2028.lo), int(a.dim_2028.lo)
    links = n_links(L, d)
    if steps <= 0:
        return dict(steps=0, t=0.0, eps=None, t_rot=None, n_rot=0.0, prims=(), pc=None, links=links)
    n_rot = rotations_per_link_step(a, gk, d, hop=False) * links * steps
    eps = eps_rot_for(n_rot, a.eps_syn.lo)
    prims, pc = link_step_primitives(a, gk, d, links, float(steps), eps)
    return dict(steps=steps, t=links * pc["per_link"] * steps, eps=eps, t_rot=t_per_rotation(eps), n_rot=n_rot,
                prims=prims, pc=pc, links=links)


def evolution_2033(a: Assumptions, steps: int) -> dict:
    """The 2033 evolution (gauge + FP hop + staggered mass) of `steps` steps at its own R-TOL tolerance (one shot)."""
    gk = a.group_2033.value
    L, d = int(a.L_2033.lo), int(a.dim_2033.lo)
    links, sites = n_links(L, d), L ** d
    n_rot = rotations_per_link_step(a, gk, d, hop=True) * links * steps
    eps = eps_rot_for(n_rot, a.eps_syn.lo)
    prims, pc = link_step_primitives(a, gk, d, links, float(steps), eps)
    prims = prims + hop_primitives(a, gk, links, float(steps), eps, sites=sites)
    per_link = pc["per_link"] + hop_link(a, gk, eps)["t"] + mass_site(a, eps)["t"] * sites / links
    return dict(steps=steps, t=links * per_link * steps, eps=eps, t_rot=t_per_rotation(eps), n_rot=n_rot,
                prims=prims, per_link=per_link, links=links)


# T-depth schedules from factory.json (Ch. 5; scratchpad factories/ch05_depth_v2.py). A Toffoli is one T layer (its
# seven pi/8 rotations commute); a synthesized rotation is t_rot layers. Toffoli and rotation counts come from groups.py;
# the numbers below are the schedule (tree depth, packing width, colour classes / hop rounds) of that analysis.
SCHED_2028 = dict(classes_lo=2, classes_hi=4, mul_tree=8, tof_width=3, inv_layers=2, tr_rot_layers=3,
                  fft_tof_width=2, fft_rot_width=5)
SCHED_2033 = dict(classes_lo=5, classes_hi=11, hop_rounds_lo=7, hop_rounds_hi=17, mul_tree=8, tof_width=7,
                  inv_width=4, tr_rot_layers_lo=2, fft_tof_width=3, rot_width=9, n_fields_parallel=3)


def step_depth_2028(a: Assumptions, t_rot: float) -> tuple[float, float]:
    """T-depth of one 2028 Trotter step, (low, high), with `workspace_2028` qubits for the U_FFT ancillas.
    low: plaquettes in 2 classes with tree products (12 -> 8 sequential U_mul) and packed Toffolis/rotations;
    high: 4 classes of serial plaquette chains. Both run ceil(64 / (w / 2)) rounds of U_FFT (2 ancilla each).
    At w = 64 this is factory.json's low 628 / 'mid' 3882 per step. factory.json's printed-LQ case (w = 2, one
    primitive at a time, F* = 1.5) is its own schedule and is not reproduced here."""
    g = GROUPS[a.group_2028.value].primitives
    mul, inv, tr, fft = g["U_mul"], g["U_inv"], g["U_Tr"], g["U_FFT"]
    L, d = int(a.L_2028.lo), int(a.dim_2028.lo)
    n_fft = PRIMCOST[a.hamiltonian.value]["U_F"](d) * n_links(L, d)
    rounds = math.ceil(n_fft / (int(a.workspace_2028.lo) // fft.ancilla))
    S = SCHED_2028
    plaq_lo = S["mul_tree"] * mul.toffoli / S["tof_width"] + S["inv_layers"] * inv.toffoli + S["tr_rot_layers"] * t_rot
    fft_lo = fft.toffoli / S["fft_tof_width"] + fft.n_rot / S["fft_rot_width"] * t_rot
    plaq_hi = 12 * mul.toffoli + 6 * inv.toffoli + tr.n_rot * t_rot
    fft_hi = fft.toffoli + fft.n_rot * t_rot
    return (S["classes_lo"] * plaq_lo + rounds * fft_lo, S["classes_hi"] * plaq_hi + rounds * fft_hi)


def step_depth_2033(a: Assumptions, t_rot: float) -> tuple[float, float]:
    """T-depth of one 2033 Trotter step, (low, high): factory.json's Ch. 5 schedule. low (200 ancilla): 5 plaquette
    classes, 7 hop matchings, packed rotations; high (100 ancilla): 11 classes, 17 hop rounds, serial chains, the 3
    fields of a hop in parallel. 4.68e4 / 6.05e5 per step at t_rot = 28.54 (factory.json)."""
    gk = a.group_2033.value
    g = GROUPS[gk].primitives
    mul, inv, tr, fft, phi = (g[k] for k in ("U_mul", "U_inv", "U_Tr", "U_FFT", "U_phi"))
    # s648: the factory.json schedule was built for H_KS (12 U_mul, 6 U_inv, 1 U_Tr, 2 U_F, 1 U_phi per link at d = 3);
    # at H_I each primitive's layers are scaled by its multiplicity ratio (DERIVED HERE, not a factory.json schedule)
    ham, d3 = ham_for(a, gk), int(a.dim_2033.lo)
    m, mk = PRIMCOST[ham], PRIMCOST["KS"]
    r_mul, r_inv, r_tr = (m[p](d3) / mk[p](d3) for p in ("U_mul", "U_inv", "U_Tr"))
    n_f, n_ph = m["U_F"](d3), phi_mult(gk, ham)
    it = {x["name"]: x for x in hop_link(a, a.group_2033.value, 1e-3)["items"]}   # counts only (eps-free)
    sq, rot, fl, dg, ph = (it[k] for k in ("colour_squish_and_parity", "colour_rotations", "hop_squish_and_flags",
                                          "hop_diagonalizers", "hop_phasing_hwp"))
    nf = SCHED_2033["n_fields_parallel"]
    rot_per_frame = rot["n_rot"] / (4 * nf)                 # 184 colour rotations per frame application per field
    n_diag = dg["t_direct"] // 2                            # 180 controlled diagonalizers (2 T + 4 rotations each)
    ms = staggered_mass_site(int(a.n_c.lo), int(a.n_stag.lo), eps=1e-3)
    S = SCHED_2033
    plaq_lo = (r_mul * S["mul_tree"] * mul.toffoli / S["tof_width"] + r_inv * 2 * inv.toffoli / S["inv_width"]
               + r_tr * (tr.toffoli / S["tof_width"] + S["tr_rot_layers_lo"] * t_rot))
    elec_lo = (n_f * (fft.toffoli / S["fft_tof_width"] + fft.n_rot / S["rot_width"] * t_rot)
               + n_ph * phi.n_rot / S["rot_width"] * t_rot)
    hop_lo = ((sq["toffoli"] + fl["toffoli"]) / nf + (2 * rot_per_frame / nf) * t_rot
              + n_diag / (3 * nf) * (2 + 2 * t_rot) + (5 + t_rot))
    mass_lo = 4 + t_rot
    plaq_hi = r_mul * 12 * mul.toffoli + r_inv * 6 * inv.toffoli + r_tr * (tr.toffoli + tr.n_rot * t_rot)
    elec_hi = n_f * (fft.toffoli + fft.n_rot * t_rot) + n_ph * (phi.toffoli + phi.n_rot * t_rot)
    hop_hi = (sq["toffoli"] + fl["toffoli"] + ph["toffoli"] + 4 * rot_per_frame * t_rot
              + n_diag / nf * (2 + 4 * t_rot) + ph["n_rot"] * t_rot)
    mass_hi = ms["n_toffoli"] + ms["n_rot"] * t_rot
    return (elec_lo + S["classes_lo"] * plaq_lo + S["hop_rounds_lo"] * hop_lo + mass_lo,
            elec_hi + S["classes_hi"] * plaq_hi + S["hop_rounds_hi"] * hop_hi + mass_hi)


def fisher_shots(cbar: float, xs, fs, target: float) -> float:
    """Shots for relative error `target` on Gamma from +/-1 outcomes with mean C(t) = cbar exp(-Gamma t), per-shot
    variance 1 - C^2, a fraction fs[i] of the shots at Gamma t = xs[i]; c0 and Gamma both free (Fisher, 2x2)."""
    fcc = fcg = fgg = 0.0
    for x, f in zip(xs, fs):
        mu = cbar * math.exp(-x)
        v = 1.0 - mu * mu
        dc, dg = math.exp(-x), -cbar * x * math.exp(-x)
        fcc += f * dc * dc / v
        fcg += f * dc * dg / v
        fgg += f * dg * dg / v
    return fcc / (fcc * fgg - fcg * fcg) / target ** 2


def tiers_2028(a: Assumptions) -> dict:
    """R1/R4/R10: the 2028 benchmark's timeslices t = 0, a_t, ..., 6 a_t, each shot evolved to its own time at its own
    R-TOL tolerance; first result = slices 0..first_result_steps_2028, full benchmark = all seven. T-depth per shot
    from step_depth_2028 at the circuit's t_rot with `workspace_2028` qubits of primitive workspace."""
    tg, t0, tr = a.t_gate_s.lo, a.shot_overhead_s.lo, REACTION_TIME_S
    n_full, n_first = int(a.trotter_steps_2028.lo), int(a.first_result_steps_2028.lo)
    ns = a.shots_per_timeslice.lo
    slices = []
    for j in range(n_full + 1):
        e = evolution_2028(a, j)
        dep = tuple(j * x for x in step_depth_2028(a, e["t_rot"])) if j else (0.0, 0.0)
        wall = ns * (max(e["t"] * tg, dep[1] * tr) + t0)          # one machine, 10 factories, high schedule
        slices.append(dict(steps=j, t=e["t"], eps=e["eps"], depth=dep, wall_s=wall, shots=ns))

    def run(js):
        sl = [slices[j] for j in js]
        shots = sum(x["shots"] for x in sl)
        t_avg = sum(x["shots"] * x["t"] for x in sl) / shots
        d_avg = tuple(sum(x["shots"] * x["depth"][k] for x in sl) / shots for k in (0, 1))
        return dict(shots=shots, t_avg=t_avg, d_avg=d_avg, wall_s=sum(x["wall_s"] for x in sl),
                    t_max=max(x["t"] for x in sl))
    return dict(slices=slices, first=run(range(n_first + 1)), full=run(range(n_full + 1)))


def second_spacing_assumptions(a: Assumptions) -> Assumptions:
    """T3: the campaign's second spacing on the same 3^3 lattice at its own step (a_t scaled with a_s at fixed
    anisotropy, so the a_t/10 fit grid and Delta t ||H|| are held)."""
    from dataclasses import replace as _replace
    a2 = a.second_spacing_a_s_fm.lo
    return _replace(a, a_s_fm=Stated(a2, "the campaign's second spacing (second_spacing_a_s_fm)"),
                    a_t_fm=Stated(a.a_t_fm.lo * a2 / a.a_s_fm.lo, "fixed anisotropy at the second spacing"))


def quench_steps(a: Assumptions, dt_fm: float, T_MeV: float, c: float) -> dict:
    """Quench preparation (author ruling 2026-10-05): ramp of `quench_ramp_a_s` lattice time units (a floor) plus a
    thermalization time t_th = c / T, both at the step dt_fm of the same Trotter circuit. Steps rounded up."""
    ramp = math.ceil(round(a.quench_ramp_a_s.lo * a.a_s_fm.lo / dt_fm, 9))
    t_th = c * a.hbarc.lo / T_MeV
    th = math.ceil(round(t_th / dt_fm, 9))
    return dict(c=c, t_th_fm=t_th, ramp=ramp, th=th, prep=ramp + th)


def own_depth_2033(a: Assumptions, _nested: bool = False) -> dict:
    """R1/R4 (r25) with the quench thermal route (2026-10-05): Gamma(k_min) from a two-point fit after the early
    transient (fit_x, fit_frac_early), each shot ONE circuit = quench preparation (ramp + c / T at the a_t/10 grid
    step) + fit evolution to its own fit time (step rule), at its own R-TOL tolerance. Two corners of the band:
    (tau = 1/k, c = quench_c.lo) and (tau = 1/(2 pi T), c = quench_c.hi); the two mixed corners are records.
    First result: one point (T = 1.5 T_c^lat, mu_B = 0) at target_first. Campaign: grid_T x grid_muB points at k_min,
    target_campaign, n_spacings; every point costs the same (the Fock-state sector fixes N_B, no sampler).
    Temperature (author ruling 2026-10-05): T = 1.5 T_c^lat with T_c^lat = 150-200 MeV is a third band coordinate.
    The cheap corner is (T hi = 300 MeV, tau = 1/k, c = 2), the expensive one (T lo = 225 MeV, tau = 1/(2 pi T), c =
    2 pi); the other six corners are records. The in-situ T_c^lat diagnostic (tc_diag_*) is priced per corner as
    prep-only shots (ramp + t_th, no fit evolution) read out in the computational basis.
    Walls: one machine, 1 us per T and 0.1 ms per shot, each shot max(N_T t_gate, D_T t_r) (contract rule 12)."""
    tg, t0, tr = a.t_gate_s.lo, a.shot_overhead_s.lo, REACTION_TIME_S
    hb = a.hbarc.lo
    T_lo, T_hi = a.T_2033()
    L, d = int(a.L_2033.lo), int(a.dim_2033.lo)
    kmin = 2 * math.pi * hb / (L * a.a_s_fm.lo)
    dt = a.dt_over_at_2033.lo * a.a_t_fm.lo                   # 0.01 fm/c: the fit-time grid and the prep step

    def lam_sd_at(T):
        # step rule (2026-10-04): the state-dependent estimate at the shot's temperature, Lambda_sd in a_s^-3
        return trotter_state_dependent((L,) * d, int(a.trotter_n_adj.lo), T * a.a_s_fm.lo / hb, 1.0, 1.0,
                                       a.trotter_eps.lo, ham=a.trotter_dispersion_2033.value)["lambda"]
    lam_sd = {T: lam_sd_at(T) for T in (T_lo, T_hi)}

    def n_rule(n_grid, T):
        if a.trotter_rule_2033.value != "state":
            return n_grid
        t_lat = n_grid * dt / a.a_s_fm.lo
        return max(n_grid, math.ceil(round(math.sqrt(lam_sd[T] * t_lat ** 3 / a.trotter_eps.lo), 9)))
    # T1 (2026-10-05): fast end = 1/k_min (free-streaming dephasing); the hydrodynamic T/((eta/s) k^2) is a record

    def tau_of(kind, T):
        return a.tau_fast_1k.lo * hb / kmin if kind == "k" else a.tau_slow_2piT.lo * hb / (2 * math.pi * T)
    taus = (tau_of("k", T_hi), tau_of("2piT", T_lo))          # the box band: 1/k, 1/(2 pi T) at T lo
    cs = (a.quench_c.lo, a.quench_c.hi)
    fs = (a.fit_frac_early.lo, 1.0 - a.fit_frac_early.lo)
    n_pts = a.grid_T.lo * a.grid_muB.lo * a.campaign_k.lo
    n_mu0 = a.mu_b_zero_fraction.lo * n_pts
    n_mup = n_pts - n_mu0
    yr = a.seconds_per_year.lo
    n_diag = a.tc_diag_points.lo * a.tc_diag_shots_per_point.lo

    def corner(T, kind, c):
        tau = tau_of(kind, T)
        q = quench_steps(a, dt, T, c)
        grid = [math.ceil(round(x * tau / dt, 9)) for x in (a.fit_x.lo, a.fit_x.hi)]
        xr = [n * dt / tau for n in grid]
        n_first = fisher_shots(a.cbar.lo, xr, fs, a.target_first.lo)
        n_camp = fisher_shots(a.cbar.lo, xr, fs, a.target_campaign.lo)
        steps = [n_rule(n, T) for n in grid]               # step rule on the fit evolution (binds on the late shots)
        pts = []
        for n, f in zip(steps, fs):
            e = evolution_2033(a, q["prep"] + n)           # ONE circuit: prep + fit, at its own R-TOL eps
            t_prep = e["per_link"] * e["links"] * q["prep"]
            dep = tuple((q["prep"] + n) * x for x in step_depth_2033(a, e["t_rot"]))
            pts.append(dict(steps=n, prep_steps=q["prep"], frac=f, t=e["t"], t_prep=t_prep, t_fit=e["t"] - t_prep,
                            eps=e["eps"], t_rot=e["t_rot"], depth=dep, evo=e))

        def wall(n_shots, base=False):
            return sum(n_shots * p["frac"] * ((p["t"] * tg if base else max(p["t"] * tg, p["depth"][1] * tr)) + t0)
                       for p in pts)

        first_s = wall(n_first)
        pt_s = wall(n_camp)
        per_sp = n_pts * pt_s
        rows = [(n_camp * p["frac"] * n_pts, p) for p in pts]
        shots_sp = sum(w_ for w_, _ in rows)
        t_avg = sum(w_ * p["t"] for w_, p in rows) / shots_sp
        d_avg = tuple(sum(w_ * p["depth"][k] for w_, p in rows) / shots_sp for k in (0, 1))
        last = pts[-1]
        # in-situ T_c^lat diagnostic: prep-only circuit (ramp + t_th at this corner's c and T) at its own R-TOL eps
        e_diag = evolution_2033(a, q["prep"])
        d_diag = tuple(q["prep"] * x for x in step_depth_2033(a, e_diag["t_rot"]))
        diag_s = n_diag * (max(e_diag["t"] * tg, d_diag[1] * tr) + t0)
        return dict(
            T=T, tau=tau, tau_kind=kind, c=c, quench=q, steps=steps, grid_steps=grid, x=xr, n_first=n_first,
            n_camp=n_camp, pts=pts, lam_sd=lam_sd[T],
            t_shot_max=last["t"], t_shot_min=pts[0]["t"], t_prep=last["t_prep"], t_fit_max=last["t_fit"],
            t_fit_min=pts[0]["t_fit"], prep_over_fit=last["t_prep"] / last["t_fit"],
            first_shots=n_first, first_s=first_s, first_s_baseline=wall(n_first, base=True),
            pt_s=pt_s, per_spacing_s=per_sp, campaign_s=a.n_spacings.lo * per_sp,
            campaign_s_baseline=a.n_spacings.lo * n_pts * wall(n_camp, base=True),
            shots_per_spacing=shots_sp, shots_campaign=a.n_spacings.lo * shots_sp,
            t_avg=t_avg, d_avg=d_avg, f_star_min_shot_hi=min(p["t"] / p["depth"][1] for p in pts),
            eps_l=a.faults_per_shot.lo / last["t"], campaign_yr=a.n_spacings.lo * per_sp / yr,
            diag_shots=n_diag, diag_t_shot=e_diag["t"], diag_s=diag_s, diag_over_first=diag_s / first_s,
        )
    out = [corner(T_hi, "k", cs[0]), corner(T_lo, "2piT", cs[1])]
    # records: the other six corners of the (T, tau, c) band
    mixed = tuple(corner(T, kind, c) for T in (T_lo, T_hi) for kind in ("k", "2piT") for c in cs
                  if (T, kind, c) not in ((T_hi, "k", cs[0]), (T_lo, "2piT", cs[1])))
    for r in out:                                               # deepest shot of the campaign (first spacing so far)
        r.update(deep_pts=r["pts"], t_shot_deep=r["t_shot_max"], t_prep_deep=r["t_prep"],
                 prep_steps_deep=r["quench"]["prep"], deep_spacing=a.a_s_fm.lo)
    second = None
    if a.second_spacing_own_step.value and not _nested:
        # T3 (2026-10-05): the second spacing at its own step; the campaign is the sum of the two spacings
        if a.n_spacings.lo != 2:
            raise ValueError("second_spacing_own_step prices exactly two spacings")
        second = own_depth_2033(second_spacing_assumptions(a), _nested=True)
        for r, r2 in zip(out, second["r"]):
            r["campaign_s"] = r["per_spacing_s"] + r2["per_spacing_s"]
            r["campaign_s_baseline"] = (r["campaign_s_baseline"] + r2["campaign_s_baseline"]) / 2
            s1, s2 = r["shots_per_spacing"], r2["shots_per_spacing"]
            r["t_avg"] = (s1 * r["t_avg"] + s2 * r2["t_avg"]) / (s1 + s2)          # exports over both spacings
            r["d_avg"] = tuple((s1 * x + s2 * y) / (s1 + s2) for x, y in zip(r["d_avg"], r2["d_avg"]))
            r["f_star_min_shot_hi"] = min(r["f_star_min_shot_hi"], r2["f_star_min_shot_hi"])
            r["shots_campaign"] = r["shots_per_spacing"] + r2["shots_per_spacing"]
            r["campaign_yr"] = r["campaign_s"] / yr
            if r2["t_shot_max"] > r["t_shot_max"]:
                r.update(deep_pts=r2["pts"], t_shot_deep=r2["t_shot_max"], t_prep_deep=r2["t_prep"],
                         prep_steps_deep=r2["quench"]["prep"], deep_spacing=second_spacing_assumptions(a).a_s_fm.lo)
            r["eps_l_deep"] = a.faults_per_shot.lo / r["t_shot_deep"]
    else:
        for r in out:
            r["eps_l_deep"] = r["eps_l"]
    return dict(taus=taus, cs=cs, Ts=(T_lo, T_hi), dt=dt, kmin=kmin, n_pts=n_pts, n_mu0=n_mu0, n_mup=n_mup,
                lam_sd=lam_sd, r=tuple(out), mixed=mixed, second=second)


# --------------------------------------------------------------------------- #
# Model
# --------------------------------------------------------------------------- #

def _model_2028(a: Assumptions) -> Result:
    gk = a.group_2028.value
    g = GROUPS[gk]
    L, d = int(a.L_2028.lo), int(a.dim_2028.lo)
    links = n_links(L, d)
    plaq = n_plaquettes(L, d)
    gq = gauge_qubits(L, d, g.link_qubits)
    lq = lq_range(L, d, g.link_qubits, 0, 0, a.anc_2028)

    steps = float(a.trotter_steps_2028.lo)
    dt = a.dt_over_at_2028.lo * a.a_t_fm.lo
    # R-TOL: this circuit's rotation count sets its tolerance
    n_rot = rotations_per_link_step(a, gk, d, hop=False) * links * steps
    eps = eps_rot_for(n_rot, a.eps_syn.lo)
    breakdown, pc = link_step_primitives(a, gk, d, links, steps, eps)
    t_link_step = pc["per_link"]
    t_per_step = links * t_link_step
    t_shot = t_per_step * steps
    t_max = steps * dt
    steps_in_rfi = a.rfi_2028_t.lo / t_per_step

    Tc = a.tc_over_sqrt_sigma.lo * a.sqrt_sigma_MeV.lo
    T = a.T_over_Tc_2028.lo * Tc
    L_perp = L * a.a_s_fm.lo
    k = 2 * math.pi * a.hbarc.lo / L_perp

    gate_s = a.t_gate_s.lo
    t0 = a.shot_overhead_s.lo
    shot_time_s = t_shot * gate_s + t0             # deepest shot on one machine: 1 us per T + ~0.1 ms per shot (E27)
    l2 = math.log2(1.0 / eps)
    # r25 (R1, R4, R10): timeslices t = 0..6 a_t, each shot at its own depth; first result = t = 0 and one step;
    # `workspace_2028` qubits of primitive workspace so the 10-factory baseline is fed (F* >= 10)
    tiers = tiers_2028(a)
    first, full = tiers["first"], tiers["full"]
    c0, dl = a.cbar.lo, a.first_result_falloff_2028.lo
    n30 = 2 * (1 - c0 * c0) / (c0 * c0 * a.target_first.lo ** 2 * dl * dl)
    s1 = tiers["slices"][1]
    w1 = s1["wall_s"] / s1["shots"]                  # one-step shot, charged as in tiers_2028
    w30 = n30 / 2 * t0 + n30 / 2 * w1
    w30_one = n30 * w1
    w = int(a.workspace_2028.lo)
    lq = (lq[0] + w, lq[1] + w)
    shots = full["shots"]
    wall_s = full["wall_s"]
    exp = depth_exports(full["t_avg"], full["d_avg"], shots, a.t_gate_s, a.shot_overhead_s)
    from dataclasses import replace as _replace
    d_w2 = step_depth_2028(_replace(a, workspace_2028=Assumed(2, "factory.json printed-LQ case")), t_per_rotation(eps))

    inter = {
        "n_links_2028": links,                          # 32
        "n_plaq_2028": plaq,                            # 16
        "qubits_per_link_2028": g.link_qubits,          # 5
        "gauge_qubits_2028": gq,                        # 160
        "lq_2028_lo": lq[0], "lq_2028_hi": lq[1],       # 177, 178 (160 + 1-2 Hadamard + 16 workspace, R10)
        "n_rot_per_link_step_2028": n_rot / (links * steps),     # 89.5 = 6 x 0 + 3 x 0 + 11/2 + 2 x 42
        "n_rot_per_shot_2028": n_rot,                   # 17,184 synthesized rotations in one shot (R-TOL)
        "eps_rot_2028": eps,                            # 7.628e-4 = sqrt(1e-2 / 17184)
        "synthesis_error_per_shot_2028": n_rot * eps ** 2,       # 1e-2 (R-TOL budget, randomized synthesis)
        "t_per_rotation_papers": 1.15 * l2,             # 11.91 at 7.628e-4 (papers' slope-only, legacy record)
        "t_per_rotation_report": t_per_rotation(eps),   # 21.11 (full fit, E20): "21.1 T"
        "magnetic_t_per_link_step_2028": pc["magnetic_per_link"],       # 1.124e3 -> "1.1e3"
        "electric_t_per_link_step_2028": pc["electric_per_link"],       # 1.969e3 -> "2.0e3"
        "electric_over_magnetic_2028": pc["electric_per_link"] / pc["magnetic_per_link"],   # 1.75
        "t_per_link_step_2028": t_link_step,            # 3.093e3 -> "3.1e3"
        "fft_compiled_2028": pc["fft_compiled"],        # True (arxiv_2408_00075)
        "magnetic_report_over_papers_2028": pc["magnetic_report"] / pc["magnetic_papers"],   # full fit over the papers' slope-only
        "electric_report_over_papers_2028": pc["electric_report"] / pc["electric_papers"],   # full fit over the papers' slope-only
        "t_per_step_2028": t_per_step,                  # 9.90e4
        "trotter_steps_2028": steps,                    # 6
        "t_per_shot_2028": t_shot,                      # 5.939e5 -> box "5.9e5"
        "t_per_shot_2028_over_rfi": t_shot / a.rfi_2028_t.lo,   # 5.94: "5.9x the 1e5-T envelope"
        "per_step_cut_needed_2028": t_shot / a.rfi_2028_t.lo,   # 5.94: ">= 5.9x" (ruling (b), at E20)
        "levers_close_2028_lo": a.per_step_cut_levers.lo >= t_shot / a.rfi_2028_t.lo,   # False: 3x does not
        "levers_close_2028_hi": a.per_step_cut_levers.hi >= t_shot / a.rfi_2028_t.lo,   # True: 10x does
        "steps_within_rfi_2028": steps_in_rfi,          # 1.01: "1e5 T buys a single step" (app03:74)
        "t_max_within_rfi_fm": math.floor(steps_in_rfi) * dt,                        # 0.025 fm/c (rule step a_t/4)
        "t_max_within_rfi_at_dt_over_10_fm": math.floor(steps_in_rfi) * a.dt_over_at_2033.lo * a.a_t_fm.lo,  # 0.01
        "dt_2028_fm": dt,                               # 0.025 (step rule; was 0.1)
        "t_max_2028_fm": t_max,                         # 0.15 (step rule; was 0.6)
        "Tc_2028_MeV": Tc,                              # 492.8
        "T_2028_MeV": T,                                # 739.2 -> "~740"
        "L_perp_fm": L_perp,                            # 0.8
        "k_2028_GeV": k / 1e3,                          # 1.55
        "k_over_T_2028": k / T,                         # 2.1
        "xi_anisotropy": a.a_s_fm.lo / a.a_t_fm.lo,     # 2
        "shots_2028": shots,                            # 7e3 = 1e3 x 7 timeslices (t = 0..6 a_t)
        "timeslices_2028": len(tiers["slices"]),        # 7 ("six steps give seven times"; the r23 box said 10)
        "gate_time_from_factories_s": REACTION_TIME_S / 10,        # 1e-6 = t_gate_s: 10 us cycle / ~10 factories
        "shot_time_2028_s": shot_time_s,                # 0.594 -> "~0.6 s" (deepest shot, requirements row)
        "wall_time_2028_s": wall_s,                     # 2059.0 s: each slice at its own depth
        "wall_time_2028_min": wall_s / 60,              # 34.3 -> box "34 min"
        "shot_overhead_s": t0,                          # 1e-4
        # r25 tiers (R1, R4) and workspace (R10)
        "t_per_slice_2028": tuple(x["t"] for x in tiers["slices"]),          # 0, 9.473e4, ..., 5.939e5
        "eps_per_slice_2028": tuple(x["eps"] for x in tiers["slices"]),
        "t_first_result_shot_2028": first["t_max"],     # 9.473e4 -> "9.5e4", inside 1e5
        "t_first_result_shot_2028_over_rfi": first["t_max"] / a.rfi_2028_t.lo,   # 0.947
        "shots_first_result_2028": first["shots"],      # 2e3 (t = 0 and one step)
        "wall_first_result_2028_s": first["wall_s"],    # 94.93 s -> "1.6 min"
        "wall_first_result_2028_min": first["wall_s"] / 60,
        # verifier r25 item 1: the first result at R1's 30%, priced in the prose only (delta open, box keeps 2e3)
        "shots_first_result_2028_at_30pct": n30,        # 8.458e4 -> "~8.5e4"
        "wall_first_result_2028_at_30pct_h": w30 / 3600,   # 1.115 h -> "~1.1 h" (half the shots at t = 0)
        "wall_first_result_2028_at_30pct_all_one_step_h": w30_one / 3600,   # 2.23 h
        "workspace_2028": w,                            # 16
        "lq_2028_without_workspace": (lq[0] - w, lq[1] - w),   # 161-162 (the r23 box)
        "f_star_2028_per_step": (t_per_step / step_depth_2028(a, t_per_rotation(eps))[1],
                                 t_per_step / step_depth_2028(a, t_per_rotation(eps))[0]),   # 10.66 / 57.1 at 6 steps' eps
        "f_star_2028_at_2_workspace": t_per_step / d_w2[1],   # ~1.7: one U_FFT at a time, the 10 factories idle
        "wall_time_2028_min_r23": a.shots_per_timeslice.lo * 10 * shot_time_s / 60,   # 99.0: the r23 box (1e4 x 6 steps)
        **{k: exp[k] for k in ("t_per_shot", "t_depth_per_shot", "f_star", "floor_wall_s", "factories_for_1yr")},
        "wall_serial_s": exp["wall_serial_s"], "baseline_ok": exp["baseline_ok"],
        "wall_first_result_s": first["wall_s"],
        "wall_campaign_s": full["wall_s"],
        # the retired working figure, for the record (groups.py 2T.plaquette_t, Stated 1e3)
        "legacy_t_per_shot_2028_at_1e3_per_plaq": plaq * g.plaquette_t.lo * steps,   # 9.6e4 (old box "1e5")
    }
    inter.update(trotter_check_2028(a))
    inter["dt_first_2028_fm"] = a.dt_first_2028_fm.lo
    # Ruling (C) 2026-10-05: the 2028 box stays a non-thermal placeholder benchmark; the quench route's cost here is
    # a prose sentence only (one circuit: ramp + c / T thermalization + the one-step first result, at its own eps)
    qs = [quench_steps(a, dt, T, cc) for cc in (a.quench_c.lo, a.quench_c.hi)]
    tq = [evolution_2028(a, q["prep"] + int(a.first_result_steps_2028.lo))["t"] for q in qs]
    inter.update({
        "quench_inv_T_steps_2028": a.hbarc.lo / T / dt,                      # 10.7 steps per 1/T
        "quench_ramp_steps_2028": qs[0]["ramp"],                              # 8
        "quench_th_steps_2028": tuple(q["th"] for q in qs),                   # 22, 68
        "quench_prep_steps_2028": tuple(q["prep"] for q in qs),               # 30, 76
        "quench_first_result_shot_2028": tuple(tq),                           # ~3.2e6, ~8.1e6 T (prep + one step)
        "quench_first_result_over_rfi_2028": tuple(x / a.rfi_2028_t.lo for x in tq),   # ~32, ~81
    })
    notes = (
        "Per-link-step cost from the papers' tables at E20 (7 T/Toffoli, every rotation at the full fit 1.15 log2(1/eps) "
        "+ 9.2 with n_rot = log coefficient / 1.15; eps = 7.628e-4 by R-TOL from this circuit's 17,184 rotations, 21.11 T "
        "each): magnetic 1.124e3 (6 U_mul + 3 U_inv + U_Tr/2 at d=2, arxiv_2208_12309) + electric 1.969e3 (2 U_FFT, "
        "arxiv_2408_00075) = 3.093e3 T; the chapter prints 1.1e3 + 2.0e3 = 3.1e3.",
        "Box prints 5.9e5 T = 32 links x 3.1e3 x 6; the exact product is 5.939e5, 5.9x the 1e5-T RFI envelope. "
        "The chapter keeps its 6 steps, now of a_t/4 to 0.15 fm/c under the step rule (2026-10-04; was a_t to 0.6 fm/c), "
        "and says the instance needs a >=5.9x per-step cut, which the named levers' 3-10x reaches only at its upper part "
        "(ruling (b) 2026-10-01).",
        "Box prints 178 LQ = 32 x 5 + 1-2 + 16 workspace; the range is 177-178.",
        "r25 (R1, R4): timeslices t = 0, a_t, ..., 6 a_t, 1e3 shots each, every shot evolved to its own time at its own "
        "R-TOL tolerance. First result (t = 0 and one step): 2e3 shots, deepest 9.473e4 T (inside 1e5, no per-step cut), "
        "94.9 s. Full benchmark: 7e3 shots, 2059 s = 34.3 min; the six-step shot is the 5.939e5 above. Serial on one "
        "machine (E27); the r23 box's 1e4 x 6-step pricing was 99.0 min.",
        "R10: 16 workspace qubits (8 U_FFT at once) give F* = 10.7-10.9 at the high schedule (57-60 low), so the 10-factory "
        "baseline is fed; LQ 160 + 1-2 + 16 = 177-178 (box '178'). Without them F* ~ 1.7 (factory.json: 1.5) and ~11 h.",
        "Legacy record: the papers' slope-only price (pre-E20) gave magnetic 1.07e3, electric 1.20e3, shot 4.358e5 (4.4x); "
        "E20 adds 5% to the magnetic term and 65% to the electric term.",
    )
    return Result("2028", lq, (first["t_max"], full["t_max"]), breakdown, inter, shots=(shots, shots),
                  wall_time_s=(wall_s, wall_s), notes=notes)


def _common_2033(a: Assumptions) -> dict:
    gk = a.group_2033.value
    g = GROUPS[gk]
    L, d = int(a.L_2033.lo), int(a.dim_2033.lo)
    nc, ns = int(a.n_c.lo), int(a.n_stag.lo)
    links = n_links(L, d)
    plaq = n_plaquettes(L, d)
    gq = gauge_qubits(L, d, g.link_qubits)
    fq = fermion_qubits(L, d, nc, ns)
    lq = lq_range(L, d, g.link_qubits, nc, ns, a.anc_2033)

    dt = a.dt_over_at_2033.lo * a.a_t_fm.lo
    steps = a.t_max_2033_fm.lo / dt
    sites = L ** d
    # R-TOL: the evolution circuit's rotation count (gauge + FP hop + mass) sets its tolerance
    n_rot_link_step = rotations_per_link_step(a, gk, d, hop=True)
    n_rot = n_rot_link_step * links * steps
    eps = eps_rot_for(n_rot, a.eps_syn.lo)
    evol_prims, pc = link_step_primitives(a, gk, d, links, steps, eps)
    # the staggered hop from the unpublished Fermion_Primitives gate counts (r17), app03:85, and the staggered
    # mass (one HWP group per site); 2028 is pure gauge and has neither
    hl = hop_link(a, gk, eps)
    ms = mass_site(a, eps)
    evol_prims = evol_prims + hop_primitives(a, gk, links, steps, eps, sites=sites)
    mass_per_link = ms["t"] * sites / links
    pc = dict(pc, gauge_per_link=pc["per_link"], hop_per_link=hl["t"], mass_per_link=mass_per_link,
              per_link=pc["per_link"] + hl["t"] + mass_per_link)
    t_per_step = links * pc["per_link"]
    t_evol = t_per_step * steps
    # sensitivities of the hop reading at the same eps (not the box)
    hop_share_draft = hop_link(a, gk, eps, share="draft")["t"]     # r17 headline, E21 sensitivity
    t_evol_share_draft = links * (pc["gauge_per_link"] + hop_share_draft + mass_per_link) * steps
    hop_no_undo = hop_link(a, gk, eps, undo=False)["t"]
    rhop_old = rhop_link_legacy(a, gk, eps)["t"]
    # the conjectured qubitized BE (app03 route (a) / budget split) is a different circuit: 10x fewer operations, read here as
    # 10x fewer rotations (N_rot / be_trim; flagged for the author), so its own R-TOL tolerance is looser
    n_rot_be = n_rot / a.be_trim.lo
    eps_be = eps_rot_for(n_rot_be, a.eps_syn.lo)
    per_link_be = (magnetic_per_link(gk, ham_for(a, gk), d, eps_be)
                   + electric_lower_bound_per_link(gk, ham_for(a, gk), d, eps_be) + hop_link(a, gk, eps_be)["t"]
                   + mass_site(a, eps_be)["t"] * sites / links)
    t_be = links * per_link_be * steps / a.be_trim.lo
    t_evol_dense = links * (pc["magnetic_per_link"] + pc["electric_dense_per_link"] + pc["hop_per_link"]
                            + pc["mass_per_link"]) * steps

    Lfm = L * a.a_s_fm.lo
    kmin = 2 * math.pi * a.hbarc.lo / Lfm
    Ts = a.T_2033()                                                      # (225, 300) MeV: every T-dependent record is a pair
    LT = tuple(Lfm * T / a.hbarc.lo for T in Ts)                         # 0.68, 0.91: the box is smaller than 1/T
    tau_at_eta_s = tuple(T * a.hbarc.lo / (a.eta_s_assumed.lo * kmin ** 2) for T in Ts)   # 0.069, 0.092 fm/c (record)
    n_decay_at_eta_s = tuple(a.t_max_2033_fm.lo / x for x in tau_at_eta_s)                # 58, 43 (record)

    n_grid = a.grid_T.lo * a.grid_muB.lo * a.grid_k.lo * a.grid_t.lo
    shots = a.shots_per_gridpt.lo * n_grid                       # r23 record (96 grid points x 1e3)
    n_grid_p = a.n_grid_prose.lo
    gate_s = a.t_gate_s.lo
    t0 = a.shot_overhead_s.lo                     # ~0.1 ms per shot (report rule)
    # the papers' fiducial synthesis error (app03:78 sensitivity sentence), same groups.py functions
    # linear split of eps_T over this circuit's rotations (r18 verifier fix; the papers' 7.1e12 Sigma(72x3) H_I
    # benchmark implies log2(1/eps) = 63.4, i.e. eps_T = 1e-8 split over ~1e11 rotations, not 1e-8 per rotation)
    eps_fid = a.papers_fiducial_eps_T.lo / n_rot
    mag_fid = magnetic_per_link(gk, ham_for(a, gk), d, eps_fid)
    elec_fid = electric_lower_bound_per_link(gk, ham_for(a, gk), d, eps_fid)
    hop_fid = hop_link(a, gk, eps_fid)["t"] + mass_site(a, eps_fid)["t"] * sites / links
    t_evol_fid = links * (mag_fid + elec_fid + hop_fid) * steps
    n_calls = a.shots_per_gridpt.lo * n_grid_p * a.n_spacings.lo

    # register extensions and volume requirements (app03:77,79,162)
    Lx = int(a.L_ext.lo)
    lq_12 = lq_range(Lx, d, g.link_qubits, nc, ns, a.anc_2033)
    Lx4 = int(a.L_ext_small.lo)
    lq_4 = lq_range(Lx4, d, g.link_qubits, nc, ns, a.anc_2033)
    m_pi_mid = 0.5 * (a.m_pi_MeV.lo + a.m_pi_MeV.hi)
    L_mpiL4 = a.m_pi_L_target.lo * a.hbarc.lo / m_pi_mid
    kmin_at_L_ext_fm = 2 * math.pi * a.hbarc.lo / a.L_ext_fm.lo
    kmin_at_12cubed = 2 * math.pi * a.hbarc.lo / (Lx * a.a_s_fm.lo)
    L_kmin_eq_T = tuple(2 * math.pi * a.hbarc.lo / T for T in Ts)      # 5.5, 4.1 fm
    q360_dense = math.ceil(math.log2(a.sigma360_order.lo))
    q360 = int(a.sigma360_qubits_compiled.lo)

    return dict(g=g, L=L, d=d, links=links, plaq=plaq, gq=gq, fq=fq, lq=lq, dt=dt, steps=steps,
                n_rot_link_step=n_rot_link_step, n_rot=n_rot, eps=eps, n_rot_be=n_rot_be, eps_be=eps_be,
                per_link_be=per_link_be,
                evol_prims=evol_prims, pc=pc, t_per_step=t_per_step, t_evol=t_evol, t_be=t_be,
                t_evol_dense=t_evol_dense, Lfm=Lfm, kmin=kmin, Ts=Ts, LT=LT,
                n_grid=n_grid, shots=shots,
                tau_at_eta_s=tau_at_eta_s, n_decay_at_eta_s=n_decay_at_eta_s, t0=t0,
                mag_fid=mag_fid, elec_fid=elec_fid, hop_fid=hop_fid, t_evol_fid=t_evol_fid, hl=hl, ms=ms,
                sites=sites, hop_share_draft=hop_share_draft,
                t_evol_share_draft=t_evol_share_draft, hop_no_undo=hop_no_undo, rhop_old=rhop_old,
                n_calls=n_calls, gate_s=gate_s, lq_12=lq_12, lq_4=lq_4, m_pi_mid=m_pi_mid, L_mpiL4=L_mpiL4,
                kmin_at_L_ext_fm=kmin_at_L_ext_fm, kmin_at_12cubed=kmin_at_12cubed,
                L_kmin_eq_T=L_kmin_eq_T, q360=q360, q360_dense=q360_dense)


def _intermediates_2033(a: Assumptions, c: dict) -> dict:
    g = c["g"]
    pc = c["pc"]
    fb = a.faults_per_shot.lo
    return {
        "ceil_log2_order_2033": math.ceil(math.log2(g.order)),          # 8 (dense; not used)
        "qubits_per_link_2033": g.link_qubits,                           # 9 (compiled register, R1)
        "qubits_per_link_2033_chapter_before_R1": int(g.link_qubits_chapter.lo),   # 8, the retired value
        "n_links_2033": c["links"],                                      # 81
        "n_plaq_2033": c["plaq"],                                        # 81
        "gauge_qubits_2033": c["gq"],                                    # 729
        "fermion_qubits_2033": c["fq"],                                  # 243
        "lq_2033_lo": c["lq"][0], "lq_2033_hi": c["lq"][1],              # 1072, 1172 -> box "~1120 (1072-1172)"
        "lq_2033_mid": 0.5 * (c["lq"][0] + c["lq"][1]),                  # 1122 -> "~1120"
        "lq_2033_lo_over_1000_limit": c["lq"][0] / 1000.0,               # 1.072: "slightly overfill" (app03:77)
        "log10_hilbert_dim_2033": c["links"] * math.log10(g.order) + c["fq"] * math.log10(2),  # 216^81 2^243
        "dt_2033_fm": c["dt"],                                           # 0.01
        "trotter_steps_2033": c["steps"],                                # 400
        "n_rot_per_link_step_2033": c["n_rot_link_step"],                # 4093.33 = 1159 gauge + 2933 FP hop + 4/3 mass
        "n_rot_per_shot_2033": c["n_rot"],                               # 1.567e8 synthesized rotations (R-TOL, quench in the circuit)
        "eps_rot_2033": c["eps"],                                        # 7.99e-6 = sqrt(1e-2 / 1.5658e8), chapter "8.0e-6"
        "synthesis_error_per_shot_2033": c["n_rot"] * c["eps"] ** 2,     # 1e-2
        "t_per_rotation_papers": 1.15 * math.log2(1.0 / c["eps"]),       # 19.34 (papers' slope-only, legacy record)
        "t_per_rotation_full_fit": t_per_rotation(c["eps"], a.fermion_synthesis.value),   # 28.54 (every rotation, E20)
        "n_rot_per_shot_be": c["n_rot_be"],                              # 1.32624e7 (N / 10, conjectured BE)
        "eps_rot_be": c["eps_be"],                                       # 2.746e-5
        "t_per_rotation_papers_be": 1.15 * math.log2(1.0 / c["eps_be"]), # 17.43
        "t_per_link_step_be_before_trim": c["per_link_be"],              # per-link-step at eps_be, before /10
        "magnetic_t_per_link_step_2033": pc["magnetic_per_link"],        # 1.888e4 -> "1.9e4"
        "electric_t_per_link_step_2033_lower_bound": pc["electric_per_link"],   # 3.394e4 -> "3.4e4" (floor)
        "electric_t_per_link_step_2033_dense": pc["electric_dense_per_link"],   # 9.25e6 (compiled U_F)
        "electric_dense_over_magnetic_2033": pc["electric_dense_per_link"] / pc["magnetic_per_link"],   # 490 -> "~490x"
        "electric_over_magnetic_2033": pc["electric_per_link"] / pc["magnetic_per_link"],   # 1.80
        "fft_compiled_2033": pc["fft_compiled"],                         # False: no Sigma(72x3) FFT
        "gauge_t_per_link_step_2033": pc["gauge_per_link"],              # 5.281e4 (magnetic + electric floor)
        "hop_t_per_link_step_2033": pc["hop_per_link"],                  # 1.0976e5 (FP draft counts, squish per link, full fit)
        "hop_toffoli_per_link_step_2033": c["hl"]["toffoli"],            # 3,672 (E21: 14,314 at share='draft')
        "hop_t_direct_per_link_step_2033": c["hl"]["t_direct"],          # 360
        "hop_n_rot_per_link_step_2033": c["hl"]["n_rot"],                # 2,933
        "hop_t_toffoli_per_link_step_2033": c["hl"]["t_toffoli"],        # 25,704
        "hop_t_rotations_per_link_step_2033": c["hl"]["t_rotations"],    # 83,694
        "hop_items_t_per_link_step_2033": {it["name"]: it["t"] for it in c["hl"]["items"]},
        "hop_undo_t_per_link_step_2033": pc["hop_per_link"] - c["hop_no_undo"],   # 3.150e4, the ESTIMATED V_g^dag (1,104 rotations)
        "hop_t_per_link_step_2033_share_draft": c["hop_share_draft"],    # 1.8425e5 (r17 headline, per-application squish)
        "t_per_shot_2033_mu0_evolution_share_draft": c["t_evol_share_draft"],   # 7.683e9 (the r18 box value)
        "hop_t_per_link_step_2033_no_undo": c["hop_no_undo"],            # 7.826e4 (the draft's 2 C^G, squish per link)
        "rhop_t_per_link_step_2033_retired": c["rhop_old"],              # 6972.7 (R-HOP, slope-only, same eps)
        "hop_over_rhop_2033": pc["hop_per_link"] / c["rhop_old"],        # 15.7
        "hop_phasing_hwp_ancilla_2033": hwp_ancilla(2 * int(a.n_c.lo) * int(a.n_stag.lo)),   # 16 = 18 - w(18) (E26, r22; was 21), inside 100-200
        "hop_t_per_step_2033": c["links"] * pc["hop_per_link"],          # 8.890e6
        "hop_t_per_shot_2033": c["links"] * pc["hop_per_link"] * c["steps"],   # 3.556e9
        "hop_over_gauge_2033": pc["hop_per_link"] / pc["gauge_per_link"],       # 2.08
        "hop_share_of_step_2033": pc["hop_per_link"] / pc["per_link"],          # 0.675 -> "67%"
        "mass_t_per_site_step_2033": c["ms"]["t"],                       # 163.1 = 4 x 28.54 + 7 x 7
        "mass_hwp_k_2033": c["ms"]["k"],                                 # 9
        "mass_hwp_ancilla_2033": c["ms"]["ancilla"],                     # 7 = 9 - w(9) (E26, r22; was 11)
        "mass_t_per_step_2033": c["sites"] * c["ms"]["t"],               # 4.405e3
        "mass_share_of_step_2033": c["sites"] * c["ms"]["t"] / c["t_per_step"],   # 3.3e-4
        "t_per_link_step_2033": pc["per_link"],                          # 1.6262e5 -> "1.6e5" (gauge + hop + mass/link)
        "magnetic_report_over_papers_2033": pc["magnetic_report"] / pc["magnetic_papers"],   # full fit over the papers' slope-only
        "electric_report_over_papers_2033": pc["electric_report"] / pc["electric_papers"],   # full fit over the papers' slope-only
        "t_per_step_2033": c["t_per_step"],                              # 1.317e7 -> "1.3e7"
        "t_per_shot_2033_mu0_evolution": c["t_evol"],                    # 5.269e9 -> box "5.3e9"
        "t_per_shot_2033_mu0_evolution_dense_fourier": c["t_evol_dense"],   # 3.040e11 -> "~3.0e11"
        "t_per_shot_2033_be": c["t_be"],                                 # 1.444e9 at the BE's own eps -> "~1.4e9"
        "L_2033_fm": c["Lfm"],                                           # 0.6
        "k_min_2033_GeV": c["kmin"] / 1e3,                               # 2.07 -> "~2.1"
        # temperature band (author ruling 2026-10-05): T = 1.5 T_c^lat, T_c^lat = 150-200 MeV; pairs are (T lo, T hi)
        "T_c_lat_2033_MeV": (a.Tc_lat_2033_MeV.lo, a.Tc_lat_2033_MeV.hi),   # (150, 200), ASSUMED range
        "T_2033_MeV": c["Ts"],                                           # (225, 300)
        "T_a_s_2033": tuple(T * a.a_s_fm.lo / a.hbarc.lo for T in c["Ts"]),   # 0.228, 0.304 (was 0.162)
        "LT_2033": c["LT"],                                              # 0.68, 0.91 -> "LT ~ 0.7-0.9", below 1 (was 0.49)
        "k_min_over_T_2033": tuple(c["kmin"] / T for T in c["Ts"]),      # 9.2, 6.9 -> "k_min ~ 7-9 T" (was 12.9)
        "eta_s_assumed": a.eta_s_assumed.lo,                             # 0.15 (app03:79, record)
        "tau_decay_at_eta_s_assumed_fm": c["tau_at_eta_s"],              # 0.069, 0.092 -> "~0.07-0.09 fm/c" (record)
        "n_decay_times_at_eta_s_assumed": c["n_decay_at_eta_s"],         # 58, 43 (record)
        "m_pi_L_2033_lo": a.m_pi_MeV.lo * c["Lfm"] / a.hbarc.lo,         # 0.91
        "m_pi_L_2033_hi": a.m_pi_MeV.hi * c["Lfm"] / a.hbarc.lo,         # 1.22  -> "~1"
        "L_for_mpiL4_fm": c["L_mpiL4"],                                  # 2.26 -> "2.3"
        "L_ext_sites": (c["L_mpiL4"] / a.a_s_fm.lo),                     # 11.3 -> "~12^3"
        "lq_12cubed_lo": c["lq_12"][0], "lq_12cubed_hi": c["lq_12"][1],  # 62,308 -> "~6e4"
        "factor_beyond_2033": c["lq_12"][0] / 1000.0,                    # 62 -> "~60" (vs the 1000-LQ limit)
        "factor_beyond_2033_register": c["lq_12"][0] / (0.5 * (c["lq"][0] + c["lq"][1])),   # 56 -> "~60"
        "lq_12cubed_over_lanl": c["lq_12"][0] / a.lanl_lq.lo,            # 1.7 -> "comparable"
        "lq_4cubed_lo": c["lq_4"][0], "lq_4cubed_hi": c["lq_4"][1],      # 2404-2504 -> "low thousands"
        "k_min_over_T_at_L_ext_fm": tuple(c["kmin_at_L_ext_fm"] / T for T in c["Ts"]),   # 2.4, 1.8 -> "1.8-2.4 T at L ~ 2.3 fm"
        "k_min_over_T_at_12cubed": tuple(c["kmin_at_12cubed"] / T for T in c["Ts"]),     # 2.3, 1.7 -> "1.7-2.3 T on 12^3"
        "L_for_kmin_eq_T_fm": c["L_kmin_eq_T"],                          # 5.5, 4.1 -> "L >~ 4-6 fm" (was 7.75)
        "sigma360_qubits_per_link_dense": c["q360_dense"],               # 11
        "sigma360_qubits_per_link": c["q360"],                           # 12 (compiled)
        "sigma360_extra_lq": c["links"] * (c["q360"] - g.link_qubits),   # 243
        # R3: RFI 1e9 ceiling, eps_l printed as a requirement at 0.1 expected faults per shot
        "ceiling_2033_t": a.ceiling_2033_t.lo,                           # 1e9
        "ratio_reference_to_ceiling": c["t_evol"] / a.ceiling_2033_t.lo,   # 14.74: the 4 fm/c reference shot (record)
        "eps_l_required_reference": fb / c["t_evol"],                    # 6.784e-12 (reference shot, record)
        "eps_l_required_reference_at_1_fault": 1.0 / c["t_evol"],        # 6.784e-11 (the pre-R3 convention, record)
        "n_grid_2033": c["n_grid"],                                      # 96 (r23 record)
        "shots_2033": c["shots"],                                        # 9.6e4 (r23 record)
        "mu_b_zero_fraction_min_realizable": 1.0 / a.grid_muB.lo,        # 0.25
        # D3: the papers' own fiducial eps_T = 1e-8 (linear split) instead of the R-TOL per-circuit tolerance
        "magnetic_t_per_link_step_2033_at_papers_fiducial_eps": c["mag_fid"],           # 1.917e4
        "electric_t_per_link_step_2033_at_papers_fiducial_eps": c["elec_fid"],          # 8.262e4 (floor)
        "hop_t_per_link_step_2033_at_papers_fiducial_eps": c["hop_fid"],                # 2.338e5 (FP hop + mass/link at eps=1e-8/N_rot, full fit)
        "t_per_link_step_2033_at_papers_fiducial_eps": c["mag_fid"] + c["elec_fid"] + c["hop_fid"],   # 3.356e5 -> "3.4e5" (app03:78)
        "t_per_shot_2033_mu0_evolution_at_papers_fiducial_eps": c["t_evol_fid"],        # 1.087e10 -> "1.1e10" (app03:78)
        "campaign_horizon_stated_yr": a.campaign_horizon_yr.lo,          # 5
        "shot_time_2033_reference_s": c["t_evol"] * c["gate_s"] + c["t0"],   # 1.474e4 s (reference shot, record)
        "shot_overhead_s": c["t0"],                                      # 1e-4 (report rule)
        "n_subroutine_calls": c["n_calls"],                              # 1.92e5 -> "~2e5" (app03:198)
        # utility box (app03:189): not a resource number, but a stated product
        "utility_musd_per_yr": a.heavy_ion_research_musd_per_yr.lo * a.utility_fraction.lo,   # 2.3727 -> "$2.4M/yr"
        "utility_musd_campaign": (a.heavy_ion_research_musd_per_yr.lo * a.utility_fraction.lo
                                  * a.campaign_horizon_yr.lo),           # 11.86 -> "$12M over a 5-year campaign"
        # the retired working figure, for the record (groups.py S72x3.plaquette_t, Stated 5e3)
        "legacy_t_per_shot_2033_at_5e3_per_plaq": (c["plaq"] * GROUPS["S72x3"].plaquette_t.lo * c["steps"]),   # 1.62e8 (old box, Sigma(72x3) working figure)
    }


_NOTES_2033 = (
    "s648 (author ruling, H. Lamm 2026-10-05): Sigma(216x3) (|G| = 648) with the improved Hamiltonian H_I at a_s = 0.2 "
    "fm, priced from the in-preparation draft (Gustafson_S648_inprep, groups.py 'S216x3', qubit encoding). Register: "
    "81 links x 11 qubits = 891 gauge + 243 fermion + 100-200 ancilla = 1234-1334; box '~1280 (1234-1334)', 1.23-1.33x "
    "the 1000-LQ reference (under 1.5x: optimization caveat). Was 9 qubits/link, 1072-1172, with Sigma(72x3).",
    "Per-link-step gauge cost at H_I (tab:primcost H_I row; U_Ph twice) at the reference shot's R-TOL eps = 7.989e-6 "
    "(1.567e8 rotations per 400-step shot, 4836.33 per link-step = 1902 gauge + 2933 hop + 4/3 mass; 28.67 T per "
    "rotation): multiplication / inversion / trace 2.858e5 + 4 U_FFT + 2 U_Ph 5.894e4 = 3.447e5, '3.4e5' (= the "
    "draft's closed form plus 9.2 T per implied rotation, groups.c_t_full). The draft's U_FFT contains the Sigma(72x3) "
    "FFT at the authors' re-synthesis (532 + 497.95 log), above the published lower-bound estimate; no dense / fast "
    "split remains (electric_dense == electric). H_KS for the same group would be 1.036e5 (3.3x less).",
    "Staggered hop: no Fermion_Primitives entry for Sigma(216x3); priced at the Sigma(72x3) counts (hop_group_2033, "
    "PLACEHOLDER): 3,672 Toffoli + 360 T + 2,933 rotations = 1.1016e5 T per link-step (undo 3.166e4); share='draft' "
    "1.847e5, reference shot 1.715e10. Mass one HWP(9) per site, 163.7 T. Link-step 3.447e5 + 1.1016e5 + 54.6 = "
    "4.550e5, '4.5e5'; hop 24%.",
    "Reference shot (mu_B=0, 400 steps of 0.01 fm/c, every shot to 4 fm/c; records) = 81 x 4.550e5 x 400 = 1.474e10, "
    "'1.5e10'; the step rule would need 1914 steps. BE (N_rot / 10 at its own eps) 1.444e9, '1.4e9'.",
    "Trotter step rule at H_I (trotter_dispersion_2033 = 'I'): Lambda_sd = 45.79 a_s^-3 on the tree-level improved "
    "dispersion (32.76 at H_KS, 1.40x; T a_s = 0.23-0.30 moves it by < 1e-2). Deepest fit shots (19, 32) steps on the "
    "(18, 26)-step a_t/10 grid, estimate (0.092, 0.098), free field (0.015, 0.019); worst case at g^2 = 1 with the H_I "
    "norms (e x 4/3, plaquette part of b x 23/12) 64-68. Second spacing (0.15 fm) late slow shot 48 steps (80 at the "
    "retired 160 MeV).",
    "Freezing (ruling (2), literature): Wilson 3+1D beta_f = 3.18(3) (Sigma(72x3)) and 3.80(5) (Sigma(216x3)), g^2_f = "
    "1.89, 1.58; non-freezing H_I trajectory at 0.15-0.2 fm ASSUMED by analogy with S(1080) (arxiv_1906_11213; there a "
    "different modification, one extra single-plaquette term, not O(a^2) improvement). The clause 'the Hamiltonian "
    "(anisotropic) limit moves freezing to stronger coupling' is an in-house statement, uncited (ruling 2026-10-05).",
    "THERMAL STATE BY QUENCH (author ruling 2026-10-05; replaces E-rho-OQ at mu_B = 0 and the KMS Gibbs sampler at "
    "mu_B > 0, both retired from the model): electric vacuum with the filled staggered Dirac sea, or an N_B-sector Fock "
    "state, ramp of one a_s (20 steps, a floor) + t_th = c / T at the a_t/10 grid step, c = 2 (132 steps at 300 MeV) to "
    "2 pi (552 steps at 225 MeV), then the fit evolution, one circuit at its own R-TOL eps. ETH with dynamical fermions "
    "at these couplings on 3^3 and the relaxation of the momentum density within t_th are ASSUMED (checked by exact "
    "evolution on Z2 up to 5x3, 2T up to 2x2, Z2 + fermion N = 16 only). No shot penalty, no new primitives, one "
    "thermal route for every grid point.",
    "TEMPERATURE (author ruling 2026-10-05): T = 1.5 T_c^lat, T_c^lat the crossover of the simulated theory (unknown; "
    "150-200 MeV ASSUMED, bracketed by HotQCD's chiral-limit 132 MeV and physical 156.5 MeV; heavier pions raise it, "
    "extra tastes lower it), measured in situ from prep-only diagnostic shots (10 energies x 100 shots, 0.18-0.67 yr, "
    "~6% of the first result). T = 225-300 MeV: 1/T = 65.8-87.7 steps, LT = 0.68-0.91 (< 1: a small hot box, not a "
    "medium), k_min = 6.9-9.2 T, tau slow end 1/(2 pi T) = 0.140 fm/c at 225 MeV. The cheap corner of every band "
    "sits at 300 MeV, the expensive one at 225 MeV. Was T = 160 MeV (physical QCD): 123.3 steps, LT = 0.49, k_min = 13 T.",
)


_NOTES_2033_R25 = (
    "r25 (R4): each shot is evolved only to its own fit time, at its own R-TOL tolerance. The fit takes two points after "
    "the early transient, Gamma t = 0.5 and 1.8 (shot audit design B), 20% / 80% of the shots, both rounded up onto the "
    "Delta t = 0.01 fm/c grid. Decay band (T1, 2026-10-05): tau = 1/k = 0.0955 fm/c to 1/(2 pi T) = 0.140 fm/c (225 MeV; "
    "0.196 at the retired 160 MeV); fit steps (5, 18) and (7, 26) on the grid; under the step rule "
    "(2026-10-04) the slow-end late shot runs 32 steps over its 0.26 fm/c, and at the second spacing (own step, T3) 48. "
    "The 4 fm/c (400-step) pricing above is the reference shot; its campaign keys are records.",
    "R1 first result: Gamma(k_min) at T = 1.5 T_c^lat, mu_B = 0, 30% statistical: 1.60e4 / 1.53e4 shots (Fisher, Cbar = "
    "0.16, c0 and Gamma free). Campaign: 20% at 16 (T, mu_B) points x k_min x 2 spacings, 3.60e4 / 3.44e4 (first "
    "spacing) and 3.60e4 / 3.47e4 (second) shots per point. Shots unchanged by the quench ruling; the walls carry the "
    "quench prep in every shot (see the quench note and the wall intermediates).",
    "Walls charge each shot max(N_T x 1 us, D_T x 10 us) + 0.1 ms at the high depth schedule (factory.json); every "
    "shot feeds the factories (F* = 34.7-35.0), so no shot is charged its depth.",
)


def _own_breakdown(a: Assumptions, o: dict) -> tuple[Primitive, ...]:
    """One shot at a_s = 0.2 fm at the deepest fit point of the slow corner (T = 225 MeV, tau = 1/(2 pi T), c = 2 pi):
    the quench preparation (572 steps) and the fit evolution (32 steps) are the same Trotter primitives, so the
    breakdown is the evolution primitives at 604 steps (the box's '2.2e10'). The second spacing's late shot is
    t_deepest_campaign."""
    return o["r"][1]["pts"][-1]["evo"]["prims"]


def _own_intermediates_2033(a: Assumptions, o: dict) -> dict:
    """r25 numbers with the quench prep (2026-10-05): every key here prices each shot at its own depth (R4) as one
    circuit, prep + fit; the un-prefixed reference keys of _intermediates_2033 (4 fm/c shot) are records."""
    lo, hi = o["r"]
    yr = a.seconds_per_year.lo
    hz = a.campaign_horizon_yr.lo * yr
    pair = lambda k: (lo[k], hi[k])
    exp = depth_exports(lo["t_avg"], lo["d_avg"], lo["shots_campaign"], a.t_gate_s, a.shot_overhead_s)
    q = (lo["quench"], hi["quench"])
    T2 = o["second"]["r"] if o.get("second") else None
    return {
        # corner order throughout: (cheap: T hi, tau = 1/k, c = 2), (expensive: T lo, tau = 1/(2 pi T), c = 2 pi)
        "T_corners_MeV": pair("T"),                                    # 300, 225
        "tau_band_fm": o["taus"],                                      # 0.0955 (1/k), 0.1396 (1/(2 pi T) at 225 MeV)
        "quench_c_band": o["cs"],                                      # (2, 2 pi)
        "quench_t_th_fm": tuple(x["t_th_fm"] for x in q),              # 1.316, 5.510 fm/c (c / T at 300 / 225 MeV)
        "quench_inv_T_steps": tuple(a.hbarc.lo / T / o["dt"] for T in pair("T")),   # 65.8, 87.7 steps per 1/T (was 123.3)
        "quench_ramp_steps": q[0]["ramp"],                             # 20 (one a_s at dt = 0.01 fm/c)
        "quench_th_steps": tuple(x["th"] for x in q),                  # 132, 551 (was 247, 775)
        "quench_prep_steps": tuple(x["prep"] for x in q),              # 152, 571 (was 267, 795)
        "quench_prep_t": pair("t_prep"),                               # T at the deepest shot's eps
        "quench_prep_over_fit": pair("prep_over_fit"),                 # preparation-bound shots
        "fit_steps": (tuple(lo["steps"]), tuple(hi["steps"])),         # step rule on the grid counts
        "fit_steps_grid": (tuple(lo["grid_steps"]), tuple(hi["grid_steps"])),   # t / (a_t/10)
        "lambda_sd_corners": pair("lam_sd"),                           # Lambda_sd at T hi / T lo (a_s^-3)
        # in-situ T_c^lat diagnostic (author ruling 2026-10-05): tc_diag_points x tc_diag_shots_per_point prep-only shots
        "tc_diag_shots": lo["diag_shots"],                             # 1e3
        "tc_diag_t_per_shot": pair("diag_t_shot"),                     # prep-only circuit at its own eps
        "tc_diag_wall_yr": (lo["diag_s"] / yr, hi["diag_s"] / yr),
        "tc_diag_over_first_result": pair("diag_over_first"),          # fraction of the first result's wall
        "fit_t_fm": (tuple(n * o["dt"] for n in lo["grid_steps"]), tuple(n * o["dt"] for n in hi["grid_steps"])),
        "fit_x_realized": (tuple(lo["x"]), tuple(hi["x"])),            # (0.524, 1.885), (0.509, 1.834)
        "eps_rot_fit_points": (tuple(p["eps"] for p in lo["pts"]), tuple(p["eps"] for p in hi["pts"])),
        "t_rot_fit_points": (tuple(p["t_rot"] for p in lo["pts"]), tuple(p["t_rot"] for p in hi["pts"])),   # 28.4-29.1
        "t_shot_fit_points": (tuple(p["t"] for p in lo["pts"]), tuple(p["t"] for p in hi["pts"])),
        "t_fit_fit_points": (tuple(p["t_fit"] for p in lo["pts"]), tuple(p["t_fit"] for p in hi["pts"])),
        "f_star_fit_points_hi": (tuple(p["t"] / p["depth"][1] for p in lo["pts"]),
                                 tuple(p["t"] / p["depth"][1] for p in hi["pts"])),   # ~35 (high schedule)
        "t_deepest": pair("t_shot_max"),                               # 6.255e9, 2.234e10: box "6.3e9-2.2e10" (was 1.05e10, 3.14e10)
        "t_fit_deepest": pair("t_fit_max"),                            # 6.9e8, 1.2e9: the fit evolution alone
        "t_earliest": pair("t_shot_min"),                              # the early fit point's shot
        "ratio_deepest_to_budget": tuple(x / a.ceiling_2033_t.lo for x in pair("t_shot_max")),   # 6.3, 22.3
        "eps_l_deepest": pair("eps_l"),                                # 1.6e-11, 4.5e-12
        "shot_time_deepest_s": tuple(x * a.t_gate_s.lo + a.shot_overhead_s.lo for x in pair("t_shot_max")),   # 6.3e3-2.2e4 s
        # T3 (2026-10-05): the campaign's deepest shot over both spacings (the second spacing's late shot)
        "t_deepest_campaign": pair("t_shot_deep"),                     # 7.9e9, 3.0e10 (second spacing, 755 + 48 steps)
        "quench_prep_steps_deepest_campaign": pair("prep_steps_deep"), # 196, 755
        "eps_l_deepest_campaign": pair("eps_l_deep"),                  # 1.3e-11, 3.4e-12
        "deepest_shot_spacing_fm": pair("deep_spacing"),               # 0.15, 0.15
        "second_spacing_quench_prep_steps": tuple(r["quench"]["prep"] for r in T2) if T2 else None,   # 196, 755
        "second_spacing_quench_prep_t": tuple(r["t_prep"] for r in T2) if T2 else None,           # 7.2e9, 2.8e10
        "second_spacing_t_fit_deepest": tuple(r["t_fit_max"] for r in T2) if T2 else None,        # 7.0e8, 1.8e9
        "first_result_corners": {(round(r["T"], 1), r["tau_kind"], round(r["c"], 4)): r["first_s"] / yr
                                 for r in o["r"] + o["mixed"]},        # 8 corners of (T, tau, c); the box quotes the extremes
        "t_deepest_corners": {(round(r["T"], 1), r["tau_kind"], round(r["c"], 4)): r["t_shot_max"]
                              for r in o["r"] + o["mixed"]},
        "shots_first_result": pair("n_first"),                         # 1.600e4, 1.550e4 (30%)
        "shots_per_point_campaign": pair("n_camp"),                    # 3.60e4, 3.49e4 (20%)
        "design_A_over_B": (fisher_shots(a.cbar.lo, (0.0, 1.25), (0.22, 0.78), 1.0)
                            / fisher_shots(a.cbar.lo, (a.fit_x.lo, a.fit_x.hi), (a.fit_frac_early.lo, 1 - a.fit_frac_early.lo), 1.0)),   # 0.366
        "design_C_over_B": (fisher_shots(a.cbar.lo, (1.0, 2.3), (0.2, 0.8), 1.0)
                            / fisher_shots(a.cbar.lo, (a.fit_x.lo, a.fit_x.hi), (a.fit_frac_early.lo, 1 - a.fit_frac_early.lo), 1.0)),   # 2.72
        "shots_per_point_campaign_design_nominal": fisher_shots(a.cbar.lo, (a.fit_x.lo, a.fit_x.hi),
                                                                (a.fit_frac_early.lo, 1 - a.fit_frac_early.lo),
                                                                a.target_campaign.lo),   # 3.42e4 (audit design B)
        "shots_campaign": pair("shots_campaign"),                      # 1.152e6, 1.145e6 (16 points x 2 spacings)
        "n_points_campaign": o["n_pts"], "n_points_mu0": o["n_mu0"], "n_points_mupos": o["n_mup"],   # 16, 4, 12
        "wall_first_result_d": (lo["first_s"] / 86400, hi["first_s"] / 86400),   # 1139, 3924 d (20/80 shot weights)
        "wall_first_result_yr": (lo["first_s"] / yr, hi["first_s"] / yr),        # 3.12, 10.74 yr: box "3.1-11 yr" (was 5.3-15)
        "wall_first_result_d_baseline": (lo["first_s_baseline"] / 86400, hi["first_s_baseline"] / 86400),
        "first_result_over_horizon": (lo["first_s"] / hz, hi["first_s"] / hz),   # 0.62, 2.1: the cheap corner fits 5 yr
        "campaign_per_spacing_yr_own": (lo["per_spacing_s"] / yr, hi["per_spacing_s"] / yr),   # 112, 387 (first spacing)
        "campaign_two_spacings_yr_own": (lo["campaign_s"] / yr, hi["campaign_s"] / yr),        # 254, 938 (was 435, 1325)
        "campaign_two_spacings_yr_own_baseline": (lo["campaign_s_baseline"] / yr, hi["campaign_s_baseline"] / yr),
        "campaign_reduction_needed_own": (lo["campaign_s"] / hz, hi["campaign_s"] / hz),       # 51, 188
        "grid_point_yr_own": (lo["pt_s"] / yr, hi["pt_s"] / yr),      # 7.0, 24.2 yr per 20% grid point
        "grid_points_in_horizon_own": (math.floor(hz / lo["pt_s"]), math.floor(hz / hi["pt_s"])),   # 0, 0
        "f_star_min_shot_hi": pair("f_star_min_shot_hi"),              # 34.7-35.0: no shot is charged its depth
        # R9 exports: the campaign at the (1/k, c = 2) corner (shot-averaged T and depth, low/high schedule);
        # the two walls are bands over the corners
        **{k: exp[k] for k in ("t_per_shot", "t_depth_per_shot", "f_star", "floor_wall_s", "factories_for_1yr")},
        "wall_serial_s": exp["wall_serial_s"], "baseline_ok": exp["baseline_ok"],
        "wall_first_result_s": (lo["first_s"], hi["first_s"]),
        "wall_campaign_s": (lo["campaign_s"], hi["campaign_s"]),
    }


def _referee_intermediates_2033(a: Assumptions, o: dict) -> dict:
    """Referee items T1-T3 (2026-10-04). Derived from the same functions as the box; none is a box number. The G8
    (E-rho-OQ) block is retired with the method (author ruling 2026-10-05; computed values in NEEDS_AUTHOR)."""
    from dataclasses import replace as _replace
    hb, Ts = a.hbarc.lo, a.T_2033()
    gk, d = a.group_2033.value, int(a.dim_2033.lo)
    L = int(a.L_2033.lo)
    links, sites = n_links(L, d), L ** d
    out = {}
    # T1: the primary output is C(t) at the two fit times; its relative error there for the printed shots
    # (prepared-state variance 1 - mu^2 per +/-1 outcome, mu = Cbar e^{-x})
    rel = []
    for r in o["r"]:
        row = []
        for N in (r["n_first"], r["n_camp"]):
            row.append(tuple(math.sqrt((1 - (a.cbar.lo * math.exp(-x)) ** 2) / (f * N)) / (a.cbar.lo * math.exp(-x))
                             for x, f in zip(r["x"], (a.fit_frac_early.lo, 1 - a.fit_frac_early.lo))))
        rel.append(tuple(row))
    out["corr_rel_err_first_at_fit_times"] = (rel[0][0], rel[1][0])     # (0.186, 0.314), (0.186, 0.351)
    out["corr_rel_err_campaign_at_fit_times"] = (rel[0][1], rel[1][1])  # (0.124, 0.210), (0.124, 0.234)
    # T1: non-hydrodynamic decay scales at k_min (free-streaming dephasing 1/k; thermal 1/(2 pi T))
    out["tau_free_streaming_fm"] = hb / o["kmin"]                         # 0.0955 -> "1/k ~ 0.1 fm/c"
    # T1 verifier (2026-10-05): free-streaming massless quanta give the transverse momentum-density correlator
    # 3 (sin x - x cos x) / x^3, x = k t; it reaches 1/e at x = 2.926, i.e. 0.28 fm/c, beyond the slow end 0.140 (225 MeV).
    f = lambda x: 3 * (math.sin(x) - x * math.cos(x)) / x ** 3 - math.exp(-1)
    lo_x, hi_x = 1.0, 4.0
    for _ in range(60):
        mid = 0.5 * (lo_x + hi_x)
        lo_x, hi_x = (mid, hi_x) if f(mid) > 0 else (lo_x, mid)
    out["free_streaming_1e_kt"] = 0.5 * (lo_x + hi_x)                       # 2.926 -> "kt ~ 2.9"
    out["free_streaming_1e_fm"] = out["free_streaming_1e_kt"] * hb / o["kmin"]   # 0.279 -> "0.28 fm/c"
    # T1 (2026-10-05): a third (test) time at Gamma t ~ 1.2 with a fifth more shots, both spacings, all 16 points.
    # sigma of ln C(t3) minus its two-point interpolation; added campaign wall (each test shot carries the quench prep
    # of its corner plus the evolution to t3 at the grid step count); NOT PRICED in the box.
    sig, frac = [], []
    o2 = o["second"] if o.get("second") else None
    for k, r in enumerate(o["r"]):
        x1, x2 = r["x"]
        n3 = math.ceil(round(1.2 * r["tau"] / o["dt"], 9))
        x3 = n3 * o["dt"] / r["tau"]
        def vln(x, n):
            mu = a.cbar.lo * math.exp(-x)
            return (1 - mu * mu) / (n * mu * mu)
        N = r["n_camp"]
        w2 = (x3 - x1) / (x2 - x1)
        sig.append(math.sqrt(vln(x3, 0.2 * N) + (1 - w2) ** 2 * vln(x1, 0.2 * N) + w2 ** 2 * vln(x2, 0.8 * N)))
        added = 0.0
        for rr, aa, oo in ((r, a, o), ) + (((o2["r"][k], second_spacing_assumptions(a), o2),) if o2 else ()):
            m3 = math.ceil(round(1.2 * rr["tau"] / oo["dt"], 9))
            t3 = evolution_2033(aa, rr["quench"]["prep"] + m3)["t"]
            added += 0.2 * rr["n_camp"] * o["n_pts"] * t3 * a.t_gate_s.lo
        frac.append(added / r["campaign_s"])
    out["test_time_sigma_lnC"] = tuple(sig)                                # 0.295, 0.289 -> "+-0.3 in ln C"
    out["test_time_campaign_fraction"] = tuple(frac)                       # ~0.20, ~0.20 -> "a fifth"
    # T3: the campaign's second spacing on the same 3^3 lattice shrinks the box; holding L = 0.6 fm needs 4^3
    a2 = a.second_spacing_a_s_fm.lo
    out["second_spacing_a_s_fm"] = a2
    out["second_spacing_L_fm_on_3cubed"] = L * a2                          # 0.45
    out["second_spacing_kmin_over_T_on_3cubed"] = tuple(2 * math.pi * hb / (L * a2) / T for T in Ts)   # 12.2, 9.2 (was 17.2)
    out["second_spacing_mpiL_on_3cubed"] = (a.m_pi_MeV.lo * L * a2 / hb, a.m_pi_MeV.hi * L * a2 / hb)   # 0.68, 0.91
    a4 = _replace(a, L_2033=Stated(4, "fixed-volume second spacing"), a_s_fm=Stated(a2, "fixed-volume second spacing"))
    o4 = own_depth_2033(a4, _nested=True)
    out["fixed_volume_4cubed_lq"] = lq_range(4, d, GROUPS[gk].link_qubits, int(a.n_c.lo), int(a.n_stag.lo), a.anc_2033)
    out["fixed_volume_4cubed_t_deepest"] = (o4["r"][0]["t_shot_max"], o4["r"][1]["t_shot_max"])   # with the quench prep
    # T3 (verifier): the second spacing at its own step (Delta t = a_t/10 with a_t scaled at fixed anisotropy, i.e.
    # Delta t ||H|| held): the slow-decay fit times need 4/3 the steps, and so does the quench (c / T at a finer step).
    a3 = _replace(a, a_s_fm=Stated(a2, "second spacing at its own step"),
                  a_t_fm=Stated(a.a_t_fm.lo * a2 / a.a_s_fm.lo, "fixed anisotropy at the second spacing"),
                  )
    o3 = own_depth_2033(a3, _nested=True)
    yr = a.seconds_per_year.lo
    out["second_spacing_own_step_fit_steps"] = (tuple(o3["r"][0]["steps"]), tuple(o3["r"][1]["steps"]))   # (5,19), (14,80)
    out["second_spacing_own_step_t_deepest"] = (o3["r"][0]["t_shot_max"], o3["r"][1]["t_shot_max"])     # 1.4e10, 4.2e10
    out["campaign_yr_second_spacing_own_step"] = tuple((r["per_spacing_s"] + r3["per_spacing_s"]) / yr
                                                       for r, r3 in zip(o["r"], o3["r"]))      # = the box since T3
    # record: the pre-2026-10-05 convention, the second spacing at the first spacing's per-shot cost
    out["campaign_yr_second_spacing_same_cost"] = tuple(2 * r["per_spacing_s"] / yr for r in o["r"])   # 748, 2223
    return out


# --------------------------------------------------------------------------- #
# Referee G7 (2026-10-04): Trotter error at a chosen step (shared with Ch. 6). DERIVED HERE.
# --------------------------------------------------------------------------- #
# One second-order step S = e^{-iD dt/2} e^{-iE dt} e^{-iD dt/2}, E the electric terms, D everything diagonal in the
# group basis (plaquettes and the staggered hop, which commute with each other). Lattice units of the spatial spacing.

def trotter_worst_case(n_links: int, d: int, g2: float, emax_casimir: float, hop_norm: float, t: float, dt: float,
                       eps: float, ham: str = "KS") -> dict:
    """(i) Worst-case second-order bound (childs2021theory), the same-link estimate of Ch. 9's trotter_bound_steps:
    ||S^N - U|| <= t dt^2 Lambda, Lambda = n_links (e^2 b / 3 + e b^2 / 6), with
        e = (g^2/2) C_max / 2           half-range of one link's electric term
        b = hop_norm + 2(d-1) 6/g^2     the hop on the link (n_stag N_c / 2 for staggered fields) and its 2(d-1)
                                        plaquettes, -(2/g^2) Re Tr U_p = -(1/g^2)(Tr U_p + h.c.), |Re Tr U_p| <= 3
                                        for SU(3)
    The plaquette normalization is the Wilson-matched one (it reduces to (E^2 + B^2)/2 in the continuum, the
    normalization of the free-field w^2 = 4 sum sin^2 used in (ii) and (iii), and Ch. 9's trotter_bound_steps). It
    was 3/g^2 before the ch05ch06 verifier (read off a factor-2 typo in app03 eq:Hks, now corrected; NEEDS_AUTHOR).
    Cross-link terms (hops sharing a site, plaquettes sharing two links) and the mass only add.
    ham="I" (s648, DERIVED HERE): H_I = K_KS - (1/12) sum tr (L_2 - R_1)^2 + (5/3) V_KS - (1/12) V_rect (arxiv_2203_02823,
    beta_K0 = 5/6, beta_K1 = 1/6 rewritten with its tr(R_1 L_2) identity). Per link the electric range grows by at most
    (1/12) x 4 = 1/3 (one two-link pair per link, |L_2 - R_1|^2 <= 4 C_max), e x 4/3; the link sits in 2(d-1) plaquettes
    at 5/3 and 6(d-1) rectangles at 1/12, so the plaquette part of b grows by (10/3 + 1/2) / 2 = 23/12."""
    e = 0.5 * (g2 / 2.0) * emax_casimir * (4.0 / 3.0 if ham == "I" else 1.0)
    b = hop_norm + 2 * (d - 1) * 6.0 / g2 * (23.0 / 12.0 if ham == "I" else 1.0)
    lam = n_links * (e * e * b / 3.0 + e * b * b / 6.0)
    return {"lambda": lam, "err": lam * t * dt * dt, "n_steps_at_eps": math.ceil(math.sqrt(lam * t ** 3 / eps))}


def free_gauge_frequencies(dims, ham: str = "KS") -> list:
    """Transverse lattice-photon frequencies w_k^2 = 4 sum_i sin^2(pi n_i / L_i), k != 0 (lattice units).
    ham="I" (s648): the tree-level Symanzik dispersion w^2 = sum_i (khat_i^2 + khat_i^4 / 12), khat_i^2 = 4 sin^2(pi n_i /
    L_i), the reading Ch. 6 uses (ch06 free_gauge_frequencies_h); assumption trotter_dispersion_2033."""
    import itertools
    ws = []
    for n in itertools.product(*[range(L) for L in dims]):
        kh = [4.0 * math.sin(math.pi * ni / L) ** 2 for ni, L in zip(n, dims)]
        w2 = sum(kh) + (sum(x * x for x in kh) / 12.0 if ham == "I" else 0.0)
        if w2 > 1e-12:
            ws.append(math.sqrt(w2))
    return ws


def trotter_state_dependent(dims, n_adj: int, temp: float, t: float, dt: float, eps: float, ham: str = "KS") -> dict:
    """(ii) State-dependent estimate (Alves, Lamm, Liu, App. 'Second-order Trotter approximation'; evaluation as in
    Ch. 9's trotter_state_dependent_steps, generalized to L_1 x ... x L_d and T = 0): error <= t dt^2 Lambda_sd,
    Lambda_sd = sigma(o1)/12 + sigma(o2)/24 = sigma/8 with the connected variance of o1 = -w^3 P^2, o2 = -w^3 X^2 per
    transverse oscillator, 2 [(w^3/2) coth(w/2T)]^2, summed over n_adj (d - 1) copies per k != 0. 'fast_phase' is the
    second-order frequency shift of the fastest mode over the shot, w_max^3 dt^2 t / 24."""
    d = len(dims)
    mult = n_adj * (d - 1)
    ws = free_gauge_frequencies(dims, ham)

    def coth(w):
        return 1.0 / math.tanh(w / (2.0 * temp)) if temp > 0 else 1.0
    var = mult * sum(2.0 * (0.5 * w ** 3 * coth(w)) ** 2 for w in ws)
    lam = math.sqrt(var) / 8.0
    return {"lambda": lam, "err": lam * t * dt * dt, "n_steps_at_eps": math.ceil(math.sqrt(lam * t ** 3 / eps)),
            "fast_phase": max(ws) ** 3 * dt * dt * t / 24.0}


def free_gauge_mode_overlap(w: float, dt: float, n: int, temp: float) -> float:
    """|Tr rho U^dag S^N| for one oscillator H = p^2/2 + w^2 x^2/2 in a thermal state (exact). In x~ = sqrt(w) x,
    p~ = p / sqrt(w) one step is the symplectic map M (theta = w dt), U^dag is a rotation by -w t, and W = R(-w t) M^N
    factors as a squeeze r and a rotation delta (polar decomposition; tr W^T W = 2 cosh 2r). With <n|S(r)|n> =
    sech(r)^{1/2} P_n(sech r) and the Legendre generating function, |Tr rho W| =
    (1 - q) sech(r)^{1/2} / |1 - 2 z sech r + z^2|^{1/2}, q = e^{-w/T}, z = q e^{-i delta}."""
    import cmath
    th = w * dt
    M = ((1 - th * th / 2, th), (-th + th ** 3 / 4, 1 - th * th / 2))
    P = ((1.0, 0.0), (0.0, 1.0))
    for _ in range(n):
        P = ((P[0][0] * M[0][0] + P[0][1] * M[1][0], P[0][0] * M[0][1] + P[0][1] * M[1][1]),
             (P[1][0] * M[0][0] + P[1][1] * M[1][0], P[1][0] * M[0][1] + P[1][1] * M[1][1]))
    c, s = math.cos(th * n), math.sin(th * n)
    a_, b_ = c * P[0][0] - s * P[1][0], c * P[0][1] - s * P[1][1]
    c_, d_ = s * P[0][0] + c * P[1][0], s * P[0][1] + c * P[1][1]
    r = 0.5 * math.acosh(max((a_ * a_ + b_ * b_ + c_ * c_ + d_ * d_) / 2.0, 1.0))
    delta = math.atan2(b_ - c_, a_ + d_)
    q = math.exp(-w / temp) if temp > 0 else 0.0
    x = 1.0 / math.cosh(r)
    z = q * cmath.exp(-1j * delta)
    return (1 - q) * math.sqrt(x) / abs(cmath.sqrt(1 - 2 * z * x + z * z))


def free_gauge_trotter_error(dims, n_adj: int, dt: float, n: int, temp: float, ham: str = "KS") -> float:
    """(iii) Exact free-field error min_theta ||(U - e^{i theta} S^N) rho^{1/2}||_2 = sqrt(2 (1 - prod_k |f_k|)) over the
    n_adj (d - 1) transverse copies per k != 0. Reproduces Ch. 9's 0.0044 / 0.022 / 0.078 (3^3, n_adj 3, T a = 0.5,
    t = 10a, N = 450 / 200 / 107; tests)."""
    mult = n_adj * (len(dims) - 1)
    logf = sum(mult * math.log(free_gauge_mode_overlap(w, dt, n, temp)) for w in free_gauge_frequencies(dims, ham))
    return math.sqrt(2.0 * (1.0 - math.exp(logf)))


def trotter_check(a: Assumptions, o: dict) -> dict:
    """G7 / step rule for the 2033 box. STEP RULE (2026-10-04, Claude's decision, NEEDS_AUTHOR closed item): the
    state-dependent second-order estimate (ii) at eps = 0.1 sets the step; the exact free-field error (iii) is the
    check; the worst-case bound (i) is quoted once as the worst case. Evaluated for the deepest fit shot of each corner
    at that corner's temperature (19 steps over 0.18 fm/c at the fast end, T a_s = 0.30; 32 steps over 0.26 fm/c at the
    slow end, T a_s = 0.23, the rule's count, grid 26) and, as a record, for the 4 fm/c reference shot (400 steps of
    a_t/10, a pricing basis, not a run) at both temperatures."""
    hb = a.hbarc.lo
    L, d = int(a.L_2033.lo), int(a.dim_2033.lo)
    dims = (L,) * d
    links = n_links(L, d)
    dt = o["dt"] / a.a_s_fm.lo                               # 0.05 a_s, the a_t/10 grid
    temps = tuple(r["T"] * a.a_s_fm.lo / hb for r in o["r"])  # (0.304, 0.228): each corner's T a_s (was 0.162 at 160 MeV)
    eps, g2, nadj = a.trotter_eps.lo, a.trotter_g2.lo, int(a.trotter_n_adj.lo)
    cmax = a.sigma72_emax_over_fund.lo * 4.0 / 3.0            # s648: the Sigma(72x3) spectrum top, a proxy for Sigma(216x3)
    hop = a.n_stag.lo * a.n_c.lo / 2.0                        # staggered hop on one link: n_stag N_c / 2
    ref = round(a.t_max_2033_fm.lo / o["dt"])                 # 400
    hw, hd = ham_for(a, a.group_2033.value), a.trotter_dispersion_2033.value   # s648: H_I norms, improved dispersion
    out = {"trotter_dt_over_a_s": dt, "trotter_temp_lat": temps, "trotter_ref_steps": ref}
    # (steps, step in a_s, T a_s): the deepest shot of each corner as priced, and the reference shot at each corner's T
    deep = tuple((r["steps"][-1], r["grid_steps"][-1] * dt / r["steps"][-1], tp) for r, tp in zip(o["r"], temps))
    grid = tuple((r["grid_steps"][-1], dt, tp) for r, tp in zip(o["r"], temps))
    refs = tuple((ref, dt, tp) for tp in temps)
    for key, ns in (("deepest", deep), ("ref", refs), ("deepest_grid", grid)):
        wc = [trotter_worst_case(links, d, g2, cmax, hop, n * h, h, eps, hw) for n, h, _ in ns]
        wmin = [min(trotter_worst_case(links, d, x / 100, cmax, hop, n * h, h, eps, hw)["err"] for x in range(20, 501))
                for n, h, _ in ns]
        sd = [trotter_state_dependent(dims, nadj, tp, n * h, h, eps, hd) for n, h, tp in ns]
        out[f"trotter_worst_err_{key}"] = tuple(x["err"] for x in wc)
        out[f"trotter_worst_err_min_g2_{key}"] = tuple(wmin)
        out[f"trotter_worst_steps_{key}"] = tuple(x["n_steps_at_eps"] for x in wc)
        out[f"trotter_sd_err_{key}"] = tuple(x["err"] for x in sd)
        out[f"trotter_sd_steps_{key}"] = tuple(x["n_steps_at_eps"] for x in sd)
        out[f"trotter_free_err_{key}"] = tuple(free_gauge_trotter_error(dims, nadj, h, n, tp, hd) for n, h, tp in ns)
    out["trotter_steps_deepest"] = tuple(n for n, _, _ in deep)
    out["trotter_dt_deepest_fm"] = tuple(h * a.a_s_fm.lo for _, h, _ in deep)
    out["trotter_worst_lambda"] = trotter_worst_case(links, d, g2, cmax, hop, 1.0, dt, eps, hw)["lambda"]
    out["trotter_sd_lambda"] = tuple(trotter_state_dependent(dims, nadj, tp, 1.0, dt, eps, hd)["lambda"] for tp in temps)
    out["trotter_sd_lambda_ks"] = tuple(trotter_state_dependent(dims, nadj, tp, 1.0, dt, eps, "KS")["lambda"]
                                        for tp in temps)   # record
    out["trotter_rule"] = a.trotter_rule_2033.value
    out["trotter_rule_meets_eps"] = max(out["trotter_sd_err_deepest"]) <= eps        # the rule, on every priced shot
    out["trotter_meets_eps"] = max(out["trotter_free_err_deepest"]) <= eps           # the free-field check
    return out


def trotter_check_2028(a: Assumptions) -> dict:
    """Step rule for the 2028 box (2T on 4^2, free-field model with 3 colour copies, T a_s = 0.75): (ii) and (iii) for
    the six-step full benchmark at a_t/4 (0.15 fm/c) and the one-step first result at dt_first_2028_fm; and, as a
    record, the old Delta t = a_t (one and six steps)."""
    L, d = int(a.L_2028.lo), int(a.dim_2028.lo)
    dims = (L,) * d
    hb = a.hbarc.lo
    Tc = a.tc_over_sqrt_sigma.lo * a.sqrt_sigma_MeV.lo
    temp = a.T_over_Tc_2028.lo * Tc * a.a_s_fm.lo / hb        # 0.749
    eps, nadj = a.trotter_eps.lo, 3                           # SU(2): 3 colour copies
    dt = a.dt_over_at_2028.lo * a.a_t_fm.lo / a.a_s_fm.lo     # 0.125 a_s
    d1 = a.dt_first_2028_fm.lo / a.a_s_fm.lo                  # 0.225 a_s
    n = int(a.trotter_steps_2028.lo)
    lam = trotter_state_dependent(dims, nadj, temp, 1.0, dt, eps)["lambda"]
    old = a.a_t_fm.lo / a.a_s_fm.lo                           # 0.5 a_s
    return {
        "trotter_2028_temp_lat": temp,
        "trotter_2028_sd_lambda": lam,                                                    # 7.07 a_s^-3
        "trotter_2028_sd_err_full": lam * n * dt * dt * dt,                               # 0.083
        "trotter_2028_free_err_full": free_gauge_trotter_error(dims, nadj, dt, n, temp),
        "trotter_2028_sd_err_first": lam * d1 ** 3,                                       # 0.081
        "trotter_2028_free_err_first": free_gauge_trotter_error(dims, nadj, d1, 1, temp),
        "trotter_2028_dt_first_max_fm": (eps / lam) ** (1.0 / 3.0) * a.a_s_fm.lo,         # 0.048 fm/c
        "trotter_2028_sd_err_old": (lam * old ** 3, lam * n * old ** 3),                  # (0.88, 5.3) at Delta t = a_t
        "trotter_2028_free_err_old": (free_gauge_trotter_error(dims, nadj, old, 1, temp),
                                      free_gauge_trotter_error(dims, nadj, old, n, temp)),   # (0.77, 0.78)
        "trotter_2028_sd_steps_old_tmax": trotter_state_dependent(dims, nadj, temp, n * old, old, eps)["n_steps_at_eps"],  # 44
    }


def _model_2033(a: Assumptions) -> Result:
    c = _common_2033(a)
    o = own_depth_2033(a)
    breakdown = _own_breakdown(a, o)
    lo, hi = o["r"]
    hard_ops = (lo["t_shot_max"], hi["t_shot_max"])
    inter = _intermediates_2033(a, c)
    inter.update(_own_intermediates_2033(a, o))
    inter.update(_referee_intermediates_2033(a, o))
    inter.update(trotter_check(a, o))
    notes = _NOTES_2033 + _NOTES_2033_R25 + (
        "breakdown = one shot at a_s = 0.2 fm at the deepest fit point of the slow corner (T = 225 MeV, tau = 1/(2 pi T), "
        "c = 2 pi): 572 quench steps + 32 fit steps of the same primitives. hard_ops = (deepest shot at (300 MeV, 1/k, "
        "c = 2), deepest shot at (225 MeV, 1/(2 pi T), c = 2 pi)) at a_s = 0.2 fm = (6.255e9, 2.234e10), box "
        "'6.3e9-2.2e10', 6.3-22x the 1e9 budget; the second spacing's late shot (755 + 48 steps) is 3.0e10 "
        "(t_deepest_campaign). Was (1.051e10, 3.142e10) at T = 160 MeV, and (6.814e8, 3.006e9) with E-rho-OQ / 1e8 "
        "Gibbs prep before the quench ruling (2026-10-05).",
    )
    shots = tuple(sorted(inter["shots_campaign"]))
    return Result("2033", c["lq"], hard_ops, breakdown, inter,
                  shots=shots, wall_time_s=inter["wall_campaign_s"],
                  epsilon_l=(hi["eps_l"], lo["eps_l"]), notes=notes)


def _model_codesign(a: Assumptions) -> Result:
    """The 2033 instance with the chapter's own levers: BE trims evolution 10x
    (app03 route (a) / budget split) and the Ch.10 methods cut Gibbs prep 10x (app03 gap (state-prep bullet)).
    Not a separate box; the numbers are the ones printed at those lines."""
    c = _common_2033(a)
    # the 4 fm/c reference shot with the BE lever (no thermal preparation: the route (a) number is evolution only)
    breakdown = (
        Primitive(f"be_walk_{a.group_2033.value}", c["links"] * c["steps"], c["per_link_be"] / a.be_trim.lo,
                  CircuitStatus.CONJECTURE, "low2019hamiltonian,arxiv_2204_03381",
                  "'qubitized BE trims 10x' (app03 budget split); the Sigma(216x3) BE is not yet synthesized "
                  "(app03:85,91), so this is the Trotter per-link-step cost (gauge + FP hop + mass) divided by 10, "
                  "priced at the BE circuit's own R-TOL tolerance (N_rot / 10 rotations, eps_rot 2.53e-5)"),
    )
    inter = _intermediates_2033(a, c)
    inter.update({
        "t_per_shot_codesign_be": c["t_be"],                                             # 1.444e9 (BE at its own R-TOL eps)
        "eps_l_required_codesign": a.faults_per_shot.lo / c["t_be"],                     # 6.92e-11
    })
    notes = (
        "No codesign box in Ch. 5. This era is the 4 fm/c reference shot with the BE 10x trim printed at app03 route "
        "(a) ('~1.4e9 T/shot'); hard_ops is that evolution figure, with no thermal preparation in front of it.",
    ) + _NOTES_2033[:4]
    return Result("codesign", c["lq"], (c["t_be"], c["t_be"]), breakdown, inter,
                  shots=(c["shots"], c["shots"]),
                  epsilon_l=(inter["eps_l_required_codesign"], inter["eps_l_required_codesign"]),
                  notes=notes)


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
    "2028": Published(lq=(178, 178), hard_ops=(9.5e4, 5.9e5),
                      src="app03:2028 box (r25 apply, 2026-10-02): '178 = 32 links x 5 + 1-2 Hadamard + 16 workspace' (R10); "
                          "per-shot T 'first result 9.5e4 (one step, inside 1e5); full benchmark 5.9e5 (six steps)', each "
                          "timeslice at its own depth (R4). Exact 9.473e4 / 5.939e5; tol 0.05 covers the rounding and 177-178.",
                      rel_tol=0.05),
    "2033": Published(lq=(1234, 1334), hard_ops=(6.3e9, 2.2e10),
                      src="app03:2033 box (temperature ruling 2026-10-05: T = 1.5 T_c^lat = 225-300 MeV; quench thermal route "
                          "2026-10-05; s648: Sigma(216x3), 11 qubits/link, H_I, a_s = 0.2 fm; T1/T3 2026-10-05; step rule "
                          "2026-10-04; r25 2026-10-02): '~1280 (1234-1334)'; per-shot T, deepest shot at a_s = 0.2 fm: "
                          "'6.3e9-2.2e10 (5.6e9-2.1e10 quench + 6.9e8-1.2e9 fit evolution)', one circuit per shot, corners "
                          "(300 MeV, tau = 1/k, c = 2) and (225 MeV, 1/(2 pi T), c = 2 pi). hard_ops = (6.255e9, 2.234e10). "
                          "Was (1.051e10, 3.142e10) at T = 160 MeV. tol 0.05: two figures.",
                      rel_tol=0.05),
    "codesign": Published(lq=(1234, 1334), hard_ops=(1.4e9, 1.4e9),
                          src="NOT A BOX. app03 route (a) '~1.4e9 T/shot ... if it is 10x cheaper than Trotter' on the "
                              "400-step reference shot (s648: Sigma(216x3) / H_I; was 5.0e8 with Sigma(72x3) / H_KS). "
                              "Same register as 2033.",
                          rel_tol=0.05),
}


def INSTANCE_ROWS(a: Assumptions, era: str, r: Result):
    """Rows for resources.json: [(label, (lq_lo, lq_hi), (t_lo, t_hi), {extra}), ...]."""
    i = r.intermediates
    if era == "2028":
        return [
            ("2+1D SU(2)->2T, 4^2, pure gauge, first result: t = 0 and one step (inside the 1e5-T envelope)", r.lq,
             (i["t_first_result_shot_2028"], i["t_first_result_shot_2028"]),
             {"group": "2T", "lattice": "4^2", "shots": i["shots_first_result_2028"], "tier": "first result"}),
            ("2+1D SU(2)->2T, 4^2, pure gauge, full benchmark to t_max=0.15 fm/c (5.9x the 1e5-T envelope)", r.lq,
             (r.hard_ops[1], r.hard_ops[1]),
             {"group": "2T", "lattice": "4^2", "shots": i["shots_2028"], "conditional": True, "tier": "full benchmark",
              "condition": ">=5.9x per-step cut; the named levers give 3-10x (app03 2028 box, gap bullet)"}),
        ]
    if era == "2033":
        return [
            ("Sigma(216x3)+3 stag, H_I, 3^3, quench-prepared thermal state, deepest shot at a_s = 0.2 fm (quench "
             "ramp + c/T thermalization + Trotter to its own fit time; primitives in preparation, hop a Sigma(72x3) "
             "placeholder)",
             r.lq, i["t_deepest"],
             {"group": "S216x3", "hamiltonian": "I", "lattice": "3^3", "route": "quench (ETH)",
              "electric": "Sigma(216x3) FFT (Gustafson_S648_inprep)",
              "band": "corners (T = 300 MeV, tau = 1/k, c = 2) to (T = 225 MeV, tau = 1/(2 pi T), c = 2 pi) at a_s = "
                      "0.2 fm, T = 1.5 T_c^lat with T_c^lat = 150-200 MeV assumed; the second spacing (0.15 fm, own step) "
                      "late shot is 3.0e10 (t_deepest_campaign)",
              "condition": "ETH with dynamical fermions at the chapter's couplings on 3^3 and the relaxation of the momentum "
                           "density within t_th assumed; T_c^lat measured in situ; shots at the thermal variance; "
                           "non-freezing H_I trajectory assumed",
              "exports_note": "R9 exports (t_per_shot, t_depth_per_shot, f_star, floor_wall_s, factories_for_1yr) are the "
                              "campaign at the (300 MeV, tau = 1/k = 0.0955 fm/c, c = 2) corner only, over both spacings, "
                              "not a band (floor is the high/low depth schedule); wall_first_result_s and wall_campaign_s "
                              "are bands"}),
        ]
    if era == "codesign":
        return []     # r25: the 4 fm/c BE lever is not a landscape mark (each shot now runs to its own fit time)
    return []
