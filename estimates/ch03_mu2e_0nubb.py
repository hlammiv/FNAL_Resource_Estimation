"""Ch. 3 — Mu2e and 0nubb nuclear matrix elements. Reproduces the LQ and
hard-op numbers of applications/app11_mu2e_0nubb.tex.

The chapter's cost model (app11:78-83) is
    N_q     = 4 N_orb + N_anc                                        (eq:Nq11, JW, one register; R14)
    N_gate  = N_proj (c_proj tau_gs / dt) N_orb^4 c_T,  tau_gs = 1/Delta  (eq:Ngate11)
    c_T     = 1.15 log2(1/eps_rot) + 9.2,  eps_rot = sqrt(eps_syn / N_rot)   (R-TOL)
    N_shot  = (lambda/|<O>|)^2 (eps m)^-2 N_op N_iso                 (eq:Nshot11; R14: m = 1 Hadamard
                                                                      test, m = m_max depth-capped MLAE)
R14 OPTION A (author ruling 2026-09-30): a Hadamard test on a BE gives diagonal elements only. For 0nubb
  the parent is prepared, one BE query of O applied, and the result projected onto the daughter ground
  state with the same filter/QPE projection, in the SAME Fock register; the success probability is
  |M/lambda|^2. So N_proj = 2 for the pairs (parent preparation + daughter projection), exactly the
  count the survey already priced as 'N_reg = 2': no gate count moves. The register is one, not two:
  pf 40, jj44 44, jj55 64, Mo/Ru 108 system qubits (were 80, 88, 128, 216).
Every survey cell (sd, pf, jj44, jj55, jj44 x jj55 for Mo/Ru) and the harmonic-
oscillator ladder (e_max = 3, 4, 8-10) is eq:Ngate11 with a different N_reg, N_orb
and gap band.  The 2028 entry is an oracle-loaded insertion-only Hadamard test whose
Gamma <= 128 LCU insertion is counted gate by gate here (R12 item 13, insertion_cost()).

INPUTS AND SOURCES
  Delta = 1-2 MeV                     Stated app11:87  ("nuclear gap scales")
  Delta_large = 4 MeV                 Stated app11:150 ("optimistic Delta = 4 MeV gap")
  Delta_small = 0.5 MeV               Stated app11:156 (R11 ruling moru-small-gap (A), 2026-09-30)
  dt = 0.5 GeV^-1                     Stated app11:88
  c_proj in [pi, 2pi]                 Assumed app11:87 (Heisenberg-limited vs textbook projection)
  c_T per circuit (R-TOL, H. Lamm 2026-09-29): the full RUS fit 1.15 log2(1/eps_rot) + 9.2
                                      (Cited BRS 2015 + Campbell 2017, randomized) at
                                      eps_rot = common.eps_rot_for(N_rot), N_rot = N_proj N_Trotter N_orb^4
                                      rotations in ONE shot of that circuit; every band end, gap
                                      variant and lever variant is its own circuit with its own eps.
                                      Replaces the fixed c_T = 27 (was Cited, app11:90).
  N_orb: sd 6, pf 10, jj44 11, jj55 16, Mo/Ru 27=11+16, He-4 10   Stated app11:90,131,156,122
  4 qubits per spatial orbital        Stated app11:79 (eq:Nq11, JW second quantization)
  HO ladder: (e+1)(e+2)(e+3)/6 spatial orbitals; chapter quotes 165 (286) at e_max=8 (10)
  N_anc = 50-150 (2033), 80-130 (2028)  Stated app11:85, app11:126
  eps = 0.05 (Mu2e), 0.1 (0nubb), 0.1-0.3 (2028)  Stated app11:92, app11:122
  eps_l = 0.1 / N_gate^shot (0.1 expected faults per shot)  Stated app11:93 (ruling R3, 2026-09-28)
  eps_synth_total = 1e-2 per shot     common.EPS_SYN (R-TOL, report-wide), stated app11:90
  RFI envelopes 150-250 LQ / 1e5 T (2028), 1000 LQ / 1e9 T (2033)  Cited DOE_RFI_2026
  t_gate_s = 1 us per T-gate          Assumed; report-wide wall-time convention, ruling R8
                                      (H. Lamm, 2026-09-29), stated in Ch. 1; app11:127,185.
                                      Was 10 us per logical cycle before R8. Moves time only. R9 name
                                      (was logical_cycle_s).
  shot_overhead_s = 0.1 ms per shot   Assumed; report convention (R9, 2026-10-02; was 1 ms, r17).
ROUNDING (ruling R4, 2026-09-28)
  Every product carries the exact inputs (N_Trotter = 1000pi-4000pi, c_T at its exact eps_rot) and is rounded
  once, to two significant figures, at print.  `printed(a)` returns every printed string the
  chapter quotes so the prose and the model come from one place; the tests compare them.
2028 INSERTION, DERIVED HERE (R12 item 13; R13 'one phase' ruling, 2026-09-30; logs
  apply_log/r12_ch03.md, apply_log/r13_ch03.md)
  The target is <Psi_0|O|Psi_0> on one oracle-loaded state (box app11:172), O Hermitian. A Hadamard
  test on the controlled BE U_O gives Re<0,Psi_0|U_O|0,Psi_0> = <Psi_0|O|Psi_0>/lambda exactly: the
  degree-1 polynomial, ONE query, no QSP phases (GSLW 1806.01838 def:phaseSeq with n = 1 and phi = 0;
  Martyn 2105.02859 main.tex:152 'phi = (0,0) ... U_phi = W(a) is just the unchanged signal'), and no
  real-part step (the expectation value is real). eps sets the shots, not the circuit. The pre-R13
  text multiplied T_BE by ceil(log2 1/eps) = 2-4 'QSP phases'; that O(log 1/eps) overhead belongs to
  approximating a function of O (e.g. e^{-iOt}, GSLW thm:blockHamSim), not to O itself.
  Definitions printed at app11:122: a query is one application of U_O (or U_O^dag); a phase is one
  projector-controlled rotation e^{i phi (2 Pi - 1)}; degree d = d queries and d + 1 phases.
  Per query: PREPARE + PREPARE^dag rotation tree, 2 x 127 R_y; SELECT controlled unary iteration,
  127 Toffolis (Babbush PRX 2018, L-1 ANDs); Pauli targets Clifford. Gamma = 128: 254 rotations,
  127 Toffolis; R-TOL eps_rot 6.27e-3, c_T 17.61; T = 4473.9 + 889 = 5362.9, 0.054x the cap.
  Same count: Gamma = 256 0.11x, Gamma = 512 0.23x; full one-body Gamma = 1600 (real, 2 C(40,2) + 40)
  0.93x, 3160 (complex, 4 C(40,2) + 40) x1.9.
  Sensitivity (not printed as the headline): the literal pre-R13 reading, d = 2-4 queries plus d + 1
  phase steps (k - 1 flag ANDs + 1 R_z each), 11201-22928 T (0.11-0.23x); doubled (W and W^dag, or a
  real-part LCU of Phi and -Phi), 0.23-0.47x. Intermediates qsp_reading_*.
  Insertion-stage register 40 + 7 + 7 + 1 = 55 LQ (no phase ancilla); load stage 40 + 11 address + 10
  workspace = 61 (R14: address sized to the printed D <~ 2e3; R13 used the 3-pass crossing D = 4508 -> 65);
  peak 61 at the load; 116 at the MLAE reflection (S_0 on 40 system + 11 load address + 7 LCU address + 1 test
  = 59 qubits + 57 AND ancillae; 94 before the R14 fix round omitted the load address).
  Wall per shot = insertion + load at D ~ 1e2-1e3: 0.0067-0.026 s.
R14 (author rulings 1a, 1d, 2 option A, 3; log apply_log/r14_ch03.md)
  1a  Dipole operator = proton charge monopole sum_i (1+tau3_i)/2 j0(q r_i), q = m_mu (Haxton et al.
      2208.07945). HO 0s,0p,1s0d (10 orbitals, 40 modes), b = 1.369 fm from r_ch(4He) = 1.67824 fm, r_p =
      0.8409 fm, c.m.-corrected 0s^4 radius. Gamma = 25 JW strings (I, 20 Z, 4 XX/YY), lambda = 15.39,
      <O>_0s^4 = 1.748, (lambda/<O>)^2 = 77.47; any 2-proton state in the space: 69.7-127.8.
      Hadamard-test shots 100 x 77.47 x 6 = 4.65e4; wall 314-1224 s (5-20 min).
  1d  Depth-capped MLAE: A = 2 load passes + 1 controlled query; m = 2k+1 applications; S_0 = 46 Toffolis.
      m_max = 13 (D=1e2) and 5 (D=1e3), 1 at D = 2e3 (m = 3 is 1.011x). Circuits (lambda/<O>)^2/(eps m)^2
      x 6 = 275-1859, wall 26.6-184 s. Estimator constant 1 (Suzuki Cramer-Rao, all at m_max).
R15 (author rulings 1-4, 2026-10-01; log apply_log/r15_ch03.md)
  1  Identity subtraction: c_0 = Tr(o)/2 is classically known; the LCU is O - c_0 I. 4He dipole: Gamma 24,
     lambda' = 7.79; the test estimates (<O> - c_0)/lambda' with per-shot variance <= 1, so relative eps on
     <O> needs (lambda'/(eps |<O>|))^2 = 19.85/eps^2 (17.9-32.8 any state). Shots 1.2e4 (was 4.6e4).
     R16: shots at the upper end of the any-state band (32.76, as the 27Al cell): 2.0e4, 7.9e2 MLAE circuits.
     27Al sd: lambda' = Tr/2 = 3.37, factor 1.45, 3.5e3 shots per scan (was 1.4e4).
  2  Oracle load = the complete sparse load of Fomichev et al. 2310.18410, (2L-2)D + 2^(L+1) + D Toffolis,
     L = ceil(log2 D), ONE per shot (it erases its own index; the test measures). Fits to D = 603, first crossing D = 604,
     printed D <~ 600; load-stage register 40 + 5L - 3 = 87 at L = 10; MLAE S_0 on 40 + 10 + 7 + 1 = 58,
     peak 114; m_max = 5 at D = 1e2 (to D = 128), 3 to D = 228, then 1.
  3  Mo/Ru = jj44 protons (22) + jj55 neutrons (32) = 54 modes, N_orb 13.5: x6.0-25 generic, x3.0-6.0
     large gap, x25-51 small gap; register 54 (104-204 LQ).
  R16 (author ruling, 2026-10-01): the jj44 pairs (76Ge/Se, 82Se/Kr) are co-design, not near miss: their
     prep alone reaches x10.8 at the top of the band, across the ~10x co-design line (main-overview:106).
     The near miss is the pf pair alone (x1.8-7.3). Mo/Ru (x6-25) keeps its role label 'reach target'
     and sits between the jj44 and jj55 pairs, co-design at the generic gap band, a near miss at a large
     gap. In resources.json the jj44 row moves to the codesign era with jj55 and Mo/Ru (drawn on the
     main axes by its T, < 1e12, as resources.py derives the drawing flag).
  4  2033 insertions: one controlled BE query, Gamma from the selection-rule count (M, parity, common pair
     J) by explicit JW algebra: jj44 0nubb 2501 terms, Gamma 40016, 4.0-4.1e6 T (0.04-0.15% of the pair
     prep); dense 53361 terms, 6.5-6.8e7 T. 27Al sd Mu2e (general scalar 1+2-body, charge conserving):
     Gamma 12900 / 49140 dense, 9.3e5-3.9e6 T, in the prep circuit (shared R-TOL tolerance).
REFEREE PASS (2026-10-04; editorial_review/responses/ch03.md)
  G2: Result.hard_ops for 2028 is the executed Hadamard-test shot, insertion (5362.9) + one complete load at
    D = 1e2 .. 600: 16254.9-99498.9 T, printed 1.6e4-9.9e4 (was the insertion alone, 5.4e3). Walls unchanged
    (they already carried the load). PUBLISHED 2028 (1.6e4, 9.9e4); Table 1.1 cell follows (shared edit).
  M3: load_angle_cost derives the load's amplitude-angle rotations (LKS 1812.00954 + Gidney 1709.06648), b = 13:
    +765-1017 T, D = 600 at 1.005x, fit to D = 596. Author ruling (2026-10-04, "do 3", Claude's recommended
    default): FOLDED IN. hard_ops 17020.3-100516.3 (printed 1.7e4-1.0e5, no longer a lower bound; conditional
    on D), T-depth + L(b-1) + one c_T layer, walls, eps_l, MLAE (load_angles=True everywhere), register: load
    stage 40 + 10 + 38 = 88, MLAE 114 + 13 reused phase-gradient qubits = 127. PUBLISHED 2028 (1.7e4, 1.0e5).
    The no-angle values are kept as *_no_angles records.
  Verifier (2026-10-04): mlae_2028(load_angles=True) charges the phase-gradient state once per circuit
    (catalytic): MLAE deepest 83754 -> 86893 (8.7e4, was 8.76e4 with the state per application); m_max = 5
    to D = 128 (unchanged), 3 to D = 221 (228 without), then 1 to D = 596 (mlae_m_max_last_d_with_load_angles).
  Max subroutine time: 0.10 s per 2028 shot, 1.8-7.3 min per 2033 27Al shot, ~16 min per f7p3p1 resonant-drive
    shot (max_subroutine_time_s; the requirements row no longer says N/A).
REFEREE M3 FOLLOW-UP (2026-10-05; scratch chopen/app11, exact diagonalization, interactions from KSHELL snt files)
  0+ gaps of the truncated Hamiltonians, KB3G restricted to the kept orbits (48Ca / 48Ti): f7p3 6.100 / 3.685,
    f7p3p1 5.501 / 3.604, pf 5.175 / 4.265 MeV, replacing ENSDF 4.283 / 2.997 in simplified_0nubb (GXPF1A, the
    record, gives a smaller 1/Delta sum in every space). f7p3 projection rung 1.19e8-2.44e8 -> 9.10e7-1.86e8 T,
    20-40 -> 15-31 d, eps_l 5.4e-10-1.1e-9. pf at its 0+ gaps 9.82e8-2.00e9 -> 7.35e8-1.50e9 T (x0.74-1.5, 23-68 yr). The ENSDF-gap
    costs are kept as t_at_ensdf_gaps records. 27Al: generic band kept; USDB gap 2.328 MeV (7/2+, the nearest
    state an M = 5/2 reference reaches) is a record, al27_prep_t_at_usdb_gap 8.93e7-1.83e8 T.
  D for 4He (d_he4_computed): Minnesota NN + Lawson (beta 2, 10) in the 40-mode space: fidelity 0.99 at
    D = 15-17, 0.999 at 77-103, dipole to 1% at D = 1 (at the chapter's b). At or below the priced band D = 1e2-600; headline unchanged.
WHAT IS NOT DERIVED HERE
  QROM oracle load D-1 Toffolis per pass  Cited Babbush_PRX_2018, priced at 7 T/Toffoli (R5; R11b
    ruling qrom-4T-vs-R5 (ii), was the source's 4D-4 temporary-AND count); 2-3 passes per shot
    Stated app11:124 (R11 ruling qrom-passes (A)). Priced condition D <~ 2e3 (R13, one-query insertion):
    5362.9 + 2 x 7 x 1999 = 33349 (0.33x) to 5362.9 + 3 x 7 x 1999 = 47342 (0.47x). The crossing is
    derived: D = 1 + 94637/14 = 6760.8 (2 passes) and 1 + 94637/21 = 4507.5 (3 passes), printed
    "D ~ 4.5-6.8e3" (R12 3.7-6.3e3 at 2-4 queries; R11b 2.0-3.3e3 at 5.4-5.7e4).
  (R15: the 2033 insertions are derived; the pre-R15 Stated jj44 5.4-5.6e7 / 4.4-4.6e8 is a record only.)
  Top-down T_BE ~ 1e3-1e4 for the 2028 BE      Stated app11:122.
  2028 headline vs oracle load (R11 ruling 2028-headline-oracle-load, Part 1 + (A)): hard_ops is
    the Gamma<=128 insertion, printed 5.4e3 since R13 (one query; 1.1-2.3e4 in R12); the oracle load
    is printed separately in the box, 2 x 7 x 99 = 1386 at D=1e2 to 3 x 7 x 999 = 20979 at D=1e3
    ("1.4-21e3"; R11b, was 0.8-12e3 at 4 T), and carried as oracle_load_t_at_d_plausible.
    hard_ops_with_oracle_load_at_d (3.3-4.7e4 at D=2e3, R13) is the priced-condition row.
  Shots (R11 ruling 2033-shots (A)): the 2033 box counts 27Al only, 6 x 400 (Mu2e, eps=0.05)
    = 2400 per basis scan, 4800-7200 over 2-3 scans, times the sd dipole factor 1.45 (R15).
    r17: Result.shots spans the 2-3 scans (7.0e3-1.0e4; was one scan to three) and the wall is
    8.5 days-1.8 months; t_shot carries the 1 ms per-shot overhead (shot-time model 1A). The 0nubb
    elements run on their own (out-of-envelope) circuits. The old box reading (100-400) x 20
    = 2-8e3 is kept as shots_box_reading_pre_r11 for the record.
  One machine (r23, author ruling 2026-10-02, apply_log/r23_ch03.md; numbers superseded by r25, see WALL TIME
    AND FACTORIES): the report assumes a single quantum
    machine. Every wall time here is shots x per-shot time, run one after another on that machine
    (1 us per T-gate + the per-shot overhead). The model never carried a machine count; the chapter's
    "divides across shot-parallel machines" sentence and the box's "shot-parallel" / "serialized" labels
    are retired. The 2033 27Al campaign (2-3 scans) is 7.358e5-4.612e6 s = 8.5 days-1.8 months, which is
    0.47-2.9% of the 5-year campaign horizon (campaign_horizon_yr, YEAR_S = 365.25 d): it fits on one
    machine with no reduction needed. Intermediates wall_fraction_of_horizon,
    wall_fits_horizon_on_one_machine; the chapter prints "at most 2.9%".
  Flagship (R11 ruling flagship-fits-1000 (A)): the chapter now says the register sits at the
    1000-LQ envelope, inside at e_max=8 (660) and 14% over at e_max=10 (1144) before ancilla.
  Mo/Ru small-gap x430-880 is Delta_small = 0.5 MeV, now printed (app11:156).
  KNOWN GAPS (see NEEDS_AUTHOR.md): the load-stage workspace is derived since R12 (log D - 1 = 11,
    output into the system register). The Gamma<=128 costing and the insertion-stage
    peak are derived here since R12. The QROM price follows R5 since R11b (qrom-4T-vs-R5 (ii)).
  CLOSED by the 2026-09-28 pass (R3/R4 + mechanical batch): every rounding-chain corner
    (pf 6.8e9, jj44 9.9e9, e_max=3 5.4e10, e_max=4 x130-510, jj55 x5.6, Mo/Ru x45-90 and
    x360-720, lever 0.99 and 4.5-72, flagship 2.3e15), the eps_l rows at 0.1 fault/shot,
    "about 1% of the jj44 pair", the synthesis budget 1e-2.
R-TOL (H. Lamm, 2026-09-29; log apply_log/rtol_ch03.md)
  The fixed c_T = 27 is gone. Every circuit (each band end of each cell, each gap variant, the
  x10 lever variant) carries its own eps_rot = sqrt(1e-2 / N_rot), N_rot = N_proj N_Trotter N_orb^4,
  and c_T = RUS fit at that eps: sd 25.65-26.80, pf 27.92-29.07, jj44 28.23-29.38, jj55
  29.48-30.63, Mo/Ru 31.21-32.36, flagship 36.64 (e_max=8 lo)-39.62 (e_max=10 hi). The 2028 entry
  had no rotation breakdown then (it has one since R12, priced by R-TOL) and the stated sd
  insertion has no rotation split; neither moves. The lever worst case (jj44 hi / 10, its own eps)
  is now 1.01x the limit, just over (was 0.99x).
WALL TIME AND FACTORIES (r25: rulings R1, R7, R9, R10, H. Lamm 2026-10-02; log apply_log/r25_ch03.md)
  t_gate_s = 1e-6 s per T, shot_overhead_s = 1e-4 s per shot, ONE machine, shots run one after another.
  Wall per shot = max(N_T t_gate_s, D_T x 10 us) + shot_overhead_s with the 10-factory baseline: the T-supply
  time when F* = N_T / D_T >= 10, the depth-limited (corrected) time otherwise (wall_one_machine).
  T-depth per shot from factory.json (Ch. 3), derived here:
    2028 Hadamard test (headline run), D = 1e2 .. 600: N_T 16255-99499, D_T 2382-6902 (compiled: PREPARE
      beside the load, pipelined rotation tree 2^(k-1) c_T, AND-tree erasure); ladder 2810-14702, serial
      6157-18049. F* per end 6.8 (D = 1e2, below 10: depth-limited) and 14.4. Wall 470-1958 s = 7.8-33 min
      (was 339-1975 s at 1 ms and the T-supply time). MLAE at D = 1e2: D_T 7923, F* 10.6, 786 x 0.0839 s = 66 s.
      Single tier: wall_first_result_s = None.
    2033 27Al (prep + insertion, N_T 1.05e8-4.40e8): D_T = N_rot c_T / G + insertion tail, G = 12 commuting
      rotations in flight (trotter_concurrency; 12-24 of the 50-150 ancillas): 8.9e6-3.7e7, F* 11.8 at both
      ends. Gate by gate (G = 1) F* = 1 and the campaign floor is 135 d-3.0 yr.
  Readout (shot audit): the bare one-body SI elements are pinned to +/-0.61% in any state, so they are a
  classical validation; the SD element <O_SD> = 2 <S_p> o_0d = 0.335-0.380 sets the shots. Direct readout after
  a single-particle basis rotation, var ~1 (estimate; Popoviciu worst case x8.3), eps 0.05, heralded
  projection success p = 0.5:
  2.77e3-3.56e3 / p = 5.5e3-7.1e3 shots per Hamiltonian (a Hadamard test would carry 79-101 x eps^-2).
  First result (R1, Ch. 3 ruling: SD at 5%, one Hamiltonian): 5.5e3-7.1e3 shots -> wall_first_result_s
    5.84e5-3.13e6 s = 6.8-36 d. At R1's generic 30%: 154-198 shots, 4.5-24 h (no category separation).
  Campaign (2-3 basis-cutoff scans, six categories at ~10% on CR): 1.1e4-2.1e4 shots -> wall_campaign_s
    1.17e6-9.40e6 s = 13.5 d-3.6 mo, at most 6.0% of the 5-year horizon; p over 0.1-0.9: 7.5 d-1.5 yr.
    Exports (campaign run): floor 11-92 d, factories_for_1yr 0.37-3.0.
  R7 rows (simplobs.json; gaps superseded by the referee M3 KB3G gaps, see above): 48Ca/48Ti first light in
    f7/2 p3/2 (24 modes, 74-174 LQ) at the 0+ gaps 4.283 / 2.997 MeV, projection readout, 30% on |M| (lambda/|M| = 69.1, p_ref 0.93, toy): 1.2e8-2.4e8 T, 1.4e4 shots, 20-40 d
    (G = 12 in flight: D_T 1.0e7-2.0e7, F* 11.9-12, so the T-supply wall stands).
    Then f7/2 p3/2 p1/2 (28 modes, 78-178 LQ) with the resonant-drive (Rabi) readout: 4.4e8-7.2e8 T, 250-1000
    shots. Full pf by Rabi: 2.3e9-3.7e9 T (x2.3-3.7, near miss). Rabi T-depth (R10 verifier fix):
    D_T = N_trotter_rot c_T / G + N_qdrift c_T, since the qDRIFT samples are random, mostly non-commuting
    single-Pauli rotations, one sequential layer each. f7p3p1: D_T 6.9e7-9.3e7, F* 6.4-7.7, wall 2.0-11 d
    (T supply alone 1.3-8.3 d). pf: D_T 6.5e8-9.8e8, F* 3.5-3.8, wall 19-110 d (T supply alone 6.6-43 d).
    eps_l (0.1 faults per shot): f7p3 4.1e-10-8.4e-10, f7p3p1 Rabi 1.4e-10-2.3e-10, pf Rabi 2.7e-11-4.4e-11.
    All classically tractable validations. 76Ge stays past 2033 (no simplification that keeps its physics fits).
  The 2033 27Al per-shot T keeps the BE insertion (~1% of N_T), which the direct readout does not use:
    conservative, left in.
"""

from __future__ import annotations

import itertools
import math
from dataclasses import dataclass, fields
from functools import lru_cache as _lru_cache

from estimates.common import (Assumed, Cited, Stated, Uncited, Tagged, Primitive,
                              CircuitStatus, Result, Published, ERAS, SYNTHESIS_SRC,
                              t_per_rotation, eps_per_rotation, eps_rot_for, EPS_SYN, EPS_SYN_SRC,
                              T_PER_TOFFOLI, depth_exports, REACTION_TIME_S, FACTORY_BASELINE)

TEX = "app11"
RFI = "DOE_RFI_2026"
# record only: the pre-R15 stated jj44 insertion (no derivation found; R14 log section 2)
INSERTION_JJ44_PRE_R15 = dict(symmetry_restricted=(5.4e7, 5.6e7), dense=(4.4e8, 4.6e8))
PI = math.pi
YEAR_S = 365.25 * 86400.0          # Julian year; the campaign horizon is campaign_horizon_yr x YEAR_S (r23)


@dataclass(frozen=True)
class Assumptions:
    # ---- RFI anchors -----------------------------------------------------
    rfi_lq_2028: Tagged = Cited((150, 250), RFI, "2028 register envelope (app11:169)")
    rfi_t_2028: Tagged = Cited(1e5, RFI, "2028 hard-op cap (app11:169)")
    rfi_lq_2033: Tagged = Cited(1000, RFI, "2033 '1000+ LQ' (app11:130)")
    rfi_t_2033: Tagged = Cited(1e9, RFI, "2033 hard-op envelope (app11:130)")

    # ---- state-prep depth, eq:Ngate11 (app11:82, 87-90) ---------------------
    gap_mev: Tagged = Stated((1.0, 2.0), f"{TEX}:87", "'nuclear gap scales Delta ~ 1-2 MeV'")
    gap_large_mev: Tagged = Stated(4.0, f"{TEX}:150", "'an optimistic Delta = 4 MeV gap'")
    gap_small_mev: Tagged = Stated(0.5, f"{TEX}:156", "'at a Delta = 0.5 MeV small-gap end, the regime "
                                                       "Mo's shape-coexistence physics argues for' (R11 ruling "
                                                       "moru-small-gap (A), 2026-09-30)")
    dt_gev_inv: Tagged = Stated(0.5, f"{TEX}:88", "'at the dt ~ 0.5 GeV^-1 step'")
    n_trotter_quoted: Tagged = Stated((3.1e3, 1.3e4), f"{TEX}:88",
                                      "'N_Trotter = 1000pi-4000pi ~ 3.1e3-1.3e4' as printed; every product "
                                      "below carries the exact 1000pi-4000pi (ruling R4)")
    c_proj: Tagged = Assumed((PI, 2 * PI), "'c_proj in [pi, 2pi] ... Heisenberg-limited phase "
                                           "estimation (pi) and the standard projection prescription (2pi)'",
                             f"{TEX}:87")
    ch2_response_horizon_gev_inv: Tagged = Stated(100.0, f"{TEX}:87",
                                                  "implied by 'i.e. 5-10x the response horizon of Ch. 2'")
    # R-TOL (2026-09-29): no fixed c_T. Each circuit's c_T is the RUS fit at eps_rot_for(N_rot).
    eps_synth_total: Tagged = Stated(EPS_SYN, f"{TEX}:90", "'a total synthesis error of 1e-2 per shot' "
                                                           "(R-TOL, report-wide; " + EPS_SYN_SRC + ")")

    # ---- registers, eq:Nq11 ------------------------------------------------
    qubits_per_orb: Tagged = Stated(4, f"{TEX}:79", "'2 . 4 N_orb + N_anc (JW second quantization)'")
    n_orb_he4: Tagged = Stated(10, f"{TEX}:122", "'N_orb ~ 10 active space of Ch. 2'")
    n_orb_sd: Tagged = Stated(6, f"{TEX}:90", "'27Al has N_reg=1, N_orb=6'; 24 system qubits")
    n_orb_pf: Tagged = Stated(10, f"{TEX}:131", "'~80 for the 48Ca/48Ti pair (pf shell)' = 2.4.10")
    n_orb_jj44: Tagged = Stated(11, f"{TEX}:156", "'jj44 x jj55 = 11 + 16 = 27 orbitals'")
    n_orb_jj55: Tagged = Stated(16, f"{TEX}:156", "'jj44 x jj55 = 11 + 16 = 27 orbitals'")
    emax_first: Tagged = Stated((3, 4), f"{TEX}:152", "'The first enlargement, e_max=3' ... 'At e_max=4'")
    emax_converged_al: Tagged = Stated((8, 10), f"{TEX}:135", "'An e_max = 8 (10) harmonic-oscillator basis'")
    emax_converged_ge: Tagged = Stated((10, 14), f"{TEX}:137", "'Spaces at e_max = 10-14'")
    mass_number_ge: Tagged = Cited(76, "LEGEND1000", "A = 76 for the first-quantized line of eq:Nq11 (app11:80)")
    ho_orbitals_quoted: Tagged = Stated((165, 286), f"{TEX}:135", "'carries 165 (286) spatial orbitals'")
    n_anc: Tagged = Stated((50, 150), f"{TEX}:85", "'The ancilla budget N_anc ~ 50-150'")
    n_anc_2028: Tagged = Stated((80, 130), f"{TEX}:126", "'plus ~80-130 ancilla for the BE workspace'")
    # R12: uncontrolled unary iteration over D amplitudes uses log D - 1 AND ancillae (Babbush PRX 2018,
    # arXiv:1805.03662 main_draft.tex:636 'only log L ancillae' for the controlled form; the uncontrolled form
    # drops the control AND). The lookup writes each configuration into the system register (no separate
    # output register). Was Stated(2) '+ 2 QROM work at the load stage'.
    qrom_and_deficit: Tagged = Cited(1, "Babbush_PRX_2018", "'only log L ancillae' (main_draft.tex:636) for "
                                                            "controlled unary iteration; uncontrolled: log L - 1")
    test_qubits: Tagged = Stated(1, f"{TEX}:127", "'+ 1 test qubit at the insertion stage'")

    # ---- precision and shots, eq:Nshot11 ----------------------------------
    eps_2028: Tagged = Stated((0.1, 0.3), f"{TEX}:122", "'With eps ~ 0.1-0.3 on a single matrix element'")
    eps_mu2e: Tagged = Stated(0.05, f"{TEX}:92", "'~10% on CR means eps ~ 0.05 on the matrix element'")
    eps_0nubb: Tagged = Stated(0.1, f"{TEX}:92", "'Beating the current factor-2-3 spread needs eps ~ 0.1'")
    n_op_mu2e: Tagged = Cited(6, "Cirigliano_Kitano_Okada_Tuzon", "six operator categories (app11:8)")
    n_op_0nubb: Tagged = Stated(5, f"{TEX}:233", "'5 0nubb operators'")
    n_pairs_cross_check: Tagged = Stated(3, f"{TEX}:233", "'x 3 pairs'")
    n_basis_scans: Tagged = Stated((2, 3), f"{TEX}:233", "'x 2-3 basis-cutoff scans'")
    n_survey_elements: Tagged = Stated(20, f"{TEX}:225", "'~20 operator-isotope elements of the 2033 "
                                                          "survey'; 6 + 5x3 = 21")
    # R9 (H. Lamm, 2026-10-02): unified name t_gate_s (was logical_cycle_s). 1 us per T-gate, the report-wide
    # wall-time convention of ruling R8 (2026-09-29); was Stated 10e-6 'at the 10 us logical cycle of Ch. 2'.
    t_gate_s: Tagged = Assumed(1e-6, "1 us per T-gate (hard operation), the report-wide wall-time convention "
                                     "(R8, H. Lamm, 2026-09-29; R9 name t_gate_s, was logical_cycle_s)",
                               f"{TEX}:127,185")
    # r17 (shot-time model 1A): a shot of N_T T-gates takes N_T x t_gate_s plus a fixed per-shot overhead for
    # register initialization, final readout and decoding. R9 (author approved, 2026-10-02): the report value
    # 0.1 ms on current surface-code hardware (1.1 us cycle, 63 us decoder latency, Google_QEC_below_threshold);
    # was the 1 ms bound of r17. Moves the 2028 Hadamard-test walls only (339-1975 s -> 470-1958 s with the
    # T-depth correction at D = 1e2); invisible at 2033 depth.
    shot_overhead_s: Tagged = Assumed(1e-4, "per-shot init + readout + decode, 0.1 ms, the report convention "
                                            "(R9, 2026-10-02; was the 1 ms bound of r17)",
                                      "Google_QEC_below_threshold")
    # r23 (author ruling 2026-10-02): one quantum machine everywhere. A campaign is its shots run one after
    # another on that machine and is compared with the 5-year campaign horizon of the requirements table.
    campaign_horizon_yr: Tagged = Stated(5, f"{TEX}:233", "'Campaign horizon & 5 years' (requirements table); "
                                                          "the 2033 wall is printed as a share of it, app11:149")
    fault_budget: Tagged = Stated(0.1, f"{TEX}:93", "'the quoted eps_l is 0.1/N_gate^shot: ... 0.1 faults "
                                                    "are expected per shot' (ruling R3, report-wide)")

    # ---- 2028 insertion-only benchmark (app11:122-127) ---------------------
    t_be_topdown: Tagged = Stated((1e3, 1e4), f"{TEX}:122", "'T_shot^2028 ~ T_BE ~ 1e3-1e4' (R13; was "
                                                            "'T_BE . ceil(log2 1/eps) ~ (1e3-1e4).(2-4)')")
    # R13 ('one phase' ruling): the Hadamard test on the controlled BE U_O returns
    # Re<0,Psi|U_O|0,Psi> = <Psi|O|Psi>/lambda, the degree-1 polynomial: one query, no QSP phases
    # (GSLW arXiv:1806.01838 def:phaseSeq, n = 1; Martyn arXiv:2105.02859 main.tex:152). Derived here.
    be_queries_2028: Tagged = Assumed(1, "derived here (R13): Hadamard test on one controlled BE query gives "
                                         "<Psi_0|O|Psi_0>/lambda exactly; O Hermitian, so no real-part step",
                                      f"{TEX}:122,124")
    gamma_trunc: Tagged = Stated(128, f"{TEX}:124", "'a Gamma <= 128-term insertion'")
    gamma_256: Tagged = Stated(256, f"{TEX}:124", "'0.11x at Gamma = 256' (R13; 0.23-0.47x in R12)")
    gamma_512: Tagged = Stated(512, f"{TEX}:124", "'0.23x at Gamma = 512' (R13; 0.47-0.96x in R12)")
    # R12: the full one-body Gamma is derived in _model_2028 from the 40 JW modes: real 2 C(40,2) + 40 = 1600
    # Pauli strings, complex 4 C(40,2) + 40 = 3160 (was Stated (2000, 3200), no count behind the low end).
    # ---- the Gamma-term insertion, counted gate by gate (R12 item 13; derived here, app11:124) ----
    # SELECT: controlled unary iteration over L = Gamma terms is L-1 AND computations (Toffolis), using
    # log L AND ancillae, the temporary ANDs uncomputed by measurement (Babbush PRX 2018, main_draft.tex:
    # 636-637, 716; Fig. 7 caption :701). The source's own T price is 4 per AND; R5 prices 7.
    select_and_deficit: Tagged = Cited(1, "Babbush_PRX_2018", "'always ends up with L-1 AND computations' "
                                                              "(arXiv:1805.03662 main_draft.tex:716); the L-1 "
                                                              "count includes the control input")
    insertion_t_per_toffoli: Tagged = Assumed(T_PER_TOFFOLI["textbook"], "R5: 7 T per Toffoli everywhere "
                                              "(was 4 in the source's temporary-AND count)",
                                              "TRACKED_CHANGES.md R5")
    qrom_t_per_amplitude: Tagged = Cited(T_PER_TOFFOLI["textbook"], "Babbush_PRX_2018",
                                         "'D-1 Toffolis (7(D-1) T at 7 T per Toffoli) per QROM pass'; D-1 from "
                                         "the source, 7 T/Toffoli from R5 (R11b ruling qrom-4T-vs-R5 (ii); was 4)")
    d_amplitudes: Tagged = Stated(600, f"{TEX}:126", "'the total stays under the cap for D <~ 600' (R15: the "
                                                     "complete load, fits to D = 603 (first crossing 604); R11b-R14 2e3 "
                                                     "lookup-only)")
    # d_break is no longer an input: the crossing D is derived in era_2028 (R11b ruling 1; was Stated(1e4))
    d_plausible: Tagged = Stated((1e2, 1e3), f"{TEX}:124", "was 'plausibly needs D ~ 1e2-1e3' (retired print, "
                                                           "referee M3); the low end is the priced band's lower "
                                                           "end, D = 1e3 the record of composed_fraction_at_d_hi")
    # Referee M3 (2026-10-05): D computed for a schematic 4He ground state in the 40-mode space: Minnesota NN
    # (Thompson, LeMere, Tang, NPA 286 (1977) 53; u = 1, no Coulomb) + intrinsic kinetic energy at the chapter's b,
    # Lawson c.m. term beta = 2 and 10 (without it the ground state is c.m.-contaminated: (0s)^4 weight 0.83,
    # D 125 / 299 at fidelity 0.99 / 0.999), exact diagonalization in M_L = M_S = 0 (dim 2556). D = number of
    # largest-|c| JW configurations kept. E0 = -26.9 / -26.8 MeV; (0s)^4 weight 0.970. Scratch: chopen/app11/he4D.py.
    he4_d_fid99: Tagged = Assumed((15, 17), "COMPUTED: D for fidelity 0.99 (beta = 10, 2)",
                                  "scratchpad chopen/app11/he4D_u1.0_s1.0_b{10.0,2.0}.json")
    he4_d_fid999: Tagged = Assumed((77, 103), "COMPUTED: D for fidelity 0.999 (beta = 2, 10)",
                                   "scratchpad chopen/app11/he4D_u1.0_s1.0_b{2.0,10.0}.json")
    he4_d_dipole_1pct: Tagged = Assumed(1, "COMPUTED: the single largest configuration already gives the dipole "
                                           "<O> within 1% (0.1% needs D = 11), both beta",
                                        "scratchpad chopen/app11/he4D_u1.0_s1.0_b{2.0,10.0}.json")
    qrom_passes: Tagged = Assumed((2, 3), "R11-R14 lookup-only accounting (load, uncompute, optional "
                                          "verification pass), kept as the record of the retired print; R15 "
                                          "prices one complete load per shot (loads_per_shot)",
                                  f"{TEX}:126")
    # R14 item 3: the QROM address is sized to the printed priced condition D <~ 2e3 (app11:124), not to
    # the cap crossing: ceil(log2 2000) = 11 address bits, 10 unary-iteration AND ancillae.

    # ---- R14 item 1a: the 4He dipole operator, lambda and <O> (derived here) ----------------------------
    # Dipole channel: coherent conversion through a virtual photon couples to the proton charge,
    # gamma_mu (1 + tau_3)/2 (Haxton et al. arXiv:2208.07945, mutoe_EFT_v25.tex:2800-2805), whose
    # lowest multipole is M_{00;p}(q) = sum_i (1 + tau3_i)/2 j0(q r_i) Y_00 (ibid. :2109-2114; KKO
    # hep-ph/0203110 :655 writes the same physics as the proton-density overlap D). Y_00 is a common
    # factor and cancels in lambda/<O>; it is dropped. q = m_mu (ibid. :500, 'q ~ m_mu').
    he4_charge_radius_fm: Tagged = Cited(1.67824, "Krauth_muonicHe_2021",
                                         "muonic-He Lamb shift, r_alpha = 1.67824(83) fm (Nature 589, 527)")
    proton_radius_fm: Tagged = Cited(0.8409, "CODATA2018", "r_p = 0.8409(4) fm")
    hbarc_mev_fm: Tagged = Cited(197.3269804, "CODATA2018", "hbar c")
    m_mu_mev: Tagged = Cited(105.6583755, "CODATA2018", "muon mass; q = m_mu for the dipole channel")
    m_nucleon_mev: Tagged = Cited(938.9187, "CODATA2018", "(m_p + m_n)/2, only for printing hbar omega")
    z_he4: Tagged = Stated(2, f"{TEX}:122", "4He: Z = 2 protons (A = 4)")
    mass_number_he4: Tagged = Stated(4, f"{TEX}:122", "4He")
    # ---- R14 item 1d: depth-capped MLAE on the 2028 Hadamard-test circuit ----------------------------
    # A = coherent oracle load + H + controlled BE query + H; Q = A S_0 A^dag S_chi. k Grover iterates
    # use 2k + 1 applications of A (Suzuki arXiv:1904.10246 Manuscript_v2.tex:247), each with its load.
    # A coherent load is the chapter's 'load + uncompute' pair of QROM passes (the address register must
    # be erased for the system to hold the superposition); the optional verification pass is dropped.
    # R15 ruling 2: the complete sparse load (Fomichev arXiv:2310.18410) erases its own index register, and
    # the Hadamard test ends in a measurement, so a shot makes ONE load: no uncompute pass. The optional
    # verification through the retained energy filter is a filter application, not a load, and is not
    # priced in the insertion-only benchmark. Under MLAE each application of A (or A^dag) carries one load.
    loads_per_shot: Tagged = Assumed(1, "derived here (R15): one complete load per shot; it erases its own "
                                        "index and the test measures, so nothing is uncomputed", f"{TEX}:126")
    # Referee M3/G2 (2026-10-04): the amplitude-angle rotations of the load's step 1 (Low-Kliuchnikov-Schaeffer
    # arXiv:1812.00954 LowTStatePrepQuantum.tex:321-339, 770-812): L levels, each a b-bit angle applied by a
    # controlled phase-gradient addition (CAdd = CNOT-conjugated Add, b - 1 Toffolis, Gidney arXiv:1709.06648),
    # plus a b-qubit Fourier (phase-gradient) state per shot. b from the source's bound 2 pi L / 2^b <= this.
    load_angle_error: Tagged = Assumed(1e-2, "angle-rounding error of the load, 2 pi L / 2^b (LKS :804); set "
                                             "equal to the per-shot synthesis budget", "arxiv_1709_06648")
    # Cramer-Rao bound (Suzuki Eq. Fisher_final, Manuscript_v2.tex:238): I(p) = sum_k N_k (2m_k+1)^2 /
    # (p(1-p)), reached by the ML estimate asymptotically (:255). With a = <O>/lambda = 2p - 1,
    # Var(a) = (1 - a^2) / sum_k N_k (2m_k+1)^2; all circuits at the deepest depth m gives
    # N = (1 - a^2)/(eps a m)^2 ~ (lambda/<O>)^2/(eps m)^2: constant 1, as Ch. 8 uses (app05:88).
    mlae_estimator_constant: Tagged = Cited(1.0, "arxiv_1904_10246",
                                            "Cramer-Rao bound with every circuit at m_max (Fisher_final); "
                                            "the LIS schedule k = 0..K is a sensitivity (mlae_lis_*)")
    # ---- R14 fix round: the 27Al sd dipole shot factor for the 2033 box (derived here, same construction
    # as he4_dipole). sd valence space: 0d (m_l = -2..2) + 1s = 6 spatial orbitals, 12 proton modes; 16O core
    # (0s^2 0p^6 protons) plus 5 valence protons. b from the 27Al charge radius by the same c.m.-corrected
    # HO fit as 4He: Z r_pt^2 = b^2 sum_occ (2n + l + 3/2) - Z (3/2) b^2 / A.
    al27_charge_radius_fm: Tagged = Cited(3.0610, "Angeli_Marinova_2013",
                                          "27Al rms charge radius 3.0610(31) fm (ADNDT 99, 69 (2013))")
    z_al27: Tagged = Stated(13, f"{TEX}:6", "27Al: Z = 13")
    mass_number_al27: Tagged = Stated(27, f"{TEX}:6", "27Al")
    z_core_al27: Tagged = Stated(8, f"{TEX}:44", "sd valence space: 16O core, 8 core protons (0s^2 0p^6)")
    # ---- R14 fix round: complete sparse-state load (sensitivity only; NEEDS_AUTHOR). Fomichev et al.
    # arXiv:2310.18410 initial.tex:350-351: amplitudes on the enumeration register, configuration lookup,
    # and index erasure, (2 log D - 2) D + 2^(log D + 1) + D < (2 log D + 3) D Toffolis, 5 log D - 3 ancillae.

    # ---- 2033 operator insertion: derived since R15 (ruling 4), gamma_0nubb / gamma_mu2e_two_body and
    # insertion_in_circuit. The pre-R15 Stated 5.4-5.6e7 (jj44, symmetry-restricted) and 4.4-4.6e8 (dense),
    # origin unknown, are kept as INSERTION_JJ44_PRE_R15 for the record.
    prep_lever: Tagged = Stated(10, f"{TEX}:163", "'improves on the tau_gs ~ 1/Delta projection cost by a factor of ten'")

    # ---- external comparison points (cited, not derived) -------------------
    fq_light_nucleus_t: Tagged = Cited(1e7, "arxiv_2507_22814", "'sits at the ~1e7-T scale' (app11:160)")
    ext_qubitization_sd_toffoli: Tagged = Cited(1e9, "arxiv_2607_21563", "'~1e9 Toffoli gates in the sd shell' (app11:144)")
    ext_qubitization_pb208_toffoli: Tagged = Cited(3e10, "arxiv_2607_21563", "'~3e10 ... 102-orbital 208Pb-core space'")
    # utility box (G5 open item 4, Claude's decision 2026-10-04): Mu2e-only base; no 0nubb dollar figure
    mu2e_tpc_musd: Tagged = Cited(315.7, "DOE_HEP_FY25_CJ", "'the $315.7M Mu2e total project cost' (app11:188); "
                                  "FY 2025 CJ, Science/HEP p. 245: BCP approved 2022-12-21, TPC $315,700,000")
    utility_fraction: Tagged = Stated(0.6, f"{TEX}:188", "'about 60% of the $315.7M Mu2e total project cost ... "
                                      "The fraction is not derived.' First set for 0nubb (NME sets the reach); for Mu2e the NMEs enter only the "
                                      "operator analysis, so 60% is a generous choice and the value is linear in it (box says so)")
    nme_spread: Tagged = Cited(3, "arxiv_2308_15634", "'factor-three spread across nuclear models' (app11:225)")

    # ---- r25 (rulings R1, R9, R10; H. Lamm 2026-10-02): T-depth, factories, and the two 2033 tiers ----
    # 2033 compile: G commuting Trotter rotations in flight at once (one Clifford diagonalizes a commuting
    # family; one parity ancilla, plus one for RUS, per concurrent rotation). factory.json (Ch. 3, 2033 run):
    # G = 1 is the gate-by-gate bound (F* = 1); G >= 10 feeds the 10-factory baseline; G_max ~ 10-60 at
    # N_anc = 50. G = 12 keeps F* >= 10 at both band ends with the insertion tail (G = 10 gives 9.9).
    trotter_concurrency: Tagged = Assumed(12, "commuting rotations in flight (factory.json Ch. 3 2033 run: "
                                              "G >= 10 validates the baseline; 12-24 of the 50-150 ancillas); "
                                              "app11:153 'With $12$ in flight'", "arxiv_1905_09749")
    # 2033 readout (shot audit, Ch. 3; ruling R2/R1 and the Ch. 3 rulings of 2026-10-02). The bare one-body
    # spin-independent (SI) elements in sd are pinned classically to +/-0.6% in any state (the derived
    # any-state band); the spin-dependent (SD, tensor) element carries the shots: <O_SD> = 2 <S_p> o_0d.
    sd_spin_expectation: Tagged = Assumed((0.30, 0.34), "<S_p>(27Al) in sd, shot-audit estimate; needs a "
                                                        "shell-model value (NEEDS_AUTHOR)", f"{TEX}:151")
    direct_readout_variance: Tagged = Assumed(1.0, "per-shot variance of the direct (basis-rotated Z) readout "
                                                   "of the SD one-body operator, variance ~1, estimate (shot audit); "
                                                   "Popoviciu worst case x8.3",
                                              f"{TEX}:151")
    projection_success: Tagged = Assumed(0.5, "heralded ground-state projection success p for 27Al (deformed; "
                                              "shot-audit estimate, band 0.1-0.9; computable exactly in sd)",
                                         f"{TEX}:151")
    projection_success_band: Tagged = Assumed((0.1, 0.9), "sensitivity band on p (shot audit)", f"{TEX}:151")
    first_result_scans: Tagged = Stated(1, f"{TEX}:153", "first result at one Hamiltonian (R1; shot audit)")
    eps_first_generic: Tagged = Assumed(0.3, "R1 report-wide first-result statistical error (H. Lamm, "
                                             "2026-10-02); quoted for comparison only", f"{TEX}:151")
    # ---- R7 (simplified 0nubb observables, simplobs.json Ch. 3, author 'item 7 sounds good') ----
    # O is a scalar, so each projection need resolve only the next 0+ state. Referee M3 (2026-10-05): the gaps
    # are those of the truncated Hamiltonian actually priced, by exact M-scheme diagonalization (Lanczos with a
    # J^2 shift, every state checked J = 0): KB3G (arXiv nucl-th/0012077; snt from the KSHELL distribution,
    # (42/A)^(1/3) scaling at A = 48) restricted to the kept orbits without renormalization. KB3G is priced: it
    # gives the larger sum 1/Delta_Ca + 1/Delta_Ti in both truncated spaces; GXPF1A (Honma et al., EPJA 25 s01,
    # 499) is the record: f7p3 6.0595 / 4.0233, f7p3p1 5.6371 / 3.8668, pf 5.2749 / 4.0477 (pf: 1/Delta sum
    # 0.4366 vs KB3G 0.4277, i.e. ~2% dearer; printed).
    # Scratch: chopen/app11/{sm.py, run_gaps.py, gap_*.json}.
    gap_0plus_f7p3_ca48_mev: Tagged = Assumed(6.100132, "COMPUTED: 48Ca 0+_2 - 0+_1, KB3G in f7/2 p3/2 (M-dim 57)",
                                              "scratchpad chopen/app11/gap_kb3g_f7p3_48Ca_J0.0.json")
    gap_0plus_f7p3_ti48_mev: Tagged = Assumed(3.684720, "COMPUTED: 48Ti 0+_2 - 0+_1, KB3G in f7/2 p3/2 (M-dim 5296)",
                                              "scratchpad chopen/app11/gap_kb3g_f7p3_48Ti_J0.0.json")
    gap_0plus_f7p3p1_ca48_mev: Tagged = Assumed(5.500798, "COMPUTED: 48Ca, KB3G in f7/2 p3/2 p1/2 (M-dim 325)",
                                                "scratchpad chopen/app11/gap_kb3g_f7p3p1_48Ca_J0.0.json")
    gap_0plus_f7p3p1_ti48_mev: Tagged = Assumed(3.603676, "COMPUTED: 48Ti, KB3G in f7/2 p3/2 p1/2 (M-dim 24453)",
                                                "scratchpad chopen/app11/gap_kb3g_f7p3p1_48Ti_J0.0.json")
    gap_0plus_pf_ca48_mev: Tagged = Assumed(5.174846, "COMPUTED: 48Ca, KB3G in full pf (M-dim 12022)",
                                            "scratchpad chopen/app11/gap_kb3g_pf_48Ca_J0.0.json")
    gap_0plus_pf_ti48_mev: Tagged = Assumed(4.264629, "COMPUTED: 48Ti, KB3G in full pf (M-dim 634744)",
                                            "scratchpad chopen/app11/gap_kb3g_pf_48Ti_J0.0.json")
    # 27Al: the 2033 deliverable keeps the generic band (its Hamiltonians are VS-IMSRG sd interactions, not
    # computed here); USDB (Brown & Richter, PRC 74, 034315; KSHELL usdb.snt, whose header is mislabelled USDA but
    # whose single-particle energies 2.1117 / -3.9257 / -3.2079 MeV are USDB's; (18/A)^0.3 at A = 27) is the check.
    # Exact Lanczos, M = 1/2 (dim 80115): 1/2+ 0.882, 3/2+ 1.063, 7/2+ 2.328, 5/2+ 2.699 MeV; M = 5/2 (dim 64299).
    gap_al27_usdb_mev: Tagged = Assumed(2.327866, "COMPUTED: first excited state reachable from an M = 5/2 "
                                                  "reference (J >= 5/2; the 7/2+), USDB in sd",
                                        "scratchpad chopen/app11/gap_usdb_sd_27Al_J-1.0.json")
    gap_al27_usdb_same_j_mev: Tagged = Assumed(2.698773, "COMPUTED: 5/2+_2 - 5/2+_1, USDB in sd (record)",
                                               "scratchpad chopen/app11/gap_usdb_sd_27Al_J2.5.json")
    gap_0plus_ensdf_ca48_mev: Tagged = Cited(4.283, "ENSDF", "48Ca first excited 0+ (4283 keV); record (was priced "
                                                             "until referee M3)")
    gap_0plus_ensdf_ti48_mev: Tagged = Cited(2.997, "ENSDF", "48Ti first excited 0+ (2997 keV); record")
    eps_simplified: Tagged = Assumed(0.3, "R1: 30% statistical error on |M| for the 0nubb first light",
                                     f"{TEX}:161")
    lam_over_m_f7p3: Tagged = Assumed(69.1, "lambda/|M| in f7/2 p3/2, schematic-interaction toy "
                                            "(simplobs adv/j0_big.out); not derived here", f"{TEX}:161")
    p_ref_f7p3: Tagged = Assumed(0.93, "closed-f7/2 reference overlap in f7/2 p3/2 (toy)", f"{TEX}:161")
    lam_over_m_f7p3p1: Tagged = Assumed(164.0, "lambda/|M| in f7/2 p3/2 p1/2 (toy); enters the qDRIFT count",
                                        f"{TEX}:161,163")
    lam_over_m_pf: Tagged = Assumed((600.0, 720.0), "lambda/|M| in full pf: Pauli-LCU lambda 208.4 over "
                                                    "|M| ~ 0.29-0.35 (simplobs)", f"{TEX}:161")
    rabi_drive_mev: Tagged = Assumed((0.3, 0.5), "theta|M|, the Rabi frequency of the resonant drive; kept "
                                                 "<= 0.5 MeV for a <~6% Stark bias (toy)", f"{TEX}:163")
    qdrift_eps: Tagged = Assumed(0.1, "qDRIFT error on the drive term; N_q = 2 (lambda/M)^2 (pi/2)^2 / eps_q",
                                 "Campbell_qDRIFT_2019")
    rabi_shots: Tagged = Assumed((250, 1000), "Rabi-readout shots incl. a 5-10-detuning resonance scan; 250 "
                                              "give ~5% on |M| at t_pi/2 (toy)", f"{TEX}:163")

    # fields the formulas use as scalars; a (lo, hi) here would silently misprice a cell
    _SCALAR_FIELDS = ("eps_synth_total", "qubits_per_orb", "dt_gev_inv", "gap_large_mev", "gap_small_mev",
                      "n_orb_he4", "n_orb_sd", "n_orb_pf", "n_orb_jj44", "n_orb_jj55",
                      "prep_lever", "gamma_trunc", "gamma_256", "gamma_512", "qrom_t_per_amplitude",
                      "be_queries_2028",
                      "select_and_deficit", "insertion_t_per_toffoli",
                      "d_amplitudes", "t_gate_s", "shot_overhead_s", "campaign_horizon_yr", "eps_mu2e", "eps_0nubb",
                      "n_op_mu2e", "n_op_0nubb", "n_pairs_cross_check", "n_survey_elements",
                      "utility_fraction", "mu2e_tpc_musd",
                      "rfi_t_2028", "rfi_t_2033", "rfi_lq_2033", "mass_number_ge", "fault_budget",
                      "he4_charge_radius_fm", "proton_radius_fm", "hbarc_mev_fm", "m_mu_mev",
                      "m_nucleon_mev", "z_he4", "mass_number_he4", "loads_per_shot",
                      "mlae_estimator_constant", "al27_charge_radius_fm", "z_al27", "mass_number_al27",
                      "z_core_al27", "trotter_concurrency", "direct_readout_variance", "projection_success",
                      "first_result_scans", "eps_first_generic", "gap_0plus_f7p3_ca48_mev", "gap_0plus_f7p3_ti48_mev",
                      "gap_0plus_f7p3p1_ca48_mev", "gap_0plus_f7p3p1_ti48_mev", "gap_0plus_pf_ca48_mev",
                      "gap_0plus_pf_ti48_mev", "gap_0plus_ensdf_ca48_mev", "gap_0plus_ensdf_ti48_mev",
                      "gap_al27_usdb_mev", "gap_al27_usdb_same_j_mev",
                      "eps_simplified", "lam_over_m_f7p3", "p_ref_f7p3", "lam_over_m_f7p3p1", "qdrift_eps",
                      "load_angle_error")

    def __post_init__(self):
        for f in fields(self):
            v = getattr(self, f.name)
            if not isinstance(v, Tagged):
                raise TypeError(f"{f.name} must be Tagged")
            if v.is_range and v.lo > v.hi:
                raise ValueError(f"{f.name}: lo > hi")
            if v.lo <= 0:
                raise ValueError(f"{f.name} must be positive")
        for name in self._SCALAR_FIELDS:
            if getattr(self, name).is_range:
                raise ValueError(f"{name} must be a single value, not a range")
        if not (self.gap_small_mev.value < self.gap_mev.lo <= self.gap_mev.hi < self.gap_large_mev.value):
            raise ValueError("gap bands must be ordered small < generic < large")
        if self.dt_gev_inv.value >= tau_gs(self.gap_large_mev.value):
            raise ValueError("Trotter step must be shorter than the shortest tau_gs")
        if self.gamma_trunc.value < 2:
            raise ValueError("the LCU needs at least two terms")
        if not (0 < self.fault_budget.value <= 1):
            raise ValueError("fault_budget is expected faults per shot, in (0, 1]")


# --------------------------------------------------------------------------- #
# Helpers: the chapter's formulas
# --------------------------------------------------------------------------- #

def tau_gs(gap_mev: float) -> float:
    """tau_gs ~ 1/Delta in GeV^-1 (app11:87): 1-2 MeV -> 1000-500 GeV^-1."""
    return 1e3 / gap_mev


def n_trotter(a: Assumptions, gap: tuple[float, float]) -> tuple[float, float]:
    """N_Trotter = c_proj tau_gs / dt over a gap band (app11:88): (pi, Delta_hi) .. (2pi, Delta_lo)."""
    lo = a.c_proj.lo * tau_gs(gap[1]) / a.dt_gev_inv.value
    hi = a.c_proj.hi * tau_gs(gap[0]) / a.dt_gev_inv.value
    return lo, hi


def c_T(a: Assumptions, n_rot: float) -> float:
    """R-TOL: T per rotation for a circuit of n_rot synthesized rotations per shot, the full RUS fit
    1.15 log2(1/eps_rot) + 9.2 at eps_rot = sqrt(eps_syn / n_rot) (randomized, incoherent)."""
    return t_per_rotation(eps_rot_for(n_rot, a.eps_synth_total.value))


def prep_rotations(a: Assumptions, n_proj: int, n_orb: int, gap: tuple[float, float],
                   lever: float = 1.0) -> tuple[float, float]:
    """N_rot per shot at each band end: N_proj N_Trotter N_orb^4 (dense Pauli strings), N_Trotter / lever."""
    lo, hi = n_trotter(a, gap)
    return (n_proj * lo / lever * n_orb ** 4, n_proj * hi / lever * n_orb ** 4)


def prep_c_T(a: Assumptions, n_proj: int, n_orb: int, gap: tuple[float, float],
             lever: float = 1.0) -> tuple[float, float]:
    """c_T at each band end, each end its own circuit (R-TOL)."""
    return tuple(c_T(a, n) for n in prep_rotations(a, n_proj, n_orb, gap, lever))


def prep_eps_rot(a: Assumptions, n_proj: int, n_orb: int, gap: tuple[float, float]) -> tuple[float, float]:
    """eps_rot at each band end (R-TOL); the low-T end has the looser tolerance."""
    return tuple(eps_rot_for(n, a.eps_synth_total.value) for n in prep_rotations(a, n_proj, n_orb, gap))


def prep_t(a: Assumptions, n_proj: int, n_orb: int, gap: tuple[float, float],
           lever: float = 1.0) -> tuple[float, float]:
    """eq:Ngate11: N_reg (c_proj tau_gs/dt / lever) N_orb^4 c_T(N_rot), for the coherent state prep alone.
    `lever` divides the projection depth (app11:163); the shorter circuit gets its own eps_rot."""
    n = prep_rotations(a, n_proj, n_orb, gap, lever)
    return (n[0] * c_T(a, n[0]), n[1] * c_T(a, n[1]))


def register(a: Assumptions, n_reg: int, n_orb: float) -> int:
    """System qubits, eq:Nq11 JW: n_reg . 4 N_orb. Since R14 (option A) every cell is one register:
    the 0nubb daughter is reached by projection in the parent's Fock register, so n_reg = 1."""
    n = n_reg * a.qubits_per_orb.value * n_orb
    if abs(n - round(n)) > 1e-9:
        raise ValueError(f"register {n} is not a whole number of modes")
    return int(round(n))


def n_orb_mo_ru(a: Assumptions) -> float:
    """Effective N_orb of the 100Mo/100Ru space (R15 ruling 3): protons in jj44 (2 N_orb_jj44 = 22 states),
    neutrons in jj55 (2 N_orb_jj55 = 32 states), 54 modes = 4 x 13.5."""
    modes = 2 * a.n_orb_jj44.value + 2 * a.n_orb_jj55.value
    return modes / a.qubits_per_orb.value


def ho_spatial_orbitals(emax: int) -> int:
    """Spatial (n, l, m_l) HO orbitals with 2n + l <= e_max: sum_N (N+1)(N+2)/2 = (e+1)(e+2)(e+3)/6."""
    return (emax + 1) * (emax + 2) * (emax + 3) // 6


def gap_l_band(a):
    return _band(a.gap_large_mev)


def _band(gap):
    return (gap.lo, gap.hi) if gap.is_range else (gap.value, gap.value)


def _inv(r):
    return (1.0 / r[1], 1.0 / r[0])


def eps_l(a: Assumptions, t: tuple[float, float]) -> tuple[float, float]:
    """Required per-logical-operation error over a T band: fault_budget / N (app11:93, R3)."""
    return (a.fault_budget.value / t[1], a.fault_budget.value / t[0])


def round_sig(x: float, digits: int = 2) -> float:
    """Round once at print, to `digits` significant figures (ruling R4).  127.3 -> 130, 0.9935 -> 0.99."""
    if x == 0:
        return 0.0
    e = math.floor(math.log10(abs(x)))
    return round(x, digits - 1 - e)


def print_range(r: tuple[float, float], digits: int = 2) -> str:
    """'lo--hi' with each end rounded once to `digits` figures; integers print without a decimal."""
    def one(v):
        v = round_sig(v, digits)
        return f"{int(v)}" if float(v).is_integer() and abs(v) >= 10 else f"{v:g}"
    return f"{one(r[0])}--{one(r[1])}"


def _mul(r, k):
    return (r[0] * k, r[1] * k)


def _div(r, k):
    return (r[0] / k, r[1] / k)


def synthesis_check(a: Assumptions) -> dict:
    """c_T across the survey under R-TOL: RUS at the per-rotation tolerance each cell's rotation count
    demands (incoherent, randomized; what the model prices) and the coherent alternative (app11:90).
    The survey ends are 27Al sd low (fewest rotations) and jj55 high (most; Mo/Ru before R15)."""
    gap = _band(a.gap_mev)
    nt = n_trotter(a, gap)
    n_rot_min = 1 * nt[0] * a.n_orb_sd.value ** 4                                  # 27Al sd, lo
    # R15: the largest survey cell is jj55 (N_orb 16), not Mo/Ru (13.5 since ruling 3; was 27)
    n_rot_max = 2 * nt[1] * max(a.n_orb_pf.value, a.n_orb_jj44.value, a.n_orb_jj55.value, n_orb_mo_ru(a)) ** 4
    et = a.eps_synth_total.value
    inc = tuple(t_per_rotation(eps_per_rotation(et, n, "incoherent")) for n in (n_rot_min, n_rot_max))
    coh = tuple(t_per_rotation(eps_per_rotation(et, n, "coherent")) for n in (n_rot_min, n_rot_max))
    return dict(rotations_survey=(n_rot_min, n_rot_max),
                eps_rot_survey=tuple(eps_rot_for(n, et) for n in (n_rot_min, n_rot_max)),
                c_T_rus_at_eps_1e4=t_per_rotation(1e-4),
                c_T_survey_range=inc,
                c_T_coherent_range=coh,
                c_T_coherent_over_randomized=(coh[0] / inc[0], coh[1] / inc[1]))


def insertion_cost(a: Assumptions, gamma: int, d: int, budget_rotations: float | None = None) -> dict:
    """Gate-level T-count of one shot of the 2028 operator insertion (derived here, R12 item 13, R13; app11:124).

    Hadamard test on a degree-d polynomial of a Gamma-term LCU block encoding of O_X. The printed
    entry is d = 1 (R13): the bare controlled BE, no phase steps. d >= 2 is the pre-R13 QSP reading,
    kept as a sensitivity:
      - JW maps each one-body term a_p^dag a_q + h.c. of the dipole operator to Pauli strings, so the
        term applied under the unary flag is a controlled Pauli string: Clifford, 0 T.
      - PREPARE: rotation tree over the Gamma real coefficient amplitudes on k = ceil(log2 Gamma) address
        qubits, 2^k - 1 R_y (level l is a uniformly controlled R_y with 2^l angles: 2^l R_y + 2^l CNOT);
        the coefficient signs go into SELECT as Paulis. PREPARE and PREPARE^dag per query, uncontrolled
        (with SELECT off, PREPARE^dag undoes PREPARE, so only SELECT carries the test-qubit control).
      - SELECT: unary iteration controlled on the test qubit, Gamma - 1 Toffolis (Babbush 2018), k AND
        ancillae, measurement-based uncompute.
      - QSP (d >= 2 only): a degree-d sequence has d queries and d + 1 phase steps e^{i phi (2 Pi_0 - 1)}
        on the address register (GSLW lemma:implementingPhasedSeq), each one flag = [address == 0]
        (k - 1 Toffolis, measurement uncompute) plus one synthesized R_z. With SELECT off the address
        stays in |0>, so the phase steps add a known global phase to that branch; it is removed in the
        readout, and they need no control. d = 1 with phi = 0 is U itself: no phase steps.
    Rotations are priced by R-TOL at N_rot = every rotation in the shot; Toffolis at 7 T (R5).
    `budget_rotations` (R14, MLAE) overrides the R-TOL count: a query inside a depth-capped amplitude-
    estimation circuit is synthesized to the deepest circuit's budget, N_rot = m_max x (rotations per query).
    """
    k = math.ceil(math.log2(gamma))
    rot_prepare = 2 ** k - 1
    n_rot_prep = 2 * d * rot_prepare
    n_phase = d + 1 if d >= 2 else 0
    n_rot_phase = n_phase
    n_rot = n_rot_prep + n_rot_phase
    tof_select = d * (gamma - a.select_and_deficit.value)
    tof_phase = n_phase * (k - 1)
    n_tof = tof_select + tof_phase
    eps = eps_rot_for(n_rot if budget_rotations is None else budget_rotations, a.eps_synth_total.value)
    ct = t_per_rotation(eps)
    t_tof = a.insertion_t_per_toffoli.value
    return dict(k=k, rot_per_prepare=rot_prepare, n_rot_prepare=n_rot_prep, n_rot_phase=n_rot_phase,
                n_rot=n_rot, toffoli_select=tof_select, toffoli_phase=tof_phase, n_toffoli=n_tof,
                eps_rot=eps, c_T=ct, t_rot=n_rot * ct, t_toffoli=n_tof * t_tof,
                t=n_rot * ct + n_tof * t_tof,
                select_and_ancillae=k, address_qubits=k)


# --------------------------------------------------------------------------- #
# R14 item 1a: the 4He dipole operator on the chapter's N_orb = 10 active space
# --------------------------------------------------------------------------- #

# The Ch. 2 active space (app01:125, '0s-0p-1s0d basis on A = 4'), spherical HO orbitals (n, l, m_l):
# 0s (1) + 0p (3) + 1s (1) + 0d (5) = 10 spatial orbitals, x 2 spin x 2 isospin = 40 JW modes.
HE4_SPATIAL = ((0, 0, 0),) + tuple((0, 1, m) for m in (-1, 0, 1)) + ((1, 0, 0),) + \
    tuple((0, 2, m) for m in range(-2, 3))


def he4_oscillator_length(a: Assumptions) -> dict:
    """b from the 4He charge radius: the 0s^4 point-proton radius, c.m.-corrected, (A-1)/A (3/2) b^2,
    set equal to r_pt^2 = r_ch^2 - r_p^2. Returns b [fm], hbar omega [MeV], q = m_mu [fm^-1], y = (qb/2)^2."""
    r_pt2 = a.he4_charge_radius_fm.value ** 2 - a.proton_radius_fm.value ** 2
    A = a.mass_number_he4.value
    b = math.sqrt(r_pt2 / (1.5 * (A - 1) / A))
    hw = a.hbarc_mev_fm.value ** 2 / (a.m_nucleon_mev.value * b * b)
    q = a.m_mu_mev.value / a.hbarc_mev_fm.value
    return dict(r_pt=math.sqrt(r_pt2), b=b, hbar_omega=hw, q=q, qb=q * b, y=(q * b / 2) ** 2)


def ho_radial(n: int, l: int, r, b: float):
    """Normalized HO radial function R_nl(r), int R^2 r^2 dr = 1, positive at the origin."""
    from scipy import special
    import numpy as np
    x = np.asarray(r) / b
    norm = math.sqrt(2 * math.factorial(n) / (b ** 3 * math.gamma(n + l + 1.5)))
    return norm * x ** l * np.exp(-x * x / 2) * special.eval_genlaguerre(n, l + 0.5, x * x)


def j0_radial_me(n1: int, n2: int, l: int, q: float, b: float) -> float:
    """<n1 l| j0(q r) |n2 l> by quadrature (the L = 0 multipole; only same-l orbitals connect)."""
    from scipy import integrate, special
    f = lambda r: float(ho_radial(n1, l, r, b) * ho_radial(n2, l, r, b)) * special.spherical_jn(0, q * r) * r * r
    return integrate.quad(f, 0.0, 25.0 * b, limit=400, epsabs=1e-13, epsrel=1e-12)[0]


def j0_radial_me_analytic(n1: int, n2: int, l: int, y: float) -> float | None:
    """Closed forms for the orbitals used here, y = (qb/2)^2 (checks the quadrature)."""
    e = math.exp(-y)
    table = {(0, 0, 0): e,                                        # 0s
             (0, 0, 1): (1 - 2 * y / 3) * e,                      # 0p
             (0, 0, 2): (1 - 4 * y / 3 + 4 * y * y / 15) * e,     # 0d
             (1, 1, 0): (1 - 4 * y / 3 + 2 * y * y / 3) * e,      # 1s
             (0, 1, 0): math.sqrt(2 / 3) * y * e}                 # 0s-1s (R_nl > 0 at the origin)
    key = (min(n1, n2), max(n1, n2), l)
    return table.get(key)


def he4_dipole_one_body(a: Assumptions):
    """The 40 x 40 one-body matrix o_pq of O = sum_i (1 + tau3_i)/2 j0(q r_i) on the JW modes.
    Mode order: p = 4 * orbital + 2 * isospin + spin, isospin 0 = proton. Returns (o, mode labels, osc)."""
    import numpy as np
    osc = he4_oscillator_length(a)
    n_sp = len(HE4_SPATIAL)
    o = np.zeros((4 * n_sp, 4 * n_sp))
    labels = []
    for i, (ni, li, mi) in enumerate(HE4_SPATIAL):
        for t in (0, 1):
            for s in (0, 1):
                labels.append((ni, li, mi, "p" if t == 0 else "n", "up" if s == 0 else "dn"))
    for i, (ni, li, mi) in enumerate(HE4_SPATIAL):
        for j, (nj, lj, mj) in enumerate(HE4_SPATIAL):
            if li != lj or mi != mj:
                continue                                   # j0 is a scalar: diagonal in l, m_l
            me = j0_radial_me(ni, nj, li, osc["q"], osc["b"])
            for s in (0, 1):                               # spin-independent, protons only (t = 0)
                o[4 * i + s, 4 * j + s] = me
    return o, tuple(labels), osc


def jw_one_body_paulis(o, tol: float = 1e-12) -> dict:
    """Jordan-Wigner Pauli decomposition of sum_pq o_pq a_p^dag a_q, o real symmetric.
    a_p^dag a_p = (1 - Z_p)/2; a_p^dag a_q + a_q^dag a_p = (X_p Z..Z X_q + Y_p Z..Z Y_q)/2 (p < q).
    Keys are tuples of (qubit, 'X'|'Y'|'Z'); () is the identity."""
    n = len(o)
    terms: dict = {}

    def add(key, c):
        terms[key] = terms.get(key, 0.0) + c
    for p in range(n):
        if abs(o[p][p]) > tol:
            add((), float(o[p][p]) / 2)
            add(((p, "Z"),), -float(o[p][p]) / 2)
    for p in range(n):
        for q in range(p + 1, n):
            if abs(o[p][q]) > tol:
                zs = tuple((k, "Z") for k in range(p + 1, q))
                add(((p, "X"),) + zs + ((q, "X"),), float(o[p][q]) / 2)
                add(((p, "Y"),) + zs + ((q, "Y"),), float(o[p][q]) / 2)
    return {k: v for k, v in terms.items() if abs(v) > tol}


def he4_dipole(a: Assumptions) -> dict:
    """lambda, <O>, Gamma and the Hadamard-test shot factor (lambda/|<O>|)^2 for the 2028 dipole benchmark.

    <O> in the closed-shell 0s^4 reference is sum over occupied modes of o_pp = Z <0s|j0|0s>. For ANY
    4He state in this space <O> = Tr(rho_p o_p), rho_p the proton one-body density matrix (eigenvalues in
    [0, 1], trace Z), so <O> lies between Z x the smallest and Z x the largest eigenvalue of the proton
    block of o (each is spin-doubled). That band is reported as the state-independence check.
    R15 ruling 1: the LCU is O - c_0 I (identity subtracted classically); `shot_factor` is
    (lambda'/|<O>|)^2 = 19.85, band 17.9-32.8; the R14 value with the identity kept (77.5) is
    `shot_factor_with_identity`.
    """
    import numpy as np
    o, labels, osc = he4_dipole_one_body(a)
    paulis = {k: float(v) for k, v in jw_one_body_paulis(o).items()}
    lam = float(sum(abs(c) for c in paulis.values()))
    c_id = float(paulis.get((), 0.0))
    z = a.z_he4.value
    occ = [p for p, lab in enumerate(labels) if lab[:3] == (0, 0, 0) and lab[3] == "p"]   # 0s p up/dn
    assert len(occ) == z
    exp_cs = float(sum(o[p, p] for p in occ))
    prot = [p for p, lab in enumerate(labels) if lab[3] == "p"]
    ev = np.linalg.eigvalsh(o[np.ix_(prot, prot)])
    exp_band = (float(z * ev.min()), float(z * ev.max()))      # spin-degenerate: two lowest = 2 x min
    # R14 record: the LCU with the identity string kept, (lambda/|<O>|)^2 = 77.5
    ratio_with_id = lam / abs(exp_cs)
    ratio_band_with_id = (lam / exp_band[1], lam / exp_band[0])
    # R15 ruling 1: the identity coefficient c_0 = Tr(o)/2 is a classically known constant (<I> = 1 in any
    # state). Block-encode O - c_0 I and add c_0 back: lambda' = lambda - |c_0|, Gamma - 1 strings. The
    # Hadamard test estimates (<O> - c_0)/lambda' with per-shot variance <= 1; the error on <O> is lambda'
    # times the statistical error, so relative accuracy eps on <O> needs N = (lambda'/(eps |<O>|))^2.
    # (lambda'/|<O> - c_0|)^2 = 1.77 is NOT the factor: it would set relative accuracy on <O> - c_0.
    lam_sub = lam - abs(c_id)
    exp_sub = exp_cs - c_id
    ratio = lam_sub / abs(exp_cs)
    ratio_band = (lam_sub / exp_band[1], lam_sub / exp_band[0])
    a_shift = exp_sub / lam_sub                                  # the test's mean, -0.751
    a_shift_band = ((exp_band[0] - c_id) / lam_sub, (exp_band[1] - c_id) / lam_sub)
    radial = {}
    for key in ((0, 0, 0), (0, 1, 0), (1, 1, 0), (0, 0, 1), (0, 0, 2)):
        radial[key] = (j0_radial_me(key[0], key[1], key[2], osc["q"], osc["b"]),
                       j0_radial_me_analytic(key[0], key[1], key[2], osc["y"]))
    counts = dict(identity=int(() in paulis),
                  z=sum(1 for k in paulis if len(k) == 1),
                  xx_yy=sum(1 for k in paulis if len(k) >= 2))
    return dict(osc=osc, o=o, labels=labels, paulis=paulis, gamma=len(paulis), gamma_counts=counts,
                gamma_without_identity=len(paulis) - counts["identity"],
                lam=lam, identity_coefficient=c_id, exp_closed_shell=exp_cs,
                exp_band_any_state=exp_band, proton_block_eigenvalues=tuple(float(x) for x in ev),
                ratio=ratio, shot_factor=ratio ** 2, ratio_band=ratio_band,
                shot_factor_band=(ratio_band[0] ** 2, ratio_band[1] ** 2),
                gamma_lcu=len(paulis) - counts["identity"], lam_lcu=lam_sub,
                test_mean=a_shift, test_mean_band=a_shift_band, bernoulli_variance=1 - a_shift ** 2,
                ratio_with_identity=ratio_with_id, shot_factor_with_identity=ratio_with_id ** 2,
                shot_factor_band_with_identity=(ratio_band_with_id[0] ** 2, ratio_band_with_id[1] ** 2),
                radial_me=radial, trace_o=float(np.trace(o)),
                lam_identity_subtracted=lam_sub, exp_identity_subtracted=exp_sub,
                shot_factor_over_shifted_mean=(lam_sub / abs(exp_sub)) ** 2)


# --------------------------------------------------------------------------- #
# R14 item 1d: depth-capped MLAE on the 2028 Hadamard-test circuit
# --------------------------------------------------------------------------- #

def mlae_2028(a: Assumptions, gamma: int, d_amp: float, t_load_override: float | None = None,
              load_angles: bool = False) -> dict:
    """Largest depth m = 2k + 1 (applications of A) with circuit T <= the 2028 cap, rotations synthesized
    to the deepest circuit's budget (R-TOL), solved self-consistently as Ch. 8 does (app05, R13).

    R15 ruling 2: A = one complete sparse load (Fomichev 2310.18410, 7 x sos_load_toffoli(D) T; A^dag
    carries its inverse at the same count) + one controlled BE query (insertion_cost at `gamma`,
    synthesized to m x its own rotations). Q = A S_0 A^dag S_chi: S_chi is Z on the test qubit (Clifford);
    S_0 reflects about |0> on every qubit A acts on and does not return to |0> for arbitrary input:
    system + the load's enumeration register (L = ceil(log2 D) at the printed D <~ 600: 10) + LCU address +
    test. The load's identification register is CNOT-computed from the system and CNOT-uncomputed with the
    system unchanged in between, so it returns to |0> for any input; the multi-control and unary-iteration
    workspaces are computed and uncomputed inside each step. C^{N-1}Z = N - 2 temporary ANDs at 7 T.
    T(m) = m T_A + (m-1)/2 T_S0. `t_load_override` replaces the complete load (lookup-only record).
    `load_angles` (referee verifier, 2026-10-04): adds the load's amplitude-angle rotations (load_angle_cost):
    the L (b - 1)-Toffoli phase-gradient additions in every application of A and of A^dag, and the
    phase-gradient (Fourier) state once per circuit, since it is catalytic and reused (b - 3 rotations
    synthesized with the circuit's m x 254, plus one T).
    """
    cap = a.rfi_t_2028.value
    system = register(a, 1, a.n_orb_he4.value)
    k_addr = math.ceil(math.log2(gamma))
    load_addr = math.ceil(math.log2(a.d_amplitudes.value))                  # 10 (ruling 3 rule, D <~ 600)
    n_s0 = system + load_addr + k_addr + a.test_qubits.value                # 58 at Gamma = 128
    t_s0 = (n_s0 - 2) * a.insertion_t_per_toffoli.value                     # 56 x 7 = 392
    t_load = (a.insertion_t_per_toffoli.value * sos_load_toffoli(d_amp)
              if t_load_override is None else t_load_override)
    rot_per_query = insertion_cost(a, gamma, 1)["n_rot"]
    ang = load_angle_cost(a, d_amp) if load_angles else None
    n_rot_once = ang["n_rot"] if load_angles else 0
    t_add = a.insertion_t_per_toffoli.value * ang["toffoli"] if load_angles else 0.0

    def t_app(m_budget):
        return (t_load + t_add + insertion_cost(a, gamma, 1,
                budget_rotations=m_budget * rot_per_query + n_rot_once)["t"])

    def t_once(m_budget):                                                   # phase-gradient state, once
        if not load_angles:
            return 0.0
        ct = t_per_rotation(eps_rot_for(n_rot_once + m_budget * rot_per_query, a.eps_synth_total.value))
        return n_rot_once * ct + 1

    def t_circ(m, m_budget):
        return m * t_app(m_budget) + (m - 1) // 2 * t_s0 + t_once(m_budget)

    def largest_odd_under(m_budget):
        m = 1
        while t_circ(m + 2, m_budget) <= cap:
            m += 2
        return m

    if t_circ(1, 1) > cap:
        raise ValueError("a single application of A exceeds the 2028 cap")
    # the iteration m -> largest m' fitting at budget m (record); it can 2-cycle at a depth boundary (R15:
    # m = 3 fits at the m = 1 budget but not at its own), so m_max is fixed directly as the largest odd m
    # whose circuit fits at its own budget (T(m, m) is increasing in m)
    iterations, m_cur, seen = [], 1, set()
    while m_cur not in seen:
        seen.add(m_cur)
        m_next = largest_odd_under(m_cur)
        iterations.append((m_cur, t_app(m_cur), m_next))
        m_cur = m_next
    m_max = 1
    while t_circ(m_max + 2, m_max + 2) <= cap:
        m_max += 2
    t_deep = t_circ(m_max, m_max)
    t_next = t_circ(m_max + 2, m_max + 2)                                 # the next depth, own budget
    q = insertion_cost(a, gamma, 1, budget_rotations=m_max * rot_per_query)
    return dict(m_max=m_max, grover_iterates=(m_max - 1) // 2, iterations=tuple(iterations),
                t_deepest=t_deep, t_next_depth=t_next, fraction_deepest=t_deep / cap,
                t_application=t_app(m_max), t_load_per_application=t_load, t_query=q["t"],
                eps_rot=q["eps_rot"], c_T=q["c_T"], n_rot_budget=m_max * rot_per_query,
                reflection_qubits=n_s0, reflection_ands=n_s0 - 2, t_reflection=t_s0, reflection_load_address=load_addr,
                peak_lq_reflection=n_s0 + (n_s0 - 2), gamma=gamma, d_amp=d_amp)


def sos_load_toffoli(d_amp: float) -> int:
    """Complete sparse-state load of D configurations, Fomichev et al. arXiv:2310.18410 initial.tex:350-351:
    (2L - 2) D + 2^(L+1) + D Toffolis, L = log D the enumeration-register width (ceil(log2 D) here, D not a
    power of two). Terms (initial.tex:304-337): 2^(L+1) the QROM state preparation of sum_i alpha_i |i> on
    the enumeration register (step 1); D the configuration lookup |i> -> |nu_i>|i> (step 2); (2L - 2) D the
    index erasure, D multi-controlled flips of |i> -> |0> conditioned on the (2L-1)-bit identifier b_i
    (step 5). The CNOT steps 3-4, 6 are Clifford. The source calls this the dominant Toffoli count."""
    d = int(round(d_amp))
    L = math.ceil(math.log2(d))
    return (2 * L - 2) * d + 2 ** (L + 1) + d


def sos_load_ancillae(d_amp: float) -> int:
    """5 log D - 3 ancillae for the same construction (ibid. :355): L enumeration + (2L - 1) identification
    + (2L - 2) workspace for the (2L-1)-controlled flips (the lookup's L - 1 AND ancillae fit inside it)."""
    return 5 * math.ceil(math.log2(int(round(d_amp)))) - 3


def load_angle_cost(a: Assumptions, d_amp: float, n_rot_other: int = 0) -> dict:
    """Referee M3/G2 (derived here): the amplitude-angle rotations that the 2^(L+1) step of the Fomichev load
    leaves unpriced. Step 1 is the QROM state preparation of Low, Kliuchnikov and Schaeffer (arXiv:1812.00954):
    at each of L levels a lookup writes a b-bit angle (its Toffolis are the source's 2^(L+1)) and the rotation
    is applied by adding that register into a b-qubit Fourier state, controlled by the target qubit (CAdd,
    LowTStatePrepQuantum.tex:770-792). CAdd is a b-bit adder conjugated by CNOTs from the target, b - 1
    Toffolis (Gidney 1709.06648). b is the smallest width with 2 pi L / 2^b <= load_angle_error (:804).
    The Fourier state costs b single-qubit Z rotations by pi/2^j, j = 0..b-1: Z and S are Clifford, one T,
    b - 3 synthesized at R-TOL with the shot's other rotations (`n_rot_other`). Signs of a real CI vector are
    a Clifford Z on a lookup bit. Register: step 1 needs b + b + (b - 1) = 3b - 1 extra qubits (angle, Fourier
    state, adder carries); the 19 + 18 identification and workspace qubits of the load are idle then and hold
    37 of them (folded into the 2028 peak since the M3 ruling, 2026-10-04). T-depth: the L additions are
    sequential, b - 1 AND layers each (ripple carry); the b - 3 Fourier-state rotations act on distinct qubits,
    one c_T layer (the T rides in it)."""
    d = int(round(d_amp))
    L = math.ceil(math.log2(d))
    b = math.ceil(math.log2(2 * PI * L / a.load_angle_error.value))
    tof = L * (b - 1)
    n_rot = b - 3
    eps = eps_rot_for(n_rot + n_rot_other, a.eps_synth_total.value)
    ct = t_per_rotation(eps)
    t = tof * a.insertion_t_per_toffoli.value + n_rot * ct + 1
    return dict(L=L, b=b, toffoli=tof, n_rot=n_rot, c_T=ct, eps_rot=eps, t=t,
                qubits_step1=2 * b + b - 1, idle_load_qubits=(2 * L - 1) + (2 * L - 2),
                t_depth_adds=L * (b - 1), t_depth_fourier=ct)


def al27_dipole_sd(a: Assumptions) -> dict:
    """The 2033 box's shot factor (lambda/|<O>|)^2 for the 27Al dipole in the sd valence space (derived here).

    Same operator and JW convention as he4_dipole. In the sd space (0d, 1s) j0 connects only equal l, m_l,
    so o is diagonal on the 12 valence proton modes: the operator is the identity plus 12 single-Z strings,
    lambda = Tr(o_p). R15 ruling 1: the identity (Tr(o_p)/2) is subtracted classically, lambda' = Tr(o_p)/2,
    and the factor is (lambda'/|<O>|)^2 at the upper end of the any-state band (1.45). The 16O core is a classical constant and is not credited (conservative). Any state
    with Z_val = 5 valence protons has <O> between the sum of the 5 smallest and the 5 largest mode values.
    b: Z r_pt^2 = b^2 sum_occ (2n + l + 3/2) - Z (3/2) b^2 / A, r_pt^2 = r_ch^2 - r_p^2 (the 4He fit, general A).
    """
    import numpy as np
    Z, A = a.z_al27.value, a.mass_number_al27.value
    z_val = Z - a.z_core_al27.value                                             # 5
    s_occ = 2 * 1.5 + 6 * 2.5 + z_val * 3.5                                     # 0s^2 0p^6 (sd)^5
    r_pt2 = a.al27_charge_radius_fm.value ** 2 - a.proton_radius_fm.value ** 2
    b = math.sqrt(r_pt2 / (s_occ / Z - 1.5 / A))
    q = a.m_mu_mev.value / a.hbarc_mev_fm.value
    hw = a.hbarc_mev_fm.value ** 2 / (a.m_nucleon_mev.value * b * b)
    o_0d = j0_radial_me(0, 0, 2, q, b)
    o_1s = j0_radial_me(1, 1, 0, q, b)
    modes = np.array([o_0d] * 10 + [o_1s] * 2)                                  # 12 valence proton modes
    paulis = jw_one_body_paulis(np.diag(modes))
    lam = float(sum(abs(c) for c in paulis.values()))
    c_id = float(paulis.get((), 0.0))
    lam_sub = lam - abs(c_id)                                                   # R15 ruling 1
    srt = np.sort(modes)
    exp_band = (float(srt[:z_val].sum()), float(srt[-z_val:].sum()))
    fac_band = ((lam_sub / exp_band[1]) ** 2, (lam_sub / exp_band[0]) ** 2)
    fac_band_with_id = ((lam / exp_band[1]) ** 2, (lam / exp_band[0]) ** 2)    # R14 record, 5.7-5.8
    return dict(b=b, hbar_omega=hw, q=q, qb=q * b, y=(q * b / 2) ** 2, o_0d=o_0d, o_1s=o_1s,
                gamma=len(paulis), gamma_lcu=len(paulis) - int(() in paulis), lam=lam,
                identity_coefficient=c_id, lam_lcu=lam_sub, z_valence=z_val, exp_band_any_state=exp_band,
                shot_factor_band=fac_band, shot_factor=fac_band[1],
                shot_factor_band_with_identity=fac_band_with_id)


def mlae_lis_constant(m_max: int) -> float:
    """Circuits relative to all-at-m_max for Suzuki's LIS, k = 0..K with K = (m_max - 1)/2 and equal
    N_shot per depth: (K + 1) m_max^2 / sum_k (2k + 1)^2 (Manuscript_v2.tex:271, I = N(2K+3)(2K+1)(K+1)/3)."""
    K = (m_max - 1) // 2
    return (K + 1) * m_max ** 2 / sum((2 * k + 1) ** 2 for k in range(K + 1))


# --------------------------------------------------------------------------- #
# R15 ruling 4: the 2033 operator insertions, one controlled BE query, Gamma counted here
# --------------------------------------------------------------------------- #

# Valence orbits (n, l, 2j) per species. sd, pf and jj44 are the chapter's spaces (app11:55, 135); jj55 is
# the 50-82 shell (0g7/2 1d5/2 1d3/2 2s1/2 0h11/2). States per species: sd 12, pf 20, jj44 22, jj55 32.
SHELLS = {
    "sd":   ((0, 2, 5), (0, 2, 3), (1, 0, 1)),
    "pf":   ((0, 3, 7), (0, 3, 5), (1, 1, 3), (1, 1, 1)),
    "jj44": ((0, 3, 5), (1, 1, 3), (1, 1, 1), (0, 4, 9)),
    "jj55": ((0, 4, 7), (1, 2, 5), (1, 2, 3), (2, 0, 1), (0, 5, 11)),
    # R7 truncation ladder for 48Ca -> 48Ti (simplobs): f7/2 p3/2 (24 modes), f7/2 p3/2 p1/2 (28 modes)
    "f7p3": ((0, 3, 7), (1, 1, 3)),
    "f7p3p1": ((0, 3, 7), (1, 1, 3), (1, 1, 1)),
}


def shell_states(shell: str) -> tuple:
    """M-scheme single-particle states (orbit index, 2j, 2m) of one species in `shell`."""
    return tuple((o, tj, m2) for o, (_, _, tj) in enumerate(SHELLS[shell]) for m2 in range(-tj, tj + 1, 2))


@_lru_cache(maxsize=None)
def _cg(j1, m1, j2, m2, J, M):
    """Clebsch-Gordan <j1 m1 j2 m2|J M>, all arguments doubled."""
    from sympy import Rational
    from sympy.physics.wigner import clebsch_gordan
    return float(clebsch_gordan(Rational(j1, 2), Rational(j2, 2), Rational(J, 2),
                                Rational(m1, 2), Rational(m2, 2), Rational(M, 2)))


def pair_quantum_numbers(orbits, s1, s2) -> tuple:
    """(2M, parity, allowed 2J) of the antisymmetrized pair of distinct states s1, s2 of one species:
    J in the triangle of j1, j2 with |M| <= J and a nonzero Clebsch-Gordan coefficient; J even when both
    sit in the same orbit (two identical nucleons in one j-shell)."""
    (o1, j1, m1), (o2, j2, m2) = s1, s2
    M = m1 + m2
    allowed = set()
    for J in range(abs(j1 - j2), j1 + j2 + 1, 2):
        if abs(M) > J or (o1 == o2 and (J // 2) % 2 == 1):
            continue
        if abs(_cg(j1, m1, j2, m2, J, M)) > 1e-12:
            allowed.add(J)
    return M, (orbits[o1][1] + orbits[o2][1]) % 2, frozenset(allowed)


def _ladder(p: int, dagger: bool):
    """JW a_p (a_p^dag) as Paulis X^x Z^z with coefficients: Z_<p (X_p +- i Y_p)/2, i Y = -X Z."""
    zs = (1 << p) - 1
    return ((0.5, 1 << p, zs), ((0.5 if dagger else -0.5), 1 << p, zs | (1 << p)))


def _pmul(A, B):
    out = {}
    for ca, xa, za in A:
        for cb, xb, zb in B:
            s = -1.0 if bin(za & xb).count("1") % 2 else 1.0     # Z^za X^xb = (-1)^{|za & xb|} X^xb Z^za
            k = (xa ^ xb, za ^ zb)
            out[k] = out.get(k, 0.0) + ca * cb * s
    return tuple((c, x, z) for (x, z), c in out.items() if abs(c) > 1e-14)


def jw_quartic(p: int, q: int, r: int, s: int):
    """a_p^dag a_q^dag a_r a_s under JW, as ((coef, xmask, zmask), ...)."""
    A = _ladder(p, True)
    for o, d in ((q, True), (r, False), (s, False)):
        A = _pmul(A, _ladder(o, d))
    return A


@_lru_cache(maxsize=None)
def gamma_0nubb(shell_p: str, shell_n: str, selection: bool = True) -> dict:
    """Pauli-string count of the 0nubb operator sum c_pqrs a_p^dag a_q^dag a_r a_s (p < q proton, r < s
    neutron) in one Fock register, protons in `shell_p`, neutrons in `shell_n` (derived here, R15).

    A rank-0 two-body operator (the GT, Fermi, tensor and contact pieces are all scalars) connects a proton
    pair and a neutron pair only if they share 2M, parity and some pair J allowed for both. These rules are
    exact; c_pqrs is then generically nonzero (accidental zeros of the radial integrals would only lower the
    count), so the selection-rule count is an upper bound on the nonzero c_pqrs. `selection=False` counts all
    C(n_p,2) C(n_n,2) terms (dense). Each term is 16 JW Pauli strings, and since p, q, r, s are four distinct
    modes the support {p, q, r, s} identifies the term: no two terms share a string, so Gamma = 16 N_terms.
    The selection-rule case is checked by explicit JW algebra."""
    P, N = shell_states(shell_p), shell_states(shell_n)
    pp = list(itertools.combinations(range(len(P)), 2))
    nn = list(itertools.combinations(range(len(N)), 2))
    if not selection:
        n_terms = len(pp) * len(nn)
        return dict(n_terms=n_terms, gamma=16 * n_terms, modes=(len(P), len(N)), checked=False)
    pq = [pair_quantum_numbers(SHELLS[shell_p], P[i], P[j]) for i, j in pp]
    nq = [pair_quantum_numbers(SHELLS[shell_n], N[i], N[j]) for i, j in nn]
    keys, n_terms = set(), 0
    for (i, j), A in zip(pp, pq):
        for (k, l), B in zip(nn, nq):
            if A[0] == B[0] and A[1] == B[1] and (A[2] & B[2]):
                n_terms += 1
                for _, x, z in jw_quartic(i, j, len(P) + k, len(P) + l):
                    keys.add((x, z))
    return dict(n_terms=n_terms, gamma=len(keys), modes=(len(P), len(N)), checked=True,
                dense_terms=len(pp) * len(nn))


@_lru_cache(maxsize=None)
def gamma_mu2e_two_body(shell: str, selection: bool = True, seed: int = 1) -> dict:
    """Pauli-string count of a general Hermitian, charge-conserving, rank-0 one- plus two-body operator on
    the protons and neutrons of `shell` (derived here, R15): the conservative cover for the six Mu2e
    categories with their two-body currents. Pair-pair terms need equal charge content (pp, nn, pn) and,
    with `selection`, equal 2M, parity and a common allowed J; the one-body part of a scalar connects equal
    (n, l, j, m) only, i.e. it is diagonal in these spaces. Strings are collected by explicit JW algebra
    with generic (random, fixed-seed) coefficients so that only exact cancellations remove a string; the
    identity is dropped (a classically known constant, R15 ruling 1)."""
    import random
    rnd = random.Random(seed)
    orbits = SHELLS[shell] + SHELLS[shell]
    nsp = len(SHELLS[shell])
    S = [(o, tj, m2, "p") for (o, tj, m2) in shell_states(shell)] + \
        [(o + nsp, tj, m2, "n") for (o, tj, m2) in shell_states(shell)]
    st = [s[:3] for s in S]
    pairs = list(itertools.combinations(range(len(S)), 2))
    info = {pr: pair_quantum_numbers(orbits, st[pr[0]], st[pr[1]]) + ("".join(sorted(S[pr[0]][3] + S[pr[1]][3])),)
            for pr in pairs}
    acc: dict = {}

    def add(terms, c):
        for cc, x, z in terms:
            acc[(x, z)] = acc.get((x, z), 0.0) + c * cc
    n_terms = 0
    for ia in range(len(pairs)):
        for ib in range(ia, len(pairs)):
            A, B = info[pairs[ia]], info[pairs[ib]]
            if A[3] != B[3]:
                continue
            if selection and not (A[0] == B[0] and A[1] == B[1] and (A[2] & B[2])):
                continue
            n_terms += 1
            (i, j), (k, l) = pairs[ia], pairs[ib]
            c = rnd.uniform(0.5, 1.5)
            add(jw_quartic(i, j, l, k), c)
            if ia != ib:
                add(jw_quartic(k, l, j, i), c)                      # + h.c.
    for p in range(len(S)):
        add(_pmul(_ladder(p, True), _ladder(p, False)), rnd.uniform(0.5, 1.5))
    nz = {k for k, v in acc.items() if abs(v) > 1e-12}
    return dict(n_terms=n_terms, gamma=len(nz - {(0, 0)}), identity=(0, 0) in nz, modes=len(S))


def insertion_in_circuit(a: Assumptions, gamma: int, n_rot_prep: float) -> dict:
    """One controlled BE query of a Gamma-string LCU inside a circuit that already holds n_rot_prep
    rotations: the R-TOL tolerance is set by the whole circuit's rotation count (R15 ruling 4)."""
    k = math.ceil(math.log2(gamma))
    n_rot_ins = 2 * (2 ** k - 1)
    q = insertion_cost(a, gamma, 1, budget_rotations=n_rot_prep + n_rot_ins)
    return dict(q, n_rot_circuit=n_rot_prep + n_rot_ins)


def survey(a: Assumptions) -> dict:
    """Every costed cell of the 2033 survey and the HO ladder, all from eq:Ngate11."""
    env = a.rfi_t_2033.value
    gap = _band(a.gap_mev)
    gap_l = _band(a.gap_large_mev)
    gap_s = _band(a.gap_small_mev)
    # R15 ruling 3: Mo/Ru puts protons in jj44 and neutrons in jj55, one species per shell: 22 + 32 = 54
    # modes (half of each shell's 4 N_orb), i.e. N_orb = 54/4 = 13.5 in eq:Nq11/eq:Ngate11 (was 11 + 16 = 27,
    # 108 modes: both species in both shells). The N_orb^4 string count is then (54/4)^4, 1/16 of 27^4.
    n_mo = n_orb_mo_ru(a)
    cells = {
        "al_sd":   dict(n_proj=1, n_orb=a.n_orb_sd.value,   label=r"$^{27}$Al $sd$"),
        "pf":      dict(n_proj=2, n_orb=a.n_orb_pf.value,   label=r"$^{48}$Ca/Ti ($pf$)"),
        "jj44":    dict(n_proj=2, n_orb=a.n_orb_jj44.value, label=r"$^{76}$Ge/Se, $^{82}$Se/Kr ($jj44$)"),
        "jj55":    dict(n_proj=2, n_orb=a.n_orb_jj55.value, label=r"$^{130}$Te/Xe, $^{136}$Xe/Ba ($jj55$)"),
        "mo_ru":   dict(n_proj=2, n_orb=n_mo,               label=r"$^{100}$Mo/Ru ($jj44\otimes jj55$)"),
    }
    for e in (a.emax_first.lo, a.emax_first.hi, a.emax_converged_al.lo, a.emax_converged_al.hi):
        cells[f"al_emax{e}"] = dict(n_proj=1, n_orb=ho_spatial_orbitals(e),
                                    label=rf"$^{{27}}$Al $e_{{\max}}={e}$")
    for k, c in cells.items():
        c["qubits"] = register(a, 1, c["n_orb"])     # R14 option A: one Fock register for every cell
        c["lq"] = (c["qubits"] + a.n_anc.lo, c["qubits"] + a.n_anc.hi)
        c["n_rot"] = prep_rotations(a, c["n_proj"], c["n_orb"], gap)
        c["eps_rot"] = prep_eps_rot(a, c["n_proj"], c["n_orb"], gap)
        c["c_T"] = prep_c_T(a, c["n_proj"], c["n_orb"], gap)
        c["t"] = prep_t(a, c["n_proj"], c["n_orb"], gap)
        c["multiple"] = _div(c["t"], env)
        c["eps_l"] = eps_l(a, c["t"])
        c["t_large_gap"] = prep_t(a, c["n_proj"], c["n_orb"], gap_l)
        c["multiple_large_gap"] = _div(c["t_large_gap"], env)
        c["t_small_gap"] = prep_t(a, c["n_proj"], c["n_orb"], gap_s)
        c["multiple_small_gap"] = _div(c["t_small_gap"], env)
        c["t_after_lever"] = prep_t(a, c["n_proj"], c["n_orb"], gap, a.prep_lever.value)
        c["multiple_after_lever"] = _div(c["t_after_lever"], env)
    # R15 ruling 4: one controlled BE query of the operator, inside the same circuit as the preparation
    # (R-TOL on the whole circuit), at the selection-rule string count and at the dense count.
    shells = dict(pf=("pf", "pf"), jj44=("jj44", "jj44"), jj55=("jj55", "jj55"), mo_ru=("jj44", "jj55"))
    for k, (sp, sn) in shells.items():
        c = cells[k]
        c["gamma_ins"] = gamma_0nubb(sp, sn)
        c["gamma_ins_dense"] = gamma_0nubb(sp, sn, selection=False)
    cells["al_sd"]["gamma_ins"] = gamma_mu2e_two_body("sd")
    cells["al_sd"]["gamma_ins_dense"] = gamma_mu2e_two_body("sd", selection=False)
    for k in list(shells) + ["al_sd"]:
        c = cells[k]
        for tag, gm in (("ins", c["gamma_ins"]["gamma"]), ("ins_dense", c["gamma_ins_dense"]["gamma"])):
            q = tuple(insertion_in_circuit(a, gm, n) for n in c["n_rot"])
            c[f"{tag}_counts"] = q
            c[f"t_{tag}"] = (q[0]["t"], q[1]["t"])
            # the whole circuit: prep and insertion rotations at the shared tolerance, plus the Toffolis
            c[f"t_with_{tag}"] = tuple(qq["n_rot_circuit"] * qq["c_T"] + qq["t_toffoli"] for qq in q)
        c["insertion_share"] = (c["t_ins"][0] / c["t"][0], c["t_ins"][1] / c["t"][1])
    return cells


# --------------------------------------------------------------------------- #
# r25: T-depth, single-machine walls, and the R7 simplified 0nubb rows
# --------------------------------------------------------------------------- #

def load_depth(d_amp: float, ladder: bool = False) -> int:
    """T-depth (AND layers) of one complete Fomichev load (factory.json, Ch. 3 2028 runs): the amplitude QROM
    2^(L+1) and the lookup D as sequential unary iteration, then D erasure flips, each a (2L-1)-controlled
    AND tree of depth ceil(log2(2L-1)) (compiled) or a 2L-2 Toffoli ladder (ladder)."""
    d = int(round(d_amp))
    L = math.ceil(math.log2(d))
    per_flip = (2 * L - 2) if ladder else math.ceil(math.log2(2 * L - 1))
    return 2 ** (L + 1) + d + d * per_flip


def depth_2028_hadamard(a: Assumptions, gamma: int, d_amp: float, load_angles: bool = False) -> dict:
    """T-depth of one 2028 Hadamard-test shot (factory.json). PREPARE acts on the LCU address only and runs
    alongside the load; the rotation tree pipelines to 2^(k-1) rotations deep (RUS, T-depth c_T each); SELECT
    is Gamma - 1 sequential ANDs. compiled = max(load, P) + SELECT + P^dag (AND-tree erasure); ladder = the
    same with Toffoli-ladder erasure; serial = everything in series as listed. `load_angles` (M3 ruling,
    2026-10-04): the load also carries its L sequential angle additions and the Fourier-state layer."""
    q = insertion_cost(a, gamma, 1)
    k, ct = q["k"], q["c_T"]
    p_pipe = 2 ** (k - 1) * ct
    sel = gamma - a.select_and_deficit.value
    lt, ll = load_depth(d_amp), load_depth(d_amp, ladder=True)
    if load_angles:
        ang = load_angle_cost(a, d_amp, q["n_rot"])
        extra = ang["t_depth_adds"] + ang["t_depth_fourier"]
        lt, ll = lt + extra, ll + extra
    return dict(compiled=max(lt, p_pipe) + sel + p_pipe, ladder=max(ll, p_pipe) + sel + p_pipe,
                serial=ll + 2 * (2 ** k - 1) * ct + sel, prepare_pipelined=p_pipe, load_tree=lt, load_ladder=ll)


def depth_2028_mlae(a: Assumptions, gamma: int, d_amp: float, load_angles: bool = False) -> dict:
    """T-depth of the deepest MLAE circuit (factory.json): A S_chi A^dag S_0 ... A. Across each S_chi (test
    qubit only) P^dag P cancels; S_0 is a log-depth AND tree over the reflected qubits. `load_angles` (M3
    ruling, 2026-10-04): every application of A and A^dag carries the L sequential angle additions, and the
    circuit one Fourier-state layer (the state is reused)."""
    ae = mlae_2028(a, gamma, d_amp, load_angles=load_angles)
    m, ct = ae["m_max"], ae["c_T"]
    k = math.ceil(math.log2(gamma))
    p_pipe = 2 ** (k - 1) * ct
    sel = gamma - a.select_and_deficit.value
    lt, ll = load_depth(d_amp), load_depth(d_amp, ladder=True)
    once = 0.0
    if load_angles:
        ang = load_angle_cost(a, d_amp, m * insertion_cost(a, gamma, 1)["n_rot"])
        lt, ll = lt + ang["t_depth_adds"], ll + ang["t_depth_adds"]
        once = ang["t_depth_fourier"]
    pairs = (m - 1) // 2
    s0 = math.ceil(math.log2(ae["reflection_qubits"]))
    compiled = pairs * (2 * max(lt, p_pipe) + 2 * sel) + pairs * s0 + (max(lt, p_pipe) + sel + p_pipe) + once
    no_cancel = m * (max(lt, p_pipe) + sel + p_pipe) + pairs * s0 + once
    serial = m * (ll + 2 * (2 ** k - 1) * ct + sel) + pairs * ae["reflection_ands"] + once
    return dict(compiled=compiled, no_cancellation=no_cancel, serial=serial, m_max=m, t_deepest=ae["t_deepest"])


def depth_2033(a: Assumptions, cnt: dict, gamma: int, conc: float) -> float:
    """T-depth of one 2033 shot (factory.json): N_rot,prep c_T / G for the Trotter prep with G commuting
    rotations in flight, plus the insertion tail (SELECT Gamma - 1 ANDs and one pipelined P^dag; PREPARE
    overlaps the prep). G = 1 returns the full T-count (gate by gate)."""
    t_all = cnt["n_rot_circuit"] * cnt["c_T"] + cnt["t_toffoli"]
    if conc <= 1:
        return t_all
    n_rot_prep = cnt["n_rot_circuit"] - cnt["n_rot"]
    tail = (gamma - a.select_and_deficit.value) + 2 ** (cnt["k"] - 1) * cnt["c_T"]
    return n_rot_prep * cnt["c_T"] / conc + tail


def wall_one_machine(a: Assumptions, shots, t, d) -> tuple[float, float]:
    """Serial wall on ONE machine with the 10-factory baseline, end by end: shots x (max(N_T t_gate, D_T t_r)
    + t_0). Equals the T-supply wall when F* >= 10, the depth-limited (corrected) wall otherwise."""
    sh = shots if isinstance(shots, tuple) else (shots, shots)
    return tuple(n * (max(x * a.t_gate_s.value, y * REACTION_TIME_S) + a.shot_overhead_s.value)
                 for n, x, y in zip(sh, t, d))


def _pair_prep_two_gaps(a: Assumptions, n_orb: float, gaps: tuple, shells: str) -> dict:
    """R7 (simplobs price.py): a 0nubb pair whose parent and daughter projections each run c_proj/Delta_i,
    summed, then one controlled BE query in the circuit (R-TOL on the whole circuit). Band ends: c_proj = pi
    and 2 pi."""
    dt = a.dt_gev_inv.value
    out = []
    for cp in (a.c_proj.lo, a.c_proj.hi):
        nt = sum(cp * tau_gs(g) / dt for g in gaps)
        n_rot = nt * n_orb ** 4
        gm = gamma_0nubb(shells, shells)["gamma"]
        q = insertion_in_circuit(a, gm, n_rot)
        out.append(dict(n_trotter=nt, n_rot=n_rot, gamma=gm, t=q["n_rot_circuit"] * q["c_T"] + q["t_toffoli"],
                        t_depth=depth_2033(a, q, gm, a.trotter_concurrency.value)))
    return dict(n_trotter=(out[0]["n_trotter"], out[1]["n_trotter"]), gamma=out[0]["gamma"],
                t=(out[0]["t"], out[1]["t"]), t_depth=(out[0]["t_depth"], out[1]["t_depth"]))


def _rabi_cost(a: Assumptions, n_orb: float, lam_over_m) -> dict:
    """R7 resonant-drive (Rabi) readout (simplobs price_rabi.py): evolve the closed-f7/2 reference under
    H + mu N_p + theta (O + O^dag) to t_pi/2 = (pi/2)/(theta|M|), then read proton number (Clifford). Trotter
    steps t/dt at N_orb^4 rotations each, plus the qDRIFT-sampled drive, N_q = 2 (lambda/M)^2 (pi/2)^2 / eps_q
    single-Pauli rotations. Low end: strong drive and the smaller lambda/M; high end: weak drive, larger."""
    lm = lam_over_m if isinstance(lam_over_m, tuple) else (lam_over_m, lam_over_m)
    ends = []
    for th, r in ((a.rabi_drive_mev.hi, lm[0]), (a.rabi_drive_mev.lo, lm[1])):
        steps = (PI / 2) / th * 1e3 / a.dt_gev_inv.value
        nq = 2 * r ** 2 * (PI / 2) ** 2 / a.qdrift_eps.value
        n_trot = steps * n_orb ** 4
        n_rot = n_trot + nq
        ct = c_T(a, n_rot)
        # T-depth (R10): the Trotter part runs G commuting rotations in flight; the qDRIFT samples are random,
        # mostly non-commuting single-Pauli rotations, one sequential layer each (grouping does not help).
        d = n_trot * ct / a.trotter_concurrency.value + nq * ct
        ends.append(dict(steps=steps, n_qdrift=nq, n_trotter_rot=n_trot, n_rot=n_rot, c_T=ct, t=n_rot * ct,
                         t_depth=d))
    return {k: (ends[0][k], ends[1][k]) for k in ("steps", "n_qdrift", "n_trotter_rot", "c_T", "t", "t_depth")}


def simplified_0nubb(a: Assumptions) -> dict:
    """R7 rows (simplobs.json Ch. 3): 48Ca -> 48Ti first light in f7/2 p3/2 (projection readout at the 0+
    gaps, 30% on |M|), then f7/2 p3/2 p1/2 and full pf with the Rabi readout. Every wall is serial on one
    machine at t_gate_s per T + shot_overhead_s per shot, depth-corrected (wall_one_machine): the Trotter part
    runs trotter_concurrency rotations in flight, the qDRIFT samples one at a time. Registers: 4 N_orb modes + N_anc."""
    env = a.rfi_t_2033.value
    # referee M3: each space at its own 0+ gaps (KB3G restricted to the kept orbits), not the ENSDF gaps
    gaps_f7p3 = (a.gap_0plus_f7p3_ca48_mev.value, a.gap_0plus_f7p3_ti48_mev.value)
    gaps_f7p3p1 = (a.gap_0plus_f7p3p1_ca48_mev.value, a.gap_0plus_f7p3p1_ti48_mev.value)
    gaps_pf = (a.gap_0plus_pf_ca48_mev.value, a.gap_0plus_pf_ti48_mev.value)
    gaps_ensdf = (a.gap_0plus_ensdf_ca48_mev.value, a.gap_0plus_ensdf_ti48_mev.value)
    t0, tg = a.shot_overhead_s.value, a.t_gate_s.value
    rows = {}
    # B: f7/2 p3/2, projection readout. Shots (lambda/|M|)^2 / (4 eps^2 p_ref) (app11 operator-insertion route)
    b = _pair_prep_two_gaps(a, 6, gaps_f7p3, "f7p3")
    b_ensdf = _pair_prep_two_gaps(a, 6, gaps_ensdf, "f7p3")                       # record: priced until M3
    shots_b = a.lam_over_m_f7p3.value ** 2 / (4 * a.eps_simplified.value ** 2 * a.p_ref_f7p3.value)
    rows["f7p3_projection"] = dict(modes=24, lq=(24 + a.n_anc.lo, 24 + a.n_anc.hi), t=b["t"], gamma=b["gamma"],
                                   n_trotter=b["n_trotter"], multiple=_div(b["t"], env), shots=(shots_b, shots_b),
                                   eps_l=eps_l(a, b["t"]), t_depth=b["t_depth"], f_star_per_end=(b["t"][0] / b["t_depth"][0], b["t"][1] / b["t_depth"][1]),
                                   wall_s=wall_one_machine(a, shots_b, b["t"], b["t_depth"]),
                                   wall_t_supply_s=(shots_b * (b["t"][0] * tg + t0), shots_b * (b["t"][1] * tg + t0)),
                                   gaps=gaps_f7p3, t_at_ensdf_gaps=b_ensdf["t"],
                                   wall_at_ensdf_gaps_s=wall_one_machine(a, shots_b, b_ensdf["t"], b_ensdf["t_depth"]))
    # B': f7/2 p3/2 p1/2 with the Rabi readout (and the projection readout it replaces, for the record)
    bp = _rabi_cost(a, 7, a.lam_over_m_f7p3p1.value)
    bp_proj = _pair_prep_two_gaps(a, 7, gaps_f7p3p1, "f7p3p1")
    sh = _band(a.rabi_shots)
    rows["f7p3p1_rabi"] = dict(modes=28, lq=(28 + a.n_anc.lo, 28 + a.n_anc.hi), t=bp["t"], steps=bp["steps"],
                               n_qdrift=bp["n_qdrift"], multiple=_div(bp["t"], env), shots=sh,
                               eps_l=eps_l(a, bp["t"]), t_depth=bp["t_depth"], f_star_per_end=(bp["t"][0] / bp["t_depth"][0], bp["t"][1] / bp["t_depth"][1]),
                               wall_s=wall_one_machine(a, sh, bp["t"], bp["t_depth"]),                  # 2.0-11 d
                               wall_t_supply_s=(sh[0] * (bp["t"][0] * tg + t0), sh[1] * (bp["t"][1] * tg + t0)),
                               projection_t=bp_proj["t"], projection_gaps=gaps_f7p3p1)
    # H: full pf with the Rabi readout: a near miss in T, weeks in wall
    h = _rabi_cost(a, a.n_orb_pf.value, _band(a.lam_over_m_pf))
    rows["pf_rabi"] = dict(modes=40, lq=(40 + a.n_anc.lo, 40 + a.n_anc.hi), t=h["t"], steps=h["steps"],
                           n_qdrift=h["n_qdrift"], multiple=_div(h["t"], env), shots=sh,
                           eps_l=eps_l(a, h["t"]), t_depth=h["t_depth"], f_star_per_end=(h["t"][0] / h["t_depth"][0], h["t"][1] / h["t_depth"][1]),
                           wall_s=wall_one_machine(a, sh, h["t"], h["t_depth"]),                           # 19-110 d
                           wall_t_supply_s=(sh[0] * (h["t"][0] * tg + t0), sh[1] * (h["t"][1] * tg + t0)))
    # A (reference): full pf, projection readout at the 0+ gaps: T near the envelope, shots out of reach
    a_pf = _pair_prep_two_gaps(a, a.n_orb_pf.value, gaps_pf, "pf")
    a_pf_ensdf = _pair_prep_two_gaps(a, a.n_orb_pf.value, gaps_ensdf, "pf")    # record: priced until M3
    lm = _band(a.lam_over_m_pf)
    sh_a = tuple(r ** 2 / (4 * a.eps_simplified.value ** 2) for r in lm)          # p_ref = 1: a lower bound
    rows["pf_projection_0plus"] = dict(t=a_pf["t"], multiple=_div(a_pf["t"], env), gamma=a_pf["gamma"],
                                       shots_lower_bound=sh_a,
                                       wall_s=(sh_a[0] * (a_pf["t"][0] * tg + t0), sh_a[1] * (a_pf["t"][1] * tg + t0)),
                                       gaps=gaps_pf, t_at_ensdf_gaps=a_pf_ensdf["t"],
                                       multiple_at_ensdf_gaps=_div(a_pf_ensdf["t"], env))
    return rows


# --------------------------------------------------------------------------- #
# The model
# --------------------------------------------------------------------------- #

def _model_2028(a: Assumptions) -> Result:
    cap = a.rfi_t_2028.value
    system = register(a, 1, a.n_orb_he4.value)                      # 40
    lq = (system + a.n_anc_2028.lo, system + a.n_anc_2028.hi)       # 120-170
    # R13: the Hadamard test on the controlled BE is the degree-1 polynomial, one query, no phases
    q = a.be_queries_2028.value                                                                  # 1
    topdown = (a.t_be_topdown.lo * q, a.t_be_topdown.hi * q)                                     # 1e3-1e4
    # Bottom-up, derived here: Gamma<=128 LCU insertion, one controlled query
    g = a.gamma_trunc.value
    one = insertion_cost(a, g, q)
    t_ins = (one["t"], one["t"])                                                                 # 5362.9
    frac128 = _div(t_ins, cap)                                                                   # 0.0536
    t_per_term = _div(t_ins, g)
    g256 = insertion_cost(a, a.gamma_256.value, q)["t"]
    frac256 = (g256 / cap, g256 / cap)                                                           # 0.111
    g512 = insertion_cost(a, a.gamma_512.value, q)["t"]
    frac512 = (g512 / cap, g512 / cap)                                                           # 0.228
    pairs = math.comb(system, 2)
    gamma_full = (2 * pairs + system, 4 * pairs + system)                                        # 1600, 3160
    full = (insertion_cost(a, gamma_full[0], q), insertion_cost(a, gamma_full[1], q))
    full_mult = (full[0]["t"] / cap, full[1]["t"] / cap)                                         # 0.927-1.900
    # Sensitivity: the pre-R13 literal reading, d = ceil(log2 1/eps) = 2 (eps=0.3) to 4 (eps=0.1)
    # queries plus d + 1 phase steps (R12 print), and the doubled degree (walk W and W^dag, or a
    # real-part LCU of Phi and -Phi)
    qsp = (math.ceil(math.log2(1 / a.eps_2028.hi)), math.ceil(math.log2(1 / a.eps_2028.lo)))   # 2-4
    qsp_counts = (insertion_cost(a, g, qsp[0]), insertion_cost(a, g, qsp[1]))
    qsp_t = (qsp_counts[0]["t"], qsp_counts[1]["t"])                                             # 11201-22928
    qsp_doubled_t = (insertion_cost(a, g, 2 * qsp[0])["t"], insertion_cost(a, g, 2 * qsp[1])["t"])
    frac_8q = insertion_cost(a, g, 8)["t"] / cap                                                 # 0.470
    # what the pre-R12 print 0.54-0.57 was: 8 queries x (254 rotations at a fixed eps_rot = 1e-4,
    # 24.48 T each, + 127 Toffolis at 4 T (low end) or 7 T (high end)) = 53809 / 56857
    pre_r12 = tuple(8 * (2 * (g - 1) * t_per_rotation(1e-4) + (g - 1) * tt) / cap for tt in (4, 7))
    # insertion-stage register: 40 system + k address + k AND workspace (controlled unary iteration) + 1 test;
    # no phase ancilla (no phase steps at d = 1)
    peak_ins = system + one["address_qubits"] + one["select_and_ancillae"] + a.test_qubits.value      # 55
    # ---- Oracle load (R15 ruling 2): the complete sparse load of Fomichev et al. arXiv:2310.18410, priced at
    # 7 T per Toffoli, one load per shot (it erases its own index; the test measures). The R11b-R14
    # lookup-only price (7(D-1) per QROM pass, 2-3 passes) is kept below as the record.
    t_tof = a.insertion_t_per_toffoli.value
    n_load = a.loads_per_shot.value                                                              # 1

    def load_t(d):
        return n_load * t_tof * sos_load_toffoli(d)
    d_pr = a.d_amplitudes.value                                                                  # 600 (printed)
    d_range = (a.d_plausible.lo, d_pr)                                                           # 1e2 .. 600
    # crossing without the angle rotations (record since the M3 ruling): insertion + load <= cap
    d_break_no_angles = max(d for d in range(2, 5000) if t_ins[0] + load_t(d) <= cap)            # 603 (last fit)
    d_cross_no_angles = d_break_no_angles + 1                                                    # 604 (first over)
    load_at = {d: load_t(d) for d in (a.d_plausible.lo, d_pr, a.d_plausible.hi)}                 # 10892, 94136, 147336
    load_range = (load_t(d_range[0]), load_t(d_range[1]))                                        # 1.1e4 - 9.4e4
    composed_range = ((t_ins[0] + load_range[0]) / cap, (t_ins[1] + load_range[1]) / cap)       # 0.16 - 0.99
    t_with_load = _mul(composed_range, cap)
    composed_at_d_hi = (t_ins[0] + load_at[a.d_plausible.hi]) / cap                              # 1.53 at D = 1e3
    sos_over_lookup = {d: sos_load_toffoli(d) / (d - 1) for d in load_at}                        # 15.7, 22.4, 21.1
    # Referee M3/G2 (2026-10-04): the load's amplitude-angle rotations (derived here, load_angle_cost). Author
    # ruling (M3, 2026-10-04, "do 3"; Claude's recommended default): FOLDED into the headline, the walls, the
    # T-depth and the register. The executed shot is insertion + complete load + angle rotations.
    n_rot_ins = one["n_rot"]
    ang = {d: load_angle_cost(a, d, n_rot_ins) for d in (d_range[0], d_pr)}
    ang_t = (ang[d_range[0]]["t"], ang[d_pr]["t"])                                               # 766, 1018
    with_angles = (t_ins[0] + load_range[0] + ang_t[0], t_ins[1] + load_range[1] + ang_t[1])    # 1.70e4, 1.005e5
    d_break_angles = max(d for d in range(2, 5000)
                         if t_ins[0] + load_t(d) + load_angle_cost(a, d, n_rot_ins)["t"] <= cap)  # 596
    d_break, d_cross = d_break_angles, d_break_angles + 1                                        # 596, 597
    # referee M3: the computed (schematic, Minnesota) 4He D is at or below the priced band (both fidelities) and fits
    d_he4 = dict(fid99=a.he4_d_fid99.value, fid999=a.he4_d_fid999.value, dipole_1pct=a.he4_d_dipole_1pct.value,
                 at_or_below_priced_band=max(a.he4_d_fid99.hi, a.he4_d_fid999.hi) <= d_range[1],
                 below_priced_band_low_end=a.he4_d_fid99.hi < d_range[0] and a.he4_d_fid999.lo < d_range[0],
                 fits=a.he4_d_fid999.hi <= d_break)
    t_no_angles = (t_ins[0] + load_range[0], t_ins[1] + load_range[1])                           # 16255-99499 (record)
    # register (ruling 3 rule, R15): size the load to the printed D condition, D <~ 600: L = 10 enumeration
    # qubits; the source's ancillae 5L - 3 = 47 (L enumeration + 2L-1 identification + 2L-2 multi-control
    # workspace, which also hosts the lookup's L - 1 ANDs). The lookup writes into the system register.
    addr = math.ceil(math.log2(d_pr))                                                            # 10
    load_anc = sos_load_ancillae(d_pr)                                                           # 47
    peak_load_no_angles = system + load_anc                                                      # 87 (record)
    # M3 ruling: step 1 of the load holds L enumeration + 3b - 1 angle/Fourier/carry qubits = 48 > 47
    peak_load = system + max(load_anc, addr + ang[d_pr]["qubits_step1"])                         # 88
    peak = max(peak_load, peak_ins)                                                              # 88 (load)
    # ---- record: the R11b-R14 lookup-only price (7(D-1) per pass, 2-3 passes) at the old D = 2e3 ----
    k = a.qrom_t_per_amplitude.value
    lookup_pass_2e3 = k * 2000 - k                                                               # 13993
    lookup_d_break = (1 + (cap - t_ins[1]) / (a.qrom_passes.hi * k),                             # 4507.5
                      1 + (cap - t_ins[0]) / (a.qrom_passes.lo * k))                             # 6760.8
    lookup_composed_2e3 = ((t_ins[0] + a.qrom_passes.lo * lookup_pass_2e3) / cap,                # 0.33
                           (t_ins[1] + a.qrom_passes.hi * lookup_pass_2e3) / cap)                # 0.47
    lookup_load_plausible = (a.qrom_passes.lo * (k * a.d_plausible.lo - k),                      # 1386
                             a.qrom_passes.hi * (k * a.d_plausible.hi - k))                      # 20979
    # R14 item 1a / R15 ruling 1: the Hadamard-test shot factor (lambda'/|<O>|)^2, identity subtracted
    dip = he4_dipole(a)
    # R16: the shots use the upper end of the any-state band (32.76), the convention of the 27Al cell; the
    # closed-shell 0s^4 value (19.85) is not a 4He ground state in this space and is kept as the reference.
    factor = dip["shot_factor_band"][1]                                                          # 32.76
    factor_closed_shell = dip["shot_factor"]                                                     # 19.85
    # insertion cost of the dipole operator itself (Gamma = 24 strings after the identity, no truncation)
    dip_ins = insertion_cost(a, dip["gamma_lcu"], q)                                             # 1.2e3
    # shots and wall time.  eps^-2 ~ 1e2 per operator (eps = 0.1); relative accuracy eps on <O> needs
    # (lambda'/|<O>|)^2 eps^-2 shots.  The dipole's factor is applied to all six categories.
    shots_per_op_bare = round(1 / a.eps_2028.lo ** 2)                                           # 1e2
    shots_bare = shots_per_op_bare * a.n_op_mu2e.value                                          # 6e2 (R13 print)
    shots_per_op = shots_per_op_bare * factor                                                   # 3276
    shots = shots_per_op * a.n_op_mu2e.value                                                     # 2.0e4
    shots_closed_shell = shots_per_op_bare * factor_closed_shell * a.n_op_mu2e.value             # 1.2e4 (R15)
    shots_r14 = shots_per_op_bare * dip["shot_factor_with_identity"] * a.n_op_mu2e.value         # 4.6e4 (R14)
    # a shot is the insertion plus one complete load with its angle rotations at D = 1e2 .. 600 (17020 - 100516 T)
    t_shot_ins_s = _mul(t_ins, a.t_gate_s.value)                                         # 0.0054 s
    # r17: plus the fixed per-shot overhead (init + readout + decode), 0.1 ms
    t0 = a.shot_overhead_s.value
    t_shot_s = (with_angles[0] * a.t_gate_s.value + t0, with_angles[1] * a.t_gate_s.value + t0)  # 0.017-0.101 s
    wall = _mul(t_shot_s, shots)                                                                 # 339-1975 s
    wall_no_overhead = ((t_shot_s[0] - t0) * shots, (t_shot_s[1] - t0) * shots)                 # 320-1956 s
    wall_bare = _mul(t_shot_s, shots_bare)
    wall_ins = _mul(t_shot_ins_s, shots)
    # ---- depth-capped MLAE at the box's Gamma <= 128 insertion, D = 1e2 .. 600 (R15: complete load; M3 ruling:
    # with the load-angle rotations, the phase-gradient state charged once per circuit, catalytic).
    # At D = 600 no circuit fits (1.005x the cap at m = 1); m_max = 1 is the plain test, reported as such.
    def _ae(gm, d):
        return mlae_2028(a, gm, d, load_angles=True)
    ae = (_ae(g, d_range[0]), _ae(g, d_break))                                                   # m = 5, 1 (D = 596)
    ae_steps = {}                                                                                # D where m_max drops
    prev = None
    for d in range(int(d_range[0]), d_break + 1):
        m = _ae(g, d)["m_max"]
        if m != prev:
            ae_steps[m] = d
            prev = m
    ae_last_d = {m: max(d for d in range(int(d_range[0]), d_break + 1) if _ae(g, d)["m_max"] >= m)
                 for m in ae_steps}                                                              # {5: 128, 3: 221, 1: 596}
    ae_no_angles = (mlae_2028(a, g, d_range[0]), mlae_2028(a, g, d_range[1]))                    # record: 8.4e4 at D = 1e2
    ae_last_d_no_angles = {m: max(d for d in range(int(d_range[0]), d_break_no_angles + 1)
                                  if mlae_2028(a, g, d)["m_max"] >= m) for m in (5, 3, 1)}       # {5: 128, 3: 228, 1: 603}
    ae_dipole = (_ae(dip["gamma_lcu"], d_range[0]), _ae(dip["gamma_lcu"], d_break))
    c_est = a.mlae_estimator_constant.value
    ae_circuits_per_op = tuple(c_est * shots_per_op / x["m_max"] ** 2 for x in ae)               # 131, 3276
    ae_circuits = tuple(n * a.n_op_mu2e.value for n in ae_circuits_per_op)                        # 786, 19656
    ae_t_circuit_s = tuple(x["t_deepest"] * a.t_gate_s.value + t0 for x in ae)          # r17: + t0
    ae_wall = (ae_circuits[0] * ae_t_circuit_s[0], ae_circuits[1] * ae_t_circuit_s[1])
    ae_lis = tuple(mlae_lis_constant(x["m_max"]) for x in ae)
    ae_queries_per_op = tuple(n * x["m_max"] for n, x in zip(ae_circuits_per_op, ae))
    # all-at-m_max needs no shallow disambiguation circuits: the state-independent band on the test's mean
    # (<O> - c_0)/lambda' maps to a window of m theta inside one monotone branch of sin^2(m theta),
    # sin^2 theta = (1 + mean)/2 = P(test qubit reads 0).
    def _mtheta(m, a_ratio):
        return m * math.asin(math.sqrt((1 + a_ratio) / 2))
    a_band = dip["test_mean_band"]
    ae_mtheta = tuple((_mtheta(x["m_max"], a_band[0]), _mtheta(x["m_max"], a_band[1])) for x in ae)
    ae_one_branch = tuple(math.floor(w[0] / (PI / 2)) == math.floor(w[1] / (PI / 2)) for w in ae_mtheta)
    ae_branch_margin = tuple(min(w[0] - math.floor(w[0] / (PI / 2)) * PI / 2,
                                 (math.floor(w[1] / (PI / 2)) + 1) * PI / 2 - w[1]) for w in ae_mtheta)
    # record: the R14 lookup-only MLAE (2 passes per application) at D = 1e2 (m_max = 13 then)
    ae_lookup_r14 = mlae_2028(a, g, a.d_plausible.lo, t_load_override=2 * k * (a.d_plausible.lo - 1))
    # ---- r25 (R9): T-depth per shot (factory.json, derived in depth_2028_*), F*, and the corrected walls.
    # The Hadamard-test run is the headline: N_T = insertion + one load at D = 1e2 .. 600 (16255 - 99499 T),
    # T-depth compiled 2381.6 - 6902.3 (AND-tree erasure); ladder 2810 - 14702, fully serial 6157 - 18049.
    dh = (depth_2028_hadamard(a, g, d_range[0], True), depth_2028_hadamard(a, g, d_range[1], True))
    dh_no_angles = (depth_2028_hadamard(a, g, d_range[0]), depth_2028_hadamard(a, g, d_range[1]))
    t_run = with_angles                                                                           # 17020-100516
    d_run = (dh[0]["compiled"], dh[1]["compiled"])                                                # 2382-6902
    ex = depth_exports(t_run, d_run, shots, a.t_gate_s, a.shot_overhead_s)
    f_star_end = (t_run[0] / d_run[0], t_run[1] / d_run[1])                                       # 6.8, 14.4
    wall_corr = wall_one_machine(a, shots, t_run, d_run)                                          # 470-1958 s
    wall_ladder = wall_one_machine(a, shots, t_run, (dh[0]["ladder"], dh[1]["ladder"]))
    wall_serial_compile = wall_one_machine(a, shots, t_run, (dh[0]["serial"], dh[1]["serial"]))
    dm = depth_2028_mlae(a, g, d_range[0], load_angles=True)
    ae_wall_corr = wall_one_machine(a, ae_circuits[0], (dm["t_deepest"],) * 2, (dm["compiled"],) * 2)[0]  # 66 s
    ae_f_star = dm["t_deepest"] / dm["compiled"]                                                  # 10.6
    # referee M3 (verifier, 2026-10-04): the MLAE solve with the load-angle rotations, the phase-gradient
    # state charged once per circuit (catalytic): m_max = 5 at D = 1e2, deepest 83754 + 5 x 588 + 177 = 8.7e4;
    # m_max = 5 still to D = 128 (L steps at 129), 3 to D = 221, then 1 (228 without the angle term).
    # Since the M3 ruling this is the headline MLAE (ae above); the names are kept for the tests.
    ae_ang = ae[0]
    mlae_with_angles = ae_ang["t_deepest"]                                                         # 86893
    ae_last_d_angles = ae_last_d                                                                  # {5: 128, 3: 221, 1: 596}
    # MLAE register: the reused b-qubit phase-gradient state is held through the reflection (114 + 13)
    mlae_peak = ae[0]["peak_lq_reflection"] + ang[d_range[0]]["b"]                               # 127
    note_eps = f"R-TOL eps_rot = sqrt(1e-2/{one['n_rot']}) = {one['eps_rot']:.3e} -> {one['c_T']:.3f} T"
    breakdown = (   # referee G2 (2026-10-04): the executed shot, insertion + one complete load, at D = 1e2
        Primitive("sos_load_complete", sos_load_toffoli(d_range[0]), t_tof, CircuitStatus.COMPILED,
                  "arxiv_2310_18410",
                  f"complete sparse load, (2L-2)D + 2^(L+1) + D Toffolis at D = {int(d_range[0])} (initial.tex:350); "
                  "7 T each (R5); amplitude-angle rotations not itemized by the source (next two rows)"),
        Primitive("load_angle_cadd", ang[d_range[0]]["toffoli"], t_tof, CircuitStatus.COMPILED,
                  "arxiv_1812_00954",
                  f"amplitude angles: L = {ang[d_range[0]]['L']} phase-gradient additions of b = {ang[d_range[0]]['b']} "
                  "bits, b - 1 Toffolis each (arXiv:1709.06648), 7 T each (R5); M3 ruling: in the headline"),
        Primitive("load_fourier_state", ang[d_range[0]]["n_rot"], ang[d_range[0]]["c_T"], CircuitStatus.COMPILED,
                  "arxiv_1812_00954",
                  "b-qubit phase-gradient state: b - 3 synthesized Z rotations (R-TOL with the insertion)"),
        Primitive("load_fourier_state_t", 1, 1.0, CircuitStatus.COMPILED, "arxiv_1812_00954",
                  "the pi/4 rotation of the phase-gradient state: one T"),
        Primitive("prepare_rotation_tree", one["n_rot_prepare"], one["c_T"], CircuitStatus.COMPILED, "derived here",
                  f"PREPARE + PREPARE^dag, 2^k - 1 = {one['rot_per_prepare']} R_y each (k = {one['k']}), "
                  f"one query; {note_eps}"),
        Primitive("select_unary_iteration", one["toffoli_select"], a.insertion_t_per_toffoli.value,
                  CircuitStatus.COMPILED, "Babbush_PRX_2018",
                  f"controlled unary iteration, Gamma - 1 = {g - 1} Toffolis (arXiv:1805.03662 "
                  "main_draft.tex:716), one query; 7 T each (R5)"),
        Primitive("pauli_string_select_targets", g, 0.0, CircuitStatus.COMPILED, "derived here",
                  "JW Pauli string (with its coefficient sign) under the unary flag: Clifford"),
        Primitive("hadamard_test_frame", 1, 0.0, CircuitStatus.COMPILED, "derived here",
                  "H on the test qubit and the X readout (<O> is real, no Y run): Clifford"),
    )
    return Result(
        era="2028", lq=lq, hard_ops=with_angles, breakdown=breakdown,
        intermediates=dict(
            system_qubits=system, n_anc=_band(a.n_anc_2028), lq_sum=lq,
            be_queries=q, qsp_phase_steps=one["n_rot_phase"], t_shot_topdown=topdown,
            gamma_trunc_fraction=frac128, t_insertion=t_ins, t_per_lcu_term=t_per_term,
            insertion_counts=one, gamma_256_fraction=frac256, gamma_512_fraction=frac512,
            gamma_full_multiple=full_mult, gamma_full_counts=full, gamma_full=gamma_full,
            qsp_reading_degree=qsp, qsp_reading_counts=qsp_counts, qsp_reading_t=qsp_t,
            qsp_reading_fraction=_div(qsp_t, cap), qsp_reading_doubled_fraction=_div(qsp_doubled_t, cap),
            gamma128_fraction_at_8_queries=frac_8q,
            pre_r12_printed_reconstruction=pre_r12,
            loads_per_shot=n_load, d_priced=d_pr, d_range=d_range, d_break=d_break, d_cross=d_cross,
            d_he4_computed=d_he4,
            d_break_no_angles=d_break_no_angles, d_cross_no_angles=d_cross_no_angles,
            hard_ops_no_angles=t_no_angles, mlae_no_angles=ae_no_angles,
            mlae_m_max_last_d_no_angles=ae_last_d_no_angles, t_depth_hadamard_no_angles=dh_no_angles,
            peak_lq_load_no_angles=peak_load_no_angles,
            load_t_at_d=load_at, load_t_range=load_range, composed_fraction_range=composed_range,
            hard_ops_with_oracle_load=t_with_load, composed_fraction_at_d_hi=composed_at_d_hi,
            sos_over_lookup_pass=sos_over_lookup,
            load_angle=ang, load_angle_t=ang_t, hard_ops_with_load_angles=with_angles,
            fraction_with_load_angles=_div(with_angles, cap), d_break_with_load_angles=d_break_angles,
            mlae_t_deepest_with_load_angles=mlae_with_angles, mlae_with_load_angles=ae_ang,
            mlae_m_max_last_d_with_load_angles=ae_last_d_angles,
            max_subroutine_time_s=(wall_corr[0] / shots, wall_corr[1] / shots),
            load_address_qubits=addr, load_ancillae=load_anc,
            peak_lq_load=peak_load, peak_lq_insertion=peak_ins, peak_lq=peak,
            lookup_r14=dict(pass_t_at_2e3=lookup_pass_2e3, d_break=lookup_d_break,
                            composed_fraction_at_2e3=lookup_composed_2e3, load_plausible=lookup_load_plausible,
                            passes=_band(a.qrom_passes), mlae_at_d_lo=ae_lookup_r14),
            dipole=dip, dipole_lambda=dip["lam"], dipole_lambda_lcu=dip["lam_lcu"],
            dipole_expectation=dip["exp_closed_shell"],
            dipole_gamma=dip["gamma"], dipole_gamma_lcu=dip["gamma_lcu"], dipole_shot_factor=factor,
            dipole_shot_factor_closed_shell=factor_closed_shell, shots_closed_shell=shots_closed_shell,
            dipole_ratio=dip["ratio"], dipole_shot_factor_band=dip["shot_factor_band"],
            dipole_shot_factor_with_identity=dip["shot_factor_with_identity"],
            dipole_insertion_counts=dip_ins, dipole_insertion_t=dip_ins["t"],
            shots_per_operator_bare=shots_per_op_bare, shots_bare=shots_bare, wall_time_bare_s=wall_bare,
            shot_overhead_s=t0, wall_time_no_overhead_s=wall_no_overhead,
            shots_per_operator=shots_per_op, shots=shots, shots_r14=shots_r14, t_shot_s=t_shot_s,
            t_shot_insertion_only_s=t_shot_ins_s, wall_time_insertion_only_s=wall_ins,
            mlae=ae, mlae_m_max=(ae[0]["m_max"], ae[1]["m_max"]), mlae_m_max_first_d=ae_steps,
            mlae_m_max_last_d=ae_last_d,
            mlae_dipole_gamma=ae_dipole, mlae_estimator_constant=c_est,
            mlae_circuits_per_operator=ae_circuits_per_op, mlae_circuits=ae_circuits,
            mlae_t_circuit_s=ae_t_circuit_s, mlae_wall_time_s=ae_wall, mlae_lis_constant=ae_lis,
            mlae_queries_per_operator=ae_queries_per_op, mlae_mtheta_window=ae_mtheta,
            mlae_window_in_one_branch=ae_one_branch, mlae_branch_margin_rad=ae_branch_margin,
            mlae_peak_lq_reflection=ae[0]["peak_lq_reflection"], mlae_peak_lq=mlae_peak,
            mlae_epsilon_l_deepest=(a.fault_budget.value / ae[0]["t_deepest"],
                                    a.fault_budget.value / ae[1]["t_deepest"]),
            fq_light_nucleus_over_cap=a.fq_light_nucleus_t.value / cap,
            rfi_lq_2028=_band(a.rfi_lq_2028),
            # r25 (R9) exports: the Hadamard-test run at D = 1e2 .. 600; single tier (benchmark)
            **{k: ex[k] for k in ("t_per_shot", "t_depth_per_shot", "f_star", "floor_wall_s", "factories_for_1yr")},
            wall_first_result_s=None, wall_campaign_s=wall_corr,
            depth_exports_hadamard=ex, f_star_per_end=f_star_end, t_depth_hadamard=dh,
            wall_time_serial_supply_s=wall, wall_time_corrected_s=wall_corr,
            wall_time_ladder_erasure_s=wall_ladder, wall_time_serial_compile_s=wall_serial_compile,
            mlae_depth=dm, mlae_f_star=ae_f_star, mlae_wall_time_corrected_s=ae_wall_corr,
        ),
        shots=(shots, shots), wall_time_s=wall_corr, epsilon_l=eps_l(a, t_run),
        notes=("Referee G2/M3 (2026-10-04, author ruling 'do 3'): hard_ops is the executed Hadamard-test shot, the "
               "one-query insertion plus one complete load with its amplitude-angle rotations at D = 1e2 .. 600 "
               "(1.7e4-1.0e5 T; 1.005x the cap at D = 600, last fit D = 596). Nothing in the shot is left out, "
               "so it is not a lower bound; it is conditional on D (no calculation fixes D for 4He). Without the "
               "angle term (record): 1.6e4-9.9e4, hard_ops_no_angles.",
               "State prep algorithm oracle-loaded (app11:122): the load is counted, the preparation is not.",
               "epsilon_l is 0.1/N (app11:93, R3); the 2028 box prints no eps_l row.",
               "hard_ops is the Gamma<=128 insertion counted gate by gate (derived here): one controlled BE "
               "query, 5.4e3 T (R13). R12 read the text's 2-4 'QSP phases' as 2-4 queries (1.1-2.3e4), "
               "kept as qsp_reading_*; the unsourced 5.4-5.7e4 before that.",
               "The oracle load is printed separately in the box (R11 (A)). R15: the complete sparse load of "
               "Fomichev 2310.18410, (2L-2)D + 2^(L+1) + D Toffolis, one per shot: 1.1e4 T at D = 1e2 and "
               "9.4e4 at the printed D <~ 600; with the angle rotations the shot fits up to D = 596 (603 without them). R11b-R14 priced the "
               "lookup alone (7(D-1) per pass, 2-3 passes, D <~ 2e3): lookup_r14.",
               "Insertion-stage register 40 + 7 address + 7 unary-iteration AND workspace + 1 test = 55; load "
               "stage 40 + 47 = 87 (L = 10 enumeration + 19 identification + 18 workspace at D <~ 600); during the "
               "angle step 40 + 10 + 3b - 1 = 88, the peak (M3 ruling). "
               "MLAE reflection stage 58 + 56 = 114 (S_0 on system + 10 enumeration + 7 LCU address + test; the "
               "identification register is clean for any input), plus the 13 reused phase-gradient qubits: 127.",
               "Shots (R16): Hadamard test (lambda'/|<O>|)^2 eps^-2 x 6 = 32.76 x 100 x 6 = 2.0e4 at the upper end of "
               "the any-state band (the 27Al convention; closed shell 19.85 gave 1.2e4 in R15, identity kept 4.6e4 "
               "in R14); depth-capped MLAE at D = 1e2: m_max = 5, 7.9e2 circuits; m_max = 1 above D = 221 (228 without the "
               "load-angle rotations; phase-gradient state reused).",
               "Wall time counts the insertion plus one complete load with its angle rotations at D = 1e2 .. 600.",
               "No Toffoli->T conversion is applied: the box counts 'T/Toffoli' as one hard op."),
    )


def _model_2033(a: Assumptions) -> Result:
    env = a.rfi_t_2033.value
    gap = _band(a.gap_mev)
    nt = n_trotter(a, gap)
    tau = (tau_gs(gap[1]), tau_gs(gap[0]))
    s = survey(a)
    al = s["al_sd"]
    # worked line (app11:90) carries the exact N_Trotter = 1000pi-4000pi and prints "1.1-4.4e8" (R4);
    # the round-first product (3.1e3-1.3e4 -> 1.1-4.6e8) is kept only as the record of what moved
    worked = al["t"]
    worked_round_first = tuple(n * al["n_orb"] ** 4 * c_T(a, n * al["n_orb"] ** 4)
                               for n in _band(a.n_trotter_quoted))
    # R15 ruling 4: the sd insertion is one controlled BE query of a general charge-conserving scalar
    # two-body operator in the sd register (the cover for the six Mu2e categories with two-body currents):
    # selection-rule string count at the low end, dense count at the high end, each in the prep's circuit.
    ins_sd = (al["t_ins"][0], al["t_ins_dense"][1])
    total = (al["t_with_ins"][0], al["t_with_ins_dense"][1])
    lq = (al["qubits"] + a.n_anc.lo, al["qubits"] + a.n_anc.hi)
    # shots / wall time (R11 ruling 2033-shots (A)).  The chapter assigns eps = 0.05 (400 shots/op)
    # to the six Mu2e operators and eps = 0.1 (100/op) to the 5 x 3 0nubb elements (app11:92, 233).
    # The 2033 box is the 27Al deliverable: 6 x 400 = 2400 shots per basis scan, 4800-7200 over
    # the 2-3 basis-cutoff scans; Result.shots spans two scans to three (r17; one to three in R15).  The 0nubb
    # shots run on their own circuits and are not in the box.  The pre-R11 box reading
    # "2-8e3" = (100-400) x 20 is kept as shots_box_reading_pre_r11.
    shots_per_op = (round(1 / a.eps_0nubb.value ** 2), round(1 / a.eps_mu2e.value ** 2))     # 100-400
    n_0nubb_elements = a.n_op_0nubb.value * a.n_pairs_cross_check.value                          # 15
    # eq:Nshot11 carries (lambda'/|<O>|)^2; the 27Al sd dipole factor is derived here (al27_dipole_sd,
    # upper end of the any-state band) and applied to all six categories. R15 ruling 1: identity
    # subtracted classically, 1.45 (5.82 in R14 with the identity in the LCU).
    al_dip = al27_dipole_sd(a)
    al_factor = al_dip["shot_factor"]                                                        # 1.45
    shots_mu2e_bare = shots_per_op[1] * a.n_op_mu2e.value                                    # 2400 (R11 print)
    shots_mu2e = shots_mu2e_bare * al_factor                                                  # 3.5e3
    shots_survey = shots_mu2e_bare + shots_per_op[0] * n_0nubb_elements                       # 3900 (bare)
    shots_scans = (shots_mu2e * a.n_basis_scans.lo, shots_mu2e * a.n_basis_scans.hi)          # 7.0e3-1.0e4
    # r17 (accepted small item): Result.shots and the wall span the 2-3 basis-cutoff scans, as the box's
    # Shots row does; the R15 reading (one scan to three) is kept as shots_one_to_three / wall_time_one_scan_s
    shots_one_to_three = (shots_mu2e, shots_scans[1])                                         # 3.5e3-1.0e4
    shots = shots_scans                                                                       # 7.0e3-1.0e4
    shots_r14 = shots_mu2e_bare * al_dip["shot_factor_band_with_identity"][1]                 # 1.4e4 (R14)
    shots_box = _mul(shots_per_op, a.n_survey_elements.value)                                 # 2e3-8e3 (pre-R11)
    # r17: + the per-shot overhead (1 ms on a 105-440 s shot, invisible here)
    t_shot_s = tuple(x * a.t_gate_s.value + a.shot_overhead_s.value for x in total)
    wall = (shots[0] * t_shot_s[0], shots[1] * t_shot_s[1])                                   # 8.5 d-1.8 mo
    # r23: the wall is serial on ONE machine (shots x per-shot time; no machine count anywhere in this model).
    # Against the 5-year campaign horizon it is 0.47-2.9%, so the campaign fits with no reduction needed.
    horizon_s = a.campaign_horizon_yr.value * YEAR_S                                          # 1.578e8 s
    wall_frac_horizon = (wall[0] / horizon_s, wall[1] / horizon_s)                            # 0.0047-0.029
    wall_one_scan = (shots_one_to_three[0] * t_shot_s[0], shots_one_to_three[1] * t_shot_s[1])  # 4.3 d (R15)
    wall_bare = (shots_mu2e_bare * t_shot_s[0], shots_mu2e_bare * a.n_basis_scans.hi * t_shot_s[1])   # 3.0 d-5.7 wk
    wall_scans = (shots_scans[0] * t_shot_s[0], shots_scans[1] * t_shot_s[1])
    wall_box = (shots_box[0] * t_shot_s[0], shots_box[1] * t_shot_s[1])
    # ---- r25 (rulings R1, R9, R10 and the Ch. 3 rulings of 2026-10-02; shot audit Ch. 3) ----------------
    # Readout. The bare one-body spin-independent (SI) elements are pinned classically to +/-0.6% in any state
    # (al27_dipole_sd any-state band 2.79-2.83), so they validate rather than measure. The spin-dependent (SD)
    # element <O_SD> = 2 <S_p> o_0d = 0.33-0.38 carries the shots. A Hadamard test would need (lambda'/<O_SD>)^2
    # = 79-101 x eps^-2; a direct readout after one single-particle basis rotation needs var / (eps <O_SD>)^2,
    # divided by the heralded projection success p. Per Hamiltonian at eps = 0.05, p = 0.5: 5.5e3-7.1e3 shots.
    o_sd = (2 * a.sd_spin_expectation.hi * al_dip["o_0d"], 2 * a.sd_spin_expectation.lo * al_dip["o_0d"])  # 0.380, 0.335
    sd_hadamard_factor = ((al_dip["lam_lcu"] / o_sd[0]) ** 2, (al_dip["lam_lcu"] / o_sd[1]) ** 2)        # 79-101
    var, p_succ = a.direct_readout_variance.value, a.projection_success.value
    sd_shots_p1 = tuple(var / (a.eps_mu2e.value * o) ** 2 for o in o_sd)                             # 2770-3564
    sd_shots = tuple(x / p_succ for x in sd_shots_p1)                                                 # 5540-7129
    si_pin = (al_dip["exp_band_any_state"][1] - al_dip["exp_band_any_state"][0]) / \
        (al_dip["exp_band_any_state"][1] + al_dip["exp_band_any_state"][0])                          # +/-0.61%
    # T-depth at G concurrent rotations, per band end (factory.json): F* ~ G with the insertion tail
    conc = a.trotter_concurrency.value
    cnt = (al["ins_counts"][0], al["ins_dense_counts"][1])
    gms = (al["gamma_ins"]["gamma"], al["gamma_ins_dense"]["gamma"])
    d33 = tuple(depth_2033(a, c_, g_, conc) for c_, g_ in zip(cnt, gms))                              # G = 12
    d33_g1 = tuple(depth_2033(a, c_, g_, 1) for c_, g_ in zip(cnt, gms))                              # = N_T
    d33_g10 = tuple(depth_2033(a, c_, g_, 10) for c_, g_ in zip(cnt, gms))
    f_star_end = (total[0] / d33[0], total[1] / d33[1])                                               # 11.6-11.8
    # tier 1, first result (R1): one Hamiltonian, SD at 5% (the category-disambiguation requirement), SI checked
    # classically. Tier 2, campaign: the 2-3 basis-cutoff scans (Hamiltonians), same readout.
    n1 = a.first_result_scans.value
    shots_first = (n1 * sd_shots[0], n1 * sd_shots[1])                                               # 5.5e3-7.1e3
    shots_camp = (a.n_basis_scans.lo * sd_shots[0], a.n_basis_scans.hi * sd_shots[1])               # 1.1e4-2.1e4
    wall_first = wall_one_machine(a, shots_first, total, d33)                                         # 6.8-36 d
    wall_camp = wall_one_machine(a, shots_camp, total, d33)                                           # 13.5-109 d
    ex = depth_exports(total, d33, shots_camp, a.t_gate_s, a.shot_overhead_s)
    ex_first = depth_exports(total, d33, shots_first, a.t_gate_s, a.shot_overhead_s)
    floor_g1 = tuple(n * d_ * REACTION_TIME_S for n, d_ in zip(shots_camp, d33_g1))                    # gate by gate
    camp_frac_horizon = (wall_camp[0] / (a.campaign_horizon_yr.value * YEAR_S),
                         wall_camp[1] / (a.campaign_horizon_yr.value * YEAR_S))
    # sensitivity: p over its band (x 0.56 - 5 on the shots), and R1's generic 30% on the SD element
    p_band = _band(a.projection_success_band)
    wall_camp_p_band = (wall_camp[0] * p_succ / p_band[1], wall_camp[1] * p_succ / p_band[0])
    shots_first_30 = tuple(var / (a.eps_first_generic.value * o) ** 2 / p_succ for o in o_sd)          # 154-198
    wall_first_30 = wall_one_machine(a, shots_first_30, total, d33)                                    # 4.5-24 h
    # R7: the simplified 0nubb rows (first light 48Ca in f7/2 p3/2, then f7/2 p3/2 p1/2 by Rabi; full pf Rabi)
    simp = simplified_0nubb(a)
    # pairs with the derived (selection-rule) insertion added, the sd headline's own convention
    jj44_with_ins = s["jj44"]["t_with_ins"]
    near_miss_share = (s["pf"]["insertion_share"][1], s["pf"]["insertion_share"][0])     # R16: pf alone
    valence_pairs_share = (min(s[k]["insertion_share"][1] for k in ("pf", "jj44")),
                           max(s[k]["insertion_share"][0] for k in ("pf", "jj44")))
    n_calls = ((a.n_op_mu2e.value + a.n_op_0nubb.value * a.n_pairs_cross_check.value) * a.n_basis_scans.lo,
               (a.n_op_mu2e.value + a.n_op_0nubb.value * a.n_pairs_cross_check.value) * a.n_basis_scans.hi)
    syn = synthesis_check(a)
    breakdown = (
        Primitive("trotter_rotations_prep_al_sd", al["n_rot"][0], al["ins_counts"][0]["c_T"],
                  CircuitStatus.SCALING, f"{TEX}:82,88,90",
                  "N_proj N_Trotter N_orb^4 dense two-body Pauli-string rotations at the low end "
                  "(Delta=2 MeV, c_proj=pi); RUS c_T at the circuit's R-TOL eps_rot = sqrt(1e-2/N_rot), "
                  f"N_rot = prep + insertion rotations: {al['ins_counts'][0]['eps_rot']:.3e} -> "
                  f"{al['ins_counts'][0]['c_T']:.3f} T (prep alone {al['c_T'][0]:.3f})"),
        Primitive("operator_insertion_sd_rotations", al["ins_counts"][0]["n_rot"], al["ins_counts"][0]["c_T"],
                  CircuitStatus.COMPILED, "derived here",
                  f"one controlled BE query, Gamma = {al['gamma_ins']['gamma']} (selection-rule count): PREPARE + "
                  "PREPARE^dag rotation trees at the circuit's R-TOL tolerance"),
        Primitive("operator_insertion_sd_select", al["ins_counts"][0]["n_toffoli"], a.insertion_t_per_toffoli.value,
                  CircuitStatus.COMPILED, "Babbush_PRX_2018", "controlled unary iteration, Gamma - 1 Toffolis"),
    )
    inter = dict(
        tau_gs_gev_inv=tau, tau_gs_over_ch2_horizon=_div(tau, a.ch2_response_horizon_gev_inv.value),
        n_trotter=nt, n_trotter_quoted=_band(a.n_trotter_quoted),
        system_qubits=al["qubits"], n_anc=_band(a.n_anc), lq_sum=lq,
        register_pf=s["pf"]["qubits"], register_jj44=s["jj44"]["qubits"], register_jj55=s["jj55"]["qubits"],
        register_mo_ru=s["mo_ru"]["qubits"], survey_lq_max=s["jj44"]["lq"][1],
        n_rot_al_sd=al["n_rot"], eps_rot_al_sd=al["eps_rot"], c_T_al_sd=al["c_T"],
        eps_rot_pf=s["pf"]["eps_rot"], c_T_pf=s["pf"]["c_T"],
        eps_rot_jj44=s["jj44"]["eps_rot"], c_T_jj44=s["jj44"]["c_T"],
        c_T_pf_jj44=(s["pf"]["c_T"][0], s["jj44"]["c_T"][1]),
        eps_rot_pf_jj44=(s["pf"]["eps_rot"][0], s["jj44"]["eps_rot"][1]),
        prep_al_sd=al["t"], worked_line=worked, worked_line_round_first=worked_round_first,
        insertion_sd=ins_sd, gamma_ins_sd=(al["gamma_ins"]["gamma"], al["gamma_ins_dense"]["gamma"]),
        insertion_sd_counts=(al["ins_counts"][0], al["ins_dense_counts"][1]), t_per_shot=total,
        fraction_of_envelope_prep=al["multiple"], fraction_of_envelope_total=_div(total, env),
        eps_l=eps_l(a, total), fault_budget=a.fault_budget.value,
        prep_pf=s["pf"]["t"], multiple_pf=s["pf"]["multiple"],
        prep_jj44=s["jj44"]["t"], multiple_jj44=s["jj44"]["multiple"],
        eps_l_pf_jj44=(s["jj44"]["eps_l"][0], s["pf"]["eps_l"][1]),
        gamma_ins_jj44=(s["jj44"]["gamma_ins"]["gamma"], s["jj44"]["gamma_ins_dense"]["gamma"]),
        terms_ins_jj44=(s["jj44"]["gamma_ins"]["n_terms"], s["jj44"]["gamma_ins_dense"]["n_terms"]),
        insertion_jj44=s["jj44"]["t_ins"], insertion_jj44_dense=s["jj44"]["t_ins_dense"],
        insertion_jj44_counts=s["jj44"]["ins_counts"], insertion_jj44_dense_counts=s["jj44"]["ins_dense_counts"],
        insertion_jj44_pre_r15=INSERTION_JJ44_PRE_R15,
        insertion_fraction_of_jj44_prep=(s["jj44"]["insertion_share"][1], s["jj44"]["insertion_share"][0]),
        insertion_fraction_of_pf_prep=(s["pf"]["insertion_share"][1], s["pf"]["insertion_share"][0]),
        insertion_fraction_of_jj55_prep=(s["jj55"]["insertion_share"][1], s["jj55"]["insertion_share"][0]),
        insertion_fraction_of_mo_ru_prep=(s["mo_ru"]["insertion_share"][1], s["mo_ru"]["insertion_share"][0]),
        insertion_fraction_near_miss=near_miss_share, insertion_fraction_pf_jj44=valence_pairs_share,
        insertions={k: dict(gamma=s[k]["gamma_ins"], gamma_dense=s[k]["gamma_ins_dense"], t=s[k]["t_ins"],
                            t_dense=s[k]["t_ins_dense"], share=s[k]["insertion_share"])
                    for k in ("al_sd", "pf", "jj44", "jj55", "mo_ru")},
        prep_jj55=s["jj55"]["t"], multiple_jj55=s["jj55"]["multiple"], eps_l_jj55=s["jj55"]["eps_l"],
        multiple_jj55_large_gap=s["jj55"]["multiple_large_gap"],
        n_orb_mo_ru=s["mo_ru"]["n_orb"], prep_mo_ru=s["mo_ru"]["t"],
        multiple_mo_ru=s["mo_ru"]["multiple"],
        multiple_mo_ru_large_gap=s["mo_ru"]["multiple_large_gap"],
        multiple_mo_ru_small_gap=s["mo_ru"]["multiple_small_gap"],
        eps_l_mo_ru=(eps_l(a, s["mo_ru"]["t_small_gap"])[0], eps_l(a, s["mo_ru"]["t_large_gap"])[1]),
        near_miss_multiple=s["pf"]["multiple"],                      # R16: pf alone (jj44 is co-design)
        prep_plus_insertion_jj44=jj44_with_ins, multiple_jj44_with_insertion=_div(jj44_with_ins, env),
        codesign_multiple=s["jj55"]["multiple"],                     # R15: jj55 alone (Mo/Ru moved below it)
        codesign_multiple_jj44=s["jj44"]["multiple"],                # R16: jj44 relabelled co-design
        register_mo_ru_modes=(2 * a.n_orb_jj44.value, 2 * a.n_orb_jj55.value),
        lq_mo_ru=s["mo_ru"]["lq"], lq_jj55=s["jj55"]["lq"],
        mo_ru_after_lever_generic=_div(s["mo_ru"]["t_after_lever"], env),
        jj44_worst_after_lever=s["jj44"]["multiple_after_lever"][1],
        pf_after_lever=s["pf"]["multiple_after_lever"], jj44_after_lever=s["jj44"]["multiple_after_lever"],
        mo_ru_after_lever=(prep_t(a, 2, s["mo_ru"]["n_orb"], gap_l_band(a), a.prep_lever.value)[0] / env,
                           prep_t(a, 2, s["mo_ru"]["n_orb"], _band(a.gap_small_mev), a.prep_lever.value)[1] / env),
        shots_per_operator=shots_per_op, shots_mu2e_per_scan=shots_mu2e, shots_survey_per_scan=shots_survey,
        shots=shots_camp, shots_pre_r25=shots, shots_with_basis_scans=shots_scans, shots_box_reading_pre_r11=shots_box, shots_r14=shots_r14,
        t_shot_s=t_shot_s, wall_time_s=wall, wall_time_with_basis_scans_s=wall_scans,
        campaign_horizon_s=horizon_s, wall_fraction_of_horizon=camp_frac_horizon,
        wall_fraction_of_horizon_pre_r25=wall_frac_horizon,
        wall_fits_horizon_on_one_machine=wall_camp[1] <= horizon_s,
        wall_time_one_scan_s=wall_one_scan, shots_one_to_three=shots_one_to_three,
        al27_dipole=al_dip, al27_dipole_shot_factor=al_factor, shots_mu2e_bare=shots_mu2e_bare,
        wall_time_bare_s=wall_bare,
        wall_time_box_reading_pre_r11_s=wall_box,
        # R14 option A: the 0nubb readout is the daughter-projection success probability p = |M/lambda|^2;
        # relative eps on |M| is 2 eps on p, so N = (1 - p)/(4 eps^2 p) ~ (lambda/|M|)^2 / (4 eps^2):
        # 25 (lambda/|M|)^2 per element at eps = 0.1. lambda for the 0nubb operators is not derived here.
        shots_0nubb_per_lambda_over_m_sq=1.0 / (4 * a.eps_0nubb.value ** 2),
        n_proj_pairs=2, n_registers_pairs=1,
        subroutine_calls=n_calls, survey_elements_sum=a.n_op_mu2e.value + a.n_op_0nubb.value * a.n_pairs_cross_check.value,
        ext_qubitization_sd_toffoli=a.ext_qubitization_sd_toffoli.value,
        **syn,
        # r25 (R1, R9, R10): the two tiers, T-depth and exports (campaign run), single machine, serial
        o_sd=o_sd, sd_hadamard_factor=sd_hadamard_factor, sd_shots_per_hamiltonian_p1=sd_shots_p1,
        sd_shots_per_hamiltonian=sd_shots, si_pinned_fraction=si_pin, projection_success=p_succ,
        trotter_concurrency=conc, t_depth_gate_by_gate=d33_g1, t_depth_g10=d33_g10, f_star_per_end=f_star_end,
        floor_wall_gate_by_gate_s=floor_g1,
        shots_first_result=shots_first, shots_campaign=shots_camp,
        **{k: ex[k] for k in ("t_depth_per_shot", "f_star", "floor_wall_s", "factories_for_1yr")},
        wall_first_result_s=wall_first, wall_campaign_s=wall_camp,
        depth_exports_campaign=ex, depth_exports_first_result=ex_first,
        campaign_fraction_of_horizon=camp_frac_horizon, wall_campaign_p_band_s=wall_camp_p_band,
        shots_first_result_at_30pct=shots_first_30, wall_first_result_at_30pct_s=wall_first_30,
        simplified_0nubb=simp,
        # referee M3 (2026-10-05): 27Al state prep at the USDB gap, against the generic band (record; band kept)
        al27_prep_t_at_usdb_gap=prep_t(a, s["al_sd"]["n_proj"], s["al_sd"]["n_orb"],
                                       (a.gap_al27_usdb_mev.value, a.gap_al27_usdb_mev.value)),
        # pre-r25 record (six categories at the dipole factor, Hadamard test, 1 ms): wall_time_s_pre_r25
        wall_time_s_pre_r25=wall,
    )
    return Result(
        era="2033", lq=lq, hard_ops=total, breakdown=breakdown, intermediates=inter,
        shots=shots_camp, wall_time_s=wall_camp, epsilon_l=eps_l(a, total),
        notes=("Headline is the 27Al sd cell, prep + insertion; the pairs are carried in intermediates "
               "and INSTANCE_ROWS with fits=False.",
               "epsilon_l = 0.1/N per shot (R3); every T band carries the exact N_Trotter and rounds once "
               "at print (R4): see printed(a).",
               "Breakdown is at the low end of the band (Delta=2 MeV, c_proj=pi, symmetry-restricted insertion).",
               "Box prints the register 74-174 LQ (24 system + 50-150 ancilla; r17, was '<= 250').",
               "r25: shots = the SD element by direct readout at eps 0.05 over p = 0.5, 5.5e3-7.1e3 per "
               "Hamiltonian; Result.shots and wall_time_s are the campaign (2-3 scans, 1.1e4-2.1e4 shots, "
               "13.5 d-3.6 mo, at most 6.0% of the 5-year horizon), serial on one machine at G = 12 (F* 11.8). "
               "First result (one Hamiltonian): 6.8-36 d. Pre-r25 (6 x 400 x 1.45 per scan, 8.5 d-1.8 mo) is "
               "kept as wall_time_s_pre_r25. The 0nubb elements run on their own circuits."),
    )


def _model_codesign(a: Assumptions) -> Result:
    env = a.rfi_t_2033.value
    s = survey(a)
    e_lo, e_hi = a.emax_converged_al.lo, a.emax_converged_al.hi
    lo, hi = s[f"al_emax{e_lo}"], s[f"al_emax{e_hi}"]
    nt = n_trotter(a, _band(a.gap_mev))
    lq = (lo["qubits"], hi["qubits"])                              # system only, as the chapter quotes
    t = (lo["t"][0], hi["t"][1])
    ge_orbs = (ho_spatial_orbitals(a.emax_converged_ge.lo), ho_spatial_orbitals(a.emax_converged_ge.hi))
    ge_qubits = tuple(register(a, 1, n) for n in ge_orbs)       # R14: one register (was two)
    e3, e4 = s[f"al_emax{a.emax_first.lo}"], s[f"al_emax{a.emax_first.hi}"]
    lq_anc = (lo["qubits"] + a.n_anc.lo, hi["qubits"] + a.n_anc.hi)                 # 710-1294
    cap_lq = a.rfi_lq_2033.value
    # first-quantized line of eq:Nq11 (app11:80) for the converged 76Ge pair: 2 A ceil(log2 N_orb)
    # R14 option A: one register, A ceil(log2 N_orb) = 684-760 (was 2 A ceil(log2 N_orb) = 1368-1520)
    fq_ge = tuple(a.mass_number_ge.value * math.ceil(math.log2(n)) for n in ge_orbs)       # 684-760
    fq_compression = (ge_qubits[0] / fq_ge[0], ge_qubits[1] / fq_ge[1])                     # x1.7-3.6
    breakdown = (
        Primitive(f"trotter_rotations_prep_al_emax{e_lo}", lo["n_rot"][0], lo["c_T"][0],
                  CircuitStatus.SCALING, f"{TEX}:152",
                  "converged 27Al, e_max=8 (165 spatial orbitals), low end of the gap band; "
                  "dense N_orb^4 strings per step, RUS c_T at R-TOL eps_rot = "
                  f"{lo['eps_rot'][0]:.3e} -> {lo['c_T'][0]:.3f} T"),
    )
    return Result(
        era="codesign", lq=lq, hard_ops=t, breakdown=breakdown,
        intermediates=dict(
            ho_orbitals={e: ho_spatial_orbitals(e) for e in (3, 4, 8, 10, 14)},
            ho_orbitals_converged=(lo["n_orb"], hi["n_orb"]), ho_orbitals_quoted=_band(a.ho_orbitals_quoted),
            system_qubits=lq, lq_with_ancilla=lq_anc,
            fits_1000_lq=hi["qubits"] <= cap_lq,                      # False: 1144 > 1000
            fits_1000_lq_emax8=lo["qubits"] <= cap_lq,                # True: 660
            fits_1000_lq_emax8_with_ancilla=lq_anc[0] <= cap_lq,      # True: 710
            fits_1000_lq_emax10=hi["qubits"] <= cap_lq,               # False
            t_flagship=t, multiple_flagship=_div(t, env),
            n_rot_flagship=(lo["n_rot"][0], hi["n_rot"][1]),
            eps_rot_flagship=(lo["eps_rot"][0], hi["eps_rot"][1]),
            c_T_flagship=(lo["c_T"][0], hi["c_T"][1]),
            eps_rot_jj55=s["jj55"]["eps_rot"], c_T_jj55=s["jj55"]["c_T"],
            eps_rot_mo_ru=s["mo_ru"]["eps_rot"], c_T_mo_ru=s["mo_ru"]["c_T"],
            t_jj55=s["jj55"]["t"], t_mo_ru=s["mo_ru"]["t"],
            eps_l_flagship=eps_l(a, t),
            register_emax3=e3["qubits"], t_emax3=e3["t"], multiple_emax3=e3["multiple"],
            register_emax4=e4["qubits"], t_emax4=e4["t"], multiple_emax4=e4["multiple"],
            ge_converged_orbitals=ge_orbs, ge_converged_qubits=ge_qubits,
            ge_first_quantized_qubits=fq_ge, ge_first_quantized_compression=fq_compression,
            multiple_jj55=s["jj55"]["multiple"], multiple_jj55_large_gap=s["jj55"]["multiple_large_gap"],
            multiple_mo_ru=s["mo_ru"]["multiple"], multiple_mo_ru_large_gap=s["mo_ru"]["multiple_large_gap"],
            multiple_mo_ru_small_gap=s["mo_ru"]["multiple_small_gap"],
            lq_jj55=s["jj55"]["lq"], lq_mo_ru=s["mo_ru"]["lq"],
        ),
        epsilon_l=eps_l(a, t),
        notes=("Flagship converged-basis 27Al, e_max=8-10, single register; LQ is system qubits only "
               "('~700-1100 system qubits in a single register', app11:135).",
               "R-TOL (2026-09-29): each end is its own circuit. Low end e_max=8, N_rot = 2.33e12, "
               "eps_rot = 6.55e-8, 36.64 T/rot -> 8.53e13; high end e_max=10, N_rot = 8.41e13, "
               "eps_rot = 1.09e-8, 39.62 T/rot -> 3.33e15 (was 6.29e13-2.27e15 at c_T = 27).",
               "app11:152 (R11 (A)): the register sits at the 1000-LQ envelope, inside at e_max=8 (660, 710 "
               "with ancilla) and 14% over at e_max=10 (1144 system, 1294 with ancilla) before ancilla."),
    )


def model(a: Assumptions, era: str) -> Result:
    """Return the Result for '2028' | '2033' | 'codesign'."""
    if era == "2028":
        return _model_2028(a)
    if era == "2033":
        return _model_2033(a)
    if era == "codesign":
        return _model_codesign(a)
    raise ValueError(f"unknown era {era!r}; expected one of {ERAS}")


def utility(a: Assumptions) -> dict:
    """The utility box's checkable arithmetic (app11:188-189).

    G5 open item 4 (Claude's decision, 2026-10-04): the base is the Mu2e TPC alone ($315.7M,
    DOE_HEP_FY25_CJ p. 245), at the stated ~60% (not derived), over the six Mu2e operator categories:
    0.6 x 315.7 = $189.4M ('~$190M'), /6 = $31.57M ('~$32M per instance'). No 0nubb dollar figure
    (LEGEND-1000 / nEXO no longer in the base). Superseded: 0.61 x ($411M + $406M) / 20 = $25M.
    """
    n = a.n_op_mu2e.value                                                             # 6
    total = a.utility_fraction.value * a.mu2e_tpc_musd.value                          # 189.42
    return dict(exposure_equivalent=a.nme_spread.value ** 4,                          # 3^4 = 81, '~80' (0nubb, no $)
                base_musd=a.mu2e_tpc_musd.value,                                      # 315.7
                fraction=a.utility_fraction.value,                                    # 0.6
                elements=n,
                musd_attributed_total=total,                                          # 189.42 -> '~$190M'
                musd_per_instance=total / n)                                          # 31.57 -> '~$32M'


PUBLISHED = {
    "2028": Published(lq=(120, 170), hard_ops=(1.7e4, 1.0e5),
                      src="app11 2028 box ('~120-170', '1.7e4-1.0e5 T per shot: insertion + complete load with its "
                          "angle rotations at D = 1e2-600'; M3 ruling 2026-10-04: angle rotations folded in, model "
                          "17020.3-100516.3; referee G2 1.6e4-9.9e4 without them; the insertion alone 5.4e3 (R13) "
                          "was the headline before)",
                      rel_tol=0.05),
    "2033": Published(lq=(74, 174), hard_ops=(1.1e8, 4.4e8),
                      src="app11:2033 box lines 208-209 ('<= 250 (24 system + 50-150 ancilla)', "
                          "'1.1e8-4.4e8 T, prep + insertion'); LQ is the printed component sum, "
                          "'<= 250' being the envelope. R15: model 1.054e8-4.404e8 with the derived one-query "
                          "insertion (1.092e8-4.771e8 in R14 with the scaled stated insertion)",
                      rel_tol=0.10),
    "codesign": Published(lq=(700, 1100), hard_ops=(8.5e13, 3.3e15),
                          src="app11:2033 box lines 214-215 and summary line 239 ('~700-1100 LQ', "
                              "'8.5e13-3.3e15 hard ops (co-design)'; R-TOL 2026-09-29, was 6e13-2.3e15 "
                              "at c_T = 27); model 660-1144 LQ, 8.53e13-3.33e15",
                          rel_tol=0.10),
}


def printed(a: Assumptions) -> dict:
    """Every printed string the chapter derives from eq:Ngate11, rounded once at print (R4).

    Keys name the site; values are the exact substrings in app11_mu2e_0nubb.tex. T bands use the
    factored two-figure form the boxes use ('$1.7$--$6.8\\times 10^{9}$'); multiples and eps_l
    are two-figure numbers; the flagship eps_l is one figure, as the chapter prints it.
    """
    from estimates.common import tex_range_factored
    s = survey(a)
    r33 = model(a, "2033").intermediates
    rcd = model(a, "codesign").intermediates

    def mant(x):
        # two significant figures, mantissa always with one decimal ("1.0", "9.0"), as the chapter prints
        e = math.floor(math.log10(abs(x)))
        m = round(x / 10 ** e, 1)
        if m >= 10:
            m, e = m / 10, e + 1
        return f"{m:.1f}", e

    def fmt_t(r):
        (ml, el), (mh, eh) = mant(r[0]), mant(r[1])
        if el == eh:
            return f"${ml}$--${mh}\\times 10^{{{el}}}$"
        return f"${ml}\\times 10^{{{el}}}$--${mh}\\times 10^{{{eh}}}$"

    return {
        "worked_line:90": fmt_t(r33["worked_line"]),                                  # $1.1$--$4.4\times 10^{8}$
        "prep_pf:139": fmt_t(s["pf"]["t"]),                                           # $1.7$--$6.8\times 10^{9}$
        "prep_jj44:139": fmt_t(s["jj44"]["t"]),                                       # $2.5$--$9.9\times 10^{9}$
        "multiple_pf:147": print_range(s["pf"]["multiple"]),                          # 1.7--6.8
        "multiple_jj44:147": print_range(s["jj44"]["multiple"]),                      # 2.5--9.9
        "near_miss:45,217": print_range(r33["near_miss_multiple"]),                     # 1.8--7.3 (R16: pf alone)
        "jj55_large_gap:150": print_range(s["jj55"]["multiple_large_gap"]),           # 5.6--11
        "t_emax3:152": fmt_t(rcd["t_emax3"]),                                         # $1.4$--$5.4\times 10^{10}$
        "multiple_emax3:152,156": print_range(rcd["multiple_emax3"]),                 # 14--54
        "multiple_emax4:152,156": print_range(rcd["multiple_emax4"]),                 # 130--510
        "mo_ru_large_gap:162": print_range(s["mo_ru"]["multiple_large_gap"]),         # 3.0--6.0 (R15)
        "mo_ru_small_gap:162,219": print_range(s["mo_ru"]["multiple_small_gap"]),     # 25--51 (R15)
        # overshoot as a percent, one figure (R4): 1.0109 -> "by about 1\\%" (R-TOL; was 0.99, under)
        "lever_worst:168": f"by about {round_sig(100 * (r33['jj44_worst_after_lever'] - 1), 1):g}\\%",
        "pf_after_lever:168": print_range(r33["pf_after_lever"]),                     # 0.16--0.68 (R16)
        "mo_ru_after_lever:169": print_range(r33["mo_ru_after_lever"]),               # 0.28--4.8 (R15)
        "mo_ru_after_lever_generic:169": print_range(r33["mo_ru_after_lever_generic"]),  # 0.56--2.4
        "n_rot_flagship:152": fmt_t(rcd["n_rot_flagship"]),                           # $2.3\times 10^{12}$--$8.4\times 10^{13}$
        "t_flagship:152,239": fmt_t(rcd["t_flagship"]),                               # $8.5\times 10^{13}$--$3.3\times 10^{15}$
        "multiple_flagship:152,215": fmt_t(rcd["multiple_flagship"]),                 # $8.5\times 10^{4}$--$3.3\times 10^{6}$
        "eps_l_sd:205,240": fmt_t(r33["eps_l"]),                                      # $2.1$--$8.7\times 10^{-10}$
        "eps_l_pf_jj44:240": fmt_t(r33["eps_l_pf_jj44"]),                             # $1.0$--$5.9\times 10^{-11}$
        "eps_l_jj55:240": fmt_t(r33["eps_l_jj55"]),                                   # $2.2$--$9.0\times 10^{-12}$
        "eps_l_flagship:240": tex_range_factored(*rcd["eps_l_flagship"], 0),          # $3\times 10^{-17}$--$1\times 10^{-15}$
        "prep_al_sd:139": fmt_t(s["al_sd"]["t"]),                                     # $1.0$--$4.4\times 10^{8}$
        "fraction_total_sd:147,210": print_range(r33["fraction_of_envelope_total"]),  # 0.11--0.44
        "multiple_jj55:150": print_range(s["jj55"]["multiple"]),                      # 12--50
        "mo_ru_generic:46,162,219": print_range(s["mo_ru"]["multiple"]),              # 6.0--25 (R15)
        "codesign_multiple:218": print_range(r33["codesign_multiple"]),               # 12--50 (jj55 alone, R15)
        "codesign_jj44:45,152,218,244": print_range(r33["codesign_multiple_jj44"]),   # 2.6--11 (R16: jj44 co-design)
        # R-TOL (2026-09-29): the tolerance each circuit sets from its own rotation count
        "eps_rot_sd:90": fmt_t(r33["eps_rot_al_sd"][::-1]),                           # $2.5$--$5.0\times 10^{-5}$
        "c_T_sd:90": f"{r33['c_T_al_sd'][0]:.1f}$--${r33['c_T_al_sd'][1]:.1f}",       # 25.6$--$26.8
        "eps_rot_survey:90": fmt_t(r33["eps_rot_survey"][::-1]),                      # $8.7\times 10^{-7}$--$5.0\times 10^{-5}$
        "c_T_survey:90": print_range(r33["c_T_survey_range"]),                        # 26--32
        "c_T_coherent:90": print_range(r33["c_T_coherent_range"]),                    # 42--56
        "eps_rot_flagship:152": fmt_t(rcd["eps_rot_flagship"][::-1]),                 # $1.1$--$6.6\times 10^{-8}$
        "c_T_flagship:152": print_range(rcd["c_T_flagship"]),                         # 37--40
        # R15 ruling 4: the derived insertion as a share of the pair preparation, one figure
        "insertion_share_jj44:143": "{:g}$--${:g}\\%".format(
            round_sig(100 * r33["insertion_fraction_of_jj44_prep"][0], 1),
            round_sig(100 * r33["insertion_fraction_of_jj44_prep"][1], 1)),               # 0.04$--$0.2\%
        "insertion_share_near_miss:217": "{:g}$--${:g}\\%".format(
            round_sig(100 * r33["insertion_fraction_near_miss"][0], 1),
            round_sig(100 * r33["insertion_fraction_near_miss"][1], 1)),                  # 0.06$--$0.2\% (R16: pf)
        "insertion_jj44:143": fmt_t(r33["insertion_jj44"]),                           # $4.0$--$4.1\times 10^{6}$
        "insertion_jj44_dense:143": fmt_t(r33["insertion_jj44_dense"]),               # $6.5$--$6.8\times 10^{7}$
        "insertion_sd:143": fmt_t(r33["insertion_sd"]),                               # $9.3\times 10^{5}$--$3.9\times 10^{6}$
        "t_per_shot_sd:209,243": "${}\\times 10^{{{}}}$--${}\\times 10^{{{}}}$".format(
            *mant(r33["t_per_shot"][0]), *mant(r33["t_per_shot"][1])),                  # box spelling, both exponents
        "eps_l_mo_ru_row:246": fmt_t(r33["eps_l_mo_ru"]),
        # r23: the single-machine serial campaign as a share of the 5-year horizon, upper end, two figures
        "wall_fraction_of_horizon_max:149": "at most ${:.1f}\\%$".format(
            round_sig(100 * r33["wall_fraction_of_horizon"][1])),                       # at most $6.0\%$ (r25)
        **printed_r25(a),
    }


def printed_r25(a: Assumptions) -> dict:
    """r25 printed strings (R1, R7, R9, R10): walls on one machine, T-depth, the two 2033 tiers and the
    0nubb first-light rows. Rounded once at print, two figures unless the chapter prints fewer."""
    i28 = model(a, "2028").intermediates
    i33 = model(a, "2033").intermediates
    sim = i33["simplified_0nubb"]
    day, month, yr = 86400.0, YEAR_S / 12, YEAR_S

    def g2(x):
        v = round_sig(x, 2)
        if v >= 10:
            return f"{int(round(v))}"
        e = math.floor(math.log10(abs(v)))
        return f"{v:.{max(1 - e, 0)}f}"

    def sci(x):
        e = math.floor(math.log10(round_sig(x, 2)))
        return f"${g2(x / 10 ** e)}\\times 10^{{{e}}}$"

    def fact(r, e=None):
        e = math.floor(math.log10(r[1])) if e is None else e
        return f"${g2(r[0] / 10 ** e)}$--${g2(r[1] / 10 ** e)}\\times 10^{{{e}}}$"

    w28, wae = i28["wall_campaign_s"], i28["mlae_wall_time_corrected_s"]
    tsh = (w28[0] / i28["shots"], w28[1] / i28["shots"])
    w1, wc = i33["wall_first_result_s"], i33["wall_campaign_s"]
    fl = i33["floor_wall_gate_by_gate_s"]
    pb = i33["wall_campaign_p_band_s"]
    w30 = i33["wall_first_result_at_30pct_s"]
    b, bp, h, apf = sim["f7p3_projection"], sim["f7p3p1_rabi"], sim["pf_rabi"], sim["pf_projection_0plus"]
    return {
        # 2028 (app11:127 and the 2028 box)
        "t_shot_supply_2028:127": f"${i28['t_shot_s'][0]:.3f}$--${round_sig(i28['t_shot_s'][1], 2):.2f}$\\,s",
        "t_layers_2028:127": f"${g2(i28['t_depth_per_shot'][0] / 1e3)}\\times 10^{{3}}$ sequential T layers against "
                             f"${g2(i28['t_per_shot'][0] / 1e4)}\\times 10^{{4}}$ T gates",
        "f_star_d100_2028:127": f"at most ${g2(i28['f_star_per_end'][0])}$ factories",
        "t_shot_corrected_d100_2028:127": f"the shot takes ${tsh[0]:.3f}$\\,s",
        "f_star_d600_2028:127": f"At $D=600$ the ratio is ${g2(i28['f_star_per_end'][1])}$",
        "mlae_f_star_2028:127": f"ratio of ${g2(i28['mlae_f_star'])}$ and takes "
                                f"${round_sig(i28['mlae_depth']['t_deepest'] * a.t_gate_s.value + a.shot_overhead_s.value, 2):.3f}$\\,s",
        "walls_2028:127": f"take ${round(w28[0] / 60)}$--${round(w28[1] / 60)}$ min and ${round(wae)}$ s on one machine",
        "box_t_shot_2028:187": f"$\\sim {tsh[0]:.3f}$--${round_sig(tsh[1], 2):.2f}\\,$s per Hadamard-test shot",
        "box_wall_2028:189": f"$\\sim {round(w28[0] / 60)}$--${round(w28[1] / 60)}$ min (Hadamard test) or "
                             f"$\\sim {round(wae)}$ s (MLAE",
        # 2033 deliverable (app11:149-151 and the box)
        "o_sd:149": f"= {round_sig(i33['o_sd'][1], 2):g}$--${round_sig(i33['o_sd'][0], 2):g}$",
        "sd_factor:149": f"= {round(i33['sd_hadamard_factor'][0])}$--${round(i33['sd_hadamard_factor'][1])}$",
        "sd_shots_p1:149": "= " + fact(i33["sd_shots_per_hamiltonian_p1"], 3)[1:] + " shots",
        "sd_shots:149,box": fact(i33["sd_shots_per_hamiltonian"], 3),
        "t_depth_2033:151": f"${g2(i33['t_depth_per_shot'][0] / 1e6)}\\times 10^{{6}}$--"
                            f"${g2(i33['t_depth_per_shot'][1] / 1e7)}\\times 10^{{7}}$",
        "f_star_2033:151": f"up to ${round(i33['f_star_per_end'][0])}$ factories stay busy",
        "floor_g1_2033:151": f"below ${round(fl[0] / day)}$ days to ${round_sig(fl[1] / yr, 2):.1f}$ years",
        "first_result_2033:151,box": f"${g2(w1[0] / day)}$--${g2(w1[1] / day)}$ days",
        "shots_campaign_2033:151,box": fact(i33["shots_campaign"], 4),
        "campaign_2033:151": f"${g2(wc[0] / day)}$ days to ${round_sig(wc[1] / month, 2):.1f}$ months",
        "campaign_box_2033": f"$\\approx {g2(wc[0] / day)}$ days--${round_sig(wc[1] / month, 2):.1f}$ months on one machine",
        "p_band_2033:151": f"${round_sig(pb[0] / day, 2):.1f}$ days and ${round_sig(pb[1] / yr, 2):.1f}$ years",
        "first_30pct_2033:151": f"${round_sig(w30[0] / 3600, 2):.1f}$--${g2(w30[1] / 3600)}$ hours",
        # R7 first light (app11:156-158 and the box)
        "pf_0plus_t:156": f"${round_sig(apf['t'][0] / 1e9, 2):.2f}$--${round_sig(apf['t'][1] / 1e9, 2):.1f}\\times 10^{{9}}$ T",
        "pf_0plus_shots:156": fact(apf["shots_lower_bound"], 6) + " projection shots",
        "pf_0plus_years:156": f"${g2(apf['wall_s'][0] / yr)}$--${g2(apf['wall_s'][1] / yr)}$ years",
        "f7p3_t:156,box": fact(b["t"], 8),
        "f7p3_shots:156": f"${g2(b['shots'][0] / 1e4)}\\times 10^{{4}}$ shots",
        "f7p3_wall:156,box": f"${g2(b['wall_s'][0] / day)}$--${g2(b['wall_s'][1] / day)}$ days",
        "f7p3p1_t:157,box": fact(bp["t"], 8),
        "f7p3p1_wall:157,box": f"${g2(bp['wall_s'][0] / day)}$--${g2(bp['wall_s'][1] / day)}$ days",
        "qdrift_f7p3p1:157": f"${g2(bp['n_qdrift'][0] / 1e6)}\\times 10^{{6}}$ at",
        "pf_rabi_t:157": fact(h["t"], 9),
        "pf_rabi_multiple:157,box": f"$\\times {g2(h['multiple'][0])}$--${g2(h['multiple'][1])}$",
        "pf_rabi_wall:157,box": f"${g2(h['wall_s'][0] / day)}$--${g2(h['wall_s'][1] / day)}$ days",
        "rabi_factories:163": f"keep only ${round(min(min(bp['f_star_per_end']), min(h['f_star_per_end'])))}$--"
                              f"${round(max(max(bp['f_star_per_end']), max(h['f_star_per_end'])))}$ factories busy",
        # referee M3: the f7p3 end moved to 1.1e-9 (KB3G gaps), so each end carries its own exponent
        "eps_l_first_light:254": f"{sci(bp['eps_l'][0])}--{sci(b['eps_l'][1])} ($0\\nu\\beta\\beta$ validation stage)",
    }


def INSTANCE_ROWS(a, era, r):
    """Rows for resources.json: one per costed instance the boxes quote."""
    env = a.rfi_t_2033.value
    if era == "2028":
        return [(r"$^4$He dipole Hadamard", r.lq, r.hard_ops,
                 {"note": "insertion + complete state load with its angle rotations (D = 1e2-600; 1.005x at D = 600); "
                          "conditional on D"})]
    s = survey(a)
    if era == "2033":
        rows = [(r"$^{27}$Al $sd$ (prep+insertion)", r.lq, r.hard_ops, {"fits": True})]
        c = s["pf"]
        rows.append((c["label"], c["lq"], c["t"], {"fits": False, "multiple": c["multiple"],
                                                    "note": "near miss; state prep only"}))
        return rows
    if era == "codesign":
        rows = []
        notes = {
            # R16: jj44 relabelled co-design (prep alone x10.8 at the top of its band)
            "jj44": "co-design (top of band past ~10x); state prep only",
            "jj55": "co-design; state prep only, generic gap band",
            "mo_ru": "reach target; co-design at the generic gap band, near miss at a large gap "
                     "(x3-6), x25-51 small gap; state prep only, generic gap band",
        }
        for k in ("jj44", "jj55", "mo_ru"):
            c = s[k]
            note = notes[k]
            rows.append((c["label"], c["lq"], c["t"], {"fits": False, "multiple": c["multiple"],
                                                        "multiple_large_gap": c["multiple_large_gap"],
                                                        "note": note}))
        rows.append((r"flagship $^{27}$Al ($e_{\max}$8--10)", r.lq, r.hard_ops,
                     {"fits": False, "multiple": tuple(x / env for x in r.hard_ops),
                      "note": "system qubits only; per state preparation"}))
        return rows
    return []
