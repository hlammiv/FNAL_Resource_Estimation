"""Ch. 9 — Chiral gauge theory (domain-wall / overlap fermions).
Reproduces the LQ and hard-op numbers of applications/app06_chiral_gauge.tex.

INPUTS AND SOURCES
  L, L5, D, Nc, Nf, |G|, N_rho     setup parameters the chapter declares (Assumed, box lines)
  link qubits                      GROUPS["Z3"].link_qubits = 2, GROUPS["2O"].link_qubits = 6 (R1: compiled = chapter)
  ancilla 32 / 38                  DERIVED HERE r20-r22 (peak of the schedule); the asserted 40 / 110 are records
                                   (ancilla_2028_asserted, ancilla_2033_asserted); 39 / 33 were the r21 counts
  Pauli strings per hop = 2        Stated app06 (the link-free fifth-direction hops)
  fifth-direction hops             DERIVED HERE r21 (a): open boundary, V (L5 - 1) (fifth_boundary)
  spatial hop                      R-HOP, groups.rhop_link (DERIVED HERE, Ch. 6 round D; ruled report-wide):
                                   Z3 diagonal -> no color move; 8 strings per copy; k = 12 passed by Ch. 9
  hop rotation price               full fit (r17 e; R-HOP's slope-only price kept as the record hop_synthesis_legacy)
  on-site term 2028                arXiv:2505.20419 (m + r(1+d)); 2N R_Z per site, FermionPrimitives_unpub section_mass.tex:12
  Pauli strings per electric = 3   Stated app06 ('three diagonal Z-strings per two-qubit link')
  grouping                         Stated app06 ('color+s': 2 x k=72, 64 x k=12, 3 x k=8)
  HWP chunk size k <= 32           Stated app06
  HWP cost                         Cited arxiv_1902_10673 (floor(log2 k)+1 rotations, k-HW(k) Toffolis, k-HW(k) ancilla)
  RUS ancilla = 1                  Cited bocharovRoettelerSvore2015 (common.RUS_ANCILLA; ruling E26)
  Toffoli = 7 T                    ruling R5; main-overview.tex conventions paragraph ("textbook")
  eps_syn = 1e-2, eps_rot per circuit  ruling R-TOL (common.EPS_SYN, eps_rot_for); was t_per_rot = 30 at 1e-4
  rotation price                   Cited BRS/Campbell RUS fit 1.15 log2(1/eps) + 9.2
  N_Trotter 2028 = 20, +1 probe    Stated app06 (ruling R6; was 30)
  faults per shot = 0.1            ruling R3 (TRACKED_CHANGES.md 2026-09-28): eps_l = 0.1 / T_shot
  shot time                        1 us per T + 0.1 ms per shot (report rule; shot_overhead_s)
  2028 shots                       DERIVED HERE r21 (d) from the chapter's sigma^2 ~ 10 at V = 27
  1+1D overlap query               DERIVED HERE r22: 40 Toffolis + 6 rotations (single-particle architecture), lift 190
                                   Toffolis; the r21 (e) second-quantized pricing (240 terms + 25 rotations) is a record
  kappa=30 (r24), delta_sgn=1e-3   Assumed app06 (budget row prints 0.1%); kappa = alpha/Delta since r22
  overlap architecture (2033)      Cited arxiv_2607_28524 (single-particle kernel on an index register, lifted to Fock
                                   space); the circuit and its counts DERIVED HERE r22 (single_particle_query_2O,
                                   colour_multiplexer_index_2O, lift_cost); unary iteration Babbush_PRX_2018; 2O encoding
                                   and U_x arXiv:2312.10285
  Higgs vertex                     40 SELECT strings per site (DERIVED HERE r21 c) + one U_x per site, once per
                                   application, outside sgn (r22)
  U_x = 392 T = 56 Toffoli x 7 T   GROUPS["2O"].primitives["U_mul"] (arXiv:2312.10285 tab:tgatecost); not retyped here
  U_x per Higgs site: 1            Stated app06 (ruling VERTEX, Higgs only since r17)
  N_Trotter 2033                   state-dependent second-order estimate (ruling E26 (ii); Alves, Lamm, Liu,
                                   preprint FERMILAB-PUB-26-0397-T); the weak-coupling evaluation is DERIVED HERE (156 oscillators,
                                   T a = 0.5 Assumed). Worst case: commutator bound (childs2021theory, r21 b), g^2 = 1
  applications per step            DERIVED HERE r22: Jacobi-Anger order at 2Q dt, QSP error 1e-2 over the shot (Assumed)
  beta||H|| = 2Q / (T a) = 864     DERIVED HERE (referee C1, 2026-10-04; the r25 Stated 1e2 is a record, beta_H_rule)
  8 AA rounds / 17 filter calls     Stated app06: a PLACEHOLDER (chapter-open pass 2026-10-05): 17 calls = 8 rounds of
                                   fixed-point AA or 17 post-selection attempts at success 1/17. k rounds call the filter
                                   2k + 1 = 17 times (verifier 2026-10-04; r25 and the first referee pass charged one call
                                   per round, 8, a record: n_dov_tpq_one_call_per_round)
  typical-reference TPQ (C1)       DERIVED HERE (chapter-open pass): free fermion sector of 3^3 at beta = 2a,
                                   p = prod (1 + e^{-beta|eps|})/2 = e^{-137.9}, (pi/4) p^{-1/2} = 7.0e29 rounds; filter
                                   normalized at -2Q loses e^{-beta(E0 + 2Q)} = e^{-490.7} (free_overlap_modes,
                                   typical_reference_tpq). Proxy numerics (bias of amplified basis/product references,
                                   IS cost; verifier 2026-10-05: a Gibbs-marginal product reference amplifies with bias
                                   5e-4 (Ising) and 0.04 (t-V) at n = 8, smaller than a uniformly random basis state's
                                   0.06 and 0.15, but not controlled):
                                   scratch chopen/app06_chiral_gauge/proxy_tpq.py and ch09_check/proxy_check.py, not in
                                   the model
  Delta N_CS readout (C3)          DERIVED HERE (chapter-open pass): Kubo form, gauge-invariant E^a B^a on 2O, Hadamard-test
                                   lambda >= 729 g^2/8pi^2 = 9.23/a (sum of ||E^a|| ||B^a||, a lower bound on an LCU
                                   1-norm), sigma^2 = (4/3)(lambda t)^4/<dN^2>^2. Verifier fix (2026-10-05): E^a_i(x) is the
                                   average of the generators at x on links (x,i) and (x-i,i); the forward-link pairing is
                                   not a total derivative in the free theory and grows as t^2. Free-field background of
                                   the link-averaged clover estimator 9.58e-4 at 10a (ncs_lattice_background), sigma^2 =
                                   1.06e14; forward link 2.56e-3 (1.5e13); helicity-A.B ideal 8.93e-3 (1.2e12,
                                   ncs_free_field_background, a record). sigma2_ncs = 1e2 stays as the Stated
                                   ideal-readout target the shot counts use
  units (C3)                       Cited Moore_Rummukainen_2000 (a = 1/(2 g^2 T) volume table) and
                                   DOnofrio_Rummukainen_Tranberg_2014 (T_c = 159 GeV); beta_L = 8, L g^2 T = 1.5 DERIVED HERE
  delta_beta = 1e-4                Assumed app06; sqrt(864 ln 1e4) = 89.2 -> 89 (r25: sqrt(1e2 ln 1e4) = 30.35 -> 30)
  Haar-start TPQ: 5% deficit, 1e3 q  RECORD (r25 prose '~3e7', retired 2026-10-05 for typical_reference_tpq)
  Gibbs sampler: 300 terms, 10 sweeps Stated app06; per-jump cost ESTIMATED HERE r20 (gibbs_* inputs Assumed),
                                   re-priced at the single-particle application (r22); the old 1e10-1e12 T is a retired record
  RFI envelopes 1e5 / 1e9          Cited DOE_RFI_2026
WHAT IS NOT DERIVED HERE
  the use of ONE U_x as the price of a controlled 2O Higgs action   ruling VERTEX; arXiv:2312.10285 compiles U_x
                                          but no controlled or multiplexed variant, so the line is SCALING
  the Higgs vertex's reading              r21 (c): the on-site bilinear dressed by the Higgs group element and weighted
                                          by rho (5-term binary expansion), once per application. r25 (a): the Yukawa is
                                          Ginsparg-Wilson projected, so its projector (1 - sgn(H_W))/2 costs a second
                                          degree-d_sgn QSVT sequence per application (priced); its effect on the
                                          block-encoding normalization 2Q is NOT priced (NEEDS_AUTHOR)
  the replacement of norms by state expectation values   Alves, Lamm, Liu (preprint FERMILAB-PUB-26-0397-T): cited, not re-derived; the
                                          constants were checked numerically (tests)
  the temperature T a = 0.5               Assumed; the chapter says T ~ T_c only
  a self-inverse block encoding           assumed for 'applications per step = Jacobi-Anger order' (the lift admits one at
                                          O(1) cost); otherwise twice the count (record t_per_shot_if_not_self_inverse)
  the gauge and Higgs terms of a step     NOT PRICED; the gauge part is estimated as a record (about 2% of the step)
  the amplification reflections           NOT PRICED (record: about 2e5 T per shot, arxiv_2407_17966)
  k - HW(k) Toffolis per HWP group        the count of arXiv:1709.06648 / 1902.10673, whose adders uncompute
                                          by measurement; priced here at 7 T per Toffoli as ruled (R5)
  diabatic error of the 20-step ramp      not estimated anywhere; the chapter flags it as the quantity
                                          the 2028 benchmark must measure
  8 amplification rounds                  Stated; chapter calls it the least-controlled input
  the bare coupling g^2                   Assumed 1 for the worst-case bound only; the chapter states none
  ~31 T coherent-synthesis figure         the RUS fit at eps_syn / N_rot = 1.64e-6 gives 31.30 T (reproduced)
  sigma^2_NCS ~ 1e2 at 10a               asserted (app06), an assumed readout variance; the 5a ratio 4.95 is DERIVED
                                          from the free-field signals at fixed absolute noise (2026-10-05); NEEDS_AUTHOR
  the T-depth layer counts                r25 factory analysis (SP_QUERY_DEPTH_2O, application_depth_2O,
                                          hwp_step_depth_2028); the HWP adder depth ceil(2 log2 k) is a carry-save
                                          estimate, the U_x internal depth is bounded by its Toffoli count
WALL TIME AND FACTORIES
  t_gate_s = 1e-6 s per T, shot_overhead_s = 1e-4 s per shot, ONE machine, every shot serial (r23).
  T-depth per shot from the r25 factory analysis (factory.json, Ch. 9 runs), low/high: a Toffoli or AND is one reaction
  layer of t_r = 10 us, a RUS rotation is T-depth = T-count. F* = N_T / D_T; each era's box schedule has F* >= 10, so
  the 10-factory baseline (1 us per T) holds:
    2028: three spatial-hop groups per slot, weight bits in parallel (42 ancilla, 250 LQ): F* 10.4-13.0 (r22-r24
          register, 240 LQ: 1.8; one group per slot with parallel weight bits: 5.3-6.1).
    2033: CNOT fan-out of the shared control in the link read (+5 LQ, 1015): F* 12.2-12.5 on every run (as itemized,
          1010 LQ: 5.6-6.4, walls 1.6-1.8x longer).
  First result: 2028 the condensate at L5 = 4 at 30% (125 shots); 2033 the method test (author ruling 2026-10-05):
  the early growth of the full <dN_CS^2> at one symmetric-side temperature at 5a (5497 shots) and 10a (1111 shots),
  30% each, no rate -> wall_first_result_s. Campaign: 2028 L5 = 2, 3, 4 at 10% (1125 shots each); 2033 the method test
  + 1e3 condensate shots at each of 4-8 temperatures (10%) -> wall_campaign_s (2.7-3.0 yr; the r25 rate campaign,
  36-60 yr, is a record). Exports: t_per_shot / t_depth_per_shot / f_star of the headline shot,
  floor_wall_s and factories_for_1yr summed over the campaign runs (common.depth_exports per run).

r23 (ruling H. Lamm 2026-10-02, "we don't want to anywhere assume we have multiple machines"; apply_log/r23_ch09.md):
one machine, every wall time serial at 1 us per T + ~0.1 ms per shot. n_machines (6.4e2), n_machines_needed (636),
wall_yr_parallel and fits_campaign_horizon are retired. In their place, per observable arm at the box's shot
(2.008e5 s): condensate, 1e3 shots, 6.36 yr per coupling point; sphaleron rate, 1e4 shots, 63.6 yr per point; the
50-point scan 3181 yr. None fits the 5-year horizon: it holds 786 shots (the condensate at one coupling to 11%); one
sphaleron point needs a 12.7x cost reduction and the scan 636x (a 3.156e8-T shot at the same shot count). Levers,
emitted as horizon_levers and not applied: low end of the step range 1.775, kappa halved 1.980, evolution halved 2.440,
link-select lines held through the copy 1.141 (78 Toffolis per query, 75 more qubits, estimated here); all four 8.93
(condensate point 0.71 yr, sphaleron point 7.1 yr). Both arms are priced at the one shot the chapter prices; the
condensate arm's own preparation is not priced separately. No per-shot T or LQ number moves.

r24 (H. Lamm 2026-10-02, after asking why the wall time is 64 yr per point and whether it can get under a year: "I think
it's ok to run for a shorter amount of time, and lower kappa is ok"; apply_log/r24_ch09.md):
  kappa 1e2 -> 30 (d_sgn 691 -> 208, same delta_sgn 1e-3); evolution 10 fm = 50a -> 2 fm = 10a.
  Sphaleron arm: Lambda_sd unchanged (20.25), N = ceil(sqrt(Lambda 10^3 / 0.1)) = 450 (range 107-1274; N ~ t^{3/2});
      x = 2Q dt = 9.6, 19 applications per step; 450 x 19 + 240 = 8790 applications (floor 2Q t = 4320);
      R-TOL 1261 rotations per application x 8790 = 11,084,190 -> eps 3.004e-5, 26.476 T; c_sp 4368.9;
      application 209 x c_sp + lift 6204.9 + reflection 152.5 + Higgs 18,144 = 937,592.7 T; shot 8.2414e9 T
      (8.2x the 1e9 limit), eps_l 1.213e-11, 8241 s per shot, 2.611 yr per point (1e4 shots), 130.6 yr for 50 points.
  Condensate arm (static observable, no real-time evolution): the TPQ preparation alone, 240 applications at its own
      R-TOL (302,640 rotations, eps 1.818e-4, 23.489 T): 2.2412e8 T, eps_l 4.46e-10, 224 s per shot, 2.59 days per
      coupling point (1e3 shots), 0.355 yr for 50 points.
  Levers for one sphaleron point under a year (not applied): low end of the step range 1.46x (1.79 yr), a 15%
      statistical target 2.25x (1.16 yr), both 3.29x (0.79 yr); spectral amplification on the 2Q t floor not priced.
      The r23 levers (kappa halved, evolution halved) are now applied in part; the hold lever is a record.
  Records: kappa_2033_r22 = 1e2 and t_evol_fm_r22 = 10 keep the r22 headline (2.00789e11), the worst-case bound at
      50a (83,195) and every second-quantized record reproducible. PUBLISHED 2033 hard_ops 2.0e11 -> 8.2e9.

r25 (H. Lamm 2026-10-02: rulings R1, R4, R9, R10 and the queued Ch. 9 rulings (a)-(f); apply_log/r25_ch09.md):
  (a) Yukawa inside (H. Singh via H. Lamm): a second degree-208 sequence per application, 417 queries (r24 209);
      R-TOL 2509 rotations per application x 8790 = 22,054,110 -> eps 2.129e-5, 27.047 T; c_sp 4372.3; application
      1,847,747.1 T (1.97x); sphaleron shot (10a) 1.6242e10 T (16.2x the limit), eps_l 6.16e-12; condensate shot
      4.4166e8 T, eps_l 2.26e-10. Records: r24 8.2414e9 / 2.2412e8 (yukawa_placement = 'outside').
  (R4) the 5a shot at its own depth: 160 steps x 23 + 240 = 3920 applications, 7.2366e9 T.
  (R10) 2033 fan-out +5 LQ: 1015 LQ, F* 12.2; 2028 three spatial-hop groups per slot: 42 ancilla, 250 LQ, F* 10.4-13.0.
  (R1, e) first result 1.59 yr (1111 x 10a + 4444 x 5a at 30%); campaign 4-8 temperatures at 10%: 29.8-50.5 yr,
      6.0-10.1x the 5-yr horizon (at 30%: 3.3-5.6 yr). 2028: first result 36 s, campaign (L5 = 2, 3, 4) 13.2 min (L5 = 2 at its corrected wall, F* 8.9).
  The 50-point scan (n_coupling, 257 yr) and the per-point levers are records. PUBLISHED: 2028 LQ 240 -> 250;
  2033 LQ 1010 -> 1015, hard_ops 8.2e9 -> 1.6e10.

ROUND HISTORY (rulings applied, by round; records)
Round B rulings (H. Lamm, 2026-09-28) applied here:
  R5      Toffoli -> T is 7 T per Toffoli everywhere in the report, Ch. 9 included. No
          measurement-assisted 4-T Toffoli, no unitemized ancilla.
  VERTEX  each controlled 2O group action is priced at the sourced primitive U_x = 392 T
          (56 Toffoli at 7 T, 4 clean ancilla; arXiv:2312.10285 tab:tgatecost), ONE per
          vertex.  The unsourced '~190 T' is retired (groups.py keeps it for the record).
  R6      the 2028 adiabatic ramp is shortened to N_Trotter = 20.
  (R3: eps_l = 0.1 expected faults per shot.  R4: carry exact, round the printed result once.)
  R-HOP   (report-wide, H. Lamm 2026-09-29) the spatial hops on Z3 links are priced by the shared
          rule groups.rhop_link: C_W = 0 (Z3 diagonal), 8 Pauli strings per copy, one HWP group
          of k = Nc L5 = 12 per (link, string), rotations at the rule's 1.15 log2(1/eps) = 15.28 T
          (author's choice 2026-09-29, 'do 1': rhop_shared.md sec. 5 item 1, first option).
          The link-free fifth-direction hops and the electric terms keep the chapter's 30 T.
          2033: no change (C_W is ruling VERTEX already; an LCU query has no HWP term).
          [Price and 2033 pricing superseded by r17 below.]
  R-TOL   (report-wide, H. Lamm 2026-09-29) total synthesis error per shot eps_syn = 1e-2; each circuit
          sets eps_rot = sqrt(eps_syn / N_rot) (common.eps_rot_for), N_rot = its synthesized rotations in
          ONE shot.  The chapter keeps its formulas: the full RUS fit 1.15 log2(1/eps) + 9.2 for its own
          rotations, R-HOP's slope-only 1.15 log2(1/eps) for the hop rotations.  Replaces the fixed
          '~30 T at eps_gs ~ 1e-4' (t_per_rot, eps_rot are kept as retired records).
          2028 unchunked N_rot = 282 x 21 = 5922 -> eps 1.2995e-3, 20.226 / 11.026 T;
               chunked   N_rot = 304 x 21 = 6384 -> eps 1.2516e-3, 20.288 / 11.088 T.
          2033 N_rot = 25 x 390 x 691 = 6.73725e6 -> eps 3.853e-5, 26.063 T; c_rot 750 -> 651.6.
          [Counts superseded by r17 below; the rule itself stands.]

  r17     (H. Lamm 2026-10-01, rulings (e) and (f); apply_log/r17_ch09.md)
          2028: the Z3 spatial-hop rotations at the full fit 1.15 log2(1/eps) + 9.2 (R-HOP's slope-only
               price retired; Z3 has no Fermion_Primitives entry, so the R-HOP structure stays), and the
               missing on-site domain-wall term m + r(1+d) as one lattice-wide phasing group of k = 192.
          2033: the 81 link terms per BE query priced by groups.hop_link_cost('2O', Wilson, d=3): the
               authors' unpublished gate counts, color frame applied and undone (ESTIMATED), n_spin = 2
               (ESTIMATED), MBU ladders, no Trotter phasing; Higgs one U_x per site. Band top = the draft's W2.
               [Hop sharing superseded by r19 below.]

  r19     (H. Lamm 2026-10-01, rulings E21 (1)-(3); apply_log/r19_core.md, r19_ch09.md)
          (1) "make the switch": the color squish and parity flags are computed once per link and held through
              V(x), V(y), hop, V^dag(x), V^dag(y) (share="link", the new groups.hop_link_cost default), undo kept.
              The undo stays OUR estimate: the draft draws V_C^dag (PD/section_hamiltonian.tex:86-91) but its
              2 C^G counts V_C on the two sites only (PD/section_su2_diag.tex:238, section_resources.tex:15).
              The r17 per-application recompute (share="draft") is the conservative sensitivity record.
          (2) n_spin = 2 at Wilson d=3: author-confirmed.  (3) SU(3) color rotations = 2 N_angles: confirmed
              (no Ch. 9 number depends on it; the 2O hop is SU(2)).
          2028 (Z3): no change.  2033: per link 1182 -> 384 Toffoli (144 T, 704 rotations unchanged), so the
              rotation count and eps_rot do not move; c_be 2,557,184.8 -> 2,104,718.8, center 6.891e11 -> 5.672e11
              T/shot. Floor and W2 top do not move (band 3.383e9-1.2952e12 unchanged).

  r20     (H. Lamm 2026-10-01, rulings E23 (1)-(5); apply_log/r20_ch09.md)
          (1) "2O just is Clifford, that's a math fact": 2O is the lift to SU(2) of the 24-element octahedral rotation
              group, which is the single-qubit Clifford group modulo phases. Stated as a standard fact (no cite).
          (2) the band floor, investigated (DERIVED HERE): the old floor priced every link and Higgs term at zero T.
              Not justified: (a) the color index is two Jordan-Wigner modes, and U_g acts on them as the matchgate
              1 + U_g + 1, a two-qubit Clifford only for the 8 elements of Q8; (b) the link is a quantum register, so the
              frame is the multiplexer sum_g |g><g| (x) U_g, a controlled Clifford; (c) group multiplication on the 6-bit
              encoding is not affine, so U_x cannot be Clifford. What IS justified: 2O entries lie in Z[1/sqrt2, i], so no
              rotation synthesis. Explicit multiplexer (colour_multiplexer_2O; verified on all 48 elements in the tests):
              3 Toffolis + 18 T, 1 ancilla, 0 rotations. Floor link = 2 (W, W^dag) x n_spin x mux + the draft's spinor
              frame = 52 Toffolis + 104 T = 468 T; floor c_be 12,551.6 -> 61,043.6; band 3.383e9-1.2952e12 ->
              1.645e10-1.2952e12. The center stays the draft's frame (ruling e); promoting the multiplexer is Henry's call.
          (3) Gibbs sampler priced from first principles (ESTIMATED HERE): jumps x D_ov queries per jump x T per D_ov.
          (4) on-site term basis: kept the draft's diagonal-gamma^0 encoding (no number moves).
          (5) ancilla itemized, peak concurrency: 2028 37 (HWP k=32) + 2 (RUS) = 39 (was 40 asserted): 248 -> 247 LQ;
              2033 43 (was ~110 asserted): 1082 -> 1015 LQ (verifier: the parity qubit is held under share='link',
              so it is counted with p, q and the hop register, not as a transient; 42 -> 43).

  r21     (H. Lamm 2026-10-01: "1) promote 2) do TPQ 3) settle them"; apply_log/r21_ch09.md)
          (1) PROMOTE: the derived 2O color multiplexer (r20) is the 2033 headline (link_frame = "multiplexer"): 468 T
              per link per query, no rotations in the link. The draft's general diagonalizing frame is the stated
              alternative (link_frame = "draft"); the draft's W2 is a record. One headline, no band. The link-term
              workspace is now the multiplexer's: ancilla 43 -> 33, register 1015 -> 1005 LQ (the draft frame needs 43).
          (2) TPQ adopted as the thermal route; the Gibbs sampler is re-priced at the promoted c_be (each query one D_ov).
          (3) five items settled, each DERIVED HERE:
              (a) fifth-direction hops: open boundary, V (L5 - 1) = 24 (two groups of k = 72), not V L5 = 32;
              (b) N_Trotter (2033) from the second-order commutator bound (childs2021theory), Lambda from the term
                  norms: 83,195 steps at g^2 = 1 (trotter_bound_steps), not the stated ~30;
              (c) LCU terms per site (lcu_terms_per_site): 24 hop + 8 on-site + 40 Higgs = 72, 1944 at V = 27 (was the
                  absolute 1.7e3, 63 per site implied);
              (d) 2028 shot variance: sigma^2 = 10 x 27 / 24 = 11.25 for the 24 site-color copies taken as independent,
                  1125 shots (was 1e3 borrowed from V = 27);
              (e) 1+1D Z3 Wilson-kernel block encoding: 30 strings per site, 240 terms, c_be = 2.2e3 T; one sgn(H_W)
                  at kappa = 30 is 3.0e5 T, a ground-state preparation 3.1e6 T, 11x the domain-wall shot.

  r22     (H. Lamm 2026-10-01, rulings E26; apply_log/r22_core.md, r22_ch09.md)
          (i)   "for 3, do whatever is easiest and cheaper": the 2033 sign function is built in the single-particle
                architecture of arxiv_2607_28524 (Lamm, Roggero, Singh, Spagnoli), the only valid one. H = psi^dag h_ov psi
                with h_ov = gamma^0 + sgn(H_W[U]); sgn acts on the Q x Q kernel on a 9-qubit index register and the link
                registers, and two selected-Majorana unitaries lift it to the Fock register. Compiled here: one query is
                598 Toffolis + 24 T + 6 rotations (561 Toffolis read the addressed link; the 2O multiplexer on a color
                qubit is 10 T); one application is 692 queries + lift (862 Toffolis) + reflection + the Higgs vertex once
                (outside sgn) = 3.06e6 T. The r17-r21 construction (QSVT of sgn on a second-quantized block encoding)
                computes the sign of the many-body operator; it is kept as a record (62,896.2; 1.809e13).
                kappa = alpha / Delta (the paper's convention). Applications per Trotter step are derived: the
                Jacobi-Anger order at 2Q dt (block-encoding normalization 2Q = 432).
          (ii)  "i dont think we can do such a wildly off result": N_Trotter (2033) is a state-dependent estimate
                (Alves, Lamm, Liu, in preparation, FERMILAB-PUB-26-0397-T; weak-coupling evaluation derived here):
                5031 steps, range 1186-14,242; 13 applications per step; 2.008e11 T per shot. The worst-case commutator
                bound (83,195 steps, 1.276e12 T) is quoted as the worst case; trotter_rule_2033 = "state".
          (iii) "do 1" and E26 core ("do 2"): 2028 ancilla 31 (HWP k=32, k - w(k)) + 1 (RUS) = 32, register 240 LQ;
                single color copy 112; unchunked 190 + 1, 399 LQ.
          (iv)  2033 ancilla itemized for the single-particle architecture: 5 + 14 + 5 + 14 = 38, register 1010 LQ.
                The amplification reflection (arxiv_2407_17966) and the gauge and Higgs terms of a step are not priced
                (records).

Two chains, one per era.

2028 benchmark (1+1D Z3 domain-wall, L=8, L5=4; app06 2028 derivation and box)
    LQ  = 2 Dirac comp x Nc=3 x V=8 x L5=4 (=192) + D*V*ceil(log2 3) (=16) + 32 ancilla (31 HWP + 1 RUS, r22) = 240
    T   = (N_Trotter=20 + 1 probe step) x T/step, one step = 1128 Pauli rotations grouped by Hamming-weight
          phasing (HWP): k equal-angle rotations -> floor(log2 k)+1 synthesized R_Z + (k - HW(k)) Toffolis at 7 T.
              fifth-direction hops   2 groups of k=72   ->  14 synth, 140 Toffoli  (no link, 2 strings; r21: 24 hops)
              spatial hops          64 groups of k=12   -> 256 synth, 640 Toffoli  (8 links x 8 strings)
              electric terms         3 groups of k=8    ->  12 synth,  21 Toffoli
              on-site m + r(1+d)     1 group  of k=192  ->   8 synth, 190 Toffoli  (r17 f)
          R-TOL: 290 x 21 = 6090 rotations/shot -> eps 1.2814e-3, 20.249 T each (full fit, hops included)
          = 290 x 20.249 + 991 x 7 = 12,809.3 T/step, x 21 = 268,994.9 T/shot.
          Chunked to k <= 32 (31 ancilla; the k=192 group alone needs 190): fifth hops 2 x (2 x 32 + 8), on-site 6 x 32
          -> 336 synth, 985 Toffoli; 7056/shot -> eps 1.1905e-3, 20.371 T: 13,739.8 T/step, 288,535.5 T/shot.
          Records: with the pre-r21 32 fifth-direction hops 276,050.9-297,372.8; without the on-site term (2028-A, 32
          hops) 244,581.6-254,029.9; pre-r17 (slope-only hops) 195,122-204,571; before R-TOL 223,333-236,899.
          hard_ops = (unchunked, chunked); only the chunked end fits the 240-LQ register.
    shots = sigma^2 / eps^2 = 11.25 / 0.01 = 1125 (r21 d); wall = shots x (T x 1 us + 0.1 ms).

2033 target (3D 3^3 SU(2)->2O overlap + Higgs, TPQ thermal state; app06 2033 derivation and box)
    LQ  = 4x1x2x27 (=216) + 3x27x6 (=486) + 27x(4+6) (=270) + 38 ancilla (itemized; r22) = 1010
    T   = N_Dov x T_Dov, with
          d_sgn = ceil(kappa ln(1/delta_sgn)) = ceil(30 x 6.91) = 208   (r24; r17-r23 kappa 1e2: 691; kappa = alpha/Delta)
          [r24: the numbers below are the r22-r23 chain at kappa 1e2, t = 50a, kept as records; the r24 chain is in
          the r24 entry at the end of this docstring.]
          c_sp  = 598 Toffolis x 7 T + 24 T + 6 rotations: one query of the single-particle H_W block encoding
                  (flags 11, controlled shifts 12, read pass 78, link copy 405, M/M^dag ANDs 10, fix-up pass 78,
                  projector 4); R-TOL 4159 rotations per application x 65,643 = 273,009,237 -> eps 6.052e-6, 29.134 T:
                  4186 + 24 + 174.8 = 4384.8 (~4.4e3)
          T_Dov = 692 c_sp (3,034,285) + lift 6220.8 (862 Toffolis + 12 T + 6 rotations) + reflection 155.1
                  (18 Toffolis + 1 rotation) + Higgs 18,144 (27 U_x at 392 T + 1080 SELECT Toffolis) = 3,058,805.3
          N_Trotter = 5031: state-dependent second-order estimate, Lambda_sd = sigma/8 = 20.25 a^-3 from 156 transverse
                  oscillators at T a = 0.5, N = ceil(sqrt(Lambda t^3 / eps)), t = 50 a, eps = 0.1 (r22 ii)
          N_Dov = N_Trotter x 13 (Jacobi-Anger order at 2Q dt = 4.29, QSP error 1e-2 over the shot)
                  + d_beta(=30) x 8 amplification rounds = 65,403 + 240 = 65,643
          per shot = 65,643 x 3,058,805.3 = 2.008e11  (printed ~2.0e11); hard_ops = (that, that).
    Stated in the prose, not the box: the step-count range 1186-14,242 (1.131e11-3.495e11 T); the worst-case commutator
          bound, 83,195 steps at 5 applications per step (1.276e12); the fixed step of Ch. 5, 1000 steps at 35
          (1.077e11); the floor 2Q t = 21,600 applications (6.6e10).
    Records: the retired second-quantized construction at the worst-case step (r21 box: c_be 62,896.2, 4.346e7 per D_ov,
          1.809e13 per shot, 33 ancilla, 1005 LQ), and its earlier prints (T/shot at the pre-r21 390 x 691 queries and
          the 1.7e3-term LCU line): r20 center (draft frame) 5.672e11; r20 floor (multiplexer) 1.645e10; W2 top 1.2952e12;
          zero-T floor 3.383e9; pre-r17 one U_x per vertex 1.479e10; R-HOP two U_x per link 2.335e10; the draft's W1
          without the undo 4.342e11; r17 center (share="draft") 6.891e11; shared squish without the undo 3.935e11;
          n_spin = 1 3.120e11; 2n-3 ladders 5.981e11.
"""

from __future__ import annotations

import itertools
import math
from dataclasses import dataclass, fields, replace

from estimates.common import (Assumed, Cited, Stated, Uncited, Tagged, Primitive,
                              CircuitStatus, Result, Published, ERAS, T_PER_TOFFOLI,
                              SYNTHESIS_SRC, t_per_rotation, eps_per_rotation, toffoli_t,
                              EPS_SYN, EPS_SYN_SRC, eps_rot_for, RUS_ANCILLA, depth_exports,
                              REACTION_TIME_S, FACTORY_BASELINE)
from estimates.groups import (GROUPS, PRIMCOST, RHOP_SRC, RHOP_N_PAULI_DIAGONAL, rhop_link, hop_link_cost,
                              wilson_mass_rotations, FERMION_DRAFT, FP_KEY, FP_DIR)

TEX = "app06"
SEC_PER_YR = 3.156e7
_GROUPINGS = ("color", "color+s", "color+site")


# --------------------------------------------------------------------------- #
# small helpers (no synthesis formulas here; those live in common.py)
# --------------------------------------------------------------------------- #

def _popcount(k: int) -> int:
    return bin(int(k)).count("1")


def hwp_synth_rotations(k: int) -> int:
    """Hamming-weight phasing: k equal-angle rotations need floor(log2 k)+1 synthesized R_Z."""
    return math.floor(math.log2(k)) + 1


def hwp_toffolis(k: int) -> int:
    """Adder-tree Toffolis to compute the Hamming weight of k bits: k - HW(k) (arXiv:1902.10673),
    the count the chapter quotes as 'k - w(k)' (app06:108); at most k."""
    return int(k) - _popcount(k)


def hwp_ancilla(k: int) -> int:
    """Ancilla of one Hamming-weight-phasing group: k - HW(k), one carry ancilla per adder Toffoli
    (arXiv:1902.10673 App. A; arXiv:1709.06648). The weight bits are those carries and the lowest bit sits on a
    data qubit, so no weight register is added: 31 at k=32, 190 at k=192, 7 at k=9. Author ruling E26 (r22,
    2026-10-01); until then this returned Toffolis + floor(log2 k)+1 (37 at k=32), counting the weight register
    twice. Same value as common.hwp_ancilla (tests/test_rhop.py)."""
    return hwp_toffolis(k)


def _chunk(k: int, n_groups: int, chunk_k: int | None) -> list[tuple[int, int]]:
    """[(k, n_groups), ...] after splitting groups larger than chunk_k into groups of chunk_k
    (plus one remainder group when chunk_k does not divide k)."""
    if chunk_k is None or k <= chunk_k:
        return [(k, n_groups)]
    q, r = divmod(k, chunk_k)
    out = [(chunk_k, n_groups * q)]
    if r:
        out.append((r, n_groups))
    return out


# --------------------------------------------------------------------------- #
# E23 (2), r20: the 2O color multiplexer (DERIVED HERE)
# --------------------------------------------------------------------------- #

def colour_multiplexer_2O() -> dict:
    """W = sum_g |g><g| (x) M(U_g) on one site's two color modes (one spinor component), DERIVED HERE (r20).

    M(U) is the Jordan-Wigner image of the color rotation psi -> U psi: identity on |00>, |11> and U on the singly
    occupied pair (|01>, |10>) (FermionPrimitives_unpub section_su2_diag.tex:27-60). The 2O register uses the ordered
    product g = (-1)^x1 j^x2 k^x3 u^(2x4+x5) t^x6 of arXiv:2312.10285, so M(U_g) is a product of six singly controlled
    generator matchgates. A CNOT (q2 -> q1) turns q1 into the odd-occupancy flag and q2 into the color qubit of the
    singly occupied sector, so each factor is the generator on q2 doubly controlled by (x_i, flag):
        x1  -1                 CZ(x1, flag)                                   0 Toffoli, 0 T
        x2  j = iY             exp(i pi/2 x2 flag Y) = 4 pi/8 Pauli rotations 0 Toffoli, 4 T
        x3  k = iZ             exp(i pi/2 x3 flag Z) = 4 pi/8 Pauli rotations 0 Toffoli, 4 T
        x4  u^2 = iZ R(X) R(-Z)   f = AND(x4, flag); two controlled R(P) = exp(-i pi/4 f P), 2 T each; C_f-(iZ) Clifford
        x5  u = -R(X) R(-Y)    same: 1 Toffoli + 4 T
        x6  t = R(X)           1 Toffoli + 2 T
    R(P) = exp(-i pi/4 P); exp(-i pi/4 f P) = exp(-i pi/8 P) exp(i pi/8 Z_f P), two T-equivalent rotations. Every f is
    uncomputed by measurement (MBU, 0 T). Verified against all 48 elements in tests/test_ch09_chiral_gauge.py.
    Total 3 Toffolis + 18 T, one reusable ancilla, no synthesized rotations. Not proven optimal: the floor is the
    cheapest known circuit. It is not zero: a controlled non-Pauli Clifford is not Clifford.
    """
    per_bit = {"x1": (0, 0), "x2": (0, 4), "x3": (0, 4), "x4": (1, 4), "x5": (1, 4), "x6": (1, 2)}
    return {"per_bit": per_bit, "toffoli": sum(v[0] for v in per_bit.values()),
            "t_direct": sum(v[1] for v in per_bit.values()), "n_rot": 0, "ancilla": 1}


def multiplexer_link_2O(n_spin: int, spinor: dict) -> dict:
    """The link term per BE query with the derived multiplexer (r20 floor, r21 headline): W and W^dag on the far site
    for each of the n_spin hopping spinor components, plus the draft's spinor frame (groups.spinor_counts). The hop
    between the near site and the rotated far site is color-diagonal and g-independent, so it enters SELECT as fixed
    Pauli strings (counted in lcu_terms_per_site): no hop squish, no diagonalizers, no color squish."""
    mx = colour_multiplexer_2O()
    n_mux = 2 * n_spin
    return {"n_mux": n_mux, "toffoli": n_mux * mx["toffoli"] + spinor["toffoli"],
            "t_direct": n_mux * mx["t_direct"] + spinor["t_direct"], "n_rot": 0,
            "mux_toffoli": n_mux * mx["toffoli"], "mux_t_direct": n_mux * mx["t_direct"],
            "spinor_toffoli": spinor["toffoli"], "spinor_t_direct": spinor["t_direct"]}


# --------------------------------------------------------------------------- #
# r21 (H. Lamm 2026-10-01, "settle them"): LCU term count and Trotter step count (DERIVED HERE)
# --------------------------------------------------------------------------- #

def lcu_terms_per_site(d: int, n_spin: int, n_c: int, n_f: int, strings_per_bilinear: int,
                       radial_terms: int = 0) -> dict:
    """Pauli strings per site of the block-encoded Wilson kernel in the Jordan-Wigner encoding (DERIVED HERE, r21 c, e).

    hop     d links per site x n_spin hopping spinor components (r = 1: the projector 1 - gamma_k keeps half of them, and
            the spinor frame makes the hop spin-diagonal) x n_c n_f color-flavor copies x strings per bilinear.
            A color-diagonal bilinear a^dag b + h.c. is XX + YY with a Z string: 2 strings (2033, after the color
            multiplexer). On a Z3 link the phase omega^g multiplies the bilinear and enters as Z strings on the link
            register: 8 strings per copy (groups.RHOP_N_PAULI_DIAGONAL, derived in Ch. 6 round D).
    onsite  the Wilson on-site term psi^dag gamma^0 (m + r d) psi: 2 N (d <= 2) or 4 N (d = 3) number operators, one
            Z string each (groups.wilson_mass_rotations; FermionPrimitives_unpub section_mass.tex:12).
    higgs   the Higgs vertex: the same on-site bilinear, dressed by the Higgs group element (the frame the U_x line
            prices) and weighted by the radial mode rho = c_0 - sum_j c_j Z_j on its log2 N_rho qubits:
            onsite x radial_terms strings (radial_terms = log2 N_rho + 1; 0 where there is no Higgs).
    """
    hop = d * n_spin * n_c * n_f * strings_per_bilinear
    onsite = wilson_mass_rotations(n_c, n_f, d)
    higgs = onsite * radial_terms
    return {"hop": hop, "hop_per_link": hop // d, "onsite": onsite, "higgs": higgs, "total": hop + onsite + higgs}


def trotter_bound_steps(n_links: int, t_lat: float, eps: float, g2: float, casimir_max: float, hop_norm: float,
                        n_plaq_per_link: int) -> dict:
    """Second-order Trotter step count from the commutator bound (DERIVED HERE, r21 b).

    One step is S(dt) = e^{-i D dt/2} e^{-i E dt} e^{-i D dt/2}: E the electric terms, D everything diagonal in the
    group basis (the plaquettes and the fermion kernel, which commute with each other; adjacent D half-steps merge, so
    the kernel is applied once per step). childs2021theory, second-order formula for two terms:
        ||S(dt) - e^{-i H dt}|| <= dt^3 ( ||[E,[E,D]]|| / 12 + ||[D,[D,E]]|| / 24 ) = dt^3 Lambda,
    so N steps over a time t give at most (t^3 / N^2) Lambda, and N = sqrt(Lambda t^3 / eps).

    Lambda from the term norms in lattice units, same-link terms only, ||[X,[Y,Z]]|| <= 4 ||X|| ||Y|| ||Z|| with
    half-ranges:
        e = (g^2 / 2) C_max / 2      one link's electric term (g^2/2) E^2, spectrum [0, (g^2/2) C_max]
        b = hop_norm + n_plaq 4/g^2  the terms that depend on that link: the Wilson hop (one unit-norm bilinear per
                                     hopping (spinor, color) mode at r = 1) and the plaquettes containing it
                                     (Kogut-Susskind -(2/g^2) Re Tr U_p, |Re Tr U_p| <= 2 for SU(2))
        Lambda = n_links (e^2 b / 3 + e b^2 / 6).
    The neglected cross-link terms (hops sharing a site, plaquettes sharing two links) and the Higgs terms only add.
    The norm inequality is NOT tight (r21 verifier): on one link, in the class-function sector j <= 3/2 with the four
    plaquettes aligned (D = (8/g^2) chi_{1/2}, E = (g^2/2) C), the exact double commutators give
    ||[E,[E,D]]||/12 + ||[D,[D,E]]||/24 = 0.57 + 4.74 = 5.31 at g^2 = 1, against 4.69 + 40 = 44.69 from the inequality
    (8.4x in Lambda, 2.9x in N for that part). So the step count is good to a factor of a few, not to 1.5x; the
    chapter says so (tests: test_r21_norm_inequality_is_loose_on_one_link).
    The order with D outside is the better of the two (E outside gives e^2 b / 6 + e b^2 / 3).
    """
    e = 0.5 * (g2 / 2.0) * casimir_max
    b = hop_norm + n_plaq_per_link * 4.0 / g2
    lam_link = e * e * b / 3.0 + e * b * b / 6.0
    lam = n_links * lam_link
    n_exact = math.sqrt(lam * t_lat ** 3 / eps)
    return {"e_half_range": e, "b_norm": b, "lambda_link": lam_link, "lambda": lam, "n_exact": n_exact,
            "n_steps": math.ceil(n_exact), "dt_lat": t_lat / math.ceil(n_exact),
            "lambda_link_other_order": e * e * b / 6.0 + e * b * b / 3.0}


# --------------------------------------------------------------------------- #
# r22 (E26 (i)): the overlap kernel in the single-particle architecture of arxiv_2607_28524, compiled here
# --------------------------------------------------------------------------- #

def colour_multiplexer_index_2O() -> dict:
    """M = sum_g |g><g| (x) U_g on the ONE color qubit of the index register (DERIVED HERE, r22).

    On the index register color is a qubit, not two Jordan-Wigner modes, so no occupancy flag is needed and the
    six bits of g = (-1)^x1 j^x2 k^x3 u^(2x4+x5) t^x6 (arXiv:2312.10285) each control one fixed generator:
        x1  -1                  a Z on the bit (as a sign of the read link it is a CZ with the address leaf)   0 T
        x2  j = iY, x3  k = iZ  controlled Pauli times S on the control: Clifford                             0 T
        x4  u^2, x5  u          phase x Pauli x R(A) R(B), R(P) = exp(-i pi/4 P); each controlled R is two
                                pi/8 Pauli rotations: 4 T per bit
        x6  t = R(X)            one controlled R: 2 T
    10 T, no Toffoli, no ancilla, no synthesized rotation. Checked on all 48 elements, with the inverse, in
    tests/test_ch09_chiral_gauge.py. The Fock-space version (colour_multiplexer_2O, r20) is 3 Toffolis + 18 T."""
    per_bit = {"x1": (0, 0), "x2": (0, 0), "x3": (0, 0), "x4": (0, 4), "x5": (0, 4), "x6": (0, 2)}
    return {"per_bit": per_bit, "toffoli": sum(v[0] for v in per_bit.values()),
            "t_direct": sum(v[1] for v in per_bit.values()), "n_rot": 0, "ancilla": 0}


def _ctrl_shift_toffolis(L: int) -> int:
    """Toffolis of one controlled shift x -> x + 1 mod L on a ceil(log2 L)-bit coordinate (ANDs uncomputed by
    measurement). L = 3 on two bits, valid subspace {0, 1, 2}: a controlled SWAP and a doubly controlled X, 2.
    L = 2^n: bit j flips under the control and the j lower bits, a C^(j+1)X at j Toffolis: n (n - 1) / 2."""
    if L == 3:
        return 2
    n = L.bit_length() - 1
    if L >= 2 and 1 << n == L:
        return n * (n - 1) // 2
    raise ValueError(f"controlled shift mod {L} is not derived here (L = 3 or a power of two)")


def single_particle_query_2O(D: int, L: int, link_qubits: int) -> dict:
    """One block-encoding query U_h = P^dag S P of the single-particle Wilson kernel H_W[U] on the index register
    (color, spin, D coordinates) and the link registers (DERIVED HERE, r22; construction of arxiv_2607_28524).

    LCU of 4D + 1 unitaries: the on-site term and, for each direction and orientation, a Dirac and a Wilson hop.
    Term register: mass/hop, direction, orientation, Dirac/Wilson (5 qubits at D = 3), so the +-i and -1
    coefficients are Clifford phases on flags. SELECT, in time order:
        direction flags     hop AND direction = i; the last one by CNOT                              D - 1
        backward flags      flag AND orientation, for the shift A_i^dag before the read              D
        controlled shifts   A_i^dag (backward, before the read) and A_i (forward, after it)          2 D x shift(L)
        read pass           unary iteration over the L^D sites below each direction flag             D (L^D - 1)
                            (Babbush_PRX_2018)
        link copy           bits x2..x6 of the addressed link into a scratch register: the link is   (b - 1) D L^D
                            quantum data, one Toffoli per bit; x1 is a sign and enters as CZ(leaf, x1)
        M or M^dag          the color multiplexer or its inverse, chosen by the orientation qubit   2 (b - 1)
        fix-up pass         the copy is measured out in X; its phase fix-ups walk the addresses again  D (L^D - 1)
        forward flags       flag AND not-orientation, recomputed (holding the backward flags through  D
                            the read would save these and cost D more qubits)
        Dirac flags         flag AND Dirac, for the controlled gamma matrices (Clifford)              D
        projector           all-zero flag on the term register for the QSP reflection               n_term - 1
    Direct T: M and M^dag at 10 T each, and a controlled-H in the uniform-over-D PREPARE and in its inverse (2 T each).
    Rotations: PREPARE and its inverse (the mass/hop split and the uniform-over-D state: 2 each) and the QSP phase,
    controlled on the LCU qubit of h_ov (2). Gate-level simulation on the full 3^3 lattice with a random 2O
    configuration: scratchpad ch9ov/recheck_select.py; ported to the tests (sampled columns)."""
    b = int(link_qubits)
    n_sites, n_links = L ** D, D * L ** D
    n_dir = math.ceil(math.log2(D)) if D > 1 else 0
    n_term = 1 + n_dir + 1 + 1
    mx = colour_multiplexer_index_2O()
    tof = {"direction_flags": D - 1, "backward_flags": D, "forward_flags": D,
           "controlled_shifts": 2 * D * _ctrl_shift_toffolis(L), "dirac_flags": D,
           "read_iteration": D * (n_sites - 1), "link_copy": (b - 1) * n_links, "mux_orientation": 2 * (b - 1),
           "fixup_iteration": D * (n_sites - 1), "projector": n_term - 1}
    t_direct = {"mux_and_inverse": 2 * mx["t_direct"], "prepare_controlled_h": 4}
    rot = {"prepare_and_inverse": 4, "qsp_phase_controlled": 2}
    link_read = tof["read_iteration"] + tof["link_copy"] + tof["fixup_iteration"]
    radix_lines = 2 if L == 3 else 1                  # AND lines per level of the iteration tree (ternary / binary)
    return {"toffoli_items": tof, "toffoli": sum(tof.values()), "t_direct_items": t_direct,
            "t_direct": sum(t_direct.values()), "rot_items": rot, "n_rot": sum(rot.values()),
            "link_read_toffoli": link_read, "n_lcu_unitaries": 4 * D + 1, "term_register": n_term,
            "anc_direction_flags": D, "anc_iteration": radix_lines * D * (1 if L == 3 else math.ceil(math.log2(L))),
            "anc_link_copy": b - 1}


def single_particle_query_Z3_1p1d(L: int, link_qubits: int) -> dict:
    """The same query for the 1+1D Z3 kernel of the 2028 comparison (DERIVED HERE, r22): index register of color
    (2 qubits), spin (1) and a log2 L-bit coordinate; 4D + 1 = 5 LCU terms on a 3-qubit term register (hop,
    orientation, Dirac/Wilson). Z3 is diagonal on the color index, so the link enters as the phase omega^(+-g):
    one synthesized rotation on each copied bit."""
    n = L.bit_length() - 1
    tof = {"backward_flag": 1, "controlled_shifts": 2 * _ctrl_shift_toffolis(L), "dirac_flag": 1,
           "read_iteration": L - 1, "link_copy": link_qubits * L, "fixup_iteration": L - 1, "projector": 3 - 1}
    rot = {"prepare_and_inverse": 2, "link_phase": link_qubits, "qsp_phase_controlled": 2}
    return {"toffoli_items": tof, "toffoli": sum(tof.values()), "t_direct": 0, "rot_items": rot,
            "n_rot": sum(rot.values()), "n_lcu_unitaries": 5, "term_register": 3, "index_register": 2 + 1 + n}


def lift_cost(q_modes: int, n_uniform3: int) -> dict:
    """The lift of the kernel to the Fock register (arxiv_2607_28524, Supplement theorem):
    <0| U_L h U_R |0> = psi^dag h psi / Q, with U_L, U_R each one of 2Q Majorana strings (type X or Y, mode a)
    selected by unary iteration with an accumulator for the Jordan-Wigner string: 2Q - 1 Toffolis each. Q, not
    2^q, because the index register is prepared uniform over the valid modes: each of its n_uniform3 three-valued
    sub-registers (the D coordinates at L = 3; the color index of the 1+1D Z3 kernel) is prepared uniform over 3
    values, in and out (one rotation and one controlled-H, 2 T, each way); two-valued and 2^n-valued ones by Hadamards."""
    return {"toffoli": 2 * (2 * q_modes - 1), "t_direct": 2 * n_uniform3 * 2, "n_rot": 2 * n_uniform3, "leaves": 2 * q_modes,
            "anc_iteration": math.ceil(math.log2(2 * q_modes)) + 1}


def bessel_j_list(x: float, kmax: int) -> list:
    """J_0..J_kmax(x) by Miller's backward recurrence, normalized with J_0 + 2 sum_k J_2k = 1 (stdlib only)."""
    if x == 0:
        return [1.0] + [0.0] * kmax
    start = int(max(kmax, x) + 40 + 2 * math.sqrt(max(kmax, x)))
    start += start % 2
    jp, j = 0.0, 1e-300
    vals = [0.0] * (start + 1)
    vals[start] = j
    for k in range(start, 0, -1):
        jp, j = j, 2.0 * k / x * j - jp
        vals[k - 1] = j
        if abs(j) > 1e250:
            vals = [v * 1e-250 for v in vals]
            j *= 1e-250
            jp *= 1e-250
    norm = vals[0] + 2.0 * sum(vals[2::2])
    return [v / norm for v in vals[:kmax + 1]]


def jacobi_anger_order(x: float, eps: float) -> int:
    """Smallest R with 2 sum_{k>R} |J_k(x)| <= eps: the degree of the Jacobi-Anger polynomial that evolves a
    block-encoded H / alpha for a time with alpha dt = x to error eps (low2019hamiltonian). It equals the number of
    applications of the block encoding when the block encoding is self-inverse, which the lifted kernel can be made
    at O(1) cost; otherwise the count is twice the order."""
    kmax = int(x + 60 + 12 * x ** (1.0 / 3.0))
    J = bessel_j_list(x, kmax)
    tail, R = 0.0, kmax
    while R > 0 and 2.0 * (tail + abs(J[R])) <= eps:
        tail += abs(J[R])
        R -= 1
    return R


def dov_per_step(n_steps: int, t_lat: float, q_modes: int, eps_qsp: float) -> int:
    """Applications of the block-encoded overlap Hamiltonian per Trotter step: normalization 2Q in lattice units,
    QSP error eps_qsp / N per step (eps_qsp over the shot)."""
    return jacobi_anger_order(2.0 * q_modes * t_lat / n_steps, eps_qsp / n_steps)


def trotter_state_dependent_steps(L: int, n_adj: int, t_lat: float, eps: float, temp_lat: float, d: int = 3) -> dict:
    """Second-order Trotter step count from the state-dependent estimate (r22, E26 (ii)). The replacement of operator
    norms by expectation values in the simulated state is that of Alves, Lamm and Liu (preprint,
    FERMILAB-PUB-26-0397-T, App. 'Second-order Trotter approximation'); the evaluation below is DERIVED HERE.

    Step S(dt) = e^{-iD dt/2} e^{-iE dt} e^{-iD dt/2}. For a state rho that commutes with H, to leading order in dt,
        || (U - e^{i theta} S^N) rho^{1/2} ||  <=  (t^3 / N^2) ( sigma(o1)/12 + sigma(o2)/24 ),
    o1 = [E,[E,D]], o2 = [D,[D,E]], sigma(O)^2 = <O^dag O> - |<O>|^2 in rho. <o1> = <o2> in any such state, and the
    common mean shifts only the global phase theta, which cancels in the measured correlator. Keeping the mean (the
    second moment, the draft's formula as written) gives 'n_high'.
    Weak-coupling evaluation: each transverse gauge mode is an oscillator with w_k^2 = 4 sum_i sin^2(k_i/2);
    o1 = -w^3 P^2, o2 = -w^3 X^2; thermal mean -(w^3/2) coth(w/2T), variance 2 [(w^3/2) coth(w/2T)]^2 for each.
    n_adj colors x (d - 1) polarizations per non-zero momentum; the variances add over modes, so sigma grows as the
    square root of the volume. At this order there is no dependence on g or on the electric cutoff.
    'n_low': the phase of the fastest mode right to eps over the whole evolution, w^3 dt^2 t / 24 = eps (the
    frequency shift of the second-order step)."""
    ws = []
    for n in itertools.product(range(L), repeat=d):
        w2 = 4.0 * sum(math.sin(math.pi * ni / L) ** 2 for ni in n)
        if w2 > 1e-12:
            ws.append(math.sqrt(w2))
    mult = n_adj * (d - 1)

    def coth(w):
        return 1.0 / math.tanh(w / (2.0 * temp_lat))
    mean = mult * sum(0.5 * w ** 3 * coth(w) for w in ws)
    var = mult * sum(2.0 * (0.5 * w ** 3 * coth(w)) ** 2 for w in ws)
    lam_var = math.sqrt(var) / 8.0                      # 1/12 + 1/24
    lam_full = math.sqrt(mean ** 2 + var) / 8.0
    w_max = max(ws)

    def n_of(lam):
        return math.ceil(math.sqrt(lam * t_lat ** 3 / eps))
    return {"n_modes": mult * len(ws), "n_momenta": len(ws), "w_max": w_max, "mean": mean, "std": math.sqrt(var),
            "lambda_var": lam_var, "lambda_full": lam_full, "n_central": n_of(lam_var), "n_high": n_of(lam_full),
            "n_low": math.ceil(t_lat * math.sqrt(w_max ** 3 * t_lat / (24.0 * eps)))}


# --------------------------------------------------------------------------- #
# Referee C1 / C3 follow-up (chapter-open pass, 2026-10-05): thermal reference state and Delta N_CS readout.
# Scratch: chopen/app06_chiral_gauge/{free_overlap_tpq,proxy_tpq,ncs_freefield}.py (+ .out).
# --------------------------------------------------------------------------- #

def free_overlap_modes(L: int, d: int, m: float, r: float) -> list:
    """Single-particle spectrum of h_ov = Gamma^0 + sgn(h_W) at U = 1 on L^d (periodic), with the Wilson kernel of
    arXiv:2607.28524 (h_W(p) = sum_i Gamma^i sin p_i + Gamma^0 c(p), c = m + r sum_i (1 - cos p_i); m = -1.5 gives the
    chapter's alpha = 7.5). DERIVED HERE: h_ov(p) = Gamma^0 (1 + c/N) + (b.Gamma)/N, N = sqrt(b^2 + c^2), so its
    eigenvalues are +-sqrt(2 + 2c/N), each twice (4 spinor components per momentum); zero at p = 0 (c < 0, b = 0).
    Checked against a dense diagonalization on 3^3 (scratch free_overlap_tpq.py). Returns one color's 4 L^d modes."""
    out = []
    for n in itertools.product(range(L), repeat=d):
        p = [2.0 * math.pi * ni / L for ni in n]
        b2 = sum(math.sin(pi) ** 2 for pi in p)
        c = m + r * sum(1.0 - math.cos(pi) for pi in p)
        lam = math.sqrt(max(0.0, 2.0 + 2.0 * c / math.sqrt(b2 + c * c)))
        out += [lam, lam, -lam, -lam]
    return out


def typical_reference_tpq(modes: list, beta: float, two_q: float) -> dict:
    """Referee C1 (DERIVED HERE). For a reference ensemble whose average is I/D (Haar, random phase, or a uniformly
    random basis state), post-selecting on the filter e^{-beta(H - E_lb)/2} gives exact thermal averages, and the mean
    success probability is p = Tr e^{-beta(H - E_lb)} / D. For free fermions it factorizes over modes: at E_lb = E0,
    p = prod (1 + e^{-beta|eps|}) / 2. Amplitude amplification needs (pi/4) p^{-1/2} rounds. Normalized at -2Q instead
    of E0 (the filter as priced), p loses a further e^{-beta (E0 + 2Q)}."""
    e0 = sum(e for e in modes if e < 0)
    ln_p = sum(math.log((1.0 + math.exp(-beta * abs(e))) / 2.0) for e in modes)
    lnz_shift = sum(math.log1p(math.exp(-beta * abs(e))) for e in modes)      # ln Z + beta E0
    e_th = e0 + sum(abs(e) / (math.exp(beta * abs(e)) + 1.0) for e in modes)   # thermal energy
    s_th = beta * (e_th - e0) + lnz_shift
    return {"n_modes": len(modes), "E0": e0, "E_minus_E0": e_th - e0, "S": s_th, "lnD": len(modes) * math.log(2),
            "ln_p": ln_p, "log10_p": ln_p / math.log(10), "aa_rounds": math.pi / 4.0 * math.exp(-ln_p / 2.0),
            "n_frozen": sum(1 for e in modes if beta * abs(e) > 2.0),
            "ln_norm_offset": -beta * (e0 + two_q)}


def ncs_free_field_background(L: int, n_adj: int, temp_lat: float, g2: float, t_lat: float, d: int = 3) -> float:
    """Referee C3 (DERIVED HERE). Non-topological <(Delta N_CS)^2>(t) of the free (abelianized) gauge field on L^d in
    the quantum thermal state. With canonical fields int E.B = -d/dt (1/2) int A.B, so N = (kappa/2) sum w (|A_+|^2 -
    |A_-|^2), kappa = g^2 / 8 pi^2, each helicity a real scalar field; Wick with the Wightman function gives
    <(Delta N)^2>(t) = (kappa^2/4) sum_{k, hel, color} ((n+1)^2 + n^2)(1 - cos 2wt), n = 1/(e^{w/T} - 1).
    Bounded in t (no diffusion in the free theory) and dominated by zero-point fluctuations at T a = 0.5.
    VERIFIER 2026-10-05: this is the IDEAL helicity-A.B operator (curl eigenvalues +-w), not the chapter's lattice
    estimator; kept as a record (8.93e-3 at 10a). The chapter's estimator is ncs_lattice_background."""
    kappa = g2 / (8.0 * math.pi ** 2)
    tot = 0.0
    for n in itertools.product(range(L), repeat=d):
        w2 = 4.0 * sum(math.sin(math.pi * ni / L) ** 2 for ni in n)
        if w2 < 1e-12:
            continue
        w = math.sqrt(w2)
        nb = 1.0 / math.expm1(w / temp_lat)
        tot += n_adj * (d - 1) * ((nb + 1.0) ** 2 + nb ** 2) * (1.0 - math.cos(2.0 * w * t_lat))
    return kappa ** 2 / 4.0 * tot


def ncs_lattice_background(L: int, n_adj: int, temp_lat: float, g2: float, t_lat: float, pairing: str = "avg",
                           n_quad: int = 0, statistics: str = "quantum") -> float:
    """Referee C3, verifier fix (DERIVED HERE, 2026-10-05). <(Delta N_CS)^2>(t) of the chapter's lattice estimator in
    the free abelian theory on L^3 (quantum thermal state, n_adj colors): Qdot = kappa sum_{x,k} E_k(x) B_k(x), B_k(x)
    the clover (average of the four k-plaquettes at x), E_k(x) the field on link (x,k) ("fwd") or the average over links
    (x,k) and (x-k,k) ("avg"). With z = (q, p) in the transverse normal modes, Delta N = z0^T A z0, A = kappa int_0^t
    S(tau)^T W S(tau) dtau (Gauss-Legendre), and <(Delta N)^2> = 2 Re Tr(A G0 A G0^T), G0 = <z z^T>. "avg" gives a
    symmetric mode matrix (a total derivative, bounded, 9.58e-4 at 10a, 95% zero point); "fwd" has an antisymmetric
    part (norm 0.43 vs 1.14) whose degenerate-mode piece is conserved, so it grows as t^2 (2.56e-3 at 10a, 2.4e-2 at
    40a). Agrees with the verifier's independent Kubo double integral (scratch ch09_check/ncs_check.py) to 2e-4.
    Scratch: chopen/app06_chiral_gauge/ncs_lattice.py.
    statistics = "classical" (author ruling 2026-10-05, method test): the classical-statistical (Rayleigh-Jeans) value,
    <q^2> = T/w^2, <p^2> = T, no commutator, i.e. what a classical real-time simulation of the free lattice gives:
    1.34e-4 at 10a, 14% of the quantum 9.58e-4 (no zero-point part). Scratch rul2/ch09/classical_bg.py."""
    import numpy as np
    sites = list(itertools.product(range(L), repeat=3))
    idx = {s: k for k, s in enumerate(sites)}
    nl = 3 * len(sites)

    def sh(x, i, dd=1):
        return tuple((x[k] + (dd if k == i else 0)) % L for k in range(3))

    def plaq(x, i, j):
        r = np.zeros(nl)
        r[3 * idx[x] + i] += 1
        r[3 * idx[sh(x, i)] + j] += 1
        r[3 * idx[sh(x, j)] + i] -= 1
        r[3 * idx[x] + j] -= 1
        return r

    K = np.zeros((nl, nl))
    for x in sites:
        for i, j in ((0, 1), (1, 2), (2, 0)):
            r = plaq(x, i, j)
            K += np.outer(r, r)
    M = np.zeros((nl, nl))
    for x in sites:
        for k in range(3):
            i, j = (k + 1) % 3, (k + 2) % 3
            B = (plaq(x, i, j) + plaq(sh(x, i, -1), i, j) + plaq(sh(x, j, -1), i, j)
                 + plaq(sh(sh(x, i, -1), j, -1), i, j)) / 4.0
            if pairing == "fwd":
                M[3 * idx[x] + k] += B
            else:
                M[3 * idx[x] + k] += B / 2.0
                M[3 * idx[sh(x, k, -1)] + k] += B / 2.0
    w2, Vm = np.linalg.eigh(K)
    ph = w2 > 1e-9
    Vp, w = Vm[:, ph], np.sqrt(w2[ph])
    N = len(w)
    nb = 1.0 / np.expm1(w / temp_lat)
    W = np.zeros((2 * N, 2 * N))
    W[N:, :N] = Vp.T @ M @ Vp
    W = (W + W.T) / 2.0
    xg, wg = np.polynomial.legendre.leggauss(n_quad or max(400, int(40 * t_lat)))
    A = np.zeros((2 * N, 2 * N))
    for tt, ww in zip((xg + 1.0) * t_lat / 2.0, wg * t_lat / 2.0):
        c, s = np.cos(w * tt), np.sin(w * tt)
        S = np.block([[np.diag(c), np.diag(s / w)], [np.diag(-w * s), np.diag(c)]])
        A += ww * (S.T @ W @ S)
    A *= g2 / (8.0 * math.pi ** 2)
    G = np.zeros((2 * N, 2 * N), complex)
    if statistics == "classical":
        G[:N, :N] = np.diag(temp_lat / w ** 2)
        G[N:, N:] = np.diag(temp_lat * np.ones(N))
    else:
        G[:N, :N] = np.diag((2 * nb + 1) / (2 * w))
        G[N:, N:] = np.diag(w * (2 * nb + 1) / 2)
        G[:N, N:] = 0.5j * np.eye(N)
        G[N:, :N] = -0.5j * np.eye(N)
    return float(n_adj * 2.0 * np.real(np.trace(A @ G @ A @ G.T)))


def ncs_hadamard_lambda(V: int, d: int, n_adj: int, g2: float, e_max: float, b_max: float) -> float:
    """Referee C3 (DERIVED HERE): sum of ||E^a|| ||B^a|| over the terms of Qdot = (g^2/8pi^2) sum_{x,i,a} E^a_i B^a_i,
    ||E^a|| = j_max = 3/2 on the 2O irreps restricting SU(2) j <= 3/2, ||B^a|| = max |Tr(sigma^a U)| = 2. A LOWER bound
    on the block-encoding normalization (verifier 2026-10-05): a Pauli/LCU decomposition of the spin-3/2 generator is
    larger. Chapter: 'lambda >= 9.2/a'."""
    return g2 / (8.0 * math.pi ** 2) * V * d * n_adj * e_max * b_max


def ncs_hadamard_rel_var(lam: float, t_lat: float, signal: float) -> float:
    """Referee C3 (DERIVED HERE). <(Delta N)^2>(t) = 2 int_0^t (t - tau) Re C(tau) dtau from C(tau) = <Qdot(tau)Qdot(0)>
    at N_tau points; each +-1 Hadamard-test outcome estimates C/lambda^2, so with N shots spread evenly over the points
    Var = (4/3) lambda^4 t^4 / N, and the relative variance per shot is (4/3) (lambda t)^4 / <(Delta N)^2>^2."""
    return 4.0 / 3.0 * (lam * t_lat) ** 4 / signal ** 2


# --------------------------------------------------------------------------- #
# r25 (ruling R9/R10, H. Lamm 2026-10-02): T-depth of one shot (factory analysis, apply_log/r25_ch09.md)
# --------------------------------------------------------------------------- #
# Reaction-limited depth: a Toffoli or AND is one reaction layer; a repeat-until-success rotation is a sequential
# Clifford+T chain, so its T-depth equals its T-count. Shots are serial on one machine.

# One query of the single-particle H_W (single_particle_query_2O), reaction layers per item. 'register': on the
# 1010-LQ register as itemized; all 38 ancilla are live during a query, so Toffolis that share a control (the 5 copy
# Toffolis of a leaf, the flags that share the orientation or Dirac qubit, the mux ANDs) are serial. The controlled
# shifts act on disjoint (flag, coordinate) pairs and run in parallel across directions (12 -> 2 + 2). 'fanout': with
# a CNOT fan-out of the shared control onto FANOUT_ANCILLA_2033 extra qubits, those groups run in parallel. The read
# and fix-up walks (78 + 78) are serial in every variant (unary iteration is a sequential walk and the directions share
# the iteration temporaries); M and M^dag are 20 serial pi/8 rotations on one color qubit.
SP_QUERY_DEPTH_2O = {
    "register": {"direction_flags": 2, "backward_flags": 3, "shift_backward": 2, "read_iteration": 78,
                 "link_copy": 405, "mux_orientation": 10, "mux_t": 20, "fixup_iteration": 78, "shift_forward": 2,
                 "forward_flags": 3, "dirac_flags": 3, "projector": 4, "prepare_controlled_h": 4},
    "fanout": {"direction_flags": 1, "backward_flags": 1, "shift_backward": 2, "read_iteration": 78,
               "link_copy": 81, "mux_orientation": 2, "mux_t": 20, "fixup_iteration": 78, "shift_forward": 2,
               "forward_flags": 1, "dirac_flags": 1, "projector": 3, "prepare_controlled_h": 4},
}
FANOUT_ANCILLA_2033 = 5          # CNOT fan-out of the shared leaf / orientation / Dirac control (factory analysis)
DEPTH_SRC_2033 = ("r25 factory analysis (apply_log/r25_ch09.md, section 'T-depth'): reaction-limited layers per "
                  "item of single_particle_query_2O, lift_cost, the reflection and the Higgs vertex")


def sp_query_depth_2O(t_rot: float, variant: str) -> tuple[float, float]:
    """(lo, hi) T-depth of one query. 'register' lo: the itemized layers with 3 rotation layers (PREPARE's two rotations
    and the two commuting halves of the controlled QSP phase run in parallel on idle select ancilla); hi: every Toffoli,
    direct T and rotation serial (598 + 24 + 6 t_rot). 'fanout': the fan-out layers with 3 rotation layers (one value)."""
    if variant == "fanout":
        d = sum(SP_QUERY_DEPTH_2O["fanout"].values()) + 3 * t_rot          # 274 + 3 t_rot = 353.4 at 26.476
        return d, d
    lo = sum(SP_QUERY_DEPTH_2O["register"].values()) + 3 * t_rot           # 614 + 3 t_rot = 693.4
    hi = 598 + 24 + 6 * t_rot                                              # 780.9
    return lo, hi


def application_depth_2O(n_queries: int, t_rot: float, variant: str, lift: dict, outer: dict,
                         n_higgs_sites: int, higgs_select_toffolis: float, umult_toffolis: int,
                         umult_concurrent: int) -> tuple[float, float]:
    """(lo, hi) T-depth of one application: n_queries sequential queries (QSP is a sequential product) + the lift +
    the QSP reflection + the Higgs vertex. Lift: two serial unary iterations (lift['toffoli']) plus the uniform-over-3
    preparations, in parallel across coordinates (lo: 4 T-layers + 2 rotations) or serial (hi: 12 + 6 rotations).
    Reflection: an n-controlled flag on idle ancilla, log-tree (lo 5) or ladder (hi = its Toffoli count). Higgs: the
    SELECT walk is serial (higgs_select_toffolis); the U_x run umult_concurrent at a time on the clean ancilla free
    outside a query (lo: ceil(sites / concurrent) rounds of umult_toffolis) or all serial (hi)."""
    q_lo, q_hi = sp_query_depth_2O(t_rot, variant)
    lift_lo = lift["toffoli"] + 4 + 2 * t_rot
    lift_hi = lift["toffoli"] + lift["t_direct"] + lift["n_rot"] * t_rot
    ref_lo, ref_hi = 5 + outer["n_rot"] * t_rot, outer["toffoli"] + outer["n_rot"] * t_rot
    rounds = math.ceil(n_higgs_sites / umult_concurrent)
    hig_lo = higgs_select_toffolis + rounds * umult_toffolis
    hig_hi = higgs_select_toffolis + n_higgs_sites * umult_toffolis
    return (n_queries * q_lo + lift_lo + ref_lo + hig_lo, n_queries * q_hi + lift_hi + ref_hi + hig_hi)


def hwp_step_depth_2028(groups: list[dict], spatial_concurrency: int, rus_ancilla: int) -> dict:
    """(lo, hi) T-depth and peak ancilla of one 2028 Trotter step under Hamming-weight phasing (r25, R10).

    Slots run one after another. A slot is one phasing group, or up to `spatial_concurrency` spatial-hop groups on
    mutually non-adjacent links (disjoint sites, so their strings commute and one Clifford maps all of them in place to
    single-qubit Z's). In a slot all weight-bit rotations run in parallel, each on its own RUS ancilla (one rotation
    layer, t_rot deep), and the adders of the slot's groups run side by side: Toffoli depth ceil(2 log2 k) for a
    carry-save tree (lo) to k - w(k) serial (hi). Peak ancilla = max over slots of the adders' k - w(k) plus one RUS
    ancilla per concurrent rotation."""
    lo = hi = 0.0
    peak = 0
    n_slots = 0
    for g in groups:
        k, n = g["k"], g["n_groups"]
        c = spatial_concurrency if g["name"] == "spatial_hops" else 1
        slots = math.ceil(n / c)
        lo += slots * (g["t_rot"] + min(math.ceil(2 * math.log2(k)), hwp_toffolis(k)))
        hi += slots * (g["t_rot"] + hwp_toffolis(k))
        peak = max(peak, min(c, n) * (hwp_ancilla(k) + rus_ancilla * hwp_synth_rotations(k)))
        n_slots += slots
    return {"lo": lo, "hi": hi, "peak_ancilla": peak, "n_slots": n_slots}


# --------------------------------------------------------------------------- #
# Assumptions
# --------------------------------------------------------------------------- #

@dataclass(frozen=True)
class Assumptions:
    # ---- shared -----------------------------------------------------------
    # R-TOL (H. Lamm 2026-09-29): total synthesis error per shot, report-wide; each circuit sets its
    # own tolerance eps_rot = sqrt(eps_syn / N_rot) (common.eps_rot_for), N_rot = its synthesized
    # rotations in one shot.  The chapter keeps its own formula (full RUS fit); only eps changes.
    eps_syn: Tagged = Cited(EPS_SYN, EPS_SYN_SRC, "app06:108 'the total synthesis error per shot is 1e-2, and the "
                                                  "tolerance is set from the circuit's rotation count'")
    rot_synthesis: Tagged = Cited("rus", SYNTHESIS_SRC, "app06:108 'from the fit, 1.15 log2(1/eps) + 9.2' "
                                                        "(repeat-until-success, randomized); the chapter's rotations "
                                                        "other than the R-HOP hop")
    # RETIRED by R-TOL (not used by this model).  Kept as records: Ch. 6's ch9_rules_step and its
    # R-HOP cross-check read them (estimates/ch06_collider.py), and the pre-round-B / pre-R-HOP
    # records below are priced at them.
    t_per_rot: Tagged = Cited(30, SYNTHESIS_SRC,
                              "RETIRED (R-TOL 2026-09-29): app06:108 was '~30 T-gates per rotation (~25 from the "
                              "fit, rounded up)' at eps_gs ~ 1e-4; record only")
    eps_rot: Tagged = Assumed(1e-4, "RETIRED (R-TOL 2026-09-29): app06:108 was 'eps_gs ~ 1e-4'; record only")
    toffoli_convention: Tagged = Cited("textbook", "main-overview:47 (ruling R5, H. Lamm 2026-09-28)",
                                       "'it uses 7 T per Toffoli throughout'; app06:99 'with every Toffoli priced "
                                       "at 7 T'; was the 4-T measurement-assisted Toffoli of arXiv:1709.06648")
    t_gate_s: Tagged = Assumed(1e-6, "app06 boxes 'at 1 us/T-gate'")
    shot_overhead_s: Tagged = Assumed(1e-4, "per-shot overhead ~0.1 ms (register initialization, final readout, decode), "
                                            "report rule as Chs. 8 and 10; moves no printed number here")
    envelope_2028: Tagged = Cited(1e5, "DOE_RFI_2026", "app06:105 '<= 1e5 T-gates'")
    envelope_2033: Tagged = Cited(1e9, "DOE_RFI_2026", "app06:102 '~1e9 2033 resource limit'")
    faults_per_shot: Tagged = Assumed(0.1, "ruling R3 (TRACKED_CHANGES.md 2026-09-28): eps_l = 0.1 expected "
                                           "faults per shot, report-wide; app06:103 'at 0.1 expected faults per shot'")

    # ---- 2028: 1+1D Z3 domain-wall ---------------------------------------
    gauge_group_2028: Tagged = Assumed("Z3", "app06:105 'Z3-digitized 1+1D toy'")
    D_2028: Tagged = Assumed(1, "app06:105 '1+1D'")
    L_2028: Tagged = Assumed(8, "app06:105 'L=8'")
    L5_2028: Tagged = Assumed(4, "app06:105 'truncated L5=4'")
    Nf_2028: Tagged = Assumed(1, "app06:105 'N_f=1'")
    Nc_2028: Tagged = Assumed(3, "app06:105 'three colors ... decoupled identical copies'")
    ancilla_2028_asserted: Tagged = Stated(40, "app06:145 (pre-r20)", "RECORD (E23 (5)): the asserted '40 ancilla for the "
                                                                   "adiabatic preparation, measurement and Hamming-weight "
                                                                   "phasing'; replaced by the itemized peak (39 at r20, "
                                                                   "32 at r22)")
    anc_rus: Tagged = Cited(RUS_ANCILLA, "bocharovRoettelerSvore2015; ruling E26 (r22, H. Lamm 2026-10-01); common.RUS_ANCILLA",
                            "repeat-until-success synthesis holds ONE ancilla while a rotation is in flight ('Clifford+T "
                            "basis, plus one ancilla qubit and measurement'); the weight-bit rotations of a phasing group "
                            "are synthesized one at a time. Was 2 in both eras (r20)")
    n_pauli_per_hop: Tagged = Stated(2, "app06:108", "'two Pauli strings per hopping bilinear'")
    n_pauli_per_electric: Tagged = Stated(
        3, "app06:108", "'three diagonal Z-strings per two-qubit link, each grouped across the eight links'")
    fifth_boundary: Tagged = Cited(
        "open", "Kaplan_chiral_lattice_fermions, arxiv_2505_20419; r21 (3a), H. Lamm 2026-10-01",
        "DERIVED HERE: domain-wall fermions have open boundaries in the fifth direction. The chiral modes sit on the "
        "walls s = 1 and s = L5, and the only term joining the walls is the quark mass m_f, taken as zero in this benchmark (the chapter states no m_f; "
        "at m_f != 0 the V wall-joining terms are two more groups of k = 24, +510 T per step). "
        "So there are V (L5 - 1) = 24 fifth-direction hops (two groups of k = 72), not the V L5 = 32 the chapter "
        "printed as 'the same number of fifth-direction hops'. 'periodic' reproduces the pre-r21 count")
    spatial_hop_grouping: Tagged = Stated(
        "color+s", "app06:108",
        "'The twelve color and s copies of each spatial hop share one s-independent link and form one group "
        "per link and Pauli string: $64$ groups of $k=12$' (ruling R6, proposal ch09-2028-grouping 'grouping C'). "
        "'color' is the pre-R6 sentence (k=Nc per hop), 'color+site' also groups the sites (k=Nc*V*L5)")
    hop_synthesis: Tagged = Cited("rus", "r17 ruling (e), H. Lamm 2026-10-01 (apply_log/r17_core.md, r17_ch09.md); "
                                         "FermionPrimitives_unpub section_resources.tex:7-10",
                                  "the Z3 spatial hop's rotations at the full fit 1.15 log2(1/eps) + 9.2, the chapter's own "
                                  "price and the draft's convention for every fermion rotation; replaces R-HOP's slope-only "
                                  "1.15 log2(1/eps) (author's 'do 1', 2026-09-29), now a legacy record. The 8 strings per "
                                  "copy and C_W = 0 still come from the R-HOP structure (Z3 has no draft entry; "
                                  "groups.hop_link_cost('Z3', ...) uses the same structure)")
    hop_synthesis_legacy: Tagged = Cited("rus-slope", "TRACKED_CHANGES.md 'RULING, R-HOP report-wide' (H. Lamm 2026-09-29)",
                                         "RETIRED by r17 ruling (e): R-HOP's slope-only price, kept so the 2.0e5 record "
                                         "stays reproducible")
    onsite_term_2028: Tagged = Cited(
        "wilson", "arxiv_2505_20419; FermionPrimitives_unpub section_mass.tex:12",
        "r17 ruling (f), H. Lamm 2026-10-01: the on-site part of h_DW = gamma^0 D_W, coefficient m + r(1+d) "
        "(= m + 2 at d=1, r=1; non-zero across the topological window -2 < m < 0), missing before r17. "
        "2N R_Z per site at d=1 (groups.wilson_mass_rotations; N = Nc Nf), on each of the V L5 sites: 192 equal-|angle| "
        "Z rotations in the draft's diagonal-gamma^0 encoding, one lattice-wide phasing group (signs are Clifford). "
        "'none' reproduces the pre-r17 step")
    hwp_chunk_k: Tagged = Stated(32, "app06 2028 derivation", "'Splitting the k=72 and k=192 groups into groups of at "
                                                              "most k=32 needs 31 and fits' (E26: one ancilla per adder, "
                                                              "no separate weight register; was '37 ancilla')")
    n_trotter_2028: Tagged = Stated(20, "app06:109", "'We shorten the adiabatic ground-state ramp from N_Trotter ~ 30 "
                                                     "to N_Trotter = 20 steps' (ruling R6)")
    n_trotter_2028_before_r6: Tagged = Stated(30, "app06:109", "'from N_Trotter ~ 30'; the ramp length before R6, "
                                                               "emitted so the trade is visible")
    n_probe_steps: Tagged = Stated(1, "app06:109", "'With a single-time-step probe ... the shot is 21 steps'")
    t_per_step_2028_stated: Tagged = Stated((1.28e4, 1.37e4), "app06 2028 derivation",
                                            "'~1.28e4 T per step' unchunked, '~1.37e4 T per step' chunked to k <= 32 "
                                            "(r21 a: 24 fifth-direction hops; r17 1.31e4 / 1.42e4)")
    t_per_shot_2028_stated: Tagged = Stated((2.69e5, 2.89e5), "app06 2028 derivation",
                                            "'~2.69e5 T-gates unchunked and ~2.89e5 chunked' (r21 a; r17 2.76-2.97e5)")
    ratio_to_envelope_2028_stated: Tagged = Stated((2.7, 2.9), "app06 2028 derivation",
                                                   "'a factor 2.7-2.9 above the first-generation 1e5 envelope'")
    t_per_rot_coherent_stated: Tagged = Stated(31, "app06:109",
                                              "'per-rotation cost would rise to ~31 T' (R-TOL: the same 1e-2 "
                                              "total summed linearly; was ~40 at 1e-4)")
    t_per_step_coherent_stated: Tagged = Stated(1.60e4, "app06 2028 derivation", "'the unchunked step to 290 x 31.3 + "
                                                                                  "991 x 7 ~ 1.60e4 T' (r21 a; r17 1.64e4)")
    t_per_shot_coherent_stated: Tagged = Stated(3.4e5, "app06 2028 derivation", "'and the shot to ~3.4e5 T' (3.363e5)")
    hwp_ancilla_unchunked_stated: Tagged = Stated(190, "app06 2028 derivation", "'The k=192 group holds 190' (E26, r22: "
                                                                                "k - w(k); was 'about 200' with the weight "
                                                                                "register counted twice)")
    hwp_ancilla_chunked_stated: Tagged = Stated(31, "app06 2028 derivation", "'needs 31 and fits' (E26, r22; was 37)")
    lq_unchunked_stated: Tagged = Stated(400, "app06 2028 derivation", "'the unchunked end needs 191, a register of about "
                                                                       "400 LQ' (192 + 16 + 190 + 1 = 399; was ~410)")
    lq_single_copy_stated: Tagged = Stated(110, "app06 2028 derivation", "'rather than the ~110 LQ a single copy would need' "
                                                                         "(64 + 16 + 32 = 112; was ~120)")
    eps_l_2028_stated: Tagged = Stated(3.5e-7, "app06 2028 derivation and box", "'eps_l <~ 3.5e-7' at 0.1 expected faults "
                                                                                "per shot (R3), from the chunked 2.885e5 T")
    # r21 (3d), DERIVED HERE: the 2028 shot count. The chapter's sigma^2 ~ 10 is for the condensate averaged over the
    # 27 sites of the 2033 lattice. The 2028 estimator averages N_c V = 24 site-color copies, the three color copies
    # being decoupled and TAKEN AS independent (they share the link, so this is an assumption), so sigma^2 = 10 x 27 / 24 = 11.25 and shots = sigma^2 / eps^2 = 1125.
    # Operator-norm bound: each copy's bilinear has eigenvalues in {-1, 0, 1}, so its variance is at most 1; if the
    # sites of one color copy were fully correlated only the 3 color copies would average: sigma^2 <= 10 x 27 / 3 = 90.
    shots_2028_stated: Tagged = Stated(1.1e3, "app06 2028 box", "'1.1e3 (psi-bar psi at one coupling)'; 1125 from "
                                                                "sigma^2 = 11.25 (r21 d; was 1e3, borrowed from V = 27)")
    shots_2028_worst_stated: Tagged = Stated(9e3, "app06 shot-count paragraph", "'at most sigma^2 ~ 90, 9e3 shots' if the "
                                                                                "sites of a color copy are fully correlated")
    wall_per_shot_2028_stated: Tagged = Stated((0.27, 0.29), "app06 2028 derivation", "'0.27-0.29 s per shot at 1 us per "
                                                                                      "T-gate' (r21 a; r17 0.28-0.30)")
    # the 1+1D overlap figures that motivate the domain-wall choice. r22 (E26 (i)): priced in the valid single-particle
    # architecture, as the 2033 query is (single_particle_query_Z3_1p1d, lift_cost); the r21 (3e) second-quantized
    # pricing (240 Pauli strings + 25 rotations) is a record.
    kappa_1p1d_a: Tagged = Assumed(30, "app06 2028 derivation 'even at kappa ~ 30'")
    delta_1p1d_a: Tagged = Assumed(1e-2, "app06 2028 derivation 'delta_sgn ~ 1e-2'")
    c_sp_1p1d_stated: Tagged = Stated(4.0e2, "app06 2028 derivation", "'a query is 40 Toffolis and 6 synthesized rotations, "
                                                                      "~4.0e2 T' (393-406 over the four circuits)")
    sgn_T_1p1d_a: Tagged = Stated(5.6e4, "app06 2028 derivation", "'one application of the block-encoded overlap "
                                                                  "Hamiltonian costs 5.6e4 T' (140 queries + lift; r21 "
                                                                  "3.0e5 in the second-quantized pricing)")
    gs_prep_1p1d_stated: Tagged = Stated(5.8e5, "app06 2028 derivation", "'at least O(10) of them, ~5.8e5 T' (r21 3.1e6)")
    gs_prep_over_dwf_stated: Tagged = Stated(2, "app06 2028 derivation", "'about six times the budget and twice the "
                                                                         "domain-wall shot' (5.78e5 / 2.885e5 = 2.0; r21 "
                                                                         "'eleven')")
    gs_prep_over_envelope_stated: Tagged = Stated(6, "app06 2028 derivation", "'about six times the budget' (5.8)")
    kappa_1p1d_b: Tagged = Assumed(50, "app06 algorithm requirements 'at kappa >~ 50'")
    sgn_T_1p1d_b: Tagged = Stated(9e4, "app06 algorithm requirements", "'already costs ~9e4 T' (9.28e4; r21 ~5e5)")
    gs_prep_1p1d_b_stated: Tagged = Stated(1e6, "app06 algorithm requirements", "'~1e6 T, an order of magnitude past the "
                                                                                "1e5 envelope' (9.55e5)")
    n_sgn_gs_prep: Tagged = Stated(10, "app06 algorithm requirements", "'ground-state preparation needs at least O(10) of "
                                                                       "them' (a floor under the lifted normalization "
                                                                       "2Q = 96)")
    c_be_1p1d_sq_stated: Tagged = Stated(2.2e3, "app06 (r21)", "RECORD (second-quantized pricing, r21 e): 'c_be ~ 2.2e3 T per "
                                                               "query' (240 Toffolis + 25 rotations)")

    # ---- 2033: 3D 3^3 SU(2)/2O overlap + Higgs ---------------------------
    gauge_group_2033: Tagged = Assumed("2O", "app06 '|G|=48 [2O binary octahedral]'")
    D_2033: Tagged = Assumed(3, "app06 'D=3'")
    L_2033: Tagged = Assumed(3, "app06 'V=3^3=27'")
    Nf_2033: Tagged = Assumed(1, "app06 'N_f=1' (Witten anomaly minimum)")
    Nc_2033: Tagged = Assumed(2, "app06 'N_c=2 [SU(2)]'")
    N_rho: Tagged = Assumed(16, "app06 'N_rho=16'")
    ancilla_2033_asserted: Tagged = Stated(110, "app06 (pre-r20)", "RECORD (E23 (5)): '+ ~110 ancilla', asserted; "
                                                                   "replaced by the itemized peak (43 at r20, 33 at r21, "
                                                                   "38 at r22 in the single-particle architecture)")
    # r22 (E26 (i)): the sign function acts on the single-particle kernel (arxiv_2607_28524); compiled here
    sgn_architecture: Tagged = Cited(
        "single_particle_lift", "arxiv_2607_28524; ruling E26 (i), H. Lamm 2026-10-01 ('do whatever is easiest and cheaper')",
        "the Hamiltonian is the bilinear psi^dag h_ov[U] psi with h_ov = gamma^0 + sgn(H_W[U]); sgn acts on the "
        "single-particle kernel, a Q x Q matrix of link operators, promoted to a 9-qubit index register and the link "
        "registers, and two selected-Majorana unitaries lift it to the Fock register. The r21 construction (QSVT of sgn "
        "on a second-quantized Jordan-Wigner block encoding of H_W) is INVALID: it computes the sign of the many-body "
        "operator psi^dag H_W psi, a different Fock-space operator with an exact zero mode (undefined kappa). It is kept "
        "as a record only (intermediates 'second_quantized_record')")
    m0_wilson: Tagged = Assumed(1.5, "r22: the negative Wilson mass -m_0 of the overlap kernel, inside the doubler window "
                                     "0 < m_0 < 2 (the chapter states none). Used only for the normalization "
                                     "alpha = 2D + |D - m_0| = 7.5 of the block encoding of H_W at r = 1 and for the "
                                     "record of the other kappa convention; moves no priced number")
    eps_qsp: Tagged = Assumed(1e-2, "r22: total QSP (Jacobi-Anger truncation) error of the kernel evolution over one "
                                    "shot; each Trotter step gets eps_qsp / N. Same size as the synthesis budget")
    # E23 (5), r20; r22: 2033 ancilla itemized from the circuits as priced, peak of the schedule.
    # r22 (single-particle architecture): 'stack' = allocated for the whole shot; 'application' = held through one
    # application of the block-encoded psi^dag h_ov psi; 'query' = the term register of one H_W query; 'select' = the
    # link read inside SELECT.
    anc_readout: Tagged = Assumed(1, "E23 (5), r20; r22 wording (validation): control qubit of the Ramsey / Hadamard-test "
                                     "readout of the unequal-time Chern-Simons correlator (arxiv_2208_13112). It acts on "
                                     "the two inserted operators only and is used after the state is prepared, so it is "
                                     "allocated for the whole shot but idle during the TPQ filter (one spare at the peak)")
    anc_qsp_signal: Tagged = Assumed(1, "E23 (5), r20; r22: one QSP signal qubit, shared by the TPQ filter and the kernel "
                                        "evolution (a real polynomial needs one signal qubit, gilyen2018QSingValTransfArXiv)")
    anc_fpaa: Tagged = Assumed(1, "E23 (5), r20: marker qubit of the fixed-point amplitude amplification (used at the "
                                  "amplification reflections only)")
    anc_h_lcu: Tagged = Assumed(2, "E23 (5), r20: LCU index over the parts of H (overlap fermion part, gauge, Higgs)")
    anc_dov_lcu: Tagged = Assumed(1, "E23 (5), r20: LCU qubit of h_ov = gamma^0 + sgn(H_W) (the chapter's "
                                     "D_ov = 1 + gamma5 sgn(H_W))")
    anc_sgn_signal: Tagged = Assumed(1, "E23 (5), r20: QSP signal qubit of sgn(H_W)")
    anc_sgn_real_part: Tagged = Assumed(1, "r22: the qubit that selects the two phase sequences +-Phi whose average is the "
                                           "real sign polynomial (gilyen2018QSingValTransfArXiv)")
    anc_lift_type: Tagged = Cited(2, "arxiv_2607_28524 (Supplement theorem)", "the two-qubit register A of the lift: four "
                                                                             "Majorana types (X or Y on each side)")
    fpaa_reflection_ancilla: Tagged = Cited(2, "arxiv_2407_17966", "r22 (validation tier A): the amplification reflects about "
                                                                   "the reference state on the whole register; an "
                                                                   "n-controlled NOT costs 2n Toffolis on two clean "
                                                                   "ancilla, borrowed from qubits idle between queries. "
                                                                   "NOT PRICED in the totals (record: about 2e5 T per shot)")
    # RECORDS of the retired second-quantized construction (r20-r21): its stack, block-encoding and link-term registers
    anc_unary_control: Tagged = Assumed(1, "RECORD (second-quantized construction, r20): the control qubit of SELECT's unary "
                                           "iteration (Babbush_PRX_2018); the AND ladder holds ceil(log2 L) - 1 temporaries")
    anc_link_pq: Tagged = Cited(6, "FermionPrimitives_unpub section_su2_diag.tex:188 (BO row, 'Ancilla 6')",
                                "RECORD (second-quantized construction, draft frame): the p, q angle registers of the 2O "
                                "color squish, held through V, hop, V^dag (share='link')")
    anc_link_hop_register: Tagged = Cited(3, "FermionPrimitives_unpub section_hopping.tex:154",
                                          "RECORD (second-quantized construction, draft frame): 'at most a 3-qubit ancilla "
                                          "register for the eigenvalue computations', held with p, q")
    anc_link_parity: Tagged = Cited(1, "FermionPrimitives_unpub section_su2_diag.tex:206",
                                    "RECORD (second-quantized construction, draft frame): 'compute the parity of the p,q "
                                    "register on to another clean ancilla'")
    anc_spinor: Tagged = Cited(2, "FermionPrimitives_unpub section_spin_diag.tex:280",
                               "RECORD (second-quantized construction): 'requires 2 spare clean ancilla' (the d=3 spinor frame)")
    L5_backup: Tagged = Assumed((8, 16), "app06 route selection 'L5 factor (8-16 for the back-up here, 16-32 in production)'")
    kappa_2033: Tagged = Assumed(30, "r24 (H. Lamm 2026-10-02, 'lower kappa is ok'): app06 'kappa = 30 (a heavier, "
                                     "coarse-mass regime)'; was 1e2 (r17-r23). Since r22 kappa = alpha / Delta, alpha the "
                                     "normalization of the block encoding of H_W (the convention of arxiv_2607_28524); "
                                     "the r21 text had ||H_W|| / Delta (record d_sgn_if_kappa_is_norm_over_gap)")
    kappa_2033_r22: Tagged = Assumed(1e2, "RECORD (r17-r23 headline kappa ~ 1e2): prices the retired second-quantized "
                                          "construction and the r22 headline (2.008e11), so every earlier print stays "
                                          "reproducible; also the kappa ~ 1e2 sensitivity of the prose")
    delta_sgn: Tagged = Assumed(1e-3, "app06 'delta_sgn ~ 1e-3'; the box budget row prints '~0.1% (1e-3)' (ch09-delta-sgn)")
    d_sgn_stated: Tagged = Stated(208, "app06 2033 derivation and box", "'d_sgn = ceil(kappa ln delta^-1) = 208' at kappa = 30 "
                                                                        "(r24; r17-r23 691, printed ~700)")
    # the Higgs vertex's SELECT strings (r21 c), and the LCU records of the second-quantized construction
    lcu_strings_per_bilinear: Tagged = Cited(2, "Jordan-Wigner encoding (standard)",
                                             "a color-diagonal bilinear a^dag b + h.c. is (XX + YY)/2 with a Z string: two "
                                             "Pauli strings (record: the hop strings of the second-quantized LCU)")
    lcu_toffoli_per_term: Tagged = Stated(1, "app06 2033 derivation; Babbush_PRX_2018",
                                          "one SELECT Toffoli per LCU term (unary iteration over L terms costs L - 1): the "
                                          "40 Higgs strings per site of the headline, and the second-quantized record")
    lcu_terms_stated_pre_r21: Tagged = Stated(1.7e3, "app06 (pre-r21)", "RECORD: '~1.7e3-term LCU', an absolute at V=27. "
                                                                        "8Q = 1728 is exactly the 4Q + 4Q selected strings "
                                                                        "of U_L and U_R in arxiv_2607_28524: the old line "
                                                                        "was the lift, paid once per application")
    lcu_terms_per_site_stated: Tagged = Stated(72, "app06 (r21)", "RECORD (second-quantized construction): '72 Pauli strings "
                                                                  "per site': 24 hop + 8 on-site + 40 Higgs")
    lcu_terms_stated: Tagged = Stated(1.9e3, "app06 (r21)", "RECORD (second-quantized construction): '1944 at V = 27'")
    higgs_strings_per_site_stated: Tagged = Stated(40, "app06 2033 derivation", "'$40$ select strings per site' of the Higgs "
                                                                                "vertex (8 on-site Z strings x the 5-term "
                                                                                "binary expansion of rho)")
    n_umult_per_higgs_site: Tagged = Stated(
        1, "app06 2033 derivation", "'27 compiled group multiplications U_x' (ruling VERTEX, kept for the Higgs; the cost is "
                                    "read from GROUPS['2O'].primitives['U_mul'], 56 Toffolis = 392 T). r22: the Higgs vertex "
                                    "is not part of H_W; it sits outside the sign function, once per application")
    # RECORD parameters of the second-quantized construction (r17-r21): the link frame and the draft's compiled hop
    link_frame: Tagged = Cited("multiplexer", "r21 ruling (1), H. Lamm 2026-10-01 ('promote'); DERIVED HERE r20 "
                                              "(colour_multiplexer_2O; arxiv_2312_10285 ordered-product encoding)",
                               "RECORD (second-quantized construction): the link-dressed Fock-space hop as the 2O color "
                               "multiplexer W and W^dag on the far site plus the draft's spinor frame, 468 T per link per "
                               "query; 'draft' = the general diagonalizing frame of the unpublished counts")
    hop_fermion_2033: Tagged = Stated("wilson", "app06 2033 derivation", "the overlap kernel H_W is the Wilson operator")
    hop_frame_undo: Tagged = Cited(True, "ruling E21 (1), H. Lamm 2026-10-01 (apply_log/r19_core.md); "
                                         "FermionPrimitives_unpub section_hamiltonian.tex:86-91, su2_diag.tex:238, "
                                         "section_resources.tex:15",
                                   "RECORD (second-quantized, draft frame). OUR ESTIMATE, ruled in: the color frame V_g is "
                                   "applied AND undone on both sites of every link in every BE query. The draft draws V_C^dag "
                                   "(fig. da_diagonalizer) but does not count it: its 2 C^G is V_C on the two sites only")
    hop_share: Tagged = Cited("link", "ruling E21 (1), H. Lamm 2026-10-01 ('make the switch'; apply_log/r19_core.md)",
                              "RECORD (second-quantized, draft frame): the color squish and parity flags are computed once per "
                              "link and held through V(x), V(y), hop, V^dag(x), V^dag(y); 'draft' recomputes them")
    hop_mcx: Tagged = Cited("mbu", "groups.mcx_toffolis (r17 core)", "C^nX ladders at n-1 Toffolis: the n-2 AND temporaries "
                                                                    "uncomputed by measurement (the report's rule)")
    hop_n_spin: Tagged = Cited(2, "ruling E21 (2), H. Lamm 2026-10-01 (author-confirmed); FermionPrimitives_unpub "
                                  "section_hamiltonian.tex:77-80",
                               "spinor components that hop after the spinor frame, n_D/2 = 2 at d=3, r=1. Sets the hop norm "
                               "N_c n_spin = 4 of the worst-case Trotter bound, and the second-quantized record's link term")
    hop_phasing_be: Tagged = Stated("none", "app06 (r21)", "RECORD (second-quantized, draft frame): a block-encoding query "
                                                           "carries no Trotter phasing")
    w2_cs_multiplicity: Tagged = Stated(4, "FermionPrimitives_unpub section_resources.tex:50-53",
                                        "RECORD (the r17-r20 band top): the draft's 3d/4d Wilson hop 4 C^S + 148 N per link")
    w2_spinor_t_per_colour: Tagged = Stated(148, "FermionPrimitives_unpub section_resources.tex:50-53", "RECORD: '+148N'")
    n_rot_be: Tagged = Stated(25, "app06 (r21)", "RECORD (second-quantized construction): '~25 synthesized rotations in "
                                                 "PREPARE and the QSP phases' per query")
    c_be_stated: Tagged = Stated(6.3e4, "app06 (r21)", "RECORD (second-quantized construction): 'c_be ~ 6.3e4 T-gates' (62,896)")
    c_be_draft_stated: Tagged = Stated(2.4e6, "app06 (r21)", "RECORD: 'c_be ~ 2.4e6' with the draft's diagonalizing frame")
    t_shot_draft_stated: Tagged = Stated(7.0e14, "app06 (r21)", "RECORD: '7.0e14 T per shot' with the draft's frame (7.007e14)")
    draft_over_mux_stated: Tagged = Stated(39, "app06 (r21)", "RECORD: '39 times the estimate we quote' (38.7)")
    draft_over_mux_link_stated: Tagged = Stated(64, "app06 (r21)", "RECORD: '64 times the multiplexer's 468 T' (63.6)")
    clifford_2O: Tagged = Cited("clifford", "standard group theory (E23 (1), H. Lamm 2026-10-01: 'a math fact')",
                                "2O is the lift to SU(2) of the octahedral rotation group O = single-qubit Clifford group "
                                "mod phases (24 elements); every 2O element normalizes the Pauli group and has entries in "
                                "Z[1/sqrt2, i] (checked in the tests)")
    t_per_dov_sq_stated: Tagged = Stated(4.3e7, "app06 (r21)", "RECORD (second-quantized construction): '691 x 6.3e4 ~ 4.3e7 T'")
    t_per_shot_sq_stated: Tagged = Stated(1.8e13, "app06 (r21 box)", "RECORD (second-quantized construction at the worst-case "
                                                                     "step): '~1.8e13' (1.809e13), 1005 LQ")
    # the single-particle query and one application, as the prose states them (r22)
    sp_query_toffoli_stated: Tagged = Stated(598, "app06 2033 derivation", "'561 of the query's 598 Toffolis'; the query "
                                                                           "is 598 Toffolis, 24 T and 6 rotations")
    c_sp_stated: Tagged = Stated(4.4e3, "app06 2033 derivation", "'598 x 7 + 24 + 6 x 29.1 ~ 4.4e3 T' per query (4384.8)")
    t_per_dov_stated: Tagged = Stated(1.8e6, "app06 2033 derivation", "'417 queries per application and 1.8e6 T' (r25 (a), "
                                                                      "Yukawa inside: 1,847,747.1; r24 9.4e5; r22 3.1e6)")
    t_evol_fm: Tagged = Assumed(2, "r24 (H. Lamm 2026-10-02, 'ok to run for a shorter amount of time'): app06 't = 2 fm "
                                   "= 10a'; was 10 fm = 50a (r17-r23)")
    t_evol_fm_r22: Tagged = Assumed(10, "RECORD (r17-r23 evolution time 10 fm = 50a): prices the retired records and the "
                                        "r22 headline")
    a_fm: Tagged = Assumed(0.2, "app06 2033 box 'a ~ 0.2 fm'")
    eps_trotter: Tagged = Assumed(0.1, "app06 'at eps ~ 0.1'")
    # r22 (E26 (ii)): the step count is a state-dependent estimate; the worst-case bound (r21 3b) is quoted as the worst case
    trotter_rule_2033: Tagged = Cited(
        "state", "ruling E26 (ii), H. Lamm 2026-10-01 ('i dont think we can do such a wildly off result'); "
                 "AlvesLammLiu_inprep (FERMILAB-PUB-26-0397-T), evaluation DERIVED HERE",
        "N_Trotter from the state-dependent second-order estimate in the thermal state (trotter_state_dependent_steps): "
        "norms of the nested commutators replaced by their standard deviations, weak-coupling evaluation, 5031 steps. "
        "Records: 'bound' = the worst-case commutator bound (childs2021theory; 83,195, the r21 headline), 'step' = the "
        "fixed step of Ch. 5 (1000), 'stated' = the pre-r21 ~30")
    temp_lat_2033: Tagged = Assumed(0.5, "r22: temperature in lattice units, T a, for the thermal-state variances. The chapter "
                                         "says T ~ T_c only; 0.15-1 moves the central count over 5008-5342")
    g2_bare: Tagged = Assumed(1.0, "r21 (3b): bare gauge coupling for the WORST-CASE bound. The chapter states none; an "
                                   "order-one coupling. The state-dependent estimate does not depend on it at leading order")
    casimir_max_2O: Tagged = Assumed(3.75, "r21 (3b): largest SU(2) Casimir j(j+1) kept by the 2O electric term, j = 3/2; "
                                           "enters the worst-case bound only")
    n_trotter_2033_pre_r21: Tagged = Stated(30, "app06 (pre-r21)", "RECORD: 'N_Trotter ~ 30 steps', dt = 1.7 a. With a valid "
                                                                   "block encoding it would need 752 applications per step")
    n_trotter_2033_stated: Tagged = Stated(450, "app06 2033 derivation and box", "'N_Trotter ~ 450 steps of dt ~ 0.022 a' "
                                                                                   "(r24 at t = 10a; r22 5031 at 50a)")
    n_trotter_range_stated: Tagged = Stated((1.1e2, 1.3e3), "app06 2033 derivation", "'a range of 1.1e2-1.3e3 steps' (r24: 107, 1274) "
                                                                                     "(1186 fastest-mode phase; 14,242 "
                                                                                     "with the mean kept)")
    lambda_sd_stated: Tagged = Stated(20, "app06 2033 derivation", "'Lambda_sd ~ 20 a^-3' at T a = 0.5 (20.25)")
    n_modes_stated: Tagged = Stated(156, "app06 2033 derivation", "'156 transverse oscillators (three colors, two "
                                                                  "polarizations, 26 nonzero momenta)'")
    n_trotter_bound_stated: Tagged = Stated(7.4e3, "app06 2033 derivation", "'7.4e3 steps' of the worst-case bound (r24: 7442 "
                                                                            "at t = 10a; r21-r23 83,195 at 50a)")
    lambda_2033_stated: Tagged = Stated(5.5e3, "app06 2033 derivation", "'Lambda ~ 5.5e3 a^-3' of the worst-case bound "
                                                                        "(5537 at g^2 = 1)")
    dt_empirical_fm: Tagged = Cited(0.01, "app03 2033 box ('2nd-order Trotter at Delta t = a_t/10', a_t = 0.1 fm/c)",
                                    "the fixed step of Ch. 5's 3^3 evolution, the printed comparison")
    t_shot_empirical_step_stated: Tagged = Stated(1.5e10, "app06 2033 derivation", "'200 steps ... the shot would be "
                                                                                    "1.4e10 T' at the fixed step of Ch. 5 "
                                                                                    "(17 filter calls: 1.4989e10; C1: 1.3508e10; r25: 1.2635e10; r24 6.4e9; 33 per step)")
    n_dov_per_step_stated: Tagged = Stated(19, "app06 2033 derivation", "'19 applications per step': the Jacobi-Anger order "
                                                                        "at 2Q dt = 9.6 (r24; r22 13 at 4.3)")
    n_dov_per_step_worst_stated: Tagged = Stated(6, "app06 2033 derivation", "'7.4e3 steps of 6 applications each' at the "
                                                                             "worst-case step (r24)")
    n_dov_per_step_bound_stated: Tagged = Stated(5, "app06 (r22-r23)", "RECORD: 'with 5 applications per step' at the r22 worst-case "
                                                                             "step (derived there too); also the r21 "
                                                                             "'~5 D_ov applications per step' that the "
                                                                             "second-quantized record keeps")
    n_dov_evolution_stated: Tagged = Stated(8.6e3, "app06 2033 derivation", "'The evolution is then 8.6e3 applications' "
                                                                            "(r24: 8550; r22 65,403)")
    n_dov_floor_stated: Tagged = Stated(4.3e3, "app06 2033 derivation", "'no step size takes the evolution below 2Qt = 4320 "
                                                                        "applications, 4.1e9 T' (r24; r22 21,600)")
    t_evolution_stated: Tagged = Stated(1.6e10, "app06 2033 derivation", "'or 1.6e10 T-gates' (r25: 1.5798e10; r24 8.0e9)")
    t_shot_range_stated: Tagged = Stated((1.3e10, 2.9e10), "app06 2033 derivation", "'Over the range above it is "
                                                                                   "1.1-2.6e10 T' (17 filter calls: 1.3469e10, 2.8708e10; C1 1.2-2.7e10; r25: 1.1115e10, 2.6352e10; "
                                                                                   "r24 5.6e9-1.3e10)")
    t_shot_bound_stated: Tagged = Stated(8.5e10, "app06 2033 derivation", "'at the worst-case step ... 8.3e10 T' "
                                                                          "(17 filter calls: 8.546e10; C1 8.4e10; r25: 8.310e10; r24 4.2e10; r22 1.276e12)")
    beta_H: Tagged = Assumed(1e2, "RECORD (r25 and earlier): app06 'beta||H|| ~ 1e2'. Referee C1 (2026-10-04): inconsistent "
                                  "with the chapter's own block encoding, whose normalization is 2Q = 432 at beta = 2a")
    beta_H_rule: Tagged = Assumed("normalization", "referee C1 (2026-10-04), DERIVED HERE: the filter e^{-beta H/2} is a "
                                  "polynomial in H / 2Q, so beta||H|| = 2Q / (T a) = 864 ('normalization'); 'stated' "
                                  "keeps the r25 record 1e2. A direct Chebyshev fit of e^{-432 (1 + x)} on [-1, 1] to "
                                  "1e-4 needs degree 84 (scratch referee/ch09/deg.py); the sqrt formula gives 89.2")
    delta_beta: Tagged = Assumed(1e-4, "app06 'and delta_beta ~ 1e-4' (ch09-delta-beta); "
                                       "sqrt(1e2 ln 1e4) = 30.35 -> d_beta = 30; 1e-3 would give 26.3")
    d_beta_stated: Tagged = Stated(89, "app06 2033 derivation", "referee C1 (2026-10-04): 'degree d_beta = 89' at beta||H|| = 2Q / (T a) = 864 (89.2); r25 '~30' at the record 1e2")
    aa_rounds: Tagged = Stated(8, "app06 2033 derivation, route (iv)", "'~8 rounds of fixed-point amplitude amplification'; "
                                                                       "'an estimate rather than a derivation'")
    aa_call_rule: Tagged = Assumed("2k+1", "verifier 2026-10-04: k rounds of amplitude amplification call the filter "
                                   "2k + 1 times (U and U^dag once per round, plus the first call); 'per_round' keeps "
                                   "the r24-r25 and first-referee-pass record of one call per round")
    n_dov_tpq_stated: Tagged = Stated(1.5e3, "app06 2033 derivation", "'1.5e3 applications' (89 x 17 = 1513: 2k + 1 filter calls in k = 8 rounds, verifier 2026-10-04; first referee pass 712 = 89 x 8; r25 240)")
    t_tpq_stated: Tagged = Stated(2.8e9, "app06 2033 derivation", "'2.8e9 T' inside the sphaleron shot "
                                                                         "(2.796e9 at 17 calls; first referee pass 1.3e9; r25: 4.4346e8; r24 2.3e8; r22 7.3e8)")
    t_per_shot_2033_stated: Tagged = Stated(1.9e10, "app06 2033 derivation and box", "'per-shot total over the N_Dov = 8.8e3 "
                                                                                    "applications is 1.6e10 T-gates' (17 filter calls: 1.8597e10; C1 1.7115e10; r25: "
                                                                                    "1.6242e10; r24 8.2e9; r22 2.008e11)")
    ratio_representative_stated: Tagged = Stated(19, "app06 2033 derivation", "'a factor 19 above the budget' (17 filter calls: 18.60; C1 2026-10-04: 17.11; r25: 16.24; r24 8.2; r22 2.0e2)")
    gibbs_T_stated: Tagged = Stated((4.8e12, 1.5e14), "app06 route (iv)", "'The range is 5.5e11-1.7e13' (r25, re-priced at the "
                                                                          "Yukawa-inside application; r24 2.8e11-8.9e12; "
                                                                          "r22 9.2e11-2.9e13)")
    gibbs_T_central_stated: Tagged = Stated(2.5e13, "app06 route (iv)", "'2.5e13 T' (C1 2026-10-04, beta||H|| = 864: 2.518e13; r25: 2.914e12; r24 1.5e12; r22 4.8e12)")
    gibbs_over_tpq_stated: Tagged = Stated(9.0e3, "app06 route (iv)", "'9.0e3 times' (17 filter calls: 9005; C1 2026-10-04: 1.914e4; r25 6.6e3; both scale with "
                                                                      "the application price)")
    kappa_physical: Tagged = Stated((1e3, 1e4), "app06 2033 derivation", "'kappa ~ 1e3-1e4 ... near-physical-mass'")
    t_physical_stated: Tagged = Stated((6.1e11, 6.1e12), "app06 2033 derivation", "'6.1e11-6.1e12 T' (17 filter calls: 6.107e11, 6.120e12; C1: 5.620e11, 5.632e12; r25: 5.333e11, 5.345e12; "
                                                                                  "r24 2.7e11-2.7e12; r22 2.0e12-2.0e13)")
    eps_l_2033_stated: Tagged = Stated(5.4e-12, "app06 2033 derivation and box", "'eps_l <~ 5.4e-12' (17 filter calls: 5.377e-12; C1: 5.843e-12; r25 '6.2e-12') at 0.1 faults per shot "
                                                                                 "for the 10a sphaleron shot (R3; r25 6.157e-12; "
                                                                                 "r24 1.2e-11; r22 5.0e-13)")
    # thermal-state route alternatives (app06 route (iv)), emitted so the prose numbers are checkable
    haar_deficit_frac: Tagged = Assumed(0.05, "app06 route (iv) 'even a 5% deficit'")
    haar_register_qubits: Tagged = Assumed(1e3, "app06 route (iv) 'on a 1e3-qubit register'")
    haar_rounds_stated: Tagged = Stated(3e7, "app06 route (iv), RECORD", "'~3e7 amplification rounds' (ch09-haar-rounds; was ~1e7); "
                                        "retired from the prose 2026-10-05 (typical_reference_tpq: 7.0e29 at 3^3)")
    gibbs_terms: Tagged = Stated(300, "app06 route (iv)", "'approximately 300 local jump operators'")
    gibbs_sweeps: Tagged = Stated(10, "app06 route (iv)", "'an optimistic ten sweeps'")
    gibbs_sampler_T: Tagged = Stated((1e10, 1e12), "app06 (pre-r20)", "RETIRED from the prose (E23 (3)): the "
                                                                      "unrepriced '1e10 ... toward 1e12'; record only")
    # E23 (3), r20: per-jump cost of a Chen-Kastoryano-Gilyen-type sampler (arxiv_2311_09207), ESTIMATED HERE.
    # One jump = operator Fourier transform of a local jump operator: controlled e^{-iHt} over a time register with
    # window T_OFT = w beta, energy resolution ~ 1/beta, plus the Metropolis weight on the frequency register.
    gibbs_n_evol: Tagged = Assumed((1, 4), "E23 (3), r20: controlled evolutions of length T_OFT per jump; low 1 (single "
                                           "evolution), high 4 (OFT e^{iHt} A e^{-iHt}, its uncompute, and the coherent term)")
    gibbs_n_evol_central: Tagged = Assumed(2, "E23 (3), r20: the OFT's forward and backward evolutions")
    gibbs_eps_oft: Tagged = Assumed(1e-3, "E23 (3), r20: Gaussian time window truncated at 1e-3, T_OFT = beta "
                                          "sqrt(ln 1/eps) (central and high; the low end takes T_OFT = beta)")
    gibbs_sweeps_high: Tagged = Assumed(30, "E23 (3), r20: the chapter calls ten sweeps optimistic; high end 30")
    gibbs_accept_T: Tagged = Assumed(1e3, "E23 (3), r20: Metropolis weight on the frequency register and the local jump "
                                          "LCU, an upper allowance per jump; negligible next to the evolution")
    sigma2_condensate: Tagged = Stated(10, "app06 shot-count paragraph", "'sigma^2 ~ 10 for the average over the 27 sites of "
                                                                         "the 2033 lattice'")
    sigma2_ncs: Tagged = Stated(1e2, "app06 shot-count paragraph", "'sigma^2_NCS ~ 1e2 at the longest t'. Chapter-open "
                                "pass (referee C3, 2026-10-05): an IDEAL-readout assumption; the Hadamard-test readout of "
                                "the Kubo form gives ncs_hadamard_rel_var (1.06e14 at 10a on the link-averaged free-field "
                                "background; >= 1.2e12 for any pairing)")
    # chapter-open pass (referee C1 / C3, 2026-10-05)
    ncs_e_max: Tagged = Assumed(1.5, "C3: ||E^a|| = j_max on the 2O irreps restricting SU(2) j <= 3/2 (the chapter's "
                                     "'Casimir cut at j = 3/2'); E^a = 0 on the other four irreps")
    ncs_b_max: Tagged = Assumed(2.0, "C3: B^a = -i Tr(sigma^a U_clover), |B^a| <= 2")
    tc_ew_gev: Tagged = Cited(159.0, "DOnofrio_Rummukainen_Tranberg_2014", "abstract: 'cross-over ... T_c = (159 +- 1) GeV'")
    mr_rate_box3: Tagged = Cited(0.0023, "Moore_Rummukainen_2000", "Table (Vol_table), a = 1/(2 g^2 T) (beta = 8): "
                                 "L g^2 T = 3.0 -> Gamma/alpha^4 T^4 = 0.0023 +- 0.0016")
    mr_rate_box10: Tagged = Cited(1.68, "Moore_Rummukainen_2000", "same table: L g^2 T = 10 -> 1.68 +- 0.03; 'large "
                                  "volume behavior is obtained by L = 8/g^2T'")
    mr_rate_box3_err: Tagged = Cited(0.0016, "Moore_Rummukainen_2000", "same table: the +-0.0016 on 0.0023 at L g^2 T = 3.0 "
                                     "(suppression 431-2400, quoted ~1e3)")
    eps_obs: Tagged = Assumed(0.1, "app06 'at eps=0.1'; box: '10% statistical'")
    shots_2033: Tagged = Stated(1e4, "app06 2033 box", "'1e4 (psi-bar psi: 1e3/coupling; Gamma_sph: 1e4/coupling)'")
    wall_per_shot_2033_stated: Tagged = Stated(1.9e4, "app06 requirements and wall-time paragraph", "'1.9e4 s' per 10a shot (18,597; C1 17,115) "
                                                                                                    "(r25: 16,241.7; r24 8.2e3; r22 2.008e5)")
    # r23 (H. Lamm 2026-10-02, apply_log/r23_ch09.md): the report assumes one machine. The Stated n_machines = 6.4e2
    # ("inside the 5-year horizon only on ~6.4e2 machines") is retired with n_machines_needed, wall_yr_parallel and
    # fits_campaign_horizon. Every wall time is serial on one machine: shots x (T x 1 us + ~0.1 ms).
    wall_yr_per_point_stated: Tagged = Stated(2.6, "RECORD (r24 box)",
                                              "'2.6 yr per point' for Gamma_sph (r24: 1e4 x 8241 s = 2.61 yr; r23 64); r25: 5.15 yr "
                                              "per temperature at 10a, inside the campaign")
    wall_yr_serial_stated: Tagged = Stated(1.3e2, "RECORD (r24 wall-time paragraph)", "'1.3e2 yr for the ~50-point scan' "
                                                                               "(r24: 130.6; r23 3.2e3); retired by r25 (e)")
    n_coupling: Tagged = Stated(50, "RECORD (r24 Eq. (Nshot), box, requirements)", "'N_coupling ~ 50'; r25 (e): 4-8 temperatures")
    campaign_yr: Tagged = Stated(5, "app06 requirements", "'Campaign horizon 5 years'")
    reduction_scan_stated: Tagged = Stated(26, "RECORD (r24 wall-time paragraph)", "'the scan needs its cost cut by 26' "
                                                                            "(r24: 130.6 / 5 = 26.1; r23 6.4e2)")
    t_shot_to_fit_scan_stated: Tagged = Stated(3.2e8, "RECORD (r24 wall-time paragraph)", "'a shot of 3.2e8 T' (3.156e8; "
                                                                                   "unchanged: set by the shot count)")
    # r24 (1): the condensate arm at its own shot (the preparation alone)
    t_shot_condensate_stated: Tagged = Stated(2.8e9, "app06 2033 box and wall-time paragraph",
                                              "'condensate shot 2.8e9 T' (TPQ alone at its own R-TOL, 17 filter calls: 2.7901e9; C1 1.3119e9; r25 4.4166e8, "
                                              "r24 2.2e8)")
    eps_l_condensate_stated: Tagged = Stated(3.6e-11, "app06 2033 box", "'eps_l <~ 3.6e-11' for the condensate (17 filter calls: 3.584e-11; C1: 7.623e-11; r25 2.264e-10; r24 4.5e-10)")
    wall_day_condensate_point_stated: Tagged = Stated(32, "app06 wall-time paragraph",
                                                      "'32 days per temperature' (1e3 condensate shots; 32.29; C1 15.18; r25 5.112; "
                                                      "r24 2.6 days per coupling point)")
    wall_yr_condensate_scan_stated: Tagged = Stated(0.36, "RECORD (r24 box)", "'0.36 yr for the ~50-point scan'; retired by r25 (e)")
    t_shot_kappa_r22_stated: Tagged = Stated(6.1e10, "app06 2033 derivation", "'at kappa ~ 1e2 the same circuit costs 5.3e10 T' "
                                                                              "(17 filter calls: 6.119e10; C1 5.631e10; r25: 5.344e10; r24 2.7e10)")
    # r24 levers for one sphaleron point under a year; sizes DERIVED HERE; none is applied to the box
    lever_eps_obs: Tagged = Assumed(0.15, "r24 lever: a 15% statistical target on Gamma_sph (shots sigma^2 / eps^2)")
    lever_low_step_stated: Tagged = Stated(1.4, "app06 wall-time paragraph", "'the low end of the step range (1.5)' (r24)")
    lever_eps15_stated: Tagged = Stated(2.25, "app06 wall-time paragraph", "'a 15% statistical target (2.25)' (r24)")
    wall_yr_lever_low_step_stated: Tagged = Stated(1.8, "app06 wall-time paragraph", "'1.8 yr' at the low end (r24)")
    wall_yr_lever_eps15_stated: Tagged = Stated(1.2, "app06 wall-time paragraph", "'1.2 yr' at 15% (r24)")
    wall_yr_lever_both_stated: Tagged = Stated(0.8, "app06 wall-time paragraph", "'together 0.8 yr' (r24)")
    # ---- r25 (H. Lamm 2026-10-02: rulings R1, R4, R9, R10 and the queued Ch. 9 rulings (a), (e), (f)) ----------------
    yukawa_placement: Tagged = Stated("inside", "author ruling r25 (a), H. Lamm 2026-10-02, relaying H. Singh",
                                      "'the Ginsparg-Wilson-projected Yukawa is INSIDE': the Yukawa couples the "
                                      "projected fields, whose projector (1 - sgn(H_W))/2 is a second QSVT polynomial of "
                                      "degree d_sgn on the same block encoding, once per application ('outside': r22-r24)")
    depth_variant_2033: Tagged = Assumed("fanout", "ruling R10 (H. Lamm 2026-10-02): where F* < 10 the box quotes the "
                                         "workspace fix; factory analysis: a CNOT fan-out of the shared control on "
                                         "5 extra LQ gives F* = 12.3 ('register' = as itemized, F* 5.6-6.4)")
    anc_fanout_2033: Tagged = Assumed(FANOUT_ANCILLA_2033, "r25 R10, factory analysis: ~5 LQ for the CNOT fan-out of "
                                      "the leaf flag (5 copy Toffolis), the orientation and Dirac flags and the mux ANDs")
    spatial_concurrency_2028: Tagged = Assumed(3, "r25 R10, DERIVED HERE: three spatial-hop phasing groups per slot on "
                                               "mutually non-adjacent links (disjoint sites, commuting strings); 3 x 10 "
                                               "adders fit inside the 31 of a k=32 group; restores F* >= 10 at 2028")
    eps_first: Tagged = Assumed(0.3, "ruling R1 (H. Lamm 2026-10-02): the minimal first result at 30% statistical error "
                                     "('for some of this physics even 30% statistical error would be interesting')")
    t_short_fm: Tagged = Assumed(1.0, "ruling r25 (e): the first result reads one temperature at 5a and 10a; 5a = 1 fm, "
                                      "priced at its own evolution depth (ruling R4)")
    ncs_var_ratio_short: Tagged = Assumed(4, "RECORD (shot audit, r25): under linear growth the 5a signal is half the 10a "
                                             "signal, so at the same absolute noise its relative variance is 4x. Replaced "
                                             "2026-10-05 (method-test ruling, no rate) by the derived (S10/S5)^2 = 4.95 of "
                                             "the free-field signals (method_test_var_ratio_5a)")
    n_temperatures: Tagged = Stated((4, 8), "author ruling r25 (e), H. Lamm 2026-10-02",
                                    "'campaign 4-8 temperatures (not 50 couplings)' on the symmetric / near-T_c side; "
                                    "the two-time linearity check at one of them, 10a only at the rest")
    t_shot_short_stated: Tagged = Stated(9.6e9, "app06 wall-time paragraph and box", "'160 steps of 23 applications ... "
                                                                                       "9.6e9 T' at 5a (17 filter calls: 9.5897e9; C1 8.109e9; r25: 7.2366e9)")
    eps_l_short_stated: Tagged = Stated(1.0e-11, "app06 2033 box", "'1.0e-11 at 5a' (17 filter calls: 1.043e-11; C1: 1.233e-11; r25: 1.382e-11)")
    wall_s_short_stated: Tagged = Stated(9.6e3, "app06 wall-time paragraph", "'or 9.6e3 s' per 5a shot (9589.7; C1 8109; r25 7236.6)")
    wall_s_condensate_stated: Tagged = Stated(2.8e3, "app06 wall-time paragraph", "'2.8e3 s' per condensate shot (17 filter calls: 2790.1; C1: 1311.9; r25 441.7)")
    f_star_register_stated: Tagged = Stated((5.6, 6.3), "app06 wall-time paragraph", "'only 5.6-6.3 T gates per reaction "
                                                                                  "layer' as itemized (10a shot: 5.59-6.32)")
    slowdown_register_stated: Tagged = Stated((1.6, 1.8), "app06 wall-time paragraph", "'1.6-1.8 times slower' (10 / F*)")
    f_star_2033_stated: Tagged = Stated(12, "app06 wall-time paragraph", "'raises this to 12' (12.17-12.54 over the runs)")
    shots_first_stated: Tagged = Stated((1.1e3, 5.5e3), "app06 box and wall-time paragraph", "'1.1e3 shots at 10a and 5.5e3 "
                                                                                           "at 5a' (1111, 5497; r25 4444 under linear growth)")
    wall_yr_first_stated: Tagged = Stated(2.3, "app06 box and wall-time paragraph", "first result (method test) '2.3 yr' (2.325; "
                                          "before the 2026-10-05 ruling 2.005 with the 5a ratio 4; C1 1.744; r25 1.591)")
    wall_yr_campaign_stated: Tagged = Stated((2.7, 3.0), "app06 box and wall-time paragraph", "campaign '2.7-3.0 yr' (method "
                                             "test + condensate at 4-8 temperatures: 2.679-3.032; the rate campaign before "
                                             "the 2026-10-05 ruling 36.08-60.00, record wall_yr_campaign_r25_stated)")
    campaign_over_horizon_stated: Tagged = Stated((0.54, 0.61), "app06 box and wall-time paragraph", "'inside the 5-year "
                                                  "horizon' (0.536-0.606; rate campaign before 2026-10-05: 7-12x)")
    wall_yr_campaign_r25_stated: Tagged = Stated((36, 60), "RECORD (app06 before the 2026-10-05 method-test ruling)",
                                                 "'36-60 yr for 4-8 temperatures, 7-12x the 5-year horizon' (rate campaign; "
                                                 "at 30%: '4.0-6.7 yr')")
    wall_s_first_2028_stated: Tagged = Stated(36, "app06 2028 box", "first result '36 s' (125 shots; 36.08)")
    wall_min_campaign_2028_stated: Tagged = Stated(13, "app06 2028 box", "campaign '13 min' (L5 = 2, 3, 4; 13.24)")
    lq_2028_stated: Tagged = Stated(250, "app06 2028 derivation and box", "'250 = 192 + 16 + 42' (r25 R10)")
    depth_step_2028_stated: Tagged = Stated((1.1e3, 1.3e3), "app06 2028 derivation", "'a T-depth of 1.1-1.3e3' per step "
                                                                                    "(1059.7-1318.7)")
    L5_campaign_2028: Tagged = Assumed((2, 4), "shot audit (r25): residual-mass scaling needs at least three L5 values; "
                                               "the 2028 campaign is L5 = 2, 3, 4 at the box's 1125 shots each")

    def __post_init__(self):
        for f in fields(self):
            v = getattr(self, f.name)
            if not isinstance(v, Tagged):
                raise TypeError(f"{f.name} must be Tagged")
            if v.is_range:
                if not (v.lo <= v.hi):
                    raise ValueError(f"{f.name}: range lo > hi")
                if v.lo <= 0:
                    raise ValueError(f"{f.name}: non-positive range")
            elif isinstance(v.value, (int, float)) and not isinstance(v.value, bool):
                if v.value <= 0:
                    raise ValueError(f"{f.name}: must be positive")
        if self.aa_call_rule.value not in ("2k+1", "per_round"):
            raise ValueError("aa_call_rule must be '2k+1' or 'per_round'")
        if self.spatial_hop_grouping.value not in _GROUPINGS:
            raise ValueError(f"spatial_hop_grouping must be one of {_GROUPINGS}")
        if self.toffoli_convention.value not in T_PER_TOFFOLI:
            raise ValueError(f"unknown Toffoli convention {self.toffoli_convention.value!r}")
        if self.n_umult_per_higgs_site.value < 1:
            raise ValueError("U_x per Higgs site: need >= 1")
        if self.onsite_term_2028.value not in ("wilson", "none"):
            raise ValueError("onsite_term_2028 must be 'wilson' or 'none'")
        if self.link_frame.value not in ("multiplexer", "draft"):
            raise ValueError("link_frame must be 'multiplexer' or 'draft'")
        if self.trotter_rule_2033.value not in ("state", "bound", "stated", "step"):
            raise ValueError("trotter_rule_2033 must be 'state', 'bound', 'stated' or 'step'")
        if self.sgn_architecture.value != "single_particle_lift":
            raise ValueError("sgn_architecture must be 'single_particle_lift': the second-quantized construction computes "
                             "the sign of the many-body operator and is kept as a record only")
        if not (0 < self.m0_wilson.value < 2):
            raise ValueError("m0_wilson must lie in the doubler window (0, 2)")
        if self.fifth_boundary.value not in ("open", "periodic"):
            raise ValueError("fifth_boundary must be 'open' or 'periodic'")
        if self.clifford_2O.value != "clifford":
            raise ValueError("clifford_2O is a fact, not a switch")
        if self.hop_share.value not in ("link", "draft"):
            raise ValueError("hop_share must be 'link' or 'draft'")
        if self.hop_fermion_2033.value != "wilson":
            raise ValueError("the 2033 overlap kernel is the Wilson operator")
        if self.yukawa_placement.value not in ("inside", "outside"):
            raise ValueError("yukawa_placement must be 'inside' (r25) or 'outside' (r22-r24 record)")
        if self.depth_variant_2033.value not in ("fanout", "register"):
            raise ValueError("depth_variant_2033 must be 'fanout' or 'register'")
        if not (0 < self.eps_first.value < 1):
            raise ValueError("eps_first must lie in (0, 1)")
        if self.hwp_chunk_k.value < 2:
            raise ValueError("hwp_chunk_k must be >= 2")
        for k in (self.gauge_group_2028.value, self.gauge_group_2033.value):
            if k not in GROUPS:
                raise ValueError(f"unknown gauge group {k!r}")
        if self.rot_synthesis.value not in ("rus", "watson", "rus-slope", "rs"):
            raise ValueError(f"unknown synthesis model {self.rot_synthesis.value!r}")
        if not (0 < self.eps_syn.value < 1 and 0 < self.eps_rot.value < 1 and 0 < self.delta_sgn.value < 1
                and 0 < self.delta_beta.value < 1 and 0 < self.haar_deficit_frac.value < 1
                and 0 < self.eps_qsp.value < 1 and 0 < self.faults_per_shot.value <= 1):
            raise ValueError("tolerances and fractions must lie in (0,1)")


# --------------------------------------------------------------------------- #
# 2028: Z3 domain-wall Trotter step under Hamming-weight phasing
# --------------------------------------------------------------------------- #

def _spatial_hop_rule(a: Assumptions, k: int, tc: str, eps: float, synthesis: str | None = None) -> dict:
    """Rule R-HOP (groups.rhop_link) for one Z3 spatial link shared by the Nc L5 color-and-s copies.

    The L5 fifth-direction slices play the role of the N_stag fields on one link (so the rule's
    default k = Nc L5 = 12 is Ch. 9's grouping ruling R6); the grouping passes its own k.  Z3 is
    diagonal on the register, so the rule prices no color move and 8 Pauli strings per copy.
    `eps` is the circuit's own R-TOL tolerance (_step_rtol_2028); the string count does not depend on it.
    r17: rotations at the full fit (a.hop_synthesis = "rus"); `synthesis` overrides it for the legacy record.
    """
    return rhop_link(a.gauge_group_2028.value, n_stag=int(a.L5_2028.value), n_c=int(a.Nc_2028.value),
                     eps=eps, k=k, synthesis=synthesis or a.hop_synthesis.value, toffoli_convention=tc)


def _hwp_groups_2028(a: Assumptions, grouping: str) -> list[tuple[str, int, int, CircuitStatus, str, bool]]:
    """(name, k, n_groups, status, note, rhop) for one Trotter step.

    Each same-angle group of k rotations costs floor(log2 k)+1 synthesized R_Z plus
    k-HW(k) Toffolis.  The Pauli strings of a hopping bilinear are grouped separately
    (app06:108).  rhop=True marks the spatial hops, priced by rule R-HOP (_step_cost_2028):
    8 strings per copy on a Z3 link (groups.rhop_link), 2 on the link-free fifth-direction hops.
    """
    D, L, L5 = int(a.D_2028.value), int(a.L_2028.value), int(a.L5_2028.value)
    V = L ** D
    Nc = int(a.Nc_2028.value)
    n_sp = V * L5                                     # spatial hops (V L5 = 32)
    n_5 = V * (L5 - 1 if a.fifth_boundary.value == "open" else L5)   # r21 (a): open fifth boundary, V (L5 - 1) = 24
    n_el = D * V                                      # electric terms (D V = 8)
    n_links = D * V                                   # s-independent links (8)
    ps_hop = int(a.n_pauli_per_hop.value)             # 2: the link-free fifth-direction hop
    ps_sp = _spatial_hop_rule(a, Nc * L5, a.toffoli_convention.value, float(a.eps_syn.value))["n_pauli"]   # 8 (R-HOP, Z3)
    ps_el = int(a.n_pauli_per_electric.value)
    groups = [
        ("fifth_hops", n_5 * Nc, ps_hop, CircuitStatus.SCALING,
         f"'The {n_5 * Nc} link-free fifth-direction hops ({n_5} hops in three colors) form one group per Pauli string' "
         f"(r21 a: fifth-direction boundary {a.fifth_boundary.value})", False),
    ]
    if grouping == "color+s":
        groups.append(("spatial_hops", Nc * n_sp // n_links, n_links * ps_sp, CircuitStatus.SCALING,
                       "'The twelve color and s copies of each spatial hop share one s-independent link': "
                       "k=Nc*L5 per (link, Pauli string); 8 strings per copy by rule R-HOP (derived here)", True))
    elif grouping == "color":
        groups.append(("spatial_hops", Nc, n_sp * ps_sp, CircuitStatus.SCALING,
                       "PRE-R6 sentence ('the three color copies of each spatial hop'): k=Nc per (hop, Pauli string); "
                       "8 strings per copy by rule R-HOP", True))
    else:
        groups.append(("spatial_hops", Nc * n_sp, ps_sp, CircuitStatus.SCALING,
                       "ALTERNATIVE (not the chapter's sentence; ignores that the dressed strings differ per link): "
                       "k=Nc*V*L5", True))
    groups.append(("electric", n_el, ps_el, CircuitStatus.SCALING,
                   "'three diagonal Z-strings per two-qubit link, each grouped across the eight links'", False))
    if a.onsite_term_2028.value == "wilson":
        # r17 (f): on-site m + r(1+d) of h_DW, 2N R_Z per site (FP section_mass.tex:12) on V L5 sites, one group
        k_os = wilson_mass_rotations(Nc, int(a.Nf_2028.value), D) * V * L5     # 6 x 32 = 192
        groups.append(("onsite_dw", k_os, 1, CircuitStatus.SCALING,
                       "on-site psi^dag gamma^0 (m + r(1+d)) psi of h_DW (arxiv_2505_20419), 2N Z per site "
                       f"({FP_KEY} section_mass.tex:12), diagonal gamma^0: one lattice-wide equal-|angle| group "
                       "(signs Clifford); DERIVED HERE r17", False))
    return groups


def _step_cost_2028(a: Assumptions, grouping: str, eps: float, chunk_k: int | None = None) -> dict:
    """One Trotter step at per-rotation tolerance `eps`.  `chunk_k` splits every phasing group larger
    than chunk_k (app06:108).

    The chapter's own rotations (fifth-direction hops, electric terms) are priced by its formula
    (a.rot_synthesis, the full RUS fit) at `eps`.  The spatial hop groups are priced by groups.rhop_link
    (rule R-HOP) at the rule's rotation price at the same `eps`; each rule group of n_pauli strings on
    one link is one call, so their T is n_units x rhop_link(k)["t"].  The counts do not depend on eps.
    """
    tc = a.toffoli_convention.value
    t_rot = t_per_rotation(eps, a.rot_synthesis.value)
    out = {"eps_rot": eps, "t_rot_chapter": t_rot, "groups": [], "n_rot": 0, "n_synth": 0, "n_synth_chapter": 0, "n_synth_rhop": 0, "n_toffoli": 0,
           "t_synth": 0.0, "t_synth_rhop": 0.0, "t_hwp": 0.0, "t_hwp_bound_k": 0.0, "max_ancilla": 0,
           "t_rhop": 0.0, "rhop": {}}
    for name, k0, n0, status, note, rule in _hwp_groups_2028(a, grouping):
        for k, n_groups in _chunk(k0, n0, chunk_k):
            if rule:
                rh = _spatial_hop_rule(a, k, tc, eps)
                n_units, rem = divmod(n_groups, rh["n_pauli"])
                assert rem == 0, (name, k, n_groups)
                t_rot_g = rh["hwp"]["t_rot"]
                n_synth = n_groups * rh["hwp"]["n_rot"]
                n_toff = n_groups * rh["hwp"]["n_toffoli"]
                out["t_rhop"] += n_units * rh["t"]
                out["rhop"][k] = rh
                out["n_synth_rhop"] += n_synth
                out["t_synth_rhop"] += n_synth * t_rot_g
            else:
                t_rot_g = t_rot
                n_synth = n_groups * hwp_synth_rotations(k)
                n_toff = n_groups * hwp_toffolis(k)
                out["n_synth_chapter"] += n_synth
            out["groups"].append(dict(name=name, k=k, n_groups=n_groups, n_rot=k * n_groups, n_synth=n_synth,
                                      n_toffoli=n_toff, t_rot=t_rot_g, rhop=rule, status=status, note=note))
            out["n_rot"] += k * n_groups
            out["n_synth"] += n_synth
            out["n_toffoli"] += n_toff
            out["t_synth"] += n_synth * t_rot_g
            out["t_hwp"] += toffoli_t(n_toff, tc)
            out["t_hwp_bound_k"] += T_PER_TOFFOLI[tc] * k * n_groups   # 'k - w(k)' is at most k Toffolis
            out["max_ancilla"] = max(out["max_ancilla"], hwp_ancilla(k))
    out["t_step"] = out["t_synth"] + out["t_hwp"]
    return out


def _step_rtol_2028(a: Assumptions, grouping: str, n_steps: int, chunk_k: int | None = None) -> dict:
    """One Trotter step priced under R-TOL: eps_rot = eps_rot_for(N_rot), N_rot = the synthesized
    rotations of ONE shot of this circuit (n_synth per step x n_steps).  The first call only counts:
    the counts are fixed by the grouping, not by eps, so there is no iteration (asserted).
    """
    eps_syn = float(a.eps_syn.value)
    n_synth = _step_cost_2028(a, grouping, eps_syn, chunk_k)["n_synth"]
    n_rot_shot = n_synth * n_steps
    out = _step_cost_2028(a, grouping, eps_rot_for(n_rot_shot, eps_syn), chunk_k)
    assert out["n_synth"] == n_synth
    out["n_rot_shot"] = n_rot_shot
    return out


def _model_2028(a: Assumptions) -> Result:
    g = GROUPS[a.gauge_group_2028.value]
    D, L, L5 = int(a.D_2028.value), int(a.L_2028.value), int(a.L5_2028.value)
    Nf, Nc = int(a.Nf_2028.value), int(a.Nc_2028.value)
    V = L ** D
    dirac = 2 ** ((D + 1) // 2)                       # app06:90: 2 in 1+1D
    lq_f = dirac * Nf * Nc * V * L5                   # 192
    lq_g = D * V * g.link_qubits                      # 16

    n_sp = V * L5
    n_5 = V * (L5 - 1 if a.fifth_boundary.value == "open" else L5)   # 24 (r21 a; 32 before)
    n_el = D * V
    n_terms = n_sp + n_5 + n_el                       # 64 hopping and electric terms
    n_steps = int(a.n_trotter_2028.value) + int(a.n_probe_steps.value)               # 21
    n_steps_before = int(a.n_trotter_2028_before_r6.value) + int(a.n_probe_steps.value)   # 31
    grouping = a.spatial_hop_grouping.value
    chunk_k = int(a.hwp_chunk_k.value)
    # R-TOL: one tolerance per circuit, from its own rotation count (5922 unchunked, 6384 chunked)
    step = _step_rtol_2028(a, grouping, n_steps)                     # unchunked default (13,145.2 T)
    chunked = _step_rtol_2028(a, grouping, n_steps, chunk_k)         # k <= 32 (14,160.6 T)
    eps_rot = step["eps_rot"]                                        # 1.2814e-3 (6090 rotations/shot)
    rh_sp = _spatial_hop_rule(a, Nc * L5, a.toffoli_convention.value, eps_rot)   # Z3 link (k=12), full fit: 1201.9 T
    hop_synthesis = a.hop_synthesis.value                              # 'rus' (r17 ruling (e))
    # the central r17 hop function prices the same Z3 link identically (R-HOP structure, full fit)
    hl_z3 = hop_link_cost(a.gauge_group_2028.value, n_stag=L5, eps=eps_rot, n_c=Nc, k=Nc * L5,
                          synthesis=hop_synthesis, toffoli_convention=a.toffoli_convention.value, legacy=False)
    assert math.isclose(hl_z3["t"], rh_sp["t"]), (hl_z3["t"], rh_sp["t"])
    n_hop_rot = n_5 * Nc * int(a.n_pauli_per_hop.value) + n_sp * Nc * rh_sp["n_pauli"]   # 192 + 768 = 960
    n_hop_rot_before_rhop = (n_sp + n_5) * Nc * int(a.n_pauli_per_hop.value)            # 384 (2 strings per hop)
    n_el_rot = n_el * int(a.n_pauli_per_electric.value)            # 24
    n_onsite_rot = sum(g["n_rot"] for g in step["groups"] if g["name"] == "onsite_dw")   # 192 (r17 f)
    # pre-r17 records: the same lattice without the on-site term, hop rotations at the full fit (2028-A) and at
    # R-HOP's retired slope-only price (the printed 1.95-2.05e5)
    # (records at the pre-r21 32 fifth-direction hops, so the old prints stay reproducible)
    a_periodic = replace(a, fifth_boundary=Cited("periodic", "record: the pre-r21 count V L5 = 32"))
    a_pre = replace(a_periodic, onsite_term_2028=Cited("none", "record: the step before r17 (f)"))
    pre = (_step_rtol_2028(a_pre, grouping, n_steps), _step_rtol_2028(a_pre, grouping, n_steps, chunk_k))
    pre_r21 = (_step_rtol_2028(a_periodic, grouping, n_steps), _step_rtol_2028(a_periodic, grouping, n_steps, chunk_k))
    t_shot_pre_r21 = tuple(n_steps * c["t_step"] for c in pre_r21)             # 276,050.9, 297,372.8 (r17-r20 print)
    # the on-site term's share of the current shot (open boundary, with and without the term)
    a_no_onsite = replace(a, onsite_term_2028=Cited("none", "sensitivity: the current step without the on-site term"))
    no_onsite = (_step_rtol_2028(a_no_onsite, grouping, n_steps), _step_rtol_2028(a_no_onsite, grouping, n_steps, chunk_k))

    def _slope_step(c):
        return c["t_step"] - c["t_synth_rhop"] + c["n_synth_rhop"] * t_per_rotation(
            c["eps_rot"], a.hop_synthesis_legacy.value)
    t_shot_no_onsite = tuple(n_steps * c["t_step"] for c in pre)               # 244,581.6, 254,029.9 (2028-A)
    t_shot_pre_r17 = tuple(n_steps * _slope_step(c) for c in pre)              # 195,122.4, 204,570.7 (old print)

    tc = a.toffoli_convention.value
    t_rot = step["t_rot_chapter"]                                    # 20.226 T (RUS fit at eps_rot)
    t_rot_retired = float(a.t_per_rot.value)                         # 30: the retired constant, records only
    eps_syn = float(a.eps_syn.value)
    others = {gname: _step_rtol_2028(a, gname, n_steps) for gname in _GROUPINGS if gname != grouping}
    t_shot = n_steps * step["t_step"]                 # 268,994.9
    t_shot_chunked = n_steps * chunked["t_step"]      # 288,535.5
    onsite_adds = (t_shot - n_steps * no_onsite[0]["t_step"], t_shot_chunked - n_steps * no_onsite[1]["t_step"])
    env = float(a.envelope_2028.value)

    # ancilla: which end fits the box's register (R5: no unitemized ancilla)
    anc_unchunked = step["max_ancilla"]               # 190 at k=192 (E26: k - w(k); 198 before r22)
    anc_chunked = chunked["max_ancilla"]              # 31 at k=32 (37 before r22)
    # E23 (5), r20; E26, r22: itemized peak. Schedule: the phasing groups run one at a time; during a group its adder
    # ancilla (one per adder; the weight bits are those ancilla and one data qubit) are live while its rotations are
    # synthesized one at a time (one RUS ancilla on top). The strings of a group are commuting and independent, mapped
    # in place to single-qubit Z's by a Clifford (no ancilla). The adiabatic ramp is the same step with ramped angles;
    # the initial state is a product state; the probe of psi-bar psi is read by measuring the fermion register. So the
    # peak is one phasing group plus RUS.
    anc_rus = int(a.anc_rus.value)                    # 1 (E26; was 2)
    lq_anc_serial = anc_chunked + anc_rus             # 32: one group at a time, rotations one at a time (r22-r24)
    # r25 (R10): that schedule has F* ~ 1.8 (all 336 rotations of a step serial on one RUS ancilla). Each slot now runs
    # its weight-bit rotations in parallel (one RUS ancilla each) and three spatial-hop groups at once on non-adjacent
    # links (hwp_step_depth_2028): peak 3 x (10 + 4) = 42 ancilla, F* >= 10.
    conc = int(a.spatial_concurrency_2028.value)      # 3
    dstep = hwp_step_depth_2028(chunked["groups"], conc, anc_rus)
    lq_anc = dstep["peak_ancilla"]                    # 42 (r22-r24 32, r20-r21 39, asserted 40)
    lq_anc_asserted = int(a.ancilla_2028_asserted.value)
    lq = lq_f + lq_g + lq_anc                         # 250 (r22-r24 240, r20-r21 247)
    lq_serial = lq_f + lq_g + lq_anc_serial           # 240 (record: the r22-r24 register)
    lq_single_copy = dirac * Nf * 1 * V * L5 + lq_g + lq_anc_serial   # '~110 LQ a single copy' (112)
    lq_unchunked = lq_f + lq_g + anc_unchunked + anc_rus      # 399 ('about 400'), serial schedule

    # sensitivity the chapter argues (app06:109): coherent synthesis.  Hold the total synthesis
    # error at R-TOL's eps_syn (incoherent: eps_syn = N eps_rot^2) and ask what per-rotation
    # tolerance a coherent (linear) sum would need for the same total: eps_syn / N.
    n_rot_shot = step["n_rot_shot"]                                  # 5922
    eps_tot = n_rot_shot * eps_rot ** 2                              # = eps_syn by construction
    eps_inc = eps_per_rotation(eps_tot, n_rot_shot, "incoherent")    # = eps_rot by construction
    eps_coh = eps_per_rotation(eps_syn, n_rot_shot, "coherent")      # 1.689e-6
    t_rot_coh_model = t_per_rotation(eps_coh, a.rot_synthesis.value)   # 31.252
    t_rot_coh_stated = float(a.t_per_rot_coherent_stated.value)
    # r17: the hop rotations are at the full fit too, so they move to the same 31.3 T
    t_rot_hop_coh = t_per_rotation(eps_coh, hop_synthesis)

    def _coh_step(c, t_rot_chapter, t_rot_hop=t_rot_hop_coh):
        return c["n_synth_chapter"] * t_rot_chapter + c["n_synth_rhop"] * t_rot_hop + c["t_hwp"]
    t_step_coh_stated = _coh_step(step, t_rot_coh_stated, t_rot_coh_stated)   # 290 x 31 + 1039 x 7 = 16263
    t_shot_coh_stated = n_steps * t_step_coh_stated
    t_step_coh_model = _coh_step(step, t_rot_coh_model)              # 290 x 31.298 + 1039 x 7 = 16349.5
    t_shot_coh_model = n_steps * t_step_coh_model                    # 343,339
    # chunked: its own count (6384 per shot) sets its own coherent tolerance
    eps_coh_ch = eps_per_rotation(eps_syn, chunked["n_rot_shot"], "coherent")
    t_shot_coh_chunked = n_steps * _coh_step(chunked, t_per_rotation(eps_coh_ch, a.rot_synthesis.value),
                                             t_per_rotation(eps_coh_ch, hop_synthesis))

    # r22 (E26 (i)), DERIVED HERE: the 1+1D overlap route the chapter rejects, in the valid single-particle architecture
    # (arxiv_2607_28524), priced as the 2033 query is. Index register: color (2 qubits), spin (1), coordinate (3); five
    # LCU terms. A query reads one of the eight two-qubit links and applies the phase omega^(+-g). One application of the
    # block-encoded overlap Hamiltonian is d_sgn + 1 queries plus the lift. Each circuit (one application, or a
    # ground-state preparation of n_sgn of them) sets its own R-TOL tolerance.
    q_modes_1p1 = dirac * Nf * Nc * V                                  # 48
    spq1 = single_particle_query_Z3_1p1d(L, g.link_qubits)             # 40 Toffolis + 6 rotations
    lift1 = lift_cost(q_modes_1p1, 1)                                  # 190 Toffolis; the color index is uniform over 3

    def _overlap_1p1d(kappa, n_sgn):
        d = math.ceil(kappa * math.log(1 / float(a.delta_1p1d_a.value)))
        nq = d + 1
        e = eps_rot_for(n_sgn * (nq * spq1["n_rot"] + lift1["n_rot"]), eps_syn)
        tr = t_per_rotation(e, a.rot_synthesis.value)
        c_sp1 = toffoli_t(spq1["toffoli"], tc) + spq1["n_rot"] * tr
        c_lift1 = toffoli_t(lift1["toffoli"], tc) + lift1["t_direct"] + lift1["n_rot"] * tr
        return {"d_sgn": d, "eps_rot": e, "c_sp": c_sp1, "t_app": nq * c_sp1 + c_lift1,
                "t": n_sgn * (nq * c_sp1 + c_lift1)}
    n_gs = int(a.n_sgn_gs_prep.value)
    ov_a1, ov_a10 = _overlap_1p1d(float(a.kappa_1p1d_a.value), 1), _overlap_1p1d(float(a.kappa_1p1d_a.value), n_gs)
    ov_b1, ov_b10 = _overlap_1p1d(float(a.kappa_1p1d_b.value), 1), _overlap_1p1d(float(a.kappa_1p1d_b.value), n_gs)
    # RECORD (r21 3e): the retired second-quantized pricing of the same kernel: an LCU of 30 Pauli strings per site at one
    # SELECT Toffoli each plus 25 rotations per query, d_sgn queries per sgn. It prices the sign of the many-body operator.
    n_spin_1p1 = dirac // 2                                            # 1
    lcu1 = lcu_terms_per_site(D, n_spin_1p1, Nc, Nf, rh_sp["n_pauli"])   # 24 hop + 6 on-site = 30 per site
    assert rh_sp["n_pauli"] == RHOP_N_PAULI_DIAGONAL[a.gauge_group_2028.value]
    lcu1_terms = V * lcu1["total"]                                     # 240
    c_toff_1p1 = toffoli_t(lcu1_terms * float(a.lcu_toffoli_per_term.value), tc)   # 1680
    n_rot_be = float(a.n_rot_be.value)                                 # 25

    def _overlap_1p1d_sq(kappa, n_sgn):
        d = math.ceil(kappa * math.log(1 / float(a.delta_1p1d_a.value)))
        e = eps_rot_for(n_rot_be * d * n_sgn, eps_syn)
        c_be = c_toff_1p1 + n_rot_be * t_per_rotation(e, a.rot_synthesis.value)
        return {"d_sgn": d, "c_be": c_be, "t": n_sgn * d * c_be}
    ov_sq = {"lcu_per_site": lcu1, "lcu_terms": lcu1_terms, "c_toffoli_T": c_toff_1p1,
             "c_be": _overlap_1p1d_sq(float(a.kappa_1p1d_a.value), 1)["c_be"],
             "c_be_stated": float(a.c_be_1p1d_sq_stated.value),
             "sgn_T_kappa30": _overlap_1p1d_sq(float(a.kappa_1p1d_a.value), 1)["t"],
             "gs_prep_T": _overlap_1p1d_sq(float(a.kappa_1p1d_a.value), n_gs)["t"],
             "sgn_T_kappa50": _overlap_1p1d_sq(float(a.kappa_1p1d_b.value), 1)["t"]}

    # record: the pre-r17 counts at the retired fixed tolerance (30 T per rotation, R-HOP hop at 1e-4)
    t_rot_hop_retired = t_per_rotation(float(a.eps_rot.value), a.hop_synthesis_legacy.value)   # 15.281

    def _retired_step(c):
        return c["n_synth_chapter"] * t_rot_retired + c["n_synth_rhop"] * t_rot_hop_retired + c["t_hwp"]
    t_shot_before_rtol = (n_steps * _retired_step(pre[0]), n_steps * _retired_step(pre[1]))   # 223333, 236899

    # r21 (3d), DERIVED HERE: shots from the estimator variance. The chapter's sigma^2 ~ 10 is the condensate averaged
    # over the 27 sites of the 2033 lattice; here N_c V = 24 site-color copies average (the color copies are decoupled,
    # hence independent). If the sites of one color copy were fully correlated only the N_c copies would average.
    eps_obs = float(a.eps_obs.value)
    copies_2033 = int(a.L_2033.value) ** int(a.D_2033.value)          # 27
    copies_2028 = Nc * V                                              # 24
    sigma2_2028 = float(a.sigma2_condensate.value) * copies_2033 / copies_2028   # 11.25
    shots = sigma2_2028 / eps_obs ** 2                                # 1125
    sigma2_worst = float(a.sigma2_condensate.value) * copies_2033 / Nc    # 90
    shots_worst = sigma2_worst / eps_obs ** 2                         # 9000
    # operator-norm bound: one copy's bilinear a^dag b + b^dag a has eigenvalues {-1, 0, 0, 1}, so its variance is at
    # most 1; the chapter's sigma^2 (variance over squared mean) per copy, 10 x 27 = 270, then needs |<o>| >= 0.061
    var_copy_max = 1.0
    condensate_min = math.sqrt(var_copy_max / (float(a.sigma2_condensate.value) * copies_2033))   # 0.0609
    t_gate = float(a.t_gate_s.value)
    t0 = float(a.shot_overhead_s.value)                               # 0.1 ms per shot
    wall_shot = (t_shot * t_gate + t0, t_shot_chunked * t_gate + t0)
    faults = float(a.faults_per_shot.value)
    eps_l = (faults / t_shot_chunked, faults / t_shot)

    # ---- r25: tiers (R1), T-depth and factories (R9, R10). The first result is the condensate at L5 = 4 at 30%; the
    # campaign is residual-mass scaling, L5 = 2, 3, 4 at the box's 10% (1125 shots each), each L5 priced on its own
    # register at its own R-TOL tolerance. Depth: the chunked schedule above (the one that fits the register).
    eps_first = float(a.eps_first.value)
    shots_first = sigma2_2028 / eps_first ** 2                        # 125
    depth_shot = (n_steps * dstep["lo"], n_steps * dstep["hi"])
    inter_fstar_serial = chunked["t_step"] / (chunked["n_synth"] * chunked["t_rot_chapter"]
                                              + sum(g["n_groups"] * math.ceil(2 * math.log2(g["k"]))
                                                    for g in chunked["groups"]))     # 1.84: r22-r24 schedule
    d_serial = hwp_step_depth_2028(chunked["groups"], 1, anc_rus)     # record: rotations parallel, one group per slot
    l5_runs = {}
    for L5_ in range(int(a.L5_campaign_2028.lo), int(a.L5_campaign_2028.hi) + 1):
        a_ = replace(a, L5_2028=Assumed(L5_, "r25 campaign point"))
        ch_ = _step_rtol_2028(a_, grouping, n_steps, chunk_k)
        d_ = hwp_step_depth_2028(ch_["groups"], conc, anc_rus)
        t_ = n_steps * ch_["t_step"]
        # R9: at the 10-factory baseline a shot takes max(N_T t_gate, D_T t_r) (the corrected wall where F* < 10)
        l5_runs[L5_] = {"t_shot": t_, "t_depth": (n_steps * d_["lo"], n_steps * d_["hi"]),
                        "lq": dirac * Nf * Nc * V * L5_ + lq_g + d_["peak_ancilla"], "shots": shots,
                        "wall_s": shots * (max(t_ * t_gate, n_steps * d_["hi"] * REACTION_TIME_S) + t0),
                        "wall_s_baseline": shots * (t_ * t_gate + t0)}
    if L5 in l5_runs:
        assert math.isclose(l5_runs[L5]["t_shot"], t_shot_chunked)
    wall_first_2028 = shots_first * (t_shot_chunked * t_gate + t0)
    wall_campaign_2028 = sum(v["wall_s"] for v in l5_runs.values())
    ex_l5 = {k: depth_exports(v["t_shot"], v["t_depth"], v["shots"], a.t_gate_s, a.shot_overhead_s)
             for k, v in l5_runs.items()}
    ex28 = depth_exports(t_shot_chunked, depth_shot, shots, a.t_gate_s, a.shot_overhead_s)

    breakdown = []
    for grp in step["groups"]:
        breakdown.append(Primitive(
            f"rot_synth_{grp['name']}", grp["n_synth"] * n_steps, grp["t_rot"], grp["status"],
            RHOP_SRC if grp["rhop"] else SYNTHESIS_SRC,
            (f"Z3 hop, R-HOP structure (groups.rhop_link = hop_link_cost('Z3'), DERIVED HERE): full fit at "
             f"eps_rot={eps_rot:.4g} (R-TOL: sqrt({eps_syn:g}/{n_rot_shot})) = {grp['t_rot']:.3f} T; " if grp["rhop"] else
             f"eps_rot={eps_rot:.4g} (R-TOL: sqrt({eps_syn:g}/{n_rot_shot})), RUS fit {grp['t_rot']:.3f} T; ")
            + f"{grp['n_groups']} group(s) of k={grp['k']} -> {grp['n_synth']}/step; {grp['note']}"))
        breakdown.append(Primitive(
            f"hwp_toffoli_{grp['name']}", grp["n_toffoli"] * n_steps, T_PER_TOFFOLI[tc], CircuitStatus.COMPILED,
            "arxiv_1902_10673,arxiv_1709_06648",
            f"Hamming-weight adders, k-HW(k) Toffolis per group at {T_PER_TOFFOLI[tc]} T ({tc}, ruling R5); "
            f"the count is the sources', whose adders uncompute by measurement"))

    def _summary(c):
        return dict(groups=tuple((x["name"], x["k"], x["n_groups"]) for x in c["groups"]),
                    n_synth_per_step=c["n_synth"], n_hwp_toffoli_per_step=c["n_toffoli"],
                    t_synth_per_step=c["t_synth"], t_hwp_per_step=c["t_hwp"], t_per_step=c["t_step"],
                    t_per_shot=n_steps * c["t_step"], synth_fraction_of_step=c["t_synth"] / c["t_step"],
                    max_ancilla=c["max_ancilla"], n_rot_shot=c["n_rot_shot"], eps_rot=c["eps_rot"])

    inter = {
        "lq_fermion": lq_f, "lq_gauge": lq_g, "lq_ancilla": lq_anc, "lq_total": lq,
        "lq_single_color_copy": lq_single_copy, "lq_single_color_copy_stated": float(a.lq_single_copy_stated.value),
        "link_qubits": g.link_qubits,
        "n_spatial_hops": n_sp, "n_fifth_hops": n_5, "n_electric_terms": n_el,
        "n_terms_per_step": n_terms,
        "n_hop_rotations_per_step": n_hop_rot, "n_electric_rotations_per_step": n_el_rot,
        "n_hop_rotations_per_step_before_rhop": n_hop_rot_before_rhop,
        # rule R-HOP on the spatial Z3 links (groups.rhop_link; derived here, ruled report-wide)
        "hop_synthesis": hop_synthesis,
        "rhop_n_pauli_spatial": rh_sp["n_pauli"], "rhop_k_spatial": rh_sp["k"],
        "rhop_t_move_spatial": rh_sp["t_moves"],
        "rhop_t_rot": rh_sp["hwp"]["t_rot"], "rhop_hwp_t": rh_sp["hwp"]["t"],
        "rhop_hwp_ancilla": rh_sp["hwp"]["ancilla"],
        "rhop_t_per_link": rh_sp["t"],
        "rhop_n_links": D * V,
        "t_hop_rhop_per_step": step["t_rhop"],
        "chunked_t_hop_rhop_per_step": chunked["t_rhop"],
        "n_synth_chapter_per_step": step["n_synth_chapter"], "n_synth_rhop_per_step": step["n_synth_rhop"],
        "t_synth_rhop_per_step": step["t_synth_rhop"],
        "hop_link_cost_z3_t": hl_z3["t"],                                # = rhop_t_per_link (central r17 function)
        "n_onsite_rotations_per_step": n_onsite_rot,                     # 192 (r17 f)
        "onsite_term": a.onsite_term_2028.value,
        "t_per_shot_no_onsite": t_shot_no_onsite,                        # 244,581.6-254,029.9 (2028-A alone, 32 hops)
        "t_per_shot_pre_r17": t_shot_pre_r17,                            # 195,122.4-204,570.7 (the old print)
        "t_per_shot_pre_r21": t_shot_pre_r21,                            # 276,050.9-297,372.8 (32 fifth-direction hops)
        "fifth_boundary": a.fifth_boundary.value,
        "onsite_term_adds_per_shot": onsite_adds,                        # the on-site term's share of the current shot
        "hop_synthesis_legacy": a.hop_synthesis_legacy.value,
        # R-TOL: one tolerance per circuit from its own rotation count (common.eps_rot_for)
        "eps_syn": eps_syn,
        "n_rot_shot": n_rot_shot, "chunked_n_rot_shot": chunked["n_rot_shot"],             # 5922, 6384
        "eps_rot": eps_rot, "chunked_eps_rot": chunked["eps_rot"],                          # 1.2995e-3, 1.2516e-3
        "t_per_rot_chapter": t_rot, "chunked_t_per_rot_chapter": chunked["t_rot_chapter"],  # 20.226, 20.288
        "t_per_rot_rhop": step["rhop"][Nc * L5]["hwp"]["t_rot"],                           # 20.249 (r17: full fit)
        "chunked_t_per_rot_rhop": chunked["rhop"][Nc * L5]["hwp"]["t_rot"],                # 20.381
        "t_per_rot_retired": t_rot_retired, "t_per_rot_hop_retired": t_rot_hop_retired,
        "t_per_shot_before_rtol": t_shot_before_rtol,
        "rotations_per_step": step["n_rot"],
        "t_per_toffoli": T_PER_TOFFOLI[tc],
        "grouping": grouping,
        "groups": tuple((x["name"], x["k"], x["n_groups"]) for x in step["groups"]),
        "n_synth_per_step": step["n_synth"],
        "n_hwp_toffoli_per_step": step["n_toffoli"],
        "t_synth_per_step": step["t_synth"],
        "t_hwp_per_step": step["t_hwp"],
        "t_hwp_bound_k_per_step": step["t_hwp_bound_k"],
        "t_per_step": step["t_step"],
        "t_per_step_stated": (float(a.t_per_step_2028_stated.lo), float(a.t_per_step_2028_stated.hi)),
        "synth_fraction_of_step": step["t_synth"] / step["t_step"],
        "t_per_rot_rus_fit": t_per_rotation(eps_rot),
        "n_trotter": int(a.n_trotter_2028.value),
        "n_steps_total": n_steps,
        "t_per_shot_unchunked": t_shot,
        # r25 (R9): t_per_shot is the shot the 250-LQ register runs, the chunked one (the unchunked end needs ~400 LQ)
        "t_per_shot": t_shot_chunked,
        "t_per_shot_stated": (float(a.t_per_shot_2028_stated.lo), float(a.t_per_shot_2028_stated.hi)),
        "fits_2028_envelope": t_shot <= env,
        "ratio_to_envelope": (t_shot / env, t_shot_chunked / env),
        "ratio_to_envelope_stated": (float(a.ratio_to_envelope_2028_stated.lo),
                                     float(a.ratio_to_envelope_2028_stated.hi)),
        # chunked to k <= hwp_chunk_k so that phasing fits the box's ancilla (app06:108, 148)
        "hwp_chunk_k": chunk_k,
        "chunked_groups": tuple((x["name"], x["k"], x["n_groups"]) for x in chunked["groups"]),
        "chunked_n_synth_per_step": chunked["n_synth"],
        "chunked_n_hwp_toffoli_per_step": chunked["n_toffoli"],
        "chunked_t_synth_per_step": chunked["t_synth"],
        "chunked_t_hwp_per_step": chunked["t_hwp"],
        "chunked_t_per_step": chunked["t_step"],
        "chunked_t_per_shot": t_shot_chunked,
        "chunked_synth_fraction_of_step": chunked["t_synth"] / chunked["t_step"],
        "chunking_delta_t_per_step": chunked["t_step"] - step["t_step"],
        "chunking_delta_n_synth": chunked["n_synth"] - step["n_synth"],
        "chunking_delta_n_toffoli": chunked["n_toffoli"] - step["n_toffoli"],
        # ancilla itemization (R5)
        "ancilla_itemized": (("hwp_chunked_k32", anc_chunked), ("rus_synthesis", anc_rus)),
        "lq_ancilla_asserted": lq_anc_asserted, "anc_rus": anc_rus,
        "hwp_ancilla_unchunked": anc_unchunked,
        "hwp_ancilla_unchunked_stated": float(a.hwp_ancilla_unchunked_stated.value),
        "hwp_ancilla_chunked": anc_chunked,
        "hwp_ancilla_chunked_stated": float(a.hwp_ancilla_chunked_stated.value),
        "unchunked_fits_box_ancilla": anc_unchunked + anc_rus <= lq_anc,
        "chunked_fits_box_ancilla": anc_chunked + anc_rus <= lq_anc,
        "lq_unchunked": lq_unchunked,
        "lq_unchunked_stated": float(a.lq_unchunked_stated.value),
        # the ramp before R6, so the trade (T-count against diabatic error) is visible
        "n_steps_total_before_r6": n_steps_before,
        "t_per_shot_at_30_step_ramp": (n_steps_before * step["t_step"], n_steps_before * chunked["t_step"]),
        # the other groupings, for the record
        "other_groupings": {k: _summary(v) for k, v in others.items()},
        # coherent-synthesis sensitivity (app06:109)
        "n_synth_per_shot": n_rot_shot,
        "eps_total_synthesis_incoherent": eps_tot,
        "eps_rot_incoherent": eps_inc,
        "eps_rot_coherent_same_total": eps_coh,
        "t_per_rot_coherent_model": t_rot_coh_model,
        "t_per_rot_rs_at_eps_rot": t_per_rotation(eps_rot, "rs"),          # 28.8: Ross-Selinger at the same eps
        "t_per_rot_coherent_stated": t_rot_coh_stated,
        "t_per_rot_hop_coherent": t_rot_hop_coh,
        "t_per_step_coherent_at_stated": t_step_coh_stated,
        "t_per_step_coherent_stated": float(a.t_per_step_coherent_stated.value),
        "t_per_step_coherent_model": t_step_coh_model,
        "t_per_shot_coherent_at_stated": t_shot_coh_stated,
        "t_per_shot_coherent_stated": float(a.t_per_shot_coherent_stated.value),
        "t_per_shot_coherent_model": t_shot_coh_model,
        "t_per_shot_coherent_chunked_at_stated": t_shot_coh_chunked,
        # 1+1D overlap (rejected route), r22: the single-particle architecture; the r21 (3e) pricing is a record
        "overlap_1p1d_q_modes": q_modes_1p1, "overlap_1p1d_index_register": spq1["index_register"],
        "overlap_1p1d_query": spq1, "overlap_1p1d_query_toffoli": spq1["toffoli"], "overlap_1p1d_query_n_rot": spq1["n_rot"],
        "overlap_1p1d_lift": lift1,
        "overlap_1p1d_d_sgn_kappa30": ov_a1["d_sgn"], "overlap_1p1d_d_sgn_kappa50": ov_b1["d_sgn"],
        "overlap_1p1d_c_sp": ov_a1["c_sp"], "overlap_1p1d_c_sp_stated": float(a.c_sp_1p1d_stated.value),
        "overlap_1p1d_c_sp_all": (ov_a1["c_sp"], ov_a10["c_sp"], ov_b1["c_sp"], ov_b10["c_sp"]),
        "overlap_1p1d_sgn_T_kappa30": ov_a1["t"], "overlap_1p1d_sgn_T_kappa30_stated": float(a.sgn_T_1p1d_a.value),
        "overlap_1p1d_sgn_T_kappa50": ov_b1["t"], "overlap_1p1d_sgn_T_kappa50_stated": float(a.sgn_T_1p1d_b.value),
        "overlap_1p1d_gs_prep_T": ov_a10["t"], "overlap_1p1d_gs_prep_T_stated": float(a.gs_prep_1p1d_stated.value),
        "overlap_1p1d_gs_prep_T_kappa50": ov_b10["t"],
        "overlap_1p1d_gs_prep_T_kappa50_stated": float(a.gs_prep_1p1d_b_stated.value),
        "overlap_1p1d_gs_prep_over_dwf": ov_a10["t"] / t_shot_chunked,
        "overlap_1p1d_gs_prep_over_dwf_stated": float(a.gs_prep_over_dwf_stated.value),
        "overlap_1p1d_gs_prep_over_envelope": ov_a10["t"] / env,
        "overlap_1p1d_gs_prep_over_envelope_stated": float(a.gs_prep_over_envelope_stated.value),
        "overlap_1p1d_sgn_over_envelope": ov_a1["t"] / env,
        "overlap_1p1d_sgn_over_dwf": ov_a1["t"] / t_shot_chunked,
        "overlap_1p1d_orders_past_envelope": math.log10(ov_b10["t"] / env),
        "overlap_1p1d_d_sgn_cut_to_fit": (ov_a10["t"] / env, ov_b10["t"] / env),
        "overlap_1p1d_second_quantized_record": ov_sq,
        # shots / wall time / implied logical error (R3: 0.1 faults per shot)
        "sigma2_condensate_2028": sigma2_2028, "condensate_copies": (copies_2033, copies_2028),
        "sigma2_condensate_2028_worst": sigma2_worst, "shots_worst": shots_worst,
        "condensate_variance_per_copy_max": var_copy_max, "condensate_per_copy_min_for_sigma2": condensate_min,
        "shots_worst_stated": float(a.shots_2028_worst_stated.value),
        "shots": shots, "shots_stated": float(a.shots_2028_stated.value),
        "shot_overhead_s": t0,
        "wall_s_per_shot": wall_shot,
        "wall_s_per_shot_stated": (float(a.wall_per_shot_2028_stated.lo), float(a.wall_per_shot_2028_stated.hi)),
        "wall_s_total": (wall_shot[0] * shots, wall_shot[1] * shots),
        "faults_per_shot": faults,
        "eps_l_implied": eps_l,
        "eps_l_stated": float(a.eps_l_2028_stated.value),
        # r25: T-depth and factories (R9, R10), the two tiers (R1)
        "lq_ancilla_serial_schedule": lq_anc_serial, "lq_serial_schedule": lq_serial,
        "spatial_concurrency": conc, "depth_per_step": (dstep["lo"], dstep["hi"]), "slots_per_step": dstep["n_slots"],
        "depth_per_step_one_group_per_slot": (d_serial["lo"], d_serial["hi"]),
        "f_star_one_group_per_slot": (chunked["t_step"] / d_serial["hi"], chunked["t_step"] / d_serial["lo"]),
        "f_star_serial_rotations": inter_fstar_serial,
        "shots_first": shots_first, "l5_campaign": l5_runs, "exports_by_l5": ex_l5,
        "t_per_shot_exports": ex28["t_per_shot"],
        "t_depth_per_shot": ex28["t_depth_per_shot"], "f_star": ex28["f_star"],
        "floor_wall_s": (sum(e["floor_wall_s"][0] for e in ex_l5.values()), sum(e["floor_wall_s"][1] for e in ex_l5.values())),
        "factories_for_1yr": (sum(e["factories_for_1yr"][0] for e in ex_l5.values()),
                              sum(e["factories_for_1yr"][1] for e in ex_l5.values())),
        "wall_first_result_s": (wall_first_2028, wall_first_2028),
        "wall_campaign_s": (wall_campaign_2028, wall_campaign_2028),
        "wall_s_first_stated": float(a.wall_s_first_2028_stated.value),
        "wall_min_campaign_stated": float(a.wall_min_campaign_2028_stated.value),
        "lq_stated": float(a.lq_2028_stated.value),
        "depth_per_step_stated": (float(a.depth_step_2028_stated.lo), float(a.depth_step_2028_stated.hi)),
    }
    notes = (
        f"2028 grouping = {grouping!r} (app06:108, ruling R6) plus the on-site term (r17 f): {step['n_synth']} "
        f"synthesized rotations + {step['n_toffoli']} HWP Toffolis at {T_PER_TOFFOLI[tc]} T (R5) per step -> "
        f"{step['t_step']:.0f} T/step; x {n_steps} steps (N_Trotter = {int(a.n_trotter_2028.value)} + probe) = "
        f"{t_shot:.0f} T/shot, {t_shot / env:.2f}x the 1e5 envelope.",
        f"Spatial hops (r17 ruling (e)): Z3 has no Fermion_Primitives entry, so the R-HOP structure stays (C_W = "
        f"{rh_sp['t_moves']:.0f}, {rh_sp['n_pauli']} Pauli strings per copy, one HWP group of k={rh_sp['k']} each), now "
        f"priced at the full fit: {rh_sp['hwp']['n_rot']} rotations at {rh_sp['hwp']['t_rot']:.3f} T + "
        f"{rh_sp['hwp']['n_toffoli']} Toffolis = {rh_sp['hwp']['t']:.2f} T per string, {rh_sp['t']:.2f} T per link "
        f"(= groups.hop_link_cost('Z3')), {step['t_rhop']:.1f} T per step on {D * V} links. At R-HOP's retired "
        f"slope-only price the pre-r17 step gave {t_shot_pre_r17[0]:.0f}-{t_shot_pre_r17[1]:.0f} T/shot.",
        f"On-site term (r17 f): m + r(1+d) of h_DW (arXiv:2505.20419), 2N R_Z per site (FP section_mass.tex:12) on "
        f"{V * L5} sites = {n_onsite_rot} equal-|angle| Z rotations, one HWP group (8 synthesized + 190 Toffolis; "
        f"chunked 6 x k=32). Without it the shot is {t_shot_no_onsite[0]:.0f}-{t_shot_no_onsite[1]:.0f} T.",
        f"R-TOL: total synthesis error per shot {eps_syn:g}, tolerance from the circuit's own count. "
        f"Unchunked: {n_rot_shot} rotations/shot ({step['n_synth']} x {n_steps}) -> eps_rot = {eps_rot:.4g}, "
        f"{t_rot:.3f} T. Chunked: {chunked['n_rot_shot']} -> {chunked['eps_rot']:.4g}, {chunked['t_rot_chapter']:.3f} T. "
        f"Before R-TOL (pre-r17 counts at 1e-4, 30 / 15.281 T): {t_shot_before_rtol[0]:.0f}-{t_shot_before_rtol[1]:.0f}.",
        f"Chunked to k <= {chunk_k}: {chunked['n_synth']} synthesized rotations + {chunked['n_toffoli']} Toffolis -> "
        f"{chunked['t_step']:.0f} T/step, {t_shot_chunked:.0f} T/shot, {t_shot_chunked / env:.2f}x the envelope. "
        f"Chunking adds {chunked['n_synth'] - step['n_synth']} synthesized rotations and "
        f"{chunked['n_toffoli'] - step['n_toffoli']:+d} Toffolis ({chunked['t_step'] - step['t_step']:+.0f} T/step).",
        f"Ancilla (E23 (5), E26; r25 R10): {dstep['n_slots']} slots per step, three spatial-hop groups per slot on "
        f"non-adjacent links, each weight bit on its own RUS ancilla: peak {lq_anc} (3 x (10 + 4)); register {lq} LQ. "
        f"The r22-r24 schedule (one group, one rotation at a time) held {lq_anc_serial} ({lq_serial} LQ) at F* "
        f"{inter_fstar_serial:.2f}. The unchunked k={max(g['k'] for g in step['groups'])} group holds {anc_unchunked}: "
        f"register {lq_unchunked} LQ. hard_ops = (unchunked, chunked); only the chunked end fits the {lq}-LQ register.",
        f"T-depth (r25, R9/R10): {dstep['lo']:.1f}-{dstep['hi']:.1f} per step, {depth_shot[0]:.0f}-{depth_shot[1]:.0f} per "
        f"shot; F* = {ex28['f_star'][0]:.2f}-{ex28['f_star'][1]:.2f} (one group per slot {chunked['t_step'] / d_serial['hi']:.2f}-"
        f"{chunked['t_step'] / d_serial['lo']:.2f}). First result: L5 = {L5} at {eps_first:.0%}, {shots_first:.0f} shots, "
        f"{wall_first_2028:.1f} s. Campaign: L5 = " + ", ".join(f"{k} ({v['lq']} LQ, {v['t_shot']:.4g} T)" for k, v in l5_runs.items())
        + f" at {shots:.0f} shots each, {wall_campaign_2028:.1f} s.",
        f"At the pre-R6 30-step ramp the same step costs give {n_steps_before * step['t_step']:.3g}-"
        f"{n_steps_before * chunked['t_step']:.3g} T/shot. The diabatic error of the 20-step ramp is not estimated "
        "anywhere; the chapter flags it as the quantity the benchmark must measure.",
        "The k-HW(k) Toffoli count is that of arXiv:1709.06648 / 1902.10673, whose adders uncompute by measurement; "
        "it is priced at 7 T per Toffoli as ruled (R5). A unitary uncompute would double the Toffoli count.",
        f"Fifth-direction hops (r21 a, derived here): open boundary, V (L5 - 1) = {n_5} hops, two groups of k = {n_5 * Nc}. "
        f"With the pre-r21 V L5 = 32 the shot was {t_shot_pre_r21[0]:.0f}-{t_shot_pre_r21[1]:.0f} T.",
        "On-site term basis: kept the draft's diagonal-gamma^0 encoding (E23 (4)). In an off-diagonal gamma^0 basis it "
        "would be two k=96 groups (14 rotations + 188 Toffolis); the chunked end is the same in either basis.",
        "First-order product formula assumed for the 2028 step (chapter states no Trotter order).",
        f"Coherent synthesis (app06:109): the same {eps_syn:g} total summed linearly needs {eps_coh:.3g} per rotation; "
        f"the RUS fit gives {t_rot_coh_model:.2f} T ('~31') for every rotation; step "
        f"{t_step_coh_model:.0f} T, shot {t_shot_coh_model:.4g} T unchunked ({t_shot_coh_chunked:.4g} chunked).",
        f"1+1D overlap (r22, derived here; single-particle architecture of arxiv_2607_28524): index register "
        f"{spq1['index_register']} qubits, 5 LCU terms, query {spq1['toffoli']} Toffolis + {spq1['n_rot']} rotations = "
        f"c_sp {ov_a1['c_sp']:.0f} T, lift {lift1['toffoli']} Toffolis. One application at kappa 30: "
        f"{ov_a1['d_sgn'] + 1} queries + lift = {ov_a1['t']:.4g} T ({ov_a1['t'] / env:.2f}x the cap); a ground-state "
        f"preparation of {n_gs}: {ov_a10['t']:.4g} T, {ov_a10['t'] / env:.1f}x the cap and "
        f"{ov_a10['t'] / t_shot_chunked:.1f}x the domain-wall shot. At kappa 50: {ov_b1['t']:.4g} and {ov_b10['t']:.4g} T. "
        f"Record, second-quantized pricing (r21 e): c_be {ov_sq['c_be']:.0f}, {ov_sq['sgn_T_kappa30']:.4g} and "
        f"{ov_sq['gs_prep_T']:.4g} T.",
        f"Shots (r21 d, derived here): sigma^2 = {float(a.sigma2_condensate.value):g} x {copies_2033} / {copies_2028} = "
        f"{sigma2_2028:.4g} -> {shots:.0f} shots; fully correlated sites would give {sigma2_worst:.0f} -> {shots_worst:.0f}.",
        f"eps_l at {faults} faults per shot (R3): {eps_l[0]:.2g} (chunked) - {eps_l[1]:.2g} (unchunked).",
    )
    return Result(era="2028", lq=(lq, lq), hard_ops=(t_shot, t_shot_chunked), breakdown=tuple(breakdown),
                  intermediates=inter, shots=(shots, shots),
                  wall_time_s=inter["wall_campaign_s"],   # r25: the campaign wall, as in Chs. 3-5, 7, 10
                  epsilon_l=eps_l, notes=notes)


# --------------------------------------------------------------------------- #
# 2033: SU(2)/2O overlap, QSP sign function, TPQ thermal state
# --------------------------------------------------------------------------- #

def _model_2033(a: Assumptions) -> Result:
    g = GROUPS[a.gauge_group_2033.value]
    D, L = int(a.D_2033.value), int(a.L_2033.value)
    Nf, Nc = int(a.Nf_2033.value), int(a.Nc_2033.value)
    V = L ** D                                        # 27
    dirac = 2 ** ((D + 1) // 2)                       # 4 in 3+1D
    n_rho_q = math.ceil(math.log2(int(a.N_rho.value)))   # 4
    lq_f = dirac * Nf * Nc * V                        # 216
    lq_g = D * V * g.link_qubits                      # 486
    lq_h = V * (n_rho_q + g.link_qubits)              # 270
    tc = a.toffoli_convention.value
    eps_syn = float(a.eps_syn.value)
    n_links = D * V                                   # 81
    n_spin = int(a.hop_n_spin.value)                  # 2
    frame = a.link_frame.value                        # record parameter of the second-quantized construction
    arch = a.sgn_architecture.value                   # 'single_particle_lift' (E26 (i))
    fd = FERMION_DRAFT[a.gauge_group_2033.value]

    # r22 (E26 (i)): the single-particle kernel. H = psi^dag h_ov[U] psi, h_ov = gamma^0 + sgn(H_W[U]); H_W[U] is a
    # Q x Q matrix of link operators, promoted to an index register (color, spin, coordinates) and the link registers.
    q_modes = lq_f                                    # Q = 216 single-particle modes = the Jordan-Wigner register
    q_index = {"color": math.ceil(math.log2(Nc * Nf)), "spin": math.ceil(math.log2(dirac)),
               "coordinates": D * math.ceil(math.log2(L))}                # 1 + 2 + 6 = 9
    n_index = sum(q_index.values())
    m0 = float(a.m0_wilson.value)
    alpha_hw = 2.0 * D + abs(D - m0)                  # r = 1: hops D (Dirac) + D (Wilson), on-site |D - m_0|: 7.5
    norm_hw_free = 2.0 * D - m0                       # the free-field ||H_W|| (all momenta at pi): 4.5
    be_norm = 2 * q_modes                             # normalization of the lifted block encoding: Q (lift) x 2 (h_ov)

    # The Higgs vertex's SELECT strings (r21 c): the on-site bilinear (4 N Z strings at d = 3) weighted by rho, whose
    # binary expansion on log2 N_rho qubits has log2 N_rho + 1 terms. The hop and on-site counts of the same function
    # are the LCU of the retired second-quantized construction (record).
    lcu = lcu_terms_per_site(D, n_spin, Nc, Nf, int(a.lcu_strings_per_bilinear.value), radial_terms=n_rho_q + 1)
    lcu_terms = V * lcu["total"]                      # 27 x 72 = 1944 (record)
    lcu_toff = lcu_terms * float(a.lcu_toffoli_per_term.value)
    c_toff = toffoli_t(lcu_toff, tc)                  # 13,608 (record)
    lcu_terms_pre = float(a.lcu_terms_stated_pre_r21.value)        # 1.7e3 (record)
    c_toff_pre = toffoli_t(lcu_terms_pre, tc)         # 11,900 (record)
    higgs_strings = V * lcu["higgs"]                  # 27 x 40 = 1080

    # QSP sign-function degree: natural log (log2 would give ~1e3). A Chebyshev degree is an
    # integer, so ceil; the chapter rounds 691 -> 700, carried as d_sgn_stated.
    # r22: kappa = alpha / Delta, alpha the normalization of the block encoding of H_W (arxiv_2607_28524).
    kappa = float(a.kappa_2033.value)
    delta = float(a.delta_sgn.value)
    d_sgn_exact = kappa * math.log(1 / delta)
    d_sgn_log2 = kappa * math.log2(1 / delta)
    d_sgn = math.ceil(d_sgn_exact)                    # 208 at kappa = 30 (r24); 691 at kappa = 1e2 (r22-r23)
    # r24 record inputs: the r17-r23 kappa (1e2) and evolution time (50a) price the retired records and the r22 headline
    kappa_r22 = float(a.kappa_2033_r22.value)
    d_sgn_r22 = math.ceil(kappa_r22 * math.log(1 / delta))       # 691
    d_sgn_chapter = float(a.d_sgn_stated.value)       # 700
    d_sgn_norm_convention = math.ceil(d_sgn_exact * alpha_hw / norm_hw_free)   # 1152: kappa read as ||H_W|| / Delta

    # ---- the Trotter step count -------------------------------------------------------------------------------
    t_over_a = float(a.t_evol_fm.value) / float(a.a_fm.value)    # 10 (r24; was 50)
    t_over_a_r22 = float(a.t_evol_fm_r22.value) / float(a.a_fm.value)   # 50 (record)
    eps_tr = float(a.eps_trotter.value)               # 0.1
    # r22 (E26 (ii)), headline: the state-dependent second-order estimate in the thermal state (weak coupling)
    n_adj = Nc * Nc - 1                               # 3 adjoint colors
    sd = trotter_state_dependent_steps(L, n_adj, t_over_a, eps_tr, float(a.temp_lat_2033.value), d=D)
    # r21 (3b), now the quoted worst case: the second-order commutator bound with operator norms
    hop_norm = float(Nc * Nf * n_spin)                # 4: one unit-norm bilinear per hopping (spinor, color) mode, r = 1
    n_plaq_per_link = 2 * (D - 1)                     # 4
    g2 = float(a.g2_bare.value)
    cas = float(a.casimir_max_2O.value)
    tb = trotter_bound_steps(n_links, t_over_a, eps_tr, g2, cas, hop_norm, n_plaq_per_link)
    tb_r22 = trotter_bound_steps(n_links, t_over_a_r22, eps_tr, g2, cas, hop_norm, n_plaq_per_link)   # 83,195 (record)
    sd_r22 = trotter_state_dependent_steps(L, n_adj, t_over_a_r22, eps_tr, float(a.temp_lat_2033.value), d=D)   # 5031
    g2_scan = [(x, trotter_bound_steps(n_links, t_over_a, eps_tr, x, cas, hop_norm, n_plaq_per_link))
               for x in (0.2 + 0.01 * i for i in range(581))]              # g^2 in [0.2, 6]
    g2_at_min, tb_min = min(g2_scan, key=lambda xr: xr[1]["lambda"])       # g^2 = 1.86: Lambda 4773, 77,244 steps
    n_trot_pre = int(a.n_trotter_2033_pre_r21.value)  # 30 (record)
    lambda_for_pre = eps_tr * n_trot_pre ** 2 / t_over_a_r22 ** 3   # 7.2e-4 a^-3: what a bound would need for 30 steps
    n_trot_step = math.ceil(round(float(a.t_evol_fm.value) / float(a.dt_empirical_fm.value), 9))   # 1000 (Ch. 5's step)
    n_trot_by_rule = {"state": sd["n_central"], "bound": tb["n_steps"], "step": n_trot_step, "stated": n_trot_pre}
    n_trot = n_trot_by_rule[a.trotter_rule_2033.value]            # 5031
    n_trot_range = (sd["n_low"], sd["n_high"])        # 1186 (fastest-mode phase), 14,242 (mean kept)
    eps_qsp = float(a.eps_qsp.value)

    # TPQ thermal state: QSVT filter of degree d_beta, carried by aa_rounds rounds.  The
    # sqrt formula is a scaling estimate ('~30'), so nearest integer: sqrt(1e2 ln 1e4) = 30.35 -> 30.
    # Referee C1 (2026-10-04): beta||H|| from the block encoding the filter uses, 2Q / (T a) = 864 (r25 record 1e2).
    beta_H_used = (be_norm / float(a.temp_lat_2033.value) if a.beta_H_rule.value == "normalization"
                   else float(a.beta_H.value))       # 864
    d_beta_exact = math.sqrt(beta_H_used * math.log(1 / float(a.delta_beta.value)))
    d_beta = round(d_beta_exact)                      # 89 (r25: 30)
    d_beta_at_1e3 = math.sqrt(beta_H_used * math.log(1e3))   # 77.3, the other candidate (r25 26.3)
    d_beta_record_r25 = round(math.sqrt(float(a.beta_H.value) * math.log(1 / float(a.delta_beta.value))))   # 30
    # k rounds of (fixed-point) amplitude amplification call the filter and its inverse once each per round, plus the
    # initial call: 2k + 1 = 17 (verifier 2026-10-04). r25 and the first referee pass charged one call per round.
    n_filter_calls = (2 * int(a.aa_rounds.value) + 1 if a.aa_call_rule.value == "2k+1"
                      else int(a.aa_rounds.value))    # 17 (record 'per_round': 8)
    n_dov_tpq = d_beta * n_filter_calls               # 1513 (first referee pass 712; r25 240)
    n_dov_tpq_one_call_per_round = d_beta * int(a.aa_rounds.value)   # 712: RECORD of the first referee pass
    n_dov_tpq_r25 = d_beta_record_r25 * int(a.aa_rounds.value)   # 240: the records below keep the r25 filter

    # ---- one application of the block-encoded psi^dag h_ov psi / (2Q) (r22, DERIVED HERE) ----------------------
    spq = single_particle_query_2O(D, L, g.link_qubits)            # 598 Toffolis + 24 T + 6 rotations
    mux_ix = colour_multiplexer_index_2O()                         # 10 T, no Toffoli, no ancilla
    lift = lift_cost(q_modes, D if L == 3 else 0)                  # 862 Toffolis + 12 T + 6 rotations (3 coordinates mod 3)
    # r25 (a): the Ginsparg-Wilson-projected Yukawa is inside. Its projector (1 - sgn(H_W))/2 is a second QSVT
    # polynomial of degree d_sgn on the same block encoding (same signal and real-part qubits, run after the first),
    # once per application; the lift and the Higgs vertex stay once per application.
    n_sgn = 2 if a.yukawa_placement.value == "inside" else 1
    n_q_per_dov = (d_sgn + 1) + (n_sgn - 1) * d_sgn  # 417 (r24 209): d_sgn + 1 for h_ov, d_sgn for the projector
    # the QSP that evolves or filters with the block encoding reflects about the all-zero state of its ancilla
    n_be_anc = (int(a.anc_lift_type.value) + n_index + int(a.anc_dov_lcu.value) + int(a.anc_sgn_signal.value)
                + int(a.anc_sgn_real_part.value) + spq["term_register"])       # 2 + 9 + 1 + 1 + 1 + 5 = 19
    outer = {"toffoli": n_be_anc - 1, "n_rot": 1}                  # 18 Toffolis + the QSP phase
    # The Higgs vertex is not part of H_W: outside the sign function, once per application.
    umult = g.primitives["U_mul"]
    t_umult = umult.t(1e-4)                           # 392 (no rotations: eps-free)
    umult_toffolis = umult.toffoli                    # 56
    umult_ancilla = umult.ancilla                     # 4 clean
    t_umult_report = umult.t_report(1e-4, tc)         # 392: same number in the report convention
    n_um = int(a.n_umult_per_higgs_site.value)        # 1
    c_higgs = V * n_um * t_umult                      # 10,584
    higgs_sel_toff = higgs_strings * float(a.lcu_toffoli_per_term.value)       # 1080
    c_higgs_app = c_higgs + toffoli_t(higgs_sel_toff, tc)          # 18,144

    def _price_sp(n_trot_, d_=d_sgn, per_step_=None, t_=t_over_a, toff_saved_=0, n_sgn_=n_sgn, tpq_=n_dov_tpq):
        """One shot in the single-particle architecture: n_trot_ steps of r applications (r = the Jacobi-Anger order
        at 2Q dt unless given) + the TPQ applications, each (d_ + 1) queries + lift + reflection + Higgs, at this
        circuit's own R-TOL tolerance. r23 levers: t_ = the evolution time in lattice units, toff_saved_ = Toffolis
        removed from one query."""
        if n_trot_ == 0:                              # r24: preparation only (the condensate arm), no evolution
            r_ = 0
        else:
            r_ = dov_per_step(n_trot_, t_, q_modes, eps_qsp) if per_step_ is None else int(per_step_)
        n_dov_ = n_trot_ * r_ + tpq_
        nq_ = (d_ + 1) + (n_sgn_ - 1) * d_
        n_rot_app_ = nq_ * spq["n_rot"] + lift["n_rot"] + outer["n_rot"]
        n_rot_shot_ = n_dov_ * n_rot_app_
        e = eps_rot_for(n_rot_shot_, eps_syn)
        tr = t_per_rotation(e, a.rot_synthesis.value)
        c_sp_ = toffoli_t(spq["toffoli"] - toff_saved_, tc) + spq["t_direct"] + spq["n_rot"] * tr
        c_lift_ = toffoli_t(lift["toffoli"], tc) + lift["t_direct"] + lift["n_rot"] * tr
        c_outer_ = toffoli_t(outer["toffoli"], tc) + outer["n_rot"] * tr
        t_dov_ = nq_ * c_sp_ + c_lift_ + c_outer_ + c_higgs_app
        return {"n_trot": n_trot_, "per_step": r_, "x": be_norm * t_ / n_trot_ if n_trot_ else 0.0, "n_dov": n_dov_,
                "n_q_per_dov": nq_, "t_lat": t_,
                "n_dov_evol": n_trot_ * r_, "n_queries": n_dov_ * nq_, "n_rot_app": n_rot_app_,
                "n_rot_shot": n_rot_shot_, "eps_rot": e, "t_rot": tr, "c_sp": c_sp_, "c_lift": c_lift_,
                "c_outer": c_outer_, "t_sgn": nq_ * c_sp_, "t_dov": t_dov_, "t_shot": n_dov_ * t_dov_}

    P = _price_sp(n_trot)                             # the headline
    by_rule = {k: _price_sp(v) for k, v in n_trot_by_rule.items()}       # the step count's alternatives
    P_low, P_high = _price_sp(n_trot_range[0]), _price_sp(n_trot_range[1])
    per_step = P["per_step"]                          # 13
    n_dov_evol, n_dov, n_queries = P["n_dov_evol"], P["n_dov"], P["n_queries"]
    n_rot_shot, eps_rot, t_rot = P["n_rot_shot"], P["eps_rot"], P["t_rot"]
    c_sp, t_dov, t_shot = P["c_sp"], P["t_dov"], P["t_shot"]
    t_rot_retired = float(a.t_per_rot.value)          # 30: retired constant, records only
    n_dov_floor = be_norm * t_over_a                  # 2Q t = 4320 applications (r22: 21,600): no step size goes below
    t_shot_floor = n_dov_floor * t_dov                # 4.1e9 T at the headline application price
    # r24 records: the r22-r23 headline (kappa 1e2, t = 50a, 5031 steps) and the worst-case step at those inputs
    P_r22 = _price_sp(sd_r22["n_central"], d_=d_sgn_r22, t_=t_over_a_r22, n_sgn_=1, tpq_=n_dov_tpq_r25)      # 2.00789e11
    P_r22_bound = _price_sp(tb_r22["n_steps"], d_=d_sgn_r22, t_=t_over_a_r22, n_sgn_=1, tpq_=n_dov_tpq_r25)  # 1.276e12
    # r25 records: the r24 headline (Yukawa outside sgn) and its condensate shot
    P_r24 = _price_sp(n_trot, n_sgn_=1, tpq_=n_dov_tpq_r25)               # 8.2414e9
    C_r24 = _price_sp(0, n_sgn_=1, tpq_=n_dov_tpq_r25)                    # 2.2412e8
    P_r25 = _price_sp(n_trot, tpq_=n_dov_tpq_r25)     # 1.6242e10: the r25 headline (filter at the record 1e2)
    C_r25 = _price_sp(0, tpq_=n_dov_tpq_r25)          # 4.4166e8: the r25 condensate shot
    P_kappa_r22 = _price_sp(n_trot, d_=d_sgn_r22)     # kappa ~ 1e2 at t = 10a: the sensitivity the prose prints
    # r24 (1): the condensate arm. <psibar psi> is static: measured on the prepared state, no real-time evolution, so its
    # shot is the preparation alone (the TPQ filter, n_dov_tpq applications) at that circuit's own R-TOL tolerance.
    C = _price_sp(0)
    P_not_self_inverse = _price_sp(n_trot, per_step_=2 * per_step)       # twice the order without a self-inverse BE
    P_kappa_norm = _price_sp(n_trot, d_=d_sgn_norm_convention)           # kappa read as ||H_W|| / Delta
    link_read_share = toffoli_t(spq["link_read_toffoli"], tc) / c_sp     # 90% of a query
    rot_share = spq["n_rot"] * t_rot / c_sp                              # 4% of a query
    # NOT PRICED (ESTIMATED HERE, record): the gauge terms of one Trotter step, per link the Kogut-Susskind
    # multiplicities of arXiv:2312.10285 (groups.PRIMCOST) with the group Fourier transform of arXiv:2408.00075 in place
    # of U_F, at the headline tolerance. The Higgs terms are not estimated.
    t_gauge_link = (PRIMCOST["KS"]["U_F"](D) * g.primitives["U_FFT"].t(eps_rot)
                    + sum(PRIMCOST["KS"][k](D) * g.primitives[k].t(eps_rot) for k in ("U_Tr", "U_inv", "U_mul")))
    t_gauge_step = n_links * t_gauge_link             # 7.1e5 T
    gauge_step_share = t_gauge_step / (per_step * t_dov)                 # about 2% of the step
    # NOT PRICED (record): the amplification's reflection about the reference state on the whole register, as an
    # n-controlled NOT at 2n Toffolis on two clean ancilla (arxiv_2407_17966, abstract), twice per round.
    n_reflect_controls = lq_f + lq_g + lq_h           # 972
    n_reflections = 2 * int(a.aa_rounds.value)        # 16
    t_fpaa_reflections = n_reflections * toffoli_t(2 * n_reflect_controls, tc)       # 2.2e5 T per shot

    # ---- RECORD: the retired second-quantized construction (r17-r21). QSVT of sgn on a Jordan-Wigner block encoding of
    # H_W computes the sign of the many-body operator, not the kernel of the chapter's Hamiltonian; kept so every earlier
    # print stays reproducible. Each of the 81 link terms is the Fock-space 2O multiplexer (r20) plus the draft's spinor
    # frame, 468 T; the alternative frame is the draft's compiled hop (groups.hop_link_cost). 5 D_ov per step (Stated).
    per_step_sq = int(a.n_dov_per_step_bound_stated.value)        # 5
    n_rot_be = float(a.n_rot_be.value)                # 25 (PREPARE + QSP phases)
    hop_kw = dict(fermion=a.hop_fermion_2033.value, undo=bool(a.hop_frame_undo.value), share=a.hop_share.value,
                  mcx=a.hop_mcx.value, n_spin=n_spin, phasing=a.hop_phasing_be.value,
                  synthesis=a.rot_synthesis.value, toffoli_convention=tc, legacy=False)

    def _hop(eps, **over):
        return hop_link_cost(a.gauge_group_2033.value, Nf, D, eps=eps, **{**hop_kw, **over})
    hop0 = _hop(1e-4)                                 # counts only (eps-free)
    n_rot_link_draft = hop0["n_rot"]                  # 704 = 480 color + 224 diagonalizer
    spinor_item = next(x for x in hop0["items"] if x["name"] == "spinor_frame")
    mux = colour_multiplexer_2O()
    ml = multiplexer_link_2O(n_spin, spinor_item)     # 52 Toffoli + 104 T
    t_link_mux = toffoli_t(ml["toffoli"], tc) + ml["t_direct"]     # 468
    # record: the draft's own Wilson hop W2 = 4 C^S + 148 N per link (section_resources.tex:50-53), the r17-r20 band top
    cs_a, cs_b = fd.printed["C_S"]
    w2_mult = int(a.w2_cs_multiplicity.value)
    w2_rot_link = round(w2_mult * cs_b / 1.15)                                         # 928
    w2_t_link = w2_mult * (cs_a - 9.2 * cs_b / 1.15) + float(a.w2_spinor_t_per_colour.value) * Nc   # 28,696

    def _price_sq(n_trot_, frame_, d_=d_sgn_r22, c_toff_=c_toff, **hop_over):
        """RECORD. One shot of the second-quantized construction: N_Dov = n_trot x 5 + TPQ queries, each d_ BE queries,
        with the link frame `frame_` ('multiplexer' | 'draft' | 'w2' | 'zero'), at this circuit's own R-TOL tolerance."""
        n_dov_ = n_trot_ * per_step_sq + n_dov_tpq_r25
        n_q = n_dov_ * d_
        n_rot_link_ = {"multiplexer": 0, "zero": 0, "w2": w2_rot_link}.get(frame_)
        if n_rot_link_ is None:
            n_rot_link_ = _hop(1e-4, **hop_over)["n_rot"]
        n_rot_q = n_rot_be + n_links * n_rot_link_
        e = eps_rot_for(n_rot_q * n_q, eps_syn)
        tr = t_per_rotation(e, a.rot_synthesis.value)
        hop_ = _hop(e, **hop_over) if frame_ == "draft" else None
        t_link_ = {"multiplexer": t_link_mux, "zero": 0.0, "w2": w2_t_link + w2_rot_link * tr}.get(frame_)
        if t_link_ is None:
            t_link_ = hop_["t"]
        c_higgs_ = 0.0 if frame_ == "zero" else c_higgs
        c_be_ = c_toff_ + n_links * t_link_ + c_higgs_ + n_rot_be * tr
        return {"n_trot": n_trot_, "frame": frame_, "n_dov": n_dov_, "n_dov_evol": n_trot_ * per_step_sq, "n_queries": n_q,
                "n_rot_link": n_rot_link_, "n_rot_query": n_rot_q, "n_rot_shot": n_rot_q * n_q, "eps_rot": e, "t_rot": tr,
                "t_link": t_link_, "hop": hop_, "c_hop": n_links * t_link_, "c_rot_be": n_rot_be * tr, "c_be": c_be_,
                "t_dov": d_ * c_be_, "t_shot": n_q * c_be_}

    SQ = _price_sq(tb_r22["n_steps"], frame)              # the r21 headline: c_be 62,896.2, 1.809e13 T per shot
    alt_frame = "draft" if frame == "multiplexer" else "multiplexer"
    SQ_alt = _price_sq(tb_r22["n_steps"], alt_frame)
    SQ_draft = SQ if frame == "draft" else SQ_alt
    SQ_mux = SQ if frame == "multiplexer" else SQ_alt
    hop = SQ_draft["hop"]                             # the draft hop at its own tolerance
    t_link_clifford_t = hop["t_toffoli"] + hop["t_direct"]         # 2832 = 7 x 384 + 144
    n_trot_by_rule_r22 = {"bound": tb_r22["n_steps"],
                          "step": math.ceil(round(float(a.t_evol_fm_r22.value) / float(a.dt_empirical_fm.value), 9)),
                          "stated": n_trot_pre}
    sq_by_rule = {k: _price_sq(n_trot_by_rule_r22[k], frame) for k in ("bound", "step", "stated")}

    # ---- ancilla, itemized at the peak of the schedule (E23 (5), r20; r22 for the single-particle architecture) ----
    # Allocated for the whole shot: the readout control (used only after the state is prepared), one QSP signal qubit
    # shared by the TPQ filter and the evolution, the amplification marker, the LCU index over the parts of H.
    # Held through one application: the lift's type register A and the index register B, the LCU qubit of h_ov, the
    # signal and real-part qubits of sgn(H_W). One query of H_W adds its term register. Inside SELECT the peak is the
    # link read: direction flags, the unary iteration's AND lines, the copy of the addressed link.
    # Not concurrent with the link read: the orientation / Dirac flag (one at a time), the M ANDs, the lift's own
    # iteration (outside the sgn sequence), the Higgs U_x, repeat-until-success synthesis.
    anc_rus = int(a.anc_rus.value)                    # 1 (E26)
    anc_stack = {"readout": int(a.anc_readout.value), "qsp_signal": int(a.anc_qsp_signal.value),
                 "fpaa": int(a.anc_fpaa.value), "h_lcu": int(a.anc_h_lcu.value)}                  # 5
    anc_app = {"lift_type_register": int(a.anc_lift_type.value), "index_register": n_index,
               "hov_lcu": int(a.anc_dov_lcu.value), "sgn_signal": int(a.anc_sgn_signal.value),
               "sgn_real_part": int(a.anc_sgn_real_part.value)}                                   # 14
    anc_query = {"term_register": spq["term_register"]}                                           # 5
    anc_select = {"direction_flags": spq["anc_direction_flags"], "iteration_temporaries": spq["anc_iteration"],
                  "link_copy": spq["anc_link_copy"]}                                              # 14
    anc_not_concurrent = {"lift_iteration": lift["anc_iteration"], "higgs_umult": umult_ancilla, "rus": anc_rus,
                          "fpaa_reflection": int(a.fpaa_reflection_ancilla.value)}                # 10, 4, 1, 2
    assert max(anc_not_concurrent.values()) <= sum(anc_query.values()) + sum(anc_select.values())
    lq_anc_itemized = (sum(anc_stack.values()) + sum(anc_app.values()) + sum(anc_query.values())
                       + sum(anc_select.values()))                                  # 38
    # r25 (R10): the CNOT fan-out of the shared control during the link read restores F* >= 10 (factory analysis)
    depth_variant = a.depth_variant_2033.value
    anc_fanout = int(a.anc_fanout_2033.value) if depth_variant == "fanout" else 0     # 5
    lq_anc = lq_anc_itemized + anc_fanout             # 43 (r22-r24 38)
    lq_anc_asserted = int(a.ancilla_2033_asserted.value)                         # 110 (record)
    lq = lq_f + lq_g + lq_h + lq_anc                  # 1015 (r22-r24 1010, r21 1005, r20 1015, r19 1082)
    lq_dwf_backup = tuple(lq_f * L5 + lq_g + lq_h + lq_anc for L5 in (int(a.L5_backup.lo), int(a.L5_backup.hi)))
    # record: the second-quantized construction's itemization (r21: 7 + 22 + 4 = 33; draft frame 43)
    n_addr = math.ceil(math.log2(lcu_terms))                                     # 11
    sq_anc_stack = {"readout": int(a.anc_readout.value), "tpq_signal": int(a.anc_qsp_signal.value),
                    "fpaa": int(a.anc_fpaa.value), "h_lcu": int(a.anc_h_lcu.value), "dov_lcu": int(a.anc_dov_lcu.value),
                    "sgn_signal": int(a.anc_sgn_signal.value)}                    # 7
    sq_anc_be = {"lcu_address": n_addr, "unary_and_ladder": n_addr - 1,
                 "unary_control": int(a.anc_unary_control.value)}                 # 22
    sq_anc_term_mux = {"link_multiplexer": mux["ancilla"], "spinor_frame": int(a.anc_spinor.value), "rus": anc_rus,
                       "higgs_umult": umult_ancilla}                              # 1, 2, 1, 4 -> 4
    mcx_tmp = max(max(fd.col_squish_mcx), max(fd.hop_squish_mcx)) - 2              # C^6X under MBU: 4 temporaries
    transient = {"mcx_ladder": mcx_tmp, "rus": anc_rus, "spinor_frame": int(a.anc_spinor.value)}
    anc_link_draft = {"pq_registers": int(a.anc_link_pq.value), "hop_eigen_register": int(a.anc_link_hop_register.value),
                      "parity": int(a.anc_link_parity.value), "largest_transient": max(transient.values())}   # 14
    sq_anc_base = sum(sq_anc_stack.values()) + sum(sq_anc_be.values())           # 29
    sq_lq_anc_mux = sq_anc_base + max(sq_anc_term_mux.values())                  # 33
    sq_lq_anc_draft = sq_anc_base + max(sum(anc_link_draft.values()), umult_ancilla)   # 43
    sq_lq_anc = sq_lq_anc_mux if frame == "multiplexer" else sq_lq_anc_draft
    sq_lq = lq_f + lq_g + lq_h + sq_lq_anc            # 1005 (r21 box)
    sq_lq_draft = lq_f + lq_g + lq_h + sq_lq_anc_draft   # 1015
    qubits_per_site = D * g.link_qubits + (n_rho_q + g.link_qubits) + dirac * Nf * Nc   # 36

    # Records at the pre-r21 query count (30 x 5 + 240 = 390 D_ov, 691 each) and the pre-r21 1.7e3-term LCU line, so
    # every earlier print stays reproducible: the r20 center (draft frame), the r20 floor (multiplexer), the W2 top,
    # the zero-T floor, and the older centers. All of them price the retired second-quantized construction.
    R_draft = _price_sq(n_trot_pre, "draft", c_toff_=c_toff_pre)       # c_be 2,104,718.8; 5.672e11
    R_mux = _price_sq(n_trot_pre, "multiplexer", c_toff_=c_toff_pre)   # 61,043.6; 1.645e10
    R_w2 = _price_sq(n_trot_pre, "w2", c_toff_=c_toff_pre)             # 4,806,210.4; 1.2952e12
    R_zero = _price_sq(n_trot_pre, "zero", c_toff_=c_toff_pre)         # 12,551.6; 3.383e9
    n_queries_pre = R_draft["n_queries"]              # 269,490

    def _c_be_vertex(n_per_link):
        e = eps_rot_for(n_rot_be * n_queries_pre, eps_syn)
        return c_toff_pre + (n_links * n_per_link + V) * t_umult + n_rot_be * t_per_rotation(e, a.rot_synthesis.value)
    c_be_records = {"vertex_1_umult_pre_r17": _c_be_vertex(1),          # 54,887.6
                    "rhop_2_umult_per_link": _c_be_vertex(2),            # 86,639.6
                    "draft_w1_no_undo": _price_sq(n_trot_pre, "draft", c_toff_=c_toff_pre, undo=False, share="draft")["c_be"],
                    "share_draft_r17_centre": _price_sq(n_trot_pre, "draft", c_toff_=c_toff_pre, share="draft")["c_be"],
                    "share_link_no_undo": _price_sq(n_trot_pre, "draft", c_toff_=c_toff_pre, undo=False, share="link")["c_be"],
                    "w1u_n_spin_1": _price_sq(n_trot_pre, "draft", c_toff_=c_toff_pre, n_spin=1)["c_be"],
                    "w1u_mcx_2n_3": _price_sq(n_trot_pre, "draft", c_toff_=c_toff_pre, mcx="2n-3")["c_be"]}
    t_shot_records = {k: n_queries_pre * v for k, v in c_be_records.items()}
    t_shot_records.update({"r20_centre_draft_frame": R_draft["t_shot"], "r20_floor_multiplexer": R_mux["t_shot"],
                           "r20_top_draft_w2": R_w2["t_shot"], "pre_r20_floor_zero_T": R_zero["t_shot"]})
    # RECORD (r19 prose, retired r20): the draft's frame with its rotation synthesis removed. Not a realizable circuit:
    # the draft diagonalizes U_g, and 2O eigenphases include e^{+-i pi/3}, e^{+-2i pi/3} (outside Z[1/sqrt2, i]).
    c_be_rotation_free = R_zero["c_be"] + n_links * t_link_clifford_t + c_higgs
    # For the record: the pre-round-B center (4-T Toffolis, the retired unsourced 190-T vertex),
    # read from groups.py's UNSOURCED record so the retired number lives in one place.
    retired = g.primitives["controlled_group_action"]
    n_ga = n_links + V                                # 108 vertices of the retired count
    c_be_before_round_b = (toffoli_t(lcu_terms_pre, "jones") + n_ga * retired.t_const
                           + n_rot_be * t_rot_retired)                                     # 28070
    t_shot_before_round_b = n_queries_pre * c_be_before_round_b     # 7.56e9, the old '~8e9'
    second_quantized_record = {
        "valid": False, "c_be": SQ["c_be"], "t_per_dov": SQ["t_dov"], "t_per_shot": SQ["t_shot"],
        "n_dov": SQ["n_dov"], "n_queries": SQ["n_queries"], "n_rot_shot": SQ["n_rot_shot"], "t_per_rot": SQ["t_rot"],
        "lq_ancilla": sq_lq_anc, "lq_total": sq_lq, "per_step": per_step_sq,
        "t_per_shot_by_trotter_rule": {k: v["t_shot"] for k, v in sq_by_rule.items()},
        "over_single_particle_per_application": SQ["t_dov"] / P_r22_bound["t_dov"],     # 14.2 at equal step count
        "over_single_particle_per_shot": SQ["t_shot"] / P_r22_bound["t_shot"]}

    # coherent-synthesis shift of one application: the same eps_syn summed linearly over this circuit's own rotations
    eps_coh = eps_per_rotation(eps_syn, n_rot_shot, "coherent")
    t_rot_coh = t_per_rotation(eps_coh, a.rot_synthesis.value)
    t_dov_shift_coh = P["n_rot_app"] * (t_rot_coh - t_rot) / t_dov  # 2%

    t_evol = n_dov_evol * t_dov                       # 2.0006e11
    t_tpq = n_dov_tpq * t_dov                         # 7.341e8
    env = float(a.envelope_2033.value)
    ratio_env = t_shot / env                          # 200.8

    # R-TOL: each kappa is its own circuit, so its own rotation count sets its own tolerance
    t_phys = tuple(_price_sp(n_trot, d_=math.ceil(k * math.log(1 / delta)))["t_shot"]
                   for k in (float(a.kappa_physical.lo), float(a.kappa_physical.hi)))
    faults = float(a.faults_per_shot.value)           # R3: 0.1 expected faults per shot
    eps_l_implied = faults / t_shot                   # 4.98e-13

    # thermal-route alternatives quoted in the prose (app06 route (iv))
    dS_haar = float(a.haar_deficit_frac.value) * float(a.haar_register_qubits.value) * math.log(2)
    haar_rounds = math.exp(dS_haar / 2)               # RECORD (r25 prose '~3e7', retired by the chapter-open pass)
    # chapter-open pass (referee C1, DERIVED HERE): typical (1-design) reference on the free fermion sector of 3^3
    temp = float(a.temp_lat_2033.value)
    ov_modes = free_overlap_modes(L, D, -m0, 1.0) * (Nc * Nf)                  # 216 modes
    typ = typical_reference_tpq(ov_modes, 1.0 / temp, float(be_norm))         # ln p = -137.9, rounds 7.0e29
    # referee C3 (DERIVED HERE): Delta N_CS by the Kubo form, Hadamard-test readout, free-field background
    g2_ncs = float(a.g2_bare.value)
    ncs_bg_helicity = {tt: ncs_free_field_background(L, n_adj, temp, g2_ncs, tt, d=D) for tt in (5.0, 10.0)}  # 3.7e-3, 8.9e-3 (record)
    ncs_bg = {tt: ncs_lattice_background(L, n_adj, temp, g2_ncs, tt, "avg") for tt in (5.0, 10.0)}  # 4.3e-4, 9.58e-4
    ncs_bg_fwd = {tt: ncs_lattice_background(L, n_adj, temp, g2_ncs, tt, "fwd") for tt in (10.0, 40.0)}  # 2.56e-3, 2.42e-2
    ncs_lam = ncs_hadamard_lambda(V, D, n_adj, g2_ncs, float(a.ncs_e_max.value), float(a.ncs_b_max.value))  # >= 9.23
    ncs_relvar = {tt: ncs_hadamard_rel_var(ncs_lam, tt, ncs_bg[tt]) for tt in ncs_bg}               # 3.3e13, 1.06e14
    ncs_relvar_helicity = {tt: ncs_hadamard_rel_var(ncs_lam, tt, ncs_bg_helicity[tt]) for tt in ncs_bg_helicity}  # 1.2e12 at 10a
    ncs_relvar_fwd10 = ncs_hadamard_rel_var(ncs_lam, 10.0, ncs_bg_fwd[10.0])                       # 1.5e13
    ncs_relvar_unit_signal = ncs_hadamard_rel_var(ncs_lam, 10.0, 1.0)                              # 9.7e7
    # author ruling 2026-10-05 (method test, no rate at 3^3): what the full <dN^2> is compared with
    ncs_bg_classical = {tt: ncs_lattice_background(L, n_adj, temp, g2_ncs, tt, "avg", statistics="classical")
                        for tt in (5.0, 10.0)}                                                     # 6.2e-5, 1.34e-4
    ncs_bg_zero_point10 = ncs_lattice_background(L, n_adj, 0.01, g2_ncs, 10.0, "avg")               # 9.07e-4 (95%)
    ncs_bg_temp_span10 = (ncs_lattice_background(L, n_adj, 0.4, g2_ncs, 10.0, "avg") / ncs_bg[10.0],
                          ncs_lattice_background(L, n_adj, 0.6, g2_ncs, 10.0, "avg") / ncs_bg[10.0])  # 0.968, 1.047
    # referee C3 (DERIVED HERE): dimensionless lattice units and the classical small-box suppression
    beta_lattice = 4.0 / (g2_ncs * temp)              # 8: Moore-Rummukainen's a = 1/(2 g^2 T) lattice
    box_lg2t = L * temp * g2_ncs                      # 1.5 (0.63 at g^2 = 0.42)
    mr_suppression = float(a.mr_rate_box10.value) / float(a.mr_rate_box3.value)   # 730 at L g^2 T = 3 (intermediate)
    mr_suppression_range = (float(a.mr_rate_box10.value) / (float(a.mr_rate_box3.value) + float(a.mr_rate_box3_err.value)),
                            float(a.mr_rate_box10.value) / (float(a.mr_rate_box3.value) - float(a.mr_rate_box3_err.value)))  # 431-2400: '~1e3'
    inv_a_gev = float(a.tc_ew_gev.value) / temp       # 318 GeV if T = T_c
    l_large_volume = math.ceil(8.0 / (g2_ncs * temp))  # 16 sites for L g^2 T >= 8
    gibbs_jumps = float(a.gibbs_terms.value) * float(a.gibbs_sweeps.value)   # 3000
    gibbs_T = (float(a.gibbs_sampler_T.lo), float(a.gibbs_sampler_T.hi))
    gibbs_T_per_jump_implied = gibbs_T[0] / gibbs_jumps                       # 3.3e6: not a D_ov, not a step (record)
    # E23 (3), r20; r21 (2); r22: the sampler from first principles (ESTIMATED HERE), at the single-particle application
    # price. Per jump: n_evol controlled evolutions of length T_OFT = w beta, each beta||H|| x w queries to H; one H query =
    # one application (the chapter's own unit: the TPQ filter's d_beta counts applications); plus an allowance for the
    # Metropolis weight and the jump LCU. The sampler replaces the TPQ preparation, so the comparison is against the TPQ line.
    bH = beta_H_used                                  # 864 (r25 1e2): the same normalization as the filter
    w_oft = math.sqrt(math.log(1 / float(a.gibbs_eps_oft.value)))           # 2.628
    sweeps_hi = float(a.gibbs_sweeps_high.value)
    gibbs_cases = {"low": (int(a.gibbs_n_evol.lo), 1.0, float(a.gibbs_sweeps.value)),
                   "central": (int(a.gibbs_n_evol_central.value), w_oft, float(a.gibbs_sweeps.value)),
                   "high": (int(a.gibbs_n_evol.hi), w_oft, sweeps_hi)}
    gibbs_est = {}
    for name, (n_ev, w, sw) in gibbs_cases.items():
        q_jump = n_ev * w * bH                                                # applications per jump
        t_jump = q_jump * t_dov + float(a.gibbs_accept_T.value)
        n_jumps = float(a.gibbs_terms.value) * sw
        gibbs_est[name] = {"n_evol": n_ev, "window": w, "sweeps": sw, "n_jumps": n_jumps, "dov_per_jump": q_jump,
                           "t_per_jump": t_jump, "t_total": n_jumps * t_jump, "over_tpq": n_jumps * t_jump / t_tpq}

    # shots and wall time (R4: carry exact, round the printed result once); 1 us per T + 0.1 ms per shot
    eps = float(a.eps_obs.value)
    shots_cond = float(a.sigma2_condensate.value) / eps ** 2       # 1e3
    shots_sph = float(a.sigma2_ncs.value) / eps ** 2               # 1e4; r25 record (rate campaign), not current
    shots = float(a.shots_2033.value)                              # 1e4; r25 record, current shots are Result.shots
    wall_shot = t_shot * float(a.t_gate_s.value) + float(a.shot_overhead_s.value)   # 2.008e5 s
    # r23 (H. Lamm 2026-10-02): one machine, serial. r24: each arm at its own shot; a point is one coupling.
    horizon_yr = float(a.campaign_yr.value)                        # 5
    n_pts = float(a.n_coupling.value)                              # 50
    t_gate, t_over = float(a.t_gate_s.value), float(a.shot_overhead_s.value)
    # condensate arm (r24 (1)): the preparation-only shot, 1e3 shots per coupling
    t_shot_cond = C["t_shot"]                                      # 2.2e8 T
    eps_l_cond = faults / t_shot_cond
    wall_shot_cond = t_shot_cond * t_gate + t_over                 # 225 s
    wall_cond_point_yr = shots_cond * wall_shot_cond / SEC_PER_YR  # 2.6 days
    wall_cond_scan_yr = wall_cond_point_yr * n_pts                 # 0.36 yr
    # sphaleron arm: the evolution shot, 1e4 shots per point
    wall_point_yr = shots * wall_shot / SEC_PER_YR                 # 2.6 yr
    wall_serial_yr = wall_point_yr * n_pts                         # 131 yr (the scan)
    shots_in_horizon = horizon_yr * SEC_PER_YR / wall_shot         # sphaleron shots in 5 years
    reduction_needed = {"condensate_point": wall_cond_point_yr / horizon_yr,
                        "condensate_scan": wall_cond_scan_yr / horizon_yr,
                        "sphaleron_point": wall_point_yr / horizon_yr,
                        "sphaleron_point_one_year": wall_point_yr / 1.0,
                        "scan": wall_serial_yr / horizon_yr}
    # the sphaleron shot that puts the scan inside the horizon at the same shot count (the ~0.1 ms overhead is kept)
    t_shot_to_fit_scan = (horizon_yr * SEC_PER_YR / (shots * n_pts) - t_over) / t_gate
    # r24 levers for one sphaleron point under a year (none applied to the box):
    #   low_step: the low end of the step range (fastest-mode phase), priced through the chain at its own R-TOL
    #   eps15:    a 15% statistical target, shots x (eps_obs / 0.15)^2
    #   hold:     RECORD (r23): the link-select lines held through the copy (78 Toffolis per query, ~75 more qubits)
    #   spectral amplification (the 2Q t floor): NOT PRICED, research
    eps_lever = float(a.lever_eps_obs.value)                       # 0.15
    shots_sph_lever = float(a.sigma2_ncs.value) / eps_lever ** 2   # 4444
    hold_saved = spq["toffoli_items"]["fixup_iteration"]           # 78
    hold_qubits = n_links - spq["anc_iteration"]                   # 75
    P_hold = _price_sp(n_trot, toff_saved_=hold_saved)
    horizon_levers = {"low_step": t_shot / P_low["t_shot"], "eps15": shots / shots_sph_lever,
                      "hold": t_shot / P_hold["t_shot"]}
    horizon_levers["low_step_and_eps15"] = horizon_levers["low_step"] * horizon_levers["eps15"]
    wall_shot_low = P_low["t_shot"] * t_gate + t_over
    wall_yr_levers = {"low_step": shots * wall_shot_low / SEC_PER_YR,
                      "eps15": shots_sph_lever * wall_shot / SEC_PER_YR,
                      "low_step_and_eps15": shots_sph_lever * wall_shot_low / SEC_PER_YR}

    # ---- r25: T-depth, factories and the two tiers (rulings R1, R4, R9, R10; Ch. 9 rulings (a), (e), (f)) ----------
    # Free ancilla outside a query (the term register and the link-read ancilla): the Higgs U_x run 4 at a time.
    umult_concurrent = (sum(anc_query.values()) + sum(anc_select.values())) // umult_ancilla        # 19 // 4 = 4

    def _depth(P_, variant_=depth_variant):
        """(lo, hi) T-depth of one shot of the circuit P_ (applications are sequential)."""
        lo_, hi_ = application_depth_2O(P_["n_q_per_dov"], P_["t_rot"], variant_, lift, outer, V, higgs_sel_toff,
                                        umult_toffolis, umult_concurrent)
        return (P_["n_dov"] * lo_, P_["n_dov"] * hi_)
    t_short = float(a.t_short_fm.value) / float(a.a_fm.value)                 # 5a
    sd_short = trotter_state_dependent_steps(L, n_adj, t_short, eps_tr, float(a.temp_lat_2033.value), d=D)
    P5 = _price_sp(sd_short["n_central"], t_=t_short)                        # the 5a shot at its own depth (R4)
    runs = {"sph10": P, "sph5": P5, "cond": C}
    depth = {k: _depth(v) for k, v in runs.items()}
    depth_register = {k: _depth(v, "register") for k, v in runs.items()}    # record: as itemized, no fan-out
    f_star_runs = {k: (runs[k]["t_shot"] / depth[k][1], runs[k]["t_shot"] / depth[k][0]) for k in runs}
    f_star_register = {k: (runs[k]["t_shot"] / depth_register[k][1], runs[k]["t_shot"] / depth_register[k][0])
                       for k in runs}
    wall_shot_run = {k: v["t_shot"] * t_gate + t_over for k, v in runs.items()}

    # shots. Author ruling 2026-10-05: the 2033 real-time measurement is a METHOD TEST on 3^3, the early growth of the
    # full <dN_CS^2> (background included) at 5a and 10a, one symmetric-side temperature, 30% each, compared with the
    # free-field quantum value and classical real-time simulation; no rate. First result = that test. Campaign = the
    # test + the condensate at 4-8 temperatures at 10% (a temperature scan of the growth is not resolvable: the free
    # signal moves 3-5% over T a = 0.4-0.6). The 5a relative variance is the 10a one times (S10/S5)^2 = 4.95 from the
    # derived free-field signals at fixed absolute readout noise (r25 record: 4, under linear growth).
    # r25 (R1, e) record: campaign 4-8 temperatures at 10% at 10a, the 5a check at one, the condensate at each.
    eps_first = float(a.eps_first.value)
    var_ratio5_r25 = float(a.ncs_var_ratio_short.value)                                          # 4 (record)
    var_ratio5 = (ncs_bg[10.0] / ncs_bg[5.0]) ** 2                                               # 4.948
    s2 = float(a.sigma2_ncs.value)
    first_shots = {"sph10": s2 / eps_first ** 2, "sph5": var_ratio5 * s2 / eps_first ** 2}       # 1111, 5498
    shots_first_hadamard = {tt: ncs_relvar[tt] / eps_first ** 2 for tt in ncs_relvar}            # 3.7e14, 1.2e15
    n_temp = (int(a.n_temperatures.lo), int(a.n_temperatures.hi))
    camp_shots = tuple(dict(first_shots, cond=nt * float(a.sigma2_condensate.value) / eps ** 2) for nt in n_temp)
    camp_shots_r25 = tuple({"sph10": nt * s2 / eps ** 2, "sph5": var_ratio5_r25 * s2 / eps ** 2,
                            "cond": nt * float(a.sigma2_condensate.value) / eps ** 2} for nt in n_temp)
    wall_first = sum(n * wall_shot_run[k] for k, n in first_shots.items())
    wall_first_parts = {k: n * wall_shot_run[k] for k, n in first_shots.items()}
    wall_camp = tuple(sum(n * wall_shot_run[k] for k, n in cs.items()) for cs in camp_shots)
    wall_camp_parts = tuple({k: n * wall_shot_run[k] for k, n in cs.items()} for cs in camp_shots)
    # R9 exports: the headline shot (10a) for t_per_shot / depth / F*; floor and factories summed over the campaign runs
    ex_head = depth_exports(P["t_shot"], depth["sph10"], camp_shots[0]["sph10"], a.t_gate_s, a.shot_overhead_s)
    ex_runs = tuple({k: depth_exports(runs[k]["t_shot"], depth[k], n, a.t_gate_s, a.shot_overhead_s)
                     for k, n in cs.items()} for cs in camp_shots)
    floor_camp = (sum(e["floor_wall_s"][0] for e in ex_runs[0].values()),
                  sum(e["floor_wall_s"][1] for e in ex_runs[1].values()))
    fact_camp = (sum(e["factories_for_1yr"][0] for e in ex_runs[0].values()),
                 sum(e["factories_for_1yr"][1] for e in ex_runs[1].values()))
    baseline_ok = all(f[0] >= FACTORY_BASELINE for f in f_star_runs.values())
    # the campaign against the 5-year horizon on one machine: what fits, and the cut it needs
    camp_yr = (wall_camp[0] / SEC_PER_YR, wall_camp[1] / SEC_PER_YR)
    camp_reduction = (camp_yr[0] / horizon_yr, camp_yr[1] / horizon_yr)
    first_yr = wall_first / SEC_PER_YR
    # r25 record: the rate campaign (4-8 temperatures at 10%) and its 30% variant
    wall_camp_r25_yr = tuple(sum(n * wall_shot_run[k] for k, n in cs.items()) / SEC_PER_YR for cs in camp_shots_r25)
    # the method test at 10% instead of 30% (a lever, not applied)
    wall_first_at_eps_yr = wall_first * (eps_first / eps) ** 2 / SEC_PER_YR
    tiers = {
        "first_result": {"temperatures": 1, "eps": eps_first, "shots": first_shots, "wall_s": wall_first,
                         "wall_s_parts": wall_first_parts, "wall_yr": first_yr},
        "campaign": {"temperatures": n_temp, "eps": eps, "eps_growth": eps_first, "shots": camp_shots,
                     "wall_s": wall_camp, "wall_s_parts": wall_camp_parts, "wall_yr": camp_yr,
                     "over_horizon": camp_reduction},
    }

    src_sp = ("DERIVED HERE r22 (single_particle_query_2O; construction of arxiv_2607_28524; unary iteration "
              "Babbush_PRX_2018; 2O encoding arxiv_2312_10285)")
    src_lift = "arxiv_2607_28524 (Supplement theorem); counts DERIVED HERE r22 (lift_cost)"
    rot_note = (f"R-TOL eps_rot = sqrt({eps_syn:g}/{n_rot_shot:.6g}) = {eps_rot:.4g}, RUS fit {t_rot:.3f} T")
    breakdown = (
        Primitive("sp_query_toffoli", spq["toffoli"] * n_queries, T_PER_TOFFOLI[tc], CircuitStatus.COMPILED, src_sp,
                  f"{spq['toffoli']} Toffolis per query of the single-particle H_W ({spq['link_read_toffoli']} of them the "
                  f"link read), x {n_dov} applications x {n_q_per_dov} queries; {T_PER_TOFFOLI[tc]} T per Toffoli (R5)"),
        Primitive("sp_query_direct_T", spq["t_direct"] * n_queries, 1.0, CircuitStatus.COMPILED, src_sp,
                  f"M and M^dag at {mux_ix['t_direct']} T each, controlled-H of PREPARE and its inverse"),
        Primitive("sp_query_rotations", spq["n_rot"] * n_queries, t_rot, CircuitStatus.SCALING, SYNTHESIS_SRC,
                  f"{spq['n_rot']} per query (PREPARE pair 4, controlled QSP phase 2); " + rot_note),
        Primitive("lift_toffoli", lift["toffoli"] * n_dov, T_PER_TOFFOLI[tc], CircuitStatus.SCALING, src_lift,
                  f"U_L and U_R: {lift['leaves']} selected Majorana strings each, {lift['leaves'] - 1} Toffolis each"),
        Primitive("lift_direct_T", lift["t_direct"] * n_dov, 1.0, CircuitStatus.SCALING, src_lift,
                  "controlled-H in the uniform-over-3 preparation of the coordinates, in and out"),
        Primitive("lift_rotations", lift["n_rot"] * n_dov, t_rot, CircuitStatus.SCALING, SYNTHESIS_SRC,
                  "uniform-over-3 preparation of the D coordinates, in and out; " + rot_note),
        Primitive("outer_reflection_toffoli", outer["toffoli"] * n_dov, T_PER_TOFFOLI[tc], CircuitStatus.SCALING,
                  "DERIVED HERE r22", f"all-zero flag on the {n_be_anc} block-encoding ancilla for the QSP that evolves or "
                                      "filters with the application"),
        Primitive("outer_reflection_rotation", outer["n_rot"] * n_dov, t_rot, CircuitStatus.SCALING, SYNTHESIS_SRC,
                  "the QSP phase of the evolution / filter polynomial; " + rot_note),
        Primitive("higgs_2O_group_action", V * n_um * n_dov, t_umult, CircuitStatus.SCALING, "arxiv_2312_10285",
                  f"{V} Higgs actions per application at {n_um} U_x = {t_umult:.0f} T ({umult_toffolis:.0f} Toffoli at 7 T, "
                  f"{umult_ancilla} clean ancilla; tab:tgatecost, COMPILED; no controlled variant, so SCALING); outside sgn"),
        Primitive("higgs_select_toffoli", higgs_sel_toff * n_dov, T_PER_TOFFOLI[tc], CircuitStatus.SCALING,
                  "app06 2033 derivation; Babbush_PRX_2018",
                  f"{lcu['higgs']} Higgs strings per site (r21 c) at one SELECT Toffoli each, once per application"),
    )

    inter = {
        "sgn_architecture": arch,
        "lq_fermion": lq_f, "lq_gauge": lq_g, "lq_higgs": lq_h, "lq_ancilla": lq_anc, "lq_total": lq,
        "link_qubits": g.link_qubits, "higgs_qubits_per_site": n_rho_q + g.link_qubits,
        "qubits_per_site": qubits_per_site,
        "lq_dwf_backup_L5_8": lq_dwf_backup[0], "lq_dwf_backup_L5_16": lq_dwf_backup[1],
        "d_sgn_exact": d_sgn_exact, "d_sgn_log2_would_be": d_sgn_log2, "d_sgn": d_sgn,
        "d_sgn_chapter": d_sgn_chapter, "d_sgn_stated": d_sgn_chapter,
        "n_links": n_links,
        "t_per_toffoli": T_PER_TOFFOLI[tc],
        # r22: the single-particle kernel and its block encoding
        "q_modes": q_modes, "index_register": q_index, "index_register_qubits": n_index,
        "n_lcu_unitaries": spq["n_lcu_unitaries"], "term_register_qubits": spq["term_register"],
        "m0_wilson": m0, "alpha_hw": alpha_hw, "norm_hw_free": norm_hw_free, "be_normalization": be_norm,
        "d_sgn_if_kappa_is_norm_over_gap": d_sgn_norm_convention,
        "t_per_shot_if_kappa_is_norm_over_gap": P_kappa_norm["t_shot"],
        "sp_query": spq, "sp_query_toffoli": spq["toffoli"], "sp_query_toffoli_stated": float(a.sp_query_toffoli_stated.value),
        "sp_query_t_direct": spq["t_direct"], "sp_query_n_rot": spq["n_rot"],
        "sp_link_read_toffoli": spq["link_read_toffoli"], "sp_link_read_share_of_query": link_read_share,
        "sp_rotation_share_of_query": rot_share,
        "mux_index_per_bit": mux_ix["per_bit"], "mux_index_t_direct": mux_ix["t_direct"],
        "mux_index_toffoli": mux_ix["toffoli"], "mux_index_ancilla": mux_ix["ancilla"],
        "c_sp": c_sp, "c_sp_stated": float(a.c_sp_stated.value),
        "n_queries_per_dov": n_q_per_dov, "t_sgn_per_dov": P["t_sgn"],
        "lift": lift, "c_lift": P["c_lift"], "outer_reflection": outer, "c_outer": P["c_outer"],
        "n_be_ancilla": n_be_anc,
        "higgs_strings_per_site": lcu["higgs"], "higgs_strings_per_site_stated": float(a.higgs_strings_per_site_stated.value),
        "higgs_select_toffolis": higgs_sel_toff, "c_higgs_umult_T": c_higgs, "c_higgs_per_dov": c_higgs_app,
        "n_rot_per_dov": P["n_rot_app"],
        "sgn_share_of_dov": P["t_sgn"] / t_dov, "lift_share_of_dov": P["c_lift"] / t_dov,
        "higgs_share_of_dov": c_higgs_app / t_dov,
        "t_dov_shift_coherent": t_dov_shift_coh,
        "t_per_rot_coherent": t_rot_coh, "eps_rot_coherent_same_total": eps_coh,
        "clifford_2O": a.clifford_2O.value,
        # R-TOL: the circuit's own tolerance
        "eps_syn": eps_syn, "n_rot_shot": n_rot_shot, "eps_rot": eps_rot, "t_per_rot": t_rot,
        "t_per_rot_retired": t_rot_retired,
        "t_umult": t_umult, "t_umult_report_convention": t_umult_report,
        "umult_toffolis": umult_toffolis, "umult_clean_ancilla": umult_ancilla,
        "umult_ancilla_fit_in_box": umult_ancilla <= lq_anc,
        "n_umult_per_higgs_site": n_um,
        # ancilla itemization (peak of the schedule)
        "ancilla_stack": anc_stack, "ancilla_application": anc_app, "ancilla_query": anc_query,
        "ancilla_select": anc_select, "ancilla_not_concurrent": anc_not_concurrent,
        "lq_ancilla_asserted": lq_anc_asserted,
        # the Trotter step count: the state-dependent estimate (headline), its range, and the worst-case bound
        "t_over_a_2033": t_over_a, "eps_trotter": eps_tr,
        "trotter_rule": a.trotter_rule_2033.value, "trotter_state_dependent": sd,
        "temp_lat": float(a.temp_lat_2033.value),
        "trotter_lambda_sd": sd["lambda_var"], "trotter_lambda_sd_stated": float(a.lambda_sd_stated.value),
        "trotter_n_modes": sd["n_modes"], "trotter_n_modes_stated": float(a.n_modes_stated.value),
        "n_trotter_range": n_trot_range,
        "n_trotter_range_stated": (float(a.n_trotter_range_stated.lo), float(a.n_trotter_range_stated.hi)),
        "trotter_bound": tb, "g2_bare": g2,
        "trotter_lambda": tb["lambda"], "trotter_lambda_stated": float(a.lambda_2033_stated.value),
        "trotter_lambda_min": tb_min["lambda"], "trotter_g2_at_min": g2_at_min, "n_trotter_at_lambda_min": tb_min["n_steps"],
        "n_trotter_bound": tb["n_steps"], "n_trotter_bound_stated": float(a.n_trotter_bound_stated.value),
        "n_trotter_2033": n_trot, "n_trotter_2033_stated": float(a.n_trotter_2033_stated.value),
        "n_trotter_by_rule": n_trot_by_rule, "n_trotter_pre_r21": n_trot_pre,
        "lambda_needed_for_pre_r21_steps": lambda_for_pre,
        "trotter_dt_lattice": t_over_a / n_trot, "trotter_dt_fm": float(a.t_evol_fm.value) / n_trot,
        "trotter_dt_lattice_bound": tb["dt_lat"],
        # applications per step: the Jacobi-Anger order at x = 2Q dt
        "eps_qsp": eps_qsp, "jacobi_anger_x": P["x"], "n_dov_per_step": per_step,
        "n_dov_per_step_stated": float(a.n_dov_per_step_stated.value),
        "n_dov_per_step_by_rule": {k: v["per_step"] for k, v in by_rule.items()},
        "n_dov_per_step_range": (P_low["per_step"], P_high["per_step"]),
        "n_dov_per_step_bound_stated": float(a.n_dov_per_step_bound_stated.value),
        "n_dov_per_step_worst_stated": float(a.n_dov_per_step_worst_stated.value),
        "n_dov_by_rule": {k: v["n_dov"] for k, v in by_rule.items()},
        "n_dov_range": (P_low["n_dov"], P_high["n_dov"]),
        "n_dov_floor": n_dov_floor, "n_dov_floor_stated": float(a.n_dov_floor_stated.value),
        "t_per_shot_floor": t_shot_floor, "floor_ratio_to_envelope": t_shot_floor / env,
        "t_per_shot_by_trotter_rule": {k: v["t_shot"] for k, v in by_rule.items()},
        "t_per_shot_range": (P_low["t_shot"], P_high["t_shot"]),
        "t_per_shot_range_stated": (float(a.t_shot_range_stated.lo), float(a.t_shot_range_stated.hi)),
        "t_per_shot_bound": by_rule["bound"]["t_shot"], "t_per_shot_bound_stated": float(a.t_shot_bound_stated.value),
        "t_per_shot_empirical_step": by_rule["step"]["t_shot"],
        "t_per_shot_empirical_step_stated": float(a.t_shot_empirical_step_stated.value),
        "t_per_shot_at_pre_r21_steps": by_rule["stated"]["t_shot"],
        "t_per_shot_if_not_self_inverse": P_not_self_inverse["t_shot"],
        # not priced (records)
        "unpriced_gauge_T_per_link_step": t_gauge_link, "unpriced_gauge_T_per_step": t_gauge_step,
        "unpriced_gauge_share_of_step": gauge_step_share,
        "unpriced_fpaa_reflections": n_reflections, "unpriced_fpaa_reflection_controls": n_reflect_controls,
        "unpriced_fpaa_reflection_T": t_fpaa_reflections,
        # the shot
        "t_per_dov": t_dov, "t_per_dov_stated": float(a.t_per_dov_stated.value),
        "n_dov_evolution": n_dov_evol, "n_dov_evolution_stated": float(a.n_dov_evolution_stated.value),
        "t_evolution": t_evol, "t_evolution_stated": float(a.t_evolution_stated.value),
        "delta_beta": float(a.delta_beta.value),
        "d_beta_exact": d_beta_exact, "d_beta": d_beta, "d_beta_stated": float(a.d_beta_stated.value),
        "beta_H_used": beta_H_used, "d_beta_record_r25": d_beta_record_r25, "n_dov_tpq_r25": n_dov_tpq_r25,
        "r25_headline": {"n_dov": P_r25["n_dov"], "t_per_shot": P_r25["t_shot"], "t_per_shot_condensate": C_r25["t_shot"]},
        "d_beta_at_delta_beta_1e-3": d_beta_at_1e3,
        "n_dov_tpq": n_dov_tpq, "n_filter_calls": n_filter_calls, "n_dov_tpq_one_call_per_round": n_dov_tpq_one_call_per_round, "n_dov_tpq_stated": float(a.n_dov_tpq_stated.value),
        "t_tpq": t_tpq, "t_tpq_stated": float(a.t_tpq_stated.value),
        "tpq_dominates": t_tpq > t_evol, "tpq_share_of_shot": t_tpq / t_shot, "evolution_share_of_shot": t_evol / t_shot,
        "n_dov_total": n_dov, "n_queries_total": n_queries,
        "t_per_shot": t_shot,
        "t_per_shot_stated": float(a.t_per_shot_2033_stated.value),
        "ratio_to_envelope": ratio_env, "ratio_to_envelope_stated": float(a.ratio_representative_stated.value),
        "query_count_cut_needed": ratio_env,
        "inside_2033_envelope": t_shot <= env,
        "t_per_shot_kappa_physical": t_phys,
        "t_per_shot_kappa_physical_stated": (float(a.t_physical_stated.lo), float(a.t_physical_stated.hi)),
        "faults_per_shot": faults,
        "eps_l_implied": eps_l_implied, "eps_l_stated": float(a.eps_l_2033_stated.value),
        "eps_l_from_envelope": faults / env,
        # RECORD: the retired second-quantized construction (r17-r21) and everything it printed
        "second_quantized_record": second_quantized_record,
        "link_frame": frame,
        "lcu_per_site": lcu, "lcu_terms_per_site": lcu["total"],
        "lcu_terms_per_site_stated": float(a.lcu_terms_per_site_stated.value),
        "lcu_terms": lcu_terms, "lcu_terms_stated": float(a.lcu_terms_stated.value), "lcu_toffolis": lcu_toff,
        "lcu_terms_pre_r21": lcu_terms_pre, "lcu_terms_per_site_implied_pre_r21": lcu_terms_pre / V,
        "lcu_toffolis_per_term": lcu_toff / lcu_terms, "lcu_address_qubits": n_addr,
        "mux_per_bit": mux["per_bit"], "mux_toffoli": mux["toffoli"], "mux_t_direct": mux["t_direct"],
        "mux_ancilla": mux["ancilla"], "mux_t_report": toffoli_t(mux["toffoli"], tc) + mux["t_direct"],
        "mux_link_toffoli": ml["toffoli"], "mux_link_t_direct": ml["t_direct"], "mux_link_n_mux": ml["n_mux"],
        "mux_link_t": t_link_mux,
        "hop_items": tuple((x["name"], x["toffoli"], x["t_direct"], x["n_rot"]) for x in hop["items"]),
        "hop_toffoli_per_link": hop["toffoli"], "hop_t_direct_per_link": hop["t_direct"],
        "hop_n_rot_per_link": n_rot_link_draft, "hop_clifford_t_per_link": t_link_clifford_t,
        "hop_t_per_link": SQ_draft["t_link"], "hop_n_spin": hop["n_spin"], "hop_undo": hop["undo"], "hop_mcx": hop["mcx"],
        "hop_share": a.hop_share.value,
        "draft_over_mux_per_link": SQ_draft["t_link"] / t_link_mux,
        "draft_over_mux_per_link_stated": float(a.draft_over_mux_link_stated.value),
        "sq_c_be_draft_frame": SQ_draft["c_be"], "sq_c_be_draft_frame_stated": float(a.c_be_draft_stated.value),
        "sq_t_per_shot_draft_frame": SQ_draft["t_shot"], "sq_t_per_shot_draft_frame_stated": float(a.t_shot_draft_stated.value),
        "sq_draft_frame_eps_rot": SQ_draft["eps_rot"], "sq_draft_frame_t_per_rot": SQ_draft["t_rot"],
        "sq_draft_frame_n_rot_shot": SQ_draft["n_rot_shot"],
        "draft_over_mux_shot": SQ_draft["t_shot"] / SQ_mux["t_shot"],
        "draft_over_mux_shot_stated": float(a.draft_over_mux_stated.value),
        "sq_c_be": SQ["c_be"], "sq_c_be_stated": float(a.c_be_stated.value),
        "sq_c_be_toffoli_T": c_toff, "sq_c_be_hop_T": SQ["c_hop"], "sq_c_be_higgs_T": c_higgs, "sq_c_be_rot_T": SQ["c_rot_be"],
        "sq_n_rot_per_query": SQ["n_rot_query"], "sq_n_rot_shot": SQ["n_rot_shot"], "sq_t_per_rot": SQ["t_rot"],
        "sq_t_per_dov": SQ["t_dov"], "sq_t_per_dov_stated": float(a.t_per_dov_sq_stated.value),
        "sq_n_dov": SQ["n_dov"], "sq_n_queries": SQ["n_queries"],
        "sq_t_per_shot": SQ["t_shot"], "sq_t_per_shot_stated": float(a.t_per_shot_sq_stated.value),
        "sq_ancilla_stack": sq_anc_stack, "sq_ancilla_be": sq_anc_be, "sq_ancilla_term": sq_anc_term_mux,
        "sq_ancilla_link_draft_frame": anc_link_draft, "sq_ancilla_link_transients": transient,
        "sq_lq_ancilla": sq_lq_anc, "sq_lq_total": sq_lq,
        "sq_lq_ancilla_draft_frame": sq_lq_anc_draft, "sq_lq_total_draft_frame": sq_lq_draft,
        "n_queries_pre_r21": n_queries_pre,
        "c_be_r20_centre_draft_frame": R_draft["c_be"], "c_be_r20_floor_multiplexer": R_mux["c_be"],
        "c_be_r20_top_draft_w2": R_w2["c_be"], "c_be_floor_zero_record": R_zero["c_be"],
        "w2_t_per_link": w2_t_link, "w2_rot_per_link": w2_rot_link, "w2_eps_rot": R_w2["eps_rot"],
        "w2_t_per_rot": R_w2["t_rot"],
        "c_be_rotation_free": c_be_rotation_free,
        "ratio_centre_to_rotation_free": R_draft["c_be"] / c_be_rotation_free,   # record (r19 prose 8.3; retired r20)
        "c_be_records": c_be_records, "t_per_shot_records": t_shot_records,
        "t_group_action_retired": retired.t_const,
        "c_be_before_round_b": c_be_before_round_b,
        "t_per_shot_before_round_b": t_shot_before_round_b,
        # thermal-route alternatives (app06 route (iv))
        "haar_entropy_deficit_nats": dS_haar, "haar_rounds": haar_rounds,
        "haar_rounds_stated": float(a.haar_rounds_stated.value),
        "gibbs_jumps": gibbs_jumps, "gibbs_sampler_T_stated": gibbs_T,
        "gibbs_T_per_jump_implied": gibbs_T_per_jump_implied,
        # E23 (3), r20; r21 (2); r22: the sampler estimated here at the application price (low / central / high)
        "gibbs_estimate": gibbs_est, "gibbs_window_factor": w_oft,
        "gibbs_T_estimate": (gibbs_est["low"]["t_total"], gibbs_est["central"]["t_total"], gibbs_est["high"]["t_total"]),
        "gibbs_T_stated": (float(a.gibbs_T_stated.lo), float(a.gibbs_T_stated.hi)),
        "gibbs_T_central_stated": float(a.gibbs_T_central_stated.value),
        "gibbs_over_tpq": gibbs_est["central"]["over_tpq"],
        "gibbs_over_tpq_stated": float(a.gibbs_over_tpq_stated.value),
        "thermal_route": "TPQ",
        # chapter-open pass (referee C1 / C3, 2026-10-05), DERIVED HERE
        "typical_ref": typ, "typical_ref_ln_p": typ["ln_p"], "typical_ref_aa_rounds": typ["aa_rounds"],
        "typical_ref_ln_norm_offset": typ["ln_norm_offset"], "typical_ref_n_frozen": typ["n_frozen"],
        "heralded_accept_min": 1.0 / n_filter_calls,
        "ncs_free_background": ncs_bg, "ncs_lambda": ncs_lam, "ncs_hadamard_rel_var": ncs_relvar,
        "ncs_free_background_helicity": ncs_bg_helicity, "ncs_hadamard_rel_var_helicity": ncs_relvar_helicity,
        "ncs_free_background_fwd": ncs_bg_fwd, "ncs_hadamard_rel_var_fwd10": ncs_relvar_fwd10,
        "mr_small_box_suppression_range": mr_suppression_range,
        "ncs_hadamard_rel_var_unit_signal": ncs_relvar_unit_signal,
        "beta_lattice_classical": beta_lattice, "box_L_g2T": box_lg2t, "box_L_g2T_ew": L * temp * 0.42,
        "mr_small_box_suppression": mr_suppression, "inv_a_gev_at_tc": inv_a_gev, "L_large_volume": l_large_volume,
        # shots / wall time. "shots_sphaleron_per_coupling" and "shots" are r25 records of the rate campaign (1e4);
        # the current 2033 shots are in tiers / Result.shots
        "shots_condensate_per_coupling": shots_cond, "shots_sphaleron_per_coupling": shots_sph,
        "shots": shots,
        "wall_s_per_shot": wall_shot, "wall_s_per_shot_stated": float(a.wall_per_shot_2033_stated.value),
        "wall_yr_per_point": wall_point_yr, "wall_yr_per_point_stated": float(a.wall_yr_per_point_stated.value),
        "wall_yr_serial_scan": wall_serial_yr,
        "wall_yr_serial_scan_stated": float(a.wall_yr_serial_stated.value),
        # r24 (1): the condensate arm, priced at the preparation alone (static observable, no real-time evolution)
        "condensate_shot": C, "t_per_shot_condensate": t_shot_cond,
        "t_per_shot_condensate_stated": float(a.t_shot_condensate_stated.value),
        "n_dov_condensate": C["n_dov"], "n_rot_shot_condensate": C["n_rot_shot"], "eps_rot_condensate": C["eps_rot"],
        "t_per_rot_condensate": C["t_rot"],
        "eps_l_condensate": eps_l_cond, "eps_l_condensate_stated": float(a.eps_l_condensate_stated.value),
        "wall_s_per_shot_condensate": wall_shot_cond,
        "wall_day_condensate_point": wall_cond_point_yr * SEC_PER_YR / 86400.0,
        "wall_day_condensate_point_stated": float(a.wall_day_condensate_point_stated.value),
        "condensate_inside_2033_envelope": t_shot_cond <= env,
        # r23/r24: one machine. Per-arm serial wall times, what 5 years hold, and the cost reduction that would fit
        "campaign_horizon_yr": horizon_yr,
        "wall_yr_condensate_point": wall_cond_point_yr,
        "wall_yr_sphaleron_point": wall_point_yr,
        "wall_yr_condensate_scan": wall_cond_scan_yr,
        "wall_yr_condensate_scan_stated": float(a.wall_yr_condensate_scan_stated.value),
        "fits_horizon": {"condensate_point": wall_cond_point_yr <= horizon_yr,
                         "condensate_scan": wall_cond_scan_yr <= horizon_yr,
                         "sphaleron_point": wall_point_yr <= horizon_yr, "scan": wall_serial_yr <= horizon_yr},
        "sphaleron_point_under_one_year": wall_point_yr <= 1.0,
        "shots_in_horizon": shots_in_horizon,
        "reduction_needed": reduction_needed,
        "reduction_scan_stated": float(a.reduction_scan_stated.value),
        "t_per_shot_to_fit_scan": t_shot_to_fit_scan,
        "t_per_shot_to_fit_scan_stated": float(a.t_shot_to_fit_scan_stated.value),
        "horizon_levers": horizon_levers,
        "horizon_levers_stated": {"low_step": float(a.lever_low_step_stated.value),
                                  "eps15": float(a.lever_eps15_stated.value)},
        "wall_yr_levers": wall_yr_levers,
        "wall_yr_levers_stated": {"low_step": float(a.wall_yr_lever_low_step_stated.value),
                                  "eps15": float(a.wall_yr_lever_eps15_stated.value),
                                  "low_step_and_eps15": float(a.wall_yr_lever_both_stated.value)},
        "lever_eps_obs": eps_lever, "lever_shots_sphaleron": shots_sph_lever,
        "low_step_over_floor": P_low["t_shot"] / t_shot_floor,     # 'gain at most 1.4 more' below the low end
        "lever_hold_toffoli_saved": hold_saved, "lever_hold_qubits": hold_qubits,
        "lever_hold_t_per_shot": P_hold["t_shot"],
        # r24 records: the r22-r23 headline (kappa 1e2, t = 50a) and the kappa ~ 1e2 sensitivity at t = 10a
        "kappa_2033": kappa, "kappa_r22": kappa_r22, "d_sgn_r22": d_sgn_r22, "t_over_a_r22": t_over_a_r22,
        "r22_headline": {"n_trot": P_r22["n_trot"], "per_step": P_r22["per_step"], "n_dov": P_r22["n_dov"],
                         "t_per_dov": P_r22["t_dov"], "t_per_shot": P_r22["t_shot"], "c_sp": P_r22["c_sp"],
                         "n_rot_shot": P_r22["n_rot_shot"], "t_per_rot": P_r22["t_rot"],
                         "t_tpq": n_dov_tpq_r25 * P_r22["t_dov"],
                         "wall_yr_per_point": shots * (P_r22["t_shot"] * t_gate + t_over) / SEC_PER_YR},
        "r22_bound": {"n_trot": tb_r22["n_steps"], "per_step": P_r22_bound["per_step"], "t_per_shot": P_r22_bound["t_shot"]},
        "trotter_state_dependent_r22": sd_r22,
        "t_per_shot_kappa_r22": P_kappa_r22["t_shot"],
        "t_per_shot_kappa_r22_stated": float(a.t_shot_kappa_r22_stated.value),
        # r25: Yukawa inside (a), the tiers (R1, (e)), each shot at its own depth (R4), T-depth and factories (R9, R10)
        "yukawa_placement": a.yukawa_placement.value, "n_sgn_per_application": n_sgn,
        "r24_t_per_shot": P_r24["t_shot"], "r24_t_per_dov": P_r24["t_dov"], "r24_t_per_shot_condensate": C_r24["t_shot"],
        "yukawa_over_r24_per_application": t_dov / P_r24["t_dov"], "yukawa_over_r24_per_shot": t_shot / P_r24["t_shot"],
        "lq_ancilla_itemized": lq_anc_itemized, "lq_ancilla_fanout": anc_fanout, "depth_variant": depth_variant,
        "umult_concurrent": umult_concurrent,
        "short_time_lat": t_short, "trotter_state_dependent_short": sd_short,
        "short_shot": {k: P5[k] for k in ("n_trot", "per_step", "x", "n_dov", "n_rot_shot", "eps_rot", "t_rot", "t_dov",
                                          "t_shot")},
        "t_per_shot_short": P5["t_shot"], "eps_l_short": faults / P5["t_shot"],
        "wall_s_per_shot_short": wall_shot_run["sph5"],
        "t_depth_by_run": depth, "t_depth_by_run_register": depth_register,
        "f_star_by_run": f_star_runs, "f_star_by_run_register": f_star_register,
        "wall_s_per_shot_by_run": wall_shot_run, "baseline_ok_all_runs": baseline_ok,
        "tiers": tiers,
        "exports_by_run": ex_runs,
        "t_per_shot_exports": ex_head["t_per_shot"],
        "t_depth_per_shot": ex_head["t_depth_per_shot"],
        "f_star": ex_head["f_star"],
        "floor_wall_s": floor_camp,
        "factories_for_1yr": fact_camp,
        "wall_first_result_s": (wall_first, wall_first),
        "wall_campaign_s": wall_camp,
        "method_test_var_ratio_5a": var_ratio5, "method_test_var_ratio_5a_r25": var_ratio5_r25,
        "shots_first_hadamard": shots_first_hadamard,
        "ncs_free_background_classical": ncs_bg_classical, "ncs_free_background_zero_point_10a": ncs_bg_zero_point10,
        "ncs_free_background_temp_span_10a": ncs_bg_temp_span10,
        "wall_yr_campaign_r25_rate": wall_camp_r25_yr, "wall_yr_method_test_at_campaign_eps": wall_first_at_eps_yr,
        "wall_day_condensate_temperature": shots_cond * wall_shot_run["cond"] / 86400.0,
        "slowdown_register": (FACTORY_BASELINE / f_star_register["sph10"][1], FACTORY_BASELINE / f_star_register["sph10"][0]),
        "t_per_shot_short_stated": float(a.t_shot_short_stated.value), "eps_l_short_stated": float(a.eps_l_short_stated.value),
        "wall_s_short_stated": float(a.wall_s_short_stated.value),
        "wall_s_condensate_stated": float(a.wall_s_condensate_stated.value),
        "f_star_register_stated": (float(a.f_star_register_stated.lo), float(a.f_star_register_stated.hi)),
        "slowdown_register_stated": (float(a.slowdown_register_stated.lo), float(a.slowdown_register_stated.hi)),
        "f_star_stated": float(a.f_star_2033_stated.value),
        "shots_first_stated": (float(a.shots_first_stated.lo), float(a.shots_first_stated.hi)),
        "wall_yr_first_stated": float(a.wall_yr_first_stated.value),
        "wall_yr_campaign_stated": (float(a.wall_yr_campaign_stated.lo), float(a.wall_yr_campaign_stated.hi)),
        "campaign_over_horizon_stated": (float(a.campaign_over_horizon_stated.lo), float(a.campaign_over_horizon_stated.hi)),
        "wall_yr_campaign_r25_rate_stated": (float(a.wall_yr_campaign_r25_stated.lo), float(a.wall_yr_campaign_r25_stated.hi)),
    }
    notes = (
        f"Architecture (E26 (i)): {arch}. sgn acts on the single-particle kernel H_W[U] (Q = {q_modes}) on a {n_index}-qubit "
        f"index register and the link registers; two selected-Majorana unitaries lift it to the Fock register "
        "(arxiv_2607_28524; compiled here). The second-quantized construction of r17-r21 is a record, not a valid circuit.",
        f"d_sgn = ceil(kappa ln(1/delta)) = {d_sgn} (natural log; log2 gives {d_sgn_log2:.0f}); chapter rounds to 700. "
        f"kappa = alpha/Delta, alpha = 2D + |D - m_0| = {alpha_hw:g} at m_0 = {m0:g}; read as ||H_W||/Delta the degree "
        f"would be {d_sgn_norm_convention} and the shot {P_kappa_norm['t_shot']:.3g} T.",
        f"Query: {spq['toffoli']} Toffolis ({spq['link_read_toffoli']} the link read) + {spq['t_direct']} T + {spq['n_rot']} "
        f"rotations = c_sp {c_sp:.1f} T (link read {link_read_share:.0%}, rotations {rot_share:.0%}).",
        f"Application: {n_q_per_dov} x c_sp = {P['t_sgn']:.0f} + lift {P['c_lift']:.1f} + reflection {P['c_outer']:.1f} + "
        f"Higgs {c_higgs_app:.0f} ({V} U_x + {higgs_sel_toff:.0f} SELECT Toffolis, outside sgn) = {t_dov:.1f} T; "
        f"normalization 2Q = {be_norm}.",
        f"R-TOL: {P['n_rot_app']} rotations per application x {n_dov} = {n_rot_shot:.6g} per shot -> eps_rot = "
        f"{eps_rot:.4g}, {t_rot:.3f} T per rotation.",
        f"N_Trotter (E26 (ii)): state-dependent second-order estimate (Alves, Lamm, Liu, preprint FERMILAB-PUB-26-0397-T), evaluated here at "
        f"weak coupling: {sd['n_modes']} transverse modes, T a = {float(a.temp_lat_2033.value):g}, Lambda_sd = "
        f"{sd['lambda_var']:.2f} a^-3 -> {sd['n_central']} steps (dt = {t_over_a / sd['n_central']:.4f} a); range "
        f"{sd['n_low']} (fastest-mode phase) to {sd['n_high']} (mean kept, Lambda {sd['lambda_full']:.1f}). Worst case "
        f"(childs2021theory, operator norms): Lambda = {tb['lambda']:.1f} a^-3 at g^2 = {g2:g}, {tb['n_steps']} steps "
        f"(dt = {tb['dt_lat']:.2e} a); minimum over g^2 {tb_min['n_steps']} at g^2 = {g2_at_min:.2f}.",
        f"Applications per step = Jacobi-Anger order at x = 2Q dt = {P['x']:.3f} with QSP error {eps_qsp:g} over the shot: "
        f"{per_step} (worst-case step {by_rule['bound']['per_step']}, Ch. 5 step {by_rule['step']['per_step']}, 30 steps "
        f"{by_rule['stated']['per_step']}); floor 2Q t = {n_dov_floor:.0f} applications = {t_shot_floor:.3g} T. Assumes a "
        f"self-inverse block encoding; otherwise {P_not_self_inverse['t_shot']:.3g} T.",
        f"Shot: {n_dov_evol} evolution + {n_dov_tpq} TPQ = {n_dov} applications x {t_dov:.1f} = {t_shot:.6g} T, "
        f"{ratio_env:.3g}x the 1e9 envelope; evolution {t_evol / t_shot:.2%}, TPQ {t_tpq / t_shot:.2%}. Range "
        f"{P_low['t_shot']:.4g}-{P_high['t_shot']:.4g}; worst-case step {by_rule['bound']['t_shot']:.4g}; Ch. 5 step "
        f"{by_rule['step']['t_shot']:.4g}.",
        f"Not priced: the gauge terms of a step, about {t_gauge_step:.3g} T ({gauge_step_share:.1%} of the step; KS "
        f"multiplicities with the group FFT, estimated here), the Higgs terms, and the amplification reflections "
        f"({n_reflections} x 2 x {n_reflect_controls} Toffolis = {t_fpaa_reflections:.3g} T per shot, arxiv_2407_17966).",
        f"Record, second-quantized construction at the worst-case step (r21 box): c_be {SQ['c_be']:.1f}, "
        f"{SQ['t_dov']:.4g} T per D_ov, {SQ['t_shot']:.4g} T per shot, {sq_lq_anc} ancilla, {sq_lq} LQ; "
        f"{SQ['t_dov'] / P_r22_bound['t_dov']:.1f}x the valid application (r22 inputs).",
        "Records (T/shot, pre-r21 query count and LCU line): " + ", ".join(f"{k} {v:.4g}" for k, v in t_shot_records.items()) + ".",
        f"At kappa 1e3-1e4 the same chain gives {t_phys[0]:.3g}-{t_phys[1]:.3g} T.",
        f"eps_l at {faults} faults per shot (R3) = {eps_l_implied:.3g}.",
        f"Coherent synthesis would move one application by {t_dov_shift_coh:.1%} (rotations are {rot_share:.1%} of a query).",
        f"Haar-start TPQ rounds = exp(0.05 x 1e3 ln2 / 2) = {haar_rounds:.1g} (app06 '~3e7').",
        "Thermal route: TPQ (r21 ruling (2)). Gibbs sampler (E23 (3), estimated here, re-priced r22): " + "; ".join(
            f"{k} {v['n_jumps']:.0f} jumps x {v['dov_per_jump']:.0f} applications = {v['t_total']:.3g} T ({v['over_tpq']:.2g}x TPQ)"
            for k, v in gibbs_est.items()) + ".",
        f"Ancilla: whole shot {sum(anc_stack.values())} + per application {sum(anc_app.values())} + term register "
        f"{sum(anc_query.values())} + link read {sum(anc_select.values())} = {lq_anc}; asserted {lq_anc_asserted} before "
        f"r20. Register {lq} LQ. The readout control is idle during the TPQ filter (one spare at the peak).",
        f"2O = single-qubit Clifford group mod phases ({a.clifford_2O.value}): the link needs no rotation synthesis; on a "
        f"color qubit the multiplexer is {mux_ix['t_direct']} T.",
        f"Condensate arm (r24 (1)): static, measured on the prepared state; shot = {C['n_dov']} TPQ applications at its "
        f"own R-TOL ({C['n_rot_shot']} rotations, eps_rot {C['eps_rot']:.4g}, {C['t_rot']:.3f} T) = {t_shot_cond:.4g} T, "
        f"eps_l {eps_l_cond:.3g}, {wall_shot_cond:.4g} s/shot.",
        f"Wall time (R4; one machine, serial): condensate {shots_cond:.0f} shots, {wall_cond_point_yr * 365.25:.2f} days/point, "
        f"{wall_cond_scan_yr:.3f} yr for {n_pts:.0f} points; sphaleron {wall_shot:.4g} s/shot, {shots_sph:.0f} shots, "
        f"{wall_point_yr:.3f} yr/point, scan {wall_serial_yr:.4g} yr (reduction {reduction_needed['scan']:.1f}x to fit "
        f"{horizon_yr:g} yr; {t_shot_to_fit_scan:.3g} T per shot).",
        "Levers for one sphaleron point under a year (r24, not applied): " + ", ".join(
            f"{k} {v:.3g}x ({wall_yr_levers[k]:.2f} yr)" for k, v in horizon_levers.items() if k in wall_yr_levers)
        + f"; record: link-select lines held {horizon_levers['hold']:.3g}x ({hold_saved} Toffolis per query, {hold_qubits} "
        "more qubits, estimated here). Spectral amplification against the 2Q t floor: not priced.",
        f"Records (r24): r22 headline (kappa {kappa_r22:g}, t = {t_over_a_r22:g}a) {P_r22['t_shot']:.6g} T; kappa {kappa_r22:g} "
        f"at t = {t_over_a:g}a {P_kappa_r22['t_shot']:.4g} T.",
        f"r25 (a): Yukawa {a.yukawa_placement.value}: {n_sgn} sign sequences, {n_q_per_dov} queries per application, "
        f"{t_dov / P_r24['t_dov']:.3f}x the r24 application; r24 shot {P_r24['t_shot']:.5g} T, condensate {C_r24['t_shot']:.5g}.",
        f"r25 (R4): 5a shot {sd_short['n_central']} steps x {P5['per_step']} + {n_dov_tpq} = {P5['n_dov']} applications, "
        f"{P5['t_shot']:.5g} T, {wall_shot_run['sph5']:.5g} s.",
        f"r25 depth ({depth_variant}, +{anc_fanout} LQ): F* " + ", ".join(
            f"{k} {v[0]:.2f}-{v[1]:.2f}" for k, v in f_star_runs.items()) + "; as itemized " + ", ".join(
            f"{k} {v[0]:.2f}-{v[1]:.2f}" for k, v in f_star_register.items()) + ".",
        f"Method test (author ruling 2026-10-05, no rate at 3^3): first result {first_shots['sph10']:.0f} x 10a + "
        f"{first_shots['sph5']:.0f} x 5a (5a variance ratio (S10/S5)^2 = {var_ratio5:.3f}) at {eps_first:.0%} = "
        f"{first_yr:.3f} yr; campaign = the test + condensate at {n_temp[0]}-{n_temp[1]} temperatures at {eps:.0%}: "
        f"{camp_yr[0]:.3f}-{camp_yr[1]:.3f} yr ({camp_reduction[0]:.2f}-{camp_reduction[1]:.2f}x the {horizon_yr:g}-yr "
        f"horizon). Hadamard-test readout instead: {shots_first_hadamard[10.0]:.3g} shots at 10a. r25 rate campaign "
        f"(record): {wall_camp_r25_yr[0]:.2f}-{wall_camp_r25_yr[1]:.2f} yr.",
    )
    eps_l = (eps_l_implied, eps_l_implied)
    return Result(era="2033", lq=(lq, lq), hard_ops=(t_shot, t_shot), breakdown=breakdown, intermediates=inter,
                  shots=(first_shots["sph10"], first_shots["sph10"]),   # 2026-10-05: the method test's 10a shots (r24-r25: 1e4)
                  wall_time_s=inter["wall_campaign_s"],   # r25: the campaign wall
                  epsilon_l=eps_l, notes=notes)


def model(a: Assumptions, era: str) -> Result:
    """Return the Result for '2028' | '2033'. Ch. 9 has no codesign box."""
    if era == "2028":
        return _model_2028(a)
    if era == "2033":
        return _model_2033(a)
    if era in ERAS:
        raise ValueError(f"Ch. 9 has no {era!r} box")
    raise ValueError(f"unknown era {era!r}")


# Eras where the chapter's own stated inputs do not reproduce the printed box.
# Empty since round B (2026-09-28).  Both entries were removed because the model now
# reproduces the printed box within PUBLISHED.rel_tol:
#   2028 (was R6): grouping 'color+s', Toffolis at 7 T, N_Trotter = 20 -> 1.11e5-1.245e5 T/shot
#                  vs the box's 1.1-1.25e5.  Since R-HOP (2026-09-29): 2.233e5-2.369e5 vs 2.2-2.4e5.
#                  Since R-TOL (2026-09-29): 1.951e5-2.046e5 vs 1.95-2.05e5.
#   2033 (was R5): 7 T throughout, one U_x per vertex at the center; band floor = vertex at
#                  0 T, top = two U_x -> 3.41e9-2.62e10 vs the box's 3.4e9-2.6e10
#                  (R-TOL: 3.383e9-2.620e10, same print).
#   r17 (2026-10-01, rulings e/f): 2028 2.76e5-2.97e5 vs the box's 2.8-3.0e5; 2033 band 3.383e9-1.2952e12
#                  vs 3.4e9-1.3e12, center 6.891e11 vs '~6.9e11'.
#   r19 (2026-10-01, E21): 2028 unchanged; 2033 band unchanged, center 5.672e11 vs '~5.7e11'.
#   r20 (2026-10-01, E23): 2028 LQ 247; 2033 LQ 1015, band 1.645e10-1.2952e12 vs 1.6e10-1.3e12.
#   r21 (2026-10-01, "promote / do TPQ / settle them"): 2028 268,995-288,535 vs the box's 2.7-2.9e5; 2033 one
#                  headline, 1.809e13 vs '~1.8e13', LQ 1005.
#   r22 (2026-10-01, E26): 2028 LQ 240 (T unchanged); 2033 single-particle architecture and the state-dependent step
#                  count, 2.008e11 vs '~2.0e11', LQ 1010.
DISPUTED: dict[str, str] = {}

PUBLISHED = {
    "2028": Published(lq=(250, 250), hard_ops=(2.7e5, 2.9e5),
                      src="app06 2028 box ('Logical qubits & 250 = 192 + 16 + 42', r25 R10: three spatial-hop groups per "
                          "slot, weight bits in parallel, F* 10.4-13.0; r22-r24 '240 = 192 + 16 + 32'). 'Per-shot T & ~2.8e5 T-gates "
                          "(2.7-2.9e5)'); lower end unchunked (268,995), upper end chunked to k <= 32 (288,535), the end that "
                          "fits the 32 ancilla. r22 (E26): 31 phasing ancilla + 1 RUS (was 37 + 2 = 39, 247 LQ); the T "
                          "count does not move",
                      rel_tol=0.10),
    "2033": Published(lq=(1015, 1015), hard_ops=(1.9e10, 1.9e10),
                      src="app06 2033 box ('1015 = 216 + 486 + 270 + 43', sphaleron shot at 10a '1.9e10'; verifier 2026-10-04: "
                          "17 filter calls in 8 rounds (2k + 1), 1513 applications, 10,063 x t_dov = 1.8597e10; referee C1 "
                          "(2026-10-04): the thermal filter at beta||H|| = 2Q / (T a) = 864, d_beta 89, 712 applications, "
                          "9262 x 1,847,860 = 1.7115e10 '1.7e10'; r25 '1.6e10': Yukawa "
                          "inside (a second sgn sequence, 417 queries per application, 8790 x 1,847,747.1 = 1.6242e10) "
                          "and the +5-LQ fan-out (R10). r24 box '1010 = 216 + 486 + 270 + 38', '8.2e9'; r24: kappa = 30, "
                          "t = 10a, 450 steps x 19 + 240 = 8790 applications x 937,592.7 T = 8.2414e9; the condensate arm "
                          "row '2.2e8' is the preparation alone). r22-r23 box: '~2.0e11': 65,643 applications x 3,058,805 T, "
                          "carried exact and rounded once (R4). r22 (E26): the sign function on the single-particle kernel "
                          "(arxiv_2607_28524; 598 Toffolis + 24 T + 6 rotations per query, 692 queries + lift + Higgs per "
                          "application), N_Trotter from the state-dependent estimate (5031; range 1186-14,242), 13 "
                          "applications per step, 38 ancilla. In the prose, not the box: the range 1.1-3.5e11, the "
                          "worst-case bound 8.3e4 steps / 1.3e12, Ch. 5's fixed step 1.1e11. r21 box was 1005 LQ, ~1.8e13",
                      rel_tol=0.10),
}


def INSTANCE_ROWS(a: Assumptions, era: str, r: Result):
    """Rows for resources.json: [(label, (lq_lo, lq_hi), (t_lo, t_hi), {extra}), ...]."""
    if era == "2028":
        return [("1+1D Z3 domain-wall, L=8, L5=4 (condensate)", r.lq, r.hard_ops,
                 {"group": "Z3", "shots": r.shots[0], "route": "domain-wall",
                  "f_star": r.intermediates["f_star"],
                  "t_box": PUBLISHED["2028"].hard_ops})]
    if era == "2033":
        i = r.intermediates
        return [("3D 3^3 SU(2)/2O overlap + Higgs, TPQ (Chern-Simons growth, method test)", r.lq, r.hard_ops,
                 {"group": "2O", "shots": r.shots[0], "route": "overlap+QSP",
                  "sgn_architecture": i["sgn_architecture"], "trotter_rule": i["trotter_rule"],
                  "t_representative": i["t_per_shot"],
                  "t_condensate_arm": i["t_per_shot_condensate"],
                  "kappa": i["kappa_2033"], "t_evol_lattice": i["t_over_a_2033"],
                  "t_range_step_count": i["t_per_shot_range"],
                  "t_alternative_fixed_step": i["t_per_shot_empirical_step"],
                  "t_alternative_worst_case_bound": i["t_per_shot_bound"],
                  "t_record_second_quantized": i["second_quantized_record"]["t_per_shot"],
                  "t_short_time": i["t_per_shot_short"],
                  "wall_first_result_s": i["wall_first_result_s"], "wall_campaign_s": i["wall_campaign_s"],
                  "f_star": i["f_star"],
                  "t_box": PUBLISHED["2033"].hard_ops})]
    return []
