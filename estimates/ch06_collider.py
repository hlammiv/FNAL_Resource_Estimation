"""Ch. 6 — QCD for colliders (label ch:collider). Reproduces the LQ and hard-op numbers of
applications/app10_collider_physics.tex, the 2026-08-28/29 rewrite (in-tree since 2026-09-29; author
rulings on its framing and instances in ch6_rewrite/CHANGE_MAP.md).

q1 (author rulings A and B, 2026-10-05; "Drop the EPOQ unless it really adds something"): E-rho-OQ is dropped from the
chapter and the model (the r27 G8 penalty exports, eroq_* assumptions and the 1e4-T thermal_prep_t placeholder are
retired). The two 2033 dipole media are quench-prepared: electric vacuum, coupling ramp over one lattice unit (20 steps,
priced as a floor), then c / T of evolution with c = 2-2 pi (90-282 steps; Hayata-Hidaka arxiv_2011_09814 for 2 pi, the
exact small-lattice checks in scratchpad/quench for 2) before the sources enter; ETH / dynamical typicality assumed at
the chapter couplings (eth_assumed). Each circuit is priced whole with its own R-TOL eps (pricing basis
scratchpad/quench/cost/price2.py). Static Im V: 2.3 / 6.7e9 -> 5.4e9-1.5e10 T per shot (5.4-15x the 1e9 limit, stated;
eps_l 1.5e-11 -> 6.6e-12), 653 shots, one temperature 34 -> 57-98 d, three 102 d -> 172-293 d. Light-like q-hat:
7.2e7 / 2.0e8 -> 2.0e8-5.5e8 T (fits; eps_l 4.9e-10 -> 1.8e-10), 405 shots, 0.65 -> 1.2-2.3 d, three 1.9 -> 3.7-6.8 d.
Shots and LQ unchanged; the species row does not move. The coarse-step thermalization lever (2-5x) is unpriced.
STAGE stays r26 for the species headline; resources.json moves at Final (dipole rows).

s648 (author ruling 2026-10-05, (5); NEEDS_AUTHOR 'Ch. 6 s648'): the static dipole row follows Ch. 5's 2033 box to
pure-gauge Sigma(216x3) (11 qubits/link, Gustafson_S648_inprep) with H_I: 833 -> 995 LQ, 2.9 / 8.4e8 -> 2.3 / 6.7e9 T
per shot (2.3-6.7x the 1e9 limit), 72 / 203 -> 85 / 240 steps (state-dependent rule on the tree-level improved free
dispersion, an assumption about H_I's free limit), eps_l 1.2e-10 -> 1.5e-11, 4.3 -> 34 d first result, 13 -> 102 d
for three temperatures; naive E-rho-OQ 378 -> 455 (|G| = 648; the printed D4 range does not move). The species /
fragmentation row and the 3+1D stretch stay Sigma(72x3) / H_KS: at 11 qubits/link the box is 1092 LQ > 1000 (would
need the ~1.1k caveat), and a switch would reuse Ch. 5's Sigma(72x3) placeholder hop (hop_group_2033), 72% of the step,
so most of the row's price would be a placeholder. The 2O light-like row does not move.

r27 additions (chapter-open pass 2026-10-05; no headline number moves, so STAGE stays r26 and resources.json does
not change). Derived here, exported as intermediates: G8 the E-rho-OQ penalty on the dipole lattices
(eroq_log10_penalty_dipole_*: naive 378 / 454; the 4-link D4 value of arxiv_2001_11490 carried per link, 40.5-44.9 /
67.5-74.8; eroq_shot_penalty Assumed 1, as Ch. 5; ALL RETIRED by q1). J2 the box term in the trend (Lund yo-yo, arxiv_2203_11601:
full extension at p / kappa, 3.9 / 7.7 / 11.6 windings of the 0.8 fm axis) and the two-box check (L_par = 16 holds
all three modes, 1828 LQ; L_par = 4 shares 3.1 GeV, 532 LQ, 2.4e9-1.25e10 T, 0.24-3.7 yr; in the campaign since
smaller ruling (a), 2026-10-05: campaign 1.5-22 -> 1.7-26 yr, 9.5e3-1.9e4 -> 1.3e4-2.5e4 shots, 4.4x -> 5.2x the horizon;
the per-shot headline does not move, so STAGE stays r26; Result.shots / wall_time_s and resources.json move).
J3 the N_h insertion: no extra T (the backward ramp is already in the shot); momentum sum rule gives
>= (sum sqrt z_i)^2 / N_z = 4.70x the shots if p_i = g_i, (sum z_i^(2/3))^3 / N_z = 23.1x at the variance bound
(ff_nh_shot_factor; not applied). J1 V/PS and calibration: no number (statements in the chapter).
STAGE: r26 (referee report 2026-10-04, items J1-J4, G7, G8; editorial_review/responses/ch06.md).
  J1  the species readout is the map N_h = U Pi_h U^dagger (arxiv_2406_05683; adiabatic readout of arxiv_1111_3633):
      each 2033 shot ends with the preparation ramp run backwards (readout_ramps_2033 = 1), then one Fourier layer to
      the electric basis (t_readout_fourier_2033, < 1e-3 of the shot, not added). 2033: 4.90e9-2.54e10 T (was
      3.90e9-1.51e10, kept as t_total_2033_r25). V/PS needs a symmetry-projected readout (open; no number).
  J2  the source is a gauge-invariant wave packet sum_d f_p(d) psibar_x W_{x,x+d} psi_{x+d} over separations up to
      L_par links (<= L_par bilocal insertions, t_source_wavepacket_2033, < 2% of a step; not added; r26 v1 priced a
      single nearest-neighbor bilinear, which has no relative-momentum profile). Box modes are multiples of
      2 pi / L_par = 1.55 GeV (source_momenta_box_modes_GeV); the 0.8 fm axis recollides at 0.4 fm/c
      (recollision_time_fm); separation through 1/Lambda needs L_par = 20, 2260 LQ.
  J4  the dipole rows measure P_s = <|w|^2>, not |W|^2. The shots are priced on the Markovian color-floor curve
      P_s = P0 (1/N^2 + (1 - 1/N^2) e^{-Gamma t}) and corrected with that model (r26 v2: 653 and 405 per temperature;
      the |W|^2-shape counts 615 and 250 are records, shots_per_T_W2). q-hat is the fundamental-source q-hat,
      <W_F> = e^{-q L r^2/4}; if the Assumed 2 GeV^2/fm is an adjoint value the SU(2) exponent takes C_F/C_A = 3/8
      (shots_qhat_if_adjoint, 695).
  G8  the dipole shot counts are binomial at an exact thermal state (thermal_variance_factor = (1, 1)); the E-rho-OQ
      factor (~<delta_ij>^-2) is open. The r25 x40 is kept as a record (thermal_variance_factor_r25).
Before r26: r25 (author rulings R1-R10, H. Lamm 2026-10-02; apply_log/r25_ch06.md).
  R1  two tiers at 2033: a minimal first result at 30% statistical error (one sampling configuration; the two
      dipole rows of R7) and the full campaign at the physics target, both with one-machine walls.
  R2  the campaign target is the 3-momentum trend at N_sigma = 3 (a 25% trend, Assumed), and the ratio shots carry
      conserved-pair clustering: N = (c / n_s + 1 / n_M) / delta^2, c = 1 (pair rate) to 2 (B + Bbar multiplicity).
  R4  the evolution window is 30-50 a (3-5 fm/c, 3-5 traversals of the 0.8 fm box; shot audit), not 100 a; the ramp
      range 1e2-1e3 is unchanged. The dipole rows price each shot at its own evolution time.
  R7  2033 first-result rows: the static color-singlet dipole (Im V(r = a, T)) on pure-gauge Sigma(216x3) / H_I 3^3
      (s648; was Sigma(72x3) / H_KS), and the
      light-like dipole q-hat on pure-gauge 2O 3x3x5 (simplobs.json). The full q-hat (the 2404-LQ 3+1D lattice) is
      past 2033.
  R9  t_gate_s (was gate_time_s), shot_overhead_s; the seven depth exports (common.depth_exports).
  R10 2028: two hop links run at once (+14 workspace qubits, 216 -> 230 LQ) so F* >= 10 (was 1.65 at the 8-ancilla
      register). 2033: F* < 10 at the worst depth estimate, so the walls are the corrected ones,
      shots x (max(N_T t_gate, D_T t_r) + t0).
Before r25: r23 (E27, H. Lamm 2026-10-02, "we don't want to anywhere assume we have multiple machines"; apply_log/r23_ch06.md).
Every wall time is serial on ONE machine at the report convention: 1 us per T-gate plus ~0.1 ms per shot (Ch. 1,
main-overview:48; shot_overhead_s, which moves no printed number). The machine count 'machines_for_horizon_ff' is retired.
Where a campaign runs past the 5-year horizon the model gives what fits in 5 years on one machine and the cost
reduction needed to fit (serial years / 5), an algorithmic / co-design reduction: ff_cost_reduction_for_horizon,
ff_shots_in_horizon, ff_calls_in_horizon, campaign_sampling_yr_at_10pct / _at_5pct, sampling_configs_in_horizon_at_5pct.
No per-shot T or LQ moves.
r22 (E26, H. Lamm 2026-10-01, "do 2"; apply_log/r22_hwp_chapters.md): a Hamming-weight-phasing group of k holds k - w(k)
ancilla, the same number as its Toffolis, and repeat-until-success synthesis holds 1. The 2028 register is
64 + 144 + 8 = 216 LQ (8 = 7 for the k = 9 phasing group + 1 for synthesis; was 11 and 219, printed '~220'). No T moves.

r19 (apply_log/r19_ch06.md). Author rulings E21 (H. Lamm 2026-10-01; apply_log/r19_core.md): (1) "make the
switch": the Sigma(72x3) hop's color squish and parity flags are computed once per link and held through V(x), V(y),
hop, V^dag(x), V^dag(y) (groups share="link", the new default), the frame undo kept (undo=True); the r17 per-application
recompute (share="draft") is the conservative sensitivity. (2) n_spin = 2 for Wilson d=3 (no effect here: staggered).
(3) SU(3) color rotations = 2 N_angles (unchanged count). The rotation count, so every R-TOL eps, is unchanged; only the
g-only Toffolis move (14,314 -> 3,672 per link). The Z3 hop has no squish and does not move, so 2028 is unchanged.
2033: 1.622e10-2.974e10 -> 1.097e10-2.021e10 T; 3+1D: 4.669e10 -> 3.238e10 T.
Before r19: r18 (apply_log/r18_ch06.md). E20 (H. Lamm 2026-10-01, "do 1"): every rotation in every circuit at the full fit
1.15 log2(1/eps) + 9.2, the gauge tables included (groups.PrimitiveCost.t: the paper's constant + (log coefficient / 1.15)
rotations at the full fit); the papers' slope-only price is the *_papers record. 2028 steps (H. Lamm 2026-10-01, "you can
do the 1 ramp+2 evolutions if its only 1.27"): 1 ramp + 2 evolution steps, 132,174.9 T = 1.32x the 1e5 reference at the
E20 price, accepted as an overshoot; the E4 step-size cross-check (t = 0.2 a at two step sizes) is restored.
Before r18: r17 (apply_log/r17_ch06.md): ruling (e) (H. Lamm 2026-10-01, "adopt the fermion primitive numbers ... estimate
the missing pieces"). The hop is groups.hop_link_cost (Sigma(72x3): the authors' unpublished color-diagonalized hop,
FermionPrimitives_unpub, with the pieces it leaves out estimated in groups.py; Z3: the R-HOP structure), the mass
groups.staggered_mass_site, and their rotations at the full fit 1.15 log2(1/eps) + 9.2 at the circuit's R-TOL eps.
r17 fitted 2028 to 1 ramp + 1 evolution step (E3, 'as many as fit'; 84,027.9 T at the r17 price).
Before it: R-TOL and the R11 rulings (apply_log/rtol_ch06.md, r11_ch06.md; earlier ch06_roundE.md,
ch06_twosteps.md, ch06_roundD.md, ch06_rewrite.md).
Report-wide rulings:
  R1  registers as the published circuits state them: Z3 2 qubits/link, Sigma(72x3) 9 (groups.link_qubits)
  R2  every term priced per LINK (the mass per SITE) per Trotter step, papers' tables (rotations at the full fit
      since E20), d = 2 (3 for the
      3+1D extension); the Sigma(72x3) electric term at the paper's stated lower bound (a floor)
  R3  epsilon_l = 0.1 expected logical faults per shot
  R4  carry exact, round the printed result once
  R5  7 T per Toffoli; Hamming-weight phasing with measurement-based uncompute (Ch. 9 ruling)
  R8  wall time at 1 us per T-gate plus ~0.1 ms per shot, serial on one machine (convention of Ch. 1,
      main-overview:48; E27: no machine count anywhere)
Round-D rulings (H. Lamm, 2026-09-29):
  (1) the 2028 benchmark has ONE Trotter step of adiabatic ramp, priced at the per-step cost
  (2) SUPERSEDED 2026-09-29 by the ruling 'Ch. 6 2028 steps' (H. Lamm): the 2028 benchmark is 2 steps in all,
      1 ramp + 1 evolution, so that it fits the 1e5 first-generation reference (99,509 T, 0.995x). The earlier
      round-D ruling (4 steps: 1 ramp + 3 evolution, 1.990e5 T, 2.0x) is kept only as the round_d_* record
  (3) the hopping term priced by the same sources and rules as every other term (rule R-HOP below)
  (4) Z3 gate costs from a source: the dihedral primitives of arxiv_2108_13305 with the Z2 dropped, on the
      circuits of arxiv_2408_00075 (groups.py)
  (5) the chapter back to 8 body pages by removing redundancy
Round-E rulings (H. Lamm, 2026-09-29, TRACKED_CHANGES.md 'RULINGS, Ch. 6 round E'):
  (E1) the structured Z3 Fourier transform is adopted: groups.GROUPS['Z3'].primitives['U_FFT'], 4 T + 2 rotations,
       explicit circuit derived here and verified; structure as arxiv_2409_17349 Eqs. (18)-(20). The transpiler
       U_F (14 rotations, arxiv_2408_00075) stays in groups.py for the record
  (E2) [tolerance SUPERSEDED by R-TOL below] every synthesized rotation at eps = 1e-3 (was 1e-4), both eras; randomized synthesis, so the
       per-shot synthesis error is N eps^2 (common.eps_per_rotation 'incoherent'), stated beside the 0.1-fault budget
  (E3) the 2028 step count: as many as fit (expected 3 = 1 ramp + 2 evolution). AT eps = 1e-3 THREE STEPS ARE
       100,338.8 T, 1.0034x THE 1e5 REFERENCE; TWO FIT (66,892.5 T). Three are priced, as the author's request of
       round E assumes, and the 0.3% overrun is printed (steps_in_cap_2028 = 2 records it). UNDER R-TOL THREE FIT
  (E4) later times by combining runs at different Trotter step sizes (dt_multi_over_a)
Ruling R-TOL (H. Lamm, 2026-09-29, TRACKED_CHANGES.md 'RULINGS, R-TOL and overshoot'; supersedes the tolerance of E2):
  total synthesis error eps_syn = 1e-2 per shot, report-wide; each circuit sets eps_rot = sqrt(eps_syn / N_rot)
  (common.eps_rot_for), N_rot = its synthesized rotations per shot (synth_rotations_per_step x steps; the 4 exact T
  of the Z3 U_FFT are not rotations). The slope-only 1.15 log2(1/eps) is kept. 2028: N = 3,792, eps = 1.624e-3,
  97,288 T (0.973x; the E3 overrun is gone). 2033: each end of the ramp range is its own circuit, N = 8.21e7 / 1.49e8,
  eps = 1.104e-5 / 8.187e-6, 2.768e9-5.106e9 T. 3+1D: N = 2.24e8, eps = 6.677e-6, 9.537e9 T.
  [r17 retired the slope-only price for the fermion terms, E20 (r18) for the gauge tables; the R-TOL rule itself
  stands. r18: 2028: N = 3,792 (3 steps), eps = 1.624e-3, 132,174.9 T (1.32x, accepted). 2033: N = 2.88e8 / 5.24e8,
  eps = 5.89e-6 / 4.37e-6, 1.622e10-2.974e10 T. 3+1D: N = 7.86e8, eps = 3.57e-6, 4.669e10 T.]

HOP SINCE r17: groups.hop_link_cost. For Sigma(72x3) the per-link counts come from the authors' unpublished
  Fermion_Primitives draft (color frame V_g on both sites and its undo, eigen-class squish, controlled diagonalizers,
  one HWP phasing group; since r19 (E21) the color squish and parity held once per link: 3,672 Toffoli + 360 T +
  2,933 rotations per link for the 3 fields; r17-r18 recomputed them per frame application, 14,314 Toffoli), with the
  frame undo (drawn in the draft, not counted there), the SU(3) squish uncompute, the MBU ladders and the phasing
  estimated there (status SCALING). For Z3 (no draft
  entry, diagonal link) the R-HOP structure below at the full fit: 56 Toffoli + 32 rotations per link.
  R-HOP is kept as the record (hop_link_rhop).

RULE R-HOP, RETIRED BY r17 (derived here; no published circuit priced a gauge-covariant staggered hop on a discrete-group link)
  Per link per Trotter step:  T_hop = n_moves x C_move(G) + n_P(G) x HWP_T(k),   k = N_c N_stag
    - group-element basis: the link operator is diagonal, so the hop needs no Fourier transform
    - C_move = 0 where D(g) is diagonal on the register (Z3 as the center of SU(3): D(g) = omega^g on all
      colors); otherwise each field is moved across the link and back, n_moves = 2 N_stag, each move priced
      at one U_mul (Ch. 9 ruling VERTEX), status SCALING (no compiled circuit; not a floor)
    - n_P = Pauli strings per color-flavor copy: 8 for the Z3-dressed hop on the arxiv_2408_00075 register
      after one Clifford CNOT (derived here, checked in tests), 2 (XX+YY) after a move
    - HWP_T(k) = (floor(log2 k)+1) rotations + (k - w(k)) Toffolis, one group per (link, Pauli string)
      [arxiv_1709_06648, arxiv_1902_10673; Ch. 9's rule, common.hwp_group], same convention
      as the gauge terms (7 T/Toffoli, 1.15 log2(1/eps) per rotation)
  Shared code (report-wide ruling, H. Lamm 2026-09-29): the rule is groups.rhop_link, the phasing
  common.hwp_group; hop_link and hwp below call them with this chapter's inputs (apply_log/rhop_shared.md).
  Staggered mass: one Z per copy, one HWP group of k per site.
  Current insertion: NOT a per-step term. Readout (a) injects by quench (none); readout (b) inserts the
  bilocal current once or twice per shot (derived here: 16 U_mul for an 8-link Wilson line, its uncompute
  and the moves, plus a controlled bilinear), recorded as an intermediate.

Instances
  2028      species-count pipeline, 2+1D Z3 subset SU(3), 4^2, 3 staggered fields     (box, app10:114-132)
  2033      one circuit, two readouts: (a) hadron-species ratios by direct sampling,
            (b) FCC-clock fragmentation D^q_pi(z); 2+1D Sigma(72x3), 4x8 cut from 4x16,
            3 staggered fields                                                      (box, app10:136-161)
  codesign  NOT A BOX: the 3+1D L_s=4 extension priced in the plug-in prose (app10:95) and
            reserved for the heavy-ion / q-hat stretch. Filed under era 'codesign' because
            common.ERAS has no 'stretch'.

The chapter's cost model (app10:69-86)
    N_q          = q_G V n_link + N_c N_stag V + N_anc           (eq:Nq_collider)
    N_gate^shot  = N_Trot V n_link polylog(1/eps)                (eq:Ngate_collider)
    N_shot^samp  = (c / n_s + 1 / n_M) / delta_rel^2             (eq:Nshot_samp_collider; R2)
    N_shot^FF    = A^-2 eps^-2 N_z N_tens N_boost                (eq:Nshot_collider)
and per Trotter step
    T_step = n_links x (magnetic + electric + T_hop per link) + V x T_mass per site

INPUTS AND SOURCES
  q_G(Z3) = 2                      GROUPS["Z3"].link_qubits; register of arxiv_2408_00075 = ceil(log2 3) of arxiv_2108_13305
  q_G(Sigma(72x3)) = 9             GROUPS["S72x3"].link_qubits (compiled register, arxiv_2511_17437); app10:69,90,140
  convention                       7 T per Toffoli, 1.15 log2(1/eps) + 9.2 T per rotation (gauge tables E20, hop and
                                   mass r17; slope-only before),
                                   eps from R-TOL per circuit (app10:88; was 1e-3 in round E, 1e-4 before); H_KS
  Z3 primitives                    groups.py: D_N primitives of arxiv_2108_13305 with the Z2 dropped, circuits of
                                   arxiv_2408_00075; the Fourier transform is U_FFT (4 T + 2 rotations, derived here,
                                   round E), not the transpiler U_F (14 R_Z); checked in tests/test_ch06_collider.py
  Sigma(72x3) primitives           groups.py, arxiv_2511_17437 tab:tgatecost; U_FFT is the stated LOWER BOUND (:905)
  multiplicities                   groups.PRIMCOST (tab:primcost, group-independent)
  hopping and mass                 groups.hop_link_cost, groups.staggered_mass_site (r17), app10:88,95;
                                   FermionPrimitives_unpub for Sigma(72x3)
  2028 lattice                     D=2, L=4 (V=16, 32 links), N_c=3, N_stag=3, N_anc=22 (2 hop links x (7 HWP + 4 RUS); R10)
  2028 steps                       1 ramp + 2 evolution (round-D ruling 1; r18 ruling: 1.32x accepted), app10:88,124
  2028 shots                       1e3, ~10%/channel at n_s ~ 0.1 (app10:128; R11 item 3 (a))
  2033 lattice                     D=2, 4x16 uncut, cut to 4x8 (V=32): register and T both priced there (R11 (A));
                                   L_par=10 (V=40) a sensitivity parenthetical; N_anc=100   app10:95,139,141
  2033 steps                       300-500 evolution (t = 30-50 a; R4) + ramp 1e2-1e3 steps 'at the same cost'   app10:90,141,143
  3+1D extension                   L_s=4 (V=64, n_link=3, 192 links); 1e3 steps; no ramp  app10:90
  sampling shots                   n_s ~ 0.1 (Assumed); 30% first result; 5.9% per point for the 3-sigma trend  app10:80,146-147
  fragmentation shots              A = 0.05-0.2 (Assumed), eps = 0.2, N_z=10, N_tens=4, N_boost=3      app10:86,148-149
  campaign                         3 source momenta, one spacing; 120 FF calls; 2 dipole rows x 3 T; 5 yr app10:168-169
  systematics                      O(a^2) 10-20%, volume 30%, subgroup 5%; 20% stat per z-bin         app10:48,155-156,176
  utility                          5% of $0.4B over 10 instances                                       app10:183
  Ch. 9 counting rules             estimates.ch09_chiral_gauge: the HWP helpers of R-HOP, and the ch9_rules_*
                                   record (Ch. 9's PRE-R-HOP 2-string, 30-T accounting on this lattice, the basis
                                   of the retired 5e3). Since R-HOP Ch. 9 prices its spatial Z3 hop by the same
                                   rule as this chapter (groups.rhop_link, 8 strings, 15.28 T per hop rotation);
                                   ch9_current_rule_hop_per_link_k9 checks it gives this chapter's hop at Ch. 9's
                                   eps = 1e-4 (881 T; Ch. 6 2028 under R-TOL: 733 T)
WHAT IS NOT DERIVED HERE, OR RESTS ON A COMPARISON
  The Sigma(72x3) hop: the unpublished draft's gate tables, with the frame undo, SU(3) squish uncompute, MBU ladders
    and HWP phasing estimated (groups.py, NEEDS_AUTHOR 'Cross-chapter: fermion hop'). 16.0x the retired R-HOP.
    Recorded variants per link at the 2033 low-end eps (central 1.116e5, share='link' + undo, E21): share='draft'
    (the r17-r18 headline) 1.861e5, no undo 7.94e4, the draft's 2(n-2) ladders 1.205e5.
  The Z3 Re/Im split of the dressed hop into two commuting blocks is a Trotter split (as Kan-Nam's U(1) split).
  The Sigma(72x3) electric term is a FLOOR: no fast Fourier transform exists for the group.
  n_s ~ 0.1, A ~ 0.05-0.2: the chapter's declared working assumptions.
  ~100 linear-response ancilla, 10x qubitized trim, the 1e3 per-step target of the algorithm gap: Stated.
  RETIRED in round D: the working figures 5e3 (2028 hopping), 3.4e5 (2033 hopping and current insertion) and
    2.04e6 (3+1D); kept as Uncited retired_* inputs only to reproduce the legacy_* record.
  The dipole rows' physics inputs (q-hat, kappa/T^3, P_d, the color factor) are Assumed; the estimator-variance
    factor 40 and the E-rho-OQ prep are Ch. 5's (Stated there).
WALL TIME AND FACTORIES
  t_gate_s = 1e-6 s per T, shot_overhead_s = 1e-4 s per shot, ONE machine, serial (report convention).
  T-depth per shot from factory.json (Ch. 6 runs, 2026-10-02): 2028 by the HWP/RUS layer count with
  workspace_links_2028 = 2 hop links at once (22 ancilla): F* = 11.0-12.4 (the 8-ancilla register gave 1.65).
  2033 (P) schedule, the priced register, hop links one at a time: per-step depth 1,907 + 205 a (low) to
  3,792 + 613 a (high) per link, a = T per rotation; F* = 7.0-18.5 per circuit, so the walls are
  shots x (max(N_T t_gate, D_T t_r) + t0) with t_r = 10 us (up to 1.43x the serial wall).
  The value to trust is f_star_per_circuit_2033 (7.0-18.5). The 2033 'f_star' export (1.8-71, baseline_ok
  (False, True)) is common.depth_exports' cross-pairing of the two ramp ends, i.e. of two different circuits; it is
  kept for the CONTRACT test and is printed nowhere.
  First result: one sampling configuration at 30% -> wall_first_result_s (dipole rows recorded beside it).
  Campaign: the 3-momentum trend at 3 sigma (3 configurations, one spacing) -> wall_campaign_s.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, fields, replace as _dc_replace

from estimates.common import (Assumed, Cited, Stated, Uncited, Tagged, Primitive,
                              CircuitStatus, Result, Published, ERAS,
                              t_per_rotation, eps_per_rotation, toffoli_t,
                              EPS_SYN, EPS_SYN_SRC, eps_rot_for, depth_exports,
                              REACTION_TIME_S, FACTORY_BASELINE)
from estimates.common import hwp_group, RUS_ANCILLA
from estimates.groups import (GROUPS, PRIMCOST, MAGNETIC, magnetic_per_link, electric_per_link,
                              DIAGONAL_LINK, rhop_link, staggered_mass_site, hop_link_cost, FERMION_DRAFT, fp_printed)
from estimates import ch09_chiral_gauge as ch9

TEX = "app10"
RFI = "DOE_RFI_2026"
S216 = "arxiv_2511_17437"
FFT = "arxiv_2408_00075"
SECONDS_PER_YEAR = 3.156e7
HBARC_GEV_FM = 0.197327
STAGE = "r26"
# DIAGONAL_LINK (D(g) diagonal on the published register: the hop needs no color move) is imported from
# groups.py, where rule R-HOP now lives as the shared function rhop_link (report-wide ruling 2026-09-29).


@dataclass(frozen=True)
class Assumptions:
    # ---- RFI anchors ---------------------------------------------------------
    rfi_lq_2028: Tagged = Cited((150, 250), RFI, "'150--250 envelope' (app10:123)")
    rfi_t_2028: Tagged = Cited(1e5, RFI, "'the 10^5 first-generation limit' (app10:88,125)")
    rfi_lq_2033: Tagged = Cited(1000, RFI, "'up to 1.2x the 1000-LQ marker' (app10:90); box 'the 2033 1000-LQ marker' (app10:144)")
    rfi_t_2033: Tagged = Cited(1e9, RFI, "'the 10^9-T resource limit' (app10:90); box 'the 10^9 limit' (app10:145)")
    rfi_eps_l: Tagged = Cited(1e-8, RFI, "'At eps_l = 10^-8 the 2+1D shot incurs 23--43 logical faults' (app10:90)")
    clean_shot_faults: Tagged = Stated(0.1, f"{TEX}:90,128,147", "'the report's criterion of 0.1 expected faults per shot' (ruling R3)")
    seconds_per_year: Tagged = Cited(SECONDS_PER_YEAR, "PDG", "Julian year")
    hbarc_GeV_fm: Tagged = Cited(HBARC_GEV_FM, "PDG", "hbar c")

    # ---- conventions (rulings R1, R2, R5, R8), stated at app10:88 ----------------
    link_width_source: Tagged = Stated("compiled", f"{TEX}:69,90",
                                       "'q_G qubits per link (the register the published circuits are built on: 2 for Z3, 9 "
                                       "for Sigma(72x3))' (ruling R1). 'compiled' reads GROUPS[k].link_qubits; 'chapter' reads "
                                       "link_qubits_chapter (8 = ceil(log2 216), the pre-ruling text) for the legacy check")
    t_gate_s: Tagged = Stated(1e-6, f"{TEX}:127,143", "'at 1 us per T-gate; Ch. 1' in both boxes (ruling R8; main-overview:48). "
                                                    "R9 (2026-10-02): the unified name, was gate_time_s")
    shot_overhead_s: Tagged = Assumed(1e-4, "per-shot overhead t0 ~ 0.1 ms (register initialization, final readout, decoding), "
                                            "the report convention of Ch. 1 (main-overview:48); as Ch. 8 and Ch. 10 "
                                            "shot_overhead_s. It moves no printed number here (2028: 0.1322 -> 0.1323 s per shot)",
                                      src="main-overview:48")
    eps_syn: Tagged = Stated(EPS_SYN, f"{TEX}:88", f"ruling R-TOL ({EPS_SYN_SRC}): 'total synthesis error of 10^-2 per "
                                                    "shot'; each circuit sets eps_rot = sqrt(eps_syn / N_rot) from its own "
                                                    "rotation count (common.eps_rot_for)")
    eps_rot: Tagged = Uncited(1e-3, "SUPERSEDED 2026-09-29 by R-TOL (was round-E ruling E2, 'eps = 10^-3 per rotation'). "
                                    "The tolerance the pricing helpers read: every era REPLACES it with "
                                    "eps_rot_for(N_rot, eps_syn) of its own circuit (_at_eps) before pricing, so this "
                                    "default reaches only the *_round_e record intermediates")
    synthesis_errors: Tagged = Stated("incoherent", f"{TEX}:88", "'randomized synthesis ... errors add as N eps^2' (ruling E2); "
                                                                 "common.eps_per_rotation(errors='incoherent')")
    eps_rot_before_round_e: Tagged = Uncited(1e-4, "SUPERSEDED 2026-09-29 (round-E ruling E2): the report's working tolerance, "
                                                   "used for the record intermediates *_at_1e-4")
    dt_multi_over_a: Tagged = Stated((0.1, 0.2, 0.3), f"{TEX}:88", "round-E ruling E4: 'runs at Delta t = 0.1 a, 0.2 a and 0.3 a'")
    synthesis_model: Tagged = Stated("rus", f"{TEX}:88", "gauge terms: the subgroup papers' tables with every rotation at the "
                                                         "full fit '1.15 log2(1/eps) + 9.2' (E20 ruling, 2026-10-01, supersedes R2's "
                                                         "papers' convention; n_rot = the paper's log coefficient / 1.15, the "
                                                         "papers' Toffoli constants kept). The papers' slope-only price is kept "
                                                         "as the *_papers record")
    fermion_synthesis: Tagged = Stated("rus", f"{TEX}:88", "r17 (2026-10-01): hop and mass rotations at the full fit "
                                                           "'1.15 log2(1/eps) + 9.2', as the Fermion_Primitives draft prices them "
                                                           "and as Chs. 5 and 10 price their fermion terms (common.t_per_rotation 'rus')")
    toffoli_convention: Tagged = Stated("textbook", f"{TEX}:88", "'7 T per Toffoli' (ruling R5)")
    hamiltonian: Tagged = Stated("KS", f"{TEX}:30,56,88", "'per link per Trotter step of H_KS'; groups.PRIMCOST key")

    # ---- shared lattice content ----------------------------------------------
    n_c: Tagged = Stated(3, f"{TEX}:69", "'one Jordan--Wigner mode per color per site, N_c = 3'")
    n_stag: Tagged = Stated(3, f"{TEX}:88,90,118,142", "'3 staggered fields' in both boxes")

    # ---- 2028 benchmark (plug-in app10:88; box app10:114-132) -----------------
    group_2028: Tagged = Stated("Z3", f"{TEX}:88,118", "'q_G = 2 (Z3 subset SU(3))'; box '2 qubits/link'")
    dim_2028: Tagged = Stated(2, f"{TEX}:118", "'2+1D'; n_link = D")
    L_2028: Tagged = Stated(4, f"{TEX}:88,118", "'on 4^2'")
    n_anc_2028: Tagged = Stated(8, f"{TEX}:88,123", "SUPERSEDED RECORD (r25 R10: the chapter now prints 230 = 64 + 144 + 22, "
                                                    "lq_anc_workspace_2028; this field is the per-link 7 + 1 and the r22 quote below is gone). "
                                                    "'+8 = 216 LQ, the 8 ancilla being one phasing group's 7 (below) and one "
                                                    "for rotation synthesis'; box '216'. = ch9.hwp_ancilla(9) + "
                                                    "common.RUS_ANCILLA = 7 + 1 (E26, r22: a group of k holds k - w(k) "
                                                    "ancilla, repeat-until-success synthesis 1; the chapter printed 11 = 7 "
                                                    "+ a 4-qubit weight register, and 219 -> '~220')")
    n_trot_2028: Tagged = Stated(2, f"{TEX}:88,124", "'One ramp step and two evolution steps'; box '1 ramp + 2 evolution steps'. "
                                                     "Ruling (H. Lamm 2026-10-01): 'you can do the 1 ramp+2 evolutions if its "
                                                     "only 1.27'. At the E20 price the 3 steps are 1.32x the 1e5 reference, "
                                                     "accepted as an overshoot (2 steps would be 0.87x). r17 had 1 evolution "
                                                     "step; R-TOL had 2 (0.973x at the slope-only price)")
    ramp_steps_2028: Tagged = Stated(1, f"{TEX}:88,119,124", "round-D ruling 1: 'one Trotter step's worth of adiabatic ramp, "
                                                             "priced at the per-step cost'; box 'one adiabatic ramp step'. The "
                                                             "chapter says one step does not reach the interacting vacuum")
    delta_rel_2028: Tagged = Stated((0.10, 0.10), f"{TEX}:128", "box 'Shots ~1e3 (binomial counts at n_s ~ 0.1, ~10%/channel)' "
                                                               "(ruling R11 item 3 (a): print the accuracy the shots buy; "
                                                               "was 20--30%). Kept as a degenerate range for the guard")
    shots_2028: Tagged = Stated(1e3, f"{TEX}:128", "'Shots ~10^3 (binomial counts at n_s ~ 0.1, ~10%/channel)'")
    t_step_gap_target: Tagged = Stated(1e3, f"{TEX}:94", "'pushing the per-step T-count below ~10^3 (a 50x cut at 2028)'")

    # ---- 2033 target (plug-in app10:90; box app10:138-165) --------------------
    group_2033: Tagged = Stated("S72x3", f"{TEX}:90,142", "'Sigma(72x3)' for SU(3)")
    dim_2033: Tagged = Stated(2, f"{TEX}:142", "'2+1D SU(3)'; the lattice is L_perp x L_par, so only 2 is accepted")
    L_perp: Tagged = Stated(4, f"{TEX}:90,142", "'4x16'; 'four-site-wide ladder' (app10:63)")
    L_par_uncut: Tagged = Stated(16, f"{TEX}:90", "'on 4x16 (V=64)'")
    L_par_cut: Tagged = Stated((8, 10), f"{TEX}:95", "sensitivity range: '(1180 LQ and 3.5--6.4e9 T at L_par=10)'. Since "
                                                     "ruling R11 option (A) the box quotes only the priced volume")
    L_par_priced: Tagged = Stated(8, f"{TEX}:95,139,141", "'cutting L_par to 8 (V=32) gives the box's 964 LQ'; box '4x8 (cut "
                                                          "from 4x16)', '964 (V=32, the volume the T-count prices)'")
    n_anc_2033: Tagged = Stated(100, f"{TEX}:86,90", "'N_anc is ~100 with the linear-response register'; '+100' in the plug-in")
    n_trot_2033: Tagged = Stated(1e3, f"{TEX}:90,143,145", "'Across N_Trot ~ 10^3 steps --- the full window for a ~2-GeV parton'")
    n_trot_eq: Tagged = Stated((1e2, 1e3), f"{TEX}:72", "eq:Ngate_collider 'N_Trot ~ t_had/Delta t ~ 10^2--10^3'")
    ramp_steps_2033: Tagged = Stated((1e2, 1e3), f"{TEX}:90,145", "'the adiabatic ramp (10^2--10^3 steps at the same cost)'. "
                                                                  "(0, 0) is accepted: a product-state start")
    readout_ramps_2033: Tagged = Stated(1, f"{TEX}:76", "referee J1 (r26): the species-number operators are "
                                        "N_h = U Pi_h U^dagger (arxiv_2406_05683 Eq. proqu; adiabatic readout of "
                                        "arxiv_1111_3633), U the preparation ramp from strong coupling, so each shot ends "
                                        "with the ramp run backwards: one more ramp per shot (derived here). 0 restores r25")
    t_had_over_a: Tagged = Stated(100, f"{TEX}:74,143", "SUPERSEDED RECORD (r25 R4: the window is 30-50 a, window_over_a; "
                                                        "the quote below is gone and this fixes only Delta t = 0.1 a and the r23 record). "
                                                        "'t_had ~ 100 a ~ 10 fm/c'; box 't = 100 a'. With n_trot_2033 this FIXES "
                                                        "the step, Delta t = t_had / N_Trot; the chapter prints no Delta t")
    a_fm: Tagged = Stated(0.1, f"{TEX}:74,142", "'at a ~ 0.1 fm'")
    source_momentum_GeV: Tagged = Stated(2.0, f"{TEX}:90,143", "'a ~2-GeV parton'; box 'p ~ 2 GeV'")
    qubitized_trim: Tagged = Stated(10, f"{TEX}:90", "'A qubitized block-encoding of H_KS trims per-shot cost a further ~10x'")
    lambda_qcd_GeV: Tagged = Assumed((0.2, 0.3), "NOT IN THE CHAPTER. Conventional range, used only to put a number on workflow "
                                                 "step 4's 't ~ 1/Lambda_QCD' (app10:56) beside the priced 10 fm/c window",
                                     src=f"{TEX}:56")

    # ---- 3+1D extension (app10:90) -------------------------------------------
    L_s_3d: Tagged = Stated(4, f"{TEX}:90", "'A 3+1D extension at L_s = 4 (V=64, n_link=3)'")
    n_trot_3d: Tagged = Stated(1e3, f"{TEX}:90", "'9.2e9 T over 10^3 steps with no ramp priced'")
    ch10_stretch_lq_printed: Tagged = Cited((2600, 3100), "app07:197", "Ch. 10's post-2033 stretch on the same lattice: "
                                                                      "'~2600--3100 LQ (1728 gauge at 9 qubits/link + 576 fermion + "
                                                                      "the 250--750 BE/bath ancilla of the 2033 route)'. Ch. 10 "
                                                                      "prints no T-count for it. Was '~2400 (+ ~50 ancilla)' before "
                                                                      "Ch. 10's R11 ruling ch10-stretch-ancilla (a)")
    ch10_stretch_n_anc: Tagged = Cited((250, 750), "app07:197", "'+ the 250--750 BE/bath ancilla of the 2033 route' in the same row")

    # ---- shots: inclusive counts (eq:Nshot_samp_collider) ---------------------
    n_s: Tagged = Assumed(0.1, "'n_s ~ O(0.1) --- the baryon/kaon scale on the chosen lattice, a working assumption' (app10:80); "
                               "CHANGE_MAP ruling 3", src=f"{TEX}:80")
    delta_rel_2033: Tagged = Stated((0.05, 0.10), f"{TEX}:43-45,80,148", "'a 5--10% determination'")
    shots_samp_printed: Tagged = Stated((1e3, 4e3), f"{TEX}:80,149", "'N ~ 10^3--4x10^3 direct samples'; box 'shots 10^3--4x10^3' "
                                                                    "(R4; the pre-ruling text printed 10^3--10^4)")

    # ---- shots: Hadamard test (eq:Nshot_collider) -----------------------------
    hadamard_amplitude: Tagged = Assumed((0.05, 0.2), "'A ~ 0.05--0.2; the resulting 25--400 dominates the budget and is measured "
                                                      "first, in the opening shots of the 2033 fragmentation instance' (app10:86; "
                                                      "CHANGE_MAP ruling 6)", src=f"{TEX}:86")
    eps_stat: Tagged = Stated(0.2, f"{TEX}:48,86", "'At eps ~ 0.2'; Objective 6 '~20% statistical accuracy'")
    n_z: Tagged = Stated(10, f"{TEX}:86,172", "'N_z ~ 10 z-bins'")
    n_tens: Tagged = Stated(4, f"{TEX}:86,172", "'N_tens ~ 4 tensor components / spin projections'")
    n_boost: Tagged = Stated(3, f"{TEX}:86,172", "'N_boost ~ 3 source momenta for LaMET extrapolation'")
    shots_ff_printed: Tagged = Stated((7.5e4, 1.2e6), f"{TEX}:86,152", "'7.5e4--1.2e6 in all'; box 'shots 7.5e4--1.2e6 (3e3 x A^-2)' "
                                                                      "(R4; the pre-ruling text printed ~10^5--10^6)")

    # ---- campaign (requirements table, app10:171-179) -------------------------
    n_spacings: Tagged = Stated(2, f"{TEX}:172", "'2 lattice spacings'")
    n_source_momenta: Tagged = Stated(3, f"{TEX}:172", "'x 3 source momenta' (CHANGE_MAP ruling 4)")
    n_qhat_T: Tagged = Stated(3, f"{TEX}:172", "'q-hat: 3 temperature points' (CHANGE_MAP ruling 5)")
    n_species_channels: Tagged = Stated(5, f"{TEX}:172", "'all ~5 species channels ... read from the same event records'")
    campaign_horizon_yr: Tagged = Stated(5, f"{TEX}:171", "'Campaign horizon 5 years'")

    # ---- systematics (2033 box, app10:158-161; Objective 6) -------------------
    syst_discretization: Tagged = Stated((0.10, 0.20), f"{TEX}:158", "'O(a^2) ~ 10--20% at a = 0.1 fm'")
    syst_volume: Tagged = Stated(0.30, f"{TEX}:158", "'~30% at L_perp = 4a'")
    syst_subgroup: Tagged = Stated(0.05, f"{TEX}:108,159", "'~5% at Sigma(72x3)'")

    # ---- clocks (overview, app10:10) and utility box (app10:186) --------------
    n_z_fcc: Tagged = Cited(6e12, "arxiv_2505_00272", "'6e12 Z bosons' at Tera-Z (app10:10)")
    n_z_lep: Tagged = Cited(1.7e7, "LEPEWWG_Zresonance_2006", "'1.7e7 at LEP' (app10:10)")
    lep_error_divisor: Tagged = Cited(300, "arxiv_2407_09593", "'the LEP uncertainty divided by approximately 300' (app10:10)")
    utility_share: Tagged = Assumed(0.05, "'a ~5% share ... the share is the stated, rescalable assumption' (app10:186)", src=f"{TEX}:186")
    cms_capital_musd: Tagged = Stated(400, f"{TEX}:186", "'~$0.4B U.S. CMS detector capital (the CMS share of the $331M ... plus the $200M "
                                                        "HL-LHC CMS upgrade)'; the CMS share of the $331M is not printed")
    n_instances_utility: Tagged = Stated(10, f"{TEX}:186", "'spread over ~10 application instances in a 5-year campaign'")

    # ---- hopping and mass: rule R-HOP (round-D ruling 3), derived here; app10:88,90 ---------------
    n_pauli_hop_diagonal: Tagged = Stated(8, f"{TEX}:88", "'the dressed hop is 8 Pauli strings per copy' (Z3). DERIVED HERE: Re and "
                                                          "Im of omega^g on the arxiv_2408_00075 register after one Clifford CNOT "
                                                          "(g0 -> g1), the forbidden state free; checked in test_z3_hop_*")
    n_pauli_hop_moved: Tagged = Stated(2, f"{TEX}:88,90", "RECORD ONLY since r17 (the retired R-HOP comparison, hop_link_rhop). "
                                                          "After a color move the hop is the plain XX+YY bilinear, 2 strings "
                                                          "(as Ch. 9's link-free hops, app06:108 'the two Pauli strings per "
                                                          "hopping bilinear (eight for a spatial hop dressed by its link, below)')")
    moves_per_field: Tagged = Stated(2, f"{TEX}:95", "since r17 used by the current insertion only ('16 group multiplications U_x "
                                                     "for an 8-link Wilson line, its uncompute and its action on one field': 2 of "
                                                     "the 16) and by the retired R-HOP record (W and W^dag per field)")
    n_pauli_mass: Tagged = Stated(1, f"{TEX}:88", "'the staggered mass is one group per site': one Z per copy")
    hop_share: Tagged = Stated("link", f"{TEX}:95", "author ruling E21 (1), 2026-10-01, 'make the switch': 'the color "
                                                    "register map and parity flags are computed once per link and held "
                                                    "through both frames, the hop and both undos' (groups share='link', "
                                                    "the shared default since r19). 'draft' (per-application recompute, "
                                                    "the r17-r18 headline) is the conservative sensitivity")
    insertion_links: Tagged = Stated(8, f"{TEX}:95", "'16 group multiplications U_x for an 8-link Wilson line, its uncompute "
                                                     "and its action on one field' = 2 x (8 - 1) + 2")
    insertions_per_shot_b: Tagged = Stated((1, 2), f"{TEX}:90", "'readout (b)'s one or two insertions of the bilocal current'")
    insertion_bilinear_rotations: Tagged = Assumed(4, "DERIVED HERE (round-D hopping research, rule_numbers.py): the controlled "
                                                      "bilinear of one insertion, 2 x (4 R_Z + 1 Toffoli); 136 T of the ~1.8e4, "
                                                      "not printed separately")
    insertion_bilinear_toffolis: Tagged = Assumed(1, "same derivation: 1 Toffoli per half, x 2")

    # ---- RETIRED in round D (ruling 3): the working figures, kept only for the legacy_* record ------
    retired_hop_2028: Tagged = Uncited(5e3, "RETIRED 2026-09-29 (round-D ruling 3): 'Staggered hopping ... carried at ~5e3 T/step, "
                                            "the chapter's working figure ... not derived here' (app10:88 before round D)")
    retired_hop_insertion_2033: Tagged = Uncited(3.4e5, "RETIRED 2026-09-29: 'staggered hopping and current insertion carried at "
                                                        "3.4e5 T/step, the chapter's working figure (no compiled circuit)'")
    retired_hop_insertion_3d: Tagged = Uncited(2.04e6, "RETIRED 2026-09-29: 'a 2.0e6 working figure for hopping and current "
                                                       "insertion' (3+1D)")
    retired_n_trot_2028: Tagged = Uncited(20, "RETIRED 2026-09-29 (round-D ruling 2): 'Twenty steps of 2.5e4 T cost 5.0e5 T'")
    round_d_n_trot_2028: Tagged = Uncited(3, "SUPERSEDED 2026-09-29 (ruling 'Ch. 6 2028 steps'): round D's 'One ramp step and "
                                             "three evolution steps ... give 2.0e5 T, 2.0x' (app10:88 before this ruling)")
    twosteps_n_trot_2028: Tagged = Uncited(1, "SUPERSEDED 2026-09-29 (round-E ruling E3): 'One ramp step and one evolution step "
                                              "... give 1.0e5 T, 0.995x' at eps = 1e-4 with the transpiler U_F")

    # ---- r25 rulings R1, R2, R4 (H. Lamm 2026-10-02): tiers, ratio shots, window ----------------
    window_over_a: Tagged = Stated((30, 50), f"{TEX}:74", "R4 (shot audit): 'the priced window is 30--50a (3--5 fm/c), three to "
                                                          "five traversals of the 0.8 fm periodic box'; was t_had = 100 a. "
                                                          "Steps = window / Delta t at the chapter's Delta t = 0.1 a")
    delta_first_2033: Tagged = Stated(0.30, f"{TEX}:135", "R1: the minimal first result at 30% statistical error (one "
                                                          "sampling configuration)")
    pair_clustering: Tagged = Stated((1, 2), f"{TEX}:80", "R2: baryon number and strangeness are conserved, so B and Bbar "
                                                          "arrive in pairs: c = 1 if n_s is the pair rate, 2 if it is the "
                                                          "B + Bbar multiplicity (shot audit)")
    n_meson: Tagged = Assumed((0.5, 1.0), "mesons per event in the ratio's denominator; the shot audit's range, giving "
                                          "(c / n_s + 1 / n_M) / delta^2 = 1.1e3-2.2e3 at 10%", src=f"{TEX}:80")
    trend_size: Tagged = Assumed(0.25, "the size of the B/M change across the three source momenta the campaign must "
                                       "resolve (shot audit's reference trend)", src=f"{TEX}:80")
    trend_sigma: Tagged = Stated(3, f"{TEX}:80", "R2: the campaign states N-sigma on the 3-momentum trend; N = 3. The "
                                                 "end-to-end change carries sqrt(2) x the per-point error")

    # ---- r25 R7: 2033 first-result dipole rows (simplobs.json, Ch. 6 candidates C' and E) ----------
    dipole_a_fm: Tagged = Stated(0.2, f"{TEX}:97", "r = a = 0.2 fm, Ch. 5's 2033 spacing; aT ~ 0.45")
    dipole_dt_fm: Tagged = Stated(0.01, f"{TEX}:97", "the 0.01 fm/c grid (Ch. 5's a_t/10); since the step rule "
                                                      "(2026-10-04) each dipole shot runs max(t / 0.01 fm/c, N_sd) "
                                                      "steps: 72 / 203 (static), 64 / 181 (light-like), dipole_rule_steps")
    dipole_times_fm: Tagged = Stated((0.5, 1.0), f"{TEX}:97", "the two times / path lengths of the log slope")
    dipole_n_anc: Tagged = Stated(100, f"{TEX}:97", "the 2033 ancilla allowance, as the box")
    static_dims: Tagged = Stated((3, 3, 3), f"{TEX}:98", "pure-gauge Sigma(216x3) on 3^3 (Ch. 5's 2033 lattice; the "
                                                         "medium here is at aT ~ 0.45, not Ch. 5's 160 MeV)")
    # s648 (author ruling 2026-10-05, (5)): the static dipole row follows Ch. 5's 2033 box to Sigma(216x3) with H_I
    # (995 LQ <= 1000). The species / fragmentation row stays Sigma(72x3) / H_KS: 1092 LQ > 1000 at 11 qubits/link
    # (would need the ~1.1k caveat), and a switch would reuse Ch. 5's Sigma(72x3) placeholder hop, 72% of the step
    # (NEEDS_AUTHOR 'Ch. 6 s648').
    group_dipole_static: Tagged = Cited("S216x3", "Gustafson_S648_inprep", "s648 ruling (5): 'pure-gauge "
                                        "Sigma(216x3) 3^3', qubit encoding, 11 qubits/link (groups.py)")
    hamiltonian_dipole_static: Tagged = Stated("I", f"{TEX}:98", "s648 ruling (1)/(5): the improved H_I, as Ch. 5's "
                                               "2033 box; tab:primcost H_I multiplicities, U_phi twice (phi_mult)")
    static_src_qubits: Tagged = Stated(4, f"{TEX}:97", "two static color sources, one qutrit each in two qubits")
    lightlike_dims: Tagged = Stated((3, 3, 5), f"{TEX}:97", "pure-gauge 2O on 3x3x5: the 5-site (1 fm) axis delays the "
                                                            "self-wake to L ~ 0.5 fm")
    lightlike_src_qubits: Tagged = Stated(2, f"{TEX}:97", "two SU(2) fundamental sources, one qubit each")
    qhat_GeV2_fm: Tagged = Assumed(2.0, "pure-glue q-hat near 1.5 T_c (simplobs.json; range 1-3)", src=f"{TEX}:97")
    kappa_over_T3: Tagged = Assumed(2.5, "quenched heavy-quark kappa/T^3 (simplobs.json; range 1.5-3.5)", src=f"{TEX}:97")
    dipole_T_GeV: Tagged = Assumed(0.44, "aT ~ 0.45 at a = 0.2 fm", src=f"{TEX}:97")
    dipole_P0: Tagged = Assumed(0.9, "singlet survival after the string build, P_d (simplobs.json)", src=f"{TEX}:97")
    static_colour_factor: Tagged = Assumed(1 / 3, "Gamma = kappa r^2 / 3 at small r (the O(1) factor, simplobs.json)",
                                           src=f"{TEX}:97")
    thermal_variance_factor: Tagged = Stated((1, 1), f"{TEX}:98", "binomial counts on the measured P_s. q1 (ruling B, "
                                             "2026-10-05): the medium is a quench-prepared pure state at the target energy "
                                             "density (ETH / typicality), so the counts carry no estimator penalty; the "
                                             "E-rho-OQ factor of r26 G8 is retired with the method (ruling A)")
    # ---- q1 (author rulings A and B, 2026-10-05): quench / ETH thermal route for the dipole media ----------------
    # Pricing basis: scratchpad/quench/cost/price2.py on the live model (same dipole_price, R-TOL eps recomputed for the
    # longer circuit). Checks (scratchpad/quench/res): exact evolution on Z2 up to 5x3, 2T up to 2x2, fermion toy N = 16.
    quench_ramp_a: Tagged = Stated(1.0, f"{TEX}:98", "coupling ramp from the electric vacuum over one lattice unit a "
                                                     "(priced as a floor): 20 steps at 0.01 fm/c")
    quench_c_range: Tagged = Cited((2.0, 2 * math.pi), "arxiv_2011_09814", "thermalization time t_th = c / T before the "
                                   "sources enter: c = 2 from the exact small-lattice checks above, c = 2 pi from "
                                   "Hayata-Hidaka's 3+1D SU(2) small-lattice study; 90-282 steps at the 0.01 fm/c "
                                   "grid, not the row's rule step (record t_shot_prep_at_rule_step, unpriced)")
    eth_assumed: Tagged = Assumed("ETH", "eigenstate thermalization / dynamical typicality hold for Sigma(216x3) H_I and 2O "
                                         "at aT ~ 0.45 (Srednicki cond-mat/9403051, Rigol 0708.1324, D'Alessio 1509.06411, "
                                         "Bartsch-Gemmer 0902.0927; SU(2) LGT checks Yao 2303.14264, Ebner 2308.16202), so "
                                         "late-time correlators are microcanonical up to O(1/V); checked here only on the "
                                         "small lattices above. The coarse-step thermalization lever (2-5x) is unpriced "
                                         "(NEEDS_AUTHOR)", src=f"{TEX}:98")
    # ---- r27 (chapter-open pass 2026-10-05: J1 V/PS + calibration, J2 box in the trend, J3 N_h) ----
    lund_kappa_GeV_fm: Tagged = Cited(1.0, "arxiv_2203_11601", "physics/hadronization.tex:21 'kappa ~ 1 GeV/fm'; :31 'the "
                                      "string reaches its maximal extension at t = sqrt(s)/2kappa' (yo-yo of a q-qbar "
                                      "pair). A 3+1D value; the 2+1D model's tension in physical units is set by step 1")
    ff_z_range: Tagged = Stated((0.1, 0.9), f"{TEX}:47", "'over z in [0.1, 0.9]'; the N_z bins split it evenly")
    n_dipole_T: Tagged = Stated(3, f"{TEX}:159", "three temperatures for the scan (the requirements row)")

    # ---- referee G7 (2026-10-04): Trotter error at the chosen Delta t = 0.1 a (trotter_check_2033) -------------
    trotter_eps: Tagged = Assumed(0.1, "G7: target Trotter error of a shot, ||(U - e^{i theta} S^N) rho^{1/2}||_2 in the "
                                       "simulated state, the measure and value of Ch. 9 (app06 'at eps ~ 0.1')", src="app06")
    trotter_g2: Tagged = Assumed(1.0, "G7: KS coupling g^2 for the worst-case bound only; the bound is also minimized "
                                      "over g^2 (trotter_check_2033)")
    trotter_temp_lat: Tagged = Assumed(0.0, "G7: the evolved state is the vacuum plus the injected pair, so the free-field "
                                            "and state-dependent evaluations take T = 0; the pair enters through the "
                                            "dispersion phase of its quanta (trotter_check_2033)")

    trotter_rule: Tagged = Assumed("state", "STEP RULE (2026-10-04, Claude's decision on the G7 item, NEEDS_AUTHOR "
                                            "closed; as Ch. 5 and Ch. 9): the state-dependent second-order estimate at "
                                            "eps = 0.1 sets the step, the exact free-field error is the check, the "
                                            "worst-case bound is quoted once. The 2033 window evolves the ramp-prepared "
                                            "vacuum, an eigenstate of H, whose connected error does not accumulate (it is "
                                            "the bounded vacuum piece the exact free-field check gives, 0.04), so the rule "
                                            "is applied to the injected excitation: its dispersion phase. The 4.65 GeV "
                                            "point at 50 a takes 575 steps (0.132 -> 0.0997); the other points keep "
                                            "Delta t = 0.1 a. The thermal dipole rows take the literal estimate (as Ch. 5), "
                                            "which in a mixed state is an upper estimate of the part that grows with t: "
                                            "each shot max(t / 0.01 fm/c, N_sd) steps. 'grid' restores the chosen steps")

    # ---- r25 R9/R10: T-depth per shot (factory.json, Ch. 6 runs, 2026-10-02) ------------------------
    workspace_links_2028: Tagged = Stated(2, f"{TEX}:88", "R10: two hop links at once, each with one k = 9 phasing "
                                                          "group's 7 ancilla and 4 repeat-until-success ancilla (its 4 "
                                                          "rotations in parallel): 22 ancilla, 230 LQ")
    hwp_adder_depth_2028: Tagged = Assumed((4, 7), "Toffoli layers of the k = 9 Hamming-weight adder, log tree to serial "
                                                   "(factory.json Ch. 6 2028)", src="factory.json")
    depth_link_const_2033: Tagged = Assumed((1907, 3792), "Sigma(72x3) hop, Toffoli layers per link per step, (P) schedule: "
                                                          "2 x 830 squish (log tree) + 90 flags + 112 + 40 + 5 (low) to "
                                                          "2 x 1219 + 892 + 326 + 120 + 16 (high); factory.json Ch. 6 2033",
                                            src="factory.json")
    depth_link_rot_2033: Tagged = Assumed((205, 613), "rotation layers per link per step (each a = T per rotation deep), "
                                                      "same schedule", src="factory.json")
    depth_step_const_2033: Tagged = Assumed((874, 5495), "per-step Toffoli layers outside the hop: electric 152, magnetic "
                                                         "718 to 5,336, mass 4 to 7", src="factory.json")
    depth_step_rot_2033: Tagged = Assumed((1155, 1170), "per-step rotation layers outside the hop: electric 1,152 (serial "
                                                        "rotations, 9-way U_phi needs 576 RUS ancilla), magnetic 2 to 14, "
                                                        "mass 1 to 4", src="factory.json")

    # ---- NOT USED IN ANY HEADLINE: alternatives the author may rule on ------------------
    draft_hop_S72x3_t_const: Tagged = Uncited(46243.6, "UNPUBLISHED DRAFT's printed staggered C^S for Sigma(72x3), per link per field "
                                                       "(FermionPrimitives_unpub, publication_draft/section_resources.tex:40). Since r17 "
                                                       "the hop is built from the draft's gate tables (groups.hop_link_cost); this "
                                                       "printed constant feeds only the draft_* record intermediates")
    draft_hop_S72x3_t_log: Tagged = Uncited(690.2, "same draft, same line: coefficient of log2(1/eps)")

    def __post_init__(self):
        may_be_zero = {"ramp_steps_2033", "n_pauli_mass", "trotter_temp_lat"}
        for name in ("window_over_a", "pair_clustering", "n_meson", "dipole_times_fm", "thermal_variance_factor"):
            v = getattr(self, name)
            if not (v.is_range and 0 < v.lo <= v.hi):
                raise ValueError(f"{name} must be a positive (lo, hi) range")
        for name in ("delta_first_2033", "trend_size", "dipole_P0"):
            if not (0 < getattr(self, name).value < 1):
                raise ValueError(f"{name} must lie in (0, 1)")
        if int(self.workspace_links_2028.value) < 1:
            raise ValueError("workspace_links_2028: at least one hop link at a time")
        for f in fields(self):
            v = getattr(self, f.name)
            if not isinstance(v, Tagged):
                raise TypeError(f"{f.name} must be Tagged")
            if f.name == "ramp_steps_2028":
                if isinstance(v.value, bool) or not isinstance(v.value, (int, float)) or v.value < 0:
                    raise ValueError("ramp_steps_2028 is a non-negative number of steps (round-D ruling 1: 1)")
                continue
            if v.is_range and v.lo > v.hi:
                raise ValueError(f"{f.name}: lo > hi")
            if isinstance(v.lo, (int, float)) and not isinstance(v.lo, bool):
                if v.lo < 0 or (v.lo == 0 and f.name not in may_be_zero):
                    raise ValueError(f"{f.name} must be positive")
        if self.synthesis_errors.value not in ("incoherent", "coherent"):
            raise ValueError("synthesis_errors must be 'incoherent' or 'coherent' (common.eps_per_rotation)")
        if not all(x > 0 for x in self.dt_multi_over_a.value):
            raise ValueError("dt_multi_over_a: positive step sizes")
        if self.link_width_source.value not in ("chapter", "compiled"):
            raise ValueError("link_width_source must be 'chapter' or 'compiled'")
        if self.hamiltonian.value not in PRIMCOST:
            raise ValueError(f"hamiltonian must be one of {tuple(PRIMCOST)}")
        if self.synthesis_model.value != "rus":
            raise ValueError("E20: the gauge tables of groups.py are priced at the full fit ('rus', PrimitiveCost.t); the "
                             "papers' slope-only price is the legacy record PrimitiveCost.t_papers, not this switch")
        if self.fermion_synthesis.value != "rus":
            raise ValueError("r17 rule: fermion-term rotations at the full fit 1.15 log2(1/eps) + 9.2 ('rus'), never slope-only")
        if self.hop_share.value not in ("link", "draft"):
            raise ValueError("hop_share must be 'link' (E21 headline) or 'draft' (sensitivity); groups.fermion_hop_counts")
        if int(self.n_pauli_hop_diagonal.value) != 8:
            raise ValueError("the Z3 hop structure in groups.fermion_hop_counts is 8 strings per copy (RHOP_N_PAULI_DIAGONAL)")
        if self.toffoli_convention.value != "textbook":
            raise ValueError("ruling R5: 7 T per Toffoli ('textbook'); the groups.py tables are in that convention")
        for name in ("group_2028", "group_2033", "group_dipole_static"):
            g = GROUPS.get(getattr(self, name).value)
            if g is None:
                raise ValueError(f"{name}: unknown gauge group")
            need = set(MAGNETIC) | {"U_phi"}
            if not need <= set(g.primitives) or not ({"U_FFT", "U_F"} & set(g.primitives)):
                raise ValueError(f"{name}: groups.py has no complete primitive table (U_inv, U_mul, U_Tr, U_phi, a Fourier "
                                 "transform) for this group")
        for name in ("dim_2028", "L_2028", "n_anc_2028", "n_c", "n_stag", "dim_2033", "L_perp",
                     "L_par_uncut", "L_par_priced", "n_anc_2033", "L_s_3d", "n_z", "n_tens", "n_boost",
                     "n_spacings", "n_source_momenta", "n_qhat_T", "n_species_channels", "n_instances_utility",
                     "n_trot_2028", "ramp_steps_2028", "twosteps_n_trot_2028", "n_pauli_hop_diagonal", "n_pauli_hop_moved", "moves_per_field",
                     "n_pauli_mass", "insertion_links"):
            if not float(getattr(self, name).value).is_integer():
                raise ValueError(f"{name} must be an integer")
        if not all(float(x).is_integer() for x in self.L_par_cut.value):
            raise ValueError("L_par_cut must be integers")
        if not (1 <= self.L_par_cut.lo <= self.L_par_cut.hi <= self.L_par_uncut.value):
            raise ValueError("cut L_par must lie inside the uncut extent")
        if not (self.L_par_cut.lo <= self.L_par_priced.value <= self.L_par_cut.hi):
            raise ValueError("L_par_priced must lie inside the cut range")
        if int(self.dim_2028.value) not in (2, 3):
            raise ValueError("dim_2028 must be 2 or 3 (1+1D is excluded by the chapter, app10:106)")
        if int(self.dim_2033.value) != 2:
            raise ValueError("dim_2033 must be 2: the 2033 lattice is L_perp x L_par. The 3+1D lattice is the separate "
                             "extension (L_s_3d, era 'codesign')")
        if not (0 < self.hadamard_amplitude.lo <= self.hadamard_amplitude.hi <= 1):
            raise ValueError("Hadamard amplitude must lie in (0, 1]")
        for name in ("eps_stat", "n_s", "utility_share", "syst_volume", "syst_subgroup", "clean_shot_faults"):
            if not (0 < getattr(self, name).value <= 1):
                raise ValueError(f"{name} must lie in (0, 1]")
        for name in ("delta_rel_2028", "delta_rel_2033", "syst_discretization"):
            v = getattr(self, name)
            if not (v.is_range and 0 < v.lo <= v.hi <= 1):
                raise ValueError(f"{name} must be a (lo, hi) range in (0, 1]")
        for name in ("eps_rot", "eps_syn", "rfi_eps_l", "eps_rot_before_round_e"):
            if not (0 < getattr(self, name).value < 1):
                raise ValueError(f"{name} must lie in (0, 1)")
        if int(self.n_trot_2028.value) < 1:
            raise ValueError("n_trot_2028: at least one evolution step")
        if not (1e-9 <= self.t_gate_s.value <= 1e-3):
            raise ValueError("gate time must be 1 ns - 1 ms")


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #

def link_width(a: Assumptions, group: str) -> int:
    """Qubits per link: the compiled register ('compiled', ruling R1) or the pre-ruling print ('chapter')."""
    g = GROUPS[group]
    if a.link_width_source.value == "chapter":
        return int(g.link_qubits_chapter.lo)
    return g.link_qubits


def n_plaquettes(V: int, D: int) -> int:
    """Periodic lattice: V plaquettes per pair of directions, D(D-1)/2 pairs. 2D -> V, 3D -> 3V."""
    return V * D * (D - 1) // 2


def register(q_g: int, V: int, n_link: int, n_c: int, n_stag: int, n_anc: int) -> tuple[int, int, int, int]:
    """eq:Nq_collider: (gauge, fermion, ancilla, total)."""
    gauge = q_g * V * n_link
    fermion = n_c * n_stag * V
    return gauge, fermion, n_anc, gauge + fermion + n_anc


def shots_binomial(n_s: float, delta_rel: float) -> float:
    """eq:Nshot_samp_collider: N = 1 / (n_s delta_rel^2)."""
    return 1.0 / (n_s * delta_rel ** 2)


def delta_binomial(n_s: float, n_shots: float) -> float:
    """The relative error N shots buy: (n_s N)^(-1/2) (app10:80)."""
    return (n_s * n_shots) ** -0.5


def shots_hadamard(a_inv2: float, eps: float, n_z: float, n_tens: float, n_boost: float) -> float:
    """eq:Nshot_collider."""
    return a_inv2 * eps ** -2 * n_z * n_tens * n_boost


def quadrature(*xs: float) -> float:
    return math.sqrt(sum(x * x for x in xs))


def fourier_key(group: str) -> str:
    """The Fourier-transform entry ruling R2 prices: the fast transform where groups.py lists one (COMPILED for
    2T, 2O, Sigma(36x3); the paper's stated LOWER BOUND for Sigma(72x3), where none has been constructed), else
    the group's only transform (Z3: the 14-rotation H_3 of arxiv_2408_00075)."""
    return "U_FFT" if "U_FFT" in GROUPS[group].primitives else "U_F"


def electric_floor_per_link(group: str, ham: str, d: int, eps: float, papers: bool = False) -> float:
    """n_F x (Fourier transform of fourier_key) + U_phi, per link per Trotter step. For Sigma(72x3) this is a
    floor (ruling R2); groups.electric_per_link would fall back to the dense U_F there. Rotations at the full fit
    (E20); papers=True gives the papers' slope-only price (legacy record)."""
    g = GROUPS[group]
    k_phi = g.phi_mult.get(ham, 1)                                        # s648: U_phi twice under H_I for Sigma(216x3)
    if papers:
        return (PRIMCOST[ham]["U_F"](d) * g.primitives[fourier_key(group)].t_papers(eps)
                + k_phi * g.primitives["U_phi"].t_papers(eps))
    return PRIMCOST[ham]["U_F"](d) * g.primitives[fourier_key(group)].t(eps) + k_phi * g.primitives["U_phi"].t(eps)


def gauge_step(a: Assumptions, group: str, d: int, ham: str | None = None) -> dict:
    """Ruling R2 structure, E20 price: magnetic and electric T per link per Trotter step from groups.py with every
    rotation at the full fit ('magnetic', 'electric', 'per_link'; the same as *_report). The papers' slope-only
    price is kept as the *_papers record. ham defaults to a.hamiltonian (H_KS); s648: the static dipole row passes
    'I'. U_phi enters phi_mult times (1 except Sigma(216x3) under H_I, 2); a group with no dense U_F (Sigma(216x3))
    has no electric_dense record (None)."""
    ham = a.hamiltonian.value if ham is None else ham
    eps = float(a.eps_rot)
    g = GROUPS[group]
    conv = a.toffoli_convention.value
    m = PRIMCOST[ham]
    fk = fourier_key(group)
    mag_papers = magnetic_per_link(group, ham, d, eps, papers=True)
    el_papers = electric_floor_per_link(group, ham, d, eps, papers=True)
    mag = sum(m[p](d) * g.primitives[p].t_report(eps, conv) for p in MAGNETIC)
    k_phi = g.phi_mult.get(ham, 1)
    el = m["U_F"](d) * g.primitives[fk].t_report(eps, conv) + k_phi * g.primitives["U_phi"].t_report(eps, conv)
    has_dense = "U_F" in g.primitives
    el_dense = (m["U_F"](d) * g.primitives["U_F"].t_report(eps, conv) + k_phi * g.primitives["U_phi"].t_report(eps, conv)
                if has_dense else None)
    # Toffoli and rotation counts per link per step, so the conventions are checked, not asserted
    mult = {p: m[p](d) for p in MAGNETIC}
    mult[fk] = m["U_F"](d)
    mult["U_phi"] = k_phi
    n_tof = sum(k * g.primitives[p].toffoli for p, k in mult.items())
    n_rot = sum(k * g.primitives[p].rot for p, k in mult.items())
    return {"magnetic": mag, "electric": el, "per_link": mag + el, "fourier_key": fk,
            "fourier_status": g.primitives[fk].status, "fft_compiled": g.has_fft,
            "electric_dense": electric_per_link(group, ham, d, eps, fft=False) if has_dense else None,
            "electric_dense_report": el_dense, "hamiltonian": ham,
            "magnetic_report": mag, "electric_report": el, "per_link_report": mag + el,
            "magnetic_papers": mag_papers, "electric_papers": el_papers, "per_link_papers": mag_papers + el_papers,
            "toffoli_per_link": n_tof, "rotations_per_link": n_rot, "multiplicity": mult,
            "t_per_rotation": t_per_rotation(eps, a.synthesis_model.value),
            "t_per_rotation_papers": t_per_rotation(eps, "rus-slope"),       # legacy record (R2)
            "t_per_rotation_fermion": t_per_rotation(eps, a.fermion_synthesis.value),
            "t_per_toffoli": toffoli_t(1, conv)}


def gauge_primitives(a: Assumptions, group: str, d: int, n_links: int, steps: float, tag: str) -> tuple[Primitive, ...]:
    """The per-shot Primitives of the gauge terms: each groups.py primitive with its tab:primcost multiplicity
    x links x steps, t_each = the table's a + (b / 1.15)(1.15 log2(1/eps) + 9.2) (E20, full fit)."""
    gs = gauge_step(a, group, d)
    g = GROUPS[group]
    eps = float(a.eps_rot)
    ham = a.hamiltonian.value
    conv = a.toffoli_convention.value
    out = []
    for p, k in gs["multiplicity"].items():
        pc = g.primitives[p]
        t_each = pc.t(eps)
        if t_each == 0:
            what = "Clifford, 0 T"
        else:
            what = (f"{pc.t_const:g} + {pc.n_rot:g} x (1.15 log2(1/eps) + 9.2) T at eps={eps:.4g} "
                    f"(paper: {pc.t_const:g} + {pc.t_log:g} log2(1/eps))")
        term = "magnetic" if p in MAGNETIC else "electric"
        extra = ""
        if p == "U_FFT" and pc.status is not CircuitStatus.COMPILED:
            extra = "; NOT CONSTRUCTED: the paper's stated lower bound, a FLOOR (ruling R2)"
        out.append(Primitive(f"{p}_{group}_{tag}", k * n_links * steps, t_each, pc.status, pc.src,
                             f"{term} term, {k:g} per link-step (tab:primcost H_{ham}) x {n_links} links x {steps:g} steps; "
                             f"{what}{extra}" + (f"; {pc.note}" if pc.note.startswith(("DERIVED", "STATED")) else "")))
    return tuple(out)


def hwp(a: Assumptions, k: int) -> dict:
    """One Hamming-weight-phasing group of k equal-angle rotations (Ch. 9's rule, app06:108; arxiv_1709_06648,
    arxiv_1902_10673): floor(log2 k)+1 rotations at the full fit and k - w(k) Toffolis at 7 T, measurement-based
    uncompute (no T, no extra qubits). Shared: common.hwp_group."""
    return hwp_group(k, float(a.eps_rot), a.fermion_synthesis.value, a.toffoli_convention.value)


def hop_link_rhop(a: Assumptions, group: str, synthesis: str = "rus-slope") -> dict:
    """LEGACY rule R-HOP (round D to R-TOL), per link per Trotter step, for the record only: groups.rhop_link with
    this chapter's Stated inputs, by default at the slope-only price the chapter printed before r17."""
    n_p = int(a.n_pauli_hop_diagonal.value) if group in DIAGONAL_LINK else int(a.n_pauli_hop_moved.value)
    return rhop_link(group, int(a.n_stag.value), int(a.n_c.value), float(a.eps_rot), n_pauli=n_p,
                     moves_per_field=int(a.moves_per_field.value), synthesis=synthesis,
                     toffoli_convention=a.toffoli_convention.value)


def hop_link(a: Assumptions, group: str) -> dict:
    """The hop per link per Trotter step since r17 (ruling (e), 2026-10-01): groups.hop_link_cost at this circuit's
    R-TOL eps and the full fit. Z3 (diagonal link, no draft entry): the R-HOP structure, 8 strings x HWP(k = N_c
    N_stag), 56 Toffoli + 32 rotations per link. Sigma(72x3): the authors' unpublished color-diagonalized staggered
    hop (FermionPrimitives_unpub) for the N_stag fields on the link, with the pieces the draft leaves out estimated
    in groups.py (frame undo, SU(3) squish uncompute, MBU ladders, HWP phasing). Since r19 (ruling E21 (1)) the color
    squish and parity are held once per link (a.hop_share = 'link'); the Z3 hop has no squish, so share does not move
    it. Returns the hop_link_cost dict plus
    "rhop" (the legacy R-HOP at the slope-only price, the pre-r17 print) and "rhop_full" (R-HOP at the full fit)."""
    n_stag, n_c = int(a.n_stag.value), int(a.n_c.value)
    eps = float(a.eps_rot)
    kw = {"n_c": n_c, "k": n_c * n_stag} if group in DIAGONAL_LINK else {}
    h = hop_link_cost(group, n_stag, eps=eps, synthesis=a.fermion_synthesis.value,
                      toffoli_convention=a.toffoli_convention.value, legacy=False, share=a.hop_share.value, **kw)
    h["rhop"] = hop_link_rhop(a, group, "rus-slope")
    h["rhop_full"] = hop_link_rhop(a, group, "rus")
    return h


def mass_site(a: Assumptions) -> dict:
    """Staggered mass, per site per Trotter step: one HWP group of the N_c N_stag copies per Z string
    (shared helper groups.staggered_mass_site, r17), at the full fit."""
    h = staggered_mass_site(int(a.n_c.value), int(a.n_stag.value), float(a.eps_rot),
                            synthesis=a.fermion_synthesis.value, toffoli_convention=a.toffoli_convention.value)
    n = int(a.n_pauli_mass.value)
    return {"hwp": h, "n_groups": n, "t": n * h["t"]}


def fermion_primitives(a: Assumptions, group: str, n_links: int, V: int, steps: float, tag: str) -> tuple[Primitive, ...]:
    """Per-shot Primitives of the hop (groups.hop_link_cost items) and the mass over `steps` Trotter steps."""
    hl, ms = hop_link(a, group), mass_site(a)
    conv = a.toffoli_convention.value
    hwp_src = "arxiv_1709_06648,arxiv_1902_10673"
    mult = n_links * steps
    out = []
    for it in hl["items"]:
        base = f"hop_{it['name']}_{group}_{tag}"
        why = f"{it['note']} [{it['src']}]"
        if it["toffoli"]:
            out.append(Primitive(f"{base}_toffoli", it["toffoli"] * mult, toffoli_t(1, conv), it["status"], it["src"],
                                 f"{it['toffoli']} Toffolis per link-step x {n_links} links x {steps:g} steps; {why}"))
        if it["t_direct"]:
            out.append(Primitive(f"{base}_t", it["t_direct"] * mult, 1.0, it["status"], it["src"],
                                 f"{it['t_direct']} direct T per link-step; {why}"))
        if it["n_rot"]:
            out.append(Primitive(f"{base}_rotations", it["n_rot"] * mult, hl["t_rot"], it["status"], it["src"],
                                 f"{it['n_rot']} rotations per link-step at {hl['t_rot']:.4f} T (full fit); {why}"))
    if ms["n_groups"]:
        m_groups = ms["n_groups"] * V * steps
        mh = ms["hwp"]
        out.append(Primitive(f"mass_hwp_rotations_{tag}", m_groups * mh["n_rot"], mh["t_rot"], CircuitStatus.SCALING,
                             f"FermionPrimitives_unpub; {hwp_src}",
                             "staggered mass: one Z per copy, one k-group per site (groups.staggered_mass_site)"))
        out.append(Primitive(f"mass_hwp_toffolis_{tag}", m_groups * mh["n_toffoli"], toffoli_t(1, conv),
                             CircuitStatus.COMPILED, hwp_src, "staggered mass, Hamming-weight adders"))
    return tuple(out)


def synth_rotations_per_step(a: Assumptions, group: str, n_links: int, V: int) -> float:
    """Synthesized rotations in one Trotter step (gauge + hop + mass), the N of the randomized-synthesis budget
    N eps^2 (R-TOL). Exact T gates (the 4 of the Z3 U_FFT, the hop diagonalizers' 2 each) and Toffolis are not
    rotations. Independent of eps."""
    gs = gauge_step(a, group, int(round(n_links / V)))
    hl, ms = hop_link(a, group), mass_site(a)
    return n_links * (gs["rotations_per_link"] + hl["n_rot"]) + V * ms["n_groups"] * ms["hwp"]["n_rot"]


def mass_site_slope_only(a: Assumptions) -> dict:
    """RECORD: the staggered mass at the retired slope-only price (the pre-r17 print)."""
    h = staggered_mass_site(int(a.n_c.value), int(a.n_stag.value), float(a.eps_rot), synthesis="rus-slope",
                            toffoli_convention=a.toffoli_convention.value)
    return {"hwp": h, "n_groups": int(a.n_pauli_mass.value), "t": int(a.n_pauli_mass.value) * h["t"]}


def _step_slope_only(a: Assumptions, group: str, d: int, n_links: int, V: int) -> float:
    """RECORD: one Trotter step as the chapter priced it before r17 (gauge tables at the papers' slope-only price,
    legacy R-HOP hop, slope-only mass), at a's eps."""
    return (n_links * (gauge_step(a, group, d)["per_link_papers"] + hop_link_rhop(a, group)["t"])
            + V * mass_site_slope_only(a)["t"])


def _n_rot_step_pre_r17(a: Assumptions, group: str, n_links: int, V: int) -> float:
    """RECORD: synthesized rotations per step under the retired R-HOP (U_mul moves + HWP strings), the pre-r17 N."""
    gs = gauge_step(a, group, int(round(n_links / V)))
    hr, ms = hop_link_rhop(a, group), mass_site(a)
    per_link = gs["rotations_per_link"] + hr["n_pauli"] * hr["hwp"]["n_rot"] + hr["n_moves"] * GROUPS[group].primitives["U_mul"].rot
    return n_links * per_link + V * ms["n_groups"] * ms["hwp"]["n_rot"]


def rtol_eps(a: Assumptions, n_rot: float) -> float:
    """Ruling R-TOL: the circuit's own rotation tolerance from its per-shot rotation count, holding the total
    synthesis error to eps_syn. Randomized synthesis (the chapter's): common.eps_rot_for = sqrt(eps_syn / N).
    n_rot is fixed by the circuit (synth_rotations_per_step does not read eps), so no iteration."""
    if a.synthesis_errors.value == "incoherent":
        return eps_rot_for(n_rot, float(a.eps_syn))
    return eps_per_rotation(float(a.eps_syn), n_rot, a.synthesis_errors.value)


def _at_eps(a: Assumptions, eps: float, n_rot: float) -> Assumptions:
    """The Assumptions every pricing helper reads, at this circuit's R-TOL tolerance."""
    return _dc_replace(a, eps_rot=Stated(eps, f"{TEX}:88", f"R-TOL: eps_rot_for(N_rot = {n_rot:.6g}, "
                                                           f"eps_syn = {float(a.eps_syn):g})"))


def synthesis_error(a: Assumptions, n_rot: float) -> float:
    """Per-shot synthesis error of n_rot rotations at eps_rot: N eps^2 (randomized) or N eps (deterministic)."""
    e = float(a.eps_rot)
    return n_rot * e * e if a.synthesis_errors.value == "incoherent" else n_rot * e


def insertion_t(a: Assumptions, group: str) -> float:
    """One controlled insertion of the bilocal current (readout (b)), derived here: 2 (L - 1) U_mul for the
    L-link Wilson line and its uncompute, plus the moves of one field, plus a controlled bilinear."""
    eps = float(a.eps_rot)
    L = int(a.insertion_links.value)
    n_umul = 2 * (L - 1) + int(a.moves_per_field.value)
    t_umul = GROUPS[group].primitives["U_mul"].t(eps)                     # gauge primitive, no rotations
    bil = 2 * (float(a.insertion_bilinear_rotations) * t_per_rotation(eps, a.fermion_synthesis.value)
               + toffoli_t(float(a.insertion_bilinear_toffolis), a.toffoli_convention.value))
    return n_umul * t_umul + bil


def ch9_rules_step(a: Assumptions, V: int, D: int) -> dict:
    """Ch. 9's 2028 Z3 counting rules (app06:108; estimates.ch09_chiral_gauge) applied to THIS lattice.

    NOT the chapter's number and not used in any headline. Before round D the chapter carried 5e3 T/step of hopping
    'based on the Ch. 9 Z3 accounting'; this is what Ch. 9's PRE-R-HOP accounting gives (2 strings per hop, no Z3
    dressing, 30 T per rotation). Ch. 9 no longer counts this way: since R-HOP (2026-09-29) its spatial Z3 hop is
    groups.rhop_link, 8 strings at 15.28 T per hop rotation, the rule hop_link applies here. Only its link-free
    fifth-direction hops keep 2 strings (c9.n_pauli_per_hop), which is what this record reads.
      hopping : one same-angle group per (link, Pauli string) over the N_c N_stag copies that share
                the link: V D x 2 groups of k = 9
      electric: 3 diagonal Z-strings per two-qubit link, each grouped across the V D links (Ch. 9 works
                in the electric basis; the chapter's electric term is priced from groups.py instead)
    Each group of k costs floor(log2 k)+1 synthesized rotations and k - w(k) Toffolis (Hamming-weight
    phasing, measurement-based uncompute, allowed by ruling).
    """
    c9 = ch9.Assumptions()
    n_links = V * D
    k_hop = int(a.n_c.value) * int(a.n_stag.value)
    ps_hop = int(c9.n_pauli_per_hop.value)
    ps_el = int(c9.n_pauli_per_electric.value)
    t_rot = float(c9.t_per_rot.value)
    conv = c9.toffoli_convention.value
    groups = (("hopping", k_hop, n_links * ps_hop), ("electric", n_links, ps_el))
    out = {"groups": groups, "n_rot": 0, "n_synth": 0, "n_toffoli": 0, "max_ancilla": 0,
           "t_per_rot": t_rot, "toffoli_convention": conv}
    for name, k, n in groups:
        out["n_rot"] += k * n
        out["n_synth"] += n * ch9.hwp_synth_rotations(k)
        out["n_toffoli"] += n * ch9.hwp_toffolis(k)
        out["max_ancilla"] = max(out["max_ancilla"], ch9.hwp_ancilla(k))
        out[f"t_{name}"] = n * (ch9.hwp_synth_rotations(k) * t_rot + toffoli_t(ch9.hwp_toffolis(k), conv))
    out["ancilla_hop_group"] = ch9.hwp_ancilla(k_hop)
    out["t_step"] = out["n_synth"] * t_rot + toffoli_t(out["n_toffoli"], conv)
    return out


# --------------------------------------------------------------------------- #
# r25 helpers: T-depth (R9/R10), ratio shots (R2), corrected walls, dipole rows (R7)
# --------------------------------------------------------------------------- #

def depth_2028(a: Assumptions, t_rot: float, n_links: int, V: int, n_steps: float, n_conc: int, est: int) -> float:
    """T-depth per 2028 shot (non-Clifford layers; a Toffoli is one layer with MBU, a rotation t_rot layers), the
    factory.json Ch. 6 2028 structure with n_conc hop links at once, each holding a k = 9 phasing group (7 ancilla)
    and 4 RUS ancilla, so a group is one adder (4..7 layers) plus one rotation deep. est = 0 (low) / 1 (high).
      electric : U_FFT (4 T + 2 rotations), U_phi, U_FFT^dag on every link; rotations bounded by the RUS-capable
                 workspace, ceil(n_links / (11 n_conc)) batches of 2 (4 + 2 a) + a
      magnetic : 2 checkerboard classes, 6 (low) to 12 (high) U_mul of 4 Toffolis plus one U_Tr rotation
      hop      : ceil(n_links / n_conc) batches x 8 strings x (adder + a)
      mass     : ceil(V / n_conc) batches x (adder + a)"""
    adder = (a.hwp_adder_depth_2028.lo, a.hwp_adder_depth_2028.hi)[est]
    anc = 11 * n_conc
    e = math.ceil(n_links / anc) * (2 * (4 + 2 * t_rot) + t_rot)
    mag = 2 * ((6, 12)[est] * 4 + t_rot)
    hop = math.ceil(n_links / n_conc) * int(a.n_pauli_hop_diagonal.value) * (adder + t_rot)
    mass = math.ceil(V / n_conc) * (adder + t_rot)
    return n_steps * (e + mag + hop + mass)


def depth_2033_step(a: Assumptions, t_rot: float, n_links: int, est: int) -> float:
    """T-depth of one 2033 Sigma(72x3) Trotter step on the priced register, factory.json schedule (P): the hop
    links run one at a time (each holds ~92 color-angle flags, ~110 LQ, against 100 ancilla), the electric
    rotations in series. est = 0 (low: log-tree squish) / 1 (high: serial ladders)."""
    pick = lambda v: (v.lo, v.hi)[est]
    return (pick(a.depth_step_const_2033) + pick(a.depth_step_rot_2033) * t_rot
            + n_links * (pick(a.depth_link_const_2033) + pick(a.depth_link_rot_2033) * t_rot))


def wall_corrected(shots: float, n_t: float, d_t: float, t_gate: float, t0: float) -> float:
    """One machine, report baseline of ~10 factories: each shot takes max(N_T t_gate, D_T t_r) + t0 (t_r = 10 us)."""
    return shots * (max(n_t * t_gate, d_t * REACTION_TIME_S) + t0)


def ratio_shots(a: Assumptions, delta: float) -> tuple[float, float]:
    """R2: shots for a species ratio at relative error delta, conserved pairs counted:
    (c / n_s + 1 / n_M) / delta^2, low end c = 1 and n_M = 1, high end c = 2 and n_M = 0.5."""
    n_s = float(a.n_s)
    return ((a.pair_clustering.lo / n_s + 1 / a.n_meson.hi) / delta ** 2,
            (a.pair_clustering.hi / n_s + 1 / a.n_meson.lo) / delta ** 2)


def trend_delta(a: Assumptions) -> float:
    """R2: per-point error at which the 3-momentum trend shows at trend_sigma: trend / (N sqrt(2))."""
    return float(a.trend_size) / (float(a.trend_sigma) * math.sqrt(2))


def pure_gauge_step(a: Assumptions, group: str, d: int, ham: str | None = None) -> dict:
    """Per-link T and rotations per step for a pure-gauge circuit. Groups with U_phi go through gauge_step; 2O (U_phi
    folded into its compiled FFT) through groups.magnetic_per_link + electric_per_link (fft=True). ham defaults to
    a.hamiltonian (H_KS)."""
    ham = a.hamiltonian.value if ham is None else ham
    if "U_phi" in GROUPS[group].primitives:
        gs = gauge_step(a, group, d, ham)
        return {"per_link": gs["per_link"], "rotations_per_link": gs["rotations_per_link"]}
    eps, m, g = float(a.eps_rot), PRIMCOST[ham], GROUPS[group]
    mult = {p: m[p](d) for p in MAGNETIC}
    mult["U_FFT"] = m["U_F"](d)
    return {"per_link": magnetic_per_link(group, ham, d, eps)
                        + electric_per_link(group, ham, d, eps, fft=True),
            "rotations_per_link": sum(k * g.primitives[p].rot for p, k in mult.items())}


def dipole_price(a: Assumptions, group: str, dims: tuple, n_steps: int, n_src_actions: int, src_qubits: int,
                 ham: str | None = None) -> dict:
    """R7: one dipole shot on a pure-gauge 3+1D lattice, its own R-TOL eps; source moves are U_mul actions."""
    d = len(dims)
    V = math.prod(dims)
    n_links = V * d
    n_rot = n_steps * n_links * pure_gauge_step(a, group, d, ham)["rotations_per_link"]
    eps = rtol_eps(a, n_rot)
    ax = _at_eps(a, eps, n_rot)
    t_step = n_links * pure_gauge_step(ax, group, d, ham)["per_link"]
    t_src = n_src_actions * GROUPS[group].primitives["U_mul"].t(eps)
    lq = GROUPS[group].link_qubits * n_links + int(a.dipole_n_anc.value) + src_qubits
    return {"V": V, "links": n_links, "lq": lq, "eps": eps, "t_step": t_step, "t_shot": n_steps * t_step + t_src,
            "n_rot": n_rot, "group": group, "hamiltonian": ham or a.hamiltonian.value}


def two_point_log_shots(p1: float, p2: float, delta: float) -> float:
    """Shots (both points) for the log slope ln(P1/P2) at relative error delta, binomial at each point."""
    var = (1 - p1) / p1 + (1 - p2) / p2
    return 2 * var / (delta * math.log(p1 / p2)) ** 2


# --------------------------------------------------------------------------- #
# 2028: species-count pipeline, 2+1D Z3, 4^2, 3 staggered fields
# --------------------------------------------------------------------------- #

def _model_2028(a: Assumptions) -> Result:
    gk = a.group_2028.value
    g = GROUPS[gk]
    q_g = link_width(a, gk)                                               # 2
    D, L = int(a.dim_2028.value), int(a.L_2028.value)
    V = L ** D                                                            # 16
    n_links = V * D                                                       # 32
    n_plaq = n_plaquettes(V, D)                                           # 16
    n_c, n_stag = int(a.n_c.value), int(a.n_stag.value)
    gauge, fermion, anc_r22, lq_r22 = register(q_g, V, D, n_c, n_stag, int(a.n_anc_2028.value))   # 64, 144, 8, 216 (r22)

    cap = float(a.rfi_t_2028)
    n_evo = float(a.n_trot_2028)                                          # 2 (r18 ruling: 1.32x accepted)
    n_ramp = float(a.ramp_steps_2028)                                     # 1 (round-D ruling 1)
    n_steps = n_evo + n_ramp                                              # 3

    # ruling R-TOL: the circuit's rotation count (independent of eps) sets its tolerance
    n_rot_step = synth_rotations_per_step(a, gk, n_links, V)              # 1264
    n_rot_shot = n_steps * n_rot_step                                     # 3792
    a_round_e = a                                                         # eps = 1e-3, round E: record only
    eps = rtol_eps(a, n_rot_shot)                                         # sqrt(1e-2 / 3792) = 1.624e-3
    a = _at_eps(a, eps, n_rot_shot)
    synth_err = synthesis_error(a, n_rot_shot)                            # 1e-2 = eps_syn

    gs = gauge_step(a, gk, D)
    hl, ms = hop_link(a, gk), mass_site(a)
    t_gauge_step = n_links * gs["per_link"]                               # 32 x 285.2 = 9126.7
    t_hop_step = n_links * hl["t"]                                        # 32 x 1027.4 = 32876.8
    t_mass_step = V * ms["t"]                                             # 16 x 128.4 = 2054.8
    t_step = t_gauge_step + t_hop_step + t_mass_step                      # 44,058.3
    t_evo = n_evo * t_step
    t_ramp = n_ramp * t_step
    t_total = t_evo + t_ramp                                              # 132,174.9 -> '1.3e5', 1.32x (accepted)
    steps_in_cap = math.floor(cap / t_step)                               # 2 (2.27): the third step overshoots

    def _shot_at(n):
        """The 2028 shot with n steps in all: its own circuit, its own R-TOL eps."""
        ax = _at_eps(a_round_e, rtol_eps(a_round_e, n * n_rot_step), n * n_rot_step)
        return n * (n_links * (gauge_step(ax, gk, D)["per_link"] + hop_link(ax, gk)["t"]) + V * mass_site(ax)["t"])

    t_total_two = _shot_at(2)                                             # 87,266.3 (0.873x), 1 ramp + 1 evolution
    t_total_four = _shot_at(4)                                            # 1.78e5 (1.78x), record
    # the r17 print: 1 ramp + 1 evolution, gauge tables at the papers' slope-only price, hop and mass at the full fit
    a_r17 = _at_eps(a_round_e, rtol_eps(a_round_e, 2 * n_rot_step), 2 * n_rot_step)
    t_total_r17 = 2 * (n_links * (gauge_step(a_r17, gk, D)["per_link_papers"] + hop_link(a_r17, gk)["t"])
                       + V * mass_site(a_r17)["t"])                       # 84,027.9 (0.84x)
    # the pre-r17 print (R-TOL, slope-only rotations, 1 ramp + 2 evolution), for the record
    n_pre = 3
    a_pre = _at_eps(a_round_e, rtol_eps(a_round_e, n_pre * n_rot_step), n_pre * n_rot_step)
    t_step_pre = _step_slope_only(a_pre, gk, D, n_links, V)               # 32,429.5
    t_total_pre = n_pre * t_step_pre                                      # 97,288.5 (0.973x)
    # round E's fixed eps = 1e-3 (slope-only, 3 steps), for the record
    t_step_round_e = _step_slope_only(a_round_e, gk, D, n_links, V)       # 33,446.26
    t_total_round_e = n_pre * t_step_round_e                              # 100,338.8 (1.003x)

    c9 = ch9_rules_step(a, V, D)
    gap = float(a.t_step_gap_target)

    # records and alternatives that are not the headline
    ham = a.hamiltonian.value
    t_rot = gs["t_per_rotation"]
    eps_old = float(a.eps_rot_before_round_e)
    ft_new, ft_old = g.primitives["U_FFT"], g.primitives["U_F"]
    # the transpiler U_F (14 rotations) in place of U_FFT, at this eps: what round E removed
    el_transpiler = PRIMCOST[ham]["U_F"](D) * ft_old.t(eps) + g.primitives["U_phi"].t(eps)
    t_step_transpiler = n_links * (gs["magnetic"] + el_transpiler) + t_hop_step + t_mass_step           # 4.73e4
    # the whole step at the superseded eps = 1e-4 with U_FFT (the E1 change alone)
    a_old = _dc_replace(a, eps_rot=Stated(eps_old, "record"))
    gs_old, hl_old, ms_old = gauge_step(a_old, gk, D), hop_link(a_old, gk), mass_site(a_old)
    t_step_old_eps_fft = n_links * (gs_old["per_link"] + hl_old["t"]) + V * ms_old["t"]                 # full fit at 1e-4
    # a plaquette written for an abelian group: fold three links into the fourth, phase, unfold (6 U_mul + U_Tr)
    plaq_abelian = 6 * g.primitives["U_mul"].t(eps) + g.primitives["U_Tr"].t(eps)
    t_gauge_step_report = n_links * gs["per_link_report"]                 # gauge step at the full fit, record
    # the retired U_mul-style price applied to Z3, at the full fit (record of the R-HOP comparison)
    umul_style_z3 = (int(a.moves_per_field.value) * n_stag * g.primitives["U_mul"].t(eps)
                     + int(a.n_pauli_hop_moved.value) * hl["rhop_full"]["hwp"]["t"])
    # Ch. 9's current rule (R-HOP, groups.rhop_link with Ch. 9's hop_synthesis) at this chapter's k = 9, at CH. 9's eps
    c9a = ch9.Assumptions()
    c9_now = rhop_link(gk, n_stag=n_stag, n_c=n_c, eps=float(c9a.eps_rot.value), k=n_c * n_stag,
                       synthesis=c9a.hop_synthesis.value, toffoli_convention="textbook")
    # Ch. 9's pre-R-HOP accounting with the Z3 dressing: 8 strings instead of 2, at Ch. 9's 30 T per rotation
    hr = hl["rhop"]["hwp"]
    c9_dressed_hop = n_links * hl["rhop"]["n_pauli"] * (hr["n_rot"] * c9["t_per_rot"] + toffoli_t(hr["n_toffoli"], "textbook"))

    dt_over_a = float(a.t_had_over_a) / float(a.n_trot_2033)             # 0.1, the chapter's step (2033 box)
    dts = tuple(float(x) for x in a.dt_multi_over_a.value)                 # (0.1, 0.2, 0.3), ruling E4
    t_reached = tuple(n_evo * d for d in dts)                             # (0.2, 0.4, 0.6) a
    t_sampled = tuple(sorted({round(k * d, 12) for d in dts for k in range(1, int(n_evo) + 1)}))   # 0.1 ... 0.6
    t_overlap = tuple(sorted({round(k * d, 12) for d in dts for k in range(1, int(n_evo) + 1)
                              if sum(any(abs(j * d2 - k * d) < 1e-9 for j in range(1, int(n_evo) + 1)) for d2 in dts) >= 2}))
    t_gauge_step_round_d = n_links * (gs_old["magnetic_papers"] + PRIMCOST[ham]["U_F"](D) * ft_old.t_papers(eps_old)
                                      + g.primitives["U_phi"].t_papers(eps_old))   # 19801, the round-D gauge step (record)
    retired_step = t_gauge_step_round_d + float(a.retired_hop_2028)       # 2.48e4

    shots = float(a.shots_2028)                                           # 1e3, Stated
    n_s = float(a.n_s)
    delta_at_shots = delta_binomial(n_s, shots)                           # 0.10
    shots_for_delta = (shots_binomial(n_s, a.delta_rel_2028.hi), shots_binomial(n_s, a.delta_rel_2028.lo))   # 1000, 1000
    gate = float(a.t_gate_s)
    t0 = float(a.shot_overhead_s)                                         # ~0.1 ms per shot (Ch. 1 convention)
    wall_shot = t_total * gate + t0                                       # 0.1323 s (0.1322 of T-gates + 0.0001)
    wall = shots * wall_shot                                              # 132.3 s, serial on one machine
    eps_l_req = float(a.clean_shot_faults) / t_total                      # 7.57e-7 -> '8e-7'

    # R10 (r25): the register with two hop links at once. Each link holds one k = 9 phasing group (7 ancilla) and
    # 4 RUS ancilla, so a group's 4 rotations run together: 22 ancilla, 230 LQ (was 8 and 216, rotations serial)
    n_conc = int(a.workspace_links_2028.value)                            # 2
    hwp9 = hl["rhop_full"]["hwp"]
    anc = n_conc * (hwp9["ancilla"] + hwp9["n_rot"])                      # 2 x (7 + 4) = 22
    lq = gauge + fermion + anc                                            # 230
    t_rot_f = gs["t_per_rotation_fermion"]                                # 19.856
    d_shot = (depth_2028(a, t_rot_f, n_links, V, n_steps, n_conc, 0),
              depth_2028(a, t_rot_f, n_links, V, n_steps, n_conc, 1))     # (1.06e4, 1.20e4)
    # the r22 register (8 ancilla: one phasing group, ONE RUS ancilla, so every rotation in series), factory.json:
    # D = N_rot a + exact T + half the Toffolis
    exact_t = ft_new.t_gates * PRIMCOST[ham]["U_F"](D) * n_links * n_steps                 # 4 x 2 x 32 x 3 = 768
    n_tof_shot = (t_total - n_rot_shot * t_rot_f - exact_t) / toffoli_t(1, a.toffoli_convention.value)
    d_shot_r22 = n_rot_shot * t_rot_f + exact_t + 0.5 * n_tof_shot        # 8.0e4
    dx = depth_exports(t_total, d_shot, shots, gate, t0)
    wall_c = wall_corrected(shots, t_total, d_shot[1], gate, t0)          # = wall: F* >= 10

    breakdown = (gauge_primitives(a, gk, D, n_links, n_evo, "evolution")
                 + fermion_primitives(a, gk, n_links, V, n_evo, "evolution")
                 + gauge_primitives(a, gk, D, n_links, n_ramp, "ramp")
                 + fermion_primitives(a, gk, n_links, V, n_ramp, "ramp"))

    inter = {
        "stage": STAGE,
        "q_G_Z3": q_g, "q_G_Z3_compiled": g.link_qubits,                  # 2, 2
        "V_2028": V, "n_links_2028": n_links, "n_plaq_2028": n_plaq,      # 16, 32, 16
        "lq_gauge_2028": gauge, "lq_fermion_2028": fermion, "lq_anc_2028": anc,   # 64, 144, 22
        "lq_anc_hwp_2028": hl["rhop_full"]["hwp"]["ancilla"],             # 7 = 9 - w(9): one phasing group (E26)
        "lq_anc_rus_2028": RUS_ANCILLA,                                   # 1: repeat-until-success synthesis (E26)
        "lq_2028": lq,                                                    # 230 (R10, r25), printed exact
        "lq_anc_workspace_2028": anc,                                     # 22 = 2 x (7 + 4)
        "lq_2028_r22": lq_r22, "lq_anc_2028_r22": anc_r22,                # 216, 8: the register before R10
        "lq_headroom_2028": a.rfi_lq_2028.hi - lq,                        # 20
        "workspace_links_2028": n_conc,                                   # 2
        # conventions
        "eps_syn": float(a.eps_syn),                                      # 1e-2 (R-TOL)
        "eps_rot": eps,                                                   # 1.624e-3 = sqrt(1e-2 / 3792) (R-TOL)
        "eps_rot_round_e": float(a_round_e.eps_rot),                      # 1e-3, superseded
        "t_per_rotation_papers": gs["t_per_rotation_papers"],             # legacy record, 1.15 log2(1/eps)
        "t_per_rotation_report": t_per_rotation(eps),                     # every rotation, full fit (E20)
        "t_per_toffoli": gs["t_per_toffoli"],                             # 7
        # Z3 primitives (groups.py: dihedral basis, Z2 dropped; U_FFT round E), at the R-TOL eps
        "z3_t_U_inv": g.primitives["U_inv"].t(eps),                       # 0 (SWAP)
        "z3_t_U_mul": g.primitives["U_mul"].t(eps),                       # 28 (4 Toffoli)
        "z3_t_U_Tr": g.primitives["U_Tr"].t(eps),                         # 1 rotation
        "z3_t_U_phi": g.primitives["U_phi"].t(eps),                       # 1 rotation
        "z3_t_U_FFT": ft_new.t(eps),                                      # 4 T + 2 rotations (derived here)
        "z3_t_U_FFT_report": ft_new.t_report(eps),                        # full fit, record
        "z3_t_U_FFT_round_e": ft_new.t_papers(float(a_round_e.eps_rot)),  # 26.92 at 1e-3 (slope-only print, record)
        "z3_t_U_FFT_at_1e-4": ft_new.t_papers(eps_old),                   # 34.56 (slope-only print, record)
        "z3_fft_exact_t": ft_new.t_gates, "z3_fft_rotations": ft_new.rot, # 4, 2
        "z3_t_U_F_transpiler": ft_old.t(eps),                             # 14 rotations (arxiv_2408_00075), record
        "z3_t_U_F_transpiler_at_1e-4": ft_old.t_papers(eps_old),          # 213.9, the round-D headline (record)
        "z3_fft_reduction_at_1e-4": ft_old.t_papers(eps_old) / ft_new.t_papers(eps_old),   # 6.19 (record)
        "z3_fft_reduction": ft_old.t(eps) / ft_new.t(eps),
        "z3_fourier_rotations": ft_new.rot,                               # 2
        "z3_fourier_key": gs["fourier_key"],                              # 'U_FFT'
        "toffoli_per_link_step_2028": gs["toffoli_per_link"],             # 24
        "rotations_per_link_step_2028": gs["rotations_per_link"],         # 5.5
        # gauge terms per link, per step (ruling R2)
        "magnetic_t_per_link_2028": gs["magnetic"],                       # 177.9 -> '178'
        "electric_t_per_link_2028": gs["electric"],                       # 107.3 -> '107'
        "t_per_link_2028": gs["per_link"],                                # 285.2
        "electric_share_of_gauge_2028": gs["electric"] / gs["per_link"],  # 0.38
        "t_gauge_per_step_2028": t_gauge_step,                            # 9126.7 -> '9.1e3'
        # hopping and mass (r17: groups.hop_link_cost, Z3 = R-HOP structure at the full fit)
        "hop_k_2028": hl["k"],                                            # 9
        "hwp_rotations_k9": hl["rhop_full"]["hwp"]["n_rot"],              # 4
        "hwp_toffolis_k9": hl["rhop_full"]["hwp"]["n_toffoli"],           # 7
        "hwp_t_k9": hl["rhop_full"]["hwp"]["t"],                          # 127.1 -> '127 T'
        "hwp_ancilla_k9": hl["rhop_full"]["hwp"]["ancilla"],              # 7 (E26, r22; was 11)
        "hop_pauli_strings_2028": hl["rhop_full"]["n_pauli"],             # 8
        "hop_moves_per_link_2028": hl["rhop_full"]["n_moves"],            # 0 (diagonal link)
        "hop_toffoli_per_link_2028": hl["toffoli"],                       # 56 = 8 x 7
        "hop_rotations_per_link_2028": hl["n_rot"],                       # 32 = 8 x 4
        "t_hop_per_link_2028": hl["t"],                                   # 1016.6 -> '1017'
        "t_hop_per_link_2028_rhop_print": hl["rhop"]["t"],                # slope-only R-HOP at this eps, record
        "t_hop_per_link_2028_at_1e-4": hl_old["rhop"]["t"],               # 881.0, the round-D print
        "t_hop_per_link_2028_round_e": hop_link_rhop(a_round_e, gk)["t"], # 758.7, the round-E print
        "t_hop_per_link_2028_full_fit_at_1e-4": hl_old["t"],              # 8 x HWP(9) at 1e-4, full fit
        "t_hop_per_step_2028": t_hop_step,                                # 23456 -> '2.3e4'
        "t_mass_per_step_2028": t_mass_step,                              # 1466 -> '1.5e3'
        "t_per_step_2028": t_step,                                        # 32429 -> '3.24e4'
        "hopping_share_of_step_2028": t_hop_step / t_step,                # 0.72
        # synthesis (ruling R-TOL): N rotations per shot, N eps^2 = eps_syn
        "synth_rotations_per_step_2028": n_rot_step,                      # 1264
        "synth_rotations_per_shot_2028": n_rot_shot,                      # 3792 -> '3.8e3'
        "synthesis_error_2028": synth_err,                                # 1e-2 = eps_syn, by construction
        "eps_rot_for_budget_2028": eps_per_rotation(float(a.clean_shot_faults), n_rot_shot, a.synthesis_errors.value),   # 5.1e-3 (record: the whole 0.1)
        # steps (round-D ruling 1; round-E ruling E3)
        "n_trot_2028": n_evo,                                             # 2 (r18 ruling, overshoot accepted)
        "ramp_steps_2028": n_ramp,                                        # 1
        "n_steps_2028": n_steps,                                          # 3
        "t_evolution_2028": t_evo,                                        # 8.81e4
        "t_ramp_2028": t_ramp,                                            # 4.41e4
        "t_total_2028": t_total,                                          # 132,174.9 -> '1.3e5'
        "t_total_2028_over_cap": t_total / cap,                           # 1.322 -> '1.3x', an accepted overshoot
        "overshoot_accepted_2028": t_total > cap,                         # True (ruling 2026-10-01)
        "steps_in_cap_2028": steps_in_cap,                                # 2 < n_steps: the third step overshoots
        "steps_in_cap_exact_2028": cap / t_step,                          # 2.27
        "headroom_2028": cap - t_total,                                   # -32,175 T
        # the same circuit with fewer or more steps, each at its own R-TOL eps
        "t_total_2028_two_steps": t_total_two,                            # 87,266.3 (0.873x), the 1+1 alternative
        "t_total_2028_four_steps": t_total_four,                          # record
        "t_per_step_2028_three_steps": t_step,                            # 44,058.3 = the headline step
        "t_total_2028_three_steps": t_total,                              # 132,174.9 (1.32x) = the headline
        # the r17 print (1 ramp + 1 evolution, gauge tables slope-only), record
        "t_total_2028_r17": t_total_r17,                                  # 84,027.9 (0.84x)
        # the pre-r17 print (R-TOL, slope-only, 1 ramp + 2 evolution), record
        "t_per_step_2028_pre_r17": t_step_pre,                            # 32,429.5
        "t_total_2028_pre_r17": t_total_pre,                              # 97,288.5 (0.973x)
        # round E's fixed eps = 1e-3 (superseded by R-TOL), record
        "t_per_step_2028_round_e": t_step_round_e,                        # 33,446.26
        "t_total_2028_round_e": t_total_round_e,                          # 100,338.8, 1.003x
        # records: the step before round E and the E1 change alone
        "t_step_transpiler_fourier_2028": t_step_transpiler,              # transpiler U_F at the R-TOL eps
        "t_step_fft_at_1e-4_2028": t_step_old_eps_fft,                    # 3.83e4: U_FFT at eps = 1e-4
        "t_gauge_per_step_2028_round_d": t_gauge_step_round_d,            # 19801 (transpiler U_F, eps = 1e-4)
        "twosteps_t_total_2028": float(a.ramp_steps_2028 .value + a.twosteps_n_trot_2028.value)
                                 * (t_gauge_step_round_d + n_links * hl_old["rhop"]["t"]
                                    + V * mass_site_slope_only(a_old)["t"]),   # 99,509
        # the superseded round-D benchmark (1 ramp + 3 evolution), record only, at the round-E step
        "round_d_n_steps_2028": n_ramp + float(a.round_d_n_trot_2028),    # 4
        "round_d_t_total_2028": (n_ramp + float(a.round_d_n_trot_2028)) * t_step_round_e,   # 1.34e5 at the round-E step
        # physical window (the chapter's step, Delta t = t_had / N_Trot of the 2033 box)
        "dt_over_a": dt_over_a,                                           # 0.1
        "t_window_2028_over_a": n_evo * dt_over_a,                        # 0.2 -> 't = 0.2 a' (two evolution steps)
        # combined runs at several step sizes (ruling E4)
        "dt_multi_over_a": dts,                                           # (0.1, 0.2, 0.3)
        "t_reached_multi_over_a": t_reached,                              # (0.2, 0.4, 0.6)
        "t_sampled_multi_over_a": t_sampled,                              # (0.1, 0.2, 0.3, 0.4, 0.6)
        "t_overlap_multi_over_a": t_overlap,                              # (0.2,): 2 x 0.1 a and 1 x 0.2 a, the E4 check
        "t_max_multi_over_a": max(t_reached),                             # 0.6 = 3x the single-run window at 0.1 a
        # the algorithm-development gap (app10:94)
        "t_step_gap_target": gap,                                         # 1e3
        "gap_cut_factor_2028": t_step / gap,                              # 44.1 -> 'a 44x cut'
        # alternatives and checks, NOT the headline
        "transpiler_electric_t_per_link_2028": el_transpiler,             # 346.2
        "alt_abelian_plaquette_t": plaq_abelian,                          # 183 per plaquette, vs 351 by tab:primcost
        "tab_primcost_plaquette_t": gs["magnetic"] * n_links / n_plaq,    # 351
        "report_convention_t_gauge_per_step_2028": t_gauge_step_report,   # gauge step at the full fit, record
        "umul_style_hop_per_link_z3": umul_style_z3,                      # 357.7: 2.1x below the explicit 759
        "explicit_over_umul_style_z3": hl["t"] / umul_style_z3,           # 2.27
        # Ch. 9's current (R-HOP) hop rule at this chapter's k = 9: the same 8 x HWP(9) as hop_link
        "ch9_current_rule_hop_per_link_k9": c9_now["t"],                  # 880.99 = t_hop_per_link_2028_at_1e-4 (Ch. 9 at 1e-4)
        "ch9_current_rule_hop_strings": c9_now["n_pauli"],                # 8
        # Ch. 9's printed step AS IT NOW READS: read from ch09's Assumptions (Ch. 9 is re-priced under R-TOL separately)
        "ch9_t_per_step_printed": (float(c9a.t_per_step_2028_stated.lo),
                                   float(c9a.t_per_step_2028_stated.hi)),   # app06:108, whatever Ch. 9 now prints
        # Ch. 9's PRE-R-HOP counting rules on this lattice (2 strings per hop, 30 T per rotation): record only
        "ch9_rules_hop_group_k": c9["groups"][0][1], "ch9_rules_hop_groups": c9["groups"][0][2],      # 9, 64
        "ch9_rules_electric_group_k": c9["groups"][1][1], "ch9_rules_electric_groups": c9["groups"][1][2],   # 32, 3
        "ch9_rules_rotations_per_step": c9["n_rot"],                      # 672
        "ch9_rules_synth_per_step": c9["n_synth"],                        # 274
        "ch9_rules_toffoli_per_step": c9["n_toffoli"],                    # 541
        "ch9_rules_t_per_rot": c9["t_per_rot"],                           # 30
        "ch9_rules_t_hopping_per_step": c9["t_hopping"],                  # 10816 (undressed: 2 strings)
        "ch9_rules_t_electric_per_step": c9["t_electric"],                # 1191
        "ch9_rules_t_per_step": c9["t_step"],                             # 12007
        "ch9_rules_dressed_hop_per_step": c9_dressed_hop,                 # 43264: 8 strings at Ch. 9's 30 T
        "ch9_rules_ancilla_hop_group": c9["ancilla_hop_group"],           # 7 (E26, r22; was 11)
        "ch9_rules_ancilla_max_group": c9["max_ancilla"],                 # 31 (electric, k=32; was 37)
        # the retired working figure, for the record
        "retired_hop_2028": float(a.retired_hop_2028),                    # 5e3
        "hop_over_retired_2028": t_hop_step / float(a.retired_hop_2028),  # 4.86
        "legacy_t_per_step_2028": retired_step,                           # 2.48e4
        "legacy_t_total_2028": float(a.retired_n_trot_2028) * retired_step,   # 4.96e5 (the APPLY-stage print)
        "legacy_t_total_2028_derive": float(a.retired_n_trot_2028) * float(a.retired_hop_2028),   # 1e5 (the DERIVE-stage print)
        # shots and wall time
        "shots_2028": shots,                                              # 1e3
        "delta_rel_at_shots_2028": delta_at_shots,                        # 0.10, the printed '~10%/channel'
        "shots_for_delta_rel_2028": shots_for_delta,                      # (1000, 1000); at the old 20-30%: (111, 250)
        "wall_per_shot_2028_s": wall_shot,                                # 0.13227 -> '~0.13 s'
        "wall_2028_s": wall,                                              # 132.27, serial on one machine
        "wall_2028_min": wall / 60.0,                                     # 2.20 -> '~2 min'
        "shot_overhead_s": t0,                                            # 1e-4
        "eps_l_required_2028": eps_l_req,                                 # 1.028e-6 -> '1e-6'
        "faults_per_shot_2028_at_rfi_eps_l": t_total * float(a.rfi_eps_l),   # 9.7e-4
        # R9/R10 depth exports (common.depth_exports; factory.json Ch. 6 2028 structure)
        **{k: dx[k] for k in ("t_per_shot", "t_depth_per_shot", "f_star", "floor_wall_s", "factories_for_1yr")},
        "wall_serial_s": dx["wall_serial_s"],
        "baseline_ok": dx["baseline_ok"],
        "wall_first_result_s": None,                                      # one tier (R1 does not apply)
        "wall_campaign_s": (wall_c, wall_c),                              # 132.3 s
        "t_depth_per_shot_r22": d_shot_r22,                               # 8.0e4 at the 8-ancilla register
        "f_star_r22": t_total / d_shot_r22,                               # 1.65 (factory.json)
        "toffolis_per_shot_2028": n_tof_shot,                             # 8,016
        "exact_t_per_shot_2028": exact_t,                                 # 768
    }
    notes = (
        f"Register (eq:Nq_collider): {q_g} x {V} x {D} = {gauge} gauge + {n_c} x {n_stag} x {V} = {fermion} fermion + "
        f"{anc} ancilla ({n_conc} hop links at once, each one k={hl['k']} phasing group's {hwp9['ancilla']} + "
        f"{hwp9['n_rot']} repeat-until-success ancilla; R10) = {lq}; the box prints '230'. With the r22 register's "
        f"{anc_r22} ancilla ({lq_r22} LQ) every rotation runs in series, D_T = {d_shot_r22:.3g} and F* = "
        f"{t_total / d_shot_r22:.2f}; with {anc}, D_T = {d_shot[0]:.4g}-{d_shot[1]:.4g} and F* = {dx['f_star'][0]:.1f}-"
        f"{dx['f_star'][1]:.1f}.",
        f"Gauge terms (R2 tables, E20 full fit, eps = {eps:.4g}): {gs['magnetic']:.1f} "
        f"magnetic + {gs['electric']:.1f} electric = {gs['per_link']:.1f} T per link per step, {t_gauge_step:.5g} T per "
        "step on {n_links} links. Z3 primitives: dihedral basis of arxiv_2108_13305 with the Z2 dropped, circuits of "
        f"arxiv_2408_00075; the Fourier transform is U_FFT, 4 T + 2 rotations = {ft_new.t(eps):.2f} T (derived "
        f"here, round E; the transpiler U_F would be {ft_old.t(eps):.1f}).".replace("{n_links}", str(n_links)),
        f"Hopping (r17, groups.hop_link_cost; Z3 has no draft entry, R-HOP structure): {hl['toffoli']} Toffoli + "
        f"{hl['n_rot']} rotations = {hl['t']:.2f} T per link, {t_hop_step:.5g} per step; mass {V} x {ms['t']:.4f} = "
        f"{t_mass_step:.5g}.",
        f"Per shot: ({n_ramp:.0f} ramp + {n_evo:.0f} evolution) x {t_step:.6g} = {t_total:.6g} T, {t_total / cap:.4f}x "
        f"the 1e5 first-generation reference, an overshoot accepted by ruling (2026-10-01); the reference holds "
        f"{steps_in_cap} steps ({cap / t_step:.3f}). With 1 ramp + 1 evolution step, at its own eps: {t_total_two:.6g} T "
        f"({t_total_two / cap:.3f}x). r17 (1 + 1, gauge slope-only): {t_total_r17:.6g} T. Pre-r17 (slope-only, 3 steps): "
        f"{t_total_pre:.6g} T.",
        f"Synthesis (R-TOL, randomized): {n_rot_shot:.0f} rotations per shot, eps_rot = sqrt({float(a.eps_syn):g} / "
        f"{n_rot_shot:.0f}) = {eps:.4g}, {t_rot:.4f} T per gauge rotation (full fit, E20), "
        f"{gs['t_per_rotation_fermion']:.4f} per hop/mass rotation (full fit); N eps^2 = {synth_err:.2g}.",
        f"Window: {n_evo:.0f} evolution step(s) at Delta t = {dt_over_a:g} a reach t = {n_evo * dt_over_a:g} a; runs at "
        f"Delta t = {', '.join(f'{d:g}' for d in dts)} a reach {', '.join(f'{x:g}' for x in t_reached)} a (ruling E4)"
        + (f", overlapping at t = {', '.join(f'{x:g}' for x in t_overlap)} a." if t_overlap else "; no time is reached twice."),
        f"Shots: {shots:.0e} at n_s = {n_s} is a {100 * delta_at_shots:.0f}% binomial error, the printed '~10%/channel' "
        f"({shots_for_delta[0]:.0f} shots buy it; ruling R11 item 3 (a), was 20-30%, which needed 111-250).",
        f"Wall time {wall:.0f} s = {wall / 60:.1f} min on one machine at 1 us per T plus {t0:g} s per shot "
        "(box '~2 min', '~0.13 s/shot').",
    )
    return Result(era="2028", lq=(lq, lq), hard_ops=(t_total, t_total), breakdown=breakdown, intermediates=inter,
                  shots=(shots, shots), wall_time_s=(wall, wall), epsilon_l=(eps_l_req, eps_l_req), notes=notes)


# --------------------------------------------------------------------------- #
# 2033: one circuit, two readouts. 2+1D Sigma(72x3), 4x8 cut from 4x16, 3 staggered fields
# --------------------------------------------------------------------------- #

def pair_phase(w: float, dt: float, n: int) -> float:
    """Second-order dispersion phase of one free gauge quantum of frequency w after n steps of dt (lattice units)."""
    return n * (2 * math.asin(w * dt / 2) - w * dt)


def rule_steps_top_mode(a: Assumptions) -> dict:
    """Step rule at the top box mode (n = 3, 4.65 GeV): the fewest steps over the window's upper end (50 a) at which
    the injected gluon's dispersion phase is <= eps; never fewer than the chosen 0.1 a grid gives."""
    L_par = int(a.L_par_priced.value)
    w = 2 * math.sin(math.pi * 3 / L_par)
    t = float(a.window_over_a.hi)
    n0 = int(round(t / 0.1))
    n = n0
    if a.trotter_rule.value == "state":
        while abs(pair_phase(w, t / n, n)) > float(a.trotter_eps.value):
            n += 1
    return {"w": w, "t_over_a": t, "n_grid": n0, "n": n, "phase": abs(pair_phase(w, t / n, n)),
            "phase_grid": abs(pair_phase(w, 0.1, n0))}


def free_gauge_frequencies_h(dims, ham: str = "KS") -> list:
    """Free transverse lattice-gluon frequencies (lattice units, k != 0). KS: w^2 = sum_i khat_i^2, khat_i =
    2 sin(pi n_i / L_i) (= ch05.free_gauge_frequencies). I (s648, DERIVED HERE, an assumption about H_I's free limit):
    the tree-level Symanzik dispersion w^2 = sum_i (khat_i^2 + khat_i^4 / 12), O(a^2)-improved. On 3^3 every nonzero
    component has khat^2 = 3, so w^2 rises by 1.25 and Lambda_sd by about 1.4."""
    import itertools
    ws = []
    for n in itertools.product(*[range(L) for L in dims]):
        kh = [4.0 * math.sin(math.pi * ni / L) ** 2 for ni, L in zip(n, dims)]
        w2 = sum(kh) + (sum(x * x for x in kh) / 12.0 if ham == "I" else 0.0)
        if w2 > 1e-12:
            ws.append(math.sqrt(w2))
    return ws


def sd_lambda_h(dims, n_adj: int, temp: float, ham: str = "KS") -> float:
    """Ch. 5's state-dependent Lambda_sd (ch05.trotter_state_dependent) on the free frequencies of `ham`."""
    mult = n_adj * (len(dims) - 1)
    var = mult * sum(2.0 * (0.5 * w ** 3 / math.tanh(w / (2.0 * temp))) ** 2 for w in free_gauge_frequencies_h(dims, ham))
    return math.sqrt(var) / 8.0


def free_err_h(dims, n_adj: int, dt: float, n: int, temp: float, ham: str = "KS") -> float:
    """Ch. 5's exact free-field error (ch05.free_gauge_trotter_error) on the free frequencies of `ham`."""
    from estimates import ch05_qgp_transport as ch5
    mult = n_adj * (len(dims) - 1)
    logf = sum(mult * math.log(ch5.free_gauge_mode_overlap(w, dt, n, temp)) for w in free_gauge_frequencies_h(dims, ham))
    return math.sqrt(2.0 * (1.0 - math.exp(logf)))


def dipole_rule_steps(a: Assumptions) -> dict:
    """Step rule for the thermal dipole rows (literal state-dependent estimate, as Ch. 5): each shot max(t / Delta t,
    N_sd) steps, N_sd = ceil(sqrt(Lambda_sd t^3 / eps)), t in lattice units, a T = 0.446. s648: the static row is H_I,
    evaluated on the improved free dispersion (free_gauge_frequencies_h); the light-like 2O row stays H_KS."""
    from estimates import ch05_qgp_transport as ch5
    da, ddt = float(a.dipole_a_fm), float(a.dipole_dt_fm)
    temp = float(a.dipole_T_GeV) * da / 0.1973269804
    eps = float(a.trotter_eps.value)
    nadj8 = int(ch5.Assumptions().trotter_n_adj.value)
    out = {}
    for key, dims, nadj, ham in (("static", tuple(a.static_dims.value), nadj8, a.hamiltonian_dipole_static.value),
                                 ("light", tuple(a.lightlike_dims.value), 3, a.hamiltonian.value)):
        lam = sd_lambda_h(dims, nadj, temp, ham)
        rows = []
        for t in a.dipole_times_fm.value:
            n0 = int(round(t / ddt))
            tl = t / da
            n = max(n0, math.ceil(round(math.sqrt(lam * tl ** 3 / eps), 9))) if a.trotter_rule.value == "state" else n0
            rows.append((n0, n, lam * tl * (tl / n) ** 2))
        out[key] = {"lambda": lam, "rows": tuple(rows), "hamiltonian": ham}
    return out


def trotter_check_2033(a: Assumptions) -> dict:
    """Referee G7 (2026-10-04, DERIVED HERE): the three Trotter-error estimates of Ch. 5's trotter_check at the chosen
    Delta t = 0.1 a over the 30-50 a window (300-500 steps) on the priced 4 x 8 lattice (64 links), and for the dipole
    rows (Delta t = 0.05 a, a T = 0.45, 50 and 100 steps; 3^3 Sigma(72x3) with 8 color copies, 3x3x5 2O with 3).
    (i) worst case (same-link nested commutators, hop n_stag N_c / 2 per link); (ii) state-dependent connected variance
    in the free-field vacuum; (iii) exact free-field error in the vacuum. The injected pair: the second-order dispersion
    phase N (2 arcsin(w dt / 2) - w dt) of one gauge quantum at the box modes k = 2 pi n / L_par, n = 1, 2, 3
    (1.55, 3.1, 4.65 GeV). Fermions: exact free staggered numerics in the tests (under 10% of the vacuum infidelity)."""
    from estimates import ch05_qgp_transport as ch5
    D = int(a.dim_2033.value)
    dims = (int(a.L_perp.value), int(a.L_par_priced.value))
    links = dims[0] * dims[1] * D
    dt = 0.1                                                              # Delta t / a, the chapter's step
    ns = tuple(int(round(x / dt)) for x in (a.window_over_a.lo, a.window_over_a.hi))   # (300, 500)
    eps, g2, temp = float(a.trotter_eps.value), float(a.trotter_g2.value), float(a.trotter_temp_lat.value)
    a5 = ch5.Assumptions()
    cmax = float(a5.sigma72_emax_over_fund.value) * 4.0 / 3.0
    nadj = int(a5.trotter_n_adj.value)
    hop = float(a.n_stag.value) * float(a.n_c.value) / 2.0
    wc = [ch5.trotter_worst_case(links, D, g2, cmax, hop, n * dt, dt, eps) for n in ns]
    wmin = [min(ch5.trotter_worst_case(links, D, x / 100, cmax, hop, n * dt, dt, eps)["err"] for x in range(20, 501))
            for n in ns]
    sd = [ch5.trotter_state_dependent(dims, nadj, temp, n * dt, dt, eps) for n in ns]
    disp = {}
    for nk in (1, 2, 3):
        w = 2 * math.sin(math.pi * nk / dims[1])
        disp[nk] = tuple(n * (2 * math.asin(w * dt / 2) - w * dt) for n in ns)
    da = float(a.dipole_a_fm)
    ddt = float(a.dipole_dt_fm) / da                                      # 0.05 a
    dtemp = float(a.dipole_T_GeV) * da / 0.1973269804                     # 0.446
    steps_d = tuple(int(round(t / float(a.dipole_dt_fm))) for t in a.dipole_times_fm.value)   # (50, 100)
    dip_cases = ((tuple(a.static_dims.value), nadj, a.hamiltonian_dipole_static.value),
                 (tuple(a.lightlike_dims.value), 3, a.hamiltonian.value))   # s648: static on H_I
    out = {
        "trotter_steps_window": ns,
        "trotter_worst_lambda_2033": wc[0]["lambda"],                     # 2092 a^-3 (Wilson-matched plaquette, 6/g^2)
        "trotter_worst_err_2033": tuple(x["err"] for x in wc),            # (628, 1046)
        "trotter_worst_err_min_g2_2033": tuple(wmin),
        "trotter_worst_steps_2033": tuple(x["n_steps_at_eps"] for x in wc),   # (23769, 51143)
        "trotter_sd_lambda_2033": sd[0]["lambda"],                        # 15.0 a^-3
        "trotter_sd_err_2033": tuple(x["err"] for x in sd),               # (4.5, 7.5)
        "trotter_sd_steps_2033": tuple(x["n_steps_at_eps"] for x in sd),  # (2011, 4326)
        "trotter_fast_phase_2033": tuple(x["fast_phase"] for x in sd),    # (0.28, 0.47): w_max = 2 sqrt 2
        "trotter_free_err_2033": tuple(ch5.free_gauge_trotter_error(dims, nadj, dt, n, temp) for n in ns),   # (0.046, 0.038)
        "trotter_pair_phase_by_box_mode": disp,                           # n=3: (0.079, 0.132)
        "trotter_dipole_free_err": tuple(free_err_h(dd, na, ddt, n, dtemp, hm)
                                         for dd, na, hm in dip_cases for n in steps_d),
        "trotter_dipole_sd_err": tuple(sd_lambda_h(dd, na, dtemp, hm) * (n * ddt) * ddt * ddt
                                       for dd, na, hm in dip_cases for n in steps_d),   # the 0.01 fm/c grid fails
    }
    # step rule: the dipole rows at their rule steps (free-field check at the same steps)
    dr = dipole_rule_steps(a)
    out["trotter_dipole_rule_steps"] = tuple(x[1] for k in ("static", "light") for x in dr[k]["rows"])   # (72, 203, 64, 181)
    out["trotter_dipole_rule_sd_err"] = tuple(x[2] for k in ("static", "light") for x in dr[k]["rows"])
    out["trotter_dipole_rule_free_err"] = tuple(
        free_err_h(dd, na, (t / da) / n, n, dtemp, hm)
        for (dd, na, hm), k in zip(dip_cases, ("static", "light"))
        for t, (_, n, _) in zip(a.dipole_times_fm.value, dr[k]["rows"]))
    top = rule_steps_top_mode(a)
    out["trotter_rule"] = a.trotter_rule.value
    out["trotter_top_mode_steps_50a"] = top["n"]                          # 575
    out["trotter_top_mode_phase_50a"] = top["phase"]                      # 0.0997 (0.132 at 500 steps)
    out["trotter_top_mode_dt_over_a"] = top["t_over_a"] / top["n"]        # 0.0870
    out["trotter_free_err_top_50a"] = ch5.free_gauge_trotter_error(dims, nadj, top["t_over_a"] / top["n"], top["n"], temp)
    return out


def _model_2033(a: Assumptions) -> Result:
    gk = a.group_2033.value
    g = GROUPS[gk]
    q_g = link_width(a, gk)                                               # 9 (R1)
    q_g_printed_before = int(g.link_qubits_chapter.lo)                    # 8, the pre-ruling print
    D = int(a.dim_2033.value)                                             # 2 (guarded)
    n_c, n_stag, n_anc = int(a.n_c.value), int(a.n_stag.value), int(a.n_anc_2033.value)
    L_perp = int(a.L_perp.value)

    V_uncut = L_perp * int(a.L_par_uncut.value)                           # 64
    gauge_u, fermion_u, _, lq_uncut = register(q_g, V_uncut, D, n_c, n_stag, n_anc)     # 1152, 576, 1828
    V_cut = (L_perp * int(a.L_par_cut.lo), L_perp * int(a.L_par_cut.hi))  # (32, 40)
    reg_cut = tuple(register(q_g, v, D, n_c, n_stag, n_anc) for v in V_cut)
    lq_cut = tuple(r[3] for r in reg_cut)                                 # (964, 1180): sensitivity; the box quotes lq below
    fermion_cut = tuple(r[1] for r in reg_cut)                            # (288, 360)
    gauge_cut = tuple(r[0] for r in reg_cut)                              # (576, 720)
    n_anc_samp = int(a.n_anc_2028.value)                                  # 'the phasing workspace plus quench workspace': the 2028 value, 8
    lq_cut_samp = tuple(register(q_g, v, D, n_c, n_stag, n_anc_samp)[3] for v in V_cut)   # (872, 1088)
    lq_cut_legacy = tuple(register(q_g_printed_before, v, D, n_c, n_stag, n_anc)[3] for v in V_cut)   # (900, 1100)

    # T-count, priced at V = 32 (app10:95); the box quotes the register at this volume only (ruling R11 option (A))
    V = L_perp * int(a.L_par_priced.value)
    lq = register(q_g, V, D, n_c, n_stag, n_anc)[3]                       # 964
    n_links = V * D                                                       # 64
    n_plaq = n_plaquettes(V, D)                                           # 32
    n_evo_full = float(a.n_trot_2033)                                     # 1e3: the t = 100 a window (record since R4)
    ramp = (float(a.ramp_steps_2033.lo), float(a.ramp_steps_2033.hi))
    dt_over_a = float(a.t_had_over_a) / n_evo_full                        # 0.1, the chapter's step
    # R4 (r25): the window is 30-50 a, so 300-500 evolution steps at Delta t = 0.1 a; the ramp range is unchanged
    n_evo_end = (float(round(a.window_over_a.lo / dt_over_a)), float(round(a.window_over_a.hi / dt_over_a)))   # (300, 500)
    n_evo = n_evo_end[0]
    # r26 (referee J1): prep ramp + the same ramp backwards for the species readout map
    n_ramps = 1 + int(a.readout_ramps_2033.value)                         # 2
    n_steps = (n_evo_end[0] + n_ramps * ramp[0], n_evo_end[1] + n_ramps * ramp[1])   # (500, 2500); r25 (400, 1500)
    n_steps_r25 = (n_evo_end[0] + ramp[0], n_evo_end[1] + ramp[1])        # (400, 1500): the r25 headline circuit (record)
    n_steps_full = (n_evo_full + ramp[0], n_evo_full + ramp[1])           # (1100, 2000): the r23 circuit (record)

    # ruling R-TOL: each end of the ramp range is its own circuit, with its own rotation count and tolerance
    n_rot_step = synth_rotations_per_step(a, gk, n_links, V)              # 261,792 (independent of eps)
    n_rot_shot = (n_steps[0] * n_rot_step, n_steps[1] * n_rot_step)       # (2.88e8, 5.24e8)
    a_round_e = a                                                         # eps = 1e-3, round E: record only
    eps_end = (rtol_eps(a, n_rot_shot[0]), rtol_eps(a, n_rot_shot[1]))    # (5.893e-6, 4.370e-6)
    a_lo, a_hi = _at_eps(a, eps_end[0], n_rot_shot[0]), _at_eps(a, eps_end[1], n_rot_shot[1])
    a = a_lo                                  # the single-valued per-link intermediates are at the low end's eps
    synth_err = (synthesis_error(a_lo, n_rot_shot[0]), synthesis_error(a_hi, n_rot_shot[1]))   # (1e-2, 1e-2)
    eps_budget = (eps_per_rotation(float(a.clean_shot_faults), n_rot_shot[1], a.synthesis_errors.value),
                  eps_per_rotation(float(a.clean_shot_faults), n_rot_shot[0], a.synthesis_errors.value))   # (2.6e-5, 3.5e-5)

    def _step(ax, v):
        g_, h_, m_ = gauge_step(ax, gk, D), hop_link(ax, gk), mass_site(ax)
        return v * D * (g_["per_link"] + h_["t"]) + v * m_["t"], g_, h_, m_

    gs = gauge_step(a, gk, D)
    hl, ms = hop_link(a, gk), mass_site(a)
    t_gauge_step = n_links * gs["per_link"]
    t_hop_step = n_links * hl["t"]
    t_mass_step = V * ms["t"]
    t_rest_step = t_hop_step + t_mass_step
    t_step = _step(a, V)[0]                                               # 9.754e6 (low end's eps); same sum as the hi end
    assert math.isclose(t_step, t_gauge_step + t_rest_step, rel_tol=1e-12)
    t_step_hi_end, gs_hi, hl_hi, ms_hi = _step(a_hi, V)                   # 1.0104e7 (high end's eps)
    t_step_end = (t_step, t_step_hi_end)
    t_evo_end = (n_evo_end[0] * t_step, n_evo_end[1] * t_step_hi_end)
    t_evo = t_evo_end[0]
    # the r23 circuit (t = 100 a, 1100-2000 steps), each end at its own R-TOL eps: the record of the old headline
    t_total_full = tuple(n * _step(_at_eps(a, rtol_eps(a, n * n_rot_step), n * n_rot_step), V)[0] for n in n_steps_full)
    t_ramp = (ramp[0] * t_step, ramp[1] * t_step_hi_end)
    t_total = (n_steps[0] * t_step, n_steps[1] * t_step_hi_end)           # r26: (4.9e9, 2.5e10); r25 (3.902e9, 1.506e10)
    t_ramps = (n_ramps * t_ramp[0], n_ramps * t_ramp[1])                  # prep + backward readout ramp
    # the r25 headline (no readout map), each end at its own R-TOL eps: record
    t_total_r25 = tuple(n * _step(_at_eps(a, rtol_eps(a, n * n_rot_step), n * n_rot_step), V)[0] for n in n_steps_r25)
    # J1: the electric-basis measurement after the backward ramp, one Fourier transform per link (not added; the
    # Sigma(72x3) transform is the paper's lower bound)
    t_readout_fourier = n_links * GROUPS[gk].primitives[fourier_key(gk)].t(float(a.eps_rot))
    # J2: the gauge-invariant source and the box geometry. r26 v1 priced a single nearest-neighbor bilinear (one
    # link's hop, record); v2 is the wave packet sum_d f_p(d) psibar_x W_{x,x+d} psi_{x+d}, d = 1..L_par, priced
    # at L_par bilocal insertions of the full L_par-link Wilson line (an upper bound per separation)
    t_source_per_step = hl["t"]
    t_source_wavepacket = int(a.L_par_priced.value) * insertion_t(a, gk)  # 8 x 1.817e4 = 1.45e5
    a_fm_ = float(a.a_fm)
    recollision_fm = int(a.L_par_priced.value) * a_fm_ / 2                # 0.4 fm/c: back-to-back pair meets again
    inv_lam_hi = float(a.hbarc_GeV_fm) / float(a.lambda_qcd_GeV.lo)       # 0.99 fm/c
    L_par_sep = math.ceil(2 * inv_lam_hi / a_fm_ - 1e-9)                  # 20 sites: separation 2 c t through 1/Lambda
    lq_sep = register(q_g, L_perp * L_par_sep, D, n_c, n_stag, n_anc)[3]  # 2260
    cap = float(a.rfi_t_2033)
    steps_in_cap = cap / t_step                                           # 102.5 at the low end's eps
    # round E's fixed eps = 1e-3, for the record
    t_step_round_e = _step_slope_only(a_round_e, gk, D, n_links, V)       # 1.95852e6 (slope-only, record)
    t_total_round_e = (n_steps_full[0] * t_step_round_e, n_steps_full[1] * t_step_round_e)   # (2.154e9, 3.917e9), r23 window
    # the pre-r17 print (R-TOL at the R-HOP rotation count, slope-only), record
    n_rot_step_pre = _n_rot_step_pre_r17(a_round_e, gk, n_links, V)       # 74,592
    t_total_pre = tuple(n * _step_slope_only(_at_eps(a_round_e, rtol_eps(a_round_e, n * n_rot_step_pre),
                                                     n * n_rot_step_pre), gk, D, n_links, V)
                        for n in n_steps_full)                            # (2.768e9, 5.106e9), r23 window
    # the hop variants (groups.fermion_hop_counts readings), at the low end's eps: sensitivities, not the headline.
    # r19 (E21): the headline is share='link' + undo; 'share_draft' is the r17-r18 per-application recompute
    hop_var = {name: hop_link_cost(gk, n_stag, eps=float(a.eps_rot), legacy=False, **kw)
               for name, kw in (("share_draft", {"share": "draft"}), ("no_undo", {"share": "link", "undo": False}),
                                ("draft_mcx", {"share": "link", "mcx": "draft-color"}),
                                ("share_draft_no_undo", {"share": "draft", "undo": False}))}

    # the same pricing at the upper end of the register range (V = 40): not in the chapter. A different circuit,
    # so its own rotation count and R-TOL tolerance at each end of the ramp range
    V_hi = V_cut[1]
    n_rot_step_vhi = synth_rotations_per_step(a, gk, V_hi * D, V_hi)
    t_total_hi = tuple(n * _step(_at_eps(a, rtol_eps(a, n * n_rot_step_vhi), n * n_rot_step_vhi), V_hi)[0]
                       for n in n_steps)                                  # (1.378e10, 2.538e10)
    t_step_hi = t_total_hi[0] / n_steps[0]

    # readout (b): one or two controlled insertions of the bilocal current per shot (derived here), not per step
    t_ins = insertion_t(a, gk)                                            # 1.817e4
    n_ins = (float(a.insertions_per_shot_b.lo), float(a.insertions_per_shot_b.hi))

    # the pre-ruling chain, for the record: 32 plaquettes x 5e3 + 3.4e5 = 5e5 per step
    t_plaq_legacy = float(g.plaquette_t)                                  # 5e3 (groups.py, the chapters' retired figure)
    retired_rest = float(a.retired_hop_insertion_2033)                    # 3.4e5 (RETIRED)
    legacy_step = n_plaq * t_plaq_legacy + retired_rest                   # 5e5
    applied_step = t_gauge_step + retired_rest                            # 2.136e6, the APPLY-stage step
    legacy_total = (n_steps_full[0] * legacy_step, n_steps_full[1] * legacy_step)   # (5.5e8, 1e9)

    # the unpublished draft's staggered hopping, NOT used in the headline
    l2 = math.log2(1.0 / float(a.eps_rot))
    draft_hop_link = float(a.draft_hop_S72x3_t_const) + float(a.draft_hop_S72x3_t_log) * l2   # printed C^S, 5.76e4
    draft_hop_step = n_links * draft_hop_link                             # 3.7e6 (one field)
    draft_hop_step_3fields = n_stag * draft_hop_step                      # 1.1e7: the printed C^S is per staggered field
    draft_step = t_gauge_step + draft_hop_step                            # one field, as the APPLY stage used it
    draft_total = (n_steps_full[0] * draft_step, n_steps_full[1] * draft_step)

    # kinematics of the priced window
    a_fm = float(a.a_fm)
    t_had_fm = float(a.t_had_over_a) * a_fm                               # 10 fm/c (the r23 window)
    window_fm = (a.window_over_a.lo * a_fm, a.window_over_a.hi * a_fm)    # (3, 5) fm/c (R4)
    hbarc = float(a.hbarc_GeV_fm)
    p_unit = tuple(2 * math.pi * hbarc / (L * a_fm) for L in
                   (int(a.L_par_cut.lo), int(a.L_par_cut.hi), int(a.L_par_uncut.value)))   # 1.55, 1.24, 0.77 GeV
    crossings = (float(a.t_had_over_a) / a.L_par_cut.hi, float(a.t_had_over_a) / a.L_par_cut.lo)   # (10, 12.5), r23 window
    crossings_window = (a.window_over_a.lo / a.L_par_priced.value, a.window_over_a.hi / a.L_par_priced.value)   # (3.75, 6.25)
    inv_lambda_fm = (hbarc / a.lambda_qcd_GeV.hi, hbarc / a.lambda_qcd_GeV.lo)                     # (0.66, 0.99) fm/c

    # step rule (2026-10-04): the 4.65 GeV point at 50 a runs n_top evolution steps (575, was 500); its own circuit,
    # rotation count and R-TOL tolerance. The 30 a end meets the rule at 0.1 a (phase 0.079), so only this shot moves
    top = rule_steps_top_mode(a)
    n_steps_top = top["n"] + n_ramps * ramp[1]                            # 2575
    n_rot_top = n_steps_top * n_rot_step
    a_top = _at_eps(a, rtol_eps(a, n_rot_top), n_rot_top)
    t_step_top, gs_top, _, _ = _step(a_top, V)
    t_top = n_steps_top * t_step_top                                      # 2.61e10 -> '2.6e10'
    hard_ops_2033 = (t_total[0], max(t_total[1], t_top))

    # logical faults (app10:90)
    eps_l = float(a.rfi_eps_l)
    faults_rfi = (t_total[0] * eps_l, hard_ops_2033[1] * eps_l)           # (49, 261)
    eps_l_req = (float(a.clean_shot_faults) / hard_ops_2033[1], float(a.clean_shot_faults) / t_total[0])   # (3.8e-12, 2.0e-11)

    # shots: the r23 accounting (1/(n_s delta^2), 5-10%), kept as the record
    n_s = float(a.n_s)
    shots_a = (shots_binomial(n_s, a.delta_rel_2033.hi), shots_binomial(n_s, a.delta_rel_2033.lo))   # (1e3, 4e3), r23
    shots_a_printed = (float(a.shots_samp_printed.lo), float(a.shots_samp_printed.hi))                 # (1e3, 4e3), r23
    a_inv2 = (a.hadamard_amplitude.hi ** -2, a.hadamard_amplitude.lo ** -2)                            # (25, 400)
    shots_b_base = shots_hadamard(1.0, float(a.eps_stat), float(a.n_z), float(a.n_tens), float(a.n_boost))   # 3e3
    shots_b = (shots_b_base * a_inv2[0], shots_b_base * a_inv2[1])                                     # (7.5e4, 1.2e6)
    shots_b_printed = (float(a.shots_ff_printed.lo), float(a.shots_ff_printed.hi))                     # (7.5e4, 1.2e6)

    gate = float(a.t_gate_s)
    t0 = float(a.shot_overhead_s)                                         # ~0.1 ms per shot (Ch. 1 convention)
    yr = float(a.seconds_per_year)
    day = 86400.0
    horizon = float(a.campaign_horizon_yr)
    horizon_s = horizon * yr                                              # 1.578e8 s
    wall_shot = (t_total[0] * gate + t0, t_total[1] * gate + t0)          # serial: (3.99e3, 1.52e4) s
    wall_shot_full = (t_total_full[0] * gate + t0, t_total_full[1] * gate + t0)   # (1.10e4, 2.02e4) s, r23 window

    # R9/R10: T-depth, schedule (P) of factory.json; each end is its own circuit at its own estimate
    t_rot_end = (gs["t_per_rotation_fermion"], gs_hi["t_per_rotation_fermion"])
    d_step = {(e, k): depth_2033_step(a, t_rot_end[e], n_links, k) for e in (0, 1) for k in (0, 1)}
    d_shot = (n_steps[0] * d_step[(0, 0)], n_steps[1] * d_step[(1, 1)])  # low circuit low estimate, high high
    f_star_circuit = ((t_total[0] / (n_steps[0] * d_step[(0, 1)]), t_total[0] / (n_steps[0] * d_step[(0, 0)])),
                      (t_total[1] / (n_steps[1] * d_step[(1, 1)]), t_total[1] / (n_steps[1] * d_step[(1, 0)])))   # 7.0-18.5
    # per-shot time with ~10 factories on one machine: max(N_T t_gate, D_T t_r) + t0 (the corrected wall)
    shot_s = (max(t_total[0] * gate, d_shot[0] * REACTION_TIME_S) + t0,
              max(t_total[1] * gate, d_shot[1] * REACTION_TIME_S) + t0)   # (3.99e3, 2.17e4) s
    depth_penalty = (shot_s[0] / wall_shot[0], shot_s[1] / wall_shot[1])  # (1.0, 1.43)
    # the top point's shot time (step rule), and the campaign's mean shot at the high end (one of three momenta)
    d_top = n_steps_top * depth_2033_step(a_top, gs_top["t_per_rotation_fermion"], n_links, 1)
    shot_s_top = max(t_top * gate, d_top * REACTION_TIME_S) + t0
    n_mom_ = int(a.n_source_momenta.value)
    shot_s_campaign_hi = (n_mom_ - 1) * shot_s[1] / n_mom_ + shot_s_top / n_mom_

    # R1 first result: one configuration (one spacing, one source momentum), B/M, V/PS, kaons at 30%
    delta_first = float(a.delta_first_2033)
    shots_first = ratio_shots(a, delta_first)                             # (122, 244)
    wall_first = (shots_first[0] * shot_s[0], shots_first[1] * shot_s[1])
    # R2 campaign: the 3-momentum trend at trend_sigma, one spacing (whether box effects cancel in the trend is open, J2)
    n_mom = int(a.n_source_momenta.value)                                 # 3
    delta_tr = trend_delta(a)                                             # 0.0589 per point
    shots_pt = ratio_shots(a, delta_tr)                                   # (3.17e3, 6.34e3) per configuration
    shots_trend = (n_mom * shots_pt[0], n_mom * shots_pt[1])              # (9.5e3, 1.9e4)
    wall_trend = (shots_trend[0] * shot_s[0], shots_trend[1] * shot_s_campaign_hi)   # step rule: the top point's shot
    trend_over_horizon = (wall_trend[0] / horizon_s, wall_trend[1] / horizon_s)   # (0.30, 4.4)
    configs_in_horizon = (horizon_s / (shots_pt[1] * shot_s_campaign_hi), horizon_s / (shots_pt[0] * shot_s[0]))
    sigma_at = {d: float(a.trend_size) / (math.sqrt(2) * d) for d in (0.05, 0.10, delta_first)}   # 3.5, 1.8, 0.59 sigma
    ratio_10 = ratio_shots(a, 0.10)                                       # (1.1e3, 2.2e3): the shot audit's range
    # a second spacing (requirements row): on the same 4x8 lattice it doubles the steps (box 0.2 x 0.4 fm) and adds
    # 2 x the one-spacing campaign; at fixed physical volume it needs 4x the sites
    wall_camp_2sp = (3 * wall_trend[0], 3 * wall_trend[1])               # 3x the trend (the L_par = 4 check aside)
    lq_fixed_volume_2sp = register(q_g, 4 * V, D, n_c, n_stag, n_anc)[3]  # 3556

    # readout (b), FCC clock: same circuit, corrected per-shot time
    wall_b = (shots_b[0] * shot_s[0], shots_b[1] * shot_s[1])
    ff_reduction = (wall_b[0] / horizon_s, wall_b[1] / horizon_s)         # serial years / 5
    ff_shots_in_horizon = (horizon_s / shot_s[1], horizon_s / shot_s[0])
    shots_per_ff_call = (float(a.eps_stat) ** -2 * a_inv2[0], float(a.eps_stat) ** -2 * a_inv2[1])   # (625, 1e4)
    ff_calls_in_horizon = (ff_shots_in_horizon[0] / shots_per_ff_call[1],
                           ff_shots_in_horizon[1] / shots_per_ff_call[0])

    # the r23 sampling record: 1e3-4e3 shots per configuration at the 100 a window, serial
    wall_a = (shots_a[0] * wall_shot_full[0], shots_a[1] * wall_shot_full[1])   # (1.1e7, 8.1e7), r23
    n_config = int(a.n_spacings.value) * int(a.n_source_momenta.value)    # 6 (requirements row)
    camp_a_10 = (n_config * shots_a[0] * wall_shot_full[0] / yr, n_config * shots_a[0] * wall_shot_full[1] / yr)   # r23
    camp_a_5 = (n_config * shots_a[1] * wall_shot_full[0] / yr, n_config * shots_a[1] * wall_shot_full[1] / yr)    # r23

    # R7 dipole rows (simplobs.json C' and E): pure gauge, each shot at its own time (R4)
    da, ddt = float(a.dipole_a_fm), float(a.dipole_dt_fm)
    times = (float(a.dipole_times_fm.lo), float(a.dipole_times_fm.hi))
    steps_d = tuple(int(round(t / ddt)) for t in times)                   # (50, 100): the 0.01 fm/c grid
    dr = dipole_rule_steps(a)                                             # step rule (2026-10-04): own steps per shot
    steps_static = tuple(x[1] for x in dr["static"]["rows"])              # (72, 203)
    steps_light = tuple(x[1] for x in dr["light"]["rows"])                # (64, 181)
    gk_dip, ham_dip = a.group_dipole_static.value, a.hamiltonian_dipole_static.value   # s648: S216x3, H_I
    # q1 (ruling B, 2026-10-05): the medium is the electric vacuum, coupling-ramped over one lattice unit and evolved
    # for c / T before the sources enter, both at the 0.01 fm/c grid (the thermalization leg needs no observable-level
    # Trotter accuracy; the dipole steps follow the rule); the circuit is priced whole (R-TOL eps on its length)
    invT_fm = hbarc / float(a.dipole_T_GeV)                               # 0.4485 fm/c
    n_ramp_q = math.ceil(float(a.quench_ramp_a) * da / ddt - 1e-9)        # 20
    c_range = (float(a.quench_c_range.lo), float(a.quench_c_range.hi))    # (2, 2 pi)
    n_th_q = tuple(math.ceil(c * invT_fm / ddt) for c in c_range)          # (90, 282)
    n_prep_q = tuple(n_ramp_q + n for n in n_th_q)                        # (110, 302)

    def _rows(n_prep):
        st = tuple(dipole_price(a, gk_dip, tuple(a.static_dims.value), n + n_prep, 2, int(a.static_src_qubits.value),
                                ham_dip) for n in steps_static)
        li = tuple(dipole_price(a, "2O", tuple(a.lightlike_dims.value), n + n_prep, 2 * (int(t / da) + 1),
                                int(a.lightlike_src_qubits.value)) for n, t in zip(steps_light, times))
        return st, li
    static_q, light_q = zip(*(_rows(n) for n in n_prep_q))               # rows at c = 2 and c = 2 pi
    static_evo, light_evo = _rows(0)                                      # record: the dipole evolution alone (pre-q1 print)

    def _rows_rule_step():
        # record only (verifier q1-1, unpriced): the prep leg at each row's own rule step instead of the 0.01 grid
        prep_fm = tuple(float(a.quench_ramp_a) * da + c * invT_fm for c in c_range)   # (1.10, 3.02) fm/c
        st = tuple(dipole_price(a, gk_dip, tuple(a.static_dims.value), n + math.ceil(pf / (t / n) - 1e-9), 2,
                                int(a.static_src_qubits.value), ham_dip)
                   for (n, t), pf in zip(zip(steps_static, times), prep_fm))
        li = tuple(dipole_price(a, "2O", tuple(a.lightlike_dims.value), n + math.ceil(pf / (t / n) - 1e-9),
                                2 * (int(t / da) + 1), int(a.lightlike_src_qubits.value))
                   for (n, t), pf in zip(zip(steps_light, times), prep_fm))
        n_prep = tuple(tuple(math.ceil(pf / (t / n) - 1e-9) for (n, t), pf in zip(zip(steps, times), prep_fm))
                       for steps in (steps_static, steps_light))
        return st, li, n_prep
    static_rule, light_rule, n_prep_rule = _rows_rule_step()              # record: 7.55e9 / 2.69e10, 2.3e8 / 8.3e8
    r2 = (da / hbarc) ** 2                                                # GeV^-2
    P0 = float(a.dipole_P0)
    G = float(a.kappa_over_T3) * float(a.dipole_T_GeV) ** 3 * r2 * float(a.static_colour_factor) / hbarc   # fm^-1
    p_static = tuple(P0 * math.exp(-G * t) for t in times)
    p_light = tuple(P0 * math.exp(-float(a.qhat_GeV2_fm) * t * r2 / 2) for t in times)
    vf = (float(a.thermal_variance_factor.lo), float(a.thermal_variance_factor.hi))   # (1, 1)
    n_T = int(a.n_dipole_T.value)

    def _floor_slope_ratio(rate, n_c):
        """J4: two-point log slope of P_s in the Markovian color description, P_s = P0 (1/N^2 + (1 - 1/N^2)
        e^{-Gamma t}), Gamma = rate N^2 / (N^2 - 1), over the |W|^2 slope 'rate' (= 2|Im V|, or q-hat r^2 / 2)."""
        f = 1 / n_c ** 2
        gam = rate / (1 - f)
        ps = [f + (1 - f) * math.exp(-gam * t) for t in times]
        return math.log(ps[0] / ps[1]) / (rate * (times[1] - times[0]))

    def _floor_curve(rate, n_c):
        """J4: P_s at the two times in the Markovian color description (|W|^2 slope 'rate')."""
        f = 1 / n_c ** 2
        gam = rate / (1 - f)
        return tuple(P0 * (f + (1 - f) * math.exp(-gam * t)) for t in times)

    def _dipole(rows_c, rows_evo, rows_rule, p, rate, n_c):
        # r26 v2 (J4): binomial on the measured P_s (color-floor curve); the rate is the P_s slope over
        # floor_slope_ratio, a fixed model factor, so the relative error of the rate is that of the slope.
        # q1: rows_c = (rows at c = 2, rows at c = 2 pi), each the two times; walls and t_shot span the c range
        p_s = _floor_curve(rate, n_c)
        n_bin = two_point_log_shots(p_s[0], p_s[1], delta_first)          # both times, one temperature, binomial
        n_w2 = two_point_log_shots(p[0], p[1], delta_first)               # record: r26 v1, the |W|^2 shape
        shots_T = (vf[0] * n_bin, vf[1] * n_bin)                          # (1, 1): no estimator penalty on a pure state
        per_pair = tuple(sum(0.5 * (r["t_shot"] * gate + t0) for r in rows) for rows in rows_c)   # mean of the two times
        wall_T = (n_bin * per_pair[0], n_bin * per_pair[1])               # one temperature at c = 2 .. 2 pi
        t_shot = (rows_c[0][0]["t_shot"], rows_c[1][-1]["t_shot"])        # (c = 2, 0.5 fm/c) .. (c = 2 pi, 1 fm/c)
        longest = rows_c[1][-1]
        return {"lq": longest["lq"], "t_shot": t_shot, "eps": (rows_c[0][0]["eps"], longest["eps"]),
                "t_shot_by_c": tuple(tuple(r["t_shot"] for r in rows) for rows in rows_c),
                "t_shot_evolution_only": tuple(r["t_shot"] for r in rows_evo),   # record: 2.35 / 6.66e9, 7.2e7 / 2.0e8
                "t_shot_prep_at_rule_step": tuple(r["t_shot"] for r in rows_rule),   # record, unpriced: 7.55e9 / 2.69e10, 2.3e8 / 8.3e8
                "n_prep_steps_at_rule_step": n_prep_rule[0] if n_c == 3 else n_prep_rule[1],   # record: (187, 725), (141, 547)
                "n_prep_steps": n_prep_q, "n_ramp_steps": n_ramp_q, "n_thermalization_steps": n_th_q,
                "t_step": longest["t_step"], "p": p, "p_s": p_s, "shots_per_T_binomial": n_bin, "shots_per_T": shots_T,
                "shots_per_T_W2": n_w2,                                                       # record: 615 / 250 (r26 v1)
                "wall_first_s": wall_T, "wall_scan_s": (n_T * wall_T[0], n_T * wall_T[1]),   # n_T = 3 temperatures, printed
                "eps_l_required": float(a.clean_shot_faults) / t_shot[1],   # 0.1 faults on the longest circuit, printed
                "t_shot_over_cap": (t_shot[0] / cap, t_shot[1] / cap),
                "t_shot_max_over_cap": t_shot[1] / cap, "lq_over_marker": longest["lq"] / float(a.rfi_lq_2033)}
    rate_light = float(a.qhat_GeV2_fm) * r2 / 2
    dip_static = _dipole(static_q, static_evo, static_rule, p_static, G, 3)
    dip_light = _dipole(light_q, light_evo, light_rule, p_light, rate_light, 2)
    dip_static["floor_slope_ratio"] = _floor_slope_ratio(G, 3)                        # J4, SU(3)
    dip_light["floor_slope_ratio"] = _floor_slope_ratio(rate_light, 2)               # J4, SU(2)
    # J4 v2: if 2 GeV^2/fm is the adjoint q-hat, the fundamental SU(2) source sees C_F/C_A = 3/8 of it
    cf_ca_su2 = (3 / 4) / 2
    p_s_adj = _floor_curve(rate_light * cf_ca_su2, 2)
    dip_light["shots_qhat_if_adjoint"] = two_point_log_shots(p_s_adj[0], p_s_adj[1], delta_first)   # 695
    dip_light["shots_qhat_if_adjoint_W2"] = two_point_log_shots(
        *(P0 * math.exp(-rate_light * cf_ca_su2 * t) for t in times), delta_first)                  # 587

    # ---- r27 (chapter-open pass 2026-10-05), derived here --------------------------------------------------
    # (q1: the r27 G8 E-rho-OQ penalty exports are retired with the method, ruling A, 2026-10-05.)
    links_static = 3 * math.prod(int(x) for x in a.static_dims.value)       # 81

    # J2: box effects in the trend. A Lund yo-yo (arxiv_2203_11601) with ends at +-p reaches full extension at
    # t = p / kappa with string length 2p / kappa; on the L_par a axis it wraps 2p / (kappa L_par a) times, so the
    # windings grow with p, so cancellation in the trend is not expected (box-effect size on B/M not derived).
    kap = float(a.lund_kappa_GeV_fm)
    L_ax_fm = int(a.L_par_priced.value) * a_fm_
    box_modes = tuple(k * p_unit[0] for k in range(1, n_mom + 1))          # 1.55, 3.10, 4.65 GeV
    yoyo_t_fm = tuple(p / kap for p in box_modes)                          # 1.55, 3.10, 4.65 fm/c
    yoyo_windings = tuple(2 * p / (kap * L_ax_fm) for p in box_modes)      # 3.9, 7.7 (7.749), 11.6
    # a second box in the same reduced model shares box mode k of L_par = 8 iff k L / 8 is an integer
    L8 = int(a.L_par_priced.value)
    shared = {L: tuple(k for k in range(1, n_mom + 1) if (k * L) % L8 == 0) for L in range(1, 2 * L8 + 1) if L != L8}
    L_all = min(L for L, ks in shared.items() if len(ks) == n_mom)          # 16: every trend point
    L_one = min(L for L, ks in shared.items() if ks)                        # 4: the 3.1 GeV point only
    lq_box_all = register(q_g, L_perp * L_all, D, n_c, n_stag, n_anc)[3]   # 1828
    lq_box_one = register(q_g, L_perp * L_one, D, n_c, n_stag, n_anc)[3]   # 532
    # price the one-point check at L_par = 4 as the headline circuit: same steps, own R-TOL eps, corrected shot time
    V_one = L_perp * L_one
    n_rot_step_one = synth_rotations_per_step(a, gk, V_one * D, V_one)
    shot_one, t_one = [], []
    for e, n in enumerate(n_steps):
        ax = _at_eps(a, rtol_eps(a, n * n_rot_step_one), n * n_rot_step_one)
        tt = n * _step(ax, V_one)[0]
        dd = n * depth_2033_step(ax, gauge_step(ax, gk, D)["t_per_rotation_fermion"], V_one * D, e)
        t_one.append(tt)
        shot_one.append(max(tt * gate, dd * REACTION_TIME_S) + t0)
    wall_box_one = (shots_pt[0] * shot_one[0], shots_pt[1] * shot_one[1])  # one configuration at the trend precision
    # smaller ruling (a) (2026-10-05, Henry: "do ... the smaller things"; Claude's recommended default): the L_par = 4
    # box-dependence check is part of the campaign. Campaign = the 3-momentum trend + the one-point check.
    shots_camp = (shots_trend[0] + shots_pt[0], shots_trend[1] + shots_pt[1])     # (1.27e4, 2.53e4)
    wall_camp = (wall_trend[0] + wall_box_one[0], wall_trend[1] + wall_box_one[1])   # 1.72-25.8 yr (was 1.48-22.1)
    camp_over_horizon = (wall_camp[0] / horizon_s, wall_camp[1] / horizon_s)     # (0.34, 5.2) (was 4.4)

    # J3: the N_h insertion. Hadamard test with the projector on both branches: (|0>|a> + |1>|b>)/sqrt2, then the
    # backward ramp (uncontrolled, already in every 2033 shot) and the label measurement; the estimator X 1_h has mean
    # Re<b|N_h|a> and variance ~ p_i (the bin's occupation). Per bin N_i = p_i / (eps A g_i)^2, g_i the bin's share of
    # the bilocal amplitude. Momentum sum rule: sum_h sum_i z_i g_{h,i} = 1, so sum_i z_i g_{pi,i} <= 1. Relative to
    # the printed N_z / (eps A)^2: p_i = g_i gives >= (sum sqrt z_i)^2 / N_z; p_i = 1 (variance bound) >= (sum z_i^(2/3))^3 / N_z.
    nz = int(a.n_z.value)
    z0, z1 = float(a.ff_z_range.lo), float(a.ff_z_range.hi)
    zc = tuple(z0 + (z1 - z0) * (i + 0.5) / nz for i in range(nz))
    ff_nh_factor = (sum(math.sqrt(z) for z in zc) ** 2 / nz,              # 4.70 (p_i = g_i)
                    sum(z ** (2 / 3) for z in zc) ** 3 / nz)                # 23.1 (variance bound)
    ff_nh_extra_t = 0.0                                                    # 0 for the projector; momentum transform unpriced

    # requirements row: 3 sampling configurations + 120 fragmentation calls + 2 dipole rows x 3 temperatures
    n_ff_calls = int(a.n_z.value) * int(a.n_tens.value) * int(a.n_boost.value)   # 120
    n_calls = n_mom + n_ff_calls + 2 * n_T                                # 129 -> '~130'
    n_calls_r23 = n_config + n_ff_calls + int(a.n_qhat_T.value)           # 129, r23

    # R9 exports for the era: the campaign run (one machine)
    # step rule: at the high end one of the three momenta runs the 575-step shot, so the exports take the campaign's
    # shot-averaged T and depth there (the low end, 30 a, is one circuit)
    t_camp_mean = (t_total[0], ((n_mom_ - 1) * t_total[1] + t_top) / n_mom_)
    d_camp_mean = (d_shot[0], ((n_mom_ - 1) * d_shot[1] + d_top) / n_mom_)
    dx = depth_exports(t_camp_mean, d_camp_mean, shots_trend, gate, t0)   # the trend run on 4x8 (the check is its own circuit)

    # systematics in quadrature
    syst = (quadrature(a.syst_discretization.lo, float(a.syst_volume), float(a.syst_subgroup)),
            quadrature(a.syst_discretization.hi, float(a.syst_volume), float(a.syst_subgroup)))   # (0.320, 0.364)
    total_err = (quadrature(float(a.eps_stat), syst[0]), quadrature(float(a.eps_stat), syst[1]))  # (0.377, 0.415)

    ramp_lo = ramp[0]      # breakdown at the low end of the ramp range (the first quoted total, 4.9e9)
    breakdown = (gauge_primitives(a, gk, D, n_links, n_evo, "evolution")
                 + fermion_primitives(a, gk, n_links, V, n_evo, "evolution"))
    if ramp_lo > 0:
        breakdown += (gauge_primitives(a, gk, D, n_links, ramp_lo, "ramp")
                      + fermion_primitives(a, gk, n_links, V, ramp_lo, "ramp"))
        if n_ramps > 1:    # r26 (J1): the backward ramp of the species readout map
            breakdown += (gauge_primitives(a, gk, D, n_links, (n_ramps - 1) * ramp_lo, "readout_ramp")
                          + fermion_primitives(a, gk, n_links, V, (n_ramps - 1) * ramp_lo, "readout_ramp"))
    inter = {
        "stage": STAGE,
        "q_G_S72x3": q_g, "q_G_S72x3_compiled": g.link_qubits,            # 9, 9
        "q_G_S72x3_printed_before_R1": q_g_printed_before,                # 8
        "ceil_log2_order_S72x3": math.ceil(math.log2(g.order)),           # 8: 'the dense floor, but no circuit has been built on it'
        "V_uncut": V_uncut, "lq_gauge_uncut": gauge_u, "lq_fermion_uncut": fermion_u, "lq_anc_2033": n_anc,   # 64, 1152, 576, 100
        "lq_uncut_2033": lq_uncut,                                        # 1828
        "lq_uncut_over_marker": lq_uncut / float(a.rfi_lq_2033),          # 1.83 -> '1.8x the register limit'
        "V_cut": V_cut, "lq_gauge_cut": gauge_cut, "lq_fermion_cut": fermion_cut,   # (32,40), (576,720), (288,360)
        "lq_cut_2033": lq_cut,                                            # (964, 1180): sensitivity; '1180 LQ ... at L_par=10'
        "lq_priced_2033": lq,                                             # 964, the box (V = 32)
        "lq_priced_2033_over_marker": lq / float(a.rfi_lq_2033),          # 0.964 -> 'at the 1000-LQ marker'
        "lq_cut_2033_over_marker": (lq_cut[0] / float(a.rfi_lq_2033), lq_cut[1] / float(a.rfi_lq_2033)),   # (0.96, 1.18)
        "lq_cut_2033_mid": 0.5 * (lq_cut[0] + lq_cut[1]),                 # 1072
        "lq_cut_2033_sampling_ancilla_only": lq_cut_samp,                 # (872, 1088): readout (a) without the 100 (record; not printed)
        "legacy_lq_cut_2033": lq_cut_legacy,                              # (900, 1100), the pre-ruling print at q_G = 8
        "V_priced_2033": V, "n_links_2033": n_links, "n_plaq_2033": n_plaq,   # 32, 64, 32
        # conventions
        "eps_syn": float(a.eps_syn),                                      # 1e-2 (R-TOL)
        "eps_rot": eps_end,                                               # (1.104e-5, 8.187e-6): each ramp end its own
        "eps_rot_round_e": float(a_round_e.eps_rot),                      # 1e-3, superseded
        "t_per_rotation_papers": (gs["t_per_rotation_papers"], gs_hi["t_per_rotation_papers"]),   # legacy record
        "t_per_rotation_report": (gs["t_per_rotation_fermion"], gs_hi["t_per_rotation_fermion"]),   # hop and mass
        "t_per_toffoli": gs["t_per_toffoli"],                             # 7
        # per link, per step (ruling R2)
        "magnetic_t_per_link_2033": gs["magnetic"],                       # 9.44e3 -> '9.4e3'
        "electric_floor_t_per_link_2033": gs["electric"],                 # 3.47e4 -> '3.5e4' (a floor: no Sigma(72x3) FFT)
        "electric_dense_t_per_link_2033": gs["electric_dense"],           # 9.5e6 (full fit)
        "fft_compiled_2033": gs["fft_compiled"],                          # False
        "fourier_key_2033": gs["fourier_key"],                            # 'U_FFT' (the stated lower bound)
        "t_per_link_2033": gs["per_link"],                                # 4.412e4 (hi end 4.469e4) -> '4.4--4.5e4'
        "electric_floor_share_of_step_2033": n_links * gs["electric"] / t_step,   # 0.15
        "t_gauge_per_step_2033": t_gauge_step,                            # 2.82e6 -> '2.8e6'
        # hop (r17: FermionPrimitives_unpub color-diagonalized staggered hop; r19: share=link + undo, E21)
        "hop_toffoli_per_link_2033": hl["toffoli"],                       # 3,672 (r17-r18, share=draft: 14,314)
        "hop_t_direct_per_link_2033": hl["t_direct"],                     # 360
        "hop_rotations_per_link_2033": hl["n_rot"],                       # 2933
        "hop_items_2033": {it["name"]: (it["toffoli"], it["t_direct"], it["n_rot"], it["status"].name)
                           for it in hl["items"]},
        "hop_estimated_items_2033": tuple(it["name"] for it in hl["items"]
                                          if it["status"] is not CircuitStatus.COMPILED),
        "t_hop_per_link_2033": hl["t"],                                   # 1.116e5 -> '1.1e5'
        "t_hop_per_link_2033_hi_end": hl_hi["t"],                         # 1.131e5
        "hop_rotation_share_2033": hl["t_rotations"] / hl["t"],           # 0.77
        "hop_variant_t_per_link_2033": {k: v["t"] for k, v in hop_var.items()},   # share_draft, no_undo, draft_mcx, share_draft_no_undo
        "rhop_t_per_link_2033": hl["rhop"]["t"],                          # retired R-HOP (slope-only) at this eps
        "hop_over_rhop_2033": hl["t"] / hl["rhop"]["t"],                  # 16.0
        "pre_r17_n_rot_per_step_2033": n_rot_step_pre,                    # 74,592
        "pre_r17_t_total_2033": t_total_pre,                              # (2.768e9, 5.106e9), the R-TOL print
        "t_hop_per_step_2033": t_hop_step,                                # 7.15e6
        "t_mass_per_step_2033": t_mass_step,                              # 5303
        "hopping_share_of_step_2033": t_hop_step / t_step,                # 0.716 -> '72%'
        "t_per_step_2033": t_step_end,                                    # (9.754e6, 1.004e7) -> '1.0e7'
        "t_per_step_2033_round_e": t_step_round_e,                        # 1.95852e6 at 1e-3 (record)
        "retired_hop_insertion_2033": retired_rest,                       # 3.4e5
        "hop_over_retired_2033": t_hop_step / retired_rest,               # 21.0
        "applied_stage_t_per_step_2033": applied_step,                    # 2.136e6, before round D
        "insertion_t_S72x3": t_ins,                                       # 1.806e4 -> '~2e4'
        "insertion_umul_count": 2 * (int(a.insertion_links.value) - 1) + int(a.moves_per_field.value),   # 16
        "insertion_per_shot_b": (n_ins[0] * t_ins, n_ins[1] * t_ins),     # (1.8e4, 3.6e4)
        "insertion_share_of_shot_b": (n_ins[0] * t_ins / t_total[1], n_ins[1] * t_ins / t_total[0]),   # (4e-6, 1.5e-5)
        "synth_rotations_per_step_2033": n_rot_step,                      # 7.46e4
        "synth_rotations_per_shot_2033": n_rot_shot,                      # (8.2e7, 1.5e8)
        "synthesis_error_2033": synth_err,                                # (1e-2, 1e-2) = eps_syn (R-TOL)
        "synthesis_error_2033_round_e": (synthesis_error(a_round_e, n_rot_shot[0]),
                                         synthesis_error(a_round_e, n_rot_shot[1])),   # (82, 149) at 1e-3, record
        "eps_rot_for_budget_2033": eps_budget,                            # (2.6e-5, 3.5e-5) record: the whole 0.1
        "n_trot_2033": n_evo_full,                                        # 1e3: the r23 t = 100 a window (record)
        "n_evo_2033": n_evo_end,                                          # (300, 500): R4, 30-50 a at 0.1 a
        "n_steps_2033": n_steps,                                          # (400, 1500) with the 1e2-1e3 ramp
        "n_steps_2033_r23": n_steps_full,                                 # (1100, 2000)
        "t_total_2033_r23": t_total_full,                                 # (1.097e10, 2.021e10): the r23 headline
        "window_over_a": (float(a.window_over_a.lo), float(a.window_over_a.hi)),   # (30, 50)
        "window_fm": window_fm,                                           # (3, 5) fm/c
        "window_over_L_par": crossings_window,                            # (3.75, 6.25) box lengths
        "t_evolution_2033": t_evo_end,                                    # (2.926e9, 5.021e9) -> '2.9--5.0e9'
        "t_evolution_2033_over_cap": t_evo / cap,                         # 2.93: the evolution alone is over the limit
        "ramp_steps_2033": ramp,                                          # (1e2, 1e3)
        "t_ramp_2033": t_ramp,                                            # (9.75e8, 1.004e10) -> '9.8e8--1.0e10'
        "t_total_2033": t_total,                                          # (3.902e9, 1.506e10) -> '3.9e9--1.5e10'
        "t_total_2033_over_cap": (t_total[0] / cap, t_total[1] / cap),    # (3.90, 15.1) -> '3.9--15x'
        # r26 referee fixes (J1 readout map, J2 source and volume)
        "n_ramps_2033": n_ramps,                                          # 2: prep + backward readout ramp (J1)
        "n_steps_2033_r25": n_steps_r25,                                  # (400, 1500), record
        "t_total_2033_r25": t_total_r25,                                  # (3.902e9, 1.506e10), the r25 headline (record)
        "t_ramps_2033": t_ramps,                                          # both ramps
        "t_readout_fourier_2033": t_readout_fourier,                      # one transform per link, not added
        "t_readout_fourier_share": t_readout_fourier / t_total[0],        # < 1e-3
        "t_source_per_step_2033": t_source_per_step,                      # one link's hop while the source is on
        "t_source_share_of_step": t_source_per_step / t_step,             # ~1% (r26 v1 record)
        "t_source_wavepacket_2033": t_source_wavepacket,                  # 1.45e5 T per insertion (J2 v2)
        "t_source_wavepacket_share_of_step": t_source_wavepacket / t_step,   # 0.015 -> '< 2% of a step'
        "recollision_time_fm": recollision_fm,                            # 0.4 fm/c on the 0.8 fm axis
        "L_par_for_separation": L_par_sep,                                # 20 sites (2 fm)
        "lq_separation_volume": lq_sep,                                   # 2260
        # r27 (2026-10-05)
        "source_box_modes_GeV": box_modes,                                # (1.55, 3.10, 4.65)
        "yoyo_full_extension_fm": yoyo_t_fm,                              # (1.55, 3.10, 4.65) fm/c at kappa = 1 GeV/fm
        "yoyo_windings_at_full_extension": yoyo_windings,                 # (3.9, 7.7, 11.6)
        "box_check_L_par_all_modes": L_all,                               # 16
        "box_check_L_par_one_mode": L_one,                                # 4 (shares 3.1 GeV)
        "box_check_shared_modes_one": shared[L_one],                      # (2,)
        "lq_box_check_all_modes": lq_box_all,                             # 1828
        "lq_box_check_one_mode": lq_box_one,                              # 532
        "t_total_box_check_one_mode": tuple(t_one),                       # (2.4e9, 1.25e10)
        "wall_per_shot_box_check_one_mode_s": tuple(shot_one),            # (2.4e3, 1.8e4)
        "wall_box_check_one_mode_s": wall_box_one,
        "wall_box_check_one_mode_yr": (wall_box_one[0] / yr, wall_box_one[1] / yr),   # (0.24, 3.7)
        "box_check_over_campaign": (wall_box_one[0] / wall_trend[0], wall_box_one[1] / wall_trend[1]),   # vs the trend
        "box_check_in_campaign": True,                                    # smaller ruling (a), 2026-10-05
        "ff_z_bin_centres": zc,
        "ff_nh_shot_factor": ff_nh_factor,                                # (4.70, 23.1): lower bounds under each p_i reading
        "ff_nh_extra_t_per_shot": ff_nh_extra_t,                          # 0 for the projector; momentum transform unpriced
        "t_total_2033_round_e": t_total_round_e,                          # (2.154e9, 3.917e9) at 1e-3, record
        "steps_in_cap_2033": steps_in_cap,                                # 100 steps (low end's eps; not printed)
        "t_per_step_2033_at_V_hi": t_step_hi,                             # at V = 40, low ramp end (not in the chapter)
        "t_total_2033_at_V_hi": t_total_hi,                               # each end at its own R-TOL eps
        "t_total_2033_qubitized": (t_total[0] / float(a.qubitized_trim), t_total[1] / float(a.qubitized_trim)),   # (3.90e8, 1.51e9)
        "gap_cut_factor_2033": t_step / float(a.t_step_gap_target),       # 1.0e4
        "report_convention_t_gauge_per_step_2033": n_links * gs["per_link_report"],
        "magnetic_report_over_papers_2033": gs["magnetic_report"] / gs["magnetic"],     # 1.0
        "electric_report_over_papers_2033": gs["electric_report"] / gs["electric"],     # 1.4-1.6
        # the pre-ruling text, for the record
        "legacy_t_per_plaquette_S72x3": t_plaq_legacy,                    # 5e3
        "legacy_t_per_step_2033": legacy_step,                            # 5e5
        "legacy_t_total_2033": legacy_total,                              # (5.5e8, 1e9)
        "gauge_step_over_legacy_plaquette_step": t_gauge_step / (n_plaq * t_plaq_legacy),   # 11.2
        # the unpublished draft's hopping, NOT the headline
        "draft_hop_t_per_link_S72x3": draft_hop_link,                     # 5.54e4
        "draft_hop_t_per_step_2033": draft_hop_step,                      # 3.55e6 (one field)
        "draft_hop_t_per_step_2033_3fields": draft_hop_step_3fields,      # 1.06e7
        "draft_hop_per_field_over_rule": draft_hop_link / (hl["rhop"]["t"] / n_stag),   # printed C^S / R-HOP per field
        "draft_printed_3fields_over_model_2033": draft_hop_step_3fields / t_hop_step,   # 1.56: 3 x printed C^S vs the model
        "draft_t_per_step_2033": draft_step,                              # 5.34e6
        "draft_t_total_2033": draft_total,                                # (5.9e9, 1.07e10)
        # kinematics
        "dt_over_a_2033": dt_over_a,                                      # 0.1 = t_had / N_Trot (derived, not printed)
        "t_had_fm": t_had_fm,                                             # 10
        "n_trot_eq_range": (float(a.n_trot_eq.lo), float(a.n_trot_eq.hi)),   # (1e2, 1e3)
        "momentum_unit_GeV_at_Lpar_8_10_16": p_unit,                      # (1.55, 1.24, 0.77)
        "source_momentum_in_units_Lpar_8": float(a.source_momentum_GeV) / p_unit[0],   # 1.29: 2 GeV is not a mode at L_par = 8
        "source_momentum_times_a": float(a.source_momentum_GeV) * a_fm / hbarc,        # 1.01
        "source_momenta_box_modes_GeV": tuple(k * p_unit[0] for k in (1, 2, 3)),       # (1.55, 3.10, 4.65), J2 v2
        "top_box_mode_times_a": 3 * p_unit[0] * a_fm / hbarc,                          # 2.36: lattice-artifact region
        "t_had_over_L_par_cut": crossings,                                # (10, 12.5) lattice lengths in the window
        "inverse_lambda_qcd_fm": inv_lambda_fm,                           # (0.66, 0.99): workflow step 4's 't ~ 1/Lambda_QCD'
        "t_had_over_inverse_lambda_qcd": (t_had_fm / inv_lambda_fm[1], t_had_fm / inv_lambda_fm[0]),   # (10, 15)
        # logical faults
        "faults_per_shot_2033_at_rfi_eps_l": faults_rfi,                  # (49.0, 261.6) -> '49--260' (step rule; r26 250)
        "faults_per_shot_2033_at_1e-9": (t_total[0] * 1e-9, t_total[1] * 1e-9),   # (3.90, 15.1)
        "eps_l_required_2033": eps_l_req,                                 # (3.82e-12, 2.04e-11) -> '3.8e-12--2.0e-11' (step rule)
        # R9/R10 depth (factory.json schedule P) and the corrected per-shot time
        "t_depth_per_step_2033": {f"{('lo', 'hi')[e]}_circuit_{('lo', 'hi')[k]}_est": v for (e, k), v in d_step.items()},
        "f_star_per_circuit_2033": f_star_circuit,                        # ((7.0, 18.5), (7.0, 18.5)) per ramp end
        "wall_per_shot_serial_2033_s": wall_shot,                         # N_T x 1 us + 0.1 ms
        "wall_per_shot_2033_s": shot_s,                                   # max(N_T t_gate, D_T t_r) + t0 -> printed
        "depth_penalty_2033": depth_penalty,                              # (1.0, 1.43)
        "shot_overhead_s": t0,                                            # 1e-4
        # R1 first result: one configuration at 30%
        "delta_first_2033": delta_first,                                  # 0.30
        "shots_first_2033": shots_first,                                  # (122, 244)
        "wall_first_2033_s": wall_first,
        "wall_first_2033_days": (wall_first[0] / day, wall_first[1] / day),
        # R2 campaign: the 3-momentum trend at 3 sigma, one spacing
        "ratio_shots_at_10pct": ratio_10,                                 # (1.1e3, 2.2e3), pair clustering (R2)
        "trend_delta_per_point": delta_tr,                                # 0.0589
        "trend_sigma_at_delta": sigma_at,                                 # {0.05: 3.5, 0.10: 1.8, 0.30: 0.59}
        "shots_per_configuration_campaign": shots_pt,                     # (3.17e3, 6.34e3)
        "shots_trend_2033": shots_trend,                                  # (9.5e3, 1.9e4)
        "wall_trend_2033_s": wall_trend,
        "wall_trend_2033_yr": (wall_trend[0] / yr, wall_trend[1] / yr),    # (1.48, 22.1), the r26/r27 campaign
        "trend_over_horizon": trend_over_horizon,                         # (0.30, 4.4)
        "shots_campaign_2033": shots_camp,                                # (1.27e4, 2.53e4): trend + L_par = 4 check
        "wall_campaign_2033_s": wall_camp,
        "wall_campaign_2033_yr": (wall_camp[0] / yr, wall_camp[1] / yr),  # (1.72, 25.8) -> '1.7--26 yr'
        "campaign_over_horizon": camp_over_horizon,                       # (0.34, 5.16) -> '5.2x'
        "configurations_in_horizon": configs_in_horizon,
        "wall_campaign_two_spacings_yr": (wall_camp_2sp[0] / yr, wall_camp_2sp[1] / yr),
        "lq_second_spacing_fixed_volume": lq_fixed_volume_2sp,           # 3556
        # R7 dipole rows (first results at 30%, one temperature; the 3-temperature scan beside)
        "dipole_static": dip_static,
        # s648 ruling (5): the static row on Sigma(216x3) / H_I; the species row kept on Sigma(72x3) / H_KS because at
        # 11 qubits/link its register is 1092 LQ > 1000 (V = 32; ~1.1k caveat band), and a switch would reuse Ch. 5's
        # Sigma(72x3) placeholder hop (hop_group_2033; FermionPrimitives_unpub stops at Sigma(72x3)), 72% of the step
        "dipole_static_group": (gk_dip, ham_dip),
        "s648_species_lq_if_switched": register(GROUPS["S216x3"].link_qubits, V, D, n_c, n_stag, n_anc)[3],   # 1092
        "s648_static_lq_if_kept_S72x3": GROUPS["S72x3"].link_qubits * links_static + int(a.dipole_n_anc.value)
                                         + int(a.static_src_qubits.value),                                     # 833
        "dipole_lightlike": dip_light,
        "dipole_steps": steps_d,                                          # (50, 100): the 0.01 fm/c grid
        "dipole_steps_rule": {"static": steps_static, "light": steps_light},   # (72, 203), (64, 181): step rule
        "dipole_sd_err_rule": {k: tuple(x[2] for x in v["rows"]) for k, v in dr.items()},   # each <= 0.1
        "dipole_sd_lambda": {k: v["lambda"] for k, v in dr.items()},
        "step_rule_top_mode": top,                                        # n 575 (was 500), phase 0.0997 (was 0.132)
        "t_shot_top_mode_2033": t_top,                                    # 2.61e10
        "wall_shot_top_mode_2033_s": shot_s_top,
        "wall_shot_campaign_mean_hi_2033_s": shot_s_campaign_hi,
        "t_campaign_mean_2033": t_camp_mean, "d_campaign_mean_2033": d_camp_mean,   # R9 export inputs (step rule)
        # the r23 sampling accounting, record
        "shots_sampling_from_eq": shots_a,                                # (1e3, 4e3), r23
        "shots_sampling_printed": shots_a_printed,                        # (1e3, 4e3), r23
        "wall_per_shot_2033_r23_s": wall_shot_full,                       # (1.097e4, 2.021e4)
        "wall_sampling_s_r23": wall_a,                                    # (1.097e7, 8.083e7), one configuration
        "campaign_sampling_yr_at_10pct_r23": camp_a_10,                   # (2.09, 3.84)
        "campaign_sampling_yr_at_5pct_r23": camp_a_5,                     # (8.34, 15.37)
        # readout (b), FCC clock
        "a_inv2": a_inv2,                                                 # (25, 400)
        "shots_ff_base": shots_b_base,                                    # 3e3
        "shots_ff_from_eq": shots_b,                                      # (7.5e4, 1.2e6)
        "shots_ff_printed": shots_b_printed,                              # (7.5e4, 1.2e6)
        "wall_ff_s": wall_b,
        "wall_ff_yr": (wall_b[0] / yr, wall_b[1] / yr),
        "ff_cost_reduction_for_horizon": ff_reduction,                    # serial years / 5
        "ff_shots_in_horizon": ff_shots_in_horizon,
        "ff_shots_per_call": shots_per_ff_call,                           # (625, 1e4) = eps^-2 A^-2 (not printed)
        "ff_calls_in_horizon": ff_calls_in_horizon,
        # campaign bookkeeping
        "n_sampling_configurations": n_mom,                               # 3 (one spacing)
        "n_sampling_configurations_r23": n_config,                        # 6
        "n_ff_calls": n_ff_calls,                                         # 120
        "n_subroutine_calls": n_calls,                                    # 129 -> '~130'
        "n_subroutine_calls_r23": n_calls_r23,                            # 129
        "campaign_horizon_yr": horizon,                                   # 5
        # R9 exports (common.depth_exports, the campaign run)
        **{k: dx[k] for k in ("t_per_shot", "t_depth_per_shot", "f_star", "floor_wall_s", "factories_for_1yr")},
        "wall_serial_s": dx["wall_serial_s"],
        "baseline_ok": dx["baseline_ok"],
        "wall_first_result_s": wall_first,
        "wall_campaign_s": wall_camp,
        # systematics
        "syst_quadrature": syst,                                          # (0.320, 0.364) -> '32--36%'
        "total_error_quadrature": total_err,                              # (0.377, 0.415) -> '38--42%' (R4)
        # clocks and utility (not resource numbers)
        "z_sample_ratio_fcc_over_lep": float(a.n_z_fcc) / float(a.n_z_lep),                  # 3.5e5
        "sqrt_z_sample_ratio": math.sqrt(float(a.n_z_fcc) / float(a.n_z_lep)),               # 594; the chapter cites ~300
        "utility_musd_per_instance": float(a.utility_share) * float(a.cms_capital_musd) / float(a.n_instances_utility),   # 2
    }
    notes = (
        "Step rule (2026-10-04, Claude's decision on G7, NEEDS_AUTHOR closed): the state-dependent estimate at eps = 0.1 "
        "sets the step. The window evolves the ramp-prepared vacuum (an eigenstate; its error is the bounded free-field "
        f"0.04-0.05), so the rule binds on the injected pair: the 4.65 GeV point at 50 a runs {top['n']} steps (phase "
        f"{top['phase']:.4f}), {t_top:.4g} T per shot, the box's upper end. The thermal dipole rows take the literal "
        f"estimate: {steps_static} and {steps_light} steps.",
        f"Register (eq:Nq_collider) at q_G = {q_g} (R1): uncut 4x16 {gauge_u} + {fermion_u} + {n_anc} = {lq_uncut}; "
        f"cut to V = {V_cut[0]}-{V_cut[1]}: {lq_cut[0]}-{lq_cut[1]}, {lq_cut[0] / float(a.rfi_lq_2033):.2f}-"
        f"{lq_cut[1] / float(a.rfi_lq_2033):.2f}x the 1000-LQ marker. The box quotes {lq} at the priced V = {V} "
        f"(ruling R11 option (A)). Before R1 the chapter printed {lq_cut_legacy[0]}-{lq_cut_legacy[1]}.",
        f"Synthesis (R-TOL, randomized): {n_rot_step:.0f} rotations per step x {n_steps[0]:.0f}-{n_steps[1]:.0f} steps = "
        f"{n_rot_shot[0]:.4g}-{n_rot_shot[1]:.4g} per shot; eps_rot = {eps_end[0]:.4g}-{eps_end[1]:.4g}, "
        f"{gs['t_per_rotation']:.3f}-{gs_hi['t_per_rotation']:.3f} T per gauge rotation, {gs['t_per_rotation_fermion']:.3f}-"
        f"{gs_hi['t_per_rotation_fermion']:.3f} per hop/mass rotation (full fit); N eps^2 = {float(a.eps_syn):g} at each end. "
        f"Pre-r17 (R-HOP, slope-only): {t_total_pre[0]:.4g}-{t_total_pre[1]:.4g} T.",
        f"Gauge terms (R2 tables, E20 full fit, low end's eps = {float(a.eps_rot):.4g}): {gs['magnetic']:.4g} magnetic + "
        f"{gs['electric']:.4g} electric = {gs['per_link']:.4g} T per link per step, {t_gauge_step:.4g} T per step on "
        f"{n_links} links. The electric term is the paper's stated lower bound for a fast transform: a FLOOR "
        f"(the dense transform is {gs['electric_dense']:.3g} T per link).",
        f"Hopping (r17, FermionPrimitives_unpub color-diagonalized hop for {n_stag} fields, groups.hop_link_cost; "
        f"share='{a.hop_share.value}', E21): {hl['toffoli']} Toffoli + {hl['t_direct']} T + {hl['n_rot']} rotations = {hl['t']:.4g} T per link, "
        f"{t_hop_step:.4g} per step ({100 * t_hop_step / t_step:.0f}% of the step), {hl['t'] / hl['rhop']['t']:.1f}x "
        f"the retired R-HOP; estimated pieces {', '.join(it['name'] for it in hl['items'] if it['status'] is not CircuitStatus.COMPILED)}. "
        f"Mass {t_mass_step:.4g}. Sensitivity: share='draft' (per-application recompute, the r17-r18 headline) "
        f"{hop_var['share_draft']['t']:.4g} T per link.",
        f"Current insertion is not a per-step term: readout (a) has none; readout (b) adds {n_ins[0]:.0f}-{n_ins[1]:.0f} x "
        f"{t_ins:.4g} T per shot (derived here).",
        f"Per shot (R4: window {a.window_over_a.lo:g}-{a.window_over_a.hi:g} a = {n_evo_end[0]:.0f}-{n_evo_end[1]:.0f} "
        f"evolution steps, plus the {ramp[0]:.0f}-{ramp[1]:.0f}-step ramp x {n_ramps}): {t_total[0]:.4g}-{t_total[1]:.4g} T, "
        f"{t_total[0] / cap:.1f}-{t_total[1] / cap:.1f}x the 1e9 limit. The r23 window (t = 100 a) gave "
        f"{t_total_full[0]:.4g}-{t_total_full[1]:.4g} T.",
        f"T and the box register are at V = {V}; at V = {V_hi} ({lq_cut[1]} LQ) the same rules give "
        f"{t_total_hi[0]:.3g}-{t_total_hi[1]:.3g} T (app10:95 parenthetical, not a box number).",
        f"T-depth (factory.json schedule P, hop links one at a time): {d_step[(0, 0)]:.3g}-{d_step[(1, 1)]:.3g} layers per "
        f"step, F* = {f_star_circuit[0][0]:.1f}-{f_star_circuit[0][1]:.1f} per circuit. Below 10 at the high estimate, so "
        f"each shot takes max(N_T t_gate, D_T t_r) + t0 = {shot_s[0]:.4g}-{shot_s[1]:.4g} s ({depth_penalty[1]:.2f}x "
        "the serial time at the high end).",
        f"First result (R1): one configuration at {100 * delta_first:.0f}%, {shots_first[0]:.0f}-{shots_first[1]:.0f} shots "
        f"(R2 pair clustering), {wall_first[0] / day:.1f}-{wall_first[1] / day:.0f} d on one machine.",
        f"Campaign (R2): the {n_mom}-momentum trend at {float(a.trend_sigma):g} sigma ({100 * delta_tr:.1f}% per point "
        f"for a {100 * float(a.trend_size):.0f}% trend), {shots_pt[0]:.3g}-{shots_pt[1]:.3g} shots per configuration, "
        f"{wall_trend[0] / yr:.2f}-{wall_trend[1] / yr:.1f} yr on one machine; with the L_par = 4 box check "
        f"(ruling 2026-10-05) {wall_camp[0] / yr:.2f}-{wall_camp[1] / yr:.1f} yr, {camp_over_horizon[1]:.1f}x the "
        f"horizon at the top; {configs_in_horizon[0]:.2f}-"
        f"{configs_in_horizon[1]:.1f} trend configurations fit in {horizon:.0f} yr. A second spacing on the same lattice "
        f"triples the trend; at fixed volume it needs {lq_fixed_volume_2sp} LQ.",
        f"Dipole rows (R7; q1 quench medium: {n_ramp_q} ramp + {n_th_q[0]}-{n_th_q[1]} thermalization steps, c = "
        f"{c_range[0]:g}-{c_range[1]:.3g}): static Sigma(216x3) H_I 3^3 {dip_static['lq']} LQ, {dip_static['t_shot'][0]:.3g}-"
        f"{dip_static['t_shot'][1]:.3g} T per shot ({dip_static['t_shot_evolution_only'][0]:.3g}/"
        f"{dip_static['t_shot_evolution_only'][1]:.3g} of it the dipole evolution at t = {times[0]:g}/{times[1]:g} fm/c), "
        f"{dip_static['t_shot_over_cap'][0]:.1f}-{dip_static['t_shot_over_cap'][1]:.1f}x the 1e9-T limit, "
        f"{dip_static['shots_per_T_binomial']:.0f} shots per temperature at 30% (binomial), {dip_static['wall_first_s'][0] / day:.3g}-"
        f"{dip_static['wall_first_s'][1] / day:.3g} d; light-like 2O 3x3x5 {dip_light['lq']} LQ, "
        f"{dip_light['t_shot'][0]:.3g}-{dip_light['t_shot'][1]:.3g} T (fits), {dip_light['shots_per_T_binomial']:.0f} shots, "
        f"{dip_light['wall_first_s'][0] / day:.3g}-{dip_light['wall_first_s'][1] / day:.3g} d; three temperatures "
        f"{dip_static['wall_scan_s'][0] / day:.3g}-{dip_static['wall_scan_s'][1] / day:.3g} and "
        f"{dip_light['wall_scan_s'][0] / day:.3g}-{dip_light['wall_scan_s'][1] / day:.3g} d.",
        f"Readout (b): {shots_b_base:.0f} x A^-2 = {shots_b[0]:.2g}-{shots_b[1]:.2g} shots; wall time "
        f"{wall_b[0]:.2g}-{wall_b[1]:.2g} s = {wall_b[0] / yr:.1f}-{wall_b[1] / yr:.0f} yr on one machine. "
        f"{horizon:.0f} years hold {ff_shots_in_horizon[0]:.2g}-{ff_shots_in_horizon[1]:.2g} shots, so the instance needs a "
        f"{ff_reduction[0]:.1f}-{ff_reduction[1]:.0f}x cut in shots x per-shot T (algorithmic / co-design; not a machine count).",
        "Result.shots and Result.wall_time_s are the campaign (R2).",
        f"Record: the draft's printed C^S (FP resources.tex:40) is {draft_hop_link:.3g} T per link per field, "
        f"{draft_hop_step_3fields:.3g} T per step for 3 fields, {draft_hop_step_3fields / t_hop_step:.2f}x the model "
        "(which adds the frame undo, shares the g-only pieces, prices MBU ladders and the full fit).",
        f"Readout map (r26, J1): N_h = U Pi_h U^dagger, the {ramp[0]:.0f}-{ramp[1]:.0f}-step ramp run backwards on every "
        f"shot ({n_ramps} ramps, {t_ramps[0]:.3g}-{t_ramps[1]:.3g} T); one Fourier layer to the electric basis "
        f"{t_readout_fourier:.3g} T ({t_readout_fourier / t_total[0]:.1e} of the shot, not added). The r25 circuit "
        f"without it: {t_total_r25[0]:.4g}-{t_total_r25[1]:.4g} T.",
        f"Source (r26, J2): gauge-invariant bilinear, {t_source_per_step:.3g} T per step while on "
        f"({100 * t_source_per_step / t_step:.1f}% of a step, not added); the pair recollides after {recollision_fm:g} fm/c; "
        f"separation through 1/Lambda needs L_par = {L_par_sep} ({lq_sep} LQ).",
        "Breakdown is at the low end of the ramp range (1e2 steps, prep and readout ramps).",
    )
    inter.update(trotter_check_2033(a))
    return Result(era="2033", lq=(lq, lq), hard_ops=hard_ops_2033, breakdown=breakdown, intermediates=inter,
                  shots=shots_camp, wall_time_s=wall_camp, epsilon_l=eps_l_req, notes=notes)


# --------------------------------------------------------------------------- #
# codesign row: the 3+1D extension of the plug-in prose (app10:90), not a box
# --------------------------------------------------------------------------- #

def _model_3d_stretch(a: Assumptions) -> Result:
    gk = a.group_2033.value
    g = GROUPS[gk]
    q_g = link_width(a, gk)                                               # 9 (R1)
    D = 3
    L = int(a.L_s_3d.value)
    V = L ** D                                                            # 64
    n_links = V * D                                                       # 192
    n_plaq = n_plaquettes(V, D)                                           # 192 = 3V
    n_c, n_stag, n_anc = int(a.n_c.value), int(a.n_stag.value), int(a.n_anc_2033.value)
    gauge, fermion, anc, lq = register(q_g, V, D, n_c, n_stag, n_anc)     # 1728, 576, 100, 2404
    lq_legacy = register(int(g.link_qubits_chapter.lo), V, D, n_c, n_stag, n_anc)[3]   # 2212

    n_evo = float(a.n_trot_3d)
    # ruling R-TOL: this circuit's rotation count (independent of eps) sets its tolerance
    n_rot_shot = n_evo * synth_rotations_per_step(a, gk, n_links, V)     # 7.859e8
    a_round_e = a                                                         # eps = 1e-3, round E: record only
    eps = rtol_eps(a, n_rot_shot)                                         # 3.567e-6
    a = _at_eps(a, eps, n_rot_shot)
    gs = gauge_step(a, gk, D)
    hl, ms = hop_link(a, gk), mass_site(a)
    t_gauge_step = n_links * gs["per_link"]                               # 1.047e7
    t_hop_step = n_links * hl["t"]                                        # 192 x 1.141e5 = 2.19e7
    t_mass_step = V * ms["t"]                                             # 64 x 169.0 = 1.08e4
    t_step = t_gauge_step + t_hop_step + t_mass_step                      # 3.238e7
    t_total = n_evo * t_step                                              # 3.238e10
    t_total_round_e = n_evo * _step_slope_only(a_round_e, gk, D, n_links, V)   # 7.673e9 at 1e-3, record
    n_rot_pre = n_evo * _n_rot_step_pre_r17(a_round_e, gk, n_links, V)  # 2.2432e8
    t_total_pre = n_evo * _step_slope_only(_at_eps(a_round_e, rtol_eps(a_round_e, n_rot_pre), n_rot_pre),
                                           gk, D, n_links, V)            # 9.537e9, the R-TOL print
    cap = float(a.rfi_t_2033)

    retired_rest = float(a.retired_hop_insertion_3d)                                                # 2.04e6 (RETIRED)
    legacy_step = n_plaq * float(g.plaquette_t) + retired_rest                                      # 3e6
    applied_step = t_gauge_step + retired_rest                                                      # 9.23e6 (APPLY stage)

    l2 = math.log2(1.0 / float(a.eps_rot))
    draft_hop_link = float(a.draft_hop_S72x3_t_const) + float(a.draft_hop_S72x3_t_log) * l2
    draft_hop_step = n_links * draft_hop_link                             # 1.06e7

    breakdown = (gauge_primitives(a, gk, D, n_links, n_evo, "evolution")
                 + fermion_primitives(a, gk, n_links, V, n_evo, "evolution"))
    eps_l_req = float(a.clean_shot_faults) / t_total
    inter = {
        "stage": STAGE,
        "V_3d": V, "n_links_3d": n_links, "n_plaq_3d": n_plaq,            # 64, 192, 192
        "lq_gauge_3d": gauge, "lq_fermion_3d": fermion, "lq_anc_3d": anc, # 1728, 576, 100
        "lq_3d": lq,                                                      # 2404 -> '~2400'
        "lq_3d_over_marker": lq / float(a.rfi_lq_2033),                   # 2.4
        "legacy_lq_3d": lq_legacy,                                        # 2212, the pre-ruling print at q_G = 8
        "ch10_lq_same_lattice_printed": a.ch10_stretch_lq_printed.value,   # (2600, 3100) (app07:197)
        "ch10_lq_same_lattice_sum": (gauge + fermion + int(a.ch10_stretch_n_anc.lo),
                                     gauge + fermion + int(a.ch10_stretch_n_anc.hi)),   # (2554, 3054): 250-750 ancilla, not 100
        "gauge_fermion_same_lattice_3d": gauge + fermion,                 # 2304, shared with Ch. 10's stretch
        "magnetic_t_per_link_3d": gs["magnetic"],                         # 1.889e4
        "electric_floor_t_per_link_3d": gs["electric"],                   # 3.564e4
        "t_per_link_3d": gs["per_link"],                                  # 5.452e4 -> '5.5e4'
        "t_gauge_per_step_3d": t_gauge_step,                              # 1.047e7
        "t_hop_per_link_3d": hl["t"],                                     # 1.141e5 -> '1.1e5'
        "t_hop_per_step_3d": t_hop_step,                                  # 2.190e7
        "t_mass_per_step_3d": t_mass_step,                                # 1.082e4
        "t_per_step_3d": t_step,                                          # 3.238e7 -> '3.2e7'
        "t_total_3d": t_total,                                            # 3.238e10 -> '3.2e10'
        "t_3d_over_cap": t_total / cap,                                   # 32.4 -> '32x its depth limit'
        "retired_hop_insertion_3d": retired_rest,                         # 2.04e6
        "hop_over_retired_3d": t_hop_step / retired_rest,                 # 10.7
        "applied_stage_t_total_3d": n_evo * applied_step,                 # 9.23e9, before round D
        "ramp_priced_3d": False,
        "eps_syn": float(a.eps_syn),                                      # 1e-2 (R-TOL)
        "eps_rot": eps,                                                   # 6.677e-6 = sqrt(1e-2 / 2.2432e8)
        "t_per_rotation_papers": gs["t_per_rotation_papers"],             # legacy record
        "t_per_rotation_report": gs["t_per_rotation_fermion"],            # hop and mass (r17)
        "synth_rotations_per_shot_3d": n_rot_shot,                        # 2.2e8
        "synthesis_error_3d": synthesis_error(a, n_rot_shot),             # 1e-2 = eps_syn
        "synthesis_error_3d_round_e": synthesis_error(a_round_e, n_rot_shot),   # 224 at 1e-3, record
        "t_total_3d_round_e": t_total_round_e,                            # 7.673e9 at 1e-3, record
        "eps_l_required_3d": eps_l_req,                                   # 3.09e-12
        "legacy_t_per_step_3d": legacy_step,                              # 3e6
        "legacy_t_total_3d": n_evo * legacy_step,                         # 3e9
        "draft_hop_t_per_step_3d": draft_hop_step,                        # 1.06e7
        "draft_hop_t_per_step_3d_3fields": n_stag * draft_hop_step,       # 3 x printed C^S, record
        "hop_toffoli_per_link_3d": hl["toffoli"], "hop_rotations_per_link_3d": hl["n_rot"],   # 3,672, 2933
        "rhop_t_per_link_3d": hl["rhop"]["t"],                            # retired R-HOP at this eps
        "pre_r17_t_total_3d": t_total_pre,                                # 9.537e9
        # R7 (r25): the full q-hat on this lattice is past 2033 (2.4x the register limit); at Ch. 5's 1e3 shots per
        # temperature, 3 temperatures would take this long on one machine of that size
        "past_2033": True,
        "qhat_full_wall_yr_at_1e3_shots_per_T": 3 * 1e3 * (t_total * float(a.t_gate_s) + float(a.shot_overhead_s))
                                                 / float(a.seconds_per_year),   # 3.1
    }
    notes = (
        f"NOT A BOX. app10:90: '1728+576+100 ~ 2400 LQ' (exact {lq}); {n_links} links x ({gs['per_link']:.4g} gauge + "
        f"{hl['t']:.4g} hop) + {V} sites x {ms['t']:.4g} mass = {t_step:.4g} T/step; {t_total:.4g} T over {n_evo:.0f} "
        f"steps. {lq / float(a.rfi_lq_2033):.1f}x the register marker, {t_total / cap:.1f}x the depth limit.",
        "No ramp is priced for the extension in the chapter ('with no ramp priced'); none is added here.",
        f"Synthesis (R-TOL, randomized): {n_rot_shot:.4g} rotations per shot, eps_rot = {eps:.4g}, "
        f"{gs['t_per_rotation']:.3f} T per gauge rotation, {gs['t_per_rotation_fermion']:.3f} per hop/mass rotation. At round E's fixed 1e-3: {t_total_round_e:.4g} T.",
        "Cross-reference to Ch. 10 (app10:95, 'with the 1728+576 gauge and fermion register of Ch. 10's stretch'): app07:197 prints "
        "~2600-3100 LQ for Sigma(72x3) on 4^3 with 3 staggered fields (1728 gauge at 9 qubits/link + 576 fermion + "
        "250-750 BE/bath ancilla = 2554-3054). Same gauge and fermion registers (2304); the ancilla differ (100 here). "
        "Ch. 10 prints NO T-count for that lattice, so the chapter claims a match for the gauge and fermion register only.",
        "Filed under era 'codesign' because common.ERAS has no stretch era; the chapter reserves it for the heavy-ion / q-hat extension.",
    )
    return Result(era="codesign", lq=(lq, lq), hard_ops=(t_total, t_total), breakdown=breakdown,
                  intermediates=inter, epsilon_l=(eps_l_req, eps_l_req), notes=notes)


def model(a: Assumptions, era: str) -> Result:
    """Return the Result for '2028' | '2033' | 'codesign' (= the 3+1D extension of app10:90)."""
    if era == "2028":
        return _model_2028(a)
    if era == "2033":
        return _model_2033(a)
    if era == "codesign":
        return _model_3d_stretch(a)
    raise ValueError(f"era must be one of {ERAS}, got {era!r}")


PUBLISHED = {
    "2028": Published(lq=(230, 230), hard_ops=(1.3e5, 1.3e5),
                      src="app10:2028 box: '230 (150--250 envelope)' (exact: 64 + 144 + 22; r25 R10: two hop links at once "
                          "for F* >= 10, was 216 = 64 + 144 + 8 at r22 E26); "
                          "'Per-shot T 1.3e5 (1 ramp + 2 evolution "
                          "steps; derived here), 1.3x the 1e5 limit' (exact 132,174.9; r18: every rotation at the full fit "
                          "(E20), R-TOL eps = 1.6e-3, the overshoot accepted by ruling 2026-10-01). tol 0.02 covers "
                          "132,175 -> 1.3e5",
                      rel_tol=0.02),
    "2033": Published(lq=(964, 964), hard_ops=(4.9e9, 2.6e10),
                      src="app10:2033 box: '964 (at the 1000-LQ marker)'; 'Per-shot T 4.9e9--2.6e10 (window + 2 x 1e2--1e3 "
                          "ramp steps; electric term a lower bound)' (step rule 2026-10-04: the 4.65 GeV point at 50 a runs "
                          "575 steps, exact 4.9014e9 and 2.6159e10; r26 2.5381e10 at 500 steps; referee J1: the readout "
                          "map runs the ramp backwards on every shot; r25 was 3.9e9--1.5e10 without it; r25 R4: window "
                          "30--50 a). Two-figure prints: tol 0.035",
                      rel_tol=0.035),
    "codesign": Published(lq=(2400, 2400), hard_ops=(3.2e10, 3.2e10),
                          src="NOT A BOX. app10:95 plug-in prose: '1728+576+100 ~ 2400 LQ' (exact 2404), '3.2e10 T over 10^3 "
                              "steps with no ramp priced' (exact 3.238e10; r19 E21, E20 full fit, eps = 3.6e-6). The 3+1D extension "
                              "reserved for the heavy-ion / q-hat stretch. tol 0.02",
                          rel_tol=0.02),
}


def INSTANCE_ROWS(a, era, r):
    """Rows for resources.json: one per priced circuit. 2028 is one row. 2033 is three rows: the box's sampling /
    fragmentation circuit (ONE circuit read out two ways, same register and per-shot T; the shot budgets are extras)
    and the two R7 dipole first-result rows (static Im V on Sigma(216x3) / H_I 3^3 since s648, light-like q-hat on
    2O 3x3x5)."""
    i = r.intermediates
    if era == "2028":
        return [(r"readout pipeline, 2+1D $\mathbb{Z}_3$ $4^2$, 3 stag.", r.lq, r.hard_ops,
                 {"group": "Z3", "lattice": "4^2", "observable": "computational-basis records vs emulation (species map "
                                                                 "checked on emulators; r26, J1)",
                  "shots": i["shots_2028"], "over_reference_t": round(i["t_total_2028_over_cap"], 3),
                  "note": "3 steps (1 ramp + 2 evolution), eps = 1.6e-3 (R-TOL: sqrt(1e-2 / 3792)), every rotation "
                          "at 1.15 log2(1/eps) + 9.2 (E20); Z3 primitives from the dihedral basis with the Z2 dropped, Fourier "
                          "transform 4 T + 2 rotations (derived here, round E); hop 8 strings x HWP(9) per link (derived "
                          "here; Z3 has no draft entry); 1.3x the 1e5 reference, an accepted overshoot"})]
    if era == "2033":
        ds, dl = i["dipole_static"], i["dipole_lightlike"]
        return [(r"species ratios / fragmentation, $\Sigma(72{\times}3)$ $4{\times}8$, 3 stag.", r.lq, r.hard_ops,
                 {"group": "S72x3", "lattice": "4x8 (cut from 4x16)",
                  "observable": "B/M, kaon rates (sampling; V/PS conditional); D_pi^q(z) (FCC era)",
                  "T_priced_at_V": i["V_priced_2033"],
                  "shots_first_result": [round(x) for x in i["shots_first_2033"]],
                  "shots_campaign": [round(x) for x in i["shots_campaign_2033"]],
                  "shots_fragmentation": list(i["shots_ff_printed"]),
                  "wall_time_s_first_result": [round(x) for x in i["wall_first_result_s"]],
                  "wall_time_s_campaign": [round(x) for x in i["wall_campaign_s"]],
                  "wall_time_s_fragmentation": [round(x) for x in i["wall_ff_s"]],
                  "note": "window 30-50 a plus a 1e2-1e3-step ramp and the same ramp backwards for the species "
                          "readout map (R4; r26 J1); electric term at the Sigma(72x3) fast-transform "
                          "lower bound, a floor; hop from the authors' unpublished color-diagonalized gate counts "
                          "(FermionPrimitives_unpub), color squish and parity held once per link (E21), 1.1e5 T per "
                          "link, 72% of the step; walls on one "
                          "machine with the T-depth correction (F* 7.0-18.5); eps at the two ramp ends from R-TOL; every "
                          "rotation at the full fit (E20)"}),
                (r"static dipole $\mathrm{Im}\,V(r,T)$, pure-gauge $\Sigma(216{\times}3)$ $H_I$ $3^3$", (ds["lq"], ds["lq"]),
                 (ds["t_shot"][0], ds["t_shot"][1]),
                 {"group": "S216x3", "hamiltonian": "I", "lattice": "3^3",
                  "observable": "Im V(r = a, T) from singlet survival (kappa as byproduct)",
                  "shots_per_T": [round(x) for x in ds["shots_per_T"]],
                  "wall_time_s_first_result": [round(x) for x in ds["wall_first_s"]],
                  "note": "R7 first-result row; s648: Sigma(216x3) (11 qubits/link, Gustafson_S648_inprep) with H_I as "
                          "Ch. 5's 2033 box, all primitives compiled in the draft; q1 (ruling B): medium quench-prepared "
                          "from the electric vacuum (20-step coupling ramp + c/T thermalization, c = 2-2 pi, 90-282 steps) "
                          "before the sources enter, ETH / typicality assumed; 5.4-15x the 1e9-T limit, stated; step rule "
                          "on the improved free dispersion (85 / 240 dipole steps); t = 0.5 and 1 fm/c, each shot at its "
                          "own time; shots binomial on P_s, no estimator penalty; measures P_s, not |W|^2 (J4)"}),
                (r"light-like dipole $\hat q$, pure-gauge $2O$ $3{\times}3{\times}5$", (dl["lq"], dl["lq"]),
                 (dl["t_shot"][0], dl["t_shot"][1]),
                 {"group": "2O", "lattice": "3x3x5", "observable": "q-hat from singlet survival at L = 0.5 and 1 fm",
                  "shots_per_T": [round(x) for x in dl["shots_per_T"]],
                  "wall_time_s_first_result": [round(x) for x in dl["wall_first_s"]],
                  "note": "R7 first-result row; compiled 2O transform; q1 (ruling B): quench-prepared medium (20-step "
                          "ramp + c/T, c = 2-2 pi) before the sources enter, ETH / typicality assumed; fits the 1e9-T "
                          "limit; shots binomial on P_s, no estimator penalty; O(1/N_c^2) color correction (J4)"})]
    if era == "codesign":
        return [(r"3+1D stretch, $\Sigma(72{\times}3)$ $4^3$, 3 stag.", r.lq, r.hard_ops,
                 {"group": "S72x3", "lattice": "4^3",
                  "note": "plug-in prose, not a box; the full-QCD q-hat vehicle, past 2033 (2.4x the register limit; R7); no ramp priced; "
                          "hop from unpublished gate counts; eps from R-TOL; every rotation at the full fit (E20), gauge tables "
                          "from the subgroup papers"})]
    return []
