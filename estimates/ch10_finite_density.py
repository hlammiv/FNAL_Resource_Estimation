"""Ch. 10 — Finite-density EoS of QCD. Reproduces the LQ and hard-op numbers of
applications/app07_finite_density.tex.

INPUTS AND SOURCES
  d, L, N_stag (2028: 2, 2, 1; 2033: 3, 2, 3)   Stated, app07:138-139 / 162-163
  N_c = 3                                        Stated, app07:91
  q_G = 8 (Sigma(36x3)), 9 (Sigma(72x3))         GROUPS[group].link_qubits (compiled, R1)
  primitive T-costs, multiplicities              GROUPS["S36x3"].primitives, groups.PRIMCOST["KS"]
                                                 (arXiv:2405.05973 tab:tgatecost/primcost,
                                                 arXiv:2408.00075 FFT), full fit (E20)
  eps_syn = 1e-2 total synthesis error per shot  R-TOL ruling (common.EPS_SYN); eps_rot per circuit
                                                 = eps_rot_for(N_rot) (was Stated 1e-4, app07:96)
  N_anc: ~100 (2028), 250-750 (2033), 250-750 (stretch)   Stated, app07:93,146 / 93,168 / 93,197
  T cap 1e5 (2028)                               Cited DOE_RFI_2026
  N_Trotter = 1 (2028)                           Stated app07:103,148 (= floor(cap / c_BE))
  c_BE = 5e5 (2033) T/step                       Stated, app07:96,105,169 BEFORE R-HOP; superseded
                                                 by the derived 3.91e5 (kept as a record only)
  hop: Fermion_Primitives draft counts, undo,    r17/r19 (groups.hop_link_cost; FermionPrimitives_unpub),
       link sharing, MBU, HWP phasing, full fit   estimated pieces labelled there
  mass: HWP(N_c N_stag) per site, mu_B folded in  r17 (groups.staggered_mass_site; FP section_mass.tex:3-10)
  Gibbs at V=2^2: >~14 sampler steps             Stated, app07:103 (restated from '>~1e6 T' = 1e6/7.25e4); record
                                                 since referee F2 (the V=2^2 sample is now priced, 1.3e12 T)
  Gibbs sampler unit cost (F2)                   ESTIMATED HERE from arXiv:2311.09207 (OFT eq. OpOFT, coherent
                                                 term eqs. mainBDef/b1/b2, Thm. I.2): gibbs_sample(); jump set,
                                                 sweeps (10/20/40), channel error 1e-2, control by 1 Toffoli per
                                                 rotation, beta||H|| at V=2^2 Assumed
  beta||H|| ~ 1e2, c ~ 20                        Stated, app07:105 (c called least controlled)
  delta = 0.1 absolute on chi4/chi2, 5x5 grid    Assumed (referee F1); grid points evenly spaced
  5% on p, 10% on chi4/chi2                      Stated, app07:56,167,224 (not connected to shots)
  t_gate = 1 us, 5-yr horizon                    Stated, app07:151,167,201 (r23: one machine, serial; no machine count)
  shot overhead ~0.1 ms                          report rule (as Ch. 8 shot_overhead_s)
  QSVT/TPQ (2033 headline, author ruling       delta_beta 1e-4, success prob 1e-2 (amplitude 0.1), 2k+1 calls: Assumed,
    2026-10-05): Ch. 9's construction          as Ch. 9 (ch09_chiral_gauge.py); query 2-4 c_step Stated app07:107;
                                               product reference 1 rotation/qubit Assumed; reflections derived here
  fault budget 0.1 expected faults per shot      R3 ruling (TRACKED_CHANGES.md)
  eps_l floor 1e-8 (RFI)                         Cited DOE_RFI_2026
  T-depth per shot (2028, 2033)                  Assumed, factory.json r25 (factories/ch10_depth_v2.py)
  2028 shots 2e3 per sector at sigma^2 <= 1      Stated app07:145; 2.2% rms (author approved 2026-10-02, R1)
  shots (grid and first result)                  DERIVED (referee F1, 2026-10-04): delta-method variance of the
                                                 k-statistic ratio k4/k2 on toy N_B distributions (baryon or quark
                                                 Skellam, kappa_2(0) Assumed), reweighting by self-normalized
                                                 importance weights; ensemble pairs Assumed from scratch searches
WHAT IS NOT DERIVED HERE
  The staggered hop comes from an UNPUBLISHED draft's gate tables; the frame undo, the SU(3)
    squish uncompute, the per-field colour rotations and the HWP phasing are estimated here
    (groups.py, NEEDS_AUTHOR 'Cross-chapter: fermion hop'). SCALING where estimated.
  Gauge primitives: the papers' tables with every rotation at the full fit (E20); the papers' own
    slope-only price is the record c_be_gauge_t_per_step_papers_convention.
  The Gibbs bound at V=2^2 is carried in sampler steps (>~14), not derived (item 10 open).
  N_anc at every rung (Stated). 2028: ~100 is a sized budget (ruling ch10-ancilla-contents A):
    workspace 29 (serial reuse; since r19 the draft's 24-qubit colour-squish scratch is held through the hop,
    so it sits beside the hop's largest piece, the C^7X ladder 5; the hop's HWP(6) holds 4 = 6 - w(6) (E26, r22;
    was 7) and with the 1 repeat-until-success synthesis ancilla is also 5; the other pieces U_FFT 8/U_Tr 7/
    U_inv 4/U_mul 2/mass HWP(3) 1 reuse it) to 55 (no reuse) plus 1 Hadamard-test qubit is accounted, the rest
    is margin; algorithmic register 106-132 of the ~180 sized (was 31-60 and 108-137 with the weight register
    added on top of the adder carries). The synthesis ancilla is not a line of the no-reuse sum; it is in the
    margin. The draft's separate flag and parity registers are not sized.
    Stretch: carries the 2033 250-750 band (ruling ch10-stretch-ancilla a).
  Shots at 2028 (~2e3, Stated): 2.2% rms on each sector's E_0 at a per-shot variance bound sigma^2 <= 1
    (author approved 2026-10-02); the term variances of the <H> family-basis estimator are not sized.
  Shots come from toy N_B distributions (Skellam in baryons or quarks), not from this Hamiltonian: kappa2 of
    N_B at V=2^3 and its T dependence are open (NEEDS_AUTHOR). The 5% pressure target is not
    separately costed (ruling ch10-pressure-5pct a): the chi4 budget at eps=0.1 binds, p is
    inherited from the same n_B shots, its mu_B quadrature error is a systematic.
  2033 hard ops (author ruling 2026-10-05) are the QSVT/TPQ thermal state, CONDITIONAL on a reference state whose
    filtered amplitude is ~0.1 and whose eigenstate weights are thermal (referee C1 of Ch. 9: a typical reference has
    success probability Z e^{beta E_0}/D, exponentially small; a product reference is biased). The query cost 2-4 c_step
    and the use of beta||H|| (not the BE normalization lambda) in the degree are Stated, not derived. The Gibbs-sampler
    estimate below is kept as a record.
  2028 hard ops are a LOWER BOUND on the executed circuit (referee G2): one ramp step (ramp length and energy
    readout open). 2033 hard ops were an ESTIMATE from referee F2 (2026-10-04) to 2026-10-05: the Chen-Kastoryano-Gilyen sampler
    priced per step (jump, OFT and inverse, coherent term and inverse, control) on the chapter's Trotter step at
    step 1/||H||; low/central/high 7.6e12 / 8.39e13 / 9.5e14. Printed range (open item 1(a), 2026-10-04): central
    8.39e13 to x K at 20 sweeps 4.42e14 (t_per_shot_range). The mixing time (c = 20 sweeps of 288 jumps) is
    still assumed; the old '>~5.5e9' (2e3 uncontrolled steps) is t_hsim_lower_bound.
WALL TIME AND FACTORIES
  t_gate_s = 1e-6 s per T, shot_overhead_s = 1e-4 s per shot, ONE machine, serial (r23, R9).
  T-depth per shot from factory.json (Ch. 10 runs, scratchpad factories/ch10_depth_v2.py), low/high:
    2028 one Trotter step 2.67e4 (T-par colour frame, 8 link rounds: its ~60 parity ancilla per hop fit the
    ~100 budget only at 8 rounds; 4 rounds, 1.43e4 / F* 28, is over budget) - 4.05e4 (as compiled, 2d = 4
    colour classes); F* = 9.9-15. F* < 10 at the as-compiled end, so the shot is depth-limited there: 0.405 s, not
    0.400 s (corrected wall, R9 rule 12). Fully serial 1.84e5 (F* 2.2) is the pessimistic bound, not priced.
    2033 Gibbs shot 9.85e7 (T-par, 12 link rounds: ~51 + ~170 parity ancilla per hop fit the 750 band only at
    ~3 concurrent hops; 6 rounds, 5.23e7 / F* 106, would need ~880) - 3.16e8 (as compiled, 12 rounds at the
    250-ancilla end); F* = 17-56, so the 10-factory baseline holds. F2: the priced sampler repeats the same step,
    controlled, so the band is scaled by the T ratio (x15,198, F* carried, Assumed); floor of the grid 1.2-3.9e5 yr.
  2028 (one tier): E_0(N_B) in 2-3 sectors at 2e3 shots each (2.2% rms at sigma^2 <= 1, author approved
    2026-10-02); wall_first_result_s = None; wall_campaign_s = the 2-3 sectors.
  2033 since the author ruling of 2026-10-05: QSVT shot 2.77e9-5.63e9 T, 2.8e3-5.6e3 s; first result 2,926-15,861 samples
    -> 0.26-2.8 yr (paired low/low, high/high); campaign 2.58e5 shots -> 23-46 yr (4.5-9.2x the horizon); kappa_2 0.05
    80-162 yr, 0.1 197-400 yr. T-depth: factory.json band scaled by each end's own T ratio (F* 17.5-56 at either end,
    f_star_matched); the cross-paired export is 8.6-114 and its (T_lo, D_hi) corner is not a circuit, so the walls use
    the matched ends (serial at both). Depth limit of the grid 4.0-26 yr.
  Records below (F2, Gibbs sampler):
  2033 first result (R1): one T row (T = 150 MeV, mu_B = 0-600 MeV), n_B (mu_B >= 150) and chi4/chi2 at 30%,
    Delta p from the same samples, each toy at min(direct, two reweighted ensembles): 2,926-15,861 samples
    -> wall_first_result_s 7.8e3-4.2e4 yr at 2.7 yr per shot (F2; 0.51-2.8 yr at the old bound).
    Campaign: the 5x5 grid at |delta(chi4/chi2)| <= 0.1, 257,603 shots at kappa_2(0) = 0.01 -> 6.8e5 yr (F2; 45.1 yr
    at the old bound). Records:
    9.08e5 (0.05), 2.24e6 (0.1), quark-like 3.88e5 (10% relative) / 4.8e3 (0.1 absolute); reweighting saves 11%
    at 0.01 and loses at 0.05.
    Since F2 the 2033 walls are estimates at the central sampler (band 7.6e12-9.5e14 carried as
    wall_campaign_band_yr / wall_first_result_band_yr); the 5-year horizon holds 1.9 shots, the grid needs a
    1.4e5x cut and the first result 1.6e3-8.4e3x.

Author ruling (Henry, 2026-10-05, "Make the QSVT/TPQ route in Ch 10 as well"): the 2033 headline is the QSVT/TPQ thermal
state of Ch. 9 (qsvt_thermal). Checked against the old rail: the degree used tolerance 1e-3 (now Ch. 9's 1e-4: 26.3 ->
30), the call count was R x 2 = 20 filter-equivalents with R = overlap^(-1/2) (now 2k + 1 = 17, k = ceil(pi/(4 arcsin
0.1)) = 8, Ch. 9's rule), and the reference preparation, signal phases and reflections were not priced (now 1e-4 of the
shot); readout of N_B is computational-basis. 2.86e9-5.81e9 -> 2.77e9-5.63e9. The 2028 benchmark stays the ground-state
ramp: the same filter on V=2^2 is 4.7e8-9.6e8 T (qsvt_2x2_t). Gibbs sampler kept as a record (gibbs_*).

Referee round (2026-10-04, editorial_review/REFEREE_REPORT.md F1-F4, G2): Eq. Nshot_fd's kappa_8/kappa_4^2 ~ 1e2
is replaced by the delta-method variance of the k-statistic ratio (F1). The shot-audit records (first result
798-2,873; reweighted grid 3.9e4-1.3e5) could not be reproduced and are replaced by the derivation here; the
direct-grid audit values (2.58e5, 2.24e6, 3.88e5) are reproduced. Campaign 2.5e5 -> 2.58e5 shots, 43.7 -> 45.1 yr,
8.7x -> 9.0x; first result 51 d-0.50 yr -> 0.51-2.8 yr. No per-shot T or LQ moves.

HISTORY AND COST MODEL
The chapter's cost model is Eqs. (Nq_fd)-(Nshot_fd), app07:87-89:
    N_q     = d L^d q_G + N_stag N_c L^d + N_anc        (q_G = compiled per-link register)
    N_gate  = t_mix * c_BE * polylog(1/eps),   t_mix = c * beta ||H||
    N_shot  = eps^-2 (kappa_8/kappa_4^2) N_T N_mu
2028: one Trotter step is priced; c_BE(2028) is DERIVED from the compiled Sigma(36x3)
primitives (groups.py: magnetic_per_link + electric_per_link with the arXiv:2408.00075 FFT,
rotations at the full fit since E20) for the gauge terms on 8 links, plus the 8 staggered hops from the
Fermion_Primitives draft (groups.hop_link_cost, r17/r19) and the staggered mass on 4 sites. Since
r17 the step is over the 1e5 cap (4.0x since r19), so the cap buys no step.
2033: [c(20) x beta||H||(1e2)] x c_BE, with c_BE = gauge (24 links at d=3) + draft hop
(24 links, N_stag = 3) + staggered mass (8 sites), derived here; it replaces the Stated 5e5,
which is kept only as a record (superseded).

r17 (ruling (e) H. Lamm 2026-10-01; apply_log/r17_core.md, apply_log/r17_ch10.md): the hop is
groups.hop_link_cost("S36x3", N_stag): the unpublished Fermion_Primitives gate tables
(FermionPrimitives_unpub) under the report rules, with the pieces the draft leaves out estimated
(frame undo, SU(3) squish uncompute, per-field colour rotations, MBU ladders, HWP phasing). The
staggered mass is groups.staggered_mass_site: one HWP group of N_c N_stag per site, mu_B N_B folded
into the same angle (it commutes with H). Hop and mass rotations are at the full fit
1.15 log2(1/eps) + 9.2. Since E20 (r18, 2026-10-01; apply_log/r18_core.md) the gauge primitives are too
(n_rot = the paper's log coefficient / 1.15, Toffoli constants kept). R-HOP is kept as a legacy comparison
intermediate only.

r19 (E21 rulings, H. Lamm 2026-10-01; apply_log/r19_core.md, apply_log/r19_ch10.md): (1) the hop's colour
squish and parity flags are computed once per link and held through V(x), V(y), hop, V^dag(x), V^dag(y)
(share="link"), with the frame undo kept (undo=True; the draft draws V_C^dag but does not count it,
PD/section_hamiltonian.tex:86-91 vs PD/section_resources.tex:15). The per-application recompute
(share="draft") is the conservative sensitivity. (3) SU(3) colour rotations = 2 N_angles (su3_diag.tex:120
"6 N_angles" is a draft slip). E20 is applied to the chapter in the same round. Since r19 each shot also
carries ~0.1 ms of overhead (report rule; it moves no printed number).
2028: N_rot 10,908 (3,708 gauge + 7,192 hop + 8 mass), eps 9.575e-4, step 399,956.8 T (4.0x the cap).
2033: N_rot 1.517e8 per shot, eps 8.118e-6, c_BE 2.760e6, 5.521e9 T/shot (5.5x the 1e9 reference).

r23 (ruling H. Lamm 2026-10-02, "we don't want to anywhere assume we have multiple machines";
apply_log/r23_ch10.md): one machine, every wall time serial at 1 us per T + ~0.1 ms per shot. n_machines (4),
wall_parallel_yr (10.9) and machines_for_horizon (9) are retired. In their place: the 5-year horizon holds
28,580 shots = 2.86 grid points (two full points of 25), and the full 5x5 grid (43.7 yr) needs an 8.75x cut
in shots x per-shot cost (6.31e8 T/shot at 2.5e5 shots). No per-shot T or LQ number moves.

R-TOL report-wide (ruling H. Lamm 2026-09-29; apply_log/rtol_ch10.md): total synthesis error
per shot eps_syn = 1e-2 (common.EPS_SYN); each circuit sets its own rotation tolerance
eps_rot = common.eps_rot_for(N_rot) = sqrt(eps_syn / N_rot), N_rot = synthesized rotations in
ONE shot of that circuit (_rot_per_step x steps). The papers' slope-only 1.15 log2(1/eps) is
kept; only eps changes (was a fixed 1e-4). 2028: N_rot = 3,740, eps_rot = 1.635e-3, 10.645 T
per rotation, 72,515.5 T/shot (was 89,854.5). 2033: N_rot = 2.28e7, eps_rot = 2.094e-5,
17.875 T per rotation, 7.818e8 T/shot (was 7.227e8). The QSVT rail and the V=3^3 promotion
are their own circuits and set their own eps from their own N_rot.

R-HOP report-wide (ruling H. Lamm 2026-09-29; apply_log/rhop_ch10.md; RETIRED by r17, kept as the
legacy comparison): the hop was priced per link per Trotter step by the shared groups.rhop_link: 2 N_stag link-controlled
colour moves at one U_mul (308 T, arXiv:2405.05973 tab:tgatecost; ruling VERTEX) plus
2 Pauli strings x HWP(N_c N_stag) (common.hwp_group, papers' convention, eps_rot):
(691.12 T/link at N_stag = 1, 2,068.25 T/link at N_stag = 3 at the old fixed 1e-4). The cap-filling UNSOURCED
allowances (1.57e4 at 2028, 1.88e5 at 2033) are gone; totals follow.

Applied 2026-09-28 (rulings R1-R4, apply_log/ch10.md): register widths are the
compiled 8 (Sigma(36x3)) / 9 (Sigma(72x3)) through GROUPS[..].link_qubits; the 2028
step is repriced from the papers (ch10-cbe-2028 option A) so the cap buys ONE step,
not ten; the QSVT rail lower end is 5e8 (ch10-qsvt-lower-end); epsilon_l is 0.1
expected faults per shot (R3); carry exact, round once (R4).

"""

from __future__ import annotations

import math
from dataclasses import dataclass

from estimates.common import (Assumed, Cited, Stated, Uncited, Tagged, Primitive,
                              CircuitStatus, Result, Published, ERAS,
                              t_per_rotation, eps_per_rotation, toffoli_t,
                              EPS_SYN, EPS_SYN_SRC, eps_rot_for, depth_exports,
                              REACTION_TIME_S, FACTORY_BASELINE)
from estimates.groups import (GROUPS, PRIMCOST, MAGNETIC, magnetic_per_link,
                              electric_per_link, hop_link_cost, staggered_mass_site,
                              fp_printed, FP_KEY, mcx_toffolis)
from estimates.common import HWP_SRC, RUS_ANCILLA, hwp_ancilla

TEX = "app07"
SECONDS_PER_YEAR = 365.25 * 86400.0
SECONDS_PER_MONTH = SECONDS_PER_YEAR / 12.0

# The electric term per link per step of H_KS: nF = PRIMCOST["KS"]["U_F"] Fourier
# transforms (the arXiv:2408.00075 FFT) plus one diagonal phase U_phi (arXiv:2405.05973
# :588). groups.electric_per_link prices exactly this; the tuple names the pieces for
# the breakdown.
# Ancilla column of the Fermion_Primitives draft's squish table for Sigma(3x36) (unpublished; su3_diag.tex:103).
SQUISH_ANCILLA_S36 = 24
ELECTRIC_PIECES = (("U_FFT", lambda d: PRIMCOST["KS"]["U_F"](d)), ("U_phi", lambda d: 1))


def _rng(t: Tagged) -> tuple[float, float]:
    return (t.lo, t.hi)


@dataclass(frozen=True)
class Assumptions:
    # ---- lattice geometry and encoding --------------------------------------
    nc: Tagged = Stated(3, f"{TEX}:91", "N_c=3 is the color count")
    group: Tagged = Stated("S36x3", f"{TEX}:138,162", "Sigma(36x3), 8 qubits/link (compiled encoding, R1); both boxes")
    ham: Tagged = Cited("KS", "arxiv_2405_05973 tab:primcost",
                        "Kogut-Susskind: plaquette + E^2 are the chapter's families (a), (b) at app07:96")
    # 2028 benchmark: 2+1D, V=2^2, one staggered field
    d_2028: Tagged = Stated(2, f"{TEX}:134,138", "2+1D SU(3) ground state")
    L_2028: Tagged = Stated(2, f"{TEX}:139", "V=2^2")
    n_stag_2028: Tagged = Stated(1, f"{TEX}:139", "1 staggered field (2 degenerate light tastes)")
    n_anc_2028: Tagged = Stated(100, f"{TEX}:93,146",
                                "~100 ancilla as a sized budget (ruling ch10-ancilla-contents A): primitive "
                                "workspace 29-55 (groups.py S36x3 + the r17 hop, its 24-qubit squish scratch held per link (r19), and the mass; "
                                "phasing groups at k - w(k) ancilla, E26) + 1 Hadamard-test qubit "
                                "accounted, rest margin")
    # 2033 target: 3D, V=2^3, three staggered fields
    d_2033: Tagged = Stated(3, f"{TEX}:158,162", "3D SU(3)")
    L_2033: Tagged = Stated(2, f"{TEX}:162", "V=2^3")
    n_stag_2033: Tagged = Stated(3, f"{TEX}:163", "3 staggered fields (N_f=4+2 tastes)")
    n_anc_2033: Tagged = Stated((250, 750), f"{TEX}:93,168",
                                "250-750 for the 2033 route: BE workspace and the parity ancilla of parallel hops (was the Gibbs route's BE/bath register; relabelled 2026-10-05)")
    # post-2033 stretch (codesign era): Sigma(72x3) on V=4^3
    group_stretch: Tagged = Stated("S72x3", f"{TEX}:197", "Sigma(72x3), 9 qubits/link (compiled encoding, R1)")
    d_stretch: Tagged = Stated(3, f"{TEX}:197", "production 3D")
    L_stretch: Tagged = Stated(4, f"{TEX}:197", "V=4^3")
    n_stag_stretch: Tagged = Stated(3, f"{TEX}:197", "3 staggered fields")
    n_anc_stretch: Tagged = Stated((250, 750), f"{TEX}:93,197",
                                   "the 250-750 BE/bath ancilla of the 2033 route, carried V-independently "
                                   "(ruling ch10-stretch-ancilla a)")

    # ---- per-step cost -------------------------------------------------------
    # c_BE(2028) is NOT a field: it is derived from GROUPS[group] through
    # groups.magnetic_per_link / electric_per_link at eps_rot (contract rule 6).
    eps_syn: Tagged = Assumed(EPS_SYN, "R-TOL: total synthesis error per shot, report-wide; each circuit sets "
                                       "eps_rot = sqrt(eps_syn / N_rot) from its own rotation count (replaces the "
                                       "Stated fixed 1e-4 of app07:96)", EPS_SYN_SRC)
    rot_synthesis_paper: Tagged = Cited("rus", "arxiv_2405_05973,arxiv_2408_00075; apply_log/r18_core.md",
                                        "E20 ruling (2026-10-01, supersedes R2's papers' convention): every gauge rotation at the "
                                        "full fit '1.15 log2(1/eps) + 9.2', n_rot = the paper's log coefficient / 1.15, the "
                                        "papers' Toffoli constants kept; eps from the R-TOL incoherent (randomized) split. "
                                        "groups.PrimitiveCost.t() prices in it; t_papers() is the slope-only legacy record")
    toffoli_convention: Tagged = Cited("textbook", "arxiv_2405_05973",
                                       "7 T per Toffoli, as in the papers' tables (PAPER_COSTS.md Sec. 4)")
    t_cap_2028: Tagged = Cited(1e5, "DOE_RFI_2026", "2028 per-shot hard-op cap; app07:134,147")
    t_ref_2033: Tagged = Cited(1e9, "DOE_RFI_2026", "2033 per-shot hard-op reference; app07:105,158,169")
    n_trotter_2028: Tagged = Stated(1, f"{TEX}:103,143",
                                    "one second-order Trotter step is priced; it is 4.0x the 1e5 cap (r19), "
                                    "so floor(cap / c_BE) = 0 and no step fits")
    n_trotter_2028_original: Tagged = Stated(10, f"{TEX}:103,148",
                                             "'the ten-step adiabatic ramp originally scoped ... does not' fit; "
                                             "re-scoping the validation target is the authors' call")
    c_be_2033_superseded: Tagged = Stated(5e5, f"{TEX}:96,105,169 (before R-HOP)",
                                          "~5e5 T/step as printed before R-HOP (gauge 3.1e5 + an unsourced hop "
                                          "allowance); superseded by the derived c_BE, kept as a record only")
    # ---- staggered hop and mass from the Fermion_Primitives draft (r17, ruling (e) 2026-10-01) ----
    hop_synthesis: Tagged = Assumed("rus", "every hop and mass rotation at the full fit 1.15 log2(1/eps) + 9.2 "
                                           "(r17 rule 7; the draft prices 9.2 + 1.15 L too)",
                                    "apply_log/r17_core.md; FermionPrimitives_unpub resources.tex:7-10")
    hop_undo: Tagged = Assumed(True, "V_g^dag on both sites after the diagonal hop (ESTIMATED; the draft draws it, "
                                     "PD/section_hamiltonian.tex:86-91, but its 2 C^G counts only the 2 forward "
                                     "applications, PD/section_su2_diag.tex:238, section_resources.tex:15)",
                               "groups.fermion_hop_counts; E21 (1), apply_log/r19_core.md")
    hop_share: Tagged = Assumed("link", "colour squish + parity computed once per link and held through V(x), V(y), "
                                        "hop, V^dag(x), V^dag(y) (E21 (1) 'make the switch'); 'draft' (recompute per "
                                        "frame application) is the conservative sensitivity, reported beside it",
                                "groups.fermion_hop_counts; apply_log/r19_core.md")
    hop_mcx: Tagged = Assumed("mbu", "C^nX ladders at n-1 Toffolis (measurement-based uncompute, report rule)",
                              "groups.mcx_toffolis")
    hop_phasing: Tagged = Assumed("hwp", "the 2 N_c N_stag equal-angle Rz of the diagonal hop as one HWP group "
                                         "(derived here)", "groups.fermion_hop_counts; arxiv_1709_06648")
    gibbs_2x2_min_steps: Tagged = Stated(14, f"{TEX}:103",
                                         "'Gibbs sampling on V=2^2 needs >~14 sampler steps': the earlier '>~1e6 T per "
                                         "shot' read at the then 7.25e4 T/step (1e6 / 72,515.5 = 13.8); lower bound, "
                                         "no beta||H|| given (item 10)")
    gibbs_2x2_t_per_shot_superseded: Tagged = Stated(1e6, f"{TEX}:103 (before r17)",
                                                     "'>~1e6 T-gates per shot' as printed before r17; record only")

    # ---- step counts ---------------------------------------------------------
    beta_h_2033: Tagged = Stated(1e2, f"{TEX}:105", "beta||H|| ~ 1e2 for the KMS Gibbs preparation")
    c_mix: Tagged = Stated(20, f"{TEX}:105",
                           "t_mix = c beta||H|| with c ~ 20; 'least controlled constant in the 2033 costing'. Referee F2 "
                           "(2026-10-04): read as c sweeps of the jump set (t_mix ~ c at per-jump rate 1)")
    # ---- referee F2 (2026-10-04): one Chen-Kastoryano-Gilyen sampler step, ESTIMATED HERE --------------------
    # A step is one unit of Lindbladian time in the normalization of arxiv_2311_09207 (||sum_a A^a† A^a|| <= 1):
    # one local jump a, drawn at random, with its operator Fourier transform (OFT) and coherent term B.
    gibbs_sweeps: Tagged = Assumed((10, 40), "sweeps of the jump set per Gibbs sample, low/high; central = c_mix "
                                             "(20, Stated). One sweep = |A| sampler steps", "referee F2; app07:105")
    gibbs_eps_channel: Tagged = Assumed(1e-2, "channel-approximation error per Gibbs sample, split evenly over its "
                                              "sampler steps; sets the OFT and B time truncations",
                                        "referee F2; same budget as common.EPS_SYN")
    gibbs_beta_h_2028: Tagged = Assumed(1e2, "beta||H|| for a Gibbs sample at V=2^2: the 2033 value, none is stated "
                                             "for the 2028 instance", "referee F2; app07:105")
    gibbs_ctrl_toffoli_per_rot: Tagged = Assumed(1, "controlled Trotter step: each synthesized rotation is controlled "
                                                    "by one Toffoli onto an AND ancilla (7 T, MBU uncompute), the control "
                                                    "qubit's own phase aggregated into one extra rotation per step; "
                                                    "Toffoli/Clifford compute-uncompute pairs need no control",
                                                 "referee F2, derived here")
    gibbs_structure: Tagged = Assumed(("one OFT, no inverse, no coherent term",
                                       "OFT and its inverse + coherent term B and its inverse",
                                       "central x K, the per-unit-time segment order of black-box Lindbladian "
                                       "simulation, K = ln(t/eps)/lnln(t/eps)"),
                                      "controlled evolution per sampler step, low/central/high (Gaussian weight, "
                                      "sigma_E = 1/beta; OFT window e^{-t^2/beta^2}, |t| <= beta sqrt(ln 1/eps_u); B: "
                                      "|t| <= ln(1/eps_u)/(2 pi), |t'| <= sqrt(ln(1/eps_u)/4)). Time register in "
                                      "sign-magnitude form: the palindromic second-order step has S(-dt) = S(dt) with "
                                      "negated angles, so the sign is CNOT conjugation (Clifford) and a window [-T, T] "
                                      "costs controlled evolution T per side (offset binary: 2T, central +40%). Low drops "
                                      "B, so it is a floor without exact KMS detailed balance. Central = one query per "
                                      "unit time: not given by Thm. I.2 (which gives x K); supported up to an unpriced "
                                      "O~(1) constant by the discrete-time channel of arXiv:2405.20322",
                                      "arxiv_2311_09207 eqs. OpOFT, mainBDef, b1, b2; Thm. I.2")

    # ---- shots, accuracy targets and wall time --------------------------------
    shots_2028: Tagged = Stated(2e3, f"{TEX}:150", "~2e3 shots, single (T,mu_B) point; no eps given")
    eps_stat: Tagged = Assumed(0.1, "absolute target |delta(chi4/chi2)| <= 0.1 per grid point, 10% at the hadron-gas "
                                    "value 1 and well defined near zeros of chi4 (referee F1)", f"{TEX}:Nshot_fd")
    grid_T_mev: Tagged = Assumed((100, 125, 150, 175, 200), "5 evenly spaced temperatures in [100,200] MeV",
                                 f"{TEX}:objectives, fig:phase_diagram")
    grid_mu_mev: Tagged = Assumed((0, 150, 300, 450, 600), "5 evenly spaced mu_B in [0,600] MeV",
                                  f"{TEX}:objectives, fig:phase_diagram")
    nb_toy_k2_priced: Tagged = Assumed(0.01, "kappa_2(N_B) at mu_B = 0 of the priced toy: Skellam in baryons, "
                                             "the same at every T (the true kappa_2 at V = 2^3 is open)",
                                       "referee F1; scratchpad referee/ch10_f1.py")
    nb_toy_k2_band: Tagged = Assumed((0.01, 0.1), "baryon-Skellam kappa_2(0) band carried beside the priced value",
                                     "referee F1")
    nb_toy_k2_mid: Tagged = Assumed(0.05, "middle of the band", "referee F1")
    nb_toy_quark_k2: Tagged = Assumed(0.05, "quark-Skellam toy, kappa_2(N_B) at mu_B = 0 (chi4/chi2 = 1/9)",
                                      "referee F1")
    fr_rw_pairs: Tagged = Assumed((("baryon", 0.01, (125, 600)), ("baryon", 0.05, (175, 600)),
                                   ("baryon", 0.1, (175, 600)), ("quark", 0.05, (200, 600))),
                                  "best two-ensemble pair (mu_B in MeV) for the T = 150 MeV first-result row, found on a "
                                  "25-MeV grid (scratchpad referee/run2.py); each toy takes min(direct, reweighted)",
                                  "referee F1")
    grid_rw_pairs_k2lo: Tagged = Assumed(((100, (150, 600)), (125, (125, 600)), (150, (125, 600)), (175, (125, 600)),
                                          (200, (150, 600))),
                                         "best two-ensemble pair per T for the grid at kappa_2(0) = 0.01 (record)",
                                         "referee F1; scratchpad referee/run2.py")
    grid_rw_pairs_k2mid: Tagged = Assumed(((100, (200, 600)), (125, (175, 600)), (150, (175, 600)), (175, (175, 600)),
                                           (200, (175, 600))),
                                          "best two-ensemble pair per T for the grid at kappa_2(0) = 0.05 (record)",
                                          "referee F1; scratchpad referee/run2.py")
    n_T: Tagged = Stated(5, f"{TEX}:56,165", "5 temperatures, T in [100,200] MeV")
    n_mu: Tagged = Stated(5, f"{TEX}:56,165", "5 chemical potentials, mu_B in [0,600] MeV")
    eps_pressure_target: Tagged = Stated(0.05, f"{TEX}:56,167,224",
                                         "'5% on p', not separately costed: chi4 at eps=0.1 sets the shots, p "
                                         "inherited (ruling ch10-pressure-5pct a)")
    eps_chi4_chi2_target: Tagged = Stated(0.10, f"{TEX}:167", "'10% on chi4/chi2'")
    t_gate_s: Tagged = Stated(1e-6, f"{TEX}:151,175", "1 us per T-gate")
    shot_overhead_s: Tagged = Assumed(1e-4, "per-shot overhead ~0.1 ms (register init, final readout, decode), "
                                            "report rule; as Ch. 8 shot_overhead_s", "R17 g_rate; ch08_baryogenesis.py")
    # r23 (H. Lamm 2026-10-02, apply_log/r23_ch10.md): the report assumes one machine. The Stated
    # n_machines = 4 ("four shot-parallel machines") is retired, with the derived 10.9 yr on four and the
    # 9 machines for the horizon. Every wall time is serial: shots x (T x 1 us + ~0.1 ms).
    campaign_horizon_yr: Tagged = Stated(5, f"{TEX}:201", "'Campaign horizon 5 years'; the 5x5 grid is 45.1 yr "
                                         "serial, so 5 years hold 2.86e4 shots and the full grid needs a 9.0x "
                                         "algorithmic cost reduction")

    # ---- r25 (rulings R1, R9, H. Lamm 2026-10-02): depth, first result, kappa2 dependence -------------
    sigma2_2028: Tagged = Assumed(1.0, "per-shot variance bound sigma^2 <= 1 of the normalized E_0 estimator; "
                                       "2e3 shots -> 2.2% rms per sector (author approved 2026-10-02)",
                                  "r25 task rulings; app07:145")
    n_sectors_2028: Tagged = Stated((2, 3), f"{TEX}:141", "'ground-state energies E_0(N_B) in 2--3 sectors'")
    t_depth_2028: Tagged = Assumed((26719.289742577264, 40488.90205167576),
                                   "T-depth of the one-step 2028 shot: low = T-par colour frame (commuting Euler-factor "
                                   "rotations in one layer), 8 link rounds (its ~60 parity ancilla per hop fit the ~100 "
                                   "budget only at 8 rounds; 4 rounds, 1.43e4, is over budget); high = as compiled, "
                                   "2d = 4 colour classes. "
                                   "Fully serial 183,854 (F* 2.2) is the pessimistic bound, not priced",
                                   "factory.json Ch. 10 run 1 (scratchpad factories/ch10_depth_v2.py), r25")
    t_depth_2033: Tagged = Assumed((98454105.32846098, 315713175.1043646),
                                   "T-depth of the 2e3-step Gibbs shot: low = T-par colour frame, 12 link rounds (its "
                                   "~170 parity ancilla per hop fit the 750 band only at ~3 concurrent hops; 6 rounds, "
                                   "5.23e7, needs ~880); high = as compiled, 12 link rounds (250-ancilla end). Fully serial "
                                   "3.89e9 (F* 1.4) is the pessimistic bound, not priced",
                                   "factory.json Ch. 10 run 2 (scratchpad factories/ch10_depth_v2.py), r25")
    eps_first_result: Tagged = Assumed(0.30, "R1: the minimal first result at 30% statistical error",
                                       "ruling R1, H. Lamm 2026-10-02")
    chain_reuse_gain: Tagged = Assumed((3, 30), "continuing one Gibbs chain with a non-demolition N_B readout: cost per "
                                                "sample ~ tau_int ~ 2 t_rel instead of t_mix; ~10x central, 3-30x "
                                                "(conditional on fault dissipation and mid-circuit Gauss checks)",
                                       "shot audit Ch. 10 (readout_levers 4)")
    chain_reuse_steps_per_tau: Tagged = Assumed(200, "Gibbs steps per N_B autocorrelation time at the 10x central gain "
                                                     "(2e3 / 10); 0.1 faults in that window sets eps_l",
                                                "shot audit Ch. 10 (readout_levers 4)")

    # ---- 2033 headline (author ruling 2026-10-05, "Make the QSVT/TPQ route in Ch 10 as well"): QSVT filter ----
    # The construction of Ch. 9 (app06 route (iv), ch09_chiral_gauge.py): a QSVT polynomial of e^{-beta(H - mu_B N_B)/2}
    # of degree d_beta = sqrt(beta||H|| ln(1/delta_beta)) applied to a reference state, k rounds of amplitude
    # amplification calling the filter 2k + 1 times. Each block-encoding query is priced at 2-4 Trotter steps of the
    # chapter's own step (c_step, at the circuit's own R-TOL eps).
    delta_beta: Tagged = Assumed(1e-4, "QSVT filter tolerance delta_beta, as Ch. 9 (app06 delta_beta = 1e-4): "
                                       "d_beta = round(sqrt(1e2 ln 1e4)) = round(30.35) = 30",
                                 "ch09_chiral_gauge.py delta_beta; author ruling 2026-10-05")
    filter_success_prob: Tagged = Assumed(1e-2, "success probability ||e^{-beta H/2}|psi>||^2 of the filter on the "
                                                "reference (amplitude 0.1, the value Ch. 9 presumes); a typical (Haar) "
                                                "reference gives Z e^{beta E_0}/D instead, exponentially small (referee C1)",
                                          f"{TEX}:107,193 ('overlap >~1e-2'); app06 route (iv)")
    aa_call_rule: Tagged = Assumed("2k+1", "k rounds of amplitude amplification call the filter 2k + 1 times (as Ch. 9, "
                                           "verifier 2026-10-04); k = ceil(pi / (4 arcsin a)) = 8 at a = 0.1",
                                   "ch09_chiral_gauge.py aa_call_rule")
    ref_rot_per_qubit: Tagged = Assumed(1, "reference preparation: a product state, one synthesized rotation per system "
                                           "qubit (one-link / one-site marginals), applied with each filter call",
                                        "derived here; app06 route (iv)")
    qsvt_query_cost: Tagged = Stated((2, 4), f"{TEX}:107", "each block-encoding query at 2-4x the Trotter-step cost c_step")
    # records of the pre-ruling QSVT comparison (app07:104-107 before 2026-10-05)
    eps_qsp: Tagged = Stated(1e-3, f"{TEX}:107 (before 2026-10-05)", "record: d_QSP was evaluated at eps=1e-3")
    R_warm_stated: Tagged = Stated(10, f"{TEX}:107,192 (before 2026-10-05)",
                                   "record: 'R ~ 10' rounds with 2 queries per degree (20 filter-equivalents; now 2k + 1 = 17)")
    qsvt_rail_stated: Tagged = Stated((3e9, 6e9), f"{TEX}:107 (before 2026-10-05)",
                                      "record: '~3-6e9 T' (2.86e9-5.81e9 derived at R = 10, 2 queries, d = 26.3)")
    d_qsp_stated: Tagged = Stated(26, f"{TEX}:106 (before 2026-10-05)", "record: 'd_QSP ~= 26' (26.3 at eps = 1e-3)")

    # ---- logical error ---------------------------------------------------------
    faults_per_shot: Tagged = Assumed(0.1, "R3 ruling: required eps_l = 0.1 expected faults per shot, report-wide",
                                      "TRACKED_CHANGES.md RULINGS round A")
    eps_l_floor: Tagged = Cited(1e-8, "DOE_RFI_2026", "RFI floor, app07:171")

    def __post_init__(self):
        for name in ("nc", "d_2028", "L_2028", "n_stag_2028", "n_anc_2028", "d_2033", "L_2033",
                     "n_stag_2033", "n_anc_2033", "d_stretch", "L_stretch", "n_stag_stretch", "n_anc_stretch",
                     "eps_syn", "t_cap_2028", "t_ref_2033", "n_trotter_2028", "n_trotter_2028_original", "c_be_2033_superseded",
                     "gibbs_2x2_min_steps", "gibbs_2x2_t_per_shot_superseded", "beta_h_2033", "c_mix", "shots_2028", "eps_stat",
                     "nb_toy_k2_priced", "nb_toy_k2_band", "nb_toy_k2_mid", "nb_toy_quark_k2", "n_T", "n_mu", "eps_pressure_target", "eps_chi4_chi2_target",
                     "t_gate_s", "shot_overhead_s", "campaign_horizon_yr", "eps_qsp", "filter_success_prob",
                     "R_warm_stated", "qsvt_query_cost", "qsvt_rail_stated", "d_qsp_stated", "delta_beta", "ref_rot_per_qubit",
                     "faults_per_shot", "eps_l_floor", "sigma2_2028", "n_sectors_2028", "t_depth_2028",
                     "t_depth_2033", "eps_first_result", "chain_reuse_gain",
                     "chain_reuse_steps_per_tau", "gibbs_sweeps", "gibbs_eps_channel", "gibbs_beta_h_2028",
                     "gibbs_ctrl_toffoli_per_rot"):
            v = getattr(self, name)
            if v.lo <= 0:
                raise ValueError(f"{name} must be positive, got {v.value}")
            if v.lo > v.hi:
                raise ValueError(f"{name} range reversed: {v.value}")
        for name in ("group", "group_stretch"):
            if getattr(self, name).value not in GROUPS:
                raise ValueError(f"{name}={getattr(self, name).value!r} not in GROUPS")
        if self.ham.value not in PRIMCOST:
            raise ValueError(f"ham={self.ham.value!r} not in groups.PRIMCOST")
        for name in ("eps_stat", "eps_qsp", "eps_pressure_target", "eps_chi4_chi2_target", "eps_syn", "eps_first_result",
                     "delta_beta"):
            if not (0 < getattr(self, name).lo < 1):
                raise ValueError(f"{name} must be in (0,1)")
        if not (0 < self.filter_success_prob.lo < 1):
            raise ValueError(f"filter_success_prob must be in (0,1), got {self.filter_success_prob.value}")
        if self.aa_call_rule.value not in ("2k+1", "per_round"):
            raise ValueError(f"aa_call_rule={self.aa_call_rule.value!r}")
        if self.toffoli_convention.value not in ("textbook", "jones"):
            raise ValueError("toffoli_convention must name a common.T_PER_TOFFOLI key")
        for name in ("d_2028", "d_2033", "d_stretch"):
            if getattr(self, name).value not in (1, 2, 3):
                raise ValueError(f"{name} must be a spatial dimension 1-3")
        g = GROUPS[self.group.value]
        for p in MAGNETIC + ("U_FFT", "U_phi"):
            if p not in g.primitives:
                raise ValueError(f"GROUPS[{self.group.value!r}] has no compiled {p}; the chapter's per-step cost "
                                 f"cannot be built from the papers")
        if self.hop_synthesis.value not in ("rus", "watson", "rus-slope", "rs"):
            raise ValueError(f"hop_synthesis={self.hop_synthesis.value!r} is not a common.t_per_rotation model")
        if self.hop_share.value not in ("draft", "link"):
            raise ValueError(f"hop_share={self.hop_share.value!r}")
        if self.hop_mcx.value not in ("mbu", "2n-3", "draft-colour"):
            raise ValueError(f"hop_mcx={self.hop_mcx.value!r}")
        if self.hop_phasing.value not in ("hwp", "plain", "none"):
            raise ValueError(f"hop_phasing={self.hop_phasing.value!r}")
        if len(self.grid_T_mev.value) != self.n_T.lo or len(self.grid_mu_mev.value) != self.n_mu.lo:
            raise ValueError("grid_T_mev / grid_mu_mev must have n_T / n_mu points")
        if self.n_trotter_2028.lo != int(self.n_trotter_2028.lo):
            raise ValueError("n_trotter_2028 must be a whole number of steps")


# --------------------------------------------------------------------------- #
# Geometry (Eq. Nq_fd, app07:87) and the per-link costs from the papers
# --------------------------------------------------------------------------- #

def _geometry(d: int, L: int, n_stag: int, nc: int, group: str) -> dict:
    g = GROUPS[group]
    n_sites = L ** d
    n_links = d * n_sites                      # d L^d, periodic
    n_plaq = (d * (d - 1) // 2) * n_sites      # d(d-1)/2 per site
    return dict(
        n_sites=n_sites, n_links=n_links, n_plaquettes=n_plaq,
        link_qubits=g.link_qubits,                          # ACTIVE = compiled (8 for S36x3, 9 for S72x3)
        link_qubits_chapter_old=int(g.link_qubits_chapter.lo),  # what the chapter printed before R1
        link_qubits_log2=math.ceil(math.log2(g.order)),     # ceil(log2|G|), the dense minimum
        gauge_lq=n_links * g.link_qubits,
        fermion_lq=n_stag * nc * n_sites,
        hops=n_stag * n_links,                              # one hop term per link per field
    )


def _per_link(a: Assumptions, d: int, e: float, papers: bool = False) -> dict:
    """T per link per Trotter step from the papers' tables at tolerance e, rotations at the full fit (E20);
    papers=True gives the papers' slope-only price (legacy record)."""
    g, h = a.group.value, a.ham.value
    mag = magnetic_per_link(g, h, d, e, papers=papers)
    ele = electric_per_link(g, h, d, e, fft=True, papers=papers)
    ele_naive = electric_per_link(g, h, d, e, fft=False, papers=papers)
    return dict(magnetic_t_per_link=mag, electric_t_per_link=ele, gauge_t_per_link=mag + ele,
                electric_share=ele / (mag + ele), electric_naive_t_per_link=ele_naive,
                naive_over_fft=ele_naive / ele)


def _gauge_t_per_step(a: Assumptions, d: int, L: int, e: float, papers: bool = False) -> float:
    pl = _per_link(a, d, e, papers=papers)
    return d * L ** d * pl["gauge_t_per_link"]


def _hop_link(a: Assumptions, n_stag: int, e: float, **over) -> dict:
    """The staggered hop per link per Trotter step (r17): groups.hop_link_cost on the chapter's group with the
    chapter's readings (undo, sharing, MBU, HWP phasing, full-fit rotations); `over` overrides one reading."""
    kw = dict(undo=bool(a.hop_undo.value), share=a.hop_share.value, mcx=a.hop_mcx.value,
              phasing=a.hop_phasing.value)
    kw.update(over)
    return hop_link_cost(a.group.value, int(n_stag), eps=e, synthesis=a.hop_synthesis.value,
                         toffoli_convention=a.toffoli_convention.value, **kw)


def _mass_site(a: Assumptions, n_stag: int, e: float) -> dict:
    """Staggered mass per site per Trotter step (r17): one HWP group of N_c N_stag, mu_B N_B folded in."""
    return staggered_mass_site(int(a.nc.value), int(n_stag), eps=e, synthesis=a.hop_synthesis.value,
                               toffoli_convention=a.toffoli_convention.value)


def _t_per_step(a: Assumptions, d: int, L: int, n_stag: int, e: float, **over) -> float:
    """c_BE: gauge terms (papers) + staggered hop (draft, r17) on d L^d links + staggered mass on L^d sites."""
    return (_gauge_t_per_step(a, d, L, e) + d * L ** d * _hop_link(a, n_stag, e, **over)["t"]
            + L ** d * _mass_site(a, n_stag, e)["t"])


def _rot_per_step(a: Assumptions, d: int, L: int, n_stag: int) -> dict:
    """R-TOL: synthesized (arbitrary-angle) rotations per Trotter step.

    Gauge: PrimitiveCost.rot (= t_log / 1.15) x multiplicity per link (U_Tr 7, U_FFT 102, U_phi 256,
    U_inv and U_mul 0 for Sigma(36x3)). Hop: groups.fermion_hop_counts n_rot per link (colour rotations,
    diagonalizer rotations, HWP phasing). Mass: the HWP(N_c N_stag) rotations per site. The counts do not
    depend on eps, so there is no iteration.
    """
    g, h = GROUPS[a.group.value], a.ham.value
    n_links, n_sites = d * L ** d, L ** d
    gauge_link = sum(PRIMCOST[h][p](d) * g.primitives[p].rot for p in MAGNETIC)
    gauge_link += sum(m(d) * g.primitives[p].rot for p, m in ELECTRIC_PIECES)
    hop_link = _hop_link(a, n_stag, 1e-3)["n_rot"]          # counts only; eps does not enter them
    mass_site = _mass_site(a, n_stag, 1e-3)["n_rot"]
    return dict(gauge=n_links * gauge_link, hop=n_links * hop_link, mass=n_sites * mass_site,
                total=n_links * (gauge_link + hop_link) + n_sites * mass_site,
                gauge_per_link=gauge_link, hop_per_link=hop_link, mass_per_site=mass_site)


def gibbs_jump_set(d: int, L: int, n_stag: int) -> dict:
    """Referee F2: a gauge-invariant local jump set, closed under adjoint, that connects every conserved sector.
    Plaquette Tr U and its adjoint (moves electric flux); staggered hop psi†_f U psi_f and adjoint, per link and
    field; local baryon eps_abc psi_a psi_b psi_c and adjoint, per site and field (changes N_B, which H conserves);
    on-site colour-singlet psi†_f psi_f' and adjoint per field pair (changes the separately conserved field numbers)."""
    v = L ** d
    fam = dict(plaquette=2 * (d * (d - 1) // 2) * v, hop=2 * d * v * n_stag, baryon=2 * v * n_stag,
               field_mixing=2 * v * (n_stag * (n_stag - 1) // 2))
    return dict(families=fam, n_jumps=sum(fam.values()))


def gibbs_sample(a: Assumptions, d: int, L: int, n_stag: int, beta_h: float, sweeps: float, case: str) -> dict:
    """Referee F2 (2026-10-04), ESTIMATED HERE: T per Gibbs sample of the Chen-Kastoryano-Gilyen sampler
    (arxiv_2311_09207) built on the chapter's own Trotter step.

    A sampler step is one unit of Lindbladian time at ||sum_a A^a† A^a|| <= 1: one local jump drawn at random. Its
    controlled Hamiltonian-simulation time, in units of beta, is h_beta (gibbs_structure):
      OFT: e^{-iHt} before and e^{iHt} after the jump, |t| <= w beta, w = sqrt(ln 1/eps_u)       -> 2w
      B:   outer e^{-/+ i beta H t}, |t| <= t1 = ln(1/eps_u)/(2 pi); inner e^{i beta H t'} A† e^{-2i beta H t'} A
           e^{i beta H t'}, |t'| <= t2 = sqrt(ln(1/eps_u)/4)                                       -> 2 t1 + 4 t2
      low 2w (no coherent term: a floor, not exactly KMS-detailed-balanced); central 2(2w) + 2(2 t1 + 4 t2) (each
      block encoding and its inverse, one query per unit time); high central x K (the Lindbladian-simulation theorem).
    Windows are one-sided per sign: sign-magnitude time register, sign applied by negating the step's angles.
    eps_u = gibbs_eps_channel / (steps per sample). Trotter steps per beta of evolution = beta||H|| (the chapter's
    step 1/||H||, app07:92). Each step controlled: c_step + 7 T per synthesized rotation (one Toffoli), plus one
    aggregated control rotation. Steps per sample = sweeps x |A|. Rotations at the sample's own R-TOL eps.
    Per-step extras (jump operator and its inverse at <= one hop link, QFT on the time register, Gaussian filter
    state, transition-weight rotation) are carried and are < 1e-4 of the step."""
    js = gibbs_jump_set(d, L, n_stag)
    n_a = js["n_jumps"]
    n_units = sweeps * n_a
    eps_u = a.gibbs_eps_channel.lo / n_units
    lu = math.log(1.0 / eps_u)
    w, t1, t2 = math.sqrt(lu), lu / (2.0 * math.pi), math.sqrt(lu / 4.0)
    x = n_units / a.gibbs_eps_channel.lo
    k_seg = math.log(x) / math.log(math.log(x))
    h_oft, h_b = 2.0 * w, 2.0 * t1 + 4.0 * t2
    h_beta = {"low": h_oft, "central": 2.0 * h_oft + 2.0 * h_b, "high": (2.0 * h_oft + 2.0 * h_b) * k_seg}[case]
    trotter_per_step = h_beta * beta_h
    rot_step = _rot_per_step(a, d, L, n_stag)["total"] + 1                  # + the aggregated control rotation
    n_trotter = trotter_per_step * n_units
    n_rot = n_trotter * rot_step
    e = _eps(a, n_rot)
    c_step = _t_per_step(a, d, L, n_stag, e)
    t_rot = t_per_rotation(e, a.rot_synthesis_paper.value)
    c_ctrl = c_step + t_rot + a.gibbs_ctrl_toffoli_per_rot.lo * toffoli_t(1, a.toffoli_convention.value) * (rot_step - 1)
    n_bits = math.ceil(math.log2(2.0 * w * beta_h))                          # time register over [-w beta, w beta]
    aux = (2.0 * _hop_link(a, n_stag, e)["t"] + 2.0 * (n_bits * (n_bits - 1) / 2.0) * 2.0 * t_rot
           + n_bits * t_rot + 2.0 * n_bits * t_rot)
    t_step = trotter_per_step * c_ctrl + aux
    return dict(case=case, jump_families=js["families"], n_jumps=n_a, sweeps=sweeps, sampler_steps=n_units,
                eps_u=eps_u, w_oft=w, t1_b=t1, t2_b=t2, k_seg=k_seg, h_beta=h_beta,
                trotter_per_sampler_step=trotter_per_step, trotter_per_sample=n_trotter, n_rot=n_rot, eps_rot=e,
                c_step=c_step, c_step_controlled=c_ctrl, control_factor=c_ctrl / c_step, aux_t_per_step=aux,
                aux_share=aux / t_step, t_per_sampler_step=t_step, t_sample=n_units * t_step)


def gibbs_band(a: Assumptions, d: int, L: int, n_stag: int, beta_h: float, mix_cut: float = 1.0) -> dict:
    """low / central / high Gibbs sample: sweeps (gibbs_sweeps lo, c_mix, gibbs_sweeps hi) / mix_cut."""
    sw = {"low": a.gibbs_sweeps.lo, "central": a.c_mix.lo, "high": a.gibbs_sweeps.hi}
    return {k: gibbs_sample(a, d, L, n_stag, beta_h, sw[k] / mix_cut, k) for k in ("low", "central", "high")}


def qsvt_thermal(a: Assumptions, d: int, L: int, n_stag: int, beta_h: float, n_sys: int,
                 query_cost: float, n_refl: int | None = None) -> dict:
    """Author ruling 2026-10-05: the QSVT/TPQ thermal state, the construction of Ch. 9 (app06 route (iv)).

    Filter: prepares e^{-beta(H - mu_B N_B)/2}|psi> with a QSVT polynomial of degree
    d_beta = round(sqrt(beta||H|| ln(1/delta_beta))) (Ch. 9 rounds the same way), one query per degree. N_B commutes
    with H, so e^{beta mu_B N_B/2} is carried by the reference's sector weights and the polynomial is of e^{-beta H/2}
    alone; its degree does not depend on mu_B (verifier 2026-10-05: a polynomial of H - mu_B N_B would widen the
    spectrum by beta mu_B Delta N_B, ~144 at (100, 600) MeV, and d_beta would be ~47). Queries are Trotter steps, so
    the polynomial is in e^{-iH tau} with tau||H|| <~ pi and ||H||, not an LCU lambda, sets d_beta. Amplification: amplitude a = sqrt(filter_success_prob),
    k = ceil(pi / (4 arcsin a)) rounds, 2k + 1 filter calls (Ch. 9's rule). A query is query_cost x c_step (Stated
    2-4, the chapter's own block-encoding cost per query) at the circuit's own R-TOL eps. Per filter call also: the
    product reference (one rotation per system qubit) and d_beta + 1 QSVT signal rotations; per round one fixed-point
    phase rotation and one reflection about the start state A S_0 A^dag. S_0 acts on every qubit A touches (system,
    BE and parity ancilla, flag): C^{n_refl + 1}Z at n_refl Toffolis (MBU), n_refl = the full register (default
    n_sys; the 2033 headline passes the high end, 1014). The reflection about
    the success flag is a Z (Clifford). Readout: N_B is diagonal in the Jordan-Wigner occupation basis, so it is read by
    measuring the fermion register with the flag (no T)."""
    p = a.filter_success_prob.lo
    amp = math.sqrt(p)
    theta = math.asin(amp)
    k = math.ceil(math.pi / (4.0 * theta))
    calls = 2 * k + 1 if a.aa_call_rule.value == "2k+1" else k
    d_exact = math.sqrt(beta_h * math.log(1.0 / a.delta_beta.lo))
    d_b = round(d_exact)
    queries = calls * d_b
    step_eq = queries * query_cost
    rot_step = _rot_per_step(a, d, L, n_stag)["total"]
    n_rot_ref = calls * n_sys * a.ref_rot_per_qubit.lo
    n_rot_sig = calls * (d_b + 1) + k
    n_rot = step_eq * rot_step + n_rot_ref + n_rot_sig
    e = _eps(a, n_rot)
    c_step = _t_per_step(a, d, L, n_stag, e)
    t_rot = t_per_rotation(e, a.rot_synthesis_paper.value)
    refl_toff = mcx_toffolis((n_sys if n_refl is None else n_refl) + 1, "mbu")
    t_filter = step_eq * c_step
    t_ref = n_rot_ref * t_rot
    t_sig = n_rot_sig * t_rot
    t_refl = k * toffoli_t(refl_toff, a.toffoli_convention.value)
    t = t_filter + t_ref + t_sig + t_refl
    return dict(success_prob=p, amplitude=amp, theta=theta, aa_rounds=k, filter_calls=calls, d_beta_exact=d_exact,
                d_beta=d_b, queries=queries, query_cost=query_cost, step_equivalents=step_eq, rot_per_step=rot_step,
                n_rot_ref=n_rot_ref, n_rot_signal=n_rot_sig, n_rot=n_rot, eps_rot=e, c_step=c_step, t_per_rotation=t_rot,
                reflection_toffolis=refl_toff, t_filter=t_filter, t_reference=t_ref, t_signal=t_sig, t_reflections=t_refl,
                t_readout=0.0, overhead_share=(t_ref + t_sig + t_refl) / t, t_shot=t)


def qsvt_band(a: Assumptions, d: int, L: int, n_stag: int, beta_h: float, n_sys: int,
              n_refl: int | None = None) -> dict:
    """low / high QSVT shot at the two ends of the stated query cost (2, 4 c_step)."""
    lo, hi = a.qsvt_query_cost.value
    return {"low": qsvt_thermal(a, d, L, n_stag, beta_h, n_sys, lo, n_refl),
            "high": qsvt_thermal(a, d, L, n_stag, beta_h, n_sys, hi, n_refl)}


def _wall_shot_at_baseline(a: Assumptions, t_shot: float, d_hi: float) -> float:
    """R9: one shot on ONE machine at the 10-factory baseline (1 us per T): max(N_T t_gate, D_T t_r) + overhead.
    Equals the serial N_T t_gate + overhead whenever F* = N_T / D_T >= 10; otherwise the shot is depth-limited
    (the corrected wall of CONTRACT rule 12). D_T is the high end of the factory.json band."""
    return max(t_shot * a.t_gate_s.lo, d_hi * REACTION_TIME_S) + a.shot_overhead_s.lo


def _eps(a: Assumptions, n_rot: float) -> float:
    """R-TOL: the circuit's own per-rotation tolerance, common.eps_rot_for(N_rot, eps_syn)."""
    return eps_rot_for(n_rot, a.eps_syn.lo)


def _hop_primitives(a: Assumptions, n_stag: int, n_links: int, n_sites: int, n_steps: float,
                    e: float) -> tuple[Primitive, ...]:
    """The r17 hop breakdown over n_links x n_steps (one Primitive per groups.fermion_hop_counts item, t_each =
    its T per link) and the staggered mass over n_sites x n_steps (HWP rotations and Toffolis)."""
    hl = _hop_link(a, n_stag, e)
    prims = []
    for it in hl["items"]:
        name = it["name"] if it["name"].startswith("hop_") else f"hop_{it['name']}"
        prims.append(Primitive(name, n_links * n_steps, it["t"], it["status"], it["src"],
                               f"per link-step: {it['toffoli']} Toffoli + {it['t_direct']} T + {it['n_rot']} rotations "
                               f"at {hl['t_rot']:.4f} T (full fit, eps={e:.4g}, R-TOL); {it['note']}"))
    ms = _mass_site(a, n_stag, e)
    prims.append(Primitive("mass_hwp_rotations", n_sites * n_steps * ms["n_rot"], ms["t_rot"], CircuitStatus.SCALING,
                           ms["src"], f"staggered mass: one HWP group of k={ms['k']} per site-step, {ms['n_rot']} "
                           f"rotations at the full fit (eps={e:.4g}); mu_B N_B folded into the same angle"))
    prims.append(Primitive("mass_hwp_toffolis", n_sites * n_steps * ms["n_toffoli"],
                           toffoli_t(1, a.toffoli_convention.value), CircuitStatus.COMPILED, HWP_SRC,
                           f"staggered mass: {ms['n_toffoli']} adder Toffolis per k={ms['k']} group (MBU)"))
    return tuple(prims)


def _hop_record(a: Assumptions, n_stag: int, e: float, n_links: int) -> dict:
    """Hop intermediates shared by both eras: the priced hop, its pieces, the conservative (share="draft"),
    no-undo and draft-literal variants at the same eps, and the retired R-HOP price for comparison."""
    hl = _hop_link(a, n_stag, e)
    by = {it["name"]: it["t"] for it in hl["items"]}
    drf = _hop_link(a, n_stag, e, share="draft")["t"]        # r17 per-application recompute (sensitivity)
    no_undo = _hop_link(a, n_stag, e, undo=False)["t"]
    fp_lit = 2 * fp_printed(a.group.value, "C_G", e) + n_stag * fp_printed(a.group.value, "C_hop_text", e)
    return dict(
        hop_t_per_link=hl["t"], hop_toffoli_per_link=hl["toffoli"], hop_t_direct_per_link=hl["t_direct"],
        hop_n_rot_per_link=hl["n_rot"], hop_t_per_rotation=hl["t_rot"], hop_pieces_t_per_link=by,
        hop_colour_frame_share=(by["colour_squish_and_parity"] + by["colour_rotations"]) / hl["t"],
        hop_t_per_link_share_draft=drf, hop_t_per_step_share_draft=n_links * drf,
        hop_t_per_link_no_undo=no_undo,
        hop_t_per_link_fp_printed=fp_lit,
        hop_t_per_link_rhop_legacy=hl["legacy"]["rhop_rus_slope"],
        hop_over_rhop_legacy=hl["legacy"]["ratio_new_over_rhop_rus_slope"],
    )


def _gauge_t_per_step_report_convention(a: Assumptions, d: int, L: int, e: float) -> float:
    """The same count through PrimitiveCost.t_report (T4 fit, named Toffoli convention); equal to
    _gauge_t_per_step at 7 T per Toffoli since E20."""
    g = GROUPS[a.group.value]
    h, conv = a.ham.value, a.toffoli_convention.value
    per_link = sum(PRIMCOST[h][p](d) * g.primitives[p].t_report(e, conv) for p in MAGNETIC)
    per_link += sum(m(d) * g.primitives[p].t_report(e, conv) for p, m in ELECTRIC_PIECES)
    return d * L ** d * per_link


def _gauge_primitives(a: Assumptions, d: int, n_links: int, n_steps: float,
                      e: float) -> tuple[Primitive, ...]:
    """One Primitive per compiled gauge primitive, multiplied through links and steps."""
    g = GROUPS[a.group.value]
    h = a.ham.value
    prims = []
    for p in MAGNETIC:
        pc = g.primitives[p]
        prims.append(Primitive(f"magnetic_{p}", n_steps * n_links * PRIMCOST[h][p](d), pc.t(e), pc.status, pc.src,
                               f"{PRIMCOST[h][p](d):g} per link per step (tab:primcost); {pc.t_const:g} + {pc.n_rot:g} x (1.15 log2(1/eps) + 9.2) "
                               f"at eps={e:.4g} (R-TOL, E20 full fit)"))
    for p, m in ELECTRIC_PIECES:
        pc = g.primitives[p]
        prims.append(Primitive(f"electric_{p}", n_steps * n_links * m(d), pc.t(e), pc.status, pc.src,
                               f"{m(d):g} per link per step; {pc.t_const:g} + {pc.n_rot:g} x (1.15 log2(1/eps) + 9.2) at eps={e:.4g} (R-TOL, E20)"))
    return tuple(prims)


# --------------------------------------------------------------------------- #
# Shots for chi4/chi2 (referee F1, 2026-10-04): delta-method variance of the k-statistic ratio
# --------------------------------------------------------------------------- #
# The estimator is r = k4/k2 (k-statistics of the N_B samples). Its asymptotic variance from N independent
# samples is V1/N with V1 = E[IF_r^2], IF_r = (IF_k4 - r IF_k2)/k2, IF_k2 = (x-m)^2 - mu2,
# IF_k4 = (x-m)^4 - mu4 - 4 mu3 (x-m) - 6 mu2 IF_k2. This is the delta method on the moments up to order 8 and
# equals the Kendall-Stuart cumulant form (kappa_8 + 16 k2 k6 + 48 k3 k5 + 34 k4^2 + 72 k2^2 k4 + 144 k2 k3^2
# + 24 k2^4 for k4, with the k2 variance and the k2-k4 covariance). A Gaussian has kappa_8 = 0 but
# Var(k4) = 24 kappa_2^4 / N. Reweighted samples (self-normalized, balance heuristic over an equal
# mixture of ensembles): V1 = E_target[(p_target / p_mixture) IF^2], which carries the effective-sample-size loss.
# Each shot is an independent preparation, so there is no autocorrelation factor.
# Toy N_B distributions (the true one at V = 2^3 is open): Skellam in baryons (kappa_2(0) = k, every even cumulant
# k cosh(mu/T), odd k sinh(mu/T), chi4/chi2 = 1) or in quarks (N_B = (n - nbar)/3, mu_q = mu_B/3, ratio 1/9),
# with kappa_2(0) the same at all five temperatures. Checked against Monte Carlo of the k-statistic ratio
# (scratchpad referee/check.py, mcrw.py).

_NMAX = 150
_SERIES_K = 60


def _log_bessel_i(n: int, z: float) -> float:
    """log I_|n|(z) by its power series (z <= ~1 here)."""
    n = abs(n)
    terms = [(2 * k + n) * math.log(z / 2.0) - math.lgamma(k + 1) - math.lgamma(k + n + 1) for k in range(_SERIES_K)]
    mx = max(terms)
    return mx + math.log(sum(math.exp(t - mx) for t in terms))


_BASE_CACHE: dict = {}


def _nb_base(kind: str, k2: float):
    key = (kind, k2)
    if key not in _BASE_CACHE:
        b = k2 / 2.0 if kind == "baryon" else 9.0 * k2 / 2.0
        ns = list(range(-_NMAX, _NMAX + 1))
        _BASE_CACHE[key] = (ns, [_log_bessel_i(n, 2.0 * b) for n in ns], b)
    return _BASE_CACHE[key]


def nb_toy_logpmf(kind: str, k2: float, u: float) -> tuple[list[float], list[float]]:
    """(N_B values, log-probabilities) of the toy at mu_B/T = u. kind 'baryon' or 'quark'."""
    if kind not in ("baryon", "quark"):
        raise ValueError(kind)
    ns, li, b = _nb_base(kind, k2)
    uu = u if kind == "baryon" else u / 3.0
    lp = [-2.0 * b * math.cosh(uu) + n * uu + l for n, l in zip(ns, li)]
    mx = max(lp)
    lz = mx + math.log(sum(math.exp(x - mx) for x in lp))
    xs = [float(n) if kind == "baryon" else n / 3.0 for n in ns]
    return xs, [x - lz for x in lp]


def ratio_variance_per_sample(kind: str, k2: float, u: float, ens_u: tuple | None = None) -> dict:
    """V1 for r = k4/k2 (and the relative per-sample variance of n_B) at mu_B/T = u, sampled directly
    (ens_u None) or reweighted from an equal mixture of ensembles at mu_B/T in ens_u."""
    xs, lt = nb_toy_logpmf(kind, k2, u)
    pt = [math.exp(x) for x in lt]
    m = sum(p * x for p, x in zip(pt, xs))
    y = [x - m for x in xs]
    mu = [sum(p * yy ** k for p, yy in zip(pt, y)) for k in range(5)]
    kap2, kap4 = mu[2], mu[4] - 3 * mu[2] ** 2
    r = kap4 / kap2
    # pw[i] = p_target x (p_target / p_sample), formed in logs (the weight alone can overflow far in the tail)
    if ens_u is None:
        pw = pt
    else:
        lens = [nb_toy_logpmf(kind, k2, ue)[1] for ue in ens_u]
        pw = []
        for i, l_t in enumerate(lt):
            ls = [le[i] for le in lens]
            mx = max(ls)
            lmix = mx + math.log(sum(math.exp(v - mx) for v in ls)) - math.log(len(ls))
            e = 2.0 * l_t - lmix
            pw.append(math.exp(e) if e > -700.0 else 0.0)
    v1 = 0.0
    vnb = 0.0
    for q, yy in zip(pw, y):
        if q == 0.0:
            continue
        i2 = yy ** 2 - mu[2]
        i4 = yy ** 4 - mu[4] - 4 * mu[3] * yy - 6 * mu[2] * i2
        v1 += q * ((i4 - r * i2) / kap2) ** 2
        vnb += q * yy ** 2
    return dict(v1=v1, r=r, kappa2=kap2, kappa4=kap4, mean=m,
                v1_nb_rel=(vnb / m ** 2) if m != 0 else math.inf)


def _grid_mu_T(a: "Assumptions"):
    return [(T, mu) for T in a.grid_T_mev.value for mu in a.grid_mu_mev.value]


def grid_shots_direct(a: "Assumptions", kind: str, k2: float, delta_abs: float | None = None,
                      rel: float | None = None) -> dict:
    """Shots to reach |delta(chi4/chi2)| <= delta_abs (or rel x |r|) at every grid point, sampled directly."""
    per = {}
    for T, mu in _grid_mu_T(a):
        v = ratio_variance_per_sample(kind, k2, mu / T)
        d = delta_abs if delta_abs is not None else rel * abs(v["r"])
        per[(T, mu)] = v["v1"] / d ** 2
    return dict(per_point=per, total=sum(per.values()), min=min(per.values()), max=max(per.values()))


def row_samples(a: "Assumptions", kind: str, k2: float, T: float, delta_rel: float, ens_mu: tuple | None,
                with_nb: bool = True) -> float:
    """Samples for one T row: chi4/chi2 at delta_rel at every mu point and (with_nb) n_B at delta_rel for
    mu_B >= the first nonzero point; direct = sum over points, reweighted = max over points."""
    need = []
    for mu in a.grid_mu_mev.value:
        v = ratio_variance_per_sample(kind, k2, mu / T, None if ens_mu is None else tuple(e / T for e in ens_mu))
        n = v["v1"] / (delta_rel * abs(v["r"])) ** 2
        if with_nb and mu > 0:
            n = max(n, v["v1_nb_rel"] / delta_rel ** 2)
        need.append(n)
    return sum(need) if ens_mu is None else max(need)


# --------------------------------------------------------------------------- #
# Eras
# --------------------------------------------------------------------------- #

def _model_2028(a: Assumptions) -> Result:
    d, L = a.d_2028.value, a.L_2028.value
    geo = _geometry(d, L, a.n_stag_2028.value, a.nc.value, a.group.value)
    anc = _rng(a.n_anc_2028)
    lq_system = geo["gauge_lq"] + geo["fermion_lq"]
    lq = (lq_system + anc[0], lq_system + anc[1])

    cap = a.t_cap_2028.lo
    n_steps = a.n_trotter_2028.lo                                  # 1
    ns28 = a.n_stag_2028.value
    rot = _rot_per_step(a, d, L, ns28)                             # 3,708 gauge + 7,192 hop + 8 mass = 10,908
    n_rot_shot = n_steps * rot["total"]                            # R-TOL: 10,908 synthesized rotations per shot
    e = _eps(a, n_rot_shot)                                        # sqrt(1e-2 / 10,908) = 9.575e-4
    pl = _per_link(a, d, e)
    gauge_step = geo["n_links"] * pl["gauge_t_per_link"]          # the derived c_BE for the gauge terms
    hl = _hop_link(a, ns28, e)                                     # r19 hop (share="link"), 36,902.7 T/link
    hop_step = geo["n_links"] * hl["t"]
    ms = _mass_site(a, ns28, e)                                    # HWP(3) per site
    mass_step = geo["n_sites"] * ms["t"]
    step = gauge_step + hop_step + mass_step                       # 399,956.8 (4.0x the cap)
    step_fixed_1e4 = _t_per_step(a, d, L, ns28, 1e-4)              # record at a fixed 1e-4
    steps_under_cap = math.floor(cap / step)                       # 0: no step fits
    t_step_total = n_steps * step
    t_shot = (t_step_total, t_step_total)
    shots = _rng(a.shots_2028)
    wall_shot = t_step_total * a.t_gate_s.lo + a.shot_overhead_s.lo   # 1 us per T + ~0.1 ms per shot
    wall = (shots[0] * wall_shot, shots[1] * wall_shot)

    # the ten-step ramp as originally scoped, and the Gibbs bound in sampler steps (item 10).
    # R-TOL: one shot of the ramp is its own circuit, N_rot = 10 x 10,908 = 109,080, priced at its own
    # eps_rot = sqrt(1e-2 / 109,080) = 3.028e-4 (420,792.3 T/step, 4.21e6 T in total).
    ramp_steps = a.n_trotter_2028_original.lo
    ramp_n_rot = ramp_steps * rot["total"]
    ramp_eps = _eps(a, ramp_n_rot)
    ramp_step_t = _t_per_step(a, d, L, ns28, ramp_eps)
    ramp_t = ramp_steps * ramp_step_t
    gibbs_steps = a.gibbs_2x2_min_steps.lo
    # referee F2 (2026-10-04): the Gibbs sample at V=2^2 priced as for 2033 (32 jumps, beta||H|| = 1e2 Assumed)
    gb28 = gibbs_band(a, d, L, ns28, a.gibbs_beta_h_2028.lo)
    # author ruling 2026-10-05: the 2033 route (QSVT/TPQ filter) at V=2^2, beta||H|| = 1e2 (Assumed, as for the sampler):
    # why the 2028 benchmark stays a ground-state ramp
    qb28 = qsvt_band(a, d, L, ns28, a.gibbs_beta_h_2028.lo, geo["gauge_lq"] + geo["fermion_lq"], lq[1])

    # ancilla budget (ruling ch10-ancilla-contents A), re-sized for the r17 hop: the paper workspace plus the
    # hop's HWP(2 N_c N_stag) and largest MBU ladder (C^7X in the Sigma(36x3) colour squish: 5 AND ancillas),
    # the draft's 24-qubit Sigma(36x3) squish scratch (su3_diag.tex:103), and the mass HWP(N_c N_stag). The draft's
    # separate flag and parity registers are not sized (NEEDS_AUTHOR); the 24 may already hold them.
    # r19 (share="link"): the squish is held through the hop, so under serial reuse the 24 sit beside the
    # largest piece the hop runs while it is held (its HWP(2 N_c N_stag) or its MBU ladder).
    # r22 (E26): a phasing group of k holds k - w(k) ancilla (HWP(6) 4, HWP(3) 1; were 7 and 3), and
    # repeat-until-success synthesis holds common.RUS_ANCILLA = 1 while a group's rotations are synthesized.
    # HWP(6) + RUS = 5 = the C^7X ladder, so the serial floor is 24 + 5 = 29 with or without the synthesis ancilla.
    prim = GROUPS[a.group.value].primitives
    k_hop = 2 * int(a.nc.value) * ns28
    ws = {"U_FFT": prim["U_FFT"].ancilla, "U_Tr": prim["U_Tr"].ancilla, "U_inv": prim["U_inv"].ancilla,
          "U_mul": prim["U_mul"].ancilla, f"hop HWP({k_hop})": hwp_ancilla(k_hop),
          "hop C^7X ladder (MBU)": mcx_toffolis(7, "mbu") - 1,
          "hop colour-squish scratch (FP su3_diag tab:su3squishcosts)": SQUISH_ANCILLA_S36,
          f"mass HWP({ms['k']})": ms["ancilla"]}
    lq_sys = geo["gauge_lq"] + geo["fermion_lq"]
    ws_serial = max(ws.values())
    if a.hop_share.value == "link":
        ws_serial = max(ws_serial, SQUISH_ANCILLA_S36 + max(hwp_ancilla(k_hop), mcx_toffolis(7, "mbu") - 1))
    ws_serial_rus = max(ws_serial, hwp_ancilla(k_hop) + RUS_ANCILLA)
    if a.hop_share.value == "link":
        ws_serial_rus = max(ws_serial_rus, SQUISH_ANCILLA_S36 + hwp_ancilla(k_hop) + RUS_ANCILLA)
    alg_reg = (lq_sys + ws_serial + 1, lq_sys + sum(ws.values()) + 1)

    prims = _gauge_primitives(a, d, geo["n_links"], n_steps, e) + _hop_primitives(
        a, ns28, geo["n_links"], geo["n_sites"], n_steps, e)
    eps_l_req = a.faults_per_shot.lo / t_step_total

    # r25 (R1, R9): depth exports, the depth-corrected wall (F* = 9.9 at the as-compiled end), the 2-3 sectors.
    dx = depth_exports(t_shot, a.t_depth_2028.value, shots, a.t_gate_s, a.shot_overhead_s)
    wall_shot_bl = _wall_shot_at_baseline(a, t_step_total, a.t_depth_2028.hi)    # 0.40499 s (depth-limited)
    wall_sector = (shots[0] * wall_shot_bl, shots[1] * wall_shot_bl)               # 810.0 s = 13.5 min
    n_sec = _rng(a.n_sectors_2028)
    wall_campaign = (n_sec[0] * wall_sector[0], n_sec[1] * wall_sector[1])        # 1,620-2,430 s = 27-40 min
    eps_sector = (math.sqrt(a.sigma2_2028.lo / shots[1]), math.sqrt(a.sigma2_2028.lo / shots[0]))   # 2.236%
    inter = dict(
        n_sites=geo["n_sites"], n_links=geo["n_links"], n_plaquettes=geo["n_plaquettes"],
        link_qubits=geo["link_qubits"], link_qubits_chapter_old=geo["link_qubits_chapter_old"],
        link_qubits_log2=geo["link_qubits_log2"],
        gauge_lq=geo["gauge_lq"], fermion_lq=geo["fermion_lq"], ancilla_lq=anc,
        lq_system=lq_system, lq_total=lq,
        hops=geo["hops"],
        **pl,
        c_be_gauge_t_per_step=gauge_step,
        eps_syn=a.eps_syn.lo, n_rot_per_step=rot["total"], n_rot_gauge_per_step=rot["gauge"],
        n_rot_hop_per_step=rot["hop"], n_rot_per_shot=n_rot_shot, eps_rot=e,
        t_per_rotation=t_per_rotation(e, a.rot_synthesis_paper.value),
        synthesis_error_per_shot=n_rot_shot * e ** 2,
        t_per_shot_at_fixed_1em4_superseded=n_steps * step_fixed_1e4,
        c_be_gauge_t_per_step_report_convention=_gauge_t_per_step_report_convention(a, d, L, e),
        c_be_gauge_t_per_step_at_paper_fiducial_eps=geo["n_links"] * (
            magnetic_per_link(a.group.value, a.ham.value, d, 1e-8)
            + electric_per_link(a.group.value, a.ham.value, d, 1e-8)),
        **_hop_record(a, ns28, e, geo["n_links"]),
        hop_t_per_step=hop_step, hop_share_of_step=hop_step / step,
        hop_t_per_step_rhop_legacy=geo["n_links"] * hl["legacy"]["rhop_rus_slope"],
        n_rot_mass_per_step=rot["mass"], mass_hwp_k=ms["k"], mass_t_per_site=ms["t"], mass_t_per_step=mass_step,
        hop_allowance_t_per_step_superseded=cap - _gauge_t_per_step(a, d, L, 1e-4, papers=True),   # record, old 1e-4, slope-only
        c_be_gauge_t_per_step_papers_convention=_gauge_t_per_step(a, d, L, e, papers=True),   # legacy record (R2)
        c_be_t_per_step=step, cap_fraction=t_step_total / cap, fits_cap=t_step_total <= cap,
        c_be_t_per_step_share_draft=gauge_step + geo["n_links"] * _hop_link(a, ns28, e, share="draft")["t"] + mass_step,
        t_cap=cap, n_trotter=n_steps, steps_under_cap=steps_under_cap,
        ancilla_workspace=ws, ancilla_workspace_max=max(ws.values()), ancilla_workspace_serial=ws_serial,
        ancilla_workspace_sum=sum(ws.values()),
        ancilla_rus=RUS_ANCILLA, ancilla_workspace_serial_with_rus=ws_serial_rus,   # 1, 29 (record: the floor does not move)
        algorithmic_register=alg_reg,
        t_per_shot_value=t_step_total, ramp_original_steps=a.n_trotter_2028_original.lo, ramp_original_t=ramp_t,
        ramp_original_n_rot=ramp_n_rot, ramp_original_eps_rot=ramp_eps, ramp_original_t_per_step=ramp_step_t,
        ramp_over_cap=ramp_t / cap,
        shots=shots, wall_per_shot_s=wall_shot, wall_total_s=wall, shot_overhead_s=a.shot_overhead_s.lo,
        faults_per_shot_budget=a.faults_per_shot.lo, eps_l_required=eps_l_req,
        faults_at_rfi_floor=t_step_total * a.eps_l_floor.lo,
        gibbs_2x2_min_steps=gibbs_steps, gibbs_2x2_t_lower_bound=gibbs_steps * step,
        gibbs_2x2_t_per_shot_superseded=a.gibbs_2x2_t_per_shot_superseded.lo,
        gibbs_2x2_sampler=gb28, gibbs_2x2_n_jumps=gb28["central"]["n_jumps"],
        gibbs_2x2_t_sample=(gb28["low"]["t_sample"], gb28["central"]["t_sample"], gb28["high"]["t_sample"]),
        gibbs_2x2_t_per_sampler_step=gb28["central"]["t_per_sampler_step"],
        gibbs_2x2_over_cap=gb28["central"]["t_sample"] / cap,
        qsvt_2x2=qb28, qsvt_2x2_t=(qb28["low"]["t_shot"], qb28["high"]["t_shot"]),
        qsvt_2x2_over_cap=(qb28["low"]["t_shot"] / cap, qb28["high"]["t_shot"] / cap),
        # r25 R9 exports (names fixed by CONTRACT rule 12) and records
        t_per_shot=dx["t_per_shot"], t_depth_per_shot=dx["t_depth_per_shot"], f_star=dx["f_star"],
        floor_wall_s=dx["floor_wall_s"], factories_for_1yr=dx["factories_for_1yr"],
        wall_first_result_s=None, wall_campaign_s=wall_campaign,
        wall_serial_s=dx["wall_serial_s"], baseline_ok=dx["baseline_ok"], fits_1yr=dx["fits_1yr"],
        wall_per_shot_baseline_s=wall_shot_bl, wall_sector_s=wall_sector, n_sectors=n_sec,
        sigma2_per_shot=a.sigma2_2028.lo, eps_per_sector=eps_sector,
        # referee F3 / G2: the 2.2% is absolute, in units of the normalization lambda of the measured terms; the
        # onset mu_c = E0(N+1) - E0(N) carries sqrt(2) of it. The one step is a lower bound on the executed circuit
        # (an N-step ramp plus the grouped energy readout; ramp length open).
        eps_mu_c_over_lambda=(math.sqrt(2.0) * eps_sector[0], math.sqrt(2.0) * eps_sector[1]),
        hard_ops_is_lower_bound=True,
    )
    return Result(
        era="2028", lq=lq, hard_ops=t_shot, breakdown=prims,
        intermediates=inter, shots=shots, wall_time_s=wall_sector, epsilon_l=(eps_l_req, eps_l_req),
        notes=(
            "LQ = 64 gauge (8 links x 8, compiled Sigma(36x3) register, R1) + 12 fermion (1 x 3 x 4) + ~100 ancilla "
            "= 176, box '~180'.",
            "R-TOL: N_rot = 3,708 gauge (U_Tr 7 x 4, U_FFT 102 x 16, U_phi 256 x 8) + 7,192 hop (899 per link x 8) "
            "+ 8 mass (HWP(3) 2 x 4) = 10,908 per shot; eps_rot = sqrt(1e-2 / 10,908) = 9.575e-4: 20.73 T per "
            "rotation, gauge, hop and mass alike (full fit, E20).",
            "Gauge terms (groups.py, arXiv:2405.05973 tab:primcost + arXiv:2408.00075 FFT, full fit): 13,068 T per "
            "link x 8 = 104,541 T/step; electric share 81%. Papers' slope-only record 70,427.",
            "Hop (r19, groups.hop_link_cost, share='link' + undo): 2,592 Toffoli + 120 T + 899 rotations = 36,902.7 "
            "T/link; colour frame 80% of it; x 8 = 295,222 T/step. Sensitivity share='draft' 85,118.7 T/link "
            "(step 785,685, 7.9x); no undo 30,102; the draft's printed 2 C^G + C^hop 3.2e4; R-HOP was 676.",
            "Mass: HWP(3) per site, 2 rotations + 1 Toffoli = 48.47 T; x 4 = 193.9 T/step; mu_B N_B folded in.",
            "Step = 104,541 + 295,222 + 194 = 399,957 T = 4.0x the 1e5 cap: no step fits (floor = 0). hard_ops is "
            "the one priced step, a LOWER BOUND on the executed circuit (referee G2): box '>~4.0e5'. A ten-step ramp "
            "(its own circuit, N_rot 109,080, eps 3.03e-4) is 4.2e6 T; the ramp length and the energy readout are open.",
            "Thermal state at V=2^2 (author ruling 2026-10-05): the 2033 QSVT/TPQ filter at beta||H|| = 1e2 (Assumed), "
            "17 calls of degree 30, 2-4 steps per query, is 4.7-9.6e8 T, 4.7-9.6e3x the cap, so the benchmark stays a "
            "ground-state ramp. Records: the Gibbs sample (referee F2) 1.3e12 T (1.2e11-1.3e13), 1.3e7x the cap; the old "
            "'>~14 steps' (5.6e6 T).",
            "Shots ~2e3 are Stated with no derivation. The ~100 ancilla are a sized budget (ruling "
            "ch10-ancilla-contents A): workspace 29 (serial reuse; the 24-qubit squish scratch held beside the "
            "C^7X ladder 5; the hop HWP(6) 4 with 1 synthesis ancilla is also 5) to 55 (none: squish scratch 24, "
            "U_FFT 8, U_Tr 7, U_inv 4, U_mul 2, hop HWP(6) 4, C^7X ladder 5, mass HWP(3) 1; phasing groups at "
            "k - w(k), E26) + 1 Hadamard-test qubit; algorithmic register 106-132 of 176. The "
            "draft's separate flag and parity registers are not sized.",
            "epsilon_l (R3): 0.1 expected faults in 399,957 ops -> 2.50e-7; the RFI floor 1e-8 suffices. "
            "Wall 0.40 s/shot serial (incl. 0.1 ms overhead), 800 s for 2e3 shots. r25 (R9): T-depth 2.67e4-4.05e4, "
            "F* 9.9-15; at the as-compiled end the shot is depth-limited, 0.405 s, so 810 s = 13.5 min per sector "
            "(box '~13.5 min') and 27-40 min for the 2-3 sectors. 2e3 shots at sigma^2 <= 1 give 2.2% of lambda "
            "(absolute) per sector, 3.2% of lambda on mu_c (referee F3).",
        ),
    )


def _model_2033(a: Assumptions) -> Result:
    d, L = a.d_2033.value, a.L_2033.value
    geo = _geometry(d, L, a.n_stag_2033.value, a.nc.value, a.group.value)
    anc_lo, anc_hi = _rng(a.n_anc_2033)
    lq_system = geo["gauge_lq"] + geo["fermion_lq"]                # 192 + 72 = 264
    lq = (lq_system + anc_lo, lq_system + anc_hi)                  # 514 - 1014

    beta_h = a.beta_h_2033.lo
    n_steps = a.c_mix.lo * beta_h                                  # 20 x 1e2 = 2e3
    ns33 = a.n_stag_2033.value
    rot = _rot_per_step(a, d, L, ns33)                             # 11,208 gauge + 64,632 hop + 32 mass = 75,872
    n_rot_shot = n_steps * rot["total"]                            # R-TOL: 1.517e8 per shot
    e = _eps(a, n_rot_shot)                                        # sqrt(1e-2 / 1.517e8) = 8.118e-6
    pl = _per_link(a, d, e)
    gauge_step = geo["n_links"] * pl["gauge_t_per_link"]          # compiled gauge terms
    hl = _hop_link(a, ns33, e)                                     # r19 hop (share="link"), 95,734 T/link
    hop_step = geo["n_links"] * hl["t"]
    ms = _mass_site(a, ns33, e)                                    # HWP(9) per site
    mass_step = geo["n_sites"] * ms["t"]
    c_be = gauge_step + hop_step + mass_step                       # 2.760e6, derived here (was Stated 5e5)
    c_be_old = a.c_be_2033_superseded.lo
    t_hsim_bound = n_steps * c_be                                  # 5.521e9: the pre-F2 '>~5.5e9' (record)
    # referee F2 (2026-10-04), RECORD since the author ruling of 2026-10-05: the Chen-Kastoryano-Gilyen sampler priced
    # per step on the chapter's Trotter step. 288 jumps x 20 sweeps = 5,760 sampler steps, each 3,758 controlled steps.
    gb = gibbs_band(a, d, L, ns33, beta_h)
    gc = gb["central"]
    g_shot = gc["t_sample"]                                        # 8.39e13 central (7.6e12 - 9.5e14)
    g_band = (gb["low"]["t_sample"], gb["high"]["t_sample"])
    gxk = gibbs_sample(a, d, L, ns33, beta_h, a.c_mix.lo, "high")
    g_xk = gxk["t_sample"]                                         # 4.42e14 (K = 5.13 at 20 sweeps)
    g_range = (g_shot, g_xk)                                       # the pre-ruling box '8.4e13--4.4e14'
    c_be_fixed_1e4 = _t_per_step(a, d, L, ns33, 1e-4)              # record at a fixed 1e-4

    # HEADLINE (author ruling 2026-10-05): the QSVT/TPQ thermal state, Ch. 9's construction. 17 filter calls (8 rounds at
    # amplitude 0.1) of degree 30, 510 queries at 2-4 c_step each, plus reference, signal rotations and reflections.
    qb = qsvt_band(a, d, L, ns33, beta_h, lq_system, lq[1])        # reflections on the full register (1014)
    q_lo, q_hi = qb["low"], qb["high"]
    t_range = (q_lo["t_shot"], q_hi["t_shot"])                     # 2.77e9 - 5.63e9
    t_shot = t_range[0]                                            # the low end: breakdown and scalar records

    # shots (Eq. Nshot_fd, referee F1): sum over the 25 points of V1(T, mu) / delta^2, V1 the delta-method
    # variance of k4/k2 for the priced toy (baryon Skellam, kappa_2(0) = 0.01), delta = 0.1 absolute
    n_grid = a.n_T.lo * a.n_mu.lo                                  # 25
    k2_lo, k2_hi = _rng(a.nb_toy_k2_band)
    g_pr = grid_shots_direct(a, "baryon", a.nb_toy_k2_priced.lo, delta_abs=a.eps_stat.lo)
    shots = g_pr["total"]                                          # 2.576e5
    shots_point = shots / n_grid                                   # 1.03e4 average per point
    g_mid = grid_shots_direct(a, "baryon", a.nb_toy_k2_mid.lo, delta_abs=a.eps_stat.lo)     # 9.08e5
    g_hi = grid_shots_direct(a, "baryon", k2_hi, delta_abs=a.eps_stat.lo)                 # 2.24e6
    g_q_rel = grid_shots_direct(a, "quark", a.nb_toy_quark_k2.lo, rel=a.eps_chi4_chi2_target.lo)   # 3.88e5
    g_q_abs = grid_shots_direct(a, "quark", a.nb_toy_quark_k2.lo, delta_abs=a.eps_stat.lo)         # 4.79e3
    def _grid_rw(k2, pairs):
        return sum(row_samples(a, "baryon", k2, T, a.eps_chi4_chi2_target.lo, pr, with_nb=False) for T, pr in pairs)
    g_rw_lo = _grid_rw(a.nb_toy_k2_priced.lo, a.grid_rw_pairs_k2lo.value)                  # 2.27e5
    g_rw_mid = _grid_rw(a.nb_toy_k2_mid.lo, a.grid_rw_pairs_k2mid.value)                   # 1.29e6
    # first result: the T = 150 MeV row at 30%, each toy at min(direct, two reweighted ensembles)
    fr_by_toy = {}
    for kind, k2, pr in a.fr_rw_pairs.value:
        direct = row_samples(a, kind, k2, 150.0, a.eps_first_result.lo, None)
        rw = row_samples(a, kind, k2, 150.0, a.eps_first_result.lo, pr)
        fr_by_toy[(kind, k2)] = dict(direct=direct, reweighted=rw, best=min(direct, rw))
    fr_exact = (min(v["best"] for v in fr_by_toy.values()), max(v["best"] for v in fr_by_toy.values()))
    fr_n = tuple(float(math.ceil(x)) for x in fr_exact)            # 2,926, 15,861 whole samples

    # R9 depth: the factory.json band was for the 2e3-step circuit; the QSVT queries repeat the same step, so each end's
    # T-depth is that band scaled by its own T ratio (F* = 17.5-56 carried at either end, Assumed). The exports pair
    # the bands crosswise (common.depth_exports), which widens F* to 8.6-114; the matched band is f_star_matched.
    sc = tuple(t / t_hsim_bound for t in t_range)
    depth33 = (a.t_depth_2033.value[0] * sc[0], a.t_depth_2033.value[1] * sc[1])
    dx = depth_exports(t_range, depth33, shots, a.t_gate_s, a.shot_overhead_s)
    f_matched = (t_hsim_bound / a.t_depth_2033.value[1], t_hsim_bound / a.t_depth_2033.value[0])   # 17.5, 56.1
    # The contract exports (f_star, baseline_ok, fits_1yr) stay cross-paired (CONTRACT rule 12, test_common checks
    # f_star = t / d crosswise). Their (T_lo, D_hi) corner is not a circuit, so the chapter quotes the matched-end
    # values below: F* 17.5-56 at either end, so the 10-factory baseline holds at both (verifier 2026-10-05).
    f1yr = dx["factories_for_1yr"]
    baseline_ok_matched = (f_matched[0] >= FACTORY_BASELINE, f_matched[1] >= FACTORY_BASELINE)
    fits_1yr_matched = (f1yr[1] <= f_matched[0], f1yr[0] <= f_matched[1])
    d_end = {"low": (a.t_depth_2033.value[0] * sc[0], a.t_depth_2033.value[1] * sc[0]),
             "high": (a.t_depth_2033.value[0] * sc[1], a.t_depth_2033.value[1] * sc[1])}
    wall_shot_lo = _wall_shot_at_baseline(a, t_range[0], d_end["low"][1])     # serial: F* >= 17.5 at each end
    wall_shot_hi = _wall_shot_at_baseline(a, t_range[1], d_end["high"][1])
    wall_shot_rng = (wall_shot_lo, wall_shot_hi)                   # 2.8e3 - 5.6e3 s
    wall_campaign = (shots * wall_shot_lo, shots * wall_shot_hi)   # 23 - 46 yr
    wall_fr = (fr_n[0] * wall_shot_lo, fr_n[1] * wall_shot_hi)     # 0.26 - 2.8 yr
    wall_point = (shots_point * wall_shot_lo, shots_point * wall_shot_hi)
    floor_campaign = (shots * d_end["low"][0] * REACTION_TIME_S, shots * d_end["high"][1] * REACTION_TIME_S)
    horizon_s = a.campaign_horizon_yr.lo * SECONDS_PER_YEAR
    # r23: one machine. What the 5-year horizon holds, and the cut the full grid needs to fit it.
    shots_in_horizon = (horizon_s / wall_shot_hi, horizon_s / wall_shot_lo)
    points_in_horizon = tuple(x / shots_point for x in shots_in_horizon)
    horizon_reduction = (wall_campaign[0] / horizon_s, wall_campaign[1] / horizon_s)   # 4.5 - 9.2x
    t_shot_to_fit = (horizon_s / shots - a.shot_overhead_s.lo) / a.t_gate_s.lo   # 6.1e8 T/shot at 2.58e5 shots
    fr_reduction = (wall_fr[0] / horizon_s, wall_fr[1] / horizon_s)              # < 1: the first result fits

    eps_l_req = (a.faults_per_shot.lo / t_range[1], a.faults_per_shot.lo / t_range[0])   # 1.8e-11, 3.6e-11
    faults_floor = (t_range[0] * a.eps_l_floor.lo, t_range[1] * a.eps_l_floor.lo)        # 28 - 57 at 1e-8

    # pre-ruling QSVT comparison (record): R x d_QSP x 2 queries x (2-4) x c_step, R = overlap^(-1/2), d at eps 1e-3
    d_qsp = math.sqrt(beta_h * math.log(1.0 / a.eps_qsp.lo))       # 26.3
    d_naive = beta_h                                               # ~100: the linear degree the QSVT saves on
    r_unstructured_log2 = lq_system / 2.0                          # sqrt(2^264) -> 2^132
    r_old = 1.0 / math.sqrt(a.filter_success_prob.lo)              # 10
    def _rail(d_poly):
        steq = tuple(r_old * d_poly * 2 * f for f in a.qsvt_query_cost.value)
        nr = tuple(s_ * rot["total"] for s_ in steq)
        ee = tuple(_eps(a, n_) for n_ in nr)
        tt = tuple(s_ * _t_per_step(a, d, L, ns33, e_) for s_, e_ in zip(steq, ee))
        return steq, nr, ee, tt
    old_steq, old_nrot, old_eps, old_t = _rail(d_qsp)

    # breakdown: the low end of the QSVT shot (query = 2 c_step), at that circuit's own eps
    n_tr, e_s = q_lo["step_equivalents"], q_lo["eps_rot"]
    t_rot_s = q_lo["t_per_rotation"]
    prims = _gauge_primitives(a, d, geo["n_links"], n_tr, e_s) + _hop_primitives(
        a, ns33, geo["n_links"], geo["n_sites"], n_tr, e_s) + (
        Primitive("qsvt_reference_rotations", q_lo["n_rot_ref"], t_rot_s, CircuitStatus.SCALING,
                  "author ruling 2026-10-05; app06 route (iv)",
                  "product reference, one rotation per system qubit, once per filter call (17 calls)"),
        Primitive("qsvt_signal_rotations", q_lo["n_rot_signal"], t_rot_s, CircuitStatus.SCALING,
                  "gilyen2018QSingValTransfArXiv; author ruling 2026-10-05",
                  "QSVT phases (d_beta + 1 per filter call) and one fixed-point phase per amplification round"),
        Primitive("qsvt_reflection_toffolis", q_lo["aa_rounds"] * q_lo["reflection_toffolis"],
                  toffoli_t(1, a.toffoli_convention.value), CircuitStatus.SCALING, "author ruling 2026-10-05",
                  "one reflection about the start state per round, on the full register (high end, 1014 Toffolis, MBU)"))

    # "V=3^3" (open questions): the same filter on V=3^3 (beta||H|| extensive, 337.5), same amplitude assumed
    geo3 = _geometry(d, 3, a.n_stag_2033.value, a.nc.value, a.group.value)
    v3_lq = (geo3["gauge_lq"] + geo3["fermion_lq"] + anc_lo, geo3["gauge_lq"] + geo3["fermion_lq"] + anc_hi)
    v3_beta_h = beta_h * geo3["n_sites"] / geo["n_sites"]
    v3_qb = qsvt_band(a, d, 3, ns33, v3_beta_h, geo3["gauge_lq"] + geo3["fermion_lq"], v3_lq[1])
    v3_t = (v3_qb["low"]["t_shot"], v3_qb["high"]["t_shot"])      # 1.8e10 - 3.7e10
    # records (pre-ruling): uncontrolled 2e3-step / 10x-cut Hamiltonian simulation and the 10x-cut Gibbs sampler
    v3_steps = (n_steps * geo3["n_sites"] / geo["n_sites"]) / 10.0  # 675
    v3_n_rot = v3_steps * _rot_per_step(a, d, 3, ns33)["total"]
    v3_eps = _eps(a, v3_n_rot)
    v3_c_be = _t_per_step(a, d, 3, ns33, v3_eps)
    v3_t_hsim = v3_steps * v3_c_be
    v3_gb = gibbs_band(a, d, 3, ns33, v3_beta_h, mix_cut=10.0)
    v3_g_shot = v3_gb["central"]["t_sample"]                       # 3.1e14

    inter = dict(
        n_sites=geo["n_sites"], n_links=geo["n_links"], n_plaquettes=geo["n_plaquettes"],
        link_qubits=geo["link_qubits"], link_qubits_chapter_old=geo["link_qubits_chapter_old"],
        link_qubits_log2=geo["link_qubits_log2"],
        gauge_lq=geo["gauge_lq"], fermion_lq=geo["fermion_lq"], ancilla_lq=(anc_lo, anc_hi),
        lq_system=lq_system, lq_total=lq, hops=geo["hops"],
        **pl,
        c_be_gauge_t_per_step=gauge_step,
        eps_syn=a.eps_syn.lo, n_rot_per_step=rot["total"], n_rot_gauge_per_step=rot["gauge"],
        n_rot_hop_per_step=rot["hop"], n_rot_per_shot=n_rot_shot, eps_rot=e,
        t_per_rotation=t_per_rotation(e, a.rot_synthesis_paper.value),
        synthesis_error_per_shot=n_rot_shot * e ** 2,
        c_be_at_fixed_1em4_superseded=c_be_fixed_1e4,
        t_per_shot_at_fixed_1em4_superseded=n_steps * c_be_fixed_1e4,
        c_be_gauge_t_per_step_report_convention=_gauge_t_per_step_report_convention(a, d, L, e),
        **_hop_record(a, ns33, e, geo["n_links"]),
        hop_t_per_step=hop_step, gauge_share_of_c_be=gauge_step / c_be, hop_share_of_c_be=hop_step / c_be,
        hop_t_per_step_rhop_legacy=geo["n_links"] * hl["legacy"]["rhop_rus_slope"],
        n_rot_mass_per_step=rot["mass"], mass_hwp_k=ms["k"], mass_t_per_site=ms["t"], mass_t_per_step=mass_step,
        c_be_share_draft=gauge_step + geo["n_links"] * _hop_link(a, ns33, e, share="draft")["t"] + mass_step,
        c_be_stated_superseded=c_be_old, c_be_over_superseded=c_be / c_be_old,
        hop_allowance_t_per_step_superseded=c_be_old - _gauge_t_per_step(a, d, L, 1e-4, papers=True),   # record, old 1e-4
        c_be_gauge_t_per_step_papers_convention=_gauge_t_per_step(a, d, L, e, papers=True),   # legacy record (R2)
        t_reference_2033=a.t_ref_2033.lo, reference_fraction=(t_range[0] / a.t_ref_2033.lo, t_range[1] / a.t_ref_2033.lo),
        hard_ops_is_lower_bound=False, hard_ops_is_conditional=True,
        t_hsim_lower_bound=t_hsim_bound, hsim_bound_is_lower_bound=True,
        beta_h=beta_h, c_mix=a.c_mix.lo, n_gibbs_steps=n_steps, c_be_t_per_step=c_be,
        # ---- headline: QSVT/TPQ thermal state (author ruling 2026-10-05) ----
        thermal_route="QSVT/TPQ", qsvt=qb, t_per_shot_range=t_range, t_per_shot_value=t_shot,
        qsvt_success_prob=q_lo["success_prob"], qsvt_amplitude=q_lo["amplitude"], qsvt_aa_rounds=q_lo["aa_rounds"],
        qsvt_filter_calls=q_lo["filter_calls"], d_beta=q_lo["d_beta"], d_beta_exact=q_lo["d_beta_exact"],
        delta_beta=a.delta_beta.lo, qsvt_queries=q_lo["queries"],
        qsvt_step_equivalents=(q_lo["step_equivalents"], q_hi["step_equivalents"]),
        qsvt_n_rot_per_shot=(q_lo["n_rot"], q_hi["n_rot"]), qsvt_eps_rot=(q_lo["eps_rot"], q_hi["eps_rot"]),
        qsvt_c_step=(q_lo["c_step"], q_hi["c_step"]),
        qsvt_t_per_rotation=(q_lo["t_per_rotation"], q_hi["t_per_rotation"]),
        qsvt_t_per_shot=t_range, qsvt_overhead_share=(q_lo["overhead_share"], q_hi["overhead_share"]),
        qsvt_t_reflections=(q_lo["t_reflections"], q_hi["t_reflections"]), qsvt_t_readout=0.0,
        # pre-ruling QSVT comparison (records)
        d_qsp_old=d_qsp, d_qsp_stated=a.d_qsp_stated.lo, d_naive_linear=d_naive, R_unstructured_log2=r_unstructured_log2,
        R_warm_old=r_old, R_warm_stated=a.R_warm_stated.lo, qsvt_old_step_equivalents=old_steq,
        qsvt_old_t_per_shot=old_t, qsvt_t_per_shot_stated=a.qsvt_rail_stated.value,
        # ---- Gibbs sampler, RECORD since the ruling (referee F2; one sentence in the chapter) ----
        gibbs_sampler=gb, gibbs_n_jumps=gc["n_jumps"], gibbs_jump_families=gc["jump_families"],
        gibbs_sampler_steps=gc["sampler_steps"], gibbs_trotter_per_sampler_step=gc["trotter_per_sampler_step"],
        gibbs_c_step_controlled=gc["c_step_controlled"], gibbs_c_step=gc["c_step"], gibbs_eps_rot=gc["eps_rot"],
        gibbs_t_per_sampler_step=gc["t_per_sampler_step"], gibbs_t_sample_band=g_band,
        gibbs_t_per_shot=g_shot, gibbs_xk_sample=gxk, gibbs_t_per_shot_xk=g_xk, gibbs_t_per_shot_range=g_range,
        gibbs_xk_over_central=g_xk / g_shot, gibbs_over_hsim_bound=g_shot / t_hsim_bound,
        gibbs_over_qsvt=(g_shot / t_range[1], g_xk / t_range[0]),
        gibbs_wall_per_shot_yr=tuple((t * a.t_gate_s.lo + a.shot_overhead_s.lo) / SECONDS_PER_YEAR for t in g_range),
        gibbs_wall_campaign_yr=tuple(shots * (t * a.t_gate_s.lo + a.shot_overhead_s.lo) / SECONDS_PER_YEAR for t in g_range),
        gibbs_eps_l_required=a.faults_per_shot.lo / g_shot,
        # chain reuse (Gibbs-only lever, record): tau = one tenth of the sample at the sampler's cost
        chain_reuse_gain=_rng(a.chain_reuse_gain), chain_reuse_steps_per_tau=a.chain_reuse_steps_per_tau.lo,
        chain_reuse_eps_l_required=a.faults_per_shot.lo * (n_steps / a.chain_reuse_steps_per_tau.lo) / g_shot,
        chain_reuse_eps_l_required_hsim_superseded=a.faults_per_shot.lo / (a.chain_reuse_steps_per_tau.lo * c_be),
        # ---- walls, one machine (r23), at each end of the per-shot range ----
        t_depth_per_shot_scaled=depth33, depth_scale=sc, f_star_matched=f_matched, t_depth_by_end=d_end,
        baseline_ok_matched=baseline_ok_matched, fits_1yr_matched=fits_1yr_matched,
        wall_per_shot_s=wall_shot_rng, wall_per_shot_min=tuple(w / 60.0 for w in wall_shot_rng),
        shot_overhead_s=a.shot_overhead_s.lo,
        shots_per_point=shots_point, n_grid_points=n_grid, shots_total=shots,
        t_campaign_total=(shots * t_range[0], shots * t_range[1]),
        eps_stat_chi4=a.eps_stat.lo, target_accuracy_pressure=a.eps_pressure_target.lo,
        target_accuracy_chi4_over_chi2=a.eps_chi4_chi2_target.lo,
        wall_campaign_yr=tuple(w / SECONDS_PER_YEAR for w in wall_campaign),
        wall_single_point_yr=tuple(w / SECONDS_PER_YEAR for w in wall_point),
        floor_wall_campaign_yr=tuple(w / SECONDS_PER_YEAR for w in floor_campaign),
        campaign_horizon_yr=a.campaign_horizon_yr.lo,
        shots_in_horizon=shots_in_horizon, grid_points_in_horizon=points_in_horizon,
        horizon_reduction_needed=horizon_reduction, t_per_shot_to_fit_horizon=t_shot_to_fit,
        fr_reduction_needed=fr_reduction, fr_fits_horizon=fr_reduction[1] <= 1.0,
        faults_per_shot_budget=a.faults_per_shot.lo, eps_l_required=eps_l_req, faults_at_rfi_floor=faults_floor,
        rfi_over_eps_l=(a.eps_l_floor.lo / eps_l_req[1], a.eps_l_floor.lo / eps_l_req[0]),
        sigma360_link_qubits_min=math.ceil(math.log2(1080)),
        compression_36_vs_72=GROUPS[a.group_stretch.value].link_qubits / GROUPS[a.group.value].link_qubits,
        v3cubed_lq=v3_lq, v3cubed_beta_h=v3_beta_h, v3cubed_qsvt=v3_qb, v3cubed_d_beta=v3_qb["low"]["d_beta"],
        v3cubed_c_step=v3_qb["low"]["c_step"], v3cubed_t_per_shot=v3_t,
        v3cubed_steps=v3_steps, v3cubed_n_rot_per_shot=v3_n_rot, v3cubed_eps_rot=v3_eps,
        v3cubed_c_be=v3_c_be, v3cubed_t_hsim_with_10x_cut_superseded=v3_t_hsim,
        v3cubed_n_jumps=v3_gb["central"]["n_jumps"], v3cubed_gibbs=v3_gb, v3cubed_gibbs_t_with_10x_cut=v3_g_shot,
        # r25 R9 exports (names fixed by CONTRACT rule 12): campaign = the 5x5 grid
        t_per_shot=dx["t_per_shot"], t_depth_per_shot=dx["t_depth_per_shot"], f_star=dx["f_star"],
        floor_wall_s=dx["floor_wall_s"], factories_for_1yr=dx["factories_for_1yr"],
        wall_first_result_s=wall_fr, wall_campaign_s=wall_campaign,
        wall_serial_s=dx["wall_serial_s"], baseline_ok=dx["baseline_ok"], fits_1yr=dx["fits_1yr"],
        floor_wall_yr=tuple(x / SECONDS_PER_YEAR for x in dx["floor_wall_s"]),
        # first result (R1; re-derived for referee F1): one T row at 30%, min(direct, two reweighted ensembles)
        eps_first_result=a.eps_first_result.lo, fr_samples=fr_n, fr_samples_exact=fr_exact, fr_by_toy=fr_by_toy,
        wall_first_result_d=tuple(x / 86400.0 for x in wall_fr), wall_first_result_yr=tuple(x / SECONDS_PER_YEAR for x in wall_fr),
        # referee F1: shots from the delta-method variance of k4/k2 (toy N_B), priced at kappa_2(0) = 0.01
        delta_chi4_chi2_abs=a.eps_stat.lo, nb_toy_k2_priced=a.nb_toy_k2_priced.lo,
        shots_per_point_by_mu_T=g_pr["per_point"], shots_per_point_min=g_pr["min"], shots_per_point_max=g_pr["max"],
        grid_shots_k2_mid=g_mid["total"], grid_shots_k2_hi=g_hi["total"],
        grid_wall_k2_mid_yr=tuple(g_mid["total"] * w / SECONDS_PER_YEAR for w in wall_shot_rng),
        grid_wall_k2_hi_yr=tuple(g_hi["total"] * w / SECONDS_PER_YEAR for w in wall_shot_rng),
        grid_shots_quark_rel=g_q_rel["total"], grid_shots_quark_abs=g_q_abs["total"],
        grid_shots_reweighted_k2_lo=g_rw_lo, grid_shots_reweighted_k2_mid=g_rw_mid,
        reweighting_saving_k2_lo=1.0 - g_rw_lo / shots, reweighting_ratio_k2_mid=g_rw_mid / g_mid["total"],
    )
    return Result(
        era="2033", lq=lq, hard_ops=t_range, breakdown=prims, intermediates=inter,
        shots=(shots, shots), wall_time_s=wall_campaign,
        epsilon_l=eps_l_req,
        notes=(
            "LQ = 192 gauge (24 links x 8, R1) + 72 fermion (3 x 3 x 8) + 250-750 ancilla (BE workspace and the parity "
            "ancilla of parallel hops) = 514-1014, box '~500-1000'.",
            "T/shot (author ruling 2026-10-05): the QSVT/TPQ thermal state, Ch. 9's construction. Filter "
            "e^{-beta(H - mu_B N_B)/2} of degree d_beta = round(sqrt(1e2 ln 1e4)) = 30; amplitude 0.1 on the reference "
            "(Assumed, as Ch. 9) -> k = ceil(pi/(4 arcsin 0.1)) = 8 rounds, 2k + 1 = 17 filter calls, 510 queries; a query "
            "is 2-4 c_step (Stated) -> 1,020-2,040 step equivalents at their own eps (1.1e-5, 8.0e-6): 2.8e9-5.6e9 T. "
            "Reference (one rotation per system qubit), QSVT and amplification phases, 8 reflections (full register, "
            "<= 1014 Toffolis each): "
            "1e-4 of the shot. Readout of N_B: computational basis, no T. Conditional on the reference (referee C1).",
            "c_BE derived here: gauge 24 x 19,230 = 461,524 (17%) + hop 24 x 95,734 = 2,297,627 (83%) + mass 8 x "
            "HWP(9) 163.6 = 1,309 = 2,760,459 T/step at the 2e3-step tolerance; replaces the Stated 5e5.",
            "Gibbs sampler (referee F2), record: 8.39e13-4.42e14 T per sample (7.6e12-9.5e14 band), 1.5e4-1.6e5x the "
            "QSVT shot; one sentence in the chapter ('~8e13-4e14 T per shot').",
            "Shots (referee F1) unchanged: 257,603 over the grid at kappa_2(0) = 0.01; 9.08e5 at 0.05; 2.24e6 at 0.1.",
            "Walls, one machine: 2.8e3-5.6e3 s per shot (46-94 min); first result (one T row, 2,926-15,861 samples) "
            "0.26-2.8 yr; campaign (5x5 grid, 2.58e5 shots) 23-46 yr, a 4.5-9.2x cut to fit 5 years; depth limit 4.0-26 yr.",
            "epsilon_l (R3): 0.1 expected faults per shot -> 1.8e-11 (upper end) to 3.6e-11, 280-560x below the RFI rate.",
            "R9: the factory.json T-depth (2e3-step circuit) scaled by each end's T ratio, F* 17.5-56 at either end "
            "(f_star_matched); the cross-paired export is 8.6-114.",
            "V=3^3: beta||H|| 337.5, d_beta 56, c_step 9.6e6: 1.8-3.7e10 T at the same amplitude (records: 10x-cut Gibbs "
            "3.1e14; 675-step Hamiltonian simulation 6.3e9).",
        ),
    )


def _model_codesign(a: Assumptions) -> Result:
    """Post-2033 aspirational stretch (app07:196): Sigma(72x3) on V=4^3. Register only; T not costed."""
    geo = _geometry(a.d_stretch.value, a.L_stretch.value, a.n_stag_stretch.value, a.nc.value,
                    a.group_stretch.value)
    anc = _rng(a.n_anc_stretch)
    lq_system = geo["gauge_lq"] + geo["fermion_lq"]
    lq = (lq_system + anc[0], lq_system + anc[1])                # 1728 + 576 + 250..750 = 2554..3054
    inter = dict(
        n_sites=geo["n_sites"], n_links=geo["n_links"], n_plaquettes=geo["n_plaquettes"],
        link_qubits=geo["link_qubits"], link_qubits_chapter_old=geo["link_qubits_chapter_old"],
        link_qubits_log2=geo["link_qubits_log2"],
        gauge_lq=geo["gauge_lq"], fermion_lq=geo["fermion_lq"], ancilla_lq=anc, lq_total=lq,
        ratio_to_2033_anchor=(lq[0] / 1000.0, lq[1] / 1000.0),
        hard_ops_costed=False,
    )
    return Result(
        era="codesign", lq=lq, hard_ops=(0.0, 0.0), breakdown=(), intermediates=inter,
        notes=(
            "LQ = 1728 gauge (192 links x 9, compiled Sigma(72x3) register, R1) + 576 fermion (3 x 3 x 64) + 250-750 "
            "BE ancilla of the 2033 route = 2554-3054, box '~2600-3100'; 'exceeds 2033 register anchor by "
            "~2.6-3.1x' (2.554-3.054).",
            "No hard-op count is stated for the stretch (resources.py carries it as LQ-only); hard_ops=(0,0) here. "
            "No Sigma(72x3) FFT exists (arXiv:2511.17437 gives a lower bound only).",
            "Ancilla carried V-independently from the 2033 Gibbs route (ruling ch10-stretch-ancilla a); the "
            "earlier ~50 (2354 LQ) is retired.",
        ),
    )


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
    "2028": Published(lq=(180, 180), hard_ops=(4.0e5, 4.0e5),
                      src="app07:2028 box ('~180' rounds 176; 'Per-shot T >~4.0e5': one Trotter step, a lower bound, 4.0x the cap; "
                          "gauge 1.0e5 (E20 full fit) + hops 3.0e5 (share='link' + undo, E21) + mass 1.9e2 = 399,957 "
                          "at the R-TOL eps_rot 9.6e-4, r19); rel_tol 0.10",
                      rel_tol=0.10),
    "2033": Published(lq=(500, 1000), hard_ops=(2.8e9, 5.6e9),
                      src="app07:2033 box ('~500-1000' rounds 514-1014; '2.8e9-5.6e9': the QSVT/TPQ thermal state, author "
                          "ruling 2026-10-05, Ch. 9's construction: 17 filter calls of degree 30, 2-4 c_step per query, "
                          "2.7726e9-5.6340e9; conditional on the reference state; was '8.4e13-4.4e14', the priced "
                          "Chen-Kastoryano-Gilyen sampler, open item 1(a) 2026-10-04; before F2 '>~5.5e9'); rel_tol 0.10",
                      rel_tol=0.10),
    "codesign": Published(lq=(2554, 3054), hard_ops=None,
                          src="app07:2033 box, post-2033 stretch ('~2600-3100 LQ (1728 + 576 + 250-750)'); T not "
                              "costed; exact register carried (ruling ch10-stretch-ancilla a), box rounds once",
                          rel_tol=0.10),
}


def INSTANCE_ROWS(a, era, r):
    """Rows for resources.json: [(label, (lq_lo, lq_hi), (t_lo, t_hi), {extra}), ...]."""
    if era == "2028":
        return [(r"2+1D SU(3) $\mu_B{>}0$", r.lq, r.hard_ops,
                 {"note": "one Trotter step (lower bound on the executed ramp), 4.0x the cap; gauge terms compiled (full fit), hops from the "
                          "unpublished Fermion_Primitives counts with the colour squish held per link and an "
                          "estimated frame undo (r19)"})]
    if era == "2033":
        return [(r"3D SU(3) $\Sigma(36{\times}3)$ $2^3$, QSVT/TPQ", r.lq, r.hard_ops,
                 {"conditional": True,
                  "note": "QSVT/TPQ thermal state (author ruling 2026-10-05): 17 filter calls of degree 30, 2-4 Trotter "
                          "steps per query; conditional on a reference state of filtered amplitude ~0.1 (referee C1). "
                          "The Gibbs sampler (8.4e13-4.4e14) is a record, not plotted"})]
    if era == "codesign":
        return [(r"post-2033 $\Sigma(72{\times}3)$ $4^3$", r.lq, r.hard_ops,
                 {"note": "T not costed", "lq_only": True})]
    return []


if __name__ == "__main__":
    a = Assumptions()
    for era in ERAS:
        r = model(a, era)
        print(f"{era}: LQ {r.lq}  T {r.hard_ops}  breakdown {r.breakdown_total():.3g}")
        for k, v in r.intermediates.items():
            print(f"    {k:44s} {v}")
