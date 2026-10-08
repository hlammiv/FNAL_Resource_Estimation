"""Ch. 9 (app06_chiral_gauge.tex): every named intermediate the prose states, so that a
future edit to one number without the others fails loudly.

Known gaps between the chapter's stated inputs and its printed numbers are marked
xfail(strict=True) and filed in NEEDS_AUTHOR.md: they fail loudly in the *other*
direction (XPASS) once the authors fix the chapter, so the marker gets removed.

2026-09-28 round B (rulings R5, VERTEX, R6; H. Lamm): Toffolis are 7 T everywhere, U_x = 392 T
read from groups.py, the 2028 ramp is N_Trotter = 20 with the grouping stated in prose ('color+s').

2026-09-29 R-HOP and R-TOL (H. Lamm): the Z3 spatial hop is 8 strings per copy, C_W = 0, k = 12;
each circuit sets eps_rot = sqrt(1e-2 / N_rot) from its own rotation count.

2026-10-01 r17 (rulings e, f; apply_log/r17_ch09.md): 2028 hop rotations at the full fit and the
on-site domain-wall term m + r(1+d) as one k=192 phasing group: 1.95-2.05e5 -> 2.76-2.97e5 T/shot.
2033 link terms priced by groups.hop_link_cost (unpublished Fermion_Primitives gate counts, frame
undo and n_spin = 2 estimated, MBU): c_be 5.5e4 -> 2.56e6, 1.48e10 -> 6.89e11 T/shot, band
3.383e9-1.2952e12 (top = the draft's W2).

2026-10-01 r19 (E21 (1)-(3); apply_log/r19_ch09.md): the colour squish and parity are computed once per link and
held through V, hop, V^dag (share="link"), undo kept (our estimate, drawn but not counted in the draft); n_spin = 2
confirmed. Per link 1182 -> 384 Toffolis, rotations unchanged (704): c_be 2.56e6 -> 2.10e6, centre 6.891e11 ->
5.672e11 T/shot. 2028, the band floor and the W2 top do not move. The r17 centre is the record share="draft".

2026-10-01 r20 (E23 (1)-(5); apply_log/r20_ch09.md): (1) 2O = single-qubit Clifford group mod phases, stated as a
fact; (2) the band floor is the derived 2O colour multiplexer (3 Toffolis + 18 T, verified below on all 48 elements),
468 T per link: floor c_be 12,551.6 -> 61,043.6, band 3.383e9 -> 1.645e10 at the bottom; (3) the Gibbs sampler is
estimated here (4.4e14 / 2.3e15 / 1.4e16 T); (4) on-site basis kept; (5) ancilla itemized at peak: 2028 39 (247 LQ),
2033 43 (1015 LQ; verifier: parity held under share='link').

2026-10-01 r21 (H. Lamm: "1) promote 2) do TPQ 3) settle them"; apply_log/r21_ch09.md): (1) the derived multiplexer is
the 2033 headline link circuit (468 T per link per query), the draft's diagonalizing frame the stated alternative, one
headline and no band; the link-term workspace becomes the multiplexer's (33 ancilla, 1005 LQ). (2) TPQ is the thermal
route; the Gibbs sampler is re-priced at the promoted c_be. (3) five derived items, each tested below: (a) 24
fifth-direction hops (open boundary); (b) N_Trotter = 83,195 from the second-order commutator bound; (c) 72 LCU terms
per site; (d) 1125 shots at 2028; (e) the 1+1D Z3 block encoding, c_be 2.2e3 T. 2028 2.69-2.89e5 T/shot; 2033
1.809e13 T/shot.

2026-10-01 r22 (E26; apply_log/r22_core.md, r22_ch09.md): (i) the 2033 sign function in the single-particle architecture
of arxiv_2607_28524, compiled here: 598 Toffolis + 24 T + 6 rotations per query, 692 queries + lift + reflection + the
Higgs vertex once per application (3.06e6 T), 38 ancilla, 1010 LQ; the second-quantized construction of r17-r21 is a
record. (ii) N_Trotter from the state-dependent estimate (Alves, Lamm, Liu, in preparation): 5031 steps (1186-14,242),
13 applications per step (the Jacobi-Anger order), 2.008e11 T/shot; the commutator bound (83,195 steps, 1.276e12) is
the quoted worst case. (iii) 2028 ancilla 31 + 1 = 32, 240 LQ (E26: k - w(k) per phasing group, one RUS ancilla).
The 1+1D overlap comparison is repriced in the same architecture (5.6e4 T per application, 5.8e5 for ten).

2026-10-02 r24 (H. Lamm: "ok to run for a shorter amount of time, and lower kappa is ok"; apply_log/r24_ch09.md):
kappa 1e2 -> 30 (d_sgn 208), t = 50a -> 10a (450 steps, 19 applications per step). The sphaleron shot is 8.2414e9 T;
the condensate arm is the TPQ preparation alone, 2.2412e8 T. Headline constants below are r24; the r22 headline and
the 50a worst case are kept as records (N_TROT_R22, N_BOUND).

2026-10-02 r25 (rulings R1, R4, R9, R10 and the queued Ch. 9 rulings; apply_log/r25_ch09.md): (a) the Ginsparg-Wilson-
projected Yukawa is inside: a second degree-208 sign sequence per application, 417 queries, 1,847,747.1 T per
application, 1.6242e10 T per 10a shot, 4.4166e8 per condensate shot; (R4) the 5a shot at its own depth, 7.2366e9 T;
(R10) the 2033 fan-out (+5 LQ, 1015) and the 2028 schedule of three spatial-hop groups per slot (42 ancilla, 250 LQ),
F* >= 10; (R1, e) first result 1.591 yr, campaign 4-8 temperatures 29.8-50.5 yr. Headline constants below are r25;
the r24 numbers are reproduced by the fixture a24 (Yukawa outside, as-itemized register).
"""

import dataclasses
import inspect
import math
import re
from pathlib import Path

import pytest

from estimates import common as c
from estimates import ch09_chiral_gauge as m
from estimates.groups import GROUPS, rhop_link, hop_link_cost, fermion_hop_counts

NEEDS = "NEEDS_AUTHOR.md"


@pytest.fixture(scope="module")
def a():
    return m.Assumptions()


@pytest.fixture(scope="module")
def r28(a):
    return m.model(a, "2028")


@pytest.fixture(scope="module")
def r33(a):
    return m.model(a, "2033")


@pytest.fixture(scope="module")
def a24(a):
    """The r24 chain: Yukawa outside sgn, the 2033 register as itemized (1010 LQ). Reproduces every r24 print."""
    return dataclasses.replace(a, yukawa_placement=c.Stated("outside", "r24 record"),
                               depth_variant_2033=c.Assumed("register", "r24 record"),
                               beta_H_rule=c.Assumed("stated", "r24 record: the filter at beta||H|| = 1e2"),
                               aa_call_rule=c.Assumed("per_round", "r24 record: one filter call per round"))


@pytest.fixture(scope="module")
def r33_24(a24):
    return m.model(a24, "2033")


def rel(x, y):
    return abs(x - y) / abs(y)


# R-TOL by hand (not through the model): one tolerance per circuit from its own count
def fit(n_rot):
    """RUS fit at the R-TOL tolerance of a circuit with n_rot rotations."""
    return 1.15 * math.log2(1 / math.sqrt(1e-2 / n_rot)) + 9.2


N28, N28C = 290 * 21, 336 * 21                        # 6090 unchunked, 7056 chunked (r21 a: 24 fifth-direction hops)
EPS28, EPS28C = math.sqrt(1e-2 / N28), math.sqrt(1e-2 / N28C)
T_RUS, T_RUS_C = fit(N28), fit(N28C)                  # 20.249, 20.371
# 2033, r22: the single-particle architecture at the state-dependent step count
E_HALF, B_NORM = 0.5 * 0.5 * 3.75, 4 + 4 * 4.0        # g^2 = 1: e = 15/16, b = hop 4 + 4 plaquettes x 4
LAMBDA = 81 * (E_HALF ** 2 * B_NORM / 3 + E_HALF * B_NORM ** 2 / 6)     # 5537.1 a^-3 (the worst-case bound)
N_BOUND = math.ceil(math.sqrt(LAMBDA * 50 ** 3 / 0.1))  # 83,195 at t = 50a (r21-r23; the records keep it)
N_BOUND_10 = math.ceil(math.sqrt(LAMBDA * 10 ** 3 / 0.1))   # 7442 at t = 10a (r24)
# r24 (H. Lamm 2026-10-02, 'ok to run for a shorter amount of time, and lower kappa is ok'): kappa = 30, t = 2 fm = 10a
D_SGN = math.ceil(30 * math.log(1e3))                 # 208
N_TROT, PER_STEP = 450, 19                            # state-dependent estimate at t = 10a; Jacobi-Anger order at 2Q dt = 9.6
# referee C1 (2026-10-04): the thermal filter at beta||H|| = 2Q / (T a) = 864: d_beta = round(sqrt(864 ln 1e4)) = 89
D_BETA = round(math.sqrt(2 * 216 / 0.5 * math.log(1e4)))      # 89 (r25: 30 at the record beta||H|| = 1e2)
# verifier (2026-10-04): 8 amplification rounds call the filter 2 x 8 + 1 = 17 times (first referee pass 8: 712)
TPQ, TPQ_R25, TPQ_C1 = D_BETA * 17, 240, D_BETA * 8     # 1513 applications (C1 712; r25 240)
NDOV = PER_STEP * N_TROT + TPQ                        # 10,063 applications (C1 9262; r25 8790)
NQD = 2 * D_SGN + 1                                   # 417 queries per application (r25 (a); r24 209)
NQD_24 = D_SGN + 1                                    # 209 (r24, Yukawa outside)
NQ = NDOV * NQD                                       # 3,665,430 single-particle queries (r24 1,837,110)
NROT_APP = NQD * 6 + 6 + 1                            # 2509 rotations per application (r24 1261)
C_HIGGS = 27 * 392 + 27 * 40 * 7                      # 18,144: U_x per site + 40 SELECT strings per site
# r22-r23 headline (kappa 1e2, t = 50a), kept as a record
N_TROT_R22, PER_STEP_R22, NDOV_R22 = 5031, 13, 65643


def price_sp(n_dov, d=D_SGN, n_sgn=2):
    """One shot of n_dov applications by hand: (d + 1) + (n_sgn - 1) d queries of 598 Toffolis + 24 T + 6 rotations
    (r25: n_sgn = 2, the Yukawa projector), the lift (862 Toffolis + 12 T + 6 rotations), the reflection (18 Toffolis +
    1 rotation), the Higgs vertex; R-TOL from its own count."""
    nq = (d + 1) + (n_sgn - 1) * d
    n_rot = n_dov * (nq * 6 + 6 + 1)
    tr = fit(n_rot)
    c_sp = 598 * 7 + 24 + 6 * tr
    t_dov = nq * c_sp + (862 * 7 + 12 + 6 * tr) + (18 * 7 + tr) + C_HIGGS
    return {"n_rot": n_rot, "t_rot": tr, "c_sp": c_sp, "t_dov": t_dov, "t_shot": n_dov * t_dov}


N33 = NDOV * NROT_APP                                 # 25,248,067 rotations per shot (C1 23,238,358; r25 22,054,110)
EPS33 = math.sqrt(1e-2 / N33)
T33 = fit(N33)                                        # 27.047 (r24 26.476; r22 29.134)
C_SP = 598 * 7 + 24 + 6 * T33                         # 4372.3 (r24 4368.9; r22 4384.8)
C_LIFT, C_OUTER = 862 * 7 + 12 + 6 * T33, 18 * 7 + T33
T_DOV = NQD * C_SP + C_LIFT + C_OUTER + C_HIGGS       # 1,847,747.1 (r24 937,592.7; r22 3,058,805.3)
T_DOV_24 = price_sp(NDOV, n_sgn=1)["t_dov"]           # 937,592.7 (r24)
# RECORD: the retired second-quantized construction at the worst-case step (the r21 box)
NDOV_SQ = 5 * N_BOUND + TPQ_R25                           # 416,215
NQ_SQ = NDOV_SQ * 691                                 # 287,604,565 BE queries
LCU_SITE = 3 * 2 * 2 * 2 + 4 * 2 + (4 * 2) * (4 + 1)  # 24 hop + 8 on-site + 40 Higgs = 72
C_TOFF = 7 * 27 * LCU_SITE                            # 13,608
LINK_MUX = 7 * (4 * 3 + 40) + (4 * 18 + 32)           # 4 multiplexers (3 Tof + 18 T) + spinor frame = 468
N33_SQ = 25 * NQ_SQ                                   # 25 rotations per query
T33_SQ = fit(N33_SQ)                                  # 31.848
C_BE33 = C_TOFF + 81 * LINK_MUX + 27 * 392 + 25 * T33_SQ  # 62,896.2
# the draft's diagonalizing frame at the same query count
LINK_CT, LINK_ROT = 7 * 384 + 144, 704                # 2O Wilson d=3 hop per link per query (r19: squish per link)
LINK_CT_DRAFT = 7 * 1182 + 144                        # r17 per-application squish (share="draft"), 8418
NROTQ = 25 + 81 * LINK_ROT                            # 57,049
T33D = fit(NROTQ * NQ_SQ)                             # 38.262
EPS33D = math.sqrt(1e-2 / (NROTQ * NQ_SQ))
C_BE33D = C_TOFF + 81 * (LINK_CT + LINK_ROT * T33D) + 27 * 392 + 25 * T33D    # 2,436,414.9
# records: the pre-r21 query count (30 x 5 + 240 = 390 D_ov) and the 1.7e3-term LCU line
NQ_PRE = 390 * 691                                    # 269,490
T_PRE = fit(NROTQ * NQ_PRE)                           # 32.478
C_BE_R20 = 11900 + 81 * LINK_CT + 27 * 392 + NROTQ * T_PRE        # 2,104,718.8 (the r19-r20 centre)
T_PRE_F = fit(25 * NQ_PRE)                            # 26.063
C_FLOOR_ZERO = 11900 + 25 * T_PRE_F                   # 12,551.6 (pre-r20 zero-T floor)
C_FLOOR = 11900 + 81 * LINK_MUX + 27 * 392 + 25 * T_PRE_F          # 61,043.6 (the r20 floor)


# --- PUBLISHED is the box as printed (rule 1) -------------------------------

def test_published_is_the_box_as_printed():
    assert m.PUBLISHED["2028"].lq == (250, 250)                  # r25 R10: 42 ancilla (r22-r24 240, r20-r21 247)
    assert m.PUBLISHED["2028"].hard_ops == (2.7e5, 2.9e5)       # r21 a (r17 2.8-3.0e5)
    assert m.PUBLISHED["2033"].lq == (1015, 1015)                # r25: 38 + 5 fan-out (r22-r24 1010, r21 1005)
    assert m.PUBLISHED["2033"].hard_ops == (1.9e10, 1.9e10)     # 17 filter calls (C1 1.7e10; r25 1.6e10, r24 8.2e9)
    assert set(m.PUBLISHED) == {"2028", "2033"}
    assert m.DISPUTED == {}


def test_both_boxes_reproduce(r28, r33):
    for era, r in (("2028", r28), ("2033", r33)):
        p = m.PUBLISHED[era]
        assert c.close(r.lq, p.lq, p.rel_tol)
        assert c.close(r.hard_ops, p.hard_ops, p.rel_tol), (era, r.hard_ops, p.hard_ops)


def test_no_codesign_box(a):
    with pytest.raises(ValueError):
        m.model(a, "codesign")


def test_every_assumption_is_read_by_the_model(a):
    # No decorative inputs: each field name must appear in the model body (rule 3 / verifier finding).
    src = inspect.getsource(m._model_2028) + inspect.getsource(m._model_2033) + \
        inspect.getsource(m._hwp_groups_2028) + inspect.getsource(m._step_cost_2028) + \
        inspect.getsource(m._spatial_hop_rule) + inspect.getsource(m._step_rtol_2028)
    for f in dataclasses.fields(a):
        assert f"a.{f.name}" in src, f"Assumptions.{f.name} is never read by model()"


# --- shared conventions -----------------------------------------------------

def test_uses_shared_synthesis_and_toffoli_conventions(a, r28, r33):
    # R5: 7 T per Toffoli everywhere
    assert a.toffoli_convention.value == "textbook" and c.T_PER_TOFFOLI["textbook"] == 7
    assert a.toffoli_convention.prov is c.Provenance.CITED
    assert r28.intermediates["t_per_toffoli"] == 7 and r33.intermediates["t_per_toffoli"] == 7
    for r in (r28, r33):
        for p in r.breakdown:
            if "toffoli" in p.name:
                assert p.t_each == 7, p.name
    # R-TOL with the FULL fit everywhere (r17: no slope-only price in any current number)
    assert a.eps_syn.value == c.EPS_SYN == 1e-2 and a.eps_syn.prov is c.Provenance.CITED
    assert a.rot_synthesis.value == "rus" and a.hop_synthesis.value == "rus"
    assert a.hop_synthesis_legacy.value == "rus-slope"
    i28, i33 = r28.intermediates, r33.intermediates
    assert (i28["n_rot_shot"], i28["chunked_n_rot_shot"]) == (N28, N28C) == (6090, 7056)
    assert i28["eps_rot"] == c.eps_rot_for(6090) and i28["chunked_eps_rot"] == c.eps_rot_for(7056)
    assert math.isclose(i28["t_per_rot_chapter"], T_RUS) and math.isclose(i28["chunked_t_per_rot_chapter"], T_RUS_C)
    assert math.isclose(i28["t_per_rot_rhop"], T_RUS) and math.isclose(i28["chunked_t_per_rot_rhop"], T_RUS_C)
    assert i33["n_rot_shot"] == N33 == 25248067 and i33["eps_rot"] == c.eps_rot_for(N33)
    assert math.isclose(i33["t_per_rot"], T33) and round(T33, 1) == 27.2 and round(EPS33 * 1e5, 1) == 2.0
    assert rel(N33, 2.5e7) <= 0.02                               # 'the 2.5e7 rotations of the sphaleron shot'
    # r24: the condensate arm is its own circuit (the preparation alone) with its own tolerance
    nc = TPQ * NROT_APP
    assert i33["n_rot_shot_condensate"] == nc == 3796117 and i33["eps_rot_condensate"] == c.eps_rot_for(nc)
    assert math.isclose(i33["t_per_rot_condensate"], fit(nc))
    assert i33["sq_draft_frame_n_rot_shot"] == NROTQ * NQ_SQ and math.isclose(i33["sq_draft_frame_t_per_rot"], T33D)
    for n, e in ((N28, EPS28), (N28C, EPS28C), (N33, EPS33)):
        assert math.isclose(n * e ** 2, 1e-2)
    # the retired fixed constants are records only (Ch. 6 reads them); no current number uses 30 T
    assert a.t_per_rot.value == 30 and a.eps_rot.value == 1e-4 and "RETIRED" in a.t_per_rot.note
    assert not any(p.t_each == 30 for r in (r28, r33) for p in r.breakdown)
    assert GROUPS["Z3"].link_qubits == 2 and GROUPS["2O"].link_qubits == 6
    assert not GROUPS["2O"].link_width_conflict
    assert a.faults_per_shot.value == 0.1
    # shot time: 1 us per T + 0.1 ms per shot
    assert a.t_gate_s.value == 1e-6 and a.shot_overhead_s.value == 1e-4


def test_hwp_formulas():
    assert m.hwp_synth_rotations(96) == 7
    assert m.hwp_synth_rotations(72) == 7
    assert m.hwp_synth_rotations(32) == 6
    assert m.hwp_synth_rotations(12) == 4
    assert m.hwp_synth_rotations(192) == 8
    assert m.hwp_synth_rotations(8) == 4
    assert [m.hwp_toffolis(k) for k in (96, 72, 32, 12, 8, 192)] == [94, 70, 31, 10, 7, 190]
    for k in (3, 8, 12, 32, 72, 96, 192):
        assert m.hwp_toffolis(k) <= k
    # E26 (r22): one ancilla per adder, k - w(k); the weight bits are those ancilla and one data qubit
    assert m.hwp_ancilla(96) == 94 and m.hwp_ancilla(32) == 31 and m.hwp_ancilla(192) == 190
    assert all(m.hwp_ancilla(k) == m.hwp_toffolis(k) == c.hwp_ancilla(k) for k in (3, 8, 9, 12, 32, 72, 96, 192))
    assert all(m.hwp_synth_rotations(k) <= m.hwp_ancilla(k) + 1 for k in (3, 8, 12, 32, 72, 96, 192))
    assert m._chunk(96, 2, 32) == [(32, 6)] and m._chunk(12, 16, 32) == [(12, 16)]
    assert m._chunk(72, 2, 32) == [(32, 4), (8, 2)]
    assert m._chunk(192, 1, 32) == [(32, 6)]
    assert m._chunk(96, 2, None) == [(96, 2)] and m._chunk(40, 1, 32) == [(32, 1), (8, 1)]


# --- 2028: logical qubits ---------------------------------------------------

def test_2028_lq_components_exact(r28):
    i = r28.intermediates
    assert i["lq_fermion"] == 192
    assert i["lq_gauge"] == 16
    # r25 (R10): three spatial-hop groups per slot, one RUS ancilla per weight bit: 3 x (10 + 4) = 42
    assert i["lq_ancilla"] == 3 * (10 + 4) == 42 and i["lq_ancilla_asserted"] == 40
    assert i["lq_ancilla_serial_schedule"] == 31 + 1 == 32 and i["lq_serial_schedule"] == 240   # E26, r22-r24
    assert i["ancilla_itemized"] == (("hwp_chunked_k32", 31), ("rus_synthesis", 1)) and i["anc_rus"] == c.RUS_ANCILLA
    assert i["lq_total"] == 250 == i["lq_stated"] and r28.lq == (250, 250)
    assert i["lq_single_color_copy"] == 64 + 16 + 32 == 112
    assert rel(112, i["lq_single_color_copy_stated"]) <= 0.02                       # '~110 LQ a single copy'


# --- 2028: the Trotter step -------------------------------------------------

def test_2028_fifth_direction_hops_open_boundary(a, r28):
    # r21 (3a), derived here: domain-wall fermions have open boundaries in the fifth direction, V (L5 - 1) hops
    i = r28.intermediates
    assert a.fifth_boundary.value == "open" and a.fifth_boundary.prov is c.Provenance.CITED
    assert "DERIVED HERE" in a.fifth_boundary.note
    assert i["fifth_boundary"] == "open" and i["n_fifth_hops"] == 8 * (4 - 1) == 24
    # a chain of L5 sites has L5 - 1 bonds; a ring has L5
    for L5 in (2, 4, 8):
        b = dataclasses.replace(a, L5_2028=c.Assumed(L5, "test"))
        assert m.model(b, "2028").intermediates["n_fifth_hops"] == 8 * (L5 - 1)
    # the pre-r21 count is the 'periodic' switch and reproduces the r17-r20 print
    b = dataclasses.replace(a, fifth_boundary=c.Cited("periodic", "test"))
    rb = m.model(b, "2028")
    assert rb.intermediates["n_fifth_hops"] == 32
    assert tuple(round(x) for x in rb.hard_ops) == (276051, 297373)
    assert tuple(round(x) for x in i["t_per_shot_pre_r21"]) == (276051, 297373)
    with pytest.raises(ValueError):
        dataclasses.replace(a, fifth_boundary=c.Cited("antiperiodic", "test"))


def test_2028_term_counts(r28):
    i = r28.intermediates
    assert i["n_spatial_hops"] == 32 and i["n_fifth_hops"] == 24 and i["n_electric_terms"] == 8
    assert i["n_terms_per_step"] == 64                                  # '64 hopping and electric terms'
    assert i["n_hop_rotations_per_step_before_rhop"] == (32 + 24) * 3 * 2 == 336
    assert i["n_hop_rotations_per_step"] == 24 * 3 * 2 + 32 * 3 * 8 == 912
    assert i["n_onsite_rotations_per_step"] == 2 * 3 * 8 * 4 == 192     # 2N per site x V L5 sites (r17 f)
    assert i["rotations_per_step"] == 912 + 24 + 192 == 1128
    assert rel(i["rotations_per_step"], 1e3) <= 0.2                     # '~1e3 Pauli rotations per step'


def test_2028_grouping_is_the_chapters_sentence(r28):
    groups = dict((n, (k, ng)) for n, k, ng in r28.intermediates["groups"])
    assert groups["fifth_hops"] == (72, 2)                               # r21 a (was 2 x k=96)
    assert groups["spatial_hops"] == (12, 64)
    assert groups["electric"] == (8, 3)
    assert groups["onsite_dw"] == (192, 1)                               # one lattice-wide group (r17 f)
    assert r28.intermediates["grouping"] == "color+s" and r28.intermediates["onsite_term"] == "wilson"
    assert sum(k * ng for k, ng in groups.values()) == r28.intermediates["rotations_per_step"] == 1128


def test_2028_spatial_hop_full_fit_and_central_function(a, r28):
    # Z3: no draft entry; R-HOP structure (C_W = 0, 8 strings, k = 12) at the full fit (r17 e)
    i = r28.intermediates
    rh = rhop_link("Z3", n_stag=4, n_c=3, eps=EPS28, k=12, synthesis="rus")
    assert rh["n_moves"] == 0 and rh["t_moves"] == 0 and rh["n_pauli"] == 8 and rh["k"] == 12
    assert i["rhop_n_pauli_spatial"] == 8 and i["rhop_k_spatial"] == 12 and i["rhop_t_move_spatial"] == 0
    assert i["hop_synthesis"] == "rus" and a.hop_synthesis.prov is c.Provenance.CITED
    assert math.isclose(i["rhop_t_rot"], T_RUS)
    assert math.isclose(i["rhop_hwp_t"], 4 * T_RUS + 70)
    assert i["rhop_hwp_ancilla"] == 10 <= 31
    assert math.isclose(i["rhop_t_per_link"], 8 * (4 * T_RUS + 70)) and math.isclose(rh["t"], i["rhop_t_per_link"])
    # the central r17 function gives the same link
    hl = hop_link_cost("Z3", n_stag=4, eps=EPS28, n_c=3, k=12, legacy=False)
    assert math.isclose(hl["t"], rh["t"]) and math.isclose(i["hop_link_cost_z3_t"], hl["t"])
    assert hl["toffoli"] == 80 and hl["n_rot"] == 32
    assert i["rhop_n_links"] == 8
    assert math.isclose(i["t_hop_rhop_per_step"], 8 * rh["t"])                  # 9663.8
    assert math.isclose(i["chunked_t_hop_rhop_per_step"], 64 * (4 * T_RUS_C + 70))


def test_2028_step_arithmetic_by_hand(r28):
    # 14 + 256 + 12 + 8 = 290 synthesized, 140 + 640 + 21 + 190 = 991 Toffolis (r21 a: two groups of k = 72)
    i = r28.intermediates
    assert i["n_synth_per_step"] == 14 + 256 + 12 + 8 == 290
    assert i["n_synth_chapter_per_step"] == 34 and i["n_synth_rhop_per_step"] == 256
    assert i["n_hwp_toffoli_per_step"] == 140 + 640 + 21 + 190 == 991
    assert math.isclose(i["t_synth_per_step"], 290 * T_RUS) and i["t_hwp_per_step"] == 991 * 7 == 6937
    assert math.isclose(i["t_per_step"], 290 * T_RUS + 991 * 7)                # 12,809.3
    assert i["n_trotter"] == 20 and i["n_steps_total"] == 21
    assert math.isclose(i["t_per_shot_unchunked"], 21 * i["t_per_step"]) and round(i["t_per_shot_unchunked"]) == 268995
    assert i["t_per_shot"] == i["chunked_t_per_shot"]                          # r25 (R9): the shot the register runs
    assert rel(i["t_per_step"], i["t_per_step_stated"][0]) <= 0.01             # '~1.28e4 T per step'
    assert rel(i["t_per_shot_unchunked"], i["t_per_shot_stated"][0]) <= 0.01   # '~2.69e5 unchunked'
    assert round(T_RUS, 1) == 20.2 and round(EPS28, 4) == 1.3e-3               # '290 x 20.2 + 991 x 7'
    assert rel(290 * 20.2 + 991 * 7, 1.28e4) <= 0.01
    # records: without the on-site term (2028-A), pre-r17 (slope-only hops), before R-TOL, all at the old 32 hops
    assert tuple(round(x, 1) for x in i["t_per_shot_no_onsite"]) == (244581.6, 254029.9)
    assert tuple(round(x) for x in i["t_per_shot_pre_r17"]) == (195122, 204571)
    assert tuple(round(x) for x in i["t_per_shot_before_rtol"]) == (223333, 236899)
    # 'The on-site term adds 3.1-4.3e4 T of this'
    lo, hi = i["onsite_term_adds_per_shot"]
    assert rel(lo, 3.1e4) <= 0.02 and rel(hi, 4.3e4) <= 0.02


def test_2028_onsite_term_switch(a):
    # 'none' with the pre-r21 hop count reproduces 2028-A (full-fit hops, no on-site term) as the headline
    b = dataclasses.replace(a, onsite_term_2028=c.Cited("none", "test"), fifth_boundary=c.Cited("periodic", "test"))
    rb = m.model(b, "2028")
    assert tuple(round(x, 1) for x in rb.hard_ops) == (244581.6, 254029.9)
    with pytest.raises(ValueError):
        dataclasses.replace(a, onsite_term_2028=c.Cited("staggered", "test"))


def test_2028_chunked_variant_by_hand(r28):
    i = r28.intermediates
    assert i["hwp_chunk_k"] == 32
    groups = [(n, k, ng) for n, k, ng in i["chunked_groups"]]
    # k = 72 = 2 x 32 + 8 per Pauli string (r21 a)
    assert ("fifth_hops", 32, 4) in groups and ("fifth_hops", 8, 2) in groups
    assert ("spatial_hops", 12, 64) in groups and ("electric", 8, 3) in groups
    assert ("onsite_dw", 32, 6) in groups
    assert i["chunked_n_synth_per_step"] == (24 + 8) + 256 + 12 + 36 == 336
    assert i["chunked_n_hwp_toffoli_per_step"] == (124 + 14) + 640 + 21 + 186 == 985
    assert math.isclose(i["chunked_t_per_step"], 336 * T_RUS_C + 985 * 7)     # 13,739.8
    assert math.isclose(i["chunked_t_per_shot"], 21 * i["chunked_t_per_step"]) and round(i["chunked_t_per_shot"]) == 288535
    assert i["chunking_delta_n_synth"] == 46 and i["chunking_delta_n_toffoli"] == -6
    assert rel(i["chunked_t_per_step"], i["t_per_step_stated"][1]) <= 0.01     # '~1.37e4'
    assert rel(i["chunked_t_per_shot"], i["t_per_shot_stated"][1]) <= 0.01     # '~2.89e5 chunked'
    assert round(T_RUS_C, 1) == 20.4 and round(EPS28C, 5) == 1.19e-3           # '336 x 20.4 + 985 x 7'


def test_2028_hard_ops_is_the_printed_range(r28):
    i = r28.intermediates
    assert r28.hard_ops == (i["t_per_shot_unchunked"], i["chunked_t_per_shot"])
    assert tuple(round(x) for x in r28.hard_ops) == (268995, 288535)
    assert not i["fits_2028_envelope"]
    lo, hi = i["ratio_to_envelope"]
    slo, shi = i["ratio_to_envelope_stated"]
    assert (slo, shi) == (2.7, 2.9)                                       # 'a factor 2.7-2.9 above'
    assert rel(lo, slo) <= 0.01 and rel(hi, shi) <= 0.01
    assert math.isclose(r28.breakdown_total(), r28.hard_ops[0])


def test_2028_ancilla_itemized(r28):
    i = r28.intermediates
    assert i["hwp_ancilla_unchunked"] == 190 == i["hwp_ancilla_unchunked_stated"]     # 'The k=192 group holds 190'
    assert i["hwp_ancilla_chunked"] == 31 == i["hwp_ancilla_chunked_stated"]          # 'needs 31 and fits'
    assert i["chunked_fits_box_ancilla"] and not i["unchunked_fits_box_ancilla"]
    assert i["lq_unchunked"] == 192 + 16 + 190 + 1 == 399                 # 'the unchunked end needs 191'
    assert rel(i["lq_unchunked"], i["lq_unchunked_stated"]) <= 0.01      # 'about 400 LQ'


def test_2028_ramp_length_trade(a, r28):
    i = r28.intermediates
    assert i["n_steps_total_before_r6"] == 31
    lo, hi = i["t_per_shot_at_30_step_ramp"]
    assert math.isclose(lo, 31 * i["t_per_step"]) and math.isclose(hi, 31 * i["chunked_t_per_step"])
    b = dataclasses.replace(a, n_trotter_2028=c.Stated(30, "test"))
    ib = m.model(b, "2028").intermediates
    assert ib["n_rot_shot"] == 290 * 31 and math.isclose(ib["eps_rot"], math.sqrt(1e-2 / (290 * 31)))
    assert all(x > y for x, y in zip(m.model(b, "2028").hard_ops, (lo, hi)))


def test_2028_split_is_explicit_and_sums(r28):
    i = r28.intermediates
    assert math.isclose(i["t_synth_per_step"] + i["t_hwp_per_step"], i["t_per_step"])
    assert i["t_hwp_per_step"] <= i["t_hwp_bound_k_per_step"]
    # 'Synthesized rotations are a little under half of each step'
    assert 0.4 < i["synth_fraction_of_step"] < 0.5                        # 46%
    assert 0.4 < i["chunked_synth_fraction_of_step"] < 0.5                # 50%


def test_2028_toffoli_convention_moves_the_step(a):
    b = dataclasses.replace(a, toffoli_convention=c.Cited("jones", "test"))
    i = m.model(b, "2028").intermediates
    assert math.isclose(i["t_per_step"], 290 * T_RUS + 991 * 4)


def test_2028_other_groupings_are_emitted(r28):
    o = r28.intermediates["other_groupings"]
    assert set(o) == {"color", "color+site"}
    assert o["color"]["n_synth_per_step"] == 14 + 512 + 12 + 8 == 546
    assert o["color"]["n_hwp_toffoli_per_step"] == 140 + 256 + 21 + 190 == 607
    assert o["color"]["n_rot_shot"] == 546 * 21
    assert math.isclose(o["color"]["t_per_step"], 546 * fit(546 * 21) + 607 * 7)
    assert o["color+site"]["n_synth_per_step"] == 14 + 56 + 12 + 8 == 90
    assert o["color+site"]["n_hwp_toffoli_per_step"] == 140 + 752 + 21 + 190 == 1103
    assert math.isclose(o["color+site"]["t_per_step"], 90 * fit(90 * 21) + 1103 * 7)


def test_2028_coherent_synthesis_sensitivity(r28):
    i = r28.intermediates
    assert math.isclose(i["eps_total_synthesis_incoherent"], 1e-2)
    assert math.isclose(i["eps_rot_coherent_same_total"], 1e-2 / 6090)             # 1.64e-6
    t_coh = 1.15 * math.log2(6090 / 1e-2) + 9.2                                     # 31.30
    assert math.isclose(i["t_per_rot_coherent_model"], t_coh)
    assert math.isclose(i["t_per_rot_hop_coherent"], t_coh)                        # hops at the full fit too
    assert i["t_per_rot_coherent_stated"] == 31 and round(t_coh, 1) == 31.3
    # 'the unchunked step to 290 x 31.3 + 991 x 7 ~ 1.60e4 T, and the shot to ~3.4e5 T'
    assert math.isclose(i["t_per_step_coherent_model"], 290 * t_coh + 991 * 7)    # 16,013.6
    assert rel(i["t_per_step_coherent_model"], i["t_per_step_coherent_stated"]) <= 0.01
    assert rel(i["t_per_step_coherent_at_stated"], i["t_per_step_coherent_stated"]) <= 0.01
    assert math.isclose(i["t_per_shot_coherent_model"], 21 * i["t_per_step_coherent_model"])      # 336,285
    assert rel(i["t_per_shot_coherent_model"], i["t_per_shot_coherent_stated"]) <= 0.02
    assert rel(i["t_per_rot_coherent_model"], i["t_per_rot_coherent_stated"]) <= 0.02


def test_2028_1p1d_overlap_in_the_single_particle_architecture(a, r28):
    # r22 (E26 (i)), derived here: the 1+1D Z3 overlap route, priced as the 2033 query (arxiv_2607_28524 architecture)
    i = r28.intermediates
    q = i["overlap_1p1d_query"]
    assert i["overlap_1p1d_q_modes"] == 2 * 3 * 8 == 48 and i["overlap_1p1d_index_register"] == 2 + 1 + 3 == 6
    assert q["n_lcu_unitaries"] == 5 and q["term_register"] == 3
    assert q["toffoli_items"] == {"backward_flag": 1, "controlled_shifts": 6, "dirac_flag": 1, "read_iteration": 7,
                                  "link_copy": 16, "fixup_iteration": 7, "projector": 2}
    assert q["toffoli"] == 40 == i["overlap_1p1d_query_toffoli"]
    assert q["rot_items"] == {"prepare_and_inverse": 2, "link_phase": 2, "qsp_phase_controlled": 2}
    assert q["n_rot"] == 6 == i["overlap_1p1d_query_n_rot"]
    assert m.single_particle_query_Z3_1p1d(8, 2) == q
    assert i["overlap_1p1d_lift"]["toffoli"] == 2 * (2 * 48 - 1) == 190 and i["overlap_1p1d_lift"]["n_rot"] == 2
    # each circuit at its own R-TOL tolerance: one application, or a ground-state preparation of 10
    d30, d50 = math.ceil(30 * math.log(100)), math.ceil(50 * math.log(100))
    assert (i["overlap_1p1d_d_sgn_kappa30"], i["overlap_1p1d_d_sgn_kappa50"]) == (d30, d50) == (139, 231)

    def app(d, n_sgn):
        tr = fit(n_sgn * ((d + 1) * 6 + 2))
        c_sp = 40 * 7 + 6 * tr
        return c_sp, n_sgn * ((d + 1) * c_sp + 190 * 7 + 4 + 2 * tr)
    assert math.isclose(i["overlap_1p1d_c_sp"], app(139, 1)[0])
    assert all(rel(x, i["overlap_1p1d_c_sp_stated"]) <= 0.03 for x in i["overlap_1p1d_c_sp_all"])   # '~4.0e2 T'
    assert math.isclose(i["overlap_1p1d_sgn_T_kappa30"], app(139, 1)[1])                # 5.620e4
    assert rel(i["overlap_1p1d_sgn_T_kappa30"], i["overlap_1p1d_sgn_T_kappa30_stated"]) <= 0.01      # '5.6e4'
    assert math.isclose(i["overlap_1p1d_gs_prep_T"], app(139, 10)[1])                   # 5.781e5
    assert rel(i["overlap_1p1d_gs_prep_T"], i["overlap_1p1d_gs_prep_T_stated"]) <= 0.01              # '~5.8e5'
    assert math.isclose(i["overlap_1p1d_sgn_T_kappa50"], app(231, 1)[1])                # 9.28e4
    assert rel(i["overlap_1p1d_sgn_T_kappa50"], i["overlap_1p1d_sgn_T_kappa50_stated"]) <= 0.04      # '~9e4'
    assert math.isclose(i["overlap_1p1d_gs_prep_T_kappa50"], app(231, 10)[1])           # 9.55e5
    assert rel(i["overlap_1p1d_gs_prep_T_kappa50"], i["overlap_1p1d_gs_prep_T_kappa50_stated"]) <= 0.05   # '~1e6'
    # the comparison: one application is more than half the 2028 budget; a ground-state preparation is about six
    # times the budget and twice the domain-wall shot
    assert 0.5 < i["overlap_1p1d_sgn_over_envelope"] < 0.6 and i["overlap_1p1d_sgn_T_kappa30"] < 1e5
    assert math.isclose(i["overlap_1p1d_gs_prep_over_dwf"], i["overlap_1p1d_gs_prep_T"] / r28.hard_ops[1])
    assert rel(i["overlap_1p1d_gs_prep_over_dwf"], i["overlap_1p1d_gs_prep_over_dwf_stated"]) <= 0.01      # 2.0 -> 'twice'
    assert rel(i["overlap_1p1d_gs_prep_over_envelope"], i["overlap_1p1d_gs_prep_over_envelope_stated"]) <= 0.04   # 5.8 -> 'six'
    assert 0.95 <= i["overlap_1p1d_orders_past_envelope"] <= 1.05              # 'an order of magnitude past'
    lo, hi = i["overlap_1p1d_d_sgn_cut_to_fit"]
    assert 5 <= lo <= hi <= 10                                                 # 'd_sgn cut by an order of magnitude'
    # RECORD (r21 3e): the retired second-quantized pricing, 240 Pauli strings + 25 rotations per query
    sq = i["overlap_1p1d_second_quantized_record"]
    per = sq["lcu_per_site"]
    assert per["hop"] == 1 * 1 * 3 * 1 * 8 == 24 and per["onsite"] == 2 * 3 == 6 and per["higgs"] == 0
    assert per["total"] == 30 and sq["lcu_terms"] == 8 * 30 == 240 and m.lcu_terms_per_site(1, 1, 3, 1, 8) == per
    assert sq["c_toffoli_T"] == 240 * 7 == 1680 and a.lcu_toffoli_per_term.value == 1 and a.n_rot_be.value == 25

    def cbe(d, n_sgn):
        return 1680 + 25 * fit(25 * d * n_sgn)
    assert math.isclose(sq["c_be"], cbe(139, 1)) and round(sq["c_be"]) == 2175 and rel(sq["c_be"], sq["c_be_stated"]) <= 0.02
    assert math.isclose(sq["sgn_T_kappa30"], 139 * cbe(139, 1)) and rel(sq["sgn_T_kappa30"], 3.0e5) <= 0.01
    assert math.isclose(sq["gs_prep_T"], 10 * 139 * cbe(139, 10)) and rel(sq["gs_prep_T"], 3.1e6) <= 0.01
    assert math.isclose(sq["sgn_T_kappa50"], 231 * cbe(231, 1)) and rel(sq["sgn_T_kappa50"], 5e5) <= 0.01


def test_2028_shots_and_wall_time(a, r28):
    # r21 (3d), derived here: sigma^2 ~ 10 is the 27-site average of the 2033 lattice; the 2028 estimator averages
    # N_c V = 24 site-colour copies
    i = r28.intermediates
    assert i["condensate_copies"] == (27, 3 * 8)
    assert math.isclose(i["sigma2_condensate_2028"], 10 * 27 / 24) and math.isclose(i["sigma2_condensate_2028"], 11.25)
    assert math.isclose(i["shots"], 11.25 / 0.1 ** 2) and round(i["shots"]) == 1125
    assert rel(i["shots"], i["shots_stated"]) <= 0.03                          # '1.1e3'
    # operator-norm bound: fully correlated sites leave only the 3 independent colour copies
    assert math.isclose(i["sigma2_condensate_2028_worst"], 10 * 27 / 3) and math.isclose(i["shots_worst"], 9000)
    assert math.isclose(i["shots_worst"], i["shots_worst_stated"])
    # 'its variance is at most 1': the bilinear a^dag b + b^dag a has eigenvalues {-1, 0, 0, 1}; the chapter's
    # sigma^2 ~ 10 at 27 sites (270 per copy) then holds when the condensate per copy exceeds 6% of its maximum
    import numpy as np
    Z, sm = np.diag([1, -1]).astype(complex), np.array([[0, 1], [0, 0]], complex)
    a0, a1 = np.kron(sm, np.eye(2)), np.kron(Z, sm)
    ev = np.linalg.eigvalsh(a0.conj().T @ a1 + a1.conj().T @ a0)
    assert np.allclose(sorted(ev), [-1, 0, 0, 1]) and i["condensate_variance_per_copy_max"] == 1.0
    assert math.isclose(i["condensate_per_copy_min_for_sigma2"], 1 / math.sqrt(270)) and rel(1 / math.sqrt(270), 0.06) <= 0.02
    # it scales as 1 / copies: one colour copy would need three times the shots
    b = dataclasses.replace(a, Nc_2028=c.Assumed(1, "test"))
    assert math.isclose(m.model(b, "2028").intermediates["shots"], 3 * i["shots"])
    assert i["wall_s_per_shot_stated"] == (0.27, 0.29)
    lo, hi = i["wall_s_per_shot"]
    assert math.isclose(lo, r28.hard_ops[0] * 1e-6 + 1e-4) and math.isclose(hi, r28.hard_ops[1] * 1e-6 + 1e-4)
    assert rel(lo, 0.27) <= 0.01 and rel(hi, 0.29) <= 0.01
    assert math.isclose(i["wall_s_total"][0], i["shots"] * lo) and math.isclose(i["wall_s_total"][1], i["shots"] * hi)
    assert 5.0 * 60 < i["wall_s_total"][0] < i["wall_s_total"][1] < 5.45 * 60     # one L5 at 1125 shots (r24 box)
    assert r28.wall_time_s == i["wall_campaign_s"]                           # r25: campaign wall (L5 = 2, 3, 4)
    assert 13.0 * 60 < r28.wall_time_s[0] < 13.5 * 60
    assert r28.epsilon_l == (0.1 / r28.hard_ops[1], 0.1 / r28.hard_ops[0]) == i["eps_l_implied"]
    assert rel(r28.epsilon_l[0], i["eps_l_stated"]) <= 0.02                    # '<~ 3.5e-7'


# --- 2033: logical qubits ---------------------------------------------------

def test_2033_lq_components_exact(r33):
    i = r33.intermediates
    assert i["lq_fermion"] == 216
    assert i["lq_gauge"] == 486
    assert i["lq_higgs"] == 270
    # r25 (R10): 38 itemized + 5 for the fan-out of the shared control (r22-r24 38, r21 33, r20 43)
    assert i["lq_ancilla_itemized"] == 38 and i["lq_ancilla_fanout"] == 5 and i["depth_variant"] == "fanout"
    assert i["lq_ancilla"] == 43 and i["lq_ancilla_asserted"] == 110
    assert i["lq_total"] == 1015 and r33.lq == (1015, 1015)
    assert i["qubits_per_site"] == 36
    assert i["lq_dwf_backup_L5_8"] == 2527 and rel(2527, 2.5e3) <= 0.02
    assert i["lq_dwf_backup_L5_16"] == 4255 and rel(4255, 4e3) <= 0.10
    # record: the retired second-quantized construction
    assert i["sq_lq_ancilla"] == 33 and i["sq_lq_total"] == 1005
    assert i["sq_lq_ancilla_draft_frame"] == 43 and i["sq_lq_total_draft_frame"] == 1015


# --- 2033: the single-particle overlap kernel (r22, E26 (i)) ------------------

def test_2033_d_sgn_is_natural_log_and_an_integer_degree(a, r33):
    i = r33.intermediates
    # r24: kappa = 30 at the same delta_sgn = 1e-3 (r17-r23: kappa 1e2, 691, printed ~700)
    assert a.kappa_2033.value == 30 and a.kappa_2033_r22.value == 1e2 and a.delta_sgn.value == 1e-3
    assert math.isclose(i["d_sgn_exact"], 30 * math.log(1e3)) and round(i["d_sgn_exact"], 2) == 207.23
    assert i["d_sgn"] == D_SGN == 208 == i["d_sgn_chapter"]
    assert i["d_sgn_r22"] == math.ceil(100 * math.log(1e3)) == 691
    assert rel(i["d_sgn_log2_would_be"], 299) <= 0.01
    # r22: kappa = alpha / Delta, alpha = 2D + |D - m_0| the block-encoding normalization at r = 1
    assert i["m0_wilson"] == 1.5 and i["alpha_hw"] == 3 + 3 + 1.5 == 7.5 and i["norm_hw_free"] == 4.5
    assert "alpha / Delta" in a.kappa_2033.note
    # the other reading (kappa = ||H_W|| / Delta at the same gap): degree 346, 1.64x the shot
    assert i["d_sgn_if_kappa_is_norm_over_gap"] == math.ceil(30 * math.log(1e3) * 7.5 / 4.5) == 346
    assert math.isclose(i["t_per_shot_if_kappa_is_norm_over_gap"], price_sp(NDOV, 346)["t_shot"])
    assert rel(i["t_per_shot_if_kappa_is_norm_over_gap"], 3.076e10) <= 0.01      # 17 calls (C1 2.831e10; r25 2.686e10, r24 1.355e10)
    with pytest.raises(ValueError):
        dataclasses.replace(a, m0_wilson=c.Assumed(2.5, "test"))


def test_2033_architecture_is_the_single_particle_lift(a, r33):
    i = r33.intermediates
    assert a.sgn_architecture.value == "single_particle_lift" == i["sgn_architecture"]
    assert a.sgn_architecture.prov is c.Provenance.CITED and "arxiv_2607_28524" in a.sgn_architecture.src
    assert "INVALID" in a.sgn_architecture.note
    # the second-quantized construction is not selectable: it is a record
    with pytest.raises(ValueError):
        dataclasses.replace(a, sgn_architecture=c.Cited("second_quantized", "test"))
    # Q = N^d N_I 2^((d+1)/2) = 27 x 2 x 4; index register q = n_I + n_S + n_X = 1 + 2 + 3 x 2
    assert i["q_modes"] == 27 * 2 * 4 == 216 == i["lq_fermion"]
    assert i["index_register"] == {"colour": 1, "spin": 2, "coordinates": 6} and i["index_register_qubits"] == 9
    assert i["n_lcu_unitaries"] == 4 * 3 + 1 == 13 and i["term_register_qubits"] == 5
    assert i["be_normalization"] == 2 * 216 == 432


def test_2033_single_particle_query_itemized(a, r33):
    # r22, derived here: one query of the single-particle H_W block encoding
    i = r33.intermediates
    q = i["sp_query"]
    assert q["toffoli_items"] == {"direction_flags": 2, "backward_flags": 3, "forward_flags": 3, "controlled_shifts": 12,
                                  "dirac_flags": 3, "read_iteration": 3 * 26, "link_copy": 5 * 81, "mux_orientation": 10,
                                  "fixup_iteration": 3 * 26, "projector": 4}
    assert q["toffoli"] == 598 == i["sp_query_toffoli"] == i["sp_query_toffoli_stated"]
    assert q["t_direct_items"] == {"mux_and_inverse": 20, "prepare_controlled_h": 4} and q["t_direct"] == 24
    assert q["rot_items"] == {"prepare_and_inverse": 4, "qsp_phase_controlled": 2} and q["n_rot"] == 6
    assert (i["sp_query_t_direct"], i["sp_query_n_rot"]) == (24, 6)
    # 'this takes 78 + 405 + 78: 561 of the query's 598 Toffolis'; 'the other 37'
    assert q["link_read_toffoli"] == 78 + 405 + 78 == 561 == i["sp_link_read_toffoli"] and 598 - 561 == 37
    assert m.single_particle_query_2O(3, 3, 6) == q
    # the multiplexer on a colour qubit: 10 T, no Toffoli, no ancilla
    assert i["mux_index_per_bit"] == {"x1": (0, 0), "x2": (0, 0), "x3": (0, 0), "x4": (0, 4), "x5": (0, 4), "x6": (0, 2)}
    assert (i["mux_index_toffoli"], i["mux_index_t_direct"], i["mux_index_ancilla"]) == (0, 10, 0)
    # the price of a query at the shot's tolerance
    assert math.isclose(i["c_sp"], C_SP) and round(i["c_sp"], 1) == 4373.0
    assert rel(i["c_sp"], i["c_sp_stated"]) <= 0.01                             # '~4.4e3 T'
    assert rel(598 * 7 + 24 + 6 * 27.2, 4.4e3) <= 0.01                         # the printed arithmetic
    assert 0.89 < i["sp_link_read_share_of_query"] < 0.905                      # '90% of it the link read'
    assert 0.035 < i["sp_rotation_share_of_query"] < 0.045                      # 'about 4% of a query' (3.7%)
    # controlled shifts: mod 3 on two bits is 2 Toffolis; mod 2^n is n (n - 1) / 2; anything else is not derived
    assert m._ctrl_shift_toffolis(3) == 2 and m._ctrl_shift_toffolis(8) == 3 and m._ctrl_shift_toffolis(4) == 1
    with pytest.raises(ValueError):
        m._ctrl_shift_toffolis(5)
    # the link read scales with the number of links (the address is in superposition)
    q4 = m.single_particle_query_2O(3, 4, 6)
    assert q4["toffoli_items"]["link_copy"] == 5 * 192 and q4["toffoli_items"]["read_iteration"] == 3 * 63
    b = dataclasses.replace(a, L_2033=c.Assumed(4, "test"))
    assert m.model(b, "2033").intermediates["sp_query_toffoli"] == q4["toffoli"]


def test_2033_one_application_components(a, r33):
    # one application of the block-encoded psi^dag h_ov psi / (2Q): r25 (a) 209 + 208 queries (h_ov and the Yukawa
    # projector) + lift + reflection + Higgs once
    i = r33.intermediates
    assert i["yukawa_placement"] == "inside" and i["n_sgn_per_application"] == 2
    assert i["n_queries_per_dov"] == 208 + 1 + 208 == NQD == 417
    assert math.isclose(i["t_sgn_per_dov"], 417 * C_SP) and rel(i["t_sgn_per_dov"], 1.8e6) <= 0.02   # '1.8e6 T'
    assert i["lift"] == {"toffoli": 2 * (2 * 216 - 1), "t_direct": 12, "n_rot": 6, "leaves": 432, "anc_iteration": 10}
    assert i["lift"]["toffoli"] == 862 and m.lift_cost(216, 3) == i["lift"]
    assert math.isclose(i["c_lift"], C_LIFT) and rel(i["c_lift"], 6.2e3) <= 0.01                     # '6.2e3 T'
    assert i["n_be_ancilla"] == 2 + 9 + 1 + 1 + 1 + 5 == 19 and i["outer_reflection"] == {"toffoli": 18, "n_rot": 1}
    assert math.isclose(i["c_outer"], C_OUTER)
    # the Higgs vertex: outside sgn, once per application; 27 U_x + 40 strings per site
    assert i["n_umult_per_higgs_site"] == 1 and i["c_higgs_umult_T"] == 27 * 392 == 10584
    assert i["higgs_strings_per_site"] == 8 * (4 + 1) == 40 == i["higgs_strings_per_site_stated"]
    assert i["higgs_select_toffolis"] == 27 * 40 == 1080
    assert i["c_higgs_per_dov"] == C_HIGGS == 18144 and rel(18144, 1.8e4) <= 0.01                    # '1.8e4 T'
    assert i["n_rot_per_dov"] == NROT_APP == 2509
    assert math.isclose(i["t_per_dov"], T_DOV) and round(i["t_per_dov"], 1) == 1848028.6
    assert rel(i["t_per_dov"], i["t_per_dov_stated"]) <= 0.03                                         # '1.8e6 T'
    assert 0.98 < i["sgn_share_of_dov"] < 0.99 and 0.003 < i["lift_share_of_dov"] < 0.004
    assert 0.009 < i["higgs_share_of_dov"] < 0.01
    # r25 (a): 1.97x the r24 application (Yukawa outside), whose numbers stay reproducible
    assert rel(i["yukawa_over_r24_per_application"], 1.9710) <= 1e-4 and rel(i["r24_t_per_dov"], 937592.7) <= 1e-7
    assert rel(i["r24_t_per_shot"], 8.241440e9) <= 1e-6 and rel(i["r24_t_per_shot_condensate"], 2.24118e8) <= 1e-5
    t_coh = 1.15 * math.log2(N33 / 1e-2) + 9.2
    assert math.isclose(i["t_per_rot_coherent"], t_coh)
    assert math.isclose(i["t_dov_shift_coherent"], NROT_APP * (t_coh - T33) / T_DOV)
    assert i["t_group_action_retired"] == 190


# --- 2033: the Trotter step count (r22, E26 (ii)) ---------------------------

def test_2033_trotter_steps_state_dependent(a, r33):
    # headline: the state-dependent second-order estimate in the thermal state, weak-coupling evaluation (derived here)
    i = r33.intermediates
    assert a.trotter_rule_2033.value == "state" == i["trotter_rule"]
    assert a.trotter_rule_2033.prov is c.Provenance.CITED and "AlvesLammLiu_inprep" in a.trotter_rule_2033.src
    assert "E26" in a.trotter_rule_2033.src and a.temp_lat_2033.prov is c.Provenance.ASSUMED
    sd = i["trotter_state_dependent"]
    # 156 transverse oscillators: 3 colours x 2 polarizations x 26 non-zero momenta; w in {sqrt3, sqrt6, 3}
    assert sd["n_momenta"] == 27 - 1 == 26 and sd["n_modes"] == 3 * 2 * 26 == 156 == i["trotter_n_modes_stated"]
    assert math.isclose(sd["w_max"], 3.0)
    ws = sorted({round(math.sqrt(4 * sum(math.sin(math.pi * n / 3) ** 2 for n in k)), 9)
                 for k in ((1, 0, 0), (1, 1, 0), (1, 1, 1), (2, 1, 0))})
    assert all(math.isclose(x, y) for x, y in zip(ws, (math.sqrt(3), math.sqrt(6), 3.0)))
    # by hand: 6 momenta at sqrt3, 12 at sqrt6, 8 at 3, six modes each
    mult = {math.sqrt(3): 6, math.sqrt(6): 12, 3.0: 8}

    def half(w):
        return 0.5 * w ** 3 / math.tanh(w / (2 * 0.5))
    mean = 6 * sum(n * half(w) for w, n in mult.items())
    var = 6 * sum(n * 2 * half(w) ** 2 for w, n in mult.items())
    assert math.isclose(sd["mean"], mean) and round(sd["mean"], 1) == 1287.8
    assert math.isclose(sd["std"], math.sqrt(var)) and round(sd["std"], 1) == 162.0
    assert math.isclose(sd["lambda_var"], math.sqrt(var) / 8) and round(sd["lambda_var"], 2) == 20.25
    assert math.isclose(sd["lambda_full"], math.sqrt(mean ** 2 + var) / 8) and round(sd["lambda_full"], 1) == 162.2
    assert rel(i["trotter_lambda_sd"], i["trotter_lambda_sd_stated"]) <= 0.02     # 'Lambda_sd ~ 20 a^-3'
    # r24: t = 2 fm = 10a (was 50a); same Lambda_sd, N scales as t^{3/2}
    assert i["t_over_a_2033"] == 10 and i["t_over_a_r22"] == 50
    assert sd["n_central"] == math.ceil(math.sqrt(sd["lambda_var"] * 10 ** 3 / 0.1)) == N_TROT == 450 == i["n_trotter_2033"]
    assert sd["n_high"] == 1274 and sd["n_low"] == math.ceil(10 * math.sqrt(27 * 10 / 2.4)) == 107
    assert i["n_trotter_range"] == (107, 1274)
    assert rel(i["n_trotter_2033"], i["n_trotter_2033_stated"]) <= 0.001          # '450 steps'
    assert rel(107, i["n_trotter_range_stated"][0]) <= 0.03 and rel(1274, i["n_trotter_range_stated"][1]) <= 0.03
    assert rel(i["trotter_dt_lattice"], 0.022) <= 0.02                            # 'dt ~ 0.022 a'
    # t^{3/2}: 5031 x 5^{-3/2} = 450.0 (ceil); the r22 count is kept as a record
    sd22 = i["trotter_state_dependent_r22"]
    assert sd22["n_central"] == N_TROT_R22 == 5031 and sd22["n_low"] == 1186 and sd22["n_high"] == 14242
    assert abs(5031 * 5 ** -1.5 - 450) < 1
    # 'N = (Lambda t^3 / eps)^(1/2) = 1.0e2 (Lambda a^3)^(1/2)'
    assert rel(math.sqrt(10 ** 3 / 0.1), 1.0e2) <= 1e-12
    # scalings: sigma grows as sqrt(volume) (N as its fourth root); weak dependence on T a
    f = m.trotter_state_dependent_steps
    assert f(3, 3, 10, 0.1, 0.5) == sd and f(3, 3, 50, 0.1, 0.5) == sd22
    s6 = f(6, 3, 50, 0.1, 0.5)
    assert s6["n_modes"] == 6 * 215 and 7.5 < s6["mean"] / sd["mean"] < 8.5             # the mean is extensive (8x the volume)
    assert 2.6 < s6["std"] / sd["std"] < 3.0                                            # sigma grows as sqrt(8) = 2.83
    assert [f(3, 3, 50, 0.1, t)["n_central"] for t in (0.15, 0.25, 1.0)] == [5008, 5008, 5342]
    assert f(3, 3, 200, 0.1, 0.5)["n_central"] == math.ceil(8 * math.sqrt(sd["lambda_var"] * 50 ** 3 / 0.1))
    b = dataclasses.replace(a, temp_lat_2033=c.Assumed(1.0, "test"))
    assert m.model(b, "2033").intermediates["n_trotter_2033"] == f(3, 3, 10, 0.1, 1.0)["n_central"] == 478
    assert [f(3, 3, 10, 0.1, t)["n_central"] for t in (0.15, 0.25, 1.0)] == [448, 448, 478]


def test_2033_trotter_worst_case_commutator_bound(a, r33):
    # r21 (3b), derived here, quoted since r22 as the worst case: N = sqrt(Lambda t^3 / eps), Lambda from operator norms
    i = r33.intermediates
    tb = i["trotter_bound"]
    assert i["t_over_a_2033"] == 10 and i["eps_trotter"] == 0.1 and i["g2_bare"] == 1.0
    assert math.isclose(tb["e_half_range"], E_HALF) and E_HALF == 15 / 16       # (g^2/2) x (15/4) / 2
    assert math.isclose(tb["b_norm"], B_NORM) and B_NORM == 20                   # hop 4 + 4 plaquettes x 4/g^2
    assert math.isclose(tb["lambda_link"], E_HALF ** 2 * 20 / 3 + E_HALF * 400 / 6)        # 5.86 + 62.5 = 68.36
    assert math.isclose(tb["lambda"], LAMBDA) and round(LAMBDA, 1) == 5537.1
    assert rel(i["trotter_lambda"], i["trotter_lambda_stated"]) <= 0.01          # '5.5e3 a^-3'
    assert tb["n_steps"] == N_BOUND_10 == 7442 == i["n_trotter_bound"]
    assert rel(i["n_trotter_bound"], i["n_trotter_bound_stated"]) <= 0.01        # '7.4e3 steps'
    assert math.isclose(N_BOUND_10 ** 2 * 0.1 / 10 ** 3, LAMBDA, rel_tol=1e-3)   # the bound is saturated at eps
    assert rel(i["trotter_dt_lattice_bound"], 1.34e-3) <= 0.01
    assert i["r22_bound"]["n_trot"] == N_BOUND == 83195                         # the r21-r23 count at 50a (record)
    # the order with the group-diagonal terms outside is the better one
    assert tb["lambda_link"] < tb["lambda_link_other_order"]
    # the direct function, and its scalings: N ~ t^{3/2}, eps^{-1/2}, sqrt(n_links)
    f = m.trotter_bound_steps
    assert f(81, 50, 0.1, 1.0, 3.75, 4.0, 4)["n_steps"] == N_BOUND and f(81, 10, 0.1, 1.0, 3.75, 4.0, 4) == tb
    assert math.isclose(f(81, 40, 0.1, 1.0, 3.75, 4.0, 4)["n_exact"], 8 * tb["n_exact"])
    assert math.isclose(f(81, 10, 0.4, 1.0, 3.75, 4.0, 4)["n_exact"], tb["n_exact"] / 2)
    assert math.isclose(f(324, 10, 0.1, 1.0, 3.75, 4.0, 4)["n_exact"], 2 * tb["n_exact"])
    # no bare coupling gives fewer than ~6.9e3 steps at t = 10a
    assert 1.5 < i["trotter_g2_at_min"] < 2.2 and rel(i["trotter_lambda_min"], 4.77e3) <= 0.01
    assert 6.9e3 < i["n_trotter_at_lambda_min"] < 6.95e3 <= N_BOUND_10
    for g2 in (0.3, 0.42, 1.0, 2.0, 4.0, 6.0):
        assert f(81, 10, 0.1, g2, 3.75, 4.0, 4)["n_steps"] >= i["n_trotter_at_lambda_min"]
    # the worst case is 16.5x the state-dependent count: 9 from the linear volume sum (81 against sqrt(81)), the rest
    # from the per-link size (the ratio does not depend on t)
    assert 16 < N_BOUND_10 / N_TROT < 17 and 16 < N_BOUND / N_TROT_R22 < 17 and 270 < LAMBDA / i["trotter_lambda_sd"] < 276
    # the pre-r21 '~30 steps' would need Lambda = 7.2e-4 a^-3
    assert i["n_trotter_pre_r21"] == 30 and rel(i["lambda_needed_for_pre_r21_steps"], 7.2e-4) <= 0.01


def test_2033_trotter_rule_switch_and_records(a, r33):
    # the rule switch: the headline, the worst-case bound, the fixed step of Ch. 5 (2 fm / 0.01 fm = 200), the old ~30
    i = r33.intermediates
    assert i["n_trotter_by_rule"] == {"state": N_TROT, "bound": N_BOUND_10, "step": 200, "stated": 30}
    assert i["n_dov_per_step_by_rule"] == {"state": 19, "bound": 6, "step": 33, "stated": 163}
    assert i["n_dov_by_rule"] == {"state": 10063, "bound": 46165, "step": 8113, "stated": 6403}   # 17 calls: + 801 (C1 + 472)
    assert i["n_dov_per_step_worst_stated"] == 6                                                  # '6 applications each'
    by = i["t_per_shot_by_trotter_rule"]
    for rule, n_dov in i["n_dov_by_rule"].items():
        assert math.isclose(by[rule], price_sp(n_dov)["t_shot"]), rule
        b = dataclasses.replace(a, trotter_rule_2033=c.Cited(rule, "test"))
        rb = m.model(b, "2033")
        assert rb.intermediates["n_dov_total"] == n_dov and math.isclose(rb.hard_ops[0], by[rule])
        assert math.isclose(rb.breakdown_total(), rb.hard_ops[0])
    assert rel(by["bound"], 8.54606e10) <= 1e-5 and rel(by["step"], 1.49894e10) <= 1e-5        # 17 calls (C1 8.40e10, 1.35e10; r25 8.31e10, 1.26e10)
    assert rel(by["stated"], 1.18269e10) <= 1e-5
    assert rel(i["t_per_shot_bound"], i["t_per_shot_bound_stated"]) <= 0.01                       # '8.4e10 T'
    assert rel(i["t_per_shot_empirical_step"], i["t_per_shot_empirical_step_stated"]) <= 0.04     # '1.4e10 T' (1.3508e10)
    assert i["t_per_shot_at_pre_r21_steps"] == by["stated"]
    # r22-r23 records at 50a and kappa 1e2: the worst case 1.276e12 and the headline 2.008e11
    assert rel(i["r22_bound"]["t_per_shot"], 1.275773e12) <= 1e-5 and i["r22_bound"]["per_step"] == 5
    # the range of the step count: 54 and 11 applications per step
    assert i["n_dov_per_step_range"] == (54, 11) and i["n_dov_range"] == (7291, 15527)       # 17 calls: + 1273 (C1 6490, 14726; r25 6018, 14254)
    lo, hi = i["t_per_shot_range"]
    assert math.isclose(lo, price_sp(7291)["t_shot"]) and math.isclose(hi, price_sp(15527)["t_shot"])
    assert rel(lo, 1.346909e10) <= 1e-5 and rel(hi, 2.87084e10) <= 1e-5                           # 17 calls (C1 1.1988e10-2.7226e10; r25 1.1115e10-2.6352e10)
    slo, shi = i["t_per_shot_range_stated"]
    assert rel(lo, slo) <= 0.04 and rel(hi, shi) <= 0.02                                          # '1.3-2.9e10 T' (1.347 -> 1.3)
    # a 12x spread in the step count is a 2.3x spread in T: the floor 2Q t
    assert 11.9 < 1274 / 107 < 12 and 2.1 < hi / lo < 2.2   # 2.13 (C1 2.27)
    assert i["n_dov_floor"] == 432 * 10 == 4320 and rel(4320, i["n_dov_floor_stated"]) <= 0.01   # '4.3e3 applications'
    assert math.isclose(i["t_per_shot_floor"], 4320 * T_DOV) and rel(i["t_per_shot_floor"], 8.0e9) <= 0.01
    assert rel(i["floor_ratio_to_envelope"], 7.98) <= 0.001                                       # 'still 8.0 times the limit'
    assert all(v > i["t_per_shot_floor"] for v in by.values())
    with pytest.raises(ValueError):
        dataclasses.replace(a, trotter_rule_2033=c.Cited("empirical", "test"))


def test_2033_jacobi_anger_applications_per_step(a, r33):
    # r22, derived here: the kernel is evolved by QSP on its block encoding (normalization 2Q = 432); the degree of
    # the Jacobi-Anger polynomial at x = 2Q dt for a QSP error of 1e-2 over the shot
    i = r33.intermediates
    assert a.eps_qsp.value == 1e-2 == i["eps_qsp"]
    assert math.isclose(i["jacobi_anger_x"], 432 * 10 / 450) and round(i["jacobi_anger_x"], 1) == 9.6     # '2Q dt = 9.6'
    assert i["n_dov_per_step"] == PER_STEP == 19 == i["n_dov_per_step_stated"]
    assert [m.dov_per_step(n, 10.0, 216, 1e-2) for n in (7442, 200, 30)] == [6, 33, 163]
    assert [m.dov_per_step(n, 10.0, 216, 1e-2) for n in (107, 450, 1274)] == [54, 19, 11]
    # r22 records at 50a
    assert [m.dov_per_step(n, 50.0, 216, 1e-2) for n in (83195, 1000, 30)] == [5, 35, 752]
    assert [m.dov_per_step(n, 50.0, 216, 1e-2) for n in (1186, 5031, 14242)] == [31, 13, 8]
    assert i["n_dov_per_step_bound_stated"] == 5
    assert rel(432 * 50 / 83195, 0.26) <= 0.01 and rel(432 * 0.05, 21.6) <= 1e-12
    # the Bessel recurrence against the power series, and the definition of the order
    for x in (0.26, 4.293, 9.6, 21.6):
        J = m.bessel_j_list(x, 60)
        for k in (0, 1, 5, 20):
            ser = sum((-1) ** s / (math.factorial(s) * math.factorial(s + k)) * (x / 2) ** (2 * s + k) for s in range(80))
            assert math.isclose(J[k], ser, rel_tol=1e-6, abs_tol=1e-12), (x, k)    # the series cancels at x = 21.6
        for eps in (1e-2 / 5031, 1e-2 / 450, 1e-6):
            R = m.jacobi_anger_order(x, eps)
            tail = lambda r: 2 * sum(abs(v) for v in m.bessel_j_list(x, 200)[r + 1:])   # noqa: E731
            assert tail(R) <= eps < tail(R - 1)
    # the truncated Jacobi-Anger series is e^{-i x cos(theta)} to the stated error
    for x, R, n in ((4.293, 13, 5031), (9.6, 19, 450)):
        J = m.bessel_j_list(x, R)
        for th in (0.0, 0.7, 1.9, math.pi):
            approx = J[0] + 2 * sum((-1j) ** k * J[k] * math.cos(k * th) for k in range(1, R + 1))
            exact = complex(math.cos(x * math.cos(th)), -math.sin(x * math.cos(th)))
            assert abs(approx - exact) <= 1e-2 / n
    # without a self-inverse block encoding the count doubles
    assert math.isclose(i["t_per_shot_if_not_self_inverse"], price_sp(450 * 38 + TPQ)["t_shot"])
    assert rel(i["t_per_shot_if_not_self_inverse"], 3.442e10) <= 0.01                             # 17 calls (C1 3.29e10; r25 3.21e10, r24 1.63e10)


def _herm_exp(H, t):
    import numpy as np
    w, v = np.linalg.eigh(H)
    return (v * np.exp(-1j * w * t)) @ v.conj().T


def test_r22_state_dependent_bound_numerically():
    """The estimate behind trotter_state_dependent_steps, on seeded random two-term Hamiltonians with thermal states.
    S = e^{-iD dt/2} e^{-iE dt} e^{-iD dt/2}, o1 = [E,[E,D]], o2 = [D,[D,E]]. (a) <o1> = <o2> in any state diagonal in
    H = E + D (o1 - o2 = [H,[E,D]]); (b) the mean of K = o1/12 - o2/24 is M/24; (c) the phase-removed thermal error is
    at most (t^3/N^2)(sigma(o1)/12 + sigma(o2)/24) to leading order, and the bound with the mean kept covers the error
    with the phase kept."""
    import numpy as np
    rng = np.random.default_rng(7)

    def herm(n):
        M = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
        return (M + M.conj().T) / 2

    def comm(X, Y):
        return X @ Y - Y @ X
    n_ok = 0
    for trial in range(24):
        d = 10
        t = [3.0, 30.0][trial % 2]
        n = [200, 800][(trial // 2) % 2]
        beta = [0.1, 1.0, 3.0][(trial // 4) % 3]
        E, D = herm(d), herm(d)
        w, v = np.linalg.eigh(E + D)
        p = np.exp(-beta * (w - w[0]))
        p /= p.sum()
        rho = (v * p) @ v.conj().T
        dt = t / n
        U = (v * np.exp(-1j * w * t)) @ v.conj().T
        S = _herm_exp(D, dt / 2) @ _herm_exp(E, dt) @ _herm_exp(D, dt / 2)
        ov = np.trace(rho @ U.conj().T @ np.linalg.matrix_power(S, n))
        err_kept, err_removed = math.sqrt(max(0, 2 - 2 * ov.real)), math.sqrt(max(0, 2 - 2 * abs(ov)))
        o1, o2 = comm(E, comm(E, D)), comm(D, comm(D, E))
        m1, m2 = np.trace(rho @ o1).real, np.trace(rho @ o2).real
        assert abs(m1 - m2) <= 1e-9 * (1 + abs(m1))                                   # (a)
        K = o1 / 12 - o2 / 24
        assert math.isclose(np.trace(rho @ K).real, m1 / 24, rel_tol=1e-9, abs_tol=1e-10)   # (b)
        r1, r2 = np.trace(rho @ o1 @ o1).real, np.trace(rho @ o2 @ o2).real
        s1, s2 = math.sqrt(max(r1 - m1 ** 2, 0)), math.sqrt(max(r2 - m2 ** 2, 0))
        pref = t ** 3 / n ** 2
        assert err_kept <= pref * (math.sqrt(r1) / 12 + math.sqrt(r2) / 24)           # the draft's second moment
        n_ok += err_removed <= pref * (s1 / 12 + s2 / 24)                             # (c)
    assert n_ok == 24


def test_r22_oscillator_commutators_and_exact_thermal_error():
    """The weak-coupling evaluation: one oscillator H = w (P^2 + X^2) / 2 with E = w P^2 / 2 (middle), D = w X^2 / 2
    (outer). o1 = -w^3 P^2, o2 = -w^3 X^2; thermal mean -(w^3/2) coth(w/2T), variance 2 [(w^3/2) coth(w/2T)]^2. And the
    exact thermal Trotter error of the 156 free modes, err^2 = 2 - 2 |prod_k Tr(rho_k U_k^dag S_k^N)|: 0.081 at the
    fixed step of Ch. 5 (1000 steps), 0.003 at the headline step count, so the estimate is a bound, loose by about 5
    in N in the one model where the exact error is known."""
    import numpy as np

    def fock(nmax):
        aa = np.diag(np.sqrt(np.arange(1, nmax)), 1)
        return (aa + aa.T) / math.sqrt(2), (aa - aa.T) / (1j * math.sqrt(2))

    def comm(X, Y):
        return X @ Y - Y @ X
    w, T, nmax, k = 1.7, 0.5, 120, 60
    X, P = fock(nmax)
    E, D = w * P @ P / 2, w * X @ X / 2
    o1, o2 = comm(E, comm(E, D)), comm(D, comm(D, E))
    assert abs((o1 + w ** 3 * P @ P)[:k, :k]).max() < 1e-8 and abs((o2 + w ** 3 * X @ X)[:k, :k]).max() < 1e-8
    p = np.exp(-w * np.arange(nmax) / T)
    p /= p.sum()
    ch = 1 / math.tanh(w / (2 * T))
    for o in (o1, o2):
        mean = (p * np.diag(o).real).sum()
        second = (p * np.diag(o.conj().T @ o).real).sum()
        assert math.isclose(mean, -(w ** 3 / 2) * ch, rel_tol=1e-9)
        assert math.isclose(second - mean ** 2, 2 * (w ** 3 / 2 * ch) ** 2, rel_tol=1e-8)

    def f_mode(wk, Ta, n, t=50.0, nm=80):
        Xm, Pm = fock(nm)
        dt = t / n
        Em, Dm = wk * (Pm @ Pm).real / 2, wk * (Xm @ Xm).real / 2
        S = _herm_exp(Dm, dt / 2) @ _herm_exp(Em, dt) @ _herm_exp(Dm, dt / 2)
        Sn_diag = np.diag(np.linalg.matrix_power(S, n))
        en = wk * (np.arange(nm) + 0.5)
        pk = np.exp(-(en - en[0]) / Ta)
        pk /= pk.sum()
        return (pk * np.exp(1j * en * t) * Sn_diag).sum()

    def exact_err(n, Ta=0.5, t=50.0):
        F = 1.0 + 0j
        for wk, mult in ((math.sqrt(3), 6), (math.sqrt(6), 12), (3.0, 8)):
            F *= f_mode(wk, Ta, n, t=t) ** (mult * 6)
        return math.sqrt(max(0.0, 2 - 2 * abs(F)))
    e1000, e5031 = exact_err(1000), exact_err(5031)
    assert rel(e1000, 0.0808) <= 0.02 and e1000 < 0.1               # 'the exact error of the free-field thermal state is 0.08'
    assert rel(e5031, 0.0032) <= 0.05
    assert 0.05 < exact_err(1186) < 0.065                           # the low end of the range also meets eps = 0.1
    # r24, t = 10a: Ch. 5's fixed step (200 steps) 0.022, the low end (107) 0.078, the headline (450) 0.0044
    e200, e107, e450 = (exact_err(n, t=10.0) for n in (200, 107, 450))
    assert rel(e200, 0.0225) <= 0.02 and rel(e107, 0.078) <= 0.02 and rel(e450, 0.0044) <= 0.05 and e107 < 0.1


def test_r21_second_order_bound_constants_numerically():
    """The bound behind trotter_bound_steps: for S(dt) = e^{-iA dt/2} e^{-iB dt} e^{-iA dt/2},
    ||S - e^{-i(A+B)dt}|| <= dt^3 (||[B,[B,A]]|| / 12 + ||[A,[A,B]]|| / 24) (childs2021theory). With the 1/12 and 1/24
    exchanged it fails when the inner term dominates, so the model's assignment (D outside, E inside) matters."""
    import numpy as np
    rng = np.random.default_rng(7)

    def herm(n):
        M = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
        return (M + M.conj().T) / 2

    def comm(X, Y):
        return X @ Y - Y @ X

    def nrm(M):
        return np.linalg.norm(M, 2)
    dt, swapped_fails = 0.05, 0
    for scale in (1.0, 0.05):
        for _ in range(40):
            A, B = scale * herm(5), herm(5)
            err = nrm(_herm_exp(A, dt / 2) @ _herm_exp(B, dt) @ _herm_exp(A, dt / 2) - _herm_exp(A + B, dt))
            bba, aab = nrm(comm(B, comm(B, A))), nrm(comm(A, comm(A, B)))
            assert err <= dt ** 3 * (bba / 12 + aab / 24) * (1 + 1e-9)
            swapped_fails += err > dt ** 3 * (bba / 24 + aab / 12)
    assert swapped_fails >= 40
    # ||[X,[Y,Z]]|| <= 4 ||X|| ||Y|| ||Z||
    for _ in range(20):
        X, Y, Zm = herm(4), herm(4), herm(4)
        assert nrm(comm(X, comm(Y, Zm))) <= 4 * nrm(X) * nrm(Y) * nrm(Zm) * (1 + 1e-12)


def test_r21_norm_inequality_is_loose_on_one_link():
    """r21 verifier: how tight is 4||X|| ||Y|| ||Z|| for the plaquette terms on one link? Class-function sector of
    SU(2) cut at j <= 3/2 (the chapter's Casimir cut): the Casimir C = diag(j(j+1)) and multiplication by the
    fundamental character chi_{1/2}, which hops j -> j +- 1/2. With the link's four plaquettes aligned,
    D = (8/g^2) chi_{1/2} and E = (g^2/2) C. The exact double commutators are 8.4 times below the inequality the
    worst-case Lambda uses: one of the three worst-case choices that stack in the 8.3e4 steps."""
    import numpy as np
    C = np.diag([j / 2 * (j / 2 + 1) for j in range(4)])                       # 0, 3/4, 2, 15/4
    T = np.diag(np.ones(3), 1) + np.diag(np.ones(3), -1)

    def comm(X, Y):
        return X @ Y - Y @ X

    def nrm(M):
        return np.linalg.norm(M, 2)
    g2 = 1.0
    E, D = (g2 / 2) * C, (8 / g2) * T
    exact = nrm(comm(E, comm(E, D))) / 12 + nrm(comm(D, comm(D, E))) / 24
    e, b_plaq = (g2 / 2) * 3.75 / 2, 4 * 4 / g2                                # the model's half-ranges, plaquettes only
    bound = e * e * b_plaq / 3 + e * b_plaq ** 2 / 6
    assert math.isclose(bound, 4.6875 + 40.0) and rel(exact, 5.31) <= 0.01
    assert exact <= bound and 8.0 < bound / exact < 8.8
    assert 2.8 < math.sqrt(bound / exact) < 3.0                                # at most ~3x in the step count


def test_r21_wilson_hop_rank_norm_and_strings():
    """Inputs of the worst-case Trotter bound (b) and of the Higgs string count: at r = 1 the hop matrix
    gamma^0 (1 - gamma_k)/2 is a rank-2 partial isometry (two hopping spinor components), the link-dressed hop has
    second-quantized norm N_c n_spin = 4, the 1+1D hop is a single bilinear, a Z3-dressed bilinear is 8 Pauli strings,
    and rho on four qubits is five Z terms."""
    import functools
    import itertools
    import numpy as np
    s0, sx = np.eye(2, dtype=complex), np.array([[0, 1], [1, 0]], complex)
    sy, sz = np.array([[0, -1j], [1j, 0]]), np.diag([1, -1]).astype(complex)
    g0 = np.kron(sz, s0)
    gam = [np.kron(sx, s) for s in (sx, sy, sz)]
    for G in gam:
        assert np.allclose(G @ G, np.eye(4)) and np.allclose(G @ g0 + g0 @ G, 0)
        M = g0 @ (np.eye(4) - G) / 2
        assert np.allclose(np.linalg.svd(M, compute_uv=False), [1, 1, 0, 0])            # n_spin = 2
    # one link, any U in SU(2): single-particle eigenvalues +-1 (x4), 0 (x8); second-quantized norm = 4
    U = _herm_exp(0.3 * sx - 1.1 * sy + 0.7 * sz, 1.0)
    M = np.kron(g0 @ (np.eye(4) - gam[0]) / 2, U)
    Hsp = np.block([[np.zeros((8, 8)), M], [M.conj().T, np.zeros((8, 8))]])
    ev = np.linalg.eigvalsh(Hsp)
    assert np.allclose(sorted(ev), [-1] * 4 + [0] * 8 + [1] * 4) and math.isclose(ev[ev > 0.5].sum(), 4.0)
    # 1+1D with gamma_1 diagonal: gamma^0 (1 - gamma_1)/2 has one non-zero entry (one hopping component, no frame)
    M1 = sx @ (np.eye(2) - sz) / 2
    assert np.count_nonzero(np.abs(M1) > 1e-12) == 1
    # a bilinear dressed by the Z3 phase omega^g on a two-qubit link register: 8 Pauli strings
    om = np.exp(2j * np.pi / 3)
    sm = np.array([[0, 1], [0, 0]], complex)
    a0, a1 = np.kron(sm, s0), np.kron(sz, sm)
    H = np.kron(a0.conj().T @ a1, np.diag([1, om, om.conjugate(), 1]))
    H = H + H.conj().T
    P = {"I": s0, "X": sx, "Y": sy, "Z": sz}
    n_str = sum(abs(np.trace(functools.reduce(np.kron, [P[x] for x in lab]) @ H)) / 16 > 1e-12
                for lab in itertools.product("IXYZ", repeat=4))
    assert n_str == 8 == m.RHOP_N_PAULI_DIAGONAL["Z3"]
    # rho on four qubits: a constant plus four Z terms
    rho = np.diag(np.arange(16)).astype(complex)
    n_rho = sum(abs(np.trace(functools.reduce(np.kron, [P[x] for x in lab]) @ rho)) / 16 > 1e-12
                for lab in itertools.product("IZ", repeat=4))
    assert n_rho == 5


# --- 2033: the shot ----------------------------------------------------------

def a_beta_rule_is_normalization(i):
    return i["beta_H_used"] == 2 * 216 / 0.5 == 864


def test_2033_dov_evolution_tpq_total(r33):
    i = r33.intermediates
    assert i["n_dov_evolution"] == 19 * N_TROT == 8550
    assert rel(i["n_dov_evolution"], i["n_dov_evolution_stated"]) <= 0.01     # '8.6e3 applications'
    assert math.isclose(i["t_evolution"], 8550 * T_DOV)
    assert rel(i["t_evolution"], i["t_evolution_stated"]) <= 0.02             # 1.5798e10 -> '1.6e10' (r24 8.0e9)
    assert i["delta_beta"] == 1e-4
    # referee C1 (2026-10-04): beta||H|| = 2Q / (T a) = 864, the normalization the filter's polynomial sees
    assert a_beta_rule_is_normalization(i)
    assert rel(i["d_beta_exact"], 89.2) <= 0.01 and i["d_beta"] == 89 == D_BETA and i["d_beta_record_r25"] == 30
    assert rel(i["d_beta_at_delta_beta_1e-3"], 77.3) <= 0.01
    assert i["n_dov_tpq"] == TPQ == 1513 and rel(TPQ, i["n_dov_tpq_stated"]) <= 0.01 and i["n_dov_tpq_r25"] == 240
    assert i["n_filter_calls"] == 17 and i["n_dov_tpq_one_call_per_round"] == TPQ_C1 == 712      # verifier: 2k + 1
    assert math.isclose(i["t_tpq"], TPQ * T_DOV)
    assert rel(i["t_tpq"], i["t_tpq_stated"]) <= 0.025                        # 2.7961e9 -> '2.8e9' (C1 1.3e9; r25 4.4e8)
    # the evolution, not the thermal state, is the shot
    assert not i["tpq_dominates"] and 0.84 < i["evolution_share_of_shot"] < 0.86 and 0.145 < i["tpq_share_of_shot"] < 0.155   # 15%
    assert rel(i["tpq_share_of_shot"], TPQ / NDOV) <= 1e-9 and round(i["tpq_share_of_shot"], 2) == 0.15   # '15%' (C1 7.7%)
    assert i["n_dov_total"] == NDOV == 10063 and i["n_queries_total"] == NQ == 4196271
    assert rel(NDOV, 1.0e4) <= 0.01                                           # 'N_Dov = 1.0e4 applications' (C1 9.3e3)
    assert math.isclose(i["t_per_shot"], NDOV * T_DOV) and rel(i["t_per_shot"], 1.8596712e10) <= 1e-6
    assert rel(i["t_per_shot"], i["t_per_shot_stated"]) <= 0.025              # 1.8597e10 -> '1.9e10' (C1 1.7e10; r25 1.6e10)
    assert rel(i["r25_headline"]["t_per_shot"], 1.6241697e10) <= 1e-6 and i["r25_headline"]["n_dov"] == 8790   # record
    assert math.isclose(r33.breakdown_total(), i["t_per_shot"])
    # not priced, recorded: the gauge terms of a step and the amplification reflections
    assert math.isclose(i["unpriced_gauge_T_per_link_step"],
                        2 * (216 + 42 * T33) + 2 * (0.5 * (350 + 4 * T33) + 3 * 112 + 6 * 392))
    assert math.isclose(i["unpriced_gauge_T_per_step"], 81 * i["unpriced_gauge_T_per_link_step"])
    assert 0.015 < i["unpriced_gauge_share_of_step"] < 0.025                  # 'about 2% of the step' (r24 4%)
    assert i["unpriced_fpaa_reflection_T"] == 16 * (2 * 972) * 7 == 217728      # 2n Toffolis (arxiv_2407_17966)
    assert rel(i["unpriced_fpaa_reflection_T"], 2e5) <= 0.1                   # 'about 2e5 T per shot'
    assert i["unpriced_fpaa_reflection_T"] < 3e-5 * i["t_per_shot"]


def test_2033_hard_ops_is_one_headline_and_responds_to_inputs(a, r33):
    i = r33.intermediates
    lo, hi = r33.hard_ops
    assert lo == hi == i["t_per_shot"] and round(lo, -5) == 1.85967e10
    # the record parameters of the second-quantized construction do not move the headline
    for b in (dataclasses.replace(a, hop_frame_undo=c.Assumed(False, "test")),
              dataclasses.replace(a, hop_share=c.Assumed("draft", "test")),
              dataclasses.replace(a, link_frame=c.Cited("draft", "test"))):
        assert m.model(b, "2033").hard_ops == r33.hard_ops and m.model(b, "2033").lq == r33.lq
    # the Toffoli convention moves everything except the U_x primitive
    b = dataclasses.replace(a, toffoli_convention=c.Cited("jones", "test"))
    tb = m.model(b, "2033")
    assert tb.hard_ops[0] < lo and tb.intermediates["t_umult"] == 392
    # the evolution time moves the step count (t^(3/2)) and the applications per step
    b = dataclasses.replace(a, t_evol_fm=c.Assumed(5, "test"))
    ib = m.model(b, "2033").intermediates
    assert ib["n_trotter_2033"] == math.ceil(math.sqrt(i["trotter_lambda_sd"] * 25 ** 3 / 0.1)) == 1779
    assert ib["n_dov_floor"] == 432 * 25 and ib["t_per_shot"] > 2 * lo
    # r24: the r22 headline is recovered with kappa 1e2 and t = 10 fm (and the Yukawa outside, r25)
    b = dataclasses.replace(a, t_evol_fm=c.Assumed(10, "test"), kappa_2033=c.Assumed(1e2, "test"),
                            yukawa_placement=c.Stated("outside", "test"), beta_H_rule=c.Assumed("stated", "test"),
                            aa_call_rule=c.Assumed("per_round", "test"))
    assert math.isclose(m.model(b, "2033").hard_ops[0], i["r22_headline"]["t_per_shot"], rel_tol=1e-12)
    # a larger QSP error budget lowers the order
    b = dataclasses.replace(a, eps_qsp=c.Assumed(0.1, "test"))
    assert m.model(b, "2033").intermediates["n_dov_per_step"] == m.dov_per_step(450, 10.0, 216, 0.1) < 19


def test_2033_hard_ops_matches_box(r33):
    assert c.close(r33.hard_ops, m.PUBLISHED["2033"].hard_ops, m.PUBLISHED["2033"].rel_tol)


def test_2033_no_rounding_inside_the_chain(a):
    b = dataclasses.replace(a, delta_sgn=c.Assumed(1e-6, "test"))
    i = m.model(b, "2033").intermediates
    assert i["d_sgn"] == math.ceil(30 * math.log(1e6)) == 415 and i["n_queries_per_dov"] == 2 * 415 + 1
    b = dataclasses.replace(a, kappa_2033=c.Assumed(60, "test"))
    t2 = m.model(b, "2033").intermediates["t_per_shot"]
    t1 = m.model(a, "2033").intermediates["t_per_shot"]
    assert math.isclose(t2, price_sp(NDOV, 415)["t_shot"]) and math.isclose(t1, price_sp(NDOV)["t_shot"])


def test_2033_envelope_ratio(r33):
    i = r33.intermediates
    assert math.isclose(i["ratio_to_envelope"], i["t_per_shot"] / 1e9)
    assert rel(i["ratio_to_envelope"], i["ratio_to_envelope_stated"]) <= 0.025     # 18.60 -> '19' (C1 17; r24 8.2)
    assert i["query_count_cut_needed"] == i["ratio_to_envelope"]
    assert not i["inside_2033_envelope"]


def test_2033_kappa_physical_scaling(r33):
    i = r33.intermediates
    lo, hi = i["t_per_shot_kappa_physical"]
    slo, shi = i["t_per_shot_kappa_physical_stated"]
    assert (slo, shi) == (6.1e11, 6.1e12)
    assert rel(lo, slo) <= 0.015 and rel(hi, shi) <= 0.015
    assert math.isclose(lo, price_sp(NDOV, 6908)["t_shot"]) and math.isclose(hi, price_sp(NDOV, 69078)["t_shot"])
    assert 32 < lo / i["t_per_shot"] < 33 and 320 < hi / i["t_per_shot"] < 330
    assert 6e2 <= lo / 1e9 < 7e2 and 6e3 <= hi / 1e9 < 7e3              # '6.1e2-6.1e3 times the budget'
    # r24: kappa ~ 1e2 at t = 10a (the r17-r23 kappa)
    assert math.isclose(i["t_per_shot_kappa_r22"], price_sp(NDOV, 691)["t_shot"])
    assert rel(i["t_per_shot_kappa_r22"], i["t_per_shot_kappa_r22_stated"]) <= 0.01     # '5.3e10' (r24 2.7e10)


def test_2033_eps_l_at_0p1_faults_per_shot(r33):
    i = r33.intermediates
    assert i["faults_per_shot"] == 0.1
    assert math.isclose(i["eps_l_implied"], 0.1 / i["t_per_shot"])
    assert rel(i["eps_l_implied"], i["eps_l_stated"]) <= 0.02           # 5.377e-12 -> '5.4e-12' (C1 5.8e-12; r25 6.2e-12)
    assert math.isclose(i["eps_l_condensate"], 0.1 / i["t_per_shot_condensate"])
    assert rel(i["eps_l_condensate"], i["eps_l_condensate_stated"]) <= 0.02   # 3.584e-11 -> '3.6e-11' (C1 7.6e-11; r25 2.3e-10)
    assert math.isclose(i["eps_l_short"], 0.1 / i["t_per_shot_short"])
    assert rel(i["eps_l_short"], i["eps_l_short_stated"]) <= 0.05             # 1.043e-11 -> '1.0e-11' (C1 1.2e-11; r25 1.4e-11)
    assert r33.epsilon_l == (i["eps_l_implied"], i["eps_l_implied"])
    assert i["eps_l_from_envelope"] == 0.1 / 1e9


def test_2033_thermal_route_is_tpq_and_the_sampler_is_repriced(r33):
    # r21 ruling (2): TPQ is the route; r22: the sampler estimate (r20) is re-priced at the single-particle application
    i = r33.intermediates
    assert i["thermal_route"] == "TPQ"
    assert rel(i["haar_rounds"], 3.36e7) <= 0.01
    assert rel(i["haar_rounds"], i["haar_rounds_stated"]) <= 0.15
    assert i["gibbs_jumps"] == 3000
    assert rel(i["gibbs_T_per_jump_implied"], 3.3e6) <= 0.02                # record: the retired '1e10'
    g = i["gibbs_estimate"]
    w = math.sqrt(math.log(1e3))
    assert (g["low"]["dov_per_jump"], g["central"]["n_jumps"], g["high"]["n_jumps"]) == (864, 3000, 9000)   # C1 (r25 100)
    assert math.isclose(g["central"]["dov_per_jump"], 2 * w * 864) and math.isclose(g["high"]["dov_per_jump"], 4 * w * 864)
    for k in g:
        assert math.isclose(g[k]["t_total"], g[k]["n_jumps"] * (g[k]["dov_per_jump"] * T_DOV + 1e3))
    lo, mid, hi = i["gibbs_T_estimate"]
    assert rel(lo, i["gibbs_T_stated"][0]) <= 0.01 and rel(hi, i["gibbs_T_stated"][1]) <= 0.03    # '4.8e12-1.5e14'
    assert rel(mid, i["gibbs_T_central_stated"]) <= 0.015                                         # '2.5e13'
    assert rel(i["gibbs_over_tpq"], i["gibbs_over_tpq_stated"]) <= 0.01                          # '9.0e3' (C1 1.9e4)
    assert rel(g["central"]["dov_per_jump"], 4.54e3) <= 0.01
    assert rel(g["low"]["over_tpq"], 1.7e3) <= 0.02 and rel(g["high"]["over_tpq"], 5.4e4) <= 0.02   # C1 3.6e3, 1.1e5
    # the sampler replaces the TPQ preparation: it would multiply the shot even at its low end
    assert lo > 4 * i["t_per_shot"] and mid > 20 * i["t_per_shot"]


def test_2033_shots_and_wall_time(r33, r33_24):
    i = r33.intermediates
    assert rel(i["shots_condensate_per_coupling"], 1e3) <= 0.01
    assert rel(i["shots_sphaleron_per_coupling"], 1e4) <= 0.01
    assert i["shots"] == 1e4
    assert math.isclose(i["wall_s_per_shot"], i["t_per_shot"] * 1e-6 + 1e-4)
    assert rel(i["wall_s_per_shot"], i["wall_s_per_shot_stated"]) <= 0.025      # 18,597 s -> '1.9e4 s'
    assert rel(i["wall_yr_per_point"], 5.8925) <= 1e-3                          # one temperature at 10a, 10% (C1 5.4230; r25 5.1463)
    assert r33.wall_time_s == i["wall_campaign_s"]                           # 2026-10-05 method test: 2.68-3.03 yr (rate campaign 36-60 yr, r25 29.8-50.5)
    assert rel(r33.wall_time_s[0], 8.45391e7) <= 1e-3 and rel(r33.wall_time_s[1], 9.56995e7) <= 1e-3
    # the r24 prints (Yukawa outside): 8241 s, 2.6 yr per point, 1.3e2 yr for the retired 50-point scan
    j = r33_24.intermediates
    assert rel(j["wall_s_per_shot"], 8241.44) <= 1e-5 and rel(j["wall_yr_per_point"], 2.6114) <= 1e-3
    assert rel(j["wall_yr_per_point"], j["wall_yr_per_point_stated"]) <= 0.01   # 2.61 -> '2.6 yr' (record)
    assert rel(j["wall_yr_serial_scan"], j["wall_yr_serial_scan_stated"]) <= 0.01   # 130.6 -> '1.3e2 yr' (record)
    assert r33_24.lq == (1010, 1010) and rel(r33_24.hard_ops[0], 8.241440e9) <= 1e-6


def test_r24_condensate_arm_is_the_preparation_alone(a, a24, r33, r33_24):
    """r24 (1): <psibar psi> is static (app06 workflow step 5): measured on the prepared state with no real-time
    evolution, so its shot is the thermal preparation (712 applications since referee C1; r25 240) at that circuit's
    own R-TOL tolerance."""
    i = r33.intermediates
    cs = i["condensate_shot"]
    assert cs["n_trot"] == 0 and cs["n_dov_evol"] == 0 and cs["n_dov"] == i["n_dov_condensate"] == TPQ == 1513
    P = price_sp(TPQ)
    assert cs["n_rot_shot"] == P["n_rot"] == TPQ * 2509 and math.isclose(cs["t_rot"], P["t_rot"])
    assert math.isclose(i["t_per_shot_condensate"], P["t_shot"]) and rel(i["t_per_shot_condensate"], 2.79010e9) <= 1e-5
    assert rel(i["t_per_shot_condensate"], i["t_per_shot_condensate_stated"]) <= 0.02          # '2.8e9 T' (C1 1.3e9; r25 4.4e8)
    assert rel(i["r25_headline"]["t_per_shot_condensate"], 4.41661e8) <= 1e-5                  # record
    # its own (looser) tolerance makes it a little cheaper than the same applications inside the sphaleron shot
    assert i["t_per_shot_condensate"] < i["t_tpq"] and rel(i["t_per_shot_condensate"], i["t_tpq"]) < 0.01
    assert not i["condensate_inside_2033_envelope"]                                            # C1: 1.3x the 1e9 limit
    w = i["wall_s_per_shot_condensate"]
    assert math.isclose(w, i["t_per_shot_condensate"] * 1e-6 + 1e-4) and rel(w, 2790.10) <= 1e-4
    assert rel(w, i["wall_s_condensate_stated"]) <= 0.01                                       # '2.8e3 s' (C1 1.3e3)
    yr = 3.156e7
    assert math.isclose(i["wall_yr_condensate_point"], 1e3 * w / yr)
    assert math.isclose(i["wall_day_condensate_temperature"], 1e3 * w / 86400)
    assert rel(i["wall_day_condensate_temperature"], i["wall_day_condensate_point_stated"]) <= 0.02   # 32.29 -> '32 days' (C1 15)
    # the r24 prints (Yukawa outside): 2.2412e8 T, 224.1 s, 2.6 days per point, 0.36 yr for the retired scan
    j = r33_24.intermediates
    P24 = price_sp(240, n_sgn=1)
    assert math.isclose(j["t_per_shot_condensate"], P24["t_shot"]) and rel(j["t_per_shot_condensate"], 2.2412e8) <= 1e-4
    assert rel(j["wall_s_per_shot_condensate"], 224.1) <= 1e-3 and rel(j["wall_day_condensate_point"], 2.594) <= 1e-3
    assert math.isclose(j["wall_yr_condensate_scan"], 50 * j["wall_yr_condensate_point"])
    assert rel(j["wall_yr_condensate_scan"], j["wall_yr_condensate_scan_stated"]) <= 0.02     # 0.355 -> '0.36 yr'
    # the kappa it is priced at moves it (the application cost), the evolution time does not
    b = dataclasses.replace(a, t_evol_fm=c.Assumed(10, "test"))
    assert m.model(b, "2033").intermediates["t_per_shot_condensate"] == i["t_per_shot_condensate"]
    b = dataclasses.replace(a24, kappa_2033=c.Assumed(1e2, "test"))
    assert rel(m.model(b, "2033").intermediates["t_per_shot_condensate"], 7.3e8) <= 0.01      # the r22 '7.3e8'


def test_r24_single_machine_wall_times(r33_24):
    """One machine, serial (ruling 2026-10-02). r24 record (Yukawa outside): one sphaleron point fits the 5-year horizon,
    not one year. The r25 tiers are tested in test_r25_*."""
    i = r33_24.intermediates
    yr = 3.156e7
    w = i["wall_s_per_shot"]
    assert math.isclose(i["wall_yr_sphaleron_point"], 1e4 * w / yr) and i["wall_yr_sphaleron_point"] == i["wall_yr_per_point"]
    assert math.isclose(i["wall_yr_serial_scan"], 50 * i["wall_yr_sphaleron_point"])
    assert i["campaign_horizon_yr"] == 5
    assert i["fits_horizon"] == {"condensate_point": True, "condensate_scan": True, "sphaleron_point": True, "scan": False}
    assert not i["sphaleron_point_under_one_year"]
    red = i["reduction_needed"]
    assert math.isclose(red["scan"], i["wall_yr_serial_scan"] / 5) and rel(red["scan"], 26.11) <= 1e-3
    assert round(red["scan"]) == i["reduction_scan_stated"] == 26
    assert math.isclose(red["sphaleron_point_one_year"], i["wall_yr_sphaleron_point"])
    # the shot that fits the scan at the same shot count: 5 yr / 5e5 shots = 315.6 s
    assert math.isclose(i["t_per_shot_to_fit_scan"], (5 * yr / 5e5 - 1e-4) / 1e-6)
    assert rel(i["t_per_shot_to_fit_scan"], i["t_per_shot_to_fit_scan_stated"]) <= 0.02            # '3.2e8 T'
    assert rel(i["t_per_shot"] / i["t_per_shot_to_fit_scan"], red["scan"]) <= 1e-6
    # the r22-r23 record: 64 yr per point
    assert rel(i["r22_headline"]["wall_yr_per_point"], 63.62) <= 1e-3
    assert rel(i["r22_headline"]["t_per_shot"], 2.007892e11) <= 1e-6 and i["r22_headline"]["n_dov"] == NDOV_R22


def test_r24_levers_for_one_point_under_a_year(a, r33_24):
    """r24 record (Yukawa outside): the levers that would bring one sphaleron point under a year (not applied)."""
    i = r33_24.intermediates
    lev, w = i["horizon_levers"], i["wall_yr_levers"]
    # the low end of the step range: the range's own low end, priced through the chain
    assert math.isclose(lev["low_step"], i["t_per_shot"] / i["t_per_shot_range"][0])
    assert round(lev["low_step"], 1) == 1.5 and i["horizon_levers_stated"]["low_step"] == 1.4   # r24 1.46; C1 1.43
    # a 15% statistical target: shots sigma^2 / eps^2 = 1e2 / 0.0225 = 4444
    assert a.lever_eps_obs.value == 0.15 and rel(i["lever_shots_sphaleron"], 4444.4) <= 1e-4
    assert math.isclose(lev["eps15"], 2.25) and i["horizon_levers_stated"]["eps15"] == 2.25
    assert math.isclose(lev["low_step_and_eps15"], lev["low_step"] * lev["eps15"])
    st = i["wall_yr_levers_stated"]
    assert (round(w["low_step"], 1), round(w["eps15"], 1), round(w["low_step_and_eps15"], 1)) == \
        (st["low_step"], st["eps15"], st["low_step_and_eps15"]) == (1.8, 1.2, 0.8)
    # each alone leaves the point above a year; together it is under
    assert w["low_step"] > 1 and w["eps15"] > 1 and w["low_step_and_eps15"] < 1
    assert rel(w["eps15"], i["wall_yr_sphaleron_point"] / 2.25) <= 1e-9
    # 'Fewer steps than the low end would gain at most 1.4 more, because of the floor 2Qt'
    assert math.isclose(i["low_step_over_floor"], i["t_per_shot_range"][0] / i["t_per_shot_floor"])
    assert round(i["low_step_over_floor"], 1) == 1.4
    # record (r23): the link-select lines held through the copy, 78 of the query's 598 Toffolis
    q = m.single_particle_query_2O(3, 3, 6)
    assert i["lever_hold_toffoli_saved"] == q["toffoli_items"]["fixup_iteration"] == 78
    assert math.isclose(lev["hold"], i["t_per_dov"] / (i["t_per_dov"] - 209 * 78 * 7), rel_tol=1e-9)
    assert i["lever_hold_qubits"] == 81 - q["anc_iteration"] == 75


def test_r24_evolution_time_and_linear_regime(a, r33):
    """t = 2 fm = 10a at a = 0.2 fm; t T = 5 at T a = 0.5 (the chapter's claim of a marginal linear window)."""
    i = r33.intermediates
    assert a.t_evol_fm.value == 2 and a.a_fm.value == 0.2 and i["t_over_a_2033"] == 10
    assert i["t_over_a_2033"] * i["temp_lat"] == 5.0
    assert i["n_dov_floor"] == 432 * 10 == 4320 and rel(4320, i["n_dov_floor_stated"]) <= 0.01
    assert math.isclose(i["t_per_shot_floor"], 4320 * T_DOV) and round(i["floor_ratio_to_envelope"], 1) == 8.0


def test_r23_no_machine_count(r33):
    # ruling 2026-10-02: no multi-machine assumption anywhere in the model or the chapter
    assert not hasattr(m.Assumptions(), "n_machines")
    for k in ("n_machines_needed", "n_machines_stated", "wall_yr_parallel", "fits_campaign_horizon", "n_machines"):
        assert k not in r33.intermediates
    assert not any("machines" in n for n in r33.notes)
    tex = (Path(m.__file__).resolve().parents[2] / "applications" / "app06_chiral_gauge.tex").read_text()
    assert not re.search(r"machines|shot-parallel|machine-years|fold parallelism|serialized|embarrassingly", tex)
    assert tex.count("one machine") == 3          # the 2028 box, the 2033 box and the wall-time paragraph (r25)


# --- RECORD: the retired second-quantized construction (r17-r21) -----------------------------

def test_2033_second_quantized_construction_is_a_record(a, r33):
    # the r21 box: QSVT of sgn on a Jordan-Wigner block encoding of H_W, at the worst-case step with 5 D_ov per step
    i = r33.intermediates
    rec = i["second_quantized_record"]
    assert rec["valid"] is False and rec["per_step"] == 5
    assert rec["n_dov"] == NDOV_SQ == 416215 and rec["n_queries"] == NQ_SQ == i["sq_n_queries"]
    assert i["sq_n_rot_shot"] == N33_SQ and math.isclose(i["sq_t_per_rot"], T33_SQ) and round(T33_SQ, 3) == 31.848
    assert math.isclose(rec["c_be"], C_BE33) and round(rec["c_be"], 1) == 62896.2
    assert rel(i["sq_c_be"], i["sq_c_be_stated"]) <= 0.01                       # r21 'c_be ~ 6.3e4'
    assert math.isclose(rec["t_per_dov"], 691 * C_BE33) and rel(rec["t_per_dov"], i["sq_t_per_dov_stated"]) <= 0.02
    assert math.isclose(rec["t_per_shot"], NQ_SQ * C_BE33) and rel(rec["t_per_shot"], 1.8089e13) <= 1e-4
    assert rel(rec["t_per_shot"], i["sq_t_per_shot_stated"]) <= 0.01            # r21 box '~1.8e13'
    assert (rec["lq_ancilla"], rec["lq_total"]) == (33, 1005)
    by = rec["t_per_shot_by_trotter_rule"]
    assert math.isclose(by["bound"], rec["t_per_shot"])
    assert rel(by["step"], 2.274e11) <= 1e-3 and rel(by["stated"], 1.691e10) <= 1e-3      # the r21 comparison prints
    # at equal step count the valid construction is 14.2x cheaper per application and per shot
    assert rel(rec["over_single_particle_per_application"], 14.18) <= 1e-3
    assert math.isclose(rec["over_single_particle_per_application"], 691 * C_BE33 / price_sp(416215, 691, n_sgn=1)["t_dov"])
    assert math.isclose(rec["over_single_particle_per_shot"], rec["over_single_particle_per_application"])
    # sgn of the many-body operator is not the bilinear with kernel sgn(h): they agree on the one-particle sector only,
    # and the many-body operator has an exact zero mode (the empty state), so its condition number is undefined
    import functools
    import numpy as np
    rng = np.random.default_rng(1)
    Qn = 6
    I2, X, Y, Z = np.eye(2), np.array([[0, 1], [1, 0]]), np.array([[0, -1j], [1j, 0]]), np.diag([1.0, -1.0])
    ann = [functools.reduce(np.kron, [Z] * k + [(X + 1j * Y) / 2] + [I2] * (Qn - k - 1)) for k in range(Qn)]
    h = rng.normal(size=(Qn, Qn)) + 1j * rng.normal(size=(Qn, Qn))
    h = h + h.conj().T

    def sgnm(A):
        w, v = np.linalg.eigh(A)
        return (v * np.sign(np.where(abs(w) < 1e-10, 0, w))) @ v.conj().T

    def bilinear(kern):
        return sum(kern[p, q] * ann[p].conj().T @ ann[q] for p in range(Qn) for q in range(Qn))
    H = bilinear(h)
    assert abs(np.linalg.eigvalsh(H)).min() < 1e-10 < abs(np.linalg.eigvalsh(h)).min()
    diff = sgnm(H) - bilinear(sgnm(h))
    assert np.linalg.norm(diff, 2) > 1.0
    number = np.rint(np.real(np.diag(sum(x.conj().T @ x for x in ann)))).astype(int)
    one = np.where(number == 1)[0]
    assert np.linalg.norm(diff[np.ix_(one, one)], 2) < 1e-10
    assert set(np.rint(np.linalg.eigvalsh(sgnm(H))).astype(int)) == {-1, 0, 1}


def test_2033_record_lcu_terms_per_site(a, r33):
    # r21 (3c): the LCU of one query of the second-quantized construction, per site. Record, except the Higgs strings,
    # which the headline still uses (40 per site, once per application)
    i = r33.intermediates
    per = i["lcu_per_site"]
    assert per["hop"] == 3 * 2 * 2 * 1 * 2 == 24 and per["hop_per_link"] == 8
    assert per["onsite"] == 4 * 2 == 8
    assert per["higgs"] == 8 * (4 + 1) == 40 == i["higgs_strings_per_site"]
    assert per["total"] == LCU_SITE == 72 == i["lcu_terms_per_site"] == i["lcu_terms_per_site_stated"]
    assert m.lcu_terms_per_site(3, 2, 2, 1, 2, radial_terms=5) == per
    assert i["lcu_terms"] == 27 * 72 == 1944 and rel(i["lcu_terms"], i["lcu_terms_stated"]) <= 0.03
    assert i["lcu_toffolis"] == 1944 and i["lcu_toffolis_per_term"] == 1.0
    assert i["sq_c_be_toffoli_T"] == C_TOFF == 13608
    assert i["lcu_address_qubits"] == 11
    b = dataclasses.replace(a, L_2033=c.Assumed(4, "test"))
    assert m.model(b, "2033").intermediates["lcu_terms"] == 64 * 72
    # the pre-r21 absolute: 1.7e3 = 8Q = 1728 is the 4Q + 4Q selected strings of the lift, paid once per application
    assert i["lcu_terms_pre_r21"] == 1.7e3 and rel(8 * 216, 1.7e3) <= 0.02 and "lift" in a.lcu_terms_stated_pre_r21.note
    assert a.lcu_toffoli_per_term.prov is c.Provenance.STATED and a.lcu_strings_per_bilinear.value == 2
    import numpy as np
    X, Y, Z = np.array([[0, 1], [1, 0]], complex), np.array([[0, -1j], [1j, 0]]), np.diag([1, -1]).astype(complex)
    sm = np.array([[0, 1], [0, 0]], complex)
    a0, a1 = np.kron(sm, np.eye(2)), np.kron(Z, sm)
    hopm = a0.conj().T @ a1 + a1.conj().T @ a0
    assert np.allclose(hopm, (np.kron(X, X) + np.kron(Y, Y)) / 2)


def test_2033_record_link_term_multiplexer(a, r33):
    # r20-r21: the Fock-space 2O colour multiplexer of the second-quantized construction (3 Toffolis + 18 T)
    i = r33.intermediates
    assert a.link_frame.value == "multiplexer" == i["link_frame"] and "RECORD" in a.link_frame.note
    assert (i["mux_toffoli"], i["mux_t_direct"], i["mux_ancilla"]) == (3, 18, 1) and i["mux_t_report"] == 39
    assert i["mux_link_n_mux"] == 4 and (i["mux_link_toffoli"], i["mux_link_t_direct"]) == (12 + 40, 72 + 32)
    assert i["mux_link_t"] == LINK_MUX == 468
    assert i["sq_n_rot_per_query"] == 25
    assert math.isclose(i["sq_c_be_hop_T"], 81 * 468) and i["sq_c_be_hop_T"] == 37908
    assert i["sq_c_be_higgs_T"] == 10584 and math.isclose(i["sq_c_be_rot_T"], 25 * T33_SQ)
    with pytest.raises(ValueError):
        dataclasses.replace(a, link_frame=c.Cited("zero", "test"))


def test_2033_record_draft_frame(a, r33):
    # r17 (e) + r19 (E21): one 2O Wilson d=3 link per BE query of the second-quantized construction, frame applied and
    # undone, squish and parity held per link (share="link"), n_spin = 2, MBU, no phasing. A record since r22.
    i = r33.intermediates
    h = fermion_hop_counts("2O", 1, fermion="wilson", d=3, undo=True, share="link", mcx="mbu", n_spin=2,
                           phasing="none")
    assert (h["toffoli"], h["t_direct"], h["n_rot"]) == (384, 144, 704)
    hd = fermion_hop_counts("2O", 1, fermion="wilson", d=3, undo=True, share="draft", mcx="mbu", n_spin=2,
                            phasing="none")
    assert (hd["toffoli"], hd["t_direct"], hd["n_rot"]) == (1182, 144, 704)               # r17 sensitivity
    items = {n: (tf, td, nr) for n, tf, td, nr in i["hop_items"]}
    assert items["colour_squish_and_parity"] == (2 * 105 + 56, 0, 0) == (266, 0, 0)       # once per link, MBU squish 105
    assert items["colour_rotations"] == (0, 0, 4 * 2 * 60)                                    # 480
    assert items["hop_squish_and_flags"] == (2 * 39, 0, 0)                                    # MBU hop squish 39
    assert items["hop_diagonalizers"] == (0, 2 * 56, 4 * 56)                                  # 2 x 7 x 2 x 2 = 56
    assert items["spinor_frame"] == (2 * 4 * 5, 2 * 16, 0)                                    # N x 4 W x (5 Tof, 4 T)
    assert i["hop_clifford_t_per_link"] == LINK_CT == 2832 and i["hop_n_rot_per_link"] == LINK_ROT
    assert math.isclose(i["hop_t_per_link"], LINK_CT + LINK_ROT * T33D) and round(i["hop_t_per_link"], 1) == 29768.7
    assert i["hop_undo"] is True and i["hop_n_spin"] == 2 and i["hop_mcx"] == "mbu" and i["hop_share"] == "link"
    assert a.hop_frame_undo.prov is c.Provenance.CITED and "E21 (1)" in a.hop_frame_undo.src
    assert "OUR ESTIMATE" in a.hop_frame_undo.note and "does not count it" in a.hop_frame_undo.note
    assert a.hop_share.prov is c.Provenance.CITED and a.hop_share.value == "link" and "E21 (1)" in a.hop_share.src
    assert a.hop_n_spin.prov is c.Provenance.CITED and "E21 (2)" in a.hop_n_spin.src
    hl = hop_link_cost("2O", 1, 3, eps=EPS33D, fermion="wilson", undo=True, n_spin=2, phasing="none", legacy=False)
    assert math.isclose(hl["t"], i["hop_t_per_link"])
    assert math.isclose(i["draft_over_mux_per_link"], i["hop_t_per_link"] / 468)
    assert rel(i["draft_over_mux_per_link"], i["draft_over_mux_per_link_stated"]) <= 0.01
    assert i["sq_draft_frame_n_rot_shot"] == NROTQ * NQ_SQ and math.isclose(i["sq_draft_frame_t_per_rot"], T33D)
    assert math.isclose(i["sq_c_be_draft_frame"], C_BE33D) and round(C_BE33D, 1) == 2436414.9
    assert rel(i["sq_c_be_draft_frame"], i["sq_c_be_draft_frame_stated"]) <= 0.02             # r21 '2.4e6'
    assert math.isclose(i["sq_t_per_shot_draft_frame"], NQ_SQ * C_BE33D)
    assert rel(i["sq_t_per_shot_draft_frame"], i["sq_t_per_shot_draft_frame_stated"]) <= 0.01  # r21 '7.0e14'
    assert rel(i["draft_over_mux_shot"], i["draft_over_mux_shot_stated"]) <= 0.01             # 38.7 -> '39 times'
    # the frame switch moves the record only
    b = dataclasses.replace(a, link_frame=c.Cited("draft", "test"))
    rb = m.model(b, "2033")
    assert rb.hard_ops == r33.hard_ops and rb.lq == (1015, 1015)
    assert math.isclose(rb.intermediates["second_quantized_record"]["t_per_shot"], NQ_SQ * C_BE33D)
    assert rb.intermediates["sq_lq_total"] == 1015


def test_2033_records_of_earlier_and_alternative_centres(r33):
    # every earlier print stays reproducible: the pre-r21 query count (390 x 691) and the 1.7e3-term LCU line
    i = r33.intermediates
    rec = i["t_per_shot_records"]
    assert i["n_queries_pre_r21"] == NQ_PRE == 269490
    assert rel(rec["vertex_1_umult_pre_r17"], 1.4792e10) <= 1e-4               # the old headline
    assert rel(rec["rhop_2_umult_per_link"], 2.335e10) <= 1e-3                  # R-HOP literal (2033-B)
    assert rel(rec["draft_w1_no_undo"], 4.342e11) <= 1e-3                       # the draft's W1 without the undo
    # r19 (E21 (1)): the r17 centre is the conservative sensitivity (squish recomputed at each application)
    assert math.isclose(rec["share_draft_r17_centre"], NQ_PRE * (11900 + 81 * LINK_CT_DRAFT + 27 * 392 + NROTQ * T_PRE))
    assert rel(rec["share_draft_r17_centre"], 6.891e11) <= 1e-4
    assert rel(rec["share_link_no_undo"], 3.935e11) <= 1e-3
    assert rel(rec["w1u_n_spin_1"], 3.120e11) <= 1e-3                           # r19-r20 centre at n_spin = 1
    assert rel(rec["w1u_mcx_2n_3"], 5.981e11) <= 1e-3
    # the r20 box: centre (draft frame), floor (multiplexer), top (the draft's W2), and the pre-r20 zero-T floor
    assert math.isclose(i["c_be_r20_centre_draft_frame"], C_BE_R20) and round(C_BE_R20, 1) == 2104718.8
    assert math.isclose(rec["r20_centre_draft_frame"], NQ_PRE * C_BE_R20) and rel(rec["r20_centre_draft_frame"], 5.672e11) <= 1e-4
    assert math.isclose(i["c_be_r20_floor_multiplexer"], C_FLOOR) and round(C_FLOOR, 1) == 61043.6
    assert rel(rec["r20_floor_multiplexer"], 1.645e10) <= 1e-3
    assert math.isclose(i["c_be_floor_zero_record"], C_FLOOR_ZERO) and round(C_FLOOR_ZERO, 1) == 12551.6
    assert round(rec["pre_r20_floor_zero_T"], -6) == 3.383e9
    # W2 = 4 C^S + 148 N, C^S(BO) = 9234.4 + 266.8 L; 266.8 / 1.15 = 232 rotations
    assert i["w2_rot_per_link"] == 4 * 232 == 928
    assert math.isclose(i["w2_t_per_link"], 4 * (9234.4 - 232 * 9.2) + 148 * 2) and round(i["w2_t_per_link"]) == 28696
    nq_top = 25 + 81 * 928
    t_top = fit(nq_top * NQ_PRE)                                                 # 32.707
    assert math.isclose(i["w2_t_per_rot"], t_top)
    assert math.isclose(i["c_be_r20_top_draft_w2"], 11900 + 81 * 28696 + 10584 + nq_top * t_top)   # 4,806,210.4
    assert round(rec["r20_top_draft_w2"], -8) == 1.2952e12
    # record (r19 prose, retired r20): the draft's frame without its rotation synthesis, not a realizable circuit
    assert math.isclose(i["c_be_rotation_free"], C_FLOOR_ZERO + 81 * 2832 + 10584)
    assert rel(i["ratio_centre_to_rotation_free"], 8.33) <= 0.01
    assert i["c_be_before_round_b"] == 1700 * 4 + 108 * 190 + 25 * 30 == 28070
    assert rel(i["t_per_shot_before_round_b"], 8e9) <= 0.06


def test_2033_umult_comes_from_groups_py(a, r33):
    um = GROUPS["2O"].primitives["U_mul"]
    i = r33.intermediates
    assert i["t_umult"] == um.t(EPS33) == um.t_const == 392
    assert i["t_umult_report_convention"] == 392
    assert i["umult_toffolis"] == um.toffoli == 56
    assert i["umult_clean_ancilla"] == um.ancilla == 4 and i["umult_ancilla_fit_in_box"]
    assert um.status is c.CircuitStatus.COMPILED and "arxiv_2312_10285" in um.src
    ga = next(p for p in r33.breakdown if p.name == "higgs_2O_group_action")
    assert ga.t_each == um.t_const and "arxiv_2312_10285" in ga.src and ga.count == 27 * NDOV    # once per application
    assert not any(f.name in ("t_umult", "umult_toffolis", "t_group_action") for f in dataclasses.fields(a))


# --- provenance honesty -----------------------------------------------------

def test_uncited_and_stated_inputs_are_tagged_as_such(a):
    assert not any(v.prov is c.Provenance.UNCITED
                   for v in (getattr(a, f.name) for f in dataclasses.fields(a)))
    assert a.delta_beta.prov is c.Provenance.ASSUMED
    assert a.n_umult_per_higgs_site.prov is c.Provenance.STATED
    assert a.spatial_hop_grouping.prov is c.Provenance.STATED
    assert a.n_trotter_2028.prov is c.Provenance.STATED and a.n_trotter_2028.value == 20
    assert a.aa_rounds.prov is c.Provenance.STATED
    assert a.t_per_step_2028_stated.prov is c.Provenance.STATED
    assert a.n_trotter_2033_pre_r21.prov is c.Provenance.STATED and "RECORD" in a.n_trotter_2033_pre_r21.note
    assert a.n_dov_per_step_stated.prov is c.Provenance.STATED and a.n_dov_per_step_bound_stated.prov is c.Provenance.STATED
    assert a.g2_bare.prov is c.Provenance.ASSUMED and a.casimir_max_2O.prov is c.Provenance.ASSUMED
    assert a.m0_wilson.prov is c.Provenance.ASSUMED and a.eps_qsp.prov is c.Provenance.ASSUMED
    assert a.gibbs_sampler_T.prov is c.Provenance.STATED
    assert "FermionPrimitives_unpub" in a.onsite_term_2028.src and "arxiv_2505_20419" in a.onsite_term_2028.src
    assert a.anc_rus.prov is c.Provenance.CITED and a.anc_rus.value == c.RUS_ANCILLA == 1 and "E26" in a.anc_rus.src
    assert a.fpaa_reflection_ancilla.prov is c.Provenance.CITED and "arxiv_2407_17966" in a.fpaa_reflection_ancilla.src
    # every record of the retired second-quantized construction says so
    for name in ("link_frame", "hop_frame_undo", "hop_share", "hop_phasing_be", "n_rot_be", "c_be_stated",
                 "c_be_draft_stated", "t_shot_draft_stated", "t_per_dov_sq_stated", "t_per_shot_sq_stated",
                 "anc_unary_control", "anc_link_pq", "anc_spinor", "lcu_terms_per_site_stated", "c_be_1p1d_sq_stated"):
        assert "RECORD" in getattr(a, name).note, name


def test_breakdown_statuses(a, r33):
    st = {p.name: p.status for p in r33.breakdown}
    assert set(st) == {"sp_query_toffoli", "sp_query_direct_T", "sp_query_rotations", "lift_toffoli", "lift_direct_T",
                       "lift_rotations", "outer_reflection_toffoli", "outer_reflection_rotation",
                       "higgs_2O_group_action", "higgs_select_toffoli"}
    # the query circuit is derived here and simulated at gate level (below): COMPILED; the rest are counts
    assert st["sp_query_toffoli"] is c.CircuitStatus.COMPILED and st["sp_query_direct_T"] is c.CircuitStatus.COMPILED
    assert all(v is c.CircuitStatus.SCALING for k, v in st.items() if not k.startswith("sp_query_")
               or k == "sp_query_rotations")
    assert all("DERIVED HERE" in p.src for p in r33.breakdown if p.name in ("sp_query_toffoli", "sp_query_direct_T"))
    assert not any("FermionPrimitives_unpub" in p.src for p in r33.breakdown)     # no unpublished input at 2033
    by = {p.name: p for p in r33.breakdown}
    assert by["sp_query_toffoli"].count == 598 * NQ and by["sp_query_rotations"].count == 6 * NQ
    assert by["lift_toffoli"].count == 862 * NDOV and by["higgs_select_toffoli"].count == 1080 * NDOV
    q = sum(p.t_total for p in r33.breakdown if p.name.startswith("sp_query_"))
    assert math.isclose(q / r33.breakdown_total(), NQD * C_SP / T_DOV)            # sgn is 99% of the shot (r25)
    rec = GROUPS["2O"].primitives["controlled_group_action"]
    assert rec.status is c.CircuitStatus.UNSOURCED and rec.t_const == 190
    assert "RETIRED 2026-09-28" in rec.note


# --- r20 E23 (1), (2): 2O, the Clifford group, and the derived multiplexer ---------------------------------

def _2o():
    """2O in the draft's fundamental irrep, ordered product g = (-1)^x1 j^x2 k^x3 u^(2x4+x5) t^x6 (arXiv:2312.10285)."""
    import itertools
    import numpy as np
    j = np.array([[0, 1], [-1, 0]], complex)
    k = np.diag([1j, -1j])
    u = 0.5 * np.array([[-1 - 1j, -1 + 1j], [1 + 1j, -1 + 1j]])
    t = np.array([[1, -1j], [-1j, 1]]) / np.sqrt(2)
    mp = np.linalg.matrix_power
    valid = [x for x in itertools.product((0, 1), repeat=6) if not (x[3] and x[4])]
    els = [(-1) ** x[0] * mp(j, x[1]) @ mp(k, x[2]) @ mp(u, 2 * x[3] + x[4]) @ mp(t, x[5]) for x in valid]
    return valid, els


def _is_pauli(M, basis):
    import numpy as np
    for P in basis:
        cc = np.trace(P.conj().T @ M) / M.shape[0]
        if abs(abs(cc) - 1) < 1e-9 and np.allclose(M, cc * P):
            return True
    return False


def _paulis(n):
    import functools
    import itertools
    import numpy as np
    one = [np.eye(2), np.array([[0, 1], [1, 0]]), np.array([[0, -1j], [1j, 0]]), np.diag([1, -1])]
    return [functools.reduce(np.kron, p).astype(complex) for p in itertools.product(one, repeat=n)]


def _jw_colour(U):
    import numpy as np
    M = np.eye(4, dtype=complex)
    M[1:3, 1:3] = U
    return M


def test_2O_is_the_clifford_group_up_to_phases(a, r33):
    import numpy as np
    valid, els = _2o()
    assert len({tuple(np.round(e, 8).ravel()) for e in els}) == 48
    p1 = _paulis(1)[1:]
    adj = set()
    for e in els:
        assert all(_is_pauli(e @ P @ e.conj().T, p1) for P in p1)          # normalizes the Pauli group
        adj.add(tuple(np.round([[np.trace(P @ e @ Q @ e.conj().T).real / 2 for Q in p1] for P in p1], 6).ravel()))
    assert len(adj) == 24                                                  # = single-qubit Clifford group mod phases
    # entries in Z[1/sqrt2, i]: 2 x entry is a + b sqrt2 (+ i(...)) with integers a, b
    for e in els:
        for z in (2 * e).ravel():
            for x in (z.real, z.imag):
                assert any(abs(x - b * math.sqrt(2) - round(x - b * math.sqrt(2))) < 1e-9 for b in range(-2, 3))
    assert a.clifford_2O.value == "clifford" and r33.intermediates["clifford_2O"] == "clifford"


def test_2O_colour_matchgate_is_clifford_only_on_Q8():
    valid, els = _2o()
    p2 = _paulis(2)
    clifford = [x for x, e in zip(valid, els) if all(_is_pauli(_jw_colour(e) @ P @ _jw_colour(e).conj().T, p2) for P in p2)]
    assert len(clifford) == 8 and all(x[3] == x[4] == x[5] == 0 for x in clifford)      # exactly Q8


def test_2O_eigenphases_leave_the_clifford_t_ring():
    # why the draft's diagonalizing route needs synthesis: eigenphases +-pi/3, +-2pi/3
    import numpy as np
    _, els = _2o()
    ph = {round(float(np.angle(v) / np.pi), 6) for e in els for v in np.linalg.eigvals(e)}
    assert {0.333333, 0.666667} <= ph


def test_2O_group_multiplication_is_not_affine():
    # U_x permutes basis states; a Clifford permutation is affine over GF(2). Find a violation.
    import itertools
    import numpy as np
    valid, els = _2o()
    key = {tuple(np.round(e, 8).ravel()): i for i, e in enumerate(els)}
    idx = {x: i for i, x in enumerate(valid)}

    def mul(i, j):
        return np.array(valid[key[tuple(np.round(els[i] @ els[j], 8).ravel())]])
    found = False
    for (a1, b1), (a2, b2), (a3, b3) in itertools.islice(itertools.combinations(itertools.product(range(8), range(8)), 3), 4000):
        xa = tuple(np.array(valid[a1]) ^ valid[a2] ^ valid[a3])
        xb = tuple(np.array(valid[b1]) ^ valid[b2] ^ valid[b3])
        if xa in idx and xb in idx and not np.array_equal(mul(a1, b1) ^ mul(a2, b2) ^ mul(a3, b3), mul(idx[xa], idx[xb])):
            found = True
            break
    assert found


def test_2O_colour_multiplexer_circuit_verified():
    """The explicit circuit behind colour_multiplexer_2O, simulated on x1..x6, q1, q2, f: equals
    sum_g |g><g| (x) M(U_g) on all 48 group states (f clean in and out)."""
    import functools
    import numpy as np
    valid, els = _2o()
    Q = ["x1", "x2", "x3", "x4", "x5", "x6", "q1", "q2", "f"]
    I2, X = np.eye(2, dtype=complex), np.array([[0, 1], [1, 0]], complex)
    Y, Z = np.array([[0, -1j], [1j, 0]]), np.diag([1, -1]).astype(complex)
    lab = {"I": I2, "X": X, "Y": Y, "Z": Z}
    Id = np.eye(2 ** len(Q), dtype=complex)

    def op(d):
        return functools.reduce(np.kron, [lab[d.get(q, "I")] for q in Q])

    def pr(q, b):
        return functools.reduce(np.kron, [np.diag([1 - b, b]).astype(complex) if qq == q else I2 for qq in Q])

    def cnot(cq, tq):
        return pr(cq, 0) + pr(cq, 1) @ op({tq: "X"})

    def tof(x, y, tq):
        return Id - pr(x, 1) @ pr(y, 1) + pr(x, 1) @ pr(y, 1) @ op({tq: "X"})

    def rot(th, d):
        return math.cos(th) * Id - 1j * math.sin(th) * op(d)
    gates, n_t, n_tof = [], 0, 0

    def ctrl_r(P, s=1):                      # exp(-i s pi/4 f P): two pi/8 rotations
        nonlocal n_t
        gates.extend([rot(s * math.pi / 8, {"q2": P}), rot(-s * math.pi / 8, {"q2": P, "f": "Z"})])
        n_t += 2

    def cc_phase(x, P):                      # exp(i pi/2 x q1 P): four pi/8 rotations
        nonlocal n_t
        th = -math.pi / 8
        gates.extend([rot(th, {"q2": P}), rot(-th, {"q2": P, x: "Z"}), rot(-th, {"q2": P, "q1": "Z"}),
                      rot(th, {"q2": P, x: "Z", "q1": "Z"})])
        n_t += 4
    gates.append(cnot("q2", "q1"))
    for x, body in (("x6", lambda: ctrl_r("X")),
                    ("x5", lambda: (ctrl_r("Y", -1), ctrl_r("X"), gates.append(op({"f": "Z"})))),
                    ("x4", lambda: (ctrl_r("Z", -1), ctrl_r("X"), gates.append(Id - 2 * pr("f", 1) @ pr("q2", 1)),
                                    gates.append(pr("f", 0) + 1j * pr("f", 1))))):
        gates.append(tof(x, "q1", "f"))
        n_tof += 1
        body()
        gates.append(tof(x, "q1", "f"))      # uncompute (by measurement in the costed circuit: 0 T)
    cc_phase("x3", "Z")
    cc_phase("x2", "Y")
    gates.append(Id - 2 * pr("x1", 1) @ pr("q1", 1))
    gates.append(cnot("q2", "q1"))
    W = functools.reduce(lambda A, B: B @ A, gates, Id)

    def idx(x, q2, q1):
        v = dict(zip(Q, list(x) + [q1, q2, 0]))
        i = 0
        for q in Q:
            i = 2 * i + v[q]
        return i
    basis = [(0, 0), (0, 1), (1, 0), (1, 1)]
    for x, e in zip(valid, els):
        B = np.array([[W[idx(x, *r), idx(x, *cc)] for cc in basis] for r in basis])
        assert np.allclose(B, _jw_colour(e)), x
    mx = m.colour_multiplexer_2O()
    assert (n_tof, n_t) == (mx["toffoli"], mx["t_direct"]) == (3, 18) and mx["n_rot"] == 0 and mx["ancilla"] == 1


def test_2O_frame_conjugation_gives_the_dressed_hop():
    # W^dag (sum_a psi_x,a^dag psi_y,a + h.c.) W = sum_ab psi_x,a^dag U_ab psi_y,b + h.c. for all 48 (JW, 4 modes)
    import functools
    import numpy as np
    _, els = _2o()
    Z, I2, sm = np.diag([1, -1]).astype(complex), np.eye(2, dtype=complex), np.array([[0, 1], [0, 0]], complex)

    def ann(i):
        return functools.reduce(np.kron, [Z] * i + [sm] + [I2] * (3 - i))
    ax, ay = [ann(0), ann(1)], [ann(2), ann(3)]
    H0 = sum(ax[cc].conj().T @ ay[cc] for cc in range(2))
    H0 = H0 + H0.conj().T
    for e in els:
        Hg = sum(ax[cc].conj().T @ sum(e[cc, d] * ay[d] for d in range(2)) for cc in range(2))
        Hg = Hg + Hg.conj().T
        Mg = np.eye(4, dtype=complex)
        ix = [2, 1]
        for b in range(2):
            for cc in range(2):
                Mg[ix[cc], ix[b]] = e[cc, b]
        Wf = np.kron(np.eye(4), Mg)
        assert np.allclose(Wf.conj().T @ H0 @ Wf, Hg)


def test_2033_ancilla_itemized_at_peak(a, r33):
    # E23 (5), r20; r22: the single-particle architecture. Peak = whole shot + one application + one query + link read
    i = r33.intermediates
    assert i["ancilla_stack"] == {"readout": 1, "qsp_signal": 1, "fpaa": 1, "h_lcu": 2}
    assert i["ancilla_application"] == {"lift_type_register": 2, "index_register": 9, "hov_lcu": 1, "sgn_signal": 1,
                                        "sgn_real_part": 1}
    assert i["ancilla_query"] == {"term_register": 5}
    assert i["ancilla_select"] == {"direction_flags": 3, "iteration_temporaries": 6, "link_copy": 5}
    assert i["lq_ancilla_itemized"] == 5 + 14 + 5 + 14 == 38 < i["lq_ancilla_asserted"] == 110
    assert i["lq_ancilla"] == 38 + 5                       # r25 (R10): the fan-out of the shared control
    # not concurrent with the link read: the lift's own iteration (9 temporaries + accumulator), the Higgs U_x, RUS,
    # the amplification reflection's two clean ancilla
    nc = i["ancilla_not_concurrent"]
    assert nc == {"lift_iteration": 10, "higgs_umult": 4, "rus": 1, "fpaa_reflection": 2}
    assert max(nc.values()) <= 5 + 14 and i["umult_clean_ancilla"] == 4
    # the readout control acts on the two inserted operators only and is idle during the TPQ filter
    assert "idle during the TPQ filter" in a.anc_readout.note and "arxiv_2208_13112" in a.anc_readout.note
    # record: the second-quantized construction (r21: 7 + 22 + 4 = 33; the draft frame's link term 14, 43)
    assert sum(i["sq_ancilla_stack"].values()) == 7
    assert i["sq_ancilla_be"] == {"lcu_address": 11, "unary_and_ladder": 10, "unary_control": 1}
    assert i["sq_ancilla_term"] == {"link_multiplexer": 1, "spinor_frame": 2, "rus": 1, "higgs_umult": 4}
    assert i["sq_lq_ancilla"] == 7 + 22 + 4 == 33
    assert i["sq_ancilla_link_draft_frame"] == {"pq_registers": 6, "hop_eigen_register": 3, "parity": 1,
                                                "largest_transient": 4}
    assert i["sq_ancilla_link_transients"] == {"mcx_ladder": 4, "rus": 1, "spinor_frame": 2}
    assert i["sq_lq_ancilla_draft_frame"] == 7 + 22 + 14 == 43


# --- r22 E26 (i): the single-particle architecture, checked circuit by circuit --------------------------------

def _rp(P, s=1):
    """exp(-i s pi/4 P)."""
    import numpy as np
    return (np.eye(2) - 1j * s * P) / math.sqrt(2)


def _2o_generator_decompositions():
    """u and u^2 as phase x Pauli x R(A, sa) R(B, sb), R(P, s) = exp(-i s pi/4 P): each controlled R is two pi/8
    Pauli rotations (2 T), the controlled phase x Pauli is Clifford."""
    import itertools
    import numpy as np
    I2, X = np.eye(2, dtype=complex), np.array([[0, 1], [1, 0]], complex)
    Y, Z = np.array([[0, -1j], [1j, 0]]), np.diag([1, -1]).astype(complex)
    Ps = {"I": I2, "X": X, "Y": Y, "Z": Z}
    u = 0.5 * np.array([[-1 - 1j, -1 + 1j], [1 + 1j, -1 + 1j]])
    out = {}
    for name, target in (("u", u), ("u2", u @ u)):
        for aa, bb in itertools.product("XYZ", repeat=2):
            for sa, sb in itertools.product((1, -1), repeat=2):
                for pn in "IXYZ":
                    for ph in (1, -1, 1j, -1j):
                        if name not in out and np.allclose(target, ph * Ps[pn] @ _rp(Ps[aa], sa) @ _rp(Ps[bb], sb)):
                            out[name] = (ph, pn, aa, sa, bb, sb)
    return out, Ps


def test_2O_index_register_multiplexer_is_10_T():
    """colour_multiplexer_index_2O: on ONE colour qubit, g = (-1)^x1 j^x2 k^x3 u^(2x4+x5) t^x6 is a product of singly
    controlled generators. x1: a sign; x2, x3: controlled iY, iZ (Clifford); x4, x5: two controlled pi/2 Pauli
    rotations each (2 T each); x6: one. 10 T, no Toffoli, no ancilla. Checked on all 48 elements, with the inverse."""
    import numpy as np
    valid, els = _2o()
    dec, Ps = _2o_generator_decompositions()
    assert set(dec) == {"u", "u2"}
    X, Y, Z = Ps["X"], Ps["Y"], Ps["Z"]
    assert np.allclose(np.array([[0, 1], [-1, 0]]), 1j * Y) and np.allclose(np.diag([1j, -1j]), 1j * Z)
    assert np.allclose(np.array([[1, -1j], [-1j, 1]]) / np.sqrt(2), _rp(X))
    for x, e in zip(valid, els):
        U, n_t = np.eye(2, dtype=complex), 0
        for bit, name in ((x[5], "t"), (x[4], "u"), (x[3], "u2")):
            if name == "t":
                U, n_t = (_rp(X) @ U if bit else U), n_t + 2
            else:
                ph, pn, aa, sa, bb, sb = dec[name]
                U, n_t = ((ph * Ps[pn] @ _rp(Ps[aa], sa) @ _rp(Ps[bb], sb) @ U) if bit else U), n_t + 4
        if x[2]:
            U = 1j * Z @ U
        if x[1]:
            U = 1j * Y @ U
        if x[0]:
            U = -U
        assert np.allclose(U, e), x
        assert n_t == 10
    mx = m.colour_multiplexer_index_2O()
    assert (mx["toffoli"], mx["t_direct"], mx["n_rot"], mx["ancilla"]) == (0, 10, 0, 0)
    # a controlled pi/2 Pauli rotation is two pi/8 rotations: exp(-i pi/4 n P) = exp(-i pi/8 P) exp(+i pi/8 Z P)
    n = np.diag([0, 1]).astype(complex)

    def ex(th, M):
        w, v = np.linalg.eigh(M)
        return (v * np.exp(-1j * th * w)) @ v.conj().T
    assert np.allclose(ex(math.pi / 4, np.kron(n, X)), ex(math.pi / 8, np.kron(np.eye(2), X)) @ ex(-math.pi / 8, np.kron(Z, X)))
    # controlled-(iY) is CY times S on the control: Clifford
    CiY = np.kron(np.diag([1, 0]), np.eye(2)) + np.kron(np.diag([0, 1]), 1j * Y)
    CY = np.kron(np.diag([1, 0]), np.eye(2)) + np.kron(np.diag([0, 1]), Y)
    assert np.allclose(CiY, np.kron(np.diag([1, 1j]), np.eye(2)) @ CY)


def test_lift_identity_and_mod3_shift():
    """The lift of arxiv_2607_28524 (Supplement theorem): with V^L_{i,a}, V^R_{i,a} in {P_a X_a, P_a Y_a} (phases on the
    type register), sum_i V^L_{i,a} V^R_{i,b} = 4 psi_a^dag psi_b, so <0| U_L h U_R |0> = psi^dag h psi / Q when the
    index register is prepared uniform over the Q valid modes (2^q otherwise). Also the mod-3 shift on two bits."""
    import functools
    import numpy as np
    rng = np.random.default_rng(1)
    Qn = 6
    I2, X, Y, Z = np.eye(2), np.array([[0, 1], [1, 0]]), np.array([[0, -1j], [1j, 0]]), np.diag([1.0, -1.0])
    PX = [functools.reduce(np.kron, [Z] * k + [X] + [I2] * (Qn - k - 1)) for k in range(Qn)]
    PY = [functools.reduce(np.kron, [Z] * k + [Y] + [I2] * (Qn - k - 1)) for k in range(Qn)]
    ann = [(PX[k] + 1j * PY[k]) / 2 for k in range(Qn)]
    VL = [[PX[k], PY[k], PX[k], -1j * PY[k]] for k in range(Qn)]
    VR = [[PX[k], PY[k], 1j * PY[k], PX[k]] for k in range(Qn)]
    for p, q in ((1, 4), (3, 3), (5, 0)):
        assert np.allclose(sum(VL[p][t] @ VR[q][t] for t in range(4)), 4 * ann[p].conj().T @ ann[q])
    h = rng.normal(size=(Qn, Qn)) + 1j * rng.normal(size=(Qn, Qn))
    h = h + h.conj().T
    H = sum(h[p, q] * ann[p].conj().T @ ann[q] for p in range(Qn) for q in range(Qn))
    # |0> -> uniform over (type, mode): amplitude 1/sqrt(4Q) on each side
    W = sum(h[p, q] * VL[p][t] @ VR[q][t] for p in range(Qn) for q in range(Qn) for t in range(4)) / (4 * Qn)
    assert np.allclose(W, H / Qn)
    # the strings are unitary, and there are 2Q distinct Majorana strings on each side (type X or Y, mode a)
    assert all(np.allclose(v @ v.conj().T, np.eye(2 ** Qn)) for row in VL + VR for v in row)
    assert m.lift_cost(Qn, 0)["toffoli"] == 2 * (2 * Qn - 1) and m.lift_cost(Qn, 0)["leaves"] == 2 * Qn
    # mod-3 increment on (b0, b1), valid values 0 = 00, 1 = 10, 2 = 01: SWAP, then X on b0 if b1 = 0 (2 Toffolis controlled)
    enc = {0: (0, 0), 1: (1, 0), 2: (0, 1)}
    for v, (b0, b1) in enc.items():
        b0, b1 = b1, b0
        if b1 == 0:
            b0 ^= 1
        assert (b0, b1) == enc[(v + 1) % 3]


class _Sim:
    """Sparse basis-state simulator for the SELECT check: counts every AND (Toffoli) and every T."""

    def __init__(self):
        self.q, self.state, self.tof, self.t, self.anc, self.peak = {}, {0: 1.0 + 0j}, 0, 0, set(), 0

    def reg(self, name):
        if name not in self.q:
            self.q[name] = len(self.q)
        return self.q[name]

    def bit(self, s, name):
        return (s >> self.q[name]) & 1

    def _on(self, s, ctrls, neg):
        return all(self.bit(s, x) for x in ctrls) and not any(self.bit(s, x) for x in neg)

    def x(self, tq, ctrls=(), neg=()):
        tb, new = 1 << self.reg(tq), {}
        for s, amp in self.state.items():
            s2 = s ^ tb if self._on(s, ctrls, neg) else s
            new[s2] = new.get(s2, 0) + amp
        self.state = new

    def phase(self, ph, ctrls=(), neg=()):
        for s in list(self.state):
            if self._on(s, ctrls, neg):
                self.state[s] *= ph

    def swap(self, p, q, ctrls=()):
        new = {}
        for s, amp in self.state.items():
            if self._on(s, ctrls, ()) and self.bit(s, p) != self.bit(s, q):
                s ^= (1 << self.q[p]) | (1 << self.q[q])
            new[s] = new.get(s, 0) + amp
        self.state = new

    def u1(self, U, tq, ctrls=()):
        tb, new = 1 << self.reg(tq), {}
        for s, amp in self.state.items():
            if self._on(s, ctrls, ()):
                v, s0 = (s >> self.q[tq]) & 1, s & ~tb
                for w in (0, 1):
                    if abs(U[w][v] * amp) > 1e-15:
                        new[s0 | (w * tb)] = new.get(s0 | (w * tb), 0) + U[w][v] * amp
            else:
                new[s] = new.get(s, 0) + amp
        self.state = {s: v for s, v in new.items() if abs(v) > 1e-14}

    def AND(self, p, q, out, negq=False):
        self.reg(out)
        self.anc.add(out)
        self.peak = max(self.peak, len(self.anc))
        self.x(out, ctrls=[p] if negq else [p, q], neg=[q] if negq else [])
        self.tof += 1

    def unAND(self, p, q, out, negq=False):                    # by measurement: 0 T
        self.x(out, ctrls=[p] if negq else [p, q], neg=[q] if negq else [])
        assert not any(self.bit(s, out) for s in self.state), out
        self.anc.discard(out)


def test_single_particle_select_reproduces_the_covariant_wilson_kernel():
    """Gate-level check of single_particle_query_2O on the full 3^3 lattice with a random 2O configuration (ported from
    the r22 verifier's simulation). (a) As matrices, the 13-term LCU with the priced schedule ([backward: A^dag] ->
    link -> [forward: A]) is the gauge-covariant Wilson kernel h_w / alpha, alpha = 2D + |D - m_0| = 7.5. (b) The
    circuit: on sampled (term, index) basis states SELECT returns the column of that term's unitary, every ancilla
    returns clean, and the count is input-independent: 594 Toffolis (598 with the 4-Toffoli projector) + 20 T, with
    14 SELECT ancilla at the link read."""
    import itertools
    import numpy as np
    rng = np.random.default_rng(11)
    valid, els = _2o()
    dec, Ps = _2o_generator_decompositions()
    I2, X, Y, Z = Ps["I"], Ps["X"], Ps["Y"], Ps["Z"]
    jm, km = 1j * Y, 1j * Z
    D = L = 3
    sites = list(itertools.product(range(L), repeat=D))
    links = {(k, x): valid[rng.integers(48)] for k in range(D) for x in sites}
    elem = dict(zip(valid, els))
    m0, r = 1.5, 1.0
    Q = 2 * 4 * L ** D

    def idx(colour, spin, x):
        return ((x[0] * L + x[1]) * L + x[2]) * 8 + spin * 2 + colour
    G0 = np.kron(I2, Z)                                        # spin index = s0 + 2 s1
    Gk = [np.kron(P, X) for P in (X, Y, Z)]

    def spinop(G):
        return np.kron(np.eye(L ** D), np.kron(G, np.eye(2)))

    def Tmat(k):
        T = np.zeros((Q, Q), complex)
        for x in sites:
            y = list(x)
            y[k] = (y[k] + 1) % L
            U = elem[links[(k, x)]]
            for sp in range(4):
                for c1 in range(2):
                    for c2 in range(2):
                        T[idx(c2, sp, tuple(y)), idx(c1, sp, x)] += U[c2, c1]
        return T
    Ts = [Tmat(k) for k in range(D)]
    hw = (D * r - m0) * spinop(G0)
    for k in range(D):
        hw = hw + 0.5j * spinop(Gk[k]) @ (Ts[k] - Ts[k].conj().T) - 0.5 * r * spinop(G0) @ (Ts[k] + Ts[k].conj().T)
    assert np.allclose(hw, hw.conj().T)
    alpha = D + D * r + abs(D * r - m0)
    assert alpha == 7.5 and abs(np.linalg.eigvalsh(hw)).max() <= alpha

    # (a) the LCU: term = (hop, direction, orientation, Dirac/Wilson); 12 hop unitaries + the on-site term
    def term_unitary(h, d, o, w):
        if not h:
            return spinop(G0)
        hopU = Ts[d].conj().T if o else Ts[d]
        return ((-1j if o else 1j) * spinop(Gk[d]) if w else -spinop(G0)) @ hopU
    terms = [(h, d, o, w) for h in (0, 1) for d in range(3) for o in (0, 1) for w in (0, 1)]
    wts = {t: ((abs(D * r - m0) / alpha) if t[0] == 0 else ((D + D * r) / alpha)) / 12 for t in terms}
    assert math.isclose(sum(wts.values()), 1.0)
    assert np.allclose(sum(wts[t] * term_unitary(*t) for t in terms), hw / alpha)
    assert len({t if t[0] else (0,) for t in terms}) == 4 * D + 1 == 13

    # (b) the circuit
    def enc(v):
        return v & 1, (v >> 1) & 1

    def build(sim, term, colour, spin, x):
        h, d, o, w = term
        s = 0
        for name, val in (("h", h), ("d0", d & 1), ("d1", d >> 1), ("o", o), ("w", w), ("c", colour),
                          ("s0", spin & 1), ("s1", spin >> 1)):
            s |= val << sim.reg(name)
        for i in range(D):
            b0, b1 = enc(x[i])
            s |= (b0 << sim.reg(f"x{i}0")) | (b1 << sim.reg(f"x{i}1"))
        sim.state = {s: 1.0 + 0j}

    def shift(sim, ctrl, i, up):                              # 2 Toffolis: controlled SWAP, doubly controlled X
        ops = [lambda: sim.swap(f"x{i}0", f"x{i}1", ctrls=[ctrl]), lambda: sim.x(f"x{i}0", ctrls=[ctrl], neg=[f"x{i}1"])]
        for op in (ops if up else ops[::-1]):
            op()
            sim.tof += 1

    def iterate(sim, root, k, leaf_fn):
        def rec(parent, level, prefix):
            if level == D:
                return leaf_fn(parent, tuple(prefix))
            b0, b1 = f"x{level}0", f"x{level}1"
            c1, c2, c0 = f"it{k}_{level}_a", f"it{k}_{level}_b", f"it{k}_{level}_z"
            sim.AND(parent, b0, c1)
            sim.AND(parent, b1, c2)
            sim.reg(c0)                                       # the third child by CNOTs, in place of the parent line
            for cc in (parent, c1, c2):
                sim.x(c0, ctrls=[cc])
            rec(c0, level + 1, prefix + [0])
            rec(c1, level + 1, prefix + [1])
            rec(c2, level + 1, prefix + [2])
            for cc in (c2, c1, parent):
                sim.x(c0, ctrls=[cc])
            sim.unAND(parent, b1, c2)
            sim.unAND(parent, b0, c1)
        rec(root, 0, [])

    def ctrl_gen(sim, ctrl, d_, dagger):
        ph, pn, aa, sa, bb, sb = d_
        if not dagger:
            seq = [_rp(Ps[bb], sb), _rp(Ps[aa], sa), ph * Ps[pn]]
        else:
            seq = [(ph * Ps[pn]).conj().T, _rp(Ps[aa], -sa), _rp(Ps[bb], -sb)]
        for U in seq:
            sim.u1(U, "c", ctrls=[ctrl])
        sim.t += 4

    def select(sim):
        sim.AND("h", "d0", "f1")
        sim.AND("h", "d1", "f2")
        sim.reg("f0")
        for cc in ("h", "f1", "f2"):
            sim.x("f0", ctrls=[cc])
        for k in range(D):                                    # backward: shift down before the read
            sim.AND(f"f{k}", "o", f"b{k}")
            shift(sim, f"b{k}", k, up=False)
            sim.unAND(f"f{k}", "o", f"b{k}")
        for jb in range(1, 6):
            sim.reg(f"r{jb}")
            sim.anc.add(f"r{jb}")

        def read(k, count):
            def fn(leaf, x):
                g = links[(k, x)]
                if count and g[0]:
                    sim.phase(-1, ctrls=[leaf])               # CZ(leaf, x1): Clifford
                for jb in range(1, 6):
                    if g[jb]:
                        sim.x(f"r{jb}", ctrls=[leaf])         # Toffoli(leaf, link bit -> scratch)
                    sim.tof += count
            return fn
        for k in range(D):
            iterate(sim, f"f{k}", k, read(k, 1))
        sim.peak = max(sim.peak, len(sim.anc) + 3 + 6)         # 3 direction flags and 2 lines x 3 levels are live here
        for jb, name in ((5, "t"), (4, "u"), (3, "u2"), (2, "k"), (1, "j")):       # M, forward
            f = f"mf{jb}"
            sim.AND(f"r{jb}", "o", f, negq=True)
            if name == "t":
                sim.u1(_rp(X), "c", ctrls=[f])
                sim.t += 2
            elif name in dec:
                ctrl_gen(sim, f, dec[name], False)
            else:
                sim.u1(km if name == "k" else jm, "c", ctrls=[f])
            sim.unAND(f"r{jb}", "o", f, negq=True)
        for jb, name in ((1, "j"), (2, "k"), (3, "u2"), (4, "u"), (5, "t")):       # M^dag, backward
            f = f"mb{jb}"
            sim.AND(f"r{jb}", "o", f)
            if name == "t":
                sim.u1(_rp(X, -1), "c", ctrls=[f])
                sim.t += 2
            elif name in dec:
                ctrl_gen(sim, f, dec[name], True)
            else:
                sim.u1((km if name == "k" else jm).conj().T, "c", ctrls=[f])
            sim.unAND(f"r{jb}", "o", f)
        for k in range(D):                                    # the copy is measured out; fix-ups walk the addresses
            iterate(sim, f"f{k}", k, read(k, 0))
        for jb in range(1, 6):
            assert not any(sim.bit(s, f"r{jb}") for s in sim.state)
            sim.anc.discard(f"r{jb}")
        for k in range(D):                                    # forward: shift up after the read (flags recomputed)
            sim.AND(f"f{k}", "o", f"w{k}", negq=True)
            shift(sim, f"w{k}", k, up=True)
            sim.unAND(f"f{k}", "o", f"w{k}", negq=True)
        for k in range(D):                                    # Dirac structure: Gamma^k = Gamma^0 gamma^k
            gq = f"g{k}"
            sim.AND(f"f{k}", "w", gq)
            sim.u1([X, Y, Z][k], "s1", ctrls=[gq])
            sim.u1(1j * Y, "s0", ctrls=[gq])
            sim.phase(1j, ctrls=[gq])
            sim.phase(-1, ctrls=[gq, "o"])
            sim.unAND(f"f{k}", "w", gq)
        sim.phase(-1, ctrls=["h"], neg=["w"])                 # Wilson hop: -1
        sim.u1(Z, "s0")                                       # Gamma^0 on every term
        for cc in ("f2", "f1", "h"):
            sim.x("f0", ctrls=[cc])
        assert not any(sim.bit(s, "f0") for s in sim.state)
        sim.unAND("h", "d1", "f2")
        sim.unAND("h", "d0", "f1")

    counts = set()
    samples = [(t, sites[rng.integers(27)], int(rng.integers(4)), int(rng.integers(2))) for t in terms for _ in range(4)]
    for term, x, sp, col in samples:
        sim = _Sim()
        build(sim, term, col, sp, x)
        select(sim)
        counts.add((sim.tof, sim.t))
        assert not sim.anc
        column = np.zeros(Q, complex)
        for s, amp in sim.state.items():
            assert [sim.bit(s, n) for n in ("h", "d0", "d1", "o", "w")] == [term[0], term[1] & 1, term[1] >> 1,
                                                                             term[2], term[3]]
            xo = tuple(sim.bit(s, f"x{i}0") + 2 * sim.bit(s, f"x{i}1") for i in range(D))
            column[idx(sim.bit(s, "c"), sim.bit(s, "s0") + 2 * sim.bit(s, "s1"), xo)] += amp
        assert np.allclose(column, term_unitary(*term)[:, idx(col, sp, x)]), term
    q = m.single_particle_query_2O(3, 3, 6)
    assert counts == {(q["toffoli"] - q["toffoli_items"]["projector"], q["t_direct_items"]["mux_and_inverse"])} == {(594, 20)}
    assert q["anc_direction_flags"] + q["anc_iteration"] + q["anc_link_copy"] == 14


def test_overlap_kernel_one_link_dependence_is_not_enhanced():
    """r22 verifier: the worst-case Trotter bound takes the Wilson hop's norm (4) for the kernel's dependence on one
    link. For the overlap kernel the dependence is that of sgn(h_w[U]), which is nonlocal. On seeded random 2O
    configurations the Fock-space half-range of psi^dag sgn(h_w) psi under a change of one link stays near the Wilson
    hop's (a one-link change is a rank-8 perturbation): no inverse-gap enhancement."""
    import itertools
    import numpy as np
    rng = np.random.default_rng(99)
    _, els = _2o()
    I2, X = np.eye(2, dtype=complex), np.array([[0, 1], [1, 0]], complex)
    Y, Z = np.array([[0, -1j], [1j, 0]]), np.diag([1, -1]).astype(complex)
    D = L = 3
    sites = list(itertools.product(range(L), repeat=D))
    sidx = {x: n for n, x in enumerate(sites)}
    G0, Gk = np.kron(I2, Z), [np.kron(P, X) for P in (X, Y, Z)]

    def hw_of(links):
        h = np.zeros((216, 216), complex)
        for x in sites:
            n = sidx[x]
            h[n * 8:(n + 1) * 8, n * 8:(n + 1) * 8] += (D - 1.5) * np.kron(G0, I2)
            for k in range(D):
                y = list(x)
                y[k] = (y[k] + 1) % L
                n2 = sidx[tuple(y)]
                hop = np.kron(0.5j * Gk[k] - 0.5 * G0, links[(k, x)])
                h[n2 * 8:(n2 + 1) * 8, n * 8:(n + 1) * 8] += hop
                h[n * 8:(n + 1) * 8, n2 * 8:(n2 + 1) * 8] += hop.conj().T
        return h

    def sgn(h):
        w, v = np.linalg.eigh(h)
        return (v * np.sign(w)) @ v.conj().T

    def fock_half_range(kern):
        w = np.linalg.eigvalsh(kern)
        return max(w[w > 0].sum(), -w[w < 0].sum()) / 2
    wil, ov = [], []
    for _ in range(2):
        links = {(k, x): els[rng.integers(48)] for k in range(D) for x in sites}
        l2 = dict(links)
        l2[(int(rng.integers(3)), sites[rng.integers(27)])] = els[rng.integers(48)]
        h, h2 = hw_of(links), hw_of(l2)
        wil.append(fock_half_range(h2 - h))
        ov.append(fock_half_range(sgn(h2) - sgn(h)))
    assert max(wil) <= 4.0 + 1e-9                      # the Wilson hop: at most N_c n_spin = 4
    assert max(ov) <= 4.5 and all(o <= 1.5 * w_ for o, w_ in zip(ov, wil))    # verifier, 100 samples: mean 2.8, max 4.5


# --- r25 (H. Lamm 2026-10-02: R1, R4, R9, R10 and the queued Ch. 9 rulings) ------------------------------------

def test_r25_yukawa_inside_doubles_the_sign_sequences(a, r33, r33_24):
    """(a) The Ginsparg-Wilson-projected Yukawa is inside: its projector (1 - sgn(H_W))/2 is a second degree-d_sgn QSVT
    polynomial on the same block encoding, once per application; the lift and the Higgs vertex stay once."""
    i, j = r33.intermediates, r33_24.intermediates
    assert a.yukawa_placement.value == "inside" and a.yukawa_placement.prov is c.Provenance.STATED
    assert i["n_queries_per_dov"] == 417 and j["n_queries_per_dov"] == 209
    assert i["c_lift"] != j["c_lift"] and i["c_higgs_per_dov"] == j["c_higgs_per_dov"] == C_HIGGS   # same circuits
    assert 1.96 < i["t_per_dov"] / j["t_per_dov"] < 1.98                     # 'twice the cost'
    assert i["lq_ancilla_itemized"] == j["lq_ancilla"] == 38                 # no new qubit: same signal and real-part
    with pytest.raises(ValueError):
        dataclasses.replace(a, yukawa_placement=c.Stated("both", "test"))


def test_r25_short_time_shot_at_its_own_depth(a, r33):
    """(R4) The 5a shot: its own state-dependent step count, applications per step and R-TOL tolerance."""
    i = r33.intermediates
    sd = i["trotter_state_dependent_short"]
    assert i["short_time_lat"] == 5.0 and sd["n_central"] == math.ceil(math.sqrt(sd["lambda_var"] * 125 / 0.1)) == 160
    sh = i["short_shot"]
    assert sh["per_step"] == m.dov_per_step(160, 5.0, 216, 1e-2) == 23 and sh["n_dov"] == 160 * 23 + TPQ == 5193
    assert math.isclose(sh["x"], 432 * 5 / 160) and sh["x"] == 13.5
    P = price_sp(5193)
    assert math.isclose(sh["t_shot"], P["t_shot"]) and rel(sh["t_shot"], 9.58966e9) <= 1e-5
    assert rel(sh["t_shot"], i["t_per_shot_short_stated"]) <= 0.01             # '9.6e9 T' (C1 8.1e9; r25 7.2e9)
    assert rel(i["wall_s_per_shot_short"], i["wall_s_short_stated"]) <= 0.01   # '9.6e3 s' (C1 8.1e3)
    assert sh["t_shot"] < 0.52 * i["t_per_shot"]   # the shared thermal filter (17 calls) keeps it above half


def test_r25_depth_matches_the_factory_analysis(a, r33, r33_24):
    """(R9, R10) T-depth per query and per application reproduce factory.json (Ch. 9 runs) at the r24 tolerance."""
    tr = 26.476
    assert m.sp_query_depth_2O(tr, "register") == (614 + 3 * tr, 622 + 6 * tr)
    assert m.sp_query_depth_2O(tr, "fanout") == (274 + 3 * tr, 274 + 3 * tr)
    assert round(m.sp_query_depth_2O(tr, "register")[0]) == 693 and round(m.sp_query_depth_2O(tr, "fanout")[0]) == 353
    lift = {"toffoli": 862, "t_direct": 12, "n_rot": 6}
    outer = {"toffoli": 18, "n_rot": 1}
    lo_f, _ = m.application_depth_2O(209, tr, "fanout", lift, outer, 27, 1080, 56, 4)
    lo_r, hi_r = m.application_depth_2O(209, tr, "register", lift, outer, 27, 1080, 56, 4)
    assert rel(lo_f, 76289) <= 2e-3 and rel(lo_r, 147349) <= 1e-3 and rel(hi_r, 166868) <= 1e-3
    # the r24 chain (Yukawa outside, as itemized): F* 5.6-6.4, the factory's 'NO' (printed 2.6 yr -> 4.1-4.6 yr)
    j = r33_24.intermediates
    fr = j["f_star_by_run"]["sph10"]
    assert 5.5 < fr[0] < 5.7 and 6.3 < fr[1] < 6.5
    # r25 headline: the fan-out holds F* >= 10 on every run
    i = r33.intermediates
    assert i["depth_variant"] == "fanout" and i["baseline_ok_all_runs"]
    assert all(12.1 < f[0] and f[1] < 12.6 for f in i["f_star_by_run"].values())
    assert round(i["f_star"][0]) == i["f_star_stated"] == 12                  # 'raises this to 12'
    fr = i["f_star_by_run_register"]["sph10"]
    assert rel(fr[0], i["f_star_register_stated"][0]) <= 0.01 and rel(fr[1], i["f_star_register_stated"][1]) <= 0.01
    lo, hi = i["slowdown_register"]
    assert rel(lo, 1.6) <= 0.03 and rel(hi, 1.8) <= 0.03                       # '1.6-1.8 times slower'
    assert i["umult_concurrent"] == 19 // 4 == 4


def test_r25_exports_2033(r33):
    """(R9) the seven exports, f_star = t / depth, floor and factories summed over the campaign runs."""
    i = r33.intermediates
    for k in c.DEPTH_EXPORTS:
        assert k in i, k
    t, d = c._band(i["t_per_shot"]), c._band(i["t_depth_per_shot"])
    assert i["f_star"] == pytest.approx((t[0] / d[1], t[1] / d[0]))
    runs = i["exports_by_run"]
    assert i["floor_wall_s"] == (sum(e["floor_wall_s"][0] for e in runs[0].values()),
                                 sum(e["floor_wall_s"][1] for e in runs[1].values()))
    assert i["floor_wall_s"][1] < i["wall_campaign_s"][1]                      # F* >= 10: the serial wall is the wall
    assert 26.5 < i["factories_for_1yr"][0] < 27.1 and 30.0 < i["factories_for_1yr"][1] < 30.6   # method test (rate campaign 360, 600)
    # the 10a test (1111 shots, 0.65 yr) and the condensate (32 days per temperature) fit a year; the 5a test does not
    assert runs[0]["sph10"]["fits_1yr"][1] and runs[0]["cond"]["fits_1yr"][1] and not runs[0]["sph5"]["fits_1yr"][1]


def test_r25_two_tiers_2033(a, r33):
    """Author ruling 2026-10-05 (method test, no rate at 3^3). First result: the growth of the full <dN_CS^2> at one
    symmetric-side temperature at 5a and 10a, 30% each, the 5a variance ratio derived from the free-field signals.
    Campaign: the test + the condensate at 4-8 temperatures at 10%. The r25 rate campaign is a record."""
    i = r33.intermediates
    tf, tc_ = i["tiers"]["first_result"], i["tiers"]["campaign"]
    assert a.eps_first.value == 0.3 and a.ncs_var_ratio_short.value == 4 and a.n_temperatures.value == (4, 8)
    bg = i["ncs_free_background"]
    vr = (bg[10.0] / bg[5.0]) ** 2
    assert i["method_test_var_ratio_5a"] == vr and 4.94 < vr < 4.96                  # chapter '2.2 times smaller ... 4.9 times'
    assert round(bg[10.0] / bg[5.0], 1) == 2.2 and round(vr, 1) == 4.9
    assert tf["shots"] == pytest.approx({"sph10": 1e2 / 0.09, "sph5": vr * 1e2 / 0.09})
    assert rel(tf["shots"]["sph10"], i["shots_first_stated"][0]) <= 0.02 and rel(tf["shots"]["sph5"], i["shots_first_stated"][1]) <= 0.02
    w = i["wall_s_per_shot_by_run"]
    assert math.isclose(tf["wall_s"], 1e2 / 0.09 * w["sph10"] + vr * 1e2 / 0.09 * w["sph5"])
    assert i["wall_first_result_s"] == (tf["wall_s"], tf["wall_s"])
    yr = 3.156e7
    assert rel(tf["wall_s"] / yr, 2.3251) <= 1e-3 and round(tf["wall_s"] / yr, 1) == i["wall_yr_first_stated"] == 2.3
    for nt, cs, wc in zip((4, 8), tc_["shots"], i["wall_campaign_s"]):
        assert cs == pytest.approx({"sph10": 1e2 / 0.09, "sph5": vr * 1e2 / 0.09, "cond": nt * 1e3})
        assert math.isclose(wc, sum(n * w[k] for k, n in cs.items()))
    lo, hi = (x / yr for x in i["wall_campaign_s"])
    assert rel(lo, 2.6787) <= 1e-3 and rel(hi, 3.0323) <= 1e-3
    assert (round(lo, 1), round(hi, 1)) == i["wall_yr_campaign_stated"] == (2.7, 3.0)          # '2.7-3.0 yr'
    olo, ohi = tc_["over_horizon"]
    assert ohi < 1 and (round(olo, 2), round(ohi, 2)) == i["campaign_over_horizon_stated"]      # 'inside the 5-year horizon'
    assert round(i["wall_yr_method_test_at_campaign_eps"]) == 21                                 # 'the test at 10% would take 21 yr'
    # record: the r25 rate campaign (4-8 temperatures at 10%, 10a at each, 5a x4 at one) as printed before the ruling
    r25 = i["wall_yr_campaign_r25_rate"]
    assert rel(r25[0], 36.078) <= 1e-3 and rel(r25[1], 60.001) <= 1e-3
    assert (round(r25[0]), round(r25[1])) == i["wall_yr_campaign_r25_rate_stated"]
    assert r33.shots == (tf["shots"]["sph10"],) * 2


def test_method_test_readout_and_reference(r33):
    """Author ruling 2026-10-05: what the method test measures and compares against, and what the Hadamard test costs."""
    i = r33.intermediates
    bg, cl = i["ncs_free_background"], i["ncs_free_background_classical"]
    assert round(cl[10.0], 5) == 1.3e-4 and 0.13 < cl[10.0] / bg[10.0] < 0.15          # 'gives 1.3e-4 at 10a'
    assert (1 - 0.3) * bg[10.0] > (1 + 0.3) * cl[10.0] * 3                               # '30%: well separated'
    assert 0.94 < i["ncs_free_background_zero_point_10a"] / bg[10.0] < 0.96               # 'mostly zero-point'
    lo, hi = i["ncs_free_background_temp_span_10a"]
    assert 0.95 < lo < 1 < hi < 1.05                                                      # 'less than 5% between Ta = 0.4 and 0.6'
    h = i["shots_first_hadamard"][10.0]
    assert math.isclose(h, i["ncs_hadamard_rel_var"][10.0] / 0.09) and round(h / 1e15, 1) == 1.2   # '1.2e15 shots at 10a'
    assert round(bg[5.0], 5) == 4.3e-4 and round(bg[10.0], 5) == 9.6e-4


def test_r25_2028_schedule_tiers_and_exports(a, r28):
    """(R10) three spatial-hop groups per slot, weight bits in parallel: 42 ancilla, F* >= 10. (R1) first result L5 = 4
    at 30%; campaign L5 = 2, 3, 4 at 10%."""
    i = r28.intermediates
    assert a.spatial_concurrency_2028.value == 3 and i["slots_per_step"] == 4 + 2 + 22 + 3 + 6 == 37
    lo, hi = i["depth_per_step"]
    rt = T_RUS_C
    assert math.isclose(lo, 10 * (rt + 10) + 5 * (rt + 6) + 22 * (rt + 8))
    assert math.isclose(hi, 10 * (rt + 31) + 5 * (rt + 7) + 22 * (rt + 10))
    assert rel(lo, i["depth_per_step_stated"][0]) <= 0.04 and rel(hi, i["depth_per_step_stated"][1]) <= 0.02
    # factory analysis: one group per slot 5.3-6.1, all rotations serial 1.84; the new schedule 10.4-13.0
    assert 5.29 < i["f_star_one_group_per_slot"][0] < 5.31 and 6.09 < i["f_star_one_group_per_slot"][1] < 6.11
    assert 1.83 < i["f_star_serial_rotations"] < 1.84
    assert 10.4 < i["f_star"][0] < 10.5 and 12.9 < i["f_star"][1] < 13.0
    assert i["t_per_shot"] == i["chunked_t_per_shot"] and i["t_depth_per_shot"] == (21 * lo, 21 * hi)
    assert i["shots_first"] == 125 and rel(i["wall_first_result_s"][0], i["wall_s_first_stated"]) <= 0.01   # '36 s'
    l5 = i["l5_campaign"]
    assert sorted(l5) == [2, 3, 4] and l5[4]["lq"] == 250 and l5[4]["t_shot"] == i["chunked_t_per_shot"]
    assert l5[2]["lq"] < l5[3]["lq"] < 250
    assert math.isclose(i["wall_campaign_s"][0], sum(v["wall_s"] for v in l5.values()))
    assert rel(i["wall_campaign_s"][0] / 60, i["wall_min_campaign_stated"]) <= 0.02                       # '13 min'
    # L5 = 3, 4 feed ten factories; L5 = 2 (k = 6 groups) does not (F* 8.9-10.0) and is priced at its corrected wall
    ok = {k: e["baseline_ok"][0] for k, e in i["exports_by_l5"].items()}
    assert ok == {2: False, 3: True, 4: True}
    assert l5[2]["wall_s"] > l5[2]["wall_s_baseline"] and l5[4]["wall_s"] == l5[4]["wall_s_baseline"]
    assert math.isclose(l5[2]["wall_s"], 1125 * (l5[2]["t_depth"][1] * 1e-5 + 1e-4))


def test_r25_chapter_prints_the_model(r28, r33):
    """The r25 prints in app06 (boxes and prose) are the model's numbers, rounded once."""
    tex = (Path(m.__file__).resolve().parents[2] / "applications" / "app06_chiral_gauge.tex").read_text()
    for frag in (r"$1015 = 216$ fermion $+\,486$ gauge $+\,270$ Higgs $+\,43$ ancilla",
                 r"$250 = 192$ fermion $+\,16$ gauge $+\,42$ ancilla",
                 r"$417$ queries per application and $1.8\times 10^{6}$ T",
                 r"$1.9\times 10^{10}$ at $t=10a$; $9.6\times 10^{9}$ at $t=5a$",
                 r"first result $36$~s; campaign $13$~min",
                 r"$2.3$~yr \\", r"$2.7$--$3.0$~yr for $4$--$8$ temperatures",
                 r"$1.1\times 10^{3}$ at $10a$ $+$ $5.5\times 10^{3}$ at $5a$, one temperature",
                 r"With the Hadamard test this needs $1.2\times 10^{15}$ shots at $10a$",
                 r"no zero-point part: $1.3\times 10^{-4}$ at $10a$",
                 r"checks the evolution and the estimator, not the thermal state", r"($430$--$2400\times$)", r"The test at $10\%$ would take $21$~yr",
                 r"has $4.9$ times the relative variance", r"less than $5\%$ between $Ta=0.4$ and $0.6$",
                 r"$\lesssim 3.6\times 10^{-11}$ (condensate), $\lesssim 5.4\times 10^{-12}$ ($N_{CS}$ test, $10a$; "
                 r"$1.0\times 10^{-11}$ at $5a$)",
                 r"$17$ filter calls per thermal state, $1.5\times 10^{3}$ applications and $2.8\times 10^{9}$ T",
                 # style pass 2026-10-08: the filter-scope sentence moved to the omission sentence after the 2033 box,
                 # the Hilbert-Schmidt-norm formula went (described in words), "budget" -> "2033 target"
                 r"$2.8$ times the 2033 target", r"Thermal-state systematic & not bounded",
                 r"preprint~\cite{AlvesLammLiu_inprep}"):
        assert frag in tex, frag
    assert "in preparation" not in tex and "N_{\\rm coupling}" not in tex
    # referee pass (2026-10-04): no lattice spacing in fm at the 2033 target, no HTL-only classical baseline, no
    # bounded thermal systematic, no 'arm'
    for gone in ("0.2\\,$fm", "hard-thermal-loop", "bounded by repeating", " arm", "Why classical methods are insufficient",
                 # author ruling 2026-10-05: no rate at 3^3, no linear-growth premise, no rate campaign
                 "Chiral condensate, sphaleron rate", "a rate only if", "Under linear growth", "ballistic",
                 "$36$--$60$~yr", "method-validation result", "sphaleron shot", "Per-shot T, sphaleron"):
        assert gone not in tex, gone
    i = r33.intermediates
    assert r33.lq[0] == 1015 and round(r33.hard_ops[0] / 1e10, 1) == 1.9
    assert round(i["t_per_shot_short"] / 1e9, 1) == 9.6 and round(i["tiers"]["first_result"]["wall_yr"], 1) == 2.3
    assert r28.lq[0] == 250 and round(r28.intermediates["wall_first_result_s"][0]) == 36


def test_utility_box_has_no_dollar_figure():
    # G5 open item 3 (Claude's decision, 2026-10-04): the $10M rested on Ch. 8's removed LISA figure; the box
    # attaches no dollar value and Table 1.2 has no Ch. 9 row
    p = Path(m.__file__).resolve().parents[2] / "applications" / "app06_chiral_gauge.tex"
    box = p.read_text().split("begin{utilitybox}")[1].split("end{utilitybox}")[0]
    # style pass 2026-10-08 (rule 9): "not priced in dollars" and "must not be added" retired with the disclaimers
    assert "\\$" not in box


# --------------------------------------------------------------------------- #
# Chapter-open pass (referee C1 / C3, 2026-10-05): thermal reference state and the Delta N_CS readout
# --------------------------------------------------------------------------- #

def test_free_overlap_spectrum_matches_closed_form():
    # h_ov = Gamma^0 + sgn(h_W) at U = 1 on 3^3, m = -1.5 (alpha = 7.5): eigenvalues +-sqrt(2 + 2c/N), twice each
    modes = m.free_overlap_modes(3, 3, -1.5, 1.0)
    assert len(modes) == 108
    vals = sorted(set(round(abs(x), 6) for x in modes))
    assert vals == [0.0, round(math.sqrt(2), 6), 1.88393, 1.946498]
    assert sum(1 for x in modes if abs(x) < 1e-12) == 4                    # p = 0: the massless mode, 4 spinor states
    assert sum(1 for x in modes if abs(abs(x) - math.sqrt(2)) < 1e-6) == 24  # 6 momenta x 4
    assert abs(sum(x for x in modes) ) < 1e-9                               # symmetric spectrum


def test_typical_reference_success_probability(r33):
    # C1: 1-design reference, post-selection exact; p = prod (1 + e^{-beta|eps|})/2 over the 216 free modes at beta = 2
    i = r33.intermediates
    t = i["typical_ref"]
    assert t["n_modes"] == 216 and i["typical_ref_n_frozen"] == 208
    assert abs(t["E0"] + 186.658) < 1e-3
    assert abs(i["typical_ref_ln_p"] + 137.934) < 1e-3                    # chapter 'p = e^{-138} ~ 1e-60'
    assert round(t["log10_p"]) == -60
    assert 7.0e29 < i["typical_ref_aa_rounds"] < 7.1e29                   # chapter '~1e30 rounds'
    assert round(i["typical_ref_ln_norm_offset"]) == -491                  # chapter 'a further e^{-491}'
    # independent check: ln p = S - ln D - beta (E - E0)
    assert abs(t["S"] - t["lnD"] - 2.0 * t["E_minus_E0"] - t["ln_p"]) < 1e-9
    assert abs(i["heralded_accept_min"] - 1 / 17) < 1e-12


def test_ncs_kubo_readout_and_units(r33):
    # C3: free-field non-topological <dN^2>, Hadamard-test normalization and variance, lattice units, small-box suppression
    i = r33.intermediates
    # verifier fix: the chapter's estimator (E link-averaged, clover B) on the free lattice theory
    bg = i["ncs_free_background"]
    assert abs(bg[10.0] - 9.577e-4) < 3e-6 and abs(bg[5.0] - 4.305e-4) < 3e-6      # verifier Kubo: 9.5775e-4, 4.3052e-4
    assert round(bg[10.0], 5) == 9.6e-4                                            # chapter '9.6e-4'
    assert abs(i["ncs_lambda"] - 729 / (8 * math.pi ** 2)) < 1e-12 and round(i["ncs_lambda"], 1) == 9.2   # 'lambda >= 9.2/a'
    assert 1.05e14 <= i["ncs_hadamard_rel_var"][10.0] < 1.15e14            # chapter '>~ 1.1e14'
    # helicity-A.B record (the first pass's 8.93e-3 / 1.2e12) and the forward-link pairing (t^2 growth)
    hel = i["ncs_free_background_helicity"]
    assert abs(hel[10.0] - 8.93e-3) < 1e-5 and abs(hel[5.0] - 3.69e-3) < 1e-5
    assert 1.2e12 <= i["ncs_hadamard_rel_var_helicity"][10.0] < 1.25e12   # chapter '1.2e12 even for an ideal helicity operator'
    fwd = i["ncs_free_background_fwd"]
    assert abs(fwd[10.0] - 2.564e-3) < 1e-5 and abs(fwd[40.0] - 2.42e-2) < 2e-4  # verifier 2.5641e-3, 2.4197e-2
    assert fwd[40.0] / fwd[10.0] > 9                                           # grows ~t^2, not bounded
    assert i["ncs_hadamard_rel_var_helicity"][10.0] < i["ncs_hadamard_rel_var_fwd10"] < i["ncs_hadamard_rel_var"][10.0]
    lo, hi = i["mr_small_box_suppression_range"]
    assert round(lo) == 431 and 2390 < hi < 2410                           # '~1e3-fold (0.0023 +- 0.0016 vs 1.68)'
    assert i["beta_lattice_classical"] == 8.0 and i["box_L_g2T"] == 1.5 and round(i["box_L_g2T_ew"], 2) == 0.63
    assert round(i["mr_small_box_suppression"]) == 730 and i["L_large_volume"] == 16
    assert i["inv_a_gev_at_tc"] == 318.0
    # the free-field background vanishes at t = 0 and is bounded in t (no diffusion in the free theory)
    assert m.ncs_free_field_background(3, 3, 0.5, 1.0, 1e-9) < 1e-12
    kap2 = (1 / (8 * math.pi ** 2)) ** 2
    assert max(m.ncs_free_field_background(3, 3, 0.5, 1.0, t) for t in (1, 3, 7, 20, 50)) <= kap2 / 4 * 156 * 2 * 1.01


def test_chapter_open_prints():
    tex = (Path(m.__file__).resolve().parents[2] / "applications" / "app06_chiral_gauge.tex").read_text()
    # symbol pass 2026-10-08: single-use p, d_beta and the alpha restatement went (rule 11); their anchors retired
    for frag in (r"$\sim 10^{30}$ rounds", r"a further $e^{-491}$",
                 r"$\lambda\ge 9.2/a$ at $g^2=1$", r"$9.6\times 10^{-4}$, mostly zero-point",
                 r"$\sigma^2\gtrsim 1.1\times 10^{14}$ ($1.2\times 10^{12}$ even for an ideal helicity operator",
                 r"suppressed $\sim 10^{3}$-fold ($0.0023\pm 0.0016$ against $1.68$) already at $L=3/(g^2T)$",
                 r"average over the links $(x,i)$ and $(x-\hat\imath,i)$", r"known cost and controlled bias",
                 r"orthonormal basis of the Gauss-law sector", r"($5\times 10^{-4}$ and $0.04$ at $n=8$)",
                 # style pass 2026-10-08: the box row now points to Algorithm choice ("amplified ... references" retired)
                 r"$a^{-1}\approx 320$~GeV", r"$L=1.5/(g^2T)$ ($0.63/(g^2T)$ at the electroweak $g^2\approx 0.42$)",
                 r"17 post-selection attempts at success probability $1/17$"):
        assert frag in tex, frag
    for gone in (r"$\sim 3\times 10^{7}$ amplification rounds", r"suppressed $730$-fold", r"$\lambda=9.2/a$",
                 r"$8.9\times 10^{-3}$, mostly zero-point", r"for a product reference the bias does not shrink", r"filtered amplitude of order $0.1$",
                 r"total systematic $\sim 40\%$", r"$\gtrsim 30\%$; $\Gamma_{\rm sph}$ at $3^3$ is validation-quality"):
        assert gone not in tex, gone
