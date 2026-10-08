"""Tests for the r17 fermion hop and mass (groups.hop_link_cost, fermion_hop_counts, staggered_mass_site).

Author ruling (e), H. Lamm 2026-10-01: adopt the unpublished Fermion_Primitives gate counts report-wide in place
of R-HOP and estimate the missing pieces. Part 1 shows which printed numbers of the draft rebuild from its own
gate tables in its own conventions (so the counts in groups.FERMION_DRAFT are read correctly). Part 2 pins the
draft's internal inconsistencies as found. Part 3 pins the report-rule counts, every estimated piece labeled.
Log: apply_log/r17_core.md.

r19 (author ruling E21 (1), 2026-10-01, "make the switch"): the default is share="link" + undo=True (color squish +
parity and hop squish + flags once per link, held through V(x), V(y), hop, V^dag(x), V^dag(y)); the r17
per-application recompute is the sensitivity share="draft", still pinned below. n_spin = 2 (Wilson d=3) and 2 N_angles
color rotations are author-confirmed (E21 (2), (3)). Log: apply_log/r19_core.md.
"""
import math

import pytest

from estimates import common as c
from estimates.groups import (FERMION_DRAFT, fermion_hop_counts, fp_printed, fp_reconstruct, hop_link_cost,
                              mcx_table_toffolis, mcx_toffolis, rhop_link, spinor_counts, staggered_mass_site,
                              wilson_mass_rotations)

GROUPS4 = ("2T", "2O", "S36x3", "S72x3")
EPS = 1e-4
L = math.log2(1 / EPS)
T_FULL = 1.15 * L + 9.2                       # 24.4809 T per rotation at 1e-4


# --------------------------------------------------------------------------- #
# 1. The draft's printed numbers rebuilt from its gate tables (self-consistent ones)
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("group,printed", [("2T", 189), ("2O", 917), ("S36x3", 7854), ("S72x3", 12789)])
def test_colour_squish_printed_T_is_the_table_at_2n_minus_4(group, printed):
    """su2_diag.tex:207-211, su3_diag.tex:86-92 equal the C^nX tables at 2(n-2) Toffolis (1 per C^2X)."""
    assert fp_reconstruct(group, "col_squish_t") == (printed, 0.0)
    assert FERMION_DRAFT[group].col_squish_t_printed == printed


@pytest.mark.parametrize("group,printed", [("2T", 84), ("2O", 434)])
def test_su2_hop_squish_printed_T_is_the_table_at_2n_minus_3(group, printed):
    """hopping.tex:202-213: the T column equals the C^nX columns at 2n-3."""
    assert fp_reconstruct(group, "hop_squish_t") == (printed, 0.0)


def test_su3_hop_squish_commented_values_are_the_table_at_2n_minus_4():
    """hopping.tex:262-263 comment out 469 and 616; those are the table at 2(n-2)."""
    assert fp_reconstruct("S36x3", "hop_squish_t") == (469, 0.0)
    assert fp_reconstruct("S72x3", "hop_squish_t") == (616, 0.0)


def test_colour_frame_constants_rebuild():
    assert fp_reconstruct("2O", "C_V") == pytest.approx((2778, 69), abs=1e-9)            # su2_diag.tex:236
    assert fp_reconstruct("2O", "C_rot") == pytest.approx((944, 69), abs=1e-9)           # :229
    assert fp_reconstruct("2T", "C_V") == pytest.approx((1164, 59.8), abs=0.5)           # :235 (786.4 printed 786)
    assert fp_reconstruct("S36x3", "C_G") == pytest.approx((10918.8, 188.6), abs=1e-9)   # su3_diag.tex:132 (flags at 2 T)
    assert fp_reconstruct("S72x3", "C_G") == pytest.approx((20725.8, 211.6), abs=1e-9)   # :133
    assert fp_reconstruct("S36x3", "C_rot") == pytest.approx((1508.8, 188.6), abs=1e-9)
    assert fp_reconstruct("S72x3", "C_rot") == pytest.approx((1692.8, 211.6), abs=1e-9)


@pytest.mark.parametrize("group", GROUPS4)
def test_text_hop_cost_rebuilds_as_two_squish_plus_diagonalizers(group):
    """hopping.tex:335-347 = 2 (squish + flags) + n_diag x C(T), C(T) = 38.8 + 4.45 L as printed."""
    assert fp_reconstruct(group, "C_hop_text") == pytest.approx(FERMION_DRAFT[group].printed["C_hop_text"], abs=1e-9)


def test_diagonalizer_reading():
    """38.8 = 2 T (controlled H) + 4 x 9.2: four rotations. The printed slope 4.45 is not 4 x 1.15 = 4.6."""
    assert 2 + 4 * 9.2 == pytest.approx(38.8)
    assert 4.45 / 1.15 != pytest.approx(round(4.45 / 1.15))


# --------------------------------------------------------------------------- #
# 2. The draft's inconsistencies, pinned as found
# --------------------------------------------------------------------------- #

def test_staggered_CS_printed_is_two_CG_plus_two_table_Chop():
    """resources.tex:31-42 vs :60-67. Sigma groups: C^S = 2 C^G + 2 C^hop(tab), exactly. BT/BO constants
    match 4 C^G + 2 C^hop while their slopes match 2 + 2. The formula on :15 says 2 C^G + C^hop."""
    for g, cg in (("S36x3", (10918.8, 188.6)), ("S72x3", (20725.8, 211.6))):
        a, b = FERMION_DRAFT[g].printed["C_S"]
        ha, hb = FERMION_DRAFT[g].printed["C_hop_tab"]
        assert (a, b) == pytest.approx((2 * cg[0] + 2 * ha, 2 * cg[1] + 2 * hb), abs=1e-9)
    for g in ("2T", "2O"):
        p = FERMION_DRAFT[g].printed
        assert p["C_S"][0] == pytest.approx(4 * p["C_G_tab"][0] + 2 * p["C_hop_tab"][0], abs=1e-9)
        assert p["C_S"][1] == pytest.approx(2 * p["C_G_tab"][1] + 2 * p["C_hop_tab"][1], abs=1e-9)


def test_su2_table_CG_drops_parity_and_offsets():
    """resources.tex:61-62 C^G = 2 x squish + slope only; su2_diag.tex C_V adds 7 x parity + 9.2 per rotation."""
    for g in ("2T", "2O"):
        f = FERMION_DRAFT[g]
        assert f.printed["C_G_tab"] == pytest.approx((2 * f.col_squish_t_printed, 1.15 * f.col_rot), abs=1e-9)


def test_su2_table_Chop_is_one_squish_and_n_cls_N_diagonalizers():
    for g in ("2T", "2O"):
        f = FERMION_DRAFT[g]
        n_rot = 4 * f.n_classes * f.n
        assert f.printed["C_hop_tab"] == pytest.approx((f.hop_squish_t_printed + 9.2 * n_rot, 1.15 * n_rot), abs=1e-9)
    assert FERMION_DRAFT["2T"].n_diag_text == 10 and FERMION_DRAFT["2O"].n_diag_text == 28   # T only vs T and T^dag


def test_su3_table_Chop_constants_exceed_half_the_text():
    """Slopes halve exactly; constants exceed text/2 by 102 (S36) and 160 (S72) (the investigation's 204 / 320)."""
    for g, extra in (("S36x3", 102), ("S72x3", 160)):
        t, tab = FERMION_DRAFT[g].printed["C_hop_text"], FERMION_DRAFT[g].printed["C_hop_tab"]
        assert tab[1] == pytest.approx(t[1] / 2)
        assert tab[0] - t[0] / 2 == pytest.approx(extra)


def test_su3_printed_hop_squish_is_not_integer_toffolis():
    assert 136 % 7 and 372 % 7


def test_bo_colour_squish_under_each_convention():
    t = FERMION_DRAFT["2O"].col_squish_mcx
    assert mcx_table_toffolis(t, "draft-color") == 131      # 917 T, printed
    assert mcx_table_toffolis(t, "2n-3") == 163              # 1141 T
    assert mcx_table_toffolis(t, "mbu") == 105               # 735 T, used


# --------------------------------------------------------------------------- #
# 3. Report-rule counts (MBU, full fit, undo estimated)
# --------------------------------------------------------------------------- #

def test_mcx_conventions():
    assert [mcx_toffolis(n) for n in range(1, 8)] == [0, 1, 2, 3, 4, 5, 6]
    assert [mcx_toffolis(n, "2n-3") for n in range(1, 8)] == [0, 1, 3, 5, 7, 9, 11]
    assert [mcx_toffolis(n, "draft-color") for n in range(1, 8)] == [0, 1, 2, 4, 6, 8, 10]
    with pytest.raises(ValueError):
        mcx_toffolis(3, "bogus")


@pytest.mark.parametrize("group,col,hop", [("2T", 23, 8), ("2O", 105, 39), ("S36x3", 759, 46), ("S72x3", 1219, 63)])
def test_squish_toffolis_under_mbu(group, col, hop):
    f = FERMION_DRAFT[group]
    assert mcx_table_toffolis(f.col_squish_mcx) == col
    assert mcx_table_toffolis(f.hop_squish_mcx) == hop


@pytest.mark.parametrize("group,n_stag,tof,tdir,rot", [
    ("2T", 1, 109, 40, 291), ("2O", 1, 347, 56, 355),
    ("S36x3", 1, 2592, 120, 899), ("S36x3", 3, 2604, 360, 2693),
    ("S72x3", 1, 3660, 120, 979), ("S72x3", 3, 3672, 360, 2933)])
def test_staggered_counts_per_link(group, n_stag, tof, tdir, rot):
    """Default since r19 (E21 (1)): share="link", undo=True."""
    h = fermion_hop_counts(group, n_stag)
    assert (h["share"], h["undo"]) == ("link", True)
    assert (h["toffoli"], h["t_direct"], h["n_rot"]) == (tof, tdir, rot)
    assert h["toffoli"] == sum(i["toffoli"] for i in h["items"])


@pytest.mark.parametrize("group,n_stag,tof,tdir,rot", [
    ("2T", 1, 379, 40, 291), ("2O", 1, 1145, 56, 355),
    ("S36x3", 1, 9480, 120, 899), ("S36x3", 3, 10076, 360, 2693),
    ("S72x3", 1, 13650, 120, 979), ("S72x3", 3, 14314, 360, 2933)])
def test_staggered_counts_per_link_draft_sensitivity(group, n_stag, tof, tdir, rot):
    """The r17 headline, kept as the conservative sensitivity: share="draft" (squish + parity per frame
    application, hop squish + flags per field). Rotation counts equal the default's; only the Toffolis move."""
    h = fermion_hop_counts(group, n_stag, share="draft")
    assert (h["toffoli"], h["t_direct"], h["n_rot"]) == (tof, tdir, rot)
    assert h["n_rot"] == fermion_hop_counts(group, n_stag)["n_rot"]


def test_staggered_breakdown_2O():
    items = {i["name"]: (i["toffoli"], i["t_direct"], i["n_rot"]) for i in fermion_hop_counts("2O", 1)["items"]}
    assert items == {"colour_squish_and_parity": (2 * 105 + 56, 0, 0),           # once per link (share="link")
                     "colour_rotations": (0, 0, 4 * 60),                      # 4 frame applications incl. V_g^dag
                     "hop_squish_and_flags": (2 * 39, 0, 0),
                     "hop_diagonalizers": (0, 2 * 28, 4 * 28),                # 28 = 2 x 7 classes x N=2
                     "hop_phasing_hwp": (3, 0, 3)}                            # HWP(k=4)
    items = {i["name"]: (i["toffoli"], i["t_direct"], i["n_rot"])
             for i in fermion_hop_counts("2O", 1, share="draft")["items"]}
    assert items == {"colour_squish_and_parity": (4 * (2 * 105 + 56), 0, 0),   # per frame application (draft)
                     "colour_rotations": (0, 0, 4 * 60),
                     "hop_squish_and_flags": (2 * 39, 0, 0),
                     "hop_diagonalizers": (0, 2 * 28, 4 * 28),                # 28 = 2 x 7 classes x N=2
                     "hop_phasing_hwp": (3, 0, 3)}                            # HWP(k=4)


def test_undo_and_sharing_variants():
    # default (share="link"): the undo adds only the V_g^dag rotation layers; the g-only pieces are held
    full, no_undo = fermion_hop_counts("2O", 1), fermion_hop_counts("2O", 1, undo=False)
    assert full["toffoli"] == no_undo["toffoli"]
    assert full["n_rot"] - no_undo["n_rot"] == 2 * 60
    # draft sensitivity: the undo also re-squishes per application
    fd, nd = fermion_hop_counts("2O", 1, share="draft"), fermion_hop_counts("2O", 1, share="draft", undo=False)
    assert fd["toffoli"] - nd["toffoli"] == 2 * (2 * 105 + 56)
    assert fd["n_rot"] - nd["n_rot"] == 2 * 60
    link = fermion_hop_counts("S72x3", 3)
    assert link == fermion_hop_counts("S72x3", 3, share="link", undo=True)
    assert link["toffoli"] == (2 * 1219 + 892) + 2 * (63 + 100) + 16
    draft = fermion_hop_counts("S72x3", 3, share="draft")
    assert draft["toffoli"] - link["toffoli"] == 3 * (2 * 1219 + 892) + 2 * 2 * (63 + 100)
    assert link["n_rot"] == draft["n_rot"]
    lit = fermion_hop_counts("S72x3", 3, frame_rot_per_field=False)            # literal "shared C^G", recorded only
    assert fermion_hop_counts("S72x3", 3)["n_rot"] - lit["n_rot"] == 4 * 2 * 184


def test_phasing_options():
    hwp = fermion_hop_counts("S36x3", 3)
    plain = fermion_hop_counts("S36x3", 3, phasing="plain")
    none = fermion_hop_counts("S36x3", 3, phasing="none")
    assert plain["n_rot"] - none["n_rot"] == 18 and plain["toffoli"] == none["toffoli"]
    assert hwp["n_rot"] - none["n_rot"] == c.hwp_synth_rotations(18) == 5
    assert hwp["toffoli"] - none["toffoli"] == c.hwp_toffolis(18) == 16


def test_wilson_2O_d3():
    """Ch. 9 2033 kernel (2O, Wilson, d=3, block encoding: no phasing). n_spin = 2 at r = 1, author-confirmed
    (E21 (2)). Default share="link"; the draft sensitivity re-squishes per application."""
    w = fermion_hop_counts("2O", 1, fermion="wilson", d=3, phasing="none")
    assert w["n_spin"] == 2
    assert (w["toffoli"], w["t_direct"], w["n_rot"]) == (384, 144, 704)
    assert 384 == (2 * 105 + 56) + 2 * 39 + 40                                  # frame g-only + hop squish + spinor
    w1 = fermion_hop_counts("2O", 1, fermion="wilson", d=3, phasing="none", n_spin=1)
    assert (w1["toffoli"], w1["t_direct"], w1["n_rot"]) == (384, 88, 352)
    wd = fermion_hop_counts("2O", 1, fermion="wilson", d=3, phasing="none", share="draft")
    assert (wd["toffoli"], wd["t_direct"], wd["n_rot"]) == (1182, 144, 704)
    assert spinor_counts(3, 2) == {"toffoli": 40, "t_direct": 32, "n_rot": 0}       # 4 W x 5 x N (MBU)
    assert spinor_counts(3, 2, "2n-3") == {"toffoli": 80, "t_direct": 32, "n_rot": 0}   # FP's 2 x 296 = 592 T
    assert 7 * 80 + 32 == 592
    with pytest.raises(ValueError):
        fermion_hop_counts("2O", 1, fermion="wilson")


def test_z3_keeps_rhop_structure_at_full_fit():
    h = hop_link_cost("Z3", 3, eps=EPS)
    assert (h["toffoli"], h["n_rot"]) == (8 * 7, 8 * 4)
    assert h["t"] == pytest.approx(rhop_link("Z3", 3, 3, EPS, synthesis="rus")["t"], rel=1e-12)
    assert h["t"] == pytest.approx(8 * (4 * T_FULL + 49), rel=1e-12)
    assert hop_link_cost("Z3", 4, eps=EPS, k=12)["n_rot"] == 32


def test_cost_at_eps_and_rtol():
    h = hop_link_cost("2O", 1, eps=EPS)
    assert h["t"] == pytest.approx(7 * 347 + 56 + 355 * T_FULL, rel=1e-12)           # 11175.71
    assert h["t"] == pytest.approx(11175.71, abs=0.01)
    assert hop_link_cost("2O", 1, eps=EPS, share="draft")["t"] == pytest.approx(16761.71, abs=0.01)
    assert sum(i["t"] for i in h["items"]) == pytest.approx(h["t"], rel=1e-12)
    r = hop_link_cost("S72x3", 3, n_rot_circuit=1e6)
    assert r["eps"] == pytest.approx(1e-4, rel=1e-12)
    assert r["t"] == pytest.approx(7 * 3672 + 360 + 2933 * T_FULL, rel=1e-12)        # 97866.39
    with pytest.raises(ValueError):
        hop_link_cost("2O", 1)


def test_legacy_comparison_carried():
    h = hop_link_cost("S72x3", 3, eps=EPS)
    old = rhop_link("S72x3", 3, 3, EPS, synthesis="rus-slope")["t"]
    assert h["legacy"]["rhop_rus_slope"] == old == pytest.approx(6940.25, abs=0.01)   # app10:90 '6.9e3'
    assert h["legacy"]["ratio_new_over_rhop_rus_slope"] == pytest.approx(h["t"] / old)
    assert 14 < h["legacy"]["ratio_new_over_rhop_rus_slope"] < 14.2             # 14.10 (share="link"; draft 24.8)


def test_fp_printed_helper():
    assert fp_printed("2O", "C_S", EPS) == pytest.approx(9234.4 + 266.8 * L)


# --------------------------------------------------------------------------- #
# Mass
# --------------------------------------------------------------------------- #

def test_staggered_mass_site():
    m = staggered_mass_site(3, 3, EPS)
    assert (m["k"], m["n_rot"], m["n_toffoli"]) == (9, 4, 7)
    assert m["t"] == pytest.approx(4 * T_FULL + 49, rel=1e-12)
    assert staggered_mass_site(3, 1, n_rot_circuit=1e6)["t"] == pytest.approx(2 * T_FULL + 7, rel=1e-12)
    with pytest.raises(ValueError):
        staggered_mass_site(3, 1)


def test_wilson_mass_rotations():
    assert wilson_mass_rotations(2, 1, 1) == 4 and wilson_mass_rotations(2, 1, 2) == 4
    assert wilson_mass_rotations(2, 1, 3) == 8 and wilson_mass_rotations(3, 2, 4) == 24


# --------------------------------------------------------------------------- #
# r19: the judge's per-link values (frame-undo reading, 2026-10-01), recomputed exactly
# --------------------------------------------------------------------------- #
EPS_J = math.sqrt(1e-2 / 1.32624e8)           # 8.683e-6, the judge's eps


@pytest.mark.parametrize("args,kw,link_true,link_false,draft_true,draft_false", [
    (("S72x3", 3), {}, (3672, 360, 2933, 109758.06), 78255.08, (14314, 360, 2933, 184252.06), 106129.08),
    (("S36x3", 3), {}, (2604, 360, 2693, 95433.58), 67354.84, (10076, 360, 2693, 147737.58), 87514.84),
    (("2O", 1, 3), {"fermion": "wilson"}, (391, 144, 708, 23084.00), 16235.52, (1189, 144, 708, 28670.00),
     18097.52)])
def test_judge_per_link_values(args, kw, link_true, link_false, draft_true, draft_false):
    """Headline (share="link", undo=True): Sigma72 x3 1.10e5 T, Sigma36 x3 9.5e4 T, 2O Wilson d=3 2.3e4 T at
    eps 8.7e-6. The draft variant (share="draft") and the no-undo readings are sensitivities."""
    h = hop_link_cost(*args, eps=EPS_J, legacy=False, **kw)
    assert (h["toffoli"], h["t_direct"], h["n_rot"]) == link_true[:3]
    assert h["t"] == pytest.approx(link_true[3], abs=0.01)
    assert h["t"] == pytest.approx(7 * h["toffoli"] + h["t_direct"] + h["n_rot"] * (1.15 * math.log2(1 / EPS_J) + 9.2),
                                   rel=1e-12)
    assert hop_link_cost(*args, eps=EPS_J, legacy=False, undo=False, **kw)["t"] == pytest.approx(link_false, abs=0.01)
    d = hop_link_cost(*args, eps=EPS_J, legacy=False, share="draft", **kw)
    assert (d["toffoli"], d["t_direct"], d["n_rot"]) == draft_true[:3]
    assert d["t"] == pytest.approx(draft_true[3], abs=0.01)
    assert hop_link_cost(*args, eps=EPS_J, legacy=False, share="draft", undo=False, **kw)["t"] == \
        pytest.approx(draft_false, abs=0.01)
    assert d["t"] > h["t"]                    # draft is the conservative sensitivity


def test_su3_rotations_are_2_n_angles():
    """E21 (3): SU(3) color rotations per frame application = 2 N_angles (su3_diag.tex:194-198; the text's
    '6 N_angles' at :120 is a draft slip). Printed rotation log coefficient / 1.15 = the model's col_rot."""
    for g, n_angles in (("S36x3", 82), ("S72x3", 92)):
        assert FERMION_DRAFT[g].col_rot == 2 * n_angles
        assert FERMION_DRAFT[g].printed["C_rot"][1] / 1.15 == pytest.approx(2 * n_angles)
