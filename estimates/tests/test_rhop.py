"""Tests for the shared rule R-HOP (groups.rhop_link) and its phasing block (common.hwp_group).

r17 (2026-10-01): R-HOP is the LEGACY comparison (author ruling (e)); the default synthesis of rhop_link and
hwp_group is now the full fit "rus". The tests below pin the legacy values, so they pass "rus-slope" explicitly;
test_defaults_are_the_full_fit pins the fix.

The rule was derived for Ch. 6 in round D (apply_log/ch06_roundD.md) and ruled report-wide by
H. Lamm on 2026-09-29. These tests pin the function itself; the Ch. 6 numbers it reproduces are
pinned in test_ch06_collider.py. Refactor log: apply_log/rhop_shared.md.
"""
import dataclasses
import math

import pytest

from estimates import common as c
from estimates import ch06_collider as ch6
from estimates import ch09_chiral_gauge as ch9
from estimates.groups import (GROUPS, DIAGONAL_LINK, RHOP_MOVES_PER_FIELD, RHOP_N_PAULI_DIAGONAL,
                              RHOP_N_PAULI_MOVED, rhop_link)

EPS = 1e-4
T_ROT = 1.15 * math.log2(1e4)          # 15.28087 T, papers' convention


# --------------------------------------------------------------------------- #
# Hamming-weight phasing
# --------------------------------------------------------------------------- #

def test_hwp_counts_agree_with_ch9_helpers():
    """common.hwp_* carries the same three counts Ch. 9 defines for itself (no drift). Since r22 (E26) the
    ancilla count of both is k - w(k); the Ch. 9 helper is the one Ch. 6 imports."""
    for k in range(1, 257):
        assert c.hwp_synth_rotations(k) == ch9.hwp_synth_rotations(k)
        assert c.hwp_toffolis(k) == ch9.hwp_toffolis(k)
        assert c.hwp_ancilla(k) == ch9.hwp_ancilla(k)


def test_hwp_group_values():
    h = c.hwp_group(9, EPS, "rus-slope")
    assert (h["k"], h["n_rot"], h["n_toffoli"], h["ancilla"]) == (9, 4, 7, 7)        # E26: ancilla = k - w(k), was 11
    assert h["t_rot"] == pytest.approx(T_ROT, rel=1e-12)
    assert h["t"] == pytest.approx(4 * T_ROT + 7 * 7, rel=1e-12)          # 110.1235
    assert h["t"] == pytest.approx(110.1235, abs=1e-4)
    h3 = c.hwp_group(3, EPS, "rus-slope")
    assert (h3["n_rot"], h3["n_toffoli"], h3["ancilla"]) == (2, 1, 1)
    assert h3["t"] == pytest.approx(2 * T_ROT + 7, rel=1e-12)              # 37.56
    h12 = c.hwp_group(12, EPS, "rus-slope")
    assert (h12["n_rot"], h12["n_toffoli"], h12["ancilla"]) == (4, 10, 10)             # was 14


def test_hwp_group_conventions_are_explicit():
    assert c.hwp_group(9, EPS, "rus")["t_rot"] == pytest.approx(T_ROT + 9.2)
    assert c.hwp_group(9, EPS, "rus-slope", toffoli_convention="jones")["t"] == pytest.approx(4 * T_ROT + 7 * 4)
    with pytest.raises(ValueError):
        c.hwp_group(0, EPS, "rus-slope")


# --------------------------------------------------------------------------- #
# R-HOP
# --------------------------------------------------------------------------- #

def test_rule_constants():
    assert DIAGONAL_LINK == frozenset({"Z3"})
    assert RHOP_N_PAULI_DIAGONAL == {"Z3": 8}
    assert RHOP_N_PAULI_MOVED == 2
    assert RHOP_MOVES_PER_FIELD == 2


def test_z3_diagonal_link_has_no_move_and_eight_strings():
    h = rhop_link("Z3", 3, 3, EPS, synthesis="rus-slope")
    assert (h["k"], h["n_moves"], h["t_move"], h["t_moves"], h["n_pauli"]) == (9, 0, 0.0, 0.0, 8)
    assert h["t"] == h["t_phasing"] == 8 * c.hwp_group(9, EPS, "rus-slope")["t"]
    assert h["t"] == pytest.approx(880.99, abs=0.01)                       # app10:88 '881'


def test_sigma72_hop_is_moves_at_umul_plus_two_strings():
    h = rhop_link("S72x3", 3, 3, EPS, synthesis="rus-slope")
    assert (h["n_moves"], h["t_move"], h["n_pauli"]) == (6, 1120.0, 2)
    assert h["t_moves"] == 6720.0
    assert h["t_phasing"] == pytest.approx(220.247, abs=1e-3)
    assert h["t"] == pytest.approx(6940.25, abs=0.01)                      # app10:90 '6.9e3'


@pytest.mark.parametrize("group", ["2T", "2O", "S36x3", "S72x3"])
@pytest.mark.parametrize("n_stag", [1, 2, 3])
def test_rule_formula_for_moved_groups(group, n_stag):
    """T_hop = 2 N_stag U_mul + 2 HWP_T(N_c N_stag), U_mul from the group's published table."""
    h = rhop_link(group, n_stag, 3, EPS, synthesis="rus-slope")
    umul = GROUPS[group].primitives["U_mul"].t(EPS)
    assert h["t"] == pytest.approx(2 * n_stag * umul + 2 * c.hwp_group(3 * n_stag, EPS, "rus-slope")["t"], rel=1e-12)
    assert h["t"] == h["t_moves"] + h["t_phasing"]


def test_overrides_k_npauli_moves():
    base = rhop_link("S36x3", 1, 3, EPS, synthesis="rus-slope")
    assert base["k"] == 3
    assert rhop_link("S36x3", 1, 3, EPS, k=12, synthesis="rus-slope")["hwp"]["k"] == 12
    assert rhop_link("Z3", 1, 3, EPS, k=12, n_pauli=8, synthesis="rus-slope")["t"] == 8 * c.hwp_group(12, EPS, "rus-slope")["t"]
    assert rhop_link("S36x3", 3, 3, EPS, moves_per_field=1, synthesis="rus-slope")["n_moves"] == 3
    assert rhop_link("Z3", 3, 3, EPS, moves_per_field=5, synthesis="rus-slope")["n_moves"] == 0   # diagonal: moves never priced


def test_rejects_groups_it_cannot_price():
    with pytest.raises(KeyError):
        rhop_link("SU4", 1, 3, EPS, synthesis="rus-slope")
    with pytest.raises(ValueError):
        rhop_link("SU3", 1, 3, EPS, synthesis="rus-slope")          # no published U_mul
    with pytest.raises(ValueError):
        rhop_link("S72x3", 0, 3, EPS, synthesis="rus-slope")


# --------------------------------------------------------------------------- #
# Ch. 6 goes through the shared function, bit-identically
# --------------------------------------------------------------------------- #

def _ch6_hop_link_before_refactor(a, group):
    """The Ch. 6 code as it stood before the refactor (ch06_collider.py round D), kept verbatim in logic."""
    eps = float(a.eps_rot)
    n_stag = int(a.n_stag.value)
    k = int(a.n_c.value) * n_stag
    t_rot = c.t_per_rotation(eps, "rus-slope")      # Ch. 6 synthesis_model before E20 (legacy R-HOP record)
    n_rot, n_tof = ch9.hwp_synth_rotations(k), ch9.hwp_toffolis(k)
    h = {"k": k, "n_rot": n_rot, "n_toffoli": n_tof, "ancilla": ch9.hwp_ancilla(k), "t_rot": t_rot,
         "t": n_rot * t_rot + c.toffoli_t(n_tof, a.toffoli_convention.value)}
    if group in {"Z3"}:
        n_moves, t_move, n_p = 0, 0.0, int(a.n_pauli_hop_diagonal.value)
    else:
        n_moves = int(a.moves_per_field.value) * n_stag
        t_move = GROUPS[group].primitives["U_mul"].t(eps)
        n_p = int(a.n_pauli_hop_moved.value)
    t_move_total = n_moves * t_move
    t_phase = n_p * h["t"]
    return {"k": k, "hwp": h, "n_pauli": n_p, "n_moves": n_moves, "t_move": t_move, "t_moves": t_move_total,
            "t_phasing": t_phase, "t": t_move_total + t_phase}


@pytest.mark.parametrize("group", ["Z3", "S72x3"])
def test_ch6_hop_link_is_bit_identical_through_the_shared_function(group):
    a = ch6.Assumptions()
    assert ch6.hop_link_rhop(a, group) == _ch6_hop_link_before_refactor(a, group)     # exact ==, no tolerance
    assert ch6.hop_link_rhop(a, group) == rhop_link(group, 3, 3, float(a.eps_rot), synthesis="rus-slope")   # Ch. 6 inputs; its synthesis is the legacy slope-only
    # at the rule's eps = 1e-4 (Ch. 6 moved to 1e-3 in round E) the values are the round-D ones
    a4 = dataclasses.replace(a, eps_rot=c.Stated(EPS, "app10:88", "test"))
    assert ch6.hop_link_rhop(a4, group) == _ch6_hop_link_before_refactor(a4, group) == rhop_link(group, 3, 3, EPS, synthesis="rus-slope")


def test_ch6_stated_inputs_equal_the_rule_constants():
    a = ch6.Assumptions()
    assert int(a.n_pauli_hop_diagonal.value) == RHOP_N_PAULI_DIAGONAL["Z3"]
    assert int(a.n_pauli_hop_moved.value) == RHOP_N_PAULI_MOVED
    assert int(a.moves_per_field.value) == RHOP_MOVES_PER_FIELD


def test_ch6_headlines_unchanged():
    """Since r17 (apply_log/r17_ch06.md) Ch. 6 prices its hop with groups.hop_link_cost; the R-HOP prints are kept as
    records (ch6.hop_link_rhop, the *_pre_r17 / *_round_e / *_at_1e-4 intermediates) and still reproduce."""
    a = ch6.Assumptions()
    r28, r33 = ch6.model(a, "2028"), ch6.model(a, "2033")
    assert r28.intermediates["t_hop_per_link_2028_round_e"] * 32 == pytest.approx(24279.71, abs=0.01)
    assert r28.intermediates["t_hop_per_link_2028_at_1e-4"] * 32 == pytest.approx(28191.61, abs=0.01)
    assert r28.intermediates["t_per_step_2028_pre_r17"] == pytest.approx(32429.49, abs=0.01)
    assert r28.intermediates["t_total_2028_pre_r17"] == pytest.approx(97288.48, abs=0.1)
    assert r28.intermediates["t_per_step_2028_round_e"] == pytest.approx(33446.26, abs=0.01)
    assert r28.intermediates["t_total_2028_round_e"] == pytest.approx(100338.79, abs=0.1)
    assert r28.intermediates["twosteps_t_total_2028"] == pytest.approx(99509.45, abs=0.1)
    assert r33.intermediates["pre_r17_t_total_2033"][0] == pytest.approx(2.767815e9, rel=1e-6)
    # r18 headlines (E20, 1 ramp + 2 evolution steps, R-TOL eps 1.6e-3): Z3 hop 8 x HWP(9) at the full fit,
    # Sigma(72x3) hop from the draft. r17 was 32,532.37 / 84,027.89 (1 ramp + 1 evolution step).
    assert r28.intermediates["t_hop_per_step_2028"] == pytest.approx(32876.79, abs=0.01)
    assert r28.hard_ops[0] == pytest.approx(132174.88, abs=0.1)
    assert r28.intermediates["t_total_2028_two_steps"] == pytest.approx(87266.29, abs=0.1)
    # r19 headlines (E21: Sigma(72x3) squish held once per link; r18 share="draft" was 1.1912865e7).
    # r25 (R4, window 30-50 a): the tolerance moves with the shorter shot; r19-r24 pin was 7.1452493e6.
    # r26 (J1, backward readout ramp): more rotations per shot, tighter R-TOL eps; r25 pin was 6.9877264e6.
    assert r33.intermediates["t_hop_per_step_2033"] == pytest.approx(7.0224735e6, rel=1e-6)


# --------------------------------------------------------------------------- #
# The rule's per-link values at the other chapters' inputs (derived here; NOT applied to any chapter)
# --------------------------------------------------------------------------- #

def test_rule_values_at_other_chapters_inputs():
    """Per link per step, papers' convention, eps = 1e-4. Predictions for rhop_shared.md only."""
    # Ch. 5 2033: Sigma(72x3), N_stag = 3, N_c = 3 (app03:137-138)
    assert rhop_link("S72x3", 3, 3, EPS, synthesis="rus-slope")["t"] * 81 == pytest.approx(562160.0, abs=1.0)
    # Ch. 10 2028: Sigma(36x3), N_stag = 1 (app07:139); 8 links
    h10a = rhop_link("S36x3", 1, 3, EPS, synthesis="rus-slope")
    assert h10a["t"] == pytest.approx(2 * 308 + 2 * (2 * T_ROT + 7), rel=1e-12)   # 691.12
    assert 8 * h10a["t"] == pytest.approx(5528.99, abs=0.01)
    # Ch. 10 2033: Sigma(36x3), N_stag = 3 (app07:163); 24 links
    h10b = rhop_link("S36x3", 3, 3, EPS, synthesis="rus-slope")
    assert 24 * h10b["t"] == pytest.approx(49637.9, abs=0.1)
    # Ch. 9 2028: Z3, k = N_c L5 = 12 by its grouping ruling (R6), 8 links
    assert 8 * rhop_link("Z3", 1, 3, EPS, k=12, synthesis="rus-slope")["t"] == 64 * c.hwp_group(12, EPS, "rus-slope")["t"]


def test_defaults_are_the_full_fit():
    """r17: the slope-only default (15.3 T per rotation at 1e-4) was a latent bug; the default is now 24.48 T."""
    assert c.hwp_group(9, EPS)["t_rot"] == pytest.approx(T_ROT + 9.2, rel=1e-12)
    assert rhop_link("Z3", 3, 3, EPS)["t"] == pytest.approx(8 * (4 * (T_ROT + 9.2) + 49), rel=1e-12)   # 1175.37
    assert rhop_link("Z3", 3, 3, EPS)["t"] == rhop_link("Z3", 3, 3, EPS, synthesis="rus")["t"]
