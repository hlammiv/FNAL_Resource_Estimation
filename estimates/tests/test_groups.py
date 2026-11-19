"""Sigma(216x3) group entry (arxiv_2611_08161).

Checks the entry against the paper's printed numbers: Table 11 (primitives), Table 13 (multiplicities, U_Ph at 1/2),
Table 14 (closed forms, synthesis-error coefficients, fiducial totals 2.3e10, 6.8e10). Table 11 and Table 14 print
their constants with 9.2 T per rotation inside; PrimitiveCost.t and c_t_full must reproduce them as printed.
"""

import math

import pytest

from estimates import common as c
from estimates import groups as G

S = "S216x3"
EPS = (1e-4, 1e-8, 1e-10)


def test_register_and_primitives_as_printed():
    g = G.GROUPS[S]
    assert g.order == 648 and g.link_qubits == 11 and not g.link_width_conflict
    p = g.primitives
    # t_const = Table 11's constant less 9.2 T per rotation (rotations = printed log coefficient / 1.15)
    rows = {"U_inv": (2065, 0, 4), "U_mul": (2548, 0, 5), "U_Tr": (30786, 17.25, 10),
            "U_FFT": (9678, 489.9, 9), "U_phi": (238, 11.5, 18)}
    assert set(p) == set(rows)
    for k, (a, b, anc) in rows.items():
        assert (p[k].t_const, p[k].t_log, p[k].ancilla) == (a, b, anc), k
        assert p[k].status is c.CircuitStatus.COMPILED and "arxiv_2611_08161" in p[k].src, k
    assert {k: p[k].n_rot for k in rows} == {"U_inv": 0, "U_mul": 0, "U_Tr": 15, "U_FFT": 426, "U_phi": 10}
    assert {k: p[k].toffoli for k in rows} == {"U_inv": 295, "U_mul": 364, "U_Tr": 4398, "U_FFT": 1382, "U_phi": 34}
    assert g.has_fft and "U_F" not in p
    assert "arxiv_2611_08161" in g.note and "in preparation" not in g.note.lower()


@pytest.mark.parametrize("eps", EPS)
def test_table_11_as_printed(eps):
    p, L = G.GROUPS[S].primitives, math.log2(1 / eps)
    assert p["U_inv"].t(eps) == 2065 and p["U_mul"].t(eps) == 2548
    assert p["U_Tr"].t(eps) == pytest.approx(30924 + 17.25 * L, rel=1e-12)
    assert p["U_FFT"].t(eps) == pytest.approx(13597.2 + 489.9 * L, rel=1e-12)
    assert p["U_phi"].t(eps) == pytest.approx(330 + 11.5 * L, rel=1e-12)


def _compose(ham, d, eps, papers, table=True):
    """table=True composes with the printed U_mul multiplicity (Table 13), which the closed forms use; table=False
    with the model's (PRIMCOST, half the printed value since card CIRC1118-1)."""
    g, m = G.GROUPS[S], dict(G.PRIMCOST[ham])
    if table:
        m["U_mul"] = G.PRIMCOST_TABLE_U_MUL[ham]
    price = (lambda pc: pc.t_papers(eps)) if papers else (lambda pc: pc.t(eps))
    t = sum(m[k](d) * price(g.primitives[k]) for k in G.MAGNETIC)
    t += m["U_F"](d) * price(g.primitives["U_FFT"]) + g.phi_mult[ham] * price(g.primitives["U_phi"])
    return t


@pytest.mark.parametrize("ham", ["KS", "I"])
@pytest.mark.parametrize("d", [2, 3])
def test_closed_form_is_the_composition(ham, d):
    """Table 14 = Table 13 x Table 11, exactly, in the slope-only record and at the report convention. The model
    counts half the printed U_mul (card CIRC1118-1), so its per-link cost is the closed form less that half."""
    mul = G.GROUPS[S].primitives["U_mul"]
    dm = G.PRIMCOST_TABLE_U_MUL[ham](d) - G.PRIMCOST[ham]["U_mul"](d)
    assert dm == pytest.approx(G.PRIMCOST[ham]["U_mul"](d))
    for eps in EPS:
        assert _compose(ham, d, eps, papers=True) == pytest.approx(G.c_t_stated(S, ham, d, eps), rel=1e-12)
        split = G.magnetic_per_link(S, ham, d, eps, papers=True) + G.electric_per_link(S, ham, d, eps, papers=True)
        assert split + dm * mul.t_papers(eps) == pytest.approx(G.c_t_stated(S, ham, d, eps), rel=1e-12)
        # report convention (E20): the same plus 9.2 T per rotation
        full = G.magnetic_per_link(S, ham, d, eps) + G.electric_per_link(S, ham, d, eps)
        assert full + dm * mul.t(eps) == pytest.approx(G.c_t_full(S, ham, d, eps), rel=1e-12)
        assert full == pytest.approx(_compose(ham, d, eps, papers=False, table=False), rel=1e-12)


@pytest.mark.parametrize("eps", EPS)
@pytest.mark.parametrize("d", [2, 3])
def test_closed_forms_as_printed(d, eps):
    """c_t_full is Table 14's printed expression, with the 9.2 T per rotation already in its constants."""
    L = math.log2(1 / eps)
    assert G.c_t_full(S, "KS", d, eps) == pytest.approx(36945 * d - 9420.6 + (8.625 * d + 982.675) * L, rel=1e-12)
    assert G.c_t_full(S, "I", d, eps) == pytest.approx(135349 * d - 65978.2 + (25.875 * d + 1956.725) * L, rel=1e-12)
    # rotations per link per step = Table 14's synthesis-error coefficient
    assert G.c_t_rotations(S, "KS", d) == pytest.approx(7.5 * d + 854.5, rel=1e-12)
    assert G.c_t_rotations(S, "I", d) == pytest.approx(22.5 * d + 1701.5, rel=1e-12)


@pytest.mark.parametrize("ham,printed,unrounded", [("KS", 2.3e10, 2.3313e10), ("I", 6.8e10, 6.7659e10)])
def test_fiducial_totals(ham, printed, unrounded):
    """Table 14: d=3, L=10, N_t=50, eps_T=1e-8 split evenly over every R_Z (no 1/2 factor)."""
    d, ln = 3, 3 * 10 ** 3 * 50
    eps = 1e-8 / (G.c_t_rotations(S, ham, d) * ln)
    total = G.c_t_full(S, ham, d, eps) * ln
    assert total == pytest.approx(unrounded, rel=1e-3)
    assert float(f"{total:.1e}") == printed == G.GROUPS[S].benchmark[ham]


def test_phi_multiplicity_leaves_other_groups_unchanged():
    for k in ("2T", "2O", "S36x3", "S72x3", "Z3"):
        assert G.GROUPS[k].phi_mult == {}
    assert G.GROUPS[S].phi_mult == {"KS": 1, "I": 2}
