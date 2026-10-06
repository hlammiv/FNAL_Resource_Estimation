"""Sigma(216x3) group entry (Gustafson_S648_inprep, in preparation; draft mains648.tex, read 2026-10-05).

Checks the entry against the draft's printed numbers: tab:qubitcosts (primitives), tab:primcost (multiplicities,
U_Ph at 1/2), tab:costsummary (closed forms and synthesis-error coefficients) and tab:fiducial (2.0e10, 6.1e10).
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
    rows = {"U_inv": (2065, 0, 4), "U_mul": (2548, 0, 5), "U_Tr": (30786, 16.1, 10),
            "U_FFT": (1086, 529, 1), "U_phi": (630, 11.5, 1)}
    assert set(p) == set(rows)
    for k, (a, b, anc) in rows.items():
        assert (p[k].t_const, p[k].t_log, p[k].ancilla) == (a, b, anc), k
        assert p[k].status is c.CircuitStatus.COMPILED and "Gustafson_S648_inprep" in p[k].src, k
    assert {k: p[k].n_rot for k in rows} == {"U_inv": 0, "U_mul": 0, "U_Tr": 14, "U_FFT": 460, "U_phi": 10}
    assert {k: p[k].toffoli for k in ("U_inv", "U_mul", "U_Tr")} == {"U_inv": 295, "U_mul": 364, "U_Tr": 4398}
    assert g.has_fft and "U_F" not in p
    assert "in preparation" in g.note.lower()


def _compose(ham, d, eps, papers):
    g, m = G.GROUPS[S], G.PRIMCOST[ham]
    price = (lambda pc: pc.t_papers(eps)) if papers else (lambda pc: pc.t(eps))
    t = sum(m[k](d) * price(g.primitives[k]) for k in G.MAGNETIC)
    t += m["U_F"](d) * price(g.primitives["U_FFT"]) + g.phi_mult[ham] * price(g.primitives["U_phi"])
    return t


@pytest.mark.parametrize("ham", ["KS", "I"])
@pytest.mark.parametrize("d", [2, 3])
def test_closed_form_is_the_composition(ham, d):
    """tab:costsummary = tab:primcost x tab:qubitcosts, exactly (slope-only, the draft's convention)."""
    for eps in EPS:
        assert _compose(ham, d, eps, papers=True) == pytest.approx(G.c_t_stated(S, ham, d, eps), rel=1e-12)
        split = G.magnetic_per_link(S, ham, d, eps, papers=True) + G.electric_per_link(S, ham, d, eps, papers=True)
        assert split == pytest.approx(G.c_t_stated(S, ham, d, eps), rel=1e-12)
        # report convention (E20): the same plus 9.2 T per rotation
        full = G.magnetic_per_link(S, ham, d, eps) + G.electric_per_link(S, ham, d, eps)
        assert full == pytest.approx(G.c_t_full(S, ham, d, eps), rel=1e-12)
        assert full == pytest.approx(_compose(ham, d, eps, papers=False), rel=1e-12)


def test_closed_forms_at_d3_as_printed():
    L = math.log2(1e10)
    assert G.c_t_stated(S, "KS", 3, 1e-10) == pytest.approx(36876 * 3 - 34074 + (8.05 * 3 + 1061.45) * L, rel=1e-14)
    assert G.c_t_stated(S, "I", 3, 1e-10) == pytest.approx(135142 * 3 - 115216 + (24.15 * 3 + 2114.85) * L, rel=1e-14)
    # rotations per link per step = the draft's synthesis-error coefficient (7d + 923), (21d + 1839)
    for d in (2, 3):
        assert G.c_t_rotations(S, "KS", d) == pytest.approx(7 * d + 923, rel=1e-12)
        assert G.c_t_rotations(S, "I", d) == pytest.approx(21 * d + 1839, rel=1e-12)


@pytest.mark.parametrize("ham,printed,unrounded", [("KS", 2.0e10, 2.022e10), ("I", 6.1e10, 6.147e10)])
def test_fiducial_totals(ham, printed, unrounded):
    """tab:fiducial: d=3, L=10, N_t=50, eps_T=1e-8 split evenly over every R_Z (no 1/2 factor)."""
    d, ln = 3, 3 * 10 ** 3 * 50
    eps = 1e-8 / (G.c_t_rotations(S, ham, d) * ln)
    total = G.c_t_stated(S, ham, d, eps) * ln
    assert total == pytest.approx(unrounded, rel=1e-3)
    assert float(f"{total:.1e}") == printed == G.GROUPS[S].benchmark[ham]


def test_phi_multiplicity_leaves_other_groups_unchanged():
    for k in ("2T", "2O", "S36x3", "S72x3", "Z3"):
        assert G.GROUPS[k].phi_mult == {}
    assert G.GROUPS[S].phi_mult == {"KS": 1, "I": 2}
