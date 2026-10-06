"""Ch. 4 (hybrid LQCD): every named number the prose states, checked against the model.

Tolerances are stated per assertion. Since ruling R4 (2026-09-28) the chapter carries exact
values and rounds the printed result once. Since ruling R7 (2026-09-28) the 2028 demonstrator is
the FREE-FERMION block encoding of arXiv:2407.13080 at V = 4^4: the tests below build the
matrix with plain Python lists, independently of the model, and pin the column register, the
sparsity, the spectrum and log det W, lambda_min(W/s), the 18 / 34 register and the per-shot range.
Since ruling 'SHIFT + 1%' (2026-09-28) the QET target is the paper's shifted one and the tolerance
is the one that guarantees 1% on log det W: d = 84, 2.0e4-4.1e4 T per shot (at the fixed 1e-4), 1.4e6 shots.
Since ruling R-TOL (2026-09-29) the synthesis error is 1e-2 per shot and each end of the band sets
eps_rot = sqrt(1e-2/N_rot) from its own rotations. R-TOL does not count Clifford+T-exact rotations; at
s = 2(m0^2+2K^2) the O_A angle theta_0 = 0, so both ends hold 338 rotations: 1.8e4-3.4e4 T per shot.
Since ruling SYN-BUDGET (2026-09-29) EPS_SYN = 1e-2 stays the report-wide default, but the 2028 circuit
sets eps_syn = 0.1 x the polynomial tolerance 1.784e-3 = 1.784e-4 per shot: eps_rot = 7.3e-4, 21.2 T per
rotation, 1.9e4-3.5e4 T per shot; synthesis adds 0.1% of log det W to the 2% combined worst case.
The old strict xfail test_2028_neighbor_register_matches_box ('box 4 vs equation 6') is
superseded by R7 and removed: 4 is the free circuit's column register, 6 the gauged one's.
Since round R11 (2026-09-30, apply_log/r11_ch04.md) d_inv is the LP minimum kappa ln(1/eps_rel) at
eps_rel = e^-5 = 0.674% of the peak of 1/x: 500 at kappa = 1e2, 5e3 at kappa = 1e3 (was 1e4), so the
co-design instance is ~1e12 T (was 2e12); the per-shot reduction deliverable is ~10x (was 'cut d_inv by
1e2x') and 48^3x96 needs ~2e2x against 6-20x of levers. Every rung prices Hadamard-test sampling.
Since R11b (2026-09-30, apply_log/r11b_ch04.md) the Ch. 4 lever (ii), '~2-3x lower QSP degree', is
dropped (d = 84 is at the LP floor): the gauged 4^4 gets 6-20x against 27.5x and fits at no lever setting.
Since ruling E27 (H. Lamm, 2026-10-02, apply_log/r23_ch04.md) the report assumes ONE quantum machine: every
wall time is serial, shots x (T x 1 us + 0.1 ms). The 1000-LQ campaign is 27.0 yr = 5.4x the 5-year horizon
(18 configurations fit); the 1e4-LQ campaign is 6.3e3x the horizon after the 10.3x per-shot reduction. No
machine count, 'device count' or shot-parallel language may appear in the chapter or the model.
Since rulings R1, R3, R7, R9, R10 (H. Lamm, 2026-10-02, apply_log/r25_ch04.md): the Tr M^-1 rungs block-encode M at
D V log2 V (83 LQ, 2.1e8 T at 2033; 90 LQ, 2.6e11 T at 1e4 LQ); the 2033 readout is the block-encoding success
probability, with a first result (1 configuration, 30%) and a campaign (3 configurations, 20%); the Banks-Casher
mode number is the 2033 demonstration row; every wall is depth-bound where F* < 10 (2028: 74 ms per shot, 6.5 h).
"""

import cmath
import importlib.util
import itertools
import math
import sys
from dataclasses import replace
from pathlib import Path

import pytest

from estimates import common as c
from estimates import ch04_hybrid_lqcd as m

A = m.Assumptions()
R28 = m.model(A, "2028")
R33 = m.model(A, "2033")
RCD = m.model(A, "codesign")
I28, I33, ICD = R28.intermediates, R33.intermediates, RCD.intermediates

M0, K, D = 0.4, 1.0, 4


def near(x, y, tol):
    return math.isclose(x, y, rel_tol=tol)


# --- the free staggered matrix, built here with plain lists (arXiv:2407.13080 Eq. 16) ----------

def staggered_M(L, m0=M0, K=K, antiperiodic_t=False):
    """Sparse rows {column: value} of M_mn = (K/2) sum_mu eta_mu(m) (delta_{n,m+mu} - delta_{n,m-mu})
    + m0 delta_{nm}, eta_mu(m) = (-1)^(m_1 + ... + m_{mu-1}). Returns (sites, rows)."""
    sites = list(itertools.product(range(L), repeat=D))
    index = {s: i for i, s in enumerate(sites)}
    rows = []
    for s in sites:
        row = {index[s]: m0}
        for mu in range(D):
            eta = -1 if sum(s[:mu]) % 2 else 1
            for hop in (+1, -1):
                t = list(s)
                t[mu] += hop
                wraps = t[mu] < 0 or t[mu] >= L
                bc = -1.0 if (wraps and antiperiodic_t and mu == D - 1) else 1.0
                t[mu] %= L
                j = index[tuple(t)]
                row[j] = row.get(j, 0.0) + hop * 0.5 * K * eta * bc
        rows.append({j: v for j, v in row.items() if abs(v) > 1e-15})
    return sites, rows


def mdag_m(rows):
    """W = M^T M for real M, as sparse rows: W[a][c] = sum_b M[b][a] M[b][c]."""
    w = [dict() for _ in rows]
    for row in rows:
        for a, va in row.items():
            for cc, vc in row.items():
                w[a][cc] = w[a].get(cc, 0.0) + va * vc
    return [{j: v for j, v in r.items() if abs(v) > 1e-12} for r in w]


SITES4, M4 = staggered_M(4)
W4 = mdag_m(M4)


def test_2028_column_register_and_sparsity():
    # arXiv:2407.13080 Sec. II.2: 'only nine nonzero entries in each row/column' = 2D+1; the column
    # register l is 4 qubits (l = 0-7 hops, 8 <= l < 16 diagonal). app12:92,111,119.
    assert I28["nnz_per_column"] == 9 == 2 * D + 1 == I28["nnz_formula_2D_plus_1"]
    assert I28["column_qubits"] == 4 == math.ceil(math.log2(9)) == I28["column_qubits_formula"]
    assert 2 ** I28["column_qubits"] >= I28["nnz_per_column"]
    # the paper's 'six qubits' is the U(1)-gauged matrix, 33 nonzeros (Sec. II.3); Eq. Nq prints that one
    assert I28["column_qubits_gauged"] == 6 == math.ceil(math.log2(2 * D * D + 1))
    assert I33["neighbor_qubits_M_dag_M"] == 6 and I33["neighbor_qubits"] == 4   # R3: the 1000-LQ box encodes M


def test_2028_W_is_eq24_with_nearest_neighbours_cancelled():
    # Eq. (24): W = (m0^2 + 2K^2) delta - (K^2/4) sum_mu (delta_{a+2mu} + delta_{a-2mu}). Built from M here,
    # compared entry by entry with the model's Eq. (24) rows. On L = 4 the +2mu and -2mu sites coincide,
    # so each row holds 5 distinct sites (diagonal 2.16, four entries of -2 x K^2/4 = -0.5), not 9.
    index = {s: i for i, s in enumerate(SITES4)}
    assert len(W4) == 256 == I28["V"]
    for s, row in zip(SITES4, W4):
        eq24 = {index[t]: v for t, v in m.free_W_row(s, 4, M0, K).items()}
        assert row.keys() == eq24.keys()
        assert all(near(row[j], eq24[j], 1e-12) for j in row)
        assert len(row) == 5 == I28["nnz_distinct_sites_this_L"]
        assert near(row[index[s]], 2.16, 1e-12) and near(row[index[s]], I28["W_diagonal"], 1e-12)
        assert sorted(v for j, v in row.items() if j != index[s]) == pytest.approx([-0.5] * 4)
    assert I28["W_hop"] == -0.25
    # on a lattice where the two-hop sites are distinct the row has the paper's nine entries
    assert len(m.free_W_row((0,) * D, 8, M0, K)) == 9


def test_2028_W_spectrum_and_logdet():
    # app12:26,111,115: eigenvalues 0.16, 1.16, 2.16, 3.16, 4.16 with multiplicities 16, 64, 96, 64, 16;
    # log det W = 150.5528 (the validation target, printed 150.55). Plane waves are a complete
    # orthonormal set, so checking W v = lambda v on all 256 of them fixes the whole spectrum.
    spec = {}
    for n in itertools.product(range(4), repeat=D):
        v = [cmath.exp(2j * math.pi * sum(a * b for a, b in zip(n, s)) / 4) for s in SITES4]
        lam = M0 ** 2 + K ** 2 * sum(math.sin(2 * math.pi * a / 4) ** 2 for a in n)
        for i, row in enumerate(W4):
            wv = sum(val * v[j] for j, val in row.items())
            assert abs(wv - lam * v[i]) < 1e-12
        key = round(lam, 9)
        spec[key] = spec.get(key, 0) + 1
    assert spec == {0.16: 16, 1.16: 64, 2.16: 96, 3.16: 64, 4.16: 16}
    assert {round(k, 9): n for k, n in I28["W_spectrum"].items()} == spec
    logdet = sum(n * math.log(lam) for lam, n in spec.items())
    assert near(logdet, 150.5528, 1e-6)
    assert near(I28["logdet_W"], logdet, 1e-12)
    assert near(I28["lambda_min_W"], 0.16, 1e-12) and near(I28["lambda_max_W"], 4.16, 1e-12)
    assert near(I28["kappa_W"], 26, 1e-12)                          # app12:26 'kappa = 26'
    # the paper's bound (App. B, Eq. B11): m0^2 <= lambda <= m0^2 + 16 K^2
    assert M0 ** 2 <= I28["lambda_min_W"] + 1e-12 and I28["lambda_max_W"] <= M0 ** 2 + 16 * K ** 2


def test_2028_spectrum_dense_crosscheck():
    np = pytest.importorskip("numpy")
    w = np.zeros((256, 256))
    for i, row in enumerate(W4):
        for j, v in row.items():
            w[i, j] = v
    ev = np.linalg.eigvalsh(w)
    vals, counts = np.unique(np.round(ev, 9), return_counts=True)
    assert dict(zip(vals.tolist(), counts.tolist())) == {0.16: 16, 1.16: 64, 2.16: 96, 3.16: 64, 4.16: 16}
    assert near(np.linalg.slogdet(w)[1], I28["logdet_W"], 1e-10)


def test_2x2x2x2_is_degenerate():
    # app12:26: 'The smaller V=2^4 periodic lattice is degenerate': +mu and -mu neighbours coincide, the
    # hopping cancels, M = m0 x identity and W = 0.16 x identity. Antiperiodic in t: W = 1.16 x identity.
    _, m2 = staggered_M(2)
    assert all(row == {i: M0} for i, row in enumerate(m2))
    w2 = mdag_m(m2)
    assert all(list(row) == [i] and near(row[i], 0.16, 1e-12) for i, row in enumerate(w2))
    _, m2a = staggered_M(2, antiperiodic_t=True)
    assert any(len(row) > 1 for row in m2a)
    w2a = mdag_m(m2a)
    assert all(list(row) == [i] and near(row[i], 1.16, 1e-12) for i, row in enumerate(w2a))
    # the model says the same, and prices the paper's O_c as the identity there
    r = m.model(replace(A, L_2028=c.Assumed(2, "x")), "2028")
    assert r.intermediates["degenerate_M_is_m0_identity"]
    assert r.intermediates["W_spectrum"] == {0.16: 16}
    assert r.intermediates["oc_toffolis_per_query"] == (0, 0)
    assert not I28["degenerate_M_is_m0_identity"]


def test_2028_subnormalization_and_lambda_min_seen_by_the_qet():
    # arXiv:2407.13080 Eqs. 25-26: cos(theta_0/2) = 2(m0^2+2K^2)/s, cos(theta_1/2) = -4K^2/s, so
    # s >= 2(m0^2+2K^2) = 4.32 and s >= 4. app12:26,112: s = 4.32, lambda_min/s = 0.037.
    s = I28["be_subnormalization_s"]
    assert near(s, 2 * (M0 ** 2 + 2 * K ** 2), 1e-12) and near(s, 4.32, 1e-12)
    assert s >= 4 * K ** 2
    assert abs(2 * (M0 ** 2 + 2 * K ** 2) / s) <= 1 and abs(-4 * K ** 2 / s) <= 1      # the arccos exist
    assert near(I28["lambda_min_over_s"], 0.16 / 4.32, 1e-12)
    assert near(I28["lambda_min_over_s"], 0.037, 0.002)              # printed 0.037
    assert I28["lambda_max_over_s"] <= 1
    # log det W = sum log(lambda/s) + V log s; the V log s term is added back classically (app12:115)
    assert near(I28["V_log_s"], 256 * math.log(4.32), 1e-12)
    assert near(I28["sum_log_lambda_over_s"] + I28["V_log_s"], I28["logdet_W"], 1e-12)
    assert near(I28["log_s_over_lambda_min"], math.log(27), 1e-12)   # 3.30


def test_2028_register_18_and_34():
    # app12:119: '18: 8 site + 4 column index + 1 block-encoding flag + 1 QET signal + 3 multi-control
    # workspace + 1 Hadamard-test ancilla'; app12:120: QME upgrade path 34; app12:92: overhead 6, '18 LQ'
    parts = [I28[k] for k in ("site_qubits", "column_qubits", "be_flag_qubits", "qet_signal_qubits",
                              "mcx_workspace_qubits", "hadamard_qubits")]
    assert parts == [8, 4, 1, 1, 3, 1]
    assert sum(parts) == 18 == I28["lq_hadamard_route"]
    assert R28.lq == (18, 18)
    assert I28["overhead_V_independent"] == 6
    # QME (the paper's trace step, Sec. III): a copy of the site register, m = 8 QPE, one sign qubit a,
    # in place of the Hadamard-test ancilla
    assert I28["qme_purification_qubits"] == 8 and I28["n_qpe"] == 8 and I28["qme_sign_qubits"] == 1
    assert I28["lq_qme_route"] == 18 - 1 + 8 + 8 + 1 == 34
    # the workspace is what the 5-control projector ladder needs
    assert I28["projector_controls"] == 5
    assert I28["projector_workspace_needed"] == 3 == I28["mcx_workspace_qubits"]
    assert R28.lq[0] < 150                                           # 'below the 150-LQ 2028 floor'


# R-TOL (H. Lamm, 2026-09-29): each end of the band is its own circuit and takes eps_rot = sqrt(1e-2/N_rot)
# and O_A holds 2 synthesized rotations at both ends: R_y(theta_0 = 0) is the identity (see test below)
N_LO, N_HI = 84 * 2 + 85 * 2, 84 * 2 + 85 * 2                          # 338, 338
TOL28 = 0.01 * 150.5528171451149 / (256 * math.log(27))                  # 1.784e-3, polynomial tolerance
EPS_SYN28 = 0.1 * TOL28                                                  # 1.784e-4 (SYN-BUDGET)
E_LO, E_HI = c.eps_rot_for(N_LO, EPS_SYN28), c.eps_rot_for(N_HI, EPS_SYN28)   # 7.266e-4 both
T_LO, T_HI = c.t_per_rotation(E_LO), c.t_per_rotation(E_HI)             # 21.191 both (full RUS fit)


def test_2028_oa_theta0_is_identity():
    # arXiv:2407.13080 Eqs. 25-26 at the chapter's s = 2(m0^2+2K^2) = 4.32 (app12:111):
    # cos(theta_0/2) = 2(m0^2 + D K^2/2)/s = 1 -> theta_0 = 0; cos(theta_1/2) = -2^4 (K^2/4)/s = -4/4.32
    s_ = 2 * (0.4 ** 2 + 2 * 1.0 ** 2)
    assert near(s_, 4.32, 1e-12) and near(I28["be_subnormalization_s"], s_, 1e-12)
    th0 = 2 * math.acos(min(1.0, 2 * (0.16 + 2.0) / s_))
    th1 = 2 * math.acos(-16 * 0.25 / s_)
    assert th0 == 0.0 and I28["oa_theta0"] == 0.0
    assert near(I28["oa_theta1"], th1, 1e-12) and near(th1, 5.5086, 1e-4)
    # theta_1/2 and +-theta_1/2 are not multiples of pi/4: synthesized
    assert abs(th1 / 2 / (math.pi / 4) - round(th1 / 2 / (math.pi / 4))) > 0.1
    assert I28["oa_rotations_generic_angles"] == (2, 4)
    assert I28["oa_rotations_per_query"] == (2, 2)
    assert m.oa_synthesized_rotations(0.0, th1) == (2, 2)
    assert m.oa_synthesized_rotations(1.0, th1) == (2, 4)                # generic angles recover (2, 4)
    assert m.oa_synthesized_rotations(0.0, 0.0) == (0, 0)


def test_2028_rtol_tolerance_per_circuit():
    # app12:125,136 (SYN-BUDGET): 'total synthesis error 1.8e-4 per shot, set below the polynomial tolerance';
    # 'eps_rot = 7.3e-4 and 21.2 T per rotation'. The report-wide default is unchanged.
    assert c.EPS_SYN == 1e-2 and I28["eps_syn_report_default"] == 1e-2
    assert A.synth_to_poly_2028.lo == 0.1 and A.synth_to_poly_2028.prov is c.Provenance.STATED
    assert not hasattr(A, "eps_syn_2028")
    assert near(I28["eps_syn"], EPS_SYN28, 1e-12) and near(I28["eps_syn"], 0.1 * I28["poly_tolerance"], 1e-12)
    assert near(I28["eps_syn"], 1.784e-4, 3e-4) and near(I28["eps_syn"], 1.8e-4, 0.01)
    assert I28["rotations_per_shot"] == (N_LO, N_HI) == (338, 338)
    assert all(near(x, y, 1e-12) for x, y in zip(I28["eps_rot"], (E_LO, E_HI)))
    assert near(E_LO, 7.266e-4, 1e-3) and E_HI == E_LO
    assert near(E_LO, 7.3e-4, 0.005)
    assert all(near(x, y, 1e-12) for x, y in zip(I28["t_per_rotation"], (T_LO, T_HI)))
    assert near(T_LO, 21.191, 1e-4) and T_HI == T_LO
    assert round(T_LO, 1) == 21.2
    assert near(T_LO, 1.15 * math.log2(1 / E_LO) + 9.2, 1e-12)          # the chapter keeps the full RUS fit
    # no fixed tolerance is left in the 2028 model
    assert not hasattr(A, "rot_eps_2028")


def test_2028_per_query_t():
    # app12:124-126: '84 BE queries at 126--322 T: shifts O_c in 12--40 Toffolis, entries O_A in 2
    # synthesized rotations'; '7 T per Toffoli; 21.2 T per rotation'
    assert I28["oc_toffolis_per_query"] == (4 * c.mcx_toffoli_count(3), 8 * c.mcx_toffoli_count(4)) == (12, 40)
    assert I28["oc_t_per_query"] == (84, 280)
    assert I28["oa_rotations_per_query"] == (2, 2)
    assert I28["t_per_toffoli"] == 7
    lo, hi = I28["t_per_query"]
    assert near(lo, 84 + 2 * T_LO, 1e-12) and round(lo) == 126
    assert near(hi, 280 + 2 * T_HI, 1e-12) and round(hi) == 322
    # app12:125: '85 controlled projector phases at 98 T: 8 Toffolis + 2 synthesized rotations'
    assert I28["projector_phases"] == 85
    plo, phi = I28["t_per_projector_phase"]
    assert near(plo, 56 + 2 * T_LO, 1e-12) and near(phi, 56 + 2 * T_HI, 1e-12)
    assert round(plo) == round(phi) == 98


def test_2028_per_shot_range():
    # app12:27,123,127: '1.9--3.5e4 T-gates', '2.8--5.3x margin' against the 1e5 cap (SYN-BUDGET)
    lo = 84 * (84 + 2 * T_LO) + 85 * (56 + 2 * T_LO)                 # 18978.4
    hi = 84 * (280 + 2 * T_HI) + 85 * (56 + 2 * T_HI)                # 35442.4
    assert I28["be_queries"] == 84 == I28["d"]
    assert near(R28.hard_ops[0], lo, 1e-12) and near(R28.hard_ops[1], hi, 1e-12)
    assert near(R28.hard_ops[0], 18978.4, 1e-5) and near(R28.hard_ops[1], 35442.4, 1e-5)
    assert near(R28.hard_ops[0], 1.9e4, 0.002) and near(R28.hard_ops[1], 3.5e4, 0.013)
    assert I28["t_per_shot"] == R28.hard_ops
    assert round(I28["margin_vs_cap"][0], 1) == 2.8 and round(I28["margin_vs_cap"][1], 1) == 5.3
    assert R28.hard_ops[1] < A.t_cap_2028.lo
    # the breakdown is the literal end (the circuit as the paper draws it) and sums to it
    assert near(R28.breakdown_total(), R28.hard_ops[1], 1e-12)
    assert [p.name for p in R28.breakdown] == ["O_c_controlled_shift_toffolis", "O_A_entry_rotations",
                                               "projector_phase_toffolis", "projector_phase_rotations"]
    assert [p.count for p in R28.breakdown] == [84 * 40, 84 * 2, 85 * 8, 85 * 2]
    assert [p.status for p in R28.breakdown] == [c.CircuitStatus.COMPILED, c.CircuitStatus.COMPILED,
                                                 c.CircuitStatus.SCALING, c.CircuitStatus.SCALING]
    assert all("ASSUMED" in p.note or "does not draw" in p.note or "Not drawn" in p.note for p in R28.breakdown)


def test_2028_rotations_and_synthesis_error():
    # app12:137 (SYN-BUDGET): 'the total synthesis error is 0.1x the tolerance, 1.8e-4 per shot ... On
    # log det W this is 0.15 against the polynomial bound of 1.5, or 0.1%, and the combined worst case stays
    # at ~2%.' app12:118-119: 'synthesis error: <= 1.8e-4 on the normalized trace, 0.1x the polynomial
    # tolerance (0.15 absolute, 0.1% of log det W)'; 'combined worst case ~2% (1%+1%+0.1%)'.
    # Report convention (T4; common.eps_per_rotation, "incoherent"): N eps^2.
    assert I28["rotations_per_shot"] == (84 * 2 + 85 * 2,) * 2 == (338, 338)
    lo, hi = I28["synthesis_error_randomized"]
    assert near(lo, EPS_SYN28, 1e-12) and near(hi, EPS_SYN28, 1e-12)
    assert near(c.eps_per_rotation(hi, 338, "incoherent"), E_HI, 1e-12)
    assert near(I28["synthesis_error_over_poly_tolerance"], 0.1, 1e-12)
    assert lo < I28["poly_tolerance"]                                # set below the polynomial tolerance
    # same units as the polynomial bound: x V log(s/lambda_min) on log det W
    assert near(I28["synthesis_abs_error"], 256 * math.log(27) * EPS_SYN28, 1e-12)
    assert near(I28["synthesis_abs_error"], 0.1 * I28["poly_abs_error_bound"], 1e-12)
    assert round(I28["synthesis_abs_error"], 2) == 0.15 and round(I28["poly_abs_error_bound"], 1) == 1.5
    assert near(I28["synthesis_rel_error"], 1e-3, 1e-12)             # 0.1% of log det W
    assert near(I28["combined_worst_case_with_synthesis_rel_error"], 0.021, 1e-12)
    assert round(100 * I28["combined_worst_case_with_synthesis_rel_error"]) == 2      # printed '~2%'
    # coherent errors would add to N eps = 0.25: randomized synthesis is still the model
    assert near(I28["synthesis_error_if_coherent"][1], 338 * E_HI, 1e-12)
    assert I28["synthesis_error_if_coherent"][0] > I28["poly_tolerance"]
    assert "synthesis_error_quadrature" not in I28                   # one convention only


def test_2028_shifted_target_and_add_backs():
    # arXiv:2407.13080 Sec. IV (txt:697-703): 'A shift of 1/2 works well. This known shift can always be
    # subtracted at the end'. app12:113: QET of log|x|/log(s/lambda_min) + 1/2; app12:115: 'classically,
    # the 1/2 shift is subtracted and the V log s term is added back'.
    assert A.poly_shift_2028.lo == 0.5 and A.poly_shift_2028.prov is c.Provenance.CITED
    assert I28["poly_shift"] == 0.5
    # the shifted target stays inside [-1/2, 1/2] on [lambda_min/s, 1], so |p| <= 1 has room on both sides
    lam, ln = I28["lambda_min_over_s"], I28["log_s_over_lambda_min"]
    assert near(math.log(lam) / ln + 0.5, -0.5, 1e-12) and math.log(1.0) / ln + 0.5 == 0.5
    # what the device estimates is the mean of the shifted target over the spectrum, 0.2345
    mean = sum(n * (math.log(k / 4.32) / ln + 0.5) for k, n in I28["W_spectrum"].items()) / 256
    assert near(I28["normalized_trace_shifted"], mean, 1e-9) and abs(mean) <= 0.5
    # and log det W is recovered with BOTH classical terms
    assert near(I28["shift_term_subtracted"], 256 * ln * 0.5, 1e-12)             # 421.87
    rebuilt = 256 * ln * mean - I28["shift_term_subtracted"] + I28["V_log_s"]
    assert near(rebuilt, 150.5528, 1e-6) and near(rebuilt, I28["logdet_W"], 1e-9)


def test_2028_accuracy_statement():
    # app12:116-118, two kinds of claim. Polynomial: 'guaranteed <= 1% of log det W by the uniform
    # tolerance 1.78e-3 ... (V log(s/lambda_min) x 1.78e-3 = 1.5 absolute)'. Sampling (R16): '1% rms
    # (eps = 1.78e-3 on the normalized trace)'. 'combined ~2%'.
    assert I28["target_rel_error"] == 0.01
    assert near(I28["poly_tolerance"], 0.01 * 150.5528 / (256 * math.log(27)), 1e-6)
    assert near(I28["poly_tolerance"], 1.784e-3, 3e-4) and near(I28["poly_tolerance"], 1.78e-3, 0.003)
    assert near(I28["poly_abs_error_bound"], 256 * math.log(27) * I28["poly_tolerance"], 1e-12)
    assert near(I28["poly_abs_error_bound"], 1.5055, 1e-4) and near(I28["poly_abs_error_bound"], 1.5, 0.004)
    assert near(I28["poly_rel_error_bound"], 0.01, 1e-12)
    assert near(I28["sampling_eps"], I28["poly_tolerance"], 1e-12)
    assert near(I28["sampling_rel_error"], 0.01, 1e-12) and I28["sampling_is_rms"] and I28["shot_variance"] == 1.0
    assert near(I28["combined_worst_case_rel_error"], 0.02, 1e-12)
    # d = 84 is minimal for that tolerance: the LP optimum at 84 is inside it, at 82 outside
    assert I28["d"] == 84
    assert I28["lp_error_at_d"] <= I28["poly_tolerance"] < I28["lp_error_at_d_minus_2"]
    assert I28["d_is_for_this_instance"]
    # app12:100,136: what the shift saves. 'Without the shift the same tolerance needs d=314' and
    # '0.72--1.3e5 T per shot' (SYN-BUDGET: its own circuit, 1258 rotations at both ends, same eps_syn);
    # 'at a tolerance of 1e-2 the degrees are 48 with the shift and 112 without'
    assert I28["d_unshifted"] == 314
    assert I28["rotations_per_shot_unshifted"] == (314 * 2 + 315 * 2,) * 2 == (1258, 1258)
    assert all(near(x, c.eps_rot_for(1258, EPS_SYN28), 1e-12) for x in I28["eps_rot_unshifted"])
    lo, hi = I28["t_per_shot_unshifted"]
    tl = th = c.t_per_rotation(c.eps_rot_for(1258, EPS_SYN28))
    assert near(lo, 314 * (84 + 2 * tl) + 315 * (56 + 2 * tl), 1e-12)
    assert near(hi, 314 * (280 + 2 * th) + 315 * (56 + 2 * th), 1e-12)
    assert near(lo, 0.72e5, 0.001) and near(hi, 1.3e5, 0.03)            # 7.20e4, 1.336e5
    assert hi > A.t_cap_2028.lo > R28.hard_ops[1]                    # the shift is what brings it inside
    assert I28["d_shifted_at_tolerance_1e-2"] == 48 and I28["d_unshifted_at_tolerance_1e-2"] == 112
    # LP degrees carried for the record
    assert I28["d_lp_at_lambda_min_0.16"] == 36 and I28["d_lp_at_tolerance_1e-3"] == 510
    assert I28["d_paper_fig12"] == (64, 70)


def test_2028_lp_degree_is_minimal():
    # d = 84 is 'LP-minimal degree, this work'. The LP needs numpy and scipy, which the package may
    # not import, so it lives in apply_log/ and is loaded by path; skipped where scipy is absent.
    pytest.importorskip("numpy")
    pytest.importorskip("scipy")
    path = Path(m.__file__).resolve().parent / "apply_log" / "ch04_roundB_lp_degree.py"
    spec = importlib.util.spec_from_file_location("ch04_roundB_lp_degree", path)
    lp = importlib.util.module_from_spec(spec)
    keep, sys.dont_write_bytecode = sys.dont_write_bytecode, True     # no __pycache__ in apply_log/
    try:
        spec.loader.exec_module(lp)
    finally:
        sys.dont_write_bytecode = keep
    lam, tol, shift = I28["lambda_min_over_s"], I28["poly_tolerance"], I28["poly_shift"]
    e84, e82 = lp.min_err(lam, 84, shift=shift)[0], lp.min_err(lam, 82, shift=shift)[0]
    assert e84 <= tol < e82
    assert near(e84, A.d_log_err.lo, 1e-3) and near(e82, A.d_log_err_below.lo, 1e-3)    # 1.6993e-3, 1.8660e-3
    # without the shift the round-B degree: err(112) <= 1e-2 < err(110)
    assert lp.min_err(lam, 112)[0] <= 1e-2 < lp.min_err(lam, 110)[0]
    assert lp.min_err(lam, 112)[0] > tol                             # and 112 unshifted is far from the 1%


def test_d_inv_lp_law_odd():
    # R11: d_inv = kappa ln(1/eps_rel) is the LP-minimal odd degree for 1/x on [1/kappa, 1] with |p| <= 1
    # on [-1, 1]; apply_log/ch04_r11_lp_odd.py. Pinned at kappa = 30, where the LP is fast: the degree
    # that reaches eps_rel = e^-5 is 151 (= 30 x 5, odd), and 149 does not. The same law gives 500 (501
    # odd) at kappa = 1e2, reproduced in apply_log/r11_ch04.md.
    pytest.importorskip("numpy")
    pytest.importorskip("scipy")
    path = Path(m.__file__).resolve().parent / "apply_log" / "ch04_r11_lp_odd.py"
    spec = importlib.util.spec_from_file_location("ch04_r11_lp_odd", path)
    lp = importlib.util.module_from_spec(spec)
    keep, sys.dont_write_bytecode = sys.dont_write_bytecode, True     # no __pycache__ in apply_log/
    try:
        spec.loader.exec_module(lp)
    finally:
        sys.dont_write_bytecode = keep
    tol = 0.5 * A.d_inv_eps_rel.lo                                   # peak of f = 1/(2 kappa x) is 1/2
    e151, e149, e121 = lp.min_err(30, 151), lp.min_err(30, 149), lp.min_err(30, 121)
    assert e151 <= tol < e149
    for d, e in ((151, e151), (149, e149), (121, e121)):
        assert near(e, 0.5 * math.exp(-d / 30), 0.01)                # err(d) = 0.5 e^(-d/kappa)


def test_2028_shots_wall_time_and_eps_l():
    # app12:128-131: 'eps_l <~ 2.8e-6 ... set by the 3.5e4-T end', '19--36 ms (1 us/T-gate plus 0.1 ms per
    # shot)' (E27; was '19--35 ms' at 1 us/T alone), '~3.1e5' shots
    # (R16 rms: Var/eps^2, Var = 1, eps = 1.78e-3; was ln(1/delta)/eps^2 = 1.4e6), '6.0e9--1.1e10 T-gates;
    # 1.7--3.1 h on one machine' (was 2.7--5.1e10, 8--14 h)
    eps = 0.01 * 150.5528171451149 / (256 * math.log(27))
    assert near(I28["shots"], 1.0 / eps ** 2, 1e-9)                  # 3.1408e5
    assert near(I28["shots"], 3.1e5, 0.014)
    assert R28.shots == (I28["shots"],) * 2
    assert near(I28["shots_qme_alternative"], 1.0 / eps, 1e-9)       # 560 on the QME route
    assert near(I28["shots_99pct_interval_factor"], 6.635, 1e-3)     # prose '6.6x more shots'
    assert I28["shot_overhead_s"] == 1e-4                            # report rule: ~0.1 ms per shot
    for w, t in zip(I28["wall_per_shot_s"], R28.hard_ops):
        assert near(w, t * 1e-6 + 1e-4, 1e-12)                       # T x 1 us + t0, one machine
    assert round(1e3 * I28["wall_per_shot_s"][0]) == 19 and round(1e3 * I28["wall_per_shot_s"][1]) == 36
    lo, hi = I28["wall_total_h"]
    assert near(lo, 1.6645, 1e-3) and round(lo, 1) == 1.7            # printed '1.7' (1.6557 without t0)
    assert near(hi, 3.1008, 1e-3) and round(hi, 1) == 3.1            # printed '3.1' (3.0921 without t0)
    for w, ws in zip(I28["wall_total_s"], I28["wall_per_shot_s"]):
        assert near(w, I28["shots"] * ws, 1e-12)                     # serial: shots x per-shot time
    assert R28.wall_time_s == I28["wall_total_corrected_s"]          # R10: the box quotes the depth-bound wall
    assert near(I28["t_total"][0], 6.0e9, 0.01) and near(I28["t_total"][1], 1.1e10, 0.013)
    assert near(R28.epsilon_l[0], 0.1 / R28.hard_ops[1], 1e-12)      # the literal, 3.5e4-T end binds
    assert near(R28.epsilon_l[0], 2.8e-6, 0.008)                     # 2.821e-6 printed '<~ 2.8e-6'
    assert R28.epsilon_l[0] > A.eps_l_2028.lo                        # looser than the RFI's 1e-8


def test_2028_gauged_instance_does_not_fit():
    # app12:28,94,136: gauged V=4^4 at D^2 V log2 V = 3.3e4 units per query, d = 84, one T per unit:
    # '~2.8e6 T per shot, a factor ~28x over the cap'; the 1e7 fault budget alone would admit it
    assert I28["gauged_units_D2Vlog2V"] == 16 * 256 * 8 == 32768
    assert near(I28["gauged_units_D2Vlog2V"], 3.3e4, 0.01)
    assert I28["gauged_t_per_shot"] == 84 * 32768 == 2752512
    assert near(I28["gauged_t_per_shot"], 2.8e6, 0.02)               # 2.75e6 printed '~2.8e6'
    assert near(I28["gauged_factor_over_cap"], 27.5, 0.001) and round(I28["gauged_factor_over_cap"]) == 28
    assert I28["gauged_saving_needed"] == I28["gauged_factor_over_cap"]
    assert I28["gauged_fits_fault_budget"]
    # R11b (ch04-lever-ii dropped): 'No single lever closes the gap. (i) and (iii) together give 6--20x,
    # short of the ~28x gap ... 1.4e5 T per shot, 1.4x over the cap ... at any setting of the levers.'
    assert not hasattr(A, "qsp_degree_saving")
    assert I28["saving_levers"] == (6, 20)
    assert not I28["lcu_alone_closes_gap"]
    assert not I28["levers_close_gap_at_upper_end"] and not I28["levers_close_gap_at_lower_end"]
    assert I28["saving_levers"][1] < I28["gauged_saving_needed"]      # 20 < 27.5
    assert near(I28["gauged_after_levers"][0], 84 * 32768 / 20, 1e-12)            # 137625.6
    assert near(I28["gauged_after_levers"][0], 1.4e5, 0.02)          # printed '1.4e5'
    assert near(I28["gauged_after_levers"][1], 84 * 32768 / 6, 1e-12)             # 4.59e5
    assert near(I28["gauged_over_cap_at_lever_top"], 1.376, 0.001)
    assert round(I28["gauged_over_cap_at_lever_top"], 1) == 1.4      # printed '1.4x over'
    # the free circuit is two orders cheaper per query than the gauged one
    assert I28["gauged_units_D2Vlog2V"] / I28["t_per_query"][1] > 80


def test_2028_model_moves_under_perturbation():
    # Verifier trap: nothing is pinned. L = 8 adds 4 site qubits and turns each controlled shift into a
    # controlled incrementer; the merge of +2mu with -2mu is no longer available, so lo = hi for O_c.
    r = m.model(replace(A, L_2028=c.Assumed(8, "x")), "2028")
    assert r.lq == (22, 22) and r.intermediates["lq_qme_route"] == 42
    per_shift = c.mcx_toffoli_count(4) + c.mcx_toffoli_count(5)
    assert r.intermediates["oc_toffolis_per_query"] == (8 * per_shift, 8 * per_shift)
    assert r.hard_ops[1] > R28.hard_ops[1]
    assert r.intermediates["nnz_distinct_sites_this_L"] == 9
    # lambda_min/s does not depend on L, but the tolerance that guarantees 1% follows log det W / V,
    # which does (0.588 at L = 4): the carried LP degree is flagged as not sized for this instance
    assert near(r.intermediates["lambda_min_over_s"], I28["lambda_min_over_s"], 1e-12)
    assert not near(r.intermediates["poly_tolerance"], I28["poly_tolerance"], 1e-3)
    assert not r.intermediates["d_is_for_this_instance"]
    # moving the mass moves lambda_min/s, and the carried LP degree is flagged as stale
    r = m.model(replace(A, m0_2028=c.Stated(0.2, "x")), "2028")
    assert not r.intermediates["d_is_for_this_instance"]
    # R-TOL: a lower degree is a smaller circuit, so its own tolerance loosens (N_rot = 36 x 2 + 37 x 2 = 146)
    r = m.model(replace(A, d_log_2028=c.Assumed(36, "x")), "2028")
    t = c.t_per_rotation(c.eps_rot_for(36 * 2 + 37 * 2, EPS_SYN28))
    assert r.intermediates["rotations_per_shot"][1] == 146
    assert near(r.hard_ops[1], 36 * (280 + 2 * t) + 37 * (56 + 2 * t), 1e-12)
    assert t < T_HI
    # and a looser budget fraction loosens every tolerance (4x the budget doubles eps_rot)
    r = m.model(replace(A, synth_to_poly_2028=c.Stated(0.4, "x")), "2028")
    assert near(r.intermediates["eps_rot"][1], 2 * E_HI, 1e-12) and r.hard_ops[1] < R28.hard_ops[1]
    # the budget follows the polynomial tolerance: L = 8 moves tol, so eps_syn moves with it
    r = m.model(replace(A, L_2028=c.Assumed(8, "x")), "2028")
    assert near(r.intermediates["eps_syn"], 0.1 * r.intermediates["poly_tolerance"], 1e-12)


def test_2028_depth_and_factories():
    # R9, R10 (factory.json, scheduled): 'the T-depth is 84 x 42.4 + 85 x 45.4 = 7.4e3 ... at most 2.6--4.8
    # factories ... 74 ms rather than 19--36 ms, and the demo takes 6.5 h'
    tr = I28["t_per_rotation"][0]
    assert near(I28["t_depth_per_query"], 2 * tr, 1e-12) and round(I28["t_depth_per_query"], 1) == 42.4
    assert near(I28["t_depth_per_phase"], 3 + 2 * tr, 1e-12) and round(I28["t_depth_per_phase"], 1) == 45.4
    d = 84 * 2 * tr + 85 * (3 + 2 * tr)
    assert I28["t_depth_per_shot"] == (d, d) and near(d, 7.42e3, 0.001)
    f = I28["f_star"]
    assert round(f[0], 1) == 2.6 and round(f[1], 1) == 4.8 and I28["baseline_ok"] == (False, False)
    w = I28["wall_per_shot_corrected_s"]
    assert all(near(x, d * 1e-5 + 1e-4, 1e-12) for x in w) and round(1e3 * w[0]) == 74
    h = I28["wall_total_corrected_h"]
    assert round(h[0], 1) == round(h[1], 1) == 6.5
    assert I28["wall_campaign_s"] == I28["wall_total_corrected_s"] == R28.wall_time_s
    assert I28["wall_first_result_s"] is None
    assert max(I28["factories_for_1yr"]) < 1                        # one factory suffices for a year
    # E29 verifier: on one factory a shot is N_T x 10 us, 0.19-0.35 s, and the demo 17-31 h
    one = [t * 1e-5 * I28["shots"] / 3600 for t in R28.hard_ops]
    assert round(one[0]) == 17 and round(one[1]) == 31
    tex = (Path(__file__).resolve().parents[3] / "applications" / "app12_hybrid_lqcd.tex").read_text()
    for s_ in ("$84\\times 42.4+85\\times 45.4=7.4\\times 10^{3}$", "$2.6$--$4.8$ factories",
               "$74\\,$ms rather than $19$--$36\\,$ms", "$6.5$ h on one machine per single-config demo",
               "$6.5$ h once three to five factories feed it ($17$--$31$ h on one)"):
        assert s_ in tex, s_


# --- Eq. Nq_hybrid (app12:88) and the gauged register boxes -----------------------------

class _NoPaper(str):
    """Stand-in for the chapter text when the paper tree is absent: any use skips the test."""

    def _skip(self, *a, **k):
        pytest.skip("paper source not present; chapter-text check skipped")

    __contains__ = lower = split = count = find = index = _skip


_TEX_PATH = Path(__file__).resolve().parents[3] / "applications" / "app12_hybrid_lqcd.tex"
TEX = _TEX_PATH.read_text() if _TEX_PATH.is_file() else _NoPaper()
YR = 365.25 * 86400


def test_2033_register_addends():
    # R3 (2026-10-02): the box block-encodes M, 2D+1 = 9 entries, a 4-qubit column index:
    # '83 (13 site + 4 column index + 4 color + 12 QPE + 50 ancilla)'; incoherent route ~71
    assert I33["site_qubits"] == 13
    assert I33["neighbor_qubits"] == 4 == math.ceil(math.log2(2 * D + 1))
    assert I33["neighbor_qubits_box"] == 4
    assert I33["neighbor_qubits_M_dag_M"] == 6                     # the M^dag M index the box printed before R3
    assert I33["color_qubits"] == 4                                # ceil(log2 9)
    assert I33["overhead_V_independent"] == 62                     # '62 qubits ... (12 QPE + 50 ancilla)'
    assert I33["lq_box_sum"] == 83 == I33["lq_eq"]
    assert R33.lq == (83, 83)
    assert I33["lq_incoherent_route"] == 71
    assert I33["lq_with_W_column_index"] == 85
    assert "at 83 LQ on the $V{=}8^3{\\times}16$ ensemble (13 site, 4 column index of $M$, 4 color, 62 overhead)" in TEX
    assert "needs $\\sim 71$" in TEX


def test_register_scales_as_log_V():
    # app12 overview: '83 LQ at 8^3x16 ... 86 LQ at 16^4 and 98 LQ at 128^4: 12 qubits per 4096x'
    # app12:92: 'at V=4^4 the site index needs 8 qubits, and at V=128^4 it needs 28'
    assert I33["lq_16x4"] == 86
    assert I33["lq_128x4"] == 98
    assert I33["qubits_per_4096x_volume"] == 12
    assert I28["site_qubits"] == 8
    assert I33["site_qubits_128x4"] == 28
    assert "costs 83 LQ on the coarse" in TEX and "costs 86 LQ at $V{=}16^4$ and 98 LQ" in TEX


# --- 1000-LQ box (rulings R1, R3, R7, R10; H. Lamm 2026-10-02) ------------------------------

def test_2033_per_shot_t():
    # R3: M at D V log2 V per query: 500 x 4 x 8192 x 13 = 2.13e8, printed '2.1e8' (1.4%); the M^dag M price
    # 500 x 16 x 8192 x 13 = 8.52e8 is 4x (= D) above it
    assert near(I33["scaling_units_DVlogV"], 4 * 8192 * 13, 1e-12)
    assert near(I33["scaling_units_D2VlogV"], 16 * 8192 * 13, 1e-12)
    assert near(R33.hard_ops[0], 500 * 4 * 8192 * 13, 1e-12)              # 2.12992e8
    assert near(R33.hard_ops[0], 2.1e8, 0.015)                             # printed '2.1e8'
    assert near(I33["t_per_shot_before_R3"], 8.52e8, 1e-3)
    assert I33["reprice_factor"] == 4
    assert R33.breakdown_total() == R33.hard_ops[0]
    assert I33["fits_budget"] and R33.hard_ops[0] <= A.t_cap_2033.lo
    assert near(R33.epsilon_l[0], 0.1 / R33.hard_ops[0], 1e-12)           # 4.695e-10, printed '<~ 4.7e-10'
    assert near(R33.epsilon_l[0], 4.7e-10, 0.002)
    assert "$2.1\\times 10^{8}$ T-gates (derived here)" in TEX
    assert "$\\lesssim 4.7\\times 10^{-10}$ ($3.7\\times 10^{-10}$ for the Banks--Casher row; 0.1 expected faults per shot)" in TEX
    # E29 verifier: the box eps_l covers the Banks-Casher row at d = 640 (INSTANCE_ROWS value)
    assert round(0.1 / I33["bc_t_per_shot"][1], 11) == 3.7e-10
    assert "at $a\\!\\approx\\!0.15\\,$fm, five configurations at $30\\%$:" in TEX
    # E29 item 7: one round of amplitude amplification (k=1) = 3 QET applications per shot
    aa = 3 * R33.hard_ops[0]
    assert aa <= A.t_cap_2033.lo and round(aa, -7) == 6.4e8 and round(0.1 / aa, 11) == 1.6e-10
    assert "triples the shot to $6.4\\times 10^{8}$ T" in TEX and "$\\epsilon_l\\approx 1.6\\times 10^{-10}$" in TEX
    # E29 verifier: am ~ 0.04 is near the strange mass, so m_pi ~ 600 MeV, not 300-400 MeV. Since the cut pass
    # (2026-10-06) the ensemble parameters are stated once, in the 1000-LQ box; objective 2 points at the box.
    assert "300--400" not in TEX and TEX.count("$m_\\pi\\!\\approx\\!600\\,$MeV ($am\\!\\approx\\!0.04$, near the strange mass)") == 1
    assert "500\\times 4.3\\times 10^{5}\\approx 2.1\\times 10^{8}$ T" in TEX
    assert near(I33["scaling_units_DVlogV"], 4.3e5, 0.01)


def test_2033_d_inv_is_the_lp_minimum():
    # app12 (R11): 'd_inv = 500 ... (linear program, this work)'; Eq. Ngate: d_inv ~= kappa ln(1/eps_poly).
    assert I33["d_inv"] == 500
    assert near(I33["d_inv_eps_rel"], math.exp(-5), 1e-12)
    assert near(I33["d_inv_eps_rel"], 0.0067, 0.01)                        # printed '0.67%'
    assert near(I33["d_inv_lp_law"], 100 * 5, 1e-12)
    assert near(I33["d_inv"], I33["d_inv_lp_law"], 1e-9)
    # the CKS-class kappa log(kappa/eps), kept for comparison: 400 (log10) or 921 (ln)
    assert near(I33["d_inv_formula_log10"], 400, 1e-9)
    assert near(I33["d_inv_formula_ln"], 921, 0.01)
    assert I33["d_inv_formula_log10"] < I33["d_inv"] < I33["d_inv_formula_ln"]
    assert near(I33["d_inv_implied_prefactor_ln"], 500 / 921, 0.01)
    assert near(I33["t_per_shot_d_inv_formula_ln"], 921 * 4 * 8192 * 13, 0.01)   # 3.9e8 at D V log2 V


def test_free_condensate_matches_a_direct_solve():
    # c = (1/V_F) Re Tr M^-1 for free staggered fermions, antiperiodic in time. Checked here on 4^3 x 4 by a
    # direct CG solve of M x = e_0 on the matrix built above (the diagonal of M^-1 is the same at every site).
    m0 = 0.04
    sites, rows = staggered_M(4, m0=m0, antiperiodic_t=True)
    n = len(sites)
    mt = [dict() for _ in range(n)]
    for i, r in enumerate(rows):
        for j, v in r.items():
            mt[j][i] = v
    def mv(rs, x):
        return [sum(v * x[j] for j, v in r.items()) for r in rs]
    def wv(x):
        return mv(mt, mv(rows, x))
    b = mv(mt, [1.0] + [0.0] * (n - 1))                  # W y = M^T e_0 -> y = M^-1 e_0
    x, r = [0.0] * n, b[:]
    p, rr = r[:], sum(v * v for v in r)
    for _ in range(400):
        ap = wv(p)
        al = rr / sum(u * v for u, v in zip(p, ap))
        x = [u + al * v for u, v in zip(x, p)]
        r = [u - al * v for u, v in zip(r, ap)]
        rn = sum(v * v for v in r)
        if rn < 1e-26:
            break
        p, rr = [u + rn / rr * v for u, v in zip(r, p)], rn
    assert near(x[0], m.free_condensate(4, 4, m0), 1e-8)
    # the box lattice: 0.0275 at am = 0.04 (shot audit, recomputed)
    assert near(I33["condensate_free"], 0.0275, 0.002)
    assert near(I33["kappa_from_am"], 101, 1e-12)                          # '(4+am)/am ~ 1e2'


def test_2033_p_readout_two_tiers():
    # R1, R3: P = f^2 am c, f = 1/2, am = 0.04, c = 0.0275-0.08 -> P = 2.7-8.0e-4; shots (1-P)/(P r^2)
    P = I33["P_success"]
    assert near(P[0], 0.25 * 0.04 * I33["condensate_free"], 1e-12) and near(P[1], 0.25 * 0.04 * 0.08, 1e-12)
    assert near(P[0], 2.7e-4, 0.02) and near(P[1], 8.0e-4, 1e-9)
    s1, s2 = I33["shots_first"], I33["shots_campaign"]
    assert near(s1[0], (1 - P[1]) / (P[1] * 0.09), 1e-12) and near(s1[1], (1 - P[0]) / (P[0] * 0.09), 1e-12)
    assert near(s1[0], 1.4e4, 0.01) and near(s1[1], 4.0e4, 0.011)         # '1.4--4.0e4 shots'
    assert near(s2[0], 3 * (1 - P[1]) / (P[1] * 0.04), 1e-12)
    assert near(s2[0], 0.94e5, 0.004) and near(s2[1], 2.7e5, 0.011)        # '0.94--2.7e5 shots'
    assert R33.shots == s2
    t1, t2 = I33["t_total_first"], I33["t_total_campaign"]
    assert near(t1[0], 3.0e12, 0.015) and near(t1[1], 8.6e12, 0.001)
    assert near(t2[0], 2.0e13, 0.003) and near(t2[1], 5.8e13, 0.002)
    # the Hadamard test of the same trace: mean am c / 2 = 5.5e-4 - 1.6e-3, 4.3e6 - 3.7e7 shots at 30%, 310-910x
    mu = I33["mu_hadamard"]
    assert near(mu[0], 5.5e-4, 0.001) and near(mu[1], 1.6e-3, 1e-12)
    h = I33["shots_hadamard_first"]
    assert near(h[0], 4.3e6, 0.01) and near(h[1], 3.7e7, 0.01)
    g = I33["p_readout_gain"]
    assert round(g[0], -1) == 310 and round(g[1], -1) == 910
    # the retired plan: eps = 1e-2 was 6-18x the Hadamard mean
    e = I33["retired_eps_over_mu"]
    assert round(e[0]) == 6 and round(e[1]) == 18
    for s in ("$P=2.7$--$8.0\\times 10^{-4}$", "$1.4$--$4.0\\times 10^{4}$ shots", "$310$--$910\\times$ more",
              "$4.3\\times 10^{6}$--$3.7\\times 10^{7}$ shots",
              "$0.94$--$2.7\\times 10^{5}$ shots, $2.0$--$5.8\\times 10^{13}$ T-gates",
              "$3.0$--$8.6\\times 10^{12}$ T-gates", "$am\\!\\approx\\!0.04$"):
        assert s in TEX, s


def test_2033_depth_factories_and_walls():
    # R9, R10: D_T = N_T / F*, F* = 1-4 (factory.json); a shot is max(N_T x 1 us, D_T x 10 us) + 0.1 ms = 533-2130 s
    t = R33.hard_ops[0]
    assert I33["t_depth_per_shot"] == (t / 4, t)
    assert I33["f_star"] == (1.0, 4.0)
    assert I33["baseline_ok"] == (False, False)
    w = I33["wall_per_shot_s"]
    assert near(w[0], t * 1e-5 / 4 + 1e-4, 1e-12) and near(w[1], t * 1e-5 + 1e-4, 1e-12)
    assert round(w[0], -1) == 530 and round(w[1], -2) == 2100          # '530--2100 s'
    assert near(I33["wall_per_shot_baseline_s"], 213, 0.001)           # 'rather than 213 s'
    s1, s2 = I33["shots_first"], I33["shots_campaign"]
    assert I33["wall_first_result_s"] == (s1[0] * w[0], s1[1] * w[1])
    assert I33["wall_campaign_s"] == (s2[0] * w[0], s2[1] * w[1]) == R33.wall_time_s
    assert round(I33["wall_first_result_days"][0]) == 86 and round(I33["wall_first_result_yr"][1], 1) == 2.7
    assert round(I33["wall_campaign_yr"][0], 1) == 1.6 and round(I33["wall_campaign_yr"][1]) == 18
    assert round(I33["horizon_ratio_campaign"][1], 1) == 3.7
    assert I33["fits_horizon_one_machine"] == (True, False)
    u = I33["wall_campaign_unary_iteration_yr"]
    assert round(u[0], 1) == 1.6 and round(u[1], 1) == 4.6 and u[1] < 5      # 'brings the campaign to 1.6--4.6 yr'
    b = I33["wall_campaign_baseline_yr"]
    assert round(b[0], 2) == 0.63 and round(b[1], 1) == 1.8                    # '0.63--1.8 yr'
    assert I33["subtrees_for_baseline"] == 3 and I33["subtree_f_star"] == 12
    assert 90 <= I33["extra_lq_for_baseline"] <= 110                           # 'about 100 of the idle LQ'
    # the floor and the factory counts come from common.depth_exports
    dx = c.depth_exports(t, I33["t_depth_per_shot"], s2, 1e-6, 1e-4)
    for k in ("f_star", "floor_wall_s", "factories_for_1yr"):
        assert I33[k] == dx[k]
    assert near(I33["floor_wall_s"][1], I33["wall_campaign_s"][1] - s2[1] * 1e-4, 1e-9)
    for s in ("Per-shot wall time & $530$--$2100\\,$s \\\\", "$86\\,$d--$2.7$ yr on one machine", "$1.6$--$18$ yr on one machine",
              "is $3.7\\times$ the 5-year horizon", "$1.6$--$4.6$ yr, which fits", "$0.63$--$1.8$ yr",
              "rather than $213\\,$s"):
        assert s in TEX, s
    # the model follows its inputs: F* = 10 everywhere restores the baseline wall
    r = m.model(replace(A, parallelism_2033=c.Cited((10, 10), "x")), "2033")
    assert near(r.intermediates["wall_campaign_s"][1], r.intermediates["wall_serial_baseline_campaign_s"][1], 1e-12)


def test_2033_banks_casher_demonstration():
    # R7 (simplobs.json C2): nu(Lambda)/V_F at a = 0.15 fm, Lambda = 0.15 (d = 430) and 0.10 (d = 640), 5 cfg, 30%
    bc = I33["bc"]
    assert [b["d"] for b in bc] == [430, 640]
    assert near(bc[0]["nu"], 27, 0.01) and near(bc[1]["nu"], 18, 0.01)       # '18--27 modes'
    assert round(bc[0]["lambda_mev"]) == 197
    t = I33["bc_t_per_shot"]
    assert t == (430 * 4 * 8192 * 13, 640 * 4 * 8192 * 13)
    assert near(t[0], 1.8e8, 0.02) and near(t[1], 2.7e8, 0.01)
    n = I33["bc_shots"]
    assert near(n[0], 1.1e4, 0.002) and near(n[1], 1.7e4, 0.015)              # simplobs: 1.1e4, 1.7e4
    assert round(I33["bc_wall_days"][0]) == 58 and round(I33["bc_wall_yr"][1], 1) == 1.5
    rows = m.INSTANCE_ROWS(A, "2033", R33)
    assert rows[1][2] == t and "Banks-Casher" in rows[1][0]
    assert "$1.8$--$2.7\\times 10^{8}$ T-gates per shot, $1.1$--$1.7\\times 10^{4}$ shots; $58\\,$d--$1.5$ yr" in TEX
    assert "there are 18--27 modes below $\\Lambda$" in TEX and "$\\sin(\\pi/16)=0.195$" in TEX
    assert abs(math.sin(math.pi / 16) - 0.195) < 5e-4


def test_e27_no_multiple_machines_in_tex_or_model():
    # Ruling E27 (H. Lamm, 2026-10-02): "we don't want to anywhere assume we have multiple machines."
    import re
    low = TEX.lower()
    for banned in ("device count", "duty cycle", "shot-parallel", "machine-year", "machine-week", "machines needed",
                   "across machines", "fold parallelism", "embarrassingly parallel", "factorizes exactly",
                   "serialized device time", "parallel machines"):
        assert banned not in low, banned
    # 'machines' survives only as the successive hardware generations ('10^3- and 10^4-LQ machines')
    assert all(low[max(0, i.start() - 3):i.start()] == "lq " for i in re.finditer(r"\bmachines\b", low))
    assert low.count("machines") <= 1
    assert "The campaign runs serially on one machine, and it is depth-bound." in TEX
    assert "Per-shot wall time & $74\\,$ms \\\\" in TEX                  # boxes carry final numbers only (E29)
    assert "$1.6$--$6.3\\times 10^{4}\\times$ the 5-year horizon" in TEX
    s = m.Assumptions()
    for name in ("n_machines", "parallel_machines", "machines_for_horizon", "n_devices"):
        assert not hasattr(s, name)
    for r in (R28, R33, RCD):
        assert not any("machines" in k or ("parallel" in k and "parallelism" not in k) for k in r.intermediates)
    assert s.shot_overhead_s.lo == 1e-4 and s.campaign_horizon_yr.lo == 5
    with pytest.raises(ValueError):
        replace(A, shot_overhead_s=c.Assumed(-1e-4, "x"))
    with pytest.raises(ValueError):
        replace(A, campaign_horizon_yr=c.Stated(0, "x"))


def test_2033_qme_alternative():
    # app12 route (c): one QET (2.1e8 T) is a fifth of the envelope; QME would put ~1e2 in one circuit (~2e10 T)
    assert near(I33["t_per_shot_qme_coherent"], 100 * R33.hard_ops[0], 1e-12)
    assert near(I33["t_per_shot_qme_coherent"], 2e10, 0.07)
    assert I33["t_per_shot_qme_coherent"] > 10 * A.t_cap_2033.lo
    assert "($2.1\\times 10^{8}$ T) is a fifth of the $10^9$ envelope" in TEX
    assert near(R33.hard_ops[0] / 1e9, 0.2, 0.07)


def test_2033_uncertainty_budget():
    assert "$30\\%$ (first result) and $20\\%$ (campaign) relative on $P$" in TEX
    assert I33["rel_first"] == 0.30 and I33["rel_campaign"] == 0.20
    assert I33["n_cfg_first"] == 1 and I33["n_cfg_campaign"] == 3


def test_2033_VF_from_chapter_formula():
    # app12:48: V_F = n_s N_c N_f V with n_s=1 (staggered); the box prints 'V_F = N_c V ~= 2.5e4 per flavor'.
    assert I33["V_F"] == 1 * 3 * 3 * 8192
    assert I33["V_F_single_flavor"] == 1 * 3 * 8192
    assert I33["V_F_box"] == 2.5e4
    assert near(I33["V_F_single_flavor"], I33["V_F_box"], 0.02)      # 24576 printed '~= 2.5e4'
    assert not near(I33["V_F"], I33["V_F_box"], 0.10)


# --- depth budgets and the 1e4-LQ rung -----------------------------------------------

def test_depth_budget_rule():
    assert near(I28["depth_budget_pf_over_eps_l"], 1e7, 1e-9)
    assert near(I33["depth_budget_pf_over_eps_l"], 1e9, 1e-9)
    assert near(ICD["depth_budget_pf_over_eps_l"], 1e11, 1e-9)


def test_codesign_24x48_instance():
    # R3 + R7: post-2033 rung, M at D V log2 V: '2.6e11 T ... on a 90-LQ algorithmic register (20 site + 4 column
    # index + 4 color + 12 QPE + ~50 ancilla), 2.6x that budget'; 'cut the per-shot T by at least 2.6x (the
    # unit-prefactor count, 2.57e11, against the 1e11 budget)'; 'LCU compression alone covers' it.
    assert 24 ** 3 * 48 == 663552 and math.ceil(math.log2(663552)) == 20
    assert [ICD[k] for k in ("site_qubits", "column_qubits", "color_qubits", "n_qpe", "n_anc")] == [20, 4, 4, 12, 50]
    assert RCD.lq == (90, 90) and ICD["lq_with_W_column_index"] == 92
    assert near(ICD["scaling_units_DVlogV"], 4 * 663552 * math.log2(663552), 1e-12)
    assert near(RCD.hard_ops[0], 5000 * 4 * 663552 * math.log2(663552), 1e-12)     # 2.5666e11
    assert near(RCD.hard_ops[0], 2.6e11, 0.014) and round(RCD.hard_ops[0] / 1e9) == 257
    assert near(ICD["t_per_shot_before_R3"], 1.0266e12, 1e-4)
    assert not ICD["fits_budget"]
    assert round(ICD["factor_over_budget"], 1) == 2.6
    assert ICD["per_shot_reduction"] == 2.6 and ICD["fits_budget_after_reduction"]
    assert ICD["factor_over_budget"] <= 2.6
    assert near(ICD["t_per_shot_after_reduction"], 9.9e10, 0.003)
    assert ICD["lever_range"] == (6, 20)
    assert ICD["fits_budget_at_lever_top"] and ICD["fits_budget_at_lever_bottom"]
    assert ICD["lcu_alone_closes"] and not ICD["even_odd_alone_closes"]
    assert ICD["d_inv"] == 5e3 and near(ICD["d_inv_lp_law"], 1e3 * 5, 1e-12)
    assert ICD["d_inv_ratio_vs_2033"] == 10
    assert near(ICD["d_inv_formula_ln"], 1.15e4, 0.01)
    assert near(ICD["d_inv_ratio_vs_2033_formula_ln"], 12.5, 1e-9)
    assert RCD.breakdown_total() == RCD.hard_ops[0]
    for s in ("costs $2.6\\times 10^{11}$ T", "on a 90-LQ algorithmic register (20 site + 4 column index",
              "at least $2.6{\\times}$ (the unit-prefactor count, $2.57\\times 10^{11}$",
              "LCU compression alone covers the $24^3{\\times}48$ instance", "past 2033"):
        assert s in TEX, s


def test_codesign_production_48x96_after_reduction():
    # '48^3x96 (5.0e12 T) needs ~50x, beyond the <~20x the levers below reach'; deliverables: 'approximately
    # 50-fold reduction'
    assert ICD["site_qubits_production"] == 24 and ICD["lq_production"] == 94
    assert near(ICD["t_per_shot_production"], 5000 * 4 * 48 ** 3 * 96 * math.log2(48 ** 3 * 96), 1e-12)
    assert near(ICD["t_per_shot_production"], 5.0e12, 0.01)
    assert not ICD["production_fits_after_reduction"]
    assert round(ICD["production_reduction_needed"]) == 50
    assert ICD["production_beyond_levers"]
    assert ICD["n_cfg_production"] == 1e3
    assert near(ICD["subroutine_calls_1e4_rung"], 1 / 1e-4 * 1e3, 1e-12)
    assert ICD["subroutine_calls_1e4_rung_box"] == 1e7
    assert "($5.0\\times 10^{12}$ T) needs $\\sim 50{\\times}$" in TEX and "approximately 50-fold reduction" in TEX


def test_codesign_one_machine_campaign_horizon():
    # 'a shot holds 9.9e10 T at a T-depth of 2.5--9.9e10 ... 2.9--11 days ... 7.8e4--3.1e5 years ... 1.6--6.3e4x
    # the 5-year horizon. Five years hold 160--640 shots.'
    t = ICD["t_per_shot_after_reduction"]
    assert ICD["f_star"] == (1, 4)
    d = ICD["t_depth_after_reduction"]
    assert near(d[0], 2.5e10, 0.02) and near(d[1], 9.9e10, 0.003)
    w = ICD["wall_per_shot_after_reduction_s"]
    assert near(w[0], t * 1e-5 / 4 + 1e-4, 1e-12) and near(w[1], t * 1e-5 + 1e-4, 1e-12)
    days = ICD["wall_per_shot_after_reduction_days"]
    assert round(days[0], 1) == 2.9 and round(days[1]) == 11
    y = ICD["wall_serial_after_reduction_yr"]
    assert near(y[0], 7.8e4, 0.003) and near(y[1], 3.1e5, 0.01)
    h = ICD["horizon_ratio_serial"]
    assert near(h[0], 1.6e4, 0.03) and near(h[1], 6.3e4, 0.01)
    n = ICD["shots_in_horizon"]
    assert round(n[0], -1) == 160 and round(n[1], -1) == 640
    assert not ICD["fits_horizon_one_machine"]
    assert ICD["qme_call_saving"] == 100
    assert RCD.wall_time_s is None and RCD.shots is None                      # conditional rung: not tabulated
    for s in ("$9.9\\times 10^{10}$ T at a T-depth of $2.5$--$9.9\\times 10^{10}$", "takes $2.9$--$11$ days",
              "$7.8\\times 10^{4}$--$3.1\\times 10^{5}$ years", "Five years hold 160--640 shots"):
        assert s in TEX, s


def test_synthesis_budget_1000lq_and_1e4lq():
    for I in (I33, ICD):
        assert I["synth_to_poly"] == 0.1
        assert near(I["poly_tolerance_normalized"], 0.5 * math.exp(-5), 1e-12)        # 3.369e-3
        assert near(I["eps_syn_2033"], 3.369e-4, 1e-3)
        assert near(I["eps_syn_2033_frac_of_peak"], 6.7e-4, 0.01)                      # 0.067%
        assert near(I["synth_over_poly_at_report_default"], 2.97, 0.01)                # '3x'
        assert near(I["synth_extra_t_per_rotation"], 2.81, 0.01)                       # '2.8 T'


def test_r9_exports_2028_and_2033():
    # CONTRACT rule 12: the seven exports, with f_star the cross-paired ratio
    for era, I in (("2028", I28), ("2033", I33)):
        for k in c.DEPTH_EXPORTS:
            assert k in I, (era, k)
        t, d = c._band(I["t_per_shot"]), c._band(I["t_depth_per_shot"])
        assert c._band(I["f_star"]) == pytest.approx((t[0] / d[1], t[1] / d[0]), rel=1e-12)
    assert I28["wall_first_result_s"] is None
    assert I33["wall_first_result_s"][0] < I33["wall_campaign_s"][0]


# --- bookkeeping ---------------------------------------------------------------------

def test_published_is_the_box_as_printed():
    assert m.PUBLISHED["2028"].lq == (18, 18) and m.PUBLISHED["2028"].hard_ops == (1.9e4, 3.5e4)   # SYN-BUDGET
    assert c.close(R28.lq, m.PUBLISHED["2028"].lq, m.PUBLISHED["2028"].rel_tol)
    assert c.close(R28.hard_ops, m.PUBLISHED["2028"].hard_ops, m.PUBLISHED["2028"].rel_tol)
    assert m.PUBLISHED["2033"].lq == (83, 83) and m.PUBLISHED["2033"].hard_ops == (2.1e8, 2.1e8)   # R3
    assert c.close(R33.hard_ops, m.PUBLISHED["2033"].hard_ops, m.PUBLISHED["2033"].rel_tol)
    assert m.PUBLISHED["codesign"].lq == (90, 90) and m.PUBLISHED["codesign"].hard_ops == (2.6e11, 2.6e11)
    assert c.close(RCD.hard_ops, m.PUBLISHED["codesign"].hard_ops, m.PUBLISHED["codesign"].rel_tol)


def test_instance_rows_carry_model_numbers():
    for era, r in (("2028", R28), ("2033", R33), ("codesign", RCD)):
        rows = m.INSTANCE_ROWS(A, era, r)
        assert len(rows) == (2 if era == "2033" else 1)
        lbl, lq, t, extra = rows[0]
        assert lq == r.lq and t == r.hard_ops and "eps_l" in extra
    lbl, _, _, extra = m.INSTANCE_ROWS(A, "2028", R28)[0]
    assert "V=4^4" in lbl and extra["volume"] == "4^4" and extra["lq_qme_route"] == 34
    assert extra["degree"] == 84 and near(extra["shots"], 3.1408e5, 1e-4)
    assert "success-probability" in m.INSTANCE_ROWS(A, "2033", R33)[0][0]
    assert "d_inv~5e3" in m.INSTANCE_ROWS(A, "codesign", RCD)[0][0]


def test_no_bare_literals_in_codesign_bookkeeping():
    r = m.model(replace(A, n_cfg_production=c.Stated(2e3, "x")), "codesign")
    assert near(r.intermediates["subroutine_calls_1e4_rung"], 2 * ICD["subroutine_calls_1e4_rung"], 1e-12)
    r = m.model(replace(A, n_anc_2033=c.Stated(60, "x")), "2033")
    assert r.lq == (93, 93) and r.intermediates["lq_box_sum"] == 93
    r = m.model(replace(A, n_cfg_campaign_2033=c.Stated(1, "x")), "2033")
    assert near(r.intermediates["wall_campaign_s"][1], R33.wall_time_s[1] / 3, 1e-12)


def test_assumptions_reject_nonsense():
    with pytest.raises(ValueError):
        replace(A, L_2028=c.Assumed(0, "x"))
    with pytest.raises(ValueError):
        replace(A, L_2028=c.Assumed(6, "x"))          # coordinate registers are binary
    with pytest.raises(ValueError):
        replace(A, d_log_2028=c.Assumed(111, "x"))    # the normalized log is even
    with pytest.raises(ValueError):
        replace(A, column_qubits_2028=c.Cited(3, "x"))   # 8 slots cannot index 9 entries
    with pytest.raises(ValueError):
        replace(A, target_rel_2028=c.Assumed(2.0, "x"))
    with pytest.raises(ValueError):
        replace(A, d_log_err=c.Assumed(2e-3, "x"))    # err(d) above err(d-2)
    with pytest.raises(ValueError):
        replace(A, log_base=c.Stated(3, "x"))
    with pytest.raises(ValueError):
        replace(A, lcu_saving=c.Stated((10, 3), "x"))
    with pytest.raises(ValueError):
        replace(A, toffoli_convention=c.Cited("whatever", "x"))
    with pytest.raises(ValueError):
        replace(A, parallelism_2033=c.Cited((0.5, 4), "x"))
    with pytest.raises(ValueError):
        replace(A, L_s_2033=c.Stated(6, "x"))         # 6^3 x 16 is not V_2033
    with pytest.raises(ValueError):
        replace(A, rel_first_2033=c.Stated(1.5, "x"))
    with pytest.raises(ValueError):
        m.model(A, "2040")


def test_r16_rms_shot_convention_in_tex():
    # R16: no confidence-1-delta language survives in the chapter; Eq. Nshot is sigma^2/eps^2.
    assert "\\delta" not in TEX and "confidence" not in TEX
    assert "\\sigma^{2}/\\varepsilon^{2}\\cdot N_\\mathrm{cfg}" in TEX
    assert "$\\sim 3.1\\!\\times\\!10^5$" in TEX and "$\\gtrsim 10^7$" in TEX
    for old in ("4.6\\times 10^{4}", "4.6\\times 10^7", "4\\times 10^{15}", "1.4\\!\\times\\!10^6",
                "\\sim 9\\times 10^{14}", "85 LQ", "92-LQ", "sub-percent precision on the normalized trace"):
        assert old not in TEX, old
    s = m.Assumptions()
    assert not hasattr(s, "delta_2028") and not hasattr(s, "delta_2033")


def test_alpha_be_is_a_fiducial_subnormalization():
    # author (2026-10-02): "say we take fiducial"; Eq. B11 of arxiv_2407_13080 is the spectral bound
    assert ("We take $\\alpha_{\\rm BE} = 16 K^2 + m_0^2$ as a fiducial subnormalization, equal to the spectral "
            "bound on $\\|W\\|$ (Eq.~B11 of~\\cite{arxiv_2407_13080})") in TEX
    assert "the paper's free-circuit block encoding has $s \\ge 2(m_0^2 + 2K^2)$" in TEX
    assert near(I28["be_subnormalization_s"], 2 * (0.16 + 2), 1e-12)


def test_referee_h1_h4_scope_and_drafting_history():
    # Referee H1 (2026-10-04): log det(M^dag M) = 2 log|det M| has no phase; the primitive does not lift the
    # finite-density volume cap, and the utility box no longer counts finite-density studies.
    assert "$\\log\\det(M^\\dagger M)=2\\log|\\det M|$ carries no phase" in TEX
    assert "The primitive priced here does not lift that cap." in TEX
    assert "lifting the $V{\\sim}16^3$ classical volume cap" not in TEX and "reweighting at $\\mu_B" not in TEX
    # H3: the device is not inside the HMC force or accept/reject step
    assert "stays classical and unchanged" not in TEX
    assert "The device does not enter the HMC force or accept/reject step" in TEX
    assert "it supplies no HMC forces or accept/reject decisions" in TEX
    # H4: the post-2033 stage is defined by eps_l, not width
    assert "$10^4$-LQ" not in TEX and "10000 LQ" not in TEX
    assert "defined by its error rate rather than its width" in TEX
    # editorial: no drafting history, no 'rung', the G6 heading
    low = TEX.lower()
    for word in ("retired", "earlier plan", "rung", "lever (ii), a tighter"):
        assert word not in low, word
    assert "\\subsection*{Classical baseline and prospective quantum advantage}" in TEX


def test_referee_h2_svt_readout_identity_and_hvp_estimate():
    # Referee H2: Re Tr M^-1 = am Tr (M^dag M)^-1 needs M = am + A with A anti-Hermitian (staggered, mu_B = 0);
    # the odd SVT of M^dag / alpha with p(y) = p0 y0 / y is p0 am M^-1, and its flag probability averages to
    # p0^2 am (1/V_F) Re Tr M^-1. Checked on random U(1) links at V = 4^4 (D = 4) and am = 0.04.
    import itertools
    import numpy as np
    rng = np.random.default_rng(7)
    L, Dd, am, p0 = 4, 4, 0.04, 0.5
    sites = list(itertools.product(range(L), repeat=Dd))
    idx = {x: i for i, x in enumerate(sites)}
    V = len(sites)
    A_ = np.zeros((V, V), complex)
    for x in sites:
        for mu in range(Dd):
            eta = (-1) ** sum(x[:mu])
            y = list(x); y[mu] = (y[mu] + 1) % L; y = tuple(y)
            u = np.exp(1j * rng.uniform(0, 2 * np.pi))
            A_[idx[x], idx[y]] += 0.5 * eta * u
            A_[idx[y], idx[x]] -= 0.5 * eta * np.conj(u)
    M = am * np.eye(V) + A_
    Minv = np.linalg.inv(M)
    c = np.trace(Minv).real / V
    assert near(c, am * np.trace(np.linalg.inv(M.conj().T @ M)).real / V, 1e-10)
    alpha = 4 + am
    U_, S_, Vh = np.linalg.svd(M.conj().T / alpha)
    assert S_.min() >= am / alpha - 1e-12
    svt = U_ @ np.diag(p0 * (am / alpha) / S_) @ Vh
    assert np.allclose(svt, p0 * am * Minv)
    P = np.mean(np.sum(np.abs(svt) ** 2, axis=0))
    assert near(P, p0 ** 2 * am * c, 1e-10)
    assert "normality alone is not enough" in TEX and "nothing is postselected" in TEX
    # HVP: noise-matched loops need V_F var / (p0 am)^2 per configuration, ~8e12 at 48^3x96, kappa = 1e3
    assert near(ICD["am_production"], 4 / 999, 1e-12) and ICD["v_f_production"] == 3 * 48 ** 3 * 96
    assert near(ICD["hvp_noise_matched_shots_per_cfg"], 7.95e12, 0.001)
    assert ICD["hvp_noise_matched_over_box"] > 1e8
    assert "$1.5\\times 10^{12}$ at $24^3{\\times}48$ and $2.4\\times 10^{13}$ at $48^3{\\times}96$" in TEX
    assert "a smaller site variance raises it proportionally" in TEX
    # verifier (2026-10-04): N = V_F / (var (p0 am)^2) scales inversely with the site variance
    assert ICD["v_f_codesign"] == 3 * 24 ** 3 * 48
    assert near(ICD["hvp_noise_matched_shots_per_cfg_24c48"], 4.967e11, 0.001)
    assert 4e7 < ICD["hvp_noise_matched_over_box_24c48"] < 6e7
    assert "noise-matched HVP loops need $1.5\\times 10^{12}$ shots per cfg at $24^3{\\times}48$" in TEX
    assert ICD["hvp_current_directions"] == 3
    half = m.model(replace(A, hvp_site_variance=m.Assumed(0.5, "x")), "codesign").intermediates
    assert near(half["hvp_noise_matched_shots_per_cfg"], 2 * ICD["hvp_noise_matched_shots_per_cfg"], 1e-12)
    assert near(half["hvp_noise_matched_shots_per_cfg_24c48"], 2 * ICD["hvp_noise_matched_shots_per_cfg_24c48"], 1e-12)
    assert "A $0.3\\%$ target is not adopted (see below)." in TEX


def test_referee_h2_hvp_mapping_and_error_budget():
    # Referee H2 (2026-10-05): loops L_k(t) = Tr(Gamma_{k,t} M^-1) are purely imaginary (staggered eps symmetry);
    # Hadamard-shot variance on L_k(t) >= 3 V_s / p0^2 x a slice-diluted classical noise vector; QME vs classical
    # count; WP25 (arXiv:2505.21476) disconnected errors bound the trace-only cut at 1.3-1.4x.
    import itertools
    import math as _m
    import numpy as np
    rng = np.random.default_rng(7)
    L, Dd, am, p0 = 4, 4, 0.04, 0.5
    sites = list(itertools.product(range(L), repeat=Dd))
    idx = {x: i for i, x in enumerate(sites)}
    V = len(sites)
    A_ = np.zeros((V, V), complex)
    G = np.zeros((V, V), complex)
    for x in sites:
        for mu in range(Dd):
            eta = (-1) ** sum(x[:mu])
            y = list(x); y[mu] = (y[mu] + 1) % L; y = tuple(y)
            u = np.exp(1j * rng.uniform(0, 2 * np.pi))
            A_[idx[x], idx[y]] += 0.5 * eta * u
            A_[idx[y], idx[x]] -= 0.5 * eta * np.conj(u)
            if mu == 1:
                G[idx[x], idx[y]] += 0.5 * eta * u
                G[idx[y], idx[x]] += 0.5 * eta * np.conj(u)
    assert np.linalg.norm(G, 2) <= 1 + 1e-12
    B = G @ np.linalg.inv(am * np.eye(V) + A_)
    for t in range(L):
        sl = [idx[x] for x in sites if x[3] == t]
        tr = np.trace(B[np.ix_(sl, sl)])
        assert abs(tr.real) < 1e-10 and abs(tr.imag) > 0.1
        hutch = sum(abs(B[i, j]) ** 2 for i in sl for j in sl if i != j)
        assert hutch <= len(sl) / am ** 2
        assert (len(sl) / (p0 * am)) ** 2 / hutch >= len(sl) / p0 ** 2
    assert near(ICD["hvp_noise_matched_shots_per_cfg_all_dirs_24c48"], 1.490e12, 0.001)
    assert near(ICD["hvp_noise_matched_shots_per_cfg_all_dirs"], 2.384e13, 0.001)
    assert near(ICD["hvp_insertion_share_of_shot"], 5e-5, 1e-12)
    assert ICD["hvp_shot_vs_noise_vector_variance_ratio_lb_24c48"] == 4 * 3 * 24 ** 3
    assert near(ICD["hvp_qme_calls_per_cfg_24c48"], 1.465e7, 0.001)
    assert near(ICD["hvp_classical_vectors_per_cfg_ub_24c48"], 8.982e6, 0.001)
    assert 1 < ICD["hvp_qme_calls_per_cfg_24c48"] / ICD["hvp_classical_vectors_per_cfg_ub_24c48"] < 2   # at var = 1
    assert 0.9e16 < ICD["hvp_qme_coherent_t_per_circuit_24c48"] < 1.1e16
    assert 0.9e5 < ICD["hvp_qme_coherent_over_budget_24c48"] < 1.1e5   # "$10^{5}\\times$ the $10^{11}$-T budget"
    lo, hi = ICD["hvp_noise_matched_wall_per_cfg_yr_24c48"]
    assert 1.1e10 < lo < 1.2e10 and 4.6e10 < hi < 4.7e10
    # WP25 tab:amu_fbf (1e-10): BMW-20 -15.4(1.2)(1.4), Mainz/CLS-24 -16.1(1.2)(1.2); total LO HVP 713.2(6.1)
    for disc, st, sy in ((-15.4, 1.2, 1.4), (-16.1, 1.2, 1.2)):
        assert 1.3 <= _m.hypot(st, sy) / sy <= 1.42
        assert 24.5 <= sy / (0.003 * abs(disc)) <= 30.5    # 24.8 (Mainz), 30.3 (BMW): "25--30"
        assert 0.003 * 713.2 > _m.hypot(st, sy)
    for s_ in ("1.5\\times 10^{7}$ calls", "$\\le 9\\times 10^{6}$ classical vectors", "$1.7\\times 10^{5}$ at $24^3$",
               "=5\\times 10^{-5}$ of the per-shot T", "by at most $1.3$--$1.4{\\times}$", "$25$--$30{\\times}$ below",
               "$10^{5}\\times$ the $10^{11}$-T budget", "which the disconnected errors already meet",
               "at most a $1.3$--$1.4{\\times}$ cut", "($1$--$5\\times 10^{10}$ yr each)",
               "purely imaginary", "no $0.3\\%$ target"):
        assert s_ in TEX, s_
    assert "uncosted HVP estimator" not in TEX and "not yet costed" not in TEX


def test_utility_box_one_time_basis():
    # G5 open item 4 (Claude's decision, 2026-10-04): value and capital base on one time basis (one-time);
    # $5M/yr over the 5-year horizon = $25M, 25/146 = 17% of the $46M + ~$100M capital base; the unsourced
    # ~$5M/yr DUNE/Mu2e add-on carries no dollar figure.
    box = TEX.split("begin{utilitybox}")[1].split("end{utilitybox}")[0]
    assert "$\\sim\\$25$M, one-time, at the post-2033 stage" in box
    assert round(100 * 5 * 5 / (46 + 100)) == 17 and "about $17\\%$ of that one-time $\\sim\\$150$M capital base" in box
    # smaller ruling (c), 2026-10-05: the ~$25M premise is kept; its basis is the DOE-reported capital base
    # ($46M project, DOE_PMAA_g2_2019; $100M reused equipment, DOE_g2_breakthrough_2026), stated in the box
    assert "both as DOE reports them" in box and "DOE_PMAA_g2_2019" in box and "DOE_g2_breakthrough_2026" in box
    assert "the $17\\%$ share is a stated, rescalable assumption" in box
    for gone in ("\\$10$M/yr", "Another $\\sim\\$5$M/yr", "\\$10M/yr"):
        assert gone not in box, gone
