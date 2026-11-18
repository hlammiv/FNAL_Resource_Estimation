"""Correctness checks for the Ch. 4 2028 benchmark circuit.

    OMP_NUM_THREADS=4 python -m pytest circuits/ch04_logdet_staggered/test_circuit.py -q

Every reference below is built with numpy from the operator definitions, not from build.py:
the staggered matrix M of arXiv:2407.13080 Eq. 16 on a periodic lattice, W = M^T M, the index map
c(n, l) of Eqs. 4-8 with add -> add2, the entries of Eqs. 24-26, and f(W/s) by eigendecomposition
of the dense W. The circuits are simulated as statevectors (qiskit-aer, 4 threads, 22 qubits at most).

Tolerances. TOL = 1e-9 on every amplitude compared (the statevector runs are double precision, the
stored QSP phases reproduce the polynomial to 1e-14, and the worst deviation seen is 1e-13).
The polynomial itself approximates the shifted normalized logarithm to the model's tolerance
1.784e-3, and that is checked separately against the exact spectrum.
"""
from __future__ import annotations

import itertools
import json
import math
import sys
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
for p in (str(REPO), str(HERE)):          # HERE goes in last, so the modules of this directory are found first
    if p in sys.path:
        sys.path.remove(p)
    sys.path.insert(0, p)

import build                      # noqa: E402
import count                      # noqa: E402
import qsp                        # noqa: E402
from sim import simulate          # noqa: E402
from estimates.common import mcx_toffoli_count          # noqa: E402

TOL = 1e-9
RNG = np.random.default_rng(20261118)


# --------------------------------------------------------------------------- #
# Register layout
# --------------------------------------------------------------------------- #

def layout(inst) -> dict:
    """Bit offsets of the registers in the statevector index (qubit 0 is the lowest bit)."""
    off, out = 0, {}
    for name, size in (("site", inst.n_site), ("col", inst.n_col), ("flag", 1), ("sig", 1),
                       ("work", inst.n_work), ("had", 1)):
        out[name] = (off, size)
        off += size
    assert off == inst.n_qubits
    return out


# --------------------------------------------------------------------------- #
# Independent references (numpy, from the definitions)
# --------------------------------------------------------------------------- #

def site_index(x, L) -> int:
    return sum(xi * L ** mu for mu, xi in enumerate(x))


def ref_M(L, D, m0, K) -> np.ndarray:
    """Free staggered matrix, arXiv:2407.13080 Eq. 16, periodic:
    M_mn = (K/2) sum_mu eta_mu(m) (delta_{n,m+mu} - delta_{n,m-mu}) + m0 delta_mn, eta_mu(m) = (-1)^{sum_{nu<mu} m_nu}."""
    V = L ** D
    M = np.zeros((V, V))
    for x in itertools.product(range(L), repeat=D):
        m = site_index(x, L)
        M[m, m] += m0
        for mu in range(D):
            eta = (-1) ** sum(x[:mu])
            up, dn = list(x), list(x)
            up[mu] = (up[mu] + 1) % L
            dn[mu] = (dn[mu] - 1) % L
            M[m, site_index(up, L)] += 0.5 * K * eta
            M[m, site_index(dn, L)] -= 0.5 * K * eta
    return M


def ref_W(L, D, m0, K) -> np.ndarray:
    M = ref_M(L, D, m0, K)
    return M.T @ M


def ref_shift_permutation(L, D) -> np.ndarray:
    """c(n, l) of arXiv:2407.13080 Eqs. 4-8 with add -> add2, sub -> sub2 (Sec. II.2), for a column
    register of 4D values: l < D is +2 mu, D <= l < 2D is -2 mu, l >= 2D is the diagonal.
    Returns perm with perm[l, n] = c(n, l)."""
    V = L ** D
    perm = np.empty((4 * D, V), dtype=int)
    for x in itertools.product(range(L), repeat=D):
        n = site_index(x, L)
        for l in range(4 * D):
            y = list(x)
            if l < D:
                y[l] = (y[l] + 2) % L
            elif l < 2 * D:
                y[l - D] = (y[l - D] - 2) % L
            perm[l, n] = site_index(y, L)
    return perm


def ref_f_of_W(W, s, cheb) -> np.ndarray:
    lam, U = np.linalg.eigh(W / s)
    return (U * np.polynomial.chebyshev.chebval(lam, cheb)) @ U.T


def random_state(dim) -> np.ndarray:
    v = RNG.normal(size=dim) + 1j * RNG.normal(size=dim)
    return v / np.linalg.norm(v)


INSTANCES = [
    dict(L=4, D=4, m0=0.4, K=1.0, variant="literal"),      # the benchmark instance, 18 qubits
    dict(L=4, D=4, m0=0.4, K=1.0, variant="merged"),
    dict(L=4, D=4, m0=0.4, K=1.0, variant="compact"),
    dict(L=4, D=2, m0=0.3, K=0.8, variant="compact"),
    dict(L=4, D=4, m0=0.7, K=0.9, s=6.5, variant="literal"),   # generic s: theta_0 != 0, 4 O_A rotations
    dict(L=4, D=4, m0=0.7, K=0.9, s=6.5, variant="merged"),
    dict(L=8, D=2, m0=0.4, K=1.0, variant="literal"),      # add2 != sub2, two-bit incrementer
    dict(L=16, D=1, m0=0.5, K=1.2, variant="literal"),     # three-bit incrementer
    dict(L=8, D=4, m0=0.4, K=1.0, variant="literal"),      # 22 qubits
]
IDS = [f"L{i['L']}D{i['D']}{'s' if 's' in i else ''}-{i['variant']}" for i in INSTANCES]


# --------------------------------------------------------------------------- #
# The reference itself
# --------------------------------------------------------------------------- #

def test_reference_W_is_eq24_and_matches_the_model():
    """M^T M built from Eq. 16 equals Eq. 24, and its log det is the model's 150.5528."""
    L, D, m0, K = 4, 4, 0.4, 1.0
    W = ref_W(L, D, m0, K)
    E = np.zeros_like(W)
    perm = ref_shift_permutation(L, D)
    for n in range(L ** D):
        E[n, n] += m0 ** 2 + 2 * K ** 2
        for l in range(2 * D):
            E[perm[l, n], n] -= K ** 2 / 4
    assert np.abs(W - E).max() < 1e-13
    mdl = count.model_items()
    assert abs(np.linalg.slogdet(W)[1] - mdl["logdet_W"]) < 1e-9
    assert abs(np.linalg.eigvalsh(W)[0] - 0.16) < 1e-12 and abs(np.linalg.eigvalsh(W)[-1] - 4.16) < 1e-12


# --------------------------------------------------------------------------- #
# Primitives
# --------------------------------------------------------------------------- #

def test_ry_is_written_correctly_with_rz():
    from qiskit import QuantumCircuit
    from qiskit.quantum_info import Operator
    for th in (0.3, -1.7, 2.75427605270114):
        qc = QuantumCircuit(1)
        build.ry(qc, th, 0)
        ref = np.array([[math.cos(th / 2), -math.sin(th / 2)], [math.sin(th / 2), math.cos(th / 2)]])
        assert np.abs(Operator(qc).data - ref).max() < 1e-14
        assert [g.operation.name for g in qc.data] == ["sdg", "h", "rz", "h", "s"]


@pytest.mark.parametrize("n", [1, 2, 3, 4, 5, 6])
def test_mcx_ladder_is_a_multicontrolled_x_with_clean_ancillas(n):
    from qiskit import QuantumCircuit
    from qiskit.quantum_info import Operator
    n_work = max(n - 2, 0)
    qc = QuantumCircuit(n + 1 + n_work)
    build.mcx_ladder(qc, list(range(n)), n, list(range(n + 1, n + 1 + n_work)))
    U = Operator(qc).data
    dim = 2 ** (n + 1)
    ref = np.eye(dim)
    a, b = 2 ** n - 1, 2 ** (n + 1) - 1                      # all controls 1, target 0 / 1
    ref[[a, b]] = ref[[b, a]]
    assert np.abs(U[:dim, :dim] - ref).max() < 1e-14         # ancillas start and end in 0
    assert qc.count_ops().get("ccx", 0) == mcx_toffoli_count(n, "linear")   # the model's 2n - 3


@pytest.mark.parametrize("n", [2, 3, 4, 5, 6])
def test_and_tree_writes_the_and(n):
    from qiskit import QuantumCircuit
    from qiskit.quantum_info import Operator
    qc = QuantumCircuit(2 * n - 1)                           # n controls, target, n - 2 ancillas
    steps = build.and_tree_steps(list(range(n)), n, list(range(n + 1, 2 * n - 1)))
    for a, b, t in steps:
        qc.ccx(a, b, t)
    assert len(steps) == n - 1
    U = Operator(qc).data
    for c in range(2 ** n):                                  # every control pattern, target and ancillas 0
        out = int(np.argmax(np.abs(U[:, c])))
        assert abs(abs(U[out, c]) - 1) < 1e-14
        assert (out >> n) & 1 == int(c == 2 ** n - 1) and out & (2 ** n - 1) == c
    for a, b, t in reversed(steps):
        qc.ccx(a, b, t)
    assert np.abs(Operator(qc).data - np.eye(2 ** (2 * n - 1))).max() < 1e-14


# --------------------------------------------------------------------------- #
# O_c, O_A, the block encoding
# --------------------------------------------------------------------------- #

def _shift_error(inst, qc_oc) -> float:
    """Largest amplitude error of O_c against the reference permutation, on a random state over
    (site, col) with every other qubit in 0. The X on the top column bit undoes build.py's frame."""
    from qiskit import QuantumCircuit
    lay = layout(inst)
    top = lay["col"][0] + inst.n_col - 1
    full = QuantumCircuit(*qc_oc.qregs)
    full.x(top)
    full.compose(qc_oc, inplace=True)
    full.x(top)
    nsc = inst.n_site + inst.n_col
    amp = random_state(2 ** nsc)
    init = np.zeros(2 ** inst.n_qubits, dtype=complex)
    init[:2 ** nsc] = amp
    out = simulate(full, init)
    perm = ref_shift_permutation(inst.L, inst.D)
    ref = np.zeros_like(init)
    V = inst.V
    for l in range(4 * inst.D):
        ref[(l << inst.n_site) + perm[l]] = amp[(l << inst.n_site) + np.arange(V)]
    return float(np.abs(out - ref).max())


@pytest.mark.parametrize("kw", INSTANCES, ids=IDS)
def test_shift_oracle_is_the_index_map(kw):
    inst = build.Instance(**kw)
    assert _shift_error(inst, build.shift_oracle_circuit(inst)) < TOL


@pytest.mark.parametrize("kw", INSTANCES, ids=IDS)
def test_entry_oracle_amplitudes(kw):
    """For column value l the flag amplitude on |0> is 2^c A_l / (slots of that entry) / s:
    hop slots carry -K^2/4 each, the 2^(c-1) diagonal slots share m0^2 + D K^2 / 2."""
    from qiskit import QuantumCircuit
    inst = build.Instance(**kw)
    lay = layout(inst)
    qc_oa = build.entry_oracle_circuit(inst)
    top = lay["col"][0] + inst.n_col - 1
    full = QuantumCircuit(*qc_oa.qregs)
    full.x(top)
    full.compose(qc_oa, inplace=True)
    full.x(top)
    ncol = 2 ** inst.n_col
    amp = random_state(ncol)
    init = np.zeros(2 ** inst.n_qubits, dtype=complex)
    init[np.arange(ncol) << inst.n_site] = amp
    out = simulate(full, init)
    s = inst.sub
    diag, hop = inst.m0 ** 2 + inst.D * inst.K ** 2 / 2, -inst.K ** 2 / 4
    for l in range(ncol):
        a0 = ncol * hop / s if l < 2 * inst.D else 2 * diag / s
        a1 = math.sqrt(max(0.0, 1 - a0 * a0))
        base = l << inst.n_site
        assert abs(out[base] - amp[l] * a0) < TOL
        assert abs(abs(out[base + (1 << lay["flag"][0])]) - abs(amp[l]) * a1) < TOL


def _block_error(inst, qc_u, W) -> float:
    """Largest error of <0 0^c| U |0 0^c> against W/s, on basis and random site states."""
    V = inst.V
    vecs = [np.eye(V)[n] for n in sorted({0, 1, V // 3, V - 1})] + [random_state(V) for _ in range(3)]
    worst = 0.0
    for v in vecs:
        init = np.zeros(2 ** inst.n_qubits, dtype=complex)
        init[:V] = v
        out = simulate(qc_u, init)
        worst = max(worst, float(np.abs(out[:V] - (W / inst.sub) @ v).max()))
    return worst


@pytest.mark.parametrize("kw", INSTANCES, ids=IDS)
def test_block_encoding_is_W_over_s(kw):
    inst = build.Instance(**kw)
    W = ref_W(inst.L, inst.D, inst.m0, inst.K)
    assert _block_error(inst, build.block_encoding_circuit(inst), W) < TOL
    assert _block_error(inst, build.block_encoding_circuit(inst, inverse=True), W) < TOL


@pytest.mark.parametrize("kw", INSTANCES[:8], ids=IDS[:8])
def test_block_encoding_inverse_undoes_it_on_the_whole_register(kw):
    """U^dag U = 1 on a random state of every qubit, workspace included."""
    inst = build.Instance(**kw)
    qc = build.block_encoding_circuit(inst)
    qc.compose(build.block_encoding_circuit(inst, inverse=True), inplace=True)
    init = random_state(2 ** inst.n_qubits)
    assert np.abs(simulate(qc, init) - init).max() < TOL


def test_block_encoding_check_catches_a_wrong_sign():
    """Flip the sign of the hopping entry (theta_1 -> 2 pi - theta_1): the block check must fail."""
    inst = build.Instance(**build.BENCHMARK, variant="merged")
    W = ref_W(inst.L, inst.D, inst.m0, inst.K)
    bad = build.block_encoding_circuit(inst, theta1_override=2 * math.pi - inst.thetas[1])
    assert _block_error(inst, bad, W) > 1e-2
    assert _block_error(inst, build.block_encoding_circuit(inst), W) < TOL


@pytest.mark.parametrize("variant", ["literal", "merged", "compact"])
def test_block_encoding_check_catches_a_dropped_term(variant):
    """Remove the last Toffoli that targets the site register (one hopping term): both the O_c check
    and the block check must fail."""
    from qiskit import QuantumCircuit
    inst = build.Instance(**build.BENCHMARK, variant=variant)
    lay = layout(inst)

    def drop(qc):
        idx = [i for i, g in enumerate(qc.data)
               if g.operation.name == "ccx" and qc.find_bit(g.qubits[2]).index < lay["site"][1]]
        out = QuantumCircuit(*qc.qregs)
        for i, g in enumerate(qc.data):
            if i != idx[-1]:
                out.append(g.operation, g.qubits)
        return out

    W = ref_W(inst.L, inst.D, inst.m0, inst.K)
    assert _shift_error(inst, drop(build.shift_oracle_circuit(inst))) > 1e-3
    assert _block_error(inst, drop(build.block_encoding_circuit(inst)), W) > 1e-3


# --------------------------------------------------------------------------- #
# Projector-controlled phase
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("kw", [INSTANCES[0], INSTANCES[3], INSTANCES[7]], ids=[IDS[0], IDS[3], IDS[7]])
@pytest.mark.parametrize("controlled", [True, False])
def test_projector_phase(kw, controlled):
    """exp(i phi (2 Pi - 1)), Pi = |0><0| on column and flag, applied when had = 1 (controlled) or always.
    Random state over site, col, flag and had, with sig and work in 0, which must come back in 0."""
    inst = build.Instance(**kw)
    lay = layout(inst)
    phi = 0.7361
    idx = np.arange(2 ** inst.n_qubits)
    anc_mask = (1 << lay["sig"][0]) | (((1 << inst.n_work) - 1) << lay["work"][0])
    live = (idx & anc_mask) == 0
    init = np.zeros(2 ** inst.n_qubits, dtype=complex)
    init[live] = random_state(int(live.sum()))
    out = simulate(build.projector_phase_circuit(inst, phi, controlled), init)
    be_mask = (((1 << inst.n_col) - 1) << lay["col"][0]) | (1 << lay["flag"][0])
    in_pi = (idx & be_mask) == 0
    phase = np.where(in_pi, np.exp(1j * phi), np.exp(-1j * phi))
    if controlled:
        phase = np.where((idx >> lay["had"][0]) & 1 == 1, phase, 1.0)
    assert np.abs(out - phase * init).max() < TOL


# --------------------------------------------------------------------------- #
# QSP conventions and the benchmark polynomial
# --------------------------------------------------------------------------- #

def test_qsp_conventions_and_solver():
    cheb = np.zeros(11)
    cheb[0::2] = [0.1, -0.3, 0.2, 0.05, -0.15, 0.1]
    wx, resid = qsp.solve_symmetric_phases(cheb)
    x = np.linspace(-1, 1, 501)
    f = np.polynomial.chebyshev.chebval(x, cheb)
    assert resid < 1e-12 and np.abs(wx - wx[::-1]).max() == 0
    assert np.abs(qsp.wx_response(wx, x).real - f).max() < 1e-12
    refl = qsp.wx_to_reflection(wx)
    assert len(refl) == len(wx) - 1
    assert np.abs(qsp.reflection_response(refl, x) - qsp.wx_response(wx, x)).max() < 1e-12
    # the reflection response against explicit 2x2 matrices at one point
    x0 = 0.37
    R = np.array([[x0, math.sqrt(1 - x0 * x0)], [math.sqrt(1 - x0 * x0), -x0]])
    U = np.eye(2, dtype=complex)
    for phi in refl:
        U = U @ np.diag([np.exp(1j * phi), np.exp(-1j * phi)]) @ R
    assert abs(U[0, 0] - qsp.reflection_response(refl, x0)[0]) < 1e-13


def test_benchmark_polynomial_meets_the_models_tolerance():
    """Stored coefficients: degree 84, even, within the model's tolerance of the shifted normalized log
    on [lambda_min/s, 1], bounded by 1, reproduced by the stored phases."""
    data = build.benchmark_polynomial()
    mdl = count.model_items()
    cheb = np.array(data["cheb_coeffs"])
    assert data["d"] == mdl["d"] == len(cheb) - 1 == len(data["phases_reflection"])
    assert np.abs(cheb[1::2]).max() == 0
    lam = 0.16 / 4.32
    x = np.linspace(lam, 1, 200001)
    err = np.abs(np.polynomial.chebyshev.chebval(x, cheb) - (np.log(x) / math.log(1 / lam) + data["shift"])).max()
    assert err < mdl["poly_tolerance"]
    assert abs(err - data["lp_error"]) < 2e-5
    y = np.linspace(-1, 1, 200001)
    assert np.abs(np.polynomial.chebyshev.chebval(y, cheb)).max() < 1
    resp = qsp.reflection_response(np.array(data["phases_reflection"]), y)
    assert np.abs(resp.real - np.polynomial.chebyshev.chebval(y, cheb)).max() < 1e-12


# --------------------------------------------------------------------------- #
# QET and the Hadamard test
# --------------------------------------------------------------------------- #

def _random_even_poly(d):
    cheb = np.zeros(d + 1)
    cheb[0::2] = RNG.normal(size=d // 2 + 1)
    x = np.linspace(-1, 1, 4001)
    cheb *= 0.7 / np.abs(np.polynomial.chebyshev.chebval(x, cheb)).max()
    return cheb


def _hadamard_branches(inst, qc, site):
    """Run the Hadamard-test circuit built without its final H. Returns (had = 0 branch, had = 1 branch),
    each normalized to one, with the phase convention fixed so the had = 0 amplitude on the input is real."""
    out = simulate(qc)
    half = 2 ** (inst.n_qubits - 1)
    b0, b1 = out[:half] * math.sqrt(2), out[half:] * math.sqrt(2)
    g = b0[site] / abs(b0[site])
    return b0 / g, b1 / g


@pytest.mark.parametrize("kw", INSTANCES[3:8], ids=IDS[3:8])
@pytest.mark.parametrize("d", [2, 8])
def test_qet_on_reduced_instances(kw, d):
    """Uncontrolled QET and the Hadamard test with a random even polynomial of degree d:
    Re <0|U_Phi|0> = f(W/s) column by column."""
    inst = build.Instance(**kw)
    W = ref_W(inst.L, inst.D, inst.m0, inst.K)
    cheb = _random_even_poly(d)
    refl = qsp.wx_to_reflection(qsp.solve_symmetric_phases(cheb)[0])
    fW = ref_f_of_W(W, inst.sub, cheb)
    V = inst.V
    site = V // 3
    init = np.zeros(2 ** inst.n_qubits, dtype=complex)
    init[site] = 1.0
    out = simulate(build.qet_circuit(inst, refl), init)
    assert np.abs(out[:V].real - fW[:, site]).max() < TOL
    for fold in (False, True):
        qc = build.hadamard_test_circuit(inst, refl, site=site, fold_last_phase=fold, final_h=False)
        b0, b1 = _hadamard_branches(inst, qc, site)
        e = np.zeros_like(b0)
        e[site] = 1.0
        assert np.abs(b0 - e).max() < TOL                    # phases off: the queries cancel exactly
        assert np.abs(b1[:V].real - fW[:, site]).max() < TOL
        assert abs(b1[site].real - fW[site, site]) < TOL     # what the Hadamard test measures
        full = build.hadamard_test_circuit(inst, refl, site=site, fold_last_phase=fold)
        psi = simulate(full)
        p1 = float(np.sum(np.abs(psi[2 ** (inst.n_qubits - 1):]) ** 2))
        assert abs((1 - 2 * p1) - fW[site, site]) < TOL


@pytest.fixture(scope="module")
def benchmark_reference():
    inst = build.Instance(**build.BENCHMARK)
    W = ref_W(inst.L, inst.D, inst.m0, inst.K)
    data = build.benchmark_polynomial()
    return W, data, ref_f_of_W(W, inst.sub, data["cheb_coeffs"])


@pytest.mark.parametrize("variant,site", [("literal", 0), ("literal", 173), ("merged", 37), ("compact", 222)])
def test_full_benchmark_circuit_against_the_exact_matrix(benchmark_reference, variant, site):
    """The full 18-qubit, d = 84 circuit. On the had = 1 branch the block-encoding ancillas in 0 hold
    P(W/s)|site>, and its real part must equal the column f(W/s)|site> of the dense reference."""
    W, data, fW = benchmark_reference
    inst = build.Instance(**build.BENCHMARK, variant=variant)
    assert inst.n_qubits == 18
    V = inst.V
    qc = build.hadamard_test_circuit(inst, data["phases_reflection"], site=site, final_h=False)
    b0, b1 = _hadamard_branches(inst, qc, site)
    e = np.zeros_like(b0)
    e[site] = 1.0
    assert np.abs(b0 - e).max() < TOL
    assert np.abs(b1[:V].real - fW[:, site]).max() < TOL
    assert abs(b1[site].real - fW[site, site]) < TOL
    assert abs(np.linalg.norm(b1) - 1) < TOL

    # the measured quantity, with the final Hadamard, and the log det it reconstructs
    psi = simulate(build.hadamard_test_circuit(inst, data["phases_reflection"], site=site))
    z = 1 - 2 * float(np.sum(np.abs(psi[2 ** 17:]) ** 2))
    assert abs(z - np.trace(fW) / V) < TOL                   # free operator: every site gives the trace
    lam_min, s = np.linalg.eigvalsh(W)[0], inst.sub
    logdet = V * math.log(s / lam_min) * (z - data["shift"]) + V * math.log(s)
    exact = np.linalg.slogdet(W)[1]
    mdl = count.model_items()
    assert abs(exact - mdl["logdet_W"]) < 1e-9
    assert abs(logdet - exact) / exact < 0.01                # the chapter's 1% polynomial guarantee
    assert abs(z - mdl["normalized_trace_shifted"]) < mdl["poly_tolerance"]


def test_full_benchmark_check_catches_a_wrong_phase_and_a_dropped_gate(benchmark_reference):
    """Two deliberate faults in the full circuit, each of which the column check must see:
    one QET phase with the wrong sign, and one Toffoli of one O_c removed."""
    from qiskit import QuantumCircuit
    W, data, fW = benchmark_reference
    inst = build.Instance(**build.BENCHMARK, variant="merged")
    V, site = inst.V, 5

    def column_error(qc):
        _, b1 = _hadamard_branches(inst, qc, site)
        return float(np.abs(b1[:V].real - fW[:, site]).max())

    phases = np.array(data["phases_reflection"])
    bad = phases.copy()
    bad[40] = -bad[40]
    assert column_error(build.hadamard_test_circuit(inst, bad, site=site, final_h=False)) > 1e-3

    good = build.hadamard_test_circuit(inst, phases, site=site, final_h=False)
    marks = [m for m in build.block_marks(good) if m[0] == "O_c"]
    name, start, stop = marks[len(marks) // 2]
    victim = next(i for i in range(start, stop) if good.data[i].operation.name == "ccx")
    dropped = QuantumCircuit(*good.qregs)
    for i, g in enumerate(good.data):
        if i != victim:
            dropped.append(g.operation, g.qubits)
    assert column_error(dropped) > 1e-3
    assert column_error(good) < TOL


# --------------------------------------------------------------------------- #
# Counting
# --------------------------------------------------------------------------- #

def test_gate_set_and_counts_are_walked_from_the_circuit():
    mdl = count.model_items()
    for variant in build.VARIANTS:
        qc = build.benchmark_circuit(variant)
        assert set(qc.count_ops()) <= {"h", "x", "s", "sdg", "cx", "ccx", "rz"}
        w = count.walk_counts(qc)
        assert w["total"]["toffoli"] == qc.count_ops()["ccx"]
        assert w["total"]["rot"] + w["total"]["t"] <= qc.count_ops()["rz"]
        assert sum(b["toffoli"] for b in w["blocks"].values()) == w["total"]["toffoli"]
        assert sum(b["rot"] for b in w["blocks"].values()) == w["total"]["rot"]
        assert w["qubits_allocated"] == 18
    lit = count.walk_counts(build.benchmark_circuit("literal"))
    mer = count.walk_counts(build.benchmark_circuit("merged"))
    d = mdl["d"]
    # per-block counts that follow from the construction, with the model's ladder formula
    assert lit["blocks"]["O_c"]["toffoli"] == d * 8 * mcx_toffoli_count(4, "linear")
    assert mer["blocks"]["O_c"]["toffoli"] == d * 4 * mcx_toffoli_count(3, "linear")
    assert lit["blocks"]["proj_and"]["toffoli"] == d * 8 and lit["blocks"]["proj_rot"]["rot"] == d * 2
    assert lit["blocks"]["O_A"]["rot"] == mer["blocks"]["O_A"]["rot"] == d * 2
    # Totals pinned by formula, so they do not rest on counts.json (which this same code regenerates):
    # no Toffoli outside O_c and the projector AND, no arbitrary rotation outside O_A, the phases and
    # the folded last phase, and no exact T gate anywhere.
    com = count.walk_counts(build.benchmark_circuit("compact"))
    expect = {"literal": (d * (8 * mcx_toffoli_count(4, "linear") + 8), 4 * d),      # 4032, 336
              "merged": (d * (4 * mcx_toffoli_count(3, "linear") + 8), 4 * d),       # 1680, 336
              "compact": (d * 6 + (d - 1) * 8, 4 * d - 1)}                           # 1168, 335
    for variant, w in (("literal", lit), ("merged", mer), ("compact", com)):
        assert (w["total"]["toffoli"], w["total"]["rot"], w["total"]["t"]) == (*expect[variant], 0)
        assert all(b["toffoli"] == 0 for name, b in w["blocks"].items() if name not in ("O_c", "proj_and"))
        assert all(b["rot"] == 0 for name, b in w["blocks"].items() if name not in ("O_A", "proj_rot", "had"))


def test_low_coordinate_bits_are_never_acted_on():
    """Free W couples only sites two hops apart, so the low bit of every coordinate is a spectator:
    the benchmark circuit touches 14 of its 18 qubits. Reference side: W at V = 4^4 is 16 identical
    blocks, each (m0^2 + 2K^2) - (K^2/2) sum_mu X_mu on the four high coordinate bits."""
    for variant in build.VARIANTS:
        w = count.walk_counts(build.benchmark_circuit(variant, site=255))
        assert w["registers"]["site"] == {"allocated": 8, "touched_by_a_gate": 8}   # the prep X gates
        w = count.walk_counts(build.benchmark_circuit(variant, site=0))
        assert w["registers"]["site"] == {"allocated": 8, "touched_by_a_gate": 4}
        assert w["qubits_allocated"] == 18 and w["qubits_touched_by_a_gate"] == 14
    L, D, m0, K = 4, 4, 0.4, 1.0
    W = ref_W(L, D, m0, K)
    X, I2 = np.array([[0.0, 1.0], [1.0, 0.0]]), np.eye(2)
    H4 = (m0 ** 2 + 2 * K ** 2) * np.eye(16)
    for mu in range(D):
        term = np.array([[1.0]])
        for nu in reversed(range(D)):
            term = np.kron(term, X if nu == mu else I2)
        H4 -= 0.5 * K ** 2 * term
    high = [site_index([2 * ((h >> mu) & 1) for mu in range(D)], L) for h in range(16)]
    for low in itertools.product(range(2), repeat=D):
        idx = [n + site_index(low, L) for n in high]
        assert np.abs(W[np.ix_(idx, idx)] - H4).max() < 1e-13
        rest = [n for n in range(L ** D) if n not in idx]
        assert np.abs(W[np.ix_(idx, rest)]).max() == 0


@pytest.mark.parametrize("L,D", [(4, 4), (8, 2), (8, 4), (16, 1), (16, 2)])
def test_shift_oracle_toffolis_follow_the_models_formula(L, D):
    """Literal O_c: 2D controlled shifts, each the model's shift2_toffolis(L, column qubits)."""
    from estimates.ch04_hybrid_lqcd import shift2_toffolis
    inst = build.Instance(L=L, D=D, variant="literal")
    w = count.walk_counts(build.shift_oracle_circuit(inst))
    assert w["total"]["toffoli"] == 2 * D * shift2_toffolis(L, inst.n_col)


def test_counts_json_is_current():
    stored = json.loads((HERE / "counts.json").read_text())
    mdl = count.model_items()
    for variant in build.VARIANTS:
        fresh = count.count_variant(variant, mdl)
        for key in ("logical_qubits", "qubits_touched_by_a_gate", "n_toffoli", "n_rot", "n_t_fixed"):
            assert stored["variants"][variant][key] == fresh[key]
        a = stored["variants"][variant]["t_report_convention_model_eps_syn"]["t_total"]
        assert abs(a - fresh["t_report_convention_model_eps_syn"]["t_total"]) < 1e-6
        assert abs(stored["variants"][variant]["t_depth"]["commuting"]["t_depth"]
                   - fresh["t_depth"]["commuting"]["t_depth"]) < 1e-6


def test_rz_classification_and_generic_instance_rotation_count():
    assert count.rz_class(math.pi / 2) == "clifford" and count.rz_class(-math.pi) == "clifford"
    assert count.rz_class(math.pi / 4) == "t" and count.rz_class(-3 * math.pi / 4) == "t"
    assert count.rz_class(0.3) == "rot"
    # at generic s the literal O_A has 4 rotations per query, the multiplexed one 2 (the model's 2-4)
    for variant, n in (("literal", 4), ("merged", 2)):
        inst = build.Instance(L=4, D=4, m0=0.7, K=0.9, s=6.5, variant=variant)
        w = count.walk_counts(build.block_encoding_circuit(inst))
        assert w["blocks"]["O_A"]["rot"] == n


def test_t_depth_walker_on_small_circuits():
    from qiskit import QuantumCircuit
    qc = QuantumCircuit(5)
    qc.metadata = {"blocks": [["x", 0, 2]]}
    qc.ccx(0, 1, 2)
    qc.ccx(0, 3, 4)                                          # shares only a control with the first
    assert count.t_depth(qc, 20.0, "commuting")["t_depth"] == 1
    assert count.t_depth(qc, 20.0, "as_written")["t_depth"] == 2
    qc = QuantumCircuit(4)
    qc.ccx(0, 1, 2)
    qc.rz(0.3, 2)                                            # diagonal on the Toffoli's target: must wait
    qc.ccx(2, 1, 3)                                          # control on qubit 2 commutes with the rz
    assert count.t_depth(qc, 20.0, "commuting")["t_depth"] == 21
    assert count.t_depth(qc, 20.0, "as_written")["t_depth"] == 22
    qc.ccx(0, 1, 2)                                          # target on qubit 2 again: after the rz
    assert count.t_depth(qc, 20.0, "commuting")["t_depth"] == 22
    qc = QuantumCircuit(2)
    qc.rz(0.3, 0)
    qc.cx(0, 1)
    qc.rz(math.pi / 4, 0)                                    # exact T, commutes with the control
    qc.h(0)
    qc.rz(0.2, 0)
    assert count.t_depth(qc, 20.0, "commuting")["t_depth"] == 40


def test_qasm_export_round_trips():
    """The exported OpenQASM 3 file parses back to the same gate sequence (names, qubits, angles)."""
    qasm3 = pytest.importorskip("qiskit.qasm3")
    pytest.importorskip("qiskit_qasm3_import")
    for variant in build.VARIANTS:
        path = HERE / f"logdet_4x4x4x4_d84_{variant}.qasm"
        loaded = qasm3.loads(path.read_text())
        built = build.benchmark_circuit(variant, measure=True)
        assert loaded.num_qubits == built.num_qubits and len(loaded.data) == len(built.data)
        for a, b in zip(loaded.data, built.data):
            assert a.operation.name == b.operation.name
            assert [loaded.find_bit(q).index for q in a.qubits] == [built.find_bit(q).index for q in b.qubits]
            assert np.allclose([float(p) for p in a.operation.params], [float(p) for p in b.operation.params],
                               rtol=0, atol=1e-12)
