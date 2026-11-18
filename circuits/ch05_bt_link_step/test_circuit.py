"""Correctness checks of the 2T link-step circuits against references built with numpy from the group and
Hamiltonian definitions, not from the circuit code.

References used here:
  * the 24 group elements as 2x2 matrices and their multiplication, inverse and trace (numpy products),
  * the link electric term as the Cayley-graph Laplacian sum_{gamma} (1 - L_gamma) built from the
    multiplication table, exponentiated with scipy.linalg.expm,
  * the plaquette term exp(i theta Re Tr(g1 g2 g3^-1 g4^-1)) from matrix products.

Tolerance: TOL = 1e-9 on every amplitude after removing one global phase (double precision; the observed
errors are below 1e-13).

Sizes: primitives on 5 to 10 qubits (full unitary), the plaquette on 20 qubits, the one-plaquette Trotter
step on 22 qubits (20 data + 2 ancilla, statevector on a random gauge-register state). The magnetic layer
of the full 4x4 benchmark lattice (160 qubits) is a permutation-plus-phase circuit and is checked on random
basis states with a bit-level simulator.

Run:  python -m pytest circuits/ch05_bt_link_step/test_circuit.py -q
"""
import json
import math
import os
import sys

import numpy as np
import pytest
import scipy.linalg as sla
from qiskit import QuantumCircuit
from qiskit.quantum_info import Operator

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import build as B  # noqa: E402

TOL = 1e-9
BT = B.group()
IDX = BT.index                      # register values of the 24 elements
MATS = np.array([BT.mats[i] for i in IDX])


# ---------------------------------------------------------------------------------------------------
# Independent references
# ---------------------------------------------------------------------------------------------------

def ref_mul(a: int, b: int) -> int:
    return BT.find(BT.mats[a] @ BT.mats[b])


def ref_inv(a: int) -> int:
    return BT.find(np.linalg.inv(BT.mats[a]))


def ref_laplacian(normalized: bool = True) -> np.ndarray:
    """sum over the 8 elements with Re Tr = 1 of (1 - L_gamma), L_gamma|h> = |gamma h>, as a 24x24 matrix
    in the order IDX, built from matrix products only."""
    lap = np.zeros((24, 24))
    for gi in IDX:
        if abs(np.trace(BT.mats[gi]).real - 1.0) > 1e-9:
            continue
        for hpos, hi in enumerate(IDX):
            lap[hpos, hpos] += 1.0
            lap[BT.pos[ref_mul(gi, hi)], hpos] -= 1.0
    assert np.allclose(lap, lap.T)
    if normalized:
        # the faithful 2-dimensional irrep carries 3/4: its eigenvalue is |Gamma| - sum chi/2 = 8 - 8/2 = 4
        lap = lap * (0.75 / 4.0)
    return lap


def ref_plaquette_retr() -> np.ndarray:
    """Re Tr(g_a g_b g_c^-1 g_d^-1) for all 24^4 element positions."""
    md = MATS.conj().transpose(0, 2, 1)
    return np.einsum("aij,bjk,ckl,dli->abcd", MATS, MATS, md, md).real


def aligned_error(out: np.ndarray, ref: np.ndarray) -> float:
    """max |out - ref| after removing one global phase."""
    ph = np.vdot(ref.ravel(), out.ravel())
    assert abs(ph) > 1e-6, "no overlap with the reference"
    return float(np.abs(out * np.conj(ph / abs(ph)) - ref).max())


def embed(psi: np.ndarray, n_anc: int) -> np.ndarray:
    """psi[a0, a1, a2, a3] over element positions of links 0..3 -> full statevector (ancillas |0>)."""
    full = np.zeros((2 ** n_anc, 32, 32, 32, 32), dtype=complex)
    full[0][np.ix_(IDX, IDX, IDX, IDX)] = psi.transpose(3, 2, 1, 0)
    return full.ravel()


def extract(vec: np.ndarray, n_anc: int) -> tuple[np.ndarray, float]:
    """Inverse of embed; also returns the norm found outside the gauge register / ancilla-zero subspace."""
    full = vec.reshape((2 ** n_anc, 32, 32, 32, 32)).copy()
    psi = full[0][np.ix_(IDX, IDX, IDX, IDX)].transpose(3, 2, 1, 0).copy()
    full[0][np.ix_(IDX, IDX, IDX, IDX)] = 0.0          # what is left is outside the subspace
    return psi, float(np.linalg.norm(full))


def run_statevector(qc: QuantumCircuit, psi0: np.ndarray) -> np.ndarray:
    from qiskit_aer import AerSimulator
    sim = AerSimulator(method="statevector", max_parallel_threads=4)
    full = QuantumCircuit(qc.num_qubits)
    full.set_statevector(psi0)
    full.compose(qc, inplace=True)
    full.save_statevector()
    return np.asarray(sim.run(full).result().get_statevector(full))


def random_state(seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    psi = rng.normal(size=(24,) * 4) + 1j * rng.normal(size=(24,) * 4)
    return psi / np.linalg.norm(psi)


def ref_step(psi, theta_e, theta_b, order=1):
    e = sla.expm(-1j * theta_e * ref_laplacian())

    def apply_e(x):
        for axis in range(4):                           # exp(-i theta E^2) on each link in turn
            x = np.moveaxis(np.tensordot(e, x, axes=([1], [axis])), 0, axis)
        return x

    retr = ref_plaquette_retr()
    if order == 1:
        return np.exp(1j * theta_b * retr) * apply_e(psi)
    half = np.exp(0.5j * theta_b * retr)
    return half * apply_e(half * psi)


# ---------------------------------------------------------------------------------------------------
# Group and register
# ---------------------------------------------------------------------------------------------------

def test_register_matches_the_paper_example():
    # arXiv:2208.12309: (1/2)[[eta, eta*], [-eta, eta*]] = (-1) i j l^2 -> |10111> = |23>
    eta = 1 + 1j
    assert BT.find(0.5 * np.array([[eta, np.conj(eta)], [-eta, np.conj(eta)]])) == 23
    assert len(IDX) == 24 and all(not ((i >> 3) & 1 and (i >> 4) & 1) for i in IDX)
    for a in IDX:                       # closure
        for b in IDX:
            ref_mul(a, b)
    traces = sorted(round(np.trace(m).real) for m in MATS)
    assert [traces.count(v) for v in (-2, -1, 0, 1, 2)] == [1, 8, 6, 8, 1]


# ---------------------------------------------------------------------------------------------------
# Primitives: full unitaries
# ---------------------------------------------------------------------------------------------------

def test_u_inv_on_every_element():
    qc = QuantumCircuit(5)
    B.u_inv(qc, list(range(5)))
    u = Operator(qc).data
    for a in IDX:
        assert abs(u[ref_inv(a), a] - 1) < TOL
    assert B.count_gates(qc)["toffoli"] == 4 and B.count_gates(qc)["rot"] == 0


def test_u_mul_on_every_pair():
    qc = QuantumCircuit(10)
    B.u_mul(qc, list(range(5)), list(range(5, 10)))
    u = Operator(qc).data
    for a in IDX:
        for b in IDX:
            assert abs(u[a + 32 * ref_mul(a, b), a + 32 * b] - 1) < TOL
    assert B.count_gates(qc) == {"toffoli": 8, "rot": 0, "t": 0, "clifford": 11}


def test_u_mul_dagger_undoes_u_mul():
    qc = QuantumCircuit(10)
    B.u_mul_dagger(qc, list(range(5)), list(range(5, 10)))
    u = Operator(qc).data
    for a in IDX:
        for b in IDX:
            assert abs(u[a + 32 * b, a + 32 * ref_mul(a, b)] - 1) < TOL


@pytest.mark.parametrize("theta", [0.37, -1.9])
def test_u_tr_phase(theta):
    qc = QuantumCircuit(5)
    B.u_tr(qc, list(range(5)), theta)
    u = Operator(qc).data
    assert np.abs(u - np.diag(np.diag(u))).max() < TOL
    ref = np.exp(1j * theta * np.array([np.trace(BT.mats[i]).real for i in IDX]))
    assert aligned_error(np.diag(u)[IDX], ref) < TOL
    assert B.count_gates(qc)["rot"] == 8


def electric_block(theta_e, f=None):
    qc = QuantumCircuit(7)
    B.electric_link(qc, list(range(5)), 5, 6, theta_e, f)
    u = Operator(qc).data
    return qc, u[np.ix_(IDX, IDX)]          # ancillas |00> in and out


@pytest.mark.parametrize("theta_e", [0.83, -0.21, 2.7])
def test_electric_link_is_the_exponential_of_the_cayley_laplacian(theta_e):
    qc, sub = electric_block(theta_e)
    # the 24x24 block is unitary: nothing leaves the gauge register and both ancillas return to |0>
    assert np.abs(sub.conj().T @ sub - np.eye(24)).max() < TOL
    assert aligned_error(sub, sla.expm(-1j * theta_e * ref_laplacian())) < TOL
    assert B.count_gates(qc) == {"toffoli": 8, "rot": 9, "t": 20, "clifford": 172}


def test_electric_eigenvalues_are_the_laplacian_spectrum():
    ev = np.linalg.eigvalsh(ref_laplacian(normalized=False))
    vals, mult = np.unique(np.round(ev, 9), return_counts=True)
    f = B.electric_eigenvalues(normalize=False)
    expect = {}
    for name, v in f.items():
        dim = int(name[0])
        expect[round(v, 9)] = expect.get(round(v, 9), 0) + dim * dim
    assert {float(v): int(m) for v, m in zip(vals, mult)} == expect
    assert abs(B.electric_eigenvalues()["2_0"] - 0.75) < 1e-12


def test_u_fourier_block_diagonalizes_the_electric_term():
    qc = QuantumCircuit(6)
    B.u_fourier(qc, list(range(5)), 5)
    v = Operator(qc).data
    lap32 = np.zeros((32, 32))
    lap32[np.ix_(IDX, IDX)] = ref_laplacian()
    cols = v[:32, :32]                       # ancilla |0> in and out
    assert np.abs(cols.conj().T @ cols - np.eye(32)).max() < TOL
    d = cols @ lap32 @ cols.conj().T
    assert np.abs(d - np.diag(np.diag(d))).max() < TOL
    assert B.count_gates(qc) == {"toffoli": 2, "rot": 2, "t": 10, "clifford": 77}


# ---------------------------------------------------------------------------------------------------
# Plaquette and Trotter step: statevector on the gauge register
# ---------------------------------------------------------------------------------------------------

@pytest.mark.parametrize("uncompute", ["lly", "dagger"])
def test_plaquette_on_20_qubits(uncompute):
    theta_b = 0.61
    qc = QuantumCircuit(20)
    B.plaquette(qc, *[B.link_qubits(k) for k in range(4)], theta_b, uncompute)
    psi = random_state(11)
    out, leak = extract(run_statevector(qc, embed(psi, 0)), 0)
    assert leak < TOL
    assert aligned_error(out, np.exp(1j * theta_b * ref_plaquette_retr()) * psi) < TOL
    c = B.count_by_block(qc)
    assert c["U_mul"]["calls"] == 6 and c["U_Tr"]["calls"] == 1
    assert c["U_inv"]["calls"] == (6 if uncompute == "lly" else 4)


@pytest.mark.parametrize("order", [1, 2])
def test_single_plaquette_trotter_step_on_22_qubits(order):
    theta_e, theta_b = B.couplings(1.3, 0.225)
    qc = B.single_plaquette_step(theta_e, theta_b, steps=1, order=order)
    assert qc.num_qubits == 22
    psi = random_state(5 + order)
    out, leak = extract(run_statevector(qc, embed(psi, 2)), 2)
    assert leak < TOL
    assert aligned_error(out, ref_step(psi, theta_e, theta_b, order)) < TOL


def test_two_steps_compose():
    theta_e, theta_b = 0.4, 0.7
    qc = B.single_plaquette_step(theta_e, theta_b, steps=2, order=1)
    psi = random_state(21)
    out, leak = extract(run_statevector(qc, embed(psi, 2)), 2)
    ref = ref_step(ref_step(psi, theta_e, theta_b), theta_e, theta_b)
    assert leak < TOL and aligned_error(out, ref) < TOL


# ---------------------------------------------------------------------------------------------------
# Tests that must fail for a wrong circuit
# ---------------------------------------------------------------------------------------------------

def test_wrong_sign_of_the_magnetic_angle_is_detected():
    theta_b = 0.61
    qc = QuantumCircuit(20)
    B.plaquette(qc, *[B.link_qubits(k) for k in range(4)], -theta_b)
    psi = random_state(11)
    out, _ = extract(run_statevector(qc, embed(psi, 0)), 0)
    assert aligned_error(out, np.exp(1j * theta_b * ref_plaquette_retr()) * psi) > 1e-3


def test_dropped_inversion_is_detected():
    # the plaquette with the inversion of link 3 left out computes Tr(g2^-1 g3 g0 g1)
    theta_b = 0.61
    l = [B.link_qubits(k) for k in range(4)]
    qc = QuantumCircuit(20)
    B.u_inv(qc, l[2])
    B.u_mul(qc, l[0], l[1]); B.u_mul(qc, l[3], l[1]); B.u_mul(qc, l[2], l[1])
    B.u_tr(qc, l[1], theta_b)
    B.u_mul_dagger(qc, l[2], l[1]); B.u_mul_dagger(qc, l[3], l[1]); B.u_mul_dagger(qc, l[0], l[1])
    B.u_inv(qc, l[2])
    psi = random_state(11)
    out, leak = extract(run_statevector(qc, embed(psi, 0)), 0)
    assert leak < TOL                                   # still a valid circuit, but the wrong operator
    assert aligned_error(out, np.exp(1j * theta_b * ref_plaquette_retr()) * psi) > 1e-3


def test_dropped_multiplication_is_detected():
    qc = QuantumCircuit(10)
    g, h = list(range(5)), list(range(5, 10))
    B.u_mul(qc, g, h)
    broken = QuantumCircuit(10)
    for k, inst in enumerate(qc.data):
        if k != 2:                                      # drop the Toffoli of the phi step
            broken.append(inst.operation, inst.qubits)
    u = Operator(broken).data
    assert any(abs(u[a + 32 * ref_mul(a, b), a + 32 * b] - 1) > 0.5 for a in IDX for b in IDX)


def test_wrong_electric_eigenvalue_is_detected():
    f = dict(B.electric_eigenvalues())
    f["3"] = 0.0                                        # the three-dimensional irrep dropped
    _, sub = electric_block(0.83, f)
    assert aligned_error(sub, sla.expm(-0.83j * ref_laplacian())) > 1e-3


# ---------------------------------------------------------------------------------------------------
# Full benchmark lattice
# ---------------------------------------------------------------------------------------------------

def basis_state_run(qc: QuantumCircuit, bits: np.ndarray) -> tuple[np.ndarray, complex]:
    """Run a permutation-plus-phase circuit on one computational basis state."""
    bits = bits.copy()
    phase = 1.0 + 0j
    idx = {q: i for i, q in enumerate(qc.qubits)}
    for inst in qc.data:
        name = inst.operation.name
        q = [idx[x] for x in inst.qubits]
        if name == "x":
            bits[q[0]] ^= 1
        elif name == "cx":
            bits[q[1]] ^= bits[q[0]]
        elif name == "ccx":
            bits[q[2]] ^= bits[q[0]] & bits[q[1]]
        elif name == "swap":
            bits[q[0]], bits[q[1]] = bits[q[1]], bits[q[0]]
        elif name == "cz":
            phase *= -1 if bits[q[0]] and bits[q[1]] else 1
        elif name in ("z", "s", "sdg", "t", "tdg"):
            if bits[q[0]]:
                phase *= {"z": -1, "s": 1j, "sdg": -1j, "t": np.exp(0.25j * np.pi),
                          "tdg": np.exp(-0.25j * np.pi)}[name]
        elif name == "rz":
            th = float(inst.operation.params[0])
            phase *= np.exp(0.5j * th * (1 if bits[q[0]] else -1))
        else:
            raise ValueError(f"{name} is not a permutation or a phase")
    return bits, phase


def ref_lattice_plaquette_sum(L: int, conf: dict) -> float:
    """sum over sites of Re Tr(U_x(s) U_y(s+x) U_x(s+y)^-1 U_y(s)^-1), conf[(x, y, mu)] = register value."""
    tot = 0.0
    for x in range(L):
        for y in range(L):
            m = (BT.mats[conf[(x, y, 0)]] @ BT.mats[conf[((x + 1) % L, y, 1)]]
                 @ BT.mats[conf[(x, (y + 1) % L, 0)]].conj().T @ BT.mats[conf[(x, y, 1)]].conj().T)
            tot += np.trace(m).real
    return tot


@pytest.mark.parametrize("L,uncompute", [(2, "lly"), (4, "lly"), (4, "dagger")])
def test_magnetic_layer_of_the_periodic_lattice(L, uncompute):
    theta_b = 0.45
    qc = B.magnetic_layer(L, theta_b, uncompute)
    links = B.lattice_links(L)
    rng = np.random.default_rng(100 + L)
    phases, sums = [], []
    for _ in range(6):
        conf = {k: int(rng.choice(IDX)) for k in links}
        bits = np.zeros(qc.num_qubits, dtype=np.int8)
        for k, li in links.items():
            for b in range(5):
                bits[5 * li + b] = (conf[k] >> b) & 1
        out, ph = basis_state_run(qc, bits)
        assert np.array_equal(out, bits)                 # every link register restored
        phases.append(ph)
        sums.append(ref_lattice_plaquette_sum(L, conf))
    assert len(set(np.round(sums, 6))) > 1
    for ph, s in zip(phases[1:], sums[1:]):              # ratios remove the global phase
        assert abs(ph / phases[0] - np.exp(1j * theta_b * (s - sums[0]))) < TOL


def test_benchmark_structure():
    qc = B.benchmark_circuit("first")
    assert qc.num_qubits == 176 and B.used_qubits(qc) == 176
    c = B.count_by_block(qc)
    assert {k: v["calls"] for k, v in c.items()} == {"U_F": 64, "U_phase": 32, "U_inv": 96, "U_mul": 96,
                                                     "U_Tr": 16}
    # every link gets exactly one electric block, every plaquette of the 4x4 torus appears once
    assert sorted(B.lattice_plaquettes(4)) == sorted(set(B.lattice_plaquettes(4)))
    assert len(B.lattice_plaquettes(4)) == 16
    hits = [0] * 32
    for pl in B.lattice_plaquettes(4):
        for k in pl:
            hits[k] += 1
    assert hits == [2] * 32
    full = B.benchmark_circuit("full")
    assert all(v == 6 * c[k]["calls"] for k, v in
               {k: v["calls"] for k, v in B.count_by_block(full).items()}.items())


def test_counts_json_is_the_walk_of_the_circuit():
    path = os.path.join(HERE, "counts.json")
    if not os.path.exists(path):
        pytest.skip("counts.json not written yet (run count.py)")
    data = json.load(open(path))
    for tier in ("first", "full"):
        c = B.count_gates(B.benchmark_circuit(tier))
        row = data["circuit"][tier]
        assert (row["n_toffoli"], row["n_rot"], row["n_t_fixed"]) == (c["toffoli"], c["rot"], c["t"])


def test_qasm_export_reloads_to_the_same_gates():
    from qiskit import qasm3
    path = os.path.join(HERE, "bt_single_plaquette_1step.qasm")
    if not os.path.exists(path):
        pytest.skip("QASM not exported yet (run build.py)")
    back = qasm3.loads(open(path).read())
    te, tb = B.couplings(B.G2_PLACEHOLDER, B.BENCHMARK["dt_first_fm"] / B.BENCHMARK["a_s_fm"])
    qc = B.single_plaquette_step(te, tb)
    assert [i.operation.name for i in back.data] == [i.operation.name for i in qc.data]
    assert B.count_gates(back) == B.count_gates(qc)


# ---------------------------------------------------------------------------------------------------
# The exported circuits against a checker that shares no code with build.py
# ---------------------------------------------------------------------------------------------------
# verify_independent.py builds 2T from quaternion arithmetic, reads the QASM text with its own parser and
# runs it with its own simulators. It cuts the lattice circuit into contiguous blocks that exhaust the gate
# list and checks every block on the whole physical domain of the qubits it touches: each electric block
# as a dense unitary against exp(-i theta_e Laplacian), each plaquette block on all 24^4 configurations.
# The text checked is a fresh export of the generator, and the file on disk must be that same text.

def _fresh_export(qc, tmp_path, name):
    from qiskit import qasm3
    text = qasm3.dumps(qc)
    path = tmp_path / name
    path.write_text(text)
    exported = os.path.join(HERE, name)
    if os.path.exists(exported):
        assert open(exported).read() == text, f"{name} on disk is not the generator's output (run build.py)"
    return str(path)


@pytest.mark.parametrize("tier", ["first", "full"])
def test_benchmark_circuit_block_by_block_with_the_independent_checker(tier, tmp_path):
    import verify_independent as V
    V.check_group()
    name, dt, steps = V.EXPORTS[tier]
    qc = B.benchmark_circuit(tier)
    path = _fresh_export(qc, tmp_path, name)
    r = V.check_lattice_file(path, V.L_BENCH, *V.angles(dt), steps)
    assert r["blocks"] == steps * (32 + 16)
    c = V.counts(*V.parse(path))                       # counted from the text, not with build.count_gates
    walk = B.count_gates(qc)
    assert (c["toffoli"], c["rot"], c["t"]) == (walk["toffoli"], walk["rot"], walk["t"])
    assert c["qubits_used"] == c["qubits_declared"] == 176


def test_single_plaquette_step_with_the_independent_checker(tmp_path):
    import verify_independent as V
    te, tb = B.couplings(B.G2_PLACEHOLDER, B.BENCHMARK["dt_first_fm"] / B.BENCHMARK["a_s_fm"])
    path = _fresh_export(B.single_plaquette_step(te, tb), tmp_path, "bt_single_plaquette_1step.qasm")
    r = V.check_single_plaquette_file(path, *V.angles(0.045))
    assert r["error"] < TOL and r["norm_outside"] < TOL
