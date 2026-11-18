"""Checks of the Ch. 2 2028 circuit against references built with numpy from the Hamiltonian definition.

    pytest test_circuit.py -q

The reference (class Ref, dense_hamiltonian, propagate_basis) never calls build.py. It works on the coordinate
tensor psi[x_0, ..., x_{nc-1}] with x in {0, 1, 2}, builds the hopping matrix and the contact potential from their
definitions, and applies the same sequence of exponentials the chapter specifies.

Tolerances
    TOL = 1e-9    maximum absolute difference of any amplitude or matrix element, circuit against reference
                  (double precision; the longest simulated circuit has 10,491 gates)
    NEG = 1e-3    a deliberately wrong reference or a mutated circuit must differ by more than this
                  (2-norm of the statevector difference, or a single phase for the basis-state check)

Sizes (the full instance has 36 qubits and is not simulated densely)
    every block type at the full-instance size: coordinate basis change and kinetic block (2 to 4 qubits),
        current insertion (3 qubits), pair contact (18 qubits) and triple contact (30 qubits) on every basis state
    complete Hadamard-test circuits: 4 nucleons on 3^1 (12 qubits), 2 nucleons on 3^3 (18 qubits),
        3 nucleons on 3^2 (20 qubits), 4 nucleons on 3^2 (24 qubits)
"""
from __future__ import annotations

import itertools
import json
import math
import os
import sys
from dataclasses import replace
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")

import numpy as np
import pytest
from scipy.linalg import expm
from qiskit import QuantumCircuit
from qiskit.quantum_info import Operator

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import build as B                                      # noqa: E402
import count as CNT                                    # noqa: E402

TOL = 1e-9
NEG = 1e-3
T, C0, D0 = 10.5794, -98.2265511, 127.839693           # app01 Eq. after "Concrete instantiations"
Q = 2 * math.pi / 3


# --------------------------------------------------------------------------- #
# Independent reference
# --------------------------------------------------------------------------- #

def hop_matrix(t):
    """One periodic three-site coordinate: amplitude -t to each of the two other sites."""
    K = np.zeros((3, 3))
    for x in range(3):
        for d in (+1, -1):
            K[(x + d) % 3, x] += -t
    return K


class Ref:
    """State tensor psi[x_0, ..., x_{nc-1}], coordinate index c = particle * dim + axis."""

    def __init__(self, dim, n, t=T, c0=C0, d0=D0):
        self.dim, self.n, self.t = dim, n, t
        self.nc = n * dim
        self.shape = (3,) * self.nc
        idx = np.indices(self.shape)
        self.x = idx
        site = [sum(idx[f * dim + a] * 3 ** a for a in range(dim)) for f in range(n)]
        V = np.zeros(self.shape)
        for a, b in itertools.combinations(range(n), 2):
            V += c0 * (site[a] == site[b])
        for a, b, c in itertools.combinations(range(n), 3):
            V += d0 * ((site[a] == site[b]) & (site[b] == site[c]))
        self.V = V
        self.K1 = hop_matrix(t)
        self.hop_constant = 2 * dim * t * n

    def uniform(self):
        return np.full(self.shape, 3.0 ** (-self.nc / 2), dtype=complex)

    def kin(self, psi, dt):
        U = expm(-1j * dt * self.K1)
        for ax in range(self.nc):
            psi = np.moveaxis(np.tensordot(U, psi, axes=(1, ax)), 0, ax)
        return psi

    def pot(self, psi, dt):
        return np.exp(-1j * dt * self.V) * psi

    def prep(self, layers):
        psi = self.uniform()
        for dt_v, dt_k in layers:
            psi = self.kin(self.pot(psi, dt_v), dt_k)
        return psi

    def trotter(self, psi, tau, steps):
        if steps == 0:
            return psi
        dt = tau / steps
        psi = self.kin(psi, dt / 2)
        for j in range(steps):
            psi = self.pot(psi, dt)
            psi = self.kin(psi, dt if j < steps - 1 else dt / 2)
        return psi

    def cur(self, psi, f, q, sign=+1):
        return np.exp(1j * sign * q * self.x[f * self.dim]) * psi

    def w_psi(self, psi, a, b, tau, steps, q, e0):
        """e^{i (E0 - const) tau} U_a^dagger S(tau) U_b psi, the branch the ancilla controls."""
        out = self.cur(self.trotter(self.cur(psi, b, q), tau, steps), a, q, -1)
        return np.exp(1j * (e0 - self.hop_constant) * tau) * out

    def h_apply(self, psi):
        out = (self.hop_constant + self.V) * psi
        for ax in range(self.nc):
            out = out - self.t * (np.roll(psi, 1, ax) + np.roll(psi, -1, ax))
        return out


def dense_hamiltonian(dim, n, t=T, c0=C0, d0=D0):
    """H as a matrix from explicit matrix elements (a second construction, used to check Ref itself)."""
    nc = n * dim
    states = list(itertools.product(range(3), repeat=nc))
    index = {s: i for i, s in enumerate(states)}
    H = np.zeros((len(states),) * 2)
    for s, i in index.items():
        sites = [tuple(s[f * dim:(f + 1) * dim]) for f in range(n)]
        H[i, i] += 2 * dim * t * n
        H[i, i] += c0 * sum(sites[a] == sites[b] for a, b in itertools.combinations(range(n), 2))
        H[i, i] += d0 * sum(sites[a] == sites[b] == sites[c] for a, b, c in itertools.combinations(range(n), 3))
        for c in range(nc):
            for d in (+1, -1):
                s2 = list(s); s2[c] = (s2[c] + d) % 3
                H[index[tuple(s2)], i] += -t
    return H, states


def legal_index(nc):
    """Position-register integer of each coordinate tuple: two bits per coordinate, coordinate 0 lowest."""
    idx = np.indices((3,) * nc)
    return sum(idx[c] * 4 ** c for c in range(nc))


def embed(tensors, inst):
    """Full statevector with workspace |0>, tensors = (anc 0 branch, anc 1 branch)."""
    out = np.zeros(2 ** inst.n_qubits, dtype=complex)
    li = legal_index(inst.n_coords).ravel()
    for anc, tens in enumerate(tensors):
        if tens is not None:
            out[li + (anc << inst.anc)] = tens.ravel()
    return out


def expected_final(inst):
    """Statevector the Hadamard test must end in, before measurement."""
    ref = Ref(inst.dim, inst.n_particles, inst.t, inst.C0, inst.D0)
    psi = ref.prep(inst.prep_layers)
    a, b = inst.current_pair
    w = ref.w_psi(psi, a, b, inst.tau, inst.steps, inst.q, inst.E0)
    if inst.component == "im":
        w = -1j * w
    c_ab = np.vdot(psi, ref.w_psi(psi, a, b, inst.tau, inst.steps, inst.q, inst.E0))
    return embed(((psi + w) / 2, (psi - w) / 2), inst), c_ab


def simulate(qc):
    """Final statevector with Aer (4 threads), including the circuit's global phase."""
    from qiskit_aer import AerSimulator
    sim = AerSimulator(method="statevector", max_parallel_threads=4, fusion_enable=True)
    c = qc.copy()
    c.save_statevector()
    sv = np.asarray(sim.run(c, shots=1).result().get_statevector(c))
    return sv


def z_anc(sv, inst):
    p = np.abs(sv) ** 2
    half = 1 << inst.anc
    return p[:half].sum() - p[half:].sum()


def propagate_basis(qc, bits):
    """Push computational basis states through a circuit of x, cx, ccx and diagonal gates.

    bits: bool array (N, n_qubits). Returns (output bits, phases). Raises on any other gate, so it also proves
    the block is a permutation times a diagonal.
    """
    bits = bits.copy()
    phase = np.full(len(bits), np.exp(1j * qc.global_phase), dtype=complex)
    for ci in qc.data:
        name = ci.operation.name
        q = [qc.find_bit(x).index for x in ci.qubits]
        if name == "x":
            bits[:, q[0]] ^= True
        elif name == "cx":
            bits[:, q[1]] ^= bits[:, q[0]]
        elif name == "ccx":
            bits[:, q[2]] ^= bits[:, q[0]] & bits[:, q[1]]
        elif name == "rz":
            th = float(ci.operation.params[0])
            phase *= np.where(bits[:, q[0]], np.exp(1j * th / 2), np.exp(-1j * th / 2))
        else:
            raise AssertionError(f"gate {name} is not classical-reversible or rz")
    return bits, phase


def all_bits(n, n_qubits, where):
    """All 2^n assignments of the qubits listed in `where`, other qubits 0."""
    vals = np.arange(2 ** n)
    bits = np.zeros((len(vals), n_qubits), dtype=bool)
    for k, q in enumerate(where):
        bits[:, q] = (vals >> k) & 1
    return bits


# --------------------------------------------------------------------------- #
# Instances
# --------------------------------------------------------------------------- #

BENCH = B.benchmark_instance(6)


def reduced(dim, n, layers=8, steps=6, e0=-7.0, **kw):
    """Same couplings, angles, q and step as the benchmark on a smaller lattice or with fewer nucleons."""
    return replace(BENCH, dim=dim, n_particles=n, prep_layers=BENCH.prep_layers[:layers], steps=steps,
                   tau=0.005 * steps, E0=e0, **kw)


# --------------------------------------------------------------------------- #
# Block tests at the full-instance size
# --------------------------------------------------------------------------- #

def test_gate_set_and_register():
    built = B.build_circuit(BENCH)
    names = {ci.operation.name for ci in built.circuit.data}
    assert names <= set(B.LOGICAL_GATES)
    assert built.circuit.num_qubits == 36
    assert (BENCH.n_pos, BENCH.n_work) == (24, 11)
    used = {built.circuit.find_bit(q).index for ci in built.circuit.data for q in ci.qubits}
    assert used == set(range(36))            # every one of the 36 qubits is touched


def test_benchmark_instance_is_the_chapters():
    """The benchmark's parameters against app01 ("Concrete instantiations"), typed here and not read from build.py.

    Without this a wrong tau in benchmark_instance changes no count and passes every other test once the QASM
    export is regenerated (found by mutation in the verification pass).
    """
    assert (BENCH.dim, BENCH.n_particles, BENCH.L) == (3, 4, 3)
    assert (BENCH.t, BENCH.C0, BENCH.D0) == (T, C0, D0)
    assert abs(BENCH.q - Q) < 1e-15
    assert len(BENCH.prep_layers) == 8
    assert BENCH.current_pair == (0, 0) and BENCH.component == "re" and not BENCH.merge_ancilla_phases
    assert abs(BENCH.E0 - (-52.67220)) < 5e-6
    for steps, tau in ((2, 0.01), (4, 0.02), (6, 0.03)):          # the three rows of the chapter's table
        inst = B.benchmark_instance(steps)
        assert inst.steps == steps and abs(inst.tau - tau) < 1e-15
        assert inst.prep_layers == BENCH.prep_layers and inst.E0 == BENCH.E0


def test_uniform_state_block():
    qc = B.block_circuit(2, lambda em: B.prepare3(em, 0, 1))
    col = Operator(qc).data[:, 0]
    assert np.max(np.abs(col - np.array([1, 1, 1, 0]) / math.sqrt(3))) < TOL
    inv = B.block_circuit(2, lambda em: (B.prepare3(em, 0, 1), B.prepare3(em, 0, 1, inverse=True)))
    assert np.max(np.abs(Operator(inv).data - np.eye(4))) < TOL


@pytest.mark.parametrize("dt", [0.003, 0.0173, -0.0125])
def test_kinetic_coordinate_block(dt):
    target = expm(-1j * dt * hop_matrix(T))
    qc = B.block_circuit(2, lambda em: B.kinetic_coord(em, 0, 1, 3 * T * dt))
    U = Operator(qc).data * np.exp(-1j * T * dt)             # the t I phase is added by kinetic_layer
    assert np.max(np.abs(U[:3, :3] - target)) < TOL
    assert abs(abs(U[3, 3]) - 1) < TOL                       # the unused pattern 11 does not mix in
    # controlled version on (lo, hi, work, control) = qubits 0..3
    qc = B.block_circuit(4, lambda em: (B.kinetic_coord(em, 0, 1, 3 * T * dt, control=3, work=2),
                                        em.phase(-T * dt, 3)))
    U = Operator(qc).data
    on = np.arange(3) + 8
    assert np.max(np.abs(U[np.ix_(on, on)] - target)) < TOL
    assert np.max(np.abs(U[:3, :3] - np.eye(3))) < TOL       # control off: identity, no phase
    assert np.max(np.abs(U[4:8, :4])) < TOL and np.max(np.abs(U[12:, 8:12])) < TOL   # workspace returns to 0


def test_kinetic_layer_two_coordinates_exact_phase():
    inst = replace(BENCH, dim=1, n_particles=2)              # 4 position bits, 1 workspace, 1 ancilla
    dt = 0.0111
    U1 = expm(-1j * dt * hop_matrix(T))
    target = np.kron(U1, U1)
    li = legal_index(2).ravel()                              # row-major over (x_0, x_1)
    qc = B.block_circuit(inst.n_qubits, lambda em: B.kinetic_layer(em, inst, dt))
    U = Operator(qc).data
    assert np.max(np.abs(U[np.ix_(li, li)] - target)) < TOL  # exact, including the global phase
    qc = B.block_circuit(inst.n_qubits, lambda em: B.kinetic_layer(em, inst, dt, control=inst.anc))
    U = Operator(qc).data
    on = li + (1 << inst.anc)
    assert np.max(np.abs(U[np.ix_(on, on)] - target)) < TOL
    assert np.max(np.abs(U[np.ix_(li, li)] - np.eye(9))) < TOL


def test_current_block():
    inst = replace(BENCH, dim=1, n_particles=2)
    for sign in (+1, -1):
        qc = B.block_circuit(inst.n_qubits, lambda em: B.current(em, inst, 1, sign, inst.anc))
        U = Operator(qc).data
        d = np.diag(U)
        assert np.max(np.abs(U - np.diag(d))) < TOL
        for x0 in range(3):
            for x1 in range(3):
                i = x0 + 4 * x1
                assert abs(d[i] - 1) < TOL                                           # control off
                assert abs(d[i + (1 << inst.anc)] - np.exp(1j * sign * Q * x1)) < TOL  # control on: particle 1 only


def _contact_case(k, dim, controlled, angle=0.7391):
    """Contact of k particles in `dim` dimensions on every computational basis input, workspace clean."""
    nb = 2 * dim
    n_pos = k * nb
    n_work = (k - 1) * nb - 1
    nq = n_pos + n_work + 1
    anc = nq - 1
    bit_lists = [list(range(f * nb, (f + 1) * nb)) for f in range(k)]
    work = list(range(n_pos, n_pos + n_work))
    qc = B.block_circuit(nq, lambda em: B.contact(em, bit_lists, angle, work, anc if controlled else None))
    inputs = all_bits(n_pos + 1, nq, list(range(n_pos)) + [anc])
    out, phase = propagate_basis(qc, inputs)
    assert np.array_equal(out, inputs)                        # positions unchanged, workspace back to 0
    vals = np.arange(len(inputs))
    pos = [(vals >> (f * nb)) & (2 ** nb - 1) for f in range(k)]
    same = np.ones(len(inputs), dtype=bool)
    for f in range(1, k):
        same &= pos[f] == pos[0]
    on = same & (inputs[:, anc] if controlled else True)
    assert np.max(np.abs(phase - np.where(on, np.exp(1j * angle), 1.0))) < TOL
    return qc


@pytest.mark.parametrize("controlled", [False, True])
def test_pair_contact_full_size(controlled):
    _contact_case(2, 3, controlled)           # 12 position bits, 5 workspace, 1 ancilla: 18 qubits, 2^13 inputs


@pytest.mark.parametrize("controlled", [False, True])
def test_triple_contact_full_size(controlled):
    qc = _contact_case(3, 3, controlled)      # 18 position bits, 11 workspace, 1 ancilla: 30 qubits, 2^19 inputs
    assert qc.num_qubits == 30 and CNT.walk(qc)["toffoli"] == 22


@pytest.mark.parametrize("controlled", [False, True])
def test_contact_basis_propagation_agrees_with_dense_unitary(controlled):
    """The basis-state propagator used above, checked against a dense unitary on a 6-qubit contact."""
    angle = 0.31
    qc = B.block_circuit(6, lambda em: B.contact(em, [[0, 1], [2, 3]], angle, [4], 5 if controlled else None))
    U = Operator(qc).data
    inputs = all_bits(6, 6, list(range(6)))
    out, phase = propagate_basis(qc, inputs)
    V = np.zeros((64, 64), dtype=complex)
    V[(out * (1 << np.arange(6))).sum(axis=1), np.arange(64)] = phase
    assert np.max(np.abs(U - V)) < TOL


@pytest.mark.parametrize("controlled", [False, True])
def test_potential_layer_all_pairs_and_triples(controlled):
    """exp(-i dt V) for 4 nucleons on 3^2 (16 position bits, 7 workspace): every legal basis state."""
    inst = replace(BENCH, dim=2)
    dt = 0.0041
    qc = B.block_circuit(inst.n_qubits, lambda em: B.potential_layer(em, inst, dt, inst.anc if controlled else None))
    ref = Ref(2, 4)
    li = legal_index(inst.n_coords).ravel()
    for anc in (0, 1):
        bits = np.zeros((len(li), inst.n_qubits), dtype=bool)
        for k in range(inst.n_pos):
            bits[:, k] = (li >> k) & 1
        bits[:, inst.anc] = bool(anc)
        out, phase = propagate_basis(qc, bits)
        assert np.array_equal(out, bits)
        want = np.exp(-1j * dt * ref.V.ravel()) if (anc or not controlled) else np.ones(len(li))
        assert np.max(np.abs(phase - want)) < TOL


# --------------------------------------------------------------------------- #
# Reference self-checks
# --------------------------------------------------------------------------- #

def test_reference_product_formula_converges_to_exact_evolution():
    """Ref against exp(-i H tau) from a matrix of explicit elements, 4 nucleons on 3^1: second order."""
    H, states = dense_hamiltonian(1, 4)
    ref = Ref(1, 4)
    v = np.random.default_rng(7).normal(size=81) + 0j
    v /= np.linalg.norm(v)
    assert np.max(np.abs(ref.h_apply(v.reshape(ref.shape)).ravel() - H @ v)) < 1e-10
    tau = 0.03
    exact = expm(-1j * tau * (H - ref.hop_constant * np.eye(81))) @ v
    err = [np.linalg.norm(ref.trotter(v.reshape(ref.shape), tau, n).ravel() - exact) for n in (6, 12, 24)]
    assert 3.6 < err[0] / err[1] < 4.4 and 3.6 < err[1] / err[2] < 4.4
    assert err[2] < 2e-3


def test_reference_reproduces_chapter_numbers_on_full_3x3x3():
    """Full 27^4 coordinate space in numpy: E0, prepared energy, overlap and the 6-step correlator."""
    from scipy.sparse.linalg import LinearOperator, eigsh
    ref = Ref(3, 4)
    n = 3 ** 12
    op = LinearOperator((n, n), matvec=lambda v: ref.h_apply(v.reshape(ref.shape)).ravel(), dtype=float)
    e, vec = eigsh(op, k=1, which="SA", tol=1e-10, v0=np.ones(n))
    e0, gs = e[0], vec[:, 0]
    study = json.loads((B.STUDY / "larger_box_results.json").read_text())
    prep = json.loads((B.STUDY / "unitary_prep_results.json").read_text())[-1]
    assert abs(e0 - study["ground_energy_MeV"]) < 1e-6           # -52.67220 (app01)
    assert abs(e0 - BENCH.E0) < 1e-6
    psi = ref.prep(BENCH.prep_layers)
    energy = np.vdot(psi, ref.h_apply(psi)).real
    assert abs(energy - prep["energy_MeV"]) < 1e-7               # -52.45933 (app01)
    assert abs(energy - (-52.45933)) < 5e-6
    overlap = abs(np.vdot(gs, psi.ravel())) ** 2
    assert abs(overlap - prep["fidelity"]) < 1e-7                # 0.997943 (app01)
    assert abs(overlap - 0.997943) < 5e-7
    resp = json.loads((B.STUDY / "larger_unitary_response.json").read_text())["rows"]
    for steps, tau in ((2, 0.01), (4, 0.02), (6, 0.03)):
        c = sum(np.vdot(psi, ref.w_psi(psi, a, b, tau, steps, Q, e0)) for a in (0, 1) for b in (0, 1)) / 4
        row = next(r for r in resp if abs(r["tau_MeV_inv"] - tau) < 1e-12 and r["steps"] == steps)
        assert abs(c - complex(*row["circuit"])) < 1e-7          # product-formula value of the study
        bias = abs(c - complex(*row["exact"]))                   # 'exact' is the study's, not recomputed here
        assert abs(bias - {2: 0.00609, 4: 0.01098, 6: 0.01434}[steps]) < 5e-6    # app01 table


# --------------------------------------------------------------------------- #
# Complete circuits on reduced instances
# --------------------------------------------------------------------------- #

def check_complete(inst):
    built = B.build_circuit(inst, measure=False)
    sv = simulate(built.circuit)
    want, c_ab = expected_final(inst)
    err = np.max(np.abs(sv - want))
    z = z_anc(sv, inst)
    return err, z, c_ab


@pytest.mark.parametrize("pair", [(0, 0), (0, 1), (1, 0), (1, 1)])
@pytest.mark.parametrize("component", ["re", "im"])
def test_complete_circuit_four_nucleons_1d(pair, component):
    """12 qubits, all eight preparation layers and six steps with the benchmark angles."""
    ref = Ref(1, 4)
    e0 = np.linalg.eigvalsh(dense_hamiltonian(1, 4)[0])[0]
    inst = reduced(1, 4, e0=e0, current_pair=pair, component=component)
    assert inst.n_qubits == 12
    err, z, c_ab = check_complete(inst)
    assert err < TOL
    assert abs(z - (c_ab.real if component == "re" else c_ab.imag)) < TOL
    assert abs(c_ab) > 0.05                                   # not a trivially zero comparison


def test_complete_circuit_two_nucleons_3d():
    """18 qubits: the full 3^3 coordinates and six-bit comparisons, eight layers, six steps."""
    inst = reduced(3, 2, current_pair=(1, 0), component="im")
    assert inst.n_qubits == 18
    err, z, c_ab = check_complete(inst)
    assert err < TOL and abs(z - c_ab.imag) < TOL


def test_complete_circuit_three_nucleons_2d():
    """20 qubits: pair and triple contacts across two axes, eight layers, six steps."""
    inst = reduced(2, 3, current_pair=(0, 1), component="re")
    assert inst.n_qubits == 20
    err, z, c_ab = check_complete(inst)
    assert err < TOL and abs(z - c_ab.real) < TOL


def test_complete_circuit_four_nucleons_2d_24_qubits():
    """24 qubits, the largest reduced instance: 4 nucleons on 3^2, eight layers, six steps."""
    inst = reduced(2, 4, current_pair=(0, 1), component="im")
    assert inst.n_qubits == 24
    err, z, c_ab = check_complete(inst)
    assert err < TOL and abs(z - c_ab.imag) < TOL


def test_static_strength_circuit_steps_zero():
    inst = reduced(1, 4, steps=0, current_pair=(0, 1))
    err, z, c_ab = check_complete(inst)
    assert err < TOL and abs(z - c_ab.real) < TOL


def test_merged_ancilla_phase_variant_is_the_same_unitary():
    for inst in (reduced(1, 4, current_pair=(1, 0), component="im", merge_ancilla_phases=True),
                 reduced(2, 3, layers=2, steps=2, merge_ancilla_phases=True)):
        err, z, c_ab = check_complete(inst)
        assert err < TOL
    base = CNT.walk(B.build_circuit(BENCH).circuit)
    merged = CNT.walk(B.build_circuit(replace(BENCH, merge_ancilla_phases=True)).circuit)
    assert merged["rotations"] < base["rotations"]
    assert (merged["toffoli"], merged["fixed_T"]) == (base["toffoli"], base["fixed_T"])


# --------------------------------------------------------------------------- #
# Tests that must fail when something is wrong
# --------------------------------------------------------------------------- #

def _mutated(built, index, new_ops):
    qc = built.circuit.copy_empty_like()
    qc.global_phase = built.circuit.global_phase
    for i, ci in enumerate(built.circuit.data):
        if i == index:
            for op, qubits in new_ops:
                qc.append(op, qubits)
        else:
            qc.append(ci.operation, ci.qubits, ci.clbits)
    return qc


def test_wrong_reference_is_detected():
    """Dropped triple term, flipped pair sign, flipped q, wrong step count: each must show up."""
    inst = reduced(1, 4, e0=-30.0)
    sv = simulate(B.build_circuit(inst, measure=False).circuit)
    assert np.max(np.abs(sv - expected_final(inst)[0])) < TOL
    for wrong in (replace(inst, D0=0.0), replace(inst, C0=-inst.C0), replace(inst, q=-inst.q),
                  replace(inst, t=-inst.t), replace(inst, steps=5), replace(inst, current_pair=(0, 1)),
                  replace(inst, E0=inst.E0 + 5.0)):
        assert np.linalg.norm(sv - expected_final(wrong)[0]) > NEG


def test_mutated_circuit_is_detected():
    """Remove one Toffoli, flip one rotation, drop one T: the full-state comparison must fail each time."""
    from qiskit.circuit.library import RZGate
    inst = reduced(1, 4, e0=-30.0)
    built = B.build_circuit(inst, measure=False)
    want = expected_final(inst)[0]
    data = built.circuit.data
    start = next(s for label, s, e in built.blocks if label == "evolve/3/potential")
    i_ccx = next(i for i in range(start, len(data)) if data[i].operation.name == "ccx")
    i_rz = next(i for i in range(start, len(data)) if data[i].operation.name == "rz")
    kstart = next(s for label, s, e in built.blocks if label == "prep/5/kinetic")
    i_t = next(i for i in range(kstart, len(data)) if data[i].operation.name == "t")
    flipped = RZGate(-float(data[i_rz].operation.params[0]))
    for index, ops in ((i_ccx, []), (i_rz, [(flipped, data[i_rz].qubits)]), (i_t, [])):
        sv = simulate(_mutated(built, index, ops))
        assert np.linalg.norm(sv - want) > NEG
    # and the unmutated copy passes through the same path
    assert np.max(np.abs(simulate(_mutated(built, -1, [])) - want)) < TOL


def test_wrong_contact_sign_is_detected_at_full_size():
    """The 30-qubit triple contact with the opposite angle fails the basis-state check."""
    nb, k = 6, 3
    bit_lists = [list(range(f * nb, (f + 1) * nb)) for f in range(k)]
    work = list(range(18, 29))
    qc = B.block_circuit(30, lambda em: B.contact(em, bit_lists, -0.7391, work, 29))
    inputs = all_bits(19, 30, list(range(18)) + [29])
    out, phase = propagate_basis(qc, inputs)
    vals = np.arange(len(inputs))
    pos = [(vals >> (f * nb)) & 63 for f in range(k)]
    on = (pos[0] == pos[1]) & (pos[0] == pos[2]) & inputs[:, 29]
    assert np.max(np.abs(phase - np.where(on, np.exp(1j * 0.7391), 1.0))) > NEG


# --------------------------------------------------------------------------- #
# Counting and export
# --------------------------------------------------------------------------- #

def test_classification_of_rotation_angles():
    assert CNT.classify("rz", [math.pi / 4]) == "fixed_T"
    assert CNT.classify("rz", [-3 * math.pi / 4]) == "fixed_T"
    assert CNT.classify("rz", [math.pi / 2]) == "clifford"
    assert CNT.classify("rz", [0.3]) == "rotation"
    assert CNT.classify("ccx", []) == "toffoli"
    with pytest.raises(ValueError):
        CNT.classify("cp", [0.3])
    built = B.build_circuit(BENCH)                 # no rz in the benchmark is accidentally exact
    n_rz = sum(ci.operation.name == "rz" for ci in built.circuit.data)
    assert n_rz == CNT.walk(built.circuit)["rotations"]


def test_blocks_partition_and_counts_json_is_current():
    built = B.build_circuit(BENCH)
    blocks, fam, total = CNT.per_block(built)
    stored = json.loads((HERE / "counts.json").read_text())
    assert stored["benchmark"]["counts"] == total
    t, c_t, eps = CNT.t_total(total)
    assert abs(stored["benchmark"]["t_total_report"] - t) < 1e-6


def test_counts_against_model_line_items():
    """Not a tuning target: this records that the walked counts and the model literals agree, per step count."""
    for steps in (2, 4, 6):
        total = CNT.walk(B.build_circuit(B.benchmark_instance(steps)).circuit)
        m = CNT.model_items(steps)
        assert (total["rotations"], total["toffoli"], total["fixed_T"]) == (m["rotations"], m["toffoli"], m["fixed_T"])
        assert abs(CNT.t_total(total)[0] / m["t_total"] - 1) < 1e-12


def test_depth_walker_on_known_circuits():
    qc = QuantumCircuit(3)
    qc.t(0); qc.t(1); qc.cx(0, 1); qc.t(1); qc.rz(0.3, 2)
    assert CNT.depth(qc, 1.0) == 2.0
    assert CNT.depth(qc, 5.0) == 5.0
    assert CNT.depth_of_one_toffoli() in (3.0, 4.0)


def test_qasm3_export_reloads_to_the_same_gates():
    from qiskit import qasm3
    path = HERE / "ch02_3x3x3_steps6_re_p0p0.qasm3"
    built = B.build_circuit(BENCH)
    loaded = qasm3.loads(path.read_text())
    assert loaded.num_qubits == 36
    a = [(ci.operation.name, [built.circuit.find_bit(q).index for q in ci.qubits],
          [float(p) for p in ci.operation.params]) for ci in built.circuit.data]
    b = [(ci.operation.name, [loaded.find_bit(q).index for q in ci.qubits],
          [float(p) for p in ci.operation.params]) for ci in loaded.data]
    assert len(a) == len(b)
    for x, y in zip(a, b):
        assert x[0] == y[0] and x[1] == y[1]
        assert all(abs(p - r) < 1e-12 for p, r in zip(x[2], y[2]))
    assert CNT.walk(loaded) == CNT.walk(built.circuit)


# --------------------------------------------------------------------------- #
# python test_circuit.py: write the observed errors of the complete-circuit checks to verification.json
# --------------------------------------------------------------------------- #

def error_report():
    import time
    e0 = np.linalg.eigvalsh(dense_hamiltonian(1, 4)[0])[0]
    cases = [("4 nucleons, 3^1", reduced(1, 4, e0=e0)),
             ("2 nucleons, 3^3", reduced(3, 2, current_pair=(1, 0), component="im")),
             ("3 nucleons, 3^2", reduced(2, 3, current_pair=(0, 1))),
             ("4 nucleons, 3^2", reduced(2, 4, current_pair=(0, 1), component="im"))]
    rows = []
    for name, inst in cases:
        t0 = time.time()
        built = B.build_circuit(inst, measure=False)
        sv = simulate(built.circuit)
        want, c_ab = expected_final(inst)
        z = z_anc(sv, inst)
        target = c_ab.real if inst.component == "re" else c_ab.imag
        c = CNT.walk(built.circuit)
        rows.append(dict(instance=name, qubits=inst.n_qubits, prep_layers=len(inst.prep_layers), steps=inst.steps,
                         gates=len(built.circuit.data), rotations=c["rotations"], toffoli=c["toffoli"],
                         fixed_T=c["fixed_T"], current_pair=list(inst.current_pair), component=inst.component,
                         max_amplitude_error=float(np.max(np.abs(sv - want))),
                         ancilla_expectation=float(z), reference_value=float(target),
                         expectation_error=float(abs(z - target)), seconds=round(time.time() - t0, 1)))
        print(rows[-1])
    out = dict(tolerance=TOL, simulator="qiskit-aer statevector, double precision, 4 threads, one machine",
               reference="numpy product of exponentials built from the Hamiltonian definition (class Ref)",
               complete_circuits=rows)
    (HERE / "verification.json").write_text(json.dumps(out, indent=1) + "\n")


if __name__ == "__main__":
    error_report()
