"""Checks of the Ch. 7 2028 step circuit against references that do not use the circuit code.

Run:  OMP_NUM_THREADS=4 python -m pytest circuits/ch07_scalar_ds_step/test_circuit.py -q

References
  per-site blocks    dense matrices from the operator definitions: x = grid values, F = centered DFT formula,
                     p = F^dag x F, exp(-i u p^2/2) by scipy.linalg.expm
  ring step          the model's own definition of the map (estimates.ch07_curved_space._ds2028_setup,
                     _ds2028_blocks, and the frame algebra of ds2028_grid_check, restated in RefMap below)
  d-dimensional step RefMap: real-space covariance and Laplacian built here (not build.step_parameters,
                     which works per momentum mode), applied to a state vector with numpy
TOL = 1e-9. Unitaries: largest entry-wise deviation after removing one global phase (dist). State vectors:
2-norm of the difference of the normalized states after removing one global phase (sdist). The observed
deviations are 1e-14 or smaller.
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np
import pytest
from scipy.linalg import expm

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
for _p in (str(REPO), str(HERE)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from qiskit import QuantumCircuit                      # noqa: E402
from qiskit.quantum_info import Operator               # noqa: E402

import build as B                                      # noqa: E402
import count as C                                      # noqa: E402
from estimates import ch07_curved_space as model       # noqa: E402

TOL = 1e-9
PI = math.pi


# --------------------------------------------------------------------------- #
# Reference operators, from the definitions
# --------------------------------------------------------------------------- #

def grid(n_q):
    K = 2 ** n_q
    dx = math.sqrt(2 * PI / K)
    jj = np.arange(K) - (K - 1) / 2
    return K, dx, jj * dx, jj


def F_matrix(n_q):
    K, _, _, jj = grid(n_q)
    return np.exp(-1j * np.outer(jj, jj) * 2 * PI / K) / math.sqrt(K)


def P_matrix(n_q):
    _, _, X, _ = grid(n_q)
    F = F_matrix(n_q)
    return F.conj().T @ np.diag(X) @ F


def kinetic_site(n_q, u):
    P = P_matrix(n_q)
    return expm(-0.5j * u * (P @ P))


def dist(U, V):
    """Largest entry-wise deviation after aligning one global phase."""
    U, V = np.asarray(U), np.asarray(V)
    ov = np.vdot(V, U)
    assert abs(ov) > 1e-12, "no overlap at all"
    return float(np.max(np.abs(U - (ov / abs(ov)) * V)))


def sdist(a, b):
    """2-norm of the difference of two normalized state vectors after aligning one global phase."""
    a, b = np.asarray(a), np.asarray(b)
    ov = np.vdot(b, a)
    assert abs(ov) > 1e-12, "no overlap at all"
    return float(np.linalg.norm(a - (ov / abs(ov)) * b))


def emit(n_qubits, fn):
    em = B.Emitter(n_qubits)
    fn(em)
    return em.qc


def unitary(qc):
    return Operator(qc).data


class RefMap:
    """The step as a list of operations on a state vector, from the Hamiltonian and the frame rule.

    Lab frame: K(v): phi += v pi.  D: pi -= W phi, W = h (a^3 m^2 + a Lap).  Frame on every site
    z' = S z, S refit (dilation and shear) to the site covariance at each kinetic-block boundary.
    ops: ("x", Q) -> exp(-i x^T Q x / 2),  ("p", u) -> exp(-i u sum_x p_x^2 / 2), in the order they act.
    """

    def __init__(self, dims, m, h=0.5, n_steps=1, C0=None, lap=None):
        self.dims = tuple(dims)
        V = self.V = int(np.prod(dims))
        b = math.sqrt(12.0 / (4.0 - m * m))
        if lap is None:
            lap = np.zeros((V, V))
            for ax, L in enumerate(dims):                         # kron of ring Laplacians
                ring = 2 * np.eye(L) - np.roll(np.eye(L), 1, 0) - np.roll(np.eye(L), -1, 0)
                mats = [np.eye(l) for l in dims]
                mats[ax] = ring
                full = mats[0]
                for mm in mats[1:]:
                    full = np.kron(full, mm)
                lap += full / b ** 2
        self.lap = lap
        if C0 is None:
            C0 = self.in_state(m, b)
        self.C0 = C0
        I, Z = np.eye(V), np.zeros((V, V))
        seq = []
        for s in range(n_steps):
            am = math.exp(s * h + h / 2)
            for blk in (("K", h / 2 / am ** 3), ("D", h * (am ** 3 * m * m * I + am * lap)), ("K", h / 2 / am ** 3)):
                if seq and blk[0] == "K" and seq[-1][0] == "K":
                    seq[-1] = ("K", seq[-1][1] + blk[1])
                else:
                    seq.append(blk)
        self.seq = seq
        covs, Mtot = [C0], np.eye(2 * V)
        for typ, v in seq:
            M = np.block([[I, v * I], [Z, I]]) if typ == "K" else np.block([[I, Z], [-v, I]])
            covs.append(M @ covs[-1] @ M.T)
            Mtot = M @ Mtot
        self.covs, self.M_lab = covs, Mtot

        def bal(Cm):
            sxx, sxp, spp = Cm[0, 0], Cm[0, V], Cm[V, V]
            sc = math.sqrt(math.sqrt(sxx * spp - sxp * sxp) / sxx)
            return np.array([[sc, 0.0], [-sxp / (sc * sxx), 1 / sc]])
        frames = [bal(C0)]
        for i, (typ, _) in enumerate(seq):
            frames.append(frames[-1] if typ == "D" else bal(covs[i + 1]))
        self.frames = frames
        ops = []
        for i, (typ, v) in enumerate(seq):
            Sin, Sout = frames[i], frames[i + 1]
            if typ == "D":
                ops.append(("x", v / Sin[0, 0] ** 2))
                continue
            A = Sout @ np.array([[1.0, v], [0.0, 1.0]]) @ np.linalg.inv(Sin)
            up = A[0, 1]
            ops += [("x", (1 - A[0, 0]) / up * I), ("p", up), ("x", (1 - A[1, 1]) / up * I)]
        self.ops = ops

    def in_state(self, m, b):
        """Real-space in-state covariance, Bunch-Davies per lattice momentum (k = 0: a = 1 ground state)."""
        from scipy.special import hankel1, h1vp
        V, dims = self.V, self.dims
        coords = np.stack(np.unravel_index(np.arange(V), dims), 1)
        nu = math.sqrt(9 / 4 - m * m)
        ph = np.exp(1j * PI * (nu + 0.5) / 2)
        Cxx, Cxp, Cpp = (np.zeros((V, V)) for _ in range(3))
        for n in coords:
            kt2 = sum((2 / b * math.sin(PI * n[a] / dims[a])) ** 2 for a in range(len(dims)))
            if kt2 < 1e-14:
                sxx, spp, sxp = 1 / (2 * m), m / 2, 0.0
            else:
                z = math.sqrt(kt2)
                q = math.sqrt(PI) / 2 * hankel1(nu, z) * ph
                p = math.sqrt(PI) / 2 * ph * (-1.5 * hankel1(nu, z) - z * h1vp(nu, z))
                sxx, spp, sxp = abs(q) ** 2, abs(p) ** 2, float(np.real(q * np.conj(p)))
            phase = 2 * PI * (coords @ (n / np.array(dims)))
            wave = np.cos(phase[:, None] - phase[None, :]) / V
            Cxx += sxx * wave
            Cxp += sxp * wave
            Cpp += spp * wave
        return np.block([[Cxx, Cxp], [Cxp.T, Cpp]])

    # ---- state-vector action; tensor axis x <-> site x, converted to qiskit order at the ends -------- #

    def apply(self, psi_qiskit, n_q, ops=None):
        K, _, X, _ = grid(n_q)
        V = self.V
        psi = psi_qiskit.reshape((K,) * V).transpose(*reversed(range(V)))     # axis x = site x
        for typ, v in (self.ops if ops is None else ops):
            if typ == "p":
                U = kinetic_site(n_q, v)
                for ax in range(V):
                    psi = np.moveaxis(np.tensordot(U, psi, axes=([1], [ax])), 0, ax)
            else:
                expo = np.zeros((K,) * V)
                for x in range(V):
                    shx = [1] * V
                    shx[x] = K
                    Xx = X.reshape(shx)
                    if v[x, x] != 0:
                        expo = expo + 0.5 * v[x, x] * Xx * Xx
                    for y in range(x + 1, V):
                        if v[x, y] != 0:
                            shy = [1] * V
                            shy[y] = K
                            expo = expo + v[x, y] * Xx * X.reshape(shy)
                psi = psi * np.exp(-1j * expo)
        return psi.transpose(*reversed(range(V))).reshape(-1)

    def dense(self, n_q, ops=None):
        dim = (2 ** n_q) ** self.V
        cols = [self.apply(np.eye(dim, dtype=complex)[:, j], n_q, ops) for j in range(dim)]
        return np.stack(cols, 1)


def run(qc, psi0):
    """Final state vector of qc from psi0 (qiskit little-endian), on Aer with at most 4 threads."""
    from qiskit_aer import AerSimulator
    full = QuantumCircuit(qc.num_qubits)
    full.set_statevector(psi0)
    full.compose(qc, inplace=True)
    full.save_statevector()
    sim = AerSimulator(method="statevector", max_parallel_threads=4, precision="double")
    return np.asarray(sim.run(full, shots=1).result().get_statevector(full))


def random_state(n_qubits, seed):
    rng = np.random.default_rng(seed)
    v = rng.normal(size=2 ** n_qubits) + 1j * rng.normal(size=2 ** n_qubits)
    return v / np.linalg.norm(v)


def model_ring_ref(m, L, n_steps=1):
    """RefMap fed with the model's own ring covariance and Laplacian (estimates._ds2028_setup)."""
    b, E, C0, lap = model._ds2028_setup(m, L)
    return RefMap((L,), m, n_steps=n_steps, C0=C0, lap=lap), E


# --------------------------------------------------------------------------- #
# 1. Per-site blocks against dense matrices
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("n_q", [2, 3, 4])
def test_site_x2_phase_block(n_q):
    K, dx, X, _ = grid(n_q)
    for gamma in (0.37, -1.9):
        qc = emit(n_q, lambda em: B.site_phase(em, list(range(n_q)), gamma, dx))
        assert dist(unitary(qc), np.diag(np.exp(-0.5j * gamma * X ** 2))) < TOL


@pytest.mark.parametrize("n_q", [2, 3, 4])
def test_centered_dft_block(n_q):
    """Dc then the swap-free inverse QFT equals (bit reversal) . F, F the centered DFT of the definition."""
    K = 2 ** n_q
    reg = list(range(n_q))
    qc = emit(n_q, lambda em: (B.centering(em, reg), B.qft_swapless(em, reg, inverse=True), B.centering(em, reg)))
    U = unitary(qc)
    rev = [int(format(j, f"0{n_q}b")[::-1], 2) for j in range(K)]
    R = np.eye(K)[rev]
    # the second centering layer sits on bit-reversed qubits here, so move it through R
    D = np.diag(np.exp(2j * PI * (K - 1) / 2 * np.arange(K) / K))
    U_fixed = D @ R @ D.conj() @ U
    assert dist(U_fixed, F_matrix(n_q)) < TOL


@pytest.mark.parametrize("n_q", [2, 3, 4])
def test_kinetic_block_equals_exp_p2(n_q):
    _, dx, _, _ = grid(n_q)
    reg = list(range(n_q))
    for u in (0.118, 0.9, -0.4):
        qc = emit(n_q, lambda em: B.kinetic_block(em, reg, u, dx))
        assert dist(unitary(qc), kinetic_site(n_q, u)) < TOL


@pytest.mark.parametrize("n_q", [2, 3])
def test_link_block(n_q):
    K, dx, X, _ = grid(n_q)
    kappa = -0.147
    qc = emit(2 * n_q, lambda em: B.link_phase(em, list(range(n_q)), list(range(n_q, 2 * n_q)), kappa, dx))
    # qiskit order: site 1 is the more significant factor
    ref = np.diag(np.exp(-1j * kappa * np.kron(X, X)))
    assert dist(unitary(qc), ref) < TOL


def test_centering_layers_cancel_between_kinetic_blocks():
    """merge_centering = True (two layers per shot) is the same unitary as a pair on every kinetic block."""
    inst = B.Instance(dims=(3,), n_q=3)
    a = unitary(B.build_step(inst, merge_centering=True).qc)
    b = unitary(B.build_step(inst, merge_centering=False).qc)
    assert dist(a, b) < TOL


# --------------------------------------------------------------------------- #
# 2. Coefficients: build.step_parameters against the Hamiltonian
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("dims,n_steps", [((3,), 1), ((4,), 1), ((4,), 2), ((2, 2, 2), 1), ((3, 3), 1), ((4, 4, 4), 1)])
def test_coefficients_reproduce_the_lab_frame_step(dims, n_steps):
    """The circuit's sequence of quadratic phases, read as a linear map on (x, p), equals
    S_out (K D K) S_in^-1 with K, D straight from the Hamiltonian (RefMap.M_lab)."""
    inst = B.Instance(dims=dims, n_steps=n_steps, m=1.2)
    p = B.step_parameters(inst)
    ref = RefMap(dims, inst.m, inst.h, n_steps)
    V = inst.V
    I, Z = np.eye(V), np.zeros((V, V))
    M = np.eye(2 * V)
    J = len(p.kinetic_u)
    for j in range(J + 1):
        Q = p.site_gamma[j] * I
        if 0 < j < J:
            for (x, y), k in p.link_kappa[j - 1].items():
                Q[x, y] += k
                Q[y, x] += k
        M = np.block([[I, Z], [-Q, I]]) @ M
        if j < J:
            M = np.block([[I, p.kinetic_u[j] * I], [Z, I]]) @ M
    big = lambda S: np.kron(S, I)
    target = big(p.frames[-1]) @ ref.M_lab @ np.linalg.inv(big(p.frames[0]))
    assert np.max(np.abs(M - target)) < 1e-10
    for f_mine, f_ref in zip(p.frames, [ref.frames[0]] + [ref.frames[i + 1] for i, s in enumerate(ref.seq) if s[0] == "K"]):
        assert np.max(np.abs(f_mine - f_ref)) < 1e-10


def test_reference_in_state_matches_the_model_ring():
    """The reference built here and the model's _ds2028_setup give the same ring covariance and Laplacian."""
    for m in (1.0, 1.4):
        b, E, C0, lap = model._ds2028_setup(m, 4)
        ref = RefMap((4,), m)
        assert np.max(np.abs(ref.C0 - C0)) < 1e-12
        assert np.max(np.abs(ref.lap - lap)) < 1e-12
        assert abs(B.Instance(dims=(4,), m=m).spacing - b) < 1e-14


# --------------------------------------------------------------------------- #
# 3. Full step: dense unitary on a 3-site ring, state vectors on larger reduced lattices
# --------------------------------------------------------------------------- #

@pytest.fixture(scope="module")
def ring3():
    inst = B.Instance(dims=(3,), n_q=3, m=1.0)
    ref, _ = model_ring_ref(inst.m, 3)
    return inst, ref, ref.dense(inst.n_q)


def test_full_step_unitary_ring3(ring3):
    """9 qubits: every entry of the 512 x 512 step unitary against the model-defined map."""
    inst, ref, U_ref = ring3
    U = unitary(B.build_step(inst).qc)
    assert np.max(np.abs(U_ref.conj().T @ U_ref - np.eye(512))) < 1e-10
    assert dist(U, U_ref) < TOL


def test_two_merged_steps_ring3():
    inst = B.Instance(dims=(3,), n_q=2, m=1.4, n_steps=2)
    ref, _ = model_ring_ref(inst.m, 3, n_steps=2)
    assert dist(unitary(B.build_step(inst).qc), ref.dense(inst.n_q)) < TOL


def test_full_step_statevector_ring4():
    """12 qubits, K = 8, random input state."""
    inst = B.Instance(dims=(4,), n_q=3, m=1.0)
    ref, _ = model_ring_ref(inst.m, 4)
    psi0 = random_state(12, 1)
    assert sdist(run(B.build_step(inst).qc, psi0), ref.apply(psi0, 3)) < TOL


def test_full_step_statevector_torus_3x3_K4():
    """18 qubits: 2D torus with distinct neighbours in both directions, K = 4."""
    inst = B.Instance(dims=(3, 3), n_q=2, m=1.0)
    ref = RefMap(inst.dims, inst.m)
    psi0 = random_state(18, 2)
    assert sdist(run(B.build_step(inst).qc, psi0), ref.apply(psi0, 2)) < TOL


def test_full_step_statevector_2x2x2_K8():
    """24 qubits: the largest reduced 3D instance at the benchmark digitization (K = 8). On sides of length 2
    the two wraps of a direction coincide, so each site has 3 distinct neighbours with doubled coupling."""
    inst = B.Instance(dims=(2, 2, 2), n_q=3, m=1.0)
    ref = RefMap(inst.dims, inst.m)
    bt = B.build_step(inst)
    assert bt.qc.num_qubits == 24
    psi0 = random_state(24, 3)
    assert sdist(run(bt.qc, psi0), ref.apply(psi0, 3)) < TOL


# --------------------------------------------------------------------------- #
# 4. End to end against the model's own exact check
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("m", [1.0, 1.4])
def test_circuit_reproduces_model_grid_check(m):
    """4-site ring, K = 8, Gaussian in-state on the grid. The shell errors of the circuit's output state
    against the continuum map equal the ones estimates.ds2028_grid_check reports (its own numpy evolution),
    and stay at the chapter's 0.14% level."""
    L, n_q = 4, 3
    K, dx, X, _ = grid(n_q)
    inst = B.Instance(dims=(L,), n_q=n_q, m=m)
    bt = B.build_step(inst)
    ref, E = model_ring_ref(m, L)
    S0, S_out = ref.frames[0], ref.frames[-1]
    big = lambda S: np.kron(S, np.eye(L))
    Cr = big(S0) @ ref.C0 @ big(S0).T
    M0 = ((Cr[L:, :L] + Cr[:L, L:].T) / 2 + 0.5j * np.eye(L)) @ np.linalg.inv(Cr[:L, :L])
    xs = np.stack(np.meshgrid(*([X] * L), indexing="ij"), 0).reshape(L, -1)      # axis x = site x
    psi = np.exp(0.5j * np.einsum("in,ij,jn->n", xs, M0, xs))
    psi /= np.linalg.norm(psi)
    to_q = lambda v: v.reshape((K,) * L).transpose(*reversed(range(L))).reshape(-1)
    out = to_q(run(bt.qc, to_q(psi)))                                             # back to axis x = site x

    F, Pm = F_matrix(n_q), P_matrix(n_q)

    def site(ps, U, only=None):
        ps = ps.reshape((K,) * L)
        for ax in (range(L) if only is None else (only,)):
            ps = np.moveaxis(np.tensordot(U, ps, axes=([1], [ax])), 0, ax)
        return ps.reshape(-1)
    pr, pq = abs(out) ** 2, abs(site(out, F)) ** 2
    Cxp = np.array([[np.real(np.vdot(out, xs[i] * site(out, Pm, j))) for j in range(L)] for i in range(L)])
    cov = np.block([[(xs * pr) @ xs.T, Cxp], [Cxp.T, (xs * pq) @ xs.T]])
    Si = big(np.linalg.inv(S_out))
    err = model._shell_err(Si @ cov @ Si.T, ref.covs[-1], E, L)
    err_model = model.ds2028_grid_check(m, K, "adapted", n_steps=1, L=L, h=0.5, merged=True)[-1]
    assert np.max(np.abs(np.array(err) - np.array(err_model))) < 1e-8
    assert max(err) < 2e-3


# --------------------------------------------------------------------------- #
# 5. Tests that fail when a block is wrong
# --------------------------------------------------------------------------- #

def test_detects_wrong_sign_of_the_gradient_coupling(ring3):
    inst, ref, U_ref = ring3
    p = B.step_parameters(inst)
    p.link_kappa = [{k: -v for k, v in d.items()} for d in p.link_kappa]
    assert dist(unitary(B.build_step(inst, params=p).qc), U_ref) > 1e-2


def test_detects_a_dropped_rotation(ring3):
    inst, ref, U_ref = ring3
    bt = B.build_step(inst)
    for block in ("d1_link", "d1_site", "k1_pi2", "k2_qft", "shear_in", "centering_out"):
        idx = next(i for i, (ins, lab) in enumerate(zip(bt.qc.data, bt.labels))
                   if lab == block and ins.operation.name == "rz")
        qc = QuantumCircuit(bt.qc.num_qubits)
        for i, ins in enumerate(bt.qc.data):
            if i != idx:
                qc.append(ins.operation, ins.qubits)
        assert dist(unitary(qc), U_ref) > 1e-3, block


def test_detects_missing_bit_reversal_and_wrong_qft_direction():
    n_q = 3
    _, dx, _, _ = grid(n_q)
    reg = list(range(n_q))
    ref = kinetic_site(n_q, 0.3)

    def no_reversal(em):
        B.centering(em, reg); B.qft_swapless(em, reg, inverse=True)
        B.site_phase(em, reg, 0.3, dx, reversed_bits=False)
        B.qft_swapless(em, reg); B.centering(em, reg, inverse=True)

    def wrong_direction(em):
        B.centering(em, reg); B.qft_swapless(em, reg, inverse=False)
        B.site_phase(em, reg, 0.3, dx, reversed_bits=True)
        B.qft_swapless(em, reg, inverse=True); B.centering(em, reg, inverse=True)

    def no_centering(em):
        B.kinetic_block(em, reg, 0.3, dx, with_centering=False)
    for bad in (no_reversal, wrong_direction, no_centering):
        assert dist(unitary(emit(n_q, bad)), ref) > 1e-2, bad.__name__


def test_outer_shear_layers_are_not_identity_but_leave_field_readout_alone(ring3):
    """The frame shear at the two ends is a real part of the map (dropping it changes the unitary), and the
    last one is diagonal, so it cannot change field-basis probabilities."""
    inst, ref, U_ref = ring3
    U_no_out = unitary(B.build_step(inst, out_shear=False).qc)
    U_no_in = unitary(B.build_step(inst, in_shear=False).qc)
    assert dist(U_no_out, U_ref) > 1e-2 and dist(U_no_in, U_ref) > 1e-2
    psi0 = random_state(9, 5)
    assert np.max(np.abs(abs(U_no_out @ psi0) ** 2 - abs(U_ref @ psi0) ** 2)) < 1e-12
    # the model's accounting (neither layer) is exactly the reference map without its first and last x phase
    U_inner = ref.dense(inst.n_q, ops=ref.ops[1:-1])
    assert dist(unitary(B.build_step(inst, in_shear=False, out_shear=False).qc), U_inner) < TOL


# --------------------------------------------------------------------------- #
# 6. Recompiled variant, counts, export
# --------------------------------------------------------------------------- #

def test_rz_merge_pass_keeps_the_unitary(ring3):
    inst, ref, U_ref = ring3
    bt = B.build_step(inst)
    mg = B.merge_rz(bt)
    U0, U1 = unitary(bt.qc), unitary(mg.qc)
    assert np.max(np.abs(U0 - U1)) < TOL                 # global phase included
    assert dist(U1, U_ref) < TOL
    assert C.walk(mg)["totals"]["n_rot"] < C.walk(bt)["totals"]["n_rot"]


def test_benchmark_counts_from_the_circuit():
    """Walk the 4^3, K = 8 circuit. Gate set is {h, cx, s, sdg, z, t, tdg, rz}, every rz is a genuine
    arbitrary rotation, and the counts equal the ones stored in counts.json."""
    bt = B.build_step(B.BENCHMARK_2028)
    assert bt.qc.num_qubits == 192
    assert {i.operation.name for i in bt.qc.data} <= {"h", "cx", "s", "sdg", "z", "t", "tdg", "rz"}
    assert all(C.classify(i.operation) == "rot" for i in bt.qc.data if i.operation.name == "rz")
    w = C.walk(bt)
    # qubit count: the 192 field qubits, each acted on, and nothing else (no ancilla, no gate wider than 2)
    assert w["qubits"] == w["qubits_acted_on"] == 64 * 3 and w["ancilla"] == 0
    assert max(len(i.qubits) for i in bt.qc.data) == 2
    stored = json.loads((HERE / "counts.json").read_text())["variants"]["full"]
    assert stored["qubits"] == 192 and stored["ancilla"] == 0
    for key in ("n_rot", "n_t_fixed", "n_toffoli"):
        assert w["totals"][key] == stored[key]
    assert abs(w["totals"]["t_total"] - stored["t_total"]) < 1e-6
    assert w["depth"]["rotation_layers"] == stored["depth"]["rotation_layers"]
    # every site of the translation-invariant lattice carries the same number of rotations per block
    V, n = 64, 3
    assert w["blocks"]["d1_link"]["n_rot"] == 3 * V * n * n
    assert w["blocks"]["d1_site"]["n_rot"] == V * math.comb(n, 2)


def _local(qc, ins, reg):
    return [reg.index(qc.find_bit(q).index) for q in ins.qubits]


def test_benchmark_per_site_blocks_are_the_verified_site_blocks():
    """Full 4^3 circuit: in every per-site block, the gates on each site are, gate for gate, the single-site
    block that tests 1 checks against dense matrices, and no gate there touches two sites."""
    inst = B.BENCHMARK_2028
    bt = B.build_step(inst)
    p, n, qc = bt.params, inst.n_q, bt.qc
    reg0 = list(range(n))
    expect = {
        "shear_in": lambda em: B.site_phase(em, reg0, p.site_gamma[0], p.dx),
        "centering_in": lambda em: B.centering(em, reg0),
        "k1_qft_inv": lambda em: B.qft_swapless(em, reg0, inverse=True),
        "k1_pi2": lambda em: B.site_phase(em, reg0, p.kinetic_u[0], p.dx, reversed_bits=True),
        "k1_qft": lambda em: B.qft_swapless(em, reg0),
        "d1_site": lambda em: B.site_phase(em, reg0, p.site_gamma[1], p.dx),
        "k2_qft_inv": lambda em: B.qft_swapless(em, reg0, inverse=True),
        "k2_pi2": lambda em: B.site_phase(em, reg0, p.kinetic_u[1], p.dx, reversed_bits=True),
        "k2_qft": lambda em: B.qft_swapless(em, reg0),
        "centering_out": lambda em: B.centering(em, reg0, inverse=True),
        "shear_out": lambda em: B.site_phase(em, reg0, p.site_gamma[2], p.dx),
    }
    order = [lab for lab, _, _ in bt.block_ranges()]
    assert order == ["shear_in", "centering_in", "k1_qft_inv", "k1_pi2", "k1_qft", "d1_site", "d1_link",
                     "k2_qft_inv", "k2_pi2", "k2_qft", "centering_out", "shear_out"]
    for lab, start, stop in bt.block_ranges():
        if lab == "d1_link":
            continue
        want = [(i.operation.name, [float(x) for x in i.operation.params], [qc0.find_bit(q).index for q in i.qubits])
                for qc0 in [emit(n, expect[lab])] for i in qc0.data]
        per_site = {x: [] for x in range(inst.V)}
        for ins in qc.data[start:stop]:
            sites = {qc.find_bit(q).index // n for q in ins.qubits}
            assert len(sites) == 1, lab
            x = sites.pop()
            per_site[x].append((ins.operation.name, [float(v) for v in ins.operation.params],
                                _local(qc, ins, B.site_register(inst, x))))
        for x in range(inst.V):
            assert per_site[x] == want, (lab, x)


def test_benchmark_diagonal_block_decodes_to_the_hamiltonian():
    """Full 4^3 circuit: read every ZZ rotation of the D block back into a quadratic form in the Z operators
    and compare with exp(-i x^T Q x / 2), x = -(dx/2) sum_q 2^q Z_q, Q from the reference map (its own
    Laplacian and frames). All 1920 rotations, 192 qubits."""
    inst = B.BENCHMARK_2028
    bt = B.build_step(inst)
    qc, n, V = bt.qc, inst.n_q, inst.V
    ref = RefMap(inst.dims, inst.m, inst.h)
    xops = [v for typ, v in ref.ops if typ == "x"]
    assert len(xops) == 5                                   # shear, shear + D + shear, shear
    Q = xops[1] + xops[2] + xops[3]
    dx = math.sqrt(2 * PI / inst.K)
    want = np.zeros((V * n, V * n))                          # coefficient of Z_a Z_b, a < b, in the exponent
    for x in range(V):
        for y in range(x, V):
            for q in range(n):
                for qp in range(n):
                    a, b = x * n + q, y * n + qp
                    if a < b:
                        want[a, b] += Q[x, y] * dx * dx / 4 * 2.0 ** (q + qp)
    got = np.zeros_like(want)
    n_rot = 0
    for lab, start, stop in bt.block_ranges():
        if lab not in ("d1_site", "d1_link"):
            continue
        data = qc.data[start:stop]
        assert len(data) % 3 == 0
        for k in range(0, len(data), 3):
            c1, rz, c2 = data[k:k + 3]
            assert (c1.operation.name, rz.operation.name, c2.operation.name) == ("cx", "rz", "cx")
            a, b = (qc.find_bit(q).index for q in c1.qubits)
            assert [qc.find_bit(q).index for q in c2.qubits] == [a, b]
            assert qc.find_bit(rz.qubits[0]).index == b
            got[min(a, b), max(a, b)] += float(rz.operation.params[0]) / 2      # exp(-i theta ZZ / 2)
            n_rot += 1
    assert n_rot == 1920
    assert np.count_nonzero(want) == 1920
    assert np.max(np.abs(got - want)) < 1e-12
    # each site couples to its 6 distinct neighbours
    assert all(np.count_nonzero(Q[x]) == 7 for x in range(V))


def test_model_accounting_variant_matches_the_model_line_items():
    """Without the two outer shear layers the walked circuit has the model's rotation and T counts
    (shot_counts(a, '2028', 1)). With them it has 2 V C(n_q,2) more rotations."""
    a = model.Assumptions()
    sc = model.shot_counts(a, "2028", 1)
    w0 = C.walk(B.build_step(B.BENCHMARK_2028, in_shear=False, out_shear=False))
    assert w0["totals"]["n_rot"] == sc["n_rot"]
    assert w0["totals"]["n_t_fixed"] == sum(sc["t_exact"].values())
    assert abs(w0["totals"]["t_total"] - sc["t_total"]) < 1e-6
    w1 = C.walk(B.build_step(B.BENCHMARK_2028))
    assert w1["totals"]["n_rot"] - sc["n_rot"] == 2 * 64 * 3


def test_qasm_export_round_trip():
    """The exported OpenQASM 3 file loads back to a circuit with the same qubits and the same counts."""
    from qiskit import qasm3
    path = HERE / C.QASM_NAME
    assert path.stat().st_size < 20e6
    qc = qasm3.loads(path.read_text())
    bt = B.build_step(B.BENCHMARK_2028)
    assert qc.num_qubits == bt.qc.num_qubits == 192
    assert len(qc.data) == len(bt.qc.data)
    ref = C.walk(bt)["totals"]
    got = C.walk(B.Built(qc=qc, labels=["all"] * len(qc.data), params=bt.params))["totals"]
    for key in ("n_rot", "n_t_fixed", "n_toffoli", "cx"):
        assert got[key] == ref[key]
    for a_, b_ in zip(qc.data, bt.qc.data):                       # every instruction, not a sample
        assert a_.operation.name == b_.operation.name
        assert [qc.find_bit(q).index for q in a_.qubits] == [bt.qc.find_bit(q).index for q in b_.qubits]
        if a_.operation.params:
            assert abs(float(a_.operation.params[0]) - float(b_.operation.params[0])) < 1e-15
