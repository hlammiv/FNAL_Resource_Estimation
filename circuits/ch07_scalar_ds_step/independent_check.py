#!/usr/bin/env python3
"""Second, independent check of the Ch. 7 2028 step circuit (written by the reviewer, not the builder).

Run:  OMP_NUM_THREADS=4 python circuits/ch07_scalar_ds_step/independent_check.py      (about 3 s, 0.2 GB)

It shares no code with test_circuit.py or count.py. What it uses from build.py is only the circuit under test
(build_step and the block emitters). Differences from the builder's references:
  - the conjugate variable p is the generator of the grid shift, P = i logm(T) / dx with T |j> = |j+1> and an
    antiperiodic wrap, not F^dag x F from the DFT formula
  - the step on a ring is a product of scipy expm of dense many-site operators (kron), with the frame
    coefficients solved here from the Hamiltonian
  - the in-state covariance is summed over complex plane waves here (rings and 4^3)
  - the gate counter reads the .qasm text with a regular expression, with no qiskit and no count.py
  - the full 192-qubit circuit is checked from the .qasm text as an exact factorization
        U = (product over sites of an 8x8 unitary) . (inter-site ZZ rotations) . (product over sites of an 8x8)
    with every 8x8 factor compared with a dense reference and every ZZ angle with the Hamiltonian.

T totals use estimates.common (imported). Tolerance 1e-9 after removing one global phase.
"""
from __future__ import annotations

import itertools
import json
import math
import re
import sys
from pathlib import Path

import numpy as np
from scipy.linalg import expm, logm
from scipy.special import h1vp, hankel1

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
for _p in (str(REPO), str(HERE)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from estimates.common import EPS_SYN, eps_rot_for, t_per_rotation, toffoli_t      # noqa: E402

QASM = HERE / "ch07_scalar_ds_step_4x4x4_K8.qasm"
TOL = 1e-9
PI = math.pi


# --------------------------------------------------------------------------- #
# Gate list from the .qasm text, and the counter
# --------------------------------------------------------------------------- #

def parse_qasm(path):
    ops, nq = [], None
    pat = re.compile(r"^(\w+)(?:\(([^)]*)\))?\s+(.*);$")
    for line in Path(path).read_text().splitlines():
        line = line.strip()
        if not line or line.startswith(("OPENQASM", "include")):
            continue
        m = re.match(r"qubit\[(\d+)\]\s+\w+;", line)
        if m:
            assert nq is None, "more than one register"
            nq = int(m.group(1))
            continue
        m = pat.match(line)
        assert m, line
        name, par, args = m.groups()
        theta = float(eval(par, {"pi": PI, "__builtins__": {}})) if par else None
        ops.append((name, theta, [int(x) for x in re.findall(r"\[(\d+)\]", args)]))
    return nq, ops


def kind(name, theta):
    """'c' Clifford, 't' exact T, 'rot' arbitrary rotation, 'toff' Toffoli. Anything else is an error, so a
    multi-controlled gate or a controlled rotation cannot pass as one gate."""
    if name in ("h", "cx", "cz", "s", "sdg", "z", "x"):
        return "c"
    if name in ("t", "tdg"):
        return "t"
    if name in ("ccx", "ccz"):
        return "toff"
    if name == "rz":
        k = theta / (PI / 4)
        if abs(k - round(k)) < 1e-9:
            return "t" if round(k) % 2 else "c"
        return "rot"
    raise ValueError(f"gate {name!r} is not in the logical gate set")


def count(nq, ops):
    c = dict(c=0, t=0, toff=0, rot=0)
    touched = set()
    for name, theta, qs in ops:
        c[kind(name, theta)] += 1
        touched.update(qs)
    eps = eps_rot_for(c["rot"], EPS_SYN)
    tr = t_per_rotation(eps)

    def depth(w):
        d = [0.0] * nq
        for name, theta, qs in ops:
            v = max(d[q] for q in qs) + w[kind(name, theta)]
            for q in qs:
                d[q] = v
        return max(d)
    return dict(qubits_declared=nq, qubits_acted_on=len(touched), widest_gate=max(len(o[2]) for o in ops),
                n_toffoli=c["toff"], n_rot=c["rot"], n_t_fixed=c["t"], eps_rot=eps, t_per_rotation=tr,
                t_total=c["rot"] * tr + c["t"] + toffoli_t(c["toff"], "textbook"),
                nonclifford_layers=depth(dict(c=0, t=1, rot=1, toff=1)),
                rotation_layers=depth(dict(c=0, t=0, rot=1, toff=0)),
                t_depth_report=depth(dict(c=0, t=1, rot=tr, toff=1)))


# --------------------------------------------------------------------------- #
# Operators from their definitions
# --------------------------------------------------------------------------- #

def grid_ops(n_q):
    K = 2 ** n_q
    dx = math.sqrt(2 * PI / K)
    X = np.diag(dx * (np.arange(K) - (K - 1) / 2))
    T = np.zeros((K, K), complex)
    for j in range(K):
        T[(j + 1) % K, j] = 1.0 if j + 1 < K else -1.0
    P = 1j * logm(T) / dx                       # exp(-i dx P) = T
    P = (P + P.conj().T) / 2
    assert np.abs(np.linalg.eigvalsh(P) - np.diag(X)).max() < 1e-10
    return K, dx, X, P


def pdist(U, W):
    ov = np.vdot(W, U)
    return float(np.max(np.abs(U - ov / abs(ov) * W)))


def bunch_davies(m, k2):
    """(<phi phi>, Re<phi pi>, <pi pi>) of one lattice mode at eta_in = -1, a = 1. k = 0: a = 1 ground state."""
    if k2 < 1e-14:
        return 1 / (2 * m), 0.0, m / 2
    nu = math.sqrt(2.25 - m * m)
    ph = np.exp(1j * PI * (nu + 0.5) / 2)
    z = math.sqrt(k2)
    q = math.sqrt(PI) / 2 * hankel1(nu, z) * ph
    p = math.sqrt(PI) / 2 * ph * (-1.5 * hankel1(nu, z) - z * h1vp(nu, z))
    return abs(q) ** 2, float((q * np.conj(p)).real), abs(p) ** 2


def balance(sxx, sxp, spp):
    s = math.sqrt(math.sqrt(sxx * spp - sxp ** 2) / sxx)
    return np.array([[s, 0.0], [-sxp / (s * sxx), 1 / s]])


def kinetic_in_frame(Sin, Sout, v):
    """S_out [[1, v], [0, 1]] S_in^-1 = Lx(g2) Up(u) Lx(g1), Lx(g) = exp(-i g x^2 / 2): p -> p - g x,
    Up(u) = exp(-i u p^2 / 2): x -> x + u p. Returns g1 (acts first), u, g2."""
    A = Sout @ np.array([[1.0, v], [0.0, 1.0]]) @ np.linalg.inv(Sin)
    u = A[0, 1]
    g1, g2 = (1 - A[0, 0]) / u, (1 - A[1, 1]) / u
    assert abs(A[1, 0] - (u * g1 * g2 - g1 - g2)) < 1e-12
    return g1, u, g2


def step_coefficients(dims, m, h=0.5):
    """One step K(h/2) D(h) K(h/2) at a_m = e^{h/2} in the per-site frame, mode by mode."""
    b = math.sqrt(12 / (4 - m * m))
    am = math.exp(h / 2)
    v = h / 2 / am ** 3
    modes = []
    for nvec in itertools.product(*[range(L) for L in dims]):
        k2 = sum((2 / b * math.sin(PI * c / L)) ** 2 for c, L in zip(nvec, dims))
        sxx, sxp, spp = bunch_davies(m, k2)
        modes.append((k2, np.array([[sxx, sxp], [sxp, spp]])))

    def site(ms):
        C = sum(c for _, c in ms) / len(ms)
        return balance(C[0, 0], C[0, 1], C[1, 1])
    MK = np.array([[1.0, v], [0.0, 1.0]])
    c1 = [(k2, MK @ c @ MK.T) for k2, c in modes]
    c2 = []
    for k2, c in c1:
        MD = np.array([[1.0, 0.0], [-h * (am ** 3 * m * m + am * k2), 1.0]])
        c2.append((k2, MD @ c @ MD.T))
    c3 = [(k2, MK @ c @ MK.T) for k2, c in c2]
    S0, S1, S3 = site(modes), site(c1), site(c3)
    g1a, ua, g2a = kinetic_in_frame(S0, S1, v)
    g1b, ub, g2b = kinetic_in_frame(S1, S3, v)
    s2 = S1[0, 0] ** 2
    n_dir = sum(1 for L in dims if L > 1)
    return dict(b=b, g_in=g1a, u1=ua, g_mid=g2a + h * (am ** 3 * m * m + 2 * n_dir * am / b ** 2) / s2 + g1b,
                u2=ub, g_out=g2b, q_link=-h * am / b ** 2 / s2)


# --------------------------------------------------------------------------- #
# Checks
# --------------------------------------------------------------------------- #

def check_counts():
    nq, ops = parse_qasm(QASM)
    r = count(nq, ops)
    stored = json.loads((HERE / "counts.json").read_text())["variants"]["full"]
    assert r["qubits_declared"] == r["qubits_acted_on"] == stored["qubits"] == 192
    assert r["widest_gate"] == 2 and r["n_toffoli"] == stored["n_toffoli"] == 0
    assert (r["n_rot"], r["n_t_fixed"]) == (stored["n_rot"], stored["n_t_fixed"])
    assert abs(r["t_total"] - stored["t_total"]) < 1e-6
    assert abs(r["t_depth_report"] - stored["depth"]["t_depth_report"]) < 1e-6
    assert r["rotation_layers"] == stored["depth"]["rotation_layers"]
    assert r["nonclifford_layers"] == stored["depth"]["nonclifford_layers"]
    return r


def check_ring(L=3, n_q=3, m=1.0):
    """Circuit unitary of the full step on a ring against a product of dense matrix exponentials."""
    import build as B
    from qiskit.quantum_info import Operator
    K, dx, X, P = grid_ops(n_q)
    out = {}
    for u in (0.118, -0.73):
        em = B.Emitter(n_q)
        B.kinetic_block(em, list(range(n_q)), u, dx)
        out[f"kinetic block u={u}"] = pdist(Operator(em.qc).data, expm(-0.5j * u * P @ P))
    c = step_coefficients((L,), m)

    def at(A, x):                                       # qiskit order: site 0 is the least significant factor
        M = np.eye(1)
        for s in reversed(range(L)):
            M = np.kron(M, A if s == x else np.eye(K))
        return M
    Xs = [at(X, x) for x in range(L)]
    X2 = sum(Xx @ Xx for Xx in Xs)
    P2 = sum(at(P @ P, x) for x in range(L))
    hop = sum(Xs[x] @ Xs[(x + 1) % L] for x in range(L))            # each link once
    U = np.eye(K ** L, dtype=complex)
    for gen in (c["g_in"] * X2 / 2, c["u1"] * P2 / 2, c["g_mid"] * X2 / 2 + c["q_link"] * hop,
                c["u2"] * P2 / 2, c["g_out"] * X2 / 2):
        U = expm(-1j * gen) @ U
    bt = B.build_step(B.Instance(dims=(L,), n_q=n_q, m=m))
    out[f"full step, ring L={L}, K={K}"] = pdist(Operator(bt.qc).data, U)
    return out


def check_full_circuit(n=3, L=4, m=1.0):
    """The 192-qubit .qasm as (x) post . ZZ links . (x) pre, every factor against its reference."""
    nq, ops = parse_qasm(QASM)
    V, K = L ** 3, 2 ** n
    assert nq == V * n
    _, dx, X, P = grid_ops(n)
    one = {"h": np.array([[1, 1], [1, -1]]) / math.sqrt(2), "s": np.diag([1, 1j]), "sdg": np.diag([1, -1j]),
           "t": np.diag([1, np.exp(0.25j * PI)]), "tdg": np.diag([1, np.exp(-0.25j * PI)]), "z": np.diag([1, -1])}

    def apply(U, name, theta, loc):
        if name == "cx":
            ctl, tgt = loc
            M = np.zeros((K, K))
            for j in range(K):
                M[j ^ ((j >> ctl & 1) << tgt), j] = 1
        else:
            g = np.diag([np.exp(-0.5j * theta), np.exp(0.5j * theta)]) if name == "rz" else one[name]
            M = np.eye(1)
            for q in reversed(range(n)):
                M = np.kron(M, g if q == loc[0] else np.eye(2))
        return M @ U
    first, last = {}, {}
    for i, (_, _, qs) in enumerate(ops):
        if len({q // n for q in qs}) > 1:
            for q in qs:
                first.setdefault(q // n, i)
                last[q // n] = i
    assert len(first) == V
    pre = [np.eye(K, dtype=complex) for _ in range(V)]
    post = [np.eye(K, dtype=complex) for _ in range(V)]
    zz = {}
    i = 0
    while i < len(ops):
        name, theta, qs = ops[i]
        sites = {q // n for q in qs}
        if len(sites) == 1:
            x = sites.pop()
            loc = [q - x * n for q in qs]
            assert i < first[x] or i > last[x], f"intra-site gate inside the link block at {i}"
            if i < first[x]:
                pre[x] = apply(pre[x], name, theta, loc)
            else:
                post[x] = apply(post[x], name, theta, loc)
            i += 1
        else:                                           # exp(-i theta Z_a Z_b / 2) as cx, rz, cx
            (n1, _, q1), (n2, th, q2), (n3, _, q3) = ops[i:i + 3]
            assert (n1, n2, n3) == ("cx", "rz", "cx") and q1 == q3 and q2 == [q1[1]]
            key = (min(q1), max(q1))
            zz[key] = zz.get(key, 0.0) + th
            i += 3
    c = step_coefficients((L, L, L), m)
    Dc = np.diag(np.exp(2j * PI * (K - 1) / 2 * np.arange(K) / K))      # centering layer, commutes with x phases
    Lx = lambda g: expm(-0.5j * g * X @ X)
    Up = lambda u: expm(-0.5j * u * P @ P)
    pre_ref = Dc @ Lx(c["g_mid"]) @ Up(c["u1"]) @ Lx(c["g_in"])
    post_ref = Lx(c["g_out"]) @ Up(c["u2"]) @ Dc.conj()
    want = {}
    for x in range(V):
        cx_ = (x // (L * L), (x // L) % L, x % L)
        for ax in range(3):
            cy = list(cx_)
            cy[ax] = (cy[ax] + 1) % L
            y = cy[0] * L * L + cy[1] * L + cy[2]
            for q in range(n):
                for qp in range(n):
                    key = (min(x * n + q, y * n + qp), max(x * n + q, y * n + qp))
                    want[key] = want.get(key, 0.0) + c["q_link"] * dx * dx * 2.0 ** (q + qp) / 2
    assert set(want) == set(zz) and len(zz) == 3 * V * n * n
    return {"192 qubits: site factor before the links, worst of 64": max(pdist(pre[x], pre_ref) for x in range(V)),
            "192 qubits: site factor after the links, worst of 64": max(pdist(post[x], post_ref) for x in range(V)),
            "192 qubits: 1728 link ZZ angles, largest difference": max(abs(want[k] - zz[k]) for k in want)}


def check_gaussian(m=1.0, L=4, n_q=3, h=0.5):
    """Sampled in-state Gaussian through the circuit on a 4-site ring, covariance against the continuum map."""
    import build as B
    from qiskit.quantum_info import Statevector
    K, dx, X, P = grid_ops(n_q)
    b = math.sqrt(12 / (4 - m * m))
    am = math.exp(h / 2)
    v = h / 2 / am ** 3
    I, Z = np.eye(L), np.zeros((L, L))
    lap = (2 * I - np.roll(I, 1, 0) - np.roll(I, -1, 0)) / b ** 2
    C = np.zeros((2 * L, 2 * L))
    for nn in range(L):
        sxx, sxp, spp = bunch_davies(m, (2 / b * math.sin(PI * nn / L)) ** 2)
        w = np.exp(2j * PI * nn * np.arange(L) / L)
        pr = np.outer(w, w.conj()).real / L
        C += np.block([[sxx * pr, sxp * pr], [sxp * pr, spp * pr]])
    M = np.block([[I, v * I], [Z, I]]) @ np.block([[I, Z], [-h * (am ** 3 * m * m * I + am * lap), I]]) \
        @ np.block([[I, v * I], [Z, I]])
    C3 = M @ C @ M.T
    S0, S3 = balance(C[0, 0], C[0, L], C[L, L]), balance(C3[0, 0], C3[0, L], C3[L, L])
    big = lambda S: np.kron(S, I)
    Cr = big(S0) @ C @ big(S0).T
    A, Bm = Cr[:L, :L], Cr[:L, L:]
    G = np.linalg.inv(A) / 4 - 0.5j * np.linalg.inv(A) @ Bm        # psi ~ exp(-x^T G x)
    G = (G + G.T) / 2
    xs = np.stack([g.reshape(-1) for g in np.meshgrid(*([np.diag(X)] * L), indexing="ij")], 0)
    psi = np.exp(-np.einsum("in,ij,jn->n", xs, G, xs))
    psi /= np.linalg.norm(psi)
    flip = lambda vec: vec.reshape((K,) * L).transpose(*reversed(range(L))).reshape(-1)
    bt = B.build_step(B.Instance(dims=(L,), n_q=n_q, m=m, h=h))
    out = flip(np.asarray(Statevector(flip(psi)).evolve(bt.qc).data)).reshape((K,) * L)
    act = lambda Aop, s: np.moveaxis(np.tensordot(Aop, out, axes=([1], [s])), 0, s)
    vecs = [act(X, s) for s in range(L)] + [act(P, s) for s in range(L)]
    cov = np.array([[np.real(np.vdot(a_, b_)) for b_ in vecs] for a_ in vecs])
    lab = big(np.linalg.inv(S3)) @ cov @ big(np.linalg.inv(S3)).T
    worst = 0.0
    for nn in range(L):
        w = np.exp(2j * PI * nn * np.arange(L) / L) / math.sqrt(L)
        f = lambda Cm, a_, b_: float(np.real(w.conj() @ Cm[a_ * L:(a_ + 1) * L, b_ * L:(b_ + 1) * L] @ w))
        d = (f(lab, 0, 0), f(lab, 1, 1), f(lab, 0, 1))
        r = (f(C3, 0, 0), f(C3, 1, 1), f(C3, 0, 1))
        worst = max(worst, abs(d[0] / r[0] - 1), abs(d[1] / r[1] - 1), abs(d[2] - r[2]) / math.sqrt(r[0] * r[1]))
    return worst


def main():
    r = check_counts()
    print("counts from the .qasm text:", {k: (float(f"{v:.6g}") if isinstance(v, float) else v) for k, v in r.items()})
    dev = {}
    dev.update(check_ring())
    dev.update(check_full_circuit())
    for k, v in dev.items():
        print(f"  {v:9.2e}  {k}")
        assert v < TOL, k
    for m in (1.0, 1.4):
        e = check_gaussian(m)
        print(f"  {e:9.2e}  Gaussian through the circuit, 4-site ring, K = 8, m = {m}: largest shell error "
              "against the continuum map")
        assert e < 2e-3
    print("independent check passed")


if __name__ == "__main__":
    main()
