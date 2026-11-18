"""Independent check of the exported QASM files of the 2T link-step circuits.

Nothing here imports build.py. The group is built from Hamilton quaternion arithmetic, the gates are
read from the QASM text with a small parser, and the circuits are run with the simulators below.

What is checked:
  A. Each 4x4 lattice file is cut into contiguous blocks that exhaust the gate list. Every electric block
     (one link, up to two ancillas) is simulated as a dense unitary and compared, on the 24 group states
     with the ancillas in |0>, with exp(-i theta_e Laplacian). Every plaquette block is run on all 24^4
     configurations of its four links (it is a permutation with phases) and compared with
     exp(i theta_b Re Tr(U_x(s) U_y(s+x) U_x(s+y)^-1 U_y(s)^-1)). Each plaquette of the torus must appear
     once per step and each link must get one electric block per step. Because the blocks are contiguous,
     exhaust the file and are each checked on the whole physical domain of the qubits they touch, the
     file as a whole is the product  prod_steps [prod_plaq exp(i theta_b Re Tr) prod_links exp(-i theta_e E^2)]
     on the physical subspace, up to a global phase.
  B. The 22-qubit single-plaquette file is run as a statevector on a random state of the gauge register.
  C. Counts: qubits touched, ccx, rz that are not multiples of pi/4, t and tdg, and an as-soon-as-possible
     non-Clifford depth (Toffoli 1 layer, T 1 layer, rotation c_T layers), all from the QASM text.

Run from the repository root:   python circuits/ch05_bt_link_step/verify_independent.py
"""
from __future__ import annotations

import math
import os
import re
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
TOL = 1e-9

# ------------------------------------------------------------------------------------------------
# The group from quaternions
# ------------------------------------------------------------------------------------------------

def qmul(a, b):
    w1, x1, y1, z1 = a
    w2, x2, y2, z2 = b
    return np.array([w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2,
                     w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2,
                     w1 * y2 - x1 * z2 + y1 * w2 + z1 * x2,
                     w1 * z2 + x1 * y2 - y1 * x2 + z1 * w2])


def qpow(a, n):
    out = np.array([1.0, 0, 0, 0])
    for _ in range(n):
        out = qmul(out, a)
    return out


Q1, QI, QJ = np.array([1.0, 0, 0, 0]), np.array([0, 1.0, 0, 0]), np.array([0, 0, 1.0, 0])
QL = -0.5 * np.array([1.0, 1, 1, 1])            # l = -(1 + i + j + k)/2

# register value = m + 2n + 4o + 8p + 16q for g = (-1)^m i^n j^o l^(p+2q)  (arXiv:2208.12309 Eq. (pres))
IDX = [v for v in range(32) if not ((v >> 3) & 1 and (v >> 4) & 1)]
ELEM = {}
for v in IDX:
    m, n, o, p, q = [(v >> k) & 1 for k in range(5)]
    ELEM[v] = (-1) ** m * qmul(qmul(qpow(QI, n), qpow(QJ, o)), qpow(QL, p + 2 * q))
POS = {v: k for k, v in enumerate(IDX)}


def find(qt):
    hits = [v for v in IDX if np.abs(ELEM[v] - qt).max() < 1e-9]
    assert len(hits) == 1, "not a unique group element"
    return hits[0]


MUL = np.array([[POS[find(qmul(ELEM[a], ELEM[b]))] for b in IDX] for a in IDX])      # positions
INV = np.array([POS[find(ELEM[a] * np.array([1, -1, -1, -1]))] for a in IDX])
RETR = np.array([2.0 * ELEM[a][0] for a in IDX])                                     # Re Tr in the 2d irrep


def check_group():
    # the 24 Hurwitz units: +-1, +-i, +-j, +-k and (+-1 +-i +-j +-k)/2
    want = set()
    for k in range(4):
        for s in (1, -1):
            e = [0.0] * 4
            e[k] = s
            want.add(tuple(e))
    for s in range(16):
        want.add(tuple(0.5 * (1 if (s >> k) & 1 else -1) for k in range(4)))
    got = {tuple(np.round(ELEM[v], 9) + 0.0) for v in IDX}
    assert got == want and len(got) == 24
    assert all(MUL[INV[a], a] == POS[0] for a in range(24))
    # the paper's worked example: (1/2)[[eta, eta*], [-eta, eta*]] is stored as |23>
    # with i = diag(i, -i), j = [[0, -1], [1, 0]], k = i j this matrix is (1 + i - j + k)/2
    assert find(np.array([0.5, 0.5, -0.5, 0.5])) == 23


def laplacian():
    """sum over the 8 elements with Re Tr = 1 of (1 - L_gamma), scaled so the 2d irrep carries 3/4."""
    gam = [a for a in range(24) if abs(RETR[a] - 1.0) < 1e-9]
    assert len(gam) == 8
    lap = np.zeros((24, 24))
    lap_r = np.zeros((24, 24))
    for g in gam:
        for h in range(24):
            lap[h, h] += 1
            lap[MUL[g, h], h] -= 1
            lap_r[h, h] += 1
            lap_r[MUL[h, g], h] -= 1
    assert np.allclose(lap, lap.T) and np.allclose(lap, lap_r)       # gauge invariant on both ends
    # Casimir of the faithful irrep: |Gamma| - (sum of the quaternions) = 8 - 4
    s = sum(ELEM[IDX[g]] for g in gam)
    assert np.allclose(s, [4, 0, 0, 0])
    return lap * (0.75 / (len(gam) - s[0]))


def expm_herm(h, theta):
    w, v = np.linalg.eigh(h)
    return (v * np.exp(-1j * theta * w)) @ v.conj().T


def plaquette_retr():
    """R[a, b, c, d] = Re Tr(g_a g_b g_c^-1 g_d^-1) over positions."""
    ab = MUL                                       # [a, b]
    cd = MUL[INV][:, INV]                          # [c, d] -> c^-1 d^-1
    return RETR[MUL[ab[:, :, None, None], cd[None, None, :, :]]]


# ------------------------------------------------------------------------------------------------
# QASM
# ------------------------------------------------------------------------------------------------

LINE = re.compile(r"^(\w+)(?:\(([^)]*)\))?\s+(.*);$")
KNOWN = {"x", "z", "s", "sdg", "h", "cx", "cz", "swap", "t", "tdg", "ccx", "rz"}


def parse(path):
    nq, gates = None, []
    for raw in open(path):
        line = raw.strip()
        if not line or line.startswith(("OPENQASM", "include")):
            continue
        m = re.match(r"^qubit\[(\d+)\]\s+q;$", line)
        if m:
            nq = int(m.group(1))
            continue
        m = LINE.match(line)
        assert m, f"cannot read {line!r}"
        name, par, args = m.groups()
        assert name in KNOWN, f"gate {name!r} is outside the logical gate set"
        qs = tuple(int(x) for x in re.findall(r"q\[(\d+)\]", args))
        assert len(qs) == {"cx": 2, "cz": 2, "swap": 2, "ccx": 3}.get(name, 1) and len(set(qs)) == len(qs)
        gates.append((name, float(par) if par is not None else None, qs))
    return nq, gates


def is_quarter(theta):
    k = theta / (math.pi / 4)
    return abs(k - round(k)) < 1e-9


def counts(nq, gates):
    c = {"toffoli": 0, "rot": 0, "t": 0, "clifford": 0}
    used = set()
    for name, par, qs in gates:
        used.update(qs)
        if name == "ccx":
            c["toffoli"] += 1
        elif name in ("t", "tdg"):
            c["t"] += 1
        elif name == "rz":
            assert not is_quarter(par), "an rz by a multiple of pi/4 is not an arbitrary rotation"
            c["rot"] += 1
        else:
            c["clifford"] += 1
    c["qubits_declared"], c["qubits_used"] = nq, len(used)
    return c


def depth(nq, gates, c_t):
    free = [0.0] * nq
    for name, par, qs in gates:
        w = 1.0 if name in ("ccx", "t", "tdg") else (c_t if name == "rz" else 0.0)
        end = max(free[q] for q in qs) + w
        for q in qs:
            free[q] = end
    return max(free)


# ------------------------------------------------------------------------------------------------
# Simulators
# ------------------------------------------------------------------------------------------------

PH = {"z": -1.0, "s": 1j, "sdg": -1j, "t": np.exp(0.25j * np.pi), "tdg": np.exp(-0.25j * np.pi)}


def run_dense(psi, gates, loc):
    """psi has shape (batch, 2, ..., 2); axis 1 + loc[q] is qubit q."""
    for name, par, qs in gates:
        ax = [1 + loc[q] for q in qs]
        if name == "h":
            v = np.moveaxis(psi, ax[0], 0)
            a, b = v[0].copy(), v[1].copy()
            v[0], v[1] = (a + b) / math.sqrt(2), (a - b) / math.sqrt(2)
        elif name == "x":
            v = np.moveaxis(psi, ax[0], 0)
            a = v[0].copy()
            v[0], v[1] = v[1], a
        elif name in PH:
            np.moveaxis(psi, ax[0], 0)[1] *= PH[name]
        elif name == "rz":
            v = np.moveaxis(psi, ax[0], 0)
            v[0] *= np.exp(-0.5j * par)
            v[1] *= np.exp(0.5j * par)
        elif name == "cx":
            v = np.moveaxis(psi, ax, [0, 1])
            a = v[1, 0].copy()
            v[1, 0], v[1, 1] = v[1, 1], a
        elif name == "cz":
            np.moveaxis(psi, ax, [0, 1])[1, 1] *= -1.0
        elif name == "ccx":
            v = np.moveaxis(psi, ax, [0, 1, 2])
            a = v[1, 1, 0].copy()
            v[1, 1, 0], v[1, 1, 1] = v[1, 1, 1], a
        elif name == "swap":
            v = np.moveaxis(psi, ax, [0, 1])
            a = v[0, 1].copy()
            v[0, 1], v[1, 0] = v[1, 0], a
        else:
            raise ValueError(name)
    return psi


def run_bits(bits, angle, gates, loc):
    """Permutation-with-phase circuits on many basis states: bits[state, qubit], angle[state]."""
    for name, par, qs in gates:
        c = [loc[q] for q in qs]
        if name == "x":
            bits[:, c[0]] ^= 1
        elif name == "cx":
            bits[:, c[1]] ^= bits[:, c[0]]
        elif name == "ccx":
            bits[:, c[2]] ^= bits[:, c[0]] & bits[:, c[1]]
        elif name == "swap":
            bits[:, c] = bits[:, c[::-1]]
        elif name == "cz":
            angle += np.pi * (bits[:, c[0]] & bits[:, c[1]])
        elif name in PH:
            angle += np.angle(PH[name]) * bits[:, c[0]]
        elif name == "rz":
            angle += par * (bits[:, c[0]].astype(float) - 0.5)
        else:
            raise ValueError(f"{name} is not a permutation or a phase")
    return bits, angle


def phase_error(out, ref):
    """max |out - ref| after removing one global phase."""
    ov = np.vdot(ref.ravel(), out.ravel())
    assert abs(ov) > 1e-6
    return float(np.abs(out * np.conj(ov / abs(ov)) - ref).max())


# ------------------------------------------------------------------------------------------------
# A. lattice files, block by block
# ------------------------------------------------------------------------------------------------

def torus_plaquettes(L):
    """frozenset of links -> ordered (U_x(s), U_y(s+x), U_x(s+y), U_y(s)); link (x, y, mu) = 2(x + L y) + mu."""
    lk = lambda x, y, mu: 2 * ((x % L) + L * (y % L)) + mu
    out = {}
    for x in range(L):
        for y in range(L):
            order = (lk(x, y, 0), lk(x + 1, y, 1), lk(x, y + 1, 0), lk(x, y, 1))
            out[frozenset(order)] = order
    assert len(out) == L * L
    return out


def localize(gates, loc):
    return tuple((name, par, tuple(loc[q] for q in qs)) for name, par, qs in gates)


def electric_block_error(local_gates, n_anc, u_ref):
    """(error against u_ref, norm outside the register) of a block on qubits 0..4 = link, 5.. = ancillas."""
    n = 5 + n_anc
    psi = np.zeros((24,) + (2,) * n, dtype=complex)
    for col, v in enumerate(IDX):
        psi[(col,) + tuple((v >> k) & 1 for k in range(5)) + (0,) * n_anc] = 1.0
    psi = run_dense(psi, local_gates, {q: q for q in range(n)})
    # index [input, anc..., q, p, o, n, m] -> [input, ancilla value, register value]
    flat = np.transpose(psi, [0] + list(range(n, 0, -1))).reshape(24, 2 ** n_anc, 32).copy()
    blk = flat[:, 0, :][:, IDX].T.copy()                     # [output position, input position]
    flat[:, 0, IDX] = 0.0                                    # what is left is outside the register
    return phase_error(blk, u_ref), float(np.linalg.norm(flat))


_ALL = {}


def all_plaquette_bits():
    """All 24^4 configurations of four links as 20 bits (bit b of link k in column 5k + b)."""
    if "bits" not in _ALL:
        ar = np.array(IDX)
        grid = np.stack(np.meshgrid(ar, ar, ar, ar, indexing="ij"), axis=-1).reshape(-1, 4)
        _ALL["bits"] = ((grid[:, :, None] >> np.arange(5)[None, None, :]) & 1).reshape(-1, 20).astype(np.int8)
    return _ALL["bits"]


def plaquette_block_error(local_gates, perm, theta_b, retr):
    """Phase error of a block on four links (local link k on qubits 5k..5k+4) over all 24^4 configurations.
    perm[r] = local link that plays role r in Re Tr(g_0 g_1 g_2^-1 g_3^-1). Raises if a register is not
    restored."""
    bits0 = all_plaquette_bits()
    bits, angle = run_bits(bits0.copy(), np.zeros(len(bits0)), local_gates, {q: q for q in range(20)})
    assert np.array_equal(bits, bits0), "registers not restored"
    ref = np.transpose(retr, np.argsort(perm)).reshape(-1)   # axes in the order of the local links
    return float(np.abs(np.exp(1j * (angle - angle[0])) - np.exp(1j * theta_b * (ref - ref[0]))).max())


def check_lattice_file(path, L, theta_e, theta_b, steps):
    nq, gates = parse(path)
    n_links = 2 * L * L
    n_data = 5 * n_links
    u_e = expm_herm(laplacian(), theta_e)
    retr = plaquette_retr()
    plaqs = torus_plaquettes(L)
    rng = np.random.default_rng(7)
    ar = np.array(IDX)
    cache = {}                                   # identical gate lists on relabelled qubits are checked once
    sequence, worst_e, worst_leak, worst_b = [], 0.0, 0.0, 0.0
    i = 0
    while i < len(gates):
        if gates[i][0] == "h":
            # electric block: gates on one link and on ancillas, up to the first gate on another link
            assert gates[i][2][0] < n_data
            link = gates[i][2][0] // 5
            j, anc = i, []
            while j < len(gates):
                qs = gates[j][2]
                if any(q < n_data and q // 5 != link for q in qs):
                    break
                anc += [q for q in qs if q >= n_data and q not in anc]
                j += 1
            assert len(anc) <= 2
            loc = {5 * link + k: k for k in range(5)}
            loc.update({a: 5 + k for k, a in enumerate(anc)})
            key = ("E", localize(gates[i:j], loc), len(anc))
            if key not in cache:
                cache[key] = electric_block_error(key[1], len(anc), u_e)
            err, leak = cache[key]
            assert leak < TOL, f"electric block on link {link}: norm {leak:.2e} outside the register"
            assert err < TOL, f"electric block on link {link}: error {err:.2e}"
            worst_e, worst_leak = max(worst_e, err), max(worst_leak, leak)
            sequence.append(("E", link))
            i = j
        else:
            # plaquette block: find its end with a few random lattice configurations. The cut is where
            # every register is back and a rotation has been applied. The block is then checked in full.
            conf = rng.choice(ar, size=(48, n_links))
            tb = ((conf[:, :, None] >> np.arange(5)[None, None, :]) & 1).reshape(48, n_data).astype(np.int8)
            tb0 = tb.copy()
            ang = np.zeros(48)
            ident = {q: q for q in range(n_data)}
            j, seen_rz = i, False
            while j < len(gates):
                assert all(q < n_data for q in gates[j][2]), "ancilla used in the magnetic layer"
                run_bits(tb, ang, gates[j:j + 1], ident)
                seen_rz |= gates[j][0] == "rz"
                j += 1
                if seen_rz and np.array_equal(tb, tb0):
                    break
            links = sorted({q // 5 for g in gates[i:j] for q in g[2]})
            assert frozenset(links) in plaqs, f"block touches links {links}, not a plaquette of the torus"
            order = plaqs[frozenset(links)]
            loc = {5 * l + k: 5 * links.index(l) + k for l in links for k in range(5)}
            perm = tuple(links.index(l) for l in order)
            key = ("B", localize(gates[i:j], loc), perm)
            if key not in cache:
                cache[key] = plaquette_block_error(key[1], perm, theta_b, retr)
            err = cache[key]
            assert err < TOL, f"plaquette {order}: phase error {err:.2e}"
            worst_b = max(worst_b, err)
            sequence.append(("B", order))
            i = j

    # structure: `steps` times (one electric block per link, then each plaquette once)
    per_step = n_links + L * L
    assert len(sequence) == steps * per_step, (len(sequence), steps * per_step)
    for s in range(steps):
        chunk = sequence[s * per_step:(s + 1) * per_step]
        assert [x[0] for x in chunk] == ["E"] * n_links + ["B"] * (L * L)
        assert sorted(x[1] for x in chunk[:n_links]) == list(range(n_links))
        assert {frozenset(x[1]) for x in chunk[n_links:]} == set(plaqs)
    return dict(blocks=len(sequence), distinct_blocks=len(cache), electric_max_error=worst_e,
                electric_norm_outside=worst_leak, plaquette_max_error=worst_b)


# ------------------------------------------------------------------------------------------------
# B. the 22-qubit single-plaquette step as a statevector
# ------------------------------------------------------------------------------------------------

def check_single_plaquette_file(path, theta_e, theta_b):
    nq, gates = parse(path)
    assert nq == 22
    rng = np.random.default_rng(3)
    amp = rng.normal(size=(24,) * 4) + 1j * rng.normal(size=(24,) * 4)
    amp /= np.linalg.norm(amp)
    full = np.zeros((32,) * 4, dtype=complex)                # [v0, v1, v2, v3], register values of the links
    full[np.ix_(IDX, IDX, IDX, IDX)] = amp
    # as 20 bit axes: link k, most significant bit first. Move bit b of link k to axis 5k + b.
    order = [5 * k + (4 - b) for k in range(4) for b in range(5)]
    psi = np.zeros((1,) + (2,) * 22, dtype=complex)
    psi[0][..., 0, 0] = np.transpose(full.reshape((2,) * 20), order)
    psi = run_dense(psi, gates, {q: q for q in range(22)})
    out_full = np.transpose(psi[0][..., 0, 0], np.argsort(order)).reshape((32,) * 4)
    out = out_full[np.ix_(IDX, IDX, IDX, IDX)].copy()
    # norm outside the gauge register with the ancillas in |0>
    rest = psi[0].copy()
    rest[..., 0, 0] = 0.0
    out_full[np.ix_(IDX, IDX, IDX, IDX)] = 0.0
    outside = math.hypot(float(np.linalg.norm(rest)), float(np.linalg.norm(out_full)))
    u_e = expm_herm(laplacian(), theta_e)
    ref = amp
    for axis in range(4):
        ref = np.moveaxis(np.tensordot(u_e, ref, axes=([1], [axis])), 0, axis)
    ref = np.exp(1j * theta_b * plaquette_retr()) * ref
    return dict(error=phase_error(out, ref), norm_outside=outside)


# ------------------------------------------------------------------------------------------------

EXPORTS = {"first": ("bt_4x4_first_result_1step.qasm", 0.045, 1),
           "full": ("bt_4x4_full_benchmark_6steps.qasm", 0.025, 6)}
G2, A_S, L_BENCH = 1.0, 0.2, 4                   # the stand-in coupling and the lattice of the exported files


def angles(dt):
    """theta_e = g^2 dt / 2a and theta_b = 2 dt / (g^2 a) from H_KS."""
    return 0.5 * G2 * dt / A_S, 2.0 * dt / (G2 * A_S)


def main():
    repo = os.path.dirname(os.path.dirname(HERE))
    sys.path.insert(0, repo)
    from estimates.common import EPS_SYN, eps_rot_for, t_per_rotation, toffoli_t

    check_group()
    res = {}
    t0 = time.time()
    r = check_single_plaquette_file(os.path.join(HERE, "bt_single_plaquette_1step.qasm"), *angles(0.045))
    assert r["error"] < TOL and r["norm_outside"] < TOL, r
    res["bt_single_plaquette_1step.qasm, statevector"] = r
    for name, dt, steps in EXPORTS.values():
        path = os.path.join(HERE, name)
        res[name] = check_lattice_file(path, L_BENCH, *angles(dt), steps)
        nq, gates = parse(path)
        c = counts(nq, gates)
        eps = eps_rot_for(c["rot"], EPS_SYN)
        c_t = t_per_rotation(eps)
        c.update(eps_rot=eps, c_t=c_t,
                 t_total=toffoli_t(c["toffoli"], "textbook") + c["t"] + c["rot"] * c_t,
                 t_depth=depth(nq, gates, c_t))
        res[name + ", counts"] = c
    for k, v in res.items():
        print(k, v)
    print(f"all checks passed in {time.time() - t0:.0f} s")
    return res


if __name__ == "__main__":
    main()
