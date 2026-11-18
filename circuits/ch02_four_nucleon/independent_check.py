"""Independent check of the exported 36-qubit circuit, written in the verification pass.

    python independent_check.py --count     recount of the two QASM files only (seconds)
    python independent_check.py             recount and exact simulation of the exported (p0, p0) Re circuit
                                            (about 3 min, 1.5 GB, one machine, 4 threads)
    python independent_check.py --pairs     the same for all four current pairs (the other three are generated with
                                            build.py and read back as QASM text) and for the merged variant
                                            (about 7 min)
    python independent_check.py --steps 4   the four pairs of the 4-step circuit, tau = 0.02 (also --steps 2), all
                                            generated with build.py. Writes independent_check_steps4.json.

Writes independent_check.json. Nothing here imports count.py or test_circuit.py, and the circuit side reads the QASM
text with its own parser. build.py is imported only with --pairs or --steps, to write the circuits that are not
exported.

Circuit side. The 36-qubit file is simulated exactly, gate for gate, without a 2^36 vector. Two properties of the
gate list make this possible. Both are checked on the gate list while it runs, and neither is assumed:
  (1) Between its two Hadamards the ancilla is touched only by rz and by cx with the ancilla as the control. The run
      is split into the ancilla = 0 and ancilla = 1 branches, each a circuit on position and workspace qubits.
      Then <Z> of the Hadamard test is Re <phi_0|phi_1> (Im with the sdg), with phi_a the two branch states.
  (2) Every workspace qubit is back in |0> whenever a group of consecutive gates ends. A group is either a run with
      no h, which is a permutation times a diagonal and is evaluated on every basis state of its support, or a run on
      at most DCAP position qubits, which is multiplied out as a dense matrix whose leak into workspace != 0 is checked.
The state kept is that of the 24 position qubits, 2^24 complex doubles.

Reference side. Split-operator evolution on the 3^12 coordinate tensor with the kinetic term applied in momentum space
by FFT (eigenvalue 2t(1 - cos(2 pi k/3)) per coordinate), the contact term from site equality, the ground state by
Lanczos on a nearest-neighbour stencil, and the exact correlator by a Chebyshev expansion of exp(-i H tau). The
couplings, q, tau and the step count are typed here from app01_neutrino_nucleus.tex. The eight preparation angles are
read from the study file, because the chapter does not print them.
"""
from __future__ import annotations

import cmath
import itertools
import json
import math
import re
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
STUDY = REPO / "applications" / "ch02_neutrino_nucleus" / "neutrino_feasibility"
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from estimates import ch02_nu_nucleus as model                                   # noqa: E402
from estimates.common import (EPS_SYN, REACTION_TIME_S, RUS_ANCILLA, depth_exports,   # noqa: E402
                              eps_rot_for, t_per_rotation, toffoli_t)

BASELINE = HERE / "ch02_3x3x3_steps6_re_p0p0.qasm3"
MERGED = HERE / "ch02_3x3x3_steps6_re_p0p0_merged.qasm3"

# app01_neutrino_nucleus.tex, "Concrete instantiations"
T_HOP, C0, D0 = 10.5794, -98.2265511, 127.839693
Q = 2 * math.pi / 3
CHAPTER_ROWS = {2: (0.01, 0.00609), 4: (0.02, 0.01098), 6: (0.03, 0.01434)}    # steps: (tau in MeV^-1, printed bias)
STEPS = 6                                                                      # the exported circuit
TOL = 1e-9

NPOS, NWORK = 24, 11
WORK = set(range(NPOS, NPOS + NWORK))
ANC = NPOS + NWORK
DCAP, PCAP = 4, 19            # largest dense group (position qubits) and largest permutation group (2^19 basis states)


# --------------------------------------------------------------------------- #
# QASM 3 text -> gate list
# --------------------------------------------------------------------------- #

_DECL = re.compile(r"^qubit\[(\d+)\]\s+(\w+);$")
_GATE = re.compile(r"^(\w+)(?:\(([^)]*)\))?\s+(.*);$")
_ARG = re.compile(r"(\w+)\[(\d+)\]")


def parse_qasm3(text):
    """(number of qubits, {register: (offset, size)}, [(name, [flat qubit indices], angle or None)])."""
    regs, off, gates = {}, 0, []
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith(("OPENQASM", "include", "bit[")):
            continue
        m = _DECL.match(line)
        if m:
            regs[m.group(2)] = (off, int(m.group(1)))
            off += int(m.group(1))
            continue
        if "measure" in line:
            qs = [regs[n][0] + int(i) for n, i in _ARG.findall(line.split("measure")[1])]
            gates.append(("measure", qs, None))
            continue
        m = _GATE.match(line)
        if not m:
            raise ValueError("line not understood: " + line)
        name, par, args = m.groups()
        ang = None
        if par is not None:
            if not re.fullmatch(r"[0-9eE+\-*/. pi()]+", par):
                raise ValueError("angle not understood: " + par)
            ang = float(eval(par, {"__builtins__": {}}, {"pi": math.pi}))
        gates.append((name, [regs[n][0] + int(i) for n, i in _ARG.findall(args)], ang))
    return off, regs, gates


# --------------------------------------------------------------------------- #
# Recount
# --------------------------------------------------------------------------- #

def classify(name, ang, tol=1e-9):
    if name == "measure":
        return "skip"
    if name in ("h", "s", "sdg", "x", "z", "cx"):
        return "clifford"
    if name in ("t", "tdg"):
        return "fixed_T"
    if name == "ccx":
        return "toffoli"
    if name == "rz":
        k = ang / (math.pi / 4)
        if abs(k - round(k)) < tol:
            return "fixed_T" if round(k) % 2 else "clifford"
        return "rotations"
    raise ValueError("gate outside the logical set: " + name)


def toffoli_7t(a, b, t):
    """The textbook 7-T Toffoli, typed out here (controls a, b, target t)."""
    return [("h", [t]), ("cx", [b, t]), ("tdg", [t]), ("cx", [a, t]), ("t", [t]), ("cx", [b, t]), ("tdg", [t]),
            ("cx", [a, t]), ("t", [b]), ("t", [t]), ("h", [t]), ("cx", [a, b]), ("t", [a]), ("tdg", [b]), ("cx", [a, b])]


def depth(gates, n, w_rot, toffoli):
    """Longest path counting non-Clifford weight. toffoli = 'layer' (weight 1) or '7t' (the circuit above)."""
    d = [0.0] * n

    def put(qs, w):
        if w == 0.0 and len(qs) == 1:
            return
        new = max(d[x] for x in qs) + w
        for x in qs:
            d[x] = new

    for name, qs, ang in gates:
        kind = classify(name, ang)
        if kind == "skip":
            continue
        if kind == "toffoli" and toffoli == "7t":
            for nm, q2 in toffoli_7t(*qs):
                put(q2, 1.0 if nm in ("t", "tdg") else 0.0)
        else:
            put(qs, {"clifford": 0.0, "fixed_T": 1.0, "toffoli": 1.0, "rotations": w_rot}[kind])
    return max(d)


def recount(gates, n):
    c = dict(rotations=0, toffoli=0, fixed_T=0, clifford=0)
    for name, qs, ang in gates:
        kind = classify(name, ang)
        if kind != "skip":
            c[kind] += 1
    eps = eps_rot_for(c["rotations"], EPS_SYN)
    c_t = t_per_rotation(eps, "rus")
    t_total = c["rotations"] * c_t + toffoli_t(c["toffoli"], "textbook") + c["fixed_T"]
    busy, full, where = set(), 0, {"workspace": 0, "ancilla": 0, "position": 0}
    for name, qs, ang in gates:                         # rotations applied while every workspace qubit is in use
        if name == "ccx" and qs[2] in WORK:
            busy ^= {qs[2]}
        elif classify(name, ang) == "rotations" and len(busy) > NWORK - RUS_ANCILLA:
            full += 1
            where["ancilla" if qs[0] == ANC else "workspace" if qs[0] in WORK else "position"] += 1
    assert not busy
    return dict(
        qubits_declared=n, qubits_touched=len({x for _, qs, _ in gates for x in qs}),
        gates_excluding_measure=sum(1 for g in gates if g[0] != "measure"), gate_names=sorted({g[0] for g in gates}),
        largest_gate_arity=max(len(g[1]) for g in gates), **c,
        eps_rot=eps, t_per_rotation=c_t, t_total_report=t_total,
        depth=dict(non_clifford_layers=depth(gates, n, 1.0, "layer"),
                   t_layers_rz_cT_toffoli_7T=depth(gates, n, c_t, "7t"),
                   t_layers_rz_cT_toffoli_one_layer=depth(gates, n, c_t, "layer")),
        rotations_with_no_idle_workspace=dict(total=full, **where))


# --------------------------------------------------------------------------- #
# Exact simulation of the gate list
# --------------------------------------------------------------------------- #

_H2 = np.array([[1, 1], [1, -1]], dtype=complex) / math.sqrt(2)
_DIAG = {"s": 1j, "sdg": -1j, "t": cmath.exp(1j * math.pi / 4), "tdg": cmath.exp(-1j * math.pi / 4), "z": -1}


def small_simulate(gates, n, state):
    """Gate-by-gate statevector on n qubits, array shape (2,)*n with axis = qubit. Used for the dense groups."""
    for name, qs, ang in gates:
        if name == "cx" or name == "ccx":
            sl = [slice(None)] * n
            for c in qs[:-1]:
                sl[c] = 1
            sub = state[tuple(sl)]
            ax = qs[-1] - sum(qs[-1] > c for c in qs[:-1])
            sub[...] = np.flip(sub, axis=ax).copy()
        elif name == "x":
            state = np.flip(state, axis=qs[0]).copy()
        elif name == "h":
            state = np.moveaxis(np.tensordot(_H2, state, axes=(1, qs[0])), 0, qs[0]).copy()
        else:
            d0, d1 = (cmath.exp(-1j * ang / 2), cmath.exp(1j * ang / 2)) if name == "rz" else (1.0, _DIAG[name])
            sl = [slice(None)] * n
            sl[qs[0]] = 0
            state[tuple(sl)] *= d0
            sl[qs[0]] = 1
            state[tuple(sl)] *= d1
    return state


def split_on_ancilla(gates):
    """(preparation gates, gates between the ancilla's Hadamards, 're' or 'im'). Asserts property (1)."""
    idx = [i for i, g in enumerate(gates) if ANC in g[1]]
    assert gates[idx[0]][0] == "h", "the ancilla must open with h"
    assert gates[-1] == ("measure", [ANC], None) and gates[-2] == ("h", [ANC], None)
    stop, comp = len(gates) - 2, "re"
    if gates[stop - 1] == ("sdg", [ANC], None):
        stop, comp = stop - 1, "im"
    for i in idx[1:]:
        if i >= stop:
            break
        name, qs, _ = gates[i]
        assert (name == "rz" and qs == [ANC]) or (name == "cx" and qs[0] == ANC), f"gate {i} {gates[i]} on the ancilla"
    return gates[:idx[0]], gates[idx[0] + 1:stop], comp


def specialize(mid, a):
    """The gates between the Hadamards for ancilla value a, and the scalar from the ancilla's own rz gates."""
    out, scalar = [], 1.0 + 0j
    for name, qs, ang in mid:
        if ANC not in qs:
            out.append((name, qs, ang))
        elif name == "rz":
            scalar *= cmath.exp((1j if a else -1j) * ang / 2)
        elif a:
            out.append(("x", [qs[1]], None))
    return out, scalar


class Sim:
    def __init__(self, state=None):
        if state is None:
            state = np.zeros((2,) * NPOS, dtype=complex)
            state[(0,) * NPOS] = 1
        self.state = state
        self.dense_groups = self.perm_groups = 0
        self.leak = 0.0

    def run(self, gates):
        seg, supp, has_h = [], set(), False
        for g in gates:
            ns, nh = supp | set(g[1]), has_h or g[0] == "h"
            n_pos = sum(1 for x in ns if x not in WORK)
            fits = (n_pos <= DCAP and len(ns) <= DCAP + 1) if nh else (n_pos <= PCAP)
            if seg and not fits:
                self._flush(seg, supp, has_h)
                seg, supp, has_h = [g], set(g[1]), g[0] == "h"
            else:
                seg.append(g)
                supp, has_h = ns, nh
        if seg:
            self._flush(seg, supp, has_h)
        return self

    def _flush(self, seg, supp, has_h):
        (self._dense if has_h else self._perm_diag)(seg, sorted(supp))

    def _dense(self, seg, L):
        k = len(L)
        loc = {x: j for j, x in enumerate(L)}
        local = [(n, [loc[x] for x in qs], a) for n, qs, a in seg]
        U = np.empty((2 ** k, 2 ** k), dtype=complex)
        for col in range(2 ** k):
            e = np.zeros(2 ** k, dtype=complex)
            e[col] = 1
            U[:, col] = small_simulate(local, k, e.reshape((2,) * k)).reshape(-1)
        ids = np.arange(2 ** k)
        clean = np.ones(2 ** k, dtype=bool)
        for j in range(k):
            if L[j] in WORK:
                clean &= ((ids >> (k - 1 - j)) & 1) == 0
        good = ids[clean]
        if not clean.all():
            self.leak = max(self.leak, float(np.max(np.abs(U[np.ix_(ids[~clean], good)]))))
            assert self.leak < 1e-12, "a dense group does not return its workspace"
        axes = [x for x in L if x not in WORK]
        m = len(axes)
        Ur = U[np.ix_(good, good)].reshape((2,) * (2 * m))
        self.state = np.moveaxis(np.tensordot(Ur, self.state, axes=(list(range(m, 2 * m)), axes)), list(range(m)), axes)
        self.dense_groups += 1

    def _perm_diag(self, seg, L):
        P = [x for x in L if x not in WORK]
        k, N = len(P), 2 ** len(P)
        ids = np.arange(N)
        bits = {x: ((ids >> (k - 1 - j)) & 1).astype(bool) for j, x in enumerate(P)}
        for x in L:
            if x in WORK:
                bits[x] = np.zeros(N, dtype=bool)
        phase = np.ones(N, dtype=complex)
        for name, qs, ang in seg:
            if name == "x":
                bits[qs[0]] = ~bits[qs[0]]
            elif name == "cx":
                bits[qs[1]] = bits[qs[1]] ^ bits[qs[0]]
            elif name == "ccx":
                bits[qs[2]] = bits[qs[2]] ^ (bits[qs[0]] & bits[qs[1]])
            elif name == "rz":
                phase *= np.where(bits[qs[0]], cmath.exp(1j * ang / 2), cmath.exp(-1j * ang / 2))
            else:
                phase[bits[qs[0]]] *= _DIAG[name]
        for x in L:
            if x in WORK:
                assert not bits[x].any(), "a permutation group does not return its workspace"
        out = np.zeros(N, dtype=np.int64)
        for j, x in enumerate(P):
            out |= bits[x].astype(np.int64) << (k - 1 - j)
        if np.array_equal(out, ids):
            shape = [1] * NPOS
            for x in P:
                shape[x] = 2
            self.state *= phase.reshape(shape)                 # P is sorted, so the axes line up
        else:
            assert len(np.unique(out)) == N
            st = np.moveaxis(self.state, P, list(range(k))).reshape(N, -1)
            new = np.empty_like(st)
            new[out] = phase[:, None] * st
            self.state = np.moveaxis(new.reshape((2,) * NPOS), list(range(k)), P)
        self.perm_groups += 1


def to_coordinates(state):
    """(2,)*24 with axis = qubit -> (3,)*12 with x = bit(2c) + 2 bit(2c+1), and the weight outside x in {0, 1, 2}."""
    s4 = np.ascontiguousarray(state).reshape((4,) * 12)        # index of one coordinate = 2 bit(2c) + bit(2c+1)
    legal = s4[np.ix_(*[[0, 2, 1]] * 12)]
    return legal, float(np.sum(np.abs(s4) ** 2) - np.sum(np.abs(legal) ** 2))


# --------------------------------------------------------------------------- #
# Reference
# --------------------------------------------------------------------------- #

class Lattice:
    """Four nucleons on 3^3. Tensor axis c = 3 * particle + axis."""

    def __init__(self):
        self.shape = (3,) * 12
        idx = np.indices(self.shape)
        self.x = idx
        site = [idx[3 * f] + 3 * idx[3 * f + 1] + 9 * idx[3 * f + 2] for f in range(4)]
        V = np.zeros(self.shape)
        for a, b in itertools.combinations(range(4), 2):
            V += C0 * (site[a] == site[b])
        for a, b, c in itertools.combinations(range(4), 3):
            V += D0 * ((site[a] == site[b]) & (site[b] == site[c]))
        self.V = V
        self.ek = sum(2 * T_HOP * (1 - np.cos(2 * math.pi * idx[c] / 3)) for c in range(12))
        self.const = 24 * T_HOP

    def kin(self, psi, dt, hop_only=False):
        e = self.ek - self.const if hop_only else self.ek
        return np.fft.ifftn(np.exp(-1j * dt * e) * np.fft.fftn(psi))

    def pot(self, psi, dt):
        return np.exp(-1j * dt * self.V) * psi

    def H(self, psi):
        out = (self.const + self.V) * psi
        for c in range(12):
            out = out - T_HOP * (np.roll(psi, 1, c) + np.roll(psi, -1, c))
        return out

    def prepare(self, layers):
        psi = np.full(self.shape, 3.0 ** -6, dtype=complex)
        for dt_v, dt_k in layers:
            psi = self.kin(self.pot(psi, dt_v), dt_k, hop_only=True)
        return psi

    def product_formula(self, psi, tau, n):
        dt = tau / n
        psi = self.kin(psi, dt / 2)
        for j in range(n):
            psi = self.kin(self.pot(psi, dt), dt if j < n - 1 else dt / 2)
        return psi

    def current(self, psi, f, sign=1):
        return np.exp(1j * sign * Q * self.x[3 * f]) * psi

    def exact(self, psi, tau):
        """exp(-i H tau) psi by Chebyshev expansion. The spectrum of H lies inside [lo, hi]."""
        lo, hi = -6 * abs(C0) - 10, 36 * T_HOP + 4 * D0 + 10
        a, b = (hi - lo) / 2, (hi + lo) / 2
        z = a * tau
        from scipy.special import jv
        t0, t1 = psi, (self.H(psi) - b * psi) / a
        out = jv(0, z) * t0 - 2j * jv(1, z) * t1
        for k in range(2, int(z + 12 * z ** (1 / 3) + 30)):
            t0, t1 = t1, 2 * (self.H(t1) - b * t1) / a - t0
            out = out + 2 * (-1j) ** k * jv(k, z) * t1
        return np.exp(-1j * b * tau) * out


def reference(steps):
    from scipy.sparse.linalg import LinearOperator, eigsh
    TAU, STEPS = CHAPTER_ROWS[steps][0], steps
    lat = Lattice()
    prep = json.loads((STUDY / "unitary_prep_results.json").read_text())[-1]
    assert prep["layers"] == 8
    th = prep["theta"]
    psi = lat.prepare([(th[2 * j] / 100, th[2 * j + 1] / 100) for j in range(8)])
    n = 3 ** 12
    op = LinearOperator((n, n), matvec=lambda v: lat.H(v.reshape(lat.shape)).ravel(), dtype=float)
    e, vec = eigsh(op, k=1, which="SA", tol=1e-11, v0=np.ones(n))
    e0, gs = float(e[0]), vec[:, 0].reshape(lat.shape).astype(complex)
    phase = cmath.exp(1j * e0 * TAU)
    w, c_pf, c_exact_prep, c_exact = {}, {}, 0, 0
    for b in (0, 1):
        pf = lat.product_formula(lat.current(psi, b), TAU, STEPS)
        ex_p = lat.exact(lat.current(psi, b), TAU)
        ex_g = lat.exact(lat.current(gs, b), TAU)
        for a in (0, 1):
            w[a, b] = phase * lat.current(pf, a, -1)
            c_pf[a, b] = np.vdot(psi, w[a, b])
            c_exact_prep += phase * np.vdot(lat.current(psi, a), ex_p) / 4
            c_exact += phase * np.vdot(lat.current(gs, a), ex_g) / 4
    j_gs = (lat.current(gs, 0) + lat.current(gs, 1)) / 2
    c_avg = sum(c_pf.values()) / 4
    row = next(r for r in json.loads((STUDY / "larger_unitary_response.json").read_text())["rows"]
               if abs(r["tau_MeV_inv"] - TAU) < 1e-12 and r["steps"] == STEPS)
    study_e0 = json.loads((STUDY / "larger_box_results.json").read_text())["ground_energy_MeV"]
    numbers = dict(
        E0_MeV=e0, E0_minus_study_file=e0 - study_e0,
        prepared_energy_MeV=float(np.vdot(psi, lat.H(psi)).real), overlap_squared=float(abs(np.vdot(gs, psi)) ** 2),
        static_strength=float(np.vdot(j_gs, j_gs).real),
        C_product_formula_prepared=[c_avg.real, c_avg.imag],
        C_exact_evolution_prepared=[c_exact_prep.real, c_exact_prep.imag],
        C_exact_ground_state=[c_exact.real, c_exact.imag],
        bias_total=abs(c_avg - c_exact), bias_over_abs_C=abs(c_avg - c_exact) / abs(c_exact),
        bias_preparation_only=abs(c_exact_prep - c_exact), bias_product_formula_only=abs(c_avg - c_exact_prep),
        study_file_circuit=row["circuit"], study_file_exact=row["exact"],
        tau_MeV_inv=TAU, steps=STEPS,
        chapter_prints=dict(E0=-52.67220, prepared_energy=-52.45933, overlap_squared=0.997943,
                            bias=CHAPTER_ROWS[steps][1], static_strength=0.65942,
                            **(dict(bias_percent=2.35, preparation_only=0.00414) if steps == 6 else {})))
    return lat, psi, gs, w, c_pf, c_exact, numbers


# --------------------------------------------------------------------------- #
# Driver
# --------------------------------------------------------------------------- #

def load(path):
    n, regs, gates = parse_qasm3(Path(path).read_text())
    assert n == 36 and regs == {"pos": (0, NPOS), "work": (NPOS, NWORK), "anc": (ANC, 1)}, regs
    return n, gates


def count_report():
    a = model.Assumptions()
    res = model.model(a, "2028")
    out = {}
    for name, path in (("baseline", BASELINE), ("merged", MERGED)):
        n, gates = load(path)
        out[name] = dict(file=path.name, **recount(gates, n))
    b = out["baseline"]
    m = dict(zip(("rotations", "toffoli", "fixed_T"), model.counts("2028", STEPS)))
    m.update(t_total=model.circuit(a, "2028", STEPS)[0], lq=res.lq[0], shots=res.shots[0], wall_campaign_s=res.wall_time_s[0],
             depth_status=res.intermediates["depth_status"],
             source="estimates/ch02_nu_nucleus.py: counts('2028', 6), circuit(Assumptions(), '2028', 6), model(a, '2028')")
    out["model"] = m
    out["circuit_over_model"] = {k: b[k] / m[k] for k in ("rotations", "toffoli", "fixed_T")} | {
        "t_total": b["t_total_report"] / m["t_total"], "qubits": b["qubits_touched"] / m["lq"]}
    walls = {}
    for key, d in b["depth"].items():
        if key == "non_clifford_layers":
            continue
        dx = depth_exports(b["t_total_report"], d, res.shots[0], a.t_gate_s.lo, a.shot_overhead_s.lo)
        walls[key] = dict(t_depth=d, f_star=dx["f_star"][0], floor_wall_s=dx["floor_wall_s"][0],
                          floor_over_model_wall=dx["floor_wall_s"][0] / res.wall_time_s[0])
    out["depth_floor_every_shot_at_the_deepest_circuit"] = dict(
        source=f"estimates.common.depth_exports, reaction time {REACTION_TIME_S} s per layer, one machine", **walls)
    return out


def simulate_file(path, lat, psi, gs, w, c_pf, cache, t_start, steps):
    n, gates = load(path)
    pre, mid, comp = split_on_ancilla(gates)
    if "prep" not in cache:
        s = Sim().run(pre)
        cache.update(prep=np.ascontiguousarray(s.state), pre=pre, prep_groups=(s.dense_groups, s.perm_groups))
        print(f"  preparation: {len(pre)} gates in {s.dense_groups} dense and {s.perm_groups} permutation groups "
              f"[{time.time() - t_start:.0f} s]", flush=True)
    assert pre == cache["pre"], "the preparation differs between files"
    branch = []
    for a in (0, 1):
        g, scalar = specialize(mid, a)
        if a == 0 and cache.get("mid0") == g:
            state = cache["phi0"]                       # the ancilla = 0 branch is the same gate list in every file
        else:
            s = Sim(cache["prep"].copy()).run(g)
            state = s.state
            if a == 0:
                cache.update(mid0=g, phi0=np.ascontiguousarray(state))
            print(f"  ancilla = {a}: {len(g)} gates in {s.dense_groups} dense and {s.perm_groups} permutation groups "
                  f"[{time.time() - t_start:.0f} s]", flush=True)
        coords, leak = to_coordinates(state)
        branch.append((coords * scalar, leak))
    (phi0, leak0), (phi1, leak1) = branch
    prep_c, leak_p = to_coordinates(cache["prep"])
    g = np.vdot(psi, prep_c)
    g /= abs(g)                                         # global phase of the QASM text, not observable
    g0 = np.vdot(psi, phi0)
    g0 /= abs(g0)
    err_w = {k: float(np.max(np.abs(phi1 - g0 * v))) for k, v in w.items()}
    pair = min(err_w, key=err_w.get)
    c = np.vdot(phi0, phi1)
    counted = recount(gates, n)
    row = dict(
        file=Path(path).name, component=comp, current_pair=list(pair),
        counts=[counted["rotations"], counted["toffoli"], counted["fixed_T"]],
        model_counts=list(model.counts("2028", steps)), t_total_report=counted["t_total_report"],
        model_t_total=model.circuit(model.Assumptions(), "2028", steps)[0], qubits_touched=counted["qubits_touched"],
        depth=counted["depth"],
        weight_outside_legal_coordinates=max(leak0, leak1, leak_p),
        prepared_state_max_amplitude_error=float(np.max(np.abs(prep_c - g * psi))),
        prepared_energy_from_circuit_MeV=float(np.vdot(prep_c, lat.H(prep_c)).real),
        overlap_squared_from_circuit=float(abs(np.vdot(gs, prep_c)) ** 2),
        ancilla0_branch_max_amplitude_error=float(np.max(np.abs(phi0 - g0 * psi))),
        ancilla1_branch_max_amplitude_error=err_w[pair],
        ancilla1_branch_error_against_nearest_other_pair=min(v for k, v in err_w.items() if k != pair),
        C_ab_circuit=[c.real, c.imag], C_ab_reference=[c_pf[pair].real, c_pf[pair].imag],
        C_ab_error=abs(c - c_pf[pair]), ancilla_expectation=float(c.real if comp == "re" else c.imag))
    ok = max(row["prepared_state_max_amplitude_error"], row["ancilla0_branch_max_amplitude_error"],
             row["ancilla1_branch_max_amplitude_error"], row["C_ab_error"], row["weight_outside_legal_coordinates"]) < TOL
    row["passes"] = bool(ok)
    print(f"  {row['file']}: {comp} of pair {pair}, branch errors {row['ancilla0_branch_max_amplitude_error']:.1e} / "
          f"{row['ancilla1_branch_max_amplitude_error']:.1e}, C_ab {c:.9f}, error {row['C_ab_error']:.1e}, "
          f"{'pass' if ok else 'FAIL'}", flush=True)
    return row, pair, c


def main(argv):
    t_start = time.time()
    out = dict(tolerance=TOL, machine="one machine, 4 threads", counts=count_report())
    b, m = out["counts"]["baseline"], out["counts"]["model"]
    print(f"recount of {BASELINE.name}: {b['qubits_touched']} qubits, {b['rotations']} rz, {b['toffoli']} ccx, "
          f"{b['fixed_T']} T, T total {b['t_total_report']:.1f} (model {m['t_total']:.1f}), "
          f"largest gate on {b['largest_gate_arity']} qubits")
    print(f"  depth {b['depth']}")
    print(f"  rotations applied with no idle workspace: {b['rotations_with_no_idle_workspace']}")
    steps = int(argv[argv.index("--steps") + 1]) if "--steps" in argv else STEPS
    assert steps in CHAPTER_ROWS
    if "--count" not in argv:
        lat, psi, gs, w, c_pf, c_exact, numbers = reference(steps)
        out["reference"] = numbers
        print(f"reference: E0 {numbers['E0_MeV']:.5f}, prepared energy {numbers['prepared_energy_MeV']:.5f}, "
              f"overlap^2 {numbers['overlap_squared']:.6f}, bias {numbers['bias_total']:.5f} "
              f"({100 * numbers['bias_over_abs_C']:.2f}% of |C|) [{time.time() - t_start:.0f} s]", flush=True)
        files = [BASELINE] if steps == STEPS else []
        scratch = []
        if "--pairs" in argv or steps != STEPS:
            import tempfile
            sys.path.insert(0, str(HERE))
            import build
            tmp = Path(tempfile.mkdtemp(prefix="ch02_pairs_"))
            todo = [((0, 1), "re"), ((1, 0), "im"), ((1, 1), "im")]
            if steps != STEPS:
                todo.insert(0, ((0, 0), "re"))
            for pair, comp in todo:
                p = tmp / f"generated_steps{steps}_{comp}_p{pair[0]}p{pair[1]}.qasm3"
                inst = build.benchmark_instance(steps, current_pair=pair, component=comp)
                assert abs(inst.tau - CHAPTER_ROWS[steps][0]) < 1e-15
                build.export_qasm3(build.build_circuit(inst), p)
                files.append(p)
                scratch.append(p)
            if steps == STEPS:
                files.append(MERGED)
        cache, rows, by_pair = {}, [], {}
        for path in files:
            row, pair, c = simulate_file(path, lat, psi, gs, w, c_pf, cache, t_start, steps)
            rows.append(row)
            by_pair[pair] = c
        for p in scratch:
            p.unlink()
        if scratch:
            scratch[0].parent.rmdir()
        out["whole_instance"] = dict(method="ancilla branches, grouped exact simulation on 2^24 amplitudes", files=rows,
                                     all_pass=all(r["passes"] for r in rows))
        if len(by_pair) == 4:
            c = sum(by_pair.values()) / 4
            out["whole_instance"]["four_pair_average_from_circuits"] = dict(
                C=[c.real, c.imag], bias_against_exact=abs(c - c_exact),
                difference_from_study_file=abs(c - complex(*numbers["study_file_circuit"])))
            print(f"four-pair average from the circuits {c:.6f}, bias against the exact value {abs(c - c_exact):.5f}")
        out["seconds"] = round(time.time() - t_start)
        name = "independent_check.json" if steps == STEPS else f"independent_check_steps{steps}.json"
        (HERE / name).write_text(json.dumps(out, indent=1) + "\n")
        print(f"wrote {name}, {out['seconds']} s")
        return 0 if out["whole_instance"]["all_pass"] else 1
    print(json.dumps(out["counts"]["circuit_over_model"]))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
