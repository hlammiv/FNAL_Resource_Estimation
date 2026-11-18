"""Gate-level circuit for the Ch. 4 2028 benchmark: log det W of free staggered fermions by a
Hadamard test on a QET (QSVT) of the block encoding of arXiv:2407.13080 Sec. II.2.

Everything is built from the instance parameters, so other sizes can be generated:

    inst = Instance(L=4, D=4, m0=0.4, K=1.0, variant="literal")
    qc = hadamard_test_circuit(inst, reflection_phases)

Gate set of every circuit returned here: h, x, s, sdg, cx, ccx, rz (and measure when asked for).
Arbitrary rotations are explicit rz gates with their angles. No gate is hidden inside a box.

Blocks (names used in the block marks, see `block_marks`):
    D_s            Hadamards on the column register, plus the X that complements the top column bit
    O_A            entry oracle, y-rotations on the flag qubit controlled on the top column bit
    O_c            shift oracle, controlled add2/sub2 on the coordinate registers (Toffoli ladders)
    proj_and       compute and uncompute of the AND of the block-encoding qubits onto the signal qubit
    proj_rot       the QET phase on the signal qubit, controlled on the Hadamard-test ancilla
    had            the two Hadamards of the Hadamard test, the folded last phase if any
    prep           X gates that prepare the site basis state

Variants (the model's two ends, and one further compile that the model does not contain):
    literal   O_c as Fig. 1 draws it: 2D controlled shifts, each controlled on all column bits.
              O_A as Fig. 2 draws it: two singly-controlled R_y, each as 2 rotations and 2 CNOT.
    merged    L = 4 only, where +2 and -2 are the same shift: D shifts without the sign control.
              O_A as one multiplexed rotation (2 rotations and 2 CNOT).
    compact   merged, with the O_c controls shared between directions (unary iteration) and the last
              QET phase applied as one rz on the Hadamard-test ancilla. D in {2, 4}, L = 4.
"""
from __future__ import annotations

import importlib.util
import json
import math
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from qiskit import ClassicalRegister, QuantumCircuit, QuantumRegister

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
POLY_FILE = HERE / "polynomial_d84.json"
VARIANTS = ("literal", "merged", "compact")
ZERO_ANGLE = 1e-12          # rotations below this angle are the identity and are not emitted
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))


# --------------------------------------------------------------------------- #
# Instance
# --------------------------------------------------------------------------- #

@dataclass(frozen=True)
class Instance:
    """Free staggered W = M^dag M on a periodic L^D lattice, block-encoded as W/s."""
    L: int = 4
    D: int = 4
    m0: float = 0.4
    K: float = 1.0
    s: float | None = None          # subnormalization; None means the smallest allowed value
    variant: str = "literal"

    def __post_init__(self):
        if self.L < 4 or self.L & (self.L - 1):
            raise ValueError("L must be a power of two >= 4 (on L = 2 every two-hop shift is the identity)")
        if self.D < 1 or self.D & (self.D - 1):
            raise ValueError("D must be a power of two (the column register is top bit, sign bit, direction bits)")
        if self.variant not in VARIANTS:
            raise ValueError(f"variant must be one of {VARIANTS}")
        if self.variant in ("merged", "compact") and self.L != 4:
            raise ValueError("merged and compact need L = 4, where the +2 and -2 shifts coincide")
        if self.variant == "compact" and self.D not in (2, 4):
            raise ValueError("compact is written for D = 2 or 4")
        if self.s is not None and self.s < self.s_min * (1 - 1e-12):
            raise ValueError(f"s = {self.s} is below the smallest allowed value {self.s_min}")

    # --- sizes -----------------------------------------------------------
    @property
    def bits(self) -> int:                 # qubits per coordinate
        return self.L.bit_length() - 1

    @property
    def n_dir(self) -> int:                # column bits that hold the direction
        return self.D.bit_length() - 1

    @property
    def n_col(self) -> int:                # direction bits, sign bit, top bit (0: hop, 1: diagonal)
        return self.n_dir + 2

    @property
    def n_site(self) -> int:
        return self.D * self.bits

    @property
    def V(self) -> int:
        return self.L ** self.D

    @property
    def shift_controls(self) -> int:
        """Largest number of controls on one X of the shift oracle."""
        if self.variant == "literal":
            return self.n_col + self.bits - 2
        return self.n_col - 1

    @property
    def n_work(self) -> int:
        """Clean workspace: n - 2 for an n-control ladder, n_col - 1 for the projector AND."""
        return max(self.shift_controls - 2, self.n_col + 1 - 2, 0)

    @property
    def n_qubits(self) -> int:
        return self.n_site + self.n_col + 1 + 1 + self.n_work + 1

    # --- matrix elements and angles (arXiv:2407.13080 Eqs. 24-26) -----------
    @property
    def diag(self) -> float:
        return self.m0 ** 2 + self.D * self.K ** 2 / 2.0

    @property
    def hop(self) -> float:
        return -self.K ** 2 / 4.0

    @property
    def s_min(self) -> float:
        return max(2.0 * self.diag, 2 ** self.n_col * self.K ** 2 / 4.0)

    @property
    def sub(self) -> float:
        return self.s_min if self.s is None else self.s

    @property
    def thetas(self) -> tuple[float, float]:
        """(theta_0, theta_1): half of the 2^c column slots carry the diagonal, one slot carries each hop."""
        c0 = min(1.0, 2.0 * self.diag / self.sub)
        c1 = max(-1.0, 2 ** self.n_col * self.hop / self.sub)
        return 2.0 * math.acos(c0), 2.0 * math.acos(c1)


BENCHMARK = dict(L=4, D=4, m0=0.4, K=1.0)      # the 2028 benchmark instance, app12 2028 box
BENCHMARK_DEGREE = 84
BENCHMARK_SHIFT = 0.5


# --------------------------------------------------------------------------- #
# Registers and block marks
# --------------------------------------------------------------------------- #

def registers(inst: Instance) -> dict:
    """Qubit order (qiskit little-endian, qubit 0 first): site, col, flag, sig, work, had.

    site: coordinate mu holds bits [mu*bits, (mu+1)*bits), low bit first, so the site index is
    sum_mu x_mu L^mu. col: direction bits (low first), then the sign bit, then the top bit."""
    return {
        "site": QuantumRegister(inst.n_site, "site"),
        "col": QuantumRegister(inst.n_col, "col"),
        "flag": QuantumRegister(1, "flag"),
        "sig": QuantumRegister(1, "sig"),
        "work": QuantumRegister(inst.n_work, "work"),
        "had": QuantumRegister(1, "had"),
    }


def _new_circuit(inst: Instance, name: str) -> tuple[QuantumCircuit, dict]:
    reg = registers(inst)
    qc = QuantumCircuit(*[r for r in reg.values() if len(r)], name=name)
    qc.metadata = {"blocks": []}
    return qc, reg


class _mark:
    """Context manager that records [name, start, stop) instruction indices in qc.metadata['blocks']."""

    def __init__(self, qc: QuantumCircuit, name: str):
        self.qc, self.name = qc, name

    def __enter__(self):
        self.start = len(self.qc.data)

    def __exit__(self, *exc):
        stop = len(self.qc.data)
        if stop > self.start:
            self.qc.metadata["blocks"].append([self.name, self.start, stop])


def block_marks(qc: QuantumCircuit) -> list:
    """[[name, start, stop), ...] covering every instruction of a circuit built here."""
    return qc.metadata["blocks"]


# --------------------------------------------------------------------------- #
# Primitives
# --------------------------------------------------------------------------- #

def ry(qc: QuantumCircuit, theta: float, q) -> None:
    """R_y(theta) = S H R_z(theta) H S^dag, written with an explicit rz."""
    if abs(theta) < ZERO_ANGLE:
        return
    qc.sdg(q)
    qc.h(q)
    qc.rz(theta, q)
    qc.h(q)
    qc.s(q)


def mcx_ladder(qc: QuantumCircuit, controls: list, target, work: list) -> None:
    """C^n X with n - 2 clean ancillas returned clean: 2n - 3 Toffolis for n >= 3."""
    n = len(controls)
    if n == 0:
        qc.x(target)
    elif n == 1:
        qc.cx(controls[0], target)
    elif n == 2:
        qc.ccx(controls[0], controls[1], target)
    else:
        if len(work) < n - 2:
            raise ValueError("not enough workspace")
        steps = [(controls[0], controls[1], work[0])]
        steps += [(work[i - 2], controls[i], work[i - 1]) for i in range(2, n - 1)]
        for a, b, t in steps:
            qc.ccx(a, b, t)
        qc.ccx(work[n - 3], controls[n - 1], target)
        for a, b, t in reversed(steps):
            qc.ccx(a, b, t)


def and_tree_steps(controls: list, target, work: list) -> list:
    """Toffoli list that writes the AND of n >= 2 controls onto `target` with n - 2 ancillas and
    n - 1 Toffolis, pairing inputs as a tree. The ancillas are left holding partial ANDs. Running
    the list backwards uncomputes. Five controls give four Toffolis in three layers."""
    items = list(controls)
    if len(items) < 2:
        raise ValueError("need at least two controls")
    if len(work) < len(items) - 2:
        raise ValueError("not enough workspace")
    steps, k = [], 0
    while len(items) > 2:
        a, b = items.pop(0), items.pop(0)
        steps.append((a, b, work[k]))
        items.append(work[k])
        k += 1
    steps.append((items[0], items[1], target))
    return steps


# --------------------------------------------------------------------------- #
# The block encoding
# --------------------------------------------------------------------------- #
# Inside one query the top column bit is complemented (X after the opening Hadamards, X before the
# closing ones), so that "hop" is top = 1. Then O_A and O_c both use the top bit only as a positive
# control, they commute, and O_c can run while the O_A rotations are in flight.

def _coordinate(reg: dict, inst: Instance, mu: int) -> list:
    return [reg["site"][mu * inst.bits + j] for j in range(inst.bits)]


def _controlled_add2(qc, inst, controls, coord, work, subtract: bool) -> None:
    """Controlled x -> x +- 2 mod L: an increment (decrement) of the coordinate without its low bit,
    as a cascade of multi-controlled X. On L = 4 it is one multi-controlled X on the high bit."""
    upper = coord[1:]
    order = range(len(upper)) if subtract else reversed(range(len(upper)))
    for j in order:
        mcx_ladder(qc, list(controls) + upper[:j], upper[j], work)


def _set_pattern(qc, qubits, current: list, wanted: list) -> list:
    """X the qubits whose complement state must change so that controls on them test `wanted`."""
    for q, c, w in zip(qubits, current, wanted):
        if c != w:
            qc.x(q)
    return list(wanted)


def append_shift_oracle(qc: QuantumCircuit, reg: dict, inst: Instance) -> None:
    """O_c in the complemented frame (top = 1 means hop). Self-contained: every X it adds for a
    negative control is removed again."""
    col, work = list(reg["col"]), list(reg["work"])
    dirs, sign, top = col[:inst.n_dir], col[inst.n_dir], col[inst.n_dir + 1]
    with _mark(qc, "O_c"):
        if inst.variant == "literal":
            pat_q = dirs + [sign]
            cur = [1] * len(pat_q)
            for sg in (0, 1):
                for mu in range(inst.D):
                    want = [(mu >> j) & 1 for j in range(inst.n_dir)] + [sg]
                    cur = _set_pattern(qc, pat_q, cur, want)
                    _controlled_add2(qc, inst, [top] + pat_q, _coordinate(reg, inst, mu), work, bool(sg))
            _set_pattern(qc, pat_q, cur, [1] * len(pat_q))
        elif inst.variant == "merged":
            cur = [1] * inst.n_dir
            for mu in range(inst.D):
                want = [(mu >> j) & 1 for j in range(inst.n_dir)]
                cur = _set_pattern(qc, dirs, cur, want)
                _controlled_add2(qc, inst, [top] + dirs, _coordinate(reg, inst, mu), work, False)
            _set_pattern(qc, dirs, cur, [1] * inst.n_dir)
        else:                                                   # compact, L = 4
            hi = [_coordinate(reg, inst, mu)[1] for mu in range(inst.D)]
            if inst.D == 2:
                qc.ccx(top, dirs[0], hi[1])
                qc.x(dirs[0])
                qc.ccx(top, dirs[0], hi[0])
                qc.x(dirs[0])
            else:                                               # D = 4: 6 Toffolis
                w = work[0]
                qc.ccx(top, dirs[1], w)                         # w = top and d1
                qc.ccx(w, dirs[0], hi[3])
                qc.x(dirs[0])
                qc.ccx(w, dirs[0], hi[2])
                qc.cx(top, w)                                   # w = top and not d1
                qc.ccx(w, dirs[0], hi[0])
                qc.x(dirs[0])
                qc.ccx(w, dirs[0], hi[1])
                qc.cx(top, w)
                qc.ccx(top, dirs[1], w)


def append_entry_oracle(qc: QuantumCircuit, reg: dict, inst: Instance, theta1_override: float | None = None) -> None:
    """O_A in the complemented frame: R_y(theta_1) on the flag when top = 1 (hop), R_y(theta_0) when
    top = 0 (diagonal). `theta1_override` exists so a test can build a deliberately wrong oracle."""
    th0, th1 = inst.thetas
    if theta1_override is not None:
        th1 = theta1_override
    top, f = reg["col"][inst.n_dir + 1], reg["flag"][0]
    with _mark(qc, "O_A"):
        if inst.variant == "literal":
            if abs(th1) >= ZERO_ANGLE:                          # R_y(theta_1) controlled on top = 1
                ry(qc, th1 / 2.0, f)
                qc.cx(top, f)
                ry(qc, -th1 / 2.0, f)
                qc.cx(top, f)
            if abs(th0) >= ZERO_ANGLE:                          # R_y(theta_0) controlled on top = 0
                ry(qc, th0 / 2.0, f)
                qc.cx(top, f)
                ry(qc, th0 / 2.0, f)
                qc.cx(top, f)
        else:                                                   # multiplexed: top = 0 -> a + b, top = 1 -> a - b
            a, b = (th0 + th1) / 2.0, (th0 - th1) / 2.0
            if abs(b) >= ZERO_ANGLE:
                qc.cx(top, f)
                ry(qc, b, f)
                qc.cx(top, f)
            ry(qc, a, f)


def append_block_encoding(qc: QuantumCircuit, reg: dict, inst: Instance, inverse: bool = False,
                          theta1_override: float | None = None) -> None:
    """U = D_s O_c O_A D_s (operator order), <0 0^c| U |0 0^c> = W/s on the site register.
    `inverse=True` appends U^dag with the same block marks."""
    sub, sreg = _new_circuit(inst, "U")
    col = list(sreg["col"])
    top = col[inst.n_dir + 1]
    with _mark(sub, "D_s"):
        sub.h(col)
        sub.x(top)
    append_entry_oracle(sub, sreg, inst, theta1_override)
    append_shift_oracle(sub, sreg, inst)
    with _mark(sub, "D_s"):
        sub.x(top)
        sub.h(col)
    _compose_marked(qc, sub, inverse)


def _compose_marked(qc: QuantumCircuit, sub: QuantumCircuit, inverse: bool = False) -> None:
    """Append `sub` (or its inverse) to `qc` gate by gate, carrying the block marks over."""
    base = len(qc.data)
    n = len(sub.data)
    marks = sub.metadata["blocks"]
    if inverse:
        src = sub.inverse()
        marks = [[name, n - stop, n - start] for name, start, stop in reversed(marks)]
    else:
        src = sub
    qmap = {q: qc.qubits[i] for i, q in enumerate(sub.qubits)}
    for inst_ in src.data:
        qc.append(inst_.operation, [qmap[q] for q in inst_.qubits])
    qc.metadata["blocks"] += [[name, base + start, base + stop] for name, start, stop in marks]


def block_encoding_circuit(inst: Instance, inverse: bool = False, theta1_override: float | None = None) -> QuantumCircuit:
    qc, reg = _new_circuit(inst, "U_dag" if inverse else "U")
    append_block_encoding(qc, reg, inst, inverse, theta1_override)
    return qc


def shift_oracle_circuit(inst: Instance) -> QuantumCircuit:
    """O_c alone, in the complemented frame (top column bit = 1 means hop)."""
    qc, reg = _new_circuit(inst, "O_c")
    append_shift_oracle(qc, reg, inst)
    return qc


def entry_oracle_circuit(inst: Instance) -> QuantumCircuit:
    """O_A alone, in the complemented frame."""
    qc, reg = _new_circuit(inst, "O_A")
    append_entry_oracle(qc, reg, inst)
    return qc


# --------------------------------------------------------------------------- #
# Projector-controlled phase
# --------------------------------------------------------------------------- #

def append_projector_phase(qc: QuantumCircuit, reg: dict, inst: Instance, phi: float, controlled: bool = True) -> None:
    """exp(i phi (2 Pi - 1)), Pi the projector on flag = 0 and column = 0.

    The AND of the n_col + 1 block-encoding qubits is written onto the signal qubit by a tree of
    n_col Toffolis, the phase is an rz on the signal qubit, and the tree is run backwards.
    controlled=True applies the phase only when the Hadamard-test ancilla is 1 (2 rz and 2 CNOT),
    controlled=False always (1 rz)."""
    be = list(reg["col"]) + [reg["flag"][0]]
    sig, had = reg["sig"][0], reg["had"][0]
    steps = and_tree_steps(be, sig, list(reg["work"]))
    with _mark(qc, "proj_and"):
        qc.x(be)
        for a, b, t in steps:
            qc.ccx(a, b, t)
    with _mark(qc, "proj_rot"):
        if controlled:
            qc.rz(phi, sig)
            qc.cx(had, sig)
            qc.rz(-phi, sig)
            qc.cx(had, sig)
        else:
            qc.rz(2.0 * phi, sig)
    with _mark(qc, "proj_and"):
        for a, b, t in reversed(steps):
            qc.ccx(a, b, t)
        qc.x(be)


def projector_phase_circuit(inst: Instance, phi: float, controlled: bool = True) -> QuantumCircuit:
    qc, reg = _new_circuit(inst, "proj_phase")
    append_projector_phase(qc, reg, inst, phi, controlled)
    return qc


# --------------------------------------------------------------------------- #
# QET and the Hadamard test
# --------------------------------------------------------------------------- #

def append_qet(qc: QuantumCircuit, reg: dict, inst: Instance, phases, controlled: bool = True,
               fold_last_phase: bool = False) -> None:
    """The QSVT sequence for even d = len(phases) in the reflection convention (qsp.py):

        U_Phi = prod_{j=1..d/2} [ e^{i phi_{2j-1} (2 Pi - 1)} U^dag e^{i phi_{2j} (2 Pi - 1)} U ]

    so in time order: U, phase phi_d, U^dag, phase phi_{d-1}, ..., U^dag, phase phi_1.
    With controlled=True only the phases are controlled on the Hadamard-test ancilla. With every
    phase off the queries cancel in pairs, so this is the controlled U_Phi.
    fold_last_phase=True replaces the last phase by an rz on the Hadamard-test ancilla. That is
    exact for the Hadamard-test expectation value, because the input state lies in Pi."""
    d = len(phases)
    if d % 2:
        raise ValueError("this builder is for even degree")
    if fold_last_phase and not controlled:
        raise ValueError("fold_last_phase needs the Hadamard-test ancilla")
    for k in range(d, 0, -1):
        append_block_encoding(qc, reg, inst, inverse=(k % 2 == 1))
        if k == 1 and fold_last_phase:
            with _mark(qc, "had"):
                qc.rz(float(phases[0]), reg["had"][0])
        else:
            append_projector_phase(qc, reg, inst, float(phases[k - 1]), controlled)


def qet_circuit(inst: Instance, phases) -> QuantumCircuit:
    """Uncontrolled QET: <0|U_Phi|0> = P(W/s), Re P = f. The Hadamard-test ancilla idles."""
    qc, reg = _new_circuit(inst, "QET")
    append_qet(qc, reg, inst, phases, controlled=False)
    return qc


def hadamard_test_circuit(inst: Instance, phases, site: int = 0, fold_last_phase: bool | None = None,
                          measure: bool = False, final_h: bool = True) -> QuantumCircuit:
    """One shot of the benchmark: <Z> on the Hadamard-test ancilla is Re <site| P(W/s) |site>
    = <site| f(W/s) |site>. `site` is the basis state of the site register, drawn uniformly at
    random per shot to estimate the normalized trace (for the free operator every site gives the
    same value). fold_last_phase defaults to True for the compact variant and False otherwise."""
    if fold_last_phase is None:
        fold_last_phase = inst.variant == "compact"
    qc, reg = _new_circuit(inst, f"logdet_{inst.variant}")
    with _mark(qc, "prep"):
        for j in range(inst.n_site):
            if (site >> j) & 1:
                qc.x(reg["site"][j])
    with _mark(qc, "had"):
        qc.h(reg["had"][0])
    append_qet(qc, reg, inst, phases, controlled=True, fold_last_phase=fold_last_phase)
    if final_h:
        with _mark(qc, "had"):
            qc.h(reg["had"][0])
    if measure:
        creg = ClassicalRegister(1, "m")
        qc.add_register(creg)
        qc.measure(reg["had"][0], creg[0])
    return qc


# --------------------------------------------------------------------------- #
# The benchmark polynomial and its phases
# --------------------------------------------------------------------------- #

def compute_benchmark_polynomial(d: int = BENCHMARK_DEGREE, shift: float = BENCHMARK_SHIFT) -> dict:
    """Run the report's linear program (estimates/apply_log/ch04_roundB_lp_degree.py, read only) for
    the shifted normalized log at lambda_min/s of the benchmark, then solve for the QSP phases."""
    import qsp
    inst = Instance(**BENCHMARK)
    lam = inst.m0 ** 2 / inst.sub
    path = REPO / "estimates" / "apply_log" / "ch04_roundB_lp_degree.py"
    spec = importlib.util.spec_from_file_location("ch04_roundB_lp_degree", path)
    lp = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(lp)
    err, even = lp.min_err(lam, d, shift=shift)
    cheb = np.zeros(d + 1)
    cheb[0::2] = even
    wx, resid = qsp.solve_symmetric_phases(cheb)
    return {
        "what": "LP-minimal even polynomial for log|x|/log(1/lam) + shift on [lam, 1], |p| <= 1 on [-1, 1], "
                "and its QSP phases",
        "lp_source": "estimates/apply_log/ch04_roundB_lp_degree.py min_err",
        "lam": lam, "d": d, "shift": shift, "lp_error": float(err),
        "cheb_coeffs": [float(v) for v in cheb],
        "phases_wx": [float(v) for v in wx],
        "phases_reflection": [float(v) for v in qsp.wx_to_reflection(wx)],
        "qsp_residual": resid,
    }


def benchmark_polynomial(regenerate: bool = False) -> dict:
    """The stored polynomial and phases of the benchmark (polynomial_d84.json), computed if absent."""
    if regenerate or not POLY_FILE.exists():
        data = compute_benchmark_polynomial()
        POLY_FILE.write_text(json.dumps(data, indent=1) + "\n")
    return json.loads(POLY_FILE.read_text())


def benchmark_circuit(variant: str = "literal", site: int = 0, measure: bool = False) -> QuantumCircuit:
    """The full 2028 benchmark shot: V = 4^4, m0 = 0.4, K = 1, d = 84, Hadamard test."""
    inst = Instance(**BENCHMARK, variant=variant)
    return hadamard_test_circuit(inst, benchmark_polynomial()["phases_reflection"], site=site, measure=measure)


if __name__ == "__main__":
    data = benchmark_polynomial(regenerate="--regenerate" in sys.argv)
    print(f"polynomial d = {data['d']}, LP error {data['lp_error']:.4e}, QSP residual {data['qsp_residual']:.1e}")
    for v in VARIANTS:
        c = benchmark_circuit(v)
        print(f"{v:8s} {c.num_qubits} qubits, {len(c.data)} gates, ops {dict(c.count_ops())}")
