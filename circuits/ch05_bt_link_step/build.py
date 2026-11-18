"""Gate-level circuits for one Trotter step of 2+1D SU(2) gauge theory digitized on the binary
tetrahedral group 2T (24 elements, 5 qubits per link), Kogut-Susskind Hamiltonian, pure gauge.

Chapter 5 (QGP transport), 2028 benchmark: 4^2 periodic lattice, 32 links, 16 plaquettes.

    H = (g^2 / 2a) sum_links E^2  -  (1 / g^2 a) sum_plaq (Tr U_p + h.c.)

One step is  exp(-i theta_e sum E^2) exp(+i theta_b sum Re Tr U_p)  with
theta_e = g^2 dt / (2 a) and theta_b = 2 dt / (g^2 a)  (function `couplings`).

Register (arXiv:2208.12309, Eq. (pres)): g = (-1)^m i^n j^o l^(p+2q), l = -(1+i+j+k)/2, stored as
|q p o n m> with m the least significant bit. Qubit order inside a link is [m, n, o, p, q]; the state
p = q = 1 is not a group element (8 of the 32 basis states are unused).

Every circuit is built from  x z s sdg h cx cz swap  (Clifford),  t tdg,  ccx,  and rz with an angle
that is not a multiple of pi/4. Nothing else is emitted. Rotations by multiples of pi/4 are written as
Clifford or T gates by `rz_any`. All circuits are equal to their targets up to a global phase.

Primitives (all DERIVED HERE from the group law, not transcribed from the papers' figures):
  u_inv     |g> -> |g^-1>                         4 Toffoli, no ancilla
  u_mul     |g>|h> -> |g>|gh>                     8 Toffoli, no ancilla
  u_tr      |g> -> exp(i theta Re Tr g)|g>        8 rotations, no ancilla
  u_fourier a unitary that maps the isotypic components of the regular representation of 2T onto
            sets of computational basis states    2 Toffoli + 10 T + 2 rotations, 1 ancilla
  u_phase   exp(-i theta f(irrep)) in that basis  4 Toffoli + 5 rotations, 2 ancilla
  electric_link = u_fourier, u_phase, u_fourier^dagger
  plaquette     = 6 u_mul + 6 u_inv + u_tr  (the gate sequence of arXiv:1903.08807 Fig. 1), or
                  6 u_mul + 4 u_inv + u_tr  with uncompute="dagger"

Use `single_plaquette_step` for the 20-qubit instance that can be checked against the exact unitary and
`lattice_step` / `benchmark_circuit` for the L x L periodic lattice.
"""
from __future__ import annotations

import itertools
import math
from contextlib import contextmanager

import numpy as np
from qiskit import QuantumCircuit

CLIFFORD_GATES = {"x", "z", "s", "sdg", "h", "cx", "cz", "swap"}
T_GATES = {"t", "tdg"}
TOFFOLI_GATES = {"ccx"}
ROT_GATES = {"rz"}
LINK_QUBITS = 5


# --------------------------------------------------------------------------------------------------
# The group, from 2x2 matrices (used to build tables for the phase gates and by the tests)
# --------------------------------------------------------------------------------------------------

class BinaryTetrahedral:
    """2T in the faithful two-dimensional irrep of arXiv:2208.12309 Eq. (reps2)."""

    def __init__(self):
        eta = 1 + 1j
        self.qi = np.diag([1j, -1j])
        self.qj = np.array([[0, -1], [1, 0]], dtype=complex)
        self.ql = -0.5 * np.array([[eta, -eta], [np.conj(eta), np.conj(eta)]])
        mp = np.linalg.matrix_power
        self.index = []      # the 24 register values (0..31 with p=q=1 skipped)
        self.mats = {}
        for idx in range(32):
            m, n, o, p, q = [(idx >> k) & 1 for k in range(5)]
            if p and q:
                continue
            self.index.append(idx)
            self.mats[idx] = ((-1) ** m) * mp(self.qi, n) @ mp(self.qj, o) @ mp(self.ql, p + 2 * q)
        self.pos = {idx: k for k, idx in enumerate(self.index)}

    def find(self, mat) -> int:
        for idx, v in self.mats.items():
            if np.allclose(v, mat, atol=1e-9):
                return idx
        raise ValueError("matrix is not a group element")

    def mul(self, a: int, b: int) -> int:
        return self.find(self.mats[a] @ self.mats[b])

    def inv(self, a: int) -> int:
        return self.find(self.mats[a].conj().T)

    def retr(self, a: int) -> float:
        return float(np.trace(self.mats[a]).real)

    def mul_table(self) -> np.ndarray:
        """24 x 24 table in positions 0..23 of self.index: tab[a, b] = position of a*b."""
        t = np.zeros((24, 24), dtype=int)
        for a in self.index:
            for b in self.index:
                t[self.pos[a], self.pos[b]] = self.pos[self.mul(a, b)]
        return t

    # irreps, for the electric eigenvalues only
    def irrep_characters(self) -> dict:
        """Characters of the seven irreps on the 24 elements (register order self.index)."""
        w = np.exp(2j * np.pi / 3)
        chars = {}
        c_of = {idx: ((idx >> 3) & 1) + 2 * ((idx >> 4) & 1) for idx in self.index}
        tr2 = {idx: np.trace(self.mats[idx]) for idx in self.index}
        for k in range(3):
            chars[f"1_{k}"] = np.array([w ** (k * c_of[i]) for i in self.index])
            chars[f"2_{k}"] = np.array([tr2[i] * w ** (k * c_of[i]) for i in self.index])
        # the three-dimensional irrep is the adjoint: chi_3 = chi_2^2 - 1
        chars["3"] = np.array([tr2[i] ** 2 - 1 for i in self.index])
        return chars


_BT = None


def group() -> BinaryTetrahedral:
    global _BT
    if _BT is None:
        _BT = BinaryTetrahedral()
    return _BT


def electric_eigenvalues(gamma: str = "trace1", normalize: bool = True) -> dict:
    """Eigenvalue of the link electric term E^2 on each irrep.

    For a discrete group E^2 is the Cayley-graph Laplacian  sum_{gamma in Gamma} (1 - L_gamma)
    (chapter text after Eq. (Hks); arXiv:2511.17437), so f(rho) = |Gamma| - sum_Gamma chi_rho / dim rho.
    gamma="trace1": Gamma = the 8 elements with Re Tr = 1, the nearest neighbours of the identity
    (a CHOICE made here; the model does not fix Gamma for 2T). normalize=True rescales so that the
    faithful irrep 2_0 carries the SU(2) Casimir 3/4.
    Returns {"1_0", "1_1", "1_2", "2_0", "2_1", "2_2", "3"} -> f.
    """
    bt = group()
    if gamma != "trace1":
        raise ValueError("only gamma='trace1' is implemented")
    sel = np.array([abs(bt.retr(i) - 1.0) < 1e-9 for i in bt.index])
    out = {}
    for name, chi in bt.irrep_characters().items():
        dim = chi[0].real
        val = sel.sum() - (chi[sel].sum() / dim)
        assert abs(val.imag) < 1e-12
        out[name] = float(val.real)
    if normalize:
        scale = 0.75 / out["2_0"]
        out = {k: v * scale for k, v in out.items()}
    return out


def couplings(g2: float, dt_over_a: float) -> tuple[float, float]:
    """(theta_e, theta_b) of one Trotter step: theta_e = g^2 dt / 2a, theta_b = 2 dt / (g^2 a)."""
    return 0.5 * g2 * dt_over_a, 2.0 * dt_over_a / g2


# --------------------------------------------------------------------------------------------------
# Emission helpers
# --------------------------------------------------------------------------------------------------

def is_quarter_multiple(theta: float, tol: float = 1e-11) -> bool:
    k = theta / (math.pi / 4)
    return abs(k - round(k)) < tol


def rz_any(qc: QuantumCircuit, theta: float, q) -> None:
    """R_z(theta) up to a global phase. Multiples of pi/4 become Clifford or T gates."""
    if is_quarter_multiple(theta):
        k = int(round(theta / (math.pi / 4))) % 8
        for name in {0: (), 1: ("t",), 2: ("s",), 3: ("s", "t"), 4: ("z",), 5: ("z", "t"),
                     6: ("sdg",), 7: ("tdg",)}[k]:
            getattr(qc, name)(q)
    else:
        qc.rz(theta, q)


def ry_any(qc: QuantumCircuit, theta: float, q) -> None:
    """R_y(theta) = S H R_z(theta) H S^dagger."""
    qc.sdg(q)
    qc.h(q)
    rz_any(qc, theta, q)
    qc.h(q)
    qc.s(q)


def controlled_h(qc: QuantumCircuit, c, t) -> None:
    """Controlled Hadamard, 2 T:  H = R_y(pi/4) Z R_y(-pi/4)."""
    ry_any(qc, -math.pi / 4, t)
    qc.cz(c, t)
    ry_any(qc, math.pi / 4, t)


@contextmanager
def tag(qc: QuantumCircuit, name: str):
    """Record the instruction range of a primitive in qc.metadata['blocks'] (for count.py)."""
    if qc.metadata is None:
        qc.metadata = {}
    start = len(qc.data)
    yield
    qc.metadata.setdefault("blocks", []).append((name, start, len(qc.data)))


def z_string_rotation(qc: QuantumCircuit, qubits, angle: float) -> None:
    """exp(-i angle/2 Z...Z) on `qubits`: CNOT ladder, R_z(angle), ladder back."""
    qubits = list(qubits)
    for a in qubits[:-1]:
        qc.cx(a, qubits[-1])
    rz_any(qc, angle, qubits[-1])
    for a in reversed(qubits[:-1]):
        qc.cx(a, qubits[-1])


# --------------------------------------------------------------------------------------------------
# Group primitives
# --------------------------------------------------------------------------------------------------
# Write g = z l^c with z = (-1)^m i^n j^o in Q8 and c = p + 2q in Z3. Conjugation by l,
# phi(z) = l z l^-1, sends i -> j -> k -> i and acts on the exponents as (n, o) -> (o, n+o) with no sign.
# phi = CX(n->o) then CX(o->n);  phi^2 = CX(o->n) then CX(n->o). With p and q never both set, the
# sequence [q: A] [p xor q: B] [p: A] applies B,A for p and A,B for q (A = CX(o->n), B = CX(n->o)).

def _phi_power(qc, ctrl_first, ctrl_last, p, q, n, o) -> None:
    """Controlled power of phi on (n, o): the control named `ctrl_last` gets phi (B then A), the control
    named `ctrl_first` gets phi^2 (A then B). ctrl_first, ctrl_last are p and q in some order."""
    qc.ccx(ctrl_first, o, n)       # A under ctrl_first
    qc.cx(q, p)
    qc.ccx(p, n, o)                # B under p xor q
    qc.cx(q, p)
    qc.ccx(ctrl_last, o, n)        # A under ctrl_last


def u_inv(qc: QuantumCircuit, link) -> None:
    """|g> -> |g^-1>.  g^-1 = phi^-c(z^-1) l^-c,  z^-1 = (-1)^(m+n+o+no) i^n j^o.  4 Toffoli."""
    m, n, o, p, q = link
    with tag(qc, "U_inv"):
        qc.ccx(n, o, m)
        qc.cx(n, m)
        qc.cx(o, m)
        _phi_power(qc, p, q, p, q, n, o)     # p (c=1): phi^-1 = phi^2;  q (c=2): phi
        qc.swap(p, q)                        # c -> -c


def u_mul(qc: QuantumCircuit, g, h) -> None:
    """|g>|h> -> |g>|gh>.  gh = z_g phi^{c_g}(z_h) l^{c_g + c_h}.  8 Toffoli, no ancilla."""
    mg, ng, og, pg, qg = g
    mh, nh, oh, ph, qh = h
    with tag(qc, "U_mul"):
        # (a) z_h -> phi^{c_g}(z_h)
        _phi_power(qc, qg, pg, pg, qg, nh, oh)
        # (b) Q8 product: m' = m_g + m_h + (n_g + o_g) n_h + o_g o_h, n' = n_g + n_h, o' = o_g + o_h
        qc.cx(mg, mh)
        qc.cx(og, ng)
        qc.ccx(ng, nh, mh)
        qc.cx(og, ng)
        qc.ccx(og, oh, mh)
        qc.cx(ng, nh)
        qc.cx(og, oh)
        # (c) c_h -> c_h + c_g mod 3: the shift (p, q) -> (p+q+1, p) is affine
        qc.ccx(qg, ph, qh)
        qc.cx(qg, pg)
        qc.ccx(pg, qh, ph)
        qc.cx(qg, pg)
        qc.ccx(pg, ph, qh)
        qc.cx(pg, ph)
        qc.cx(qg, qh)


def u_mul_dagger(qc: QuantumCircuit, g, h) -> None:
    """|g>|gh> -> |g>|h>: the gates of u_mul in reverse order (all are self-inverse)."""
    tmp = QuantumCircuit(*qc.qregs)
    u_mul(tmp, g, h)
    with tag(qc, "U_mul"):
        for inst in reversed(tmp.data):
            qc.append(inst.operation, inst.qubits)


# Walsh (Z-string) expansion of Re Tr g on the register, the eight unused states left free.
# Bit masks over [m, n, o, p, q] (bit 0 = m). Found by exhaustive search over supports (the function is
# odd in m, so only strings containing Z_m can appear): no expansion with fewer than 8 strings exists.
TRACE_SUPPORTS = (1, 11, 13, 15, 17, 27, 29, 31)


def _fit_z_strings(n_bits: int, care: dict, supports) -> np.ndarray:
    """Weights w_S with  const + sum_S w_S prod_{k in S} (-1)^{bit_k} = care[state]  on the care states."""
    states = sorted(care)
    a = np.array([[1.0] + [(-1.0) ** bin(s & S).count("1") for S in supports] for s in states])
    b = np.array([care[s] for s in states])
    sol = np.linalg.lstsq(a, b, rcond=None)[0]
    if np.abs(a @ sol - b).max() > 1e-10:
        raise ValueError("the Z-string support does not reproduce the diagonal")
    return sol[1:]


def trace_weights() -> np.ndarray:
    bt = group()
    return _fit_z_strings(5, {i: bt.retr(i) for i in bt.index}, TRACE_SUPPORTS)


def u_tr(qc: QuantumCircuit, link, theta: float) -> None:
    """|g> -> exp(i theta Re Tr g)|g>, 8 Z-string rotations."""
    with tag(qc, "U_Tr"):
        for S, w in zip(TRACE_SUPPORTS, trace_weights()):
            qs = [link[k] for k in range(5) if (S >> k) & 1]
            # exp(i theta w Z_S) = exp(-i (angle/2) Z_S) with angle = -2 theta w
            z_string_rotation(qc, qs[::-1], -2.0 * theta * w)


# --------------------------------------------------------------------------------------------------
# Electric term
# --------------------------------------------------------------------------------------------------
# The electric term is a function of the irrep only, so the circuit needs a unitary V that sends each
# isotypic component of the regular representation onto a set of basis states, not the full group
# Fourier transform with a fixed basis inside every irrep block.
#   1. H on m. m' = 0: functions even under -1 (on 2T/Z2 = A4); m' = 1: odd functions.
#   2. Odd sector: the three two-dimensional irreps 2_k = 2_0 x omega^{kc}. Its states
#      sum_c omega^{kc} Tr(A sigma(n,o) u^c)|n,o,c>, u the 2x2 matrix of l, become product states after
#      B: |n,o> -> vec(sigma(n,o)^dagger)/sqrt2 (a Clifford) and u^-c on the row qubit controlled on c.
#      u is a single-qubit Clifford of order three, u = R_z(3pi/2) R_y(pi/2), and H u^-1 H = u, so the
#      doubly controlled u^-c is: w = m' AND (p xor q) in an ancilla, controlled-H from q, one
#      controlled-u^-1 from w, controlled-H from q. Exact in Clifford+T.
#   3. B^dagger, then H on n and o. Even sector: (a, b) = (0, 0) are the one-dimensional irreps,
#      everything else is the three-dimensional irrep.
#   4. A Z3 Fourier transform on (p, q) up to output phases: Fourier mode k -> basis state k.
# After V: m'=1 -> irrep 2_k with k = c';  m'=0, a=b=0 -> 1_k;  m'=0 otherwise -> 3.

_ALPHA = math.acos(1.0 / math.sqrt(3.0))


def _bell_basis(qc, n, o) -> None:
    """|n,o> -> vec(sigma(n,o)^dagger)/sqrt2, row index on the n wire. Clifford."""
    qc.sdg(n)
    qc.cz(n, o)
    qc.cx(o, n)
    qc.h(n)
    qc.cx(n, o)


def _bell_basis_dagger(qc, n, o) -> None:
    qc.cx(n, o)
    qc.h(n)
    qc.cx(o, n)
    qc.cz(n, o)
    qc.s(n)


def _controlled_u_inv(qc, c, t, merge_last: bool = False) -> None:
    """Controlled u^-1 on t, u = R_z(3pi/2) R_y(pi/2) (the 2x2 matrix of l), phase exact.
    u^-1 = R_y(-pi/2) R_z(-3pi/2) = A X B X C with A = R_y(-pi/4), B = R_y(pi/4) R_z(3pi/4),
    C = R_z(-3pi/4), A B C = 1."""
    rz_any(qc, -3 * math.pi / 4, t)        # C
    qc.cx(c, t)
    rz_any(qc, 3 * math.pi / 4, t)         # B = R_y(pi/4) R_z(3pi/4): R_z first
    ry_any(qc, math.pi / 4, t)
    qc.cx(c, t)
    # A = R_y(-pi/4); with merge_last the caller's next gate is another R_y(-pi/4), so emit R_y(-pi/2)
    ry_any(qc, -math.pi / 2 if merge_last else -math.pi / 4, t)


def _two_level_h(qc, p, q) -> None:
    """Hadamard on span{|p=1,q=0>, |p=0,q=1>}; |00> and |11> untouched. 2 T."""
    qc.cx(q, p)
    controlled_h(qc, p, q)
    qc.cx(q, p)


def _z3_fourier(qc, p, q) -> None:
    """Sends the Z3 Fourier modes on {|00>, |10>, |01>} (c = p + 2q) to basis states: mode 0 -> c = 0,
    modes 1, 2 -> c = 1, 2 in some order and up to phases. 4 T + 2 rotations. |11> is untouched."""
    _two_level_h(qc, p, q)
    # R_y(-2 alpha) on p when q = 0, alpha = arccos(1/sqrt 3)
    qc.cx(q, p)
    ry_any(qc, -_ALPHA, p)
    qc.cx(q, p)
    ry_any(qc, -_ALPHA, p)
    qc.s(q)
    _two_level_h(qc, p, q)


def _z3_fourier_dagger(qc, p, q) -> None:
    _two_level_h(qc, p, q)
    qc.sdg(q)
    ry_any(qc, _ALPHA, p)
    qc.cx(q, p)
    ry_any(qc, _ALPHA, p)
    qc.cx(q, p)
    _two_level_h(qc, p, q)


def _twiddle(qc, link, anc) -> None:
    """Step 2 above; self-contained (the ancilla is returned to |0>)."""
    m, n, o, p, q = link
    _bell_basis(qc, n, o)
    qc.cx(q, p)
    qc.ccx(m, p, anc)          # w = m' AND (p xor q)
    qc.cx(q, p)
    controlled_h(qc, q, n)
    _controlled_u_inv(qc, anc, n, merge_last=True)
    # second controlled-H; its leading R_y(-pi/4) is merged with the trailing R_y(-pi/4) of the controlled u^-1
    qc.cz(q, n)
    ry_any(qc, math.pi / 4, n)
    qc.cx(q, p)
    qc.ccx(m, p, anc)
    qc.cx(q, p)
    _bell_basis_dagger(qc, n, o)


def _twiddle_dagger(qc, link, anc) -> None:
    tmp = QuantumCircuit(*qc.qregs)
    _twiddle(tmp, link, anc)
    inv = tmp.inverse()
    for inst in inv.data:
        qc.append(inst.operation, inst.qubits)


def u_fourier(qc: QuantumCircuit, link, anc) -> None:
    """V: isotypic components -> basis-state sets (see the comment block above). 1 ancilla."""
    m, n, o, p, q = link
    with tag(qc, "U_F"):
        qc.h(m)
        _twiddle(qc, link, anc)
        qc.h(n)
        qc.h(o)
        _z3_fourier(qc, p, q)


def u_fourier_dagger(qc: QuantumCircuit, link, anc) -> None:
    m, n, o, p, q = link
    with tag(qc, "U_F"):
        _z3_fourier_dagger(qc, p, q)
        qc.h(n)
        qc.h(o)
        _twiddle_dagger(qc, link, anc)
        qc.h(m)


def phase_weights(f: dict) -> tuple[tuple, np.ndarray]:
    """Z-string weights of the electric phase on (m', y, t): y = (m'=0 and a=b=0), t = (c' != 0).
    Strings over [m', t, y] as bit masks (bit 0 = m', bit 1 = t, bit 2 = y)."""
    if abs(f["1_1"] - f["1_2"]) > 1e-12 or abs(f["2_1"] - f["2_2"]) > 1e-12:
        raise ValueError("conjugate irreps must carry the same electric eigenvalue")
    care = {0b000: f["3"], 0b010: f["3"], 0b001: f["2_0"], 0b011: f["2_1"],
            0b100: f["1_0"], 0b110: f["1_1"]}
    supports = (0b001, 0b010, 0b011, 0b100, 0b110)
    return supports, _fit_z_strings(3, care, supports)


def u_phase(qc: QuantumCircuit, link, anc0, anc1, theta_e: float, f: dict | None = None) -> None:
    """exp(-i theta_e f(irrep)) in the basis reached by u_fourier. 4 Toffoli + 5 rotations, 2 ancilla."""
    m, n, o, p, q = link
    f = f or electric_eigenvalues()
    supports, w = phase_weights(f)
    with tag(qc, "U_phase"):
        qc.x(n)
        qc.x(o)
        qc.x(m)
        qc.ccx(n, o, anc0)
        qc.ccx(anc0, m, anc1)                 # y = (not m') (not a) (not b)
        wires = {0: [m], 1: [p, q], 2: [anc1]}  # t = p xor q is the parity of p and q
        # x(m) is still applied: Z_m' = -Z_m here, so undo it before the rotations
        qc.x(m)
        for S, ws in zip(supports, w):
            qs = [x for k in range(3) if (S >> k) & 1 for x in wires[k]]
            z_string_rotation(qc, qs, 2.0 * theta_e * ws)
        qc.x(m)
        qc.ccx(anc0, m, anc1)
        qc.ccx(n, o, anc0)
        qc.x(m)
        qc.x(o)
        qc.x(n)


def electric_link(qc: QuantumCircuit, link, anc0, anc1, theta_e: float, f: dict | None = None) -> None:
    """exp(-i theta_e E^2) on one link: V, phase, V^dagger."""
    u_fourier(qc, link, anc0)
    u_phase(qc, link, anc0, anc1, theta_e, f)
    u_fourier_dagger(qc, link, anc0)


# --------------------------------------------------------------------------------------------------
# Magnetic term
# --------------------------------------------------------------------------------------------------

def plaquette(qc: QuantumCircuit, l1, l2, l3, l4, theta_b: float, uncompute: str = "lly") -> None:
    """exp(+i theta_b Re Tr(g1 g2 g3^-1 g4^-1)). The product g3^-1 g4^-1 g1 g2 (same trace) is built in
    place on link 2, phased by u_tr and uncomputed.

    uncompute="lly":    forward multiplications only, as in arXiv:1903.08807 Fig. 1: 6 u_mul + 6 u_inv.
    uncompute="dagger": u_mul run backwards:                                    6 u_mul + 4 u_inv.
    """
    u_inv(qc, l3)
    u_inv(qc, l4)
    u_mul(qc, l1, l2)
    u_mul(qc, l4, l2)
    u_mul(qc, l3, l2)
    u_tr(qc, l2, theta_b)
    if uncompute == "lly":
        u_inv(qc, l3)
        u_inv(qc, l4)
        u_inv(qc, l1)
        u_mul(qc, l3, l2)
        u_mul(qc, l4, l2)
        u_mul(qc, l1, l2)
        u_inv(qc, l1)
    elif uncompute == "dagger":
        u_mul_dagger(qc, l3, l2)
        u_mul_dagger(qc, l4, l2)
        u_mul_dagger(qc, l1, l2)
        u_inv(qc, l4)
        u_inv(qc, l3)
    else:
        raise ValueError(uncompute)


# --------------------------------------------------------------------------------------------------
# Instances
# --------------------------------------------------------------------------------------------------

def link_qubits(k: int) -> list[int]:
    return list(range(LINK_QUBITS * k, LINK_QUBITS * (k + 1)))


def lattice_links(L: int) -> dict:
    """Periodic L x L lattice: link (x, y, mu) -> index, mu = 0 (x direction), 1 (y direction)."""
    return {(x, y, mu): 2 * (x + L * y) + mu for y in range(L) for x in range(L) for mu in (0, 1)}


def lattice_plaquettes(L: int) -> list[tuple[int, int, int, int]]:
    """Plaquettes as link indices (l1, l2, l3, l4): U_x(s) U_y(s+x) U_x(s+y)^-1 U_y(s)^-1, ordered by
    checkerboard colour so that plaquettes of one colour share no link (even L)."""
    if L < 2:
        raise ValueError("a periodic lattice needs L >= 2; use single_plaquette_step for one plaquette")
    lk = lattice_links(L)
    out = []
    for colour in (0, 1):
        for y in range(L):
            for x in range(L):
                if (x + y) % 2 == colour:
                    out.append((lk[(x, y, 0)], lk[((x + 1) % L, y, 1)], lk[(x, (y + 1) % L, 0)],
                                lk[(x, y, 1)]))
    return out


def _step_layers(qc, links, plaqs, anc, theta_e, theta_b, steps, order, uncompute, f):
    n_pairs = len(anc) // 2

    def electric():
        for k, link in enumerate(links):
            a0, a1 = anc[2 * (k % n_pairs)], anc[2 * (k % n_pairs) + 1]
            electric_link(qc, link, a0, a1, theta_e, f)

    def magnetic(theta):
        for pl in plaqs:
            plaquette(qc, *[links[i] for i in pl], theta, uncompute)

    if order == 1:
        for _ in range(steps):
            electric()
            magnetic(theta_b)
    elif order == 2:
        # exp(i V/2) exp(-i K) exp(i V/2) per step, the magnetic halves of adjacent steps merged
        magnetic(theta_b / 2)
        for s in range(steps):
            electric()
            magnetic(theta_b if s < steps - 1 else theta_b / 2)
    else:
        raise ValueError("order must be 1 or 2")


def trotter_circuit(n_links: int, plaqs, theta_e: float, theta_b: float, steps: int = 1, order: int = 1,
                    n_workspace: int = 2, uncompute: str = "lly", f: dict | None = None) -> QuantumCircuit:
    """`steps` Trotter steps on `n_links` links with the plaquettes `plaqs` (tuples of link indices).
    Qubits: 5 per link, then `n_workspace` ancillas (even, >= 2), used in pairs by the electric term.
    order=1: electric layer then magnetic layer, the count the chapter model uses for a step.
    order=2: symmetric step, steps + 1 magnetic layers in all."""
    if n_workspace < 2 or n_workspace % 2:
        raise ValueError("n_workspace must be even and at least 2")
    qc = QuantumCircuit(LINK_QUBITS * n_links + n_workspace)
    qc.metadata = {"blocks": [], "n_links": n_links, "n_plaquettes": len(plaqs), "steps": steps,
                   "order": order, "n_workspace": n_workspace, "uncompute": uncompute,
                   "theta_e": theta_e, "theta_b": theta_b}
    links = [link_qubits(k) for k in range(n_links)]
    anc = list(range(LINK_QUBITS * n_links, LINK_QUBITS * n_links + n_workspace))
    _step_layers(qc, links, plaqs, anc, theta_e, theta_b, steps, order, uncompute, f or electric_eigenvalues())
    return qc


def single_plaquette_step(theta_e: float, theta_b: float, steps: int = 1, order: int = 1,
                          n_workspace: int = 2, uncompute: str = "lly", f: dict | None = None) -> QuantumCircuit:
    """One open plaquette, four links (20 data qubits + ancillas): Tr(g0 g1 g2^-1 g3^-1)."""
    return trotter_circuit(4, [(0, 1, 2, 3)], theta_e, theta_b, steps, order, n_workspace, uncompute, f)


def lattice_step(L: int, theta_e: float, theta_b: float, steps: int = 1, order: int = 1,
                 n_workspace: int = 16, uncompute: str = "lly", f: dict | None = None) -> QuantumCircuit:
    """Periodic L x L lattice: 2 L^2 links, L^2 plaquettes."""
    return trotter_circuit(2 * L * L, lattice_plaquettes(L), theta_e, theta_b, steps, order, n_workspace,
                           uncompute, f)


def magnetic_layer(L: int, theta_b: float, uncompute: str = "lly") -> QuantumCircuit:
    """Only the plaquette layer of the periodic L x L lattice (no ancilla): a diagonal unitary."""
    qc = QuantumCircuit(LINK_QUBITS * 2 * L * L)
    qc.metadata = {"blocks": []}
    links = [link_qubits(k) for k in range(2 * L * L)]
    for pl in lattice_plaquettes(L):
        plaquette(qc, *[links[i] for i in pl], theta_b, uncompute)
    return qc


# The 2028 benchmark of Chapter 5. Lattice and step sizes are read from the chapter model in count.py;
# the values here are the generator defaults and are checked against the model there.
BENCHMARK = dict(L=4, a_s_fm=0.2, dt_first_fm=0.045, dt_full_fm=0.025, steps_first=1, steps_full=6,
                 n_workspace=16)
G2_PLACEHOLDER = 1.0   # the chapter does not fix g^2 for the 2028 benchmark; counts do not depend on it


def benchmark_circuit(tier: str = "first", order: int = 1, uncompute: str = "lly",
                      n_workspace: int | None = None, g2: float = G2_PLACEHOLDER) -> QuantumCircuit:
    """The 4^2 benchmark: tier="first" is one step of 0.045 fm/c, tier="full" six steps of 0.025 fm/c."""
    b = BENCHMARK
    dt, steps = {"first": (b["dt_first_fm"], b["steps_first"]), "full": (b["dt_full_fm"], b["steps_full"])}[tier]
    theta_e, theta_b = couplings(g2, dt / b["a_s_fm"])
    return lattice_step(b["L"], theta_e, theta_b, steps, order, n_workspace or b["n_workspace"], uncompute)


# --------------------------------------------------------------------------------------------------
# Walking a circuit
# --------------------------------------------------------------------------------------------------

def classify(inst) -> str:
    """'clifford', 't', 'toffoli' or 'rot' for one instruction; anything else is an error."""
    name = inst.operation.name
    if name in CLIFFORD_GATES:
        return "clifford"
    if name in T_GATES:
        return "t"
    if name in TOFFOLI_GATES:
        return "toffoli"
    if name in ROT_GATES:
        if is_quarter_multiple(float(inst.operation.params[0])):
            raise ValueError("an rz by a multiple of pi/4 must be emitted as Clifford/T")
        return "rot"
    raise ValueError(f"gate {name!r} is outside the logical gate set")


def count_gates(qc: QuantumCircuit, start: int = 0, stop: int | None = None) -> dict:
    out = {"toffoli": 0, "rot": 0, "t": 0, "clifford": 0}
    for inst in qc.data[start:stop]:
        out[classify(inst)] += 1
    return out


def count_by_block(qc: QuantumCircuit) -> dict:
    """Counts per tagged primitive: {name: {"calls", "toffoli", "rot", "t", "clifford"}}."""
    out = {}
    covered = 0
    for name, a, b in qc.metadata["blocks"]:
        c = count_gates(qc, a, b)
        d = out.setdefault(name, {"calls": 0, "toffoli": 0, "rot": 0, "t": 0, "clifford": 0})
        d["calls"] += 1
        for k, v in c.items():
            d[k] += v
        covered += b - a
    if covered != len(qc.data):
        raise ValueError("some instructions are outside the tagged primitives")
    return out


def weighted_depth(qc: QuantumCircuit, w_toffoli: float = 1.0, w_t: float = 1.0, w_rot: float = 1.0) -> float:
    """Non-Clifford depth by as-soon-as-possible scheduling on the circuit's own qubits: a gate starts when
    all its qubits are free and lasts its weight (Clifford 0). With w_rot = c_T this is the T-depth
    convention of the chapter model (a Toffoli is one layer, a synthesized rotation c_T layers)."""
    w = {"clifford": 0.0, "t": w_t, "toffoli": w_toffoli, "rot": w_rot}
    free = [0.0] * qc.num_qubits
    idx = {q: i for i, q in enumerate(qc.qubits)}
    for inst in qc.data:
        qs = [idx[q] for q in inst.qubits]
        end = max(free[q] for q in qs) + w[classify(inst)]
        for q in qs:
            free[q] = end
    return max(free)


def used_qubits(qc: QuantumCircuit) -> int:
    seen = set()
    for inst in qc.data:
        seen.update(inst.qubits)
    return len(seen)


def export_qasm(path: str, qc: QuantumCircuit) -> int:
    """Write the circuit as OpenQASM 3; returns the file size in bytes."""
    from qiskit import qasm3
    text = qasm3.dumps(qc)
    with open(path, "w") as fh:
        fh.write(text)
    return len(text)


if __name__ == "__main__":
    import os
    here = os.path.dirname(os.path.abspath(__file__))
    te, tb = couplings(G2_PLACEHOLDER, BENCHMARK["dt_first_fm"] / BENCHMARK["a_s_fm"])
    jobs = {
        "bt_4x4_first_result_1step.qasm": benchmark_circuit("first"),
        "bt_4x4_full_benchmark_6steps.qasm": benchmark_circuit("full"),
        "bt_single_plaquette_1step.qasm": single_plaquette_step(te, tb),
    }
    for name, qc in jobs.items():
        size = export_qasm(os.path.join(here, name), qc)
        print(f"{name}: {qc.num_qubits} qubits, {len(qc.data)} gates, {size / 1e6:.2f} MB")
