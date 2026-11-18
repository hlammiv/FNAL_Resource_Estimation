"""Explicit logical circuit for the Ch. 7 2028 benchmark step: free real scalar on de Sitter.

Physics instance (paper applications/app04_curved_space.tex, model estimates/ch07_curved_space.py,
functions _ds2028_setup, _ds2028_blocks, ds2028_grid_check with frame='adapted', merged=True):

  H(t) = sum_x [ pi_x^2 / (2 a^3) + a^3 m^2 phi_x^2 / 2 ] + (a/2) phi^T Lap phi,   a = e^t, H_inf = 1,
  Lap  = periodic lattice Laplacian / b^2,  b fixed by m^2 + 12/b^2 = 4 (omega_max = 2 H_inf at a = 1).

One symmetric second-order (Stormer-Verlet) step of length h at the midpoint a_m = e^{t + h/2}:

  U_step = K(h/2) D(h) K(h/2),   K(tau) = exp(-i (tau / a_m^3) sum_x pi_x^2 / 2),
                                 D(h)   = exp(-i (h/2) phi^T (a_m^3 m^2 + a_m Lap) phi).

Digitization: each site holds n_q qubits, K = 2^n_q grid points x_j = dx (j - (K-1)/2), dx = sqrt(2 pi / K)
(arXiv:2108.10793). With qubit q of a site carrying bit q of j (little-endian),
  x = -(dx/2) sum_q 2^q Z_q,
so x^2 on a site is C(n_q,2) ZZ strings and x_x x_y on a link is n_q^2 ZZ strings (arXiv:2210.07985).
The conjugate variable is p = F^dag x F with the centered DFT F_{jj'} = exp(-2 pi i (j-c)(j'-c)/K)/sqrt(K),
c = (K-1)/2.

Grid rescale and shear ("adapted frame"): the register on every site holds z' = S z, z = (phi, pi),
S = [[s, 0], [c, 1/s]] refit to the continuum site covariance at each kinetic-block boundary. A kinetic
block in the frame is  S_out Up(v) S_in^-1 = Lx(g1) Up(u) Lx(g2):  an x^2 phase, a p^2 phase, an x^2 phase.
D in the frame is an x phase with Q / s^2.

Everything here is a plain gate list: H, CX, S, Sdg, Z, T, Tdg and Rz(theta). Rz appears only for angles that
are not multiples of pi/4. No Toffoli and no ancilla is used.

Public entry points
  Instance, BENCHMARK_2028           instance parameters
  step_parameters(inst)              classical coefficients of the step (no circuit code)
  build_step(inst, ...)              -> Built(qc, labels, params)
  block builders                     site_phase, link_phase, centering, qft_swapless, kinetic_block
  merge_rz(built)                    optional peephole pass (commute Rz through CX controls and merge)
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from itertools import combinations, product

import numpy as np
from qiskit import QuantumCircuit

PI = math.pi


# --------------------------------------------------------------------------- #
# Instance
# --------------------------------------------------------------------------- #

@dataclass(frozen=True)
class Instance:
    """Lattice, digitization and step. Units H_inf = 1."""
    dims: tuple = (4, 4, 4)        # periodic lattice, any number of dimensions
    n_q: int = 3                   # qubits per site, K = 2^n_q grid points
    m: float = 1.0                 # field mass (the chapter quotes m ~ 1-1.4 H_inf); needs m <= 1.5
    h: float = 0.5                 # step, Delta t = 0.5 / H_inf
    n_steps: int = 1               # one macro-step at 2028; kinetic half-steps merge between steps
    b: float | None = None         # lattice spacing; None -> the 2028 rule m^2 + 12/b^2 = 4

    @property
    def V(self) -> int:
        return int(np.prod(self.dims))

    @property
    def K(self) -> int:
        return 2 ** self.n_q

    @property
    def spacing(self) -> float:
        return self.b if self.b is not None else math.sqrt(12.0 / (4.0 - self.m ** 2))

    @property
    def n_qubits(self) -> int:
        return self.V * self.n_q


BENCHMARK_2028 = Instance()        # 4^3, K = 8, one step: 192 qubits


# --------------------------------------------------------------------------- #
# Lattice
# --------------------------------------------------------------------------- #

def laplacian(dims, b: float) -> np.ndarray:
    """sum over axes of (2 - T - T^T) / b^2 on the periodic lattice, site index = C-order ravel."""
    V = int(np.prod(dims))
    lap = np.zeros((V, V))
    idx = np.arange(V).reshape(dims)
    for ax in range(len(dims)):
        nb = np.roll(idx, -1, axis=ax).ravel()
        for x, y in zip(idx.ravel(), nb):
            lap[x, x] += 1.0
            lap[y, y] += 1.0
            lap[x, y] -= 1.0
            lap[y, x] -= 1.0
    return lap / b ** 2


def links(dims):
    """Distinct nearest-neighbour pairs (x, y), x < y, ordered by (axis, parity) so that a greedy edge
    colouring gives few layers. Returns [(x, y, colour)]. A side of length 2 gives one pair per two sites
    (its two wraps coincide and the Laplacian weight doubles). A side of length 1 gives none."""
    V = int(np.prod(dims))
    idx = np.arange(V).reshape(dims)
    coords = np.stack(np.unravel_index(np.arange(V), dims), 1)
    seen, out = set(), []
    for ax in range(len(dims)):
        nb = np.roll(idx, -1, axis=ax).ravel()
        for par in (0, 1, 2):
            for x in range(V):
                c = coords[x][ax]
                # parity classes: even, odd, and the wrap link of an odd side on its own
                cls = 2 if (dims[ax] % 2 == 1 and c == dims[ax] - 1) else c % 2
                if cls != par:
                    continue
                y = int(nb[x])
                key = (min(x, y), max(x, y))
                if x == y or key in seen:
                    continue
                seen.add(key)
                out.append(key)
    used = [set() for _ in range(V)]
    coloured = []
    for x, y in out:
        col = 0
        while col in used[x] or col in used[y]:
            col += 1
        used[x].add(col)
        used[y].add(col)
        coloured.append((x, y, col))
    coloured.sort(key=lambda t: t[2])
    return coloured


def lattice_k2(dims, b: float) -> np.ndarray:
    """Laplacian eigenvalues sum_d 4 sin^2(pi n_d / L_d) / b^2 for every lattice momentum."""
    k2 = np.zeros(dims)
    for ax, L in enumerate(dims):
        shape = [1] * len(dims)
        shape[ax] = L
        k2 = k2 + (4.0 * np.sin(PI * np.arange(L) / L) ** 2).reshape(shape)
    return k2.ravel() / b ** 2


def in_state_mode_cov(m: float, k2: np.ndarray):
    """Per-mode (phi, pi) covariance of the in-state: Bunch-Davies mode functions at the lattice momentum
    (eta_in = -1, a = 1) for k != 0 and the a = 1 oscillator ground state for k = 0. Same definition as
    estimates.ch07_curved_space._ds2028_setup."""
    from scipy.special import hankel1, h1vp
    if m > 1.5:
        raise ValueError("m > 1.5 H_inf (imaginary nu) is outside the chapter's benchmark window")
    nu = math.sqrt(9.0 / 4.0 - m * m)
    ph = np.exp(1j * PI * (nu + 0.5) / 2)
    sxx, sxp, spp = np.zeros(len(k2)), np.zeros(len(k2)), np.zeros(len(k2))
    for i, kk in enumerate(k2):
        if kk < 1e-14:
            sxx[i], spp[i] = 1.0 / (2.0 * m), m / 2.0
            continue
        z = math.sqrt(kk)
        q = math.sqrt(PI) / 2 * hankel1(nu, z) * ph
        p = math.sqrt(PI) / 2 * ph * (-1.5 * hankel1(nu, z) - z * h1vp(nu, z))
        sxx[i], spp[i], sxp[i] = abs(q) ** 2, abs(p) ** 2, float(np.real(q * np.conj(p)))
    return sxx, sxp, spp


def balance(sxx: float, sxp: float, spp: float) -> np.ndarray:
    """Lower-triangular site frame S (dilation and shear) that makes the site covariance S C S^T
    proportional to the identity."""
    sc = math.sqrt(math.sqrt(sxx * spp - sxp * sxp) / sxx)
    return np.array([[sc, 0.0], [-sxp / (sc * sxx), 1.0 / sc]])


# --------------------------------------------------------------------------- #
# Classical coefficients of the step
# --------------------------------------------------------------------------- #

@dataclass
class StepParams:
    """Coefficients of the circuit, in the order the blocks act.

    site_gamma[j]   x^2 phase exp(-i gamma x^2 / 2) on every site: j = 0 before the first kinetic block,
                    j = J after the last, the others between kinetic blocks j-1 and j (they absorb the site
                    part of D and the two frame shears).
    kinetic_u[j]    p^2 phase exp(-i u p^2 / 2) of kinetic block j.
    link_kappa[j]   D block j: exp(-i kappa[(x, y)] x_x x_y) on every link.
    frames[j]       site frame at kinetic-block boundary j (J + 1 of them).
    lab_seq         the lab-frame sequence [("K", v) | ("D", mass_coeff, lap_coeff)] the frames were fit to.
    """
    inst: Instance
    dx: float
    b: float
    site_gamma: list
    kinetic_u: list
    link_kappa: list
    frames: list
    lab_seq: list
    links: list = field(default_factory=list)


def lab_sequence(inst: Instance):
    """K(h/2) D(h) K(h/2) per step at a_m = exp(s h + h/2), kinetic halves of neighbouring steps merged."""
    seq = []
    for s in range(inst.n_steps):
        am = math.exp(s * inst.h + inst.h / 2)
        for blk in (("K", inst.h / 2 / am ** 3), ("D", inst.h * am ** 3 * inst.m ** 2, inst.h * am),
                    ("K", inst.h / 2 / am ** 3)):
            if seq and blk[0] == "K" and seq[-1][0] == "K":
                seq[-1] = ("K", seq[-1][1] + blk[1])
            else:
                seq.append(blk)
    return seq


def step_parameters(inst: Instance) -> StepParams:
    b = inst.spacing
    dx = math.sqrt(2 * PI / inst.K)
    k2 = lattice_k2(inst.dims, b)
    sxx, sxp, spp = in_state_mode_cov(inst.m, k2)
    lap = laplacian(inst.dims, b)
    lk = links(inst.dims)
    seq = lab_sequence(inst)

    def site_frame():
        return balance(sxx.mean(), sxp.mean(), spp.mean())

    frames = [site_frame()]
    gammas, us, kappas = [], [], []
    pending = 0.0                       # x^2 coefficient waiting to merge into the next site layer
    for blk in seq:
        if blk[0] == "K":
            v = blk[1]
            sxx, sxp = sxx + 2 * v * sxp + v * v * spp, sxp + v * spp
            s_in, s_out = frames[-1], site_frame()
            A = s_out @ np.array([[1.0, v], [0.0, 1.0]]) @ np.linalg.inv(s_in)
            u = A[0, 1]
            gammas.append(pending + (1 - A[0, 0]) / u)
            us.append(u)
            pending = (1 - A[1, 1]) / u
            frames.append(s_out)
        else:
            _, cm, cl = blk
            w = cm + cl * k2
            spp, sxp = spp - 2 * w * sxp + w * w * sxx, sxp - w * sxx
            s2 = frames[-1][0, 0] ** 2
            pending += (cm + cl * lap[0, 0]) / s2
            kappas.append({(x, y): cl * lap[x, y] / s2 for x, y, _ in lk})
    gammas.append(pending)
    return StepParams(inst=inst, dx=dx, b=b, site_gamma=gammas, kinetic_u=us, link_kappa=kappas,
                      frames=frames, lab_seq=seq, links=lk)


# --------------------------------------------------------------------------- #
# Gate emission
# --------------------------------------------------------------------------- #

class Emitter:
    """Appends gates to a QuantumCircuit and records a block label for every instruction."""

    def __init__(self, n_qubits: int):
        self.qc = QuantumCircuit(n_qubits)
        self.labels: list[str] = []
        self.block = "unlabelled"

    def begin(self, label: str):
        self.block = label

    def _tag(self, k: int = 1):
        self.labels.extend([self.block] * k)

    def h(self, q):
        self.qc.h(q); self._tag()

    def cx(self, a, b):
        self.qc.cx(a, b); self._tag()

    def rz(self, theta, q):
        self.qc.rz(theta, q); self._tag()

    def phase(self, lam: float, q: int):
        """diag(1, e^{i lam}). Multiples of pi/4 become exact Z, S, T gates. Anything else is Rz(lam)
        (equal up to the global phase lam/2, which is recorded on the circuit)."""
        k = lam / (PI / 4)
        if abs(k - round(k)) < 1e-12:
            k = int(round(k)) % 8
            seq = {0: (), 1: ("t",), 2: ("s",), 3: ("s", "t"), 4: ("z",), 5: ("z", "t"), 6: ("sdg",),
                   7: ("tdg",)}[k]
            for g in seq:
                getattr(self.qc, g)(q); self._tag()
        else:
            self.qc.rz(lam, q); self._tag()
            self.qc.global_phase += lam / 2

    def zz(self, theta: float, a: int, b: int):
        """exp(-i theta Z_a Z_b / 2) = CX(a,b) Rz_b(theta) CX(a,b)."""
        self.cx(a, b)
        k = theta / (PI / 4)
        if abs(k - round(k)) < 1e-12:
            self.phase(theta, b)
            self.qc.global_phase -= theta / 2
        else:
            self.rz(theta, b)
        self.cx(a, b)

    def cphase(self, lam: float, a: int, b: int):
        """Controlled phase diag(1,1,1,e^{i lam}) = P_a(lam/2) P_b(lam/2) CX(a,b) P_b(-lam/2) CX(a,b).
        lam = pi/2 is the controlled-S (T, T, Tdg). Smaller angles give three Rz."""
        self.phase(lam / 2, a)
        self.phase(lam / 2, b)
        self.cx(a, b)
        self.phase(-lam / 2, b)
        self.cx(a, b)


def _pair_rounds(n: int):
    """All pairs q < q' of n qubits, ordered by a greedy edge colouring (parallel pairs adjacent)."""
    used = [set() for _ in range(n)]
    out = []
    for a, b in combinations(range(n), 2):
        col = 0
        while col in used[a] or col in used[b]:
            col += 1
        used[a].add(col); used[b].add(col)
        out.append((col, a, b))
    out.sort()
    return [(a, b) for _, a, b in out]


def site_phase(em: Emitter, reg, gamma: float, dx: float, reversed_bits: bool = False):
    """exp(-i gamma x^2 / 2) on one site, up to a global phase. x = -(dx/2) sum_q 2^q Z_q, so
    x^2 = (dx^2/4) [sum_q 4^q + 2 sum_{q<q'} 2^{q+q'} Z_q Z_q'] and the gate is C(n_q,2) ZZ rotations with
    theta = gamma dx^2 2^{q+q'} / 2. reversed_bits = True reads qubit q as bit n_q-1-q (used between the
    swap-free QFTs)."""
    n = len(reg)
    w = [n - 1 - q for q in range(n)] if reversed_bits else list(range(n))
    for a, b in _pair_rounds(n):
        em.zz(gamma * dx * dx * 2.0 ** (w[a] + w[b]) / 2.0, reg[a], reg[b])


def link_phase(em: Emitter, reg_x, reg_y, kappa: float, dx: float):
    """exp(-i kappa x_x x_y) on one link: n_q^2 ZZ rotations, theta = kappa dx^2 2^{q+q'} / 2, emitted in
    n_q rounds of n_q disjoint pairs."""
    n = len(reg_x)
    for r in range(n):
        for q in range(n):
            qp = (q + r) % n
            em.zz(kappa * dx * dx * 2.0 ** (q + qp) / 2.0, reg_x[q], reg_y[qp])


def centering(em: Emitter, reg, inverse: bool = False):
    """R_z layer of the centered DFT: diag(exp(+-2 pi i c j / K)), c = (K-1)/2, i.e. P(+-pi (K-1) 2^q / K)
    on qubit q."""
    n = len(reg)
    K = 2 ** n
    sgn = -1.0 if inverse else 1.0
    for q in range(n):
        em.phase(sgn * PI * (K - 1) * 2 ** q / K, reg[q])


def qft_swapless(em: Emitter, reg, inverse: bool = False):
    """W with QFT = W R, R the bit reversal and QFT_{kj} = exp(2 pi i jk / K) / sqrt(K) (little-endian).
    W: for q = 0..n-1: H(q), then CP(pi / 2^(c-q)) from every c > q. The inverse runs it backwards with
    negated angles. No swap gates: the caller applies the p^2 phase with reversed bit weights."""
    n = len(reg)
    ops = []
    for q in range(n):
        ops.append(("h", q))
        for c in range(q + 1, n):
            ops.append(("cp", PI / 2 ** (c - q), c, q))
    if inverse:
        ops = [(o[0], -o[1], o[2], o[3]) if o[0] == "cp" else o for o in reversed(ops)]
    for o in ops:
        if o[0] == "h":
            em.h(reg[o[1]])
        else:
            em.cphase(o[1], reg[o[2]], reg[o[3]])


def kinetic_block(em: Emitter, reg, u: float, dx: float, with_centering: bool = True, labels=None):
    """exp(-i u p^2 / 2) on one site = F^dag exp(-i u x^2 / 2) F
       = Dc^dag . W . Phi^R(u) . W^dag . Dc   (rightmost first).
    with_centering = False leaves out Dc and Dc^dag (they cancel between kinetic blocks of a shot because
    everything between two kinetic blocks is diagonal)."""
    lab = labels or {}
    if with_centering:
        em.begin(lab.get("c_in", em.block)); centering(em, reg)
    em.begin(lab.get("qft_inv", em.block)); qft_swapless(em, reg, inverse=True)
    em.begin(lab.get("pi2", em.block)); site_phase(em, reg, u, dx, reversed_bits=True)
    em.begin(lab.get("qft", em.block)); qft_swapless(em, reg, inverse=False)
    if with_centering:
        em.begin(lab.get("c_out", em.block)); centering(em, reg, inverse=True)


# --------------------------------------------------------------------------- #
# Full step
# --------------------------------------------------------------------------- #

@dataclass
class Built:
    qc: QuantumCircuit
    labels: list
    params: StepParams

    def block_ranges(self):
        """[(label, first instruction index, one past last)] for the contiguous blocks, in circuit order."""
        out, start = [], 0
        for i in range(1, len(self.labels) + 1):
            if i == len(self.labels) or self.labels[i] != self.labels[start]:
                out.append((self.labels[start], start, i))
                start = i
        return out


def site_register(inst: Instance, x: int):
    return [x * inst.n_q + q for q in range(inst.n_q)]


def build_step(inst: Instance = BENCHMARK_2028, in_shear: bool = True, out_shear: bool = True,
               merge_centering: bool = True, params: StepParams | None = None) -> Built:
    """The circuit of inst.n_steps steps on the whole lattice. Qubit x*n_q + q is bit q of site x.

    in_shear / out_shear   keep the x^2 phase before the first and after the last kinetic block (the frame
                           shears at the two ends of the shot). With both True the circuit is the exact
                           digitized map of ds2028_grid_check. The model's count has neither.
    merge_centering        True: one centering layer at the start and one at the end of the shot (the model's
                           counting). False: every kinetic block carries its own pair.
    """
    p = params or step_parameters(inst)
    em = Emitter(inst.n_qubits)
    regs = [site_register(inst, x) for x in range(inst.V)]
    J = len(p.kinetic_u)

    def site_layer(label, gamma):
        em.begin(label)
        for r in regs:
            site_phase(em, r, gamma, p.dx)

    if in_shear:
        site_layer("shear_in", p.site_gamma[0])
    if merge_centering:
        em.begin("centering_in")
        for r in regs:
            centering(em, r)
    for j in range(J):
        tag = f"k{j + 1}"
        # same gates as kinetic_block on every site, emitted layer by layer so each label is one
        # contiguous range of the instruction list
        if not merge_centering:
            em.begin(f"{tag}_centering_in")
            for r in regs:
                centering(em, r)
        em.begin(f"{tag}_qft_inv")
        for r in regs:
            qft_swapless(em, r, inverse=True)
        em.begin(f"{tag}_pi2")
        for r in regs:
            site_phase(em, r, p.kinetic_u[j], p.dx, reversed_bits=True)
        em.begin(f"{tag}_qft")
        for r in regs:
            qft_swapless(em, r, inverse=False)
        if not merge_centering:
            em.begin(f"{tag}_centering_out")
            for r in regs:
                centering(em, r, inverse=True)
        if j < J - 1:
            site_layer(f"d{j + 1}_site", p.site_gamma[j + 1])
            em.begin(f"d{j + 1}_link")
            for x, y, _ in p.links:
                link_phase(em, regs[x], regs[y], p.link_kappa[j][(x, y)], p.dx)
    if merge_centering:
        em.begin("centering_out")
        for r in regs:
            centering(em, r, inverse=True)
    if out_shear:
        site_layer("shear_out", p.site_gamma[J])
    return Built(qc=em.qc, labels=em.labels, params=p)


# --------------------------------------------------------------------------- #
# Optional recompilation: merge Rz through CX controls
# --------------------------------------------------------------------------- #

_PHASE_OF = {"t": PI / 4, "tdg": -PI / 4, "s": PI / 2, "sdg": -PI / 2, "z": PI}


def merge_rz(built: Built) -> Built:
    """Peephole pass. A diagonal single-qubit gate commutes with a CX on its control, so every run of
    Z-phases on a qubit that is interrupted only by CX gates controlled on that qubit is one phase.
    The merged phase is re-emitted (exact gates if it is a multiple of pi/4, one Rz otherwise) just before
    the gate that ends the run, under the label of the last phase merged. Same unitary, global phase kept."""
    qc = built.qc
    em = Emitter(qc.num_qubits)
    em.qc.global_phase = qc.global_phase
    pending = [0.0] * qc.num_qubits
    plabel = [None] * qc.num_qubits

    def flush(q):
        if plabel[q] is not None:
            em.begin(plabel[q])
            em.phase(pending[q], q)
            pending[q], plabel[q] = 0.0, None

    for inst_, lab in zip(qc.data, built.labels):
        name = inst_.operation.name
        qs = [qc.find_bit(q).index for q in inst_.qubits]
        if name == "rz" or name in _PHASE_OF:
            if name == "rz":                      # Rz(lam) = e^{-i lam/2} P(lam)
                lam = float(inst_.operation.params[0])
                em.qc.global_phase -= lam / 2
            else:
                lam = _PHASE_OF[name]
            pending[qs[0]] += lam
            plabel[qs[0]] = lab
        elif name == "cx":
            flush(qs[1])
            em.begin(lab); em.cx(qs[0], qs[1])
        elif name == "h":
            flush(qs[0])
            em.begin(lab); em.h(qs[0])
        else:
            raise ValueError(f"unexpected gate {name}")
    for q in range(qc.num_qubits):
        flush(q)
    return Built(qc=em.qc, labels=em.labels, params=built.params)
