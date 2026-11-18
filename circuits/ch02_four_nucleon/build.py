"""Gate-level circuit for the Ch. 2 2028 benchmark: the four-nucleon charge correlator on a 3^3 lattice.

The circuit is one Hadamard test of

    C_ab(tau) = <psi| U_a(q)^dagger  exp(-i (H - E0) tau)  U_b(q) |psi>,   U_f(q) = exp(i q x_f),

for one distinguishable nucleon per species on a periodic L = 3 lattice in `dim` dimensions, with

    H = sum_f [2 dim t - t sum_(hops to the two other sites of each axis)] + C0 sum_(f<g) delta(r_f, r_g)
        + D0 sum_(f<g<h) delta(r_f, r_g) delta(r_g, r_h).

The chapter's correlator is C = (1/4) sum_(a,b in protons) C_ab (app01, "Concrete instantiations").

Everything is emitted in the logical gate set {h, s, sdg, x, cx, t, tdg, ccx, rz, measure}. Arbitrary rotations stay
as rz with their angle. The circuit's global phase is tracked in QuantumCircuit.global_phase so that each block
equals its target operator exactly, not up to a phase.

Register layout (flat qubit index):
    pos[f*dim*2 + axis*2 + bit]   coordinate x = bit0 + 2*bit1 in {0, 1, 2} (the pattern 11 is never populated)
    work[0 .. n_work-1]           clean comparison workspace, returned to |0> by every block
    anc                           Hadamard-test ancilla (last qubit)

Usage:
    from build import benchmark_instance, build_circuit
    built = build_circuit(benchmark_instance(steps=6))      # built.circuit, built.blocks
    python build.py                                         # writes the OpenQASM 3 exports
"""
from __future__ import annotations

import itertools
import json
import math
from contextlib import contextmanager
from dataclasses import dataclass, replace
from pathlib import Path

from qiskit import ClassicalRegister, QuantumCircuit, QuantumRegister

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
STUDY = REPO / "applications" / "ch02_neutrino_nucleus" / "neutrino_feasibility"

THETA3 = 2 * math.asin(1 / math.sqrt(3))      # Ry angle that puts weight 1/3 on the high bit
LOGICAL_GATES = ("h", "s", "sdg", "x", "cx", "t", "tdg", "ccx", "rz", "measure")


# --------------------------------------------------------------------------- #
# Instance
# --------------------------------------------------------------------------- #

@dataclass(frozen=True)
class Instance:
    """One Hadamard-test circuit. Energies in MeV, times in MeV^-1, q in units of 1/a."""
    dim: int = 3                      # spatial dimensions of the periodic lattice
    n_particles: int = 4              # one nucleon per species; species 0 and 1 are the protons
    L: int = 3                        # sites per axis; only L = 3 is implemented (projector form of the kinetic term)
    t: float = 10.5794                # hopping
    C0: float = -98.2265511           # pair contact
    D0: float = 127.839693            # triple contact
    prep_layers: tuple = ()           # ((dt_V, dt_K), ...): exp(-i dt_V V) then exp(-i dt_K K_hop), in order
    tau: float = 0.03                 # evolution time
    steps: int = 6                    # second-order steps K/2 V K ... V K/2
    q: float = 2 * math.pi / 3        # momentum transfer along axis 0
    E0: float = 0.0                   # reference energy of H (including the 2 dim t n_particles constant)
    current_pair: tuple = (0, 0)      # (a, b): U_a^dagger on the way out, U_b on the way in
    component: str = "re"             # "re" or "im" part of C_ab
    merge_ancilla_phases: bool = False  # variant: collect every phase gate on the ancilla into one rz

    def __post_init__(self):
        if self.L != 3:
            raise NotImplementedError("only L = 3 is implemented: the kinetic block uses t I - 3t |s><s|")
        if self.n_particles < 2 or self.dim < 1:
            raise ValueError("need at least two particles and one dimension")
        if self.component not in ("re", "im"):
            raise ValueError("component is 're' or 'im'")

    bits_per_coord = 2

    @property
    def n_coords(self):
        return self.n_particles * self.dim

    @property
    def n_pos(self):
        return self.n_coords * self.bits_per_coord

    @property
    def bits_per_particle(self):
        return self.dim * self.bits_per_coord

    @property
    def n_work(self):
        """Peak clean workspace: one AND chain over the comparison bits of the largest contact."""
        k = 3 if self.n_particles >= 3 else 2
        return (k - 1) * self.bits_per_particle - 1

    @property
    def n_qubits(self):
        return self.n_pos + self.n_work + 1

    @property
    def anc(self):
        return self.n_qubits - 1

    @property
    def work(self):
        return list(range(self.n_pos, self.n_pos + self.n_work))

    def coord_bits(self, particle, axis):
        lo = particle * self.bits_per_particle + axis * self.bits_per_coord
        return lo, lo + 1

    def particle_bits(self, particle):
        lo = particle * self.bits_per_particle
        return list(range(lo, lo + self.bits_per_particle))

    @property
    def hop_constant(self):
        """The diagonal 2 dim t per particle of H. The circuit evolves the hopping part only."""
        return 2 * self.dim * self.t * self.n_particles


def benchmark_instance(steps=6, tau=None, **overrides):
    """The 2028 benchmark: 4 nucleons, 3^3, eight preparation layers, tau = 0.005 MeV^-1 per step.

    The preparation angles are the 8-layer optimum of the study (unitary_prep_results.json, in units of
    1/(100 MeV)) and E0 is the study's exact ground energy (larger_box_results.json). Both are read, not copied.
    """
    prep = json.loads((STUDY / "unitary_prep_results.json").read_text())[-1]
    if prep["layers"] != 8:
        raise RuntimeError("expected the 8-layer preparation as the last entry")
    theta = prep["theta"]
    layers = tuple((theta[2 * j] / 100.0, theta[2 * j + 1] / 100.0) for j in range(prep["layers"]))
    e0 = json.loads((STUDY / "larger_box_results.json").read_text())["ground_energy_MeV"]
    if tau is None:
        tau = 0.005 * steps
    return replace(Instance(prep_layers=layers, tau=tau, steps=steps, E0=e0), **overrides)


# --------------------------------------------------------------------------- #
# Emitter: the only place gates are appended
# --------------------------------------------------------------------------- #

class Emitter:
    """Appends logical gates to a circuit, tracks the global phase and records named blocks."""

    def __init__(self, qc, merge_phases_on=None):
        self.qc = qc
        self.blocks = []                  # (label, start, stop) index ranges into qc.data
        self.merge_qubit = merge_phases_on
        self.deferred = 0.0               # summed phase-gate angle waiting on merge_qubit
        self._open = None

    @contextmanager
    def block(self, label):
        if self._open is not None:
            raise RuntimeError(f"block {label!r} opened inside {self._open!r}: blocks do not nest")
        self._open = label
        start = len(self.qc.data)
        yield
        self.blocks.append((label, start, len(self.qc.data)))
        self._open = None

    # Clifford and Toffoli
    def h(self, q): self.qc.h(q)
    def x(self, q): self.qc.x(q)
    def s(self, q): self.qc.s(q)
    def sdg(self, q): self.qc.sdg(q)
    def cx(self, c, t): self.qc.cx(c, t)
    def ccx(self, a, b, t): self.qc.ccx(a, b, t)

    # Rotations
    def phase(self, theta, q):
        """P(theta) = diag(1, e^{i theta}) = e^{i theta/2} Rz(theta): one rz."""
        if q == self.merge_qubit:
            self.deferred += theta
            return
        self.qc.rz(float(theta), q)
        self.qc.global_phase += theta / 2

    def flush_phases(self):
        """Emit the collected ancilla phase as one rz (merge variant only)."""
        if self.merge_qubit is None:
            return
        theta = math.remainder(self.deferred, 2 * math.pi)
        self.deferred = 0.0
        q, self.merge_qubit = self.merge_qubit, None
        self.phase(theta, q)
        self.merge_qubit = q

    def cphase(self, theta, c, t):
        """Controlled-P(theta): three rz and two cx."""
        self.phase(theta / 2, c)
        self.phase(theta / 2, t)
        self.cx(c, t)
        self.phase(-theta / 2, t)
        self.cx(c, t)

    def ry(self, theta, q):
        """Ry(theta) = S H Rz(theta) H Sdg: one rz."""
        self.sdg(q); self.h(q); self.qc.rz(float(theta), q); self.h(q); self.s(q)

    def ry_pi4(self, q, sign=+1):
        """Ry(+-pi/4) with one exact T gate: Rz(pi/4) = e^{-i pi/8} T."""
        self.sdg(q); self.h(q)
        (self.qc.t if sign > 0 else self.qc.tdg)(q)
        self.h(q); self.s(q)
        self.qc.global_phase -= sign * math.pi / 8


# --------------------------------------------------------------------------- #
# Blocks
# --------------------------------------------------------------------------- #

def prepare3(em, lo, hi, inverse=False):
    """|00> -> (|00> + |01> + |10>)/sqrt(3) on (hi, lo), i.e. the uniform state |s> of one coordinate.

    Ry(THETA3) on hi, then a Hadamard on lo controlled on hi = 0, written as Ry(pi/4) CX Ry(-pi/4).
    One arbitrary rz and two T gates. The controlled-Hadamard part is its own inverse.
    """
    if not inverse:
        em.ry(THETA3, hi)
    em.x(hi)
    em.ry_pi4(lo, +1)
    em.cx(hi, lo)
    em.ry_pi4(lo, -1)
    em.x(hi)
    if inverse:
        em.ry(-THETA3, hi)


def kinetic_coord(em, lo, hi, phi, control=None, work=None):
    """exp(i phi |s><s|) on one coordinate, optionally controlled.

    With phi = 3 t dt this is exp(-i dt K_coord) e^{+i t dt}, because K_coord = t I - 3 t |s><s| on three sites.
    Uncontrolled: basis change, controlled-phase between the two coordinate bits, basis change back (5 rz, 4 T).
    Controlled: the |s> flag is computed into a workspace bit (2 Toffolis) and the phase is a controlled-phase
    from the control onto that bit (5 rz, 4 T, 2 Toffolis).
    """
    prepare3(em, lo, hi, inverse=True)
    em.x(lo); em.x(hi)
    if control is None:
        em.cphase(phi, lo, hi)
    else:
        em.ccx(lo, hi, work)
        em.cphase(phi, control, work)
        em.ccx(lo, hi, work)
    em.x(lo); em.x(hi)
    prepare3(em, lo, hi)


def kinetic_layer(em, inst, dt, control=None):
    """exp(-i dt K_hop) on all coordinates, exactly (including the phase from the t I term)."""
    for f in range(inst.n_particles):
        for axis in range(inst.dim):
            lo, hi = inst.coord_bits(f, axis)
            kinetic_coord(em, lo, hi, 3 * inst.t * dt, control, inst.work[0])
    if control is None:
        em.qc.global_phase += -inst.n_coords * inst.t * dt
    else:
        em.phase(-inst.n_coords * inst.t * dt, control)


def contact(em, bit_lists, angle, work, control=None):
    """Phase e^{i angle} when all listed particles sit on the same site, optionally controlled.

    bit_lists[k] is the list of position bits of particle k. Each later particle is compared bitwise with the
    first (cx, x), the comparison bits are ANDed along a chain of clean workspace bits (n - 1 Toffolis for n
    comparison bits), the phase acts on the last workspace bit, and everything is uncomputed.
    """
    first = bit_lists[0]
    comparisons = [(a, b) for other in bit_lists[1:] for a, b in zip(first, other)]
    checks = [b for _, b in comparisons]
    for a, b in comparisons:
        em.cx(a, b); em.x(b)
    if len(checks) == 1:
        flag = checks[0]
        chain = []
    else:
        w = work[:len(checks) - 1]
        chain = [(checks[0], checks[1], w[0])] + [(w[i - 2], checks[i], w[i - 1]) for i in range(2, len(checks))]
        flag = w[-1]
    for a, b, c in chain:
        em.ccx(a, b, c)
    if control is None:
        em.phase(angle, flag)
    else:
        em.cphase(angle, control, flag)
    for a, b, c in reversed(chain):
        em.ccx(a, b, c)
    for a, b in reversed(comparisons):
        em.x(b); em.cx(a, b)


def potential_layer(em, inst, dt, control=None):
    """exp(-i dt V): every pair with -C0 dt, every triple with -D0 dt."""
    bits = [inst.particle_bits(f) for f in range(inst.n_particles)]
    for pair in itertools.combinations(range(inst.n_particles), 2):
        contact(em, [bits[f] for f in pair], -inst.C0 * dt, inst.work, control)
    for triple in itertools.combinations(range(inst.n_particles), 3):
        contact(em, [bits[f] for f in triple], -inst.D0 * dt, inst.work, control)


def current(em, inst, particle, sign, control):
    """Controlled exp(i sign q x) on axis 0 of one particle, x = bit0 + 2 bit1."""
    lo, hi = inst.coord_bits(particle, 0)
    em.cphase(sign * inst.q, control, lo)
    em.cphase(sign * 2 * inst.q, control, hi)


# --------------------------------------------------------------------------- #
# Whole circuit
# --------------------------------------------------------------------------- #

@dataclass
class Built:
    circuit: QuantumCircuit
    blocks: list            # (label, start, stop) ranges into circuit.data, a partition of the circuit
    instance: Instance


def new_circuit(inst, measure=True):
    regs = [QuantumRegister(inst.n_pos, "pos"), QuantumRegister(inst.n_work, "work"), QuantumRegister(1, "anc")]
    if measure:
        regs.append(ClassicalRegister(1, "m"))
    return QuantumCircuit(*regs, name="ch02_four_nucleon")


def emit_preparation(em, inst):
    with em.block("init/uniform"):
        for f in range(inst.n_particles):
            for axis in range(inst.dim):
                prepare3(em, *inst.coord_bits(f, axis))
    for j, (dt_v, dt_k) in enumerate(inst.prep_layers):
        with em.block(f"prep/{j}/potential"):
            potential_layer(em, inst, dt_v)
        with em.block(f"prep/{j}/kinetic"):
            kinetic_layer(em, inst, dt_k)


def build_preparation(inst):
    """Preparation only (no ancilla operations, no measurement): |0...0> -> |psi>."""
    qc = new_circuit(inst, measure=False)
    em = Emitter(qc)
    emit_preparation(em, inst)
    return Built(qc, em.blocks, inst)


def build_circuit(inst, measure=True):
    """The complete Hadamard-test circuit: preparation, current, controlled evolution, current, readout.

    The ancilla expectation value <Z> is Re C_ab (component 're') or Im C_ab (component 'im').
    """
    qc = new_circuit(inst, measure)
    anc = inst.anc
    em = Emitter(qc, merge_phases_on=anc if inst.merge_ancilla_phases else None)
    a, b = inst.current_pair

    emit_preparation(em, inst)
    with em.block("test/open"):
        em.h(anc)
    with em.block("current/in"):
        current(em, inst, b, +1, anc)
    if inst.steps > 0:
        dt = inst.tau / inst.steps
        with em.block("evolve/first_half/kinetic"):
            kinetic_layer(em, inst, dt / 2, anc)
        for j in range(inst.steps):
            with em.block(f"evolve/{j}/potential"):
                potential_layer(em, inst, dt, anc)
            with em.block(f"evolve/{j}/kinetic"):
                kinetic_layer(em, inst, dt if j < inst.steps - 1 else dt / 2, anc)
    with em.block("current/out"):
        current(em, inst, a, -1, anc)
    with em.block("energy_phase"):
        em.phase((inst.E0 - inst.hop_constant) * inst.tau, anc)
    with em.block("test/close"):
        em.flush_phases()
        if inst.component == "im":
            em.sdg(anc)
        em.h(anc)
        if measure:
            qc.measure(anc, 0)
    return Built(qc, em.blocks, inst)


# Standalone blocks, for tests and for recompiling one block at a time.

def block_circuit(n_qubits, emit):
    qc = QuantumCircuit(n_qubits)
    em = Emitter(qc)
    emit(em)
    return qc


# --------------------------------------------------------------------------- #
# Export
# --------------------------------------------------------------------------- #

def export_qasm3(built, path):
    from qiskit import qasm3
    qc = built.circuit.copy()
    qc.global_phase = 0.0          # a global phase of the whole shot is unobservable; keep the file to plain gates
    text = qasm3.dumps(qc)
    Path(path).write_text(text)
    return len(text)


def main():
    out = {}
    for name, inst in (("ch02_3x3x3_steps6_re_p0p0.qasm3", benchmark_instance(6)),
                       ("ch02_3x3x3_steps6_re_p0p0_merged.qasm3", benchmark_instance(6, merge_ancilla_phases=True))):
        out[name] = export_qasm3(build_circuit(inst), HERE / name)
    for name, size in out.items():
        print(f"{name}: {size / 1e6:.3f} MB")


if __name__ == "__main__":
    main()
