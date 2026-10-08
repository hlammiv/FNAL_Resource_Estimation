"""Ch. 6 (QCD for colliders), the 2026-08-28/29 rewrite: every named number the chapter
states, pinned against applications/app10_collider_physics.tex as it reads today.

q1 (author rulings A and B, 2026-10-05): E-rho-OQ dropped from the chapter and the model (G8 penalty exports retired);
the dipole media are quench-prepared (20-step coupling ramp + c/T, c = 2-2 pi, before the sources enter; ETH / typicality
assumed). Static Im V 2.3/6.7e9 -> 5.4e9-1.5e10 T per shot (5.4-15x the limit), 57-98 d; light-like 2.0-5.5e8 (fits),
1.2-2.3 d; shots and LQ unchanged.

STAGE: r26 (referee report 2026-10-04; editorial_review/responses/ch06.md): J1 the species readout map N_h = U Pi_h U^dagger
runs the preparation ramp backwards on every 2033 shot (2033: 4.9e9-2.5e10 T, was 3.9e9-1.5e10; first result 6.9-103 d,
campaign 1.5-22 yr, fragmentation >= 12-1400 yr); J2 gauge-invariant source and the 0.4 fm/c recollision; J4 P_s vs |W|^2
(Markovian color floor); G8 the dipole shot counts binomial at an exact thermal state, the E-rho-OQ factor open (the r25
x40 is a record). Before it:
r25 (apply_log/r25_ch06.md): author rulings R1-R10 (H. Lamm 2026-10-02). R1 two tiers at 2033 (first result
at 30%, campaign at the physics target); R2 conserved-pair clustering in the ratio shots and the 3-momentum trend at
3 sigma; R4 the 30-50 a window (2033: 3.9e9-1.5e10 T, was 1.1-2.0e10); R7 the two dipole first-result rows and the full
q-hat past 2033; R9 t_gate_s and the seven depth exports; R10 the 2028 register 230 LQ (two hop links at once, F* 11-12)
and the 2033 depth-corrected walls. Before it:
r23 (apply_log/r23_ch06.md): author ruling E27 (H. Lamm 2026-10-02, "we don't want to anywhere assume we have
multiple machines"). Every wall time is serial on one machine at 1 us per T-gate plus ~0.1 ms per shot. The machine count
(machines_for_horizon_ff, '5.2--150 machines', 'divides across machines') is retired; where a campaign runs past the
5-year horizon the chapter gives what fits in 5 years on one machine and the cost reduction needed (serial years / 5).
No per-shot T or LQ moves. Before it:
r19 (apply_log/r19_ch06.md): author ruling E21 (1) (H. Lamm 2026-10-01, "make the switch"): the Sigma(72x3)
hop's color squish and parity are computed once per link and held through V(x), V(y), hop, V^dag(x), V^dag(y)
(groups share='link'), the frame undo kept; share='draft' (the r17-r18 headline) is the sensitivity. Rotation counts and
every eps are unchanged; the hop is 3,672 Toffoli + 360 T + 2,933 rotations per link (was 14,314 Toffoli). 2028 unchanged.
2033: 1.0971e10-2.0208e10 T. 3+1D: 3.238e10 T. Before it:
r18 (apply_log/r18_ch06.md): E20 (H. Lamm 2026-10-01, "do 1") prices every rotation at the full fit, the gauge
tables included (n_rot = the paper's log coefficient / 1.15, the papers' constants kept). The 2028 benchmark is 1 ramp + 2
evolution steps (H. Lamm 2026-10-01, "you can do the 1 ramp+2 evolutions if its only 1.27"): 132,174.9 T, 1.32x, accepted
as an overshoot; the E4 step-size cross-check is back. 2033: 1.6216e10-2.9743e10 T. 3+1D: 4.669e10 T. Before it:
r17 (apply_log/r17_ch06.md): ruling (e) of 2026-10-01 adopts the authors' unpublished fermion gate counts.
The hop is groups.hop_link_cost (Sigma(72x3) from the Fermion_Primitives draft with the missing pieces estimated;
Z3 by the R-HOP structure), the mass groups.staggered_mass_site, both at the full fit 1.15 log2(1/eps) + 9.2; the
gauge tables stay in the papers' convention (as Chs. 5 and 10). 2028 is 1 ramp + 1 evolution step (E3 'as many as
fit'): 84,028 T, 0.840x. 2033: 1.547e10-2.838e10 T. 3+1D: 4.464e10 T. Before it:
R-TOL (apply_log/rtol_ch06.md, on top of apply_log/ch06_roundE.md, ch06_twosteps.md and ch06_roundD.md). Report rulings R1-R5
and R8, and the round-D rulings of 2026-09-29: (1) one ramp step at the per-step cost, (2) superseded
by the ruling 'Ch. 6 2028 steps' (2026-09-29): 2 steps in all (1 ramp + 1 evolution), 99,509 T, 0.995x
the 1e5 reference, (3) hopping priced by the report's own
rules (rule R-HOP, derived here), (4) Z3 primitives on the dihedral basis with the Z2 dropped,
(5) 8 body pages; and the round-E rulings (E1 the Z3 U_FFT, 4 T + 2 rotations; E2 eps = 1e-3 per rotation,
randomized synthesis; E3 1 ramp + 2 evolution steps, 100,339 T, 1.003x; E4 combined runs at several step sizes) and
ruling R-TOL (2026-09-29): total synthesis error 1e-2 per shot, each circuit's eps_rot = sqrt(1e-2 / N_rot), which
supersedes E2's 1e-3 (2028: 97,288 T, 0.973x; 2033: 2.768-5.106e9; 3+1D: 9.537e9). PUBLISHED is the boxes as printed after the .tex was edited. Earlier prints are
reproduced by the `legacy_*` intermediates.

Strict xfails are the items still open for the author, filed in apply_log/ch06_roundD.md,
section APPLY, "Needs author", under the item number quoted in the reason. If an edit makes one
pass, the strict xfail fails the run and the item can be closed.

The Z3 circuits that groups.py records as DERIVED HERE, the scope limit of the dihedral paper,
and the Pauli content of the Z3-dressed hop are checked in test_z3_* below, standard library only.
"""

import cmath
import dataclasses
import json
import math
from pathlib import Path

import pytest

from estimates import common as c
from estimates import ch06_collider as m
from estimates import ch09_chiral_gauge as ch9
from estimates.groups import GROUPS, PRIMCOST, magnetic_per_link, electric_per_link

TEX = Path(__file__).resolve().parents[3] / "applications" / "app10_collider_physics.tex"
LOG = "apply_log/ch06_roundD.md (APPLY, needs author)"
# ruling R-TOL: eps_rot = sqrt(eps_syn / N_rot) per circuit, eps_syn = 1e-2 (written out here, not read from common)
EPS_SYN = 1e-2
N28 = 3 * 1264                                   # 2028: 3 steps x 1264 synthesized rotations (r18)
N33 = (500 * 261792, 2500 * 261792)              # 2033 (r26 J1): (300 + 2 x 1e2) and (500 + 2 x 1e3) steps x 261,792
N33_R25 = (400 * 261792, 1500 * 261792)          # the r25 circuit (no readout ramp), record
N33_R23 = (1100 * 261792, 2000 * 261792)         # the r23 circuit, t = 100 a (record)
N3D = 1000 * 785920                              # 3+1D: 1e3 steps x 785,920
N28_PRE = 3 * 1264                               # the pre-r17 2028 circuit (3 steps, slope-only), record
N28_R17 = 2 * 1264                               # the r17 2028 circuit (2 steps), record
EPS = math.sqrt(EPS_SYN / N28)                   # 1.624e-3 (2028)
EPS33 = tuple(math.sqrt(EPS_SYN / n) for n in N33)   # (8.741e-6, 3.909e-6); r25 (9.772e-6, 5.046e-6)
EPS3D = math.sqrt(EPS_SYN / N3D)                 # 3.567e-6
EPS_PRE = math.sqrt(EPS_SYN / N28_PRE)           # 1.624e-3, the R-TOL 2028 print
EPS_ROUND_E = 1e-3                               # round-E ruling E2, superseded; record intermediates only
L2 = math.log2(1 / EPS)
L2_33 = math.log2(1 / EPS33[0])                  # the 2033 per-link intermediates are at the low ramp end's eps
L2_3D = math.log2(1 / EPS3D)


@pytest.fixture(scope="module")
def a():
    return m.Assumptions()


@pytest.fixture(scope="module")
def r28(a):
    return m.model(a, "2028")


@pytest.fixture(scope="module")
def r33(a):
    return m.model(a, "2033")


@pytest.fixture(scope="module")
def r3d(a):
    return m.model(a, "codesign")


def near(x, y, rel):
    assert math.isclose(x, y, rel_tol=rel), (x, y, rel)


# --- contract and conventions -----------------------------------------------------

def test_stage_and_conventions(a, r28, r33, r3d):
    assert m.STAGE == "r26"
    assert all(r.intermediates["stage"] == "r26" for r in (r28, r33, r3d))
    assert a.hop_share.value == "link"                         # E21 (1): squish held once per link
    assert a.link_width_source.value == "compiled"              # R1
    assert a.toffoli_convention.value == "textbook"             # R5: 7 T per Toffoli
    assert c.toffoli_t(1, a.toffoli_convention.value) == 7
    assert a.synthesis_model.value == "rus"                     # E20: gauge tables at the full fit too
    assert a.fermion_synthesis.value == "rus"                   # r17: hop and mass at the full fit
    assert a.eps_syn.value == EPS_SYN == c.EPS_SYN and a.eps_syn.prov is c.Provenance.STATED   # R-TOL
    assert a.eps_rot.value == EPS_ROUND_E and a.eps_rot.prov is c.Provenance.UNCITED and "SUPERSEDED" in a.eps_rot.note
    near(r28.intermediates["eps_rot"], EPS, 1e-12)
    near(r28.intermediates["eps_rot"], c.eps_rot_for(N28), 1e-12)
    near(r33.intermediates["eps_rot"][0], EPS33[0], 1e-12); near(r33.intermediates["eps_rot"][1], EPS33[1], 1e-12)
    near(r3d.intermediates["eps_rot"], EPS3D, 1e-12)
    near(EPS, 1.624e-3, 1e-3); near(EPS33[0], 8.741e-6, 1e-3); near(EPS33[1], 3.909e-6, 1e-3); near(EPS3D, 3.567e-6, 1e-3)
    assert a.t_gate_s.value == 1e-6                             # R8: 1 us per T-gate (main-overview:47); R9 name
    assert not hasattr(a, "gate_time_s")
    assert a.clean_shot_faults.value == 0.1                     # R3
    near(r28.intermediates["t_per_rotation_papers"], 1.15 * L2, 1e-12)
    near(r28.intermediates["t_per_rotation_papers"], c.t_per_rotation(EPS, "rus-slope"), 1e-12)
    near(r28.intermediates["t_per_rotation_report"], 1.15 * L2 + 9.2, 1e-12)   # every rotation (E20)
    near(r28.intermediates["t_per_rotation_papers"], 10.656, 1e-4)             # legacy record (slope-only), not printed
    near(r28.intermediates["t_per_rotation_report"], 19.856, 1e-4)             # '19.9 T per rotation'
    near(r33.intermediates["t_per_rotation_papers"][0], 19.324, 1e-4)         # legacy record
    near(r33.intermediates["t_per_rotation_papers"][1], 20.660, 1e-4)
    near(r33.intermediates["t_per_rotation_report"][0], 28.524, 1e-4)         # '28.5--29.9 T' (r26; r25 28.3--29.4)
    near(r33.intermediates["t_per_rotation_report"][1], 29.860, 1e-4)
    near(r3d.intermediates["t_per_rotation_papers"], 20.811, 1e-4)            # legacy record
    near(r3d.intermediates["t_per_rotation_report"], 30.011, 1e-4)            # '30.0 T'
    assert a.synthesis_errors.value == "incoherent"             # randomized synthesis, N eps^2 (E2, R-TOL)


def test_every_input_tagged_and_stated_inputs_carry_line_refs(a):
    for f in dataclasses.fields(a):
        v = getattr(a, f.name)
        assert isinstance(v, c.Tagged), f.name
        if v.prov is c.Provenance.STATED:
            assert v.src.startswith("app10:"), f"{f.name} is Stated but has no app10:line"
        if v.prov in (c.Provenance.ASSUMED, c.Provenance.UNCITED):
            assert v.note, f"{f.name} is {v.prov.value} with no note"
        if v.prov is c.Provenance.CITED:
            assert v.src, f"{f.name} is Cited with no source"


def test_working_figures_are_retired(a):
    """Round-D ruling 3: the unsourced 5e3 / 3.4e5 / 2.04e6 are gone from every headline; kept Uncited for the record."""
    names = {f.name for f in dataclasses.fields(a)}
    assert not {"t_hop_per_step_2028", "t_hop_insertion_per_step_2033", "t_hop_insertion_per_step_3d"} & names
    for name, val in (("retired_hop_2028", 5e3), ("retired_hop_insertion_2033", 3.4e5), ("retired_hop_insertion_3d", 2.04e6),
                      ("retired_n_trot_2028", 20)):
        v = getattr(a, name)
        assert v.prov is c.Provenance.UNCITED and v.value == val and "RETIRED" in v.note, name
    for era in c.ERAS:
        for p in m.model(a, era).breakdown:
            assert p.status is not c.CircuitStatus.UNSOURCED and "working figure" not in p.note, (era, p.name)


def test_working_assumptions_are_tagged_assumed(a):
    assert a.n_s.prov is c.Provenance.ASSUMED and a.n_s.value == 0.1            # CHANGE_MAP ruling 3
    assert a.hadamard_amplitude.prov is c.Provenance.ASSUMED and a.hadamard_amplitude.value == (0.05, 0.2)
    assert a.utility_share.prov is c.Provenance.ASSUMED
    # what is not the chapter's is tagged as such
    assert a.lambda_qcd_GeV.prov is c.Provenance.ASSUMED and "NOT IN THE CHAPTER" in a.lambda_qcd_GeV.note
    assert a.draft_hop_S72x3_t_const.prov is c.Provenance.UNCITED and "UNPUBLISHED DRAFT" in a.draft_hop_S72x3_t_const.note
    assert not hasattr(a, "z3_fourier_alt_rotations")                   # round E: the transform is the headline
    assert a.eps_rot_before_round_e.prov is c.Provenance.UNCITED and "SUPERSEDED" in a.eps_rot_before_round_e.note
    assert a.twosteps_n_trot_2028.prov is c.Provenance.UNCITED and a.twosteps_n_trot_2028.value == 1


def test_bad_era_raises(a):
    with pytest.raises(ValueError):
        m.model(a, "2040")


# --- the Z3 primitive circuits of groups.py, checked gate by gate -----------------------
# Register: |g> on two qubits (q0, q1), 0 = 00, 1 = 01, 2 = 10 (q0 the low bit), |11> forbidden
# (arxiv_2408_00075 main.tex:237-241).

ENC = {0: (0, 0), 1: (1, 0), 2: (0, 1)}       # g -> (q0, q1)
DEC = {v: k for k, v in ENC.items()}


class Counter:
    def __init__(self):
        self.toffoli = self.cnot = self.x = 0

    def X(self, b, t):
        self.x += 1
        b[t] ^= 1

    def CNOT(self, b, ctl, t):
        self.cnot += 1
        if b[ctl]:
            b[t] ^= 1

    def TOF(self, b, c1, c2, t):
        self.toffoli += 1
        if b[c1] and b[c2]:
            b[t] ^= 1


def _chi(k, b, lo, hi, ctrl=None, inverse=False):
    """The paper's increment chi = X . CNOT . CNOT (Fig. clockmatrixqubit), optionally controlled on one
    qubit: each CNOT becomes a Toffoli and the X a CNOT."""
    def cx(ctl, t):
        k.CNOT(b, ctl, t) if ctrl is None else k.TOF(b, ctrl, ctl, t)

    def x(t):
        k.X(b, t) if ctrl is None else k.CNOT(b, ctrl, t)
    if not inverse:
        cx(hi, lo); cx(lo, hi); x(lo)
    else:
        x(lo); cx(lo, hi); cx(hi, lo)


def test_z3_increment_is_the_papers_chi():
    for g in range(3):
        b = list(ENC[g]); _chi(Counter(), b, 0, 1)
        assert DEC[tuple(b)] == (g + 1) % 3
        b = list(ENC[g]); _chi(Counter(), b, 0, 1, inverse=True)
        assert DEC[tuple(b)] == (g - 1) % 3


def test_z3_multiplication_is_four_toffolis():
    """U_mul |g>|h> = |g>|g+h mod 3>: chi on h controlled on g0 (g = 1), chi^-1 controlled on g1 (g = 2)."""
    k = Counter()
    for g in range(3):
        for h in range(3):
            b = list(ENC[g]) + list(ENC[h])                 # g0 g1 h0 h1
            k.toffoli = 0
            _chi(k, b, 2, 3, ctrl=0)
            _chi(k, b, 2, 3, ctrl=1, inverse=True)
            assert DEC[tuple(b[:2])] == g and DEC[tuple(b[2:])] == (g + h) % 3
            assert k.toffoli == 4
    p = GROUPS["Z3"].primitives["U_mul"]
    assert p.toffoli == 4 and p.t(EPS) == 28 and p.rot == 0 and p.ancilla == 0
    assert "DERIVED HERE" in p.note and "arxiv_2108_13305" in p.src and p.status is c.CircuitStatus.COMPILED


def test_z3_inversion_is_a_swap():
    for g in range(3):
        q0, q1 = ENC[g]
        assert DEC[(q1, q0)] == (-g) % 3
    p = GROUPS["Z3"].primitives["U_inv"]
    assert p.t(EPS) == 0 and "DERIVED HERE" in p.note and "arxiv_2108_13305" in p.src
    assert p.status is c.CircuitStatus.COMPILED


def test_z3_trace_and_electric_phases_are_one_rotation():
    """Re Tr(omega^g 1_3) = 3, -3/2, -3/2 and -2cos(2 pi e/3) = -2, 1, 1: up to a global phase both act on
    g != 0, which on the register is q0 XOR q1. CNOT, one R_Z on q0, CNOT."""
    tr = [3 * math.cos(2 * math.pi * g / 3) for g in range(3)]
    el = [-2 * math.cos(2 * math.pi * e / 3) for e in range(3)]
    near(tr[1], tr[2], 1e-12); near(el[1], el[2], 1e-12)
    assert not math.isclose(tr[0], tr[1]) and not math.isclose(el[0], el[1])
    k = Counter()
    for g in range(3):
        b = list(ENC[g])
        k.CNOT(b, 1, 0)
        flagged = b[0]                                      # the qubit the single R_Z acts on
        k.CNOT(b, 1, 0)
        assert tuple(b) == ENC[g] and flagged == (1 if g else 0)
    for name in ("U_Tr", "U_phi"):
        p = GROUPS["Z3"].primitives[name]
        assert p.rot == 1 and p.toffoli == 0 and "DERIVED HERE" in p.note and "arxiv_2108_13305" in p.src
        assert p.status is c.CircuitStatus.SCALING           # a circuit derived here, drawn in no paper
        near(p.t(EPS), 1.15 * L2 + 9.2, 1e-12)                # E20: one rotation at the full fit
        near(p.t_papers(EPS), 1.15 * L2, 1e-12)               # legacy: the papers' slope-only price
    # the Walsh form of arxiv_2108_13305 Eqs. (21)-(22): with the never-occupied |11> set to the |00> value the
    # trace is a + b Z0 Z1, one rotation; with |11> fixed at 0 it needs Z0, Z1 and Z0 Z1 (three rotations)
    def walsh(vals):                                        # vals on (q0 q1) = 00, 10, 01, 11; returns I, Z0, Z1, Z0Z1
        sgn = lambda q0, q1, s0, s1: (-1) ** (q0 * s0 + q1 * s1)
        st = [(0, 0), (1, 0), (0, 1), (1, 1)]
        return [sum(v * sgn(q0, q1, s0, s1) for v, (q0, q1) in zip(vals, st)) / 4 for s0, s1 in st]
    w_free = walsh([tr[0], tr[1], tr[2], tr[0]])
    w_zero = walsh([tr[0], tr[1], tr[2], 0.0])
    assert sum(abs(x) > 1e-12 for x in w_free[1:]) == 1 and abs(w_free[3]) > 0
    assert sum(abs(x) > 1e-12 for x in w_zero[1:]) == 3


def test_z3_fourier_transform_transpiler_count_kept_for_the_record():
    """arxiv_2408_00075 main.tex:250, :717: H_3 on the two-qubit register is '3 CNOTs and 14 R_z gates' (a transpiler
    count). The same paper's FFT tables are built from it: 2T 42 = 3 x 14 rotations. Since round E the models price
    U_FFT; U_F stays in groups.py for the record."""
    p = GROUPS["Z3"].primitives["U_F"]
    assert p.rot == 14 and p.toffoli == 0 and p.src.startswith("arxiv_2408_00075") and "STATED" in p.note
    assert p.status is c.CircuitStatus.SCALING               # a stated count; the circuit is not printed
    near(p.t(EPS), 14 * (1.15 * L2 + 9.2), 1e-12)             # E20: 14 rotations at the full fit
    near(p.t_papers(1e-4), 213.93, 1e-4)                       # legacy: 16.1 log2(1/eps), the round-D figure
    near(p.t(1e-4), 213.93 + 14 * 9.2, 1e-4)
    near(GROUPS["2T"].primitives["U_FFT"].rot, 3 * 14, 1e-9)
    assert GROUPS["Z3"].has_fft and m.fourier_key("Z3") == "U_FFT" and m.fourier_key("S72x3") == "U_FFT"


# the round-E circuit, index 2 q1 + q0 (g = 0, 1, 2 -> |00>, |01>, |10>), gates in time order
def _mm4(A, B):
    return [[sum(A[i][k] * B[k][j] for k in range(4)) for j in range(4)] for i in range(4)]


def _lift(g, q):
    U = [[0j] * 4 for _ in range(4)]
    for i in range(4):
        for j in range(4):
            bi, bj = (i >> q) & 1, (j >> q) & 1
            if (i ^ j) & ~(1 << q) & 3:
                continue
            U[i][j] = g[bi][bj]
    return U


def _gate(name, *args):
    r2 = 1 / math.sqrt(2)
    one = {"H": [[r2, r2], [r2, -r2]], "S": [[1, 0], [0, 1j]], "Sdg": [[1, 0], [0, -1j]], "Z": [[1, 0], [0, -1]],
           "T": [[1, 0], [0, cmath.exp(1j * math.pi / 4)]], "Tdg": [[1, 0], [0, cmath.exp(-1j * math.pi / 4)]]}
    if name in one:
        return _lift(one[name], args[0])
    if name == "RZ":
        b, q = args
        return _lift([[cmath.exp(-0.5j * b), 0], [0, cmath.exp(0.5j * b)]], q)
    if name == "CZ":
        return [[(1 if i == j else 0) * (-1 if i == 3 else 1) for j in range(4)] for i in range(4)]
    if name == "CX":                                          # control q1, target q0
        perm = {0: 0, 1: 1, 2: 3, 3: 2}
        return [[1 if perm[j] == i else 0 for j in range(4)] for i in range(4)]
    raise KeyError(name)


def _ry(q, rz_gate):                                          # Ry(b) = S H Rz(b) H Sdg: time order Sdg, H, Rz, H, S
    return [("Sdg", q), ("H", q), rz_gate, ("H", q), ("S", q)]


def _z3_fft_circuit():
    th = math.acos(1 / math.sqrt(3))
    ch = _ry(1, ("Tdg", 1)) + [("CZ",)] + _ry(1, ("T", 1))   # CH(q0 -> q1) = Ry(pi/4) CZ Ry(-pi/4)
    return ([("CX",)] + ch + _ry(0, ("RZ", -th, 0)) + [("CZ",), ("Z", 0), ("S", 1)] + _ry(0, ("RZ", th, 0))
            + ch + [("CX",)])


def test_z3_fft_circuit_is_f3_plus_i():
    """Round-E ruling E1: the explicit circuit behind groups.py's Z3 U_FFT, F3 (+) i = B X B^dag, rebuilt gate by
    gate: the 3x3 block is F3 up to a global phase, nothing leaks to or from |11>, 4 T, 2 synthesized rotations."""
    circ = _z3_fft_circuit()
    U = [[1 if i == j else 0 for j in range(4)] for i in range(4)]
    for g in circ:
        U = _mm4(_gate(*g), U)
    w = cmath.exp(2j * math.pi / 3)
    F = [[w ** (j * k) / math.sqrt(3) for k in range(3)] for j in range(3)]
    ph = U[0][0] / F[0][0]
    assert abs(abs(ph) - 1) < 1e-12
    assert max(abs(U[i][j] - ph * F[i][j]) for i in range(3) for j in range(3)) < 1e-12
    assert max(abs(U[3][j]) + abs(U[j][3]) for j in range(3)) < 1e-12
    assert abs(U[3][3] / ph - 1j) < 1e-12
    names = [g[0] for g in circ]
    assert names.count("T") + names.count("Tdg") == 4 and names.count("RZ") == 2
    assert names.count("CX") + names.count("CZ") == 5
    p = GROUPS["Z3"].primitives["U_FFT"]
    assert p.status is c.CircuitStatus.COMPILED and p.t_gates == 4 and p.rot == 2 and p.toffoli == 0
    assert "arxiv_2409_17349" in p.src and "DERIVED HERE" in p.note
    near(p.t_papers(1e-4), 34.5617, 1e-5)                     # legacy: the papers' slope-only price
    near(p.t_papers(1e-3), 26.9213, 1e-5)
    near(GROUPS["Z3"].primitives["U_F"].t_papers(1e-4) / p.t_papers(1e-4), 6.19, 0.001)
    near(p.t(1e-4), 34.5617 + 2 * 9.2, 1e-5)                  # E20: 4 T + 2 rotations at the full fit
    near(p.t(1e-3), 26.9213 + 2 * 9.2, 1e-5)
    near(p.t_report(1e-3), 4 + 2 * c.t_per_rotation(1e-3), 1e-12)   # the 4 T are not re-costed as Toffolis


def test_z3_fft_in_the_2028_step(a, r28):
    i = r28.intermediates
    assert i["z3_fourier_key"] == "U_FFT" and i["z3_fourier_rotations"] == 2 and i["z3_fft_exact_t"] == 4
    near(i["z3_t_U_FFT"], 4 + 2 * (1.15 * L2 + 9.2), 1e-12)               # 43.7 -> '44 T' (E20)
    assert round(i["z3_t_U_FFT"]) == 44
    near(i["z3_t_U_FFT_report"], 4 + 2 * (1.15 * L2 + 9.2), 1e-12)       # full fit, record
    near(i["z3_t_U_FFT_round_e"], 26.9213, 1e-5)                          # the round-E '27 T'
    near(i["z3_t_U_F_transpiler"], 14 * (1.15 * L2 + 9.2), 1e-12)        # 278.0 -> '278 T'
    assert round(i["z3_t_U_F_transpiler"]) == 278
    near(i["z3_fft_reduction"], 6.3595, 0.001)
    near(i["transpiler_electric_t_per_link_2028"], 29 * (1.15 * L2 + 9.2), 1e-12)
    near(i["t_step_transpiler_fourier_2028"], 59051.89, 1e-6)             # what the transpiler count would give
    near(i["t_step_fft_at_1e-4_2028"], 49903.82, 1e-6)                    # the step at 1e-4, every rotation at the full fit
    near(i["twosteps_t_total_2028"], 99509.45, 1e-6)                      # the superseded two-step print


def test_dihedral_basis_is_z_2n_only():
    """Round-D ruling 4. arxiv_2108_13305 builds D_N only for N = 2^n (Sec. III). Its circuits with the Z2 dropped are
    Z_{2^n} circuits: on the Z3 states the 2^n inversion (1's complement, then +1 mod 4) sends g = 1 to the forbidden
    |11>, and its mod-4 adder is wrong on (1,2), (2,1), (2,2). So the Z3 circuits come from arxiv_2408_00075, and the
    dihedral paper supplies the algebra (k -> N-k, k1+k2 mod N) that they implement."""
    inv_2n = {g: ((3 - g) + 1) % 4 for g in range(3)}       # 1's complement on 2 bits, then increment mod 4
    assert inv_2n[1] == 3                                    # |11>, forbidden in the Z3 encoding
    bad = [(g, h) for g in range(3) for h in range(3) if (g + h) % 4 != (g + h) % 3]
    assert bad == [(1, 2), (2, 1), (2, 2)]
    note = GROUPS["Z3"].primitives["U_mul"].note
    assert "Dropped with the Z2" in note and "two's complement" in note
    assert GROUPS["Z3"].link_qubits_chapter.src == "arxiv_2108_13305" and GROUPS["Z3"].link_qubits == 2 == math.ceil(math.log2(3))


# the Z3-dressed hop: qubits (g1, g0, a, b), link register of arxiv_2408_00075, one color-flavor copy
_PAULI = {"I": ((1, 0), (0, 1)), "X": ((0, 1), (1, 0)), "Y": ((0, -1j), (1j, 0)), "Z": ((1, 0), (0, -1))}


def _kron(A, B):
    return [[A[i][j] * B[k][l] for j in range(len(A)) for l in range(len(B))] for i in range(len(A)) for k in range(len(B))]


def _pauli_terms(H):
    import itertools
    out = {}
    for lab in itertools.product("IXYZ", repeat=4):
        M = _PAULI[lab[0]]
        for ch in lab[1:]:
            M = _kron(M, _PAULI[ch])
        cf = sum(M[j][i] * H[i][j] for i in range(16) for j in range(16)) / 16
        if abs(cf) > 1e-12:
            out["".join(lab)] = cf
    return out


def _hop(link_phase):
    """sum_g |g><g| (x) (u_g a^dag b + h.c.) on (g1, g0, a, b); link_phase maps (g1, g0) -> u_g."""
    H = [[0j] * 16 for _ in range(16)]
    for (g1, g0), u in link_phase.items():
        i = int(f"{g1}{g0}10", 2); j = int(f"{g1}{g0}01", 2)
        H[i][j] += u; H[j][i] += u.conjugate()
    return H


def test_z3_hop_is_eight_pauli_strings_in_two_commuting_blocks(a, r28):
    """R-HOP for Z3 (derived here). D(g) = omega^g on every color, so the link is diagonal and there is no move. After
    one Clifford CNOT(g0 -> g1) the register reads g = 0, 1, 2 as 00, 11, 10 (01 never occupied); choosing the free
    value there leaves 8 strings, which split into the Re block {XX, YY} and the Im block {XY, YX}. On the raw register
    the fewest found is 10."""
    w = cmath.exp(2j * math.pi / 3)
    target = {(0, 0): 1 + 0j, (1, 1): w, (1, 0): w ** 2}
    Ht = _hop(target)
    best = None
    for re_x in [x / 4 for x in range(-12, 13)]:
        for im_x in [x / 4 for x in range(-8, 9)]:
            d = _pauli_terms(_hop({**target, (0, 1): complex(re_x, im_x)}))
            if best is None or len(d) < len(best[1]):
                best = ((re_x, im_x), d)
    (fv, terms) = best
    assert len(terms) == 8 == a.n_pauli_hop_diagonal.value == r28.intermediates["hop_pauli_strings_2028"]
    Hf = _hop({**target, (0, 1): complex(*fv)})
    allowed = [s for s in range(16) if (s >> 2) & 3 != 0b01]
    assert all(abs(Hf[i][j] - Ht[i][j]) < 1e-12 for i in allowed for j in allowed)   # same operator where the link lives
    re_block = {k for k in terms if k[2:] in ("XX", "YY")}
    im_block = {k for k in terms if k[2:] in ("XY", "YX")}
    assert len(re_block) == len(im_block) == 4 and re_block | im_block == set(terms)
    raw = {(0, 0): 1 + 0j, (0, 1): w, (1, 0): w ** 2, (1, 1): complex(-2, 0)}
    assert len(_pauli_terms(_hop(raw))) == 10


def test_2028_hop_rule_numbers(a, r28):
    """app10:93 'Hamming-weight phasing ..., 4 rotations and 7 Toffolis (128 T, 7 ancilla) per group' (E26, r22: a group
    of k holds k - w(k) ancilla, the same number as its Toffolis; the chapter printed 11); '8 x 128 ~ 1.0e3 T
    per link, 3.3e4 T per step; the staggered mass is one group per site, 2.1e3 T' (r17: groups.hop_link_cost, Z3 by the
    R-HOP structure, rotations at the full fit; r18 R-TOL eps = 1.624e-3)."""
    i = r28.intermediates
    assert (i["hop_k_2028"], i["hwp_rotations_k9"], i["hwp_toffolis_k9"], i["hwp_ancilla_k9"]) == (9, 4, 7, 7)
    assert i["hwp_ancilla_k9"] == c.hwp_ancilla(9) == i["hwp_toffolis_k9"]
    near(i["hwp_t_k9"], 4 * (1.15 * L2 + 9.2) + 7 * 7, 1e-12)
    assert round(i["hwp_t_k9"]) == 128
    assert i["hop_moves_per_link_2028"] == 0                 # Z3: diagonal link, no color move
    assert (i["hop_toffoli_per_link_2028"], i["hop_rotations_per_link_2028"]) == (8 * 7, 8 * 4)
    near(i["t_hop_per_link_2028"], 8 * i["hwp_t_k9"], 1e-12)
    near(i["t_hop_per_link_2028"], 1027.40, 1e-5)            # '8 x 128 ~ 1.0e3'
    assert round(i["t_hop_per_link_2028"] / 1e3, 1) == 1.0
    assert round(i["t_hop_per_link_2028_at_1e-4"]) == 881   # the round-D print (slope-only)
    assert round(i["t_hop_per_link_2028_round_e"]) == 759   # the round-E print
    near(i["t_hop_per_link_2028_rhop_print"], 8 * (4 * 1.15 * L2 + 49), 1e-12)   # the same structure, slope-only
    near(i["t_hop_per_step_2028"], 32876.79, 1e-6)          # '3.3e4'
    near(i["t_mass_per_step_2028"], 16 * i["hwp_t_k9"], 1e-12)
    near(i["t_mass_per_step_2028"], 2054.80, 1e-5)          # '2.1e3'
    near(i["hopping_share_of_step_2028"], 0.7462, 0.001)
    # the retired U_mul-style price on Z3, record
    near(i["umul_style_hop_per_link_z3"], 2 * 3 * 28 + 2 * i["hwp_t_k9"], 1e-12)
    near(i["explicit_over_umul_style_z3"], 2.418, 0.002)
    near(i["hop_over_retired_2028"], 6.575, 0.002)          # 6.6x the retired 5e3


# --- 2028 benchmark: plug-in app10:88, box app10:114-130 --------------------------------

def test_2028_register_from_eq_Nq(r28):
    """R10 (r25): '2.16.2 + 3.3.16 + 22 = 230 LQ. The 22 ancilla run the hop on two links at once, each with one phasing
    group's 7 (below) and 4 for synthesis; with 8 (216 LQ) every rotation runs in series, the T-count is 1.7x the
    T-depth ..., while with 22 it is 11--12x'; box '230'. The r22 register (8 = 7 + 1, 216) is the record."""
    i = r28.intermediates
    assert i["q_G_Z3"] == GROUPS["Z3"].link_qubits == 2
    assert GROUPS["Z3"].link_qubits_compiled.src == "arxiv_2408_00075"
    assert (i["V_2028"], i["n_links_2028"], i["n_plaq_2028"]) == (16, 32, 16)
    assert (i["lq_gauge_2028"], i["lq_fermion_2028"], i["lq_anc_2028_r22"]) == (64, 144, 8)
    assert (i["lq_anc_hwp_2028"], i["lq_anc_rus_2028"]) == (ch9.hwp_ancilla(9), c.RUS_ANCILLA) == (7, 1)
    assert i["lq_anc_hwp_2028"] == c.hwp_ancilla(9) == i["hwp_ancilla_k9"]
    assert i["lq_anc_2028_r22"] == i["lq_anc_hwp_2028"] + i["lq_anc_rus_2028"]   # the Stated 8 is the derived 7 + 1
    assert i["lq_2028_r22"] == 216
    assert i["workspace_links_2028"] == 2
    assert i["lq_anc_workspace_2028"] == 22 == 2 * (7 + i["hwp_rotations_k9"])
    assert i["lq_2028"] == 230 == 64 + 144 + 22 and r28.lq == (230, 230)      # printed exact: '= 230', box '230'
    assert 150 <= i["lq_2028"] <= 250                       # '150--250 envelope'
    assert i["lq_headroom_2028"] == 20


def test_2028_depth_exports(a, r28):
    """R9/R10 (factory.json Ch. 6 2028): at the r22 register (8 ancilla) D_T = 8.0e4 and F* = 1.65 ('1.7x'); with two hop
    links at once D_T = 1.06--1.20e4 and F* = 11.0--12.4 ('11--12x'). The seven exports are common.depth_exports'."""
    i = r28.intermediates
    near(i["t_depth_per_shot_r22"], 80070.88, 1e-5)
    near(i["f_star_r22"], 1.6507, 1e-4)                     # '1.7x the T-depth'
    near(i["toffolis_per_shot_2028"], 8016, 1e-9); assert i["exact_t_per_shot_2028"] == 768
    near(i["t_depth_per_shot"][0], 10640.17, 1e-6); near(i["t_depth_per_shot"][1], 12008.17, 1e-6)
    near(i["f_star"][0], 11.007, 1e-4); near(i["f_star"][1], 12.422, 1e-4)   # '11--12x'
    assert i["baseline_ok"] == (True, True)
    dx = c.depth_exports(r28.hard_ops, i["t_depth_per_shot"], i["shots_2028"], a.t_gate_s, a.shot_overhead_s)
    for k in ("t_per_shot", "t_depth_per_shot", "f_star", "floor_wall_s", "factories_for_1yr"):
        assert i[k] == dx[k], k
    near(i["floor_wall_s"][1], 120.08, 1e-4)
    near(i["factories_for_1yr"][0], 4.188e-5, 1e-3)
    assert i["wall_first_result_s"] is None                 # one tier
    near(i["wall_campaign_s"][0], 132.275, 1e-5)            # = the serial wall: F* >= 10
    # one hop link at a time would not restore F* >= 10
    r1 = m.model(dataclasses.replace(a, workspace_links_2028=c.Stated(1, "app10:88", "test")), "2028")
    assert r1.intermediates["f_star"][0] < 10 and r1.lq == (219, 219)


def test_2028_gauge_terms_come_from_groups(a, r28):
    """R2: '178 T magnetic and 107 T electric per link, 9.1e3 T per step on 32 links' (U_FFT, every rotation at the full
    fit (E20), R-TOL eps = 1.624e-3). The multiplicities are groups.PRIMCOST (tab:primcost): per link 3 U_inv, 6 U_mul, 1/2 U_Tr, 2 Fourier
    transforms, one U_phi."""
    i = r28.intermediates
    assert i["magnetic_t_per_link_2028"] == magnetic_per_link("Z3", "KS", 2, EPS)
    assert i["electric_t_per_link_2028"] == electric_per_link("Z3", "KS", 2, EPS)
    mult = PRIMCOST["KS"]
    assert (mult["U_inv"](2), mult["U_mul"](2), mult["U_Tr"](2), mult["U_F"](2)) == (3, 6, 0.5, 2)
    R = 1.15 * L2 + 9.2
    near(i["magnetic_t_per_link_2028"], 3 * 0 + 6 * 28 + 0.5 * R, 1e-12)
    near(i["electric_t_per_link_2028"], 2 * (4 + 2 * R) + R, 1e-12)
    assert round(i["magnetic_t_per_link_2028"]) == 178
    assert round(i["electric_t_per_link_2028"]) == 107
    assert round(i["t_per_link_2028"]) == 285
    assert (i["toffoli_per_link_step_2028"], i["rotations_per_link_step_2028"]) == (24, 5.5)
    near(i["t_per_link_2028"], 24 * 7 + 2 * 4 + 5.5 * R, 1e-12)
    near(i["t_gauge_per_step_2028"], 9126.70, 1e-6)         # '9.1e3'
    near(i["report_convention_t_gauge_per_step_2028"], 9126.70, 1e-6)   # E20: the same
    near(i["t_gauge_per_step_2028_round_d"], 19801.14, 1e-6)   # round-D record, papers' slope-only price
    near(i["electric_share_of_gauge_2028"], 0.3761, 0.001)
    near(i["tab_primcost_plaquette_t"], 355.86, 0.001)
    near(i["alt_abelian_plaquette_t"], 187.86, 0.001)       # an abelian-specific plaquette (not used)


def test_2028_t_count_chain(r28):
    """'The step is 4.4e4 T. One ramp step and two evolution steps give 1.3e5 T, 1.3x the 1e5 first-generation limit, an
    overshoot we accept; one evolution step would fit, at 8.7e4 T'; box '1.3e5 (1 ramp + 2 evolution steps; derived
    here), 1.3x the 1e5 limit' (r18: E20 full fit; ruling of 2026-10-01 on the 2028 steps)."""
    i = r28.intermediates
    near(i["t_per_step_2028"], i["t_gauge_per_step_2028"] + i["t_hop_per_step_2028"] + i["t_mass_per_step_2028"], 1e-12)
    near(i["t_per_step_2028"], 44058.29, 1e-6)              # '4.4e4'
    assert round(i["t_per_step_2028"] / 1e4, 1) == 4.4
    assert (i["ramp_steps_2028"], i["n_trot_2028"], i["n_steps_2028"]) == (1, 2, 3)
    near(i["t_ramp_2028"], i["t_per_step_2028"], 1e-12)     # the ramp step at the per-step cost (ruling 1)
    near(i["t_evolution_2028"], 2 * i["t_per_step_2028"], 1e-12)
    near(i["t_total_2028"], 3 * i["t_per_step_2028"], 1e-12)
    near(i["t_total_2028"], 132174.88, 1e-6)                # '1.3e5'
    assert round(i["t_total_2028"] / 1e5, 1) == 1.3
    near(i["t_total_2028_over_cap"], 1.3217, 1e-4)          # '1.3x', an accepted overshoot (<= ~1.35x by the ruling)
    assert round(i["t_total_2028_over_cap"], 1) == 1.3 and i["t_total_2028_over_cap"] <= 1.35
    assert i["overshoot_accepted_2028"] is True
    assert i["steps_in_cap_2028"] == 2 < i["n_steps_2028"]  # the third step overshoots
    near(i["steps_in_cap_exact_2028"], 2.2697, 0.001)
    near(i["headroom_2028"], -32174.88, 1e-6)
    near(i["t_total_2028_three_steps"], i["t_total_2028"], 1e-12)
    near(i["t_per_step_2028_three_steps"], i["t_per_step_2028"], 1e-12)
    # the same circuit with 2 steps (1 ramp + 1 evolution) at its own R-TOL eps: '8.7e4', fits
    near(i["t_total_2028_two_steps"], 87266.29, 1e-6)
    assert round(i["t_total_2028_two_steps"] / 1e4, 1) == 8.7 and i["t_total_2028_two_steps"] < 1e5
    near(i["t_total_2028_four_steps"], 177439.77, 1e-6)    # record
    # the r17 print (1 + 1 steps, gauge tables slope-only), the pre-r17 print (R-TOL, slope-only, 3 steps) and round E
    near(i["t_total_2028_r17"], 84027.89, 1e-6)
    near(i["t_per_step_2028_pre_r17"], 32429.49, 1e-6)
    near(i["t_total_2028_pre_r17"], 97288.48, 1e-6)
    near(i["t_per_step_2028_round_e"], 33446.26, 1e-6)
    near(i["t_total_2028_round_e"], 100338.79, 1e-6)
    assert (i["round_d_n_steps_2028"], round(i["round_d_t_total_2028"])) == (4, 133785)   # round D's count, round-E step
    assert r28.hard_ops == (i["t_total_2028"], i["t_total_2028"])
    near(r28.breakdown_total(), i["t_total_2028"], 1e-12)
    assert c.tex_hard_ops(*r28.hard_ops) == r"$1.3{\times}10^{5}$"       # the box's spelling, and Table 1.1's


def test_2028_one_ramp_step_and_two_evolution_steps(a, r28):
    """Round-D ruling 1 and the ruling of 2026-10-01: one ramp step and two evolution steps, 1.32x accepted."""
    assert a.ramp_steps_2028.value == 1 and a.n_trot_2028.value == 2
    assert "if its only 1.27" in a.n_trot_2028.note
    assert a.round_d_n_trot_2028.value == 3 and a.round_d_n_trot_2028.prov is c.Provenance.UNCITED
    rows = {p.name: p for p in r28.breakdown}
    ramp = sum(p.t_total for n, p in rows.items() if n.endswith("_ramp") or "_ramp_" in n)
    evo = sum(p.t_total for n, p in rows.items() if n.endswith("_evolution") or "_evolution_" in n)
    near(ramp, r28.intermediates["t_per_step_2028"], 1e-12)
    near(evo, 2 * r28.intermediates["t_per_step_2028"], 1e-12)
    # an N-step ramp adds N x the per-step cost; under R-TOL a longer shot is a different circuit with its own eps
    for n in (0, 2, 5):
        rb = m.model(dataclasses.replace(a, ramp_steps_2028=c.Stated(n, "app10:124", "test")), "2028")
        ib = rb.intermediates
        near(ib["eps_rot"], math.sqrt(EPS_SYN / ((2 + n) * 1264)), 1e-12)
        near(rb.hard_ops[0], (2 + n) * ib["t_per_step_2028"], 1e-12)
        assert (ib["t_per_step_2028"] < r28.intermediates["t_per_step_2028"]) == (n < 1)
        near(rb.breakdown_total(), rb.hard_ops[0], 1e-12)
    # one evolution step is the 2-step alternative the chapter quotes ('8.7e4')
    r1 = m.model(dataclasses.replace(a, n_trot_2028=c.Stated(1, "app10:124", "test")), "2028")
    near(r1.hard_ops[0], r28.intermediates["t_total_2028_two_steps"], 1e-12)


def test_2028_breakdown_statuses(r28):
    rows = {p.name: p for p in r28.breakdown}
    gauge = {f"{p}_Z3_{t}" for p in ("U_inv", "U_mul", "U_Tr", "U_FFT", "U_phi") for t in ("evolution", "ramp")}
    ferm = ({f"hop_z3_hop_hwp_Z3_{t}_{k}" for t in ("evolution", "ramp") for k in ("toffoli", "rotations")}
            | {f"{p}_{t}" for p in ("mass_hwp_rotations", "mass_hwp_toffolis") for t in ("evolution", "ramp")})
    assert set(rows) == gauge | ferm                          # no color move for Z3, no unsourced line
    assert rows["U_mul_Z3_evolution"].count == 6 * 32 * 2 and rows["U_mul_Z3_evolution"].t_each == 28
    assert rows["U_FFT_Z3_ramp"].count == 2 * 32 * 1
    assert rows["U_FFT_Z3_evolution"].status is c.CircuitStatus.COMPILED
    hr, ht = rows["hop_z3_hop_hwp_Z3_evolution_rotations"], rows["hop_z3_hop_hwp_Z3_evolution_toffoli"]
    assert hr.count == 32 * 8 * 4 * 2 and ht.count == 32 * 8 * 7 * 2 and ht.t_each == 7
    near(hr.t_each, 1.15 * L2 + 9.2, 1e-12)                  # r17: full fit
    assert rows["mass_hwp_rotations_ramp"].count == 16 * 4
    near(rows["mass_hwp_rotations_ramp"].t_each, 1.15 * L2 + 9.2, 1e-12)
    for name in ("U_inv_Z3_evolution", "U_mul_Z3_evolution", "U_Tr_Z3_evolution", "U_phi_Z3_evolution"):
        assert "DERIVED HERE" in rows[name].note
    assert "DERIVED HERE" in rows["U_FFT_Z3_evolution"].note
    assert "DERIVED HERE" in hr.note and "no FP entry" in hr.note
    assert hr.status is c.CircuitStatus.SCALING and ht.status is c.CircuitStatus.SCALING
    assert rows["mass_hwp_toffolis_evolution"].status is c.CircuitStatus.COMPILED


def test_2028_ch9_rules_on_this_lattice(a, r28):
    """Ch. 9's PRE-R-HOP accounting (2 strings per hop, 30 T per rotation) on 4^2: the record behind the retired 5e3.
    Ch. 9's printed step is read from its own model, and its spatial-hop rule (groups.rhop_link at Ch. 9's own synthesis
    setting) at this chapter's k = 9 and eps = 1e-4 gives the same 8 x HWP(9) as this chapter's hop at that eps."""
    i = r28.intermediates
    assert i["ch9_t_per_step_printed"] == (float(ch9.Assumptions().t_per_step_2028_stated.lo),
                                           float(ch9.Assumptions().t_per_step_2028_stated.hi))
    assert i["ch9_current_rule_hop_strings"] == 8
    syn = ch9.Assumptions().hop_synthesis.value
    expect = {"rus": "t_hop_per_link_2028_full_fit_at_1e-4", "rus-slope": "t_hop_per_link_2028_at_1e-4"}[syn]
    assert i["ch9_current_rule_hop_per_link_k9"] == i[expect]
    assert (i["ch9_rules_hop_group_k"], i["ch9_rules_hop_groups"]) == (9, 64)
    assert (i["ch9_rules_electric_group_k"], i["ch9_rules_electric_groups"]) == (32, 3)
    assert i["ch9_rules_t_per_rot"] == 30
    assert i["ch9_rules_t_hopping_per_step"] == 64 * (4 * 30 + 7 * 7) == 10816
    assert i["ch9_rules_t_electric_per_step"] == 3 * (6 * 30 + 31 * 7) == 1191
    assert i["ch9_rules_t_per_step"] == 12007
    assert i["ch9_rules_dressed_hop_per_step"] == 32 * 8 * (4 * 30 + 7 * 7) == 43264   # with the Z3 dressing
    assert (i["ch9_rules_ancilla_hop_group"], i["ch9_rules_ancilla_max_group"]) == (7, 31)   # k - w(k) (E26; were 11, 37)


def test_2028_gap_target_and_window(r28):
    """app10:99 'a 44x cut at 2028'; app10:88 'Two steps at the chapter's Delta t = 0.1 a reach t = 0.2 a'. Ruling E4,
    restored with the 2 evolution steps: 'one or two evolution steps at Delta t = 0.1 a, 0.2 a and 0.3 a sample
    t = 0.1--0.4 a in steps of 0.1 a and 0.6 a, three times the reach of one run'; 't = 0.2 a, reached at two step sizes
    (2 x 0.1 a and 1 x 0.2 a), checks it'."""
    i = r28.intermediates
    assert i["t_step_gap_target"] == 1e3
    near(i["gap_cut_factor_2028"], 44.06, 0.001)            # '44x'
    assert round(i["gap_cut_factor_2028"]) == 44
    near(i["dt_over_a"], 0.1, 1e-12)
    near(i["t_window_2028_over_a"], 0.2, 1e-12)             # 't = 0.2 a'
    assert i["dt_multi_over_a"] == (0.1, 0.2, 0.3)
    assert all(math.isclose(x, y) for x, y in zip(i["t_reached_multi_over_a"], (0.2, 0.4, 0.6)))
    assert i["t_sampled_multi_over_a"] == (0.1, 0.2, 0.3, 0.4, 0.6)   # '0.1--0.4 a in steps of 0.1 a and 0.6 a'
    assert i["t_overlap_multi_over_a"] == (0.2,)            # the step-size cross-check
    near(i["t_max_multi_over_a"], 0.6, 1e-12)
    near(i["t_max_multi_over_a"] / i["t_window_2028_over_a"], 3.0, 1e-12)   # 'three times the reach of one run'


def test_2028_shots_wall_time_and_fault_budget(r28):
    """Box: 'Shots ~1e3'; 'Wall time ~2 min (~0.13 s/shot at 1 us per T-gate; Ch. 1)'; 'eps_l <~ 8e-7'."""
    i = r28.intermediates
    assert i["shots_2028"] == 1e3 and r28.shots == (1e3, 1e3)
    assert i["shot_overhead_s"] == 1e-4                     # Ch. 1 convention: ~0.1 ms per shot (E27)
    near(i["wall_per_shot_2028_s"], r28.hard_ops[0] * 1e-6 + 1e-4, 1e-12)   # 1 us per T + 0.1 ms
    near(i["wall_per_shot_2028_s"], 0.132275, 1e-5)         # '~0.13 s'
    assert round(i["wall_per_shot_2028_s"], 2) == 0.13
    near(i["wall_2028_s"], 132.275, 1e-5)                   # 1e3 shots in series on one machine
    near(i["wall_2028_s"], i["shots_2028"] * i["wall_per_shot_2028_s"], 1e-12)
    assert round(i["wall_2028_min"]) == 2                   # '~2 min' (2.20)
    near(r28.wall_time_s[0], 132.275, 1e-5)
    near(i["eps_l_required_2028"], 7.5657e-7, 0.001)        # R3: 0.1 / 132,175 -> '8e-7'
    assert round(i["eps_l_required_2028"] * 1e7) == 8
    assert r28.epsilon_l == (i["eps_l_required_2028"],) * 2
    near(i["faults_per_shot_2028_at_rfi_eps_l"], 1.3217e-3, 0.001)


def test_2028_synthesis_budget(a, r28):
    """Ruling R-TOL, app10:88: 'For 2028, N = 3.8e3 gives eps = 1.6e-3 and 19.9 T per rotation' (E20: every rotation)."""
    i = r28.intermediates
    assert i["synth_rotations_per_step_2028"] == 32 * (0.5 + 2 * 2 + 1 + 8 * 4) + 16 * 4 == 1264
    assert i["synth_rotations_per_shot_2028"] == 3 * 1264 == 3792
    near(i["synthesis_error_2028"], EPS_SYN, 1e-12)         # N eps^2 = eps_syn by construction
    near(i["eps_rot"], 1.6239e-3, 1e-4)                     # '1.6e-3'
    near(i["t_per_rotation_papers"], 1.15 * L2, 1e-12)      # legacy record
    assert round(i["t_per_rotation_report"], 1) == 19.9
    near(i["eps_rot_for_budget_2028"], c.eps_per_rotation(0.1, 3792, "incoherent"), 1e-12)
    # the count agrees with the breakdown's rotation rows
    rot_rows = sum(p.count * GROUPS["Z3"].primitives[p.name.split("_Z3_")[0]].rot for p in r28.breakdown
                   if p.name.startswith("U_"))
    rot_rows += sum(p.count for p in r28.breakdown if p.name.endswith("rotations") or "rotations_" in p.name)
    assert rot_rows == 3792


# --- 2033 target: plug-in app10:90, box app10:134-157 ------------------------------------

def test_2033_register_uncut_and_cut(r33):
    """R1 and R11 option (A): '9.64.2 + 3.3.64 + 100 = 1828 LQ; cutting L_par to 8 (V=32) gives the box's 964 LQ,
    at the 1000-LQ marker (1180 LQ and 2.0--3.7e10 T at L_par=10), whereas the uncut lattice is 1.8x the register limit'."""
    i = r33.intermediates
    g = GROUPS["S72x3"]
    assert i["q_G_S72x3"] == g.link_qubits == 9
    assert i["q_G_S72x3_printed_before_R1"] == int(g.link_qubits_chapter.lo) == 8
    assert i["ceil_log2_order_S72x3"] == 8                  # 'the dense floor, but no circuit has been built on it'
    assert (i["V_uncut"], i["lq_gauge_uncut"], i["lq_fermion_uncut"], i["lq_anc_2033"]) == (64, 1152, 576, 100)
    assert i["lq_uncut_2033"] == 1828
    near(i["lq_uncut_over_marker"], 1.83, 0.002)            # '1.8x'
    assert i["V_cut"] == (32, 40)
    assert i["lq_gauge_cut"] == (576, 720)
    assert i["lq_fermion_cut"] == (288, 360)                # app10:30 '288--360 fermion qubits cut, 576 on the uncut'
    assert i["lq_cut_2033"] == (964, 1180)                  # 1180: the L_par=10 parenthetical
    assert i["lq_priced_2033"] == 964 == 9 * 32 * 2 + 3 * 3 * 32 + 100 and r33.lq == (964, 964)
    assert i["lq_priced_2033_over_marker"] == 0.964         # 'at the 1000-LQ marker'
    assert i["lq_cut_2033_over_marker"] == (0.964, 1.18)
    assert i["lq_cut_2033_mid"] == 1072
    assert i["lq_cut_2033_sampling_ancilla_only"] == (872, 1088)   # the 2028 phasing + synthesis 8, not 100 (record)


def test_2033_gauge_terms_come_from_groups(a, r33):
    """R2: nothing is retyped. 'The magnetic term is 9.4e3 T per link per step. ... the electric term is priced at
    that paper's stated lower bound for one, 3.4e4 T per link per step, and is a floor.' (E20 full fit: n_rot = the
    paper's log coefficient / 1.15 at 1.15 log2(1/eps) + 9.2; R-TOL, low ramp end's eps, the R4 window)"""
    i = r33.intermediates
    u = GROUPS["S72x3"].primitives
    E = r33.intermediates["eps_rot"][0]
    near(E, EPS33[0], 1e-12)
    assert i["magnetic_t_per_link_2033"] == magnetic_per_link("S72x3", "KS", 2, E)
    assert i["electric_floor_t_per_link_2033"] == PRIMCOST["KS"]["U_F"](2) * u["U_FFT"].t(E) + u["U_phi"].t(E)
    assert i["electric_dense_t_per_link_2033"] == electric_per_link("S72x3", "KS", 2, E, fft=False)
    assert i["fft_compiled_2033"] is False and u["U_FFT"].status is c.CircuitStatus.CONJECTURE
    assert i["fourier_key_2033"] == "U_FFT"
    F = lambda coef: (coef / 1.15) * (1.15 * L2_33 + 9.2)  # E20: the paper's log term as rotations at the full fit
    near(i["magnetic_t_per_link_2033"], 3 * 637 + 6 * 1120 + 0.5 * (1414 + F(8.05)), 1e-12)
    near(i["electric_floor_t_per_link_2033"], 2 * (532 + F(515.2)) + F(294.4), 1e-12)
    near(i["magnetic_t_per_link_2033"], 9.4372e3, 0.001)    # '9.4e3'
    near(i["electric_floor_t_per_link_2033"], 3.3924e4, 0.001)  # '3.4e4'
    near(i["electric_dense_t_per_link_2033"], 9.2492e6, 0.001)  # the compiled dense transform (not printed)
    near(i["t_per_link_2033"], 4.3362e4, 0.001)             # '4.3--4.4e4' (low end)
    near(i["t_gauge_per_step_2033"], 2.7752e6, 0.001)
    near(i["electric_floor_share_of_step_2033"], 0.2212, 0.002)
    near(i["magnetic_report_over_papers_2033"], 1.0, 1e-12)  # E20: the report price IS the headline
    near(i["electric_report_over_papers_2033"], 1.0, 1e-12)
    near(i["report_convention_t_gauge_per_step_2033"], 2.7752e6, 0.001)
    # the papers' slope-only price, legacy record
    near(m.gauge_step(dataclasses.replace(a, eps_rot=c.Stated(EPS33[0], "app10:88", "test")), "S72x3", 2)["electric_papers"],
         2.3326e4, 0.001)
    # the high ramp end, at its own eps: '4.3--4.4e4'
    hi = dataclasses.replace(a, eps_rot=c.Stated(EPS33[1], "app10:88", "test"))
    near(m.gauge_step(hi, "S72x3", 2)["per_link"], 4.4905e4, 0.001)


def test_2033_hop_from_the_draft(a, r33):
    """r19, app10:95: 'The hop uses the unpublished gate counts ... The color register map and parity flags are computed
    once per link and held through both frames, the hop and both undos. ... the draft draws it but does not count it, so
    its cost is our estimate ... The hop is 1.1e5 T per link (3672 Toffolis, 360 T, 2933 rotations), 8.6e4 T of it the
    color frame; recomputing the map for each frame would give 1.8e5 T'; 'the step is 1.0e7 T, 72% of it the hop'.
    groups.hop_link_cost defaults since r19 (E21 (1)): frame undo on, share='link', MBU ladders, HWP phasing, full fit
    at the circuit's R-TOL eps."""
    from estimates.groups import hop_link_cost
    i = r33.intermediates
    assert a.hop_share.value == "link" and a.hop_share.prov is c.Provenance.STATED
    ref = hop_link_cost("S72x3", 3, eps=EPS33[0], legacy=False)
    assert ref["share"] == "link" and ref["undo"] is True
    assert (i["hop_toffoli_per_link_2033"], i["hop_t_direct_per_link_2033"], i["hop_rotations_per_link_2033"]) \
        == (ref["toffoli"], ref["t_direct"], ref["n_rot"]) == (3672, 360, 2933)
    near(i["t_hop_per_link_2033"], 7 * 3672 + 360 + 2933 * (1.15 * L2_33 + 9.2), 1e-12)
    near(i["t_hop_per_link_2033"], 109726.15, 1e-6)          # '1.1e5'
    near(i["t_hop_per_link_2033_hi_end"], 113642.02, 1e-6)
    it = i["hop_items_2033"]
    assert it["colour_squish_and_parity"][0] == 2 * 1219 + 892 == 3330   # once per link (E21), was 13,320
    assert it["colour_rotations"][2] == 2208 == 12 * 184     # 4 frame applications x 3 fields x 2 N_angles (E21 (3))
    frame = 7 * it["colour_squish_and_parity"][0] + it["colour_rotations"][2] * (1.15 * L2_33 + 9.2)
    near(frame, 8.6292e4, 0.001)                             # '8.6e4 T of it the color frame'
    near(frame / i["t_hop_per_link_2033"], 0.787, 0.002)
    assert set(i["hop_estimated_items_2033"]) == {"colour_squish_and_parity", "colour_rotations", "hop_phasing_hwp"}
    near(i["hop_over_rhop_2033"], 15.74, 0.002)              # vs the retired R-HOP (slope-only)
    v = i["hop_variant_t_per_link_2033"]                     # sensitivities, not the headline
    near(v["share_draft"], 184220.15, 1e-6)                  # '1.8e5 T' (the r17-r18 headline)
    near(v["no_undo"], 7.8235e4, 0.001); near(v["draft_mcx"], 1.1859e5, 0.001)
    near(v["share_draft_no_undo"], 1.0611e5, 0.001)          # the draft's literal 2 C^G, per-application squish
    near(i["t_hop_per_step_2033"], 64 * i["t_hop_per_link_2033"], 1e-12)
    near(i["t_mass_per_step_2033"], 5219.13, 1e-5)          # '5.2--5.3e3' (low end)
    near(i["hopping_share_of_step_2033"], 0.7164, 0.002)     # '72% of it the hop'
    near(i["hop_over_retired_2033"], 20.65, 0.002)
    rows = {p.name: p for p in r33.breakdown}
    assert not any("move" in n for n in rows)                # the U_mul color move is retired
    cs = rows["hop_colour_squish_and_parity_S72x3_evolution_toffoli"]
    assert cs.status is c.CircuitStatus.SCALING and "FermionPrimitives_unpub" in cs.src and "ESTIMATED" in cs.note
    assert "share='link'" in cs.note
    assert (cs.count, cs.t_each) == (3330 * 64 * 300, 7)    # R4: 300 evolution steps at the low end
    # the switch is a single input: share='draft' gives back the r18 shot
    rd = m.model(dataclasses.replace(a, hop_share=c.Stated("draft", "app10:95", "test")), "2033")
    near(rd.hard_ops[0], 7.285238e9, 1e-5); near(rd.hard_ops[1], 3.8435519e10, 1e-5)   # top point at 575 steps (step rule; was 3.73e10)
    assert rd.intermediates["eps_rot"] == i["eps_rot"]       # same rotation count, same R-TOL eps
    assert rows["hop_hop_diagonalizers_S72x3_evolution_rotations"].status is c.CircuitStatus.COMPILED


def test_2028_z3_hop_does_not_move_with_e21(a, r28):
    """E21 (1) moves only the g-only squish/parity Toffolis; the Z3 hop (R-HOP structure, diagonal link) has none, so the
    2028 shot is the same under share='link' and share='draft' (132,174.9 T, 1.32x; left as is pending the author)."""
    rd = m.model(dataclasses.replace(a, hop_share=c.Stated("draft", "app10:95", "test")), "2028")
    assert rd.hard_ops == r28.hard_ops
    near(r28.hard_ops[0], 132174.88, 1e-6)


def test_2033_current_insertion_is_not_per_step(a, r33):
    """app10:95 'Readout (a) injects by quench; readout (b)'s one or two insertions of the bilocal current cost ~2e4 T
    each (derived here: 16 group multiplications ...), below 1e-4 of the shot'. The bilinear's rotations at the full fit."""
    i = r33.intermediates
    assert i["insertion_umul_count"] == 16
    near(i["insertion_t_S72x3"], 16 * 1120 + 2 * (4 * (1.15 * L2_33 + 9.2) + 7), 1e-12)
    near(i["insertion_t_S72x3"], 1.8161e4, 0.001)           # '~2e4'
    assert i["insertion_share_of_shot_b"][1] < 1e-4
    assert not any("insertion" in p.name for p in r33.breakdown)


def test_2033_t_count_chain(r33):
    """r26 (referee J1): 'The 30--50 a window is 300--500 steps, 2.9--5.1e9 T; the adiabatic ramp and the same ramp run
    backwards for the readout map (1e2--1e3 steps each, at the same cost) add 2.0e9--2.0e10 T: 4.9e9--2.5e10 T per shot,
    4.9--25x the 1e9-T resource limit'; 'the step is 1.0e7 T'. The r25 circuit (no readout ramp, 3.9e9--1.5e10) and the
    r23 circuit (t = 100 a) are records. Each ramp end is its own circuit, with its own eps (R-TOL)."""
    i = r33.intermediates
    assert (i["V_priced_2033"], i["n_links_2033"], i["n_plaq_2033"]) == (32, 64, 32)
    near(i["t_per_step_2033"][0], i["t_gauge_per_step_2033"] + i["t_hop_per_step_2033"] + i["t_mass_per_step_2033"], 1e-12)
    near(i["t_per_step_2033"][0], 9.802859e6, 1e-5)         # '1.0e7'
    near(i["t_per_step_2033"][1], 1.015238e7, 1e-5)
    assert i["n_trot_2033"] == 1e3                          # the r23 window, record
    assert i["window_over_a"] == (30, 50) and i["n_evo_2033"] == (300, 500)
    assert i["n_ramps_2033"] == 2 and i["n_steps_2033"] == (500, 2500)       # r26: prep + backward readout ramp
    assert i["n_steps_2033_r25"] == (400, 1500)
    assert i["window_fm"] == (3.0, 5.0) and i["window_over_L_par"] == (3.75, 6.25)
    near(i["t_evolution_2033"][0], 2.940858e9, 1e-5)        # '2.9e9'
    near(i["t_evolution_2033"][1], 5.076190e9, 1e-5)        # '5.1e9'
    assert i["ramp_steps_2033"] == (1e2, 1e3)
    near(i["t_ramp_2033"][0], 9.802859e8, 1e-5)             # one ramp
    near(i["t_ramp_2033"][1], 1.015238e10, 1e-5)
    near(i["t_ramps_2033"][0], 1.960572e9, 1e-5)            # '2.0e9' (both ramps)
    near(i["t_ramps_2033"][1], 2.030476e10, 1e-5)           # '2.0e10'
    near(i["t_total_2033"][0], 500 * i["t_per_step_2033"][0], 1e-12)
    near(i["t_total_2033"][1], 2500 * i["t_per_step_2033"][1], 1e-12)
    near(i["t_total_2033"][0], 4.9014296e9, 1e-5)           # '4.9e9'
    near(i["t_total_2033"][1], 2.5380950e10, 1e-5)          # the 500-step shot (1.55 and 3.1 GeV points)
    # step rule (2026-10-04): the 4.65 GeV point at 50 a runs 575 steps; hard_ops = (low end, that shot)
    assert i["step_rule_top_mode"]["n"] == 575 and i["step_rule_top_mode"]["phase"] <= 0.1
    near(i["t_shot_top_mode_2033"], 2.6158908e10, 1e-5)     # '2.6e10'
    assert r33.hard_ops == (i["t_total_2033"][0], i["t_shot_top_mode_2033"])
    near(i["t_total_2033_over_cap"][0], 4.901, 0.001)       # '4.9--26x'
    near(r33.hard_ops[1] / 1e9, 26.16, 0.001)
    near(i["t_total_2033_r25"][0], 3.9017597e9, 1e-5)       # the r25 headline, a record
    near(i["t_total_2033_r25"][1], 1.5062166e10, 1e-5)
    near(i["t_total_2033_r23"][0], 1.0971497e10, 1e-5)      # the r23 headline, a record (no longer printed)
    near(i["t_total_2033_r23"][1], 2.0207840e10, 1e-5)
    near(i["pre_r17_t_total_2033"][0], 2.76781e9, 1e-5)     # the R-TOL print (r23 window), record
    near(i["pre_r17_t_total_2033"][1], 5.10638e9, 1e-5)
    near(i["t_total_2033_round_e"][0], 2.15437e9, 1e-5)     # round E's fixed 1e-3, record
    near(i["t_total_2033_round_e"][1], 3.91704e9, 1e-5)
    near(i["applied_stage_t_per_step_2033"], 3.1152e6, 0.001)   # the pre-round-D hop figure on this gauge step
    near(i["t_total_2033_qubitized"][0], 4.9014e8, 0.001)
    assert c.tex_hard_ops(*r33.hard_ops) == r"$4.9{\times}10^{9}$--$2.6{\times}10^{10}$"   # Table 1.1 spelling
    # the same rules at V = 40 (L_par = 10, 1180 LQ), its own N_rot and eps: a sensitivity, not printed since r26
    near(i["t_total_2033_at_V_hi"][0], 6.15707e9, 1e-4)
    near(i["t_total_2033_at_V_hi"][1], 3.18776e10, 1e-4)


def test_2033_readout_map_and_source_r26(r33):
    """Referee J1/J2 (r26): the backward ramp is in the shot; the electric-basis Fourier layer (< 1e-3 of the shot) and the
    gauge-invariant source (~1% of a step while on) are recorded, not added; the 0.8 fm axis recollides at 0.4 fm/c and
    separation through 1/Lambda_QCD needs L_par = 20 sites, 2260 LQ."""
    i = r33.intermediates
    near(i["t_readout_fourier_2033"], 64 * GROUPS["S72x3"].primitives["U_FFT"].t(i["eps_rot"][0]), 1e-12)
    assert i["t_readout_fourier_share"] < 1e-3
    near(i["t_source_per_step_2033"], i["t_hop_per_link_2033"], 1e-12)
    assert 0.01 < i["t_source_share_of_step"] < 0.012                        # r26 v1 record
    # J2 v2: the wave-packet source, L_par separations at the full 8-link bilocal price: '< 2% of a step'
    near(i["t_source_wavepacket_2033"], 8 * i["insertion_t_S72x3"], 1e-12)
    assert 0.014 < i["t_source_wavepacket_share_of_step"] < 0.02
    assert i["recollision_time_fm"] == 0.4 and i["L_par_for_separation"] == 20 and i["lq_separation_volume"] == 2260
    assert 9 * 80 * 2 + 9 * 80 + 100 == 2260


def test_2033_synthesis_budget_is_reported(a, r33, r3d):
    """Ruling R-TOL, app10:88: 'for 2033, N = 1.3--6.5e8 gives 8.7--3.9e-6 and 28.5--29.9 T; for 3+1D, 7.9e8 gives
    3.6e-6 and 30.0 T' (E20: one price for every rotation; R4 window; r26 readout ramp). At round E's 1e-3 N eps^2 would
    be 131--654 (record)."""
    i = r33.intermediates
    assert i["synth_rotations_per_step_2033"] == 261792 == 64 * (i["hop_rotations_per_link_2033"] + 1155.5) + 32 * 4
    assert i["pre_r17_n_rot_per_step_2033"] == 74592
    assert i["synth_rotations_per_shot_2033"] == N33
    near(i["synthesis_error_2033"][0], EPS_SYN, 1e-12)
    near(i["synthesis_error_2033"][1], EPS_SYN, 1e-12)
    near(i["synthesis_error_2033_round_e"][0], 130.896, 1e-4)
    near(i["synthesis_error_2033_round_e"][1], 654.48, 1e-4)
    near(i["eps_rot_for_budget_2033"][0], 1.2361e-5, 0.001)  # record: the whole 0.1
    near(i["eps_rot_for_budget_2033"][1], 2.7640e-5, 0.001)
    assert [round(x * 1e6, 1) for x in i["eps_rot"]] == [8.7, 3.9]   # '8.7--3.9e-6' (r26; r25 9.8--5.0e-6)
    assert [round(x, 1) for x in i["t_per_rotation_papers"]] == [19.3, 20.7]
    assert [round(x, 1) for x in i["t_per_rotation_report"]] == [28.5, 29.9]
    assert round(i["synth_rotations_per_shot_2033"][0] / 1e8, 1) == 1.3 and round(N33[1] / 1e8, 1) == 6.5
    # the total synthesis error is a knob: a tighter budget tightens every rotation
    b = dataclasses.replace(a, eps_syn=c.Stated(1e-4, "app10:88", "test"))
    near(m.model(b, "2033").intermediates["eps_rot"][1], math.sqrt(1e-4 / N33[1]), 1e-12)
    assert m.model(b, "2033").hard_ops[0] > r33.hard_ops[0]
    assert r3d.intermediates["synth_rotations_per_shot_3d"] == N3D
    near(r3d.intermediates["synthesis_error_3d"], EPS_SYN, 1e-12)
    near(r3d.intermediates["synthesis_error_3d_round_e"], 785.92, 1e-6)
    assert round(r3d.intermediates["eps_rot"] * 1e6, 1) == 3.6 and round(r3d.intermediates["t_per_rotation_papers"], 1) == 20.8
    assert round(r3d.intermediates["t_per_rotation_report"], 1) == 30.0 and round(N3D / 1e8, 1) == 7.9


def test_2033_breakdown_sums_and_statuses(r33):
    rows = {p.name: p for p in r33.breakdown}
    tags = ("evolution", "ramp", "readout_ramp")                 # r26 (J1): the backward ramp of the readout map
    gauge = {f"{p}_S72x3_{tag}" for p in ("U_inv", "U_mul", "U_Tr", "U_FFT", "U_phi") for tag in tags}
    hop = {f"hop_{it}_S72x3_{tag}_{k}" for tag in tags for it, k in
           (("colour_squish_and_parity", "toffoli"), ("colour_rotations", "rotations"), ("hop_squish_and_flags", "toffoli"),
            ("hop_diagonalizers", "t"), ("hop_diagonalizers", "rotations"), ("hop_phasing_hwp", "toffoli"),
            ("hop_phasing_hwp", "rotations"))}
    mass = {f"{p}_{tag}" for p in ("mass_hwp_rotations", "mass_hwp_toffolis") for tag in tags}
    assert set(rows) == gauge | hop | mass
    assert (rows["U_mul_S72x3_evolution"].count, rows["U_mul_S72x3_evolution"].t_each) == (6 * 64 * 300, 1120)
    assert (rows["U_inv_S72x3_ramp"].count, rows["U_inv_S72x3_ramp"].t_each) == (3 * 64 * 1e2, 637)
    assert rows["U_inv_S72x3_readout_ramp"].count == rows["U_inv_S72x3_ramp"].count
    near(r33.breakdown_total(), r33.hard_ops[0], 1e-12)     # the low end of the ramp range
    for name, p in rows.items():
        assert p.status is not c.CircuitStatus.UNSOURCED
        if name.startswith("U_FFT"):
            assert p.status is c.CircuitStatus.CONJECTURE and "FLOOR" in p.note
        elif name.startswith("U_"):
            assert p.status is c.CircuitStatus.COMPILED
        elif name.startswith("hop_"):
            assert "FermionPrimitives_unpub" in p.src or "arxiv_1709_06648" in p.src


def test_2033_ramp_may_be_a_product_state_start(a):
    """A zero-step ramp is a legitimate input (product-state start): evolution only."""
    b = dataclasses.replace(a, ramp_steps_2033=c.Stated((0, 0), "app10:90", "test"))
    rb = m.model(b, "2033")
    near(rb.hard_ops[0], rb.intermediates["t_evolution_2033"][0], 1e-12)
    near(rb.intermediates["t_total_2033"][1], rb.intermediates["t_evolution_2033"][1], 1e-12)   # 300 and 500 steps (R4 window ends)
    near(rb.hard_ops[1], rb.intermediates["t_shot_top_mode_2033"], 1e-12)  # step rule: the 575-step top point
    assert not any(p.name.endswith("_ramp") for p in rb.breakdown)
    near(rb.breakdown_total(), rb.hard_ops[0], 1e-12)


def test_2033_kinematics(a, r33):
    """R4: 'the window is t_had = 30--50 a (3--5 fm/c at a ~ 0.1 fm), three to five traversals of the 0.8 fm periodic
    box'; 'N_Trot ~ t_had/Delta t ~ 1e2--1e3'; 'p ~ 2 GeV'. The step Delta t = 0.1 a is the r23 t_had / N_Trot."""
    i = r33.intermediates
    near(i["dt_over_a_2033"], 0.1, 1e-12)
    near(i["t_had_fm"], 10.0, 1e-12)
    assert i["n_trot_eq_range"] == (1e2, 1e3)
    near(i["momentum_unit_GeV_at_Lpar_8_10_16"][0], 1.55, 0.002)
    near(i["momentum_unit_GeV_at_Lpar_8_10_16"][1], 1.24, 0.002)
    near(i["source_momentum_in_units_Lpar_8"], 1.29, 0.002)
    near(i["source_momentum_times_a"], 1.01, 0.005)
    # J2 v2: box modes 1.55, 3.1, 4.65 GeV; the top one at pa ~ 2.4
    bm = i["source_momenta_box_modes_GeV"]
    near(bm[0], 1.55, 0.002); near(bm[1], 3.10, 0.002); near(bm[2], 4.65, 0.002)
    near(i["top_box_mode_times_a"], 2.36, 0.002)
    assert i["t_had_over_L_par_cut"] == (10.0, 12.5)
    assert 300 >= i["n_trot_eq_range"][0] and 500 <= i["n_trot_eq_range"][1]   # the R4 steps lie in the eq. range
    # workflow step 4 used to say 't ~ 1/Lambda_QCD'; the r23 window was ten times that (record)
    near(i["inverse_lambda_qcd_fm"][0], 0.66, 0.005)
    near(i["inverse_lambda_qcd_fm"][1], 0.99, 0.005)
    assert 10 <= i["t_had_over_inverse_lambda_qcd"][0] <= i["t_had_over_inverse_lambda_qcd"][1] <= 15.3
    # the step length is derived from two printed inputs; changing the window changes it, visibly
    b = dataclasses.replace(a, t_had_over_a=c.Stated(50, "app10:74", "test"))
    near(m.model(b, "2033").intermediates["dt_over_a_2033"], 0.05, 1e-12)


def test_2033_logical_faults(r33):
    """R3: 'At eps_l = 1e-8 the 2+1D shot incurs 49--250 logical faults ...; the report's criterion of 0.1 expected
    faults per shot needs eps_l <~ 3.9e-12--2.0e-11' (r26: with the readout ramp; r25 39--150, 6.6e-12--2.6e-11)."""
    i = r33.intermediates
    near(i["faults_per_shot_2033_at_rfi_eps_l"][0], 49.014, 0.001)   # '49' (r26)
    near(i["faults_per_shot_2033_at_rfi_eps_l"][1], 261.59, 0.001)   # '260' (step rule; r26 250)
    near(i["eps_l_required_2033"][0], 3.8228e-12, 0.001)    # '3.8e-12' (step rule; r26 3.9e-12)
    near(i["eps_l_required_2033"][1], 2.0402e-11, 0.001)    # '2.0e-11' (r26)
    assert r33.epsilon_l == i["eps_l_required_2033"]


def test_2033_sampling_shots_and_wall_time(a, r33):
    """R2 eq:Nshot_samp: N = (c / n_s + 1 / n_M) / delta^2, c = 1--2, n_M = 0.5--1: 'a 30% ratio needs 120--240 shots
    and a 10% ratio 1.1--2.2e3'. R10/R9: per shot max(N_T t_gate, D_T t_r) + t0 = '4.9e3--3.6e4 s' (requirements; r26), up
    to '1.43x' the T-count time. R1 box: 'First result ... one configuration at 30%: 120--240 shots, 6.9--103 d' (r26)."""
    i = r33.intermediates
    near(i["ratio_shots_at_10pct"][0], 1.1e3, 1e-9); near(i["ratio_shots_at_10pct"][1], 2.2e3, 1e-9)
    near(i["shots_first_2033"][0], 122.22, 1e-4); near(i["shots_first_2033"][1], 244.44, 1e-4)   # '120--240'
    assert i["delta_first_2033"] == 0.30
    near(i["wall_per_shot_serial_2033_s"][0], r33.hard_ops[0] * 1e-6 + 1e-4, 1e-12)
    near(i["wall_per_shot_serial_2033_s"][1], i["t_total_2033"][1] * 1e-6 + 1e-4, 1e-12)
    near(i["wall_shot_top_mode_2033_s"], 37480.6, 1e-5)    # '3.7e4' (the 575-step top point, step rule)
    near(i["wall_per_shot_2033_s"][0], 4901.43, 1e-5)       # '4.9e3' (F* >= 10 at the low end: the serial time)
    near(i["wall_per_shot_2033_s"][1], 36364.2, 1e-5)      # '3.6e4' (D_T x 10 us + 0.1 ms)
    near(i["wall_per_shot_2033_s"][1], 36364.2, 1e-5)      # the 500-step circuit, depth-bound
    near(i["depth_penalty_2033"][1], 1.4327, 1e-4)          # 'up to 1.43x'
    assert i["depth_penalty_2033"][0] == 1.0
    for k in (0, 1):
        near(i["wall_first_2033_s"][k], i["shots_first_2033"][k] * i["wall_per_shot_2033_s"][k], 1e-12)
    near(i["wall_first_2033_days"][0], 6.934, 1e-3)         # '6.9 d'
    near(i["wall_first_2033_days"][1], 102.9, 1e-3)         # '103 d'
    assert i["wall_first_result_s"] == i["wall_first_2033_s"]
    # the r23 accounting is kept as the record
    near(i["shots_sampling_from_eq"][0], 1e3, 1e-9); near(i["shots_sampling_from_eq"][1], 4e3, 1e-9)
    assert i["shots_sampling_printed"] == (1e3, 4e3)
    near(i["wall_sampling_s_r23"][0], 1.0971e7, 0.001); near(i["wall_sampling_s_r23"][1], 8.0831e7, 0.001)


def test_2033_sampling_campaign_on_one_machine(a, r33):
    """R2 campaign: 'a 25% trend shows at 1.8 sigma with 10% points and at 3 sigma with 5.9% points (3.2--6.3e3 shots
    each)'; gap: 'The 3 sigma trend over three source momenta at one spacing, with the L_par = 4 box check, takes 1.7--26 yr;
    at the high end five years hold 0.68 of the three trend configurations and the campaign needs a 5.2x cut. A second
    spacing on the same lattice ... triples the trend; at fixed volume it needs 3.6e3 LQ.' Box: '1.3e4--2.5e4 shots,
    1.7--26 yr'. Smaller ruling (a), 2026-10-05: the check is in the campaign (the trend alone: 9.5e3--1.9e4, 1.5--22 yr)."""
    i = r33.intermediates
    yr = 3.156e7
    near(i["trend_delta_per_point"], 0.25 / (3 * math.sqrt(2)), 1e-12)
    near(i["trend_delta_per_point"], 0.0589, 0.002)         # '5.9%'
    near(i["trend_sigma_at_delta"][0.10], 1.768, 1e-3)      # '1.8 sigma'
    near(i["trend_sigma_at_delta"][0.05], 3.536, 1e-3)
    near(i["shots_per_configuration_campaign"][0], 3168, 1e-9); near(i["shots_per_configuration_campaign"][1], 6336, 1e-9)
    near(i["shots_trend_2033"][0], 9504, 1e-9); near(i["shots_trend_2033"][1], 19008, 1e-9)   # the trend: 9.5e3--1.9e4
    assert i["n_sampling_configurations"] == 3 and i["campaign_horizon_yr"] == 5
    near(i["wall_trend_2033_s"][0], i["shots_trend_2033"][0] * i["wall_per_shot_2033_s"][0], 1e-12)
    # step rule: one of the three momenta (4.65 GeV) runs its 575-step shot at the high end
    near(i["wall_trend_2033_s"][1], i["shots_trend_2033"][1] * (2 * i["wall_per_shot_2033_s"][1]
                                                                  + i["wall_shot_top_mode_2033_s"]) / 3, 1e-12)
    near(i["wall_trend_2033_yr"][0], 1.476, 1e-4)          # the trend alone, 1.5--22 yr (r26/r27 campaign)
    near(i["wall_trend_2033_yr"][1], 22.126, 1e-4)
    near(i["trend_over_horizon"][1], 4.425, 1e-3)
    # smaller ruling (a), 2026-10-05: campaign = trend + one configuration at L_par = 4 at the trend precision
    assert i["box_check_in_campaign"] is True
    for k in (0, 1):
        near(i["shots_campaign_2033"][k], i["shots_trend_2033"][k] + i["shots_per_configuration_campaign"][k], 1e-12)
        near(i["wall_campaign_2033_s"][k], i["wall_trend_2033_s"][k] + i["wall_box_check_one_mode_s"][k], 1e-12)
    near(i["shots_campaign_2033"][0], 12672, 1e-9); near(i["shots_campaign_2033"][1], 25344, 1e-9)   # '1.3e4--2.5e4'
    assert (round(i["shots_campaign_2033"][0] / 1e4, 1), round(i["shots_campaign_2033"][1] / 1e4, 1)) == (1.3, 2.5)
    near(i["wall_campaign_2033_yr"][0], 1.7180, 1e-3)       # '1.7'
    near(i["wall_campaign_2033_yr"][1], 25.817, 1e-3)       # '26'
    assert (round(i["wall_campaign_2033_yr"][0], 1), round(i["wall_campaign_2033_yr"][1])) == (1.7, 26)
    near(i["campaign_over_horizon"][1], 5.163, 1e-3)        # '5.2x cut'
    assert round(i["campaign_over_horizon"][1], 1) == 5.2
    near(i["configurations_in_horizon"][0], 0.6779, 1e-3)    # "0.68 of the three" trend configurations
    near(i["wall_campaign_two_spacings_yr"][0], 3 * i["wall_trend_2033_yr"][0], 1e-12)   # 'triples the trend'
    assert i["lq_second_spacing_fixed_volume"] == 3556      # '3.6e3 LQ'
    assert r33.shots == i["shots_campaign_2033"] and r33.wall_time_s == i["wall_campaign_2033_s"]
    assert i["wall_campaign_s"] == i["wall_campaign_2033_s"]
    # the r23 six-configuration record
    near(i["campaign_sampling_yr_at_10pct_r23"][0], 2.0858, 1e-4); near(i["campaign_sampling_yr_at_5pct_r23"][1], 15.367, 1e-4)


def test_2033_depth_exports(a, r33):
    """R9/R10 (factory.json Ch. 6 2033, schedule P): 'a step is 5.3e5--1.5e6 layers deep and the T-count is 7.0--18.5x
    the T-depth'. Exports from common.depth_exports over the campaign; f_star cross-pairs the two ramp ends (1.8--71),
    the per-circuit band is 7.0--18.5."""
    i = r33.intermediates
    d = i["t_depth_per_step_2033"]
    near(d["lo_circuit_lo_est"], 530108.22, 1e-6)           # '5.3e5'
    near(d["hi_circuit_hi_est"], 1454568.03, 1e-6)          # '1.5e6' (r26; r25 1.4e6)
    near(i["f_star_per_circuit_2033"][0][0], 6.9989, 1e-4); near(i["f_star_per_circuit_2033"][0][1], 18.493, 1e-4)
    near(i["f_star_per_circuit_2033"][1][0], 6.9797, 1e-4)
    near(i["t_depth_per_shot"][0], 500 * d["lo_circuit_lo_est"], 1e-12)
    near(i["t_depth_per_shot"][1], i["d_campaign_mean_2033"][1], 1e-12)   # step rule: campaign mean (one 575-step point)
    assert i["t_depth_per_shot"][1] > 2500 * d["hi_circuit_hi_est"]
    dx = c.depth_exports(i["t_campaign_mean_2033"], i["d_campaign_mean_2033"], i["shots_trend_2033"], a.t_gate_s,
                         a.shot_overhead_s)   # step rule: campaign shot-averages at the high end (one 575-step point)
    for k in ("t_per_shot", "t_depth_per_shot", "f_star", "floor_wall_s", "factories_for_1yr"):
        assert i[k] == dx[k], k
    assert i["baseline_ok"] == (False, True)
    near(i["factories_for_1yr"][0], 14.761, 1e-4); near(i["factories_for_1yr"][1], 154.44, 1e-4)   # step rule (was 152.88)
    near(i["floor_wall_s"][1] + 19008 * 1e-4, i["wall_trend_2033_s"][1], 1e-9)   # the trend's high end at the depth floor
    near(i["wall_trend_2033_s"][0], i["wall_serial_s"][0], 1e-12)                 # the low end at the T-count


def test_2033_fragmentation_shots_and_wall_time(r33):
    """Box: '7.5e4--1.2e6 shots, 12--1400 yr (2.3--280x over 5 yr)'. Gap: 'The fragmentation instance takes 12--1400 yr:
    five years hold 4.3e3--3.2e4 of its shots, so it needs a 2.3--280x cut in shots times per-shot T' (r26: the readout
    ramp is in the shot; the N_h insertion of eq:FF_def is not priced, so these are lower bounds, J3)."""
    i = r33.intermediates
    near(i["a_inv2"][0], 25, 1e-9)
    near(i["a_inv2"][1], 400, 1e-9)
    near(i["shots_ff_base"], 3e3, 1e-9)
    near(i["shots_ff_from_eq"][0], 7.5e4, 1e-9)
    near(i["shots_ff_from_eq"][1], 1.2e6, 1e-9)
    assert i["shots_ff_printed"] == (7.5e4, 1.2e6)
    near(i["wall_ff_yr"][0], 11.648, 1e-4)                  # '12' (r26)
    near(i["wall_ff_yr"][1], 1382.7, 1e-4)                  # '1400' (r26)
    for k in (0, 1):
        near(i["wall_ff_s"][k], i["shots_ff_from_eq"][k] * i["wall_per_shot_2033_s"][k], 1e-12)
        near(i["ff_cost_reduction_for_horizon"][k], i["wall_ff_yr"][k] / i["campaign_horizon_yr"], 1e-12)
    near(i["ff_cost_reduction_for_horizon"][0], 2.3296, 1e-4)   # '2.3' (r26)
    near(i["ff_cost_reduction_for_horizon"][1], 276.53, 1e-4)   # '280' (r26)
    near(i["ff_shots_in_horizon"][0], 4339.43, 1e-5)        # '4.3e3' (r26)
    near(i["ff_shots_in_horizon"][1], 32194.7, 1e-5)        # '3.2e4' (r26)
    assert i["ff_shots_in_horizon"][1] < i["shots_ff_from_eq"][0]   # five years do not hold even the low-end budget
    near(i["ff_shots_per_call"][0], 625, 1e-9); near(i["ff_shots_per_call"][1], 1e4, 1e-9)
    near(i["ff_calls_in_horizon"][0], 0.43394, 1e-4)
    near(i["ff_calls_in_horizon"][1], 51.512, 1e-4)


def test_2033_dipole_rows(a, r33):
    """R7 dipole rows; q1 (rulings A/B, 2026-10-05): the media are quench-prepared (20-step coupling ramp + c/T, c = 2-2 pi,
    90-282 steps before the sources enter), each circuit priced whole with its own R-TOL eps (scratchpad/quench/cost/price2.py).
    Chapter: '995 LQ and 5.4e9--1.5e10 T per shot (2.3--6.7e9 of it the dipole evolution), 5.4--15x the 1e9-T limit
    (eps_l <~ 6.6e-12)'; '912 LQ and 2.0--5.5e8 T (eps_l <~ 1.8e-10), under the limit'; '653 shots per temperature for Im V
    and 405 for q-hat: 57--98 and 1.2--2.3 d, and 0.47--0.80 yr and 3.7--6.8 d for three temperatures'."""
    i = r33.intermediates
    ds, dl = i["dipole_static"], i["dipole_lightlike"]
    assert i["dipole_steps"] == (50, 100)                                   # the 0.01 fm/c grid
    # step rule (2026-10-04): thermal rows take the literal state-dependent estimate, each shot its own steps
    # s648 (2026-10-05): the static row is Sigma(216x3) / H_I; its rule runs on the improved free dispersion
    assert i["dipole_steps_rule"] == {"static": (85, 240), "light": (64, 181)}   # static was (72, 203) at S72x3 / H_KS
    assert all(x <= 0.1 for v in i["dipole_sd_err_rule"].values() for x in v)
    assert i["dipole_static_group"] == ("S216x3", "I")
    assert ds["lq"] == 995 == 11 * 81 + 100 + 4 and dl["lq"] == 912 == 6 * 135 + 100 + 2   # was 833 = 9 x 81 + 104
    assert i["s648_static_lq_if_kept_S72x3"] == 833 and i["s648_species_lq_if_switched"] == 1092 == 11 * 64 + 288 + 100
    assert GROUPS["2O"].link_qubits == 6
    # q1: the quench preparation, 1 a coupling ramp + c/T at 0.01 fm/c, aT = 0.446
    import math as _m
    assert ds["n_ramp_steps"] == 20 == _m.ceil(0.2 / 0.01 - 1e-9)
    assert ds["n_thermalization_steps"] == (90, 282) == tuple(_m.ceil(c_ * (0.1973269804 / 0.44) / 0.01) for c_ in (2.0, 2 * _m.pi))
    assert ds["n_prep_steps"] == (110, 302) == dl["n_prep_steps"]
    # verifier q1-1 (2026-10-05): the prep leg at the row's rule step instead of the grid, record only, unpriced
    assert ds["n_prep_steps_at_rule_step"] == (187, 725) and dl["n_prep_steps_at_rule_step"] == (141, 547)
    near(ds["t_shot_prep_at_rule_step"][0], 7.55e9, 1e-2); near(ds["t_shot_prep_at_rule_step"][1], 2.69e10, 1e-2)   # 7.6-27x, unpriced
    near(dl["t_shot_prep_at_rule_step"][0], 2.3e8, 2e-2); near(dl["t_shot_prep_at_rule_step"][1], 8.3e8, 2e-2)       # still fits
    assert a.quench_c_range.value == (2.0, 2 * _m.pi) and a.quench_c_range.prov is c.Provenance.CITED
    assert a.eth_assumed.prov is c.Provenance.ASSUMED and a.quench_ramp_a.value == 1.0
    assert not hasattr(a, "thermal_prep_t") and not hasattr(a, "eroq_shot_penalty")   # retired with E-rho-OQ (ruling A)
    # per-shot T spans (c = 2, 0.5 fm/c) .. (c = 2 pi, 1 fm/c); the dipole evolution alone is the pre-q1 record
    near(ds["t_shot"][0], 5.4040e9, 1e-3); near(ds["t_shot"][1], 1.5092e10, 1e-3)      # '5.4e9--1.5e10'
    near(dl["t_shot"][0], 1.9650e8, 1e-3); near(dl["t_shot"][1], 5.5030e8, 1e-3)       # '2.0--5.5e8'
    near(ds["t_shot_by_c"][0][1], 9.726e9, 1e-3); near(ds["t_shot_by_c"][1][0], 1.076e10, 1e-3)
    near(ds["t_shot_evolution_only"][0], 2.346627e9, 1e-5); near(ds["t_shot_evolution_only"][1], 6.657600e9, 1e-5)   # '2.3--6.7e9 of it'
    near(dl["t_shot_evolution_only"][0], 7.164794e7, 1e-5); near(dl["t_shot_evolution_only"][1], 2.044818e8, 1e-5)
    near(ds["t_shot_over_cap"][0], 5.404, 1e-3); near(ds["t_shot_max_over_cap"], 15.09, 1e-3)   # '5.4--15x the 10^9-T limit'
    assert dl["t_shot_max_over_cap"] < 1 and ds["lq_over_marker"] < 1                 # q-hat row 'under the limit'
    near(ds["eps_l_required"], 6.63e-12, 2e-3); near(dl["eps_l_required"], 1.82e-10, 2e-3)   # '6.6e-12', '1.8e-10'
    # shots (unchanged by q1): binomial on the color-floor P_s curve, rate corrected by the color model (r26 v2, J4)
    near(ds["shots_per_T_binomial"], 653.14, 1e-4)          # '653'
    near(dl["shots_per_T_binomial"], 405.00, 1e-4)          # '405'
    near(ds["shots_per_T_W2"], 614.74, 1e-4); near(dl["shots_per_T_W2"], 249.40, 1e-4)   # r26 v1 record (|W|^2 shape)
    near(dl["shots_qhat_if_adjoint"], 695.05, 1e-4); near(dl["shots_qhat_if_adjoint_W2"], 587.24, 1e-4)   # C_F/C_A = 3/8
    near(_m.log(ds["p_s"][0] / ds["p_s"][1]) / _m.log(ds["p"][0] / ds["p"][1]), ds["floor_slope_ratio"], 1e-9)
    near(_m.log(dl["p_s"][0] / dl["p_s"][1]) / _m.log(dl["p"][0] / dl["p"][1]), dl["floor_slope_ratio"], 1e-9)
    near(ds["floor_slope_ratio"], 0.9608, 1e-3); near(dl["floor_slope_ratio"], 0.6902, 1e-3)   # '0.96 and 0.69'
    assert a.thermal_variance_factor.value == (1, 1)
    assert ds["shots_per_T"] == (ds["shots_per_T_binomial"],) * 2 and dl["shots_per_T"] == (dl["shots_per_T_binomial"],) * 2
    # walls: one temperature at c = 2 .. 2 pi, then three temperatures; each shot at its own time (mean of the two circuits)
    near(ds["wall_first_s"][0] / 86400, 57.19, 1e-3); near(ds["wall_first_s"][1] / 86400, 97.71, 1e-3)   # '57--98 d'
    near(dl["wall_first_s"][0] / 86400, 1.234, 1e-3); near(dl["wall_first_s"][1] / 86400, 2.265, 1e-3)   # '1.2--2.3 d'
    near(ds["wall_scan_s"][0] / 3.156e7, 0.4697, 2e-3); near(ds["wall_scan_s"][1] / 3.156e7, 0.8025, 2e-3)   # '0.47--0.80 yr'
    near(dl["wall_scan_s"][0] / 86400, 3.70, 2e-3); near(dl["wall_scan_s"][1] / 86400, 6.79, 2e-3)           # '3.7--6.8 d'
    per = sum(0.5 * (t * 1e-6 + 1e-4) for t in ds["t_shot_by_c"][0])
    near(ds["wall_first_s"][0], ds["shots_per_T"][0] * per, 1e-12)
    rows = m.INSTANCE_ROWS(a, "2033", r33)
    assert [r[1] for r in rows] == [(964, 964), (995, 995), (912, 912)]
    assert [tuple(round(x / 1e8) for x in r[2]) for r in rows[1:]] == [(54, 151), (2, 6)]
    assert "quench-prepared" in rows[1][3]["note"] and "quench-prepared" in rows[2][3]["note"]


def test_dipole_markov_colour_floor_derivation():
    """J4 (r26), standard library only. In a color-covariant Markovian (Lindblad) description with Hermitian color kicks,
    the q qbar color space is 1 + 8 (1 + 3 for SU(2)), so P_s relaxes as 1/N^2 + (1 - 1/N^2) e^{-Gamma t}, while the
    averaged singlet amplitude decays at |Im V|; equal initial slopes give Gamma = 2|Im V| N^2 / (N^2 - 1). Checked
    numerically against the full 9x9 (4x4) Lindbladian at source correlations 0, 0.5, 0.9 in
    scratchpad/referee/ch06/j4_lindblad.py. Here: the large-r case, P_s from one-site depolarizing at rate
    1/(2N) + C_F per unit kick strength (Fierz), |W|^2 at 2 C_F."""
    for n in (2, 3):
        cf = (n * n - 1) / (2 * n)
        gamma_ps = 2 * (0.5 / n + cf)          # both sources depolarize
        assert math.isclose(gamma_ps / (2 * cf), n * n / (n * n - 1), rel_tol=1e-12)


def test_no_machine_count_anywhere(a, r28, r33, r3d):
    """E27 (H. Lamm 2026-10-02): the report never assumes more than one machine. No model field is a machine count,
    no note quotes one, and the chapter's only 'parallel' are L_parallel and the parallel Hadamard test."""
    import re
    assert not any("machine" in f.name for f in dataclasses.fields(a))
    for r in (r28, r33, r3d):
        assert not any("machines" in k or "n_machine" in k or "parallel" in k for k in r.intermediates)
        assert not any(re.search(r"\d\s*machines|machines for|across machines|shot-parallel|machine-years", n) for n in r.notes)
    if TEX.exists():
        text = TEX.read_text()
        assert "machines" not in text and "serialized" not in text
        assert not re.search(r"shot-parallel|machine-year|embarrassingly|fold parallel", text)
        assert text.count("one machine") == 1                  # the campaign box row (the gap paragraph's restatement cut 2026-10-06)
        rest = text.replace(r"L_\parallel", "").replace("parallel Hadamard test", "")
        assert "parallel" not in rest.lower()


def test_requirements_table(r28, r33):
    """app10:165-176."""
    i = r33.intermediates
    assert i["n_sampling_configurations"] == 3              # '3 sampling configurations (3 source momenta, one spacing)'
    assert i["n_ff_calls"] == 120
    assert i["n_subroutine_calls"] == 129 == 3 + 120 + 2 * 3   # '~130'
    near(i["n_subroutine_calls"], 130, 0.01)
    assert round(r28.intermediates["wall_per_shot_2028_s"], 2) == 0.13   # '~0.13 s (2028)'
    near(r28.hard_ops[0], 1.3e5, 0.02)                      # '1.3e5 T (2028)'
    near(r33.hard_ops[0], 4.9e9, 0.01); near(r33.hard_ops[1], 2.6e10, 0.02)   # '4.9e9--2.6e10 T (2033)' (step rule)
    assert round(i["wall_per_shot_2033_s"][0], -2) == 4900 and round(i["wall_shot_top_mode_2033_s"], -3) == 37000
    assert r33.lq == (964, 964) and r28.lq == (230, 230)    # 'Size 230 LQ (2028) / 964 LQ (2033; dipoles 833--912)'


def test_systematics_in_quadrature(r33):
    """Objective 6 and box: '20% stat, 32--36% syst (~38--42% total)' from 10--20% (+) 30% (+) 5%. R4: the
    exact 41.5% prints 42 (the pre-ruling 41 came from combining the already-rounded 36)."""
    i = r33.intermediates
    near(i["syst_quadrature"][0], 0.3202, 0.001)
    near(i["syst_quadrature"][1], 0.3640, 0.001)
    near(i["total_error_quadrature"][0], 0.3775, 0.001)
    near(i["total_error_quadrature"][1], 0.4153, 0.001)
    assert [round(100 * x) for x in i["syst_quadrature"]] == [32, 36]
    assert [round(100 * x) for x in i["total_error_quadrature"]] == [38, 42]


def test_clock_and_utility_numbers(r33):
    i = r33.intermediates
    near(i["z_sample_ratio_fcc_over_lep"], 3.5e5, 0.01)     # 6e12 / 1.7e7
    near(i["sqrt_z_sample_ratio"], 594, 0.001)              # the chapter cites '~300' (arxiv_2407_09593)
    near(i["utility_musd_per_instance"], 2.0, 1e-12)        # 5% x $0.4B / 10


def test_draft_printed_hop_is_recorded(r33, r3d):
    """The draft's printed staggered C^S for Sigma(72x3) (FP resources.tex:40): 46243.6 + 690.2 log2(1/eps) T per link
    PER FIELD, kept as a record beside the r17 model (which adds the frame undo and prices MBU ladders and the full fit)."""
    i = r33.intermediates
    near(i["draft_hop_t_per_link_S72x3"], 46243.6 + 690.2 * L2_33, 1e-12)
    near(i["draft_hop_t_per_step_2033"], 3.7019e6, 0.001)
    near(i["draft_hop_t_per_step_2033_3fields"], 1.1106e7, 0.001)
    near(i["draft_printed_3fields_over_model_2033"], 1.5814, 0.002)   # r19: the model is below 3 x printed C^S
    near(r3d.intermediates["draft_hop_t_per_step_3d_3fields"], 3.3831e7, 0.001)


# --- 3+1D extension (plug-in prose app10:90; era 'codesign') -----------------------------

def test_3d_extension_register_and_depth(r3d):
    """'its 192 links at 1.7e5 T per link per step (5.5e4 gauge, 1.1e5 hop), plus 1.1e4 T of mass, give 3.2e7 T/step and
    3.2e10 T over 1e3 steps with no ramp priced --- 2.4x the 2033 register limit and 32x its depth limit' (r19 E21 hop,
    E20, eps = 3.6e-6)."""
    i = r3d.intermediates
    assert (i["V_3d"], i["n_links_3d"], i["n_plaq_3d"]) == (64, 192, 192)
    assert (i["lq_gauge_3d"], i["lq_fermion_3d"], i["lq_anc_3d"]) == (1728, 576, 100)
    assert i["lq_3d"] == 2404 and r3d.lq == (2404, 2404)
    near(i["lq_3d_over_marker"], 2.4, 0.002)
    assert i["magnetic_t_per_link_3d"] == magnetic_per_link("S72x3", "KS", 3, EPS3D)
    near(i["t_per_link_3d"], 5.4523e4, 0.001)               # '5.5e4 gauge'
    near(i["t_per_link_3d"] + i["t_hop_per_link_3d"], 1.6861e5, 0.001)   # '1.7e5 T per link per step'
    near(i["t_hop_per_link_3d"], 114087.31, 1e-6)           # '1.1e5 T hop': the staggered hop does not depend on d
    assert (i["hop_toffoli_per_link_3d"], i["hop_rotations_per_link_3d"]) == (3672, 2933)
    near(i["t_hop_per_step_3d"], 2.19048e7, 1e-5)
    near(i["t_mass_per_step_3d"], 10818.91, 1e-5)           # '1.1e4 T of mass'
    near(i["t_per_step_3d"], 3.2384e7, 0.001)               # '3.2e7'
    near(i["t_total_3d"], 3.2384e10, 0.001)                 # '3.2e10'
    assert r3d.hard_ops == (i["t_total_3d"], i["t_total_3d"])
    near(i["t_3d_over_cap"], 32.38, 0.002)                  # '32x'
    near(i["pre_r17_t_total_3d"], 9.5374e9, 0.001)          # the R-TOL print, record
    near(i["t_total_3d_round_e"], 7.6731e9, 0.001)          # round E's fixed 1e-3, record
    assert i["ramp_priced_3d"] is False
    near(r3d.breakdown_total(), i["t_total_3d"], 1e-12)
    near(i["eps_l_required_3d"], 3.0879e-12, 0.001)
    near(i["hop_over_retired_3d"], 10.738, 0.002)
    near(i["applied_stage_t_total_3d"], 1.25084e10, 0.001)  # the retired figure on this gauge step


def test_3d_cross_reference_to_ch10(r3d):
    """app07:197 prints ~2600--3100 LQ for the same lattice (1728 + 576 + 250--750 BE/bath ancilla, Ch. 10 ruling
    ch10-stretch-ancilla (a)) and no T-count. The chapter claims the match for the gauge and fermion register only."""
    from estimates import ch10_finite_density as ch10
    i = r3d.intermediates
    assert ch10.PUBLISHED["codesign"].lq == (2554, 3054) and ch10.PUBLISHED["codesign"].hard_ops is None
    assert i["ch10_lq_same_lattice_printed"] == (2600, 3100)
    assert i["ch10_lq_same_lattice_sum"] == ch10.PUBLISHED["codesign"].lq
    assert i["gauge_fermion_same_lattice_3d"] == 1728 + 576 == i["lq_3d"] - i["lq_anc_3d"]
    assert any("NO T-count" in n for n in r3d.notes)


# --- earlier prints are still reproduced, as a record ---------------------------------------

def test_legacy_numbers(a, r28, r33, r3d):
    """What the rewrite printed before round D (APPLY stage) and before the rulings (DERIVE stage)."""
    i28 = r28.intermediates
    near(i28["legacy_t_per_step_2028"], 2.480e4, 0.001)
    near(i28["legacy_t_total_2028"], 4.96e5, 0.001)          # APPLY: '20 steps of 2.5e4 T cost 5.0e5 T'
    assert i28["legacy_t_total_2028_derive"] == 1e5          # DERIVE: '<~1e5: ~5e3/step x ~20 steps + ramp'
    i = r33.intermediates
    assert i["legacy_lq_cut_2033"] == (900, 1100)
    assert i["legacy_t_per_plaquette_S72x3"] == float(GROUPS["S72x3"].plaquette_t) == 5e3
    assert i["legacy_t_per_step_2033"] == 5e5
    assert i["legacy_t_total_2033"] == (5.5e8, 1e9)
    assert r3d.intermediates["legacy_lq_3d"] == 2212
    assert r3d.intermediates["legacy_t_per_step_3d"] == 3e6 and r3d.intermediates["legacy_t_total_3d"] == 3e9
    b = dataclasses.replace(a, link_width_source=c.Stated("chapter", "app10:90", "test"))
    rb = m.model(b, "2033")
    assert rb.intermediates["lq_cut_2033"] == (900, 1100) and rb.lq == (900, 900)   # R11 (A): the box quotes V = 32
    assert m.model(b, "codesign").lq == (2212, 2212)
    assert m.model(b, "2028").lq == (230, 230)              # Z3: 2 either way


# --- the headline T responds to the lattice, the group, the fields and the tolerance ----------

def test_headline_t_scales(a):
    S = c.Stated
    r0 = m.model(a, "2028")
    t33 = m.model(a, "2033").hard_ops
    # volume: 32 -> 72 links at L = 6 (a larger circuit: more rotations, so a tighter eps under R-TOL)
    r = m.model(dataclasses.replace(a, L_2028=S(6, "app10:88")), "2028")
    i0, i6 = r0.intermediates, r.intermediates
    near(i6["eps_rot"], math.sqrt(EPS_SYN / (3 * (72 * 37.5 + 36 * 4))), 1e-12)
    near(r.hard_ops[0], 3 * (72 * (i6["t_per_link_2028"] + i6["t_hop_per_link_2028"]) + 36 * i6["hwp_t_k9"]), 1e-12)
    assert r.hard_ops[0] > 3 * (72 * (i0["t_per_link_2028"] + i0["t_hop_per_link_2028"]) + 36 * i0["hwp_t_k9"])
    assert r.lq == (2 * 36 * 2 + 9 * 36 + 22,) * 2
    # dimension: the magnetic multiplicities double at d = 3 (at the d = 3 circuit's own eps)
    r = m.model(dataclasses.replace(a, dim_2028=S(3, "app10:118")), "2028")
    e3 = r.intermediates["eps_rot"]
    near(r.intermediates["magnetic_t_per_link_2028"], 2 * magnetic_per_link("Z3", "KS", 2, e3), 1e-12)
    # tolerance: a tighter total synthesis error makes rotations cost more (R-TOL); the fixed eps_rot is no longer read
    r = m.model(dataclasses.replace(a, eps_syn=S(1e-8, "app10:88")), "2033")
    assert r.hard_ops[0] > t33[0]
    near(r.intermediates["t_per_rotation_report"][0], 1.15 * math.log2(1 / math.sqrt(1e-8 / N33[0])) + 9.2, 1e-12)
    assert m.model(dataclasses.replace(a, eps_rot=S(1e-8, "app10:88")), "2033").hard_ops == t33
    # priced volume: V = 40
    r = m.model(dataclasses.replace(a, L_par_priced=S(10, "app10:90")), "2033")
    near(r.hard_ops[0], m.model(a, "2033").intermediates["t_total_2033_at_V_hi"][0], 1e-12)
    near(r.hard_ops[1], m.model(a, "2033").intermediates["t_total_2033_at_V_hi"][1], 1e-12)
    # group: Sigma(36x3) has a compiled fast transform, an 8-qubit register and a cheaper U_mul
    r = m.model(dataclasses.replace(a, group_2033=S("S36x3", "app10:90")), "2033")
    assert r.hard_ops[0] < t33[0] and r.lq == (8 * 64 + 288 + 100,) * 2
    from estimates.groups import hop_link_cost
    e36 = r.intermediates["eps_rot"][0]
    near(r.intermediates["t_hop_per_link_2033"], hop_link_cost("S36x3", 3, eps=e36, legacy=False)["t"], 1e-12)
    # staggered fields move the T-count: the draft hop for 2 fields (g-only squish shared, rotations per field)
    r = m.model(dataclasses.replace(a, n_stag=S(2, "app10:88")), "2033")
    h2 = hop_link_cost("S72x3", 2, eps=r.intermediates["eps_rot"][0], legacy=False)
    near(r.intermediates["t_hop_per_link_2033"], h2["t"], 1e-12)
    assert h2["toffoli"] <= 3672 and h2["n_rot"] < 2933
    assert r.hard_ops[0] < t33[0] and r.lq[0] < 964


# --- PUBLISHED, INSTANCE_ROWS, guards ------------------------------------------------------

def test_published_is_the_boxes_as_printed():
    p = m.PUBLISHED
    assert set(p) == {"2028", "2033", "codesign"}
    assert (p["2028"].lq, p["2028"].hard_ops) == ((230, 230), (1.3e5, 1.3e5))
    assert (p["2033"].lq, p["2033"].hard_ops) == ((964, 964), (4.9e9, 2.6e10))   # step rule (r26 J1 readout ramp: 2.5e10)
    assert (p["codesign"].lq, p["codesign"].hard_ops) == ((2400, 2400), (3.2e10, 3.2e10))
    assert "NOT A BOX" in p["codesign"].src
    assert all(x.rel_tol <= 0.035 for x in p.values())
    assert not getattr(m, "DISPUTED", {})


def test_instance_rows_describe_the_rewrite(a):
    rows = {era: m.INSTANCE_ROWS(a, era, m.model(a, era)) for era in c.ERAS}
    assert [len(rows[e]) for e in c.ERAS] == [1, 3, 1]     # R7: the two dipole rows
    for era, rr in rows.items():
        for label, lq, t, extra in rr:
            assert isinstance(label, str) and len(lq) == 2 and len(t) == 2 and isinstance(extra, dict)
            assert not ({"ch", "era", "label", "lq", "t", "codesign", "src"} & set(extra))
            json.dumps([label, lq, t, extra])
    l28, l33, l3d = (rows[e][0] for e in c.ERAS)
    assert "readout pipeline" in l28[0] and r"\mathbb{Z}_3" in l28[0]
    assert l28[1] == (230, 230) and l28[3]["over_reference_t"] == 1.322 and "t_bound" not in l28[3]
    near(l28[2][0], 132174.88, 1e-6)
    assert "R-TOL" in l28[3]["note"] and "R-TOL" in l3d[3]["note"] and "3 steps (1 ramp + 2 evolution)" in l28[3]["note"]
    assert "accepted overshoot" in l28[3]["note"] and "E20" in l28[3]["note"] and "E20" in l33[3]["note"]
    assert "R-HOP" not in json.dumps(rows) and "working figure" not in json.dumps(rows)
    assert "unpublished" in l33[3]["note"] and "unpublished" in l3d[3]["note"]
    assert r"\Sigma(72{\times}3)" in l33[0] and l33[1] == (964, 964)
    near(l33[2][0], 4.9014e9, 0.001); near(l33[2][1], 2.6159e10, 0.001)   # step rule (was 2.5381e10)
    assert "readout map" in l33[3]["note"] and "quench-prepared" in rows["2033"][1][3]["note"]
    assert "once per link (E21)" in l33[3]["note"] and "1.1e5 T per link, 72% of the step" in l33[3]["note"]
    r33 = m.model(a, "2033")
    assert l33[3]["shots_campaign"] == [12672, 25344] == [round(x) for x in r33.shots]   # trend + L_par = 4 check
    assert l33[3]["shots_first_result"] == [122, 244]
    assert l33[3]["shots_fragmentation"] == [7.5e4, 1.2e6]
    assert l33[3]["wall_time_s_campaign"] == [round(x) for x in r33.wall_time_s]
    assert [rows["2033"][k][1] for k in (1, 2)] == [(995, 995), (912, 912)]
    assert l3d[1] == (2404, 2404)
    near(l3d[2][0], 3.2384e10, 0.001)
    blob = json.dumps(rows)
    for retired in ("2T", "167"):
        assert retired not in blob, retired


def test_resources_rows_build_in_memory():
    """The generator accepts the rows (nothing is written: the coordinator regenerates resources.json)."""
    from estimates.__main__ import resources_json
    inst = resources_json([6])["instances"]
    assert [(r["era"], r["lq"]) for r in inst] == [("2028", [230, 230]), ("2033", [964, 964]), ("2033", [995, 995]),
                                                   ("2033", [912, 912]), ("codesign", [2404, 2404])]
    near(inst[0]["t"][0], 132174.88, 1e-6)
    near(inst[1]["t"][1], 2.6159e10, 0.001)                 # step rule (was 2.5381e10)
    near(inst[4]["t"][0], 3.2384e10, 0.001)


def test_assumption_guards(a):
    S = c.Stated
    bad = [
        dict(link_width_source=S("dense", "app10:90")),
        dict(group_2033=S("SU3x", "app10:90")),
        dict(group_2033=S("SU3", "app10:90")),
        dict(group_2028=S("Z2", "app10:88")),
        dict(L_par_cut=S((10, 8), "app10:90")),
        dict(L_par_cut=S((8, 20), "app10:90")),
        dict(L_par_priced=S(16, "app10:90")),
        dict(dim_2028=S(1, "app10:118")),
        dict(dim_2033=S(3, "app10:140")),
        dict(n_stag=S(2.5, "app10:88")),
        dict(n_trot_2028=S(0, "app10:88")),
        dict(n_trot_2028=S(2.5, "app10:88")),
        dict(ramp_steps_2028=S(-1, "app10:88")),
        dict(ramp_steps_2028=S(None, "app10:88")),          # round-D ruling 1: the ramp has a length
        dict(ramp_steps_2033=S((-1, 10), "app10:90")),
        dict(n_pauli_hop_diagonal=S(0, "app10:88")),
        dict(moves_per_field=S(1.5, "app10:88")),
        dict(n_trot_2033=S(0, "app10:90")),
        dict(hadamard_amplitude=c.Assumed((0.2, 1.5), "x")),
        dict(n_s=c.Assumed(0.0, "x")),
        dict(delta_rel_2033=S(0.05, "app10:80")),
        dict(t_gate_s=S(1.0, "app10:128")),
        dict(window_over_a=S((50, 30), "app10:74")),
        dict(pair_clustering=S(2, "app10:80")),
        dict(delta_first_2033=S(1.5, "app10:135")),
        dict(workspace_links_2028=S(0, "app10:88")),
        dict(hamiltonian=S("Wilson", "app10:56")),
        dict(toffoli_convention=S("jones", "app10:88")),    # R5
        dict(synthesis_model=S("rus-slope", "app10:88")),       # E20: never slope-only for the gauge tables
        dict(fermion_synthesis=S("rus-slope", "app10:88")),     # r17: never slope-only for the hop and mass
        dict(n_pauli_hop_diagonal=S(10, "app10:88")),           # the shared Z3 hop structure is 8 strings
        dict(hop_share=S("site", "app10:95")),                  # E21: 'link' (headline) or 'draft' (sensitivity)
        dict(eps_syn=S(0.0, "app10:88")),                   # R-TOL: a positive total synthesis error
        dict(eps_syn=S(1.5, "app10:88")),
    ]
    for kw in bad:
        with pytest.raises((ValueError, TypeError)):
            dataclasses.replace(a, **kw)
    with pytest.raises(TypeError):
        dataclasses.replace(a, n_c=3)
    # legitimate: a zero-step ramp at either era, no mass term
    dataclasses.replace(a, ramp_steps_2033=S((0, 0), "app10:90"))
    dataclasses.replace(a, ramp_steps_2028=S(0, "app10:124"))
    dataclasses.replace(a, n_pauli_mass=S(0, "app10:88"))


def test_no_breakdown_row_is_negative(a):
    for era in c.ERAS:
        for p in m.model(a, era).breakdown:
            assert p.count >= 0 and p.t_each >= 0, (era, p.name)


# --- the chapter still says what the model reads -------------------------------------------

PHRASES = [
    # style pass 2026-10-08 (editorial_review/STYLE_GUIDE.md): 72 prose anchors whose sentences were cut or reworded are retired;
    # 11 box-row anchors follow the rewritten rows (rule 10, numbers unchanged). Numeric pins live in the model tests above.
    # r26 (referee report 2026-10-04: J1-J4, G7, G8)
    # r26 v2 (verifier): J1 error term, box-confined state, ramp lower bound; J2 wave-packet source and box modes;
    # J4 floor-curve shots and the q-hat convention; G8 preparation mismatch
    r"D^h_q(z,\mu^2) = \frac{z^{d-3}}{2N_c}\int",
    # r27 (2026-10-05): J1 V/PS and calibration, J2 box term, J3 N_h, G8 dipole penalty
    r"so nothing makes the box term cancel in the trend",
    r"and global analyses combine lattice constraints with collider data~\cite{Ji_LaMET,arxiv_1908_10439,arxiv_1711_07916}",
    r"the string wraps the axis $3.9$, $7.7$ and $11.6$ times",
    r"The smallest box with all three modes is $L_\parallel=16$ ($1828$ LQ); $L_\parallel=4$ ($532$ LQ) shares $3.1$ GeV, a one-point check of box dependence ($0.24$--$3.7$ yr).",
    r"${\gtrsim}\,4.7\times$ the shots if a bin's occupation equals its signal share",
    r"the pair meets again after $0.4$ fm/$c$, before the $0.7$--$1$ fm/$c$ formation time",
    r"\mathcal{N}_h(P_h)\,\bar\psi_q(0)\,W_0^\dagger",
    r"LaMET-type matching~\cite{Ji_LaMET,arxiv_1908_10439} has no established range of control for it",
    r"the slope of $P_s$ is $0.96$ (SU(3)) and $0.69$ (SU(2)) of that of $|W|^2$",
    r"$\Gamma=2|\mathrm{Im}\,V|N_c^2/(N_c^2-1)$",
    # R1
    # R2, R5: the convention sentence and the floor sentence
    # ruling R-TOL: the rule, one sentence, and each circuit's tolerance and T per rotation
    r"The total synthesis error is $10^{-2}$ per shot, and each circuit sets its tolerance from its own rotation count $N$: randomized synthesis errors add as $N\epsilon^2$, so $\epsilon=\sqrt{10^{-2}/N}$.",
    # round-D ruling 4: the Z3 basis; round-E ruling E1: the Fourier transform, displayed, with provenance
    # r17: the hop (Z3 derived here, Sigma(72x3) from the unpublished counts)
    # round-D ruling 1; the ruling of 2026-10-01: 1 ramp + 2 evolution steps, the overshoot stated plainly
    # round-E ruling E4: combined runs and the step-size cross-check (restored r18)
    # Aug-29 framing rulings (CHANGE_MAP ruling 6), restored after the round-E length cuts
    r"is measured first, in the opening shots of the 2033 fragmentation instance",
    r"2\!\cdot\!16\!\cdot\!2+3\!\cdot\!3\!\cdot\!16+22=230",
    # 2033
    r"The magnetic term is $9.4\times 10^{3}$ T per link per step~\cite{arxiv_2511_17437}",
    # r19 (E21 (1)): the color squish held per link; the undo drawn in the draft but not counted there
    # R3
    r"Required $\epsilon_l$ & $\lesssim 8\times 10^{-7}$ \\",   # style pass 2026-10-08: box row rewritten (rule 10); numbers unchanged
    r"Required $\epsilon_l$ & $\lesssim 3.8\times 10^{-12}$--$2.0\times 10^{-11}$ \\",   # style pass 2026-10-08: box row rewritten (rule 10); numbers unchanged
    # R4
    r"$\mathcal{A}\sim 0.05$--$0.2$",
    r"$32$--$36\%$ syst.\ ($\approx 38$--$42\%$ total)",
    r"systematics-dominated at $32$--$36\%$ ($\approx 38$--$42\%$ total)",
    r"20\%/bin stat., $32$--$36\%$ syst.\ ($\approx 38$--$42\%$ total)",   # cut pass 2026-10-06: the requirements accuracy row is cut; the box FCC row carries it
    # boxes
    r"Logical qubits & $230$ \\",   # style pass 2026-10-08: box row rewritten (rule 10); numbers unchanged
    # R16 sweep: the 2028 register split and step breakdown live in the plug-in prose (pinned above)
    r"Per-shot T & $1.3\times 10^{5}$ (1 ramp $+$ 2 evolution steps), $1.3\times$ the $10^{5}$ target \\",   # style pass 2026-10-08: box row rewritten (rule 10); numbers unchanged
    r"Shots & $\sim 10^{3}$ ($\sim 10\%$/channel)",
    r"\textbf{Wall time} & $\sim 2$ min ($\sim 0.13$ s/shot at $1\,\mu$s per T-gate; Ch.~\ref{ch:overview})",
    r"quench-injected source; one adiabatic ramp step, two evolution steps",
    r"Setup & 2+1D SU(3) via $\Sigma(72{\times}3)$ ($q_G{=}9$~\cite{arxiv_2511_17437}), $4{\times}8$, $a{\approx}0.1$ fm",
    r"Logical qubits & $964$ \\",   # style pass 2026-10-08: box row rewritten (rule 10); numbers unchanged
    r"Per-shot T & $4.9\times 10^{9}$--$2.6\times 10^{10}$ (window $+$ two ramps) \\",   # style pass 2026-10-08: box row rewritten (rule 10); numbers unchanged
    # R16 sweep: the 2033 step breakdown moved from the box to the plug-in prose
    # E27 (2026-10-02): one machine, shots in series; the horizon gap is a cost reduction, not a machine count
    r"First result & $120$--$240$ shots, $6.9$--$103$ d \\",   # style pass 2026-10-08: box row rewritten (rule 10); numbers unchanged
    r"Campaign & trend $+\,L_\parallel{=}4$ check: $1.3$--$2.5\times 10^{4}$ shots, $1.7$--$26$ yr (on one machine at $1\,\mu$s per T-gate, Ch.~\ref{ch:overview}; high end $5.2\times$ the five-year window) \\",   # style pass 2026-10-08: box row rewritten (rule 10); numbers unchanged
    # cut pass 2026-10-06: the gap paragraph's one-machine sentence restated the campaign box row (pinned above) and is cut
    # requirements and utility
    r"Number of subroutine calls & $\sim 130$",
    r"Maximum subroutine time & $4.9\times 10^{3}$--$3.7\times 10^{4}$ s per shot (2033); $\sim 0.13$ s (2028)",
    # cut pass 2026-10-06: the Hard-op row restated the box Per-shot T and is cut; the Size row is back (template rows kept for cross-chapter consistency)
    r"$3^3$: $995$ LQ, $5.4\times 10^{9}$--$1.5\times 10^{10}$ T, $653$ shots",
    r"$3{\times}3{\times}5$: $912$ LQ, $2.0$--$5.5\times 10^{8}$ T, $405$ shots",
    r"Per-shot T & $1.3\times 10^{5}$ (1 ramp $+$ 2 evolution steps",
    r"Per-shot T & $4.9\times 10^{9}$--$2.6\times 10^{10}$ (window $+$ two ramps",
    r"spread over $\sim 10$ application instances in a 5-year campaign, gives $\sim\$2$M per instance",
    # r25 (R1, R2, R4, R7, R10)
    r"Baryon number and strangeness are conserved, so baryons and kaons come in pairs",
    # q1 (2026-10-05, rulings A/B): quench-prepared media, E-rho-OQ dropped
    r"$\Sigma(72\!\times\!3)$ for SU(3) at $\sim 5\%$ truncation accuracy, $\Sigma(216\!\times\!3)$ for the static dipole and $2O$ for SU(2) dipoles",
    r"First: $\mathrm{Im}\,V$, one temperature & static dipole, $\Sigma(216{\times}3)$ $H_I$ $3^3$: $995$ LQ, $5.4\times 10^{9}$--$1.5\times 10^{10}$ T, $653$ shots, $57$--$98$ d \\",   # style pass 2026-10-08: box row rewritten (rule 10); numbers unchanged
    r"First: $\hat q$, one temperature & light-like, $2O$ $3{\times}3{\times}5$: $912$ LQ, $2.0$--$5.5\times 10^{8}$ T, $405$ shots, $1.2$--$2.3$ d \\",   # style pass 2026-10-08: box row rewritten (rule 10); numbers unchanged
    r"$7.5\times 10^{4}$--$1.2\times 10^{6}$ shots, ${\geq}\,12$--$1400$ yr \\",   # style pass 2026-10-08: box row rewritten (rule 10); numbers unchanged
    r"$\sim 130$: $3$ sampling configurations ($3$ source momenta, one spacing); $120$ fragmentation ($10{\cdot}4{\cdot}3$); $2$ dipole rows $\times$ $3$ temperatures \\",
]

RETIRED = [
    # q1 (2026-10-05, rulings A/B): E-rho-OQ dropped, dipole media quench-prepared
    r"E$\rho$OQ",
    r"exact thermal state",
    r"a one-ancilla Hadamard test of $W$",
    r"model-dependent small-$r$ byproduct",
    r"$2.3$ or $6.7\times 10^{9}$ T per shot",
    r"$7.2\times 10^{7}$ or $2.0\times 10^{8}$",
    r"three temperatures take $102$ and $1.9$ d",
    r"The per-shot T assumes the E$\rho$OQ configuration load",
    r"n_{\rm eff}^{2N_\ell}",
    r"arxiv_2001_11490",
    r"arxiv_2603_22422",
    # s648 (2026-10-05): the static dipole moves to Sigma(216x3) / H_I
    r"On pure-gauge $\Sigma(72\!\times\!3)$ $3^3$ this is $729+100+4=833$ LQ",
    r"static dipole, $\Sigma(72{\times}3)$ $3^3$: $833$ LQ",
    r"dipoles $833$--$912$",
    r"three temperatures take $13$ and $1.9$ d",
    # r27 (2026-10-05): superseded r26 wording
    r"symmetry-projected measurement not constructed here",
    r"a lower bound until $\mathcal{N}_h$ is priced",
    r"exponential in the number of links and not yet estimated",
    r"($\mathcal{N}_h$ unpriced)",
    r"whether box effects cancel in the trend is not established",
    # r26 (referee J1-J4, G7, G8): the r25 prints superseded by the readout ramp and the dipole conditionality
    r"$3.9\times 10^{9}$--$1.5\times 10^{10}$",
    r"$5.5$--$61$ d",
    r"$1.2$--$13$ yr",
    r"$9.3$--$820$",
    r"can raise both by up to $40\times$",
    r"Neither is a classically tractable validation",
    r"crosswalk",
    r"the draft draws it but does not count it",
    r"controlled only for $x,z\gtrsim 0.1$",
    r"FCC clock",
    r"FCC-clock",
    r"electric term a floor",
    r"box systematics largely cancel",
    # r26 v2 (verifier): superseded r26 v1 prints
    r"binomial counts need $615$ shots",
    r"$0.73$ d",
    r"measures $W$ directly",
    r"\mathcal{O}(\bar n^2V_{hh'}",
    r"diabatic errors that the 2028 emulators calibrate",
    r"on the hadronized out-state",
    r"the trend does.",
    r"mechanism-discriminating test",
    r"$\bar\psi_x U_{xy}\psi_y$",
    r"with large relative momentum",
    r"LaMET matching in the fragmentation instance",
    r"its readout cost is subleading",
    # the superseded chapter
    r"167",
    r"4.9\times 10^{5}",
    r"binary tetrahedral",
    # the rewrite before the rulings (q1 2026-10-05: the bare "5.5\times 10^{8}" and "20$ steps" retirements are dropped;
    # both strings are live in the quench text, 2.0--5.5e8 T and the 20-step coupling ramp)
    r"q_G=8",
    r"q_G{=}8",
    r"q_G\sim\log_2|G|",
    r"is minimal",
    r"900$--$1100",
    r"\approx 1700",
    r"1536+576+100",
    r"T/plaquette/step",
    r"from plaquettes",
    r"5\times 10^{5}$/step",
    r"\lesssim 10^{5}$: ",
    r"steps $+$ ramp \\",
    r"plus the ramp",
    r"38$--$41",
    r"$10^{3}$--$10^{4}$",
    r"$10^{5}$--$10^{6}$",
    r"$1\,\mu$s/T)",
    r"matching Ch.",
    # the APPLY stage, retired by round D
    r"working figure",                        # ruling 3: no unsourced hopping figure survives
    r"3.4\times 10^{5}",
    r"2.0\times 10^{6}$ working",
    r"5\times 10^{3}$ T/step",
    r"Twenty steps",
    r"5.0\times 10^{5}",
    r"no ramp length is set",
    r"We conservatively retain the current-insertion cost",
    r"No paper tabulates $\mathbb{Z}_3$ link primitives",
    r"single Hadamard test, 2028",
    r"one to a few string breaks",
    # superseded by the ruling 'Ch. 6 2028 steps' (2026-09-29)
    r"three evolution steps",
    r"3 evolution",
    r"the four steps",
    r"short of any string break",
    r"short of the first string break",
    r"2.0\times 10^{5}",
    r"t=0.3\,a",
    r"first toy $B/M$",
    r"2.3$--$4.3",
    r"9.2\times 10^{9}",
    # superseded by the round-E rulings (2026-09-29)
    r"0.995\times",
    r"at the report's working tolerance",
    r"the Fourier transform that paper's $14$ rotations",
    r"$176$ T magnetic",
    r"$443$ T electric",
    r"881",
    r"110.1",
    r"5.0\times 10^{4}",
    r"2.5$--$4.5",
    r"8.5\times 10^{9}",
    r"2.2$--$4.1",
    r"$25$--$45$",
    r"$50\times$ cut",
    r"caps coherent evolution at one step",
    r"\epsilon{=}10^{-4}",
    # the round-E apply's E4 overclaim, cut by the round-E fix (one extrapolation point only; no Trotter order is priced)
    r"second-order Trotter",
    r"extrapolated in $\Delta t$",
    r"the two evolution steps reach",
    r"the first quantity a pilot should measure",
    # superseded by ruling R-TOL (2026-09-29): the fixed 1e-3 and every number priced at it
    r"at $\epsilon=10^{-3}$ per rotation",
    r"\epsilon{=}10^{-3}",
    r"inside the $0.1$-fault budget",
    r"(not priced here)",
    r"1.003\times",
    r"3.34\times 10^{4}",
    r"94.8",
    r"759",
    r"$174$ T magnetic",
    r"$65$ T electric",
    r"($27$ T)",
    r"($160$ T)",
    r"($95$ T",
    r"7.6\times 10^{3}$ T per step",
    r"6.9\times 10^{3}",
    r"2.0\times 10^{6}$ T",
    r"2.2$--$3.9",
    r"7.7\times",
    r"$22$--$39$",
    r"2.6$--$4.6",
    r"$33\times$ cut",
    r"$1.0\times 10^{5}$ T",
    # R11 (2026-09-30): option (A) register at the priced volume; 2028 accuracy (a); the 3+1D per-link sum
    r"$964$--$1180$",
    r"$V{=}32$--$40$",
    r"up to $1.2\times$",
    r"$20$--$30\%$/channel",
    r"T gauge and $7.0\times 10^{3}$ T hop",
    r"prints for the same lattice",
    # R16 report-wide box sweep (2026-10-01): derivation chains cut from the boxes (finals kept, breakdowns in prose)
    r"$= 64$ gauge ($2{\cdot}16{\cdot}2$)",
    r"$= 3.24\times 10^{4}$/step",
    r"$(235$ gauge",
    r"(cut from $4{\times}16$)",
    r"the volume the T-count prices",
    r"$= 2.5$--$2.6\times 10^{6}$/step",
    r"$U_\times$ comparison$)$",
    r"($\bar n_s{\sim}0.1$, binomial)",
    r"($3\times 10^{3}\times\mathcal{A}^{-2}$)",
    # r17 (2026-10-01): R-TOL prints at the slope-only price and the R-HOP color move, retired
    r"9.7\times 10^{4}",
    r"3.24\times 10^{4}",
    r"$0.97\times$",
    r"$8\times 91.6",
    r"$92$ T, $11$ ancilla",
    r"$733$",
    r"by the rule above with two moves per field",
    r"is not a floor",
    r"compiled link-controlled color move",
    r"2.8$--$5.1",
    r"2.5$--$2.6",
    r"2.0$--$3.6\times 10^{-11}",
    r"$28$--$51$",
    r"$10.7$ T per rotation",
    r"$32\times$ cut",
    r"3.5$--$6.4\times 10^{9}",
    r"9.5\times 10^{9}",
    r"$9.5\times$ its depth",
    # r18 (2026-10-01): E20 retires the slope-only gauge tables; the 2028 shot goes back to 1 ramp + 2 evolution steps
    r"no offset",
    r"per gauge rotation",
    r"8.4\times 10^{4}",
    r"$0.84\times$",
    r"one evolution step, one evolution",
    r"1 evolution step;",
    r"one adiabatic ramp step, one evolution step",
    r"no time is reached twice",
    r"$173$ T magnetic",
    r"$60$ T electric",
    r"($25$ T)",
    r"($144$ T)",
    r"7.4\times 10^{3}$ T per step",
    r"$8\times 127",
    r"4.2\times 10^{4}",
    r"$42\times$ cut",
    r"2.4\times 10^{4}$ T per link",
    r"3.3$--$3.4\times 10^{4}",
    r"$85\%$ of it",
    r"1.5$--$2.8",
    r"$15$--$28\times$",
    r"1.9$--$3.6\times 10^{10}",
    r"4.5\times 10^{10}",
    r"$45\times$ its depth",
    r"$150$--$280$",
    r"3.5$--$6.5\times 10^{-12}",
    r"$\sim 0.08$ s",
    r"3.6 yr",
    r"$7.4$--$220$",
    # r19 (2026-10-01): E21 (1) holds the color squish once per link; the r18 2033 and 3+1D prints retired
    r"second uncompute",
    r"1.4\times 10^{4}$ Toffolis",
    r"1.9\times 10^{5}$ T per link",
    r"1.9\times 10^{5}$ hop",
    r"1.6\times 10^{5}$ T of it",
    r"$81\%$ of it",
    r"1.5\times 10^{7}$ T",
    r"1.5\times 10^{9}$",
    r"1.6$--$3.0",
    r"$16$--$30\times$",
    r"2.0$--$3.7\times 10^{10}",
    r"2.4\times 10^{5}$ T per link",
    r"4.7\times 10^{7}",
    r"4.7\times 10^{10}",
    r"$47\times$ its depth",
    r"$160$--$300$",
    r"3.4$--$6.2\times 10^{-12}",
    r"6 months",
    r"3.8 yr",
    r"1.2\times 10^{8}$ s",
    r"3.6\times 10^{10}$ s",
    r"$39$--$1100$",
    r"$7.7$--$230$",
    # r22 (2026-10-01): E26, a phasing group of k holds k - w(k) ancilla and repeat-until-success synthesis 1
    r"$11$ ancilla",
    r"+11\approx 220",
    r"$\sim 220$",
    r"one phasing group's workspace",
    # r23 (2026-10-02): E27, no machine count anywhere
    r"divides across machines",
    r"machines for the 5-year horizon",
    r"s serialized",
    r"shot-parallel",
    # r25 (2026-10-02): R1, R2, R4, R7, R10
    r"+8=216",
    r"We reserve this configuration",
    r"$1.1$--$2.0\times 10^{10}$ T per shot",
    r"$11$--$20\times$",
    r"$110$--$200$",
    r"4.9$--$9.1\times 10^{-12}",
    r"$1.4$--$2.5\times 10^{10}$",
    r"$8.8\times 10^{4}$",
    r"$26$--$770$",
    r"$5.2$--$150\times$",
    r"$7.8\times 10^{3}$",
    r"The six sampling configurations",
    r"at zero extra shot cost",
    r"binning the same final state in $z$",
    r"$t_{\rm had}\sim 10$ fm/$c$",
    r"\frac{1}{\bar n_s\,\delta_{\rm rel}^2}",
    r"Readout (a) &",
    r"$1.9\times 10^{5}$ T.",
]


@pytest.mark.skipif(not TEX.exists(), reason="chapter .tex not beside the scripts tree")
@pytest.mark.parametrize("phrase", PHRASES)
def test_tex_still_states_the_inputs(phrase):
    assert phrase in TEX.read_text(), f"chapter no longer says: {phrase}"


@pytest.mark.skipif(not TEX.exists(), reason="chapter .tex not beside the scripts tree")
@pytest.mark.parametrize("phrase", RETIRED)
def test_tex_carries_no_superseded_or_pre_ruling_number(phrase):
    assert phrase not in TEX.read_text(), f"chapter still carries: {phrase}"


@pytest.mark.skipif(not TEX.exists(), reason="chapter .tex not beside the scripts tree")
def test_tex_states_the_conventions():
    """R2 and R5 each ask for one sentence; R8 asks the boxes to point to Ch. 1; ruling 4 one sentence for the basis."""
    text = TEX.read_text()
    assert text.count("7 T per Toffoli") == 1               # the plug-in sentence (R16 sweep: box tags cut)
    assert text.count(r"$1.15\log_2(1/\epsilon)$ T per rotation") == 0   # E20: no slope-only price anywhere
    assert text.count(r"$1.15\log_2(1/\epsilon)+9.2$") == 1
    assert text.count("papers' convention") == 0             # R16 sweep: the convention lives in the plug-in sentence
    assert "stated lower bound for one" in text   # r26: 'floor' -> 'lower bound' (editorial); the box tag moved to the omission sentence (style pass 2026-10-08)
    assert text.count(r"$1\,\mu$s per T-gate") == 2
    assert text.count(r"The total synthesis error is $10^{-2}$ per shot") == 1   # R-TOL: one sentence states the rule
    assert text.count(r"\epsilon=\sqrt{10^{-2}/N}") == 1
    assert text.count(r"\cite{arxiv_2409_17349}") == 1      # round-E ruling E1: the published structure
    assert r"at $1\,\mu$s per T-gate; Ch.~\ref{ch:overview})" in text
    assert r"on one machine at $1\,\mu$s per T-gate, Ch.~\ref{ch:overview};" in text   # E27: serial, one machine
    assert text.count(r"\cite{FermionPrimitives_unpub}") == 2   # r17: 2028 hop scope, 2033 hop (r25: box row cut for length)
    assert "our estimate" in text                           # the frame undo is labeled as ours
    assert text.count(r"\cite{arxiv_2108_13305}") == 1      # the basis sentence


# --- open items (strict xfail) -----------------------------------------------------------------

def test_2028_benchmark_overshoot_is_the_accepted_one(a, r28):
    """Item 1 of apply_log/ch06_roundD.md. r18: the author accepts 1 ramp + 2 evolution steps over the 1e5 reference
    ('you can do the 1 ramp+2 evolutions if its only 1.27', 2026-10-01). At the E20 price the overshoot is 1.32x; the
    r18 brief adopts it up to ~1.35x. Fails if a later change pushes it past that, so the ruling is asked again."""
    over = r28.hard_ops[1] / float(a.rfi_t_2028)
    assert 1.0 < over <= 1.35
    assert r28.intermediates["t_total_2028_two_steps"] <= float(a.rfi_t_2028)   # the fitting alternative exists


def test_2033_synthesis_error_inside_the_fault_budget(a, r33):
    """Round-E item E2 CLOSED by R-TOL (apply_log/rtol_ch06.md): each 2033 end is priced at its own eps, N eps^2 = 1e-2."""
    assert r33.intermediates["synthesis_error_2033"][1] <= float(a.clean_shot_faults)
    assert max(r33.intermediates["synthesis_error_2033"]) <= float(a.eps_syn) * (1 + 1e-12)


def test_2028_shots_give_the_printed_accuracy(a, r28):
    """Item 3 CLOSED by R11 ruling (a): the box prints the accuracy 1e3 shots buy at n_s = 0.1, '~10%/channel'."""
    d = r28.intermediates["delta_rel_at_shots_2028"]
    near(d, 0.10, 1e-12)
    assert tuple(round(x) for x in r28.intermediates["shots_for_delta_rel_2028"]) == (1000, 1000)
    assert a.delta_rel_2028.lo * 0.9 <= d <= a.delta_rel_2028.hi * 1.1


def test_2033_hop_move_is_compiled(r33):
    """Item 4 of apply_log/ch06_roundD.md CLOSED at r17 (apply_log/r17_ch06.md): the U_mul color move (ruling VERTEX) is
    gone. The authors' unpublished draft compiles the color frame V_g, the eigen-class squish and the controlled
    diagonalizers, so the hop now rests on gate tables, not on a U_mul comparison."""
    rows = {p.name: p for p in r33.breakdown}
    assert not any("move" in n for n in rows)
    assert rows["hop_hop_squish_and_flags_S72x3_evolution_toffoli"].status is c.CircuitStatus.COMPILED
    assert rows["hop_hop_diagonalizers_S72x3_evolution_rotations"].status is c.CircuitStatus.COMPILED


@pytest.mark.xfail(strict=True, reason="r17 open item (apply_log/r17_ch06.md; NEEDS_AUTHOR 'For the Fermion_Primitives "
                                       "authors'). E21 (1) (r19) fixes the structure (undo counted, squish held once per "
                                       "link), but the undo V_g^dag and the SU(3) squish uncompute are still OUR ESTIMATES: "
                                       "the draft draws V_C^dag and does not count it. Closes when the draft counts them")
def test_2033_colour_frame_count_is_the_drafts(r33):
    rows = {p.name: p for p in r33.breakdown}
    assert rows["hop_colour_squish_and_parity_S72x3_evolution_toffoli"].status is c.CircuitStatus.COMPILED


def test_2033_t_covers_the_register_range(r33):
    """Item 6 CLOSED by R11 option (A): the box register is quoted at the volume the T-count prices (V = 32)."""
    i = r33.intermediates
    assert i["V_priced_2033"] == 32 and r33.lq == (i["lq_priced_2033"],) * 2 == (964, 964)
    assert i["lq_priced_2033"] == i["lq_cut_2033"][0]
    near(i["t_total_2033_at_V_hi"][0], 6.2e9, 0.02); near(i["t_total_2033_at_V_hi"][1], 3.2e10, 0.02)   # the L_par = 10 sensitivity


# --------------------------------------------------------------------------- #
# Referee G7 (2026-10-04): Trotter error at the chosen Delta t = 0.1 a
# --------------------------------------------------------------------------- #

def _staggered_step(dims, m_lat, dt):
    """Single-particle free staggered H and its Strang step over the mass and the 2d (direction, parity) hop groups."""
    import itertools
    np = pytest.importorskip("numpy")
    sites = list(itertools.product(*[range(L) for L in dims]))
    idx = {s: i for i, s in enumerate(sites)}
    V = len(sites)
    terms = [np.diag([m_lat * (-1) ** sum(s) for s in sites]).astype(complex)]
    for mu in range(len(dims)):
        for p in (0, 1):
            G = np.zeros((V, V), complex)
            for s in sites:
                if s[mu] % 2 == p:
                    t = list(s)
                    t[mu] = (t[mu] + 1) % dims[mu]
                    G[idx[s], idx[tuple(t)]] += -0.5j * (-1) ** sum(s[:mu])
                    G[idx[tuple(t)], idx[s]] += 0.5j * (-1) ** sum(s[:mu])
            terms.append(G)

    def ex(X, tau):
        e, v = np.linalg.eigh(X)
        return (v * np.exp(-1j * e * tau)) @ v.conj().T
    S = ex(terms[-1], dt)
    for X in reversed(terms[:-1]):
        h = ex(X, dt / 2)
        S = h @ S @ h
    return sum(terms), S, ex


def test_g7_2033_step_numbers(r33):
    """Step rule (2026-10-04, Claude's decision on G7). app10:73: the window evolves the ramp-prepared vacuum, an eigenstate
    whose error does not grow (free field 0.04-0.05 with the fermions), so the state-dependent estimate applies to the
    injected pair: a gluon phase of 0.06 (0.13) at 3.1 (4.65) GeV by 50 a; the 4.65 GeV point runs 575 steps there
    (0.0997). The literal connected-variance estimate on the vacuum (4.5-7.5, 2.0-4.3e3 steps) is a record: it counts the
    vacuum's zero-point fluctuation of the error operator as if it accumulated. Worst case 6.3e2-1.0e3 at g^2 = 1.
    The injected particles are a q-qbar pair; the 0.06 quark bound below is the lattice-corner mode, not a box momentum."""
    np = pytest.importorskip("numpy")
    i = r33.intermediates
    assert i["trotter_steps_window"] == (300, 500)
    w0, w1 = i["trotter_worst_err_2033"]
    assert round(w0, -1) == 630 and round(w1, -2) == 1000
    assert min(i["trotter_worst_err_min_g2_2033"]) > 100                  # no g^2 rescues the bound
    s0, s1 = i["trotter_sd_err_2033"]
    assert round(s0, 1) == 4.5 and round(s1, 1) == 7.5
    ph = i["trotter_pair_phase_by_box_mode"]
    assert round(ph[2][1], 2) == 0.06 and round(ph[3][1], 2) == 0.13 and ph[1][1] < 0.01
    # fermions: 3 fields x 3 colors of free staggered quarks on the 4 x 8 lattice, T = 0, m a = 0.02
    H, S, ex = _staggered_step((4, 8), 0.02, 0.1)
    e, v = np.linalg.eigh(H)
    nF = (v * (e < 0)) @ v.conj().T
    tot = []
    for n, g in zip((300, 500), i["trotter_free_err_2033"]):
        w = ex(H, 0.1 * n).conj().T @ np.linalg.matrix_power(S, n)
        inf_f = 9 * (1 - abs(np.linalg.det(np.eye(32) - nF + nF @ w)))
        assert inf_f < 0.15 * g * g / 2                                    # under 15% of the gauge infidelity
        tot.append(math.sqrt(g * g + 2 * inf_f))
    assert 0.04 <= min(tot) and max(tot) <= 0.05
    # quark dispersion phase over 500 steps: at most 0.06
    for k in range(32):
        if e[k] > 0:
            ov = v[:, k].conj() @ S @ v[:, k]
            assert abs(500 * (-np.angle(ov) - e[k] * 0.1)) < 0.065
    assert max(i["trotter_dipole_free_err"]) < 0.035                      # dipole rows at 0.05 a, a T = 0.45 (s648: static
    #                                                                       on the H_I dispersion, 0.029 / 0.034; was < 0.025)
    # the rule at the top box mode and in the thermal dipole rows
    assert i["trotter_rule"] == "state" and i["trotter_top_mode_steps_50a"] == 575
    assert round(i["trotter_top_mode_phase_50a"], 3) == 0.100 and i["trotter_top_mode_phase_50a"] <= 0.1
    assert i["trotter_free_err_top_50a"] < 0.03                           # free-field check at the 575-step shot
    assert i["trotter_dipole_rule_steps"] == (85, 240, 64, 181)            # s648: static on H_I (was 72, 203)
    assert max(i["trotter_dipole_rule_sd_err"]) <= 0.1 and max(i["trotter_dipole_rule_free_err"]) < 0.012
    assert min(i["trotter_dipole_sd_err"]) > 0.1                          # the 0.01 fm/c grid fails the rule there
    tex = TEX.read_text()
    assert "$6.3\\times 10^{2}$--$1.0\\times 10^{3}$" in tex and "the $4.65$ GeV point runs $575$ steps there" in tex
    assert "an eigenstate whose error does not grow ($0.04$--$0.05$ for free fields)" in tex


# --- r27 (chapter-open pass 2026-10-05): J1, J2, J3 derivations ---------------------------------------------
# (q1, 2026-10-05: test_r27_g8_dipole_eroq_penalty retired with the E-rho-OQ route, ruling A)


def test_r27_j2_box_term_and_check(r33):
    i = r33.intermediates
    near(i["source_box_modes_GeV"][0], 1.5498, 1e-4)
    near(i["yoyo_full_extension_fm"][2], 4.6494, 1e-4)
    assert [round(w, 1) for w in i["yoyo_windings_at_full_extension"]] == [3.9, 7.7, 11.6]   # printed, rounded once
    assert i["box_check_L_par_all_modes"] == 16 and i["lq_box_check_all_modes"] == 1828
    assert i["box_check_L_par_one_mode"] == 4 and i["box_check_shared_modes_one"] == (2,) and i["lq_box_check_one_mode"] == 532
    near(i["t_total_box_check_one_mode"][0], 2.413e9, 1e-3); near(i["t_total_box_check_one_mode"][1], 1.250e10, 1e-3)
    near(i["wall_box_check_one_mode_yr"][0], 0.242, 1e-2); near(i["wall_box_check_one_mode_yr"][1], 3.69, 1e-2)
    assert 0.16 < i["box_check_over_campaign"][0] < 0.17 and 0.16 < i["box_check_over_campaign"][1] < 0.17
    near(i["wall_trend_2033_yr"][1], 22.13, 1e-3)                       # the trend alone
    near(i["wall_campaign_2033_yr"][1], 22.13 + 3.69, 2e-3)             # the campaign includes the check (ruling 2026-10-05)


def test_r27_j3_nh_insertion_sum_rule(r33):
    i = r33.intermediates
    zc = i["ff_z_bin_centres"]
    assert len(zc) == 10 and abs(zc[0] - 0.14) < 1e-12 and abs(zc[-1] - 0.86) < 1e-12
    f_eq, f_var = i["ff_nh_shot_factor"]
    near(f_eq, sum(math.sqrt(z) for z in zc) ** 2 / 10, 1e-12); near(f_eq, 4.699, 1e-3)
    near(f_var, sum(z ** (2 / 3) for z in zc) ** 3 / 10, 1e-12); near(f_var, 23.07, 1e-3)
    # the Lagrange optimum: g_i ~ z_i^-1/2 saturates sum z_i g_i = 1 and gives the minimum
    g = [1 / (math.sqrt(z) * sum(math.sqrt(x) for x in zc)) for z in zc]
    assert abs(sum(z * gi for z, gi in zip(zc, g)) - 1) < 1e-12
    near(sum(1 / gi for gi in g) / 10, f_eq, 1e-12)
    assert i["ff_nh_extra_t_per_shot"] == 0.0


# --- s648 (author ruling 2026-10-05, (5)): static dipole on Sigma(216x3) / H_I ------------------------------------
def test_s648_static_dipole_pricing(a):
    """The static row's per-link step is the draft's H_I closed form at d = 3 with every rotation at the full fit
    (+9.2 T per rotation, 21 d + 1839 = 1902 rotations), U_phi twice; Sigma(72x3) / H_KS is unchanged by phi_mult."""
    from estimates.groups import c_t_stated
    for eps in (1e-5, 2.7634e-5, 1e-8):
        ax = m._at_eps(a, eps, 1)
        gs = m.gauge_step(ax, "S216x3", 3, "I")
        assert gs["rotations_per_link"] == 1902 and gs["multiplicity"]["U_phi"] == 2
        assert gs["electric_dense"] is None and gs["hamiltonian"] == "I"
        near(gs["per_link"], c_t_stated("S216x3", "I", 3, eps) + 9.2 * 1902, 1e-9)
        near(m.pure_gauge_step(ax, "S216x3", 3, "I")["per_link"], gs["per_link"], 1e-12)
        assert m.gauge_step(ax, "S72x3", 3)["multiplicity"]["U_phi"] == 1
    assert a.group_dipole_static.value == "S216x3" and a.hamiltonian_dipole_static.value == "I"
    assert a.group_2033.value == "S72x3" and a.hamiltonian.value == "KS"   # the species row is not switched


def test_s648_improved_free_dispersion():
    """H_I free limit (assumed tree-level Symanzik): on 3^3 every nonzero momentum component has khat^2 = 3, so
    w^2 rises by exactly 1.25; KS reproduces Ch. 5's frequencies and Lambda_sd."""
    from estimates import ch05_qgp_transport as ch5
    ks, imp = m.free_gauge_frequencies_h((3, 3, 3), "KS"), m.free_gauge_frequencies_h((3, 3, 3), "I")
    assert sorted(ks) == sorted(ch5.free_gauge_frequencies((3, 3, 3)))
    for w0, w1 in zip(ks, imp):
        near(w1 * w1, 1.25 * w0 * w0, 1e-12)
    T = 0.44 * 0.2 / 0.1973269804
    near(m.sd_lambda_h((3, 3, 3), 8, T, "KS"), ch5.trotter_state_dependent((3, 3, 3), 8, T, 1, 1, 0.1)["lambda"], 1e-12)
    near(m.sd_lambda_h((3, 3, 3), 8, T, "I"), 45.9006, 1e-5)                    # 32.92 at KS
