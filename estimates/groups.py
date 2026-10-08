"""Gauge-group registers and gate costs, as the published papers state them.

Rewritten 2026-09-28 from the LaTeX sources on the author's instruction to use
the papers' values rather than any inference. Every stated number carries its
paper and line; the quotable record is PAPER_COSTS.md next to this file.

Sources (bib keys): arxiv_2208_12309 (2T), arxiv_2312_10285 (2O),
arxiv_2405_05973 (Sigma(36x3)), arxiv_2408_00075 (fast Fourier transforms),
arxiv_2511_17437 (Sigma(72x3)), Gustafson_S648_inprep (Sigma(216x3), in preparation;
draft mains648.tex and cost_notes/, read 2026-10-05).

What the papers give, and what they do not:
- Register per link: 2T 5, 2O 6, Sigma(36x3) 8, Sigma(72x3) 9. The chapters'
  7 and 8 for the Sigma groups are ceil(log2|G|), an encoding no paper builds a
  circuit on. `link_qubits_chapter` records what the chapter prints;
  `link_qubits_compiled` what the circuits use; `link_qubits` is the ACTIVE
  value the models read: the compiled register (ruling R1, 2026-09-28).
- Primitive T-costs as "a + b log2(1/eps)" (Table tab:tgatecost of each paper),
  printed in the papers' convention: 7 T per Toffoli, 1.15 log2(1/eps) T per R_Z
  with no offset, coherent (linear) synthesis-error splitting. E20 ruling
  (2026-10-01, supersedes R2): every rotation in every circuit is priced at the
  full BRS fit 1.15 log2(1/eps) + 9.2 under R-TOL. Each primitive's rotation
  count is n_rot = b / 1.15 (the paper's printed coefficient); `a` stays.
  PrimitiveCost.t is that price (one function, `rotation_t`); t_papers and
  c_t_stated keep the papers' slope-only price as the legacy record.
- Multiplicities per link per Trotter step (tab:primcost, group-independent)
  and closed-form C_T per link per step. NO paper states a per-plaquette cost.
  The chapters' "per plaquette" figure is the magnetic term's per-link cost;
  `magnetic_per_link` and `electric_per_link` price the two pieces separately
  from the papers' tables so the split is explicit.
- Fast Fourier transforms exist for 2T, 2O and Sigma(36x3) (arxiv_2408_00075).
  None exists for Sigma(72x3); its paper gives a strict lower bound only. Z3's U_FFT
  (4 T + 2 rotations) is an explicit circuit DERIVED HERE (Ch. 6 round E, 2026-09-29) on
  the published two-level structure of arxiv_2409_17349; `t_gates` carries its 4 exact T.
- Not in any paper: a controlled/multiplexed group action; "190 T" (the figure
  Ch. 9 used until 2026-09-28, retired in favor of the sourced U_mul = 392 T);
  that the 2O fundamental irrep is the single-qubit Clifford group (Ch. 9 now
  states it as the authors' own observation).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Mapping

from .common import CircuitStatus, Cited, Stated, Tagged


# --------------------------------------------------------------------------- #
# Primitive costs, as stated
# --------------------------------------------------------------------------- #

@dataclass(frozen=True)
class PrimitiveCost:
    """T-count of one primitive. The paper prints t_const + t_log*log2(1/eps) (slope-only rotations).

    E20 ruling (2026-10-01; apply_log/r18_core.md): every rotation is priced at the full
    Bocharov-Roetteler-Svore fit 1.15 log2(1/eps) + 9.2 (common.t_per_rotation, "rus"). The paper's log
    coefficient is read as a rotation count, n_rot = t_log / 1.15 (PAPER_COSTS.md sec. 6-7), and the
    paper's constant t_const (Toffolis at 7 T, plus any exact T in t_gates) stays as printed:

        t(eps) = t_const + n_rot * (1.15 log2(1/eps) + 9.2)  =  t_papers(eps) + 9.2 n_rot

    `t_papers` is the paper's own slope-only price, kept as the legacy record (rule R2, superseded).
    `toffoli` is the count implied at 7 T per Toffoli (derived, flagged as such in `note`).
    """
    t_const: float
    t_log: float
    ancilla: int = 0
    status: CircuitStatus = CircuitStatus.COMPILED
    src: str = ""
    note: str = ""
    t_gates: float = 0.0    # exact T/T^dag gates inside t_const that are NOT from Toffolis (Z3 U_FFT: 4)

    @staticmethod
    def _snap(x: float) -> float:
        """The implied counts are integers by construction; snap float residue."""
        r = round(x)
        return float(r) if abs(x - r) < 1e-6 else x

    @property
    def toffoli(self) -> float:
        return self._snap((self.t_const - self.t_gates) / 7.0)

    @property
    def n_rot(self) -> float:
        """Synthesized rotations, from the paper's printed log coefficient: t_log / 1.15. Integer for every
        primitive except the dense 2O U_F (11370.1 -> 9887.04) and Sigma(36x3) U_F (185898 -> 161650.4),
        where the printed coefficient is carried exactly (carry exact, round once)."""
        from .common import RUS_SLOPE
        return self._snap(self.t_log / RUS_SLOPE)

    @property
    def rot(self) -> float:
        """Alias of n_rot (the name the chapter models use)."""
        return self.n_rot

    def t(self, eps: float) -> float:
        """T-count at per-rotation tolerance eps, report convention (E20): rotations at the full fit."""
        return self.t_const + self.n_rot * rotation_t(eps)

    def t_papers(self, eps: float) -> float:
        """LEGACY: the paper's own slope-only price t_const + t_log log2(1/eps) (rule R2, superseded by E20)."""
        return self.t_const + self.t_log * math.log2(1.0 / eps)

    def t_report(self, eps_rot: float, toffoli_convention: str = "textbook") -> float:
        """Re-costed at a named Toffoli convention; at "textbook" (7 T) this equals t(eps_rot)."""
        from .common import T_PER_TOFFOLI
        return self.t_gates + self.toffoli * T_PER_TOFFOLI[toffoli_convention] + self.n_rot * rotation_t(eps_rot)


def rotation_t(eps: float) -> float:
    """The ONE price of a synthesized rotation for every gauge primitive (E20): the full BRS fit
    1.15 log2(1/eps) + 9.2 (common.t_per_rotation, model "rus")."""
    from .common import t_per_rotation
    return t_per_rotation(eps, "rus")


def papers_convention_t(pc: "PrimitiveCost", eps: float) -> float:
    """LEGACY record: a primitive in its paper's slope-only convention (rule R2, superseded by E20)."""
    return pc.t_papers(eps)


def _pc(const, log, anc, src, note="", status=CircuitStatus.COMPILED, t_gates=0.0):
    return PrimitiveCost(const, log, anc, status, src, note, t_gates)


# Multiplicities per link per Trotter step (tab:primcost; identical in
# arxiv_2312_10285:797-810, arxiv_2405_05973:647-658, arxiv_2511_17437:934-947;
# the 2T paper :548-561 gives the H_I row). Callables of the spatial dimension d.
PRIMCOST = {
    "KS": {"U_F": lambda d: 2,       "U_Tr": lambda d: 0.5 * (d - 1),  "U_inv": lambda d: 3 * (d - 1),      "U_mul": lambda d: 6 * (d - 1)},
    "I":  {"U_F": lambda d: 4,       "U_Tr": lambda d: 1.5 * (d - 1),  "U_inv": lambda d: 2 + 11 * (d - 1), "U_mul": lambda d: 4 + 26 * (d - 1)},
}
MAGNETIC = ("U_inv", "U_mul", "U_Tr")   # the plaquette (Wilson-loop) term
ELECTRIC = ("U_F",)                      # plus U_phi where the paper lists it separately


@dataclass(frozen=True)
class Group:
    name: str
    order: int
    link_qubits_chapter: Tagged             # what the chapter prints
    link_qubits_compiled: Tagged | None     # the register every circuit in the paper is built on
    plaquette_t: Tagged | None              # the chapter's working "per plaquette" figure (not in any paper)
    primitives: Mapping[str, PrimitiveCost] = field(default_factory=dict)
    c_t: Mapping[str, object] = field(default_factory=dict)   # {"KS": f(d, eps), "I": ..., "I_fft": ...} per link per step, as stated
    benchmark: Mapping[str, float] = field(default_factory=dict)  # stated totals, L^3=10^3, N_t=50, eps_T=1e-8
    dominant: str = ""
    used_by: tuple[str, ...] = ()
    note: str = ""
    # U_phi (electric phase) calls per link per step, by Hamiltonian. Absent -> 1 (the Sigma(36x3)/(72x3) papers
    # add it once). The Sigma(216x3) draft's tab:primcost lists U_Ph at 1 (H_KS) and 2 (H_I).
    phi_mult: Mapping[str, float] = field(default_factory=dict)

    @property
    def link_qubits(self) -> int:
        """ACTIVE width the models read: the register the published circuits are
        built on (ruling R1, 2026-09-28). Falls back to the chapter value only
        where no compiled register is recorded."""
        src = self.link_qubits_compiled or self.link_qubits_chapter
        return int(src.lo)

    @property
    def link_width_conflict(self) -> bool:
        return (self.link_qubits_compiled is not None
                and self.link_qubits_compiled.lo != self.link_qubits_chapter.lo)

    @property
    def has_fft(self) -> bool:
        p = self.primitives.get("U_FFT")
        return p is not None and p.status is CircuitStatus.COMPILED


def _log(eps):
    return math.log2(1.0 / eps)


GROUPS: dict[str, Group] = {
    "Z2": Group("Z2", 2, Cited(1, "Lamm_QC_RoundTable"), Cited(1, "Lamm_QC_RoundTable"), None,
                used_by=("fig_qubit_budget",)),
    # Z3 primitives (Ch. 6 round D, H. Lamm ruling 2026-09-29: "the dihedral-group primitive costs with the extra Z2
    # dropped (D_N = Z_N x| Z_2) are an acceptable basis"; apply_log/ch06_roundD.md).
    # BASIS. arxiv_2108_13305 (Alam, Hadfield, Lamm, Li) gives the four D_N primitives. At m = 0 (the Z2 factor
    # dropped) they reduce to Z_N primitives: inversion k -> N-k (Eq. 17), multiplication k1+k2 mod N (Eq. 18),
    # trace by the direct Walsh method (Eqs. 21-22), Fourier transform the F_{Z_N} factor of Eq. 27 / Fig. 6, on
    # ceil(log2 N) qubits (Sec. III). WHAT IS DROPPED WITH THE Z2: the m qubit, every m-controlled gate of Figs.
    # 2-4 and 12 (the controlled two's complement of the multiplier), the Lambda_{m=0} controls of the trace, and
    # Phi(omega) plus the m-Hadamard of the Fourier transform. SCOPE LIMIT: that paper builds circuits only for
    # N = 2^n and prints no T-count, so its circuits give Z_{2^n}, not Z3 (checked: its inversion sends g=1 to the
    # forbidden |11>, its mod-4 adder fails on (1,2), (2,1), (2,2)). The Z3 circuits are therefore those of
    # arxiv_2408_00075 (main.tex): register :237-241 (|0>=|00>, |1>=|01>, |2>=|10>, |11> forbidden), chi = 2 CNOT
    # + X (Fig. clockmatrixqubit :259-270), controlled-sum rule (:233, :257), single-qubit control on |01>/|10>
    # (:306), X_{1,2} = SWAP (Fig. z2qubit :283-304), H_3 = 3 CNOT + 14 R_z (:250, :717). Each circuit is checked
    # on every basis state in tests/test_ch06_collider.py (test_z3_*). Convention as the other groups: 7 T per
    # Toffoli; the rotation count is the log coefficient / 1.15, each priced at the full fit (E20). COMPILED = the paper draws the gates; SCALING = a stated count or a
    # circuit derived here that no paper draws.
    "Z3": Group(
        "Z3", 3,
        link_qubits_chapter=Cited(2, "arxiv_2108_13305", "'|k> uses ceil(log2 N) qubits' (Sec. III) = ceil(log2 3)"),
        link_qubits_compiled=Cited(2, "arxiv_2408_00075", "qutrit on two qubits: |0>=|00>, |1>=|01>, |2>=|10>, |11> forbidden (:237-241)"),
        plaquette_t=None,
        primitives={
            "U_inv": _pc(0, 0, 0, "arxiv_2108_13305 Eq. (17) at m=0; arxiv_2408_00075 Fig. z2qubit :283-304",
                         "DERIVED HERE from the two sources: k -> -k mod 3 exchanges |01> and |10>, the drawn X_{1,2} = SWAP; "
                         "Clifford, 0 T. Dropped with the Z2: the m=0 control of Figs. 2, 3, 12 of arxiv_2108_13305"),
            "U_mul": _pc(28, 0, 0, "arxiv_2108_13305 Eq. (18) at m1=m2=0; arxiv_2408_00075 Fig. clockmatrixqubit :259-270, rule :257, :306",
                         "DERIVED HERE from the drawn chi and the paper's stated control rule: (g, h) -> (g, g+h mod 3) = chi on h "
                         "controlled on g0 (g=1), chi^-1 controlled on g1 (g=2); each controlled chi is 2 Toffoli + 1 CNOT: "
                         "4 Toffoli, no ancilla. Dropped with the Z2: the m2-controlled two's complement of Fig. 4 "
                         "(arxiv_2108_13305), whose mod-2^n adder does not apply to N=3"),
            "U_Tr":  _pc(0, 1.15, 0, "arxiv_2108_13305 Eqs. (21)-(22), direct Walsh method, Lambda_{m=0} dropped",
                         "DERIVED HERE: Re Tr(omega^g 1_3) = 3, -3/2, -3/2 on |00>, |01>, |10>; with the never-occupied |11> "
                         "(arxiv_2408_00075 :306) set to the |00> value, the diagonal is one Z Z rotation: CNOT, one R_Z, CNOT "
                         "(3 rotations if |11> is fixed at 0)",
                         status=CircuitStatus.SCALING),
            "U_phi": _pc(0, 1.15, 0, "arxiv_2108_13305 Eq. (22) Walsh method on the Fourier-basis electric diagonal",
                         "DERIVED HERE: a class-function electric term is equal on the conjugate irreps e=1, 2 "
                         "(-2cos(2 pi e/3) = -2, 1, 1), so diag(a, b, b, free) is one Z Z rotation",
                         status=CircuitStatus.SCALING),
            "U_F":   _pc(0, 16.1, 0, "arxiv_2408_00075 :250, :717",
                         "STATED: H_3 on the two-qubit register 'decomposed using available transpilers ... to 3 CNOTs and "
                         "14 R_z gates'; 14 x 1.15. The circuit is not printed. arxiv_2108_13305 Eq. (27) with m dropped "
                         "reduces U_F to F_{Z_N}, which it builds only for N=2^n (the standard QFT). KEPT FOR THE RECORD: "
                         "since Ch. 6 round E (H. Lamm, 2026-09-29) the models price U_FFT below",
                         status=CircuitStatus.SCALING),
            # Ch. 6 round E ruling (H. Lamm, 2026-09-29, 'Structured Z3 Fourier transform adopted'). Explicit circuit,
            # index 2 q1 + q0, g = 0, 1, 2 -> |00>, |01>, |10>:
            #   F3 (+) i = B X B^dag,  B = CNOT(q1->q0) CH(q0->q1),  CH = Ry_q1(pi/4) CZ Ry_q1(-pi/4),
            #   X = Ry_q0(th) S_q1 Z_q0 CZ Ry_q0(-th),  th = arccos(1/sqrt 3),  Ry(b) = S H Rz(b) H S^dag.
            # Rz(+-pi/4) = T / T^dag (4 T in all), 2 synthesized Rz(+-th), 2 CNOT + 3 CZ, 0 Toffoli. Verified in plain
            # Python (round-E scratchpad z3fft/circuit/circuit.py): 3x3 block = F3 to 5.6e-16, global phase 1, no
            # leakage to or from |11> (6.8e-17), |11> -> i|11>; checked again with numpy in the round-E apply. The
            # same two-level structure (Eq. (19): V = e^{-i pi/4 Y_b} e^{-i phi Y_a} e^{-i 3pi/4 X_b}, a = {0,1}, b = {1,2},
            # cos phi = 1/sqrt 3; checked numerically to 5e-16, ../notes/check_z3_fourier.py) is published by
            # arxiv_2409_17349 (Hidaka, Yamamoto, PRD 111 014510), Sec. IV Eqs. (18)-(20), V = F3 diag(1,1,i), in the
            # |00>,|01>,|11> encoding (one CNOT from this one); that paper gives no T-count. The qubit circuit and its
            # 4 T + 2 rotations are DERIVED HERE. Ancilla-free optimality evidence (Pauli-support bound, PTM filter,
            # one-rotation lemma) is in the scratchpad; it is a filter, not an exact-synthesis proof.
            "U_FFT": _pc(4, 2.3, 0, "arxiv_2409_17349 Sec. IV Eqs. (18)-(20) (structure); qubit circuit derived here",
                         "DERIVED HERE, circuit explicit and verified: F3 (+) i = B X B^dag with 2 CH (T, T^dag each: 4 T "
                         "exact) around X = Ry(th) [S_q1 Z_q0 CZ] Ry(-th), th = arccos(1/sqrt 3): 2 synthesized rotations. "
                         "4 + 2 x 1.15 log2(1/eps) = 34.56 T at eps=1e-4, 26.92 at 1e-3 (the transpiler U_F: 213.9, 160.4). "
                         "Structure as arxiv_2409_17349 Eqs. (18)-(20), which prints no T-count", t_gates=4),
        },
        dominant="U_mul, 70% of the gauge step at H_KS, d=2, eps=1e-3; the two U_FFT 23% (derived here). With the "
                 "transpiler U_F at eps=1e-4 it was U_F, 69%",
        used_by=("Ch9 2028 domain-wall benchmark", "Ch6 2028"),
        note=("Ch. 9 uses Z3 as the 2028 toy; Ch. 6 as the first-generation species-count benchmark; 2 qubits/link. "
              "Primitive table: the D_N primitives of arxiv_2108_13305 with the Z2 dropped, on the circuits of "
              "arxiv_2408_00075 (round D, 2026-09-29). Multiplicities are tab:primcost (group-independent); an "
              "abelian-specific plaquette would be cheaper."),
    ),

    "2T": Group(
        "2T", 24,
        link_qubits_chapter=Cited(5, "arxiv_2208_12309", "'five qubits ... per gauge link' (abstract :95)"),
        link_qubits_compiled=Cited(5, "arxiv_2208_12309", "binary |qponm>, :146; = ceil(log2 24), :136"),
        plaquette_t=Stated(1e3, "app03:78 / app10", "chapters' working figure; the paper's magnetic term per link at H_KS, d=2, eps=1e-4 is 1.09e3"),
        primitives={
            "U_inv": _pc(28, 0, 0, "arxiv_2208_12309 tab:tgatecost :526-544"),
            "U_mul": _pc(154, 0, 1, "arxiv_2208_12309 tab:tgatecost"),
            "U_Tr":  _pc(0, 12.65, 0, "arxiv_2208_12309 tab:tgatecost"),
            "U_F":   _pc(0, 1150, 0, "arxiv_2208_12309 tab:tgatecost", "naive; 1025 CNOT, 2139 R_Z, 1109 R_Y (:484)"),
            "U_FFT": _pc(98, 48.3, 2, "arxiv_2408_00075 tab:tgatecostbt :871-893"),
        },
        c_t={"I": lambda d, e: 4312 * d - 3640 + (4581.03 + 18.975 * d) * _log(e)},   # :564-566
        benchmark={"I": 2.0e10, "I_superseded_by_BO_paper": 1.1e11, "I_fft": 3.4e9},
        dominant="U_F, 44% (:568)",
        used_by=("Ch5 2028",),
        note="Toffoli 7 T; C^nNOT = 4(n-1) Toffoli + (n-1) clean ancilla; R_Z 'at worst 1.15 log(1/eps)' (:546).",
    ),

    "2O": Group(
        "2O", 48,
        link_qubits_chapter=Stated(6, "app06:99", "'six-qubit links'"),
        link_qubits_compiled=Cited(6, "arxiv_2312_10285", "'a total of six -- per gauge link' (abstract :104); 16 forbidden states (:153)"),
        plaquette_t=None,
        primitives={
            "U_inv": _pc(112, 0, 1, "arxiv_2312_10285 tab:tgatecost :774-793"),
            "U_mul": _pc(392, 0, 4, "arxiv_2312_10285 tab:tgatecost", "the '392 T' Ch. 9 cites; at 7 T/Toffoli"),
            "U_Tr":  _pc(350, 4.6, 2, "arxiv_2312_10285 tab:tgatecost", "188 + 4.6 log with 6 more ancilla (footnote :790)"),
            "U_F":   _pc(0, 11370.1, 0, "arxiv_2312_10285 tab:tgatecost"),
            "U_FFT": _pc(216, 48.3, 4, "arxiv_2408_00075 tab:tgatecostbt"),
            # Historical: Ch. 9 priced a '~190 T controlled 2O group action' until the
            # round-B ruling (2026-09-28) replaced it with U_mul. Never in the paper; no
            # controlled U_mul exists there. Nearest number is 188 (ancilla-assisted U_Tr).
            "controlled_group_action": _pc(190, 0, 0, "", "NOT IN PAPER; no controlled/multiplexed U_mul in arxiv_2312_10285. RETIRED 2026-09-28 (ruling VERTEX, H. Lamm): Ch. 9 prices each controlled 2O group action at one U_mul (392 T) instead; the 190-T figure is kept UNSOURCED for the record",
                                           status=CircuitStatus.UNSOURCED),
        },
        c_t={"KS": lambda d, e: 2863 * (d - 1) + (22737.8 + 2.3 * d) * _log(e),      # :813
             "I":  lambda d, e: 11949 * d - 10157 + (45473.3 + 6.9 * d) * _log(e)},   # :821
        benchmark={"I": 4.1e11, "KS": 2.0e11, "I_fft": 5.6e9},
        dominant="U_F, 99% (H_KS) / 98% (H_I) (:827)",
        used_by=("Ch9 2033",),
        note=("Toffoli 7 T; C^nNOT = 2ceil(log2 n)-1 Toffoli, n-2 dirty ancilla; R_Z 'on average 1.15 log2(1/eps)' (:795). "
              "The paper does NOT say the 2O fundamental irrep is the single-qubit Clifford group (app06:99 claims it)."),
    ),

    "S36x3": Group(
        "Σ(36×3)", 108,
        link_qubits_chapter=Stated(7, "app07:87", "= ceil(log2 108); no paper builds circuits on this"),
        link_qubits_compiled=Cited(8, "arxiv_2405_05973", "'an eight-qubit encoding' (abstract :128); '8 qubits ... or 3 qutrits and 2 qubits' (:186-188)"),
        plaquette_t=Stated((1e3, 1e4), "app07:96", "chapters' working figure; paper magnetic per link at H_KS: 2.4e3 (d=2), 4.9e3 (d=3) at eps=1e-4"),
        primitives={
            "U_inv": _pc(119, 0, 4, "arxiv_2405_05973 tab:tgatecost :622-642", "alternatives 203 T/2 anc, 420 T (App. A)"),
            "U_mul": _pc(308, 0, 2, "arxiv_2405_05973 tab:tgatecost"),
            "U_Tr":  _pc(378, 8.05, 7, "arxiv_2405_05973 tab:tgatecost"),
            "U_phi": _pc(0, 294.4, 0, "arxiv_2405_05973 tab:tgatecost", "electric phase; 256 R_Z + 254 CNOT (:588)"),
            "U_F":   _pc(0, 185898, 0, "arxiv_2405_05973 tab:tgatecost", "naive; 30956 CNOT, 2666 R_X, 32806 R_Y, 55234 R_Z (:586)"),
            "U_FFT": _pc(532, 117.3, 8, "arxiv_2408_00075 :835, tab:tgatecostbt", "T-width 8; 34 Toffoli + implied 102 R_Z"),
        },
        c_t={"KS":    lambda d, e: 2394 * (d - 1) + (371791 + 4.025 * d) * _log(e),     # :666-681
             "I":     lambda d, e: 9884 * d - 8414 + (744167 + 12.075 * d) * _log(e),
             "I_fft": lambda d, e: 9632 * d - 6034 + (1045.93 + 12.075 * d) * _log(e)},  # arxiv_2408_00075 tab:simcost
        benchmark={"I": 7.0e12, "KS": 3.5e12, "I_fft": 1.2e10},
        dominant="U_F, 'over 99% ... regardless of Hamiltonian' (:682); 51% with the FFT",
        used_by=("Ch10 2028", "Ch10 2033"),
        note="LINK-WIDTH: chapter 7 (ceil log2), compiled 8. Toffoli 7 T; C^nNOT = 2ceil(log2 n)-1; R_Z 1.15 log2 average; fiducial eps 1e-8 (:664).",
    ),

    "S72x3": Group(
        "Σ(72×3)", 216,
        link_qubits_chapter=Stated(8, "app03:68 / app10:78", "= ceil(log2 216); no paper builds circuits on this"),
        link_qubits_compiled=Cited(9, "arxiv_2511_17437", "'a register consisting of 9 qubits' (:150); '3 qutrits and 3 qubits' hybrid (:300-306)"),
        plaquette_t=Stated(5e3, "app03:78 / app10", "chapters' working figure; paper magnetic per link at H_KS, d=3, eps=1e-4 is 1.9e4 -- the 5e3 coincides with Sigma(36x3)'s"),
        primitives={
            "U_inv": _pc(637, 0, 4, "arxiv_2511_17437 tab:tgatecost :874-893"),
            "U_mul": _pc(1120, 0, 5, "arxiv_2511_17437 tab:tgatecost"),
            "U_Tr":  _pc(1414, 8.05, 12, "arxiv_2511_17437 tab:tgatecost"),
            "U_phi": _pc(0, 294.4, 0, "arxiv_2511_17437 tab:tgatecost"),
            "U_F":   _pc(112, 186295.4, 5, "arxiv_2511_17437 tab:tgatecost"),
            # NO FFT is constructed for Sigma(72x3). The paper's 'strict lower bound' (:905).
            "U_FFT": _pc(532, 515.2, 0, "arxiv_2511_17437 :905", "NOT CONSTRUCTED; 'estimated ... a strict lower bound'",
                         status=CircuitStatus.CONJECTURE),
        },
        c_t={"KS": lambda d, e: 9338 * d - 9114 + (372881 + 4.025 * d) * _log(e),     # :949-964
             "I":  lambda d, e: 38248 * d - 32046 + (745758 + 12.075 * d) * _log(e)},
        benchmark={"I": 7.1e12, "KS": 3.5e12},
        dominant="U_F, 'over 99%' (:965); '700 times more expensive than Sigma(36x3) with the fast Fourier transform'",
        used_by=("Ch5 2033", "Ch6 2033"),
        note="LINK-WIDTH: chapter 8 (ceil log2), compiled 9. Toffoli 7 T; C^nNOT = 4(n-1) clean; R_Z 1.15 log average; eps 1e-10 in Fig., eps_T 1e-8 benchmark.",
    ),

    # Sigma(216x3) (order 648; the qutrit Clifford group). Author ruling (H. Lamm, 2026-10-05): Ch. 5 2033 moves to
    # this group with H_I. Source: E. J. Gustafson, H. Lamm, E. M. Murairi, S. Osorio Perez, "Primitive quantum gates
    # for an SU(3) discrete subgroup: Sigma(216x3)", IN PREPARATION (bib Gustafson_S648_inprep); draft
    # /home/hlamm/Desktop/QC/S648_primitive_gates/mains648.tex (D below), cost notes in cost_notes/ (CN below).
    # QUBIT ENCODING only (D tab:qubitcosts :836-846); the draft's qubit-qutrit rows (tab:qutritcosts) are not used.
    # Same conventions as the other groups: 7 T per Toffoli, the printed log coefficient read as n_rot = b / 1.15,
    # each rotation at the full fit (E20). Multiplicities: D tab:primcost :870-883 = PRIMCOST, plus U_Ph at 1 (KS) /
    # 2 (I) per link per step (phi_mult). Composition reproduces D tab:costsummary :978-979 exactly (tests).
    # Rotations per link per step: KS 923 + 7 d, I 1839 + 21 d (= D's synthesis-error column, no 1/2 factor).
    "S216x3": Group(
        "Σ(216×3)", 648,
        link_qubits_chapter=Cited(11, "Gustafson_S648_inprep", "'can be encoded into a register of 11 qubits' (D :171); ceil(log2 648) = 10 is not used"),
        link_qubits_compiled=Cited(11, "Gustafson_S648_inprep", "omega^p C^q E^r V^(2s+t) X^u D^y: 4 qutrits (p,q,r,y) as qubit pairs + 3 qubits (s,t,u) = 11 (D :321-329)"),
        plaquette_t=None,
        primitives={
            "U_inv": _pc(2065, 0, 4, "Gustafson_S648_inprep tab:qubitcosts (D :840); CN cost_U_inverse.tex :145",
                         "T-width 2; 295 Toffoli; circuit Fig. dsbqubitinverse (D :542) with Fig. mathfrakvqubit"),
            "U_mul": _pc(2548, 0, 5, "Gustafson_S648_inprep tab:qubitcosts (D :841); CN cost_U_multiplication.tex :148",
                         "T-width 3; 364 Toffoli; circuit Fig. dsbmultiplyqubit (D :609)"),
            "U_Tr":  _pc(30786, 16.1, 10, "Gustafson_S648_inprep tab:qubitcosts (D :842); squish D :693, phasing D :695",
                         "T-width 1; squish 2 x 15393 T (4398 Toffoli) + 14 R_Z syntheses (15 rotations, 14 distinct "
                         "angles, x_{0,1,2} = x_{0,1,2,3} synthesized once and applied twice); 10 clean ancilla = 4 class "
                         "qubits + 6 scratch. Circuits Figs. s216x3tracesquishdiagram, dsbphasetrace, App. tracecircuits"),
            "U_phi": _pc(630, 11.5, 1, "Gustafson_S648_inprep tab:qubitcosts (D :844), sec:uph_impl (D :810); CN cost_U_phase.tex",
                         "kinetic phase U_Ph: decode-and-rotate on a 5-qubit irrep label, 90 Toffoli + 10 controlled R_Z "
                         "(ten nonzero phase classes); T-width 1; ancilla 1 in D (CN resource_comparison.tex :81 lists 5). "
                         "Called 1 x (H_KS), 2 x (H_I) per link per step"),
            "U_FFT": _pc(1086, 529, 1, "Gustafson_S648_inprep tab:qubitcosts (D :843); CN cost_U_fourier.tex :168-170",
                         "T-width 11. = D_D twiddle 551 + 9 rot, Phi^{72x3} 0 + 8 rot (fixed cube roots), qubit H_3 3 exact T + "
                         "10 rot, inner U_FFT^{Sigma(72x3)} 532 + 433 rot: 460 rotations. The inner FFT 'enters by citation' "
                         "to arxiv_2511_17437 (D :752), which prints only a lower-bound estimate 532 + 515.2 log (:905); the "
                         "draft's 497.95 log is the authors' own re-synthesis (CN resource_comparison.tex :746, :779-782)",
                         t_gates=3),
        },
        c_t={"KS": lambda d, e: 36876 * d - 34074 + (8.05 * d + 1061.45) * _log(e),      # D tab:costsummary :978
             "I":  lambda d, e: 135142 * d - 115216 + (24.15 * d + 2114.85) * _log(e)},  # D tab:costsummary :979
        benchmark={"KS": 2.0e10, "I": 6.1e10},   # D tab:fiducial :1021 (L^3 = 10^3, N_t = 50, eps_T = 1e-8, d = 3)
        dominant="no single primitive (unlike U_F for Sigma(72x3)). Report convention (E20), d=3, eps=1e-10: H_I U_mul 38%, "
                 "U_Tr 25%, U_FFT 24%, U_inv 13%, U_phi <1%; H_KS U_FFT 38%, U_Tr 26%, U_mul 25%, U_inv 10%",
        used_by=("Ch5 2033",),
        phi_mult={"KS": 1, "I": 2},
        note=("IN PREPARATION (Gustafson, Lamm, Murairi, Osorio Perez). Qubit encoding, 11 qubits/link, no ceil(log2) "
              "conflict. Freezing beta_f^{3+1} = 3.80(5) (Wilson, D tab:subgroups :190). Toffoli 7 T; R_Z 1.15 log2 per "
              "synthesis in the draft; eps_T 1e-8 fiducial; error budget (7d+923) / (21d+1839) d L^d N_t eps with no 1/2."),
    ),

    "SU2": Group("SU(2)", 0, Cited(12, "Lamm_QC_RoundTable", "digitized; fig_qubit_budget.py"), None, None,
                 used_by=("fig_qubit_budget",), note="Continuous group; order 0 = digitization-dependent."),
    "SU3": Group("SU(3)", 0, Cited(27, "Lamm_QC_RoundTable", "digitized; fig_qubit_budget.py"), None, None,
                 used_by=("fig_qubit_budget",), note="Continuous group; order 0 = digitization-dependent."),
}


# --------------------------------------------------------------------------- #
# Per-link-per-step pricing from the papers' own tables
# --------------------------------------------------------------------------- #

def _price(pc: PrimitiveCost, eps: float, papers: bool) -> float:
    return pc.t_papers(eps) if papers else pc.t(eps)


def magnetic_per_link(group: str, ham: str, d: int, eps: float, papers: bool = False) -> float:
    """T per link per Trotter step for the plaquette term.

    Sum over tab:primcost multiplicities of the magnetic primitives (U_inv, U_mul, U_Tr). This is what the
    chapters call their 'per plaquette' figure. Rotations at the full fit (E20); papers=True gives the
    papers' slope-only price (legacy record, rule R2).
    """
    g = GROUPS[group]
    return sum(PRIMCOST[ham][p](d) * _price(g.primitives[p], eps, papers) for p in MAGNETIC)


def electric_per_link(group: str, ham: str, d: int, eps: float, fft: bool = True, papers: bool = False) -> float:
    """T per link per Trotter step for the electric term: nF Fourier transforms
    (FFT if one is COMPILED for the group and fft=True, else the naive U_F)
    plus the electric phase U_phi where the paper lists it separately. Full fit (E20) unless papers=True."""
    g = GROUPS[group]
    # A group with no naive U_F in its paper (Sigma(216x3): the draft builds only the FFT) always uses U_FFT.
    key = "U_FFT" if ((fft and g.has_fft) or "U_F" not in g.primitives) else "U_F"
    t = PRIMCOST[ham]["U_F"](d) * _price(g.primitives[key], eps, papers)
    if "U_phi" in g.primitives:
        t += g.phi_mult.get(ham, 1) * _price(g.primitives["U_phi"], eps, papers)
    return t


def c_t_stated(group: str, ham: str, d: int, eps: float) -> float | None:
    """LEGACY record: the paper's closed-form C_T per link per step, as printed (slope-only rotations)."""
    f = GROUPS[group].c_t.get(ham)
    return None if f is None else f(d, eps)


def c_t_log_coefficient(group: str, ham: str, d: int) -> float | None:
    """The printed log2(1/eps) coefficient of the paper's C_T at dimension d (C_T(eps=1/2) - C_T(eps=1))."""
    f = GROUPS[group].c_t.get(ham)
    return None if f is None else f(d, 0.5) - f(d, 1.0)


def c_t_rotations(group: str, ham: str, d: int) -> float | None:
    """Rotations per link per step implied by the paper's C_T: log coefficient / 1.15 (E20)."""
    from .common import RUS_SLOPE
    b = c_t_log_coefficient(group, ham, d)
    return None if b is None else b / RUS_SLOPE


def c_t_full(group: str, ham: str, d: int, eps: float) -> float | None:
    """The paper's closed-form C_T per link per step at the report convention (E20): the printed formula plus
    9.2 T per implied rotation, i.e. const + n_rot (1.15 log2(1/eps) + 9.2), n_rot = log coefficient / 1.15."""
    f = GROUPS[group].c_t.get(ham)
    if f is None:
        return None
    from .common import RUS_OFFSET
    return f(d, eps) + RUS_OFFSET * c_t_rotations(group, ham, d)


def plaquette_t_from_primitives(group: str, toffoli_convention: str, t_per_rot: float,
                                multiplicity: Mapping[str, float]) -> float:
    """Sum multiplicity[prim] x (toffoli x T + rot x t_per_rot). Kept for callers
    that price in the report's convention with their own multiplicities."""
    from .common import T_PER_TOFFOLI
    g = GROUPS[group]
    tt = T_PER_TOFFOLI[toffoli_convention]
    return sum(m * (g.primitives[p].t_gates + g.primitives[p].toffoli * tt + g.primitives[p].rot * t_per_rot)
               for p, m in multiplicity.items())


# --------------------------------------------------------------------------- #
# Rule R-HOP: a gauge-covariant fermion hop, per link per Trotter step
# --------------------------------------------------------------------------- #
# DERIVED HERE (Ch. 6 round D, apply_log/ch06_roundD.md; ruled report-wide by H. Lamm 2026-09-29,
# TRACKED_CHANGES.md "RULING, R-HOP report-wide"). No published circuit prices a gauge-covariant hop
# on a discrete-group link; the rule is assembled from costs the report already uses:
#
#   T_hop = moves_per_field N_stag C_W(G) + n_P(G) HWP_T(k),     k = N_c N_stag by default
#
#   C_W(G)  = 0 where D(g) is diagonal on the register (Z3 as the center of SU(3): D(g) = omega^g on
#             every color), otherwise one U_mul of the group's published table (ruling VERTEX), status
#             SCALING: a comparison, no compiled circuit, not a floor.
#   moves   = 2 per field (move across the link and back: W and W^dag).
#   n_P(G)  = Pauli strings per color-flavor copy: 2 for the plain XX+YY bilinear after a move; 8 for
#             the Z3-dressed hop in the group basis of arxiv_2408_00075 (Walsh expansion of Re/Im omega^g
#             after one Clifford CNOT; checked in tests/test_ch06_collider.py test_z3_hop_*).
#   HWP_T   = common.hwp_group: one Hamming-weight-phasing group of k per (link, Pauli string),
#             floor(log2 k)+1 rotations at the chapter's synthesis tolerance + k - w(k) Toffolis,
#             measurement-based uncompute (Ch. 9's rule).
# Default synthesis since r17 (2026-10-01) is the report's full fit 1.15 log2(1/eps) + 9.2 ("rus"); it was the
# slope-only "rus-slope" (15.3 T per rotation at 1e-4), a latent bug no ruling adopted. Every chapter passes
# `synthesis` explicitly, so the fix moves no chapter number. 7 T per Toffoli. The staggered mass is NOT part of
# R-HOP (see staggered_mass_site).
#
# LEGACY (r17, author ruling (e), 2026-10-01): R-HOP is RETIRED for the non-diagonal groups (2T, 2O, Sigma(36x3),
# Sigma(72x3)) in favor of hop_link_cost below, which prices the hop from the authors' unpublished Fermion_Primitives
# gate counts. rhop_link is kept, unchanged in logic, as the recorded legacy comparison; hop_link_cost reports it
# beside the new price. For Z3 (diagonal link, no draft entry) hop_link_cost still uses the R-HOP structure, priced
# at the full fit.

DIAGONAL_LINK = frozenset({"Z3"})        # D(g) diagonal on the published register: no color move
RHOP_N_PAULI_DIAGONAL = {"Z3": 8}        # dressed hop, strings per copy (derived here, Ch. 6 round D)
RHOP_N_PAULI_MOVED = 2                   # plain XX+YY after the move (app06:108, app10:88,90)
RHOP_MOVES_PER_FIELD = 2                 # W and W^dag per field per link (app10:88,90)
RHOP_SRC = ("DERIVED HERE: rule R-HOP (apply_log/ch06_roundD.md; TRACKED_CHANGES.md E6b and 'RULING, R-HOP "
            "report-wide'); moves at U_mul by ruling VERTEX; phasing arxiv_1709_06648,arxiv_1902_10673")


def rhop_link(group: str, n_stag: int, n_c: int, eps: float, *, k: int | None = None,
              n_pauli: int | None = None, moves_per_field: int = RHOP_MOVES_PER_FIELD,
              synthesis: str = "rus", toffoli_convention: str = "textbook") -> dict:
    """LEGACY rule R-HOP: T per link per Trotter step of the hop of `n_stag` fields of `n_c` colors.

    Superseded for non-diagonal groups by hop_link_cost (r17); kept as the recorded comparison.

    `k` is the number of copies that share one link and one Pauli string and so form one phasing
    group; the default N_c N_stag is Ch. 6's. A chapter whose own grouping ruling differs (Ch. 9's
    color-and-s copies) passes its k. `n_pauli` and `moves_per_field` default to the rule; a chapter
    passes its own Stated inputs so their provenance stays in its Assumptions.

    Returns {"k", "hwp", "n_pauli", "n_moves", "t_move", "t_moves", "t_phasing", "t"}: `hwp` is the
    common.hwp_group dict, t_moves = n_moves t_move, t_phasing = n_pauli hwp["t"], t = the sum.
    """
    from .common import hwp_group
    if group not in GROUPS:
        raise KeyError(f"unknown group {group!r}")
    if n_stag < 1 or n_c < 1:
        raise ValueError("R-HOP needs n_stag >= 1 and n_c >= 1")
    if k is None:
        k = n_c * n_stag
    h = hwp_group(k, eps, synthesis, toffoli_convention)
    if group in DIAGONAL_LINK:
        n_moves, t_move = 0, 0.0
        if n_pauli is None:
            n_pauli = RHOP_N_PAULI_DIAGONAL[group]
    else:
        if "U_mul" not in GROUPS[group].primitives:
            raise ValueError(f"R-HOP: {group} has no published U_mul to price the color move")
        n_moves = moves_per_field * n_stag
        t_move = GROUPS[group].primitives["U_mul"].t(eps)
        if n_pauli is None:
            n_pauli = RHOP_N_PAULI_MOVED
    t_move_total = n_moves * t_move
    t_phase = n_pauli * h["t"]
    return {"k": k, "hwp": h, "n_pauli": n_pauli, "n_moves": n_moves, "t_move": t_move, "t_moves": t_move_total,
            "t_phasing": t_phase, "t": t_move_total + t_phase}


# --------------------------------------------------------------------------- #
# Fermion hop and mass from the authors' unpublished Fermion_Primitives draft (r17)
# --------------------------------------------------------------------------- #
# Author ruling (e), H. Lamm 2026-10-01: "adopt the fermion primitive numbers ... estimate or approximate the
# missing pieces". Source: the unpublished draft "Primitive Quantum Gates for Lattice Gauge Theories with Fermions"
# (Hidalgo, Murairi, Lamm, Gustafson), work/extracted/Fermion_Primitives/publication_draft/ (FP below); cite as
# unpublished gate counts (private communication), bib key FermionPrimitives_unpub. Every count below is read from
# FP's gate tables; the full line-by-line record is apply_log/r17_core.md.
#
# STRUCTURE of one gauge-covariant hop on a link (FP section_hamiltonian.tex:42-62, mapping_circuits/
# FinalDiagonalizations.tex, HoppingDiag.tex):
#   color frame V_g on both sites (FP's "2 C^G", su2_diag.tex:238)  ->  [spinor frame V_S, Wilson only]
#   ->  eigen-class squish of g  ->  per (class, color[, spinor]) controlled diagonalizer T  ->  phasing Rz(-th), Rz(th)
#   ->  T^dag  ->  squish^dag  ->  V_g^dag on both sites (the undo; drawn in FP, never counted, see FRAME UNDO below).
#
# FRAME UNDO (r19, author ruling E21 (1), H. Lamm 2026-10-01 "make the switch"; frame-undo reading, two readers +
# judge). FP draws V_C^dag on both sites but prices only the forward frame:
#   - FP's "factor of 2" / "2 C^G" is two SITES, forward only: su2_diag.tex:238 ("The squish operations have to be
#     applied to neighboring sites in the hopping term, which yields an extra factor of 2"); resources.tex:15
#     (C^Stag = 2 C^G + C^hop).
#   - SU(3) C^G is one forward application on one site, one spinor register: su3_diag.tex:130-134 = rotations
#     (:121-126) + one squish (:86-91) + flags (:208-213); Sigma(72x3) 1692.8 + 12789 + 892 x 7 = 20725.8 exactly.
#   - The printed staggered totals are 2 C^G + 2 C^hop(table) exactly (resources.tex:37-41 with Tab fermion_costs
#     :59-67): Sigma(72x3) 2(20725.8) + 2(2396) = 46243.6, log 2(211.6) + 2(133.5) = 690.2; Sigma(36x3) 26041.6,
#     644.2. The 2 on C^hop is the hop T_diag / EigSquish compute + uncompute (FinalDiagonalizations.pdf,
#     hopping.tex:323-333), which IS priced; the 2 on C^G is the two sites.
#   - V_C^dag(x), V_C^dag(y) are drawn in big_fig_group_agnostic.pdf (fig da_diagonalizer, section_hamiltonian.tex:
#     86-91), in the identical Fermionic_Primitive.pdf (running_notes.tex:124-128) and written as V_g^dag(s+e)
#     V_g^dag(s) in running_notes.tex:213-225, but no count or sentence anywhere prices them.
#   So the model adds V_g^dag on both sites: 4 color-rotation layers per link per field per spinor component
#   (undo=True, the default). Adjacent links carry different g, so the frame cannot telescope into the next link.
#
# SQUISH / PARITY MULTIPLICITY (a separate assumption from the undo). The color squish and parity flags depend on
# g only, and g is untouched between V(x), V(y), hop, V^dag(x), V^dag(y). Headline since r19 (E21 (1)):
# share="link", computed once per link, held live through all four frame applications and the hop, then
# uncomputed (compute + uncompute = 2 squish + parity per link; su2_diag.tex:238 allows the memory/operations
# trade-off). Sensitivity share="draft" (the r17 headline, conservative): FP's per-site compute + uncompute extended
# to every application, 4 x (2 squish + parity) per link, and the hop squish + flags per field. FP never calls for
# re-squishing for the undo, nor for an SU(3) squish uncompute; the draft variant is the model's extension.
#
# COMPONENT COUNTS, as FP prints them, and the reading used here:
#   color squish   C^nX table (su2_diag.tex:186-191; su3_diag.tex:98-106). FP's printed T (917, 189, 7854, 12789,
#                   su2_diag.tex:207-211, su3_diag.tex:86-92) equal the table at 2(n-2) Toffolis per C^nX (n>=3), 1 per
#                   C^2X. Under the report's measurement-based uncompute (MBU) a C^nX costs n-1 Toffolis.
#   color parity   SU(2): 44 (BT), 56 (BO) Toffolis (su2_diag.tex:225). SU(3): flag Toffolis 778 (S36), 892 (S72)
#                   (su3_diag.tex:202-215); FP prices S36's at 2 T each in C^G = 10918.8 (su3_diag.tex:132), here 7.
#                   Kept as printed (FP gives totals only, no C^nX split for MBU).
#   color rotations per frame application per field: BT 52, BO 60 (su2_diag.tex:225-229, 4 R_Z per controlled
#                   rotation); S36 164, S72 184 (su3_diag.tex:121-126, 194-198: 2 N_angles,
#                   author-confirmed E21 (3); the "6 N_angles" in the text at su3_diag.tex:120 is a draft slip).
#   squish copies   SU(2) C_V counts compute + uncompute (su2_diag.tex:232-238); SU(3) C^G counts the squish once
#                   (su3_diag.tex:128-134). ESTIMATED HERE: the SU(3) uncompute is a second squish, as in SU(2).
#   hop squish      C^nX table (hopping.tex:202-213, 253-266). SU(2) printed T (84, 434) equal the table at 2n-3. SU(3)
#                   printed T (136, 372) match no convention; the table at 2(n-2) gives 469 and 616, the values FP
#                   comments out (hopping.tex:262-263). Here: the table under MBU (n-1).
#   hop flags       SU(3) 100 Toffolis (hopping.tex:268-281); SU(2) none. Kept as printed.
#   diagonalizer    C(T) = 38.8 + 4.45 log2(1/eps) (hopping.tex:309-311) = 2 T (controlled H) + 4 rotations at
#                   9.2 + 1.15 L; the printed slope 4.45 is read as 4 x 1.15 = 4.6 (38.8 = 2 + 4 x 9.2 fixes the
#                   count at 4). One T and one T^dag per (eigen-class, color) (FinalDiagonalizations.tex). FP's text
#                   C_T(hop) (hopping.tex:335-347) reproduces exactly as 2 squish (+ 2 flags) + n_diag C(T), with
#                   n_diag = 2 n_cls N for BO and the SU(3) groups (28; 60) but n_cls N for BT (10): BT counts T only.
#                   Here: 2 n_cls N for every group.
#   phasing         Rz(-th) Rz(th) per (color[, spinor]) pair (FinalDiagonalizations.tex); not in any FP constant.
#                   After the frames every pair carries the same |phi lambda| = 1 (staggered) or the same 1-r gamma^0
#                   weight 2 (Wilson, r=1), so the 2 N n_spin n_stag rotations on one link are one equal-angle
#                   Hamming-weight-phasing group (HWP; arxiv_1709_06648, arxiv_1902_10673). A block encoding drops them.
#   spinor (Wilson) per color 280 floor((d-1)/2) + 16 T (resources.tex:60-67) = 4 W gates x (4 T + 10 Toffolis) at
#                   d=3,4 (spin_diag.tex:282-297); 16 T and no Toffoli at d=1,2. The 10 Toffolis per W are 5 triply-
#                   controlled rotations at 2 compute + 2 uncompute (spin_diag.tex:280); MBU halves them to 5.
#
# PER LINK, per Trotter step, n_stag fields sharing the link (share="link", undo=True: the default since r19, E21):
#   frame:   n_app = 4 applications (2 sites x apply/undo; 2 without the undo). g-only part: 2 x squish + parity
#            Toffolis ONCE per link (share="link"); n_app x that under share="draft". Plus n_app x n_stag x n_spin
#            color-rotation sets: the rotations act on each field's color modes and are not shared
#            (correction to "one shared C^G": FP's formula is per field).
#   hop:     g-only hop squish + flags 2 (squish + flags) Toffolis once per link (n_stag x under share="draft"),
#            plus 2 n_cls N n_spin n_stag diagonalizers.
#   phasing: one HWP group of k = 2 N n_spin n_stag (or plain, or none for a block encoding).
#   spinor:  n_stag x N x C^spinor(d), Wilson only.
# AUTHOR-CONFIRMED (E21 (2), (3), 2026-10-01): n_spin = 2 for Wilson d = 3 (one color frame per hopping spinor
#   component; 1 - gamma^0 keeps 2 of 4 at r = 1); SU(3) color rotations = 2 N_angles per frame application
#   (184 = 2 x 92 for Sigma(72x3), su3_diag.tex:194-198; the "6 N_angles" at su3_diag.tex:120 is a draft slip).
# Z3 has no FP entry; D(g) is diagonal on the register: the R-HOP structure (C_W = 0, 8 strings, HWP(k)) at the full fit.

MCX_CONVENTIONS = ("mbu", "2n-3", "draft-color")


def mcx_toffolis(n_controls: int, convention: str = "mbu") -> int:
    """Toffolis in one C^nX. "mbu": n-1 (the report's rule: the n-2 AND temporaries uncomputed by measurement);
    "2n-3": clean-ancilla ladder with coherent uncompute (FP's hop squish); "draft-color": 2(n-2) for n >= 3,
    the convention that reproduces FP's printed color-squish T (917, 189, 7854, 12789)."""
    n = int(n_controls)
    if n < 1:
        raise ValueError("need >= 1 control")
    if n == 1:
        return 0
    if n == 2:
        return 1
    if convention == "mbu":
        return n - 1
    if convention == "2n-3":
        return 2 * n - 3
    if convention == "draft-color":
        return 2 * (n - 2)
    raise ValueError(f"unknown C^nX convention {convention!r}")


def mcx_table_toffolis(table: Mapping[int, int], convention: str = "mbu") -> int:
    """Toffolis of a C^nX gate table {n_controls: count} (CX entries, n=1, cost nothing)."""
    return sum(cnt * mcx_toffolis(n, convention) for n, cnt in table.items())


@dataclass(frozen=True)
class FermionDraftGroup:
    """One group's fermion-sector counts as the unpublished Fermion_Primitives draft prints them."""
    key: str
    draft_name: str
    n: int                                   # fiducial irrep dimension N
    col_squish_mcx: Mapping[int, int]        # color (V_g) squish C^nX table
    col_squish_t_printed: float
    col_parity_toffoli: int                  # SU(2) parity / SU(3) flag Toffolis per frame application
    col_parity_t_each_printed: float         # T per parity Toffoli in FP's constant (7; 2 for S36)
    col_rot: int                             # color rotations per frame application per field
    col_squish_copies_printed: int           # squish copies inside FP's C_V / C^G (2 SU(2), 1 SU(3))
    hop_squish_mcx: Mapping[int, int]
    hop_squish_t_printed: float
    hop_flag_toffoli: int                    # SU(3) eigen-class flag Toffolis (hopping.tex:268-281)
    n_classes: int                           # eigen-classes the hop distinguishes (hopping.tex:15-78)
    n_diag_text: int                         # diagonalizers inside FP's text C_T(hop)
    printed: Mapping[str, tuple[float, float]]   # FP's printed a + b log2(1/eps)
    src: Mapping[str, str]
    note: str = ""


FP_DIR = "work/extracted/Fermion_Primitives/publication_draft"
FP_KEY = "FermionPrimitives_unpub"
FP_C_T = (38.8, 4.45)        # diagonalizer, hopping.tex:311, as printed
FP_DIAG_T, FP_DIAG_ROT = 2, 4  # its reading: 2 T (controlled H) + 4 rotations (38.8 = 2 + 4 x 9.2)

FERMION_DRAFT: dict[str, FermionDraftGroup] = {
    "2T": FermionDraftGroup(
        "2T", "BT", 2,
        col_squish_mcx={1: 11, 2: 11, 3: 1, 4: 2, 5: 1}, col_squish_t_printed=189,
        col_parity_toffoli=44, col_parity_t_each_printed=7, col_rot=52, col_squish_copies_printed=2,
        hop_squish_mcx={1: 3, 2: 2, 3: 1, 5: 1}, hop_squish_t_printed=84, hop_flag_toffoli=0,
        n_classes=5, n_diag_text=10,
        printed={"C_V": (1164, 59.8), "C_rot": (786, 59.8), "C_G_tab": (378, 59.8), "C_hop_tab": (452, 46),
                 "C_hop_text": (556, 44.5), "C_S": (2416, 211.6)},
        src={"col": "su2_diag.tex:186-191 (BT row), :209, :225, :228, :235", "hop": "hopping.tex:28-30, :211, :338",
             "tab": "resources.tex:61, :32"},
        note="BT text C_T(hop) counts n_cls N = 10 diagonalizers (T only); BO and SU(3) count 2 n_cls N. "
             "Printed C_rot 786 = 308 + 478.4 = 786.4 (rounded)."),
    "2O": FermionDraftGroup(
        "2O", "BO", 2,
        col_squish_mcx={1: 10, 2: 15, 3: 14, 4: 11, 5: 6, 6: 1}, col_squish_t_printed=917,
        col_parity_toffoli=56, col_parity_t_each_printed=7, col_rot=60, col_squish_copies_printed=2,
        hop_squish_mcx={1: 2, 2: 1, 3: 9, 4: 5, 6: 1}, hop_squish_t_printed=434, hop_flag_toffoli=0,
        n_classes=7, n_diag_text=28,
        printed={"C_V": (2778, 69), "C_rot": (944, 69), "C_G_tab": (1834, 69), "C_hop_tab": (949.2, 64.4),
                 "C_hop_text": (1954.4, 124.6), "C_S": (9234.4, 266.8)},
        src={"col": "su2_diag.tex:186-191 (BO row), :210, :225, :229, :236", "hop": "hopping.tex:32-34, :212, :339",
             "tab": "resources.tex:62, :33"},
        note="Printed color squish 917 T = 131 Toffolis = the table at 2(n-2); at 2n-3 it is 163 (1141 T), "
             "under MBU 105 (735 T)."),
    "S36x3": FermionDraftGroup(
        "S36x3", "Sigma(36x3)", 3,
        col_squish_mcx={1: 38, 3: 10, 4: 60, 5: 90, 6: 29, 7: 9}, col_squish_t_printed=7854,
        col_parity_toffoli=778, col_parity_t_each_printed=2, col_rot=164, col_squish_copies_printed=1,
        hop_squish_mcx={1: 8, 2: 1, 3: 1, 4: 1, 5: 10}, hop_squish_t_printed=136, hop_flag_toffoli=100,
        n_classes=10, n_diag_text=60,
        printed={"C_rot": (1508.8, 188.6), "C_G": (10918.8, 188.6), "C_hop_tab": (2102, 133.5),
                 "C_hop_text": (4000, 267), "C_S": (26041.6, 644.2)},
        src={"col": "su3_diag.tex:89, :103, :124, :132, :196, :211", "hop": "hopping.tex:47-61, :262, :277, :345",
             "tab": "resources.tex:65, :39"},
        note="C^G prices the 778 flag Toffolis at 2 T (S72 and Delta(27) use 7). Printed hop-squish T 136 is not an "
             "integer number of Toffolis; FP comments out 469 (= the table at 2(n-2)). n_cls = 10 (Delta(27) 4 + "
             "Delta(54) 3 + Sigma(36x3) 3)."),
    "S72x3": FermionDraftGroup(
        "S72x3", "Sigma(72x3)", 3,
        col_squish_mcx={1: 360, 2: 3, 3: 12, 4: 73, 5: 123, 6: 95, 7: 1}, col_squish_t_printed=12789,
        col_parity_toffoli=892, col_parity_t_each_printed=7, col_rot=184, col_squish_copies_printed=1,
        hop_squish_mcx={1: 10, 3: 2, 4: 9, 5: 8}, hop_squish_t_printed=372, hop_flag_toffoli=100,
        n_classes=10, n_diag_text=60,
        printed={"C_rot": (1692.8, 211.6), "C_G": (20725.8, 211.6), "C_hop_tab": (2396, 133.5),
                 "C_hop_text": (4472, 267), "C_S": (46243.6, 690.2)},
        src={"col": "su3_diag.tex:90, :104, :125, :133, :197, :212", "hop": "hopping.tex:263, :278, :346",
             "tab": "resources.tex:66, :40"},
        note="FP lists no Sigma(72x3)-specific eigen-classes (hopping.tex:39-78 has none); its text C_T(hop) uses the "
             "Sigma(36x3) count, 30 = 10 x 3. Printed hop-squish T 372 is not integer Toffolis; FP comments out 616 "
             "(= the table at 2(n-2))."),
}


def fp_printed(group: str, quantity: str, eps: float) -> float:
    """FP's printed a + b log2(1/eps) for `quantity` (keys of FermionDraftGroup.printed)."""
    a, b = FERMION_DRAFT[group].printed[quantity]
    return a + b * _log(eps)


def fp_reconstruct(group: str, quantity: str) -> tuple[float, float]:
    """FP's printed constants rebuilt from its own component counts in its own conventions (no report rule):
    rotations at 9.2 + 1.15 L (the slope as FP prints it), color squish at 2(n-2), hop squish at 2n-3 (SU(2)) or
    its printed T (SU(3)). Returns (a, b). Used by the tests to show which printed numbers are self-consistent."""
    g = FERMION_DRAFT[group]
    if quantity == "col_squish_t":
        return 7.0 * mcx_table_toffolis(g.col_squish_mcx, "draft-color"), 0.0
    if quantity == "hop_squish_t":
        if g.n == 2:
            return 7.0 * mcx_table_toffolis(g.hop_squish_mcx, "2n-3"), 0.0
        return 7.0 * mcx_table_toffolis(g.hop_squish_mcx, "draft-color"), 0.0   # the commented-out 469 / 616
    sq = g.col_squish_t_printed
    rot_a, rot_b = 9.2 * g.col_rot, 1.15 * g.col_rot
    par = g.col_parity_toffoli * g.col_parity_t_each_printed
    if quantity == "C_rot":
        return par + rot_a if g.n == 2 else rot_a, rot_b
    if quantity in ("C_V", "C_G"):
        return g.col_squish_copies_printed * sq + par + rot_a, rot_b
    if quantity == "C_hop_text":
        a = 2 * (g.hop_squish_t_printed + 7 * g.hop_flag_toffoli) + g.n_diag_text * FP_C_T[0]
        return a, g.n_diag_text * FP_C_T[1]
    raise KeyError(quantity)


def spinor_counts(d: int, n: int, mcx: str = "mbu") -> dict:
    """Wilson spinor diagonalization per link per field, all N colors: FP resources.tex:60-67
    (280 floor((d-1)/2) + 16 T per color) split as spin_diag.tex:282-297. Under MBU the 10 Toffolis per W gate
    (5 triply-controlled rotations x 2 compute + 2 uncompute, spin_diag.tex:280) become 5."""
    if d <= 2:
        return {"toffoli": 0, "t_direct": 16 * n, "n_rot": 0}
    per_w = 10 if mcx != "mbu" else 5
    return {"toffoli": 4 * per_w * n, "t_direct": 16 * n, "n_rot": 0}


def _hwp_counts(k: int) -> tuple[int, int]:
    from .common import hwp_synth_rotations, hwp_toffolis
    return hwp_synth_rotations(k), hwp_toffolis(k)


def fermion_hop_counts(group: str, n_stag: int = 1, *, fermion: str = "staggered", d: int | None = None,
                       n_c: int = 3, k: int | None = None, undo: bool = True, share: str = "link",
                       mcx: str = "mbu", phasing: str = "hwp", n_spin: int | None = None,
                       frame_rot_per_field: bool = True) -> dict:
    """Gate counts (eps-free) of one gauge-covariant hop on one link per Trotter step: Toffolis, direct T gates,
    synthesized rotations, with a breakdown. See the block comment above for every reading.

    group     "2T", "2O", "S36x3", "S72x3" (FP counts) or "Z3" (R-HOP structure; n_c, k used there only).
    fermion   "staggered" or "wilson" (needs d). n_spin: spinor components that hop after V_S; default 1
              (staggered) and n_D / 2 at r = 1 (Wilson, n_D = 2^{floor((d+1)/2)}: the projector 1 - gamma^0 keeps
              half; FP section_hamiltonian.tex:77-80). FP's W1 has no spinor multiplicity; n_spin = 2 at d = 3
              (one color frame per hopping spinor component) is AUTHOR-CONFIRMED (E21 (2), 2026-10-01). V_S and V_g
              act on different indices and commute, so only the hopping components need the color frame.
    undo      count V_g^dag on both sites (default True). FP draws V_C^dag (big_fig_group_agnostic.pdf,
              section_hamiltonian.tex:86-91; running_notes.tex:124-128, :213-225) but never counts it: its
              factor 2 / "2 C^G" is two sites, forward only (su2_diag.tex:238; resources.tex:15), SU(3) C^G is one
              forward application (su3_diag.tex:130-134), and the printed totals are 2 C^G + 2 C^hop(table)
              exactly (resources.tex:37-41, e.g. 46243.6 = 2 x 20725.8 + 2 x 2396). False = FP's literal 2 C^G.
    share     SEPARATE assumption from the undo (the g-only squish / parity multiplicity). "link" (default since
              r19, author ruling E21 (1)): color squish + parity and hop squish + flags computed once per link
              and held through V(x), V(y), hop, V^dag(x), V^dag(y). "draft" (conservative sensitivity, the r17
              headline): color squish + parity per frame application and hop squish + flags per field.
    mcx       C^nX convention for squish ladders: "mbu" (report rule), "2n-3", "draft-color".
    phasing   "hwp" (one group of 2 N n_spin n_stag), "plain" (that many rotations), "none" (block encoding).
    frame_rot_per_field  False reproduces the literal "4 C^G shared by the fields" reading (rotations not
              multiplied by n_stag); recorded only.
    """
    if n_stag < 1:
        raise ValueError("n_stag >= 1")
    if share not in ("draft", "link"):
        raise ValueError(f"unknown share {share!r}")
    if phasing not in ("hwp", "plain", "none"):
        raise ValueError(f"unknown phasing {phasing!r}")
    items: list[dict] = []

    def add(name, toffoli=0, t_direct=0, n_rot=0, status=CircuitStatus.COMPILED, src="", note=""):
        items.append({"name": name, "toffoli": toffoli, "t_direct": t_direct, "n_rot": n_rot,
                      "status": status, "src": src, "note": note})

    if group == "Z3":
        if fermion != "staggered":
            raise ValueError("Z3: staggered only (no FP entry, R-HOP structure)")
        kk = n_c * n_stag if k is None else int(k)
        r, t = _hwp_counts(kk)
        npl = RHOP_N_PAULI_DIAGONAL["Z3"]
        add("z3_hop_hwp", toffoli=npl * t, n_rot=npl * r, status=CircuitStatus.SCALING, src=RHOP_SRC,
            note=f"no FP entry; R-HOP structure: C_W = 0 (diagonal), {npl} strings x HWP(k={kk}); rotations at the "
                 "full fit (r17 rule 6)")
        out = {"group": group, "fermion": fermion, "n_stag": n_stag, "k": kk, "items": items}
    else:
        if group not in FERMION_DRAFT:
            raise KeyError(f"no Fermion_Primitives counts for {group!r}")
        g = FERMION_DRAFT[group]
        if fermion == "staggered":
            ns = 1 if n_spin is None else int(n_spin)
        elif fermion == "wilson":
            if d is None:
                raise ValueError("Wilson needs d")
            n_dirac = 2 ** ((d + 1) // 2)              # 2 at d = 1, 2; 4 at d = 3
            ns = n_dirac // 2 if n_spin is None else int(n_spin)   # r = 1: 1 - gamma^0 keeps half
        else:
            raise ValueError(f"unknown fermion {fermion!r}")
        n_app = 4 if undo else 2
        src_fp = f"{FP_KEY} ({FP_DIR})"
        # color frame: g-only part
        sq = mcx_table_toffolis(g.col_squish_mcx, mcx)
        frame_shared = 2 * sq + g.col_parity_toffoli
        n_shared = n_app if share == "draft" else 1
        add("colour_squish_and_parity", toffoli=n_shared * frame_shared, src=f"{src_fp} {g.src['col']}",
            note=(f"{n_shared} x (2 x {sq} squish [{mcx}] + {g.col_parity_toffoli} parity/flag) Toffolis; "
                  + ("share='draft' sensitivity: FP's per-site compute + uncompute (su2_diag.tex:238) at every "
                     "frame application" if share == "draft"
                     else "share='link' (E21): computed once per link, held through V(x), V(y), hop, V^dag(x), "
                          "V^dag(y), then uncomputed")
                  + ("; SU(3) squish uncompute ESTIMATED (FP's C^G has one squish, su3_diag.tex:130-134)"
                     if g.n == 3 else "")),
            status=CircuitStatus.COMPILED if not undo else CircuitStatus.SCALING)
        rot_sets = n_app * (n_stag if frame_rot_per_field else 1) * ns
        add("colour_rotations", n_rot=rot_sets * g.col_rot, src=f"{src_fp} {g.src['col']}",
            note=f"{rot_sets} sets (= {n_app} frame applications x {n_stag if frame_rot_per_field else 1} fields x "
                 f"{ns} spinor) x {g.col_rot} rotations"
                 + ("; the V_g^dag half is drawn in FP (big_fig_group_agnostic.pdf) but not in its 2 C^G "
                    "(su2_diag.tex:238, resources.tex:15, :37-41)" if undo else "; FP's literal 2 C^G (no undo)"),
            status=CircuitStatus.SCALING if undo else CircuitStatus.COMPILED)
        # hop: g-only part
        hsq = mcx_table_toffolis(g.hop_squish_mcx, mcx)
        hop_shared = 2 * (hsq + g.hop_flag_toffoli)
        n_hshared = n_stag if share == "draft" else 1
        add("hop_squish_and_flags", toffoli=n_hshared * hop_shared, src=f"{src_fp} {g.src['hop']}",
            note=f"{n_hshared} x 2 x ({hsq} squish [{mcx}] + {g.hop_flag_toffoli} flags): compute + uncompute "
                 "(FinalDiagonalizations.tex)")
        n_diag = 2 * g.n_classes * g.n * ns * n_stag
        add("hop_diagonalizers", t_direct=FP_DIAG_T * n_diag, n_rot=FP_DIAG_ROT * n_diag,
            src=f"{src_fp} hopping.tex:309-311; FinalDiagonalizations.tex",
            note=f"{n_diag} = 2 (T, T^dag) x {g.n_classes} classes x N={g.n} x {ns} spinor x {n_stag} fields; "
                 "each 2 T + 4 rotations (C(T) = 38.8 + 4.45 L read as 2 + 4 x (9.2 + 1.15 L))")
        n_ph = 2 * g.n * ns * n_stag
        if phasing == "hwp":
            r, t = _hwp_counts(n_ph)
            add("hop_phasing_hwp", toffoli=t, n_rot=r, status=CircuitStatus.SCALING,
                src="FinalDiagonalizations.tex; arxiv_1709_06648,arxiv_1902_10673",
                note=f"Rz(-th), Rz(th) per pair: {n_ph} equal-angle rotations as one HWP group (derived here)")
        elif phasing == "plain":
            add("hop_phasing_plain", n_rot=n_ph, src="FinalDiagonalizations.tex", note=f"{n_ph} rotations")
        if fermion == "wilson":
            sp = spinor_counts(d, g.n, mcx)
            add("spinor_frame", toffoli=n_stag * sp["toffoli"], t_direct=n_stag * sp["t_direct"],
                src=f"{src_fp} resources.tex:60-67; spin_diag.tex:280-297",
                note=f"N C^spinor(d={d}) per field, as FP's W1 (resources.tex:15-16) counts it")
        out = {"group": group, "fermion": fermion, "n_stag": n_stag, "d": d, "n_spin": ns, "undo": undo,
               "share": share, "mcx": mcx, "phasing": phasing, "items": items}
    out["toffoli"] = sum(i["toffoli"] for i in items)
    out["t_direct"] = sum(i["t_direct"] for i in items)
    out["n_rot"] = sum(i["n_rot"] for i in items)
    return out


def hop_link_cost(group: str, n_stag: int = 1, d: int | None = None, n_rot_circuit: float | None = None, *,
                  eps: float | None = None, synthesis: str = "rus", toffoli_convention: str = "textbook",
                  legacy: bool = True, **kw) -> dict:
    """T per link per Trotter step of the gauge-covariant hop, priced under the report's rules (r17).

    Counts from fermion_hop_counts (keyword arguments pass through). The rotation tolerance is `eps` if given,
    else R-TOL at the circuit's own rotation count: eps = sqrt(1e-2 / n_rot_circuit). Rotations at the full
    fit 1.15 log2(1/eps) + 9.2 (synthesis="rus"); 7 T per Toffoli. With legacy=True the result also carries
    the retired R-HOP price at the same eps (both at the chapters' old "rus-slope" and at the full fit) and the
    ratio new / old.

    Defaults (r19, author ruling E21 (1)): undo=True (V_g^dag on both sites) and share="link" (color squish +
    parity and hop squish + flags once per link, held through V(x), V(y), hop, V^dag(x), V^dag(y)). The r17
    per-application recompute stays reachable as the conservative sensitivity: share="draft".

    Returns the counts dict plus {"eps", "t_rot", "t_toffoli", "t_rotations", "t", "legacy"}.
    """
    from .common import T_PER_TOFFOLI, eps_rot_for, t_per_rotation
    if eps is None:
        if n_rot_circuit is None:
            raise ValueError("give eps or n_rot_circuit (R-TOL)")
        eps = eps_rot_for(n_rot_circuit)
    c = fermion_hop_counts(group, n_stag, d=d, **kw)
    t_rot = t_per_rotation(eps, synthesis)
    c["eps"], c["t_rot"] = eps, t_rot
    c["t_toffoli"] = c["toffoli"] * T_PER_TOFFOLI[toffoli_convention]
    c["t_rotations"] = c["n_rot"] * t_rot
    c["t"] = c["t_toffoli"] + c["t_direct"] + c["t_rotations"]
    for i in c["items"]:
        i["t"] = i["toffoli"] * T_PER_TOFFOLI[toffoli_convention] + i["t_direct"] + i["n_rot"] * t_rot
    if legacy and kw.get("fermion", "staggered") == "staggered":
        n_c = FERMION_DRAFT[group].n if group in FERMION_DRAFT else kw.get("n_c", 3)
        old = rhop_link(group, n_stag, n_c, eps, k=kw.get("k"), synthesis="rus-slope")
        old_full = rhop_link(group, n_stag, n_c, eps, k=kw.get("k"), synthesis="rus")
        c["legacy"] = {"rhop_rus_slope": old["t"], "rhop_rus": old_full["t"],
                       "ratio_new_over_rhop_rus_slope": c["t"] / old["t"],
                       "ratio_new_over_rhop_rus": c["t"] / old_full["t"]}
    return c


def staggered_mass_site(n_c: int, n_stag: int, eps: float | None = None, n_rot_circuit: float | None = None, *,
                        synthesis: str = "rus", toffoli_convention: str = "textbook") -> dict:
    """Staggered mass per site per Trotter step: one R_Z per color-flavor copy (FP section_mass.tex:3-10,
    resources.tex:5-7), the N_c N_stag equal-angle copies grouped as one HWP group (Ch. 6's ruled form; option M-B).
    A baryon chemical potential mu_B N_B commutes with H and adds the same (-mu_B/3) dt to every copy on a site,
    so it folds into this angle at no extra cost. eps as hop_link_cost (or R-TOL from n_rot_circuit)."""
    from .common import eps_rot_for, hwp_group
    if eps is None:
        if n_rot_circuit is None:
            raise ValueError("give eps or n_rot_circuit (R-TOL)")
        eps = eps_rot_for(n_rot_circuit)
    h = hwp_group(int(n_c) * int(n_stag), eps, synthesis, toffoli_convention)
    h["src"] = f"{FP_KEY} section_mass.tex:3-10; HWP arxiv_1709_06648,arxiv_1902_10673"
    h["eps"] = eps
    return h


def wilson_mass_rotations(n_c: int, n_f: int, d: int) -> int:
    """Wilson on-site term per site per Trotter step, plain R_Z count: 2 N (d = 1, 2) or 4 N (d = 3, 4),
    N = N_c N_f (FP section_mass.tex:12). FP's resources.tex:24 total with floor(d/2) disagrees at d = 3
    (2 N instead of 4 N) and is not used. Callers may group equal angles by HWP."""
    n = int(n_c) * int(n_f)
    return 2 * n if d <= 2 else 4 * n


def conflicts() -> list[tuple[str, int, int]]:
    """(group, chapter_width, compiled_width) for every link-width disagreement."""
    return [(k, g.link_qubits, int(g.link_qubits_compiled.lo))
            for k, g in GROUPS.items() if g.link_width_conflict]
