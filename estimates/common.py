"""Shared conventions for the chapter resource-estimate models.

One home for every number that more than one chapter uses, so a change here
(a synthesis constant, a Toffoli convention) propagates and cannot drift.
Everything is standard-library only.

Provenance is mandatory. Every input to a chapter model is a `Tagged` value
carrying one of four provenance tags, and every primitive in a T-count carries
a `CircuitStatus`. The point is not decoration: the `__main__` table prints
the tag beside every input, and `CIRCUIT_STATUS.md` is generated from the
statuses, so a reader sees at a glance what is derived, what is asserted, and
what would have to be compiled before any of this could be handed to a
logical-circuit tool such as QLX.

Rotation synthesis follows T4/T4b in TRACKED_CHANGES.md: the repeat-until-
success (RUS) fit of Bocharov, Roetteler & Svore, with errors adding
incoherently under randomized synthesis (Campbell 2017) unless a model says
otherwise. Chapters 2, 3 and 9 all state this model; encode it here, not there.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Callable, Iterable, Mapping


# --------------------------------------------------------------------------- #
# Provenance
# --------------------------------------------------------------------------- #

class Provenance(str, Enum):
    CITED = "cited"                    # a published source gives this value
    ASSUMED = "assumed"                # a working assumption the chapter declares
    STATED = "stated-not-derived"      # the chapter asserts it without derivation
    UNCITED = "uncited"                # asserted with no source at all


@dataclass(frozen=True)
class Tagged:
    """A model input with its provenance. `value` may be a number or (lo, hi)."""
    value: Any
    prov: Provenance
    src: str = ""          # bib key, arXiv id, or "appNN:line"
    note: str = ""

    def __float__(self) -> float:
        return float(self.lo)

    @property
    def lo(self):
        return self.value[0] if isinstance(self.value, tuple) else self.value

    @property
    def hi(self):
        return self.value[1] if isinstance(self.value, tuple) else self.value

    @property
    def is_range(self) -> bool:
        return isinstance(self.value, tuple)


def Cited(value, src: str, note: str = "") -> Tagged:
    return Tagged(value, Provenance.CITED, src, note)


def Assumed(value, note: str = "", src: str = "") -> Tagged:
    return Tagged(value, Provenance.ASSUMED, src, note)


def Stated(value, src: str, note: str = "") -> Tagged:
    """Chapter states the number but shows no derivation. `src` = appNN:line."""
    return Tagged(value, Provenance.STATED, src, note)


def Uncited(value, note: str = "") -> Tagged:
    return Tagged(value, Provenance.UNCITED, "", note)


# --------------------------------------------------------------------------- #
# Circuit knowledge
# --------------------------------------------------------------------------- #

class CircuitStatus(str, Enum):
    """How well the gate-level circuit behind a primitive cost is known."""
    COMPILED = "COMPILED"        # explicit gate-level circuit exists in the source
    SCALING = "SCALING"          # counting/scaling argument; no circuit
    CONJECTURE = "CONJECTURE"    # order-of-magnitude guess or conjecture by analogy
    UNSOURCED = "UNSOURCED"      # a number with no stated origin


@dataclass(frozen=True)
class Primitive:
    """One named line item in a T-count breakdown.

    `count` × `t_each` = its T contribution. For synthesized rotations, set
    `t_each` from `t_per_rotation()` and record `eps_rot` in `note`.
    """
    name: str
    count: float
    t_each: float
    status: CircuitStatus
    src: str = ""
    note: str = ""

    @property
    def t_total(self) -> float:
        return self.count * self.t_each


# --------------------------------------------------------------------------- #
# Rotation synthesis (T4 / T4b)
# --------------------------------------------------------------------------- #

RUS_SLOPE = 1.15        # T per log2(1/eps), Bocharov-Roetteler-Svore fit
RUS_OFFSET = 9.2        # additive constant in the same fit
RS_SLOPE = 3.0          # Ross-Selinger ancilla-free, ~3 log2(1/eps): the coherent method
RUS_ANCILLA = 1         # ancilla held while one RUS rotation is in flight (Bocharov-Roetteler-Svore use one);
                        # author ruling E26, r22, 2026-10-01. Chapters had counted 2

SYNTHESIS_SRC = "bocharovRoettelerSvore2015,campbell2017"


def t_per_rotation(eps: float, model: str = "rus") -> float:
    """Expected T-gates to synthesize one R_Z to tolerance `eps`.

    model="rus"        : 1.15*log2(1/eps) + 9.2   (BRS's own fit; what the chapters use)
    model="watson"     : 1.15*log2(2/eps) + 9.2   (Watson et al. Eq. 80; Perdue's T4 note)
    model="rus-slope"  : 1.15*log2(1/eps)         (su3qlx / the SU(3) papers' convention; LEGACY records only:
                                                    E20, 2026-10-01, prices every rotation at "rus")
    model="rs"         : 3*log2(1/eps)            (Ross-Selinger, ancilla-free)

    At eps=1e-4: rus -> 24.5, watson -> 25.7, rus-slope -> 15.3, rs -> 39.9. The
    chapters' "~30 T at 1e-4" is the RUS figure rounded up; see T4 in TRACKED_CHANGES.md.
    """
    l2 = math.log2(1.0 / eps)
    if model == "rus":
        return RUS_SLOPE * l2 + RUS_OFFSET
    if model == "watson":
        return RUS_SLOPE * (l2 + 1.0) + RUS_OFFSET
    if model == "rus-slope":
        return RUS_SLOPE * l2
    if model == "rs":
        return RS_SLOPE * l2
    raise ValueError(f"unknown synthesis model {model!r}")


def eps_per_rotation(eps_total: float, n_rot: float, errors: str = "incoherent") -> float:
    """Per-rotation tolerance that holds the total synthesis error to `eps_total`.

    errors="incoherent": randomized synthesis, errors add in quadrature, eps ~ sqrt(eps_total/N)
    errors="coherent"  : deterministic synthesis, errors add linearly,     eps ~ eps_total/N
    """
    if n_rot <= 0:
        raise ValueError("n_rot must be positive")
    if errors == "incoherent":
        return math.sqrt(eps_total / n_rot)
    if errors == "coherent":
        return eps_total / n_rot
    raise ValueError(f"unknown error model {errors!r}")


# Ruling R-TOL (H. Lamm, 2026-09-29; TRACKED_CHANGES.md "RULINGS, R-TOL and overshoot"): the total
# synthesis error per shot is eps_syn = 1e-2, report-wide, and each circuit sets its own per-rotation
# tolerance from its size. Replaces every fixed per-rotation tolerance (1e-4, 1e-3, and the tolerances
# behind "~30 T" and c_T = 27). Each chapter keeps its own synthesis formula; only eps changes.
EPS_SYN = 1e-2
EPS_SYN_SRC = "TRACKED_CHANGES.md R-TOL (H. Lamm, 2026-09-29)"


def eps_rot_for(n_rot: float, eps_syn: float = EPS_SYN) -> float:
    """Per-rotation tolerance of a circuit under ruling R-TOL: sqrt(eps_syn / n_rot).

    `n_rot` is the number of synthesized (arbitrary-angle) rotations in ONE shot of the circuit;
    Clifford+T-exact rotations (multiples of pi/4) are not synthesized and are not counted.
    Randomized synthesis, errors add as n_rot eps^2 (eps_per_rotation, "incoherent"). n_rot is
    fixed by the circuit, not by eps, so no iteration is needed.
    """
    return eps_per_rotation(eps_syn, n_rot, "incoherent")


def synthesized_rotations(n_rot: float, eps_rot: float, model: str = "rus",
                          src: str = SYNTHESIS_SRC, note: str = "") -> Primitive:
    """Convenience: the Primitive for `n_rot` synthesized rotations at `eps_rot`."""
    return Primitive("rot_synth", n_rot, t_per_rotation(eps_rot, model),
                     CircuitStatus.COMPILED, src,
                     note or f"eps_rot={eps_rot:.1e}, model={model}")


# --------------------------------------------------------------------------- #
# Toffoli and multi-controlled gates
# --------------------------------------------------------------------------- #

T_PER_TOFFOLI = {
    "textbook": 7,      # standard Clifford+T decomposition; su3qlx and the SU(3) papers
    "jones":    4,      # Jones 2013 with a measurement/ancilla; used by Ch. 9 via arXiv:1709.06648
}


def toffoli_t(n_toffoli: float, convention: str) -> float:
    """T-count of `n_toffoli` Toffolis. A model must name the convention."""
    return n_toffoli * T_PER_TOFFOLI[convention]


def mcx_toffoli_count(n_controls: int, convention: str = "linear") -> int:
    """Toffolis in a C^nNOT. Copied from su3qlx/gates.py.

    convention="linear": 2n-3 (clean ancilla)        convention="log": 2*ceil(log2 n)-1 (dirty)
    """
    if n_controls < 1:
        raise ValueError("need >=1 control")
    if n_controls == 1:
        return 0
    if n_controls == 2:
        return 1
    if convention == "linear":
        return 2 * n_controls - 3
    if convention == "log":
        return 2 * math.ceil(math.log2(n_controls)) - 1
    raise ValueError(convention)


# --------------------------------------------------------------------------- #
# Hamming-weight phasing (HWP): the building block of rule R-HOP (groups.rhop_link)
# --------------------------------------------------------------------------- #
# k equal-angle Z rotations applied as floor(log2 k)+1 synthesized R_Z on the Hamming-weight
# register plus k - w(k) adder Toffolis, the adders uncomputed by measurement and a classically
# controlled Clifford (no T, no extra qubits) [arxiv_1709_06648, arxiv_1902_10673].
# The same three counts as estimates.ch09_chiral_gauge.hwp_* (app06:108), which Ch. 9 still
# defines for itself; tests/test_rhop.py asserts the two agree for every k it checks.
#
# ANCILLA (author ruling E26, r22, 2026-10-01; arxiv_1902_10673 App. A, arxiv_1709_06648): k - w(k), the same
# number as the Toffolis. The weight is computed in place by a tree of half and full adders. Each adder is one
# Toffoli, leaves its sum bit on one of its input qubits and writes its carry to one new ancilla. The bits of
# the weight above the lowest ARE those carry ancilla; the lowest bit sits on a data qubit. There is no separate
# weight register to add. Until r22 hwp_ancilla returned toffolis + floor(log2 k) + 1, which counted the weight
# register twice (k=32: 37 instead of 31; k=9: 11 instead of 7).

HWP_SRC = "arxiv_1709_06648,arxiv_1902_10673"


def hwp_synth_rotations(k: int) -> int:
    """Synthesized R_Z in one HWP group of k equal-angle rotations: floor(log2 k)+1."""
    if int(k) < 1:
        raise ValueError("an HWP group needs k >= 1")
    return math.floor(math.log2(k)) + 1


def hwp_toffolis(k: int) -> int:
    """Adder Toffolis to compute the Hamming weight of k bits: k - w(k), w = binary popcount."""
    if int(k) < 1:
        raise ValueError("an HWP group needs k >= 1")
    return int(k) - bin(int(k)).count("1")


def hwp_ancilla(k: int) -> int:
    """Ancilla of one HWP group of k equal-angle rotations: k - w(k), one carry ancilla per adder Toffoli
    (arxiv_1902_10673 App. A; arxiv_1709_06648). The floor(log2 k)+1-bit weight is held on those ancilla, with
    its lowest bit on a data qubit, so no weight register is added (E26, r22). k=3 -> 1, 8 -> 7, 9 -> 7,
    12 -> 10, 32 -> 31, 192 -> 190. All of them are live while the group's rotations are synthesized."""
    return hwp_toffolis(k)


def hwp_group(k: int, eps: float, synthesis: str = "rus",
              toffoli_convention: str = "textbook") -> dict:
    """T-cost of one HWP group of k equal-angle rotations synthesized to tolerance `eps`.

    Default convention is the report's (r17, 2026-10-01): the full RUS fit 1.15 log2(1/eps) + 9.2 per
    rotation ("rus") and 7 T per Toffoli ("textbook"). Until r17 the default was the slope-only
    "rus-slope" (15.3 T at 1e-4 instead of 24.5), which no ruling adopted for fermion terms; every
    chapter model passes `synthesis` explicitly, so the change moves no chapter number. Returns
    {"k", "n_rot", "n_toffoli", "ancilla", "t_rot", "t"} with t = n_rot t_rot + T(n_toffoli).
    """
    t_rot = t_per_rotation(eps, synthesis)
    n_rot, n_tof = hwp_synth_rotations(k), hwp_toffolis(k)
    return {"k": k, "n_rot": n_rot, "n_toffoli": n_tof, "ancilla": hwp_ancilla(k), "t_rot": t_rot,
            "t": n_rot * t_rot + toffoli_t(n_tof, toffoli_convention)}


# --------------------------------------------------------------------------- #
# T-depth, factory parallelism and wall-time exports (ruling R9, H. Lamm, 2026-10-02)
# --------------------------------------------------------------------------- #
# Factory model (factory.json synthesis; Gidney-Ekera, arxiv_1905_09749). One factory delivers one T state per
# logical cycle t_c = 10 us; each sequential non-Clifford layer waits one reaction time t_r = 10 us. With F
# factories a shot takes max(N_T t_c / F, D_T t_r). The report's wall-time convention (1 us per T + 0.1 ms per shot)
# is the F = 10 baseline: 10 factories at 10 us each deliver one T per us. Above F* = N_T / D_T extra factories sit
# idle, so the shortest possible run is shots x D_T x t_r ("floor"). A run fits in one year on ONE machine only if
# F_1yr = 10 x (serial years at 1 us/T) <= F*.
#
# No central override of gate time or shot overhead: each chapter keeps its own t_gate_s and shot_overhead_s in its
# Assumptions and passes them in. These helpers only fix the arithmetic and the names.

REACTION_TIME_S = 10e-6        # t_r: one reaction time per sequential non-Clifford layer, arxiv_1905_09749
FACTORY_BASELINE = 10          # factories behind the report's 1 us per T convention (10 x one T per 10 us)
FACTORY_SRC = "arxiv_1905_09749"
SECONDS_PER_YEAR = 365.25 * 86400.0   # Julian year, as ch03/ch04/ch07/ch10 (chapters may keep their own constant)

DEPTH_EXPORTS = ("t_per_shot", "t_depth_per_shot", "f_star", "floor_wall_s", "factories_for_1yr",
                 "wall_first_result_s", "wall_campaign_s")


def _band(x) -> tuple[float, float]:
    """A number, a (lo, hi) tuple, or a Tagged value -> (lo, hi) floats."""
    if isinstance(x, Tagged):
        x = x.value
    if isinstance(x, tuple):
        lo, hi = float(x[0]), float(x[1])
    else:
        lo = hi = float(x)
    if hi < lo:
        raise ValueError(f"range is reversed: {lo} > {hi}")
    return lo, hi


def depth_exports(t_per_shot, t_depth_per_shot, shots, t_gate_s, shot_overhead_s,
                  factory_baseline: float = FACTORY_BASELINE,
                  reaction_time_s: float = REACTION_TIME_S) -> dict:
    """The R9 depth/factory numbers for ONE run (one point set, one shot count) on ONE machine.

    Arguments may be numbers, (lo, hi) tuples, or Tagged values. `t_gate_s` and `shot_overhead_s` are the
    chapter's own Assumptions values (1e-6 and 1e-4 by the report convention); nothing here overrides them.
    T-depth comes from factory.json (low/high as a range).

    Returns (every value a (lo, hi) tuple of floats):
      t_per_shot, t_depth_per_shot   the inputs, as bands
      f_star             N_T / D_T, cross-paired (t_lo / d_hi, t_hi / d_lo) so the band is the widest honest one
      floor_wall_s       shots x D_T x t_r (t_r = 10 us): the shortest run any number of factories allows
      wall_serial_s      shots x (N_T x t_gate_s + shot_overhead_s): the report's baseline wall on one machine
      serial_yr          wall_serial_s in years
      factories_for_1yr  factory_baseline x (shots x N_T x t_gate_s in years). The shot overhead is left out,
                         because it does not shrink with more factories; at 0.1 ms per shot it matters only for
                         circuits under ~1e3 T.
      baseline_ok        (f_star_lo >= factory_baseline, f_star_hi >= factory_baseline): whether the 10-factory
                         baseline can be fed at each end of the band (False -> quote a corrected wall or a fix)
      fits_1yr           (F_1yr_hi <= f_star_lo, F_1yr_lo <= f_star_hi): one year reachable at the worst / best end
    """
    t, d, n = _band(t_per_shot), _band(t_depth_per_shot), _band(shots)
    tg, t0 = _band(t_gate_s)[0], _band(shot_overhead_s)[0]
    if min(t[0], d[0], n[0]) <= 0 or tg <= 0 or t0 < 0:
        raise ValueError("T-count, T-depth, shots and t_gate_s must be positive; shot_overhead_s non-negative")
    if d[1] > t[1] * (1 + 1e-12):
        raise ValueError(f"T-depth {d} exceeds T-count {t}")
    f_star = (t[0] / d[1], t[1] / d[0])
    floor = (n[0] * d[0] * reaction_time_s, n[1] * d[1] * reaction_time_s)
    serial = (n[0] * (t[0] * tg + t0), n[1] * (t[1] * tg + t0))
    f1 = tuple(factory_baseline * n_ * t_ * tg / SECONDS_PER_YEAR for n_, t_ in ((n[0], t[0]), (n[1], t[1])))
    return {
        "t_per_shot": t,
        "t_depth_per_shot": d,
        "f_star": f_star,
        "floor_wall_s": floor,
        "wall_serial_s": serial,
        "serial_yr": (serial[0] / SECONDS_PER_YEAR, serial[1] / SECONDS_PER_YEAR),
        "factories_for_1yr": f1,
        "baseline_ok": (f_star[0] >= factory_baseline, f_star[1] >= factory_baseline),
        "fits_1yr": (f1[1] <= f_star[0], f1[0] <= f_star[1]),
    }


# --------------------------------------------------------------------------- #
# Model contract
# --------------------------------------------------------------------------- #

ERAS = ("2028", "2033", "codesign")


@dataclass(frozen=True)
class Result:
    """What every chapter model returns for one era.

    lq / hard_ops are (lo, hi); lo == hi when the chapter quotes one value.
    `breakdown` is the list of Primitives whose t_total sums to hard_ops for the
    representative point; `intermediates` holds every named number the prose
    states on the way (e.g. "700 x 3e4 x 390"), so a verifier can diff each.
    """
    era: str
    lq: tuple[float, float]
    hard_ops: tuple[float, float]
    breakdown: tuple[Primitive, ...] = ()
    intermediates: Mapping[str, Any] = field(default_factory=dict)
    shots: tuple[float, float] | None = None
    wall_time_s: tuple[float, float] | None = None
    epsilon_l: tuple[float, float] | None = None
    notes: tuple[str, ...] = ()

    def breakdown_total(self) -> float:
        return sum(p.t_total for p in self.breakdown)


@dataclass(frozen=True)
class Published:
    """The number the chapter box actually prints, with where it prints it."""
    lq: tuple[float, float]
    hard_ops: tuple[float, float] | None
    src: str                      # "app01:2033 box, line 224"
    rel_tol: float = 0.10         # how close the model must land


def close(model: tuple[float, float], published: tuple[float, float], rel_tol: float) -> bool:
    """Model band vs published band.

    Range vs range: both ends within rel_tol. Single published value ("~1.3e5")
    vs a model band: passes if the value lies inside the band (with rel_tol slack
    at the ends) -- a box that rounds a derived range to one number is not a
    disagreement with the range.
    """
    mlo, mhi = model
    plo, phi = published
    if plo == phi and mlo != mhi:
        return mlo * (1 - rel_tol) <= plo <= mhi * (1 + rel_tol)
    return all(math.isclose(m, p, rel_tol=rel_tol) for m, p in zip(model, published))


# --------------------------------------------------------------------------- #
# Formatting, shared with verify_table.py so the table and the model agree
# --------------------------------------------------------------------------- #

def sci(x: float, digits: int = 1) -> tuple[str, int]:
    """Mantissa string and exponent, chapter style: 1.1e8 -> ("1.1", 8)."""
    if x == 0:
        return "0", 0
    e = math.floor(math.log10(abs(x)))
    m = round(x / 10 ** e, digits)
    if m >= 10:
        m /= 10
        e += 1
    return f"{m:g}", e


def _mantissa(x: float, digits: int) -> tuple[float, int]:
    """(mantissa rounded to `digits` decimals, exponent), renormalised if it rounds to 10."""
    e = math.floor(math.log10(abs(x)))
    m = round(x / 10 ** e, digits)
    if m >= 10:
        m /= 10
        e += 1
    return m, e


def tex_hard_ops(lo: float, hi: float | None = None, digits: int = 1, rel: str = "",
                 keep_zeros: bool = False) -> str:
    """THE spelling of a hard-operation count in Table 1.1 (ruling R10-e, 2026-09-29).

    One value, or a range whose ends round to the same number:
        $5{\\times}10^{4}$          and, when the mantissa is 1, the bare power $10^{9}$
    Range, equal exponents (the exponent is factored and written once):
        $2.0$--$4.1{\\times}10^{4}$
    Range, different exponents:
        $3.4{\\times}10^{9}$--$2.6{\\times}10^{10}$

    `digits` is the largest number of decimals a mantissa may carry; trailing zeros are
    dropped (1.25 at digits=2 prints "1.25", 1.1 prints "1.1"). In a factored range a
    whole-number mantissa beside a fractional one keeps one decimal ("2.0"--"4.1"), so the
    two ends read at the same precision. `rel` is a relation put in front of the number
    inside the math group, e.g. rel=r"\\lesssim" gives ${\\lesssim}10^{5}$ for a count the
    chapter states as a bound. The unit (" T") is the caller's.

    `keep_zeros=True` prints every mantissa with exactly `digits` decimals, trailing zeros kept,
    so a box that prints "1.0e9" or "3.0e7" is spelled $1.0{\\times}10^{9}$, not $10^{9}$ (E20 table
    pass, 2026-10-01: each Table 1.1 cell carries its box's significant figures).

    tests/test_table.py compares every Table 1.1 cell against this function's output and
    nothing else. `tex_range` and `tex_range_factored` below are the older spellings, kept
    unchanged for their existing callers (chapter-box strings, test_common).
    """
    if hi is None:
        hi = lo
    if lo <= 0 or hi <= 0:
        raise ValueError("hard-operation counts are positive")
    if hi < lo:
        raise ValueError(f"range is reversed: {lo} > {hi}")
    pre = f"{{{rel}}}" if rel else ""
    (ml, el), (mh, eh) = _mantissa(lo, digits), _mantissa(hi, digits)

    def txt(m: float) -> str:
        return f"{m:.{digits}f}" if keep_zeros else f"{m:g}"

    def one(m: float, e: int) -> str:
        return f"10^{{{e}}}" if (m == 1 and not keep_zeros) else f"{txt(m)}{{\\times}}10^{{{e}}}"

    if (ml, el) == (mh, eh):
        return f"${pre}{one(ml, el)}$"
    if el == eh:
        sl, sh = txt(ml), txt(mh)
        if "." in sl + sh:
            sl, sh = (s if "." in s else s + ".0" for s in (sl, sh))
        return f"${pre}{sl}$--${sh}{{\\times}}10^{{{el}}}$"
    return f"${pre}{one(ml, el)}$--${one(mh, eh)}$"


def tex_range(lo: float, hi: float, digits: int = 1) -> str:
    """LaTeX for a T-count range with the exponent repeated at both ends.

    Same exponent: $1.1{\\times}10^{8}$--$4.8{\\times}10^{8}$
    Superseded for Table 1.1 by `tex_hard_ops` (R10-e); kept for existing callers.
    Use `tex_range_factored` for boxes.
    """
    ml, el = sci(lo, digits)
    mh, eh = sci(hi, digits)
    if lo == hi:
        return f"${ml}{{\\times}}10^{{{el}}}$"
    return f"${ml}{{\\times}}10^{{{el}}}$--${mh}{{\\times}}10^{{{eh}}}$"


def tex_range_factored(lo: float, hi: float, digits: int = 1) -> str:
    """$1.1$--$4.8\\times 10^{8}$ when exponents match; falls back otherwise."""
    ml, el = sci(lo, digits)
    mh, eh = sci(hi, digits)
    if lo == hi:
        return f"${ml}\\times 10^{{{el}}}$"
    if el == eh:
        return f"${ml}$--${mh}\\times 10^{{{el}}}$"
    return f"${ml}\\times 10^{{{el}}}$--${mh}\\times 10^{{{eh}}}$"


def fmt(x: float) -> str:
    m, e = sci(x)
    return f"{m}e{e}"


def fmt_range(r: tuple[float, float] | None) -> str:
    if r is None:
        return "--"
    return fmt(r[0]) if r[0] == r[1] else f"{fmt(r[0])}-{fmt(r[1])}"
