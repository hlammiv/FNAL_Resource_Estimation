"""Table 1.1 vs the chapter boxes and the chapter models: every cell of both generation
columns, Chs. 2-10, prints its chapter box's headline register and per-shot hard-operation
count, with the box's own significant figures, in both copies of the table (rulings R10-a..e,
H. Lamm, 2026-09-29; author request 2026-10-01 "for table 1.1 get the figures and numbers to
match": a cell carries the box's value at the box's precision, e.g. 1.0e9 not 10^9, 1082 not 1.1k).

The two copies are doc/main-overview.tex and doc/table_standalone.tex; rows are found by
the chapter label phrases in scripts/verify_table.py ROW_KEY, and the cells are the third
(first generation, "2028") and fourth (next generation, "2033") `&`-separated fields.

What is compared. A cell reads "<instance>; <register> LQ, <count> T ...". The test takes
the FIRST "... LQ" token and the FIRST "$...$ T" token of the cell and requires
  (1) cell == box spelling. The box value is the chapter model's PUBLISHED[era] (the box as
      printed, kept in step with the box by each chapter's own test file), or the HEADLINE
      override when the box headlines one value out of a published band (Ch. 5 2033 '~1120').
      It is formatted by
        hard ops : common.tex_hard_ops at the box's digits (t_digits, keep_zeros)
        register : lq_text below (integers, a range written "lo--hi"; no k-notation)
  (2) box == model. The model's hard-op count (Result.hard_ops, or the intermediate named by
      HEADLINE.t) formatted at the box's precision must give the same spelling, and the box
      register must lie inside the model's register band within PUBLISHED.rel_tol.
Strict xfails (XFAIL) are the filed disagreements between a box headline and its model;
each names its NEEDS_AUTHOR.md item, and a fix flips it to XPASS and fails the run.
Every other "$...$ T" token in a cell is
checked for the one spelling by a parse-and-reformat round trip.
"""

import re
from dataclasses import dataclass
from pathlib import Path

import pytest

import estimates
from estimates import common as c
from verify_table import ROW_KEY, fields, norm, rows

DOC = Path(__file__).resolve().parents[3]
TABLES = {p: DOC / p for p in ("main-overview.tex", "table_standalone.tex")}
ERAS = ("2028", "2033")

# rows whose cell deliberately carries no T-count: none since R10-a (2026-09-29).
# Was {6, 7, 10} for the next-generation column, under a caption clause now deleted.
NO_T: set = set()


@dataclass(frozen=True)
class Headline:
    lq: str | None = None              # model intermediate behind the box register (None: Result.lq)
    t: str | None = None               # model intermediate behind the box count (None: Result.hard_ops)
    t_digits: int = 1                  # decimals the box prints on the mantissa
    keep_zeros: bool = False           # the box keeps a trailing zero ("1.0e9", "3.0e7")
    box_lq: tuple | None = None        # box headline register when it is not PUBLISHED.lq
    box_t: tuple | None = None         # box headline count when it is not PUBLISHED.hard_ops
    t_bound: bool = False              # the box states the count as an upper bound (\lesssim)
    t_rel: str = ""                    # relation the box prints before a one-sided count, e.g. r"\gtrsim" (lower bound, referee G2)


H = Headline
HEADLINE = {
    # ---- first generation -------------------------------------------------------------
    (2, "2028"): H(t_digits=2),                      # box '36', '3.89e4 T' (ruling 2026-10-12: report synthesis convention; study 8.03e4)
    (3, "2028"): H(keep_zeros=True),                 # box '~120-170', '1.7e4-1.0e5 T' (insertion + complete load + load-angle rotations, D = 1e2-600; M3 ruling 2026-10-04 folded the angles in, conditional on D; G2 1.6e4-9.9e4; was the 5.4e3 insertion)
    (4, "2028"): H(),                                # box '18', '1.9-3.5e4' (1.898e4-3.544e4)
    (5, "2028"): H(lq="lq_2028_hi"),                 # box '178', first result '9.5e4' / full benchmark '5.9e5' (r25: 16 workspace LQ, R10; each timeslice at its own depth, R4; was 162, '5.9e5')
    (6, "2028"): H(lq="lq_2028"),                    # box '230' (64 + 144 + 22; r25 R10, was 216 = 64 + 144 + 8), '1.3e5' (E20, 1 ramp + 2 evo: 132,174.9)
    (7, "2028"): H(),                                # box '220-240' (222-242), '9.5e4' (94,753; option (c): one step at K = 8, smaller ruling (b) 2026-10-05; was '160-180', '9.4e4')
    (8, "2028"): H(t_rel=r"\gtrsim"),                # box '174', '>~2.8e5' (2.756e5 at 20 steps, ruling 2026-10-11; preparation unpriced, assumed 3.3x cut; referee B3 2026-10-05)
    (9, "2028"): H(),                                # box '250' (r25 R10: 42 ancilla, was 240), '2.7-2.9e5'; the cell prints the band (268,995-288,535)
    (10, "2028"): H(keep_zeros=True, t_rel=r"\gtrsim"),  # box '~180' (176), '>~4.0e5' (r19 E20+E21: 399,956.8; one Trotter step, a lower bound on the ramp, referee G2)
    # ---- next generation --------------------------------------------------------------
    (2, "2033"): H(t_digits=2),                                # box '908', '2.07e8' (ruling 2026-10-12; study 7.23e8)
    (3, "2033"): H(),                                # box '74-174', '1.1e8-4.4e8' (1.054e8-4.404e8)
    (4, "2033"): H(),                                # box '83', '2.1e8' (2.12992e8; r25 R3: M at D V log2 V; was '~85', '~9e8')
    (5, "2033"): H(lq="lq_2033_mid", box_lq=(1280, 1280), keep_zeros=True),   # box '~1280 (1234-1334)'; deepest shot 6.3e9-2.2e10 (quench route, T = 225-300 MeV x c = 2-2pi, temperature ruling 2026-10-05; 3.0e10 on the second spacing's late shot; was 1.1e10-3.1e10 at 160 MeV); s648 2026-10-05: Sigma(216x3), 11 qubits/link, H_I, a_s = 0.2 fm, second spacing 0.15 fm own step, 48 steps on its slow late shot (was 1120, 1.3-3.2e8 at rulings2 g^2 = 3.5; 2.2-9.6e8 with Sigma(72x3) at 0.2 fm)
    (6, "2033"): H(),                                # box '964', '4.9e9-2.6e10' (step rule 2026-10-04: 4.65 GeV point 575 steps; 4.9014e9-2.6159e10; r26 J1 2.5381e10; r25 3.9e9-1.5e10); electric term a lower bound
    (7, "2033"): H(keep_zeros=True),                 # box '1050' (prose 1030-1050), '1.5e10' (1.486e10; referee U2 step m_phi dt = 0.05, 2026-10-04; was '1.0e9')
    (8, "2033"): H(t_rel=r"\gtrsim"),                # box '1000', deepest MLAE circuit '>~9.9e8' (33 x 3.0056e7 + 16 x 7.42e4 = 9.930e8; free preparation and reflections priced 2026-10-05; G2 2026-10-04, was the 3.0e7 base circuit; lower bound, dressing ramp unpriced)
    (9, "2033"): H(keep_zeros=True),                 # box '1015', N_CS method test '2.5e10' (Yukawa weight in the normalization 2026-10-10: 2.4534e10; before 1.9e10 / 1.8909e10; r25 '1.6e10'); the condensate (2.8e9) is not the cell
    (10, "2033"): H(t="t_per_shot_range"),           # box '~500-1000' (514-1014), '2.8-5.6e9' (author ruling 2026-10-05: QSVT/TPQ thermal state, 17 filter calls of degree 30, 2-4 c_step per query, 2.7726e9-5.6340e9, conditional on the reference; was the Gibbs sampler 8.4e13-4.4e14, open item 1(a) 2026-10-04; before that '>~5.5e9')
}

# second instances a cell prints after the headline one: (chapter, column) -> model era

# filed disagreements, (chapter, column, quantity) -> reason naming the NEEDS_AUTHOR.md item.
# "Tmodel" = the box's count against the model at the box's precision.
# none open: the E20 Ch. 10 hold closed in r19 (E21, 2026-10-01).
XFAIL: dict = {}

CELLS = {name: {era: rows(path, era) for era in ERAS} for name, path in TABLES.items()}
FIELDS = {name: fields(path) for name, path in TABLES.items()}

LQ_TOKEN = re.compile(r"(?<![\w.$\-])([0-9][0-9.]*(?:--[0-9][0-9.]*)?k?) LQ(?![\w-])")
T_TOKEN = re.compile(r"(?<![\w.\-])(\$[^\s;,]*\$) T(?![\w-])")


def _round_sig(v: float, sig: int | None) -> float:
    if sig is None:
        return float(v)
    m, e = c.sci(v, sig - 1)
    return float(m) * 10 ** e


def lq_text(lo: float, hi: float, sig: int | None = None) -> str:
    """The table's register spelling: the box's integers as printed, a range written once with
    "--" ("800--1000"), and one value when both ends round alike. No k-notation (2026-10-01:
    a cell carries the box's significant figures, so 1082 stays 1082, not 1.1k)."""
    a, b = _round_sig(lo, sig), _round_sig(hi, sig)
    sa, sb = f"{a:g}", f"{b:g}"
    return sa if sa == sb else f"{sa}--{sb}"


def _pair(v) -> tuple[float, float]:
    return (float(v[0]), float(v[1])) if isinstance(v, (tuple, list)) else (float(v), float(v))


def _model(ch: int, era: str):
    m = estimates.load(ch)
    return m, m.model(m.Assumptions(), era)


def _fmt_t(t: tuple[float, float], h: Headline) -> str:
    if h.t_bound:
        return c.tex_hard_ops(t[1], digits=h.t_digits, rel=r"\lesssim", keep_zeros=h.keep_zeros)
    if h.t_rel:
        return c.tex_hard_ops(t[0], digits=h.t_digits, rel=h.t_rel, keep_zeros=h.keep_zeros)
    return c.tex_hard_ops(*t, digits=h.t_digits, keep_zeros=h.keep_zeros)


def box(ch: int, era: str, h: Headline) -> tuple[str, str]:
    """(register, hard ops) as the chapter box headlines them, in the table's spelling."""
    m = estimates.load(ch)
    pub = m.PUBLISHED[era]
    lq = _pair(h.box_lq) if h.box_lq else _pair(pub.lq)
    t = _pair(h.box_t) if h.box_t else _pair(pub.hard_ops)
    return lq_text(*lq), _fmt_t(t, h)


def expected(ch: int, era: str, h: Headline) -> tuple[tuple[float, float], str]:
    """(register band, hard ops at the box's precision) from the chapter model."""
    m, r = _model(ch, era)
    lq = _pair(r.intermediates[h.lq]) if h.lq else _pair(r.lq)
    t = _pair(r.intermediates[h.t]) if h.t else _pair(r.hard_ops)
    tol = m.PUBLISHED[era].rel_tol
    # a headline taken from an intermediate must sit inside the band the model returns
    assert r.lq[0] * (1 - tol) <= lq[0] <= lq[1] <= r.lq[1] * (1 + tol), (ch, era, h.lq, lq, r.lq)
    assert r.hard_ops[0] * (1 - tol) <= t[0] <= t[1] <= r.hard_ops[1] * (1 + tol), \
        (ch, era, h.t, t, r.hard_ops)
    return (lq, tol), _fmt_t(t, h)


def _cell(ch: int, era: str, table: str) -> str:
    cell = CELLS[table][era].get(ch)
    assert cell is not None, f"Ch{ch} ({ROW_KEY[ch]!r}) has no row in {table}"
    return norm(cell)


def _tokens(cell: str) -> tuple[list[str], list[str]]:
    return LQ_TOKEN.findall(cell), T_TOKEN.findall(cell)


def _xfail(request, key):
    if key in XFAIL:
        request.applymarker(pytest.mark.xfail(strict=True, reason=XFAIL[key]))


PARAMS = [(ch, era) for era in ERAS for ch in sorted(ROW_KEY)]


@pytest.mark.parametrize("table", sorted(TABLES))
@pytest.mark.parametrize("ch,era", PARAMS)
def test_lq_in_table(ch, era, table, request):
    """The cell prints the box's register, as the box prints it."""
    _xfail(request, (ch, era, "LQ"))
    want, _ = box(ch, era, HEADLINE[ch, era])
    lqs, _ = _tokens(_cell(ch, era, table))
    assert lqs and lqs[0] == want, (ch, era, table, want, lqs)


@pytest.mark.parametrize("table", sorted(TABLES))
@pytest.mark.parametrize("ch,era", PARAMS)
def test_t_in_table(ch, era, table, request):
    """The cell prints the box's per-shot count, at the box's significant figures."""
    if ch in NO_T:
        pytest.skip(f"Ch{ch} row carries no T-count (NO_T)")
    _xfail(request, (ch, era, "T"))
    _, want = box(ch, era, HEADLINE[ch, era])
    _, ts = _tokens(_cell(ch, era, table))
    assert ts and ts[0] == want, (ch, era, table, want, ts)


@pytest.mark.parametrize("ch,era", PARAMS)
def test_box_register_is_the_model(ch, era, request):
    _xfail(request, (ch, era, "LQmodel"))
    h = HEADLINE[ch, era]
    (lq, tol), _ = expected(ch, era, h)
    m = estimates.load(ch)
    blq = _pair(h.box_lq) if h.box_lq else _pair(m.PUBLISHED[era].lq)
    assert lq[0] * (1 - tol) <= blq[0] <= blq[1] <= lq[1] * (1 + tol), (ch, era, blq, lq, tol)


@pytest.mark.parametrize("ch,era", PARAMS)
def test_box_count_is_the_model_at_box_precision(ch, era, request):
    _xfail(request, (ch, era, "Tmodel"))
    h = HEADLINE[ch, era]
    _, model_t = expected(ch, era, h)
    _, box_t = box(ch, era, h)
    assert model_t == box_t, (ch, era, "model", model_t, "box", box_t)


@pytest.mark.parametrize("table", sorted(TABLES))
def _parse(token: str) -> tuple[float, float, int, str]:
    """Numbers, mantissa decimals and relation of a '$...$' hard-op token, whatever its spelling."""
    rel = re.match(r"\$\{?(\\[a-z]+)\}?", token)
    body = token.replace(rel.group(1), "") if rel else token
    body = body.replace("{", "").replace("}", "").replace("$", "").replace(" ", "")
    ends = body.split("--")
    exps = [re.search(r"10\^(-?\d+)", e) for e in ends]
    assert exps[-1], token
    out, digits, zeros = [], 0, False
    for e, x in zip(ends, exps):
        x = x or exps[-1]                                  # factored: the exponent is written once
        mant = re.sub(r"(\\times)?10\^-?\d+", "", e) or "1"
        digits = max(digits, len(mant.partition(".")[2]))
        zeros = zeros or ("." in mant and mant.endswith("0"))
        out.append(float(mant) * 10 ** int(x.group(1)))
    return out[0], out[-1], digits, (rel.group(1) if rel else ""), zeros


@pytest.mark.parametrize("table", sorted(TABLES))
@pytest.mark.parametrize("ch,era", PARAMS)
def test_one_spelling(ch, era, table):
    """R10-e: every hard-op count in the cell is spelled as common.tex_hard_ops spells it."""
    _, ts = _tokens(_cell(ch, era, table))
    assert ts, (ch, era, table, "no T-count in the cell")
    for tok in ts:
        lo, hi, digits, rel, zeros = _parse(tok)
        assert tok == c.tex_hard_ops(lo, hi, digits=digits, rel=rel, keep_zeros=zeros), (ch, era, table, tok)


@pytest.mark.parametrize("ch", sorted(ROW_KEY))
def test_table_copies_agree(ch):
    """Application, both generation cells and the program column are identical in the two copies."""
    a, b = (FIELDS[t].get(ch) for t in sorted(TABLES))
    assert a is not None and b is not None, ch
    assert [norm(x) for x in a[1:5]] == [norm(x) for x in b[1:5]], (ch, a, b)


def test_every_cell_is_specified():
    assert set(HEADLINE) == {(ch, era) for ch in ROW_KEY for era in ERAS}
    assert not NO_T
    assert {k[:2] for k in XFAIL} <= set(HEADLINE)


def test_formatter_is_the_table_convention():
    f = c.tex_hard_ops
    assert f(2.0e4, 4.1e4) == r"$2.0$--$4.1{\times}10^{4}$"
    assert f(3.4e9, 2.6e10) == r"$3.4{\times}10^{9}$--$2.6{\times}10^{10}$"
    assert f(5e4) == r"$5{\times}10^{4}$"
    assert f(1e9) == r"$10^{9}$"
    assert f(8.52e8, digits=0) == r"$9{\times}10^{8}$"
    assert f(1.1e5, 1.25e5, digits=2) == r"$1.1$--$1.25{\times}10^{5}$"
    assert f(1e5, rel=r"\lesssim") == r"${\lesssim}10^{5}$"
    assert f(9.763e8, rel=r"\gtrsim") == r"${\gtrsim}9.8{\times}10^{8}$"
    assert f(6.7e8, 9.6e8) == f(6.7e8, 9.6e8, 1) == r"$6.7$--$9.6{\times}10^{8}$"
    assert f(1.008e9, keep_zeros=True) == r"$1.0{\times}10^{9}$"
    assert f(2.958e7, keep_zeros=True) == r"$3.0{\times}10^{7}$"
    assert f(1.00428e5, digits=3) == r"$1.004{\times}10^{5}$"
    assert f(1.6e10, 3.0e10, keep_zeros=True) == f(1.6e10, 3.0e10) == r"$1.6$--$3.0{\times}10^{10}$"
    assert lq_text(1082, 1082) == "1082"
    assert lq_text(500, 1000) == "500--1000"
    assert lq_text(74, 174) == "74--174"
    assert lq_text(171, 172, 2) == "170"
    # the tokenizer does not find "1k LQ" inside "1.1k LQ", nor "170 LQ" inside "120--170 LQ"
    assert LQ_TOKEN.findall("x; 1.1k LQ, y") == ["1.1k"]
    assert LQ_TOKEN.findall("x; 120--170 LQ, 85 LQ register in a $\\sim$1k-LQ machine") == ["120--170", "85"]
    assert T_TOKEN.findall(r"a $3^3$; 1k LQ, $1.2$--$2.2{\times}10^{9}$ T; pairs: T ${\times}1.7$--$9.9$ over") \
        == [r"$1.2$--$2.2{\times}10^{9}$"]
