"""Table 1.3 (tab:results, the top-line results by chapter) vs the chapter models.

Every number in every cell of the table is read from the paper, in reading order, and compared
with its source at the precision the table prints. No value is written down twice: SPEC names
where each number comes from, not what it should be, so the table cannot drift from the models.

Three things fail the run:
  - a printed number that does not round to its source (the table or a model changed);
  - a printed number with no source in SPEC (a number was added to the table by hand);
  - a "fits" cell whose model is above the target.

A source is a model value, or, for the few assumptions the models do not export, a phrase in the
chapter's own text (Stated). The words around the numbers are not checked: a qualifier such as
"at least" or "lower bound" is the chapter's statement, carried by hand.

The top-line bullets above the table are checked the same way at the end of this file.

Runs through sync_paper.sh, which copies this file beside the other tests; it needs the paper.
"""

import re
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Callable

import pytest

from estimates import TEX, load
from estimates.common import SECONDS_PER_YEAR as YEAR

DOC = Path(__file__).resolve().parents[3]
OVERVIEW = DOC / "main-overview.tex"
WINDOW_S = 5 * YEAR
TARGET = {"2028": 1e5, "2033": 1e9}
LABEL = {"nu_nucleus": 2, "mu2e_0nubb": 3, "hybrid_lqcd": 4, "qgp_transport": 5, "collider": 6,
         "curved_space": 7, "baryogenesis": 8, "chiral": 9, "finite_density": 10}
COLUMNS = ("2028", "2033", "first", "campaign", "conditions")
SECONDS = {"d": 86400.0, "yr": YEAR, "mo": YEAR / 12}

# Labels that contain digits but are not results. Removed before numbers are read.
NOT_NUMBERS = (r"\$\\Sigma\(216\\times 3\)\$", r"\$0\\nu\\beta\\beta\$", r"\$3\^3\$", r"\$H_I\$",
               r"\$N_\{CS\}\$", r"\$2\\beta\$--\$2\\pi\\beta\$", r"\b20(?:28|33)\b",
               r"\\(?:ref|cite)\{[^}]*\}")
NUMBER = re.compile(r"(\\sim\s*)?(\d+(?:\.\d+)?)(?:\\times\s*10\^\{(\d+)\})?")
UNIT = re.compile(r"\$?\s*(\\times|\\%|yr\b|d\b|mo\b|LQ\b)")


@lru_cache(maxsize=None)
def result(ch, era):
    module = load(ch)
    return module.model(module.Assumptions(), era)


def pick(value, end):
    return value[end] if isinstance(value, (tuple, list)) else value


def ops(ch, era, end=1):
    return pick(result(ch, era).hard_ops, end)


def lq(ch, era, end=1):
    return pick(result(ch, era).lq, end)


def x(ch, era, key, end=1):
    return pick(result(ch, era).intermediates[key], end)


def over(ch, era, end=1):
    """Per-shot hard operations as a multiple of the era's target."""
    return ops(ch, era, end) / TARGET[era]


def campaign_over_window(ch, end=1):
    return x(ch, "2033", "wall_campaign_s", end) / WINDOW_S


@dataclass(frozen=True)
class Source:
    """A printed number and where it comes from. value is in seconds for a wall time."""
    value: Callable[[], float]


@dataclass(frozen=True)
class AtMost:
    """A printed number that bounds its source from above ("register also fits 250 LQ")."""
    value: Callable[[], float]


@dataclass(frozen=True)
class Stated:
    """A printed number that only the chapter's text states: the phrase must be in that chapter."""
    ch: int
    phrase: str
    number: float


@dataclass(frozen=True)
class Fits:
    """The cell says "fits": no printed number, and the source must be at or below the target."""
    value: Callable[[], float]


def S(fn, *args, **kw):
    return Source(lambda: fn(*args, **kw))


def fits(ch, era):
    return Fits(lambda: over(ch, era))


def walls(ch, key="wall_first_result_s"):
    return [S(x, ch, "2033", key, 0), S(x, ch, "2033", key, 1)]


def wall(ch, key="wall_first_result_s", end=0):
    return [S(x, ch, "2033", key, end)]


SPEC = {
    2: {"2028": [fits(2, "2028")],
        "2033": [fits(2, "2033")],
        "first": wall(2),
        "campaign": wall(2, "wall_campaign_s"),
        "conditions": []},
    3: {"2028": [fits(3, "2028"), S(x, 3, "2028", "d_range", 1)],
        "2033": [fits(3, "2033")],
        "first": walls(3),
        "campaign": walls(3, "wall_campaign_s"),
        "conditions": [AtMost(lambda: lq(3, "2033", 1)),
                       S(x, 3, "2033", "multiple_pf", 0), S(x, 3, "2033", "multiple_jj55", 1),
                       S(lambda: result(3, "2033").intermediates["simplified_0nubb"]["pf_rabi"]["lq"][0]),
                       S(x, 3, "2033", "lq_jj55", 1)]},
    4: {"2028": [fits(4, "2028")],
        "2033": [fits(4, "2033")],
        "first": walls(4),
        "campaign": walls(4, "wall_campaign_s") + [S(campaign_over_window, 4)],
        "conditions": [S(lq, 4, "codesign", 0), S(ops, 4, "codesign", 0)]},
    5: {"2028": [Fits(lambda: x(5, "2028", "t_first_result_shot_2028_over_rfi")), S(over, 5, "2028")],
        "2033": [S(over, 5, "2033", 0), S(over, 5, "2033", 1)],
        "first": walls(5),
        "campaign": walls(5, "wall_campaign_s") + [S(campaign_over_window, 5, 0), S(campaign_over_window, 5, 1)],
        "conditions": [S(lambda: 100 * x(5, "2033", "quench_prep_t", 0) / x(5, "2033", "t_deepest", 0)),
                       S(lambda: 100 * x(5, "2033", "quench_prep_t", 1) / x(5, "2033", "t_deepest", 1)),
                       Stated(5, r"$3$--$10\times$", 3), Stated(5, r"$3$--$10\times$", 10)]},
    6: {"2028": [S(over, 6, "2028")],
        "2033": [S(over, 6, "2033", 0), S(over, 6, "2033", 1)],
        "first": walls(6),
        "campaign": walls(6, "wall_campaign_s") + [S(campaign_over_window, 6),
                    S(lambda: x(6, "2033", "wall_ff_yr", 0) * YEAR), S(lambda: x(6, "2033", "wall_ff_yr", 1) * YEAR),
                    S(x, 6, "2033", "ff_cost_reduction_for_horizon", 0), S(x, 6, "2033", "ff_cost_reduction_for_horizon", 1)],
        "conditions": [S(lambda: 100 * x(6, "2033", "hopping_share_of_step_2033")),
                       S(lambda: 100 * x(10, "2033", "hop_share_of_c_be"))]},
    7: {"2028": [fits(7, "2028")],
        "2033": [S(over, 7, "2033")],
        "first": walls(7),
        "campaign": walls(7, "wall_campaign_s"),
        "conditions": []},
    8: {"2028": [S(over, 8, "2028"), S(x, 8, "2028", "ct_reduction_banked"),
                 S(lambda: x(8, "2028", "t_raw_per_shot") / TARGET["2028"])],
        "2033": [fits(8, "2033")],
        "first": wall(8),
        "campaign": wall(8, "wall_campaign_s"),
        "conditions": []},
    9: {"2028": [S(over, 9, "2028", 0), S(over, 9, "2028", 1)],
        "2033": [S(over, 9, "2033")],
        "first": wall(9),
        "campaign": walls(9, "wall_campaign_s"),
        "conditions": []},
    10: {"2028": [S(over, 10, "2028")],
         "2033": [S(over, 10, "2033", 0), S(over, 10, "2033", 1)],
         "first": walls(10),
         "campaign": walls(10, "wall_campaign_s") + [S(campaign_over_window, 10, 0), S(campaign_over_window, 10, 1)],
         "conditions": []},
}


def table_rows():
    text = OVERVIEW.read_text()
    start = text.index(r"\textbf{Main conditions}")
    body = text[start:text.index(r"\end{tabularx}", start)]
    rows = {}
    for line in body.splitlines():
        m = re.match(r"Ch\.~\\ref\{ch:(\w+)\}\s*&", line)
        if m:
            cells = [c.strip() for c in re.split(r"(?<!\\)&", line)]
            assert len(cells) == 6, f"expected 6 columns: {line[:80]}"
            rows[LABEL[m.group(1)]] = dict(zip(COLUMNS, cells[1:]))
    return rows


def tokens(cell):
    """(text, value, unit, approximate) for each number in a cell, in reading order. A number
    joined by -- to the next takes that number's unit: $6.0$--$25$ d is two values in days."""
    for pattern in NOT_NUMBERS:
        cell = re.sub(pattern, " ", cell)
    out = []
    for m in NUMBER.finditer(cell):
        unit = UNIT.match(cell[m.end():])
        exponent = int(m.group(3)) if m.group(3) else 0
        out.append({"text": m.group(2), "scale": 10.0 ** exponent, "unit": unit.group(1) if unit else None,
                    "approx": bool(m.group(1)), "joined": cell[m.end():].lstrip("$ ").startswith("--")})
    for here, next_one in zip(reversed(out[:-1]), reversed(out[1:])):
        if here["unit"] is None and here["joined"]:
            here["unit"] = next_one["unit"]
    return out


def significant(text):
    digits = text.replace(".", "").lstrip("0")
    return len(digits) if "." in text else max(1, len(digits.rstrip("0")))


def agrees(token, value):
    """Does value round to the printed number at the printed significant figures?"""
    printed = float(token["text"])
    value = value / SECONDS.get(token["unit"], 1.0) / token["scale"]
    if token["approx"]:
        return abs(value / printed - 1) < 0.05
    return float(f"{value:.{significant(token['text'])}g}") == printed


def chapter_text(ch):
    return (DOC / TEX[ch]).read_text()


CELLS = [(ch, col) for ch in SPEC for col in COLUMNS]


def test_every_chapter_has_a_row():
    assert set(table_rows()) == set(SPEC)


@pytest.mark.parametrize("ch,col", CELLS)
def test_cell_matches_its_sources(ch, col):
    cell = table_rows()[ch][col]
    spec = SPEC[ch][col]
    found = tokens(cell)
    printed = iter(found)
    n_printed = sum(not isinstance(s, Fits) for s in spec)
    assert len(found) == n_printed, (
        f"Ch. {ch} {col}: {len(found)} numbers printed {[t['text'] for t in found]}, "
        f"{n_printed} have a source. Every number in Table 1.3 needs one in SPEC.")
    for source in spec:
        if isinstance(source, Fits):
            assert "fits" in cell, f"Ch. {ch} {col}: expected a 'fits' cell: {cell}"
            assert float(f"{source.value():.2g}") <= 1.0, f"Ch. {ch} {col} says fits, source {source.value():.3g}x"
            continue
        token = next(printed)
        if isinstance(source, Stated):
            assert source.phrase in chapter_text(source.ch), f"Ch. {source.ch} no longer says {source.phrase}"
            assert float(token["text"]) == source.number, f"Ch. {ch} {col}: {token['text']} vs {source.phrase}"
        elif isinstance(source, AtMost):
            assert source.value() <= float(token["text"]), f"Ch. {ch} {col}: {source.value()} exceeds {token['text']}"
        else:
            assert agrees(token, source.value()), (
                f"Ch. {ch} {col}: prints {token['text']} {token['unit'] or ''}, source gives "
                f"{source.value() / SECONDS.get(token['unit'], 1.0) / token['scale']:.4g}")


# The bullets above the table: what they count and range over comes from the same sources.

NUMBER_WORD = {4: "four", 5: "Five"}
DIGITAL = sorted(SPEC)


def top_line():
    return re.sub(r"\s+", " ", OVERVIEW.read_text())


def test_registers_bullet():
    registers = [lq(ch, "2028") for ch in DIGITAL]
    assert f"All nine 2028 benchmarks fit the 250-LQ register, at {min(registers):.0f}--{max(registers):.0f} LQ." in top_line()
    assert max(registers) <= 250
    large = [lq(ch, "2033") for ch in (5, 6, 7, 9)]
    assert f"reach ${min(large):.0f}$--${max(large):.0f}$ LQ in Chs." in top_line()


@pytest.mark.parametrize("era,target,lo_end", [("2028", r"$10^{5}$", 1), ("2033", r"$10^{9}$", 0)])
def test_hard_operations_bullet(era, target, lo_end):
    above = [ch for ch in DIGITAL if float(f"{over(ch, era):.2g}") > 1.0]
    assert len(above) == 5
    lo = min(over(ch, era, lo_end) for ch in above)
    hi = max(over(ch, era) for ch in above)
    assert (f"Five of the nine {era} " in top_line()
            and f"exceed {target}, by factors of {lo:.2g} to {hi:.2g}" in top_line())


def test_co_design_band_in_bullet():
    assert [ch for ch in DIGITAL if over(ch, "2033") >= 10] == [5, 6, 7, 9]
    assert over(6, "2033", 0) < 10 <= over(6, "2033", 1)      # "Ch. 6 at its upper end"


def test_machine_time_bullet():
    within = [ch for ch in DIGITAL if campaign_over_window(ch) <= 1]
    assert within == [2, 3, 7, 8, 9]
    text = top_line()
    assert "Five of the nine campaigns fit the five-year window on one machine." in text
    assert f"{float(f'{campaign_over_window(5):.2g}'):.0f}\\times$ (Ch." in text
    worst = x(6, "2033", "ff_cost_reduction_for_horizon", 1)
    assert f"${float(f'{worst:.2g}'):.0f}\\times$." in text


def test_far_end_bullet():
    text = top_line()
    r3, r10 = result(3, "codesign"), result(10, "codesign")
    assert f"{r3.lq[0]:.0f}--{r3.lq[1]:.0f} LQ" in text
    for value in r3.hard_ops:
        mantissa, exponent = f"{value:.1e}".split("e+")
        assert f"{mantissa}\\times 10^{{{int(exponent)}}}" in text
    assert f"{r10.lq[0]:.0f}--{r10.lq[1]:.0f}-LQ register without an estimated depth" in text


def test_wall_time_sentence():
    text = top_line()
    assert f"from ${x(9, '2028', 'wall_first_result_s'):.0f}$~s (the first result of Ch." in text
    assert f"to ${x(4, '2028', 'wall_total_corrected_h'):.1f}$~h (Ch." in text
