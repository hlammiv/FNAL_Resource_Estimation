"""The results in the chapter table (Table 1.1) and the top-line bullets vs the chapter models.

test_table.py checks each row's register and hard-operation count. This file checks what the
table adds to them: the 2033 first-result and campaign wall times, and the status symbols.

The wall times are read from the paper in reading order and compared with the models at the
printed precision, and nothing is written down twice. Each symbol must be the one its model
gives: the upper end of the range against the target, the worse of register and hard
operations, or for the campaign against the five-year window. Fits means at or below the target
at two significant figures, near means below 3 times it, and anything else is not yet fitting.

Three things fail the run: a wall time that does not round to its model, a printed number with
no source in SPEC, and a symbol that disagrees with its model. Qualifiers such as a lower-bound
mark are the chapter's statement and are not checked.

Runs through sync_paper.sh, which copies this file beside the other tests; it needs the paper.
"""

import re
from functools import lru_cache
from pathlib import Path

import pytest

from estimates import load
from estimates.common import SECONDS_PER_YEAR as YEAR

DOC = Path(__file__).resolve().parents[3]
OVERVIEW = DOC / "main-overview.tex"
WINDOW_S = 5 * YEAR
TARGET = {"2028": 1e5, "2033": 1e9}
LABEL = {"nu_nucleus": 2, "mu2e_0nubb": 3, "hybrid_lqcd": 4, "qgp_transport": 5, "collider": 6,
         "curved_space": 7, "baryogenesis": 8, "chiral": 9, "finite_density": 10}
LQ_TARGET = {"2028": 250, "2033": 1000}
NEAR = 3.0
SECONDS = {"d": 86400.0, "yr": YEAR, "mo": YEAR / 12}

# Not results: removed before numbers are read.
NOT_NUMBERS = (r"\\fits(?:yes|near|no)",)
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


def status(ratio):
    """The symbol for a value as a multiple of its target."""
    rounded = float(f"{ratio:.2g}")
    return r"\fitsyes" if rounded <= 1.0 else r"\fitsnear" if rounded < NEAR else r"\fitsno"


def machine_status(ch, era):
    """The worse of register and hard operations, at the upper end of each range."""
    return status(max(over(ch, era), lq(ch, era) / LQ_TARGET[era]))


def campaign_status(ch):
    return status(campaign_over_window(ch))


def wall_sources(ch, key):
    value = result(ch, "2033").intermediates[key]
    return [lambda: value[0], lambda: value[1]] if isinstance(value, (tuple, list)) else [lambda: value]


# Ranges printed as one number because the model gives one value at the printed precision.
ONE_VALUE = {(2, "wall_first_result_s"), (2, "wall_campaign_s"), (8, "wall_first_result_s"),
             (8, "wall_campaign_s"), (9, "wall_first_result_s")}
DIGITAL = sorted(LABEL.values())


def sources(ch, key):
    found = wall_sources(ch, key)
    return found[:1] if (ch, key) in ONE_VALUE else found


def table_rows():
    rows = {}
    for line in OVERVIEW.read_text().splitlines():
        m = re.match(r"\\ref\{ch:(\w+)\}\s*&", line)
        if m and m.group(1) in LABEL:
            cells = [c.strip() for c in re.split(r"(?<!\\)&", line)]
            assert len(cells) == 6, f"expected 6 columns: {line[:80]}"
            rows[LABEL[m.group(1)]] = dict(zip(("2028", "2033", "first", "campaign"), cells[2:]))
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


def test_every_digital_chapter_has_a_row():
    assert set(table_rows()) == set(DIGITAL)


@pytest.mark.parametrize("ch", DIGITAL)
@pytest.mark.parametrize("col,key", [("first", "wall_first_result_s"), ("campaign", "wall_campaign_s")])
def test_wall_times_match_the_models(ch, col, key):
    found = tokens(table_rows()[ch][col])
    wanted = sources(ch, key)
    assert len(found) == len(wanted), (
        f"Ch. {ch} {col}: {len(found)} numbers printed {[t['text'] for t in found]}, {len(wanted)} have a source")
    for token, value in zip(found, wanted):
        assert agrees(token, value()), (
            f"Ch. {ch} {col}: prints {token['text']} {token['unit'] or ''}, model gives "
            f"{value() / SECONDS.get(token['unit'], 1.0):.4g}")


@pytest.mark.parametrize("ch", DIGITAL)
@pytest.mark.parametrize("col", ["2028", "2033", "campaign"])
def test_status_symbols_match_the_models(ch, col):
    cell = table_rows()[ch][col]
    printed = re.findall(r"\\fits(?:yes|near|no)\b", cell)
    assert len(printed) == 1, f"Ch. {ch} {col}: expected one symbol: {cell}"
    wanted = campaign_status(ch) if col == "campaign" else machine_status(ch, col)
    assert printed[0] == wanted, f"Ch. {ch} {col}: prints {printed[0]}, model gives {wanted}"


# The bullets above the table: what they count and range over comes from the same sources.

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


def test_far_end_bullet():
    text = top_line()
    r3, r10 = result(3, "codesign"), result(10, "codesign")
    r2 = result(2, "codesign")
    assert f"$^{{12}}$C, {r2.lq[0]:.0f}--{r2.lq[1]:.0f} LQ, ${r2.hard_ops[0] / 1e12:.1f}$--${r2.hard_ops[1] / 1e12:.1f}\\times 10^{{12}}$" in text
    assert f"{r3.lq[0]:.0f}--{r3.lq[1]:.0f} LQ" in text
    for value in r3.hard_ops:
        mantissa, exponent = f"{value:.1e}".split("e+")
        assert f"{mantissa}\\times 10^{{{int(exponent)}}}" in text
    assert f"{r10.lq[0]:.0f}--{r10.lq[1]:.0f}-LQ register without an estimated depth" in text


def test_wall_time_sentence():
    text = top_line()
    assert f"from ${x(9, '2028', 'wall_first_result_s'):.0f}$~s (the first result of Ch." in text
    assert f"to ${x(4, '2028', 'wall_total_corrected_h'):.1f}$~h (Ch." in text
