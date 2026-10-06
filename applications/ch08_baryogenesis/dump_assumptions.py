"""Print the Ch. 8 assumptions table (Markdown) from the live model.

    python applications/ch08_baryogenesis/dump_assumptions.py

Every row is read from `estimates.ch08_baryogenesis.Assumptions()`: field name, value, provenance tag,
source and note, exactly as the code carries them. Rows are grouped by the section comments
(`# ---- ... ----`) of the dataclass, in source order. Nothing is computed or rounded here beyond
printing floats to ten significant figures.
"""

from __future__ import annotations

import inspect
import re
import sys
from collections import Counter
from dataclasses import fields
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from estimates import ch08_baryogenesis as m          # noqa: E402
from estimates.common import Provenance, Tagged       # noqa: E402

TAG = {Provenance.CITED: "Cited", Provenance.ASSUMED: "Assumed",
       Provenance.STATED: "Stated-not-derived", Provenance.UNCITED: "Uncited"}

# Reader-facing titles for the dataclass section comments (matched on the comment's opening words).
TITLES = (
    ("lattice and registers", "Lattice and registers"),
    ("c_T derivation inputs", "Inputs to the gate-level c_T"),
    ("Trotter", "Trotter schedule"),
    ("shots and readout", "Shots, readout and fault budget"),
    ("wall time and campaign", "Wall time and campaign"),
    ("two-tier schedule", "Two-tier schedule (first result and campaign) and T-depth"),
    ("state preparation", "State preparation, flag contrast, bath and fault model"),
    ("interim variant", "Interim scalar+Wilson variant (printed values)"),
    ("utility box", "Utility box"),
)
SECTION = re.compile(r"^\s*# -{3,}\s*(.*?)\s*-{3,}\s*$")
FIELD = re.compile(r"^\s{4}(\w+)\s*:\s*Tagged\b")


def _title(raw: str) -> str:
    for key, title in TITLES:
        if raw.startswith(key):
            return title
    return raw


def sections() -> dict[str, str]:
    """field name -> section title, from the comments in the class source."""
    out, cur = {}, "Other"
    for line in inspect.getsource(m.Assumptions).splitlines():
        s = SECTION.match(line)
        if s:
            cur = _title(s.group(1))
            continue
        f = FIELD.match(line)
        if f:
            out[f.group(1)] = cur
    return out


def fmt(x) -> str:
    if isinstance(x, bool) or isinstance(x, str):
        return str(x)
    if isinstance(x, int):
        return str(x)
    if isinstance(x, float):
        return f"{x:.10g}"
    return str(x)


def value(v: Tagged) -> str:
    if v.is_range:
        return "(" + ", ".join(fmt(x) for x in v.value) + ")"
    return fmt(v.value)


def cell(s: str) -> str:
    return s.replace("|", "\\|").replace("\n", " ")


def main() -> None:
    a = m.Assumptions()
    sec = sections()
    rows = [(f.name, getattr(a, f.name)) for f in fields(a)]
    missing = [n for n, _ in rows if n not in sec]
    if missing:
        raise SystemExit(f"fields not found in the class source: {missing}")
    counts = Counter(TAG[v.prov] for _, v in rows)
    print(f"{len(rows)} fields: " + ", ".join(f"{k} {counts[k]}" for k in TAG.values() if counts[k]))
    current = None
    for name, v in rows:
        if sec[name] != current:
            current = sec[name]
            print(f"\n**{current}**\n")
            print("| field | value | provenance | source | note |")
            print("|---|---|---|---|---|")
        print(f"| `{name}` | {cell(value(v))} | {TAG[v.prov]} | {cell(v.src) or '-'} | {cell(v.note) or '-'} |")


if __name__ == "__main__":
    main()
