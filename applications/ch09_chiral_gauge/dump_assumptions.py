"""Print the Ch. 9 assumption table (Markdown) from the live model.

    python applications/ch09_chiral_gauge/dump_assumptions.py          # from the repository root

Every row is read from estimates.ch09_chiral_gauge.Assumptions(): field name, value, provenance tag,
source and note. Nothing is typed here. Notes are quoted verbatim, whitespace collapsed and cut at
NOTE_CHARS characters; the full text is in the code.

Rows are grouped by instance and, inside each, by role:
  inputs              fields the code does not mark as a printed figure or a comparison record
  printed figures     fields named *_stated: numbers the chapter prints, which the tests compare with
                      the model's own value
  comparison records  fields whose source or note the code marks RECORD or RETIRED: inputs of the
                      alternative constructions and sensitivities the model also evaluates
"""
from __future__ import annotations

import sys
from dataclasses import fields
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from estimates import ch09_chiral_gauge as m  # noqa: E402

NOTE_CHARS = 150
PROV = {"cited": "Cited", "assumed": "Assumed", "stated-not-derived": "Stated-not-derived", "uncited": "Uncited"}
# The dataclass is laid out shared / 2028 / 2033; these names open the 2028 and 2033 blocks.
FIRST_2028, FIRST_2033 = "gauge_group_2028", "gauge_group_2033"
# Fields in the 2033 block that both eras read (the condensate variance, the target precision, the
# first-result precision); fields with '2028' or '2033' in the name belong to that instance wherever they sit.
BOTH_ERAS = {"sigma2_condensate", "eps_obs", "eps_first"}


def _cell(s: str) -> str:
    s = " ".join(str(s).split())
    return s.replace("|", "\\|")


def _num(x) -> str:
    """A float in short form (1.5e+10, not 15000000000.0); checked to round-trip exactly."""
    if isinstance(x, float) and x != 0 and (abs(x) >= 1e5 or abs(x) < 1e-3):
        mant, exp = f"{x:.6e}".split("e")
        s = f"{mant.rstrip('0').rstrip('.')}e{exp}"
        return s if float(s) == x else repr(x)
    return str(x)


def _value(v) -> str:
    if isinstance(v, tuple):
        return "(" + ", ".join(_num(x) for x in v) + ")"
    return _num(v)


def _note(s: str) -> str:
    s = " ".join(str(s).split())
    return _cell(s if len(s) <= NOTE_CHARS else s[:NOTE_CHARS].rstrip() + " ...")


def _role(name: str, v) -> str:
    text = f"{v.src} {v.note}"
    if "RECORD" in text or "RETIRED" in text:
        return "comparison records"
    if name.endswith("_stated"):
        return "printed figures"
    return "inputs"


def rows():
    a = m.Assumptions()
    block = "shared"
    out = []
    for f in fields(a):
        if f.name == FIRST_2028:
            block = "2028"
        elif f.name == FIRST_2033:
            block = "2033"
        v = getattr(a, f.name)
        inst = block
        if "2028" in f.name:
            inst = "2028"
        elif "2033" in f.name:
            inst = "2033"
        elif f.name in BOTH_ERAS:
            inst = "shared"
        out.append((inst, _role(f.name, v), f.name, v))
    return out


def main() -> None:
    rs = rows()
    counts = {}
    for _, _, _, v in rs:
        counts[PROV[v.prov.value]] = counts.get(PROV[v.prov.value], 0) + 1
    print(f"{len(rs)} fields: " + ", ".join(f"{k} {n}" for k, n in counts.items()) + ".")
    titles = {"shared": "Shared by both instances", "2028": "2028 benchmark (1+1D Z3 domain wall)",
              "2033": "2033 target (3D 2O overlap + Higgs)"}
    for inst in ("shared", "2028", "2033"):
        for role in ("inputs", "printed figures", "comparison records"):
            sel = [(n, v) for i, r, n, v in rs if i == inst and r == role]
            if not sel:
                continue
            print()
            print(f"#### {titles[inst]}: {role} ({len(sel)})")
            print()
            print("| name | value | provenance | source | note |")
            print("|---|---|---|---|---|")
            for n, v in sel:
                print(f"| `{n}` | {_cell(_value(v.value))} | {PROV[v.prov.value]} | {_cell(v.src) or '-'} | {_note(v.note) or '-'} |")


if __name__ == "__main__":
    main()
