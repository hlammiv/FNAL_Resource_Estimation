"""Print the Ch. 4 assumptions table (Markdown) from the live model.

Run from the repository root:

    python applications/ch04_hybrid_lqcd/dump_assumptions.py

Every row is read from estimates.ch04_hybrid_lqcd.Assumptions() at run time: the field name, its
value, its provenance tag, and the source and note strings exactly as the code holds them. Rows are
grouped by the section comments of the dataclass. The "used by" column lists the model functions
(_model_2028, _model_2033, _model_codesign) whose source reads the field as `a.<name>`; a field read
by none of them is a record kept for comparison or a check in __post_init__.
"""

from __future__ import annotations

import dataclasses
import inspect
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from estimates import ch04_hybrid_lqcd as m  # noqa: E402
from estimates.common import Tagged  # noqa: E402

# Display names for the section comments of the dataclass, matched on their leading words.
GROUPS = (
    ("lattice / formulation", "Shared: lattice and formulation"),
    ("2028 tech demo", "2028 benchmark: log det W, free staggered V = 4^4 "
                       "(this section also holds the shot, wall-time and campaign conventions every era uses)"),
    ("the gauged V=4^4 comparison", "2028 comparison: the gauged V = 4^4 circuit "
                                    "(its levers are reused by the co-design row)"),
    ("1000-LQ rung", "2033 target: Tr M^-1 on 8^3x16, and the Banks-Casher demonstration row"),
    ("1e4-LQ (codesign) rung", "Co-design row (post-2033): Tr M^-1 on 24^3x48 "
                               "(p_f here is the fault convention of every era)"),
    ("informational only", "Informational: outside every per-shot count"),
)
ERA_FUNCS = (("2028", m._model_2028), ("2033", m._model_2033), ("codesign", m._model_codesign))


def _groups() -> dict[str, str]:
    """field name -> display group, from the '# --- title ---' comments in the class source."""
    out, current = {}, None
    for line in inspect.getsource(m.Assumptions).splitlines():
        s = line.strip()
        if s.startswith("# ---"):
            title = s.strip("#- ").strip()
            for key, name in GROUPS:
                if title.startswith(key):
                    current = name
                    break
            else:
                raise RuntimeError(f"unmapped section comment: {title!r}")
            continue
        hit = re.match(r"(\w+)\s*:\s*Tagged\b", s)
        if hit:
            out[hit.group(1)] = current
    return out


def _used_by(name: str) -> str:
    pat = re.compile(rf"\ba\.{re.escape(name)}\b")
    eras = [era for era, f in ERA_FUNCS if pat.search(inspect.getsource(f))]
    return ", ".join(eras) if eras else "-"


def _value(v) -> str:
    if isinstance(v, tuple):
        return "(" + ", ".join(repr(x) for x in v) + ")"
    if isinstance(v, str):
        return v
    return repr(v)


def _cell(s: str) -> str:
    return s.replace("|", "\\|").replace("\n", " ")


def main() -> None:
    a = m.Assumptions()
    groups = _groups()
    fields = [f.name for f in dataclasses.fields(a)]
    missing = [n for n in fields if n not in groups]
    if missing:
        raise RuntimeError(f"fields outside any section: {missing}")
    counts: dict[str, int] = {}
    for n in fields:
        v = getattr(a, n)
        assert isinstance(v, Tagged), n
        counts[v.prov.value] = counts.get(v.prov.value, 0) + 1
    print(f"{len(fields)} fields: " + ", ".join(f"{k} {c}" for k, c in counts.items()))
    order = []
    for n in fields:
        if groups[n] not in order:
            order.append(groups[n])
    for g in order:
        print()
        print(f"#### {g}")
        print()
        print("| field | value | provenance | used by | source | note |")
        print("|---|---|---|---|---|---|")
        for n in fields:
            if groups[n] != g:
                continue
            v = getattr(a, n)
            print(f"| `{n}` | {_cell(_value(v.value))} | {v.prov.value} | {_used_by(n)} "
                  f"| {_cell(v.src) or '-'} | {_cell(v.note) or '-'} |")


if __name__ == "__main__":
    main()
