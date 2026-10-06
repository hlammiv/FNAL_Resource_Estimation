"""Print the Ch. 7 Assumptions as Markdown tables, straight from the live model.

    python applications/ch07_curved_space/dump_assumptions.py

One table per section of the `Assumptions` dataclass in estimates/ch07_curved_space.py (the sections are
the `# ---- ... ----` comment headers in the class body). Columns: field name, value as stored, provenance
tag, source, note. Nothing is rounded or rewritten, so each row matches the code exactly.
"""

from __future__ import annotations

import inspect
import re
import sys
from dataclasses import fields
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from estimates import ch07_curved_space as ch7          # noqa: E402
from estimates.common import Provenance                 # noqa: E402

TAG = {Provenance.CITED: "Cited", Provenance.ASSUMED: "Assumed",
       Provenance.STATED: "Stated-not-derived", Provenance.UNCITED: "Uncited"}


def sections() -> dict[str, str]:
    """field name -> the section header it sits under in the class body."""
    src = inspect.getsource(ch7.Assumptions)
    out, current = {}, "(no header)"
    for line in src.splitlines():
        m = re.match(r"\s*# ----\s*(.*?)\s*-*\s*$", line)
        if m:
            current = m.group(1)
            continue
        m = re.match(r"\s{4}(\w+)\s*:\s*Tagged\b", line)
        if m:
            out[m.group(1)] = current
    return out


def cell(s) -> str:
    return str(s).replace("|", "\\|").replace("\n", " ")


def main() -> None:
    a = ch7.Assumptions()
    sec = sections()
    groups: dict[str, list] = {}
    for f in fields(a):
        groups.setdefault(sec.get(f.name, "(no header)"), []).append((f.name, getattr(a, f.name)))
    counts = {}
    for f in fields(a):
        t = TAG[getattr(a, f.name).prov]
        counts[t] = counts.get(t, 0) + 1
    print(f"{len(fields(a))} fields: " + ", ".join(f"{k} {v}" for k, v in counts.items()))
    for title, rows in groups.items():
        print(f"\n#### {title}\n")
        print("| field | value | provenance | source | note |")
        print("|---|---|---|---|---|")
        for name, v in rows:
            print(f"| `{name}` | {cell(repr(v.value))} | {TAG[v.prov]} | {cell(v.src)} | {cell(v.note)} |")


if __name__ == "__main__":
    main()
