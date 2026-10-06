"""Print the Ch. 2 assumptions table (Markdown) from the live model.

Run from anywhere:  python applications/ch02_neutrino_nucleus/dump_assumptions.py

Every row comes from estimates.ch02_nu_nucleus.Assumptions(). Rows are grouped by which
code path actually reads the field: the script runs model(a, era) for each era on a proxy
that records attribute reads, so the grouping is measured, not hand-sorted. Values are the
Python repr of Tagged.value; sources are printed verbatim (newlines folded, '|' escaped
for Markdown). The note column is printed for group A only; the notes of the other fields
are in ch02_legacy.py.
"""
import dataclasses
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from estimates import ch02_nu_nucleus as m  # noqa: E402
from estimates.common import Provenance  # noqa: E402


class _Recorder:
    """Wraps an Assumptions instance and records which fields are read."""

    def __init__(self, a):
        object.__setattr__(self, "_a", a)
        object.__setattr__(self, "seen", set())

    def __getattr__(self, name):
        self.seen.add(name)
        return getattr(self._a, name)


def fields_read(era: str) -> set[str]:
    rec = _Recorder(m.Assumptions())
    m.model(rec, era)
    names = {f.name for f in dataclasses.fields(m.Assumptions)}
    return rec.seen & names


def _cell(s: str) -> str:
    return " ".join(str(s).split()).replace("|", "\\|")


def _prov(v) -> str:
    if v.prov is Provenance.CITED:
        return f"Cited ({v.src})"
    if v.prov is Provenance.ASSUMED:
        return "Assumed"
    if v.prov is Provenance.STATED:
        return "Stated-not-derived"
    return v.prov.value


def _src(v) -> str:
    # Cited rows carry their key in the provenance column; others show their source here.
    return "" if v.prov is Provenance.CITED else v.src


def table(names: list[str], a, notes: bool = False) -> str:
    head = "| name | value | provenance | source |" + (" note |" if notes else "")
    out = [head, "|---|---|---|---|" + ("---|" if notes else "")]
    for n in names:
        v = getattr(a, n)
        row = f"| `{n}` | `{v.value!r}` | {_cell(_prov(v))} | {_cell(_src(v))} |"
        out.append(row + (f" {_cell(v.note)} |" if notes else ""))
    return "\n".join(out)


def main() -> None:
    a = m.Assumptions()
    order = [f.name for f in dataclasses.fields(a)]
    bench = fields_read("2028") | fields_read("2033")
    codesign = fields_read("codesign")
    groups = [
        ("A. Read by the 2028 and 2033 benchmarks (`model`, ch02_nu_nucleus.py)",
         [n for n in order if n in bench]),
        ("B. Read by the 12C co-design scaling row (`_model_codesign`, ch02_legacy.py)",
         [n for n in order if n in codesign and n not in bench]),
        ("C. Inherited from ch02_legacy.Assumptions, read by no exported number",
         [n for n in order if n not in bench and n not in codesign]),
    ]
    counts = {}
    for n in order:
        p = _prov(getattr(a, n)).split(" ")[0]
        counts[p] = counts.get(p, 0) + 1
    print(f"{len(order)} fields: " + ", ".join(f"{k} {v}" for k, v in counts.items()))
    for i, (title, names) in enumerate(groups):
        print(f"\n#### {title}: {len(names)} fields\n")
        print(table(names, a, notes=(i == 0)))


if __name__ == "__main__":
    main()
