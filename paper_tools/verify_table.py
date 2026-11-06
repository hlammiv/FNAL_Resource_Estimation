#!/usr/bin/env python3
"""Check the chapter table (Table 1.1) against the chapter resource models.

The table lives in two places -- table_standalone.tex and a mirrored copy in
main-overview.tex -- and both restate numbers that are derived in the chapter
benchmarkbox/utilitybox environments. Nothing in LaTeX enforces that they
agree, and in August 2026 four rows had drifted (Chs. 4, 5, 8, 9).

This script re-checks the agreement:

  1. the two table copies carry the same row (application, both generation
     cells, program) for every chapter, Chs. 2-10, and
  2. each cell of both columns prints its chapter box's headline register and
     per-shot hard-operation count at the box's significant figures, and the box
     count is the chapter model's at that precision (estimates/tests/test_table.py,
     which formats with common.tex_hard_ops / lq_text and compares the cell's
     leading "... LQ" and "... T" tokens for equality; 2026-10-01).

Scope, stated plainly so the PASS is not over-read: the header row, the
captions and the two hardware rows (Chs. 11-12) are outside the check, and so
is the descriptive text of a cell (lattice, group, observable). Strict xfails
in test_table.py are filed disagreements, listed in its XFAIL table.
INSTANCES is imported only to confirm that every chapter with a 2033 instance
has a table row.

Exit status is non-zero if anything disagrees, so it can gate a build.

    python3 verify_table.py
"""

import re
import sys
from pathlib import Path

from resources import INSTANCES

DOC = Path(__file__).resolve().parent.parent

# The cell check is model-backed: estimates/tests/test_table.py formats each
# chapter model's LQ and hard-op output for both eras and asserts it is what that
# chapter's row prints in BOTH table copies (strict xfails mark the filed
# disagreements). This replaced a third hand-kept copy of the numbers on 2026-09-28;
# the first-generation column joined the check on 2026-09-29 (R10).
TABLE_TEST = Path(__file__).resolve().parent / "estimates" / "tests" / "test_table.py"

ROW_KEY = {
    2: "Decoding DUNE", 3: "Hunting LNV/LFV", 4: "Quantum-accelerated lattice fermions",
    5: "Viscosity of the perfect liquid", 6: "QCD for colliders",
    7: "After the Big Bang", 8: "Origins of matter", 9: "Chiral lattice fermions",
    10: "Dense QCD",
}


COLUMN = {"2028": 2, "2033": 3}      # `&`-separated field of each generation cell


def fields(path):
    """The `&`-separated fields of each chapter row, Chs. 2-10, of one table copy."""
    text = path.read_text()
    out = {}
    for ch, key in ROW_KEY.items():
        for line in text.splitlines():
            if key in line:
                # split on column separators only: several application names
                # contain an escaped \& ("baryogenesis \& FOPTs")
                cells = [c.strip() for c in re.split(r"(?<!\\)&", line)]
                # Ch. | application | first generation | next generation | first result | campaign \\
                if len(cells) >= 5:
                    out[ch] = cells
                break
    return out


def rows(path, era="2033"):
    """Pull one generation's cell of each chapter row out of a table copy."""
    return {ch: cells[COLUMN[era]] for ch, cells in fields(path).items()}


def norm(cell):
    return re.sub(r"\s+", " ", re.sub(r"\\\\.*$", "", cell)).strip()


def main():
    standalone = fields(DOC / "table_standalone.tex")
    overview = fields(DOC / "main-overview.tex")
    bad = []

    print("Table 1.1 vs. chapter resource models\n" + "-" * 62)

    for ch in sorted(ROW_KEY):
        s, o = standalone.get(ch), overview.get(ch)
        if s is None or o is None:
            bad.append(f"Ch{ch}: row missing from a table copy")
            print(f"  Ch{ch:<3} MISSING ROW")
            continue

        # 1. the two copies must agree with each other, field by field
        #    (field 0 is the chapter number: a \ref in the overview, a digit standalone)
        names = ("application", "first generation", "next generation", "first result", "campaign")
        diff = [n for n, a, b in zip(names, s[1:6], o[1:6]) if norm(a) != norm(b)]
        if diff:
            bad.append(f"Ch{ch}: the two table copies disagree ({', '.join(diff)})")
            for n, a, b in zip(names, s[1:6], o[1:6]):
                if norm(a) != norm(b):
                    print(f"  Ch{ch:<3} COPIES DIFFER in the {n} cell\n"
                          f"        standalone: {norm(a)[:88]}\n"
                          f"        overview:   {norm(b)[:88]}")
            continue

        print(f"  Ch{ch:<3} copies agree (both columns)")

    # 2. the cells must match the chapter MODELS (estimates/tests/test_table.py):
    #    each model's LQ and T, both eras, formatted the way the table writes them
    #    and compared with the cell in both copies. Strict xfails there are the filed
    #    disagreements (NEEDS_AUTHOR.md); an unexpected pass or a new failure
    #    fails this gate.
    import subprocess
    r = subprocess.run([sys.executable, "-m", "pytest", str(TABLE_TEST), "-q", "-p", "no:cacheprovider"],
                       cwd=str(TABLE_TEST.parents[2]), capture_output=True, text=True)
    tail = r.stdout.strip().splitlines()[-1] if r.stdout.strip() else r.stderr.strip()[-200:]
    print(f"  model-backed cell check: {tail}")
    if r.returncode != 0:
        bad.append(f"test_table.py failed: {tail}")

    # 3. every chapter with a 2033 instance should have a row
    have = {i["ch"] for i in INSTANCES if i["era"] == "2033"}
    if have - set(ROW_KEY):
        bad.append(f"chapters with 2033 instances but no table row: {have - set(ROW_KEY)}")

    print("-" * 62)
    if bad:
        print(f"FAIL: {len(bad)} problem(s)")
        for b in bad:
            print("  -", b)
        return 1
    print("PASS: both table copies agree with each other and with the chapter models, both columns.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
