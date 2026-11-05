#!/usr/bin/env python3
"""Accept/reject/strip tracked-change markup (\\tadd/\\tdel/\\trep/\\tnote).

Companion to trackchanges.tex. Every tracked macro must carry a same-line
%%TC:<ID> tag after its final closing brace; IDs key into TRACKED_CHANGES.md.

  --list                 enumerate all tracked edits (ID, file:line, kind)
  --accept-all           resolve every edit in favor of the new text
  --reject ID [ID ...]   resolve the named IDs in favor of the OLD text,
                         accept everything else
  --dry-run              print unified diffs instead of writing files

Accepting removes the markup: \\tadd -> content, \\tdel -> (gone),
\\trep{new}{old} -> new, \\tnote -> (gone); %%TC: tags and %%TC-OLD: comment
lines are dropped, and the \\input{trackchanges} hook is removed from
main.tex.  Rejecting inverts: \\tadd -> (gone), \\tdel -> content,
\\trep -> old.  An ID whose old text lives only in %%TC-OLD: comments (the
applied-clean fragile-context cases) cannot be auto-rejected; the script
keeps those lines and tells you to restore by hand.

Safety: character-scan brace matcher (escape-aware). The script refuses to
run if it finds unbalanced braces, a tracked macro without a %%TC: tag, a
nested tracked macro, or an unescaped % inside a tracked argument.
"""

import argparse
import difflib
import re
import sys
from pathlib import Path

DOC = Path(__file__).resolve().parent.parent
MACROS = ("tadd", "tdel", "trep", "tnote")
TAG_RE = re.compile(r"%%TC:([A-Za-z0-9_\-]+)")
MACRO_RE = re.compile(r"\\(tadd|tdel|trep|tnote)\b")


def tex_files():
    files = [DOC / "main.tex", DOC / "main-overview.tex",
             DOC / "table_standalone.tex", DOC / "commands.tex"]
    files += sorted((DOC / "applications").glob("*.tex"))
    return [f for f in files if f.exists()]


def scan_argument(text, i, path, lineno_at):
    """text[i] must be '{'; return (content, index_after_closing_brace)."""
    assert text[i] == "{"
    depth, j, out = 1, i + 1, []
    while j < len(text):
        c = text[j]
        if c == "\\" and j + 1 < len(text):
            out.append(text[j:j + 2])
            j += 2
            continue
        if c == "%":
            sys.exit(f"ERROR: unescaped %% inside tracked argument at "
                     f"{path}:{lineno_at(j)} -- not supported, fix the markup")
        if c == "{":
            depth += 1
        elif c == "}":
            depth -= 1
            if depth == 0:
                return "".join(out), j + 1
        out.append(c)
        j += 1
    sys.exit(f"ERROR: unbalanced braces starting at {path}:{lineno_at(i)}")


def parse_file(path):
    """Return list of edits: dict(kind, src, args, id, start, end, line)."""
    text = path.read_text()

    def lineno_at(idx):
        return text.count("\n", 0, idx) + 1

    edits = []
    for m in MACRO_RE.finditer(text):
        # skip definitions in trackchanges.tex-style contexts (newcommand)
        line_start = text.rfind("\n", 0, m.start()) + 1
        if "newcommand" in text[line_start:m.start()]:
            continue
        kind = m.group(1)
        i = m.end()
        src = "team"
        if i < len(text) and text[i] == "[":
            close = text.index("]", i)
            src = text[i + 1:close]
            i = close + 1
        nargs = 2 if kind == "trep" else 1
        args = []
        for _ in range(nargs):
            while i < len(text) and text[i] in " \t\n%":
                if text[i] == "%":  # tag or comment between args: not allowed
                    sys.exit(f"ERROR: comment between arguments at "
                             f"{path}:{lineno_at(i)}")
                i += 1
            if i >= len(text) or text[i] != "{":
                sys.exit(f"ERROR: expected '{{' after \\{kind} at "
                         f"{path}:{lineno_at(i)}")
            content, i = scan_argument(text, i, path, lineno_at)
            if MACRO_RE.search(content):
                sys.exit(f"ERROR: nested tracked macro at "
                         f"{path}:{lineno_at(m.start())} -- not supported")
            args.append(content)
        # tag must follow on the same line as the final closing brace
        eol = text.find("\n", i)
        eol = len(text) if eol == -1 else eol
        tag_m = TAG_RE.search(text[i:eol])
        if not tag_m:
            sys.exit(f"ERROR: \\{kind} at {path}:{lineno_at(m.start())} has "
                     f"no same-line %%TC:<ID> tag")
        # the tag is a TeX comment: anything after it on the line is silently
        # dropped from the build, so the tag must terminate the line
        trailing = text[i + tag_m.end():eol]
        if trailing.strip():
            sys.exit(f"ERROR: content after tag {tag_m.group(0)} at "
                     f"{path}:{lineno_at(m.start())} is commented out by it "
                     f"-- move it to the next line: {trailing.strip()[:60]!r}")
        # span covers only the macro; the %%TC: tag is stripped separately so
        # text between the closing brace and the tag survives resolution
        edits.append(dict(kind=kind, src=src, args=args, id=tag_m.group(1),
                          start=m.start(), end=i,
                          line=lineno_at(m.start())))
    return text, edits


def resolve(text, edits, reject_ids):
    """Replace each edit span with its resolution, back to front."""
    manual = []
    for e in sorted(edits, key=lambda e: e["start"], reverse=True):
        rejected = e["id"] in reject_ids
        kind = e["kind"]
        if kind == "tadd":
            repl = "" if rejected else e["args"][0]
        elif kind == "tdel":
            repl = e["args"][0] if rejected else ""
        elif kind == "trep":
            repl = e["args"][1] if rejected else e["args"][0]
        else:  # tnote
            repl = ""
            if rejected:
                manual.append(e)
        text = text[:e["start"]] + repl + text[e["end"]:]
    return text, manual


def strip_lines(text, reject_ids, path):
    """Drop %%TC-OLD: lines (unless their ID is rejected) and stray tags."""
    out, kept_old = [], []
    for line in text.split("\n"):
        m = re.match(r"\s*%%TC-OLD:([A-Za-z0-9_\-]+)", line)
        if m:
            if m.group(1) in reject_ids:
                kept_old.append((path, m.group(1)))
                out.append(line)
            continue
        out.append(TAG_RE.sub("", line))
    return "\n".join(out), kept_old


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--accept-all", action="store_true")
    ap.add_argument("--reject", nargs="+", default=[], metavar="ID")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()

    all_edits = []
    parsed = {}
    for path in tex_files():
        text, edits = parse_file(path)
        parsed[path] = (text, edits)
        all_edits += [(path, e) for e in edits]

    if a.list or not (a.accept_all or a.reject):
        for path, e in sorted(all_edits, key=lambda x: (str(x[0]), x[1]["line"])):
            print(f"{e['id']:<12} {path.relative_to(DOC)}:{e['line']:<5} "
                  f"\\{e['kind']}[{e['src']}]")
        print(f"-- {len(all_edits)} tracked edits")
        return

    reject_ids = set(a.reject)
    known = {e["id"] for _, e in all_edits}
    unknown = reject_ids - known
    if unknown:
        sys.exit(f"ERROR: unknown IDs: {', '.join(sorted(unknown))}")

    warnings = []
    for path, (text, edits) in parsed.items():
        new, manual = resolve(text, edits, reject_ids)
        new, kept = strip_lines(new, reject_ids, path)
        for p, i in kept:
            warnings.append(f"NOTE: {p.relative_to(DOC)} keeps %%TC-OLD:{i} "
                            f"lines -- restore that rejected edit by hand")
        for e in manual:
            warnings.append(f"NOTE: rejected \\tnote {e['id']} removed; if it "
                            f"marked an applied-clean edit, revert by hand")
        if path.name == "main.tex":
            new = new.replace("\\input{trackchanges}\n", "")
        if new != text:
            if a.dry_run:
                sys.stdout.writelines(difflib.unified_diff(
                    text.splitlines(True), new.splitlines(True),
                    fromfile=str(path), tofile=f"{path} (resolved)"))
            else:
                path.write_text(new)
                print(f"resolved: {path.relative_to(DOC)}")
    for w in warnings:
        print(w)
    if not a.dry_run:
        print("Done. trackchanges.tex itself can now be deleted; "
              "rebuild and re-run verify_table.py.")


if __name__ == "__main__":
    main()
