"""Costed resource instances for Chs. 2-10, transcribed from the chapter
benchmarkbox/utilitybox contents.

Single source of truth: fig_resource_landscape.py reads from here, so the
figure and any future table cannot drift apart. Every entry cites the file and
line it came from, so it can be re-checked against the prose.

Units: lq = algorithmic logical qubits (DOE RFI convention -- the algorithmic
register, exclusive of magic-state-distillation factories). t = T/Toffoli
("hard operation") count per shot.

Ranges are (lo, hi); lo == hi means the chapter quotes a single value.
`t_bound="upper"` means the chapter quotes "<=", so the true cost sits at or
below the plotted point.

Three entries formerly conflicted with Table 1.1 (table_standalone.tex,
mirrored in main-overview.tex). Resolved 2026-08-06 in favor of the chapter
boxes, and the table was edited to match; see RECONCILED below.
"""

# Frontier grouping follows Table 1.1's own section headers, so the figure and
# the table group the chapters identically.
FRONTIER = {
    2: "intensity", 3: "intensity", 4: "intensity",
    5: "energy", 6: "energy",
    7: "cosmic", 8: "cosmic", 9: "cosmic", 10: "cosmic",
}

CHAPTER_SHORT = {
    2: r"$\nu$--nucleus",
    3: r"$\mu\to e$ / $0\nu\beta\beta$",
    4: "hybrid LQCD",
    5: "QGP transport",
    6: "collider",
    7: "early universe",
    8: "baryogenesis",
    9: "chiral gauge",
    10: "finite density",
}


# Instances are GENERATED. Run `python3 -m estimates --resources` (from doc/scripts)
# to regenerate resources.json from the chapter models in estimates/; this module
# only loads it. The hand-transcribed list this replaced (26 rows, 2026-08-06 ..
# 2026-09-28) is preserved in backups/ and in git history.
_DEFAULTS = dict(t_bound=None, conditional=False, codesign=False, reconciled=False, src="")


def _load():
    """Split generated rows into costed INSTANCES and LQ_ONLY (a row whose model
    emits no T-count -- t == (0, 0) -- is a register quoted without a costed
    depth, drawn as an LQ-only rule rather than invented a depth for)."""
    import json
    from pathlib import Path
    p = Path(__file__).resolve().parent / "resources.json"
    rows = json.loads(p.read_text())["instances"]
    costed, lq_only = [], []
    for r in rows:
        d = dict(_DEFAULTS)
        d.update(r)
        d["lq"] = tuple(d["lq"])
        d["t"] = tuple(d["t"])
        # The figure draws `codesign` rows in the beyond-RFI band (>= 1e12 T).
        # The models tag `codesign` by era; a co-design-era row priced below
        # 1e12 (Ch. 5 rails, Ch. 6 stretch, Ch. 3 jj55/Mo pairs) belongs on the
        # main axes, so the drawing flag is derived from the T-count here.
        d["era_label"] = d["era"]
        d["codesign"] = bool(d["t"][0] >= 1e12)
        if d["t"] == (0, 0) or d["t"][1] <= 0:
            lq_only.append(dict(ch=d["ch"], era=d["era"], label=d["label"],
                                lq=d["lq"][1], lq_lo=d["lq"][0], note=d.get("note") or "T not costed", src=d["src"]))
        else:
            costed.append(d)
    return costed, lq_only


INSTANCES, LQ_ONLY = _load()


# DOE 2026 RFI anchors (RFI 892431-26-RFI-0001), as used throughout the report.
ENVELOPES = {
    "2028": dict(lq=(150, 250), t=(1e4, 1e5), label="2028 RFI envelope"),
    "2033": dict(lq=(250, 1000), t=(1e5, 1e9), label="2033 RFI envelope"),
}

# ---------------------------------------------------------------------------
# Resolved 2026-08-06: Table 1.1 was edited to match the chapter boxes in all
# three cases, so figure and table now agree. Kept as a record of what moved
# and why; verify_table.py re-checks the agreement on demand.
RECONCILED = {
    4: ("Table said '1k LQ' for the 2033 cell, which is the machine era, not "
        "the algorithmic register the caption declares. Table now reads "
        "'85 LQ register in a ~1k-LQ machine'; the caption's three-rung "
        "ladder is now labeled as machine-era-indexed."),
    5: ("Table said '1--1.1k LQ'; the chapter box quotes a single 1000. Table "
        "now reads '1k LQ, 1e9 T'."),
    8: ("Table gave volume 64x32 and no T-count, so the caption's default rule "
        "silently assigned it the 1e9 column anchor. Chapter says 10x8 and "
        "2e7 T -- two decades inside. Table now states both explicitly."),
    9: ("Table said SU(3) overlap on 4^3, 700 LQ, 1e9-1e10 T; chapter (then) said 3D "
        "SU(2)/2O overlap on 3^3, 1000 LQ, 2e8 T. Table now matches. "
        "(Current: 1010 LQ, 8.2e9 T on the sphaleron arm, r24 E28: kappa 30, evolution to t = 2 fm = 10a, 450 steps; condensate arm 2.2e8 T. 1010 LQ, 2.0e11 T at r22 E26: single-particle overlap kernel lifted to Fock space, Trotter steps from the state-dependent estimate; 1005 LQ, 1.8e13 T at r21 E25 with the steps from the commutator bound; 1015 LQ, 5.7e11 central at r20 E23; 1082 LQ at r19; ~1.5e10 T after rulings R5/VERTEX of 2026-09-28; ~8e9 after R-032.)"),
}


def geometry(inst):
    """Classify an instance by how well its resources are pinned down.

    This is the visual variable the figure encodes: definiteness, not size.
    """
    lo_lq, hi_lq = inst["lq"]
    lo_t, hi_t = inst["t"]
    wide_lq, wide_t = hi_lq > lo_lq, hi_t > lo_t
    if inst["codesign"]:
        return "codesign"
    if inst["conditional"]:
        return "conditional"
    if inst["t_bound"]:
        return "bound"
    if wide_lq and wide_t:
        return "blob"
    if wide_t:
        return "capsule_t"
    if wide_lq:
        return "capsule_lq"
    return "point"


def print_table():
    hdr = f"{'Ch':>3} {'era':>9} {'LQ':>13} {'T-gates':>25} {'geometry':>11}  label"
    print(hdr)
    print("-" * len(hdr))
    for i in INSTANCES:
        lq = f"{i['lq'][0]:g}" if i["lq"][0] == i["lq"][1] else f"{i['lq'][0]:g}-{i['lq'][1]:g}"
        t = f"{i['t'][0]:.3g}" if i["t"][0] == i["t"][1] else f"{i['t'][0]:.3g}-{i['t'][1]:.3g}"
        if i["t_bound"] == "upper":
            t = "<= " + t
        flag = "  "
        lbl = i["label"].replace("$", "").replace("\\", "")
        print(f"{i['ch']:>3} {i['era']:>9} {lq:>13} {t:>25} {geometry(i):>11}{flag} {lbl}")
    for i in LQ_ONLY:
        print(f"{i['ch']:>3} {i['era']:>9} {i['lq']:>13g} {i['note']:>25} {'lq_only':>11}   "
              + i["label"].replace("$", "").replace("\\", ""))
    print(f"\n{len(INSTANCES)} costed instances + {len(LQ_ONLY)} LQ-only.")
    print("\nReconciled with Table 1.1 on 2026-08-06:")
    for ch, why in RECONCILED.items():
        print(f"  Ch{ch}: {why}")


if __name__ == "__main__":
    print_table()
