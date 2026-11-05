#!/usr/bin/env python3
"""Portfolio-level resource landscape: logical qubits vs T-gates for every
costed instance in Chs. 2-10.

Marker geometry encodes how well each instance is pinned down -- a filled point
means the chapter quotes one number on both axes; a capsule or open box means it
quotes a range; an arrow means it quotes an upper bound. Nested shaded regions
are the DOE 2026 RFI "fits" envelopes, so an instance inside a box is one the
corresponding machine can actually run.

Both axes are logarithmic: the 2028 instances run from 18 to 250 LQ and
would be unreadable on a linear register axis.

Data lives in resources.py. Run that module directly to dump the table.
"""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle

import fnalfig as ff
from resources import ENVELOPES, FRONTIER, INSTANCES, LQ_ONLY, geometry

COLOR = {"intensity": ff.BLUE, "energy": ff.ORANGE, "cosmic": ff.GREEN}

# Main panel spans the real portfolio; the narrow right panel carries the
# co-design horizon at 1e12-1e15 without letting it compress everything else.
XMAIN = (3e3, 1e12)   # from 3e3 so the Ch. 3 2028 insertion (5.4e3, R13) and the inset's left edge (4e3) are on the axes; to 1e12 (E19) so every
                      # row priced below 1e12 is on the axes (E28: Ch. 9's 2033 target is 8.2e9, in this panel; 2.0e11 at E26; it was 1.8e13 at r21 E25, in the right-hand panel)
XCODE = (8e11, 1e16)   # left pad so a rung near 1e12 is not clipped at the break (the Ch. 4 rung was 1.03e12 at R11; 2.6e11 after R3, now in the main panel). A main-panel range that runs past
                       # 1e12 continues into this panel (none does after E25 or E26; Ch. 9 2033 was 1.6e10-1.3e12 at r20); the 8e11-1e12 overlap only repeats such a capsule's end
YLIM = (13, 3400)      # floor below the smallest register, 18 LQ (Ch. 4 2028, ruling R7)

# --- One placement rule per geometry, applied identically in the main panel,
# --- the co-design panel and the inset, so the three cannot drift apart.
# --- Each geometry names an ordered list of candidate slots; the first that
# --- does not collide with an already-placed label or another mark is taken.
# --- Deterministic, so a re-run reproduces the figure exactly.
ABOVE, RIGHT, LEFT, BELOW = (0, 8), (11, 0), (-11, 0), (0, -8)
# Diagonals are tried after the four cardinals: at the 2028 rung five marks sit
# within twenty logical qubits of one another and four slots are not enough.
NE, NW, SE, SW = (9, 6), (-9, 6), (9, -6), (-9, -6)
LABEL_SLOTS = {
    "point":       [ABOVE, RIGHT, LEFT, BELOW, NE, NW, SE, SW],
    "bound":       [ABOVE, RIGHT, LEFT, BELOW, NE, NW, SE, SW],
    "capsule_t":   [ABOVE, RIGHT, LEFT, BELOW, NE, NW, SE, SW],
    "capsule_lq":  [RIGHT, LEFT, ABOVE, BELOW, NE, SE, NW, SW],  # top cap busy
    "blob":        [NE, NW, SE, SW, ABOVE, RIGHT, LEFT, BELOW],  # cardinals busy
    "conditional": [ABOVE, RIGHT, LEFT, BELOW, NE, NW, SE, SW],
    "codesign":    [ABOVE, RIGHT, LEFT, BELOW, NE, NW, SE, SW],
}
# Hand-placed first choices where the automatic order lands a label between two
# marks (E21: Ch. 10's 2028 label sat between the Ch. 9 and Ch. 10 points).
LABEL_PREFER = {(10, "2028"): NE}
CHAR_PT = 3.9      # width of one label character at 6.5pt
LINE_PT = 7.0      # label height

# Upper bounds are drawn as a short fixed stub rather than a long trailing
# arrow: two 2028 bounds sit ten logical qubits apart, and long arrows overlap.
BOUND_DECADES = 0.28

MARK_MS = 4.2          # filled point
CAP_MS = 2.6           # center dot on a capsule
CAP_LW = 3.4


def draw_geometry(inst):
    """Geometry as drawn. A beyond-RFI (>= 1e12 T) row is placed in the right-hand
    panel, and the panel is what marks it; its mark still encodes how tightly it
    is pinned down (range, single value), so it is not drawn as a dashed
    'conditional' box it is not."""
    return geometry(dict(inst, codesign=False))


def mark_geometry(inst):
    """Marker anchor and the point the label hangs off, in data coords.

    Returned anchor is where the chapter number attaches; keeping this in one
    place is what lets the inset reuse the main panel's rule verbatim.
    """
    lo_t, hi_t = inst["t"]
    lo_q, hi_q = inst["lq"]
    g = draw_geometry(inst)
    mid_t, mid_q = np.sqrt(lo_t * hi_t), np.sqrt(lo_q * hi_q)
    if inst.get("t_representative") and lo_t <= inst["t_representative"] <= hi_t:
        # a chapter that prints a central value inside its band (Ch. 9 2033 at r20, center 5.7e11 of 1.6e10-1.3e12; a single value since E25, 2.0e11 at E26, 8.2e9 since E28)
        # is marked there, not at the band's geometric middle
        mid_t = float(inst["t_representative"])
    if g == "bound":
        return g, (hi_t, mid_q), (hi_t, mid_q)
    if g == "capsule_lq":
        return g, (mid_t, mid_q), (mid_t, mid_q)
    if g in ("conditional", "codesign"):
        return g, (mid_t, mid_q), (mid_t, _box(inst)[3])   # label above the box top
    return g, (mid_t, mid_q), (mid_t, mid_q)


BOX_MIN_T = 0.07     # decades: half-width of a conditional box quoted at one T
BOX_MIN_Q = 0.025    # decades: half-height of one quoted at one register


def _box(inst):
    """Drawn extent of a conditional/co-design box, with degenerate sides padded."""
    lo_t, hi_t = inst["t"]
    lo_q, hi_q = inst["lq"]
    if np.log10(hi_t / lo_t) < 2 * BOX_MIN_T:
        m = np.sqrt(lo_t * hi_t)
        lo_t, hi_t = m * 10 ** -BOX_MIN_T, m * 10 ** BOX_MIN_T
    if np.log10(hi_q / lo_q) < 2 * BOX_MIN_Q:
        m = np.sqrt(lo_q * hi_q)
        lo_q, hi_q = m * 10 ** -BOX_MIN_Q, m * 10 ** BOX_MIN_Q
    return lo_t, hi_t, lo_q, hi_q


def draw_mark(ax, inst, bound_decades=BOUND_DECADES):
    """Draw one instance's mark. Labels are added separately, after clustering."""
    c = COLOR[FRONTIER[inst["ch"]]]
    lo_t, hi_t = inst["t"]
    lo_q, hi_q = inst["lq"]
    g, (mx, my), _ = mark_geometry(inst)

    if g == "point":
        ax.plot(mx, my, "o", ms=MARK_MS, color=c, mec="white", mew=0.5, zorder=5)
    elif g == "capsule_t":
        ax.plot([lo_t, hi_t], [my, my], lw=CAP_LW, color=c, alpha=0.55,
                solid_capstyle="round", zorder=4)
        ax.plot(mx, my, "o", ms=CAP_MS, color=c, zorder=5)
    elif g == "capsule_lq":
        ax.plot([mx, mx], [lo_q, hi_q], lw=CAP_LW, color=c, alpha=0.55,
                solid_capstyle="round", zorder=4)
        ax.plot(mx, my, "o", ms=CAP_MS, color=c, zorder=5)
    elif g == "blob":
        ax.plot([lo_t, hi_t], [my, my], lw=CAP_LW, color=c, alpha=0.55,
                solid_capstyle="round", zorder=4)
        ax.plot([mx, mx], [lo_q, hi_q], lw=CAP_LW, color=c, alpha=0.55,
                solid_capstyle="round", zorder=4)
        ax.plot(mx, my, "o", ms=CAP_MS, color=c, zorder=5)
    elif g == "bound":
        ax.plot(mx, my, "o", ms=MARK_MS, color=c, mec="white", mew=0.5, zorder=5)
        ax.annotate("", xy=(mx * 10 ** -bound_decades, my), xytext=(mx, my),
                    arrowprops=dict(arrowstyle="-|>", color=c, lw=0.8,
                                    mutation_scale=6, shrinkA=1.5, shrinkB=0),
                    zorder=4)
    elif g == "codesign" and lo_t == hi_t and lo_q == hi_q:
        # a single value on both axes has no box to draw (Ch. 4 co-design, 92 LQ at 2e12 T)
        ax.plot(mx, my, "o", ms=MARK_MS, color=c, mec="white", mew=0.5, zorder=5)
    elif g in ("conditional", "codesign"):
        # A conditional row quoted as a single value on one axis (Ch. 5 2028 at
        # 4.36e5 T) would otherwise be a zero-width box and
        # vanish; give each degenerate side a fixed half-width in log space.
        bt0, bt1, bq0, bq1 = _box(inst)
        ax.add_patch(Rectangle((bt0, bq0), bt1 - bt0, bq1 - bq0,
                               facecolor=c, alpha=0.10, edgecolor=c,
                               ls="--", lw=0.9, zorder=4))


def draw_labels(ax, instances, merge_pts=9.0, lines_block=True):
    """Attach chapter numbers using the per-geometry rule.

    Marks whose label anchors land within `merge_pts` of each other in display
    space share one label ("5,6,7") rather than being scattered by hand -- at
    1e9 T three chapters sit within a few points of one another, and any
    per-instance nudging there is arbitrary and unreproducible.
    """
    ax.figure.canvas.draw()
    entries, obstacles = [], []
    pt = ax.figure.dpi / 72.0        # points -> display pixels

    # Text already on the axes -- the frontier and geometry keys, the envelope
    # labels, the panel titles -- is just as much an obstacle as a mark. Sample
    # each existing label's box so a chapter number never lands on top of one.
    # (Only pre-existing text is scanned: this runs before any number is placed.)
    rend = ax.figure.canvas.get_renderer()
    for txt in list(ax.texts):
        try:
            bb = txt.get_window_extent(renderer=rend)
        except Exception:
            continue
        if bb.width <= 0 or bb.height <= 0:
            continue
        obstacles += [(bb.x0 + fx * bb.width, bb.y0 + fy * bb.height)
                      for fx in np.linspace(0, 1, 6)
                      for fy in np.linspace(0, 1, 3)]

    for ln in (list(ax.lines) if lines_block else []):
        if ln.get_linestyle() in ("-", "solid") or len(ln.get_xdata()) < 2:
            continue
        pts = ln.get_transform().transform(
            np.column_stack([np.asarray(ln.get_xdata(), float),
                             np.asarray(ln.get_ydata(), float)]))
        for a, b in zip(pts[:-1], pts[1:]):
            obstacles += [tuple(a + f * (b - a)) for f in np.linspace(0, 1, 14)]

    for inst in instances:
        g, (mx, my), anchor = mark_geometry(inst)
        px, py = ax.transData.transform(anchor)
        entries.append(dict(inst=inst, g=g, anchor=anchor, px=px, py=py))

        # Sample the whole drawn extent, not just the anchor: a bar or an
        # arrow shaft is just as much an obstacle as the point at its center.
        lo_t, hi_t = inst["t"]
        lo_q, hi_q = inst["lq"]
        if g == "bound":
            span = [(mx * 10 ** (-BOUND_DECADES * f), my)
                    for f in np.linspace(0, 1, 7)]
        elif g == "capsule_t":
            span = [(t, my) for t in np.geomspace(lo_t, hi_t, 9)]
        elif g == "capsule_lq":
            span = [(mx, q) for q in np.geomspace(lo_q, hi_q, 9)]
        elif g == "blob":
            span = ([(t, my) for t in np.geomspace(lo_t, hi_t, 7)]
                    + [(mx, q) for q in np.geomspace(lo_q, hi_q, 7)])
        elif g in ("conditional", "codesign"):
            b = _box(inst)
            span = [(t, q) for t in np.geomspace(b[0], b[1], 5)
                    for q in np.geomspace(b[2], b[3], 3)]
        else:
            span = [(mx, my)]
        obstacles += [tuple(ax.transData.transform(p)) for p in span]

    used, groups = set(), []
    for i, e in enumerate(entries):
        if i in used:
            continue
        grp = [e]
        used.add(i)
        for j, f in enumerate(entries):
            if j in used:
                continue
            if FRONTIER[f["inst"]["ch"]] != FRONTIER[e["inst"]["ch"]]:
                continue
            if np.hypot(f["px"] - e["px"], f["py"] - e["py"]) <= merge_pts * pt:
                grp.append(f)
                used.add(j)
        groups.append(grp)

    # Place topmost-first so the order is stable across runs and panels.
    groups.sort(key=lambda g: -np.mean([e["py"] for e in g]))
    placed = []                                    # occupied label rects

    def rect(cx, cy, w, h):
        return (cx - w / 2, cy - h / 2, cx + w / 2, cy + h / 2)

    ab = ax.get_window_extent()

    def clashes(r):
        # a label that leaves the axes is clipped by the frame (Ch. 6 stretch
        # at 2404 LQ used to lose its number off the top edge)
        if r[0] < ab.x0 or r[2] > ab.x1 or r[1] < ab.y0 or r[3] > ab.y1:
            return True
        for q in placed:
            if not (r[2] < q[0] or r[0] > q[2] or r[3] < q[1] or r[1] > q[3]):
                return True
        pad = 1.5 * pt
        for mx, my in obstacles:                   # keep labels off the marks
            if r[0] - pad <= mx <= r[2] + pad and r[1] - pad <= my <= r[3] + pad:
                return True
        return False

    for grp in groups:
        chs = sorted({g["inst"]["ch"] for g in grp})
        text = ",".join(str(c) for c in chs)
        col = COLOR[FRONTIER[chs[0]]]      # merging is frontier-local
        px = float(np.mean([g["px"] for g in grp]))
        py = float(np.mean([g["py"] for g in grp]))
        w, h = (len(text) * CHAR_PT + 2.0) * pt, LINE_PT * pt

        g0 = grp[0]["g"]
        lq0 = grp[0]["inst"]["lq"]
        if g0 == "capsule_lq" and lq0[1] / lq0[0] < 1.05:
            # a 171-172 LQ "range" draws as a dot, so its top cap is free: use
            # the point order (label above) rather than beside a neighbor
            g0 = "point"
        base = LABEL_SLOTS[g0]
        key = (grp[0]["inst"]["ch"], grp[0]["inst"]["era"])
        if len(grp) == 1 and key in LABEL_PREFER:
            base = [LABEL_PREFER[key]] + [b for b in base if b != LABEL_PREFER[key]]
        for scale in (1.0, 1.7, 2.4, 3.2, 4.2):
            for bx, by in base:
                dx, dy = bx * scale, by * scale
                cx = px + dx * pt + (0 if dx == 0 else np.sign(dx) * w / 2)
                cy = py + dy * pt + (0 if dy == 0 else np.sign(dy) * h / 2)
                if not clashes(rect(cx, cy, w, h)):
                    break
            else:
                continue
            break
        else:                                # genuinely nowhere: take the top
            dx, dy = base[0][0] * 4.2, base[0][1] * 4.2
            cx = px + dx * pt
            cy = py + dy * pt + np.sign(dy) * h / 2
        placed.append(rect(cx, cy, w, h))

        ha = "center" if dx == 0 else ("left" if dx > 0 else "right")
        va = "bottom" if dy > 0 else ("top" if dy < 0 else "center")
        ax.annotate(text, xy=grp[0]["anchor"] if len(grp) == 1 else
                    (float(np.mean([g["anchor"][0] for g in grp])),
                     float(np.mean([g["anchor"][1] for g in grp]))),
                    xytext=(dx, dy), textcoords="offset points",
                    ha=ha, va=va, fontsize=6.5, color=col, zorder=6)


def envelopes(ax):
    """Nested 'fits' regions: an instance inside a corner is one that machine
    runs. Only the two binding edges are drawn, so each reads as a budget
    corner rather than as a background panel."""
    for era, alpha in (("2033", 0.045), ("2028", 0.075)):
        env = ENVELOPES[era]
        t_hi, q_hi = env["t"][1], env["lq"][1]
        ax.add_patch(Rectangle((XMAIN[0], YLIM[0]), t_hi - XMAIN[0], q_hi - YLIM[0],
                               facecolor=ff.FNAL_BLUE, alpha=alpha,
                               edgecolor="none", zorder=1))
        ax.plot([XMAIN[0], t_hi], [q_hi, q_hi], color=ff.FNAL_BLUE, lw=0.8,
                ls=(0, (4, 2)), zorder=2)
        ax.plot([t_hi, t_hi], [YLIM[0], q_hi], color=ff.FNAL_BLUE, lw=0.8,
                ls=(0, (4, 2)), zorder=2)
    # Each label rides its own vertical budget edge in the low-register band.
    # That strip is empty: the inset owns the upper left, the key the lower
    # right, and no instance except Ch. 4's 2028 point sits below 100 LQ.
    for era, t_hi in (("2028", 1e5), ("2033", 1e9)):
        ax.annotate(f"{era} RFI target", xy=(t_hi * 1.16, 45), fontsize=6.5,
                    color=ff.FNAL_BLUE, va="center", ha="left", rotation=90)


def main():
    ff.setup()
    fig, (ax, axc) = plt.subplots(
        1, 2, figsize=(6.5, 4.0), dpi=300, sharey=True,
        gridspec_kw=dict(width_ratios=[5, 1], wspace=0.06),
    )

    envelopes(ax)

    for inst in INSTANCES:
        draw_mark(axc if inst["codesign"] else ax, inst)
        if not inst["codesign"] and inst["t"][1] > XMAIN[1]:
            # a main-panel range that runs past the break continues in the right-hand panel (clipped there)
            draw_mark(axc, inst)

    # Register quoted without a costed depth: an LQ rule, not an invented point.
    # A quoted register range (Ch. 10 stretch 2554-3054, R11) is a dotted band.
    by_lq = {}
    for i in LQ_ONLY:
        by_lq.setdefault((i.get("lq_lo", i["lq"]), i["lq"]), []).append(i)
    for (lq_lo, lq), group in sorted(by_lq.items()):
        chs = sorted({g["ch"] for g in group})
        fr = {FRONTIER[c] for c in chs}
        c = COLOR[FRONTIER[chs[0]]] if len(fr) == 1 else ff.DARKGREY
        for y in sorted({lq_lo, lq}):
            axc.axhline(y, color=c, lw=0.9, ls=":", zorder=3)
        if lq_lo < lq:
            axc.axhspan(lq_lo, lq, color=c, alpha=0.08, lw=0, zorder=2)
        axc.annotate(",".join(str(x) for x in chs) + "  register only",
                     xy=(1.15e12, lq_lo / 1.06), fontsize=5.5, color=c,
                     ha="left", va="top")

    for a, xlim in ((ax, XMAIN), (axc, XCODE)):
        a.set_xscale("log")
        a.set_yscale("log")
        a.set_xlim(*xlim)
        a.set_ylim(*YLIM)

    ax.set_ylabel("algorithmic logical qubits")
    ax.set_yticks([20, 50, 100, 200, 500, 1000, 2000])
    ax.set_yticklabels(["20", "50", "100", "200", "500", "1000", "2000"])
    axc.tick_params(axis="y", which="both", left=False, right=False)
    axc.set_xticks([1e13, 1e15])
    # The main panel ends at 1e12 (E19). A label on that last tick would sit on
    # the cut marks, so the labeled ticks stop at 1e11; the minor ticks run on.
    ax.set_xticks([1e4, 1e5, 1e6, 1e7, 1e8, 1e9, 1e10, 1e11])

    # Axis break: hide the facing spines. The cut marks are drawn in main()
    # after the layout is final, since figure coords are needed.
    ax.spines["right"].set_visible(False)
    axc.spines["left"].set_visible(False)
    ax.tick_params(right=False, which="both")

    axc.annotate("beyond RFI", xy=(0.5, 1.015), xycoords="axes fraction",
                 ha="center", va="bottom", fontsize=6.5, color=ff.DARKGREY)

    fig.supxlabel("T-gates (hard operations) per shot", fontsize=8, y=0.055)

    # Frontier key and geometry key, set as text rather than boxed legends.
    # Lower middle (3e5-1e8 T, below 100 LQ): the one band no instance reaches
    # once the main axis runs to 1e11 and the 2033 labels fill the top.
    for i, (name, col) in enumerate((("intensity (Chs. 2-4)", ff.BLUE),
                                     ("energy (5-6)", ff.ORANGE),
                                     ("cosmic (7-10)", ff.GREEN))):
        ax.annotate(name, xy=(0.245, 0.300 - 0.052 * i), xycoords="axes fraction",
                    fontsize=6.5, color=col, ha="left", va="top")

    # Geometry key, in the co-design panel's empty lower half. That panel
    # carries two marks and a rule, so it has room the main panel does not,
    # and putting the key there keeps it clear of the 2033 cluster entirely.
    # Swatches are drawn by the same code paths as the data marks, so the key
    # cannot describe a style the plot does not use.
    tr = axc.transAxes
    gx, gy, dyr = 0.92, 0.235, 0.062
    sx = gx - 0.70
    # Only styles the plot actually uses get a key row.
    used = {draw_geometry(i) for i in INSTANCES}
    rows = [("single value", "point"),
            ("range", "range"),
            ("upper bound", "bound"),
            ("conditional", "conditional")]
    have = {"point": "point" in used or "codesign" in used,
            "range": bool(used & {"capsule_t", "capsule_lq", "blob"}),
            "bound": "bound" in used,
            "conditional": "conditional" in used}
    rows = [r for r in rows if have[r[1]]]
    for i, (lab, kind) in enumerate(rows):
        y = gy - i * dyr
        axc.annotate(lab, xy=(gx, y), xycoords="axes fraction",
                     fontsize=5.5, color=ff.DARKGREY, ha="right", va="center")
        if kind == "point":
            axc.plot([sx], [y], "o", ms=MARK_MS, color=ff.DARKGREY, transform=tr,
                     mec="white", mew=0.5, clip_on=False, zorder=6)
        elif kind == "range":
            axc.plot([sx - 0.10, sx + 0.10], [y] * 2, lw=CAP_LW,
                     color=ff.DARKGREY, alpha=0.55, solid_capstyle="round",
                     transform=tr, clip_on=False, zorder=6)
        elif kind == "bound":
            axc.plot([sx + 0.08], [y], "o", ms=MARK_MS, color=ff.DARKGREY,
                     transform=tr, mec="white", mew=0.5, clip_on=False, zorder=6)
            axc.annotate("", xy=(sx - 0.09, y), xytext=(sx + 0.08, y),
                         xycoords="axes fraction", textcoords="axes fraction",
                         arrowprops=dict(arrowstyle="-|>", color=ff.DARKGREY, lw=0.8,
                                         mutation_scale=6, shrinkA=1.5, shrinkB=0))
        else:
            axc.add_patch(Rectangle((sx - 0.11, y - 0.022), 0.22, 0.044,
                                    transform=tr, facecolor=ff.DARKGREY, alpha=0.10,
                                    edgecolor=ff.DARKGREY, ls="--", lw=0.9,
                                    clip_on=False, zorder=6))

    # The 2028 rung is a genuine pile-up: twelve instances (Chs. 2-3, 5-10)
    # between 113 and 286 LQ at 5.4e3-4.4e5 T. Zoom it rather than let the labels collide.
    # Sits fully above the 2028 envelope's 250-LQ edge (y ~ 0.53 in axes
    # fraction) so the inset frame never touches the box it magnifies.
    # Width 0.31 keeps the frame left of Ch. 8's 2033 point (3.0e7 T, 1000 LQ).
    # Bottom 0.605 lifts the inset's x tick labels clear of the main panel's
    # 250-LQ envelope edge, and it sits just above the zoomed window's top
    # (300 LQ, axes fraction 0.56), so the inset hides none of the marks it
    # magnifies (E18; the old inset covered the 248- and 286-LQ points).
    axz = ax.inset_axes([0.085, 0.605, 0.31, 0.335])
    # Every 2028 instance lives at 5e3-1e5 T, so the T window is tight around
    # that; the register window runs past 250 so the envelope corner is
    # unambiguously inside the frame.
    # Upper T edge 1.1e6 so Ch. 5's 2028 row (5.94e5 T, E20) is inside the zoom with room for its label.
    zx, zy = (4.0e3, 1.1e6), (105, 300)   # 4.0e3: just left of the Ch. 3 mark (5.4e3, R13 one query), inside the main axes;
                                          # 1.1e6 (E19): Ch. 5 2028 at 5.9e5 is the largest 2028 mark (r19; Ch. 10 now 4.0e5, was 7.5e5)
    axz.add_patch(Rectangle((zx[0], zy[0]), 1e5 - zx[0], 250 - zy[0],
                            facecolor=ff.FNAL_BLUE, alpha=0.075,
                            edgecolor="none", zorder=1))
    axz.plot([zx[0], 1e5], [250, 250], color=ff.FNAL_BLUE, lw=0.8,
             ls=(0, (4, 2)), zorder=2)
    axz.plot([1e5, 1e5], [zy[0], 250], color=ff.FNAL_BLUE, lw=0.8,
             ls=(0, (4, 2)), zorder=2)
    zoomed = [i for i in INSTANCES if i["era"] == "2028" and i["ch"] != 4]
    for inst in zoomed:
        # The full stub would reach back past the Ch. 3 mark at 5.4e3 in the
        # tightened window; a short one still reads as an upper bound.
        draw_mark(axz, inst, bound_decades=0.16)
    axz.set_xscale("log")
    axz.set_yscale("log")
    axz.set_xlim(*zx)
    axz.set_ylim(*zy)
    # Log minor ticks would relabel these axes over the top of the fixed ticks.
    axz.minorticks_off()
    axz.set_xticks([1e4, 1e5, 1e6])
    axz.set_xticklabels([r"$10^4$", r"$10^5$", r"$10^6$"], fontsize=5.5)
    axz.set_yticks([120, 180, 250, 300])
    axz.set_yticklabels(["120", "180", "250", "300"], fontsize=5.5)
    axz.tick_params(length=2, width=0.5, pad=1.5)
    # The x tick labels hang below the frame, over the main panel's 2028
    # cluster; a white backing keeps them legible and hides nothing in the inset.
    axz.set_zorder(8)
    # The y tick labels sit left of the frame, where the main panel's 1000-LQ
    # envelope edge runs through them; they get the same backing.
    for lab in axz.get_xticklabels() + axz.get_yticklabels():
        lab.set_bbox(dict(facecolor="white", edgecolor="none", pad=0.4))
    for s in axz.spines.values():
        s.set_linewidth(0.5)
        s.set_color(ff.GREY)
    axz.annotate("2028 rung, zoomed", xy=(0.5, 1.02), xycoords="axes fraction",
                 ha="center", va="bottom", fontsize=6, color=ff.DARKGREY)
    ax.indicate_inset_zoom(axz, edgecolor=ff.GREY, lw=0.5, alpha=0.9)

    fig.subplots_adjust(left=0.085, right=0.985, top=0.93, bottom=0.135)

    # Labels last, once every axes has its final limits: placement is computed
    # in display space, and the same call serves the main panel, the co-design
    # panel and the inset, so all three obey one rule.
    in_inset = {id(i) for i in zoomed}
    draw_labels(ax, [i for i in INSTANCES
                     if not i["codesign"] and id(i) not in in_inset])
    draw_labels(axc, [i for i in INSTANCES if i["codesign"]])
    # In the inset the 1e5 envelope edge runs through the 2028 cluster; treating
    # it as an obstacle pushed Ch. 2's 1.004e5-T label beside the Ch. 7/8 marks.
    draw_labels(axz, zoomed, merge_pts=17.0, lines_block=False)

    # Cut marks last: figure coordinates are only final after the layout is.
    fig.canvas.draw()
    kw = dict(transform=fig.transFigure, color=ff.INK, lw=0.7, clip_on=False)
    xa, xb = ax.get_position().x1, axc.get_position().x0
    ya, yb = ax.get_position().y0, ax.get_position().y1
    for xx in (xa, xb):
        for yy in (ya, yb):
            fig.add_artist(plt.Line2D([xx - 0.005, xx + 0.005],
                                      [yy - 0.011, yy + 0.011], **kw))

    ff.save(fig, "resource_landscape")


if __name__ == "__main__":
    main()
