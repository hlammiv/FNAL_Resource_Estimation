#!/usr/bin/env python3
"""Reusable collider event-display ornament (transverse view).

Dashed beam axis through a primary vertex, ~35 charged tracks with
charge-dependent curvature, three collimated jets, a segmented barrel
calorimeter with radial energy-histogram bars (leading towers orange),
two muons with chamber stubs beyond the outer ring, a displaced
secondary vertex, tracker hit ticks, and a dashed missing-ET arrow.

Everything is drawn in a local design frame and pushed through a
data-level rotate -> scale -> translate map, so the piece can be placed
anywhere (e.g. the cover sky), rotated, scaled, and faded.  All artists
are polylines / Polygons in data coordinates -- rotation is rigid by
construction.  Deterministic (fixed rng seed per call).
"""

import numpy as np
from matplotlib.patches import Polygon

WHITE = "#f2f6fb"
PALE = "#9fc4e8"
ORANGE = "#e8833a"
D2R = np.pi / 180.0

# ----------------------------------------------------- detector radii -----
R_PIPE = 0.030
LAYERS = (0.085, 0.150, 0.215, 0.285)   # silicon tracker layers
R_TRK = 0.355                            # tracker outer boundary
R_CAL0, R_CAL1 = 0.400, 0.455            # segmented ECAL ring
R_BAR = 0.475                            # calo histogram baseline
R_OUT = 0.600                            # outer (muon-system) boundary
SEG_STEP = 7.5                           # calo segmentation, degrees

PHI_A, PHI_B, PHI_C = 150.0, 38.0, 266.0   # jet axes (A = leading)
PHI_MET = -26.0                            # recoils against jet A
PHI_MU1, PHI_MU2 = 97.0, 207.0

# max extent from the design center: the dashed beam axis reaches +-0.85
# (detector rings / muon stubs stop at 0.685).
DESIGN_RADIUS = 0.85


def _helix(phi0_deg, q, R, r_stop, start=(0.0, 0.0), n=260):
    """Circular-arc track of signed curvature q/R from `start`, tangent
    direction phi0; truncated where the radius from the VERTEX reaches
    r_stop (loopers simply exhaust their turning angle)."""
    phi0 = phi0_deg * D2R
    th_cap = 0.4
    dx = dy = None
    for _ in range(26):
        th = np.linspace(0.0, th_cap, n)
        dx = start[0] + q * R * (np.sin(phi0 + q * th) - np.sin(phi0))
        dy = start[1] - q * R * (np.cos(phi0 + q * th) - np.cos(phi0))
        r = np.hypot(dx, dy)
        hit = np.nonzero(r >= r_stop)[0]
        if hit.size:
            k = hit[0] + 1
            return dx[:k], dy[:k]
        if th_cap > 3.6 * np.pi:
            break
        th_cap *= 1.6
    return dx, dy


def _spiral(phi0_deg, q, R0, th_max, shrink=0.35, n=300):
    """Low-pT curler: arc whose radius decays (dE/dx), so it spirals."""
    phi0 = phi0_deg * D2R
    th = np.linspace(0.0, th_max, n)
    Rt = R0 * (1.0 - shrink * th / th_max)
    ds = Rt * np.gradient(th)
    ang = phi0 + q * th
    return np.cumsum(np.cos(ang) * ds), np.cumsum(np.sin(ang) * ds)


def draw_event_display(ax, center, rot_deg=0.0, scale=1.0, alpha=1.0, z=1):
    """Collider event display ornament. center = data coords of the design
    center; rot_deg CCW; alpha multiplies every element's alpha; artists at
    zorder z..z+2."""
    rot = rot_deg * D2R
    cr, sr = np.cos(rot), np.sin(rot)
    cx, cy = center
    A = float(alpha)

    def xf(xs, ys):
        """Design coords -> data coords: rotate, scale, translate."""
        xs = np.asarray(xs, dtype=float)
        ys = np.asarray(ys, dtype=float)
        return (cx + scale * (cr * xs - sr * ys),
                cy + scale * (sr * xs + cr * ys))

    def pl(xs, ys, alpha=1.0, zorder=5, **kw):
        X, Y = xf(xs, ys)
        ax.plot(X, Y, alpha=min(1.0, alpha * A), zorder=z + (zorder - 5),
                **kw)

    def poly(vx, vy, facecolor, alpha, zorder):
        X, Y = xf(vx, vy)
        ax.add_patch(Polygon(np.column_stack((X, Y)), closed=True,
                             facecolor=facecolor, edgecolor="none",
                             alpha=min(1.0, alpha * A),
                             zorder=z + (zorder - 5)))

    def radial(phi_deg, r0, r1, **kw):
        a = phi_deg * D2R
        pl([r0 * np.cos(a), r1 * np.cos(a)],
           [r0 * np.sin(a), r1 * np.sin(a)], **kw)

    def hit_ticks(xs, ys, radii, color=WHITE, alpha=0.85, t=0.013):
        r = np.hypot(xs, ys)
        for R0 in radii:
            hit = np.nonzero(r >= R0)[0]
            if not hit.size or hit[0] == 0:
                continue
            i = hit[0]
            tx, ty = xs[i] - xs[i - 1], ys[i] - ys[i - 1]
            L = np.hypot(tx, ty)
            if L == 0:
                continue
            nx, ny = -ty / L, tx / L
            pl([xs[i] - nx * t, xs[i] + nx * t],
               [ys[i] - ny * t, ys[i] + ny * t],
               color=color, lw=0.4, alpha=alpha, zorder=6,
               solid_capstyle="butt")

    def calo_bar(phi_deg, h, wdeg, color, lw, alpha, zb=6, fill=False):
        a0, a1 = (phi_deg - wdeg / 2) * D2R, (phi_deg + wdeg / 2) * D2R
        r1 = R_BAR + h
        arc = np.linspace(a0, a1, 9)
        xs = np.concatenate(([R_BAR * np.cos(a0)], r1 * np.cos(arc),
                             [R_BAR * np.cos(a1)]))
        ys = np.concatenate(([R_BAR * np.sin(a0)], r1 * np.sin(arc),
                             [R_BAR * np.sin(a1)]))
        pl(xs, ys, color=color, lw=lw, alpha=alpha, zorder=zb,
           solid_capstyle="butt", solid_joinstyle="miter")
        if fill:                             # sampling rungs inside the tower
            rr = R_BAR + 0.05
            while rr < r1 - 0.018:
                pl(rr * np.cos(arc), rr * np.sin(arc), color=color, lw=0.35,
                   alpha=alpha * 0.6, zorder=zb, solid_capstyle="butt")
                rr += 0.05
        if fill:
            inner = np.linspace(a1, a0, 9)
            vx = np.concatenate((xs, R_BAR * np.cos(inner)))
            vy = np.concatenate((ys, R_BAR * np.sin(inner)))
            poly(vx, vy, color, 0.12, zb - 1)

    # ------------------------------------------------ beam axis -----------
    pl([-0.85, -R_PIPE - 0.01], [0, 0], color=PALE, lw=0.4, alpha=0.55,
       ls=(0, (2.6, 1.9)), zorder=5)
    pl([R_PIPE + 0.01, 0.85], [0, 0], color=PALE, lw=0.4, alpha=0.55,
       ls=(0, (2.6, 1.9)), zorder=5)
    for sgn in (-1, 1):                       # inward beam chevrons
        x0 = sgn * 0.80
        for dy_ in (0.020, -0.020):
            pl([x0, x0 - sgn * 0.035], [dy_, 0.0], color=PALE, lw=0.45,
               alpha=0.7, zorder=5)
    th = np.linspace(0, 2 * np.pi, 60)        # beam pipe
    pl(R_PIPE * np.cos(th), R_PIPE * np.sin(th), color=WHITE, lw=0.4,
       alpha=0.8, zorder=6)

    # ------------------------------------------------ tracker layers ------
    for i, r in enumerate(LAYERS):
        pl(r * np.cos(th), r * np.sin(th), color=PALE, lw=0.3,
           alpha=0.26 + 0.02 * i, zorder=5)
    pl(R_TRK * np.cos(th), R_TRK * np.sin(th), color=PALE, lw=0.35,
       alpha=0.38, zorder=5)

    # ------------------------------------------------ calorimeter ring ----
    thf = np.linspace(0, 2 * np.pi, 200)
    pl(R_CAL0 * np.cos(thf), R_CAL0 * np.sin(thf), color=PALE, lw=0.45,
       alpha=0.60, zorder=5)
    pl(R_CAL1 * np.cos(thf), R_CAL1 * np.sin(thf), color=PALE, lw=0.35,
       alpha=0.45, zorder=5)
    for adeg in np.arange(0.0, 360.0, SEG_STEP):
        radial(adeg, R_CAL0, R_CAL1, color=PALE, lw=0.3, alpha=0.30,
               zorder=5)

    # outer boundary (muon-system inner edge), dashed
    pl(R_OUT * np.cos(thf), R_OUT * np.sin(thf), color=PALE, lw=0.4,
       alpha=0.40, ls=(0, (3.5, 2.5)), zorder=5)
    for adeg in (0.0, 180.0):                # fiducial ticks on the beam line
        radial(adeg, 0.615, 0.655, color=PALE, lw=0.35, alpha=0.55, zorder=5)

    # ------------------------------------------------ jet cones -----------
    for phi_j, half, aa in ((PHI_A, 11, 0.08), (PHI_B, 10, 0.07),
                            (PHI_C, 12, 0.05)):
        arc = np.linspace((phi_j - half) * D2R, (phi_j + half) * D2R, 12)
        vx = np.concatenate(([0.0], R_CAL0 * np.cos(arc)))
        vy = np.concatenate(([0.0], R_CAL0 * np.sin(arc)))
        poly(vx, vy, PALE, aa, 5)

    # ------------------------------------------------ jet tracks ----------
    rng = np.random.default_rng(20260806)
    tick_layers = LAYERS + (R_TRK,)

    # jet A (leading): 4 prompt tracks ...
    for k, (dphi, q, R) in enumerate((( -7.5, +1, 3.2), (-2.0, -1, 4.5),
                                      (  3.5, +1, 2.2), ( 8.0, -1, 1.5))):
        xs, ys = _helix(PHI_A + dphi, q, R, R_CAL0)
        pl(xs, ys, color=WHITE, lw=0.6, alpha=0.95, zorder=6,
           solid_capstyle="round")
        if k in (0, 2):
            hit_ticks(xs, ys, tick_layers)
    # ... plus a secondary vertex at r = 0.105 with 3 displaced tracks
    sv_a = 147.0 * D2R
    sv = (0.105 * np.cos(sv_a), 0.105 * np.sin(sv_a))
    pl([0, sv[0]], [0, sv[1]], color=PALE, lw=0.35, alpha=0.7,
       ls=(0, (1.2, 1.2)), zorder=6)
    pl([sv[0]], [sv[1]], marker="o", ms=2.0 * scale, mfc="none", mec=WHITE,
       mew=0.5, alpha=0.95, zorder=7)
    for dphi, q, R in ((-14.0, -1, 1.1), (2.0, +1, 1.8), (16.0, -1, 0.85)):
        xs, ys = _helix(PHI_A + dphi, q, R, R_CAL0, start=sv)
        pl(xs, ys, color=WHITE, lw=0.5, alpha=0.85, zorder=6,
           solid_capstyle="round")

    # jet B: 5 tracks
    for k, (dphi, q, R) in enumerate(((-8.0, -1, 1.8), (-3.5, +1, 3.5),
                                      (0.5, +1, 1.2), (4.5, -1, 2.6),
                                      (9.0, +1, 0.9))):
        xs, ys = _helix(PHI_B + dphi, q, R, R_CAL0)
        pl(xs, ys, color=WHITE, lw=0.55, alpha=0.92, zorder=6,
           solid_capstyle="round")
        if k == 1:
            hit_ticks(xs, ys, tick_layers)

    # jet C (soft): 4 tracks
    for dphi, q, R in ((-10.0, +1, 0.55), (-3.0, -1, 0.9),
                       (4.0, +1, 0.7), (11.0, -1, 0.45)):
        xs, ys = _helix(PHI_C + dphi, q, R, R_CAL0)
        pl(xs, ys, color=WHITE, lw=0.5, alpha=0.8, zorder=6,
           solid_capstyle="round")

    # ------------------------------------------------ muons ---------------
    for phi_m, q in ((PHI_MU1, -1), (PHI_MU2, +1)):
        xs, ys = _helix(phi_m, q, 5.0, R_CAL0)
        pl(xs, ys, color=WHITE, lw=0.65, alpha=0.95, zorder=7,
           solid_capstyle="round")
        hit_ticks(xs, ys, tick_layers, alpha=0.9)
        a = (phi_m + q * 2.5) * D2R          # slight bend continues
        pl([R_CAL0 * np.cos(a), 0.585 * np.cos(a)],
           [R_CAL0 * np.sin(a), 0.585 * np.sin(a)],
           color=PALE, lw=0.45, alpha=0.6, zorder=6)
        for r0, r1 in ((0.615, 0.645), (0.657, 0.685)):
            pl([r0 * np.cos(a), r1 * np.cos(a)],
               [r0 * np.sin(a), r1 * np.sin(a)],
               color=WHITE, lw=0.7, alpha=0.95, zorder=7,
               solid_capstyle="butt")
            # muon-chamber outline: little annular box around the stub
            adeg = a / D2R
            arc0 = np.linspace((adeg - 5.0) * D2R, (adeg + 5.0) * D2R, 7)
            box_x = np.concatenate((r0 * np.cos(arc0), r1 * np.cos(arc0[::-1]),
                                    [r0 * np.cos(arc0[0])]))
            box_y = np.concatenate((r0 * np.sin(arc0), r1 * np.sin(arc0[::-1]),
                                    [r0 * np.sin(arc0[0])]))
            pl(box_x, box_y, color=PALE, lw=0.35, alpha=0.55, zorder=6,
               solid_capstyle="butt")

    # ------------------------------------------------ soft tracks ---------
    veto = (PHI_A, PHI_B, PHI_C, PHI_MU1, PHI_MU2, PHI_MET % 360.0)
    placed = 0
    while placed < 11:
        phi = rng.uniform(0.0, 360.0)
        if any(min(abs(phi - v), 360 - abs(phi - v)) < 18.0 for v in veto):
            continue
        q = int(rng.choice((-1, 1)))
        R = float(np.exp(rng.uniform(np.log(0.20), np.log(1.5))))
        r_stop = float(rng.uniform(0.30, 0.40))
        xs, ys = _helix(phi, q, R, r_stop)
        pl(xs, ys, color=PALE, lw=0.45, alpha=0.58, zorder=6,
           solid_capstyle="round")
        if placed == 4:
            hit_ticks(xs, ys, LAYERS[:3], color=PALE, alpha=0.7)
        placed += 1

    # tight curlers (very low pT, spiralling in the tracker)
    for phi, q, R0, tm in ((63.0, +1, 0.095, 3.2 * np.pi),
                           (232.0, -1, 0.080, 3.5 * np.pi),
                           (299.0, +1, 0.105, 2.9 * np.pi)):
        xs, ys = _spiral(phi, q, R0, tm)
        pl(xs, ys, color=PALE, lw=0.45, alpha=0.62, zorder=6,
           solid_capstyle="round")

    # neutral particles: dashed to the calo face, then a deposit
    for phi in (118.0, 295.0):
        radial(phi, 0.045, R_CAL0, color=PALE, lw=0.4, alpha=0.55,
               ls=(0, (1.6, 1.6)), zorder=6)
        calo_bar(phi, 0.045, 4.0, PALE, 0.5, 0.6)

    # ------------------------------------------------ calo histogram ------
    # leading-jet towers in orange, sidebands pale
    calo_bar(PHI_A + 1.5, 0.160, 5.5, ORANGE, 0.8, 0.95, zb=7, fill=True)
    calo_bar(PHI_A - 4.5, 0.115, 5.0, ORANGE, 0.75, 0.9, zb=7, fill=True)
    for dphi, h in ((-10.5, 0.055), (7.5, 0.075), (13.0, 0.035)):
        calo_bar(PHI_A + dphi, h, 4.5, PALE, 0.55, 0.75)
    calo_bar(PHI_B - 1.0, 0.125, 5.5, ORANGE, 0.75, 0.9, zb=7, fill=True)
    for dphi, h in ((-7.0, 0.065), (5.0, 0.085), (10.5, 0.040)):
        calo_bar(PHI_B + dphi, h, 4.5, PALE, 0.55, 0.75)
    for dphi, h in ((-6.0, 0.050), (1.0, 0.070), (7.0, 0.040)):
        calo_bar(PHI_C + dphi, h, 4.5, PALE, 0.55, 0.75)
    # underlying-event activity around the ring
    for _ in range(16):
        phi = float(rng.uniform(0.0, 360.0))
        if any(min(abs(phi - v), 360 - abs(phi - v)) < 16.0
               for v in (PHI_A, PHI_B, PHI_C)):
            continue
        radial(phi, R_BAR, R_BAR + float(rng.uniform(0.015, 0.048)),
               color=PALE, lw=0.5, alpha=0.45, zorder=6,
               solid_capstyle="butt")

    # ------------------------------------------------ missing ET ----------
    a = PHI_MET * D2R
    ca, sa = np.cos(a), np.sin(a)
    pl([0.055 * ca, 0.60 * ca], [0.055 * sa, 0.60 * sa], color=ORANGE,
       lw=0.7, alpha=0.95, ls=(0, (3.0, 2.0)), zorder=7,
       solid_capstyle="butt")
    for chev in (+152.0, -152.0):            # arrowhead chevron
        b = (PHI_MET + chev) * D2R
        pl([0.60 * ca, 0.60 * ca + 0.045 * np.cos(b)],
           [0.60 * sa, 0.60 * sa + 0.045 * np.sin(b)],
           color=ORANGE, lw=0.7, alpha=0.95, zorder=7)

    # ------------------------------------------------ primary vertex ------
    pl([0], [0], marker="o", ms=1.8 * scale, mfc=WHITE, mec="none",
       alpha=1.0, zorder=7)


if __name__ == "__main__":
    import pathlib

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    DIR = pathlib.Path(__file__).resolve().parent
    PLATE = "#0d2c56"

    fig = plt.figure(figsize=(8.5, 11))
    fig.patch.set_facecolor(PLATE)
    ax = fig.add_axes((0, 0, 1, 1))
    ax.set_facecolor(PLATE)
    ax.set_xlim(0, 8.5)
    ax.set_ylim(0, 11)
    ax.set_aspect("equal")
    ax.axis("off")

    draw_event_display(ax, (2.0, 8.0), rot_deg=0.0, scale=1.0)
    draw_event_display(ax, (7.0, 6.0), rot_deg=-15.0, scale=0.6, alpha=0.8)
    draw_event_display(ax, (1.0, 3.0), rot_deg=25.0, scale=0.5)

    out = DIR / "test.png"
    fig.savefig(out, dpi=150, facecolor=fig.get_facecolor())
    print(f"wrote {out}")
