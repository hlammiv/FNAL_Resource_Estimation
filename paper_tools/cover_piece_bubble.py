#!/usr/bin/env python3
"""Reusable bubble-chamber window ornament (15-foot chamber photograph).

A dense tangle of particle tracks in an elliptical chamber window: a beam
ribbon sweeps in (two tracks die in hadronic stars), secondaries spray with
both curvature signs, a V0 decays after an invisible neutral gap, delta-ray
curls bud off the beams, stale tracks from earlier expansions cross faintly,
and e+e- pairs wind up as opposite decaying spirals -- one pair orange as
the single accent.  Pure linework with fine bubble dots and fiducial crosses.

Exposes draw_bubble_chamber(ax, center, rot_deg, scale, alpha, z): all
geometry lives in a fixed design frame and is pushed through
rotate(rot_deg) -> scale -> translate(center) at data level.  Deterministic.
"""

import numpy as np
from matplotlib.patches import Ellipse

# ------------------------------------------------------ self-contained palette
WHITE = "#f2f6fb"
PALE = "#9fc4e8"
ORANGE = "#e8833a"

# ------------------------------------------------- design-frame geometry
# The ornament was authored in absolute page coordinates centered here; the
# transform helper subtracts this center so callers only supply `center`.
_DCX, _DCY = 4.25, 8.52        # design center (window center)
_EA, _EB = 0.82, 0.615          # window semi-axes in the design frame

DESIGN_RADIUS = _EA             # max half-extent of the design (data units)


# ---------------------------------------------------------- track engine
def _trk(p0, ang, S, kfun, n=240):
    """Integrate a planar track: heading theta' = kappa(s)."""
    s = np.linspace(0.0, S, n)
    ds = S / (n - 1)
    kap = kfun(s)
    th = ang + np.concatenate(
        [[0.0], np.cumsum(0.5 * (kap[1:] + kap[:-1])) * ds])
    xs = p0[0] + np.concatenate(
        [[0.0], np.cumsum(0.5 * (np.cos(th[1:]) + np.cos(th[:-1]))) * ds])
    ys = p0[1] + np.concatenate(
        [[0.0], np.cumsum(0.5 * (np.sin(th[1:]) + np.sin(th[:-1]))) * ds])
    return xs, ys, th


def _arc(p0, ang, S, k0, drift=0.0, n=240):
    """Constant-curvature track with a slow tightening term (energy loss)."""
    return _trk(p0, ang, S, lambda s: k0 + drift * s, n)


def _spiral(p0, ang, k0, S, g=3.3, n=340):
    """Decaying spiral: curvature grows exponentially (radiating electron)."""
    return _trk(p0, ang, S, lambda s: k0 * np.exp(g * s), n)


def _at(xs, ys, th, frac):
    i = int(frac * (len(xs) - 1))
    return xs[i], ys[i], th[i], i


# --------------------------------------------------------------- ornament
def draw_bubble_chamber(ax, center, rot_deg=0.0, scale=1.0, alpha=1.0, z=1):
    """Bubble-chamber window ornament. center = data coords of the design
    center; rot_deg CCW; alpha multiplies every element's alpha; artists at
    zorder z..z+2."""
    rng = np.random.default_rng(19730404)      # fixed: deterministic art
    arts = []

    # data-level transform: rotate about design center, scale, translate
    _c = np.cos(np.radians(rot_deg))
    _s = np.sin(np.radians(rot_deg))
    cx0, cy0 = center

    def T(xs, ys):
        xs = np.asarray(xs, dtype=float) - _DCX
        ys = np.asarray(ys, dtype=float) - _DCY
        return (cx0 + scale * (_c * xs - _s * ys),
                cy0 + scale * (_s * xs + _c * ys))

    def A(a):                     # alpha
        return a * alpha

    def LW(w):                    # linewidths / marker sizes follow scale
        return w * scale

    def Z(zorig):                 # original zorders 5..6.8 -> z..z+1.8
        return z + (zorig - 5.0)

    # clip ellipse: transformed center + rotated angle, added to the axes so
    # its clip path lives in data coordinates
    clip = Ellipse(center, 2 * _EA * scale, 2 * _EB * scale, angle=rot_deg,
                   facecolor="none", edgecolor="none", zorder=Z(5))
    ax.add_patch(clip)

    def P(xs, ys, **kw):
        kw.setdefault("solid_capstyle", "round")
        kw["lw"] = LW(kw["lw"])
        kw["alpha"] = A(kw["alpha"])
        kw["zorder"] = Z(kw["zorder"])
        ln, = ax.plot(*T(xs, ys), **kw)
        arts.append(ln)
        return ln

    def bubbles(xs, ys, spacing=0.024, ms=0.85, color=None,
                alpha=0.55, z=6.1):
        """Fine bubble texture: dots resampled evenly along the track."""
        seg = np.hypot(np.diff(xs), np.diff(ys))
        cum = np.concatenate([[0.0], np.cumsum(seg)])
        tgt = np.arange(0.0, cum[-1], spacing)
        bx = np.interp(tgt, cum, xs)
        by = np.interp(tgt, cum, ys)
        jx = rng.normal(0, 0.0032, len(bx))
        jy = rng.normal(0, 0.0032, len(by))
        ln, = ax.plot(*T(bx + jx, by + jy), "o", ms=LW(ms), color=color,
                      alpha=A(alpha), zorder=Z(z), markeredgewidth=0)
        arts.append(ln)

    def pair(pv, ang, k0a, k0b, Sa, Sb, g, color, lw, alpha, z, open_a):
        """e+e- conversion: two opposite decaying spirals from one point."""
        xa, ya, _ = _spiral(pv, ang + open_a, k0a, Sa, g)
        xb, yb, _ = _spiral(pv, ang - open_a, -k0b, Sb, g)
        P(xa, ya, color=color, lw=lw, alpha=alpha, zorder=z)
        P(xb, yb, color=color, lw=lw, alpha=alpha, zorder=z)

    # ---- the window: faint photographic pane + double rim -------------
    ax.add_patch(Ellipse(center, 2 * _EA * scale, 2 * _EB * scale,
                         angle=rot_deg, facecolor=WHITE, alpha=A(0.05),
                         edgecolor="none", zorder=Z(5)))
    ax.add_patch(Ellipse(center, 2 * _EA * scale, 2 * _EB * scale,
                         angle=rot_deg, facecolor="none", edgecolor=PALE,
                         lw=LW(0.7), alpha=A(0.8), zorder=Z(5.1)))
    ax.add_patch(Ellipse(center, 2 * _EA * scale * 0.955,
                         2 * _EB * scale * 0.955, angle=rot_deg,
                         facecolor="none", edgecolor=PALE, lw=LW(0.3),
                         alpha=A(0.35), zorder=Z(5.1)))

    # fiducial crosses, as on the 15-foot chamber windows; the rotated
    # 4-spoke asterisk marker keeps the "+" aligned with the window
    fx = [3.76, 4.76, 3.68, 4.84, 4.28, 4.22]
    fy = [8.93, 8.90, 8.14, 8.17, 9.02, 8.06]
    ln, = ax.plot(*T(fx, fy), linestyle="none", marker=(4, 2, rot_deg),
                  ms=LW(3.2), mew=LW(0.5), color=PALE, alpha=A(0.38),
                  zorder=Z(5.2))
    arts.append(ln)

    # ---- stale tracks from earlier expansions: faint full chords ------
    strays = [((3.48, 7.98), 0.52, -0.28, 1.95),
              ((3.40, 8.78), -0.22, 0.22, 1.85),
              ((3.98, 7.90), 1.08, 0.30, 1.55),
              ((4.92, 8.02), 2.25, -0.26, 1.60)]
    for j, (p0, a0, k, S) in enumerate(strays):
        xs, ys, _ = _arc(p0, a0, S, k)
        P(xs, ys, color=PALE, lw=0.32, alpha=0.30, zorder=5.6)
        if j == 0:
            bubbles(xs, ys, spacing=0.030, ms=0.7, color=PALE, alpha=0.30,
                    z=5.6)

    # ---- beam ribbon: six tracks bunched, sweeping gently upward ------
    beams = []
    for i in range(6):
        y0 = 8.26 + 0.048 * i + rng.normal(0, 0.004)
        k = 0.21 + rng.normal(0, 0.02)
        a0 = 0.08 + rng.normal(0, 0.012)
        beams.append(_arc((3.30, y0), a0, 2.0, k, drift=0.02))

    vAx, vAy, vAth, iA = _at(*beams[2], 0.38)
    vBx, vBy, vBth, iB = _at(*beams[4], 0.70)
    for i, (xs, ys, th) in enumerate(beams):
        if i == 2:
            P(xs[:iA + 1], ys[:iA + 1], color=WHITE, lw=0.65, alpha=0.92,
              zorder=6.4)
        elif i == 4:
            P(xs[:iB + 1], ys[:iB + 1], color=WHITE, lw=0.65, alpha=0.92,
              zorder=6.4)
        else:
            P(xs, ys, color=WHITE, lw=0.6, alpha=0.85, zorder=6.4)
    bubbles(*beams[1][:2], spacing=0.027, ms=0.8, color=WHITE, alpha=0.5)

    # ---- primary star A: eight prongs, one kinked, one bubbled --------
    prongsA = [(-2.10, 1.4, 0.60), (-1.75, -1.2, 0.80), (-1.35, -2.2, 0.65),
               (-0.85, 1.6, 0.85), (-0.42, -0.9, 1.15), (-0.10, 0.50, 1.35),
               (0.18, -0.35, 1.30), (0.55, 1.9, 0.80), (0.95, -1.3, 0.85),
               (1.45, 2.8, 0.55)]
    endsA = []
    for (da, k0, S) in prongsA:
        xs, ys, th = _arc((vAx, vAy), vAth + da, S, k0,
                          drift=0.6 * np.sign(k0))
        lw = 0.5 if abs(k0) < 1.4 else 0.42
        P(xs, ys, color=PALE, lw=lw, alpha=0.72, zorder=6)
        endsA.append((xs, ys, th))
    bubbles(*endsA[6][:2], spacing=0.022, ms=0.8, color=PALE, alpha=0.6)
    # decay-in-flight kink on prong 8
    xs, ys, th = endsA[8]
    xs2, ys2, _ = _arc((xs[-1], ys[-1]), th[-1] + 0.60, 0.36, -2.0,
                       drift=-1.2)
    P(xs2, ys2, color=PALE, lw=0.4, alpha=0.65, zorder=6)

    # ---- V0: invisible neutral, then a two-prong V --------------------
    nang = vAth + 0.95
    dvx = vAx + 0.34 * np.cos(nang)
    dvy = vAy + 0.34 * np.sin(nang)
    P([vAx, dvx], [vAy, dvy], color=PALE, lw=0.35, ls=(0, (0.8, 1.8)),
      alpha=0.3, zorder=5.5)
    xs, ys, _ = _arc((dvx, dvy), nang + 0.32, 0.62, 1.0, drift=0.6)
    P(xs, ys, color=PALE, lw=0.45, alpha=0.75, zorder=6)
    xsv, ysv, _ = _arc((dvx, dvy), nang - 0.28, 0.55, -1.5, drift=-0.8)
    P(xsv, ysv, color=PALE, lw=0.45, alpha=0.75, zorder=6)
    bubbles(xsv, ysv, spacing=0.021, ms=0.75, color=PALE, alpha=0.6)

    # ---- secondary star B: five prongs, one kinked --------------------
    prongsB = [(-0.70, 2.0, 0.60), (-0.25, -0.70, 1.00),
               (0.10, 0.60, 0.90), (0.50, -1.8, 0.65), (0.95, 1.2, 0.55)]
    endsB = []
    for (da, k0, S) in prongsB:
        xs, ys, th = _arc((vBx, vBy), vBth + da, S, k0,
                          drift=0.55 * np.sign(k0))
        P(xs, ys, color=PALE, lw=0.45, alpha=0.7, zorder=6)
        endsB.append((xs, ys, th))
    bubbles(*endsB[1][:2], spacing=0.023, ms=0.75, color=PALE, alpha=0.55)
    xs, ys, th = endsB[3]
    xs2, ys2, _ = _arc((xs[-1], ys[-1]), th[-1] - 0.70, 0.28, 2.8, drift=1.6)
    P(xs2, ys2, color=PALE, lw=0.38, alpha=0.6, zorder=6)

    # ---- THE accent: a large orange e+e- pair off star A --------------
    gang = vAth - 0.62
    gx = vAx + 0.58 * np.cos(gang)
    gy = vAy + 0.58 * np.sin(gang)
    P([vAx, gx], [vAy, gy], color=PALE, lw=0.35, ls=(0, (0.8, 1.8)),
      alpha=0.3, zorder=5.5)
    pair((gx, gy), gang, 4.2, 5.2, 0.85, 0.75, 3.3, ORANGE, 0.6, 0.95,
         6.8, open_a=0.15)

    # ---- supporting pale pairs and lone electron spirals --------------
    pair((3.85, 8.86), 0.90, 8.0, 10.0, 0.54, 0.47, 3.5, PALE, 0.45,
         0.60, 6.2, open_a=0.18)
    pair((4.88, 8.58), -2.60, 8.0, 10.0, 0.54, 0.47, 3.5, PALE, 0.45,
         0.60, 6.2, open_a=0.18)
    pair((3.84, 8.14), 0.55, 9.0, 11.0, 0.48, 0.42, 3.5, PALE, 0.42,
         0.55, 6.2, open_a=0.20)
    xs, ys, _ = _spiral((3.64, 8.28), 2.2, 8.0, 0.55, 3.4)
    P(xs, ys, color=PALE, lw=0.4, alpha=0.55, zorder=6.2)
    xs, ys, _ = _spiral((4.05, 8.97), -2.4, -9.0, 0.52, 3.4)
    P(xs, ys, color=PALE, lw=0.4, alpha=0.55, zorder=6.2)
    xs, ys, _ = _spiral((4.66, 8.97), 2.9, 10.0, 0.48, 3.4)
    P(xs, ys, color=PALE, lw=0.4, alpha=0.55, zorder=6.2)
    xs, ys, _ = _spiral((4.52, 8.08), 2.6, 11.0, 0.44, 3.4)
    P(xs, ys, color=PALE, lw=0.4, alpha=0.55, zorder=6.2)

    # ---- delta rays: tight curls budding off two beams ----------------
    dx, dy, dth, _ = _at(*beams[0], 0.24)
    xs, ys, _ = _spiral((dx, dy), dth + 1.35, 12.0, 0.42, 3.2)
    P(xs, ys, color=PALE, lw=0.35, alpha=0.6, zorder=6)
    dx, dy, dth, _ = _at(*beams[5], 0.72)
    xs, ys, _ = _spiral((dx, dy), dth - 1.25, -11.0, 0.44, 3.2)
    P(xs, ys, color=PALE, lw=0.35, alpha=0.6, zorder=6)

    # ---- background boiling: sparse isolated bubbles ------------------
    u = rng.uniform(0, 1, 14)
    aa = rng.uniform(0, 2 * np.pi, 14)
    bx = _DCX + 0.94 * _EA * np.sqrt(u) * np.cos(aa)
    by = _DCY + 0.94 * _EB * np.sqrt(u) * np.sin(aa)
    ln, = ax.plot(*T(bx, by), "o", ms=LW(0.9), color=PALE, alpha=A(0.35),
                  zorder=Z(5.4), markeredgewidth=0, linestyle="none")
    arts.append(ln)

    for a in arts:
        a.set_clip_path(clip)


# ------------------------------------------------------------- test harness
if __name__ == "__main__":
    import pathlib

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    HERE = pathlib.Path(__file__).resolve().parent
    fig = plt.figure(figsize=(8.5, 11))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 8.5)
    ax.set_ylim(0, 11)
    ax.set_axis_off()

    # plate + page frame (the corner piece should read as cropped by it)
    ax.add_patch(Rectangle((0, 0), 8.5, 11, facecolor="#0d2c56", zorder=0))
    ax.add_patch(Rectangle((0.45, 0.45), 7.6, 10.1, fill=False,
                           edgecolor="#9fc4e8", lw=1.0, alpha=0.7,
                           zorder=0.5))

    draw_bubble_chamber(ax, (2.0, 8.0), rot_deg=0.0, scale=1.0)
    draw_bubble_chamber(ax, (7.8, 10.2), rot_deg=30.0, scale=0.85,
                        alpha=0.6)
    draw_bubble_chamber(ax, (2.0, 3.0), rot_deg=-20.0, scale=0.6)

    fig.savefig(HERE / "test.png", dpi=150)
    print("wrote", HERE / "test.png")
