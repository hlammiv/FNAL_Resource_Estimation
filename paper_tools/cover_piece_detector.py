#!/usr/bin/env python3
"""Reusable detector-slice ornament: particle-detector cross-section wedge.

A CMS/ATLAS-style transverse slice: beam pipe (beam out of page), silicon
tracker layers, ECAL with crystal segmentation, HCAL with tower segmentation
and absorber laminations, solenoid winding cross-section, steel return yoke
with hatching, and two staggered muon stations drawn as segmented drift-tube
boxes.  Traversed by six particle signatures: a low-p_T looper, an electron
showering in the ECAL, a dashed photon showering in the ECAL, a charged
hadron showering in the HCAL, an orange muon punching through to the
chambers, and a faint dotted neutrino escaping.  Pure linework in the
blueprint idiom; fills at alpha <= 0.13.

All geometry is defined in *design coordinates* relative to the interaction
point (the piece's design center) and pushed through a rigid
rotate -> scale -> translate transform at draw time, so the ornament can be
placed anywhere, at any angle, at any size, and faded with a global alpha.
The tether / plumb line that tied the original variant to the building has
been removed: the ornament is self-contained.
"""

import numpy as np
from matplotlib.patches import Circle, Polygon, Wedge

WHITE = "#f2f6fb"
PALE = "#9fc4e8"
ORANGE = "#e8833a"

d2r = np.deg2rad

# ------------------------------------------------------------ geometry ----
# Design center = interaction point = (0, 0) in design coordinates.
TH0, TH1 = 50.0, 130.0     # wedge angular span (degrees, in design frame)

R_PIPE = 0.045
R_PIX = (0.085, 0.105, 0.125)          # pixel layers, tightly spaced
R_STR = (0.18, 0.23, 0.285, 0.34)      # strip layers, module-segmented
R_TRK = (0.105, 0.18, 0.23, 0.285, 0.34)   # layers that get hit marks
R_EC = (0.40, 0.54)        # ECAL band
R_HC = (0.60, 0.82)        # HCAL band
R_SOL = (0.86, 0.94)       # solenoid cryostat band
R_YOKE = (1.10, 1.16)      # return-yoke band
R_MU1 = (1.00, 1.10)       # muon station 1
R_MU2 = (1.16, 1.26)       # muon station 2
R_TOP = 1.26               # outermost detector radius

MU1_SEG = ((51, 76), (79, 104), (107, 129))
MU2_SEG = ((51, 64), (67, 90), (93, 116), (119, 129))

# Approx max extent from the design center, in design units (muon arrowhead
# and the dash-dot center line both reach ~1.30-1.32).
DESIGN_RADIUS = 1.32


def draw_detector_slice(ax, center, rot_deg=0.0, scale=1.0, alpha=1.0, z=1):
    """Detailed detector-slice ornament. center = data coords of the piece's
    design center; rot_deg CCW; alpha multiplies every element's alpha;
    artists drawn at zorder z..z+2."""
    W, P, O = WHITE, PALE, ORANGE
    cx0, cy0 = center
    cr, sr = np.cos(d2r(rot_deg)), np.sin(d2r(rot_deg))

    # -- rigid data-level transform: rotate -> scale -> translate ----------
    def T(x, y):
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        return ((x * cr - y * sr) * scale + cx0,
                (x * sr + y * cr) * scale + cy0)

    def zm(z0):
        """Map the design zorders 5..7 onto z..z+2."""
        return z + (z0 - 5)

    def plot(x, y, *fmt, **kw):
        """ax.plot in design coords: transform, fade, remap zorder, and
        scale stroke widths / marker sizes with the piece."""
        kw['alpha'] = kw.get('alpha', 1.0) * alpha
        kw['zorder'] = zm(kw.get('zorder', 6))
        if 'lw' in kw:
            kw['lw'] = kw['lw'] * scale
        if 'ms' in kw:
            kw['ms'] = kw['ms'] * scale
        xt, yt = T(x, y)
        ax.plot(xt, yt, *fmt, **kw)

    def circle(c, r, **kw):
        """Circle patch in design coords (rotation-invariant)."""
        kw['alpha'] = kw.get('alpha', 1.0) * alpha
        kw['zorder'] = zm(kw.get('zorder', 6))
        if 'lw' in kw:
            kw['lw'] = kw['lw'] * scale
        xc, yc = T(*c)
        ax.add_patch(Circle((float(xc), float(yc)), r * scale, **kw))

    def wedge_band(r_out, a0, a1, width, **kw):
        """Annular wedge centered on the design origin; rotation folds into
        the angular span."""
        kw['alpha'] = kw.get('alpha', 1.0) * alpha
        kw['zorder'] = zm(kw.get('zorder', 5))
        xc, yc = T(0.0, 0.0)
        ax.add_patch(Wedge((float(xc), float(yc)), r_out * scale,
                           a0 + rot_deg, a1 + rot_deg, width=width * scale,
                           **kw))

    def polygon(pts, **kw):
        """Polygon patch: transform every vertex."""
        kw['alpha'] = kw.get('alpha', 1.0) * alpha
        kw['zorder'] = zm(kw.get('zorder', 5))
        if 'lw' in kw:
            kw['lw'] = kw['lw'] * scale
        xs, ys = T([p[0] for p in pts], [p[1] for p in pts])
        ax.add_patch(Polygon(np.column_stack([xs, ys]), closed=True, **kw))

    # -- design-coordinate helpers (identical to the source variant, with
    #    the interaction point at the origin) ------------------------------
    def pol(r, a_deg):
        a = d2r(a_deg)
        return r * np.cos(a), r * np.sin(a)

    def arc(r, a0=TH0, a1=TH1, color=P, lw=0.35, alpha=0.8, ls='-', z=6,
            n=90):
        th = d2r(np.linspace(a0, a1, n))
        plot(r * np.cos(th), r * np.sin(th), color=color, lw=lw,
             alpha=alpha, ls=ls, zorder=z, solid_capstyle='round')

    def radial(a, r0, r1, color=P, lw=0.3, alpha=0.5, ls='-', z=6):
        (xa, ya), (xb, yb) = pol(r0, a), pol(r1, a)
        plot([xa, xb], [ya, yb], color=color, lw=lw, alpha=alpha, ls=ls,
             zorder=z, solid_capstyle='round')

    # --------------------------------------------------------- beam pipe --
    circle((0.0, 0.0), R_PIPE, facecolor='none', edgecolor=W, lw=0.5,
           zorder=6)
    plot([0.0], [0.0], 'o', ms=1.5, color=O, zorder=7)   # beam out of page

    # ------------------------------------------------------ faint fills ---
    wedge_band(R_EC[1], TH0, TH1, R_EC[1] - R_EC[0], facecolor=W,
               edgecolor='none', alpha=0.06, zorder=5)
    wedge_band(R_HC[1], TH0, TH1, R_HC[1] - R_HC[0], facecolor=W,
               edgecolor='none', alpha=0.04, zorder=5)
    wedge_band(R_YOKE[1], TH0, TH1, R_YOKE[1] - R_YOKE[0], facecolor=W,
               edgecolor='none', alpha=0.05, zorder=5)

    # -------------------------------------------------------- tracker -----
    for r in R_PIX:                                     # pixel detector
        arc(r, color=P, lw=0.3, alpha=0.7)
    for r in R_STR:                                     # strip detector
        arc(r, color=P, lw=0.3, alpha=0.8)
        for a in np.arange(53.5, 128, 7):               # module segmentation
            (xa, ya), (xb, yb) = pol(r - 0.007, a), pol(r + 0.007, a)
            plot([xa, xb], [ya, yb], color=P, lw=0.3, alpha=0.55, zorder=6)
    # dashed continuation of the wedge borders through the tracker
    for a in (TH0, TH1):
        radial(a, 0.065, 0.37, color=P, lw=0.3, alpha=0.3,
               ls=(0, (1.2, 1.8)))

    # -------------------------------------------------------- preshower ---
    arc(0.365, color=P, lw=0.3, alpha=0.6)
    arc(0.378, color=P, lw=0.3, alpha=0.6)

    # ------------------------------------------------------------ ECAL ----
    arc(R_EC[0], color=W, lw=0.5, alpha=0.95)
    arc(R_EC[1], color=W, lw=0.5, alpha=0.95)
    for a in np.arange(55, 126, 5):                    # crystal segmentation
        radial(a, R_EC[0], R_EC[1], color=P, lw=0.3, alpha=0.45)
    # construction overshoot of the outer arc past the wedge borders
    arc(R_EC[1], 43.5, 49.0, color=P, lw=0.3, alpha=0.28)
    arc(R_EC[1], 131.0, 136.5, color=P, lw=0.3, alpha=0.28)

    # ------------------------------------------------------------ HCAL ----
    arc(R_HC[0], color=W, lw=0.5, alpha=0.95)
    arc(R_HC[1], color=W, lw=0.5, alpha=0.95)
    for a in np.arange(60, 121, 10):                    # tower segmentation
        radial(a, R_HC[0], R_HC[1], color=P, lw=0.3, alpha=0.45)
    for r in (0.672, 0.745):                            # absorber laminations
        arc(r, 51.5, 128.5, color=P, lw=0.3, alpha=0.35, ls=(0, (2.6, 2.0)))
    arc(R_HC[1], 43.5, 49.0, color=P, lw=0.3, alpha=0.28)
    arc(R_HC[1], 131.0, 136.5, color=P, lw=0.3, alpha=0.28)

    # -------------------------------------------------------- solenoid ----
    arc(R_SOL[0], color=W, lw=0.5, alpha=0.95)
    arc(R_SOL[1], color=W, lw=0.5, alpha=0.95)
    for a in np.arange(54, 127, 4):                     # winding cross-section
        circle(pol(0.90, a), 0.024, facecolor='none', edgecolor=P, lw=0.3,
               alpha=0.8, zorder=6)

    # ------------------------------------------------------------ yoke ----
    arc(R_YOKE[0], color=P, lw=0.35, alpha=0.6)
    arc(R_YOKE[1], color=P, lw=0.35, alpha=0.6)
    for a in np.arange(52, 127, 5):                     # steel hatching
        (xa, ya) = pol(R_YOKE[0], a)
        (xb, yb) = pol(R_YOKE[1], a + 2.4)
        plot([xa, xb], [ya, yb], color=P, lw=0.3, alpha=0.45, zorder=6)

    # --------------------------------------------------- muon stations ----
    for (rin, rout), segs in ((R_MU1, MU1_SEG), (R_MU2, MU2_SEG)):
        for (a0, a1) in segs:
            arc(rin, a0, a1, color=W, lw=0.45, alpha=0.95)
            arc(rout, a0, a1, color=W, lw=0.45, alpha=0.95)
            radial(a0, rin, rout, color=W, lw=0.45, alpha=0.95)
            radial(a1, rin, rout, color=W, lw=0.45, alpha=0.95)
            for f in (1 / 3, 2 / 3):                    # drift-tube layers
                arc(rin + f * (rout - rin), a0 + 0.8, a1 - 0.8, color=P,
                    lw=0.3, alpha=0.5)
            for a in np.arange(a0 + 6, a1 - 2, 6):      # drift cells
                radial(a, rin + 0.012, rout - 0.012, color=P, lw=0.3,
                       alpha=0.3)

    # -------------------------------------------- wedge borders + ticks ---
    for a in (TH0, TH1):
        radial(a, 0.38, R_TOP, color=W, lw=0.5, alpha=0.9)
        for r in np.arange(0.45, 1.24, 0.1):            # scale tick marks
            (xa, ya) = pol(r, a)
            out = 1 if a == TH0 else -1
            dx, dy = -np.sin(d2r(a)) * out, np.cos(d2r(a)) * out
            plot([xa, xa - dx * 0.028], [ya, ya - dy * 0.028], color=P,
                 lw=0.3, alpha=0.6, zorder=6)

    # ================================================ particle signatures =
    rng = np.random.default_rng(77)

    def track(th_of_r, r0, r1, color=W, lw=0.5, alpha=0.95, ls='-', z=6):
        rs = np.linspace(r0, r1, 140)
        th = d2r(th_of_r(rs))
        xs, ys = rs * np.cos(th), rs * np.sin(th)
        plot(xs, ys, color=color, lw=lw, alpha=alpha, ls=ls, zorder=z,
             solid_capstyle='round')
        return xs, ys                                   # design coords

    def hits(th_of_r, r0, r1, color=W, ms=0.55, z=6):
        for rr in R_TRK:
            if r0 <= rr <= r1:
                hx, hy = pol(rr, th_of_r(np.array([rr]))[0])
                plot([hx], [hy], 'o', ms=ms, color=color, alpha=0.9,
                     zorder=z)

    def clamp(px, py, rmax):
        d = np.hypot(px, py)
        if d > rmax:
            px, py = px * rmax / d, py * rmax / d
        return px, py

    def shower(sx, sy, axis_deg, r_hi, n_seed, spread, core_len, core_hw,
               color=W, lw=0.35, z=6, dots=True):
        """Mini branching cascade confined to r < r_hi, plus a faint core."""
        a0 = d2r(axis_deg)
        ux, uy = np.cos(a0), np.sin(a0)
        nx, ny = -uy, ux
        # faint core with a whisper of an outline
        fx, fy = sx + ux * core_len * 0.95, sy + uy * core_len * 0.95
        fx, fy = clamp(fx, fy, r_hi - 0.01)
        polygon([(sx, sy),
                 (fx + nx * core_hw, fy + ny * core_hw),
                 (fx + ux * 0.02, fy + uy * 0.02),
                 (fx - nx * core_hw, fy - ny * core_hw)],
                facecolor=color, edgecolor=color, lw=0.3, alpha=0.13,
                zorder=5)
        # fluff: short strokes and specks scattered through the core body
        for _ in range(n_seed + 2):
            t = rng.uniform(0.3, 0.85)
            off = rng.normal(0, core_hw * 0.9)
            px, py = sx + ux * t * core_len + nx * off, \
                sy + uy * t * core_len + ny * off
            aa = a0 + rng.normal(0, spread * 1.8)
            LL = rng.uniform(0.028, 0.05)
            qx, qy = clamp(px + LL * np.cos(aa), py + LL * np.sin(aa),
                           r_hi - 0.008)
            plot([px, qx], [py, qy], color=color, lw=lw * 0.8, alpha=0.6,
                 zorder=z)
            if dots and rng.uniform() < 0.5:
                plot([px], [py], 'o', ms=0.5, color=color, alpha=0.6,
                     zorder=z)
        for _ in range(n_seed):
            a = a0 + rng.normal(0, spread)
            L = rng.uniform(0.45, 0.8) * core_len
            ex, ey = clamp(sx + L * np.cos(a), sy + L * np.sin(a),
                           r_hi - 0.012)
            t = np.linspace(0, 1, 24)
            bow = rng.uniform(-0.1, 0.1) * L
            bx = sx + t * (ex - sx) + np.sin(np.pi * t) * bow * (-np.sin(a))
            by = sy + t * (ey - sy) + np.sin(np.pi * t) * bow * (np.cos(a))
            plot(bx, by, color=color, lw=lw, alpha=0.8, zorder=z,
                 solid_capstyle='round')
            for _ in range(int(rng.integers(1, 3))):    # secondary branches
                aa = a + rng.normal(0, spread * 1.4)
                LL = L * rng.uniform(0.3, 0.5)
                gx, gy = clamp(ex + LL * np.cos(aa), ey + LL * np.sin(aa),
                               r_hi - 0.008)
                plot([ex, gx], [ey, gy], color=color, lw=lw * 0.85,
                     alpha=0.7, zorder=z)
                if dots:
                    plot([gx], [gy], 'o', ms=0.6, color=color, alpha=0.8,
                         zorder=z)

    # -- low-p_T looper: decaying spiral curling up inside the tracker
    lc = (0.065, 0.145)
    phi = d2r(np.linspace(-95, 460, 260))
    lr = 0.145 - 0.095 * (phi - phi[0]) / (phi[-1] - phi[0])
    plot(lc[0] + lr * np.cos(phi), lc[1] + lr * np.sin(phi), color=P,
         lw=0.4, alpha=0.75, zorder=6, solid_capstyle='round')

    # -- electron: bends clockwise, showers compactly in the ECAL
    e_th = lambda r: 73 - 13 * ((r - 0.05) / 0.35) ** 1.5
    xs, ys = track(e_th, 0.05, R_EC[0], color=W, lw=0.5)
    hits(e_th, 0.05, R_EC[0])
    e_dir = np.degrees(np.arctan2(ys[-1] - ys[-3], xs[-1] - xs[-3]))
    shower(xs[-1], ys[-1], e_dir, R_EC[1], 4, 0.11, 0.12, 0.028)

    # -- photon: dashed, no track, showers in the ECAL
    track(lambda r: 83 + 0 * r, 0.06, R_EC[0] - 0.008, color=W, lw=0.45,
          alpha=0.95, ls=(0, (2.6, 1.6)))
    gx, gy = pol(R_EC[0], 83)
    shower(gx, gy, 83, R_EC[1], 3, 0.09, 0.115, 0.024)

    # -- charged hadron: MIP through the ECAL, wide shower in the HCAL
    h_th = lambda r: 100 + 7 * ((r - 0.05) / 0.55) ** 1.5
    xs, ys = track(h_th, 0.05, R_HC[0], color=W, lw=0.5)
    hits(h_th, 0.05, R_HC[0])
    h_dir = np.degrees(np.arctan2(ys[-1] - ys[-3], xs[-1] - xs[-3]))
    shower(xs[-1] + 0.006, ys[-1] + 0.006, h_dir, R_HC[1], 6, 0.15, 0.17,
           0.05, lw=0.32)

    # -- neutral hadron: dashed, no track, showers in the HCAL
    track(lambda r: 67 + 0 * r, 0.06, R_HC[0] - 0.008, color=W, lw=0.4,
          alpha=0.85, ls=(0, (2.6, 1.6)))
    kx, ky = pol(R_HC[0], 67)
    shower(kx, ky, 67, R_HC[1], 4, 0.13, 0.155, 0.04, lw=0.32)

    # -- neutrino: faint dotted, escapes everything
    track(lambda r: 95 + 0 * r, 0.05, R_TOP + 0.03, color=P, lw=0.4,
          alpha=0.32, ls=(0, (0.8, 2.4)))

    # blueprint center line of the wedge (dash-dot symmetry axis)
    radial(90, 0.055, R_TOP + 0.035, color=P, lw=0.3, alpha=0.3,
           ls=(0, (4.0, 1.5, 0.8, 1.5)))

    # -- muon: orange, sagitta flips outside the solenoid, arrowhead tip
    def mu_th(r):
        r = np.asarray(r, dtype=float)
        t1 = 116 + 8 * np.clip((r - 0.05) / 0.85, 0, 1) ** 1.6
        s = 8 * 1.6 / 0.85
        t2 = 124 + s * (r - 0.90) - 45 * (r - 0.90) ** 2
        return np.where(r <= 0.90, t1, t2)

    xs, ys = track(mu_th, 0.05, R_TOP + 0.015, color=O, lw=0.7, z=7)
    hits(mu_th, 0.05, 0.36, color=O, ms=0.7, z=7)
    # drift-tube hit ticks in both stations
    for rr in (1.02, 1.05, 1.08, 1.18, 1.21, 1.24):
        a = mu_th(np.array([rr]))[0]
        hx, hy = pol(rr, a)
        dx, dy = -np.sin(d2r(a)), np.cos(d2r(a))
        plot([hx - dx * 0.016, hx + dx * 0.016],
             [hy - dy * 0.016, hy + dy * 0.016], color=O, lw=0.5,
             alpha=0.95, zorder=7)
    # arrowhead (design coords; rotates rigidly with the track)
    tip = np.array([xs[-1], ys[-1]])
    tdir = np.arctan2(ys[-1] - ys[-4], xs[-1] - xs[-4])
    for da in (+2.65, -2.65):
        plot([tip[0], tip[0] + 0.045 * np.cos(tdir + da)],
             [tip[1], tip[1] + 0.045 * np.sin(tdir + da)], color=O,
             lw=0.7, zorder=7, solid_capstyle='round')


# ------------------------------------------------------------ test --------
if __name__ == "__main__":
    import pathlib

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    HERE = pathlib.Path(__file__).resolve().parent
    PLATE = "#0d2c56"

    fig = plt.figure(figsize=(8.5, 11))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 8.5)
    ax.set_ylim(0, 11)
    ax.set_aspect('equal')
    ax.axis('off')
    fig.patch.set_facecolor(PLATE)
    ax.set_facecolor(PLATE)

    draw_detector_slice(ax, (2.0, 8.0), rot_deg=0.0, scale=1.0)
    draw_detector_slice(ax, (1.0, 9.8), rot_deg=-35.0, scale=0.9, alpha=0.6)
    draw_detector_slice(ax, (6.5, 4.0), rot_deg=30.0, scale=0.6)

    out = HERE / "test.png"
    fig.savefig(out, dpi=150, facecolor=PLATE)
    print(f"wrote {out}")
