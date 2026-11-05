#!/usr/bin/env python3
"""Cosmology ornament: the expansion history of the universe drawn as an
engineering blueprint.

A horn/light-cone opens from a point (the inflation burst, in orange) along a
dash-dot engineering centerline.  The comoving grid is drawn as curving
meridians and latitude cross-sections inside the horn; epoch tick-marks with
tiny labels (EW, QCD, BBN, CMB) hang from the bottom silhouette on a common
dimension baseline.  A first-order phase-transition bubble cluster (orange)
sits near the EW epoch with gravitational-wave ripple arcs radiating toward
the mouth.  The CMB is a speckled anisotropy band near the mouth; a small
cosmic web of filaments and cluster nodes spills out beyond it.  Construction
guides (opening-angle rays + dimension arc, an a(t) callout, a time arrow)
complete the blueprint idiom.

Exposes draw_cosmology(ax, center, rot_deg, scale, alpha, z): all geometry is
authored in a design frame about the origin and pushed through
rotate(rot_deg) -> scale -> translate(center) at data level.  Deterministic.
"""

import numpy as np
from matplotlib.patches import Arc, Ellipse

# ------------------------------------------------------ self-contained palette
WHITE = "#f2f6fb"
PALE = "#9fc4e8"
ORANGE = "#e8833a"

# ------------------------------------------------- design-frame geometry
_X0, _X1 = -0.78, 0.34          # horn tip / mouth axial positions
_L = _X1 - _X0                  # horn length
_RMAX = 0.335                   # mouth radius
_E = 0.17                       # cross-section foreshortening (ellipse ratio)
_BASE = -0.415                  # epoch-label dimension baseline

DESIGN_RADIUS = 0.88            # max half-extent of the design (data units)

_EPOCHS = [(0.190, "EW"), (0.345, "QCD"), (0.525, "BBN"), (0.800, "CMB")]


def _prof(u):
    """Horn radius vs. fractional axial position u in [0,1]: a fast
    inflationary rise, decelerating power-law growth, late-time flare."""
    u = np.asarray(u, dtype=float)
    infl = 0.10 * (1.0 - np.exp(-u / 0.030))
    grow = 0.62 * np.clip(u, 0.0, None) ** 0.65
    dark = 0.28 * np.clip((u - 0.78) / 0.22, 0.0, None) ** 2.2
    return _RMAX * (infl + grow + dark)


def _xu(u):
    return _X0 + _L * np.asarray(u, dtype=float)


# --------------------------------------------------------------- ornament
def draw_cosmology(ax, center, rot_deg=0.0, scale=1.0, alpha=1.0, z=1):
    """Expansion-history horn ornament.  center = data coords of the design
    origin; rot_deg CCW; alpha multiplies every element's alpha; artists at
    zorder z..z+2."""
    rng = np.random.default_rng(19670615)      # fixed: deterministic art

    _c = np.cos(np.radians(rot_deg))
    _s = np.sin(np.radians(rot_deg))
    cx0, cy0 = center

    def T(xs, ys):
        xs = np.asarray(xs, dtype=float)
        ys = np.asarray(ys, dtype=float)
        return (cx0 + scale * (_c * xs - _s * ys),
                cy0 + scale * (_s * xs + _c * ys))

    def Tp(x, y):                 # scalar point
        tx, ty = T(x, y)
        return float(tx), float(ty)

    def A(a):
        return min(1.0, a * alpha)

    def LW(w):                    # linewidths follow scale, with a floor
        return max(w * scale, 0.18)

    def MS(m):                    # marker sizes follow scale, with a floor
        return max(m * scale, 0.45)

    def FS(f):                    # font sizes follow scale, with a floor
        return max(f * scale, 2.1)

    def P(xs, ys, **kw):
        kw.setdefault("solid_capstyle", "round")
        kw["lw"] = LW(kw["lw"])
        kw["alpha"] = A(kw["alpha"])
        kw["zorder"] = z + kw.pop("zo", 1.0)
        ln, = ax.plot(*T(xs, ys), **kw)
        return ln

    def dots(xs, ys, ms, color, a, zo=1.0):
        ax.plot(*T(xs, ys), "o", ms=MS(ms), color=color, alpha=A(a),
                zorder=z + zo, markeredgewidth=0, linestyle="none")

    def sect(x, r, **kw):
        """Latitude cross-section: front half solid, back half dotted."""
        cx, cy = Tp(x, 0.0)
        w, h = 2 * _E * r * scale, 2 * r * scale
        lw, a_f, zo = kw.get("lw", 0.3), kw.get("a", 0.4), kw.get("zo", 1.0)
        col = kw.get("color", PALE)
        ax.add_patch(Arc((cx, cy), w, h, angle=rot_deg, theta1=-90.0,
                         theta2=90.0, edgecolor=col, lw=LW(lw),
                         alpha=A(a_f), zorder=z + zo))
        ax.add_patch(Arc((cx, cy), w, h, angle=rot_deg, theta1=90.0,
                         theta2=270.0, edgecolor=col, lw=LW(0.8 * lw),
                         alpha=A(0.42 * a_f), zorder=z + zo - 0.2,
                         linestyle=(0, (1.2, 1.8))))

    def txt(x, y, s, f, color=PALE, a=0.8, zo=1.6, **kw):
        kw.setdefault("ha", "center")
        kw.setdefault("va", "center")
        if "$" in s:
            kw.setdefault("math_fontfamily", "cm")
        ax.text(*Tp(x, y), s, fontsize=FS(f), color=color, alpha=A(a),
                rotation=rot_deg, rotation_mode="anchor", zorder=z + zo, **kw)

    uu = np.linspace(0.0, 1.0, 220)
    xx, rr = _xu(uu), _prof(uu)

    # ---- engineering centerline with epoch ticks and time arrow --------
    P([-0.86, 0.60], [0.0, 0.0], color=PALE, lw=0.3,
      ls=(0, (5.0, 1.6, 1.0, 1.6)), alpha=0.42, zo=0.3)
    P([0.60, 0.865], [0.0, 0.0], color=PALE, lw=0.3, alpha=0.42, zo=0.3)
    P([0.865, 0.838], [0.0, 0.014], color=PALE, lw=0.32, alpha=0.55, zo=0.3)
    P([0.865, 0.838], [0.0, -0.014], color=PALE, lw=0.32, alpha=0.55, zo=0.3)
    txt(0.862, 0.038, r"$t$", 3.4, a=0.6)
    for u_m in (0.118, 0.265, 0.435, 0.615, 0.705, 0.905):   # minor ticks
        xm = float(_xu(u_m))
        P([xm, xm], [-0.009, 0.009], color=PALE, lw=0.28, alpha=0.45, zo=0.8)

    # ---- construction guides: opening-angle rays + dimension arc -------
    for sgn in (+1.0, -1.0):
        P([_X0, _X1], [0.0, sgn * _RMAX], color=PALE, lw=0.26,
          ls=(0, (1.0, 2.2)), alpha=0.24, zo=0.4)
    half = np.degrees(np.arctan2(_RMAX, _L))
    tx0, ty0 = Tp(_X0, 0.0)
    ax.add_patch(Arc((tx0, ty0), 2 * 0.195 * scale, 2 * 0.195 * scale,
                     angle=rot_deg, theta1=-half, theta2=half, edgecolor=PALE,
                     lw=LW(0.26), alpha=A(0.42), zorder=z + 0.4))
    txt(_X0 + 0.235, -0.054, r"$2\theta$", 2.6, a=0.45)

    # faint construction-circle fragments framing the whole piece
    for th1 in (118.0, 298.0):
        ax.add_patch(Arc(Tp(0.0, 0.0), 2 * 0.84 * scale, 2 * 0.84 * scale,
                         angle=rot_deg, theta1=th1, theta2=th1 + 38.0,
                         edgecolor=PALE, lw=LW(0.24), alpha=A(0.16),
                         zorder=z + 0.2, linestyle=(0, (1.0, 2.4))))

    # ---- horn silhouette (white) and comoving meridians (pale) ---------
    P(xx, rr, color=WHITE, lw=0.62, alpha=0.9, zo=1.4)
    P(xx, -rr, color=WHITE, lw=0.62, alpha=0.9, zo=1.4)
    for phi in (32.0, 62.0, 118.0, 148.0):
        cph, sph = np.cos(np.radians(phi)), np.sin(np.radians(phi))
        um = np.linspace(0.012, 1.0, 180)
        rm = _prof(um)
        P(_xu(um) + _E * rm * sph, rm * cph, color=PALE, lw=0.3,
          alpha=0.5, zo=1.0)

    # ---- latitude cross-sections: the comoving grid --------------------
    for u, kind in [(0.058, "r"), (0.118, "f"), (0.190, "e"), (0.265, "f"),
                    (0.345, "e"), (0.435, "f"), (0.525, "e"), (0.615, "f"),
                    (0.705, "f"), (0.800, "e"), (0.905, "f")]:
        r = float(_prof(u))
        if kind == "f":
            sect(_xu(u), r, lw=0.26, a=0.32)
        elif kind == "e":
            sect(_xu(u), r, lw=0.4, a=0.55)
        else:                      # end of inflation / reheating
            sect(_xu(u), r, lw=0.42, a=0.7, color=WHITE)

    # ---- mouth: double rim + very faint pane ---------------------------
    mcx, mcy = Tp(_X1, 0.0)
    ax.add_patch(Ellipse((mcx, mcy), 2 * _E * _RMAX * scale,
                         2 * _RMAX * scale, angle=rot_deg, facecolor=WHITE,
                         alpha=A(0.05), edgecolor="none", zorder=z + 0.5))
    ax.add_patch(Ellipse((mcx, mcy), 2 * _E * _RMAX * scale,
                         2 * _RMAX * scale, angle=rot_deg, facecolor="none",
                         edgecolor=WHITE, lw=LW(0.62), alpha=A(0.9),
                         zorder=z + 1.5))
    ax.add_patch(Ellipse((mcx, mcy), 2 * _E * _RMAX * 0.9 * scale,
                         2 * _RMAX * 0.94 * scale, angle=rot_deg,
                         facecolor="none", edgecolor=PALE, lw=LW(0.26),
                         alpha=A(0.4), zorder=z + 1.4))
    # radius dimension inside the mouth: today's scale factor a0
    psi_d = np.radians(148.0)
    rdx = _X1 + _E * _RMAX * 0.92 * np.sin(psi_d)
    rdy = _RMAX * 0.92 * np.cos(psi_d)
    P([_X1, rdx], [0.0, rdy], color=PALE, lw=0.26, alpha=0.55, zo=1.55)
    dots([_X1], [0.0], 0.7, PALE, 0.7, zo=1.55)
    txt(_X1 + 0.052, -0.148, r"$a_0$", 2.9, a=0.6, zo=1.6)

    # ---- inflation burst at the tip (THE orange accent) ----------------
    for k, da in enumerate((-58.0, -38.0, -19.0, 0.0, 19.0, 38.0, 58.0)):
        ang = np.radians(180.0 + da)
        r0 = 0.028
        r1 = r0 + (0.085 if k % 2 == 0 else 0.052)
        P([_X0 + r0 * np.cos(ang), _X0 + r1 * np.cos(ang)],
          [r0 * np.sin(ang), r1 * np.sin(ang)], color=ORANGE, lw=0.46,
          alpha=0.9, zo=1.8)
    for da in (-48.0, 48.0):       # interleaved fainter pale rays
        ang = np.radians(180.0 + da)
        P([_X0 + 0.030 * np.cos(ang), _X0 + 0.105 * np.cos(ang)],
          [0.030 * np.sin(ang), 0.105 * np.sin(ang)], color=PALE, lw=0.28,
          alpha=0.45, zo=1.2)
    ax.add_patch(Arc((tx0, ty0), 2 * 0.021 * scale, 2 * 0.021 * scale,
                     angle=rot_deg, theta1=0.0, theta2=360.0,
                     edgecolor=ORANGE, lw=LW(0.5), alpha=A(0.95),
                     zorder=z + 1.9))
    dots([_X0], [0.0], 1.25, ORANGE, 0.95, zo=1.9)
    for ang in (0.0, 90.0, 180.0, 270.0):   # datum crosshair ticks
        a_r = np.radians(ang + rot_deg * 0.0)
        P([_X0 + 0.027 * np.cos(a_r), _X0 + 0.040 * np.cos(a_r)],
          [0.027 * np.sin(a_r), 0.040 * np.sin(a_r)], color=ORANGE,
          lw=0.34, alpha=0.8, zo=1.8)

    # ---- primordial fluctuations: two tiny squiggles after the tip -----
    for y0, u0 in ((0.024, 0.075), (-0.020, 0.105)):
        xs = _xu(u0) + np.linspace(0.0, 0.062, 60)
        ys = y0 + 0.0075 * np.sin(2 * np.pi * (xs - xs[0]) / 0.021)
        P(xs, ys, color=PALE, lw=0.26, alpha=0.5, zo=1.0)

    # ---- first-order transition: bubble nucleation + GW ripples --------
    bx, by = _xu(0.215), -0.030
    for (dx, dy, r) in ((0.0, 0.0, 0.020), (-0.024, 0.020, 0.011),
                        (0.020, 0.026, 0.008)):
        ax.add_patch(Arc(Tp(bx + dx, by + dy), 2 * r * scale, 2 * r * scale,
                         angle=rot_deg, theta1=0.0, theta2=360.0,
                         edgecolor=ORANGE, lw=LW(0.42), alpha=A(0.9),
                         zorder=z + 1.8))
    ax.add_patch(Arc(Tp(bx, by), 2 * 0.032 * scale, 2 * 0.032 * scale,
                     angle=rot_deg, theta1=25.0, theta2=245.0,
                     edgecolor=ORANGE, lw=LW(0.28), alpha=A(0.55),
                     zorder=z + 1.7))
    for r, a_g in ((0.060, 0.72), (0.096, 0.55), (0.132, 0.40)):
        ax.add_patch(Arc(Tp(bx, by), 2 * r * scale, 2 * r * scale,
                         angle=rot_deg, theta1=-44.0, theta2=44.0,
                         edgecolor=PALE, lw=LW(0.38), alpha=A(a_g),
                         zorder=z + 1.1))
    for r, a_g in ((0.052, 0.45), (0.082, 0.32)):   # fainter back-ripples
        ax.add_patch(Arc(Tp(bx, by), 2 * r * scale, 2 * r * scale,
                         angle=rot_deg, theta1=142.0, theta2=218.0,
                         edgecolor=PALE, lw=LW(0.3), alpha=A(a_g),
                         zorder=z + 1.1))

    # ---- BBN: a few light-nuclei glyphs after the BBN section ----------
    for (nx0, ny0, ang_n) in ((_xu(0.560), 0.052, 25.0),
                              (_xu(0.585), -0.045, -40.0),
                              (_xu(0.545), -0.098, 70.0)):
        ca, sa = np.cos(np.radians(ang_n)), np.sin(np.radians(ang_n))
        d_n = 0.0115
        p1 = (nx0 - 0.5 * d_n * ca, ny0 - 0.5 * d_n * sa)
        p2 = (nx0 + 0.5 * d_n * ca, ny0 + 0.5 * d_n * sa)
        P([p1[0], p2[0]], [p1[1], p2[1]], color=PALE, lw=0.26,
          alpha=0.55, zo=1.1)
        dots([p1[0], p2[0]], [p1[1], p2[1]], 0.8, PALE, 0.7, zo=1.2)

    # ---- CMB anisotropy: speckled band near the mouth ------------------
    n_cmb = 112
    u_c = rng.uniform(0.772, 0.842, n_cmb)
    psi = rng.uniform(0.0, 2 * np.pi, n_cmb)
    r_c = _prof(u_c)
    sx = _xu(u_c) + _E * r_c * np.sin(psi)
    sy = r_c * np.cos(psi) * rng.uniform(0.97, 1.0, n_cmb)
    front = np.sin(psi) > 0.0
    a_c = rng.uniform(0.28, 0.75, n_cmb)
    m_c = rng.uniform(0.55, 1.05, n_cmb)
    for i in range(n_cmb):
        col = WHITE if a_c[i] > 0.58 else PALE
        aa = a_c[i] if front[i] else 0.38 * a_c[i]
        dots([sx[i]], [sy[i]], m_c[i], col, aa, zo=1.3 if front[i] else 0.6)

    # ---- large-scale structure: cosmic web spilling from the mouth -----
    def rim(psi_deg):
        pr = np.radians(psi_deg)
        return (_X1 + _E * _RMAX * np.sin(pr), _RMAX * np.cos(pr))

    nodes = {
        "r0": rim(35.0), "r1": rim(72.0), "r2": rim(108.0), "r3": rim(146.0),
        "a0": (0.530, 0.235), "a1": (0.568, 0.058), "a2": (0.528, -0.118),
        "a3": (0.560, -0.262),
        "b0": (0.672, 0.150), "b1": (0.706, -0.038), "b2": (0.660, -0.205),
        "c0": (0.792, 0.070), "c1": (0.772, -0.140),
    }
    jit = {k: (v[0] + rng.normal(0, 0.007), v[1] + rng.normal(0, 0.007))
           for k, v in nodes.items()}

    def edge(k1, k2, bow, lw=0.3, a_e=0.5):
        p = np.array(jit[k1]) if isinstance(k1, str) else np.array(k1)
        q = np.array(jit[k2]) if isinstance(k2, str) else np.array(k2)
        t_ = np.linspace(0.0, 1.0, 40)
        d = q - p
        nvec = np.array([-d[1], d[0]]) / (np.linalg.norm(d) + 1e-12)
        pts = (p[None, :] + t_[:, None] * d[None, :]
               + (bow * np.sin(np.pi * t_))[:, None] * nvec[None, :])
        P(pts[:, 0], pts[:, 1], color=PALE, lw=lw, alpha=a_e, zo=1.0)

    for (k1, k2, bow) in [("r0", "a0", 0.010), ("r1", "a1", -0.008),
                          ("r2", "a2", 0.009), ("r3", "a3", -0.007),
                          ("a0", "a1", -0.006), ("a2", "a3", 0.005),
                          ("a0", "b0", -0.008), ("a1", "b0", 0.010),
                          ("a1", "b1", -0.007), ("a2", "b1", 0.008),
                          ("a3", "b2", -0.009), ("b1", "b2", 0.010),
                          ("b0", "c0", 0.007), ("b1", "c0", -0.008),
                          ("b1", "c1", 0.008), ("b2", "c1", -0.007)]:
        edge(k1, k2, bow)
    # dead-end filament hairs off the web
    edge("a0", (0.565, 0.315), 0.006, lw=0.26, a_e=0.35)
    edge("b2", (0.705, -0.285), -0.006, lw=0.26, a_e=0.35)
    edge("c0", (0.845, 0.125), 0.005, lw=0.26, a_e=0.35)
    edge("c1", (0.828, -0.208), -0.005, lw=0.26, a_e=0.35)

    for k, (nx, ny) in jit.items():
        if k in ("b0", "a2"):      # two cluster nodes: open circle + dot
            ax.add_patch(Arc((Tp(nx, ny)), 2 * 0.016 * scale,
                             2 * 0.016 * scale, angle=rot_deg, theta1=0.0,
                             theta2=360.0, edgecolor=PALE, lw=LW(0.34),
                             alpha=A(0.75), zorder=z + 1.3))
            dots([nx], [ny], 1.0, WHITE, 0.85, zo=1.3)
        elif k.startswith("r"):
            dots([nx], [ny], 0.9, PALE, 0.7, zo=1.3)
        else:
            dots([nx], [ny], 1.2, PALE, 0.75, zo=1.3)
    # faint field galaxies sprinkled around and along the filaments
    gx = rng.uniform(0.45, 0.85, 9)
    gy = rng.uniform(-0.30, 0.28, 9)
    dots(gx, gy, 0.6, PALE, 0.3, zo=0.6)
    for (k1, k2) in [("r1", "a1"), ("a1", "b1"), ("b0", "c0"), ("a3", "b2")]:
        p, q = np.array(jit[k1]), np.array(jit[k2])
        for f in (0.3, 0.62):
            m = p + f * (q - p) + rng.normal(0, 0.006, 2)
            dots([m[0]], [m[1]], 0.65, WHITE, 0.45, zo=1.1)

    # ---- epoch dimension row: ticks, leaders, tiny labels --------------
    for u, lab in _EPOCHS:
        xe, re = float(_xu(u)), float(_prof(u))
        P([xe, xe], [-0.016, 0.016], color=WHITE, lw=0.4, alpha=0.8, zo=1.5)
        P([xe, xe], [-re, _BASE + 0.035], color=PALE, lw=0.24,
          ls=(0, (1.0, 2.0)), alpha=0.35, zo=0.5)
        P([xe - 0.012, xe + 0.012], [_BASE + 0.035, _BASE + 0.035],
          color=PALE, lw=0.32, alpha=0.6, zo=1.0)
        txt(xe, _BASE - 0.003, r"$\mathrm{%s}$" % lab, 3.1, a=0.85, va="top")

    # ---- a(t) callout on the top silhouette ----------------------------
    ux = 0.635
    axx, ayy = float(_xu(ux)), float(_prof(ux))
    P([axx, axx + 0.055, axx + 0.10], [ayy, ayy + 0.085, ayy + 0.085],
      color=PALE, lw=0.26, alpha=0.5, zo=1.0)
    dots([axx], [ayy], 0.85, PALE, 0.7, zo=1.2)
    txt(axx + 0.135, ayy + 0.085, r"$a(t)$", 3.0, a=0.7)


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

    ax.add_patch(Rectangle((0, 0), 8.5, 11, facecolor="#0d2c56", zorder=0))
    ax.add_patch(Rectangle((0.45, 0.45), 7.6, 10.1, fill=False,
                           edgecolor="#9fc4e8", lw=1.0, alpha=0.7,
                           zorder=0.5))

    draw_cosmology(ax, (2.0, 8.0), rot_deg=0.0, scale=1.0)
    draw_cosmology(ax, (6.95, 4.45), rot_deg=-10.0, scale=0.75, alpha=0.8)
    draw_cosmology(ax, (2.0, 3.0), rot_deg=20.0, scale=0.55)

    fig.savefig(HERE / "test.png", dpi=150)
    print("wrote", HERE / "test.png")
