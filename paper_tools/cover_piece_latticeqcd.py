"""Lattice-QCD ornament for the FNAL report cover (blueprint idiom).

One cartoon: spacetime lattice underneath, two protons (3 quarks bound by a
Y of gluon springs) colliding, a compact burst, a color flux-tube that breaks
into a q-qbar pair (string breaking), and outgoing hadrons (mesons, baryon).

Contract:
    draw_lattice_qcd(ax, center, rot_deg=0.0, scale=1.0, alpha=1.0, z=1)

All geometry lives in design coordinates about the origin (DESIGN_RADIUS ~0.8)
and is pushed through one rotate -> scale -> translate helper.
"""

import numpy as np

# ---------------------------------------------------------------- palette --
PLATE = "#0d2c56"
WHITE = "#f2f6fb"
PALE = "#9fc4e8"
ORANGE = "#e8833a"

DESIGN_RADIUS = 0.80

# ---------------------------------------------------------- transform kit --


def _make_tf(center, rot_deg, scale):
    """Return tf(x, y) -> rotated, scaled, translated arrays."""
    th = np.deg2rad(rot_deg)
    c, s = np.cos(th), np.sin(th)
    cx, cy = center

    def tf(x, y):
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        return (c * x - s * y) * scale + cx, (s * x + c * y) * scale + cy

    return tf


def _coil(p0, p1, nloop, amp, npts=160):
    """Gluon curl (prolate cycloid) from p0 to p1; endpoints exact."""
    p0 = np.asarray(p0, float)
    p1 = np.asarray(p1, float)
    d = p1 - p0
    L = float(np.hypot(*d))
    if L < 1e-9:
        return np.array([p0[0]]), np.array([p0[1]])
    ux, uy = d / L
    px, py = -uy, ux
    t = np.linspace(0.0, 1.0, npts)
    # straight leads so the curl does not touch the anchor dots
    lead = 0.10
    ph = 2 * np.pi * nloop * np.clip((t - lead) / (1 - 2 * lead), 0.0, 1.0)
    a_l = t * L - amp * np.sin(ph)          # along-axis (loops when amp big)
    a_p = amp * (1.0 - np.cos(ph))          # one-sided curls
    return p0[0] + ux * a_l + px * a_p, p0[1] + uy * a_l + py * a_p


# ------------------------------------------------------------- the piece --


def draw_lattice_qcd(ax, center, rot_deg=0.0, scale=1.0, alpha=1.0, z=1):
    rng = np.random.default_rng(7)  # deterministic
    tf = _make_tf(center, rot_deg, scale)

    # linewidths / markers with floors
    def lw(w):
        return max(0.30, w * scale)

    def ms(m):
        return max(0.7, m * scale)

    fs = max(2.6, 4.6 * scale)

    def line(x, y, color, w, a, ls="-", zo=0):
        ax.plot(*tf(x, y), color=color, lw=lw(w), ls=ls,
                alpha=np.clip(a * alpha, 0, 1), zorder=z + zo,
                solid_capstyle="round", dash_capstyle="round")

    def dot(x, y, color, m, a, zo=0, open_face=False):
        ax.plot(*tf([x], [y]), marker="o", ms=ms(m), ls="",
                mfc=(PLATE if open_face else color), mec=color,
                mew=lw(0.55), alpha=np.clip(a * alpha, 0, 1), zorder=z + zo)

    def fill(x, y, color, a, zo=0):
        ax.fill(*tf(x, y), color=color,
                alpha=np.clip(a * alpha, 0, 0.15), lw=0, zorder=z + zo)

    def text(x, y, s, color, a, size=1.0, zo=2, style="italic"):
        X, Y = tf([x], [y])
        ax.text(X[0], Y[0], s, color=color, fontsize=fs * size,
                alpha=np.clip(a * alpha, 0, 1), rotation=rot_deg,
                rotation_mode="anchor", ha="center", va="center",
                style=style, family="serif", zorder=z + zo)

    def arrow(p0, p1, color, w, a, ls="-", zo=1, style="-|>", head=3.2):
        (x0,), (y0,) = tf([p0[0]], [p0[1]])
        (x1,), (y1,) = tf([p1[0]], [p1[1]])
        ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
                    xycoords="data", textcoords="data",
                    arrowprops=dict(arrowstyle=style, lw=lw(w), color=color,
                                    alpha=np.clip(a * alpha, 0, 1),
                                    linestyle=ls, shrinkA=0, shrinkB=0,
                                    mutation_scale=max(2.0, head * scale)),
                    zorder=z + zo)

    def circle(cx, cy, r, color, w, a, ls="-", zo=1, fill_a=None):
        t = np.linspace(0, 2 * np.pi, 90)
        if fill_a:
            fill(cx + r * np.cos(t), cy + r * np.sin(t), color, fill_a, zo=zo)
        line(cx + r * np.cos(t), cy + r * np.sin(t), color, w, a, ls=ls, zo=zo)

    # ================================================== 1. spacetime lattice
    su, sv = 0.17, 0.155
    U = np.arange(-4, 5)
    V = np.arange(-3, 4)

    def site(u, v):
        x = u * su * (1.0 - 0.045 * (v + 3))      # mild perspective recession
        y = v * sv * (1.0 - 0.022 * (v + 3))
        return x, y

    def fade(x, y):
        d = np.hypot(x / 0.72, y / 0.50)
        return float(np.clip((1.02 - d) / 0.50, 0.0, 1.0))

    for v in V:                                   # spatial links
        for u in U[:-1]:
            x0, y0 = site(u, v)
            x1, y1 = site(u + 1, v)
            f = fade((x0 + x1) / 2, (y0 + y1) / 2)
            if f > 0.03:
                line([x0, x1], [y0, y1], PALE, 0.42, 0.30 * f, zo=0)
    for u in U:                                   # temporal links
        for v in V[:-1]:
            x0, y0 = site(u, v)
            x1, y1 = site(u, v + 1)
            f = fade((x0 + x1) / 2, (y0 + y1) / 2)
            if f > 0.03:
                line([x0, x1], [y0, y1], PALE, 0.42, 0.30 * f, zo=0)
    for u in U:                                   # sites
        for v in V:
            x0, y0 = site(u, v)
            f = fade(x0, y0)
            if f > 0.06:
                dot(x0, y0, PALE, 1.05, 0.46 * f, zo=0)

    # a few links carry gluon curls (gauge field lives on links)
    for (ua, va), (ub, vb) in [((-3, -2), (-2, -2)), ((2, 2), (3, 2)),
                               ((0, -2), (1, -2)), ((3, -1), (3, -2))]:
        pa, pb = site(ua, va), site(ub, vb)
        f = fade((pa[0] + pb[0]) / 2, (pa[1] + pb[1]) / 2)
        cx, cy = _coil(pa, pb, 3, 0.014)
        line(cx, cy, PALE, 0.5, max(0.45 * f, 0.15), zo=0)

    # two shaded plaquettes with orientation loop-arrows (Wilson loops)
    for (u0, v0), r_arrow, a_ls in [((1, -2), 0.030, "-"),
                                    ((-3, 1), 0.024, (0, (2.0, 1.6)))]:
        corners = [site(u0, v0), site(u0 + 1, v0),
                   site(u0 + 1, v0 + 1), site(u0, v0 + 1)]
        ccx = np.mean([p[0] for p in corners])
        ccy = np.mean([p[1] for p in corners])
        corners = [(ccx + 0.90 * (p[0] - ccx), ccy + 0.90 * (p[1] - ccy))
                   for p in corners]
        xs = [p[0] for p in corners] + [corners[0][0]]
        ys = [p[1] for p in corners] + [corners[0][1]]
        fill(xs, ys, PALE, 0.05, zo=0)
        pcx, pcy = np.mean(xs[:4]), np.mean(ys[:4])
        t = np.linspace(np.deg2rad(20), np.deg2rad(335), 60)
        line(pcx + r_arrow * np.cos(t), pcy + r_arrow * np.sin(t),
             PALE, 0.45, 0.5, ls=a_ls, zo=0)
        # arrowhead at end of arc, along local tangent
        te = t[-1]
        hx, hy = pcx + r_arrow * np.cos(te), pcy + r_arrow * np.sin(te)
        tx, ty = -np.sin(te), np.cos(te)
        nx, ny = np.cos(te), np.sin(te)
        hl = 0.011
        line([hx - hl * (tx + 0.8 * nx), hx], [hy - hl * (ty + 0.8 * ny), hy],
             PALE, 0.45, 0.5, zo=0)
        line([hx - hl * (tx - 0.8 * nx), hx], [hy - hl * (ty - 0.8 * ny), hy],
             PALE, 0.45, 0.5, zo=0)

    # lattice-spacing dimension callout  |<- a ->|
    (ax0, ay0), (ax1, ay1) = site(-4, -3), site(-3, -3)
    yd = ay0 - 0.055
    line([ax0, ax0], [ay0 - 0.014, yd - 0.014], PALE, 0.4, 0.55, zo=0)
    line([ax1, ax1], [ay1 - 0.014, yd - 0.014], PALE, 0.4, 0.55, zo=0)
    arrow((ax0, yd), (ax1, yd), PALE, 0.45, 0.6, style="<->", head=2.6, zo=0)
    text((ax0 + ax1) / 2, yd - 0.035, "a", PALE, 0.75, size=0.85)

    # ============================================= 2. construction geometry
    C = np.array([-0.10, 0.01])                   # collision point
    dirA = np.array([np.cos(np.deg2rad(-12)), np.sin(np.deg2rad(-12))])
    dirB = np.array([np.cos(np.deg2rad(118)), np.sin(np.deg2rad(118))])
    for d, back, fwd in ((dirA, 0.68, 0.48), (dirB, 0.58, 0.50)):
        p0, p1 = C - d * back, C + d * fwd        # beam construction lines
        line([p0[0], p1[0]], [p0[1], p1[1]], PALE, 0.4, 0.22,
             ls=(0, (5.0, 2.0, 1.0, 2.0)), zo=0)

    # ==================================================== 3. hadron builder
    def baryon(cx, cy, r, ang0, a=1.0, zo=1, label=None, lab_ang=135.0):
        """3 quarks + Y of gluon springs inside a dashed boundary."""
        circle(cx, cy, r, PALE, 0.5, 0.55 * a, ls=(0, (1.6, 1.8)),
               zo=zo, fill_a=0.06 * a)
        qs = []
        for k in range(3):
            th = np.deg2rad(ang0 + 120 * k)
            qs.append((cx + 0.64 * r * np.cos(th), cy + 0.64 * r * np.sin(th)))
        for q in qs:                              # Y-junction flux strings
            gx, gy = _coil((cx, cy), q, 2, 0.005)
            line(gx, gy, WHITE, 0.5, 0.80 * a, zo=zo)
        dot(cx, cy, WHITE, 1.0, 0.75 * a, zo=zo)  # junction
        for q in qs:
            dot(q[0], q[1], WHITE, 2.05, 0.95 * a, zo=zo)
        if label:
            la = np.deg2rad(lab_ang)
            text(cx + 1.5 * r * np.cos(la), cy + 1.5 * r * np.sin(la),
                 label, PALE, 0.8 * a, zo=zo + 1)

    def meson(cx, cy, r, ang, a=1.0, zo=1, label=None):
        """quark + antiquark with a single gluon spring."""
        circle(cx, cy, r, PALE, 0.5, 0.5 * a, ls=(0, (1.6, 1.8)),
               zo=zo, fill_a=0.06 * a)
        th = np.deg2rad(ang)
        q0 = (cx - 0.58 * r * np.cos(th), cy - 0.58 * r * np.sin(th))
        q1 = (cx + 0.58 * r * np.cos(th), cy + 0.58 * r * np.sin(th))
        gx, gy = _coil(q0, q1, 2, 0.0075)
        line(gx, gy, WHITE, 0.55, 0.8 * a, zo=zo)
        dot(q0[0], q0[1], WHITE, 2.1, 0.95 * a, zo=zo)
        dot(q1[0], q1[1], WHITE, 2.1, 0.9 * a, zo=zo, open_face=True)
        if label:
            text(cx + 1.4 * r, cy + 1.1 * r, label, PALE, 0.8 * a, zo=zo + 1)

    # ====================================== 4. incoming protons + kinematics
    pA = C - dirA * 0.52                          # from the left
    pB = C - dirB * 0.47                          # from below-right
    baryon(pA[0], pA[1], 0.115, 18, zo=1, label="p", lab_ang=128)
    baryon(pB[0], pB[1], 0.100, -30, zo=1, label="p", lab_ang=-12)
    for pc, d, r in [(pA, dirA, 0.115), (pB, dirB, 0.100)]:
        tail = pc + d * (r + 0.02)
        head = C - d * 0.115
        arrow(tail, head, WHITE, 0.6, 0.8, ls=(0, (2.6, 1.8)), head=3.4)
        # speed strokes behind the hadron
        perp = np.array([-d[1], d[0]])
        for off in (-0.05, 0.0, 0.05):
            b0 = pc - d * (r + 0.030) + perp * off * 1.1
            b1 = pc - d * (r + 0.085 - 0.25 * abs(off)) + perp * off
            line([b0[0], b1[0]], [b0[1], b1[1]], PALE, 0.45, 0.5)

    # ============================================== 5. the collision burst
    nspk = 12
    for k in range(nspk):
        th = 2 * np.pi * k / nspk + 0.13 + rng.uniform(-0.09, 0.09)
        r0 = 0.030 + rng.uniform(0, 0.008)
        r1 = (0.095 if k % 3 == 0 else 0.062) + rng.uniform(0, 0.014)
        col = ORANGE if k % 3 == 0 else WHITE
        line([C[0] + r0 * np.cos(th), C[0] + r1 * np.cos(th)],
             [C[1] + r0 * np.sin(th), C[1] + r1 * np.sin(th)],
             col, 0.6, 0.9, zo=2)
    circle(C[0], C[1], 0.047, PALE, 0.45, 0.6, ls=(0, (1.4, 1.6)), zo=2)
    circle(C[0], C[1], 0.020, ORANGE, 0.55, 0.9, zo=2, fill_a=0.13)

    # ====================================== 6. string breaking (flux tube)
    sax = np.deg2rad(20)                          # vignette axis
    sdir = np.array([np.cos(sax), np.sin(sax)])
    sperp = np.array([-sdir[1], sdir[0]])

    def on_axis(s, w=0.0):
        p = C + sdir * s + sperp * w
        return p[0], p[1]

    def tube(s0, s1, hw):
        """capsule flux tube with faint fill + transverse field rungs."""
        t = np.linspace(-np.pi / 2, np.pi / 2, 30)
        top = [on_axis(s, hw) for s in np.linspace(s0 + hw, s1 - hw, 24)]
        cap1 = [on_axis(s1 - hw + hw * np.cos(u), hw * np.sin(u))
                for u in np.linspace(np.pi / 2, -np.pi / 2, 20)][::-1]
        bot = [on_axis(s, -hw) for s in np.linspace(s1 - hw, s0 + hw, 24)]
        cap0 = [on_axis(s0 + hw - hw * np.cos(u), hw * np.sin(u))
                for u in np.linspace(-np.pi / 2, np.pi / 2, 20)][::-1]
        poly = top + cap1 + bot + cap0
        xs = [p[0] for p in poly]
        ys = [p[1] for p in poly]
        fill(xs, ys, ORANGE, 0.13, zo=1)
        line(xs + [xs[0]], ys + [ys[0]], ORANGE, 0.6, 0.85, zo=2)
        for s in np.linspace(s0 + 2.0 * hw, s1 - 2.0 * hw, 3):
            x0, y0 = on_axis(s, hw * 0.62)
            x1, y1 = on_axis(s, -hw * 0.62)
            line([x0, x1], [y0, y1], WHITE, 0.4, 0.55, zo=2)

    sq0, st1a, st1b, sbrk, st2a, st2b, sq1 = \
        0.135, 0.165, 0.345, 0.395, 0.445, 0.625, 0.655
    hw = 0.026
    tube(st1a, st1b, hw)                          # tube 1: q ... qbar_new
    tube(st2a, st2b, hw)                          # tube 2: q_new ... qbar
    dot(*on_axis(sq0), WHITE, 2.4, 0.95, zo=2)                    # quark
    dot(*on_axis(st1b + 0.018), WHITE, 2.2, 0.95, zo=2, open_face=True)
    dot(*on_axis(st2a - 0.018), WHITE, 2.2, 0.95, zo=2)           # new q
    dot(*on_axis(sq1), WHITE, 2.4, 0.95, zo=2, open_face=True)    # antiquark
    for sgn in (1, -1):                           # snap sparks at the break
        for dth in (np.deg2rad(55), np.deg2rad(112)):
            vx = np.cos(sax + sgn * dth)
            vy = np.sin(sax + sgn * dth)
            b0 = C + sdir * sbrk + np.array([vx, vy]) * 0.014
            b1 = C + sdir * sbrk + np.array([vx, vy]) * 0.040
            line([b0[0], b1[0]], [b0[1], b1[1]], ORANGE, 0.55, 0.9, zo=2)
    text(*(C + sdir * sq0 + np.array([0.005, -0.055])), "q", PALE, 0.95)
    text(*(C + sdir * sq1 + np.array([0.012, 0.055])), "q̄", PALE, 0.95)

    # a radiated gluon curling away from the vertex
    gdir = np.array([np.cos(np.deg2rad(52)), np.sin(np.deg2rad(52))])
    g0, g1 = C + gdir * 0.105, C + gdir * 0.285
    gx, gy = _coil(g0, g1, 3, 0.013)
    line(gx, gy, WHITE, 0.55, 0.75, zo=2)
    text(*(g1 + gdir * 0.045 + np.array([0.028, 0.0])), "g", PALE, 0.85)

    # ================================================= 7. outgoing hadrons
    outs = [
        (np.deg2rad(78), 0.470, "meson", 0.052, "π"),   # pion, up
        (np.deg2rad(203), 0.470, "baryon", 0.075, None),     # baryon, low-L
        (np.deg2rad(133), 0.430, "meson", 0.047, None),      # meson, up-left
    ]
    for th, dist, kind, r, lab in outs:
        d = np.array([np.cos(th), np.sin(th)])
        hc = C + d * dist
        t0 = C + d * 0.115
        t1 = C + d * (dist - r - 0.035)
        arrow(t0, t1, PALE, 0.5, 0.55, ls=(0, (2.2, 2.0)), head=2.8)
        if kind == "meson":
            meson(hc[0], hc[1], r, np.rad2deg(th) + 90, zo=1, label=lab)
        else:
            baryon(hc[0], hc[1], r, np.rad2deg(th) + 35, zo=1, label=lab)


# ------------------------------------------------------------ test harness
if __name__ == "__main__":
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(8.5, 11))
    fig.patch.set_facecolor(PLATE)
    ax.set_facecolor(PLATE)
    ax.set_xlim(0, 8.5)
    ax.set_ylim(0, 11)
    ax.set_aspect("equal")
    ax.axis("off")
    fig.subplots_adjust(left=0, right=1, top=1, bottom=0)

    draw_lattice_qcd(ax, (2.0, 8.0), rot_deg=0.0, scale=1.0, alpha=1.0, z=1)
    draw_lattice_qcd(ax, (1.30, 4.70), rot_deg=8.0, scale=0.75, alpha=0.8, z=1)
    draw_lattice_qcd(ax, (6.5, 3.0), rot_deg=-20.0, scale=0.55, alpha=1.0, z=1)

    out = __file__.rsplit("/", 1)[0] + "/test.png"
    fig.savefig(out, dpi=150, facecolor=PLATE)
    print("wrote", out)
