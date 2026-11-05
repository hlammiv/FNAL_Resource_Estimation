#!/usr/bin/env python3
"""Reusable "quantum-vacuum froth" loop-cloud ornament.

A roughly spherical froth of virtual-particle processes: gluon coils,
photon wavy lines, fermion loops with direction arrows, quark-antiquark
bubbles, nested two-loop structures, all joined at vertex dots into one
seething web. Brightness and line weight fall off from the center to the
rim so the cloud reads as a ball. One orange two-loop self-energy sits
at the core with an orange Higgs-like dashed coupling.

All geometry is generated in *design coordinates* (design center at the
origin) and pushed through a rotate -> scale -> translate map into data
coordinates, so the piece can be placed anywhere, rotated, scaled and
faded. Self-contained: only numpy + matplotlib. Deterministic (fixed
seed).
"""
import numpy as np
from matplotlib.patches import Circle, FancyArrowPatch

WHITE = "#f2f6fb"
PALE = "#9fc4e8"
ORANGE = "#e8833a"

RX, RY = 0.80, 0.62          # design semi-axes of the cloud
DESIGN_RADIUS = 0.85         # max stroke extent from the design center
SEED = 77


# ------------------------------------------------------- design-space geom --

def _rho(x, y):
    """Normalized elliptical radius in design space (0 center .. 1 rim)."""
    return float(np.hypot(x / RX, y / RY))


def _style(r):
    """(color, alpha, lw) as a function of normalized radius 0..1."""
    r = min(r, 1.0)
    color = WHITE if r < 0.45 else PALE
    alpha = max(1.00 - 0.70 * r ** 1.4, 0.20)
    lw = max(0.75 - 0.40 * r, 0.30)
    return color, alpha, lw


def _seg_point_d(p, q, c):
    """Distance from point c to segment p-q (design space)."""
    p, q, c = (np.asarray(v, float) for v in (p, q, c))
    d = q - p
    L2 = float(d @ d)
    if L2 == 0:
        return float(np.hypot(*(c - p)))
    t = float(np.clip((c - p) @ d / L2, 0, 1))
    return float(np.hypot(*(c - (p + t * d))))


# ------------------------------------------------------------------ piece --

def draw_loop_cloud(ax, center, rot_deg=0.0, scale=1.0, alpha=1.0, z=1):
    """Virtual-particle loop-cloud ornament. center = data coords of the
    design center; rot_deg CCW; alpha multiplies every element's alpha;
    artists at zorder z..z+2."""
    rng = np.random.default_rng(SEED)

    # ---- design -> data transform: rotate(rot_deg) -> scale -> translate --
    th = np.deg2rad(rot_deg)
    R = np.array([[np.cos(th), -np.sin(th)],
                  [np.sin(th), np.cos(th)]])
    off = np.array([float(center[0]), float(center[1])])

    def T(pts):
        """(N,2) or (2,) design coords -> (N,2) data coords."""
        pts = np.atleast_2d(np.asarray(pts, float))
        return pts @ R.T * scale + off

    def T1(p):
        return tuple(T(p)[0])

    def A(a):                       # global alpha multiplier
        return float(np.clip(a * alpha, 0.0, 1.0))

    def LW(w):                      # stroke widths track the scale
        return max(w * scale, 0.25)

    def MS(m):                      # marker sizes track the scale
        return max(m * scale, 0.7)

    def MUT(m):                     # arrowhead sizes track the scale
        return max(m * scale, 2.2)

    # ---- primitive draw helpers (inputs in design coords) -----------------
    def line(pts, color, a, lw, zz, ls="-"):
        P = T(pts)
        ax.plot(P[:, 0], P[:, 1], color=color, lw=LW(lw), ls=ls,
                alpha=A(a), zorder=zz, solid_capstyle="round")

    def dot(p, ms, color, a, zz):
        x, y = T1(p)
        ax.plot([x], [y], "o", ms=MS(ms), color=color, alpha=A(a),
                zorder=zz)

    def head(p, ang, color, a, lw=0.35, mut=4.0, zz=z):
        """Tiny arrowhead at design point p along design angle ang."""
        p0 = (p[0] - 0.025 * np.cos(ang), p[1] - 0.025 * np.sin(ang))
        ax.add_patch(FancyArrowPatch(
            T1(p0), T1(p), arrowstyle="-|>", mutation_scale=MUT(mut),
            color=color, lw=LW(lw), alpha=A(a), zorder=zz,
            shrinkA=0, shrinkB=0))

    def circle(c, r, color, a, lw, zz):
        ax.add_patch(Circle(T1(c), r * scale, facecolor="none",
                            edgecolor=color, lw=LW(lw), alpha=A(a),
                            zorder=zz))

    def wavy(p0, p1, n, amp, color, lw, a, zz):
        p0, p1 = np.asarray(p0, float), np.asarray(p1, float)
        d = p1 - p0
        L = float(np.hypot(*d))
        if L == 0:
            return
        u = d / L
        nv = np.array([-u[1], u[0]])
        t = np.linspace(0, 1, 200)
        pts = (p0[None, :] + np.outer(t * L, u)
               + np.outer(amp * np.sin(2 * np.pi * n * t), nv))
        line(pts, color, a, lw, zz)

    def coil(p0, p1, n, amp, color, lw, a, zz):
        p0, p1 = np.asarray(p0, float), np.asarray(p1, float)
        d = p1 - p0
        L = float(np.hypot(*d))
        if L == 0:
            return
        u = d / L
        nv = np.array([-u[1], u[0]])
        t = np.linspace(0, 1, 400)
        par = t * L + 0.55 * amp * np.sin(2 * np.pi * n * t)
        perp = amp * (1 - np.cos(2 * np.pi * n * t))
        pts = p0[None, :] + np.outer(par, u) + np.outer(perp - amp, nv)
        line(pts, color, a, lw, zz)

    # ---- propagators ------------------------------------------------------
    def fermion(p, q, color, a, lw, bow, fwd, zz=z):
        p, q = np.asarray(p, float), np.asarray(q, float)
        d = q - p
        L = float(np.hypot(*d))
        u = d / L
        nv = np.array([-u[1], u[0]])
        t = np.linspace(0, 1, 40)
        pts = p[None, :] + np.outer(t * L, u) + np.outer(
            np.sin(np.pi * t) * bow * L, nv)
        line(pts, color, a, lw, zz)
        mid = p + 0.5 * L * u + nv * bow * L
        uu = u if fwd else -u
        head(tuple(mid), np.arctan2(uu[1], uu[0]), color, a,
             lw=0.9 * lw, mut=3.4 + 2.2 * lw, zz=zz)

    def photon(p, q, color, a, lw, zz=z):
        L = float(np.hypot(q[0] - p[0], q[1] - p[1]))
        wavy(p, q, int(np.clip(round(L * 13), 3, 6)), 0.021,
             color, lw, a, zz)

    def gluon(p, q, color, a, lw, zz=z):
        L = float(np.hypot(q[0] - p[0], q[1] - p[1]))
        coil(p, q, int(np.clip(round(L * 15), 3, 7)), 0.020,
             color, lw, a, zz)

    def scalar(p, q, color, a, lw, zz=z):
        line([p, q], color, a, lw, zz, ls=(0, (2.4, 1.7)))

    # ---- features ---------------------------------------------------------
    def floop(c, r, color, a, lw, arrows=(75, 255), zz=z):
        """Fermion loop: circle with tangential direction arrows."""
        circle(c, r, color, a, lw, zz)
        for adeg in arrows:
            ang = np.deg2rad(adeg)
            p = (c[0] + r * np.cos(ang), c[1] + r * np.sin(ang))
            head(p, ang + np.pi / 2, color, min(a + 0.05, 1.0),
                 lw=0.9 * lw, mut=3.2 + 2.4 * lw, zz=zz)

    def bubble(p, q, bow, color, a, lw, zz=z):
        """Quark-antiquark bubble: lens of two bowed arcs, opposite arrows."""
        p, q = np.asarray(p, float), np.asarray(q, float)
        d = q - p
        L = float(np.hypot(*d))
        u = d / L
        nv = np.array([-u[1], u[0]])
        t = np.linspace(0, 1, 50)
        for s in (+1.0, -1.0):
            pts = p[None, :] + np.outer(t * L, u) + np.outer(
                np.sin(np.pi * t) * bow * L * s, nv)
            line(pts, color, a, lw, zz)
            mid = p + 0.5 * L * u + nv * bow * L * s
            uu = u * s
            head(tuple(mid), np.arctan2(uu[1], uu[0]), color,
                 min(a + 0.05, 1.0), lw=0.9 * lw, mut=3.2 + 2.4 * lw, zz=zz)

    def sunset(c, r, zz=z + 1):
        """Central two-loop self-energy: orange fermion loop, gluon rung."""
        circle(c, r, ORANGE, 0.95, 0.80, zz)
        coil((c[0] - r + 0.012, c[1]), (c[0] + r - 0.012, c[1]),
             4, 0.024, WHITE, 0.50, 0.95, zz)
        for adeg in (90, 270):
            ang = np.deg2rad(adeg)
            p = (c[0] + r * np.cos(ang), c[1] + r * np.sin(ang))
            head(p, ang + np.pi / 2, ORANGE, 1.0, lw=0.55, mut=4.6, zz=zz)
        for sx in (-1, 1):
            dot((c[0] + sx * r, c[1]), 1.9, WHITE, 0.95, zz + 1)

    def nested_chord(c, r, color, a, lw, zz=z):
        """Two-loop structure: fermion loop with a photon chord across it."""
        floop(c, r, color, a, lw, arrows=(120, 300), zz=zz)
        wavy((c[0] - r, c[1]), (c[0] + r, c[1]), 3, 0.016,
             color, 0.85 * lw, a, zz)
        for sx in (-1, 1):
            dot((c[0] + sx * r, c[1]), 1.5, color,
                min(a + 0.1, 1.0), zz + 1)

    def concentric(c, r_out, r_in, color, a, lw, zz=z):
        """Nested loops joined by two short radial photon spokes."""
        floop(c, r_out, color, a, lw, arrows=(30, 210), zz=zz)
        circle(c, r_in, color, a, 0.85 * lw, zz)
        for adeg in (100, 280):
            ang = np.deg2rad(adeg)
            p0 = (c[0] + r_in * np.cos(ang), c[1] + r_in * np.sin(ang))
            p1 = (c[0] + r_out * np.cos(ang), c[1] + r_out * np.sin(ang))
            line([p0, p1], color, a, 0.85 * lw, zz)

    # ---- feature layout (design coords, center at origin) -----------------
    sunset_c, sunset_r = (0.00, 0.02), 0.16
    loops = [  # (dx, dy, r)
        (-0.34, +0.27, 0.085),
        (+0.33, -0.24, 0.075),
        (+0.45, +0.21, 0.058),
        (-0.30, -0.29, 0.068),
        (+0.63, +0.03, 0.048),
        (-0.63, -0.14, 0.050),
        (+0.30, +0.06, 0.048),
    ]
    nested_c, nested_r = (-0.44, +0.03), 0.105
    conc_c, conc_ro, conc_ri = (+0.17, +0.38), 0.075, 0.042
    bubbles = [  # (p, q, bow)
        ((-0.06, -0.36), (+0.15, -0.44), 0.30),
        ((+0.30, +0.34), (+0.46, +0.43), 0.32),
    ]

    # nodes: (x, y, exclusion radius); obstacles for edge routing
    nodes = [(sunset_c[0], sunset_c[1], sunset_r)]
    obstacles = [(sunset_c, sunset_r + 0.030)]
    for (dx, dy, r) in loops:
        nodes.append((dx, dy, r))
        obstacles.append(((dx, dy), r + 0.028))
    nodes.append((nested_c[0], nested_c[1], nested_r))
    obstacles.append((nested_c, nested_r + 0.028))
    nodes.append((conc_c[0], conc_c[1], conc_ro))
    obstacles.append((conc_c, conc_ro + 0.028))
    for (p, q, bow) in bubbles:
        nodes.append((p[0], p[1], 0.0))
        nodes.append((q[0], q[1], 0.0))
        mid = ((p[0] + q[0]) / 2, (p[1] + q[1]) / 2)
        L = float(np.hypot(q[0] - p[0], q[1] - p[1]))
        obstacles.append((mid, 0.5 * L + 0.015))
    n_feature = len(nodes)

    # anchor for the orange Higgs dash, top-center
    higgs_top = (+0.02, +0.50)
    nodes.append((higgs_top[0], higgs_top[1], 0.0))

    # ---- plain vertices: an even rim ring + blue-noise interior -----------
    plain = []

    def ok(x, y, dmin):
        if any(np.hypot(x - c[0], y - c[1]) < rc + 0.08
               for (c, rc) in obstacles):
            return False
        if any(np.hypot(x - px, y - py) < dmin for (px, py, _) in
               nodes + plain):
            return False
        if abs(x) < 0.10 and y < -0.42:   # keep the bottom lane clean
            return False
        return True

    for k in range(15):        # rim ring -> round silhouette
        a = 2 * np.pi * k / 15 + rng.normal(0, 0.09)
        rr = rng.uniform(0.84, 0.95)
        x = RX * rr * np.cos(a)
        y = RY * rr * np.sin(a)
        if ok(x, y, 0.13):
            plain.append((x, y, 0.0))
    # extra arc fills: bottom + right mid-rim
    for adeg in (-142, -116, -64, -38, 25, -22):
        a = np.deg2rad(adeg + rng.normal(0, 4))
        rr = rng.uniform(0.70, 0.86)
        x = RX * rr * np.cos(a)
        y = RY * rr * np.sin(a)
        if ok(x, y, 0.11):
            plain.append((x, y, 0.0))
    tries = 0
    while len(plain) < 47 and tries < 9000:
        tries += 1
        a = rng.uniform(0, 2 * np.pi)
        rr = np.sqrt(rng.uniform()) * 0.90
        x = RX * rr * np.cos(a)
        y = RY * rr * np.sin(a)
        r = _rho(x, y)
        if ok(x, y, 0.105 + 0.11 * r):
            plain.append((x, y, 0.0))
    nodes += plain
    N = len(nodes)

    # ---- Prim spanning tree so the froth is ONE connected web -------------
    def clear(i, j):
        p, q = nodes[i][:2], nodes[j][:2]
        for (c, rc) in obstacles:
            d = _seg_point_d(p, q, c)
            own = (np.hypot(p[0] - c[0], p[1] - c[1]) <= rc + 0.02 or
                   np.hypot(q[0] - c[0], q[1] - c[1]) <= rc + 0.02)
            if not own and d < rc:
                return False
        return True

    in_tree = {0}
    edges = []
    while len(in_tree) < N:
        best, bd = None, 9e9
        for i in in_tree:
            for j in range(N):
                if j in in_tree:
                    continue
                d = (np.hypot(nodes[i][0] - nodes[j][0],
                              nodes[i][1] - nodes[j][1])
                     - nodes[i][2] - nodes[j][2])
                pen = 0.0 if clear(i, j) else 1.5
                if d + pen < bd:
                    best, bd = (i, j), d + pen
        edges.append(best)
        in_tree.add(best[1])

    # extra cycle edges: nearest non-tree neighbors, sparser at the rim
    have = {tuple(sorted(e)) for e in edges}
    order = sorted(range(N), key=lambda i: _rho(nodes[i][0], nodes[i][1]))
    extra = 0
    for i in order:
        if extra >= 21:
            break
        ri = _rho(nodes[i][0], nodes[i][1])
        if ri > 0.70 and rng.uniform() < 0.65:
            continue
        cand = sorted(
            (j for j in range(N) if j != i
             and tuple(sorted((i, j))) not in have),
            key=lambda j: np.hypot(nodes[i][0] - nodes[j][0],
                                   nodes[i][1] - nodes[j][1]))
        for j in cand[:3]:
            d = (np.hypot(nodes[i][0] - nodes[j][0],
                          nodes[i][1] - nodes[j][1])
                 - nodes[i][2] - nodes[j][2])
            if 0.06 < d < 0.36 and clear(i, j):
                edges.append((i, j))
                have.add(tuple(sorted((i, j))))
                extra += 1
                break

    # ---- draw the web -----------------------------------------------------
    dots = []

    def endpoint(i, other):
        x, y, ex = nodes[i]
        if ex == 0.0:
            return np.array([x, y])
        ox, oy = nodes[other][0], nodes[other][1]
        u = np.array([ox - x, oy - y])
        u = u / np.hypot(*u)
        return np.array([x, y]) + (ex + 0.008) * u

    for (i, j) in edges:
        p = endpoint(i, j)
        q = endpoint(j, i)
        L = float(np.hypot(*(q - p)))
        if L < 0.05:
            continue
        mid = 0.5 * (p + q)
        r_mid = _rho(mid[0], mid[1])
        color, a, lw = _style(r_mid)
        u = rng.uniform()
        if L < 0.14:               # too short to wiggle legibly
            cuts = (0.62, 0.62, 0.62, 1.0)
        elif r_mid > 0.72:         # calmer types at the faint rim
            cuts = (0.46, 0.66, 0.80, 1.0)
        else:
            cuts = (0.36, 0.66, 0.90, 1.0)
        if u < cuts[0]:
            fermion(p, q, color, a, lw,
                    bow=rng.choice([-1, 1]) * rng.uniform(0.08, 0.17),
                    fwd=bool(rng.integers(0, 2)))
        elif u < cuts[1]:
            photon(p, q, color, a, lw)
        elif u < cuts[2]:
            gluon(p, q, color, a, lw)
        else:
            scalar(p, q, color, a, lw)
        dots.append(tuple(p))
        dots.append(tuple(q))

    # ---- features ---------------------------------------------------------
    sunset(sunset_c, sunset_r)
    for (dx, dy, r) in loops:
        c = (dx, dy)
        color, a, lw = _style(_rho(*c))
        floop(c, r, color, a, lw,
              arrows=tuple(rng.choice([0, 30, 60], 1)[0] + aa
                           for aa in (75, 255)))
    ncol, nal, nlw = _style(_rho(*nested_c))
    nested_chord(nested_c, nested_r, ncol, nal, nlw)
    ccol, cal, clw = _style(_rho(*conc_c))
    concentric(conc_c, conc_ro, conc_ri, ccol, cal, clw)
    for (p, q, bow) in bubbles:
        mid = ((p[0] + q[0]) / 2, (p[1] + q[1]) / 2)
        color, a, lw = _style(_rho(*mid))
        bubble(p, q, bow, color, a, lw)

    # orange Higgs-like dashed coupling: sunset rim up to the top vertex
    uH = np.array([higgs_top[0] - sunset_c[0], higgs_top[1] - sunset_c[1]])
    uH = uH / np.hypot(*uH)
    pS = (sunset_c[0] + (sunset_r + 0.008) * uH[0],
          sunset_c[1] + (sunset_r + 0.008) * uH[1])
    line([pS, higgs_top], ORANGE, 1.0, 0.70, z + 1, ls=(0, (2.2, 1.5)))
    dots += [pS, higgs_top]

    # ---- rim froth: faint outward stubs from the outermost vertices -------
    outer = [n for n in nodes[n_feature:] if _rho(n[0], n[1]) > 0.80]
    outer = sorted(outer, key=lambda n: -_rho(n[0], n[1]))[:7]
    for (x, y, _) in outer:
        a = np.arctan2(y / RY, x / RX)
        a += rng.normal(0, 0.22)
        L = rng.uniform(0.05, 0.085)
        q = (x + L * np.cos(a), y + L * np.sin(a))
        if abs(q[0]) > 0.84 or not (-0.75 < q[1] < 0.61):
            continue
        if rng.uniform() < 0.5:    # plain / dashed only: no tiny scribbles
            line([(x, y), q], PALE, 0.28, 0.32, z)
        else:
            scalar((x, y), q, PALE, 0.28, 0.32)

    # ---- vertex dots ------------------------------------------------------
    seen = []
    for (x, y) in dots:
        if any(np.hypot(x - sx, y - sy) < 0.02 for (sx, sy) in seen):
            continue
        seen.append((x, y))
        r = _rho(x, y)
        color, a, _ = _style(r)
        dot((x, y), max(2.0 - 1.1 * r, 0.8), color,
            min(a + 0.10, 1.0), z + 1)


# --------------------------------------------------------------- test rig --

if __name__ == "__main__":
    import pathlib

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    fig = plt.figure(figsize=(8.5, 11))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 8.5)
    ax.set_ylim(0, 11)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.add_patch(Rectangle((0, 0), 8.5, 11, facecolor="#0d2c56", zorder=0))

    draw_loop_cloud(ax, (2.0, 8.0), rot_deg=0.0, scale=1.0)
    draw_loop_cloud(ax, (1.5, 6.3), rot_deg=15.0, scale=0.6, alpha=0.8)
    draw_loop_cloud(ax, (6.5, 2.5), rot_deg=-30.0, scale=0.45)

    out = pathlib.Path(__file__).with_name("test.png")
    fig.savefig(out, dpi=150)
    print(f"wrote {out}")
