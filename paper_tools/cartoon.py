"""Drawing primitives shared by the chapter cartoons.

These are schematic, not to scale. Everything works in whatever data coords the
calling axes uses; the cartoon scripts set up a bare canvas with fnalfig.blank()
and then draw in arbitrary units.
"""

import numpy as np
from matplotlib.patches import (Circle, FancyArrowPatch, PathPatch, Polygon,
                                Rectangle)
from matplotlib.path import Path

import fnalfig as ff


def arrow(ax, p0, p1, color=ff.INK, lw=0.9, style="-|>", ls="-", z=5, mut=7,
          shrinkA=0, shrinkB=0, **kw):
    a = FancyArrowPatch(p0, p1, arrowstyle=style, mutation_scale=mut,
                        color=color, lw=lw, linestyle=ls, zorder=z,
                        shrinkA=shrinkA, shrinkB=shrinkB, **kw)
    ax.add_patch(a)
    return a


def wavy(ax, p0, p1, n=6, amp=0.09, color=ff.INK, lw=0.9, z=4, alpha=1.0):
    """Boson-like wavy propagator between two points."""
    p0, p1 = np.asarray(p0, float), np.asarray(p1, float)
    d = p1 - p0
    L = np.hypot(*d)
    if L == 0:
        return
    u = d / L
    nvec = np.array([-u[1], u[0]])
    t = np.linspace(0, 1, 200)
    pts = (p0[None, :] + np.outer(t * L, u)
           + np.outer(amp * np.sin(2 * np.pi * n * t), nvec))
    ax.plot(pts[:, 0], pts[:, 1], color=color, lw=lw, zorder=z,
            alpha=alpha, solid_capstyle="round")


def coil(ax, p0, p1, n=7, amp=0.11, color=ff.INK, lw=0.9, z=4, alpha=1.0):
    """Gluon-like coil: a wave with loops, drawn as a stretched cycloid."""
    p0, p1 = np.asarray(p0, float), np.asarray(p1, float)
    d = p1 - p0
    L = np.hypot(*d)
    u = d / L
    nvec = np.array([-u[1], u[0]])
    t = np.linspace(0, 1, 400)
    par = t * L + 0.55 * amp * np.sin(2 * np.pi * n * t)
    perp = amp * (1 - np.cos(2 * np.pi * n * t))
    pts = p0[None, :] + np.outer(par, u) + np.outer(perp - amp, nvec)
    ax.plot(pts[:, 0], pts[:, 1], color=color, lw=lw, zorder=z,
            alpha=alpha, solid_capstyle="round")


def nucleus(ax, xy, r=0.5, n=9, color=ff.BLUE, alt=ff.ORANGE, seed=3, z=3,
            halo=True):
    """A cluster of nucleons, two species, inside an optional halo."""
    rng = np.random.default_rng(seed)
    x0, y0 = xy
    if halo:
        ax.add_patch(Circle(xy, r * 1.32, facecolor=color, alpha=0.07,
                            edgecolor=color, lw=0.6, ls=(0, (3, 2)), zorder=z - 1))
    pts = []
    rr = r * 0.80
    while len(pts) < n:
        p = rng.uniform(-rr, rr, 2)
        if np.hypot(*p) <= rr and all(np.hypot(*(p - q)) > 0.30 * r for q in pts):
            pts.append(p)
    for i, p in enumerate(pts):
        ax.add_patch(Circle((x0 + p[0], y0 + p[1]), 0.155 * r,
                            facecolor=color if i % 2 else alt,
                            edgecolor="white", lw=0.4, zorder=z))
    return pts


def lattice(ax, x0, y0, nx, ny, dx, color=ff.GREY, lw=0.5, site_ms=2.0,
            site_color=None, z=2, alpha=1.0):
    """Square lattice of sites joined by links."""
    for i in range(nx):
        for j in range(ny):
            x, y = x0 + i * dx, y0 + j * dx
            if i < nx - 1:
                ax.plot([x, x + dx], [y, y], color=color, lw=lw, zorder=z,
                        alpha=alpha)
            if j < ny - 1:
                ax.plot([x, x], [y, y + dx], color=color, lw=lw, zorder=z,
                        alpha=alpha)
    if site_ms:
        xs = [x0 + i * dx for i in range(nx) for _ in range(ny)]
        ys = [y0 + j * dx for _ in range(nx) for j in range(ny)]
        ax.plot(xs, ys, "o", ms=site_ms, color=site_color or color,
                zorder=z + 1, alpha=alpha)


def box(ax, xy, w, h, color=ff.INK, fc="white", lw=0.8, r=0.06, z=5, ls="-",
        alpha=1.0):
    """Rounded-corner box, given by lower-left corner."""
    x, y = xy
    r = min(r, w / 2, h / 2)
    verts, codes = [], []
    corners = [((x + r, y), (x + w - r, y), (x + w, y), (x + w, y + r)),
               ((x + w, y + h - r), (x + w, y + h), (x + w - r, y + h), None),
               ((x + r, y + h), (x, y + h), (x, y + h - r), None),
               ((x, y + r), (x, y), (x + r, y), None)]
    verts.append(corners[0][0]); codes.append(Path.MOVETO)
    verts.append(corners[0][1]); codes.append(Path.LINETO)
    verts += [corners[0][2], corners[0][3]]; codes += [Path.CURVE3, Path.CURVE3]
    verts.append(corners[1][0]); codes.append(Path.LINETO)
    verts += [corners[1][1], corners[1][2]]; codes += [Path.CURVE3, Path.CURVE3]
    verts.append(corners[2][0]); codes.append(Path.LINETO)
    verts += [corners[2][1], corners[2][2]]; codes += [Path.CURVE3, Path.CURVE3]
    verts.append(corners[3][0]); codes.append(Path.LINETO)
    verts += [corners[3][1], corners[3][2]]; codes += [Path.CURVE3, Path.CURVE3]
    codes.append(Path.CLOSEPOLY); verts.append(corners[0][0])
    p = PathPatch(Path(verts, codes), facecolor=fc, edgecolor=color, lw=lw,
                  zorder=z, linestyle=ls, alpha=alpha)
    ax.add_patch(p)
    return p


def txt(ax, xy, s, size=6.5, color=ff.INK, ha="center", va="center", z=7, **kw):
    return ax.text(xy[0], xy[1], s, fontsize=size, color=color, ha=ha, va=va,
                   zorder=z, **kw)


def brace(ax, x, y0, y1, width=0.12, color=ff.GREY, lw=0.7, z=4):
    """Vertical curly-ish brace spanning y0..y1 at x, opening left."""
    ym = 0.5 * (y0 + y1)
    t = np.linspace(0, 1, 60)
    for (ya, yb) in ((y0, ym), (y1, ym)):
        yy = ya + (yb - ya) * t
        xx = x + width * np.sin(np.pi * t) ** 0.7
        ax.plot(xx, yy, color=color, lw=lw, zorder=z)


def wedge_band(ax, xs, lo, hi, color, alpha=0.13, z=1, **kw):
    ax.fill_between(xs, lo, hi, color=color, alpha=alpha, lw=0, zorder=z, **kw)


def cone(ax, apex, half_angle_slope, ymax, color=ff.GREY, alpha=0.10, z=1):
    """Light-cone / expanding-universe wedge from an apex upward."""
    x0, y0 = apex
    dy = ymax - y0
    dx = half_angle_slope * dy
    ax.add_patch(Polygon([[x0, y0], [x0 - dx, ymax], [x0 + dx, ymax]],
                         closed=True, facecolor=color, alpha=alpha, lw=0,
                         zorder=z))


# --- Richer elements, added for the second pass on the chapter cartoons ------

def hadron(ax, xy, r=0.30, quarks=3, color=ff.BLUE, seed=1, z=4, label=None,
           lw=0.9, alpha=0.13):
    """A confined blob with valence quarks inside -- proton, pion, whatever."""
    rng = np.random.default_rng(seed)
    x0, y0 = xy
    ax.add_patch(Circle(xy, r, facecolor=color, alpha=alpha, edgecolor=color,
                        lw=lw, zorder=z))
    for i in range(quarks):
        a = 2 * np.pi * i / quarks + rng.uniform(0, 1)
        rr = r * rng.uniform(0.28, 0.52)
        ax.add_patch(Circle((x0 + rr * np.cos(a), y0 + rr * np.sin(a)),
                            0.30 * r, facecolor=color, edgecolor="white",
                            lw=0.4, zorder=z + 1))
    if label:
        txt(ax, (x0, y0 - r - 0.16), label, size=6.5, color=color)


def detector_rings(ax, xy, radii, color=ff.GREY, lw=0.7, z=1, alpha=0.55):
    """Concentric detector barrels seen end-on."""
    for r in radii:
        ax.add_patch(Circle(xy, r, facecolor="none", edgecolor=color, lw=lw,
                            zorder=z, alpha=alpha))


def spray(ax, xy, n=14, r0=0.35, rmax=1.6, color=ff.ORANGE, seed=2, z=3,
          lw=0.7, spread=360, ang0=0, alpha=1.0):
    """Isotropic-ish spray of tracks from a vertex."""
    rng = np.random.default_rng(seed)
    x0, y0 = xy
    tips = []
    for i in range(n):
        a = np.deg2rad(ang0 + spread * (i / n) + rng.uniform(-6, 6))
        L = rng.uniform(0.7, 1.0) * rmax
        p0 = (x0 + r0 * np.cos(a), y0 + r0 * np.sin(a))
        p1 = (x0 + L * np.cos(a), y0 + L * np.sin(a))
        ax.plot([p0[0], p1[0]], [p0[1], p1[1]], color=color, lw=lw, zorder=z,
                alpha=alpha, solid_capstyle="round")
        tips.append(p1)
    return tips


def cube(ax, origin, size, nx=3, color=ff.GREY, lw=0.5, skew=(0.36, 0.26),
         z=2, alpha=0.9, site_ms=1.7, site_color=None):
    """Oblique-projected cubic lattice: a 3D box drawn as nx^3 sites."""
    sx, sy = skew
    ox, oy = origin
    d = size / max(nx - 1, 1)

    def P(i, j, k):
        return (ox + i * d + sx * k * d, oy + j * d + sy * k * d)

    for k in range(nx):
        for i in range(nx):
            for j in range(nx):
                if i < nx - 1:
                    a, b = P(i, j, k), P(i + 1, j, k)
                    ax.plot(*zip(a, b), color=color, lw=lw, zorder=z, alpha=alpha)
                if j < nx - 1:
                    a, b = P(i, j, k), P(i, j + 1, k)
                    ax.plot(*zip(a, b), color=color, lw=lw, zorder=z, alpha=alpha)
                if k < nx - 1:
                    a, b = P(i, j, k), P(i, j, k + 1)
                    ax.plot(*zip(a, b), color=color, lw=lw, zorder=z,
                            alpha=alpha * 0.75)
    if site_ms:
        pts = [P(i, j, k) for k in range(nx) for i in range(nx) for j in range(nx)]
        ax.plot([p[0] for p in pts], [p[1] for p in pts], "o", ms=site_ms,
                color=site_color or color, zorder=z + 1, alpha=alpha)
    return P


def matrix_box(ax, xy, w, h, color=ff.INK, lw=0.9, z=5, label=None,
               fill=None, band=None, fc="white"):
    """Square-bracket matrix with an optional sparsity band drawn inside."""
    x, y = xy
    ax.add_patch(Rectangle((x, y), w, h, facecolor=fc, edgecolor="none", zorder=z))
    t = 0.16 * w
    for xx, s in ((x, 1), (x + w, -1)):
        ax.plot([xx + s * t, xx, xx, xx + s * t], [y, y, y + h, y + h],
                color=color, lw=lw, zorder=z + 1, solid_joinstyle="miter")
    if band is not None:
        n = band
        d = min(w, h) / n
        for i in range(n):
            for j in range(n):
                if abs(i - j) <= 1:
                    ax.add_patch(Rectangle((x + 0.10 * w + j * d * 0.8,
                                            y + 0.10 * h + i * d * 0.8),
                                           d * 0.6, d * 0.6, facecolor=fill or color,
                                           alpha=0.35, edgecolor="none", zorder=z + 1))
    if label:
        txt(ax, (x + w / 2, y + h + 0.10), label, size=7, color=color)


def circuit(ax, xy, w, nwires=3, gates=None, color=ff.INK, lw=0.7, z=5,
            wire_color=None, label=None):
    """Miniature quantum circuit: horizontal wires with boxed gates.

    `gates` is a list of (frac_x, wire_index, text) or (frac_x, (w0,w1), text)
    for a multi-wire box.
    """
    x, y = xy
    dy = 0.26
    for i in range(nwires):
        ax.plot([x, x + w], [y + i * dy, y + i * dy],
                color=wire_color or ff.GREY, lw=lw, zorder=z)
    for g in gates or []:
        fx, wires, text = g
        if isinstance(wires, tuple):
            w0, w1 = wires
        else:
            w0 = w1 = wires
        gy0, gy1 = y + w0 * dy, y + w1 * dy
        bw, bh = max(0.30, 0.080 * len(text) + 0.20), (gy1 - gy0) + 0.19
        box(ax, (x + fx * w - bw / 2, gy0 - 0.085), bw, bh, color=color,
            fc="white", lw=0.8, r=0.03, z=z + 1)
        txt(ax, (x + fx * w, (gy0 + gy1) / 2), text, size=5.2, color=color,
            z=z + 2)
    if label:
        txt(ax, (x + w / 2, y + nwires * dy - 0.10), label, size=6, color=ff.MIDGREY)


def thermometer(ax, xy, h, frac, color=ff.ORANGE, w=0.13, z=5):
    """Small vertical gauge, used to indicate temperature in phase cartoons."""
    x, y = xy
    box(ax, (x, y), w, h, color=ff.GREY, fc="white", lw=0.6, r=w / 2, z=z)
    ax.add_patch(Rectangle((x + 0.02, y + 0.02), w - 0.04, (h - 0.04) * frac,
                           facecolor=color, alpha=0.75, edgecolor="none",
                           zorder=z + 1))
