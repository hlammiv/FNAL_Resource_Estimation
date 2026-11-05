"""Glyph library for the Wilson Hall cover art.

Every glyph draws at unit scale -- it fits inside a box roughly 1 wide by
1 tall centered on `xy` -- and is scaled by `s`. Styling flows in from the
caller (`c` = primary line color, `lw` = line weight, `accent` = optional
accent color), so the blueprint palette is applied uniformly and can be
changed in one place in the composition script.

Glyphs never call plt directly and never set colors of their own beyond the
passed ones; that is what keeps the sheet coherent.
"""

import numpy as np
from matplotlib.patches import Circle, Ellipse, FancyArrowPatch, Polygon, Rectangle


# ---------------------------------------------------------------- helpers --

def _line(ax, xs, ys, c, lw, z=6, alpha=1.0, ls="-"):
    ax.plot(xs, ys, color=c, lw=lw, zorder=z, alpha=alpha, ls=ls,
            solid_capstyle="round")


def _rect(ax, x, y, w, h, c, lw, z=6, fc="none", alpha=1.0):
    ax.add_patch(Rectangle((x, y), w, h, facecolor=fc, edgecolor=c, lw=lw,
                           zorder=z, alpha=alpha))


def _circ(ax, xy, r, c, lw, z=6, fc="none", alpha=1.0):
    ax.add_patch(Circle(xy, r, facecolor=fc, edgecolor=c, lw=lw, zorder=z,
                        alpha=alpha))


def _txt(ax, xy, s_, size, c, z=7, alpha=1.0):
    ax.text(xy[0], xy[1], s_, fontsize=size, color=c, ha="center",
            va="center", zorder=z, alpha=alpha)


# --------------------------------------------------- tier 1: ground floor --

def transmon(ax, xy, s, c, lw, accent=None):
    """Xmon cross: the workhorse qubit. Single 12-point outline, no seam."""
    x, y = xy
    a = 0.42 * s
    w = 0.075 * s
    pts = np.array([
        (-a, w), (-w, w), (-w, a), (w, a), (w, w), (a, w),
        (a, -w), (w, -w), (w, -a), (-w, -a), (-w, -w), (-a, -w),
    ])
    ax.add_patch(Polygon(pts + np.array([x, y]), closed=True,
                         facecolor="none", edgecolor=c, lw=lw, zorder=6))


def jjunction(ax, xy, s, c, lw, accent=None):
    """Josephson junction: boxed X in a loop with leads."""
    x, y = xy
    b = 0.26 * s
    _rect(ax, x - b, y - b, 2 * b, 2 * b, c, lw)
    _line(ax, [x - b, x + b], [y - b, y + b], c, lw)
    _line(ax, [x - b, x + b], [y + b, y - b], c, lw)
    _line(ax, [x - 0.48 * s, x - b], [y, y], c, lw)
    _line(ax, [x + b, x + 0.48 * s], [y, y], c, lw)


def gate1q(ax, xy, s, c, lw, accent=None, label="H"):
    """Boxed single-qubit gate on a wire stub."""
    x, y = xy
    b = 0.30 * s
    _line(ax, [x - 0.50 * s, x - b], [y, y], c, lw)
    _line(ax, [x + b, x + 0.50 * s], [y, y], c, lw)
    _rect(ax, x - b, y - b, 2 * b, 2 * b, accent or c, lw, fc="#0d2c56")
    # Text size tracks physical glyph size: s is in axes units (inches in the
    # cover composition), so 26*s keeps the letter ~60% of the box height.
    _txt(ax, (x, y), label, 26 * s, accent or c)


def cnot(ax, xy, s, c, lw, accent=None):
    """Control dot and target circle-plus across two wire stubs."""
    x, y = xy
    dy = 0.30 * s
    r = 0.16 * s
    for yy in (y + dy, y - dy):
        _line(ax, [x - 0.48 * s, x + 0.48 * s], [yy, yy], c, lw)
    _circ(ax, (x, y + dy), 0.06 * s, c, lw, fc=c)
    _circ(ax, (x, y - dy), r, c, lw)
    _line(ax, [x, x], [y + dy, y - dy - r], c, lw)
    _line(ax, [x - r, x + r], [y - dy, y - dy], c, lw)


def cz(ax, xy, s, c, lw, accent=None):
    """Symmetric controlled-Z: two dots joined."""
    x, y = xy
    dy = 0.28 * s
    for yy in (y + dy, y - dy):
        _line(ax, [x - 0.42 * s, x + 0.42 * s], [yy, yy], c, lw)
        _circ(ax, (x, yy), 0.06 * s, c, lw, fc=c)
    _line(ax, [x, x], [y - dy, y + dy], c, lw)


def meter(ax, xy, s, c, lw, accent=None):
    """Measurement dial."""
    x, y = xy
    b = 0.34 * s
    _rect(ax, x - b, y - 0.26 * s, 2 * b, 0.52 * s, c, lw)
    th = np.linspace(0.25 * np.pi, 0.75 * np.pi, 30)
    _line(ax, x + 0.24 * s * np.cos(th), y - 0.14 * s + 0.24 * s * np.sin(th),
          c, lw)
    _line(ax, [x, x + 0.16 * s], [y - 0.12 * s, y + 0.12 * s], c, lw)


# ------------------------------------------------- tier 2: NISQ era, small --

def bloch(ax, xy, s, c, lw, accent=None):
    """Bloch sphere with a tilted state arrow."""
    x, y = xy
    r = 0.40 * s
    _circ(ax, (x, y), r, c, lw)
    ax.add_patch(Ellipse((x, y), 2 * r, 0.66 * r, facecolor="none",
                         edgecolor=c, lw=lw * 0.8, zorder=6, alpha=0.75))
    _line(ax, [x, x], [y - r, y + r], c, lw * 0.7, alpha=0.6)
    ar = FancyArrowPatch((x, y), (x + 0.55 * r, y + 0.62 * r),
                         arrowstyle="-|>", mutation_scale=4.5 * s,
                         color=accent or c, lw=lw, zorder=7)
    ax.add_patch(ar)


def ion_chain(ax, xy, s, c, lw, accent=None, n=5):
    """A linear trap with a chain of ions."""
    x, y = xy
    hw = 0.48 * s
    for dy in (0.16 * s, -0.16 * s):
        _line(ax, [x - hw, x + hw], [y + dy, y + dy], c, lw * 0.8)
    for i in range(n):
        xx = x - 0.36 * s + i * 0.72 * s / (n - 1)
        _circ(ax, (xx, y), 0.055 * s, c, lw, fc=c)


def chip5(ax, xy, s, c, lw, accent=None):
    """Five-qubit die, bow-tie coupling: the 2016-era cloud processor."""
    x, y = xy
    b = 0.44 * s
    _rect(ax, x - b, y - b, 2 * b, 2 * b, c, lw)
    pts = [(-0.22, 0.20), (0.22, 0.20), (0, 0), (-0.22, -0.20), (0.22, -0.20)]
    for (dx, dy) in pts:
        _circ(ax, (x + dx * s, y + dy * s), 0.055 * s, c, lw, fc=c)
    for (a_, b_) in [(0, 2), (1, 2), (3, 2), (4, 2)]:
        _line(ax, [x + pts[a_][0] * s, x + pts[b_][0] * s],
              [y + pts[a_][1] * s, y + pts[b_][1] * s], c, lw * 0.8)


def resonator(ax, xy, s, c, lw, accent=None):
    """Readout resonator meander."""
    x, y = xy
    t = np.linspace(0, 1, 200)
    xs = x - 0.45 * s + 0.9 * s * t
    ys = y + 0.16 * s * np.sign(np.sin(10 * np.pi * t)) * (np.abs(
        np.sin(10 * np.pi * t)) ** 0.15)
    _line(ax, xs, ys, c, lw * 0.85)


def circuit3(ax, xy, s, c, lw, accent=None, seed=0):
    """Three-wire circuit snippet with a few boxed gates."""
    x, y = xy
    rng = np.random.default_rng(seed)
    hw = 0.52 * s
    ys = [y - 0.26 * s, y, y + 0.26 * s]
    for yy in ys:
        _line(ax, [x - hw, x + hw], [yy, yy], c, lw * 0.8)
    for k in range(3):
        gx = x - 0.30 * s + k * 0.30 * s
        wi = int(rng.integers(0, 3))
        b = 0.085 * s
        _rect(ax, gx - b, ys[wi] - b, 2 * b, 2 * b, c, lw * 0.9,
              fc="#0d2c56")
    # one two-qubit link
    gx = x + 0.36 * s
    _circ(ax, (gx, ys[2]), 0.035 * s, c, lw, fc=c)
    _circ(ax, (gx, ys[0]), 0.06 * s, c, lw)
    _line(ax, [gx, gx], [ys[2], ys[0] - 0.06 * s], c, lw * 0.9)


# ------------------------------------------- tier 3: bigger integrated HW --

def heavy_hex(ax, xy, s, c, lw, accent=None):
    """Heavy-hex qubit patch."""
    x, y = xy
    r = 0.30 * s
    th = np.linspace(0, 2 * np.pi, 7) + np.pi / 6
    hx, hy = x + r * np.cos(th), y + r * np.sin(th)
    _line(ax, hx, hy, c, lw * 0.8)
    for (xx, yy) in zip(hx[:-1], hy[:-1]):
        _circ(ax, (xx, yy), 0.05 * s, c, lw, fc=c)
    for i in range(6):
        mx, my = (hx[i] + hx[i + 1]) / 2, (hy[i] + hy[i + 1]) / 2
        _circ(ax, (mx, my), 0.038 * s, c, lw)


def srf_cavity(ax, xy, s, c, lw, accent=None, cells=4):
    """Multicell elliptical SRF cavity in profile -- the SQMS substrate."""
    x, y = xy
    hw = 0.52 * s
    beam = 0.07 * s
    cw = 2 * hw / cells
    t = np.linspace(0, np.pi, 40)
    for i in range(cells):
        cx = x - hw + (i + 0.5) * cw
        xs = cx + 0.5 * cw * np.cos(t + np.pi)
        for sgn in (1, -1):
            ys = y + sgn * (beam + 0.20 * s * np.sin(t))
            _line(ax, xs, ys, c, lw * 0.85)
    _line(ax, [x - hw - 0.10 * s, x - hw], [y + beam, y + beam], c, lw * 0.85)
    _line(ax, [x - hw - 0.10 * s, x - hw], [y - beam, y - beam], c, lw * 0.85)
    _line(ax, [x + hw, x + hw + 0.10 * s], [y + beam, y + beam], c, lw * 0.85)
    _line(ax, [x + hw, x + hw + 0.10 * s], [y - beam, y - beam], c, lw * 0.85)


def tweezer_array(ax, xy, s, c, lw, accent=None):
    """Optical-tweezer atom array: dot grid under lens cones."""
    x, y = xy
    for i in range(4):
        for j in range(2):
            _circ(ax, (x - 0.33 * s + i * 0.22 * s, y - 0.10 * s + j * 0.22 * s),
                  0.045 * s, c, lw, fc=c)
    for i in (0, 2):
        ax0 = x - 0.33 * s + i * 0.22 * s
        _line(ax, [ax0 - 0.09 * s, ax0, ax0 + 0.09 * s],
              [y + 0.42 * s, y + 0.14 * s, y + 0.42 * s], c, lw * 0.7,
              alpha=0.8)


def qpu_die(ax, xy, s, c, lw, accent=None, n=4):
    """A denser many-qubit die: n x n grid with couplers."""
    x, y = xy
    b = 0.46 * s
    _rect(ax, x - b, y - b, 2 * b, 2 * b, c, lw)
    g = np.linspace(-0.30 * s, 0.30 * s, n)
    for gx in g:
        for gy in g:
            _circ(ax, (x + gx, y + gy), 0.035 * s, c, lw, fc=c)
    for gx in g:
        _line(ax, [x + gx, x + gx], [y - 0.30 * s, y + 0.30 * s], c, lw * 0.5,
              alpha=0.7)
        _line(ax, [x - 0.30 * s, x + 0.30 * s], [y + gx, y + gx], c, lw * 0.5,
              alpha=0.7)


# --------------------------------------------- tier 4: cryo + integration --

def dilfridge(ax, xy, s, c, lw, accent=None):
    """Dilution refrigerator, properly dressed: dashed vacuum jacket, pulse
    tube on the top flange, five narrowing plates with support rods, coax
    lines dropping to a chip mount at the mixing-chamber stage."""
    x, y = xy
    col = accent or c
    # vacuum jacket
    _rect(ax, x - 0.46 * s, y - 0.52 * s, 0.92 * s, 1.00 * s, c, lw * 0.6)
    ax.add_patch(Rectangle((x - 0.40 * s, y - 0.46 * s), 0.80 * s, 0.88 * s,
                           facecolor="none", edgecolor=c, lw=lw * 0.45,
                           linestyle=(0, (2, 2)), zorder=6, alpha=0.8))
    # top flange and pulse tube
    _line(ax, [x - 0.46 * s, x + 0.46 * s], [y + 0.48 * s] * 2, c, lw)
    _rect(ax, x - 0.10 * s, y + 0.48 * s, 0.20 * s, 0.10 * s, c, lw * 0.8)
    _line(ax, [x - 0.04 * s, x - 0.04 * s], [y + 0.58 * s, y + 0.64 * s], c,
          lw * 0.7)
    _line(ax, [x + 0.04 * s, x + 0.04 * s], [y + 0.58 * s, y + 0.64 * s], c,
          lw * 0.7)
    # five plates, narrowing toward the cold end
    widths = [0.38, 0.33, 0.28, 0.23, 0.18]
    ys = [y + 0.34 * s - i * 0.155 * s for i in range(5)]
    for w, yy in zip(widths, ys):
        _line(ax, [x - w * s, x + w * s], [yy, yy], col, lw * 1.25)
    # support rods between consecutive plates
    for w0, w1, y0, y1 in zip(widths[:-1], widths[1:], ys[:-1], ys[1:]):
        for sgn in (-1, 1):
            _line(ax, [x + sgn * 0.80 * w0 * s, x + sgn * 0.80 * w1 * s],
                  [y0, y1], c, lw * 0.55)
    # coax lines wiggling down the middle
    t = np.linspace(0, 1, 60)
    for k, xoff in enumerate((-0.07, 0.0, 0.07)):
        xs = x + xoff * s + 0.018 * s * np.sin(9 * np.pi * t + k)
        _line(ax, xs, y + 0.34 * s - t * (0.34 + 0.28) * s, c, lw * 0.4,
              alpha=0.9)
    # chip mount at the mixing chamber
    _rect(ax, x - 0.07 * s, y - 0.40 * s, 0.14 * s, 0.07 * s, c, lw * 0.7,
          fc=c)


def cryostat(ax, xy, s, c, lw, accent=None):
    """Cryostat can with nested shells."""
    x, y = xy
    for k, (w, h) in enumerate([(0.42, 0.48), (0.33, 0.39), (0.24, 0.30)]):
        _rect(ax, x - w * s, y - h * s, 2 * w * s, 2 * h * s, c,
              lw * (1.0 - 0.15 * k), alpha=1.0 - 0.2 * k)
    _line(ax, [x, x], [y + 0.48 * s, y + 0.62 * s], c, lw)


def rack(ax, xy, s, c, lw, accent=None, seed=0):
    """Equipment rack: HPC or control electronics."""
    x, y = xy
    w, h = 0.34 * s, 0.52 * s
    _rect(ax, x - w, y - h, 2 * w, 2 * h, c, lw)
    rng = np.random.default_rng(seed)
    n = 6
    for i in range(n):
        yy = y - h + (i + 0.5) * 2 * h / n
        _line(ax, [x - 0.78 * w, x + 0.78 * w], [yy, yy], c, lw * 0.6,
              alpha=0.85)
        if rng.uniform() < 0.6:
            _circ(ax, (x + 0.52 * w, yy + 0.035 * s), 0.018 * s,
                  accent or c, lw * 0.5, fc=accent or c)


# --------------------------------------------------------- tier 5: QEC ----

def surface_code(ax, xy, s, c, lw, accent=None, n=4):
    """Surface-code patch: data qubits on a grid, checkered plaquettes."""
    x, y = xy
    b = 0.42 * s
    d = 2 * b / n
    for i in range(n):
        for j in range(n):
            if (i + j) % 2 == 0:
                _rect(ax, x - b + i * d, y - b + j * d, d, d, c, lw * 0.5,
                      fc=c, alpha=0.22)
    for i in range(n + 1):
        _line(ax, [x - b + i * d] * 2, [y - b, y + b], c, lw * 0.6)
        _line(ax, [x - b, x + b], [y - b + i * d] * 2, c, lw * 0.6)
    for i in range(n + 1):
        for j in range(n + 1):
            _circ(ax, (x - b + i * d, y - b + j * d), 0.030 * s, c, lw * 0.5,
                  fc=c)


def lattice_surgery(ax, xy, s, c, lw, accent=None):
    """Two code patches with a merged seam."""
    x, y = xy
    b = 0.26 * s
    for dx in (-0.30 * s, 0.30 * s):
        _rect(ax, x + dx - b, y - b, 2 * b, 2 * b, c, lw * 0.8)
    _rect(ax, x - 0.06 * s, y - b, 0.12 * s, 2 * b, accent or c, lw * 0.8,
          fc=accent or c, alpha=0.30)
    _line(ax, [x - 0.04 * s, x - 0.04 * s], [y - b, y + b], accent or c,
          lw * 0.6, ls=(0, (2, 1.6)))
    _line(ax, [x + 0.04 * s, x + 0.04 * s], [y - b, y + b], accent or c,
          lw * 0.6, ls=(0, (2, 1.6)))


def logical_ket(ax, xy, s, c, lw, accent=None):
    """A logical qubit: |0_L> written on a small shielded patch."""
    x, y = xy
    b = 0.40 * s
    _rect(ax, x - b, y - b * 0.8, 2 * b, 1.6 * b, c, lw)
    _rect(ax, x - b * 0.82, y - b * 0.62, 1.64 * b, 1.24 * b, c, lw * 0.6,
          alpha=0.7)
    _txt(ax, (x, y), r"$|0_L\rangle$", 15 * s, accent or c)


# ------------------------------------------ gates & factories, added v2 ----

def toffoli(ax, xy, s, c, lw, accent=None):
    """Toffoli: controls on two wires, target plus on the third."""
    x, y = xy
    dy = 0.26 * s
    r = 0.13 * s
    for yy in (y + dy, y, y - dy):
        _line(ax, [x - 0.44 * s, x + 0.44 * s], [yy, yy], c, lw * 0.9)
    for yy in (y + dy, y):
        _circ(ax, (x, yy), 0.05 * s, c, lw, fc=c)
    _circ(ax, (x, y - dy), r, c, lw)
    _line(ax, [x, x], [y + dy, y - dy - r], c, lw * 0.9)
    _line(ax, [x - r, x + r], [y - dy, y - dy], c, lw * 0.9)


def swapg(ax, xy, s, c, lw, accent=None):
    """SWAP: crosses joined across two wires."""
    x, y = xy
    dy = 0.26 * s
    a = 0.09 * s
    for yy in (y + dy, y - dy):
        _line(ax, [x - 0.42 * s, x + 0.42 * s], [yy, yy], c, lw * 0.9)
        _line(ax, [x - a, x + a], [yy - a, yy + a], c, lw)
        _line(ax, [x - a, x + a], [yy + a, yy - a], c, lw)
    _line(ax, [x, x], [y - dy, y + dy], c, lw * 0.9)


def msd_factory(ax, xy, s, c, lw, accent=None):
    """Magic-state distillation, 15 raw T states funnelled to one clean one:
    three tiers of inputs converging through the triangle, the distilled
    state glowing at the tip, and a legible 15:1 label."""
    x, y = xy
    col = accent or c
    top, tip = y + 0.42 * s, y - 0.40 * s
    _line(ax, [x - 0.52 * s, x + 0.52 * s], [top, top], c, lw)
    _line(ax, [x - 0.52 * s, x], [top, tip], c, lw)
    _line(ax, [x + 0.52 * s, x], [top, tip], c, lw)
    # 15 raw T states: tiers of 7, 5, 3 descending into the funnel
    for (ny, cnt) in ((0.32, 7), (0.13, 5), (-0.06, 3)):
        hw = 0.52 * s * (ny * s and (top - (y + ny * s)) or 1)
        hw = 0.52 * s * (1 - (top - (y + ny * s)) / (top - tip)) \
            if False else 0.44 * s * (cnt / 7.0)
        for i in range(cnt):
            xx = x - hw + (2 * hw) * (i / (cnt - 1) if cnt > 1 else 0.5)
            _circ(ax, (xx, y + ny * s), 0.032 * s, c, lw * 0.5, fc=c,
                  alpha=0.9)
    # convergence arrows
    for sgn in (-1, 1):
        _line(ax, [x + sgn * 0.20 * s, x + sgn * 0.05 * s],
              [y - 0.16 * s, y - 0.30 * s], c, lw * 0.55, alpha=0.8)
    # the distilled magic state, glowing
    _circ(ax, (x, tip + 0.09 * s), 0.10 * s, col, lw * 0.6, fc="none",
          alpha=0.5)
    _circ(ax, (x, tip + 0.09 * s), 0.062 * s, col, lw * 0.8, fc=col)
    # label centered under the tip, on clear plate
    _txt(ax, (x, tip - 0.16 * s), r"$15\,{:}\,1$", 14 * s, col)


def group_circuit(ax, xy, s, c, lw, accent=None):
    """Register-level primitive-gate circuit: slashed register wires, a
    Fourier box on one register and a group-multiplication box spanning
    both."""
    x, y = xy
    dy = 0.22 * s
    for yy in (y + dy, y - dy):
        _line(ax, [x - 0.50 * s, x + 0.50 * s], [yy, yy], c, lw * 0.9)
        _line(ax, [x - 0.42 * s, x - 0.36 * s], [yy - 0.05 * s, yy + 0.05 * s],
              c, lw * 0.8)
    b = 0.155 * s
    _rect(ax, x - 0.26 * s - b, y + dy - b, 2 * b, 2 * b, c, lw * 0.9,
          fc="#0d2c56")
    _txt(ax, (x - 0.26 * s, y + dy), r"$U_F$", 10 * s, accent or c)
    bx = 0.17 * s
    _rect(ax, x + 0.15 * s - bx, y - dy - b, 2 * bx, 2 * dy + 2 * b,
          accent or c, lw * 0.9, fc="#0d2c56")
    _txt(ax, (x + 0.15 * s, y), r"$U_{\times}$", 10 * s, accent or c)


# ---------------------------------- hardware & abstractions, added v3 ----

def squid(ax, xy, s, c, lw, accent=None):
    """rf-SQUID: a loop interrupted by two tiny junctions."""
    x, y = xy
    r = 0.30 * s
    th = np.linspace(0, 2 * np.pi, 60)
    _line(ax, x + r * np.cos(th), y + r * np.sin(th), c, lw * 0.9)
    for sgn in (-1, 1):
        bx = x + sgn * r
        b = 0.085 * s
        _rect(ax, bx - b, y - b, 2 * b, 2 * b, c, lw * 0.8, fc="#0d2c56")
        _line(ax, [bx - b, bx + b], [y - b, y + b], c, lw * 0.7)
        _line(ax, [bx - b, bx + b], [y + b, y - b], c, lw * 0.7)


def spin_qubit(ax, xy, s, c, lw, accent=None):
    """Spin in a dot: circle with a spin-up arrow."""
    x, y = xy
    _circ(ax, (x, y), 0.30 * s, c, lw * 0.9)
    ar = FancyArrowPatch((x, y - 0.20 * s), (x, y + 0.26 * s),
                         arrowstyle="-|>", mutation_scale=4.5 * s,
                         color=accent or c, lw=lw, zorder=7)
    ax.add_patch(ar)


def blade_trap(ax, xy, s, c, lw, accent=None):
    """A single trapped ion between blade electrodes."""
    x, y = xy
    for sgn in (-1, 1):
        pts = np.array([(-0.34, 0.38 * sgn), (0.34, 0.38 * sgn),
                        (0.12, 0.14 * sgn), (-0.12, 0.14 * sgn)]) * s
        ax.add_patch(Polygon(pts + np.array([x, y]), closed=True,
                             facecolor="none", edgecolor=c, lw=lw * 0.85,
                             zorder=6))
    _circ(ax, (x, y), 0.055 * s, accent or c, lw, fc=accent or c)


def cold_atoms(ax, xy, s, c, lw, accent=None):
    """A cold-atom cloud at the crossing of two dashed beams (a MOT)."""
    x, y = xy
    for a in (np.pi / 4, 3 * np.pi / 4):
        dx, dy = 0.46 * s * np.cos(a), 0.46 * s * np.sin(a)
        _line(ax, [x - dx, x + dx], [y - dy, y + dy], c, lw * 0.6,
              alpha=0.8, ls=(0, (2, 1.6)))
    rng = np.random.default_rng(4)
    for _ in range(7):
        a, r = rng.uniform(0, 2 * np.pi), rng.uniform(0, 0.13 * s)
        _circ(ax, (x + r * np.cos(a), y + r * np.sin(a)), 0.035 * s, c,
              lw * 0.5, fc=c)


def tensor_net(ax, xy, s, c, lw, accent=None):
    """Matrix-product state: bonded tensors with physical legs."""
    x, y = xy
    n = 4
    xs = [x - 0.36 * s + i * 0.72 * s / (n - 1) for i in range(n)]
    _line(ax, [xs[0], xs[-1]], [y + 0.08 * s, y + 0.08 * s], c, lw * 0.8)
    for xx in xs:
        _circ(ax, (xx, y + 0.08 * s), 0.075 * s, c, lw * 0.85, fc="#0d2c56")
        _line(ax, [xx, xx], [y - 0.01 * s, y - 0.22 * s], c, lw * 0.7)


def plaquette(ax, xy, s, c, lw, accent=None):
    """One stabilizer check: four data qubits around an ancilla star."""
    x, y = xy
    d = 0.30 * s
    corners = [(x - d, y - d), (x + d, y - d), (x + d, y + d), (x - d, y + d)]
    for (cx_, cy_) in corners:
        _line(ax, [x, cx_], [y, cy_], c, lw * 0.6, alpha=0.9)
        _circ(ax, (cx_, cy_), 0.055 * s, c, lw * 0.7, fc=c)
    _circ(ax, (x, y), 0.075 * s, accent or c, lw * 0.8, fc="#0d2c56")


def braid(ax, xy, s, c, lw, accent=None):
    """Anyon worldlines exchanging: a braid with an over-crossing."""
    x, y = xy
    t = np.linspace(0, 1, 60)
    yy = y - 0.38 * s + 0.76 * s * t
    xa = x - 0.16 * s + 0.32 * s * (3 * t ** 2 - 2 * t ** 3)
    xb = x + 0.16 * s - 0.32 * s * (3 * t ** 2 - 2 * t ** 3)
    _line(ax, xa, yy, c, lw * 0.85)
    mask = np.abs(t - 0.5) > 0.13          # break the under-strand
    _line(ax, np.where(mask, xb, np.nan), yy, c, lw * 0.85)
    _line(ax, [x - 0.34 * s, x - 0.34 * s], [y - 0.38 * s, y + 0.38 * s], c,
          lw * 0.85)


# --------------------------------------- hardware additions, v4 ----------

def wirebond_chip(ax, xy, s, c, lw, accent=None):
    """Die in a package, wirebonds fanning to the frame."""
    x, y = xy
    B = 0.44 * s
    _rect(ax, x - B, y - B, 2 * B, 2 * B, c, lw * 0.9)
    d = 0.20 * s
    _rect(ax, x - d, y - d, 2 * d, 2 * d, c, lw * 0.8, fc="#0d2c56")
    for k in range(3):
        f = -0.6 + 0.6 * k
        _line(ax, [x - d, x - B], [y + f * d, y + f * B * 0.9], c, lw * 0.45)
        _line(ax, [x + d, x + B], [y + f * d, y + f * B * 0.9], c, lw * 0.45)
        _line(ax, [x + f * d, x + f * B * 0.9], [y + d, y + B], c, lw * 0.45)
        _line(ax, [x + f * d, x + f * B * 0.9], [y - d, y - B], c, lw * 0.45)


def hemt_amp(ax, xy, s, c, lw, accent=None):
    """Cryogenic amplifier: the triangle symbol on a signal line."""
    x, y = xy
    _line(ax, [x - 0.48 * s, x - 0.20 * s], [y, y], c, lw * 0.9)
    _line(ax, [x + 0.20 * s, x + 0.48 * s], [y, y], c, lw * 0.9)
    pts = np.array([(-0.20, 0.26), (-0.20, -0.26), (0.20, 0.0)]) * s
    ax.add_patch(Polygon(pts + np.array([x, y]), closed=True,
                         facecolor="#0d2c56", edgecolor=accent or c, lw=lw,
                         zorder=6))
    _txt(ax, (x - 0.05 * s, y), "+", 13 * s, accent or c)


def circulator(ax, xy, s, c, lw, accent=None):
    """Microwave circulator: circle, arced arrow, three ports."""
    x, y = xy
    r = 0.28 * s
    _circ(ax, (x, y), r, c, lw * 0.9)
    th = np.linspace(0.25 * np.pi, 1.45 * np.pi, 40)
    _line(ax, x + 0.15 * s * np.cos(th), y + 0.15 * s * np.sin(th), c,
          lw * 0.7)
    ar = FancyArrowPatch((x + 0.15 * s * np.cos(th[-2]), y + 0.15 * s * np.sin(th[-2])),
                         (x + 0.15 * s * np.cos(th[-1] + 0.3),
                          y + 0.15 * s * np.sin(th[-1] + 0.3)),
                         arrowstyle="-|>", mutation_scale=3.5 * s, color=c,
                         lw=lw * 0.7, zorder=7)
    ax.add_patch(ar)
    for a in (np.pi / 2, np.pi * 7 / 6, np.pi * 11 / 6):
        _line(ax, [x + r * np.cos(a), x + (r + 0.16 * s) * np.cos(a)],
              [y + r * np.sin(a), y + (r + 0.16 * s) * np.sin(a)], c,
              lw * 0.8)


def laser(ax, xy, s, c, lw, accent=None):
    """Laser head, beam, and focusing lens -- the atoms' best friend."""
    x, y = xy
    _rect(ax, x - 0.46 * s, y - 0.13 * s, 0.34 * s, 0.26 * s, c, lw * 0.9)
    col = accent or c
    _line(ax, [x - 0.12 * s, x + 0.20 * s], [y, y], col, lw * 0.8)
    ax.add_patch(Ellipse((x + 0.24 * s, y), 0.055 * s, 0.30 * s,
                         facecolor="none", edgecolor=c, lw=lw * 0.8,
                         zorder=6))
    _line(ax, [x + 0.27 * s, x + 0.46 * s], [y + 0.10 * s, y + 0.015 * s], col,
          lw * 0.6)
    _line(ax, [x + 0.27 * s, x + 0.46 * s], [y - 0.10 * s, y - 0.015 * s], col,
          lw * 0.6)


def color_code(ax, xy, s, c, lw, accent=None):
    """Triangular color-code patch: hexagon split into three shaded cells."""
    x, y = xy
    r = 0.36 * s
    th = np.linspace(0, 2 * np.pi, 7) + np.pi / 6
    hx, hy = x + r * np.cos(th), y + r * np.sin(th)
    _line(ax, hx, hy, c, lw * 0.9)
    for k, a0 in enumerate((0, 2, 4)):
        pts = [(x, y), (hx[a0], hy[a0]), (hx[a0 + 1], hy[a0 + 1])]
        ax.add_patch(Polygon(pts, closed=True, facecolor=(accent or c),
                             alpha=0.10 + 0.09 * k, edgecolor="none",
                             zorder=5))
    for k in range(6):
        _line(ax, [x, hx[k]], [y, hy[k]], c, lw * 0.5, alpha=0.65)
    for k in range(6):
        _circ(ax, (hx[k], hy[k]), 0.045 * s, c, lw * 0.6, fc=c)


# ------------------------- physics content additions, v5 -----------------

def atom(ax, xy, s, c, lw, accent=None):
    """Bohr atom: nucleus, two tilted orbits, two electrons."""
    x, y = xy
    _circ(ax, (x, y), 0.055 * s, accent or c, lw, fc=accent or c)
    for ang in (25, -25):
        e = Ellipse((x, y), 0.92 * s, 0.34 * s, angle=ang, facecolor="none",
                    edgecolor=c, lw=lw * 0.75, zorder=6)
        ax.add_patch(e)
    for (ang, ph) in ((25, 0.3), (-25, 2.6)):
        a = np.deg2rad(ang)
        ex = 0.46 * s * np.cos(ph)
        ey = 0.17 * s * np.sin(ph)
        rx = ex * np.cos(a) - ey * np.sin(a)
        ry = ex * np.sin(a) + ey * np.cos(a)
        _circ(ax, (x + rx, y + ry), 0.035 * s, c, lw * 0.6, fc=c)


def molecule(ax, xy, s, c, lw, accent=None):
    """Ball-and-stick molecule: center atom, three neighbors, one double
    bond -- the quantum-chemistry payload."""
    x, y = xy
    ctr = (x - 0.06 * s, y)
    others = [(x + 0.34 * s, y + 0.26 * s), (x + 0.30 * s, y - 0.30 * s),
              (x - 0.42 * s, y - 0.06 * s)]
    _line(ax, [ctr[0], others[0][0]], [ctr[1], others[0][1]], c, lw * 0.8)
    _line(ax, [ctr[0] + 0.02 * s, others[1][0] + 0.02 * s],
          [ctr[1] - 0.03 * s, others[1][1] + 0.03 * s], c, lw * 0.7)
    _line(ax, [ctr[0] - 0.02 * s, others[1][0] - 0.02 * s],
          [ctr[1] + 0.01 * s, others[1][1] + 0.07 * s], c, lw * 0.7)
    _line(ax, [ctr[0], others[2][0]], [ctr[1], others[2][1]], c, lw * 0.8)
    _circ(ax, ctr, 0.10 * s, accent or c, lw * 0.8, fc="#0d2c56")
    for o in others:
        _circ(ax, o, 0.07 * s, c, lw * 0.7, fc=c)


def circuit5(ax, xy, s, c, lw, accent=None, seed=1):
    """A deeper circuit: four wires, entangling ladder, rotations, readout."""
    x, y = xy
    rng = np.random.default_rng(seed)
    hw = 0.62 * s
    ys = [y + (i - 1.5) * 0.20 * s for i in range(4)]
    for yy in ys:
        _line(ax, [x - hw, x + hw], [yy, yy], c, lw * 0.75)
    # H column
    for i in (0, 2):
        b = 0.07 * s
        _rect(ax, x - hw + 0.06 * s, ys[i] - b, 2 * b, 2 * b, c, lw * 0.7,
              fc="#0d2c56")
    # CNOT ladder
    for k in range(3):
        gx = x - 0.22 * s + k * 0.16 * s
        _circ(ax, (gx, ys[k]), 0.028 * s, c, lw, fc=c)
        _circ(ax, (gx, ys[k + 1]), 0.05 * s, c, lw * 0.7)
        _line(ax, [gx, gx], [ys[k], ys[k + 1] - 0.05 * s], c, lw * 0.7)
    # rotation boxes
    for i in (1, 3):
        b = 0.07 * s
        _rect(ax, x + 0.22 * s - b, ys[i] - b, 2 * b, 2 * b, accent or c,
              lw * 0.7, fc="#0d2c56")
    # measurement meters at line ends
    for i in range(4):
        b = 0.055 * s
        _rect(ax, x + hw - 2.4 * b, ys[i] - b, 2.4 * b, 2 * b, c, lw * 0.6,
              fc="#0d2c56")


def pulse_seq(ax, xy, s, c, lw, accent=None):
    """Control pulses on two channels: Gaussian, flat-top, and a DRAG wiggle."""
    x, y = xy
    for k, yy in enumerate((y + 0.16 * s, y - 0.22 * s)):
        _line(ax, [x - 0.55 * s, x + 0.55 * s], [yy, yy], c, lw * 0.5,
              alpha=0.75)
    t = np.linspace(-1, 1, 120)
    g = np.exp(-((t + 0.55) / 0.16) ** 2)
    flat = np.clip(1.15 * np.exp(-((t - 0.05) / 0.30) ** 4), 0, 1)
    _line(ax, x + t * 0.55 * s, y + 0.16 * s + 0.24 * s * g, accent or c,
          lw * 0.8)
    _line(ax, x + t * 0.55 * s, y + 0.16 * s + 0.24 * s * flat, c, lw * 0.7)
    drag = np.exp(-((t - 0.45) / 0.14) ** 2) * np.sin(18 * t)
    _line(ax, x + t * 0.55 * s, y - 0.22 * s + 0.16 * s * drag, c, lw * 0.7)


BRAKETS = [r"$\langle\psi|\hat{H}|\psi\rangle$",
           r"$\frac{|00\rangle+|11\rangle}{\sqrt{2}}$",
           r"$\langle\phi|\psi\rangle$",
           r"$|+\rangle$",
           r"$\frac{|000\rangle+|111\rangle}{\sqrt{2}}$"]


def braket(ax, xy, s, c, lw, accent=None, seed=0):
    """Dirac notation, floating free."""
    _txt(ax, xy, BRAKETS[seed % len(BRAKETS)], 13.5 * s, accent or c)


def shor_eq(ax, xy, s, c, lw, accent=None):
    """Shor's point, stacked so the plaque stays narrow: the semiprime, its
    factors, and the order relation that finds them."""
    x, y = xy
    col = accent or c
    _txt(ax, (x, y + 0.26 * s), r"$4\,088\,459$", 12 * s, col)
    _txt(ax, (x, y - 0.02 * s), r"$= 2017 \times 2027$", 10.5 * s, col)
    _txt(ax, (x, y - 0.30 * s), r"$a^{r}\equiv 1\ (\mathrm{mod}\ N)$",
         8 * s, c)


def qcd_ham(ax, xy, s, c, lw, accent=None):
    """The pure-gauge Kogut-Susskind Hamiltonian: what the stack is for."""
    x, y = xy
    col = accent or c
    _txt(ax, (x, y + 0.17 * s),
         r"$H_{\mathrm{KS}} = \frac{g^{2}}{2}\sum_{\ell} E_{\ell}^{2}$",
         11 * s, col)
    _txt(ax, (x, y - 0.21 * s),
         r"$-\,\frac{1}{2g^{2}}\sum_{p}\left(\mathrm{Tr}\,U_{p}+\mathrm{h.c.}\right)$",
         10 * s, col)


# ------------------------------- theory glyphs, v6 -----------------------

def chiral_scurve(ax, xy, s, c, lw, accent=None):
    """Domain-wall mass profile: the tanh kink, with the chiral zero mode
    bound where it crosses."""
    x, y = xy
    _line(ax, [x - 0.50 * s, x + 0.50 * s], [y, y], c, lw * 0.5, alpha=0.7)
    t = np.linspace(-1, 1, 100)
    _line(ax, x + 0.48 * s * t, y + 0.26 * s * np.tanh(3.2 * t), c, lw * 0.9)
    zm = np.exp(-((t / 0.28) ** 2))
    _line(ax, x + 0.48 * s * t, y + 0.20 * s * zm, accent or c, lw * 0.75)
    _txt(ax, (x + 0.36 * s, y + 0.30 * s), r"$m(s)$", 7 * s, c)


def sho(ax, xy, s, c, lw, accent=None):
    """Harmonic-oscillator potential with its evenly spaced levels."""
    x, y = xy
    t = np.linspace(-1, 1, 80)
    _line(ax, x + 0.46 * s * t, y - 0.30 * s + 0.72 * s * t * t, c, lw * 0.9)
    for n in range(3):
        e = 0.13 * s + n * 0.20 * s
        hw = 0.46 * s * np.sqrt(e / (0.72 * s))
        _line(ax, [x - hw, x + hw], [y - 0.30 * s + e, y - 0.30 * s + e],
              accent or c, lw * 0.6, alpha=0.9)
    _circ(ax, (x, y - 0.30 * s + 0.13 * s), 0.035 * s, accent or c, lw,
          fc=accent or c)


def qldpc(ax, xy, s, c, lw, accent=None):
    """QLDPC Tanner graph: sparse checks with long-range edges."""
    x, y = xy
    vx = [x + (i - 2) * 0.23 * s for i in range(5)]
    cx_ = [x + (i - 1) * 0.32 * s for i in range(3)]
    vy, cy = y - 0.28 * s, y + 0.28 * s
    edges = [(0, 0), (1, 0), (1, 1), (2, 1), (3, 1), (3, 2), (4, 2),
             (0, 2), (4, 0)]                    # last two are long-range
    for k, (i, j) in enumerate(edges):
        _line(ax, [vx[i], cx_[j]], [vy, cy], c,
              lw * (0.75 if k < 7 else 0.55),
              alpha=0.95 if k < 7 else 0.75)
    for xx in vx:
        _circ(ax, (xx, vy), 0.055 * s, c, lw * 0.7, fc=c)
    for k, xx in enumerate(cx_):
        b = 0.06 * s
        _rect(ax, xx - b, cy - b, 2 * b, 2 * b, accent or c, lw * 0.8,
              fc="#0d2c56")


def kubo(ax, xy, s, c, lw, accent=None):
    """Transport from the Kubo limit, and the KSS bound under it."""
    x, y = xy
    col = accent or c
    _txt(ax, (x, y + 0.19 * s),
         r"$\eta = \lim_{\omega\to 0}\frac{\rho_{xy}(\omega)}{\omega}$",
         13 * s, col)
    _txt(ax, (x, y - 0.24 * s), r"$\eta/s \geq \frac{1}{4\pi}$",
         14 * s, c)


def peps(ax, xy, s, c, lw, accent=None):
    """2D tensor network (PEPS): a 3x2 patch with bonds and physical legs."""
    x, y = xy
    xs = [x - 0.30 * s, x, x + 0.30 * s]
    ys_ = [y + 0.14 * s, y - 0.14 * s]
    for yy in ys_:
        _line(ax, [xs[0], xs[-1]], [yy, yy], c, lw * 0.75)
    for xx in xs:
        _line(ax, [xx, xx], [ys_[1], ys_[0]], c, lw * 0.75)
    for xx in xs:
        for yy in ys_:
            _circ(ax, (xx, yy), 0.065 * s, c, lw * 0.8, fc="#0d2c56")
            _line(ax, [xx + 0.05 * s, xx + 0.11 * s],
                  [yy - 0.05 * s, yy - 0.11 * s], c, lw * 0.5, alpha=0.8)


def mera(ax, xy, s, c, lw, accent=None):
    """MERA: the renormalization pyramid -- four sites coarse-grained to one."""
    x, y = xy
    b = [x + (i - 1.5) * 0.24 * s for i in range(4)]
    m = [x - 0.18 * s, x + 0.18 * s]
    yb, ym, yt = y - 0.30 * s, y + 0.02 * s, y + 0.32 * s
    for i, xx in enumerate(b):
        _line(ax, [xx, m[i // 2]], [yb, ym], c, lw * 0.7)
    for xx in m:
        _line(ax, [xx, x], [ym, yt], c, lw * 0.7)
    for xx in b:
        _circ(ax, (xx, yb), 0.045 * s, c, lw * 0.6, fc=c)
    for xx in m:
        _rect(ax, xx - 0.055 * s, ym - 0.055 * s, 0.11 * s, 0.11 * s, c,
              lw * 0.7, fc="#0d2c56")
    _circ(ax, (x, yt), 0.06 * s, accent or c, lw * 0.8, fc=accent or c)


def wilson_loop(ax, xy, s, c, lw, accent=None):
    """A Wilson loop: closed oriented path on a faint lattice patch."""
    x, y = xy
    for i in range(4):
        xx = x - 0.36 * s + i * 0.24 * s
        _line(ax, [xx, xx], [y - 0.28 * s, y + 0.28 * s], c, lw * 0.35,
              alpha=0.5)
    for j in range(3):
        yy = y - 0.28 * s + j * 0.28 * s
        _line(ax, [x - 0.36 * s, x + 0.36 * s], [yy, yy], c, lw * 0.35,
              alpha=0.5)
    col = accent or c
    x0, x1 = x - 0.24 * s, x + 0.24 * s
    y0, y1 = y - 0.14 * s, y + 0.14 * s
    for (xa, ya, xb, yb) in ((x0, y0, x1, y0), (x1, y0, x1, y1),
                             (x1, y1, x0, y1), (x0, y1, x0, y0)):
        _line(ax, [xa, xb], [ya, yb], col, lw * 1.1)
        mx, my = (xa + xb) / 2, (ya + yb) / 2
        dx, dy = (xb - xa), (yb - ya)
        n = max(abs(dx), abs(dy)) or 1
        ar = FancyArrowPatch((mx - 0.03 * s * dx / n, my - 0.03 * s * dy / n),
                             (mx + 0.03 * s * dx / n, my + 0.03 * s * dy / n),
                             arrowstyle="-|>", mutation_scale=3.6 * s,
                             color=col, lw=0, zorder=7)
        ax.add_patch(ar)


def trotter_block(ax, xy, s, c, lw, accent=None):
    """One Trotter step in a dashed block, repeated x N."""
    x, y = xy
    ys = [y - 0.20 * s, y, y + 0.20 * s]
    for yy in ys:
        _line(ax, [x - 0.55 * s, x + 0.48 * s], [yy, yy], c, lw * 0.75)
    ax.add_patch(Rectangle((x - 0.34 * s, y - 0.31 * s), 0.62 * s, 0.62 * s,
                           facecolor="none", edgecolor=accent or c,
                           lw=lw * 0.8, linestyle=(0, (2.4, 1.8)), zorder=6))
    b = 0.065 * s
    _rect(ax, x - 0.24 * s - b, ys[2] - b, 2 * b, 2 * b, c, lw * 0.7,
          fc="#0d2c56")
    _rect(ax, x + 0.02 * s - b, ys[1] - b, 2 * b, 2 * b, c, lw * 0.7,
          fc="#0d2c56")
    _circ(ax, (x + 0.18 * s, ys[0]), 0.028 * s, c, lw, fc=c)
    _circ(ax, (x + 0.18 * s, ys[1]), 0.05 * s, c, lw * 0.7)
    _line(ax, [x + 0.18 * s, x + 0.18 * s], [ys[0], ys[1] - 0.05 * s], c,
          lw * 0.7)
    _txt(ax, (x + 0.42 * s, y + 0.36 * s), r"$\times N$", 9 * s,
         accent or c)


def ghz(ax, xy, s, c, lw, accent=None):
    """The three-qubit GHZ state."""
    _txt(ax, xy, BRAKETS[4], 12.5 * s, accent or c)


CODES = [r"$[[7,1,3]]$", r"$[[5,1,3]]$", r"$[[144,12,12]]$",
         r"$[[4,2,2]]$", r"$[[9,1,3]]$", r"$[[15,1,3]]$",
         r"$[[23,1,7]]$", r"$[[42,14,6]]_{64}$", r"$[[n,k,d]]_{q}$"]


def code_params(ax, xy, s, c, lw, accent=None, seed=0):
    """QEC code parameters in double-bracket notation."""
    _txt(ax, xy, CODES[seed % len(CODES)], 14 * s, accent or c)


def mpo(ax, xy, s, c, lw, accent=None):
    """Matrix-product operator: tensors threaded on a horizontal bond with
    physical legs above and below."""
    x, y = xy
    _line(ax, [x - 0.46 * s, x + 0.46 * s], [y, y], c, lw * 0.8)
    xs = [x - 0.33 * s + i * 0.22 * s for i in range(4)]
    for xx in xs:
        _rect(ax, xx - 0.055 * s, y - 0.055 * s, 0.11 * s, 0.11 * s, c,
              lw * 0.8, fc="#0d2c56")
        _line(ax, [xx, xx], [y + 0.055 * s, y + 0.20 * s], c, lw * 0.55)
        _line(ax, [xx, xx], [y - 0.055 * s, y - 0.20 * s], c, lw * 0.55)


def tn_loop(ax, xy, s, c, lw, accent=None):
    """A traced contraction: the bond leaves the chain's last tensor and
    loops back over the top to the first -- Tr(ABC)."""
    x, y = xy
    _line(ax, [x - 0.40 * s, x + 0.40 * s], [y - 0.10 * s, y - 0.10 * s], c,
          lw * 0.8)
    th = np.linspace(0, np.pi, 60)
    ax.plot(x + 0.40 * s * np.cos(th), y - 0.10 * s + 0.34 * s * np.sin(th),
            color=accent or c, lw=lw * 0.7, zorder=6, alpha=0.9)
    for i in range(3):
        xx = x - 0.26 * s + i * 0.26 * s
        _circ(ax, (xx, y - 0.10 * s), 0.065 * s, c, lw * 0.8, fc="#0d2c56")


NU_KETS = [r"$|\nu_e\rangle$", r"$|\nu_\mu\rangle$",
           r"$|\nu_\tau\rangle$"]


def nu_ket(ax, xy, s, c, lw, accent=None, seed=0):
    """A neutrino flavor ket."""
    _txt(ax, xy, NU_KETS[seed % len(NU_KETS)], 13 * s, accent or c)


def nu_osc(ax, xy, s, c, lw, accent=None):
    """Flavor mixing through the PMNS matrix."""
    _txt(ax, xy,
         r"$|\nu_{\alpha}\rangle = \sum_{i} U_{\alpha i}\,|\nu_{i}\rangle$",
         12 * s, accent or c)


def u_nu(ax, xy, s, c, lw, accent=None):
    """The mixing matrix acting on a flavor ket."""
    _txt(ax, xy, r"$U_{\mathrm{PMNS}}\,|\nu_{\mu}\rangle$", 11.5 * s,
         accent or c)


# |U_alpha_i| magnitudes, roughly the measured PMNS moduli
_PMNS = [[0.82, 0.55, 0.15], [0.37, 0.57, 0.70], [0.39, 0.59, 0.69]]


def pmns_matrix(ax, xy, s, c, lw, accent=None):
    """The PMNS matrix, drawn: a 3x3 of dots sized by |U|, applied to the
    flavor column vector."""
    x, y = xy
    mx = x - 0.22 * s
    # matrix brackets
    for sgn, bx in ((1, mx - 0.34 * s), (-1, mx + 0.34 * s)):
        _line(ax, [bx + sgn * 0.06 * s, bx, bx, bx + sgn * 0.06 * s],
              [y + 0.30 * s, y + 0.30 * s, y - 0.30 * s, y - 0.30 * s], c,
              lw * 0.8)
    for i in range(3):
        for j in range(3):
            r = 0.028 * s + 0.052 * s * _PMNS[i][j]
            _circ(ax, (mx + (j - 1) * 0.20 * s, y + (1 - i) * 0.19 * s), r,
                  accent or c, lw * 0.4, fc=accent or c, alpha=0.9)
    # flavor column
    vx = x + 0.42 * s
    for sgn, bx in ((1, vx - 0.16 * s), (-1, vx + 0.16 * s)):
        _line(ax, [bx + sgn * 0.05 * s, bx, bx, bx + sgn * 0.05 * s],
              [y + 0.30 * s, y + 0.30 * s, y - 0.30 * s, y - 0.30 * s], c,
              lw * 0.8)
    for k, lab in enumerate((r"$\nu_e$", r"$\nu_\mu$", r"$\nu_\tau$")):
        _txt(ax, (vx, y + (1 - k) * 0.19 * s), lab, 7.5 * s, c)


# ------------------------------------------------------------- scenery ----

def bison(ax, xy, s, c, alpha=0.9, z=6):
    """Grazing bison silhouette, facing left.

    Identity comes from three cues: the tall hump over the FRONT shoulders,
    the head hanging low off its front face, and the slighter hindquarters.
    """
    x, y = xy
    pts = np.array([
        # nose, low and forward
        (-0.58, -0.10), (-0.60, -0.02),
        # forehead rising steeply into the hump
        (-0.55, 0.06), (-0.48, 0.14), (-0.40, 0.26), (-0.30, 0.38),
        # hump peak over the front shoulders
        (-0.18, 0.44), (-0.06, 0.44),
        # back sloping DOWN toward the rear
        (0.10, 0.34), (0.24, 0.26), (0.36, 0.20),
        # rump and tail nub
        (0.46, 0.16), (0.50, 0.10), (0.53, 0.12), (0.52, 0.06),
        # rear legs
        (0.48, -0.02), (0.46, -0.24), (0.40, -0.24), (0.39, -0.04),
        (0.20, -0.04), (0.18, -0.24), (0.12, -0.24), (0.11, -0.04),
        # belly
        (-0.06, -0.06),
        # front legs, under the hump
        (-0.12, -0.06), (-0.14, -0.24), (-0.20, -0.24), (-0.21, -0.06),
        (-0.34, -0.06), (-0.36, -0.24), (-0.42, -0.24), (-0.43, -0.08),
        # beard hanging from the low head, back up to the nose
        (-0.47, -0.10), (-0.50, -0.22), (-0.54, -0.20), (-0.56, -0.12),
    ]) * s
    ax.add_patch(Polygon(pts + np.array([x, y]), closed=True, facecolor=c,
                         edgecolor="none", alpha=alpha, zorder=z))
