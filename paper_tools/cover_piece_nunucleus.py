#!/usr/bin/env python3
"""Detailed neutrino--nucleus scattering ornament (DUNE flagship vignette).

Textbook-quality CC event with full intranuclear cascade, in the cover's
blueprint idiom: a dashed neutrino arrives at the leptonic vertex, the
charged lepton bends away arrowed, and the W (wavy, ORANGE) dives into a
packed two-shade nucleus of ~19 nucleons to strike one (ORANGE).  The
struck nucleon rescatters off two others with visible deflection vertices,
a produced pion escapes as a fine chevron-ticked track, knocked-out
nucleons leave with arrowed stubs, and a faint dashed "transparency" ghost
continues straight through the far side.  Fine touches: momentum-transfer
arrow q, Fermi-motion arrows inside spectators, residual-nucleus recoil
arrow, dashed nuclear boundary, vertex flash ticks.

Exposes draw_nu_nucleus(ax, center, rot_deg, scale, alpha, z): all
geometry lives in an origin-centered design frame and is pushed through
rotate(rot_deg) -> scale -> translate(center) at data level.  Fills stay
at alpha <= 0.15; every alpha is multiplied by the alpha parameter.
Deterministic (fixed seed).
"""

import numpy as np
from matplotlib.patches import Circle

# ------------------------------------------------------ self-contained palette
WHITE = "#f2f6fb"
PALE = "#9fc4e8"
ORANGE = "#e8833a"

DESIGN_RADIUS = 0.82          # max half-extent of the design (data units)

# ------------------------------------------------- design-frame geometry
_NC = np.array([0.00, -0.10])   # nucleus center
_RB = 0.345                     # dashed nuclear boundary radius
_RPACK = 0.235                  # nucleon-center packing radius
_RN = 0.056                     # nucleon radius

_V = np.array([-0.42, 0.18])    # leptonic vertex (nu -> lepton + W)
_NU0 = np.array([-0.74, 0.35])  # neutrino entry point
_S = np.array([-0.10, -0.04])   # struck nucleon
_N1 = np.array([0.08, -0.16])   # first rescatter partner
_N2 = np.array([0.22, -0.06])   # second rescatter partner
_LE = np.array([0.61, 0.54])    # lepton track end


def _bez(p0, p1, bow, n=90):
    """Quadratic bezier p0->p1, control at midpoint + perp * bow."""
    p0, p1 = np.asarray(p0, float), np.asarray(p1, float)
    d = p1 - p0
    L = np.hypot(*d)
    nrm = np.array([-d[1], d[0]]) / max(L, 1e-9)
    c = 0.5 * (p0 + p1) + bow * nrm
    t = np.linspace(0, 1, n)[:, None]
    pts = (1 - t) ** 2 * p0 + 2 * (1 - t) * t * c + t ** 2 * p1
    return pts[:, 0], pts[:, 1]


def _wavy(p0, p1, n_per=4, amp=0.022, n=220):
    """Photon/W wiggle: sinusoidal perpendicular offset, nodes at both ends."""
    p0, p1 = np.asarray(p0, float), np.asarray(p1, float)
    d = p1 - p0
    L = np.hypot(*d)
    u = d / max(L, 1e-9)
    nrm = np.array([-u[1], u[0]])
    t = np.linspace(0, 1, n)
    off = amp * np.sin(2 * np.pi * n_per * t) * np.sin(np.pi * t) ** 0.25
    xs = p0[0] + t * d[0] + off * nrm[0]
    ys = p0[1] + t * d[1] + off * nrm[1]
    return xs, ys


def _pack_nucleons(rng):
    """Blue-noise disc packing, seeded with the three cascade nucleons."""
    pts = [_S.copy(), _N1.copy(), _N2.copy()]
    tries = 0
    while len(pts) < 21 and tries < 8000:
        tries += 1
        a = rng.uniform(0, 2 * np.pi)
        r = _RPACK * np.sqrt(rng.uniform())
        p = _NC + r * np.array([np.cos(a), np.sin(a)])
        if all(np.hypot(*(p - q)) > 0.088 for q in pts):
            pts.append(p)
    return pts


# --------------------------------------------------------------- ornament
def draw_nu_nucleus(ax, center, rot_deg=0.0, scale=1.0, alpha=1.0, z=1):
    """Neutrino--nucleus CC event with intranuclear cascade. center = data
    coords of the design center; rot_deg CCW; alpha multiplies every
    element's alpha; artists at zorder z..z+2."""
    rng = np.random.default_rng(20260806)      # fixed: deterministic art

    _c = np.cos(np.radians(rot_deg))
    _s = np.sin(np.radians(rot_deg))
    cx0, cy0 = center

    def T(xs, ys):
        xs = np.asarray(xs, dtype=float)
        ys = np.asarray(ys, dtype=float)
        return (cx0 + scale * (_c * xs - _s * ys),
                cy0 + scale * (_s * xs + _c * ys))

    def T1(p):
        x, y = T(p[0], p[1])
        return float(x), float(y)

    def A(a):
        return a * alpha

    def LW(w):                     # linewidths follow scale, with a floor
        return max(w * scale, 0.28)

    def MS(m):                     # marker/mutation sizes, with a floor
        return max(m * scale, 2.2)

    def P(xs, ys, **kw):
        kw.setdefault("solid_capstyle", "round")
        kw["lw"] = LW(kw["lw"])
        kw["alpha"] = A(kw["alpha"])
        kw["zorder"] = z + kw.pop("dz", 1.0)
        ax.plot(*T(xs, ys), **kw)

    def ARR(p0, p1, color, lw, mut, a, dz=1.0, style="-|>"):
        ax.annotate("", xy=T1(p1), xytext=T1(p0),
                    arrowprops=dict(arrowstyle=style, color=color,
                                    lw=LW(lw), mutation_scale=MS(mut),
                                    shrinkA=0, shrinkB=0, alpha=A(a)),
                    zorder=z + dz)

    def head(xs, ys, color, lw, mut, a, dz=1.0):
        """Arrowhead aligned with the end tangent of a plotted curve."""
        ARR((xs[-4], ys[-4]), (xs[-1], ys[-1]), color, lw, mut, a, dz)

    def dot(p, ms, color, a, dz=1.2):
        ax.plot(*T1(p), "o", ms=MS(ms), color=color, alpha=A(a),
                markeredgewidth=0, zorder=z + dz)

    # ---- nuclear boundary: faint dashed ring + whisper of a halo ------
    ax.add_patch(Circle(T1(_NC), _RB * scale, facecolor=WHITE,
                        edgecolor="none", alpha=A(0.04), zorder=z))
    ax.add_patch(Circle(T1(_NC), _RB * scale, facecolor="none",
                        edgecolor=PALE, lw=LW(0.45), ls=(0, (1.6, 2.4)),
                        alpha=A(0.55), zorder=z + 0.1))

    # ---- nucleons: packed, two shades (protons WHITE / neutrons PALE) --
    pts = _pack_nucleons(rng)
    kinds = rng.uniform(size=len(pts)) < 0.5
    kinds[0] = True                                  # struck one: a proton
    spect = [i for i in range(3, len(pts))]
    fermi = rng.choice(spect, size=3, replace=False)  # Fermi-motion arrows
    for i, p in enumerate(pts):
        if i == 0:
            continue                                 # struck: drawn later
        col = WHITE if kinds[i] else PALE
        ax.add_patch(Circle(T1(p), _RN * scale, facecolor=col,
                            edgecolor="none",
                            alpha=A(0.11 if kinds[i] else 0.07),
                            zorder=z + 0.4))
        ax.add_patch(Circle(T1(p), _RN * scale, facecolor="none",
                            edgecolor=col, lw=LW(0.5),
                            alpha=A(0.80 if kinds[i] else 0.55),
                            zorder=z + 0.5))
    for i in fermi:
        a0 = rng.uniform(0, 2 * np.pi)
        u = 0.024 * np.array([np.cos(a0), np.sin(a0)])
        ARR(pts[i] - u, pts[i] + u, PALE, 0.38, 2.4, 0.35, dz=0.6)

    # ---- incoming neutrino: dashed, with a mid-flight direction tick --
    P([_NU0[0], _V[0]], [_NU0[1], _V[1]], color=PALE, lw=0.55,
      ls=(0, (2.6, 2.0)), alpha=0.85, dz=1.0)
    m0 = _NU0 + 0.54 * (_V - _NU0)
    m1 = _NU0 + 0.60 * (_V - _NU0)
    ARR(m0, m1, PALE, 0.55, 4.0, 0.85, dz=1.0)

    # ---- leptonic vertex ----------------------------------------------
    dot(_V, 1.7, WHITE, 0.95, dz=1.4)

    # ---- outgoing charged lepton: solid, arrowed, slight curve --------
    lx, ly = _bez(_V, _LE, bow=0.055)
    P(lx, ly, color=WHITE, lw=0.7, alpha=0.95, dz=1.2)
    head(lx, ly, WHITE, 0.7, 5.0, 0.95, dz=1.2)

    # ---- W boson: wavy ORANGE, vertex to the struck nucleon's edge ----
    uW = (_S - _V) / np.hypot(*(_S - _V))
    S_edge = _S - (_RN + 0.010) * uW
    wx, wy = _wavy(_V, S_edge, n_per=4, amp=0.021)
    P(wx, wy, color=ORANGE, lw=0.6, alpha=0.95, dz=1.8)

    # ---- momentum transfer q: small offset arrow beside the W ---------
    nW = np.array([uW[1], -uW[0]])                   # lower-left side
    q0 = _V + 0.10 * (S_edge - _V) + 0.075 * nW
    q1 = _V + 0.44 * (S_edge - _V) + 0.075 * nW
    ARR(q0, q1, PALE, 0.5, 3.6, 0.75, dz=1.0)
    qm = 0.5 * (q0 + q1) + 0.052 * nW
    ax.text(*T1(qm), r"$q$", fontsize=max(4.2 * scale, 2.6), color=PALE,
            alpha=A(0.8), ha="center", va="center", zorder=z + 1.0)

    # ---- struck nucleon: THE orange accent, with impact flash ---------
    ax.add_patch(Circle(T1(_S), _RN * scale, facecolor=WHITE,
                        edgecolor="none", alpha=A(0.10), zorder=z + 1.5))
    ax.add_patch(Circle(T1(_S), _RN * scale, facecolor="none",
                        edgecolor=ORANGE, lw=LW(0.75), alpha=A(0.95),
                        zorder=z + 1.6))
    dot(_S, 1.3, ORANGE, 0.9, dz=1.6)
    for af in np.deg2rad([112, 146, 180]):           # burst at the W impact
        r0, r1 = _RN + 0.014, _RN + 0.044
        ca, sa = np.cos(af), np.sin(af)
        P([_S[0] + r0 * ca, _S[0] + r1 * ca],
          [_S[1] + r0 * sa, _S[1] + r1 * sa],
          color=ORANGE, lw=0.5, alpha=0.85, dz=1.6)

    # ---- transparency ghost: straight-through path, very faint --------
    g0 = _S + 0.10 * uW
    g1 = _S + 0.62 * uW
    P([g0[0], g1[0]], [g0[1], g1[1]], color=PALE, lw=0.4,
      ls=(0, (1.4, 2.2)), alpha=0.32, dz=0.7)

    # ---- intranuclear cascade: struck nucleon rescatters twice --------
    u1 = (_N1 - _S) / np.hypot(*(_N1 - _S))
    c1x, c1y = _bez(_S + 0.9 * _RN * u1, _N1, bow=0.014)
    P(c1x, c1y, color=WHITE, lw=0.6, alpha=0.9, dz=1.2)
    c2x, c2y = _bez(_N1, _N2, bow=-0.016)
    P(c2x, c2y, color=WHITE, lw=0.6, alpha=0.9, dz=1.2)
    dot(_N1, 1.3, WHITE, 0.9)
    dot(_N2, 1.3, WHITE, 0.9)

    # leading knocked-out nucleon: N2 -> out through the boundary
    exit1 = np.array([0.55, 0.05])
    e1x, e1y = _bez(_N2, exit1, bow=0.020)
    P(e1x, e1y, color=WHITE, lw=0.6, alpha=0.9, dz=1.2)
    head(e1x, e1y, WHITE, 0.6, 4.6, 0.9, dz=1.2)

    # recoil partner from N1: short stub escaping down-right
    p1 = np.array([0.235, -0.53])
    s1x, s1y = _bez(_N1, p1, bow=-0.020)
    P(s1x, s1y, color=PALE, lw=0.5, alpha=0.8, dz=1.1)
    head(s1x, s1y, PALE, 0.5, 4.0, 0.8, dz=1.1)

    # recoil partner from N2: short stub escaping right-down
    p2 = np.array([0.52, -0.29])
    s2x, s2y = _bez(_N2, p2, bow=0.016)
    P(s2x, s2y, color=PALE, lw=0.5, alpha=0.8, dz=1.1)
    head(s2x, s2y, PALE, 0.5, 4.0, 0.8, dz=1.1)

    # ---- produced pion: fine track with direction chevrons ------------
    pi_end = np.array([-0.27, -0.65])
    px, py = _bez(_N1, pi_end, bow=-0.075, n=140)
    P(px, py, color=PALE, lw=0.42, alpha=0.85, dz=1.1)
    for frac in (0.38, 0.62, 0.86):
        i = int(frac * (len(px) - 1))
        th = np.arctan2(py[i + 1] - py[i - 1], px[i + 1] - px[i - 1])
        for da in (2.55, -2.55):                     # backward chevron legs
            P([px[i] + 0.030 * np.cos(th + da), px[i]],
              [py[i] + 0.030 * np.sin(th + da), py[i]],
              color=PALE, lw=0.42, alpha=0.85, dz=1.1)

    # ---- residual nucleus recoil: small faint arrow, lower left -------
    r0 = _NC + 1.10 * _RB * np.array([np.cos(3.67), np.sin(3.67)])
    r1 = _NC + 1.48 * _RB * np.array([np.cos(3.72), np.sin(3.72)])
    ARR(r0, r1, PALE, 0.45, 3.2, 0.5, dz=0.8)


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

    draw_nu_nucleus(ax, (2.0, 8.0), rot_deg=0.0, scale=1.0)
    draw_nu_nucleus(ax, (1.38, 6.05), rot_deg=0.0, scale=1.0, alpha=0.8)
    draw_nu_nucleus(ax, (6.5, 3.0), rot_deg=-15.0, scale=0.6)

    fig.savefig(HERE / "test.png", dpi=150)
    print("wrote", HERE / "test.png")
