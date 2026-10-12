"""Schematic lambda/|M| for the 48Ca -> 48Ti 0nubb-like two-body operator, Pauli-LCU 1-norm.
Shell-model toy: HO radial functions (b from hw = 41 A^-1/3), delta-force (SDI-like) Hamiltonian with
single-particle energies, pn formalism. Operators: short-range contact delta(r12) and long-range 1/r12
(Fermi-like proxy of the light-Majorana neutrino potential H(r) ~ R/r). Read-only scratch script."""
import itertools, math, sys
import numpy as np
from scipy.special import sph_harm_y, eval_genlaguerre
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import eigsh
from sympy import Rational
from sympy.physics.wigner import clebsch_gordan

A_MASS = 48
KINDS = ("gt_coulomb",)
HW = 41.0 * A_MASS ** (-1 / 3)
B = 197.327 / math.sqrt(938.92 * HW)

ORB = {"f7": (0, 3, 7), "p3": (1, 1, 3), "f5": (0, 3, 5), "p1": (1, 1, 1)}
SPE = {"f7": 0.0, "p3": 2.0, "p1": 4.0, "f5": 6.5}   # 41Ca-like single-particle energies (MeV)

# quadrature grids
xg, wg = np.polynomial.legendre.leggauss(24)
NPHI = 48
TH = np.arccos(xg)
PH = np.arange(NPHI) * 2 * np.pi / NPHI
THg, PHg = np.meshgrid(TH, PH, indexing="ij")
WANG = np.outer(wg, np.full(NPHI, 2 * np.pi / NPHI))
rg = np.linspace(1e-4, 8 * B, 1600)
dr = rg[1] - rg[0]


def radial(n, l, r):
    x = r / B
    norm = math.sqrt(2 * math.factorial(n) / (B ** 3 * math.gamma(n + l + 1.5)))
    return norm * x ** l * np.exp(-x * x / 2) * eval_genlaguerre(n, l + 0.5, x * x)


def cg(j1, m1, j2, m2, J, M):
    return float(clebsch_gordan(Rational(j1, 2), Rational(j2, 2), Rational(J, 2), Rational(m1, 2), Rational(m2, 2), Rational(M, 2)))


def spinor(orb, m2):
    n, l, tj = ORB[orb]
    out = np.zeros((2,) + THg.shape, complex)   # index 0: s=+1/2, 1: s=-1/2
    for si, s2 in enumerate((1, -1)):
        ml2 = m2 - s2
        if abs(ml2) > 2 * l:
            continue
        c = cg(2 * l, ml2, 1, s2, tj, m2)
        if c != 0:
            out[si] = c * sph_harm_y(l, ml2 // 2, THg, PHg)
    return out


def states(orbs):
    return [(o, m2) for o in orbs for m2 in range(-ORB[o][2], ORB[o][2] + 1, 2)]


class TBME:
    def __init__(self, orbs, kind):
        self.st = states(orbs)
        self.kind = kind
        self.sp = {s: spinor(*s) for s in self.st}
        self.R = {o: radial(ORB[o][0], ORB[o][1], rg) for o in orbs}
        kmax = 2 * max(ORB[o][1] for o in orbs)
        self.Ykq = {(k, q): sph_harm_y(k, q, THg, PHg) for k in range(kmax + 1) for q in range(-k, k + 1)}
        self.cache = {}
        self.rcache = {}

    def ang(self, a, c, k, q):
        rho = (np.conj(self.sp[a]) * self.sp[c]).sum(0)
        return (rho * self.Ykq[(k, q)] * WANG).sum()

    def rad(self, oa, ob, oc, od, k):
        key = (oa, ob, oc, od, k)
        if key in self.rcache:
            return self.rcache[key]
        f1 = self.R[oa] * self.R[oc] * rg ** 2
        f2 = self.R[ob] * self.R[od] * rg ** 2
        if self.kind in ("delta", "gt_delta"):
            v = (self.R[oa] * self.R[ob] * self.R[oc] * self.R[od] * rg ** 2).sum() * dr
        else:  # 1/r12 multipole: 4pi/(2k+1) r<^k / r>^(k+1)
            r1 = rg[:, None]; r2 = rg[None, :]
            g = np.where(r1 < r2, r1 ** k / r2 ** (k + 1), r2 ** k / r1 ** (k + 1)) * 4 * np.pi / (2 * k + 1)
            v = f1 @ g @ f2 * dr * dr
        self.rcache[key] = v
        return v

    def dens(self, a, c, mu):
        """sum_{s s'} Phi*_a(s) S_mu[s,s'] Phi_c(s'); mu=None identity, else spherical sigma_mu."""
        A_, C_ = np.conj(self.sp[a]), self.sp[c]
        if mu is None:
            return (A_ * C_).sum(0)
        r2 = math.sqrt(2)
        if mu == 0:
            return A_[0] * C_[0] - A_[1] * C_[1]
        if mu == 1:
            return -r2 * A_[0] * C_[1]
        return r2 * A_[1] * C_[0]

    def me(self, a, b, c, d):
        if not self.kind.startswith("gt"):
            return self.me0(a, b, c, d)
        key = (a, b, c, d)
        if key in self.cache:
            return self.cache[key]
        if a[1] + b[1] != c[1] + d[1]:
            return 0.0
        la, lb, lc, ld = (ORB[x[0]][1] for x in (a, b, c, d))
        if (la + lb + lc + ld) % 2:
            return 0.0
        tot = 0.0
        for mu in (-1, 0, 1):
            t1 = self.dens(a, c, mu); t2 = self.dens(b, d, -mu)
            ph = (-1) ** mu
            if self.kind == "gt_delta":
                tot += ph * self.rad(a[0], b[0], c[0], d[0], 0) * (t1 * t2 * WANG).sum()
                continue
            for k in range(min(la + lc, lb + ld) + 1):
                if (la + lc + k) % 2 or (lb + ld + k) % 2:
                    continue
                rk = self.rad(a[0], b[0], c[0], d[0], k)
                for q in range(-k, k + 1):
                    i1 = (t1 * np.conj(self.Ykq[(k, q)]) * WANG).sum()
                    i2 = (t2 * self.Ykq[(k, q)] * WANG).sum()
                    tot += ph * rk * i1 * i2
        tot = tot.real
        self.cache[key] = tot
        return tot

    def me0(self, a, b, c, d):
        """<ab|f|cd> direct, distinguishable-particle (no species check)."""
        key = (a, b, c, d)
        if key in self.cache:
            return self.cache[key]
        if a[1] + b[1] != c[1] + d[1]:
            return 0.0
        la, lb, lc, ld = (ORB[x[0]][1] for x in (a, b, c, d))
        if (la + lb + lc + ld) % 2:
            return 0.0
        kmax = min(la + lc, lb + ld)
        tot = 0.0
        for k in range(kmax + 1):
            if (la + lc + k) % 2 or (lb + ld + k) % 2:
                continue
            rk = self.rad(a[0], b[0], c[0], d[0], k)
            for q in range(-k, k + 1):
                # f = sum_k g_k sum_q Y*_kq(1) Y_kq(2): <a|Y*_kq|c> <b|Y_kq|d>
                t1 = (np.conj(self.sp[a]) * self.sp[c]).sum(0)
                t2 = (np.conj(self.sp[b]) * self.sp[d]).sum(0)
                i1 = (t1 * np.conj(self.Ykq[(k, q)]) * WANG).sum()
                i2 = (t2 * self.Ykq[(k, q)] * WANG).sum()
                tot += rk * (i1 * i2).real if self.kind != "delta" else 0.0
        if self.kind == "delta":
            # delta: g_k = delta(r1-r2)/r^2, sum_kq Y*Y = delta(Omega1-Omega2): direct angular overlap
            t1 = (np.conj(self.sp[a]) * self.sp[c]).sum(0)
            t2 = (np.conj(self.sp[b]) * self.sp[d]).sum(0)
            tot = self.rad(a[0], b[0], c[0], d[0], 0) * (t1 * t2 * WANG).sum().real
        self.cache[key] = tot
        return tot


def popcount(x):
    return bin(x).count("1")


def apply_2b(state, a, b, c, d):
    """a_a^+ a_b^+ a_d a_c |state>, modes as bit positions; returns (sign, newstate) or None."""
    s = state
    sign = 1
    for op, mode in ((-1, c), (-1, d), (1, b), (1, a)):  # rightmost first: a_c, then a_d, then a_b^+, then a_a^+
        occ = (s >> mode) & 1
        if op == -1 and not occ or op == 1 and occ:
            return None
        sign *= -1 if popcount(s & ((1 << mode) - 1)) % 2 else 1
        s ^= 1 << mode
    return sign, s


def run(orbs, kind_H="delta", V0=None, Z=2, N=6, label=""):
    st = states(orbs)
    ns = len(st)
    # modes: protons 0..ns-1, neutrons ns..2ns-1
    mm = [s[1] for s in st] * 2
    spe = [SPE[s[0]] for s in st] * 2
    Hme = TBME(orbs, "delta")

    def basis(zp, nn):
        out = []
        for pc in itertools.combinations(range(ns), zp):
            for nc in itertools.combinations(range(ns, 2 * ns), nn):
                if sum(mm[i] for i in pc + nc) == 0:
                    out.append(sum(1 << i for i in pc + nc))
        return out

    # two-body list: (a<b, c<d) with V_AS, species-aware, delta force strength -V0
    def vas(a, b, c, d):
        sa, sb, sc, sd = (x >= ns for x in (a, b, c, d))
        A_, Bb, C_, D_ = (st[x % ns] for x in (a, b, c, d))
        v = 0.0
        if sa == sc and sb == sd:
            v += Hme.me(A_, Bb, C_, D_)
        if sa == sd and sb == sc:
            v -= Hme.me(A_, Bb, D_, C_)
        return -V0 * v

    pairs = list(itertools.combinations(range(2 * ns), 2))
    byM = {}
    for p in pairs:
        byM.setdefault((mm[p[0]] + mm[p[1]], (p[0] >= ns) + (p[1] >= ns)), []).append(p)
    V = {}
    for key, pl in byM.items():
        for (a, b) in pl:
            for (c, d) in pl:
                v = vas(a, b, c, d)
                if abs(v) > 1e-10:
                    V[(a, b, c, d)] = v

    def ham(bas):
        idx = {s: i for i, s in enumerate(bas)}
        rows, cols, vals = [], [], []
        for i, s in enumerate(bas):
            e = sum(spe[k] for k in range(2 * ns) if (s >> k) & 1)
            rows.append(i); cols.append(i); vals.append(e)
        # group V by (c,d)
        byCD = {}
        for (a, b, c, d), v in V.items():
            byCD.setdefault((c, d), []).append((a, b, v))
        for i, s in enumerate(bas):
            occ = [k for k in range(2 * ns) if (s >> k) & 1]
            for c, d in itertools.combinations(occ, 2):
                for a, b, v in byCD.get((c, d), ()):
                    r = apply_2b(s, a, b, c, d)
                    if r is None:
                        continue
                    j = idx.get(r[1])
                    if j is not None:
                        rows.append(j); cols.append(i); vals.append(v * r[0])
        n = len(bas)
        return coo_matrix((vals, (rows, cols)), shape=(n, n)).tocsr(), idx

    def gs(bas, nev=2):
        H, idx = ham(bas)
        if H.shape[0] <= 400:
            w, v = np.linalg.eigh(H.toarray())
        else:
            w, v = eigsh(H, k=min(nev + 4, H.shape[0] - 1), which="SA")
            o = np.argsort(w); w, v = w[o], v[:, o]
        return w, v, idx

    bi = basis(0, Z + N)          # 48Ca: 8 neutrons
    bf = basis(Z, N)              # 48Ti: 2p 6n
    wi, vi, idxi = gs(bi)
    wf, vf, idxf = gs(bf)
    closed = sum(1 << (ns + k) for k in range(ns) if st[k][0] == "f7")
    out = dict(label=label, dim_i=len(bi), dim_f=len(bf), p_ref=float(vi[idxi[closed], 0] ** 2))
    for kind in KINDS:
        Om = TBME(orbs, kind)
        terms = {}
        for a, b in itertools.combinations(range(ns), 2):
            for c, d in itertools.combinations(range(ns), 2):
                if mm[a] + mm[b] != mm[c] + mm[d]:
                    continue
                v = Om.me(st[a], st[b], st[c], st[d]) - Om.me(st[a], st[b], st[d], st[c])
                if abs(v) > 1e-10:
                    terms[(a, b, c + ns, d + ns)] = v
        lam = sum(abs(v) for v in terms.values())
        # O|i>
        psi = vi[:, 0]
        res = {}
        for s, i in idxi.items():
            amp = psi[i]
            if abs(amp) < 1e-14:
                continue
            for (a, b, c, d), v in terms.items():
                r = apply_2b(s, a, b, c, d)
                if r is None:
                    continue
                res[r[1]] = res.get(r[1], 0.0) + amp * v * r[0]
        vec = np.zeros(len(bf))
        for s, x in res.items():
            vec[idxf[s]] += x
        M = float(vf[:, 0] @ vec)
        S = float(vec @ vec)
        out[kind] = dict(n_terms=len(terms), lam=lam, M=M, S=S, lam_over_M=lam / abs(M) if M else float("inf"),
                         p_gs=(M / lam) ** 2, p_sum=S / lam ** 2, gs_fraction=M * M / S)
    out["Ei"] = wi[:3].tolist(); out["Ef"] = wf[:3].tolist()
    return out


def fit_v0(orbs):
    """V0 so that the two-neutron f7/2^2 J=0 - J=2 splitting is ~1.5 MeV (42Ca: 0+ -> 2+ 1.525 MeV)."""
    H = TBME(["f7"], "delta")
    st = states(["f7"])
    # pairing energy of J=0 pair vs J=2 via m-scheme: use M=0 two-neutron states
    pairs = [(i, j) for i, j in itertools.combinations(range(len(st)), 2) if st[i][1] + st[j][1] == 0]
    Hm = np.zeros((len(pairs), len(pairs)))
    for x, (a, b) in enumerate(pairs):
        for y, (c, d) in enumerate(pairs):
            Hm[x, y] = -(H.me(st[a], st[b], st[c], st[d]) - H.me(st[a], st[b], st[d], st[c]))
    w = np.linalg.eigvalsh(Hm)
    w = np.unique(np.round(w, 9))
    return 1.525 / (w[1] - w[0]), w


if __name__ == "__main__":
    v0, w = fit_v0(["f7"])
    print("b =", B, "fm; V0 =", v0, "MeV fm^3; f7^2 levels per unit V0:", w)
    for orbs in (["f7"], ["f7", "p3"]) + ((["f7", "p3", "p1"],) if "big" in sys.argv else ()):
        o = run(orbs, V0=v0, label="+".join(orbs))
        print(o["label"], "p_ref(closed f7)", round(o["p_ref"], 4), "dims", o["dim_i"], o["dim_f"], "Ei", np.round(o["Ei"], 3), "Ef", np.round(o["Ef"], 3))
        for kind in KINDS:
            d = o[kind]
            print(f"   {kind:8s} terms {d['n_terms']:5d} lam {d['lam']:.4g} M {d['M']:.4g} lam/|M| {d['lam_over_M']:.3g}"
                  f"  p_gs {d['p_gs']:.3g}  p_sum {d['p_sum']:.3g}  gs share of sum {d['gs_fraction']:.3g}")
