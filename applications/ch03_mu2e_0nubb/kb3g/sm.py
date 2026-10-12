"""Minimal M-scheme shell-model code (KSHELL .snt interactions), pn-factorized, for the Ch. 3 M3 gap checks.

H Psi = Hp Psi + Psi Hn^T + sum_{beta,delta} C_{beta delta} Psi (a+_beta a_delta)_n^T, Psi[p_det, n_det].
Interactions are used as published (no renormalization when orbits are dropped; the truncated Hamiltonian is
the full-space interaction restricted to the kept orbits).
"""
import itertools, math, functools, sys, time
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh, LinearOperator
from sympy.physics.wigner import clebsch_gordan
from sympy import Rational


@functools.lru_cache(maxsize=None)
def cg(j1, m1, j2, m2, J, M):          # all args doubled (2j, 2m) except J,M plain-doubled too
    if m1 + m2 != M or abs(M) > J:
        return 0.0
    return float(clebsch_gordan(Rational(j1, 2), Rational(j2, 2), Rational(J, 2),
                                Rational(m1, 2), Rational(m2, 2), Rational(M, 2)))


def read_snt(path):
    lines = [l.split('!')[0].split() for l in open(path)]
    lines = [l for l in lines if l]
    it = iter(lines)
    npo, nno, _, _ = map(int, next(it))
    orbs = {}
    for _ in range(npo + nno):
        i, n, l, j2, tz = map(int, next(it))
        orbs[i] = (n, l, j2, tz)
    n1, _ = map(int, next(it)[:2])
    spe = {}
    for _ in range(n1):
        i, j, e = next(it)
        assert i == j
        spe[int(i)] = float(e)
    h = next(it)
    n2, method = int(h[0]), int(h[1])
    a0, pw = (float(h[2]), float(h[3])) if method == 1 else (1.0, 0.0)
    tb = {}
    for _ in range(n2):
        i, j, k, l, J, v = next(it)
        key = (int(i), int(j), int(k), int(l), int(J))
        tb[key] = float(v)
    return orbs, spe, tb, (a0, pw)


class Model:
    def __init__(self, snt, keep_orbit_labels, Z, N, A):
        """keep_orbit_labels: set of (n, l, 2j) to keep; Z, N valence particles; A for mass scaling."""
        orbs, spe, tb, (a0, pw) = read_snt(snt)
        self.scale = (A / a0) ** pw
        keep = {i for i, o in orbs.items() if o[:3] in keep_orbit_labels}
        self.orbs = {i: orbs[i] for i in keep}
        # single-particle states per species: list of (orbit, 2m)
        self.sps = {-1: [], 1: []}
        for i in sorted(keep):
            n, l, j2, tz = orbs[i]
            for m2 in range(-j2, j2 + 1, 2):
                self.sps[tz].append((i, m2))
        self.idx = {tz: {s: k for k, s in enumerate(v)} for tz, v in self.sps.items()}
        self.m2 = {tz: np.array([s[1] for s in v]) for tz, v in self.sps.items()}
        self.Z, self.N = Z, N
        self.spe = spe
        # coupled TBME, symmetrized (file lists each unordered (ab),(cd) once)
        V = {}
        for (i, j, k, l, J), v in tb.items():
            if not all(x in keep for x in (i, j, k, l)):
                continue
            for key in ((i, j, k, l, J), (k, l, i, j, J)):
                if key in V:
                    assert abs(V[key] - v) < 1e-9, key
                V[key] = v
        self.V = V
        self._terms()
        self._dets()

    def _terms(self):
        """Uncoupled two-body terms. pp/nn: dict (alpha<beta, gamma<delta) -> c for a+a a+b a_d a_c.
        pn: dict (alpha_p, gamma_p, beta_n, delta_n) -> c for (a+alpha a gamma)_p (a+beta a delta)_n."""
        o = self.orbs
        same = {-1: {}, 1: {}}
        pn = {}
        for (a, b, c, d, J), v in self.V.items():
            v = v * self.scale
            ta, tb_, tc, td = o[a][3], o[b][3], o[c][3], o[d][3]
            ja, jb, jc, jd = o[a][2], o[b][2], o[c][2], o[d][2]
            Nab = 1 / math.sqrt(2) if a == b else 1.0
            Ncd = 1 / math.sqrt(2) if c == d else 1.0
            for M in range(-2 * J, 2 * J + 1, 2):
                for ma in range(-ja, ja + 1, 2):
                    mb = M - ma
                    if abs(mb) > jb:
                        continue
                    c1 = cg(ja, ma, jb, mb, 2 * J, M)
                    if c1 == 0:
                        continue
                    for mc in range(-jc, jc + 1, 2):
                        md = M - mc
                        if abs(md) > jd:
                            continue
                        c2 = cg(jc, mc, jd, md, 2 * J, M)
                        if c2 == 0:
                            continue
                        coef = v * Nab * Ncd * c1 * c2
                        if ta == tb_:
                            t = ta
                            al, be = self.idx[t][(a, ma)], self.idx[t][(b, mb)]
                            ga, de = self.idx[t][(c, mc)], self.idx[t][(d, md)]
                            if al == be or ga == de:
                                continue
                            s = 1
                            if al > be:
                                al, be, s = be, al, -s
                            if ga > de:
                                ga, de, s = de, ga, -s
                            k = (al, be, ga, de)
                            same[t][k] = same[t].get(k, 0.0) + s * coef
                        else:
                            assert ta == -1 and tb_ == 1 and tc == -1 and td == 1
                            k = (self.idx[-1][(a, ma)], self.idx[-1][(c, mc)],
                                 self.idx[1][(b, mb)], self.idx[1][(d, md)])
                            pn[k] = pn.get(k, 0.0) + coef
        self.same = {t: {k: x for k, x in d.items() if abs(x) > 1e-12} for t, d in same.items()}
        self.pn = {k: x for k, x in pn.items() if abs(x) > 1e-12}

    def _dets(self):
        self.dets, self.detM, self.lookup = {}, {}, {}
        for t, n in ((-1, self.Z), (1, self.N)):
            ns = len(self.sps[t])
            ds = np.array([sum(1 << i for i in c) for c in itertools.combinations(range(ns), n)], dtype=np.uint64)
            ds.sort()
            self.dets[t] = ds
            self.detM[t] = np.array([sum(self.m2[t][i] for i in range(ns) if (int(d) >> i) & 1) for d in ds])
        self.Hs = {t: self._species_H(t) for t in (-1, 1)}

    # ---- vectorized fermion strings on det arrays ----
    @staticmethod
    def _below(d, i):
        return np.bitwise_count(d & np.uint64((1 << i) - 1)).astype(np.int64)

    def _apply(self, t, ops, src=None):
        """ops: list of (mode, create?) applied right to left (ops[-1] first). Returns (src_idx, dst_idx, sign)."""
        ds = self.dets[t] if src is None else src
        idx = np.arange(len(ds))
        d = ds.copy()
        sgn = np.ones(len(ds), dtype=np.int64)
        for mode, create in reversed(ops):
            bit = np.uint64(1 << mode)
            occ = (d & bit) != 0
            ok = ~occ if create else occ
            d, idx, sgn = d[ok], idx[ok], sgn[ok]
            sgn = sgn * (1 - 2 * (self._below(d, mode) & 1))
            d = d ^ bit
        dst = np.searchsorted(self.dets[t], d)
        return idx, dst, sgn

    def one_body(self, t, i, j):
        """(a+_i a_j) on species t as sparse matrix."""
        n = len(self.dets[t])
        s, dd, g = self._apply(t, [(i, True), (j, False)])
        return sp.csr_matrix((g.astype(float), (dd, s)), shape=(n, n))

    def _species_H(self, t):
        n = len(self.dets[t])
        rows, cols, vals = [], [], []
        # one-body SPE
        e = np.array([self.spe[o] for (o, m) in self.sps[t]])
        diag = np.zeros(n)
        for i, ei in enumerate(e):
            diag += ei * ((self.dets[t] >> np.uint64(i)) & np.uint64(1)).astype(float)
        rows.append(np.arange(n)); cols.append(np.arange(n)); vals.append(diag)
        if (self.Z if t == -1 else self.N) >= 2:
            for (al, be, ga, de), c in self.same[t].items():
                s, d, g = self._apply(t, [(al, True), (be, True), (de, False), (ga, False)])
                rows.append(d); cols.append(s); vals.append(c * g)
        return sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(n, n))

    def setup_sector(self, M2):
        mp, mn = self.detM[-1], self.detM[1]
        self.mask = (mp[:, None] + mn[None, :]) == M2
        self.M2 = M2
        # pn: group by neutron (beta, delta)
        by_n = {}
        for (al, ga, be, de), c in self.pn.items():
            by_n.setdefault((be, de), []).append((al, ga, c))
        self.pn_ops = []
        Np = len(self.dets[-1])
        if self.Z == 0 or self.N == 0:
            return
        cache = {}
        for (be, de), lst in by_n.items():
            Cp = sp.csr_matrix((Np, Np))
            for al, ga, c in lst:
                if (al, ga) not in cache:
                    cache[(al, ga)] = self.one_body(-1, al, ga)
                Cp = Cp + c * cache[(al, ga)]
            s, d, g = self._apply(1, [(be, True), (de, False)])
            self.pn_ops.append((Cp.tocsr(), s, d, g.astype(float)))

    def matvec_full(self, Psi):
        out = self.Hs[-1] @ Psi + (self.Hs[1] @ Psi.T).T
        for Cp, s, d, g in self.pn_ops:
            Phi = Psi[:, s] * g
            out[:, d] += Cp @ Phi
        return out

    def lowest(self, M2, k=8, tol=1e-9):
        self.setup_sector(M2)
        dim = int(self.mask.sum())
        shape = self.mask.shape

        def mv(x):
            Psi = np.zeros(shape)
            Psi[self.mask] = x
            return self.matvec_full(Psi)[self.mask]
        op = LinearOperator((dim, dim), matvec=mv, dtype=float)
        k = min(k, dim - 1) if dim > 1 else 1
        if dim <= 1500:
            H = np.column_stack([mv(c) for c in np.eye(dim)])
            w, v = np.linalg.eigh(H)
            w, v = w[:k], v[:, :k]
        else:
            w, v = eigsh(op, k=k, which='SA', tol=tol)
            o = np.argsort(w); w, v = w[o], v[:, o]
        return dim, w, v

    def J2(self, vec):
        """<J^2> = |J+ psi|^2 + M^2 + M (M = M2/2)."""
        Psi = np.zeros(self.mask.shape)
        Psi[self.mask] = vec
        out = np.zeros_like(Psi)
        for t in (-1, 1):
            for (o, m2), i in self.idx[t].items():
                if (o, m2 + 2) in self.idx[t]:
                    j2 = self.orbs[o][2]
                    f = math.sqrt((j2 - m2) / 2 * ((j2 + m2) / 2 + 1))
                    j = self.idx[t][(o, m2 + 2)]
                    s, d, g = self._apply(t, [(j, True), (i, False)])
                    if t == -1:
                        out[d, :] += f * g[:, None] * Psi[s, :]
                    else:
                        out[:, d] += f * g[None, :] * Psi[:, s]
        M = self.M2 / 2
        return float(np.sum(out ** 2)) + M * M + M


def spectrum(model, M2, k=8, label=""):
    t0 = time.time()
    dim, w, v = model.lowest(M2, k)
    Js = []
    for i in range(len(w)):
        j2 = model.J2(v[:, i])
        J = (-1 + math.sqrt(1 + 4 * j2)) / 2
        Js.append(J)
    print(f"{label}: M-dim {dim}, {time.time()-t0:.1f}s")
    for e, J in zip(w, Js):
        print(f"   E = {e:12.5f}  Ex = {e - w[0]:8.4f}  J = {J:6.3f}")
    sys.stdout.flush()
    return dim, w, Js
