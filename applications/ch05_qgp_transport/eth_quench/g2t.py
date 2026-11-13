"""Binary tetrahedral (2T) pure gauge in 2+1D on an open Lx x Ly plaquette patch, maximal-tree gauge fixed
(all horizontal links + the x=0 vertical column fixed to e; root (0,0)); physical states = singlets of the residual
global conjugation, handled exactly in the orbit basis.  H = (g^2/2)(3/16) sum_l sum_{h in N}(1 - L_h)
+ (1/g^2) sum_p (2 - tr_F U_p), N = the 8 elements with tr_F = 1 (KS-matched on the fundamental: E^2 = 3/4)."""
import numpy as np, scipy.sparse as sp, sys, itertools, time, os
from qa import analyze, dump, canon, T_of_energy

Lx, Ly, g2 = int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3])
# ---- group
el = []
for s in itertools.product([1, -1], repeat=1):
    pass
for i in range(4):
    for sgn in (1, -1):
        q = [0, 0, 0, 0]; q[i] = sgn; el.append(tuple(q))
for sg in itertools.product([0.5, -0.5], repeat=4):
    el.append(tuple(sg))
el = [np.array(e, float) for e in el]
G = len(el); assert G == 24


def qmul(a, b):
    return np.array([a[0]*b[0]-a[1]*b[1]-a[2]*b[2]-a[3]*b[3], a[0]*b[1]+a[1]*b[0]+a[2]*b[3]-a[3]*b[2],
                     a[0]*b[2]-a[1]*b[3]+a[2]*b[0]+a[3]*b[1], a[0]*b[3]+a[1]*b[2]-a[2]*b[1]+a[3]*b[0]])


def find(q):
    for i, e in enumerate(el):
        if np.allclose(e, q):
            return i
    raise ValueError
mul = np.array([[find(qmul(a, b)) for b in el] for a in el])
inv = np.array([find(np.array([e[0], -e[1], -e[2], -e[3]])) for e in el])
q0 = np.array([e[0] for e in el])
ident = find(np.array([1, 0, 0, 0.]))
Nset = [i for i in range(G) if abs(q0[i] - 0.5) < 1e-9]; assert len(Nset) == 8
# ---- lattice
sid = lambda x, y: x + (Lx + 1) * y
links = []  # (a, b, kind, x, y)
for y in range(Ly + 1):
    for x in range(Lx):
        links.append((sid(x, y), sid(x + 1, y), 'h', x, y))
for y in range(Ly):
    for x in range(Lx + 1):
        links.append((sid(x, y), sid(x, y + 1), 'v', x, y))
lid = {(k, x, y): i for i, (_, _, k, x, y) in enumerate(links)}
tree = [i for i, l in enumerate(links) if l[2] == 'h' or (l[2] == 'v' and l[3] == 0)]
nontree = [i for i in range(len(links)) if i not in tree]
n = len(nontree); D = G ** n
pos = {l: j for j, l in enumerate(nontree)}
plaqs = [(x, y) for y in range(Ly) for x in range(Lx)]
cfg = np.array(list(np.ndindex(*(G,) * n)), dtype=np.int64) if n else np.zeros((1, 0), int)
wts = G ** np.arange(n - 1, -1, -1)
index = lambda c: c @ wts


def U(c, l):
    return c[:, pos[l]] if l in pos else np.full(len(c), ident)


def holonomy(c, x, y):
    a = U(c, lid[('h', x, y)]); b = U(c, lid[('v', x + 1, y)])
    cc = inv[U(c, lid[('h', x, y + 1)])]; d = inv[U(c, lid[('v', x, y)])]
    return mul[mul[mul[a, b], cc], d]
hol = {p: holonomy(cfg, *p) for p in plaqs}
HB = {p: sp.diags((2 - 2 * q0[hol[p]]) / g2) for p in plaqs}
# subtree sides
nsite = (Lx + 1) * (Ly + 1)
adj = {s: [] for s in range(nsite)}
for i in tree:
    a, b = links[i][:2]; adj[a].append((b, i)); adj[b].append((a, i))


def side(cut):
    seen = {0}; st = [0]
    while st:
        u = st.pop()
        for v, i in adj[u]:
            if i != cut and v not in seen:
                seen.add(v); st.append(v)
    return set(range(nsite)) - seen   # component without the root


def act(l, h):
    c = cfg.copy()
    if l in pos:
        c[:, pos[l]] = mul[h, cfg[:, pos[l]]]
        return index(c)
    a, b = links[l][:2]; S = side(l)
    w = h if b in S else inv[h]; wi = inv[w]
    for j in nontree:
        x, y = links[j][:2]; g = cfg[:, pos[j]]
        if x in S and y in S: c[:, pos[j]] = mul[mul[w, g], wi]
        elif y in S: c[:, pos[j]] = mul[g, wi]
        elif x in S: c[:, pos[j]] = mul[w, g]
    return index(c)
coef = 0.5 * g2 * 3 / 16
HE = {}
for l in range(len(links)):
    M = sp.csr_matrix((D, D))
    for hh in Nset:
        tgt = act(l, hh)
        M = M + sp.csr_matrix((np.ones(D), (tgt, np.arange(D))), shape=(D, D))
    HE[l] = coef * (8 * sp.identity(D, format='csr') - M)
# orbit (singlet) basis
lab = np.arange(D)
for w in range(G):
    lab = np.minimum(lab, index(mul[mul[w, cfg], inv[w]]))
uniq, orb = np.unique(lab, return_inverse=True)
size = np.bincount(orb)
P = sp.csr_matrix((1 / np.sqrt(size[orb]), (np.arange(D), orb)), shape=(D, len(uniq)))
proj = lambda M: (P.T @ M @ P).tocsr()
print(f"2T {Lx}x{Ly}: D_fixed={D}, singlet={len(uniq)}", flush=True)
H = proj(sum(HE.values()) + sum(HB.values()))
# local energy density per plaquette
share = {l: 0 for l in range(len(links))}
links_of = {}
for (x, y) in plaqs:
    ls = [lid[('h', x, y)], lid[('v', x + 1, y)], lid[('h', x, y + 1)], lid[('v', x, y)]]
    links_of[(x, y)] = ls
    for l in ls: share[l] += 1
hden = {p: HB[p] + sum(HE[l] / share[l] for l in links_of[p]) for p in plaqs}
cosm = {p: np.cos(np.pi * (p[0] + 0.5) / Lx) for p in plaqs}   # lowest open-boundary mode along x
ops_local = {"plaq_trF/2": proj(sp.diags(sum(q0[hol[p]] for p in plaqs) / len(plaqs))),
             "E2_link": proj(sum(HE.values()) / (coef * len(links)))}
ops_corr = {"plaq_kmin": proj(sp.diags(sum(cosm[p] * q0[hol[p]] for p in plaqs))),
            "eps_kmin": proj(sum(cosm[p] * hden[p] for p in plaqs))}
# initial states: vacuum (E=0) and fundamental Wilson loops on plaquette subsets (class functions of U_p)
vac = P.T @ np.ones(D)


def loops(S, rep):
    f = np.ones(D)
    for p in S:
        f = f * (2 * q0[hol[p]] if rep == 'F' else 4 * q0[hol[p]] ** 2 - 1)
    return P.T @ f
NS = int(sys.argv[4]) if len(sys.argv) > 4 else 12
rng = np.random.default_rng(7)
Hcomp = sum(HE.values()) + sum(HB.values())


def ensembles(E):
    Ev = vac @ (H @ vac) / (vac @ vac)
    Tstar = T_of_energy(E, Ev)
    ens = {"vac": [(1.0, vac)]}
    Ts = [Tstar, 0.749] if not np.isnan(Tstar) else [0.749]
    subsets = [S for r in range(len(plaqs) + 1) for S in itertools.combinations(plaqs, r)]
    vecs = {S: loops(S, 'F') for S in subsets}
    En = {S: v @ (H @ v) / (v @ v) for S, v in vecs.items()}
    TL = [float(x) for x in os.environ.get("TLIST", "").split(",") if x] or sorted(set([0.749, 1.5 * Tstar, 2.5 * Tstar]))
    for T in TL:
        Et = canon(E, T) @ E
        if Et <= Ev:
            print(f"  T={T:.3f}: E(T)={Et:.3f} below vacuum-state energy {Ev:.3f}: unreachable by E-basis states")
            continue
        lo, hi = 0.0, 1.0
        for _ in range(60):
            q = 0.5 * (lo + hi)
            m = sum(q ** len(S) * (1 - q) ** (len(plaqs) - len(S)) * En[S] for S in subsets)
            if m < Et: lo = q
            else: hi = q
        if abs(m - Et) > 1e-3:
            print(f"  T={T:.3f}: loop ensemble cannot reach E(T)={Et:.3f} (max {m:.3f})"); continue
        pw = np.array([q ** len(S) * (1 - q) ** (len(plaqs) - len(S)) for S in subsets])
        NSAMP = int(os.environ.get("NSAMP", "0"))
        if NSAMP:
            pick = rng.choice(len(subsets), size=NSAMP, p=pw / pw.sum())
            ens[f"loopsF_T{T:.3f}"] = [(1.0, vecs[subsets[i]]) for i in pick]
        else:
            ens[f"loopsF_T{T:.3f}"] = [(w, vecs[S]) for w, S in zip(pw, subsets) if w > 1e-6]
        Ts.append(T)
    return ens, Ts


if "--check" in sys.argv:
    Hc = Hcomp
    err = np.abs(P @ (P.T @ (Hc @ P)) - Hc @ P).max()
    v1 = loops([plaqs[0]], 'F'); v0 = vac
    print("singlet closure err", err, " E_E(fund loop) - vac:", (v1 @ (proj(sum(HE.values())) @ v1)) / (v1 @ v1), "expected", 1.5 * g2)
    sys.exit()
tl = np.linspace(0, 30, 241)
tc = np.linspace(0, 20, 241)
tw = [float(x) for x in os.environ.get('TWS', '0,2,8,16').split(',')]
res = analyze(f"2T_{Lx}x{Ly}_g2{g2}", H, ops_local, ops_corr, ensembles, [], tl, tc, tw,
              extra={"model": "2T open patch, gauge fixed", "Lx": Lx, "Ly": Ly, "g2": g2, "D_singlet": len(uniq)})
dump(res, f"out_2t_{Lx}x{Ly}_g{g2}.json")
