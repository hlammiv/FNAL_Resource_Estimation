"""1+1D Z2 gauge theory + one staggered fermion, open chain of N sites, Gauss law solved (left boundary flux +1),
neutral sector (N/2 fermions). H = -w sum (c^dag_x c_{x+1} + h.c.) + m sum (-1)^x n_x + hE sum_l (1 - sigma^x_l)/2,
sigma^x_l = prod_{y<=l} (-1)^{n_y - s_y}, s_y = y mod 2 (odd sites filled in the strong-coupling vacuum)."""
import numpy as np, scipy.sparse as sp, sys, itertools
from qa import analyze, dump, canon, T_of_energy

N, w, m, hE = int(sys.argv[1]), 1.0, float(sys.argv[2]), float(sys.argv[3])
NS = int(sys.argv[4]) if len(sys.argv) > 4 else 12
states = np.array([sum(1 << i for i in c) for c in itertools.combinations(range(N), N // 2)], dtype=np.int64)
D = len(states)
idx = {s: i for i, s in enumerate(states)}
occ = ((states[:, None] >> np.arange(N)[None, :]) & 1).astype(float)
s_bg = (np.arange(N) % 2).astype(float)
par = np.cumprod(1 - 2 * ((occ - s_bg) % 2), axis=1)[:, :N - 1]          # sigma^x on links 0..N-2
Emass = m * occ @ ((-1.0) ** np.arange(N))
Eflux = hE * (1 - par) / 2
Ecl = Emass + Eflux.sum(1)
rows, cols, vals = [np.arange(D)], [np.arange(D)], [Ecl]
for x in range(N - 1):
    a = (states >> x) & 1; b = (states >> (x + 1)) & 1
    sel = np.where(a != b)[0]
    tgt = np.array([idx[s ^ (3 << x)] for s in states[sel]])
    rows.append(tgt); cols.append(sel); vals.append(-w * np.ones(len(sel)))
H = sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(D, D))
H = 0.5 * (H + H.T)  # each hop counted from both ends above: symmetrize (each (i,j) pair appears twice with -w)
H = H.tolil(); H.setdiag(Ecl); H = H.tocsr()
# local energy density: site mass + half of adjacent link fluxes + half of adjacent hops
cosm = np.cos(np.pi * (np.arange(N) + 0.5) / N)
lc = np.cos(np.pi * (np.arange(N - 1) + 1.0) / N)
hop = {}
for x in range(N - 1):
    a = (states >> x) & 1; b = (states >> (x + 1)) & 1
    sel = np.where(a != b)[0]
    tgt = np.array([idx[s ^ (3 << x)] for s in states[sel]])
    hop[x] = sp.csr_matrix((-w * np.ones(len(sel)), (tgt, sel)), shape=(D, D))
ops_local = {"flux_density": sp.diags((1 - par).mean(1) / 2), "pair_density": sp.diags(np.abs(occ - s_bg).mean(1))}
ops_corr = {"charge_kmin": sp.diags((occ - s_bg) @ cosm),
            "eps_kmin": sp.diags((occ * ((-1.0) ** np.arange(N)) * m) @ cosm + Eflux @ lc) + sum(lc[x] * hop[x] for x in hop)}


def unit(i):
    v = np.zeros(D); v[i] = 1.0; return v
rng = np.random.default_rng(3)


def ensembles(E):
    vac = idx[int(sum(1 << i for i in range(1, N, 2)))]
    Tstar = T_of_energy(E, Ecl[vac])
    ens = {"vac": [(1.0, unit(vac))]}
    Ts = [Tstar]
    for fac in (1.5, 2.5):
        T = Tstar * fac; Et = canon(E, T) @ E
        lo, hi = -3, 3
        for _ in range(100):
            mm = 0.5 * (lo + hi); p = np.exp(-(Ecl - Ecl.min()) * 10 ** mm); p /= p.sum()
            if p @ Ecl > Et: lo = mm
            else: hi = mm
        pick = rng.choice(D, size=NS, p=p)
        ens[f"fockE_{fac}Tstar"] = [(1.0, unit(i)) for i in pick]
        Ts.append(T)
    return ens, Ts


tl = np.linspace(0, 30, 241); tc = np.linspace(0, 20, 241); tw = [0.0, 2.0, 8.0, 16.0]
res = analyze(f"Z2F_N{N}_m{m}_hE{hE}", H, ops_local, ops_corr, ensembles, [], tl, tc, tw,
              extra={"model": "1+1D Z2 + staggered fermion, OBC", "N": N, "m": m, "hE": hE, "w": w})
dump(res, f"out_z2f_N{N}_m{m}_h{hE}.json")
