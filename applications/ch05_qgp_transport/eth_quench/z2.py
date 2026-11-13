"""Z2 lattice gauge theory in 2+1D on an Lx x Ly torus, zero-winding sector, via its exact dual:
H = -J sum_p B_p - h sum_l sigma^x_l  <->  H_dual = -J sum_p X_p - h sum_<pq> Z_p Z_q on the Z2-even sector
(prod_p B_p = 1). Electric-basis (sigma^x) Gauss-law states <-> dual Z-configurations (pairs s, ~s)."""
import numpy as np, scipy.sparse as sp, sys, json, time, os
from qa import analyze, dump, canon, T_of_energy

Lx, Ly, h = int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3])
J = 1.0
N = Lx * Ly
site = lambda x, y: (x % Lx) + Lx * (y % Ly)
bonds = [(site(x, y), site(x + 1, y), 'x') for y in range(Ly) for x in range(Lx)] + \
        [(site(x, y), site(x, y + 1), 'y') for y in range(Ly) for x in range(Lx)]
full = np.arange(2 ** N, dtype=np.int64)
reps = full[(full & 1) == 0]
D = len(reps)
mask = 2 ** N - 1
idx = -np.ones(2 ** N, dtype=np.int64); idx[reps] = np.arange(D); idx[reps ^ mask] = np.arange(D)
spin = lambda s, p: 1 - 2 * ((s >> p) & 1)
zz = {b: spin(reps, b[0]) * spin(reps, b[1]) for b in bonds}
Ecl = -h * sum(zz.values()).astype(float)
rows, cols, vals = [np.arange(D)], [np.arange(D)], [Ecl]
Xop = {}
for p in range(N):
    tgt = idx[reps ^ (1 << p)]
    Xop[p] = sp.csr_matrix((np.ones(D), (tgt, np.arange(D))), shape=(D, D))
    rows.append(tgt); cols.append(np.arange(D)); vals.append(-J * np.ones(D))
H = sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(D, D))
xp = np.array([p % Lx for p in range(N)]); kx = 2 * np.pi / Lx
ops_local = {"B": sum(Xop.values()) / N, "Elink": sp.diags(sum(zz.values()) / len(bonds))}
ops_corr = {"B_kmin": sum(np.cos(kx * xp[p]) * Xop[p] for p in range(N)),
            "eps_kmin": sum(np.cos(kx * xp[p]) * (-J) * Xop[p] for p in range(N)) +
                        sp.diags(sum(0.5 * (np.cos(kx * xp[b[0]]) + np.cos(kx * xp[b[1]])) * (-h) * zz[b] for b in bonds)),
            "Eaniso_kmin": sp.diags(sum(np.cos(kx * xp[b[0]]) * (1 if b[2] == 'x' else -1) * zz[b] for b in bonds))}
rng = np.random.default_rng(12345)
NS = int(sys.argv[4]) if len(sys.argv) > 4 else 12


def unit(i):
    v = np.zeros(D); v[i] = 1.0; return v


def ensembles(E):
    vac = idx[0]
    Tstar = T_of_energy(E, Ecl[vac])
    ens = {"vac": [(1.0, unit(vac))]}
    Ts = [Tstar] if not np.isnan(Tstar) else []
    for fac in (1.5, 2.5):
        T = Tstar * fac
        Et = canon(E, T) @ E
        if Et >= 0:
            continue
        # classical Boltzmann over electric configurations, beta' tuned to <E_cl> = E(T)
        lo, hi = -3, 3
        for _ in range(100):
            m = 0.5 * (lo + hi); w = np.exp(-(Ecl - Ecl.min()) * 10 ** m); w /= w.sum()
            if w @ Ecl > Et: lo = m
            else: hi = m
        w = np.exp(-(Ecl - Ecl.min()) * 10 ** m); w /= w.sum()
        pick = rng.choice(D, size=NS, p=w)
        ens[f"gibbsE_{fac}Tstar"] = [(1.0, unit(i)) for i in pick]
        dd = np.abs(Ecl - Et); win = np.where(dd <= max(2 * h, dd.min()) + 1e-9)[0]
        pick = rng.choice(win, size=NS)
        ens[f"microE_{fac}Tstar"] = [(1.0, unit(i)) for i in pick]
        Ts.append(T)
    keep = os.environ.get("ENSKEEP")
    if keep:
        ks = keep.split(","); ens = {k: v for k, v in ens.items() if k in ks}
    return ens, Ts


tl = np.linspace(0, 20, 201)
tc = np.linspace(0, 20, 241)
tw = [float(x) for x in os.environ.get('TWS', '0,2,8,16').split(',')]
res = analyze(f"Z2_{Lx}x{Ly}_h{h}", H, ops_local, ops_corr, ensembles, [], tl, tc, tw,
              extra={"model": "Z2 dual TFIM", "Lx": Lx, "Ly": Ly, "J": J, "h": h})
dump(res, f"out_z2_{Lx}x{Ly}_h{h}.json")
