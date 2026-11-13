import numpy as np, scipy.sparse as sp, sys, json
from qk import analyze_k
Lx, Ly, h = int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3]); J = 1.0
NS = int(sys.argv[4]) if len(sys.argv) > 4 else 6
N = Lx * Ly
site = lambda x, y: (x % Lx) + Lx * (y % Ly)
bonds = [(site(x, y), site(x + 1, y), 'x') for y in range(Ly) for x in range(Lx)] + [(site(x, y), site(x, y + 1), 'y') for y in range(Ly) for x in range(Lx)]
full = np.arange(2 ** N, dtype=np.int64); reps = full[(full & 1) == 0]; D = len(reps); mask = 2 ** N - 1
idx = -np.ones(2 ** N, dtype=np.int64); idx[reps] = np.arange(D); idx[reps ^ mask] = np.arange(D)
spin = lambda s, p: 1 - 2 * ((s >> p) & 1)
zz = {b: spin(reps, b[0]) * spin(reps, b[1]) for b in bonds}
Ecl = -h * sum(zz.values()).astype(float)
Xop = {p: sp.csr_matrix((np.ones(D), (idx[reps ^ (1 << p)], np.arange(D))), shape=(D, D)) for p in range(N)}
H = (sp.diags(Ecl) - J * sum(Xop.values())).tocsr()
xp = np.array([p % Lx for p in range(N)]); kx = 2 * np.pi / Lx
ops_local = {"B": (sum(Xop.values()) / N).tocsr(), "Elink": sp.diags(sum(zz.values()) / len(bonds)).tocsr()}
ops_corr = {"B_kmin": sum(np.cos(kx * xp[p]) * Xop[p] for p in range(N)).tocsr(),
            "eps_kmin": (sum(np.cos(kx * xp[p]) * (-J) * Xop[p] for p in range(N)) + sp.diags(sum(0.5 * (np.cos(kx * xp[b[0]]) + np.cos(kx * xp[b[1]])) * (-h) * zz[b] for b in bonds))).tocsr(),
            "Eaniso_kmin": sp.diags(sum(np.cos(kx * xp[b[0]]) * (1 if b[2] == 'x' else -1) * zz[b] for b in bonds)).tocsr()}
rng = np.random.default_rng(12345)
def unit(i):
    v = np.zeros(D); v[i] = 1.0; return v
def ens(T_of, E_of):
    vac = idx[0]; Ts_ = T_of(Ecl[vac]); out = {"vac": [(1.0, unit(vac))]}; Ts = [Ts_]
    for fac in (1.5, 2.5):
        T = Ts_ * fac; Et = E_of(T)
        if Et >= 0: continue
        lo, hi = -3, 3
        for _ in range(100):
            m = 0.5 * (lo + hi); w = np.exp(-(Ecl - Ecl.min()) * 10 ** m); w /= w.sum()
            if w @ Ecl > Et: lo = m
            else: hi = m
        out[f"gibbsE_{fac}Tstar"] = [(1.0, unit(i)) for i in rng.choice(D, size=NS, p=w)]
        Ts.append(T)
    return out, Ts
res = analyze_k(f"Z2K_{Lx}x{Ly}_h{h}", H, ops_local, ops_corr, ens, 20.0, 201, 20.0, 241, [0.0, 4.0, 16.0],
                extra={"model": "Z2 dual TFIM (Krylov)", "Lx": Lx, "Ly": Ly, "J": J, "h": h})
json.dump(res, open(f"out_z2k_{Lx}x{Ly}_h{h}.json", "w"))
