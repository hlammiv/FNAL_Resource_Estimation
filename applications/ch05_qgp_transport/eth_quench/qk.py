"""Krylov / typicality version of qa.analyze for sectors too large for dense ED.
Canonical values by dynamical typicality (R random Gaussian vectors in the sector, e^{-beta H/2}), real time by
scipy expm_multiply on a uniform grid."""
import numpy as np, scipy.sparse as sp, time
from scipy.sparse.linalg import expm_multiply, eigsh


def evolve_grid(H, v, t1, n):
    """v(t) on t = linspace(0, t1, n), returns D x n."""
    return expm_multiply(-1j * H, v.astype(complex), start=0, stop=t1, num=n, endpoint=True).T


def imag_grid(H, v, betas_half):
    out = []; cur = v.copy(); b0 = 0.0
    for b in betas_half:
        cur = expm_multiply(-(b - b0) * H, cur); b0 = b
        out.append(cur.copy())
    return out


def corr_from(H, A, v, t1, n):
    """connected Re<v|A(t)A|v> - <A(t)><A>, unnormalized by <v|v> (caller normalizes)."""
    Av = A @ v
    P1 = evolve_grid(H, Av, t1, n); P2 = evolve_grid(H, v, t1, n)
    full = np.sum(np.conj(P2) * (A @ P1), axis=0)
    nv = np.vdot(v, v).real
    At = np.real(np.sum(np.conj(P2) * (A @ P2), axis=0)) / nv
    A0 = np.real(np.vdot(v, Av)) / nv
    return np.real(full) / nv, At, A0


def t_cross(times, y, level):
    for i in range(1, len(times)):
        if y[i] <= level:
            f = (y[i - 1] - level) / (y[i - 1] - y[i]) if y[i - 1] != y[i] else 1
            return float(times[i - 1] + f * (times[i] - times[i - 1]))
    return float("nan")


def analyze_k(name, H, ops_local, ops_corr, make_ensembles, t_local, n_local, t_corr, n_corr, t_waits, R=6,
              Tgrid=None, seed=0, extra=None):
    rng = np.random.default_rng(seed)
    D = H.shape[0]
    res = {"name": name, "D": D, "method": "Krylov+typicality", "R": R}
    if extra: res.update(extra)
    E0 = float(eigsh(H, k=1, which="SA")[0][0]); res["E0"] = E0
    Hs = (H - E0 * sp.identity(D, format="csr")).tocsr()
    Tgrid = np.array(Tgrid if Tgrid is not None else np.geomspace(0.3, 20, 120))
    bh = np.sort(1 / (2 * Tgrid))            # half-betas ascending
    t0 = time.time()
    rs = [rng.standard_normal(D) for _ in range(R)]
    states = [imag_grid(Hs, r, bh) for r in rs]   # list over r of list over beta
    Ecurve = []
    for k, b in enumerate(bh):
        num = sum(np.vdot(st[k], H @ st[k]).real for st in states); den = sum(np.vdot(st[k], st[k]).real for st in states)
        Ecurve.append(num / den)
    Ts = 1 / (2 * bh)
    order = np.argsort(Ts); TsS, EcS = Ts[order], np.array(Ecurve)[order]
    res["E_of_T"] = {"T": TsS.tolist(), "E": EcS.tolist()}
    print(f"  typicality E(T) {time.time()-t0:.1f}s", flush=True)

    def T_of(Et):
        if Et < EcS[0] or Et > EcS[-1]: return float("nan")
        return float(np.interp(Et, EcS, TsS))

    def thermal_state(T):
        k = int(np.argmin(np.abs(Ts - T)))
        return [st[k] for st in states], float(Ts[k])
    E_of = lambda T: float(np.interp(T, TsS, EcS))
    ensembles, Textra = make_ensembles(T_of, E_of)
    can = {}
    for T in Textra:
        vs, Tg = thermal_state(T)
        den = sum(np.vdot(v, v).real for v in vs)
        d = {"T": Tg, "E": float(sum(np.vdot(v, H @ v).real for v in vs) / den),
             "local": {k: float(sum(np.vdot(v, O @ v).real for v in vs) / den) for k, O in ops_local.items()}, "corr": {}}
        for k, A in ops_corr.items():
            parts = [corr_from(H, A, v, t_corr, n_corr) for v in vs]; ws_ = [np.vdot(v, v).real / den for v in vs]
            F = sum(w * p[0] for w, p in zip(ws_, parts)); Atm = sum(w * p[1] for w, p in zip(ws_, parts)); A0m = sum(w * p[2] for w, p in zip(ws_, parts))
            C = F - Atm * A0m
            tt = np.linspace(0, t_corr, n_corr)
            d["corr"][k] = {"C": C.tolist(), "t_e": t_cross(tt, C / C[0], np.exp(-1))}
        can[f"{Tg:.4g}"] = d
        print(f"  canonical T={Tg:.3f} {time.time()-t0:.1f}s", flush=True)
    res["canonical"] = can
    q = {}
    for en, samples in ensembles.items():
        W = np.array([w for w, _ in samples]); W = W / W.sum()
        Es = np.array([np.vdot(v, H @ v).real / np.vdot(v, v).real for _, v in samples])
        Em = float(W @ Es)
        d = {"n": len(samples), "E_mean": Em, "T_matched": T_of(Em), "local": {}, "corr": {}}
        loc = {k: np.zeros(n_local) for k in ops_local}
        for w, (_, v) in zip(W, samples):
            Vt = evolve_grid(H, v / np.linalg.norm(v), t_local, n_local)
            for k, O in ops_local.items():
                loc[k] += w * np.real(np.sum(np.conj(Vt) * (O @ Vt), axis=0))
        d["local"] = {k: x.tolist() for k, x in loc.items()}
        for k, A in ops_corr.items():
            rows, rows_ens = [], []
            for tw in t_waits:
                C = np.zeros(n_corr); F = np.zeros(n_corr); Atm = np.zeros(n_corr); A0m = 0.0
                for w, (_, v) in zip(W, samples):
                    vw = v / np.linalg.norm(v)
                    if tw > 0: vw = expm_multiply(-1j * tw * H, vw.astype(complex))
                    full, At, A0 = corr_from(H, A, vw, t_corr, n_corr)
                    C += w * (full - At * A0); F += w * full; Atm += w * At; A0m += w * A0
                rows.append(C.tolist()); rows_ens.append((F - Atm * A0m).tolist())
            tt = np.linspace(0, t_corr, n_corr)
            d["corr"][k] = {"C_tw": rows, "C_tw_ens": rows_ens, "t_e_tw_ens": [t_cross(tt, np.array(r) / r[0], np.exp(-1)) for r in rows_ens], "t_e_tw": [t_cross(tt, np.array(r) / r[0], np.exp(-1)) for r in rows]}
        q[en] = d
        print(f"  ensemble {en} E={Em:.3f} T*={d['T_matched']:.3f} {time.time()-t0:.1f}s", flush=True)
    res["quench"] = q
    res["t_local"] = np.linspace(0, t_local, n_local).tolist(); res["t_corr"] = np.linspace(0, t_corr, n_corr).tolist()
    res["t_waits"] = list(t_waits)
    return res
