"""Generic exact-diagonalization analysis of quench thermalization vs canonical ensemble.

Inputs: dense real-symmetric H on a physical (gauge-invariant) sector, Hermitian real ops (sparse or dense,
same basis), initial-state ensembles (lists of real/complex vectors with weights).
Outputs (all lattice units, a = 1):
  canonical E(T), local <O>_T, connected two-time correlators C_T(t) = Re<A(t)A(0)>_T - <A>^2,
  quench: per-ensemble energy mean/width, local O(t), two-time correlator at waiting times t_w,
  diagonal-ensemble correlator, thermalization times.
"""
import numpy as np, scipy.linalg as sla, scipy.sparse as sp, time, json


def rc(M, X):
    """real matrix times complex matrix without promoting M."""
    if np.iscomplexobj(X):
        return M @ X.real + 1j * (M @ X.imag)
    return M @ X


def eig(H):
    t0 = time.time()
    E, V = sla.eigh(H, driver="evd", overwrite_a=True)
    print(f"  eigh D={len(E)} {time.time()-t0:.1f}s", flush=True)
    return E, V


def to_eig(V, A):
    AV = A @ V if sp.issparse(A) else A @ V
    return V.T @ AV


def canon(E, T):
    w = np.exp(-(E - E[0]) / T)
    return w / w.sum()


def T_of_energy(E, Etarget, lo=1e-3, hi=1e3):
    """canonical temperature with <H>_T = Etarget (bisection on log T)."""
    if Etarget <= E[0] or Etarget >= E.mean():
        return np.nan
    a, b = np.log(lo), np.log(hi)
    for _ in range(200):
        m = 0.5 * (a + b)
        if canon(E, np.exp(m)) @ E < Etarget:
            a = m
        else:
            b = m
    return float(np.exp(0.5 * (a + b)))


def spectral_corr(E, Aeig, p, times, block=2048):
    """C(t) = sum_nm p_n |A_nm|^2 cos((E_m-E_n) t) - (sum_n p_n A_nn)^2, by exact accumulation over omega bins."""
    D = len(E)
    span = E[-1] - E[0]
    nb = 400001
    edges = np.linspace(-span - 1e-9, span + 1e-9, nb + 1)
    hist = np.zeros(nb)
    for i in range(0, D, block):
        W = p[i:i + block, None] * Aeig[i:i + block, :] ** 2
        om = E[None, :] - E[i:i + block, None]
        idx = ((om.ravel() - edges[0]) * (nb / (edges[-1] - edges[0]))).astype(np.int64)
        np.clip(idx, 0, nb - 1, out=idx)
        hist += np.bincount(idx, weights=W.ravel(), minlength=nb)
    cen = 0.5 * (edges[1:] + edges[:-1])
    nz = hist != 0
    C = np.cos(np.outer(times, cen[nz])) @ hist[nz]
    mean = p @ np.diag(Aeig)
    return C - mean ** 2


def t_cross(times, y, level):
    """first time y(t) <= level (linear interpolation); nan if never."""
    for i in range(1, len(times)):
        if y[i] <= level < y[i - 1] or y[i] <= level:
            if y[i - 1] == y[i]:
                return float(times[i])
            f = (y[i - 1] - level) / (y[i - 1] - y[i])
            return float(times[i - 1] + f * (times[i] - times[i - 1]))
    return float("nan")


def quench_local(E, V, c, ops_comp, times):
    """<psi(t)|O|psi(t)> for psi(0) with eigen-coefficients c; ops in computational basis (sparse)."""
    ph = np.exp(-1j * np.outer(E, times)) * c[:, None]
    Psi = rc(V, ph)                             # D x Nt
    out = {}
    for k, O in ops_comp.items():
        out[k] = np.real(np.sum(np.conj(Psi) * (O @ Psi), axis=0))
    return out


def quench_corr(E, Aeig, a, times):
    """connected C(t) = Re <a| e^{iHt} A e^{-iHt} A |a> - <a|A(t)|a><a|A|a>, a in eigenbasis."""
    b = rc(Aeig, a)
    ph = np.exp(-1j * np.outer(E, times))
    Mb = rc(Aeig, ph * b[:, None])
    Ma = rc(Aeig, ph * a[:, None])
    left = np.conj(ph * a[:, None])
    full = np.sum(left * Mb, axis=0)
    At = np.real(np.sum(left * Ma, axis=0))
    A0 = np.real(np.vdot(a, b))
    return np.real(full), At, A0


def analyze(name, H, ops_local, ops_corr, ensembles, T_targets, t_local, t_corr, t_waits, extra=None):
    """ensembles: dict ename -> list of (weight, vector in sector basis). ops_corr: dict -> Hermitian real op."""
    res = {"name": name, "D": H.shape[0]}
    if extra:
        res.update(extra)
    Hd = H.toarray() if sp.issparse(H) else np.array(H, dtype=float)
    E, V = eig(Hd)
    del Hd
    res["E0"] = float(E[0]); res["Emax"] = float(E[-1]); res["Emid"] = float(E.mean())
    if callable(ensembles):
        ensembles, T_extra = ensembles(E)
        T_targets = sorted(set(list(T_targets) + list(T_extra)))
    loc_eig = {k: np.einsum("in,in->n", V, (O @ V)) for k, O in ops_local.items()}
    Aeig = {k: to_eig(V, A) for k, A in ops_corr.items()}
    # canonical
    can = {}
    for T in T_targets:
        p = canon(E, T)
        d = {"T": T, "E": float(p @ E), "sigE": float(np.sqrt(p @ E ** 2 - (p @ E) ** 2))}
        d["local"] = {k: float(p @ v) for k, v in loc_eig.items()}
        d["corr"] = {}
        for k, A in Aeig.items():
            C = spectral_corr(E, A, p, t_corr)
            d["corr"][k] = {"C": C.tolist(), "t_e": t_cross(t_corr, C / C[0], np.exp(-1)),
                            "t_half": t_cross(t_corr, C / C[0], 0.5)}
        can[f"{T:.4g}"] = d
    res["canonical"] = can
    Tg = np.geomspace(0.1, 10, 200)
    res["curves"] = {"T": Tg.tolist(), "E": [float(canon(E, T) @ E) for T in Tg],
                     "local": {k: [float(canon(E, T) @ v) for T in Tg] for k, v in loc_eig.items()}}
    res["t_local"] = list(map(float, t_local)); res["t_corr"] = list(map(float, t_corr))
    res["t_waits"] = list(map(float, t_waits))
    # quench ensembles
    qres = {}
    for en, samples in ensembles.items():
        t0 = time.time()
        W = np.array([w for w, _ in samples]); W = W / W.sum()
        cs = [rc(V.T, v) for _, v in samples]
        cs = [c / np.linalg.norm(c) for c in cs]
        Es = np.array([np.real(np.vdot(c, E * c)) for c in cs])
        sE = np.array([np.sqrt(max(np.real(np.vdot(c, E ** 2 * c)) - e ** 2, 0)) for c, e in zip(cs, Es)])
        Emean = float(W @ Es)
        Ttgt = T_of_energy(E, Emean)
        pdiag = sum(w * np.abs(c) ** 2 for w, c in zip(W, cs))
        d = {"n": len(samples), "E_mean": Emean, "E_spread_samples": float(np.sqrt(W @ (Es - Emean) ** 2)),
             "sigE_within_mean": float(W @ sE),
             "sigE_total": float(np.sqrt(pdiag @ E ** 2 - (pdiag @ E) ** 2)),
             "T_matched": Ttgt}
        # local relaxation (ensemble-averaged and per-sample spread)
        loc = {k: np.zeros(len(t_local)) for k in ops_local}
        loc_s = {k: [] for k in ops_local}
        for w, c in zip(W, cs):
            o = quench_local(E, V, c, ops_local, t_local)
            for k in o:
                loc[k] += w * o[k]; loc_s[k].append(o[k])
        d["local"] = {k: loc[k].tolist() for k in loc}
        d["local_sample_sd"] = {k: np.std(np.array(loc_s[k]), axis=0).tolist() for k in loc}
        d["local_diag"] = {k: float(pdiag @ v) for k, v in loc_eig.items()}
        if not np.isnan(Ttgt):
            pc = canon(E, Ttgt)
            d["local_canon_at_T"] = {k: float(pc @ v) for k, v in loc_eig.items()}
        # two-time correlators
        d["corr"] = {}
        for k, A in Aeig.items():
            Cdiag = spectral_corr(E, A, pdiag, t_corr)
            rows, rows_ens = [], []
            for tw in t_waits:
                C = np.zeros(len(t_corr)); F = np.zeros(len(t_corr)); At_m = np.zeros(len(t_corr)); A0_m = 0.0
                for w, c in zip(W, cs):
                    a = np.exp(-1j * E * tw) * c
                    full, At, A0 = quench_corr(E, A, a, t_corr)
                    C += w * (full - At * A0); F += w * full; At_m += w * At; A0_m += w * A0
                rows.append(C.tolist()); rows_ens.append((F - At_m * A0_m).tolist())
            ent = {"C_tw": rows, "C_tw_ens": rows_ens, "C_diag": Cdiag.tolist(),
                   "t_e_tw_ens": [t_cross(t_corr, np.array(r) / r[0], np.exp(-1)) for r in rows_ens],
                   "t_e_diag": t_cross(t_corr, Cdiag / Cdiag[0], np.exp(-1)),
                   "t_e_tw": [t_cross(t_corr, np.array(r) / r[0], np.exp(-1)) for r in rows]}
            if not np.isnan(Ttgt):
                Cc = spectral_corr(E, A, canon(E, Ttgt), t_corr)
                ent["C_canon_at_T"] = Cc.tolist()
                ent["t_e_canon"] = t_cross(t_corr, Cc / Cc[0], np.exp(-1))
            d["corr"][k] = ent
        qres[en] = d
        print(f"  ensemble {en}: n={len(samples)} E={Emean:.4f} T*={Ttgt:.4g} ({time.time()-t0:.1f}s)", flush=True)
    res["quench"] = qres
    return res


def dump(res, fn):
    def conv(o):
        if isinstance(o, (np.floating, np.integer)):
            return o.item()
        if isinstance(o, np.ndarray):
            return o.tolist()
        raise TypeError(type(o))
    with open(fn, "w") as f:
        json.dump(res, f, default=conv)
