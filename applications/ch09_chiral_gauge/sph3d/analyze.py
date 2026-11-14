"""Analyze runs/c_*.json -> kappa_table.json.
kappa = Gamma/(alpha T)^4, Gamma = d<dN^2>/dt / V. Two conversions: tree (g^2 a T = 4/beta_L) and improved
(Moore-Rummukainen hep-ph/9906259 Eq. 12, a = 4/((beta_L - 0.6) g^2 T)). Errors: jackknife over replicas.
Definitions: 'naive' = Kubo integral of E.B along the real-time path (what the quantum counting-field readout
measures); 'cool_k' = cooled N_CS at flow depth k (topological-style, Moore-Turok)."""
import json, glob, math, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


def diffusion(rows, dtm, lo, hi):
    """rows: list of 1D arrays (replicas, any length). Pooled MSD slope over lags [lo, hi], jackknife over replicas."""
    lags = np.arange(max(1, int(lo / dtm)), int(hi / dtm) + 1)
    S = np.array([[((x[l:] - x[:-l]) ** 2).sum() if l < len(x) else 0.0 for l in lags] for x in rows])
    C = np.array([[max(len(x) - l, 0) for l in lags] for x in rows], float)
    s_ = lags * dtm
    A = np.vstack([s_, np.ones_like(s_)]).T

    def slope(Ssum, Csum):
        return np.linalg.lstsq(A, Ssum / Csum, rcond=None)[0][0]
    full = slope(S.sum(0), C.sum(0))
    B = len(rows)
    jk = np.array([slope(S.sum(0) - S[b], C.sum(0) - C[b]) for b in range(B)])
    err = math.sqrt((B - 1) / B * ((jk - jk.mean()) ** 2).sum())
    return full, err

HERE = os.path.dirname(os.path.abspath(__file__))
FAC = (4 * math.pi) ** 4


def lag_window(beta, t_max):
    lo, hi = (10.0, 60.0) if beta <= 4 else (20.0, 120.0)
    hi = min(hi, t_max / 8)
    lo = min(lo, hi / 3)
    return lo, hi


def kappa(slope, N, beta, improved):
    ga = 4.0 / (beta - 0.6) if improved else 4.0 / beta
    return slope * FAC / (ga ** 4 * N ** 3)


def load():
    groups = {}
    for f in sorted(glob.glob(os.path.join(HERE, "runs", "c_*.json"))):
        d = json.load(open(f))
        import re as _re
        tag = _re.match(r".*_s\d+(.*)\.json$", os.path.basename(f)).group(1).lstrip("_")
        groups.setdefault((d["N"], d["beta"], tag), []).append(d)
    return groups


def analyze(groups):
    rows = []
    for (N, beta, tag), ds in sorted(groups.items()):
        dtm = ds[0]["dt_meas"]
        tmax = max(d["t_max"] for d in ds)
        lo, hi = lag_window(beta, min(d["t_max"] for d in ds))
        out = {"N": N, "beta_L": beta, "tag": tag, "replicas": sum(d["B"] for d in ds), "t_max": tmax,
               "lag_window": (lo, hi), "L_g2T_tree": 4 * N / beta, "L_g2T_improved": 4 * N / (beta - 0.6) if beta > 0.6 else None,
               "flow_time_a2": ds[0]["flow_time_a2"], "Ecool_mean": ds[0]["Ecool_mean"],
               "energy_drift_rel": max(d["energy_drift_rel"] for d in ds), "gauss_max": max(d["gauss_max"] for d in ds),
               "cpu_s": sum(d.get("cpu_s", 0) for d in ds)}
        series = {"naive": [np.asarray(x) for d in ds for x in d["Nnaive"]]}
        for k in range(len(ds[0]["Ncool"])):
            series[f"cool{k}"] = [np.asarray(x) for d in ds for x in d["Ncool"][k]]
        for key, X in series.items():
            s, e = diffusion(X, dtm, lo, hi)
            out[key] = {"slope": s, "slope_err": e, "kappa_tree": kappa(s, N, beta, False), "kappa_tree_err": kappa(e, N, beta, False),
                        "kappa_imp": kappa(s, N, beta, True), "kappa_imp_err": kappa(e, N, beta, True)}
            # distribution of dN over windows of length t for the counting-field estimate (naive and cooled)
            dist = {}
            for tw in (2.0, 5.0, 10.0, 20.0, 40.0, 80.0, 160.0):
                l = int(round(tw / dtm))
                if l >= min(len(x) for x in X) // 2:
                    continue
                dN = np.concatenate([(x[l:] - x[:-l])[::max(1, l // 2)] for x in X])
                m2 = float((dN ** 2).mean()); m4 = float((dN ** 4).mean())
                dist[tw] = {"m2": m2, "m4": m4, "n": int(dN.size),
                            "G": {c: float(np.cos(c * dN).mean()) for c in (0.25, 0.5, 1.0, 1.5, 2.0, math.pi)},
                            "frac_ge_half": float((np.abs(dN) >= 0.5).mean())}
            out[key]["dist"] = dist
        rows.append(out)
    return rows


if __name__ == "__main__":
    rows = analyze(load())
    json.dump(rows, open(os.path.join(HERE, "kappa_table.json"), "w"), indent=1)
    print(f"{'N':>3} {'bL':>4} {'tag':>4} {'Lg2T_t':>6} {'Lg2T_i':>6} {'R':>4} {'k_naive_tree':>14} {'k_cool1_tree':>14} {'k_cool1_imp':>14} {'k_cool0':>8} {'k_cool2':>8} drift")
    for r in rows:
        c0, c1, c2, nv = r["cool0"], r["cool1"], r["cool2"], r["naive"]
        print(f"{r['N']:3d} {r['beta_L']:4g} {r['tag']:>4} {r['L_g2T_tree']:6.2f} {r['L_g2T_improved']:6.2f} {r['replicas']:4d} "
              f"{nv['kappa_tree']:8.3f}+-{nv['kappa_tree_err']:.3f} {c1['kappa_tree']:8.3f}+-{c1['kappa_tree_err']:.3f} "
              f"{c1['kappa_imp']:8.3f}+-{c1['kappa_imp_err']:.3f} {c0['kappa_tree']:8.3f} {c2['kappa_tree']:8.3f} {r['energy_drift_rel']:.1e}")
