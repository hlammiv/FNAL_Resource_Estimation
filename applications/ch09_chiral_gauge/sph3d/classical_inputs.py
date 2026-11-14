"""Classical inputs for the quantum cost model at one (N, beta_L): naive (Kubo) <dN^2>(t) and Re G(chi, t) =
<cos(chi dN)> from all replicas of runs/c_N{N}_b{beta}_s*.json (plain tag), with jackknife errors."""
import json, glob, math, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))


def inputs(N, beta, times=(5.0, 10.0, 20.0, 40.0), chis=(1.5, 2.5, 3.5), key="Nnaive"):
    rows = []
    for f in sorted(glob.glob(os.path.join(HERE, "runs", f"c_N{N}_b{beta:g}_s*.json"))):
        b = os.path.basename(f)
        if not b.split("_s")[1].split(".json")[0].isdigit():
            continue
        d = json.load(open(f))
        X = d[key] if key == "Nnaive" else d["Ncool"][2]
        rows += [np.asarray(x) for x in X]
    out = {"N": N, "beta_L": beta, "replicas": len(rows), "S": {}, "S_err": {}, "G": {}}
    for t in times:
        l = int(round(t))
        per = [(x[l:] - x[:-l])[::max(1, l // 2)] for x in rows]
        allv = np.concatenate(per)
        S = float((allv ** 2).mean())
        jk = [float(np.concatenate(per[:i] + per[i + 1:]) .__pow__(2).mean()) for i in range(len(per))]
        n = len(per)
        out["S"][t] = S
        out["S_err"][t] = math.sqrt((n - 1) / n * sum((j - np.mean(jk)) ** 2 for j in jk))
        out["G"][t] = {c: float(np.cos(c * allv).mean()) for c in chis}
    return out


if __name__ == "__main__":
    N, beta = int(sys.argv[1]), float(sys.argv[2])
    r = inputs(N, beta)
    print(json.dumps(r, indent=1))
