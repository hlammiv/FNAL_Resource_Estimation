"""Response of the naive (Kubo) N_CS to the cooled one: regression slope Z of dN_naive on dN_cool over windows t,
and the share of the naive <dN^2> it explains, Z^2 <dN_cool^2> / <dN_naive^2>."""
import json, glob, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))


def zfac(N, beta, t=40, k=2):
    a, b = [], []
    for f in sorted(glob.glob(os.path.join(HERE, "runs", f"c_N{N}_b{beta:g}_s*.json"))):
        if not os.path.basename(f).split("_s")[1].split(".json")[0].isdigit():
            continue
        d = json.load(open(f))
        for xn, xc in zip(d["Nnaive"], d["Ncool"][k]):
            xn, xc = np.asarray(xn), np.asarray(xc)
            a.append((xn[t:] - xn[:-t])[::t // 2]); b.append((xc[t:] - xc[:-t])[::t // 2])
    if not a:
        return None
    a, b = np.concatenate(a), np.concatenate(b)
    Z = float((a * b).mean() / (b * b).mean())
    return {"Z": Z, "explained": Z * Z * float((b * b).mean() / (a * a).mean()), "corr": float(np.corrcoef(a, b)[0, 1])}


if __name__ == "__main__":
    for N, beta in ((4, 2.5), (5, 2.5), (5, 3), (6, 2), (6, 2.5), (6, 3), (8, 2.5), (12, 2.5)):
        print(N, beta, {t: zfac(N, beta, t) for t in (20, 80)})
