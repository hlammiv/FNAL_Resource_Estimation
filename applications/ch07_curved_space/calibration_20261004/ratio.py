"""The 2026-10-04 calibration exactly as computed then (recovered from the session record, 2026-11-11):
per coupling set, the largest exact step error over the read times (chi register, summary.json from analyze.py)
over the largest ring-proxy step error (proxy_red.json, Richardson from dt, dt/2, dt/4), and the full 5^3 proxy
error (proxy_full.json, against dt = 0.001, shells with n_k > 0.1). exact_over_proxy_step_2033 is the largest ratio.
Recorded output: ratios 0.973, 1.230, 0.518; full 0.0146, 0.0033, 0.0163."""
import json, sys
from pathlib import Path
import numpy as np
D = Path(__file__).resolve().parent
s = json.load(open(D / "summary.json")); p = json.load(open(D / "proxy_red.json")); f = json.load(open(D / "proxy_full.json"))
H = ["0.05", "0.025", "0.0125"]
out = {}
for lam, g in [(1, 1), (3, 0.5), (1, 2)]:
    ex = max(v["rel"] for k, v in s.items() if k.startswith(f"chi,{lam},{g},"))
    pr = []
    for i in range(len(p[f"{lam},{g},0.05"])):
        o = [p[f"{lam},{g},{h}"][i]["n_ferm"][1] for h in H]; e = o[2] + (o[2] - o[1]) / 3; pr.append(abs(o[0] - e) / e)
    ref = np.array(f[f"{lam},{g},0.001"]["nf"]); a = np.array(f[f"{lam},{g},0.05"]["nf"]); m = ref > 0.1
    full = float(np.max(abs(a / ref - 1)[m]))
    out[f"{lam},{g}"] = dict(exact=ex, proxy_ring=max(pr), ratio=ex / max(pr), proxy_full=full, extrap=full * ex / max(pr))
    print(lam, g, "exact %.2e proxy ring %.2e ratio %.3f full %.4f extrap %.4f" % (ex, max(pr), ex / max(pr), full, full * ex / max(pr)))
r = max(v["ratio"] for v in out.values()); pf = max(v["proxy_full"] for v in out.values())
print("largest ratio %.3f; largest full proxy %.4f; model's product %.4f; largest per-set extrap %.4f"
      % (r, pf, pf * r, max(v["extrap"] for v in out.values())))
json.dump(out, open(D / "ratio.json", "w"), indent=1)
