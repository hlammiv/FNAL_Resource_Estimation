import sys, json, numpy as np
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[2]))   # repository root
from estimates import ch07_curved_space as c
res = {}
for lam, g in [(1, 1), (3, 0.5), (1, 2)]:
    for h in [0.05, 0.025, 0.0125, 0.001]:
        s, f = c.verlet_proxy(h, lam, 0.31, g)
        res[f"{lam},{g},{h}"] = dict(ns=s.tolist(), nf=f.tolist())
        print(lam, g, h, np.round(f, 4), flush=True)
json.dump(res, open(HERE / "proxy_full.json", "w"))
