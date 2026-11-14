"""Summary numbers for result.json and the model's records: finite-volume ratios at each beta_L (4^3 and 5^3 against
the weighted mean of N >= 6), topological share at 4^3, spacing comparison with Moore-Rummukainen hep-ph/9906259."""
import json, math, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from zfactor import zfac
R = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "kappa_table.json")))
MR = {"table1_beta8": {3.0: (0.0023, 0.0016), 3.5: (0.029, 0.008), 4.0: (0.112, 0.012), 4.5: (0.276, 0.015), 5.0: (0.545, 0.035),
                       6.0: (1.18, 0.04), 7.0: (1.40, 0.05), 8.0: (1.61, 0.07), 10.0: (1.68, 0.03)},
      "table2_large_volume_imp": {4: (2.123, 0.075), 6: (1.800, 0.055), 8: (1.683, 0.032), 12: (1.336, 0.036), 16: (1.123, 0.046),
                                  24: (0.853, 0.031), 32: (0.626, 0.030)},
      "table2_large_volume_tree": {4: (3.78, 0.12), 6: (2.66, 0.08), 8: (2.25, 0.04), 12: (1.62, 0.04), 16: (1.30, 0.05)}}


def row(N, b, tag=""):
    for r in R:
        if r["N"] == N and r["beta_L"] == b and r["tag"] == tag:
            return r


def wmean(vals):
    w = [1 / e ** 2 for v, e in vals if e > 0]
    m = sum(v / e ** 2 for v, e in vals if e > 0) / sum(w)
    return m, 1 / math.sqrt(sum(w))


out = {"finite_volume": {}, "topological_share_4cube": {}, "validation_beta_L_8.6": {}, "spacing": {}}
for b in (2.0, 2.5, 3.0):
    big = [row(N, b) for N in (6, 8, 10, 12, 16) if row(N, b)]
    d = {}
    for key, k in (("naive", "naive"), ("cooled", "cool2")):
        m, e = wmean([(r[k]["kappa_tree"], r[k]["kappa_tree_err"]) for r in big])
        d[key] = {"large_volume_kappa_tree": (m, e), "large_volume_kappa_imp": (m * ((b - 0.6) / b) ** 4, e * ((b - 0.6) / b) ** 4)}
        for N in (3, 4, 5):
            r = row(N, b)
            if r:
                v, ve = r[k]["kappa_tree"], r[k]["kappa_tree_err"]
                d[key][f"ratio_N{N}"] = (v / m, (v / m) * math.sqrt((ve / v) ** 2 + (e / m) ** 2) if v > 0 else None)
                d[key][f"kappa_tree_N{N}"] = (v, ve)
                d[key][f"kappa_imp_N{N}"] = (r[k]["kappa_imp"], r[k]["kappa_imp_err"])
    out["finite_volume"][b] = d
    out["topological_share_4cube"][b] = zfac(4, b, 40)
for N, L in ((8, 4.0), (10, 5.0), (12, 6.0), (16, 8.0)):
    r = row(N, 8.6)
    if r:
        out["validation_beta_L_8.6"][N] = {"L_g2T_imp": L, "ours_cool1_imp": (r["cool1"]["kappa_imp"], r["cool1"]["kappa_imp_err"]),
                                          "ours_cool2_imp": (r["cool2"]["kappa_imp"], r["cool2"]["kappa_imp_err"]),
                                          "MR_table1": MR["table1_beta8"][L], "replicas": r["replicas"], "t_max": r["t_max"]}
r4 = row(4, 3.0)
out["spacing"] = {"kappa_imp_cool2_4cube_beta3": r4["cool2"]["kappa_imp"], "kappa_tree_cool2_4cube_beta3": r4["cool2"]["kappa_tree"],
                  "ratio_imp_to_MR_beta8": r4["cool2"]["kappa_imp"] / 1.683, "ratio_tree_to_MR_beta8_tree": r4["cool2"]["kappa_tree"] / 2.25,
                  "note": "classical lattice rates have no continuum limit (Gamma ~ a, Moore-Rummukainen); the reference is MR's a = 1/(2 g^2 T) value"}
out["MR"] = MR
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "summary.json"), "w"), indent=1)
print(json.dumps(out["finite_volume"], indent=0)[:3000]); print(json.dumps(out["topological_share_4cube"])); print(json.dumps(out["validation_beta_L_8.6"])); print(json.dumps(out["spacing"]))
