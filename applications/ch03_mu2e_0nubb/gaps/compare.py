"""Compare this folder's regenerated outputs with the values the Ch. 3 model holds (estimates/ch03_mu2e_0nubb.py).
Exits nonzero if any output is missing or differs by more than its tolerance."""
import json, sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[2]))
from estimates import ch03_mu2e_0nubb as ch3
a = ch3.Assumptions()
v = lambda f: getattr(a, f).value if hasattr(getattr(a, f), "value") else getattr(a, f)
load = lambda n: json.loads((HERE / n).read_text())
gap = lambda n: (lambda d: d["E"][1] - d["E"][0])(load(n))
rows = []
for sp in ("f7p3", "f7p3p1", "pf"):
    for nuc, tag in (("48Ca", "ca48"), ("48Ti", "ti48")):
        rows.append((f"gap_0plus_{sp}_{tag}_mev", lambda sp=sp, nuc=nuc: gap(f"gap_kb3g_{sp}_{nuc}_J0.0.json"), 1e-4))
for (sp, nuc), rec in {("f7p3", "48Ca"): 6.0595, ("f7p3", "48Ti"): 4.0233, ("f7p3p1", "48Ca"): 5.6371,
                       ("f7p3p1", "48Ti"): 3.8668, ("pf", "48Ca"): 5.2749}.items():
    rows.append((f"GXPF1A {sp} {nuc} (record in the model comment)", lambda sp=sp, nuc=nuc, rec=rec: (gap(f"gap_gxpf1a_{sp}_{nuc}_J0.0.json"), rec), 1e-3))
def al27_first_reachable():
    d = load("gap_usdb_sd_27Al_J-1.0.json")
    return min(e - d["E"][0] for e, j in zip(d["E"][1:], d["J"][1:]) if j >= 2.49)
rows += [("gap_al27_usdb_mev", al27_first_reachable, 1e-4),
         ("gap_al27_usdb_same_j_mev", lambda: gap("gap_usdb_sd_27Al_J2.5.json"), 1e-4),
         ("sd_spin_expectation", lambda: load("al27_sd_readout.json")["S_p"], 1e-4),
         ("sd_element_expectation", lambda: load("al27_sd_readout.json")["O_SD"]["expectation"], 1e-4),
         ("direct_readout_variance", lambda: load("al27_sd_readout.json")["O_SD"]["variance"], 1e-4),
         ("projection_success", lambda: load("al27_sd_readout.json")["p_largest_component"], 1e-4),
         ("he4_d_fid99", lambda: (load("he4D_u1.0_s1.0_b10.0.json")["D_infid_0.01"], load("he4D_u1.0_s1.0_b2.0.json")["D_infid_0.01"]), 0),
         ("he4_d_fid999", lambda: (load("he4D_u1.0_s1.0_b2.0.json")["D_infid_0.001"], load("he4D_u1.0_s1.0_b10.0.json")["D_infid_0.001"]), 0),
         ("he4_d_dipole_1pct", lambda: max(load(f"he4D_u1.0_s1.0_b{b}.json")["D_Oerr_0.01"] for b in ("2.0", "10.0")), 0)]
bad = 0
for name, f, tol in rows:
    try:
        got = f()
    except (FileNotFoundError, KeyError) as e:
        print(f"MISSING  {name}: {e}"); bad += 1; continue
    if isinstance(got, tuple) and isinstance(got[1], float) and name.startswith("GXPF1A"):
        got, want = got
    else:
        want = v(name)
    ok = (abs(got - want) <= tol * max(1, abs(want))) if not isinstance(want, tuple) else tuple(got) == tuple(want)
    bad += not ok
    print(f"{'ok  ' if ok else 'DIFF'}     {name}: regenerated {got}, model {want}")
print("all match" if not bad else f"{bad} differ or missing")
sys.exit(1 if bad else 0)
