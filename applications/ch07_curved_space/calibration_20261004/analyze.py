import json, glob, numpy as np
D = str(__import__("pathlib").Path(__file__).resolve().parent) + "/"
H = ["0.05", "0.025", "0.0125"]
SETS = [(1, 1), (3, 0.5), (1, 2)]


def rich(o1, o2, o4):
    ex = o4 + (o4 - o2) / 3
    ratio = (o1 - o2) / (o2 - o4) if abs(o2 - o4) > 1e-15 else float("nan")
    return ex, ratio


out = {}
pr = json.load(open(D + "proxy_red.json"))
for enc in ["fixed", "chi"]:
    for lam, g in SETS:
        f = D + f"exact_{enc}_l{float(lam)}_g{float(g)}.json"
        try:
            r = json.load(open(f))
        except FileNotFoundError:
            continue
        if not all(h in r["runs"] for h in H):
            continue
        print(f"\n== exact enc={enc} lam={lam} g={g}")
        print("  t      a    n_f(k1)@dt   rel_err(dt)  rel_err(dt/2) ratio | n_f(k0) relerr | n_s(k1) relerr(dt) | phibar")
        for i, rec in enumerate(r["runs"][H[0]]):
            o = [r["runs"][h][i] for h in H]
            nf = [x["n_ferm"][1] for x in o]
            ex, ra = rich(*nf)
            n0 = [x["n_ferm"][0] for x in o]
            ex0, _ = rich(*n0)
            ns = [x["n_scalar"] for x in o]
            exs, _ = rich(*ns)
            print(f"{rec['t']:7.1f} {rec['a']:5.2f}  {nf[0]:.5f}   {abs(nf[0]-ex)/ex:.2e}   {abs(nf[1]-ex)/ex:.2e}  {ra:5.2f} | "
                  f"{n0[0]:.4f} {abs(n0[0]-ex0)/max(ex0,1e-9):.1e} | {ns[0]:.4f} {abs(ns[0]-exs)/max(abs(exs),1e-9):.1e} | {rec['phibar']:+.4f}")
            out[f"{enc},{lam},{g},{rec['t']}"] = dict(nf=nf, ex=ex, rel=abs(nf[0] - ex) / ex, ratio=ra, n0=n0, ns=ns)
    print()
for lam, g in SETS:
    print(f"== proxy (reduced lattice) lam={lam} g={g}")
    ref = pr[f"{lam},{g},0.001"]
    for i, rec in enumerate(pr[f"{lam},{g},0.05"]):
        o = [pr[f"{lam},{g},{h}"][i]["n_ferm"][1] for h in H]
        ex, ra = rich(*o)
        print(f"{rec['t']:7.1f}  n_f={o[0]:.5f} rel_err(dt) vs 1e-3: {abs(o[0]-ref[i]['n_ferm'][1])/ref[i]['n_ferm'][1]:.2e} "
              f"rich {abs(o[0]-ex)/ex:.2e} ratio {ra:5.2f}  n_s={rec['n_scalar']:.4f}")
json.dump(out, open(D + "summary.json", "w"), indent=1)
