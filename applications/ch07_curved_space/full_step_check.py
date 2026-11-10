#!/usr/bin/env python3
"""Full-window rerun of the Ch. 7 exact step check on the reduced instance (U2, 2026-10-04).

The chapter's 2033 step error, '~2% (scaled from a reduced instance)', is the 5^3 mean-field proxy error
(Assumptions.proxy_step_err_2033 = 0.016, verlet_proxy) times the largest exact-to-proxy ratio on a 3-site ring
(exact_over_proxy_step_2033 = 1.23, exact_reduced_step): intermediate step_err_extrapolated = 0.020. The test
test_exact_reduced_step_check reruns one coupling set to t = 12.5 only. This driver reruns the whole study:
the three coupling sets (lambda, g) = (1, 1), (3, 0.5), (1, 2) at Phi0 = 0.31, m_psi = 0.1, the chapter's a(t),
K = 16, comoving register, at m_phi dt = 0.05, 0.025, 0.0125, read at eight times from t = 12.5 to 1000/m_phi.
At each read time it reports, for the fermion shell k = 2 pi/(3b):
  n_exact at the three steps, the Richardson ratio (n1 - n2)/(n2 - n4) (4.0 = clean second order),
  the exact step error |n(dt)/n_Richardson - 1| at m_phi dt = 0.05,
  the proxy step error on the same ring (verlet_proxy at the ring's two |k|, Richardson from the same three steps),
  and their ratio.
The summary compares the largest exact error, the largest ratio and the extrapolated 5^3 error with the constants
recorded in estimates/ch07_curved_space.py (0.0031, 1.23, 0.020). The read grid of the 2026-10-04 run is recorded
only as 'eight times from 12.5 to 1000'; this driver uses a doubling grid plus t = 1000.
Cost: about 13 ms per exact step on one core, so 2e4 + 4e4 + 8e4 steps per set, ~35 min per set, ~1.7 h in all
(OMP_NUM_THREADS=2). --quick runs t = 12.5 only (~1 min) and reproduces the test's numbers.
Usage: python applications/ch07_curved_space/full_step_check.py [--quick] [--out full_step_check.json]"""
import argparse, json, math, sys, time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from estimates import ch07_curved_space as c   # noqa: E402

SETS = ((1.0, 1.0), (3.0, 0.5), (1.0, 2.0))
HS = (0.05, 0.025, 0.0125)
REC = (12.5, 25.0, 50.0, 100.0, 200.0, 400.0, 800.0, 1000.0)
L, B, PHI0, M_PSI = 3, 3.46, 0.31, 0.1
KS = [0.0, 2 * math.pi / (L * B)]              # the ring's two |k| shells; index 1 is the shell the check reads


def richardson(n1, n2, n4):
    return n4 + (n4 - n2) / 3


def main(rec, out):
    a = c.Assumptions()
    rows = []
    for lam, g in SETS:
        t0 = time.time()
        o = c.exact_reduced_step(lam, g, list(HS), "chi", T_m=max(rec), rec=rec)
        for i, t in enumerate(rec):
            n = [o[h][i]["n_ferm"][1] for h in HS]
            ex = richardson(*n)
            ratio = (n[0] - n[1]) / (n[1] - n[2]) if n[1] != n[2] else float("nan")
            err_exact = abs(n[0] / ex - 1)
            p = [float(c.verlet_proxy(h, lam, PHI0, g, m_psi=M_PSI, ks=KS, T_m=t)[1][1]) for h in HS]
            err_proxy = abs(p[0] / richardson(*p) - 1)
            rows.append(dict(lam=lam, g=g, t=t, n_exact=n, n_proxy=p, richardson_ratio=ratio,
                             step_err_exact=err_exact, step_err_proxy=err_proxy,
                             exact_over_proxy=(err_exact / err_proxy if err_proxy else float("nan")),
                             n_scalar=o[HS[0]][i]["n_scalar"], phibar=o[HS[0]][i]["phibar"]))
            r = rows[-1]
            print(f"set ({lam:g}, {g:g}) t = {t:6.1f}: n_k {n[0]:.5f} ratio {ratio:5.2f} exact err {err_exact:.2e} "
                  f"proxy err {err_proxy:.2e} exact/proxy {r['exact_over_proxy']:.2f}", flush=True)
        print(f"set ({lam:g}, {g:g}) done in {time.time() - t0:.0f} s", flush=True)
    worst = max(rows, key=lambda r: r["step_err_exact"])
    worst_ratio = max(rows, key=lambda r: r["exact_over_proxy"])
    extrap = float(a.proxy_step_err_2033.lo) * worst_ratio["exact_over_proxy"]
    summary = dict(max_step_err_exact=worst["step_err_exact"], at=(worst["lam"], worst["g"], worst["t"]),
                   recorded_exact_step_err=float(a.exact_step_err_reduced_2033.lo),
                   max_exact_over_proxy=worst_ratio["exact_over_proxy"],
                   at_ratio=(worst_ratio["lam"], worst_ratio["g"], worst_ratio["t"]),
                   recorded_exact_over_proxy=float(a.exact_over_proxy_step_2033.lo),
                   step_err_extrapolated=extrap,
                   recorded_step_err_extrapolated=float(a.proxy_step_err_2033.lo) * float(a.exact_over_proxy_step_2033.lo))
    print("summary:", json.dumps(summary, indent=1))
    Path(out).write_text(json.dumps(dict(rec=rec, hs=HS, sets=SETS, rows=rows, summary=summary), indent=1))
    return summary


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true", help="t = 12.5 only (the test's window)")
    ap.add_argument("--out", default="full_step_check.json")
    args = ap.parse_args()
    main((12.5,) if args.quick else REC, args.out)
