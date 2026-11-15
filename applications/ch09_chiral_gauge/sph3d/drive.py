"""python3 drive.py N beta B t_max dt_meas seed [depth_base] [tag] -> runs/c_N{N}_b{beta}_s{seed}{tag}.json
depths = depth_base x (1/2, 1, 2), default depth_base = max(4, beta^2/8) (flow time max(1 a^2, 0.5/(g^2T)^2))."""
import sys, json, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import su2rt
N, beta, B, t_max, dtm, seed = int(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5]), int(sys.argv[6])
db = float(sys.argv[7]) if len(sys.argv) > 7 and sys.argv[7] != "-" else max(4.0, beta**2 / 8)
tag = sys.argv[8] if len(sys.argv) > 8 else ""
t0 = time.time()
r = su2rt.run(N, beta, B=B, t_max=t_max, dt_meas=dtm, depths=(db / 2, db, 2 * db), seed=seed)
r["cpu_s"] = time.time() - t0
RUNS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "runs"); os.makedirs(RUNS, exist_ok=True)
json.dump(r, open(os.path.join(RUNS, f"c_N{N}_b{beta:g}_s{seed}{tag}.json"), "w"))
print(N, beta, seed, tag, round(r["cpu_s"]), r["energy_drift_rel"], r["gauss_max"], flush=True)
