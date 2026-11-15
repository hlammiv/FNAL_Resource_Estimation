"""Classical quench check: U = 1, e Gaussian at beta_init = beta/2 projected on Gauss; record the kinetic temperature
and plaquette energy vs t (lattice units) against the thermal values at beta."""
import sys, json, math, numpy as np
sys.path.insert(0, '.')
import su2rt as s
N, beta, B = int(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3])
rng = np.random.default_rng(5)
U = np.zeros((B, 3, N, N, N, 4)); U[..., 0] = 1
e = s.project_gauss(U, rng.normal(size=(B, 3, N, N, N, 3)) / math.sqrt(beta / 2))
F = None; out = []
for t in range(0, 41):
    kin = float((0.5 * (e * e).sum((1, 2, 3, 4, 5)) / (3 * N ** 3)).mean() * beta * 1.5)  # compared with the same quantity in the thermal ensemble ('thermal')
    pl = float((s.plaq_energy(U) / (3 * N ** 3)).mean())
    out.append((t, kin, pl))
    U, e, F = s.leapfrog(U, e, 0.2, 10, F)
Ut, et, Ft = s.thermalize(B, N, beta, rng, n_ref=30)
ref = (float((0.5 * (et * et).sum((1, 2, 3, 4, 5)) / (3 * N ** 3)).mean() * beta * 1.5), float((s.plaq_energy(Ut) / (3 * N ** 3)).mean()))
json.dump({"N": N, "beta": beta, "series": out, "thermal": ref}, open(f"quench_N{N}_b{beta:g}.json", "w"))
for r in out[::2]: print(r)
print("thermal", ref)
