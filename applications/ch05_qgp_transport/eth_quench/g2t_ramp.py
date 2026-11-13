"""2T open patch: cool below the E=0 product-state temperature T* by a finite-speed coupling ramp g^2: g2i -> g2f
(linear in 1/g^2 over time tau, starting from the E=0 product state = ground state at g^2 -> inf), then hold at g2f.
Reports energy after the ramp, matched canonical T (typicality / ED), and relaxation of local observables."""
import numpy as np, scipy.sparse as sp, sys, json, os
from scipy.sparse.linalg import expm_multiply, eigsh
Lx, Ly, g2f = int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3])
src = open("g2t.py").read().split("NS = int")[0]
sys.argv = ["g2t.py", str(Lx), str(Ly), "1.0"]
exec(src)                                  # builds HE, HB at g2 = 1.0 (HE ~ g2, HB ~ 1/g2)
HE1 = proj(sum(HE.values())); HB1 = proj(sum(HB.values()))
Hf = (g2f * HE1 + HB1 / g2f).tocsr()
D = Hf.shape[0]
plaq = ops_local["plaq_trF/2"]
v0 = vac / np.linalg.norm(vac)
# canonical E(T), plaquette(T): exact for small D, else typicality
if D <= 3000:
    E, V = np.linalg.eigh(Hf.toarray()); pl = np.einsum("in,in->n", V, plaq @ V)
    def can(T):
        w = np.exp(-(E - E[0]) / T); w /= w.sum(); return w @ E, w @ pl
else:
    rng = np.random.default_rng(0); E0 = eigsh(Hf, k=1, which="SA")[0][0]
    Hs = Hf - E0 * sp.identity(D)
    Ts = np.geomspace(0.2, 5, 80)[::-1]; rs = [rng.standard_normal(D) for _ in range(8)]
    cur = [r.copy() for r in rs]; b0 = 0; tab = []
    for T in Ts:
        b = 1 / (2 * T); cur = [expm_multiply(-(b - b0) * Hs, c) for c in cur]; b0 = b
        den = sum(c @ c for c in cur)
        tab.append((T, sum(c @ (Hf @ c) for c in cur) / den, sum(c @ (plaq @ c) for c in cur) / den))
    tab = np.array(tab[::-1])
    def can(T):
        return np.interp(T, tab[:, 0], tab[:, 1]), np.interp(T, tab[:, 0], tab[:, 2])
Tg = np.geomspace(0.2, 5, 400); Ecan = np.array([can(T)[0] for T in Tg])
Tof = lambda e: float(np.interp(e, Ecan, Tg)) if Ecan[0] <= e <= Ecan[-1] else float("nan")
out = {"Lx": Lx, "Ly": Ly, "g2f": g2f, "D": D, "E_vac_state": float(v0 @ (Hf @ v0)), "T_star": Tof(v0 @ (Hf @ v0)), "ramps": {}}
print("T* (sudden) =", out["T_star"], flush=True)
g2i = 8.0
for tau in [0.0, 1.0, 2.0, 4.0, 8.0, 16.0]:
    v = v0.astype(complex); n = max(int(tau / 0.05), 1) if tau > 0 else 0
    for k in range(n):
        s = (k + 0.5) / n; ig = (1 - s) / g2i + s / g2f; g2 = 1 / ig
        v = expm_multiply(-1j * (tau / n) * (g2 * HE1 + HB1 / g2), v)
    Er = float(np.real(np.vdot(v, Hf @ v))); T = Tof(Er)
    ts = np.linspace(0, 20, 81); vt = expm_multiply(-1j * Hf, v, start=0, stop=20, num=81, endpoint=True)
    pl_t = [float(np.real(np.vdot(x, plaq @ x))) for x in vt]
    out["ramps"][str(tau)] = {"E": Er, "T_eff": T, "plaq_canon": float(can(T)[1]) if T == T else None, "plaq_t": pl_t, "t": ts.tolist()}
    print(f"tau={tau}: E={Er:.4f} T_eff={T:.4f} plaq(0)={pl_t[0]:.4f} plaq late={np.mean(pl_t[40:]):.4f} canon={out['ramps'][str(tau)]['plaq_canon']}", flush=True)
json.dump(out, open(f"out_2tramp_{Lx}x{Ly}_g{g2f}.json", "w"))
