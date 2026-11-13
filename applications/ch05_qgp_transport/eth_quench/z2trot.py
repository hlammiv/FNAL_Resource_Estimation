"""Trotterized vs exact quench thermalization, Z2 dual on Lx x Ly (full 2^N space; initial states Z2-symmetric).
Second-order step: e^{-i dt Hx/2} e^{-i dt Hzz} e^{-i dt Hx/2}. Tracks <H>, <B>, <E_link> at the Trotter step dt
used for the thermalization stage.  Output: json with time series."""
import numpy as np, scipy.sparse as sp, sys, json
from scipy.sparse.linalg import expm_multiply
Lx, Ly, h = int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3]); J = 1.0
N = Lx * Ly; D = 2 ** N
site = lambda x, y: (x % Lx) + Lx * (y % Ly)
bonds = [(site(x, y), site(x + 1, y)) for y in range(Ly) for x in range(Lx)] + [(site(x, y), site(x, y + 1)) for y in range(Ly) for x in range(Lx)]
s = np.arange(D)
spin = lambda p: 1 - 2 * ((s >> p) & 1)
zzsum = sum(spin(a) * spin(b) for a, b in bonds).astype(float)
Hzz = -h * zzsum


def Xall(v):
    out = np.zeros_like(v)
    for p in range(N):
        out += v[s ^ (1 << p)]
    return out
def Hx(v): return -J * Xall(v)
def H(v): return Hx(v) + Hzz * v
from scipy.sparse.linalg import LinearOperator
Hop = LinearOperator((D, D), matvec=H, dtype=complex)


def expX(v, th):
    # exp(+i th sum_p X_p) applied site by site (Hx = -J sum X  => exp(-i dt Hx) = prod exp(i J dt X))
    c, sn = np.cos(th), 1j * np.sin(th)
    for p in range(N):
        v = c * v + sn * v[s ^ (1 << p)]
    return v
obs = lambda v: (np.vdot(v, H(v)).real / N, np.vdot(v, Xall(v)).real / N, np.vdot(v, zzsum * v).real / len(bonds))
v0s = {"vac": np.zeros(D, complex)}
v0s["vac"][0] = v0s["vac"][-1] = 1 / np.sqrt(2)
rng = np.random.default_rng(1)
tgt = np.abs(zzsum - 0.5 * len(bonds)); cand = np.where(tgt == tgt.min())[0]; cfg = int(rng.choice(cand))
v = np.zeros(D, complex); v[cfg] = v[cfg ^ (D - 1)] = 1 / np.sqrt(2); v0s["rand_cfg"] = v
tmax = float(sys.argv[4]) if len(sys.argv) > 4 else 20.0
out = {"Lx": Lx, "Ly": Ly, "h": h, "tmax": tmax, "runs": {}}
for name, v0 in v0s.items():
    # use small-step 4th-order-accurate reference: Trotter dt=0.02 as "exact"
    for dt in (0.02, 0.1, 0.25, 0.5, 1.0):
        v = v0.copy(); rec = [obs(v)]; ph = np.exp(-1j * dt * Hzz); n_per = int(round(0.5 / dt)) if dt <= 0.5 else 1
        tt = [0.0]; t = 0.0
        while t < tmax - 1e-9:
            v = expX(v, J * dt / 2); v = ph * v; v = expX(v, J * dt / 2); t += dt
            if abs((t / 0.5) - round(t / 0.5)) < 1e-6:
                rec.append(obs(v)); tt.append(t)
        out["runs"][f"{name}_dt{dt}"] = {"t": tt, "E_per_site": [r[0] for r in rec], "B": [r[1] for r in rec], "Elink": [r[2] for r in rec]}
        print(name, dt, rec[-1], flush=True)
json.dump(out, open(f"out_z2trot_{Lx}x{Ly}_h{h}.json", "w"))
