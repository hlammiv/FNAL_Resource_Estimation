import sys, math, numpy as np, json
sys.dont_write_bytecode = True
import scatter as S
# 1. convergence: N = 1000, sigma = 40, measurement time x1.5 and x2
S.N = 1000
wall, dth, wd = 0.5, 0.1, 2.0
rows = {}
for ny in (0, 1, 4):
    ky = 2 * math.pi * ny / S.LY
    H = S.chain_h(ky, wall, dth, wd, 500); ev = np.linalg.eigh(H)
    for kz in [2 * math.pi * n / S.LZ for n in (1, 3, 5, 7)]:
        for sig in (25.0, 40.0):
            RL, E, vg = S.reflect(ev, None, ky, kz, 0, 250, sig, 500)
            RR, _, _ = S.reflect(ev, None, ky, kz, 1, 250, sig, 500)
            rows[(ny, round(kz, 3), sig)] = (RL, RR, RL - RR)
    # hole (antiparticle) of flavour L: negative-energy packet; C' predicts R(pbar_L) = R(p_R)
    for kz in [2 * math.pi * n / S.LZ for n in (3, 5)]:
        hk = lambda q: S.s3 * (S.AM + 2 * S.RW - S.RW * math.cos(q) - S.RW * math.cos(ky)) + S.s1 * math.sin(q) + S.s2 * math.sin(ky)
        e, v = np.linalg.eigh(hk(kz)); u = v[:, 0]
        dq = 1e-4; vg = (np.linalg.eigvalsh(hk(kz + dq))[0] - np.linalg.eigvalsh(hk(kz - dq))[0]) / (2 * dq)
        # negative-energy packet moving toward +z needs vg > 0: use -kz if needed
        if vg < 0:
            e, v = np.linalg.eigh(hk(-kz)); u = v[:, 0]; kk = -kz; vg = -vg
        else:
            kk = kz
        z = np.arange(S.N); g = np.exp(-(z - 250) ** 2 / (4 * 25.0 ** 2)) * np.exp(1j * kk * z)
        psi = np.zeros(S.N * 4, complex)
        for s in range(2): psi[(z * 2 + 0) * 2 + s] = g * u[s]
        psi /= np.linalg.norm(psi); t = 2.0 * 250 / vg
        Ev, V = ev; c = V.conj().T @ psi; pt = V @ (np.exp(-1j * Ev * t) * c)
        Rhole = float((np.abs(pt) ** 2).reshape(S.N, 4)[:500].sum())
        print(f"hole check ny={ny} kz={kz:.3f}: R(hole_L) {Rhole:.4f} vs R_R {rows[(ny, round(kz,3), 25.0)][1]:.4f} vs R_L {rows[(ny, round(kz,3), 25.0)][0]:.4f}")
for k, v in sorted(rows.items()):
    print(k, [round(x, 5) for x in v])
