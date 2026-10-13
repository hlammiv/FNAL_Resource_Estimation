"""S_int for Ch. 8: free-fermion reflection asymmetry of the chapter's two-flavour 2-component Wilson fermions
(r = 1, a m_f = 0.5) off the C-violating tanh wall |y| phi(z) e^{i theta_C(z)} psi_L^dag beta psi_R + h.c.,
prof(z) = (1 + tanh((z - z0)/width))/2, theta_C = dtheta x prof, on a long open chain per k_y sector of the
L_y = 8 cylinder. R(p) = weight of a Gaussian packet of the symmetric-phase positive-energy eigenstate
(momentum +k_z, flavour f) left on the incident side after scattering. A_R(k) = R_L(k) - R_R(k) (C' maps the
flavour-L antiparticle onto the flavour-R particle). S_int = sum_k f_eq(k) A_R(k) / sum_k f_eq(k) over the
incident modes (k_z > 0), f_eq Fermi-Dirac at T = 1/beta, beta m_phi = 1, a m_phi = 0.5 -> a T = 0.5.
Usage: python3 scatter.py [N] ; writes sint_results.json."""
import sys, json, math, time
import numpy as np
N = int(sys.argv[1]) if len(sys.argv) > 1 else 800
AM, RW, NF, LY, LZ = 0.5, 1.0, 2, 8, 15
aT = 0.5
s1 = np.array([[0, 1], [1, 0]], complex); s2 = np.array([[0, -1j], [1j, 0]]); s3 = np.diag([1.0, -1.0]).astype(complex)
hop = lambda alpha: -0.5 * (RW * s3 - 1j * alpha)

def chain_h(ky, wall, dtheta, width, z0):
    n = N * NF * 2
    h = np.zeros((n, n), complex)
    onsite = s3 * (AM + 2 * RW) + hop(s2) * np.exp(1j * ky) + hop(s2).conj().T * np.exp(-1j * ky)
    prof = lambda z: 0.5 * (1 + math.tanh((z - z0) / width))
    for z in range(N):
        for f in range(NF):
            i = (z * NF + f) * 2
            h[i:i + 2, i:i + 2] += onsite
            if z + 1 < N:
                j = ((z + 1) * NF + f) * 2
                h[i:i + 2, j:j + 2] += hop(s1); h[j:j + 2, i:i + 2] += hop(s1).conj().T
        if wall:
            i, j = (z * NF) * 2, (z * NF + 1) * 2
            blk = wall * prof(z) * np.exp(1j * dtheta * prof(z)) * s3
            h[i:i + 2, j:j + 2] += blk; h[j:j + 2, i:i + 2] += blk.conj().T
    return h

def bloch(kz, ky):
    """symmetric-phase 2x2 Bloch Hamiltonian of chain_h (hop(alpha) e^{iq} + h.c. = -s3 cos q - alpha sin q), its
    positive-energy spinor, energy and z group velocity. The sign of the sin terms was wrong until 2026-10-12 (verification
    pass), which launched s3 u instead of u and mixed in a negative-energy component; S_int was then 1.33e-2."""
    hk = lambda q: s3 * (AM + 2 * RW - RW * math.cos(q) - RW * math.cos(ky)) - s1 * math.sin(q) - s2 * math.sin(ky)
    e, v = np.linalg.eigh(hk(kz)); u = v[:, 1]; E = e[1]
    dq = 1e-4; vg = (np.linalg.eigvalsh(hk(kz + dq))[1] - np.linalg.eigvalsh(hk(kz - dq))[1]) / (2 * dq)
    return u, E, vg

def reflect(ev, vec, ky, kz, f, zpk, sig, zwall):
    """R for a packet of flavour f at (kz, ky); ev = (E, V) eigendecomposition of the chain Hamiltonian."""
    u, E, vg = bloch(kz, ky)
    if vg <= 0.02:
        return None, E, vg
    z = np.arange(N)
    g = np.exp(-(z - zpk) ** 2 / (4 * sig ** 2)) * np.exp(1j * kz * z)
    psi = np.zeros(N * NF * 2, complex)
    for s in range(2):
        psi[(z * NF + f) * 2 + s] = g * u[s]
    psi /= np.linalg.norm(psi)
    t = 2.0 * (zwall - zpk) / vg                       # incident packet would be 2x the wall distance past its start
    Ev, V = ev
    c = V.conj().T @ psi
    psi_t = V @ (np.exp(-1j * Ev * t) * c)
    w = np.abs(psi_t) ** 2
    left = w.reshape(N, NF * 2)[:zwall].sum()
    return float(left), float(E), float(vg)

def run(wall, dtheta, width, zwall=None, zpk=None, sig=25.0):
    zwall = zwall or N // 2; zpk = zpk or N // 4
    rows = []
    kzs = [2 * math.pi * n / LZ for n in range(1, LZ // 2 + 1)]       # incident |k_z| bins on the 15-site axis
    for ny in range(LY):
        ky = 2 * math.pi * ny / LY
        H = chain_h(ky, wall, dtheta, width, zwall)
        ev = np.linalg.eigh(H)
        for kz in kzs:
            RL, E, vg = reflect(ev, None, ky, kz, 0, zpk, sig, zwall)
            RR, _, _ = reflect(ev, None, ky, kz, 1, zpk, sig, zwall)
            rows.append(dict(ny=ny, ky=ky, kz=kz, E=E, vg=vg, R_L=RL, R_R=RR,
                             A_R=(None if RL is None else RL - RR), f_eq=1 / (math.exp(E / aT) + 1)))
    ok = [r for r in rows if r["A_R"] is not None]
    wsum = sum(r["f_eq"] for r in ok)
    S = sum(r["f_eq"] * r["A_R"] for r in ok) / wsum
    Sabs = sum(r["f_eq"] * abs(r["A_R"]) for r in ok) / wsum
    return dict(wall=wall, dtheta=dtheta, width=width, S_int=S, S_int_abs=Sabs, n_modes=len(ok), n_slow=len(rows) - len(ok),
                max_abs_A=max(abs(r["A_R"]) for r in ok), mean_R=float(np.mean([r["R_L"] for r in ok])), rows=rows)

if __name__ == "__main__":
    t0 = time.time(); out = []
    configs = [(0.5, 0.1, 1.0), (0.5, 0.1, 2.0), (0.3, 0.1, 2.0), (1.0, 0.1, 2.0), (0.5, 0.2, 2.0), (0.5, 0.0, 2.0)]
    for wall, dth, wd in configs:
        r = run(wall, dth, wd)
        out.append(r)
        print(f"wall {wall} dtheta {dth} width {wd}: S_int {r['S_int']:+.3e} (|A| weighted {r['S_int_abs']:.3e}), max|A_R| {r['max_abs_A']:.3e}, "
              f"mean R_L {r['mean_R']:.3f}, modes {r['n_modes']} (+{r['n_slow']} too slow)  [{time.time()-t0:.0f}s]", flush=True)
    json.dump(out, open("sint_results.json", "w"), indent=1)
