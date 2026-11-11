"""Mean-field/linear proxy (the chapter's verlet_proxy construction) on the SAME reduced lattice as exact_step.py:
classical condensate, linear scalar fluctuation at k1 = 2 pi/(3b) with mass^2 1 + 3 lam phibar^2 and lattice k^2,
free Wilson fermion modes h(k) with mass m_psi + g phibar. Same splitting, same readout, no fermion backreaction,
condensate starts at x = 0 (no Yukawa tadpole), as in verlet_proxy."""
import math, json, sys
import numpy as np
import exact_step as E


def proxy(lam, g, h, rec=E.REC_T, x0=0.0):
    V, B, R = E.V, E.B, E.R
    N = int(round(E.TEND / h))
    recn = {int(round(tr / h)): tr for tr in rec}
    xb, pb = x0, math.sqrt(V) * E.PHI0
    k1 = 2 * math.pi / (3 * B)
    kl2 = 4 * math.sin(k1 * B / 2) ** 2 / B ** 2
    W0 = math.sqrt(1 + kl2)
    u, pu = 1 / math.sqrt(2 * W0) + 0j, -1j * math.sqrt(W0 / 2)
    ks = E.ks
    vs = [np.linalg.eigh(E.hk(k, 1.0, E.MPSI))[1][:, 0].astype(complex) for k in ks]
    out = []
    t = 0.0
    for n in range(1, N + 1):
        am = E.a_of(t + h / 2); a3 = am ** 3
        for half in (0, 1):
            pb -= a3 * (xb + lam * xb ** 3 / V) * h / 2
            pu = pu - (a3 * (1 + 3 * lam * xb ** 2 / V) + am * kl2) * u * h / 2
            M = E.MPSI + g * xb / math.sqrt(V) + R / (am * B)
            ph = np.array([np.exp(-1j * M * h / 2), np.exp(1j * M * h / 2)])
            vs = [ph * v for v in vs]
            if half == 0:
                xb += pb / a3 * h
                u = u + pu / a3 * h
                newv = []
                for k, v in zip(ks, vs):
                    hh = (-(R / (am * B)) * math.cos(k * B)) * E.s3 + math.sin(k * B) / (am * B) * E.s1
                    w, U = np.linalg.eigh(hh)
                    newv.append(U @ (np.exp(-1j * w * h) * (U.conj().T @ v)))
                vs = newv
        t = n * h
        if n in recn:
            a = E.a_of(t); a3 = a ** 3
            phibar = xb / math.sqrt(V)
            W = math.sqrt(1 + 3 * lam * phibar ** 2 + kl2 / a ** 2)
            ns = (abs(pu) ** 2 / a3 + a3 * W ** 2 * abs(u) ** 2) / (2 * W) - 0.5
            Mo = E.MPSI + g * phibar
            nf = {}
            for k, v in zip(ks, vs):
                up, _ = E.pos_mode(k, a, Mo)
                nf.setdefault(round(abs(k), 8), []).append(float(abs(np.vdot(up, v)) ** 2))
            out.append(dict(t=recn[n], a=a, phibar=phibar, n_scalar=float(ns),
                            n_ferm=[float(np.mean(v)) for kk, v in sorted(nf.items())]))
    return out


if __name__ == "__main__":
    res = {}
    for lam, g in [(1, 1), (3, 0.5), (1, 2)]:
        for h in [0.05, 0.025, 0.0125, 0.001]:
            res[f"{lam},{g},{h}"] = proxy(lam, g, h)
            print(lam, g, h, res[f"{lam},{g},{h}"][-1], flush=True)
    json.dump(res, open("proxy_red.json", "w"), indent=1)
