"""S_int on the chapter's own 15 x 8 geometry: the comparator of scatter.py (same Hamiltonian, wall and packets) on the
open 15-site scattering axis of the 2033 instance, read at the times the instance can reach (t <= t_max = 10/m_phi =
20 a at a m_phi = 0.5). scatter.py scatters on an 800-site chain and reads each bin at 2 (z_wall - z_pk)/v_g, which
the 15-site axis cannot do: the wall sits at z0 = 7 (sites 0-6 incident), the packet starts at z_pk = 2.5 with width
sigma_z = 1.5 sites (2.0 as a check), and the transmitted part reaches the open far end and returns inside t_max.
Here every bin is read on the same time grid, A_R(t) = R_L(t) - R_R(t) with R the incident-side weight, and
S_int(t) is the f_eq-weighted mean over the 7 x 8 incident (k_z, k_y) bins (k_z = 2 pi n/15, n = 1..7; k_z = 0 has
v_g = 0 and is not an incident packet) at T_c (a T = 0.5) and 0.9 T_c (a T = 0.45). The free chain (wall = 0) gives
R0(t), the incident-side weight a non-scattering packet would leave, from which each bin's readout window is read:
first grid time with R0 < 0.2 (cleared), last time before R0 rises again (far-end return).
Usage: python3 scatter15.py ; writes sint15_results.json."""
import sys, json, math
import numpy as np
sys.dont_write_bytecode = True
import scatter as S

N = 15
T_MAX = 20.0                                  # 10/m_phi at a m_phi = 0.5, the priced evolution (100 steps of 0.1/m_phi)
TIMES = list(range(0, 21))                     # lattice units a; t_max = 20 a
TEMPS = {"Tc": 0.5, "0.9Tc": 0.45}


def chain_h(ky, wall, dtheta, width, z0):
    """scatter.chain_h on N = 15 sites (scatter.N is the module-level chain length)."""
    S.N = N
    return S.chain_h(ky, wall, dtheta, width, z0)


def packet(ky, kz, f, zpk, sig):
    u, E, vg = S.bloch(kz, ky)
    z = np.arange(N)
    g = np.exp(-(z - zpk) ** 2 / (4 * sig ** 2)) * np.exp(1j * kz * z)
    psi = np.zeros(N * S.NF * 2, complex)
    for s in range(2):
        psi[(z * S.NF + f) * 2 + s] = g * u[s]
    psi /= np.linalg.norm(psi)
    return psi, E, vg


def left_weight(ev, psi, zcut, times):
    Ev, V = ev
    c = V.conj().T @ psi
    out = []
    for t in times:
        w = np.abs(V @ (np.exp(-1j * Ev * t) * c)) ** 2
        out.append(float(w.reshape(N, S.NF * 2)[:zcut].sum()))
    return out


def readout_window(r0, times, cleared=0.2):
    """(t_clear, t_last): first time the free packet has left the incident side (R0 < cleared) and the last time before
    R0 rises again (the transmitted part returning from the open end). None if the bin never clears inside t_max."""
    t_clear = next((t for t, w in zip(times, r0) if w < cleared), None)
    if t_clear is None:
        return None, None
    i0 = times.index(t_clear)
    i_min = i0 + int(np.argmin(r0[i0:]))
    return t_clear, times[i_min]


def run(wall, dtheta, width, sig=1.5, zpk=2.5, z0=7.0, times=TIMES):
    zcut = int(math.ceil(z0))
    kzs = [2 * math.pi * n / S.LZ for n in range(1, S.LZ // 2 + 1)]
    rows = []
    for ny in range(S.LY):
        ky = 2 * math.pi * ny / S.LY
        ev = np.linalg.eigh(chain_h(ky, wall, dtheta, width, z0))
        ev0 = np.linalg.eigh(chain_h(ky, 0.0, 0.0, width, z0))
        for kz in kzs:
            pL, E, vg = packet(ky, kz, 0, zpk, sig)
            pR, _, _ = packet(ky, kz, 1, zpk, sig)
            RL = left_weight(ev, pL, zcut, times)
            RR = left_weight(ev, pR, zcut, times)
            R0 = left_weight(ev0, pL, zcut, times)
            t_clear, t_last = readout_window(R0, times)
            rows.append(dict(ny=ny, ky=ky, kz=kz, E=E, vg=vg, R_L=RL, R_R=RR, R0=R0,
                             A_R_t=[x - y for x, y in zip(RL, RR)], t_clear=t_clear, t_last=t_last,
                             f_eq={k: 1 / (math.exp(E / aT) + 1) for k, aT in TEMPS.items()}))
    s_int = {}
    for k in TEMPS:
        wsum = sum(r["f_eq"][k] for r in rows)
        s_int[k] = [sum(r["f_eq"][k] * r["A_R_t"][i] for r in rows) / wsum for i in range(len(times))]
    i_max = times.index(T_MAX)
    for r in rows:                                 # the per-bin value the chapter would read at t_max
        r["A_R"] = r["A_R_t"][i_max]
        r["f_eq_Tc"] = r["f_eq"]["Tc"]
    return dict(wall=wall, dtheta=dtheta, width=width, sigma=sig, zpk=zpk, z0=z0, n_sites=N, times=times,
                S_int_t=s_int, S_int=s_int["Tc"][i_max], S_int_09=s_int["0.9Tc"][i_max],
                max_abs_A=max(abs(r["A_R"]) for r in rows), n_modes=len(rows), rows=rows)


if __name__ == "__main__":
    configs = [(0.5, 0.1, 2.0, 1.5), (0.5, 0.1, 2.0, 2.0), (0.5, 0.1, 1.0, 1.5), (0.3, 0.1, 2.0, 1.5),
               (1.0, 0.1, 2.0, 1.5), (0.5, 0.0, 2.0, 1.5)]
    out = []
    for wall, dth, wd, sig in configs:
        r = run(wall, dth, wd, sig)
        out.append(r)
        print(f"wall {wall} dtheta {dth} width {wd} sigma_z {sig}: S_int(t_max) {r['S_int']:+.3e} at T_c, "
              f"{r['S_int_09']:+.3e} at 0.9 T_c; max|A_R| {r['max_abs_A']:.3e}", flush=True)
    main = out[0]
    print("S_int(t) at T_c for the comparator (|y| phi 0.5, width 2 a, sigma_z 1.5):")
    for t, s in zip(main["times"], main["S_int_t"]["Tc"]):
        if t % 2 == 0:
            print(f"  t = {t:2d} a ({t / 2:4.1f}/m_phi): {s:+.3e}")
    print("per-bin readout window at k_y = 0 (t_clear, t_last) and A_R at t_last and t_max:")
    for r in main["rows"]:
        if r["ny"] == 0:
            tl = r["t_last"]
            a_last = r["A_R_t"][main["times"].index(tl)] if tl is not None else float("nan")
            print(f"  k_z = {r['kz']:.3f} v_g = {r['vg']:.2f} f_eq = {r['f_eq_Tc']:.3f}: window ({r['t_clear']}, {tl}) a, "
                  f"A_R {a_last:+.4f} at t_last, {r['A_R']:+.4f} at t_max")
    json.dump(out, open("sint15_results.json", "w"), indent=1)
