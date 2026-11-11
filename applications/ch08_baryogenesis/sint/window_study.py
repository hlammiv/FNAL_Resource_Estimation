"""Ch. 8 readout study (referee follow-up F2, 2026-11-11), on the comparator of scatter15.py (free Wilson fermions,
|y| phi = 0.5 tanh wall of width 2 a, Delta theta_C = 0.1, 15 x 8 cylinder, packets sigma_z = 1.5 at z = 2.5, wall at
z0 = 7, incident side z < 7). Three questions:

1. Window test. What does the 15-site incident-side asymmetry read, and what part of it is reflection? The same wall
   and packets on an axis padded with PAD sites at both ends (nothing returns inside the times read) separate the
   far-end return (15-site minus padded at the same t) from the finite-packet momentum band (padded, late time, minus
   the 800-site sharp-k value of scatter.py).
2. Separate calibrations. The first result reads the f_eq-weighted superposition once (at t_max = 20 a, or at the
   weighted window); the campaign reads each bin at its own t_last (bins with no window at t_max). Each signal is
   reported with its query consequence: first result ~ 9 / (p^2 S^2 m), campaign ~ 1 / (p^2 eps_C^2 m) per bin.
3. Flag efficiency. The model's a = 1/2 + p0 A_R / 2 takes the vacuum charge fluctuation on the incident half as
   independent of the packet (p0 = 0.658 free, 0.646 in the wall vacuum). Here the state is the exact Slater
   determinant: the wall vacuum (all negative-energy orbitals of every k_y sector) plus the packet, projected on the
   positive-energy subspace, evolved to the readout time. Q_inc = sum over sectors of (N_inc - 14). Particle branch
   [Q_inc >= 1], antiparticle branch [Q_inc >= 0]; by C' the flavor-L antiparticle equals the flavor-R particle with
   Q -> -Q, so a_anti = P(Q_inc <= 0 | flavor-R particle). The effective efficiency is
   p_eff = (a_part + a_anti - 1) / A_R(t), with A_R(t) the one-body incident weight difference at the same t.
Usage: python3 window_study.py ; writes window_study.json. About a minute."""
import sys, json, math
import numpy as np
sys.dont_write_bytecode = True
import scatter as S

WALL, DTH, WIDTH, SIG, ZPK, Z0 = 0.5, 0.1, 2.0, 1.5, 2.5, 7.0
N15, PAD = 15, 45
TIMES = list(range(0, 21))
LATE = 60.0                                                    # padded axis: every packet has left the wall by then
TEMPS = {"Tc": 0.5, "0.9Tc": 0.45}
DIM = S.NF * 2                                                 # components per site


def hamiltonian(n, ky, z0, wall=WALL):
    S.N = n
    return S.chain_h(ky, wall, DTH, WIDTH, z0)


def packet(n, ky, kz, f, zpk):
    u, E, vg = S.bloch(kz, ky)
    z = np.arange(n)
    g = np.exp(-(z - zpk) ** 2 / (4 * SIG ** 2)) * np.exp(1j * kz * z)
    psi = np.zeros(n * DIM, complex)
    for s in range(2):
        psi[(z * S.NF + f) * 2 + s] = g * u[s]
    return psi / np.linalg.norm(psi), E, vg


def incident_weight(ev, psi, n, zcut, times):
    Ev, V = ev
    c = V.conj().T @ psi
    return [float((np.abs(V @ (np.exp(-1j * Ev * t) * c)) ** 2).reshape(n, DIM)[:zcut].sum()) for t in times]


def fcs(occ_A):
    """Charge distribution on region A of a Slater determinant whose occupied orbitals restricted to A are the columns
    of occ_A: independent Bernoulli(nu), nu the eigenvalues of occ_A occ_A^dag."""
    nus = np.clip(np.linalg.eigvalsh(occ_A @ occ_A.conj().T), 0, 1)
    d = np.array([1.0])
    for x in nus:
        d = np.convolve(d, [1 - x, x])
    return d                                                    # index = particle number on A


def sector_vacuum(ky, nA):
    e, v = np.linalg.eigh(hamiltonian(N15, ky, Z0))
    return e, v


def run():
    zcut = int(math.ceil(Z0))
    nA = zcut * DIM                                            # 28 modes on the incident half per sector
    kzs = [2 * math.pi * k / S.LZ for k in range(1, S.LZ // 2 + 1)]
    vac = {}                                                   # vacuum charge distribution per sector, offset nA/2
    eig = {}
    for ny in range(S.LY):
        ky = 2 * math.pi * ny / S.LY
        e, v = sector_vacuum(ky, nA)
        eig[ny] = (e, v)
        vac[ny] = fcs(v[:nA, e < 0])
    rows = []
    for ny in range(S.LY):
        ky = 2 * math.pi * ny / S.LY
        e, v = eig[ny]
        neg, pos = v[:, e < 0], v[:, e >= 0]
        ev15 = (e, v)
        evp = np.linalg.eigh(hamiltonian(N15 + 2 * PAD, ky, Z0 + PAD))
        others = np.array([1.0])                               # convolution of the other sectors' vacua
        for m in range(S.LY):
            if m != ny:
                others = np.convolve(others, vac[m])
        off_others = (S.LY - 1) * nA // 2
        for kz in kzs:
            rec = dict(ny=ny, ky=ky, kz=kz)
            for f, tag in ((0, "L"), (1, "R")):
                psi, E, vg = packet(N15, ky, kz, f, ZPK)
                rec.update(E=E, vg=vg)
                rec[f"R15_{tag}"] = incident_weight(ev15, psi, N15, zcut, TIMES)
                pp, _, _ = packet(N15 + 2 * PAD, ky, kz, f, ZPK + PAD)
                rec[f"Rpad_{tag}"] = incident_weight(evp, pp, N15 + 2 * PAD, zcut + PAD, TIMES + [LATE])
                # exact many-body readout: sea + positive-energy part of the packet
                cp = pos @ (pos.conj().T @ psi)
                rec[f"pos_weight_{tag}"] = float(np.vdot(cp, cp).real)
                cp /= np.linalg.norm(cp)
                ct = pos.conj().T @ cp
                probs = []
                for t in TIMES:
                    phi = pos @ (np.exp(-1j * e[e >= 0] * t) * ct)
                    occ = np.column_stack([neg, phi])
                    d_sec = fcs(occ[:nA])                       # N_inc in the packet sector (nA/2 + 1 particles above sea... )
                    tot = np.convolve(d_sec, others)
                    q = np.arange(len(tot)) - (nA // 2 + off_others)
                    probs.append(dict(ge1=float(tot[q >= 1].sum()), le0=float(tot[q <= 0].sum())))
                rec[f"P_{tag}"] = probs
            rec["A15"] = [a - b for a, b in zip(rec["R15_L"], rec["R15_R"])]
            rec["Apad"] = [a - b for a, b in zip(rec["Rpad_L"], rec["Rpad_R"])]
            rec["a_sum"] = [pl["ge1"] + pr["le0"] - 1 for pl, pr in zip(rec["P_L"], rec["P_R"])]
            r0 = incident_weight(np.linalg.eigh(hamiltonian(N15, ky, Z0, wall=0.0)), packet(N15, ky, kz, 0, ZPK)[0], N15, zcut, TIMES)
            t_clear = next((t for t, w in zip(TIMES, r0) if w < 0.2), None)
            rec["t_last"] = None if t_clear is None else TIMES[TIMES.index(t_clear) + int(np.argmin(r0[TIMES.index(t_clear):]))]
            rec["R0"] = r0
            rec["f_eq"] = {k: 1 / (math.exp(rec["E"] / aT) + 1) for k, aT in TEMPS.items()}
            rows.append(rec)
    return rows


def summarize(rows, vac_p0=0.6464):
    it = TIMES.index(20)
    out = {}
    for k in TEMPS:
        w = np.array([r["f_eq"][k] for r in rows]); w /= w.sum()
        S15 = [float(w @ np.array([r["A15"][i] for r in rows])) for i in range(len(TIMES))]
        Spad = [float(w @ np.array([r["Apad"][i] for r in rows])) for i in range(len(TIMES) + 1)]
        Wsum = [float(w @ np.array([r["a_sum"][i] for r in rows])) for i in range(len(TIMES))]
        camp_t = [r["t_last"] if r["t_last"] is not None else 20 for r in rows]
        S_camp = float(sum(wi * r["A15"][TIMES.index(t)] for wi, r, t in zip(w, rows, camp_t)))
        W_camp = float(sum(wi * r["a_sum"][TIMES.index(t)] for wi, r, t in zip(w, rows, camp_t)))
        R0w = [float(w @ np.array([r["R0"][i] for r in rows])) for i in range(len(TIMES))]
        out[k] = dict(S15_t=S15, Spad_t=Spad[:-1], Spad_late=Spad[-1], S15_tmax=S15[it], Spad_tmax=Spad[it],
                      return_part_tmax=S15[it] - Spad[it], S_campaign_windowed=S_camp,
                      flag_signal_tmax=Wsum[it], p_eff_first=Wsum[it] / S15[it] if S15[it] else None,
                      flag_signal_campaign=W_camp, p_eff_campaign=W_camp / S_camp if S_camp else None,
                      weighted_R0=R0w)
        # per-bin efficiency where the bin's signal is not tiny
        pe = [r["a_sum"][TIMES.index(t)] / r["A15"][TIMES.index(t)] for r, t in zip(rows, camp_t) if abs(r["A15"][TIMES.index(t)]) > 5e-3]
        out[k]["p_eff_bins_range"] = (float(min(pe)), float(max(pe))) if pe else None
    out["pos_weight_min"] = min(min(r["pos_weight_L"], r["pos_weight_R"]) for r in rows)
    return out


if __name__ == "__main__":
    rows = run()
    s = summarize(rows)
    for k in TEMPS:
        o = s[k]
        print(f"[{k}] S15(t_max) {o['S15_tmax']:+.3e}  padded at t_max {o['Spad_tmax']:+.3e}  return part {o['return_part_tmax']:+.3e}  "
              f"padded late {o['Spad_late']:+.3e}  campaign windowed {o['S_campaign_windowed']:+.3e}")
        print(f"     flag signal (a_part + a_anti - 1) at t_max {o['flag_signal_tmax']:+.3e} -> p_eff first {o['p_eff_first']:.3f}; "
              f"campaign {o['flag_signal_campaign']:+.3e} -> p_eff {o['p_eff_campaign']:.3f}; per-bin p_eff {o['p_eff_bins_range']}")
    print("weighted free R0 (incident weight of a non-scattering packet) vs t:", [round(x, 3) for x in s["Tc"]["weighted_R0"]])
    print("S15(t) at T_c:", [f"{x:+.2e}" for x in s["Tc"]["S15_t"][8:]])
    print("Spad(t) at T_c:", [f"{x:+.2e}" for x in s["Tc"]["Spad_t"][8:]])
    print("min positive-energy weight of a packet:", round(s["pos_weight_min"], 4))
    json.dump(dict(summary=s, rows=rows), open("window_study.json", "w"))
