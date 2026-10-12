"""lambda/|M| and p_ref for the Ch. 3 toy 0nubb-like operator (gt_coulomb: sigma.sigma tau+tau+ with a 1/r12
radial form, HO b from hw = 41 A^-1/3) between KB3G ground states (exact M-scheme Lanczos, kept orbits, no
renormalization) of 48Ca and 48Ti.  The operator and its Pauli-LCU 1-norm lambda are the toy's; only the
states change.  Usage: python3 kb3g_lam.py <snt> <space> [k] [tol]."""
import sys, time, json, itertools, math
import numpy as np
sys.dont_write_bytecode = True
import toyop as L
from sm import Model

SPACES = {'f7p3': ['f7', 'p3'], 'f7p3p1': ['f7', 'p3', 'p1'], 'pf': ['f7', 'p3', 'p1', 'f5']}
snt, space = sys.argv[1], sys.argv[2]
k = int(sys.argv[3]) if len(sys.argv) > 3 else 4
tol = float(sys.argv[4]) if len(sys.argv) > 4 else 1e-9
orbs = SPACES[space]
keep = {L.ORB[o] for o in orbs}            # (n, l, 2j) labels
t0 = time.time()

# ---- toy operator terms (unchanged) ----
st = L.states(orbs); ns = len(st); mm = [s[1] for s in st]
Om = L.TBME(orbs, "gt_coulomb")
terms = {}
for a, b in itertools.combinations(range(ns), 2):
    for c, d in itertools.combinations(range(ns), 2):
        if mm[a] + mm[b] != mm[c] + mm[d]:
            continue
        v = Om.me(st[a], st[b], st[c], st[d]) - Om.me(st[a], st[b], st[d], st[c])
        if abs(v) > 1e-10:
            terms[(a, b, c, d)] = v
lam = sum(abs(v) for v in terms.values())
print(f"{space}: {len(terms)} terms, lambda {lam:.4f}, operator built in {time.time()-t0:.1f}s", flush=True)

# ---- KB3G ground states ----
Ca = Model(f'snt/{snt}.snt', keep, 0, 8, 48)
Ti = Model(f'snt/{snt}.snt', keep, 2, 6, 48)
def ground(m, label):
    t = time.time()
    dim, w, v = m.lowest(0, k=k, tol=tol)
    Js = [(-1 + math.sqrt(1 + 4 * m.J2(v[:, i]))) / 2 for i in range(len(w))]
    print(f"  {label}: M-dim {dim}, E {np.round(w, 4)}, Ex {np.round(w - w[0], 4)}, J {np.round(Js, 3)}, {time.time()-t:.1f}s", flush=True)
    assert abs(Js[0]) < 0.05, "lowest state is not J = 0"
    Psi = np.zeros(m.mask.shape); Psi[m.mask] = v[:, 0]
    return Psi, w, Js
PsiCa, wCa, JCa = ground(Ca, '48Ca')
PsiTi, wTi, JTi = ground(Ti, '48Ti')

# ---- mode map toy -> sm ----
def smmode(m, t, kk):
    o, m2 = st[kk]
    lab = L.ORB[o]
    for i, (oi, mi) in enumerate(m.sps[t]):
        if m.orbs[oi][:3] == lab and mi == m2:
            return i
    raise KeyError((o, m2))
p_of = [smmode(Ti, -1, kk) for kk in range(ns)]
n_of = [smmode(Ca, 1, kk) for kk in range(ns)]
assert Ca.sps[1] == Ti.sps[1] and Ca.sps[-1] == Ti.sps[-1]

def apply_ops(src_dets, dst_dets, ops):
    """ops applied right to left on src_dets; returns (src_idx, dst_idx, sign) for survivors found in dst_dets."""
    d = src_dets.copy(); idx = np.arange(len(d)); sgn = np.ones(len(d), dtype=np.int64)
    for mode, create in reversed(ops):
        bit = np.uint64(1 << mode)
        occ = (d & bit) != 0
        ok = ~occ if create else occ
        d, idx, sgn = d[ok], idx[ok], sgn[ok]
        sgn = sgn * (1 - 2 * (np.bitwise_count(d & np.uint64((1 << mode) - 1)).astype(np.int64) & 1))
        d = d ^ bit
    pos = np.searchsorted(dst_dets, d)
    pos = np.minimum(pos, len(dst_dets) - 1)
    found = dst_dets[pos] == d
    return idx[found], pos[found], sgn[found]

# ---- O |Ca> in the Ti det space (protons-first full ordering, as the toy's apply_2b) ----
OPsi = np.zeros(PsiTi.shape)
vac = Ca.dets[-1]; assert len(vac) == 1 and int(vac[0]) == 0
ca_n = PsiCa[0, :]
for (a, b, c, d), v in terms.items():
    si, pi, sp_ = apply_ops(vac, Ti.dets[-1], [(p_of[a], True), (p_of[b], True)])   # a+_a a+_b |0>: b first
    assert len(pi) == 1
    sn_src, sn_dst, sn = apply_ops(Ca.dets[1], Ti.dets[1], [(n_of[d], False), (n_of[c], False)])  # a_d a_c: c first
    np.add.at(OPsi[pi[0]], sn_dst, v * sp_[0] * sn * ca_n[sn_src])
M = float(np.sum(PsiTi * OPsi)); S = float(np.sum(OPsi ** 2))
closed = sum(1 << n_of[kk] for kk in range(ns) if st[kk][0] == 'f7')
p_ref = float(ca_n[np.searchsorted(Ca.dets[1], np.uint64(closed))] ** 2)
out = dict(snt=snt, space=space, n_terms=len(terms), lam=lam, M=M, S=S, lam_over_M=lam / abs(M), p_ref=p_ref,
           gs_share=M * M / S, dim_Ca=int(Ca.mask.sum()), dim_Ti=int(Ti.mask.sum()),
           E_Ca=wCa.tolist(), J_Ca=JCa, E_Ti=wTi.tolist(), J_Ti=JTi, k=k, tol=tol, t_s=time.time() - t0)
print(f"RESULT {snt} {space}: lambda {lam:.4f} M {M:+.5f} lambda/|M| {lam/abs(M):.2f} p_ref {p_ref:.4f} gs share {M*M/S:.4f} "
      f"||O|Ca>||^2 {S:.5f}  ({time.time()-t0:.1f}s)", flush=True)

# ---- cross-check by the toy's own dict-based application (small spaces) ----
if Ca.mask.sum() <= 20000:
    amp = {}
    for j, dn in enumerate(Ca.dets[1]):
        if abs(ca_n[j]) > 1e-14:
            # toy bitstring: protons low bits (none), neutron toy mode kk at bit ns+kk
            s_ = 0
            for kk in range(ns):
                if (int(dn) >> n_of[kk]) & 1:
                    s_ |= 1 << (ns + kk)
            amp[s_] = ca_n[j]
    res = {}
    for s_, x in amp.items():
        for (a, b, c, d), v in terms.items():
            r = L.apply_2b(s_, a, b, c + ns, d + ns)
            if r is None:
                continue
            res[r[1]] = res.get(r[1], 0.0) + x * v * r[0]
    M2 = 0.0
    for s_, x in res.items():
        pbits = sum(1 << p_of[kk] for kk in range(ns) if (s_ >> kk) & 1)
        nbits = sum(1 << n_of[kk] for kk in range(ns) if (s_ >> (ns + kk)) & 1)
        ip = np.searchsorted(Ti.dets[-1], np.uint64(pbits)); jn = np.searchsorted(Ti.dets[1], np.uint64(nbits))
        M2 += x * PsiTi[ip, jn]
    print(f"  cross-check (toy apply_2b route): M {M2:+.5f}  (vectorized {M:+.5f})", flush=True)
    out['M_crosscheck'] = M2
json.dump(out, open(f'lam_{snt}_{space}.json', 'w'), indent=1)
