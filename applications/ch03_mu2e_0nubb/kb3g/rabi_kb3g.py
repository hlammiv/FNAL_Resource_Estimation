"""The Ch. 3 resonant-drive (Rabi) readout toy (simplobs adv/rabi.py), rerun with KB3G ground states and
spectrum in a truncated pf space.  Same protocol: evolve the closed-f7/2 reference under H_p - dE + det,
H_n, theta (O + O^dag) to t_pi/2 = pi/(2 theta |M|), scan the detuning for the resonance, then fit
P_daughter(t) = A sin^2(W t) over 0..2 t_pi/2 and compare W with theta |M|."""
import sys, time, math, itertools, numpy as np
sys.dont_write_bytecode = True
import toyop as L
from sm import Model
from scipy.sparse import coo_matrix, bmat, identity, csr_matrix
from scipy.sparse.linalg import expm_multiply, eigsh
SPACES = {'f7p3': ['f7', 'p3'], 'f7p3p1': ['f7', 'p3', 'p1']}
snt, space = sys.argv[1], sys.argv[2]
thMs = tuple(float(x) for x in (sys.argv[3] if len(sys.argv) > 3 else "0.3,0.5").split(","))
orbs = SPACES[space]; keep = {L.ORB[o] for o in orbs}
st = L.states(orbs); ns = len(st); mm = [s[1] for s in st]
Om = L.TBME(orbs, "gt_coulomb"); terms = {}
for a, b in itertools.combinations(range(ns), 2):
    for c, d in itertools.combinations(range(ns), 2):
        if mm[a] + mm[b] == mm[c] + mm[d]:
            v = Om.me(st[a], st[b], st[c], st[d]) - Om.me(st[a], st[b], st[d], st[c])
            if abs(v) > 1e-10: terms[(a, b, c, d)] = v
Ca = Model(f'snt/{snt}.snt', keep, 0, 8, 48); Ti = Model(f'snt/{snt}.snt', keep, 2, 6, 48)
Ca.setup_sector(0); Ti.setup_sector(0)
def smmode(m, t, kk):
    lab = L.ORB[st[kk][0]]
    return next(i for i, (oi, mi) in enumerate(m.sps[t]) if m.orbs[oi][:3] == lab and mi == st[kk][1])
p_of = [smmode(Ti, -1, kk) for kk in range(ns)]; n_of = [smmode(Ca, 1, kk) for kk in range(ns)]
def apply_ops(src, dst, ops):
    d = src.copy(); idx = np.arange(len(d)); sgn = np.ones(len(d), dtype=np.int64)
    for mode, create in reversed(ops):
        bit = np.uint64(1 << mode); occ = (d & bit) != 0; ok = ~occ if create else occ
        d, idx, sgn = d[ok], idx[ok], sgn[ok]
        sgn = sgn * (1 - 2 * (np.bitwise_count(d & np.uint64((1 << mode) - 1)).astype(np.int64) & 1)); d = d ^ bit
    pos = np.minimum(np.searchsorted(dst, d), len(dst) - 1); f = dst[pos] == d
    return idx[f], pos[f], sgn[f]
def dense_H(m):
    """H in the M = 0 sector: dense below 6000, else sparse (columns from matvec)."""
    dim = int(m.mask.sum())
    if dim <= 6000:
        H = np.zeros((dim, dim))
        for j in range(dim):
            P = np.zeros(m.mask.shape); x = np.zeros(dim); x[j] = 1; P[m.mask] = x
            H[:, j] = m.matvec_full(P)[m.mask]
        return H
    R, C, V = [], [], []
    for j in range(dim):
        P = np.zeros(m.mask.shape); x = np.zeros(dim); x[j] = 1; P[m.mask] = x
        col = m.matvec_full(P)[m.mask]; nz = np.flatnonzero(np.abs(col) > 1e-13)
        R.append(nz); C.append(np.full(len(nz), j)); V.append(col[nz])
    return csr_matrix((np.concatenate(V), (np.concatenate(R), np.concatenate(C))), shape=(dim, dim))
t0 = time.time()
Hi = dense_H(Ca); Hf = dense_H(Ti)
ni, nf = Hi.shape[0], Hf.shape[0]
print(f"{snt} {space}: dims {ni} {nf}, H built {time.time()-t0:.1f}s", flush=True)
# O as nf x ni in mask-flat indices
flatCa = -np.ones(Ca.mask.shape, dtype=np.int64); flatCa[Ca.mask] = np.arange(ni)
flatTi = -np.ones(Ti.mask.shape, dtype=np.int64); flatTi[Ti.mask] = np.arange(nf)
R, C, V = [], [], []
for (a, b, c, d), v in terms.items():
    _, pi, sp_ = apply_ops(Ca.dets[-1], Ti.dets[-1], [(p_of[a], True), (p_of[b], True)])
    s_, dd, sn = apply_ops(Ca.dets[1], Ti.dets[1], [(n_of[d], False), (n_of[c], False)])
    src = flatCa[0, s_]; dst = flatTi[pi[0], dd]; ok = (src >= 0) & (dst >= 0)
    R.append(dst[ok]); C.append(src[ok]); V.append(v * sp_[0] * sn[ok])
O = coo_matrix((np.concatenate(V), (np.concatenate(R), np.concatenate(C))), shape=(nf, ni)).tocsr()
def low(H, k=6):
    if isinstance(H, np.ndarray):
        return np.linalg.eigh(H)
    w, v = eigsh(H, k=k, which='SA', tol=1e-10); o = np.argsort(w); return w[o], v[:, o]
wi, vi = low(Hi); wf, vf = low(Hf)
M = float(vf[:, 0] @ (O @ vi[:, 0]))
closed = sum(1 << n_of[kk] for kk in range(ns) if st[kk][0] == 'f7')
ref = np.zeros(ni); ref[flatCa[0, np.searchsorted(Ca.dets[1], np.uint64(closed))]] = 1.0
pref = float(vi[:, 0] @ ref) ** 2
Sfull = float(np.linalg.norm(O @ vi[:, 0]) ** 2)
print(f"  M {M:+.5f} lambda/|M| {sum(abs(v) for v in terms.values())/abs(M):.1f} p_ref {pref:.4f} gs share {M*M/Sfull:.4f} "
      f"Ei {wi[:3].round(3)} Ef {wf[:4].round(3)} gap_f(0+) {wf[np.argmax(np.array([abs(vf[:,i] @ np.ones(nf)) for i in range(1)]))] if False else ''}", flush=True)
dE = wi[0] - wf[0]
Hp = csr_matrix(Hi - dE * np.eye(ni)) if isinstance(Hi, np.ndarray) else (Hi - dE * identity(ni)).tocsr(); Hfs = csr_matrix(Hf); Os = O
for thM in thMs:
    th = thM / abs(M)
    Htot = bmat([[Hp, th * Os.T], [th * Os, Hfs]]).tocsr()
    tq = math.pi / 2 / (th * abs(M))
    psi0 = np.concatenate([ref, np.zeros(nf)]).astype(complex)
    best = None
    for det in np.linspace(-1.0, 1.0, 41) * thM / 0.1:
        Hd = Htot + bmat([[det * identity(ni), None], [None, csr_matrix((nf, nf))]]).tocsr()
        psi = expm_multiply(-1j * Hd, psi0, start=0, stop=tq, num=2, endpoint=True)[-1]
        pd = float(np.linalg.norm(psi[ni:]) ** 2); pgs = float(abs(vf[:, 0] @ psi[ni:]) ** 2)
        if best is None or pd > best[1]: best = (det, pd, pgs)
    det = best[0]
    Hd = Htot + bmat([[det * identity(ni), None], [None, csr_matrix((nf, nf))]]).tocsr()
    ts = np.linspace(0, 2 * tq, 41)
    traj = expm_multiply(-1j * Hd, psi0, start=0, stop=2 * tq, num=41, endpoint=True)
    pd = np.array([np.linalg.norm(p[ni:]) ** 2 for p in traj])
    best2 = None
    for W in np.linspace(0.3, 2.0, 1701) * th * abs(M):
        s2 = np.sin(W * ts) ** 2; A = (s2 @ pd) / (s2 @ s2); r = ((pd - A * s2) ** 2).sum()
        if best2 is None or r < best2[0]: best2 = (r, W, A)
    print(f"  theta|M| = {thM:.1f} MeV: resonance offset {det:+.3f} MeV, P_daughter(t_pi/2) {best[1]:.3f} (gs {best[2]:.3f});"
          f" Rabi fit |M|_fit/|M| = {best2[1]/(th*abs(M)):.3f}, amplitude {best2[2]:.3f}; t_pi/2 = {tq:.2f} MeV^-1  ({time.time()-t0:.0f}s)", flush=True)
