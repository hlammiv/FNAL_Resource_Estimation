"""27Al in sd (USDB): <S_p>, <O_SD>, Var(O_SD) in the J = M = 5/2 ground state, and the overlap of
single-determinant references with it (projection success p). Ch. 3 open items (author ruling 2026-10-09, item 8).

O_SD = sum over valence protons of o_l(q) sigma_z, with o_l the j0(qr) radial matrix element of the orbit
(o_0d for 0d5/2, 0d3/2; o_1s for 1s1/2), from estimates.ch03_mu2e_0nubb.al27_dipole_sd (the chapter's b and q).
sigma_z connects j = l +- 1/2 at the same l, m, so <O_SD> is not simply 2 <S_p> o_0d; both are printed.
"""
import sys, json, math, time
import numpy as np
import scipy.sparse as sp
HERE = str(__import__('pathlib').Path(__file__).resolve().parent)
sys.path.insert(0, HERE)
sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[3]))   # repository root
from sm import Model, cg
from estimates import ch03_mu2e_0nubb as ch3

d = ch3.al27_dipole_sd(ch3.Assumptions())
o_l = {2: d['o_0d'], 0: d['o_1s']}
print("o_0d", d['o_0d'], "o_1s", d['o_1s'], "lam_lcu", d['lam_lcu'], flush=True)

sd = {(0, 2, 5), (1, 0, 1), (0, 2, 3)}
t0 = time.time()
m = Model(HERE + '/usdb.snt', sd, 5, 6, 27)
dim, w, v = m.lowest(5, k=3)
print(f"M=5/2 dim {dim}, E {w}, {time.time()-t0:.0f}s", flush=True)
vec = v[:, 0]
Psi = np.zeros(m.mask.shape); Psi[m.mask] = vec
J = (-1 + math.sqrt(1 + 4 * m.J2(vec))) / 2
print("ground state J =", J, flush=True)


def sigz_me(l, j2a, j2b, m2):
    """<l j_a m | sigma_z | l j_b m>, all j, m doubled."""
    tot = 0.0
    for ms2 in (-1, 1):
        ml2 = m2 - ms2
        if abs(ml2) > 2 * l:
            continue
        tot += cg(2 * l, ml2, 1, ms2, j2a, m2) * cg(2 * l, ml2, 1, ms2, j2b, m2) * ms2
    return tot


Np = len(m.dets[-1])
op_spin = sp.csr_matrix((Np, Np)); op_sd = sp.csr_matrix((Np, Np))
for (oa, ma), ia in m.idx[-1].items():
    na, la, j2a, _ = m.orbs[oa]
    for (ob, mb), ib in m.idx[-1].items():
        nb, lb, j2b, _ = m.orbs[ob]
        if (na, la, ma) != (nb, lb, mb):
            continue
        c = sigz_me(la, j2a, j2b, ma)
        if abs(c) < 1e-14:
            continue
        ob_ = m.one_body(-1, ia, ib)
        op_spin = op_spin + c * ob_
        op_sd = op_sd + c * o_l[la] * ob_
op_spin = op_spin.tocsr(); op_sd = op_sd.tocsr()

res = {"dim": dim, "E0": float(w[0]), "J": J, "o_0d": d['o_0d'], "o_1s": d['o_1s'], "lam_lcu": d['lam_lcu']}
for name, op in (("sigma_z_p", op_spin), ("O_SD", op_sd)):
    OPsi = op @ Psi
    ev = float(np.sum(Psi * OPsi)); ev2 = float(np.sum(OPsi * OPsi))
    res[name] = {"expectation": ev, "second_moment": ev2, "variance": ev2 - ev * ev}
    print(name, res[name], flush=True)
res["S_p"] = res["sigma_z_p"]["expectation"] / 2
res["O_SD_over_2Sp_o0d"] = res["O_SD"]["expectation"] / (2 * res["S_p"] * d['o_0d'])
res["hadamard_factor"] = (d['lam_lcu'] / res["O_SD"]["expectation"]) ** 2
res["shots_p1_eps005"] = res["O_SD"]["variance"] / (0.05 * res["O_SD"]["expectation"]) ** 2

# references: naive 0d5/2 filling (protons m = 5/2..-3/2, neutrons all six), and the largest component
od5 = {t: [i for i, o in m.orbs.items() if o[:3] == (0, 2, 5) and o[3] == t][0] for t in (-1, 1)}
pdet = sum(1 << m.idx[-1][(od5[-1], m2)] for m2 in (5, 3, 1, -1, -3))
ndet = sum(1 << m.idx[1][(od5[1], m2)] for m2 in (5, 3, 1, -1, -3, -5))
ip = int(np.searchsorted(m.dets[-1], np.uint64(pdet))); in_ = int(np.searchsorted(m.dets[1], np.uint64(ndet)))
assert m.dets[-1][ip] == pdet and m.dets[1][in_] == ndet
res["p_naive_0d52"] = float(Psi[ip, in_] ** 2)
flat = np.argmax(Psi ** 2); ipm, inm = np.unravel_index(flat, Psi.shape)
res["p_largest_component"] = float(Psi[ipm, inm] ** 2)
res["largest_is_naive"] = bool(ipm == ip and inm == in_)
# occupation of the largest-component determinant, for the record
occ = lambda t, det: [(m.orbs[o][:3], m2) for (o, m2), i in m.idx[t].items() if (int(det) >> i) & 1]
res["largest_component_occ"] = {"p": occ(-1, m.dets[-1][ipm]), "n": occ(1, m.dets[1][inm])}
# cumulative weight of the k largest components (how much a few-determinant reference could reach)
srt = np.sort((Psi ** 2).ravel())[::-1]
res["cum_weight_top"] = {k: float(srt[:k].sum()) for k in (1, 2, 5, 10, 100)}
print(json.dumps(res, indent=1, default=str), flush=True)
json.dump(res, open(HERE + '/al27_sd_readout.json', 'w'), indent=1, default=str)
