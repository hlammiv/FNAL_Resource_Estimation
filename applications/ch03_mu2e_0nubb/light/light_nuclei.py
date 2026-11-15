"""Light-nucleus alternatives to the 2028 4He benchmark (Roni P1-01 follow-up).

Same space and Hamiltonian as gaps/he4D.py: the 40-mode 0s,0p,1s0d JW register, b = 1.369 fm (the 4He fit,
kept fixed for every nucleus), Minnesota NN (u = 1, no Coulomb), intrinsic kinetic energy for mass A, Lawson
c.m. term beta (H_cm - 3/2 hw). The Minnesota force is central (Wigner/Majorana/Bartlett/Heisenberg mixture),
so total L and S commute with H; we diagonalize in the stretched sector M_L = L, M_S = S, whose lowest state
is |J = L + S, M_J = J>, after checking that the lower sectors give the same energy.

For each ground state: E0, L, S, D (largest-|c| determinants for fidelity 0.99 / 0.999), and for one-body
proxies of the six Mu2e responses (Haxton et al. 2208.07945), proton and neutron separately:
  M0     j0(qr)                         (charge, J = 0)
  M2     j2(qr) Y20                     (charge, J = 2)
  S1L0   j0(qr) sigma_z                 (L = 0 part of Sigma'_1 and Sigma''_1)
  S1L2   j2(qr) [Y2 x sigma]^1_0        (L = 2 part; separates Sigma'_1 from Sigma''_1)
  D1     j0(qr) l_z                     (Delta_1, long-wavelength proxy)
  P0     j0(qr) l.s                     (Phi''_0, long-wavelength proxy)
  P2     j0(qr) [l x s]^2_0             (Phi~'_2 / Phi''_2 class, long-wavelength proxy)
<O> in |J, M_J = J>, the JW one-norm lambda' after identity subtraction, Gamma, and (lambda'/|<O>|)^2.
Usage: python light_nuclei.py NAME Z N L 2S [beta] [nmax_extra]
"""
import sys, math, itertools, json, time
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh
from numpy.polynomial.hermite import hermgauss
from scipy.special import eval_hermite, spherical_jn
try:
    from scipy.special import sph_harm_y
except ImportError:                      # scipy < 1.15
    from scipy.special import sph_harm
    sph_harm_y = lambda l, m, th, ph: sph_harm(m, l, ph, th)
from scipy import integrate
from sympy.physics.wigner import gaunt, clebsch_gordan
from sympy import Rational
sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[3]))
from estimates import ch03_mu2e_0nubb as c3

NAME = sys.argv[1]; Z = int(sys.argv[2]); N = int(sys.argv[3]); LGS = int(sys.argv[4]); S2 = int(sys.argv[5])
BETA = float(sys.argv[6]) if len(sys.argv) > 6 else 10.0
NEXTRA = int(sys.argv[7]) if len(sys.argv) > 7 else 99          # HO quanta above the minimum kept
A = Z + N
a = c3.Assumptions()
osc = c3.he4_oscillator_length(a)
b, hw, q = osc['b'], osc['hbar_omega'], osc['q']
MINN = [(200.0, 1.487, 'R'), (-178.0, 0.639, 't'), (-91.85, 0.465, 's')]
U_EXCH = 1.0
t0 = time.time()

# ---------------- one- and two-body elements (as gaps/he4D.py, mass A) ----------------
cart = [(nx, ny, nz) for Nq in range(3) for nx in range(Nq + 1) for ny in range(Nq + 1 - nx) for nz in [Nq - nx - ny]]
def phi1(n, u):
    return eval_hermite(n, u) * np.exp(-u * u / 2) / math.sqrt(2 ** n * math.factorial(n) * math.sqrt(math.pi))
xg, wg = hermgauss(80)
U1, U2 = np.meshgrid(xg, xg, indexing='ij'); W = np.outer(wg, wg)
P = np.array([phi1(n, xg) * np.exp(xg * xg / 2) for n in range(3)])
def G1(kap):
    g = np.exp(-kap * b * b * (U1 - U2) ** 2) * W
    return np.einsum('ai,bj,ci,dj,ij->abcd', P, P, P, P, g)
def V3(kap):
    g = G1(kap); n = len(cart); out = np.empty((n, n, n, n))
    for i, ci in enumerate(cart):
        for j, cj in enumerate(cart):
            for k, ck in enumerate(cart):
                for l, cl in enumerate(cart):
                    out[i, j, k, l] = np.prod([g[ci[d], cj[d], ck[d], cl[d]] for d in range(3)])
    return out
pt = np.zeros((3, 3), complex)
for n in range(3):
    if n + 1 < 3: pt[n + 1, n] = 1j * math.sqrt((n + 1) / 2)
    if n - 1 >= 0: pt[n - 1, n] = -1j * math.sqrt(n / 2)
def p2half(m, n):
    if m == n: return (2 * n + 1) / 4
    if abs(m - n) == 2: return -math.sqrt(max(m, n) * (max(m, n) - 1)) / 4
    return 0.0
nc = len(cart)
T1 = np.zeros((nc, nc)); Pd = []; Rd = []
for i, ci in enumerate(cart):
    for j, cj in enumerate(cart):
        for d in range(3):
            if all(ci[e] == cj[e] for e in range(3) if e != d):
                T1[i, j] += p2half(ci[d], cj[d])
for d in range(3):
    m = np.zeros((nc, nc), complex); r_ = np.zeros((nc, nc))
    for i, ci in enumerate(cart):
        for j, cj in enumerate(cart):
            if all(ci[e] == cj[e] for e in range(3) if e != d):
                m[i, j] = pt[ci[d], cj[d]]
                if ci[d] == cj[d] + 1: r_[i, j] = math.sqrt(ci[d] / 2)
                if ci[d] == cj[d] - 1: r_[i, j] = math.sqrt(cj[d] / 2)
    Pd.append(m); Rd.append(r_)
xs, ws = hermgauss(12)
X, Y, Zg = np.meshgrid(xs, xs, xs, indexing='ij')
Wt = (ws[:, None, None] * ws[None, :, None] * ws[None, None, :]).ravel()
X, Y, Zg = X.ravel(), Y.ravel(), Zg.ravel()
R = np.sqrt(X * X + Y * Y + Zg * Zg); TH = np.arccos(Zg / R); PH = np.arctan2(Y, X)
cart_vals = np.array([phi1(c[0], X) * phi1(c[1], Y) * phi1(c[2], Zg) for c in cart]) * np.exp(R * R / 2)
sph_vals = np.array([c3.ho_radial(n, l, R, 1.0) * sph_harm_y(l, m, TH, PH) * np.exp(R * R / 2)
                     for (n, l, m) in c3.HE4_SPATIAL])
Umat = np.einsum('sk,ck,k->sc', sph_vals, cart_vals, Wt)
assert np.allclose(Umat @ Umat.conj().T, np.eye(10), atol=1e-10)
to1 = lambda h: Umat.conj() @ h @ Umat.T
to2 = lambda v: np.einsum('ia,jb,kc,ld,abcd->ijkl', Umat.conj(), Umat.conj(), Umat, Umat, v)
Ncart = np.diag([sum(c) for c in cart]).astype(float)
t_sph = to1(T1 * hw * (A - 1) / A + BETA * hw / A * (Ncart + 1.5 * np.eye(nc)))
K2 = sum(np.einsum('ac,bd->abcd', Pd[d], Pd[d]) for d in range(3)) * (-hw / A)
K2 = K2 + BETA * hw / A * sum(np.einsum('ac,bd->abcd', Pd[d], Pd[d]) + np.einsum('ac,bd->abcd', Rd[d], Rd[d]) for d in range(3))
K2s = to2(K2)
Vg = {lab: V0 * to2(V3(kap)) for V0, kap, lab in MINN}

NM = 40
orb = np.arange(NM) // 4; iso = (np.arange(NM) // 2) % 2; spn = np.arange(NM) % 2
ml = np.array([c3.HE4_SPATIAL[o][2] for o in orb]); ms2 = 1 - 2 * spn
nq = np.array([2 * c3.HE4_SPATIAL[o][0] + c3.HE4_SPATIAL[o][1] for o in orb])
def direct(p, q_, r, s):
    if iso[p] != iso[r] or iso[q_] != iso[s]: return 0.0
    op, oq, orr, os_ = orb[p], orb[q_], orb[r], orb[s]
    sd = float(spn[p] == spn[r] and spn[q_] == spn[s]); sx = float(spn[p] == spn[s] and spn[q_] == spn[r])
    val = 0.0
    for wr, (i, j, k, l) in ((U_EXCH / 2, (op, oq, orr, os_)), ((2 - U_EXCH) / 2, (op, oq, os_, orr))):
        vr = Vg['R'][i, j, k, l]; vt = Vg['t'][i, j, k, l]; vs = Vg['s'][i, j, k, l]
        val += wr * (vr * sd + 0.5 * vt * (sd + sx) + 0.5 * vs * (sd - sx))
    if sd: val += K2s[op, oq, orr, os_]
    return val

# ---------------- many-body basis: Z protons, N neutrons, fixed M_L, M_S ----------------
prot = [m for m in range(NM) if iso[m] == 0]; neut = [m for m in range(NM) if iso[m] == 1]
def species(modes, k):
    out = {}
    for c in itertools.combinations(modes, k):
        key = (int(ml[list(c)].sum()), int(ms2[list(c)].sum()), int(nq[list(c)].sum()))
        out.setdefault(key, []).append(sum(1 << m for m in c))
    return out
SP, SN = species(prot, Z), species(neut, N)
NMIN = min(k[2] for k in SP) + min(k[2] for k in SN)
def basis_for(ML, MS2):
    st = []
    for (a1, b1, n1), lp in SP.items():
        for (a2, b2, n2), ln in SN.items():
            if a1 + a2 == ML and b1 + b2 == MS2 and n1 + n2 <= NMIN + NEXTRA:
                st.extend(x | y for x in lp for y in ln)
    return np.array(sorted(st), dtype=np.int64)

TAB = np.array([bin(i).count('1') for i in range(256)], dtype=np.int64)
def popcount(v):
    return TAB[v.view(np.uint8).reshape(-1, 8)].sum(1)
def below(st, m):
    return popcount(st & np.int64((1 << m) - 1))

# antisymmetrized TBMEs, p<q, r<s, conserving iso, ml, ms
pairs = [(p, q_) for p in range(NM) for q_ in range(p + 1, NM)]
key = lambda p, q_: (iso[p] + iso[q_], ml[p] + ml[q_], ms2[p] + ms2[q_])
groups = {}
for pq in pairs: groups.setdefault(key(*pq), []).append(pq)
TB = {}
for g in groups.values():
    for (r, s) in g:
        lst = []
        for (p, q_) in g:
            v = direct(p, q_, r, s) - direct(p, q_, s, r)
            if abs(v) > 1e-11: lst.append((p, q_, v))
        TB[(r, s)] = lst
OB = [(p, r, t_sph[orb[p], orb[r]]) for p in range(NM) for r in range(NM)
      if iso[p] == iso[r] and spn[p] == spn[r] and abs(t_sph[orb[p], orb[r]]) > 1e-11]
print(f'{NAME}: TBME setup {time.time()-t0:.1f}s', flush=True)

def apply_ann(st, m):
    sg = np.where(below(st, m) % 2 == 0, 1, -1); return st & ~np.int64(1 << m), sg
def apply_cre(st, m):
    sg = np.where(below(st, m) % 2 == 0, 1, -1); return st | np.int64(1 << m), sg

def build_H(basis):
    dim = len(basis); rows, cols, vals = [], [], []
    idx_all = np.arange(dim)
    for (p, r, v) in OB:
        occ = ((basis >> r) & 1).astype(bool)
        st = basis[occ]; ix = idx_all[occ]
        st1, s1 = apply_ann(st, r)
        ok = ((st1 >> p) & 1) == 0
        st2, s2 = apply_cre(st1[ok], p)
        j = np.searchsorted(basis, st2); good = (j < dim)
        good[good] &= basis[j[good]] == st2[good]
        rows.append(j[good]); cols.append(ix[ok][good]); vals.append(v * (s1[ok] * s2)[good])
    for (r, s), lst in TB.items():
        if not lst: continue
        occ = (((basis >> r) & 1) & ((basis >> s) & 1)).astype(bool)
        if not occ.any(): continue
        st = basis[occ]; ix = idx_all[occ]
        st1, s1 = apply_ann(st, r); st1, s2 = apply_ann(st1, s); sg0 = s1 * s2
        for (p, q_, v) in lst:
            ok = (((st1 >> p) & 1) | ((st1 >> q_) & 1)) == 0
            if not ok.any(): continue
            st2, s3 = apply_cre(st1[ok], q_); st2, s4 = apply_cre(st2, p)
            j = np.searchsorted(basis, st2); good = j < dim
            good[good] &= basis[j[good]] == st2[good]
            rows.append(j[good]); cols.append(ix[ok][good]); vals.append(v * (sg0[ok] * s3 * s4)[good])
    H = sp.coo_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(dim, dim)).tocsr()
    return H

def lowest(ML, MS2, k=3):
    basis = basis_for(ML, MS2)
    if len(basis) == 0: return basis, None, None
    H = build_H(basis)
    herm = abs(H - H.conj().T).max() if H.nnz else 0
    if len(basis) < 600:
        w, v = np.linalg.eigh(H.toarray()); return basis, w[:k] - 1.5 * BETA * hw, v[:, 0]
    w, v = eigsh(H, k=k, which='SA', tol=1e-10)
    o = np.argsort(w)
    print(f'  sector ML={ML} 2MS={MS2} dim={len(basis)} nnz={H.nnz} herm={herm:.1e} E={w[o]-1.5*BETA*hw}', flush=True)
    return basis, w[o] - 1.5 * BETA * hw, v[:, o[0]]

# ---------------- ground state; check L, S by lower sectors ----------------
res = dict(name=NAME, Z=Z, N=N, A=A, L=LGS, S2=S2, J2=2 * LGS + S2, beta=BETA, nextra=NEXTRA, b=b, hw=hw)
basis, E, psi = lowest(LGS, S2)
res.update(dim=len(basis), E0=float(E[0]), E_sector=[float(x) for x in E])
checks = {}
s2min = S2 % 2
for (ML, MS2) in {(0, s2min), (LGS, s2min), (0, S2), (LGS + 1, S2), (LGS, S2 + 2)}:
    if (ML, MS2) == (LGS, S2): continue
    bb = basis_for(ML, MS2)
    if len(bb) == 0 or len(bb) > 1_500_000: continue
    _, Ek, _ = lowest(ML, MS2, k=2)
    checks[f'{ML},{MS2}'] = float(Ek[0])
res['sector_E0'] = checks
print(NAME, 'E0', E[0], 'checks', checks, flush=True)

# ---------------- one-body proxies ----------------
nsp = len(c3.HE4_SPATIAL)
def rad(n1, l1, n2, l2, L):
    f = lambda r: float(c3.ho_radial(n1, l1, r, b) * c3.ho_radial(n2, l2, r, b)) * spherical_jn(L, q * r) * r * r
    return integrate.quad(f, 0, 25 * b, limit=400, epsabs=1e-13, epsrel=1e-12)[0]
def Ysp(L, mu):           # <i| j_L(qr) Y_L mu |j> on spatial orbitals
    M_ = np.zeros((nsp, nsp), complex)
    for i, (n1, l1, m1) in enumerate(c3.HE4_SPATIAL):
        for j, (n2, l2, m2) in enumerate(c3.HE4_SPATIAL):
            if m1 != mu + m2: continue
            g = float(gaunt(l1, L, l2, -m1, mu, m2)) * (-1) ** m1
            if g: M_[i, j] = rad(n1, l1, n2, l2, L) * g
    return M_
def lsp(mu):              # j0(qr) l_mu (spherical components)
    M_ = np.zeros((nsp, nsp), complex)
    for i, (n1, l1, m1) in enumerate(c3.HE4_SPATIAL):
        for j, (n2, l2, m2) in enumerate(c3.HE4_SPATIAL):
            if l1 != l2 or l1 == 0: continue
            r0 = rad(n1, l1, n2, l2, 0)
            if mu == 0 and m1 == m2: M_[i, j] = r0 * m1
            if mu == 1 and m1 == m2 + 1: M_[i, j] = -r0 * math.sqrt(l1 * (l1 + 1) - m2 * (m2 + 1)) / math.sqrt(2)
            if mu == -1 and m1 == m2 - 1: M_[i, j] = r0 * math.sqrt(l1 * (l1 + 1) - m2 * (m2 - 1)) / math.sqrt(2)
    return M_
sz = np.array([[0.5, 0], [0, -0.5]], complex)       # spin index 0 = up
splus = np.array([[0, 1], [0, 0]], complex)
s_mu = {0: sz, 1: -splus / math.sqrt(2), -1: splus.T / math.sqrt(2)}
def cgf(j1, m1, j2, m2, J, M):
    return float(clebsch_gordan(Rational(j1), Rational(j2), Rational(J), Rational(m1), Rational(m2), Rational(M)))
def couple(Aop, Bop, k):  # [A^(1) x B^(1)]^k_0, A spatial dict mu->10x10, B spin dict mu->2x2
    return sum(cgf(ra, mu, 1, -mu, k, 0) * np.kron(Aop[mu], Bop[-mu]) for mu in (-1, 0, 1) for ra in [Aop['rank']]
               if abs(mu) <= ra)
def to_modes(sp_spin, species_t):     # 20x20 (orb x spin) one-body -> 40 modes, isospin species_t
    o = np.zeros((NM, NM), complex)
    for i in range(nsp):
        for si in range(2):
            for j in range(nsp):
                for sj in range(2):
                    o[4 * i + 2 * species_t + si, 4 * j + 2 * species_t + sj] = sp_spin[2 * i + si, 2 * j + sj]
    return o
I2 = np.eye(2)
Y2 = {mu: Ysp(2, mu) for mu in range(-2, 3)}
L1 = {mu: lsp(mu) for mu in (-1, 0, 1)}; L1['rank'] = 1
S1 = dict(s_mu)
# [Y2 x sigma]^1_0
Y2s1 = sum(cgf(2, mu, 1, -mu, 1, 0) * np.kron(Y2[mu], 2 * S1[-mu]) for mu in (-1, 0, 1))
ops_sp = {
    'M0': np.kron(Ysp(0, 0) * math.sqrt(4 * math.pi), I2),
    'M2': np.kron(Y2[0], I2),
    'S1L0': np.kron(Ysp(0, 0) * math.sqrt(4 * math.pi), 2 * sz),
    'S1L2': Y2s1,
    'D1': np.kron(L1[0], I2),
    'P0': couple(L1, S1, 0),
    'P2': couple(L1, S1, 2),
}
for k_, v in list(ops_sp.items()):
    if abs(v - v.conj().T).max() > 1e-10:
        if abs(v + v.conj().T).max() < 1e-10: ops_sp[k_] = 1j * v
        else: raise SystemExit(f'{k_} not (anti)Hermitian')

def jw_lambda(o):
    """JW one-norm of sum o_pq a+_p a_q, complex Hermitian: identity Tr/2, Z_p -o_pp/2, and for p<q
    Re o_pq (XX+YY)/2 and Im o_pq (XY-YX)/2 strings. Returns lambda (identity excluded), c0, Gamma."""
    c0 = float(np.trace(o).real) / 2; lam = 0.0; gam = 0
    for p in range(NM):
        if abs(o[p, p]) > 1e-12: lam += abs(o[p, p].real) / 2; gam += 1
        for q_ in range(p + 1, NM):
            re, im = o[p, q_].real, o[p, q_].imag
            if abs(re) > 1e-12: lam += abs(re); gam += 2
            if abs(im) > 1e-12: lam += abs(im); gam += 2
    return lam, c0, gam

def density(vec, basis):
    """rho[p, r] = <a+_p a_r> in vec."""
    rho = np.zeros((NM, NM), complex); dim = len(basis); idx = np.arange(dim)
    for r in range(NM):
        occ = ((basis >> r) & 1).astype(bool); st = basis[occ]; ix = idx[occ]
        st1, s1 = apply_ann(st, r)
        for p in range(NM):
            if iso[p] != iso[r]: continue
            ok = ((st1 >> p) & 1) == 0
            st2, s2 = apply_cre(st1[ok], p)
            j = np.searchsorted(basis, st2); good = j < dim
            good[good] &= basis[j[good]] == st2[good]
            rho[p, r] = np.vdot(vec[j[good]], vec[ix[ok][good]] * (s1[ok] * s2)[good])
    return rho / np.vdot(vec, vec).real

rho = density(psi, basis)
res['Zcheck'] = float(sum(rho[p, p].real for p in prot)); res['Ncheck'] = float(sum(rho[p, p].real for p in neut))
order = np.argsort(-np.abs(psi) ** 2); cum = np.cumsum(np.abs(psi[order]) ** 2)
Dfor = lambda thr: int(np.searchsorted(1 - cum, thr, side='right') == 0) or int(np.argmax(1 - cum <= thr) + 1)
res['w0'] = float(cum[0])
res['D_infid_0.01'] = int(np.argmax(1 - cum <= 1e-2) + 1); res['D_infid_0.001'] = int(np.argmax(1 - cum <= 1e-3) + 1)
res['nonzero'] = int(np.sum(np.abs(psi) ** 2 > 1e-14))
rhoD = {}
DSCAN = sorted({res['D_infid_0.01'], res['D_infid_0.001']} | {d for d in (50, 100, 150, 200, 300, 400, 500, 596, 800, 1000, 1500)
                                                          if d < res['nonzero']})
res['infid_at_D'] = {str(D): float(1 - cum[D - 1]) for D in DSCAN}
for D in DSCAN:
    v = np.zeros_like(psi); v[order[:D]] = psi[order[:D]]; rhoD[D] = density(v, basis)
ops = {}
for k_, o_sp in ops_sp.items():
    for t, tn in ((0, 'p'), (1, 'n')):
        o = to_modes(o_sp, t)
        lam, c0, gam = jw_lambda(o)
        ex = float(np.sum(o * rho.T).real)           # sum_pr o_pr <a+_p a_r>
        exD = {D: float(np.sum(o * rr.T).real) for D, rr in rhoD.items()}
        ops[f'{k_}_{tn}'] = dict(lam=lam, c0=c0, gamma=gam, exp=ex,
                                 shot_factor=(lam / abs(ex)) ** 2 if abs(ex) > 1e-9 else None,
                                 rel_err_at_D={str(D): (abs(e - ex) / abs(ex) if abs(ex) > 1e-9 else None) for D, e in exD.items()})
res['ops'] = ops
res['seconds'] = time.time() - t0
print(json.dumps(res, indent=1))
json.dump(res, open(__import__('pathlib').Path(__file__).resolve().parent / f'light_{NAME}_b{BETA}_n{NEXTRA}.json', 'w'), indent=1)
