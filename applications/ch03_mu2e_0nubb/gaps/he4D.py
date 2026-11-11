"""M3: data size D of the 4He ground state in the Ch. 2/3 active space (0s, 0p, 1s, 0d; 40 JW modes),
with the Minnesota NN potential (Thompson, LeMere, Tang, NPA 286 (1977) 53; u = 1, no Coulomb), stated as a
schematic stand-in for the chiral interaction, intrinsic kinetic energy, b from the chapter's model (1.369 fm).
D = number of largest-|c| Slater determinants (JW occupation configurations) the oracle load must write.
"""
import sys, math, itertools, json
import numpy as np
from numpy.polynomial.hermite import hermgauss
from scipy.special import eval_hermite, sph_harm_y
sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[3]))   # repository root
from estimates import ch03_mu2e_0nubb as c3

a = c3.Assumptions()
osc = c3.he4_oscillator_length(a)
b, hw = osc['b'], osc['hbar_omega']
A = 4
HBARC = 197.3269804
U_EXCH = float(sys.argv[1]) if len(sys.argv) > 1 else 1.0
SCALE = float(sys.argv[2]) if len(sys.argv) > 2 else 1.0     # overall strength factor (sensitivity)
BETA = float(sys.argv[3]) if len(sys.argv) > 3 else 0.0      # Lawson c.m. term beta (H_cm - 3/2 hw)
MINN = [(200.0, 1.487, 'R'), (-178.0, 0.639, 't'), (-91.85, 0.465, 's')]

# --- Cartesian HO 1D functions in u = x/b, normalized w.r.t. du (b drops out of dimensionless overlaps)
cart = [(nx, ny, nz) for N in range(3) for nx in range(N + 1) for ny in range(N + 1 - nx) for nz in [N - nx - ny]]
assert len(cart) == 10
def phi1(n, u):
    return eval_hermite(n, u) * np.exp(-u * u / 2) / math.sqrt(2 ** n * math.factorial(n) * math.sqrt(math.pi))

# 1D Gaussian two-body integrals G[a,b,c,d] = int phi_a(u1) phi_b(u2) e^{-k b^2 (u1-u2)^2} phi_c(u1) phi_d(u2)
xg, wg = hermgauss(80)
U1, U2 = np.meshgrid(xg, xg, indexing='ij')
W = np.outer(wg, wg)
P = np.array([phi1(n, xg) * np.exp(xg * xg / 2) for n in range(3)])   # strip weight
def G1(kap):
    g = np.exp(-kap * b * b * (U1 - U2) ** 2) * W
    return np.einsum('ai,bj,ci,dj,ij->abcd', P, P, P, P, g)
# check: <00|e^{-k r^2}|00> 1D = 1/sqrt(1 + 2 k b^2)
for V0, kap, _ in MINN:
    assert abs(G1(kap)[0, 0, 0, 0] - 1 / math.sqrt(1 + 2 * kap * b * b)) < 1e-12

def V3(kap):
    g = G1(kap)
    n = len(cart)
    out = np.empty((n, n, n, n))
    for i, ci in enumerate(cart):
        for j, cj in enumerate(cart):
            for k, ck in enumerate(cart):
                for l, cl in enumerate(cart):
                    out[i, j, k, l] = np.prod([g[ci[d], cj[d], ck[d], cl[d]] for d in range(3)])
    return out

# 1D dimensionless p~ = -i d/du (p = hbar/b p~), and p~^2/2 matrices on n = 0..2 (exact one-body elements)
pt = np.zeros((3, 3), complex)
for n in range(3):
    if n + 1 < 3: pt[n + 1, n] = 1j * math.sqrt((n + 1) / 2)
    if n - 1 >= 0: pt[n - 1, n] = -1j * math.sqrt(n / 2)
def p2half(m, n):
    if m == n: return (2 * n + 1) / 4
    if abs(m - n) == 2: return -math.sqrt(max(m, n) * (max(m, n) - 1)) / 4
    return 0.0
ncart = len(cart)
T1 = np.zeros((ncart, ncart))          # sum_d p~_d^2/2 in units of hbar omega
for i, ci in enumerate(cart):
    for j, cj in enumerate(cart):
        for d in range(3):
            if all(ci[e] == cj[e] for e in range(3) if e != d):
                T1[i, j] += p2half(ci[d], cj[d])
Pd = []
for d in range(3):
    m = np.zeros((ncart, ncart), complex)
    for i, ci in enumerate(cart):
        for j, cj in enumerate(cart):
            if all(ci[e] == cj[e] for e in range(3) if e != d):
                m[i, j] = pt[ci[d], cj[d]]
    Pd.append(m)

# --- Cartesian -> spherical (chapter's HE4_SPATIAL order, R_nl positive at origin, Condon-Shortley Y_lm)
xs, ws = hermgauss(12)
X, Y, Z = np.meshgrid(xs, xs, xs, indexing='ij')
Wt = (ws[:, None, None] * ws[None, :, None] * ws[None, None, :]).ravel()
X, Y, Z = X.ravel(), Y.ravel(), Z.ravel()
R = np.sqrt(X * X + Y * Y + Z * Z)
TH = np.arccos(Z / R); PH = np.arctan2(Y, X)
cart_vals = np.array([phi1(c[0], X) * phi1(c[1], Y) * phi1(c[2], Z) for c in cart]) * np.exp(R * R / 2)
sph_vals = []
for (n, l, m) in c3.HE4_SPATIAL:
    rad = c3.ho_radial(n, l, R, 1.0)
    sph_vals.append(rad * sph_harm_y(l, m, TH, PH) * np.exp(R * R / 2))
sph_vals = np.array(sph_vals)
Umat = np.einsum('sk,ck,k->sc', sph_vals, cart_vals, Wt)        # psi_s = sum_c U[s,c] phi_c (phi real)
assert np.allclose(Umat @ Umat.conj().T, np.eye(10), atol=1e-10)

def to_sph1(h):
    return Umat.conj() @ h @ Umat.T
def to_sph2(v):
    return np.einsum('ia,jb,kc,ld,abcd->ijkl', Umat.conj(), Umat.conj(), Umat, Umat, v)

Rd = []
for d in range(3):
    m_ = np.zeros((ncart, ncart))
    for i, ci in enumerate(cart):
        for j, cj in enumerate(cart):
            if all(ci[e] == cj[e] for e in range(3) if e != d):
                if ci[d] == cj[d] + 1: m_[i, j] = math.sqrt(ci[d] / 2)
                if ci[d] == cj[d] - 1: m_[i, j] = math.sqrt(cj[d] / 2)
    Rd.append(m_)
Ncart = np.diag([sum(c) for c in cart]).astype(float)
t_sph = to_sph1(T1 * hw * (A - 1) / A + BETA * hw / A * (Ncart + 1.5 * np.eye(ncart)))
K2 = sum(np.einsum('ac,bd->abcd', Pd[d], Pd[d]) for d in range(3)) * (-hw / A)   # -(1/(mA)) p1.p2
K2 = K2 + BETA * hw / A * sum(np.einsum('ac,bd->abcd', Pd[d], Pd[d]) + np.einsum('ac,bd->abcd', Rd[d], Rd[d]) for d in range(3))
K2s = to_sph2(K2)
Vg = {lab: SCALE * V0 * to_sph2(V3(kap)) for V0, kap, lab in MINN}

# --- modes p = 4 orb + 2 iso + spin (iso 0 = proton)
NM = 40
orb = np.arange(NM) // 4; iso = (np.arange(NM) // 2) % 2; spn = np.arange(NM) % 2
ml = np.array([c3.HE4_SPATIAL[o][2] for o in orb]); ms = 1 - 2 * spn
def direct(p, q, r, s):
    """<pq|V|rs> (unsymmetrized), V = [V_R + (1+Ps)/2 V_t + (1-Ps)/2 V_s] (u + (2-u) P^r)/2 + kinetic 2-body."""
    if iso[p] != iso[r] or iso[q] != iso[s]:
        return 0.0
    op, oq, orr, os_ = orb[p], orb[q], orb[r], orb[s]
    sd = float((spn[p] == spn[r]) and (spn[q] == spn[s]))
    sx = float((spn[p] == spn[s]) and (spn[q] == spn[r]))
    val = 0.0
    # spatial direct and spatial-exchanged (P^r) parts
    for wr, (i, j, k, l) in ((U_EXCH / 2, (op, oq, orr, os_)), ((2 - U_EXCH) / 2, (op, oq, os_, orr))):
        vr = Vg['R'][i, j, k, l]; vt = Vg['t'][i, j, k, l]; vs = Vg['s'][i, j, k, l]
        val += wr * (vr * sd + 0.5 * vt * (sd + sx) + 0.5 * vs * (sd - sx))
    if sd:
        val += K2s[op, oq, orr, os_]
    return val

# --- 2p2n basis restricted to M_L = 0, M_S = 0
prot = [m for m in range(NM) if iso[m] == 0]; neut = [m for m in range(NM) if iso[m] == 1]
pp = list(itertools.combinations(prot, 2)); nn = list(itertools.combinations(neut, 2))
basis = [(P, Nn) for P in pp for Nn in nn
         if sum(ml[list(P + Nn)]) == 0 and sum(ms[list(P + Nn)]) == 0]
dim = len(basis)
print('b', b, 'hw', hw, 'u', U_EXCH, 'scale', SCALE, 'sector dim', dim)
index = {bb: i for i, bb in enumerate(basis)}

def sign_sort(lst):
    lst = list(lst); s = 1
    for i in range(len(lst)):
        for j in range(len(lst) - 1 - i):
            if lst[j] > lst[j + 1]:
                lst[j], lst[j + 1] = lst[j + 1], lst[j]; s = -s
    return s, tuple(lst)

# antisymmetrized elements cache
from functools import lru_cache
@lru_cache(maxsize=None)
def vas(p, q, r, s):
    return direct(p, q, r, s) - direct(p, q, s, r)

H = np.zeros((dim, dim), complex)
for i, (P, Nn) in enumerate(basis):
    occ = tuple(sorted(P + Nn))
    for j, (P2, N2) in enumerate(basis):
        occ2 = tuple(sorted(P2 + N2))
        sa, sb = set(occ), set(occ2)
        diff_a = sorted(sa - sb); diff_b = sorted(sb - sa)
        nd = len(diff_a)
        if nd > 2: continue
        # Slater-Condon with explicit phases via second quantization on small tuples
        val = 0.0
        if nd == 0:
            for p in occ:
                val += t_sph[orb[p], orb[p]] if (iso[p] == iso[p]) else 0
            for x, y in itertools.combinations(occ, 2):
                val += vas(x, y, x, y)
        else:
            # generic: apply a+ a (nd=1) or a+a+aa (nd=2) to |occ2>, compare with <occ|
            def amp(cre, ann):
                st = list(occ2); s = 1
                for m in ann:
                    k = st.index(m); s *= (-1) ** k; st.pop(k)
                for m in reversed(cre):
                    k = sum(1 for x in st if x < m); s *= (-1) ** k; st.insert(k, m)
                return s if tuple(st) == tuple(sorted(occ)) else 0
            if nd == 1:
                p, r = diff_a[0], diff_b[0]
                if iso[p] != iso[r] or spn[p] != spn[r]:
                    pass
                else:
                    s1 = amp([p], [r])
                    if s1:
                        val += s1 * t_sph[orb[p], orb[r]]
                for x in sorted(sa & sb):
                    s2 = amp([p, x], [r, x])      # a+_p a+_x a_x a_r
                    if s2:
                        val += s2 * vas(p, x, r, x)
            else:
                p, q = diff_a; r, s_ = diff_b
                s2 = amp([p, q], [r, s_])         # a+_p a+_q a_s a_r (a_r first)
                if s2:
                    val += s2 * vas(p, q, r, s_)
        H[i, j] = val
print('herm err', np.abs(H - H.conj().T).max(), 'imag max', np.abs(H.imag).max())
w, v = np.linalg.eigh(H)
psi = v[:, 0]
print('E0', w[0] - 1.5 * BETA * hw, 'E1', w[1] - 1.5 * BETA * hw, 'gap', w[1] - w[0], 'beta', BETA)

# dipole: one-body proton operator from the chapter's model, same mode order
o, labels, _ = c3.he4_dipole_one_body(a)
def expO(vec):
    tot = 0.0
    nrm = np.vdot(vec, vec).real
    for i, (P, Nn) in enumerate(basis):
        if vec[i] == 0: continue
        for j, (P2, N2) in enumerate(basis):
            if vec[j] == 0 or Nn != N2: continue
            sa, sb = set(P), set(P2)
            if sa == sb:
                tot += (vec[i].conjugate() * vec[j]).real * sum(o[p, p] for p in P)
            elif len(sa - sb) == 1:
                p = (sa - sb).pop(); r = (sb - sa).pop()
                if o[p, r] == 0: continue
                st = sorted(P2 + N2); k = st.index(r); s = (-1) ** k; st.pop(k)
                k2 = sum(1 for x in st if x < p); s *= (-1) ** k2; st.insert(k2, p)
                tot += (vec[i].conjugate() * vec[j]).real * s * o[p, r]
    return tot / nrm
O_ex = expO(psi)
order = np.argsort(-np.abs(psi) ** 2)
cum = np.cumsum(np.abs(psi[order]) ** 2)
print('<O> exact', O_ex, ' leading weight', cum[0], ' (0s)^4 weight check')
res = dict(b=b, hw=hw, u=U_EXCH, scale=SCALE, beta=BETA, dim=dim, E0=float(w[0] - 1.5 * BETA * hw), gap=float(w[1]-w[0]), O_exact=float(O_ex), w0=float(cum[0]))
def D_for(crit, thr):
    for D in range(1, dim + 1):
        if crit(D) <= thr:
            return D
    return dim
nz = int(np.sum(np.abs(psi) ** 2 > 1e-14))
res['nonzero'] = nz
infid = lambda D: 1 - cum[D - 1]
def oerr(D):
    vec = np.zeros_like(psi); vec[order[:D]] = psi[order[:D]]
    return abs(expO(vec) - O_ex) / abs(O_ex)
for thr in (1e-2, 1e-3, 1e-4):
    res[f'D_infid_{thr:g}'] = D_for(infid, thr)
cache = {}
for thr in (1e-2, 1e-3):
    # scan D in steps (expO is O(dim^2)); report first D meeting the threshold
    Dfound = None
    for D in sorted(set([1, 2, 4, 8, 16, 32, 64] + list(range(1, 400)))):
        if D > nz: break
        e = oerr(D); cache[D] = e
        if e <= thr:
            Dfound = D; break
    res[f'D_Oerr_{thr:g}'] = Dfound
print('nonzero amplitudes', nz)
print({k: v for k, v in res.items()})
print('infidelity at D=1,10,100,600:', [float(1 - cum[min(D, dim) - 1]) for D in (1, 10, 100, 600)])
json.dump(res, open(f'he4D_u{U_EXCH}_s{SCALE}_b{BETA}.json', 'w'), indent=1)
