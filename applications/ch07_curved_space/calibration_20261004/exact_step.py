"""U2 exact step-convergence check for Ch.7 (2033 preheating instance), reduced instance.

Units m_phi = 1. Chapter values: b = 3.46, H_inf = 0.01, T = 1000, a = (1 + 1.5 H t)^(2/3), K = 16 uniform grid
(arXiv:2108.10793, spacing sqrt(2 pi / K) in units of the a = 1 site vacuum width), Phi0 = 0.31, m_psi = 0.1, r = 1,
(lam, g) in {(1,1), (3,0.5), (1,2)}. Site normalisation v = b^3 (the 3D site volume, so quanta per site, the
digitisation and the per-site couplings match the 5^3 instance).

Reduced lattice: 1D periodic ring of L = 3 sites (k = 0, +-2 pi/(3b)), real scalar on K = 16 per site, one
2-component Wilson fermion (gamma0 = s3, gamma0 gamma1 = s1), sector at half filling (C(6,3) = 20).
Hilbert dim 20 * 16^3 = 81920. Exact many-body evolution under the chapter's splitting:
  S(h) = D(h/2) [Kin(h) Hop(h)] D(h/2), coefficients at a(t + h/2);
D = everything diagonal in the field/occupation basis (mass, phi^4, gradient, Wilson on-site, bare mass, Yukawa),
Kin = per-site pi^2 (exact via the centred DFT), Hop = fermion hop (exact, 20 x 20). Kin and Hop act on different
registers and commute.

enc = "fixed": phi on the fixed grid set at a = 1 (the chapter's circuit as written).
enc = "chi":   chi = a^{3/2} phi on the grid (width fixed); adds the dilation term (3/2)(adot/a)(x p + p x)/2, applied
               as an extra palindromic block D/2 Dil/2 [Kin Hop] Dil/2 D/2 (not in the chapter's per-step count).

Initial state: ground state of the digitised interacting H at a = 1 (the ramp target; the map vacuum differs at
O((w h)^4) ~ 1e-7 in n_k), then the boost exp(i pbar sum_x x_x), pbar = sqrt(v) Phi0.
Readout (same for every h): scalar quasiparticle n_k with W^2 = 1 + 3 lam phibar^2 + 4 sin^2(kb/2)/(a b)^2;
fermion particle number in the positive-energy out-modes of h(k) at M = m_psi + g phibar.
"""
import sys, json, math, time, itertools
import numpy as np
from scipy.sparse.linalg import eigsh, LinearOperator

B = 3.46; V = B ** 3; R = 1.0; MPSI = 0.1; HINF = 0.01; TEND = 1000.0; PHI0 = 0.31
L = 3; K = 16
REC_T = [12.5, 25, 50, 100, 200, 400, 700, 1000]


def a_of(t):
    return (1 + 1.5 * HINF * t) ** (2 / 3)


def adot_over_a(t):
    return HINF / (1 + 1.5 * HINF * t)


# ---------------- scalar grid --------------------------------------------------
dx = math.sqrt(2 * math.pi / K)
jj = np.arange(K) - (K - 1) / 2
X = jj * dx
F = np.exp(-1j * np.outer(jj, jj) * 2 * np.pi / K) / np.sqrt(K)   # x -> p basis (p grid = x grid)
P2 = F.conj().T @ np.diag(X ** 2) @ F
Pm = F.conj().T @ np.diag(X) @ F
G = (np.diag(X) @ Pm + Pm @ np.diag(X)) / 2
G = (G + G.conj().T) / 2
gw, gV = np.linalg.eigh(G)

# ---------------- fermions -----------------------------------------------------
NM = 2 * L
def mode(x, s):
    return 2 * x + s

states = [s for s in range(2 ** NM) if bin(s).count("1") == L]
sidx = {s: i for i, s in enumerate(states)}
NF = len(states)

def cdag_c(i, j):
    """many-body matrix of c_i^dag c_j in the half-filled sector (JW, mode 0 first)."""
    M = np.zeros((NF, NF))
    for s in states:
        if not (s >> j) & 1:
            continue
        s1 = s ^ (1 << j)
        sg = (-1) ** bin(s & ((1 << j) - 1)).count("1")
        if (s1 >> i) & 1:
            continue
        s2 = s1 | (1 << i)
        sg *= (-1) ** bin(s1 & ((1 << i) - 1)).count("1")
        M[sidx[s2], sidx[s]] += sg
    return M

CC = [[cdag_c(i, j) for j in range(NM)] for i in range(NM)]
nmode = np.array([[(s >> i) & 1 for i in range(NM)] for s in states], float)   # (NF, NM)
cbc = np.stack([nmode[:, mode(x, 0)] - nmode[:, mode(x, 1)] for x in range(L)], 1)  # (NF, L): psibar psi per site

s1 = np.array([[0, 1], [1, 0.]]); s3 = np.array([[1, 0], [0, -1.]])
T = np.roll(np.eye(L), 1, axis=1)          # (T c)_x = c_{x+1}
h_hop1 = (-(R / 2) * np.kron(T + T.T, s3) + np.kron((T - T.T) / 2j, s1)) / B   # times 1/a
H_hop1 = sum(h_hop1[i, j] * CC[i][j] for i in range(NM) for j in range(NM))
H_hop1 = (H_hop1 + H_hop1.conj().T) / 2
hw, hV = np.linalg.eigh(H_hop1)

ks = 2 * np.pi * np.arange(L) / (L * B)
ks[ks > np.pi / B] -= 2 * np.pi / B


def hk(k, a, M):
    return (M + R / (a * B) * (1 - np.cos(k * B))) * s3 + np.sin(k * B) / (a * B) * s1


def pos_mode(k, a, M):
    w, U = np.linalg.eigh(hk(k, a, M))
    return U[:, 1], U[:, 0]


# ---------------- diagonal pieces (shape NF, K^L) ------------------------------
xs = np.stack(np.meshgrid(*([X] * L), indexing="ij"), 0).reshape(L, -1)   # (L, K^L)
grad = sum((xs[x] - xs[(x + 1) % L]) ** 2 for x in range(L)) / (2 * B ** 2)
quad = (xs ** 2).sum(0) / 2
quart = (xs ** 4).sum(0) / (4 * V)
wil = cbc.sum(1)                                # sum_x psibar psi  (NF,)
yuk = cbc @ xs / math.sqrt(V)                   # (NF, K^L): sum_x psibar psi x_x / sqrt(v)


def diag(a, lam, g, enc):
    if enc == "fixed":
        sc = a ** 3 * (quad + lam * quart) + a * grad
        return sc[None, :] + ((MPSI + R / (a * B)) * wil)[:, None] + g * yuk
    a3 = a ** 3
    sc = quad + lam * quart / a3 + grad / a ** 2
    return sc[None, :] + ((MPSI + R / (a * B)) * wil)[:, None] + g * yuk / a ** 1.5


def apply_site(psi, U):
    """apply single-site K x K matrix U on every scalar site. psi shape (NF, K, K, K)."""
    for ax in range(1, L + 1):
        psi = np.moveaxis(np.tensordot(U, psi, axes=([1], [ax])), 0, ax)
    return psi


def sum_site(psi, O):
    """sum_x O_x psi (one-site operator O summed over sites)."""
    out = np.zeros_like(psi)
    for ax in range(1, L + 1):
        out = out + np.moveaxis(np.tensordot(O, psi, axes=([1], [ax])), 0, ax)
    return out


def kin_U(theta):
    return F.conj().T @ (np.exp(-1j * theta * X ** 2)[:, None] * F)


def dil_U(theta):
    return gV @ (np.exp(-1j * theta * gw)[:, None] * gV.conj().T)


def hop_U(theta):
    return hV @ (np.exp(-1j * theta * hw)[:, None] * hV.conj().T)


def ground(lam, g, enc):
    d = diag(1.0, lam, g, enc).ravel()
    shape = (NF,) + (K,) * L
    Kin = P2 / 2

    def mv(v):
        psi = v.reshape(shape)
        out = d * v
        out = out + sum_site(psi, Kin).ravel()
        out = out + np.tensordot(H_hop1, psi, axes=([1], [0])).ravel()
        return out
    op = LinearOperator((d.size, d.size), matvec=mv, dtype=complex)
    w, vv = eigsh(op, k=1, which="SA", tol=1e-12)
    return w[0], vv[:, 0].reshape(shape)


def observe(psi, t, lam, g, enc):
    a = a_of(t)
    shape = psi.shape
    prob = abs(psi) ** 2
    ps = prob.sum(0).reshape(-1)
    xm = xs @ ps                                     # <x_x>
    Cx = (xs * ps) @ xs.T - np.outer(xm, xm)
    # p basis
    q = apply_site(psi, F)
    pq = (abs(q) ** 2).sum(0).reshape(-1)
    pm = xs @ pq
    Cp = (xs * pq) @ xs.T - np.outer(pm, pm)
    xbar = xm.mean()
    if enc == "fixed":
        phibar = xbar / math.sqrt(V)
        sx, sp = a ** 3, 1 / a ** 3
    else:
        phibar = xbar / math.sqrt(V) / a ** 1.5
        sx, sp = 1.0, 1.0
    ns = []
    for k in ks:
        if abs(k) < 1e-12:
            continue
        ph = np.exp(-1j * k * B * np.arange(L))
        xk = np.real(ph.conj() @ Cx @ ph) / L
        pk = np.real(ph.conj() @ Cp @ ph) / L
        W = math.sqrt(1 + 3 * lam * phibar ** 2 + 4 * math.sin(k * B / 2) ** 2 / (a * B) ** 2)
        ns.append((sp * pk + sx * W ** 2 * xk) / (2 * W) - 0.5)
    # fermion one-body density
    Mf = psi.reshape(NF, -1)
    RHO = Mf @ Mf.conj().T                           # RHO_{b a} = sum psi_b psi_a^*
    rho = np.array([[np.trace(CC[i][j] @ RHO) for j in range(NM)] for i in range(NM)])  # <c_i^dag c_j>
    M = MPSI + g * phibar
    nf = {}
    for k in ks:
        up, _ = pos_mode(k, a, M)
        w = np.kron(np.exp(1j * k * B * np.arange(L)) / math.sqrt(L), up)
        nf.setdefault(round(abs(k), 8), []).append(float(np.real(w @ rho @ w.conj())))
    nfs = [float(np.mean(v)) for kk, v in sorted(nf.items())]
    return dict(t=t, a=a, phibar=float(phibar), xbar=float(xbar), n_scalar=float(np.mean(ns)), n_ferm=nfs,
                norm=float(prob.sum()))


def run(lam, g, h, enc, psi0):
    psi = psi0.copy()
    N = int(round(TEND / h))
    rec = {int(round(tr / h)): tr for tr in REC_T}
    out = []
    t = 0.0
    hwop = None
    for n in range(1, N + 1):
        am = a_of(t + h / 2)
        ph = np.exp(-1j * (h / 2) * diag(am, lam, g, enc)).reshape(psi.shape)
        psi = psi * ph
        if enc == "chi":
            Ud = dil_U(1.5 * adot_over_a(t + h / 2) * h / 2)
            psi = apply_site(psi, Ud)
            psi = apply_site(psi, kin_U(h / 2))
        else:
            psi = apply_site(psi, kin_U(h / (2 * am ** 3)))
        psi = np.tensordot(hop_U(h / am), psi, axes=([1], [0]))
        if enc == "chi":
            psi = apply_site(psi, Ud)
        psi = psi * ph
        t = n * h
        if n in rec:
            out.append(observe(psi, rec[n], lam, g, enc))
    return out


if __name__ == "__main__":
    lam, g, enc = float(sys.argv[1]), float(sys.argv[2]), sys.argv[3]
    hs = [float(x) for x in sys.argv[4].split(",")]
    e0, psi = ground(lam, g, enc)
    pbar = math.sqrt(V) * PHI0
    boost = np.exp(1j * pbar * xs.sum(0)).reshape((1,) + (K,) * L)
    psi0 = psi * boost
    res = {"lam": lam, "g": g, "enc": enc, "E0": e0, "init": observe(psi0, 0.0, lam, g, enc), "runs": {}}
    for h in hs:
        t0 = time.time()
        res["runs"][str(h)] = run(lam, g, h, enc, psi0)
        res["runs"][str(h) + "_sec"] = time.time() - t0
        json.dump(res, open(f"exact_{enc}_l{lam}_g{g}.json", "w"), indent=1)
        print(h, time.time() - t0, res["runs"][str(h)][-1], flush=True)
