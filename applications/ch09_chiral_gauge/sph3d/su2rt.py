"""Classical real-time SU(2) lattice gauge theory (Kogut-Susskind / Ambjorn-Krasnitz), pure gauge, with a cooled
(topological) Chern-Simons number after Moore-Turok / Moore-Rummukainen (hep-ph/9906259), and the naive (Kubo)
N_CS = integral of E.B along the real-time path, which is what a quantum counting-field readout measures.

Conventions
  U = u0 + i u.sigma (unit quaternion), U = exp(i theta.sigma/2). Left electric generator e on each link.
  Classical canonical ensemble exp(-beta_L H_L), H_L = sum_links e^2/2 + sum_plaq (1 - Tr U_p / 2),
  beta_L = 4/(g^2 a T) at tree level.
  EOM in internal time tau: dU/dtau = i (e.sigma/2) U, de/dtau = -dV/d(left), tau = 2 t / a.
  Gauss: G_x = sum_i [ e_i(x) - R(U_i(x-i))^T e_i(x-i) ] = 0.
  dN_CS = (1/8 pi^2) sum_links dtheta_i(x) . Bbar_i(x), Bbar = clover plaquette angle (2 x vector part),
  averaged over the two ends of the link (the far end transported back). A winding-1 vacuum reached along the
  link-wise geodesic gives |dN| = 0.81, 0.92, 0.95 on 8^3, 12^3, 16^3 (tests_su2.py).
Array layout (component axis first, for speed): U (4, B, 3, N, N, N), e (3, B, 3, N, N, N). Site arrays
(c, B, N, N, N).
"""
import numpy as np, math


def qmul(a, b):
    a0, a1, a2, a3 = a[0], a[1], a[2], a[3]
    b0, b1, b2, b3 = b[0], b[1], b[2], b[3]
    return np.stack([a0 * b0 - a1 * b1 - a2 * b2 - a3 * b3,
                     a0 * b1 + b0 * a1 - (a2 * b3 - a3 * b2),
                     a0 * b2 + b0 * a2 - (a3 * b1 - a1 * b3),
                     a0 * b3 + b0 * a3 - (a1 * b2 - a2 * b1)])


def qdag(a):
    return np.stack([a[0], -a[1], -a[2], -a[3]])


def qexp(v):
    """exp(i v.sigma/2) for algebra vector v (3, ...)."""
    n = np.sqrt(v[0] ** 2 + v[1] ** 2 + v[2] ** 2)
    h = 0.5 * n
    s = np.where(n > 1e-14, np.sin(h) / np.where(n > 1e-14, n, 1.0), 0.5)
    return np.stack([np.cos(h), s * v[0], s * v[1], s * v[2]])


def qlog(q):
    w = q[1:]
    sn = np.sqrt(w[0] ** 2 + w[1] ** 2 + w[2] ** 2)
    ang = np.arctan2(sn, q[0])
    f = np.where(sn > 1e-14, 2 * ang / np.where(sn > 1e-14, sn, 1.0), 2.0)
    return f * w


def rot(U, v):
    """R(U) v with U (v.sigma) U^dag = (R v).sigma."""
    q0, q1, q2, q3 = U[0], U[1], U[2], U[3]
    x, y, z = v[0], v[1], v[2]
    return np.stack([(q0 * q0 + q1 * q1 - q2 * q2 - q3 * q3) * x + 2 * (q1 * q2 + q0 * q3) * y + 2 * (q1 * q3 - q0 * q2) * z,
                     2 * (q1 * q2 - q0 * q3) * x + (q0 * q0 - q1 * q1 + q2 * q2 - q3 * q3) * y + 2 * (q2 * q3 + q0 * q1) * z,
                     2 * (q1 * q3 + q0 * q2) * x + 2 * (q2 * q3 - q0 * q1) * y + (q0 * q0 - q1 * q1 - q2 * q2 + q3 * q3) * z])


def rotT(U, v):
    """R(U)^T v."""
    q0, q1, q2, q3 = U[0], U[1], U[2], U[3]
    x, y, z = v[0], v[1], v[2]
    return np.stack([(q0 * q0 + q1 * q1 - q2 * q2 - q3 * q3) * x + 2 * (q1 * q2 - q0 * q3) * y + 2 * (q1 * q3 + q0 * q2) * z,
                     2 * (q1 * q2 + q0 * q3) * x + (q0 * q0 - q1 * q1 + q2 * q2 - q3 * q3) * y + 2 * (q2 * q3 - q0 * q1) * z,
                     2 * (q1 * q3 - q0 * q2) * x + 2 * (q2 * q3 + q0 * q1) * y + (q0 * q0 - q1 * q1 - q2 * q2 + q3 * q3) * z])


def sh(X, mu, s):
    """X(x + s mu) for site arrays (c, B, N, N, N)."""
    return np.roll(X, -s, axis=2 + mu)


def staples(U):
    """W_i(x) = U_i(x) * sum of staples, shape (4, B, 3, N, N, N)."""
    W = np.empty_like(U)
    for i in range(3):
        Ui = U[:, :, i]
        S = 0.0
        for j in range(3):
            if j == i:
                continue
            Uj = U[:, :, j]
            S = S + qmul(qmul(sh(Uj, i, 1), qdag(sh(Ui, j, 1))), qdag(Uj))
            Ujm = sh(Uj, j, -1)
            S = S + qmul(qmul(qdag(sh(Ujm, i, 1)), qdag(sh(Ui, j, -1))), Ujm)
        W[:, :, i] = qmul(Ui, S)
    return W


def force(U):
    return 0.5 * staples(U)[1:]


def plaq_energy(U):
    W = staples(U)
    return (4 * 3 * np.prod(U.shape[3:6]) - W[0].sum((1, 2, 3, 4))) / 4.0


def energy(U, e):
    return 0.5 * (e * e).sum((0, 2, 3, 4, 5)) + plaq_energy(U)


def Cop(U, e):
    G = 0.0
    for i in range(3):
        G = G + e[:, :, i] - sh(rotT(U[:, :, i], e[:, :, i]), i, -1)
    return G


def CT(U, chi):
    out = np.empty((3,) + U.shape[1:])
    for i in range(3):
        out[:, :, i] = chi - rot(U[:, :, i], sh(chi, i, 1))
    return out


def gauss(U, e):
    return Cop(U, e)


def project_gauss(U, e, tol=1e-10, maxit=500):
    b = Cop(U, e)
    x = np.zeros_like(b)
    r = b.copy(); p = r.copy()
    rr = (r * r).sum((0, 2, 3, 4))
    b2 = rr.copy() + 1e-300
    for it in range(maxit):
        Ap = Cop(U, CT(U, p))
        pAp = (p * Ap).sum((0, 2, 3, 4)) + 1e-300
        al = (rr / pAp)[None, :, None, None, None]
        x += al * p
        r -= al * Ap
        rn = (r * r).sum((0, 2, 3, 4))
        if np.all(rn < tol * tol * b2):
            break
        p = r + (rn / rr)[None, :, None, None, None] * p
        rr = rn
    return e - CT(U, x)


def leapfrog(U, e, dtau, nsteps, F=None):
    if F is None:
        F = force(U)
    for _ in range(nsteps):
        e = e - 0.5 * dtau * F
        U = qmul(qexp(dtau * e), U)
        F = force(U)
        e = e - 0.5 * dtau * F
    U = U / np.sqrt((U * U).sum(0, keepdims=True))
    return U, e, F


def cool(U, depth, eps=0.2):
    """Gradient flow dU/ds = -grad V (left); flow time in a^2 = s/4."""
    n = max(1, int(round(depth / eps)))
    for _ in range(n):
        U = qmul(qexp(-eps * force(U)), U)
    return U / np.sqrt((U * U).sum(0, keepdims=True))


def clover_B(U):
    B = np.empty((3,) + U.shape[1:])
    for i in range(3):
        j, k = (i + 1) % 3, (i + 2) % 3
        Uj, Uk = U[:, :, j], U[:, :, k]
        P = qmul(qmul(Uj, sh(Uk, j, 1)), qmul(qdag(sh(Uj, k, 1)), qdag(Uk)))
        acc = P[1:].copy()
        Ujm = sh(Uj, j, -1); Ukm = sh(Uk, j, -1)
        P = qmul(qmul(Uk, qdag(sh(Ujm, k, 1))), qmul(qdag(Ukm), Ujm))
        acc += P[1:]
        Ujmk = sh(Ujm, k, -1); Ukmk = sh(Ukm, k, -1); Ukk = sh(Uk, k, -1)
        P = qmul(qmul(qdag(Ujm), qdag(Ukmk)), qmul(Ujmk, Ukk))
        acc += P[1:]
        Ujk = sh(Uj, k, -1)
        P = qmul(qmul(qdag(Ukk), Ujk), qmul(sh(Ukk, j, 1), qdag(Uj)))
        acc += P[1:]
        B[:, :, i] = 0.5 * acc            # (2 w) averaged over the 4 leaves
    return B


def link_B(U):
    B = clover_B(U)
    out = np.empty_like(B)
    for i in range(3):
        out[:, :, i] = 0.5 * (B[:, :, i] + rot(U[:, :, i], sh(B[:, :, i], i, 1)))
    return out


def dNcs_path(Ua, Ub, nsub=4):
    th = qlog(qmul(Ub, qdag(Ua)))
    tot = 0.0
    for s in range(nsub):
        Um = qmul(qexp((s + 0.5) / nsub * th), Ua)
        tot = tot + (th / nsub * link_B(Um)).sum((0, 2, 3, 4, 5))
    return tot / (8 * math.pi ** 2)


def leapfrog_naive(U, e, dtau, nsteps, F=None, Bl=None):
    """Leapfrog that integrates the naive dN_CS = dtau e.(B(U)+B(U'))/2 / 8pi^2 along the trajectory."""
    if F is None:
        F = force(U)
    if Bl is None:
        Bl = link_B(U)
    dN = 0.0
    for _ in range(nsteps):
        e = e - 0.5 * dtau * F
        Un = qmul(qexp(dtau * e), U)
        Bn = link_B(Un)
        dN = dN + (dtau * e * 0.5 * (Bl + Bn)).sum((0, 2, 3, 4, 5))
        U, Bl = Un, Bn
        F = force(U)
        e = e - 0.5 * dtau * F
    U = U / np.sqrt((U * U).sum(0, keepdims=True))
    return U, e, F, Bl, dN / (8 * math.pi ** 2)


def thermalize(B, N, beta, rng, n_ref=30, tau_ref=4.0, dtau=0.2):
    U = np.zeros((4, B, 3, N, N, N)); U[0] = 1.0
    F = None
    for r in range(n_ref):
        e = rng.normal(size=(3, B, 3, N, N, N)) / math.sqrt(beta)
        e = project_gauss(U, e)
        U, e, F = leapfrog(U, e, dtau, int(round(tau_ref / dtau)), F)
    return U, e, F


def run(N, beta, B=8, t_max=200.0, dt_meas=1.0, dtau=0.2, depths=(2.0, 4.0, 8.0), seed=1, n_ref=30, nsub=4,
        cool_eps=0.2):
    """t in lattice units t/a (tau = 2 t). Cooled N_CS at several flow depths (internal units; flow time in a^2 =
    depth/4) and the naive Kubo N_CS integrated along the real-time trajectory, every dt_meas."""
    rng = np.random.default_rng(seed)
    depths = sorted(depths)
    U, e, F = thermalize(B, N, beta, rng, n_ref=n_ref)
    E0 = energy(U, e)
    ekin = 0.5 * (e * e).sum((0, 2, 3, 4, 5)) / (3 * N ** 3)
    nmeas = int(round(t_max / dt_meas))
    nstep = int(round(2 * dt_meas / dtau))
    Nc = np.zeros((len(depths), nmeas + 1, B)); Nn = np.zeros((nmeas + 1, B))
    Ecool = np.zeros((len(depths), nmeas + 1, B))

    def cooled_set(U):
        out, Uc, done = [], U, 0.0
        for d in depths:
            Uc = cool(Uc, d - done, eps=cool_eps); done = d
            out.append(Uc)
        return out
    prev = cooled_set(U)
    for k in range(len(depths)):
        Ecool[k, 0] = plaq_energy(prev[k])
    Gmax, Bl = 0.0, None
    for m in range(1, nmeas + 1):
        U, e, F, Bl, dN = leapfrog_naive(U, e, dtau, nstep, F, Bl)
        Nn[m] = Nn[m - 1] + dN
        cur = cooled_set(U)
        for k in range(len(depths)):
            Nc[k, m] = Nc[k, m - 1] + dNcs_path(prev[k], cur[k], nsub=nsub)
            Ecool[k, m] = plaq_energy(cur[k])
        prev = cur
        if m % 100 == 0:
            Gmax = max(Gmax, float(np.abs(gauss(U, e)).max()))
    E1 = energy(U, e)
    return {"N": N, "beta": beta, "B": B, "t_max": t_max, "dt_meas": dt_meas, "dtau": dtau, "depths_internal": depths,
            "flow_time_a2": [d / 4 for d in depths],
            "Ncool": [Nc[k].T.tolist() for k in range(len(depths))], "Nnaive": Nn.T.tolist(),
            "Ecool_mean": [float(Ecool[k].mean()) for k in range(len(depths))],
            "energy_drift_rel": float(np.max(np.abs(E1 - E0) / np.abs(E0))),
            "ekin_per_dof_x_beta": float(ekin.mean() * beta), "plaq_per_p": float((plaq_energy(U) / (3 * N ** 3)).mean()),
            "gauss_max": Gmax}
