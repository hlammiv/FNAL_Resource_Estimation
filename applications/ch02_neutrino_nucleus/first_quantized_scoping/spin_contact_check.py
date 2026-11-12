"""Numerical checks behind the spin-dependent contact construction of README.md, section 3(i).

Claim. On antisymmetric states the chapter's on-site interaction (C_t in the spin-triplet pair channel,
C_s in the spin-singlet pair channel, D for every on-site triple) equals, in first quantization,

    V_full = sum_{i<j} delta(r_i,r_j) [Cbar + dC SWAP^sigma_ij] + D sum_{i<j<k} delta_ij delta_jk ,
    Cbar = (C_t + C_s)/2 ,  dC = (C_t - C_s)/2 ,

whose pair terms do not commute (App. E of arXiv:2507.22814), and ALSO equals

    V_excl = sum_i g(n_i) + dC sum_{i<j} e_ij SWAP^sigma_ij ,
    n_i = number of other particles on particle i's site,  g(n) = Cbar n/2 + D n(n-1)/6 ,
    e_ij = delta_ij [n_i = 1]   (the pair is alone on its site),

whose terms all commute, so that one potential layer exp(-i V_excl dt) is an exact product of gadgets.

Checks (small systems, dense tensors, one thread, a few seconds):
 1. V_full and V_excl agree on random antisymmetric states (eta=3 on 3 sites; eta=4 on 2 sites).
 2. The exclusive pair terms commute on the FULL space; the plain pair terms do not.
 3. The spectrum of V on antisymmetric states is the occupation-model spectrum {C_t, C_s, 3Cbar+D, 6Cbar+4D}.
 4. The chapter's 22 on-site Pauli terms (2033_model_study/triton_resources.json) give the same on-site energies.
 5. Two-qubit gadget: exp(-i th SWAP) = e^{-i th} B^dag diag(1,1,1,e^{2 i th}) B with B = (H x 1) CNOT.
 6. Species-labeled encoding (protons and neutrons in fixed registers, antisymmetrized within each species, no
    isospin qubit): the spectrum of V_excl equals the one-proton sector of the fully antisymmetrized spectrum.
 7. Walsh expansion of the hopping dispersion cos(2 pi p/8) on a 3-bit momentum register: number of terms.
"""
import itertools, json, os
import numpy as np

CT, CS, DD = -93.62604, -82.19588, 126.95409          # chapter Eq. for the 2033 model (app01, line 139)
CBAR, DC = (CT + CS) / 2, (CT - CS) / 2
rng = np.random.default_rng(7)


def shape(M, eta):
    return (M, 2, 2) * eta          # per particle: site, spin, isospin


def antisym(psi, eta):
    out = np.zeros_like(psi)
    for perm in itertools.permutations(range(eta)):
        sign = np.linalg.det(np.eye(eta)[list(perm)])
        axes = [3 * p + k for p in perm for k in range(3)]
        out += sign * np.transpose(psi, axes)
    return out


def delta(M, eta, i, j):
    """Boolean tensor delta(r_i, r_j), broadcastable against psi."""
    idx = np.indices((M,) * eta)
    d = (idx[i] == idx[j])
    return d.reshape(tuple(s for p in range(eta) for s in (M, 1, 1)))


def counts(M, eta):
    idx = np.indices((M,) * eta)
    n = [sum((idx[i] == idx[j]).astype(int) for j in range(eta) if j != i) for i in range(eta)]
    return [x.reshape(tuple(s for p in range(eta) for s in (M, 1, 1))) for x in n]


def swap_spin(psi, i, j):
    return np.swapaxes(psi, 3 * i + 1, 3 * j + 1)


def V_full(psi, M, eta):
    out = np.zeros_like(psi)
    for i, j in itertools.combinations(range(eta), 2):
        dij = delta(M, eta, i, j)
        out += dij * (CBAR * psi + DC * swap_spin(psi, i, j))
    for i, j, k in itertools.combinations(range(eta), 3):
        out += DD * (delta(M, eta, i, j) & delta(M, eta, j, k)) * psi
    return out


def excl_term(psi, M, eta, i, j, n):
    return (delta(M, eta, i, j) & (n[i] == 1)) * swap_spin(psi, i, j)


def plain_term(psi, M, eta, i, j):
    return delta(M, eta, i, j) * swap_spin(psi, i, j)


def V_excl(psi, M, eta):
    n = counts(M, eta)
    g = sum(CBAR * x / 2 + DD * x * (x - 1) / 6 for x in n)
    out = g * psi
    for i, j in itertools.combinations(range(eta), 2):
        out += DC * excl_term(psi, M, eta, i, j, n)
    return out


def check_equal(M, eta, trials=4):
    worst = 0.0
    for _ in range(trials):
        psi = antisym(rng.standard_normal(shape(M, eta)), eta)
        psi /= np.linalg.norm(psi)
        worst = max(worst, np.linalg.norm(V_full(psi, M, eta) - V_excl(psi, M, eta)))
    return worst


def check_commute(M, eta):
    n = counts(M, eta)
    psi = rng.standard_normal(shape(M, eta))
    psi /= np.linalg.norm(psi)
    pairs = list(itertools.combinations(range(eta), 2))
    we = wp = 0.0
    for a, b in itertools.combinations(pairs, 2):
        c1 = excl_term(excl_term(psi, M, eta, *b, n), M, eta, *a, n) - excl_term(excl_term(psi, M, eta, *a, n), M, eta, *b, n)
        c2 = plain_term(plain_term(psi, M, eta, *b), M, eta, *a) - plain_term(plain_term(psi, M, eta, *a), M, eta, *b)
        we, wp = max(we, np.linalg.norm(c1)), max(wp, np.linalg.norm(c2))
    return we, wp


def antisym_spectrum(M, eta):
    """Eigenvalues of V_excl restricted to the antisymmetric subspace (dense, small)."""
    dim = (4 * M) ** eta
    basis = []
    for _ in range(600):
        v = antisym(rng.standard_normal(shape(M, eta)), eta).reshape(dim)
        for b in basis:
            v -= b * (b @ v)
        nv = np.linalg.norm(v)
        if nv > 1e-8:
            v /= nv
            for b in basis:               # second pass
                v -= b * (b @ v)
            v /= np.linalg.norm(v)
            basis.append(v)
        if len(basis) >= comb(4 * M, eta):
            break
    B = np.array(basis)
    Vb = np.array([V_excl(b.reshape(shape(M, eta)), M, eta).reshape(dim) for b in B])
    return np.linalg.eigvalsh((B @ Vb.T + Vb @ B.T) / 2), len(basis)


def sector_spectrum(M, labeled):
    """eta = 3, one proton and two neutrons. labeled=True: isospin pattern fixed to (p, n, n), antisymmetric in the two
    neutrons only. labeled=False: fully antisymmetric states projected on the one-proton sector."""
    eta, dim = 3, (4 * M) ** 3
    basis = []
    target = 2 * M * comb(2 * M, 2)
    for _ in range(40 * target):
        psi = rng.standard_normal(shape(M, eta))
        if labeled:
            mask = np.zeros(shape(M, eta))
            mask[:, :, 0, :, :, 1, :, :, 1] = 1.0
            psi = psi * mask
            psi = psi - np.transpose(psi, (0, 1, 2, 6, 7, 8, 3, 4, 5))
        else:
            iso = np.indices(shape(M, eta))
            nprot = (iso[2] == 0).astype(int) + (iso[5] == 0) + (iso[8] == 0)
            psi = antisym(psi * (nprot == 1), eta)
        v = psi.reshape(dim)
        for _pass in range(2):
            for b in basis:
                v = v - b * (b @ v)
        nv = np.linalg.norm(v)
        if nv > 1e-8:
            basis.append(v / nv)
        if len(basis) >= target:
            break
    B = np.array(basis)
    Vb = np.array([V_excl(b.reshape(shape(M, eta)), M, eta).reshape(dim) for b in B])
    return np.sort(np.linalg.eigvalsh((B @ Vb.T + Vb @ B.T) / 2)), len(basis)


def comb(n, k):
    from math import comb as c
    return c(n, k)


def chapter_onsite():
    here = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(here, "..", "2033_model_study", "triton_resources.json")) as f:
        terms = json.load(f)["onsite_nonidentity_paulis"]
    P = dict(I=np.eye(2), X=np.array([[0, 1], [1, 0]]), Y=np.array([[0, -1j], [1j, 0]]), Z=np.diag([1., -1.]))
    H = np.zeros((16, 16), complex)
    N = np.zeros((16, 16))
    for q in range(4):
        ops = [np.eye(2)] * 4
        ops[q] = (np.eye(2) - P["Z"]) / 2
        N += np.kron(np.kron(ops[0], ops[1]), np.kron(ops[2], ops[3]))
    for s, c in terms:
        H += c * np.kron(np.kron(P[s[0]], P[s[1]]), np.kron(P[s[2]], P[s[3]]))
    w, v = np.linalg.eigh(H)
    nn = np.real(np.einsum("ij,jk,ki->i", v.conj().T, N, v))
    by_n = {}
    for e, k in zip(w, np.rint(nn).astype(int)):
        by_n.setdefault(int(k), []).append(float(e))
    e0 = by_n[0][0]                               # identity term is not in the list: shift so the empty site is 0
    return {k: sorted(round(e - e0, 5) for e in v_) for k, v_ in by_n.items()}


def gadget_error(th=0.37):
    SW = np.eye(4)[[0, 2, 1, 3]]
    Hd = np.array([[1, 1], [1, -1]]) / np.sqrt(2)
    CN = np.eye(4)[[0, 1, 3, 2]]                 # control = first qubit
    B = np.kron(Hd, np.eye(2)) @ CN
    w, v = np.linalg.eigh(SW)
    U = (v * np.exp(-1j * th * w)) @ v.T
    G = np.exp(-1j * th) * B.conj().T @ np.diag([1, 1, 1, np.exp(2j * th)]) @ B
    return float(np.abs(U - G).max())


def walsh_terms():
    p = np.arange(8)
    f = np.cos(2 * np.pi * p / 8)
    coeffs = {}
    for s in range(8):
        chi = np.array([(-1) ** bin(s & x).count("1") for x in p])
        c = float(f @ chi) / 8
        if abs(c) > 1e-12:
            coeffs[format(s, "03b")] = round(c, 6)
    return coeffs


if __name__ == "__main__":
    out = {}
    out["V_full_minus_V_excl_on_antisym_eta3_M3"] = check_equal(3, 3)
    out["V_full_minus_V_excl_on_antisym_eta4_M2"] = check_equal(2, 4)
    we, wp = check_commute(3, 3)
    out["max_commutator_exclusive_pair_terms_full_space"] = we
    out["max_commutator_plain_pair_terms_full_space"] = wp
    ev, nb = antisym_spectrum(2, 3)
    out["antisym_dim_eta3_M2"] = nb
    out["antisym_spectrum_eta3_M2_distinct"] = sorted(set(np.round(ev, 5)))
    out["expected_levels"] = dict(zero=0.0, C_t=CT, C_s=CS, triple=round(3 * CBAR + DD, 5), quad=round(6 * CBAR + 4 * DD, 5))
    out["chapter_onsite_energies_by_occupation"] = chapter_onsite()
    ea, na = sector_spectrum(2, labeled=False)
    el, nl = sector_spectrum(2, labeled=True)
    out["one_proton_sector_dim_full_vs_labeled"] = [na, nl]
    out["one_proton_sector_spectrum_max_difference"] = float(np.abs(ea - el).max()) if na == nl else None
    out["two_qubit_gadget_error"] = gadget_error()
    out["hopping_dispersion_walsh_terms_m3"] = walsh_terms()
    here = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(here, "spin_contact_check.json"), "w") as f:
        json.dump(out, f, indent=1)
    print(json.dumps(out, indent=1))
