"""Published T-count and qubit formulas of Spagnoli, Lissoni, Roggero, Quantum 10, 2200 (2026),
arXiv:2507.22814, in THAT PAPER'S OWN convention (4 T per Toffoli, T_ROT = 0.57 log2(1/eps) + 8.83,
phase-kickback kinetic term with q^2 dispersion, SU(4)-symmetric contacts C and G).

Nothing here is in the report's convention. The only purpose is to reproduce Table 5 of the paper,
so that the report-convention model in scoping.py starts from a checked implementation.

Equation numbers refer to the published version (v3).
"""
import math

HBARC = 197.3269804
# Table 4 of the paper (a = 1.4 fm): hbar^2/(2 mu a^2), C, G in MeV
KIN, C0, G0 = 10.58, -98.23, 127.84
D_DIM = 3


def log2(x):
    return math.log2(x)


# ---- Appendix A subroutines -------------------------------------------------------------------
def T_SQU(n):                       # Thm A.1, Eq. (142)
    return 4 * n * n - 4 * n


def T_QFT(n, eps):                  # Thm A.2, Eq. (143)
    s = 7 * n - 11
    for k in range(3, n):
        s += 8 * min(math.ceil(log2(n / eps)), k) - 15
    return s


def T_ROT(eps):                     # Thm A.3, Eq. (144)
    return 0.57 * log2(1 / eps) + 8.83


def T_MCX(n):                       # Thm A.4, Eq. (145)
    return 4 * (n - 1)


# ---- norms (Appendix C) -----------------------------------------------------------------------
def K_coeff(m):                     # Eq. (10)
    return KIN * (2 * math.pi / 2 ** m) ** 2


def lam_T(eta, m, d=D_DIM):         # Eq. (22): d K eta 2^(2m-2) = KIN d eta pi^2
    return d * K_coeff(m) * eta * 2 ** (2 * m - 2)


def normV(eta, C=C0, G=G0):         # Eqs. (203)-(204)
    a0 = abs(C) * (eta // 2)
    a1 = abs(3 * C + G) * (eta // 3) + (abs(C) if eta % 3 == 2 else 0)
    a2 = abs(6 * C + 4 * G) * (eta // 4)
    if eta % 4 == 2:
        a2 += abs(C)
    elif eta % 4 == 3:
        a2 += max(abs(C), abs(3 * C + G))
    return max(a0, a1, a2)


def lam_H(eta, m, d=D_DIM, C=C0, G=G0):   # Eq. (90)
    return eta * ((3 * abs(C) + 4 * G) / 2 + d * K_coeff(m) * 2 ** (2 * (m - 1)))


def alpha2(eta, m, d=D_DIM, C=C0, G=G0, tmax=None):   # Eq. (217)
    """tmax = max_i ||T_i||; the paper uses d K 2^(2m-2). A caller may pass another single-particle bound."""
    if tmax is None:
        tmax = d * K_coeff(m) * 2 ** (2 * m - 2)
    a = max(abs(2 * C), abs(4 * C + 3 * G), abs(6 * C + 9 * G))
    b = max(abs(C), abs(2 * C + G), abs(3 * C + 3 * G))
    return 4 * eta / 24 * tmax ** 2 * a + 4 * eta / 12 * tmax * b ** 2


# ---- times (Sec. 4) ---------------------------------------------------------------------------
def t_cross(m, a_fm=1.4, mu=939.0, E=10.0):          # Eq. (137)
    return a_fm * 2 ** m / HBARC * math.sqrt(mu / (2 * E))


def t_response(eta, m, domega=100.0):                # Eq. (138)
    dH = lam_T(eta, m) + normV(eta) + 18 * eta
    return (dH / domega - 1) * 2 * math.pi / dH


# ---- Lemma 1, Lemma 3 -------------------------------------------------------------------------
def b_diag_pf2(eta, m, t, eps, d=D_DIM):             # Eq. (224)
    x = 24 * d * eta * t * lam_T(eta, m, d) / eps
    return math.ceil(log2(x * log2(x)))


def b_qft_pf2(eta, m, r, eps, d=D_DIM):              # Eq. (225)
    return min(m - 1, math.ceil(log2(12 * d * eta * m * r / eps))) + 1


def T_kin(eta, m, b_diag, eps_qft, d=D_DIM, variant="eq19"):
    """Eq. (14). variant="eq19": T_DIAG = 8 b (2m-1)/2 as printed in Eq. (19) (the min taken at its first argument).
    variant="lemma6_wH1": T_DIAG = 4 b w_H with w_H = 1, the smallest value Lemma 6, Eq. (161) allows."""
    t_diag = 8 * b_diag * (2 * m - 1) / 2 if variant == "eq19" else 4 * b_diag
    return d * eta * (2 * T_SQU(m - 1) + t_diag + 2 * T_QFT(m, eps_qft))


def Gamma(eta):                                      # Eq. (30)
    return eta * (eta - 1) * (eta - 2) // 6 + eta * (eta - 1) // 2


def T_pot(eta, m, eps, d=D_DIM):                     # Eq. (29)
    g = Gamma(eta)
    return 4 * eta * (eta - 1) * (eta - 2) / 3 + g * (2 * T_MCX(d * m) + T_ROT(eps / g))


def pf2(eta, m, t, eps, d=D_DIM, variant="eq19"):
    """Theorem 2 / Lemma 7: second-order product formula. Returns dict with T count, steps, qubits."""
    r = math.ceil(math.sqrt(4 * t ** 3 * alpha2(eta, m, d) / eps))
    bd = b_diag_pf2(eta, m, t, eps, d)
    bq = b_qft_pf2(eta, m, r, eps, d)
    tv = T_pot(eta, m, eps / (4 * r), d)
    tt = T_kin(eta, m, bd, eps / (4 * r) / (3 * d * eta), d, variant)
    total = r * tv + (r + 1) * tt
    anc = bq + bd + max(2 * bq - 1, bd, m * (m - 1), d * m - 1)        # Eq. (223)
    return dict(T=total, T_potential_only=r * tv, steps=r, T_V=tv, T_T=tt, b_diag=bd, b_qft=bq, ancilla=anc,
                qubits_dm=eta * d * m + anc, qubits_eq5=eta * (d * m + 2) + anc)


# ---- Theorem 3-5: block encoding, QSP, GQSP ---------------------------------------------------
def b_r_qsp(eta, m, t, eps):                         # Eq. (110)
    return math.ceil(0.5 * log2(72 * math.pi ** 2 * lam_H(eta, m) * abs(t) / eps))


def b_qft_qsp(eta, m, t, eps, d=D_DIM):              # Eq. (110)
    return min(m - 1, math.ceil(log2(16 * d * lam_H(eta, m) * abs(t) * m / eps))) + 1


def T_prep(eta, m, eps):                             # Eq. (89), b_r of Eq. (92)
    br = math.ceil(0.5 * log2(18 * math.pi ** 2 / eps))
    return 12 * math.ceil(log2(eta)) + 16 * br + 8 * m - 64 + 4 * T_ROT(eps / 32)


def T_sel(eta, m, eps, d=D_DIM):                     # Eq. (89)
    return 24 * (eta - 1) * d * m + 4 * (d + 2) * (m - 1) + 16 + 2 * d * T_QFT(m, eps / (4 * d)) - 8 * eta + 52


def Q_qsp(eta, m, t, eps_t):                         # Eq. (108)
    return math.ceil(2 * lam_H(eta, m) * abs(t) + 3 * math.log(6 / eps_t))


def qsp_qubits(eta, m, t, eps, d=D_DIM, extra=12):   # Eq. (109) [extra=12] / Eq. (132) [extra=13]
    br, bq = b_r_qsp(eta, m, t, eps), b_qft_qsp(eta, m, t, eps, d)
    anc = eta + (2 + d) * m + extra + br + bq + max(2 * bq - 1, d * m + 6, eta + 2 * m + 9)
    return dict(ancilla=anc, qubits_dm=eta * d * m + anc, qubits_eq5=eta * (d * m + 2) + anc, b_r=br, b_qft=bq)


def qsp(eta, m, t, eps, d=D_DIM):                    # Eq. (107)
    Q = Q_qsp(eta, m, t, eps / 2)
    e = eps / (4 * lam_H(eta, m) * abs(t))
    T = 2 * Q * (T_sel(eta, m, e, d) + 2 * T_prep(eta, m, e) + 3) + T_prep(eta, m, e) \
        + 2 * (Q + 1) * T_ROT(eps / (4 * (Q + 1)))
    return dict(T=T, Q=Q, **qsp_qubits(eta, m, t, eps, d, 12))


def gqsp(eta, m, t, eps, d=D_DIM):                   # Eq. (130)
    Q = Q_qsp(eta, m, t, eps / 2)
    e = eps / (4 * lam_H(eta, m) * abs(t))
    T = Q * T_sel(eta, m, e, d) + 2 * (Q + 1) * T_prep(eta, m, e) \
        + 3 * (Q + 1) * T_ROT(eps / (6 * (Q + 1))) + 2 * Q * T_MCX(eta + 2 * m + 10)
    return dict(T=T, Q=Q, T_sel=T_sel(eta, m, e, d), T_prep=T_prep(eta, m, e),
                **qsp_qubits(eta, m, t, eps, d, 13))


# ---- Table 5 ----------------------------------------------------------------------------------
TABLE5 = {  # (method, eta): (T[t_cross], qubits[t_cross], T[t_dw], qubits[t_dw])
    ("PF2", 16): (1.3e9, 196, 8.2e7, 190), ("PF2", 40): (3.3e10, 416, 2.1e9, 410),
    ("QSP", 16): (7e7, 235, 1.12e7, 234), ("QSP", 40): (4e8, 500, 6.4e7, 498),
    ("GQSP", 16): (3.75e7, 235, 5.99e6, 234), ("GQSP", 40): (2.11e8, 500, 3.38e7, 498),
}


def check_table5(m=3, eps=0.1):
    rows = []
    for (meth, eta), pub in TABLE5.items():
        f = dict(PF2=pf2, QSP=qsp, GQSP=gqsp)[meth]
        tc, tr = t_cross(m), t_response(eta, m)
        a, b = f(eta, m, tc, eps), f(eta, m, tr, eps)
        extra = {}
        if meth == "PF2":
            a2, b2 = pf2(eta, m, tc, eps, variant="lemma6_wH1"), pf2(eta, m, tr, eps, variant="lemma6_wH1")
            extra = dict(ratio_cross_wH1=a2["T"] / pub[0], ratio_resp_wH1=b2["T"] / pub[2],
                         ratio_cross_potential_only=a["T_potential_only"] / pub[0],
                         ratio_resp_potential_only=b["T_potential_only"] / pub[2])
        rows.append(dict(method=meth, eta=eta, t_cross=tc, t_response=tr, **extra,
                         T_cross=a["T"], T_cross_pub=pub[0], ratio_cross=a["T"] / pub[0],
                         q_cross=a["qubits_dm"], q_cross_pub=pub[1],
                         T_resp=b["T"], T_resp_pub=pub[2], ratio_resp=b["T"] / pub[2],
                         q_resp=b["qubits_dm"], q_resp_pub=pub[3],
                         q_resp_eq5=b["qubits_eq5"]))
    return rows


if __name__ == "__main__":
    print(f"{'method':5s} {'eta':>3s} | {'T[t_cross]':>10s} {'pub':>8s} {'ratio':>5s} {'Q':>4s} {'pub':>4s} |"
          f" {'T[t_dw]':>10s} {'pub':>8s} {'ratio':>5s} {'Q':>4s} {'pub':>4s}")
    for r in check_table5():
        print(f"{r['method']:5s} {r['eta']:3d} | {r['T_cross']:10.3e} {r['T_cross_pub']:8.2e} {r['ratio_cross']:5.2f} "
              f"{r['q_cross']:4d} {r['q_cross_pub']:4d} | {r['T_resp']:10.3e} {r['T_resp_pub']:8.2e} "
              f"{r['ratio_resp']:5.2f} {r['q_resp']:4d} {r['q_resp_pub']:4d}")
    for r in check_table5():
        if r["method"] == "PF2":
            print(f"PF2 eta={r['eta']}: ratio to published with T_DIAG=4b (w_H=1): t_cross {r['ratio_cross_wH1']:.3f}, "
                  f"t_dw {r['ratio_resp_wH1']:.3f}; potential term alone: {r['ratio_cross_potential_only']:.3f}, "
                  f"{r['ratio_resp_potential_only']:.3f}")
    print("t_cross =", t_cross(3), " t_response(16), (40) =", t_response(16, 3), t_response(40, 3))
