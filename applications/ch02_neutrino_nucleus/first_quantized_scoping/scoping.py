"""First-quantized A >= 12 scoping for the Ch. 2 2033 benchmark, in the report's conventions.

Run from anywhere:  python3 scoping.py            (writes results.json, prints the tables of README.md)

Conventions (imported from estimates/common.py, never copied): 7 T per Toffoli; rotations by randomized
repeat-until-success synthesis at eps_syn = 1e-2 per shot, eps_rot = sqrt(eps_syn / N_rot),
c_T = 1.15 log2(1/eps_rot) + 9.2; 0.1 expected faults per shot; 1 us per T plus 0.1 ms per shot; ONE machine,
serial shots; 2033 targets 1000 LQ and 1e9 T; five-year window.

Every input carries a provenance tag in INPUTS: Cited / Computed / Assumed. See README.md section 6.
"""
import json, math, os, sys
from dataclasses import dataclass

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "..", "..")))
from estimates.common import (t_per_rotation, eps_rot_for, toffoli_t, EPS_SYN, RFI_LQ_2033, RFI_T_2033,  # noqa: E402
                              SECONDS_PER_YEAR)
import spagnoli_lemmas as sp  # noqa: E402

# ----------------------------------------------------------------------------------------------
# Inputs
# ----------------------------------------------------------------------------------------------
T_HOP, C_T, C_S, D3 = 10.5794, -93.62604, -82.19588, 126.95409      # chapter 2033 model (MeV)
CBAR, DC = (C_T + C_S) / 2, (C_T - C_S) / 2
A_FM, HBARC = 1.4, 197.3269804
TAUS = (0.05, 0.10, 0.20)                                            # chapter time grid (MeV^-1)
T_GATE_S, T0_S, FAULTS = 1e-6, 1e-4, 0.1                             # report time model and fault budget
WINDOW_YR = 5.0
EPS_EVOL = 0.00731       # Assumed: evolution error target at tau = 0.20, set equal to the chapter's measured triton bias
INFID = 1e-3             # Assumed: residual excited-state weight after the filter
TRITON_REF = dict(steps=100, tau=0.20, bias=0.00731)                 # chapter: measured second-order bias
SHOTS_FIRST = 2 * math.ceil(2 * math.log(4 / .05) / .05 ** 2)        # 7,012  (chapter Eq. Nshot, m=2)
SHOTS_PER_EST = math.ceil(2 * math.log(12 / .05) / .02 ** 2)         # 27,404 (m=6)
CURRENT = dict(lq=908, t_shot=2.07e8, first_d=9.4, campaign_d=250.0)

NUCLEI = {  # A, Z, gap proxy (MeV), overlap^2
    "3H":   dict(A=3, Z=1, gap=None, gamma2=1.0),
    "12C":  dict(A=12, Z=6, gap=4.44, gamma2=0.1),
    "16O":  dict(A=16, Z=8, gap=6.05, gamma2=0.1),
    "40Ar": dict(A=40, Z=18, gap=1.46, gamma2=0.1),
}

INPUTS = [
    # name, value, tag, source, note / plausible range
    ("t_hop, C_t, C_s, D", f"{T_HOP}, {C_T}, {C_S}, {D3} MeV", "Cited", "app01_neutrino_nucleus.tex 2033 model; TRITON.md",
     "Chapter couplings, fitted for A<=3 in the 6^3 box. Not refitted for A>=12 or for 8^3."),
    ("lattice spacing", "1.4 fm", "Cited", "app01; arXiv:1911.06368", ""),
    ("time grid", "0.05, 0.10, 0.20 MeV^-1", "Cited", "app01 2033 box", "Kept. q = 2 pi/(8a) = 110.7 MeV on 8^3."),
    ("shots", f"{SHOTS_FIRST} first result; {SHOTS_PER_EST} x 6 campaign", "Cited", "app01 Eq. (Nshot)",
     "Hoeffding plus union bound; +-0.05 (2 estimates) and +-0.02 (6 estimates), 95% simultaneous."),
    ("system register", "(3m+2) A qubits, m=3", "Cited", "arXiv:2507.22814 Eq. (5)", "Labeled triton: (3m+1) A = 30."),
    ("Lemma 1, Lemma 3, Thm 2-5, App. A, Eqs. (90), (108), (217)", "formulas", "Cited", "arXiv:2507.22814",
     "Implemented in spagnoli_lemmas.py; Table 5 reproduced (README section 2)."),
    ("7 T/Toffoli, RUS synthesis, eps_syn=1e-2", "common.py", "Cited", "main-overview.tex conventions; estimates/common.py", ""),
    ("1 us per T, 0.1 ms per shot, 0.1 faults per shot", "report", "Cited", "main-overview.tex", "Serial T count on one machine; no T-depth derived."),
    ("spin-dependent contact, exclusive-pair form", "V_excl", "Computed", "spin_contact_check.py",
     "Equal to the chapter's on-site interaction on antisymmetric states; all terms commute. Circuit is a construction-level count, not compiled."),
    ("potential layer counts", "P(6(dm-1)+14) Toffoli, P+6A rotations", "Computed", "scoping.py potential_layer",
     "Three equality sweeps; compute and uncompute both counted at 7 T per Toffoli."),
    ("kinetic layer counts (8^3)", "14 rotations + 12 T per coordinate", "Computed", "scoping.py kinetic_layer; Walsh count in spin_contact_check.py",
     "Exact QFT_8 (2 CS + 1 CT each way) and 4 controlled Walsh phases for the hopping dispersion."),
    ("kinetic layer counts (6^3)", "<= 40 rotations + 16 Toffoli per coordinate", "Computed (bound)", "scoping.py kinetic_layer",
     "CRT relabeling, generic two-qubit QFT_3; no circuit constructed."),
    ("generalized alpha_2 for C_t != C_s", "Eq. (217) with C -> Cbar and +|dC| terms", "Computed", "scoping.py alpha2_chapter",
     "Our extension of the paper's semi-norm bound; hopping dispersion uses max single-particle kinetic energy 12 t."),
    ("block-encoding norm", "lambda = A(6t + (3|Cbar|+4D)/2) + |dC| A(A-1)/2", "Computed (bound)", "scoping.py lam_chapter",
     "Paper Eq. (90) form with hopping kinetic LCU and a pair-SWAP LCU. The paper leaves the SWAP block encoding to future work."),
    ("spin SELECT cost", "2(A-1)(dm+2)+2(dm-1)+5 Toffoli plus a second PREPARE", "Computed (bound)", "scoping.py walk_step",
     "Second swap network and comparator; construction-level, self-inverse property not verified."),
    ("Slater-determinant preparation", "A^2 orbital unitaries, A(A-1)/2 swap+C^kX, A^2-4A+3 rotations", "Cited",
     "arXiv:2509.07279 Sec. 2-3", "Each orbital unitary costed here as a 511-rotation real loader on 9 spatial qubits (Computed)."),
    ("eigenstate filter degree", "2 l, l = ln(2/eps_f) lambda sqrt(2)/gap", "Cited", "arXiv:1910.14596 Lemma 2",
     "Needs E_0 known to much better than the gap. That knowledge is Assumed."),
    ("evolution error target", f"{EPS_EVOL}", "Assumed", "this study", "Half the +-0.02 tolerance; chapter measured 0.0073. T scales as eps^-1/2 (PF) or log (QSP)."),
    ("PF step count, calibrated", "100 sqrt(A/3) (tau/0.2)^1.5 sqrt(0.00731/eps)", "Assumed", "chapter triton bias, extensive scaling",
     "No A>=12 numerics exist. Range: this value up to the rigorous bound (60-110 times more steps)."),
    ("PF step count, bound", "sqrt(4 tau^3 alpha_2/eps)", "Cited", "arXiv:2507.22814 Lemma 7", "Worst-case semi-norm; 100 times above the triton's measured need."),
    ("ground-state gap", "4.44 (12C), 6.05 (16O), 1.46 (40Ar) MeV", "Assumed", "physical first excited states as proxies",
     "The model's finite-volume gap is unknown. Range 1-8 MeV; LO pionless nuclei may sit near cluster thresholds."),
    ("mean-field overlap squared", "0.1", "Assumed", "no source", "Range 1e-3 to 0.5. Sets attempts per accepted shot."),
    ("filter infidelity", f"{INFID}", "Assumed", "this study", "Enters logarithmically."),
    ("rodeo filter time", "N_cyc x 0.8 x 3/gap", "Assumed", "arXiv:2009.04092 cycle factor [1+exp(-gap^2 sigma^2/2)]/2, sigma=3/gap",
     "Product-formula route only. Constants from the published cycle factor; schedule not optimized."),
    ("charge-operator LCU norm", "(A/Z)^2 on the correlator, (A/Z)^4 on shots", "Computed", "scoping.py shot_factor",
     "Antisymmetrized registers carry no proton label. Isoscalar density has factor 1."),
    ("qubitized-route extra qubits", "A + 6", "Assumed", "this study", "Second unary register and flags; not laid out."),
    ("current triton numbers", "908 LQ, 2.07e8 T, 9.4 d, 250 d", "Cited", "app01 2033 box", ""),
]


# ----------------------------------------------------------------------------------------------
# Gate bookkeeping in the report convention
# ----------------------------------------------------------------------------------------------
@dataclass
class G:
    rot: float = 0.0     # synthesized rotations
    tof: float = 0.0     # Toffolis
    fix: float = 0.0     # exact T gates

    def __add__(s, o):
        return G(s.rot + o.rot, s.tof + o.tof, s.fix + o.fix)

    def __mul__(s, k):
        return G(s.rot * k, s.tof * k, s.fix * k)
    __rmul__ = __mul__


def c_T(n_rot):
    return t_per_rotation(eps_rot_for(max(n_rot, 1.0), EPS_SYN), "rus")


def T_of(g, n_rot_circuit):
    """T count of gate bundle g inside a circuit whose total rotation count is n_rot_circuit."""
    return g.rot * c_T(n_rot_circuit) + toffoli_t(g.tof, "textbook") + g.fix


# ----------------------------------------------------------------------------------------------
# Product-formula route: exact kinetic layer, exact (exclusive-pair) potential layer
# ----------------------------------------------------------------------------------------------
def kinetic_layer(A, L=8):
    """Controlled exp(-i K dt) for the chapter's hopping dispersion, per layer.
    L=8: QFT_8 and inverse per coordinate (each 2 controlled-S = 6 T and one controlled-T = 3 rotations),
         4 Walsh phases, each controlled by the Hadamard qubit = 2 rotations.
    L=6: bound, see INPUTS."""
    ncoord = 3 * A
    if L == 8:
        return G(rot=ncoord * (2 * 3 + 2 * 4), fix=ncoord * 12)
    if L == 6:
        return G(rot=ncoord * (2 * 15 + 2 * 5), tof=ncoord * 16)
    raise ValueError(L)


def potential_layer(A, dm=9, store_flags=False, labeled_triton=False):
    """Controlled exp(-i V dt), V = sum_i g(n_i) + dC sum_{i<j} e_ij SWAP_ij (spin_contact_check.py).
    Per pair, one equality flag costs (dm-1) Toffolis to compute and the same to uncompute (as in Lemma 3)."""
    P = A * (A - 1) // 2
    if labeled_triton:
        # three flags kept; diagonal: 3 controlled pair phases + 1 triple phase; spin: 3 exclusive-pair phases
        return G(rot=6 + 1 + 3, tof=3 * 2 * (dm - 1) + 2 + 3 * 4)
    sweeps = 1 if store_flags else 3
    tof = P * (sweeps * 2 * (dm - 1)      # equality flags, compute + uncompute
               + 2 * 2                    # counter increments (compute sweep and uncompute sweep, 2 counters)
               + 2 * 5)                   # AND of (flag, n_lo, not n_hi, s_i, s_j, control), compute + uncompute
    rot = P + 6 * A                       # one singlet phase per pair; 3 Walsh phases per counter, controlled
    return G(rot=rot, tof=tof)


def pf_evolution(A, n_steps, L=8, **kw):
    return n_steps * kinetic_layer(A, L) + (n_steps + 1) * potential_layer(A, **kw) + G(rot=1)


def alpha2_chapter(A, tmax=12 * T_HOP):
    """Eq. (217) of arXiv:2507.22814 generalized to C_t != C_s (README 3(i)); tmax = max single-particle kinetic energy."""
    a = max(abs(2 * CBAR) + 2 * abs(DC), abs(4 * CBAR + 3 * D3) + 2 * abs(DC), abs(6 * CBAR + 9 * D3))
    b = max(abs(CBAR) + abs(DC), abs(2 * CBAR + D3) + abs(DC), abs(3 * CBAR + 3 * D3))
    return 4 * A / 24 * tmax ** 2 * a + 4 * A / 12 * tmax * b ** 2


def dt_bound(A, eps=EPS_EVOL, tau=0.20):
    """Step from the rigorous second-order bound (Lemma 7 of arXiv:2507.22814) at the longest time."""
    return tau / math.ceil(math.sqrt(4 * tau ** 3 * alpha2_chapter(A) / eps))


def dt_calibrated(A, eps=EPS_EVOL):
    """ASSUMED: the triton's measured bias (0.00731 at 100 steps of 0.002) with an error that grows linearly in A
    at fixed step, so dt = 0.002 sqrt(3/A) sqrt(eps/0.00731). Fixed step for all times, as in the chapter."""
    r = TRITON_REF
    return r["tau"] / r["steps"] * math.sqrt(3 / A * eps / r["bias"])


def steps_bound(A, tau, eps=EPS_EVOL):
    return math.ceil(tau / dt_bound(A, eps) - 1e-9)


def steps_calibrated(A, tau, eps=EPS_EVOL):
    return math.ceil(tau / dt_calibrated(A, eps) - 1e-9)


RODEO_S = (1 + math.exp(-4.5)) / 2                    # per-cycle survival of an excited state at sigma = 3/gap


def rodeo(gap, gamma2, infid=INFID):
    """Cycles and total controlled-evolution time (MeV^-1) of a rodeo filter, arXiv:2009.04092.
    Gaussian times of rms sigma = 3/gap, mean |t| = 0.8 sigma. Residual excited weight (1-g2) s^n / g2 <= infid."""
    n = max(1, math.ceil(math.log(gamma2 * infid / (1 - gamma2)) / math.log(RODEO_S)))
    return n, n * 0.8 * 3 / gap


def rodeo_cycles_before_abort(n):
    """Expected cycles run by an attempt that starts in an excited state: it is rejected at the first '1' outcome."""
    s = RODEO_S
    return sum(k * s ** (k - 1) * (1 - s) for k in range(1, n + 1)) + n * s ** n


# ----------------------------------------------------------------------------------------------
# Qubitization route
# ----------------------------------------------------------------------------------------------
def lam_chapter(A, dispersion="hopping", spin=True):
    lt = 6 * T_HOP * A if dispersion == "hopping" else sp.lam_T(A, 3)
    return lt + A * (3 * abs(CBAR) + 4 * D3) / 2 + (abs(DC) * A * (A - 1) / 2 if spin else 0.0)


def walk_step(A, m=3, d=3, spin=True, eps_be=1e-6):
    """One controlled walk step, as Toffolis and rotations. Published SELECT/PREPARE T counts (4 T per Toffoli)
    are converted to Toffolis by /4; their rotations are synthesized in the report convention."""
    dm = d * m
    sel = sp.T_sel(A, m, eps_be, d) / 4
    prep_tof = (12 * math.ceil(math.log2(A)) + 16 * math.ceil(0.5 * math.log2(18 * math.pi ** 2 / eps_be)) + 8 * m - 64) / 4
    prep_tof = max(prep_tof, 0)
    n_prep = 2                                            # PREPARE and its inverse per step
    refl = A + 2 * m + 10 - 1
    tof = sel + n_prep * prep_tof + refl + 4              # +4: Hadamard-test control (analogy with the paper's +3)
    rot = n_prep * 4 + 3                                  # 4 PREPARE rotations each way, 3 signal rotations (GQSP)
    if spin:
        tof += 2 * (A - 1) * (dm + 2) + 2 * (dm - 1) + 5  # second swap network, comparator, controlled spin SWAP
        tof += n_prep * prep_tof + (A + 4)                # second unary PREPARE; wider reflection
        rot += n_prep * 3
    return G(rot=rot, tof=tof)


def Q_evol(lam, tau, eps=EPS_EVOL):
    return math.ceil(2 * lam * tau + 3 * math.log(6 / (eps / 2)))     # Eq. (108) with eps_t = eps/2


def filter_degree(lam, gap, gamma2, infid=INFID):
    eps_f = math.sqrt(infid * gamma2 / (1 - gamma2))
    ell = math.ceil(math.log(2 / eps_f) / (math.sqrt(2) * gap / (2 * lam)))    # arXiv:1910.14596 Lemma 2(ii)
    return 2 * ell


def qubitized_ancilla(A, lam, tmax, m=3, d=3):
    br = math.ceil(0.5 * math.log2(72 * math.pi ** 2 * lam * tmax / EPS_EVOL))
    bq = min(m - 1, math.ceil(math.log2(16 * d * lam * tmax * m / EPS_EVOL))) + 1
    return A + (2 + d) * m + 13 + br + bq + max(2 * bq - 1, d * m + 6, A + 2 * m + 9)    # Eq. (132)


# ----------------------------------------------------------------------------------------------
# Preparation
# ----------------------------------------------------------------------------------------------
def slater(groups, k_spatial=9, k=10):
    """arXiv:2509.07279, for each separately antisymmetrized group of n particles: n(n+1)/2 uses of U_j and n(n-1)/2 of
    U_j^dagger; n(n-1)/2 x (k cswap + C^kX); n^2-4n+3 rotations for the Y states.
    U_j: real-amplitude loader on the spatial register, 2^k_spatial - 1 rotations (spin set by X gates)."""
    g = G()
    for n in groups:
        g = g + G(rot=n * n * (2 ** k_spatial - 1) + max(n * n - 4 * n + 3, 0), tof=n * (n - 1) / 2 * (2 * k - 1))
    return g


def triton_loader(L=8):
    """Chapter construction without the occupation write: 2^20-1 rotation loader, centre of mass, translations, spin decode."""
    if L == 8:
        return G(rot=2 ** 20 - 1, tof=6 * 2 + 32)
    return G(rot=2 ** 20 - 1 + 21, tof=6 * 8 + 32)


def insertion():
    """Controlled e^{i q x} on one coordinate register, q = 2 pi/L: controlled T, S, Z. Term sampled classically."""
    return G(rot=3, fix=3)


def shot_factor(A, Z, encoding):
    """Shots multiplier for the proton-charge correlator normalized as in the chapter.
    labeled: protons occupy registers 0..Z-1, J = (1/Z) sum_p exp(i q x_p) has unit LCU norm -> 1.
    eq5: fully antisymmetrized register with isospin qubits (Eq. (5) of arXiv:2507.22814). J needs the projector
         (1+tau_z)/2 on every particle, LCU norm A/Z per insertion, (A/Z)^2 on the correlator, (A/Z)^4 on shots."""
    return 1.0 if encoding == "labeled" else (A / Z) ** 4


# ----------------------------------------------------------------------------------------------
# One instance
# ----------------------------------------------------------------------------------------------
def instance(nuc, route, L=8, steps="calibrated", encoding="labeled", gap=None, gamma2=None, store_flags=False,
             dispersion="hopping"):
    p = NUCLEI[nuc]
    A, Z = p["A"], p["Z"]
    gap = p["gap"] if gap is None else gap
    gamma2 = p["gamma2"] if gamma2 is None else gamma2
    triton = (nuc == "3H")
    dm = 9
    out = dict(nucleus=nuc, A=A, lattice=f"{L}^3", route=route, steps_mode=steps, encoding=encoding,
               gap_MeV=gap, gamma2=gamma2)
    spin_bits = 1 if encoding == "labeled" else 2
    stepf = steps_calibrated if steps == "calibrated" else steps_bound

    # --- preparation ---------------------------------------------------------------------------
    if triton:
        prep_det, filt, attempts = triton_loader(L), G(), 1.0
    else:
        prep_det = slater((Z, A - Z), k=dm + 1) if encoding == "labeled" else slater((A,), k=dm + 2)
        attempts = 1.0 / gamma2
    fail_g = None            # gates spent by one rejected attempt

    # --- evolution and filter ------------------------------------------------------------------
    evol = {}
    if route == "PF2":
        kw = dict(labeled_triton=triton, store_flags=store_flags)
        for tau in TAUS:
            evol[tau] = pf_evolution(A, stepf(A, tau), L, **kw)
        out["steps"] = {str(t): stepf(A, t) for t in TAUS}
        g1 = pf_evolution(A, 1, L, **kw) + (-1) * pf_evolution(A, 0, L, **kw)
        out["per_step"] = dict(rot=g1.rot, toffoli=g1.tof, fixed_T=g1.fix)
        if not triton:
            ncyc, tf = rodeo(gap, gamma2)
            dt = TAUS[-1] / stepf(A, TAUS[-1])
            nf = math.ceil(tf / dt)
            filt = pf_evolution(A, nf, L, **kw) + G(rot=ncyc)
            nab = rodeo_cycles_before_abort(ncyc)
            fail_g = prep_det + (nab / ncyc) * filt
            out["filter"] = dict(kind="rodeo, rejected attempts stop at the first failed cycle", cycles=ncyc,
                                 time_MeV_inv=tf, steps=nf, cycles_before_abort=nab)
        sys_q = (dm + spin_bits) * A
        anc = (3 + 5 + (dm - 1) + 2) if triton else (2 * A + (dm - 1) + 5 + 2 + 1)
        if store_flags and not triton:
            anc += A * (A - 1) // 2
        out["register"] = dict(system=sys_q, counters=0 if triton else 2 * A, flag_ladder=dm - 1, and_chain=5,
                               pair_flags=(3 if triton else (A * (A - 1) // 2 if store_flags else 0)),
                               rus_hadamard_filter=2 if triton else 3, total=sys_q + anc)
    elif route == "GQSP":
        lam = lam_chapter(A, dispersion)
        w = walk_step(A)
        for tau in TAUS:
            Q = Q_evol(lam, tau)
            evol[tau] = Q * w + G(rot=3)
        out["lambda_MeV"] = lam
        out["walk_queries"] = {str(t): Q_evol(lam, t) for t in TAUS}
        out["per_step"] = dict(rot=w.rot, toffoli=w.tof, fixed_T=0)
        tmax = TAUS[-1]
        if not triton:
            deg = filter_degree(lam, gap, gamma2)
            filt = deg * w + G(rot=deg + 1)
            fail_g = prep_det + filt
            out["filter"] = dict(kind="eigenstate filter (QSVT), one flag at the end, no early stop", degree=deg)
            tmax = max(tmax, deg / (2 * lam))
        sys_q = (dm + spin_bits) * A
        anc = qubitized_ancilla(A, lam, tmax) + A + 6
        out["register"] = dict(system=sys_q, block_encoding_and_phase=anc - A - 6, second_selection_and_flags=A + 6,
                               total=sys_q + anc)
    else:
        raise ValueError(route)

    # --- totals --------------------------------------------------------------------------------
    deepest = prep_det + filt + insertion() + evol[TAUS[-1]]
    n_rot = deepest.rot
    T = lambda g: T_of(g, n_rot)          # one synthesis tolerance, set by the deepest circuit  # noqa: E731
    out["c_T"] = c_T(n_rot)
    out["rotations_deepest"] = n_rot
    out["toffolis_deepest"] = deepest.tof
    out["T_prep_deterministic"] = T(prep_det)
    out["T_filter_per_attempt"] = T(filt)
    out["T_evolution"] = {str(t): T(evol[t]) for t in TAUS}
    out["T_per_step"] = T(G(**{k: v for k, v in zip(("rot", "tof", "fix"), (out["per_step"]["rot"], out["per_step"]["toffoli"], out["per_step"]["fixed_T"]))}))
    out["T_deepest_circuit"] = T(deepest)
    out["attempts_per_accepted_shot"] = attempts
    T_fail = 0.0 if fail_g is None else T(fail_g)
    out["T_rejected_attempt"] = T_fail
    exp_T = {t: (attempts - 1) * T_fail + T(prep_det + filt + insertion() + evol[t]) for t in TAUS}
    out["T_expected_per_accepted_shot"] = {str(t): exp_T[t] for t in TAUS}
    F = 1.0 if triton else shot_factor(A, Z, encoding)
    out["shot_factor"] = F
    n_first = SHOTS_FIRST * F
    n_camp = 6 * SHOTS_PER_EST * F
    T_first = n_first * exp_T[0.10]
    T_camp = 2 * SHOTS_PER_EST * F * sum(exp_T.values())
    out["shots_first"], out["shots_campaign"] = n_first, n_camp
    out["T_first_result"], out["T_campaign"] = T_first, T_camp
    out["wall_shot_s"] = out["T_deepest_circuit"] * T_GATE_S + T0_S
    out["wall_first_d"] = (T_first * T_GATE_S + n_first * attempts * T0_S) / 86400
    out["wall_campaign_d"] = (T_camp * T_GATE_S + n_camp * attempts * T0_S) / 86400
    out["wall_campaign_yr"] = out["wall_campaign_d"] * 86400 / SECONDS_PER_YEAR
    out["eps_l"] = FAULTS / out["T_deepest_circuit"]
    out["fits_LQ"] = out["register"]["total"] <= RFI_LQ_2033
    out["fits_T"] = "yes" if out["T_deepest_circuit"] <= 1.5 * RFI_T_2033 else ("near" if out["T_deepest_circuit"] <= 3 * RFI_T_2033 else "no")
    out["fits_window"] = "yes" if out["wall_campaign_yr"] <= WINDOW_YR else ("near" if out["wall_campaign_yr"] <= 3 * WINDOW_YR else "no")
    return out


# ----------------------------------------------------------------------------------------------
# Published lemmas converted to the report convention (SU(4) Hamiltonian the paper costs; evolution only)
# ----------------------------------------------------------------------------------------------
def su4_converted(A, tau, eps, m=3, d=3):
    """Lemma 1, Lemma 3, Thm 2 and Thm 5 of arXiv:2507.22814 with the paper's couplings and q^2 dispersion.
    Conversion: Toffoli-derived T counts x 7/4; T_ROT -> report c_T. Uncontrolled, as published."""
    # --- PF2 ---
    r = math.ceil(math.sqrt(4 * tau ** 3 * sp.alpha2(A, m, d) / eps))
    gam = sp.Gamma(A)
    tofV = (4 * A * (A - 1) * (A - 2) / 3 + gam * 2 * sp.T_MCX(d * m)) / 4
    bd = sp.b_diag_pf2(A, m, tau, eps, d)
    tofT = sp.T_kin(A, m, bd, eps / (4 * r) / (3 * d * A), d) / 4
    nrot = r * gam
    ct = c_T(nrot)
    TV, TT = 7 * tofV + gam * ct, 7 * tofT
    pf = dict(steps=r, T_V_step=TV, T_T_step=TT, T_step=TV + TT, T_total=r * TV + (r + 1) * TT,
              T_paper_convention=sp.pf2(A, m, tau, eps, d)["T"],
              qubits=(d * m + 2) * A + sp.pf2(A, m, tau, eps, d)["ancilla"])
    # --- GQSP ---
    g = sp.gqsp(A, m, tau, eps, d)
    Q = g["Q"]
    e = eps / (4 * sp.lam_H(A, m) * tau)
    tof = (Q * sp.T_sel(A, m, e, d) + 2 * (Q + 1) * (sp.T_prep(A, m, e) - 4 * sp.T_ROT(e / 32)) + 2 * Q * sp.T_MCX(A + 2 * m + 10)) / 4
    nrot = 2 * (Q + 1) * 4 + 3 * (Q + 1)
    Tg = 7 * tof + nrot * c_T(nrot)
    gq = dict(Q=Q, lambda_MeV=sp.lam_H(A, m), T_step=Tg / Q, T_total=Tg, T_paper_convention=g["T"],
              qubits=(d * m + 2) * A + g["ancilla"])
    return dict(A=A, tau=tau, eps=eps, PF2=pf, GQSP=gq)


# ----------------------------------------------------------------------------------------------
def fmt(x):
    return f"{x:.2e}" if (abs(x) >= 1e4 or abs(x) < 1e-2) else (f"{x:.0f}" if abs(x) >= 100 else f"{x:.3g}")


def days(x):
    return f"{x:.1f} d" if x < 730 else f"{x / 365.25:.1f} yr"


def main():
    res = dict(conventions=dict(T_per_toffoli=7, eps_syn=EPS_SYN, t_gate_s=T_GATE_S, shot_overhead_s=T0_S,
                                faults_per_shot=FAULTS, target_LQ=RFI_LQ_2033, target_T=RFI_T_2033, window_yr=WINDOW_YR),
               current_triton=CURRENT, inputs=[dict(zip(("name", "value", "tag", "source", "note"), i)) for i in INPUTS])

    # 0. Table 5 reproduction
    res["table5_check"] = sp.check_table5()

    # 2. published lemmas, converted
    res["su4_converted"] = [su4_converted(A, t, e) for A in (12, 16, 40) for (t, e) in ((0.20, EPS_EVOL),)]
    res["su4_converted_at_paper_point"] = [su4_converted(A, sp.t_response(A, 3), 0.1) for A in (16, 40)]
    print("\n== Published SU(4) lemmas converted to the report convention (evolution only, tau=0.20, eps=0.00731) ==")
    print("| A | PF2 steps (bound) | PF2 T/step | PF2 T | PF2 LQ | GQSP queries | GQSP T/query | GQSP T | GQSP LQ |")
    print("|---|---|---|---|---|---|---|---|---|")
    for r in res["su4_converted"]:
        p, g = r["PF2"], r["GQSP"]
        print(f"| {r['A']} | {p['steps']} | {fmt(p['T_step'])} | {fmt(p['T_total'])} | {p['qubits']} | {g['Q']} | {fmt(g['T_step'])} | {fmt(g['T_total'])} | {g['qubits']} |")
    for r in res["su4_converted_at_paper_point"]:
        print(f"  paper point A={r['A']} t={r['tau']:.4f} eps=0.1: PF2 {fmt(r['PF2']['T_total'])} (paper conv. {fmt(r['PF2']['T_paper_convention'])}),"
              f" GQSP {fmt(r['GQSP']['T_total'])} (paper conv. {fmt(r['GQSP']['T_paper_convention'])})")

    # 1,3,4. chapter Hamiltonian, complete circuits
    cases = []
    for nuc in ("12C", "16O", "40Ar"):
        for route in ("PF2", "GQSP"):
            for steps in (("calibrated", "bound") if route == "PF2" else ("calibrated",)):
                cases.append(instance(nuc, route, steps=steps))
    eq5 = [instance(nuc, "GQSP", encoding="eq5") for nuc in ("12C", "16O", "40Ar")]
    res["fully_antisymmetrized_eq5_variant"] = eq5
    tritons = [instance("3H", "PF2", L=8), instance("3H", "PF2", L=6), instance("3H", "PF2", L=8, steps="bound"),
               instance("3H", "GQSP", L=8)]
    c12_6 = [instance("12C", "PF2", L=6)]
    res["instances"] = cases
    res["triton_first_quantized"] = tritons
    res["12C_6cubed"] = c12_6

    def row(c):
        s = c["route"] + ("" if c["route"] != "PF2" else (" cal." if c["steps_mode"] == "calibrated" else " bound"))
        return (f"| {c['nucleus']} {c['lattice']} | {s} | {c['encoding']} | {c['register']['total']} | {fmt(c['T_per_step'])} | "
                f"{fmt(c['T_evolution']['0.2'])} | {fmt(c['T_filter_per_attempt'])} | {fmt(c['T_deepest_circuit'])} | "
                f"{fmt(c['T_expected_per_accepted_shot']['0.2'])} | {fmt(c['shots_campaign'])} | {days(c['wall_first_d'])} | "
                f"{days(c['wall_campaign_d'])} | {c['eps_l']:.1e} | {c['fits_T']} | {c['fits_window']} |")
    hdr = ("| instance | method | encoding | LQ | T/step | T evol (0.20) | T filter/attempt | T deepest circuit | E[T]/accepted shot | "
           "campaign shots | first result | campaign | eps_l | fits 1e9 T | fits 5 yr |\n|" + "---|" * 15)
    print("\n== Complete circuits, chapter Hamiltonian (C_t != C_s, D), baseline assumptions ==")
    print(hdr)
    for c in tritons + c12_6 + cases + eq5:
        print(row(c))
    print(f"| 3H 6^3 occupation (current) | PF2 | charge | 908 | 1.7e6 | 1.8e8 | none | 2.07e8 | 2.07e8 | 1.64e+05 | 9.4 d | 250.0 d | 4.8e-10 | yes | yes |")

    # registers
    print("\n== Registers ==")
    for c in tritons[:1] + tritons[3:] + [x for x in cases if x["steps_mode"] == "calibrated"]:
        print(c["nucleus"], c["route"], c["register"])

    # 5/6. sensitivity to overlap and gap
    G2S = (0.9, 0.5, 0.2, 0.1, 0.01, 0.001)
    sens = []
    for nuc in ("12C", "16O", "40Ar"):
        for route in ("PF2", "GQSP"):
            for g2 in G2S:
                for gap in (1.0, 2.0, NUCLEI[nuc]["gap"], 8.0):
                    c = instance(nuc, route, gap=gap, gamma2=g2)
                    sens.append(dict(nucleus=nuc, route=route, gamma2=g2, gap=gap, T_deepest=c["T_deepest_circuit"],
                                     campaign_yr=c["wall_campaign_yr"], first_d=c["wall_first_d"]))
    res["sensitivity_overlap_gap"] = sens
    print("\n== Sensitivity: campaign wall in years (deepest-circuit T); rows overlap^2, columns gap (MeV) ==")
    for nuc in ("12C", "16O", "40Ar"):
        for route in ("PF2", "GQSP"):
            gaps = (1.0, 2.0, NUCLEI[nuc]["gap"], 8.0)
            print(f"\n{nuc} {route}" + (" (calibrated steps)" if route == "PF2" else "") + "\n| overlap^2 | "
                  + " | ".join(f"gap {g}" for g in gaps) + " |\n|" + "---|" * 5)
            for g2 in G2S:
                cells = []
                for gap in gaps:
                    x = next(x for x in sens if x["nucleus"] == nuc and x["route"] == route and x["gamma2"] == g2 and x["gap"] == gap)
                    cells.append(f"{x['campaign_yr']:.2g} ({fmt(x['T_deepest'])})")
                print(f"| {g2:g} | " + " | ".join(cells) + " |")

    # overlap needed for a five-year campaign at the baseline gap
    thr = []
    for nuc in ("12C", "16O", "40Ar"):
        for route, steps in (("PF2", "calibrated"), ("PF2", "bound"), ("GQSP", "calibrated")):
            lo, hi = 1e-6, 0.999
            f = lambda g: instance(nuc, route, steps=steps, gamma2=g)["wall_campaign_yr"] - WINDOW_YR  # noqa: E731
            if f(hi) > 0:
                g_need = None
            else:
                for _ in range(60):
                    mid = math.sqrt(lo * hi)
                    lo, hi = (mid, hi) if f(mid) > 0 else (lo, mid)
                g_need = hi
            thr.append(dict(nucleus=nuc, route=route, steps=steps, overlap2_needed_for_5yr=g_need))
    res["overlap_needed_for_window"] = thr
    print("\n== Overlap^2 needed for a five-year campaign at the baseline gap ==")
    for x in thr:
        print(x)

    # step-count sensitivity for PF
    stepsens = []
    for nuc in ("3H", "12C", "16O", "40Ar"):
        A = NUCLEI[nuc]["A"]
        stepsens.append(dict(nucleus=nuc, calibrated=steps_calibrated(A, 0.2), bound=steps_bound(A, 0.2),
                             ratio=steps_bound(A, 0.2) / steps_calibrated(A, 0.2)))
    res["pf_steps_tau0.2"] = stepsens
    print("\n== PF2 steps at tau=0.20, eps=0.00731 ==")
    for s in stepsens:
        print(s)

    # what the first-quantized triton could do with the 1e9 budget
    t3 = tritons[0]
    spare = (RFI_T_2033 - t3["T_prep_deterministic"]) / t3["T_per_step"]
    res["triton_8cubed_steps_within_1e9"] = spare
    res["triton_8cubed_tau_within_1e9_MeV_inv"] = spare * 0.002
    print(f"\nFirst-quantized triton 8^3: {spare:.3g} steps of 0.002 fit in 1e9 T (tau = {spare * 0.002:.0f} MeV^-1, "
          f"energy resolution 2 pi/tau = {2 * math.pi / (spare * 0.002) * 1e3:.1f} keV)")

    with open(os.path.join(HERE, "results.json"), "w") as f:
        json.dump(res, f, indent=1, default=float)


if __name__ == "__main__":
    main()
