"""Ch. 3 (Mu2e / 0nubb): every named number in app11_mu2e_0nubb.tex, as the chapter
states it, so that editing one without the others fails loudly.

The generic contract tests (test_common.py) already check provenance, the breakdown
sum, and closeness to PUBLISHED. These pin the prose intermediates. Tolerances are
stated per assertion: 1e-9 where the chapter quotes an exact product, 5% where it
prints two significant figures, 10% where one. Numbers the chapter's own inputs do
not reproduce are strict xfails pointing at NEEDS_AUTHOR.md; if an edit makes
them pass, the XPASS fails the suite and the item must be re-examined.

2026-09-28 pass (rulings R3, R4 + the mechanical batch; apply_log/ch03.md): epsilon_l is
0.1/N per shot; every band carries the exact N_Trotter = 1000pi-4000pi and rounds once at
print. The rounding-chain xfails are gone; `test_printed_strings_are_in_the_tex` pins the
model's printed strings to the chapter verbatim.

2026-09-30 pass (R11 rulings; apply_log/r11_ch03.md): QROM passes and Delta_small are stated,
the 2028 headline prints 5.4-5.7e4 with the oracle load separate, the 2033 box counts 27Al shots
only, the flagship register is 'at the envelope', the utility fraction is 0.8%, the full-Gamma
multiple is x8.4-14. The strict xfails those rulings resolve are now passing asserts.

2026-09-30 follow-up (R11b; apply_log/r11b_ch03.md): the QROM load is priced at R5's 7 T per
Toffoli, 7(D-1) per pass; the priced condition is D <~ 2e3 (0.82-0.99x) and the cap crossing
D = 2.0-3.3e3 is derived in the model (d_break is no longer a Stated input).

2026-09-30 R12 item 13 (apply_log/r12_ch03.md): the Gamma<=128 insertion is counted gate by gate
(rotation-tree PREPARE, controlled unary-iteration SELECT, d = 2-4 QSP queries, R-TOL, 7 T/Toffoli):
1.1-2.3e4 T (was the unsourced 5.4-5.7e4). Gamma = 256 / 512 / full, the QROM composed fractions,
the D crossing (3.7-6.3e3), the wall times and the 55-LQ insertion peak follow.

2026-09-30 R13 'one phase' ruling (apply_log/r13_ch03.md): the benchmark's <Psi_0|O|Psi_0> is the
degree-1 polynomial, one controlled BE query and no QSP phases (GSLW 1806.01838, Martyn 2105.02859).
Headline 5.4e3 T (0.054x); Gamma = 256 / 512 at 0.11x / 0.23x; full real one-body 0.93x, complex x1.9;
composed 0.33-0.47x, D crossing 4.5-6.8e3, load-stage peak 65 LQ, wall 4-16 s. The R12 reading
(d = 2-4 queries) is kept as the qsp_reading_* sensitivity (0.11-0.23x, doubled 0.23-0.47x).

2026-09-30 R14 (apply_log/r14_ch03.md): (1a) the 4He dipole shot factor (lambda/|<O>|)^2 = 77 is derived
(Gamma = 25, lambda = 15.4, <O> = 1.75; 70-130 for any state in the space), shots 4.6e4; (1d) depth-capped
MLAE, m_max = 13-5 at D ~ 1e2-1e3, 2.8e2-1.9e3 circuits; (2) option A: 0nubb by daughter projection in one
register (N_proj = 2 unchanged in T, registers halved); (3) QROM address sized to D <~ 2e3, peak 61 LQ.

2026-10-01 R15 (apply_log/r15_ch03.md): (1) the identity term is subtracted classically: Gamma = 24,
lambda' = 7.79, shot factor (lambda'/|<O>|)^2 = 20 (18-33 any state), 2.0e4 Hadamard-test shots at the band top (R16); 27Al sd 1.45,
3.5e3 shots per scan. (2) The complete sparse load of Fomichev et al. 2310.18410, one per shot: D <~ 600
(fits to 603, crosses at 604), load-stage peak 87 LQ, MLAE m_max = 5 at D = 1e2 (S_0 on 58 qubits, peak 114). (3) Mo/Ru
is jj44 protons + jj55 neutrons = 54 modes: x6-25 (x3-6 large gap, x25-51 small). (4) The 2033 insertions
are one controlled BE query with Gamma counted by selection rules: jj44 0nubb Gamma = 40016 (2501 terms),
4.0-4.1e6 T, 0.04-0.2% of the pair prep; 27Al sd Mu2e two-body 12900 strings, 9.3e5-3.9e6 T.

2026-10-02 r23 (author ruling "we don't want to anywhere assume we have multiple machines";
apply_log/r23_ch03.md): every wall time is shots x per-shot time, serial on one machine. The chapter no
longer says "shot-parallel", "serialized" or "divides across ... machines"; the 2033 27Al campaign
(8.5 days-1.8 months) is printed as at most 2.9% of the 5-year campaign horizon. No T, LQ or shot count moves.
"""

import math
import pathlib

import pytest

from estimates import common as c
from estimates import ch03_mu2e_0nubb as ch3

NA = "see NEEDS_AUTHOR.md"
TEX = pathlib.Path(__file__).resolve().parents[3] / "applications" / "app11_mu2e_0nubb.tex"


@pytest.fixture(scope="module")
def a():
    return ch3.Assumptions()


@pytest.fixture(scope="module")
def tex():
    return TEX.read_text() if TEX.exists() else None


@pytest.fixture(scope="module")
def r28(a):
    return ch3.model(a, "2028")


@pytest.fixture(scope="module")
def r33(a):
    return ch3.model(a, "2033")


@pytest.fixture(scope="module")
def rcd(a):
    return ch3.model(a, "codesign")


def _is(x, y, rel=1e-9):
    assert math.isclose(x, y, rel_tol=rel), (x, y)


def _div(r, k):
    return (r[0] / k, r[1] / k)


def _rng(model, quoted, rel):
    assert len(model) == 2 and len(quoted) == 2
    for m, q in zip(model, quoted):
        assert math.isclose(m, q, rel_tol=rel), (model, quoted, rel)


# ---------------------------------------------------------------- provenance

def test_every_input_is_tagged(a):
    import dataclasses
    for f in dataclasses.fields(a):
        assert isinstance(getattr(a, f.name), c.Tagged), f.name


def test_stated_and_uncited_inputs_are_tagged_honestly(a):
    # the jj44 insertion and the gap bands are asserted, not derived; the Gamma<=128 insertion is
    # derived since R12 (no stated fraction input remains)
    assert not hasattr(a, "gamma_trunc_fraction")
    # R13: one BE query is derived here (Hadamard test on the controlled BE), declared as an assumption
    assert a.be_queries_2028.prov is c.Provenance.ASSUMED and a.be_queries_2028.value == 1
    assert a.select_and_deficit.prov is c.Provenance.CITED and a.select_and_deficit.value == 1
    assert a.insertion_t_per_toffoli.value == c.T_PER_TOFFOLI["textbook"] == 7
    # R15 ruling 4: the 2033 insertions are derived; the stated 5.4-5.6e7 / 4.4-4.6e8 survive only as a record
    assert not hasattr(a, "insertion_jj44_sym") and not hasattr(a, "insertion_jj44_dense")
    assert ch3.INSERTION_JJ44_PRE_R15["symmetry_restricted"] == (5.4e7, 5.6e7)
    assert a.gap_mev.prov is c.Provenance.STATED
    assert a.gap_large_mev.prov is c.Provenance.STATED
    # R11 (2026-09-30): the small gap is printed (app11:162); R15: the lookup pass count is a record only,
    # the load is one complete load per shot (derived, declared as an assumption)
    assert a.gap_small_mev.prov is c.Provenance.STATED and a.gap_small_mev.value == 0.5
    assert a.qrom_passes.prov is c.Provenance.ASSUMED and a.qrom_passes.value == (2, 3)
    assert a.loads_per_shot.prov is c.Provenance.ASSUMED and a.loads_per_shot.value == 1
    assert a.d_amplitudes.prov is c.Provenance.STATED and a.d_amplitudes.value == 600
    # R11b: D-1 Toffolis from the source, priced at 7 T/Toffoli (R5); the crossing D is derived
    assert a.qrom_t_per_amplitude.prov is c.Provenance.CITED and a.qrom_t_per_amplitude.value == 7
    assert a.qrom_t_per_amplitude.value == c.T_PER_TOFFOLI["textbook"]
    assert not hasattr(a, "d_break")
    assert a.utility_fraction.prov is c.Provenance.STATED and a.utility_fraction.value == 0.6
    # the synthesis budget (R-TOL, report-wide) and the fault budget are stated in the prose (app11:90, :93)
    assert a.eps_synth_total.prov is c.Provenance.STATED and a.eps_synth_total.value == 0.01
    assert a.fault_budget.prov is c.Provenance.STATED and a.fault_budget.value == 0.1


def test_fault_budget_is_applied_everywhere(a, r28, r33, rcd):
    # "the quoted eps_l is 0.1/N_gate^shot" (app11:93, ruling R3)
    for r in (r28, r33, rcd):
        _is(r.epsilon_l[0], 0.1 / r.hard_ops[1])
        _is(r.epsilon_l[1], 0.1 / r.hard_ops[0])
    with pytest.raises(ValueError, match="fault_budget"):
        ch3.Assumptions(fault_budget=c.Stated(2.0, "x"))


def test_c_T_goes_through_common_synthesis(a, r33):
    # R-TOL (2026-09-29): no fixed c_T; each circuit's c_T is the RUS fit at eps_rot_for(N_rot)  app11:90
    assert not hasattr(a, "c_T")
    assert a.eps_synth_total.value == c.EPS_SYN == 1e-2
    i = r33.intermediates
    _is(i["c_T_rus_at_eps_1e4"], c.t_per_rotation(1e-4))
    # "giving eps_rot = 2.5e-6-5.0e-5 and c_T = 26-31 T per rotation across the survey's cells" (R15: the
    # top end is the jj55 pair, the largest cell once Mo/Ru is 54 modes)
    lo, hi = i["c_T_survey_range"]
    assert 25.6 <= lo <= 25.7 and 30.6 <= hi <= 30.7, (lo, hi)
    _rng(i["eps_rot_survey"], (4.956e-5, 2.464e-6), 1e-3)
    # "Were the errors instead to add coherently, c_T would rise to 42-52 ... by ~70%"
    clo, chi = i["c_T_coherent_range"]
    assert 42 <= clo <= 42.2 and 52.0 <= chi <= 52.1, (clo, chi)
    for x in i["c_T_coherent_over_randomized"]:
        assert 1.6 <= x <= 1.8                                      # "~70%"


def test_rtol_each_band_end_is_its_own_circuit(a, r33, rcd):
    # every band end: N_rot = N_reg N_Trotter N_orb^4, eps_rot = sqrt(1e-2/N_rot), c_T = 1.15 log2(1/eps) + 9.2
    nt = r33.intermediates["n_trotter"]
    for key, n_reg, n_orb in (("prep_al_sd", 1, 6), ("prep_pf", 2, 10), ("prep_jj44", 2, 11),
                              ("prep_jj55", 2, 16), ("prep_mo_ru", 2, 13.5)):
        for end in (0, 1):
            n = n_reg * nt[end] * n_orb ** 4
            eps = math.sqrt(1e-2 / n)
            _is(r33.intermediates[key][end], n * (1.15 * math.log2(1 / eps) + 9.2))
    # inventory values (apply_log/rtol_shared.md)
    _rng(r33.intermediates["eps_rot_al_sd"], (4.956e-5, 2.478e-5), 1e-3)
    _rng(r33.intermediates["c_T_al_sd"], (25.646, 26.796), 1e-4)
    _rng(r33.intermediates["c_T_pf"], (27.916, 29.066), 1e-4)
    _rng(r33.intermediates["c_T_jj44"], (28.232, 29.382), 1e-4)
    _rng(rcd.intermediates["c_T_jj55"], (29.475, 30.625), 1e-4)
    _rng(rcd.intermediates["c_T_mo_ru"], (28.911, 30.061), 1e-4)              # R15: 54 modes (was 31.2-32.4)
    _rng(rcd.intermediates["c_T_flagship"], (36.643, 39.618), 1e-4)
    _rng(rcd.intermediates["eps_rot_flagship"], (6.553e-8, 1.091e-8), 1e-3)
    # the 2028 insertion (R13): one query, 254 rotations and 127 Toffolis, its own eps_rot;
    # the R12 sensitivity ends (511/272, 1021/538) are their own circuits too
    r28 = ch3.model(a, "2028")
    for end in (0, 1):
        eps = math.sqrt(1e-2 / 254)
        _is(r28.intermediates["t_insertion"][end], 254 * (1.15 * math.log2(1 / eps) + 9.2) + 7 * 127)
    for end, n_rot, n_tof in ((0, 511, 272), (1, 1021, 538)):
        eps = math.sqrt(1e-2 / n_rot)
        _is(r28.intermediates["qsp_reading_t"][end], n_rot * (1.15 * math.log2(1 / eps) + 9.2) + 7 * n_tof)
    # a gap variant is its own circuit: Delta = 4 MeV carries its own eps
    n = 2 * math.pi * 250 / 0.5 * 16 ** 4
    _is(ch3.prep_t(a, 2, 16, (4.0, 4.0))[0], n * c.t_per_rotation(c.eps_rot_for(n)))


def test_c_proj_is_the_declared_working_assumption(a):
    assert a.c_proj.prov is c.Provenance.ASSUMED
    _is(a.c_proj.lo, math.pi)
    _is(a.c_proj.hi, 2 * math.pi)


def test_no_breakdown_primitive_claims_compiled(r28, r33, rcd):
    # the 2028 insertion is an explicit gate-level count since R12: COMPILED, 'derived here' except
    # SELECT, whose L-1 Toffolis are Babbush's circuit; the state-prep lines stay SCALING. R15: the 2033
    # sd insertion is a gate-level count too (one controlled query at the derived Gamma)
    for p in r28.breakdown:
        assert p.status is c.CircuitStatus.COMPILED, p
        assert p.src in ("derived here", "Babbush_PRX_2018", "arxiv_2310_18410", "arxiv_1812_00954"), p   # G2 + M3
    for r in (r33, rcd):
        for p in r.breakdown:
            if p.name.startswith("operator_insertion"):
                assert p.status is c.CircuitStatus.COMPILED and p.src in ("derived here", "Babbush_PRX_2018"), p
            else:
                assert p.status is c.CircuitStatus.SCALING, p


def test_bad_era_rejected(a):
    with pytest.raises(ValueError):
        ch3.model(a, "2040")


def test_gap_ordering_enforced():
    with pytest.raises(ValueError):
        ch3.Assumptions(gap_small_mev=c.Uncited(3.0))


def test_scalar_inputs_reject_ranges():
    # the formulas use eps_synth_total and qubits_per_orb as scalars; a range must fail loudly, not TypeError
    with pytest.raises(ValueError, match="eps_synth_total"):
        ch3.Assumptions(eps_synth_total=c.Stated((1e-3, 1e-2), "x"))
    with pytest.raises(ValueError, match="qubits_per_orb"):
        ch3.Assumptions(qubits_per_orb=c.Stated((2, 4), "x"))


# ------------------------------------------------------- the formula itself

def test_ho_ladder_orbital_counts(rcd):
    # "An e_max = 8 (10) harmonic-oscillator basis carries 165 (286) spatial orbitals"  app11:135
    i = rcd.intermediates
    assert i["ho_orbitals"] == {3: 20, 4: 35, 8: 165, 10: 286, 14: 680}
    assert i["ho_orbitals_converged"] == i["ho_orbitals_quoted"] == (165, 286)


def test_state_prep_horizon(r33):
    i = r33.intermediates
    _is(i["tau_gs_gev_inv"][0], 500)              # "tau_gs ~ 1/Delta ~ 500-1000 GeV^-1"   app11:87
    _is(i["tau_gs_gev_inv"][1], 1000)
    _rng(i["tau_gs_over_ch2_horizon"], (5, 10), 1e-9)   # "5-10x the response horizon of Ch. 2"
    # "N_Trotter = c_proj tau_gs / dt = 1000pi-4000pi ~ 3.1e3-1.3e4"                   app11:88
    _is(i["n_trotter"][0], 1000 * math.pi)
    _is(i["n_trotter"][1], 4000 * math.pi)
    _rng(i["n_trotter"], (3.1e3, 1.3e4), 0.05)
    assert i["n_trotter_quoted"] == (3.1e3, 1.3e4)


def test_worked_line_27al(r33):
    # "27Al has N_reg=1, N_orb=6, N_Trotter = 1000pi-4000pi, so eps_rot = 2.5-5.0e-5 and
    #  c_T = 25.6-26.8, and Eq. gives ~1.0-4.4e8 T"   app11:90 (R-TOL)
    i = r33.intermediates
    assert i["worked_line"] == i["prep_al_sd"]                    # exact carry (R4)
    n = (1000 * math.pi * 6 ** 4, 4000 * math.pi * 6 ** 4)
    _is(i["worked_line"][0], n[0] * c.t_per_rotation(c.eps_rot_for(n[0])))
    _is(i["worked_line"][1], n[1] * c.t_per_rotation(c.eps_rot_for(n[1])))
    _rng(i["worked_line"], (1.044e8, 4.364e8), 1e-3)                                    # app11:139
    # the record of what rounding N_Trotter first would give (each end at its own eps)
    _rng(i["worked_line_round_first"], (1.03e8, 4.52e8), 0.01)


# ------------------------------------------------------------- 2028 register

def test_2028_register(r28):
    i = r28.intermediates
    assert i["system_qubits"] == 40               # "4 N_orb = 40 system qubits"          app11:130
    assert i["n_anc"] == (80, 130)                # "~80-130 ancilla"                      app11:130
    assert i["lq_sum"] == (120, 170)              # "~120-170"                             app11:185
    assert r28.lq == (120, 170)
    # R15 ruling 2 + the R14 sizing rule: the complete load is sized to the printed D <~ 600: L = 10
    # enumeration qubits, 5L - 3 = 47 ancillae (Fomichev initial.tex:355): 10 + 19 + 18
    assert i["load_address_qubits"] == 10 and 2 ** 9 < 600 <= 2 ** 10
    assert i["load_ancillae"] == 47 == 10 + (2 * 10 - 1) + (2 * 10 - 2)
    assert i["peak_lq_load_no_angles"] == 40 + 47 == 87
    # M3 ruling (2026-10-04): the angle step holds 10 enumeration + 3b - 1 = 38 qubits (37 of them idle load ancillae)
    assert i["peak_lq_load"] == 40 + 10 + 38 == 88
    # insertion stage 40 system + 7 LCU address + 7 unary-iteration workspace + 1 test = 55; no phase ancilla
    one = i["insertion_counts"]
    assert one["address_qubits"] == 7 and one["select_and_ancillae"] == 7
    assert i["peak_lq_insertion"] == 40 + 7 + 7 + 1 == 55
    assert i["peak_lq"] == 88                     # "a peak of 88 LQ at the load stage"
    # MLAE: S_0 on 40 system + 10 enumeration + 7 LCU address + 1 test = 58 (the identification register is
    # CNOT-computed and CNOT-uncomputed from the unchanged system, clean for any input), 56 AND ancillae
    assert i["mlae_peak_lq_reflection"] == 58 + 56 == 114
    assert i["mlae_peak_lq"] == 114 + 13 == 127   # + the reused phase-gradient state (M3 ruling)
    assert i["rfi_lq_2028"] == (150, 250)
    assert r28.lq[1] <= i["rfi_lq_2028"][1]     # "sits far inside a [150,250]-LQ machine"  app11:130
    assert i["peak_lq"] < i["rfi_lq_2028"][0]   # "the envelope sets the machine size, not a floor"

# ------------------------------------------------------------ 2028 hard ops

def test_2028_topdown_estimate(r28):
    # R13: the Hadamard test on the controlled BE is one query, so T_shot ~ T_BE ~ 1e3-1e4    app11:122
    i = r28.intermediates
    assert i["be_queries"] == 1 and i["qsp_phase_steps"] == 0
    _rng(i["t_shot_topdown"], (1e3, 1e4), 1e-9)


def test_2028_bottom_up_costing_is_the_box(r28):
    # counted gate by gate, "derived here"; one controlled query (R13)               app11:124
    i = r28.intermediates
    one = i["insertion_counts"]
    assert (one["k"], one["rot_per_prepare"]) == (7, 127)            # k = ceil(log2 128), 2^k - 1
    assert one["n_rot"] == 254 and one["n_rot_phase"] == 0           # "2 x 127 = 254 rotations"
    assert one["n_toffoli"] == 127 and one["toffoli_phase"] == 0     # "127 Toffolis"
    assert round(one["c_T"], 1) == 17.6                              # "(17.6 T each)"
    _is(one["eps_rot"], math.sqrt(1e-2 / 254))
    _rng(i["t_insertion"], (5362.88, 5362.88), 1e-6)
    assert ch3.round_sig(i["t_insertion"][0] / 1e3) == 5.4                   # "5.4e3"   app11:122,124,181,236
    assert ch3.round_sig(i["gamma_trunc_fraction"][0]) == 0.054              # "0.054x the cap"  app11:124,182
    # referee G2 + M3 ruling (2026-10-04): the headline is the executed shot, insertion + one complete load with
    # its angle rotations at D = 1e2 .. 600; 1.005x the cap at D = 600, last fit D = 596
    _rng(r28.hard_ops, (i["t_insertion"][0] + i["load_t_range"][0] + i["load_angle_t"][0],
                        i["t_insertion"][1] + i["load_t_range"][1] + i["load_angle_t"][1]), 1e-12)
    _is(r28.breakdown_total(), r28.hard_ops[0])
    assert r28.hard_ops[1] < 1.01e5
    # the top-down band 1e3-1e4 (app11:122) contains the derived insertion
    assert i["t_shot_topdown"][0] <= i["t_insertion"][0] <= i["t_insertion"][1] <= i["t_shot_topdown"][1]


def test_2028_qsp_reading_sensitivity(r28):
    # the R12 reading, kept as a sensitivity: d = 2-4 queries plus d + 1 phase steps (printed at app11:124
    # as the counterfactual '1.1-2.3e4 T (0.11-0.23x)'); doubled degree 0.23-0.47x
    i = r28.intermediates
    assert i["qsp_reading_degree"] == (2, 4)
    lo, hi = i["qsp_reading_counts"]
    assert (lo["n_rot"], hi["n_rot"]) == (511, 1021) and (lo["n_toffoli"], hi["n_toffoli"]) == (272, 538)
    _rng(i["qsp_reading_t"], (11200.93, 22927.91), 1e-6)
    assert ch3.print_range(_div(i["qsp_reading_t"], 1e4)) == "1.1--2.3"
    assert ch3.print_range(i["qsp_reading_fraction"]) == "0.11--0.23"
    assert ch3.print_range(i["qsp_reading_doubled_fraction"]) == "0.23--0.47"


def test_2028_pre_r12_print_is_reconstructed(r28):
    # the unsourced 0.54-0.57x is 8 queries x (254 rotations at a fixed eps_rot = 1e-4 + 127 Toffolis
    # at 4 T / 7 T); under the current rules the same 8 queries give ~0.47x (triage reconstruction)
    i = r28.intermediates
    assert tuple(round(x, 2) for x in i["pre_r12_printed_reconstruction"]) == (0.54, 0.57)
    assert round(i["gamma128_fraction_at_8_queries"], 2) == 0.47


def test_2028_gamma_scaling(r28):
    # the same one-query count at larger Gamma (not linear: tree padding, c_T grows with N_rot)
    i = r28.intermediates
    assert ch3.round_sig(i["gamma_256_fraction"][0]) == 0.11     # "0.11x at Gamma = 256"  app11:124
    assert ch3.round_sig(i["gamma_512_fraction"][0]) == 0.23     # "0.23x at Gamma = 512"  app11:124,183
    # full one-body Gamma from the 40 JW modes: real 2 C(40,2) + 40, complex 4 C(40,2) + 40
    assert i["gamma_full"] == (1600, 3160)
    lo, hi = i["gamma_full_counts"]
    assert (lo["k"], lo["n_rot"], lo["n_toffoli"]) == (11, 4094, 1599)     # Gamma = 1600, one query
    assert (hi["k"], hi["n_rot"], hi["n_toffoli"]) == (12, 8190, 3159)     # Gamma = 3160, one query


def test_2028_gamma_full_upper_multiple(r28):
    # "fits at 0.93x" (real) and "exceeds the cap by x1.9" (complex)                     app11:124,183
    lo, hi = r28.intermediates["gamma_full_multiple"]
    assert (round(lo, 2), round(hi, 1)) == (0.93, 1.9), (lo, hi)
    assert lo < 1.0 < hi


def test_2028_oracle_load(r28):
    # R15 ruling 2: the complete sparse load (Fomichev 2310.18410), (2L-2)D + 2^(L+1) + D Toffolis at 7 T,
    # one load per shot
    i = r28.intermediates
    assert i["loads_per_shot"] == 1
    assert i["load_t_at_d"][100.0] == 7 * (12 * 100 + 256 + 100) == 10892      # "1.1e4 T at D = 1e2"
    assert i["load_t_at_d"][600] == 7 * (18 * 600 + 2048 + 600) == 94136
    assert ch3.round_sig(i["load_t_at_d"][100.0]) == 1.1e4
    assert round(i["sos_over_lookup_pass"][100.0]) == 16                      # "the price of 16 bare lookups"
    # composed with the 5362.9-T insertion (record without the angle term): 0.16 at D = 1e2 to 0.99 at D = 600;
    # fits to 603, crosses at 604. With the angle term (M3 ruling, the headline): fits to 596, crosses at 597
    t0, t1 = i["t_insertion"]
    _rng(i["hard_ops_with_oracle_load"], (t0 + 10892, t1 + 94136), 1e-12)
    assert i["hard_ops_no_angles"] == i["hard_ops_with_oracle_load"]
    assert r28.hard_ops == i["hard_ops_with_load_angles"]                     # M3 ruling: the headline
    assert ch3.print_range(i["composed_fraction_range"]) == "0.16--0.99"
    assert i["d_break_no_angles"] == 603 and i["d_cross_no_angles"] == 604 and i["d_priced"] == 600
    assert t0 + 7 * ch3.sos_load_toffoli(603) <= 1e5 < t0 + 7 * ch3.sos_load_toffoli(604)
    assert i["d_break"] == 596 and i["d_cross"] == 597
    # record (retired print since referee M3): at D = 1e3 the load alone is 1.5x the cap
    assert round(i["load_t_at_d"][1000.0] / 1e5, 1) == 1.5 and i["composed_fraction_at_d_hi"] > 1
    # the R11b-R14 lookup-only record
    lk = i["lookup_r14"]
    assert lk["passes"] == (2, 3) and lk["load_plausible"] == (1386, 20979)
    assert ch3.print_range(lk["composed_fraction_at_2e3"]) == "0.33--0.47"
    assert (ch3.round_sig(lk["d_break"][0] / 1e3), ch3.round_sig(lk["d_break"][1] / 1e3)) == (4.5, 6.8)



def test_2028_oracle_load_is_printed_separately(r28, tex):
    # R11 ruling 2028-headline-oracle-load (A): the headline is the insertion and the box prints the
    # oracle load on its own. R15: one complete load, 1.1e4 at D = 1e2 to 9.4e4 at D = 600
    i = r28.intermediates
    assert tuple(ch3.round_sig(x) for x in i["load_t_range"]) == (1.1e4, 9.4e4)
    # referee G2 (2026-10-04): the load IS in the headline now (was the insertion alone); M3 ruling: + angles
    assert tuple(ch3.round_sig(x) for x in r28.hard_ops) == (1.7e4, 1.0e5)
    assert tuple(ch3.round_sig(x) for x in i["hard_ops_no_angles"]) == (1.6e4, 9.9e4)
    if tex is None:
        pytest.skip("chapter .tex not on disk")
    # R16 sweep: box shows final numbers only
    # style pass 2026-10-08: box wording follows rule 10/12 (numbers unchanged)
    assert ("$1.7\\times 10^{4}$--$1.0\\times 10^{5}$: $5.4\\times 10^{3}$ (insertion, one "
            "query) $+$ $1.1\\times 10^{4}$--$9.4\\times 10^{4}$ (state preparation, $D=10^{2}$--$600$) $+$ "
            "$7.7\\times 10^{2}$--$1.0\\times 10^{3}$ (amplitude rotations)") in tex
    # cut pass 2026-10-06: the requirements rows that restated the boxes are gone; one row points at the boxes
    # style pass 2026-10-08: retired, the sentence was cut (rule 1/2); the value stays pinned on the model
    assert "peak $88$ LQ, $127$ under amplitude estimation" in tex
    assert "lower bound" not in tex.split("The 2028 benchmark}")[1].split("The 2033 tiers")[0]
    assert "$61$ LQ" not in tex and "$65$ LQ" not in tex
    assert "$63$ LQ" not in tex and "QROM address" not in tex
    assert "the register needs $55$ LQ during the insertion" in tex and "a peak of $116$" not in tex
    assert "2 QROM work" not in tex and "$2$ work at load" not in tex
    assert "$5.4$--$5.7" not in tex and "0.54" not in tex and "walk-PREPARE" not in tex
    # R13: 'query' and 'phase' defined; one query, no phases
    assert "A \\emph{query} is one application of the BE unitary" in tex
    # style pass 2026-10-08: retired, the sentence was cut (rule 1/2); the value stays pinned on the model
    assert "with one query and no QSP phases" in tex
    assert "wrapped in $\\mathcal{O}(\\log 1/\\epsilon)$ QSP phases" not in tex
    assert "T_{\\rm BE}\\cdot\\lceil\\log_2(1/\\epsilon)\\rceil" not in tex
    # style pass 2026-10-08: retired, the sentence was cut (rule 1/2); the value stays pinned on the model
    assert "T_{\\rm shot}^{2028}\\sim T_{\\rm BE}" not in tex
    # style pass 2026-10-08: retired, the sentence was cut (rule 1/2); the value stays pinned on the model
    assert "the insertion is $0.054\\times$ the 2028 target" in tex
    # E30 cap pass: the Gamma = 256/512/1600/3160 what-ifs and the d = 2-4 QSP counterfactual are no longer
    # printed; their values stay pinned on the model above (gamma_*_fraction, qsp_reading_*)
    assert "2000$--$3200" not in tex
    assert "$254d+d+1" not in tex
    assert "\\sim 5\\times 10^{4}$ T/Toffoli" not in tex and "unqualified $5\\times 10^{4}$" not in tex
    # R15: the complete load prices the box; the lookup-only wording is gone
    # style pass 2026-10-08: wording follows rules 7/12 (oracle load -> state preparation, cap -> 2028 target)
    assert "state preparation $D\\lesssim 600$ configurations" in tex
    assert "The executed shot is then $0.17\\times$ the 2028 target at $D=10^{2}$" in tex
    assert ("and $1.005\\times$ at $D=600$, above it by less than the precision of the count, "
            "and it fits up to $D=596$") in tex and "$D=603$" not in tex and "$D=604$" not in tex
    assert "lookup passes" not in tex and "lookup-only price" not in tex and "2--3" not in tex.split("Sizes &")[1][:400]
    assert "4D-4" not in tex and "D\\sim 10^{4}" not in tex
    # R16 sweep: the Gamma = 512 / 1600 what-ifs and the 0.054x fraction live in the prose only (pinned above)
    assert "insertion LCU $\\Gamma\\lesssim 128$ (dipole $\\Gamma{=}24$)" in tex



def test_2028_shots_and_wall_time(r28, tex):
    # R15 ruling 1: shots carry (lambda'/|<O>|)^2 with the identity subtracted classically; R16: at the upper
    # end of the any-state band (32.76), the 27Al convention, not the closed-shell 19.85
    i = r28.intermediates
    _is(i["dipole_shot_factor"], i["dipole_shot_factor_band"][1])
    assert round(i["dipole_shot_factor"], 2) == 32.76 and round(i["dipole_shot_factor_closed_shell"], 2) == 19.85
    assert ch3.round_sig(i["shots_closed_shell"]) == 12000                       # R15 print, record
    assert i["shots_per_operator_bare"] == 100 and i["shots_bare"] == 600        # the R13 eps^-2 x 6
    _is(i["shots_per_operator"], 100 * i["dipole_shot_factor"])
    _is(i["shots"], 600 * i["dipole_shot_factor"])
    assert ch3.round_sig(i["shots_per_operator"]) == 3300                        # "3.3e3 shots per operator"
    assert ch3.round_sig(i["shots"]) == 20000                                    # "2.0e4"   app11:123,box
    assert ch3.round_sig(i["shots_r14"]) == 46000                                # the R14 print, record
    assert r28.shots == (i["shots"], i["shots"])
    # R8: 1 us per T-gate. A shot is the insertion plus one complete load at D = 1e2 .. 600
    load = i["load_t_range"]
    # r25 (R9, author approved): 0.1 ms per shot (init + readout + decode), was the 1 ms bound of r17
    assert i["shot_overhead_s"] == 1e-4 and a_default().t_gate_s.value == 1e-6
    ins = i["t_insertion"]
    ang = i["load_angle_t"]                                                  # M3 ruling: in the shot
    _rng(i["t_shot_s"], ((ins[0] + load[0] + ang[0]) / 1e6 + 1e-4, (ins[1] + load[1] + ang[1]) / 1e6 + 1e-4), 1e-12)
    _rng(i["t_shot_s"], (r28.hard_ops[0] / 1e6 + 1e-4, r28.hard_ops[1] / 1e6 + 1e-4), 1e-12)   # referee G2
    _rng(i["t_shot_insertion_only_s"], _div(ins, 1e6), 1e-12)
    # T supply alone: 0.017-0.10 s per shot, 337-1978 s (wall_time_serial_supply_s; M3 ruling, was 321-1958)
    lo, hi = i["wall_time_serial_supply_s"]
    assert round(lo) == 337 and round(hi) == 1978, (lo, hi)
    # r25 (R9/R10): at D = 1e2 F* = 7.1 < 10, so the T-depth sets the shot (2382 x 10 us + 0.1 ms = 0.024 s);
    # at D = 600 F* = 14.3 and the T supply sets it (0.1006 s). Wall 470-1978 s = 7.8-33 min ("~8-33 min").
    lo, hi = r28.wall_time_s
    assert r28.wall_time_s == i["wall_campaign_s"] == i["wall_time_corrected_s"]
    assert round(lo) == 470 and round(hi) == 1978 and round(lo / 60) == 8 and round(hi / 60) == 33, (lo, hi)
    assert (round(i["f_star_per_end"][0], 1), round(i["f_star_per_end"][1], 1)) == (7.1, 14.3)
    nlo, nhi = i["wall_time_no_overhead_s"]
    assert round(nlo / 60) == 6 and round(nhi / 60) == 33
    if tex is not None:
        assert "$\\sim 0.024$--$0.10\\,$s per Hadamard-test shot" in tex
        assert "$\\sim 8$--$33$ min (Hadamard test) or $\\sim 68$ s (MLAE at $D=10^{2}$) on one machine" in tex
        assert "(at $1\\,\\mu$s per T-gate and $0.1$\\,ms per shot, convention in Ch.~\\ref{ch:overview}" in tex
        # the 0.1 ms is sourced in the chapter (r17 verifier); the 1 ms bound is retired (r25)
        assert "$0.1$\\,ms per shot for initialization, readout and decoding, as on current superconducting hardware~\\cite{Google_QEC_below_threshold}" in tex
        assert "which we bound at $1$\\,ms" not in tex and "1$\\,ms per-shot bound" not in tex
        assert ("$2.0\\times 10^{4}$ Hadamard-test shots, or $7.9\\times 10^{2}$ MLAE circuits at "
                "$D=10^{2}$ ($m_{\\max}=5$, $8.7\\times 10^{4}$ T each)") in tex
        assert "not estimated here" not in tex and "$4.6\\times 10^{4}$" not in tex
    _is(i["fq_light_nucleus_over_cap"], 100)              # "a factor ~1e2 over the 1e5 cap"
    _rng(r28.epsilon_l, (0.1 / r28.hard_ops[1], 0.1 / r28.hard_ops[0]), 1e-12)  # the 2028 box prints no row

# ------------------------------------------------- R14 item 1a: the dipole factor

def test_dipole_oscillator_length(a):
    # b from r_ch(4He) = 1.67824 fm, r_p = 0.8409 fm: (A-1)/A (3/2) b^2 = r_ch^2 - r_p^2      app11:124
    o = ch3.he4_oscillator_length(a)
    _is(o["b"] ** 2 * 1.5 * 3 / 4, 1.67824 ** 2 - 0.8409 ** 2)
    assert round(o["b"], 2) == 1.37 and round(o["hbar_omega"]) == 22 and round(o["qb"], 2) == 0.73
    _is(o["q"], 105.6583755 / 197.3269804)


def test_dipole_radial_matrix_elements_match_closed_forms(a):
    d = ch3.he4_dipole(a)
    printed = {(0, 0, 0): 0.874, (0, 0, 1): 0.796, (1, 1, 0): 0.728, (0, 0, 2): 0.722, (0, 1, 0): 0.096}
    for key, (quad, closed) in d["radial_me"].items():
        assert math.isclose(quad, closed, rel_tol=1e-9), (key, quad, closed)
        assert round(quad, 3) == printed[key], key
    # q -> 0: j0 -> 1, orthonormality
    for n1, n2, l, want in ((0, 0, 0, 1), (0, 1, 0, 0), (1, 1, 0, 1), (0, 0, 2, 1)):
        assert abs(ch3.j0_radial_me(n1, n2, l, 1e-6, 1.37) - want) < 1e-9


def test_dipole_operator_structure(a):
    import numpy as np
    o, labels, _ = ch3.he4_dipole_one_body(a)
    assert o.shape == (40, 40) and len(labels) == 40 and len(ch3.HE4_SPATIAL) == 10
    assert np.allclose(o, o.T)
    neutrons = [p for p, lab in enumerate(labels) if lab[3] == "n"]
    assert not o[neutrons, :].any() and not o[:, neutrons].any()      # protons only
    off = [(p, q) for p in range(40) for q in range(p + 1, 40) if abs(o[p, q]) > 1e-12]
    assert len(off) == 2                                              # 0s-1s, spin up and down


def test_jw_decomposition_is_exact_on_a_small_instance():
    # the decomposition used for Gamma and lambda reproduces sum o_pq a_p^dag a_q as a matrix
    import numpy as np
    rng = np.random.default_rng(7)
    n = 4
    o = rng.normal(size=(n, n)); o = o + o.T
    X = np.array([[0, 1], [1, 0]]); Y = np.array([[0, -1j], [1j, 0]]); Z = np.diag([1, -1]); I2 = np.eye(2)
    def kron(ops):
        out = np.array([[1.0]])
        for m in ops:
            out = np.kron(out, m)
        return out
    # a_j = Z^{<j} |0><1|_j (occupied = |1>, Z = -1)
    a_ops = [kron([Z] * j + [np.array([[0, 1], [0, 0]])] + [I2] * (n - j - 1)) for j in range(n)]
    target = sum(o[p, q] * a_ops[p].conj().T @ a_ops[q] for p in range(n) for q in range(n))
    P = {"X": X, "Y": Y, "Z": Z}
    built = np.zeros((2 ** n, 2 ** n), dtype=complex)
    for key, c in ch3.jw_one_body_paulis(o).items():
        ops = [I2] * n
        for qb, l in key:
            ops[qb] = P[l]
        built += c * kron(ops)
    assert np.allclose(built, target)


def test_dipole_lambda_expectation_and_gamma(a, r28):
    d = ch3.he4_dipole(a)
    # JW strings: identity + 20 Z (proton modes) + 4 XX/YY (0s-1s, two spins)
    assert d["gamma"] == 25 and d["gamma_counts"] == dict(identity=1, z=20, xx_yy=4)
    # R15 ruling 1: the identity is added back classically; the LCU has 24 strings
    assert d["gamma_without_identity"] == d["gamma_lcu"] == 24 and d["gamma_lcu"] <= 128
    # lambda (all strings) = Tr(o) + 2 |o_{0s,1s}|; the identity coefficient c_0 = Tr(o)/2 = 7.60
    _is(d["lam"], d["trace_o"] + 2 * abs(d["radial_me"][(0, 1, 0)][0]))
    _is(d["identity_coefficient"], d["trace_o"] / 2)
    assert round(d["identity_coefficient"], 2) == 7.60
    _is(d["lam_lcu"], d["lam"] - d["identity_coefficient"])
    assert round(d["lam"], 1) == 15.4 and round(d["lam_lcu"], 2) == 7.79 and round(d["exp_closed_shell"], 2) == 1.75
    _is(d["exp_closed_shell"], 2 * d["radial_me"][(0, 0, 0)][0])
    # the estimator: the test returns (<O> - c_0)/lambda' with per-shot variance <= 1; relative accuracy eps
    # on <O> needs (lambda'/(eps |<O>|))^2 shots
    _is(d["ratio"], d["lam_lcu"] / d["exp_closed_shell"])
    assert round(d["ratio"], 1) == 4.5 and ch3.round_sig(d["shot_factor"]) == 20
    # the rejected reading lambda'/|<O> - c_0| (the R14 log's 1.77) sets accuracy on <O> - c_0, not on <O>
    assert round(d["shot_factor_over_shifted_mean"], 2) == 1.77
    # any 2-proton state: <O> in [2 e_min, 2 e_max] of the proton block -> factor 18-33
    lo, hi = d["exp_band_any_state"]
    assert (round(lo, 2), round(hi, 2)) == (1.36, 1.84) and lo <= d["exp_closed_shell"] <= hi
    assert tuple(ch3.round_sig(x) for x in d["shot_factor_band"]) == (18, 33)
    # the test's mean -0.75, Bernoulli variance 0.44 (not credited)
    assert round(d["test_mean"], 2) == -0.75 and round(d["bernoulli_variance"], 2) == 0.44
    # the R14 record with the identity kept: 77 (70-130)
    assert ch3.round_sig(d["shot_factor_with_identity"]) == 77
    assert tuple(ch3.round_sig(x) for x in d["shot_factor_band_with_identity"]) == (70, 130)
    # the dipole's own insertion: k = 5, 62 rotations, 23 Toffolis, 1.2e3 T
    c = r28.intermediates["dipole_insertion_counts"]
    assert (c["k"], c["n_rot"], c["n_toffoli"]) == (5, 62, 23) and ch3.round_sig(c["t"]) == 1200



def test_dipole_factor_as_printed(tex):
    if tex is None:
        pytest.skip("chapter .tex not on disk")
    # E30 cap pass: the oscillator-length, c_0, closed-shell 1.75 / 4.5 / 20 and 0.44-variance steps of the
    # derivation are no longer printed; the model pins them (al27/he4 dipole tests above)
    # style pass 2026-10-08: retired, the sentence was cut (rule 1/2); the value stays pinned on the model
    for frag in ("it is the identity plus $\\Gamma=24$ Pauli strings with $\\lambda=7.79$",
                 "is $18$--$33$, and we count shots at $33$",
                 "$3.3\\times 10^{3}$ shots per operator, or $2.0\\times 10^{4}$ over the six categories",
                 "after its identity term, a classically known constant, has been subtracted",
                 "the dipole's own $24$ strings cost $1.2\\times 10^{3}$ T"):
        assert frag in tex, frag
    assert "=77$" not in tex and "$70$--$130$" not in tex

# ------------------------------------------------- R14 item 1d: depth-capped MLAE

def test_mlae_depth_is_self_consistent(a, r28):
    i = r28.intermediates
    lo, hi = i["mlae"]                                        # D = 1e2 and the last fit, 596 (M3 ruling)
    assert (lo["d_amp"], hi["d_amp"]) == (100.0, 596)
    assert i["mlae_m_max"] == (5, 1)
    for x in (lo, hi):
        m = x["m_max"]
        assert m % 2 == 1 and x["grover_iterates"] == (m - 1) // 2
        assert x["t_deepest"] <= 1e5 < x["t_next_depth"]     # m fits, m + 2 (its own budget) does not
        # rotations synthesized to the deepest circuit's budget: N_rot = m x 254
        assert x["n_rot_budget"] == m * 254
        _is(x["eps_rot"], math.sqrt(1e-2 / (m * 254)))
        # T(m) = m (7 x (SOS(D) + angle adds) + T_query) + (m-1)/2 x 56 x 7 + the phase-gradient state once
        # (one complete load with its angle rotations per application; M3 ruling)
        ang = ch3.load_angle_cost(a, x["d_amp"])
        ct = c.t_per_rotation(c.eps_rot_for(ang["n_rot"] + m * 254, 1e-2))
        _is(x["eps_rot"], math.sqrt(1e-2 / (m * 254)))
        tq = 254 * ct + 7 * 127
        assert x["reflection_qubits"] == 40 + 10 + 7 + 1 and x["reflection_ands"] == 56
        _is(x["t_deepest"], m * (7 * (ch3.sos_load_toffoli(x["d_amp"]) + ang["toffoli"]) + tq)
            + (m - 1) // 2 * 56 * 7 + ang["n_rot"] * ct + 1, 1e-4)
    assert lo["iterations"][-1][0] == lo["iterations"][-1][2] == 5          # fixed point at D = 1e2
    assert ch3.round_sig(lo["t_deepest"]) == 8.7e4
    assert ch3.round_sig(i["mlae_no_angles"][0]["t_deepest"]) == 8.4e4      # record without the angle term
    # m_max = 5 up to D = 128, 3 from 129 to 221, 1 from 222 ("3 above D = 128, the plain test above D = 221")
    assert i["mlae_m_max_first_d"] == {5: 100, 3: 129, 1: 222}
    assert i["mlae_m_max_last_d"] == {5: 128, 3: 221, 1: 596}
    assert i["mlae_m_max_last_d_no_angles"] == {5: 128, 3: 228, 1: 603}
    # a boundary where the naive iteration 2-cycles (m = 3 fits at the m = 1 budget, not its own): m_max = 1
    x = ch3.mlae_2028(a, 128, 229)
    assert x["m_max"] == 1 and x["t_next_depth"] > 1e5
    # with the dipole's own Gamma = 24 the depth at D = 1e2 would be 7 (sensitivity, not printed)
    assert tuple(x["m_max"] for x in i["mlae_dipole_gamma"]) == (7, 1)
    # the R14 lookup-only record at D = 1e2: m_max = 13
    assert i["lookup_r14"]["mlae_at_d_lo"]["m_max"] == 13



def test_mlae_circuits_and_wall(a, r28):
    i = r28.intermediates
    assert i["mlae_estimator_constant"] == 1.0
    for n, x in zip(i["mlae_circuits_per_operator"], i["mlae"]):
        _is(n, i["dipole_shot_factor"] * 100 / x["m_max"] ** 2)
    lo, hi = i["mlae_circuits"]
    assert ch3.round_sig(lo) == 790 and ch3.round_sig(hi) == 20000      # D = 600: the plain test
    w = i["mlae_wall_time_s"]
    assert round(w[0]) == 68                                             # "~68 s (MLAE at D = 1e2)" (M3 ruling; was 66)
    assert round(i["mlae_t_circuit_s"][0], 3) == 0.087                   # "~0.087 s per amplitude-estimation circuit"
    # r25: F* = 10.9 >= 10 (compiled depth 7942), so the T supply sets the circuit time and the wall stands
    assert round(i["mlae_f_star"], 1) == 10.9 and i["mlae_wall_time_corrected_s"] == w[0]
    # the aliasing window stays in one monotone branch of sin^2(m theta), so no shallow circuits
    assert i["mlae_window_in_one_branch"] == (True, True)
    assert round(i["mlae_branch_margin_rad"][0], 3) == 0.035             # "0.035 rad from its edge"
    # Suzuki LIS k = 0..K would need 2.1x more circuits at m = 5 (sensitivity)
    assert round(i["mlae_lis_constant"][0], 1) == 2.1
    # eps_l at the deepest circuit (the 2028 box prints no eps_l row)
    assert all(0.99e-6 < e < 1.25e-6 for e in i["mlae_epsilon_l_deepest"])



def test_mlae_as_printed(tex):
    if tex is None:
        pytest.skip("chapter .tex not on disk")
    # E30 cap pass: the Grover-iterate and reflection-register itemization is condensed to its totals
    # style pass 2026-10-08: retired, the sentence was cut (rule 1/2); the value stays pinned on the model
    for frag in ("give $m_{\\max}=5$ at "
                 "$D=10^{2}$ (longest circuit $8.7\\times 10^{4}$ T), $3$ above $D=128$, and the plain test above $D=221$",
                 "keeps it there at $m_{\\max}=5$, $0.035$ rad from the edge",
                 "$(\\lambda/|\\langle\\hat{\\mathcal{O}}\\rangle|)^2(\\epsilon\\,m_{\\max})^{-2}$ circuits per operator",
                 "\\cite{arxiv_1904_10246,arxiv_2012_03348}"):
        assert frag in tex, frag
    assert "$m_{\\max}=13$" not in tex and "0.005 rad" not in tex

# ------------------------------------------------- R14 item 2: option A, 0nubb by projection

def test_0nubb_projection_route(a, r33, rcd, tex):
    i = r33.intermediates
    # N_proj = 2 (parent preparation + daughter projection) is what the survey already prices:
    # each pair's T equals 2 x N_Trotter x N_orb^4 x c_T at that doubled rotation count
    assert i["n_proj_pairs"] == 2 and i["n_registers_pairs"] == 1
    nt = i["n_trotter"]
    n = 2 * nt[1] * 11 ** 4
    _is(i["prep_jj44"][1], n * c.t_per_rotation(c.eps_rot_for(n)))
    # the projection readout: 25 (lambda/|M|)^2 runs per element at eps = 0.1
    _is(i["shots_0nubb_per_lambda_over_m_sq"], 25.0)
    if tex is not None:
        assert "QPE-based overlap estimator" not in tex and "algorithmically equivalent alternative" not in tex
        assert "with probability $|M/\\lambda|^2$" in tex
        assert "$25\\,(\\lambda/|M|)^2$ at $\\epsilon=0.1$" in tex
        assert "$N_{\\rm proj}=2$ for the $\\beta\\beta$ pairs and $1$ for $\\mu\\to e$" in tex     # style pass: wording
        assert "N_{\\rm reg}" not in tex and "two registers" not in tex and "held at once" not in tex
        assert "4\\,N_{\\rm orb} \\;+\\; N_{\\rm anc} & \\text{(JW second quantization, one register" in tex


# ------------------------------------------------------------- 2033 register

def test_2033_registers(a, r33):
    i = r33.intermediates
    assert i["system_qubits"] == 24               # "~24 system qubits for 27Al (sd shell)"  app11:134
    # R14 option A: one register per pair (was 80, 88, 128, 216 for two)
    assert i["register_pf"] == 40                 # "~40 for the 48Ca/48Ti pair"             app11:134
    assert i["register_jj44"] == 44               # "~44 for the 76Ge/76Se pair"             app11:134
    assert i["register_jj55"] == 64
    # R15 ruling 3: Mo/Ru = jj44 protons (22) + jj55 neutrons (32) = 54 modes, N_orb 13.5 (was 27, 108)
    assert i["register_mo_ru_modes"] == (22, 32)
    assert i["n_orb_mo_ru"] == 13.5 and i["register_mo_ru"] == 54
    # per-species state counts behind every single-register cell: sd 12, pf 20, jj44 22, jj55 32
    assert {k: len(ch3.shell_states(k)) for k in ch3.SHELLS} == dict(sd=12, pf=20, jj44=22, jj55=32,
                                                                     f7p3=12, f7p3p1=14)    # r25 R7 ladder
    for key, n_orb in (("sd", a.n_orb_sd.value), ("pf", a.n_orb_pf.value), ("jj44", a.n_orb_jj44.value),
                       ("jj55", a.n_orb_jj55.value)):
        assert 2 * len(ch3.shell_states(key)) == 4 * n_orb, key
    assert i["n_anc"] == (50, 150)                # "N_anc ~ 50-150"                         app11:85
    assert i["lq_sum"] == (74, 174)               # box "24 system + 50-150 ancilla"
    assert r33.lq == (74, 174) and r33.lq[1] <= 250
    assert i["survey_lq_max"] == 194 <= 200       # "These three cells fit below ~200 LQ"     app11:134



def test_2033_heavier_pairs_lq(rcd, tex):
    # "The heavier jj55 pairs (64 system qubits) and the 100Mo/100Ru pair (54 ...) sit at ~100-210 LQ"
    i = rcd.intermediates
    assert i["lq_jj55"] == (114, 214) and i["lq_mo_ru"] == (104, 204)
    assert min(i["lq_mo_ru"][0], i["lq_jj55"][0]) == 104 and max(i["lq_mo_ru"][1], i["lq_jj55"][1]) == 214
    if tex is not None:
        assert "($54$: protons in $jj44$, neutrons in $jj55$) sit at $\\sim 100$--$210$ LQ" in tex
        assert "11 + 16 = 27" not in tex and "$108$ system" not in tex

# ------------------------------------------------------------ 2033 hard ops

def test_2033_al_sd_cell(r33):
    i = r33.intermediates
    _rng(i["prep_al_sd"], (1.0e8, 4.4e8), 0.05)                    # "1.0-4.4e8 for 27Al (sd)"
    _rng(i["fraction_of_envelope_prep"], (0.10, 0.44), 0.05)       # "0.10-0.44x ... on state prep"
    # R15 ruling 4: one controlled query, Gamma = 12900 (selection rules) / 49140 (dense), in the prep circuit
    assert i["gamma_ins_sd"] == (12900, 49140)
    lo, hi = i["insertion_sd_counts"]
    assert (lo["k"], lo["n_rot"], lo["n_toffoli"]) == (14, 2 * 16383, 12899)
    assert (hi["k"], hi["n_rot"], hi["n_toffoli"]) == (16, 2 * 65535, 49139)
    _is(lo["n_rot_circuit"], ch3.survey(ch3.Assumptions())["al_sd"]["n_rot"][0] + 2 * 16383)
    _is(lo["eps_rot"], math.sqrt(1e-2 / lo["n_rot_circuit"]))
    _rng(i["insertion_sd"], (9.3e5, 3.9e6), 0.02)                  # "9.3e5-3.9e6 T"
    _rng(i["t_per_shot"], (1.1e8, 4.4e8), 0.05)                    # box "1.1e8-4.4e8 T, prep + insertion"
    _rng(i["t_per_shot"], (1.0537e8, 4.4036e8), 1e-3)              # R14 1.092e8-4.771e8
    _rng(i["fraction_of_envelope_total"], (0.11, 0.44), 0.05)      # "= 0.11-0.44x the box"
    assert r33.hard_ops == i["t_per_shot"]
    _is(r33.breakdown_total(), i["t_per_shot"][0])
    _rng(i["eps_l"], (2.3e-10, 9.5e-10), 0.02)                     # "eps_l required 2.3-9.5e-10"
    assert r33.epsilon_l == i["eps_l"]



def test_2033_valence_pairs(r33):
    i = r33.intermediates
    _rng(i["prep_pf"], (1.754e9, 7.305e9), 1e-3)         # "1.8-7.3e9 for 48Ca/48Ti (pf)"      app11:139
    _rng(i["prep_jj44"], (2.597e9, 1.081e10), 1e-3)      # "2.6e9-1.1e10 for 76Ge/76Se"        app11:139
    _rng(i["multiple_pf"], (1.8, 7.3), 0.03)             # "x1.8-7.3 for 48Ca/48Ti"            app11:147
    _rng(i["multiple_jj44"], (2.6, 11), 0.03)            # "x2.6-11 for 76Ge/76Se and 82Se/82Kr"
    _rng(i["near_miss_multiple"], (1.8, 7.3), 0.03)      # R16: pf alone, box near-miss row    app11:45,217
    assert i["near_miss_multiple"] == i["multiple_pf"]
    assert i["codesign_multiple_jj44"] == i["multiple_jj44"]   # R16: jj44 relabelled co-design
    _rng(i["eps_l_pf_jj44"], (9.2e-12, 5.7e-11), 0.02)   # "9.2e-12-5.7e-11 (pf/jj44)"         app11:240
    assert i["multiple_pf"][0] > 1                        # "one cell fits across its full band and the rest exceed it"
    # the pair carries "the same factor of two" as the others (same N_rot per register, doubled)
    n = i["n_trotter"][0] * 10 ** 4
    _is(i["prep_pf"][0], 2 * n * c.t_per_rotation(c.eps_rot_for(2 * n)))
    # the exact products, each end at its own eps: 7.305e9 and 1.0812e10, printed 7.3e9 / 1.1e10
    n = 2 * 4000 * math.pi * 11 ** 4
    _is(i["prep_jj44"][1], n * c.t_per_rotation(c.eps_rot_for(n)))
    assert ch3.round_sig(i["prep_pf"][1]) == 7.3e9 and ch3.round_sig(i["prep_jj44"][1]) == 1.1e10


def test_2033_jj44_with_insertion_crosses_ten(r33):
    # the sd headline is prep + insertion; the same convention on the jj44 pair (derived insertion, R15)
    i = r33.intermediates
    s = ch3.survey(ch3.Assumptions())["jj44"]
    for end in (0, 1):
        q = s["ins_counts"][end]
        _is(i["prep_plus_insertion_jj44"][end], q["n_rot_circuit"] * q["c_T"] + q["t_toffoli"])
        assert i["prep_plus_insertion_jj44"][end] > i["prep_jj44"][end]
    assert i["multiple_jj44_with_insertion"][1] > 10       # crosses the >~10x co-design line of main-overview:106
    assert i["multiple_jj44"][1] > 10                        # R-TOL: prep alone 10.8x; R16 ruling: co-design



def test_2033_operator_insertion_is_subdominant(r33):
    # R15 ruling 4: one controlled BE query of the 0nubb operator in the jj44 register (derived here)
    i = r33.intermediates
    assert i["terms_ins_jj44"] == (2501, 53361) and 53361 == math.comb(22, 2) ** 2
    assert i["gamma_ins_jj44"] == (16 * 2501, 16 * 53361) == (40016, 853776)
    lo = i["insertion_jj44_counts"][0]
    assert (lo["k"], lo["n_rot"], lo["n_toffoli"]) == (16, 2 * 65535, 40015)
    _rng(i["insertion_jj44"], (3.98e6, 4.13e6), 1e-3)            # "4.0-4.1e6 T"
    _rng(i["insertion_jj44_dense"], (6.52e7, 6.76e7), 1e-3)      # "6.5-6.8e7 T"
    assert i["insertion_fraction_of_jj44_prep"][1] < 0.05        # "sub-dominant but not zero"
    assert i["insertion_sd"][1] < i["prep_al_sd"][0]
    lo, hi = i["insertion_fraction_of_jj44_prep"]                # "0.04-0.2% of the jj44 pair preparation"
    assert (ch3.round_sig(100 * lo, 1), ch3.round_sig(100 * hi, 1)) == (0.04, 0.2)
    lo, hi = i["insertion_fraction_near_miss"]                   # box "(the insertion adds 0.06-0.2%)", pf alone (R16)
    assert (ch3.round_sig(100 * lo, 1), ch3.round_sig(100 * hi, 1)) == (0.06, 0.2)
    lo, hi = i["insertion_fraction_pf_jj44"]                     # main-overview:71 "all but 0.04-0.2%"
    assert (ch3.round_sig(100 * lo, 1), ch3.round_sig(100 * hi, 1)) == (0.04, 0.2)
    # the other pairs (not printed): pf 62208, jj55 158592, Mo/Ru 78176 strings
    ins = i["insertions"]
    assert (ins["pf"]["gamma"]["gamma"], ins["jj55"]["gamma"]["gamma"], ins["mo_ru"]["gamma"]["gamma"]) == \
        (62208, 158592, 78176)
    assert all(v["share"][0] < 0.01 for v in ins.values())



def test_2033_jj55_codesign_pairs(r33):
    i = r33.intermediates
    _rng(i["multiple_jj55"], (12, 50), 0.02)                 # "exceed the envelope by x12-50"      app11:150
    _rng(i["multiple_jj55_large_gap"], (6, 12), 0.02)        # "still by x6-12 at ... Delta = 4 MeV"
    n = 2 * (math.pi * 250 / 0.5) * 16 ** 4
    _is(i["multiple_jj55_large_gap"][0], n * c.t_per_rotation(c.eps_rot_for(n)) / 1e9)   # 5.95 exact
    _rng(i["eps_l_jj55"], (2.0e-12, 8.2e-12), 0.03)          # "2.0-8.2e-12 (jj55)"                 app11:240


def test_2033_mo_ru_generic_gap(r33):
    # R15 ruling 3: 54 modes, N_orb^4 = 13.5^4 = 27^4 / 16
    i = r33.intermediates
    _rng(i["multiple_mo_ru"], (6.0, 25), 0.02)                # "x6-25 over the envelope at the generic gap band"
    _rng(i["codesign_multiple"], (12, 50), 0.02)              # box co-design row, jj55 alone: "x12-50"
    _rng(i["multiple_mo_ru_large_gap"], (3.0, 6.0), 0.02)     # "x3-6 for a large gap"
    assert i["multiple_mo_ru"][1] < i["multiple_jj55"][1]                       # cheaper than a jj55 pair
    assert i["multiple_mo_ru_large_gap"][1] < i["near_miss_multiple"][1]        # large gap: a near miss
    assert i["multiple_mo_ru"][0] < 10 < i["multiple_mo_ru"][1]                 # straddles the ~10x line
    _rng(ch3.prep_t(a_default(), 2, 13.5, (4.0, 4.0)), (2.957e9, 6.034e9), 1e-3)



def a_default():
    return ch3.Assumptions()


def test_2033_mo_ru_small_gap(r33):
    # "x25-51 at a Delta = 0.5 MeV small-gap end" (box "up to x51"); R11 (A) Stated gap
    i = r33.intermediates
    _rng(i["multiple_mo_ru_small_gap"], (25.1, 51.1), 0.01)
    assert i["multiple_mo_ru_small_gap"][0] > 10                  # small gap: co-design
    _rng(i["eps_l_mo_ru"], (2.0e-12, 3.4e-11), 0.03)          # "2.0e-12-3.4e-11 (100Mo/Ru)"



def test_2033_prep_lever(r33):
    i = r33.intermediates
    # "A preparation algorithm that improves ... by a factor of ten ... the top of the jj44 band ...
    #  sits just over it, by about 1%"  app11:163. The x10 shorter circuit gets its own eps (R-TOL).
    n = 2 * 4000 * math.pi / 10 * 11 ** 4
    _is(i["jj44_worst_after_lever"], n * c.t_per_rotation(c.eps_rot_for(n)) / 1e9)
    assert 1 < i["jj44_worst_after_lever"] < 1.02                     # just over (was 0.99 at c_T = 27)
    assert ch3.round_sig(100 * (i["jj44_worst_after_lever"] - 1), 1) == 1   # "by about 1%", two-figure rule kept
    assert i["jj44_worst_after_lever"] - 1 > 0.003                    # outside the ~0.3% overshoot ruling: an overrun
    assert i["jj44_worst_after_lever"] < i["multiple_jj44"][1] / 10   # "slightly more than that factor"
    # "The 100Mo estimate would fall to 0.28-4.8 times the limit ... (0.56-2.4 at the generic band)" (R15)
    _rng(i["mo_ru_after_lever"], (0.2758, 4.796), 1e-3)
    _rng(i["mo_ru_after_lever_generic"], (0.5635, 2.350), 1e-3)


def test_2033_shots_wall_time_campaign(r33):
    i = r33.intermediates
    assert i["shots_per_operator"] == (100, 400)           # "eps^-2 ~ 100-400 shots per operator"  app11:92
    assert i["shots_mu2e_bare"] == 6 * 400 == 2400
    assert i["shots_survey_per_scan"] == 6 * 400 + 15 * 100 == 3900     # record: one pass over the survey
    # times the derived 27Al sd dipole factor, identity subtracted (R15), upper end of the any-state band
    f = i["al27_dipole_shot_factor"]
    assert round(f, 2) == 1.45
    _is(i["shots_mu2e_per_scan"], 2400 * f)
    _rng(i["shots_with_basis_scans"], (4800 * f, 7200 * f), 1e-12)
    # r17: Result.shots spans the 2-3 basis-cutoff scans (R15: one scan to three, kept as a record)
    assert r33.shots == i["shots"] == i["shots_campaign"] and i["shots_pre_r25"] == i["shots_with_basis_scans"]
    _rng(i["shots_one_to_three"], (2400 * f, 7200 * f), 1e-12)
    assert ch3.round_sig(i["shots_r14"]) == 1.4e4                      # R14 record
    assert i["shots_box_reading_pre_r11"] == (2000, 8000)   # the retired "~2-8e3" = (100-400) x 20
    assert i["survey_elements_sum"] == 21                   # 6 + 5x3, quoted as "~20"
    assert 42 <= i["subroutine_calls"][0] <= 50 <= i["subroutine_calls"][1] <= 63   # "~50"
    # 1 us per T-gate on 1.054e8-4.404e8 T -> 105.4-440.4 s: "~1.8-7.3 min"
    assert round(i["t_shot_s"][0] / 60, 1) == 1.8 and round(i["t_shot_s"][1] / 60, 1) == 7.3
    _rng(i["t_shot_s"], (r33.hard_ops[0] / 1e6 + 1e-4, r33.hard_ops[1] / 1e6 + 1e-4), 1e-12)   # r25 + 0.1 ms
    # pre-r25 record: six categories at the dipole factor, 2-3 scans (8.5 d-1.8 mo)
    _rng(i["wall_time_s_pre_r25"], (4800 * f * i["t_shot_s"][0], 7200 * f * i["t_shot_s"][1]), 1e-12)
    _rng(i["wall_time_one_scan_s"], (2400 * f * i["t_shot_s"][0], 7200 * f * i["t_shot_s"][1]), 1e-12)
    # r25: Result.shots / wall_time_s are the campaign (SD element, direct readout, over p, 2-3 scans)
    assert r33.shots == i["shots_campaign"] and r33.wall_time_s == i["wall_campaign_s"]



def test_2033_shots_as_printed(a, r33, tex):
    """r25 (shot audit, R1/R2 and the Ch. 3 rulings of 2026-10-02): the SD element sets the shots."""
    i = r33.intermediates
    d = ch3.al27_dipole_sd(a)
    # SI elements pinned classically: the any-state band is +/-0.61%
    assert round(100 * i["si_pinned_fraction"], 2) == 0.61
    # <O_SD> = 2 <S_p> o_0d at <S_p> = 0.34 .. 0.30
    _rng(i["o_sd"], (2 * 0.34 * d["o_0d"], 2 * 0.30 * d["o_0d"]), 1e-12)
    assert tuple(round(x, 2) for x in i["o_sd"]) == (0.38, 0.34)
    # Hadamard-test factor (lambda'/<O_SD>)^2 = 79-101 (the audit's 78.6-101.0)
    assert tuple(round(x) for x in i["sd_hadamard_factor"]) == (79, 101)
    # direct readout: var / (eps <O_SD>)^2 = 2.77e3-3.56e3, over p = 0.5: 5.54e3-7.11e3 per Hamiltonian
    _rng(i["sd_shots_per_hamiltonian_p1"], tuple(1 / (0.05 * o) ** 2 for o in i["o_sd"]), 1e-12)
    assert tuple(round(x) for x in i["sd_shots_per_hamiltonian_p1"]) == (2769, 3557)
    _rng(i["sd_shots_per_hamiltonian"], _div(i["sd_shots_per_hamiltonian_p1"], 0.5), 1e-12)
    _rng(i["shots_first_result"], i["sd_shots_per_hamiltonian"], 1e-12)
    _rng(i["shots_campaign"], (2 * i["sd_shots_per_hamiltonian"][0], 3 * i["sd_shots_per_hamiltonian"][1]), 1e-12)
    assert tuple(ch3.round_sig(x) for x in i["shots_campaign"]) == (1.1e4, 2.1e4)
    if tex is None:
        pytest.skip("chapter .tex not on disk")
    # style pass 2026-10-08: box rows hold quantity and value only (rule 10); numbers unchanged
    assert ("$5.5$--$7.1\\times 10^{3}$ per Hamiltonian; "
            "$1.1$--$2.1\\times 10^{4}$ over the $2$--$3$ basis-cutoff scans") in tex
    assert "Per-shot wall time & $\\sim 1.8$--$7.3$ min" in tex and "T-depth & $8.9\\times 10^{6}$--$3.7\\times 10^{7}$" in tex
    assert "Flagship $^{27}$Al & $e_{\\max}{=}8$--$10$: $\\approx 700$--$1100$ LQ, $\\times 8.5\\times 10^{4}$--$3.3\\times 10^{6}$ the 2033 target" in tex   # box row; multiple_flagship is in retired
    assert "(all six categories at the dipole's shot factor)" not in tex and "estimand;" not in tex
    assert "per basis scan" not in tex
    assert "2$--$8\\times 10^{3}$ across the survey" not in tex
    assert "fixes the dipole operator to within $\\pm 0.6\\%$" in tex


def test_2033_tiers_depth_and_walls(a, r33, tex):
    """r25 (R1, R9, R10): first result and campaign on one machine; T-depth at G = 12 (factory.json)."""
    i = r33.intermediates
    day, month, yr = 86400, ch3.YEAR_S / 12, ch3.YEAR_S
    assert a.trotter_concurrency.value == 12
    # T-depth = N_rot c_T / G + insertion tail; G = 1 returns the T-count; factory.json G = 10: 1.07e7-4.5e7
    _rng(i["t_depth_gate_by_gate"], r33.hard_ops, 1e-9)
    assert tuple(ch3.round_sig(x, 3) for x in i["t_depth_g10"]) == (1.07e7, 4.46e7)
    assert tuple(ch3.round_sig(x) for x in i["t_depth_per_shot"]) == (8.9e6, 3.7e7)
    assert all(11.7 < f < 11.9 for f in i["f_star_per_end"])         # >= 10: the baseline holds
    # walls: shots x (max(N_T t_g, D_T t_r) + 0.1 ms); here the T supply binds
    for key, shots in (("wall_first_result_s", i["shots_first_result"]), ("wall_campaign_s", i["shots_campaign"])):
        _rng(i[key], (shots[0] * i["t_shot_s"][0], shots[1] * i["t_shot_s"][1]), 1e-12)
    assert (round(i["wall_first_result_s"][0] / day, 1), round(i["wall_first_result_s"][1] / day)) == (6.8, 36)
    assert round(i["wall_campaign_s"][0] / day) == 14 and round(i["wall_campaign_s"][1] / month, 1) == 3.6
    assert round(100 * i["wall_fraction_of_horizon"][1], 1) == 6.0 and i["wall_fits_horizon_on_one_machine"]
    # gate by gate (G = 1): campaign floor 135 d - 3.0 yr
    assert round(i["floor_wall_gate_by_gate_s"][0] / day) == 135
    assert round(i["floor_wall_gate_by_gate_s"][1] / yr, 1) == 3.0
    # p over 0.1-0.9: 7.5 d - 1.5 yr; R1's generic 30% on SD: 154-198 shots, 4.5-24 h
    assert round(i["wall_campaign_p_band_s"][0] / day, 1) == 7.5 and round(i["wall_campaign_p_band_s"][1] / yr, 1) == 1.5
    assert tuple(round(x) for x in i["shots_first_result_at_30pct"]) == (154, 198)
    assert round(i["wall_first_result_at_30pct_s"][0] / 3600, 1) == 4.5
    assert round(i["wall_first_result_at_30pct_s"][1] / 3600) == 24
    # exports (campaign run), via common.depth_exports
    ex = c.depth_exports(r33.hard_ops, i["t_depth_per_shot"], i["shots_campaign"], a.t_gate_s, a.shot_overhead_s)
    for k in ("t_depth_per_shot", "f_star", "floor_wall_s", "factories_for_1yr"):
        assert i[k] == ex[k], k
    assert tuple(round(x, 2) for x in i["factories_for_1yr"]) == (0.37, 2.98)
    if tex is None:
        pytest.skip("chapter .tex not on disk")
    for frag in ("one Hamiltonian, $5\\%$ on the spin-dependent element: $\\approx 6.8$--$36$ days on one machine",
                 "$\\approx 14$ days--$3.6$ months on one machine ($2$--$3$ scans)",   # referee M2: SD element only
                 "takes at most $6.0\\%$ of the five-year window on one machine",     # style pass: rule 7 wording
                 "no number of factories brings the campaign below $135$ days to $3.0$ years",
                 "between $7.5$ days and $1.5$ years"):     # style pass: the 30% alternative (4.5-24 h) is cut (rule 3)
        assert frag in tex, frag
    assert "8.5$ days" not in tex and "$2.9\\%$" not in tex


def test_simplified_0nubb_rows(a, r33, tex):
    """r25 R7 (simplobs.json Ch. 3, reproduced): 48Ca first light in f7/2 p3/2, then f7/2 p3/2 p1/2 by Rabi,
    full pf by Rabi a near miss; 76Ge past 2033."""
    sim = r33.intermediates["simplified_0nubb"]
    day, yr = 86400, ch3.YEAR_S
    b, bp, h, apf = sim["f7p3_projection"], sim["f7p3p1_rabi"], sim["pf_rabi"], sim["pf_projection_0plus"]
    assert b["lq"] == (74, 174) and bp["lq"] == (78, 178) and h["lq"] == (90, 190)
    assert b["gamma"] == 7104 and apf["gamma"] == 62208
    # referee M3 (2026-10-05): each space at its own KB3G 0+ gaps (was ENSDF: 3563.3-7126.7 steps, 1.194-2.436e8 T,
    # 20-40 d; kept as t_at_ensdf_gaps)
    assert b["gaps"] == (6.100132, 3.684720) and apf["gaps"] == (5.174846, 4.264629)
    _rng(b["n_trotter"], (2735.2, 5470.4), 1e-4)
    _rng(b["t"], (9.098e7, 1.856e8), 2e-3)
    _rng(b["t_at_ensdf_gaps"], (1.194e8, 2.436e8), 2e-3)
    assert round(b["shots"][0]) == 14262 and round(b["wall_s"][0] / day) == 15 and round(b["wall_s"][1] / day) == 31
    _rng(bp["t"], (4.399e8, 7.200e8), 2e-3)
    assert round(bp["n_qdrift"][0], -4) == 1.33e6
    # E29 (verifier, R10): the qDRIFT samples run one at a time, D_T = N_trot c_T / 12 + N_q c_T
    for r_ in (b, bp, h):
        for e in (0, 1):
            assert r_["wall_s"][e] >= r_["wall_t_supply_s"][e] * (1 - 1e-12)
            assert r_["f_star_per_end"][e] == pytest.approx(r_["t"][e] / r_["t_depth"][e])
    _rng(bp["t_depth"], (6.9e7, 9.3e7), 1e-2)
    assert (round(bp["f_star_per_end"][0], 1), round(bp["f_star_per_end"][1], 1)) == (6.4, 7.7)
    assert (round(bp["wall_s"][0] / day, 1), round(bp["wall_s"][1] / day)) == (2.0, 11)
    assert (round(bp["wall_t_supply_s"][0] / day, 1), round(bp["wall_t_supply_s"][1] / day, 1)) == (1.3, 8.3)
    _rng(h["t"], (2.267e9, 3.716e9), 2e-3)
    _rng(h["t_depth"], (6.5e8, 9.8e8), 1e-2)
    assert (round(h["f_star_per_end"][0], 1), round(h["f_star_per_end"][1], 1)) == (3.5, 3.8)
    assert (round(h["wall_s"][0] / day), round(h["wall_s"][1] / day, -1)) == (19, 110)
    assert b["f_star_per_end"][0] > 10 and b["wall_s"] == b["wall_t_supply_s"]   # 12 in flight, F* ~ 12
    _rng(b["eps_l"], (5.39e-10, 1.10e-9), 1e-2)
    _rng(bp["eps_l"], (1.4e-10, 2.3e-10), 2e-2)
    _rng(h["eps_l"], (2.7e-11, 4.4e-11), 1e-2)
    assert 1.5 < h["multiple"][0] and h["multiple"][1] < 10        # a near miss, not co-design
    _rng(apf["t"], (7.354e8, 1.498e9), 2e-3)
    _rng(apf["t_at_ensdf_gaps"], (9.821e8, 2.001e9), 2e-3)
    assert apf["shots_lower_bound"] == (1.0e6, 1.44e6)
    assert (round(apf["wall_s"][0] / yr), round(apf["wall_s"][1] / yr)) == (23, 68)
    assert all(g < e for g, e in zip(bp["projection_gaps"], (6.100132, 3.684720)))   # p1/2 lowers both gaps
    if tex is None:
        pytest.skip("chapter .tex not on disk")
    # style pass 2026-10-08: the truncated-space gaps and the generic-band comparison are cut (rule 2); box rows
    # follow the vocabulary of rule 12 (readout -> estimate, resonant drive -> resonant transition); numbers unchanged
    for frag in ("are $5.17$ and $4.26$ MeV ($4.28$ and $3.00$ MeV measured~\\cite{ENSDF})",
                 "so KB3G is the conservative choice", "\\cite{arxiv_nucl-th_0012077}",
                 "$74$--$174$ LQ, $0.91$--$1.9\\times 10^{8}$ T, projection estimate: $15$--$31$ days",
                 "$78$--$178$ LQ, $4.4$--$7.2\\times 10^{8}$ T, resonant-transition estimate: $2.0$--$11$ days",
                 "full $pf$ by resonant transition $\\times 2.3$--$3.7$, $19$--$110$ days",
                 "($\\times 0.74$--$1.5$ at the model $0^+$ gaps)",
                 "The qDRIFT samples run one at a time, so these circuits keep only $4$--$8$ factories busy.",
                 "accidental resonances between excited parent and daughter $0^+$ states",
                 "$1.4\\times 10^{-10}$--$1.1\\times 10^{-9}$ ($0\\nu\\beta\\beta$ validation stage)",
                 "with per-shot variance about one (our estimate; the worst case is about 8 times larger)",
                 "\\cite{Campbell_qDRIFT_2019}", "It stays past 2033.",
                 "We take the daughter projection at the parent's cost; a daughter with a smaller gap raises it"):
        assert frag in tex, frag


def test_single_machine_serial_campaign(a, r28, r33, tex):
    """r23 (author ruling 2026-10-02): one quantum machine everywhere. Every wall is shots x per-shot time."""
    import dataclasses
    # the model carries no machine count, in the inputs or in any era's intermediates
    banned = ("machine", "parallel")
    assert not [f.name for f in dataclasses.fields(a) if any(b in f.name for b in banned)]
    for r in (r28, r33):
        hits = [k for k in r.intermediates if any(b in k for b in ("n_machine", "machines", "parallel"))]
        assert not hits, hits
    # serial wall = shots x per-shot time, band end by band end, both eras
    # r25: per shot max(N_T t_g, D_T t_r) + t_0, so the wall is at least the T-supply wall, equal where F* >= 10
    for r in (r28, r33):
        t = r.intermediates["t_shot_s"]
        assert r.wall_time_s[0] >= r.shots[0] * t[0] * (1 - 1e-12) and r.wall_time_s[1] >= r.shots[1] * t[1] * (1 - 1e-12)
    _rng(r33.wall_time_s, (r33.shots[0] * r33.intermediates["t_shot_s"][0],
                           r33.shots[1] * r33.intermediates["t_shot_s"][1]), 1e-12)
    i28 = r28.intermediates
    _is(i28["mlae_wall_time_s"][0], i28["mlae_circuits"][0] * i28["mlae_t_circuit_s"][0])
    # the 2033 campaign against the 5-year horizon: 0.74-6.0%, fits on one machine, no reduction needed
    i = r33.intermediates
    assert a.campaign_horizon_yr.prov is c.Provenance.STATED and a.campaign_horizon_yr.value == 5
    _is(ch3.YEAR_S, 365.25 * 86400)
    _is(i["campaign_horizon_s"], 5 * 365.25 * 86400)
    _rng(i["wall_fraction_of_horizon"], _div(r33.wall_time_s, i["campaign_horizon_s"]), 1e-12)
    assert tuple(ch3.round_sig(100 * x) for x in i["wall_fraction_of_horizon"]) == (0.74, 6.0)
    assert tuple(ch3.round_sig(100 * x) for x in i["wall_fraction_of_horizon_pre_r25"]) == (0.47, 2.9)   # record
    assert i["wall_fits_horizon_on_one_machine"] is True and r33.wall_time_s[1] < i["campaign_horizon_s"]
    with pytest.raises(ValueError, match="campaign_horizon_yr"):
        ch3.Assumptions(campaign_horizon_yr=c.Stated((5, 10), "x"))
    if tex is None:
        pytest.skip("chapter .tex not on disk")
    assert "Campaign horizon & 5 years" in tex
    low = tex.lower()
    for s in ("shot-parallel", "shot parallel", "embarrassingly", "across machines", "machine-year",
              "machine year", "-fold parallelism", "serialized", "divides across", "machines needed",
              "every shot is independent", "shot throughput"):
        assert s not in low, s
    # the only 'parallel' left is the parallel Hadamard test (one circuit, one machine); cut pass 2026-10-06
    # reworded "closely paralleling Ch. 2" to "as in Ch. 2"
    import re
    assert sorted(re.findall(r"parallel\w*", low)) == ["parallel"]
    assert "the parallel Hadamard test~\\cite{arxiv_2607_12344}" in tex
    # 'machines' (plural) appears nowhere; 'machine' only as one machine / machine time / the RFI machine
    assert not re.search(r"\bmachines\b", low)
    assert "shot-parallel" not in low and "basis-cutoff scans \\\\" in tex



def test_2033_external_comparison_is_cited_not_derived(a, r33):
    # "~1e9 Toffoli gates in the sd shell" (arXiv:2607.21563) is quoted, not reproduced      app11:144
    assert a.ext_qubitization_sd_toffoli.prov is c.Provenance.CITED
    assert r33.intermediates["ext_qubitization_sd_toffoli"] == 1e9


# --------------------------------------------------------- co-design ladder

def test_codesign_emax3(rcd):
    i = rcd.intermediates
    assert i["register_emax3"] == 80                        # "approximately 80 system qubits"    app11:152
    _rng(i["t_emax3"], (1.49e10, 6.19e10), 0.005)           # "1.5-6.2e10 T operations" (1.490e10 exact)
    _rng(i["multiple_emax3"], (15, 62), 0.01)               # "exceeding the limit by a factor of 15-62"
    n = 4000 * math.pi * 20 ** 4
    _is(i["t_emax3"][1], n * c.t_per_rotation(c.eps_rot_for(n)))     # exact 6.19e10
    assert ch3.round_sig(i["t_emax3"][1]) == 6.2e10


def test_codesign_emax4(rcd):
    # 35 spatial orbitals at e_max = 4: x148.5-615.6 under R-TOL, printed "150-620" (app11:152, 156).
    # (x127-509 at c_T = 27; the older "150-600" was the e_max=3 pair scaled with N_orb = 36.)
    i = rcd.intermediates
    assert i["ho_orbitals"][4] == 35
    _rng(i["multiple_emax4"], (148.5, 615.6), 0.005)
    assert (ch3.round_sig(i["multiple_emax4"][0]), ch3.round_sig(i["multiple_emax4"][1])) == (150, 620)


def test_codesign_flagship(rcd):
    i = rcd.intermediates
    assert i["system_qubits"] == (660, 1144)                # "~700-1100 system qubits"          app11:135
    _rng(rcd.lq, (700, 1100), 0.10)
    assert i["lq_with_ancilla"] == (710, 1294)              # + N_anc 50-150
    assert i["fits_1000_lq_emax8"] and i["fits_1000_lq_emax8_with_ancilla"]   # 660 / 710 fit
    assert not i["fits_1000_lq_emax10"]                     # 1144 does not
    _rng(i["t_flagship"], (8.532e13, 3.331e15), 1e-3)       # "8.5e13-3.3e15 T operations per state preparation"
    n = 4000 * math.pi * 286 ** 4
    _is(i["t_flagship"][1], n * c.t_per_rotation(c.eps_rot_for(n)))   # exact 3.331e15 at 39.62 T/rot
    assert ch3.round_sig(i["t_flagship"][1]) == 3.3e15
    _rng(i["multiple_flagship"], (8.5e4, 3.3e6), 0.01)      # "x8.5e4-3.3e6"                     app11:152,215
    _rng(i["eps_l_flagship"], (3e-17, 1e-15), 0.25)        # "3e-17-1e-15 (flagship)", one figure  app11:240
    assert rcd.hard_ops == i["t_flagship"]
    _is(rcd.breakdown_total(), i["t_flagship"][0])


def test_codesign_flagship_fits_1000_lq_as_printed(rcd, tex):
    # R11 (A): "sits at the 1000-LQ envelope (inside it at e_max=8, ~14% over at e_max=10 before ancilla)"
    i = rcd.intermediates
    assert i["fits_1000_lq_emax8"] is True
    assert i["fits_1000_lq_emax10"] is False
    assert round(100 * (i["system_qubits"][1] / 1000 - 1)) == 14          # 1144 / 1000
    if tex is not None:
        assert "fits within approximately 1000 logical qubits" not in tex
        assert "$\\sim 14\\%$ over at $e_{\\max}=10$ before ancilla" in tex


def test_codesign_ge_converged_registers(rcd):
    i = rcd.intermediates
    assert i["ge_converged_orbitals"] == (286, 680)
    # R14 option A: one register, 1144-2720 ("~1100-2700 system qubits in the single register")  app11:141
    assert i["ge_converged_qubits"] == (1144, 2720)
    _rng(i["ge_converged_qubits"], (1100, 2700), 0.05)
    assert i["ge_converged_qubits"][0] > 1000               # "exceeds the 2033 LQ envelope" (JW)
    # first-quantized line of eq:Nq11: A ceil(log2 N_orb) = 684-760; "compresses ... by factors of a few"
    assert i["ge_first_quantized_qubits"] == (76 * 9, 76 * 10)
    lo, hi = i["ge_first_quantized_compression"]
    assert 1.6 <= lo <= 1.7 and 3.5 <= hi <= 3.6, (lo, hi)
    assert i["ge_first_quantized_qubits"][1] < 1000         # "fits the envelope in width"


# ------------------------------------------------------------- cross-checks

def test_all_cells_share_one_formula(a, r33, rcd):
    # every cell is eq:Ngate11 with N_reg, N_orb, gap changed; nothing else
    nt = r33.intermediates["n_trotter"]

    def f(n):          # N_rot x c_T(N_rot), R-TOL
        return n * c.t_per_rotation(c.eps_rot_for(n))
    for key, n_reg, n_orb in (("prep_al_sd", 1, 6), ("prep_pf", 2, 10), ("prep_jj44", 2, 11),
                              ("prep_jj55", 2, 16), ("prep_mo_ru", 2, 13.5)):
        _is(r33.intermediates[key][0], f(n_reg * nt[0] * n_orb ** 4))
        _is(r33.intermediates[key][1], f(n_reg * nt[1] * n_orb ** 4))
    _is(rcd.intermediates["t_emax3"][0], f(nt[0] * 20 ** 4))
    _is(rcd.intermediates["t_flagship"][0], f(nt[0] * 165 ** 4))


def test_band_ratio_is_four(r33):
    # Delta 1-2 MeV (x2) times c_proj pi-2pi (x2) = 4 in rotations across every prep band; under
    # R-TOL the top end also synthesizes tighter (eps / 2), +1.15 T per rotation
    for key in ("prep_al_sd", "prep_pf", "prep_jj44", "prep_jj55", "prep_mo_ru"):
        lo, hi = r33.intermediates[key]
        assert 4 < hi / lo < 4 * 1.05


def test_published_is_the_box_as_printed():
    assert ch3.PUBLISHED["2028"].lq == (120, 170)
    assert ch3.PUBLISHED["2028"].hard_ops == (1.7e4, 1.0e5)          # M3 ruling: + angle rotations; G2 1.6e4-9.9e4; R13 5.4e3
    assert ch3.PUBLISHED["2033"].lq == (74, 174)
    assert ch3.PUBLISHED["2033"].hard_ops == (1.1e8, 4.4e8)                # R15 derived insertion
    assert ch3.PUBLISHED["codesign"].lq == (700, 1100)
    assert ch3.PUBLISHED["codesign"].hard_ops == (8.5e13, 3.3e15)    # R-TOL (2026-09-29); was 6e13-2.3e15


# ---------------------------------------------------- printed strings (R4)

def test_round_once_at_print():
    assert ch3.round_sig(127.29) == 130 and ch3.round_sig(509.15) == 510
    assert ch3.round_sig(0.99351) == 0.99 and ch3.round_sig(6.7858e9) == 6.8e9
    assert ch3.print_range((45.078, 90.157)) == "45--90"
    assert ch3.print_range((5.559, 11.118)) == "5.6--11"
    assert ch3.print_range((5.95, 12.14)) == "6--12"
    assert ch3.print_range((4.508, 72.125)) == "4.5--72"


def test_printed_strings_are_in_the_tex(a, tex):
    """Every string printed(a) returns must appear verbatim in the chapter (exact carry, one rounding)."""
    p = ch3.printed(a)
    expected = {
        "worked_line:90": "$1.0$--$4.4\\times 10^{8}$",
        "prep_al_sd:139": "$1.0$--$4.4\\times 10^{8}$",
        "prep_pf:139": "$1.8$--$7.3\\times 10^{9}$",
        "prep_jj44:139": "$2.6\\times 10^{9}$--$1.1\\times 10^{10}$",
        "fraction_total_sd:147,210": "0.11--0.44",
        "multiple_pf:147": "1.8--7.3",
        "multiple_jj44:147": "2.6--11",
        "near_miss:45,217": "1.8--7.3",
        "multiple_jj55:150": "12--50",
        "jj55_large_gap:150": "6--12",
        "t_emax3:152": "$1.5$--$6.2\\times 10^{10}$",
        "multiple_emax3:152,156": "15--62",
        "multiple_emax4:152,156": "150--620",
        # R15 ruling 3: Mo/Ru on 54 modes
        "mo_ru_generic:46,162,219": "6--25",
        "mo_ru_large_gap:162": "3--6",
        "mo_ru_small_gap:162,219": "25--51",
        "codesign_multiple:218": "12--50",
        "codesign_jj44:45,152,218,244": "2.6--11",
        "lever_worst:168": "by about 1\\%",
        "pf_after_lever:168": "0.16--0.68",
        "mo_ru_after_lever:169": "0.28--4.8",
        "mo_ru_after_lever_generic:169": "0.56--2.4",
        "n_rot_flagship:152": "$2.3\\times 10^{12}$--$8.4\\times 10^{13}$",
        "t_flagship:152,239": "$8.5\\times 10^{13}$--$3.3\\times 10^{15}$",
        "multiple_flagship:152,215": "$8.5\\times 10^{4}$--$3.3\\times 10^{6}$",
        "eps_l_sd:205,240": "$2.3$--$9.5\\times 10^{-10}$",
        "eps_l_pf_jj44:240": "$9.2\\times 10^{-12}$--$5.7\\times 10^{-11}$",
        "eps_l_jj55:240": "$2.0$--$8.2\\times 10^{-12}$",
        "eps_l_mo_ru_row:246": "$2.0\\times 10^{-12}$--$3.4\\times 10^{-11}$",
        "eps_l_flagship:240": "$3\\times 10^{-17}$--$1\\times 10^{-15}$",
        # R-TOL (2026-09-29): the rule sentence and the worked line, app11:90; flagship, app11:152
        "eps_rot_sd:90": "$2.5$--$5.0\\times 10^{-5}$",
        "c_T_sd:90": "25.6$--$26.8",
        "eps_rot_survey:90": "$2.5\\times 10^{-6}$--$5.0\\times 10^{-5}$",
        "c_T_survey:90": "26--31",
        "c_T_coherent:90": "42--52",
        "eps_rot_flagship:152": "$1.1$--$6.6\\times 10^{-8}$",
        "c_T_flagship:152": "37--40",
        # R15 ruling 4: the derived insertions
        "insertion_share_jj44:143": "0.04$--$0.2\\%",
        "insertion_share_near_miss:217": "0.06$--$0.2\\%",
        "insertion_jj44:143": "$4.0$--$4.1\\times 10^{6}$",
        "insertion_jj44_dense:143": "$6.5$--$6.8\\times 10^{7}$",
        "insertion_sd:143": "$9.3\\times 10^{5}$--$3.9\\times 10^{6}$",
        "t_per_shot_sd:209,243": "$1.1\\times 10^{8}$--$4.4\\times 10^{8}$",
        # r23: the single-machine serial campaign as a share of the 5-year horizon
        "wall_fraction_of_horizon_max:149": "at most $6.0\\%$",
    }
    r25 = ch3.printed_r25(a)
    assert set(expected) | set(r25) == set(p), set(expected) ^ set(p)
    # style pass 2026-10-08 (rules 1, 2, 5): these strings are intermediates or restate a box row, and their
    # sentences were cut from the prose. Their values stay pinned above (p[k] == v) and in the model; only the
    # check that the chapter prints them is retired.
    retired = {"t_shot_supply_2028:127", "t_layers_2028:127", "f_star_d100_2028:127",
               "t_shot_corrected_d100_2028:127", "f_star_d600_2028:127", "mlae_f_star_2028:127", "walls_2028:127",
               "first_result_2033:151,box", "campaign_2033:151", "first_30pct_2033:151", "pf_0plus_t:156",
               "qdrift_f7p3p1:157", "pf_rabi_t:157",
               "mo_ru_after_lever:169", "mo_ru_after_lever_generic:169", "n_rot_flagship:152",
               "multiple_flagship:152,215", "eps_rot_sd:90", "eps_rot_survey:90", "eps_rot_flagship:152",
               "c_T_flagship:152", "insertion_share_near_miss:217", "insertion_jj44:143", "insertion_jj44_dense:143"}
    assert retired <= set(p), retired - set(p)
    if tex is not None:
        for k, v in r25.items():
            if k in retired:
                continue
            assert v in tex, (k, v)
    for k, v in expected.items():
        assert p[k] == v, (k, p[k], v)
    if tex is None:
        pytest.skip("chapter .tex not on disk")
    for k, v in expected.items():
        if k in retired:
            continue
        # plain multiples are set as '$\times 1.7$--$6.8$' or '14--54' in the prose; accept either spelling
        assert v in tex or v.replace("--", "$--$") in tex, f"{k}: {v!r} not in app11_mu2e_0nubb.tex"
    assert "10^{-2}/N_{\\rm rot}" in tex and "c_T\\approx 27" not in tex     # R-TOL rule sentence, app11:90
    assert "0.1/N_{\\rm gate}^{\\rm shot}" in tex                    # R3 definition, app11:93


def test_utility_box(a):
    # G5 open item 4 (Claude's decision, 2026-10-04): Mu2e-only base, six operator categories, no 0nubb dollars
    u = ch3.utility(a)
    assert u["exposure_equivalent"] == 81                   # "3^4 ~ 80" (0nubb exposure range; no dollar value)
    assert a.mu2e_tpc_musd.prov is c.Provenance.CITED and a.mu2e_tpc_musd.src == "DOE_HEP_FY25_CJ"
    assert u["base_musd"] == 315.7 and u["fraction"] == 0.6 and u["elements"] == 6
    _is(u["musd_attributed_total"], 0.6 * 315.7)            # 189.42 -> "~$190M"
    _is(u["musd_per_instance"], 0.6 * 315.7 / 6)            # 31.57 -> "~$32M"
    assert round(u["musd_attributed_total"], -1) == 190 and round(u["musd_per_instance"]) == 32
    assert not hasattr(a, "nexo_tpc_musd") and not hasattr(a, "utility_leverage_fraction")


def test_utility_box_text_matches_model(a, tex):
    if tex is None:
        pytest.skip("chapter .tex not on disk")
    box = tex.split("begin{utilitybox}")[1].split("end{utilitybox}")[0]
    assert "$\\sim\\$32$M per application instance" in box
    assert ("an attribution of $\\approx\\$190$M, about $60\\%$ of the \\$315.7M Mu2e total project cost"
            "~\\cite{DOE_HEP_FY25_CJ}") in box and "The fraction is not derived." in box
    # style pass 2026-10-08 (rule 9): the "no dollar figure" sentence is cut; the box still carries none
    assert "We attach no dollar figure" not in box and "M per instance" not in box.split("0\\nu\\beta\\beta$ programs")[-1]
    for gone in ("nEXO", "\\$820", "\\$406", "\\$25$M", "0.5$B"):
        assert gone not in box, gone


def test_instance_rows(a, r28, r33, rcd):
    rows28 = ch3.INSTANCE_ROWS(a, "2028", r28)
    rows33 = ch3.INSTANCE_ROWS(a, "2033", r33)
    rowscd = ch3.INSTANCE_ROWS(a, "codesign", rcd)
    assert len(rows28) == 1 and rows28[0][1] == r28.lq and rows28[0][2] == r28.hard_ops
    assert len(rows33) == 2 and rows33[0][3]["fits"] is True      # sd + pf (R16: jj44 moved to co-design)
    assert all(not row[3]["fits"] for row in rows33[1:])
    assert all(row[1][1] <= 250 for row in rows33)          # "<= 250 LQ" for the near-miss cells   app11:238
    assert len(rowscd) == 4 and rowscd[-1][1] == rcd.lq and rowscd[-1][2] == rcd.hard_ops
    assert "jj44" in rowscd[0][0] and rowscd[0][3]["note"].startswith("co-design")     # R16
    assert rowscd[2][3]["note"].startswith("reach target; co-design at the generic gap band")
    assert all(row[1][1] <= 250 for row in rowscd[:3])
    assert ch3.INSTANCE_ROWS(a, "2040", r33) == []


# ------------------------------------------------- R14 fix round (verifier findings)

def test_al27_dipole_shot_factor(a, r33, tex):
    d = ch3.al27_dipole_sd(a)
    # b from r_ch(27Al) = 3.0610 fm, same c.m.-corrected HO fit as 4He, 0s^2 0p^6 (sd)^5 protons
    S = 2 * 1.5 + 6 * 2.5 + 5 * 3.5
    _is(d["b"] ** 2 * (S / 13 - 1.5 / 27), 3.0610 ** 2 - 0.8409 ** 2)
    assert round(d["b"], 2) == 1.80
    y = d["y"]
    _is(d["o_0d"], (1 - 4 * y / 3 + 4 * y * y / 15) * math.exp(-y), 1e-8)
    _is(d["o_1s"], (1 - 4 * y / 3 + 2 * y * y / 3) * math.exp(-y), 1e-8)
    assert (round(d["o_0d"], 3), round(d["o_1s"], 3)) == (0.559, 0.576)
    assert d["gamma"] == 13 and d["gamma_lcu"] == 12           # identity + 12 single-Z strings
    _is(d["lam"], 10 * d["o_0d"] + 2 * d["o_1s"])
    # R15 ruling 1: the identity (half the trace) is added back classically, lambda' = Tr/2 = 3.37
    _is(d["lam_lcu"], d["lam"] / 2)
    assert round(d["lam"], 2) == 6.74 and round(d["lam_lcu"], 2) == 3.37
    lo, hi = d["exp_band_any_state"]
    assert (round(lo, 2), round(hi, 2)) == (2.79, 2.83)
    assert tuple(round(x, 2) for x in d["shot_factor_band"]) == (1.42, 1.45)
    assert tuple(round(x, 1) for x in d["shot_factor_band_with_identity"]) == (5.7, 5.8)     # R14 record
    assert r33.intermediates["al27_dipole_shot_factor"] == d["shot_factor_band"][1]
    if tex is not None:
        # E30 cap pass: b = 1.80 fm and the 0.559 / 0.576 diagonal elements are pinned on the model above
        # style pass 2026-10-08: the lambda = 3.37, 2.79-2.83 and 1.42-1.45 steps are cut (rule 2); pinned above
        for frag in ("fixes the dipole operator to within $\\pm 0.6\\%$",
                     "\\cite{Angeli_Marinova_2013}",
                     "before the factor $(\\lambda/|\\langle\\hat{\\mathcal{O}}\\rangle|)^2$ of Eq.~\\eqref{eq:Nshot11}"):
            assert frag in tex, frag



def test_complete_load_counts(a, r28, tex):
    # Fomichev arXiv:2310.18410 initial.tex:350: (2L-2)D + 2^(L+1) + D Toffolis, L = ceil(log2 D); :355 5L - 3
    assert ch3.sos_load_toffoli(100) == 12 * 100 + 256 + 100 == 1556
    assert ch3.sos_load_toffoli(600) == 18 * 600 + 2048 + 600 == 13448
    assert ch3.sos_load_toffoli(1000) == 18 * 1000 + 2048 + 1000
    assert all(ch3.sos_load_toffoli(d) < (2 * math.ceil(math.log2(d)) + 3) * d for d in (100, 600, 1000, 2000))
    assert ch3.sos_load_ancillae(600) == 47 and ch3.sos_load_ancillae(100) == 32
    if tex is not None:
        # E30 cap pass: the per-term itemization of the load is condensed; 1.1e4 at D = 1e2 is in the box row
        # style pass 2026-10-08: the Toffoli formula is cut (rule 1); the counts stay pinned above
        for frag in ("The state-preparation circuit is the leading Toffoli count of Ref.~\\cite{arxiv_2310_18410}, once per shot.",
                     "needs about the low end of the range counted (box), so the count is conservative.",   # referee M3
                     "\\cite{arxiv_2310_18410}"):
            assert frag in tex, frag
        assert "The box keeps the lookup-only price." not in tex



def test_insertion_share_and_moru_wording(a, tex):
    p = ch3.printed(a)
    assert p["insertion_share_jj44:143"] == "0.04$--$0.2\\%"
    assert p["insertion_share_near_miss:217"] == "0.06$--$0.2\\%"                # R16: pf alone
    if tex is not None:
        assert "$0.04$--$0.2\\%$ of the $jj44$ pair preparation" in tex
        # style pass 2026-10-08: retired, the sentence was cut (rule 1/2); the value stays pinned on the model
        assert "0.5$--$2\\%" not in tex
        assert "$\\Gamma=40016$" in tex
        assert "$12900$ Pauli strings besides the identity" in tex
        assert "5.4$--$5.6" not in tex and "4.4$--$4.6" not in tex and "scaling by $N_{\\rm orb}^{4}$" not in tex
        assert "same factor of two in Eq." not in tex and "each register spans" not in tex
        assert "register asymmetry" not in tex
        # R15 ruling 3: Mo/Ru classification and register
        assert "the $22$ proton states of $jj44$ and the $32$ neutron states of $jj55$, $54$ modes in all" in tex
        assert "falls between the $jj44$ and $jj55$ pairs, both co-design" in tex     # R16
        assert "near-miss cross-check" not in tex and "near-miss pairs" not in tex
        assert "$0\\nu\\beta\\beta$ co-design & $^{76}$Ge/Se, $^{82}$Se/Kr: $\\times 2.6$--$11$" in tex
        assert "$0\\nu\\beta\\beta$ near miss & $^{48}$Ca/Ti: $\\times 1.8$--$7.3$" in tex
        assert "$0\\nu\\beta\\beta$ reach target & $^{100}$Mo/Ru: $\\times 6$--$25$ on state prep" in tex
        assert "co-design target at every plausible spectral gap" not in tex and "$\\times 880$" not in tex



def test_jw_quartic_is_exact_on_a_small_instance():
    # R15 ruling 4: a_p^dag a_q^dag a_r a_s as 16 Pauli strings, checked as an operator identity on 5 modes
    import numpy as np
    n = 5
    I2, X, Y, Z = np.eye(2), np.array([[0, 1], [1, 0]]), np.array([[0, -1j], [1j, 0]]), np.diag([1.0, -1.0])
    lower = np.array([[0, 1], [0, 0]])                         # |1> -> |0>, occupation = |1>

    def kron(ops):
        out = np.array([[1.0]])
        for o in ops:
            out = np.kron(out, o)
        return out

    def a_op(p):                                               # qubit 0 is the leftmost factor
        return kron([Z] * p + [lower] + [I2] * (n - p - 1))

    def pauli(x, z):
        ops = []
        for q in range(n):
            xb, zb = (x >> q) & 1, (z >> q) & 1
            ops.append(X @ Z if xb and zb else X if xb else Z if zb else I2)   # X^x Z^z per qubit
        return kron(ops)
    for p, q, r, s in ((0, 2, 3, 4), (1, 3, 0, 4), (4, 0, 2, 1)):
        target = a_op(p).conj().T @ a_op(q).conj().T @ a_op(r) @ a_op(s)
        terms = ch3.jw_quartic(p, q, r, s)
        built = sum(c * pauli(x, z) for c, x, z in terms)
        assert np.allclose(built, target), (p, q, r, s)
        assert len(terms) == 16 and all(abs(abs(c) - 1 / 16) < 1e-12 for c, _, _ in terms)


def test_insertion_counting_rules():
    # pair quantum numbers: two nucleons in one j-shell couple to even J only
    st = ch3.shell_states("jj44")
    g = [s for s in st if s[0] == 3]                          # 0g9/2, 10 states
    M, par, Js = ch3.pair_quantum_numbers(ch3.SHELLS["jj44"], g[0], g[-1])   # m = -9/2, +9/2: M = 0
    assert M == 0 and par == 0 and Js == frozenset({0, 4, 8, 12, 16})        # 2J = 0, 4, ..., 16
    # the selection-rule count is explicit-JW checked (Gamma = 16 x terms) and well under the dense count
    for sp, sn in (("pf", "pf"), ("jj44", "jj44"), ("jj55", "jj55"), ("jj44", "jj55")):
        d = ch3.gamma_0nubb(sp, sn)
        assert d["checked"] and d["gamma"] == 16 * d["n_terms"] < 16 * d["dense_terms"]
    assert ch3.gamma_0nubb("jj44", "jj44")["n_terms"] == 2501
    # the sd Mu2e operator: identity present and dropped from Gamma
    sd = ch3.gamma_mu2e_two_body("sd")
    assert sd["identity"] and sd["gamma"] == 12900 and sd["modes"] == 24


def test_r17_registers_printed_and_load_count(r33, tex):
    # r17 (accepted small items): the 2033 box and requirements table print the registers, not '<= 250';
    # the complete load is Fomichev's leading Toffoli count, the angle loading not itemized
    s = ch3.survey(ch3.Assumptions())
    assert r33.lq == (74, 174) and s["pf"]["lq"] == (90, 190) and s["jj44"]["lq"] == (94, 194)
    assert (s["pf"]["qubits"], s["jj44"]["qubits"]) == (40, 44)
    if tex is None:
        pytest.skip("chapter .tex not on disk")
    assert "\\le 250" not in tex
    assert "Logical qubits & $74$--$174$ \\\\" in tex
    assert "$N_{\\rm anc}\\sim 50$--$150$ ancillae, $74$--$174$ LQ" in tex
    assert "$90$--$190$ and $94$--$194$ LQ" in tex
    # referee M3 (2026-10-04): the angle rotations are now estimated in the text (load_angle_cost)
    assert ("That source does not itemize the rotations that apply the amplitude angles;") in tex   # style pass: L dropped
    assert "and we do not price them" not in tex


# ------------------------------------------------- referee report 2026-10-04 (M1-M4, G2, G5, editorial)

def test_referee_g2_m3_load_angle_rotations(a, r28, tex):
    # G2: the 2028 cost object is the executed shot (insertion + one complete load), 1.6e4-9.9e4 T.
    # M3: the load's amplitude-angle rotations, derived from LKS 1812.00954 (phase-gradient CAdd, b - 1 Toffolis
    # per level, Gidney 1709.06648) with b = 13 (2 pi L / 2^b <= 1e-2), plus b - 3 synthesized Fourier-state rotations
    i = r28.intermediates
    lo, hi = i["load_angle"][100.0], i["load_angle"][600]
    assert (lo["L"], lo["b"], lo["toffoli"], lo["n_rot"]) == (7, 13, 84, 10)
    assert (hi["L"], hi["b"], hi["toffoli"], hi["n_rot"]) == (10, 13, 120, 10)
    assert 2 * math.pi * 10 / 2 ** 13 <= 1e-2 < 2 * math.pi * 10 / 2 ** 12
    assert tuple(ch3.round_sig(x) for x in i["load_angle_t"]) == (770, 1000)       # "7.7e2-1.0e3"
    assert round(i["fraction_with_load_angles"][1], 3) == 1.005                      # "1.005x the cap"
    assert i["d_break_with_load_angles"] == 596                                       # "the fit then ends at D = 596"
    # verifier (2026-10-04): the phase-gradient state is catalytic, charged once per circuit, not per application
    x = i["mlae_with_load_angles"]
    assert x["m_max"] == 5 and x["t_deepest"] <= 1e5 < x["t_next_depth"]
    _is(x["t_deepest"], i["mlae_no_angles"][0]["t_deepest"] + 5 * 7 * lo["toffoli"]
        + lo["n_rot"] * c.t_per_rotation(c.eps_rot_for(10 + 5 * 254, 1e-2)) + 1, 1e-4)
    assert ch3.round_sig(i["mlae_t_deepest_with_load_angles"]) == 8.7e4               # "8.7e4 with the angle rotations"
    # m_max with the angle term: 5 to D = 128 (unchanged, L steps at 129), 3 to D = 221, 1 to the D = 596 fit
    assert i["mlae_m_max_last_d_with_load_angles"] == {5: 128, 3: 221, 1: 596}
    assert ch3.round_sig(i["mlae_no_angles"][0]["t_deepest"]) == 8.4e4
    # M3 ruling (2026-10-04, "do 3"): the angle term is folded into the headline, the walls, depth and register
    assert r28.hard_ops == i["hard_ops_with_load_angles"] and i["hard_ops_no_angles"][1] < 1e5 < r28.hard_ops[1]
    assert i["mlae"][0] is x
    # maximum subroutine time: 0.10 s per 2028 shot (D = 600)
    assert round(i["max_subroutine_time_s"][1], 2) == 0.10
    if tex is None:
        pytest.skip("chapter .tex not on disk")
    # style pass 2026-10-08: retired, the sentence was cut (rule 1/2); the value stays pinned on the model
    for frag in ("we count them as phase-gradient additions~\\cite{arxiv_1812_00954,arxiv_1709_06648}",
                 "(fits to $D=596$; schematic $^4$He needs $77$--$103$ at fidelity $0.999$)",
                 "$+$ $7.7\\times 10^{2}$--$1.0\\times 10^{3}$ (amplitude rotations)",
                 "by Hadamard test on a prepared ground state]"):
        assert frag in tex, frag
    assert "Size (2028 benchmark)" not in tex          # cut pass 2026-10-06: row restated the box


def test_referee_max_subroutine_time(a, tex):
    # editorial: the requirements row reports the deepest-circuit durations instead of N/A
    i33 = ch3.model(a, "2033").intermediates
    rabi = i33["simplified_0nubb"]["f7p3p1_rabi"]
    per_shot = [w / n for w, n in zip(rabi["wall_s"], rabi["shots"])]
    assert round(max(per_shot) / 60) == 16                                          # "~16 min per shot"
    assert (round(i33["t_shot_s"][0] / 60, 1), round(i33["t_shot_s"][1] / 60, 1)) == (1.8, 7.3)
    if tex is None:
        pytest.skip("chapter .tex not on disk")
    assert "Maximum subroutine time & N/A" not in tex
    assert ("Maximum subroutine time & $\\sim 16$ min per shot ($0\\nu\\beta\\beta$ validation stage, resonant "
            "transition); $1.8$--$7.3$ min ($^{27}$Al, 2033); $0.10$\\,s (2028)") in tex


def test_referee_text_corrections(a, tex):
    if tex is None:
        pytest.skip("chapter .tex not on disk")
    # M4
    assert "zero-background" not in tex
    assert "An established conversion signal would demonstrate charged LFV." in tex
    # M1: the half-life formula is the light-exchange mechanism; other mechanisms enter separately
    assert "For $0\\nu\\beta\\beta$ mediated by light-Majorana-neutrino exchange alone," in tex
    assert "at comparable weight" not in tex
    assert "enters separately, with its own coupling, NME and phase-space factor" in tex   # style pass: eta_X one use
    # M2: the priced 2033 readout is the bare SD element; one rate fixes one combination
    assert "for all six operator categories (dipole primary)" not in tex
    assert "One conversion rate on one target fixes one combination of the coefficients." in tex
    assert "spin-independent categories validated classically" in tex
    # G5: no 0nubb spread transferred to Mu2e; the 3^4 range is not a realized gain
    assert "degraded by a factor $2$--$3$ by NME uncertainty alone" not in tex
    assert "is equivalent to a factor $3^{4}" not in tex and "more leverage than any feasible" not in tex
    assert "A precise NME does not supply that exposure." in tex
    # editorial: no glossary clash (qudit error, qudit dimension); terminology
    assert "\\epsilon_q" not in tex and "degree-$d$" not in tex
    assert "first light" not in tex and "rung" not in tex
    # verifier (2026-10-04): M2 box rows, M4 sourcing, SD-element motivation
    assert "on $\\mathrm{CR}(\\mu^-\\to e^-,{}^{27}\\mathrm{Al})$; $5\\%$" not in tex
    # style pass 2026-10-08 (rules 9, 10): box rows carry no scope disclaimers or "not derived" tags; the
    # five-category assumption is stated once in the 2028 prose
    assert "Target observable & spin-dependent element to $5\\%$" in tex
    assert "all six categories (other five at the dipole's shot factor and $\\Gamma\\le 128$)" in tex
    assert "whose string counts are not derived here" in tex
    assert "An expression of interest for Mu2e-II was submitted in 2018" in tex
    assert "nEXO was not selected" not in tex and "We do not date them here." not in tex
    assert "pre-loaded" not in tex and "---" not in tex
    assert "so it isolates the operators that couple to nuclear spin~\\cite{arxiv_2208_07945}" in tex


def test_referee_m3_he4_data_size(a, r28, tex):
    """Referee M3 (2026-10-05): D for the 4He ground state, computed for a schematic stand-in (Minnesota NN,
    Lawson beta = 2 and 10, exact diagonalization in the 40-mode space; scratch chopen/app11/he4D.py). The
    computed D is at or below the priced band D = 1e2-600 and below the fit's end (596); the headline is unchanged."""
    i = r28.intermediates["d_he4_computed"]
    assert i["fid99"] == (15, 17) and i["fid999"] == (77, 103) and i["dipole_1pct"] == 1
    assert i["at_or_below_priced_band"] and i["below_priced_band_low_end"] and i["fits"]
    assert a.he4_d_fid999.hi <= r28.intermediates["d_break"] == 596
    assert tuple(ch3.round_sig(x) for x in r28.hard_ops) == (1.7e4, 1.0e5)          # headline unchanged
    if tex is None:
        pytest.skip("chapter .tex not on disk")
    # style pass 2026-10-08: the 0.99 / D = 15-17 and 1% steps are cut (rule 2); 77-103 is in the box
    for frag in ("the Minnesota potential~\\cite{Minnesota_NN_1977} with a Lawson term",
                 "schematic $^4$He needs $77$--$103$ at fidelity $0.999$",
                 "A chiral ground state may need a larger $D$ and is not checked here."):
        assert frag in tex, frag
    assert "plausibly needs" not in tex and "$D$ for $^4$He not computed" not in tex


def test_referee_m3_al27_usdb_gap(a, r33, tex):
    """Referee M3 (2026-10-05): 27Al keeps the generic 1-2 MeV band; the USDB gap (2.33 MeV, the 7/2+, nearest state
    an M = 5/2 reference reaches; exact Lanczos in sd) is printed as the check that the band is conservative."""
    assert a.gap_al27_usdb_mev.value > a.gap_mev.hi
    lo, hi = r33.intermediates["al27_prep_t_at_usdb_gap"]
    assert (ch3.round_sig(lo), ch3.round_sig(hi)) == (8.9e7, 1.8e8)
    assert tuple(ch3.round_sig(x) for x in r33.hard_ops) == (1.1e8, 4.4e8)          # 2033 headline unchanged
    if tex is None:
        pytest.skip("chapter .tex not on disk")
    assert "in full $pf$ GXPF1A costs about $2\\%$ more" in tex
    assert ("$2.33$ MeV above the ground state (the $7/2^+$), which would lower the preparation to "
            "$0.89$--$1.8\\times 10^{8}$ T") in tex and "\\cite{Brown_Richter_USDB_2006}" in tex

