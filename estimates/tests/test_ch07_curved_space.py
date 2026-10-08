"""Ch. 7 (curved spacetime): every named number in app04_curved_space.tex, as the
chapter states it, so that editing one without the others fails loudly.

The generic contract tests (test_common.py) already check provenance, the
breakdown sum, and closeness to PUBLISHED. These pin the prose intermediates.
"""

import math
from pathlib import Path

import pytest

from estimates import common as c
from estimates import ch07_curved_space as ch7

TEX = Path(__file__).resolve().parents[3] / "applications" / "app04_curved_space.tex"


@pytest.fixture(scope="module")
def a():
    return ch7.Assumptions()


@pytest.fixture(scope="module")
def r28(a):
    return ch7.model(a, "2028")


@pytest.fixture(scope="module")
def r33(a):
    return ch7.model(a, "2033")


def _is(x, y, rel=1e-9):
    assert math.isclose(x, y, rel_tol=rel), (x, y)


# ---------------------------------------------------------------- provenance

def test_every_input_is_tagged(a):
    import dataclasses
    for f in dataclasses.fields(a):
        assert isinstance(getattr(a, f.name), c.Tagged), f.name


def test_per_step_costs_are_derived_not_stated(a):
    # R12: app04:84 derives Gamma_step term by term; no Stated per-step constant remains
    assert not hasattr(a, "t_per_step_2028") and not hasattr(a, "t_per_step_2033")
    assert not hasattr(a, "t_per_prep_step_2033")
    # the 2028 Gaussian-network prep is still Stated, and the chapter says so (app04:78)
    assert a.t_prep_2028.prov is c.Provenance.STATED


def test_yukawa_form_is_author_confirmed(a):
    # R13: the author confirmed g phi psibar psi (2026-09-30); the chapter writes it at app04:84
    assert a.yukawa_on.prov is c.Provenance.STATED
    assert a.yukawa_on.src == "app04:84" and a.yukawa_on.lo == 1
    assert "author-confirmed 2026-09-30" in a.yukawa_on.note


def test_gamma_step_is_exercisable(a):
    # changing V, K or N_f now moves the per-step T, not only the register
    base = ch7.model(a, "2033").intermediates["t_per_step"]
    for kw in (dict(L_2033=c.Stated(6, "test")), dict(K_2033=c.Stated(8, "test")),
               dict(N_f_2033=c.Stated(2, "test"))):
        b = ch7.Assumptions(**kw)
        assert ch7.model(b, "2033").intermediates["t_per_step"] != base, kw
    with pytest.raises(ValueError):
        ch7.Assumptions(K_2033=c.Stated(12, "test"))    # not a power of two


def test_term_counts_uniform_grid():
    # phi linear in Z (arXiv:2210.07985 Eq. phiqub): C(n,2)+C(n,4) on-site, n^2 per link
    a = ch7.Assumptions()
    s28, s33 = ch7.step_terms(a, "2028"), ch7.step_terms(a, "2033")
    assert (s28["scalar_site_strings"], s28["gradient_link_strings"], s28["pi2_site_strings"]) == (3, 9, 3)   # K = 8
    assert (s33["scalar_site_strings"], s33["gradient_link_strings"], s33["pi2_site_strings"]) == (7, 16, 6)
    assert s33["fermion_site_strings"] == 4 + 16          # Wilson mass + Yukawa 4 n_q
    assert s33["hop_link_strings"] == 16                  # 8 bilinears x 2 strings
    assert ch7.wilson_hop_nonzeros() == [8, 8, 8]
    assert (s28["cs_per_qft"], s28["cr_per_qft"]) == (2, 1)   # K = 8 (option (c))
    assert (s33["cs_per_qft"], s33["cr_per_qft"]) == (3, 3)
    assert ch7.centering_layer(2) == (0, 1)               # 3pi/2 Clifford, 3pi/4 one T
    assert ch7.centering_layer(4) == (2, 1)
    assert ch7.centering_layer(3) == (1, 1)
    assert s28["diag_block"] == 1920 and s28["pi2_block"] == 192 and s28["hop_block"] == 0
    assert s28["qft_pair_rot"] == 384 and s28["qft_pair_cs"] == 256
    assert s33["diag_block"] == 9375 and s33["pi2_block"] == 750 and s33["hop_block"] == 6000
    assert s33["qft_pair_rot"] == 2250 and s33["qft_pair_cs"] == 750


def test_no_commutator_lever(a):
    # R11 ch07-assumed-lever (a): app04:81 "We therefore price the worst case"
    assert a.commutator_lever.prov is c.Provenance.STATED
    assert a.commutator_lever.lo == 1
    _is(a.stability_omega_dt.lo, 2)


def test_x10_lever_is_rejected_by_accuracy():
    # U2 (2026-10-04): any lever that lengthens the 2033 step past m_phi dt = 0.05 fails the accuracy criterion
    with pytest.raises(ValueError, match="accuracy"):
        ch7.Assumptions(commutator_lever=c.Assumed(10, "test"))
    with pytest.raises(ValueError, match="accuracy"):
        ch7.Assumptions(dt_Hinf_2033=c.Stated(0.01, "test"))     # the retired m_phi dt = 1 step


def test_fault_budget_is_the_r3_convention(a):
    # R3 (2026-09-28): 0.1 expected logical faults per shot, report-wide; Ch. 7's "~90% survival"
    # is the same criterion (exp(-0.1) = 0.905). A named, Assumed input, not a hidden 0.9.
    assert a.fault_budget_per_shot.prov is c.Provenance.ASSUMED
    _is(a.fault_budget_per_shot.lo, 0.1)
    assert 0.90 <= math.exp(-a.fault_budget_per_shot.lo) <= 0.91


def test_breakdown_is_primitive_level(r28, r33):
    # R12: every Trotter primitive is an explicit gate-level count (COMPILED, 'derived here');
    # only the 2028 Gaussian-network prep stays SCALING (Stated, app04:78)
    for r in (r28, r33):
        for p in r.breakdown:
            if p.name == "gaussian_network_prep":
                assert p.status is c.CircuitStatus.SCALING
            else:
                assert p.status is c.CircuitStatus.COMPILED, p
                assert "derived here" in p.src or "2210.07985" in p.src, p
    assert [p.name for p in r28.breakdown].count("gaussian_network_prep") == 1


def test_no_codesign_box(a):
    with pytest.raises(ValueError):
        ch7.model(a, "codesign")


def test_dt_stability_rule_enforced():
    # Delta t <~ 1/omega_max (app04:73) is validated at construction
    with pytest.raises(ValueError):
        ch7.Assumptions(dt_Hinf_2028=c.Stated(1.0, "test"), omega_max_Hinf_2028=c.Stated(2, "test"))


# ------------------------------------------------------------- 2028 register

def test_2028_register(r28):
    i = r28.intermediates
    assert i["V"] == 64                      # "V = 4^3 (64 sites)"        app04:112
    # smaller ruling (b), 2026-10-05: option (c), one step at K = 8 (was K = 4, 128 system, 158-178)
    assert i["qubits_per_site"] == 3         # "K = 8, 3 qubits/site"      box
    assert i["system_qubits"] == 192         # "192 system qubits (64 x 3)"
    assert i["n_anc"] == (30, 50)            # "30-50 ramp/out-mode/..."   app04:120
    assert i["n_hadamard"] == (0, 0)         # R6 (r25): no Hadamard test
    assert i["lq_sum"] == (222, 242)         # "222--242 in all"; box prints 220-240 (2 sig figs)
    assert r28.lq == (222, 242)
    assert i["log2_hilbert_dim"] == 192      # "dimension ~ K^V" = 8^64 = 2^192


# ------------------------------------------------------------ 2028 hard ops

def test_2028_time_step(r28):
    i = r28.intermediates
    _is(i["dt_Hinf"], 0.5)                   # "Delta t = 0.5/H_inf"       app04:73,123
    _is(i["omega_max_Hinf"], 2)              # "omega_max ~ 2 H_inf"       app04:123
    _is(i["dt_times_omega_max"], 1.0)        # Delta t ~ 1/omega_max       app04:65
    _is(i["expansion_time_over_step"], 2)    # "which is >= 2 ... times the step"  app04:73
    assert i["n_trotter_steps"] == 1         # option (c): "as many ... as the cap holds ... At K=8 that is one"
    _is(i["T_phys_Hinf"], 0.5)               # "The step spans T_phys = 0.5/H_inf"
    _is(i["time_ordering_over_trotter_error"], 0.5)   # "smaller ... by H_inf/omega_max"  app04:73


def test_2028_hard_ops_derived(r28):
    # smaller ruling (b), 2026-10-05, option (c): one macro-step at K = 8; B/2 D B/2 (app04:84)
    i = r28.intermediates
    assert i["diag_block_rotations"] == 64 * 3 + 192 * 9 == 1920
    assert i["pi2_block_rotations"] == 64 * 3 and i["qft_pair_rotations"] == 64 * 2 * 1 * 3
    assert i["n_rot_per_shot"] == 1920 + 2 * 192 + 2 * 384 + 128 == 3200   # "3200 rotations" (128 centering)
    _is(i["t_exact_per_shot"], 1664)                      # 1536 cS + 128 centering
    _is(i["t_controlled_s"], 1536)
    _is(i["eps_rot"], math.sqrt(1e-2 / 3200))             # 1.77e-3 -> "1.8e-3"
    assert round(i["eps_rot"], 4) == 1.8e-3
    _is(i["t_per_rotation"], 1.15 * math.log2(1 / math.sqrt(1e-2 / 3200)) + 9.2)
    assert 19.65 <= i["t_per_rotation"] < 19.75           # "19.7 T per rotation"
    _is(i["t_evolution"], 3200 * i["t_per_rotation"] + 1664)
    assert 6.45e4 <= i["t_evolution"] < 6.55e4            # 64,753 -> "6.5e4"
    _is(i["t_per_step"], i["t_evolution"])
    _is(i["t_single_step_alone"], i["t_evolution"])
    _is(i["t_prep"], 3e4)                                 # Stated 2e4 scaled to 192 system qubits
    _is(i["t_prep_per_system_qubit"], 156.25)             # "~150 T per system qubit"
    _is(i["t_per_shot"][0], i["t_evolution"] + 3e4)
    assert 9.45e4 <= i["t_per_shot"][0] < 9.5e4           # 94,753 -> "9.5e4"   box
    assert 0.945 <= i["ratio_to_cap"] < 0.955             # "0.95 of the cap"
    assert i["macro_steps_that_fit"] == 1                 # "At K=8 that is one"
    _is(i["shot_t_at_n_steps"][1], i["t_per_shot"][0])
    assert 1.45e5 <= i["shot_t_at_n_steps"][2] < 1.5e5    # "a second would bring the shot to 1.5e5 T"
    with pytest.raises(ValueError, match="exceed"):
        ch7.model(ch7.Assumptions(n_steps_2028=c.Assumed(2, "test")), "2028")
    _is(r28.breakdown_total(), i["t_per_shot"][0])
    assert r28.hard_ops == (i["t_per_shot"][0], i["t_per_shot"][0])
    # the retired four-step K = 4 benchmark reproduces (record)
    old = ch7.model(ch7.Assumptions(K_2028=c.Stated(4, "test"), n_steps_2028=c.Stated(4, "test"),
                                    t_prep_2028=c.Stated(2e4, "test")), "2028")
    _is(old.hard_ops[0], ch7.RETIRED_2028_K4["t_per_shot_2028"], rel=1e-6)
    assert old.lq == ch7.RETIRED_2028_K4["lq_2028"]


def test_2028_shots_and_wall_time(r28):
    # R6 (r25): three computational-basis settings, every mode read from the same shots
    i = r28.intermediates
    assert i["n_settings"] == 3                           # "three settings"
    _is(i["samples_for_precision"], 400)                  # 5% -> 400 samples
    assert i["worst_shell_modes"] == 3                    # lowest shell: 6 momenta, 3 complex modes
    assert i["shots_min_per_setting"] == 134              # 400 / 3
    assert i["shots_per_setting"] == 150                  # "150 shots per setting"
    assert i["shots"] == 450                              # "450 shots"
    assert 0.045 <= i["stat_precision_worst_shell"] <= 0.05
    _is(i["readout_extra_t"], 0.0)                        # settings reuse the last pi^2 half-step
    # report rule: 1 us per T plus 0.1 ms per shot, serial on one machine (E27, r23)
    _is(i["t_gate_s"], 1e-6)
    _is(i["shot_overhead_s"], 1e-4)
    _is(i["t_shot_s"], i["t_per_shot"][0] * 1e-6 + 1e-4)
    _is(i["wall_time_s"], 450 * (i["t_per_shot"][0] * 1e-6 + 1e-4))   # 42.68 s
    assert 42.5 <= i["wall_time_s"] < 43.5                # "43 s"
    assert i["wall_first_result_s"] is None               # one tier
    _is(i["wall_campaign_s"], i["wall_time_s"])
    with pytest.raises(ValueError, match="worst-shell"):
        ch7.model(ch7.Assumptions(shots_per_setting_2028=c.Assumed(100, "test")), "2028")


def test_2028_depth_exports(r28):
    # R9 (r25): T-depth from factory.json (346 low, 1308 high with parallel prep)
    i = r28.intermediates
    assert i["t_depth_per_shot"] == (346.0, 1308.0)
    assert 72 <= i["f_star"][0] < 73 and 273 <= i["f_star"][1] < 275   # K = 4 depth kept for option (c)
    assert all(i["baseline_ok"])
    _is(i["floor_wall_s"][1], 450 * 1308 * 1e-5)


def test_2028_envelope_claim_and_rfi_survival(r28):
    i = r28.intermediates
    # app04:90 "O(10)-step preheating run within the same hard-op envelope": what that needs
    assert i["envelope_steps"] == 10
    _is(i["t_per_step_for_envelope"], 7e3)           # (1e5 - 3e4) / 10 at K = 8
    assert 9.0 <= i["step_reduction_for_envelope"] < 9.5     # "from 6.5e4 to <~ 7e3 T at K=8, about 9x"
    assert i["shot_t_at_envelope_steps"] > 1e5
    # the 2028 box states no eps_l requirement, only that the RFI 1e-8 suffices (app04:122)
    assert r28.epsilon_l is None
    t = i["t_per_shot"][0]
    _is(i["expected_faults_at_rfi"], t * 1e-8)
    assert 9.45e-4 <= i["expected_faults_at_rfi"] < 9.55e-4   # "9.5e-4 expected faults"
    assert 0.99905 <= i["shot_survival_at_rfi"] < 0.99915     # "99.91% shot survival"
    _is(i["fault_budget_per_shot"], 0.1)
    _is(i["eps_l_for_fault_budget"], 0.1 / t)
    _is(i["eps_l_for_survival_target"], -math.log(0.9) / t)


# ------------------------------------------------------------- 2033 register

def test_2033_register(r33):
    i = r33.intermediates
    assert i["V"] == 125                       # "V = 5^3 (125 sites)"     app04:137
    assert i["qubits_per_site_scalar"] == 4    # "K = 16, 4 qubits/site"   app04:138
    assert i["qubits_per_site_fermion"] == 4   # "1 Wilson fermion (4 qubits/site)"  app04:138
    assert i["system_qubits"] == 1000          # "1000 system (2 x 125 x 4)"  app04:144
    assert i["lq_sum"] == (1030, 1050)         # "1030--1050 in all"; box prints 1050
    assert r33.lq == (1030, 1050)
    assert i["log2_hilbert_dim_scalar"] == 500     # K^V = 16^125 = 2^500      app04:53
    assert i["log2_hilbert_dim_fermion"] == 500    # 2^(4 N_f V) = 2^500       app04:53
    assert i["log2_hilbert_dim"] == 1000
    _is(i["H_inf_over_Mpl"], 1e-5)                 # "H_inf = 1e-5 M_Pl"       app04:139,159


# ------------------------------------------------------------ 2033 hard ops

def test_2033_time_step(r33):
    i = r33.intermediates
    _is(i["T_phys_Hinf"], 10)                  # "T_phys H_inf = 10"
    assert 1.5 <= i["e_folds_matter_like"] <= 2.5   # "~2 e-folds of the matter-like background"
    _is(i["dt_Hinf"], 5e-4)                    # "Delta t = 5e-4/H_inf"  (U2)
    _is(i["omega_max_Hinf"], 1e2)              # m_phi ~ 1e2 H_inf
    _is(i["dt_times_omega_max"], 0.05)         # "m_phi Delta t = 0.05"
    _is(i["expansion_time_over_step"], 2e3)    # "2 to 2e3 times the step"
    _is(i["n_trotter_steps_raw"], 2e4)         # "2e4 Trotter steps"
    _is(i["time_ordering_over_trotter_error"], 0.01)   # H_inf/omega_max


def test_2033_per_step_derived(r33):
    i = r33.intermediates
    assert i["per_step_rotations"] == 9375 + 750 + 2 * 6000 + 2250    # 24,375 -> "2.44e4"
    _is(i["per_step_exact_t"], 2250)                                   # "2250 T"
    # U2: 2e4 evolution + 500 ramp steps in one product, + the 500-rotation boost (U1)
    assert i["n_rot_per_shot"] == 20501 * 9375 + 20500 * (750 + 12000 + 2250) + 500 + 500
    _is(i["eps_rot"], math.sqrt(1e-2 / i["n_rot_per_shot"]))
    assert 4.45e-6 <= i["eps_rot"] < 4.55e-6                           # "4.5e-6"
    assert 29.55 <= i["t_per_rotation"] < 29.65                        # "29.6 T"
    _is(i["t_per_step"], 24375 * i["t_per_rotation"] + 2250)
    assert 7.2e5 <= i["t_per_step"] < 7.3e5                            # "7.2e5 T"


def test_2033_evolution_has_no_lever(r33):
    i = r33.intermediates
    _is(i["commutator_lever"], 1)
    _is(i["n_trotter_steps_net"], 2e4)
    _is(i["t_evolution_net"], 2e4 * i["t_per_step"])
    assert 1.4e10 <= i["t_evolution_net"] < 1.46e10   # "1.4e10 T of evolution", "1.45e10"


def test_hwp_lever_not_taken(r28, r33):
    """app04:84 'Hamming-weight phasing ... would roughly halve the rotation cost at 31 workspace qubits per
    32-rotation group. The 2028 shot fits without it, the 2033 register has no room for it' (E26, r22: a group of k
    holds k - w(k) ancilla, the same number as its Toffolis; the chapter printed 37 = 31 + a 6-qubit weight register)."""
    for r in (r28, r33):
        i = r.intermediates
        assert i["hwp32_ancilla"] == c.hwp_ancilla(32) == c.hwp_toffolis(32) == 31
        assert 0.4 < i["hwp32_t_over_plain"] < 0.55        # "roughly halve": 0.53 at 2028, 0.43 at 2033
    assert r33.lq[0] + 31 > 1000                           # "the 2033 register has no room for it"
    if TEX.exists():
        tex = TEX.read_text()
        assert r"at $31$ workspace qubits per $32$-rotation group" in tex
        assert r"$37$ workspace qubits" not in tex


def test_2033_step_accuracy_argument(r33):
    # U2: Stormer-Verlet cos(w~ dt) = 1 - (w dt)^2/2; at w dt = 1 a 4.72% frequency error, 47 rad over 1e3
    i = r33.intermediates
    _is(i["omega_dt_taken"], 0.05)
    assert i["omega_dt_taken"] < i["stability_omega_dt"]
    _is(i["verlet_freq_err_at_omega_dt_1"], math.pi / 3 - 1)          # "4.7%"
    assert 47.0 <= i["verlet_phase_err_at_omega_dt_1"] < 47.5         # "47 rad"
    assert i["verlet_freq_err_taken"] < 1.1e-4
    _is(i["omega_dt_ramp"], 1.0)
    assert 1.4 * i["omega_dt_ramp"] < i["stability_omega_dt"]         # "fastest mode at omega dt ~ 1.4"
    lo, hi = i["switch_occupation"]                                   # sudden switch, map-vacuum aspect ratio
    assert 5e-3 <= lo < 5.5e-3 and 2.9e-2 <= hi < 3.1e-2               # "5e-3 to 3e-2 quanta per mode"
    assert i["switch_occupation_swept"] < 2e-5                        # "below 2e-5" after the 50-step sweep
    assert i["dt_sweep_steps"] == 50
    assert 3.5e7 <= i["t_dt_sweep"] < 3.7e7                           # "3.6e7 T"
    assert 2e-3 <= i["dt_sweep_fraction"] < 2.5e-3                    # "0.2% of the shot", not priced


def test_2033_ramp_window_ruling(a, r33):
    # R17 (d): window 5/H_inf (500/m_phi), 500 steps at 0.01/H_inf; U2 reprice: 3.6e8 T, 2.5% of the evolution
    i = r33.intermediates
    assert a.ramp_window_Hinf_2033.prov is c.Provenance.STATED
    _is(a.ramp_window_Hinf_2033.lo, 5)
    assert i["n_prep_steps"] == 500
    _is(i["dt_ramp_Hinf"], 0.01)
    _is(i["ramp_window_Hinf"], 5)
    _is(i["ramp_window_over_m_phi_period"], 500)
    _is(i["t_prep"], 500 * i["t_per_step"])    # same Hamiltonian, same per-step count
    assert 3.55e8 <= i["t_prep"] < 3.65e8      # "3.6e8 T"
    _is(i["prep_fraction"], 0.025)             # ramp/evolution
    assert 0.0235 <= i["t_prep"] / i["t_per_shot"][0] < 0.0245           # "2.4% of the shot"
    assert 0.0235 <= i["vacuum_reference_fraction"] < 0.0245             # prep-only run, "2.4% of a shot"


def test_2033_de_sitter_variant(r33):
    # app04:103: "a de Sitter run over the same T_phys (... dt ~ 0.5/H_inf as at 2028, ~20 steps, ~1.2e7 T)
    # is the cheaper variant"
    i = r33.intermediates
    _is(i["de_sitter_variant_steps"], 20)
    assert 1.15e7 <= i["de_sitter_variant_t"] < 1.25e7   # "~1.2e7 T"   app04:103


def test_2033_per_shot_total_is_sum_of_box_lines(r33):
    i = r33.intermediates
    # U2: 2e4 x 7.25e5 + 500 x 7.25e5 + boost + closing block = 1.486e10, box "1.5e10"
    _is(i["t_per_shot"][0], i["t_evolution_net"] + i["t_prep"] + i["t_condensate_boost"] + i["t_closing"])
    assert 2.8e5 <= i["t_closing"] < 3.0e5
    assert 1.45e4 <= i["t_condensate_boost"] < 1.55e4      # "1.5e4-T boost"
    assert i["condensate_boost_rotations"] == 125 * 4      # "500 single-qubit rotations"
    assert 1.48e10 <= i["t_per_shot"][0] < 1.49e10
    assert r33.hard_ops == (i["t_per_shot"][0], i["t_per_shot"][0])
    _is(r33.breakdown_total(), i["t_per_shot"][0])
    assert c.close(r33.hard_ops, (1.5e10, 1.5e10), 0.01)
    # "15 times the 1e9 reference budget, far past the 1.5x overshoot we accept"
    assert 14.5 <= i["ratio_to_ceiling"] < 15.0
    assert i["ratio_to_ceiling"] > 1.5


def test_2033_shots_wall_time_campaign(r33):
    # R1/R6 (r25): two settings, N = N_set ((n+1/2)/n)^2 / (m' eps^2) per profile, x1-2 margin
    i = r33.intermediates
    t = i["t_per_shot"][0]
    assert i["n_settings"] == 2 and i["occupation_band"] == 1 and i["worst_shell_modes"] == 3
    _is(i["shots_per_profile_campaign_strict"], 150)       # 2 x 2.25 / (3 x 0.01)
    _is(i["shots_per_profile_first_strict"], 150 / 9)
    assert i["shots_per_profile_first"] == (17, 34)        # "17--34 shots"
    assert i["shots_per_profile_campaign"] == (150, 300)   # "150--300 per profile"
    assert i["profiles_first"] == 1 and i["profiles_campaign"] == (10, 20)
    assert i["shots_first"] == (17, 34)
    assert i["shots_campaign"] == (1500, 6000)             # "1.5--6.0e3"
    _is(i["t_shot_s"], t * 1e-6 + 1e-4)        # "1.5e4 s", max subroutine time
    assert 1.45e4 <= i["t_shot_s"] < 1.55e4
    assert 4.05 <= i["t_shot_hours"] < 4.15     # "4.1 h"
    lo, hi = i["wall_first_result_d"]
    assert 2.85 <= lo < 2.95 and 5.75 <= hi < 5.85          # "2.9--5.8 days"
    lo, hi = i["wall_campaign_yr"]
    assert 0.705 <= lo < 0.715 and 2.75 <= hi < 2.85        # "0.71--2.8 yr"
    _is(i["wall_campaign_s"][1], 6000 * (t * 1e-6 + 1e-4))
    _is(i["wall_first_result_s"][0], 17 * (t * 1e-6 + 1e-4))
    assert i["campaign_fits_horizon"] is True              # "inside the 5-year horizon"
    # Result carries the campaign band as flat (lo, hi) floats (common.Result contract)
    assert r33.shots == (1500, 6000)
    _is(r33.wall_time_s[0], 1500 * (t * 1e-6 + 1e-4))
    _is(r33.wall_time_s[1], 6000 * (t * 1e-6 + 1e-4))
    row = ch7.INSTANCE_ROWS(ch7.Assumptions(), "2033", r33)[0][3]
    assert row["shots"] == 1500 and row["shots_hi"] == 6000


def test_2033_depth_exports(r33):
    # R9 (r25): factory.json 27 / 92 rotation layers per step, now over 20,500 steps: T-depth 1.6-5.6e7
    i = r33.intermediates
    lo, hi = i["t_depth_per_shot"]
    _is(lo, 27 * i["t_per_rotation"] * 20500)
    _is(hi, 92 * i["t_per_rotation"] * 20500)
    assert 1.6e7 <= lo < 1.7e7 and 5.55e7 <= hi < 5.65e7                # "1.6--5.6e7"
    assert 265 <= i["f_star"][0] < 267 and 905 <= i["f_star"][1] < 907    # "270--910"
    assert all(i["baseline_ok"]) and all(i["fits_1yr"])
    _is(i["floor_wall_s"][0], 1500 * lo * 1e-5)


def test_2033_readout_circuits(r33):
    # R6: scalar conjugate setting adds one per-site QFT; fermion out-basis bounded by N(N-1)/2 Givens
    i = r33.intermediates
    assert 3.5e4 <= i["t_qft_readout"] < 4.5e4                # "about 4e4 T"
    assert 2.5e-6 <= i["t_qft_readout"] / i["t_per_shot"][0] < 3.5e-6   # "3e-6 of the shot"
    assert i["n_fermion_modes"] == 500 and i["n_givens_bound"] == 124750
    assert 7.35e6 <= i["t_fermion_readout_bound"] < 7.45e6    # "7.4e6 T"
    assert 4.5e-4 <= i["fermion_readout_fraction_bound"] < 5.5e-4   # "5e-4 of the shot"


def test_retired_r25_record():
    assert ch7.RETIRED_R25["shots_2028"] == 4000 and ch7.RETIRED_R25["n_profiles"] == 50
    a = ch7.Assumptions()
    for k in ("shots_per_mode_2028", "n_modes_2028", "shots_per_mode_2033", "n_modes_2033", "n_profiles",
              "gate_time_s"):
        assert not hasattr(a, k), k


def test_no_multi_machine_quantity_remains(a, r28, r33):
    # E27: the retired fields are a record only; nothing in the model is a machine count
    assert not hasattr(a, "n_machines")
    for r in (r28, r33):
        for k in r.intermediates:
            assert "machine" not in k and k != "scan_wall_time_weeks", k
        for n in r.notes:
            assert "machine-" not in n and "machines" not in n, n
    assert set(ch7.RETIRED_R23) == {
        "n_machines", "scan_machine_weeks", "scan_machine_years", "scan_wall_time_weeks",
        "machines_for_campaign_horizon", "machines_for_campaign_horizon_whole", "machines_for_one_year"}
    assert not set(ch7.RETIRED_R23) & set(r33.intermediates)


def test_chapter_text_assumes_one_machine():
    # E27 + R1/R6 (r25): serial single-machine times, two tiers, no machine count, no Hadamard test
    if not TEX.exists():
        pytest.skip("chapter tex not present")
    tex = TEX.read_text()
    for banned in ("machine-years", "machine-weeks", "shot-parallel", "machines over", "machines fit",
                   "would finish it in a year", "embarrassingly parallel", "-fold parallelism",
                   "50-profile", "Hadamard-test", "10 IR modes", "60 momentum bins",
                   "split into", "re-preparation", "shot survival", "rung", " arm "):
        assert banned not in tex, banned
    for s in (r"$2.9$--$5.8$ days", r"$0.71$--$2.8$~yr", r"$17$--$34$", r"$450$ shots",
              r"$1.5$--$6.0\times 10^{3}$", r"$\epsilon_l\lesssim 7\times 10^{-12}$", r"$222$--$242$",
              r"$1030$--$1050$", r"\approx 270$--$910$", r"$1.5\times 10^{10}$ T-gates",
              r"$m_\phi\Delta t=0.05$", r"$1.6$--$5.6\times 10^{7}$"):
        assert s in tex, s


def test_2033_error_rate_argument(r33):
    # "eps_l <~ 7e-12 under the report-wide fault budget" (box and prose); "three orders below" the RFI 1e-8
    i = r33.intermediates
    t = i["t_per_shot"][0]
    _is(i["eps_l"], 0.1 / t)                   # 6.73e-12 -> "<~7e-12"
    assert 6.7e-12 <= i["eps_l"] < 7e-12
    _is(i["eps_l_printed"], 7e-12)
    _is(i["fault_budget_per_shot"], 0.1)
    assert 0.1 <= i["expected_faults_at_printed_eps"] < 0.105
    assert 0.90 <= i["shot_survival"] <= 0.91  # "about 90% of shots are fault-free"
    _is(i["eps_l_rfi"], 1e-8)
    assert 1e3 <= i["eps_l_rfi_over_requirement"] < 1e4          # "three orders below"
    assert 145 <= i["expected_faults_at_rfi"] < 155               # "~150 expected faults"
    assert r33.epsilon_l == (0.1 / t, 0.1 / t)


def test_no_split_circuit_route(r33):
    # U4: splitting the shot does not reset accumulated faults; the model carries no split-circuit quantity
    for k in r33.intermediates:
        assert "circuit" not in k or k == "t_qft_readout", k
    if TEX.exists():
        tex = TEX.read_text()
        assert "Splitting the shot does not help" in tex


def test_fault_budget_observable_bias(a, r33):
    # G4: the T-only count excludes idle storage (1.6e12 qubit-cycles, 1e2 x N_T); one-site-fault bias <= 3%
    i = r33.intermediates
    assert 1.55e12 <= i["idle_qubit_cycles"] < 1.6e12            # "1.6e12 idle qubit-cycles"
    assert 100 <= i["idle_qubit_cycles_over_t"] < 110            # "1e2 times N_T"
    assert 0.15 <= i["fault_bias_dn_per_fault"] < 0.16           # "at most 0.16 per mode" (scalar)
    _is(i["fault_bias_dn_fermion"], a.fermion_fault_dn_max.lo)   # "0.12 per mode in the worst shell"
    b = i["fault_bias_fermion"]
    assert 0.022 <= b[0.5] < 0.023 and 0.11 <= b[0.1] < 0.115    # "about 2% at 0.5 and 11% at 0.1"
    assert "fault_bias_bound" not in i                           # the old "<= 3%" claim is withdrawn
    assert 0.094 <= 1 - math.exp(-0.1) < 0.096                   # "9.5% of shots faulted"
    assert 19.3 <= a.site_energy_max_quanta.lo < 19.5            # "the 19 quanta of the highest K = 16 level"


def test_2033_prep_label_is_checked(a):
    # with a 10/H_inf window the ramp is 5% of the evolution; the model must flag the printed 2.5%
    b = ch7.Assumptions(ramp_window_Hinf_2033=c.Assumed(10, "test"))
    r = ch7.model(b, "2033")
    _is(r.intermediates["prep_fraction"], 0.05)
    assert any(n.startswith("WARNING") for n in r.notes)
    assert not any(n.startswith("WARNING") for n in ch7.model(a, "2033").notes)


# ------------------------------------------------------------- cross-checks

def test_headline_ratio_2033_over_2028(r28, r33):
    # requirements summary: 2028 9.5e4 T, 2033 1.5e10 T (both headlines include prep)
    assert 1.55e5 <= r33.hard_ops[0] / r28.hard_ops[0] < 1.6e5


def test_published_is_the_box_as_printed():
    # 2028 option (c), one step at K = 8: 220-240 LQ, 9.5e4 (smaller ruling (b), 2026-10-05; was 160-180, 9.4e4)
    assert ch7.PUBLISHED["2028"].lq == (220, 240)
    assert ch7.PUBLISHED["2028"].hard_ops == (9.5e4, 9.5e4)
    assert ch7.PUBLISHED["2033"].lq == (1050, 1050)
    assert ch7.PUBLISHED["2033"].hard_ops == (1.5e10, 1.5e10)
    assert "codesign" not in ch7.PUBLISHED


def test_utility_box(a):
    u = ch7.utility(a)
    _is(u["utility_musd_per_yr"], 3.27)        # 3% x $109M/yr             app04:164
    assert u["utility_musd_per_yr_rounded"] == 3   # "~$3M/yr"              app04:159
    assert u["utility_musd_campaign"] == 15        # "~$15M over a 5-year campaign"


def test_instance_rows(a, r28, r33):
    rows28 = ch7.INSTANCE_ROWS(a, "2028", r28)
    rows33 = ch7.INSTANCE_ROWS(a, "2033", r33)
    assert len(rows28) == 1 and len(rows33) == 1
    assert rows28[0][1] == r28.lq and rows28[0][2] == r28.hard_ops
    assert "preheating" in rows33[0][0]                     # box title (app04:133), Z-D7-M06
    assert rows33[0][3]["t_bound"] is None                  # a point, no lever
    assert rows33[0][3]["conditional"] is False             # U2: step checked; the comoving-register squeeze is free
    assert ch7.INSTANCE_ROWS(a, "codesign", r33) == []


def test_ramp_step_is_stability_checked():
    # R12: the ramp Trotterizes the same H (m_phi^2 phi^2 included), so a coarse ramp step
    # at 0.5/H_inf (omega dt = 50) must be refused like a levered evolution step
    with pytest.raises(ValueError, match="ramp step"):
        ch7.Assumptions(dt_ramp_Hinf_2033=c.Assumed(0.5, "test"))


# ------------------------------------------------------------- referee checks U1-U3 (2026-10-04)

def test_digitized_site_k16_and_k4():
    # U3: "K = 16 ... about two quanta per site"; "At K = 4 the site gap is 15% low"; "K = 8 (gap error 0.1%)"
    d16, d8, d4 = ch7.digitized_site(16), ch7.digitized_site(8), ch7.digitized_site(4)
    assert abs(d16["gap"] - 1) < 1e-6 and abs(d16["x2"] - 0.5) < 1e-6
    assert 0.84 <= d4["gap"] < 0.86                          # 15% low
    assert abs(d8["gap"] - 1) < 1.5e-3                        # 0.1%
    assert 19.3 <= d16["e_max_above_ground"] < 19.5


def test_coherent_condensate_fits_k16():
    # U1/U3: "fidelity above 0.98 ... over 100 periods; 0.86 at three quanta"
    assert ch7.coherent_fidelity(16, 2.0) > 0.978
    assert 0.85 <= ch7.coherent_fidelity(16, 3.0) < 0.87


def test_step_proxy_rejects_omega_dt_1_and_accepts_0p05(a):
    # U2 (verifier rerun): set 3 (lam = 1, g = 2) at the register-consistent Phi0 = 0.31 m_phi. m_phi dt = 1
    # gives a ~27% error; 0.05 is within 2% of a fine step on n_k >~ 0.1, all of them fermion shells (reference
    # at 0.01 here to keep the test short; 0.001 in the response gives the same verdict)
    import numpy as np
    P = a.condensate_amplitude_2033.lo
    ref_s, ref_f = ch7.verlet_proxy(0.01, 1.0, P, 2.0)
    assert ref_s.max() < 0.06                                    # "scalar shells stay at n_k <~ 0.06"
    assert 0.5 <= ref_f.min() and ref_f.max() <= 1.0
    for mdt, ok in ((1.0, False), (0.05, True)):
        s, f = ch7.verlet_proxy(mdt, 1.0, P, 2.0)
        err = np.max(abs(f / ref_f - 1)[ref_f > 0.1])
        assert bool(err <= 0.02) is ok, (mdt, err)


def test_condensate_amplitude_fits_register(r33):
    # U1 (verifier): per-site condensate quanta m b^3 Phi0^2 / 2 at b = 3.46/m_phi; K = 16 holds about two
    i = r33.intermediates
    assert 1.9 <= i["condensate_quanta_per_site"] <= 2.0
    _is(i["condensate_amplitude"], 0.31)


def test_switch_occupation_formula():
    # U2 verifier: sudden switch n = (r + 1/r)/4 - 1/2, r = w_eff(1)/w_eff(0.05), w_eff = w sqrt(1 - (w h)^2/4)
    r = math.sqrt(1 - 0.25) / math.sqrt(1 - 0.05 ** 2 / 4)
    _is(ch7.switch_occupation(1.0, 1.0, 0.05), (r + 1 / r) / 4 - 0.5)
    assert ch7.switch_occupation(1.0, 0.05, 0.05) == 0.0
    assert abs(ch7.switch_occupation(1.0, 0.05, 0.05, 10)) < 1e-12   # no step change, stays in the map vacuum


def test_fermion_shots_and_scalar_precision(a, r33):
    # U2 verifier: the payoff is the fermion (binomial, m' = 12); scalar shells at n <~ 0.06 are read to ~50%
    i = r33.intermediates
    assert a.fermion_shell_modes_2033.lo == 12
    assert 0.052 <= i["fermion_n_min_priced"] < 0.054             # "n_k >~ 0.05"
    _is(i["fermion_n_min_priced_correlated"], 0.4)               # "m' = 1 ... n_k >~ 0.4"
    assert 6.3e3 <= i["scalar_shots_10pct_per_profile"] < 6.5e3  # "6e3 shots per profile"
    lo, hi = i["scalar_eps_at_priced_shots"]
    assert 0.45 <= hi < 0.47 and 0.64 <= lo < 0.66                # "45--65%"


def test_fermion_fault_shift_worst_position(a):
    # G4 verifier: an X fault mid-way along the JW order (j = 248) is the worst case on 5^3
    assert abs(ch7.fermion_fault_shift(js=[248]) - a.fermion_fault_dn_max.lo) < 5e-4
    assert ch7.fermion_fault_shift(js=[0]) < 0.08


def test_exact_reduced_step_check(a, r33):
    # U2 exact check (2026-10-04): the chapter's splitting run exactly on a 3-site ring at the actual couplings,
    # K = 16, one Wilson fermion, comoving register. Short window here (t = 12.5, set 1); the full run (t <= 1000,
    # three sets) gives <= 0.3%, Richardson ratio 4.0, and at most 1.23x the proxy on the same lattice
    o = ch7.exact_reduced_step(1.0, 1.0, [0.05, 0.025, 0.0125], "chi", T_m=12.5, rec=(12.5,))
    n1, n2, n4 = (o[h][0]["n_ferm"][1] for h in (0.05, 0.025, 0.0125))
    ex = n4 + (n4 - n2) / 3
    assert 3.9 < (n1 - n2) / (n2 - n4) < 4.1                   # second order
    assert 1.3e-3 < abs(n1 / ex - 1) < 1.5e-3                  # 1.38e-3 in the full run at t = 12.5
    assert 0.31 < n1 < 0.33
    i = r33.intermediates
    assert i["step_err_exact_reduced"] <= 3.1e-3               # "at most 0.3%"
    assert 0.0195 <= i["step_err_extrapolated"] < 0.0200       # "about 2% scaled to the 5^3 proxy"
    assert i["step_err_extrapolated"] < a.precision_campaign_2033.lo / 2
    if TEX.exists():
        tex = TEX.read_text()
        assert r"at most $0.3\%$ step error on the fermion" in tex
        assert r"at most $1.2$ times the proxy's error" in tex
        assert "still needs a convergence study" not in tex


def test_fixed_grid_fails_late(a):
    # U2 exact check: the K = 16 grid fixed at a = 1 loses the site vacuum as it narrows like a^{-3/2}
    assert 0.930 <= ch7.fixed_grid_gap(16, 4.0) < 0.935       # "by a^3 = 4 ... 7% low"
    assert ch7.fixed_grid_gap(16, 1.0) > 0.9999
    assert ch7.fixed_grid_gap(16, 256.0) < 0.03                # t = 1000
    t4 = (math.sqrt(4.0) - 1) / (1.5 * 0.01)                   # a^3 = (1 + 1.5 H t)^2 = 4
    assert 66 < t4 < 68                                        # "m_phi t ~ 70"
    assert a.fixed_grid_err_t50.lo == 0.30                    # 0.391 vs 0.300 (set 2) at m_phi t = 50
    a3_50 = (1 + 1.5 * 0.01 * 50) ** 2
    assert 3.0 < a3_50 < 3.1                                   # "a^3 \approx 3"
    if TEX.exists():
        tex = TEX.read_text()
        assert r"off by up to $30\%$" in tex
        assert "conditional on the comoving-register squeeze" not in tex
        assert "does not price" not in tex


def test_shear_encoding_matches_explicit_squeeze():
    # U2 verifier: p^2/2 + g(xp+px)/2 = U (p^2/2) U^dag - g^2 x^2/2, U = exp(-i g x^2/2), so the comoving squeeze is
    # a conjugation of the pi^2 block by ZZ phases that merge into the diagonal block (no added rotation)
    o_c = ch7.exact_reduced_step(1.0, 1.0, [0.05], "chi", T_m=12.5, rec=(12.5,))[0.05][0]
    o_s = ch7.exact_reduced_step(1.0, 1.0, [0.05], "shear", T_m=12.5, rec=(12.5,))[0.05][0]
    assert max(abs(x - y) for x, y in zip(o_c["n_ferm"], o_s["n_ferm"])) < 1e-5
    assert abs(o_c["n_scalar"] - o_s["n_scalar"]) < 1e-4
    if TEX.exists():
        assert r"-\gamma^2\chi^2/2" in TEX.read_text()


def test_2028_grid_check(a):
    # U2 for 2028 (2026-10-05): exact digitized evolution on a 4-site ring at the 2028 spacing, against the
    # continuum Gaussian through the same four-step map. A grid fixed at t_in fails within one step; refitting the
    # grid every step (rescale + shear, angle changes only) holds K = 4 for two steps and K = 8 for four at
    # m ~ 1.4 H_inf. Accuracy numbers are from the priced merged circuit (merged = True)
    fixed = [ch7.ds2028_grid_check(m, 4, "fixed", n_steps=1)[0][1] for m in (0.5, 1.5)]
    _is(fixed[0], 0.399, rel=0.01); _is(fixed[1], 2.49, rel=0.01)
    assert a.grid_2028_fixed_err_step1.value == (0.40, 2.49)          # "40--250% after one step"
    m2 = ch7.ds2028_grid_check(1.4, 4, "adapted", n_steps=2, merged=True)[-1][1]
    assert m2 < 0.05 and abs(m2 - a.grid_2028_k4_two_steps_ir.hi) < 1e-3           # "within 5% for two steps"
    assert abs(ch7.ds2028_grid_check(1.0, 4, "adapted", n_steps=2, merged=True)[-1][1]
               - a.grid_2028_k4_two_steps_ir.lo) < 1e-3
    m4 = ch7.ds2028_grid_check(1.4, 4, "adapted", merged=True)[-1][1]
    assert abs(m4 - a.grid_2028_k4_adapted_ir_end.lo) < 1e-3                         # "11--24%"
    _is(ch7.ds2028_grid_check(1.5, 4, "adapted", merged=True)[-1][1], a.grid_2028_k4_adapted_ir_end.hi, rel=0.01)
    s4 = ch7.ds2028_grid_check(1.4, 4, "adapted")                                    # split circuit, per step
    assert s4[1][1] < 0.05 and abs(s4[3][1] - m4) < 0.01                             # split ~ merged at m >= 1.4
    for m in (1.4, 1.5):
        assert max(ch7.ds2028_grid_check(m, 8, "adapted", merged=True)[-1]) < 1.5e-3 # "0.15% over four steps"
    for m, ref in ((1.0, 1.42e-3), (1.4, 1.05e-3), (1.5, 1.04e-3)):                 # option (c), the priced circuit
        e = max(ch7.ds2028_grid_check(m, 8, "adapted", n_steps=1, merged=True)[-1])
        assert e < 1.5e-3 and abs(e - ref) < 1e-5                                    # "0.15% for m ~ 1-1.4 H_inf"
    assert a.grid_2028_k8_one_step.value == (1.04e-3, 1.42e-3)
    _is(ch7.ds2028_map_vs_exact(1.0, n_steps=1), a.map_vs_desitter_2028_one_step.lo, rel=0.01)   # "13--25%"
    _is(ch7.ds2028_map_vs_exact(1.4, n_steps=1), 0.246, rel=0.01)
    assert ch7.ds2028_grid_check(1.0, 8, "adapted", merged=True)[-1][1] > 0.1        # light fields fail at 4 steps
    _is(ch7.ds2028_map_vs_exact(0.5), 0.188, rel=0.01)
    _is(ch7.ds2028_map_vs_exact(1.5), 0.868, rel=0.01)                               # "20--90% from de Sitter"
    assert ch7.ds2028_map_vs_exact(0.5, h=0.125, n_steps=16) < 0.015                 # second order: the map converges
    if TEX.exists():
        tex = TEX.read_text()
        assert r"by $40$--$250\%$ after one step" in tex
        assert r"within $5\%$ for two steps" in tex
        assert r"so the 2028 benchmark runs one step at $K=8$, which holds every shell to $0.15\%$" in tex
        assert "is an open choice" not in tex
        assert "Re-centering" not in tex
        assert r"that map is $13$--$25\%$ from continuum de Sitter" in tex
        assert "needs a reference for the digitized circuit" not in tex


def test_twopi_leading_order_comparator(a):
    # G6 (2026-10-05): Gaussian truncation (Hartree + mean-field Yukawa with fermion back-reaction, below the two-loop
    # 2PI of hep-ph/0212404) on the U2 reduced
    # instance, against the exact K = 16 runs (m_phi dt = 0.0125): exact fermion n_k at t = 12.5 is 0.3187 (set 1),
    # 0.6262 (set 3); at t = 1000, 0.9044 (set 1). Late error <= 5.3% (all sets), early up to 18%
    o1 = ch7.twopi_lo(1.0, 1.0, rec=(12.5,))[0]["n_ferm"][1]
    o3 = ch7.twopi_lo(1.0, 2.0, rec=(12.5,))[0]["n_ferm"][1]
    _is(o1, 0.2988, rel=2e-3)
    _is(o3, 0.5136, rel=2e-3)
    assert abs(o3 / 0.6262 - 1) <= a.twopi_lo_fermion_err_early.lo + 1e-3          # "18% early"
    assert a.twopi_lo_fermion_err_late.lo <= 0.055                                  # "to 5%" late
    assert a.twopi_lo_fermion_err_late.lo < a.precision_campaign_2033.lo             # inside the 10% target
    assert a.twopi_lo_scalar_err_late.lo == 0.66                                    # "up to 66% late"
    assert a.twopi_lo_scalar_err_early.lo == 0.96                                   # "96% earlier"
    ns3 = ch7.twopi_lo(1.0, 2.0, rec=(25.0,))[0]["n_scalar"]                         # exact K = 16: 0.0339
    assert abs(1 - ns3 / 0.0339) >= a.twopi_lo_scalar_err_early.lo - 0.005
    if TEX.exists():
        tex = TEX.read_text()
        assert r"we therefore treat the $5^3$ instance as a cross-check, not an advantage" in tex
        assert "Neither the two-loop truncation nor a 3D lattice was run" in tex
        assert "within reach of 2PI" not in tex
        assert "contingent on that error exceeding" not in tex
