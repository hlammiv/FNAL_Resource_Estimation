"""Ch. 8 baryogenesis: every named intermediate the prose states, so that editing one
number without the others fails loudly. Generic contract tests live in test_common.py."""

import math
from dataclasses import replace

import pytest

from estimates import ch08_baryogenesis as m
from estimates.common import CircuitStatus, Provenance, Stated, close, sci

TEX = __import__("pathlib").Path(__file__).resolve().parents[3] / "applications" / "app05_baryogenesis.tex"


@pytest.fixture(scope="module")
def a():
    return m.Assumptions()


@pytest.fixture(scope="module")
def r28(a):
    return m.model(a, "2028")


@pytest.fixture(scope="module")
def r33(a):
    return m.model(a, "2033")


def _close(x, y, tol):
    assert math.isclose(x, y, rel_tol=tol), (x, y, tol)


def _in(x, rng):
    assert rng[0] <= x <= rng[1], (x, rng)


# --- c_T derivation (R12 item 14; app05:85) ---------------------------------------

def test_scalar_counts_are_the_papers_term_counts():
    # arxiv_2407_13819 Eq. rzTrotter (:1756-1762): |Omega| (C1+C2+C3+C4) + E_D n_q^2 R_Z;
    # arxiv_2210_07985 Eq. trpi2 (:397-407): C(n_q,2) ZZ for pi^2 on the centered grid
    s16 = m.scalar_counts(4, 2, 7)
    assert s16["rot"] == {"onsite_polynomial": 15, "gradient_links": 32, "kinetic_pi2_zz": 6, "qft_controlled_phase": 2}
    assert s16["t_exact"] == {"qft_exact": 64}               # 2 x (3x3 CS + 2x8 CT + 1x7 Toffoli)
    s4 = m.scalar_counts(2, 2, 7)
    assert s4["rot"] == {"onsite_polynomial": 3, "gradient_links": 8, "kinetic_pi2_zz": 1, "qft_controlled_phase": 0}
    assert s4["t_exact"] == {"qft_exact": 6}                 # 2 x one controlled-S
    assert m.qft_cost(4, 7) == (1, 32.0)
    # a 4-qubit AQFT at the paper's formula (Eq. approxQFT_ampTrot, :1769) would cost ~5e2 T, not 32
    eps_qft = 1e-3
    l = math.log2(4 / eps_qft)
    assert 8 * 4 * l + l * math.log2(l / eps_qft) > 500


def test_fermion_counts_are_derived_here():
    f16 = m.fermion_counts(4, 2, 2, 2, 2)
    assert f16["rot"] == {"wilson_hopping": 8, "mass_wilson_onsite": 4, "yukawa_phi_bilinear": 16, "yukawa_phase_frame": 4}
    assert f16["t_exact"] == {"hopping_basis_change": 16}
    assert sum(f16["rot"].values()) == 32
    f4 = m.fermion_counts(2, 2, 2, 2, 2)
    assert sum(f4["rot"].values()) == 24 and f4["t_exact"]["hopping_basis_change"] == 16
    # 4-component spinor doubles every fermion count (rank-2 Wilson projector, 4 bilinears)
    f16x4 = m.fermion_counts(4, 2, 4, 2, 2)
    assert sum(f16x4["rot"].values()) == 64 and f16x4["t_exact"]["hopping_basis_change"] == 32


def test_ct_is_priced_by_rtol_per_instance(a, r28, r33):
    from estimates.common import t_per_rotation, eps_rot_for
    i = r28.intermediates
    assert i["rot_per_site_step"] == 55 and i["t_exact_per_site_step"] == 64
    assert i["n_rot_shot"] == 55 * 432 == 23760
    _close(i["eps_rot"], math.sqrt(1e-2 / 23760), 1e-12)
    _close(i["t_per_rotation"], 1.15 * math.log2(1 / i["eps_rot"]) + 9.2, 1e-12)
    _close(i["c_t"], 55 * i["t_per_rotation"] + 64, 1e-12)
    _close(i["c_t"], 1239.82, 1e-5)
    _close(i["c_t"], i["c_t_quoted"], 0.005)                   # printed '1.24e3'
    j = r33.intermediates
    assert j["rot_per_site_step"] == 87 and j["t_exact_per_site_step"] == 80
    assert j["n_rot_shot"] == 87 * 120 * 100
    # R13: synthesized to the deepest MLAE circuit's budget, N_rot = m_max x (N_rot(base) + priced preparation, B2)
    assert j["n_rot_budget"] == 33 * (1_044_000 + 1800 + 14400 + 472)
    _close(j["t_per_rotation"], t_per_rotation(eps_rot_for(33 * (1_044_000 + 16672))), 1e-12)
    _close(j["c_t"], 2466.43, 1e-5)
    assert sci(j["c_t"]) == ("2.5", 3)
    assert sci(j["c_t_scalar"]) == ("1.6", 3) and sci(j["c_t_fermion"]) == ("8.9", 2)   # '1.6e3 scalar, 8.9e2 fermion'
    _close(j["c_t"], j["c_t_quoted"], 0.02)
    # the base circuit's own budget (superseded) is recorded: c_T 2212.9, 2.656e7, m_max 37
    _close(j["base_budget_c_t"], 2212.94, 1e-5)
    assert j["base_budget_m_max"] == 37


def test_ct_range_is_the_instances(a, r28, r33):
    # app05:85: 'c_T runs from 2.5e2 to 2.5e3 across the instances below'
    i = r28.intermediates
    instances = (i["c_t"], i["family_8x8_nphi4_c_t"], i["family_8x8_nphi16_c_t"], i["with_fermions_6x6_c_t"],
                 i["interim_c_t"], r33.intermediates["c_t"])
    assert sci(min(instances)) == ("2.5", 2) and sci(max(instances)) == ("2.5", 3)
    assert (a.c_t_range_quoted.lo, a.c_t_range_quoted.hi) == (2.5e2, 2.5e3)


# --- 2028 box (app05:117-150) and T-count summary (app05:91) ---------------------

def test_2028_register_is_a_component_sum(r28):
    i = r28.intermediates
    assert i["qubits_per_site"] == 4
    assert i["system_qubits"] == 144
    assert i["n_anc"] == 30
    assert i["lq"] == 174
    assert r28.lq == (174, 174)


def test_2028_raw_t_is_the_derived_product(r28):
    i = r28.intermediates
    assert i["site_steps"] == 36 * 12
    _close(i["t_raw_per_shot"], i["c_t"] * 432, 1e-12)
    _close(i["t_raw_per_shot"], 5.356e5, 1e-3)                 # '= 5.4e5 raw'
    assert sci(i["t_raw_per_shot"]) == ("5.4", 5)
    assert sci(1.24e3 * 432) == ("5.4", 5)                     # printed c_T reproduces the printed raw


def test_2028_step_count_and_window(a, r28):
    i = r28.intermediates
    assert i["n_trot"] == 12 == int(a.n_trot_2028.lo)
    _close(i["t_window_over_mphi"], 1.2, 1e-9)
    assert i["n_trot_from_prose_window"] == (10.0, 20.0)
    assert i["t_window_within_prose"]


def test_2028_banked_factor_is_assumed_and_carried_exactly(a, r28):
    i = r28.intermediates
    assert a.ct_reduction_banked.prov is Provenance.ASSUMED   # R17 (c): assumed development target, not derived or cited
    assert i["ct_reduction_banked"] == 3.3
    t = r28.hard_ops[0]
    assert r28.hard_ops == (t, t)
    _close(t, i["t_raw_per_shot"] / 3.3, 1e-12)                 # 1.623e5
    assert sci(t) == ("1.6", 5)                                 # '~1.6e5'
    lo, hi = i["t_banked_range"]                                # model-internal 3-4x sanity span (not printed)
    assert (sci(lo), sci(hi)) == (("1.3", 5), ("1.8", 5))
    assert i["t_box_quoted"] == 1.6e5 and i["t_box_within_banked_range"]
    assert sci(i["t_over_cap"], 1) == ("1.6", 0)                # T-count summary (r16: moved from box): '1.6 times the 1e5 cap'


def test_2028_box_row_closes(a, r28):
    pub = m.PUBLISHED["2028"]
    assert pub.rel_tol == 0.05
    assert close(r28.hard_ops, pub.hard_ops, pub.rel_tol)
    assert close(r28.lq, pub.lq, pub.rel_tol)


def test_2028_breakdown_is_primitive_level_and_compiled(r28):
    tot = r28.breakdown_total()
    _close(tot, r28.hard_ops[0], 1e-9)
    priced = [p for p in r28.breakdown if p.t_each > 0]
    assert {p.name for p in priced} == {"scalar_onsite_polynomial", "scalar_gradient_links", "scalar_kinetic_pi2_zz",
                                        "scalar_qft_controlled_phase", "scalar_qft_exact"}
    assert all(p.status is CircuitStatus.COMPILED for p in priced)
    assert all("derived here" in p.src or "arxiv_" in p.src for p in priced)
    holes = [p for p in r28.breakdown if p.t_each == 0]
    assert holes and all(p.status is CircuitStatus.UNSOURCED for p in holes)


def test_2028_shots_and_wall_time(r28):
    i = r28.intermediates
    assert i["shots_within_n_rare_range"]
    # R17 g_rate: wall = shots x (T x 1 us + t0), t0 = 0.1 ms; gate-limited, no 1 s^-1 shot rate
    assert i["shot_overhead_s"] == 1e-4
    # r23 single-time record (no longer printed): 5e3 shots at one time, 13.5 min (44.6 raw)
    _close(i["wall_min"] * 60, 5e3 * (1.6230e5 * 1e-6 + 1e-4), 1e-4)   # 812.0 s
    assert round(i["wall_min"]) == 14 and round(i["wall_min_raw"]) == 45
    assert i["shot_overhead_fraction"] < 1e-3
    assert round(i["gate_time_min"]) == 14 and round(i["gate_time_min_raw"]) == 45


def test_2028_slope_two_times(a, r28):
    # R25 (R2 + shot audit): P_FV at 4 and 12 steps, 3 x 5e3 = 1.5e4 shots at each; each circuit at its own depth
    i = r28.intermediates
    assert i["slope_early_steps"] == 4
    assert i["slope_shots_per_time"] == 1.5e4 and r28.shots == (3e4, 3e4)
    early = m.price_instance(a, 36, 4, 16, 0, 2, reduction=3.3)
    _close(i["slope_t_early"], early["t_raw"] / 3.3, 1e-12)
    assert sci(i["slope_t_early"]) == ("5.2", 4)                # '5.2e4 T'
    wall = 1.5e4 * ((r28.hard_ops[0] + i["slope_t_early"]) * 1e-6 + 2e-4)
    _close(i["slope_wall_min"] * 60, wall, 1e-12)
    _close(r28.wall_time_s[0], wall, 1e-12)
    assert round(i["slope_wall_min"]) == 54                     # '~54 min per parameter point'
    assert round(i["slope_wall_min_raw"] / 60, 1) == 2.9        # '~2.9 h at the unreduced c_T'


def test_2028_depth_exports(a, r28):
    # R9: one tier; T-depth from factory.json (banked 5.4e3-7.0e3), F* 23-30 >= 10
    i = r28.intermediates
    assert i["t_per_shot"] == r28.hard_ops
    assert i["t_depth_per_shot"] == (5.4e3, 7.0e3)
    assert (round(i["f_star"][0]), round(i["f_star"][1])) == (23, 30)   # 'F* = 23-30'
    assert i["baseline_ok"] == (True, True)
    _close(i["floor_wall_s"][0], i["slope_shots_full_equiv"] * 5.4e3 * 1e-5, 1e-12)
    assert i["wall_first_result_s"] is None
    _close(i["wall_campaign_s"], r28.wall_time_s[0], 1e-12)
    assert i["factories_for_1yr"][0] < 1


def test_wilson_mass_matches_fermion_primitives():
    # R17 (e): Ch.8's mass + Wilson on-site count equals FP section_mass.tex 2N at d = 2 (N = N_c N_f = 2)
    from estimates.groups import wilson_mass_rotations
    fc = m.fermion_counts(n_q=4, n_f=2, components=2, links=2, dims=2)
    assert fc["rot"]["mass_wilson_onsite"] == wilson_mass_rotations(1, 2, 2) == 4


def test_2028_family_rows(r28):
    # r16 sweep: the family numbers moved from the 2028 box to the T-count summary prose
    i = r28.intermediates
    assert i["family_8x8_nphi4_lq"] == 158
    _close(i["family_8x8_nphi4_t_banked"], i["family_8x8_nphi4_t_raw"] / 3.3, 1e-12)
    assert sci(i["family_8x8_nphi4_t_banked"]) == ("5.9", 4)    # '5.9e4 T'
    assert i["family_8x8_nphi16_lq"] == 286
    _close(i["family_8x8_nphi16_t_banked"], i["family_8x8_nphi16_t_raw"] / 3.3, 1e-12)
    assert sci(i["family_8x8_nphi16_t_banked"]) == ("2.9", 5)   # '2.9e5 T'


def test_2028_fermion_arm_numbers_are_facts(r28):
    # app05:91 (R12 light touch): the fermion arm enters at N_phi=16 and waits for 2033 (resolution);
    # '288 system qubits and c_T ~ 2.0e3, or 8.5e5 T raw and ~2.6e5 T with the same reduction applied to every sector'
    i = r28.intermediates
    assert i["with_fermions_6x6_system_qubits"] == 288
    assert i["with_fermions_6x6_rot_per_site_step"] == 87
    assert round(i["with_fermions_6x6_c_t"], -2) == 2000
    assert sci(i["with_fermions_6x6_t_raw"]) == ("8.5", 5) and i["with_fermions_6x6_t_raw_quoted"] == 8.5e5
    _close(i["with_fermions_6x6_t_banked_all_sectors_point"], i["with_fermions_6x6_t_raw"] / 3.3, 1e-12)
    assert sci(i["with_fermions_6x6_t_banked_all_sectors_point"]) == ("2.6", 5)
    assert i["with_fermions_6x6_t_quoted"] == 2.6e5


def test_interim_variant(r28):
    i = r28.intermediates
    assert i["interim_qubits_per_site"] == 6
    assert i["interim_system_qubits"] == 96
    assert i["interim_phase_register_bits"] == 13
    assert i["interim_workspace"] == 17 == i["interim_workspace_from_prior"]
    assert i["interim_lq"] == 113 and i["interim_lq_quoted"] == 115
    assert i["interim_lq_with_2028_ancilla"] == 143
    assert i["interim_two_register_lq_quoted"] == 190
    assert i["interim_two_register_lq_candidates"] == (192, 130)
    # c_T ~ 7.5e2: 2.5e2 scalar + 5.0e2 fermions; x 16 x 12 = 1.449e5, printed '~1.4e5'; gap '~1.4'
    assert i["interim_site_steps"] == 192 and i["interim_rot_per_site_step"] == 36
    assert sci(i["interim_c_t"]) == ("7.5", 2)
    assert sci(i["interim_c_t_scalar"]) == ("2.5", 2) and round(i["interim_c_t_fermion"], -1) == 500
    _close(i["interim_t_from_inputs"], i["interim_c_t"] * 192, 1e-12)
    assert sci(i["interim_t_from_inputs"]) == ("1.4", 5) and i["interim_t_quoted"] == 1.4e5
    assert sci(i["gap_factor_vs_cap"]) == ("1.4", 0) and i["gap_factor_vs_cap_quoted"] == 1.4
    assert sci(i["interim_t_banked_all_sectors"]) == ("4.4", 4) and i["interim_fits_cap_with_banked_reduction"]
    assert sci(i["interim_t_banked_scalar_only"]) == ("1.1", 5)
    assert not i["interim_fits_cap_scalar_share_only"]
    assert i["interim_scalar_share_le_benchmark"]               # 'below the benchmark's 1.6e5'


# --- 2033 box (app05:152-174) and T-count summary (app05:91) ---------------------

def a_generic_range():
    g = m.Assumptions().n_anc_generic
    return (g.lo, g.hi)


def test_2033_register_is_a_component_sum(r33):
    i = r33.intermediates
    assert i["volume"] == 120
    assert i["qubits_per_site"] == 8
    assert i["system_qubits"] == 960
    assert i["n_anc"] == 40
    assert i["n_anc_within_generic_range"]
    assert a_generic_range() == (10, 40)
    assert i["workspace_anc_implied"] == 10
    assert i["lq"] == 1000 and r33.lq == (1000, 1000)
    assert i["four_component_system_qubits"] == 1440
    _close(i["lz_over_mphi"], 7.5, 1e-9)


def test_2033_t_is_the_derived_product(r33):
    i = r33.intermediates
    assert i["n_trot"] == 100
    _close(i["t_per_step"], i["c_t"] * 120, 1e-12)
    assert sci(i["t_per_step"]) == ("3", 5)                   # '3.0e5/step'
    t = i["t_per_shot"][0]
    assert i["t_per_shot"] == (i["base_circuit_t"], i["base_circuit_t"])
    _close(t, i["t_per_step"] * 100, 1e-12)
    assert sci(t) == ("3", 7)                                 # 'Base circuit 3.0e7 T'
    _close(t, 2.9597e7, 1e-4)                                 # evolution
    assert 0.02 < i["fraction_of_envelope"] < 0.03              # base circuit / 1e9 (no longer printed)
    # B2 (2026-10-05): one application of A = evolution + priced free preparation, 'Base circuit 3.0e7'
    _close(i["t_A"], t + i["free_prep_t"], 1e-12)
    assert sci(i["t_A"]) == ("3", 7)
    # G2: the headline is the deepest executed circuit, 33 applications of A + 16 reflection pairs, '>~9.9e8'
    assert r33.hard_ops == (i["deepest_mlae_circuit_t"],) * 2 == (i["headline_t"],) * 2
    _close(r33.hard_ops[0], 33 * i["t_A"] + 16 * i["mlae_reflection_t_per_iterate"], 1e-12)
    assert sci(r33.hard_ops[0]) == ("9.9", 8) and r33.hard_ops[0] <= 1e9
    _close(r33.breakdown_total(), r33.hard_ops[0], 1e-9)
    assert close(r33.hard_ops, m.PUBLISHED["2033"].hard_ops, m.PUBLISHED["2033"].rel_tol)
    assert m.PUBLISHED["2033"].hard_ops == (9.9e8, 9.9e8)
    # shots are depth-33 circuits, so shots x headline = (queries / 33) x headline
    _close(r33.shots[0] * 33, i["queries_total"], 1e-9)
    # 'Maximum circuit time ... ~17 min (2033 deepest MLAE circuit)'
    assert round(i["deepest_circuit_s"] / 60) == 17


def test_2033_breakdown_shows_the_unpriced_pieces(r33):
    names = {p.name: p for p in r33.breakdown}
    for hole in ("gibbs_state_prep", "lindblad_dissipator_step", "vacuum_prep_dressing"):
        assert names[hole].t_each == 0 and names[hole].status is CircuitStatus.UNSOURCED
    # B2 (2026-10-05): the free preparation and the reflections are priced in the deepest circuit
    for priced in ("vacuum_prep_scalar_product", "vacuum_prep_slater_ky", "vacuum_prep_fft_transverse", "packet_prep",
                   "mlae_reflections"):
        assert names[priced].t_each > 0 and names[priced].status is CircuitStatus.COMPILED
    assert names["mlae_reflections"].count == 16 and names["vacuum_prep_slater_ky"].count == 33 * 14400
    # the deepest circuit is 33 base circuits (G2)
    fermion = sum(p.t_total for p in r33.breakdown if p.name.startswith("fermion_"))
    scalar = sum(p.t_total for p in r33.breakdown if p.name.startswith("scalar_"))
    _close(fermion, r33.intermediates["c_t_fermion"] * 12000 * 33, 1e-9)
    _close(scalar, r33.intermediates["c_t_scalar"] * 12000 * 33, 1e-9)
    assert all(p.status is CircuitStatus.COMPILED for p in r33.breakdown if p.t_each > 0)


def test_2033_m_max_is_the_fixed_point(r33):
    # R13: m_max = largest m with m T_A(m) + (m-1)/2 refl <= 1e9, each depth priced at its own budget; downward
    # scan from m0 = floor((1e9 + refl/2) / (T_A(1) + refl/2)) = 37 (T_A(1) = 2.6980e7)
    i = r33.intermediates
    it = i["m_max_iterations"]
    assert [x[0] for x in it] == [37, 36, 35, 34, 33] and [x[2] for x in it] == [False] * 4 + [True]
    _close(it[0][1] / 37 - 18 / 37 * i["mlae_reflection_t_per_iterate"], 3.0157e7, 1e-4)   # T_A(37)
    _close(it[-1][1], i["deepest_mlae_circuit_t"], 1e-12)
    t, refl = i["t_A"], i["mlae_reflection_t_per_iterate"]
    _close(t, 3.0056e7, 1e-4)
    assert math.floor((1e9 + refl / 2) / (t + refl / 2)) == i["m_max"]
    assert 34 * t + 16 * refl > 1e9 >= i["deepest_mlae_circuit_t"]
    # B3: the depth-33 circuit leaves less than one Trotter step of room per application for the dressing ramp
    assert 0.5 < i["ramp_steps_left_at_depth33"] < 1
    assert round(i["prep_r_depth33_exceeds_envelope"], 3) == 0.007          # 'any r > 0.007' 


def test_2033_mlae_schedule(r33):
    i = r33.intermediates
    assert i["m_max"] == 33 == i["m_max_quoted"]                # 'm_max = floor(1e9/3.0e7) = 33'
    assert math.floor(1e9 / 2.9583e7) == 33
    assert i["mlae_depths"] == (33,)                           # R16: every circuit at depth 33
    assert sci(i["per_bin_mlae"]) == ("7", 4)                 # '~7.0e4 base-circuit queries per kinematic bin'
    _close(i["per_bin_mlae"], i["per_bin_mlae_quoted"], 0.02)
    _close(i["per_bin_mlae"], 1e6 / (33 * i["flag_contrast_p0"] ** 2), 1e-12)   # B2 flag contrast
    _close(i["per_bin_unit_contrast_r26"], 1e6 / 33, 1e-12)
    assert sci(i["queries_total"]) == ("7", 6)                # box and :83: 7.0e6
    assert sci(i["queries_total"], 0) == ("7", 6)               # Eq. (Nshot_baryo): 7e6
    assert sci(i["mlae_circuits_per_bin"]) == ("2.1", 3)      # prose: '~2.1e3 circuits' per bin
    assert sci(i["mlae_circuits_per_bin"] * i_bins(i)) == ("2.1", 5)   # box: '2.1e5 circuits'
    assert round(i["mlae_overhead_vs_heisenberg"]) == 70
    _close(i["per_bin_shot_noise"], 1e6, 1e-9)
    _close(i["per_bin_heisenberg"], 1e3, 1e-9)
    _close(i["queries_heisenberg_total"], 1e5, 1e-9)
    _close(r33.shots[0], i["mlae_circuits_total"], 1e-9)          # G2: shots are depth-33 circuits
    _close(i["queries_total"] * i["flag_contrast_p0"] ** 2, i["queries_total_single_depth_r13"], 1e-12)
    # LIS sensitivity (not priced): 1.46x the queries
    _close(i["mlae_lis_queries_total"] / i["queries_total"], 33 * 289 / 6545, 1e-12)
    assert f"{i['mlae_lis_query_constant']:.2f}" == "1.46"    # prose: '1.46 times as many queries'


def i_bins(i):
    return i["queries_total"] / i["per_bin_mlae"]


# --- MLAE estimator constant (R15) -------------------------------------------------------

def test_lis_constant_is_suzukis_fisher_ratio():
    # Suzuki :271: LIS with M = 16 iterates: I = N (2M+3)(2M+1)(M+1)/(3 a(1-a)), N_q = N (M+1)^2
    big_m = 16
    fisher = (2 * big_m + 3) * (2 * big_m + 1) * (big_m + 1) // 3
    nq = (big_m + 1) ** 2
    assert fisher == 6545 and nq == 289
    _close(m.lis_query_constant(33), 33 * nq / fisher, 1e-15)   # 1.45714, printed 1.46
    assert f"{m.lis_query_constant(33):.2f}" == "1.46"
    _close(m.lis_circuit_constant(33), 17 * 33 ** 2 / 6545, 1e-15)   # 2.83 (circuits)
    assert m.lis_query_constant(1) == 1.0
    # Ch. 3's depths (r14), in queries, for the record: 13 -> 1.4, 5 -> 1.286
    _close(m.lis_query_constant(13), 13 * 49 / 455, 1e-15)
    _close(m.lis_query_constant(5), 5 * 9 / 35, 1e-15)
    _close(m.lis_circuit_constant(13), 2.6, 1e-12)


def test_balanced_encoding_gives_constant_one_at_a_single_depth():
    # a = (1+D)/2, p(D) = sin^2(m theta), theta = arcsin sqrt(a): Fisher per circuit for D is m^2/(1-D^2)
    for mm in (1, 5, 33):
        for d in (1e-3, 0.02, -0.03):
            h = 1e-7
            p = lambda x: math.sin(mm * math.asin(math.sqrt((1 + x) / 2))) ** 2
            dp = (p(d + h) - p(d - h)) / (2 * h)
            fisher = dp ** 2 / (p(d) * (1 - p(d)))
            _close(fisher, mm ** 2 / (1 - d ** 2), 1e-5)


def test_branch_condition_and_author_prior(a, r33):
    i = r33.intermediates
    b = m.balanced_branch_bound(33)
    _close(b, math.sin(math.pi / 66), 1e-15)
    assert f"{b:.3f}" == "0.048"                               # printed 'sin(pi/66) = 0.048'
    # m theta just below / above the bound stays inside / leaves the branch [8 pi, 8.5 pi]
    for d, inside in ((0.999 * b, True), (1.001 * b, False)):
        th = math.pi / 4 + 0.5 * math.asin(d)
        assert (8 * math.pi < 33 * th < 8.5 * math.pi) is inside
    # R15: no derived bound; R16: the authors adopt |A_R| < sin(pi/66) = 0.048 as the prior
    assert a.dn_state_independent_bound.prov is Provenance.ASSUMED and a.dn_state_independent_bound.lo == 1.0
    assert i["mlae_derived_bound_ok"] is False and i["mlae_author_prior_ok"] is True
    assert a.c_est.prov is Provenance.STATED and a.c_est.lo == 1.0
    assert a.dn_branch_prior.prov is Provenance.STATED and f"{b:.3f}" == f"{a.dn_branch_prior.lo:.3f}"
    assert i["mlae_single_depth"] is True and i["mlae_estimator_constant"] == 1.0
    assert round(i["dn_estimate_headroom"]) == 48              # 'sits 48 times inside it'
    # without the prior: 'below Delta theta_C' taken as a bound (0.1) is ~2x the branch condition, so LIS
    no_prior = replace(a, mlae_schedule=m.Cited("LIS", "test"), c_est=Stated(m.lis_query_constant(33), "test"))
    r01 = m.model(replace(no_prior, dn_state_independent_bound=Stated(0.1, "test", "test")), "2033").intermediates
    assert r01["mlae_single_depth"] is False and round(0.1 / b) == 2
    _close(r01["mlae_estimator_constant"], 33 * 289 / 6545, 1e-12)
    _close(r01["per_bin_mlae"], i["mlae_lis_per_bin"], 1e-12)
    assert r01["mlae_depths"] == tuple(range(1, 34, 2))
    # a derived bound under 0.048 would also give constant one
    r1 = m.model(replace(no_prior, dn_state_independent_bound=Stated(0.01, "test", "test")), "2033").intermediates
    assert r1["mlae_single_depth"] is True and r1["mlae_estimator_constant"] == 1.0
    _close(r1["per_bin_mlae"], 1e6 / (33 * i["flag_contrast_p0"] ** 2), 1e-12)


def test_lis_resolves_the_branch(r33):
    i = r33.intermediates
    # LIS sensitivity (R15 record; 6.5x before the B2 flag contrast): the ladder resolves the branch with ~9.8x margin
    assert round(i["mlae_lis_branch_margin_sigma"], 1) == 9.8
    # the binding step is the first one: m = 1 shots, next depth 3
    n = i["mlae_lis_circuits_per_depth_per_bin"]
    assert round(n, -1) == 350
    _close(i["mlae_lis_branch_margin_sigma"], (math.pi / 12) * 2 * math.sqrt(n), 1e-12)


def test_2033_wall_time_and_campaign(r33):
    # E27 (r23): one machine everywhere; every wall time is serial, shots x (T x 1 us + t0)
    i = r33.intermediates
    assert round(i["base_circuit_s"]) == 30                     # base circuit '~30 s'
    assert i["shot_overhead_s"] == 1e-4
    _close(i["mlae_circuits_total"], i["mlae_circuits_per_bin"] * 100, 1e-12)   # 2.1e5 circuits
    wall = i["mlae_circuits_total"] * (i["deepest_mlae_circuit_t"] * 1e-6 + 1e-4)
    _close(r33.wall_time_s[0], wall, 1e-12)
    _close(i["dn_arm_yr"] * 365.25 * 86400, wall, 1e-12)
    assert i["mlae_circuits_total"] * 1e-4 / wall < 1e-6        # 'the ~0.1 ms per-shot overhead is negligible'
    assert round(i["dn_arm_yr"], 1) == 6.7                      # '~6.7 yr per point'
    # r23 single-time nucleation record (no longer printed): 1e4 x 29.6 s = 3.4 d
    _close(i["nucleation_point_days_single_time_r23"] * 86400, 1e4 * (i["base_circuit_t"] * 1e-6 + 1e-4), 1e-12)
    # R25 slope: 3e4 shots at 33 and at 100 steps, each circuit at its own synthesis budget
    assert i["nucleation_shots_per_time"] == 3e4 and i["nucleation_early_steps"] == 33
    _close(i["nucleation_t_late_own_budget"], i["base_budget_t_per_shot"], 1e-12)      # 2.656e7
    assert sci(i["nucleation_t_early"]) == ("8.4", 6)            # '8.4e6 T'
    _close(i["nucleation_point_days"] * 86400,
           3e4 * ((i["nucleation_t_late_own_budget"] + i["nucleation_t_early"]) * 1e-6 + 2e-4), 1e-12)
    assert round(i["nucleation_point_days"]) == 12              # '~12 days per nucleation point'
    _close(i["dn_campaign_serial_yr"], 10 * i["dn_arm_yr"], 1e-12)
    assert round(i["dn_campaign_serial_yr"]) == 67              # 'the ten Delta n points ~67 yr'
    _close(i["nucleation_scan_serial_yr"], 100 * i["nucleation_point_days"] / 365.25, 1e-12)
    assert round(i["nucleation_scan_serial_yr"], 1) == 3.3      # 'the nucleation scan ~3.3 yr'
    _close(i["campaign_serial_yr"], i["dn_campaign_serial_yr"] + i["nucleation_scan_serial_yr"], 1e-12)
    assert round(i["campaign_serial_yr"]) == 70                 # '~70 yr on one machine'
    assert i["horizon_yr"] == 5
    _close(i["campaign_over_horizon"], i["campaign_serial_yr"] / 5, 1e-12)
    assert round(i["campaign_over_horizon"]) == 14              # '14 times the 5-year horizon'
    # what fits in 5 years on one machine: the nucleation scan and the first result, not a campaign Delta-n point
    assert i["dn_points_in_horizon"] == 0
    _close(i["fits_in_horizon_yr"], i["nucleation_scan_serial_yr"] + 2 * i["first_dn_point_days"] / 365.25, 1e-12)
    assert round(i["fits_in_horizon_yr"], 1) == 4.7             # '(~4.7 yr)'
    assert i["nucleation_scan_serial_yr"] + i["dn_arm_yr"] > 5
    # LIS sensitivity (not printed): 9.7 yr per point, 100 yr serial, 20x the horizon
    assert round(i["mlae_lis_dn_arm_yr"], 1) == 9.7
    assert round(i["mlae_lis_campaign_serial_yr"]) == 100
    assert round(i["mlae_lis_campaign_over_horizon"]) == 20


def test_2033_first_result_tier(a, r33):
    # R1/R2: 3-sigma detection of Delta n at T_c and 0.9 T_c, plus the nucleation slope at both temperatures
    i = r33.intermediates
    assert i["first_points"] == 2 and i["first_detection_sigma"] == 3
    assert i["s_int_taken"] == 1e-3 and a.s_int.prov is Provenance.ASSUMED      # not computed (NEEDS_AUTHOR)
    # B2: label-controlled packets lengthen A, so the first-result circuits run at depth 31; queries x 1/p0^2
    assert i["first_m_top"] == 31 and sci(i["label_packets_t"]) == ("1.7", 6)
    _close(i["first_queries_per_point"], 9 / (1e-6 * 31 * i["flag_contrast_p0"] ** 2), 1e-12)
    assert sci(i["first_queries_per_point"]) == ("6.7", 5)      # '6.7e5 queries'
    assert sci(i["first_circuits_per_point"]) == ("2.2", 4)
    assert round(i["first_dn_point_days"], -1) == 250           # '~250 days'
    _close(i["first_result_days"], 2 * (i["first_dn_point_days"] + i["nucleation_point_days"]), 1e-12)
    assert round(i["first_result_yr"], 1) == 1.4                # '~1.4 yr'
    _close(i["wall_first_result_s"], i["first_result_days"] * 86400, 1e-12)
    # the Delta-n wall scales as (1e-3 / S_int)^2
    half = m.model(replace(a, s_int=Stated(5e-4, "test")), "2033").intermediates
    _close(half["first_dn_point_days"], 4 * i["first_dn_point_days"], 1e-6)
    assert i["first_label_register"] == 7 and i["first_workspace_left"] == 3


def test_2033_depth_exports(r33):
    # R9: T-depth per base circuit from factory.json, 7.2e5 (A = 40) to 2.87e6 (A = 10); F* 10-41
    i = r33.intermediates
    assert i["t_depth_per_shot"] == (7.2e5, 2.87e6)
    assert (round(i["f_star"][0]), round(i["f_star"][1])) == (10, 41)   # 'F* = 10-41'
    assert i["baseline_ok"] == (True, True)
    _close(i["wall_campaign_s"], i["campaign_serial_yr"] * 365.25 * 86400, 1e-12)
    _close(i["campaign_serial_yr_check"][0], i["campaign_serial_yr"], 1e-4)   # full-circuit equivalents agree
    _close(i["floor_wall_s"][1], i["campaign_shots_full_equiv"] * 2.87e6 * 1e-5, 1e-12)
    assert round(i["factories_for_1yr"][0]) == 700 > i["f_star"][1]          # a 1-yr campaign is out of reach
    # factory.json had 28 and 2.76 yr at A = 10 before the B2 contrast and preparation (x 2.35)
    assert round(i["dn_point_factories_for_1yr"][0]) == 67
    _close(i["dn_point_floor_wall_s"][1] / (365.25 * 86400), 6.46, 0.01)


def test_2033_first_result_depth(r33):
    # verifier r25 issue 1: the 7-qubit label leaves A = 3; borrowing the bath gives A = 33
    i = r33.intermediates
    assert i["first_t_depth_per_shot"] == (8.7e5, 9.6e6)
    _close(8.7e5, 7.2e5 * 317 / 261, 0.01)                                # ceil(10440/33) / ceil(10440/40)
    _close(9.6e6, 2.87e6 * 3480 / 1044, 0.01)                             # ceil(10440/3) / ceil(10440/10)
    assert round(i["first_f_star"][1]) == 34                              # 'F* ~ 34'
    assert round(i["first_f_star"][0], 1) == 3.1 and i["first_baseline_ok"] == (False, True)
    assert round(i["first_result_yr_unborrowed"], 1) == 4.4              # '~4.4 yr' without borrowing
    import pathlib
    tex = (pathlib.Path(__file__).resolve().parents[3] / "applications" / "app05_baryogenesis.tex").read_text()
    assert r"T-depth $\sim 8.7\times 10^{5}$, $F^*\approx 34$" in tex
    assert r"the first result takes $\sim 4.4$ yr" in tex
    assert r"(1.5\text{--}3)\times 10^{4}\text{ per time}" in tex
    assert "rare-event nucleation budget" not in tex


def test_no_machine_count_anywhere(r28, r33):
    # E27 ruling (H. Lamm 2026-10-02): the report never assumes more than one machine
    for r in (r28, r33):
        assert not [k for k in r.intermediates if "machine" in k]
    import pathlib, re
    tex = (pathlib.Path(__file__).resolve().parents[3] / "applications" / "app05_baryogenesis.tex").read_text()
    assert not re.search(r"machine-yr|machine-year|\d+\s+machines|shot-parallel|machines (inside|needed|over)", tex)
    assert "on one machine" in tex


def test_2033_unbanked_reduction_is_visible(r33):
    _close(r33.intermediates["t_per_shot_if_reduction_banked_point"], r33.intermediates["base_circuit_t"] / 3.3, 1e-12)


# --- R3: 0.1 expected faults per shot --------------------------------------------

def test_epsilon_l_follows_r3(a, r28, r33):
    assert float(a.fault_budget_per_shot.lo) == 0.1 and a.fault_budget_per_shot.prov is Provenance.ASSUMED
    e28 = r28.intermediates["epsilon_l_required"]
    _close(e28, 0.1 / r28.hard_ops[0], 1e-12)
    assert r28.epsilon_l == (e28, e28)
    assert sci(e28, 0) == ("6", -7) and r28.intermediates["epsilon_l_quoted"] == 6e-7
    i = r33.intermediates
    _close(i["epsilon_l_required_base_circuit"], 0.1 / i["base_circuit_t"], 1e-12)
    _close(i["deepest_mlae_circuit_t"], 33 * i["t_A"] + 16 * i["mlae_reflection_t_per_iterate"], 1e-12)
    _close(i["epsilon_l_required_deepest_mlae_circuit"], 0.1 / r33.hard_ops[0], 1e-12)
    assert i["deepest_mlae_circuit_t"] <= 1e9
    assert r33.epsilon_l == (i["epsilon_l_required_deepest_mlae_circuit"], i["epsilon_l_required_base_circuit"])
    assert sci(i["epsilon_l_required_deepest_mlae_circuit"], 0) == ("1", -10)
    assert sci(i["epsilon_l_required_base_circuit"], 0) == ("3", -9)
    # r16 sweep: the box prints only the deepest-circuit eps_l (R-TOL); the base-circuit value is computed above, no longer quoted
    assert i["epsilon_l_quoted_deepest_mlae_circuit"] == 1e-10


def test_utility_box_has_no_unsourced_dollar_figure(a, r33):
    # referee G5 ruling (2026-10-04): the unsourced LISA dollar figure is gone from the model and the box
    assert "utility_usd" not in r33.intermediates and not hasattr(a, "lisa_us_contribution_usd")
    from pathlib import Path
    p = Path(__file__).resolve().parents[3] / "applications" / "app05_baryogenesis.tex"
    if not p.exists():
        pytest.skip("chapter .tex not on disk")
    box = p.read_text().split("begin{utilitybox}")[1].split("end{utilitybox}")[0]
    assert "\\$" not in box and "several-hundred-million" not in box
    assert "LISA" in box


# --- scaling hooks ---------------------------------------------------------------

def test_2033_t_grows_with_volume_near_linearly(a, r33):
    # 'T-gates scale with volume, linearly apart from the slow log2 growth of the per-rotation cost'
    big = m.model(replace(a, lz_sites=Stated(30, "x")), "2033")
    assert big.intermediates["volume"] == 240
    ratio = big.intermediates["base_circuit_t"] / r33.intermediates["base_circuit_t"]
    # at the deepest-circuit budget m_max halves (33 -> 16), which offsets the log2 growth: 1.998
    assert big.intermediates["m_max"] == 16
    assert 1.95 < ratio < 2.05
    assert big.lq == (2 * 960 + 40, 2 * 960 + 40)


def test_2033_four_component_spinor(a, r33):
    four = m.model(replace(a, spinor_components=Stated(4, "x")), "2033")
    assert four.intermediates["system_qubits"] == 1440 and four.lq == (1480, 1480)
    assert four.intermediates["rot_per_site_step"] == 55 + 64
    assert four.intermediates["c_t_fermion"] > 1.9 * r33.intermediates["c_t_fermion"]


def test_envelope_below_per_shot_t_is_refused(a):
    with pytest.raises(ValueError):
        m.model(replace(a, envelope_2033=Stated(1e7, "x")), "2033")


# --- provenance and structure ------------------------------------------------------

def test_uses_the_report_pricing_helpers():
    import inspect
    src = inspect.getsource(m)
    assert "t_per_rotation(" in src and "eps_rot_for(" in src and "toffoli_t(" in src


def test_model_does_not_read_published():
    import inspect
    for fn in (m._model_2028, m._model_2033):
        assert "PUBLISHED" not in inspect.getsource(fn)


def test_provenance_is_honest(a):
    assert a.t_2028_box.prov is Provenance.STATED
    assert a.ct_reduction.prov is Provenance.ASSUMED
    assert a.ct_reduction_banked.prov is Provenance.ASSUMED   # R17 (c): assumed development target, not derived or cited
    assert a.toffoli_convention.prov is Provenance.CITED and a.toffoli_convention.value == "textbook"
    assert a.eps_syn.prov is Provenance.CITED and a.eps_syn.value == 1e-2
    assert a.links_per_site.prov is Provenance.ASSUMED
    assert a.envelope_2033.prov is Provenance.CITED


def test_codesign_is_refused(a):
    with pytest.raises(ValueError):
        m.model(a, "codesign")


def test_instance_rows(a, r28, r33):
    rows28 = m.INSTANCE_ROWS(a, "2028", r28)
    assert len(rows28) == 4
    for _, lq, t, extra in rows28:
        assert lq[0] <= lq[1] and t[0] <= t[1]
        assert {"t_bound", "conditional", "codesign", "reconciled"} <= set(extra)
    main = rows28[0]
    assert main[1] == (174, 174) and main[2] == r28.hard_ops and main[3]["reconciled"]
    interim = [row for row in rows28 if "interim" in row[0]][0]
    assert interim[1] == (113, 113)
    assert interim[2] == (r28.intermediates["interim_t_from_inputs"],) * 2
    assert interim[3]["conditional"]
    rows33 = m.INSTANCE_ROWS(a, "2033", r33)
    assert len(rows33) == 1
    lbl, lq, t, extra = rows33[0]
    assert (lbl, lq, t) == (r"wall scattering $15{\times}8$", (1000, 1000), r33.hard_ops)
    # G2/B3: the deepest circuit is a lower bound while preparation is unpriced -> conditional mark
    assert {k: extra[k] for k in ("t_bound", "conditional", "codesign", "reconciled")} == \
        dict(t_bound=None, conditional=True, codesign=False, reconciled=True)


def test_assumption_validation_rejects_nonsense():
    for bad in (dict(n_phi=Stated(12, "x")), dict(n_phi=Stated(1, "x")), dict(ct_reduction=Stated((4, 3), "x")),
                dict(ct_reduction_banked=Stated(2.5, "x")), dict(fault_budget_per_shot=Stated(1.5, "x")),
                dict(shots_2028=Stated(0, "x")), dict(v_2028=Stated(0, "x")),
                dict(n_bath_anc_2033=Stated(45, "x")), dict(spinor_components=Stated(3, "x")),
                dict(wilson_r=Stated(0.5, "x")), dict(links_per_site=Stated(0, "x"))):
        with pytest.raises(ValueError):
            replace(m.Assumptions(), **bad)


# --- referee revision 2026-10-04: B2 construction, B3 preparation sensitivity, editorial max circuit time ---

def test_b2_construction_priced(r33):
    # app05:83 (2026-10-05): S_chi counts <= 480 occupations into a 9-bit register (<= 10 Toffolis each), computed
    # and uncomputed, + S_0 on the 1000-LQ input: ~7e4 T per iterate, 0.25%. U_vac: 15 rotations per scalar site;
    # 8 k_y sectors x 30 x 30 Givens (arxiv_1711_05395); rotation-free 8-point FFFT; packet <= 59 Givens per charge.
    i = r33.intermediates
    assert i["n_fermion_modes"] == 480 == 120 * 2 * 2 and i["sector_modes"] == 60
    assert (480).bit_length() == 9 and m.COUNTER_TOF_PER_INPUT == 10
    assert i["mlae_reflection_t_per_iterate"] == (2 * 480 * 10 + 998) * 7
    assert sci(i["mlae_reflection_t_per_iterate"], 0) == ("7", 4)
    assert round(100 * i["mlae_reflection_fraction"], 2) == 0.25
    assert m.slater_rotations_ky(15, 8, 2, 2) == 8 * 30 * 30 * 2 == 14400
    assert m.scalar_product_rotations(120, 16) == 1800 and m.packet_rotations(60) == 472
    assert m.fft_transverse_t(8, 60) == 60 * (12 * 2 + 8) == 1920
    _close(i["slater_t"], 14400 * i["t_per_rotation"], 1e-12)
    _close(i["free_prep_t"], (1800 + 14400 + 472) * i["t_per_rotation"] + 1920, 1e-12)
    assert sci(i["free_prep_t"]) == ("4.6", 5) and round(100 * i["r_free"], 1) == 1.6     # '4.6e5 T, 1.6%'
    # the superseded dense bound (480 x 479 / 2 Givens) was 6.3e6 T, 16x the k_y construction
    assert sci(i["slater_dense_bound_t_r26"]) == ("6.3", 6)
    assert i["first_m_top"] == 31 and i["first_t_A"] > i["t_A"]


def test_b2_flag_contrast_and_c_prime(a, r33):
    # a = 1/2 + p0 A_R / 2; p0 = P(q = 0) of the free-vacuum incident-half charge; delta = 0 by C'
    i = r33.intermediates
    assert round(i["flag_contrast_p0"], 2) == 0.66 and round(i["query_factor_contrast"], 1) == 2.3   # '0.66', '2.3'
    assert abs(i["flag_vacuum_background"]) < 1e-12
    lo, hi = (m.flag_contrast(a, x)["p0"] for x in (0.25, 1.0))
    assert (round(lo, 2), round(hi, 2)) == (0.57, 0.77)                     # '0.57--0.77 for am_f = 0.25--1'
    # p0 >= 1 - var (q integer, mean 0): the FCS respects the Chebyshev bound
    assert i["flag_contrast_p0"] >= 1 - i["flag_contrast_var"]
    # C' = charge conjugation x flavor swap x spinor swap commutes with H on a C-violating wall; plain C does not
    assert i["c_prime_residual"] == 0.0
    h = m._wilson_h(a, 0.5, 0.5, 0.1)
    assert abs(h.conj() + h).max() > 0.1
    assert a.am_f.prov is Provenance.ASSUMED


def test_b3_preparation_sensitivity(a, r33):
    # r = 0 reproduces the printed schedule; the dressing ramp r lowers the depth: r = 1 -> depth 15, ~29 yr per
    # campaign Delta-n point instead of 6.7; r = 0.2 -> depth 27
    i = r33.intermediates
    s0 = m.prep_sensitivity(a, 0.0)
    assert s0["m_max"] == 33 == s0["m_top"] and s0["first_m_top"] == 31
    _close(s0["dn_point_yr"], i["dn_arm_yr"], 1e-12)
    _close(s0["first_dn_point_days"], i["first_dn_point_days"], 1e-12)
    _close(s0["deepest_t"], r33.hard_ops[0], 1e-12)
    assert (i["prep_r1_m_max"], i["prep_r1_m_top"]) == (16, 15)
    assert round(i["prep_r1_dn_point_yr"]) == 29
    assert i["prep_r02_m_top"] == 27
    # walls grow roughly as (1 + r)^2
    assert 3.5 < i["prep_r1_dn_wall_ratio"] < 5
    with pytest.raises(ValueError):
        m.prep_sensitivity(a, -0.5)


def test_dn_schedule_depth_never_rises_with_r(a):
    # verifier 2026-10-05: the old fixed-point loop cycled (33 -> 32 -> 33) for r near 0.007, 0.0175, 0.039, 0.054.
    # The downward scan gives a depth for every r, nonincreasing in r, and every reported circuit fits.
    prev = prev1 = 10 ** 6
    for k in range(121):
        r = 0.0005 * k
        s = m.dn_schedule(a, m.free_prep(a), r)
        s1 = m.dn_schedule(a, m.free_prep(a, int(a.label_register_first.lo)), r)
        for x in (s, s1):
            assert x["m_top"] % 2 == 1 and x["deepest"] <= 1e9
            assert x["m_top"] in (x["m_max"], x["m_max"] - 1)
        assert s["m_top"] <= prev and s1["m_top"] <= prev1
        prev, prev1 = s["m_top"], s1["m_top"]
    assert m.dn_schedule(a, m.free_prep(a), 0.0065)["m_top"] == 33
    assert [m.dn_schedule(a, m.free_prep(a), r)["m_top"] for r in (0.0072, 0.0075, 0.01)] == [31, 31, 31]


def test_flag_contrast_in_wall_background(a, r33):
    # verifier 2026-10-05: U_vac prepares the wall-background vacuum; p0 0.60-0.65 for |y| phi <= 2 m_f, so the free
    # p0^-2 = 2.3 is a lower bound (up to 2.75). wall -> 0 reproduces the free k_y-sector value.
    i = r33.intermediates
    _close(m.flag_contrast(a, wall=1e-9)["p0"], i["flag_contrast_p0"], 1e-9)
    w = i["flag_contrast_p0_wall"]
    assert i["flag_contrast_p0"] > w[0.3] > w[0.5] > w[1.0]
    assert (round(w[1.0], 2), round(w[0.3], 2)) == (0.60, 0.65)
    assert round(i["query_factor_contrast_wall_max"], 2) == 2.75
    assert abs(m.flag_contrast(a, wall=0.5, dtheta=1.0)["delta"]) < 1e-12        # C' holds on the wall too
    tex = TEX.read_text()
    assert "$0.60$--$0.65$ for $|y|\\phi\\le 2m_f$" in tex and "up to $2.75$" in tex
    assert "\\arcsin(p_0A_R)" in tex and "|p_0A_R|<\\sin(\\pi/66)" in tex
    assert "flavor-resolved charge asymmetry" in tex


def test_g4_fault_bias_model(a, r33):
    # app05:91: 0.1 expected faults -> p_f = 0.095; scale 1 - p_f; shift <= p_f / (m p0) = 4e-3; bias below the
    # first result's 3.3e-4 needs p_f <= 6.8e-3, eps_l <~ 7e-12 per T; idle qubit-cycles 9.9e10 at 10 us
    i = r33.intermediates
    _close(i["fault_prob_per_circuit"], 1 - math.exp(-0.1), 1e-12)
    _close(i["fault_bias_worst_A_R"], i["fault_prob_per_circuit"] / (33 * i["flag_contrast_p0"]), 1e-12)
    assert sci(i["fault_bias_worst_A_R"], 0) == ("4", -3)
    _close(i["fault_prob_for_bias_below_first_sigma"], 31 * i["flag_contrast_p0"] * 1e-3 / 3, 1e-12)
    assert sci(i["epsilon_l_bias_first"], 0) == ("7", -12)
    assert sci(i["epsilon_l_bias_bin"], 0) == ("2", -11)
    _close(i["idle_qubit_cycles_deepest"], 1000 * r33.hard_ops[0] * 1e-6 / 1e-5, 1e-12)
    assert sci(i["idle_qubit_cycles_deepest"]) == ("9.9", 10)
    assert round(i["idle_faults_at_box_epsilon_l"]) == 10
    assert sci(i["epsilon_l_with_idle"], 0) == ("1", -12)
    assert sci(i["epsilon_l_with_idle_first_bias"], 0) == ("7", -14)
    assert i["null_control_query_factor"] == 4.0
    # null route costed (verifier 2026-10-05): first result ~5.5 yr, ~8.7 yr with the scan, past the horizon
    assert round(i["first_result_null_yr"], 1) == 5.5 and round(i["scan_plus_first_null_yr"], 1) == 8.7
    assert i["scan_plus_first_null_yr"] > i["horizon_yr"]
    tex = TEX.read_text()
    assert "$\\lesssim 7\\times 10^{-14}$ per location" in tex and "7\\times 10^{-12}" not in tex
    assert "$\\sim 5.5$ yr, and $\\sim 8.7$ yr with the scan" in tex
    box = tex.split("2033 target: wall scattering")[1].split("end{benchmarkbox}")[0]
    assert "idle qubit-cycle" not in box and "10^{-14}" not in box


def test_b5_filtered_bath_estimate(a, r33):
    # app05:99: ~beta of patch evolution per unit Lindblad time (arxiv_2311_09207): 25 sites x 10 steps x c_T;
    # >= 2V = 240 scalar jumps per sweep: 1.3e8 T, five nucleation circuits; about seven sweeps fit 1e9
    i = r33.intermediates
    _close(i["filtered_jump_t"], 25 * 10 * i["base_budget_c_t"], 1e-12)
    assert sci(i["filtered_jump_t"]) == ("5.5", 5)
    assert i["bath_jumps_min"] == 240 and sci(i["bath_sweep_t"]) == ("1.3", 8)
    assert round(i["bath_sweep_over_nucleation_circuit"]) == 5
    assert round(i["bath_sweeps_in_envelope"]) == 7
    assert a.gibbs_cost_model.prov is Provenance.CITED and a.lr_patch_radius.prov is Provenance.ASSUMED


def test_maximum_circuit_time(r28, r33):
    # requirements table: '0.16 s (2028 shot); ~17 min (2033 deepest MLAE circuit)'
    assert round(r28.intermediates["shot_s"], 2) == 0.16
    assert round(r33.intermediates["deepest_circuit_s"] / 60) == 17


def test_referee_text_hooks():
    import pathlib
    tex = (pathlib.Path(__file__).resolve().parents[3] / "applications" / "app05_baryogenesis.tex").read_text()
    assert r"${\gtrsim}\,9.9\times 10^{8}$ T at depth $33$" in tex              # G2 box headline
    assert r"any ramp ($r>0.007$) forces a shallower schedule" in tex              # B3: bound at the priced schedule
    assert r"$p_0=0.66$" in tex and r"$p_0^{-2}=2.3$" in tex                      # B2 flag contrast
    assert "flavor-summed asymmetry vanishes" in tex                               # C'
    assert "unquantified" not in tex and "is not priced. We write $r$ for the full" not in tex
    assert r"$\gtrsim 5.5\times 10^{5}$ T" in tex and r"$\gtrsim 1.3\times 10^{8}$ T" in tex   # B5
    assert r"$\lesssim 7\times 10^{-14}$ per location" in tex                     # G4 (idles counted)
    assert "break down for the fermions" not in tex                                 # B1
    assert "where no classical method is established" not in tex                    # B1
    assert "imprint a bubble profile" not in tex                                    # B2
    assert "476 Toffolis" not in tex                                                # B2 ancilla budget
    assert r"\beta/H = -T\,d\ln(\Gamma/V)/dT" in tex                           # B4 sign
    assert r"e^{-S_3(T)/T}" in tex and r"/\hbar}" not in tex                      # B4 convention
    assert "Maximum subroutine time & N/A" not in tex                               # editorial
    assert "only known method" not in tex                                           # B1
