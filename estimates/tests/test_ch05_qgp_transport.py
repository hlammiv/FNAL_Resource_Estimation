"""Ch. 5 (QGP transport): every named number the prose states, pinned.

Rewritten 2026-09-28 for the round-A apply pass (rulings R1-R4, TRACKED_CHANGES.md):
the register is the compiled one (2T 5, Sigma(72x3) 9; Sigma(216x3) 11 at H_I since s648, 2026-10-05), the per-link-step cost is
the cited papers' magnetic + electric terms at H_KS in the papers' convention
(groups.py), eps_l is 0.1 expected faults per shot against the RFI 1e9 ceiling,
and printed numbers are exact products rounded once.

The generic contract tests in test_common.py check shape and the box totals.
These tests pin each intermediate the chapter states on the way, so that a
future edit to one number without the others fails loudly. Tolerances are
stated per assertion: `rel` is the rounding slack the chapter's own wording
implies ("~", "≈", two significant figures in boxes, one in prose).

Chapter statements the stated inputs do not reproduce are xfail(strict=True)
and filed in NEEDS_AUTHOR.md; if an author edit makes one pass, the strict
xfail fails and the item can be closed.
"""

import dataclasses
import math
from pathlib import Path

import pytest

from estimates import common as c
from estimates import ch05_qgp_transport as m
from estimates.groups import GROUPS, PRIMCOST, magnetic_per_link, electric_per_link

TEX = Path(__file__).resolve().parents[3] / "applications" / "app03_qgp_transport.tex"
NA = "NEEDS_AUTHOR.md"
EPS_28 = math.sqrt(1e-2 / 17184)          # R-TOL: 2028 shot, 17,184 synthesized rotations -> 7.628e-4
# card RF3-ch5-s216 (2026-11-12): the hop is the draft's Sigma(216x3) row, 4925 rotations per link-step (was the
# Sigma(72x3) placeholder's 2933: 4702 x 81 x 400 + mass = 1.52388e8, eps 8.101e-6, BE 2.562e-5)
N_ROT_33 = 6694 * 24 * 400 + 4 * 8 * 400   # S648: (1769 H_I gauge + 4925 hop) x 24 x 400 + mass = 6.42752e7 (3^3: 81 links, 27 sites, 2.16929e8)
EPS_33 = math.sqrt(1e-2 / N_ROT_33)       # R-TOL: 2033 evolution -> 6.790e-6 (arxiv_2611_08161 costs)
EPS_BE = math.sqrt(1e-2 / (N_ROT_33 / 10))   # R-TOL: conjectured BE, N_rot / 10 -> 2.147e-5


@pytest.fixture(scope="module")
def a():
    return m.Assumptions()


@pytest.fixture(scope="module")
def r28(a):
    return m.model(a, "2028")


@pytest.fixture(scope="module")
def r33(a):
    return m.model(a, "2033")


@pytest.fixture(scope="module")
def rcd(a):
    return m.model(a, "codesign")


def near(x, y, rel):
    assert math.isclose(x, y, rel_tol=rel), (x, y, rel)


# --- per-link-step pricing comes from groups.py (rule 3 of the apply pass) ------

def test_per_link_step_costs_are_the_groups_py_numbers(a, r28, r33):
    """Nothing is retyped: the box's gauge terms ARE groups.py's. s648: Sigma(216x3) at H_I (U_Ph twice)."""
    i28, i33 = r28.intermediates, r33.intermediates
    assert a.eps_syn.lo == c.EPS_SYN == 1e-2 and a.hamiltonian.value == "KS" and a.hamiltonian_2033.value == "I"
    e28, e33 = i28["eps_rot_2028"], i33["eps_rot_2033"]
    assert i28["magnetic_t_per_link_step_2028"] == magnetic_per_link("2T", "KS", 2, e28)
    assert i28["electric_t_per_link_step_2028"] == electric_per_link("2T", "KS", 2, e28)   # BT FFT is COMPILED
    assert i33["magnetic_t_per_link_step_2033"] == magnetic_per_link("S216x3", "I", 3, e33)
    g = GROUPS["S216x3"]
    el = 4 * g.primitives["U_FFT"].t(e33) + 2 * g.primitives["U_phi"].t(e33)
    near(i33["electric_t_per_link_step_2033_lower_bound"], el, 1e-12)
    near(i33["electric_t_per_link_step_2033_dense"], electric_per_link("S216x3", "I", 3, e33), 1e-12)
    near(i33["electric_t_per_link_step_2033_dense"], el, 1e-12)                  # no dense / fast split remains
    from estimates.groups import c_t_full
    near(i33["gauge_t_per_link_step_2033"], c_t_full("S216x3", "I", 3, e33), 1e-12)   # Table 14's closed form
    assert g.has_fft and i33["fft_compiled_2033"] is True and i28["fft_compiled_2028"] is True
    assert PRIMCOST["I"]["U_F"](3) == 4 and m.phi_mult("S216x3", "I") == 2 and m.phi_mult("S72x3", "KS") == 1


def test_papers_convention_and_report_fit_sentence(r28, r33):
    """E20 (r18): every rotation at the full fit; app03:78 '21.1 T per rotation' (2028) and '28.5 T' (2033).
    The papers' slope-only price is a legacy record (11.91 / 19.34 T); the full fit adds <~10% to the magnetic
    terms and 1.45-1.65x to the electric terms over it (model check only, no longer printed)."""
    i28, i33 = r28.intermediates, r33.intermediates
    near(i28["t_per_rotation_papers"], 11.91, 0.001)
    near(i33["t_per_rotation_papers"], 18.73, 0.001)   # 2^3 (staggered ruling 2026-11-13): was 19.74
    near(i33["t_per_rotation_full_fit"], 27.93, 0.001)       # s648 + the Sigma(216x3) hop (28.65 with the placeholder)   # 2^3 (staggered ruling 2026-11-13): was 28.94
    near(i28["t_per_rotation_report"], 21.11, 0.001)
    assert 1.0 < i28["magnetic_report_over_papers_2028"] < 1.10
    assert 1.0 < i33["magnetic_report_over_papers_2033"] < 1.10
    assert 1.45 < i28["electric_report_over_papers_2028"] < 1.65      # 1.646
    assert 1.20 < i33["electric_report_over_papers_2033"] < 1.65      # 1.217 (S648: the larger non-rotation T in U_FFT)



# --- ruling R-TOL (H. Lamm 2026-09-29): each circuit sets its own tolerance ----------

def test_rtol_rotation_counts_and_tolerances(a, r28, r33, rcd):
    """eps_rot = sqrt(1e-2 / N_rot) with N_rot this circuit's synthesized rotations per shot:
    2028 89.5/link-step (U_Tr 11/2 + 2 U_FFT x 42) x 32 x 6 = 17,184; 2033 (S648) 6695.33/link-step
    (H_I: 3 U_Tr x 15 + 4 U_FFT x 426 + 2 U_Ph x 10 = 1769 gauge = 22.5 d + 1701.5, + 4925 hop + 4/3 mass)
    x 81 x 400 = 2.16929e8; BE N/10."""
    i28, i33 = r28.intermediates, r33.intermediates
    assert i28["n_rot_per_link_step_2028"] == 89.5 and i28["n_rot_per_shot_2028"] == 17184
    near(i33["n_rot_per_link_step_2033"], 1769 + 4925 + 4 / 3, 1e-12)
    near(i33["n_rot_per_shot_2033"], N_ROT_33, 1e-12)
    near(i33["n_rot_per_shot_be"], N_ROT_33 / 10, 1e-12)
    near(i28["eps_rot_2028"], c.eps_rot_for(17184), 1e-15)
    near(i28["eps_rot_2028"], 7.628e-4, 1e-3)
    near(i33["eps_rot_2033"], 1.247e-5, 1e-3)   # 2^3 (staggered ruling 2026-11-13): was 6.790e-6
    near(i33["eps_rot_be"], 3.944e-5, 1e-3)   # 2^3 (staggered ruling 2026-11-13): was 2.147e-5
    near(i28["synthesis_error_per_shot_2028"], 1e-2, 1e-12)
    near(i33["synthesis_error_per_shot_2033"], 1e-2, 1e-12)
    assert rcd.intermediates["eps_rot_be"] == i33["eps_rot_be"]


@pytest.mark.parametrize("era,key", [("2028", "n_rot_per_shot_2028"), ("2033", "n_rot_per_shot_2033")])
def test_rtol_rotation_count_is_complete_by_perturbation(a, era, key):
    """Every eps-dependent T in the evolution is counted in N_rot: re-price the same primitives at
    two tolerances; the shot changes by exactly N_rot x 1.15 x log2(eps1/eps2)."""
    r = m.model(a, era)
    n = r.intermediates[key]
    if era == "2028":
        gk, d, links, steps, hop = "2T", 2, 32, 6.0, False
    else:
        gk, d, links, steps, hop = "S216x3", 3, 24, 400.0, True # 2^3 (staggered ruling 2026-11-13): was 3^3
    def shot(eps):
        prims, _ = m.link_step_primitives(a, gk, d, links, steps, eps)
        if hop:
            prims = prims + m.hop_primitives(a, gk, links, steps, eps, sites=8) # 2^3 (staggered ruling 2026-11-13): was 3^3
        return sum(p.count * p.t_each for p in prims)
    e1, e2 = 1e-4, 1e-6
    near(shot(e2) - shot(e1), n * 1.15 * math.log2(e1 / e2), 1e-9)
    near(n, m.rotations_per_link_step(a, gk, d, hop) * links * steps, 1e-12)

# --- 2028 box (app03:107-129) ----------------------------------------------

def test_2028_register(r28):
    """R10 (r25): 16 workspace qubits on top of 160 link + 1-2 Hadamard; box '178'."""
    i = r28.intermediates
    assert i["n_links_2028"] == 32                  # "32 links"
    assert i["n_plaq_2028"] == 16
    assert i["qubits_per_link_2028"] == GROUPS["2T"].link_qubits == 5
    assert i["gauge_qubits_2028"] == 160            # 32 x 5
    assert i["lq_2028_without_workspace"] == (161, 162)   # "+ 1--2 Hadamard ancilla" (the r23 box '162')
    assert i["workspace_2028"] == 16
    assert (i["lq_2028_lo"], i["lq_2028_hi"]) == (177, 178)   # box "178"
    assert r28.lq == (177, 178)


def test_2028_t_count_is_the_printed_product(r28):
    i = r28.intermediates
    near(i["magnetic_t_per_link_step_2028"], 1124.1, 0.001)   # "1.1e3 magnetic" (E20; 1073.5 slope-only)
    near(i["electric_t_per_link_step_2028"], 1969.2, 0.001)   # "2.0e3 electric": 2 x (98 + 42 x (1.15 log + 9.2))
    near(i["electric_over_magnetic_2028"], 1.75, 0.01)
    near(i["t_per_link_step_2028"], 3093.3, 0.001)            # "3.1e3 T/link-step"
    near(i["t_per_link_step_2028"], 3.1e3, 0.01)
    near(i["t_per_step_2028"], 9.899e4, 0.002)
    assert i["trotter_steps_2028"] == 6
    near(i["t_per_shot_2028"], 5.939e5, 0.001)                # box "5.9e5" (E20; was 4.358e5 slope-only)
    near(i["t_per_shot_2028"], 5.9e5, 0.01)
    # r25: hard_ops = (first-result deepest shot, full-benchmark deepest shot) = (9.473e4, 5.939e5)
    assert r28.hard_ops == (i["t_first_result_shot_2028"], i["t_per_shot_2028"])
    near(r28.hard_ops[0], 9.4729e4, 1e-4)                    # box "9.5e4" (one step at its own eps)
    near(r28.breakdown_total(), i["t_per_shot_2028"], 1e-9)
    assert {p.name for p in r28.breakdown} == {"U_inv_2T", "U_mul_2T", "U_Tr_2T", "U_FFT_2T"}
    assert all(p.status is c.CircuitStatus.COMPILED for p in r28.breakdown)


def test_2028_over_the_rfi_envelope(r28):
    """app03:78,120: '5.9x the 1e5-T first-generation envelope'; app03:74: '1e5 T buys a single step'."""
    i = r28.intermediates
    near(i["t_per_shot_2028_over_rfi"], 5.9, 0.01)
    assert 1.0 < i["steps_within_rfi_2028"] < 2.0
    assert i["t_max_within_rfi_fm"] == pytest.approx(0.025)     # one step of a_t/4 (step rule, 2026-10-04)
    # ruling (b) 2026-10-01: the benchmark needs >= 5.9x (E20); the levers' 3-10x closes it only at the upper part
    near(i["per_step_cut_needed_2028"], 5.939, 0.002)
    assert i["levers_close_2028_lo"] is False and i["levers_close_2028_hi"] is True
    assert i["t_max_within_rfi_at_dt_over_10_fm"] == pytest.approx(0.01)
    # the retired working figure (1e3 T/plaquette) gave the old box's 9.6e4
    assert i["legacy_t_per_shot_2028_at_1e3_per_plaq"] == 9.6e4


def test_2028_trotter_conventions(r28):
    i = r28.intermediates
    assert i["dt_2028_fm"] == pytest.approx(0.025)         # Delta t = a_t/4 (step rule, 2026-10-04; was a_t)
    assert i["t_max_2028_fm"] == pytest.approx(0.15)       # "t_max = 0.15 fm/c" (was 0.6)
    assert i["dt_first_2028_fm"] == pytest.approx(0.045)   # "one step of 0.045 fm/c"
    assert i["xi_anisotropy"] == pytest.approx(2.0)        # "xi = a_s/a_t = 2"


def test_2028_temperature_and_momentum(r28):
    i = r28.intermediates
    near(i["Tc_2028_MeV"], 492.8, 1e-6)             # 1.12 x 440
    near(i["T_2028_MeV"], 740, 0.01)                # "T = 1.5 T_c ~ 740 MeV"
    near(i["L_perp_fm"], 0.8, 1e-9)                 # "L_perp = 0.8 fm"
    near(i["k_2028_GeV"], 1.55, 0.01)               # "k = 2 pi / L_perp ~ 1.55 GeV"
    near(i["k_over_T_2028"], 2.1, 0.02)             # "~ 2.1 T"


def test_2028_shots_and_wall_time(r28):
    """r25 (R1, R4): timeslices t = 0..6 a_t (seven), 1e3 each, every shot at its own depth and tolerance;
    first result t = 0 and one step. One machine, 1 us per T + 0.1 ms per shot."""
    i = r28.intermediates
    assert i["timeslices_2028"] == 7 and i["shots_2028"] == 7e3          # box "7e3 (1e3/timeslice x 7)"
    assert r28.shots == (7e3, 7e3)
    assert i["shots_first_result_2028"] == 2e3                            # box "2e3"
    assert i["gate_time_from_factories_s"] == pytest.approx(1e-6)         # 10 us / 10 factories
    near(i["shot_time_2028_s"], 0.594, 0.002)       # requirements "~0.6 s at the 2028 benchmark" (deepest shot)
    t = i["t_per_slice_2028"]
    assert t[0] == 0.0 and t[6] == i["t_per_shot_2028"]
    assert all(t[j] < t[j + 1] for j in range(6))
    for j in range(1, 7):                                                 # each slice is its own R-TOL circuit
        near(i["eps_per_slice_2028"][j], c.eps_rot_for(17184 * j / 6), 1e-12)
    near(i["t_first_result_shot_2028"], 9.473e4, 1e-3)
    assert i["t_first_result_shot_2028_over_rfi"] < 1.0                    # "fits the 1e5-T envelope as it stands"
    near(i["wall_first_result_2028_s"], 1e3 * (t[1] * 1e-6 + 1e-4) + 1e3 * 1e-4, 1e-12)
    near(i["wall_first_result_2028_min"], 1.58, 0.005)                    # "1.6 min"
    near(r28.wall_time_s[0], sum(1e3 * (x * 1e-6 + 1e-4) for x in t), 1e-12)   # F* >= 10: no depth charge
    near(i["wall_time_2028_min"], 34.32, 0.001)                           # box "34 min"
    near(i["wall_time_2028_min_r23"], 99.0, 0.001)                        # the r23 box (1e4 x 6 steps), record
    assert i["shot_overhead_s"] == 1e-4


def test_2028_workspace_restores_factory_feed(a, r28):
    """R10: with 16 workspace qubits 8 U_FFT run at once and F* = N_T / D_T >= 10 at the high schedule; with the
    2 workspace qubits needed to run any primitive F* ~ 1.7 (factory.json: 1.5). At 64 workspace the schedule is
    factory.json's 628 / 3882 per step."""
    i = r28.intermediates
    assert 10 <= i["f_star_2028_per_step"][0] < 11 and i["f_star_2028_per_step"][1] > 50
    assert i["f_star_2028_at_2_workspace"] < 2
    b = dataclasses.replace(a, workspace_2028=c.Assumed(64, "test"))
    lo, hi = m.step_depth_2028(b, 21.1098)
    near(lo, 628, 0.002)
    near(hi, 3882, 0.001)
    b = dataclasses.replace(a, workspace_2028=c.Assumed(8, "test"))
    assert m.model(b, "2028").intermediates["f_star"][0] < 10           # 8 is too few
    # the R9 exports: full benchmark, shot-averaged
    e = c.depth_exports(i["t_per_shot"], i["t_depth_per_shot"], 7e3, 1e-6, 1e-4)
    for k in c.DEPTH_EXPORTS[:5]:
        assert i[k] == e[k]
    assert i["baseline_ok"] == (True, True)
    assert i["wall_first_result_s"] == i["wall_first_result_2028_s"]
    assert i["wall_campaign_s"] == r28.wall_time_s[0]
    near(i["wall_serial_s"][0], i["wall_campaign_s"], 1e-12)


# --- 2033 box (app03:131-171) ----------------------------------------------

def test_2033_register_from_eq_Nq(r33):
    i = r33.intermediates
    g = GROUPS["S216x3"]
    assert i["qubits_per_link_2033"] == g.link_qubits == 11      # s648: the paper's qubit register
    assert i["ceil_log2_order_2033"] == 10                       # "a dense ceil(log2 648) = 10 encoding would save 81 LQ"
    assert i["n_links_2033"] == 24 and i["n_plaq_2033"] == 24 # 2^3 (staggered ruling 2026-11-13): was 3^3
    assert i["gauge_qubits_2033"] == 264            # "264 gauge (3 . 8 . 11)", was 891 # 2^3 (staggered ruling 2026-11-13): was 3^3
    assert i["fermion_qubits_2033"] == 72          # "72 fermion (3 . 3 . 8)", was 243 # 2^3 (staggered ruling 2026-11-13): was 3^3
    assert (i["lq_2033_lo"], i["lq_2033_hi"]) == (436, 536)   # "+ 100--200 ancilla" -> box "~490 (436-536)"; was (1234, 1334) on 3^3
    assert r33.lq == (436, 536)   # 2^3; was (1234, 1334)
    assert i["lq_2033_mid"] == 486   # 2^3; was 1284
    near(i["lq_2033_lo_over_1000_limit"], 0.436, 1e-9)   # "overfill it by a quarter to a third"; box "1.2-1.3x"   # 2^3 (staggered ruling 2026-11-13): was 1.234
    near(i["log10_hilbert_dim_2033"], 24 * math.log10(648) + 72 * math.log10(2), 1e-12)  # 648^24 2^72 # 2^3 (staggered ruling 2026-11-13): was 3^3
    assert round(i["log10_hilbert_dim_2033"]) == 89                                       # "~10^89" (2^3; was 301)


def test_2033_evolution_t_count(r33):
    i = r33.intermediates
    assert i["dt_2033_fm"] == pytest.approx(0.01)   # Delta t = a_t/10
    assert i["trotter_steps_2033"] == pytest.approx(400)   # "4 x 10^2 2nd-order Trotter steps"
    near(i["magnetic_t_per_link_step_2033"], 2.8590e5, 0.001)   # "2.9e5": 56 U_mul + 24 U_inv + 3 U_Tr (H_I)
    near(i["electric_t_per_link_step_2033_lower_bound"], 87347, 0.001)   # "8.9e4": 4 U_FFT + 2 U_Ph   # 2^3 (staggered ruling 2026-11-13): was 8.9087e4
    near(i["gauge_t_per_link_step_2033"], 3.7321e5, 0.001)   # "3.7e5 T in all"   # 2^3 (staggered ruling 2026-11-13): was 3.7500e5
    # card RF3-ch5-s216 (2026-11-12): the draft's Sigma(216x3) hop row (was the Sigma(72x3) placeholder: hop 1.10095e5,
    # link-step 4.8463e5, step 3.9255e7, hop 22.7%, reference shot 1.5702e10, share='draft' 1.8116e10, BE 1.5411e9)
    near(i["hop_t_per_link_step_2033"], 278845, 1e-4)     # Sigma(216x3) draft counts: "2.8e5 T"   # 2^3 (staggered ruling 2026-11-13): was 2.83815e5
    near(i["mass_t_per_site_step_2033"], 160.74, 1e-4)       # one HWP(9) -> "160 T per site-step"   # 2^3 (staggered ruling 2026-11-13): was 164.77
    near(i["t_per_link_step_2033"], 6.5211e5, 0.001)         # "6.6e5 per link-step"   # 2^3 (staggered ruling 2026-11-13): was 6.5886e5
    near(i["t_per_step_2033"], 1.5651e7, 0.001)              # "5.3e7 T"   # 2^3 (staggered ruling 2026-11-13): was 5.3368e7
    near(i["hop_share_of_step_2033"], 0.4276, 0.002)          # "43% of the step"   # 2^3 (staggered ruling 2026-11-13): was 0.4308
    assert i["mass_share_of_step_2033"] < 1e-3               # "under 0.1% of the step"
    near(i["t_per_shot_2033_mu0_evolution"], 6.2602e9, 0.001)   # reference shot to 4 fm/c "2.1e10"   # 2^3 (staggered ruling 2026-11-13): was 2.1347e10
    near(i["t_per_shot_2033_mu0_evolution_share_draft"], 1.0109e10, 0.001)   # gap bullet "3.4e10"   # 2^3 (staggered ruling 2026-11-13): was 3.4336e10
    near(i["t_per_shot_2033_be"], 6.1375e8, 0.001)  # route (a) "~2.1e9 T/shot" with BE, at the BE's own eps   # 2^3 (staggered ruling 2026-11-13): was 2.0933e9
    # quench ruling (2026-10-05): the box prints the deepest shot at a_s = 0.2 fm, one circuit of quench prep + fit
    # evolution, at the corners (T = 300 MeV, 1/k, c = 2) and (225 MeV, 1/(2 pi T), c = 2 pi) (temperature ruling
    # 2026-10-05: T = 1.5 T_c^lat, T_c^lat = 150-200 MeV; was 1.05086e10, 3.14203e10 at T = 160 MeV)
    assert r33.hard_ops == (pytest.approx(2.54721e9, rel=1e-3), pytest.approx(9.40696e9, rel=1e-3))   # 2^3; were 9.06056e9, 3.23463e10 on 3^3;   # were 6.66665e9, 2.37886e10 (placeholder hop)
    # the retired working figure (5e3 T/plaquette) gave the old box's 1.62e8
    assert i["legacy_t_per_shot_2033_at_5e3_per_plaq"] == pytest.approx(4.8e7)   # 2^3; was 1.62e8 on 3^3


def test_2033_hop_is_the_fermion_primitives_draft(a, r33, r28):
    """r17 ruling (e), H. Lamm 2026-10-01: the 2033 hop is groups.hop_link_cost from the unpublished
    Fermion_Primitives gate counts with the estimated pieces (undo, SU(3) squish uncompute, MBU, HWP phasing),
    rotations at the full fit; the mass is one HWP(9) per site. 2028 is pure gauge and has neither.
    E21 (1) (2026-10-01, "make the switch"): squish + parity held once per link (share='link'), undo kept;
    share='draft' is the conservative sensitivity. Card RF3-ch5-s216 (2026-11-12): the hop is the draft's own
    Sigma(216x3) row (groups.FERMION_DRAFT['S216x3']); the Sigma(72x3) placeholder (3672, 360, 2933) is retired."""
    from estimates.groups import FERMION_DRAFT, FP_SQUISH_ANCILLA, hop_link_cost, staggered_mass_site, rhop_link
    assert a.hop_group_2033.value == "S216x3" and a.hop_group_2033.prov is c.Provenance.CITED
    assert a.hop_group_2033.src == "FermionPrimitives_unpub" and "RF3-ch5-s216" in a.hop_group_2033.note
    assert m.hop_group(a, "S216x3") == "S216x3" and m.hop_group(a, "S36x3") == "S36x3" and m.hop_group(a, "S72x3") == "S72x3"
    assert "S216x3" in FERMION_DRAFT and a.group_2033.value == "S216x3"   # the hop group IS the gauge group: no stand-in
    i = r33.intermediates
    near(i["eps_rot_2033"], EPS_33, 1e-12)
    h = hop_link_cost("S216x3", 3, eps=EPS_33, synthesis="rus")
    assert h["t"] > 2.5 * hop_link_cost("S72x3", 3, eps=EPS_33, synthesis="rus")["t"]   # 2.57x the retired placeholder
    assert a.fermion_synthesis.value == "rus"
    near(i["hop_t_per_link_step_2033"], h["t"], 1e-12)
    assert (i["hop_toffoli_per_link_step_2033"], i["hop_t_direct_per_link_step_2033"],
            i["hop_n_rot_per_link_step_2033"]) == (20068, 792, 4925)          # 2 x 7881 + 1402 + 2 x (1224 + 220) + 16; 2 x 396; 12 x 278 + 4 x 396 + 5
    t_rot = 1.15 * math.log2(1 / EPS_33) + 9.2
    near(i["hop_t_per_link_step_2033"], 7 * 20068 + 792 + 4925 * t_rot, 1e-12)
    near(i["hop_undo_t_per_link_step_2033"], 1668 * t_rot, 1e-9)             # the undo is 2 x 3 x 278 = 1,668 rotations
    near(i["hop_undo_t_per_link_step_2033"], 46595, 0.002)                 # "inverse color basis change (4.8e4 T per link-step)"   # 2^3 (staggered ruling 2026-11-13): was 4.8278e4
    near(i["hop_t_per_link_step_2033_share_draft"], 7 * 77336 + 792 + 4925 * t_rot, 1e-12)   # r17 headline
    near(i["hop_t_per_link_step_2033_share_draft"], 6.7972e5, 0.001)   # 2^3 (staggered ruling 2026-11-13): was 6.8469e5
    near(i["hop_t_per_link_step_2033_no_undo"], 2.3225e5, 0.002)   # 2^3 (staggered ruling 2026-11-13): was 2.3554e5
    assert hop_link_cost("S216x3", 3, eps=EPS_33, share="link")["t"] == h["t"]   # the default IS share='link'
    near(i["hop_over_rhop_2033"], 17.95, 0.005)                               # model-only comparison (sentence dropped from the chapter, r17 fix)   # 2^3 (staggered ruling 2026-11-13): was 18.26
    near(i["rhop_t_per_link_step_2033_retired"], rhop_link("S216x3", 3, 3, EPS_33, synthesis="rus-slope")["t"], 1e-12)
    # the draft's squish scratch (su3_diag.tex:105) is recorded and sits inside the box's ancilla allowance; LQ unchanged
    assert i["hop_squish_ancilla_2033"] == FP_SQUISH_ANCILLA["S216x3"] == 29
    assert 29 + c.hwp_ancilla(18) < a.anc_2033.lo and r33.lq == (436, 536)   # 2^3; was (1234, 1334)
    ms = staggered_mass_site(3, 3, eps=EPS_33)
    near(i["mass_t_per_site_step_2033"], ms["t"], 1e-12)
    near(i["mass_t_per_site_step_2033"], 4 * t_rot + 49, 1e-12)
    # E26 (r22): a phasing group of k holds k - w(k) ancilla (was 21 and 11 with the weight register added on top)
    assert i["hop_phasing_hwp_ancilla_2033"] == c.hwp_ancilla(18) == 16 < a.anc_2033.lo   # LQ unchanged
    assert i["mass_hwp_ancilla_2033"] == c.hwp_ancilla(9) == 7
    near(i["t_per_link_step_2033"], i["gauge_t_per_link_step_2033"] + i["hop_t_per_link_step_2033"]
         + i["mass_t_per_site_step_2033"] / 3, 1e-12)
    names = {p.name for p in r33.breakdown}
    assert {"hop_colour_squish_and_parity_S216x3", "hop_colour_rotations_S216x3", "hop_hop_squish_and_flags_S216x3",
            "hop_hop_diagonalizers_S216x3", "hop_hop_phasing_hwp_S216x3", "mass_hwp_S216x3"} <= names
    assert all("unpublished gate counts" in p.note for p in r33.breakdown if p.name.startswith("hop_"))
    near(r33.breakdown_total(), r33.hard_ops[1], 1e-12)
    assert not any(p.name.startswith(("hop_", "mass_")) for p in r28.breakdown)
    near(i["t_per_link_step_2033_at_papers_fiducial_eps"], 9.3101e5, 0.001)   # record (eps_T=1e-8 split over N_rot)   # 2^3 (staggered ruling 2026-11-13): was 9.4452e5
    near(i["t_per_shot_2033_mu0_evolution_at_papers_fiducial_eps"], 8.9377e9, 0.001)   # 2^3 (staggered ruling 2026-11-13): was 3.0602e10
    with pytest.raises(ValueError):
        dataclasses.replace(a, fermion_synthesis=c.Cited("rus-slope", "x"))   # never slope-only for new prices


def test_2033_electric_is_the_compiled_s216_fft(r33):
    """s648: the Sigma(216x3) paper constructs the fast transform (COMPILED in groups.py), so there is no
    floor / dense split; the old Sigma(72x3) dense fallback (~3.0e11 T/shot) is retired."""
    i = r33.intermediates
    near(i["electric_dense_over_magnetic_2033"], i["electric_over_magnetic_2033"], 1e-12)
    near(i["t_per_shot_2033_mu0_evolution_dense_fourier"], i["t_per_shot_2033_mu0_evolution"], 1e-12)
    rows = {p.name: p for p in r33.breakdown}
    for name in ("U_inv_S216x3", "U_mul_S216x3", "U_Tr_S216x3", "U_FFT_S216x3", "U_phi_S216x3"):
        assert rows[name].status is c.CircuitStatus.COMPILED
    assert "COMPILED fast transform" in rows["U_FFT_S216x3"].note


def test_2033_breakdown_is_one_quench_plus_fit_shot(a, r33):
    """Contract rule 2: the breakdown sums to a per-shot total: ONE shot at a_s = 0.2 fm at the slow corner (tau =
    1/(2 pi T), c = 2 pi, T = 225 MeV), 572 quench steps + 32 fit steps of the same primitives (quench ruling 2026-10-05; the Gibbs
    and E-rho-OQ preparation rows are retired). s648: H_I multiplicities, 56 U_mul, 24 U_inv, 3 U_Tr, 4 U_FFT, 2 U_Ph
    per link-step at d = 3."""
    i = r33.intermediates
    rows = {p.name: p for p in r33.breakdown}
    hop_rows = {"hop_colour_squish_and_parity_S216x3", "hop_colour_rotations_S216x3", "hop_hop_squish_and_flags_S216x3",
                "hop_hop_diagonalizers_S216x3", "hop_hop_phasing_hwp_S216x3"}
    assert set(rows) == {"U_inv_S216x3", "U_mul_S216x3", "U_Tr_S216x3", "U_FFT_S216x3", "U_phi_S216x3",
                         "mass_hwp_S216x3"} | hop_rows
    steps = i["quench_prep_steps"][1] + i["fit_steps"][1][1]
    assert steps == 572 + 27 == 599   # 2^3; was 572 + 32 = 604 on 3^3
    for n in hop_rows:
        assert rows[n].count == pytest.approx(24 * steps) # 2^3 (staggered ruling 2026-11-13): was 3^3
    assert rows["mass_hwp_S216x3"].count == pytest.approx(8 * steps) # 2^3 (staggered ruling 2026-11-13): was 3^3
    near(r33.breakdown_total(), i["t_deepest"][1], 1e-9)
    near(r33.breakdown_total(), r33.hard_ops[1], 1e-9)
    for name, mult in (("U_mul_S216x3", 56), ("U_inv_S216x3", 24), ("U_Tr_S216x3", 3), ("U_FFT_S216x3", 4),
                       ("U_phi_S216x3", 2)):
        assert rows[name].count == pytest.approx(mult * 24 * steps) # 2^3 (staggered ruling 2026-11-13): was 3^3
    assert rows["U_mul_S216x3"].t_each == GROUPS["S216x3"].primitives["U_mul"].t(EPS_33) == 2548   # Toffoli-only
    eps604 = c.eps_rot_for(N_ROT_33 * 599 / 400)   # 2^3: 599 steps (was 604)
                                  # its own R-TOL tolerance (604 steps)
    near(rows["U_Tr_S216x3"].t_each, GROUPS["S216x3"].primitives["U_Tr"].t(eps604), 1e-12)
    assert not any("gibbs" in n or "eroq" in n for n in rows)


def test_2033_momentum_and_decay_times(r33):
    i = r33.intermediates
    near(i["L_2033_fm"], 0.4, 1e-9)                 # "L = 0.4 fm" (2^3; was 0.6)
    near(i["k_min_2033_GeV"], 3.1, 0.02)            # "k_min = 2 pi / L ~ 3.1 GeV" (2^3; was 2.1)
    # temperature ruling (2026-10-05): T = 1.5 T_c^lat, T_c^lat = 150-200 MeV -> T = 225-300 MeV; pairs are (T lo, T hi)
    assert i["T_c_lat_2033_MeV"] == (150, 200) and i["T_2033_MeV"] == (225, 300)
    near(i["T_a_s_2033"][0], 0.228, 0.002) and near(i["T_a_s_2033"][1], 0.304, 0.002)
    near(i["LT_2033"][0], 0.456, 0.002)                 # "LT ~ 0.5-0.6" (2^3; was 0.684, "0.7-0.9"), below 1: a small hot box, not a medium
    near(i["LT_2033"][1], 0.608, 0.002)   # was 0.912
    assert max(i["LT_2033"]) < 1.0
    near(i["k_min_over_T_2033"][0], 13.78, 0.002)        # "k_min ~ 10-14 T" (2^3; was 9.18, "7-9 T") (was ~13 T at 160 MeV)
    near(i["k_min_over_T_2033"][1], 10.33, 0.002)   # was 6.89
    # the hydrodynamic tau = T/((eta/s) k^2) at eta/s = 0.15 is a record the chapter says has no basis at this k
    assert i["eta_s_assumed"] == 0.15
    near(i["tau_decay_at_eta_s_assumed_fm"][0], 0.0308, 0.002)   # "~0.07-0.09 fm/c" (was 0.049)   # 2^3 (staggered ruling 2026-11-13): was 0.0693
    near(i["tau_decay_at_eta_s_assumed_fm"][1], 0.0411, 0.002)   # was 0.0924
    near(i["n_decay_times_at_eta_s_assumed"][0], 129.8, 0.002)   # was 57.7
    near(i["n_decay_times_at_eta_s_assumed"][1], 97.4, 0.002)   # was 43.3


def test_temperature_is_1p5_tc_lat(a):
    """Author ruling 2026-10-05: the 2033 target is T = 1.5 T_c^lat (the 2028 convention), T_c^lat the crossover of the
    simulated lattice theory, unknown, bracketed 150-200 MeV from the HotQCD chiral-limit and physical crossovers and
    measured in situ; no T_2033_MeV input remains."""
    assert a.eta_s_assumed.prov is c.Provenance.STATED and a.eta_s_assumed.lo == 0.15
    assert not hasattr(a, "T_2033_MeV") and not hasattr(a, "tau_decay_fm")
    assert a.Tc_lat_2033_MeV.prov is c.Provenance.ASSUMED and (a.Tc_lat_2033_MeV.lo, a.Tc_lat_2033_MeV.hi) == (150, 200)
    assert "arxiv_1903_04801" in a.Tc_lat_2033_MeV.src and "arxiv_1812_08235" in a.Tc_lat_2033_MeV.src
    assert a.T_over_Tc_2033.prov is c.Provenance.STATED and a.T_over_Tc_2033.lo == 1.5 == a.T_over_Tc_2028.lo
    assert a.T_2033() == (225, 300)
    for bad in (dict(Tc_lat_2033_MeV=c.Assumed((200.0, 150.0), "x")), dict(T_over_Tc_2033=c.Stated(0.0, "x")),
                dict(tc_diag_points=c.Assumed(-1, "x"))):
        with pytest.raises(ValueError):
            dataclasses.replace(a, **bad)
    tex = TEX.read_text()
    assert r"$T=1.5\,T_c^{\rm lat}$" in tex and r"$T_c^{\rm lat}$" in tex
    assert r"\cite{arxiv_1903_04801}" in tex and r"\cite{arxiv_1812_08235}" in tex
    assert "160" not in tex.replace("$160$ T per site-step", "")      # no physical-QCD 160 MeV left in the chapter
    assert r"$LT\approx 0.5$--$0.6$" in tex   # 2^3; was 0.7--0.9
    assert True and "With $LT<1$ the box is smaller than the thermal wavelength, so it is not a medium" in tex
    # re-pinned 2026-11-09 (referee cards, Q5-6): the decay-time sentence became Fig. (eta/s = 0.15 is in the plot label); T stays in the prose
    assert r"The counts are for $T=225$--$300$ MeV ($1.5\,T_c^{\rm lat}$ over that range)" in tex


def test_2033_volume_requirements(r33):
    i = r33.intermediates
    assert 0.6 <= i["m_pi_L_2033_lo"] <= 0.61 and 0.81 <= i["m_pi_L_2033_hi"] <= 0.82   # "m_pi L ~ 0.6-0.8" (2^3; was 0.9-1.3, "~1")
    near(i["L_for_mpiL4_fm"], 2.3, 0.03)            # "L >~ 2.3 fm" (at m_pi = 350 MeV)
    near(i["L_ext_sites"], 12, 0.07)                # "~12^3 lattice"
    assert i["lq_12cubed_lo"] == 3 * 1728 * 11 + 9 * 1728 + 100 == 72676  # Eq. (Nq) at 11 qubits/link (s648)
    near(i["lq_12cubed_lo"], 7e4, 0.04)             # "~7e4 LQ"
    near(i["factor_beyond_2033"], 70, 0.04)         # "factor of ~70" beyond the 1000-LQ limit
    near(i["lq_12cubed_over_lanl"], 2.0, 0.02)      # "about twice the 36,000-LQ LANL estimate"
    assert (i["lq_4cubed_lo"], i["lq_4cubed_hi"]) == (2788, 2888)   # "4^3, 2788--2888 LQ"
    near(i["k_min_over_T_at_L_ext_fm"][0], 2.40, 0.005)   # "k_min ~ 1.8-2.4 T" at L ~ 2.3 fm (T = 225-300 MeV)
    near(i["k_min_over_T_at_L_ext_fm"][1], 1.80, 0.005)
    near(i["k_min_over_T_at_12cubed"][0], 2.30, 0.005)    # "(1.7-2.3 T on the 12^3 lattice)"
    near(i["k_min_over_T_at_12cubed"][1], 1.72, 0.005)
    near(i["L_for_kmin_eq_T_fm"][0], 5.51, 0.005)         # "requires L >~ 4-6 fm" (was 7.75 at 160 MeV)
    near(i["L_for_kmin_eq_T_fm"][1], 4.13, 0.005)
    assert i["sigma360_qubits_per_link_dense"] == 11   # "ceil(log2 1080) = 11 dense"
    assert i["sigma360_qubits_per_link"] == 12      # "12 qubits/link in its compiled encoding"
    assert i["sigma360_extra_lq"] == 24             # "+24 LQ on the 2^3 lattice": 24 x (12 - 11); was 81 on 3^3


def test_2033_eps_l_requirement_r3(r33):
    """R3: RFI 1e9 is the ceiling; eps_l is a requirement at 0.1 expected faults per shot. These keys are the 4 fm/c
    reference shot (the r23 box), kept as model records; r25 numbers are in test_2033_own_depth_*. s648 values."""
    i = r33.intermediates
    assert i["ceiling_2033_t"] == 1e9
    near(i["ratio_reference_to_ceiling"], 6.26, 0.001)   # 2^3 (staggered ruling 2026-11-13): was 21.35
    near(i["eps_l_required_reference"], 1.597e-11, 0.001)   # 2^3 (staggered ruling 2026-11-13): was 4.684e-12
    # the reference-shot (4 fm/c) keys above are records; the box's eps_l are the deepest quench + fit shots
    assert r33.epsilon_l == (i["eps_l_deepest"][1], i["eps_l_deepest"][0])
    near(0.1 / i["t_per_shot_2033_be"], 1.629e-10, 0.001)   # 2^3 (staggered ruling 2026-11-13): was 4.777e-11


def test_2033_shots_campaign_and_wall_time_reference_records(r33):
    """The r23 pricing (every shot to 4 fm/c, 1e3 shots per grid point, 96 points), kept as model records only; the
    chapter no longer prints it (r25, R4). s648 values."""
    i = r33.intermediates
    assert i["n_grid_2033"] == 96 and i["shots_2033"] == pytest.approx(9.6e4)
    assert not any("machine" in k for k in i)
    near(i["n_subroutine_calls"], 1.92e5, 1e-9)


def test_2033_own_depth_design_and_shots(a, r33):
    """r25 (R1, R4): two fit points at Gamma t = 0.5 and 1.8 (20/80), rounded up onto the 0.01 fm/c grid, over the
    decay band 1/k = 0.0955 .. 1/(2 pi T) = 0.140 fm/c at T = 225 MeV (temperature ruling 2026-10-05; was 0.196 at
    160 MeV); shots from the Fisher information (c0, Gamma free, Cbar = 0.16): 30% first result, 20% campaign.
    s648: the H_I step rule binds both deepest shots."""
    i = r33.intermediates
    assert i["T_corners_MeV"] == (300, 225)                                   # cheap corner at T hi, expensive at T lo
    near(i["tau_band_fm"][0], 197.327 / (2 * math.pi * 197.327 / 0.4), 1e-12)   # T1: fast end 1/k = 0.0637 fm/c (2^3; was 0.0955 on 3^3)
    near(i["tau_band_fm"][0], 0.0637, 1e-3)                                    # "1/k ~ 0.06 fm/c" (was 0.0955)
    near(i["tau_band_fm"][1], 197.327 / (2 * math.pi * 225), 1e-12)          # "1/(2 pi T) ~ 0.14 fm/c"
    assert i["fit_steps_grid"] == ((4, 12), (7, 26))                          # the a_t/10 grid (2^3; was ((5, 18), (7, 26)))
    assert i["fit_steps"] == ((4, 12), (7, 27))                               # step rule (H_I): 26 -> 27 (2^3; was ((5, 19), (7, 32)))
    near(i["fit_t_fm"][1][1], 0.26, 1e-9)                                      # "0.18-0.26 fm/c"
    near(i["fit_t_fm"][0][1], 0.12, 1e-9)   # 2^3 (staggered ruling 2026-11-13): was 0.18
    assert all(x >= 0.5 for pair in i["fit_x_realized"] for x in pair)        # the transient is excluded
    near(i["shots_per_point_campaign_design_nominal"], 3.42e4, 0.001)         # shot audit design B
    near(i["shots_per_point_campaign"][0], 4.43e4, 0.001)                     # "3.4-3.6e4 per point" (a_s = 0.2 fm)   # 2^3 (staggered ruling 2026-11-13): was 3.599e4
    near(i["shots_per_point_campaign"][1], 3.443e4, 0.001)
    near(i["shots_first_result"][0], 1.969e4, 0.001)                           # "1.5-2.0e4" (2^3; was 1.600e4, "1.5-1.6e4")
    near(i["shots_first_result"][1], 1.530e4, 0.001)
    for k in (0, 1):                                                          # 30% vs 20%: Fisher scales as 1/target^2
        near(i["shots_first_result"][k] / i["shots_per_point_campaign"][k], (0.2 / 0.3) ** 2, 1e-12)
    near(i["design_A_over_B"], 0.366, 0.002)                                   # "0.37x"
    near(i["design_C_over_B"], 2.72, 0.002)                                    # "2.7x"
    assert (i["n_points_campaign"], i["n_points_mu0"], i["n_points_mupos"]) == (16, 4, 12)
    near(i["shots_campaign"][0], 1.4175e6, 0.001)                              # "1.1e6" (both spacings)   # 2^3 (staggered ruling 2026-11-13): was 1.1518e6
    near(i["shots_campaign"][1], 1.1406e6, 0.001)
    assert r33.shots == (i["shots_campaign"][1], i["shots_campaign"][0])
    # each shot is its own R-TOL circuit: quench prep + fit steps (quench ruling 2026-10-05)
    for q, (n1, n2), (e1, e2) in zip(i["quench_prep_steps"], i["fit_steps"], i["eps_rot_fit_points"]):
        near(e1, c.eps_rot_for(N_ROT_33 * (q + n1) / 400), 1e-12)
        near(e2, c.eps_rot_for(N_ROT_33 * (q + n2) / 400), 1e-12)
    near(i["eps_rot_fit_points"][0][0], 1.997e-5, 0.01)                        # 152 + 5 steps   # 2^3 (staggered ruling 2026-11-13): was 1.084e-5
    near(i["eps_rot_fit_points"][1][1], 1.019e-5, 0.01)                        # 572 + 32 steps   # 2^3 (staggered ruling 2026-11-13): was 5.525e-6
    near(i["t_rot_fit_points"][0][0], 27.15, 0.002)                            # 28.2-29.3 T per rotation   # 2^3 (staggered ruling 2026-11-13): was 28.17
    near(i["t_rot_fit_points"][1][1], 28.27, 0.002)   # 2^3 (staggered ruling 2026-11-13): was 29.29
    o = m.own_depth_2033(m.Assumptions())
    p2 = o["second"]["r"][1]["pts"][-1]                                        # the second spacing's late shot
    assert p2["steps"] == 40 and p2["prep_steps"] == 755   # 2^3; was 48
    near(p2["eps"], c.eps_rot_for(N_ROT_33 * 795 / 400), 1e-12)   # 2^3: 755 + 40 (was 803)


def test_2033_own_depth_per_shot_and_eps_l(r33):
    """Box per-shot row and eps_l (R3 at 0.1 faults per shot) with the quench preparation in every shot (author ruling
    2026-10-05): ramp 20 steps + c / T at the a_t/10 grid step, c = 2 at T = 300 MeV (132 steps) to 2 pi at 225 MeV
    (552), then the fit evolution
    (19 / 52 steps), one circuit at its own eps. Shots unchanged."""
    i = r33.intermediates
    near(i["quench_inv_T_steps"][0], 65.776, 1e-4)                 # 1/T = 0.658 fm/c at 300 MeV, dt = 0.01 fm/c
    near(i["quench_inv_T_steps"][1], 87.701, 1e-4)                 # 0.877 fm/c at 225 MeV (was 123.3 at 160 MeV)
    assert i["quench_ramp_steps"] == 20                            # one a_s = 0.2 fm/c
    assert i["quench_th_steps"] == (132, 552)                      # ceil(2 x 65.78), ceil(2 pi x 87.70)
    assert i["quench_prep_steps"] == (152, 572)                    # box "152-572 steps"
    near(i["quench_t_th_fm"][0], 1.3155, 1e-3)
    near(i["quench_t_th_fm"][1], 5.5104, 1e-3)
    # card RF3-ch5-s216 (2026-11-12), Sigma(216x3) hop row; with the placeholder: 5.9259e9, 2.2528e10; 7.4074e8,
    # 1.2603e9; 6.66665e9, 2.37886e10 ("6.7e9-2.4e10")
    near(i["quench_prep_t"][0], 2.3608e9, 1e-3)                    # "8.1e9-3.1e10 quench"   # 2^3 (staggered ruling 2026-11-13): was 8.0538e9
    near(i["quench_prep_t"][1], 8.9829e9, 1e-3)   # 2^3 (staggered ruling 2026-11-13): was 3.0633e10
    near(i["t_fit_deepest"][0], 1.86381e8, 1e-3)                   # "1.0e9-1.7e9 fit evolution"   # 2^3 (staggered ruling 2026-11-13): was 1.00673e9
    near(i["t_fit_deepest"][1], 4.2402e8, 1e-3)   # 2^3 (staggered ruling 2026-11-13): was 1.71371e9
    near(i["t_deepest"][0], 2.54721e9, 1e-3)                       # box "9.1e9-3.2e10" at a_s = 0.2 fm   # 2^3 (staggered ruling 2026-11-13): was 9.06056e9
    near(i["t_deepest"][1], 9.40696e9, 1e-3)   # 2^3 (staggered ruling 2026-11-13): was 3.23463e10
    for k in (0, 1):
        near(i["t_deepest"][k], i["quench_prep_t"][k] + i["t_fit_deepest"][k], 1e-12)
    near(i["ratio_deepest_to_budget"][0], 2.547, 0.002)            # "9-32 times the 2033 target"   # 2^3 (staggered ruling 2026-11-13): was 9.061
    near(i["ratio_deepest_to_budget"][1], 9.407, 0.002)   # 2^3 (staggered ruling 2026-11-13): was 32.35
    near(i["quench_prep_over_fit"][0], 12.67, 0.002)                # "13-21 times the fit evolution" (2^3; was 8.00, "8-18")
    near(i["quench_prep_over_fit"][1], 21.19, 0.002)   # was 17.88
    assert i["t_deepest"][0] > 1.5 * 1e9                           # both corners are past the 1.5x line: a co-design gap
    near(i["eps_l_deepest"][1], 1.0630e-11, 0.002)                 # box "<~ 1.1e-11" (2^3; was 3.0915e-12)
    near(i["eps_l_deepest"][0], 3.9259e-11, 0.002)                 # "3.9e-11 at the fast corner" (was 1.1037e-11)
    near(i["shot_time_deepest_s"][0] / 3600, 0.7076, 0.001)         # box "shot 0.7-2.6 h" (2^3; was 2.517, "2.5-9.0 h")
    near(i["shot_time_deepest_s"][1], 9407.1, 0.001)                # 9.4e3 s = 2.6 h at a_s = 0.2 fm (was 32346)
    # T3: the second spacing's late shot (own step, dt = 0.0075 fm/c): 755 quench steps + 48 fit steps
    assert i["second_spacing_quench_prep_steps"] == (196, 755)
    assert i["quench_prep_steps_deepest_campaign"] == (196, 755)
    near(i["t_deepest_campaign"][1], 1.2515e10, 1e-3)              # "1.3e10 on the second spacing's late shot" (2^3; was 4.3106e10)
    near(i["eps_l_deepest_campaign"][1], 7.990e-12, 0.002)        # "8.0e-12 on the second spacing" (2^3; was 2.3198e-12)
    assert i["deepest_shot_spacing_fm"] == (0.15, 0.15)
    near(i["second_spacing_own_step_t_deepest"][1], i["t_deepest_campaign"][1], 1e-12)
    # the eight corners of the (T, tau, c) band: the box quotes (300, 1/k, 2) and (225, 1/(2 pi T), 2 pi)
    tc = i["t_deepest_corners"]
    assert len(tc) == 8 and min(tc.values()) == i["t_deepest"][0] and max(tc.values()) == i["t_deepest"][1]
    near(tc[(300.0, "2piT", round(2 * math.pi, 4))], 7.097e9, 0.002)   # 2^3; was 2.4255e10
      # the ruling's "~17x" is this corner (T hi only)
    # in-situ T_c^lat diagnostic: 10 energies x 100 prep-only shots, ~6% of the first result's wall
    assert i["tc_diag_shots"] == 1000
    near(i["tc_diag_t_per_shot"][0], 2.3593e9, 0.002)                  # prep-only circuit at its own eps   # 2^3 (staggered ruling 2026-11-13): was 8.0458e9
    near(i["tc_diag_t_per_shot"][1], 8.9794e9, 0.002)   # 2^3 (staggered ruling 2026-11-13): was 3.0619e10
    near(i["tc_diag_wall_yr"][0], 0.07791, 0.005)                        # 0.25-0.97 yr (model record)   # 2^3 (staggered ruling 2026-11-13): was 0.2550
    near(i["tc_diag_wall_yr"][1], 0.2992, 0.005)   # 2^3 (staggered ruling 2026-11-13): was 0.9702
    assert 0.045 < min(i["tc_diag_over_first_result"]) and max(i["tc_diag_over_first_result"]) < 0.07   # 2^3: 4.7-6.3% (was > 5%)


def test_2033_own_depth_walls(r33):
    """One machine, each shot max(N_T x 1 us, D_T x 10 us) + 0.1 ms at the high depth schedule; every shot feeds the
    factories (F* >= 14), so no shot is charged its depth. Quench ruling: the walls carry the preparation."""
    i = r33.intermediates
    # card RF3-ch5-s216 (2026-11-12): were 3.324, 11.44 ("3.3-11 yr"); 270.9, 998.5; 54.18, 199.7; 7.478, 25.74; 0.6647
    near(i["wall_first_result_yr"][0], 1.64, 0.002)              # "4.5-16 yr"   # 2^3 (staggered ruling 2026-11-13): was 4.517
    near(i["wall_first_result_yr"][1], 4.766, 0.002)   # 2^3 (staggered ruling 2026-11-13): was 15.556
    # 2^3: the fast-end shot is depth-limited (F* = 9.6 < 10 at the high schedule), so the wall exceeds the serial
    # baseline there (was equal on 3^3)
    assert i["wall_first_result_d"][0] > i["wall_first_result_d_baseline"][0]
    assert i["wall_first_result_d"][1] > i["wall_first_result_d_baseline"][1]
    fc = i["first_result_corners"]
    assert len(fc) == 8 and 1.52 < min(fc.values()) < 1.53 and 5.99 < max(fc.values()) < 6.01   # 2^3; was 4.5-4.55, 15.95-16.0
      # (225, 1/k, 2 pi) 15.96 is not a box corner; the box quotes (225, 1/(2 pi T), 2 pi), 15.56
    near(i["campaign_two_spacings_yr_own"][0], 134.4, 0.001)      # "130-420 yr" (2^3; was 368.2, "370-1400 yr")
    near(i["campaign_two_spacings_yr_own"][1], 416.1, 0.001)   # was 1357.8
    # 2^3: the fast-end campaign is depth-limited (F* 9.5-9.6 < 10); on 3^3 it equalled the baseline with F* > 10
    assert i["campaign_two_spacings_yr_own"][0] > i["campaign_two_spacings_yr_own_baseline"][0]
    assert 9 < min(i["f_star_min_shot_hi"]) < 10
    near(i["campaign_reduction_needed_own"][0], 26.88, 0.002)     # "27-83x the five-year window" (2^3; was 73.64, "74-270x")
    near(i["campaign_reduction_needed_own"][1], 83.23, 0.002)   # was 271.6
    near(i["grid_point_yr_own"][0], 3.691, 0.002)                 # one 20% grid point takes 3.7-11 yr (2^3; was 10.16)
    near(i["grid_point_yr_own"][1], 10.72, 0.002)   # was 35.00
    assert i["grid_points_in_horizon_own"] == (1, 0)             # 2^3: one point at the fast end (was (0, 0))
    near(i["first_result_over_horizon"][0], 0.3281, 0.002)   # 2^3; was 0.9034
            # the fast corner's first result still fits the window (4.5 yr)
    assert i["first_result_over_horizon"][1] < 1   # 2^3: the first result fits at every corner (was only at the fast end)
    assert r33.wall_time_s == i["wall_campaign_s"]
    yr = m.SECONDS_PER_YEAR          # common.SECONDS_PER_YEAR since 2026-10-09 (was 3.156e7 here)
    near(i["wall_campaign_s"][1] / yr, i["campaign_two_spacings_yr_own"][1], 1e-12)


def test_2033_r9_exports(r33):
    """R9: the five depth exports are the campaign at the chapter's tau = 0.049 fm/c, shot-averaged (low/high
    schedule); the two walls are bands over the decay band."""
    i = r33.intermediates
    e = c.depth_exports(i["t_per_shot"], i["t_depth_per_shot"], i["shots_campaign"][0], 1e-6, 1e-4)
    for k in c.DEPTH_EXPORTS[:5]:
        assert i[k] == e[k]
    near(i["t_per_shot"][0], 2.86712e9, 0.001)              # over both spacings, quench prep included (was 7.422e9)   # 2^3 (staggered ruling 2026-11-13): was 1.00876e10
    near(i["f_star"][0], 9.583, 0.002)                       # were 37.26, 548.5   # 2^3 (staggered ruling 2026-11-13): was 32.07
    near(i["f_star"][1], 119.7, 0.002)   # 2^3 (staggered ruling 2026-11-13): was 401.5
    near(i["wall_serial_s"][0], i["campaign_two_spacings_yr_own_baseline"][0] * m.SECONDS_PER_YEAR, 1e-9)
    near(i["factories_for_1yr"][0], 1288, 0.002)             # was 2709   # 2^3 (staggered ruling 2026-11-13): was 3682
    near(i["floor_wall_s"][1] / 3.156e7, 134.4, 0.002)       # was 72.70   # 2^3 (staggered ruling 2026-11-13): was 114.8
    assert i["baseline_ok"] == (False, True)   # 2^3: F* 9.58 at the high schedule (was (True, True))
    assert i["wall_first_result_s"][0] < i["wall_first_result_s"][1]


def test_2033_depth_schedule_reproduces_factory_json(a):
    """step_depth_2033 at the Sigma(72x3) / H_KS reference t_rot (28.54) gives factory.json's 4.68e4 / 6.05e5 per step;
    s648: at H_I the layers scale with the multiplicities; with the Sigma(216x3) hop row (card RF3-ch5-s216) the
    step is 1.334e5 / 1.672e6 at 28.94 (7.20e4 / 1.06e6 at 28.65 with the Sigma(72x3) placeholder hop)."""
    b = dataclasses.replace(a, group_2033=c.Stated("S72x3", "t"), hamiltonian_2033=c.Stated("KS", "t"))
    lo, hi = m.step_depth_2033(b, 1.15 * math.log2(1 / math.sqrt(1e-2 / 132624000)) + 9.2)
    near(lo, 4.68e4, 0.002)
    near(hi, 6.05e5, 0.002)
    lo, hi = m.step_depth_2033(a, 1.15 * math.log2(1 / EPS_33) + 9.2)
    near(lo, 1.3123e5, 0.002)   # 2^3 (staggered ruling 2026-11-13): was 1.3339e5
    near(hi, 1.6418e6, 0.002)   # 2^3 (staggered ruling 2026-11-13): was 1.6722e6


def test_2033_campaign_split_is_stated(a):
    """ch05-mub-split (R11, option A): mu_B in {0,200,400,600} MeV, 1 of 4 at zero, Stated at app03:150. Since the
    quench ruling (2026-10-05) one thermal route serves every point, so the split moves no cost."""
    assert a.mu_b_zero_fraction.prov is c.Provenance.STATED
    assert r"$4$ points at $N_B=0$ and $12$ at $N_B=1$--$3$" in TEX.read_text()   # SR-ch5-1 (2026-11-10): canonical N_B sectors
    assert a.mu_b_zero_fraction.lo == 1.0 / a.grid_muB.lo == 0.25
    r0 = m.model(a, "2033")
    for f in (0.0, 1.0):
        b = dataclasses.replace(a, mu_b_zero_fraction=c.Uncited(f))
        r = m.model(b, "2033")
        assert r.hard_ops == r0.hard_ops and r.wall_time_s == pytest.approx(r0.wall_time_s)




# --- codesign: the box's own BE-trim line, not a separate box -----------------

def test_codesign_is_the_be_trim(rcd, r33):
    i = rcd.intermediates
    assert rcd.lq == r33.lq
    near(rcd.hard_ops[0], 6.1375e8, 0.001)                   # R-TOL at the BE's own N_rot/10 (route (a) "~2.1e9"; 1.5411e9 with the placeholder hop)   # 2^3 (staggered ruling 2026-11-13): was 2.0933e9
    near(i["eps_rot_be"], EPS_BE, 1e-12)
    near(i["t_per_rotation_papers_be"], 16.82, 0.001)   # 2^3 (staggered ruling 2026-11-13): was 17.83
    assert rcd.hard_ops[0] == rcd.hard_ops[1]
    near(rcd.breakdown_total(), rcd.hard_ops[0], 1e-9)        # the BE walk alone (no thermal preparation in route (a))
    assert {p.name for p in rcd.breakdown} == {"be_walk_S216x3"}
    near(i["eps_l_required_codesign"], 1.629e-10, 0.01)   # 2^3 (staggered ruling 2026-11-13): was 4.777e-11
    assert rcd.epsilon_l == (i["eps_l_required_codesign"], i["eps_l_required_codesign"])
    assert any("No codesign box" in n for n in rcd.notes)
    assert "NOT A BOX" in m.PUBLISHED["codesign"].src


# --- provenance and conflict surfacing ---------------------------------------

def test_register_and_costs_follow_the_group(a):
    """Changing the group changes the register and the per-link-step cost together, from groups.py."""
    b = dataclasses.replace(a, group_2033=c.Stated("S36x3", "test"), hamiltonian_2033=c.Stated("KS", "test"))
    r = m.model(b, "2033")
    assert r.intermediates["qubits_per_link_2033"] == GROUPS["S36x3"].link_qubits == 8
    assert r.intermediates["gauge_qubits_2033"] == 24 * 8 # 2^3 (staggered ruling 2026-11-13): was 3^3
    near(r.intermediates["magnetic_t_per_link_step_2033"],
         magnetic_per_link("S36x3", "KS", 3, r.intermediates["eps_rot_2033"]), 1e-12)   # its own R-TOL eps
    assert r.intermediates["fft_compiled_2033"] is True     # Sigma(36x3) has a COMPILED FFT
    near(r.intermediates["magnetic_t_per_link_step_2033"], 4.9e3, 0.02)   # the number the old '5e3' was
    # a group with no primitive table is refused, not silently priced (Z3 has a complete table since Ch. 6
    # round E added its U_FFT, 2026-09-29; Z2 has none)
    with pytest.raises(ValueError):
        dataclasses.replace(a, group_2033=c.Stated("Z2", "test"))
    with pytest.raises(ValueError):
        dataclasses.replace(a, group_2028=c.Stated("Z2", "test"))


def test_r3_fault_budget_is_a_cited_input(a):
    assert a.faults_per_shot.lo == 0.1 and a.faults_per_shot.prov is c.Provenance.CITED
    assert "R3" in a.faults_per_shot.src
    assert a.ceiling_2033_t.lo == 1e9 and a.rfi_eps_l.lo == 1e-8
    b = dataclasses.replace(a, faults_per_shot=c.Cited(1.0, "test"))
    r = m.model(b, "2033")
    near(r.epsilon_l[1], 1 / 2.54721e9, 0.001)   # 2^3; was 1 / 9.06056e9
   
    near(r.intermediates["eps_l_required_reference"], 1.597e-10, 0.01)   # 2^3 (staggered ruling 2026-11-13): was 4.684e-11


def test_s648_assumptions_are_noted(a, r33, rcd):
    """s648: the group, Hamiltonian, hop source and freezing assumption are named in the notes and Assumptions.
    Card RF3-ch5-s216 (2026-11-12): the hop note names the draft's Sigma(216x3) row; no note calls the hop a
    PLACEHOLDER any more (the retired Sigma(72x3) stand-in is mentioned as history only)."""
    for r in (r33, rcd):
        assert any("s648" in n and "Sigma(216x3)" in n for n in r.notes)
    assert any("RF3-ch5-s216" in n and "FERMION_DRAFT['S216x3']" in n and "20,068 Toffoli" in n for n in r33.notes)
    assert not any("PLACEHOLDER" in n for n in r33.notes + rcd.notes)
    assert a.nonfreezing_HI_assumed.prov is c.Provenance.ASSUMED and "arxiv_1906_11213" in a.nonfreezing_HI_assumed.src
    assert (a.beta_f_S72x3.lo, a.beta_f_S216x3.lo) == (3.18, 3.80)
    near(6 / a.beta_f_S72x3.lo, 1.89, 0.002)
    near(6 / a.beta_f_S216x3.lo, 1.58, 0.002)
    tex = TEX.read_text()
    assert r"$\beta_f=3.18(3)$~\cite{arxiv_2511_17437} and $\Sigma(216\times 3)$ at $\beta_f=3.80(5)$" in tex
    assert r"i.e.\ $g^2_f=1.89$ and $1.58$" in tex
    assert r"one extra single-plaquette term" in tex and r"$H_I$ is a different modification, an $\mathcal O(a^2)$ improvement." in tex
    assert r"\cite{arxiv_1906_11213,arxiv_2203_02330}" in tex and r"\cite{arxiv_2611_08161}" in tex


def test_instance_rows(a, r28, r33, rcd):
    rows28 = m.INSTANCE_ROWS(a, "2028", r28)
    assert len(rows28) == 2
    assert rows28[0][2][0] == pytest.approx(9.473e4, rel=1e-3) and not rows28[0][3].get("conditional")
    assert rows28[1][3]["conditional"] is True and rows28[1][3]["condition"].startswith(">=5.9x")
    assert all(x[1] == (177, 178) for x in rows28)
    rows33 = m.INSTANCE_ROWS(a, "2033", r33)
    assert len(rows33) == 1 and rows33[0][2] == r33.intermediates["t_deepest"] == r33.hard_ops
    assert not rows33[0][3].get("conditional") and "ETH" in rows33[0][3]["condition"]
    assert rows33[0][3]["route"].startswith("quench")
    assert m.INSTANCE_ROWS(a, "codesign", rcd) == []


def test_bad_inputs_rejected(a):
    with pytest.raises(ValueError):
        dataclasses.replace(a, mu_b_zero_fraction=c.Uncited(1.5))
    with pytest.raises(ValueError):
        dataclasses.replace(a, anc_2033=c.Stated(100, "x"))
    with pytest.raises(ValueError):
        dataclasses.replace(a, quench_c=c.Assumed((2.0, 1.0), "x"))
    with pytest.raises(ValueError):
        dataclasses.replace(a, quench_ramp_a_s=c.Assumed(-1.0, "x"))
    with pytest.raises(ValueError):
        dataclasses.replace(a, hamiltonian=c.Stated("XX", "x"))
    with pytest.raises(ValueError):
        dataclasses.replace(a, faults_per_shot=c.Cited(0.0, "x"))
    with pytest.raises(ValueError):
        m.model(a, "2040")


# --- the .tex still says what the model encodes -----------------------------

@pytest.mark.skipif(not TEX.exists(), reason="chapter .tex not beside the scripts tree")
@pytest.mark.parametrize("phrase", [
    # r16 sweep: box derivation rows moved to the prose; pins follow them there
    r"Logical qubits & $178$ \\",
    r"Per-shot T & first result: $9.5\times 10^{4}$ (one step of $0.045$ fm/$c$)",   # style pass 2026-10-08: box row reworded, numbers unchanged
    r"full benchmark: $5.9\times 10^{5}$ ($6$ Trotter steps of $a_t/4$ to $t_{\max} = 0.15$ fm/$c$)",   # style pass 2026-10-08: box row reworded, numbers unchanged
    r"Shots & first result $2\times 10^{3}$ ($8.5\times 10^{4}$ for $30\%$); full benchmark $7\times 10^{3}$ ($7$ measurement times) \\",   # re-pinned 2026-11-12 (batch 6): the 30% count added (shots_first_result_2028_at_30pct = 8.458e4);   # style pass 2026-10-08: box row reworded, numbers unchanged
    r"first result $1.6$ min; full benchmark $34$ min per $k$ on one machine; longest shot $\sim 0.6$ s \\",   # style pass 2026-10-08: box row reworded, numbers unchanged
    r"Logical qubits & $\sim 490$ ($436$--$536$), within the 2033 target \\",   # staggered ruling 2026-11-13 (2^3): was the 3^3 numbers   # Q5-9 (markup board 2026-11-09): box trimmed to the other chapters' rows, numbers unchanged
    r"Per-shot T & $2.5\times 10^{9}$--$9.4\times 10^{9}$, mostly the quench (Fig.~\ref{fig:budget_split}); second spacing $1.3\times 10^{10}$ \\",   # staggered ruling 2026-11-13 (2^3): was the 3^3 numbers   # re-pinned 2026-11-12 (card RF3-ch5-s216: Sigma(216x3) hop row; was 6.7e9-2.4e10, 3.2e10);   # Q5-9 (markup board 2026-11-09): box trimmed to the other chapters' rows, numbers unchanged
    r"The hydrodynamic $\tau$ has no basis at this $k$",   # re-pinned 2026-11-09 (referee cards, Q5-6): sentence became the measurement figure; the caption keeps the claim
    # retired 2026-11-13: sentence cut in the accepted cut pass (CUT3) -- "with $3.4$--$3.6\times 10^{4}$ shots per point"
    r"$4$ points at $N_B=0$ and $12$ at $N_B=1$--$3$, at two spacings, each estimated at $T=1.5\,T_c^{\rm lat}$.",   # re-pinned 2026-11-09 (referee cards)
    r"Required $\epsilon_l$ & $\lesssim 1.1\times 10^{-11}$ ($8.0\times 10^{-12}$ at the second spacing) \\",   # staggered ruling 2026-11-13 (2^3): was the 3^3 numbers   # re-pinned 2026-11-12 (card RF3-ch5-s216; was 4.2e-12, 3.2e-12);   # Q5-9 (markup board 2026-11-09): box trimmed to the other chapters' rows, numbers unchanged
    r"Shots & first result $1.5$--$2.0\times 10^{4}$; campaign $1.1$--$1.4\times 10^{6}$ \\",   # staggered ruling 2026-11-13 (2^3): was the 3^3 numbers   # Q5-9 (markup board 2026-11-09): box trimmed to the other chapters' rows, numbers unchanged
    # verifier r25 items 1, 3, 8, 9, 10
    r"the eventual first-principles (continuum, larger-volume) $\eta/s$ and $\zeta/s$ at $20\%$",
    r"so it adds up to $1\times$ the $\eta/s$ campaign's shots, and the subtraction for $\zeta$ is not counted.",
    r"\textbf{Wall time} & shot $0.7$--$2.6$ h, first result $1.6$--$4.8$ yr on one machine;",   # staggered ruling 2026-11-13 (2^3): was the 3^3 numbers   # re-pinned 2026-11-12 (card RF3-ch5-s216; was 1.9-6.6 h, 3.3-11 yr);   # re-pinned 2026-11-12 (batch 6): per-shot wall added to the row;   # style pass 2026-10-08: box row reworded, numbers unchanged
    r"campaign $130$--$420$ yr, $27$--$83\times$ the five-year window \\",   # staggered ruling 2026-11-13 (2^3): was the 3^3 numbers   # re-pinned 2026-11-12 (card RF3-ch5-s216; was 270-1000 yr, 54-200x);   # Q5-9 (markup board 2026-11-09): box trimmed to the other chapters' rows, numbers unchanged
    # R-TOL (H. Lamm 2026-09-29): the rule and each circuit's tolerance and T per rotation
    # card RF3-ch5-s216 (2026-11-12): the hop sentence, its share, the reference circuit, the BE line and the gap bullet
    r"The $\Sigma(216\times 3)$ hopping term uses unpublished gate counts~\cite{FermionPrimitives_unpub}.",
    r"At 2033 the hopping term is $43\%$ of the step.",
    r"at $1.6\times 10^{7}$ T ($6.5\times 10^{5}$ per link-step), $6.3\times 10^{9}$ T.",   # staggered ruling 2026-11-13 (2^3): was the 3^3 numbers
    r"so the longest shot is $2.5$--$9.4$ times the 2033 target of $10^{9}$ T",   # staggered ruling 2026-11-13 (2^3): was the 3^3 numbers
    r"puts every 2033 shot $2.5$--$9.4$ times above the 2033 target.",   # staggered ruling 2026-11-13 (2^3): was the 3^3 numbers
    r"the campaign is $27$--$83$ times the five-year window",   # staggered ruling 2026-11-13 (2^3): was the 3^3 numbers
    r"the 4 fm/$c$ reference circuit would be $\sim 6.1\times 10^{8}$ T.",   # staggered ruling 2026-11-13 (2^3): was the 3^3 numbers
    r"The $\Sigma(216\times 3)$ hop is $43\%$ of the 2033 step.  Recomputing its register map at both color basis changes gives $1.0\times 10^{10}$ T at 4 fm/$c$.  The inverse color basis change ($4.7\times 10^{4}$ T per link-step) may be unnecessary.",   # staggered ruling 2026-11-13 (2^3): was the 3^3 numbers
    r"(the second spacing's late circuit) is $7.6\times 10^{9}$ T.",   # staggered ruling 2026-11-13 (2^3): was the 3^3 numbers
    r"on $\Sigma(216\times 3)$ (11 qubits/link)",   # Q5-9 (markup board 2026-11-09): box trimmed to the other chapters' rows, numbers unchanged

    # E27 (2026-10-02): one machine, serial; the excess over the horizon is an algorithmic reduction, never a machine count
    # referee verifier 2026-10-04 (T3, G8)
    r"so $\Sigma(72\times 3)$ under $H_{\rm KS}$ would sit on a coarse, confined lattice.",
    r"The Hamiltonian (anisotropic) limit moves freezing to stronger coupling, so",
    # quench thermal route (author ruling 2026-10-05)
    r"\cite{arxiv_cond-mat_9403051,arxiv_0708_1324,arxiv_1509_06411,arxiv_0902_0927}",
    r"\cite{arxiv_2303_14264,arxiv_2308_16202}",
    r"$Z_2$ gauge theory up to $5\times 3$, $2T$ up to $2\times 2$ and a $Z_2$-plus-staggered-fermion chain of $16$ sites",
    r"reproduces canonical local observables to $1$--$4\%$ ($6\%$ on $3\times 3$) within $t\approx 0.6$--$2.5/T$",
])
def test_tex_still_states_the_inputs(phrase):
    text = TEX.read_text()
    assert phrase in text, f"chapter no longer says: {phrase}"


@pytest.mark.skipif(not TEX.exists(), reason="chapter .tex not beside the scripts tree")
@pytest.mark.parametrize("phrase", [
    r"T/plaq.",                                   # the retired per-plaquette figures
    r"$\lceil\log_2 216\rceil=8$ qubits/link, and",
    r"$= 648$ gauge",
    r"$20$ min per",
    r"$\sim 8\times 10^{13}$",
    r"Gibbs prep-dominated",
    r"both arms therefore require either",                       # D2: false for the mu_B>0 arm
    r"need the $10\times$ BE trim or $\gtrsim 3$ machines",    # D1: BE alone leaves 5.4 yr
    r"$1.2\times 10^{9}$ T/shot sits at the RFI",              # R-HOP: the mu_B=0 rail is 1.4x over
    r"The $\epsilon=10^{-4}$ per rotation is the report's working tolerance",   # R-TOL: no fixed tolerance
    r"$\epsilon=10^{-4}$",                                     # R-TOL: no fixed 1e-4 anywhere in the chapter
    r"contingent on the $3$--$10\times$ per-step cut",          # ruling (b): needs >= 4.4x
    r"Any one of these by 2028 brings",                        # ruling (b)
    r"by the rule of Ch.~\ref{ch:collider}: two link-controlled color moves",   # r17: R-HOP retired here
    r"split roughly evenly",                                   # r17: evolution dominates both rails
    r"only halves the on-device-thermal rail",                 # r17: BE cuts it 4.9x
    r"with no offset",                                          # E20: no slope-only gauge pricing
    r"in their convention (7 T per Toffoli",                    # E20
    r"Fermion-term rotations are priced at the full fit",      # E20: one sentence covers every rotation
    r"$4.4\times$",                                           # E20: 2028 needs >= 5.9x
    r"$14{,}314$ Toffolis",                                     # E21 (1): squish held per link, 3,672
    r"The draft counts the frame change once per site.",       # E21: replaced by drawn-but-not-counted
    r"the second uncompute of the $\Sigma(72\times 3)$ register map that the draft omits",   # E21
    r"$7.8\times 10^{4}$ T",                                  # E21: the undo is 3.2e4
    r"cuts it from $1.8\times 10^{5}$ to $1.1\times 10^{5}$",   # E21: share='link' is now the price
    r"just inside",                                            # E21: 3.5 yr, comfortably inside
    r"$8.7\times 10^{9}$",                                    # E21: mu_B>0 shot 6.3e9
    r"machines",                                               # E27: no machine count anywhere (one machine, serial)
    r"serialized",                                             # E27: "serial on one machine"
    r"to fit the 5-year horizon; the BE trim alone leaves",    # E27: the old requirements row
    # r25 (R1, R4, R5, R10, 2026-10-02): each shot at its own depth; two tiers; Gibbs cut required; 2028 workspace
    r"Logical qubits & $162$",
    r"$99$ min",
    r"$\sim 18$ yr",
    r"$\sim 37$ yr",
    r"$7.3\times$",
    r"carry $78\%$",
    r"$\mu_B{=}0$ rail: $5.3\times 10^{9}$ T",
    r"$t_{\max} = 4$ fm/$c$ \\",
    r"ample for the decay-rate fit",
    r"$10^3$ shots per grid point",
    r"estimator-variance factor of $\sim 40$",
    r"$\lesssim 1.9\times 10^{-11}$",
    r"$\eta/s$ to 20\%]",
    r"$10^4$ ($10^3$/timeslice $\times\, 10$ timeslices)",
    r"evolution still dominates",
    # verifier r25
    r"bounded against the $k\to 0$ Kubo channel",
    r"rides the $\eta/s$ campaign",
    r"a superconducting-modality assumption",
    r"Shear and bulk viscosities $\eta$, $\zeta$ \\",
    r"$10^3$/timeslice",
    r"per point $\times$ $16$ points",
    r"the 2028 shot's",
    # referee report 2026-10-04 (T1-T3, G7, G8, editorial)
    r"rail",                                                    # terminology: route
    r"ceiling",                                                 # G3: reference budget
    r"a floor",                                                 # lower bound
    r"The draft draws",                                         # drafting history
    r"We retire",                                               # drafting history
    r"$5\%$ assumed by analogy",                                # T3: no subgroup estimate
    r"$20$--$40\%$ at $L=0.6$ fm",                              # T3: finite volume uncontrolled
    r"$30$--$50\%$ at $a_s=0.2$ fm",                            # T3: discretization not quantified
    r"(the $12^3$ extension) \\",                               # T1: 12^3 still has k_min ~ 3.2 T
    r"This finite-$k$ law is what is measured here",            # T1
    r"licenses $3$--$10\times$ larger $a_t$ for free",          # G7: sqrt of the slack at 2nd order
    r"to hold Trotter error below",                             # G7: chosen step
    r"stochastic / E$\rho$OQ",                                  # G8: the estimator is stated
    r"would mean $100\times$",                                   # G8 verifier: source estimate n_eff^-V, not 0.1
    r"assumed improvements, estimated at",                      # G7 verifier: assumed, not derived
    r"(method validation at fixed volume;",                     # T3: the two spacings differ in volume
    # 2026-10-05: T1 fast end, T3 second spacing, G6 phase bracket, G8 computed penalty
    r"the fast end is optimistic",                              # T1: fast end now 1/k
    r"fast end hydrodynamic",                                   # T1
    r"the second priced at the first's per-shot cost",          # T3: own step
    r"$\geq 0.17$",                                             # G6: hop-free value is not a bound (as Ch. 10)
    r"average phase of $1$ at zero coupling",                   # G6: Gauss's law keeps a singlet projection
    r"conditional on signal-to-noise mitigation",               # G8: penalty computed, not a mitigation away
    # quench ruling (2026-10-05): E-rho-OQ as a route and the Gibbs sampler as a priced route are gone
    r"eq:erhooq",
    r"\langle\delta_{ij}\rangle_\rho",
    r"Gibbs prep",
    r"Gibbs-prep",
    r"heat-reset",
    r"$10\times$ Gibbs",
    r"$\mu_B{>}0$ route",
    r"$\mu_B{=}0$ route",
    r"Thermal-pure-quantum",
    r"causal",
    r"moves freezing to stronger coupling~\cite{Gustafson_S648_inprep}",   # ruling (E): clause uncited
    r"n_{\rm eff}^{-V}",                                        # G8: replaced by the computed value
    # s648 (2026-10-05): Sigma(216x3) / H_I at a_s = 0.2 fm; the rulings2 g^2 = 3.5 / 0.69 fm instance retired
    r"$0.69$ fm",
    r"$g^2=3.5$",
    r"hot starts",
    r"conditional on a fast $\Sigma(72\times 3)$ transform",
    r"compiled dense transform",
    r"($|G|=216$, 9 qubits/link)",
    r"$1072$--$1172$",
    r"we charge them their circuit depth",
    # card RF3-ch5-s216 (2026-11-12): the Sigma(72x3) placeholder hop and its numbers are gone from the chapter
    r"No hopping-term circuit exists",
    r"placeholder",
    r"likely low",
    r"$23\%$",
    r"$6.7\times 10^{9}$--$2.4\times 10^{10}$",
    r"$3.3$--$11$ yr",
    r"$270$--$1000$ yr",
    r"$54$--$200",
    r"$7$--$24$ times",
    r"$4.8\times 10^{5}$ per link-step",
])
def test_tex_no_longer_states_the_retired_inputs(phrase):
    assert phrase not in TEX.read_text(), f"chapter still says: {phrase}"


@pytest.mark.skipif(not TEX.exists(), reason="chapter .tex not beside the scripts tree")
def test_tex_assumes_one_machine():
    """E27 (H. Lamm 2026-10-02): the report never assumes more than one quantum machine. The only 'parallel' left in
    the chapter is the magic-state factories inside the one machine; every 'machine' is singular."""
    import re
    tex = TEX.read_text()
    for pat in (r"\bmachines\b", r"machine-years", r"shot-parallel", r"embarrassingly", r"fold parallelism",
                r"across machines", r"per machine", r"\d\s*\$?\s*machines"):
        assert not re.search(pat, tex), f"multi-machine wording: {pat}"
    assert not re.search(r"parallel(?! magic-state factories)", tex) and tex.count("parallel") <= 1   # CUT-pinned 2026-11-11: the schedule sentence with the factories went
    assert tex.count("on one machine") >= 3   # Q5-9: the box campaign row no longer repeats it


def test_shots_from_fisher_not_the_old_gloss(a):
    """r25: shots come from the Fisher information of the two-point fit at Cbar = 0.16 (Stated; 1/Cbar^2 = 39, the
    old '~40' factor). The r23 gloss '1e3 = eps^-2 x ~40' is retired from the chapter; shots_per_gridpt is a record."""
    assert a.cbar.prov is c.Provenance.STATED and a.cbar.lo == 0.16
    near(1 / a.cbar.lo ** 2, 39.06, 0.001)
    assert a.shots_per_gridpt.lo == 1e3                      # record (r23 reference campaign)
    assert a.target_first.lo == 0.3 and a.target_campaign.lo == 0.2
    tex = TEX.read_text()
    assert r"\epsilon^{-2}=25$ at $\epsilon=20\%$" not in tex
    assert r"$\bar C\approx 0.16$" in tex


def test_r25_verifier_2028_first_result_at_30pct():
    """Verifier r25 item 1: the 2028 first result at R1's 30% on a relative fall-off delta = 0.1 (prose only)."""
    i = m.model(m.Assumptions(), "2028").intermediates
    assert i["shots_first_result_2028_at_30pct"] == pytest.approx(2 * (1 - 0.16**2) / (0.16**2 * 0.3**2 * 0.1**2))
    assert i["shots_first_result_2028_at_30pct"] == pytest.approx(8.458e4, rel=1e-3)      # "~8.5e4"
    assert i["wall_first_result_2028_at_30pct_h"] == pytest.approx(1.115, rel=1e-3)        # "~1.1 h"
    assert i["wall_first_result_2028_at_30pct_all_one_step_h"] == pytest.approx(2.228, rel=1e-3)
    assert i["shots_first_result_2028"] == 2e3                                             # box unchanged


def test_r25_verifier_r9_naming_and_exports_note():
    a = m.Assumptions()
    assert not hasattr(a, "logical_cycle_us") and not hasattr(a, "logical_cycle_s")    # R9: legacy names dropped
    r = m.model(a, "2033")
    rows = m.INSTANCE_ROWS(a, "2033", r)
    assert "tau = 1/k = 0.0637" in rows[0][3]["exports_note"]   # 2^3; was 0.0955


def test_referee_t1_t3_g8_numbers(r33):
    i = r33.intermediates
    (f0, f1), (s0, s1) = i["corr_rel_err_first_at_fit_times"]
    near(f0, 0.19, 0.03)                                                 # "19% at the early time"
    near(s0, 0.19, 0.03)
    assert 0.360 <= s1 <= 0.365 and 0.325 <= f1 <= 0.330                 # 33-36% at the late one (2^3; fast corner was 36%)
    near(i["tau_free_streaming_fm"], 0.0637, 0.001)                       # "1/k ~ 0.06 fm/c" (2^3; was 0.0955)
    near(i["second_spacing_kmin_over_T_on_3cubed"][0], 18.4, 0.005)       # "k_min ~ 14-18 T" on 2^3 (key name kept; was 12.2 on 3^3)
    near(i["second_spacing_kmin_over_T_on_3cubed"][1], 13.78, 0.005)   # was 9.18
    lo, hi = i["second_spacing_mpiL_on_3cubed"]
    assert 0.45 <= lo <= 0.46 and 0.60 <= hi <= 0.61                    # "m_pi L ~ 0.5--0.6" (2^3; was 0.7--0.9)
    assert i["fixed_volume_4cubed_lq"] == (2788, 2888)                   # "2788--2888 LQ" (11 qubits/link)
    assert i["fixed_volume_4cubed_a_s_fm"] == 0.1                        # 2^3: the fixed 0.4 fm box is 4^3 at 0.1 fm (was 0.15 on 3^3)


def test_referee_verifier_t3_g8_numbers(r33):
    """Verifier (2026-10-04): the second spacing is an Assumed example; priced at its own step (Delta t ||H|| held,
    step rule at that spacing) the slow late shot needs 48 fit steps at H_I after a 755-step quench (80 after 1054 at
    the retired 160 MeV), and the campaign's slow end is ~1400 yr (1120 yr with the second spacing at the first's
    per-shot cost; ~1000 and 820 yr with the Sigma(72x3) placeholder hop, before card RF3-ch5-s216)."""
    a = m.Assumptions()
    assert a.second_spacing_a_s_fm.prov is c.Provenance.ASSUMED and a.second_spacing_a_s_fm.lo == 0.15
    i = r33.intermediates
    assert i["second_spacing_own_step_fit_steps"][1] == (10, 40)          # grid 35 = 4/3 x 26; the H_I rule needs 40 (2^3; was 48)
    near(i["campaign_yr_second_spacing_own_step"][1], 416.1, 0.01)   # 2^3; was 1357.8 (998.5 placeholder hop)
    assert i["campaign_yr_second_spacing_own_step"][1] > i["campaign_per_spacing_yr_own"][1] * 2


# --------------------------------------------------------------------------- #
# Referee G7 (2026-10-04): Trotter error at the chosen 2033 step
# --------------------------------------------------------------------------- #

def _free_staggered_infidelity(dims, m_lat, dt, n, temp):
    """1 - |Tr rho U^dag S^N| for ONE free staggered field (U = 1), exact: single-particle Strang step over the mass
    and the 2d hop groups (direction, parity of x_mu), many-body overlap det((1 - n_F) + n_F w)."""
    import itertools
    np = pytest.importorskip("numpy")
    sites = list(itertools.product(*[range(L) for L in dims]))
    idx = {s: i for i, s in enumerate(sites)}
    V = len(sites)
    M = np.diag([m_lat * (-1) ** sum(s) for s in sites]).astype(complex)
    groups = []
    for mu in range(len(dims)):
        for p in (0, 1):
            G = np.zeros((V, V), complex)
            for s in sites:
                if s[mu] % 2 != p:
                    continue
                t = list(s)
                t[mu] = (t[mu] + 1) % dims[mu]
                eta = (-1) ** sum(s[:mu])
                G[idx[s], idx[tuple(t)]] += -0.5j * eta
                G[idx[tuple(t)], idx[s]] += 0.5j * eta
            groups.append(G)

    def ex(X, tau):
        e, v = np.linalg.eigh(X)
        return (v * np.exp(-1j * e * tau)) @ v.conj().T
    terms = [M] + groups
    S = ex(terms[-1], dt)
    for X in reversed(terms[:-1]):
        h = ex(X, dt / 2)
        S = h @ S @ h
    H = sum(terms)
    e, v = np.linalg.eigh(H)
    f = 1 / (1 + np.exp(e / temp)) if temp > 0 else (e < 0).astype(float)
    nF = (v * f) @ v.conj().T
    w = ex(H, n * dt).conj().T @ np.linalg.matrix_power(S, n)
    return 1 - abs(np.linalg.det(np.eye(V) - nF + nF @ w))


def test_g7_free_field_reproduces_ch9():
    """The exact free-gauge-field error reproduces Ch. 9's r22 numbers (3^3, 3 colors, T a = 0.5, t = 10 a)."""
    for n, err in ((450, 0.0044), (200, 0.022), (107, 0.078)):
        near(m.free_gauge_trotter_error((3, 3, 3), 3, 10 / n, n, 0.5), err, 0.03)


def test_g7_free_field_matches_fock_space():
    """The closed form of one mode against a truncated Fock-space product of the three exponentials."""
    np = pytest.importorskip("numpy")
    w, dt, n, temp, D = 2.4, 0.3, 7, 1.1, 90
    a_ = np.diag(np.sqrt(np.arange(1, D)), 1)
    x, p = (a_ + a_.T) / np.sqrt(2 * w), 1j * np.sqrt(w / 2) * (a_.T - a_)

    def ex(X, tau):
        e, v = np.linalg.eigh(X)
        return (v * np.exp(-1j * e * tau)) @ v.conj().T
    S = ex(w * w * x @ x / 2, dt / 2) @ ex((p @ p).real / 2, dt) @ ex(w * w * x @ x / 2, dt / 2)
    U = np.diag(np.exp(-1j * w * (np.arange(D) + 0.5) * n * dt))
    rho = np.diag(np.exp(-w * np.arange(D) / temp))
    rho /= rho.trace()
    K = 60                                       # keep the trace away from the truncation edge
    f = abs(np.trace((rho @ U.conj().T @ np.linalg.matrix_power(S, n))[:K, :K]))
    near(m.free_gauge_mode_overlap(w, dt, n, temp), f, 1e-6)


def test_g7_2033_step_numbers(r33):
    """Step rule (2026-10-04): the state-dependent estimate sets the step at eps = 0.1, the exact free field is the
    check, the worst case is quoted once. s648: at H_I both use the tree-level improved dispersion (Lambda_sd 45.8 vs
    32.8 at H_KS, 1.4x); temperature ruling 2026-10-05 (T a_s = 0.30 / 0.23 at the two corners): the deepest fast
    shot runs 19 steps (0.092), the deepest slow shot 32 (0.098); free field 0.015-0.019; worst case 64-68 at g^2 = 1
    with the H_I norms."""
    i = r33.intermediates
    near(i["trotter_dt_over_a_s"], 0.05, 1e-9)
    near(i["trotter_temp_lat"][0], 0.304, 0.002) and near(i["trotter_temp_lat"][1], 0.228, 0.002)
    assert i["trotter_rule"] == "state" and i["trotter_steps_deepest"] == (12, 27)   # 2^3; was (19, 32)
    s0, s1 = i["trotter_sd_err_deepest"]
    assert round(s0, 3) == 0.048 and round(s1, 3) == 0.096 and i["trotter_rule_meets_eps"]   # 2^3; was 0.092, 0.098
    g0, g1 = i["trotter_sd_err_deepest_grid"]
    assert g0 < 0.1 < g1                                                      # 2^3: the a_t/10 grid fails at the slow end only (was both)
    assert i["trotter_sd_steps_deepest_grid"] == (9, 27)                      # the fewest steps that meet it (2^3; was (19, 32))
    w0, w1 = i["trotter_worst_err_deepest"]
    assert round(w0) == 14 and round(w1) == 28   # 2^3; was 64, 68
    assert min(i["trotter_worst_err_min_g2_deepest"]) > 1                   # no g^2 rescues the bound
    f0, f1 = i["trotter_free_err_deepest"]
    assert round(f0, 3) == 0.020 and round(f1, 3) == 0.018   # 2^3; was 0.015, 0.019
    assert i["trotter_meets_eps"]
    for lam, lamk in zip(i["trotter_sd_lambda"], i["trotter_sd_lambda_ks"]):
        near(lam / lamk, 1.54, 0.003)                                          # "raises the estimate 1.5x over H_KS" (2^3; was 1.40)
        near(lamk, 20.79, 0.001)   # 2^3; was 32.76
    tex = TEX.read_text()
    assert r"\frac{1}{g^2 a}\sum_{\Box}" in tex
    assert "raises the estimate $1.5\\times$ over $H_{\\rm KS}$" in tex   # staggered ruling 2026-11-13 (2^3): was the 3^3 numbers
    assert "the fundamental carries the SU(3) Casimir $4/3$" in tex
    # the 4 fm/c reference shot is a pricing basis: the rule would need 1601 steps ("1.6e3"; 1914, "1.9e3", on 3^3)
    assert i["trotter_sd_steps_ref"][0] == 1601   # 2^3 (staggered ruling 2026-11-13)
    # Lambda_sd at T a_s = 0.23-0.30 is almost all zero-point: the T = 0 value differs by < 1e-2 (verifier 2026-10-04)
    lam0 = m.trotter_state_dependent((2, 2, 2), 8, 0.0, 1.0, 0.05, 0.1, "I")["lambda"]   # 2^3 (was (3, 3, 3), 45.8)
    assert max(abs(x - lam0) for x in i["trotter_sd_lambda"]) < 1e-2 and round(lam0, 1) == 32.0
    assert round(m.trotter_state_dependent((3, 3, 3), 8, 0.0, 1.0, 0.05, 0.1, "I")["lambda"], 1) == 45.8   # the 3^3 record
    lam0k = m.trotter_state_dependent((3, 3, 3), 8, 0.0, 1.0, 0.05, 0.1)["lambda"]
    assert round(lam0k, 1) == 32.8                                            # the H_KS default is unchanged (Ch. 6 uses it)
    assert "with the estimate taken in the thermal state" in tex   # Q5-5: the rule is defined once in Ch. 1


def test_g7_improved_dispersion_matches_ch6():
    """s648: the H_I free dispersion is the same tree-level Symanzik reading Ch. 6 uses."""
    from estimates import ch06_collider as ch6
    for dims in ((3, 3, 3), (4, 4, 6)):
        assert m.free_gauge_frequencies(dims, "I") == pytest.approx(ch6.free_gauge_frequencies_h(dims, "I"))
        assert m.free_gauge_frequencies(dims) == pytest.approx(ch6.free_gauge_frequencies_h(dims, "KS"))
    # worst case at H_I: electric half-range x 4/3, plaquette part of b x 23/12
    ks = m.trotter_worst_case(81, 3, 1.0, 8 / 3, 4.5, 1.0, 0.05, 0.1)
    hi = m.trotter_worst_case(81, 3, 1.0, 8 / 3, 4.5, 1.0, 0.05, 0.1, "I")
    e, b = 0.5 * 0.5 * 8 / 3 * 4 / 3, 4.5 + 4 * 6.0 * 23 / 12
    near(hi["lambda"], 81 * (e * e * b / 3 + e * b * b / 6), 1e-12)
    assert hi["lambda"] > ks["lambda"]


def test_g7_2028_step_rule(r28):
    """Step rule at 2028 (2T on 4^2, 3 color copies, T a_s = 0.75): Delta t = a_t gave 0.88 (state-dependent) / 0.77
    (free field) after one step; six steps of a_t/4 reach 0.15 fm/c at 0.083 (free field 0.044); the one-step first
    result at 0.045 fm/c is 0.081 (largest one-step Delta t 0.048 fm/c). The T per shot does not change."""
    i = r28.intermediates
    assert round(m.free_gauge_trotter_error((4, 4), 3, 0.5, 1, 0.75), 1) == 0.8
    assert round(i["trotter_2028_sd_err_old"][0], 2) == 0.88 and i["trotter_2028_sd_steps_old_tmax"] == 44
    assert round(i["trotter_2028_sd_err_full"], 3) == 0.083 and round(i["trotter_2028_free_err_full"], 3) == 0.044
    assert round(i["trotter_2028_sd_err_first"], 3) == 0.081 and i["trotter_2028_free_err_first"] < 0.1
    assert round(i["trotter_2028_dt_first_max_fm"], 3) == 0.048
    near(i["t_per_shot_2028"], 5.939e5, 1e-3) and near(i["t_first_result_shot_2028"], 9.473e4, 1e-3)
    tex = TEX.read_text()
    # re-pinned 2026-11-09 (referee cards, Q5-6): both sentences became rows of Table tab:ch05_trotter
    assert r"2028 benchmark, six steps of $a_t/4$ & $6$ & $0.08$ & $0.04$ & -- \\" in tex
    assert r"2028 first result, one step of $0.045$ fm/$c$ & $1$ & $0.08$ & $0.08$ & -- \\" in tex
    assert r"One 2028 step may be up to $0.048$ fm/$c$" in tex


def test_g7_fermions_negligible_in_2033_check():
    """Three free staggered fields x 3 colors change the free-field error by about 1% (27 sites, density from 4^3)."""
    dens = _free_staggered_infidelity((4, 4, 4), 0.023, 0.05, 36, 0.162) / 64
    gauge = m.free_gauge_trotter_error((3, 3, 3), 8, 0.05, 36, 0.162) ** 2 / 2
    assert 9 * 27 * dens < 0.05 * gauge


def test_eham_normalization():
    """eq:Hks normalization (2026-10-04, Claude's decision): for Sigma(72x3), E^a E^a is the Cayley-graph Laplacian
    f(rho) = |Gamma| - Re sum_Gamma chi_rho / dim rho of arxiv_2511_17437 (Gamma = elements of maximal Re Tr), scaled
    so the 3 carries C_F = 4/3. Built from the paper's generators: |G| = 216, Gamma = 54 elements of trace 1,
    f(3) = 36 and f(8) = 54 as in tab:eham, so f(8)/f(3) = 1.5 (not the Casimir 9/4) is a property of the group.
    The worst-case bound uses the top of the spectrum, f = 72 -> (4/3) 72 / 36 = 8/3."""
    np = pytest.importorskip("numpy")
    w = np.exp(2j * np.pi / 3)
    gens = [np.diag([1, w, w * w]), np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0]], complex),
            np.array([[1, 1, 1], [1, w, w * w], [1, w * w, w]]) / (np.sqrt(3) * 1j),
            np.array([[1, 1, w * w], [1, w, w], [w, 1, w]]) / (np.sqrt(3) * 1j)]
    G = {tuple(np.round(np.eye(3).flatten(), 6)): np.eye(3, dtype=complex)}
    front = list(G.values())
    while front:
        new = []
        for A in front:
            for g in gens:
                B = A @ g
                k = tuple(np.round(B.flatten(), 6))
                if k not in G:
                    G[k] = B
                    new.append(B)
        front = new
    els = [M for M in G.values() if not np.allclose(M, np.eye(3))]
    assert len(G) == 216
    mx = max(np.trace(M).real for M in els)
    gam = [M for M in els if np.isclose(np.trace(M).real, mx)]
    assert len(gam) == 54 and np.isclose(mx, 1.0) and all(np.isclose(np.trace(M), 1.0) for M in gam)
    f3 = len(gam) - sum(np.trace(M) for M in gam).real / 3
    f8 = len(gam) - sum(abs(np.trace(M)) ** 2 - 1 for M in gam) / 8
    assert np.isclose(f3, 36) and np.isclose(f8, 54)
    a = m.Assumptions()
    near(a.sigma72_emax_over_fund.lo * 4 / 3, 8 / 3, 1e-12)               # 72/36 x C_F (s648: a proxy for Sigma(216x3))
    assert "is its Cayley-graph Laplacian~\\cite{arxiv_2511_17437}, scaled so that the fundamental carries the SU(3) Casimir $4/3$." in TEX.read_text()


# --------------------------------------------------------------------------- #
# 2026-10-05 (referee T1, T3, G6, G8 follow-up)
# --------------------------------------------------------------------------- #

def test_t1_free_streaming_1e_time(r33):
    """T1 verifier (2026-10-05): 1/k is a dephasing scale; free streaming reaches 1/e at kt = 2.93 (0.28 fm/c),
    beyond the slow end 1/(2 pi T) = 0.140 fm/c (225 MeV; 0.196 at the retired 160 MeV)."""
    i = r33.intermediates
    near(i["free_streaming_1e_kt"], 2.926, 0.001)
    near(i["free_streaming_1e_fm"], 0.186, 0.005)   # 2^3 (staggered ruling 2026-11-13): was 0.279
    assert i["free_streaming_1e_fm"] > i["tau_band_fm"][1]


def test_t1_fast_end_is_free_streaming(a, r33):
    """T1: the fast end of the decay band is the collisionless dephasing time 1/k_min, not the hydrodynamic
    T/((eta/s) k^2), which stays a record."""
    assert a.tau_fast_1k.prov is c.Provenance.ASSUMED and a.tau_fast_1k.lo == 1.0
    i = r33.intermediates
    near(i["tau_band_fm"][0], i["tau_free_streaming_fm"], 1e-12)
    near(i["tau_decay_at_eta_s_assumed_fm"][0], 0.0308, 0.002)           # record: "0.07-0.09 fm/c ... has no basis"   # 2^3 (staggered ruling 2026-11-13): was 0.0693
    (s0, s1), (f0, f1) = i["test_time_sigma_lnC"], i["test_time_campaign_fraction"]
    assert 0.26 <= min(s0, s1) and max(s0, s1) <= 0.30                   # "+-0.3 in ln C" (2^3: 0.26-0.29; was 0.28-0.30)
    assert round(100 * f1) == 19 and round(100 * f0) == 19                # "add a fifth to the campaign" (2^3: 19%; was 20%)


def test_t3_second_spacing_at_its_own_step(a, r33):
    """T3: the campaign is the first spacing plus the second at its own step (Delta t ||H|| held); the old convention
    (second at the first's per-shot cost) is a record."""
    assert a.second_spacing_own_step.prov is c.Provenance.ASSUMED and a.second_spacing_own_step.value is True
    i = r33.intermediates
    o = m.own_depth_2033(a)
    yr = a.seconds_per_year.lo
    for k in (0, 1):
        near(i["campaign_two_spacings_yr_own"][k],
             (o["r"][k]["per_spacing_s"] + o["second"]["r"][k]["per_spacing_s"]) / yr, 1e-12)
    assert i["second_spacing_own_step_fit_steps"] == ((4, 12), (10, 40))   # 2^3; was ((5, 19), (10, 48))
    near(i["campaign_yr_second_spacing_same_cost"][1], 343.15, 0.001)    # the pre-T3 convention (record; 2^3, was 1120.0)
    b = dataclasses.replace(a, second_spacing_own_step=c.Assumed(False, "x"))
    rb = m.model(b, "2033")
    near(rb.intermediates["campaign_two_spacings_yr_own"][1], 343.15, 0.001)   # 2^3; was 1120.0
    assert rb.hard_ops == r33.hard_ops                                   # the box shot is the first spacing's either way
    near(rb.intermediates["t_deepest_campaign"][1], r33.hard_ops[1], 1e-12)
    assert i["deepest_shot_spacing_fm"] == (0.15, 0.15)


# --------------------------------------------------------------------------- #
# Thermal state by quench (author ruling 2026-10-05): E-rho-OQ and the Gibbs sampler retired as priced routes
# --------------------------------------------------------------------------- #

def test_quench_assumptions_and_2028_record(a, r28):
    """quench_c = (2, 2 pi) and a one-a_s ramp are Assumed with the ETH sources; the 2028 box stays non-thermal and
    the prose prices the quench there as a record: 8 + 22-68 = 30-76 prep steps in front of the one-step first
    result, 3.2e6-8.1e6 T, 32-81x the 1e5-T target (ruling (C))."""
    assert a.quench_c.prov is c.Provenance.ASSUMED and a.quench_c.lo == 2.0 and a.quench_c.hi == pytest.approx(2 * math.pi)
    for k in ("arxiv_cond-mat_9403051", "arxiv_0708_1324", "arxiv_1509_06411", "arxiv_0902_0927", "arxiv_2303_14264",
              "arxiv_2308_16202", "arxiv_2011_09814"):
        assert k in a.quench_c.src
    assert a.quench_ramp_a_s.prov is c.Provenance.ASSUMED and a.quench_ramp_a_s.lo == 1.0
    assert not hasattr(a, "gibbs_prep_t") and not hasattr(a, "eroq_prep_t") and not hasattr(a, "eroq_shot_penalty")
    q = m.quench_steps(a, 0.01, 160.0, 2.0)                               # the retired 160 MeV record
    assert (q["ramp"], q["th"], q["prep"]) == (20, 247, 267)
    assert m.quench_steps(a, 0.01, 300.0, 2.0)["prep"] == 152 and m.quench_steps(a, 0.01, 225.0, 2 * math.pi)["prep"] == 572
    i = r28.intermediates
    near(i["quench_inv_T_steps_2028"], 10.68, 1e-3)
    assert i["quench_ramp_steps_2028"] == 8 and i["quench_th_steps_2028"] == (22, 68)
    assert i["quench_prep_steps_2028"] == (30, 76)
    near(i["quench_first_result_shot_2028"][0], 3.1895e6, 1e-3)          # "3.2e6-8.1e6 T"
    near(i["quench_first_result_shot_2028"][1], 8.0888e6, 1e-3)
    near(i["quench_first_result_over_rfi_2028"][0], 31.9, 0.003)         # "32-81x the 1e5-T target"
    near(i["quench_first_result_over_rfi_2028"][1], 80.9, 0.003)
    assert r28.hard_ops == (pytest.approx(9.473e4, rel=1e-3), pytest.approx(5.939e5, rel=1e-3))   # box unchanged
    assert i["shots_first_result_2028"] == 2e3 and i["shots_2028"] == 7e3


def test_quench_shot_is_one_circuit(a, r33):
    """Each 2033 shot is one circuit, prep + fit, at its own R-TOL eps: the prep T is the per-link-step price at that
    eps times the prep steps, and the fit T the rest; the shots and the register do not move with the quench."""
    i = r33.intermediates
    o = m.own_depth_2033(a)
    for r in o["r"]:
        for p in r["pts"]:
            e = m.evolution_2033(a, p["prep_steps"] + p["steps"])
            near(p["t"], e["t"], 1e-12)
            near(p["t_prep"], e["per_link"] * 24 * p["prep_steps"], 1e-12) # 2^3 (staggered ruling 2026-11-13): was 3^3
            near(p["t_prep"] + p["t_fit"], p["t"], 1e-9)
    b = dataclasses.replace(a, quench_c=c.Assumed((1.0, 1.0), "x"))
    rb = m.model(b, "2033")
    assert rb.shots == r33.shots and rb.lq == r33.lq
    assert rb.intermediates["quench_th_steps"] == (66, 88)              # ceil(1/T) at 300 / 225 MeV
    assert rb.hard_ops[1] < r33.hard_ops[0]                               # c = 1 at both corners: cheaper than c = 2
    tex = TEX.read_text()
    assert "E$\\rho$OQ" in tex and tex.count("E$\\rho$OQ") == 1       # one sentence names the method (ruling (A))
    assert tex.count("KMS") == 1


def test_systematics_row_2026_10_12(a):
    """Author decision 2026-10-12 (option 1): the three systematics the referee pass removed come back as model
    differences scaled from Euclidean lattice results, in the 2033 box and the text, with citations."""
    import math
    i = m.model(a, "2033").intermediates
    s = i["systematics_2033"]
    assert s["discretization_0p2fm"] == (0.5, 1.0) and s["discretization_0p15fm"] == (0.35, 0.7)
    assert s["discretization_improved_0p2fm"] == (0.3, 0.6)
    assert s["finite_volume_floor"] == 0.5 and s["subgroup"] == (0.05, 0.1)
    lo, hi = i["pi_a_s_T_sq_2033"]
    assert math.isclose(lo, 0.51, rel_tol=0.02) and math.isclose(hi, 0.91, rel_tol=0.02)     # "(pi a_s T)^2 = 0.5-0.9"
    for lo_, hi_ in (s["discretization_0p2fm"], s["discretization_0p15fm"], s["subgroup"]):
        assert 0 < lo_ < hi_ <= 1
    if not TEX.exists():
        pytest.skip("chapter .tex not on disk")
    tex = TEX.read_text()
    for frag in ("Systematics & model differences and not error bars (Table~\\ref{tab:ch05_systematics}); finite volume not "
                 "quantified at $LT<1$",   # Q5-9: the percentages stay in Table 5.2
                 "$50$--$100\\%$ at $0.2$ fm, $35$--$70\\%$ at $0.15$ fm; $30$--$60\\%$ and $20$--$40\\%$ with an improved taste sector",
                 "not quantified",   # SR-ch5-3 (2026-11-10)
                 "Neither spacing is in the $a^2$ scaling regime, since $(\\pi a_s T)^2=0.5$--$0.9$ at $0.2$ fm, so the two give neither an $a^2$ slope nor a continuum limit.",
                 # retired 2026-11-13: sentence cut in the accepted cut pass (CUT3) -- ", and we assign it no percentage" (prose clause; the Table 5.2 cell stays)
                 "whose relation to infinite-volume $\\eta/s$ is uncontrolled at order one.",
                 "a few percent or less & no percentage \\\\",
                 "\\cite{arxiv_1802_07718}", "\\cite{arxiv_2401_01169}", "\\cite{Bhanot_1982_S648,arxiv_1906_11213}"):
        assert frag in tex, frag

