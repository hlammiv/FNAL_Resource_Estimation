"""Ch. 10 (app07_finite_density.tex): every named intermediate the prose states.

The generic contract and box-reproduction tests live in test_common.py. These
pin each intermediate number the chapter prints on the way to the box, so an
edit to one number without the others fails loudly. Tolerances are the
chapter's own rounding ('~180' for 176, '~26' for 26.3, '~8 yr' for 7.9).

Updated 2026-09-28 for the round-A rulings (apply_log/ch10.md): R1 compiled
widths 8/9, R2 the 2028 step repriced from the papers (one step under the cap),
R3 0.1 faults per shot, R4 carry exact; plus the mechanical batch
(qsvt-lower-end, horizon, v3cubed-width, fft-citation, rounding, link-hop-count).
Chapter statements that the stated inputs still do NOT reproduce are marked
xfail(strict=True) with the NEEDS_AUTHOR.md item.

Updated 2026-09-29 for R-HOP report-wide (apply_log/rhop_ch10.md): the staggered hop is
priced by the shared groups.rhop_link; the cap-filling allowances are gone and the totals
follow (2028 89,854.5 T; 2033 c_BE 361,353.9, 7.227e8 T/shot).

Updated 2026-09-29 for R-TOL (apply_log/rtol_ch10.md): total synthesis error 1e-2 per shot,
each circuit's eps_rot = common.eps_rot_for(N_rot). 2028: N_rot 3,740, eps 1.635e-3, 72,515.5 T;
2033: N_rot 2.28e7, eps 2.094e-5, c_BE 390,923.3, 7.818e8 T/shot; QSVT rail 4.046e8-8.230e8.

Updated 2026-09-29 for R11-R15 (apply_log/r11_ch10.md): stretch ancilla carries the 2033
250-750 band (2554-3054 LQ, box '~2600-3100', '~2.6-3.1x'); the 2028 ~100 ancilla is a sized
budget with 8-24 paper workspace + 1 Hadamard qubit accounted (register 85-101 of 176);
the 5% on p is not separately costed (chi4 budget binds).

Updated 2026-10-01 for r17 (apply_log/r17_ch10.md): the hop is groups.hop_link_cost on Sigma(36x3)
(the unpublished Fermion_Primitives counts with the estimated undo, squish uncompute, per-field color
rotations, MBU, HWP phasing, full-fit rotations); the staggered mass is groups.staggered_mass_site
with mu_B N_B folded in. 2028: N_rot 10,908, eps 9.575e-4, step 751,571 T = 7.5x the cap (no step
fits). 2033: N_rot 1.517e8, eps 8.118e-6, c_BE 3.913e6, 7.825e9 T/shot; QSVT rail 4.07e9-8.23e9;
62 yr serialized, 15.5 yr on four machines, 13 machines for the 5-yr horizon.

Updated 2026-10-01 for r19 (E20 applied to the chapter; E21 rulings; apply_log/r19_ch10.md): every gauge
rotation at the full fit; the hop's color squish and parity flags held per link (share="link") with the
estimated frame undo; share="draft" kept as the conservative sensitivity; ~0.1 ms per-shot overhead.
2028: step 399,956.8 T = 4.0x the cap, eps_l 2.50e-7, 0.40 s/shot, 13 min. 2033: c_BE 2.760e6,
5.521e9 T/shot (5.5x), eps_l 1.81e-11, ~55 faults at the floor, 43.7 yr serialized, 10.9 yr on four,
9 machines; QSVT 2.86-5.81e9; V=3^3 6.3e9, ~63x cut. Algorithmic register 108-137.

Updated 2026-10-01 for r22 (E26; apply_log/r22_hwp_chapters.md): a Hamming-weight-phasing group of k holds
k - w(k) ancilla (HWP(6) 4, HWP(3) 1; were 7 and 3). 2028 primitive workspace 29-55 (was 31-60), algorithmic
register 106-132 (was 108-137). The ~100 sized budget, the ~180 box and every T count are unchanged.

Updated 2026-10-02 for r23 (apply_log/r23_ch10.md; ruling "we don't want to anywhere assume we have multiple
machines"): one machine, serial wall times only. n_machines, wall_parallel_yr, machines_for_horizon and
horizon_ratio_parallel are retired (were 4, 10.9 yr, 9). New: 5 years hold 2.86e4 shots = two full grid
points; the 43.7-yr grid needs an 8.7x cut in shots x per-shot cost (6.3e8 T/shot at 2.5e5 shots).

Updated 2026-10-02 for r25 (rulings R1, R9; apply_log/r25_ch10.md): the seven R9 exports in both eras, T-depth from
factory.json. 2028: F* 9.9-28 (E29: 9.9-15), depth-limited wall 810 s = ~13.5 min per sector, 27-40 min for 2-3 sectors, 2.2% rms
per sector at sigma^2 <= 1; t_per_shot is now the (lo, hi) band and the scalar is t_per_shot_value. 2033: F* 17-106 (E29: 17-56),
grid 43.7 yr (floor 4.1-25 yr); first result one T row from two reweighted ensembles at 30%, 798-2,873 samples,
51 d - 0.50 yr; kappa2 records (2.5-3.9e5 / 2.2e6 / reweighted 3.9e4-1.3e5); chain-reuse eps_l 1.8e-10.

Updated 2026-10-02 for E29 (Ch. 10 verifier fix 1; author approved open item 7): the T-par low ends sit inside the
priced ancilla budget. 2028 low end 8 link rounds 2.67e4 (was 4 rounds 1.43e4), F* 9.9-15; 2033 low end 12 rounds
9.85e7 (was 6 rounds 5.23e7, ~880 ancilla > 750), F* 17-56, grid floor 7.8-25 yr. Walls and box values unchanged.

Updated 2026-10-05 for the author ruling "Make the QSVT/TPQ route in Ch 10 as well": the 2033 headline is the QSVT/TPQ
thermal state of Ch. 9 (qsvt_thermal): d_beta 30 at delta_beta 1e-4, 8 rounds at amplitude 0.1 -> 17 filter calls,
510 queries at 2-4 c_step, reference / phases / reflections < 1e-4: 2.77e9-5.63e9 T ('2.8-5.6e9'); walls 46-94 min
per shot, first result 0.26-2.8 yr, grid 23-46 yr (4.5-9.2x), eps_l 1.8e-11. The Gibbs sampler is a record (one
sentence). 2028 unchanged; the same filter at V=2^2 is 4.7-9.6e8 T.
"""

import math
import re
from pathlib import Path

import pytest

from estimates import common as c
from estimates import ch10_finite_density as m
from estimates.groups import (GROUPS, magnetic_per_link, electric_per_link, rhop_link,
                              hop_link_cost, staggered_mass_site, fp_printed)


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


def _close(x, y, rel=0.0, abs_=0.0):
    return math.isclose(x, y, rel_tol=rel, abs_tol=abs_)


# --- provenance and validation --------------------------------------------

def test_every_assumption_is_tagged(a):
    import dataclasses
    for f in dataclasses.fields(a):
        assert isinstance(getattr(a, f.name), c.Tagged), f.name


def test_stated_inputs_carry_line_refs(a):
    # the numbers the chapter asserts without derivation must point at a line
    for name in ("c_be_2033_superseded", "n_anc_2028", "n_anc_2033", "c_mix", "beta_h_2033",
                 "n_trotter_2028", "shots_2028"):
        v = getattr(a, name)
        assert v.prov is c.Provenance.STATED, name
        assert v.src.startswith("app07:"), (name, v.src)
    # and the paper-convention inputs are Cited to the papers; the hop and mass rotations use the full fit (r17)
    assert a.rot_synthesis_paper.prov is c.Provenance.CITED
    assert a.hop_synthesis.value == "rus" and a.hop_synthesis.prov is c.Provenance.ASSUMED
    assert a.toffoli_convention.prov is c.Provenance.CITED
    assert a.t_cap_2028.src == "DOE_RFI_2026"
    # R-TOL: the per-shot synthesis budget is the report-wide ruling, not a chapter number
    assert a.eps_syn.lo == c.EPS_SYN and a.eps_syn.src == c.EPS_SYN_SRC


def test_post_init_rejects_nonsense():
    with pytest.raises(ValueError):
        m.Assumptions(c_mix=c.Stated(-1, "x"))
    with pytest.raises(ValueError):
        m.Assumptions(n_anc_2033=c.Stated((750, 250), "x"))
    with pytest.raises(ValueError):
        m.Assumptions(group=c.Stated("SU7", "x"))
    with pytest.raises(ValueError):
        m.Assumptions(toffoli_convention=c.Stated("whatever", "x"))
    with pytest.raises(ValueError):
        m.model(m.Assumptions(), "2040")
    with pytest.raises(ValueError):
        m.Assumptions(filter_success_prob=c.Assumed(5, "x"))           # probability > 1
    with pytest.raises(ValueError):
        m.Assumptions(aa_call_rule=c.Assumed("3k", "x"))
    with pytest.raises(ValueError):
        m.Assumptions(n_anc_2028=c.Stated((200, 100), "x"))            # reversed range
    with pytest.raises(ValueError):
        m.Assumptions(group=c.Stated("SU3", "x"))                      # no compiled primitives in GROUPS
    with pytest.raises(ValueError):
        m.Assumptions(ham=c.Stated("XX", "x"))                         # not a PRIMCOST Hamiltonian
    # r17: the 2028 step no longer fits, so the cap check is gone; the hop readings are validated instead
    with pytest.raises(ValueError):
        m.Assumptions(hop_synthesis=c.Assumed("magic", "x"))
    with pytest.raises(ValueError):
        m.Assumptions(hop_share=c.Assumed("site", "x"))
    with pytest.raises(ValueError):
        m.Assumptions(hop_mcx=c.Assumed("4n", "x"))
    with pytest.raises(ValueError):
        m.Assumptions(n_trotter_2028=c.Stated(1.5, "x"))


def test_ranges_on_2028_inputs_propagate():
    a = m.Assumptions(n_anc_2028=c.Stated((100, 200), "x"), shots_2028=c.Stated((1e3, 3e3), "x"))
    r = m.model(a, "2028")
    assert r.lq == (176, 276)
    assert r.shots == (1e3, 3e3)
    # r25: per sector at the 10-factory baseline, depth-limited at the as-compiled end (D_T 40,488.9 x 10 us)
    w = max(399956.82855117304 * 1e-6, 40488.90205167576 * 1e-5) + 1e-4
    assert _close(r.wall_time_s[0], 1e3 * w, rel=1e-9) and _close(r.wall_time_s[1], 3e3 * w, rel=1e-9)


# --- gauge group: compiled width (R1) ---------------------------------------

def test_link_width_is_the_compiled_value(r28, r33, rcd):
    g = GROUPS["S36x3"]
    assert g.link_qubits == 8                              # app07:91,115,138,161 "8 qubits/link"
    assert int(g.link_qubits_compiled.lo) == 8
    for r in (r28, r33):
        assert r.intermediates["link_qubits"] == 8
        assert r.intermediates["link_qubits_chapter_old"] == 7   # what was printed before R1
        assert r.intermediates["link_qubits_log2"] == 7    # ceil(log2 108) = 7, the dense minimum
    assert rcd.intermediates["link_qubits"] == 9           # Sigma(72x3), app07:115,196


def test_group_arithmetic_quoted_in_prose(r33):
    assert math.ceil(math.log2(216)) == 8                  # dense minimum; compiled is 9
    assert r33.intermediates["sigma360_link_qubits_min"] == 11  # "at least ceil(log2 1080) = 11", app07:127
    assert _close(r33.intermediates["compression_36_vs_72"], 1.125, rel=1e-9)  # app07:115 "~1.13"
    assert _close(r33.intermediates["compression_36_vs_72"], 1.13, rel=0.01)


# --- per-link costs come from the papers through groups.py (R2) ----------------

def test_per_link_costs_are_the_papers_numbers(a, r28, r33):
    g, h = a.group.value, a.ham.value
    i, j = r28.intermediates, r33.intermediates
    e = i["eps_rot"]
    assert i["magnetic_t_per_link"] == magnetic_per_link(g, h, 2, e)
    assert i["electric_t_per_link"] == electric_per_link(g, h, 2, e)
    # E20 (r19 apply): every gauge rotation at the full fit 1.15 log2(1/eps) + 9.2
    assert _close(i["magnetic_t_per_link"], 2.5e3, rel=0.02)      # app07:96 "~2.5e3 T per link at d=2" (2,466.6)
    assert _close(i["electric_t_per_link"], 1.1e4, rel=0.04)      # "~1.1e4 T per link" (2028 eps, 10,601.1)
    assert _close(i["electric_share"], 0.81, rel=0.01)            # app07:103 "~81%"
    assert i["electric_t_per_link"] > i["magnetic_t_per_link"]    # "the electric term, not the plaquette, dominates"
    assert round(i["naive_over_fft"], -1) == 630                 # "~630-650x more" without the FFT (632.8)
    assert round(j["naive_over_fft"], -1) == 650                 # (650.8)
    assert j["magnetic_t_per_link"] == magnetic_per_link(g, h, 3, j["eps_rot"])
    assert _close(j["magnetic_t_per_link"], 5.0e3, rel=0.01)      # "5.0e3 at d=3" (2033 eps, 4,988.5)
    assert _close(j["electric_t_per_link"], 1.4e4, rel=0.02)      # "1.4e4" (2033 eps, 14,241.6)
    # the printed full-fit forms: FFT 532 + 102 x 9.2 = 1470 constant, U_phi 256 x 9.2 = 2355 constant
    pr = GROUPS[g].primitives
    assert round(pr["U_FFT"].t_const + 9.2 * pr["U_FFT"].n_rot) == 1470 and round(1.15 * pr["U_FFT"].n_rot, 1) == 117.3
    assert round(pr["U_phi"].t_const + 9.2 * pr["U_phi"].n_rot) == 2355 and round(1.15 * pr["U_phi"].n_rot, 1) == 294.4
    # the electric term is d-independent at fixed eps; the eras differ only through their own eps
    assert electric_per_link(g, h, 3, e) == i["electric_t_per_link"]


def test_electric_pieces_match_groups():
    # the breakdown's 2 U_FFT + 1 U_phi must be what groups.electric_per_link prices
    g = GROUPS["S36x3"]
    e = 1e-4
    assert _close(2 * g.primitives["U_FFT"].t(e) + g.primitives["U_phi"].t(e),
                  electric_per_link("S36x3", "KS", 2, e), rel=1e-12)


# --- 2028 benchmark box (app07:134-156) ------------------------------------

def test_2028_register_components(r28):
    i = r28.intermediates
    assert i["n_links"] == 8 and i["n_sites"] == 4 and i["n_plaquettes"] == 4
    assert i["gauge_lq"] == 64          # "64 gauge (d L^d . q_G = 8 . 8)"
    assert i["fermion_lq"] == 12        # "+ 12 fermion (1 staggered field, JW on V=2^2)"
    assert i["ancilla_lq"] == (100, 100)  # "+ ~100 ancilla budget" (ruling ch10-ancilla-contents A)
    assert i["lq_system"] == 76
    assert i["lq_total"] == (176, 176)
    assert _close(i["lq_total"][0], 180, rel=0.03)   # "~180"
    assert r28.lq == (176, 176)


def test_2028_ancilla_budget_accounting(r28):
    # app07:93 "primitive workspace 29--55 ..., Hadamard test 1; ... algorithmic register is 106--132 LQ
    # of the ~180 sized" (ruling ch10-ancilla-contents A; r17 adds the hop HWP(6), the C^7X ladder, the draft's
    # 24-qubit Sigma(36x3) squish scratch (su3_diag.tex:103), mass HWP(3); r19 holds the squish through the hop.
    # r22 (E26): a phasing group of k holds k - w(k) ancilla, so HWP(6) is 4 (was 7) and HWP(3) is 1 (was 3); the
    # serial-reuse floor is 24 + the C^7X ladder 5 = 29 (was 24 + HWP(6) 7 = 31) and the no-reuse sum 55 (was 60))
    prim = GROUPS["S36x3"].primitives
    assert (c.hwp_ancilla(6), c.hwp_ancilla(3)) == (4, 1)
    ws = [prim[k].ancilla for k in ("U_FFT", "U_Tr", "U_inv", "U_mul")] + [c.hwp_ancilla(6), 5, 24, c.hwp_ancilla(3)]
    assert sorted(ws) == [1, 2, 4, 4, 5, 7, 8, 24]
    assert sorted(r28.intermediates["ancilla_workspace"].values()) == sorted(ws)
    assert (max(ws), sum(ws)) == (24, 55)
    assert r28.intermediates["ancilla_workspace_serial"] == 24 + max(c.hwp_ancilla(6), 5) == 29
    # the repeat-until-success synthesis ancilla (1, E26) does not move the floor: HWP(6) 4 + 1 = the ladder's 5
    assert r28.intermediates["ancilla_rus"] == c.RUS_ANCILLA == 1
    assert r28.intermediates["ancilla_workspace_serial_with_rus"] == 24 + c.hwp_ancilla(6) + c.RUS_ANCILLA == 29
    lq_sys = r28.intermediates["lq_system"]
    assert r28.intermediates["algorithmic_register"] == (lq_sys + 29 + 1, lq_sys + sum(ws) + 1) == (106, 132)
    assert lq_sys + sum(ws) + 1 <= r28.intermediates["lq_total"][0]   # accounted part fits the budget
    tex = Path(__file__).resolve().parents[3] / "applications" / "app07_finite_density.tex"
    if tex.exists():
        text = tex.read_text()
        assert r"primitive workspace $29$--$55$" in text and r"algorithmic register is $106$--$132$ LQ" in text
        assert r"$31$--$60$" not in text and r"$108$--$137$" not in text


def test_2028_per_step_cost_derived_from_the_papers(r28):
    i = r28.intermediates
    assert _close(i["c_be_gauge_t_per_step"], 104541.03412795652, rel=1e-12)  # "1.0e5 T for the gauge terms" (E20)
    assert _close(i["c_be_gauge_t_per_step"], 8 * i["gauge_t_per_link"], rel=1e-12)
    assert _close(i["c_be_gauge_t_per_step_at_paper_fiducial_eps"], 1.8e5, rel=0.03)   # "1.8e5 at eps=1e-8" (175,100.5)
    assert _close(i["c_be_gauge_t_per_step_papers_convention"], 7.0e4, rel=0.01)       # "the papers' own convention 7.0e4"
    assert _close(i["c_be_gauge_t_per_step_report_convention"], i["c_be_gauge_t_per_step"], rel=1e-12)
    # r19 hop: Sigma(36x3), 1 field, share="link" + undo: 2,592 Toffoli + 120 T + 899 rotations at 20.733 T
    # = 36,902.7 T/link ("3.7e4")
    assert i["hops"] == 8
    assert (i["hop_toffoli_per_link"], i["hop_t_direct_per_link"], i["hop_n_rot_per_link"]) == (2592, 120, 899)
    assert _close(i["hop_t_per_link"], 36902.74155367662, rel=1e-12)
    assert _close(i["hop_t_per_rotation"], c.t_per_rotation(i["eps_rot"], "rus"), rel=1e-12)   # full fit
    assert _close(i["hop_colour_frame_share"], 0.80, rel=0.01)                  # "the color frame is 80%"
    assert _close(i["hop_t_per_step"], 3.0e5, rel=0.02)                         # "the 8 hops are 3.0e5"
    assert _close(i["hop_t_per_link_fp_printed"], 3.2e4, rel=0.01)              # "the draft's own formula 3.2e4"
    assert _close(i["hop_t_per_link_share_draft"], 85118.74155367662, rel=1e-12)   # sensitivity "8.5e4" (r17 value)
    assert _close(i["c_be_t_per_step_share_draft"], 785684.828551173, rel=1e-12)
    assert _close(i["hop_t_per_link_rhop_legacy"], 676.13, rel=1e-4)            # the retired R-HOP price
    # mass: HWP(3) per site, 2 rotations + 1 Toffoli at the full fit ("48 T per site")
    assert i["mass_hwp_k"] == 3 and i["n_rot_mass_per_step"] == 8
    assert _close(i["mass_t_per_site"], 48.465498450893485, rel=1e-12)
    assert _close(i["mass_t_per_step"], 1.9e2, rel=0.03)                        # "1.9e2 for the mass"
    assert _close(i["c_be_t_per_step"], 399956.82855117304, rel=1e-12)         # "4.0e5 T/step"
    assert _close(i["c_be_t_per_step"],
                  i["c_be_gauge_t_per_step"] + i["hop_t_per_step"] + i["mass_t_per_step"], rel=1e-12)
    assert round(i["cap_fraction"], 1) == 4.0 and i["fits_cap"] is False      # "4.0x the cap"
    assert i["t_cap"] == 1e5 and i["n_trotter"] == 1 and i["steps_under_cap"] == 0   # no step fits
    assert i["ramp_original_steps"] == 10
    # R-TOL: the ramp is its own circuit, N_rot = 10 x 10,908, eps_rot = sqrt(1e-2/109,080) = 3.028e-4
    assert i["ramp_original_n_rot"] == 109080
    assert _close(i["ramp_original_eps_rot"], (1e-2 / 109080) ** 0.5, rel=1e-12)
    assert _close(i["ramp_original_t"], 4.207922937551e6, rel=1e-9)     # "a ten-step ramp costs 4.2e6 T"
    assert round(i["ramp_original_t"], -5) == 4.2e6
    assert i["ramp_over_cap"] > 1


def test_2028_hard_ops_band_and_breakdown(r28):
    lo, hi = r28.hard_ops
    assert lo == hi and _close(lo, 399956.82855117304, rel=1e-12)     # one priced step; box "~4.0e5"
    assert _close(r28.breakdown_total(), lo, rel=1e-9)                # gauge + hop + mass, no allowance
    by = {p.name: p for p in r28.breakdown}
    assert set(by) == {"magnetic_U_inv", "magnetic_U_mul", "magnetic_U_Tr", "electric_U_FFT", "electric_U_phi",
                       "hop_colour_squish_and_parity", "hop_colour_rotations", "hop_squish_and_flags",
                       "hop_diagonalizers", "hop_phasing_hwp", "mass_hwp_rotations", "mass_hwp_toffolis"}
    assert by["hop_colour_squish_and_parity"].count == 8 and by["hop_colour_squish_and_parity"].t_each == 7 * 2296
    assert by["hop_colour_squish_and_parity"].status is c.CircuitStatus.SCALING        # undo estimated
    assert by["mass_hwp_rotations"].count == 8 and by["mass_hwp_toffolis"].count == 4
    assert by["magnetic_U_mul"].count == 48 and by["magnetic_U_mul"].t_each == 308     # 6 per link x 8 links
    assert by["magnetic_U_inv"].count == 24 and by["magnetic_U_inv"].t_each == 119
    assert by["magnetic_U_Tr"].count == 4                                              # (d-1)/2 per link
    assert by["electric_U_FFT"].count == 16 and by["electric_U_phi"].count == 8
    for n in ("magnetic_U_inv", "magnetic_U_mul", "magnetic_U_Tr", "electric_U_FFT", "electric_U_phi"):
        assert by[n].status is c.CircuitStatus.COMPILED and by[n].src
    assert all(p.status is not c.CircuitStatus.UNSOURCED for p in r28.breakdown)
    assert c.close(r28.hard_ops, m.PUBLISHED["2028"].hard_ops, 0.10)


def test_2028_shots_wall_and_eps_l(r28):
    i = r28.intermediates
    assert i["shots"] == (2e3, 2e3)                           # "~2e3 (single point)"
    assert i["shot_overhead_s"] == 1e-4                       # ~0.1 ms per shot (report rule)
    assert round(i["wall_per_shot_s"], 2) == 0.40             # "~0.40 s (at 1 us/T-gate)"
    assert round(i["wall_total_s"][0] / 60) == 13             # serial record, 800.1 s
    # r25 (R9): F* = 9.9 at the as-compiled depth, so the box quotes the depth-limited 810.0 s = 13.5 min ("~13.5 min")
    assert r28.shots == (2e3, 2e3) and _close(r28.wall_time_s[0], 809.9780410335, rel=1e-9)
    assert round(r28.wall_time_s[0] / 60, 1) == 13.5 and round(i["wall_per_shot_baseline_s"], 2) == 0.40
    # R3: 0.1 expected faults in 4.0e5 ops -> eps_l <~ 2.5e-7; the RFI floor suffices
    assert i["faults_per_shot_budget"] == 0.1
    assert _close(i["eps_l_required"], 2.50027e-7, rel=1e-5)
    assert r28.epsilon_l == (i["eps_l_required"], i["eps_l_required"])
    assert i["eps_l_required"] > 1e-8


def test_f2_2028_gibbs_sample_priced(a, r28):
    # app07 2028 paragraph: "a Gibbs sample on V=2^2 ... 32 jumps ... 1.3e12 T (1.2e11--1.3e13), 1e7 times the cap"
    i = r28.intermediates
    assert i["gibbs_2x2_n_jumps"] == 32
    lo, mid, hi = i["gibbs_2x2_t_sample"]
    assert (round(lo, -10), round(mid, -11), round(hi, -12)) == (1.2e11, 1.3e12, 1.3e13)
    assert round(i["gibbs_2x2_over_cap"], -6) == 1.3e7
    assert a.gibbs_beta_h_2028.prov is c.Provenance.ASSUMED


def test_f2_sampler_step_derivation(a, r33):
    """Referee F2 record: one Chen-Kastoryano-Gilyen step on the chapter's Trotter step. Since the 2026-10-05 ruling
    the chapter keeps one sentence: '~8e13--4e14 T per shot'."""
    i = r33.intermediates
    js = m.gibbs_jump_set(3, 2, 3)
    assert js["families"] == dict(plaquette=48, hop=144, baryon=48, field_mixing=48) and js["n_jumps"] == 288
    g = i["gibbs_sampler"]["central"]
    assert g["sampler_steps"] == 20 * 288 == 5760
    lu = math.log(1 / g["eps_u"])
    assert _close(g["eps_u"], 1e-2 / 5760, rel=1e-12)
    assert _close(g["h_beta"], 8 * math.sqrt(lu) + 2 * lu / math.pi, rel=1e-12) and round(g["h_beta"]) == 38
    assert round(g["control_factor"], 2) == 1.16 and g["aux_share"] < 1e-4
    assert _close(i["gibbs_t_per_shot"], 8.390685037940064e13, rel=1e-9)
    lo, hi = i["gibbs_t_per_shot_range"]
    assert (round(lo, -13), round(hi, -14)) == (8e13, 4e14)                 # "~8e13--4e14 T per shot"
    olo, ohi = i["gibbs_over_qsvt"]
    assert (round(olo, -3), round(ohi, -4)) == (1.5e4, 1.6e5)               # "four to five orders of magnitude"
    assert 1e4 < olo and ohi < 1e6
    assert i["n_gibbs_steps"] == a.c_mix.lo * a.beta_h_2033.lo


def test_2028_gibbs_lower_bound_is_carried_not_derived(a, r28):
    # app07:103 "Gibbs sampling on a V=2^2 SU(3) system needs >~14 sampler steps at 4.0e5 T each": Stated (item 10),
    # restated from the earlier ">~1e6 T per shot" at the then 72,515.5 T/step (1e6 / 72,515.5 = 13.8)
    i = r28.intermediates
    assert a.gibbs_2x2_min_steps.prov is c.Provenance.STATED
    assert i["gibbs_2x2_min_steps"] == 14 and round(1e6 / 72515.51564015418) == 14
    assert i["gibbs_2x2_t_per_shot_superseded"] == 1e6
    assert _close(i["gibbs_2x2_t_lower_bound"], 14 * i["c_be_t_per_step"], rel=1e-12)


# --- 2033 utility box (app07:158-201) --------------------------------------

def test_2033_register_components(r33):
    i = r33.intermediates
    assert i["n_links"] == 24 and i["n_sites"] == 8 and i["n_plaquettes"] == 24
    assert i["gauge_lq"] == 192        # "192 gauge"
    assert i["fermion_lq"] == 72       # "72 fermion"
    assert i["ancilla_lq"] == (250, 750)
    assert i["lq_system"] == 264       # sqrt(2^264) in the QSVT paragraph
    assert i["lq_total"] == (514, 1014)  # "= 514--1014"
    assert r33.lq == (514, 1014)
    assert c.close(r33.lq, (500, 1000), 0.05)   # "~500-1000"


def test_2033_hard_ops_product(r33):
    """Author ruling 2026-10-05: the QSVT/TPQ thermal state, Ch. 9's construction."""
    i = r33.intermediates
    assert i["beta_h"] == 1e2                      # "beta||H|| = 1e2"
    assert _close(i["c_be_t_per_step"], 2760459.402694799, rel=1e-12)  # "2.8e6 T/step" at the 2e3-step tolerance
    assert _close(i["t_hsim_lower_bound"], 5.520918805e9, rel=1e-9)   # record (pre-F2 '>~5.5e9')
    assert i["thermal_route"] == "QSVT/TPQ" and i["hard_ops_is_conditional"] is True
    assert i["delta_beta"] == 1e-4 and _close(i["d_beta_exact"], math.sqrt(1e2 * math.log(1e4)), rel=1e-12)
    assert i["d_beta"] == 30                                            # "d_beta = 30"
    assert _close(i["qsvt_amplitude"], 0.1, rel=1e-12)                 # "a = 0.1"
    assert i["qsvt_aa_rounds"] == 8 == math.ceil(math.pi / (4 * math.asin(0.1)))   # "k = 8"
    assert i["qsvt_filter_calls"] == 17 and i["qsvt_queries"] == 510   # "17 times, or 510 queries"
    assert i["qsvt_step_equivalents"] == (1020, 2040)                  # "1.0--2.0e3 step equivalents"
    elo, ehi = i["qsvt_eps_rot"]
    assert (round(elo, 6), round(ehi, 7)) == (1.1e-5, 8.0e-6)          # "1.1e-5--8.0e-6"
    lo, hi = i["t_per_shot_range"]
    assert _close(lo, 2.7726e9, rel=1e-4) and _close(hi, 5.6340e9, rel=1e-4)
    assert (round(lo, -8), round(hi, -8)) == (2.8e9, 5.6e9)            # "2.8--5.6e9"
    assert r33.hard_ops == (lo, hi) and c.close(r33.hard_ops, m.PUBLISHED["2033"].hard_ops, 0.10)
    assert all(x < 1e-4 for x in i["qsvt_overhead_share"])            # "less than 1e-4 of the shot"
    assert i["qsvt"]["low"]["reflection_toffolis"] == 1014 and i["qsvt_t_readout"] == 0.0  # "at most 1014 Toffolis each"
    rlo, rhi = i["reference_fraction"]
    assert (round(rlo, 1), round(rhi, 1)) == (2.8, 5.6)                # "2.8--5.6 times the 1e9 reference"
    assert _close(r33.breakdown_total(), lo, rel=1e-9)                 # breakdown = the low end
    # the pre-ruling rail (R = 10 x 2 queries x d 26.3 at eps 1e-3, no reference / phases / reflections) is a record
    olo, ohi = i["qsvt_old_t_per_shot"]
    assert _close(olo, 2.8595e9, rel=1e-3) and _close(ohi, 5.8108e9, rel=1e-3)
    assert i["qsvt_t_per_shot_stated"] == (3e9, 6e9) and i["d_qsp_stated"] == 26 and i["R_warm_stated"] == 10


def test_2033_c_be_consistent_with_compiled_floor(r33):
    i = r33.intermediates
    assert _close(i["c_be_gauge_t_per_step"], 4.6e5, rel=0.01)       # app07:96 "4.6e5 at V=2^3" (E20, 461,524)
    assert i["hops"] == 72
    # r19 hop: Sigma(36x3), 3 fields, share="link" + undo: 2,604 Toffoli + 360 T + 2,693 rotations at 28.647 T
    # = 95,734 T/link ("9.6e4")
    assert (i["hop_toffoli_per_link"], i["hop_t_direct_per_link"], i["hop_n_rot_per_link"]) == (2604, 360, 2693)
    assert _close(i["hop_t_per_link"], 95734.44975033074, rel=1e-12)
    assert _close(i["hop_t_per_link"], 9.6e4, rel=0.01)
    assert _close(i["hop_colour_frame_share"], 0.76, rel=0.01)                  # "76%"
    assert _close(i["hop_t_per_step"], 2.3e6, rel=0.01)                         # "the 72 hops ... 2.3e6"
    assert round(100 * i["hop_share_of_c_be"]) == 83                           # "(83%)"
    assert _close(i["hop_t_per_link_fp_printed"], 5.4e4, rel=0.01)              # draft formula, 2 C^G + 3 C^hop
    assert _close(i["hop_t_per_link_fp_printed"],
                  2 * fp_printed("S36x3", "C_G", i["eps_rot"]) + 3 * fp_printed("S36x3", "C_hop_text", i["eps_rot"]),
                  rel=1e-12)
    # sensitivity share="draft" (r17 value): "1.5e5" per link, "an 8.0e9-T shot at 2033"
    assert _close(i["hop_t_per_link_share_draft"], 148038.44975033074, rel=1e-12)
    assert _close(i["hop_t_per_link_share_draft"], 1.5e5, rel=0.02)
    assert round(2e3 * i["c_be_share_draft"] / 1e9, 1) == 8.0
    assert i["mass_hwp_k"] == 9 and i["n_rot_mass_per_step"] == 32
    assert _close(i["mass_t_per_site"], 164, rel=0.01)                          # "164 T at 3"
    assert _close(i["mass_t_per_step"], 1.3e3, rel=0.01)                        # "1.3e3 for the mass"
    assert _close(i["c_be_gauge_t_per_step"] + i["hop_t_per_step"] + i["mass_t_per_step"], i["c_be_t_per_step"],
                  rel=1e-12)
    assert i["c_be_stated_superseded"] == 5e5
    assert round(i["c_be_over_superseded"], 1) == 5.5                          # the derived value replaces 5e5
    by = {p.name: p for p in r33.breakdown}
    assert all(p.status is not c.CircuitStatus.UNSOURCED for p in r33.breakdown)
    n_tr = i["qsvt_step_equivalents"][0]                             # the QSVT shot's step equivalents (low end)
    assert _close(by["hop_colour_rotations"].count, n_tr * 24, rel=1e-12)
    assert _close(by["mass_hwp_rotations"].count, n_tr * 8 * 4, rel=1e-12)
    assert _close(by["magnetic_U_mul"].count, n_tr * 24 * 12, rel=1e-12)   # 6(d-1) = 12 per link at d=3
    assert _close(by["electric_U_FFT"].count, n_tr * 24 * 2, rel=1e-12)


def test_2033_shots_from_eq_nshot(r33):
    # referee F1: Eq. Nshot_fd = sum over points of V1/delta^2, V1 the delta-method variance of k4/k2
    i = r33.intermediates
    assert i["n_grid_points"] == 25                        # "5x5 points"
    assert _close(i["shots_total"], 257603.4856, rel=1e-8)  # "2.6e5 over the 5x5 grid"
    assert round(i["shots_total"], -4) == 2.6e5
    assert (round(i["shots_per_point_min"], -2), round(i["shots_per_point_max"], -3)) == (3.7e3, 6.1e4)
    assert i["shots_per_point_by_mu_T"][(100, 600)] == i["shots_per_point_max"]
    assert _close(i["shots_per_point"] * 25, i["shots_total"], rel=1e-12)
    assert c.close(r33.shots, (i["shots_total"],) * 2, 1e-12)


def test_f1_delta_method_variance():
    # Skellam at mu = 0 has every even cumulant k and odd zero: V1 = 36 + 72 k + 24 k^2 exactly
    for k in (0.01, 0.1, 1.0):
        v = m.ratio_variance_per_sample("baryon", k, 0.0)
        assert _close(v["r"], 1.0, rel=1e-9) and _close(v["v1"], 36 + 72 * k + 24 * k * k, rel=1e-6), (k, v)
    # a near-Gaussian (large-k Skellam) gives Var(k4/k2) ~ 24 k^2 per sample, not kappa_8 / kappa_4^2 = 1/k
    v = m.ratio_variance_per_sample("baryon", 30.0, 0.0)
    assert _close(v["v1"] / (24 * 30.0 ** 2), 1.0, rel=0.15)
    # the quark-like toy has chi4/chi2 = 1/9 at mu = 0
    assert _close(m.ratio_variance_per_sample("quark", 0.05, 0.0)["r"], 1 / 9, rel=1e-9)
    # reweighting from the target itself changes nothing
    d = m.ratio_variance_per_sample("baryon", 0.05, 1.0)
    w = m.ratio_variance_per_sample("baryon", 0.05, 1.0, (1.0, 1.0))
    assert _close(d["v1"], w["v1"], rel=1e-9)


def test_f1_kappa2_records(r33):
    i = r33.intermediates
    assert round(i["grid_shots_k2_mid"], -4) == 9.1e5                    # "9.1e5 at kappa_2 = 0.05"
    assert round(i["grid_shots_k2_hi"], -5) == 2.2e6                     # "2.2e6 at 0.1"
    lo, hi = i["grid_wall_k2_mid_yr"]
    assert (round(lo, -1), round(hi, -1)) == (80, 160)                   # "(80--160 yr)"
    lo, hi = i["grid_wall_k2_hi_yr"]
    assert (round(lo, -2), round(hi, -2)) == (200, 400)                  # "(200--400 yr)"
    assert round(i["grid_shots_quark_abs"], -2) == 4.8e3                 # "4.8e3 at delta = 0.1"
    assert round(i["grid_shots_quark_rel"], -4) == 3.9e5                 # "3.9e5 for 10% relative"
    assert round(100 * i["reweighting_saving_k2_lo"]) == 11               # "saves 11% at kappa_2 = 0.01"
    assert i["reweighting_ratio_k2_mid"] > 1                              # "and loses at 0.05"
    assert i["hard_ops_is_lower_bound"] is False and i["hsim_bound_is_lower_bound"] is True


def test_2033_wall_time_chain(r33):
    i = r33.intermediates
    assert i["shot_overhead_s"] == 1e-4
    w = i["wall_per_shot_s"]
    assert _close(w[0], i["t_per_shot_range"][0] * 1e-6 + 1e-4, rel=1e-12)   # serial at both ends (F* >= 17.5)
    assert _close(w[1], i["t_per_shot_range"][1] * 1e-6 + 1e-4, rel=1e-12)
    assert (round(w[0], -2), round(w[1], -2)) == (2.8e3, 5.6e3)          # "2.8--5.6e3 s"
    assert [round(x) for x in i["wall_per_shot_min"]] == [46, 94]        # "(46--94 min)"
    lo, hi = i["wall_campaign_yr"]
    assert (round(lo), round(hi)) == (23, 46)                            # "23--46 years of serial running"
    assert _close(r33.wall_time_s[0], i["shots_total"] * w[0], rel=1e-12)
    plo, phi = i["wall_single_point_yr"]
    assert (round(plo, 2), round(phi, 1)) == (0.91, 1.8)                 # "0.91--1.8 yr for an average point"
    flo, fhi = i["floor_wall_campaign_yr"]
    assert (round(flo, 1), round(fhi)) == (4.0, 26)                      # "depth limit is 4.0--26 yr"
    assert i["campaign_horizon_yr"] == 5
    hlo, hhi = i["horizon_reduction_needed"]
    assert (round(hlo, 1), round(hhi, 1)) == (4.5, 9.2)                  # "falls by 4.5--9.2"
    assert _close(i["t_per_shot_to_fit_horizon"], 6.1e8, rel=0.01)       # "6.1e8-T shot at the same shot count"
    assert i["fr_fits_horizon"] is True                                  # "Five years hold the first result"


def test_r23_no_machine_count(r33):
    # ruling 2026-10-02: no multi-machine assumption anywhere in the model or the chapter
    assert not hasattr(m.Assumptions(), "n_machines")
    for k in ("wall_parallel_yr", "machines_for_horizon", "horizon_ratio_parallel", "n_machines"):
        assert k not in r33.intermediates
    tex = (Path(m.__file__).resolve().parents[2] / "applications" / "app07_finite_density.tex").read_text()
    assert not re.search(r"shot-parallel|parallel machines|across machines|machine-years|fold parallelism|on four", tex)


def test_2033_accuracy_targets_are_carried(r33):
    # app07:167 "10% on chi4/chi2 (sets the shots); 5% on p is not separately costed"
    # (ruling ch10-pressure-5pct a): only the chi4 estimator at eps=0.1 is costed
    i = r33.intermediates
    assert i["target_accuracy_pressure"] == 0.05
    assert i["target_accuracy_chi4_over_chi2"] == 0.10
    assert i["eps_stat_chi4"] == 0.1


def test_2033_logical_error_row(r33):
    # R3: 0.1 expected faults per shot -> box "<~ 1.8e-11 ... at the upper end; 3.6e-11 at the lower", "280--560x below"
    i = r33.intermediates
    assert i["faults_per_shot_budget"] == 0.1
    lo, hi = i["eps_l_required"]
    assert (round(lo, 12), round(hi, 12)) == (1.8e-11, 3.6e-11)
    rlo, rhi = i["rfi_over_eps_l"]
    assert (round(rlo, -1), round(rhi, -1)) == (280, 560)
    flo, fhi = i["faults_at_rfi_floor"]
    assert (round(flo), round(fhi)) == (28, 56)                         # "would accrue 28--56 faults"
    assert r33.epsilon_l == (lo, hi)


def test_2033_qsvt_matches_ch09_construction(a, r33):
    """The same construction as Ch. 9 (ch09_chiral_gauge.py): degree from the tolerance with Python round(), 2k + 1
    filter calls, k from the amplitude; each end at its own R-TOL eps from its own N_rot."""
    from estimates import ch09_chiral_gauge as ch9
    a9 = ch9.Assumptions()
    assert a.delta_beta.lo == float(a9.delta_beta.value) and a.aa_call_rule.value == a9.aa_call_rule.value
    assert int(a9.aa_rounds.value) == r33.intermediates["qsvt_aa_rounds"]       # 8 in both
    i = r33.intermediates
    for end in ("low", "high"):
        q = i["qsvt"][end]
        assert q["n_rot"] == q["step_equivalents"] * 75872 + 17 * 264 + (17 * 31 + 8)
        assert q["eps_rot"] == c.eps_rot_for(q["n_rot"])
        assert _close(q["t_shot"], q["t_filter"] + q["t_reference"] + q["t_signal"] + q["t_reflections"], rel=1e-12)
    assert i["R_unstructured_log2"] == 132                              # "up to sqrt(2^264) rounds"
    assert i["d_naive_linear"] == 1e2


def test_qsvt_rounds_follow_the_success_probability():
    r = m.model(m.Assumptions(filter_success_prob=c.Assumed(1e-4, "x")), "2033")
    i = r.intermediates
    assert i["qsvt_aa_rounds"] == math.ceil(math.pi / (4 * math.asin(1e-2))) == 79
    assert i["qsvt_filter_calls"] == 159
    assert 159 / 17 * 2.7726e9 < r.hard_ops[0] < 1.1 * 159 / 17 * 2.7726e9


def test_2033_v3cubed_promotion(r33):
    # open questions: "At V=3^3 the same filter (d_beta = 56, a 9.6e6-T step) is 1.8--3.7e10 T ... on 1141--1641 LQ"
    i = r33.intermediates
    assert i["v3cubed_lq"] == (1141, 1641)
    assert i["v3cubed_beta_h"] == 337.5 and i["v3cubed_d_beta"] == 56
    assert round(i["v3cubed_c_step"], -5) == 9.6e6
    lo, hi = i["v3cubed_t_per_shot"]
    assert (round(lo, -9), round(hi, -9)) == (1.8e10, 3.7e10)
    # records
    assert i["v3cubed_steps"] == 675 and i["v3cubed_n_jumps"] == 972
    assert i["v3cubed_eps_rot"] == c.eps_rot_for(i["v3cubed_n_rot_per_shot"])
    assert round(i["v3cubed_gibbs_t_with_10x_cut"], -14) == 3e14


# --- post-2033 stretch (app07:196) ------------------------------------------

def test_codesign_register(rcd):
    i = rcd.intermediates
    assert i["link_qubits"] == 9 and i["n_links"] == 192 and i["n_sites"] == 64
    assert i["gauge_lq"] == 1728       # "~1728 gauge at 9 qubits/link"
    assert i["fermion_lq"] == 576      # "~576 fermion"
    assert i["ancilla_lq"] == (250, 750)  # "the 250--750 BE/bath ancilla of the 2033 route"
    assert i["lq_total"] == (2554, 3054)
    assert rcd.lq == (2554, 3054)
    assert _close(i["lq_total"][0], 2600, rel=0.02) and _close(i["lq_total"][1], 3100, rel=0.02)  # "~2600--3100 LQ"
    assert i["link_qubits_chapter_old"] == 8          # ceil(log2 216), what was printed before R1
    lo, hi = i["ratio_to_2033_anchor"]
    assert _close(lo, 2.554, rel=1e-12) and _close(hi, 3.054, rel=1e-12)
    assert round(lo, 1) == 2.6 and round(hi, 1) == 3.1   # "exceeds 2033 register anchor by ~2.6--3.1x"
    assert rcd.hard_ops == (0.0, 0.0) and rcd.breakdown == ()
    assert m.PUBLISHED["codesign"].hard_ops is None


# --- resources rows and requirements summary --------------------------------

def test_instance_rows(a, r28, r33, rcd):
    rows28 = m.INSTANCE_ROWS(a, "2028", r28)
    assert len(rows28) == 1 and rows28[0][1] == (176, 176) and "t_bound" not in rows28[0][3]
    rows33 = m.INSTANCE_ROWS(a, "2033", r33)
    assert len(rows33) == 1 and rows33[0][3]["conditional"] is True and rows33[0][2] == r33.hard_ops
    rowscd = m.INSTANCE_ROWS(a, "codesign", rcd)
    assert rowscd[0][3]["lq_only"] is True


def test_requirements_summary_row(r33):
    # requirements: "2.6e5 over the 25 grid points"; "2.8--5.6e3 s (46--94 min) per QSVT shot"
    i = r33.intermediates
    assert round(i["n_grid_points"] * i["shots_per_point"], -4) == 2.6e5
    assert [round(x, -2) for x in i["wall_per_shot_s"]] == [2.8e3, 5.6e3]


# --- r17: the hop and mass through the shared functions -----------------------

def test_hop_and_mass_are_the_shared_functions(a, r28, r33):
    for r, n_stag in ((r28, 1), (r33, 3)):
        e = r.intermediates["eps_rot"]
        want = hop_link_cost("S36x3", n_stag, eps=e)            # defaults = the chapter's readings (r19: link + undo)
        assert r.intermediates["hop_t_per_link"] == want["t"]
        assert r.intermediates["hop_t_per_step"] == r.intermediates["n_links"] * want["t"]
        assert r.intermediates["mass_t_per_site"] == staggered_mass_site(3, n_stag, eps=e)["t"]
        # the retired R-HOP is carried as a comparison only, at the chapter's old convention
        assert r.intermediates["hop_t_per_link_rhop_legacy"] == rhop_link("S36x3", n_stag, 3, e, synthesis="rus-slope")["t"]


def test_published_matches_the_boxes_as_printed():
    assert m.PUBLISHED["2028"].lq == (180, 180)
    assert m.PUBLISHED["2033"].lq == (500, 1000)
    assert m.PUBLISHED["2028"].hard_ops == (4.0e5, 4.0e5) and m.PUBLISHED["2033"].hard_ops == (2.8e9, 5.6e9)
    assert m.PUBLISHED["codesign"].lq == (2554, 3054)


# --- ruling R-TOL: each circuit's tolerance from its own rotation count ---------

@pytest.mark.parametrize("era,n_step,n_shot,eps,t_rot", [
    ("2028", 10908, 10908, 9.575e-4, 20.7327),     # every rotation at the full fit (E20)
    ("2033", 75872, 1.51744e8, 8.118e-6, 28.6470),
])
def test_rtol_rotation_count_and_tolerance(era, n_step, n_shot, eps, t_rot):
    i = m.model(m.Assumptions(), era).intermediates
    assert i["n_rot_per_step"] == n_step and i["n_rot_per_shot"] == n_shot
    assert i["eps_rot"] == c.eps_rot_for(n_shot)                      # the shared function, eps_syn = 1e-2
    assert _close(i["eps_rot"], eps, rel=1e-3) and _close(i["t_per_rotation"], t_rot, rel=1e-4)
    assert _close(i["synthesis_error_per_shot"], 1e-2, rel=1e-12)     # N eps^2 = eps_syn


def test_rtol_rotation_count_matches_the_breakdown():
    """Perturbation: dT/d log2(1/eps) must be 1.15 x N_rot, so no synthesized rotation is missing."""
    a = m.Assumptions()
    for era, (d, L, ns) in (("2028", (2, 2, 1)), ("2033", (3, 2, 3))):
        i = m.model(a, era).intermediates
        n_steps = i["n_trotter"] if era == "2028" else i["n_gibbs_steps"]
        e1, e2 = 1e-3, 1e-5
        dT = n_steps * (m._t_per_step(a, d, L, ns, e2) - m._t_per_step(a, d, L, ns, e1))
        assert _close(dT, 1.15 * i["n_rot_per_shot"] * math.log2(e1 / e2), rel=1e-9)


def test_r17_2028_no_step_fits_the_cap(r28):
    i = r28.intermediates
    assert i["steps_under_cap"] == 0 and i["n_trotter"] == 1        # one step priced, 4.0x the cap


def test_r19_hop_reading_is_e21(a, r28, r33):
    """E21 (1): squish and flags held per link with the frame undo kept; share="draft" is the sensitivity."""
    assert a.hop_share.value == "link" and a.hop_undo.value is True
    for r, n_stag in ((r28, 1), (r33, 3)):
        e = r.intermediates["eps_rot"]
        assert r.intermediates["hop_t_per_link"] == hop_link_cost("S36x3", n_stag, eps=e, share="link", undo=True)["t"]
        assert r.intermediates["hop_t_per_link_share_draft"] == hop_link_cost("S36x3", n_stag, eps=e, share="draft")["t"]
        assert r.intermediates["hop_t_per_link_no_undo"] == hop_link_cost("S36x3", n_stag, eps=e, undo=False)["t"]
        assert r.intermediates["hop_t_per_link_share_draft"] > r.intermediates["hop_t_per_link"]


# --- r25: rulings R1 (two tiers), R9 (depth and factory exports), 2028 shot precision -------------

def test_r25_2028_exports_and_sectors(a, r28):
    i = r28.intermediates
    assert i["t_per_shot"] == (399956.82855117304,) * 2
    assert i["t_depth_per_shot"] == a.t_depth_2028.value == (26719.289742577264, 40488.90205167576)
    lo, hi = i["f_star"]
    assert round(lo, 1) == 9.9 and round(hi) == 15                       # "9.9--15 factories"
    assert i["baseline_ok"] == (False, True)                              # depth-limited at the as-compiled end
    assert i["wall_first_result_s"] is None                               # one-tier 2028 benchmark
    assert i["n_sectors"] == (2, 3)
    w = i["wall_campaign_s"]
    assert _close(w[0], 2 * 809.9780410335, rel=1e-9) and _close(w[1], 3 * 809.9780410335, rel=1e-9)
    assert (round(w[0] / 60), math.floor(w[1] / 60)) == (27, 40)           # "27--40 min for 2--3 sectors"
    assert i["eps_per_sector"][0] == pytest.approx(0.02236, rel=1e-3)     # "2.2% rms at sigma^2 <= 1"
    assert round(100 * i["eps_per_sector"][0], 1) == 2.2
    assert i["factories_for_1yr"][0] < 1e-3                               # never a long run


def test_r25_2033_exports(a, r33):
    i = r33.intermediates
    # the factory.json band (2e3-step circuit) scaled by each end's own T ratio; F* 17.5-56 at either end
    assert i["t_depth_per_shot"] == i["t_depth_per_shot_scaled"]
    s_lo, s_hi = i["depth_scale"]
    assert _close(i["t_depth_per_shot"][0], a.t_depth_2033.value[0] * s_lo, rel=1e-12)
    assert _close(i["t_depth_per_shot"][1], a.t_depth_2033.value[1] * s_hi, rel=1e-12)
    flo, fhi = i["f_star_matched"]
    assert round(flo) == 17 and round(fhi) == 56                          # "17--56 factories"
    xlo, xhi = i["f_star"]                                                # cross-paired export (record)
    assert round(xlo, 1) == 8.6 and round(xhi) == 114
    assert i["wall_campaign_s"] == r33.wall_time_s
    assert i["fits_1yr"] == (False, False)
    # matched ends (what the chapter quotes): the 10-factory baseline holds at both ends
    assert i["baseline_ok_matched"] == (True, True) and i["baseline_ok"] == (False, True)
    assert i["fits_1yr_matched"] == (False, False)


def test_r25_2033_first_result(a, r33):
    i = r33.intermediates
    assert i["eps_first_result"] == 0.30
    assert i["fr_samples"] == (2926.0, 15861.0)                           # "2.9e3 ... and 1.6e4"
    by = i["fr_by_toy"]
    assert by[("baryon", 0.01)]["best"] == by[("baryon", 0.01)]["reweighted"]
    assert by[("baryon", 0.1)]["best"] == by[("baryon", 0.1)]["direct"]
    y_lo, y_hi = i["wall_first_result_yr"]
    assert (round(y_lo, 2), round(y_hi, 1)) == (0.26, 2.8)               # "0.26--2.8 yr on one machine"
    assert _close(i["wall_first_result_s"][0], 2926 * i["wall_per_shot_s"][0], rel=1e-12)
    assert _close(i["wall_first_result_s"][1], 15861 * i["wall_per_shot_s"][1], rel=1e-12)
    assert i["wall_first_result_yr"][1] < 5


def test_r25_tex_prints_the_new_numbers():
    tex = (Path(m.__file__).resolve().parents[2] / "applications" / "app07_finite_density.tex").read_text()
    for s in (r"$\sim 13.5\,$min per sector", r"$27$--$40\,$min", r"$2.2\%$", r"$2.9\times 10^{3}$--$1.6\times 10^{4}$",
              r"$9.9$--$15$", r"$2.7$--$4.0\times 10^{4}$", r"$9.1\times 10^{5}$",
              r"$2.2\times 10^{6}$", r"$2.6\times 10^{5}$",
              r"$3.7\times 10^{3}$", r"$6.1\times 10^{4}$", r"$6.1\times 10^{8}$",
              r"$\gtrsim 4.0\times 10^{5}$", r"$4.2\times 10^{6}$", r"$3.2\%$",
              # author ruling 2026-10-05: the QSVT/TPQ headline
              r"$2.8$--$5.6\times 10^{9}$", r"$d_\beta=30$", r"$2m+1=17$", r"$510$ block-encoding queries",
              r"$0.26$--$2.8$~yr", r"$23$--$46$", r"$0.91$--$1.8$~yr", r"$4.0$--$26$~yr", r"$80$--$160$~yr",
              r"$200$--$400$~yr", r"$4.5$--$9.2", r"$17$--$56$", r"$46$--$94$~min", r"\lesssim 1.8\times 10^{-11}$",
              r"$3.6\times 10^{-11}$", r"$280$--$560", r"$28$--$56$ faults", r"$4.7$--$9.6\times 10^{8}$ T",
              r"$1.8$--$3.7\times 10^{10}$", r"$d_\beta=56$",
              r"$9.6\times 10^{6}$-T step", r"$\sim 8\times 10^{13}$--$4\times 10^{14}$ T per shot",
              r"\cite{arxiv_2311_09207,arxiv_2405_20322}", r"lattice cutoff"):
        assert s in tex, s
    # retired numbers: the referee round and the Gibbs-sampler headline (2026-10-04 to 2026-10-05)
    for old in (r"\kappa_8/\kappa_4^2\sim", r"$51$~days", r"$8.7\times$", r"$\sim 44$", r"$3.9\times 10^{4}$--$1.3\times 10^{5}$",
                r"$2.5\times 10^{5}$", r"originally scoped", r"draft does not count", r"rung", r"RFI floor",
                r"$\gtrsim 5.5\times 10^{9}$", r"$\gtrsim 14$", r"$\sim 390$~yr", r"$0.51$--$2.8$~yr",
                r"$8.4\times 10^{13}$--$4.4\times 10^{14}$", r"$7.6\times 10^{12}$", r"$1.3\times 10^{12}$",
                r"$7.8\times 10^{3}$--$4.2\times 10^{4}$~yr", r"$1.2$--$3.9\times 10^{5}$~yr", r"$1.4\times 10^{5}$",
                r"$5{,}760$ sampler steps", r"h_\beta", r"BE/bath", r"warm-start", r"$1.2\times 10^{-15}$",
                r"Gibbs-prep shot"):
        assert old not in tex, old


def test_r25_chain_reuse_record(r33):
    # Gibbs-only lever, a record since 2026-10-05 (no longer in the chapter)
    i = r33.intermediates
    assert i["chain_reuse_gain"] == (3, 30) and i["chain_reuse_steps_per_tau"] == 200
    assert _close(i["chain_reuse_eps_l_required"], 1.2e-14, rel=0.01)
    assert _close(i["chain_reuse_eps_l_required_hsim_superseded"], 1.8e-10, rel=0.01)
    assert _close(i["chain_reuse_eps_l_required"], 10 * i["gibbs_eps_l_required"], rel=1e-9)


def test_r25_depth_bands_and_t_per_shot_scalar(r28, r33):
    t28, d28 = r28.intermediates["t_per_shot"], r28.intermediates["t_depth_per_shot"]
    assert (round(d28[0], -2), round(d28[1], -2)) == (2.67e4, 4.05e4)
    d33 = r33.intermediates["t_depth_per_shot"]
    assert _close(d33[0] / r33.intermediates["depth_scale"][0], 98454105.32846098, rel=1e-12)
    assert t28[0] == r28.intermediates["t_per_shot_value"]
    assert r33.intermediates["t_per_shot"] == r33.intermediates["t_per_shot_range"]
    assert r33.intermediates["t_per_shot_value"] == r33.intermediates["t_per_shot_range"][0]


def test_e29_tex_wording_fixes():
    tex = (Path(m.__file__).resolve().parents[2] / "applications" / "app07_finite_density.tex").read_text()
    for s in (r"$n_B$ for $\mu_B\ge 150$~MeV and $\chi_4/\chi_2$ at all five $\mu_B$",
              r"$\Delta p=p(T,\mu_B)-p(T,0)$ from the same samples",
              r"($2.2\%$ of $\lambda$ per sector, absolute)", r"(depth-limited at the as-compiled end)"):
        assert s in tex, s
    for old in (r"$1.4$--$4.0\times 10^{4}$", r"$17$--$106$", r"$4.1$--$25$~yr", r"$9.9$--$28$"):
        assert old not in tex, old


def test_e29_chain_reuse_ratio_to_floor(r33):
    assert round(1e-8 / r33.intermediates["chain_reuse_eps_l_required"], -4) == 8.4e5   # record


def test_utility_box_heavy_ion_research_base():
    # G5 open item 4 (Claude's decision, 2026-10-04): RHIC operations end for EIC construction; the base is the
    # NP Heavy Ion research line ($47.454M/yr FY24-25 enacted, DOE_NP_FY26_budget); 2% disjoint from Ch. 5's 5%
    p = Path(__file__).resolve().parents[3] / "applications" / "app07_finite_density.tex"
    if not p.exists():
        pytest.skip("chapter .tex not on disk")
    box = p.read_text().split("begin{utilitybox}")[1].split("end{utilitybox}")[0]
    assert round(0.02 * 47.454, 2) == 0.95
    assert "$\\sim\\$0.95$M/yr, for a future production EoS" in box and "($\\sim\\$0.95$M/yr)" in box
    assert "$\\$47.5$M/yr DOE NP Heavy Ion research line" in box
    assert "187" not in box and "3.7$M" not in box
