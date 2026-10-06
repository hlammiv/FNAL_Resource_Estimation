"""Shared conventions, and the contract every chapter model must satisfy."""

import math

import pytest

import estimates
from estimates import common as c
from estimates.groups import GROUPS, conflicts


# --- synthesis --------------------------------------------------------------

def test_rus_matches_bocharov_fit():
    # BRS: 1.15 log2(1/eps) + 9.2. At 1e-4 this is 24.5; the chapters round to "~30".
    assert math.isclose(c.t_per_rotation(1e-4), 1.15 * math.log2(1e4) + 9.2)
    assert 24 < c.t_per_rotation(1e-4) < 25
    assert 28 < c.t_per_rotation(1e-5) < 29


def test_ross_selinger_costs_more():
    assert c.t_per_rotation(1e-4, "rs") > c.t_per_rotation(1e-4, "rus")


def test_eps_splitting_incoherent_is_looser():
    n = 5e3
    inc = c.eps_per_rotation(1e-2, n, "incoherent")
    coh = c.eps_per_rotation(1e-2, n, "coherent")
    assert inc > coh
    assert math.isclose(inc, math.sqrt(1e-2 / n))
    assert math.isclose(coh, 1e-2 / n)


def test_t4_note_numbers():
    # hank-synthesis-note.txt, which used Watson et al. Eq. 80 (log2(2/eps)):
    # 5e3 rotations, total synthesis error 0.01-0.1 -> coherent 29-33 T, randomized 20-22 T.
    coh = sorted(c.t_per_rotation(c.eps_per_rotation(b, 5e3, "coherent"), "watson") for b in (0.1, 0.01))
    rnd = sorted(c.t_per_rotation(c.eps_per_rotation(b, 5e3, "incoherent"), "watson") for b in (0.1, 0.01))
    assert 28 <= coh[0] and coh[1] <= 34, coh
    assert 19 <= rnd[0] and rnd[1] <= 23, rnd
    # and the converged 12C row: ~1e11 rotations -> coherent 56-61, randomized 33-36
    coh = sorted(c.t_per_rotation(c.eps_per_rotation(b, 1e11, "coherent"), "watson") for b in (0.1, 0.01))
    rnd = sorted(c.t_per_rotation(c.eps_per_rotation(b, 1e11, "incoherent"), "watson") for b in (0.1, 0.01))
    assert 55 <= coh[0] and coh[1] <= 62, coh
    assert 32 <= rnd[0] and rnd[1] <= 37, rnd


# --- toffoli ----------------------------------------------------------------

def test_toffoli_conventions_explicit():
    assert c.toffoli_t(10, "textbook") == 70
    assert c.toffoli_t(10, "jones") == 40
    with pytest.raises(KeyError):
        c.toffoli_t(10, "whatever")


@pytest.mark.parametrize("k, anc", [(1, 0), (2, 1), (3, 1), (4, 3), (8, 7), (9, 7), (12, 10), (18, 16), (21, 18),
                                    (32, 31), (72, 70), (96, 94), (192, 190)])
def test_hwp_ancilla_is_k_minus_popcount(k, anc):
    """Author ruling E26 (r22, 2026-10-01), arxiv_1902_10673 App. A / arxiv_1709_06648: a group of k equal-angle
    rotations needs k - w(k) Toffolis and k - w(k) ancilla. The weight bits are the adders' carry ancilla and the
    lowest bit sits on a data qubit, so the weight register is not added on top (the pre-r22 count did)."""
    assert c.hwp_ancilla(k) == anc == k - bin(k).count("1")
    assert c.hwp_ancilla(k) == c.hwp_toffolis(k)
    assert c.hwp_group(k, 1e-4)["ancilla"] == anc


def test_hwp_weight_register_fits_in_the_ancilla():
    """The floor(log2 k)+1 weight bits live on the carry ancilla plus one data qubit, for every k."""
    for k in range(1, 1025):
        assert c.hwp_synth_rotations(k) <= c.hwp_ancilla(k) + 1
        assert c.hwp_ancilla(k) <= k - 1
    assert c.hwp_ancilla(32) + c.hwp_synth_rotations(32) == 37       # the retired double count, for the record
    with pytest.raises(ValueError):
        c.hwp_ancilla(0)


def test_rus_synthesis_holds_one_ancilla():
    """E26: repeat-until-success synthesis (bocharovRoettelerSvore2015) needs 1 ancilla, not 2."""
    assert c.RUS_ANCILLA == 1


def test_mcx():
    assert c.mcx_toffoli_count(2) == 1
    assert c.mcx_toffoli_count(4, "linear") == 5
    assert c.mcx_toffoli_count(4, "log") == 3


# --- formatting matches the table convention -------------------------------

def test_tex_range_table_style():
    assert c.tex_range(1.1e8, 4.8e8) == r"$1.1{\times}10^{8}$--$4.8{\times}10^{8}$"
    assert c.tex_range(5e4, 5e4) == r"$5{\times}10^{4}$"


def test_tex_range_factored():
    assert c.tex_range_factored(1.1e8, 4.8e8) == r"$1.1$--$4.8\times 10^{8}$"


# --- groups -----------------------------------------------------------------

def test_link_width_conflicts_are_recorded_not_hidden():
    ks = {k for k, _, _ in conflicts()}
    assert {"S36x3", "S72x3"} <= ks


def test_group_primitives_have_sources():
    for g in GROUPS.values():
        for name, p in g.primitives.items():
            if p.status is c.CircuitStatus.COMPILED:
                assert p.src, f"{g.name}.{name} is COMPILED but has no source"


# --- the contract, for every chapter model that exists ---------------------

@pytest.mark.parametrize("ch", estimates.available())
def test_model_contract(ch):
    m = estimates.load(ch)
    a = m.Assumptions()
    assert m.PUBLISHED, "no PUBLISHED numbers"
    for era, pub in m.PUBLISHED.items():
        assert era in c.ERAS
        r = m.model(a, era)
        assert r.era == era
        assert r.lq[0] <= r.lq[1] and r.hard_ops[0] <= r.hard_ops[1]
        # every T-count must be a sum of named primitives with a status
        if r.breakdown:
            for p in r.breakdown:
                assert isinstance(p.status, c.CircuitStatus)
                assert p.count >= 0 and p.t_each >= 0
            # the breakdown must land inside the quoted hard-op range (representative point)
            tot = r.breakdown_total()
            assert r.hard_ops[0] * 0.5 <= tot <= r.hard_ops[1] * 2.0, (tot, r.hard_ops)
        # every Assumptions field must be Tagged (provenance is mandatory)
    import dataclasses
    for f in dataclasses.fields(a):
        v = getattr(a, f.name)
        assert isinstance(v, c.Tagged), f"Ch{ch}.Assumptions.{f.name} has no provenance tag"


@pytest.mark.parametrize("ch", estimates.available())
def test_model_reproduces_published(ch, request):
    """Model vs the box as printed. A chapter may declare DISPUTED = {era: reason}
    for eras where its own stated inputs do not reproduce the box and the item is
    filed in NEEDS_AUTHOR.md; those become strict xfails, so the suite stays green
    while the dispute stays visible, and any silent resolution flips the test."""
    m = estimates.load(ch)
    a = m.Assumptions()
    disputed = getattr(m, "DISPUTED", {})
    for era, pub in m.PUBLISHED.items():
        r = m.model(a, era)
        if era in disputed:
            ok = c.close(r.lq, pub.lq, pub.rel_tol) and (pub.hard_ops is None or c.close(r.hard_ops, pub.hard_ops, pub.rel_tol))
            assert not ok, f"Ch{ch} {era} is marked DISPUTED but now reproduces: {disputed[era]}"
            continue
        assert c.close(r.lq, pub.lq, pub.rel_tol), (era, "lq", r.lq, pub.lq)
        if pub.hard_ops is not None:
            assert c.close(r.hard_ops, pub.hard_ops, pub.rel_tol), (era, "T", r.hard_ops, pub.hard_ops)


# --- ruling R-TOL: per-circuit rotation tolerance ----------------------------

def test_rtol_eps_syn_is_report_wide_1e2():
    assert c.EPS_SYN == 1e-2


def test_rtol_eps_rot_for_is_sqrt_budget_over_n():
    for n in (1, 338, 3792, 5e3, 3.78108e7, 5.0625e10):
        e = c.eps_rot_for(n)
        assert math.isclose(e, math.sqrt(1e-2 / n))
        assert math.isclose(n * e ** 2, c.EPS_SYN)            # randomized: errors add as N eps^2
        assert math.isclose(e, c.eps_per_rotation(1e-2, n, "incoherent"))
    assert math.isclose(c.eps_rot_for(1e4), 1e-3)
    assert math.isclose(c.eps_rot_for(1e6), 1e-4)
    assert math.isclose(c.eps_rot_for(100, eps_syn=1e-4), 1e-3)


def test_rtol_bigger_circuit_tighter_tolerance():
    assert c.eps_rot_for(1e8) < c.eps_rot_for(1e4) < c.eps_rot_for(1e2)
    # and the chapter formulas cost more per rotation as the circuit grows
    for model in ("rus", "rus-slope"):
        assert c.t_per_rotation(c.eps_rot_for(1e8), model) > c.t_per_rotation(c.eps_rot_for(1e4), model)


def test_rtol_rejects_empty_circuit():
    with pytest.raises(ValueError):
        c.eps_rot_for(0)
    with pytest.raises(ValueError):
        c.eps_rot_for(-5)


def test_rtol_matches_ch3_survey_check():
    # Ch. 3 already prices its survey at a 1e-2 per-shot synthesis budget (app11:90, the '25-32 T');
    # eps_rot_for must reproduce its incoherent RUS numbers at the same rotation counts.
    m = estimates.load(3)
    i = m.model(m.Assumptions(), "2033").intermediates
    for n, t in zip(i["rotations_survey"], i["c_T_survey_range"]):
        assert math.isclose(c.t_per_rotation(c.eps_rot_for(n)), t)


# --------------------------------------------------------------------------- #
# E20 (r18, 2026-10-01): every gauge-primitive rotation at the full fit 1.15 log2(1/eps) + 9.2
# --------------------------------------------------------------------------- #

from estimates import groups as G   # noqa: E402

# rotations per primitive = the paper's printed log coefficient / 1.15 (PAPER_COSTS.md sec. 6-7)
E20_N_ROT = {
    ("Z3", "U_Tr"): 1, ("Z3", "U_phi"): 1, ("Z3", "U_F"): 14, ("Z3", "U_FFT"): 2,
    ("2T", "U_Tr"): 11, ("2T", "U_F"): 1000, ("2T", "U_FFT"): 42,
    ("2O", "U_Tr"): 4, ("2O", "U_F"): 11370.1 / 1.15, ("2O", "U_FFT"): 42,
    ("S36x3", "U_Tr"): 7, ("S36x3", "U_phi"): 256, ("S36x3", "U_F"): 185898 / 1.15, ("S36x3", "U_FFT"): 102,
    ("S72x3", "U_Tr"): 7, ("S72x3", "U_phi"): 256, ("S72x3", "U_F"): 161996, ("S72x3", "U_FFT"): 448,
}


def test_e20_rotation_counts_per_primitive():
    for (g, p), n in E20_N_ROT.items():
        assert math.isclose(G.GROUPS[g].primitives[p].n_rot, n, rel_tol=1e-12), (g, p)
    for g in ("Z3", "2T", "2O", "S36x3", "S72x3"):
        for p in ("U_inv", "U_mul"):
            assert G.GROUPS[g].primitives[p].n_rot == 0          # Toffoli-only primitives do not move


@pytest.mark.parametrize("eps", [1e-3, 8.683e-6, 1e-8])
def test_e20_every_primitive_priced_at_the_full_fit(eps):
    L = math.log2(1 / eps)
    assert G.rotation_t(eps) == c.t_per_rotation(eps, "rus") == pytest.approx(1.15 * L + 9.2, rel=1e-15)
    for g, grp in G.GROUPS.items():
        for p, pc in grp.primitives.items():
            assert pc.t(eps) == pytest.approx(pc.t_const + pc.n_rot * (1.15 * L + 9.2), rel=1e-12), (g, p)
            assert pc.t(eps) == pytest.approx(pc.t_papers(eps) + 9.2 * pc.n_rot, rel=1e-12), (g, p)
            assert pc.t_report(eps) == pytest.approx(pc.t(eps), rel=1e-12), (g, p)   # 7 T per Toffoli: one price


def test_e20_examples_from_the_ruling():
    s36 = G.GROUPS["S36x3"].primitives
    assert s36["U_Tr"].n_rot == 7 and s36["U_phi"].n_rot == 256 and s36["U_FFT"].n_rot == 102
    assert s36["U_FFT"].t(1e-4) == pytest.approx(532 + 102 * (1.15 * math.log2(1e4) + 9.2), rel=1e-12)


@pytest.mark.parametrize("group,ham", [("2T", "I"), ("2O", "KS"), ("2O", "I"), ("S36x3", "KS"),
                                       ("S36x3", "I"), ("S36x3", "I_fft"), ("S72x3", "KS"), ("S72x3", "I")])
def test_e20_closed_form_c_t(group, ham):
    for d in (2, 3):
        n = G.c_t_rotations(group, ham, d)
        assert n == pytest.approx(G.c_t_log_coefficient(group, ham, d) / 1.15, rel=1e-12)
        for eps in (1e-4, 1e-8):
            assert G.c_t_full(group, ham, d, eps) == pytest.approx(G.c_t_stated(group, ham, d, eps) + 9.2 * n, rel=1e-12)


def test_e20_per_link_terms_and_papers_switch():
    eps = 1e-5
    for g, ham, d in (("2T", "KS", 2), ("S36x3", "KS", 2), ("S72x3", "KS", 3), ("Z3", "KS", 2)):
        grp = G.GROUPS[g]
        rot_mag = sum(G.PRIMCOST[ham][p](d) * grp.primitives[p].n_rot for p in G.MAGNETIC)
        assert G.magnetic_per_link(g, ham, d, eps) == pytest.approx(
            G.magnetic_per_link(g, ham, d, eps, papers=True) + 9.2 * rot_mag, rel=1e-12)
        key = "U_FFT" if grp.has_fft else "U_F"
        rot_el = G.PRIMCOST[ham]["U_F"](d) * grp.primitives[key].n_rot + (grp.primitives["U_phi"].n_rot if "U_phi" in grp.primitives else 0)
        assert G.electric_per_link(g, ham, d, eps) == pytest.approx(
            G.electric_per_link(g, ham, d, eps, papers=True) + 9.2 * rot_el, rel=1e-12)


# ---- legacy reproduction: one per paper, the papers' own slope-only convention (rule R2, superseded) ----

def _compose_papers(group, ham, d, eps, with_phi):
    grp, m = G.GROUPS[group], G.PRIMCOST[ham]
    t = sum(m[p](d) * grp.primitives[p].t_papers(eps) for p in G.MAGNETIC) + m["U_F"](d) * grp.primitives["U_F"].t_papers(eps)
    return t + (grp.primitives["U_phi"].t_papers(eps) if with_phi else 0.0)


@pytest.mark.parametrize("group,ham,with_phi,paper", [
    ("2T", "I", False, "arxiv_2208_12309 :564-566"),
    ("2O", "KS", False, "arxiv_2312_10285 :813"),
    ("S36x3", "KS", False, "arxiv_2405_05973 :666-681"),
    ("S72x3", "KS", True, "arxiv_2511_17437 :949-964"),
])
def test_legacy_papers_closed_form_reproduced(group, ham, with_phi, paper):
    """The paper's printed C_T is the tab:primcost composition of its tab:tgatecost at 1.15 log2(1/eps), no offset."""
    for d in (2, 3):
        for eps in (1e-4, 1e-8):
            assert _compose_papers(group, ham, d, eps, with_phi) == pytest.approx(
                G.c_t_stated(group, ham, d, eps), rel=1e-5), (paper, d, eps)


def test_legacy_fft_paper_reproduced():
    """arxiv_2408_00075: the transpiler H_3 (14 R_z) at 1e-4 is 213.93 T and the FFT rows are a + b log2(1/eps)
    in the papers' convention (tab:tgatecostbt; Sigma(36x3) 532 + 117.3 log, :835)."""
    assert G.GROUPS["Z3"].primitives["U_F"].t_papers(1e-4) == pytest.approx(213.93, abs=0.005)
    L = math.log2(1e4)
    for g, a, b in (("2T", 98, 48.3), ("2O", 216, 48.3), ("S36x3", 532, 117.3)):
        assert G.GROUPS[g].primitives["U_FFT"].t_papers(1e-4) == pytest.approx(a + b * L, rel=1e-12)
        assert G.papers_convention_t(G.GROUPS[g].primitives["U_FFT"], 1e-4) == pytest.approx(a + b * L, rel=1e-12)


# --------------------------------------------------------------------------- #
# R9 (H. Lamm, 2026-10-02): T-depth, factory parallelism and unified wall-time names
# --------------------------------------------------------------------------- #

def test_r9_factory_constants():
    assert c.REACTION_TIME_S == 10e-6
    assert c.FACTORY_BASELINE == 10
    assert c.FACTORY_SRC == "arxiv_1905_09749"
    assert c.DEPTH_EXPORTS == ("t_per_shot", "t_depth_per_shot", "f_star", "floor_wall_s",
                               "factories_for_1yr", "wall_first_result_s", "wall_campaign_s")


def test_r9_depth_exports_arithmetic():
    # factory.json, Ch. 4 2033 Tr M^-1 campaign: 1e6 shots, 8.5e8 T, D_T 2.1e8-8.5e8 -> F* 1-4, F_1yr ~270,
    # floor 67-270 yr, printed (serial) wall 27 yr.
    x = c.depth_exports(8.5e8, (2.1e8, 8.5e8), 1e6, 1e-6, 1e-4)
    assert x["f_star"] == pytest.approx((1.0, 8.5e8 / 2.1e8))
    assert x["floor_wall_s"] == pytest.approx((1e6 * 2.1e8 * 1e-5, 1e6 * 8.5e8 * 1e-5))
    assert x["floor_wall_s"][1] / c.SECONDS_PER_YEAR == pytest.approx(269.4, rel=1e-3)
    assert x["wall_serial_s"][0] == pytest.approx(1e6 * (8.5e8 * 1e-6 + 1e-4))
    assert x["serial_yr"][0] == pytest.approx(26.94, rel=1e-3)
    assert x["factories_for_1yr"][0] == pytest.approx(10 * 8.5e8 / c.SECONDS_PER_YEAR * 1e6 * 1e-6)
    assert x["baseline_ok"] == (False, False)
    assert x["fits_1yr"] == (False, False)


def test_r9_depth_exports_ranges_and_pairing():
    # Ch. 2 2033: T 6.7e8-9.6e8, D_T 4.7e7-1.31e8 -> F* cross-paired 5.1-20.4; shots paired lo with lo.
    x = c.depth_exports((6.7e8, 9.6e8), (4.7e7, 1.31e8), (5.7e4, 3.1e5), 1e-6, 1e-4)
    assert x["f_star"] == pytest.approx((6.7e8 / 1.31e8, 9.6e8 / 4.7e7))
    assert x["floor_wall_s"] == pytest.approx((5.7e4 * 4.7e7 * 1e-5, 3.1e5 * 1.31e8 * 1e-5))
    assert x["wall_serial_s"] == pytest.approx((5.7e4 * (6.7e8 * 1e-6 + 1e-4), 3.1e5 * (9.6e8 * 1e-6 + 1e-4)))
    # small benchmark: F* comfortably above the baseline, one year trivially reachable
    y = c.depth_exports(1.0e5, (1.26e3, 1.0e4), 2e3, 1e-6, 1e-4)
    assert y["baseline_ok"] == (True, True) and y["fits_1yr"] == (True, True)
    assert y["f_star"] == pytest.approx((10.0, 1e5 / 1.26e3))


def test_r9_depth_exports_takes_tagged_and_keeps_chapter_values():
    tg, t0 = c.Assumed(2e-6, "a slower machine"), c.Assumed(1e-3, "1 ms per shot")
    x = c.depth_exports(1e6, 1e5, 10, tg, t0)
    assert x["wall_serial_s"][0] == pytest.approx(10 * (1e6 * 2e-6 + 1e-3))       # no central override
    assert x["factories_for_1yr"][0] == pytest.approx(10 * 10 * 1e6 * 2e-6 / c.SECONDS_PER_YEAR)


def test_r9_depth_exports_rejects_nonsense():
    with pytest.raises(ValueError):
        c.depth_exports(1e3, 1e4, 1, 1e-6, 1e-4)          # depth > count
    with pytest.raises(ValueError):
        c.depth_exports(1e4, (2e3, 1e3), 1, 1e-6, 1e-4)   # reversed range
    with pytest.raises(ValueError):
        c.depth_exports(1e4, 1e3, 0, 1e-6, 1e-4)
    with pytest.raises(ValueError):
        c.depth_exports(1e4, 1e3, 1, 0, 1e-4)


_R9_LEGACY_NAMES = ("logical_cycle_s", "gate_time_s", "gate_time_us")


# R9 exports: xfail marker removed in E29 (2026-10-02); all nine chapters pass.
@pytest.mark.parametrize("ch", estimates.available())
def test_model_exports_depth(ch):
    """CONTRACT.md rules 11-13: unified Assumptions names, docstring header, and the R9 exports in every
    2028 / 2033 Result.intermediates, with f_star = t_per_shot / t_depth_per_shot (cross-paired)."""
    import dataclasses
    import re
    m = estimates.load(ch)
    a = m.Assumptions()
    names = {f.name for f in dataclasses.fields(a)}
    assert {"t_gate_s", "shot_overhead_s"} <= names, f"Ch{ch}: missing t_gate_s / shot_overhead_s"
    assert not names & set(_R9_LEGACY_NAMES), f"Ch{ch}: legacy names {names & set(_R9_LEGACY_NAMES)}"
    doc = m.__doc__ or ""
    assert re.match(rf"Ch\. {ch} — \S", doc), f"Ch{ch}: docstring must start 'Ch. {ch} — <title>.'"
    assert estimates.TEX[ch] in doc, f"Ch{ch}: docstring must name {estimates.TEX[ch]}"
    heads = [doc.find(h) for h in ("INPUTS AND SOURCES", "WHAT IS NOT DERIVED HERE", "WALL TIME AND FACTORIES")]
    assert all(i >= 0 for i in heads) and heads == sorted(heads), f"Ch{ch}: docstring sections missing or out of order"
    for era in m.PUBLISHED:
        if era not in ("2028", "2033"):
            continue
        i = m.model(a, era).intermediates
        missing = [k for k in c.DEPTH_EXPORTS if k not in i]
        assert not missing, f"Ch{ch} {era}: missing exports {missing}"
        if i.get("depth_status") == "not_established":
            # Accepted Ch. 2 studies establish counts, not a device schedule.
            assert ch == 2
            assert m.model(a, era).wall_time_s is None
            assert i["t_per_shot"][0] > 0
            assert all(i[k] is None for k in c.DEPTH_EXPORTS if k != "t_per_shot")
            continue
        for k in c.DEPTH_EXPORTS:
            if k == "wall_first_result_s" and era == "2028" and i[k] is None:
                continue
            lo, hi = c._band(i[k])
            assert lo > 0, (ch, era, k, i[k])
        t, d = c._band(i["t_per_shot"]), c._band(i["t_depth_per_shot"])
        assert c._band(i["f_star"]) == pytest.approx((t[0] / d[1], t[1] / d[0]), rel=1e-6), (ch, era)
