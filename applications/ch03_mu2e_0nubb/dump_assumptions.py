"""Print the Ch. 3 assumptions table from the live model.

Run from the repository root:

    python applications/ch03_mu2e_0nubb/dump_assumptions.py          # Markdown tables, grouped
    python applications/ch03_mu2e_0nubb/dump_assumptions.py --raw    # notes exactly as in the code

Every field of estimates.ch03_mu2e_0nubb.Assumptions is printed once, with its value, provenance tag
and source exactly as the code holds them. The note column is the code's note with internal
decision labels and edit history removed (see `clean`); --raw prints the notes verbatim. The script
fails if a field is missing from the grouping or listed twice, so the table cannot drift from the code.
"""

from __future__ import annotations

import re
import sys
from dataclasses import fields
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from estimates import ch03_mu2e_0nubb as m          # noqa: E402
from estimates.common import Provenance              # noqa: E402

# Groups follow the order of the Assumptions dataclass. Every field must appear exactly once.
GROUPS = [
    ("DOE reference budgets", [
        "rfi_lq_2028", "rfi_t_2028", "rfi_lq_2033", "rfi_t_2033"]),
    ("State preparation by projection, Eq. (Ngate11): 2033 and co-design cells", [
        "gap_mev", "gap_large_mev", "gap_small_mev", "dt_gev_inv", "n_trotter_quoted", "c_proj",
        "ch2_response_horizon_gev_inv", "eps_synth_total", "prep_lever"]),
    ("Registers, Eq. (Nq11)", [
        "qubits_per_orb", "n_orb_he4", "n_orb_sd", "n_orb_pf", "n_orb_jj44", "n_orb_jj55", "emax_first",
        "emax_converged_al", "emax_converged_ge", "mass_number_ge", "ho_orbitals_quoted", "n_anc",
        "n_anc_2028", "qrom_and_deficit", "test_qubits"]),
    ("Precision, shots, fault budget and wall-time conventions", [
        "eps_2028", "eps_mu2e", "eps_0nubb", "n_op_mu2e", "n_op_0nubb", "n_pairs_cross_check",
        "n_basis_scans", "n_survey_elements", "t_gate_s", "shot_overhead_s", "campaign_horizon_yr",
        "fault_budget"]),
    ("2028 benchmark: operator insertion, state load and amplitude estimation", [
        "t_be_topdown", "be_queries_2028", "gamma_trunc", "gamma_256", "gamma_512", "select_and_deficit",
        "insertion_t_per_toffoli", "qrom_t_per_amplitude", "d_amplitudes", "d_plausible", "he4_d_fid99",
        "he4_d_fid999", "he4_d_dipole_1pct", "qrom_passes", "loads_per_shot", "load_angle_error",
        "mlae_estimator_constant"]),
    ("2028 benchmark: the 4He dipole operator", [
        "he4_charge_radius_fm", "proton_radius_fm", "hbarc_mev_fm", "m_mu_mev", "m_nucleon_mev", "z_he4",
        "mass_number_he4"]),
    ("2033 deliverable: 27Al dipole factor, readout and T-depth", [
        "al27_charge_radius_fm", "z_al27", "mass_number_al27", "z_core_al27", "trotter_concurrency",
        "sd_spin_expectation", "direct_readout_variance", "projection_success", "projection_success_band",
        "first_result_scans", "eps_first_generic", "gap_al27_usdb_mev", "gap_al27_usdb_same_j_mev"]),
    ("2033 0nubb validation stage: 48Ca/48Ti truncation ladder and resonant-drive readout", [
        "gap_0plus_f7p3_ca48_mev", "gap_0plus_f7p3_ti48_mev", "gap_0plus_f7p3p1_ca48_mev",
        "gap_0plus_f7p3p1_ti48_mev", "gap_0plus_pf_ca48_mev", "gap_0plus_pf_ti48_mev",
        "gap_0plus_ensdf_ca48_mev", "gap_0plus_ensdf_ti48_mev", "eps_simplified", "lam_over_m_f7p3",
        "p_ref_f7p3", "lam_over_m_f7p3p1", "lam_over_m_pf", "rabi_drive_mev", "qdrift_eps", "rabi_shots"]),
    ("External comparison points and the utility box (not used in any T-count)", [
        "fq_light_nucleus_t", "ext_qubitization_sd_toffoli", "ext_qubitization_pb208_toffoli",
        "mu2e_tpc_musd", "utility_fraction", "nme_spread"]),
]

# Notes whose text in the code is mostly edit history; these one-line readings replace them in the
# cleaned table (--raw shows the originals).
NOTE_OVERRIDES = {
    "load_angle_error": "angle-rounding error of the load, 2 pi L / 2^b (Low, Kliuchnikov and Schaeffer); "
                        "set equal to the per-shot synthesis budget",
    "insertion_t_per_toffoli": "7 T per Toffoli, the report-wide convention (the source's temporary-AND count would give 4)",
    "qrom_t_per_amplitude": "'D-1 Toffolis per QROM pass' from the source, at the report's 7 T per Toffoli; "
                            "used only in the lookup-only comparison",
    "d_plausible": "lower end of the priced D band (1e2); D = 1e3 is kept only as a comparison point "
                   "(composed_fraction_at_d_hi)",
    "qrom_passes": "lookup-only accounting (load, uncompute, optional verification pass), kept for comparison; "
                   "the priced shot makes one complete load (loads_per_shot)",
    "trotter_concurrency": "commuting Trotter rotations in flight at once (12-24 of the 50-150 ancillas); "
                           "G >= 10 keeps the 10-factory baseline fed",
    "sd_spin_expectation": "<S_p>(27Al) in sd, an estimate pending a shell-model value",
    "direct_readout_variance": "per-shot variance of the direct (basis-rotated Z) readout of the SD one-body "
                               "operator, ~1 (estimate); worst case x8.3 (Popoviciu bound)",
    "projection_success": "heralded ground-state projection success p for deformed 27Al (estimate; computable "
                          "exactly in sd); sensitivity band 0.1-0.9",
    "projection_success_band": "sensitivity band on p",
    "first_result_scans": "the first result uses one Hamiltonian (one basis-cutoff scan)",
    "eps_first_generic": "report-wide first-result statistical error, quoted for comparison only",
    "eps_simplified": "30% statistical error on |M| for the 0nubb first light",
    "gap_0plus_ensdf_ca48_mev": "48Ca first excited 0+ (4283 keV); comparison only, the priced gaps are KB3G",
    "gap_0plus_ensdf_ti48_mev": "48Ti first excited 0+ (2997 keV); comparison only, the priced gaps are KB3G",
    "gap_al27_usdb_same_j_mev": "COMPUTED: 5/2+_2 - 5/2+_1, USDB in sd (comparison only)",
    "lam_over_m_f7p3": "lambda/|M| in f7/2 p3/2 from a schematic interaction; not derived here",
    "lam_over_m_f7p3p1": "lambda/|M| in f7/2 p3/2 p1/2 (schematic interaction); enters the qDRIFT count",
    "lam_over_m_pf": "lambda/|M| in full pf: Pauli-LCU lambda 208.4 over |M| ~ 0.29-0.35 (schematic interaction)",
    "p_ref_f7p3": "closed-f7/2 reference overlap in f7/2 p3/2 (schematic interaction)",
    "rabi_drive_mev": "theta|M|, the Rabi frequency of the resonant drive; kept <= 0.5 MeV for a <~6% Stark "
                      "bias (schematic interaction)",
    "rabi_shots": "Rabi-readout shots including a 5-10-detuning resonance scan; 250 give ~5% on |M| at "
                  "t_pi/2 (schematic interaction)",
    "utility_fraction": "'about 60% of the $315.7M Mu2e total project cost ... The fraction is not derived.' "
                        "The dollar value is linear in it",
}

# Parentheticals that carry only decision labels, dates or the previous value of a field.
_HISTORY = re.compile(r"\bR\d|\br\d\d\b|R-TOL|ruling|referee|retired|H\. Lamm|\bwas\b|NEEDS_AUTHOR|"
                      r"TRACKED_CHANGES|shot audit")
_FORBIDDEN = re.compile(r"\bR\d+\b|ruling|referee|retired|\bwas\b|NEEDS_AUTHOR|TRACKED_CHANGES|"
                        r"\bround\b|\bagent\b|verifier|workflow|earlier plan")


def _strip_parens(text: str) -> str:
    """Drop every balanced (...) group whose content matches _HISTORY."""
    out, i = [], 0
    while i < len(text):
        if text[i] == "(":
            depth, j = 0, i
            while j < len(text):
                depth += text[j] == "("
                depth -= text[j] == ")"
                if depth == 0:
                    break
                j += 1
            group = text[i:j + 1]
            if _HISTORY.search(group):
                while out and out[-1] == " ":
                    out.pop()
            else:
                out.append(group)
            i = j + 1
        else:
            out.append(text[i])
            i += 1
    return "".join(out)


def clean(name: str, note: str) -> str:
    if name in NOTE_OVERRIDES:
        return NOTE_OVERRIDES[name]
    s = _strip_parens(note)
    s = re.sub(r"\(app11:[\d,]+\)", "", s)                  # chapter line refs repeat the source column
    s = s.strip(" ;")
    return re.sub(r"\s+;", ";", re.sub(r"\s{2,}", " ", s))


def fmt_value(v) -> str:
    if isinstance(v, tuple):
        return "(" + ", ".join(repr(x) for x in v) + ")"
    return repr(v)


def fmt_prov(t) -> str:
    return {Provenance.CITED: "Cited", Provenance.ASSUMED: "Assumed",
            Provenance.STATED: "Stated-not-derived", Provenance.UNCITED: "Uncited"}[t.prov]


def main(raw: bool = False) -> int:
    a = m.Assumptions()
    names = [f.name for f in fields(a)]
    listed = [n for _, g in GROUPS for n in g]
    missing = [n for n in names if n not in listed]
    extra = [n for n in listed if n not in names]
    dup = sorted({n for n in listed if listed.count(n) > 1})
    if missing or extra or dup:
        print(f"grouping out of date: missing {missing}, unknown {extra}, duplicated {dup}", file=sys.stderr)
        return 1
    counts = {}
    leaks = []
    for title, group in GROUPS:
        print(f"#### {title}\n")
        print("| field | value | provenance | source | note |")
        print("|---|---|---|---|---|")
        for n in group:
            t = getattr(a, n)
            counts[fmt_prov(t)] = counts.get(fmt_prov(t), 0) + 1
            note = t.note if raw else clean(n, t.note)
            if not raw and _FORBIDDEN.search(note):
                leaks.append(n)
            cells = (f"`{n}`", fmt_value(t.value), fmt_prov(t), t.src, note)
            print("| " + " | ".join(c.replace("|", "\\|") for c in cells) + " |")
        print()
    print(f"{len(names)} fields: " + ", ".join(f"{k} {v}" for k, v in counts.items()))
    if leaks:
        print(f"note text still carries edit history: {leaks}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main(raw="--raw" in sys.argv[1:]))
