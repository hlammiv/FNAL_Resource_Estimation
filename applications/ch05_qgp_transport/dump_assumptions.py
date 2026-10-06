"""Print the Ch. 5 assumptions table (Markdown) from the live model.

Run from the repository root:

    python applications/ch05_qgp_transport/dump_assumptions.py          # Markdown tables, grouped
    python applications/ch05_qgp_transport/dump_assumptions.py --raw    # source and note verbatim

Every row is read from estimates.ch05_qgp_transport.Assumptions() at run time: the field name, its
value and its provenance tag exactly as the code holds them. The source and note columns are the
code's strings with internal decision labels and edit history removed (see `clean`); --raw prints
them verbatim. Every edit is an exact substring replacement that must find its target, and the run
fails if a label survives, so the cleaned table cannot drift silently from the code. Rows are
grouped by the section comments of the dataclass. The "used by" column lists the eras whose model
(_model_2028, _model_2033, _model_codesign, followed through every module function they call) reads
the field as `a.<name>`; a field read by none of them is a record kept for comparison, or a value the
chapter text quotes that no exported number depends on.
"""

from __future__ import annotations

import ast
import dataclasses
import inspect
import re
import sys
import textwrap
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from estimates import ch05_qgp_transport as m  # noqa: E402
from estimates.common import Tagged  # noqa: E402

# Display names for the section comments of the dataclass, matched on their leading words.
GROUPS = (
    ("shared physical constants", "Shared: physical constants"),
    ("per-link-step pricing", "Shared: gate-pricing conventions (both eras)"),
    ("2028 box", "2028 benchmark: lattice, group, Trotter step, shots and machine"),
    ("2028 tiers and workspace", "2028 benchmark: first result and primitive workspace"),
    ("2033 box", "2033 target: lattice, group, Hamiltonian and fermions"),
    ("2033 temperature", "2033 target: temperature and the in-situ T_c diagnostic"),
    ("thermal state by quench", "2033 target: quench preparation, budgets and campaign grid"),
    ("referee G7", "2033 target: Trotter step rule"),
    ("r25", "Shots, fit design, decay band, volume extensions and utility"),
)
ERA_ROOTS = (("2028", "_model_2028"), ("2033", "_model_2033"), ("codesign", "_model_codesign"))


# Notes replaced whole. Three of them, and the be_trim edit below, correct code notes that lag the
# model: the late-shot step counts are the model's fit_steps_grid -> fit_steps with
# trotter_sd_err_deepest_grid -> trotter_sd_err_deepest; the second spacing is priced at its own step
# (second_spacing_own_step); F* at two workspace qubits is f_star_2028_at_2_workspace; and the 2033
# group is Sigma(216x3).
NOTE_OVERRIDES = {
    "trotter_rule_2033": "each fit shot runs max(t / (a_t/10), N_sd) steps, N_sd = ceil(sqrt(Lambda_sd t^3 / eps)), "
                         "the fewest steps at which the state-dependent second-order estimate (connected variance, "
                         "Ch. 9) meets eps = 0.1. Only the late shots bind: 18 -> 19 and 26 -> 32 steps (estimate "
                         "0.149 -> 0.098 on the slow-end shot). The estimate is an upper estimate of the part that "
                         "grows with t. 'grid' restores the a_t/10 step everywhere",
    "second_spacing_a_s_fm": "an example finer spacing for the campaign's second spacing ('e.g. a_s = 0.15 fm') on the "
                             "same 3^3 lattice, priced at its own step (second_spacing_own_step); the choice is open",
    "workspace_2028": "primitive workspace on top of the 160 link qubits and 1-2 Hadamard-test ancillas, so 8 U_FFT "
                      "(2 ancillas each) run at once and F* >= 10 (f_star_2028_per_step); with only the 2 qubits any "
                      "primitive needs, F* = 1.66 (f_star_2028_at_2_workspace)",
    "hamiltonian_2033": "the 2033 box is priced at the improved Hamiltonian H_I (tab:primcost H_I row: 4 U_F, 2 U_Ph, "
                        "3(d-1)/2 U_Tr, 2 + 11(d-1) U_inv, 4 + 26(d-1) U_mul per link-step; Gustafson_S648_inprep). "
                        "'KS' restores the H_KS multiplicities for Sigma(216x3) (sensitivity)",
    "eps_syn": "total synthesis error per shot, report-wide; each circuit sets its own per-rotation tolerance "
               "eps_rot = sqrt(eps_syn / N_rot) (common.eps_rot_for, randomized synthesis), N_rot = synthesized "
               "rotations in one shot of that circuit. The papers' own fiducial is eps_T = 1e-8 total with a linear "
               "split (PAPER_COSTS.md sec. 3-4); app03:78 prices that as a sensitivity",
    "synthesis_convention": "the papers' tables with every rotation at the full fit 1.15 log2(1/eps) + 9.2 "
                            "(n_rot = log coefficient / 1.15), 7 T/Toffoli. The papers' slope-only price is kept "
                            "as the *_papers record",
    "dt_over_at_2028": "the state-dependent second-order estimate (trotter_check_2028) at eps = 0.1 sets the step; "
                       "six steps of a_t/4 reach t_max = 0.15 fm/c at 0.083. Same T per shot: the gate count of a "
                       "step does not depend on Delta t",
    "t_gate_s": "1 us per T gate, i.e. ~10 factories at the 10 us logical cycle",
    "shot_overhead_s": "per-shot overhead ~0.1 ms (register init, final readout, decode), report rule, as Ch. 8 "
                       "shot_overhead_s. Every wall time is serial on one machine: shots x (T x 1 us + t0)",
    "per_step_cut_levers": "the named per-step levers give 3-10x; the 2028 full benchmark needs >= 5.9x, so only "
                           "the upper part of the range closes it",
    "first_result_steps_2028": "the 2028 first result is the fall-off from t = 0 to one step (two timeslices), each "
                               "shot evolved to its own time; one step at its own eps is 9.47e4 T, inside 1e5",
    "first_result_falloff_2028": "NOT A BOX NUMBER. Relative fall-off delta of the normalized correlator between "
                                 "t = 0 and one step (0.045 fm/c; for an even correlator delta ~ t^2, so 0.1 at a_t "
                                 "would be ~0.02 here), used only to price the first result at 30% in the prose: "
                                 "N = 2 (1 - c0^2) / (c0^2 r^2 delta^2), c0 = cbar. delta is open",
    "group_2033": "Sigma(216x3) subset SU(3), |G| = 648, 11 qubits/link in the draft's qubit encoding "
                  "(Gustafson_S648_inprep, in preparation)",
    "T_over_Tc_2033": "T = 1.5 T_c^lat, the same convention as the 2028 benchmark",
    "ceiling_2033_t": "'the RFI's ~1e9-hard-op 2033 envelope' (app03:149), the 2033 hard-operation ceiling",
    "shots_per_gridpt": "RECORD: 10^3/grid-pt (eps^-2 = 25 x ~40) of the 4 fm/c reference campaign; the priced "
                        "shots come from fisher_shots",
    "grid_k": "RECORD: x 3 k of the 4 fm/c reference campaign; the priced campaign uses campaign_k",
    "grid_t": "RECORD: x 2 t-windows of the reference campaign; the priced campaign uses the two fit points of fit_x",
    "mu_b_zero_fraction": "mu_B in {0,200,400,600} MeV: 1 of 4 at mu_B=0; the campaign is '4 points at mu_B=0 and "
                          "12 at mu_B>0'. With the quench preparation every point costs the same per shot (one "
                          "thermal route)",
    "cbar": "normalized correlator Cbar ~ 0.16 of the +/-1 Hadamard-test estimator (variance factor ~40 = "
            "1/Cbar^2); not derived",
    "target_first": "the first result, Gamma(k_min) at T = 1.5 T_c^lat, mu_B = 0, to 30% statistical",
    "tau_fast_1k": "fast end of the decay band = factor x 1/k_min, the collisionless (free-streaming) dephasing "
                   "time of the transverse momentum density at k_min ~ 7-9 T (0.0955 fm/c on 3^3 at a_s = 0.2 fm). "
                   "The hydrodynamic T/((eta/s) k^2) = 0.049 fm/c, which the chapter says has no basis at this k, "
                   "is kept for comparison (tau_decay_at_eta_s_assumed_fm)",
    "second_spacing_own_step": "the campaign's second spacing is priced at its own step under the step rule (a_t "
                               "scaled with a_s at fixed anisotropy, so Delta t ||H|| is held). Pricing it at the "
                               "first spacing's per-shot cost would let Delta t ||H|| grow by 4/3 and break the step "
                               "rule; False gives that comparison",
}

# Exact substring edits on the remaining notes: (old, new). Each old string must be present.
NOTE_EDITS = {
    "be_trim": [("Sigma(72x3) BE not yet synthesized", "Sigma(216x3) BE is not yet synthesized")],
    "shots_per_timeslice": [(" (NEEDS_AUTHOR)", "")],
    "hop_group_2033": [("s648: ", ""), ("NEEDS_AUTHOR", "Open.")],
    "nonfreezing_HI_assumed": [("s648 ruling (2): ", ""), ("; author ruling 2026-10-05)", ")")],
    "hop_rule": [(" (r17 ruling (e), H. Lamm 2026-10-01)", ""), ("(E21 (1): ", "(")],
    "fermion_synthesis": [(" (r17)", "")],
    "quench_c": [(", NEEDS_AUTHOR)", ")")],
    "trotter_eps": [("G7: ", "")],
    "trotter_g2": [("G7: ", "")],
    "trotter_n_adj": [("G7: ", "")],
    "trotter_dispersion_2033": [("s648: ", "")],
    "sigma72_emax_over_fund": [(" (G7). NORMALIZATION (2026-10-04, Claude's decision, NEEDS_AUTHOR closed; stated at "
                                "app03 eq:Hks)", ". Normalization (stated at app03 eq:Hks)"),
                               ("the workflow's string-tension tuning", "the string-tension tuning")],
    "fit_x": [(" (NEEDS_AUTHOR: which transient the fit excludes)", "; which transient the fit excludes is open")],
    "campaign_k": [("R4 / shot audit: ", "")],
    "tau_slow_2piT": [(" (0.196 at the retired 160 MeV)", "")],
}

# Source strings: decision labels and dates removed, file and paper references kept.
SRC_EDITS = {
    "eps_syn": [("TRACKED_CHANGES.md R-TOL (H. Lamm, 2026-09-29)", "report-wide convention (common.eps_rot_for)")],
    "first_result_steps_2028": [(": author rulings R1-R10 (H. Lamm 2026-10-02)", " (author's choice)")],
    "workspace_2028": [(": author rulings R1-R10 (H. Lamm 2026-10-02)", " (author's choice)")],
    "target_first": [(": author rulings R1-R10 (H. Lamm 2026-10-02)", " (author's choice)")],
    "faults_per_shot": [("TRACKED_CHANGES.md 'RULINGS, NEEDS_AUTHOR round A' R3 (2026-09-28)",
                         "report-wide convention")],
    "shots_per_gridpt": [("app03 r23 box (retired r25)", "4 fm/c reference campaign")],
    "grid_k": [("app03 r23 box (retired r25)", "4 fm/c reference campaign")],
    "grid_t": [("app03 r23 box (retired r25)", "4 fm/c reference campaign")],
    "n_grid_prose": [("app03 r23 box (retired r25)", "4 fm/c reference campaign")],
}

_FORBIDDEN = re.compile(r"\bR\d+\b|\br\d\d\b|R-TOL|\b[Rr]uling|referee|retired|\b[Ww]as\b|NEEDS_AUTHOR|"
                        r"TRACKED_CHANGES|\bround\b|\bagent\b|verifier|workflow|earlier plan|Claude|supersede|"
                        r"H\. Lamm|\bs648\b|\bG7\b|\bE\d\d\b")


def _edit(name: str, text: str, edits: dict) -> str:
    for old, new in edits.get(name, ()):
        if old not in text:
            raise RuntimeError(f"{name}: expected text not found, update the edit: {old!r}")
        text = text.replace(old, new)
    return text


def clean(name: str, note: str) -> str:
    if NOTE_OVERRIDES.get(name):
        return NOTE_OVERRIDES[name]
    return _edit(name, note, NOTE_EDITS)


def clean_src(name: str, src: str) -> str:
    return _edit(name, src, SRC_EDITS)


def _groups() -> dict[str, str]:
    """field name -> display group, from the '# ---- title ----' comments in the class source."""
    out, current = {}, None
    for line in inspect.getsource(m.Assumptions).splitlines():
        s = line.strip()
        if s.startswith("# ---"):
            title = s.strip("#- ").strip()
            for key, name in GROUPS:
                if title.startswith(key):
                    current = name
                    break
            else:
                raise RuntimeError(f"unmapped section comment: {title!r}")
            continue
        hit = re.match(r"(\w+)\s*:\s*Tagged\b", s)
        if hit:
            out[hit.group(1)] = current
    return out


def _closure() -> dict[str, str]:
    """era -> concatenated source of every module-level function its model reaches (plus the
    Assumptions.T_2033 method where it is called), so a field read in a helper counts."""
    funcs = {n: f for n, f in vars(m).items()
             if inspect.isfunction(f) and f.__module__ == m.__name__}
    src = {n: textwrap.dedent(inspect.getsource(f)) for n, f in funcs.items()}
    calls = {}
    for n, s in src.items():
        names = {x.id for x in ast.walk(ast.parse(s)) if isinstance(x, ast.Name)}
        calls[n] = {c for c in names if c in funcs and c != n}
    out = {}
    for era, root in ERA_ROOTS:
        seen, todo = set(), [root]
        while todo:
            f = todo.pop()
            if f not in seen:
                seen.add(f)
                todo.extend(calls[f])
        text = "\n".join(src[f] for f in sorted(seen))
        if ".T_2033(" in text:
            text += "\n" + inspect.getsource(m.Assumptions.T_2033).replace("self.", "a.")
        out[era] = text
    return out


def _used_by(name: str, closure: dict[str, str]) -> str:
    pat = re.compile(rf"\ba\.{re.escape(name)}\b")
    eras = [era for era, text in closure.items() if pat.search(text)]
    return ", ".join(eras) if eras else "-"


def _value(v) -> str:
    if isinstance(v, tuple):
        return "(" + ", ".join(repr(x) for x in v) + ")"
    if isinstance(v, str):
        return v
    return repr(v)


def _cell(s: str) -> str:
    return s.replace("|", "\\|").replace("\n", " ")


def main(raw: bool = False) -> int:
    a = m.Assumptions()
    leaks: list[str] = []
    groups = _groups()
    closure = _closure()
    fields = [f.name for f in dataclasses.fields(a)]
    missing = [n for n in fields if n not in groups]
    if missing:
        raise RuntimeError(f"fields outside any section: {missing}")
    counts: dict[str, int] = {}
    for n in fields:
        v = getattr(a, n)
        assert isinstance(v, Tagged), n
        counts[v.prov.value] = counts.get(v.prov.value, 0) + 1
    print(f"{len(fields)} fields: " + ", ".join(f"{k} {c}" for k, c in counts.items()))
    order = []
    for n in fields:
        if groups[n] not in order:
            order.append(groups[n])
    for g in order:
        print()
        print(f"#### {g}")
        print()
        print("| field | value | provenance | used by | source | note |")
        print("|---|---|---|---|---|---|")
        for n in fields:
            if groups[n] != g:
                continue
            v = getattr(a, n)
            src = v.src if raw else clean_src(n, v.src)
            note = v.note if raw else clean(n, v.note)
            if not raw and (_FORBIDDEN.search(src) or _FORBIDDEN.search(note)):
                leaks.append(n)
            print(f"| `{n}` | {_cell(_value(v.value))} | {v.prov.value} | {_used_by(n, closure)} "
                  f"| {_cell(src) or '-'} | {_cell(note) or '-'} |")
    if leaks:
        print(f"source or note text still carries internal labels: {leaks}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main(raw="--raw" in sys.argv[1:]))
