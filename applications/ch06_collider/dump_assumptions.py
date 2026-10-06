"""Print the Ch. 6 assumption table from the live model, as Markdown.

    python applications/ch06_collider/dump_assumptions.py               # the README table
    python applications/ch06_collider/dump_assumptions.py --code-notes  # the model's own annotation per field

Name, value, provenance tag and source key are read from estimates.ch06_collider.Assumptions()
at run time; nothing in those columns is typed here. The "export" column is also computed: the
script perturbs each field (numbers by 7% or +1, a string to its documented alternative), reruns
the three eras and reports whether any export moves. The exports are what the chapter report of
`python -m estimates --chapter ch06` lists: logical qubits, per-shot T, shots, wall time and
epsilon_l of each era, the depth and factory keys, and every number in the chapter's
resources.json rows. "yes" means the field moves one of them, "no" means it feeds only
intermediates or cross-checks, "guarded" means the model accepts only the value shown.

The note column is a one-line reader's summary written for this table; `--code-notes` prints the
annotation the model carries instead. The script fails if a field is added to or removed from the
dataclass without a matching entry below, so the table cannot silently drift from the code.
"""

from __future__ import annotations

import json
import math
import sys
from dataclasses import fields, replace
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from estimates import ch06_collider as m  # noqa: E402
from estimates.common import Provenance, Tagged  # noqa: E402

TAG = {Provenance.CITED: "Cited", Provenance.ASSUMED: "Assumed",
       Provenance.STATED: "Stated-not-derived", Provenance.UNCITED: "Uncited"}

# Valid alternatives for the string-valued switches (the guards in Assumptions.__post_init__ list them).
ALT = {"link_width_source": "chapter", "synthesis_errors": "coherent", "hop_share": "draft",
       "trotter_rule": "grid", "hamiltonian_dipole_static": "KS", "group_2028": "S72x3",
       "group_2033": "S216x3", "group_dipole_static": "S72x3", "hamiltonian": "I",
       "synthesis_model": "rus-slope", "fermion_synthesis": "rus-slope", "toffoli_convention": "other",
       "eth_assumed": "none"}

SECTIONS = [
    ("Reference budgets and report conventions", [
        ("rfi_lq_2028", "DOE 2028 register envelope; comparison only"),
        ("rfi_t_2028", "DOE 2028 per-shot T reference; the 2028 row reports its ratio to it (over_reference_t)"),
        ("rfi_lq_2033", "DOE 2033 register marker; comparison only"),
        ("rfi_t_2033", "DOE 2033 per-shot T reference; comparison only"),
        ("rfi_eps_l", "DOE logical error rate; used to count expected faults per shot at that rate"),
        ("clean_shot_faults", "expected logical faults allowed per shot; required eps_l = 0.1 / (T per shot)"),
        ("seconds_per_year", "converts walls to years"),
        ("hbarc_GeV_fm", "converts lattice momenta, times and temperatures to GeV and fm"),
        ("t_gate_s", "1 us per T gate on one machine"),
        ("shot_overhead_s", "per-shot initialization, readout and decoding time"),
        ("eps_syn", "total rotation-synthesis error per shot; each circuit sets eps_rot = sqrt(eps_syn / N_rot)"),
        ("eps_rot", "placeholder tolerance; every circuit replaces it with its own sqrt(eps_syn / N_rot) before pricing"),
        ("eps_rot_before_round_e", "fixed tolerance read only by comparison intermediates (*_at_1e-4)"),
        ("synthesis_errors", "randomized synthesis: rotation errors add as N eps^2"),
        ("synthesis_model", "gauge-table rotations at the full fit 1.15 log2(1/eps) + 9.2 T"),
        ("fermion_synthesis", "hop and mass rotations at the same full fit"),
        ("toffoli_convention", "7 T per Toffoli"),
        ("link_width_source", "qubits per link from the compiled registers in groups.py (2 for Z3, 9 for Sigma(72x3))"),
        ("hamiltonian", "Kogut-Susskind H_KS; selects the per-link primitive multiplicities (groups.PRIMCOST)"),
        ("n_c", "3 colours, one Jordan-Wigner mode per colour per site"),
        ("n_stag", "3 staggered fields (u, d, s-like)"),
    ]),
    ("2028 benchmark: Z3, 2+1D, 4^2", [
        ("group_2028", "Z3 subgroup of SU(3), 2 qubits per link"),
        ("dim_2028", "spatial dimension; links per site = D"),
        ("L_2028", "4^2 lattice: V = 16 sites, 32 links"),
        ("n_trot_2028", "evolution steps per shot"),
        ("ramp_steps_2028", "adiabatic ramp steps per shot, each at the per-step cost"),
        ("shots_2028", "shots for the benchmark"),
        ("delta_rel_2028", "per-channel accuracy 1e3 shots buy at n_s = 0.1; a check, not an input to the shots"),
        ("workspace_links_2028", "hop links run at once, each with 7 phasing + 4 synthesis ancilla: 22 ancilla, 230 LQ"),
        ("n_anc_2028", "7 + 1 ancilla of a register that runs every rotation in series (216 LQ); comparison only"),
        ("n_pauli_hop_diagonal", "Pauli strings per colour-flavour copy in the Z3-dressed hop (derived and tested)"),
        ("n_pauli_mass", "one Z string per copy: one phasing group per site for the staggered mass"),
        ("hwp_adder_depth_2028", "Toffoli layers of the k = 9 Hamming-weight adder (log tree to serial); sets the 2028 T-depth"),
        ("dt_multi_over_a", "step sizes of the combined 2028 runs; sets which times are sampled, no T change"),
        ("t_step_gap_target", "per-step T target of the algorithm gap; gives the 44x cut"),
    ]),
    ("2033 species / fragmentation circuit: Sigma(72x3), 2+1D, 4x8", [
        ("group_2033", "Sigma(72x3) subgroup of SU(3), 9 qubits per link"),
        ("dim_2033", "the 2033 lattice is L_perp x L_par"),
        ("L_perp", "transverse extent, 4 sites"),
        ("L_par_uncut", "full axis length 16 (1828 LQ); comparison"),
        ("L_par_cut", "cut axis lengths 8-10 (964-1180 LQ); sensitivity"),
        ("L_par_priced", "axis length at which register and T are priced (V = 32)"),
        ("n_anc_2033", "ancilla allowance (linear-response register)"),
        ("a_fm", "lattice spacing"),
        ("t_had_over_a", "with n_trot_2033 fixes the step: Delta t = 100 a / 1000 = 0.1 a"),
        ("n_trot_2033", "with t_had_over_a fixes Delta t = 0.1 a"),
        ("n_trot_eq", "the 1e2-1e3 range of eq. Ngate_collider; a check"),
        ("window_over_a", "evolution window 30-50 a: 300-500 steps at 0.1 a"),
        ("ramp_steps_2033", "adiabatic ramp steps; each end of the range is its own circuit"),
        ("readout_ramps_2033", "the ramp is run backwards once per shot for the species map N_h = U Pi_h U^dag"),
        ("hop_share", "colour register map and parity flags computed once per link and held ('draft' recomputes them)"),
        ("n_pauli_hop_moved", "2 strings (XX+YY) of an undressed hop; read only by a comparison hop rule"),
        ("moves_per_field", "group multiplications per field; read by the current-insertion count"),
        ("insertion_links", "8-link Wilson line of the bilocal insertion: 2 x 7 + 2 = 16 U_mul"),
        ("insertions_per_shot_b", "bilocal-current insertions per fragmentation shot; reported, not added to the shot"),
        ("insertion_bilinear_rotations", "R_Z per half of the controlled bilinear in one insertion"),
        ("insertion_bilinear_toffolis", "Toffolis per half of the controlled bilinear"),
        ("source_momentum_GeV", "nominal parton momentum; kinematic intermediates only"),
        ("lambda_qcd_GeV", "puts 1/Lambda_QCD in fm/c beside the window; intermediates only"),
        ("qubitized_trim", "assumed per-shot cut from a block encoding; reported, not applied"),
        ("lund_kappa_GeV_fm", "string tension for the yo-yo windings on the 0.8 fm axis"),
        ("trotter_rule", "the state-dependent Trotter estimate sets the steps ('grid' keeps 0.1 a and 0.01 fm/c)"),
        ("trotter_eps", "Trotter error target per shot; sets the 575-step top-mode shot and the dipole steps"),
        ("trotter_g2", "coupling for the worst-case Trotter bound only"),
        ("trotter_temp_lat", "temperature of the free-field Trotter checks of the 2033 window (vacuum)"),
        ("depth_link_const_2033", "hop Toffoli layers per link per step, hop links one at a time"),
        ("depth_link_rot_2033", "hop rotation layers per link per step (each one rotation deep)"),
        ("depth_step_const_2033", "Toffoli layers per step outside the hop (electric, magnetic, mass)"),
        ("depth_step_rot_2033", "rotation layers per step outside the hop"),
    ]),
    ("2033 shots and campaign", [
        ("n_s", "mean multiplicity of a species class per event"),
        ("n_meson", "mesons per event in the ratio's denominator"),
        ("pair_clustering", "c in (c / n_s + 1 / n_M) / delta^2: 1 if n_s counts pairs, 2 if B + Bbar"),
        ("delta_first_2033", "relative error of the first result"),
        ("trend_size", "end-to-end change of B/M across the three source momenta"),
        ("trend_sigma", "significance required of that trend"),
        ("n_source_momenta", "source momenta in the trend"),
        ("delta_rel_2033", "5-10% accuracy range of the inclusive-count formula; comparison"),
        ("shots_samp_printed", "1/(n_s delta^2) at 5-10%; checked against the formula"),
        ("hadamard_amplitude", "bilocal-current amplitude A of the Hadamard test; shots scale as A^-2"),
        ("eps_stat", "statistical accuracy per z-bin of D_pi^q(z)"),
        ("n_z", "z-bins"),
        ("n_tens", "tensor components / spin projections of the bilocal current"),
        ("n_boost", "source momenta for the LaMET extrapolation"),
        ("shots_ff_printed", "fragmentation shots 3e3 x A^-2; written to resources.json (shots_fragmentation)"),
        ("ff_z_range", "z range split by the bins; sets the N_h shot factor (4.7-23x, not applied)"),
        ("n_spacings", "lattice spacings in the requirements row; bookkeeping only"),
        ("n_qhat_T", "q-hat temperatures in the requirements row; bookkeeping only"),
        ("n_species_channels", "species channels read from the same records; not used in a number"),
        ("campaign_horizon_yr", "campaign horizon; sets the cut factors"),
    ]),
    ("2033 dipole rows: pure gauge, quench-prepared medium", [
        ("group_dipole_static", "static row on Sigma(216x3), 11 qubits per link"),
        ("hamiltonian_dipole_static", "static row uses the improved Hamiltonian H_I (U_phi applied twice)"),
        ("static_dims", "static row lattice 3^3 (81 links)"),
        ("static_src_qubits", "two static colour sources, one qutrit each in two qubits"),
        ("lightlike_dims", "light-like row lattice 3x3x5 on 2O (135 links)"),
        ("lightlike_src_qubits", "two SU(2) fundamental sources, one qubit each"),
        ("dipole_n_anc", "ancilla allowance of each dipole row"),
        ("dipole_a_fm", "dipole lattice spacing; r = a"),
        ("dipole_dt_fm", "time grid; each shot runs max(t / 0.01 fm/c, N_sd) steps"),
        ("dipole_times_fm", "the two times (path lengths) of the log slope"),
        ("dipole_T_GeV", "medium temperature, aT ~ 0.45"),
        ("quench_ramp_a", "coupling ramp from the electric vacuum over one lattice unit: 20 steps"),
        ("quench_c_range", "thermalization time c / T before the sources enter: 90-282 steps"),
        ("eth_assumed", "eigenstate thermalization / typicality at these couplings; label only"),
        ("thermal_variance_factor", "estimator penalty on the binomial counts (none on a pure state)"),
        ("kappa_over_T3", "quenched heavy-quark kappa / T^3; sets the static decay rate"),
        ("static_colour_factor", "Gamma = kappa r^2 / 3 at small r"),
        ("qhat_GeV2_fm", "pure-glue q-hat near 1.5 T_c; sets the light-like decay rate"),
        ("dipole_P0", "singlet survival after the string build"),
        ("n_dipole_T", "temperatures in the scan; the scan wall is 3x the one-temperature wall"),
    ]),
    ("3+1D extension (not a box)", [
        ("L_s_3d", "4^3 lattice: V = 64, 192 links"),
        ("n_trot_3d", "steps, no ramp"),
        ("ch10_stretch_lq_printed", "Ch. 10's printed register on the same lattice; cross-reference"),
        ("ch10_stretch_n_anc", "Ch. 10's ancilla on that lattice; cross-reference"),
    ]),
    ("Systematics, clocks and utility (context numbers)", [
        ("syst_discretization", "O(a^2) at a = 0.1 fm"),
        ("syst_volume", "finite volume at L_perp = 4a"),
        ("syst_subgroup", "Sigma(72x3) truncation"),
        ("n_z_fcc", "Z bosons at FCC-ee Tera-Z"),
        ("n_z_lep", "Z bosons at LEP"),
        ("lep_error_divisor", "FCC-ee statistical gain over LEP"),
        ("utility_share", "share of the capital anchor assigned to the program"),
        ("cms_capital_musd", "U.S. CMS detector capital, $M"),
        ("n_instances_utility", "application instances sharing the value"),
    ]),
    ("Comparison values (read only by comparison intermediates)", [
        ("retired_hop_2028", "a hand-carried 2028 hop figure per step that the derived hop replaces"),
        ("retired_hop_insertion_2033", "a hand-carried 2033 hop + insertion figure per step"),
        ("retired_hop_insertion_3d", "a hand-carried 3+1D hop + insertion figure per step"),
        ("retired_n_trot_2028", "a 20-step 2028 variant"),
        ("round_d_n_trot_2028", "a 1 ramp + 3 evolution 2028 variant"),
        ("twosteps_n_trot_2028", "a 1 ramp + 1 evolution 2028 variant at eps = 1e-4"),
        ("draft_hop_S72x3_t_const", "printed per-link, per-field hop constant of the unpublished gate-count draft"),
        ("draft_hop_S72x3_t_log", "its log2(1/eps) coefficient; the model builds the hop from the draft's gate tables instead"),
    ]),
]

EXPORT_KEYS = ("t_per_shot", "t_depth_per_shot", "f_star", "floor_wall_s", "factories_for_1yr",
               "wall_first_result_s", "wall_campaign_s")


def fingerprint(a) -> str:
    out = {}
    for era in m.PUBLISHED:
        r = m.model(a, era)
        out[era] = {"lq": r.lq, "t": r.hard_ops, "shots": r.shots, "wall": r.wall_time_s, "eps_l": r.epsilon_l,
                    **{k: r.intermediates.get(k) for k in EXPORT_KEYS}}
        for n, (_, lq, t, extra) in enumerate(m.INSTANCE_ROWS(a, era, r)):
            out[f"{era}:{n}"] = {"lq": lq, "t": t, **{k: v for k, v in extra.items() if not isinstance(v, str)}}
    return json.dumps(out, default=str, sort_keys=True)


def _bump(x):
    if isinstance(x, bool) or not isinstance(x, (int, float)):
        return None
    return x + 1 if isinstance(x, int) else x * 1.07


def moves_export(a0, name: str, base: str) -> str:
    t = getattr(a0, name)
    v = t.value
    if isinstance(v, str):
        trials = [ALT[name]] if name in ALT else []
    elif isinstance(v, tuple):
        trials = [tuple(_bump(x) for x in v), v[:-1] + (_bump(v[-1]),)]
    else:
        trials = [_bump(v)]
    trials = [x for x in trials if x is not None and not (isinstance(x, tuple) and None in x)]
    if not trials:
        return "no"
    guarded = True
    for nv in trials:
        try:
            a = replace(a0, **{name: Tagged(nv, t.prov, t.src, t.note)})
            fp = fingerprint(a)
        except (ValueError, KeyError, TypeError):
            continue
        guarded = False
        if fp != base:
            return "yes"
    return "guarded" if guarded else "no"


def fmt(v) -> str:
    if isinstance(v, tuple):
        return "(" + ", ".join(fmt(x) for x in v) + ")"
    if isinstance(v, float):
        short = f"{v:.6g}"
        return short if float(short) == v else repr(v)
    return str(v)


def main(code_notes: bool = False) -> None:
    a = m.Assumptions()
    names = [f.name for f in fields(a)]
    listed = [n for _, rows in SECTIONS for n, _ in rows]
    missing, extra = set(names) - set(listed), set(listed) - set(names)
    if missing or extra or len(listed) != len(set(listed)):
        sys.exit(f"table out of step with the model: missing {sorted(missing)}, unknown {sorted(extra)}")
    base = fingerprint(a)
    counts = {}
    for n in names:
        p = TAG[getattr(a, n).prov]
        counts[p] = counts.get(p, 0) + 1
    print(f"{len(names)} fields: " + ", ".join(f"{k} {v}" for k, v in counts.items()) + ".")
    for title, rows in SECTIONS:
        print(f"\n**{title}**\n")
        print("| name | value | tag | source | export | note |")
        print("|---|---|---|---|---|---|")
        for n, note in rows:
            t = getattr(a, n)
            text = t.note.replace("|", "/") if code_notes else note
            print(f"| `{n}` | {fmt(t.value)} | {TAG[t.prov]} | {t.src or '-'} | {moves_export(a, n, base)} | {text} |")


if __name__ == "__main__":
    main(code_notes="--code-notes" in sys.argv[1:])
