"""Print the Ch. 10 assumptions table (Markdown) from the live model.

    python applications/ch10_finite_density/dump_assumptions.py

Every row is read from estimates.ch10_finite_density.Assumptions(): field name, value, provenance tag,
source and note, exactly as the code holds them. Only the grouping into instances is added here; the
script fails if a field is missing from the grouping or listed twice, so the table cannot drift from
the dataclass.
"""

from __future__ import annotations

import dataclasses
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from estimates.ch10_finite_density import Assumptions  # noqa: E402
from estimates.common import Provenance  # noqa: E402

GROUPS = (
    ("Shared: lattice theory and encoding",
     ("nc", "group", "ham")),
    ("Shared: synthesis, Toffoli and hop readings",
     ("eps_syn", "rot_synthesis_paper", "toffoli_convention", "hop_synthesis", "hop_undo", "hop_share",
      "hop_mcx", "hop_phasing")),
    ("Shared: wall time and logical error",
     ("t_gate_s", "shot_overhead_s", "faults_per_shot", "eps_l_floor")),
    ("2028 benchmark: 2+1D SU(3), V = 2^2, ground state at mu_B > 0",
     ("d_2028", "L_2028", "n_stag_2028", "n_anc_2028", "t_cap_2028", "n_trotter_2028", "n_trotter_2028_original",
      "shots_2028", "sigma2_2028", "n_sectors_2028", "t_depth_2028")),
    ("2033 target: 3D SU(3), V = 2^3, register and step count",
     ("d_2033", "L_2033", "n_stag_2033", "n_anc_2033", "t_ref_2033", "beta_h_2033", "c_mix", "t_depth_2033")),
    ("2033 target: QSVT/TPQ thermal state (the priced route)",
     ("delta_beta", "filter_success_prob", "aa_call_rule", "ref_rot_per_qubit", "qsvt_query_cost")),
    ("2033 target: grid, accuracy targets and shots",
     ("grid_T_mev", "grid_mu_mev", "n_T", "n_mu", "eps_stat", "eps_chi4_chi2_target", "eps_pressure_target",
      "eps_first_result", "nb_toy_k2_priced", "nb_toy_k2_band", "nb_toy_k2_mid", "nb_toy_quark_k2", "fr_rw_pairs",
      "grid_rw_pairs_k2lo", "grid_rw_pairs_k2mid", "campaign_horizon_yr")),
    ("Post-2033 stretch: Sigma(72x3), V = 4^3 (register only)",
     ("group_stretch", "d_stretch", "L_stretch", "n_stag_stretch", "n_anc_stretch")),
    ("Comparison routes (computed, not in the headline numbers)",
     ("c_be_2033_superseded", "gibbs_2x2_min_steps", "gibbs_2x2_t_per_shot_superseded", "gibbs_beta_h_2028",
      "gibbs_sweeps", "gibbs_eps_channel", "gibbs_ctrl_toffoli_per_rot", "gibbs_structure", "chain_reuse_gain",
      "chain_reuse_steps_per_tau", "eps_qsp", "R_warm_stated", "qsvt_rail_stated", "d_qsp_stated")),
)

TAG = {Provenance.CITED: "Cited", Provenance.ASSUMED: "Assumed",
       Provenance.STATED: "Stated-not-derived", Provenance.UNCITED: "Uncited"}


def _cell(x) -> str:
    return str(x).replace("|", "\\|").replace("\n", " ")


def main() -> int:
    a = Assumptions()
    names = [f.name for f in dataclasses.fields(a)]
    listed = [n for _, ns in GROUPS for n in ns]
    missing = [n for n in names if n not in listed]
    extra = [n for n in listed if n not in names]
    dup = sorted({n for n in listed if listed.count(n) > 1})
    if missing or extra or dup:
        raise SystemExit(f"grouping out of date: missing {missing}, unknown {extra}, duplicated {dup}")

    counts = {}
    for n in names:
        p = TAG[getattr(a, n).prov]
        counts[p] = counts.get(p, 0) + 1
    print(f"{len(names)} fields: " + ", ".join(f"{k} {v}" for k, v in sorted(counts.items())))
    for title, ns in GROUPS:
        print()
        print(f"#### {title}")
        print()
        print("| name | value | provenance | source | note |")
        print("|---|---|---|---|---|")
        for n in ns:
            v = getattr(a, n)
            print(f"| `{n}` | {_cell(v.value)} | {TAG[v.prov]} | {_cell(v.src)} | {_cell(v.note)} |")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
