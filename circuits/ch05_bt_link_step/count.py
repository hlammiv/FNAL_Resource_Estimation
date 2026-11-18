"""Counts of the 2T link-step circuits, obtained by walking the circuit objects of build.py, next to the
line items of the chapter model (estimates/ch05_qgp_transport.py and estimates/groups.py).

    python circuits/ch05_bt_link_step/count.py          # prints the tables and writes counts.json

Report convention (estimates/common.py, imported, not copied): 7 T per Toffoli (toffoli_t, "textbook"),
rotations at c_T = t_per_rotation(eps_rot) with eps_rot = eps_rot_for(N_rot), eps_syn = 1e-2 per shot.
T-depth convention: the chapter model's (ch05 SCHED comment): a Toffoli is one layer, a T gate one layer,
a synthesized rotation c_T layers, scheduled as soon as its qubits are free (build.weighted_depth).
"""
from __future__ import annotations

import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, HERE)
sys.path.insert(0, REPO)

import build as B  # noqa: E402
from estimates import ch05_qgp_transport as M  # noqa: E402
from estimates.common import EPS_SYN, eps_rot_for, t_per_rotation, toffoli_t  # noqa: E402
from estimates.groups import GROUPS, PRIMCOST  # noqa: E402

MAGNETIC_BLOCKS = ("U_inv", "U_mul", "U_Tr")
ELECTRIC_BLOCKS = ("U_F", "U_phase")


def report_t(n_toffoli: float, n_rot: float, n_t: float, c_t: float) -> float:
    return toffoli_t(n_toffoli, "textbook") + n_t + n_rot * c_t


def circuit_row(qc, n_links: int, steps: int) -> dict:
    """Everything counted from one circuit object."""
    tot = B.count_gates(qc)
    blocks = B.count_by_block(qc)
    eps = eps_rot_for(tot["rot"], EPS_SYN)
    c_t = t_per_rotation(eps)
    per = {}
    for name, c in blocks.items():
        per[name] = dict(calls=c["calls"], n_toffoli=c["toffoli"], n_rot=c["rot"], n_t_fixed=c["t"],
                         clifford=c["clifford"],
                         per_call=dict(n_toffoli=c["toffoli"] / c["calls"], n_rot=c["rot"] / c["calls"],
                                       n_t_fixed=c["t"] / c["calls"]),
                         calls_per_link_step=c["calls"] / (n_links * steps),
                         t=report_t(c["toffoli"], c["rot"], c["t"], c_t))
    t_total = report_t(tot["toffoli"], tot["rot"], tot["t"], c_t)
    mag = sum(per[b]["t"] for b in MAGNETIC_BLOCKS)
    ele = sum(per[b]["t"] for b in ELECTRIC_BLOCKS)
    assert abs(mag + ele - t_total) < 1e-6 * t_total
    return dict(
        qubits_declared=qc.num_qubits, logical_qubits=B.used_qubits(qc),
        data_qubits=B.LINK_QUBITS * n_links, ancilla_qubits=B.used_qubits(qc) - B.LINK_QUBITS * n_links,
        n_links=n_links, steps=steps, n_gates=len(qc.data),
        n_toffoli=tot["toffoli"], n_rot=tot["rot"], n_t_fixed=tot["t"], n_clifford=tot["clifford"],
        eps_rot=eps, c_t=c_t, t_total=t_total,
        t_per_link_step=t_total / (n_links * steps),
        t_magnetic_per_link_step=mag / (n_links * steps),
        t_electric_per_link_step=ele / (n_links * steps),
        n_toffoli_per_link_step=tot["toffoli"] / (n_links * steps),
        n_rot_per_link_step=tot["rot"] / (n_links * steps),
        n_t_fixed_per_link_step=tot["t"] / (n_links * steps),
        t_depth=B.weighted_depth(qc, 1.0, 1.0, c_t),
        non_clifford_layers_unit_rotation=B.weighted_depth(qc, 1.0, 1.0, 1.0),
        blocks=per,
    )


def model_values() -> dict:
    """The chapter model's line items, each from a named function or record of the model code."""
    a = M.Assumptions()
    gk = a.group_2028.value
    L, d = int(a.L_2028.lo), int(a.dim_2028.lo)
    prim = GROUPS[gk].primitives
    mult = {p: PRIMCOST[a.hamiltonian.value][p](d) for p in ("U_mul", "U_inv", "U_Tr", "U_F")}
    r = M.model(a, "2028")
    it = r.intermediates
    out = dict(
        source="estimates/ch05_qgp_transport.py (Assumptions, evolution_2028, step_depth_2028, model(a,'2028')"
               ".intermediates) and estimates/groups.py (GROUPS['2T'].primitives, PRIMCOST['KS'])",
        L=L, dim=d, n_links=M.n_links(L, d), n_plaquettes=M.n_plaquettes(L, d),
        dt_first_fm=a.dt_first_2028_fm.lo, dt_full_fm=a.dt_over_at_2028.lo * a.a_t_fm.lo, a_s_fm=a.a_s_fm.lo,
        lq=list(r.lq), gauge_qubits=it["gauge_qubits_2028"], workspace=it["workspace_2028"],
        hadamard_ancilla=list(a.anc_2028.value),
        primitives={p: dict(n_toffoli=prim[k].toffoli, n_rot=prim[k].n_rot, ancilla=prim[k].ancilla,
                            calls_per_link_step=mult[p], src=prim[k].src)
                    for p, k in (("U_mul", "U_mul"), ("U_inv", "U_inv"), ("U_Tr", "U_Tr"), ("U_F", "U_FFT"))},
        n_rot_per_link_step=it["n_rot_per_link_step_2028"],
        box_t_per_rotation=it["t_per_rotation_report"],
        box_magnetic_t_per_link_step=it["magnetic_t_per_link_step_2028"],
        box_electric_t_per_link_step=it["electric_t_per_link_step_2028"],
        box_t_per_link_step=it["t_per_link_step_2028"],
    )
    out["n_toffoli_per_link_step"] = sum(mult[p] * out["primitives"][p]["n_toffoli"] for p in mult)
    for tier, steps in (("first", int(a.first_result_steps_2028.lo)), ("full", int(a.trotter_steps_2028.lo))):
        e = M.evolution_2028(a, steps)
        pc = e["pc"]
        depth = M.step_depth_2028(a, e["t_rot"])
        out[tier] = dict(steps=steps, n_rot=e["n_rot"], eps_rot=e["eps"], c_t=e["t_rot"], t_total=e["t"],
                         n_toffoli=out["n_toffoli_per_link_step"] * e["links"] * steps,
                         t_per_link_step=pc["per_link"], t_magnetic_per_link_step=pc["magnetic_per_link"],
                         t_electric_per_link_step=pc["electric_per_link"],
                         t_depth_low=steps * depth[0], t_depth_high=steps * depth[1])
    return out


def ratio(x, y):
    return None if not y else x / y


def main() -> dict:
    model = model_values()
    b = B.BENCHMARK
    # the generator defaults must be the model's instance
    assert (b["L"], b["a_s_fm"], b["dt_first_fm"], b["n_workspace"]) == (
        model["L"], model["a_s_fm"], model["dt_first_fm"], model["workspace"])
    assert abs(b["dt_full_fm"] - model["dt_full_fm"]) < 1e-12
    assert (b["steps_first"], b["steps_full"]) == (model["first"]["steps"], model["full"]["steps"])
    n_links = model["n_links"]

    circuit = {}
    for tier in ("first", "full"):
        circuit[tier] = circuit_row(B.benchmark_circuit(tier), n_links, model[tier]["steps"])
    variants = {
        "first, uncompute=dagger": circuit_row(B.benchmark_circuit("first", uncompute="dagger"), n_links, 1),
        "first, order=2 (two magnetic half layers)": circuit_row(B.benchmark_circuit("first", order=2), n_links, 1),
        "first, 2 workspace qubits": circuit_row(B.benchmark_circuit("first", n_workspace=2), n_links, 1),
        "first, 64 workspace qubits": circuit_row(B.benchmark_circuit("first", n_workspace=64), n_links, 1),
        "full, uncompute=dagger": circuit_row(B.benchmark_circuit("full", uncompute="dagger"), n_links, 6),
    }
    te, tb = B.couplings(B.G2_PLACEHOLDER, b["dt_first_fm"] / b["a_s_fm"])
    plaq = circuit_row(B.single_plaquette_step(te, tb), 4, 1)

    # the circuit counted at the model's own T per rotation (isolates the gate counts from the tolerance)
    at_model_ct = {}
    for tier in ("first", "full"):
        c, ct = circuit[tier], model[tier]["c_t"]
        mag = sum(report_t(c["blocks"][k]["n_toffoli"], c["blocks"][k]["n_rot"], c["blocks"][k]["n_t_fixed"], ct)
                  for k in MAGNETIC_BLOCKS)
        ele = sum(report_t(c["blocks"][k]["n_toffoli"], c["blocks"][k]["n_rot"], c["blocks"][k]["n_t_fixed"], ct)
                  for k in ELECTRIC_BLOCKS)
        at_model_ct[tier] = dict(c_t=ct, t_total=mag + ele, t_magnetic_per_link_step=mag / (n_links * c["steps"]),
                                 t_electric_per_link_step=ele / (n_links * c["steps"]))

    ratios = {}
    for tier in ("first", "full"):
        c, m = circuit[tier], model[tier]
        ratios[tier] = dict(
            t_total=ratio(c["t_total"], m["t_total"]),
            t_per_link_step=ratio(c["t_per_link_step"], m["t_per_link_step"]),
            t_magnetic_per_link_step=ratio(c["t_magnetic_per_link_step"], m["t_magnetic_per_link_step"]),
            t_electric_per_link_step=ratio(c["t_electric_per_link_step"], m["t_electric_per_link_step"]),
            n_toffoli=ratio(c["n_toffoli"], m["n_toffoli"]),
            n_rot=ratio(c["n_rot"], m["n_rot"]),
            c_t=ratio(c["c_t"], m["c_t"]),
            t_total_at_model_c_t=ratio(at_model_ct[tier]["t_total"], m["t_total"]),
            t_depth_vs_model_low=ratio(c["t_depth"], m["t_depth_low"]),
            t_depth_vs_model_high=ratio(c["t_depth"], m["t_depth_high"]),
            logical_qubits_vs_model_lq_hi=ratio(c["logical_qubits"], model["lq"][1]),
        )
    prim_ratio = {}
    cb = circuit["first"]["blocks"]
    for p in ("U_mul", "U_inv", "U_Tr", "U_F"):
        mp, cp = model["primitives"][p], cb[p]
        prim_ratio[p] = dict(
            n_toffoli=ratio(cp["per_call"]["n_toffoli"], mp["n_toffoli"]),
            n_rot=ratio(cp["per_call"]["n_rot"], mp["n_rot"]),
            calls_per_link_step=ratio(cp["calls_per_link_step"], mp["calls_per_link_step"]))

    out = dict(
        instance="Ch. 5 2028 benchmark: 2+1D SU(2) -> 2T, 4x4 periodic lattice, H_KS, pure gauge; "
                 "first = one step, full = six steps; order-1 step (one electric + one magnetic layer), "
                 "plaquette uncompute as arXiv:1903.08807 Fig. 1, 16 workspace qubits",
        convention=dict(t_per_toffoli=toffoli_t(1, "textbook"), eps_syn=EPS_SYN,
                        c_t="estimates.common.t_per_rotation(eps_rot_for(N_rot))",
                        t_depth="Toffoli 1 layer, T 1 layer, rotation c_T layers, ASAP on the circuit's qubits"),
        model=model, circuit=circuit, circuit_at_model_c_t=at_model_ct, ratios_circuit_over_model=ratios,
        primitive_ratios_circuit_over_model=prim_ratio, variants=variants, single_plaquette_step=plaq,
        g2_placeholder=B.G2_PLACEHOLDER,
    )
    with open(os.path.join(HERE, "counts.json"), "w") as fh:
        json.dump(out, fh, indent=1, default=float)
    return out


def _fmt(x):
    if x is None:
        return "-"
    if isinstance(x, float) and x != int(x):
        return f"{x:.4g}"
    return f"{int(x)}"


if __name__ == "__main__":
    o = main()
    m, c, r = o["model"], o["circuit"], o["ratios_circuit_over_model"]
    print("PRIMITIVES (per call; calls per link-step)")
    print(f"{'':10s}{'model Tof':>10s}{'rot':>6s}{'anc':>5s}{'calls':>7s} | {'circ Tof':>9s}{'rot':>6s}{'T':>5s}{'calls':>7s}")
    cb = c["first"]["blocks"]
    for p in ("U_mul", "U_inv", "U_Tr", "U_F", "U_phase"):
        mp = m["primitives"].get(p)
        cp = cb[p]
        left = (f"{_fmt(mp['n_toffoli']):>10s}{_fmt(mp['n_rot']):>6s}{_fmt(mp['ancilla']):>5s}"
                f"{_fmt(mp['calls_per_link_step']):>7s}") if mp else f"{'-':>10s}{'-':>6s}{'-':>5s}{'-':>7s}"
        print(f"{p:10s}{left} | {_fmt(cp['per_call']['n_toffoli']):>9s}{_fmt(cp['per_call']['n_rot']):>6s}"
              f"{_fmt(cp['per_call']['n_t_fixed']):>5s}{_fmt(cp['calls_per_link_step']):>7s}")
    for tier in ("first", "full"):
        print(f"\n{tier.upper()} ({c[tier]['steps']} step(s), {c[tier]['n_links']} links)")
        rows = [("logical qubits", c[tier]["logical_qubits"], m["lq"][1], r[tier]["logical_qubits_vs_model_lq_hi"]),
                ("N_Toffoli", c[tier]["n_toffoli"], m[tier]["n_toffoli"], r[tier]["n_toffoli"]),
                ("N_rot", c[tier]["n_rot"], m[tier]["n_rot"], r[tier]["n_rot"]),
                ("N_T_fixed", c[tier]["n_t_fixed"], 0, None),
                ("eps_rot", c[tier]["eps_rot"], m[tier]["eps_rot"], None),
                ("c_T per rotation", c[tier]["c_t"], m[tier]["c_t"], r[tier]["c_t"]),
                ("T magnetic / link-step", c[tier]["t_magnetic_per_link_step"], m[tier]["t_magnetic_per_link_step"],
                 r[tier]["t_magnetic_per_link_step"]),
                ("T electric / link-step", c[tier]["t_electric_per_link_step"], m[tier]["t_electric_per_link_step"],
                 r[tier]["t_electric_per_link_step"]),
                ("T / link-step", c[tier]["t_per_link_step"], m[tier]["t_per_link_step"], r[tier]["t_per_link_step"]),
                ("T total", c[tier]["t_total"], m[tier]["t_total"], r[tier]["t_total"]),
                ("T total at model c_T", o["circuit_at_model_c_t"][tier]["t_total"], m[tier]["t_total"],
                 r[tier]["t_total_at_model_c_t"]),
                ("T-depth (model low)", c[tier]["t_depth"], m[tier]["t_depth_low"], r[tier]["t_depth_vs_model_low"]),
                ("T-depth (model high)", c[tier]["t_depth"], m[tier]["t_depth_high"], r[tier]["t_depth_vs_model_high"])]
        print(f"{'':26s}{'circuit':>12s}{'model':>12s}{'ratio':>8s}")
        for name, x, y, q in rows:
            print(f"{name:26s}{_fmt(float(x)):>12s}{_fmt(float(y)):>12s}{_fmt(q):>8s}")
    print("\nVARIANTS")
    for name, v in o["variants"].items():
        print(f"{name:44s} LQ {v['logical_qubits']:4d}  Tof {v['n_toffoli']:6d}  rot {v['n_rot']:5d}  "
              f"T {v['n_t_fixed']:5d}  T_total {v['t_total']:.4g}  T-depth {v['t_depth']:.4g}")
    v = o["single_plaquette_step"]
    print(f"{'single plaquette step (22 qubits)':44s} LQ {v['logical_qubits']:4d}  Tof {v['n_toffoli']:6d}  "
          f"rot {v['n_rot']:5d}  T {v['n_t_fixed']:5d}  T_total {v['t_total']:.4g}  T-depth {v['t_depth']:.4g}")
