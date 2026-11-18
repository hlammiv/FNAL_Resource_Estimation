#!/usr/bin/env python3
"""Walk the Ch. 7 2028 step circuit and count it, per block and in total, next to the model's line items.

Usage:  python count.py            writes counts.json and the OpenQASM 3 export next to this file
        python count.py --no-qasm  counts only

Counting rules (estimates/common.py, imported, not copied):
  Toffoli            T_PER_TOFFOLI["textbook"] = 7 T (none occur in this circuit)
  arbitrary Rz       t_per_rotation(eps_rot_for(N_rot, EPS_SYN)): 1.15 log2(1/eps) + 9.2 T, eps = sqrt(1e-2/N_rot)
  Rz at k pi/4       exact: k even is Clifford, k odd is one T (the builder emits these as S/T gates already)
  T, Tdg             one T each
Depth: longest path through the circuit counting only non-Clifford gates. Three weightings are reported:
  nonclifford_layers   every non-Clifford gate counts 1
  rotation_layers      arbitrary rotations count 1, exact T gates 0
  t_depth_report       arbitrary rotations count t_per_rotation, exact T gates 1, Toffoli 1
The last is the model's convention where it defines one (estimates/ch07_curved_space.py, 2033 note:
"rotation layers per step x T/rot x steps"). The model's 2028 T-depth (346-1308) is an input taken from an
audit file that is not in the repository, so it is quoted for reference and not reproduced.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from collections import OrderedDict
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
for p in (str(REPO), str(HERE)):
    if p not in sys.path:
        sys.path.insert(0, p)

from estimates.common import EPS_SYN, T_PER_TOFFOLI, eps_rot_for, t_per_rotation, toffoli_t   # noqa: E402
from estimates import ch07_curved_space as model                                               # noqa: E402
import build as B                                                                              # noqa: E402

CLIFFORD = {"h", "cx", "cz", "s", "sdg", "z", "x", "y", "swap", "id", "barrier"}
T_GATES = {"t", "tdg"}
TOFFOLI = {"ccx", "ccz"}
QASM_NAME = "ch07_scalar_ds_step_4x4x4_K8.qasm"
WIDTH_CAPS = (10, 30, 50)        # 10 factories of the report baseline; N_anc = 30-50 of the chapter


def classify(op) -> str:
    """'clifford' | 't' | 'rot' | 'toffoli' for one instruction."""
    name = op.name
    if name in CLIFFORD:
        return "clifford"
    if name in T_GATES:
        return "t"
    if name in TOFFOLI:
        return "toffoli"
    if name in ("rz", "p"):
        k = float(op.params[0]) / (math.pi / 4)
        if abs(k - round(k)) < 1e-9:
            return "t" if int(round(k)) % 2 else "clifford"
        return "rot"
    raise ValueError(f"gate {name!r} is outside the logical gate set of this circuit")


def t_total(n_rot: int, n_t_fixed: int, n_toffoli: int) -> dict:
    """T count under the report convention, from estimates.common."""
    eps = eps_rot_for(n_rot, EPS_SYN) if n_rot else float("nan")
    t_rot = t_per_rotation(eps) if n_rot else 0.0
    return dict(eps_rot=eps, t_per_rotation=t_rot,
                t_total=n_rot * t_rot + n_t_fixed + toffoli_t(n_toffoli, "textbook"))


def walk(built: B.Built) -> dict:
    """Counts by block, totals, and depths, from the instruction list."""
    qc = built.qc
    blocks: "OrderedDict[str, dict]" = OrderedDict()
    tot = dict(gates=0, clifford=0, cx=0, n_rot=0, n_t_fixed=0, n_toffoli=0)
    kinds = []
    for ins, lab in zip(qc.data, built.labels):
        kind = classify(ins.operation)
        kinds.append(kind)
        b = blocks.setdefault(lab, dict(gates=0, clifford=0, cx=0, n_rot=0, n_t_fixed=0, n_toffoli=0))
        for d in (b, tot):
            d["gates"] += 1
            d["cx"] += ins.operation.name == "cx"
            d["clifford"] += kind == "clifford"
            d["n_rot"] += kind == "rot"
            d["n_t_fixed"] += kind == "t"
            d["n_toffoli"] += kind == "toffoli"
    conv = t_total(tot["n_rot"], tot["n_t_fixed"], tot["n_toffoli"])
    t_rot = conv["t_per_rotation"]
    for b in blocks.values():
        b["t_report"] = b["n_rot"] * t_rot + b["n_t_fixed"] + toffoli_t(b["n_toffoli"], "textbook")
    weights = {"nonclifford_layers": dict(rot=1, t=1, toffoli=1),
               "rotation_layers": dict(rot=1, t=0, toffoli=0),
               "t_depth_report": dict(rot=t_rot, t=1, toffoli=1)}
    depth = {}
    busiest = None
    for key, w in weights.items():
        d = [0.0] * qc.num_qubits
        per_layer = {}
        for ins, kind in zip(qc.data, kinds):
            qs = [qc.find_bit(q).index for q in ins.qubits]
            v = max(d[q] for q in qs) + w.get(kind, 0)
            for q in qs:
                d[q] = v
            if key == "rotation_layers" and kind == "rot":
                per_layer[v] = per_layer.get(v, 0) + 1
        depth[key] = max(d)
        if key == "rotation_layers":
            depth["peak_rotations_in_one_layer"] = max(per_layer.values()) if per_layer else 0
    # rotation layers when at most `cap` rotations are in flight at once (one synthesis ancilla each,
    # estimates.common.RUS_ANCILLA): greedy list schedule in circuit order, a valid schedule, not an optimum
    depth["rotation_layers_width_capped"] = {}
    for cap in WIDTH_CAPS:
        d = [0] * qc.num_qubits
        load = {}
        for ins, kind in zip(qc.data, kinds):
            qs = [qc.find_bit(q).index for q in ins.qubits]
            v = max(d[q] for q in qs)
            if kind == "rot":
                v += 1
                while load.get(v, 0) >= cap:
                    v += 1
                load[v] = load.get(v, 0) + 1
            for q in qs:
                d[q] = v
        depth["rotation_layers_width_capped"][str(cap)] = max(d)
    # rotations acting on the busiest qubit (the model's "rotation layers on its busiest qubit" counts the
    # rotations of every string a qubit belongs to)
    per_qubit = [0] * qc.num_qubits
    for i, (ins, kind) in enumerate(zip(qc.data, kinds)):
        if kind != "rot":
            continue
        q = qc.find_bit(ins.qubits[0]).index
        per_qubit[q] += 1
        # a ZZ rotation sits between two CX gates with the same control: charge the control too
        prev = qc.data[i - 1] if i else None
        if prev is not None and prev.operation.name == "cx" and qc.find_bit(prev.qubits[1]).index == q \
                and i + 1 < len(qc.data) and qc.data[i + 1].operation.name == "cx" \
                and [qc.find_bit(x).index for x in qc.data[i + 1].qubits] == [qc.find_bit(x).index for x in prev.qubits]:
            per_qubit[qc.find_bit(prev.qubits[0]).index] += 1
    busiest = max(per_qubit)
    # ancilla = qubits of the circuit beyond the field register of the instance. Counted, not assumed.
    n_data = built.params.inst.n_qubits
    touched = {qc.find_bit(q).index for ins in qc.data for q in ins.qubits}
    return dict(qubits=qc.num_qubits, ancilla=qc.num_qubits - n_data, qubits_acted_on=len(touched),
                blocks=blocks, totals={**tot, **conv}, depth=depth,
                rotation_strings_on_busiest_qubit=busiest)


def model_items() -> dict:
    """The model's 2028 line items, each with the function or key it was read from."""
    a = model.Assumptions()
    sc = model.shot_counts(a, "2028", a.n_steps_2028.lo)
    r = model.model(a, "2028")
    reg = model.register(a, "2028")
    it = r.intermediates
    return dict(
        source="estimates/ch07_curved_space.py: shot_counts(Assumptions(), '2028', 1), register(), model(a, '2028')",
        system_qubits=reg["system_qubits"], n_anc=list(reg["n_anc"]), lq=list(reg["lq"]),
        rot_diag=sc["rot"]["diag"], rot_pi2=sc["rot"]["pi2"], rot_qft=sc["rot"]["qft"],
        rot_centering=sc["rot"]["centering"], rot_outer_shear=0,
        t_controlled_s=sc["t_exact"]["qft_cs"], t_centering=sc["t_exact"]["centering"],
        n_rot=sc["n_rot"], n_t_fixed=sum(sc["t_exact"].values()), n_toffoli=0,
        eps_rot=sc["eps_rot"], t_per_rotation=sc["t_per_rot"],
        t_step=sc["t_total"], t_evolution=it["t_evolution"], t_prep_stated=it["t_prep"],
        t_per_shot=r.hard_ops[0], hard_op_cap=it["hard_op_cap"],
        t_depth_per_shot=list(a.t_depth_2028.value), t_depth_source=a.t_depth_2028.src,
    )


def _sum(blocks, names, key):
    return sum(blocks[n][key] for n in names if n in blocks)


def line_items(w: dict, m: dict) -> list:
    """Circuit count next to the model line item, with the ratio (None where the model has 0)."""
    bl = w["blocks"]
    kin = [n for n in bl if n.endswith("_pi2")]
    qft = [n for n in bl if n.endswith("_qft") or n.endswith("_qft_inv")]
    cen = [n for n in bl if "centering" in n]
    diag = [n for n in bl if n.startswith("d") and (n.endswith("_site") or n.endswith("_link"))]
    shear = [n for n in bl if n.startswith("shear_")]
    tot = w["totals"]
    rows = [
        ("logical qubits used by the circuit (data + ancilla)", w["qubits"], m["system_qubits"],
         "model: system qubits, register()['system_qubits']"),
        ("logical qubits, model total with N_anc = 30-50", w["qubits"], m["lq"][0],
         "model: register()['lq'] low end; high end %d" % m["lq"][1]),
        ("rotations: diagonal block D (site x^2 + link)", _sum(bl, diag, "n_rot"), m["rot_diag"],
         "shot_counts()['rot']['diag']"),
        ("rotations: pi^2 phases (two kinetic blocks)", _sum(bl, kin, "n_rot"), m["rot_pi2"],
         "shot_counts()['rot']['pi2']"),
        ("rotations: QFT controlled-R_3 (4 QFTs)", _sum(bl, qft, "n_rot"), m["rot_qft"],
         "shot_counts()['rot']['qft']"),
        ("rotations: centering layers", _sum(bl, cen, "n_rot"), m["rot_centering"],
         "shot_counts()['rot']['centering']"),
        ("rotations: frame shear at the two ends of the shot", _sum(bl, shear, "n_rot"), m["rot_outer_shear"],
         "not a model line item"),
        ("N_rot (arbitrary rotations)", tot["n_rot"], m["n_rot"], "shot_counts()['n_rot']"),
        ("T: QFT controlled-S", _sum(bl, qft, "n_t_fixed"), m["t_controlled_s"],
         "shot_counts()['t_exact']['qft_cs']"),
        ("T: centering layers", _sum(bl, cen, "n_t_fixed"), m["t_centering"],
         "shot_counts()['t_exact']['centering']"),
        ("N_T_fixed (exact T gates)", tot["n_t_fixed"], m["n_t_fixed"], "sum of shot_counts()['t_exact']"),
        ("N_Toff", tot["n_toffoli"], m["n_toffoli"], "the model has no Toffoli"),
        ("eps_rot", tot["eps_rot"], m["eps_rot"], "shot_counts()['eps_rot']"),
        ("T per rotation", tot["t_per_rotation"], m["t_per_rotation"], "shot_counts()['t_per_rot']"),
        ("T, one step (report convention)", tot["t_total"], m["t_step"], "shot_counts()['t_total']"),
    ]
    out = []
    for name, c, mv, src in rows:
        out.append(dict(item=name, circuit=c, model=mv, ratio=(c / mv if mv else None), model_source=src))
    return out


VARIANTS = OrderedDict([
    ("full", dict(in_shear=True, out_shear=True, merged_rz=False,
                  note="headline: the exact digitized one-step map, textbook gate decompositions")),
    ("in_shear_only", dict(in_shear=True, out_shear=False, merged_rz=False,
                           note="last shear layer dropped (diagonal before a field-basis readout), first kept")),
    ("no_outer_shear", dict(in_shear=False, out_shear=False, merged_rz=False,
                            note="frame shear at the two ends left out (the model's accounting)")),
    ("full_merged_rz", dict(in_shear=True, out_shear=True, merged_rz=True,
                            note="headline circuit after the Rz merge pass build.merge_rz")),
    ("no_outer_shear_merged_rz", dict(in_shear=False, out_shear=False, merged_rz=True,
                                      note="no outer shear, after the Rz merge pass")),
])


def build_variant(inst: B.Instance, name: str) -> B.Built:
    v = VARIANTS[name]
    bt = B.build_step(inst, in_shear=v["in_shear"], out_shear=v["out_shear"])
    return B.merge_rz(bt) if v["merged_rz"] else bt


def main(write_qasm: bool = True) -> dict:
    inst = B.BENCHMARK_2028
    m = model_items()
    out = OrderedDict()
    out["instance"] = dict(dims=list(inst.dims), n_q=inst.n_q, K=inst.K, m_over_Hinf=inst.m, h_Hinf=inst.h,
                           n_steps=inst.n_steps, lattice_spacing_b=inst.spacing,
                           generator="build.build_step(build.BENCHMARK_2028)")
    out["conventions"] = dict(t_per_toffoli=T_PER_TOFFOLI["textbook"], eps_syn=EPS_SYN,
                              rotation_cost="estimates.common.t_per_rotation(eps_rot_for(N_rot))",
                              depth="longest non-Clifford path; t_depth_report weights rotations by T per rotation")
    out["model"] = m
    out["variants"] = OrderedDict()
    for name, v in VARIANTS.items():
        bt = build_variant(inst, name)
        w = walk(bt)
        tot = w["totals"]
        out["variants"][name] = dict(
            note=v["note"], qubits=w["qubits"], ancilla=w["ancilla"],
            n_toffoli=tot["n_toffoli"], n_rot=tot["n_rot"], n_t_fixed=tot["n_t_fixed"],
            cx=tot["cx"], gates=tot["gates"], eps_rot=tot["eps_rot"], t_per_rotation=tot["t_per_rotation"],
            t_total=tot["t_total"], t_total_over_model_step=tot["t_total"] / m["t_step"],
            t_total_plus_stated_prep=tot["t_total"] + m["t_prep_stated"],
            fits_hard_op_cap_with_stated_prep=bool(tot["t_total"] + m["t_prep_stated"] <= m["hard_op_cap"]),
            depth=w["depth"],
            t_depth_width_capped_estimate={k: v * tot["t_per_rotation"]
                                           for k, v in w["depth"]["rotation_layers_width_capped"].items()},
            rotation_strings_on_busiest_qubit=w["rotation_strings_on_busiest_qubit"],
            blocks=w["blocks"], line_items=line_items(w, m))
        if name == "full":
            out["block_ranges_full"] = [dict(label=l, start=s, stop=e) for l, s, e in bt.block_ranges()]
            if write_qasm:
                from qiskit import qasm3
                text = qasm3.dumps(bt.qc)
                (HERE / QASM_NAME).write_text(text)
                out["qasm"] = dict(file=QASM_NAME, bytes=len(text.encode()), instructions=len(bt.qc.data))
    (HERE / "counts.json").write_text(json.dumps(out, indent=1, default=float) + "\n")
    return out


def _print(out: dict):
    m = out["model"]
    for name, v in out["variants"].items():
        print(f"\n[{name}] {v['note']}")
        print(f"  qubits {v['qubits']}  N_Toff {v['n_toffoli']}  N_rot {v['n_rot']}  N_T_fixed {v['n_t_fixed']}  "
              f"T/rot {v['t_per_rotation']:.3f}  T step {v['t_total']:.1f}  (model {m['t_step']:.1f}, "
              f"ratio {v['t_total_over_model_step']:.3f})")
        print(f"  depth {v['depth']}  busiest qubit: {v['rotation_strings_on_busiest_qubit']} rotation strings")
        for row in v["line_items"]:
            r = "  n/a" if row["ratio"] is None else f"{row['ratio']:5.3f}"
            print(f"    {row['item']:<55s} circuit {row['circuit']:<12.6g} model {row['model']:<12.6g} ratio {r}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-qasm", action="store_true")
    args = ap.parse_args()
    _print(main(write_qasm=not args.no_qasm))
