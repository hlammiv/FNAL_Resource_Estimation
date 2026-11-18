"""Count the Ch. 4 2028 benchmark circuit by walking the circuit object, and set the counts next to the
line items of the chapter model (estimates/ch04_hybrid_lqcd.py, run here, never copied).

    python count.py            # prints the tables and writes counts.json
    python count.py --qasm     # also writes the OpenQASM 3 files

Conventions (estimates/common.py, imported):
    Toffoli            7 T                      common.toffoli_t(n, "textbook")
    arbitrary rz       1.15 log2(1/eps) + 9.2   common.t_per_rotation(eps)
    eps_rot            sqrt(eps_syn / N_rot)    common.eps_rot_for(N_rot, eps_syn), N_rot from the circuit
    eps_syn            the chapter model's departure, 0.1 x the polynomial tolerance (model intermediate
                       'eps_syn'), and the report default common.EPS_SYN for reference
An rz whose angle is a multiple of pi/4 is exact: Clifford for even multiples, one T for odd ones.

T-depth. Weighted longest path through the circuit. Weights follow the model (its 't_depth_per_shot'):
one layer per Toffoli, t_rot layers per arbitrary rotation, one layer per exact T. Two orderings:
    commuting    a gate waits only for earlier gates it does not commute with. Two gates commute here when,
                 on every qubit they share, both are diagonal (controls, rz, s, t) or both are X-type
                 (targets of cx/ccx, x). This is the model's convention: O_c runs under the O_A rotations.
    as_written   a gate waits for every earlier gate on any of its qubits.
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
for p in (str(REPO), str(HERE)):          # HERE goes in last, so the modules of this directory are found first
    if p in sys.path:
        sys.path.remove(p)
    sys.path.insert(0, p)

import build                                                                        # noqa: E402
from estimates.common import EPS_SYN, eps_rot_for, t_per_rotation, toffoli_t        # noqa: E402
from estimates.ch04_hybrid_lqcd import Assumptions, PUBLISHED, model                # noqa: E402

TOFFOLI_CONVENTION = "textbook"
CLIFFORD_1Q = {"h", "x", "y", "z", "s", "sdg"}
IGNORED = {"measure", "barrier"}
BLOCK_ORDER = ("prep", "had", "D_s", "O_A", "O_c", "proj_and", "proj_rot")


def rz_class(theta: float) -> str:
    """'clifford', 't' or 'rot' for an rz angle."""
    k = theta / (math.pi / 4.0)
    if abs(k - round(k)) < 1e-9:
        return "clifford" if round(k) % 2 == 0 else "t"
    return "rot"


def gate_class(inst) -> str:
    """One of 'toffoli', 'rot', 't', 'cx', 'clifford_1q', 'ignored'. Anything else is an error, so
    no gate can go uncounted."""
    name = inst.operation.name
    if name == "ccx":
        return "toffoli"
    if name == "rz":
        c = rz_class(float(inst.operation.params[0]))
        return {"rot": "rot", "t": "t", "clifford": "clifford_1q"}[c]
    if name in ("t", "tdg"):
        return "t"
    if name == "cx":
        return "cx"
    if name in CLIFFORD_1Q:
        return "clifford_1q"
    if name in IGNORED:
        return "ignored"
    raise ValueError(f"gate {name!r} is outside the counted gate set")


def walk_counts(qc) -> dict:
    """Totals and per-block counts, from the circuit data and the block marks."""
    keys = ("toffoli", "rot", "t", "cx", "clifford_1q")
    owner = [None] * len(qc.data)
    instances: dict = {}
    for name, start, stop in build.block_marks(qc):
        instances[name] = instances.get(name, 0) + 1
        for i in range(start, stop):
            if owner[i] is not None:
                raise ValueError("overlapping block marks")
            owner[i] = name
    blocks = {name: dict.fromkeys(keys, 0) for name in instances}
    total = dict.fromkeys(keys, 0)
    used = set()
    for i, inst in enumerate(qc.data):
        c = gate_class(inst)
        if c == "ignored":
            continue
        if owner[i] is None:
            raise ValueError(f"instruction {i} ({inst.operation.name}) is in no block")
        blocks[owner[i]][c] += 1
        total[c] += 1
        used.update(qc.find_bit(q).index for q in inst.qubits)
    for name in blocks:
        blocks[name]["marked_segments"] = instances[name]
    per_reg = {}
    for reg in qc.qregs:
        idx = [qc.find_bit(q).index for q in reg]
        per_reg[reg.name] = {"allocated": len(idx), "touched_by_a_gate": sum(i in used for i in idx)}
    return {"total": total, "blocks": blocks, "qubits_allocated": qc.num_qubits,
            "qubits_touched_by_a_gate": len(used), "registers": per_reg}


# qubit roles for the commutation rule: 'z' diagonal on the qubit, 'x' X-type on it, 'o' anything else
def _roles(inst) -> list:
    name = inst.operation.name
    n = len(inst.qubits)
    if name in ("rz", "s", "sdg", "t", "tdg", "z"):
        return ["z"]
    if name == "x":
        return ["x"]
    if name == "cx":
        return ["z", "x"]
    if name == "ccx":
        return ["z", "z", "x"]
    return ["o"] * n


def t_depth(qc, t_rot: float, ordering: str = "commuting") -> dict:
    """Weighted longest path. Returns the depth and what lies on one longest path."""
    if ordering not in ("commuting", "as_written"):
        raise ValueError(ordering)
    weight = {"toffoli": 1.0, "rot": float(t_rot), "t": 1.0}
    n = qc.num_qubits
    zero = (0.0, 0, 0, 0)                                    # (time, toffolis, rotations, exact T) on the path
    fin = {r: [zero] * n for r in "zxo"}
    best = zero
    for inst in qc.data:
        c = gate_class(inst)
        if c == "ignored":
            continue
        roles = _roles(inst) if ordering == "commuting" else ["o"] * len(inst.qubits)
        qs = [qc.find_bit(q).index for q in inst.qubits]
        start = zero
        for q, r in zip(qs, roles):
            waits = {"z": "xo", "x": "zo", "o": "zxo"}[r]
            for w in waits:
                if fin[w][q][0] > start[0]:
                    start = fin[w][q]
        end = (start[0] + weight.get(c, 0.0), start[1] + (c == "toffoli"), start[2] + (c == "rot"),
               start[3] + (c == "t"))
        for q, r in zip(qs, roles):
            if end[0] >= fin[r][q][0]:
                fin[r][q] = end
        if end[0] > best[0]:
            best = end
    return {"t_depth": best[0], "toffoli_layers_on_path": best[1], "rotations_on_path": best[2],
            "exact_t_on_path": best[3]}


def t_total(counts: dict, eps_syn: float) -> dict:
    """T under the report convention for one circuit, with its own eps_rot (common.eps_rot_for)."""
    n_rot, n_tof, n_t = counts["rot"], counts["toffoli"], counts["t"]
    eps = eps_rot_for(n_rot, eps_syn) if n_rot else float("nan")
    t_rot = t_per_rotation(eps) if n_rot else 0.0
    return {"eps_syn": eps_syn, "eps_rot": eps, "t_per_rotation": t_rot,
            "t_toffoli": toffoli_t(n_tof, TOFFOLI_CONVENTION), "t_rotations": n_rot * t_rot, "t_fixed": n_t,
            "t_total": toffoli_t(n_tof, TOFFOLI_CONVENTION) + n_rot * t_rot + n_t}


def model_items() -> dict:
    """The chapter model's 2028 line items, read from a live run of estimates.ch04_hybrid_lqcd.model."""
    r = model(Assumptions(), "2028")
    i = r.intermediates
    d, n_phase = i["be_queries"], i["projector_phases"]
    src = "estimates.ch04_hybrid_lqcd.model(Assumptions(), '2028')"
    ends = {}
    for k, variant in enumerate(("merged", "literal")):       # the model's band: index 0 merged, 1 literal
        ends[variant] = {
            "logical_qubits": (r.lq[k], f"{src}.lq"),
            "O_c_toffoli": (d * i["oc_toffolis_per_query"][k], "intermediates be_queries x oc_toffolis_per_query"),
            "O_A_rotations": (d * i["oa_rotations_per_query"][k], "intermediates be_queries x oa_rotations_per_query"),
            "projector_toffoli": (n_phase * i["projector_toffolis_each"],
                                  "intermediates projector_phases x projector_toffolis_each"),
            "projector_rotations": (n_phase * i["projector_rotations_each"],
                                    "intermediates projector_phases x projector_rotations_each"),
            "n_toffoli": (d * i["oc_toffolis_per_query"][k] + n_phase * i["projector_toffolis_each"],
                          "sum of the two Toffoli items"),
            "n_rot": (i["rotations_per_shot"][k], "intermediates rotations_per_shot"),
            "t_per_rotation": (i["t_per_rotation"][k], "intermediates t_per_rotation"),
            "t_total": (r.hard_ops[k], f"{src}.hard_ops"),
            "t_depth": (i["t_depth_per_shot"][k], "intermediates t_depth_per_shot"),
        }
    return {"ends": ends, "eps_syn": i["eps_syn"], "d": d, "projector_phases": n_phase,
            "published_lq": PUBLISHED["2028"].lq, "published_hard_ops": PUBLISHED["2028"].hard_ops,
            "poly_tolerance": i["poly_tolerance"], "logdet_W": i["logdet_W"],
            "normalized_trace_shifted": i["normalized_trace_shifted"],
            "breakdown_literal": [{"name": p.name, "count": p.count, "t_each": p.t_each, "t_total": p.t_total,
                                   "status": p.status.value} for p in r.breakdown]}


def count_variant(variant: str, mdl: dict) -> dict:
    qc = build.benchmark_circuit(variant)
    w = walk_counts(qc)
    tot, blk = w["total"], w["blocks"]
    t_model_budget = t_total(tot, mdl["eps_syn"])
    t_default_budget = t_total(tot, EPS_SYN)
    depth = {o: t_depth(qc, t_model_budget["t_per_rotation"], o) for o in ("commuting", "as_written")}
    out = {
        "gates": len(qc.data), "ops": {k: int(v) for k, v in qc.count_ops().items()},
        "logical_qubits": w["qubits_allocated"], "qubits_touched_by_a_gate": w["qubits_touched_by_a_gate"],
        "registers": w["registers"],
        "n_toffoli": tot["toffoli"], "n_rot": tot["rot"], "n_t_fixed": tot["t"], "n_cx": tot["cx"],
        "n_clifford_1q": tot["clifford_1q"],
        "t_report_convention_model_eps_syn": t_model_budget,
        "t_report_convention_default_eps_syn": t_default_budget,
        "t_depth": depth,
        "blocks": {name: blk[name] for name in BLOCK_ORDER if name in blk},
    }
    if variant in mdl["ends"]:
        m = mdl["ends"][variant]
        circ = {
            "logical_qubits": w["qubits_allocated"],
            "O_c_toffoli": blk["O_c"]["toffoli"],
            "O_A_rotations": blk["O_A"]["rot"],
            "projector_toffoli": blk["proj_and"]["toffoli"],
            "projector_rotations": blk["proj_rot"]["rot"],
            "n_toffoli": tot["toffoli"],
            "n_rot": tot["rot"],
            "t_per_rotation": t_model_budget["t_per_rotation"],
            "t_total": t_model_budget["t_total"],
            "t_depth": depth["commuting"]["t_depth"],
        }
        out["comparison"] = [{"item": k, "circuit": circ[k], "model": m[k][0], "ratio": circ[k] / m[k][0],
                              "model_source": m[k][1]} for k in circ]
    return out


def main(write_qasm: bool = False) -> dict:
    mdl = model_items()
    inst = build.Instance(**build.BENCHMARK)
    poly = build.benchmark_polynomial()
    result = {
        "chapter": "Ch. 4, hybrid lattice QCD, 2028 benchmark",
        "instance": {**build.BENCHMARK, "V": inst.V, "s": inst.sub, "theta_0": inst.thetas[0],
                     "theta_1": inst.thetas[1], "d": poly["d"], "shift": poly["shift"],
                     "lp_error": poly["lp_error"], "qsp_residual": poly["qsp_residual"]},
        "conventions": {
            "toffoli": "7 T, common.toffoli_t(n, 'textbook')",
            "rotation": "common.t_per_rotation(common.eps_rot_for(N_rot, eps_syn)), N_rot counted from the circuit",
            "eps_syn_model": mdl["eps_syn"], "eps_syn_default": EPS_SYN,
            "t_depth": "weighted longest path: 1 per Toffoli, t_per_rotation per arbitrary rz, 1 per exact T. "
                       "'commuting' lets gates that share only same-type qubits overlap (the model's convention), "
                       "'as_written' orders every pair of gates that share a qubit",
            "single_machine": "all counts are per shot of one circuit on one machine",
        },
        "model": mdl,
        "variants": {v: count_variant(v, mdl) for v in build.VARIANTS},
    }
    (HERE / "counts.json").write_text(json.dumps(result, indent=1) + "\n")
    if write_qasm:
        from qiskit import qasm3
        for v in build.VARIANTS:
            qc = build.benchmark_circuit(v, measure=True)
            (HERE / f"logdet_4x4x4x4_d84_{v}.qasm").write_text(qasm3.dumps(qc))
    return result


def _print(result: dict) -> None:
    for v, c in result["variants"].items():
        t = c["t_report_convention_model_eps_syn"]
        print(f"\n== {v}: {c['logical_qubits']} qubits ({c['qubits_touched_by_a_gate']} touched by a gate), "
              f"N_Toff {c['n_toffoli']}, N_rot {c['n_rot']}, N_T_fixed {c['n_t_fixed']}, "
              f"T {t['t_total']:.1f} at {t['t_per_rotation']:.3f} T per rotation "
              f"(default eps_syn: {c['t_report_convention_default_eps_syn']['t_total']:.1f})")
        for o, dd in c["t_depth"].items():
            print(f"   T-depth {o:10s} {dd['t_depth']:.1f}  ({dd['toffoli_layers_on_path']} Toffoli layers, "
                  f"{dd['rotations_on_path']} rotations on the path)")
        for name, b in c["blocks"].items():
            print(f"   {name:9s} Toffoli {b['toffoli']:5d}  rot {b['rot']:4d}  T {b['t']:3d}  cx {b['cx']:4d}  "
                  f"1q Clifford {b['clifford_1q']:5d}")
        for row in c.get("comparison", []):
            print(f"   {row['item']:20s} circuit {row['circuit']:10.4g}  model {row['model']:10.4g}  "
                  f"ratio {row['ratio']:.4f}")


if __name__ == "__main__":
    _print(main(write_qasm="--qasm" in sys.argv))
