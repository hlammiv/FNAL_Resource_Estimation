"""Count the Ch. 2 2028 circuit by walking the circuit object, and compare with the chapter model.

    python count.py        # writes counts.json and prints the comparison table

Classification of every instruction (anything else raises):
    Clifford       h, s, sdg, x, z, cx, and rz with an angle that is an even multiple of pi/4
    fixed T        t, tdg, and rz with an angle that is an odd multiple of pi/4
    Toffoli        ccx
    rotation       rz with any other angle (kept in the circuit with its angle)
T totals use estimates/common.py: 7 T per Toffoli (toffoli_t, "textbook") and t_per_rotation(eps_rot_for(N_rot)).
"""
from __future__ import annotations

import json
import math
import sys
from collections import OrderedDict
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from estimates import ch02_nu_nucleus as model                      # noqa: E402
from estimates.common import (EPS_SYN, RUS_ANCILLA, depth_exports,   # noqa: E402
                              eps_rot_for, t_per_rotation, toffoli_t)
from build import LOGICAL_GATES, benchmark_instance, build_circuit   # noqa: E402

ANGLE_TOL = 1e-9
CLIFFORD = {"h", "s", "sdg", "x", "z", "cx"}
KEYS = ("rotations", "toffoli", "fixed_T", "clifford")


def classify(name, params):
    """'clifford', 'fixed_T', 'toffoli', 'rotation' or 'skip' for one instruction."""
    if name in ("measure", "barrier"):
        return "skip"
    if name in CLIFFORD:
        return "clifford"
    if name in ("t", "tdg"):
        return "fixed_T"
    if name == "ccx":
        return "toffoli"
    if name == "rz":
        k = float(params[0]) / (math.pi / 4)
        if abs(k - round(k)) < ANGLE_TOL:
            return "fixed_T" if round(k) % 2 else "clifford"
        return "rotation"
    raise ValueError(f"gate {name!r} is outside the logical gate set {LOGICAL_GATES}")


def walk(circuit, start=0, stop=None):
    """Counts over circuit.data[start:stop]."""
    out = dict(rotations=0, toffoli=0, fixed_T=0, clifford=0)
    key = {"rotation": "rotations", "toffoli": "toffoli", "fixed_T": "fixed_T", "clifford": "clifford"}
    for ci in circuit.data[start:stop]:
        kind = classify(ci.operation.name, ci.operation.params)
        if kind != "skip":
            out[key[kind]] += 1
    return out


def t_total(c, n_rot_for_eps=None):
    """(total T, T per rotation, eps_rot) in the report convention, from estimates/common.py."""
    n_rot = c["rotations"] if n_rot_for_eps is None else n_rot_for_eps
    eps = eps_rot_for(n_rot, EPS_SYN) if n_rot else None
    c_t = t_per_rotation(eps, "rus") if n_rot else 0.0
    return c["rotations"] * c_t + toffoli_t(c["toffoli"], "textbook") + c["fixed_T"], c_t, eps


def block_family(label):
    """Collapse 'prep/3/kinetic' to 'prep/kinetic' and so on."""
    parts = label.split("/")
    if parts[0] == "prep":
        return f"prep/{parts[2]}"
    if parts[0] == "evolve":
        return f"evolve/{parts[-1]}"
    return label


def per_block(built):
    """Counts per block and per block family. Asserts that the blocks partition the circuit."""
    qc = built.circuit
    pos = 0
    blocks = []
    for label, start, stop in built.blocks:
        assert start == pos, f"gap or overlap before block {label}"
        pos = stop
        blocks.append(dict(block=label, **walk(qc, start, stop)))
    assert pos == len(qc.data), "instructions after the last block"
    fam = OrderedDict()
    for b in blocks:
        f = fam.setdefault(block_family(b["block"]), dict(blocks=0, rotations=0, toffoli=0, fixed_T=0, clifford=0))
        f["blocks"] += 1
        for k in KEYS:
            f[k] += b[k]
    total = walk(qc)
    for k in KEYS:
        assert total[k] == sum(b[k] for b in blocks)
    return blocks, fam, total


# --------------------------------------------------------------------------- #
# Depth
# --------------------------------------------------------------------------- #

_CCX_DEF = None


def _ccx_definition():
    """Qiskit's 7-T decomposition of ccx as [(name, local qubit indices)], used only for the weighted depth."""
    global _CCX_DEF
    if _CCX_DEF is None:
        from qiskit.circuit.library import CCXGate
        d = CCXGate().definition
        _CCX_DEF = [(ci.operation.name, [d.find_bit(q).index for q in ci.qubits]) for ci in d.data]
        assert sum(n in ("t", "tdg") for n, _ in _CCX_DEF) == 7
        assert all(n in ("h", "cx", "t", "tdg") for n, _ in _CCX_DEF)
    return _CCX_DEF


def depth(circuit, rotation_weight=1.0, decompose_toffoli=False):
    """Longest path through the circuit counting only non-Clifford gates.

    Each T counts 1 and each arbitrary rotation counts `rotation_weight`. A Toffoli counts 1 layer when
    decompose_toffoli is False, and is replaced by Qiskit's 7-T decomposition when True. Clifford gates are free
    but a two-qubit Clifford synchronizes its qubits, as in any DAG longest-path depth.
    """
    d = [0.0] * circuit.num_qubits

    def apply(name, params, qs):
        kind = classify(name, params)
        if kind == "skip":
            return
        if kind == "toffoli" and decompose_toffoli:
            for n, loc in _ccx_definition():
                apply(n, (), [qs[i] for i in loc])
            return
        w = {"clifford": 0.0, "fixed_T": 1.0, "toffoli": 1.0, "rotation": rotation_weight}[kind]
        if len(qs) == 1 and w == 0.0:
            return
        new = max(d[q] for q in qs) + w
        for q in qs:
            d[q] = new

    for ci in circuit.data:
        apply(ci.operation.name, ci.operation.params, [circuit.find_bit(q).index for q in ci.qubits])
    return max(d)


def depths(circuit, c_t):
    """The two depths reported: non-Clifford gate layers, and T layers at the report's synthesis cost."""
    return dict(
        non_clifford_layers=depth(circuit, 1.0, False),
        t_depth_report=depth(circuit, c_t, True),
        t_depth_toffoli_one_layer=depth(circuit, c_t, False),
        toffoli_t_depth_used=depth_of_one_toffoli(),
    )


def workspace_at_rotations(built):
    """How many clean workspace qubits are idle when each arbitrary rotation is applied.

    A workspace qubit is busy between the Toffoli that computes into it and the Toffoli that uncomputes it.
    Repeat-until-success synthesis holds estimates.common.RUS_ANCILLA clean qubits per rotation in flight, so a
    rotation applied with no idle workspace needs a qubit beyond the register.
    """
    inst, qc = built.instance, built.circuit
    busy = {w: False for w in inst.work}
    none_free = on_work = on_anc = 0
    min_free = len(busy)
    for ci in qc.data:
        qs = [qc.find_bit(q).index for q in ci.qubits]
        name = ci.operation.name
        if name == "ccx" and qs[2] in busy:
            busy[qs[2]] = not busy[qs[2]]
        elif classify(name, ci.operation.params) == "rotation":
            free = sum(not b for b in busy.values())
            min_free = min(min_free, free)
            if free < RUS_ANCILLA:
                none_free += 1
                on_work += qs[0] in busy
                on_anc += qs[0] == inst.anc
    assert not any(busy.values()), "workspace not returned"
    return dict(rus_ancilla=RUS_ANCILLA, min_idle_workspace_at_a_rotation=min_free,
                rotations_with_too_few_idle=none_free, of_which_on_workspace_flag=on_work,
                of_which_on_hadamard_ancilla=on_anc,
                qubits_if_rus_needs_its_own_ancilla=qc.num_qubits + (RUS_ANCILLA if none_free else 0))


def depth_of_one_toffoli():
    from qiskit import QuantumCircuit
    qc = QuantumCircuit(3)
    qc.ccx(0, 1, 2)
    return depth(qc, 1.0, True)


# --------------------------------------------------------------------------- #
# Model line items
# --------------------------------------------------------------------------- #

def model_items(steps):
    """What estimates/ch02_nu_nucleus.py says for the 2028 circuit with `steps` evolution steps."""
    a = model.Assumptions()
    nr, nt, nf = model.counts("2028", steps)
    total, rate = model.circuit(a, "2028", steps)
    return dict(rotations=nr, toffoli=nt, fixed_T=nf, t_total=total, t_per_rotation=rate,
                source="estimates/ch02_nu_nucleus.py counts('2028', steps) and circuit(a, '2028', steps)")


def ratio(x, y):
    return None if not y else x / y


def summarize(inst):
    built = build_circuit(inst)
    blocks, fam, total = per_block(built)
    t, c_t, eps = t_total(total)
    qc = built.circuit
    for f in fam.values():       # T of each family at the whole circuit's per-rotation cost
        f["t_report"] = f["rotations"] * c_t + toffoli_t(f["toffoli"], "textbook") + f["fixed_T"]
    used = sorted({qc.find_bit(q).index for ci in qc.data for q in ci.qubits})
    out = dict(
        qubits=dict(total=qc.num_qubits, used=len(used), position=inst.n_pos, workspace=inst.n_work, ancilla=1),
        counts=total, t_total_report=t, t_per_rotation=c_t, eps_rot=eps,
        depth=depths(qc, c_t), rus_workspace=workspace_at_rotations(built), families=fam, blocks=blocks,
        gates_total=len(qc.data),
    )
    return built, out


def main():
    a = model.Assumptions()
    res = model.model(a, "2028")
    out = OrderedDict()
    out["conventions"] = dict(
        t_per_toffoli=toffoli_t(1, "textbook"), eps_synth=EPS_SYN,
        rotation_cost="estimates.common.t_per_rotation(eps_rot_for(N_rot), 'rus')",
        fixed_T="t, tdg and rz at odd multiples of pi/4", angle_tolerance=ANGLE_TOL,
        depth="non_clifford_layers: each T, Toffoli or arbitrary rz is one layer. t_depth_report: T layers with "
              "each arbitrary rz counted as c_T sequential T and each Toffoli replaced by Qiskit's 7-T circuit. "
              "The chapter model defines no T-depth (depth_status = serial_t_count_bound).")
    out["model"] = dict(
        lq=res.lq[0], hard_ops=res.hard_ops[0],
        breakdown=[dict(name=p.name, count=p.count, t_each=p.t_each, t_total=p.t_total) for p in res.breakdown],
        counts_formula="counts('2028', n) = (646 + 91 n, 1208 + 172 n, 456 + 48 n)",
        per_step=dict(zip(("rotations", "toffoli", "fixed_T"),
                          (y - x for x, y in zip(model.counts("2028", 5), model.counts("2028", 6))))),
        source="estimates/ch02_nu_nucleus.py model(Assumptions(), '2028'), counts, circuit")

    out["by_steps"] = OrderedDict()
    summaries = {}
    for steps in (2, 4, 6):
        built, s = summarize(benchmark_instance(steps))
        summaries[steps] = s
        m = model_items(steps)
        row = dict(circuit=dict(s["counts"], t_total=s["t_total_report"], t_per_rotation=s["t_per_rotation"],
                                qubits=s["qubits"]["used"]),
                   model=m,
                   ratio_circuit_over_model={k: ratio(s["counts"][k], m[k]) for k in ("rotations", "toffoli", "fixed_T")})
        row["ratio_circuit_over_model"]["t_total"] = ratio(s["t_total_report"], m["t_total"])
        row["ratio_circuit_over_model"]["qubits"] = ratio(s["qubits"]["used"], res.lq[0])
        out["by_steps"][str(steps)] = row

    deep = summaries[6]
    out["benchmark"] = deep

    # Per-step increment measured from the circuit (difference of the 6- and 4-step circuits, per step).
    out["per_step_from_circuit"] = {k: (summaries[6]["counts"][k] - summaries[4]["counts"][k]) / 2
                                    for k in ("rotations", "toffoli", "fixed_T")}

    # Line items of the study's own recount file (the source the model's literals were taken from).
    study = json.loads((REPO / "applications/ch02_neutrino_nucleus/neutrino_feasibility/"
                        "larger_circuit_counts_steps6.json").read_text())
    fam = deep["families"]

    def fam_sum(names):
        return {k: sum(fam[n][k] for n in names) for k in ("rotations", "toffoli", "fixed_T")}

    mine_init = fam_sum(["init/uniform"])
    mine_prep = fam_sum(["init/uniform", "prep/potential", "prep/kinetic"])
    rename = dict(rotations="rotations", toffoli="Toffoli", fixed_T="fixed_T")
    out["study_line_items"] = dict(
        source="applications/ch02_neutrino_nucleus/neutrino_feasibility/larger_circuit_counts_steps6.json",
        initial=dict(circuit=mine_init, study={k: study["initial_counts"][v] for k, v in rename.items()}),
        after_preparation=dict(circuit=mine_prep,
                               study={k: study["after_preparation"][v] for k, v in rename.items()}),
        complete=dict(circuit={k: deep["counts"][k] for k in rename},
                      study={k: study["complete_counts"][v] for k, v in rename.items()}))

    # Every one of the eight Hadamard tests (four current pairs, two components) has the same count.
    variants = {}
    for pair in ((0, 0), (0, 1), (1, 0), (1, 1)):
        for comp in ("re", "im"):
            b = build_circuit(benchmark_instance(6, current_pair=pair, component=comp))
            variants[f"{comp}_p{pair[0]}p{pair[1]}"] = {k: walk(b.circuit)[k] for k in ("rotations", "toffoli", "fixed_T")}
    out["all_eight_tests_same_count"] = len({json.dumps(v, sort_keys=True) for v in variants.values()}) == 1

    # Variant: the phase gates on the Hadamard ancilla all commute (the ancilla is only ever a control), so they
    # merge into one rz. Same unitary, checked in test_circuit.py.
    _, merged = summarize(benchmark_instance(6, merge_ancilla_phases=True))
    out["variant_merged_ancilla_phases"] = dict(
        counts=merged["counts"], t_total_report=merged["t_total_report"], t_per_rotation=merged["t_per_rotation"],
        eps_rot=merged["eps_rot"], depth=merged["depth"], qubits=merged["qubits"]["used"],
        rus_workspace=merged["rus_workspace"], families=merged["families"],
        ratio_to_model={k: ratio(merged["counts"][k], model_items(6)[k]) for k in ("rotations", "toffoli", "fixed_T")}
        | {"t_total": ratio(merged["t_total_report"], model_items(6)["t_total"])})

    # Depth and wall-time exports from estimates/common.py, one machine, serial.
    shots = res.shots[0]
    dx = depth_exports(deep["t_total_report"], deep["depth"]["t_depth_report"], shots, a.t_gate_s.lo, a.shot_overhead_s.lo)
    dx1 = depth_exports(deep["t_total_report"], deep["depth"]["t_depth_toffoli_one_layer"], shots, a.t_gate_s.lo,
                        a.shot_overhead_s.lo)
    out["depth_exports"] = dict(
        source="estimates.common.depth_exports(t_total_report, T-depth, campaign shots, t_gate_s, shot_overhead_s); "
               "campaign shots and the two times from estimates/ch02_nu_nucleus.py; every shot at the deepest circuit",
        shots=shots, model_wall_campaign_s=res.wall_time_s[0],
        t_depth_report=dict(t_depth=deep["depth"]["t_depth_report"], f_star=dx["f_star"][0],
                            floor_wall_s=dx["floor_wall_s"][0], wall_serial_s=dx["wall_serial_s"][0],
                            baseline_ok=dx["baseline_ok"][0]),
        t_depth_toffoli_one_layer=dict(t_depth=deep["depth"]["t_depth_toffoli_one_layer"], f_star=dx1["f_star"][0],
                                       floor_wall_s=dx1["floor_wall_s"][0], baseline_ok=dx1["baseline_ok"][0]))

    (HERE / "counts.json").write_text(json.dumps(out, indent=1) + "\n")

    # ---- printed table -------------------------------------------------------------------------------------
    print("2028 four-nucleon circuit, counts by walking the circuit (circuit / model / ratio)")
    for steps, row in out["by_steps"].items():
        c, m, r = row["circuit"], row["model"], row["ratio_circuit_over_model"]
        print(f" steps {steps}: rot {c['rotations']}/{m['rotations']} ({r['rotations']:.3f})  "
              f"Toff {c['toffoli']}/{m['toffoli']} ({r['toffoli']:.3f})  "
              f"T_fixed {c['fixed_T']}/{m['fixed_T']} ({r['fixed_T']:.3f})  "
              f"T {c['t_total']:.1f}/{m['t_total']:.1f} ({r['t_total']:.4f})  qubits {c['qubits']}/{res.lq[0]}")
    print(" block families, 6 steps:")
    for name, f in deep["families"].items():
        print(f"  {name:28s} x{f['blocks']:<3d} rot {f['rotations']:5d}  Toff {f['toffoli']:5d}  T_fixed {f['fixed_T']:4d}"
              f"  Clifford {f['clifford']:6d}  T {f['t_report']:9.1f} ({100 * f['t_report'] / deep['t_total_report']:.1f}%)")
    print(f"  total gates {deep['gates_total']}, Clifford {deep['counts']['clifford']}")
    print(" depth, 6 steps:", deep["depth"])
    mv = out["variant_merged_ancilla_phases"]
    print(f" merged-ancilla-phase variant: rot {mv['counts']['rotations']}  Toff {mv['counts']['toffoli']}  "
          f"T_fixed {mv['counts']['fixed_T']}  T {mv['t_total_report']:.1f}  ratio to model {mv['ratio_to_model']['t_total']:.4f}")
    print(" depth exports:", json.dumps(out["depth_exports"], indent=1))
    print(" workspace at rotations:", deep["rus_workspace"])
    print(" all eight Hadamard tests have the same count:", out["all_eight_tests_same_count"])


if __name__ == "__main__":
    main()
