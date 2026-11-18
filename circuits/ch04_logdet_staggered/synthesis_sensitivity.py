"""How much the Hadamard-test expectation value moves when every arbitrary rz of the benchmark circuit
is off by a small angle. Not a test, and not a synthesis: the rz gates are kept, and each angle is
shifted by +-2 eps_rot with a random sign (an rz off by delta is 2 sin(delta/4), about delta/2, from
the exact one in operator norm, so this is an error of eps_rot per rotation).

    OMP_NUM_THREADS=4 python synthesis_sensitivity.py      # about one minute, writes synthesis_sensitivity.json

Two numbers come out, for z = <Z> on the Hadamard-test ancilla (the normalized shifted trace):
    first_order_rms    spread of (z(+delta) - z(-delta))/2 over sign patterns: the shift of z for ONE fixed
                       set of rotation errors. Randomized synthesis redraws the errors, so this averages out.
    second_order_mean  mean of (z(+delta) + z(-delta))/2 - z(0): the bias that survives the redrawing.
The model budgets eps_syn = N_rot eps_rot^2 for the second and assumes the first averages out.
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
for p in (str(HERE.parents[1]), str(HERE)):   # HERE goes in last, so the modules of this directory are found first
    if p in sys.path:
        sys.path.remove(p)
    sys.path.insert(0, p)

import build                      # noqa: E402
import count                      # noqa: E402
from sim import simulate          # noqa: E402


def main(variant: str = "merged", pairs: int = 24, seed: int = 1) -> dict:
    from qiskit import QuantumCircuit
    mdl = count.model_items()
    qc = build.benchmark_circuit(variant)
    tot = count.walk_counts(qc)["total"]
    eps = count.t_total(tot, mdl["eps_syn"])["eps_rot"]
    n_rot = tot["rot"]

    def z_of(deltas) -> float:
        out = QuantumCircuit(*qc.qregs)
        k = 0
        for g in qc.data:
            if g.operation.name == "rz":
                out.rz(float(g.operation.params[0]) + float(deltas[k]), g.qubits[0])
                k += 1
            else:
                out.append(g.operation, g.qubits)
        psi = simulate(out)
        return 1.0 - 2.0 * float(np.sum(np.abs(psi[2 ** (qc.num_qubits - 1):]) ** 2))

    rng = np.random.default_rng(seed)
    z0 = z_of(np.zeros(n_rot))
    first, second = [], []
    for _ in range(pairs):
        d = 2.0 * eps * rng.choice([-1.0, 1.0], size=n_rot)
        zp, zm = z_of(d), z_of(-d)
        first.append((zp - zm) / 2.0)
        second.append((zp + zm) / 2.0 - z0)
    out = {
        "variant": variant, "n_rot": n_rot, "eps_rot": eps, "angle_shift": 2.0 * eps, "sign_patterns": pairs,
        "seed": seed, "z_exact_angles": z0,
        "first_order_rms": float(np.std(first)),
        "second_order_mean": float(np.mean(second)),
        "second_order_stderr": float(np.std(second) / math.sqrt(pairs)),
        "model_eps_syn": mdl["eps_syn"], "model_poly_tolerance": mdl["poly_tolerance"],
        "first_order_rms_over_poly_tolerance": float(np.std(first) / mdl["poly_tolerance"]),
        "second_order_mean_over_eps_syn": float(abs(np.mean(second)) / mdl["eps_syn"]),
    }
    (HERE / "synthesis_sensitivity.json").write_text(json.dumps(out, indent=1) + "\n")
    return out


if __name__ == "__main__":
    for k, v in main().items():
        print(f"{k:38s} {v}")
