"""Statevector simulation helper shared by the tests and synthesis_sensitivity.py."""
from __future__ import annotations

import numpy as np


def simulate(qc, init=None, threads: int = 4) -> np.ndarray:
    """Final statevector of `qc` from |0...0> or from `init` (qiskit-aer, at most `threads` threads;
    falls back to qiskit.quantum_info.Statevector when qiskit-aer is absent)."""
    n = qc.num_qubits
    if init is None:
        init = np.zeros(2 ** n, dtype=complex)
        init[0] = 1.0
    init = np.asarray(init, dtype=complex)
    try:
        from qiskit_aer import AerSimulator
    except ImportError:
        from qiskit.quantum_info import Statevector
        return np.asarray(Statevector(init).evolve(qc))
    from qiskit import QuantumCircuit
    run = QuantumCircuit(*qc.qregs)
    run.set_statevector(init)
    run.compose(qc, inplace=True)
    run.save_statevector()
    sim = AerSimulator(method="statevector", max_parallel_threads=threads)
    return np.asarray(sim.run(run, shots=1).result().get_statevector())
