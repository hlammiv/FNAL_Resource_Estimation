"""Pauli strings per Trotter step of a generic rank-0 one- plus two-body nuclear Hamiltonian, per space.

Runs estimates.ch03_mu2e_0nubb.hamiltonian_strings_count (the structured count; checked against explicit
JW algebra in sd, f7p3, f7p3p1, pf, Mo/Ru and e_max = 1, 2) and writes strings_<space>.json here. The
oscillator spaces e_max = 8 and 10 take minutes; the model reads their counts from HAMILTONIAN_STRINGS_STORED.

Usage (repository root): python3 applications/ch03_mu2e_0nubb/pauli/count_strings.py emax8 emax10
"""
import json
import pathlib
import sys
import time

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[3]))
from estimates import ch03_mu2e_0nubb as ch3  # noqa: E402

OUT = pathlib.Path(__file__).resolve().parent

for key in sys.argv[1:]:
    t0 = time.time()
    r = ch3.hamiltonian_strings_count(*ch3.step_space(key))
    r.update(space=key, seconds=round(time.time() - t0, 1))
    (OUT / f"strings_{key}.json").write_text(json.dumps(r, indent=1) + "\n")
    print(json.dumps(r), flush=True)
