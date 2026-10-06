# FNAL Resource Estimation

This repository holds the resource models behind the Fermilab whitepaper *A Program for Scientific Quantum Utility at Fermilab* (2026). Each application chapter of the paper (Chapters 2 to 10) has one Python model here. The model takes the chapter's physics inputs, builds the circuit cost from named primitives, and returns the logical-qubit count, T-gate count, shot count and single-machine wall time that the chapter prints in its benchmark boxes and that the paper collects in its summary table. Every input carries its provenance, and a test suite pins every published number to the model that produces it. The paper itself is not included; the only thing taken from it is this code and the `resources.json` file the code generates.

## Conventions every model shares

These are set once in `estimates/common.py` or stated once in the paper's overview and carried by each chapter's `Assumptions`.

- **Logical qubits** count the algorithmic register only. Magic-state factories and distillation are excluded, as in the DOE reference budgets below.
- **Toffoli gates** cost 7 T (`T_PER_TOFFOLI["textbook"]`).
- **Rotations** are synthesized by repeat-until-success at 1.15 log2(1/ε_rot) + 9.2 T each (`RUS_SLOPE`, `RUS_OFFSET`; Bocharov, Roetteler and Svore). The per-rotation tolerance is set per circuit from that circuit's own rotation count, ε_rot = sqrt(10⁻² / N_rot) (`eps_rot_for`, `EPS_SYN = 1e-2`), so that the synthesis errors of one circuit, added incoherently, total 1%. Chapter 2 is the one exception: it prices rotations at the leading coherent Ross-Selinger term, 3 log2(N_R / 0.005) T, as its chapter text does (`coherent_synthesis_allowance` in `ch02_nu_nucleus.py`).
- **Fault budget**: 0.1 expected logical faults per circuit. This fixes the logical error rate each chapter quotes, ε_l ≈ 0.1 / (hard operations per circuit).
- **Wall time** is for one machine, run serially: 1 µs per T gate plus 0.1 ms per shot. The 1 µs corresponds to about ten magic-state factories each delivering one T state per 10 µs logical cycle (`FACTORY_BASELINE = 10`). The models also report the T-depth, the factory count F* = N_T / D_T above which extra factories sit idle, and the shortest possible run at one 10 µs reaction time per sequential non-Clifford layer (`REACTION_TIME_S`, `depth_exports`). Nothing assumes several machines.
- **Two tiers** per instance where the chapter gives them: a first result at 30% precision, and a campaign at the chapter's target precision.
- **DOE reference budgets**: 150-250 logical qubits and 10⁵ T (or Toffoli) gates for 2028, and about 1000 logical qubits and 10⁹ hard gates for 2033. They are reference points for comparison, not limits. Several instances exceed them, and the chapters say by how much.

## Provenance

Every model input is a `Tagged` value (`estimates/common.py`) with one of these tags:

- **Cited**: a published source gives the value; the source is recorded (an arXiv identifier or bibliography key, often with a file and line).
- **Assumed**: a working assumption the chapter declares, with a note on why.
- **Stated** (printed as `stated-not-derived`): a value the chapter uses without deriving it, recorded with the chapter line where it appears.

A fourth tag, **Uncited**, exists for a value with no source at all. Separately, each primitive in a T-count carries a circuit status (compiled, scaling argument, conjecture, unsourced), and `docs/CIRCUIT_STATUS.md` is generated from those statuses. `python -m estimates --inputs` prints every input with its tag and source.

## Layout

```
estimates/              the package
  common.py             rotation synthesis, Toffoli and wall-time conventions, provenance tags
  groups.py             gauge-group primitives (qubits per link, plaquette and group-action T-counts)
  chNN_*.py             one model per application chapter (Chapters 2-10)
  ch02_legacy.py        the longer-term 12C scaling row of Chapter 2
  apply_log/            the polynomial-degree fits behind Chapter 4's QET degrees (need scipy;
                        not imported by the package)
  tests/                pytest suite pinning every published number
  __main__.py           command-line entry point
applications/chNN_*/    one README per chapter (what is estimated, every assumption with its
                        provenance, how the numbers are built, open items) and the
                        dump_assumptions.py script that prints its assumptions table
docs/                   CONTRACT.md (what every chapter model must expose), CIRCUIT_STATUS.md
                        (generated), PAPER_COSTS.md, OPEN_ITEMS.md (what the models cannot decide)
resources.json          the generated numbers, in the form the paper's figures and table read
sync_paper.sh           copies estimates/ into the paper's working tree and runs its suite there
```

## Running

Python 3.11 (the version CI tests) with `numpy`, `scipy` and `sympy` (`pip install -r requirements.txt`).

```
python -m pytest estimates/tests -q        # the full suite
python -m estimates                        # every instance: qubits and T-count against the printed values
python -m estimates --inputs               # the same, with every assumption and its provenance
python -m estimates --chapter ch03         # one chapter: assumptions, resources.json rows, exports
python -m estimates --chapter ch03 --intermediates   # every intermediate the model records
python -m estimates --resources            # regenerate resources.json
python -m estimates --status               # regenerate docs/CIRCUIT_STATUS.md
```

Some tests check that a chapter's text quotes the model's numbers. They need the paper source alongside and skip without it; section 2 ("How to run") of each chapter README says which tests these are and what they check before they skip.

## Applications

| Ch. | Application | What the model prices |
|---|---|---|
| 2 | [Neutrino-nucleus response for DUNE](applications/ch02_neutrino_nucleus/README.md) | Real-time charge correlators of finite-volume nuclear models by Hadamard tests, plus a longer-term 12C scaling row |
| 3 | [Mu2e and neutrinoless double-beta decay](applications/ch03_mu2e_0nubb/README.md) | Nuclear matrix elements of muon-to-electron conversion and 0νββ operators in shell-model spaces, ground state by Trotterized projection |
| 4 | [Hybrid lattice QCD](applications/ch04_hybrid_lqcd/README.md) | Fermion-matrix traces (log det, Tr M⁻¹) by block encoding and QSVT on classically generated gauge configurations |
| 5 | [Transport in the quark-gluon plasma](applications/ch05_qgp_transport/README.md) | The shear-channel stress-tensor autocorrelator and its decay rate on discrete-subgroup gauge theories |
| 6 | [QCD for colliders](applications/ch06_collider/README.md) | Real-time hadronization of a quark-antiquark pair, hadron-species ratios, fragmentation, and dipole survival in a medium |
| 7 | [Quantum fields in curved spacetime](applications/ch07_curved_space/README.md) | Particle production during reheating on an expanding FLRW lattice |
| 8 | [Baryogenesis](applications/ch08_baryogenesis/README.md) | False-vacuum decay rate and C-violating fermion reflection off a bubble wall |
| 9 | [Chiral lattice fermions](applications/ch09_chiral_gauge/README.md) | Domain-wall and overlap fermions coupled to gauge fields for real-time chiral gauge dynamics |
| 10 | [QCD at finite baryon density](applications/ch10_finite_density/README.md) | The equation of state and baryon-number susceptibilities at nonzero chemical potential |

Chapters 11 (qudit platforms) and 12 (dark-matter sensing) of the paper are parametric models stated in the paper itself and are not part of this repository.

## How to cite

For the estimates, cite the paper:

> N. Holtkamp *et al.*, *A Program for Scientific Quantum Utility at Fermilab*, Fermilab whitepaper (2026).

For the code itself:

> H. Lamm, *FNAL_Resource_Estimation*, https://github.com/hlammiv/FNAL_Resource_Estimation (2026).
