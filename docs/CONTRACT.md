# Chapter model contract

Every `chNN_<name>.py` must satisfy `tests/test_common.py::test_model_contract`. Read
`common.py` and `groups.py` first; do not redefine anything they provide.

## Shape

```python
"""Ch. N — <title>. Reproduces the LQ and hard-op numbers of applications/appNN_*.tex.

INPUTS AND SOURCES
  <one line per assumption: symbol, value, where it comes from>
WHAT IS NOT DERIVED HERE
  <every number carried as Stated/Uncited, and why>
WALL TIME AND FACTORIES
  t_gate_s = <value> s per T, shot_overhead_s = <value> s per shot, ONE machine, serial.
  T-depth per shot from factory.json (<run name>), low/high. F* = N_T / D_T; F* >= 10 or the fix.
  First result: <point set, precision> -> wall_first_result_s. Campaign: <target> -> wall_campaign_s.
"""
from dataclasses import dataclass
from estimates.common import (Assumed, Cited, Stated, Uncited, Tagged, Primitive,
                              CircuitStatus, Result, Published, ERAS,
                              t_per_rotation, eps_per_rotation, toffoli_t, depth_exports)
from estimates.groups import GROUPS

@dataclass(frozen=True)
class Assumptions:
    # EVERY field is a Tagged value. Ranges are (lo, hi) tuples.
    n_orb_2028: Tagged = Cited(10, "app01:116", "active 0s-0p-1s0d on A=4")
    t_per_rot:  Tagged = Cited(30, "bocharovRoettelerSvore2015,campbell2017",
                               "RUS, randomized; chapter rounds 24.5 -> 30")
    t_gate_s:        Tagged = Assumed(1e-6, "1 us per T, report convention")      # seconds
    shot_overhead_s: Tagged = Assumed(1e-4, "0.1 ms per shot: init, readout, decode")
    ...
    def __post_init__(self): ...   # validate ranges, forbid nonsense

def model(a: Assumptions, era: str) -> Result:
    """Return the Result for '2028' | '2033' | 'codesign'."""
    ...
    return Result(era=era, lq=(lo, hi), hard_ops=(lo, hi),
                  breakdown=(Primitive(...), ...),
                  intermediates={"rotations_per_step": 5e3, "t_per_step": 1.5e5, ...,
                                 # R9 exports (2028 and 2033 eras), names fixed:
                                 "t_per_shot": ..., "t_depth_per_shot": (lo, hi), "f_star": ...,
                                 "floor_wall_s": ..., "factories_for_1yr": ...,
                                 "wall_first_result_s": ..., "wall_campaign_s": ...},
                  shots=..., wall_time_s=..., epsilon_l=..., notes=(...))

PUBLISHED = {
    "2028": Published(lq=(170, 170), hard_ops=(9e4, 1.5e5), src="app01:2028 box lines 201-205", rel_tol=0.10),
    "2033": Published(...),
    "codesign": Published(...),
}

def INSTANCE_ROWS(a, era, r):
    """Rows for resources.json: [(label, (lq_lo, lq_hi), (t_lo, t_hi), {extra}), ...].
    One per costed instance the chapter's box quotes (a chapter may have several per era)."""
    return [...]
```

## Rules

1. **Provenance on everything.** A bare number in `Assumptions` fails the test. Use
   `Cited(v, bibkey)` when a published source gives it, `Assumed(v, note)` when the
   chapter declares it a working assumption, `Stated(v, "appNN:line")` when the chapter
   asserts it without derivation, `Uncited(v)` when there is no source at all.
   Be honest: `Stated` and `Uncited` are the useful tags.

2. **`breakdown` is primitive-level, and sums to `hard_ops`.** Not "3e3 T/step" but the
   Toffolis, the synthesized rotations (count, eps_rot, t_each), the grouped rotations,
   the group actions, each a `Primitive` with a `CircuitStatus`. Multiply through the
   Trotter steps / queries / rounds so the sum is the per-shot total. The test checks the
   sum lands within the quoted range.

3. **`intermediates` holds every named number the prose states on the way** to the
   headline — the "5e3 rotations", the "700 × 3e4", the "0.11–0.48× the box". The
   verifier diffs each one. A hidden factor is a finding.

4. **Rotation synthesis goes through `common.t_per_rotation` / `eps_per_rotation`.**
   If the chapter uses "~30 T", set `t_per_rot = Cited(30, ...)` and record in the note
   that the fit gives 24.5 at 1e-4. Never hard-code a second synthesis formula.

5. **Toffoli→T names its convention** (`toffoli_t(n, "textbook"|"jones")`). Ch. 9 uses
   `"jones"` (4 T) via arXiv:1709.06648; most others `"textbook"` (7 T). If the chapter
   doesn't say, that is a `Stated` input and a note.

6. **Gauge-group widths and primitive costs come from `groups.py`.** Use
   `GROUPS[k].link_qubits` (the chapter's value) so the recorded conflict surfaces; do not
   type 7 or 8 inline.

7. **`PUBLISHED` is the box as printed.** You may not change it to make a test pass. If
   the model cannot reproduce it, the model returns what the physics gives, the test
   fails, and you write the discrepancy into `NEEDS_AUTHOR.md` under the chapter heading
   with both numbers and the reason.

8. **`rel_tol`** defaults to 0.10; tighten it when the chapter quotes an exact product,
   loosen (up to 0.5) only for order-of-magnitude rows, and say why in `Published.src`.

9. Standard library only. No numpy, no imports from outside `estimates/`.

10. **Line references** in `src` strings are to the current `.tex` (run `grep -n`); they
    will drift, so also quote a few words of the sentence in `note`.

11. **Unified wall-time names (ruling R9, H. Lamm, 2026-10-02).** Every `Assumptions` carries
    `t_gate_s` (seconds per T; replaces `logical_cycle_s`, `gate_time_s`, `gate_time_us`) and
    `shot_overhead_s` (seconds per shot; 1e-4 by the report convention). The old names are not
    allowed. There is no central override: each chapter keeps its own values and tags.

12. **Depth and factory exports (R9).** For each of the `2028` and `2033` eras in `PUBLISHED`,
    `Result.intermediates` holds these keys, with exactly these names:

    | key | meaning |
    |---|---|
    | `t_per_shot` | T (hard ops) per shot, (lo, hi) |
    | `t_depth_per_shot` | T-depth per shot from factory.json, (lo, hi) |
    | `f_star` | t_per_shot / t_depth_per_shot, the factory count above which factories idle |
    | `floor_wall_s` | shots x t_depth x 10 us (`REACTION_TIME_S`): shortest possible run |
    | `factories_for_1yr` | 10 (`FACTORY_BASELINE`) x serial years at t_gate_s per T |
    | `wall_first_result_s` | single-machine wall of the minimal first result (ruling R1); `None` only for a 2028 benchmark with one tier |
    | `wall_campaign_s` | single-machine wall of the full campaign at the chapter's physics target |

    Compute the first five with `common.depth_exports(t_per_shot, t_depth_per_shot, shots,
    a.t_gate_s, a.shot_overhead_s)` and copy its keys; do not re-derive them. Where `f_star < 10`
    the box quotes the corrected wall (or the workspace fix, with the LQ updated), never the
    10-factory baseline. All walls are serial on ONE machine. `tests/test_common.py::
    test_model_exports_depth` checks the keys, the names in `Assumptions`, the docstring header,
    and that `f_star` equals the ratio.

13. **Module docstring header.** The docstring starts `Ch. N — <title>.`, names its
    `applications/appNN_*.tex`, and has the three section headings shown above, in that order:
    `INPUTS AND SOURCES`, `WHAT IS NOT DERIVED HERE`, `WALL TIME AND FACTORIES`. Free text may follow.

## What "reproduce" means

The model regenerates the chapter's box numbers *from the chapter's own stated inputs*.
It does not improve the physics. Where the chapter's inputs are wrong, or its arithmetic
does not close, the model exposes that — it does not paper over it. Improvements are the
authors' job and go through `NEEDS_AUTHOR.md`.

## Chapter 2 revision, 2026-10-04

The accepted finite-volume benchmarks establish complete-circuit counts but no device schedule. Their depth/factory/wall-time exports are explicitly `None` with `depth_status="not_established"`; gate volumes are reported instead. Do not substitute legacy timing estimates. The carbon codesign comparator is evolution scaling, not a complete shot.
