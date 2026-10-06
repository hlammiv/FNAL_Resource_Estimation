# Chapter 4: hybrid lattice QCD

Model: [`estimates/ch04_hybrid_lqcd.py`](../../estimates/ch04_hybrid_lqcd.py). Tests: [`estimates/tests/test_ch04_hybrid_lqcd.py`](../../estimates/tests/test_ch04_hybrid_lqcd.py).

## 1. What is estimated

Classical Monte Carlo produces the gauge configurations; the quantum computer evaluates one trace of the fermion matrix per configuration. It block-encodes the sparse staggered matrix, applies a QET/QSVT polynomial ($\log|x|$ or $1/x$), and reads out a trace. The register grows as $\log_2 V$, so the per-shot T-count, not the register width, is what binds. The model prices three instances.

- **2028 benchmark, $\log\det W$ on free staggered fermions at $V = 4^4$** (`_model_2028`). $W = M^\dagger M$, $m_0 = 0.4$, $K = 1$, periodic. The circuit is the free-fermion block encoding of arXiv:2407.13080 Sec. II.2 (9 nonzeros per column, a 4-qubit column register, one flag qubit), a QET of the shifted normalized logarithm at degree $d = 84$, and a Hadamard-test trace. The validation target is the exact answer $\log\det W = 150.55$ (`logdet_W`, computed from the plane-wave spectrum). The model also prices the gauged $V = 4^4$ circuit as a comparison, and it does not fit the budget.
- **2033 target, $\mathrm{Tr}\,M^{-1}$ on SU(3) staggered $V = 8^3\times 16$** (`_model_2033`). $N_f = 2+1$, $am = 0.04$, $\kappa \approx 10^2$. The block encoding is of $M$ itself at $D V\log_2 V$ T per query. QSVT applies $1/x$ at $d_\mathrm{inv} = 500$, and the readout is the block-encoding success probability $P = p_0^2\,am\,c$, $c = (1/V_F)\,\mathrm{Re\,Tr}\,M^{-1}$. A second row on the same register, the Banks-Casher mode number $\nu(\Lambda)/V_F$ at $a = 0.15$ fm, is the physics demonstration.
- **Co-design row (post-2033), $\mathrm{Tr}\,M^{-1}$ at $V = 24^3\times 48$, $\kappa \approx 10^3$** (`_model_codesign`). The same circuit at $d_\mathrm{inv} = 5\times 10^3$. It is the stage at $\epsilon_l \sim 10^{-12}$ whose target is the disconnected HVP. Its per-shot T is quoted before and after the required $2.6\times$ reduction, and the production lattice $48^3\times 96$ is costed alongside.

## 2. How to run

From the repository root:

```
python -m estimates --chapter ch04                    # assumptions, resources.json rows, per-era exports
python -m estimates --chapter ch04 --intermediates    # also every intermediate the model records
OMP_NUM_THREADS=2 python -m pytest estimates/tests/test_ch04_hybrid_lqcd.py -q
python applications/ch04_hybrid_lqcd/dump_assumptions.py   # the table of Section 3
```

`--chapter ch04` recomputes the model and writes nothing. It prints the 95 assumptions with their provenance, then the four rows this chapter contributes to `resources.json`. For each era it prints the headline LQ and T beside the values the paper prints (with an `ok` flag at the 2% tolerance), followed by shots, wall time, the required $\epsilon_l$, the T-depth and factory exports, the primitive breakdown and the model's notes. The intermediates (121, 96 and 77 further keys for 2028, 2033 and the co-design row) appear with `--intermediates`.

The tests need `numpy` and `scipy` for the dense cross-check of the spectrum and for the two linear programs that reproduce $d = 84$ and the $d_\mathrm{inv}$ law (the two LP scripts in `estimates/apply_log/`). In this repository the run gives 29 passed and 19 skipped. The skipped tests first check the model's numbers and then check that the chapter prints them. Here the model checks run and pass, and the test is reported as skipped when it reaches the first chapter-text check. Section 5 says which numbers are pinned in each layout.

## 3. Assumptions

`dump_assumptions.py` reads `Assumptions()` from the live model and prints every field. The provenance tags are those of `estimates/common.py`. `cited` means a published source gives the value. `assumed` is a working assumption of this chapter. `stated-not-derived` is a number the chapter asserts without derivation. `uncited` is asserted with no source. The source and note columns are the code's own strings, verbatim. Short labels in them (R3, R7, SYN-BUDGET, E27 and the like) are the identifiers of author decisions recorded in the module docstring (H2, in the `hvp_site_variance` note, is a referee question; the comments of `_model_codesign` treat it). `appNN:line` is the code's pointer into the chapter. The line numbers are not kept in step with the .tex, so the quoted phrase in the note is the reliable locator. Finally, `txt:NNN` points into the text of arXiv:2407.13080. Rows are grouped by the section comments of the dataclass. "Used by" lists the era functions that read the field, so the conventions stored in one section but used everywhere (`p_f`, `t_gate_s`, `shot_overhead_s`) are visible.

Output of `python applications/ch04_hybrid_lqcd/dump_assumptions.py`:

95 fields: stated-not-derived 45, uncited 2, cited 18, assumed 30

#### Shared: lattice and formulation

| field | value | provenance | used by | source | note |
|---|---|---|---|---|---|
| `D` | 4 | stated-not-derived | 2028, 2033, codesign | app12:48 | D-dimensional Euclidean lattice; 4D (app12:200) |
| `n_s` | 1 | stated-not-derived | 2033 | app12:48 | n_s=1 for staggered fermions |
| `log_base` | 2 | stated-not-derived | 2028, 2033, codesign | app12:55,89 | 'we evaluate the gauged cost as D^2 V log_2 V' (R7); 1.7e6 at V=8192 = 16*8192*13 |
| `scaling_prefactor` | 1.0 | uncited | 2028, 2033, codesign | - | GAUGED circuit only: app12:55 'a working prefactor of one T per unit' of D^2 V log2 V; no compiled circuit and no citation for the prefactor |
| `toffoli_convention` | textbook | cited | 2028 | main-overview:47 (ruling R5, H. Lamm 2026-09-28) | 7 T per Toffoli; app12:125 '7 T per Toffoli' |

#### 2028 benchmark: log det W, free staggered V = 4^4 (this section also holds the shot, wall-time and campaign conventions every era uses)

| field | value | provenance | used by | source | note |
|---|---|---|---|---|---|
| `L_2028` | 4 | assumed | 2028 | app12:26,110 | linear extent, periodic; V = L^D = 4^4 = 256 (app12:110). Ruling R7: the old V = 2^4 instance is degenerate (M = m0 x identity) |
| `m0_2028` | 0.4 | stated-not-derived | 2028 | app12:110 | bare mass m_0 = 0.4 |
| `K_2028` | 1.0 | cited | 2028 | arxiv_2407_13080 | txt:266-267 'The K coupling is typically set to one'; app12:110 'K=1' |
| `nnz_free_2028` | 9 | cited | 2028 | arxiv_2407_13080 | Sec. II.2, txt:323-324 'only nine nonzero entries in each row/column' = 2D+1; app12:111 |
| `column_qubits_2028` | 4 | cited | 2028 | arxiv_2407_13080 | Sec. II.1, txt:103-105 'l = 0-7 ... 8 <= l < 16', Fig. 1 draws l_0..l_3; reused for free fermions (txt:317-319). The paper's 'six qubits' (txt:349) is the U(1)-gauged matrix, 33 nonzeros |
| `be_flag_qubits_2028` | 1 | cited | 2028 | arxiv_2407_13080 | Figs. 2-3: the \|0> qubit the O_A y-rotations act on |
| `qet_signal_qubits_2028` | 1 | assumed | 2028 | - | QET signal qubit carrying the projector-controlled phases; the standard construction, not drawn in the paper |
| `mcx_workspace_2028` | 3 | assumed | 2028 | - | clean workspace for the Toffoli ladders: the 5-control projector AND needs 5-2 = 3; the 4-control shifts of O_c need 2 of the same 3 |
| `hadamard_qubits_2028` | 1 | assumed | 2028 | - | Hadamard-test ancilla. The Hadamard test is this chapter's choice; the paper's trace step is QME (Sec. III) |
| `n_qpe_2028` | 8 | assumed | 2028 | - | m = 8 QPE register on the QME upgrade path; the paper's example runs m = 6, 7, 8, 9 (txt:554) |
| `qme_sign_qubits_2028` | 1 | cited | 2028 | arxiv_2407_13080 | txt:592-594 'QME ... requires one additional qubit, which we shall label a' |
| `poly_shift_2028` | 0.5 | cited | 2028 | arxiv_2407_13080 | Sec. IV, txt:697-703 'it helps to shift the entire polynomial away from y = -1 ... A shift of 1/2 works well. This known shift can always be subtracted at the end of the calculation'; app12:100,113 |
| `target_rel_2028` | 0.01 | assumed | 2028 | - | relative error on log det W: guaranteed for the polynomial by the uniform tolerance, and resolved by the sampling as an rms error (app12:116-117; ruling 'SHIFT + 1%'; rms since R16) |
| `d_log_2028` | 84 | assumed | 2028 | - | LP-minimal degree, this work: smallest even d with uniform error <= 1.784e-3 on the shifted target log\|x\|/log(s/lambda_min) + 1/2 over [lambda_min/s, 1] and \|p\| <= 1 on [-1, 1], at lambda_min/s = 0.0370. Not in the paper |
| `d_log_lam_norm` | 0.037037037037037035 | assumed | 2028 | - | the lambda_min/s the LP degrees here were computed for; the model flags d as stale if the instance moves off it |
| `d_log_err` | 0.0016993 | assumed | 2028 | - | LP optimum at d = 84 with the shift: the uniform error reached |
| `d_log_err_below` | 0.001866 | assumed | 2028 | - | LP optimum at d - 2 = 82 with the shift: above the tolerance, so 84 is minimal |
| `d_shifted_tol_1e2` | 48 | assumed | 2028 | - | LP-minimal degree with the shift at tolerance 1e-2 (err(48) = 9.939e-3, err(46) = 1.104e-2); app12:100 |
| `d_unshifted_tol_1e2` | 112 | assumed | 2028 | - | LP-minimal degree WITHOUT the shift at tolerance 1e-2 (err(112) = 9.907e-3, err(110) = 1.026e-2): the round-B box before this ruling; app12:100 |
| `d_unshifted_2028` | 314 | assumed | 2028 | - | LP-minimal degree WITHOUT the shift at the tolerance 1.784e-3 (err(314) = 1.7704e-3, err(312) = 1.7945e-3): what the shift saves; app12:100,136 |
| `d_lp_at_lam016` | 36 | assumed | 2028 | - | record: LP-minimal degree at lambda_min = 0.16, no shift, tolerance 1e-2 (err 8.99e-3; err(34) = 1.011e-2): what the old 16+12 = 28 box was sized for |
| `d_lp_at_tol_1e3` | 510 | assumed | 2028 | - | record: LP-minimal degree at lambda_min/s = 0.0370, no shift, tolerance 1e-3: err(510) = 9.985e-4, err(508) = 1.0061e-3 |
| `d_fig12_paper` | (64, 70) | cited | 2028 | arxiv_2407_13080 | Fig. 12 caption, txt:644-647: total degree d = 64, 66, 68, 70 with d_r = 60, d = d_f + d_r (txt:686); the paper states NO lambda_min for it |
| `oc_shift_gates_2028` | (4, 8) | assumed | 2028 | - | controlled shifts in O_c: 8 as Fig. 1 draws them (add2 for +mu, sub2 for -mu); 4 when the +2mu and -2mu shifts, which coincide on L = 4, are merged |
| `oc_controls_2028` | (3, 4) | assumed | 2028 | - | controls per shift: the 4 column qubits as drawn; 3 when the sign bit is dropped by the merge |
| `oa_rotations_2028` | (2, 4) | assumed | 2028 | - | O_A at GENERIC angles: Fig. 2 draws two singly-controlled R_y. 4 synthesized rotations if each controlled R_y is 2 rotations + 2 CNOT; 2 if the pair is compiled as one multiplexed rotation. The model drops rotations whose angle is Clifford+T exact (R-TOL): at s = 2(m0^2+2K^2) theta_0 = 0, so both ends are 2 (oa_synthesized_rotations) |
| `projector_toffolis_2028` | 8 | assumed | 2028 | - | one controlled projector phase: the AND of the 5 block-encoding qubits (4 column + 1 flag) is computed onto the signal qubit by a 4-Toffoli ladder and uncomputed by 4 more |
| `and_tree_depth_2028` | 3 | cited | 2028 | factory.json (T-depth audit, 2026-10-02) ch04 2028 run | the 5-control AND onto the signal qubit is a 3-layer tree on the 3 workspace qubits; measurement uncompute is Clifford (T-depth 0) |
| `projector_rotations_2028` | 2 | assumed | 2028 | - | the phase rotation controlled on the Hadamard-test ancilla: 2 synthesized rotations + 2 CNOT. For even d the block-encoding queries need no control |
| `synth_to_poly_2028` | 0.1 | stated-not-derived | 2028 | ruling SYN-BUDGET (H. Lamm, 2026-09-29): EPS_SYN = 1e-2 stays the report-wide default; a circuit may set a tighter budget when its accuracy target needs one. Ch. 4 2028: eps_syn = 0.1 x the polynomial tolerance 1.784e-3 (app12:125) | synthesis budget of the 2028 circuit as a fraction of the polynomial tolerance: eps_syn = 0.1 x 1.784e-3 = 1.784e-4 per shot, tighter than the report-wide EPS_SYN = 1e-2 because the 1% guarantee on log det W rests on the polynomial tolerance; eps_rot = sqrt(eps_syn / N_rot) from the circuit's own rotation count (common.eps_rot_for). Replaces eps_syn = 1e-2 |
| `shot_variance` | 1.0 | stated-not-derived | 2028, codesign | ruling R16 (H. Lamm, 2026-10-01); app12:92 | per-shot variance of the +-1 Hadamard-test outcome, 1 - mu^2 <= 1; the bound 1 is used. Shots = Var/eps^2 at rms additive eps (the convention of Chs. 3 and 8). Replaces confidence delta = 1e-2 and ln(1/delta)/eps^2 |
| `eps_l_2028` | 1e-08 | cited | 2028 | DOE_RFI_2026 | app12:26: the RFI's eps_l=1e-8. The box prints the requirement 0.1/T instead (R3) |
| `t_cap_2028` | 100000.0 | cited | 2028 | DOE_RFI_2026 | app12:26,126: <=1e5 T per circuit |
| `t_gate_s` | 1e-06 | cited | 2028, 2033, codesign | app12:128 | 1 us per T-gate; convention in Ch. overview |
| `shot_overhead_s` | 0.0001 | assumed | 2028, 2033, codesign | - | per-shot overhead t0 ~ 0.1 ms (register init, final readout, decode), report rule, main-overview:48; as Ch. 8 shot_overhead_s. Wall = shots x (T x 1 us + t0), serial on ONE machine (ruling E27, 2026-10-02) |
| `campaign_horizon_yr` | 5 | stated-not-derived | 2033, codesign | app12:Requirements | 'Campaign horizon 5 years'. Serial time on one machine is compared with it; no machine count enters (ruling E27) |

#### 2028 comparison: the gauged V = 4^4 circuit (its levers are reused by the co-design row)

| field | value | provenance | used by | source | note |
|---|---|---|---|---|---|
| `lcu_saving` | (3, 10) | stated-not-derived | 2028, codesign | app12:137 | LCU-based BE compression '~3-10x saving' |
| `sublattice_saving` | 2 | cited | 2028, codesign | arxiv_2407_13080 | Sec. V, txt:755-760 'only half the matrix need be block-encoded, since the even and odd sub-lattices are completely independent' |

#### 2033 target: Tr M^-1 on 8^3x16, and the Banks-Casher demonstration row

| field | value | provenance | used by | source | note |
|---|---|---|---|---|---|
| `V_2033` | 8192 | stated-not-derived | 2033 | app12:144 | V=8^3x16 |
| `VF_box_2033` | 25000.0 | stated-not-derived | 2033 | app12:144 | 'V_F = N_c V ~= 2.5e4 per flavor' as printed (was 1.6e5; ch04-VF-8x8x8x16) |
| `Nc_2033` | 3 | stated-not-derived | 2033, codesign | app12:145 | SU(3) staggered |
| `Nf_2033` | 3 | stated-not-derived | 2033 | app12:145 | N_f=2+1 |
| `kappa_2033` | 100.0 | stated-not-derived | 2033, codesign | app12:146 | kappa ~ 1e2 (am ~ 0.04, m_pi ~ 600 MeV, near the strange mass; was m_pi ~ 300-400 MeV) |
| `d_inv_eps_rel` | 0.006737946999085467 | assumed | 2033, codesign | - | uniform error of the 1/x polynomial relative to the peak of 1/x on [1/kappa, 1]: e^-5 = 0.674% (app12:147,169; ruling R11 ch04-dinv-500-1e4 B') |
| `d_inv_2033` | 500 | assumed | 2033, codesign | - | LP-minimal odd degree, this work: kappa ln(1/eps_rel) = 1e2 x 5 (apply_log/ch04_r11_lp_odd.py: err(501) = 3.334e-3 = 0.5 e^-5.01); app12:147 'd_inv = 500, the smallest degree ... (linear program, this work)' |
| `neighbor_box_2033` | 4 | stated-not-derived | 2033 | app12 1000-LQ box | '13 site + 4 column index + 4 color' (R3, 2026-10-02: the box block-encodes M, 2D+1 = 9 entries; was 6, the M^dag M count) |
| `n_qpe_2033` | 12 | stated-not-derived | 2033 | app12:152 | 12 QPE |
| `n_anc_2033` | 50 | stated-not-derived | 2033 | app12:152 | ~50 ancilla (BE workspace + phase loader + QME) |
| `eps_2033` | 0.01 | stated-not-derived | 2033, codesign | app12 1e4-LQ paragraph | eps=1e-2 rms per configuration (R16): the Hadamard-test target the 1000-LQ box used before R3; kept as the 1e4-LQ rung's inherited shot plan (1e4 per configuration) |
| `n_cfg_2033` | 100 | stated-not-derived | - | retired 2033 plan (R3) | ~1e2 configurations; informational since R3 |
| `be_D_power_M` | 1 | stated-not-derived | 2033, codesign | ruling R3 (H. Lamm, 2026-10-02); shot audit ch04 needs_author 4 | the Tr M^-1 rungs block-encode M, 2D+1 = 9 entries per column: D V log2 V per query. D^2 V log2 V (33 entries of M^dag M) stays for the gauged 4^4 log det |
| `L_s_2033` | 8 | stated-not-derived | 2033 | app12 1000-LQ box | V = 8^3 x 16, spatial extent |
| `L_t_2033` | 16 | stated-not-derived | 2033 | app12 1000-LQ box | V = 8^3 x 16, temporal extent (antiperiodic for the free c) |
| `am_2033` | 0.04 | stated-not-derived | 2033 | ruling R3; shot audit ch04 change (1) | bare mass am ~ 0.04: kappa = (4 + am)/am ~ 1e2 at subnormalization 4 + am (am = 4/99) |
| `condensate_interacting_2033` | 0.08 | assumed | 2033 | - | shot audit ch04: interacting estimate of c = (1/V_F) Re Tr M^-1 at am = 0.04, 'up to 0.08'; the free value is computed here |
| `peak_norm_2033` | 0.5 | assumed | 2033, codesign | - | QET polynomial normalized to 1/2 at sigma_min (the 0.5 e^-5 tolerance of the uncertainty row); P = f^2 am c (shot audit readout lever 1) |
| `rel_first_2033` | 0.3 | stated-not-derived | 2033 | ruling R1 (H. Lamm, 2026-10-02) | first result at 30% statistical error |
| `n_cfg_first_2033` | 1 | stated-not-derived | 2033 | ruling R1; shot audit ch04 minimal_first_result | one configuration |
| `rel_campaign_2033` | 0.2 | stated-not-derived | 2033 | ruling R3 | campaign: 20% relative against the classical value |
| `n_cfg_campaign_2033` | 3 | stated-not-derived | 2033 | ruling R3 '1-3 configurations' | three priced, the top of 1-3 |
| `parallelism_2033` | (1, 4) | cited | 2033 | factory.json (T-depth audit, 2026-10-02) ch04 1000-LQ runs | F* = N_T / D_T: 1 for a multiplexed rotation on the flag qubit (every T on the critical path), 4 for the unary-iteration (QROM) link load. Ratio-preserving: D_T = N_T / F*, carried to the D V log2 V count |
| `subtree_lq` | 34 | cited | 2033 | factory.json (T-depth audit, 2026-10-02) ch04 1000-LQ campaign | ~2 log2 L LQ per parallel unary-iteration subtree, F* ~ 4k for k subtrees; unpriced co-design lever |
| `bc_a_fm` | 0.15 | stated-not-derived | 2033 | simplobs.json ch04 C2; ruling R7 | coarse a = 0.15 fm variant (MILC HISQ class) |
| `bc_hbarc_gev_fm` | 0.1973 | cited | 2033 | PDG | hbar c = 0.1973 GeV fm |
| `bc_sigma_third_gev` | 0.27 | assumed | 2033 | - | simplobs.json ch04 C2: Sigma^(1/3) = 270 MeV; P uncertain by ~2x in the small box |
| `bc_tastes` | 4 | stated-not-derived | 2033 | simplobs.json ch04 C2 | 4 staggered tastes |
| `bc_lambda_lo` | 0.1 | stated-not-derived | 2033 | simplobs.json ch04 C2 | Lambda = 0.10 (lattice units), d = 640 |
| `bc_lambda_hi` | 0.15 | stated-not-derived | 2033 | simplobs.json ch04 C2 | Lambda = 0.15 (lattice units, 197 MeV), d = 430 |
| `bc_d_lambda_lo` | 640 | assumed | 2033 | - | simplobs.json ch04 C2: step-filter LP degree at Lambda = 0.10 |
| `bc_d_lambda_hi` | 430 | assumed | 2033 | - | simplobs.json ch04 C2: step-filter LP degree at Lambda = 0.15 |
| `bc_n_cfg` | 5 | stated-not-derived | 2033 | simplobs.json ch04 C2 | 5 configurations; gauge variance of nu ~ nu per config |
| `eps_l_2033` | 1e-10 | stated-not-derived | 2033 | app12:154 | eps_l ~ 1e-10 |
| `t_cap_2033` | 1000000000.0 | cited | 2033, codesign | DOE_RFI_2026 | app12:156: ~1e9 hard-op envelope |
| `V_16x4` | 65536 | stated-not-derived | 2033 | app12:22 | '88 LQ at V=16^4' |
| `V_128x4` | 268435456 | stated-not-derived | 2033 | app12:22 | '100 LQ at V=128^4' |

#### Co-design row (post-2033): Tr M^-1 on 24^3x48 (p_f here is the fault convention of every era)

| field | value | provenance | used by | source | note |
|---|---|---|---|---|---|
| `p_f` | 0.1 | stated-not-derived | 2028, 2033, codesign | app12:180 | <=0.1 expected fault per shot (app12:94: p_f/eps_l); report-wide convention by ruling R3 (2026-09-28) |
| `eps_l_codesign` | 1e-12 | stated-not-derived | codesign | app12:180 | eps_l ~ 1e-12 |
| `V_codesign` | 663552 | stated-not-derived | codesign | app12:180 | V ~ 24^3x48 |
| `kappa_codesign` | 1000.0 | stated-not-derived | codesign | app12:180 | kappa ~ 1e3 at physical mass |
| `d_inv_codesign` | 5000.0 | assumed | codesign | - | kappa ln(1/eps_rel) = 1e3 x 5 at the 2033 tolerance: the LP law, extrapolated (the LP runs out of memory near d ~ 5e3); app12:33,182 'd_inv ~= 5e3' (R11; was 1e4) |
| `per_shot_reduction` | 2.6 | stated-not-derived | codesign | app12 1e4-LQ paragraph | 'cut the per-shot T by at least 2.6x' (R3 reprice at D V log2 V, 2026-10-02: 2.567e11 / 1e11; was 10.3x at D^2 V log2 V, E19 P1) |
| `parallelism_codesign` | (1, 4) | cited | codesign | factory.json (T-depth audit, 2026-10-02) ch04 1e4-LQ run | F* = 1-4 on the 92-LQ register; LCU keeps the PREP+SELECT unary-iteration form, so F* does not change |
| `synth_to_poly_2033` | 0.1 | stated-not-derived | 2033, codesign | ruling SYN-BUDGET (H. Lamm, 2026-09-29): EPS_SYN = 1e-2 stays the report-wide default; a circuit may set a tighter budget when its accuracy target needs one. Ch. 4 2028: eps_syn = 0.1 x the polynomial tolerance 1.784e-3 (app12:125) extended to the 1000-LQ and 1e4-LQ rungs (author, E19 P1) | total synthesis error = 0.1 x the polynomial error 0.5 e^-5 = 3.369e-3 on the normalized trace (peak of 1/x = 1/2), i.e. 3.369e-4 = 0.067% of the peak; app12:155,162. Informational: the SCALING count has a unit prefactor, so no hard_ops change |
| `V_production` | 10616832 | stated-not-derived | codesign | app12:180 | physical-mass 48^3x96 |
| `n_cfg_production` | 1000.0 | stated-not-derived | codesign | app12:196 | 'x 1e3 production cfg' |
| `subroutine_calls_box_1e4` | 10000000.0 | stated-not-derived | codesign | app12:197 | '~1e7' (1e4 shots/cfg x 1e3 production cfg); R16 rms, was 4.6e7 |
| `hvp_site_variance` | 1.0 | assumed | codesign | - | referee H2 (2026-10-04): gauge variance per site of the local one-link current term <x\|Gamma M^-1\|x>, order one; not measured (NEEDS_AUTHOR, George). Sets the noise-matched disconnected-HVP shot count |
| `n_qpe_codesign` | 12 | stated-not-derived | codesign | app12:180 | '12 QPE', carried from the 1000-LQ box (app12:152) |
| `n_anc_codesign` | 50 | stated-not-derived | codesign | app12:180 | '~50 ancilla', carried from the 1000-LQ box (app12:152) |

#### Informational: outside every per-shot count

| field | value | provenance | used by | source | note |
|---|---|---|---|---|---|
| `qsp_phase_eps` | 0.0001 | uncited | 2033 | - | no synthesis tolerance stated for the phase loader (app12:100); 1e-4 is the value this informational figure carried before R-TOL (there is no longer a fixed report-wide tolerance; R-TOL would set it from the loader's own rotation count); NOT in the chapter's per-shot T and not printed |

## 4. How the numbers are built

All line numbers refer to `estimates/ch04_hybrid_lqcd.py`. Values below are the model's own outputs (`python -m estimates --chapter ch04 --intermediates`), shown to four figures.

### 4.1 2028: $\log\det W$, free staggered, $V = 4^4$ (`_model_2028`, lines 664-928)

**The matrix (lines 673-683).** `free_W_spectrum` (605-614) diagonalizes Eq. 24 of arXiv:2407.13080 by plane waves, $\lambda = m_0^2 + K^2\sum_\mu \sin^2(2\pi n_\mu/L)$. It gives $\lambda_\mathrm{min} = 0.16$, $\lambda_\mathrm{max} = 4.16$, $\kappa = 26$ and $\log\det W = \sum \log\lambda = 150.5528$. The tests rebuild $M$ and $W$ from plain lists, independently of the model, and check the sparsity, the spectrum and this value (with a dense numpy cross-check). The circuit block-encodes $W/s$, and `be_subnormalization` (617-623) takes the smallest $s$ for which the $O_A$ angles of Eqs. 25-26 exist, $s = \max(2(m_0^2 + DK^2/2),\ 2^c K^2/4) = 4.32$. The QET therefore sees $\lambda_\mathrm{min}/s = 0.0370$. Two terms are restored classically: $V\log s = 374.59$ is added back, and the shift term $V\log(s/\lambda_\mathrm{min})\times\tfrac12 = 421.87$ is subtracted.

**Register (685-692).** 8 site + 4 column + 1 flag + 1 QET signal + 3 workspace + 1 Hadamard-test ancilla = **18 LQ**. The workspace is what the 5-control projector AND needs ($5 - 2 = 3$). On the QME route the Hadamard ancilla is replaced by a copy of the site register (8), an 8-qubit QPE register and one sign qubit, which gives 34.

**Accuracy target (695-700, 744-756).** The polynomial tolerance is the one that guarantees 1% on $\log\det W$:

$$\mathrm{tol} = \frac{0.01\,|\log\det W|}{V\log(s/\lambda_\mathrm{min})} = 1.784\times10^{-3}.$$

$d = 84$ is the smallest even degree whose minimax error on the shifted target $\log|x|/\log(s/\lambda_\mathrm{min}) + \tfrac12$ is within it: the LP error is $1.6993\times10^{-3}$ at 84 and $1.8660\times10^{-3}$ at 82. The test `test_2028_lp_degree_is_minimal` reruns the LP. The resulting error budget on $\log\det W$ is polynomial $\le 1.5055$ (1%), sampling 1% rms at $\varepsilon = \mathrm{tol}$, and synthesis $0.1506$ (0.1%), a combined worst case of 2.1%. Without the shift the same tolerance needs $d = 314$ and $7.20\times10^4$-$1.34\times10^5$ T per shot, above the cap.

**Per-shot T (701-742; the inner function `circuit` is 713-723).** The two ends of the band are two circuits, and each sets its own synthesis tolerance from its own rotation count.

| piece | merged / multiplexed end | literal end (as Figs. 1-2 draw it) | where |
|---|---|---|---|
| $O_c$ controlled shifts per query | 4 shifts x $C^3X$ (3 Toffolis) = 12 Toffolis | 8 shifts x $C^4X$ (5 Toffolis) = 40 Toffolis | `shift2_toffolis` 649-657, `common.mcx_toffoli_count` ($2n-3$) |
| $O_A$ synthesized rotations per query | 2 | 2 | `oa_angles` 626-631, `oa_synthesized_rotations` 640-646 |
| projector phase, $d+1 = 85$ of them | 8 Toffolis + 2 rotations | same | assumptions `projector_toffolis_2028`, `projector_rotations_2028` |

On L = 4 the $+2\hat\mu$ and $-2\hat\mu$ shifts land on the same site, which is what allows the merge. At $s = 4.32$, Eq. 25 gives $\theta_0 = 0$, so the controlled $R_y(\theta_0)$ is the identity and is not synthesized. Only the $\theta_1 = 5.509$ rotations are, which is why both ends carry 2 rotations per query and not the 2-4 that generic angles would need. Each shot then holds $N_\mathrm{rot} = 2d + 2(d+1) = 338$ synthesized rotations, and

$$\varepsilon_\mathrm{syn} = 0.1\,\mathrm{tol} = 1.784\times10^{-4},\quad \varepsilon_\mathrm{rot} = \sqrt{\varepsilon_\mathrm{syn}/N_\mathrm{rot}} = 7.266\times10^{-4},\quad t_\mathrm{rot} = 1.15\log_2(1/\varepsilon_\mathrm{rot}) + 9.2 = 21.19\ \mathrm{T}$$

(`common.eps_rot_for`, `common.t_per_rotation`). With 7 T per Toffoli, a query costs 126.4 or 322.4 T and a projector phase costs 98.38 T, so

$$N_T = d\,T_\mathrm{query} + (d+1)\,T_\mathrm{phase} = 1.898\times10^4\ \text{to}\ 3.544\times10^4,$$

which is 2.82-5.27x below the $10^5$ cap. The primitive breakdown is the literal end and sums to its total:

| primitive | count | T each | T total | status |
|---|---|---|---|---|
| `O_c_controlled_shift_toffolis` | 3360 | 7 | 23520 | COMPILED |
| `O_A_entry_rotations` | 168 | 21.19 | 3560 | COMPILED |
| `projector_phase_toffolis` | 680 | 7 | 4760 | SCALING |
| `projector_phase_rotations` | 170 | 21.19 | 3602 | SCALING |

**Shots, $\epsilon_l$ and wall time (758-771).** The shot count is $\mathrm{Var}/\varepsilon^2$ at Var = 1 and $\varepsilon = \mathrm{tol}$, giving $3.1408\times10^5$ shots (`shots_incoherent`). The required logical error rate is $p_f/N_T$ with $p_f = 0.1$, i.e. $2.821\times10^{-6}$ at the literal end, which binds, and $5.269\times10^{-6}$ at the other. For the T-depth, the $O_A$ and phase rotations form one sequential chain on the flag and signal qubits, and $O_c$ runs underneath:

$$D_T = d\,(2t_\mathrm{rot}) + (d+1)(3 + 2t_\mathrm{rot}) = 84\times42.38 + 85\times45.38 = 7417,$$

where 3 is the depth of the AND tree. $F^* = N_T/D_T = 2.56$-$4.78 < 10$, so the 10 factories behind 1 µs per T cannot be fed. Each shot is charged `corrected_shot_s` (538-540), $\max(N_T\cdot 1\,\mu\mathrm s,\ D_T\cdot 10\,\mu\mathrm s) + t_0 = 74.27$ ms. The run is $2.333\times10^4$ s = 6.48 h, against 1.66-3.10 h at the 1 µs-per-T baseline. The total T is $5.961\times10^9$-$1.113\times10^{10}$. The 6.48 h assumes the 2.6-4.8 factories a shot can feed. On one factory a shot is $N_T\times10\,\mu$s, and the demo takes 16.6-30.9 h (`test_2028_depth_and_factories`). Either way it finishes far inside a year (`factories_for_1yr` < $4\times10^{-3}$).

**The gauged $V = 4^4$ comparison (773-781).** This is a different circuit. Loading links costs $D^2V\log_2V = 32768$ units per query at one T per unit (`scaling_prefactor`, uncited), so the shot is $84\times 32768 = 2.753\times10^6$ T, 27.5x the cap. The two levers are LCU compression (3-10x) and block-encoding one sublattice (2x), together 6-20x. At the top of that range the shot is still $1.376\times10^5$ T, 1.38x the cap.

### 4.2 2033: $\mathrm{Tr}\,M^{-1}$ on $8^3\times16$, and the Banks-Casher row (`_model_2033`, lines 931-1109)

**Register (943-952).** `register_M` (510-512): $\lceil\log_2 V\rceil$ site + $\lceil\log_2(2D+1)\rceil$ column + $\lceil\log_2 N_c^2\rceil$ color + 12 QPE + 50 ancilla = 13 + 4 + 4 + 12 + 50 = **83 LQ**. Without the QPE register the count is 71. Because the register grows as $\log_2 V$, the same formula gives 86 at $16^4$ and 98 at $128^4$.

**Per-shot T (954-965).** $N_T = d_\mathrm{inv}\cdot D V\log_2 V$ with a unit prefactor (`scaling_units_M`, 515-518). $d_\mathrm{inv} = \kappa\ln(1/\varepsilon_\mathrm{rel}) = 100\times5 = 500$ at $\varepsilon_\mathrm{rel} = e^{-5}$ of the peak of $1/x$ (`d_inv_lp`, 569-572; the LP law is checked at $\kappa = 30$ by `test_d_inv_lp_law_odd`). One query is $4\times8192\times13 = 4.260\times10^5$ T, and the shot is $500\times4.260\times10^5 = 2.130\times10^8$ T. The required $\epsilon_l = 0.1/N_T = 4.695\times10^{-10}$.

**Readout (967-982).** The observable is the block-encoding success probability $P = p_0^2\,am\,c$, with $p_0 = 1/2$ the polynomial's value at $\sigma_\mathrm{min}$ (`peak_norm_2033`) and $am = 0.04$. $c$ runs from the free value 0.02749 (`free_condensate`, 521-530, $8^3\times16$, antiperiodic in time; the function is checked against a direct CG solve on $4^3\times4$, and the $8^3\times16$ value is pinned at 0.0275 by `test_free_condensate_matches_a_direct_solve`) to the interacting estimate 0.08. That gives $P = 2.749\times10^{-4}$-$8.0\times10^{-4}$. The shots per configuration at relative error $r$ are $(1-P)/(P r^2)$ (`shots_success_probability`, 533-535):

| tier | configurations, $r$ | shots | T total |
|---|---|---|---|
| first result | 1 at 30% | $1.388\times10^4$-$4.040\times10^4$ | $2.956\times10^{12}$-$8.605\times10^{12}$ |
| campaign | 3 at 20% | $9.367\times10^4$-$2.727\times10^5$ | $1.995\times10^{13}$-$5.808\times10^{13}$ |

A Hadamard test of the same trace has mean $am\,c/2$ and would need $4.340\times10^6$-$3.675\times10^7$ shots for the first result.

**Depth and walls (984-997).** No link-loading circuit is compiled, so the T-depth is carried as a ratio, $D_T = N_T/F^*$. $F^* = 1$ for a multiplexed rotation on the flag qubit and 4 for a unary-iteration load (`parallelism_2033`). A shot takes $\max(N_T\cdot1\,\mu\mathrm s,\ D_T\cdot10\,\mu\mathrm s) + t_0$ = 532.5 s ($F^* = 4$) to 2129.9 s ($F^* = 1$), against 213.0 s at the baseline. The band pairs the cheap ends ($c = 0.08$, $F^* = 4$) and the expensive ends ($c = 0.0275$, $F^* = 1$). The first result takes 85.5 d to 2.727 yr. The campaign takes 1.581-18.41 yr, which is 0.32-3.68x the 5-year horizon; with $F^* = 4$ throughout it is 1.581-4.601 yr. Three parallel unary-iteration subtrees on about 102 idle LQ would give $F^* = 12$ and restore the baseline, 0.632-1.841 yr. This lever is computed but not adopted in the box.

**Banks-Casher row (999-1015).** With $a^{-1} = 0.1973/0.15 = 1.315$ GeV, $\Sigma^{1/3} = 270$ MeV and 4 tastes, $P = 2\Lambda\,(4\,a^3\Sigma)/(\pi N_c)$. The modes per configuration are $\nu = P\cdot V_F = 27.1$ ($\Lambda = 0.15$) and 18.0 ($\Lambda = 0.10$), with $V_F = 24576$ per flavor. The shots are $(1-P)/\big(P\,(0.3^2 - 1/(5\nu))\big)$ over 5 configurations: $1.098\times10^4$-$1.725\times10^4$. The T per shot is the filter degree times one query: $430\times4.26\times10^5 = 1.832\times10^8$ to $640\times4.26\times10^5 = 2.726\times10^8$. The wall time is 58.2 d to 1.490 yr, and the required $\epsilon_l$ is $3.668\times10^{-10}$.

### 4.3 Co-design row: $24^3\times48$, $\kappa \approx 10^3$ (`_model_codesign`, lines 1112-1273)

The same formulas give 20 + 4 + 4 + 12 + 50 = **90 LQ** and $N_T = 5000\times 4\times663552\times\log_2 663552 = 5000\times5.133\times10^7 = 2.567\times10^{11}$ T. At $\epsilon_l = 10^{-12}$ the per-shot budget is $p_f/\epsilon_l = 10^{11}$, so the row is 2.57x over it. A 2.6x reduction (`per_shot_reduction`) leaves $9.872\times10^{10}$. LCU alone (3x) closes the gap and the even-odd split alone (2x) does not. The production lattice $48^3\times96$ needs 94 LQ and $4.956\times10^{12}$ T, 49.6x the budget, beyond the 20x the levers reach. The campaign figures are intermediates only; the row exports no shots or wall time. $10^4$ shots per configuration ($\mathrm{Var}/\varepsilon^2$ at $\varepsilon = 10^{-2}$) on $10^3$ configurations make $10^7$ shots. At 2.86-11.4 days per shot that is $7.820\times10^4$-$3.128\times10^5$ yr, and five years hold 160-640 shots. Lines 1153-1179 cost the disconnected HVP itself. Noise-matching the current loops needs $4.967\times10^{11}$ shots per configuration per flavor and current direction at $24^3\times48$ ($1.490\times10^{12}$ for the three spatial currents), and $7.947\times10^{12}$ per direction at $48^3\times96$ ($2.384\times10^{13}$ for the three currents), at order-one site variance.

## 5. What the paper prints

`resources.json` holds four rows for this chapter (written by `INSTANCE_ROWS`, lines 1310-1330). The paper's `scripts/resources.py` reads them.

| era | label | LQ | T per shot | other keys |
|---|---|---|---|---|
| 2028 | log det W, free staggered V=4^4 (Hadamard test) | 18 | 18978.4-35442.4 | eps_l 2.821e-6, shots 314075.1, lq_qme_route 34, degree 84 |
| 2033 | Tr M^-1, SU(3) staggered V=8^3x16, kappa~1e2 (success-probability readout) | 83 | 212992000 | eps_l 4.695e-10, shots_first 13877.8-40401.1, shots_campaign 93675.0-272707.3, n_cfg_campaign 3 |
| 2033 | Banks-Casher mode number, V=8^3x16, a=0.15 fm (demonstration) | 83 | 183173120-272629760 | eps_l 3.668e-10, shots 10979.8-17247.1 |
| codesign | Tr M^-1, V=24^3x48, kappa~1e3, d_inv~5e3 (before reduction) | 90 | 256659922982.3 | eps_l 1e-12, codesign true |

Table 1.1 of the paper (the chapter list in the overview) prints "18 LQ, 1.9-3.5e4 T" for 2028 and "83 LQ, 2.1e8 T" for 2033. `PUBLISHED` (lines 1287-1307) holds the box values the model is checked against at a 2% tolerance: (18, 1.9e4-3.5e4), (83, 2.1e8) and (90, 2.6e11).

Two sets of tests pin these numbers. The first runs in every layout:

- `test_published_is_the_box_as_printed`: LQ and T of all three eras against `PUBLISHED`.
- `test_instance_rows_carry_model_numbers`: the `resources.json` rows equal the model (18 LQ, QME route 34, degree 84, $3.1408\times10^5$ shots).
- `test_2028_register_18_and_34`, `test_2028_per_query_t`, `test_2028_per_shot_range`: the register, the 126-322 T query, the 98 T phase, $1.898\times10^4$/$3.544\times10^4$ per shot, the 2.8/5.3 margins and the breakdown.
- `test_2028_accuracy_statement`, `test_2028_rotations_and_synthesis_error`, `test_2028_rtol_tolerance_per_circuit`: tol, $d$, the 1%/1%/0.1% budget and $\varepsilon_\mathrm{rot}$.
- `test_2028_shots_wall_time_and_eps_l`: $3.1\times10^5$ shots, $\epsilon_l \lesssim 2.8\times10^{-6}$, $6.0\times10^9$-$1.1\times10^{10}$ T.
- `test_2028_gauged_instance_does_not_fit`: the gauged 4^4 comparison.
- `test_r9_exports_2028_and_2033`: the seven depth and factory exports are present, and $F^* = N_T/D_T$ holds.
- `test_depth_budget_rule`: $p_f/\epsilon_l = 10^7, 10^9, 10^{11}$.
- `test_synthesis_budget_1000lq_and_1e4lq`: the 2033 and co-design synthesis budget.
- `test_2033_d_inv_is_the_lp_minimum`: $d_\mathrm{inv} = 500$.
- `test_free_condensate_matches_a_direct_solve`: $c = 0.0275$.
- `test_2028_model_moves_under_perturbation`: changing $L$, $m_0$, $d$ or the synthesis fraction moves the outputs as it should.

The second set checks the model's numbers and then the chapter's printed strings. In this repository the model checks run and pass, and pytest reports the test as skipped when it reaches the chapter text. A few model checks come after the first chapter-text check and run only in the paper tree: the Banks-Casher $\epsilon_l = 3.7\times10^{-10}$, the amplitude-amplification $6.4\times10^8$ T and $1.6\times10^{-10}$, and the $4.3\times10^5$ per-query T (`test_2033_per_shot_t`), and the $F^* = 10$ perturbation (`test_2033_depth_factories_and_walls`).

- `test_2028_depth_and_factories`: 74 ms, 6.5 h, $F^* = 2.6$-4.8.
- `test_2033_register_addends`, `test_register_scales_as_log_V`: 83, 71, 86, 98.
- `test_2033_per_shot_t`: $2.1\times10^8$, $\epsilon_l \lesssim 4.7\times10^{-10}$ and $3.7\times10^{-10}$.
- `test_2033_p_readout_two_tiers`: $P$, $1.4$-$4.0\times10^4$ and $0.94$-$2.7\times10^5$ shots, $3.0$-$8.6\times10^{12}$ and $2.0$-$5.8\times10^{13}$ T.
- `test_2033_depth_factories_and_walls`: 530-2100 s, 86 d-2.7 yr, 1.6-18 yr.
- `test_2033_banks_casher_demonstration`: $1.8$-$2.7\times10^8$ T, $1.1$-$1.7\times10^4$ shots, 58 d-1.5 yr.
- `test_codesign_24x48_instance`, `test_codesign_production_48x96_after_reduction`, `test_codesign_one_machine_campaign_horizon`: 90 LQ, $2.6\times10^{11}$, 2.6x, $5.0\times10^{12}$, ~50x, 2.9-11 days per shot, $7.8\times10^4$-$3.1\times10^5$ yr.

## 6. Circuit status

From `docs/CIRCUIT_STATUS.md`: COMPILED 2, SCALING 4, CONJECTURE 0, UNSOURCED 0.

| era | primitive | status | what that means here |
|---|---|---|---|
| 2028 | `O_c_controlled_shift_toffolis` | COMPILED | The source draws the circuit (arXiv:2407.13080 Fig. 1). The adder realization (one X on the high bit) and the clean-ancilla ladders are this model's, since the paper gives no T-count. |
| 2028 | `O_A_entry_rotations` | COMPILED | Fig. 2 and Eqs. 25-26 of the same paper. The split of each controlled $R_y$ into 2 rotations + 2 CNOT and the synthesis cost are this model's. |
| 2028 | `projector_phase_toffolis` | SCALING | The standard projector-controlled phase of QET, which the paper does not draw. |
| 2028 | `projector_phase_rotations` | SCALING | The QET phases controlled on the Hadamard-test ancilla, which the paper does not draw. |
| 2033 | `BE_query_M_SU3_staggered_8^3x16` | SCALING | $D V\log_2 V$ per query at a unit prefactor. No compiled gauged block encoding. |
| codesign | `BE_query_M_SU3_staggered_24^3x48` | SCALING | Same. |

The 2028 count is therefore gate-level throughout, with assumptions stated per primitive. The 2033 and co-design counts rest entirely on the scaling of the gauged block encoding.

## 7. Open items and limitations

- **Polynomial degrees are inputs.** $d = 84$ (and 314, 48, 112) and $d_\mathrm{inv} = 500$ come from linear programs carried as `Assumed`. The tests rerun them with scipy. $d_\mathrm{inv} = 5\times10^3$ at $\kappa = 10^3$ is the law $\kappa\ln(1/\varepsilon_\mathrm{rel})$ extrapolated, because the LP runs out of memory near that degree. If the instance moves off the $\lambda_\mathrm{min}/s$ the degrees were computed for, the model flags $d$ as stale (`d_is_for_this_instance`).
- **The gauged block encoding has no circuit.** One T per unit of $D V\log_2 V$ (or $D^2V\log_2 V$ for $M^\dagger M$) is uncited (`scaling_prefactor`). Every gauged query is SCALING, and $F^* = 1$-4 is a structural estimate, not a schedule.
- **The 2028 ancillas and Toffoli ladders are this model's compilation** of the structure arXiv:2407.13080 draws. The paper specifies no adder and no T-counts.
- **The condensate $c$ sets the 2033 shot band.** It lies between the free value 0.0275 and an estimate of 0.08. An ensemble value at $am = 0.04$ on $8^3\times16$ would pin it.
- **What the 2033 result validates.** $P$ is compared with a classical evaluation of the same polynomial on the same configurations. Calling the result $\mathrm{Tr}\,M^{-1}$ itself would need $d_\mathrm{inv} \approx 700$ and 1.4x the per-shot T (`docs/OPEN_ITEMS.md`). Continuum and finite-volume effects are not addressed at this stage.
- **The 2033 walls are depth-bound** at 530-2100 s per shot. The three-subtree layout ($F^* \approx 12$ on about 102 idle LQ) and one step of amplitude amplification are levers that are computed but not adopted. Amplitude amplification would make a $6.4\times10^8$-T shot at $\epsilon_l \approx 1.6\times10^{-10}$ (`test_2033_per_shot_t`, checked in the paper tree).
- **The Banks-Casher inputs** ($\Sigma^{1/3} = 270$ MeV, 4 tastes, 5 configurations, filter degrees 430/640) come from a simplified-observable study. In a 1.2 fm box $P$ is uncertain by about 2x, and the physics is in the $\epsilon$ regime.
- **The post-2033 row is priced at $D V\log_2 V$** with the 2033 block encoding. If link loading needs the extra factor $D$, the cost is $D^2 V\log_2 V$ instead.
- **The disconnected HVP needs far more shots than the row prices.** The $10^4$ shots per configuration price $\mathrm{Tr}\,M^{-1}$ and are a lower bound for the HVP. Noise-matched loops need $1.5\times10^{12}$ per configuration at $24^3\times48$. That count assumes `hvp_site_variance` = 1, which is not conservative, and no coherent variance-reducing estimator is costed.
- **The co-design campaign does not fit the horizon on one machine** even after the 2.6x reduction: $7.8\times10^4$-$3.1\times10^5$ yr. QME would cut the calls 100x, but each circuit would then be 100x deeper.
- **Finite-density reweighting is not a deliverable.** A phase estimator for $\arg\det M$ is not costed.
- **Synthesis of the 2033 phase loader** (`qsp_phase_eps`) is informational and outside every per-shot T.

## 8. Conventions used

The shared conventions are set out in the top-level [README](../../README.md) and `docs/CONTRACT.md`. Logical qubits count the algorithmic register only, without factories. A Toffoli is 7 T (`toffoli_convention = "textbook"`). Rotations use the full synthesis fit $1.15\log_2(1/\varepsilon) + 9.2$, with $\varepsilon_\mathrm{rot} = \sqrt{\varepsilon_\mathrm{syn}/N_\mathrm{rot}}$ from the circuit's own rotation count (R-TOL). The fault budget is $p_f = 0.1$ expected faults per shot. Wall times are serial on one machine at 1 µs per T plus 0.1 ms per shot, with a $10\,\mu$s reaction time per T-depth layer wherever $F^* < 10$. The chapter departs from the report default in one place: its synthesis budget is $0.1\times$ the polynomial tolerance (`synth_to_poly_2028`, `synth_to_poly_2033`) rather than `EPS_SYN = 1e-2`, because the accuracy claim rests on that tolerance.
