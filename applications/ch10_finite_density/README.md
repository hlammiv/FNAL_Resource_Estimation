# Chapter 10: QCD at finite baryon density

Model: `estimates/ch10_finite_density.py`. Tests: `estimates/tests/test_ch10_finite_density.py`. Paper chapter: `applications/app07_finite_density.tex` (not part of this repository). Line numbers below refer to the model file unless another file is named.

## 1. What is estimated

The chapter targets the Kogut-Susskind Hamiltonian of lattice QCD with staggered quarks at nonzero baryon chemical potential, H(mu_B) = H_QCD - mu_B N_B. Gauge links are digitized on the discrete subgroup Sigma(36x3) of SU(3) (108 elements, 8 qubits per link in the compiled encoding of arXiv:2405.05973) and the quarks are mapped to qubits by Jordan-Wigner. The physics goal is the equation of state on a 5 x 5 grid, T in [100, 200] MeV and mu_B in [0, 600] MeV: the pressure p, the baryon density n_B, the susceptibilities chi_2 to chi_4 and the ratio chi_4/chi_2, the lattice counterpart of the net-proton kappa sigma^2 measured in the RHIC Beam Energy Scan.

The model prices three instances, one `resources.json` row each.

- **2028 benchmark.** 2+1D SU(3), V = 2^2, one staggered field (two degenerate light tastes). The validation observable is the ground-state energy E_0(N_B) in two or three baryon-number sectors at T -> 0, and from it the onset chemical potential mu_c = E_0(N_B+1) - E_0(N_B). The state is prepared by an adiabatic ramp; every Trotter term conserves N_B, so the ramp stays in its sector. The model prices one second-order Trotter step, which is a lower bound on the executed ramp.
- **2033 target.** 3D SU(3), V = 2^3, three staggered fields (N_f = 4+2 tastes, since a Hamiltonian has no fourth root). The thermal state proportional to exp(-(H - mu_B N_B)/T) is prepared by a QSVT imaginary-time filter applied to a reference state and boosted by amplitude amplification, the same construction as the thermal route of Chapter 9 (`ch09_chiral_gauge.py`). Each shot reads N_B; its moments give n_B, chi_2, chi_3, chi_4. Two tiers are costed: a first result (one temperature row, 30% on n_B and chi_4/chi_2) and the campaign (the full grid at an absolute 0.1 on chi_4/chi_2).
- **Post-2033 stretch (co-design row).** Sigma(72x3) (9 qubits per link) on V = 4^3 with three staggered fields. Only the register is computed; no T-count is given because no Sigma(72x3) Fourier transform circuit exists yet.

## 2. How to run

From the repository root:

```
python -m estimates --chapter ch10                     # this chapter only
python -m estimates --chapter ch10 --intermediates     # also every named intermediate of each era
python -m pytest estimates/tests/test_ch10_finite_density.py -q
python applications/ch10_finite_density/dump_assumptions.py   # the table of section 3, as Markdown
```

`--chapter ch10` recomputes the model and writes nothing. It prints (i) the 74 assumptions with value, provenance tag and source (6 Cited, 33 Assumed, 35 Stated-not-derived); (ii) the three rows the chapter writes to `resources.json`; (iii) for each era, the model's logical qubits and per-shot T beside the numbers the chapter box prints, with `[ok]` when they agree within the 10% tolerance, then shots, wall time, required logical error rate, the depth and factory exports, and the primitive-by-primitive T breakdown with each primitive's circuit status. The test file runs in under a second: 44 pass and 4 skip when the paper source is not next to the package (the skipped tests check that the chapter text quotes the model's numbers).

## 3. Assumptions

Every input of the model is a field of the `Assumptions` dataclass (lines 215-484) and carries a provenance tag:

- **Cited**: a published source gives the value; the source column holds its bibliography key.
- **Assumed**: a working assumption of this package or of the chapter, with its reasoning in the note.
- **Stated-not-derived**: the chapter states the number and the model uses it without deriving it; the source column gives the chapter location as the field records it (`app07:line`). Those line numbers may not match the current chapter file line for line; the quoted phrase in the note is the stable anchor. This holds for the rows that feed the headline. The comparison-route records and `n_trotter_2028_original` and `c_mix` quote wording the current chapter no longer contains; they are kept as records and do not enter the boxes.

The table below is the output of `dump_assumptions.py`, which reads the live dataclass; only the grouping is added by the script. Values, sources and notes are the code's own strings, printed unedited; some sources point to internal change-log files and some notes carry internal labels, which are not needed to follow the numbers. The last group holds inputs of comparison routes that the model computes and reports but that do not enter the headline numbers (section 4.9).

74 fields: Assumed 33, Cited 6, Stated-not-derived 35

#### Shared: lattice theory and encoding

| name | value | provenance | source | note |
|---|---|---|---|---|
| `nc` | 3 | Stated-not-derived | app07:91 | N_c=3 is the color count |
| `group` | S36x3 | Stated-not-derived | app07:138,162 | Sigma(36x3), 8 qubits/link (compiled encoding, R1); both boxes |
| `ham` | KS | Cited | arxiv_2405_05973 tab:primcost | Kogut-Susskind: plaquette + E^2 are the chapter's families (a), (b) at app07:96 |

#### Shared: synthesis, Toffoli and hop readings

| name | value | provenance | source | note |
|---|---|---|---|---|
| `eps_syn` | 0.01 | Assumed | TRACKED_CHANGES.md R-TOL (H. Lamm, 2026-09-29) | R-TOL: total synthesis error per shot, report-wide; each circuit sets eps_rot = sqrt(eps_syn / N_rot) from its own rotation count (replaces the Stated fixed 1e-4 of app07:96) |
| `rot_synthesis_paper` | rus | Cited | arxiv_2405_05973,arxiv_2408_00075; apply_log/r18_core.md | E20 ruling (2026-10-01, supersedes R2's papers' convention): every gauge rotation at the full fit '1.15 log2(1/eps) + 9.2', n_rot = the paper's log coefficient / 1.15, the papers' Toffoli constants kept; eps from the R-TOL incoherent (randomized) split. groups.PrimitiveCost.t() prices in it; t_papers() is the slope-only legacy record |
| `toffoli_convention` | textbook | Cited | arxiv_2405_05973 | 7 T per Toffoli, as in the papers' tables (PAPER_COSTS.md Sec. 4) |
| `hop_synthesis` | rus | Assumed | apply_log/r17_core.md; FermionPrimitives_unpub resources.tex:7-10 | every hop and mass rotation at the full fit 1.15 log2(1/eps) + 9.2 (r17 rule 7; the draft prices 9.2 + 1.15 L too) |
| `hop_undo` | True | Assumed | groups.fermion_hop_counts; E21 (1), apply_log/r19_core.md | V_g^dag on both sites after the diagonal hop (ESTIMATED; the draft draws it, PD/section_hamiltonian.tex:86-91, but its 2 C^G counts only the 2 forward applications, PD/section_su2_diag.tex:238, section_resources.tex:15) |
| `hop_share` | link | Assumed | groups.fermion_hop_counts; apply_log/r19_core.md | color squish + parity computed once per link and held through V(x), V(y), hop, V^dag(x), V^dag(y) (E21 (1) 'make the switch'); 'draft' (recompute per frame application) is the conservative sensitivity, reported beside it |
| `hop_mcx` | mbu | Assumed | groups.mcx_toffolis | C^nX ladders at n-1 Toffolis (measurement-based uncompute, report rule) |
| `hop_phasing` | hwp | Assumed | groups.fermion_hop_counts; arxiv_1709_06648 | the 2 N_c N_stag equal-angle Rz of the diagonal hop as one HWP group (derived here) |

#### Shared: wall time and logical error

| name | value | provenance | source | note |
|---|---|---|---|---|
| `t_gate_s` | 1e-06 | Stated-not-derived | app07:151,175 | 1 us per T-gate |
| `shot_overhead_s` | 0.0001 | Assumed | R17 g_rate; ch08_baryogenesis.py | per-shot overhead ~0.1 ms (register init, final readout, decode), report rule; as Ch. 8 shot_overhead_s |
| `faults_per_shot` | 0.1 | Assumed | TRACKED_CHANGES.md RULINGS round A | R3 ruling: required eps_l = 0.1 expected faults per shot, report-wide |
| `eps_l_floor` | 1e-08 | Cited | DOE_RFI_2026 | RFI floor, app07:171 |

#### 2028 benchmark: 2+1D SU(3), V = 2^2, ground state at mu_B > 0

| name | value | provenance | source | note |
|---|---|---|---|---|
| `d_2028` | 2 | Stated-not-derived | app07:134,138 | 2+1D SU(3) ground state |
| `L_2028` | 2 | Stated-not-derived | app07:139 | V=2^2 |
| `n_stag_2028` | 1 | Stated-not-derived | app07:139 | 1 staggered field (2 degenerate light tastes) |
| `n_anc_2028` | 100 | Stated-not-derived | app07:93,146 | ~100 ancilla as a sized budget (ruling ch10-ancilla-contents A): primitive workspace 29-55 (groups.py S36x3 + the r17 hop, its 24-qubit squish scratch held per link (r19), and the mass; phasing groups at k - w(k) ancilla, E26) + 1 Hadamard-test qubit accounted, rest margin |
| `t_cap_2028` | 100000.0 | Cited | DOE_RFI_2026 | 2028 per-shot hard-op cap; app07:134,147 |
| `n_trotter_2028` | 1 | Stated-not-derived | app07:103,143 | one second-order Trotter step is priced; it is 4.0x the 1e5 cap (r19), so floor(cap / c_BE) = 0 and no step fits |
| `n_trotter_2028_original` | 10 | Stated-not-derived | app07:103,148 | 'the ten-step adiabatic ramp originally scoped ... does not' fit; re-scoping the validation target is the authors' call |
| `shots_2028` | 2000.0 | Stated-not-derived | app07:150 | ~2e3 shots, single (T,mu_B) point; no eps given |
| `sigma2_2028` | 1.0 | Assumed | r25 task rulings; app07:145 | per-shot variance bound sigma^2 <= 1 of the normalized E_0 estimator; 2e3 shots -> 2.2% rms per sector (author approved 2026-10-02) |
| `n_sectors_2028` | (2, 3) | Stated-not-derived | app07:141 | 'ground-state energies E_0(N_B) in 2--3 sectors' |
| `t_depth_2028` | (26719.289742577264, 40488.90205167576) | Assumed | factory.json Ch. 10 run 1 (scratchpad factories/ch10_depth_v2.py), r25 | T-depth of the one-step 2028 shot: low = T-par color frame (commuting Euler-factor rotations in one layer), 8 link rounds (its ~60 parity ancilla per hop fit the ~100 budget only at 8 rounds; 4 rounds, 1.43e4, is over budget); high = as compiled, 2d = 4 color classes. Fully serial 183,854 (F* 2.2) is the pessimistic bound, not priced |

#### 2033 target: 3D SU(3), V = 2^3, register and step count

| name | value | provenance | source | note |
|---|---|---|---|---|
| `d_2033` | 3 | Stated-not-derived | app07:158,162 | 3D SU(3) |
| `L_2033` | 2 | Stated-not-derived | app07:162 | V=2^3 |
| `n_stag_2033` | 3 | Stated-not-derived | app07:163 | 3 staggered fields (N_f=4+2 tastes) |
| `n_anc_2033` | (250, 750) | Stated-not-derived | app07:93,168 | 250-750 for the 2033 route: BE workspace and the parity ancilla of parallel hops (was the Gibbs route's BE/bath register; relabelled 2026-10-05) |
| `t_ref_2033` | 1000000000.0 | Cited | DOE_RFI_2026 | 2033 per-shot hard-op reference; app07:105,158,169 |
| `beta_h_2033` | 100.0 | Stated-not-derived | app07:105 | beta\|\|H\|\| ~ 1e2 for the KMS Gibbs preparation |
| `c_mix` | 20 | Stated-not-derived | app07:105 | t_mix = c beta\|\|H\|\| with c ~ 20; 'least controlled constant in the 2033 costing'. Referee F2 (2026-10-04): read as c sweeps of the jump set (t_mix ~ c at per-jump rate 1) |
| `t_depth_2033` | (98454105.32846098, 315713175.1043646) | Assumed | factory.json Ch. 10 run 2 (scratchpad factories/ch10_depth_v2.py), r25 | T-depth of the 2e3-step Gibbs shot: low = T-par color frame, 12 link rounds (its ~170 parity ancilla per hop fit the 750 band only at ~3 concurrent hops; 6 rounds, 5.23e7, needs ~880); high = as compiled, 12 link rounds (250-ancilla end). Fully serial 3.89e9 (F* 1.4) is the pessimistic bound, not priced |

#### 2033 target: QSVT/TPQ thermal state (the priced route)

| name | value | provenance | source | note |
|---|---|---|---|---|
| `delta_beta` | 0.0001 | Assumed | ch09_chiral_gauge.py delta_beta; author ruling 2026-10-05 | QSVT filter tolerance delta_beta, as Ch. 9 (app06 delta_beta = 1e-4): d_beta = round(sqrt(1e2 ln 1e4)) = round(30.35) = 30 |
| `filter_success_prob` | 0.01 | Assumed | app07:107,193 ('overlap >~1e-2'); app06 route (iv) | success probability \|\|e^{-beta H/2}\|psi>\|\|^2 of the filter on the reference (amplitude 0.1, the value Ch. 9 presumes); a typical (Haar) reference gives Z e^{beta E_0}/D instead, exponentially small (referee C1) |
| `aa_call_rule` | 2k+1 | Assumed | ch09_chiral_gauge.py aa_call_rule | k rounds of amplitude amplification call the filter 2k + 1 times (as Ch. 9, verifier 2026-10-04); k = ceil(pi / (4 arcsin a)) = 8 at a = 0.1 |
| `ref_rot_per_qubit` | 1 | Assumed | derived here; app06 route (iv) | reference preparation: a product state, one synthesized rotation per system qubit (one-link / one-site marginals), applied with each filter call |
| `qsvt_query_cost` | (2, 4) | Stated-not-derived | app07:107 | each block-encoding query at 2-4x the Trotter-step cost c_step |

#### 2033 target: grid, accuracy targets and shots

| name | value | provenance | source | note |
|---|---|---|---|---|
| `grid_T_mev` | (100, 125, 150, 175, 200) | Assumed | app07:objectives, fig:phase_diagram | 5 evenly spaced temperatures in [100,200] MeV |
| `grid_mu_mev` | (0, 150, 300, 450, 600) | Assumed | app07:objectives, fig:phase_diagram | 5 evenly spaced mu_B in [0,600] MeV |
| `n_T` | 5 | Stated-not-derived | app07:56,165 | 5 temperatures, T in [100,200] MeV |
| `n_mu` | 5 | Stated-not-derived | app07:56,165 | 5 chemical potentials, mu_B in [0,600] MeV |
| `eps_stat` | 0.1 | Assumed | app07:Nshot_fd | absolute target \|delta(chi4/chi2)\| <= 0.1 per grid point, 10% at the hadron-gas value 1 and well defined near zeros of chi4 (referee F1) |
| `eps_chi4_chi2_target` | 0.1 | Stated-not-derived | app07:167 | '10% on chi4/chi2' |
| `eps_pressure_target` | 0.05 | Stated-not-derived | app07:56,167,224 | '5% on p', not separately costed: chi4 at eps=0.1 sets the shots, p inherited (ruling ch10-pressure-5pct a) |
| `eps_first_result` | 0.3 | Assumed | ruling R1, H. Lamm 2026-10-02 | R1: the minimal first result at 30% statistical error |
| `nb_toy_k2_priced` | 0.01 | Assumed | referee F1; scratchpad referee/ch10_f1.py | kappa_2(N_B) at mu_B = 0 of the priced toy: Skellam in baryons, the same at every T (the true kappa_2 at V = 2^3 is open) |
| `nb_toy_k2_band` | (0.01, 0.1) | Assumed | referee F1 | baryon-Skellam kappa_2(0) band carried beside the priced value |
| `nb_toy_k2_mid` | 0.05 | Assumed | referee F1 | middle of the band |
| `nb_toy_quark_k2` | 0.05 | Assumed | referee F1 | quark-Skellam toy, kappa_2(N_B) at mu_B = 0 (chi4/chi2 = 1/9) |
| `fr_rw_pairs` | (('baryon', 0.01, (125, 600)), ('baryon', 0.05, (175, 600)), ('baryon', 0.1, (175, 600)), ('quark', 0.05, (200, 600))) | Assumed | referee F1 | best two-ensemble pair (mu_B in MeV) for the T = 150 MeV first-result row, found on a 25-MeV grid (scratchpad referee/run2.py); each toy takes min(direct, reweighted) |
| `grid_rw_pairs_k2lo` | ((100, (150, 600)), (125, (125, 600)), (150, (125, 600)), (175, (125, 600)), (200, (150, 600))) | Assumed | referee F1; scratchpad referee/run2.py | best two-ensemble pair per T for the grid at kappa_2(0) = 0.01 (record) |
| `grid_rw_pairs_k2mid` | ((100, (200, 600)), (125, (175, 600)), (150, (175, 600)), (175, (175, 600)), (200, (175, 600))) | Assumed | referee F1; scratchpad referee/run2.py | best two-ensemble pair per T for the grid at kappa_2(0) = 0.05 (record) |
| `campaign_horizon_yr` | 5 | Stated-not-derived | app07:201 | 'Campaign horizon 5 years'; the 5x5 grid is 45.1 yr serial, so 5 years hold 2.86e4 shots and the full grid needs a 9.0x algorithmic cost reduction |

#### Post-2033 stretch: Sigma(72x3), V = 4^3 (register only)

| name | value | provenance | source | note |
|---|---|---|---|---|
| `group_stretch` | S72x3 | Stated-not-derived | app07:197 | Sigma(72x3), 9 qubits/link (compiled encoding, R1) |
| `d_stretch` | 3 | Stated-not-derived | app07:197 | production 3D |
| `L_stretch` | 4 | Stated-not-derived | app07:197 | V=4^3 |
| `n_stag_stretch` | 3 | Stated-not-derived | app07:197 | 3 staggered fields |
| `n_anc_stretch` | (250, 750) | Stated-not-derived | app07:93,197 | the 250-750 BE/bath ancilla of the 2033 route, carried V-independently (ruling ch10-stretch-ancilla a) |

#### Comparison routes (computed, not in the headline numbers)

| name | value | provenance | source | note |
|---|---|---|---|---|
| `c_be_2033_superseded` | 500000.0 | Stated-not-derived | app07:96,105,169 (before R-HOP) | ~5e5 T/step as printed before R-HOP (gauge 3.1e5 + an unsourced hop allowance); superseded by the derived c_BE, kept as a record only |
| `gibbs_2x2_min_steps` | 14 | Stated-not-derived | app07:103 | 'Gibbs sampling on V=2^2 needs >~14 sampler steps': the earlier '>~1e6 T per shot' read at the then 7.25e4 T/step (1e6 / 72,515.5 = 13.8); lower bound, no beta\|\|H\|\| given (item 10) |
| `gibbs_2x2_t_per_shot_superseded` | 1000000.0 | Stated-not-derived | app07:103 (before r17) | '>~1e6 T-gates per shot' as printed before r17; record only |
| `gibbs_beta_h_2028` | 100.0 | Assumed | referee F2; app07:105 | beta\|\|H\|\| for a Gibbs sample at V=2^2: the 2033 value, none is stated for the 2028 instance |
| `gibbs_sweeps` | (10, 40) | Assumed | referee F2; app07:105 | sweeps of the jump set per Gibbs sample, low/high; central = c_mix (20, Stated). One sweep = \|A\| sampler steps |
| `gibbs_eps_channel` | 0.01 | Assumed | referee F2; same budget as common.EPS_SYN | channel-approximation error per Gibbs sample, split evenly over its sampler steps; sets the OFT and B time truncations |
| `gibbs_ctrl_toffoli_per_rot` | 1 | Assumed | referee F2, derived here | controlled Trotter step: each synthesized rotation is controlled by one Toffoli onto an AND ancilla (7 T, MBU uncompute), the control qubit's own phase aggregated into one extra rotation per step; Toffoli/Clifford compute-uncompute pairs need no control |
| `gibbs_structure` | ('one OFT, no inverse, no coherent term', 'OFT and its inverse + coherent term B and its inverse', 'central x K, the per-unit-time segment order of black-box Lindbladian simulation, K = ln(t/eps)/lnln(t/eps)') | Assumed | arxiv_2311_09207 eqs. OpOFT, mainBDef, b1, b2; Thm. I.2 | controlled evolution per sampler step, low/central/high (Gaussian weight, sigma_E = 1/beta; OFT window e^{-t^2/beta^2}, \|t\| <= beta sqrt(ln 1/eps_u); B: \|t\| <= ln(1/eps_u)/(2 pi), \|t'\| <= sqrt(ln(1/eps_u)/4)). Time register in sign-magnitude form: the palindromic second-order step has S(-dt) = S(dt) with negated angles, so the sign is CNOT conjugation (Clifford) and a window [-T, T] costs controlled evolution T per side (offset binary: 2T, central +40%). Low drops B, so it is a floor without exact KMS detailed balance. Central = one query per unit time: not given by Thm. I.2 (which gives x K); supported up to an unpriced O~(1) constant by the discrete-time channel of arXiv:2405.20322 |
| `chain_reuse_gain` | (3, 30) | Assumed | shot audit Ch. 10 (readout_levers 4) | continuing one Gibbs chain with a non-demolition N_B readout: cost per sample ~ tau_int ~ 2 t_rel instead of t_mix; ~10x central, 3-30x (conditional on fault dissipation and mid-circuit Gauss checks) |
| `chain_reuse_steps_per_tau` | 200 | Assumed | shot audit Ch. 10 (readout_levers 4) | Gibbs steps per N_B autocorrelation time at the 10x central gain (2e3 / 10); 0.1 faults in that window sets eps_l |
| `eps_qsp` | 0.001 | Stated-not-derived | app07:107 (before 2026-10-05) | record: d_QSP was evaluated at eps=1e-3 |
| `R_warm_stated` | 10 | Stated-not-derived | app07:107,192 (before 2026-10-05) | record: 'R ~ 10' rounds with 2 queries per degree (20 filter-equivalents; now 2k + 1 = 17) |
| `qsvt_rail_stated` | (3000000000.0, 6000000000.0) | Stated-not-derived | app07:107 (before 2026-10-05) | record: '~3-6e9 T' (2.86e9-5.81e9 derived at R = 10, 2 queries, d = 26.3) |
| `d_qsp_stated` | 26 | Stated-not-derived | app07:106 (before 2026-10-05) | record: 'd_QSP ~= 26' (26.3 at eps = 1e-3) |

## 4. How the numbers are built

The chain is register, then rotation tolerance, then the cost of one Trotter step, then the cost of one shot, then shots, then wall time. Every number quoted here is a key of `model(a, era).intermediates` (shown in backticks) or a line of the primitive breakdown, and, where a test pins it, the test is named.

### 4.1 Logical qubits (`_geometry`, lines 491-504)

N_q = d L^d q_G + N_stag N_c L^d + N_anc: links times qubits per link, plus one qubit per color per field per site, plus an ancilla budget.

| instance | gauge | fermion | ancilla | total (`lq_total`) | test |
|---|---|---|---|---|---|
| 2028 | 8 links x 8 = 64 | 1 x 3 x 4 = 12 | 100 | 176 | `test_2028_register_components` |
| 2033 | 24 links x 8 = 192 | 3 x 3 x 8 = 72 | 250-750 | 514-1014 | `test_2033_register_components` |
| stretch | 192 links x 9 = 1728 | 3 x 3 x 64 = 576 | 250-750 | 2554-3054 | `test_codesign_register` |

The 2028 budget of 100 ancillas is itemized (lines 937-951; `ancilla_workspace`): U_FFT 8, U_Tr 7, U_inv 4, U_mul 2, the hop's Hamming-weight-phasing group HWP(6) 4, its C^7X ladder 5, the 24-qubit color-squish scratch of the hop, and the mass HWP(3) 1. With serial reuse the workspace is 29 (the squish scratch is held through the hop, so it sits beside the ladder); without reuse it is 55. With one Hadamard-test qubit the algorithmic register is 106-132 of the 176 (`algorithmic_register`; `test_2028_ancilla_budget_accounting`); the rest is margin. The 2033 band of 250-750 (block-encoding workspace and the parity ancillas of hops run in parallel) is Stated and not itemized.

### 4.2 Rotation tolerance (`_rot_per_step`, lines 546-562; `_eps`, lines 690-692)

Synthesized rotations are counted from the gate tables, independent of the tolerance. Each circuit then sets its own tolerance eps_rot = sqrt(0.01 / N_rot) (`common.eps_rot_for`), so the total synthesis error per shot is 0.01, and every rotation costs 1.15 log2(1/eps_rot) + 9.2 T (`common.t_per_rotation`). No iteration is needed because N_rot does not depend on eps.

Rotations per Trotter step (`n_rot_per_step`): 10,908 at 2028 (3,708 gauge + 7,192 hop + 8 mass) and 75,872 at 2033 (11,208 + 64,632 + 32). The gauge part per link is U_Tr (d-1)/2 x 7 + 2 U_FFT x 102 + U_phi x 256: 463.5 at d = 2 and 467 at d = 3.

| circuit | N_rot per shot | eps_rot | T per rotation | test |
|---|---|---|---|---|
| 2028, one step | 10,908 | 9.575e-4 | 20.733 | `test_rtol_rotation_count_and_tolerance` |
| 2033, 2e3 plain steps (step price quoted in the chapter) | 1.517e8 | 8.118e-6 | 28.647 | same |
| 2033 QSVT shot, low end | 77,394,463 | 1.137e-5 | 28.089 | `test_2033_hard_ops_product`, `test_2033_qsvt_matches_ch09_construction` |
| 2033 QSVT shot, high end | 154,783,903 | 8.038e-6 | 28.663 | same |

`test_rtol_rotation_count_matches_the_breakdown` checks that no synthesized rotation is missed: the change of the T-count between two tolerances equals 1.15 N_rot log2(eps_1/eps_2) to 1e-9.

### 4.3 One Trotter step: c_step = gauge + hop + mass (`_t_per_step`, lines 540-543)

**Gauge terms** (`_per_link`, lines 507-516; `groups.magnetic_per_link` and `groups.electric_per_link`, `groups.py` lines 407-428). Per link per step the Kogut-Susskind multiplicities of arXiv:2405.05973 (tab:primcost; `groups.PRIMCOST`, `groups.py` lines 132-135) are U_inv 3(d-1), U_mul 6(d-1), U_Tr (d-1)/2 for the plaquette, and two Fourier transforms (the FFT of arXiv:2408.00075) plus one diagonal phase U_phi for the electric term. Each primitive costs t_const + n_rot (1.15 log2(1/eps) + 9.2) (`groups.PrimitiveCost.t`), with the papers' Sigma(36x3) values U_inv 119 + 0, U_mul 308 + 0, U_Tr 378 + 7, U_FFT 532 + 102, U_phi 0 + 256 (Toffolis at 7 T in t_const; n_rot read as the papers' log coefficient over 1.15).

| | magnetic T/link | electric T/link | gauge T/step | test |
|---|---|---|---|---|
| 2028 (eps 9.575e-4) | 2,466.6 | 10,601.1 | 8 x 13,067.6 = 104,541.0 | `test_per_link_costs_are_the_papers_numbers`, `test_2028_per_step_cost_derived_from_the_papers` |
| 2033 (eps 8.118e-6) | 4,988.5 | 14,241.6 | 24 x 19,230.2 = 461,523.9 | `test_2033_c_be_consistent_with_compiled_floor` |

The electric term is 81% of the 2028 gauge cost (`electric_share`). Without the FFT it would be 633 (2028) to 651 (2033) times larger (`naive_over_fft`).

**Staggered hop** (`_hop_link`, lines 524-531; `groups.fermion_hop_counts` and `groups.hop_link_cost`, `groups.py` lines 787-941). The gate counts come from an unpublished draft on gauge-covariant fermion primitives (bibliography key `FermionPrimitives_unpub`). A hop on a link applies a color frame V on both sites, a diagonal hop, and V^dag on both sites. Per link per step, at the chapter's readings (`hop_undo` True, `hop_share` "link", `hop_mcx` "mbu", `hop_phasing` "hwp"):

| piece | 2028 (N_stag = 1) | 2033 (N_stag = 3) | how it is counted |
|---|---|---|---|
| color squish and parity | 2,296 Toffoli = 16,072 T | same | 2 x 759 squish + 778 parity/flag Toffolis, computed once per link and held through the hop |
| color rotations | 656 rotations = 13,600.7 T | 1,968 = 56,377.4 T | 164 per frame application x 4 applications (V and V^dag on both sites) x N_stag |
| hop squish and flags | 292 Toffoli = 2,044 T | same | 2 x (46 + 100), compute and uncompute |
| diagonalizers | 60 x (2 T + 4 rot) = 5,095.9 T | 180 x (2 T + 4 rot) = 20,985.9 T | 2 x 10 classes x 3 colors x N_stag |
| phasing | HWP(6): 3 rot + 4 Toffoli = 90.2 T | HWP(18): 5 rot + 16 Toffoli = 255.2 T | the 2 N_c N_stag equal-angle Rz as one Hamming-weight-phasing group |
| **per link** (`hop_t_per_link`) | 2,592 Toffoli + 120 T + 899 rot = **36,902.7** | 2,604 + 360 + 2,693 = **95,734.4** | |
| per step (`hop_t_per_step`) | x 8 = 295,221.9 | x 24 = 2,297,626.8 | |

The color frame (first two rows) is 80% of the 2028 hop and 76% of the 2033 hop (`hop_colour_frame_share`). The model also carries the alternative readings at the same tolerance: recomputing the squish at every frame application (`hop_t_per_link_share_draft`, 85,118.7 and 148,038.4 T/link), dropping the frame undo (`hop_t_per_link_no_undo`, 30,102.4 at 2028), and the draft's own printed formula (`hop_t_per_link_fp_printed`, 32,297.9 at 2028). These are pinned in `test_2028_per_step_cost_derived_from_the_papers`, `test_2033_c_be_consistent_with_compiled_floor` and `test_r19_hop_reading_is_e21`.

**Staggered mass** (`_mass_site`, lines 534-537; `groups.staggered_mass_site`, `groups.py` lines 944-958). One Rz per color and field on each site, the N_c N_stag equal angles grouped as one HWP group. Because N_B commutes with H, the mu_B N_B term adds the same angle to every copy on a site and folds in at no cost. HWP(3) is 2 rotations + 1 Toffoli = 48.47 T per site (x 4 sites = 193.9 T); HWP(9) is 4 rotations + 7 Toffolis = 163.6 T per site (x 8 = 1,308.7 T).

**The step** (`c_be_t_per_step`): at 2028, 104,541.0 + 295,221.9 + 193.9 = **399,956.8 T**, 4.0 times the 1e5 cap, so no step fits (`steps_under_cap` = 0). At 2033 and the tolerance of 2e3 plain steps, 461,523.9 + 2,297,626.8 + 1,308.7 = 2,760,459.4 T, of which the hop is 83% (`hop_share_of_c_be`). Inside the QSVT shot the same step is re-priced at that circuit's tolerance: 2,718,083 T at the low end and 2,761,708 T at the high end (`qsvt_c_step`).

### 4.4 2028 per-shot T and breakdown (`_model_2028`, lines 885-1056)

The per-shot count is the one priced step, 399,956.8 T (`hard_ops`; `test_2028_hard_ops_band_and_breakdown`). Breakdown as printed by `--chapter ch10`:

```
magnetic_U_inv                  24 x      119 =       2856  COMPILED
magnetic_U_mul                  48 x      308 =  1.478e+04  COMPILED
magnetic_U_Tr                    4 x    523.1 =       2093  COMPILED
electric_U_FFT                  16 x     2647 =  4.235e+04  COMPILED
electric_U_phi                   8 x     5308 =  4.246e+04  COMPILED
hop_colour_squish_and_parity     8 x 1.607e+04 =  1.286e+05  SCALING
hop_colour_rotations             8 x 1.36e+04 =  1.088e+05  SCALING
hop_squish_and_flags             8 x     2044 =  1.635e+04  COMPILED
hop_diagonalizers                8 x     5096 =  4.077e+04  COMPILED
hop_phasing_hwp                  8 x     90.2 =      721.6  SCALING
mass_hwp_rotations               8 x    20.73 =      165.9  SCALING
mass_hwp_toffolis                4 x        7 =         28  COMPILED
```

The executed circuit is the whole ramp plus the energy readout. A ten-step ramp priced as its own circuit (N_rot 109,080, eps_rot 3.03e-4) is 4.21e6 T (`ramp_original_t`; `test_2028_per_step_cost_derived_from_the_papers`). The ramp length, the preparation error and the grouped readout of H are not priced.

### 4.5 2033 per-shot T: the QSVT thermal filter (`qsvt_thermal`, lines 628-672; `qsvt_band`, lines 675-680)

- Degree: d_beta = round(sqrt(beta||H|| ln(1/delta_beta))) = round(sqrt(100 ln 1e4)) = round(30.35) = 30 (`d_beta`). The polynomial is in e^{-beta H/2} alone: N_B commutes with H, so the e^{beta mu_B N_B/2} weight rides on the reference's sector weights and the degree does not depend on mu_B.
- Amplification: amplitude a = sqrt(0.01) = 0.1, k = ceil(pi / (4 arcsin a)) = 8 iterations, 2k + 1 = 17 filter calls (`qsvt_filter_calls`), 17 x 30 = 510 block-encoding queries (`qsvt_queries`).
- A query costs 2-4 Trotter steps (Stated), so a shot is 1,020-2,040 step equivalents (`qsvt_step_equivalents`).
- N_rot = step equivalents x 75,872 + 17 x 264 (product reference, one rotation per system qubit per call) + (17 x 31 + 8) (QSVT phases and one phase per amplification iteration), which sets each end's tolerance (section 4.2).
- T = step equivalents x c_step(eps) + (reference + phase rotations) x T per rotation + 8 reflections x 1,014 Toffolis x 7 T. Each reflection about the start state is priced on the full high-end register of 1,014 qubits at both ends (n - 1 Toffolis per C^nZ; an upper bound at the 514-LQ end, as the chapter's 'at most 1014 Toffolis' says). N_B is read in the computational basis with no T.

| end | filter | reference | phases | reflections | **T per shot** |
|---|---|---|---|---|---|
| low (2 steps per query) | 2,772,445,126 | 126,061 | 15,027 | 56,784 | **2,772,642,999** |
| high (4 steps per query) | 5,633,883,942 | 128,642 | 15,335 | 56,784 | **5,634,084,702** |

Everything besides the filter is under 1e-4 of the shot (`qsvt_overhead_share`, 7.1e-5 and 3.6e-5). The shot is 2.8-5.6 times the 1e9 reference (`reference_fraction`). Pinned by `test_2033_hard_ops_product` and `test_2033_qsvt_matches_ch09_construction`; `test_qsvt_rounds_follow_the_success_probability` checks that a smaller success probability raises the cost through k. The breakdown (low end, 15 primitives, as printed by `--chapter ch10`):

```
magnetic_U_inv               1.469e+05 x      119 =  1.748e+07  COMPILED
magnetic_U_mul               2.938e+05 x      308 =  9.048e+07  COMPILED
magnetic_U_Tr                2.448e+04 x    574.6 =  1.407e+07  COMPILED
electric_U_FFT               4.896e+04 x     3397 =  1.663e+08  COMPILED
electric_U_phi               2.448e+04 x     7191 =   1.76e+08  COMPILED
hop_colour_squish_and_parity 2.448e+04 x 1.607e+04 =  3.934e+08  SCALING
hop_colour_rotations         2.448e+04 x 5.528e+04 =  1.353e+09  SCALING
hop_squish_and_flags         2.448e+04 x     2044 =  5.004e+07  COMPILED
hop_diagonalizers            2.448e+04 x 2.058e+04 =  5.039e+08  COMPILED
hop_phasing_hwp              2.448e+04 x    252.4 =   6.18e+06  SCALING
mass_hwp_rotations           3.264e+04 x    28.09 =  9.168e+05  SCALING
mass_hwp_toffolis            5.712e+04 x        7 =  3.998e+05  COMPILED
qsvt_reference_rotations          4488 x    28.09 =  1.261e+05  SCALING
qsvt_signal_rotations              535 x    28.09 =  1.503e+04  SCALING
qsvt_reflection_toffolis          8112 x        7 =  5.678e+04  SCALING
```

### 4.6 Shots

**2028.** 2e3 shots per sector are Stated. With H normalized by the total norm lambda of its measured terms and a per-shot variance sigma^2 <= 1 (Assumed), that is sqrt(1/2000) = 2.24% of lambda per sector and sqrt(2) x 2.24% = 3.16% of lambda on mu_c (`eps_per_sector`, `eps_mu_c_over_lambda`; the per-sector value is pinned by `test_r25_2028_exports_and_sectors`).

**2033 campaign** (`ratio_variance_per_sample`, lines 817-849; `grid_shots_direct`, lines 856-864). The estimator of chi_4/chi_2 is the ratio of k-statistics r = k_4/k_2 of the N_B samples. Each shot is an independent preparation, so the variance is V_1/N with V_1 the delta-method (influence-function) variance per sample; shots per grid point are V_1/delta^2 at delta = 0.1 absolute. The true N_B distribution at V = 2^3 is unknown, so V_1 is evaluated on a toy: a Skellam distribution in baryon number (`nb_toy_logpmf`, lines 804-814) with kappa_2(N_B) = 0.01 at mu_B = 0, the same at every T. `test_f1_delta_method_variance` checks the closed form V_1 = 36 + 72 kappa + 24 kappa^2 at mu_B = 0 and the quark-toy ratio 1/9.

| quantity | value | key | test |
|---|---|---|---|
| shots over the 25 points, kappa_2(0) = 0.01 | 257,603 | `shots_total` | `test_2033_shots_from_eq_nshot` |
| per point, min (mu_B = 0) to max (T = 100, mu_B = 600 MeV) | 3,672 - 60,535 | `shots_per_point_min`, `_max` | same |
| same grid at kappa_2(0) = 0.05 / 0.1 | 907,698 / 2,240,611 | `grid_shots_k2_mid`, `_hi` | `test_f1_kappa2_records` |
| quark-Skellam toy (kappa_2 = 0.05), 0.1 absolute / 10% relative | 4,787 / 387,747 | `grid_shots_quark_abs`, `_rel` | same |
| reweighting from two ensembles per T row | saves 11% at 0.01, costs 42% more at 0.05 | `reweighting_saving_k2_lo`, `reweighting_ratio_k2_mid` | same |

**2033 first result** (`row_samples`, lines 867-878). The T = 150 MeV row at 30% relative on chi_4/chi_2 at all five mu_B and on n_B for mu_B >= 150 MeV. For each of four toys (baryon Skellam at kappa_2 = 0.01, 0.05, 0.1; quark Skellam at 0.05) the model takes the smaller of direct sampling and reweighting from the best pair of ensembles (`fr_rw_pairs`). The range over the toys is 2,926 - 15,861 samples (`fr_samples`; `test_r25_2033_first_result`). Delta p is taken from the same samples.

### 4.7 Wall time, depth and factories (`_wall_shot_at_baseline`, lines 683-687; `common.depth_exports`)

One machine, serial: wall per shot = max(N_T x 1 us, D_T x 10 us) + 0.1 ms, where D_T is the T-depth and 10 us the reaction time per sequential non-Clifford layer. The first term wins whenever F* = N_T / D_T >= 10, the number of factories behind the 1 us per T convention.

- **2028.** D_T = 2.67e4 - 4.05e4 (Assumed, from a factory schedule of the step; `t_depth_2028`), F* = 9.88 - 14.97 (`f_star`). At the as-compiled end F* < 10, so the shot is depth-limited: 0.40499 s instead of 0.40006 s serial. 2e3 shots give 810.0 s = 13.5 min per sector (`wall_time_s`) and 1,620 - 2,430 s = 27 - 40 min for 2-3 sectors (`wall_campaign_s`). Tests: `test_2028_shots_wall_and_eps_l`, `test_r25_2028_exports_and_sectors`.
- **2033.** The factory band was computed for the 2e3-step circuit; the QSVT shot repeats the same step, so each end's depth is that band scaled by its own T ratio (0.502 and 1.020, `depth_scale`), giving F* = 17.5 - 56.1 (`f_star_matched`). The shot is serial at both ends: 2,772.6 - 5,634.1 s, 46 - 94 min (`wall_per_shot_s`). The contract export `f_star` pairs the bands crosswise (8.6 - 114) and its low corner is not a real circuit; the walls use the matched ends. Tests: `test_2033_wall_time_chain`, `test_r25_2033_exports`.

| 2033 quantity | value | key |
|---|---|---|
| first result (2,926 shots at the low end, 15,861 at the high) | 0.26 - 2.83 yr | `wall_first_result_yr` |
| campaign, 257,603 shots | 22.6 - 46.0 yr | `wall_campaign_yr` |
| one average grid point | 0.91 - 1.84 yr | `wall_single_point_yr` |
| depth floor of the campaign (unlimited factories) | 4.0 - 26.3 yr | `floor_wall_campaign_yr` |
| cut needed to fit the 5-year horizon | 4.5 - 9.2x | `horizon_reduction_needed` |
| per-shot T that would fit 5 years at the same shots | 6.1e8 | `t_per_shot_to_fit_horizon` |
| campaign at kappa_2(0) = 0.05 / 0.1 | 80 - 162 / 197 - 400 yr | `grid_wall_k2_mid_yr`, `grid_wall_k2_hi_yr` |

### 4.8 Required logical error rate

eps_l = 0.1 / N_T, so that a shot accrues 0.1 expected faults. 2028: 0.1 / 399,957 = 2.50e-7, so the RFI rate of 1e-8 suffices. 2033: 1.77e-11 at the upper end of the shot and 3.61e-11 at the lower, 277 - 563 times below the RFI rate; at 1e-8 a shot would accrue 28 - 56 faults (`eps_l_required`, `rfi_over_eps_l`, `faults_at_rfi_floor`; `test_2033_logical_error_row`).

### 4.9 Comparison routes the model also computes

These are reported in the chapter text but are not the headline.

- **The same filter at V = 2^2** (beta||H|| = 100, Assumed): 4.72e8 - 9.57e8 T, 4.7e3 - 9.6e3 times the 2028 cap (`qsvt_2x2_t`). This is why the 2028 benchmark is a ground-state ramp and not a thermal state.
- **V = 3^3 at 2033 parameters**: beta||H|| scaled with the volume to 337.5, d_beta = 56, a 9.56e6-T step, 1.82e10 - 3.70e10 T per shot on 1,141 - 1,641 LQ (`v3cubed_*`; `test_2033_v3cubed_promotion`).
- **A Lindbladian Gibbs sampler** (Chen-Kastoryano-Gilyen, arXiv:2311.09207; `gibbs_sample`, lines 576-619) on the same Trotter step: 288 local jumps, 20 sweeps, 8.39e13 T per sample at the central structure and 4.42e14 with the x K segment factor of Thm. I.2 (`gibbs_t_per_shot_range`, which the chapter prints as ~8e13-4e14), 1.5e4 - 1.6e5 times the QSVT shot (`gibbs_over_qsvt`; `test_f2_sampler_step_derivation`). Over the sweep band of 10-40 the central structure spans 7.6e12 - 9.5e14 (`gibbs_t_sample_band`). The same sampler at V = 2^2 is 1.3e12 T (`gibbs_2x2_t_sample`; `test_f2_2028_gibbs_sample_priced`).
- **2e3 plain Trotter steps** (c_mix x beta||H|| = 20 x 100) at 2.76e6 T each: 5.52e9 T (`t_hsim_lower_bound`). This circuit also sets the reference for the 2033 depth scaling.

## 5. What the paper prints

The rows below are what the chapter writes to `resources.json`, which feeds the paper's resource-landscape figure. The two benchmark boxes and the summary table (Table 1.1) print rounded values of the same numbers, recorded in `PUBLISHED`. The current rows:

| era | label | `lq` | `t` | flags |
|---|---|---|---|---|
| 2028 | 2+1D SU(3) mu_B > 0 | [176, 176] | [399956.82855117304, 399956.82855117304] | |
| 2033 | 3D SU(3) Sigma(36x3) 2^3, QSVT/TPQ | [514, 1014] | [2772642998.6754518, 5634084702.228161] | `conditional: true` |
| codesign | post-2033 Sigma(72x3) 4^3 | [2554, 3054] | [0.0, 0.0] | `lq_only: true`, `codesign: true` |

The `PUBLISHED` table (lines 1370-1386) holds what the boxes print: ~180 LQ and >~4.0e5 T (2028), ~500-1000 LQ and 2.8-5.6e9 T (2033), 2554-3054 LQ (stretch). `test_common.py::test_model_reproduces_published` requires the model to land within 10% of each, and `test_published_matches_the_boxes_as_printed` pins the printed values themselves. The other rows of the boxes, and where they come from:

| box row | 2028 | 2033 | tests |
|---|---|---|---|
| required eps_l | <~2.5e-7 | <~1.8e-11 (3.6e-11 at the lower end), 280-560x below RFI | `test_2028_shots_wall_and_eps_l`, `test_2033_logical_error_row` |
| shots | ~2e3 per sector, 2.2% of lambda | first result 2.9e3-1.6e4; campaign 2.6e5, 3.7e3-6.1e4 per point | `test_r25_2028_exports_and_sectors`, `test_r25_2033_first_result`, `test_2033_shots_from_eq_nshot` |
| per-shot wall | ~0.40 s (depth-limited) | 2.8-5.6e3 s (46-94 min) | `test_2028_shots_wall_and_eps_l`, `test_2033_wall_time_chain` |
| wall | ~13.5 min per sector, 27-40 min for 2-3 sectors | first result 0.26-2.8 yr; campaign 23-46 yr (0.91-1.8 yr per point) | `test_r25_2028_exports_and_sectors`, `test_r25_2033_first_result`, `test_2033_wall_time_chain` |

The intermediate numbers the chapter prose quotes (per-link costs, hop pieces, d_beta, call counts, factory counts, the kappa_2 sensitivity) are each pinned by the test named beside them in section 4. When the paper source sits next to the package, `test_r25_tex_prints_the_new_numbers` and `test_e29_tex_wording_fixes` also check that the chapter text contains them.

## 6. Circuit status

From `docs/CIRCUIT_STATUS.md` (generated by `python -m estimates --status`). COMPILED means an explicit gate-level circuit exists for the primitive; SCALING means the count rests on a counting argument. No Ch. 10 primitive is CONJECTURE or UNSOURCED.

| primitive | status | why |
|---|---|---|
| magnetic_U_inv, magnetic_U_mul, magnetic_U_Tr | COMPILED | arXiv:2405.05973 gate tables |
| electric_U_FFT | COMPILED | arXiv:2408.00075 |
| electric_U_phi | COMPILED | arXiv:2405.05973 |
| hop_squish_and_flags, hop_diagonalizers | COMPILED | gate tables of the fermion-primitives draft (unpublished) |
| mass_hwp_toffolis | COMPILED | Hamming-weight-phasing adders, arXiv:1709.06648, arXiv:1902.10673 |
| hop_colour_squish_and_parity | SCALING | the SU(3) squish uncompute is estimated here; the draft counts one squish |
| hop_colour_rotations | SCALING | the V^dag half of the frame is drawn in the draft but not counted there; added here |
| hop_phasing_hwp, mass_hwp_rotations | SCALING | grouping of equal-angle rotations into one HWP group, derived here |
| qsvt_reference_rotations, qsvt_signal_rotations, qsvt_reflection_toffolis (2033 only) | SCALING | the QSVT wrapper, priced by counting |

Across both eras the file lists 16 COMPILED and 11 SCALING rows. By T, SCALING primitives carry 59.6% of the 2028 step and 63.3% of the 2033 low-end shot, almost all of it the hop's color frame.

## 7. Open items and limitations

From `docs/OPEN_ITEMS.md` and the chapter's own gap list.

- **The 2033 price is conditional on the reference state.** It assumes a reference whose filtered amplitude is about 0.1 with thermal weights correctly spread over the N_B and field-number sectors. A Haar-random reference succeeds with probability Z e^{beta E_0}/D (up to sqrt(2^264) amplification iterations); a product reference is biased, with amplitude falling exponentially in V. No candidate closes this.
- **The QSVT degree and query cost are not derived.** d_beta is set by ||H||, not a block-encoding normalization; at the Chapter 9 value beta lambda = 864 (lambda/||H|| = 8.64) it would be round(sqrt(864 ln 1e4)) = 89 (`ch09_chiral_gauge.py`, `d_beta_stated`), and the shot about 3 times larger. The 2-4 Trotter steps per query are Stated. The 17 calls are the standard amplification count at amplitude 0.1; a fixed-point sequence would need more calls, which is not priced (`docs/OPEN_ITEMS.md`).
- **kappa_2 of N_B at V = 2^3 sets the campaign.** The grid needs 2.6e5 shots at 0.01, 9.1e5 at 0.05 and 2.2e6 at 0.1. A classical exact-diagonalization or strong-coupling estimate of kappa_2 and its T dependence would settle the band.
- **The chi_4/chi_2 target may be the wrong one.** It is an absolute 0.1; at small kappa_2 the ratio is pinned near 1/9 or 1 by quark or baryon content, so the precision on the deviation from that baseline may be what matters.
- **The 2028 count is a lower bound.** One ramp step is priced; the executed circuit is the N-step ramp plus the grouped energy readout (a ten-step ramp is 4.2e6 T). lambda, mu_c/lambda, the readout grouping, the ramp length and the term variances behind the 2e3 shots are open. The step itself is 4.0x the cap and needs a substantially cheaper color frame for the hop before it can run.
- **The 5% pressure target is not separately costed.** It references the classical mu_B = 0 pressure of this exact Hamiltonian; if that is not available, p needs an energy channel that is not priced.
- **The hop comes from an unpublished draft.** The frame undo, the SU(3) squish uncompute, the per-field color rotations, the C^nX convention (n - 1 Toffolis) and 7 T per flag Toffoli are this package's readings where the draft is silent; the recompute-per-application reading would raise the 2028 step to 7.9x the cap.
- **2033 depth is scaled, not scheduled.** The T-depth is the 2e3-step factory band scaled by each end's T ratio (F* 17.5-56, Assumed).
- **Ancillas.** The 2033 band of 250-750 is Stated, not sized. The 2028 itemization does not list the synthesis ancilla as its own line, and the draft's separate flag and parity registers are not sized.
- **Physics systematics outside the cost model.** Sigma(36x3) truncation (5-10%, a conjecture), O(a^2) discretization at one spacing (~20-30%), uncontrolled finite volume (m_pi L <~ 1), and the 4+2 taste content, which shifts the universality class and the critical point. The digitized Sigma(36x3) path integral has its own sign problem at mu_B = 0, so the classical-comparator verdict is not established.
- **Fault budget.** 0.1 expected faults per shot counts T gates only; Clifford, idle, measurement and injection faults are excluded.
- **Stretch.** The Sigma(72x3) instance has a register count only; no Sigma(72x3) Fourier transform circuit exists (arXiv:2511.17437 gives a lower bound only).

## 8. Conventions used

This chapter follows the package-wide conventions described in the top-level `README.md`: logical qubits count the algorithmic register only (magic-state factories excluded); a Toffoli costs 7 T; every synthesized rotation costs 1.15 log2(1/eps) + 9.2 T, with eps set per circuit from its own rotation count so the total synthesis error per shot is 0.01; the fault budget is 0.1 expected faults per shot; and wall times are for one machine, 1 us per T plus 0.1 ms per shot, with the depth correction of section 4.7 when fewer than ten factories can be kept busy.
