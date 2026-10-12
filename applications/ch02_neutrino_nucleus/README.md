# Chapter 2: neutrino-nucleus response for DUNE

Model: `estimates/ch02_nu_nucleus.py` (benchmarks) and `estimates/ch02_legacy.py` (the longer-term 12C scaling row).
Tests: `estimates/tests/test_ch02_nu_nucleus.py`, plus the per-chapter contract tests in `estimates/tests/test_common.py`.
Chapter text: `applications/app01_neutrino_nucleus.tex` in the report.

## 1. What is estimated

The observable is the real-time charge correlator
C(τ, q) = ⟨Ψ0| J†(q) e^{-i(H-E0)τ} J(q) |Ψ0⟩, measured by Hadamard tests at three specified times, real and imaginary parts.
A sampled time series of C determines a smeared inclusive response; the two costed deliverables are finite-volume model
benchmarks that can be checked classically, not physical cross sections.

| era | instance | lattice / model | register (LQ = logical qubits) |
|---|---|---|---|
| 2028 | four-nucleon charge correlator | one nucleon of each spin-isospin species, periodic 3^3 lattice, a = 1.4 fm, SU(4)-symmetric hopping t and contact couplings C0, D0 of arXiv:1911.06368; τ = 0.01, 0.02, 0.03 MeV^-1 (2, 4, 6 steps) | 36 LQ: 24 position, 11 clean workspace, 1 Hadamard |
| 2033 | fully mobile three-nucleon (triton) charge response | one proton and two neutrons, periodic 6^3 lattice, singlet/triplet two-body contacts plus a three-body contact fitted in the 6^3 box; τ = 0.05, 0.10, 0.20 MeV^-1 (25, 50, 100 steps of Δτ = 0.002) | 908 LQ: 864 occupation modes plus 44 workspace/measurement |
| codesign | converged-basis 12C evolution scaling | second-quantized chiral NN+3N in N_orb = 150-200 orbitals, N_orb^4/2 rotations per step, 200 steps | 800-1000 LQ |

Each benchmark headline is the deepest complete executed circuit: fresh deterministic preparation, current insertion,
controlled evolution and measurement. The 12C row is an evolution-only scaling comparison, not a complete response shot.
The register descriptions are from the chapter text; the model carries the totals (36, 908, and 4 N_orb + 200).

## 2. How to run

From the repository root:

```
python -m estimates --chapter ch02                  # assumptions, resources.json rows, per-era exports
python -m estimates --chapter ch02 --intermediates  # also every intermediate the model records
pytest estimates/tests/test_ch02_nu_nucleus.py -q
python applications/ch02_neutrino_nucleus/dump_assumptions.py   # the tables in section 3
```

`--chapter ch02` prints three blocks. ASSUMPTIONS lists all 142 fields of `Assumptions` with value, provenance tag and
source. INSTANCES prints the three rows this chapter writes to `resources.json`. EXPORTS prints, per era, the model's
register and per-shot T against the value printed in the chapter's benchmark box and its tolerance (`[ok]` when they
agree). It then prints the shot count, wall time, ε_l (the per-operation hard-fault rate, section 4), the wall time (a serial T-count bound since the 2026-10-12 ruling) and the depth/factory keys (`-` for the two benchmarks; see section 4), and the primitive breakdown whose sum is the
headline. Nothing is written to disk.

The test file has eight tests, and all run here: the study data they compare against (`neutrino_feasibility/`,
`2033_model_study/`) is mirrored in this folder since 2026-10-12, so the tests read it from either the report tree or here. The register and per-shot T of each era are also pinned by `test_model_reproduces_published[2]` in
`test_common.py` (section 5).

## 3. Assumptions

`Assumptions` in `ch02_nu_nucleus.py` (lines 31-34) adds one field, `coherent_synthesis_allowance` (study convention only), to the 141 fields it
inherits from `ch02_legacy.Assumptions` (lines 190-441). The table below is the output of `dump_assumptions.py`, pasted
unedited. It groups the fields by which code path reads them, measured by running `model(a, era)` for each era on an object
that records attribute reads. Values are the Python value of each field; source strings are the code's own, verbatim; the
code's notes are printed for group A and are otherwise in `ch02_legacy.py`. Six Assumed fields cite an internal
decision record rather than a publication (`R1`, `R2`, `R3`, `E27`, `TRACKED_CHANGES.md`); read them as working
assumptions. Provenance: Cited (key) means a published source gives the value; Assumed is a working assumption; Stated-not-derived
means the chapter asserts the number without derivation.

Two things to know before reading it:

- Only two fields enter the 2028 and 2033 numbers. Everything else those benchmarks use (step counts, gate counts per
  step, the register sizes, the shot tolerances) is written as literals in `counts` and `model`; they are listed after the
  generated table.
- Group C describes 4He, 12C and 40Ar instances that the chapter does not cost. Nothing exported reads these
  fields; they remain because the 12C scaling row shares the dataclass. Their `app01:NNN` line references point at the
  chapter text the fields were written against and do not line up with the current `app01_neutrino_nucleus.tex`.

<!-- begin dump_assumptions.py output -->

142 fields: Cited 43, Assumed 21, Stated-not-derived 78

#### A. Read by the 2028 and 2033 benchmarks (`model`, ch02_nu_nucleus.py): 5 fields

| name | value | provenance | source | note |
|---|---|---|---|---|
| `eps_syn` | `0.01` | Assumed | TRACKED_CHANGES.md R-TOL (H. Lamm, 2026-09-29) | ruling R-TOL: total synthesis error per shot 1e-2, report-wide; eps_rot = sqrt(eps_syn/N_rot) per circuit (randomized synthesis) |
| `t_gate_s` | `1e-06` | Assumed | app01:2028,2033 boxes | 1 us per T-gate (hard operation), the report-wide wall-time convention of ruling R8 (H. Lamm, 2026-09-29; convention in Ch. 1); was 10e-6. R9 (2026-10-02) renamed it from logical_cycle_s. app01 2028 box '~0.10 s'; 2033 box '11-16 min' |
| `shot_overhead_s` | `0.0001` | Assumed | E27; ch08_baryogenesis.py | per-shot overhead t0 ~ 0.1 ms (register initialization, final readout, decode), the report rule stated in Ch. 1; as Ch. 8 shot_overhead_s. 0.1% of the 2028 shot, < 2e-7 of a 2033 shot: no printed number moves (E27) |
| `fault_budget` | `0.1` | Stated-not-derived | app01:100,137 | 'expected number of logical faults per circuit to <~ 0.1'; '<= 0.1-expected-fault budget per shot' |
| `coherent_synthesis_allowance` | `0.005` | Assumed | app01:eq:nu_synthesis (pre-2026-10-12) | The chapter authors' study convention (sensitivity only since the 2026-10-12 ruling): coherent leading Ross-Selinger synthesis, delta_syn = 0.005 shared among N_R rotations; reproduces the study files' T_RS_leading_estimate / leading_T |

#### B. Read by the 12C co-design scaling row (`_model_codesign`, ch02_legacy.py): 33 fields

| name | value | provenance | source |
|---|---|---|---|
| `pauli_fraction` | `0.5` | Cited (arxiv_2312_05344) |  |
| `eps_target` | `0.05` | Stated-not-derived | app01:57,289 |
| `n_comp` | `5` | Stated-not-derived | app01:90 |
| `rfi_2033_t` | `1000000000.0` | Cited (DOE_RFI_2026) |  |
| `eps_l_2028_floor` | `1e-08` | Cited (DOE_RFI_2026) |  |
| `fault_budget_unit` | `1.0` | Assumed | R3 |
| `sigma_dune_MeV` | `10` | Stated-not-derived | app01:90,145 |
| `n_orb_explore` | `20` | Stated-not-derived | app01:166 |
| `n_orb_conv` | `(150, 200)` | Stated-not-derived | app01:65,166 |
| `n_orb_16o` | `(200, 300)` | Stated-not-derived | app01:94 |
| `n_orb_ar_route_i` | `(250, 500)` | Stated-not-derived | app01:96,119 |
| `t_ar_route_i` | `(10000000000000.0, 500000000000000.0)` | Stated-not-derived | app01:119 |
| `tau_max_GeVinv` | `100` | Stated-not-derived | app01:90,166 |
| `n_trotter_codesign` | `200` | Stated-not-derived | app01:166 |
| `window_trunc_factor` | `2` | Stated-not-derived | app01:90,168 |
| `n_anc` | `200` | Stated-not-derived | app01:94 |
| `n_anc_bracket` | `(100, 300)` | Stated-not-derived | app01:90 |
| `A_c` | `12` | Stated-not-derived | app01:58 |
| `shots_codesign` | `(300000.0, 600000.0)` | Stated-not-derived | app01:172,253 |
| `q0_max_GeV` | `1.0` | Stated-not-derived | app01:166,172 |
| `horizon_yr` | `5` | Stated-not-derived | app01:139,185,273 |
| `codesign_serial_yr` | `(17000.0, 110000.0)` | Stated-not-derived | app01:185 |
| `codesign_shots_in_horizon` | `(28, 93)` | Stated-not-derived | app01:185,249 |
| `codesign_reduction` | `(3000.0, 20000.0)` | Stated-not-derived | app01:185,249 |
| `qubitized_toffoli` | `(400000000000000.0, 7000000000000000.0)` | Cited (arxiv_2607_21563) |  |
| `qubitized_precision_MeV` | `0.01` | Cited (arxiv_2607_21563) |  |
| `rescale_bin_MeV` | `30` | Stated-not-derived | app01:171 |
| `eps_l_codesign` | `(2e-14, 6e-14)` | Stated-not-derived | app01:185 |
| `eps_l_codesign_shorthand` | `1e-13` | Stated-not-derived | app01:185 |
| `configs_momentum_transfers` | `(5, 10)` | Stated-not-derived | app01:277 |
| `n_nuclei` | `3` | Stated-not-derived | app01:277 |
| `configs_stated` | `100.0` | Stated-not-derived | app01:277 |
| `shots_per_config_stated` | `13000.0` | Stated-not-derived | app01:280 |

#### C. Inherited from ch02_legacy.Assumptions, read by no exported number: 104 fields

| name | value | provenance | source |
|---|---|---|---|
| `n_orb_2028` | `10` | Stated-not-derived | app01:125,196 |
| `n_trotter_2028` | `1` | Stated-not-derived | app01:126 |
| `pionless_gain` | `5` | Cited (arxiv_2312_05344) |  |
| `n_trotter_pionless` | `3` | Stated-not-derived | app01:126,207 |
| `anc_qpe_window` | `80` | Stated-not-derived | app01:203 |
| `anc_lcu_be` | `50` | Stated-not-derived | app01:204 |
| `anc_hadamard` | `(1, 2)` | Stated-not-derived | app01:204 |
| `lq_2028_prose` | `170` | Stated-not-derived | app01:94 |
| `n_bin_2028` | `1` | Stated-not-derived | app01:208 |
| `alpha_2028` | `10` | Assumed | shot_audit Ch. 2 |
| `jpsi_norm_sq_2028` | `1.5` | Assumed | shot_audit Ch. 2 |
| `n_comp_2028_full` | `2` | Assumed | shot_audit Ch. 2 |
| `n_parts_2028` | `2` | Assumed | shot_audit Ch. 2 |
| `eps_first` | `0.3` | Assumed | R1 |
| `strings_per_excitation_2028` | `8` | Assumed | factory.json |
| `sigma_2028_MeV` | `50` | Stated-not-derived | app01:128,197 |
| `rfi_2028_t` | `100000.0` | Cited (DOE_RFI_2026) |  |
| `rfi_2028_lq` | `(150, 250)` | Cited (DOE_RFI_2026) |  |
| `rfi_2033_lq` | `1000` | Cited (DOE_RFI_2026) |  |
| `eps_l_2028` | `1e-06` | Stated-not-derived | app01:2028 box |
| `t_4he` | `(670000000.0, 960000000.0)` | Stated-not-derived | app01:137,224 |
| `lq_4he` | `141` | Cited (arxiv_1911_06368) |  |
| `lattice_4he` | `3` | Stated-not-derived | app01:137,221 |
| `box_4he_fm` | `4.2` | Stated-not-derived | app01:137 |
| `hbar_c_MeV_fm` | `197.327` | Cited (PDG) |  |
| `q_coverage_3cube_GeV` | `(0.3, 0.51)` | Stated-not-derived | app01:246 |
| `q_coverage_4cube_GeV` | `(0.22, 0.77)` | Stated-not-derived | app01:246 |
| `w_q` | `12` | Stated-not-derived | app01:137,222 |
| `hop_t_MeV` | `10.5794` | Cited (arxiv_1911_06368) |  |
| `c0_abs_MeV` | `98.2265511` | Cited (arxiv_1911_06368) |  |
| `d0_MeV` | `127.839693` | Cited (arxiv_1911_06368) |  |
| `lambda_GeV` | `19.1` | Cited (arxiv_1911_06368) |  |
| `bin_4he_MeV` | `29.3` | Stated-not-derived | app01:144,223 |
| `bin_target_MeV` | `5` | Stated-not-derived | app01:145 |
| `eps_l_4he` | `1e-10` | Stated-not-derived | app01:137,226 |
| `block_frac` | `(0.68, 0.7)` | Stated-not-derived | app01:137 |
| `gs_prep_frac_block` | `0.015` | Stated-not-derived | app01:139 |
| `gs_prep_frac_ci` | `(0.156, 0.285)` | Stated-not-derived | app01:139 |
| `walk_frac` | `0.35` | Stated-not-derived | app01:141 |
| `walk_terms_per_site` | `38` | Cited (arxiv_1911_06368) |  |
| `walk_rot_per_u2` | `1` | Assumed | arxiv_1911_06368 |
| `p_bin_min` | `0.1` | Assumed | R2 |
| `n_comp_first` | `1` | Assumed | R1 |
| `eps_norm_campaign` | `0.02` | Assumed | shot_audit Ch. 2 |
| `eps_norm_first` | `0.1` | Assumed | r25 |
| `norm_grouped_factor` | `(1, 3)` | Assumed | shot_audit Ch. 2 |
| `norm_snapshot_var` | `0.6666666666666666` | Assumed | shot_audit Ch. 2 |
| `prepare_ancilla_4he` | `21` | Assumed | factory.json |
| `select_depth_per_index` | `(1, 3)` | Assumed | factory.json |
| `plus_bit_factor` | `(1.34, 1.38)` | Stated-not-derived | app01:141 |
| `block_plus_bit_frac` | `(0.91, 0.97)` | Stated-not-derived | app01:141 |
| `plus_cube_factor` | `(5.4, 5.6)` | Stated-not-derived | app01:141 |
| `tau_grid_frac` | `(0.63, 0.65)` | Stated-not-derived | app01:142 |
| `tau_grid_lq` | `130` | Cited (arxiv_1911_06368) |  |
| `modes_per_site` | `4` | Assumed | app01:148 |
| `spin_isospin_bits` | `2` | Cited (arxiv_2507_22814) |  |
| `lattice_12c` | `4` | Stated-not-derived | app01:150,233 |
| `box_12c_fm` | `5.6` | Stated-not-derived | app01:150 |
| `box_5cube_fm` | `7.0` | Stated-not-derived | app01:150 |
| `t_12c` | `(2700000000.0, 2800000000.0)` | Stated-not-derived | app01:152,234 |
| `lq_12c` | `292` | Cited (arxiv_1911_06368) |  |
| `bin_12c_MeV` | `34.7` | Stated-not-derived | app01:152,234 |
| `w_q_12c` | `13` | Stated-not-derived | app01:152,234 |
| `eps_l_12c` | `3.6e-11` | Stated-not-derived | app01:152,235 |
| `five_cube_factor` | `(10, 11)` | Stated-not-derived | app01:150 |
| `w_q_5cube` | `14` | Stated-not-derived | app01:150 |
| `te_pe_gain` | `8` | Cited (arxiv_1911_06368) |  |
| `t_onramp` | `(20000000.0, 50000000.0)` | Stated-not-derived | app01:154,239 |
| `lq_onramp` | `280` | Cited (arxiv_1911_06368) |  |
| `t_onramp_eps_ref` | `30000000.0` | Stated-not-derived | app01:154,240,288 |
| `eps_l_onramp` | `3e-09` | Stated-not-derived | app01:154,240,288 |
| `euclid_factor` | `(1.5, 3.1)` | Stated-not-derived | app01:154 |
| `A_ar` | `40` | Stated-not-derived | app01:44 |
| `lattice_ar` | `8` | Cited (arxiv_2507_22814) |  |
| `fq_dim` | `3` | Cited (arxiv_2507_22814) |  |
| `fq_eps` | `0.1` | Cited (arxiv_2507_22814) |  |
| `fq_kin_MeV` | `10.58` | Cited (arxiv_2507_22814) |  |
| `fq_c_abs_MeV` | `98.23` | Cited (arxiv_2507_22814) |  |
| `fq_g_MeV` | `127.84` | Cited (arxiv_2507_22814) |  |
| `fq_spacing_fm` | `1.4` | Cited (arxiv_2507_22814) |  |
| `fq_mass_MeV` | `939` | Cited (arxiv_2507_22814) |  |
| `fq_cross_E_MeV` | `10` | Cited (arxiv_2507_22814) |  |
| `fq_heavy_H_per_nucleon_MeV` | `18` | Cited (arxiv_2507_22814) |  |
| `fq_dw_MeV` | `100` | Cited (arxiv_2507_22814) |  |
| `fq_dw_extrap_MeV` | `10` | Stated-not-derived | app01:156,249 |
| `fq_t_src_per_toffoli` | `4` | Cited (arxiv_2507_22814) |  |
| `tab5_gqsp_t_dw` | `33800000.0` | Cited (arxiv_2507_22814) |  |
| `tab5_gqsp_q_dw` | `498` | Cited (arxiv_2507_22814) |  |
| `tab5_gqsp_t_cross` | `211000000.0` | Cited (arxiv_2507_22814) |  |
| `tab5_gqsp_q_cross` | `500` | Cited (arxiv_2507_22814) |  |
| `tab5_gqsp16_t_dw` | `5990000.0` | Cited (arxiv_2507_22814) |  |
| `tab5_gqsp16_t_cross` | `37500000.0` | Cited (arxiv_2507_22814) |  |
| `tab5_gqsp16_q` | `(234, 235)` | Cited (arxiv_2507_22814) |  |
| `tab5_trotter2q_t_dw` | `347000000.0` | Cited (arxiv_2507_22814) |  |
| `tab5_trotter2q_t_cross` | `5890000000.0` | Cited (arxiv_2507_22814) |  |
| `tab5_trotter2q_q` | `3072` | Cited (arxiv_2507_22814) |  |
| `lq_ar` | `(580, 601)` | Stated-not-derived | app01:156,159,243 |
| `lq_ar_requirements` | `(580, 601)` | Stated-not-derived | app01:286 |
| `ar_amplified_factor` | `(0.15, 9.3)` | Stated-not-derived | app01:158,244 |
| `ar_postselected_factor` | `1.6` | Stated-not-derived | app01:159,244 |
| `ar_repeats` | `1000.0` | Stated-not-derived | app01:159 |
| `eps_l_ar` | `(1.1e-11, 6.7e-10)` | Stated-not-derived | app01:159,245,288 |
| `eps_l_ar_postselected` | `6.3e-11` | Stated-not-derived | app01:159 |
| `ar_window_MeV` | `(10, 100)` | Stated-not-derived | app01:156,247 |

## 4. How the numbers are built

### Benchmarks (2028, 2033): `ch02_nu_nucleus.py`

**Gate counts.** `counts(era, n)` (lines 37-44) returns the number of arbitrary rotations N_R, Toffolis N_Toff and
other exact T gates N_T^fixed of the complete circuit with n evolution steps, from the literals above. For the headline
circuits:

| era | n | N_R | N_Toff | N_T^fixed |
|---|---|---|---|---|
| 2028 | 6 | 1,192 | 2,240 | 744 |
| 2033 | 100 | 7,884,133 | 60,066 | 0 |

**T per shot.** `circuit(a, era, n, convention='report')` (lines 47-66) prices each circuit with the chapter's Eq. (nu_synthesis)
in the report convention (ruling 2026-10-12):

  N_T = N_R [1.15 log2(1/ε_R) + 9.2] + 7 N_Toff + N_T^fixed,  ε_R = sqrt(ε_syn / N_R),  ε_syn = `eps_syn` = 0.01,

i.e. randomized synthesis, errors adding incoherently, at the Bocharov-Roetteler-Svore fit (`t_per_rotation`,
`eps_per_rotation` in `common.py`). The tolerance is set per circuit, so the deeper circuits pay slightly more per
rotation. `model` (lines 69-112) evaluates all three times and returns the deepest as the headline. Per-time costs
(`intermediates['costs_by_time']`):

| era | steps | T per rotation | N_T |
|---|---|---|---|
| 2028 | 2, 4, 6 | 18.59, 18.76, 18.90 | 26,812; 32,866; 38,948 |
| 2033 | 25, 50, 100 | 25.32, 25.72, 26.19 | 7.0441e+07; 1.1544e+08; 2.0694e+08 |

The headline breakdown (`Result.breakdown`, three `Primitive` rows that sum exactly to the headline; the test checks this):

| era | rotations | Toffoli | fixed T | total |
|---|---|---|---|---|
| 2028 | 1,192 × 18.896 = 22,524 | 2,240 × 7 = 15,680 | 744 | 38,948.32 |
| 2033 | 7,884,133 × 26.194 = 2.0652e+08 | 60,066 × 7 = 420,462 | 0 | 2.06936e+08 |

At 2033 the rotations are 99.8% of the T count. Of the 7,884,133 rotations, 6,835,104 (87%) are in the evolution term
58,752n + 9,504(n+1) and 1,048,596 (13%) in the loader.

Two sensitivities are exported next to the headline. `convention='study'` (`study_T_per_shot`, `costs_by_time_study`,
`study_T_per_rotation`) is the chapter authors' own pricing, the coherent leading Ross-Selinger term 3 log2(N_R/δ) at
δ = 0.005: 54,482; 67,321; 80,302 T (2028) and 2.4132e+08; 3.9930e+08; 7.2310e+08 T (2033), 2.06x and
3.49x the headline; these are the counts the study files under this folder record, and the chapter quotes them in
parentheses. `convention='coherent'` (`coherent_T_per_shot`) is the BRS fit with coherently adding errors, ε_R = ε_syn/N_R:
5.051e+04 and 3.409e+08 T. Neither enters a headline.

**Shots.** Eq. (Nshot) gives ⌈(2/η²) ln(2m/α)⌉ repetitions per bounded estimate, for m estimates at absolute
tolerance η and simultaneous failure probability α. With m = 6, η = 0.02 and α = 0.05 that is 27,404 per estimate
(`shots_per_estimate`) and 164,424 for the six (`campaign_shots`, `Result.shots`) in both eras. The 2033 first result (Re and Im
at τ = 0.10, η = 0.05) is 2 × 3,506 = 7,012 shots (`first_result_shots`).

**Campaign gate volume.** 2028: every shot is priced at the deepest circuit, 164,424 × 38,948 = 6.4040e+09 T
(`campaign_T`, scope "deepest-shot upper estimate"). 2033: each time gets 2 × 27,404 shots at its own cost,
2 × 27,404 × (7.0441e+07 + 1.1544e+08 + 2.0694e+08) = 2.1529e+13 T ("sum of time-specific costs"). First result:
7,012 × 1.1544e+08 = 8.0945e+11 T (`first_result_T`).

**Logical error rate.** `epsilon_l` = `fault_budget` / N_T = 0.1 / N_T, the rate at which the deepest shot expects 0.1
hard faults: 2.5675e-06 (2028) and 4.8324e-10 (2033). This is a fault allowance, not a total observable-error estimate.

**Wall time.** Serial T-count bounds since the 2026-10-12 ruling (line 86): the report's 1 µs per T plus 0.1 ms per shot on one
machine, with no T-depth, Clifford depth, routing or factory schedule derived for either benchmark. So `t_depth_per_shot`,
`f_star`, `floor_wall_s` and `factories_for_1yr` are `None`, `depth_status = "serial_t_count_bound"`, and
`wall_per_shot_s` = 0.0390 s (2028) / 206.9 s = 3.4 min (2033); `wall_campaign_s` (= `Result.wall_time_s`) =
6420 s = 1.8 h (2028) / 2.1529e+07 s = 249 d = 0.68 yr (2033); `wall_first_result_s` = 8.0945e+05 s = 9.4 d (2033 only). They are lower
bounds: a depth-bound circuit on the 44 workspace qubits could run several times longer. `test_common.py::test_model_exports_depth`
allows this status for Ch. 2 alone.

### 12C scaling row (codesign): `ch02_legacy._model_codesign`

`model(a, 'codesign')` hands off to `ch02_legacy.model` (line 1386), which calls `_model_codesign` (lines 1227-1383).
The cost chain uses `_route_i` (lines 448-457) and the report-wide helpers `eps_rot_for` and `t_per_rotation` in
`common.py`:

- Pauli strings per step = `pauli_fraction` × N_orb^4 = 0.5 N_orb^4: 2.53125e8 at N_orb = 150, 8e8 at N_orb = 200.
- Rotations per shot N_R = strings × `n_trotter_codesign` (200): 5.0625e10 and 1.6e11.
- Per-rotation tolerance, randomized synthesis: ε_R = sqrt(`eps_syn` / N_R) with `eps_syn` = 0.01: 4.444e-7 and 2.5e-7.
- T per rotation, Bocharov-Roetteler-Svore fit: 1.15 log2(1/ε_R) + 9.2 = 33.467 and 34.421.
- T per shot = N_R × T per rotation = 1.6943e12 and 5.5074e12 (one `rot_synth` primitive, status SCALING).
- Logical qubits = 4 N_orb + `n_anc` (200) = 800 and 1000.

The function also exports quantities the chapter does not print:

- shots 3.2e5-6.4e5 (ε_target^-2 × N_bin × N_comp with ε_target = 0.05, N_comp = 5, N_bin = (5-10) × 32 Nyquist τ
  points out to τ_max = 1/σ = 100 GeV^-1);
- per-shot wall time T × 1 µs + 0.1 ms = 1.694e6-5.507e6 s (19.6-63.7 days on one machine);
- ε_l = 0.1/T = 1.82e-14 to 5.90e-14;
- the remaining 75 intermediate keys (`--intermediates` lists them).

None of these enter `resources.json` beyond the register and T band.

Since the 2026-10-12 ruling the benchmarks and the 12C row use the same report-wide randomized convention; the chapter
quotes the study's coherent counts in parentheses next to each headline and says so next to Eq. (nu_synthesis).

## 5. What the paper prints

Rows in `resources.json` (`ch: 2`), with the values as stored:

| era | label | lq | t | codesign |
|---|---|---|---|---|
| 2028 | Four-nucleon charge correlator, 3^3, complete | [36, 36] | [38948.317136396356, 38948.317136396356] | false |
| 2033 | Mobile triton charge response, 6^3, complete | [908, 908] | [206935602.44923684, 206935602.44923684] | false |
| codesign | 12C converged-basis evolution scaling (incomplete circuit) | [800, 1000] | [1694252578823.0613, 5507408616755.647] | true |

How each printed number is backed:

| printed (chapter box / Table 1.1) | model value | pinned by |
|---|---|---|
| 2028: 36 LQ; 3.89e4 T | `lq` 36; `hard_ops` 38,948.32 | `test_common.py::test_model_reproduces_published[2]` against `PUBLISHED['2028']` = (36, 3.89e4) at rel_tol 0.005; `test_ruling_2026_10_12_report_synthesis_convention_and_serial_walls` (3.8948e4); in the report tree also `test_table.py` (box at its printed precision, both copies of Table 1.1) |
| 2028: 8.03e4 in the study's coherent synthesis convention | `study_T_per_shot` 80,302.18 | the ruling test (8.0302e4, 2.06x) and `test_2028_matches_independent_circuit_recount` against the study files |
| 2028: 27,404 shots per estimate; 164,424 shots | `shots_per_estimate`, `Result.shots` | `test_complete_counts_shots_and_no_invented_schedule` (exact) |
| 2028: at most 6.40e9 T; ≈ 1.8 h | `campaign_T` 6.4040e+09; `wall_campaign_s` 6420 s | the ruling test |
| 2028: ε_l ≲ 2.6e-6 | `epsilon_l` 2.5675e-06 | the ruling test |
| 2028 table: 2.68e4, 3.29e4, 3.89e4 T | `costs_by_time` | the ruling test; `test_2028_matches_independent_circuit_recount` (n = 2, 6) compares the counts and the study-convention T with the study's recount files |
| 2033: 908 LQ; 2.07e8 T | `lq` 908; `hard_ops` 2.0694e+08 | `test_model_reproduces_published[2]` against `PUBLISHED['2033']` = (908, 2.07e8) at rel_tol 0.005; the ruling test (2.0694e8); `test_table.py` in the report tree |
| 2033: 7.23e8 in the study's coherent synthesis convention; 3.4 min per shot | `study_T_per_shot` 7.2310e8; `wall_per_shot_s` 206.9 s | the ruling test (3.49x) and `test_2033_matches_construction_and_loader_checks` against the study files |
| 2033: first result 7,012 shots, 8.09e11 T, ≈ 9.4 d | `first_result_shots`, `first_result_T` 8.0945e+11, `wall_first_result_s` 8.095e+05 s | `test_complete_counts_shots_and_no_invented_schedule` and the ruling test |
| 2033: campaign 164,424 shots, 2.15e13 T, ≈ 250 d | `Result.shots`, `campaign_T` 2.1529e+13, `wall_campaign_s` 2.153e+07 s | same tests |
| 2033: ε_l ≲ 4.8e-10 | `epsilon_l` 4.8324e-10 | the ruling test |
| 2033 table: 7.04e7, 1.15e8, 2.07e8 T; 1,048,596 rotations; 60,066 Toffolis | `costs_by_time`, `counts` | the ruling test; `test_2033_matches_construction_and_loader_checks` (counts and study-convention T against the study files) |
| 12C: 800-1000 LQ; 1.7-5.5e12 T | `lq` (800, 1000); `hard_ops` (1.6943e12, 5.5074e12) | `test_model_reproduces_published[2]` against `PUBLISHED['codesign']` at rel_tol 0.1; the second-instance check in `test_table.py` |

Other tests. `test_fermionic_reference_and_observable_checks` checks the 2033 study's classical
cross-checks: energy and correlator agreement below 1e-9, second-order step convergence, and an 8^3 ground energy that
differs from the 6^3 one by more than 2 MeV, so the volume limitation stays visible.
`test_tighter_synthesis_budget_increases_count` checks that a tighter `eps_syn` raises the 2033 count and that
`coherent_synthesis_allowance` moves only the study-convention sensitivity. The schedule test checks that `t_depth_per_shot`
and `f_star` stay `None` and that `Result.wall_time_s` is the serial campaign bound.

## 6. Circuit status

From `docs/CIRCUIT_STATUS.md` (generated by `python -m estimates --status`):

| era | primitive | count | T each | T total | status | source |
|---|---|---|---|---|---|---|
| 2028 | rotations | 1.19e+03 | 18.9 | 2.25e+04 | COMPILED | editorial_review/neutrino_feasibility/LARGER_BOX.md |
| 2028 | Toffoli | 2.24e+03 | 7 | 1.57e+04 | COMPILED | same |
| 2028 | fixed_T | 744 | 1 | 744 | COMPILED | same |
| 2033 | rotations | 7.88e+06 | 26.2 | 2.07e+08 | SCALING | editorial_review/2033_model_study/TRITON.md |
| 2033 | Toffoli | 6.01e+04 | 7 | 4.2e+05 | SCALING | same |
| 2033 | fixed_T | 0 | 1 | 0 | SCALING | same |
| codesign | rot_synth | 5.06e+10 | 33.5 | 1.69e+12 | SCALING | arxiv_2312_05344; bocharovRoettelerSvore2015, campbell2017 |

COMPILED 3, SCALING 4, CONJECTURE 0, UNSOURCED 0. COMPILED for 2028 means explicit unsynthesized circuits and an
independent gate recount exist; the rotations are still priced by the synthesis fit, not synthesized. The
2033 counts are a constructive allowance for the 908-qubit circuit, which has not been exported as a full-register
circuit. The benchmark statuses are set in `model` (line 104); the 12C status is set in
`ch02_legacy._model_codesign`.

## 7. Open items and limitations

From `docs/OPEN_ITEMS.md`:

- No T-depth or device schedule is derived for either benchmark: the depth and factory exports are `None`
  (`depth_status = "serial_t_count_bound"`), and the printed walls (1.8 h; 9.4 d and 250 d) are serial T-count lower
  bounds that a depth-bound circuit on 44 workspace qubits could exceed several times.
- The 2033 primitive counts are constructive allowances, not a compiled full-register circuit.
- Rotation synthesis follows the report convention since 2026-10-12 (randomized BRS fit at ε_syn = 1e-2,
  `exact_synthesis = False`); the chapter authors' study files price the same rotations coherently at the leading
  Ross-Selinger term (`convention='study'`, 2.06x and 3.49x the headline), and coherent accumulation under the BRS fit
  would give 5.1e4 and 3.4e8 T (`coherent_T_per_shot`).
- Physical accuracy of the finite-volume correlators and the logical-noise behavior of the circuits are not derived.
- The 12C row is an evolution-scaling estimate on the legacy assumptions, not a complete response shot, with no
  step-error validation.

From the chapter text (physics limits of the benchmarks as defined):

- The 12C estimate's step accuracy, three-body implementation, preparation and current/readout costs need their own
  audit.

- 2028: the couplings are not refitted to the 4.2 fm box or to physical helium; q a = 2π/3 gives large dispersion
  artifacts and q ≈ 295 MeV is beyond controlled pionless EFT. Validation is against the finite-volume Hamiltonian, not
  (e,e') data.
- 2033: the three-body coupling is fitted to the 6^3 ground energy, and the same couplings give -15.12, -8.48 and -6.45 MeV
  on 4^3, 6^3 and 8^3, so the model is not volume converged; q ≈ 148 MeV is near the pionless breakdown scale; effective-range
  terms, tensor forces and two-body currents are omitted. The elastic pole is subtracted using a classical calculation whose
  quantum determination is not costed, and the inclusive correlator does not separate breakup channels.
- Neither benchmark establishes quantum advantage, a physical cross section, or a reduction of DUNE's systematic error.

Code-level points a reader should know:

- Only `eps_syn`, `fault_budget`, `t_gate_s`, `shot_overhead_s` and (study convention only) `coherent_synthesis_allowance` are adjustable inputs to the benchmarks. Step counts, gate counts
  and shot tolerances are literals, so sensitivity studies on them mean editing `counts` and `model`.
- `Assumptions` carries 104 fields (group C) that no exported number reads, with source line numbers that refer to a
  different version of the chapter text.
- The study data the tests compare against is mirrored in this folder (`neutrino_feasibility/`, `2033_model_study/`).

## 8. Conventions used

The report-wide conventions are described in the top-level `README.md`:

- the per-circuit rotation tolerance rule (ε_R = sqrt(10^-2 / N_R), randomized synthesis);
- 7 T per Toffoli;
- the logical-qubit count as the algorithmic register, excluding magic-state factories;
- the wall-time model (1 µs per T gate plus 0.1 ms per shot, serial on one machine).

Chapter 2 uses all four since the 2026-10-12 ruling (before it the two benchmarks priced rotations with coherent
Ross-Selinger synthesis at δ = 0.005 per circuit and quoted no wall time). The walls are the serial T-count form of the
wall-time model only, since no T-depth has been derived.
