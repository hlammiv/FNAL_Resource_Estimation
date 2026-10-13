# Chapter 8: baryogenesis (bubble nucleation and C-violating wall scattering)

Model: `estimates/ch08_baryogenesis.py`. Tests: `estimates/tests/test_ch08_baryogenesis.py`. Chapter source: `applications/app05_baryogenesis.tex` in the paper. Unless a line says otherwise, every number below is a model intermediate (printed by `python -m estimates --chapter ch08 --intermediates`), a value in `resources.json`, or an assertion in the test file. Each number names its source.

## 1. What is estimated

The chapter has two real-time calculations. The first is the false-vacuum decay rate Γ/V, measured as the slope of ln P_FV(t). The second is the reflection asymmetry of fermions off a bubble wall, A_R = |R_p|² − |R_p̄|², which builds the charge asymmetry Δn. The model prices two instances and no co-design instance. `model(a, "codesign")` raises, and `test_codesign_is_refused` checks that it does.

- **2028 benchmark: pure-scalar bubble nucleation on 6² (2+1)D.** One real φ⁴ scalar with a biased double well, N_φ = 16 field values per site (4 qubits). The false vacuum is prepared adiabatically and evolved in closed real time with 20 first-order Trotter steps of δt = 0.1/m_φ. P_FV is measured at two times, 10 and 20 steps (1.0 and 2.0/m_φ, both inside the chapter's 1-2/m_φ exponential window; ruling 2026-10-11, was 4 and 12 steps), and the slope gives Γ_vac/V to 20% at one parameter point. The same function also records three side instances. Two are 8² family rows kept at the r16 record's 12 steps (`n_trot_family`; N_φ = 4 and N_φ = 16). The third is an interim 4² scalar plus two Wilson flavors at N_φ = 4, flagged conditional in `resources.json`. The model also records the cost the fermion arm would have at N_φ = 16 on 6²; that arm waits for 2033.
- **2033 target: wall scattering on a 15×8 cylinder.** One scalar (N_φ = 16) plus two 2-component Wilson fermion flavors (r = 1), coupled by the C-violating Yukawa term |y|φ[e^{iθ_C(z)} ψ̄_L ψ_R + h.c.], at a m_φ = 0.5. The 15-site axis is the scattering axis and is open; the 8-site transverse axis is periodic. The instance has two arms.
  - *Nucleation:* a Lindblad bath on 30 shared ancillas prepares the thermal state. P_FV is measured at t/3 and t (33 and 100 steps) on a ~100-point (T, Δθ_C) scan.
  - *Δn readout:* one base circuit A = U_evo U_pk U_vac prepares the free vacuum of the wall background (U_vac), loads a flavor-L particle or antiparticle packet chosen by a control qubit (U_pk), and evolves for 100 steps (U_evo). The amplitude of a flag on the incident-side charge is then read by depth-capped maximum-likelihood amplitude estimation (MLAE).

  Two tiers are priced. The first result is a 3σ detection of Δn at T_c and 0.9 T_c, together with both nucleation slopes. The campaign is Δn to ±10⁻³ in 8 momentum bins (the distinct |k_z| on the 15-site axis; ruling 2026-10-11, was 100) at 10 points, plus the nucleation slope at all 100 scan points.

## 2. How to run

From the repository root:

```
python -m estimates --chapter ch08                  # assumptions, resources.json rows, per-era exports
python -m estimates --chapter ch08 --intermediates  # also every named intermediate (89 for 2028, 164 for 2033)
python -m pytest estimates/tests/test_ch08_baryogenesis.py -q
python applications/ch08_baryogenesis/dump_assumptions.py   # the table in section 3
```

`--chapter ch08` prints four things. First, the 83 assumptions with their provenance. Second, the five rows this chapter writes to `resources.json`. Third, for each era: the headline against the printed box (`[ok]` when it falls within the box's 5% tolerance), shots, wall time, ε_l, the T-depth and factory keys, and the primitive-by-primitive T breakdown, whose sum equals the headline. Fourth, the model's notes, which spell out the arithmetic.

`flag_contrast` diagonalizes 480×480 one-body Hamiltonians with numpy, so numpy is required. The whole chapter runs in under a second.

Test status in this repository: 48 passed, 4 skipped. The four skipped tests also read the chapter `.tex`, and they skip when the paper source is absent: `test_no_machine_count_anywhere`, `test_utility_box_has_no_unsourced_dollar_figure`, `test_g4_fault_bias_model` and `test_referee_text_hooks`. In `test_no_machine_count_anywhere` and `test_g4_fault_bias_model` the numeric assertions run before the file is read, so a wrong number still fails. Only the text check ends in a skip.

## 3. Assumptions

Output of `dump_assumptions.py`, generated from the live `Assumptions` dataclass (lines 243-460). The model has a single dataclass. The rows are grouped by its section comments, and each row states which instance it serves. The provenance tags are the code's own. **Cited** means a published source or report-wide ruling gives the value. **Assumed** means a working assumption the model declares. **Stated-not-derived** means a chapter value taken without a derivation; many of these are printed values the tests compare against the model. Source tags of the form `app05:NN` name the chapter; their line numbers do not match the current file, so find the passage by the quoted phrase. Some quoted phrases (for example in `c_t_2028_quoted`, `eps_l_2033_base_quoted`, `lq_two_register_quoted`, `shots_2028`, `t_gate_s`) are not in the current chapter text; those values are model inputs the chapter uses without printing them.

92 fields: Cited 9, Assumed 27, Stated-not-derived 56

**Lattice and registers**

| field | value | provenance | source | note |
|---|---|---|---|---|
| `v_2028` | 36 | Stated-not-derived | app05:49,122 | 6^2 (2+1)D lattice, 36 sites |
| `n_phi` | 16 | Stated-not-derived | app05:122,157 | N_phi=16 (4 qubits/site); resolves both minima and the barrier |
| `n_phi_small` | 4 | Stated-not-derived | app05:134,176 | N_phi=4 (2 qubits/site) in the family row and the interim variant |
| `n_f_2033` | 2 | Stated-not-derived | app05:156,163 | 2 Wilson-fermion flavors; 2N_f=4 fermion qubits/site |
| `spinor_components` | 2 | Stated-not-derived | app05:95 | 2^floor((D+1)/2)=2-component Dirac spinor in 2+1D |
| `n_anc_2028` | 30 | Stated-not-derived | app05:93,126-127 | ~30 ancilla: Trotter control, basin projector |
| `n_anc_2033` | 40 | Stated-not-derived | app05:95,164 | 40 ancilla (30 bath + 10 workspace) |
| `n_bath_anc_2033` | 30 | Stated-not-derived | app05:159,164 | Lindblad bath, 30 shared ancillas (reset each step) |
| `n_anc_generic` | (10, 40) | Stated-not-derived | app05:86 | 'N_anc ~ 10-40' in the scaling section |
| `lz_sites` | 15 | Stated-not-derived | app05:114,156 | 15-site scattering axis |
| `lperp_sites` | 8 | Stated-not-derived | app05:114,156 | 8-site transverse direction |
| `a_mphi` | 0.5 | Stated-not-derived | app05:157 | a m_phi = 0.5 |
| `v_family` | 64 | Stated-not-derived | app05:134 | 8^2 family rows |
| `v_interim` | 16 | Stated-not-derived | app05:176 | 4^2 interim scalar+Wilson-fermion variant |
| `workspace_interim_prior` | 30 | Stated-not-derived | app05:177 | '~30-qubit amplitude-estimation workspace' of a phase-register estimator, of which ~13 qubits are the phase register; MLAE needs no phase register |
| `workspace_interim` | 17 | Stated-not-derived | app05:177 | '~17-qubit amplitude-estimation workspace ... a workspace budget, not an itemized count' = 30 - 13 (MLAE needs no phase register) |
| `lq_interim_quoted` | 115 | Stated-not-derived | app05:177 | 'The total is ~115 LQ'; 96+17=113 |
| `lq_two_register_quoted` | 190 | Stated-not-derived | app05:177 | 'two-register protocol ... ~190': 96+96=192 |

**Inputs to the gate-level c_T**

| field | value | provenance | source | note |
|---|---|---|---|---|
| `lattice_dims` | 2 | Stated-not-derived | app05:49,156 | 2+1D: d = 2 hopping directions |
| `links_per_site` | 2 | Assumed | derived here | d = 2 links per site, as for periodic boundaries (arxiv_2407_13819 E_D); periodic on the square 2028 lattices; the 2033 15x8 cylinder is open along its 15-site axis, 8 of 240 links fewer (~1% of c_T), so the 2033 count is slightly conservative |
| `wilson_r` | 1 | Assumed | derived here | Wilson parameter r = 1: the hopping matrix (r beta - i alpha_j)/2 has rank 1 |
| `eps_syn` | 0.01 | Cited | TRACKED_CHANGES.md R-TOL (H. Lamm, 2026-09-29) | report-wide: total synthesis error per circuit 1e-2; eps_rot = sqrt(eps_syn/N_rot) |
| `rot_synthesis` | rus | Cited | bocharovRoettelerSvore2015,campbell2017 | report-wide: 1.15 log2(1/eps_rot) + 9.2 T per rotation |
| `toffoli_convention` | textbook | Cited | main-overview:47 (report-wide) | 7 T per Toffoli everywhere |
| `c_t_2028_quoted` | 1260 | Stated-not-derived | r12 record | c_T 1263.1 T/site/step at 20 steps (the rotations of a longer shot are synthesized tighter; 1239.8 at 12 steps); derived here, no longer printed since the Oct 6 cut |
| `c_t_2033_quoted` | 2500 | Stated-not-derived | app05:91,165 | 'c_T ~ 2.5e3 (1.6e3 scalar, 8.9e2 fermion)' at the deepest circuit's budget |
| `c_t_range_quoted` | (250, 2500) | Stated-not-derived | app05:85 | 'c_T runs from 2.5e2 to 2.5e3 across the instances below' |
| `ct_reduction` | (3, 4) | Assumed | - | model-internal sanity span for the banked ~3.3x (not printed in the chapter); assumed development target, not derived or cited |
| `ct_reduction_banked` | 3.3 | Assumed | app05:91,126,162,163,179 | assumed development target, not derived or cited: 'divided by an assumed ~3.3x reduction that is a development target, not derived or cited'; carried exactly, printed results rounded once |
| `t_2028_box` | 280000 | Stated-not-derived | app05:94 | '2.8e5 (9.1e5 before the assumed 3.3x step reduction)': 9.095e5/3.3 = 2.756e5 at 20 steps (ruling 2026-10-11; was 1.6e5 at 12 steps) |
| `t_with_fermions_6x6_raw_quoted` | 850000 | Stated-not-derived | app05:91 | '8.5e5 T raw': 1973.0 x 36 x 12 = 8.523e5 |
| `t_with_fermions_6x6_quoted` | 260000 | Stated-not-derived | app05:91 | '~2.6e5 T with the same reduction applied to every sector' |
| `c_t_with_fermions_6x6_quoted` | 2000 | Stated-not-derived | app05:91 | 'c_T ~ 2.0e3' |
| `system_qubits_with_fermions_6x6_quoted` | 288 | Stated-not-derived | app05:91 | '288 system qubits' = 36 x (4 + 4) |

**Trotter schedule**

| field | value | provenance | source | note |
|---|---|---|---|---|
| `n_trot_2028` | 20 | Stated-not-derived | app05:59,97 | '20 steps at the same delta t', to t = 2/m_phi: both slope times (10 and 20 steps) inside the chapter's 1-2/m_phi exponential window (ruling 2026-10-11, was 12 steps with the early time 0.4/m_phi in the Zeno regime) |
| `n_trot_family` | 12 | Stated-not-derived | r16 record | the 8^2 family rows and the interim 4^2 variant (no longer printed) were set at 12 steps to t = 1.2/m_phi; kept there as records |
| `t_max_2033` | 10 | Stated-not-derived | app05:82,100 | t_max = 10/m_phi (units of 1/m_phi) |
| `dt` | 0.1 | Stated-not-derived | app05:82,91 | 'delta t = 0.1/m_phi' |
| `t_window_2028` | (1, 2) | Stated-not-derived | app05:91,101,142 | 't ~ 1-2/m_phi', Zeno-to-exponential crossover |

**Shots, readout and fault budget**

| field | value | provenance | source | note |
|---|---|---|---|---|
| `shots_2028` | 5000 | Stated-not-derived | r23 record | single-time 20% count at t = 1.2/m_phi (r23 box '5e3'); it fixes the ln P_FV decrement there, delta = 1/(N eps^2) = 5e-3, from which the slope shots follow |
| `decrement_ref_time_2028` | 1.2 | Assumed | derived here 2026-10-11 | time (1/m_phi) at which shots_2028 set the decrement; in the exponential regime the decrement grows linearly with t |
| `n_rare_range` | (1000, 10000) | Stated-not-derived | app05:83 | N_rare ~ 1e3-1e4 in Eq. Nshot_baryo |
| `shots_nucleation_2033` | 10000 | Stated-not-derived | app05:169 | nucleation arm '1e4 shots' per grid point |
| `eps_c` | 0.003 | Stated-not-derived | app05:63,115 | target precision on the per-bin asymmetry: a tenth of the comparator's typical bin (\|A_R\| 0.03-0.05; author decision 2026-10-12, option 2; was 1e-3, the size the asymmetry itself was then estimated at) |
| `n_bin` | 8 | Stated-not-derived | app05:63,115 | 'the 8 bins the 15-site axis resolves': floor(15/2) + 1 distinct \|k_z\| (ruling 2026-10-11, batch 4 item 3b; was 'N_bin ~ 1e2', which the lattice cannot resolve) |
| `mlae_queries_per_bin_quoted` | 7800 | Stated-not-derived | app05:63 | '~8e3 queries per momentum bin' at eps_c = 3e-3 (C_est = 1, x 1/p0^2 flag contrast; 2026-10-12, was '~7.0e4' at 1e-3) |
| `c_est` | 1 | Stated-not-derived | app05:63 | author-adopted branch prior \|p0 A_R\| < sin(pi/66) = 0.048: every circuit at m_max, Cramer-Rao constant one (Suzuki Fisher_final) |
| `dn_branch_prior` | 0.048 | Stated-not-derived | app05:63 | 'which we adopt as a prior rather than derive': \|p0 A_R\| < sin(pi/66) (author choice); the comparator's largest bin gives p0 A_R = 0.033 (2026-10-12) |
| `mlae_schedule` | single | Cited | arxiv_1904_10246 Manuscript_v2.tex:238 | all circuits at m_max; LIS m = 1, 3, ..., m_max (:268, :271) kept as sensitivity |
| `dn_estimate` | 0.056 | Assumed | sint/sint_results.json | COMPUTED (2026-10-12, sint/scatter.py, after the Bloch-sign fix of the verification pass): the largest per-bin \|A_R\| of the chapter's free-fermion comparator (two Wilson flavours, a m_f = 0.5, static tanh wall \|y\| phi = 0.5 of width 1/m_phi, Delta theta_C = 0.1) over the 8 x 8 incident (k_z, k_y) bins (0.059 at width 1 a); it changes sign across the axis. The 1+1D estimate of arxiv_2412_10365 gave ~1e-3 |
| `dn_estimate_1p1d_record` | 0.001 | Stated-not-derived | app05:61 | '~1e-3 at Delta theta_C = 0.1' (arxiv_2412_10365), the printed estimate until 2026-10-12 (record) |
| `dn_state_independent_bound` | 1 | Assumed | derived here | no derived bound on \|Delta\| follows from the chapter's inputs (the 'below Delta theta_C' scaling, taken as a bound, gives 0.1 > 0.048); recorded only: the single depth rests on the author prior (c_est), not on a derived bound |
| `m_max_quoted` | 33 | Stated-not-derived | app05:88,166 | 'm_max = floor(1e9/3.0e7) = 33' (fixed point) |
| `envelope_2033` | 1000000000 | Cited | DOE RFI 2026 | 1e9 hard-op envelope, app05:88,91,166 |
| `cap_2028` | 100000 | Cited | DOE RFI 2026 | 1e5 T cap, app05:128,181 |
| `accuracy` | 0.2 | Stated-not-derived | app05:49,160,201 | 20% relative accuracy on Gamma/V |
| `fault_budget_per_shot` | 0.1 | Assumed | main:272 (report-wide) | 0.1 expected faults per shot, report-wide |
| `eps_l_2028_quoted` | 4e-07 | Stated-not-derived | app05:96 | 'Required epsilon_l <~ 4e-7' (0.1 expected faults per 2.7e5-T shot; was 6e-7 at 12 steps) |
| `eps_l_2033_deepest_quoted` | 1e-10 | Stated-not-derived | app05:167 | '<~ 1e-10 for the deepest MLAE circuit (33 x 3.0e7 T)' |
| `eps_l_2033_base_quoted` | 3e-09 | Stated-not-derived | app05:167 | '<~ 3e-9 for the base circuit' |

**Wall time and campaign**

| field | value | provenance | source | note |
|---|---|---|---|---|
| `t_gate_s` | 1e-06 | Stated-not-derived | app05:89,91 | '1 us per T-gate' (seconds per T) |
| `shot_overhead_s` | 0.0001 | Assumed | app05:130 | per-shot overhead t0 ~ 0.1 ms: register initialization (~1 logical cycle, 10 us), final transversal readout and decode (63 us decoder latency, Google_QEC_below_threshold); the wall is gate-limited, with no separate shot rate |
| `n_scan_points` | 100 | Stated-not-derived | app05:89,151,176 | '~100-point (T, Delta theta_C) scan' |
| `n_dn_points` | 10 | Stated-not-derived | app05:89 | Delta n arm 'at ~10 representative points' |
| `horizon_yr` | 5 | Stated-not-derived | app05:89,151,175 | '5-year campaign horizon', one machine, serial (no machine count is formed) |

**Two-tier schedule (first result and campaign) and T-depth**

| field | value | provenance | source | note |
|---|---|---|---|---|
| `slope_times` | 2 | Assumed | shot_audit.json ch8 (1) | Gamma/V is the slope of ln P_FV past the Zeno time, so two times (t/3 and t), not one |
| `slope_early_fraction_2028` | 0.5 | Assumed | ruling 2026-10-11 | 2028 early time t/2: 10 of 20 steps (1.0 and 2.0/m_phi), both in the chapter's exponential window; shots per time follow from slope_shots (was t/3, 4 of 12 steps, with 0.4/m_phi in the Zeno regime) |
| `slope_early_fraction_2033` | 0.3333333333 | Assumed | shot_audit.json ch8 lever (5) | 2033 early time t/3: 33 of 100 steps; a long lever arm is ~5x cheaper than t/1.5 |
| `detection_sigma` | 3 | Assumed | H. Lamm 2026-10-02 | the 2033 first result is a 3-sigma detection of Delta n |
| `n_points_first` | 2 | Stated-not-derived | app05:143 | first result at T = T_c and 0.9 T_c (the box's two temperatures) |
| `s_int` | 0.0214 | Assumed | sint/sint_results.json | COMPUTED (2026-10-12, sint/scatter.py, after the Bloch-sign fix of the verification pass; 1.33e-2 before it): f_eq-weighted (Fermi-Dirac at a T = 0.5) mean of A_R over the incident bins for the comparator at \|y\| phi = 0.5, width 2a, Delta theta_C = 0.1; 2.4e-2 at 0.9 T_c; 1.9-2.7e-2 over wall amplitudes 0.3-1.0 and widths 1-2 a at both temperatures; linear in Delta theta_C, zero at 0; the C' relation R(hole_L) = R(p_R) holds to 1e-4. Converged in chain length and packet width. The first-result Delta-n wall scales as S_int^-2 |
| `s_int_band` | (0.019, 0.027) | Assumed | sint/sint_results.json | COMPUTED record: S_int over the comparator's wall amplitudes 0.3-1.0 and widths 1-2 a at T_c (1.9-2.4e-2) and 0.9 T_c (2.1-2.7e-2) |
| `label_register_first` | 7 | Assumed | shot_audit.json ch8 lever (1) | packet-label register of the f_eq-weighted packet superposition, taken from the 10-qubit workspace |
| `t_depth_2028` | (5400, 7000) | Assumed | factory.json ch08 run 1 | T-depth per 12-step shot, banked: ceil(1980/29) x 21.38 + ... per step (A = 29-30 RUS ancillas) to the serial structural schedule, x 12 / 3.3; the model scales it by n_trot_2028 / t_depth_2028_ref_steps |
| `t_depth_2028_ref_steps` | 12 | Assumed | factory.json ch08 run 1 | step count of the factory.json run behind t_depth_2028; depth is proportional to the step count |
| `t_depth_2033` | (720000, 2870000) | Assumed | factory.json ch08 runs 2-5 | T-depth per base circuit: A = 40 (bath ancillas borrowed during the closed scattering segments) to A = 10 (workspace only) |
| `t_depth_2033_first` | (870000, 9600000) | Assumed | factory.json ch08 runs 2-5, rescaled | T-depth per first-result Delta-n base circuit: the 7-qubit label leaves A = 33 (bath borrowed) to A = 3 (workspace only). Scaled from t_depth_2033 by the rotation-layer count per step, ceil(10440/A): 7.2e5 x 317/261 = 8.74e5 and 2.87e6 x 3480/1044 = 9.57e6 |

**State preparation, flag contrast, bath and fault model**

| field | value | provenance | source | note |
|---|---|---|---|---|
| `am_f` | 0.5 | Assumed | NEEDS_AUTHOR (flag contrast) | fermion Wilson mass a m_f in the symmetric phase: not fixed by the chapter, taken equal to a m_phi; flag contrast p0 = 0.57, 0.66, 0.77 at a m_f = 0.25, 0.5, 1 (flag_contrast) |
| `slater_givens` | (N - eta) eta | Cited | arxiv_1711_05395 Eq. n_gates (fhm_arXiv4.tex:346-350) | Givens rotations for one Slater determinant of eta particles in N modes |
| `fft8_rotation_free` | True | Cited | arxiv_1902_10673 main.tex:1102 | FFFT of side 8 needs no arbitrary rotations; F_2 = two pi/8 Pauli rotations = 2 T, twiddles powers of T (<= 1 T each), derived here |
| `beta_mphi` | 1 | Stated-not-derived | app05:99 | 'plasma correlation time ~1/T ~ 1/m_phi' at T ~ T_c: beta = 1/m_phi |
| `lr_patch_radius` | 2 | Assumed | derived here | Lieb-Robinson radius v_LR beta in sites: beta = 1/m_phi = 2a at a m_phi = 0.5, v_LR ~ 1 site per unit a; patch (2R+1)^2 = 25 sites |
| `gibbs_cost_model` | t beta | Cited | arxiv_2311_09207 Thm. L_cost and lattice estimate (main.tex:398-436) | Hamiltonian simulation time ~beta per unit Lindblad time on a patch of radius v_LR beta; jumps normalized \|\|sum_a A^a+ A^a\|\| <= 1 (main.tex:346), adjoints included |
| `idle_cycle_s` | 1e-05 | Assumed | derived here | one logical cycle ~10 us (as in shot_overhead_s); idle qubit-cycles = LQ x circuit time / cycle |
| `bath_mixing_sweeps` | (3, 20) | Assumed | estimate 2026-10-10 (scratch ch8_heur) | sweeps of the filtered dissipator to mix the Gibbs state, a dynamic-exponent heuristic (not measured): relaxation time tau = xi^z = 4-4.5 sweeps at xi = 1/m_phi = 2a and z = 2-2.17 (model A); warm start from a mean-field thermal state 2-3 tau = 3-14 sweeps, cold start to 1% local accuracy 7-21; printed as words only, the nucleation walls stay unpriced lower bounds (ruling 2026-10-11) |
| `bath_mixing_sweeps_central` | 8 | Assumed | estimate 2026-10-10 | central value of bath_mixing_sweeps (record; nothing printed uses it) |

**Finite-volume saddle behind the nucleation decrement**

| field | value | provenance | source | note |
|---|---|---|---|---|
| `sigma_over_mphi2` | 1.32 | Assumed | derived here 2026-10-10 (scratch ch8_heur) | wall tension sigma / m_phi^2 of the 2033 scalar potential, chosen so that the wall-pair (slab) saddle on the 15x8 cylinder reproduces the decrement 2.5e-3 at t = 10/m_phi behind the printed 1e4 single-time shots; thin-wall phi^4: sigma = m^3/(3 lambda), so lambda/m_phi = 0.25. The chapter never fixed sigma or lambda; this is the stated assumption |
| `latent_heat_ell` | 1 | Assumed | assumed 2026-10-10 | free-energy difference Delta F = ell m_phi^3 (1 - T/T_c) between the wells (thin wall, T_c = m_phi from beta_mphi); lowers the slab action by ell (1 - T/T_c) L_perp over a strip of width 1/m_phi |

**Interim scalar+Wilson variant (printed values)**

| field | value | provenance | source | note |
|---|---|---|---|---|
| `t_interim_quoted` | 140000 | Stated-not-derived | app05:178,181 | '~1.4e5 T raw (c_T ~ 7.5e2: 2.5e2 + 5.0e2, x 16 x 12)' |
| `gap_interim_quoted` | 1.4 | Stated-not-derived | app05:178,181 | '~1.4x above the 2028 hard-op cap' |
| `t_interim_banked_quoted` | 44000 | Stated-not-derived | app05:181 | 'brings it to ~4.4e4 if applied to every sector' |
| `t_interim_scalar_share_quoted` | 110000 | Stated-not-derived | app05:181 | '~1.1e5 if applied to the scalar share alone, below the benchmark's 1.6e5' |

## 4. How the numbers are built

The module docstring (lines 1-210) narrates the same derivation. Its closing paragraph (lines 200-209) gives walls at unit flag contrast, without the priced preparation and reflections. The computed intermediates and the tests include both and carry the numbers the paper prints; those are quoted here.

### 4.1 Qubits

`qubits_per_site(n_phi, n_f, components)` (lines 467-469) gives log₂N_φ + N_f × components per site. LQ is then V × that count + N_anc. This is the chapter's Eq. (Nq_baryo), with the 2N_f term as the 2-component case.

| instance | per site | system | ancilla | LQ | test |
|---|---|---|---|---|---|
| 2028, 6², N_φ=16 | 4 | 144 | 30 (Trotter control, basin projector) | **174** | `test_2028_register_is_a_component_sum` |
| 8², N_φ=4 / N_φ=16 | 2 / 4 | 128 / 256 | 30 | 158 / 286 | `test_2028_family_rows` |
| interim 4², N_φ=4, N_f=2 | 6 | 96 | 17 (MLAE workspace, a budget, not itemized) | 113 (printed ~115) | `test_interim_variant` |
| 2033, 15×8, N_φ=16, N_f=2 | 8 | 960 | 40 (30 bath + 10 workspace) | **1000** | `test_2033_register_is_a_component_sum` |

At N_φ = 16 the fermion arm on 6² needs 288 system qubits. A 4-component spinor would put the 2033 system register at 1440 qubits (`four_component_system_qubits`; `test_2033_four_component_spinor`).

### 4.2 Gate content of one site per Trotter step

Every count is for one site and one Trotter step. Each term is either a synthesized rotation or exact T.

**Scalar**, `scalar_counts` (lines 491-500). n_q = log₂N_φ qubits, with φ linear in Z on the centered grid.

- On-site polynomial (mass, φ⁴, bias, and the gradient's φ² share): one R_Z per Z-string of weight 1 to min(4, n_q), so Σ_{w≤4} C(n_q, w). That is 15 at n_q = 4 and 3 at n_q = 2 (arXiv:2407.13819 Eq. rzTrotter; arXiv:2210.07985 Eqs. trdispl to trphi4).
- Gradient φ_iφ_j: n_q² ZZ rotations per link, with 2 links per site. That is 32 at n_q = 4.
- Kinetic π²: conjugating by the centered QFT turns it into C(n_q, 2) ZZ rotations, so 6 at n_q = 4 (arXiv:2210.07985 Eqs. fftqft, trpi2). The diagonal phase matrices around the QFT commute through or merge into the weight-1 on-site rotations, so they add nothing.
- Two exact QFTs per site per step, from `qft_cost` (lines 472-488), counted by controlled phase π/2^k at distance k. A k = 1 phase is a controlled-S at 3 T. A k = 2 phase is a Toffoli into an ancilla plus one T, 8 T. Each k ≥ 3 phase is a Toffoli plus one synthesized rotation. At n_q = 4 one QFT is 3×3 + 2×8 + 1×7 = 32 T plus 1 rotation. The AQFT formula of arXiv:2407.13819 would give more than 500 T per 4-qubit QFT and is not used (`test_scalar_counts_are_the_papers_term_counts`).
- Scalar total at N_φ = 16: **55 rotations + 64 exact T**. At N_φ = 4: 12 rotations + 6 T.

**Fermions**, `fermion_counts` (lines 503-516). Two Wilson flavors at r = 1 under Jordan-Wigner. The JW strings are parity ladders and cost CNOTs only.

- Hopping: at r = 1 the hopping matrix (σ₃ − iσ_j)/2 has rank one. A π/4 Givens basis change per site per direction (2 T, in and out) leaves one c†c + h.c. per link per flavor, which is 2 R_Z. Per site this is 8 rotations, plus 16 exact T for the basis change.
- Mass plus Wilson term: 2 R_Z per flavor, so 4. This equals the Fermion-Primitives count 2N at d = 2 (`test_wilson_mass_matches_fermion_primitives`).
- Yukawa: per spinor component, φ ⊗ (XX+YY)/2 gives 2n_q rotations, plus 2 for the θ_C(z) phase frame. Per site that is 16 + 4 at n_q = 4 and 8 + 4 at n_q = 2.
- Fermion total: **32 rotations + 16 T** at N_φ = 16, and 24 rotations + 16 T at N_φ = 4. A 4-component spinor doubles every fermion count (`test_fermion_counts_are_derived_here`).

### 4.3 From counts to c_T and T per circuit

`price_instance` (lines 519-572) multiplies the per-site-step counts by V × N_Trot. It then applies the report's synthesis rule from `estimates/common.py`. The per-rotation tolerance is ε_rot = √(ε_syn / N_rot), with ε_syn = 10⁻² per circuit (`eps_rot_for`, line 172). The cost per rotation is t_rot = 1.15 log₂(1/ε_rot) + 9.2 T (`t_per_rotation`, line 125), and a Toffoli is 7 T (`toffoli_t`, line 201). Then

  c_T = n_rot × t_rot + n_T,  T_circuit = c_T × V × N_Trot.

N_rot counts the synthesized rotations in the whole circuit, so c_T varies between instances of the same Hamiltonian through t_rot. The `budget_circuits` argument covers the deep MLAE circuits: it synthesizes rotations to the budget of a circuit that contains m applications of the base circuit, plus `extra_rot` preparation rotations per application.

| instance | site-steps | N_rot | ε_rot | t_rot | c_T | T raw |
|---|---|---|---|---|---|---|
| 2028 6², N_φ=16 (20 steps) | 720 | 39,600 | 5.025e-4 | 21.802 | 1263.1 | 9.095e5 |
| 2028 6², N_φ=16 at 12 steps (r12 record) | 432 | 23,760 | 6.487e-4 | 21.379 | 1239.8 | 5.356e5 |
| 8², N_φ=4 | 768 | | | | 253.1 | 1.944e5 |
| 8², N_φ=16 | 768 | | | | 1266 | 9.723e5 |
| 6², N_φ=16 + fermions | 432 | | | | 1973 | 8.523e5 |
| interim 4², N_φ=4 + fermions | 192 | | | | 754.8 (250.3 scalar + 504.5 fermion) | 1.449e5 |
| 2033 nucleation, own budget | 12,000 | 1,044,000 | 9.787e-5 | | 2212.9 | 2.656e7 |
| 2033 Δn, depth-33 budget | 12,000 | 33 × 1,060,672 = 35,002,176 | 1.690e-5 | 27.430 | 2466.4 (1572.7 scalar + 893.8 fermion) | 2.960e7 |

All values are from the `--intermediates` output. They are pinned by `test_ct_is_priced_by_rtol_per_instance` (c_T 1263.1 at 20 steps, 1239.82 at 12, and 2466.43; N_rot budget 33 × (1,044,000 + 1800 + 14400 + 472)), by `test_ct_range_is_the_instances` (the printed span 2.5e2 to 2.5e3), and by the family, interim and fermion-arm tests.

### 4.4 The 2028 box (`_model_2028`, lines 645-822)

- **Per-shot T.** The raw count is c_T × 720 = 9.095e5 (20 steps; ruling 2026-10-11, was 5.356e5 at 12). The box divides it by an assumed 3.3× c_T reduction (`ct_reduction_banked`). That factor is a development target, neither derived nor cited, and it is carried exactly: 9.095e5 / 3.3 = **2.756e5** (`resources.json`: 275591.4), printed 2.8e5. In the breakdown, every primitive's T-per-unit is divided by the same 3.3 (rotations 21.802/3.3 = 6.607 T; QFT exact T 64/3.3 = 19.39 per site-step):

  | primitive | count | T each | T total |
  |---|---|---|---|
  | scalar_onsite_polynomial | 10,800 | 6.607 | 7.135e4 |
  | scalar_gradient_links | 23,040 | 6.607 | 1.522e5 |
  | scalar_kinetic_pi2_zz | 4320 | 6.607 | 2.854e4 |
  | scalar_qft_controlled_phase | 1440 | 6.607 | 9514 |
  | scalar_qft_exact | 720 | 19.39 | 1.396e4 |
  | adiabatic_state_prep | 1 | 0 | 0 (UNSOURCED) |
  | basin_projector_and_readout | 1 | 0 | 0 (UNSOURCED) |

  The sum is 2.756e5 (`test_2028_breakdown_is_primitive_level_and_compiled`).
- **ε_l.** `eps_l_required` (lines 632-634) gives ε_l = 0.1 / T = **3.63e-7**, printed ≲ 4e-7.
- **Shots.** The slope of ln P_FV comes from two times, t/2 (10 steps, 1.343e5 T after the 3.3) and t (20 steps), both inside the chapter's 1-2/m_φ exponential window. The shots follow from the lever arm, `slope_shots`: N = (1 + r)/(ε² δ₂ (1 − r)²) at each time, with r = t₁/t₂ and δ₂ the decrement at the late time (binomial statistics, decrement growing linearly with t). The r23 single-time count 5e3 at 20% fixed δ(1.2/m_φ) = 5e-3 (`shots_2028`, `decrement_ref_time_2028`), so δ(2.0/m_φ) = 8.33e-3 and r = 1/2 gives 6 × the single-time count: **1.8e4 shots at each time**, `Result.shots` = 3.6e4 (`slope_shot_factor` 6; the r25 schedule's factor 3 is the r = 1/3 case of the same formula).
- **Wall time.** 1.8e4 × [(2.756e5 + 1.343e5) × 1 µs + 2 × 0.1 ms] = 7382 s = **2.05 h** per parameter point, printed ~2.1 h. At the unreduced c_T it is 6.8 h. One shot takes 0.276 s (`shot_s`).
- **Depth.** The T-depth band 5.4e3 to 7.0e3 is an Assumed input from the factory analysis at 12 steps (`t_depth_2028_ref_steps`); the model scales it by 20/12 to 9.0e3 to 1.17e4 (`t_depth_scaled`). `depth_exports` (common.py line 316) turns it into F* = N_T/D_T = 23.6 to 30.6, printed 24-31. Since F* ≥ 10, the ten-factory baseline (section 8) behind 1 µs per T is fed. The depth floor is 2.4e3 to 3.1e3 s.
- **Record.** The decrement 8.33e-3 at 2.0/m_φ on 6² (area 9/m_φ²) corresponds to a vacuum bounce action of about 8.0 (`bounce_action_implied_2028`, `bounce_action_implied`), the 2028 counterpart of the 2033 saddle below.

### 4.5 The 2033 target (`_model_2033`, lines 1134-1530)

**One application of A.** Per site-step the circuit has 87 rotations and 80 exact T, over 120 × 100 = 12,000 site-steps. The evolution is therefore c_T × 12,000 = 2.960e7 T, or 2.96e5 T per step. On top of that, `free_prep` (lines 930-940) prices the free part of the preparation for every application:

- scalar product state: N_φ − 1 = 15 rotations per site, 1800 rotations, 4.937e4 T (`scalar_product_rotations`);
- free-fermion vacuum of the wall background: the transverse axis is periodic, so the problem is block-diagonal in k_y. That gives 8 sectors of N = 60 modes at half filling, each needing (N − η)η = 900 Givens rotations (arXiv:1711.05395), at 2 rotations each: 14,400 rotations, 3.95e5 T (`slater_rotations_ky`);
- inverse 8-point fermionic FFT on 60 chains: rotation-free (arXiv:1902.10673), 32 T per chain, 1920 T (`fft_transverse_t`);
- selector-controlled packet in one k_y sector, ≤ 59 Givens per charge: 472 rotations, 1.295e4 T (`packet_rotations`).

The free preparation totals 4.592e5 T, or 1.55% of the evolution (`r_free`). One application then costs t_A = **3.006e7** T, printed "Base circuit 3.0e7".

**Reflections.** `mlae_reflection_t` (lines 871-877) prices each Grover iterate. S_χ counts ≤ 480 fermion occupations into a 9-bit register using ≤ 10 Toffolis per increment, computed and uncomputed. S₀ is a controlled phase on the 1000-qubit input, at LQ − 2 Toffolis. Together that is (2 × 480 × 10 + 998) × 7 = **74,186 T** per iterate, 0.25% of the evolution (`test_b2_construction_priced`).

**Depth cap and the headline.** `dn_schedule` (lines 943-983) finds the deepest circuit that fits in 10⁹ T. The cost of a depth-m circuit is m t_A(m) + (m − 1)/2 × refl, where t_A(m) is priced with every rotation synthesized to that circuit's own budget. Because t_A rises with m, the depths that fit form an interval. The scan starts at m₀ = ⌊(10⁹ + refl/2)/(t_A(1) + refl/2)⌋ = 37 and steps down: 37, 36, 35 and 34 exceed 10⁹ (1.1171e9, 1.0861e9, 1.0550e9, 1.0240e9), and 33 fits at 9.9305e8. The headline, `Result.hard_ops`, is the deepest executed circuit:

  33 × 3.0056e7 + 16 × 74,186 = **9.930e8 T** (`resources.json`: 993049073.3), printed ≳ 9.9e8.

It is a lower bound because the interacting dressing ramp is not priced (section 7). The primitive breakdown multiplies each per-application count by 33. The largest pieces are the scalar gradient (3.476e8), the Yukawa bilinear (1.738e8) and the on-site polynomial (1.629e8). Preparation adds the Slater determinant (1.303e7), scalar loading (1.629e6), packet (4.27e5) and FFT (6.3e4), and the reflections add 1.187e6. The sum equals the headline (`test_2033_t_is_the_derived_product`, `test_2033_m_max_is_the_fixed_point`, `test_2033_breakdown_shows_the_unpriced_pieces`).

**ε_l.** At 0.1 expected faults, ε_l = 0.1/9.93e8 = **1.01e-10** for the deepest circuit and 0.1/2.96e7 = 3.38e-9 for one evolution (`test_epsilon_l_follows_r3`).

**Δn queries.** The flag fires on a threshold of the incident-side charge. The vacuum has no sharp value of that charge, so the good-outcome amplitude is a = ½ + ½ p₀ A_R. `flag_contrast` (lines 1025-1055) computes p₀ = P(q = 0) from the full counting statistics of the free vacuum restricted to the incident half, sector by sector in k_y. At a m_f = 0.5 it gives **p₀ = 0.6584**, which multiplies the query count by 1/p₀² = 2.31.

A second check covers the background. `c_prime_residual` (lines 1093-1102) confirms that charge conjugation combined with the flavor and spinor swap commutes with H on a C-violating wall; the residual is exactly 0. The vacuum background δ therefore vanishes, the flavor-summed Δn is zero, and the readout uses flavor-L packets.

With every circuit at depth 33 the estimator constant is 1. Writing a = (1 + Δ)/2 with Δ = p₀A_R, this is the Cramér-Rao bound from Suzuki's Fisher information m²/(1 − Δ²) per circuit, and it requires |Δ| < sin(π/66) = 0.0476 (`balanced_branch_bound`, line 836; intermediate `mlae_branch_bound_needed`). The authors adopt |A_R| < 0.048, which implies it, as a prior; the ~10⁻³ estimate sits 47.6 times inside it (`dn_estimate_headroom`). The model checks that no derived bound supplies it (`mlae_derived_bound_ok = False`, `mlae_author_prior_ok = True`; `test_branch_condition_and_author_prior`). Then:

- queries per bin = 1/(p₀² ε_C² m) = 1/(0.6584² × 9×10⁻⁶ × 33) = **7.77e3**, printed ~8e3 (ε_C = 3e-3 since 2026-10-12, a tenth of the comparator's typical bin; at the earlier 1e-3 it was 6.99e4);
- 8 bins (`n_bin`, the distinct |k_z| on the 15-site axis, floor(15/2) + 1 = 8, checked in `__post_init__`; ruling 2026-10-11, batch 4 item 3b, was 100) give 6.21e4 queries in **1.88e3 circuits** of depth 33 (`Result.shots`);
- this is 23.3 times the Heisenberg count ε_C⁻¹ N_bin = 2.67e3 (`test_2033_mlae_schedule`).

**Walls.** All walls are serial on one machine. Each circuit costs T × 1 µs + 0.1 ms; `circuit_s` is at lines 986-988.

- One depth-33 circuit takes 993 s (~17 min; `deepest_circuit_s`; with the 2028 `shot_s` pinned by `test_maximum_circuit_time`), and one application of A takes 30.06 s (`base_circuit_s`).
- One campaign Δn point is 1.88e3 × 993 s = **0.059 yr** (`dn_arm_yr`; `Result.wall_time_s` = 1.87e6 s).
- The nucleation circuits are single circuits, so each is synthesized to its own budget. The shots follow `slope_shots` at the decrement the box's 1e4 single-time count fixes, δ(10/m_φ) = 2.5e-3 (`nucleation_decrement_late`): 2.96e4 at 33 steps (8.446e6 T) and 2.96e4 at 100 steps (2.656e7 T), printed 3e4, which comes to **12.0 d** at T_c. The 100-point scan takes 3.29 yr.
- **Finite-volume saddle** (`slab_saddle`; ruling 2026-10-11, printed as one sentence). On the 15×8 cylinder no critical bubble fits: a thin-wall bubble of tension s = σ/m_φ² against ΔF = ℓ m_φ³ (1 − T/T_c) has a minimum S/T = 4πs²/ℓ = 21.9, above the wall pair's 2sL_⊥ = 10.6, so the saddle is a pair of walls wrapping the 8-site direction and Γ/V is finite at T_c. σ = 1.32 m_φ² (`sigma_over_mphi2`, with `latent_heat_ell` = 1) reproduces the assumed decrement (2.52e-3 against 2.5e-3, `slab_reproduces_assumed_decrement`); in thin-wall φ⁴ it means λ/m_φ = 0.25. At 0.9 T_c the action over T grows to 11.3 and the decrement halves (`decrement_ratio_0p9tc` 0.50), so that nucleation point costs 2 × the shots, **24.1 d** (`nucleation_point_days_0p9tc`).
- The campaign is 10 × 0.059 + 3.29 = **3.88 yr**, 0.78 of the 5-year horizon (`wall_campaign_s` 1.22e8 s; 8.62 yr at ε_C = 1e-3).
- The first result runs at depth 31, not 33. MLAE depths are odd (`dn_schedule` returns the largest odd depth that fits as `m_top`). The 7-qubit label register drives 2⁷ label-controlled packets, which add 1.657e6 T per application. With the label packets t_A = 3.169e7, so depth 33 no longer fits in 10⁹ and the deepest odd depth that fits is 31 (`first_m_top`). Here S_int is the f_eq-weighted integral of A_R in the chapter's Eq. (deltan), taken at 10⁻³ (`s_int`, not computed). Each temperature then needs 9/(p₀² S_int² × 31) = 6.70e5 queries in 2.16e4 circuits, or 246.0 d. The first result is 2 × 246.0 + 12.0 + 24.1 d = 528 d = **1.45 yr** (`wall_first_result_s` 4.56e7 s; the 0.9 T_c nucleation point carries the slab factor 2). The nucleation scan plus the first result fits in **4.63 yr** (the two first-result nucleation points are already in the scan) (`test_2033_wall_time_and_campaign`, `test_2033_first_result_tier`).

**Depth and factories.** The T-depth per evolution, 7.2e5 to 2.87e6, is an Assumed input from the factory analysis. The low end is A = 40, which borrows the 30 bath ancillas as rotation workspace; no Δn circuit needs the bath. The high end is A = 10, the workspace alone. This gives F* = 10.3 to 41.1, so the ten-factory baseline holds at both ends. A 1-year campaign would need about 86 factories (`factories_for_1yr`; 700 at 100 bins).

The first-result label leaves A = 3 of the workspace. Rescaled by rotation layers per step, its T-depth band is 8.7e5 to 9.6e6, with F* = 3.1 to 34. Without borrowing the bath, the first result is depth-bound at 4.47 yr (`test_2033_depth_exports`, `test_2033_first_result_depth`).

### 4.6 Sensitivities the model computes but the boxes do not bank

- **Suzuki's linearly incremental schedule** (LIS, m = 1, 3, …, 33) instead of the single depth: query constant 33 × 289/6545 = 1.457. That gives 1.019e5 per bin, 0.78 yr per Δn point and 11.0 yr for the campaign (9.70 and 100.3 at 100 bins) (`lis_query_constant`, `test_lis_constant_is_suzukis_fisher_ratio`).
- **Bath mixing band** (`bath_mixing_sweeps` 3-20, `bath_mixing_sweeps_central` 8, `bath_mixing_sensitivity`; ruling 2026-10-11, words only). A dynamic-exponent heuristic (relaxation time τ = ξ^z = 4-4.5 sweeps at ξ = 1/m_φ = 2a and z = 2-2.17; warm start from a mean-field thermal state 2-3 τ = 3-14 sweeps, cold start to 1% local accuracy 7-21) puts the mixing time at 3-20 sweeps of the filtered dissipator. Every shot ends in a projective measurement, so every nucleation shot would pay the preparation: at 3 / 8 / 20 sweeps a nucleation shot is 4.2e8 / 1.1e9 / 2.7e9 T (the central value no longer fits 10⁹), a nucleation point 290 d / 2.0 yr / 5.1 yr, the 100-point scan 78 / 203 / 502 yr and the first result 3.7 / 7.4 / 16 yr. The chapter prints the band and the low-to-high consequence in words; the nucleation walls stay unpriced lower bounds (`bath_prep_dominates_nucleation_shot`, `test_bath_mixing_band_is_words_only`).
- **Dressing ramp** costing r × the evolution (`prep_sensitivity`, lines 1115-1128). At r = 0.2 the depth drops to 27 and a Δn point takes 0.78 yr. At r = 1 the depth is 15 and a point takes 2.32 yr (9.74 and 29.0 yr at 100 bins). Depth 33 leaves 0.71 Trotter steps of headroom per application, so any r > 0.007 lowers the depth (`test_b3_preparation_sensitivity`, `test_dn_schedule_depth_never_rises_with_r`).
- **Flag contrast in the wall background:** p₀ = 0.654, 0.646 and 0.603 at |y|φ_wall = 0.3, 0.5 and 1.0. The free p₀⁻² = 2.31 is therefore a lower bound on the query factor, which can reach 2.75.
- **c_T reduction applied to 2033:** banking the 3.3× would give 8.97e6 T per evolution. It is recorded but not banked.

## 5. What the paper prints

The chapter's numbers come from `PUBLISHED` (lines 1545-1555), the box values the model is checked against within a 5% tolerance, and from `INSTANCE_ROWS` (lines 1558-1582), which writes `resources.json`. The tests check the model against `PUBLISHED`; the model never reads it (`test_model_does_not_read_published`).

**`resources.json` rows for chapter 8** (current values):

| era | label | LQ | T | conditional | reconciled |
|---|---|---|---|---|---|
| 2028 | scalar bubble 6² | 174 | 162303.809916798 | false | true |
| 2028 | family 8², N_φ=4 | 158 | 58906.78463287797 | false | false |
| 2028 | family 8², N_φ=16 | 286 | 294649.4585081366 | false | false |
| 2028 | interim scalar+Wilson 4² | 113 | 144912.76692923915 (raw, no reduction) | true | false |
| 2033 | wall scattering 15×8 | 1000 | 993049073.279399 | true | true |

Table 1.1 of the paper prints the two headline rows as "Scalar bubble, 6²; 174 LQ, ≳1.6×10⁵ T (conditional)" and "Scalar+Wilson+Lindblad, 15×8; 1000 LQ, ≳9.9×10⁸ T". The check of Table 1.1 against `resources.json` (`scripts/estimates/tests/test_table.py` in the paper tree) lives with the paper source, not in this repository. It checks the ≳ relation and the value for (8, "2028"); it does not check the conditional mark.

**Box values and the tests that pin them**

| printed | model value | test |
|---|---|---|
| 2028 LQ 174 | `lq` 174 | `test_2028_register_is_a_component_sum`, `test_2028_box_row_closes`, `test_instance_rows` |
| 2028 per-shot T 2.8e5 (c_T derived, /3.3; 20 steps) | 2.756e5; c_T 1263.1; raw 9.095e5 | `test_2028_banked_factor_is_assumed_and_carried_exactly`, `test_2028_raw_t_is_the_derived_product`, `test_ct_is_priced_by_rtol_per_instance` |
| 2028 T-depth 9.0e3-1.2e4, F* 24-31 | `t_depth_scaled` 9.0e3-1.17e4; `f_star` 23.6-30.6 | `test_2028_depth_exports` |
| 2028 ε_l ≲ 4e-7 | 3.63e-7 | `test_epsilon_l_follows_r3` |
| 2028 shots 1.8e4 at each of 10 and 20 steps | `slope_shots` at r = 1/2, δ₂ = 8.33e-3; 3.6e4 total | `test_2028_slope_two_times` |
| 2028 wall ~2.1 h (6.8 h unreduced) | 123.0 min; 405.9 min | `test_2028_slope_two_times` |
| 8² rows 158 LQ / 5.9e4 T and 286 LQ / 2.9e5 T | as in the table above | `test_2028_family_rows` |
| interim ~115 LQ, ~1.4e5 T raw, 4.4e4 / 1.1e5 with the reduction | 113; 1.449e5; 4.391e4 / 1.114e5 | `test_interim_variant` |
| fermion arm 288 qubits, c_T ~2.0e3, 8.5e5 raw, ~2.6e5 | 288; 1973; 8.523e5; 2.583e5 | `test_2028_fermion_arm_numbers_are_facts` |
| 2033 LQ 1000 | 960 + 40 | `test_2033_register_is_a_component_sum` |
| 2033 deepest circuit ≳9.9e8 T at depth 33 | 9.930e8; m_max 33 | `test_2033_t_is_the_derived_product`, `test_2033_m_max_is_the_fixed_point` |
| 2033 base circuit 3.0e7 T; c_T ≈ 2.5e3 | t_A 3.006e7; 2466.4 (1572.7 + 893.8) | `test_2033_t_is_the_derived_product`, `test_ct_is_priced_by_rtol_per_instance` |
| 2033 T-depth 0.72-2.9e6, F* 10-41 | 10.3-41.1 | `test_2033_depth_exports` |
| 2033 ε_l ≲ 1e-10 | 1.01e-10 | `test_epsilon_l_follows_r3` |
| flag contrast p₀ = 0.66 | 0.6584 | `test_b2_flag_contrast_and_c_prime` |
| Δn shots: 1.5e3 queries at depth 31 (first result); 6.2e4 over 8 bins in 1.9e3 circuits at depth 33 (campaign) | 1.462e3; 6.214e4; 1.883e3 | `test_2033_first_result_tier`, `test_2033_mlae_schedule` |
| nucleation 3e4 shots at each of two times; decrement 2.5e-3 at 10/m_φ from a wall-pair saddle with σ ≈ 1.3 m_φ²; ×2 at 0.9 T_c | 2.963e4 at 33 and 100 steps; `slab_decrement_tc` 2.52e-3; 24.1 d | `test_2033_wall_time_and_campaign` |
| walls: first result ≳37 d; campaign ≳3.9 yr, inside the horizon; scan + first result ≳3.3 yr | 37.1 d; 3.88, 0.78; 3.29 yr | `test_2033_first_result_tier`, `test_2033_wall_time_and_campaign` |
| any ramp r > 0.007 forces a shallower schedule | 0.007 | `test_2033_m_max_is_the_fixed_point` |
| bath sweep ≳1.3e8 T; mixing 3-20 sweeps (words): nucleation shot 0.4-2.7e9 T, scan 78-500 yr | 1.328e8; `bath_mixing_sensitivity` | `test_b5_filtered_bath_estimate`, `test_bath_mixing_band_is_words_only` |
| ε_l ≲ 1e-12 per location with idles; ≲ 7e-14 to keep the fault bias below the first-result error | 9.97e-13; 6.87e-14 | `test_g4_fault_bias_model` |

`test_referee_text_hooks` and parts of the tests above also check that the quoted strings appear in the chapter text. Those checks run only with the paper source present.

## 6. Circuit status

From `docs/CIRCUIT_STATUS.md` (Ch. 8: COMPILED 20, SCALING 0, CONJECTURE 0, UNSOURCED 5).

- **COMPILED, 2028 (5):** `scalar_onsite_polynomial`, `scalar_gradient_links` and `scalar_kinetic_pi2_zz`, all three term-level counts from arXiv:2407.13819 and arXiv:2210.07985, plus `scalar_qft_controlled_phase` and `scalar_qft_exact`, both derived here.
- **COMPILED, 2033 (15):** the same five scalar terms; the five fermion terms derived here (`fermion_wilson_hopping`, `fermion_mass_wilson_onsite`, `fermion_yukawa_phi_bilinear`, `fermion_yukawa_phase_frame`, `fermion_hopping_basis_change`); `vacuum_prep_scalar_product`; `vacuum_prep_slater_ky` (arXiv:1711.05395); `vacuum_prep_fft_transverse` (arXiv:1902.10673); `packet_prep`; and `mlae_reflections`.
- **UNSOURCED, priced at 0 T (5):**
  - 2028: `adiabatic_state_prep` and `basin_projector_and_readout`.
  - 2033: `vacuum_prep_dressing` (the interacting ramp of the Δn preparation), plus `lindblad_dissipator_step` and `gibbs_state_prep`. The last two belong to the nucleation circuits only.
- No primitive is SCALING or CONJECTURE.

"Compiled" here means a gate-level count, of rotations and exact T, of the Trotterized Hamiltonian and the free preparation. It does not mean a full compiled circuit has been verified.

## 7. Open items and limitations

Taken from `docs/OPEN_ITEMS.md` and the chapter's algorithm-development gap. None of them changes a printed number.

- **Unpriced preparation makes the 2033 numbers lower bounds.** The free part of the Δn preparation is priced (4.6e5 T per application). The interacting dressing ramp is not, and any ramp lowers the MLAE depth: r = 0.2 gives depth 27 and 9.7 yr per campaign point, and r = 1 gives depth 15 and 29 yr (`prep_sensitivity`). A static wall needs a degenerate well or a pinning term.
- **The 2028 count is neither a lower nor an upper bound.** It includes the 3.3× c_T reduction, an assumed development target (`ct_reduction_banked`). It also omits the adiabatic preparation and the basin projector. If the reduction cuts T-count but not T-depth, F* falls to 7-9 and the 2028 wall becomes depth-bound at 1.1-1.4 times the quoted 2.1 h.
- **The 2028 run is scalar only.** N_φ = 4 cannot resolve the biased well. At N_φ = 16 the fermion arm needs 288 system qubits on 6², above the 2028 register reference, at c_T ≈ 2.0e3. It waits for 2033.
- **The flag contrast depends on a fermion mass the chapter does not fix.** p₀ = 0.66 at a m_f = 0.5 (taken equal to a m_φ) and 0.57-0.77 for a m_f = 0.25-1. On the wall background p₀ = 0.60-0.65, so the query factor 2.3 is a lower bound and can reach 2.75.
- **Δn is flavor-resolved by necessity.** An exact symmetry (charge conjugation with the flavor and spinor swap, `c_prime_residual` = 0) makes the flavor-summed asymmetry vanish identically.
- **The bath is estimated, not priced.** One unit of Lindblad time is about 5.5e5 T (`filtered_jump_t`). One application of each of the 240 scalar jumps is 1.3e8 T, five nucleation circuits, and 7.3 such sweeps fit in 10⁹. A dynamic-exponent heuristic puts the sweeps needed to mix at 3-20 (`bath_mixing_sweeps`, printed as words); at any value in the band the preparation dominates the nucleation shot (section 4.6). The dissipator rate during evolution is open. The Gibbs preparation and the basin projector are carried at 0 T.
- **S_int comes from the free-fermion comparator.** `sint/scatter.py` (2026-10-12, after the Bloch-sign fix of the verification pass) gives S_int = 2.14e-2 at T_c (2.4e-2 at 0.9 T_c; 1.9-2.7e-2 over wall amplitudes 0.3-1 and widths 1-2 a at both temperatures) and per-bin |A_R| up to 0.056 with sign changes, for a static prescribed tanh wall in the free theory; the interacting, plasma-dressed wall is what the calculation determines. The first-result Δn wall scales as S_int⁻². The 1+1D estimate 10⁻³ is kept as a record (`dn_estimate_1p1d_record`).
- **The Δn readout rests on an adopted prior.** The single-depth schedule (estimator constant 1) uses the adopted prior |p₀A_R| < 0.048, and no bound is derived from the chapter's inputs; the comparator's largest bin gives p₀A_R = 0.037, headroom 1.3. Without the prior, Suzuki's LIS schedule costs 1.46 times the queries.
- **Fault bias.** At the 0.1-fault budget the deepest circuit has p_f = 0.095, which can shift A_R by up to 4.4e-3, above the per-bin 3×10⁻³. The bias is below the first-result error (7e-3) but above the campaign's per-bin 3e-3; keeping it below the per-bin target needs ε_l ≲ 7e-11 per T, or ≲ 7e-13 per location once 9.9e10 idle qubit-cycles are counted. A null circuit (the flavor-R antiparticle branch, a = ½ exactly) quadruples the first-result Δn queries: 0.11 yr for the first result, 8.7 yr with the scan. The fault budget counts T gates only; idle, Clifford and measurement faults are outside it.
- **Physics of the schedule.**
  - The 2028 slope now uses 1.0 and 2.0/m_φ (ruling 2026-10-11), both in the chapter's exponential window; the expected decrement there (5e-3 at 1.2/m_φ, from the r23 count) is still the chapter author's to confirm.
  - The decrement 2.5e-3 at T_c rests on the stated assumption σ = 1.32 m_φ² for the wall-pair saddle (`slab_saddle`); Γ/V is finite at T_c on 15×8 because no critical bubble fits. The 0.9 T_c point carries the factor 2 that follows.
  - The campaign uses the 8 momentum bins the 15-site axis resolves (ruling 3b). A scan over temperature alone (10 points) would cut the nucleation scan from 3.3 yr to about 0.3 yr and is the chapter author's option; if the Gibbs preparation were priced at the mixing band, the nucleation scan (78-500 yr), not the Δn arm, would set the campaign.
- **First-result workspace.** The first-result walls assume the 30 bath ancillas serve as rotation workspace during the closed Δn circuits. Without them, F* is 3.1 and the first result takes 4.5 yr.
- **Boundaries.** c_T counts 2 links per site, as for periodic boundaries. The 15×8 cylinder is open along z and has 8 of 240 links fewer, so the 2033 count is slightly conservative (~1% of c_T).
- **Not estimated:** the tensor-network cost of the classical comparator for the 2+1D quantum scalar.

## 8. Conventions used

This chapter follows the conventions shared by every model, described in the top-level `README.md` and implemented in `estimates/common.py`.

- **R-TOL.** Each circuit's synthesized rotations share a total error of 10⁻². The per-rotation tolerance is √(10⁻²/N_rot), priced at 1.15 log₂(1/ε) + 9.2 T.
- **Toffoli** = 7 T.
- **Logical qubits** count the algorithmic register only, not the magic-state factories.
- **Fault budget:** 0.1 expected faults per circuit, counted on T gates.
- **Wall-time model:** one machine, serial, 1 µs per T plus 0.1 ms per circuit, backed by about ten factories. F* = N_T/D_T ≥ 10 is the check that the baseline can be fed.
