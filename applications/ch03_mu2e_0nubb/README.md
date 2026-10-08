# Ch. 3: Mu2e and neutrinoless double-beta-decay nuclear matrix elements

Model: [`estimates/ch03_mu2e_0nubb.py`](../../estimates/ch03_mu2e_0nubb.py). Tests: [`estimates/tests/test_ch03_mu2e_0nubb.py`](../../estimates/tests/test_ch03_mu2e_0nubb.py). Paper chapter: `applications/app11_mu2e_0nubb.tex` (not in this repository).

Every number below is printed by the model, read from `resources.json`, or asserted in the test file; the text says which. Numbers from the model are quoted to two or three figures; the code carries them unrounded.

## 1. What is estimated

The chapter asks what a fault-tolerant machine costs to compute nuclear matrix elements of beyond-Standard-Model operators. Two families are covered: the muon-to-electron conversion operators that Mu2e constrains on 27Al, and the 0νββ transition operators between parent and daughter ground states (48Ca/48Ti, 76Ge/76Se, 82Se/82Kr, 130Te/130Xe, 136Xe/136Ba, 100Mo/100Ru). Nuclei live in a harmonic-oscillator or valence shell-model single-particle basis of N_orb spatial orbitals, mapped to 4 N_orb qubits by Jordan-Wigner (spin times isospin). The ground state is prepared by a Trotterized projection whose length is set by the gap Δ. The operator is block-encoded as a sum of Pauli strings (an LCU) and read out by a Hadamard test, by projection onto the daughter, or by direct measurement.

The model has three eras (`model(a, era)`, lines 1996-2004):

- **2028 benchmark** (`_model_2028`): the 4He dipole operator, the proton charge monopole Σ_i ½(1+τ3,i) j0(q r_i) at q = m_μ, in the 10-orbital 0s-0p-1s0d space (40 qubits). The ground state is not prepared: it is loaded from D classically computed configurations, and the load is counted. The observable is ⟨Ψ0|O|Ψ0⟩ by a Hadamard test on one controlled query of the block encoding. Result: 120-170 LQ, 1.7e4-1.0e5 T per shot over D = 100-600.
- **2033 target** (`_model_2033`): the 27Al μ→e element in the sd valence space (6 orbitals, 24 qubits), state preparation by projection plus one controlled operator query per shot. Result: 74-174 LQ, 1.05e8-4.40e8 T per shot. The 2033 era also carries the 48Ca/48Ti pf pair (1.75-7.30x the 1e9 budget, a near miss) and a 48Ca 0νββ validation ladder in truncated pf spaces.
- **Co-design** (`_model_codesign`): the converged 27Al calculation at e_max = 8-10 (165-286 orbitals, 660-1144 system qubits, 8.5e13-3.3e15 T per preparation), plus the jj44, jj55 and 100Mo/100Ru 0νββ pairs at 2.6-50x the 2033 budget.

## 2. How to run

From the repository root:

```
python -m estimates --chapter ch03                    # assumptions, resources.json rows, per-era exports
python -m estimates --chapter ch03 --intermediates    # also every intermediate the model records
python -m pytest estimates/tests/test_ch03_mu2e_0nubb.py -q
python applications/ch03_mu2e_0nubb/dump_assumptions.py   # the table in section 3 (--raw: notes verbatim)
```

The model needs numpy and scipy (oscillator integrals in `he4_dipole` and `al27_dipole_sd`) and sympy (Clebsch-Gordan coefficients in `_cg`, used by the 2033 operator string counts). The `--chapter ch03` report has three blocks. ASSUMPTIONS lists the 99 inputs with value, provenance and source. INSTANCES lists the seven rows the chapter writes to `resources.json`. EXPORTS gives, per era, the headline LQ and T against the published box (`[ok]` when inside its tolerance), shots, wall time, ε_l, T-depth, F* (= N_T/D_T, the number of factories the circuit can keep busy), the factory-limited floor (EXPORTS' `f_star` pairs the low-T end with the high-depth end and vice versa, so it is a wider band than the per-end F* quoted below, which the model keeps as `f_star_per_end`: 7.1-14.3 in 2028, 11.8 in 2033), and the primitive breakdown with its T sum. The run takes a few seconds.

The test file runs 58 tests and skips 15 in this repository. The skips are the parts of tests that compare against the chapter's `.tex`, which is not shipped here; numeric assertions placed before the skip still run. Next to the paper source all 73 pass.

## 3. Assumptions

`Assumptions` (lines 233-546) is one frozen dataclass of 99 `Tagged` fields:

- 23 Cited: a published source gives the value, and the source column is its bibliography key.
- 33 Assumed: a working assumption the chapter declares.
- 43 Stated-not-derived: the chapter prints the value without derivation. The source is `app11:NN`, a reference into the chapter source; the tests match the printed strings by content, not by line.

Fields marked COMPUTED in the note are Assumed because they come from exact diagonalizations run outside this repository; the source names the output file. `insertion_t_per_toffoli` (7 T per Toffoli) is the report-wide convention stated in the top-level README; its source string names an entry in the paper's convention log. Neither that log nor the diagonalization output files named by the COMPUTED fields are shipped in this repository. `__post_init__` (lines 527-546) rejects a non-positive value in any field, a range in any of the 69 fields listed in `_SCALAR_FIELDS`, a mis-ordered gap band, a Trotter step longer than the shortest τ_gs, an LCU of fewer than two terms, and a fault budget outside (0, 1].

The table below is the output of `dump_assumptions.py`, which reads the live dataclass. Field, value, provenance and source are printed exactly as the code holds them. The note column drops internal decision labels and previous values from the code's notes; `--raw` prints the notes verbatim. The script exits with an error if a field is missing from its grouping, so the table cannot fall behind the code.

#### DOE reference budgets

| field | value | provenance | source | note |
|---|---|---|---|---|
| `rfi_lq_2028` | (150, 250) | Cited | DOE_RFI_2026 | 2028 register envelope |
| `rfi_t_2028` | 100000.0 | Cited | DOE_RFI_2026 | 2028 hard-op cap |
| `rfi_lq_2033` | 1000 | Cited | DOE_RFI_2026 | 2033 '1000+ LQ' |
| `rfi_t_2033` | 1000000000.0 | Cited | DOE_RFI_2026 | 2033 hard-op envelope |

#### State preparation by projection, Eq. (Ngate11): 2033 and co-design cells

| field | value | provenance | source | note |
|---|---|---|---|---|
| `gap_mev` | (1.0, 2.0) | Stated-not-derived | app11:87 | 'nuclear gap scales Delta ~ 1-2 MeV' |
| `gap_large_mev` | 4.0 | Stated-not-derived | app11:150 | 'an optimistic Delta = 4 MeV gap' |
| `gap_small_mev` | 0.5 | Stated-not-derived | app11:156 | 'at a Delta = 0.5 MeV small-gap end, the regime Mo's shape-coexistence physics argues for' |
| `dt_gev_inv` | 0.5 | Stated-not-derived | app11:88 | 'at the dt ~ 0.5 GeV^-1 step' |
| `n_trotter_quoted` | (3100.0, 13000.0) | Stated-not-derived | app11:88 | 'N_Trotter = 1000pi-4000pi ~ 3.1e3-1.3e4' as printed; every product below carries the exact 1000pi-4000pi |
| `c_proj` | (3.141592653589793, 6.283185307179586) | Assumed | app11:87 | 'c_proj in [pi, 2pi] ... Heisenberg-limited phase estimation (pi) and the standard projection prescription (2pi)' |
| `ch2_response_horizon_gev_inv` | 100.0 | Stated-not-derived | app11:87 | implied by 'i.e. 5-10x the response horizon of Ch. 2' |
| `eps_synth_total` | 0.01 | Stated-not-derived | app11:90 | 'a total synthesis error of 1e-2 per shot' |
| `prep_lever` | 10 | Stated-not-derived | app11:163 | 'improves on the tau_gs ~ 1/Delta projection cost by a factor of ten' |

#### Registers, Eq. (Nq11)

| field | value | provenance | source | note |
|---|---|---|---|---|
| `qubits_per_orb` | 4 | Stated-not-derived | app11:79 | '2 . 4 N_orb + N_anc (JW second quantization)' |
| `n_orb_he4` | 10 | Stated-not-derived | app11:122 | 'N_orb ~ 10 active space of Ch. 2' |
| `n_orb_sd` | 6 | Stated-not-derived | app11:90 | '27Al has N_reg=1, N_orb=6'; 24 system qubits |
| `n_orb_pf` | 10 | Stated-not-derived | app11:131 | '~80 for the 48Ca/48Ti pair (pf shell)' = 2.4.10 |
| `n_orb_jj44` | 11 | Stated-not-derived | app11:156 | 'jj44 x jj55 = 11 + 16 = 27 orbitals' |
| `n_orb_jj55` | 16 | Stated-not-derived | app11:156 | 'jj44 x jj55 = 11 + 16 = 27 orbitals' |
| `emax_first` | (3, 4) | Stated-not-derived | app11:152 | 'The first enlargement, e_max=3' ... 'At e_max=4' |
| `emax_converged_al` | (8, 10) | Stated-not-derived | app11:135 | 'An e_max = 8 (10) harmonic-oscillator basis' |
| `emax_converged_ge` | (10, 14) | Stated-not-derived | app11:137 | 'Spaces at e_max = 10-14' |
| `mass_number_ge` | 76 | Cited | LEGEND1000 | A = 76 for the first-quantized line of eq:Nq11 |
| `ho_orbitals_quoted` | (165, 286) | Stated-not-derived | app11:135 | 'carries 165 (286) spatial orbitals' |
| `n_anc` | (50, 150) | Stated-not-derived | app11:85 | 'The ancilla budget N_anc ~ 50-150' |
| `n_anc_2028` | (80, 130) | Stated-not-derived | app11:126 | 'plus ~80-130 ancilla for the BE workspace' |
| `qrom_and_deficit` | 1 | Cited | Babbush_PRX_2018 | 'only log L ancillae' (main_draft.tex:636) for controlled unary iteration; uncontrolled: log L - 1 |
| `test_qubits` | 1 | Stated-not-derived | app11:127 | '+ 1 test qubit at the insertion stage' |

#### Precision, shots, fault budget and wall-time conventions

| field | value | provenance | source | note |
|---|---|---|---|---|
| `eps_2028` | (0.1, 0.3) | Stated-not-derived | app11:122 | 'With eps ~ 0.1-0.3 on a single matrix element' |
| `eps_mu2e` | 0.05 | Stated-not-derived | app11:92 | '~10% on CR means eps ~ 0.05 on the matrix element' |
| `eps_0nubb` | 0.1 | Stated-not-derived | app11:92 | 'Beating the current factor-2-3 spread needs eps ~ 0.1' |
| `n_op_mu2e` | 6 | Cited | Cirigliano_Kitano_Okada_Tuzon | six operator categories |
| `n_op_0nubb` | 5 | Stated-not-derived | app11:233 | '5 0nubb operators' |
| `n_pairs_cross_check` | 3 | Stated-not-derived | app11:233 | 'x 3 pairs' |
| `n_basis_scans` | (2, 3) | Stated-not-derived | app11:233 | 'x 2-3 basis-cutoff scans' |
| `n_survey_elements` | 20 | Stated-not-derived | app11:225 | '~20 operator-isotope elements of the 2033 survey'; 6 + 5x3 = 21 |
| `t_gate_s` | 1e-06 | Assumed | app11:127,185 | 1 us per T-gate (hard operation), the report-wide wall-time convention |
| `shot_overhead_s` | 0.0001 | Assumed | Google_QEC_below_threshold | per-shot init + readout + decode, 0.1 ms, the report convention |
| `campaign_horizon_yr` | 5 | Stated-not-derived | app11:233 | 'Campaign horizon & 5 years' (requirements table); the 2033 wall is printed as a share of it, app11:149 |
| `fault_budget` | 0.1 | Stated-not-derived | app11:93 | 'the quoted eps_l is 0.1/N_gate^shot: ... 0.1 faults are expected per shot' |

#### 2028 benchmark: operator insertion, state load and amplitude estimation

| field | value | provenance | source | note |
|---|---|---|---|---|
| `t_be_topdown` | (1000.0, 10000.0) | Stated-not-derived | app11:122 | 'T_shot^2028 ~ T_BE ~ 1e3-1e4' |
| `be_queries_2028` | 1 | Assumed | app11:122,124 | derived here: Hadamard test on one controlled BE query gives <Psi_0\|O\|Psi_0>/lambda exactly; O Hermitian, so no real-part step |
| `gamma_trunc` | 128 | Stated-not-derived | app11:124 | 'a Gamma <= 128-term insertion' |
| `gamma_256` | 256 | Stated-not-derived | app11:124 | '0.11x at Gamma = 256' |
| `gamma_512` | 512 | Stated-not-derived | app11:124 | '0.23x at Gamma = 512' |
| `select_and_deficit` | 1 | Cited | Babbush_PRX_2018 | 'always ends up with L-1 AND computations' (arXiv:1805.03662 main_draft.tex:716); the L-1 count includes the control input |
| `insertion_t_per_toffoli` | 7 | Assumed | TRACKED_CHANGES.md R5 | 7 T per Toffoli, the report-wide convention (the source's temporary-AND count would give 4) |
| `qrom_t_per_amplitude` | 7 | Cited | Babbush_PRX_2018 | 'D-1 Toffolis per QROM pass' from the source, at the report's 7 T per Toffoli; used only in the lookup-only comparison |
| `d_amplitudes` | 600 | Stated-not-derived | app11:126 | 'the total stays under the cap for D <~ 600' |
| `d_plausible` | (100.0, 1000.0) | Stated-not-derived | app11:124 | lower end of the priced D band (1e2); D = 1e3 is kept only as a comparison point (composed_fraction_at_d_hi) |
| `he4_d_fid99` | (15, 17) | Assumed | scratchpad chopen/app11/he4D_u1.0_s1.0_b{10.0,2.0}.json | COMPUTED: D for fidelity 0.99 (beta = 10, 2) |
| `he4_d_fid999` | (77, 103) | Assumed | scratchpad chopen/app11/he4D_u1.0_s1.0_b{2.0,10.0}.json | COMPUTED: D for fidelity 0.999 (beta = 2, 10) |
| `he4_d_dipole_1pct` | 1 | Assumed | scratchpad chopen/app11/he4D_u1.0_s1.0_b{2.0,10.0}.json | COMPUTED: the single largest configuration already gives the dipole <O> within 1% (0.1% needs D = 11), both beta |
| `qrom_passes` | (2, 3) | Assumed | app11:126 | lookup-only accounting (load, uncompute, optional verification pass), kept for comparison; the priced shot makes one complete load (loads_per_shot) |
| `loads_per_shot` | 1 | Assumed | app11:126 | derived here: one complete load per shot; it erases its own index and the test measures, so nothing is uncomputed |
| `load_angle_error` | 0.01 | Assumed | arxiv_1709_06648 | angle-rounding error of the load, 2 pi L / 2^b (Low, Kliuchnikov and Schaeffer); set equal to the per-shot synthesis budget |
| `mlae_estimator_constant` | 1.0 | Cited | arxiv_1904_10246 | Cramer-Rao bound with every circuit at m_max (Fisher_final); the LIS schedule k = 0..K is a sensitivity (mlae_lis_*) |

#### 2028 benchmark: the 4He dipole operator

| field | value | provenance | source | note |
|---|---|---|---|---|
| `he4_charge_radius_fm` | 1.67824 | Cited | Krauth_muonicHe_2021 | muonic-He Lamb shift, r_alpha = 1.67824(83) fm (Nature 589, 527) |
| `proton_radius_fm` | 0.8409 | Cited | CODATA2018 | r_p = 0.8409(4) fm |
| `hbarc_mev_fm` | 197.3269804 | Cited | CODATA2018 | hbar c |
| `m_mu_mev` | 105.6583755 | Cited | CODATA2018 | muon mass; q = m_mu for the dipole channel |
| `m_nucleon_mev` | 938.9187 | Cited | CODATA2018 | (m_p + m_n)/2, only for printing hbar omega |
| `z_he4` | 2 | Stated-not-derived | app11:122 | 4He: Z = 2 protons (A = 4) |
| `mass_number_he4` | 4 | Stated-not-derived | app11:122 | 4He |

#### 2033 deliverable: 27Al dipole factor, readout and T-depth

| field | value | provenance | source | note |
|---|---|---|---|---|
| `al27_charge_radius_fm` | 3.061 | Cited | Angeli_Marinova_2013 | 27Al rms charge radius 3.0610(31) fm (ADNDT 99, 69 (2013)) |
| `z_al27` | 13 | Stated-not-derived | app11:6 | 27Al: Z = 13 |
| `mass_number_al27` | 27 | Stated-not-derived | app11:6 | 27Al |
| `z_core_al27` | 8 | Stated-not-derived | app11:44 | sd valence space: 16O core, 8 core protons (0s^2 0p^6) |
| `trotter_concurrency` | 12 | Assumed | arxiv_1905_09749 | commuting Trotter rotations in flight at once (12-24 of the 50-150 ancillas); G >= 10 keeps the 10-factory baseline fed |
| `sd_spin_expectation` | (0.3, 0.34) | Assumed | app11:151 | <S_p>(27Al) in sd, an estimate pending a shell-model value |
| `direct_readout_variance` | 1.0 | Assumed | app11:151 | per-shot variance of the direct (basis-rotated Z) readout of the SD one-body operator, ~1 (estimate); worst case x8.3 (Popoviciu bound) |
| `projection_success` | 0.5 | Assumed | app11:151 | heralded ground-state projection success p for deformed 27Al (estimate; computable exactly in sd); sensitivity band 0.1-0.9 |
| `projection_success_band` | (0.1, 0.9) | Assumed | app11:151 | sensitivity band on p |
| `first_result_scans` | 1 | Stated-not-derived | app11:153 | the first result uses one Hamiltonian (one basis-cutoff scan) |
| `eps_first_generic` | 0.3 | Assumed | app11:151 | report-wide first-result statistical error, quoted for comparison only |
| `gap_al27_usdb_mev` | 2.327866 | Assumed | scratchpad chopen/app11/gap_usdb_sd_27Al_J-1.0.json | COMPUTED: first excited state reachable from an M = 5/2 reference (J >= 5/2; the 7/2+), USDB in sd |
| `gap_al27_usdb_same_j_mev` | 2.698773 | Assumed | scratchpad chopen/app11/gap_usdb_sd_27Al_J2.5.json | COMPUTED: 5/2+_2 - 5/2+_1, USDB in sd (comparison only) |

#### 2033 0nubb validation stage: 48Ca/48Ti truncation ladder and resonant-drive readout

| field | value | provenance | source | note |
|---|---|---|---|---|
| `gap_0plus_f7p3_ca48_mev` | 6.100132 | Assumed | scratchpad chopen/app11/gap_kb3g_f7p3_48Ca_J0.0.json | COMPUTED: 48Ca 0+_2 - 0+_1, KB3G in f7/2 p3/2 (M-dim 57) |
| `gap_0plus_f7p3_ti48_mev` | 3.68472 | Assumed | scratchpad chopen/app11/gap_kb3g_f7p3_48Ti_J0.0.json | COMPUTED: 48Ti 0+_2 - 0+_1, KB3G in f7/2 p3/2 (M-dim 5296) |
| `gap_0plus_f7p3p1_ca48_mev` | 5.500798 | Assumed | scratchpad chopen/app11/gap_kb3g_f7p3p1_48Ca_J0.0.json | COMPUTED: 48Ca, KB3G in f7/2 p3/2 p1/2 (M-dim 325) |
| `gap_0plus_f7p3p1_ti48_mev` | 3.603676 | Assumed | scratchpad chopen/app11/gap_kb3g_f7p3p1_48Ti_J0.0.json | COMPUTED: 48Ti, KB3G in f7/2 p3/2 p1/2 (M-dim 24453) |
| `gap_0plus_pf_ca48_mev` | 5.174846 | Assumed | scratchpad chopen/app11/gap_kb3g_pf_48Ca_J0.0.json | COMPUTED: 48Ca, KB3G in full pf (M-dim 12022) |
| `gap_0plus_pf_ti48_mev` | 4.264629 | Assumed | scratchpad chopen/app11/gap_kb3g_pf_48Ti_J0.0.json | COMPUTED: 48Ti, KB3G in full pf (M-dim 634744) |
| `gap_0plus_ensdf_ca48_mev` | 4.283 | Cited | ENSDF | 48Ca first excited 0+ (4283 keV); comparison only, the priced gaps are KB3G |
| `gap_0plus_ensdf_ti48_mev` | 2.997 | Cited | ENSDF | 48Ti first excited 0+ (2997 keV); comparison only, the priced gaps are KB3G |
| `eps_simplified` | 0.3 | Assumed | app11:161 | 30% statistical error on \|M\| for the 0nubb first light |
| `lam_over_m_f7p3` | 69.1 | Assumed | app11:161 | lambda/\|M\| in f7/2 p3/2 from a schematic interaction; not derived here |
| `p_ref_f7p3` | 0.93 | Assumed | app11:161 | closed-f7/2 reference overlap in f7/2 p3/2 (schematic interaction) |
| `lam_over_m_f7p3p1` | 164.0 | Assumed | app11:161,163 | lambda/\|M\| in f7/2 p3/2 p1/2 (schematic interaction); enters the qDRIFT count |
| `lam_over_m_pf` | (600.0, 720.0) | Assumed | app11:161 | lambda/\|M\| in full pf: Pauli-LCU lambda 208.4 over \|M\| ~ 0.29-0.35 (schematic interaction) |
| `rabi_drive_mev` | (0.3, 0.5) | Assumed | app11:163 | theta\|M\|, the Rabi frequency of the resonant drive; kept <= 0.5 MeV for a <~6% Stark bias (schematic interaction) |
| `qdrift_eps` | 0.1 | Assumed | Campbell_qDRIFT_2019 | qDRIFT error on the drive term; N_q = 2 (lambda/M)^2 (pi/2)^2 / eps_q |
| `rabi_shots` | (250, 1000) | Assumed | app11:163 | Rabi-readout shots including a 5-10-detuning resonance scan; 250 give ~5% on \|M\| at t_pi/2 (schematic interaction) |

#### External comparison points and the utility box (not used in any T-count)

| field | value | provenance | source | note |
|---|---|---|---|---|
| `fq_light_nucleus_t` | 10000000.0 | Cited | arxiv_2507_22814 | 'sits at the ~1e7-T scale' |
| `ext_qubitization_sd_toffoli` | 1000000000.0 | Cited | arxiv_2607_21563 | '~1e9 Toffoli gates in the sd shell' |
| `ext_qubitization_pb208_toffoli` | 30000000000.0 | Cited | arxiv_2607_21563 | '~3e10 ... 102-orbital 208Pb-core space' |
| `mu2e_tpc_musd` | 315.7 | Cited | DOE_HEP_FY25_CJ | 'the $315.7M Mu2e total project cost'; FY 2025 CJ, Science/HEP p. 245: BCP approved 2022-12-21, TPC $315,700,000 |
| `utility_fraction` | 0.6 | Stated-not-derived | app11:188 | 'about 60% of the $315.7M Mu2e total project cost ... The fraction is not derived.' The dollar value is linear in it |
| `nme_spread` | 3 | Cited | arxiv_2308_15634 | 'factor-three spread across nuclear models' |

99 fields: Cited 23, Stated-not-derived 43, Assumed 33

## 4. How the numbers are built

Two conventions from `estimates/common.py` run through everything. A synthesized rotation costs c_T = 1.15 log2(1/ε_rot) + 9.2 T (`t_per_rotation`, line 125). Every circuit sets its own tolerance from the total synthesis budget of 1e-2 per shot: ε_rot = √(1e-2 / N_rot), N_rot the number of rotations in one shot of that circuit (`eps_rot_for`, line 172; Ch. 3 wraps it as `c_T`, lines 565-568). Each end of each band is therefore its own circuit with its own c_T. A Toffoli costs 7 T. The required logical error rate is ε_l = 0.1 / N_T (`eps_l`, lines 630-632).

### 4.1 Registers (Eq. Nq11)

`register` (lines 597-603) returns 4 N_orb system qubits for every cell. Parent and daughter of a 0νββ pair share one Fock register, since the operator only changes how many protons and neutrons fill the same orbitals. LQ = 4 N_orb + N_anc with N_anc = 50-150 (2033) or 80-130 (2028). The cells: 4He 40, 27Al sd 24, pf 40, jj44 44, jj55 64. 100Mo/100Ru puts protons in jj44 and neutrons in jj55, 22 + 32 = 54 modes, so N_orb = 13.5 (`n_orb_mo_ru`, lines 606-610). The oscillator ladder has (e+1)(e+2)(e+3)/6 spatial orbitals (`ho_spatial_orbitals`, lines 613-615): 20, 35, 165, 286 at e_max = 3, 4, 8, 10. The co-design flagship LQ is the system register alone, 660-1144 (710-1294 with ancillae).

The 2028 register is counted stage by stage (`_model_2028`, lines 1494-1502; pinned by `test_2028_register`). The insertion stage holds 40 system + 7 LCU address + 7 unary-iteration AND ancillae + 1 test qubit = 55. The load needs 5L - 3 = 47 ancillae at L = ⌈log2 600⌉ = 10 (`sos_load_ancillae`, lines 965-968), 87 in all. While the amplitude angles are applied it holds 40 + 10 + (3b - 1) = 88 with b = 13, the peak. Under amplitude estimation the reflection about |0⟩ acts on 40 + 10 + 7 + 1 = 58 qubits with 56 AND ancillae, 114, plus the 13-qubit phase-gradient state: 127. All sit inside the 120-170 LQ budgeted register and the 150-250 LQ reference.

### 4.2 2028 benchmark: T per shot

The shot is one controlled block-encoding query plus one complete state load (`_model_2028`, lines 1424-1714).

**Insertion** (`insertion_cost`, lines 679-718). The LCU is truncated at Γ = 128 strings, so the address has k = ⌈log2 Γ⌉ = 7 qubits. PREPARE is a rotation tree of 2^k - 1 = 127 R_y; PREPARE and PREPARE† give 254 rotations. SELECT is unary iteration controlled on the test qubit, Γ - 1 = 127 Toffolis (Babbush et al. 2018); the Pauli targets are Clifford. A Hadamard test on the controlled query returns ⟨O⟩/λ exactly, so there is one query and no signal-processing phases (`be_queries_2028` = 1). With N_rot = 254, ε_rot = 6.27e-3 and c_T = 17.61 (the insertion is synthesized at its own 254 rotations; the 10 load rotations below use the shot's 264, c_T = 17.65, a difference of about 10 T), so T_ins = 254 × 17.61 + 127 × 7 = 4474 + 889 = 5363 T, 0.054 of the 1e5 cap. The same count at Γ = 256 and 512 gives 0.11 and 0.23 of the cap; the full one-body operator on 40 modes (Γ = 1600 real, 3160 complex) gives 0.93-1.9. The dipole's own 24 strings would cost 1181 T; the 128-string cap leaves room for the other five operator categories, whose string counts are not derived.

**Load** (`sos_load_toffoli`, lines 953-962). The complete sparse-state load of Fomichev et al. (arXiv:2310.18410) writes D configurations and their amplitudes and erases the index in (2L - 2)D + 2^(L+1) + D Toffolis, with L = ⌈log2 D⌉. With one load per shot that is 1556 Toffolis (10892 T) at D = 100 and 13448 (94136 T) at D = 600. The source leaves the amplitude-angle rotations of its 2^(L+1) step unpriced; `load_angle_cost` (lines 971-995) adds them following Low, Kliuchnikov and Schaeffer. It uses L phase-gradient additions of b bits, b - 1 Toffolis each, with b = ⌈log2(2πL / 0.01)⌉ = 13, plus b - 3 = 10 synthesized rotations and one T for the phase-gradient state. At D = 100 (L = 7) that is 84 Toffolis + 10 rotations at c_T = 17.65 + 1 T = 765 T; at D = 600 (L = 10) it is 120 Toffolis and the same 10 rotations, 1017 T.

**Total.** hard_ops = 5363 + 10892 + 765 = 17020 T at D = 100 and 5363 + 94136 + 1017 = 100516 T at D = 600, 1.005x the cap; the largest D that fits is 596 (`d_break`). The primitive breakdown at D = 100, as `--chapter ch03` prints it:

| primitive | count | T each | T |
|---|---|---|---|
| `sos_load_complete` | 1556 | 7 | 10892 |
| `load_angle_cadd` | 84 | 7 | 588 |
| `load_fourier_state` | 10 | 17.65 | 176.5 |
| `load_fourier_state_t` | 1 | 1 | 1 |
| `prepare_rotation_tree` | 254 | 17.61 | 4474 |
| `select_unary_iteration` | 127 | 7 | 889 |
| `pauli_string_select_targets` | 128 | 0 | 0 |
| `hadamard_test_frame` | 1 | 0 | 0 |
| **sum** | | | **17020** |

D is the one physics input that moves this headline. The model carries an exact diagonalization of a schematic stand-in (Minnesota NN interaction plus a Lawson center-of-mass term, in the same 40 modes). It gives fidelity 0.99 at D = 15-17 and 0.999 at D = 77-103, and the largest single configuration already gives the dipole to 1% (`he4_d_*`, `d_he4_computed`). These sit at or below the priced band, so the band is conservative for that Hamiltonian; a chiral ground state is not checked.

**Shots** (`he4_dipole`, lines 814-870). The operator matrix o_pq on the 40 modes uses oscillator radial integrals of j0(qr) at q = m_μ. The oscillator length b = 1.369 fm is fit to the 4He charge radius after removing the proton size and the center-of-mass term (`he4_oscillator_length`, lines 731-739). Its Jordan-Wigner form is the identity plus 24 strings (`jw_one_body_paulis`, lines 792-811), λ = 15.39. The identity coefficient c_0 = Tr(o)/2 is a known constant, so the encoded operator is O - c_0 with λ' = 7.79. For any two-proton state ⟨O⟩ lies in 1.361-1.844 (Z times the extreme eigenvalues of the proton block), so the shot factor (λ'/|⟨O⟩|)² lies in 17.9-32.8; the model takes the upper end, 32.76. At ε = 0.1 that is 100 × 32.76 = 3276 shots per operator and 19656 over the six Mu2e operator categories, each assumed to carry the dipole's factor.

**Depth-capped amplitude estimation** (`mlae_2028`, lines 877-950). Each application of A is one load with its angle additions plus one query; the reflection S_0 on 58 qubits costs 56 Toffolis. The circuit of m = 2k+1 applications costs T(m) = m T_A + (m-1)/2 × 392 + (phase-gradient state, once), with the query rotations synthesized to the deepest circuit's budget, m × 254. The model solves for the largest odd m whose circuit fits under 1e5 at its own budget: m_max = 5 up to D = 128, 3 up to D = 221, 1 beyond. At D = 100 the deepest circuit is 86893 T, and the Cramér-Rao bound with every circuit at m_max gives 3276 / 25 = 131 circuits per operator, 786 in all. The state-independent band on ⟨O⟩ keeps m_max θ inside one monotone branch of sin²(mθ), 0.035 rad from its edge, so no shallow disambiguation circuits are needed.

**Wall time** (`depth_2028_hadamard`, lines 1265-1281; `wall_one_machine`, lines 1320-1325). Per shot the time is max(N_T × 1 μs, D_T × 10 μs) + 0.1 ms: the T supply of the ten-factory baseline, or the T-depth at one reaction time per layer, whichever binds. The compiled T-depth runs PREPARE alongside the load, pipelines the rotation tree to 2^(k-1) c_T = 1127 layers, and erases the load index with an AND tree: D_T = 2382 at D = 100 and 7040 at D = 600. F* = N_T / D_T is 7.1 at D = 100 (below 10, so depth binds: 0.0239 s per shot) and 14.3 at D = 600 (T supply binds: 0.1006 s). The Hadamard-test run is 19656 shots, 470-1978 s (8-33 min) on one machine. The MLAE route at D = 100 has D_T = 7942, F* = 10.9, 0.087 s per circuit, 68 s in all. ε_l = 0.1 / N_T = 9.9e-7 to 5.9e-6.

### 4.3 2033 and co-design: state preparation (Eq. Ngate11)

`prep_t` (lines 589-594) is the chapter's formula for every survey cell:

N_T(shot) = N_rot × c_T(N_rot), N_rot = N_proj × N_Trotter × N_orb^4, N_Trotter = c_proj τ_gs / δt, τ_gs = 1/Δ.

With Δ = 1-2 MeV (τ_gs = 500-1000 GeV^-1), c_proj = π to 2π and δt = 0.5 GeV^-1, N_Trotter = 1000π to 4000π (3142-12566), the low end pairing c_proj = π with Δ = 2 MeV (`n_trotter`, lines 558-562). Each Trotter step applies the dense N_orb^4 two-body Pauli strings as rotations. N_proj = 1 for 27Al (one preparation) and 2 for a 0νββ pair (parent preparation plus projection onto the daughter, priced at the parent's gap). `survey` (lines 1195-1248) evaluates every cell at the generic band, at Δ = 4 MeV, at Δ = 0.5 MeV, and with the preparation depth divided by 10 (`prep_lever`):

| cell | N_proj | N_orb | system qubits | N_rot | c_T | T per shot | x 1e9 |
|---|---|---|---|---|---|---|---|
| 27Al sd | 1 | 6 | 24 | 4.07e6-1.63e7 | 25.6-26.8 | 1.04e8-4.36e8 | 0.10-0.44 |
| 48Ca/Ti pf | 2 | 10 | 40 | 6.28e7-2.51e8 | 27.9-29.1 | 1.75e9-7.30e9 | 1.75-7.30 |
| jj44 pairs | 2 | 11 | 44 | 9.20e7-3.68e8 | 28.2-29.4 | 2.60e9-1.08e10 | 2.60-10.8 |
| 100Mo/Ru | 2 | 13.5 | 54 | 2.09e8-8.35e8 | 28.9-30.1 | 6.03e9-2.51e10 | 6.03-25.1 |
| jj55 pairs | 2 | 16 | 64 | 4.12e8-1.65e9 | 29.5-30.6 | 1.21e10-5.04e10 | 12.1-50.4 |
| 27Al e_max = 3 | 1 | 20 | 80 | 5.03e8-2.01e9 | 29.6-30.8 | 1.49e10-6.19e10 | 14.9-61.9 |
| 27Al e_max = 4 | 1 | 35 | 140 | 4.71e9-1.89e10 | 31.5-32.6 | 1.48e11-6.16e11 | 148-616 |
| 27Al e_max = 8 | 1 | 165 | 660 | 2.33e12-9.31e12 | 36.6-37.8 | 8.53e13-3.52e14 | 8.5e4-3.5e5 |
| 27Al e_max = 10 | 1 | 286 | 1144 | 2.10e13-8.41e13 | 38.5-39.6 | 8.09e14-3.33e15 | 8.1e5-3.3e6 |

The co-design flagship takes the low end of e_max = 8 and the high end of e_max = 10: 8.53e13-3.33e15 T, ε_l = 3.0e-17 to 1.2e-15. Gap variants (`multiple_*_large_gap`, `multiple_mo_ru_small_gap`): jj55 5.95-12.1 at Δ = 4 MeV; 100Mo/Ru 2.96-6.03 at 4 MeV and 25.1-51.1 at 0.5 MeV. With the preparation depth cut tenfold, pf falls to 0.16-0.68 of the budget, jj44 to 0.24-1.01, 100Mo/Ru to 0.28-4.8 across the gap variants (0.56-2.35 at the generic band).

### 4.4 2033: the operator insertion

The insertion is one controlled query of a Γ-string LCU inside the preparation circuit, so its rotations share the circuit's tolerance (`insertion_in_circuit`, lines 1186-1192). Γ is counted by explicit Jordan-Wigner algebra. For 0νββ (`gamma_0nubb`, lines 1112-1139) every piece of the operator is a scalar, so a proton pair and a neutron pair connect only if they share 2M, parity and an allowed pair J. Each surviving term a†a†aa on four distinct modes is 16 strings. In jj44 this keeps 2501 of 53361 terms, Γ = 40016 (853776 dense), and costs 3.98e6-4.13e6 T (6.52e7-6.76e7 dense), 0.04-0.15% of the pair's preparation. For 27Al (`gamma_mu2e_two_body`, lines 1143-1183) the operator is a general charge-conserving scalar one- plus two-body operator on the 24 sd modes, the conservative cover for the six Mu2e categories with two-body currents: Γ = 12900 with the selection rules, 49140 dense. The 2033 headline uses the selection-rule count at its low end and the dense count at its high end: 9.31e5-3.86e6 T.

The 27Al shot is then the whole circuit at the shared tolerance: N_T = (N_rot,prep + N_rot,ins) × c_T + 7 × (Γ - 1) = 1.054e8-4.404e8 (`_model_2033`, lines 1729-1733), 0.105-0.440 of the 1e9 budget, ε_l = 2.3e-10 to 9.5e-10. At the low end:

| primitive | count | T each | T |
|---|---|---|---|
| `trotter_rotations_prep_al_sd` | 4.072e6 | 25.65 | 1.044e8 |
| `operator_insertion_sd_rotations` | 32766 | 25.65 | 8.41e5 |
| `operator_insertion_sd_select` | 12899 | 7 | 9.03e4 |
| **sum** | | | **1.054e8** |

### 4.5 2033: shots, T-depth and wall time for 27Al

**What is measured** (`al27_dipole_sd`, lines 998-1031; `_model_2033`, lines 1769-1806). The sd dipole is diagonal on the 12 valence proton modes, λ' = 3.37, and any state of five valence protons has ⟨O⟩ in 2.794-2.829, ±0.61% (`si_pinned_fraction`). The bare spin-independent elements are therefore fixed classically and serve as validation. The shots go to the spin-dependent element, ⟨O_SD⟩ = 2⟨S_p⟩ o_0d = 0.335-0.380 for ⟨S_p⟩ = 0.30-0.34 and o_0d = 0.559. A Hadamard test would need (λ'/⟨O_SD⟩)² = 79-101 times ε^-2. The model instead reads it directly after a single-particle basis rotation: var / (ε ⟨O_SD⟩)² = 2769-3557 shots at ε = 0.05 and variance 1, divided by the heralded projection success p = 0.5, giving 5539-7114 shots per Hamiltonian. The first result is one Hamiltonian; the campaign is 2-3 basis-cutoff scans, 11077-21342 shots.

**T-depth** (`depth_2033`, lines 1308-1317). Commuting Trotter rotations share one diagonalizing Clifford, so G = 12 run at once: D_T = N_rot,prep c_T / 12 + (Γ - 1) + 2^(k-1) c_T = 8.93e6-3.73e7, F* = 11.8 at both ends. The ten-factory T supply binds, and a shot takes 105-440 s (1.8-7.3 min).

**Wall time.** First result 5.84e5-3.13e6 s (6.8-36 days); campaign 1.17e6-9.40e6 s (13.5 days to 3.6 months), at most 6.0% of the 5-year horizon. Sensitivities: p over 0.1-0.9 moves the campaign to 7.5 days-1.5 years. Gate by gate (G = 1), no number of factories brings it below 135 days-3.0 years. A 30% target on the spin-dependent element would cut the first result to 154-198 shots, 4.5-24 hours, without separating the operator categories.

### 4.6 2033: the 0νββ validation stage

`simplified_0nubb` (lines 1366-1417) prices 48Ca → 48Ti in truncated pf spaces at the 0+ gaps of the KB3G interaction restricted to the kept orbits, which the model carries as computed inputs (section 3). For a scalar operator each projection only has to resolve the next 0+ state, so the pair cost is Σ_i c_proj (1/Δ_i)/δt Trotter steps over parent and daughter (`_pair_prep_two_gaps`, lines 1328-1342).

- f7/2 p3/2 (24 modes, 74-174 LQ), projection readout: 9.10e7-1.86e8 T per shot; (λ/|M|)²/(4ε²p_ref) = 69.1² / (4 × 0.09 × 0.93) = 14262 shots at 30% on |M|; 15-31 days, F* 11.9.
- f7/2 p3/2 p1/2 (28 modes, 78-178 LQ), resonant-drive readout (`_rabi_cost`, lines 1345-1363): evolve under H + μN_p + θ(O + O†) for a quarter Rabi period and measure proton number. Steps (π/2)/(θ|M|)/δt at N_orb^4 rotations each, plus the qDRIFT-sampled drive, 2(λ/|M|)²(π/2)²/ε_qD = 1.33e6 single-Pauli rotations. 4.40e8-7.20e8 T, 250-1000 shots, 2.0-11 days. The qDRIFT samples run one at a time, so D_T = N_Trotter-rot c_T/12 + N_qD c_T and F* = 6.4-7.7.
- Full pf by the same readout: 2.27e9-3.72e9 T (2.3-3.7x the budget), F* 3.5-3.8, 19-110 days.
- Full pf by projection at its 0+ gaps: 7.35e8-1.50e9 T, but at λ/|M| = 600-720 at least 1.0e6-1.44e6 shots, 23-68 years.

### 4.7 What goes where

`PUBLISHED` (lines 2025-2043) holds the box values with their tolerances; `INSTANCE_ROWS` (lines 2210-2243) writes the `resources.json` rows; `printed` (lines 2046-2128, with its 2033 and 2028 wall-time companion at lines 2131-2207) builds every number the chapter prints from the model, each carried at full precision and cut to two figures only at print; `utility` (lines 2007-2022) is the utility box.

## 5. What the paper prints

`resources.json` rows for Ch. 3 (written by `INSTANCE_ROWS`; they feed the paper's resource-landscape figure and Table 1.1):

| era | label | LQ | T per shot | x budget | role |
|---|---|---|---|---|---|
| 2028 | 4He dipole Hadamard | 120-170 | 17020.3-100516.3 | | conditional on D (100-600) |
| 2033 | 27Al sd (prep+insertion) | 74-174 | 1.0537e8-4.4036e8 | | fits |
| 2033 | 48Ca/Ti (pf) | 90-190 | 1.7540e9-7.3050e9 | 1.75-7.30 | near miss; state prep only |
| codesign | 76Ge/Se, 82Se/Kr (jj44) | 94-194 | 2.5971e9-1.0812e10 | 2.60-10.8 | co-design; state prep only |
| codesign | 130Te/Xe, 136Xe/Ba (jj55) | 114-214 | 1.2137e10-5.0443e10 | 12.1-50.4 | co-design; state prep only |
| codesign | 100Mo/Ru (jj44 x jj55) | 104-204 | 6.0337e9-2.5095e10 | 6.03-25.1 | reach target |
| codesign | flagship 27Al (e_max 8-10) | 660-1144 | 8.5324e13-3.3309e15 | 8.5e4-3.3e6 | system qubits only, per preparation |

Table 1.1's Ch. 3 row quotes the 2028 and 2033 headlines (120-170 LQ, 1.7e4-1.0e5 T; 74-174 LQ, 1.1-4.4e8 T) and the 0νββ pairs at 1.8-11x, the pf low end to the jj44 high end of the rows above. The check of the table text lives with the paper.

The boxes, as `PUBLISHED` holds them and as the model computes them:

| box | printed | model | tests that pin it |
|---|---|---|---|
| 2028 LQ | 120-170 | 120-170 (peak used 88; 127 under MLAE) | `test_published_is_the_box_as_printed`, `test_2028_register` |
| 2028 T per shot | 1.7e4-1.0e5 | 17020-100516 | `test_2028_bottom_up_costing_is_the_box`, `test_2028_oracle_load_is_printed_separately`, `test_complete_load_counts` |
| 2028 shots | 2.0e4 (or 7.9e2 MLAE circuits) | 19656 (786) | `test_2028_shots_and_wall_time`, `test_mlae_circuits_and_wall`, `test_mlae_depth_is_self_consistent` |
| 2028 wall | 8-33 min (or 68 s MLAE) | 470-1978 s (68.4 s) | `test_2028_shots_and_wall_time`, `test_mlae_circuits_and_wall` |
| 2028 time per shot | 0.024-0.10 s | 0.0239-0.1006 s | `test_2028_shots_and_wall_time` |
| 2033 LQ | 74-174 | 74-174 | `test_2033_registers` |
| 2033 T per shot | 1.1e8-4.4e8 | 1.054e8-4.404e8 | `test_2033_al_sd_cell` |
| 2033 ε_l | 2.3-9.5e-10 | 2.27e-10-9.49e-10 | `test_2033_al_sd_cell`, `test_fault_budget_is_applied_everywhere` |
| 2033 shots | 5.5-7.1e3 per Hamiltonian; 1.1-2.1e4 campaign | 5539-7114; 11077-21342 | `test_2033_shots_as_printed` |
| 2033 wall | 6.8-36 days first result; 14 days-3.6 months campaign | 5.84e5-3.13e6 s; 1.17e6-9.40e6 s | `test_2033_tiers_depth_and_walls` |
| 2033 time per shot | 1.8-7.3 min | 105.4-440.4 s | `test_2033_shots_wall_time_campaign` |
| 0νββ near miss (pf) | x1.8-7.3 | 1.75-7.30 | `test_2033_valence_pairs` |
| co-design pairs | x2.6-11 (jj44), x12-50 (jj55), x6-25 (Mo/Ru) | as in section 4.3 | `test_2033_valence_pairs`, `test_2033_jj55_codesign_pairs`, `test_2033_mo_ru_generic_gap` |
| validation ladder | f7p3 0.91-1.9e8 T, 15-31 days; f7p3p1 4.4-7.2e8 T, 2.0-11 days; pf x2.3-3.7, 19-110 days | as in section 4.6 | `test_simplified_0nubb_rows` |
| co-design LQ | ~700-1100 | 660-1144 | `test_codesign_flagship`, `test_codesign_flagship_fits_1000_lq_as_printed` |
| co-design T | 8.5e13-3.3e15 | 8.532e13-3.331e15 | `test_codesign_flagship` |
| utility | ~$32M per instance (~$190M over six) | 31.57 (189.42) $M = 0.6 × 315.7 / 6 | `test_utility_box` |

`test_instance_rows` pins the `resources.json` rows and their roles, and `test_printed_strings_are_in_the_tex` checks every string `printed(a)` returns against a fixed expected list (and, next to the paper, against the chapter text). The load-angle add-on (L, b, 765-1017 T, the D = 596 fit) is pinned by the load-angle test (test file, line 1440).

## 6. Circuit status

From `docs/CIRCUIT_STATUS.md` (generated from the breakdowns by `python -m estimates --status`): 10 COMPILED, 2 SCALING, 0 CONJECTURE, 0 UNSOURCED. COMPILED here means the count follows from an explicit gate-level construction written out in the model; no circuit file is emitted.

| era | primitive | status | source |
|---|---|---|---|
| 2028 | `sos_load_complete` | COMPILED | arXiv:2310.18410 |
| 2028 | `load_angle_cadd` | COMPILED | arXiv:1812.00954 (adder: arXiv:1709.06648) |
| 2028 | `load_fourier_state` | COMPILED | arXiv:1812.00954 |
| 2028 | `load_fourier_state_t` | COMPILED | arXiv:1812.00954 |
| 2028 | `prepare_rotation_tree` | COMPILED | derived in the model |
| 2028 | `select_unary_iteration` | COMPILED | Babbush et al. 2018 |
| 2028 | `pauli_string_select_targets` | COMPILED (Clifford) | derived in the model |
| 2028 | `hadamard_test_frame` | COMPILED (Clifford) | derived in the model |
| 2033 | `trotter_rotations_prep_al_sd` | SCALING | chapter, Eq. Ngate11 |
| 2033 | `operator_insertion_sd_rotations` | COMPILED | derived in the model |
| 2033 | `operator_insertion_sd_select` | COMPILED | Babbush et al. 2018 |
| codesign | `trotter_rotations_prep_al_emax8` | SCALING | chapter, Eq. Ngate11 |

The two SCALING rows carry essentially all of the 2033 and co-design T: the preparation is a counting argument (dense N_orb^4 strings per step, a projection length c_proj/Δ), not a constructed circuit. `test_no_breakdown_primitive_claims_compiled` checks that the state-preparation rows stay SCALING and the insertion rows are the ones marked COMPILED.

## 7. Open items and limitations

From `docs/OPEN_ITEMS.md` and the model's own notes:

- **The 2028 headline depends on D.** It is priced over D = 100-600. A schematic Minnesota ground state reaches fidelity 0.999 at D = 77-103, but a chiral (SRG-evolved) ground state could need more. A narrower band would move the headline toward the insertion alone (5363 T).
- **The other five Mu2e categories** are assumed to carry the 4He dipole's shot factor and to fit in the 128-string insertion; their string counts are not derived.
- **The 2033 readout covers only the spin-dependent element.** If the deliverable is the VS-IMSRG-evolved operators, the readout of their two-body parts is not priced. ⟨S_p⟩ = 0.30-0.34 and the per-shot variance of about one are estimates (the worst-case variance is about 8 times larger); a shell-model value would replace both.
- **Projection success p = 0.5 for 27Al** is assumed. It is computable exactly in sd; p over 0.1-0.9 spans 7.5 days to 1.5 years of campaign.
- **Gaps.** The 27Al and survey cells use a generic 1-2 MeV band, while the 48Ca/48Ti ladder uses the KB3G 0+ gaps of each truncated space. The USDB gap that an M = 5/2 reference can reach in 27Al is 2.33 MeV and would lower the 27Al preparation to 8.9e7-1.8e8 T (`al27_prep_t_at_usdb_gap`); the band is kept. The daughter projection is priced at the parent's gap; a daughter with a smaller gap raises the cost in proportion to 1/Δ.
- **The resonant-drive readout is checked only on a schematic interaction.** λ/|M| = 69, 164 and 600-720 and p_ref = 0.93 come from that interaction; Stark bias, resonance identification and the Trotter and qDRIFT errors are not simulated, so its shot counts and walls are conditional.
- **Depth.** The resonant-drive walls are depth-limited (qDRIFT samples run one at a time, F* 3.5-7.7). The 2033 27Al walls assume 12 commuting rotations in flight; gate by gate the campaign floor is 135 days to 3.0 years.
- **Cancellations not taken.** PREPARE trees that cancel or act on unused amplitudes in the 2028 shot are kept as priced, so the insertion count is an upper bound for its construction. The 2033 27Al shot also keeps the operator query (about 1% of its T), which the direct readout does not use.
- **Category separation** needs targets with different Z, N/Z and spin; a second target such as 48Ti in pf is not priced.
- **Fault model.** ε_l = 0.1/N_T counts T gates only; Clifford, idle, measurement and injection faults are not modeled.
- **Utility.** No Mu2e-specific spread of the nuclear matrix elements is quoted, and the 60% utility fraction is a stated choice; the dollar figure is linear in it.
- **Beyond 2033.** 76Ge (jj44 and converged bases) has no simplification that keeps its physics inside the 2033 budget, and the converged 27Al flagship is 8.5e4-3.3e6 times over it.

## 8. Conventions used

The report-wide conventions are logical-qubit counts, 7 T per Toffoli, the per-circuit rotation-synthesis rule (1.15 log2(1/ε_rot) + 9.2 T at ε_rot = √(1e-2/N_rot)), the fault budget of 0.1 per shot and the single-machine wall-time model (1 μs per T from about ten factories, 10 μs per sequential T layer (`REACTION_TIME_S`, `estimates/common.py` line 294), 0.1 ms per shot, nothing run in parallel across machines). These are shared by every chapter and described in the [top-level README](../../README.md) and `estimates/common.py`. Two Ch. 3 specifics: the 2033 first result is quoted at 5% on the spin-dependent element, because 30% would not separate the operator categories (the 30% figure is computed for comparison), and the 2028 benchmark has a single tier.
