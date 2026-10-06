# Open items

What the chapter models cannot decide on their own. Each line is either a physics judgment the
model carries as an input (tagged `Assumed` or `Stated` in the chapter's `Assumptions`), or a piece
of the circuit or the campaign that is bounded, estimated or left at zero rather than priced. Where
a line names a model field or function (in backticks), the number comes from there; where it names
a sensitivity, the model computes both ends and the chapter prints the one stated in its box. None
of these items changes a printed number until an author decides it; the tests pin the boxes as
printed.

Chapters 11 and 12 have no model in this package and are not listed.

## Ch. 2, neutrino-nucleus response (`ch02_nu_nucleus.py`)

- No T-depth or device schedule is derived for the two benchmarks: the depth, factory and wall-time exports are `None` with `depth_status = "not_established"`; shot counts and gate volumes are reported instead.
- The 2033 primitive counts (`counts("2033", steps)`) are constructive allowances for the 908-qubit construction, not a compiled full-register circuit.
- Rotation synthesis is priced at the leading Ross-Selinger term, 3 log2(N_R / 0.005) T per rotation (`coherent_synthesis_allowance`, `exact_synthesis = False`); the subleading term is pending.
- Physical accuracy of the finite-volume correlators and the logical-noise behaviour of the circuits are not derived.
- The longer-term 12C row is an evolution-scaling estimate on the legacy assumptions, not a complete response shot, with no step-error validation.

## Ch. 3, Mu2e and neutrinoless double beta decay (`ch03_mu2e_0nubb.py`)

- The 4He load dimension D is priced over 1e2-600 (`d_range`); a schematic Minnesota Hamiltonian reaches fidelity 0.999 at D = 77-103 (`d_he4_computed`), but a chiral (SRG-evolved) ground state could need more, and a narrower band would lower the 2028 headline toward the insertion alone.
- The 2033 campaign measures only the spin-dependent (tensor) element because the bare one-body SI elements in sd are pinned classically to 0.61% (`si_pinned_fraction`); if the deliverable is the VS-IMSRG-evolved operators, the readout of their two-body parts is not priced.
- <S_p>(27Al) = 0.30-0.34 and the per-shot variance of about one are estimates (`o_sd`, `sd_shots_per_hamiltonian`); the worst-case variance is about 8 times larger, and a shell-model value would replace both.
- Projection success p = 0.5 for 27Al is assumed (`projection_success`); it is computable exactly in sd, and the printed sensitivity p = 0.1-0.9 spans 7.5 d to 1.5 yr of campaign (`wall_campaign_p_band_s`).
- The 27Al survey cells use a generic 1-2 MeV gap band while the 48Ca/48Ti cells use the KB3G 0+ gaps of each truncated space; per-nucleus model gaps everywhere would move the 2033 cell.
- The daughter projection is priced at the parent's gap band; a daughter with a smaller gap raises the cost in proportion to 1/Delta.
- The resonant-drive (Rabi) readout is toy-checked only: lambda/|M| = 69, 164 and 600-720 and p_ref = 0.93 come from a schematic interaction, and Stark bias and Trotter/qDRIFT errors are not simulated.
- The Rabi walls are depth-limited (qDRIFT samples run one at a time, F* 3.5-7.7, `f_star_per_end`); a drive that groups commuting samples would recover the factories.
- The 2033 walls assume 12 commuting Trotter rotations in flight (`trotter_concurrency`, F* = 11.8); gate by gate the campaign floor is 135 d to 3.0 yr (`floor_wall_gate_by_gate_s`).
- Dead and cancelling PREPARE trees would remove about 2.2e3 T from the 2028 Hadamard-test shot (42% of the insertion); the count is kept as priced.
- Separating operator classes needs targets with different Z, N/Z and spin; a second target such as 48Ti in pf is not priced.
- No Mu2e-specific nuclear-matrix-element spread is quoted, and the 60% utility fraction (`utility_fraction`) is a stated choice first made for the 0nubb case.

## Ch. 4, hybrid lattice QCD (`ch04_hybrid_lqcd.py`)

- The polynomial degrees d (84 and the other rungs) and d_inv come from a linear program that needs scipy; they are carried `Assumed` and reproduced by `apply_log/ch04_roundB_lp_degree.py` and `apply_log/ch04_r11_lp_odd.py`, and at kappa = 1e3 the degree 5e3 is the fitted law extrapolated.
- The gauged-circuit prefactor of one T per unit of D^2 V log2 V is uncited (`scaling_prefactor`); every gauged block-encoding query is a SCALING primitive.
- The interacting condensate c lies between the free value 0.0275 (computed, `condensate_free`) and an estimate 0.08 (`condensate_interacting`); it sets the 2033 shot band, and an ensemble value at am = 0.04 on 8^3x16 would pin it.
- The box checks P against a classical evaluation of the same polynomial (d_inv = 500); calling the result Tr M^-1 needs d_inv about 700 and 1.4 times the per-shot T.
- The 2033 walls are depth-bound (F* = 1-4, 530-2100 s per shot, `wall_per_shot_s`); three parallel unary-iteration subtrees on about 100 idle logical qubits would give F* about 12 (`subtree_f_star`, `wall_campaign_unary_iteration_yr`) and are not adopted in the box.
- The Banks-Casher inputs (Sigma^(1/3) = 270 MeV, 4 tastes, 5 configurations, step-filter degrees 430/640) come from a simplified-observable pass; P is uncertain by about 2x in a 1.2 fm box.
- The post-2033 rung is priced at D V log2 V with the 2033 block encoding; D^2 V log2 V applies only if the D^2 there is meant to carry link loading.
- The disconnected-HVP plan needs 1.5e12 noise-matched shots per configuration at 24^3x48 (`hvp_noise_matched_shots_per_cfg_24c48`); `hvp_site_variance = 1` is assumed and is not conservative, and no coherent variance-reducing estimator is costed.
- Finite-density reweighting is not a deliverable: a phase estimator for arg det M to absolute error well below pi on an extensive quantity is not costed.
- The 2028 sampling target is 1% (3.1e5 shots); 2% rms would be 7.4e4 shots, and MLAE with k = 2 fits the cap only at the low end of the T band.

## Ch. 5, quark-gluon-plasma transport (`ch05_qgp_transport.py`)

- The crossover temperature of the simulated theory (Sigma(216x3) with H_I and three staggered fields) is unknown; it is priced at 150-200 MeV (`Tc_lat_2033_MeV`, Assumed) with the target T = 1.5 T_c^lat, so every 2033 range carries a cheap corner at 300 MeV and an expensive corner at 225 MeV (`T_corners_MeV`); below 150 MeV the slow end grows again.
- The in-situ diagnostic (`tc_diag_*`, about 1e3 preparation-only shots) locates the crossover in ramp energy, not in temperature; placing the campaign at 1.5 T_c^lat still needs an energy-to-temperature map for the simulated theory.
- Thermal preparation by a quench assumes eigenstate thermalization at the chapter's couplings on 3^3, a thermalization time c/T with c in [2, 2 pi] (`quench_c`), a one-spacing ramp as a floor (`quench_ramp_a_s`), relaxation of the momentum density within that time, and the thermal variance; all are checked only on small exactly solvable lattices.
- The preparation runs at the a_t/10 grid step; it needs energy conservation rather than a faithful unitary, so a 2-5 times coarser preparation step is an unpriced lever.
- The non-freezing H_I trajectory for Sigma(216x3) at 0.2 fm and at the second spacing is assumed by analogy with S(1080); no hop circuit exists for Sigma(216x3), so the hop is a Sigma(72x3) placeholder (`hop_group_2033`), and with 24 classes against 16 the true hop is likely larger (at twice the hop the step grows about 24%).
- The Sigma(72x3) electric term is a lower bound for a fast group transform that has not been constructed (CONJECTURE in `groups.py`); the compiled dense transform is about 270 times larger (`electric_t_per_link_step_2033_dense`).
- The fit excludes the early-time transient on a design that is not derivable from the chapter's inputs; fitting from t = 0 takes 0.37 times the shots and dropping Gamma t < 1 takes 2.7 times (`design_A_over_B`, `design_C_over_B`), moving every 2033 shot count and wall in proportion.
- The normalized correlator Cbar = 0.16 is carried without the block encoding of A_k that would give it; at 0.05 the shots grow about 10 times, at 0.3 they fall about 3.5 times.
- The decay band runs from 1/k to 1/(2 pi T) (`tau_band_fm`); the slow end is not conservative, since free streaming reaches 1/e only at kt = 2.93 (`free_streaming_1e_fm`) and strong-coupling holography makes the shear mode decay slower still.
- The campaign is 254-938 yr, 51-188 times the horizon (`campaign_two_spacings_yr_own`), a stated gap with no priced lever; a cost-optimal shot split (about 20%), a qubitized block encoding (the 10x trim is Stated and no block encoding exists), a coarser second spacing or a single spacing (12-29 yr) are the unpriced options.
- The second spacing is priced on the same 3^3 lattice (`second_spacing_a_s_fm`), so volume and k_min change with a and the two spacings do not give a continuum limit at fixed physics; holding L = 0.6 fm needs 4^3 at 2788-2888 logical qubits (`fixed_volume_4cubed_*`).
- A rate from two fit times is conditional on a third time confirming a single exponential; one more point at Gamma t about 1.2 adds 10-15% to the campaign and tests only a 30% departure (`test_time_*`), and is not priced.
- The 2028 shot count of 1e3 per time slice has no set precision; at 30% with an assumed relative fall-off of 0.1 over the 0.045 fm/c step it is 8.46e4 shots and 1.1 h (`shots_first_result_2028_at_30pct`), and because the correlator is even in t the fall-off at that step is closer to 0.02, which would raise the count about 25 times.
- The O(a^2), finite-volume and subgroup-truncation systematics have no derivation and are not printed as error bars.
- The average sign of the fermion determinant on the 0.6 fm box is not computed, and mu_B > 0 may be classically sampleable through the canonical ensemble, so the classical-comparator verdict is prospective.

## Ch. 6, collider physics (`ch06_collider.py`)

- The species readout N_h = U Pi_h U^dag uses strong-coupling local labels and the preparation ramp run backwards (`readout_ramps_2033`); the relation of label populations to hadron yields is not established, the backward ramp is not shown to be adiabatic far above the vacuum (it is priced as a lower bound), and the misidentification rate has no estimate at Sigma(72x3) 4x8.
- Vector and pseudoscalar mesons are degenerate at strong coupling; a per-meson V/PS readout needs a Fourier transform of hard-core link excitations in the label register, which is neither constructed nor priced, and the same object is needed for momentum-resolved N_h(P_h).
- On the 0.8 fm axis the back-to-back pair recollides after 0.4 fm/c (`recollision_time_fm`), so the 2033 counts are the hadron content of a box-confined state; the three-momentum trend uses the box modes 1.55, 3.1 and 4.65 GeV (`source_box_modes_GeV`), the top one at pa about 2.4, and whether box effects cancel in the trend is not established (`L_par_for_separation` = 20 would need 2260 logical qubits).
- The fragmentation instance prices the bilocal current without N_h, so its shots are lower bounds; the per-bin variance raises them 4.7 times (occupations equal to amplitude shares) to 23 times (`ff_nh_shot_factor`), and the momentum transform is unpriced.
- The dipole rows measure <|w|^2> rather than |<w>|^2 and correct the slope with a Markovian colour model (`floor_slope_ratio`); q-hat = 2 GeV^2/fm is taken as the fundamental value (695 shots if adjoint, `shots_qhat_if_adjoint`), and kappa/T^3 = 2.5, P_d = 0.9 and the 1/3 colour factor are assumed.
- The quench-prepared dipole media assume eigenstate thermalization at aT about 0.45 for Sigma(216x3) H_I and 2O; the ramp is not tuned to the target temperature, a coarser thermalization step (2-5x) is an unpriced lever, and the static Im V row is 5.4-15 times the 1e9-T budget.
- The species row runs Sigma(72x3) with H_KS at a about 0.1 fm in 2+1D, where the freezing point is unchecked; switching it to Sigma(216x3) would need a fermion hop that does not exist and 1092 logical qubits.
- The dipole circuits' T-depth is not analysed; their walls are serial T-count walls.
- The 10x per-shot cut from a qubitized Sigma(72x3) block encoding is an assumption; no block encoding or normalization is computed.
- n_s = 0.1, with c = 1 or 2 depending on whether n_s is the pair rate or the multiplicity, n_M = 0.5-1 and a 25% end-to-end trend are the declared working assumptions; the ratio shots scale as 1/trend^2 (`trend_delta_per_point`, `ratio_shots_at_10pct`).
- The adiabatic ramp of 1e2-1e3 steps runs twice per shot (preparation and readout); at 1e2 steps the campaign high end would fall from about 26 yr to about 6 yr.
- Source-free calibration circuits that measure the vacuum fake rate b multiply the shots by (sqrt(1 + r) + sqrt(r))^2 with r = b/n_s; they are in no printed number.
- At L_par = 8 the momentum quantum is 1.55 GeV, so about 3 of the 10 z-bins of readout (b) are resolvable.
- The 2033 hop runs one link at a time because its colour-angle flags fill the workspace (F* 7.0-18.5, `f_star_per_circuit_2033`); a streamed-flag schedule would give F* 26-61 at about 5% more hop T.

## Ch. 7, quantum fields in curved space (`ch07_curved_space.py`)

- The 2033 step m_phi dt = 0.05 comes from a mean-field convergence proxy (`verlet_proxy`) and is checked exactly on a 3-site ring with one Wilson fermion (`exact_reduced_step`, error at most 0.3%, at most 2% scaled to 5^3); the 3D case and the Pauli-string split of the hop are not tested.
- The condensate amplitude Phi0 about 0.3 m_phi is set by what the K = 16 register holds (`condensate_quanta_per_site`); the scalar shells then sit at n_k below 0.06 and are read to only 45-65%, so the 10% target is on the fermion; the alternatives (b about 1.5/m_phi at twice the cost, or K = 32 at +125 logical qubits) are priced and not taken.
- At the couplings (lambda, g) = (1, 1), (3, 0.5), (1, 2) a Gaussian truncation reaches the late fermion n_k to within 5.3% (`twopi_lo`), so the 5^3 instance is a cross-check, not an advantage; an advantage needs stronger per-site coupling or the scalar shells, which the truncation misses by up to 66%.
- The 2028 run is one step at K = 8 (`ds2028_grid_check`); its one-step map is 13-25% from exact de Sitter, light fields m below H fail beyond two steps at any K up to 16, and the chapter should state m about 1-1.4 H_inf.
- The 50-step dt sweep between ramp and evolution (3.6e7 T, 0.24% of the shot, `t_dt_sweep`) is stated and not priced.
- The lattice spacing b about 3.5/m_phi is implied by omega_max about m_phi and is to be confirmed.
- The background a(t) proportional to t^(2/3) is prescribed; self-consistent backreaction is not priced.
- The fault budget counts T gates only; idle storage (1.6e12 qubit-cycles per 2033 shot, `idle_qubit_cycles`) would dominate at the same rate per location, a late fault shifts the fermion n_k by at most 0.12 per mode (`fermion_fault_shift`), and early faults are not bounded.
- The shot count N_set ((n + 1/2)/n)^2 / (m' eps^2) at n = 1, m' = 3 with a 1-2x margin (`shots_per_profile_campaign`) covers fermion shells with n_k above 0.05; whether the target is the fermion spectrum or a scalar precision is open.
- The scalar energy density, pressure and Yukawa readout are a first equation-of-state proxy with unestimated statistical error.
- The campaign is 10-20 profiles (`profiles_campaign`); which parameters are scanned is not fixed.
- The preparation-only vacuum reference run (2.4% of a shot, `vacuum_reference_fraction`) and the fermion out-basis readout (at most 7.4e6 T, `t_fermion_readout_bound`) are not priced, and the 2028 Gaussian-network preparation (3e4 T, `t_prep_2028`) is stated, not derived.

## Ch. 8, electroweak baryogenesis (`ch08_baryogenesis.py`)

- The free part of the Delta-n state preparation is priced (4.6e5 T per application, 1.6% of the evolution, `free_prep`); the interacting dressing ramp is not, and any ramp lowers the MLAE depth (`prep_sensitivity`: r = 0.2 gives depth 27 and 9.7 yr per campaign point, r = 1 gives depth 15 and 29 yr); a static wall needs a degenerate well or a pinning term.
- The flag contrast p0 = 0.66 at a m_f = 0.5 on the free vacuum, 0.60-0.65 on the wall background (`flag_contrast`), multiplies the queries by about 2.3; the chapter does not fix a m_f.
- By an exact charge-conjugation-plus-flavour-swap symmetry (`c_prime_residual` = 0) the flavour-summed asymmetry vanishes identically, so Delta n is flavour-resolved.
- The bath for the nucleation circuits is estimated, not priced: one unit of Lindblad time is about 5.5e5 T (`bath_sweep_t`), one application of each of the 240 scalar jumps is 1.3e8 T (five nucleation circuits), and the sweeps needed to mix and the dissipator rate during the evolution are open; the Gibbs preparation and the basin projector are carried at 0 T, UNSOURCED.
- The 2028 count 1.6e5 T includes an unbanked 3.3x reduction of c_T that is an assumed development target (`ct_reduction_banked`) and omits the adiabatic preparation and basin projector, so it is neither a lower nor an upper bound.
- The integrated asymmetry S_int is not computed; the first result takes 1e-3 (`s_int_taken`) and its wall scales as (1e-3/S_int)^2.
- The 2028 slope uses times 0.4 and 1.2/m_phi, and 0.4/m_phi may sit in the Zeno regime; a shorter lever arm costs 5-6 times the shots.
- At T = T_c the nucleation rate vanishes in infinite volume; the finite-volume decrement on 15x8 may be far below 2.5e-3, so that point may cost more than the 12 d priced.
- The scan over N_bin = 100 exceeds what a 15x8 lattice resolves; a scan over temperature alone would cut it from 3.3 yr to about 0.3 yr.
- At the 0.1-fault budget the deepest circuit has p_f = 0.095 (`fault_prob_per_circuit`), which shifts A_R by up to 4.4e-3 (`fault_bias_worst_A_R`), above the per-bin 1e-3; a null circuit quadruples the first-result queries (`null_control_query_factor`), and idle locations would need about 1e-12 per location.
- The tensor-network cost of the 2+1D quantum-scalar comparator is not estimated.
- If the 3.3x cut applies to T-count but not depth, F* = 7-9 and the 2028 wall is depth-bound at 1.1-1.4 times the quoted 54 min.

## Ch. 9, chiral lattice fermions for real-time gauge dynamics (`ch09_chiral_gauge.py`)

- A typical reference state for the thermal filter succeeds with probability Tr e^(-beta(H - E0))/D, about 1e-60 at 3^3 and beta = 2a (`typical_reference_tpq`); an amplified random basis-state reference is biased, a Gibbs(H0) reference does better but uncontrolled, and the priced 17 filter calls (`n_filter_calls`) are a placeholder.
- The thermal filter is priced on the fermion block encoding only; the electric, magnetic and Higgs terms inside e^(-beta H/2) raise the normalization above 2Q = 432 (`be_normalization`) and add their own cost, neither priced.
- The Chern-Simons estimator is the Kubo form with link-averaged E.B (`ncs_lattice_background`); its Hadamard-test variance against the free-field background is at least 1.1e14 (`ncs_hadamard_rel_var`), so the assumed readout variance of about 1e2 needs a readout whose noise is set by delta N_CS itself, which is open together with the multiplicative renormalization of the rate operator.
- At g^2 = 1 and T a = 0.5 the box has L g^2 T = 1.5 (`box_L_g2T`), where the classical sphaleron rate is suppressed about 1e3 times relative to large volume (`mr_small_box_suppression`); no rate is claimed at 3^3, and the Higgs-mass matching and the location of T_c are open.
- The Ginsparg-Wilson projected Yukawa is priced as a second sign sequence per application (`n_sgn_per_application`), but its weight y rho_max in the block-encoding normalization is not.
- The condensate arm is priced on the prepared thermal state (`condensate_shot`); if the ground state is meant, its preparation has no price.
- The step count is the state-dependent estimate, 450 steps at t = 10a (`n_trotter_2033`, `trotter_rule`), against a worst-case bound of 7.4e3; a criterion that keeps the common phase would need about 7e2 steps, T a = 0.5 and g^2 = 1 are choices, and the Higgs modes would add 15-20%.
- A self-inverse block encoding is assumed so that applications per step equal the Jacobi-Anger order; otherwise the count doubles (`t_per_shot_if_not_self_inverse`).
- The gauge and Higgs terms of each Trotter step (gauge estimated at about 2% of the step, `unpriced_gauge_share_of_step`) and the amplification reflections (about 2e5 T per shot, `unpriced_fpaa_reflection_T`) are not priced.
- The window t = 10a is marginal for the linear regime of Chern-Simons diffusion, which sets in after a time of order 1/(g^4 T).
- The diabatic error of the 20-step ramp is not estimated; the 2028 benchmark is meant to measure it.
- The 2028 benchmark takes m_f = 0 (no wall-joining term); m_f different from zero adds about 4% to the T count.
- The 2028 register holds the chunked circuit at 2.9e5 T on 250 logical qubits (`chunked_t_per_shot`); the unchunked 2.7e5 needs about 400 (`lq_unchunked`), and F* of at least 10 needs three concurrent spatial-hop groups per slot.
- The 2033 register of 1015 includes five qubits for a CNOT fan-out of the shared control (`lq_ancilla_fanout`); without it F* is 5.6-6.3 and every wall is 1.6-1.8 times longer.

## Ch. 10, QCD at finite baryon density (`ch10_finite_density.py`)

- The 2033 QSVT/TPQ price is conditional on a reference whose filtered amplitude is about 0.1 (`qsvt_amplitude`) with thermal eigenstate weights spread correctly over the N_B and field-number sectors; a Haar reference succeeds with probability Z e^(beta E0)/D, a product reference is biased with amplitude falling exponentially in V, and no candidate closes it.
- The polynomial degree is set in ||H|| rather than a block-encoding normalization (`d_beta`); at the Ch. 9 ratio lambda/||H|| = 8.64 it would be 89 and the shot 3 times larger, and the query cost of 2-4 Trotter steps (`qsvt_c_step`) is stated.
- The 17 filter calls (`qsvt_filter_calls`) are the standard amplitude-amplification count at amplitude 0.1; a true fixed-point sequence needs about 31 calls, 1.8 times (the same applies to Ch. 9).
- kappa_2 of N_B at V = 2^3 and its temperature dependence set the grid shots: 2.6e5 at 0.01 (`nb_toy_k2_priced`), 9.1e5 at 0.05 and 2.2e6 at 0.1 (`grid_shots_k2_mid`, `grid_shots_k2_hi`); a classical exact-diagonalization or strong-coupling estimate would settle the band.
- The target on chi4/chi2 is an absolute 0.1 (`delta_chi4_chi2_abs`); at small kappa_2 the ratio is pinned near 1/9 or 1 by the quark or baryon content, so the precision on the deviation from that baseline may be the real target.
- The 2028 box prices one ramp step as a lower bound (`hard_ops_is_lower_bound`, `ramp_original_*`); the executed circuit is the N-step ramp plus the grouped energy readout, about 4.2e6 T, and lambda, mu_c/lambda, the readout grouping and the ramp length are open, as are the term variances of the energy estimator behind the 2e3 shots.
- The 5% pressure target references the classical mu_B = 0 pressure of this exact Hamiltonian; if it is not available, p needs the uncosted energy channel.
- The 2033 T-depth is the 2e3-step factory band scaled by each end's T ratio (`f_star_matched` = 17.5-56, Assumed); the shared `depth_exports` pairs the bands crosswise (F* 8.6-114), so the walls are taken serial at the matched ends.
- The digitized Sigma(36x3) path integral has its own sign problem at mu_B = 0, so it does not measure the mu_B obstruction; the classical-comparator verdict of "probably a validation" is not established, and the Euclidean staggered phase at the tuned g^2 would close it.
- The Gibbs-sampler record (8.4e13-4.4e14 T per shot, `gibbs_t_per_shot_range`) assumes a mixing time of 20 sweeps (`c_mix`).
- The 2028 workspace itemization does not list the synthesis ancilla as its own line (29 with serial reuse, 55 without; `ancilla_workspace_*`).

## Shared across chapters

- The staggered and Wilson hops of Chs. 5, 6, 9 and 10 come from an unpublished gate-count draft (`groups.hop_link_cost`); the frame undo, the SU(3) squish uncompute, the per-field colour rotations, the C^nX convention (n - 1 Toffolis), the diagonalizer slope (4 rotations), the Sigma(72x3) class count (taken from Sigma(36x3)) and 7 T per flag Toffoli are this package's readings where the draft is silent or inconsistent; Ch. 6 keeps a strict xfail on the colour-frame count as a tripwire.
- The fault budget of 0.1 expected faults per shot counts T gates only (`eps_l = 0.1 / N_T`); Clifford, idle, measurement and magic-state-injection faults are excluded, only Ch. 7 bounds the bias on its observable, and a per-instance detection model is not derived.
- The Hamming-weight-phasing ancilla count k - w(k) plus one repeat-until-success ancilla (`common.hwp_ancilla`, `common.RUS_ANCILLA`) is applied everywhere; Ch. 10 does not itemize the synthesis ancilla.
- The reference budgets (1e5 and 1e9 hard ops, 150-250 and 1000 logical qubits, eps_l 1e-8) and the seconds per year (3.156e7 in some modules, 3.15576e7 in others) are typed per module; no printed number depends on the difference.
- The illustrative dollar figures are stated shares of cited budget lines (Ch. 3: 60% of the Mu2e project cost; Ch. 4: 17% of a one-time capital base of about $150M; Chs. 5 and 10: 5% and 2% of the DOE NP Heavy Ion research line), rescalable and not derived; Chs. 8 and 9 attach none.
