# Chapter 9: chiral lattice fermions for real-time gauge dynamics

Model: `estimates/ch09_chiral_gauge.py`. Tests: `estimates/tests/test_ch09_chiral_gauge.py`. Chapter text: `applications/app06_chiral_gauge.tex` in the paper source (not in this repository).

Every number below was produced by the model as it stands, through `python -m estimates --chapter ch09`, `resources.json`, `model(Assumptions(), era)` or a direct call of one of the model's functions (named where used). Where a number is pinned by a test, the test is named.

## 1. What is estimated

The chapter has two instances and no co-design row (`model(a, era)` raises for any era other than `2028` and `2033`).

**2028 benchmark: 1+1D domain-wall fermions with a Z3 gauge field.** One Dirac flavor, three colors, L = 8 sites, fifth-direction extent L5 = 4, gauge group Z3 inside SU(3) at 2 qubits per link. The state is prepared by an adiabatic ramp of 20 Trotter steps, and one probe step follows. The observable is the chiral condensate ⟨ψ̄ψ⟩ and its residual-mass scaling in L5, compared with a classical Z3 domain-wall calculation at the same L5. The benchmark exercises the register layout and the fifth-direction hopping. It does not exercise the overlap sign function, non-Abelian digitization, the Higgs sector or thermal-state preparation. The first result is L5 = 4 at 30%; the campaign is L5 = 2, 3, 4 at 10%.

**2033 target: SU(2)-Higgs with overlap fermions on a 3^3 spatial lattice.** SU(2) is digitized to the binary octahedral group 2O (48 elements, 6 qubits per link). One Dirac flavor of overlap fermions, a Higgs field Φ = ρU with the radial mode on N_ρ = 16 levels, and a Ginsparg-Wilson-projected Yukawa. The Hamiltonian is ψ†h_ov[U]ψ with h_ov = γ⁰ + sgn(H_W[U]); the sign function is a quantum-signal-processing (QSP) polynomial of the single-particle Wilson kernel, lifted to Fock space by the construction of arXiv:2607.28524. A thermal state at T a = 0.5 is prepared by a quantum-singular-value-transformation (QSVT) filter on a reference state and evolved to t = 10a. The observable is a method test: the early-time growth of the full ⟨(ΔN_CS)²⟩ (background included) at t = 5a and 10a at one temperature to 30%, compared with the free-field value and with classical real-time simulation. The campaign adds the condensate at 4-8 temperatures to 10%. No sphaleron rate is claimed at 3^3.

Headlines (`resources.json`): 2028, 250 logical qubits (LQ) and 2.690e5-2.885e5 T per shot; 2033, 1015 logical qubits and 1.860e10 T per shot.

## 2. How to run

From the repository root (Python 3 with `numpy`; see `requirements.txt`):

```
python -m estimates --chapter ch09                    # assumptions, resources.json rows, exports per era
python -m estimates --chapter ch09 --intermediates    # also every intermediate the model records
python -m pytest estimates/tests/test_ch09_chiral_gauge.py -q
python applications/ch09_chiral_gauge/dump_assumptions.py   # the table in section 3
```

`--chapter ch09` recomputes the model and writes nothing. It prints:

- the 217 assumption fields with value, provenance and source;
- the two rows this chapter writes to `resources.json`;
- for each era, the headline against the printed box, shots, wall time, ε_l, and the depth and factory exports;
- the primitive breakdown (count × T each = T total, with its circuit status);
- the model's own notes.

The 2033 model computes free-field backgrounds with `numpy` and takes a few seconds.

The test file has 83 tests: 79 pass and 4 skip here (`test_r23_no_machine_count`, `test_r25_chapter_prints_the_model`, `test_utility_box_has_no_dollar_figure`, `test_chapter_open_prints`), because they read the chapter's `.tex` from `../applications/` (the repository's parent directory, `conftest.PAPER_ROOT`); with the paper source there, all four pass. `estimates/tests/test_rhop.py` (26 tests) also concerns this chapter: it pins the Z3 spatial-hop rule `groups.rhop_link` priced in the 2028 step and checks that the Hamming-weight-phasing helpers of this module agree with `common.py` for k = 1-256. `test_fermion_primitives.py` covers `groups.hop_link_cost`, which this chapter uses only in a consistency check of the Z3 hop (lines 1582-1584) and in the second-quantized comparison of section 4.3.

## 3. Assumptions

The table is printed by `dump_assumptions.py` from `Assumptions()` (lines 874-1432 of the model). Name, value and provenance are the code's. Sources and notes are quoted from the code, whitespace collapsed, notes cut at 150 characters (the full text is in the model). Notes use the code's own labels for report-wide conventions (R-TOL for the synthesis tolerance, R3 for the fault budget, R5 for 7 T per Toffoli, and so on); the conventions themselves are summarized in the top-level README.

Rows are grouped by instance and by role. *Inputs* are the fields the code does not mark otherwise. Most feed the headline; perturbing each by 20% shows that some feed only the comparisons and sensitivities of section 4.4. These are: the two envelopes, the 1+1D overlap comparison (`kappa_1p1d_*`, `delta_1p1d_a`, `sgn_T_1p1d_*`, `n_sgn_gs_prep`), the worst-case Trotter bound (`casimir_max_2O`, `hop_n_spin`), the Chern-Simons readout and unit checks (`ncs_e_max`, `ncs_b_max`, `tc_ew_gev`, `mr_rate_*`), the thermal-route alternatives (`haar_*`, `gibbs_*`), and `m0_wilson`, `L5_backup`, `kappa_physical`, `dt_empirical_fm`, `shots_2033`, `campaign_yr`, `lever_eps_obs`, `fpaa_reflection_ancilla`, `lcu_strings_per_bilinear`, `hop_mcx`. *Printed figures* are fields named `*_stated`: numbers the chapter prints, which the tests compare with the model's own value (their notes give the printed text and the model value it rounds). *Comparison records* are fields whose source or note carries the code's record marker: inputs of alternative constructions the model evaluates so that their costs can be compared. Perturbing each of them leaves every headline (qubits, T, shots, wall, ε_l) unchanged.

217 fields: Cited 31, Assumed 62, Stated-not-derived 124.

#### Shared by both instances: inputs (9)

| name | value | provenance | source | note |
|---|---|---|---|---|
| `eps_syn` | 0.01 | Cited | TRACKED_CHANGES.md R-TOL (H. Lamm, 2026-09-29) | app06:108 'the total synthesis error per shot is 1e-2, and the tolerance is set from the circuit's rotation count' |
| `rot_synthesis` | rus | Cited | bocharovRoettelerSvore2015,campbell2017 | app06:108 'from the fit, 1.15 log2(1/eps) + 9.2' (repeat-until-success, randomized); the chapter's rotations other than the R-HOP hop |
| `toffoli_convention` | textbook | Cited | main-overview:47 (ruling R5, H. Lamm 2026-09-28) | 'it uses 7 T per Toffoli throughout'; app06:99 'with every Toffoli priced at 7 T'; was the 4-T measurement-assisted Toffoli of arXiv:1709.06648 |
| `t_gate_s` | 1e-06 | Assumed | - | app06 boxes 'at 1 us/T-gate' |
| `shot_overhead_s` | 1e-04 | Assumed | - | per-shot overhead ~0.1 ms (register initialization, final readout, decode), report rule as Chs. 8 and 10; moves no printed number here |
| `faults_per_shot` | 0.1 | Assumed | - | ruling R3 (TRACKED_CHANGES.md 2026-09-28): eps_l = 0.1 expected faults per shot, report-wide; app06:103 'at 0.1 expected faults per shot' |
| `sigma2_condensate` | 10 | Stated-not-derived | app06 shot-count paragraph | 'sigma^2 ~ 10 for the average over the 27 sites of the 2033 lattice' |
| `eps_obs` | 0.1 | Assumed | - | app06 'at eps=0.1'; box: '10% statistical' |
| `eps_first` | 0.3 | Assumed | - | ruling R1 (H. Lamm 2026-10-02): the minimal first result at 30% statistical error ('for some of this physics even 30% statistical error would be inter ... |

#### Shared by both instances: comparison records (2)

| name | value | provenance | source | note |
|---|---|---|---|---|
| `t_per_rot` | 30 | Cited | bocharovRoettelerSvore2015,campbell2017 | RETIRED (R-TOL 2026-09-29): app06:108 was '~30 T-gates per rotation (~25 from the fit, rounded up)' at eps_gs ~ 1e-4; record only |
| `eps_rot` | 1e-04 | Assumed | - | RETIRED (R-TOL 2026-09-29): app06:108 was 'eps_gs ~ 1e-4'; record only |

#### 2028 benchmark (1+1D Z3 domain wall): inputs (26)

| name | value | provenance | source | note |
|---|---|---|---|---|
| `envelope_2028` | 1e+05 | Cited | DOE_RFI_2026 | app06:105 '<= 1e5 T-gates' |
| `gauge_group_2028` | Z3 | Assumed | - | app06:105 'Z3-digitized 1+1D toy' |
| `D_2028` | 1 | Assumed | - | app06:105 '1+1D' |
| `L_2028` | 8 | Assumed | - | app06:105 'L=8' |
| `L5_2028` | 4 | Assumed | - | app06:105 'truncated L5=4' |
| `Nf_2028` | 1 | Assumed | - | app06:105 'N_f=1' |
| `Nc_2028` | 3 | Assumed | - | app06:105 'three colors ... decoupled identical copies' |
| `anc_rus` | 1 | Cited | bocharovRoettelerSvore2015; ruling E26 (r22, H. Lamm 2026-10-01); common.RUS_ANCILLA | repeat-until-success synthesis holds ONE ancilla while a rotation is in flight ('Clifford+T basis, plus one ancilla qubit and measurement'); the weigh ... |
| `n_pauli_per_hop` | 2 | Stated-not-derived | app06:108 | 'two Pauli strings per hopping bilinear' |
| `n_pauli_per_electric` | 3 | Stated-not-derived | app06:108 | 'three diagonal Z-strings per two-qubit link, each grouped across the eight links' |
| `fifth_boundary` | open | Cited | Kaplan_chiral_lattice_fermions, arxiv_2505_20419; r21 (3a), H. Lamm 2026-10-01 | DERIVED HERE: domain-wall fermions have open boundaries in the fifth direction. The chiral modes sit on the walls s = 1 and s = L5, and the only term ... |
| `spatial_hop_grouping` | color+s | Stated-not-derived | app06:108 | 'The twelve color and s copies of each spatial hop share one s-independent link and form one group per link and Pauli string: $64$ groups of $k=12$' ( ... |
| `hop_synthesis` | rus | Cited | r17 ruling (e), H. Lamm 2026-10-01 (apply_log/r17_core.md, r17_ch09.md); FermionPrimitives_unpub section_resources.tex:7-10 | the Z3 spatial hop's rotations at the full fit 1.15 log2(1/eps) + 9.2, the chapter's own price and the draft's convention for every fermion rotation; ... |
| `onsite_term_2028` | wilson | Cited | arxiv_2505_20419; FermionPrimitives_unpub section_mass.tex:12 | r17 ruling (f), H. Lamm 2026-10-01: the on-site part of h_DW = gamma^0 D_W, coefficient m + r(1+d) (= m + 2 at d=1, r=1; non-zero across the topologic ... |
| `hwp_chunk_k` | 32 | Stated-not-derived | app06 2028 derivation | 'Splitting the k=72 and k=192 groups into groups of at most k=32 needs 31 and fits' (E26: one ancilla per adder, no separate weight register; was '37 ... |
| `n_trotter_2028` | 20 | Stated-not-derived | app06:109 | 'We shorten the adiabatic ground-state ramp from N_Trotter ~ 30 to N_Trotter = 20 steps' (ruling R6) |
| `n_trotter_2028_before_r6` | 30 | Stated-not-derived | app06:109 | 'from N_Trotter ~ 30'; the ramp length before R6, emitted so the trade is visible |
| `n_probe_steps` | 1 | Stated-not-derived | app06:109 | 'With a single-time-step probe ... the shot is 21 steps' |
| `kappa_1p1d_a` | 30 | Assumed | - | app06 2028 derivation 'even at kappa ~ 30' |
| `delta_1p1d_a` | 0.01 | Assumed | - | app06 2028 derivation 'delta_sgn ~ 1e-2' |
| `sgn_T_1p1d_a` | 56000.0 | Stated-not-derived | app06 2028 derivation | 'one application of the block-encoded overlap Hamiltonian costs 5.6e4 T' (140 queries + lift; r21 3.0e5 in the second-quantized pricing) |
| `kappa_1p1d_b` | 50 | Assumed | - | app06 algorithm requirements 'at kappa >~ 50' |
| `sgn_T_1p1d_b` | 90000.0 | Stated-not-derived | app06 algorithm requirements | 'already costs ~9e4 T' (9.28e4; r21 ~5e5) |
| `n_sgn_gs_prep` | 10 | Stated-not-derived | app06 algorithm requirements | 'ground-state preparation needs at least O(10) of them' (a floor under the lifted normalization 2Q = 96) |
| `spatial_concurrency_2028` | 3 | Assumed | - | r25 R10, DERIVED HERE: three spatial-hop phasing groups per slot on mutually non-adjacent links (disjoint sites, commuting strings); 3 x 10 adders fit ... |
| `L5_campaign_2028` | (2, 4) | Assumed | - | shot audit (r25): residual-mass scaling needs at least three L5 values; the 2028 campaign is L5 = 2, 3, 4 at the box's 1125 shots each |

#### 2028 benchmark (1+1D Z3 domain wall): printed figures (23)

| name | value | provenance | source | note |
|---|---|---|---|---|
| `t_per_step_2028_stated` | (12800.0, 13700.0) | Stated-not-derived | app06 2028 derivation | '~1.28e4 T per step' unchunked, '~1.37e4 T per step' chunked to k <= 32 (r21 a: 24 fifth-direction hops; r17 1.31e4 / 1.42e4) |
| `t_per_shot_2028_stated` | (2.69e+05, 2.89e+05) | Stated-not-derived | app06 2028 derivation | '~2.69e5 T-gates unchunked and ~2.89e5 chunked' (r21 a; r17 2.76-2.97e5) |
| `ratio_to_envelope_2028_stated` | (2.7, 2.9) | Stated-not-derived | app06 2028 derivation | 'a factor 2.7-2.9 above the first-generation 1e5 envelope' |
| `t_per_rot_coherent_stated` | 31 | Stated-not-derived | app06:109 | 'per-rotation cost would rise to ~31 T' (R-TOL: the same 1e-2 total summed linearly; was ~40 at 1e-4) |
| `t_per_step_coherent_stated` | 16000.0 | Stated-not-derived | app06 2028 derivation | 'the unchunked step to 290 x 31.3 + 991 x 7 ~ 1.60e4 T' (r21 a; r17 1.64e4) |
| `t_per_shot_coherent_stated` | 3.4e+05 | Stated-not-derived | app06 2028 derivation | 'and the shot to ~3.4e5 T' (3.363e5) |
| `hwp_ancilla_unchunked_stated` | 190 | Stated-not-derived | app06 2028 derivation | 'The k=192 group holds 190' (E26, r22: k - w(k); was 'about 200' with the weight register counted twice) |
| `hwp_ancilla_chunked_stated` | 31 | Stated-not-derived | app06 2028 derivation | 'needs 31 and fits' (E26, r22; was 37) |
| `lq_unchunked_stated` | 400 | Stated-not-derived | app06 2028 derivation | 'the unchunked end needs 191, a register of about 400 LQ' (192 + 16 + 190 + 1 = 399; was ~410) |
| `lq_single_copy_stated` | 110 | Stated-not-derived | app06 2028 derivation | 'rather than the ~110 LQ a single copy would need' (64 + 16 + 32 = 112; was ~120) |
| `eps_l_2028_stated` | 3.5e-07 | Stated-not-derived | app06 2028 derivation and box | 'eps_l <~ 3.5e-7' at 0.1 expected faults per shot (R3), from the chunked 2.885e5 T |
| `shots_2028_stated` | 1100.0 | Stated-not-derived | app06 2028 box | '1.1e3 (psi-bar psi at one coupling)'; 1125 from sigma^2 = 11.25 (r21 d; was 1e3, borrowed from V = 27) |
| `shots_2028_worst_stated` | 9000.0 | Stated-not-derived | app06 shot-count paragraph | 'at most sigma^2 ~ 90, 9e3 shots' if the sites of a color copy are fully correlated |
| `wall_per_shot_2028_stated` | (0.27, 0.29) | Stated-not-derived | app06 2028 derivation | '0.27-0.29 s per shot at 1 us per T-gate' (r21 a; r17 0.28-0.30) |
| `c_sp_1p1d_stated` | 400.0 | Stated-not-derived | app06 2028 derivation | 'a query is 40 Toffolis and 6 synthesized rotations, ~4.0e2 T' (393-406 over the four circuits) |
| `gs_prep_1p1d_stated` | 5.8e+05 | Stated-not-derived | app06 2028 derivation | 'at least O(10) of them, ~5.8e5 T' (r21 3.1e6) |
| `gs_prep_over_dwf_stated` | 2 | Stated-not-derived | app06 2028 derivation | 'about six times the budget and twice the domain-wall shot' (5.78e5 / 2.885e5 = 2.0; r21 'eleven') |
| `gs_prep_over_envelope_stated` | 6 | Stated-not-derived | app06 2028 derivation | 'about six times the budget' (5.8) |
| `gs_prep_1p1d_b_stated` | 1e+06 | Stated-not-derived | app06 algorithm requirements | '~1e6 T, an order of magnitude past the 1e5 envelope' (9.55e5) |
| `wall_s_first_2028_stated` | 36 | Stated-not-derived | app06 2028 box | first result '36 s' (125 shots; 36.08) |
| `wall_min_campaign_2028_stated` | 13 | Stated-not-derived | app06 2028 box | campaign '13 min' (L5 = 2, 3, 4; 13.24) |
| `lq_2028_stated` | 250 | Stated-not-derived | app06 2028 derivation and box | '250 = 192 + 16 + 42' (r25 R10) |
| `depth_step_2028_stated` | (1100.0, 1300.0) | Stated-not-derived | app06 2028 derivation | 'a T-depth of 1.1-1.3e3' per step (1059.7-1318.7) |

#### 2028 benchmark (1+1D Z3 domain wall): comparison records (3)

| name | value | provenance | source | note |
|---|---|---|---|---|
| `ancilla_2028_asserted` | 40 | Stated-not-derived | app06:145 (pre-r20) | RECORD (E23 (5)): the asserted '40 ancilla for the adiabatic preparation, measurement and Hamming-weight phasing'; replaced by the itemized peak (39 a ... |
| `hop_synthesis_legacy` | rus-slope | Cited | TRACKED_CHANGES.md 'RULING, R-HOP report-wide' (H. Lamm 2026-09-29) | RETIRED by r17 ruling (e): R-HOP's slope-only price, kept so the 2.0e5 record stays reproducible |
| `c_be_1p1d_sq_stated` | 2200.0 | Stated-not-derived | app06 (r21) | RECORD (second-quantized pricing, r21 e): 'c_be ~ 2.2e3 T per query' (240 Toffolis + 25 rotations) |

#### 2033 target (3D 2O overlap + Higgs): inputs (66)

| name | value | provenance | source | note |
|---|---|---|---|---|
| `envelope_2033` | 1e+09 | Cited | DOE_RFI_2026 | app06:102 '~1e9 2033 resource limit' |
| `gauge_group_2033` | 2O | Assumed | - | app06 '\|G\|=48 [2O binary octahedral]' |
| `D_2033` | 3 | Assumed | - | app06 'D=3' |
| `L_2033` | 3 | Assumed | - | app06 'V=3^3=27' |
| `Nf_2033` | 1 | Assumed | - | app06 'N_f=1' (Witten anomaly minimum) |
| `Nc_2033` | 2 | Assumed | - | app06 'N_c=2 [SU(2)]' |
| `N_rho` | 16 | Assumed | - | app06 'N_rho=16' |
| `sgn_architecture` | single_particle_lift | Cited | arxiv_2607_28524; ruling E26 (i), H. Lamm 2026-10-01 ('do whatever is easiest and cheaper') | the Hamiltonian is the bilinear psi^dag h_ov[U] psi with h_ov = gamma^0 + sgn(H_W[U]); sgn acts on the single-particle kernel, a Q x Q matrix of link ... |
| `m0_wilson` | 1.5 | Assumed | - | r22: the negative Wilson mass -m_0 of the overlap kernel, inside the doubler window 0 < m_0 < 2 (the chapter states none). Used only for the normaliza ... |
| `eps_qsp` | 0.01 | Assumed | - | r22: total QSP (Jacobi-Anger truncation) error of the kernel evolution over one shot; each Trotter step gets eps_qsp / N. Same size as the synthesis b ... |
| `anc_readout` | 1 | Assumed | - | E23 (5), r20; r22 wording (validation): control qubit of the Ramsey / Hadamard-test readout of the unequal-time Chern-Simons correlator (arxiv_2208_13 ... |
| `anc_qsp_signal` | 1 | Assumed | - | E23 (5), r20; r22: one QSP signal qubit, shared by the TPQ filter and the kernel evolution (a real polynomial needs one signal qubit, gilyen2018QSingV ... |
| `anc_fpaa` | 1 | Assumed | - | E23 (5), r20: marker qubit of the fixed-point amplitude amplification (used at the amplification reflections only) |
| `anc_h_lcu` | 2 | Assumed | - | E23 (5), r20: LCU index over the parts of H (overlap fermion part, gauge, Higgs) |
| `anc_dov_lcu` | 1 | Assumed | - | E23 (5), r20: LCU qubit of h_ov = gamma^0 + sgn(H_W) (the chapter's D_ov = 1 + gamma5 sgn(H_W)) |
| `anc_sgn_signal` | 1 | Assumed | - | E23 (5), r20: QSP signal qubit of sgn(H_W) |
| `anc_sgn_real_part` | 1 | Assumed | - | r22: the qubit that selects the two phase sequences +-Phi whose average is the real sign polynomial (gilyen2018QSingValTransfArXiv) |
| `anc_lift_type` | 2 | Cited | arxiv_2607_28524 (Supplement theorem) | the two-qubit register A of the lift: four Majorana types (X or Y on each side) |
| `fpaa_reflection_ancilla` | 2 | Cited | arxiv_2407_17966 | r22 (validation tier A): the amplification reflects about the reference state on the whole register; an n-controlled NOT costs 2n Toffolis on two clea ... |
| `L5_backup` | (8, 16) | Assumed | - | app06 route selection 'L5 factor (8-16 for the back-up here, 16-32 in production)' |
| `kappa_2033` | 30 | Assumed | - | r24 (H. Lamm 2026-10-02, 'lower kappa is ok'): app06 'kappa = 30 (a heavier, coarse-mass regime)'; was 1e2 (r17-r23). Since r22 kappa = alpha / Delta, ... |
| `delta_sgn` | 0.001 | Assumed | - | app06 'delta_sgn ~ 1e-3'; the box budget row prints '~0.1% (1e-3)' (ch09-delta-sgn) |
| `lcu_strings_per_bilinear` | 2 | Cited | Jordan-Wigner encoding (standard) | a color-diagonal bilinear a^dag b + h.c. is (XX + YY)/2 with a Z string: two Pauli strings (record: the hop strings of the second-quantized LCU) |
| `lcu_toffoli_per_term` | 1 | Stated-not-derived | app06 2033 derivation; Babbush_PRX_2018 | one SELECT Toffoli per LCU term (unary iteration over L terms costs L - 1): the 40 Higgs strings per site of the headline, and the second-quantized re ... |
| `n_umult_per_higgs_site` | 1 | Stated-not-derived | app06 2033 derivation | '27 compiled group multiplications U_x' (ruling VERTEX, kept for the Higgs; the cost is read from GROUPS['2O'].primitives['U_mul'], 56 Toffolis = 392 ... |
| `hop_fermion_2033` | wilson | Stated-not-derived | app06 2033 derivation | the overlap kernel H_W is the Wilson operator |
| `hop_mcx` | mbu | Cited | groups.mcx_toffolis (r17 core) | C^nX ladders at n-1 Toffolis: the n-2 AND temporaries uncomputed by measurement (the report's rule) |
| `hop_n_spin` | 2 | Cited | ruling E21 (2), H. Lamm 2026-10-01 (author-confirmed); FermionPrimitives_unpub section_hamiltonian.tex:77-80 | spinor components that hop after the spinor frame, n_D/2 = 2 at d=3, r=1. Sets the hop norm N_c n_spin = 4 of the worst-case Trotter bound, and the se ... |
| `clifford_2O` | clifford | Cited | standard group theory (E23 (1), H. Lamm 2026-10-01: 'a math fact') | 2O is the lift to SU(2) of the octahedral rotation group O = single-qubit Clifford group mod phases (24 elements); every 2O element normalizes the Pau ... |
| `t_evol_fm` | 2 | Assumed | - | r24 (H. Lamm 2026-10-02, 'ok to run for a shorter amount of time'): app06 't = 2 fm = 10a'; was 10 fm = 50a (r17-r23) |
| `a_fm` | 0.2 | Assumed | - | app06 2033 box 'a ~ 0.2 fm' |
| `eps_trotter` | 0.1 | Assumed | - | app06 'at eps ~ 0.1' |
| `trotter_rule_2033` | state | Cited | ruling E26 (ii), H. Lamm 2026-10-01 ('i dont think we can do such a wildly off result'); AlvesLammLiu_inprep (FERMILAB-PUB-26-0397-T), evaluation DERIVED HERE | N_Trotter from the state-dependent second-order estimate in the thermal state (trotter_state_dependent_steps): norms of the nested commutators replace ... |
| `temp_lat_2033` | 0.5 | Assumed | - | r22: temperature in lattice units, T a, for the thermal-state variances. The chapter says T ~ T_c only; 0.15-1 moves the central count over 5008-5342 |
| `g2_bare` | 1.0 | Assumed | - | r21 (3b): bare gauge coupling for the WORST-CASE bound. The chapter states none; an order-one coupling. The state-dependent estimate does not depend o ... |
| `casimir_max_2O` | 3.75 | Assumed | - | r21 (3b): largest SU(2) Casimir j(j+1) kept by the 2O electric term, j = 3/2; enters the worst-case bound only |
| `dt_empirical_fm` | 0.01 | Cited | app03 2033 box ('2nd-order Trotter at Delta t = a_t/10', a_t = 0.1 fm/c) | the fixed step of Ch. 5's 3^3 evolution, the printed comparison |
| `beta_H_rule` | normalization | Assumed | - | referee C1 (2026-10-04), DERIVED HERE: the filter e^{-beta H/2} is a polynomial in H / 2Q, so beta\|\|H\|\| = 2Q / (T a) = 864 ('normalization'); 'stated' ... |
| `delta_beta` | 1e-04 | Assumed | - | app06 'and delta_beta ~ 1e-4' (ch09-delta-beta); sqrt(1e2 ln 1e4) = 30.35 -> d_beta = 30; 1e-3 would give 26.3 |
| `aa_rounds` | 8 | Stated-not-derived | app06 2033 derivation, route (iv) | '~8 rounds of fixed-point amplitude amplification'; 'an estimate rather than a derivation' |
| `aa_call_rule` | 2k+1 | Assumed | - | verifier 2026-10-04: k rounds of amplitude amplification call the filter 2k + 1 times (U and U^dag once per round, plus the first call); 'per_round' k ... |
| `kappa_physical` | (1000.0, 10000.0) | Stated-not-derived | app06 2033 derivation | 'kappa ~ 1e3-1e4 ... near-physical-mass' |
| `haar_deficit_frac` | 0.05 | Assumed | - | app06 route (iv) 'even a 5% deficit' |
| `haar_register_qubits` | 1000.0 | Assumed | - | app06 route (iv) 'on a 1e3-qubit register' |
| `gibbs_terms` | 300 | Stated-not-derived | app06 route (iv) | 'approximately 300 local jump operators' |
| `gibbs_sweeps` | 10 | Stated-not-derived | app06 route (iv) | 'an optimistic ten sweeps' |
| `gibbs_n_evol` | (1, 4) | Assumed | - | E23 (3), r20: controlled evolutions of length T_OFT per jump; low 1 (single evolution), high 4 (OFT e^{iHt} A e^{-iHt}, its uncompute, and the coheren ... |
| `gibbs_n_evol_central` | 2 | Assumed | - | E23 (3), r20: the OFT's forward and backward evolutions |
| `gibbs_eps_oft` | 0.001 | Assumed | - | E23 (3), r20: Gaussian time window truncated at 1e-3, T_OFT = beta sqrt(ln 1/eps) (central and high; the low end takes T_OFT = beta) |
| `gibbs_sweeps_high` | 30 | Assumed | - | E23 (3), r20: the chapter calls ten sweeps optimistic; high end 30 |
| `gibbs_accept_T` | 1000.0 | Assumed | - | E23 (3), r20: Metropolis weight on the frequency register and the local jump LCU, an upper allowance per jump; negligible next to the evolution |
| `sigma2_ncs` | 100.0 | Stated-not-derived | app06 shot-count paragraph | 'sigma^2_NCS ~ 1e2 at the longest t'. Chapter-open pass (referee C3, 2026-10-05): an IDEAL-readout assumption; the Hadamard-test readout of the Kubo f ... |
| `ncs_e_max` | 1.5 | Assumed | - | C3: \|\|E^a\|\| = j_max on the 2O irreps restricting SU(2) j <= 3/2 (the chapter's 'Casimir cut at j = 3/2'); E^a = 0 on the other four irreps |
| `ncs_b_max` | 2.0 | Assumed | - | C3: B^a = -i Tr(sigma^a U_clover), \|B^a\| <= 2 |
| `tc_ew_gev` | 159.0 | Cited | DOnofrio_Rummukainen_Tranberg_2014 | abstract: 'cross-over ... T_c = (159 +- 1) GeV' |
| `mr_rate_box3` | 0.0023 | Cited | Moore_Rummukainen_2000 | Table (Vol_table), a = 1/(2 g^2 T) (beta = 8): L g^2 T = 3.0 -> Gamma/alpha^4 T^4 = 0.0023 +- 0.0016 |
| `mr_rate_box10` | 1.68 | Cited | Moore_Rummukainen_2000 | same table: L g^2 T = 10 -> 1.68 +- 0.03; 'large volume behavior is obtained by L = 8/g^2T' |
| `mr_rate_box3_err` | 0.0016 | Cited | Moore_Rummukainen_2000 | same table: the +-0.0016 on 0.0023 at L g^2 T = 3.0 (suppression 431-2400, quoted ~1e3) |
| `shots_2033` | 10000.0 | Stated-not-derived | app06 2033 box | '1e4 (psi-bar psi: 1e3/coupling; Gamma_sph: 1e4/coupling)' |
| `campaign_yr` | 5 | Stated-not-derived | app06 requirements | 'Campaign horizon 5 years' |
| `lever_eps_obs` | 0.15 | Assumed | - | r24 lever: a 15% statistical target on Gamma_sph (shots sigma^2 / eps^2) |
| `yukawa_placement` | inside | Stated-not-derived | author ruling r25 (a), H. Lamm 2026-10-02, relaying H. Singh | 'the Ginsparg-Wilson-projected Yukawa is INSIDE': the Yukawa couples the projected fields, whose projector (1 - sgn(H_W))/2 is a second QSVT polynomia ... |
| `depth_variant_2033` | fanout | Assumed | - | ruling R10 (H. Lamm 2026-10-02): where F* < 10 the box quotes the workspace fix; factory analysis: a CNOT fan-out of the shared control on 5 extra LQ ... |
| `anc_fanout_2033` | 5 | Assumed | - | r25 R10, factory analysis: ~5 LQ for the CNOT fan-out of the leaf flag (5 copy Toffolis), the orientation and Dirac flags and the mux ANDs |
| `t_short_fm` | 1.0 | Assumed | - | ruling r25 (e): the first result reads one temperature at 5a and 10a; 5a = 1 fm, priced at its own evolution depth (ruling R4) |
| `n_temperatures` | (4, 8) | Stated-not-derived | author ruling r25 (e), H. Lamm 2026-10-02 | 'campaign 4-8 temperatures (not 50 couplings)' on the symmetric / near-T_c side; the two-time linearity check at one of them, 10a only at the rest |

#### 2033 target (3D 2O overlap + Higgs): printed figures (50)

| name | value | provenance | source | note |
|---|---|---|---|---|
| `d_sgn_stated` | 208 | Stated-not-derived | app06 2033 derivation and box | 'd_sgn = ceil(kappa ln delta^-1) = 208' at kappa = 30 (r24; r17-r23 691, printed ~700) |
| `higgs_strings_per_site_stated` | 40 | Stated-not-derived | app06 2033 derivation | '$40$ select strings per site' of the Higgs vertex (8 on-site Z strings x the 5-term binary expansion of rho) |
| `sp_query_toffoli_stated` | 598 | Stated-not-derived | app06 2033 derivation | '561 of the query's 598 Toffolis'; the query is 598 Toffolis, 24 T and 6 rotations |
| `c_sp_stated` | 4400.0 | Stated-not-derived | app06 2033 derivation | '598 x 7 + 24 + 6 x 29.1 ~ 4.4e3 T' per query (4384.8) |
| `t_per_dov_stated` | 1.8e+06 | Stated-not-derived | app06 2033 derivation | '417 queries per application and 1.8e6 T' (r25 (a), Yukawa inside: 1,847,747.1; r24 9.4e5; r22 3.1e6) |
| `n_trotter_2033_stated` | 450 | Stated-not-derived | app06 2033 derivation and box | 'N_Trotter ~ 450 steps of dt ~ 0.022 a' (r24 at t = 10a; r22 5031 at 50a) |
| `n_trotter_range_stated` | (110.0, 1300.0) | Stated-not-derived | app06 2033 derivation | 'a range of 1.1e2-1.3e3 steps' (r24: 107, 1274) (1186 fastest-mode phase; 14,242 with the mean kept) |
| `lambda_sd_stated` | 20 | Stated-not-derived | app06 2033 derivation | 'Lambda_sd ~ 20 a^-3' at T a = 0.5 (20.25) |
| `n_modes_stated` | 156 | Stated-not-derived | app06 2033 derivation | '156 transverse oscillators (three colors, two polarizations, 26 nonzero momenta)' |
| `n_trotter_bound_stated` | 7400.0 | Stated-not-derived | app06 2033 derivation | '7.4e3 steps' of the worst-case bound (r24: 7442 at t = 10a; r21-r23 83,195 at 50a) |
| `lambda_2033_stated` | 5500.0 | Stated-not-derived | app06 2033 derivation | 'Lambda ~ 5.5e3 a^-3' of the worst-case bound (5537 at g^2 = 1) |
| `t_shot_empirical_step_stated` | 1.5e+10 | Stated-not-derived | app06 2033 derivation | '200 steps ... the shot would be 1.4e10 T' at the fixed step of Ch. 5 (17 filter calls: 1.4989e10; C1: 1.3508e10; r25: 1.2635e10; r24 6.4e9; 33 per st ... |
| `n_dov_per_step_stated` | 19 | Stated-not-derived | app06 2033 derivation | '19 applications per step': the Jacobi-Anger order at 2Q dt = 9.6 (r24; r22 13 at 4.3) |
| `n_dov_per_step_worst_stated` | 6 | Stated-not-derived | app06 2033 derivation | '7.4e3 steps of 6 applications each' at the worst-case step (r24) |
| `n_dov_evolution_stated` | 8600.0 | Stated-not-derived | app06 2033 derivation | 'The evolution is then 8.6e3 applications' (r24: 8550; r22 65,403) |
| `n_dov_floor_stated` | 4300.0 | Stated-not-derived | app06 2033 derivation | 'no step size takes the evolution below 2Qt = 4320 applications, 4.1e9 T' (r24; r22 21,600) |
| `t_evolution_stated` | 1.6e+10 | Stated-not-derived | app06 2033 derivation | 'or 1.6e10 T-gates' (r25: 1.5798e10; r24 8.0e9) |
| `t_shot_range_stated` | (1.3e+10, 2.9e+10) | Stated-not-derived | app06 2033 derivation | 'Over the range above it is 1.1-2.6e10 T' (17 filter calls: 1.3469e10, 2.8708e10; C1 1.2-2.7e10; r25: 1.1115e10, 2.6352e10; r24 5.6e9-1.3e10) |
| `t_shot_bound_stated` | 8.5e+10 | Stated-not-derived | app06 2033 derivation | 'at the worst-case step ... 8.3e10 T' (17 filter calls: 8.546e10; C1 8.4e10; r25: 8.310e10; r24 4.2e10; r22 1.276e12) |
| `d_beta_stated` | 89 | Stated-not-derived | app06 2033 derivation | referee C1 (2026-10-04): 'degree d_beta = 89' at beta\|\|H\|\| = 2Q / (T a) = 864 (89.2); r25 '~30' at the record 1e2 |
| `n_dov_tpq_stated` | 1500.0 | Stated-not-derived | app06 2033 derivation | '1.5e3 applications' (89 x 17 = 1513: 2k + 1 filter calls in k = 8 rounds, verifier 2026-10-04; first referee pass 712 = 89 x 8; r25 240) |
| `t_tpq_stated` | 2.8e+09 | Stated-not-derived | app06 2033 derivation | '2.8e9 T' inside the sphaleron shot (2.796e9 at 17 calls; first referee pass 1.3e9; r25: 4.4346e8; r24 2.3e8; r22 7.3e8) |
| `t_per_shot_2033_stated` | 1.9e+10 | Stated-not-derived | app06 2033 derivation and box | 'per-shot total over the N_Dov = 8.8e3 applications is 1.6e10 T-gates' (17 filter calls: 1.8597e10; C1 1.7115e10; r25: 1.6242e10; r24 8.2e9; r22 2.008 ... |
| `ratio_representative_stated` | 19 | Stated-not-derived | app06 2033 derivation | 'a factor 19 above the budget' (17 filter calls: 18.60; C1 2026-10-04: 17.11; r25: 16.24; r24 8.2; r22 2.0e2) |
| `gibbs_T_stated` | (4.8e+12, 1.5e+14) | Stated-not-derived | app06 route (iv) | 'The range is 5.5e11-1.7e13' (r25, re-priced at the Yukawa-inside application; r24 2.8e11-8.9e12; r22 9.2e11-2.9e13) |
| `gibbs_T_central_stated` | 2.5e+13 | Stated-not-derived | app06 route (iv) | '2.5e13 T' (C1 2026-10-04, beta\|\|H\|\| = 864: 2.518e13; r25: 2.914e12; r24 1.5e12; r22 4.8e12) |
| `gibbs_over_tpq_stated` | 9000.0 | Stated-not-derived | app06 route (iv) | '9.0e3 times' (17 filter calls: 9005; C1 2026-10-04: 1.914e4; r25 6.6e3; both scale with the application price) |
| `t_physical_stated` | (6.1e+11, 6.1e+12) | Stated-not-derived | app06 2033 derivation | '6.1e11-6.1e12 T' (17 filter calls: 6.107e11, 6.120e12; C1: 5.620e11, 5.632e12; r25: 5.333e11, 5.345e12; r24 2.7e11-2.7e12; r22 2.0e12-2.0e13) |
| `eps_l_2033_stated` | 5.4e-12 | Stated-not-derived | app06 2033 derivation and box | 'eps_l <~ 5.4e-12' (17 filter calls: 5.377e-12; C1: 5.843e-12; r25 '6.2e-12') at 0.1 faults per shot for the 10a sphaleron shot (R3; r25 6.157e-12; r2 ... |
| `wall_per_shot_2033_stated` | 19000.0 | Stated-not-derived | app06 requirements and wall-time paragraph | '1.9e4 s' per 10a shot (18,597; C1 17,115) (r25: 16,241.7; r24 8.2e3; r22 2.008e5) |
| `t_shot_condensate_stated` | 2.8e+09 | Stated-not-derived | app06 2033 box and wall-time paragraph | 'condensate shot 2.8e9 T' (TPQ alone at its own R-TOL, 17 filter calls: 2.7901e9; C1 1.3119e9; r25 4.4166e8, r24 2.2e8) |
| `eps_l_condensate_stated` | 3.6e-11 | Stated-not-derived | app06 2033 box | 'eps_l <~ 3.6e-11' for the condensate (17 filter calls: 3.584e-11; C1: 7.623e-11; r25 2.264e-10; r24 4.5e-10) |
| `wall_day_condensate_point_stated` | 32 | Stated-not-derived | app06 wall-time paragraph | '32 days per temperature' (1e3 condensate shots; 32.29; C1 15.18; r25 5.112; r24 2.6 days per coupling point) |
| `t_shot_kappa_r22_stated` | 6.1e+10 | Stated-not-derived | app06 2033 derivation | 'at kappa ~ 1e2 the same circuit costs 5.3e10 T' (17 filter calls: 6.119e10; C1 5.631e10; r25: 5.344e10; r24 2.7e10) |
| `lever_low_step_stated` | 1.4 | Stated-not-derived | app06 wall-time paragraph | 'the low end of the step range (1.5)' (r24) |
| `lever_eps15_stated` | 2.25 | Stated-not-derived | app06 wall-time paragraph | 'a 15% statistical target (2.25)' (r24) |
| `wall_yr_lever_low_step_stated` | 1.8 | Stated-not-derived | app06 wall-time paragraph | '1.8 yr' at the low end (r24) |
| `wall_yr_lever_eps15_stated` | 1.2 | Stated-not-derived | app06 wall-time paragraph | '1.2 yr' at 15% (r24) |
| `wall_yr_lever_both_stated` | 0.8 | Stated-not-derived | app06 wall-time paragraph | 'together 0.8 yr' (r24) |
| `t_shot_short_stated` | 9.6e+09 | Stated-not-derived | app06 wall-time paragraph and box | '160 steps of 23 applications ... 9.6e9 T' at 5a (17 filter calls: 9.5897e9; C1 8.109e9; r25: 7.2366e9) |
| `eps_l_short_stated` | 1e-11 | Stated-not-derived | app06 2033 box | '1.0e-11 at 5a' (17 filter calls: 1.043e-11; C1: 1.233e-11; r25: 1.382e-11) |
| `wall_s_short_stated` | 9600.0 | Stated-not-derived | app06 wall-time paragraph | 'or 9.6e3 s' per 5a shot (9589.7; C1 8109; r25 7236.6) |
| `wall_s_condensate_stated` | 2800.0 | Stated-not-derived | app06 wall-time paragraph | '2.8e3 s' per condensate shot (17 filter calls: 2790.1; C1: 1311.9; r25 441.7) |
| `f_star_register_stated` | (5.6, 6.3) | Stated-not-derived | app06 wall-time paragraph | 'only 5.6-6.3 T gates per reaction layer' as itemized (10a shot: 5.59-6.32) |
| `slowdown_register_stated` | (1.6, 1.8) | Stated-not-derived | app06 wall-time paragraph | '1.6-1.8 times slower' (10 / F*) |
| `f_star_2033_stated` | 12 | Stated-not-derived | app06 wall-time paragraph | 'raises this to 12' (12.17-12.54 over the runs) |
| `shots_first_stated` | (1100.0, 5500.0) | Stated-not-derived | app06 box and wall-time paragraph | '1.1e3 shots at 10a and 5.5e3 at 5a' (1111, 5497; r25 4444 under linear growth) |
| `wall_yr_first_stated` | 2.3 | Stated-not-derived | app06 box and wall-time paragraph | first result (method test) '2.3 yr' (2.325; before the 2026-10-05 ruling 2.005 with the 5a ratio 4; C1 1.744; r25 1.591) |
| `wall_yr_campaign_stated` | (2.7, 3.0) | Stated-not-derived | app06 box and wall-time paragraph | campaign '2.7-3.0 yr' (method test + condensate at 4-8 temperatures: 2.679-3.032; the rate campaign before the 2026-10-05 ruling 36.08-60.00, record w ... |
| `campaign_over_horizon_stated` | (0.54, 0.61) | Stated-not-derived | app06 box and wall-time paragraph | 'inside the 5-year horizon' (0.536-0.606; rate campaign before 2026-10-05: 7-12x) |

#### 2033 target (3D 2O overlap + Higgs): comparison records (38)

| name | value | provenance | source | note |
|---|---|---|---|---|
| `ancilla_2033_asserted` | 110 | Stated-not-derived | app06 (pre-r20) | RECORD (E23 (5)): '+ ~110 ancilla', asserted; replaced by the itemized peak (43 at r20, 33 at r21, 38 at r22 in the single-particle architecture) |
| `anc_unary_control` | 1 | Assumed | - | RECORD (second-quantized construction, r20): the control qubit of SELECT's unary iteration (Babbush_PRX_2018); the AND ladder holds ceil(log2 L) - 1 t ... |
| `anc_link_pq` | 6 | Cited | FermionPrimitives_unpub section_su2_diag.tex:188 (BO row, 'Ancilla 6') | RECORD (second-quantized construction, draft frame): the p, q angle registers of the 2O color squish, held through V, hop, V^dag (share='link') |
| `anc_link_hop_register` | 3 | Cited | FermionPrimitives_unpub section_hopping.tex:154 | RECORD (second-quantized construction, draft frame): 'at most a 3-qubit ancilla register for the eigenvalue computations', held with p, q |
| `anc_link_parity` | 1 | Cited | FermionPrimitives_unpub section_su2_diag.tex:206 | RECORD (second-quantized construction, draft frame): 'compute the parity of the p,q register on to another clean ancilla' |
| `anc_spinor` | 2 | Cited | FermionPrimitives_unpub section_spin_diag.tex:280 | RECORD (second-quantized construction): 'requires 2 spare clean ancilla' (the d=3 spinor frame) |
| `kappa_2033_r22` | 100.0 | Assumed | - | RECORD (r17-r23 headline kappa ~ 1e2): prices the retired second-quantized construction and the r22 headline (2.008e11), so every earlier print stays ... |
| `lcu_terms_stated_pre_r21` | 1700.0 | Stated-not-derived | app06 (pre-r21) | RECORD: '~1.7e3-term LCU', an absolute at V=27. 8Q = 1728 is exactly the 4Q + 4Q selected strings of U_L and U_R in arxiv_2607_28524: the old line was ... |
| `lcu_terms_per_site_stated` | 72 | Stated-not-derived | app06 (r21) | RECORD (second-quantized construction): '72 Pauli strings per site': 24 hop + 8 on-site + 40 Higgs |
| `lcu_terms_stated` | 1900.0 | Stated-not-derived | app06 (r21) | RECORD (second-quantized construction): '1944 at V = 27' |
| `link_frame` | multiplexer | Cited | r21 ruling (1), H. Lamm 2026-10-01 ('promote'); DERIVED HERE r20 (colour_multiplexer_2O; arxiv_2312_10285 ordered-product encoding) | RECORD (second-quantized construction): the link-dressed Fock-space hop as the 2O color multiplexer W and W^dag on the far site plus the draft's spin ... |
| `hop_frame_undo` | True | Cited | ruling E21 (1), H. Lamm 2026-10-01 (apply_log/r19_core.md); FermionPrimitives_unpub section_hamiltonian.tex:86-91, su2_diag.tex:238, section_resources.tex:15 | RECORD (second-quantized, draft frame). OUR ESTIMATE, ruled in: the color frame V_g is applied AND undone on both sites of every link in every BE que ... |
| `hop_share` | link | Cited | ruling E21 (1), H. Lamm 2026-10-01 ('make the switch'; apply_log/r19_core.md) | RECORD (second-quantized, draft frame): the color squish and parity flags are computed once per link and held through V(x), V(y), hop, V^dag(x), V^da ... |
| `hop_phasing_be` | none | Stated-not-derived | app06 (r21) | RECORD (second-quantized, draft frame): a block-encoding query carries no Trotter phasing |
| `w2_cs_multiplicity` | 4 | Stated-not-derived | FermionPrimitives_unpub section_resources.tex:50-53 | RECORD (the r17-r20 band top): the draft's 3d/4d Wilson hop 4 C^S + 148 N per link |
| `w2_spinor_t_per_colour` | 148 | Stated-not-derived | FermionPrimitives_unpub section_resources.tex:50-53 | RECORD: '+148N' |
| `n_rot_be` | 25 | Stated-not-derived | app06 (r21) | RECORD (second-quantized construction): '~25 synthesized rotations in PREPARE and the QSP phases' per query |
| `c_be_stated` | 63000.0 | Stated-not-derived | app06 (r21) | RECORD (second-quantized construction): 'c_be ~ 6.3e4 T-gates' (62,896) |
| `c_be_draft_stated` | 2.4e+06 | Stated-not-derived | app06 (r21) | RECORD: 'c_be ~ 2.4e6' with the draft's diagonalizing frame |
| `t_shot_draft_stated` | 7e+14 | Stated-not-derived | app06 (r21) | RECORD: '7.0e14 T per shot' with the draft's frame (7.007e14) |
| `draft_over_mux_stated` | 39 | Stated-not-derived | app06 (r21) | RECORD: '39 times the estimate we quote' (38.7) |
| `draft_over_mux_link_stated` | 64 | Stated-not-derived | app06 (r21) | RECORD: '64 times the multiplexer's 468 T' (63.6) |
| `t_per_dov_sq_stated` | 4.3e+07 | Stated-not-derived | app06 (r21) | RECORD (second-quantized construction): '691 x 6.3e4 ~ 4.3e7 T' |
| `t_per_shot_sq_stated` | 1.8e+13 | Stated-not-derived | app06 (r21 box) | RECORD (second-quantized construction at the worst-case step): '~1.8e13' (1.809e13), 1005 LQ |
| `t_evol_fm_r22` | 10 | Assumed | - | RECORD (r17-r23 evolution time 10 fm = 50a): prices the retired records and the r22 headline |
| `n_trotter_2033_pre_r21` | 30 | Stated-not-derived | app06 (pre-r21) | RECORD: 'N_Trotter ~ 30 steps', dt = 1.7 a. With a valid block encoding it would need 752 applications per step |
| `n_dov_per_step_bound_stated` | 5 | Stated-not-derived | app06 (r22-r23) | RECORD: 'with 5 applications per step' at the r22 worst-case step (derived there too); also the r21 '~5 D_ov applications per step' that the second-qu ... |
| `beta_H` | 100.0 | Assumed | - | RECORD (r25 and earlier): app06 'beta\|\|H\|\| ~ 1e2'. Referee C1 (2026-10-04): inconsistent with the chapter's own block encoding, whose normalization is ... |
| `haar_rounds_stated` | 3e+07 | Stated-not-derived | app06 route (iv), RECORD | '~3e7 amplification rounds' (ch09-haar-rounds; was ~1e7); retired from the prose 2026-10-05 (typical_reference_tpq: 7.0e29 at 3^3) |
| `gibbs_sampler_T` | (1e+10, 1e+12) | Stated-not-derived | app06 (pre-r20) | RETIRED from the prose (E23 (3)): the unrepriced '1e10 ... toward 1e12'; record only |
| `wall_yr_per_point_stated` | 2.6 | Stated-not-derived | RECORD (r24 box) | '2.6 yr per point' for Gamma_sph (r24: 1e4 x 8241 s = 2.61 yr; r23 64); r25: 5.15 yr per temperature at 10a, inside the campaign |
| `wall_yr_serial_stated` | 130.0 | Stated-not-derived | RECORD (r24 wall-time paragraph) | '1.3e2 yr for the ~50-point scan' (r24: 130.6; r23 3.2e3); retired by r25 (e) |
| `n_coupling` | 50 | Stated-not-derived | RECORD (r24 Eq. (Nshot), box, requirements) | 'N_coupling ~ 50'; r25 (e): 4-8 temperatures |
| `reduction_scan_stated` | 26 | Stated-not-derived | RECORD (r24 wall-time paragraph) | 'the scan needs its cost cut by 26' (r24: 130.6 / 5 = 26.1; r23 6.4e2) |
| `t_shot_to_fit_scan_stated` | 3.2e+08 | Stated-not-derived | RECORD (r24 wall-time paragraph) | 'a shot of 3.2e8 T' (3.156e8; unchanged: set by the shot count) |
| `wall_yr_condensate_scan_stated` | 0.36 | Stated-not-derived | RECORD (r24 box) | '0.36 yr for the ~50-point scan'; retired by r25 (e) |
| `ncs_var_ratio_short` | 4 | Assumed | - | RECORD (shot audit, r25): under linear growth the 5a signal is half the 10a signal, so at the same absolute noise its relative variance is 4x. Replace ... |
| `wall_yr_campaign_r25_stated` | (36, 60) | Stated-not-derived | RECORD (app06 before the 2026-10-05 method-test ruling) | '36-60 yr for 4-8 temperatures, 7-12x the 5-year horizon' (rate campaign; at 30%: '4.0-6.7 yr') |

## 4. How the numbers are built

### 4.1 Rules both instances use

- **Rotations.** Each circuit sets its own per-rotation tolerance from the number N_rot of synthesized rotations in one shot, ε_rot = sqrt(ε_syn / N_rot) with ε_syn = 10⁻² (`common.eps_rot_for`), and prices each rotation at the repeat-until-success (RUS) fit 1.15 log2(1/ε_rot) + 9.2 T (`common.t_per_rotation`). The tolerance depends on the count and the count does not depend on the tolerance, so there is no iteration (`_step_rtol_2028`, lines 1544-1555, asserts this).
- **Toffolis** cost 7 T. **Logical error rate** ε_l = 0.1 / (T per shot).
- **Wall time** on one machine: shots × (T per shot × 1 µs + 0.1 ms). Where a circuit cannot feed ten factories (F* = N_T / D_T < 10, with D_T the T-depth and a 10 µs reaction time per layer), each shot is charged the longer of N_T × 1 µs and D_T × 10 µs. The depth exports come from `common.depth_exports`.

### 4.2 2028: Trotter step under Hamming-weight phasing (`_model_2028`, lines 1558-2002)

**Qubits** (lines 1560-1565, 1617-1638). Fermions 2 Dirac components × 3 colors × V = 8 × L5 = 4 = 192; gauge D × V × 2 = 16; ancilla 42; total 250.

**Terms and grouping** (`_hwp_groups_2028`, lines 1452-1496). The Pauli rotations of one first-order step are grouped into sets of k equal-angle rotations. Hamming-weight phasing (`hwp_synth_rotations`, `hwp_toffolis`, `hwp_ancilla`, lines 332-349) does a group with floor(log2 k) + 1 synthesized rotations and k − w(k) adder Toffolis, where w(k) is the binary weight of k, and holds k − w(k) ancilla.

| term | rotations per step | groups | synthesized per step | Toffolis per step |
|---|---|---|---|---|
| fifth-direction hops: open boundary, V(L5 − 1) = 24 hops × 3 colors × 2 strings | 144 | 2 × (k = 72) | 14 | 140 |
| spatial hops: V L5 = 32 hops × 3 colors × 8 strings per copy on a Z3 link (`groups.rhop_link`); the 12 color and s copies share one link | 768 | 64 × (k = 12) | 256 | 640 |
| electric: 8 links × 3 diagonal Z strings | 24 | 3 × (k = 8) | 12 | 21 |
| on-site m + r(1 + d) of the domain-wall kernel: 2N Z rotations per site on 32 sites | 192 | 1 × (k = 192) | 8 | 190 |
| total | 1128 | | 290 | 991 |

**Step and shot** (`_step_cost_2028`, lines 1499-1541). One shot is 21 steps (20 ramp + 1 probe), so N_rot = 290 × 21 = 6090, ε_rot = 1.2814e-3 and 20.249 T per rotation. One step is 290 × 20.249 + 991 × 7 = 12,809.3 T and the shot 268,994.9 T. The spatial hops are 75% of it (9,663.8 T per step).

**Chunking.** The k = 192 group needs 190 adder ancilla, a register of 399 LQ. Splitting every group to k ≤ 32 (`hwp_chunk_k`; fifth-direction hops 4 × 32 + 2 × 8, on-site 6 × 32) gives 336 synthesized rotations and 985 Toffolis per step; N_rot = 7056, ε_rot = 1.1905e-3, 20.371 T per rotation, 13,739.8 T per step, **288,535.5 T per shot**. This is the circuit the 250-LQ register runs. `hard_ops` is the pair (268,994.9, 288,535.5) and `t_per_shot` is the chunked value.

**Ancilla and depth** (`hwp_step_depth_2028`, lines 847-867). Groups run one slot at a time. Inside a slot every weight-bit rotation runs at once on its own repeat-until-success ancilla, and three spatial-hop groups on mutually non-adjacent links share a slot. The peak is a spatial-hop slot, 3 × (10 adder + 4 RUS) = 42 ancilla. A slot costs one rotation layer plus the adder depth, ceil(2 log2 k) for a carry-save tree (low) or k − w(k) (high). One step is 37 slots, T-depth 1059.7-1318.7; one shot 22,255-27,694; F* = 10.42-12.97.

**Shots and wall** (lines 1712-1759). The chapter's σ² ≈ 10 is the condensate averaged over the 27 sites of the 2033 lattice. The 2028 estimator averages 3 colors × 8 sites = 24 copies, taken as independent: σ² = 10 × 27 / 24 = 11.25, so 1125 shots at 10% and 125 at 30%. If the sites of a color were fully correlated, σ² = 90 and 9000 shots. First result: 125 × (288,535.5 × 1 µs + 0.1 ms) = 36.1 s. Campaign: L5 = 2, 3, 4 at 1125 shots each, each on its own register (149, 197, 250 LQ) at its own tolerance (1.586e5, 2.391e5, 2.885e5 T). L5 = 2 has F* below 10 and is charged its depth-limited wall. Total 794.2 s. ε_l = 0.1 / T = 3.47e-7 (chunked) to 3.72e-7.

**Why not the overlap operator in 1+1D** (lines 1663-1683). The same single-particle architecture as 2033, built for the 1+1D Z3 kernel (`single_particle_query_Z3_1p1d`, `lift_cost`): one query is 40 Toffolis and 6 rotations, 391.6 T. At κ = 30 and δ = 10⁻² (κ the condition number α/Δ of the Wilson kernel, Δ its spectral gap), one application of the overlap Hamiltonian is 140 queries plus the lift, 5.62e4 T, and a ground-state preparation of 10 applications is 5.78e5 T, 2.0 times the domain-wall shot and 5.8 times the 10⁵ budget.

### 4.3 2033: overlap fermions by QSP, thermal state by a filter (`_model_2033`, lines 2009-2906)

**Qubits** (lines 2013-2018, 2252-2279). Fermions 4 × 1 × 2 × 27 = 216 (Q = 216 single-particle modes); gauge 3 × 27 × 6 = 486; Higgs 27 × (4 + 6) = 270; ancilla 43. The ancilla are itemized at the peak of the schedule:

- 5 held for the whole shot: readout control, QSP signal, amplification marker, and 2 linear-combination-of-unitaries (LCU) qubits over the parts of H;
- 14 held through one application: 2 lift-type qubits, the 9-qubit index register of color, spin and coordinates, the LCU qubit of h_ov, and the signal and real-part qubits of the sign polynomial;
- 5 for the term register of one query;
- 14 inside the link read: 3 direction flags, 6 unary-iteration temporaries, 5 for the copied link;
- 5 for a CNOT fan-out of the shared control, which lets the link read feed ten factories.

**Sign-function degree** (lines 2049-2061). d_sgn = ceil(κ ln(1/δ_sgn)) = ceil(30 × ln 10³) = 208, with κ = α/Δ and α = 2D + |D − m_0| = 7.5 the normalization of the H_W block encoding.

**One query of the single-particle Wilson kernel** (`single_particle_query_2O`, lines 498-539). An LCU of 4D + 1 = 13 unitaries (the on-site term, and a Dirac and a Wilson hop per direction and orientation) on a 5-qubit term register (PREPARE loads the term weights into that register; SELECT applies the chosen term). SELECT walks the sites by unary iteration (`Babbush_PRX_2018`), copies the addressed link, and applies the 2O color rotation to the color qubit of the index register. On a single qubit that multiplexer is 10 T and no Toffoli (`colour_multiplexer_index_2O`, lines 469-483; checked on all 48 elements in `test_2O_index_register_multiplexer_is_10_T`).

| item | Toffolis |
|---|---|
| direction flags (D − 1) | 2 |
| backward / forward / Dirac flags (D each) | 3 + 3 + 3 |
| controlled shifts mod 3 (2D × 2) | 12 |
| read pass, unary iteration D(L^D − 1) | 78 |
| link copy (b − 1) × D L^D, b = 6 | 405 |
| multiplexer ANDs on the orientation qubit 2(b − 1) | 10 |
| fix-up pass D(L^D − 1) | 78 |
| projector on the term register | 4 |
| total | 598 |

Plus 24 direct T (multiplexer and inverse 20, controlled-H in PREPARE 4) and 6 rotations (PREPARE and inverse 4, controlled QSP phase 2). The link read (read, copy, fix-up) is 561 of the 598 Toffolis and 90% of the query's T. A gate-level check of SELECT against the covariant Wilson kernel is `test_single_particle_select_reproduces_the_covariant_wilson_kernel`.

**One application** of ψ†h_ov ψ / 2Q (lines 2106-2128, `_price_sp` lines 2130-2153):
- 417 queries: 209 for sgn(H_W) in h_ov and 208 for the projector (1 − sgn(H_W))/2 of the projected Yukawa;
- the lift to Fock space (`lift_cost`, lines 555-563): two selected-Majorana unitaries over 2Q = 432 strings each, 2 × 431 = 862 Toffolis, 12 T and 6 rotations for the uniform preparation over the three coordinates mod 3;
- the reflection of the outer QSP on the 19 block-encoding ancilla: 18 Toffolis and 1 rotation;
- the Higgs vertex, once per application and outside the sign function: 27 group multiplications U_x at 392 T each (56 Toffolis, `GROUPS["2O"].primitives["U_mul"]`, arXiv:2312.10285) and 1080 SELECT Toffolis (40 strings per site: 8 on-site Z strings × the 5-term binary expansion of ρ).

**Step count** (`trotter_state_dependent_steps`, lines 606-641). Second-order Trotter steps. The nested-commutator norms are replaced by their standard deviations in the thermal state (Alves, Lamm, Liu, FERMILAB-PUB-26-0397-T). These are evaluated at weak coupling on the 156 transverse oscillators of 3^3 (3 colors × 2 polarizations × 26 momenta) at T a = 0.5, giving Λ_sd = 20.25 a⁻³ and N = ceil(sqrt(Λ_sd t³ / ε)) = 450 steps for t = 10a, ε = 0.1.

**Applications per step** (`jacobi_anger_order`, `dov_per_step`, lines 586-603). The Jacobi-Anger order at x = 2Q dt = 432 × 10 / 450 = 9.6 with truncation error 10⁻² / 450 per step: 19. Evolution: 450 × 19 = 8550 applications. This assumes a self-inverse block encoding; otherwise the count doubles.

**Thermal state** (lines 2089-2104). The filter e^(−βH/2) is a polynomial in H/2Q, so β‖H‖ = 2Q / (T a) = 864 and d_β = round(sqrt(864 ln 10⁴)) = 89 applications per call. Eight rounds of amplitude amplification call the filter 2 × 8 + 1 = 17 times: 1513 applications.

**Tolerance and shot.** Each application has 417 × 6 + 6 + 1 = 2509 rotations; over 8550 + 1513 = 10,063 applications that is 25,248,067, so ε_rot = 1.990e-5 and 27.159 T per rotation. Then

- c_sp = 598 × 7 + 24 + 6 × 27.159 = 4,373.0 T per query,
- application = 417 × 4,373.0 + 6,209.0 (lift) + 153.2 (reflection) + 18,144 (Higgs) = 1,848,028.6 T,
- **shot = 10,063 × 1,848,028.6 = 1.85967e10 T**, 18.6 times the 10⁹ budget; 85% evolution, 15% thermal state.

In the primitive breakdown the query Toffolis are 1.757e10 T (94%), query rotations 6.84e8, the Higgs U_x 1.07e8, direct T 1.01e8, Higgs SELECT 7.6e7, the lift 6.2e7.

**The other two shots.** At t = 5a the same estimate gives 160 steps of 23 applications plus the thermal state, 5193 applications and 9.590e9 T. The condensate is read from the prepared state, so its shot is the 1513 filter applications alone, at its own tolerance (25.587 T per rotation): 2.790e9 T.

**Depth** (`sp_query_depth_2O`, `application_depth_2O`, lines 803-844). T-depth counted in 10 µs reaction layers per item of the query, with the read and fix-up walks serial. With the fan-out the 10a shot has T-depth 1.516e9-1.529e9 and F* = 12.17-12.27; as itemized without it, F* = 5.58-6.32.

**Shots and wall** (lines 2475-2519). The readout variance σ²_NCS = 10² is assumed, so the 10a test takes 10² / 0.3² = 1111 shots. The 5a signal of the free lattice estimator is smaller, and at fixed absolute readout noise its relative variance is (S_10a / S_5a)² = 4.947 times larger (`ncs_lattice_background`, lines 702-772): 5497 shots. Wall per shot: 18,596.7 s (10a), 9,589.7 s (5a), 2,790.1 s (condensate). First result 7.338e7 s = 2.325 yr (3.156e7 s per year in this module). Campaign: the first result plus 10³ condensate shots at each of 4-8 temperatures, 8.454e7-9.570e7 s = 2.679-3.032 yr, 0.54-0.61 of the five-year horizon. ε_l = 5.38e-12 (10a), 1.04e-11 (5a), 3.58e-11 (condensate).

### 4.4 Sensitivities and comparisons the model also computes

All are in `--intermediates` and none enters a headline.

- Step-count range 107 steps (only the phase of the fastest gauge mode held to ε) to 1274 (mean kept): 1.347e10-2.871e10 T. Worst-case commutator bound (`trotter_bound_steps`, lines 430-462; Λ = 5537 a⁻³ at g² = 1): 7442 steps, 8.546e10 T. Ch. 5's fixed step: 1.499e10 T. No step count brings the evolution below 2Q t = 4320 applications, 7.98e9 T (`n_dov_floor`; `test_2033_trotter_rule_switch_and_records`).
- Not self-inverse: 3.44e10 T. κ read as ‖H_W‖/Δ: d_sgn = 346, 3.08e10 T. κ = 10²: 6.12e10 T. κ = 10³-10⁴: 6.11e11-6.12e12 T.
- Thermal route: a Gibbs sampler estimated at the same application price, 4.79e12 / 2.52e13 / 1.51e14 T (low / central / high), 9.0e3 times the filter at the center. A typical reference state succeeds with probability e^(−137.9) on the free fermion sector, which would need 7.0e29 amplification rounds (`typical_reference_tpq`).
- A second-quantized construction (the sign of the many-body operator ψ†H_Wψ rather than of the kernel, so not the chapter's Hamiltonian): 1.809e13 T per shot, recorded as `t_record_second_quantized`.
- Readout: the Hadamard test of the Kubo form needs relative variance 1.06e14 at 10a on the free-field background, so 1.17e15 shots at 30%, against the assumed 10².

## 5. What the paper prints

`resources.json` rows (the chapter's boxes and Table 1.1 read these):

| key | 2028 | 2033 |
|---|---|---|
| `lq` | (250, 250) | (1015, 1015) |
| `t` | (268,994.9, 288,535.5) | (1.85967e10, 1.85967e10) |
| `t_box` (the box as printed) | (2.7e5, 2.9e5) | (1.9e10, 1.9e10) |
| `shots` | 1125 | 1111.1 (10a, first result) |
| `f_star` | (10.42, 12.97) | (12.17, 12.27) |
| `t_condensate_arm` | | 2.790e9 |
| `t_short_time` | | 9.590e9 |
| `t_range_step_count` | | (1.347e10, 2.871e10) |
| `t_alternative_fixed_step` | | 1.499e10 |
| `t_alternative_worst_case_bound` | | 8.546e10 |
| `t_record_second_quantized` | | 1.809e13 |
| `wall_first_result_s` | | 7.338e7 |
| `wall_campaign_s` | | (8.454e7, 9.570e7) |

Other exports (`python -m estimates --chapter ch09`): 2028 ε_l 3.47e-7-3.72e-7, wall 36.1 s first result and 794.2 s campaign, T-depth per shot 22,255-27,694, campaign floor wall 653-771 s; 2033 ε_l 5.38e-12, T-depth per 10a shot 1.516e9-1.529e9, campaign floor wall 6.87e7-7.83e7 s, and 26.8-30.3 factories to run the campaign in one year.

Table 1.1 prints, for Ch. 9: 250 LQ and 2.7-2.9e5 T (2028); 1015 LQ and 1.9e10 T (2033). The 2028 box prints 250 = 192 + 16 + 42, 2.7-2.9e5 T, ε_l ≲ 3.5e-7, 1.3e2 and 1.1e3 shots, 36 s and 13 min. The 2033 box prints 1015 = 216 + 486 + 270 + 43, 2.8e9 T (condensate), 1.9e10 T (10a) and 9.6e9 T (5a), ε_l ≲ 3.6e-11 / 5.4e-12 / 1.0e-11, 1.1e3 + 5.5e3 shots, 2.3 yr and 2.7-3.0 yr.

Tests that pin them:

- `PUBLISHED` (lines 2939-2962) is the box as printed: `test_published_is_the_box_as_printed`; the model lies within 10% of it: `test_both_boxes_reproduce`, `test_2033_hard_ops_matches_box`.
- 2028: `test_2028_lq_components_exact`, `test_2028_step_arithmetic_by_hand`, `test_2028_chunked_variant_by_hand`, `test_2028_hard_ops_is_the_printed_range` (268,995 and 288,535 exactly), `test_2028_ancilla_itemized`, `test_2028_shots_and_wall_time`, `test_r25_2028_schedule_tiers_and_exports` (depth, F*, 36 s, 13 min).
- 2033: `test_2033_lq_components_exact`, `test_2033_single_particle_query_itemized`, `test_2033_one_application_components`, `test_2033_trotter_steps_state_dependent`, `test_2033_jacobi_anger_applications_per_step`, `test_2033_dov_evolution_tpq_total` (1.8596712e10 to 1e-6), `test_2033_eps_l_at_0p1_faults_per_shot`, `test_r24_condensate_arm_is_the_preparation_alone` (2.79010e9), `test_r25_short_time_shot_at_its_own_depth`, `test_r25_exports_2033`, `test_r25_two_tiers_2033` (2.325 yr, 2.679-3.032 yr), `test_2033_trotter_rule_switch_and_records` (the 8.0e9 floor).
- With the paper source present, `test_r25_chapter_prints_the_model` checks that the chapter's boxes and prose quote these numbers.

## 6. Circuit status

From `docs/CIRCUIT_STATUS.md`: 6 COMPILED, 12 SCALING, no CONJECTURE, no UNSOURCED.

- **COMPILED.** 2028: the Hamming-weight adders of all four term groups (k − w(k) Toffolis per group, arXiv:1709.06648, 1902.10673). 2033: the 598 Toffolis and 24 direct T of the single-particle query, compiled here; the tests simulate SELECT at gate level on the 3^3 lattice and compare sampled columns with the covariant Wilson kernel.
- **SCALING, 2028.** The synthesized rotations of the four groups (the RUS fit at the circuit's tolerance; the spatial-hop string count comes from the derived Z3 hop rule).
- **SCALING, 2033.** The query rotations; the lift's Toffolis, direct T and rotations (counts derived here from the construction of arXiv:2607.28524); the outer QSP reflection and its rotation; the Higgs U_x (compiled in arXiv:2312.10285, but no controlled variant is compiled, so one U_x is the price of a controlled 2O action); and the Higgs SELECT Toffolis.

## 7. Open items and limitations

- **Thermal reference state.** The 17 filter calls are a placeholder. A typical reference succeeds with probability about e^(−138) at 3^3 and β = 2a; an amplified random basis state is biased; a Gibbs-like product reference does better but without control. The 2033 thermal preparation is conditional.
- **Filter normalization.** The filter is priced on the fermion block encoding only. The electric, magnetic and Higgs terms inside e^(−βH/2) raise the normalization above 2Q = 432 and add their own cost; neither is priced. The weight of the projected Yukawa in the normalization is not priced either.
- **Chern-Simons readout.** The assumed σ²_NCS ≈ 10² needs a readout whose noise is set by ΔN_CS itself. A Hadamard test of the Kubo form against the free-field background has relative variance 1.06e14 at 10a (`ncs_hadamard_rel_var`, with the block-encoding normalization taken at its lower bound 9.23/a). The multiplicative renormalization of the rate operator is open.
- **Box size.** At g² = 1 and T a = 0.5 the box has L g² T = 1.5, where the classical sphaleron rate is suppressed about 730 times (431-2400 within the quoted errors) relative to large volume. No rate is claimed at 3^3; Higgs-mass matching and the location of T_c are open. The window t = 10a is marginal for the linear regime of Chern-Simons diffusion.
- **Step count.** 450 steps is the state-dependent estimate; the worst-case bound is 7442. Keeping the common mean of the nested commutators, which only shifts a global phase, gives 1274 (`n_trotter_range`). T a = 0.5 is a choice (the chapter says T ~ T_c; `trotter_state_dependent_steps(3, 3, 10.0, 0.1, T)` gives 448-478 steps for T a = 0.15-1), g² = 1 enters only the worst-case bound, and the Higgs modes are not in the oscillator sum. The chapter estimates that they would add 15–20% and that the fermion terms change the count by about 1%; the model computes neither.
- **Not priced.** The gauge and Higgs terms of each Trotter step (gauge about 2.0% of a step) and the amplification reflections (2.18e5 T per shot). A self-inverse block encoding is assumed; otherwise the shot is 3.44e10 T.
- **2028 ramp.** The diabatic error of the 20-step ramp is not estimated; the benchmark is meant to measure it. The step is taken first order (the chapter states no order). The quark mass m_f = 0 removes the wall-joining term; m_f ≠ 0 adds two groups of k = 24, +510 T per step (the `fifth_boundary` note), 4% of the unchunked step.
- **2028 shots.** 1125 shots treat the 24 site-color copies as independent (they share the link); full correlation within a color gives 9000.
- **Registers.** The 2028 register holds the chunked circuit; the unchunked one needs 399 LQ. F* ≥ 10 at 2028 needs three concurrent spatial-hop groups per slot, and at 2033 the five fan-out qubits; without them F* is 5.6-6.3 and every 2033 wall is 1.6-1.8 times longer.
- **Primitive counts.** The k − w(k) Toffoli count assumes adders uncomputed by measurement; a unitary uncompute doubles it. The adder depth ceil(2 log2 k) is a carry-save estimate. The 1+1D overlap comparison and the 2028 Z3 hop use string counts derived in this package, not compiled circuits.
- **Shared items** (`docs/OPEN_ITEMS.md`, last section): the fault budget counts T gates only; Clifford, idle, measurement and injection faults are excluded.

## 8. Conventions used

R-TOL is the 1% total synthesis error per shot and the per-circuit tolerance. The other conventions are Toffoli = 7 T, the fault budget of 0.1 expected faults per shot, and the one-machine wall-time model (1 µs per T, 0.1 ms per shot, ten factories at a 10 µs reaction time, depth correction where F* < 10). All are report-wide conventions described in the top-level `README.md` and implemented in `estimates/common.py`. This module uses 3.156e7 s per year (`SEC_PER_YR`).
