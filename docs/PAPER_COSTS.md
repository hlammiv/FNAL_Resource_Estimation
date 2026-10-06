# Gauge-group registers and gate costs, as the papers state them

Extracted 2026-09-28 from the LaTeX sources on disk, on the author's instruction to
use the published values rather than any proposal. Every number below is quoted
from the paper at the location given; nothing is computed except in the last
section, which says so. `log` means log2(1/eps) throughout.

Sources: **P36** = arXiv:2405.05973v4 `mainV2.tex` (Sigma(36x3)); **FFT** =
arXiv:2408.00075v2 `main.tex` (fast Fourier transforms); **P72** = arXiv:2511.17437v1
`main_s216.tex` (Sigma(72x3)); **BT** = arXiv:2208.12309 `apssamp.tex` (2T);
**BO** = arXiv:2312.10285v1 `main.tex` (2O). Line numbers are into those files.

## 1. Register per link — stated

| group | qubits | quoted statement |
|---|---|---|
| 2T (|G|=24) | **5** | "requiring five qubits or one quicosotetrit per gauge link" (BT abstract, :95); "cannot go below ... ceil(log2 24) = 5" (BT :136) |
| 2O (|G|=48) | **6** | "one additional qubit -- for a total of six -- per gauge link" (BO abstract, :104); 16 of the 64 states are forbidden (BO :153) |
| Sigma(36x3) (|G|=108) | **8** (or 3 qutrits + 2 qubits) | "presented for both an eight-qubit encoding and a heterogeneous three-qutrit plus two-qubit register" (P36 abstract, :128); "either 8 qubits ... or 3 qutrits and 2 qubits will be required to store the group register" (P36 :186-188) |
| Sigma(72x3) (|G|=216) | **9** (or 3 qutrits + 3 qubits) | "The elements of Sigma(72x3) can be encoded into a register consisting of 9 qubits" (P72 :150); "a Sigma(72x3)-register requires 9 qubits" (P72 :300-306) |

The whitepaper's 7 (Sigma(36x3)) and 8 (Sigma(72x3)) are ceil(log2|G|). **No paper
describes or builds circuits on that encoding**; the strings `\lceil` and `log_2`
never occur in a register-size statement. All T-counts below are for the
registers in this table. For 2T and 2O the whitepaper's 5 and 6 are correct.

## 2. Primitive gate costs — stated (Table `tab:tgatecost` in each paper)

T-gates / clean ancilla, in the paper's own convention (Sec. 4).

| primitive | 2T (BT :526-544) | 2O (BO :774-793) | Sigma(36x3) (P36 :622-642) | Sigma(72x3) (P72 :874-893) |
|---|---|---|---|---|
| U_inv | 28 / 0 | 112 / 1 | 119 / 4 | 637 / 4 |
| U_mult | 154 / 1 | **392 / 4** | 308 / 2 | 1120 / 5 |
| U_Tr | 12.65 log / 0 | 350 + 4.6 log / 2 (188 with 6 more ancilla) | 378 + 8.05 log / 7 | 1414 + 8.05 log / 12 |
| U_phi (electric phase) | (U_Tr plays this role) | (U_Tr plays this role) | 294.4 log / 0 | 294.4 log / 0 |
| U_F (naive Fourier) | 1150 log / 0 | 11370.1 log / 0 | 185898 log / 0 | 186295.4 log + 112 / 5 |
| U_FFT | 98 + 48.3 log / 2 (FFT :871-893) | 216 + 48.3 log / 4 (FFT) | 532 + 117.3 log / 8 (FFT :835) | **not constructed**; lower bound 515.2 log + 532 (P72 :905) |

Not stated anywhere: a controlled or multiplexed group action; "190 T"; that the
2O fundamental irrep is the single-qubit Clifford group. The nearest number to
190 is 188, the ancilla-assisted 2O trace gate.

## 3. Per link per Trotter step — stated

Multiplicities (`tab:primcost`, identical in BO :797-810, P36 :647-658, P72 :934-947;
BT :548-561 gives the H_I row):

| Hamiltonian | U_F | U_Tr | U_inv | U_mult |
|---|---|---|---|---|
| H_KS | 2 | (d-1)/2 | 3(d-1) | 6(d-1) |
| H_I (improved) | 4 | 3(d-1)/2 | 2 + 11(d-1) | 4 + 26(d-1) |

Closed forms, N_T = C_T x d L^d N_t (T per link per step):

| group | C_T^KS | C_T^I | source |
|---|---|---|---|
| 2T | -- | 4312d - 3640 + (4581.03 + 18.975d) log | BT :564-566 |
| 2O | 2863(d-1) + (22737.8 + 2.3d) log | 11949d - 10157 + (45473.3 + 6.9d) log | BO :813, :821 |
| Sigma(36x3) | 2394(d-1) + (371791 + 4.025d) log | 9884d - 8414 + (744167 + 12.075d) log | P36 :666-681 |
| Sigma(36x3), with FFT | -- | 9632d - 6034 + (1045.93 + 12.075d) log | FFT `tab:simcost` :895-908 |
| Sigma(72x3) | 9338d - 9114 + (372881 + 4.025d) log | 38248d - 32046 + (745758 + 12.075d) log | P72 :949-964 |

Synthesis-error budgeting the papers use with these: eps_T = (coefficient) x d L^d N_t x eps
(P36 :669, :681; P72 :952, :964) -- i.e. a coherent, linear split of the total
synthesis error over all rotations.

## 4. Conventions — stated

| | 2T | 2O | Sigma(36x3) | Sigma(72x3) |
|---|---|---|---|---|
| T per Toffoli | 7 | 7 | 7 | 7 |
| C^nNOT | 4(n-1) Toffoli, n-1 clean anc | 2 ceil(log2 n) - 1 Toffoli, n-2 dirty | 2 ceil(log2 n) - 1 Toffoli, n-2 dirty | 4(n-1) Toffoli, n-1 clean |
| T per R_Z | "at worst 1.15 log(1/eps)" RUS | "on average 1.15 log2(1/eps) (worst -9 + 4 log2)" | same as 2O | "on average 1.15 log(1/eps)" |
| R_X, R_Y | -- | <= 3 R_Z | <= 3 R_Z | <= 1 R_Z + 2 Cliffords |
| fiducial eps | 1e-8 (total) | 1e-8 (total) | 1e-8 (P36 :664) | 1e-10 in the figure, eps_T = 1e-8 in the benchmark |

The rotation figure is the slope-only 1.15 log2(1/eps). The whitepaper's report-wide
convention (T4) is the full Bocharov-Roetteler-Svore fit 1.15 log2(1/eps) + 9.2 with
randomized (incoherent) error splitting. Since E20 (section 7) the models price these
tables in the report's convention; the papers' own price is kept as a legacy record.

## 5. What dominates, and benchmark totals — stated

| group | dominant | benchmark (L^3 = 10^3, N_t = 50, eps_T = 1e-8) |
|---|---|---|
| 2T | U_F, 44% (BT :568) | 2.0e10 (H_I; BT :568), superseded by 1.1e11 (BO :827); with FFT 3.4e9 (FFT `tab:simcost`) |
| 2O | U_F, 99% (H_KS) / 98% (H_I) (BO :827) | 4.1e11 (H_I), 2.0e11 (H_KS); with FFT 5.6e9 |
| Sigma(36x3) | U_F, "over 99% ... regardless of Hamiltonian" (P36 :682) | 7.0e12 (H_I), 3.5e12 (H_KS); with FFT 1.2e10, "a factor of 580" (FFT :870), Fourier share 51% |
| Sigma(72x3) | U_F, "over 99%" (P72 :965) | 7.1e12 (H_I), 3.5e12 (H_KS); "700 times more expensive than Sigma(36x3) with the fast Fourier transform" |

"Per plaquette" occurs in none of the five papers.

## 6. Derived, not in the papers (for checking only)

- Implied Toffoli counts at 7 T: Sigma(36x3) U_inv 17, U_mult 44, U_Tr const 54;
  Sigma(72x3) U_inv 91, U_mult 160; 2T U_mult 22; 2O U_mult 56.
- Implied rotation counts at 1.15 T/log: U_phi 256; Sigma(36x3) U_F 161650; U_Tr 7;
  Sigma(36x3) FFT 102 (the text's "36 R_z" increment does not match its own 36.8 log;
  internal inconsistency at FFT :834-835); 2O U_F 9887; 2T U_Tr 11.
- A naive "4 U_mult + U_Tr" magnetic plaquette would be 1.6e3 + 8 log (Sigma(36x3))
  and 5.9e3 + 8 log (Sigma(72x3)). The papers do not compose it that way: the
  magnetic term is priced per link via `tab:primcost` (3(d-1) U_inv, 6(d-1) U_mult,
  (d-1)/2 U_Tr for H_KS).
- P36 and FFT disagree on Sigma(36x3)'s C_T^I constant (9884d - 8414 vs 9632d - 8192)
  and on the with-FFT total (5.7e9 estimated in P36 vs 1.2e10 in FFT).

## 7. Conversion to the report convention (E20, 2026-10-01) — derived

Ruling E20 (r18; apply_log/r18_core.md): every rotation in every circuit is priced at the
full Bocharov-Roetteler-Svore fit 1.15 log2(1/eps) + 9.2, with eps per circuit from R-TOL.
This supersedes ruling R2's use of the papers' slope-only convention for the gauge tables.
The tables above are converted as follows, and nothing else changes:

- The paper's log coefficient b is read as a rotation count, n_rot = b / 1.15, from the
  printed coefficient (section 6 lists the implied counts).
- The paper's constant a (Toffolis at 7 T, plus any exact T) stays as printed.
- Each primitive costs a + n_rot (1.15 log2(1/eps) + 9.2) = (a + b log2(1/eps)) + 9.2 n_rot.
- The closed forms C_T become C_T(printed) + 9.2 x (log coefficient / 1.15).

Rotation counts per primitive (n_rot):

| primitive | Z3 (derived) | 2T | 2O | Sigma(36x3) | Sigma(72x3) |
|---|---|---|---|---|---|
| U_inv, U_mult | 0 | 0 | 0 | 0 | 0 |
| U_Tr | 1 | 11 | 4 | 7 | 7 |
| U_phi | 1 | -- | -- | 256 | 256 |
| U_F (naive) | 14 | 1000 | 9887.04 | 161650.43 | 161996 |
| U_FFT | 2 (+4 exact T) | 42 | 42 | 102 | 448 (lower bound) |

Two printed coefficients are not multiples of 1.15 (2O U_F 11370.1, Sigma(36x3) U_F 185898).
They are carried exactly as b / 1.15. The FFT text's "36 R_z" increment (FFT :834-835) is
not used; the printed 117.3 log gives 102. The P36/FFT disagreement on the C_T^I constant is
untouched (constants stay as printed).

Code: `groups.PrimitiveCost.t` (report convention, via `groups.rotation_t` = 
`common.t_per_rotation(eps, "rus")`), `groups.c_t_full`, `groups.c_t_rotations`. Legacy:
`PrimitiveCost.t_papers`, `groups.papers_convention_t`, `groups.c_t_stated`, and
`papers=True` on `magnetic_per_link` / `electric_per_link`. Tests: tests/test_common.py
(`test_e20_*` pin the conversion; `test_legacy_*` reproduce each paper's printed formulas).
