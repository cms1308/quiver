# Reference and Anchor Map

**Analysis Date:** 2026-04-14

## Active Anchor Registry

| Anchor ID | Anchor | Type | Source / Locator | Why It Matters | Contract Subject IDs | Required Action | Carry Forward To |
|-----------|--------|------|------------------|----------------|----------------------|-----------------|------------------|
| REF-001 | Intriligator-Wecht a-maximization | method | hep-th/0304128 (Intriligator & Wecht, 2003) | Defines the primary computational method: trial a-function maximization subject to anomaly-free constraints determines exact IR R-charges. The entire Module 3 pipeline is built on this. | -- | use | execution, verification |
| REF-002 | IR fixed point boundary analysis | method | hep-th/0502049 | Provides the criterion B_a < 0 for existence of non-trivial IR fixed points in multi-node quivers. Module 4 is entirely based on this method. The project strategy ("follows hep-th/0502049") is stated in README.md line 10. | -- | use | planning, execution |
| REF-003 | NSVZ exact beta function | method | Novikov, Shifman, Vainshtein, Zakharov (1983); see reviews | Provides the exact beta function for N=1 SUSY gauge theories. Used to derive the anomaly-free R-symmetry condition (Tr[R G^2] = 0 at fixed point) and to connect gamma_i = 3R_i - 2. Referenced in `module1_beta_functions.md` (Sec. 6). | -- | cite | execution, verification |
| REF-004 | Single-node classification (rank-2 matter) | benchmark | arXiv:2007.16165 | Complete classification of single-node 4D N=1 gauge theories with rank-2 matter. The project's 19 single-node theories should reproduce these results. Referenced in `module2_quiver_generation.md` (Sec. 7) and `module5_classification.md` (Sec. 2a). | -- | compare | verification |
| REF-005 | Single-node classification (extended) | benchmark | arXiv:2510.19136 | Further single-node classification results. Must be consistent with project's single-node table. Referenced alongside REF-004 in same locations. | -- | compare | verification |
| REF-006 | Seiberg SQCD conformal window | benchmark | Seiberg (1994-1995); Intriligator & Seiberg reviews | SU(N) SQCD conformal window: 3N/2 < N_f < 3N. R^* = 1 - N/N_f. Used as primary sanity check for a-maximization code. Verified in `module3_a_maximization.md` (Sec. 10). | -- | compare | verification |
| REF-007 | Kutasov-Schwimmer models | benchmark | Kutasov & Schwimmer; Intriligator & Wecht | SU(N) with one antisymmetric + flavors: known conformal windows. Referenced in `module5_classification.md` (Sec. 4b). Project's single-node results should reproduce these. | -- | compare | verification |
| REF-008 | Kutasov-Schwimmer decoupling procedure | method | Kutasov & Schwimmer; Intriligator & Wecht | a-maximization with iterative unitarity decoupling: when R[O] < 2/3, pin R[O] = 2/3 and re-maximize. Implemented in `a_maximization.py` (`a_maximize_with_decoupling`) and `a_maximization_large_N.py`. Referenced in `module3_a_maximization.md` (Sec. 9). | -- | use | execution |
| REF-009 | Douglas-Moore quiver gauge theories | background | Douglas & Moore (1996) | Original construction of quiver gauge theories from D-branes at singularities. Provides the mathematical framework of quiver diagrams used throughout. Referenced in `README.md` (Key References). | -- | cite | writing |
| REF-010 | Komargodski-Schwimmer a-theorem | benchmark | Komargodski & Schwimmer (2011) | a-theorem: a_UV >= a_IR for any RG flow between conformal fixed points in 4D. Provides a consistency check on all classified theories. Referenced in `module5_classification.md` (Sec. 10). | -- | compare | verification |
| REF-011 | Banks-Zaks fixed point | benchmark | Banks & Zaks (1982) | Perturbative IR fixed point near the upper boundary of the conformal window (N_f near 3N for SU(N) SQCD). Provides a weak-coupling cross-check for a-maximization results in the Banks-Zaks limit. Referenced in `module3_a_maximization.md` (Sec. 10) and `module4_ir_fixed_points.md` (Sec. 9). | -- | compare | verification |

## Benchmarks and Comparison Targets

- **SU(N) SQCD R-charges**: R^* = 1 - N/N_f for N_f fundamentals + N_f anti-fundamentals. Conformal window 3N/2 < N_f < 3N.
  - Source: REF-006 (Seiberg); verified in `module3_a_maximization.md` (Sec. 10)
  - Compared in: `a_maximization.py` (validation in `__main__` block, expects R = 0.75 for N=5, N_f=10)
  - Status: matched

- **Meson unitarity bound for SQCD**: R[M] = 2R^* = 2 - 2N/N_f >= 2/3, giving N_f >= 3N/2 as lower conformal window boundary.
  - Source: REF-006
  - Compared in: `module3_a_maximization.md` (Sec. 10)
  - Status: matched

- **Circular SU(N)^2 pure bifundamental**: R_Q = R_Qtilde = 0, no interacting SCFT (B_a > 0).
  - Source: Direct computation in `module3_a_maximization.md` (Sec. 11)
  - Compared in: `a_maximization.py` (validation block)
  - Status: matched (confirms no fixed point for pure bifundamental without flavors)

- **Single-node classification (19 theories)**: exact large-N R-charges and a/N^2 for all 19 single-node quivers.
  - Source: REF-004, REF-005 (arXiv:2007.16165, arXiv:2510.19136)
  - Compared in: `a_maximization_large_N.py` (`_classify_single_node()`)
  - Status: pending explicit comparison against literature values

- **Two-node universality classes (135 classes, 4560 theories)**: exact symbolic results for 76 classes, numerical for 59.
  - Source: novel result (no prior literature benchmark)
  - Compared in: `quivers.db`; `two_node_db.py`
  - Status: internal consistency (exact vs numerical agree to 0.001)

- **a-theorem consistency (a_UV >= a_IR)**: must hold for all classified theories.
  - Source: REF-010 (Komargodski-Schwimmer)
  - Compared in: not yet checked
  - Status: pending

- **Seiberg duality closure**: Seiberg duals of classified theories should also appear in the classification.
  - Source: REF-006
  - Compared in: `module5_classification.md` (Sec. 10, mentioned as cross-check)
  - Status: pending

## Prior Artifacts and Baselines

- `quivers.db`: SQLite database containing 4560 theories in 135 universality classes across all six gauge pair types (SU-SU, SU-SO, SU-Sp, SO-SO, SO-Sp, Sp-Sp). Contains exact symbolic R-charges and a/N^2 for 76 classes and numerical approximations for 59 classes. Built by `two_node_db.py`. This is the primary deliverable artifact and must be preserved and extended (not overwritten) as Module 4 filtering is applied.

- `beta_functions.py`: Validated implementation of Dynkin indices, anomaly coefficients, and one-loop beta function computation. All downstream modules depend on this. Changes to Dynkin index values would propagate to every result in the project.

- `quiver_generation.py`: Enumeration engine producing 19 single-node and 6099 two-node quiver candidates. Includes anomaly cancellation, Witten anomaly checks, and charge conjugation deduplication. The 6099 -> 4560 reduction (after filtering a/N^2 > 2 outliers) is a baseline count.

- `a_maximization.py`: Numerical a-maximization at fixed (N, N_f) with iterative unitarity decoupling. Validated against SU(N) SQCD analytic result.

- `a_maximization_large_N.py`: Exact symbolic large-N a-maximization using sympy. Produces algebraic R-charges (rational or involving sqrt). 76/135 classes solved exactly within 30s timeout.

## Open Reference Questions

- **hep-th/0502049 applicability to single-node theories**: Module 4 notes (Sec. 9) identify a subtlety where the boundary analysis B < 0 for single-node SQCD gives N_f > 3N, which contradicts the known conformal window 3N/2 < N_f < 3N. The resolution involves Seiberg duality and the free magnetic phase. The boundary analysis is designed for multi-node quivers; its applicability to single-node theories needs careful treatment.

- **Superpotential terms**: The entire classification assumes zero superpotential (W = 0). Adding superpotential terms (e.g., W = Tr Phi Q Qtilde for adjoint + flavor theories) would add linear constraints R[Phi] + R[Q] + R[Qtilde] = 2 and change the R-charge landscape. This is explicitly noted as an open question in `plan.md` (Sec. "Open questions", item 1).

- **59 numerically-solved classes**: 59 of 135 universality classes have only numerical R-charges (sympy timeout at 30s). These involve 4+ free parameters and quartic/higher a-functions. Whether exact algebraic solutions exist for all of them, or whether some require higher-degree algebraic extensions, is unknown.

- **Module 4 implementation**: The boundary analysis (B_a < 0) is fully specified in `module4_ir_fixed_points.md` but NOT implemented in code. All 4560 theories in the database have NOT been filtered for IR fixed point existence. This is the critical missing step before the classification is complete.

- **k >= 3 node theories**: Enumeration infrastructure supports arbitrary node count, but 3+ node quivers have not been systematically enumerated or classified. The combinatorial explosion is acknowledged in `module5_classification.md` (Sec. 6).

- **Explicit comparison with arXiv:2007.16165 and arXiv:2510.19136**: The project's 19 single-node theories should reproduce the exact conformal windows and R-charges from these papers. This comparison has not been performed explicitly.

## Background Reading

- **Quiver gauge theory reviews** (various): General background on quiver gauge theories from string theory constructions (D-branes at singularities, orbifolds). Referenced in `README.md` as "Douglas & Moore (1996); various reviews." Provides physical motivation but is not directly used in the classification pipeline.

- **NSVZ beta function reviews**: The exact beta function is cited but only the one-loop coefficient and fixed-point condition are used. A deeper understanding of the NSVZ formula (e.g., scheme dependence, higher-order corrections) is background context.

- **Superconformal index / partition function methods**: Not used in this project, but provide alternative approaches to studying the same IR fixed points. Could serve as independent verification of a-maximization results.

- **Large N techniques (planar diagrams, 't Hooft limit)**: The project uses large N only at the level of leading-order coefficients (b_0/N, a/N^2). The full 't Hooft expansion (planar vs non-planar) is background context.

---

_Reference map: 2026-04-14_
