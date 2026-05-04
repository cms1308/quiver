# Module 6 — Marginal Operator Search

Enumerate gauge-invariant chiral operators in two-node quiver gauge theories,
evaluate their R-charges via finite-N a-maximization, and identify those that
are marginal (R = 2) at every tested N and flavor singlets.

## Physical motivation

A flavor-singlet chiral operator with R = 2 at the IR fixed point is a
candidate for an **exactly marginal deformation** in the Leigh–Strassler
sense: turning it on preserves all global symmetries, so the only running
of its coupling comes from gauge-coupling beta functions. The dimension of
the conformal manifold at any N=1 SCFT is

$$
\dim \mathcal{M}_{c} = \#\{\text{marginal flavor-singlet ops}\} - \#\{\text{independent gauge anomaly conditions}\}
$$

(Green–Komargodski–Seiberg–Tachikawa–Wecht).

Restricting the search to operators that are marginal **at every** $N$ in
$\{10, 20, 30\}$ filters out coincidental marginality at any single $N$:
any operator passing this test has an N-independent R-charge equal to 2,
indicating that R-charges of its constituents are pinned by exact rational
relations (e.g. R = 1 for an adjoint coupled to balanced bifundamentals,
R = 2/3 for class-44 matter).

## Pipeline

```
  quivers.db row ──▶ quiver_from_row ──▶ Quiver
                                          │
                                          ▼
              a_maximize at N ∈ {10,20,30} ──▶ R_per_N: {N → {label → R}}
                                          │
                                          ▼
              enumerate_candidates(quiver) ──▶ list[CandidateOp]
                                          │
                                          ▼
                  for each candidate:
                    R(O; N) = Σ mult · R[label] (per N)
                    keep iff |R-2| < 1e-6 at every N AND flavor singlet
                                          │
                                          ▼
                          marginal_operators.tex
```

## Operator enumeration

A `CandidateOp` is a multiset of field labels (matching the labels used by
`a_maximization.build_fields`: `node{i}_{rep}` for intra-node fields,
`edge_{src}_{dst}_{rep}` for bifundamentals). Three sub-enumerators:

1. **single-node**: products of fields at one gauge node. Index balance:
   - SU node: number of upper indices = number of lower indices (mesonic;
     baryons subleading at large N and exhibit N-dependent dimension, so
     excluded by construction).
   - SO/Sp node: total index count is even (δ or Ω contracts any pair).

2. **bifund-loop**: closed walks in the bifundamental multigraph. Each edge
   carries an orientation given by its rep (`+-` directed; `++`/`--`
   symmetric same-chirality pairs; `std` for SO/Sp-SO/Sp). Cyclic-word
   canonical form (lex-min over rotations + reflection) deduplicates
   physically equivalent operators.

3. **mixed**: rank-2 tensors at endpoints of bifundamental loops, e.g.
   $\mathrm{tr}(\bar{Q}_{12}\,Q_{12}\,\mathrm{adj}_1^2)$.

Operators of degree ≥ 2 are emitted (a single trace of one field is either
identically zero, e.g. $\mathrm{tr}(\Phi^a T^a) = 0$, or non-singlet).

## Flavor singlet check (deferred)

A naïve "multiplets of size $\ge 2$ carry flavor charge, multiplicity-1
fields are flavor-trivial" rule is implemented in `is_flavor_singlet` but
**is not used by default** because it is incomplete:

- Every chiral field with multiplicity 1 still carries a $U(1)$ flavor (a
  phase rotation), and only specific operator combinations are
  $U(1)$-neutral.
- For multiplets of size $n \ge 2$, the operator decomposes under $U(n)$
  into reps; only the singlet rep has zero flavor charge. Different cyclic
  orderings of the same multiset can sit in different $U(n)$ reps.

Identifying the genuine flavor-singlet subset requires explicit
$U(1)^k \times \prod_m U(n_m)$ rep-decomposition (Young tableau
bookkeeping). This is left as future work.

The `find_marginal_ops` helper and `scripts/dump_marginal_operators.py`
report **all** always-marginal operators; users analyse flavor charges by
hand. The `is_flavor_singlet` heuristic is retained in the module but
unused by the default pipeline.

## Trivial gauge-anomaly marginals

For each gauge node $a$, the anomaly-free condition

$$
T(\mathrm{adj}_a) + \sum_i T_a(\rho_i)\,(R_i - 1) = 0
$$

can be re-stated as: an operator with multiplicities $n_i \propto T_a(\rho_i)$,
restricted to node $a$, has $R = 2$ identically. These "Konishi-like"
operators are always among the marginal singlets and contribute trivially
to the conformal manifold count — they are *not* genuine exactly marginal
deformations beyond the gauge anomaly. Subtracting them from the per-theory
count gives the dimension of the conformal manifold.

## Usage

```bash
# Single theory
python3 scripts/dump_marginal_operators.py --theory-id 215 --max-degree 4

# All theories in one universality class
python3 scripts/dump_marginal_operators.py --class-id 9 --max-degree 6

# Full DB sweep + TeX output → paper/sections/generated/marginal_operators.tex
python3 scripts/dump_marginal_operators.py --all --max-degree 6 --tex

# Acceptance tests
python3 tests/test_marginal_operators.py
```

## Results summary (max-degree = 6, N ∈ {10, 20, 30})

Of 21,315 theories in `quivers.db`:
- 15,743 processed (5,572 skipped where finite-N a-max diverges or class is null)
- 350 theories have at least one always-marginal operator

Class 44 (R_bif = 2/3 leading-N marginal class) shows $\mathrm{tr}(\mathrm{adj}^3)$
marginal (R_adj = 2/3 at every N). Class 9 features bifundamentals with
R_bif = 1/2 yielding $\mathrm{tr}((Q\bar{Q})^2)$ and $\mathrm{tr}(\mathrm{adj}^2)$
when adj is at R = 1.

To inspect per-theory results inline, run:

```
python3 two_node_db.py show <class_id>
```

— a `Marginal ops` column is appended to the theory table.

## Scope cuts

- Baryons and determinants ($\epsilon^{a_1\dots a_N}$ contractions) excluded
  — their dimension scales with $N$, so they are never marginal at all $N$.
- F-term equivalences from a chosen superpotential are not modeled; the
  count is "candidates for exactly marginal" not the final manifold dimension.
- Multi-trace operators are not enumerated (subleading at large $N$).

## Files

- `marginal_operators.py` — module: data model, enumerators, R evaluation,
  filters, DB row parsing
- `scripts/dump_marginal_operators.py` — CLI sweep, TeX output
- `tests/test_marginal_operators.py` — acceptance tests (13 tests)
- `paper/sections/generated/marginal_operators.tex` — auto-generated table
