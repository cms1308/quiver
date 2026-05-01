# Module 5: Final Classification

## Goal
Combine all preceding modules to produce the complete list of 4D N=1 quiver gauge theories that are:
1. UV asymptotically free (Module 1, 2)
2. Admit a large N limit (Module 2)
3. Satisfy all consistency conditions (anomaly cancellation, Witten anomaly) (Module 2)
4. Flow to a non-trivial IR superconformal fixed point (Modules 3, 4)

---

## 1. Classification Pipeline

```
┌─────────────────────────────────────────┐
│  Enumerate quivers (Module 2)           │
│  - Graph topologies                     │
│  - Gauge group assignments              │
│  - Matter content                       │
└────────────────┬────────────────────────┘
                 │
                 ▼
┌─────────────────────────────────────────┐
│  Filter: UV asymptotic freedom (Mod 1)  │
│  b_0^(a) > 0 for all nodes a           │
│  B_0^(a) > 0 at large N               │
└────────────────┬────────────────────────┘
                 │
                 ▼
┌─────────────────────────────────────────┐
│  Filter: Consistency conditions (Mod 2) │
│  - SU(N) gauge anomaly = 0             │
│  - Sp(N) Witten condition              │
└────────────────┬────────────────────────┘
                 │
                 ▼
┌─────────────────────────────────────────┐
│  a-maximization (Module 3)              │
│  - Compute R_i*, a*, c*                │
│  - Check unitarity bounds              │
└────────────────┬────────────────────────┘
                 │
                 ▼
┌─────────────────────────────────────────┐
│  IR fixed point check (Module 4)        │
│  - Boundary analysis B_a < 0 for all a │
└────────────────┬────────────────────────┘
                 │
                 ▼
┌─────────────────────────────────────────┐
│  CLASSIFIED THEORIES                    │
└─────────────────────────────────────────┘
```

---

## 2. Classification by Quiver Topology

We organize results by the structure of the quiver graph.

### 2a. Single-Node Theories

The complete classification of single-node theories (SU(N), SO(N), Sp(N) with matter in fundamental, adjoint, rank-2 symmetric, and rank-2 antisymmetric representations) that are asymptotically free and have non-trivial IR fixed points is carried out in:

- **arXiv:2007.16165** — classification of 4D $\mathcal{N}=1$ gauge theories with rank-2 matter
- **arXiv:2510.19136** — further classification results

These references should be consulted for the exact conformal windows, anomaly-free conditions, and IR R-charges of all single-node theories.

### 2b. Two-Node Quivers (Classified)

The full two-node classification is stored in `quivers.db` (built by `two_node_db.py`). It contains **4560 theories in 135 universality classes** across all six gauge pair types (SU-SU, SU-SO, SU-Sp, SO-SO, SO-Sp, Sp-Sp).

Each class has:
- Exact symbolic R-charges and $a/N^2$ (76 classes) or numerical approximation (59 classes)
- Chiral excess $\delta(N) = aN + b$ per node
- AF bound $N_f < \alpha N + \gamma$ per node

**Key structural features:**
- SU-SU quivers dominate: 79 classes, 3366 theories
- Many classes contain dozens of theories that differ only in their rank-2 matter composition but share the same large-N physics
- Edge configurations include multiple bifundamentals, mixed chirality ($++$, $--$, $+-$), and multi-edges up to multiplicity 2

### 2c. Linear Quivers (Chain)

$G_1 - G_2 - \cdots - G_k$ with bifundamentals on each link.

**Pure SU(N)$^k$ linear quiver:**
- End nodes: $b_0 = 5N/2$ ✓
- Interior nodes: $b_0 = 2N$ ✓
- Anomaly: each interior node has equal in/out arrows ✓
- Valid for any $k$ at large N

**With additional node matter:**
- Add adjoints, symmetric, antisymmetric at specific nodes to tune conformal window

### 2d. Circular (Necklace) Quivers

$G_1 - G_2 - \cdots - G_k - G_1$ (cyclic).

**Pure SU(N)$^k$ circular quiver (affine $A_{k-1}$ quiver):**
- All nodes: $b_0 = 2N > 0$ ✓
- Balanced at each node ✓
- Valid for any $k$ ≥ 2

At $k=1$ with 1 adjoint: $b_0 = 3N - N = 2N$ (pure $\mathcal{N}=2$ SYM, marginal not AF).

**SU(N)$^k$ with one adjoint per node (affine $\hat{A}_{k-1}$, $\mathcal{N}=2$ quiver):**
- Each node: $b_0 = 3N - N - 2 \cdot N/2 = N > 0$... wait, this needs careful counting.

### 2e. Star (Hub-and-Spoke) Quivers

One central node connected to several outer nodes.

**SU(N) hub with $k$ SU(N) spokes:**
- Hub node: $b_0 = 3N - k \cdot N/2$; requires $k \leq 5$
- Spoke nodes: $b_0 = 3N - N/2 = 5N/2 > 0$ ✓

For $k = 5$: hub has $b_0 = N/2 > 0$ ✓ (marginally AF at large N; finite N corrections matter)

---

## 3. Large N Admissible Quivers — Catalog

We organize by the **large N leading beta function coefficient** $B_0^{(a)} = \lim_{N\to\infty} b_0^{(a)}/N$.

### Requirement: $B_0^{(a)} > 0$ for all nodes

For **SU(N) nodes** with $n_b$ bifundamental neighbors (all SU-SU edges):
$$B_0 = 3 - n_{\mathrm{adj}} - \frac{n_S + n_A}{2} - \frac{n_b}{2}$$

| $n_{\mathrm{adj}}$ | $n_S + n_A$ | $n_b^{\max}$ |
|-------------------|-------------|--------------|
| 0 | 0 | 5 |
| 0 | 1 | 3 |
| 0 | 2 | 1 |
| 1 | 0 | 3 |
| 1 | 1 | 1 |
| 2 | 0 | 1 |

---

## 4. Theories Known to Have Non-Trivial IR Fixed Points

### 4a. $\mathcal{N}=1$ SQCD (SU, SO, Sp variants)
- Seiberg's conformal window: $\frac{3}{2}N < N_f < 3N$ for SU(N)

### 4b. Kutasov–Schwimmer Models
- SU(N) with one $A$ (or $\bar{A}$) + $N_f$ fundamentals (+ anti-fundamentals for anomaly cancellation)
- Conformal window: $N+2 \leq N_f \leq \frac{3}{2}(N-2)$ (depends on exact conventions)

### 4c. Circular SU(N)$^k$ Quivers
- All nodes balanced, always UV-AF for any $k$
- IR analysis: $k=2$ (two nodes) is equivalent to SU(N) SQCD-like analysis

### 4d. Linear SU(N)$^k$ Quivers with Flavors
- Adding $N_f$ fundamental flavors at end nodes extends the conformal window

### 4e. Mixed Gauge Group Quivers
- SO(N)–Sp(N) pairs: special care needed due to large Dynkin indices for bifundamentals
- SU(N)–SO(N) and SU(N)–Sp(N) quivers: more natural candidates

---

## 5. Output Format for Each Classified Theory

For each valid theory, record:

```
Theory ID: [unique identifier]
Quiver type: [linear / circular / star / other]
Nodes: [(G_1, type_1), (G_2, type_2), ...]
Edges: [(a, b, rep_ab), ...]
Node matter: {a: [rep_1, rep_2, ...], ...}
Superpotential: [W = ...]

UV properties:
  b_0^(a) = [...]  for each a
  B_0^(a) = [...]  (large N coefficient)
  Anomaly cancellation: [✓/✗]

IR properties:
  R_charges: {field_i: R_i*, ...}
  a* = [...]
  c* = [...]
  a/c = [...]
  B_a = [...]  (boundary analysis values, all < 0 ✓)

Notes: [special properties, known dual, etc.]
```

---

## 6. Results by Quiver Size

### 1-node theories (complete)
19 quivers parametrized by $N_f$. Exact large-N R-charges and $a/N^2$ computed symbolically for all. These reproduce known SQCD and Kutasov-type conformal windows.

### 2-node theories (complete, pending Module 4 IR filter)
4560 theories in 135 universality classes. Exact symbolic results for 76 classes; numerical approximation for 59 classes (sympy solver timeout at 30s due to 4+ free parameters in the a-function). The Module 4 boundary analysis ($\mathcal{B}_a < 0$) has not yet been applied to filter for non-trivial IR fixed points.

### $k \geq 3$-node theories (not yet enumerated)
The enumeration infrastructure (`enumerate_quivers(n_nodes)`) supports arbitrary $n_\text{nodes}$, but the combinatorial explosion grows rapidly. The 3-node enumeration has been tested and produces a large number of candidates.

---

## 7. Special Cases and Marginal Theories

Theories with $B_0^{(a)} = 0$ for some node $a$ are **marginally AF at large N**: they are asymptotically free for any finite $N$ but the AF condition becomes marginal as $N \to \infty$.

These include:
- SU(N) with 6 bifundamental SU-SU neighbors: $B_0 = 3 - 3 = 0$
- SO(N) node with 3 SO-SU bifundamentals: $B_0 = 3 - 3 \cdot 1 = 0$ (needs careful $\mathcal{O}(1)$ analysis)

For these, the finite-N corrections to $b_0$ determine asymptotic freedom. They may or may not be included depending on whether one requires strict $b_0 > 0$ or $b_0 \geq 0$.

---

## 8. Generalization Directions

The current classification assumes all gauge groups have the same N. Natural generalizations:

1. **Different ranks**: $G_a = SU(N_a)$ with $N_a$ varying — allows for "quiver tails" with decreasing ranks
2. **Framing**: adding decoupled flavor nodes (non-gauged SU($N_f$)) at various nodes
3. **Superpotential deformations**: adding relevant/marginal operators that modify R-charges and conformal window
4. **Global symmetry structure**: classifying the flavor symmetries and anomalies of the classified theories

---

## 9. Implementation: `two_node_db.py` and `a_maximization_large_N.py`

### Single-Node Classification

`a_maximization_large_N.py` contains `_classify_single_node()` which prints a table of all 19 single-node theories with their matter content, chiral excess $\delta(N)$, AF bound, leading-order R-charges, and $a/N^2$.

### Two-Node Classification: `two_node_db.py`

SQLite database + CLI for browsing two-node quiver universality classes.

**Database schema:**
- `universality_class` — `class_id`, `gauge_pair`, `a_over_N2` (numerical), `n_theories`, `a_exact`, `c_exact`, `R_exact` (symbolic strings from sympy)
- `theory` — per-theory metadata including matter content, edges, delta per node ($\delta_0$, $\delta_1$ as $aN+b$ strings), AF bounds
- `build_info` — provenance (timestamp, parameters)

**Build pipeline (3 phases):**
1. **Enumerate + fast scan** (~3 min): generate all 6099 two-node quivers via `enumerate_quivers(2)`, compute numerical $a/N^2$ for each with `a_maximize_large_N_fast`
2. **Cluster into universality classes**: group by gauge pair, cluster by numerical $a/N^2$ within tolerance $10^{-5}$. Produces 135 classes across 4560 valid theories (theories with $a/N^2 > 2$ are filtered out).
3. **Exact symbolic solve**: for each of the 135 classes, pick the representative with fewest fields, run `a_maximize_large_N` with 30s timeout. Validate exact result against numerical (reject if $|\text{exact} - \text{numerical}| > 0.001$). Result: 76 exact classes, 59 numerical-only.

**CLI commands:**
```bash
python two_node_db.py build                    # build database
python two_node_db.py classes                  # list all 135 classes
python two_node_db.py classes --pair SU-SU     # filter by gauge pair
python two_node_db.py show 7                   # show all theories in class 7
python two_node_db.py search --matter0 adj     # search by node matter
python two_node_db.py search --delta0 N        # search by chiral excess
python two_node_db.py stats                    # database summary
```

**Universality classes** arise because at large N, all rank-2 representations ($S$, $\bar{S}$, $A$, $\bar{A}$) of SU(N) have identical leading-order Dynkin index ($T_\text{lead} = 1/2$) and dimension ($\dim_\text{lead} = 1/2$). Theories with different rank-2 content but the same total counts share the same $a/N^2$ and R-charges at leading order.

### Classification Numbers

| Gauge Pair | Classes | Theories |
|------------|---------|----------|
| SU–SU | 79 | 3366 |
| SU–SO | 27 | 588 |
| SU–Sp | 18 | 354 |
| SO–SO | 4 | 56 |
| SO–Sp | 4 | 112 |
| Sp–Sp | 3 | 84 |
| **Total** | **135** | **4560** |

---

## 10. Cross-Checks and Validation

- **$\mathcal{N}=2$ theories**: Setting the adjoint chiral R-charge to $R = 2/3$ should recover $\mathcal{N}=2$ results
- **Known dualities**: Seiberg duals of classified theories should also appear in the classification (or relate to theories outside the classification via relevant deformations)
- **Unitarity**: All gauge-invariant chiral operators $\mathcal{O}$ must satisfy $R[\mathcal{O}] \geq 2/3$; elementary (gauge-charged) fields are not directly constrained and may have $R < 2/3$
- **$a$-theorem consistency**: $a_{\mathrm{UV}} \geq a_{\mathrm{IR}}$ for all theories in the classification (Komargodski-Schwimmer)
