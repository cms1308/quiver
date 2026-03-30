# Module 2: Quiver Diagram Generation and Constraints

## Goal
Systematically enumerate all quiver gauge theories satisfying:
1. UV asymptotic freedom (from Module 1)
2. Large N admissibility
3. Gauge anomaly cancellation (SU(N) nodes)
4. Sp(N) Witten anomaly cancellation
5. Consistency of bifundamental representations across mixed-gauge-group edges

---

## 1. Quiver Data Structure

A quiver $Q = (V, E, \ell_V, \ell_E, \ell_{\mathrm{node}})$ consists of:

- **$V$** — finite set of gauge nodes
- **$E \subseteq V \times V$** — directed edges (bifundamental matter)
- **$\ell_V : V \to \{\mathrm{SU}, \mathrm{SO}, \mathrm{Sp}\}$** — gauge group type at each node
- **$\ell_E : E \to \mathcal{R}_{\mathrm{bif}}$** — bifundamental representation on each edge
- **$\ell_{\mathrm{node}} : V \to \mathrm{MultiSet}(\mathcal{R}_{\mathrm{node}})$** — non-bifundamental matter at each node (adjoints, symmetric, antisymmetric)

For **undirected** edges between nodes of the same gauge group type, orientation encodes which copy is fundamental vs. anti-fundamental (for SU(N)); for SO and Sp the orientation is unphysical.

---

## 2. Allowed Gauge Groups

| Type | Group | Fundamental dim | Large N scaling of $T(\mathrm{adj})$ |
|------|-------|-----------------|--------------------------------------|
| SU | SU(N) | $N$ | $N$ |
| SO | SO(N) | $N$ | $N$ |
| Sp | USp(2N) | $2N$ | $N$ |

All nodes share the **same parameter N**. The theory admits a large N limit (planar limit) only when:
- The 't Hooft coupling $\lambda_a = g_a^2 N$ is kept fixed as $N \to \infty$
- $b_0^{(a)} > 0$ for all nodes even as $N \to \infty$ (see Module 1)

---

## 3. Allowed Matter Representations

### Node matter (at a single gauge node)

| Gauge group | Representation | Notation | Large N scaling of $T$ |
|-------------|---------------|---------|------------------------|
| SU(N) | Fundamental | $\square$ | $1/2$ |
| SU(N) | Anti-fundamental | $\bar{\square}$ | $1/2$ |
| SU(N) | Adjoint | $\mathrm{adj}$ | $N$ |
| SU(N) | Rank-2 symmetric | $S$ | $N/2$ |
| SU(N) | Rank-2 sym conjugate | $\bar{S}$ | $N/2$ |
| SU(N) | Rank-2 antisymmetric | $A$ | $N/2$ |
| SU(N) | Rank-2 antisym conjugate | $\bar{A}$ | $N/2$ |
| SO(N) | Vector | $V$ | $1$ |
| SO(N) | Adjoint | $\mathrm{adj}$ | $N$ |
| SO(N) | Rank-2 symmetric (traceless) | $S$ | $N$ |
| Sp(N) | Fundamental | $f$ | $1/2$ |
| Sp(N) | Adjoint | $\mathrm{adj}$ | $N$ |
| Sp(N) | Rank-2 antisymmetric (traceless) | $A$ | $N$ |

> Spinors of SO(N) are excluded at large N (exponential growth of $T$).

### Bifundamental matter (edge between two nodes)

| Node $a$ | Node $b$ | Bifundamental representation | $T_{G_a}$ | $T_{G_b}$ |
|----------|----------|------------------------------|-----------|-----------|
| SU(N) | SU(N) | $(\square, \bar{\square})$ | $N/2$ | $N/2$ |
| SU(N) | SU(N) | $(\square, \square)$ | $N/2$ | $N/2$ |
| SU(N) | SO(N) | $(\square, V)$ | $N/2$ | $N$ |
| SU(N) | Sp(N) | $(\square, f)$ | $N$ | $N/2$ |
| SO(N) | SO(N) | $(V, V)$ | $N$ | $N$ |
| SO(N) | Sp(N) | $(V, f)$ | $2N$ | $N/2$ |
| Sp(N) | Sp(N) | $(f, f)$ | $N$ | $N$ |

> Note: SO(N)–Sp(N) bifundamentals contribute $T_{SO} = 2N$ to the SO node, which is very costly for asymptotic freedom.

---

## 4. Anomaly Cancellation

### 4a. SU(N) Gauge Anomaly

For each SU(N) node $a$, the cubic gauge anomaly must vanish:
$$\sum_{\text{chiral multiplets charged under } G_a} A(r_i) \cdot \dim(\text{other groups}) = 0$$

where $A(r)$ is the anomaly coefficient (third-order Dynkin index):

| SU(N) rep | $A(r)$ |
|-----------|--------|
| $\square$ | $+1$ |
| $\bar{\square}$ | $-1$ |
| adjoint | $0$ |
| $S$ | $+(N+4)$ |
| $\bar{S}$ | $-(N+4)$ |
| $A$ | $+(N-4)$ |
| $\bar{A}$ | $-(N-4)$ |

> Note: $A(A)$ changes sign for $N < 4$. For large N both $S$ and $A$ have $|A| \sim N$. The conjugates $\bar{S}$ and $\bar{A}$ carry opposite anomaly coefficients.

**For a purely bifundamental SU(N) quiver:**

A bifundamental $(\square_a, \bar{\square}_b)$:
- Contributes $A = +N$ to node $a$ (there are $N$ fundamentals of $G_a$)
- Contributes $A = -N$ to node $b$ (there are $N$ anti-fundamentals of $G_b$)

A bifundamental $(\square_a, \square_b)$ (both fundamentals):
- Contributes $A = +N$ to node $a$
- Contributes $A = +N$ to node $b$

**Anomaly-free condition for SU(N)$^k$ quiver:**

Assigning each edge a direction $a \to b$ meaning $(\square_a, \bar{\square}_b)$:
$$\sum_{\text{edges } e: a \to b} N - \sum_{\text{edges } e: b \to a} N + \sum_{\text{node reps at } a} A(r) \cdot d_{\text{other}} = 0$$

For pure bifundamental quivers with all $(\square, \bar{\square})$ edges: the quiver must be **balanced** at each node:
$$\deg^+(a) = \deg^-(a) \quad \forall\, a \in V$$

This means the underlying directed graph must have equal in-degree and out-degree at every SU(N) node.

**With rank-2 tensors $S, \bar{S}, A, \bar{A}$ at node $a$:**
$$N\big(\deg^+(a) - \deg^-(a)\big) + (n_S - n_{\bar{S}})(N+4) + (n_A - n_{\bar{A}})(N-4) = 0$$

For large N, the leading condition is:
$$\deg^+(a) - \deg^-(a) + (n_S - n_{\bar{S}}) + (n_A - n_{\bar{A}}) = 0$$

### 4b. SO(N) and Sp(N) Anomalies

**SO(N):** No cubic gauge anomaly (gauge group is real/pseudo-real for odd N; for SO(N) the third-order Dynkin index vanishes). No constraint from anomaly cancellation.

**Sp(N) = USp(2N):** No cubic gauge anomaly (group is pseudo-real). However, there is a **Witten global anomaly** for odd numbers of fundamental half-hypermultiplets. The condition:
$$n_f^{Sp(N)} \equiv 0 \pmod{2}$$

i.e., the total number of fundamental chiral multiplets at each Sp(N) node must be **even**. (A bifundamental edge from node $b$ to an Sp(N) node contributes a fundamental of Sp(N), counted in this total.)

---

## 5. Large N Constraints Summary

For the quiver to admit a large N limit, at each node $a$ the leading-N contribution to $b_0$ must be positive. Define:

$$B_0^{(a)} \equiv \lim_{N\to\infty} \frac{b_0^{(a)}}{N}$$

Requirements:
- **SU(N) node:** $B_0 = 3 - n_{\mathrm{adj}} - \frac{(n_S + n_{\bar{S}}) + (n_A + n_{\bar{A}})}{2} - \sum_{b \sim a} c_{ab} > 0$
  - where $c_{ab} = 1/2$ for each SU-SU, SU-SO edge at the SU end; $c_{ab} = 1$ for SU-Sp edge at the SU end
- **SO(N) node:** $B_0 = 3 - n_{\mathrm{adj}} - n_S - \sum_{b\sim a}c_{ab} > 0$
  - $c_{ab} = 1$ for SO-SU and SO-SO edges; $c_{ab} = 2$ for SO-Sp edge
- **Sp(N) node:** $B_0 = 3 - n_{\mathrm{adj}} - n_A - \sum_{b\sim a}c_{ab} > 0$
  - $c_{ab} = 1/2$ for Sp-SU and Sp-SO edges; $c_{ab} = 1$ for Sp-Sp edge

---

## 6. Enumeration Strategy

### Step 1: Bound on quiver size
From large N asymptotics, $B_0 > 0$ gives an upper bound on the number of edges and node matter at each node. For example, an SU(N) node with no node matter can have at most **5 bifundamental neighbors** connected via SU-SU edges (since $B_0 = 3 - 5/2 = 1/2 > 0$ but $3 - 6/2 = 0$).

**Upper bounds on neighbor count per node type (pure bifundamental case):**

| Node type | Edge type to neighbor | Max neighbors for $B_0 > 0$ |
|-----------|----------------------|------------------------------|
| SU(N) | all SU-SU | 5 |
| SU(N) | all SU-Sp | 2 (since $c=1$ each) |
| SO(N) | all SO-SU | 2 |
| Sp(N) | all Sp-SU | 5 |
| Sp(N) | all Sp-Sp | 2 |

### Step 2: Enumerate graphs
Enumerate all connected directed multigraphs $G = (V, E)$ up to isomorphism with $|V| \leq V_{\max}$ and $\deg(v) \leq d_{\max}(v)$ for each node.

### Step 3: Assign gauge group types
For each graph, assign gauge group types (SU, SO, Sp) to each node.

### Step 4: Assign edge representations
Determine the allowed bifundamental representations for each edge given the gauge group types at its endpoints.

### Step 5: Assign node matter
Enumerate additional matter at each node: adjoints, symmetric tensors, antisymmetric tensors.

### Step 6: Apply all constraints
For each candidate quiver:
- Check large N $b_0 > 0$ at each node (Module 1)
- Check exact $b_0 > 0$ at each node (Module 1)
- Check SU(N) anomaly cancellation
- Check Sp(N) Witten anomaly condition
- (Optionally) check theory is connected and irreducible

---

## 7. Examples of Valid Quivers

### Example A: SU(N) SQCD with $N_f$ fundamentals
- 1 node (SU(N)), $N_f$ flavor fundamentals + $N_f$ anti-fundamentals
- AF: $b_0 = 3N - N_f > 0 \Rightarrow N_f < 3N$
- Anomaly: $N_f \cdot 1 + N_f \cdot(-1) = 0$ ✓
- Large N: $N_f = \kappa N$ with $\kappa < 3$

### Example B: $\mathcal{N}=2$ SQCD — SU(N) with adjoint + $N_f$ fund
- 1 node, 1 adjoint + $N_f$ fundamentals + $N_f$ anti-fundamentals
- AF: $b_0 = 3N - N - N_f = 2N - N_f > 0 \Rightarrow N_f < 2N$

### Example C: Linear SU(N)$^k$ quiver (necklace tail)
- $k$ nodes, bifundamentals between adjacent nodes
- All nodes: $b_0 = 2N$ (interior) or $5N/2$ (ends)
- Anomaly: each interior node balanced ✓

### Example D: Circular (affine $A_{k-1}$) quiver SU(N)$^k$
- $k$ nodes in a circle, one bifundamental per adjacent pair
- All nodes: $b_0 = 2N > 0$ ✓, balanced ✓

### Example E: Sp(N)–SU(2N)–SO(N) tail
- Three nodes: Sp(N)–SU(2N)–SO(N)
- Need to check all constraints carefully

---

## 8. Known Classification for Special Cases

### Single-node theories (quivers with one gauge group)
- **SU(N):** $(n_S + n_{\bar{S}}) + (n_A + n_{\bar{A}}) < 6$ (rank-2 tensors), $n_f < 6N - \frac{N+2}{2}(n_S+n_{\bar{S}}) - \frac{N-2}{2}(n_A+n_{\bar{A}})$ (fundamentals)
- **SO(N):** only vectors and symmetric tensors; $n_S = 0$ for AF in large N (since $T(S) = N+2 \approx N$)
- **Sp(N):** similar analysis; adjoints very costly

### Two-node theories
The full analysis must also pass Module 4 (non-trivial IR fixed point check).

---

## 9. Implementation Notes

```python
# Pseudocode for quiver generation

def enumerate_quivers(max_nodes, max_edges_per_node):
    graphs = enumerate_directed_graphs(max_nodes, max_edges_per_node)
    valid_quivers = []

    for G in graphs:
        for gauge_assignment in assign_gauge_groups(G):
            for edge_reps in assign_edge_representations(G, gauge_assignment):
                for node_matter in assign_node_matter(G, gauge_assignment):
                    q = Quiver(G, gauge_assignment, edge_reps, node_matter)
                    if check_large_N_AF(q) and check_anomalies(q):
                        valid_quivers.append(q)

    return valid_quivers


def check_large_N_AF(quiver):
    for node in quiver.nodes:
        B0 = compute_large_N_beta(quiver, node)
        if B0 <= 0:
            return False
    return True


def check_anomalies(quiver):
    for node in quiver.nodes:
        if node.gauge_group == "SU":
            if su_anomaly(quiver, node) != 0:
                return False
        elif node.gauge_group == "Sp":
            if sp_witten_anomaly(quiver, node) != 0:
                return False
    return True
```
