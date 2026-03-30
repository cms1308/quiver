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

| Type | Group |
|------|-------|
| SU | SU(N) |
| SO | SO(N) |
| Sp | USp(2N) |

All nodes share the **same parameter N**. $b_0^{(a)} > 0$ must hold for all valid N at every node (see Module 1).

---

## 3. Allowed Matter Representations

### Node matter (at a single gauge node)

| Gauge group | Representation | Notation | Dynkin index $T$ |
|-------------|---------------|----------|------------------|
| SU(N) | Fundamental | $\square$ | $\frac{1}{2}$ |
| SU(N) | Anti-fundamental | $\bar{\square}$ | $\frac{1}{2}$ |
| SU(N) | Adjoint | $\mathrm{adj}$ | $N$ |
| SU(N) | Rank-2 symmetric | $S$ | $\frac{N+2}{2}$ |
| SU(N) | Rank-2 sym conjugate | $\bar{S}$ | $\frac{N+2}{2}$ |
| SU(N) | Rank-2 antisymmetric | $A$ | $\frac{N-2}{2}$ |
| SU(N) | Rank-2 antisym conjugate | $\bar{A}$ | $\frac{N-2}{2}$ |
| SO(N) | Vector | $V$ | $1$ |
| SO(N) | Adjoint | $\mathrm{adj}$ | $N-2$ |
| SO(N) | Rank-2 symmetric (traceless) | $S$ | $N+2$ |
| Sp(N) | Fundamental | $f$ | $\frac{1}{2}$ |
| Sp(N) | Adjoint | $\mathrm{adj}$ | $N+1$ |
| Sp(N) | Rank-2 antisymmetric (traceless) | $A$ | $N-1$ |

> Spinors of SO(N) are excluded because their Dynkin index grows exponentially with N.

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

**For a purely bifundamental SU(N) quiver:**

A bifundamental $(\square_a, \bar{\square}_b)$:
- Contributes $A = +N$ to node $a$ (there are $N$ fundamentals of $G_a$)
- Contributes $A = -N$ to node $b$ (there are $N$ anti-fundamentals of $G_b$)

A bifundamental $(\square_a, \square_b)$ (both fundamentals):
- Contributes $A = +N$ to node $a$
- Contributes $A = +N$ to node $b$

**Anomaly-free condition for SU(N)$^k$ quiver:**

Denote:
- $n_{ab}^{+-}$: number of $(\square_a, \bar{\square}_b)$ bifundamentals from node $a$ to node $b$ (directed)
- $n_{ab}^{++}$: number of $(\square_a, \square_b)$ bifundamentals between nodes $a$ and $b$ (undirected; $n_{ab}^{++} = n_{ba}^{++}$)

Each $(\square_a, \bar{\square}_b)$ contributes $+N$ to node $a$ and $-N$ to node $b$.
Each $(\square_a, \square_b)$ contributes $+N$ to **both** node $a$ and node $b$.

The general anomaly-free condition at node $a$:
$$N\sum_{b \neq a} \Big(n_{ab}^{+-} - n_{ba}^{+-} + n_{ab}^{++}\Big) + (n_f - n_{\bar{f}}) + (n_S - n_{\bar{S}})(N+4) + (n_A - n_{\bar{A}})(N-4) = 0$$

where $n_f, n_{\bar{f}}$ are the numbers of standalone fundamental and anti-fundamental chirals at node $a$ (node-level matter, not bifundamentals).

Note: $(\square, \square)$ edges contribute $+N$ to both endpoints and must be compensated by net incoming $(\square, \bar{\square})$ arrows, anti-fundamentals $\bar{\square}$, or $\bar{S}/\bar{A}$ at the node.

### 4b. SO(N) and Sp(N) Anomalies

**SO(N):** No cubic gauge anomaly (gauge group is real/pseudo-real for odd N; for SO(N) the third-order Dynkin index vanishes). No constraint from anomaly cancellation.

**Sp(N) = USp(2N):** No cubic gauge anomaly (group is pseudo-real). However, there is a **Witten global anomaly** for odd numbers of fundamental half-hypermultiplets. The condition:
$$n_f^{Sp(N)} \equiv 0 \pmod{2}$$

i.e., the total number of fundamental chiral multiplets at each Sp(N) node must be **even**. (A bifundamental edge from node $b$ to an Sp(N) node contributes a fundamental of Sp(N), counted in this total.)

---

## 5. Asymptotic Freedom Conditions

$b_0^{(a)} > 0$ must hold for all valid $N$ at every node ($N \geq 2$ for SU(N), $N \geq 3$ for SO(N), $N \geq 1$ for Sp(N)). The exact formulas, using the bifundamental contributions from the table in §3, are:

**SU(N):**
$$b_0 = 3N - \frac{n_f + n_{\bar{f}}}{2} - N n_{\mathrm{adj}} - \frac{N+2}{2}(n_S + n_{\bar{S}}) - \frac{N-2}{2}(n_A + n_{\bar{A}}) - \sum_{\text{bif. edges}} T_{G_a}$$

**SO(N):**
$$b_0 = 3(N-2) - n_V - (N-2)n_{\mathrm{adj}} - (N+2)n_S - \sum_{\text{bif. edges}} T_{G_a}$$

**Sp(N):**
$$b_0 = 3(N+1) - \frac{n_f}{2} - (N+1)n_{\mathrm{adj}} - (N-1)n_A - \sum_{\text{bif. edges}} T_{G_a}$$

where the bifundamental contribution $T_{G_a}$ per edge is read from the table in §3.

---

## 6. Enumeration Strategy

### Step 1: Enumerate graphs
Enumerate all connected directed multigraphs up to isomorphism with $|V| \leq V_{\max}$. The maximum degree per node is bounded by requiring $b_0 > 0$ for some $N$: since each bifundamental edge contributes at least $T_{G_a} \geq N/2$ (from §3), the degree is finite for any fixed matter content.

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
- Check $b_0 > 0$ at each node (Module 1)
- Check SU(N) anomaly cancellation
- Check Sp(N) Witten anomaly condition
- (Optionally) check theory is connected and irreducible

---

## 7. Examples of Valid Quivers

### Example A: SU(N) SQCD with $N_f$ fundamentals
- 1 node (SU(N)), $N_f$ flavor fundamentals + $N_f$ anti-fundamentals
- AF: $b_0 = 3N - N_f > 0 \Rightarrow N_f < 3N$
- Anomaly: $N_f \cdot 1 + N_f \cdot(-1) = 0$ ✓

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
The complete classification of single-node 4D $\mathcal{N}=1$ theories (with gauge anomaly cancellation, UV asymptotic freedom, and non-trivial IR fixed point) is given in arXiv:2007.16165 and arXiv:2510.19136.

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
                    if check_af(q) and check_anomalies(q):
                        valid_quivers.append(q)

    return valid_quivers


def check_af(quiver):
    for node in quiver.nodes:
        b0 = compute_beta(quiver, node)
        if b0 <= 0:
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
