# Module 1: One-Loop Beta Function Coefficients

## Goal
For each gauge node $G_a$ in a quiver, compute the one-loop beta function coefficient $b_0^{(a)}$ and determine whether the node (and the full theory) is asymptotically free at UV.

---

## 1. The One-Loop Beta Function

For a 4D N=1 supersymmetric gauge theory with gauge group $G$ and chiral superfields in representations $r_i$, the one-loop beta function is:

$$\mu \frac{dg}{d\mu} = -\frac{g^3}{16\pi^2} b_0 + \mathcal{O}(g^5)$$

$$\boxed{b_0 = 3\,T(\mathrm{adj}) - \sum_i T(r_i)}$$

where $T(r)$ is the **Dynkin index** of representation $r$, defined by:
$$\mathrm{tr}_r[T^a T^b] = T(r)\,\delta^{ab}$$

**Asymptotic freedom** requires $b_0 > 0$.

The sum runs over all chiral superfields charged under $G$. For a bifundamental field $(\square_a, \bar{\square}_b)$ shared between two nodes, it contributes to **both** $b_0^{(a)}$ and $b_0^{(b)}$.

---

## 2. Dynkin Indices

### SU(N) — normalization: $T(\square) = \frac{1}{2}$

| Representation | Symbol | Dimension | Dynkin Index $T$ | Anomaly $A$ |
|----------------|--------|-----------|-------------------|-------------|
| Fundamental | $\square$ | $N$ | $\frac{1}{2}$ | $+1$ |
| Anti-fundamental | $\bar{\square}$ | $N$ | $\frac{1}{2}$ | $-1$ |
| Adjoint | $\mathrm{adj}$ | $N^2 - 1$ | $N$ | $0$ |
| Rank-2 symmetric | $S$ | $\frac{N(N+1)}{2}$ | $\frac{N+2}{2}$ | $+(N+4)$ |
| Rank-2 sym conjugate | $\bar{S}$ | $\frac{N(N+1)}{2}$ | $\frac{N+2}{2}$ | $-(N+4)$ |
| Rank-2 antisymmetric | $A$ | $\frac{N(N-1)}{2}$ | $\frac{N-2}{2}$ | $+(N-4)$ |
| Rank-2 antisym conjugate | $\bar{A}$ | $\frac{N(N-1)}{2}$ | $\frac{N-2}{2}$ | $-(N-4)$ |

The anomaly coefficient $A(r)$ enters gauge anomaly cancellation (see Module 2).

**Beta function for SU(N):**
$$b_0^{SU(N)} = 3N - \frac{n_f + n_{\bar{f}}}{2} - N\,n_{\mathrm{adj}} - \frac{N+2}{2}\,(n_S + n_{\bar{S}}) - \frac{N-2}{2}\,(n_A + n_{\bar{A}}) - [\text{bifundamental contributions}]$$

where $n_f, n_{\bar{f}}, n_{\mathrm{adj}}, n_S, n_{\bar{S}}, n_A, n_{\bar{A}}$ are the numbers of fundamentals, anti-fundamentals, adjoints, symmetric, conjugate symmetric, antisymmetric, and conjugate antisymmetric chiral multiplets. Note that $\bar{S}$ and $\bar{A}$ have the same Dynkin index as $S$ and $A$ respectively, so their contributions to $b_0$ are identical.

---

### SO(N) — normalization: $T(V) = 1$

| Representation | Symbol | Dimension | Dynkin Index $T$ |
|----------------|--------|-----------|-------------------|
| Vector (fundamental) | $V$ | $N$ | $1$ |
| Adjoint (= rank-2 antisym) | $\mathrm{adj}$ | $\frac{N(N-1)}{2}$ | $N-2$ |
| Rank-2 symmetric (traceless) | $S^2$ | $\frac{N(N+1)}{2}-1$ | $N+2$ |
| Spinor | $\mathbf{s}$ | $2^{\lfloor N/2 \rfloor - 1}$ | grows exponentially |

> **Note:** Spinor representations grow exponentially with $N$ and are **excluded** in the large N classification since they violate large N factorization.

SO(N) is automatically free of gauge anomalies (no $A(r)$ constraint needed).

**Beta function for SO(N):**
$$b_0^{SO(N)} = 3(N-2) - n_V \cdot 1 - n_{\mathrm{adj}}(N-2) - n_S(N+2) - [\text{bifundamental contributions}]$$

---

### USp(2N) = Sp(N) — normalization: $T(f) = \frac{1}{2}$

The group USp(2N) has rank $N$; its fundamental representation has dimension $2N$.

| Representation | Symbol | Dimension | Dynkin Index $T$ |
|----------------|--------|-----------|-------------------|
| Fundamental | $f$ | $2N$ | $\frac{1}{2}$ |
| Adjoint (= rank-2 symmetric) | $\mathrm{adj}$ | $N(2N+1)$ | $N+1$ |
| Rank-2 antisymmetric (traceless) | $A$ | $N(2N-1)-1$ | $N-1$ |

> **Note:** For USp(2) $\cong$ SU(2), the antisymmetric rank-2 is a singlet and plays no role.

For Sp(N), the number of fundamental representations must be **even** to avoid the Witten anomaly (a global gauge anomaly).

**Beta function for Sp(N):**
$$b_0^{Sp(N)} = 3(N+1) - \frac{n_f}{2} - n_{\mathrm{adj}}(N+1) - n_A(N-1) - [\text{bifundamental contributions}]$$

---

## 3. Restriction to Bifundamental Inter-Node Matter

### Why only bifundamentals?

In principle, a chiral superfield can transform in any representation of the product gauge group $G_a \times G_b$. For example, a field in $(\square_a, \mathrm{adj}_b)$ is gauge-invariant and consistent as a quantum field theory. However, **for theories admitting a large N limit, only bifundamental representations are allowed as inter-node matter.** The reason is a strict constraint from asymptotic freedom at large N.

The Dynkin index contribution of a field in representation $(r_a, r_b)$ of $G_a \times G_b$ to node $G_a$'s beta function is:
$$\Delta b_0^{(a)} = T_{G_a}(r_a) \cdot \dim(r_b)$$

The factor $\dim(r_b)$ counts the multiplicity: for each color state of $G_b$, there is an independent copy of $r_a$ under $G_a$.

**Scaling with N for different choices of $r_b$:**

| $r_b$ (under $G_b = \mathrm{SU}(N)$) | $\dim(r_b)$ | $\Delta b_0^{(a)}$ scaling |
|--------------------------------------|-------------|----------------------------|
| Fundamental $\square$ | $N$ | $\sim N$ |
| Adjoint $\mathrm{adj}$ | $N^2 - 1$ | $\sim N^2$ |
| Rank-2 symmetric $S$ | $N(N+1)/2$ | $\sim N^2$ |
| Rank-2 antisymmetric $A$ | $N(N-1)/2$ | $\sim N^2$ |

Since $3T(\mathrm{adj}_a) \sim 3N$ at large N, any inter-node matter with $\dim(r_b) \sim N^2$ would contribute $\Delta b_0^{(a)} \sim N^2$ to node $a$'s beta function — overwhelming the $3N$ gauge contribution and making $b_0^{(a)} \to -\infty$. Such theories are **never asymptotically free at large N**.

> **Conclusion:** For a theory to remain asymptotically free as $N \to \infty$, every inter-node matter field must transform in representations $r_a$, $r_b$ such that both $\dim(r_a)$ and $\dim(r_b)$ scale at most as $N^1$. This restricts inter-node matter to:
> - $r_a, r_b \in \{\square, \bar{\square}, V, f\}$ (fundamental/vector representations of each group)

This is precisely the **bifundamental** matter: the field transforms in the fundamental (or anti-fundamental/vector) of each of the two gauge groups it connects. Higher representations — e.g. $(\square_a, \mathrm{adj}_b)$, $(\mathrm{adj}_a, \mathrm{adj}_b)$, $(S_a, \square_b)$ — are excluded by large N asymptotic freedom.

**Node matter** (a field charged under only one gauge group $G_a$) is not subject to this cross-node constraint. Adjoints, symmetric, and antisymmetric representations at a single node are allowed because their Dynkin index contribution scales at most as $N$ (see Section 2 and 4).

---

## 4. Bifundamental Contributions (Large N Allowed Cases)

When two gauge nodes $G_a$ and $G_b$ are connected by a bifundamental edge, the matter field is charged under both groups. The contribution to each node's beta function:

### SU(N) × SU(N) bifundamental $(\square_a, \bar{\square}_b)$:
- To node $a$: $T_{G_a} = N \cdot \frac{1}{2} = \frac{N}{2}$  (there are $N$ fundamentals of $G_a$, one per color of $G_b$)
- To node $b$: $T_{G_b} = N \cdot \frac{1}{2} = \frac{N}{2}$

### SU(N) × SU(N) bifundamental $(\square_a, \square_b)$:
- To node $a$: $T_{G_a} = N \cdot \frac{1}{2} = \frac{N}{2}$  (there are $N$ fundamentals of $G_a$, one per color of $G_b$)
- To node $b$: $T_{G_b} = N \cdot \frac{1}{2} = \frac{N}{2}$

Note: the Dynkin index contributions to $b_0$ are identical for $(\square, \bar{\square})$ and $(\square, \square)$. They differ in their **gauge anomaly** contributions (see Module 2) and give inequivalent chiral theories.

### SU(N) × SO(N) bifundamental $(\square_{SU}, V_{SO})$:
- To SU(N) node: $T_{SU} = N \cdot \frac{1}{2} = \frac{N}{2}$ (N fundamentals, one per SO(N) vector index)
- To SO(N) node: $T_{SO} = N \cdot 1 = N$ (N vectors, one per SU(N) color)

### SU(N) × Sp(N) bifundamental $(\square_{SU}, f_{Sp})$:
- To SU(N) node: $T_{SU} = 2N \cdot \frac{1}{2} = N$
- To Sp(N) node: $T_{Sp} = N \cdot \frac{1}{2} = \frac{N}{2}$

### SO(N) × Sp(N) bifundamental $(V_{SO}, f_{Sp})$:
- To SO(N) node: $T_{SO} = 2N \cdot 1 = 2N$
- To Sp(N) node: $T_{Sp} = N \cdot \frac{1}{2} = \frac{N}{2}$

---

## 5. Large N Asymptotics

At large N, the leading-N behavior of $b_0$ determines whether the theory remains AF.

| Group | $3T(\mathrm{adj})$ at large $N$ | Node matter contribution per field |
|-------|-------------------------------|----------------------------|
| SU(N) | $3N$ | $\square, \bar{\square}$: $\frac{1}{2} = \mathcal{O}(1)$; $S, \bar{S}, A, \bar{A}$: $\sim \frac{N}{2}$; adj: $\sim N$ |
| SO(N) | $3N$ | vector $V$: $1 = \mathcal{O}(1)$; sym $S$: $\sim N$; adj: $\sim N$ |
| Sp(N) | $3N$ | fund $f$: $\frac{1}{2} = \mathcal{O}(1)$; antisym $A$: $\sim N$; adj: $\sim N$ |

**Large N asymptotic freedom requires that the $\mathcal{O}(N)$ part of $b_0$ is positive:**

For **SU(N)** with $n_{\mathrm{adj}}$ adjoints, $n_S + n_{\bar{S}}$ symmetric (plus conjugates), $n_A + n_{\bar{A}}$ antisymmetric (plus conjugates), and $k$ bifundamental neighbors:
$$b_0 \approx N \left(3 - n_{\mathrm{adj}} - \frac{(n_S + n_{\bar{S}}) + (n_A + n_{\bar{A}})}{2} - \frac{k}{2}\right) + \mathcal{O}(1)$$

Large N AF condition:
$$n_{\mathrm{adj}} + \frac{(n_S + n_{\bar{S}}) + (n_A + n_{\bar{A}})}{2} + \frac{k}{2} < 3$$

where $k$ = number of bifundamental edges to other nodes (each contributing $N/2$ at large N regardless of the neighbor's gauge group).

---

## 6. Algorithm

```
for each quiver candidate Q:
    for each gauge node a with group G_a:
        b0[a] = 3 * T_adj(G_a)
        for each matter field i charged under G_a:
            b0[a] -= T_{G_a}(r_i)
        if b0[a] <= 0:
            reject Q  # not asymptotically free
    accept Q as UV-AF
```

### Input
- Quiver graph: nodes labeled by gauge group type, edges labeled by bifundamental type
- Additional node matter: adjoints, symmetric, antisymmetric at each node

### Output
- $b_0^{(a)}$ for each node
- Boolean: is the full quiver UV-asymptotically free?
- Boolean: does it remain AF at large N?


---

## 8. Notes on the NSVZ Exact Beta Function

At higher loops (and exactly via NSVZ):
$$\beta(g) = -\frac{g^3}{16\pi^2}\,\frac{3T(\mathrm{adj}) - \sum_i T(r_i)(1 - \gamma_i)}{1 - T(\mathrm{adj})\,g^2/(8\pi^2)}$$

where $\gamma_i$ is the anomalous dimension of field $i$. At the IR fixed point, $\beta(g) = 0$ (numerator = 0) and $\gamma_i = 3 R_i - 2$ (from superconformal algebra), recovering the condition:
$$3T(\mathrm{adj}) = \sum_i T(r_i)(3 - 3R_i) \quad \Leftrightarrow \quad \mathrm{Tr}[R\,G^2] = 0$$

This links the beta function vanishing to the anomaly-free R-symmetry condition used in a-maximization (Module 3).
