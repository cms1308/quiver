# Phase 2: Two-Node Summary (No-Veneziano Cut) - Research

**Researched:** 2026-04-21
**Focus:** Two-node quivers with $m \leq 3$ and NO Veneziano limit (`ven_all=0` or `class_id IS NULL`).
**Confidence:** HIGH

## Summary
The "no-Veneziano" cut represents a physically distinct sector of the two-node quiver landscape where the matter content does not scale with $N$ in a way that allows for a traditional $N_f/N_c$ fixed ratio. In these theories, the leading-order one-loop beta function coefficient $\alpha$ is identically zero ($\alpha = 0$), making the theories marginally asymptotically free at leading order ($b_0 = O(1)$).

**Striking Physical Universality:** All 105 universality classes in the no-Veneziano cut exhibit **exactly $a/c = 1.000$** at leading order in large $N$. This is a major physical signal, likely indicating $\mathcal{N}=2$ supersymmetry enhancement or a rigid $R$-charge configuration where the vector multiplet contribution dominates or forces $\text{Tr} R = 0$.

## Cut Statistics (m ≤ 3, no-ven = 0 OR class_id IS NULL)
- **Total theories:** 1369 (1264 classified with `ven=0` + 105 `class-NULL`).
- **Total classes:** 119 (105 with `ven_all=0` + 14 "mixed" classes where `ven_any=1, ven_all=0`).
- **$a/c$ Ratio:** 1.000 (universal across all no-Veneziano classes).

### Breakdown by Gauge Pair and Rank Multipliers:
- **SO-SO**: (1,1): 7 theories, (2,1): 6, (3,2): 18.
- **SO-Sp**: (1,1): 6, (2,1): 21, (3,1): 21, (3,2): 22.
- **SU-SO**: (1,1): 54, (1,2): 45, (1,3): 8, (2,3): 17.
- **SU-SU**: (1,1): 432, (2,1): 177, (3,1): 56, (3,2): 199.
- **SU-Sp**: (1,1): 54, (2,1): 126, (2,3): 8, (3,1): 15, (3,2): 17.
- **Sp-Sp**: (1,1): 28, (2,1): 18, (3,1): 6, (3,2): 21.

## Central Charge Spectrum
The $a/N^2$ values range from negative (non-unitary) to large positive values, despite the fixed $a/c = 1$.

- **Largest $a/N^2$**: Sp-Sp(3,2) Class 597 ($a/N^2 = 6.24$), SU-Sp(2,3) Class 278 ($a/N^2 = 5.34$).
- **Unitary Classes ($a/N^2 > 0$)**: The majority of theories flow to unitary SCFTs with $a=c$.
- **Non-Unitary Classes ($a/N^2 < 0$)**: Class 111 (SU-SU(2,1)) has $a/N^2 = -0.6328$. These represent non-unitary fixed points where the Hofman-Maldacena bounds are violated but the $a$-function stationary point exists.
- **Trivial/Type-I Classes ($a/N^2 = 0$)**: Class 1 (SU-SU(1,1)) and Class 117 (SU-SU(2,1)) represent theories that are likely IR-free or "Type-I" at leading order.

## The 14 "Mystery" Classes (877-899)
These classes have associated theories (ranging from 1 to 238 members) but **NULL central charges**. 
- **Interpretation**: These correspond to theories where a-maximization failed to find a solution (e.g., diverged or no stationary point in the unitary region).
- **Classification**: These are categorized as **Unclassified / Non-Conformal**.

## Structural Theorem Candidate
> "4D $\mathcal{N}=1$ two-node quivers without a Veneziano limit ($\alpha=0$) flow to $a=c$ SCFTs (when unitary) at large $N$ because the vanishing of the leading-order beta function forces the $R$-charges to satisfy $\text{Tr} R = 0$ at the fixed point."

## Conventions
- Same as Phase 1 & 2 original.
- Ordering: Alphabetical by gauge pair.

## Data Sources
- `quivers.db`: `theory` and `universality_class` tables.
- `results/two_node_classification.json`: Source for existing superconformal data.
