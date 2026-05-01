# Phase 2: Two-Node Summary (No-Veneziano Cut) - Summary

**Focus:** Two-node quivers with $m \leq 3$ and NO Veneziano limit (`ven_all=0` or `class_id IS NULL`).
**Status:** Complete
**Key Finding:** $a/c = 1.000$ Universality

## Performance
- **Duration:** Redo (2026-04-21)
- **Theories Analysed:** 1369
- **Classes Analysed:** 119

## Key Results
1.  **Universality of $a/c = 1.0$**: Every single class in the no-Veneziano cut exhibits exactly $a=c$ at leading order in large $N$. This confirms that theories with rigid matter content ($O(1)$ fields) flow to $a=c$ fixed points, likely due to $\mathcal{N}=2$ enhancement or vanishing leading-order beta functions.
2.  **Central Charge Spectrum**:
    - **Unitary range**: up to $a/N^2 = 6.24$ (Sp-Sp).
    - **Non-Unitary**: 31 classes with $a/N^2 < 0$. Class 111 (SU-SU) is a representative non-unitary fixed point ($a/N^2 = -0.6328$).
    - **Trivial**: Classes with $a=c=0$ identified as Type-I / IR-free.
3.  **Failed a-maximization**: 14 "Mystery" classes (877-899) identified as non-conformal or pathological flows where a-max failed to converge.

## Deliverables
- `results/no_veneziano_classification.json`: Data for all 1369 theories.
- `results/tables/no_veneziano_summary.md`: Detailed classification tables grouped by gauge pair.

## Physics Interpretation
The "no-Veneziano" condition corresponds to $\alpha = 0$ in the one-loop beta function $b_0 = \alpha N + \beta$. The vanishing of the $O(N)$ term forces the $R$-charges to satisfy the $\text{Tr} R = 0$ condition at the IR fixed point, which in turn leads to the exact equality of the central charges $a$ and $c$ at large $N$.
