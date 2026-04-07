"""
a_maximization.py — Module 3

IR R-charges and central charges via Intriligator-Wecht a-maximization.

For a quiver at specific (N, N_f), builds the full chiral field content,
imposes anomaly-free constraints, and maximizes the trial a-function.

Reference: module3_a_maximization.md
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction

import numpy as np
from scipy.linalg import null_space
from scipy.optimize import minimize

from beta_functions import GaugeType, T_adj, T_rep, T_bifund
from quiver_generation import Quiver, chiral_excess_coeffs


# ── Dimension helpers ──────────────────────────────────────────────────────────

def dim_group(gauge_type: GaugeType, N: int) -> int:
    """Dimension of gauge group (= number of gaugino Weyl fermions)."""
    if gauge_type == "SU":
        return N * N - 1
    if gauge_type == "SO":
        return N * (N - 1) // 2
    if gauge_type == "Sp":
        return N * (2 * N + 1)
    raise ValueError(f"Unknown gauge type: {gauge_type!r}")


def dim_rep(gauge_type: GaugeType, rep: str, N: int) -> int:
    """Dimension of a single-node representation."""
    if gauge_type == "SU":
        if rep in ("fund", "antifund"):
            return N
        if rep == "adj":
            return N * N - 1
        if rep in ("S", "Sbar"):
            return N * (N + 1) // 2
        if rep in ("A", "Abar"):
            return N * (N - 1) // 2
    elif gauge_type == "SO":
        if rep == "V":
            return N
        if rep == "adj":
            return N * (N - 1) // 2
        if rep == "S":
            return N * (N + 1) // 2 - 1  # traceless symmetric
    elif gauge_type == "Sp":
        if rep == "fund":
            return 2 * N
        if rep == "adj":
            return N * (2 * N + 1)
        if rep == "A":
            return N * (2 * N - 1) - 1  # traceless antisymmetric
    raise ValueError(f"Unknown rep {rep!r} for gauge type {gauge_type!r}")


# ── ChiralField ────────────────────────────────────────────────────────────────

@dataclass
class ChiralField:
    """
    One type of chiral multiplet (all copies share the same R-charge).

    Attributes
    ----------
    label : str
        Human-readable identifier, e.g. "node0_adj", "node1_fund", "edge_0_1_pm".
    R_index : int
        Index into the R-charge vector (unique per field type).
    dim : int
        Total dimension = product of dim(r_a) over all gauge groups the field
        is charged under, times the multiplicity (number of copies).
    T_contributions : dict[int, Fraction]
        node_index → effective Dynkin index T_{G_a} at that node
        (already multiplied by multiplicity).
    """
    label: str
    R_index: int
    dim: int
    T_contributions: dict[int, Fraction]


# ── build_fields ───────────────────────────────────────────────────────────────

def build_fields(quiver: Quiver, N: int, N_f: int) -> list[ChiralField]:
    """
    Expand a Module 2 quiver at (N, N_f) into the full list of chiral fields.

    Each node i has effective rank N_i = rank_multipliers[i] * N.

    Three sources of fields:
    1. Node-level rank-2/adj matter (from quiver.node_matter)
    2. Fund-like matter derived from N_f and chiral excess δ
    3. Bifundamental matter (one ChiralField per edge)
    """
    fields: list[ChiralField] = []
    idx = 0  # R_index counter
    mults = quiver.rank_multipliers

    # 1. Rank-2/adj node matter
    for i, (g, matter) in enumerate(zip(quiver.gauge_types, quiver.node_matter)):
        N_i = mults[i] * N
        for rep, count in matter.items():
            if count == 0:
                continue
            d = count * dim_rep(g, rep, N_i)
            T = count * T_rep(g, rep, N_i)
            fields.append(ChiralField(
                label=f"node{i}_{rep}",
                R_index=idx,
                dim=d,
                T_contributions={i: T},
            ))
            idx += 1

    # 2. Fund-like matter
    for i, g in enumerate(quiver.gauge_types):
        N_i = mults[i] * N
        if g == "SU":
            a, b = chiral_excess_coeffs(quiver, i)
            delta = a * N + b  # chiral excess = n_f - n_fbar (N is the base rank)
            n_f    = N_f + max(delta, 0)
            n_fbar = N_f + max(-delta, 0)
            if n_f > 0:
                fields.append(ChiralField(
                    label=f"node{i}_fund",
                    R_index=idx,
                    dim=n_f * dim_rep(g, "fund", N_i),
                    T_contributions={i: n_f * T_rep(g, "fund", N_i)},
                ))
                idx += 1
            if n_fbar > 0:
                fields.append(ChiralField(
                    label=f"node{i}_antifund",
                    R_index=idx,
                    dim=n_fbar * dim_rep(g, "antifund", N_i),
                    T_contributions={i: n_fbar * T_rep(g, "antifund", N_i)},
                ))
                idx += 1
        elif g == "SO":
            if N_f > 0:
                fields.append(ChiralField(
                    label=f"node{i}_V",
                    R_index=idx,
                    dim=N_f * dim_rep(g, "V", N_i),
                    T_contributions={i: N_f * T_rep(g, "V", N_i)},
                ))
                idx += 1
        elif g == "Sp":
            n_fund = 2 * N_f
            if n_fund > 0:
                fields.append(ChiralField(
                    label=f"node{i}_fund",
                    R_index=idx,
                    dim=n_fund * dim_rep(g, "fund", N_i),
                    T_contributions={i: n_fund * T_rep(g, "fund", N_i)},
                ))
                idx += 1

    # 3. Bifundamental edges
    for e in quiver.edges:
        i, j = e.src, e.dst
        gi, gj = quiver.gauge_types[i], quiver.gauge_types[j]
        N_i, N_j = mults[i] * N, mults[j] * N
        T_i, T_j = T_bifund(gi, gj, N, N_i, N_j)
        d_i = _bifund_dim_at_node(gi, e.rep, N_i, side="src")
        d_j = _bifund_dim_at_node(gj, e.rep, N_j, side="dst")
        rep_label = e.rep.replace("+", "p").replace("-", "m")
        fields.append(ChiralField(
            label=f"edge_{i}_{j}_{rep_label}",
            R_index=idx,
            dim=d_i * d_j,
            T_contributions={i: T_i, j: T_j},
        ))
        idx += 1

    return fields


def _bifund_dim_at_node(gauge_type: GaugeType, rep: str, N: int, side: str) -> int:
    """Dimension of the fundamental-like representation at one endpoint of a bifundamental edge.
    N here is the effective rank (already multiplied by rank_mult)."""
    if gauge_type == "SU":
        return N  # both □ and □̄ have dim N
    if gauge_type == "SO":
        return N  # vector V
    if gauge_type == "Sp":
        return 2 * N  # fundamental f of USp(2N)
    raise ValueError(f"Unknown gauge type: {gauge_type!r}")


# ── Anomaly constraint matrix ──────────────────────────────────────────────────

def anomaly_matrix(
    fields: list[ChiralField],
    quiver: Quiver,
    N: int,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Build (A, b) encoding the anomaly-free condition at each gauge node:

        T(adj_a) + Σ_i T_a(field_i) * (R_i - 1) = 0

    Rearranged to A @ R = b:
        A[a, R_index] = T_a(field)   (summed if multiple fields share R_index)
        b[a] = Σ_i T_a(field_i) - T(adj_a)

    Each node a uses effective rank N_a = rank_multipliers[a] * N.
    """
    n_nodes = quiver.n_nodes
    n_R = len(fields)
    A = np.zeros((n_nodes, n_R))
    b = np.zeros(n_nodes)

    mults = quiver.rank_multipliers
    for a, g in enumerate(quiver.gauge_types):
        N_a = mults[a] * N
        b[a] = float(-T_adj(g, N_a))
        for f in fields:
            T_a = f.T_contributions.get(a)
            if T_a is not None:
                A[a, f.R_index] += float(T_a)
                b[a] += float(T_a)

    return A, b


# ── Trial central charges ──────────────────────────────────────────────────────

def _traces(R_values: np.ndarray, fields: list[ChiralField], quiver: Quiver, N: int) -> tuple[float, float]:
    """Return (Tr R, Tr R³) including gaugino contributions.
    Each node a uses effective rank N_a = rank_multipliers[a] * N."""
    mults = quiver.rank_multipliers
    tr_R = sum(dim_group(g, m * N) for g, m in zip(quiver.gauge_types, mults))
    tr_R3 = float(tr_R)
    for f in fields:
        r = R_values[f.R_index]
        tr_R  += f.dim * (r - 1)
        tr_R3 += f.dim * (r - 1) ** 3
    return float(tr_R), float(tr_R3)


def a_trial(R_values: np.ndarray, fields: list[ChiralField], quiver: Quiver, N: int) -> float:
    """a = (3/32)(3 Tr R³ - Tr R)"""
    tr_R, tr_R3 = _traces(R_values, fields, quiver, N)
    return (3 / 32) * (3 * tr_R3 - tr_R)


def c_trial(R_values: np.ndarray, fields: list[ChiralField], quiver: Quiver, N: int) -> float:
    """c = (1/32)(9 Tr R³ - 5 Tr R)"""
    tr_R, tr_R3 = _traces(R_values, fields, quiver, N)
    return (1 / 32) * (9 * tr_R3 - 5 * tr_R)


# ── Gauge-invariant operator checks ───────────────────────────────────────────

def gauge_invariant_ops(
    fields: list[ChiralField],
    quiver: Quiver,
    R_values: np.ndarray,
) -> dict[str, tuple[list[int], float]]:
    """
    Return the simple gauge-invariant operators and their R-charges.

    Each entry: label → (list_of_R_indices, R_value_of_operator)

    Operators checked:
      - SU node: meson Q*Qtilde (fund + antifund), R = R_fund + R_antifund
      - SU bifundamental pair Q, Qtilde (opposite directions): R = R_Q + R_Qtilde
      - Adjoint bilinear Tr(Φ²): R = 2*R_adj
      - S+Sbar, A+Abar bilinears: R = R_S + R_Sbar, R_A + R_Abar
      - SO vector bilinear V·V: R = 2*R_V
      - Sp fundamental bilinear f·Ω·f: R = 2*R_f

    The unitarity bound R[O] ≥ 2/3 applies to these operators, not to the
    elementary fields themselves.
    """
    by_label: dict[str, ChiralField] = {f.label: f for f in fields}
    ops: dict[str, tuple[list[int], float]] = {}

    for i, g in enumerate(quiver.gauge_types):
        # Fund + antifund meson
        fl, afl = f"node{i}_fund", f"node{i}_antifund"
        if g == "SU" and fl in by_label and afl in by_label:
            r_idx = [by_label[fl].R_index, by_label[afl].R_index]
            ops[f"meson_{i}"] = (r_idx, sum(R_values[k] for k in r_idx))

        # Self-conjugate bilinears (Tr Φ², S·Sbar, A·Abar, V·V, f·Ω·f)
        bilinear_pairs = [
            ("adj",  "adj",  f"adj_bilinear_{i}"),
            ("S",    "Sbar", f"SSbar_{i}"),
            ("A",    "Abar", f"AAbar_{i}"),
            ("V",    "V",    f"VV_{i}"),   # SO
        ]
        if g == "Sp":
            bilinear_pairs.append(("fund", "fund", f"Sp_ff_{i}"))

        for rep1, rep2, op_label in bilinear_pairs:
            l1, l2 = f"node{i}_{rep1}", f"node{i}_{rep2}"
            if rep1 == rep2:
                if l1 in by_label:
                    r_idx = [by_label[l1].R_index, by_label[l1].R_index]
                    ops[op_label] = (r_idx, 2 * R_values[by_label[l1].R_index])
            else:
                if l1 in by_label and l2 in by_label:
                    r_idx = [by_label[l1].R_index, by_label[l2].R_index]
                    ops[op_label] = (r_idx, R_values[r_idx[0]] + R_values[r_idx[1]])

    # Bifundamental pairs: edge (i,j,"+-") paired with edge (j,i,"+-")
    edge_fields = [f for f in fields if f.label.startswith("edge_")]
    for f1 in edge_fields:
        parts1 = f1.label.split("_")
        src1, dst1, rep1 = parts1[1], parts1[2], "_".join(parts1[3:])
        if rep1 != "pm":
            continue
        for f2 in edge_fields:
            parts2 = f2.label.split("_")
            src2, dst2, rep2 = parts2[1], parts2[2], "_".join(parts2[3:])
            if rep2 != "pm" or src2 != dst1 or dst2 != src1:
                continue
            key = f"bif_meson_{min(src1,src2)}_{max(src1,src2)}"
            if key not in ops:
                r_idx = [f1.R_index, f2.R_index]
                ops[key] = (r_idx, R_values[f1.R_index] + R_values[f2.R_index])

    return ops


# ── AMaxResult ─────────────────────────────────────────────────────────────────

@dataclass
class AMaxResult:
    fields: list[ChiralField]        # all chiral fields used in the computation
    R_charges: dict[str, float]      # field.label → R-charge value
    a: float                         # central charge a*
    c: float                         # central charge c*
    unitarity_ok: bool               # all gauge-invariant operators have R ≥ 2/3


# ── Core a-maximization ────────────────────────────────────────────────────────

def a_maximize(quiver: Quiver, N: int, N_f: int) -> AMaxResult:
    """
    Determine IR R-charges by maximizing the trial a-function subject to
    anomaly-free constraints, at specific (N, N_f).

    Algorithm:
    1. Build chiral fields at (N, N_f)
    2. Construct anomaly constraint matrix A @ R = b
    3. Particular solution R0 via lstsq; null space F via scipy
    4. Parameterize R = R0 + F @ s and maximize a_trial over s
    5. Compute a*, c*, check unitarity
    """
    fields = build_fields(quiver, N, N_f)
    A, b = anomaly_matrix(fields, quiver, N)

    # Particular solution
    R0, *_ = np.linalg.lstsq(A, b, rcond=None)

    # Null space (flavor directions)
    F = null_space(A)  # shape (n_R, n_free)
    n_free = F.shape[1]

    def neg_a(s: np.ndarray) -> float:
        R = R0 + F @ s
        return -a_trial(R, fields, quiver, N)

    if n_free == 0:
        s_opt = np.zeros(0)
    else:
        res = minimize(neg_a, x0=np.zeros(n_free), method="BFGS")
        s_opt = res.x

    R_opt = R0 + F @ s_opt
    a_val = a_trial(R_opt, fields, quiver, N)
    c_val = c_trial(R_opt, fields, quiver, N)

    R_charges = {f.label: float(R_opt[f.R_index]) for f in fields}
    gi_ops = gauge_invariant_ops(fields, quiver, R_opt)
    unitarity_ok = all(r >= 2 / 3 - 1e-9 for _, r in gi_ops.values())

    return AMaxResult(fields=fields, R_charges=R_charges, a=a_val, c=c_val, unitarity_ok=unitarity_ok)


def a_maximize_with_decoupling(quiver: Quiver, N: int, N_f: int) -> AMaxResult:
    """
    a-maximization with iterative unitarity decoupling.

    The unitarity bound R[O] ≥ 2/3 applies to gauge-invariant operators O,
    not to elementary (gauge-charged) fields.

    When a gauge-invariant operator O = Φ_1...Φ_k has R[O] < 2/3, it decouples
    as a free field. We enforce this by adding the linear constraint
        Σ R[Φ_i for i in O] = 2/3
    and redo a-maximization. Iterate until all gauge-invariant operators satisfy
    R[O] ≥ 2/3.
    """
    fields = build_fields(quiver, N, N_f)
    n_R = len(fields)

    # extra_constraints: list of (weight_vector, rhs) to add to anomaly system
    # Each encodes Σ_i w_i * R_i = rhs for a decoupled operator.
    extra_constraints: list[tuple[np.ndarray, float]] = []
    R_opt = np.zeros(n_R)

    for _ in range(n_R + 1):
        A, b = anomaly_matrix(fields, quiver, N)

        if extra_constraints:
            extra_A = np.array([w for w, _ in extra_constraints])
            extra_b = np.array([rhs for _, rhs in extra_constraints])
            A = np.vstack([A, extra_A])
            b = np.concatenate([b, extra_b])

        R0, *_ = np.linalg.lstsq(A, b, rcond=None)
        F = null_space(A)
        n_free = F.shape[1]

        def neg_a(s: np.ndarray) -> float:
            return -a_trial(R0 + F @ s, fields, quiver, N)

        s_opt = minimize(neg_a, x0=np.zeros(n_free), method="BFGS").x if n_free > 0 else np.zeros(0)
        R_opt = R0 + F @ s_opt

        gi_ops = gauge_invariant_ops(fields, quiver, R_opt)
        violating = [(label, r_idx) for label, (r_idx, r_val) in gi_ops.items()
                     if r_val < 2 / 3 - 1e-9]
        if not violating:
            break

        # Add one constraint per violating operator: Σ R_i = 2/3
        already = {tuple(sorted(w.nonzero()[0])) for w, _ in extra_constraints}
        added = False
        for label, r_idx in violating:
            key = tuple(sorted(set(r_idx)))
            if key in already:
                continue
            w = np.zeros(n_R)
            for k in r_idx:
                w[k] += 1.0
            extra_constraints.append((w, 2 / 3))
            already.add(key)
            added = True
        if not added:
            break  # no new constraints can be added; stop

    a_val = a_trial(R_opt, fields, quiver, N)
    c_val = c_trial(R_opt, fields, quiver, N)
    R_charges = {f.label: float(R_opt[f.R_index]) for f in fields}
    gi_ops = gauge_invariant_ops(fields, quiver, R_opt)
    unitarity_ok = all(r >= 2 / 3 - 1e-9 for _, r in gi_ops.values())

    return AMaxResult(fields=fields, R_charges=R_charges, a=a_val, c=c_val, unitarity_ok=unitarity_ok)


# ── Validation ─────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    from quiver_generation import Edge

    # ── Example 1: SU(N) SQCD ─────────────────────────────────────────────────
    # With T(fund) = 1/2 convention (used throughout), the anomaly-free condition gives:
    #   T(adj) + N_f*T(fund)*(R_f-1) + N_f*T(antifund)*(R_fbar-1) = 0
    #   N + N_f*(R-1) = 0  (symmetric R_f = R_fbar = R)
    #   R = 1 - N/N_f
    # Conformal window: N_f > 3N/2 (lower bound from meson R_M = 2R ≥ 2/3 → R ≥ 1/3).
    # Upper bound: N_f < 3N (AF). As N_f → 3N, R → 2/3 (free limit).
    # Note: individual quarks have R < 2/3; unitarity bound applies to gauge-invariant operators.
    print("=== SU(5) SQCD, N_f=10 ===")
    q = Quiver(["SU"])
    res = a_maximize(q, N=5, N_f=10)
    for label, R in res.R_charges.items():
        expect = 1 - 5 / 10
        print(f"  {label}: R = {R:.6f}  (expect {expect:.4f} = 1 - N/N_f)")
    print(f"  a = {res.a:.6f},  c = {res.c:.6f}")

    print()

    # Scan: as N_f → 3N=15, R → 2/3; at N_f=3N/2=7.5, meson hits unitarity bound
    print("=== SU(5) SQCD conformal window scan ===")
    print(f"  {'N_f':>4}  {'R_quark':>8}  {'R_meson':>8}")
    for N_f in [8, 10, 12, 14]:
        r = a_maximize(q, N=5, N_f=N_f)
        R_q = list(r.R_charges.values())[0]
        print(f"  {N_f:>4}  {R_q:>8.4f}  {2*R_q:>8.4f}  (expect R = {1-5/N_f:.4f})")

    print()

    # ── Example 2: SU(5) SQCD below conformal window — decoupling ─────────────
    # N_f=4 < 3N/2=7.5: R_quark = 1 - 5/4 = -0.25, meson R = -0.5 < 2/3.
    # Decoupling adds constraint R_fund + R_antifund = 2/3 (pins meson).
    # This is inconsistent with the anomaly constraint (theory is not an SCFT);
    # lstsq finds the best-fit compromise and unitarity_ok = False signals no SCFT.
    print("=== SU(5) SQCD, N_f=4 (below conformal window) — with decoupling ===")
    res2 = a_maximize_with_decoupling(q, N=5, N_f=4)
    for label, R in res2.R_charges.items():
        print(f"  {label}: R = {R:.6f}")
    print(f"  unitarity_ok = {res2.unitarity_ok}  (False: no SCFT in this phase)")

    print()

    # ── Example 3: SU(N)² circular quiver with extra flavors ───────────────────
    # Balanced +- (one each way) + N_f=3 flavors per node.
    # Anomaly at each node: b = N_f (from fund+antifund). One free parameter per quiver.
    print("=== SU(5)² circular quiver (balanced +-), N_f=3 ===")
    q2 = Quiver(["SU", "SU"], [Edge(0, 1, "+-"), Edge(1, 0, "+-")])
    res3 = a_maximize(q2, N=5, N_f=3)
    for label, R in res3.R_charges.items():
        print(f"  {label}: R = {R:.6f}")
    print(f"  a = {res3.a:.6f},  c = {res3.c:.6f}")

    print()

    # ── Example 4: SU(N) + 1 adj + N_f flavors ─────────────────────────────────
    # Anomaly: N + N*(R_adj-1) + N_f*(R_f-1) + N_f*(R_fbar-1) = 0
    # One constraint, three unknowns → two free parameters.
    # By flavor symmetry R_f = R_fbar: one remaining free parameter → nontrivial a-max.
    print("=== SU(5) + 1 adj + N_f=8 flavors ===")
    q3 = Quiver(["SU"], node_matter=[{"adj": 1}])
    res4 = a_maximize(q3, N=5, N_f=8)
    for label, R in res4.R_charges.items():
        print(f"  {label}: R = {R:.6f}")
    print(f"  a = {res4.a:.6f},  c = {res4.c:.6f},  unitarity_ok = {res4.unitarity_ok}")
